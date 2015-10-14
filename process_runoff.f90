program process_runoff
! Take the weights file and the model file and write out our new file
! i.e. spread the runoffs according to a cunning plan!
!
!
use iso_fortran_env
use netcdf
implicit none

type :: nc_info_type
   character(len=128)  :: fname
   integer(kind=int32) :: ncid
   integer(kind=int32) :: vid
   integer(kind=int32) :: tid
   integer(kind=int32) :: idim
   integer(kind=int32) :: jdim
   integer(kind=int32) :: nrecs
end type nc_info_type

type :: model_type
   integer(int32):: npts
   integer(int32),allocatable,dimension(:) :: idx  !index of source array
   integer(int32),allocatable,dimension(:) :: i    !i index of 2d model grid
   integer(int32),allocatable,dimension(:) :: j    !j index of model grid
   real(real64),allocatable,dimension(:)   :: x    !x of model grid
   real(real64),allocatable,dimension(:)   :: y    !y of model grid
   real(real64),allocatable,dimension(:)   :: weight    ! factor to multiply source  runoff
end type model_type
   
type :: source_type
   integer(int32):: npts
   integer(int32),allocatable,dimension(:) :: i
   integer(int32),allocatable,dimension(:) :: j
   real(real64),allocatable,dimension(:)   :: x
   real(real64),allocatable,dimension(:)   :: y
   real(real64),allocatable,dimension(:)   :: area
   real(real64),allocatable,dimension(:)   :: dis
   integer(int32),allocatable,dimension(:) :: idx
end type source_type

type(model_type) :: model
type(source_type) :: source
type(nc_info_type) :: source_info,model_info

integer(kind=int32) :: nargs

integer(kind=int32) :: i,irecs


real(kind=real32),allocatable,dimension(:,:)   :: runoff_source,runoff_model         ! 
real(kind=real64),allocatable,dimension(:,:)   :: area_source, area_model
real(kind=real64),allocatable,dimension(:)     :: lon_model, lat_model

character(len=128)   :: carg,carg1

nargs=command_argument_count()
if ( nargs < 4 ) then
   write(*,*) 'ERROR: Usage :: process_runoff -i infile -o outfile'
   stop 1
endif
do i = 1,nargs,2
   call get_command_argument(i,carg)
   call get_command_argument(i+1,carg1)
   select case(trim(carg))
   case('-i')
       source_info%fname=trim(carg1)
   case('-o')
       model_info%fname=trim(carg1)
! Placeholder
   case DEFAULT
        write(*,*) 'Unknown argument',trim(carg)
        stop
   end select
enddo

!call read_connection_file_nn(source,coast,coast_nn,qc_nn)
write(*,*) 'Reading weights'
call read_weights_file(source,model)
write(*,*) 'Reading source info'
call get_source_info(source_info)
write(*,*) 'Reading model info'
call get_model_info(model_info)

allocate(runoff_source(source_info%idim,source_info%jdim))
allocate(runoff_model(model_info%idim,model_info%jdim))
! Equal all points weights by their area 

write(*,*) 'Setting up output file'
call setup_model_file(source_info,model_info)

do irecs = 1,source_info%nrecs
    write(*,*) 'Record',irecs
    call process_record(source_info,model_info,source,model,runoff_source,runoff_model,irecs)
    write(*,*) 'global area sum of source runoff = ', sum(runoff_source*area_source)
    write(*,*) 'global area sum of target runoff = ', sum(runoff_model*area_model)
enddo
call handle_error(nf90_close(source_info%ncid))
call handle_error(nf90_close(model_info%ncid))

write(*,*) 'Done'

contains

subroutine process_record(source_info,model_info,source,model,runoff_source,runoff_model,record)
   type(nc_info_type),intent(in) :: source_info
   type(nc_info_type),intent(in) :: model_info
   type(source_type),intent(in) :: source
   type(model_type),intent(in) :: model
   real(kind=real32),intent(out),dimension(:,:)   :: runoff_source,runoff_model
   integer(kind=int32),intent(in)  :: record
   integer(kind=int32)  :: i, i_s, j_s
   real(kind=real64)    :: time

   runoff_model = 0.0_real32
   call handle_error(nf90_get_var(source_info%ncid,source_info%vid,runoff_source,start=(/1,1,record/), &
                                  count=(/source_info%idim,source_info%jdim,1/)))
   call handle_error(nf90_get_var(source_info%ncid,source_info%tid,time,start=(/record/)))
   do i = 1,model%npts
      i_s = source%i(model%idx(i))
      j_s = source%j(model%idx(i))
      runoff_model(model%i(i),model%j(i)) = runoff_model(model%i(i),model%j(i)) + runoff_source(i_s,j_s)*model%weight(i)
   enddo
   call handle_error(nf90_put_var(model_info%ncid,model_info%vid,runoff_model,start=(/1,1,record/), &
                                  count=(/model_info%idim,model_info%jdim,1/)))
   call handle_error(nf90_put_var(model_info%ncid,model_info%tid,time,start=(/record/)))
end subroutine process_record

subroutine setup_model_file(source_info,model_info)
   type(nc_info_type),intent(inout) :: source_info
   type(nc_info_type),intent(inout) :: model_info
   integer(kind=int32)              :: id_x,id_y,id_t, id_lon, id_lat, id_area
   character(len=32)::dummy
   integer:: nats,i
   call handle_error(nf90_create(trim(model_info%fname),ior(nf90_netcdf4,nf90_clobber),model_info%ncid))
   call handle_error(nf90_def_dim(model_info%ncid,'nx',model_info%idim,id_x))
   call handle_error(nf90_def_dim(model_info%ncid,'ny',model_info%jdim,id_y))
   call handle_error(nf90_def_dim(model_info%ncid,'time',nf90_unlimited,id_t))
   call handle_error(nf90_def_var(model_info%ncid,'nx',nf90_double,(/ id_x /), id_lon))
   call handle_error(nf90_put_att(model_info%ncid, id_lon, "long_name", "longitude"))
   call handle_error(nf90_put_att(model_info%ncid, id_lon, "units", "degrees_east"))
   call handle_error(nf90_put_att(model_info%ncid, id_lon, "cartesian_axis", "X"))
   call handle_error(nf90_def_var(model_info%ncid,'ny',nf90_double,(/ id_y /), id_lat))
   call handle_error(nf90_put_att(model_info%ncid, id_lat, "long_name", "latitude"))
   call handle_error(nf90_put_att(model_info%ncid, id_lat, "units", "degrees_north"))
   call handle_error(nf90_put_att(model_info%ncid, id_lat, "cartesian_axis", "Y"))
   call handle_error(nf90_def_var(model_info%ncid,'time',nf90_double,(/ id_t /), model_info%tid))
   call handle_error(nf90_def_var(model_info%ncid,'area',nf90_float,(/ id_x,id_y/), id_area, &
                                  chunksizes=(/ 200,200/),deflate_level=1,shuffle=.true.))
   call handle_error(nf90_def_var(model_info%ncid,'runoff',nf90_float,(/ id_x,id_y,id_t /), model_info%vid, &
                                  chunksizes=(/ 200,200,1 /),deflate_level=1,shuffle=.true.))
   call handle_error(nf90_copy_att(source_info%ncid,source_info%vid,'units',model_info%ncid,model_info%vid))
   call handle_error(nf90_copy_att(source_info%ncid,source_info%tid,'units',model_info%ncid,model_info%tid))

   call handle_error(nf90_inquire_variable(source_info%ncid,source_info%tid,name=dummy,natts=nats))
   print *,trim(dummy),nats
   do i=1,nats
      call handle_error(nf90_inq_attname(source_info%ncid,source_info%tid,i,dummy))
      print *,trim(dummy),len_trim(dummy)
   enddo
   call handle_error(nf90_copy_att(source_info%ncid,source_info%tid,'calendar',model_info%ncid,model_info%tid))
   call handle_error(nf90_enddef(model_info%ncid,h_minfree=4096))

   call handle_error(nf90_put_var(model_info%ncid,id_lon,lon_model))
   call handle_error(nf90_put_var(model_info%ncid,id_lat,lat_model))
   call handle_error(nf90_put_var(model_info%ncid,id_area,area_model))   



end subroutine setup_model_file



subroutine get_source_info(info)
   type(nc_info_type),intent(inout) :: info
   integer                          :: dids(3), vid
   character(len=32)                :: time_name
   call handle_error(nf90_open(trim(info%fname),nf90_nowrite,info%ncid))
   call handle_error(nf90_inq_varid(info%ncid,'runoff',info%vid))
   call handle_error(nf90_inquire_variable(info%ncid,info%vid,dimids=dids))
   call handle_error(nf90_inquire_dimension(info%ncid,dids(1),len=info%idim))
   call handle_error(nf90_inquire_dimension(info%ncid,dids(2),len=info%jdim))
   call handle_error(nf90_inquire_dimension(info%ncid,dids(3),len=info%nrecs,name=time_name))
   call handle_error(nf90_inq_varid(info%ncid,trim(time_name),info%tid))
   allocate(area_source(info%idim,info%jdim))
   call handle_error(nf90_inq_varid(info%ncid,'area',vid))
   call handle_error(nf90_get_var(info%ncid, vid, area_source ))

end subroutine get_source_info

subroutine get_model_info(info)
   type(nc_info_type),intent(inout) :: info
   integer                          :: ncid,vid,dids(2),i,j,ind
   character(len=128)  :: odir='',ofile=''  ! for ocean_mosaic.nc
   character(len=128)  :: gdir='',gfile=''  ! for hgrid file
   character(len=256)  :: dirfile=''        ! concatenation
   logical             :: fexist
   real(kind=real64),allocatable :: wrk_super(:,:)

   call handle_error(nf90_open('topog.nc',nf90_nowrite,ncid))
   call handle_error(nf90_inq_varid(ncid,'depth',vid))
   call handle_error(nf90_inquire_variable(ncid,vid,dimids=dids))
   call handle_error(nf90_inquire_dimension(ncid,dids(1),len=info%idim))
   call handle_error(nf90_inquire_dimension(ncid,dids(2),len=info%jdim))
   call handle_error(nf90_close(ncid))

   allocate(area_model(info%idim,info%jdim))
   call handle_error(nf90_open('mosaic.nc',nf90_nowrite,ncid))
   call handle_error(nf90_inq_varid(ncid,'ocn_mosaic_dir',vid))
   call handle_error(nf90_get_var(ncid,vid,odir))
   call handle_error(nf90_inq_varid(ncid,'ocn_mosaic_file',vid))
   call handle_error(nf90_get_var(ncid,vid,ofile))
   dirfile=odir(1:scan(odir,'/',back=.true.)) // ofile(1:scan(ofile,'c',back=.true.))
   inquire(file=trim(dirfile),exist=fexist)
   if ( .not. fexist ) then
      write(*,*) 'ocn_mosaic_dir/ocn_mosaic_file =',trim(dirfile), ' does not exist. Bailing out' 
      stop 1
   endif
   call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
   call handle_error(nf90_inq_varid(ncid,'gridlocation',vid))
   call handle_error(nf90_get_var(ncid,vid,gdir))
   call handle_error(nf90_inq_varid(ncid,'gridfiles',vid))
   call handle_error(nf90_get_var(ncid,vid,gfile))
   call handle_error(nf90_close(ncid))
   dirfile=gdir(1:scan(gdir,'/',back=.true.)) // gfile(1:scan(gfile,'c',back=.true.))
   inquire(file=trim(dirfile),exist=fexist)
   if ( .not. fexist ) then
      write(*,*) 'gridlocation/gridfiles =',trim(dirfile), ' does not exist. Bailing out' 
      stop 1
   endif
   call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
   allocate(wrk_super(2*info%idim,2*info%jdim))
   call handle_error(nf90_inq_varid(ncid,'area',vid))
   call handle_error(nf90_get_var(ncid,vid,wrk_super))

   do j = 1, info%jdim
      do i = 1,info%idim
         area_model(i,j) = wrk_super(2*i-1,2*j-1)+wrk_super(2*i,2*j-1)+wrk_super(2*i-1,2*j)+wrk_super(2*i,2*j)
      enddo
   enddo
   deallocate(wrk_super)
   allocate(lon_model(info%idim))
   allocate(lat_model(info%jdim))
   allocate(wrk_super(2*info%idim+1,2*info%jdim+1))
   call handle_error(nf90_inq_varid(ncid,'x',vid))
   call handle_error(nf90_get_var(ncid,vid,wrk_super))
   do i = 1, info%idim
      lon_model(i) = wrk_super(2*i,2)
   enddo
   call handle_error(nf90_inq_varid(ncid,'y',vid))
   call handle_error(nf90_get_var(ncid,vid,wrk_super))
   ind = info%idim/4;
   do j = 1, info%jdim
      lat_model(j) = wrk_super(2*ind+1,2*j)
   enddo
   deallocate(wrk_super)
   call handle_error(nf90_close(ncid))

end subroutine get_model_info

!subroutine read_connection_file_nn(source,coast,coast_nn,qc)

subroutine read_weights_file(source,model)
   type(source_type), intent(out) :: source
   type(model_type), intent(out) :: model
   integer(int32) :: ncid, did, vid

   call handle_error(nf90_open('runoff_weights.nc',nf90_nowrite,ncid))

   call handle_error(nf90_inq_dimid(ncid,'ir',did))
   call handle_error(nf90_inquire_dimension(ncid,did,len=model%npts))

   call handle_error(nf90_inq_dimid(ncid,'is',did))
   call handle_error(nf90_inquire_dimension(ncid,did,len=source%npts))

   allocate(model%i(model%npts),model%j(model%npts),model%weight(model%npts),model%idx(model%npts))
   allocate(source%i(source%npts),source%j(source%npts))

   call handle_error(nf90_inq_varid(ncid,'runoff_i',vid))
   call handle_error(nf90_get_var(ncid,vid,model%i))
   call handle_error(nf90_inq_varid(ncid,'runoff_j',vid))
   call handle_error(nf90_get_var(ncid,vid,model%j))
   call handle_error(nf90_inq_varid(ncid,'runoff_weight',vid))
   call handle_error(nf90_get_var(ncid,vid,model%weight))
   call handle_error(nf90_inq_varid(ncid,'source_index',vid))
   call handle_error(nf90_get_var(ncid,vid,model%idx))


   call handle_error(nf90_inq_varid(ncid,'source_i',vid))
   call handle_error(nf90_get_var(ncid,vid,source%i))
   call handle_error(nf90_inq_varid(ncid,'source_j',vid))
   call handle_error(nf90_get_var(ncid,vid,source%j))
   call handle_error(nf90_close(ncid))

end subroutine read_weights_file

subroutine handle_error(error_flag,isfatal,err_string)
! Simple error handle for NetCDF
integer,intent(in) :: error_flag
logical, intent(in),optional :: isfatal
character(*), intent(in),optional :: err_string
logical            :: fatal
fatal = .true.
if(present(isfatal)) fatal=isfatal
if ( error_flag  /= nf90_noerr ) then
   if ( fatal ) then
      write(*,*) 'FATAL ERROR:',nf90_strerror(error_flag)
      if (present(err_string)) write(*,*) trim(err_string)
      stop
   endif
endif
end subroutine handle_error
    
end program process_runoff
