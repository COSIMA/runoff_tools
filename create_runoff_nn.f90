program create_runoff_nn
use iso_fortran_env
use netcdf
use kdtree2_precision_module
use kdtree2_module
implicit none

type :: coast_type
   integer(kind=int32)                     :: npts
   integer(kind=int32),allocatable,dimension(:) :: i
   integer(kind=int32),allocatable,dimension(:) :: j
   real(kind=kdkind),allocatable,dimension(:)   :: x
   real(kind=kdkind),allocatable,dimension(:)   :: y
   real(kind=kdkind),allocatable,dimension(:)   :: area
end type coast_type
   
type :: source_type
   integer(kind=int32)                     :: npts
   integer(kind=int32),allocatable,dimension(:) :: i
   integer(kind=int32),allocatable,dimension(:) :: j
   real(kind=kdkind),allocatable,dimension(:)   :: x
   real(kind=kdkind),allocatable,dimension(:)   :: y
   real(kind=kdkind),allocatable,dimension(:)   :: area
   real(kind=kdkind),allocatable,dimension(:)   :: dis
   integer(kind=int32),allocatable,dimension(:) :: idx
end type source_type

real(kind=real64),parameter             :: DEG2RAD = asin(1.0_real64)/90.0_real64  ! PI/180
real(kind=real64),parameter             :: REARTH = 6.371d3  !km


real(kind=kdkind), allocatable,dimension(:,:) :: pos_source,pos_coast
type(kdtree2_result) :: results(1)
type(kdtree2),pointer    :: tree

type(coast_type) :: coast
type(source_type) :: source

real(kind=real64)   :: x,y

integer(kind=int32) :: i,j,k
integer(kind=int32) :: nx,ny                  ! Size of model grid
integer(kind=int32) :: nxp,nyp                ! Size of model supergrid
integer(kind=int32) :: ni_source,nj_source    ! Size of runoff input variable

integer(kind=int32) :: ncid,vid,did           ! NetCDF ids
integer(kind=int32) :: dids2d(2)              ! 2D dimension ids array

real(kind=real64),allocatable,dimension(:)     :: x1d_source,y1d_source ! x,y of source (spherical source)
real(kind=real64),allocatable,dimension(:,:)   :: wrk,wrk_super
integer(kind=int16),allocatable,dimension(:,:) :: mask
logical,allocatable,dimension(:,:)       ::wet

character(len=128)  :: odir='',ofile=''  ! for ocean_mosaic.nc
character(len=128)  :: gdir='',gfile=''  ! for hgrid file
character(len=128)  :: tdir='',tfile=''  ! for topgraphy file
character(len=256)  :: dirfile=''        ! concatenation

logical             :: tripolar=.false.,x_cyclic=.false.,south_open=.false.,east_open=.false., &
                       north_open=.false.,west_open=.false. ! BC types

integer(kind=int32)      :: nargs
character(len=12)   :: carg,dimname

! Get info on the grid from input
nargs=command_argument_count()
do i = 1,nargs
   call get_command_argument(i,carg)
   select case(trim(carg))
   case ('tripolar')
        tripolar = .true.
        x_cyclic = .true.
        south_open=.false.
   case ('x_cyclic')
        x_cyclic = .true.
   case ('north_open')
        north_open=.true.
   case ('south_open')
        south_open=.true.
   case ('east_open')
        south_open=.true.
   case ('west_open')
        west_open=.true.
   case DEFAULT
        write(*,*) 'Unknown argument',trim(carg)
        stop
   end select
enddo
if(tripolar .and. south_open) then
  write(*,*) 'Illegal combination: tripolar=.true. and south_open=true.'
  stop
endif

write(*,*) 'Getting model grid info'
! Get mosaic info
call handle_error(nf90_open('mosaic.nc',nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'ocn_mosaic_dir',vid))
call handle_error(nf90_get_var(ncid,vid,odir))
call handle_error(nf90_inq_varid(ncid,'ocn_mosaic_file',vid))
call handle_error(nf90_get_var(ncid,vid,ofile))
call handle_error(nf90_inq_varid(ncid,'ocn_topog_dir',vid))
call handle_error(nf90_get_var(ncid,vid,tdir))
call handle_error(nf90_inq_varid(ncid,'ocn_topog_file',vid))
call handle_error(nf90_get_var(ncid,vid,tfile))
call handle_error(nf90_close(ncid))
! Get horizontal grid
dirfile=odir(1:scan(odir,'/',back=.true.)) // ofile(1:scan(ofile,'c',back=.true.))
write(*,*) len_trim(dirfile),dirfile
call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'gridlocation',vid))
call handle_error(nf90_get_var(ncid,vid,gdir))
call handle_error(nf90_inq_varid(ncid,'gridfiles',vid))
call handle_error(nf90_get_var(ncid,vid,gfile))
call handle_error(nf90_close(ncid))

! Read wet via wrk
write(*,*) 'Creating Coastal mask'

dirfile=tdir(1:scan(tdir,'/',back=.true.)) // tfile(1:scan(tfile,'c',back=.true.))
call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
call handle_error(nf90_inq_dimid(ncid,'nx',did))
call handle_error(nf90_inquire_dimension(ncid,did,len=nx))
call handle_error(nf90_inq_dimid(ncid,'ny',did))
call handle_error(nf90_inquire_dimension(ncid,did,len=ny))
call handle_error(nf90_inq_varid(ncid,'depth',vid))
allocate(wrk(nx,ny))
call handle_error(nf90_get_var(ncid,vid,wrk))
call handle_error(nf90_close(ncid))

! Extend cyclic/tripolar/OB/solid

allocate(wet(0:nx+1,0:ny+1))
do j=1,ny
   do i=1,nx
      wet(i,j) = wrk(i,j)>0.0_real64
   enddo
enddo
if ( tripolar ) then
   do i = 0, nx+1
      wet(i,0) = .false.
   enddo
   do i = 1,nx/2
      wet(i,ny+1) = wet(nx-i+1,ny)
   enddo
   do j = 1, ny+1
      wet(0,j) = wet(nx,j)
      wet(nx+1,j) = wet(1,j)
   enddo
else
   if ( south_open ) then
      do i = 1,nx
         wet(i,0) = wet(i,1)
      enddo
   else
      do i = 1,nx
         wet(i,0) = .false.
      enddo
   endif
   if ( north_open ) then
      do i = 1,nx
         wet(i,ny+1) = wet(i,ny)
      enddo
   else
      do i = 1,nx
         wet(i,ny+1) = .false.
      enddo
   endif
   if ( x_cyclic ) then
      do j = 0, ny+1
         wet(0,j) = wet(nx,j)
         wet(nx+1,j) = wet(1,j)
      enddo
   else
      if ( west_open ) then
         do j = 0, ny+1
            wet(0,j) = wet(1,j)
         enddo
      else
         do j = 0, ny+1
            wet(0,j) = .false.
         enddo
      endif
      if ( east_open ) then
         do j = 0, ny+1
            wet(nx+1,j) = wet(nx,j)
         enddo
      else
         do j = 0, ny+1
            wet(nx+1,j) = .false.
         enddo
      endif
   endif
endif

! Test for coast
k = 0
do j = 1, ny
   do i = 1, nx
      if( wet(i,j)  .and.  .not. all(wet(i-1:i+1,j-1:j+1))) then
         k = k + 1
      endif
   enddo
enddo
coast%npts=k
write(*,*) coast%npts, ' model coastal tiles'
allocate(coast%i(coast%npts),coast%j(coast%npts),coast%x(coast%npts),coast%y(coast%npts),coast%area(coast%npts))
k = 0
do j = 1, ny
   do i = 1, nx
      if( wet(i,j)  .and.  .not. all(wet(i-1:i+1,j-1:j+1))) then
         k = k + 1
         coast%i(k) = i
         coast%j(k) = j
      endif
   enddo
enddo

deallocate(wet)
!
! On mosaic "supergrid" we need to get every second point
!
write(*,*) 'Reading supergrid info'
nxp = 2*nx+1
nyp = 2*ny+1
allocate(wrk_super(nxp,nyp))
! Read xt
dirfile=gdir(1:scan(gdir,'/',back=.true.)) // gfile(1:scan(gfile,'c',back=.true.))
call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
call handle_error(nf90_inq_varid(ncid,'x',vid))
call handle_error(nf90_get_var(ncid,vid,wrk_super))
do j=1,ny
   do i = 1,nx
      wrk(i,j)= wrk_super(2*i,2*j)
   enddo
enddo


do i = 1,coast%npts
   coast%x(i)=wrk(coast%i(i),coast%j(i))
enddo

! Read yt
call handle_error(nf90_inq_varid(ncid,'y',vid))
call handle_error(nf90_get_var(ncid,vid,wrk_super))
do j=1,ny
   do i = 1,nx
      wrk(i,j)= wrk_super(2*i,2*j)
   enddo
enddo
do i = 1,coast%npts
   coast%y(i)=wrk(coast%i(i),coast%j(i))
enddo

! Area, probably should do dxt*dyt correctly but I think this is ok.
deallocate(wrk_super)
allocate(wrk_super(nxp-1,nyp-1))
call handle_error(nf90_inq_varid(ncid,'area',vid))
call handle_error(nf90_get_var(ncid,vid,wrk_super))
call handle_error(nf90_close(ncid))
do j = 1, ny
   do i = 1,nx
      wrk(i,j) = wrk_super(2*i-1,2*j-1)+wrk_super(2*i,2*j-1)+wrk_super(2*i-1,2*j)+wrk_super(2*i,2*j)
   enddo
enddo
do i = 1,coast%npts
   coast%area(i)=wrk(coast%i(i),coast%j(i))
enddo

deallocate(wrk,wrk_super)

write(*,*) 'Creating kdtree'

! Points on surface of sphere radius 1
allocate(pos_coast(3,coast%npts))
do i = 1,coast%npts
   pos_coast(1,i) =  cos(coast%x(i)*DEG2RAD)*cos(coast%y(i)*DEG2RAD)
   pos_coast(2,i) =  sin(coast%x(i)*DEG2RAD)*cos(coast%y(i)*DEG2RAD)
   pos_coast(3,i) =  sin(coast%y(i)*DEG2RAD)
enddo
! Create tree
tree => kdtree2_create(pos_coast,sort=.false.,rearrange=.true.)

! Read wet mask for runoff, MUST have coordinate variables in source.
write(*,*) 'Getting coastal mask for input'
call handle_error(nf90_open('runoff_mask.nc',nf90_nowrite,ncid))

! Get logical mask and find dimensions

call handle_error(nf90_inq_varid(ncid,'mask',vid))
call handle_error(nf90_inquire_variable(ncid,vid,dimids=dids2d))
call handle_error(nf90_inquire_dimension(ncid,dids2d(1),len=ni_source))
call handle_error(nf90_inquire_dimension(ncid,dids2d(2),len=nj_source))
allocate(mask(ni_source,nj_source),x1d_source(ni_source),y1d_source(nj_source),wrk(ni_source,nj_source))

! Get variables and names of coordinate variables

call handle_error(nf90_get_var(ncid, vid, mask))

call handle_error(nf90_inquire_dimension(ncid,dids2d(1),name=dimname))
call handle_error(nf90_inq_varid(ncid,trim(dimname),vid))
call handle_error(nf90_get_var(ncid, vid, x1d_source ))

call handle_error(nf90_inquire_dimension(ncid,dids2d(2),name=dimname))
call handle_error(nf90_inq_varid(ncid,trim(dimname),vid))
call handle_error(nf90_get_var(ncid, vid, y1d_source ))

call handle_error(nf90_inq_varid(ncid,'area',vid))
call handle_error(nf90_get_var(ncid, vid, wrk ))
call handle_error(nf90_close(ncid))


source%npts = count(mask>0_int16)
write(*,*) source%npts,' coastal tiles'
allocate(source%i(source%npts),source%j(source%npts))
allocate(source%x(source%npts),source%y(source%npts))
allocate(source%area(source%npts))
k = 0
do j = 1, nj_source
   do i = 1, ni_source
      if ( mask(i,j)>0_int16 ) then
         k = k + 1
         source%i(k) = i
         source%j(k) = j
      endif
   enddo
enddo

allocate(pos_source(3,source%npts))
do i = 1, source%npts
   x = x1d_source(source%i(i))*DEG2RAD
   y = y1d_source(source%j(i))*DEG2RAD
   source%area(i) = wrk(source%i(i),source%j(i))
   source%x(i) = x1d_source(source%i(i))
   source%y(i) = y1d_source(source%j(i))
   pos_source(1,i) =  cos(x)*cos(y)
   pos_source(2,i) =  sin(x)*cos(y)
   pos_source(3,i) =  sin(y)
enddo
deallocate(wrk,mask,x1d_source,y1d_source)

allocate(source%dis(source%npts),source%idx(source%npts))
! Just get nearest neighbour

write(*,*) 'Finding nearest neighbours'

do i = 1, source%npts
   call kdtree2_n_nearest(tp=tree,qv=pos_source(:,i),nn=1,results=results)
   source%dis(i) = 2.0*asin(sqrt(results(1)%dis)/2.0_real64)*REARTH
   source%idx(i)  = results(1)%idx
enddo
print *,minval(source%dis),maxval(source%dis),DEG2RAD,1./DEG2RAD

! We now have the nearest neighbour to our runoff points.
! Write these points out to a file.
! We also want to apply some sort of QC

write(*,*) 'Writing out connection file'

call create_connection_file(source,coast)
write(*,*) 'Done'

contains

subroutine create_connection_file(source,coast)
   type(coast_type), intent(in) :: coast
   type(source_type), intent(in) :: source
   integer(kind=int32),dimension(:),allocatable :: itarget
   real(kind=real64),dimension(:),allocatable :: rtarget
   integer(kind=int32) :: ncid, did_ic, did_is
   integer(kind=int32) :: coast_i_id,coast_j_id
   integer(kind=int32) :: coast_x_id,coast_y_id
   integer(kind=int32) :: coast_area_id
   integer(kind=int32) :: source_i_id,source_j_id
   integer(kind=int32) :: source_x_id,source_y_id
   integer(kind=int32) :: source_area_id
   integer(kind=int32) :: target_i_id,target_j_id
   integer(kind=int32) :: target_x_id,target_y_id
   integer(kind=int32) :: target_area_id
   integer(kind=int32) :: target_idx_id
   integer(kind=int32) :: dist_id
   integer(kind=int32) :: qc_id


   call handle_error(nf90_create('runoff_connection_nn.nc',ior(NF90_CLOBBER,NF90_NETCDF4),ncid))
   call handle_error(nf90_def_dim(ncid,'ic',coast%npts,did_ic))
   call handle_error(nf90_def_dim(ncid,'is',source%npts,did_is))
!
! These are for ALL the model coastal points
!
   call handle_error(nf90_def_var(ncid,'coast_i',nf90_int,did_ic,coast_i_id))
   call handle_error(nf90_put_att(ncid,coast_i_id,'long_name','model i index'))
   call handle_error(nf90_def_var(ncid,'coast_j',nf90_int,did_ic,coast_j_id))
   call handle_error(nf90_put_att(ncid,coast_j_id,'long_name','model j index'))
   call handle_error(nf90_def_var(ncid,'coast_x',nf90_double,did_ic,coast_x_id))
   call handle_error(nf90_put_att(ncid,coast_x_id,'long_name','model longitude'))
   call handle_error(nf90_put_att(ncid,coast_x_id,'units','degrees_E'))
   call handle_error(nf90_def_var(ncid,'coast_y',nf90_double,did_ic,coast_y_id))
   call handle_error(nf90_put_att(ncid,coast_y_id,'long_name','model latitude'))
   call handle_error(nf90_put_att(ncid,coast_y_id,'units','degrees_N'))
   call handle_error(nf90_def_var(ncid,'coast_area',nf90_double,did_ic,coast_area_id))
   call handle_error(nf90_put_att(ncid,coast_area_id,'long_name','model area'))
   call handle_error(nf90_put_att(ncid,coast_area_id,'units','m^2'))

!
! Source
!

   call handle_error(nf90_def_var(ncid,'source_i',nf90_int,did_is,source_i_id))
   call handle_error(nf90_put_att(ncid,source_i_id,'long_name','runoff source i index'))
   call handle_error(nf90_def_var(ncid,'source_j',nf90_int,did_is,source_j_id))
   call handle_error(nf90_put_att(ncid,source_j_id,'long_name','runoff source j index'))
   call handle_error(nf90_def_var(ncid,'source_x',nf90_double,did_is,source_x_id))
   call handle_error(nf90_put_att(ncid,source_x_id,'long_name','runoff source longitude'))
   call handle_error(nf90_put_att(ncid,source_x_id,'units','degrees_E'))
   call handle_error(nf90_def_var(ncid,'source_y',nf90_double,did_is,source_y_id))
   call handle_error(nf90_put_att(ncid,source_y_id,'long_name','runoff source latitude'))
   call handle_error(nf90_put_att(ncid,source_y_id,'units','degrees_N'))
   call handle_error(nf90_def_var(ncid,'source_area',nf90_double,did_is,source_area_id))
   call handle_error(nf90_put_att(ncid,source_area_id,'long_name','runoff source area'))
   call handle_error(nf90_put_att(ncid,source_area_id,'units','m^2'))
!
! Define nearest neighbour  subset
!
   call handle_error(nf90_def_var(ncid,'target_i',nf90_int,did_is,target_i_id))
   call handle_error(nf90_put_att(ncid,target_i_id,'long_name','model target i index'))
   call handle_error(nf90_def_var(ncid,'target_j',nf90_int,did_is,target_j_id))
   call handle_error(nf90_put_att(ncid,target_j_id,'long_name','model target j index'))
   call handle_error(nf90_def_var(ncid,'target_x',nf90_double,did_is,target_x_id))
   call handle_error(nf90_put_att(ncid,target_x_id,'long_name','model target longitude'))
   call handle_error(nf90_put_att(ncid,target_x_id,'units','degrees_E'))
   call handle_error(nf90_def_var(ncid,'target_y',nf90_double,did_is,target_y_id))
   call handle_error(nf90_put_att(ncid,target_y_id,'long_name','model target latitude'))
   call handle_error(nf90_put_att(ncid,target_y_id,'units','degrees_N'))
   call handle_error(nf90_def_var(ncid,'target_area',nf90_double,did_is,target_area_id))
   call handle_error(nf90_put_att(ncid,target_area_id,'long_name','model target area'))
   call handle_error(nf90_put_att(ncid,target_area_id,'units','m^2'))
   call handle_error(nf90_def_var(ncid,'target_idx',nf90_int,did_is,target_idx_id))
   call handle_error(nf90_put_att(ncid,target_idx_id,'long_name','packed target index'))
   call handle_error(nf90_def_var(ncid,'dist',nf90_double,did_is,dist_id))
   call handle_error(nf90_put_att(ncid,dist_id,'long_name','distance to target'))
   call handle_error(nf90_put_att(ncid,dist_id,'units','km'))
   call handle_error(nf90_def_var(ncid,'source_qc',nf90_int,did_is,qc_id))
   call handle_error(nf90_put_att(ncid,qc_id,'long_name','QC flag'))
   call handle_error(nf90_enddef(ncid,h_minfree=4096))

! Put it there
   call handle_error(nf90_put_var(ncid,coast_i_id,coast%i))
   call handle_error(nf90_put_var(ncid,coast_j_id,coast%j))
   call handle_error(nf90_put_var(ncid,coast_x_id,coast%x))
   call handle_error(nf90_put_var(ncid,coast_y_id,coast%y))
   call handle_error(nf90_put_var(ncid,coast_area_id,coast%area))

   call handle_error(nf90_put_var(ncid,source_i_id,source%i))
   call handle_error(nf90_put_var(ncid,source_j_id,source%j))
   call handle_error(nf90_put_var(ncid,source_x_id,source%x))
   call handle_error(nf90_put_var(ncid,source_y_id,source%y))
   call handle_error(nf90_put_var(ncid,source_area_id,source%area))
!
! Target arrays
!
   allocate(itarget(source%npts),rtarget(source%npts))

   itarget = coast%i(source%idx(:))
   call handle_error(nf90_put_var(ncid,target_i_id,itarget))
   itarget = coast%j(source%idx(:))
   call handle_error(nf90_put_var(ncid,target_j_id,itarget))
   rtarget = coast%x(source%idx(:))
   call handle_error(nf90_put_var(ncid,target_x_id,rtarget))
   rtarget = coast%y(source%idx(:))
   call handle_error(nf90_put_var(ncid,target_y_id,rtarget))
   rtarget = coast%area(source%idx(:))
   call handle_error(nf90_put_var(ncid,target_area_id,rtarget))
   call handle_error(nf90_put_var(ncid,target_idx_id,source%idx))
   call handle_error(nf90_put_var(ncid,dist_id,source%dis))
   deallocate(rtarget)
   itarget = 0
! Simple QC. If further than 100km then flag.
   do i = 1, source%npts
      if ( source%dis(i) > 100.0_real64 ) itarget(i) = 1
   enddo

   
   call handle_error(nf90_put_var(ncid,qc_id,itarget))
   call handle_error(nf90_close(ncid))

end subroutine create_connection_file

subroutine handle_error(error_flag,isfatal,err_string)
! Simple error handle for NetCDF
integer(kind=int32),intent(in) :: error_flag
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
    
end program create_runoff_nn
