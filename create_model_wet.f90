program create_model_wet
  !
  ! Create a file with all the wet points that we wish to map to.
  !
  ! Mosaic input.
  !
  ! Usage: create_model_wet
  !
  ! Assumes there's a file 'mosaic.nc' lurking about....
  !
  ! Output. A file called model_wet.nc
  !
  !
  use iso_fortran_env
  use netcdf
  implicit none

  type :: wet_type
     integer(kind=int32)                     :: npts
     integer(kind=int32),allocatable,dimension(:) :: i
     integer(kind=int32),allocatable,dimension(:) :: j
     real(kind=real64),allocatable,dimension(:)   :: x
     real(kind=real64),allocatable,dimension(:)   :: y
     real(kind=real64),allocatable,dimension(:)   :: area
  end type wet_type

  type(wet_type) :: wet

  integer(kind=int32) :: i,j,k
  integer(kind=int32) :: nx,ny                  ! Size of model grid
  integer(kind=int32) :: nxp,nyp                ! Size of model supergrid

  integer(kind=int32) :: ncid,vid,did           ! NetCDF ids

  real(kind=real64),allocatable,dimension(:,:)   :: wrk,wrk_super
  logical,allocatable,dimension(:,:)       ::iswet

  character(len=128)  :: odir='',ofile=''  ! for ocean_mosaic.nc
  character(len=128)  :: gdir='',gfile=''  ! for hgrid file
  character(len=128)  :: tdir='',tfile=''  ! for topgraphy file
  character(len=256)  :: dirfile=''        ! concatenation

  logical             :: fexist = .false.

  ! Get info on the grid from input

  write(*,*) 'Getting model grid info'
  ! Get mosaic info
  inquire(file=trim('mosaic.nc'),exist=fexist)
  if ( .not. fexist ) then
     write(*,*) 'mosaic.nc does not exist. Bailing out' 
     stop 1
  endif
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

  ! Read wet via wrk
  write(*,*) 'Water mask'

  dirfile=tdir(1:scan(tdir,'/',back=.true.)) // tfile(1:scan(tfile,'c',back=.true.))
  inquire(file=trim(dirfile),exist=fexist)
  if ( .not. fexist ) then
     write(*,*) 'ocn_topog_dir/ocn_topog_file =',trim(dirfile), ' does not exist. Bailing out' 
     stop 1
  endif
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

  allocate(iswet(1:nx,ny+1))
  do j=1,ny
     do i=1,nx
        iswet(i,j) = wrk(i,j)>0.0_real64
     enddo
  enddo

  ! Test for coast
  k = 0
  do j = 1, ny
     do i = 1, nx
        if( iswet(i,j) ) then
           k = k + 1
        endif
     enddo
  enddo
  wet%npts=k
  write(*,*) wet%npts, ' model wet tiles'
  allocate(wet%i(wet%npts),wet%j(wet%npts),wet%x(wet%npts),wet%y(wet%npts),wet%area(wet%npts))
  k = 0
  do j = 1, ny
     do i = 1, nx
        if( iswet(i,j)) then
           k = k + 1
           wet%i(k) = i
           wet%j(k) = j
        endif
     enddo
  enddo

  deallocate(iswet)
  !
  ! On mosaic "supergrid" we need to get every second point
  !
  write(*,*) 'Reading supergrid info'
  nxp = 2*nx+1
  nyp = 2*ny+1
  allocate(wrk_super(nxp,nyp))
  ! Read xt
  dirfile=gdir(1:scan(gdir,'/',back=.true.)) // gfile(1:scan(gfile,'c',back=.true.))
  inquire(file=trim(dirfile),exist=fexist)
  if ( .not. fexist ) then
     write(*,*) 'gridlocation/gridfiles =',trim(dirfile), ' does not exist. Bailing out' 
     stop 1
  endif
  call handle_error(nf90_open(trim(dirfile),nf90_nowrite,ncid))
  call handle_error(nf90_inq_varid(ncid,'x',vid))
  call handle_error(nf90_get_var(ncid,vid,wrk_super))
  do j=1,ny
     do i = 1,nx
        wrk(i,j)= wrk_super(2*i,2*j)
     enddo
  enddo


  do i = 1,wet%npts
     wet%x(i)=wrk(wet%i(i),wet%j(i))
  enddo

  ! Read yt
  call handle_error(nf90_inq_varid(ncid,'y',vid))
  call handle_error(nf90_get_var(ncid,vid,wrk_super))
  do j=1,ny
     do i = 1,nx
        wrk(i,j)= wrk_super(2*i,2*j)
     enddo
  enddo
  do i = 1,wet%npts
     wet%y(i)=wrk(wet%i(i),wet%j(i))
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
  do i = 1,wet%npts
     wet%area(i)=wrk(wet%i(i),wet%j(i))
  enddo

  deallocate(wrk,wrk_super)

  ! Write out coastal mask file

  write(*,*) 'Writing out connection file'

  call create_wet_file(wet)
  write(*,*) 'Done'

contains

  subroutine create_wet_file(wet)
    type(wet_type), intent(in) :: wet
    integer(kind=int32) :: ncid, did_ic
    integer(kind=int32) :: wet_i_id,wet_j_id
    integer(kind=int32) :: wet_x_id,wet_y_id
    integer(kind=int32) :: wet_area_id

    call handle_error(nf90_create('model_wet.nc',ior(NF90_CLOBBER,NF90_NETCDF4),ncid))
    call handle_error(nf90_def_dim(ncid,'iw',wet%npts,did_ic))
    !
    ! These are for ALL the model wetal points
    !
    call handle_error(nf90_def_var(ncid,'wet_i',nf90_int,did_ic,wet_i_id))
    call handle_error(nf90_put_att(ncid,wet_i_id,'long_name','model i index'))
    call handle_error(nf90_def_var(ncid,'wet_j',nf90_int,did_ic,wet_j_id))
    call handle_error(nf90_put_att(ncid,wet_j_id,'long_name','model j index'))
    call handle_error(nf90_def_var(ncid,'wet_x',nf90_double,did_ic,wet_x_id))
    call handle_error(nf90_put_att(ncid,wet_x_id,'long_name','model longitude'))
    call handle_error(nf90_put_att(ncid,wet_x_id,'units','degrees_E'))
    call handle_error(nf90_def_var(ncid,'wet_y',nf90_double,did_ic,wet_y_id))
    call handle_error(nf90_put_att(ncid,wet_y_id,'long_name','model latitude'))
    call handle_error(nf90_put_att(ncid,wet_y_id,'units','degrees_N'))
    call handle_error(nf90_def_var(ncid,'wet_area',nf90_double,did_ic,wet_area_id))
    call handle_error(nf90_put_att(ncid,wet_area_id,'long_name','model area'))
    call handle_error(nf90_put_att(ncid,wet_area_id,'units','m^2'))

    call handle_error(nf90_enddef(ncid,h_minfree=4096))

    ! Put it there
    call handle_error(nf90_put_var(ncid,wet_i_id,wet%i))
    call handle_error(nf90_put_var(ncid,wet_j_id,wet%j))
    call handle_error(nf90_put_var(ncid,wet_x_id,wet%x))
    call handle_error(nf90_put_var(ncid,wet_y_id,wet%y))
    call handle_error(nf90_put_var(ncid,wet_area_id,wet%area))

    call handle_error(nf90_close(ncid))

  end subroutine create_wet_file

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

end program create_model_wet
