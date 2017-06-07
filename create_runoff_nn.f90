program create_runoff_nn

  use iso_fortran_env
  use netcdf
  use kdtree2_precision_module
  use kdtree2_module
  use runoff_modules

  implicit none

  real(kind=real64),parameter             :: DEG2RAD = asin(1.0_real64)/90.0_real64  ! PI/180
  real(kind=real64),parameter             :: REARTH = 6.371d3  !km


  real(kind=kdkind), allocatable,dimension(:,:) :: pos_source,pos_coast
  type(kdtree2_result) :: results(1)
  type(kdtree2),pointer    :: tree

  type(coast_type) :: coast
  type(source_type) :: source

  real(kind=real64)   :: x,y

  integer(kind=int32) :: i,j,k
  integer(kind=int32) :: ni_source,nj_source    ! Size of runoff input variable

  integer(kind=int32) :: ncid,vid           ! NetCDF ids
  integer(kind=int32) :: dids2d(2)              ! 2D dimension ids array

  real(kind=real64),allocatable,dimension(:)     :: x1d_source,y1d_source ! x,y of source (spherical source)
  real(kind=real64),allocatable,dimension(:,:)   :: wrk
  integer(kind=int16),allocatable,dimension(:,:) :: mask

  character(len=128)  :: coastfile

  integer(kind=int32)      :: nargs
  character(len=32)   :: carg,dimname

  ! Get info on the grid from input
  nargs=command_argument_count()
  do i = 1,nargs
     call get_command_argument(i,carg)
     select case(trim(carg))
     case ('-f')
        call get_command_argument(i+1,carg)
        coastfile = trim(carg)
        exit
     case DEFAULT
        write(*,*) 'Unknown argument',trim(carg)
        stop
     end select
  enddo

  call read_coast_mask(coast,coastfile)

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

  ! Note that results%dis is the square of the distance. Need to take sqare root
  ! and convert chord to distance over great circle.
  do i = 1, source%npts
     call kdtree2_n_nearest(tp=tree,qv=pos_source(:,i),nn=1,results=results)
     source%dis(i) = 2.0*asin(sqrt(results(1)%dis)/2.0_real64)*REARTH
     source%idx(i)  = results(1)%idx
  enddo

  ! We now have the nearest neighbour to our runoff points.
  ! Write these points out to a file.
  ! We may also want to apply some sort of QC

  write(*,*) 'Writing out connection file'

  call create_connection_file(source,coast)
  write(*,*) 'Done'

end program create_runoff_nn
