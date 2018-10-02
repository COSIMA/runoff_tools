program create_model_coast
  !
  ! Create a file with all the candidate coast points that we wish to map to.
  ! At the moment this program only finds points immediately adjacent to land
  !
  ! Mosaic input.
  !
  ! Usage: create_model_coast boundary_types
  !
  ! boundary_types may be
  !     tripolar  (sets x_cyclic and solid southern boundary)
  ! OR
  !     a combination north_open, south_open, east_open, west_open. No argument
  !     means this boundary is solid
  !AND
  !     x_cyclic for periodic boundaries.
  !
  ! OUtput: A file called model_coast.nc
  !
  !
  use iso_fortran_env
  use netcdf
  use runoff_modules
  implicit none

  type :: att_type
     character(len=64)                            :: command
     character(len=256)                           :: mosaic_file
     character(len=256)                           :: grid_file
     character(len=256)                           :: topog_file
     character(len=24)                            :: creation_time
     character(len=24)                            :: creator = 'Dr V. Frankenstein'
  end type att_type

  type(coast_type) :: coast
  type(att_type)   :: attributes

  integer(kind=int32) :: i,j,k
  integer(kind=int32) :: nx,ny                  ! Size of model grid
  integer(kind=int32) :: nxp,nyp                ! Size of model supergrid

  integer(kind=int32) :: ncid,vid,did           ! NetCDF ids
  integer(kind=int32) :: istat                  ! 

  real(kind=real64),allocatable,dimension(:,:)   :: wrk,wrk_super
  logical,allocatable,dimension(:,:)       ::wet

  character(len=128)  :: odir='',ofile=''  ! for ocean_mosaic.nc
  character(len=128)  :: gdir='',gfile=''  ! for hgrid file
  character(len=128)  :: tdir='',tfile=''  ! for topgraphy file
  character(len=256)  :: dirfile=''        ! concatenation

  logical             :: fexist = .false.

  logical             :: tripolar=.false.,x_cyclic=.false.,south_open=.false.,east_open=.false., &
       north_open=.false.,west_open=.false. ! BC types

  integer(kind=int32)      :: nargs
  character(len=12)   :: carg


  call get_command(attributes%command)
  call date_and_time(carg)
  attributes%creation_time=carg(1:4) // '-' // carg(5:6) // '-' // carg(7:8)
  call date_and_time(time=carg)
  attributes%creation_time=trim(attributes%creation_time) // ' ' // carg(1:2) // ':' // carg(3:4) // ':' // carg(5:8)
  call get_environment_variable('USER',value=attributes%creator,status=istat)
  if(istat /= 0) attributes%creator = 'A mysterious entity'
  ! Get info on the grid from input
  nargs=command_argument_count()
  if ( nargs == 0 ) then
     write(*,*) 'No arguments, all solid boundaries'
  else
     do i = 1,nargs
        call get_command_argument(i,carg)
        select case(trim(carg))
        case ('tripolar')
           tripolar = .true.
           x_cyclic = .true.
           south_open=.false.
           write(*,*) 'tripolar'
        case ('x_cyclic')
           x_cyclic = .true.
           write(*,*) 'cyclic x boundaries'
        case ('north_open')
           north_open=.true.
           write(*,*) 'Open boundary to the North'
        case ('south_open')
           south_open=.true.
           write(*,*) 'Open boundary to the South'
        case ('east_open')
           south_open=.true.
           write(*,*) 'Open boundary to the East'
        case ('west_open')
           west_open=.true.
           write(*,*) 'Open boundary to the West'
        case DEFAULT
           write(*,*) 'Unknown argument',trim(carg)
           stop
        end select
     enddo
     if(tripolar .and. south_open) then
        write(*,*) 'Illegal combination: tripolar=.true. and south_open=true.'
        stop
     endif
  endif

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
  attributes%mosaic_file=dirfile
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
  write(*,*) 'Creating Coastal mask'

  dirfile=tdir(1:scan(tdir,'/',back=.true.)) // tfile(1:scan(tfile,'c',back=.true.))
  attributes%topog_file=trim(dirfile)
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
     do i = 1,nx
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
  attributes%grid_file=trim(dirfile)
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

  ! Write out coastal mask file

  write(*,*) 'Writing out connection file'

  call create_coast_file(coast,attributes)
  write(*,*) 'Done'

contains

  subroutine create_coast_file(coast,attributes)
    type(coast_type), intent(in) :: coast
    type(att_type), intent(in) ::  attributes
    integer(kind=int32) :: ncid, did_ic
    integer(kind=int32) :: coast_i_id,coast_j_id
    integer(kind=int32) :: coast_x_id,coast_y_id
    integer(kind=int32) :: coast_area_id

    call handle_error(nf90_create('model_coast.nc',ior(NF90_CLOBBER,NF90_NETCDF4),ncid))
    call handle_error(nf90_def_dim(ncid,'ic',coast%npts,did_ic))
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

    call handle_error(nf90_put_att(ncid,nf90_global,'command',trim(attributes%command)))
    call handle_error(nf90_put_att(ncid,nf90_global,'mosaic_file',trim(attributes%mosaic_file)))
    call handle_error(nf90_put_att(ncid,nf90_global,'grid_file',trim(attributes%grid_file)))
    call handle_error(nf90_put_att(ncid,nf90_global,'topog_file',trim(attributes%topog_file)))
    call handle_error(nf90_put_att(ncid,nf90_global,'creation_time',trim(attributes%creation_time)))
    call handle_error(nf90_put_att(ncid,nf90_global,'creator',trim(attributes%creator)))
    call handle_error(nf90_enddef(ncid,h_minfree=4096))

    ! Put it there
    call handle_error(nf90_put_var(ncid,coast_i_id,coast%i))
    call handle_error(nf90_put_var(ncid,coast_j_id,coast%j))
    call handle_error(nf90_put_var(ncid,coast_x_id,coast%x))
    call handle_error(nf90_put_var(ncid,coast_y_id,coast%y))
    call handle_error(nf90_put_var(ncid,coast_area_id,coast%area))

    call handle_error(nf90_close(ncid))

  end subroutine create_coast_file

end program create_model_coast
