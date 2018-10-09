module runoff_modules

  use iso_fortran_env
  use netcdf
  use kdtree2_precision_module
  use kdtree2_module

  real(kind=real64),parameter             :: DEG2RAD = asin(1.0_real64)/90.0_real64  ! PI/180
  real(kind=real64),parameter             :: REARTH = 6.371d3  !km

  type :: runoff_type
     integer(int32),allocatable,dimension(:) :: idx  !index of source array
     integer(int32),allocatable,dimension(:) :: i    !i index of 2d model grid
     integer(int32),allocatable,dimension(:) :: j    !j index of model grid
     real(kdkind),allocatable,dimension(:)   :: x    !x of model grid
     real(kdkind),allocatable,dimension(:)   :: y    !y of model grid
     real(kdkind),allocatable,dimension(:)   :: weight    ! factor to multiply source  runoff
  end type runoff_type

  type :: coast_type
     integer(kind=int32)                     :: npts
     integer(kind=int32),allocatable,dimension(:) :: i
     integer(kind=int32),allocatable,dimension(:) :: j
     real(kind=real64),allocatable,dimension(:)   :: x
     real(kind=real64),allocatable,dimension(:)   :: y
     real(kind=real64),allocatable,dimension(:)   :: area
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

  type :: wet_type
     integer(kind=int32)                     :: npts
     integer(kind=int32),allocatable,dimension(:) :: i
     integer(kind=int32),allocatable,dimension(:) :: j
     real(kind=real64),allocatable,dimension(:)   :: x
     real(kind=real64),allocatable,dimension(:)   :: y
     real(kind=real64),allocatable,dimension(:)   :: area
  end type wet_type

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

contains

  subroutine read_coast_mask(coast,coastfile)
    type(coast_type), intent(out) :: coast
    character(*)                  :: coastfile
    integer(kind=int32) :: ncid, did_ic,vid
    call handle_error(nf90_open(trim(coastfile),nf90_nowrite,ncid))
    call handle_error(nf90_inq_dimid(ncid,'ic',did_ic))
    call handle_error(nf90_inquire_dimension(ncid,did_ic,len=coast%npts))
    allocate(coast%i(coast%npts),coast%j(coast%npts),coast%x(coast%npts),coast%y(coast%npts),coast%area(coast%npts))
    call handle_error(nf90_inq_varid(ncid,'coast_i',vid))
    call handle_error(nf90_get_var(ncid,vid,coast%i))
    call handle_error(nf90_inq_varid(ncid,'coast_j',vid))
    call handle_error(nf90_get_var(ncid,vid,coast%j))
    call handle_error(nf90_inq_varid(ncid,'coast_x',vid))
    call handle_error(nf90_get_var(ncid,vid,coast%x))
    call handle_error(nf90_inq_varid(ncid,'coast_y',vid))
    call handle_error(nf90_get_var(ncid,vid,coast%y))
    call handle_error(nf90_inq_varid(ncid,'coast_area',vid))
    call handle_error(nf90_get_var(ncid,vid,coast%area))
    call handle_error(nf90_close(ncid))
  end subroutine  read_coast_mask

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
    integer(kind=int32) :: qc_id, i


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

  subroutine get_source_info(info,var_name)
    type(nc_info_type),intent(inout) :: info
    character(len=32),intent(in)     :: var_name
    integer                          :: dids(3)
    character(len=32)                :: time_name
    call handle_error(nf90_open(trim(info%fname),nf90_nowrite,info%ncid),err_string='opening '//info%fname)
    call handle_error(nf90_inq_varid(info%ncid,var_name,info%vid))
    call handle_error(nf90_inquire_variable(info%ncid,info%vid,dimids=dids))
    call handle_error(nf90_inquire_dimension(info%ncid,dids(1),len=info%idim))
    call handle_error(nf90_inquire_dimension(info%ncid,dids(2),len=info%jdim))
    call handle_error(nf90_inquire_dimension(info%ncid,dids(3),len=info%nrecs,name=time_name))
    call handle_error(nf90_inq_varid(info%ncid,trim(time_name),info%tid))
  end subroutine get_source_info

  subroutine get_model_info(info)
    type(nc_info_type),intent(inout) :: info
    integer                          :: ncid,vid,dids(2)
    call handle_error(nf90_open('topog.nc',nf90_nowrite,ncid))
    call handle_error(nf90_inq_varid(ncid,'depth',vid))
    call handle_error(nf90_inquire_variable(ncid,vid,dimids=dids))
    call handle_error(nf90_inquire_dimension(ncid,dids(1),len=info%idim))
    call handle_error(nf90_inquire_dimension(ncid,dids(2),len=info%jdim))
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

  !subroutine read_connection_file_nn(source,coast,coast_nn,qc)
  subroutine read_connection_file_nn(source,coast,coast_nn)
    type(coast_type), intent(out) :: coast
    type(coast_type), intent(out) :: coast_nn
    type(source_type), intent(out) :: source
    !   integer(kind=int32),allocatable,intent(out) :: qc
    integer(int32) :: ic,is
    integer(int32) :: ncid, did_ic, did_is
    integer(int32) :: coast_i_id,coast_j_id
    integer(int32) :: coast_x_id,coast_y_id
    integer(int32) :: coast_area_id
    integer(int32) :: source_i_id,source_j_id
    integer(int32) :: source_x_id,source_y_id
    integer(int32) :: source_area_id
    integer(int32) :: target_i_id,target_j_id
    integer(int32) :: target_x_id,target_y_id
    integer(int32) :: target_area_id
    integer(int32) :: target_idx_id
    integer(int32) :: dist_id
    !   integer(int32) :: qc_id


    call handle_error(nf90_open('runoff_connection_nn.nc',NF90_NOWRITE,ncid))
    call handle_error(nf90_inq_dimid(ncid,'ic',did_ic))
    call handle_error(nf90_inquire_dimension(ncid,did_ic,len=ic))
    call handle_error(nf90_inq_dimid(ncid,'is',did_is))
    call handle_error(nf90_inquire_dimension(ncid,did_is,len=is))
    !
    ! Do allocations
    !
    allocate(coast%i(ic),coast%j(ic),coast%x(ic),coast%y(ic),coast%area(ic))
    allocate(coast_nn%i(is),coast_nn%j(is),coast_nn%x(is),coast_nn%y(is),coast_nn%area(is))
    allocate(source%i(is),source%j(is),source%x(is),source%y(is),source%area(is))
    allocate(source%dis(is),source%idx(is))
    !   allocate(qc(is))

    call handle_error(nf90_inq_varid(ncid,'coast_i',coast_i_id))
    call handle_error(nf90_inq_varid(ncid,'coast_j',coast_j_id))
    call handle_error(nf90_inq_varid(ncid,'coast_x',coast_x_id))
    call handle_error(nf90_inq_varid(ncid,'coast_y',coast_y_id))
    call handle_error(nf90_inq_varid(ncid,'coast_area',coast_area_id))

    call handle_error(nf90_inq_varid(ncid,'source_i',source_i_id))
    call handle_error(nf90_inq_varid(ncid,'source_j',source_j_id))
    call handle_error(nf90_inq_varid(ncid,'source_x',source_x_id))
    call handle_error(nf90_inq_varid(ncid,'source_y',source_y_id))
    call handle_error(nf90_inq_varid(ncid,'source_area',source_area_id))

    call handle_error(nf90_inq_varid(ncid,'target_i',target_i_id))
    call handle_error(nf90_inq_varid(ncid,'target_j',target_j_id))
    call handle_error(nf90_inq_varid(ncid,'target_x',target_x_id))
    call handle_error(nf90_inq_varid(ncid,'target_y',target_y_id))
    call handle_error(nf90_inq_varid(ncid,'target_area',target_area_id))
    call handle_error(nf90_inq_varid(ncid,'target_idx',target_idx_id))
    call handle_error(nf90_inq_varid(ncid,'dist',dist_id))
    !   call handle_error(nf90_inq_varid(ncid,'source_qc',qc_id))

    ! Get it here
    call handle_error(nf90_get_var(ncid,coast_i_id,coast%i))
    call handle_error(nf90_get_var(ncid,coast_j_id,coast%j))
    call handle_error(nf90_get_var(ncid,coast_x_id,coast%x))
    call handle_error(nf90_get_var(ncid,coast_y_id,coast%y))
    call handle_error(nf90_get_var(ncid,coast_area_id,coast%area))

    call handle_error(nf90_get_var(ncid,source_i_id,source%i))
    call handle_error(nf90_get_var(ncid,source_j_id,source%j))
    call handle_error(nf90_get_var(ncid,source_x_id,source%x))
    call handle_error(nf90_get_var(ncid,source_y_id,source%y))
    call handle_error(nf90_get_var(ncid,source_area_id,source%area))
    !
    ! Target arrays
    !

    call handle_error(nf90_get_var(ncid,target_i_id,coast_nn%i))
    call handle_error(nf90_get_var(ncid,target_j_id,coast_nn%j))
    call handle_error(nf90_get_var(ncid,target_x_id,coast_nn%x))
    call handle_error(nf90_get_var(ncid,target_y_id,coast_nn%y))
    call handle_error(nf90_get_var(ncid,target_area_id,coast_nn%area))
    call handle_error(nf90_get_var(ncid,target_idx_id,source%idx))
    call handle_error(nf90_get_var(ncid,dist_id,source%dis))

    !   call handle_error(nf90_get_var(ncid,qc_id,qc))
    call handle_error(nf90_close(ncid))

  end subroutine read_connection_file_nn

  subroutine setup_model_file(source_info,model_info)
    type(nc_info_type),intent(inout) :: source_info
    type(nc_info_type),intent(inout) :: model_info
    integer(kind=int32)              :: id_x,id_y,id_t
    character(len=32)::dummy
    integer:: nats,i
    call handle_error(nf90_create(trim(model_info%fname),ior(nf90_netcdf4,nf90_clobber),model_info%ncid))
    call handle_error(nf90_def_dim(model_info%ncid,'lon',model_info%idim,id_x))
    call handle_error(nf90_def_dim(model_info%ncid,'lat',model_info%jdim,id_y))
    call handle_error(nf90_def_dim(model_info%ncid,'time',nf90_unlimited,id_t))
    call handle_error(nf90_def_var(model_info%ncid,'time',nf90_double,(/ id_t /), model_info%tid))
    call handle_error(nf90_def_var(model_info%ncid,'runoff',nf90_float,(/ id_x,id_y,id_t /), model_info%vid, &
         chunksizes=(/ 200,200,1 /),deflate_level=1,shuffle=.true.))

    ! Copy time attributes
    call handle_error(nf90_inquire_variable(source_info%ncid,source_info%tid,name=dummy,natts=nats))
    print *,trim(dummy),nats
    do i=1,nats
       call handle_error(nf90_inq_attname(source_info%ncid,source_info%tid,i,dummy))
       call handle_error(nf90_copy_att(source_info%ncid,source_info%tid,dummy,model_info%ncid,model_info%tid))
       print *,trim(dummy),len_trim(dummy)
    enddo
    ! Copy variable attributes
    call handle_error(nf90_inquire_variable(source_info%ncid,source_info%vid,name=dummy,natts=nats))
    do i=1,nats
       call handle_error(nf90_inq_attname(source_info%ncid,source_info%vid,i,dummy))
       call handle_error(nf90_copy_att(source_info%ncid,source_info%vid,dummy,model_info%ncid,model_info%vid))
       print *,trim(dummy),len_trim(dummy)
    enddo
    ! Copy global attributes
    call handle_error(nf90_inquire(source_info%ncid,nAttributes=nats))
    do i=1,nats
       call handle_error(nf90_inq_attname(source_info%ncid,NF90_GLOBAL,i,dummy))
       call handle_error(nf90_copy_att(source_info%ncid,NF90_GLOBAL,dummy,model_info%ncid,NF90_GLOBAL))
       print *,trim(dummy),len_trim(dummy)
    enddo

    call handle_error(nf90_enddef(model_info%ncid,h_minfree=4096))
  end subroutine setup_model_file

  subroutine write_mask_file(outfile, mask)

    ! Interface variables
    logical, dimension(:,:), intent(in) :: mask
    character(len=*), intent(in)  :: outfile

    ! Local variables
    integer :: ncid, vid
    integer :: id_x,id_y

    call handle_error(nf90_create(outfile,ior(nf90_netcdf4,nf90_clobber),ncid))
    call handle_error(nf90_def_dim(ncid,'longitude',size(mask,1),id_x))
    call handle_error(nf90_def_dim(ncid,'latitude',size(mask,2),id_y))
    call handle_error(nf90_def_var(ncid,'mask',nf90_short,(/ id_x,id_y /), vid, &
         deflate_level=4,shuffle=.true.))
    call handle_error(nf90_enddef(ncid))

    call handle_error(nf90_put_var(ncid,vid,merge(1,0,mask)))

    call handle_error(nf90_close(ncid))

  end subroutine write_mask_file

  function read_record(source_info, irecord) result(record)

    ! Interface variables
    type(nc_info_type),intent(in) :: source_info
    integer, intent(in)           :: irecord

    ! Return value
    real(kind=real32), dimension(source_info%idim,source_info%jdim) :: record

    call handle_error(nf90_get_var(source_info%ncid,source_info%vid,record,start=(/1,1,irecord/), &
         count=(/source_info%idim,source_info%jdim,1/)))

  end function read_record

end module runoff_modules
