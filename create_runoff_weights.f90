program create_runoff_weights
! Take the nearest neighbours file and create a file with multiple targets.
! i.e. spread the runoff according to a cunning plan!
!
! We can either 
!
!   Distribute to N neighbours of the nearest neighbour
!   Distribute to  radius of the nearest neighbour

!   Weight evenly
!   Weight inverse distance from the source
!   Weight inverse distance fron the nearest neighbour
!
!
use iso_fortran_env
use netcdf
use kdtree2_precision_module
use kdtree2_module
implicit none

type :: runoff_type
   integer(int32),allocatable,dimension(:) :: idx  !index of source array
   integer(int32),allocatable,dimension(:) :: i    !i index of 2d model grid
   integer(int32),allocatable,dimension(:) :: j    !j index of model grid
   real(kdkind),allocatable,dimension(:)   :: x    !x of model grid
   real(kdkind),allocatable,dimension(:)   :: y    !y of model grid
   real(kdkind),allocatable,dimension(:)   :: weight    ! factor to multiply source  runoff
end type runoff_type
   
type :: coast_type
   integer(int32),allocatable,dimension(:) :: i
   integer(int32),allocatable,dimension(:) :: j
   real(kdkind),allocatable,dimension(:)   :: x
   real(kdkind),allocatable,dimension(:)   :: y
   real(kdkind),allocatable,dimension(:)   :: area
end type coast_type
   
type :: source_type
   integer(int32),allocatable,dimension(:) :: i
   integer(int32),allocatable,dimension(:) :: j
   real(kdkind),allocatable,dimension(:)   :: x
   real(kdkind),allocatable,dimension(:)   :: y
   real(kdkind),allocatable,dimension(:)   :: area
   real(kdkind),allocatable,dimension(:)   :: dis
   integer(int32),allocatable,dimension(:) :: idx
end type source_type




real(kind=real64),parameter             :: DEG2RAD = asin(1.0_real64)/90.0_real64  ! PI/180
real(kind=real64),parameter             :: REARTH = 6.371d3  !km



real(kdkind), allocatable,dimension(:,:) :: pos_source,pos_coast,pos_nn
type(kdtree2_result),allocatable,dimension(:) :: results
type(kdtree2),pointer    :: tree

type(runoff_type) :: runoff
type(coast_type) :: coast,coast_nn
type(source_type) :: source

integer(kind=int32) :: nargs
integer(kind=int32) :: method, weights

integer(kind=int32) :: i,j,n
integer(kind=int32) :: i_cyc=-1
integer(kind=int32) :: iptr,iptr_tmp
integer(kind=int32) :: num_max=20
integer(kind=int32) :: num_found
integer(kind=int32) :: num_nn=20
integer(kind=int32) :: num_new
integer(kind=int32) :: n_coast                 ! number of model coastal points
integer(kind=int32) :: n_source               ! Number of source runoff tiles

real(kind=real64)   :: x,y
real(kind=real64)   :: escale,area_source,dist,x_source,y_source,z_source
real(kind=real64)   :: R2_max=(50.0/REARTH)**2 ! 50km
real(kind=real64)   :: R_max=50.0 ! 50km

real(kind=real64),allocatable,dimension(:)     :: wgt_tmp       ! temporary array for holding local weights.
integer(kind=int32),allocatable,dimension(:)   :: n_distribute  ! Number of points source is distibuted to.
integer(kind=int32),allocatable,dimension(:)   :: i_tmp
integer(kind=int32),allocatable,dimension(:)   :: idx_tmp
integer(kind=int32),allocatable,dimension(:)   :: idx_i,idx_j

logical,allocatable,dimension(:)             :: good

character(len=12)   :: carg,carg1

! Get info on the grid from input
nargs=command_argument_count()
do i = 1,nargs,2
   call get_command_argument(i,carg)
   call get_command_argument(i+1,carg1)
   select case(trim(carg))
! Placeholder
   case('-m')
      read(carg1,'(I)') method
   case('-w')
      read(carg1,'(I)') weights
   case('-n')
      read(carg1,'(I)') num_nn
      num_max=num_nn
   case('-r')
      read(carg1,*) r_max
      R2_max=(R_max/REARTH)**2
   case DEFAULT
        write(*,*) 'Unknown argument',trim(carg)
        stop
   end select
enddo

! Write summary of methods to be used

select case(method)
  case(1) 
     write(*,*) 'Selecting ',num_nn, 'nearest neighbours up to',r_max,'km distant'
  case(2) 
     write(*,*) 'Selecting nearest neighbours up to',r_max,'km distant up to a maximum of',num_max
  case DEFAULT
     write(*,*) 'Unknown method. Exiting'
     stop 1
end select
select case(weights)
  case(1) 
     write(*,*) 'Equal weighting'
  case(2) 
     write(*,*) 'Gaussian weighting by distance from original source'
  case(3) 
     write(*,*) 'Gaussian weighting by distance from nearest neighbour to source'
  case DEFAULT
     write(*,*) 'Unknown weighting method. Exiting'
     stop 1
end select

!call read_connection_file_nn(source,coast,coast_nn,qc_nn)
call read_connection_file_nn(source,coast,coast_nn)

n_coast = size(coast%x)
n_source = size(source%x) ! Same size for coast_nn and qc_nn


write(*,*) 'Creating kdtree'

! Points on surface of sphere radius 1
allocate(pos_coast(3,n_coast))
do i = 1,n_coast
   pos_coast(1,i) =  cos(coast%x(i)*DEG2RAD)*cos(coast%y(i)*DEG2RAD)
   pos_coast(2,i) =  sin(coast%x(i)*DEG2RAD)*cos(coast%y(i)*DEG2RAD)
   pos_coast(3,i) =  sin(coast%y(i)*DEG2RAD)
enddo
! Create tree we would like results to be sorted.

tree => kdtree2_create(pos_coast,sort=.true.,rearrange=.true.)  !?? on the rearrange

allocate(pos_source(3,n_source),pos_nn(3,n_source))

do i = 1, n_source
   x = source%x(i)*DEG2RAD
   y = source%y(i)*DEG2RAD
   pos_source(1,i) =  cos(x)*cos(y)
   pos_source(2,i) =  sin(x)*cos(y)
   pos_source(3,i) =  sin(y)
   x = coast_nn%x(i)*DEG2RAD
   y = coast_nn%y(i)*DEG2RAD
   pos_nn(1,i) =  cos(x)*cos(y)
   pos_nn(2,i) =  sin(x)*cos(y)
   pos_nn(3,i) =  sin(y)
enddo

write(*,*) 'Finding targets and weights'

allocate(n_distribute(n_source))
n_distribute = 0

select case (method)
   case(1)  ! Nearest neighbours of target
      allocate(idx_tmp(num_nn*n_source))
      allocate(i_tmp(num_nn*n_source))
      allocate(results(num_nn))
      allocate(idx_i(num_nn),idx_j(num_nn))
      allocate(good(num_nn))

      iptr=0
      do i = 1, n_source
         call kdtree2_n_nearest(tp=tree,qv=pos_nn(:,i),nn=num_nn,results=results)
         do n = 1,num_nn
            idx_i(n) = coast%i(results(n)%idx)
            idx_j(n) = coast%j(results(n)%idx)
         enddo

         call test_adjacency(idx_i,idx_j,i_cyc,num_nn,good)  ! test if we can fid a path fro the closest

         do n = 1, num_nn
            if ( idx_i(n) > 0 .and. results(n)%dis < R2_max .and. good(n)) then
               iptr = iptr + 1
               idx_tmp(iptr) = results(n)%idx
               i_tmp(iptr) = i
               n_distribute(i)   = n_distribute(i) + 1
            endif
         enddo
      enddo

   case(2)  ! Nearest neighbours of nearest neighbour(!) within R^2
      allocate(idx_tmp(num_max*n_source))
      allocate(results(num_max))
      allocate(idx_i(num_max),idx_j(num_max))
      allocate(good(num_max))

      iptr=0

      do i = 1, n_source
         call kdtree2_r_nearest(tp=tree,qv=pos_nn(:,i),r2=R2_max,nfound=num_found, &
                                             nalloc=num_max,results=results)
         if( num_found > num_max ) then
           write(*,*) 'Extending number of results' , num_max, ' to ', num_found
           num_new = num_found
           deallocate(results)
           allocate(results(num_found))
           call kdtree2_r_nearest(tp=tree,qv=pos_nn(:,i),r2=R2_max,nfound=num_found, &
                                             nalloc=num_new,results=results)
         endif

         do n = 1,min(num_max,num_found)
            idx_i(n) = coast%i(results(n)%idx)
            idx_j(n) = coast%j(results(n)%idx)
         enddo
            
         call test_adjacency(idx_i,idx_j,i_cyc,min(num_max,num_found),good)  ! test if we can fid a path fro the closest
         do n = 1, min(num_max,num_found)
            if ( idx_i(n) > 0 .and. good(n) ) then
               iptr = iptr + 1
               idx_tmp(iptr) = results(n)%idx
               n_distribute(i)   = n_distribute(i) + 1
            endif
         enddo
      enddo

   case DEFAULT
      write(*,*) 'Invalid method'
      stop 1
end select
               
! We now  have our source indices source%XXX(1:n_source) and ! coast%XXX(idx_tmp(1:iptr))
! We now need to get weights

deallocate(coast_nn%i,coast_nn%j,coast_nn%x,coast_nn%y,coast_nn%area)
deallocate(results)
deallocate(idx_i,idx_j)

write(*,*) 'Number of coastal points to be mapped to',iptr

allocate(runoff%weight(iptr),runoff%idx(iptr),runoff%i(iptr),runoff%j(iptr),runoff%x(iptr),runoff%y(iptr))
do i=1,iptr
   runoff%i(i)=coast%i(idx_tmp(i))
   runoff%j(i)=coast%j(idx_tmp(i))
   runoff%x(i)=coast%x(idx_tmp(i))
   runoff%y(i)=coast%y(idx_tmp(i))
enddo
do i=1,iptr
   runoff%idx(i)=i_tmp(i)
enddo
   

! Equal all points weights by their area 

allocate(wgt_tmp(maxval(n_distribute)))
select case(weights)
   case(1)
      iptr = 0
      do i = 1, n_source
         area_source = source%area(i)
         iptr_tmp = iptr
         do j = 1, n_distribute(i)
              wgt_tmp(j) = area_source/n_distribute(i)
         enddo
         iptr = iptr_tmp
         do j = 1, n_distribute(i)
              iptr = iptr + 1
              runoff%weight(iptr)=wgt_tmp(j)/coast%area(idx_tmp(iptr))
         enddo
      enddo
   case(2)
      iptr = 0
      escale = (REARTH/50.0)**2
      do i = 1, n_source
         area_source = source%area(i)
         iptr_tmp = iptr
         x_source = pos_source(1,i)
         y_source = pos_source(2,i)
         z_source = pos_source(3,i)
         do j = 1, n_distribute(i)
              iptr = iptr + 1
              dist = (pos_coast(1,idx_tmp(iptr))-x_source)**2 + &
                     (pos_coast(2,idx_tmp(iptr))-y_source)**2 + &
                     (pos_coast(3,idx_tmp(iptr))-z_source)**2 
              wgt_tmp(j) = exp(-dist/escale)
         enddo
         wgt_tmp(1:n_distribute(i)) =  wgt_tmp(1:n_distribute(i))*area_source/sum(wgt_tmp(1:n_distribute(i)))
         iptr = iptr_tmp
         do j = 1, n_distribute(i)
              iptr = iptr + 1
              runoff%weight(iptr)=wgt_tmp(j)/coast%area(idx_tmp(iptr))
              !runoff%weight(idx_tmp(iptr))=wgt_tmp(j)
         enddo
      enddo
   case(3)
      iptr = 0
      escale = (REARTH/50.0)**2
      do i = 1, n_source
         area_source = source%area(i)
         iptr_tmp = iptr
         x_source = pos_nn(1,i)
         y_source = pos_nn(2,i)
         z_source = pos_nn(3,i)
         do j = 1, n_distribute(i)
              iptr = iptr + 1
              dist = (pos_coast(1,idx_tmp(iptr))-x_source)**2 + &
                     (pos_coast(2,idx_tmp(iptr))-y_source)**2 + &
                     (pos_coast(3,idx_tmp(iptr))-z_source)**2 
              wgt_tmp(j) = exp(-dist/escale)
         enddo
         wgt_tmp(1:n_distribute(i)) =  wgt_tmp(1:n_distribute(i))*area_source/sum(wgt_tmp(1:n_distribute(i)))
         iptr = iptr_tmp
         do j = 1, n_distribute(i)
              iptr = iptr + 1
              runoff%weight(iptr)=wgt_tmp(j)/coast%area(idx_tmp(iptr))
         enddo
      enddo
end select   
         
call write_weights_file(source,runoff)
write(*,*) 'Done'

contains

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

subroutine write_weights_file(source,runoff)
   type(runoff_type), intent(in) :: runoff
   type(source_type), intent(in) :: source
   integer(int32) :: ir,is
   integer(int32) :: ncid, did_ir, did_is
   integer(int32) :: runoff_i_id,runoff_j_id
   integer(int32) :: runoff_x_id,runoff_y_id
   integer(int32) :: runoff_weight_id, runoff_idx_id
   integer(int32) :: source_i_id,source_j_id

   ir=size(runoff%x)
   is=size(source%x)

   call handle_error(nf90_create(path='runoff_weights.nc',cmode=ior(nf90_netcdf4,nf90_clobber),ncid=ncid))
   call handle_error(nf90_def_dim(ncid,'ir',ir,did_ir))
   call handle_error(nf90_def_dim(ncid,'is',is,did_is))
   call handle_error(nf90_def_var(ncid,'runoff_i',nf90_int,did_ir,runoff_i_id))
   call handle_error(nf90_put_att(ncid,runoff_i_id,'long_name','model i index'))
   call handle_error(nf90_def_var(ncid,'runoff_j',nf90_int,did_ir,runoff_j_id))
   call handle_error(nf90_put_att(ncid,runoff_j_id,'long_name','model j index'))
   call handle_error(nf90_def_var(ncid,'runoff_x',nf90_double,did_ir,runoff_x_id))
   call handle_error(nf90_put_att(ncid,runoff_x_id,'long_name','model longitude'))
   call handle_error(nf90_put_att(ncid,runoff_x_id,'units','degrees_E'))
   call handle_error(nf90_def_var(ncid,'runoff_y',nf90_double,did_ir,runoff_y_id))
   call handle_error(nf90_put_att(ncid,runoff_y_id,'long_name','model latitude'))
   call handle_error(nf90_put_att(ncid,runoff_y_id,'units','degrees_N'))
   call handle_error(nf90_def_var(ncid,'runoff_weight',nf90_double,did_ir,runoff_weight_id))
   call handle_error(nf90_put_att(ncid,runoff_weight_id,'long_name','multiplication factor'))
   call handle_error(nf90_put_att(ncid,runoff_weight_id,'units','none'))
   call handle_error(nf90_def_var(ncid,'source_index',nf90_int,did_ir,runoff_idx_id))
   call handle_error(nf90_put_att(ncid,runoff_idx_id,'long_name','linear index of source array'))

   call handle_error(nf90_def_var(ncid,'source_i',nf90_int,did_is,source_i_id))
   call handle_error(nf90_put_att(ncid,source_i_id,'long_name','runoff source i index'))
   call handle_error(nf90_def_var(ncid,'source_j',nf90_int,did_is,source_j_id))
   call handle_error(nf90_put_att(ncid,source_j_id,'long_name','runoff source j index'))

   call handle_error(nf90_enddef(ncid,h_minfree=4096))

   call handle_error(nf90_put_var(ncid,runoff_i_id,runoff%i))
   call handle_error(nf90_put_var(ncid,runoff_j_id,runoff%j))
   call handle_error(nf90_put_var(ncid,runoff_x_id,runoff%x))
   call handle_error(nf90_put_var(ncid,runoff_y_id,runoff%y))
   call handle_error(nf90_put_var(ncid,runoff_weight_id,runoff%weight))
   call handle_error(nf90_put_var(ncid,runoff_idx_id,runoff%idx))

   call handle_error(nf90_put_var(ncid,source_i_id,source%i))
   call handle_error(nf90_put_var(ncid,source_j_id,source%j))
   call handle_error(nf90_close(ncid))

end subroutine write_weights_file

subroutine test_adjacency(idx_i,idx_j,i_cyc,n,good)  ! test if we can fid a path fro the closest
integer(kind=int32),dimension(:),intent(in) :: idx_i,idx_j
integer(kind=int32),intent(in)              :: i_cyc  ! If positive i_cyc is length of cycle in i
integer(kind=int32),intent(in)              :: n     ! number of points to test
logical,dimension(:),intent(out)            :: good

good = .true. ! placeholder
end subroutine test_adjacency

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
    
end program create_runoff_weights
