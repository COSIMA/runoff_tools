program process_runoff
  ! Take the weights file and the model file and write out our new file
  ! i.e. spread the runoffs according to a cunning plan!
  !
  !
  use iso_fortran_env
  use netcdf
  use runoff_modules
  implicit none

  type(model_type) :: model
  type(source_type) :: source
  type(nc_info_type) :: source_info,model_info

  integer(kind=int32) :: nargs

  integer(kind=int32) :: i,irecs


  real(kind=real32),allocatable,dimension(:,:)   :: runoff_source,runoff_model         ! 

  character(len=128)   :: carg,carg1
  character(len=32)   :: var_name

  nargs=command_argument_count()
  if ( nargs < 4 ) then
     write(*,*) 'ERROR: Usage :: process_runoff -i infile -o outfile'
     stop 1
  endif
  var_name = 'runoff'
  do i = 1,nargs,2
     call get_command_argument(i,carg)
     call get_command_argument(i+1,carg1)
     select case(trim(carg))
     case('-i')
        source_info%fname=trim(carg1)
     case('-o')
        model_info%fname=trim(carg1)
        ! Placeholder
     case('-v')
        var_name=trim(carg1)
     case DEFAULT
        write(*,*) 'Unknown argument',trim(carg)
        stop
     end select
  enddo

  !call read_connection_file_nn(source,coast,coast_nn,qc_nn)
  write(*,*) 'Reading weights'
  call read_weights_file(source,model)
  write(*,*) 'Reading source info'
  call get_source_info(source_info,var_name)
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
    integer(kind=int32)  :: skipped

    runoff_model = 0.0_real32
    skipped = 0
    call handle_error(nf90_get_var(source_info%ncid,source_info%vid,runoff_source,start=(/1,1,record/), &
         count=(/source_info%idim,source_info%jdim,1/)))
    call handle_error(nf90_get_var(source_info%ncid,source_info%tid,time,start=(/record/)))
    do i = 1,model%npts
       i_s = source%i(model%idx(i))
       j_s = source%j(model%idx(i))

       if (runoff_source(i_s,j_s) > 0) then
          runoff_model(model%i(i),model%j(i)) = runoff_model(model%i(i),model%j(i)) + runoff_source(i_s,j_s)*model%weight(i)
       else
          skipped = skipped + 1
       end if
    enddo
    print '(A,I0,A)','Skipped ',skipped,' zero records'
    call handle_error(nf90_put_var(model_info%ncid,model_info%vid,runoff_model,start=(/1,1,record/), &
         count=(/model_info%idim,model_info%jdim,1/)))
    call handle_error(nf90_put_var(model_info%ncid,model_info%tid,time,start=(/record/)))
  end subroutine process_record
end program process_runoff
