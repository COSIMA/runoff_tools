program make_mask

  use iso_fortran_env
  use netcdf
  use runoff_modules

  implicit none

  type(nc_info_type) :: source_info

  character(len=128)   :: carg, carg1, outfile
  character(len=32)   :: var_name

  integer :: i, nargs

  var_name = 'runoff'

  nargs=command_argument_count()
  if ( nargs < 2 ) then
     write(*,*) 'ERROR: Usage :: process_runoff [-v varname] [-o outfile] -f runoff_file'
     stop 1
  endif

  var_name = 'runoff'
  outfile = 'runoff_mask.nc'
  do i = 1,nargs,2
     call get_command_argument(i,carg)
     call get_command_argument(i+1,carg1)
     select case(trim(carg))
     case('-f')
        source_info%fname=trim(carg1)
     case('-v')
        var_name=trim(carg1)
     case('-o')
        outfile=trim(carg1)
     case DEFAULT
        write(*,*) 'Unknown argument',trim(carg)
        stop
     end select
  enddo

  write(*,*) 'Reading source info'
  call get_source_info(source_info, var_name)

  call write_mask(source_info, outfile)

  call handle_error(nf90_close(source_info%ncid))

contains

  subroutine write_mask(source_info, outfile)

    ! Interface variables
    type(nc_info_type),intent(in) :: source_info
    character(len=*) :: outfile

    ! Local variables
    logical :: mask(source_info%idim,source_info%jdim)
    integer :: i

    mask = .false.

    do i = 1,source_info%nrecs
       mask = mask .or. read_record(source_info,i) > 0
    end do

    print *,'Number of runoff pixels: ',count(mask)

    write(*,*) 'Setting up output file'
    call write_mask_file(outfile, mask)

  end subroutine write_mask

end program make_mask
