!> \brief provides limited exception handling
!!
!! To extend with more exception codes, for stopping execution, code should be
!! less than zero, otherwise greater than zero. Zero itself is reserved for
!! the case of no error.
!!
module bexception
 use iso_fortran_env
 implicit none

 integer, parameter :: BEXCEPTION_FATAL = -1   !< use for fatal exception
 integer, parameter :: BEXCEPTION_WARNING = 1  !< use for warning, continuing execution

contains

 !> throwException should be called if the program detects an exception.
 !!
 subroutine throwException(message, errlevel)
  implicit none

  character(*), intent(in) :: message  !< string to print to standard error

  !> level of error:
  !! - < 1 : stop program execution
  !! - > 1 : print message, continue
  integer,      intent(in) :: errlevel

  character(19) :: prepend

  if(errlevel>0) prepend='Warning: demagoque:'
  if(errlevel<0) prepend='Error: demagoque:'

  write(*,*)prepend//message

  if(errlevel<0) then
   write(*,*)
   write(*,*)'Exiting...'
   stop 1
  endif

 end subroutine throwException

end module bexception

