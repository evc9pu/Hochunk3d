module messages

  implicit none
  save

  logical :: show_warnings = .true.
  integer :: count_warnings

  logical :: show_errors = .true.
  integer :: count_errors

  logical :: show_debug = .true.

contains

  subroutine set_debug(flag)
    implicit none
    integer,intent(in) :: flag
    show_debug = flag==1
  end subroutine set_debug

  subroutine set_warnings(flag)
    implicit none
    integer,intent(in) :: flag
    show_warnings = flag==1
  end subroutine set_warnings

  subroutine warning(message,from)
    implicit none
    character(len=*),intent(in) :: message,from
    if(show_warnings) then
       print *,''
       print *,'>WARNING : ',message
       print *,'>ORIGIN  : ',from
       print *,''
    end if
    count_warnings = count_warnings + 1
  end subroutine warning

  subroutine set_errors(flag)
    implicit none
    integer,intent(in) :: flag
    show_errors = flag==1
  end subroutine set_errors

  subroutine error(message,from)
    implicit none
    character(len=*),intent(in) :: message,from
    if(show_errors) then
       print *,''
       print *,'>ERROR  : ',message
       print *,'>ORIGIN : ',from
       print *,''
    end if
    count_errors = count_errors + 1
  end subroutine error

end module messages
