module log_mod

  implicit none
  save

  interface diskdat_write
     module procedure diskdat_write_integer
     module procedure diskdat_write_real4
     module procedure diskdat_write_real8
     module procedure diskdat_write_logical
  end interface

contains

  subroutine diskdat_section(title)
    implicit none
    character(len=*),intent(in) :: title
    write(12,'(A20)') repeat("=",20)
    write(12,'("Section : ")',advance='no')
    write(12,*) title
    write(12,'(A20)') repeat("=",20)
  end subroutine diskdat_section

  subroutine diskdat_write_logical(name,value,description)
    implicit none
    character(len=*),intent(in) :: name,description
    logical,intent(in) :: value
    character(len=16) :: string
    if(value) then
       string = "T"
    else
       string = "F"
    end if
    call diskdat_write_str(name,string,description)
  end subroutine diskdat_write_logical

  subroutine diskdat_write_real4(name,value,description)
    implicit none
    character(len=*),intent(in) :: name,description
    real,intent(in) :: value
    character(len=16) :: string
    write(string,'(ES16.7)') value
    call diskdat_write_str(name,string,description)
  end subroutine diskdat_write_real4

  subroutine diskdat_write_real8(name,value,description)
    implicit none
    character(len=*),intent(in) :: name,description
    real*8,intent(in) :: value
    character(len=16) :: string
    write(string,'(ES16.7)') value
    call diskdat_write_str(name,string,description)
  end subroutine diskdat_write_real8

  subroutine diskdat_write_integer(name,value,description)
    implicit none
    character(len=*),intent(in) :: name,description
    integer,intent(in) :: value
    character(len=16) :: string
    write(string,'(I16)') value
    call diskdat_write_str(name,string,description)
  end subroutine diskdat_write_integer

  subroutine diskdat_write_str(keyword,string,comment)
    implicit none
    character(len=*),intent(in) :: keyword,comment
    character(len=16),intent(in) :: string
    character(len=80) :: card
    write(card,'(A8," = ",A16," / ",A50)') keyword,string,comment
    write(12,'(A80)') card
  end subroutine diskdat_write_str

end module log_mod
