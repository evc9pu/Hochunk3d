module configuration

  implicit none
  save

  private
  public :: load_parameter_file
  public :: get_parameter_value

  integer :: n_par

  character(len=100),allocatable :: names(:)
  character(len=1000),allocatable :: values(:)

  interface get_parameter_value
     module procedure get_parameter_value_logical
     module procedure get_parameter_value_int
     module procedure get_parameter_value_char
     module procedure get_parameter_value_r4
     module procedure get_parameter_value_r8
     module procedure get_parameter_value_r4arr
     module procedure get_parameter_value_intarr
     module procedure get_parameter_value_r8arr
  end interface

contains

  subroutine load_parameter_file(filename)
    implicit none
    character(len=*),intent(in) :: filename
    integer :: ioerr,pos1,pos2,iter,i
    character(len=1000) :: line
    do iter=1,2
       open(unit=10,file=filename,status='old')
       i = 0
       do
          read(10,'(A1000)',iostat=ioerr) line
          if(ioerr.ne.0) exit
          pos1 = index(line,'=')
          pos2 = index(line(pos1+1:),'=') + pos1
          if(pos1.gt.0.and.pos2.gt.0) then
             i = i + 1
             if(iter==2) then
                values(i) = trim(line(1:pos1-1))
                names(i) = trim(adjustl(line(pos1+1:pos2-1)))
             end if
          end if
       end do
       close(unit=10)
       if(iter==1) then
          n_par = i
          allocate(names(n_par),values(n_par))
       end if
    end do
    write(*,'(I0," parameter(s) read from ",3X,A)') n_par,filename
  end subroutine load_parameter_file

  integer function key(name)
    implicit none
    character(len=*),intent(in) :: name
    integer :: i
    do i=1,n_par
       if(trim(name)==trim(names(i))) then
          key = i
          return
       end if
    end do
    key = 0
  end function key

  subroutine get_parameter_value_int(name,value,comment)
    implicit none
    character(len=*),intent(in) :: name,comment
    character(len=3) :: cvalue
    integer,intent(out) :: value
    integer :: ioerr
    if(key(name) > 0) then
       read(values(key(name)),*,iostat=ioerr) value
       if(ioerr.ne.0) then
          read(values(key(name)),*,iostat=ioerr) cvalue
          if(cvalue=='YES'.or.cvalue=='NO') then
             if(cvalue=='YES') value = 1
             if(cvalue=='NO') value = 0
          else
             print *,"Parameter "//trim(name)//" is not an integer"
             stop
          end if
       end if
    else
       write(*,*) "ERROR: parameter "//trim(name)//" not found"
       stop
    end if
    write(*,'(A10,3X,I10,3X,A)') name,value,comment
  end subroutine get_parameter_value_int

!subroutine get_parameter_value_intarr(name,value,comment)
!  implicit none
!  character(len=*),intent(in) :: name,comment
!  integer,intent(out) :: value(:)
!  real*8 :: value_r8(size(value))
!  call get_parameter_value_r8arr(name,value_r8,comment)
!  value = nint(value_r8)
! end subroutine get_parameter_value_intarr

subroutine get_parameter_value_intarr(name,value,comment)
  implicit none
  character(len=*),intent(in) :: name,comment
  integer,intent(out) :: value(:)
  integer :: ioerr,i
  if(key(name) > 0) then
     read(values(key(name)),*,iostat=ioerr) value(:)
     if(ioerr.ne.0) then
        print *,"Parameter "//trim(name)//" is not an array of real numbers"
        stop
     end if
  else
     write(*,*) "ERROR: parameter "//trim(name)//" not found"
     stop
  end if
  write(*,'(A10,1X,I2,3X,I10,3X,A)') name,1,value(1),comment
  do i=2,size(value)
     write(*,'(A10,1X,I2,3X,I10,3X,A)') "",i,value(i),""
  end do
end subroutine get_parameter_value_intarr


  subroutine get_parameter_value_r4(name,value,comment)
    implicit none
    character(len=*),intent(in) :: name,comment
    real,intent(out) :: value
    real*8 :: value_r8
    call get_parameter_value_r8(name,value_r8,comment)
    value = real(value_r8)
  end subroutine get_parameter_value_r4

  subroutine get_parameter_value_r8(name,value,comment)
    implicit none
    character(len=*),intent(in) :: name,comment
    real*8,intent(out) :: value
    integer :: ioerr
    if(key(name) > 0) then
       read(values(key(name)),*,iostat=ioerr) value
       if(ioerr.ne.0) then
          print *,"Parameter "//trim(name)//" is not a real number"
          stop
       end if
    else
       write(*,*) "ERROR: parameter "//trim(name)//" not found"
       stop
    end if
    write(*,'(A10,3X,ES10.4,3X,A)') name,value,comment
  end subroutine get_parameter_value_r8

  subroutine get_parameter_value_r4arr(name,value,comment)
    implicit none
    character(len=*),intent(in) :: name,comment
    real,intent(out) :: value(:)
    real*8 :: value_r8(size(value))
    call get_parameter_value_r8arr(name,value_r8,comment)
    value = real(value_r8)
  end subroutine get_parameter_value_r4arr

  subroutine get_parameter_value_r8arr(name,value,comment)
    implicit none
    character(len=*),intent(in) :: name,comment
    real*8,intent(out) :: value(:)
    integer :: ioerr,i
    if(key(name) > 0) then
       read(values(key(name)),*,iostat=ioerr) value(:)
       if(ioerr.ne.0) then
          print *,"Parameter "//trim(name)//" is not an array of real numbers"
          stop
       end if
    else
       write(*,*) "ERROR: parameter "//trim(name)//" not found"
       stop
    end if
    write(*,'(A10,1X,I2,3X,ES10.4,3X,A)') name,1,value(1),comment
    do i=2,size(value)
       write(*,'(A10,1X,I2,3X,ES10.4,3X,A)') "",i,value(i),""
    end do
  end subroutine get_parameter_value_r8arr

  subroutine get_parameter_value_char(name,value,comment)
    implicit none
    character(len=*),intent(in) :: name,comment
    character(len=*),intent(out) :: value
    integer :: ioerr
    if(key(name) > 0) then
       read(values(key(name)),*,iostat=ioerr) value
       if(ioerr.ne.0) then
          print *,"Parameter "//trim(name)//" is not a string"
          stop
       end if
    else
       write(*,*) "ERROR: parameter "//trim(name)//" not found"
       stop
    end if
    write(*,'(A10,3X,A10,3X,A)') name,value,comment
  end subroutine get_parameter_value_char

  subroutine get_parameter_value_logical(name,value,comment)
    implicit none
    character(len=*),intent(in) :: name,comment
    character(len=3) :: cvalue
    logical,intent(out) :: value
    integer :: ioerr
    if(key(name) > 0) then
       read(values(key(name)),*,iostat=ioerr) cvalue
       if(ioerr.ne.0.or.(trim(cvalue).ne.'YES'.and.trim(cvalue).ne.'NO')) then
          print *,"Parameter "//trim(name)//" should be 'YES' or 'NO'"
          stop
       end if
       value = cvalue == 'YES'
    else
       write(*,*) "ERROR: parameter "//trim(name)//" not found"
       stop
    end if
    write(*,'(A10,3X,L10,3X,A)') name,value,comment
  end subroutine get_parameter_value_logical

end module configuration


