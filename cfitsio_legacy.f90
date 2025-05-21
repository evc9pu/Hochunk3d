subroutine ftexist(filename, status, exists)
	implicit none
	character(len=*),intent(in) :: filename
	integer :: status,exists
	call ftexest(filename, status, exists)
end subroutine ftexist