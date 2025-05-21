real*8 function check(x, routine)

  ! x is cos(theta) or sin(theta)
  ! check to see if it's greater than 1

  implicit none
  real*8 :: x
  character(len=20) :: routine
  character(len=50) :: cmsger

  if (x .gt. 1.d0) then
     if (x .gt. 1.01d0) then
        write(cmsger,'(a,f14.11)')'cos,sin > 1',x
        call ERRMSG('WARNING',routine,cmsger)
        print*,'hi'
     end if
     x = 1.d0
  else if (x .lt. -1.d0) then
     if (x .lt. -1.01d0) then
        print*,'hi'
        write(cmsger,'(a,f14.11)')'cos,sin < -1',x
        call ERRMSG('WARNING',routine,cmsger)
     end if
     x = -1.d0
  end if

  check = x

  return
end function check



