! ********************************************************************
! Function to Calculate the Normalized Planck Function (B_nu/B)*(kT/h)
! ********************************************************************
real*8 function Bnu(nu,T)

  implicit none

  real*8,parameter :: norm = 0.1539897d0    ! 15/pi**4

  real*8 :: nu,T,x
  real*8 :: expm1

  x=nu/T

  if (x.eq.0.d0) then
     Bnu=0.d0
  else
     Bnu=norm*x*x*x/expm1(x)  !Bx = (B_nu/B)*(kT/h)
  end if

  return
end function Bnu

! ********************************************************************
! Function to Calculate the temperature derivative of the Normalized
! Planck Function (dB_nu/dT)/(dB/dT) * (kT/h)
! ********************************************************************

real*8 function dBnu(nu,T)

  implicit none

  real*8,parameter :: norm = 0.1539897d0    ! 15/pi**4

  real*8 :: nu,T,x
  real*8 :: expm1

  x=nu/T

  if (x.eq.0.d0) then
     dBnu=0.d0
  else
     dBnu=-norm*x*x*x*x/(4.d0*expm1(x)*expm1(-x))  !x*Bx/(1-exp(-x))/4
  end if

  return
end function dBnu

real*8 function ierfc(x)

  ! ***************************************************************************
  ! Inverse error function complent.  It uses the ration func approx in
  ! Abramowitz and Stegun (AMS-55), eq. 26.2.23, page 933. (Note sqrt(2) diff
  ! between arguments of Q and erfc)
  !
  ! ***************************************************************************

  implicit none

  real*8 :: x,t

  real*8,parameter :: sqhalf = 0.70710678d0

  real*8,parameter :: c0 = 2.515517d0
  real*8,parameter :: c1 = 0.802853d0
  real*8,parameter :: c2 = 0.010328d0

  real*8,parameter :: d1 = 1.432788d0
  real*8,parameter :: d2 = 0.189269d0
  real*8,parameter :: d3 = 0.001308d0

  t=sqrt(2.d0*log(2.d0/x))

  ierfc=(t-(c0+t*(c1+t*c2))/(1.d0+t*(d1+t*(d2+t*d3))))*sqhalf

  return
end function ierfc

real*8 function expm1(x)

  ! ********************************************************************
  ! this function computes the value of exp(x)-1 using the
  ! expansion in Abramowitz and Stegun (sec. 4.2.45).  this
  ! is used to get around the accuracy problems inherent
  ! when subtracting exp(x) - 1 when x is very
  ! small.
  ! ********************************************************************

  real*8 :: x
  real*8 :: arg,sum
  integer :: i

  real*8,dimension(7),parameter :: coef = (/0.9999999995d0,&
       &0.4999999206d0,0.1666653019d0,0.0416573475d0,&
       &0.0083013598d0,0.0013298820d0,0.0001413161d0/)

  real*8,parameter :: maxarg = 0.6931472d0

  arg=x

  if (abs(x) .gt. maxarg) then
     expm1 = exp(arg)-1.0d0
  else
     sum = 0.d0
     do i=7,1,-1
        sum = arg*(coef(i)+sum)
     end do
     expm1 = sum
  end if

  return

end function expm1



