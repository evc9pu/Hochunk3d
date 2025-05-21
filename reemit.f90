subroutine reemit(T,idust2)

  use tts_mod
  use stokesvar_mod
  use random
  use constants

  implicit none

  real*8 :: T
  real*8, external :: dustfreq,lucydustfreq
  integer :: idust2

  ! incident radiation chosen to be unpolarized.
  ! can put in a desired polarization in the limb darkening
  ! subroutine.
  sip=1.d0
  sqp=0.d0
  sup=0.d0
  svp=0.d0

  cost=1.d0-2.d0*ran()
  sint=sqrt(1.d0-cost**2.d0)

  ! sample phi
  phi=r2p*ran()
  cosp=cos(phi)
  sinp=sin(phi)

  if (ilucy.eq.0) then
     nub=dustfreq(T,idust2)
  else
     nub=lucydustfreq(T,idust2)
  end if

  return

end subroutine reemit
