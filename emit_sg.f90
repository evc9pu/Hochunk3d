subroutine emit_sg

  use tts_mod
  use stokesvar_mod
  use random
  use constants

  implicit none

  ! incident radiation chosen to be unpolarized.
  ! can put in a desired polarization in the limb darkening
  ! subroutine.
  sip=1.d0
  sqp=0.d0
  sup=0.d0
  svp=0.d0

  ! sample theta
  cost=1.d0-2.d0*ran()
  sint=sqrt(1.d0-cost**2)

  ! sample phi
  phi=r2p*ran()
  cosp=cos(phi)
  sinp=sin(phi)

  return

end subroutine emit_sg
