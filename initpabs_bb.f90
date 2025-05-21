subroutine initpabs

  ! reemit thermal photon from surface of star after absorption
  ! by star.  differs from initp because we know where it is
  ! on surface.  otherwise, it's the same

  use tts_mod
  use stokesvar_mod
  use taunum_mod
  use random
  implicit none

  ! incident radiation chosen to be unpolarized.
  ! can put in a desired polarization in the limb darkening
  ! subroutine.
  sip=1.d0
  sqp=0.d0
  sup=0.d0
  svp=0.d0

  cosb=zp/rtot
  sinb=sqrt(1.d0-cosb**2)

  lp=atan2(yp,xp)

  call initp_sample()
  call initp_transform()
  
  ux=sint*cosp
  uy=sint*sinp
  uz=cost

  ! nub=atmosfreq()
  nub=bbfreq(tstar)
  ! wave=1.2398d0/nub

  return
end subroutine initpabs


! **********************************************************
