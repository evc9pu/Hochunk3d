! ********************************************************

subroutine initp()

  ! called by subroutine disk
  ! initialize variables for new photon which scatters in disk

  use tts_mod
  use stokesvar_mod
  use random
  use constants

  implicit none

  real*8 :: xran

  ! incident radiation chosen to be unpolarized.
  ! can put in a desired polarization in the limb darkening
  ! subroutine.
  sip=1.d0
  sqp=0.d0
  sup=0.d0
  svp=0.d0

  ! sample cosb,sinb, position on star in latitude.
  do
     xran=ran()
     cosb=1.d0-2.d0*xran
     sinb=sqrt(1.d0-cosb**2)
     if(cosb.ne.0.d0) exit
     print *,'WARNING: cosb is zero, resampling'
  end do

  ! sample azimuthal coord (longitude)
  xran=ran()
  lp=r2p*xran

  call initp_sample()
  call initp_transform()
  
  print*,'tstar',tstar
  nub=bbfreq(tstar)
  ! wave=1.2398d0/nub

  return
end subroutine initp


! **********************************************************
