! ********************************************************

subroutine initpout()

  ! called by subroutine disk
  ! initialize variables for new photon which scatters in disk
  !
  ! history:
  ! 00/03/19 (mjw):  set initial stokes vector using Stokes vector
  ! from REAPAR if iveeg.eq.1 (CVEEJ='YES')
  !

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
  sup=0.0d0
  svp=0.d0

  ! sample cosb,sinb, position on star in latitude.
  xran=ran()
  cosb=1.d0-2.d0*xran
  sinb=sqrt(1.d0-cosb**2)

  ! isotropic intensity, sample from n=mu
  xran=ran()
  ! is it isotropic intensity from a surface or isotropic in a volume?
  ! spherical surface
  ! cost=-sqrt(xran)    !emitting towards center
  ! volume
  cost=-sqrt(xran)
  sint=sqrt(1.d0-cost*cost)

  ! sample phis
  xran=ran()
  phi=r2p*xran
  cosp=cos(phi)
  sinp=sin(phi)

  ! sample azimuthal coord (longitude)
  xran=ran()
  lp=r2p*xran

  call initp_transform()

  call rep_isrf(nub)
  ! print*,'sampled freq ',nub
  ! do this for testing of peeling off, easier to differentiate stellar
  ! from envelope spectrum.
  ! nub=atmosfreq()
  ! nub=bbfreq(dble(tstar))
  ! wave=1.2398d0/nub

  return
end subroutine initpout


! **********************************************************
