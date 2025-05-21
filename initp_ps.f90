! ********************************************************

subroutine initp_ps(ispotph)

  ! called by subroutine disk
  ! initialize variables for new photon which scatters in disk
  ! for really large envelope, makes approximation of point source.
  ! program still set up for emitting from a star, but assume direction
  ! is radial.

  use tts_mod
  use stokesvar_mod
  use random
  use constants

  implicit none

  real*8 :: xran,check
  integer :: ispotph
  character(len=20) :: routine
  real*8,external :: atmosfreq

  routine='initp_ps'

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
     cosb=check(cosb,routine)
     sinb=sqrt(1.d0-cosb**2)
     sinb=check(sinb,routine)
     if(cosb.ne.0.d0) exit
     print *,'WARNING: cosb is zero, resampling'
  end do

  ! assume direction is same as position.  same as assuming it's
  ! a point source.
  cost=cosb
  sint=sinb

  ! sample azimuthal coord (longitude)
  xran=ran()
  lp=r2p*xran

  phi=lp
  sinp=sin(phi)
  cosp=cos(phi)
  sinp=check(sinp,routine)
  cosp=check(cosp,routine)

  ! print*,'ispotph,iplanckst',ispotph,iplanckst
  if (ispotph.eq.0) then
     ! print*,'iplanckst',iplanckst
     if (iplanckst.eq.0) then
        nub=atmosfreq()
     else
        nub=bbfreq(tstar)
        ! print*,'hi,tstar,wave',tstar,1.2398d0/nub
     end if
  else
     nub=bbfreq(tshock)
  end if
  ! nub=atmosfreq()
  ! nub=bbfreq(real(tstar))
  ! wave=1.2398d0/nub

  return
end subroutine initp_ps


! **********************************************************
