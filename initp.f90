! ********************************************************

subroutine initp(ispotph)

  ! called by subroutine disk
  ! initialize variables for new photon which scatters in disk

  use tts_mod
  use stokesvar_mod
  use spot_mod
  use random
  use constants

  implicit none

  real*8 :: xran,dellon,costhp
  integer :: ispotph
  real*8,external :: atmosfreq

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
  
  !if photon is in spot region, set ispotph=1 to emit at tshock temperature
  if (ispot.eq.1) then
    if (fspot.eq.1.d0) then
      ispotph=0
    else
      ! spot 1
      if (lp.lt.pi) then
          dellon=lp
      else
          dellon=r2p-lp
      endif
      costhp=cspot1*cosb+sspot1*sinb*cos(dellon)
      if (costhp.lt.csprad) ispotph=1
      ! spot 2
      if (nspot.eq.2) then
        if (lp.lt.pi) then
          dellon=pi-lp
        else
          dellon=lp-pi
        endif
        costhp=cspot2*cosb+sspot2*sinb*cos(dellon)
        if (costhp.lt.csprad) ispotph=1
      endif
    endif
  endif

   if (ispotph.eq.0) then
     if (iplanckst.eq.0) then
        nub=atmosfreq()
     else
        nub=bbfreq(tstar)
     end if
    else
     nub=bbfreq(tshock)
    end if
   wave=1.2398d0/nub

  return
end subroutine initp


! **********************************************************
