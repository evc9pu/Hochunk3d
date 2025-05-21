! ********************************************************

subroutine initpspot(ispotph)

  ! called by subroutine disk
  ! initialize variables for new photon which scatters in disk

  use tts_mod
  use stokesvar_mod
  use random
  use spot_mod
  use constants

  implicit none

  real*8 :: xran
  real*8 :: sthet,cthet,sphi,cphi,dotsp1,dotsp2,dellon,costhp

  integer :: ispotph
  real*8,external :: atmosfreq

  ! print*,'ispotph, nspot, npspot ',ispotph,nspot,npspot

  ! incident radiation chosen to be unpolarized.
  ! can put in a desired polarization in the limb darkening
  ! subroutine.
  sip=1.d0
  sqp=0.d0
  sup=0.d0
  svp=0.d0

  ! sample cosb,sinb, position on star in latitude.
100 continue
  do
     xran=ran()
     cosb=1.d0-2.d0*xran
     sinb=sqrt(1.d0-cosb**2)
     if(cosb.ne.0.d0.and.abs(cosb).lt.1.d0) exit
     print *,'WARNING: cosb is zero or +/-1, resampling'
  end do

  ! sample azimuthal coord (longitude)
  xran=ran()
  lp=r2p*xran
  cphi=cos(lp)
  sphi=sin(lp)
  sthet=sinb
  cthet=cosb

  ! **************SPOT STUFF***************************
  
  if (ispot.eq.1) then
    
    !ispotph=1
    
    dotsp1=scspot1*sthet*cphi+ &
       & ssspot1*sthet*sphi+ &
       & cspot1*cthet
       dotsp2=scspot2*sthet*cphi+ &
       & ssspot2*sthet*sphi+ &
       & cspot2*cthet
       
  ! if((dotsp1.lt.csprad).and.    ! for two spots
  ! $       (dotsp2.lt.csprad))then    ! that are filled in

  ! if( ((dotsp1.lt.csprad).or.(dotsp1.gt.cspradin)) .and.
  ! $         ((dotsp2.lt.csprad).or.(dotsp2.gt.cspradin)) )
  ! $        then

  ! if(dotsp1.lt.csprad)then  ! for one spot filled in
  ! probspot=ampl
  ! else
  ! probspot=1.d0
  ! end if
  ! xran=ran()
  ! if(xran.gt.probspot) then
  ! goto 100
  ! end if

  ! baw, 6/27/00, doing it differently;
  ! use equation 3 from Wood et al 2000.  If iphot < N_spot, put photon
  ! source at spot.  otherwise, put it outside spot.
  ! BAW, 20080828 change to, if ispotph=1...  

    if (nspot.eq.1) then
     ! print*,'dotsp1,csprad ',dotsp1,csprad
     ! this is for a single spot.
     ! if (iphot.lt.npspot) then
      if (ispotph.eq.1) then
        ! note , csprad goes from 1 to 0 as theta goes from 0 to 90
        ! selecting hotspot photon
          if (dotsp1.lt.csprad) go to 100
      !else
        !2010 feb 20, don't reject on stellar photons, let them be in hotspot if they fall there
       ! if (dotsp1.gt.csprad) go to 100
      end if
    else
     ! this is for 2 spots
     ! if (iphot.lt.npspot) then
      if (ispotph.eq.1) then
        if (dotsp1.lt.csprad.and.dotsp2.lt.csprad) go to 100
      !else
      !2010 feb 20, don't reject on stellar photons, let them be in hotspot if they fall there
      !if (dotsp1.gt.csprad.or.dotsp2.gt.csprad) go to 100
      end if
    end if
    
!   print*,'dotsp1,csprad',dotsp1,csprad    
  
  !if photon is in spot region, set ispotph=1 to emit at tshock temperature  
  !if fspot=1, set ispotph=0 so photon is emitted from stellar spectrum.
    if (fspot.lt.1.d0) then
    ! spot 1
      if (lp.lt.pi) then
        dellon=lp
      else
        dellon=r2p-lp
      endif
      costhp=cspot1*cosb+sspot1*sinb*cos(dellon)
      if (costhp.gt.csprad) ispotph=1
    ! spot 2
      if (nspot.eq.2) then
        if (lp.lt.pi) then
          dellon=pi-lp
        else
          dellon=lp-pi
        endif
        costhp=cspot2*cosb+sspot2*sinb*cos(dellon)
        if (costhp.gt.csprad) ispotph=1
      endif
    endif

  ! else
    
  !  ispotph=0
    
  endif

  ! ***************************************************

  call initp_sample()
  call initp_transform()
  
  ! nub=bbfreq(real(tspot))
  
  if (fspot.eq.1.) then
      if (iplanckst.eq.0) then
         nub=atmosfreq() ! this atmosphere file would need to have T=tshock
      else
         nub=bbfreq(tshock)
      end if
   else if(ispotph.eq.0) then
      if (iplanckst.eq.0) then
         nub=atmosfreq()
      else
         nub=bbfreq(tstar)
      end if
   else
      nub=bbfreq(tshock)
   end if

  ! if (ispotph.eq.0) then
  !   if (iplanckst.eq.0) then
  !      nub=atmosfreq()
  !   else
  !      nub=bbfreq(tstar)
  !   end if
  ! else
  !   nub=bbfreq(tshock)
  ! end if
  
  ! print*,'ispotph',ispotph,'tshock',tshock*11605.

  return
end subroutine initpspot


! **********************************************************
