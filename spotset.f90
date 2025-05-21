subroutine spotset()

  use stokesvar_mod
  use spot_mod
  use constants

  implicit none

  thspot=thspot*deg2rad
  ! 20080829 BAW, not doing ring
  thspotin=0.d0
  ! thspotin=thspotin*deg2rad
  th0sp1=(90.d0-spotlat)*deg2rad
  phi0sp1=spotlon*deg2rad    !note that spotlon has been set to 0 in reapar which is fine, leave that wayispo.  
  th0sp2=pi-th0sp1
  phi0sp2=phi0sp1+pi

  scspot1=sin(th0sp1)*cos(phi0sp1)
  ssspot1=sin(th0sp1)*sin(phi0sp1)
  cspot1=cos(th0sp1)
  sspot1=sin(th0sp1)
  scspot2=sin(th0sp2)*cos(phi0sp2)
  ssspot2=sin(th0sp2)*sin(phi0sp2)
  cspot2=cos(th0sp2)
  sspot2=sin(th0sp2)

  csprad=cos(thspot)
  cspradin=cos(thspotin)

  spotflag=0

  ! TESTING 20080830
  ! tspot=10000.
  ! thspot=10.
  ! npspot=

  ! determine ratio of spot photons to non-spot based on size and
  ! temperatures
  ! 20080829 don't do this in radeq code, already calculated this
  ! stuff in gridset.f
  ! print*,'wave',wave
  ! nu=2.9979d14/wave    !in microns
  ! print*,'nu',nu
  ! call plancknu(tstar,nu,bnustar)
  ! call plancknu(tspot,nu,bnuspot)

  ! one spot
  ! f=(1.d0-csprad)/(1.d0+csprad)*bnuspot/bnustar
  ! two spots
  ! if (nspot.eq.2) f=2.d0*f
  ! haven't done formula for ring....

  ! npspot=real(np)*f/(1.d0+f)

  ! print*,'ratio of Bspot/Bstar ',bnuspot/bnustar
  ! print*,'ratio of spot photons to star photons ',f
  ! print*,'spot photons, star photons, total ',npspot, np-npspot, np
  ! print*, ' '

  ! stop

  return
end subroutine spotset
