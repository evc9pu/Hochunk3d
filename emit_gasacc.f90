subroutine emit_gasacc()

  use tts_mod
  use stokesvar_mod, only : sint, cost, phi, sinp, cosp
  use taunum_mod, only : rtot
  use opacin_mod, only : massc
  use random, only : ran, bbfreq
  use constants, only : msol, pi, rsol, r2p, sigt, gn

  implicit none

  real*8 :: tgas, t0

  ! Use Equation (3.23) from Pringle (1981) to find the accretion 
  ! disk temperature:
  !
  !       /    3*G*M*Mdot                         \ 1/4
  ! T_s = |  -------------- * [1 - sqrt(Rstar/R)] |
  !       \  8*pi*R^3*sigma                       /


  ! Find the disk accretion rate in cgs, already calculated
  ! mdot=mdotdisk/3600.d0/24.d0/365.25d0*msol

  ! Use Pringle (1981)'s equation to find the gas temperature
  t0 = 3.d0 * gn * (msol*massc) * mdot / &
       &    (8.d0 * pi * sigt * (rtot*rstar*rsol)**3.d0) 

  tgas = ( t0 * (1.d0 - sqrt(1.d0/rtot)) ) **0.25d0 / 11605.d0

  ! sample photon frequency
  nub = bbfreq(tgas)
  call opacset(nub)

  ! sample theta
  cost=1.d0-2.d0*ran()
  sint=sqrt(1.d0-cost**2.d0)

  ! sample phi
  phi=r2p*ran()
  cosp=cos(phi)
  sinp=sin(phi)

end subroutine emit_gasacc
