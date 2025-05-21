subroutine plancknu(t,nu,bnu)

  ! calculate planck function,
  ! for given temperature, frequency

  use constants
  implicit none

  real*8 :: nu,bnu,t,xpn,hnkt

  hnkt=h/k*nu/t
  if(hnkt.gt.170.d0) hnkt=170.d0
  xpn=exp(hnkt)
  bnu=2.d0*h*nu*nu/c*nu/c/(xpn-1.d0)

  return
end subroutine plancknu

