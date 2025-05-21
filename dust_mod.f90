module dust_mod

  use grid_mod, only: ndg
  implicit none
  save

  real*8,dimension(ndg) :: kappa,kappav,albedo,g,pl,g2,onemg2,g2p1,twog,p1maxinv
  real*8 :: pc,sc,Tcond,Tsub

  integer,parameter :: nlammax = 1001
  integer :: nlambda(ndg)
  real*8,dimension(nlammax,ndg) :: lamdust,nudust,kappad,adust,gdust,pldust

  integer,parameter :: nT=1000
  integer,parameter :: nnu=1025
  real*8 :: lTint(nT),Tint(nT),nuint(nnu),kappap(nT,ndg),kappar(nT,ndg)
  real*8 :: kapBint(nT,nnu,ndg),kapdBint(nT,nnu,ndg),dBint(nT,nnu,ndg),kappaf(ndg)

end module dust_mod
