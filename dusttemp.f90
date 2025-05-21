! ********************************************************************
! Function to calculate dust temperature from condition of radiative
! equilibrium
! ********************************************************************

     real*8 function dusttemp(ir,it,ip,na,T0,ltot,nphot,idust2)

  use grid_mod
  use dust_mod
  use constants
  implicit none

  integer :: ir,it,ip             ! grid cell index
  integer :: nphot,idust2
  real*8 :: T0,Td,Told
  real*8 :: ltot,na

  real*8 :: KapP

  real*8,parameter :: eps = 1.d-5

  if (T0 .ne. 0.d0) then
     Td=T0
  else
     Td=0.1d0/11605.d0
  end if

  Told=0.5d0*Td

  do while( abs(Td-Told) .gt. eps*Td)

     Told=Td

  !   Td=Tstar*( (na*pi*Rstar*Rstar) / (nphot*kappav(idust2)* KapP(Td,idust2)) )**0.25d0
      Td=( (na*ltot)/(4.d0*sigt*nphot*kappav(idust2)*KapP(Td,idust2) ) )**0.25d0 / 11605.d0

     if(massarr(ir,it,ip,idust2) > 0.d0) then
     !  Td=Td*263818.1192d0/massarr(ir,it,ip,idust2)
       Td=Td/massarr(ir,it,ip,idust2)
     else
        Td=0.1d0/11605.d0
     end if
     
  !   print*,'in dusttemp, td,ir,it',td * 11605.d0,ir,it

     if (Td.gt.1600.d0/11605.d0) Td=1600.d0/11605.d0

     Td=max(Td,0.1d0/11605.d0)

  end do

  if (Td.ge.Tsub/11605.d0) nodust(ir,it,ip,idust2)=.true.

  dusttemp=Td

  return

end function dusttemp

