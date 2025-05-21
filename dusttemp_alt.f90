! ********************************************************************
! Function to calculate dust temperature from condition of radiative
! equilibrium
! ********************************************************************

real function dusttemp(i,j,k,na,T0,rstar,tstar,nphot,idust2)

  use grid_mod
  use dust_mod
  use constants
  implicit none

  integer :: i,j,k             ! grid cell index
  integer :: nphot,idust2
  real :: T0,Td,Told,eps
  real*8 :: rstar,tstar,na

  real :: KapP

  data eps /1.e-5/

  integer,parameter :: n_iter = 100
  integer :: iter
  real :: required,current,tmax,tmin

  if (massarr(i,j,k) .gt. 0.d0) then

     required = na * pi * (Rstar*rsun)**2. * Tstar**4 / (nphot * massarr(i,j,k)**4.)

     if (required == 0.d0) then

        dusttemp =  0.1 / 11605

     else

        tmin = 0. / 11605.d0
        tmax = 100000. / 11605.d0

        do iter=1,n_iter

           Td = (tmin + tmax) / 2.
           current = Td**4 * kappav(idust2)* KapP(Td,idust2)

           if(abs(current-required)/required.lt.1.e-2) exit

           if(current < required) then
              tmin = Td
           else
              tmax = Td
           end if

        end do

        if(iter==101) then
           print *,'did not converge'
           stop
        end if

        if(Td.ge.Tsub/11605.d0) nodust(i,j,k)=.true.

        Td=max(Td,0.1/11605.d0)

        dusttemp = Td

     end if

  else
     dusttemp = 0.1 / 11605
  end if

end function dusttemp

