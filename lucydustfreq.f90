! ********************************************************************
! Temperature Correction procedure that calculates photon frequency
! for dust emission from the distribution
!
! kappa_nu*B_nu(T1)      kappa_nu*B_nu(T0)
! N1*----------------- - N0*-----------------
! dP          kappa_P(T1)*B(T1)      kappa_P(T0)*B(T0)
! --- =   -------------------------------------------
! dnu                        N1 - N0
!
! where N0 is the number of photons previously emitted at temperature
! T0 and N1 is the total number to be emitted at temperature T1

! 20070804 BAW, try to modify subroutine to just sample from

! dp/dnu = kappa_nu*B_nu(T1)

! ********************************************************************
real*8 function lucydustfreq(T0,idust2)

  use dust_mod
  use random

  implicit none

  real*8 :: T0,lT0

  integer :: j,i0,i1,i,idust2
  real*8 :: xi,f,P0,P1,P

!  if (T0 .ge. Tsub/11605.d0) then
!     lucydustfreq=bbfreq(T0)
!     return
!  end if
  ! if (T0 .le. 11./11605.d0) then
  ! dustfreq=bbfreq(T0)
  ! return
  ! end if

  lT0=log(T0)

  xi=ran()

  j=int(dble(nT-1)*(lT0-lTint(1))/(lTint(nT)-lTint(1)))+1
  if ( (j .lt. 1) .or. (j .gt. (nT-1)) ) then
     write(*,*) 'index j out of bounds'
     write(*,*) 'T0 =',11605.d0*T0
  end if
  f=(lT0-lTint(j))/(lTint(j+1)-lTint(j))
  if (f.lt.0.d0) then
     f=0.d0
  end if
  if (f.gt.1.d0) then
     f=1.d0
  end if

  i0=1
  i1=nnu

  P0=  (1.d0-f)*kapBint(j,i0,idust2) + f*kapBint(j+1,i0,idust2)
  if(xi.lt.P0) then
     if (P0 .gt. 0.d0) then
        lucydustfreq=T0*xi*nuint(1)/P0
     else
        lucydustfreq=T0*nuint(1)
     end if
     return
  end if

  P1=  (1.d0-f)*kapBint(j,i1,idust2) + f*kapBint(j+1,i1,idust2)
  if (P1 .gt. 1.d0) then
     P1=1.d0
  end if
  if (xi.ge.P1) then
     if (P1 .lt. 1.d0) then
        if (xi.lt.1.d0) then
           lucydustfreq=T0*nuint(nnu)*(1.d0-P1)/(1.d0-xi)
        else
           lucydustfreq=T0*nuint(nnu)
        end if
     else
        lucydustfreq=T0*nuint(nnu)
     end if
     return
  end if

  do while((i1-i0).gt.1)

     i=(i0+i1)/2
     P=  (1.d0-f)*kapBint(j,i,idust2) + f*kapBint(j+1,i,idust2)

     if (P .gt. 1.d0) then
        if (i.ge.nnu-2) then
           write(*,*)'P>1; i,P = ',i,P
        end if
        P=1.d0
     end if

     if (P .lt. P0) then
        write(*,*) &
             & 'P<P0; P0,P,P1 = ',P0,P,P1
        P0=P
        i0=i
     end if
     if (P .gt. P1) then
        write(*,*) &
             & 'P>P0; P0,P,P1 = ',P0,P,P1
        P1=P
        i1=1
     end if

     if (P.gt.xi) then
        i1=i
        P1=P
     else
        i0=i
        P0=P
     end if

  end do

  if (P1.lt.P0) then
     write(*,*) 'P1<P0; P0, P1 = ',P0,P1
     write(*,*) &
          & 'xi,f,P0,P1,P',xi,f,P0,P1,P
  end if

  ! f=(xi-P0)/(P1-P0)
  if (xi.lt.P0) then
     write(*,*)'xi<P0, xi,P0 = ',xi,P0
     f=0.d0
  else if (xi.gt.P1) then
     write(*,*)'xi>P1, xi,P1 = ',xi,P1
     f=1.d0
  else if (P0.le.0.d0) then
     write(*,*) &
          & 'P0<=0, P0,i0,i1 = ',P0,i0,i1
     f=0.5d0
  else if (P1.le.0.d0) then
     write(*,*)'P1<=0, P1 = ',P1
  else if (xi.le.0.d0) then
     write(*,*)'xi<=0, xi = ',xi
  else
     f=log(xi/P0)/log(P1/P0)
  end if

  lucydustfreq= T0 * nuint(i0) * (nuint(i1)/nuint(i0))**f

  if (f.lt.0.d0) then
     f=0.d0
  end if
  if (f.gt.1.d0) then
     f=1.d0
  end if

  return
end function lucydustfreq
