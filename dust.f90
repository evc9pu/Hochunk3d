module kappa_B_int

  implicit none

  real*8 :: nu0,nu1,dnu,kappa0,kappa1,chi0,chi1,Teff

end module kappa_B_int

subroutine opacset(nu)

  use dust_mod
  implicit none

  integer :: i,i0,i1,id
  real*8 :: nu,f,fm1

  do id=1,ndg

     if (nu.ge.nudust(1,id)) then

        g(id)=gdust(1,id)
        pl(id)=pldust(1,id)
        kappa(id)=kappad(1,id)
        albedo(id)=adust(1,id)

     else if (nu.le.nudust(nlambda(id),id)) then

        g(id)=gdust(nlambda(id),id)
        pl(id)=pldust(nlambda(id),id)
        kappa(id)=kappad(nlambda(id),id)
        albedo(id)=adust(nlambda(id),id)

     else

        i0=1
        i1=nlambda(id)

        do while(i1-i0 .gt. 1)
           i=(i1+i0)/2
           if (nu.le.nudust(i,id)) then
              i0=i
           else
              i1=i
           end if
        end do

        f=(nu-nudust(i1,id))/(nudust(i0,id)-nudust(i1,id))
        fm1=1.d0-f
        g(id)=fm1*gdust(i1,id)+f*gdust(i0,id)
        pl(id)=fm1*pldust(i1,id)+f*pldust(i0,id)
        kappa(id)=fm1*kappad(i1,id)+f*kappad(i0,id)
        albedo(id)=fm1*adust(i1,id)+f*adust(i0,id)

        ! f=log(nu/nudust(i1))/log(nudust(i0)/nudust(i1))
        ! g=gdust(i1)*(gdust(i0)/gdust(i1))**f
        ! pl=pldust(i1)*(pldust(i0)/pldust(i1))**f
        ! kappa=kappad(i1)*(kappad(i0)/kappad(i1))**f
        ! albedo=adust(i1)*(adust(i0)/adust(i1))**f

     end if

     ! pc=0.39
     ! sc=1.0
     ! if (pl+pc .gt. 1.d0) then
     ! write(*,*) 'Error - pl+pc > 1'
     ! pc=0.99-pl
     ! end if
     pc=0.d0
     sc=0.d0

     g2(id)=g(id)*g(id)
     g2p1(id)=g2(id)+1.d0
     onemg2(id)=1.d0-g2(id)
     twog(id)=2.d0*g(id)
     p1maxinv(id)=0.9d0*min((1.d0-g(id))**2.d0/(1.d0+g(id)),0.5d0)

  end do

  return
end subroutine opacset


real*8 function chiR(T,idust2)

  use dust_mod
  implicit none

  integer :: k,idust2
  real*8 :: T,lT,f

  lT=log(T)
  k=int(dble(nT-1)*(lT-lTint(1))/(lTint(nT)-lTint(1)))+1

  if ( (k .lt. 1) .or. (k .gt. (nT-1)) ) then
     write(*,*) 'chiR: index k out of bounds'
     write(*,*) '          k = ',k,'  T = ',11605.d0*T
     if (k.lt.1) then
        chiR=kappar(1,idust2)
     else
        chiR=kappar(nT,idust2)
     end if
     return
  end if

  f=(lT-lTint(k))/(lTint(k+1)-lTint(k))

  chiR = (1.d0-f)*kappar(k,idust2)+f*kappar(k+1,idust2)

  return

end function chiR


real*8 function KapP(T,idust2)

  use dust_mod
  implicit none

  integer :: k,idust2
  real*8 :: T,lT,f

  lT=log(T)
  k=int(dble(nT-1)*(lT-lTint(1))/(lTint(nT)-lTint(1)))+1

  if ( (k .lt. 1) .or. (k .gt. (nT-1)) ) then
     write(*,*) 'KapP: index k out of bounds'
     write(*,*) '          k = ',k,'  T = ',11605.d0*T
     if (k.lt.1) then
        KapP=kappap(1,idust2)
     else
        KapP=kappap(nT,idust2)
     end if
     return
  end if

  f=(lT-lTint(k))/(lTint(k+1)-lTint(k))

  KapP = (1.d0-f)*kappap(k,idust2)+f*kappap(k+1,idust2)

  return

end function KapP


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
! ********************************************************************
real*8 function dustfreq(T0,idust2)

  use dust_mod
  use random
  implicit none

  real*8 :: T0,lT0

  integer :: j,i0,i1,i,idust2
  real*8 :: xi,f,P0,P1,P

  if (T0 .ge. Tsub/11605.d0) then
     dustfreq=bbfreq(T0)
     return
  end if

  lT0=log(T0)

  xi=ran()

  j=int(dble(nT-1)*(lT0-lTint(1))/(lTint(nT)-lTint(1)))+1
  if ( (j .lt. 1) .or. (j .gt. (nT-1)) ) then
     write(*,*) 'dustfreq: index j out of bounds'
     write(*,*) 'T0 =',11605.d0*T0
     stop
  end if
  f=(lT0-lTint(j))/(lTint(j+1)-lTint(j))

  if (f.lt.0.d0) then
     if (f.lt.-0.001d0) then
        write(*,*) 'dustfreq: Pf = ',f,' out of bounds [0,1]'
        write(*,*) '          xi = ',xi,'freq = ',1.24d0/dustfreq
     end if
     f=0.d0
  end if

  if (f.gt.1.d0) then
     if (f.gt.1.001d0) then
        write(*,*) 'dustfreq: Pf = ',f,' out of bounds [0,1]'
        write(*,*) '          xi = ',xi,'freq = ',1.24d0/dustfreq
     end if
     f=1.d0
  end if

  i0=1
  i1=nnu

  P0=(1.d0-f)*kapdBint(j,i0,idust2) + f*kapdBint(j+1,i0,idust2)
  if(xi.lt.P0) then
     if (P0 .gt. 0.d0) then
        dustfreq=T0*xi*nuint(1)/P0
     else
        dustfreq=T0*nuint(1)
     end if
     return
  end if

  P1=(1.d0-f)*kapdBint(j,i1,idust2) + f*kapdBint(j+1,i1,idust2)
  if (P1 .gt. 1.d0) then
     P1=1.d0
  end if

  if (xi.ge.P1) then
     if (P1 .lt. 1.d0) then
        if (xi.lt.1.d0) then
           dustfreq=T0*nuint(nnu)*(1.d0-P1)/(1.d0-xi)
        else
           dustfreq=T0*nuint(nnu)
        end if
     else
        dustfreq=T0*nuint(nnu)
     end if
     return
  end if

  do while((i1-i0).gt.1)

     i=(i0+i1)/2
     P=  (1.d0-f)*kapdBint(j,i,idust2) + f*kapdBint(j+1,i,idust2)

     if (P .gt. 1.d0) then
        if (i.ge.nnu-2) then
           write(*,*)'dustfreq: P>1; i,P = ',i,P
        end if
        P=1.d0
     end if

     if (P .lt. P0) then
        write(*,*) 'dustfreq: P<P0; P0,P,P1 = ',P0,P,P1
        P0=P
        i0=i
     end if
     if (P .gt. P1) then
        write(*,*) 'dustfreq: P>P1; P0,P,P1 = ',P0,P,P1
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
     write(*,*) 'dustfreq: P1<P0; P0, P1 = ',P0,P1
     print*,'xi,f,P0,P1,P',xi,f,P0,P1,P
     stop
  end if

  if (xi.lt.P0) then
     write(*,*) 'dustfreq: xi<P0, xi,P0 = ',xi,P0
     f=0.d0
  else if (xi.gt.P1) then
     write(*,*) 'dustfreq: xi>P1, xi,P1 = ',xi,P1
     f=1.d0
  else if (P0.le.0.d0) then
     write(*,*) 'dustfreq: P0<=0, P0,i0,i1 = ',P0,i0,i1
     f=0.5d0
  else if (P1.le.0.d0) then
     write(*,*) 'dustfreq: P1<=0, P1 = ',P1
     stop
  else if (xi.le.0.d0) then
     write(*,*) 'dustfreq: xi<=0, xi = ',xi
     stop
  else
     f=log(xi/P0)/log(P1/P0)
  end if

  dustfreq= T0 * nuint(i0) * (nuint(i1)/nuint(i0))**f

  if (f.lt.0.d0) then
     if (f.lt.-0.001d0) then
        write(*,*) 'dustfreq: Pf = ',f,' out of bounds [0,1]'
        write(*,*) '          xi = ',xi,'freq = ',1.24d0/dustfreq
     end if
     f=0.d0
  end if

  if (f.gt.1.d0) then
     if (f.gt.1.001d0) then
        write(*,*) 'dustfreq: Pf = ',f,' out of bounds [0,1]'
        write(*,*) '          xi = ',xi,'freq = ',1.24d0/dustfreq
     end if
     f=1.d0
  end if

  return

end function dustfreq


subroutine dustinit(dustfile)

  use dust_mod
  implicit none

  integer :: i,iv,id
  real*8 :: lam,f,hgg,kap,pmax,thet,cext,csca
  real*8, parameter :: Tmin = 0.1d0, Tmax = 2100.d0 !Tmax>Tcond*(5/2)**0.25
  real*8 :: numin,numax
  character(len=80) :: dustfile(ndg)

  Tcond=1600.d0/11605.d0

  numax=32.d0                 ! in units of T; i.e., numax=32.k*T
  numin=numax/1024.d0
  

  do id=1,ndg

     open(7,file=dustfile(id),status='old')

     do i=1,nlammax

        read(7,*,end=100) lam,cext,csca,kap,hgg,pmax,thet
        if (lam .le. 0.55d0) iv=i
        lamdust(i,id)=10000.d0*lam
        nudust(i,id)=12398.d0/lamdust(i,id)

        kappad(i,id)=kap
        adust(i,id)=csca/cext
        gdust(i,id)=hgg
        pldust(i,id)=pmax
        
        ! test
!       if (id.ge.5.and.lam.gt.0.4.and.lam.lt.0.1) kappad(i,id)=0.00001*kappad(i,id)

!   for VSG/PAH set first element of wavelength array to small opacity
!   smallest wavelength is 1000 A.  don't know opacity below that.  in opacset, we
!   use this value for all lower wavelengths.  I'd rather use a small number than one I don't know.
        if (i.eq.1) kappad(i,id) = 0.00001*kappad(i,id)

     end do

     write(*,*) 'error, too many wavelengths in dust opacity file'
     stop

100  continue

     close(7)
     nlambda(id)=i-1

     ! normalize opacity array

     f=(5500.d0-lamdust(iv,id))/(lamdust(iv+1,id)-lamdust(iv,id))
     kappav(id)=(1.d0-f)*kappad(iv,id)+f*kappad(iv+1,id)

     do i=1,nlambda(id)
        kappad(i,id)=kappad(i,id)/kappav(id)
     end do

  end do

  ! set up temperature and frequency arrays for Planck mean opacity
  do i=1,nT
     lTint(i)=log(Tmin/11605.d0)+dble(i-1)*log(Tmax/Tmin)/dble(nT-1)
     Tint(i)=exp(lTint(i))
     ! write(*,*) i,Tint(i),lTint(i)
  end do

  do i=1,nnu
     nuint(i)=numin+dble(i-1)*(numax-numin)/dble(nnu-1)
  end do

  call setemiss()

  return

end subroutine dustinit


subroutine setemiss()

  use dust_mod
  use kappa_B_int
  implicit none

  integer :: i,j,id
  real*8 :: dx,chidBint

  real*8, external :: kappaB,kappadB,dB,chidB

  do id=1,ndg

     do i=1,nT

        Teff=Tint(i)

        nu1=0.d0                 ! be sure to multiply by Teff if nu1 .ne. 0
        call opacset(nu1)
        chi1=kappa(id)
        kappa1=(1.d0-albedo(id))*kappa(id) !kappa1 is absorptive opacity

        do j=1,nnu

           nu0=nu1
           nu1=nuint(j)*Teff   !frequency array is scalled to Teff
           dnu=nu1-nu0
           dx=dnu/Teff         !x=h*nu/kT
           chi0=chi1
           kappa0=kappa1
           call opacset(nu1)
           chi1=kappa(id)
           kappa1=(1.d0-albedo(id))*kappa(id)

           if ( (nu0/Teff) .lt. 50.d0) then
              call QSIMP(kappaB,0.d0,1.d0,kapBint(i,j,id))
              call QSIMP(kappadB,0.d0,1.d0,kapdBint(i,j,id))
              call QSIMP(dB,0.d0,1.d0,dBint(i,j,id))
              call QSIMP(chidB,0.d0,1.d0,chidBint)
           else
              kapBint(i,j,id)=0.d0
              kapdBint(i,j,id)=0.d0
              dBint(i,j,id)=0.d0
              chidBint=0.d0
           end if

           ! sum integration step into cumulative integral

           if (j.eq.1) then
              kapBint(i,j,id)=dx*kapBint(i,j,id)
              kapdBint(i,j,id)=dx*kapdBint(i,j,id)
              dBint(i,j,id)=dx*dBint(i,j,id)
              kappar(i,id)=dx*chidBint
           else
              kapBint(i,j,id)=dx*kapBint(i,j,id)+kapBint(i,j-1,id)
              kapdBint(i,j,id)=dx*kapdBint(i,j,id)+kapdBint(i,j-1,id)
              dBint(i,j,id)=dx*dBint(i,j,id)+dBint(i,j-1,id)
              kappar(i,id)=dx*chidBint+kappar(i,id)
           end if

        end do

        kappar(i,id)=1.d0/kappar(i,id)
        kappap(i,id)=kapBint(i,nnu,id) ! + integral to infinity

        do j=1,nnu-1
           kapdBint(i,j,id)=kapdBint(i,j,id)/kapdBint(i,nnu,id)
           ! BAW, 20090317, normalize kapBint also
           kapBint(i,j,id)=kapBint(i,j,id)/kapBint(i,nnu,id)
           dBint(i,j,id)=dBint(i,j,id)/dBint(i,nnu,id)
        end do
        kapdBint(i,nnu,id)=1.d0
        ! BAW 20090317, normalize kapBint also
        kapBint(i,nnu,id)=1.d0
        dBint(i,nnu,id)=1.d0

     end do
  end do

  return
end subroutine setemiss


! ********************************************************************
! Function to Calculate LTE Emissivity j_nu=kappa_nu*B_nu
! ********************************************************************
real*8 function kappaB(f)

  use kappa_B_int
  implicit none

  real*8 :: f,nu
  real*8 :: Bnu

  nu=nu0+f*dnu
  kappaB=((1.d0-f)*kappa0+f*kappa1)*Bnu(nu,Teff)

  return
end function kappaB


! ********************************************************************
! Function to Calculatekappa_nu*dB_nu/dT
! ********************************************************************
real*8 function kappadB(f)

  use kappa_B_int
  implicit none

  real*8 :: f,nu
  real*8 :: dBnu

  nu=nu0+f*dnu
  kappadB=((1.d0-f)*kappa0+f*kappa1)*dBnu(nu,Teff)

  return
end function kappadB

! ********************************************************************
! Function to Calculate (1/chi_nu)*dB_nu/dT
! ********************************************************************
real*8 function chidB(f)

  use kappa_B_int
  implicit none

  real*8 :: f,nu
  real*8 :: dBnu

  nu=nu0+f*dnu
  chidB=dBnu(nu,Teff)/((1.d0-f)*chi0+f*chi1)

  return
end function chidB

! ********************************************************************
! Function to Calculate dB_nu/dT for integration
! ********************************************************************

real*8 function dB(f)

  use kappa_B_int
  implicit none

  real*8 :: f,nu
  real*8 :: dBnu

  nu=nu0+f*dnu
  dB=dBnu(nu,Teff)

  return
end function dB


real*8 function dBBfreq(T0,idust2)

  use dust_mod
  use random
  implicit none

  real*8 :: T0,lT0

  integer :: j,i0,i1,i,idust2
  real*8 :: xi,f,P0,P1,P

  lT0=log(T0)

  xi=ran()

  j=int(dble(nT-1)*(lT0-lTint(1))/(lTint(nT)-lTint(1)))+1

  if ( (j .lt. 1) .or. (j .gt. (nT-1)) ) then
     write(*,*) 'dBBfreq: index j out of bounds'
     print*,'T0',T0*11605.d0
     stop
  end if

  f=(lT0-lTint(j))/(lTint(j+1)-lTint(j))
  if ((f.lt.-0.0001d0).or.(f.ge.1.0001d0)) then
     write(*,*) 'dBBfreq: Tf = ',f,' out of bounds'
     write(*,*) '          T0  = ',T0*11606.d0
  end if

  i0=1
  i1=nnu

  P0=  (1.d0-f)*dBint(j,i0,idust2) + f*dBint(j+1,i0,idust2)
  if(xi.lt.P0) then
     if (P0 .gt. 0.d0) then
        dBBfreq=T0*xi*nuint(1)/P0
     else
        dBBfreq=T0*nuint(1)
     end if
     return
  end if

  P1=  (1.d0-f)*dBint(j,i1,idust2) + f*dBint(j+1,i1,idust2)
  if (P1 .gt. 1.d0) P1=1.d0

  if (xi.ge.P1) then
     if (P1 .lt. 1.d0) then
        if (xi.lt.1.d0) then
           dBBfreq=T0*nuint(nnu)*(1.d0-P1)/(1.d0-xi)
        else
           dBBfreq=T0*nuint(nnu)
        end if
     else
        dBBfreq=T0*nuint(nnu)
     end if
     return
  end if

  do while((i1-i0).gt.1)

     i=(i0+i1)/2
     P=  (1.d0-f)*dBint(j,i,idust2) + f*dBint(j+1,i,idust2)

     if (P .gt. 1.d0) then
        if (i.ge.nnu-2) then
           write(*,*)'dBBfreq: P>1; i,P = ',i,P
        end if
        P=1.d0
     end if

     if (P .lt. P0) then
        write(*,*) 'dBBfreq: P<P0; P0,P,P1 = ',P0,P,P1
        P0=P
        i0=i
     end if
     if (P .gt. P1) then
        write(*,*) 'dBBfreq: P>P1; P0,P,P1 = ',P0,P,P1
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
     write(*,*) 'dBBfreq: P1<=P0; P0, P1 = ',P0,P1
     stop
  end if

  if (xi.lt.P0) then
     write(*,*) 'dBBfreq: xi<P0, xi,P0 = ',xi,P0
     f=0.d0
  else if (xi.gt.P1) then
     write(*,*) 'dBBfreq: xi>P1, xi,P1 = ',xi,P1
     f=1.d0
  else if (P0.le.0.d0) then
     write(*,*) 'dBBfreq: P0<=0, P0,i0,i1 = ',P0,i0,i1
     f=0.5d0
  else if (P1.le.0.d0) then
     write(*,*) 'dBBfreq: P1<=0, P1 = ',P1
     stop
  else if (xi.le.0.d0) then
     write(*,*) 'dBBfreq: xi<=0, xi = ',xi
     stop
  else
     f=log(xi/P0)/log(P1/P0)
  end if

  dBBfreq= T0 * nuint(i0) * (nuint(i1)/nuint(i0))**f

  if ((f.lt.-0.0001d0).or.(f.ge.1.0001d0)) then
     write(*,*) 'dBBfreq: Pf = ',f,' out of bounds'
     write(*,*) '          xi = ',xi,'freq = ',1.24d0/dBBfreq
  end if

  return

end function dBBfreq



