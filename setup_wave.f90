! *******************************************************

subroutine setup_wave()
  !
  ! history:
  ! 00/03/19 (mjw):  set vger position and outside source using
  ! values from parameter file (e.g. REAPAR)
  ! 00/07/27 (mjw):  fix formatting issue for lf95/gcc (a**-1.5 type)
  !
  !

  use tts_mod
  use grid_mod
  use stokesvar_mod
  use taunum_mod
  use log_mod
  use opacin_mod
  use out_mod
  use dust_mod
  use random
  use tabl_mod
  use constants

  implicit none

  real*8 :: const,a1,c1,c0
  real*8 :: z1cgs,rmincgs,rmindcgs(ndg),rmaxdcgs(ndg),rddustcgs(ndg),p
  real*8 :: hnoam,taunoam,c0r2h
  real*8 :: rho1,sigz1,gmax,rhodust,tauzdust
  real*8 :: zmr

  integer :: i,id

  call diskdat_section('setup_wave')

  ! for calculating density in cgs
  rmincgs=rstar*rsol
  rmindcgs=rmind*rmincgs
  rmaxdcgs=rmaxd*rmincgs
  rddustcgs=rddust*rmincgs

  ! scaling factor for computing scattered light in peeled-off image
  sfactor=1.d0/(4.d0*pi)

  ! calculate quantities for disk midplane dust properties, id=1
  do id=1,ndg

     if(is_disk(id)) then
     print*,''
     print*,'dust grain index',id
     kapd(id)=kappav(id)
     z1cgs=rstar*rsol*z1(id)


     gmax=0.d0
     do i=1,nlambda(id)
        if (gdust(i,id).gt.gmax) gmax=gdust(i,id)
     end do
     ! peak=1.5*(1.d0-gmax**2)/(2.+gmax**2)*2./(1.d0+gmax**2-2.*gmax)**1.5

     ! calculate density from massd=mass of disk.
     ! = 2*pi*(integral(density*r*dr*dz)).
     a1=1.d0-a(id)
     p=b(id)-a(id)


     ! 11/19/02 BAW, change density formula to Whitney et al. 2002
     if(rmaxd(id).gt.0.d0.and.fmass(id).gt.0.d0) then
        ! c stays same
        const=r2p**1.5d0*rmincgs**2*z1cgs
        if(p.ne.-2.d0.and.p.ne.-1.5d0) then
           ! rho=rho0(id) at r=rmin but mass is integrated from rmind, hence
           ! print*,'yo!'
           ! print*,massd,msol,c,p,rmaxd,rmind
           rho0(id)=massd*msol/const/ &
                & (1.d0/(p+2.d0)*(rmaxd(id)**(p+2.d0)-rmind(id)**(p+2.d0)) &
                & -1.d0/(p+1.5d0)*(rmaxd(id)**(p+1.5d0)-rmind(id)**(p+1.5d0)))
           ! print*,rho0
           ! stop
        else if (p.eq.-2.d0) then
           rho0(id)=massd*msol/const/ &
                & (log(rmaxd(id)/rmind(id)) &
                & -1.d0/(p+1.5d0)*(rmaxd(id)**(p+1.5d0)-rmind(id)**(p+1.5d0)))
        else if (p.eq.-1.5d0) then
           rho0(id)=massd*msol/const/ &
                & (1.d0/(p+2.d0)*(rmaxd(id)**(p+2.d0)-rmind(id)**(p+2.d0)) &
                & -log(rmaxd(id)/rmind(id)))
        end if
        rho0(id)=rho0(id)*fmass(id)
        write(6,*)'rho0 of disk',rho0(id)
        ! opacity in units of rstar**-1 is kappa*rho*rmin.  c0 is multiplier
        ! used in radiative transfer.
        c0=kapd(id)*rho0(id)*rmincgs
        c1=kapd(id)*rho0(id)*rmincgs !for optical depth calcs here
        print*,'kapd',kapd(id)
        if(a(id).ne.1.0d0.and.a(id).ne.0.5d0) then
           taur(id)=c1/a1*(rmaxd(id)**a1-rddust(id)**a1) &
                & -c1/(0.5d0-a(id))*(rmaxd(id)**(0.5d0-a(id))- &
                & rddust(id)**(0.5d0-a(id)))
        else if (a(id).eq.1.0d0) then
           print*,'hi'
           taur(id)=c1*log(rmaxd(id)/rddust(id)) &
                & -c1/(0.5d0-a(id))*(rmaxd(id)**(0.5d0-a(id))- &
                & rddust(id)**(0.5d0-a(id)))
        else if (a(id).eq.0.5d0) then
           taur(id)=c1/a1*(rmaxd(id)**a1-rddust(id)**a1) &
                & -c1*log(rmaxd(id)/rddust(id))
        end if
     else
        const=0.d0
        c0=0.d0
        c1=0.d0
        rho0(id)=0.d0
        taur(id)=0.d0
     end if
     ! rhod0=rho0
     write(6,*)'c0,c1,taur',c0,c1,taur(id)
     ! stop
     if(rmaxd(id).gt.0.d0.and.fmass(id).gt.0.d0) then
        tauzmin=kapd(id)*sqrt(r2p)*z1cgs*rho0(id)*rddust(id)**p &
             & *(1.d0-sqrt(1.d0/rmind(id)))
        tauzdust=kapd(id)*sqrt(r2p)*z1cgs*rho0(id)*rddust(id)**p &
             & *(1.d0-sqrt(1.d0/rddust(id)))
        tauzmax=kapd(id)*sqrt(r2p)*z1cgs*rho0(id)*rmaxd(id)**p &
             & *(1.d0-sqrt(1.d0/rmaxd(id)))
        tauz30=kapd(id)*sqrt(r2p)*z1cgs*rho0(id)*(30.d0*autors)**p &
             & *(1.d0-sqrt(1.d0/(30.d0*autors)))
        write(6,*) 'tauzmin,tauzdust,tauzmax,tauz30' &
             & ,tauzmin,tauzdust,tauzmax,tauz30             
        ! density at rmind,30AU, rmaxd
        !  huh?  I don't think this is right. depends on b, no?  20111203.  fortunately it's not used, but still...don't print this out
  !      rhomin=rho0(id)*(rmind/rmin)**(-a(id))* &
  !           & (1.d0-sqrt(rmin/rmind))
  !      rhodust=rho0(id)*(rddust/rmin)**(-a(id))* &
  !           & (1.d0-sqrt(1.d0/rddust))
  !      rho30=rho0(id)*(30.d0*autors)**(-a(id))* &
  !           & (1.d0-sqrt(1.d0/(30.d0*autors)))
  !      rhomax=rho0(id)*(rmaxd/rmin)**(-a(id))* &
  !           & (1.d0-sqrt(rmin/rmaxd))
             rhomin=rho0(id)*(rmind(id)/rmin)**(-a(id))
             rhodust=rho0(id)*(rddust(id)/rmin)**(-a(id))
             rho30=rho0(id)*(30.d0*autors)**(-a(id))
             rhomax=rho0(id)*(rmaxd(id)/rmin)**(-a(id))
        ! write(6,*)'density at rmind,rddust,30au,rmaxd', &
        !     & rhomin,rhodust,rho30,rhomax
        sigz1=sqrt(r2p)*z1cgs*rho0(id)*(rmind(id)/rmin)**p* &
             & (1.d0-sqrt(rmin/rmind(id)))
        write(6,*) 'sigma of disk at rmind',sigz1
        ! density at 1 AU, sigma at .1 AU
        rho1=rho0(id)*(1.d0*autors/rmin)**(-a(id))* &
             & (1.d0-sqrt(1.d0/autors))
        sigz1=sqrt(r2p)*z1cgs*rho0(id)*(0.1d0*autors)**p* &
             & (1.d0-sqrt(0.1d0/autors))
        print*,'Sigma at 0.1 AU ',sigz1
        ! scale height at 100 AU
        zmr=z1(id)*(100.d0*autors/rmin)**b(id)
        print*,'h/r at 100 AU ',zmr/autors/100.d0
        ! scale height at rmaxd
        zmr=z1(id)*(rmaxd(id)/rmin)**b(id)
        print*,'h/r at rmaxd ',zmr/rmaxd(id)
        ! scale height at rmind
        zmr=z1(id)*(rddust(id)/rmin)**b(id)
        print*,'h/r at rddustd, h at rddust ',zmr/rddust(id),zmr
        print*,'id,z1,rddust,rmin,b',id,z1(id),rddust(id),rmin,b(id)
     end if

     ! make a file of disk height and constants vs radius for use with
     ! error function to calculate taus.
     ! if(rmaxd.gt.0.d0) then
     ! open(unit=14,file='tau.dat',status='unknown')
     ! write(14,*)
     ! 1        'radius (r),scale height (h),c0*r**-a*2h,tau at 3h,F(x1)'
     ! dnoam=rmaxd/10.d0
     ! do inoam=1,10
     ! rnoam=(dble(inoam)-0.5d0)*dnoam
     ! hnoam=z1*(rnoam/rmin)**b
     ! c0r2h=c0*(rnoam/rmin)**(-a)*2.*hnoam*(1.d0-sqrt(rmin/rnoam))
     ! if(c0.gt.0.d0) then
     ! taunoam=c0r2h*0.0013d0
     ! xnoam=1.d0-1.d0/c0r2h
     ! end if
     ! write(14,*) rnoam,hnoam,c0r2h,taunoam,xnoam
     ! end do
     ! close(14)
     ! end if

     ! find tau for H=3.5 at rmin
     hnoam=z1(id)*rmind(id)**b(id)
!     c0r2h=c0*(rmind)**(-a(id))*2.*hnoam*(1.d0-sqrt(rmin/rmind))
     c0r2h=c0*(rmind(id))**(-a(id))*2.*hnoam
     if(c0.gt.0.d0) then
        taunoam=c0r2h*0.0002d0
     end if

     ! call diskdat_write('taunoam',taunoam,
     ! 1 'optical depth to H=3.5 at rmind')
     ! print*, 'optical depth to H=3.5 at rmind',taunoam


     call diskdat_write('b',b(id),'z=z1*r**b.')
     call diskdat_write('z1',z1(id),'')
     call diskdat_write('a',a(id),'rho(z=0)=c0*r**(-a)')
     ! call diskdat_write('c0',c0,'')
     !call diskdat_write('rhomin',rhomin,'density at rmind')
     !call diskdat_write('rho30',rho30,'density at 30 AU')
     !call diskdat_write('rhomax',rhomax,'density at rmaxd')
     ! call diskdat_write('oadeg',oadeg,'oa (tan(oa)=zmax/rmaxd)')
     call diskdat_write('massd',massd,'mass in disk')
     call diskdat_write('rho0(id)',rho0(id), &
          & 'density at midplane of disk at stellar radius')
     call diskdat_write('tauzmin',tauzmin,'tau-z at rmind')
     call diskdat_write('tauz30',tauz30,'tau-z at 30 AU')
     call diskdat_write('tauzmax',tauzmax,'tau-z at rmaxd')

     end if
  end do
  call diskdat_write('massd',massd,'mass in disk')
  print*,''


  return
end subroutine setup_wave

