subroutine gridset(iter,i1)

! set up grid:  3-D spherical grid with variable spacing in r and
! theta (set by exponents, rexp, texp)
! make a bunch of arrays necessary for find_wall subroutine

! history:
! 2000/09/06 (baw): write subroutine
! 2012/03/22 (mjw): replace random number stuff for fractal density
! 			including resetting generator to previous value
!			after using generator FRACTAL_GEN for the density
!			stuff
! 2012/03/28 (mjw): enabled old RNG use for fractal generation (sort of undo
! 			changes of 03/22 but leave code in, but commented out)
!

  use output_mod
  use log_mod

  use tts_mod
  use grid_mod
  use opacin_mod
  use out_mod
  use dust_mod
  use spot_mod
  use constants
  use random
  use wrimsg_mod

  implicit none

  real*8 :: dr,dt,dp,rad,phi,thet,densd,dense,vol,densd2,densd1
  real*8 :: cost,sint,eps,tau,rming1,tiny
  real*8 :: rmincgs,taud,taue,taufmin,mu,thetatm
  real*8 :: thetmin,sigmu,rsonr,Tdisk
  real*8 :: mudisk,a1,tauw,rtop,rfac,tauRr,taumid,tauRth
  real*8 :: res,thettop,maxdens,tauRmu,radamb,rbeg2
  real*8 :: tauave,lstar,lacc,lacc2,tave,Vc,lacc_missing
  real*8 :: rminsub,fluxratio
  real*8 :: fplay,rexpnew,ltot2,diskscale,rmintest
  real*8 :: dr3,dcost,r,masssg,x,y,z,masstot,mu0tmp
  real*8 :: aave,tmp,taup,dphi,totvol,kappafave
  real*8 :: fl,const,fm1,fm2,fm3,fm4,densav,densav2,rntot,std
  real*8 :: maxdensvar,mindensvar,mdottotcgs,rhoave,mdotratio,tant
  real*8 :: adens,bdens,zscale,rhodens,sigpw,sigpwtot,sigtot,sigscale
  real*8 :: delz,tanb,rcyl
  real*8 :: tauzrat,taudisk,tauRos,tauv,mu0
  real*8 :: sigtmp(2),r1,r2,th1,th2,tmpmass(8),diskscl(8)
  real*8 :: inc,cosi,sini,x1,y1,z1hi,rsq1,cosb1,ph1,dens1au,rdduste

  ! real :: xyarr(200,200),x,y,xmax,dx,r,dens1,z
  real*8 :: chiR,ierfc,erfc,m
  integer :: ir,it,ip,nttmp,nr0,ntatm,nrwall,irbeg,itmid,ircyl,ir1,it1,ip1
  integer :: dcount,id,i,iter,igas,inhole,j,i1,icent,nf,pexp,iskip,irgap

  integer :: peelid
  character :: cmsgnm*70

  call diskdat_section('gridset')

  tiny=1.d-15
  rmincgs=rstar*rsol
  radamb=rmax
  iskip=0

  ! I have some gridding below in case there's a gas disk, but that's a future update
  ! for now, set gas disk to no.
  igas=0

  print*,'idiskacc',idiskacc

  if (idiskacc.eq.1) then

     ! set disk accretion fraction based on alpha disk paramters
     ! see eqns 5 and 6 in Whitney et al. 2003, apj, 591, 1049
       lstar=4.d0*pi*sigt*tstar**4*rmincgs**2
     ! 2011 july 6. modification, only counting non-hotspot photons
   !  if (fspot.lt.1.d0) then
     !  lstar=4.d0*pi*sigt*tstar**4*rmincgs**2*(1-fspot)
   !  else
  
     if (ialpha.eq.1) then
        ! if specifying alphad in input instead of mdotdisk
        print*,'calculating disk accretion rate from alphad)'
        print*,'ialpha,alphad',ialpha,alphad
        Vc=dsqrt(gn*msol*massc/rmincgs)
        ! I think this will work...20090731
        mdot=0.d0
        do id=1,ndg
           if(is_disk(id)) then
              mdotcgs(id)=dsqrt(18.d0*pi**3.d0)*alphad*Vc*rho0(id)*(z1(id)*rmincgs)**3.d0/rmincgs
              mdot=mdot+mdotcgs(id)
              ! units of cgs
           end if
        end do
        mdotdisk=mdot*3600.d0*24.d0*365.25d0/msol   !units of Msol/yr
     else
        ! if specifying mdotdisk in input
        Vc=dsqrt(gn*msol*massc/rmincgs)
        ! mdot=mdotdisk/3600.d0/24.d0/365.25d0*msol
        mdottotcgs=mdotdisk/(3600.d0*24.d0*365.25d0/msol)
        mdotratio=rho0(1)*z1(1)**3/rho0(2)/z1(2)**3
        mdotcgs=0.
        mdotcgs(2)=mdottotcgs/(1.d0+mdotratio)
        mdotcgs(1)=mdottotcgs-mdotcgs(2)
        mdotcgs(3:6)=0.d0
        print*,'rho0',rho0
        print*,'z1',z1
        do id=1,ndg	      
           if(is_disk(id)) then
             if (mdotcgs(id).gt.0.) alphad=mdotcgs(id)/(dsqrt(18.d0*pi**3.d0)*Vc*rho0(id)*(z1(id)*rmincgs)**3.d0/rmincgs)
             print*,'alphad,id',alphad,id
           end if
        enddo
     end if
     print*,'alphad',alphad
     print*,'disk accretion rate',mdotcgs*3600.d0*24.d0*365.25d0/msol
     print*,'total disk accretion rate, input',sum(mdotcgs)*3600.d0*24.d0*365.25d0/msol,mdotdisk
     mdotdisk=sum(mdotcgs)*3600.d0*24.d0*365.25d0/msol
  !   stop

     ! accretion in disk
     lacc=gn*(msol*massc)*sum(mdotcgs)/2.d0/(rddustmin*rmincgs)

     ! accretion in shock, the optically thick part from the heated
     ! atmosphere  (calvet & gullbring 1998).
     ! we will split it into half blackbody at Tshock, and half x-rays.
     ! assumes co-Rotation radius is 5 Rstar.
     if (massd.gt.0.d0) then
        lacc2=gn*(msol*massc)*sum(mdotcgs)*(1.d0/rmincgs-1.d0/(rmincgs*rtrunc))
        if (rtrunc<rddustmin) &
		! missing accretion
		   &  lacc_missing=gn*(msol*massc)*sum(mdotcgs)/2.d0/(rtrunc*rmincgs)  &
		    &             -gn*(msol*massc)*sum(mdotcgs)/2.d0/(rddustmin*rmincgs)                 
     else
        lacc2=0.d0
        lacc_missing=0.
        lacc=0.
     end if

          
     ! total luminosity
     ltot=lacc+lacc2+lstar+l_isrf
     print*,'lacc,lacc2,lstar,l_isrf',lacc/lsun,lacc2/lsun,lstar/lsun,l_isrf/lsun
!     print*,'ltot',ltot/lsun
     print*,'missing accretion between rtrunc and rddust',lacc_missing/lsun


     accfrac=lacc/ltot
     accfrac2=lacc2/ltot
     isrf_frac=l_isrf/ltot
 !    mdot=sum(mdotcgs)*3600.d0*24.d0*365.25d0/msol
     print*, 'including disk accretion luminosity '
     print*, 'disk accretion rate (msol/yr) ',mdotdisk
     print*, 'fraction of acc luminosity from disk ',accfrac
     print*, 'fraction of lum. on stellar hotspot ',accfrac2
     print*, 'total system luminosity (Lsun)',ltot/lsun

     call diskdat_write('isrf_frac',isrf_frac, 'fractional ISRF lum')
     call diskdat_write('diskacc',.true.,'whether disk luminosity is included')
     call diskdat_write('mdot',mdotdisk,'disk accretion rate (msol/yr)')
     call diskdat_write('accfrac',accfrac,'fraction of lum. from disk')
     call diskdat_write('accfrac2',accfrac2,'fraction of lum. on stellar hotspot')
     call diskdat_write('ltot',ltot/lsun,'total system lum. (Lsun)')
     call diskdat_write('lmissing',lacc_missing/lsun,'missing acc. lum. bet. rtrunc and rddust')

     ! f=calvet's filling factor, ranges from 0.01 to 0.001
     ! median is 0.007
     ! if (fspot.gt.0.d0.and.lacc2.gt.0.d0) then  !eek, 20130612
	 if (lacc2.gt.0.d0) then
	    if (itemp.eq.0) then
       	  fluxratio=0.5d0*lacc2/lstar/fspot
       	  tshock=tstar*(1.d0+fluxratio)**0.25d0
       	  tspot=tshock
        else
	      tshock=tspot
		  fluxratio=(tshock/tstar)**4-1.d0
		  fspot=0.5d0*lacc2/lstar/fluxratio
 	    endif
        if (nspot.eq.1) then
           thspot=acos(1.d0-2.d0*fspot)*rad2deg
        else if (nspot.eq.2) then
           thspot=acos(1.d0-fspot)*rad2deg
        else
           print*,'error!  nspot not set,stopping program'
           stop
        end if

        print*,'fspot',fspot
        print*,'number and size of spots ',nspot,thspot
        print*,'thermal shock (spot) T ',tshock
 !       stop

        call diskdat_write('fspot',fspot,'fractional stellar hotspot(s) size')
        call diskdat_write('tshock',tshock,'thermal shock (spot) T')

        fplay=0.5d0*lacc2/fspot/4.d0/pi/rmincgs**2
        print*,'log flux shock ',log10(fplay)
     else
        thspot = 0.d0
     end if
     
     if(partial_peeloff) call output_accretion(lacc, lacc2)

  else

     lstar=4.d0*pi*sigt*tstar**4*rmincgs**2
     ltot=lstar+l_isrf
     isrf_frac=l_isrf/ltot
     mdot=0.d0
     accfrac=0.d0
     lacc=0.d0
     lacc2=0.d0

     call diskdat_write('diskacc',.false., 'whether disk luminosity is included')
     call diskdat_write('ltot',lstar/lsun, 'total stellar lum. (Lsun)')
     call diskdat_write('ltot',ltot/lsun, 'total system lum. (Lsun)')
     call diskdat_write('isrf_frac',isrf_frac,'fractional ISRF lum')

     print*, 'disk accretion rate (msol/yr) ',mdot
     print*, 'fraction of luminosity from disk ',accfrac
     print*,'lstar/lsun,l_isrf/lsun',lstar/lsun,l_isrf/lsun
     print*, 'total stellar luminosity (Lsun)',lstar/lsun
     print*, 'total system luminosity (Lsun)',ltot/lsun
     print*, 'isrf_frac',l_isrf/ltot

  end if
  
  !BAW TESTING!!!!   2010 apr 28.   delete all these lines after testing
!  lacc2=0.5*lstar
!  accfrac2=lacc2/ltot
!  tspot=tstar*2.5
!  tshock=tspot
!  ! total luminosity
!  ltot=ltot+lacc2
!  print*,'lacc,lacc2,lstar,l_isrf',lacc,lacc2,lstar,l_isrf
!  print*,'ltot',ltot
!  if (nspot.eq.1) then
!     thspot=acos(1.d0-2.d0*fspot)*rad2deg
!  else if (nspot.eq.2) then
!     thspot=acos(1.d0-fspot)*rad2deg
!  else
!     print*,'error!  nspot not set,stopping program'
!     stop
!  end if
!  print*,'number and size of spots ',nspot,thspot
!  print*,'calculated thermal shock T ',tshock
!!  stop 

  ! calculate an average stellar temp based on luminosity
  ! of star+accretion
  tave=tstar*(1+0.5*lacc2/lstar)**0.25
  print*,'average stellar temp is now',tave

  print*,'tshock,tave',tshock,tave
  print*,'note: should be same if fspot=1'
 ! if (fspot.eq.1) tshock=tave

  if (fspot.eq.1.0.and.iplanckst.eq.0) then
    print*,'WARNING: due to accretion, emitting from tstar = ',tave
    print*, 'you may want to change stellar atmosphere file to agree with this!'
  endif

  rminsub=(tsub/tave)**(-2.1d0)
  call diskdat_write('Tave',Tave,'averave Tstar including hotspot')

  ! reset tstar to tave, since it is used for temperature
  ! calculation (sets the luminosity scale).
  tstarave=tave
  ! tstar=tave

  if ((rminsub.gt.rmaxd(1).or.rminsub.gt.rmaxd(2)).and.massd.gt.0.0d0) then
     print*,''
     print*,'****************'
     print*,'ERROR.  dust sublimation radius larger than RMAXD'
     print*,'probably because you chose a hot star and'
     print*,'dust sublimation radius is large.'
     print*,'Increase RMAXD.'
     print*,'RMIND, RMAXD in AU',rminsub/autors,rmaxd(1)/autors,rmaxd(2)/autors
     print*,'Stopping program'
     stop
  end if

  if (irminsub.eq.1) then
     rmine=rminsub*rmine_in
     rmind=rminsub*rmind_in
     rddust=rminsub*rmind_in
    ! rmind2=rmind2*rminsub
     rmin_sg=rminsub*rmin_sg
     print*,'resetting Rsub to ',rminsub
     print*,'inner envelope radius ',rmine
     print*,'inner disk radius ',rmind
     call diskdat_write('rdust',rmine,'updated dust destruction radius')
     call diskdat_write('rmine',rmine,'updated inner envelope radius')
     call diskdat_write('rmind',rmind(1),'updated inner disk radius')
     call diskdat_write('rmind',rmind(2),'updated inner disk radius')
  end if

  if (rddust(1).gt.rmaxd(1)) then
     print*,'ERROR, rddust > rmaxd'
     print*,'rddust,rmaxd (AU)',rddust(1)/autors,rmaxd(1)/autors
     print*,'stopping program'
     stop
  end if

if (rddust(2).gt.rmaxd(2)) then
   print*,'ERROR, rddust > rmaxd'
   print*,'rddust,rmaxd (AU)',rddust(1)/autors,rmaxd(1)/autors
   print*,'stopping program'
   stop
end if

  ! calculate scale height based on Tsub at Rsub
  ! can't do this because we already set disk parameters, would
  ! have to iterate.  save this for hseq code.
  ! testing shows it only increases scale height slightly.
  ! k=1.38e-16
  ! muc=2.3
  ! mH=1.67e-24
  ! rcgs=rminsub*rstar*rsol
  ! mcgs=massc*msol
  ! honr=sqrt(k*Tsub/Gn/mcgs*rcgs/muc/mH)
  ! h=honr*rminsub
  ! print*,'h on r at Rsub',honr
  ! print*,'h (rsub) ',h
  ! zmin=honr*rminsub/rminsub**b
  ! print*,'h_0 based on h(r_sub)',zmin
  ! print*,'assumes Tsub=',Tsub,'    beta =',b

  ! convert opacity to units of cm^2/gm*rstar(cgs) because distance
  ! is in units if 1/rstar, and dtau=kapd*rho*ds

  do id=1,ndg
     kapd(id)=kappav(id)*rmincgs
  end do

  ! 20070119 BAW, only recalcalate rminsub after first iteration
  ! note we aren't doing this since we aren't calling this subroutine after the first iteration
  if (iter.gt.1) then
     it=(ntg+1)/2
     ir=3
     do while (tdust2(ir,it,1,1).ge.1600.d0*11605.d0)
        ir=ir+1
        print*,tdust2(ir,it,1,1)*11605.d0
     end do
     rmintest=rarr(ir-1)

     ! ip=1
     ! ir=1
     ! rmintest=0.d0
     ! do it=1,ntg
     ! do while (densarr(ir,it,ip).eq.0.d0)
     ! ir=ir+1
     ! if(rarr(ir).gt.rmintest) rmintest=rarr(ir)
     ! end do
     ! end do

     ! new rminsub
     rminsub=rmintest
     if ((rminsub.gt.rmaxd(1).or.rminsub.gt.rmaxd(2)).and.massd.gt.0.0d0) then
        print*,''
        print*,'****************'
        print*,'ERROR.  dust sublimation radius larger than RMAXD'
        print*,'probably because you chose a hot star and'
        print*,'dust sublimation radius is large.'
        print*,'Increase RMAXD.'
        print*,'RMIND, RMAXD in AU',rminsub/autors,rmaxd(1)/autors
        print*,'RMIND, RMAXD in AU',rminsub/autors,rmaxd(2)/autors
        print*,'Stopping program'
        stop
     end if

     if (irminsub.eq.1) then
        rmine=rminsub*rmine_in
        rmind=rminsub*rmind_in
        rddust=rminsub*rmind_in
        print*,'resetting Rsub to ',rminsub
        print*,'inner envelope radius ',rmine
        print*,'inner disk radius ',rmind
        write(12,*) rmine, ' updated dust destruction radius'
        write(12,*) rmine, ' updated inner envelope radius'
        write(12,*) rmind, ' updated inner disk radius'
     end if
  end if


  if (rddust(1).gt.rmaxd(1)) then
     print*,'ERROR, rddust > rmaxd(1)'
     print*,'rddust,rmaxd (AU)',rddust(1)/autors,rmaxd(1)/autors
     print*,'stopping program'
     stop
  end if
  if (rddust(2).gt.rmaxd(2)) then
     print*,'ERROR, rddust > rmaxd(2)'
     print*,'rddust,rmaxd (AU)',rddust(2)/autors,rmaxd(2)/autors
     print*,'stopping program'
     stop
  end if

  rddustmin=min(rddust(1),rddust(2))


  ! calculate some constants for the TSC envelope density
  call envset()

  ! make grid include minimum and maximum values.

  print*, ' '
  print*, 'grid setup'

! set up rarr(nrg), radial grid
  print*,'rddust,rmine,rmind,rmin',rddust,rmine,rmind,rmin
  ! rgrid
  ! 20050628, change rexp for disks with large inner holes
  if (rddust(1).gt.0.2d0*rmaxd(1).and.rddust(2).gt.0.2d0*rmaxd(2)) then
     rexpnew=2.d0
     print*,'disk has large inner hole'
     print*,'changing radial density exponent to rexp=2'
  else
     rexpnew=rexp
  end if
  ! gflag=0   !commented out 20090731
  rarr(1)=rmin

  if (massd.eq.0.d0) then
  !   rarr(2)=rddust
  !   dr=(rmax-rmine)/(dble(nrg-2))**rexpnew
  !   print*,'rexp,dr,rmin',rexpnew,dr/autors,rmin/autors
  !   do ir=3,nrg
  !      rarr(ir)=rmine+dr*(dble(ir-2))**rexpnew
  !   end do

     dr=(rmax-rmine)/(dble(nrg-1))**rexpnew
     print*,'rexp,dr,rmin',rexpnew,dr/autors,rmin/autors
     do ir=2,nrg
       rarr(ir)=rmine+dr*(dble(ir-1))**rexpnew
     end do

  else

     irbeg=2
   !  rarr(irbeg)=rddust
     rdduste=min(rddust(1),rddust(2),rmine)
     rarr(irbeg)=rdduste
          
     !I'm going to try log spacing through entire r grid.  !3/19/12   got rid of old version with log spacing in center, power law outside
     rfac=(rmax/rdduste)**(1./dble(nrg-2))
     print*,'rfac',rfac
     do ir=2,nrg
       rarr(ir)=rdduste*rfac**(ir-2)
     !  print*,'rarr',rarr
     enddo
     print*,'rarr(irbeg),rarr(nrg),rddust,rmax',rarr(irbeg),rarr(nrg),rddust,rmax
     
  end if

  ! set rmine to a grid location
  ! if (gflag.eq.1) then
  ! changed condition, 20090731
  if (rddustmin.lt.rmine) then
     call locate(rarr,nrg,rmine,ir)
     rarr(ir+1)=rmine
  end if

! same with rddust
if (rmine.lt.rddustmin) then
   call locate(rarr,nrg,rddustmin,ir)
   rarr(ir+1)=rddustmin
end if


  ! set rgapd1 and rgapd2 to a grid location
  if (igapd.eq.1) then
	if (rgapd1.gt.rmin.and.rgapd1.lt.rmax) then
	  call locate(rarr,nrg,rgapd1,ir)
	  rarr(ir+1)=rgapd1
    endif
    if (rgapd2.gt.rmin.and.rgapd2.lt.rmax) then
	  call locate(rarr,nrg,rgapd2,ir)
	  rarr(ir+1)=rgapd2
    endif
  endif

  ! r-squared array
  do ir=1,nrg
     r2arr(ir)=rarr(ir)**2.d0
     ! print first few points of rarr
     if (ir.lt.10) then
        print*,'rarr/rmin,ir ',rarr(ir)/rmin,ir
     end if
  end do

  ! r-ave array
  do ir=1,nrg-1
     r1=log10(rarr(ir))
     r2=log10(rarr(ir+1))
     ravearr(ir)=10**(0.5d0*(r1+r2))
  !   ravearr(ir)=0.5d0*(rarr(ir)+rarr(ir+1))
  end do

  open(unit=15,file='rwalls.dat',status='unknown')
  write(15,*) nrg,' = number of grid points in r'
  write(15,*) &
       & '       index      r/rstar         r(rsub)         r(au)'
  do ir=1,nrg
     write(15,900) ir,rarr(ir), &
          & rarr(ir)/rminsub, &
          & rarr(ir)/autors
  end do
  close(15)
900 format(i10,3(1x,f15.5))

  open(unit=15,file='rave.dat',status='unknown')
  write(15,*) nrg-1,' = number of r cells'
  write(15,*) &
       & '       index      r/rstar         r(rsub)         r(au)'
  do ir=1,nrg-1
     write(15,920) ir,ravearr(ir), &
          & ravearr(ir)/rminsub, &
          & ravearr(ir)/autors
  end do
  close(15)
920 format(i10,3(1x,f15.5))

  print*,'rarr(nrg),r2arr(nrg)',rarr(nrg),r2arr(nrg)

  !  theta array

    nttmp=(ntg)/2

!!  log spacing
!     tmptharr(1)=pi/dble(nttmp)
!     m=(pihalf/tmptharr(1))**(1.d0/dble(nttmp-1.d0))
!      do it=2,nttmp
!        tmptharr(it)=tmptharr(it-1)*m
!     enddo
!  !  redo with new value for tmptharr(1)
!     tmptharr(1)=(tmptharr(2)-tmptharr(1))
!     m=(pihalf/tmptharr(1))**(1.d0/dble(nttmp-1.d0))
!      do it=2,nttmp
!        tmptharr(it)=tmptharr(it-1)*m
!     enddo

 !  power law
    pexp=2
    m=pi/2.d0/dble(nttmp)**pexp
    do it=1,nttmp
	  tmptharr(it)=m*dble(it)**pexp
    enddo

   do it=1,nttmp
	   thetarr(nttmp+1+it)=tmptharr(it)+pihalf
!     print*,'tharr ',thetarr(it)*rad2deg
   enddo
	 thetarr(nttmp+1)=pihalf
   do it=1,nttmp
 	   thetarr(it)=pi-thetarr(ntg-it+1)
!   print*,'tharr ',thetarr(it)*rad2deg
	 enddo
	
!  thetarr(1)=0.d0
!  thetarr(ntg)=pi

!  print*,''
!  do it=1,ntg
!    print*,'thetarr final, ',thetarr(it)*rad2deg
!  enddo

  do it=1,size(thete_arr)
     if(minval(abs(thete_arr(it) - thetarr))==0.) then
        print *,'ERROR: peeloff angle is exactly the same as one of the wall angles'
        print *,'       and this is not recommended. The angle is ',thete_arr(it) * rad2deg
        print *,'change thete to a slightly different value and restart!'
        stop
     end if
  end do

  ! this is not the eps used in the rest of the code!  see
  ! newdisktherm for that
  eps=1.d-8
  do it=1,ntg
     if (it.eq.1) then
        thetarr(it)=0.d0
        costarr(it)=1.d0
        sintarr(it)=0.d0
        tan2arr(it)=0.d0
     else if (it.eq.ntg) then
        thetarr(it)=pi
        costarr(it)=-1.d0
        sintarr(it)=0.d0
        tan2arr(it)=0.d0
     else if (it.eq.nttmp+1) then
        thetarr(it)=pihalf
        costarr(it)=0.d0
        sintarr(it)=1.d0
        tan2arr(it)=-1.d0
     else
        costarr(it)=cos(thetarr(it))
        sintarr(it)=sin(thetarr(it))
        tan2arr(it)=tan(thetarr(it))**2.d0
     end if
     ! print*,'thetarr,costarr,tan2arr '
     ! print*,'thetarr',thetarr(it)*rad2deg
  end do

  ! theta-ave array
  do it=1,ntg-1

     ! The following check is needed because
     ! ifort -m32 produces a discontinuity
     if(thetarr(it+1)<thetarr(it)) then
        stop "Error - thetarr is not monotonically increasing"
     end if

     th1=thetarr(it)**pexp
     th2=thetarr(it+1)**pexp

     thetavearr(it)=(0.5d0*(th1+th2))**(1./pexp)
 !    thetavearr(it)=0.5d0*(thetarr(it)+thetarr(it+1))

  end do

  open(unit=15,file='thwalls.dat',status='unknown')
  write(15,'(I6," / number of grid points in theta")') ntg
  write(15,'(" index    thet(rad)    thet(deg)      '// &
       & 'cost      tan**2(thet)")')
  do it=1,ntg
     write(15,'(i5,4(1x,f15.5))') it,thetarr(it), &
          & thetarr(it)*rad2deg,costarr(it),tan2arr(it)
  end do
  close(15)

  open(unit=15,file='thave.dat',status='unknown')
  write(15,'(I6," / number theta cells")') ntg-1
  write(15,'(" index    thet(rad)    thet(deg)      ")')
  do it=1,ntg-1
     write(15,'(i5,2(1x,f15.5))') it,thetavearr(it), &
          & thetavearr(it)*rad2deg
  end do
  close(15)
  
  dp=r2p/dble(npg-1)
  do ip=1,npg
     phiarr(ip)=dp*(dble(ip-1))
     aarr(ip)=dsin(phiarr(ip))
     barr(ip)=-dcos(phiarr(ip))
     ! print*,'phiarr, aarr, barr ',phiarr(ip),aarr(ip),barr(ip)
     ! carr(ip)=
     ! darr(ip)=
  end do

  ! phi-ave array
  do ip=1,npg-1
     phiavearr(ip)=0.5d0*(phiarr(ip)+phiarr(ip+1))
  end do

  open(unit=15,file='phiwalls.dat',status='unknown')
  write(15,'(I6," / number of grid points in phi")') npg
  write(15,'(" index    phi(rad)    phi(deg)")')
  do ip=1,npg
     write(15,'(i5,2(1x,f12.5))') ip,phiarr(ip), &
          & phiarr(ip)*rad2deg
  end do
  close(15)

  if(diffusion) then
     ! set up diffusion grid
     diffus = .false.
     diffdir = 0
     dcount=0
	if (massd.gt.0.d0) then
		!define some constants	
		mu0=0.d0
		tauv=0.d0
	    tauzrat=0.d0
	    tauRos=0.d0
		do id=1,2
          mu0tmp=fmass(id)*z1(id)*Rddust(id)**b(id)/Rddust(id)
          mu0=mu0+mu0tmp
          tauv=tauv+taur(id)
          a1=1.0d0-a(id)
          tauzrat=tauzrat+sqrt(pihalf)*a1*mu0tmp/ &
               & ((rmaxd(id)/rddust(id))**a1-1.d0)
          tauRos=tauRos+taur(id)*chiR(1600.d0/11605.d0,id)
	   enddo
	   tauRth=tauzrat*tauRos
	   taudisk=min(10.d0,tauRth)

        if (rddust(1).eq.rmin.and.rddust(2).eq.rmin) then
           nr0=1
        else
           nr0=3
           diffus(1,:,:) = .false.
           diffdir(1,:,:) = 0
        end if

        do ir=nr0,nrg-1

           r=ravearr(ir)/rddustmin
 		   if (r*rddustmin.gt.rmaxdmax) exit

             Rsonr=1.d0/(rddustmin*r)
             Tdisk=min(1600.d0/11605.d0,max(3.d0/11605.d0, &
                & (Tstar/11605.d0) &
                & *(max(2.d0/3.d0*(Rsonr)**3, &
                & (asin(Rsonr)-(Rsonr)*sqrt(1.d0-(Rsonr)**2)))/pi &
                & )**0.25d0))

	 	 ! print*,''
         ! print*,'r/rddust,tdisk,',r,Tdisk*11605.d0

           sigmu=0.d0
           taumid=0.d0

  !         do id=1,ndg
  !            if(is_disk(id)) then
  !!               sigmu = sigmu + fmass(id)*mu0*r**(b(id)-1.d0)
   !              taumid = taumid + fmass(id)*tauzrat*tauv*chiR(Tdisk,id)/r**(a(id)-b(id))
   !           end if
   !        end do

!  just define diffusion based on large-grain disk.  has to be within that disk to work.

                do id=1,1
                sigmu = sigmu + fmass(id)*mu0*r**(b(id)-1.d0)
                taumid = taumid + fmass(id)*tauzrat*tauv*chiR(Tdisk,id)/r**(a(id)-b(id))
                a1=1.0d0-a(id)
                enddo

 
!	print*,'sigmu,taumid,a1',sigmu,taumid,a1
!	print*,'taudisk',taudisk   !taudisk may need to be defined only with dust 1!

           if (taumid.gt.taudisk) then
              mudisk=sigmu*sqrt(2.d0)*ierfc(taudisk/taumid)
           else
              mudisk=0.d0
           end if

	 !     print*,'mudisk',mudisk

           do it=1,ntg-1

              mu=abs(cos(thetavearr(it)))
              tauRmu=taumid*erfc(mu/(sigmu*sqrt(2.0d0)))
              tauRr=tauRos*exp(-0.5d0*mu*mu/(sigmu*sigmu))* &
                   & (1.d0-r**a1)/(1.d0-(rmaxdmax/rddustmin)**a1)

	
	!		if (it.eq.ntg/2+1) then
	!			print*,'mu,mudisk,tauRr,taudisk,',mu,mudisk,tauRr,taudisk
	!		endif

              do ip=1,npg-1
	
	           if (igapd.eq.1.and.r*rddustmin.lt.rgapd2.and.r*rddustmin.gt.rgapd1) then
	!	         if (it.eq.1.and.ip.eq.1) print*,'in gap region, setting diffusion =false, r,r_au,ir',r*rddust,r*rddust/autors,ir
		         diffus(ir,it,ip)=.false.
		         diffdir(ir,it,ip)=0
		       else
                if ((mu.lt.mudisk).and.(tauRr.gt.taudisk)) then
  !               if ((mu.lt.0.25*mudisk).and.(tauRr.gt.4.*taudisk)) then
                    diffus(ir,it,ip) = .true.
                    dcount=dcount+1
                 else
                    diffus(ir,it,ip) = .false.
                 end if

    			  if (tauRr.lt.tauRmu) then
                    diffdir(ir,it,ip) = -1
                 else
                    if (thetarr(it).lt.pihalf) then
                       diffdir(ir,it,ip) = -2
                    else
                       diffdir(ir,it,ip) = 2
                    end if
                 end if

               endif

              end do
           end do
        end do
     end if
  else
     diffus = .false.
     diffdir = 0
     dcount = 0
  end if
  
  print*,'number of diffusion cells in grid',dcount
  call diskdat_write('dcount',dcount,'number of diffusion cells in grid')

 !stop

  !calculate fractal density variations
    if (ifractal.eq.1) then
      call diskdat_write('ifractal',ifractal,'adding fractal density variations')

      write(cmsgnm,'(a)')  'computing fractal density variations'
      call wrimsg('GRIDSET',cmsgnm)
  !    call get_current_generator(current_gen)
  !    call set_current_generator(fractal_gen)
  !    write(cmsgnm,'(a,i2,a,i2)')  'Switching from RNG stream ',current_gen,&
  !        ' to fractal generator ',fractal_gen
      write(cmsgnm,'(a)')  'Using old RNG for fractal density'
      call wrimsg('GRIDSET',cmsgnm)
	  call set_random_seed(ifseed)
      nf=32
      fl=fractl
  !   fl=3.792
      icent=0
      const=0.
      fm1=1.
      fm2=1.
      fm3=1.
      fm4=1.
 !   for more info, see ../fractal/fractal_sph.f
      call fractal_sph(nf,fl,const,fm1,fm2,fm3,fm4,icent)
  !modify densvararr to have average value of 1, std dev = densratio
	  densav=0.d0
	  densav2=0.d0
	  rntot=(nrg-1)*(ntg-1)*(npg-1)
      do ir=1,nrg-1
	  do it=1,ntg-1
	  do ip=1,npg-1
		 densav=densvararr(ir,it,ip)+densav
		 densav2=densvararr(ir,it,ip)**2+densav2
	  enddo
      enddo
      enddo
      print*,'rntot',rntot
      densav=densav/rntot
      densav2=densav2/rntot
      print*,'densav2,densav**2',densav2,densav**2
      print*,'std^2',densav2-densav**2
      std=sqrt(densav2-densav**2)
      print*,'densav,std',densav,std
      print*,'std/densav',std/densav
      print*,'std',std
  !    call diskdat_write('std/densav',std/densav,'std dev of dens. vars normalized to ave')
      do ir=1,nrg-1
    do it=1,ntg-1
    do ip=1,npg-1
!      densvararr(ir,it,ip)=densvararr(ir,it,ip)/densav
         densvararr(ir,it,ip)=densvararr(ir,it,ip)/densav*densratio+1-(densratio)
!   densvararr(ir,it,ip)=(densvararr(ir,it,ip)-densav)*densratio/std+1.
!   densvararr(ir,it,ip)=densvararr(ir,it,ip)*densratio/std
    enddo
    enddo
    enddo
   ! test; can comment this out when it works
  densav=0.d0
	  densav2=0.d0
	  mindensvar=0.d0
	  maxdensvar=0.d0
  	  do ir=1,nrg-1
	    do it=1,ntg-1
   	  do ip=1,npg-1
  		 densav=densvararr(ir,it,ip)+densav
  		 densav2=densvararr(ir,it,ip)**2+densav2
  		 if (densvararr(ir,it,ip).lt.0) print*,'densvararr lt 0',densvararr(ir,it,ip)
  		 if (densvararr(ir,it,ip).gt.maxdensvar) maxdensvar=densvararr(ir,it,ip)
  		 if (densvararr(ir,it,ip).lt.mindensvar) mindensvar=densvararr(ir,it,ip)
      enddo
      enddo
      enddo
      densav=densav/rntot
       densav2=densav2/rntot
     std=sqrt(densav2-densav**2)
	print*,'fractal densav,std, min,max',densav,std,mindensvar,maxdensvar
    call diskdat_write('fractal dens ave',densav,'should be 1')
    call diskdat_write('fractal std',std,'standard deviation')
  !  write(cmsgnm,'(a,i2)')  'Switching RNG stream back to ',&
!	current_gen
 !   call wrimsg('GRIDSET',cmsgnm)
 !   call set_current_generator(current_gen)
!    call set_random_seed(i1)
    write(IDVOUT,*) ' '
  endif

!!scale surface density inside and outside gap to input scale factor
!! this uses the z-integration for sigma; still deciding which to use
!  if (iskip.eq.1) then
!  if (igapd.eq.1) then
!	  ! set scale factor for disk mass.  involves looping through and
!	  ! calculating mass twice.  :)
!    do id=1,2
!      rhogap0(id)=rhoscale_gap(id)*rho0(id)*rgapd2**(b(id)-a(id))/rgapd2**(bgap(id)-agap(id))*z1(id)/z1gap(id)
!      print*,'rhoscale_gap,rho0',rhoscale_gap(id)*rho0(id)
!    end do
!    rhogap0(5)=fmass(5)/fmass(1)*rhogap0(1)
!    rhogap0(6)=fmass(6)/fmass(2)*rhogap0(2)
!    rhogap0(3:4)=0.d0
!    print*,'rho0',rho0
!    print*,'rhogap0',rhogap0
!    print*,'rhoscale_gap',rhoscale_gap
!    print*,'z1gap',z1gap
!  endif
!  endif

print*,''

!this integrates along theta for sigma calculation
!first assume rhogap0 is same as rho0
rhogap0=rho0
!print*,'rhogap0',rhogap0
!now calculate disk surface density along theta at gap
do id=1,2
  call locate(rarr,nrg,rgapd2,irgap)
!  print*,'rarr(irgap),rarr(irgap+1),rgapd2',rarr(irgap),rarr(irgap+1),rgapd2
 ! rarr(irgap+1)=rgapd2   !WHAT?   only if igapd=1
  do ir=irgap,irgap+1   !want just inside gap and just outside gap.
  print*,'rarr(ir),rgapd2',rarr(ir),rgapd2
  i=ir-irgap+1  !should range from 1-2
  print*,'i',i
  if (i.eq.1) then
	rad=rgapd2*0.9d0
  else
	rad=rgapd2*1.1d0
  endif
  sigtmp(i)=0.d0
  do it=1,ntg-1
	  thet=thetavearr(it)
    cost=cos(thet)
    sint=sin(thet)
    dcost=costarr(it)-costarr(it+1)
 !   print*,'dcost',dcost
	do ip=1,npg-1 !3-D atmosphere
        phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
   !     print*,'phi',phi
        call densdisk(rad,sint,cost,phi,densd,id)
        sigtmp(i)=sigtmp(i)+abs(densd*rad*dcost*rmincgs)/dble(npg-1)
  !      print*,'densd,sigtmp(i)',densd,sigtmp(i)
     enddo  !phi loop
  enddo  !theta loop
  print*,'i,sigtmp(i)',i,sigtmp(i)
  enddo  !rad loop
  print*,'sigtmp',sigtmp
  print*,'rhoscale_gap(id)',rhoscale_gap(id)
  rhogap0(id)=rho0(id)*rhoscale_gap(id)*sigtmp(2)/sigtmp(1)
enddo
print*,'rho0',rho0
print*,'rhogap0',rhogap0
print*,'rhoscale_gap',rhoscale_gap
print*,'z1',z1
print*,'z1gap',z1gap
print*,'fmassd',fmass
!rhogap0(2)=1.351142946984419E-005
rhogap0(5)=fmass(5)/fmass(1)*rhogap0(1)
rhogap0(6)=fmass(6)/fmass(2)*rhogap0(2)
rhogap0(3:4)=0.d0

  if (massd.gt.0.d0) then
     massdisk=0.d0
     tmpmass=0.d0
     do id=1,ndg
        if(is_disk(id)) then
           do ir=1,nrg-1
              rad=ravearr(ir)
              dr=rarr(ir+1)-rarr(ir)
              dr3=rarr(ir+1)**3.d0-rarr(ir)**3.d0
              do it=1,ntg-1
                 thet=thetavearr(it)
                 cost=cos(thet)
                 sint=sin(thet)
                 dcost=cos(thetarr(it+1))-cos(thetarr(it))
                 dt=thetarr(it+1)-thetarr(it)
                 do ip=1,npg-1 !3-D atmosphere
                    phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
                    dp=phiarr(ip+1)-phiarr(ip)
                    vol = - dr3 * dcost * dp / 3.d0 * rmincgs**3.d0
                    ! vol=rad**2*sint*dt*dp*dr*rmincgs**3
                    call densdisk(rad,sint,cost,phi,densd,id)
                    massdisk=massdisk+densd*vol
                    tmpmass(id)=tmpmass(id)+densd*vol
                 end do
              end do
           end do
        end if
     end do
     massdisk=massdisk/msun
     tmpmass=tmpmass/msun
     diskscale=massd/massdisk
     do id=1,ndg
      if(tmpmass(id).gt.0.) then
       diskscl(id)=massd*fmass(id)/tmpmass(id)
      endif
     end do 
  else
     diskscale=1.d0
     diskscl=1.d0
  end if

  print*,'massdisk before scaling (all grains)',massdisk
  print*,'mass of each disk before scaling',tmpmass
  print*,'diskscale',diskscale
  print*,'diskscl',diskscl
  ! calculate density in grid
  massenv=0.d0
  massdisk=0.d0
  maxdens=0.d0
  masssg=0.d0
  totvol=0.d0
  tmpmass=0.d0

  do ir=1,nrg-1

     rad=ravearr(ir)
     dr=rarr(ir+1)-rarr(ir)
     dr3=rarr(ir+1)**3.d0-rarr(ir)**3.d0
          
     do it=1,ntg-1

        thet=thetavearr(it)
        cost=cos(thet)
        sint=sin(thet)
        tant=tan(thet)
        dt=thetarr(it+1)-thetarr(it)
        dcost=cos(thetarr(it+1))-cos(thetarr(it))

        do ip=1,npg-1 ! 3-D atmosphere

           phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
           dp=phiarr(ip+1)-phiarr(ip)

           vol = - dr3 * dcost * dp / 3.d0 * rmincgs**3.d0
           totvol=totvol+vol
           vcell(ir,it,ip)=vol
           
           ! Disk grains
           do id=1,ndg
             
              ir1=ir
              it1=it
              ip1=ip
 
              if(is_disk(id)) then
				 if (imisalign.eq.1.and.rad.lt.radmisalign) then
                   x=rad*sint*cos(phi)
                   y=rad*sint*sin(phi)
                   z=rad*cost
                   cosi=cos(incmisalign*pi/180.d0)
                   sini=sin(incmisalign*pi/180.d0)
                   print*,'cosi',cosi
				   print*,'sini',sini
                   x1=x*cosi-z*sini
                   y1=y
                   z1hi=x*sini+z*cosi
                   rsq1=x1**2+y1**2+z1hi**2
                   cosb1=z1hi/sqrt(rsq1)
                   ph1=atan2(y1,x1)
                   if (ph1.lt.0) ph1=ph1+r2p  
				   print*,'ph1',ph1
                   call locate(costarr,ntg,cosb1,it1)
                   call locate(r2arr,nrg,rsq1,ir1)
                   call locate(phiarr,npg,ph1,ip1)
                 endif
                 call densdisk(rad,sint,cost,phi,densd,id)
                 densarr(ir1,it1,ip1,id)=densd*diskscl(id) ! densd already includes fmassd scale factor in rho0
                 massdisk=massdisk+densd*diskscl(id)*vcell(ir1,it1,ip1)
                 tmpmass(id)=tmpmass(id)+densd*diskscl(id)*vcell(ir1,it1,ip1)
      !           sigarr(ir1,id)=sigarr(ir1,id)+abs(densarr(ir1,it1,ip1,id)*rad*dcost*rmincgs)/dble(npg-1)  !rad*dcost=dz
                 sigarr(ir1,ip1,id)=sigarr(ir1,ip1,id)+abs(densarr(ir1,it1,ip1,id)*rad*dcost*rmincgs)/dble(npg-1)  !rad*dcost=dz
      !           sigarr(ir,id)=sigarr(ir,id)+abs(densarr(ir,it,ip,id)*rad*dt*rmincgs)/dble(npg-1)  !try rad*dt instead; how come this makes no difference?
              end if

              if(is_envelope(id) .or. is_cavity(id)) then
				         if (ienvtype.eq.1) then
                  call densenv(rad,thet,cost,sint,phi,dense,radamb,id,inhole)
				         else
				          call densenvpw(rad,thet,cost,sint,phi,dense,radamb,id,inhole,rmincgs)
				         endif
                 densarr(ir1,it1,ip1,id)=densarr(ir1,it1,ip1,id)+dense*fmass(id)
                 massenv=massenv+dense*fmass(id)*vcell(ir1,it1,ip1)
              end if
              
              if (ifractal.eq.1.and.fractmask(id).eq.1) then
                densarr(ir1,it1,ip1,id)=densarr(ir1,it1,ip1,id)*(densvararr(ir1,it1,ip1))
     !       densarr(ir,it,ip,id)=densarr(ir,it,ip,id)*(densvararr(ir,it,ip))
              endif

              massarr(ir1,it1,ip1,id)=(densarr(ir1,it1,ip1,id)*vcell(ir1,it1,ip1))**0.25d0 ! bad programming:  massarr is really mass**0.25
              if(is_sg(id)) masssg=masssg+densarr(ir1,it1,ip1,id)*vcell(ir1,it1,ip1)


           end do

           tmp = 0.d0
           do id=1,5
              tmp=tmp+densarr(ir1,it1,ip1,id)
           end do

           if (tmp.gt.maxdens) maxdens=tmp

        end do
     end do
  end do

  maxdens = maxval(sum(densarr,dim=4))
  
! calculate calculate scale factor for sigarr (theta integration); also do an analytic comparison first (along z)
 print*,''
 print*,'calculating surface density of disk'
  if (massd.gt.0.d0) then
 do id=1,ndg
   if(is_disk(id)) then    
      ir=1   
      rad=ravearr(ir+1)
  !    rad=0.5d0*(rarr(ir+1)+rarr(ir+2))  
      do while (ir.lt.nrg-1)
  !       ip=1       
  !       phi=0.d0
  !       dp=r2p
         ir=ir+1   !want to start at ir=2
         rad=ravearr(ir)
         dr=rarr(ir+1)-rarr(ir)
         do ip=1,npg-1
	      dp=phiarr(ip+1)-phiarr(ip)
	      phi=0.5*(phiarr(ip+1)+phiarr(ip))
          if (igapd.eq.1.and.rad.lt.rgapd2.and.rad.gt.rgapd1) then
            adens=agap(id)
            bdens=bgap(id)
            zscale=z1gap(id)
            rhodens=rhogap0(id)
         else
            adens=a(id)
            bdens=b(id)
            zscale=z1(id)
            rhodens=rho0(id)
         endif
         sigpw=dsqrt(r2p)*zscale*rmincgs*rhodens*rad**(bdens-adens)*(1.d0-dsqrt(1./rad))
         if (sigpw.lt.0.d0) then
            print*,'huh?',sigpw
            stop
         endif
         sigpwtot=sigpwtot+sigpw*rad*dr*rmincgs**2*dp
         sigtot=sigtot+sigarr(ir,ip,id)*rad*dr*rmincgs**2*dp
       !  print*,'rad,ir,rmaxd',rad,ir,rmaxd
       !  print*,'id,ir,sigpw,sigarr',id,ir,sigpw,sigarr(ir,id)
       !  print*,sigpwtot,sigtot
      enddo
     enddo
   endif
 enddo
 print*,'sigpwtot,sigtot,input mass',sigpwtot/msun,sigtot/msun,massd

 sigscale=massd/sigtot*msun

! scale surface density; print it out also
!scale surface density and write out to file
  do id=1,ndg
	if (is_disk(id)) then
	open(unit=15,file='sig_theta'//char(48+id)//'.dat',status='unknown')
	write(15,*) nrg,' = number of grid points in r'
	write(15,*) &
	     & 'r/rstar    r(au)    sigma'
  do ir=1,nrg
  do ip=1,npg
    sigarr(ir,ip,id)=sigarr(ir,ip,id)*sigscale  
    write(15,911)  ravearr(ir), &
         & ravearr(ir)/autors, &
         & sigarr(ir,ip,id)
  enddo
  enddo
  close(15)
  endif
  enddo
911 format(3(1x,f15.5))


 endif

  if (iskip.eq.1) then !skip this, does z integration and we are using theta.
  if (massd.gt.0.d0) then
  sigpwtot=0.d0
  sigtot=0.d0
!  print*,'rmaxd',rmaxd
  itmid=ntg/2+1
  do id=1,ndg
    if(is_disk(id)) then    
	  do ircyl=2,nrg-1
	     rcyl=ravearr(ircyl)
		do ip=1,npg-1
		   sigarr(ircyl,ip,id)=0.d0
		   z=0.d0
		   it=itmid
		   delz=costarr(itmid-1)*rcyl*1.01
		!   print*,'delz',delz

		
		 do while (z.lt.rmax.and.it.gt.2)
			z=z+delz
		!	cosb=z/rcyl
			tanb=(rcyl/z)
			cosb=cos(atan(tanb))
		    call locate(costarr,ntg,cosb,it)
		    rad=z/cosb
		    call locate(rarr,nrg,rad,ir)
		    sigarr(ircyl,ip,id)=sigarr(ircyl,ip,id)+2.d0*densarr(ir,it,ip,id)*delz*rmincgs   !rmincgs converts delz to cm
		   ! print*,costarr(it-1)*rcyl-costarr(it)*rcyl,delz
		    delz=max(costarr(it-1)*rcyl-costarr(it)*rcyl,delz)
		   ! print*,'it,delz,z',it,delz,z
		   ! print*,'z,delz,densarr',z,delz, densarr(ir,it,ip,id)
		   ! print*,'sigarr',sigarr(ir,id)
		 enddo
		
	    enddo
	    sigarr(ircyl,ip,id)=sigarr(ircyl,ip,id)/dble(npg-1)
	
	!   for fun, do it analytically 	
		 dr=rarr(ircyl+1)-rarr(ircyl)
         if (igapd.eq.1.and.rcyl.lt.rgapd2.and.rcyl.gt.rgapd1) then
           adens=agap(id)
           bdens=bgap(id)
           zscale=z1gap(id)
           rhodens=rhogap0(id)
        else
           adens=a(id)
           bdens=b(id)
           zscale=z1(id)
           rhodens=rho0(id)
        endif
        sigpw=dsqrt(r2p)*zscale*rmincgs*rhodens*rcyl**(bdens-adens)*(1.d0-dsqrt(1./rcyl))
        if (sigpw.lt.0.d0) then
           print*,'huh?',sigpw
           stop
        endif
          sigpwtot=sigpwtot+sigpw*rcyl*dr*rmincgs**2*r2p
          sigtot=sigtot+sigarr(ircyl,ip,id)*rcyl*dr*rmincgs**2*r2p
   !       print*,'rcyl,ir,rmaxd',rcyl,ir,rmaxd
   !       print*,'id,ir,sigpw,sigarr',id,ir,sigpw,sigarr(ir,id)
   !       print*,sigpwtot,sigtot 		
	  enddo
      endif
    enddo
  print*,'sigpwtot,sigtot,input mass',sigpwtot/msun,sigtot/msun,massd

  sigscale=massd/sigtot*msun

! scale surface density; print it out also
!scale surface density and write out to file
  do id=1,ndg
	if (is_disk(id)) then
	open(unit=15,file='sig_z'//char(48+id)//'.dat',status='unknown')
	write(15,*) nrg,' = number of grid points in r'
	write(15,*) &
	     & 'r/rstar    r(au)    sigma'
  do ir=1,nrg
  do ip=1,npg
    sigarr(ir,ip,id)=sigarr(ir,ip,id)*sigscale  
    write(15,912)  ravearr(ir), &
         & ravearr(ir)/autors, &
         & sigarr(ir,ip,id)
  enddo
  enddo
  close(15)
  endif
  enddo
912 format(3(1x,f15.5))

  endif
  endif

!!testing
!this doesn't work. 
!  call output_grid('diffus_old',diffus)
! ! testing.  redo diffusion region based on disk density
! if(diffusion) then
!   ! set up diffusion grid
!    dcount=0
!    do ir=nr0,nrg-1
!    do it=1,ntg-1
!    do ip=1,npg-1 ! 3-D atmosphere
!      !only do big grain region
!	  !  rhoave=0.d0
!	!	do id=1,ndg
!	!	   if(is_disk(id)) then
!	!	      rhoave=rhoave+densarr(ir,it,ip,id)
!	!	   end if
!	!	end do
!	  rhoave=densarr(ir,it,ip,1) 
!	   if (rhoave.gt.1.e-8) then
!		  diffus(ir,it,ip) = .true.
!		else
!		  diffus(ir,it,ip) = .false.
!		endif
!!		if (igapd.eq.1.and.r.lt.rgapd2.and.r.gt.rgapd1) then
!1	         diffus(ir,it,ip)=.false.
!!	         diffdir(ir,it,ip)=0
!!	     endif
!	enddo
!	enddo
!	enddo
!  endif
!  call output_grid('diffus_new',diffus)

!this doesn't work
! test diffusion region
! if (0) then
! if(diffusion) then
!   ! set up diffusion grid
!    dcount=0
!    do ir=nr0,nrg-1
!    do it=1,ntg-1
!    do ip=1,npg-1 ! 3-D atmosphere
!		do id=1,ndg
!		   if(is_disk(id)) then
!			  if (diffus(ir,it,ip)) then
!			  if (densarr(ir,it,ip,id).lt.1.d-9) then
!				print*,'density low, setting diffusion to false'
!				diffus(ir,it,ip)=.false.
!		      end if
!	          endif
!		   endif
!		end do
!	enddo
!	enddo
!	enddo
!  endif
!  endif
 


 ! print*,'z1,b,rho0',z1,b,rho0
  massenv=massenv/msun
  massdisk=massdisk/msun
  masssg=masssg/msun
  tmpmass=tmpmass/msun
  ! important note:  densarr(ir,it,ip) is the density halfway between
  ! ir and ir+1, it and it+1, ip and ip+1; it is NOT the density at
  ! ir,it,ip.
  print*, 'massenv (all grains) ',massenv
  print*, 'massdisk from grid, compared to input ',massdisk, &
       & massd
  print*,'disk mass of each dust type',tmpmass
  ! masssg is actually mass sg+20-200 A grains
  print*, 'mass very small grains (non-thermal) (solar) ', masssg
  print*,'massvsgs/masstot',masssg/(massdisk+massenv)
  print*,'max density in disk/envelope',maxdens
  print*,'min radius where rho=rhoamb ',radamb/autors
  print*, ' '
  print*,'totvol',totvol

  !calculate envelope density at 1 AU and thet=60 degrees.  average over phi if clumpy envelope
  call locate(thetarr,ntg,90.d0*deg2rad,it)
  call locate(rarr,nrg,autors,ir)
 ! print*,'it,ir',it,ir    !correct...
  dens1au=0.d0
  do id=1,ndg
  do ip=1,npg-1
	if(is_envelope(id)) then
		dens1au=dens1au+densarr(ir,it,ip,id)
	!	print*,'hi, ip,id',ip,id
	!	print*,'densarr(ir,it,ip,id)',densarr(ir,it,ip,id)
	endif
  enddo
  enddo
  dens1au=dens1au/float(npg-1)
  print*,'envelope density at 1 AU in midplane',dens1au,' gm/cm^3'
  call diskdat_write('dens1au',dens1au, &
       & 'envelope density at 1 AU in midplane')

  open(unit=15,file='Av_view.dat',status='unknown')
  ! calculate average optical depth (A_v) along all theta & phi directions
  write(15,*) &
       & 'Av along viewing directions'
  write(15,*) &
       & 'theta       phi      Av_env    Av_disk      Av_sg     Av_tot'
  tauave=0.d0
  dcount=0
  dphi=r2p/dble(nph)
  do i=1,nmu
     thet=acos(u(i))
     if (thet.eq.0.d0) thet=1.d-5
     sint=sin(thet)
     call locate(thetarr,ntg,thet,it)
     if (ntg.eq.1) it=1
     ! print*,'it',it
     do j=1,nph
        phi=dble(j)*dphi-0.5d0*dphi
        ! print*,'phi',phi
        if (phi.gt.r2p) print*,'oops!  phi',phi
        call locate(phiarr,npg,phi,ip)
        tau=0.d0
        taud=0.d0
        taue=0.d0
        taup=0.d0
        cost=cos(thet)
        do ir=1,nrg-1
           dr=rarr(ir+1)-rarr(ir)
           rad=ravearr(ir)
           ! print*,'ir,it,ip',ir,it,ip
           do id=1,2
              taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr
           end do
           do id=3,4
              taue=taue+kapd(id)*densarr(ir,it,ip,id)*dr
           end do
           do id=5,8
              taup=taup+kapd(id)*densarr(ir,it,ip,id)*dr
           end do
        end do
        write(15,901) thet*rad2deg,phi*rad2deg &
             & ,taue*1.086d0,taud*1.086d0,taup*1.086d0 &
             & ,(taue+taud+taup)*1.086d0
901     format((f10.5,1x,f10.5),5(1x,1pe12.5))
        tauave=tauave+taud+taue+taup
        dcount=dcount+1
     end do
  end do
  tauave=tauave/dble(dcount)
  print*,'dcount',dcount
  print*,'A_V average over all directions ',tauave*1.086d0
  write(15,*) 'A_V average over all directions ',tauave*1.086d0
  call diskdat_write('Av_ave',tauave*1.086d0, &
       & 'A_V ave over all directions')
  close(15)

  ! stop

  ! integrate optical depth along all the grid angle directions
  tauave=0.d0
  dcount=0
  open(unit=15,file='Av_grid.dat',status='unknown')
  write(15,*) 'Av through each grid theta bin '
  write(15,*) &
       & 'theta    phi    AV_env       AV_disk    AV_sg      AV_tot'
  do ip=1,npg-1
     phi=phiavearr(ip)
     do it=1,ntg-1
        tau=0.d0
        taud=0.d0
        taue=0.d0
        taup=0.d0
        thet=thetavearr(it)
        if (thet.gt.pi/2.d0-eps.and.thet.lt.pi/2.d0+eps) then
           ! TSC blows up at thet=90
           print*,'thet',thet*rad2deg
           thet=89.999d0*deg2rad
        end if
        cost=cos(thet)
        sint=sin(thet)
        phi=0.d0
        do ir=1,nrg-1
           dr=rarr(ir+1)-rarr(ir)
           rad=ravearr(ir)
           do id=1,2
              taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
           end do
           do id=3,4
              taue=taue+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
           end do
           do id=5,8
              taup=taup+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
           end do
        end do
        write(15,901) thet*rad2deg,phi*rad2deg &
             & ,taue &
             & ,taud,taup,taue+taud+taup
        tauave=tauave+taud+taue+taup
        dcount=dcount+1
     end do
  end do
  tauave=tauave/dble(dcount)
  print*,'dcount',dcount
  print*,'A_V average over all directions ',tauave
  call diskdat_write('Av_ave',tauave, &
       & 'A_V ave over all directions')
  write(15,*) 'A_V average over all directions ',tauave
  close(15)
  ! 902     format(f10.5,4(1x,f15.5))

  ! calculate Av along ethet,ephi directions (peeled image)
  open(unit=15,file='Av_peel.dat',status='unknown')
  write(15,*) 'Av along peeled direction ethet'
  write(15,*) &
       & 'theta     phi      Av_env     Av_disk    AV_sg     Av_tot'
  do peelid=1,npeel
     thet=thete_arr(peelid)
     cost=coste_arr(peelid)
     sint=sinte_arr(peelid)
     call locate(thetarr,ntg,thet,it)
     phi=phie
     call locate(phiarr,npg,phi,ip)
     taud=0.d0
     taue=0.d0
     taup=0.d0
     do ir=1,nrg-1
        dr=rarr(ir+1)-rarr(ir)
        rad=ravearr(ir)
        do id=1,2
           taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
        do id=3,4
           taue=taue+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
        do id=5,8
           taup=taup+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
     end do
     write(15,901) thet*rad2deg,phi*rad2deg &
          & ,taue &
          & ,taud,taup,taue+taud+taup
     write(15,*) 'average A_v along *all* directions ', &
          & tauave
  end do
  close(15)
  
  ! calculate Av along the disk midplane
  open(unit=15,file='Av_90.dat',status='unknown')
  write(15,*) 'Av in the disk midplane'
  write(15,*) &
       & 'theta     phi      Av_env     Av_disk    AV_sg     Av_tot'
     thet=pi/2.d0
     cost=0.d0
     sint=1.d0
     call locate(thetarr,ntg,thet,it)
     phi=0.d0
     call locate(phiarr,npg,phi,ip)
     taud=0.d0
     taue=0.d0
     taup=0.d0
     do ir=1,nrg-1
        dr=rarr(ir+1)-rarr(ir)
        rad=ravearr(ir)
        do id=1,2
           taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
        do id=3,4
           taue=taue+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
        do id=5,8
           taup=taup+kapd(id)*densarr(ir,it,ip,id)*dr*1.086d0
        end do
     end do
     write(15,901) thet*rad2deg,phi*rad2deg &
          & ,taue &
          & ,taud,taup,taue+taud+taup
     write(15,*) 'average A_v along *all* directions ', &
          & tauave
  close(15)


  ! calculate tau_V=1 surf from origin (not quite the stellar surface!).
  open(unit=15,file='taudsurfin.dat',status='unknown')
  write(15,*) 'cyl r, z values for tau=1 surface'
  do it=1,ntg-1
     ! thet=(90.-(45./2000.*(it-1)))*deg2rad
     ! call locate(thetarr,ntg,thet,it)
     taud=0.d0
     cost=cos(thet)
     sint=sin(thet)
     phi=0.d0
     ip=1
     ir=0
     do while (taud.lt.1.d0.and.ir.lt.(nrg-1))
        ir=ir+1
        dr=rarr(ir+1)-rarr(ir)
        rad=ravearr(ir)
        do id=1,2
           taud=taud+kapd(id)*densarr(ir,it,ip,id)*dr
        end do
     end do
     x=rad*sint
     z=rad*cost
     write(15,*) x/autors,z/autors
  end do
  close(15)

  ! r=1.d0*autors
  ! call locate(rarr,nrg,dble(r),ir)
  ! dens1=densarr(ir,1,1)
  ! print*,'r,dens1',r,dens1

  ! write out massarr
  call output_grid('marr',massarr**4)

  ! write out densarr
  call output_grid('darr',densarr)
  call output_grid('darrtot',sum(densarr,dim=4))
  call output_grid('tarr',tdust2)
  call output_grid('darr_000',densarr)

  call output_grid('diffus',diffus)
  call output_grid('diffdir',diffdir)


  ! oops, not f77 compatible
  ! masstot = sum(densarr) / msun
  ! massdisk = sum(densarr(:,:,:,1)*vcell + densarr(:,:,:,2)*vcell) / msun
  ! massenv = sum(densarr(:,:,:,3)*vcell) / msun
  ! masssgreg = sum(densarr(:,:,:,4)*vcell) / msun

  masstot=massenv+massdisk+masssg

  call diskdat_write('masstot',masstot,'total circumstellar mass')
  call diskdat_write('massenv',massenv,'envelope mass without sgs')
  call diskdat_write('massdisk',massdisk,'disk mass without sgs')
  call diskdat_write('masssgs',masssg,'total  mass, small grains')

  print*,'done with gridset'

end subroutine gridset










