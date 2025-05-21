subroutine disk(converge,iter,i1)

! monte carlo radiative transfer for disk surrounding star.
! uses cartesian coordinates, arbitrary density distribution,
! arbitrary disk structure.
! This program works for extremely (geometrically) thin disks of large
! radial extent plus tenuous envelope.  In order to do the
! radiative transfer in the thin disk, the opacity is calculated
! at each step, and the photon path integration is variable.
! In the disk, the step size is small for directions perpendicular
! to the disk, say zmindisk/10.  Outside the
! disk, z gt. zmaxdisk, the step size is much larger, zmax/100
! or zmax/200, something on that order.
! Note that zmindisk is probably about .1 stellar radius.
! zmax is about 10000 stellar radii, that is, about 100 AU.
! The envelope extend is even larger for protostars--10^4 AU,
! so the step size increases with distance from the source.

! calls these subroutines:
! stokes
! initp
! opacdisk
! opacinfall
! dels
!
  ! history:
  ! 99/04/12 (mjw): add cputime for photon loop update output
  ! (this breaks the overall CPU clock in newtts)
  ! 99/11/20 (baw):  big changes..
! 2012/03/08-10 (mjw):  MPICH stuff; convert output statements to unit IDVOUT
!                    add use of wrimsg_mod
!                    random cleanup of diagnostic output
!
!
! *************************************************************

  use tts_mod
  use grid_mod
  use stokesvar_mod
  use filt_mod
  use taunum_mod
  use log_mod
  use opacin_mod
  use out_mod
  use peeloff_mod
  use dust_mod
  use random
  use spot_mod
  use constants
  use output_mod
  use wrimsg_mod
#ifdef MPI
  use ttsre_mpi_mod
#endif

  implicit none

  real*8 :: sini(nmu),xmax
  real*8 :: thetb
  real*8 :: tsum,coslp,sinlp
  real*8 :: cosnorm,sipnew
  real*8 :: rmincgs,tcount,tau_therm,rn,lstar

  integer :: ii(3),ns,icount,ithet,iphot,ir
  integer :: it,ip,npacc,ips,idust2,id,npaccstar
  integer :: ispotph,iter,nsum,i1,i_orig_init,inub

  character(len=4) :: suffix

  integer :: peelid

  ! cpu time variables
  real*8 :: cpusec,rjunk
  character(len=70) :: cmsgnm
  character(len=3) :: converge

  write(6,*) 'check, pi',pi
  icount=0

  ! set up density grid
  if (iter.lt.2) then
     write(IDVOUT,*) 'iteration #',iter,'calling gridset'
     if (ihydro.eq.1) then
     call gridsethydro(iter,i1)
     else
	 call gridset(iter,i1)
	 endif
  end if

  print*,''
  if (converge.eq.'yes') write(IDVOUT,*) 'final iteration'
  print*,'iter',iter

  ! set up filter functions if making peeled images
  if (ipeel.eq.1.and.output_imfilt) call filt

  ! The flux is summed to image(ix,iy,it)
  ! polarization to imagei(ixp,iyp,it),imageq...,imageu...

  ! inclination arrays.
  do ithet=1,nmu
     sini(ithet)=sqrt(1.d0-u(ithet)**2)
  end do

  ! if the inner hole is large, star is effectively point source
  ! and theta,phi grid cells are really close together at rmin,
  ! so emit photons into radial direction so grid can handle it.
  if (rddustmin.gt.100.d0) then
     print*,'emitting from star as point source'
     ips=1
  else
     ips=0
  end if

  ! call this after gridset since some variables are set there.
  if (ispot.eq.1.and.spotflag.eq.1) then
     print*,'calling spotset'
     call spotset()
  end if
  ! only call spotset once, it will set spotflag=0 after first call

  nscat=0.d0
  tot=0.d0

  ! origin of photons, set to 0 for now...
  i_orig=0

  ! x,y arrays for peeled-off image
  fractxh=2.d0*rmaxi/dble(nxhst)
  xmax=rmaxi

  rmincgs=rstar*rsol

  ! np=npsav
  print*,'np',np
  flux=0.d0
  sflux=0.d0
  aflux=0.d0
  abflux=0.d0
  scount=0.d0
  dscount=0.d0

  npout=int(isrf_frac*dble(np))
  rn=isrf_frac*dble(np)
  print*,'rn',rn

  if ((rn-npout).gt.0.5d0) npout=npout+1
! number of stellar photons
  ns=int(dble(np)*(1.d0-accfrac-isrf_frac))
  rn=real(np)*(1.d0-accfrac-isrf_frac)
  print*,'rn',rn

  if ((rn-ns).gt.0.5d0) ns=ns+1
  npacc=int(dble(np)*accfrac)
  rn=int(dble(np)*accfrac)
  print*,'rn',rn

  if ((rn-npacc).gt.0.5d0) npacc=npacc+1
!  number of hotspot photons
  npaccstar=int(accfrac2*dble(np))
  rn=int(dble(ns)*accfrac2)
  print*,'rn',rn

  if ((rn-npaccstar).gt.0.5d0) npaccstar=npaccstar+1
  nsum=ns+npacc+npout

  print*,'ns (star),npaccstar, npacc (disk),npout (outside),sum,should=np'
  print*, ns,npaccstar,npacc,npout,nsum,np

  if (nsum.ne.np) then
     print*,'resetting ns'
     ns=ns+np-nsum
     print*,'ns (star),npacc (disk),npout (outside),sum,should=np'
     nsum=ns+npacc+npout
     print*, ns,npacc,npout,nsum,np
  end if
  print*,'np_disk,np_shock',npacc,npaccstar
  print*,'np_spot = np_shock = npaccstar'
  
  lstar=4.d0*pi*sigt*tstar**4*rmincgs**2
  print*,'ltot/np*(ns-npaccstar),lstar',ltot/lsun/dble(np)*dble(ns-npaccstar),lstar/lsun

  ! convert Tstar to Jon's units
  if (iter.eq.1) then
     tstar=tstar/11605.d0
     tstarave=tstarave/11605.d0
     tshock=tshock/11605.d0
  end if

  ! ********  do loop over stellar photons  *******
  print*,'ns',ns
  print*,'npout',npout

  do iphot=1,np

     on_wall = .false.

     if(mod(iphot,iwrite).eq.0) then
        ! cpu time
        call cpu_time(cpusec)
        ! use explicit format to keep g77 happy
        write(cmsgnm,'(i12,a20,f11.2,a)') iphot,' photons completed. ',cpusec,' = CPU time (sec)'
        call WRIMSG('MAIN',cmsgnm)
     end if



     ! if (ran().lt.accfrac) then
     if (iphot.lt.(npacc)) then ! DISK ACCRETION PHOTONS

        if (iphot.eq.1) print*,'doing disk accretion photons now'

        call initpacc(ii,np,ns,iter)
        call opacset(nub)
        ! bin photon into initial emitted spectrum array
        inub=min(max(int(dble(nfreq)*log(nub/numin)/lnurat)+1,1),nfreq)
        i_orig_init=4
        si_init(inub,i_orig_init)=si_init(inub,i_orig_init)+1.
        si_init(inub,1)=si_init(inub,1)+1.

        wave=1.2398d0/nub

        ir=ii(1)
        it=ii(2)
        ip=ii(3)

        do id=1,ndg
           kapd(id)=kappa(id)*rmincgs*kappav(id)
           ! rlam(id)=albedo(id)
        end do

        iflag=0

        ! origin of photon, 2 for disk
        i_orig=2
        ! else
        ! i_orig=3
        ! end if
        ! i_orig=2

        if (ipeel.eq.1) then
           do peelid=1,npeel
              ! 20090730, idust is not defined, should be idust2, but should be okay, recalculated in peeloff
              sipnew = sip
              call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp &
                   & ,cost,sint,cosp,sinp,phi &
                   & ,hit,htot,rsq,rtot &
                   & ,tsum,ii(1),ii(2),ii(3),idust,idust2,iflag,iphot &
                   & ,peelid)
           end do
        end if

        ! ******  check peeling off of accretion photon!!!!!!!  071902
        ! set normalization!!!!!

     else if (iphot.le.(npacc+ns)) then ! STELLAR AND HOT SPOT PHOTON

        if (iphot.lt.(npacc+npaccstar)) then ! HOT SPOT PHOTON
           if (iphot.eq.npacc+1) print*,'doing hot spot photons now'

           ! 20080828 BAW selecting region inside hotspot
           iplanckst=1
           ispotph=1
           if (ips.eq.0) then
              call initpspot(ispotph)
           else
              call initp_ps(ispotph)
           end if

           ! emit half the photons as x-rays
           if (ran().lt.0.5d0) then
              nub=1.2398d0/0.06d0+ran()*(1.2398d0/0.015d0-1.2398d0/0.06d0)
           end if

          ! bin photon into initial emitted spectrum array
          inub=min(max(int(dble(nfreq)*log(nub/numin)/lnurat)+1,1),nfreq)
          i_orig_init=3
          si_init(inub,i_orig_init)=si_init(inub,i_orig_init)+1.
          si_init(inub,1)=si_init(inub,1)+1.

        else ! STELLAR PHOTONS

           if (iphot.eq.npacc+npaccstar+1) print*,'now doing stellar photons'

           ! 20080828 BAW, selecting region outside hotspot
           iplanckst=0
           ispotph=0
           if (ips.eq.0) then
              call initpspot(ispotph)
           else
              call initp_ps(ispotph)
           end if

          ! bin photon into initial emitted spectrum array
          inub=min(max(int(dble(nfreq)*log(nub/numin)/lnurat)+1,1),nfreq)
          i_orig_init=2
          si_init(inub,i_orig_init)=si_init(inub,i_orig_init)+1.
          si_init(inub,1)=si_init(inub,1)+1.

        end if

        ! origin of photon i_orig = 1 for star
        i_orig=1

        call opacset(nub)
        wave=1.2398d0/nub


        do id=1,ndg
           kapd(id)=kappa(id)*rmincgs*kappav(id)
        end do

        coslp=cos(lp)
        sinlp=sin(lp)

        zp=cosb*rmin
        rp=sinb*rmin
        xp=rp*coslp
        yp=rp*sinlp

        rsq=zp**2+rp**2
        rtot=sqrt(rsq)

        ux=sint*cosp
        uy=sint*sinp
        uz=cost

        iflag=0
        ii(1)=1 ! index of radial grid (stellar surface is 1)
        on_wall = .true.

        thetb=acos(cosb)
        call locate(thetarr,ntg,thetb,it)
        ii(2)=it

        call locate(phiarr,npg,lp,ip)
        ii(3)=ip

        ! first, peel off direct flux
        ! weight photon intensity by angle between normal and
        ! photon direction, since emitted from a surface.
        if (ipeel.eq.1.and..not.partial_peeloff) then
           do peelid=1,npeel

              if (ips.eq.0) then

                 cosnorm=cosb*coste_arr(peelid)+(sinb*sinte_arr(peelid)*(cospe_arr(peelid)*coslp+sinpe_arr(peelid)*sinlp))

                 if (limb.eq.0) then
                    ! intensity constant, energy/Sr proportional to mu
                    sipnew=4.d0*sip*cosnorm !normalization = 2
                 else
                    ! intensity goes as (1+mu), energy/Sr has another factor of mu.
                    sipnew=12.d0/5.d0*sip*(cosnorm+cosnorm*cosnorm)
                 end if

                 if (cosnorm.gt.0.d0) then
                    call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp &
                         & ,cost,sint,cosp,sinp,phi &
                         & ,hit,htot,rsq,rtot &
                         & ,tsum,ii(1),ii(2),ii(3),idust,idust2,iflag &
                         & ,iphot,peelid)
                 end if

              else

                 sipnew=sip*2.d0
                 call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp &
                      & ,cost,sint,cosp,sinp,phi &
                      & ,hit,htot,rsq,rtot &
                      & ,tsum,ii(1),ii(2),ii(3),idust,idust2,iflag,iphot &
                      & ,peelid)

              end if

           end do
        end if

     else ! EXTERNAL ILLUMINATION

        if (iphot.eq.npacc+ns+1) print*,'now doing external photons'

        call initpout()
        on_wall=.true.
        i_orig=4
        call opacset(nub)
        wave=1.2398d0/nub
        
		! bin photon into initial emitted spectrum array
        inub=min(max(int(dble(nfreq)*log(nub/numin)/lnurat)+1,1),nfreq)
        i_orig_init=5
        si_init(inub,i_orig_init)=si_init(inub,i_orig_init)+1.
        si_init(inub,1)=si_init(inub,1)+1.
 
        do id=1,ndg
           kapd(id)=kappa(id)*rmincgs*kappav(id)
        end do

        coslp=cos(lp)
        sinlp=sin(lp)
        zp=cosb*rmax*0.999999
        rp=sinb*rmax*0.999999
        xp=rp*coslp
        yp=rp*sinlp

        rsq=zp**2+rp**2
        rtot=sqrt(rsq)

        ux=sint*cosp
        uy=sint*sinp
        uz=cost

        iflag=0
        ii(1)=nrg-1          !index of radial grid at rmax

        thetb=acos(cosb)
        call locate(thetarr,ntg,thetb,it)
        ii(2)=it

        call locate(phiarr,npg,lp,ip)
        ii(3)=ip

        ! first, peel off direct flux
        ! I think we weight photon simply by 1/4pi, or just one?  no, 2 because
        ! we only count half the flux.
        ! however, we don't want to include these photons unless they
        ! interact. anlogous to observing, it's like subtracting off the sky
        ! background. so don't call peel here.

! uncomment these if you want to include background radiation that doesn't interact
!        if (ipeel.eq.1) then
!          do peelid=1,npeel
!            cosnorm=cosb*coste_arr(peelid)+(sinb*sinte_arr(peelid)*(cospe_arr(peelid)*coslp+sinpe_arr(peelid)*sinlp))
!            sipnew=sip*2.d0
!            if (cosnorm.lt.0.d0) then
!              call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp &
!                     & ,cost,sint,cosp,sinp,phi &
!                     & ,hit,htot,rsq,rtot &
!                     & ,tsum,ii(1),ii(2),ii(3),idust,idust2,iflag,iphot &
!                     & ,peelid)
!            end if
!          enddo
!        end if

     end if

     call propagate(iphot,sini,ii,np,ns,idust2,iter)

 end do

  ! ***** end of loop over stellar photons  ******

  ! np=ns

  write(IDVOUT,*) 'photons killed: ', sum(killed)
  write(IDVOUT,*) 'SED photons used: ', sum(nums(:,:,:,:,1))
  write(IDVOUT,*) 'Number of photons: ',np 
#ifdef MPI
  write(IDVOUT,*) 'at Reduction Barrier...'
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)



! some scalars
  call MPI_ALLREDUCE(MPI_IN_PLACE,icount,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,np,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dscount,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,scount,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  write(IDVOUT,*) 'REDUCE: Number of photons: ',np 

! quantities which may need to be summed, but are not explicitly in the routine
! KILLED (reset before each call to DISK)
  icountbuf = nrg*ntg*npg
  call MPI_ALLREDUCE(MPI_IN_PLACE,killed,icountbuf,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  write(IDVOUT,*) 'REDUCE: photons killed: ', sum(killed)

! quantities which may need to be summed, but are not explicitly in the routine
! NABS,DBIG,DSG (reset before each call to DISK)  
! JMEAN not needed as it is calculated from DBIG and DSG
  icountbuf = nrg*ntg*npg*ndg
  call MPI_ALLREDUCE(MPI_IN_PLACE,nabs,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dtauabs,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  icountbuf = nrg*ntg*npg
  call MPI_ALLREDUCE(MPI_IN_PLACE,dbig,icountbuf,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dsg,icountbuf,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,jmean,icountbuf,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

! standard SED quantities
! SI,SQ,SU,SV,SI2,SQ2,SU2,SV2,NUMS (these are reset before each call to DISK)
  icountbuf = nfreq*nmu*nph*nap*no
  call MPI_ALLREDUCE(MPI_IN_PLACE,si,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sq,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,su,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sv,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,si2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sq2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,su2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,sv2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(MPI_IN_PLACE,nums,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  write(IDVOUT,*) 'REDUCE: SED photons used: ', sum(nums(:,:,:,:,1))

! do reduction for arrays available ONLY when temperature is converged
  if (converge == "yes" ) then 
     write(IDVOUT,*) 'Reducing image arrays...'

! peeloff SED quantities
! TI,TQ,TU,TV,TI2,TQ2,TU2,TV2,NUMT (these should be reset before each call to DISK, and
!  at present, this is the case in TTS)
!
     icountbuf = nfreq*npeel*nap*no
     call MPI_ALLREDUCE(MPI_IN_PLACE,ti,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tq,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tu,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tv,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,ti2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tq2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tu2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,tv2,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,numt,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     write(IDVOUT,*) 'REDUCE: peeled SED photons used: ', sum(numt(:,:,:,1))

! broadband peeloff IMAGES
   if(output_imfilt) then
     icountbuf = npeel*nxhst*nxhst*nbnd
     call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_i,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_q,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_u,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,image_b_v,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
  endif

 if(output_imcube) then
! monochromatic peeloff IMAGES
     icountbuf = npeel*nxhst*nxhst*nfreq*no
     call MPI_ALLREDUCE(MPI_IN_PLACE,image_m_i,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     icountbuf = npeel*nxhst*nxhst*nfreq
     call MPI_ALLREDUCE(MPI_IN_PLACE,image_m_q,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,image_m_u,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
     call MPI_ALLREDUCE(MPI_IN_PLACE,image_m_v,icountbuf,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
 endif

  endif
! the ALLREDUCE calls are synchronous, but let's just be safe
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  write(cmsgnm,'(a,i2)') 'Done with reduction in process ',myid
  call wrimsg('DISK',cmsgnm)

#endif



  call wrimsg('DISK','done with photons')

  ! calculate final temp using lucy method
  call wrimsg('DISK','calling tfinal')

  print*,'np before tfinal',np
  call tfinal(np,iter,converge)

  if (ihseq.eq.1.and.iter.ge.iter_hseq) then
	  call hseq(iter)
  endif

  ! not sure I want to do this yet, at least while testing
  ! if (iter.lt.3.and.converge.eq.'yes') then
  ! print*,'resetting converge=no since iter < 3'
  ! print*,'i.e., requiring at least 3 iterations'
  ! converge='no'
  ! end if

  call int_mean(ns+npout+npacc)

  write(IDVOUT,*) 'fraction of phots scattered outside image'
  write(IDVOUT,*) real(icount)/real(np)
  write(IDVOUT,*) np,scount,dscount
  tcount=(dble(np)-scount-dscount)/dble(np)
  tau_therm=-log(1.d0-tcount)

  call diskdat_section('newdisktherm')
  call diskdat_write('tcount',tcount,'')
  call diskdat_write('tautherm',tau_therm,'')

  if (sum(fmass,mask=is_sg).gt.1.d-6) then
     rjunk=real(imave)/real(imct)
  else
     rjunk=0.
  end if

  write(IDVOUT,*) 'imave,imct',imave,imct
  write(IDVOUT,*) 'average imean index for PAH emissivity',rjunk
  call diskdat_write('imave/imct',rjunk,'average imean index for PAH emissivity')

! 20120220  taking this out.  need a diffdir to go along, plus doesn't make sense in 2-D diffusion.  can put in with 3-D.
 ! if (diffusion) then
!  do ir=3,nrg-1
 !    do it=1,ntg-1
!        do ip=1,npg-1
!           if(killed(ir,it,ip) > 1 .and. sum(densarr(ir,it,ip,:), mask=is_thermal) > 0.) then
!              diffus(ir,it,ip) = .true.
!           end if
!        end do
!     end do
!  end do
!  endif

  write(suffix,'("_",I3.3)') iter

!  call output_grid('denstotal'//suffix,sum(densarr(:,:,:,1:4), dim=4))
  call output_grid('diffus' //suffix,diffus)
  call output_grid('diffdir' //suffix,diffdir)
  call output_grid('killed' //suffix,killed)
  
  write(IDVOUT,*) 'iter,massdisk,massenv',iter,massdisk,massenv

  return
end subroutine disk

! *********************************************************
