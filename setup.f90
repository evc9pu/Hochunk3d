! *******************************************************

subroutine setup(atname,dustname,i1)

  use tts_mod
  use grid_mod
  use stokesvar_mod
  use log_mod
  use opacin_mod
  use out_mod
  use random
  use tabl_mod
  use dust_mod, only : Tsub
  use constants
  implicit none

  real*8 :: rmine1,muc,rcgs,mcgs,honr,dela,rminsub,zmin1,zmin2

  character(len=80) :: atname,dustname(ndg)

  integer :: i,id,i1

  call diskdat_section('setup')

  ! data
  rmin = 1.d0
  aflux = 0.d0
  flux  = 0.d0
  sflux = 0.d0

  ! read in parameters
  call reapar('mctherm.par',atname,dustname,i1)
  ! occult not in current documentation
  occult = 1

  ! wavelength stuff
  numin=1.2398d0/5000.d0        ! lambda = 5000microns
  ! numax=min(32.*tstar/11605.d0,24.6d0)
  numax=1.2398d0/0.01d0         ! lambda = 0.01microns
  print*,'wavelength range (microns)',1.2398d0/numax,1.2398d0/numin
  nurat=numax/numin
  lnurat=log(numax/numin)
  nunorm=nurat**(0.5d0/dble(nfreq))-nurat**(-0.5d0/dble(nfreq))
  print*,'nunorm',nunorm
  ! wavelength stuff for image cube
  ! numincube=1.2398d0/wavemaxcube        
  ! numaxcube=1.2398d0/wavemincube        
  ! print*,'image cube wavelength range (microns)',1.2398d0/numaxcube,1.2398d0/numincube
  ! nuratcube=numaxcube/numincube
  ! lnuratcube=log(numaxcube/numincube)
  ! nunormcube=nuratcube**(0.5d0/dble(nfreqcube))-nuratcube**(-0.5d0/dble(nfreqcube))

  ! convert units from AU to rstar.
  ! autors=214.94d0/rstar
  autors=au2cm/rsol/rstar
  print*,'AU to Rstar conversion',autors
  ! first calculate rmine based on dust sublimation temperature
  ! 2004/05/13. change definition of rmine, rmind = minimum radius
  ! in units of dust destruction radius.  (so you don't have to
  ! calculate it each time).
  tsub=1600.d0
  ! tsub=300.d0
  ! tsub=150.d0
  ! rminsub=0.5d0*(tsub/tstar)**(-2.4d0)
  ! empirical fit
  ! rminsub=(tsub/tstar)**(-2.085d0)
  rminsub=(tsub/tstar)**(-2.10d0)
  print*,'dust destruction radius ',rminsub
  ! including finite source effect (but not limb darkening)
  ! optically thin approximation
  rmine1=sqrt(1.d0/(1.d0-(1.d0-2.d0*(tsub/tstar)**5)**2))
  ! print*,'rmine including finite source',rmine1
  ! rmine=rmine1
  ! if (rmine.lt.1.d0) rmine=1.
  print*,'dust destruction radius (au) ',rminsub/autors
  ! print*,'including finite source effect',rmine1/autors

  if (irminsub.eq.1) then
     rmine=rminsub*rmine
     rmind=rminsub*rmind_in
     rddust=rminsub*rddust
     rmin_sg=rminsub*rmin_sg
     ! rmaxd=rmaxd+rminsub
  else if (irminsub.eq.2) then
     rmine=rmine*autors
     rmind=rmind*autors
     rddust=rddust*autors
     rmin_sg=rmin_sg*autors
  end if
  rspiral1=rspiral1*autors
  rspiral2=rspiral2*autors

  print*,'rmine,rmind,rddust ',rmine,rmind,rddust

  ! IZMIN=1
  if (czmin.eq.'RSUB') then
     ! calculate scale height based on Tsub at Rsub
     muc=2.3d0
     rcgs=rminsub*rstar*rsol
     mcgs=massc*msol
     honr=sqrt(k*Tsub/Gn/mcgs*rcgs/muc/mH)
     print*,'h on r at Rsub',honr
     print*,'h (rsub) ',honr*rminsub
     ! zmin=honr*rminsub/rminsub**b
     ! new, 20050514 scale zmin
     do id=1,ndg
        if(is_disk(id))  then
			z1(id)=z1(id)*honr*rminsub/rminsub**b(id)
			z1gap(id)=z1gap(id)*honr*rminsub/rminsub**bgap(id)
		endif
     end do
     print*,'h_0 based on h(r_sub)',z1
     print*,'assumes Tsub=',Tsub,'    beta =',b
  else if (czmin.eq.'R100') then
     do id=1,ndg
        if(is_disk(id)) then
			z1(id) = z1(id) * autors / (100.*autors)**b(id)
			z1gap(id)=z1gap(id)*autors/(100.*autors)**bgap(id)
		endif
     end do
  end if

  rmax=rmax*autors
  rmaxd=rmaxd*autors
  rgapd1=rgapd1*autors
  if (rgapd1.lt.rmin) rgapd1=rmin
  rgapd2=rgapd2*autors
  rgape1=rgape1*autors
  rgape2=rgape2*autors
  rmaxi=rmaxi*autors        !image size
  warplength=warplength*autors
  gaplength=gaplength*autors
  rimlength=rimlength*autors
  rimcurvelength=rimcurvelength*autors
  gapcurvelength=gapcurvelength*autors
  rchole=rchole*autors
  rc=rc*autors
  zbub1=zbub1*autors
  zbub2=zbub2*autors
  radexp=radexp*autors
  radmisalign=radmisalign*autors
  ! note: aperture in radius, not diameter

  call diskdat_write('nap',nap,'number of apertures')
  write(12,*) 'aperture sizes in AU, Rstar'
  print*, 'aperture sizes in AU, Rstar'
  if (nap.eq.1) then
     dela=0.d0
     aperture2(1) = apmin
  else
     forall(i=1:nap) aperture2(i) = 10.**(dble(i-1)/dble(nap-1) * (log10(apmax) - log10(apmin)) + log10(apmin))
  end if
  do i=1,nap
     write(12,*) aperture2(i),aperture2(i)*autors
     print*, aperture2(i),aperture2(i)*autors
     aperture2(i)=aperture2(i)*autors
     ! square the aperture
     aperture2(i)=aperture2(i)**2.d0
  end do
  z01=z01*autors
  z02=z02*autors
  ! zflowmin=zflowmin*autors
  npsav=np+npout
  print*,'npsav',npsav

  do id=1,2
  if (rddust(id).gt.rmaxd(id)) then
     print*,'ERROR, rddust > rmaxd'
     print*,'rmind,rmaxd (AU)',rddust(id)/autors,rmaxd(id)/autors
     print*,'stopping program'
     stop
  end if

  if (rmaxd(id).gt.rmax) then
     print*,'ERROR, rmaxd > rmax'
     print*,'resetting rmax=rmaxd = (AU)',rmaxd(id)/autors
     rmax=rmaxd(id)*1.001
  end if
  enddo

  rmaxdmax=max(rmaxd(1),rmaxd(2))
  rddustmin=min(rddust(1),rddust(2))

  ! for now, calculate scattering in old way, h-g
  idust=1

  call diskdat_write('nxhst',nxhst,'image size')
  call diskdat_write('nfreq',nfreq,'number of frequencies')
  call diskdat_write('nmu',nmu,'number of theta angles')
  call diskdat_write('nph',nph,'number of phi angles')

  call diskdat_write('rate',rate,'')
  call diskdat_write('ihole',ihole,'')
  call diskdat_write('thet1',thet1,'cavity opening angle')
  call diskdat_write('thet2',thet2,'cavity opening angle')
  call diskdat_write('idust',idust,'')
  call diskdat_write('limb',limb,'')
  call diskdat_write('occult',occult,'')
  call diskdat_write('massc',massc,'mass of stellar core')


  ! write(12,*) 'initial envelope inner dust radius, rmine ',rmine
  ! write(12,*) 'disk inner radius, rmind ',rmind
  ! write(12,*) 'disk inner dust radius, rddust ',rddust

  ! don't do this anymore, doing all wavelengths
  ! call waveset(cwav,cfldust,ctherm,nwav,MXWAV)
  ! read in from wavein.dat

  ! checks
  if (ihole.eq.1.and.istream.eq.1.and.ipoly.eq.1) then
     write(6,*) 'ERROR, cant have both a streamline and'
     write(6,*) 'polynomial hole.  pick one (istream,  ipoly)'
     write(6,*) ' '
     go to 666
  end if
  go to 777
666 stop
777 continue

  write(6,*) 'pi,r2p',pi,r2p

  ! modified for temporary multiple-angle peeloff
  thete_arr = thete_arr * r2p/360.d0
  phie_arr  = phie_arr  * r2p/360.d0
  sinte_arr=sin(thete_arr)
  coste_arr=cos(thete_arr)
  sinpe_arr=sin(phie_arr)
  cospe_arr=cos(phie_arr)

  ! g2=g**2
  ! if(g.eq.0.d0) isot=1
  hit=0.d0
  htot=0.d0

  ! limb darkening = 0 for constant, 1 for eddington.
  ! make new table if you want something else
  if (limb.eq.1) call table

  ! calculate density based on input mass
  ! convert r to solar units.

  ! z1=zmin
  ! if you want, input h/r instead of z1
  ! z1=rmaxd*z1/rmaxd**b
  ! print*,'scale height at rstar ',z1
  ! write(12,*) 'scale height at rstar ',z1
  if(rmaxd(1).gt.0.d0.or.rmaxd(2).gt.0.d0) then
     zmax=0.d0
     zmin1=0.d0
     zmin2=0.d0
     do id=1,ndg
        if(is_disk(id)) then
           zmax=zmax+fmass(id)*z1(id)*rmaxd(id)**b(id)
           zmin1=zmin1+fmass(id)*z1(id)*rmin**b(id)
           zmin2=zmin2+fmass(id)*z1(id)*rddust(id)**b(id)
        end if
     end do
  else
     zmax=0.d0
  end if

  if (ienvtype.eq.0) then
 ! rescale fiducial density to 1 rstar
  	rhodens1=rhodens1*(autors)**(envexp)
  	print*,'power law envelope density at 1 rstar is ',rhodens1
  endif
  call diskdat_write &
       & ('zminstar',zmin1,'weighted scale height at rstar')
  call diskdat_write('zmin',zmin2,'rddust')
  call diskdat_write('zmax',zmax,'rmaxd')
  call diskdat_write('zminstar',zmin1/autors,'units of AU at rstar')
  call diskdat_write('zmin',zmin2/autors,'rddust')
  call diskdat_write('zmax',zmax/autors,'rmaxd')

  !call diskdat_write('Fsg',Fsg,'mass fraction of PAH grains')

  !spiral arm.  
  pitch=pitch*pi/180.d0

  !misaligned disk. required diffusion=no
  if (diffusion.and.imisalign.eq.1) then
	print*,'*****'
	print*,'you have diffusion set to yes and misaligned disk set to yes'
	print*,'for misaligned disk, set diffusion to no'
	print*,'I am stopping the code now, just go set diffusion to no and rerun'
	stop
  end if

  return
end subroutine setup
