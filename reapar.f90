subroutine reapar(filename,atname,dustname,i1)

! history:
! 2012/03/22 (mjw):  remove call to SET_RANDOM_SEED.  RNG initialization is
!			now done in MAIN (TTS)
!

  use tts_mod
  use opacin_mod
  use out_mod
  use filt_mod, only : filter_dir
  use spot_mod
  use messages, only : set_debug, set_warnings, set_errors
  use random
  use grid_mod, only : nrg, ntg, npg, is_disk, is_envelope, is_cavity, fractmask
  use constants, only : msol

  use configuration

  implicit none

  character(len=*),intent(in) :: filename
  character(len=*),intent(out) :: atname,dustname(ndg)

  character(len=50) :: cshape,cgapd,cgape,cgapedens,cgapddens,envtype,crmin,crimcurve,cgapcurve,cdiskcurv,cmisalign,cspotpar

  integer :: iwarn,idisk_hseq

  real*8 :: fmassd1,fsg(4)
  integer :: i1,id,fractmasktmp(4)

  call load_parameter_file(filename)

  ! set these arrays to zero before reading in 2 of the elements
  rmind_in=0.d0
  rmaxd=0.d0

  ! Preliminaries

  call get_parameter_value('NP',np,'Total number of photons')
  call get_parameter_value('NPMIN',npmin,'Minimum number of photons for Lucy temperature calculation')
  call get_parameter_value('NIMIN',n_iter_min,'Minimum number of iterations for Lucy temperature calculation')
  call get_parameter_value('NIMAX',n_iter_max,'Maxmimum number of iterations for Lucy temperature calculation')
  call get_parameter_value('CPEEL',ipeel,'Use peel-off algorithm for SEDs and images')
  call get_parameter_value('CLUCY',ilucy,'Use the Lucy method for temperature calculation')
  call get_parameter_value('I1',i1,'Random number seed')
  i1 = -abs(i1)
!  call set_random_seed(i1)
  call get_parameter_value('CWARN',iwarn,'Whether to print out all warnings and errors')
  call set_warnings(iwarn)
  call set_errors(iwarn)
  call set_debug(iwarn)
  call get_parameter_value('IWRITE',iwrite,'How often to print statistics')
  call get_parameter_value('DENSINPUT',ihydro,'YES for reading in density grid, NO for calculating analytic')

  npout=0

  call get_parameter_value('PARDIR',par_dir,'Directory containing parameter files')
  filter_dir = par_dir

  ! External illumination
  call get_parameter_value('COUTILLUM',iout,'Use outside illumination')
  call get_parameter_value('ISRF_SCL',isrf_scl,'Scale factor for the ISRF')
  call get_parameter_value('ISRF_AV',isrf_av,'Extinction of the ISRF')

  ! SGs
  call get_parameter_value('RMIN_SG',rmin_sg,'minimum radius of SGs/vsgs')

  ! Central source properties
  call get_parameter_value('CLIMB',limb,'Limb-darkening for central source')
  call get_parameter_value('CPLANCKST',iplanckst,'Planck function for stellar emission')
  call get_parameter_value('ATNAME',atname,'atmosphere file')
  atname = trim(par_dir)//'/'//trim(atname)
  call get_parameter_value('RSTAR',rstar,'Stellar radius in Solar radii')
  call get_parameter_value('TSTAR',tstar,'Blackbody temperature of central star (K)')
  call get_parameter_value('MASSC',massc,' Mass of central star (for TSC properties)')

  ! Disk properties
  call get_parameter_value('MASSD',massd,'Disk mass in solar masses')
  call get_parameter_value('CRMIN',crmin,'Units of Rmin:  Rstar, Rsub, or AU')

  
  select case(trim(crmin))
  case('RSTAR')
    write(*,'("RMIND, RMINE, RMINSG are in units of Rstar")')
    irminsub=0
  case('RSUB')
      write(*,'("RMIND, RMINE, RMINSG are in units of Rsub")')
      irminsub=1
    case('AU')
      write(*,'("RMIND, RMINE, RMINSG are in units of AU")')
      irminsub=2
  case default
     stop "CRMIN should be one of RSTAR/RSUB/AU"
  end select
       
  call get_parameter_value('CZMIN',czmin,'Calculate disk scale height based on HSEQ at Rsub')

  select case(trim(czmin))
  case('RSTAR')
     write(*,'("ZSCALE is the scaleheight in R_star at R_star")')
  case('RSUB')
     write(*,'("ZSCALE is the fraction of HSEQ scaleheight at Rsub")')
  case('R100')
     write(*,'("ZSCALE is the scaleheight in AU at 100AU")')
  case default
     stop "CZMIN should be one of RSTAR/RSUB/R100"
  end select

  call get_parameter_value('FMASSD1',fmassd1,' fraction of mass in disk 1 (Large grains settled disk)')

  fmass = -huge(1.d0)
  z1 = -huge(1.d0)
  a = -huge(1.d0)
  b = -huge(1.d0)

  call get_parameter_value('CDISKCURV',cdiskcurv,'define disk geometry using spherical or cylindrical radius')
  select case(trim(cdiskcurv))
  case('SPH')
    idiskcurve=1
  case('CYL')
    idiskcurve=0
  case default
    stop "CDISKCURV should be either SPH or CYL"
  end select

  call get_parameter_value('CHSEQ',ihseq,'solve for hydrostatic equilibrium?')
  call get_parameter_value('ITER_HSEQ',iter_hseq,'start HSEQ solution after this # of iterations')
!  if (ihseq.eq.1) then
!	if (n_iter_max.lt.8) then
!	  print*,'Since HSEQ=yes, we will do 4 iterations for temperature first'
!	  print*,'I am resetting NIMAX=8'
!	  Print*,'go get a drink of coffee'
!	  n_iter_max=8
 !  endif
 ! endif
! let user choose n_iter_max for now
  call get_parameter_value('DISK_HSEQ',idisk_hseq,'which disk to apply hseq; 1=big grain; 2=small grain; 3=both')
  
  select case(idisk_hseq)
  case(1)
     nhseq1=1
     nhseq2=1
  case(2)
     nhseq1=2
     nhseq2=2
  case(3)
     nhseq1=1
     nhseq2=2
  case default
     stop "DISK_HSEQ should be either 1, 2 or 3"
  end select


  ! gap in disk.  
  call get_parameter_value('CGAPD',igapd,'gap in disk?')
  call get_parameter_value('RGAPD1',rgapd1,'Inner gap radius (AU)')
  call get_parameter_value('RGAPD2',rgapd2,'Outer gap radius (AU)')

 ! spiral arms in disk
  call get_parameter_value('CSPIRAL',ispiral,'spiral warps in disk (YES or NO)')
  call get_parameter_value('PITCH',pitch,'pitch angle of spiral arms (degrees)')
  call get_parameter_value('SN',SN,'parameter that determines width of the arms')
  call get_parameter_value('SW',sw,'fraction of dust mass entrained in the arms')
  call get_parameter_value('RSPIRAL1',RSPIRAL1,'radius where spiral arms start')
  call get_parameter_value('RSPIRAL2',RSPIRAL2,'radius where spiral arms end')
 
  ! Disk warp

  call get_parameter_value('CDISKWARP',idiskwarp,'Include disk warp')
  call get_parameter_value('WARPHEIGHT',warpheight,'Scale for disk warp (units of scale height)')
  call get_parameter_value('WEXP',wexp,'Exponent for azimuthal disk warp (cos**wexp)')
  call get_parameter_value('WARPLENGTH',warplength,'radial falloff factor for warp')

! curved inner rim

  call get_parameter_value('CRIMCURVE',irimcurve,'curved inner rim wall')
  call get_parameter_value('RIMCURVEHEIGHT',rimcurveheight,'Scale for rim wall curve (units of scale height)')
!  rimcurveheight=1.d0
  call get_parameter_value('RIMCURVELENGTH',rimcurvelength,'radial falloff factor for curved rim wall')
  call get_parameter_value('RIMCURVEEXP',rimcurveexp,'affects curvature of rim wall')

  ! puffed up inner rim

  call get_parameter_value('CRIMPUFF',irimpuff,'Include disk puffed-up inner rim')
  call get_parameter_value('RIMHEIGHT',rimheight,'Scale for disk puff (units of scale height)')
  call get_parameter_value('RIMLENGTH',rimlength,'radial falloff factor for puffed inner rim')

! curved gap wall

  call get_parameter_value('CGAPCURVE',igapcurve,'Include curved inner gap wall')
  call get_parameter_value('GAPCURVEHEIGHT',gapcurveheight,'Scale for gap wall curve (units of scale height)')
  call get_parameter_value('GAPCURVELENGTH',gapcurvelength,'AU, radial falloff factor for curved gap wall')
  call get_parameter_value('GAPCURVEEXP',gapcurveexp,'affects curvature of gap wall')
!  igapcurve=0
!  gapcurveheight=1.d0
!  gapcurvelength=0.1


! puffed up outer gap

  call get_parameter_value('CGAPPUFF',igappuff,'Include disk puffed-up outer gap')
  call get_parameter_value('GAPHEIGHT',gapheight,'Scale for gap puff (units of scale height)')
  call get_parameter_value('GAPLENGTH',gaplength,'AU, radial falloff factor for puffed outer gap')

! misaligned inner disk

  call get_parameter_value('CMISALIGN',imisalign,'misaligned inner disk')
  call get_parameter_value('RADMISALIGN',radmisalign,'outer radius of misaligned inner portion of disk (AU)')
  call get_parameter_value('INCMISALIGN',incmisalign,'inclination of the misaligned disk (degrees)')

  ! Large grains settled disk
  call get_parameter_value('DUSTNAME(1)',dustname(1),'dust file, thermal grains [disk 1]')
  call get_parameter_value('DUSTNAME(5)',dustname(5),'dust file, <200A grains [disk 1]')
  call get_parameter_value('FSG(1)',fsg(1),'fraction of the mass in <200A grains')
  call get_parameter_value('ZSCALE(1)',z1(1),'Scale height of disk at Rstar')
  call get_parameter_value('RMIND(1)',rmind_in(1),'Disk outer radius in Rstar or Rsub')
  call get_parameter_value('RMAXD(1)',rmaxd(1),'Disk outer radius in AU')
  call get_parameter_value('A(1)',a(1),'Disk density exponent (~ r^(-a))')
  call get_parameter_value('B(1)',b(1),'Disk scale height exponent (~ r^(b))')
  call get_parameter_value('ZSCALE_GAP(1)',z1gap(1),'Scale height of gap (same units as above)')
  call get_parameter_value('A_GAP(1)',agap(1),'gap density exponent (~ r^(-a))')
  call get_parameter_value('B_GAP(1)',bgap(1),'gap scale height exponent (~ r^(b))')
  call get_parameter_value('RHOSCALE_GAP(1)',rhoscale_gap(1),'ratio of gap to disk density at RGAPD2')
  rhoscale_gap(5)=rhoscale_gap(1)

  ! Small grains higher scale height disk
  call get_parameter_value('DUSTNAME(2)',dustname(2),'dust file, thermal grains [disk 2]')
  call get_parameter_value('DUSTNAME(6)',dustname(6),'dust file, <200A grains [disk 2]')
  call get_parameter_value('FSG(2)',fsg(2),'fraction of the mass in <200A grains')
  call get_parameter_value('ZSCALE(2)',z1(2),'Scale height of disk at Rstar')
  call get_parameter_value('RMIND(2)',rmind_in(2),'Disk outer radius in Rstar or Rsub')
  call get_parameter_value('RMAXD(2)',rmaxd(2),'Disk outer radius in AU')
  call get_parameter_value('A(2)',a(2),'Disk density exponent (~ r^(-a))')
  call get_parameter_value('B(2)',b(2),'Disk scale height exponent (~ r^(b))')
  call get_parameter_value('ZSCALE_GAP(2)',z1gap(2),'Scale height of gap (same units as above)')
  call get_parameter_value('A_GAP(2)',agap(2),'gap density exponent (~ r^(-a))')
  call get_parameter_value('B_GAP(2)',bgap(2),'gap scale height exponent (~ r^(b))')
  call get_parameter_value('RHOSCALE_GAP(2)',rhoscale_gap(2),'ratio of gap to disk density at RGAPD2')
  rhoscale_gap(6)=rhoscale_gap(2)

  rmind=rmind_in
  rddust=rmind_in

  ! SG disks
  z1(5:6) = z1(1:2)
  z1(7:8)=0.d0
  a(5:6) = a(1:2)
  b(5:6) = b(1:2)
  agap(5:6)=agap(1:2)
  bgap(5:6)=bgap(1:2)
  z1gap(5:6)=z1gap(1:2)

  z1(3:4)=0
  a(3:4)=0
  b(3:4)=0
  rho0(3:4)=0

  fmass(1) = fmassd1 * (1. - fsg(1))
  fmass(5) = fmassd1 * fsg(1)
  fmass(2) = (1. - fmassd1) * (1. - fsg(2))
  fmass(6) = (1. - fmassd1) * fsg(2)

  rmind_in(5)=rmind_in(1)
  rmind_in(6)=rmind_in(2)
  rmaxd(5)=rmaxd(1)
  rmaxd(6)=rmaxd(2)

 !  Exponential cutoff to radial disk density profile
  call get_parameter_value('CRADEXP',iradexp,'add radial exponential cutoff to disk density?')
  call get_parameter_value('RADEXP',radexp,'radial exponential cutoff to disk density')

  ! Disk accretion

  call get_parameter_value('CDISKACC',idiskacc,'Include disk accretion luminosity')
  if(idiskacc==1) then
     call get_parameter_value('CALPHA',ialpha,'using alpha parameter to calculate disk Mdot')
     if (ialpha.eq.1) then
        call get_parameter_value('ALPHA/MDOT',alphad,'Disk alpha parameter')
     else
        call get_parameter_value('ALPHA/MDOT',mdotdisk,'Disk accretion rate')
     end if
  end if

  call get_parameter_value('RTRUNC',rtrunc,'Magnetosphere co-rotation radius')
!  print*,'ispot',ispot

  call get_parameter_value('CSPOT',ispot,'emission from star+spot')
  spotflag = ispot
  call get_parameter_value('SPOTLAT',spotlat,'latitude of spot (degrees)')
  spotlon=0.d0
  call get_parameter_value('NSPOT',nspot,'number of spots')
  !  if(nspot>2) stop "ERROR: nspot cannot be larger than 2"

  call get_parameter_value('CTEMP',itemp,'YES to read spot temperature, NO to read Fractional spot size')
  if (itemp.eq.1) then
   call get_parameter_value('TSPOT/FSPOT',tspot,'Temperature of spot')
  else
   call get_parameter_value('TSPOT/FSPOT',fspot,'Fractional area of spot')
  end if

  if (ispot.eq.0) then
    fspot=1.d0
    print*,'no stellar hotspots, so fspot is reset to 1'
  endif

  if (itemp.eq.0.and.fspot.eq.0.d0.and.idiskacc.eq.1) then
     print*,'ERROR: fspot=0 and CDISKACC=yes - 0 area hotspot is a singularity'
     print*,'       maybe you meant fspot=1? '
     stop
  end if


  ! reset some variables if disk mass is 0
  if (massd.eq.0.d0) then
     czmin='RSUB'
     z1=1.d0
     idiskacc=0
     idiskwarp=0
     igapd=0
     ispiral=0
     igappuff=0
     irimpuff=0
     print*,'since disk mass is zero, setting:'
     print*,'CZMIN=RSUB, zscale=1.0, CDISKACC=NO '
  end if

  if (ispiral.eq.1.and.idiskwarp.eq.1) then
    print*,'WARNING:  you have both spiral warps and an inner warp.  do you want that?'
  endif

  ! Envelope properties

  call get_parameter_value('DUSTNAME(3)',dustname(3),'dust file, thermal grains [envelope]')
  call get_parameter_value('FSG(3)',fsg(3),'fraction of the mass in <200A grains [envelope]')
  call get_parameter_value('DUSTNAME(7)',dustname(7),'dust file, <200A grains [envelope]')

  fmass(3) = 1. - fsg(3)
  fmass(7) = fsg(3)

  call get_parameter_value('RMAX',rmax,'Envelope outer radius')
  call get_parameter_value('RMINE',rmine_in,'Envelope inner radius in Rstar or Rsub')
  ! gap in envelope.  

  if (rmax.lt.rmaxd(1).or.rmax.lt.rmaxd(2)) then
	print*,'ERROR:  RMAX must be larger to or equal to RMAXD'
	print*,'STOPPING PROGRAM'
	stop
  endif


  rmine=rmine_in

  call get_parameter_value('ENVTYPE',envtype,'ULRICH or POWLAW')
  call get_parameter_value('RATE',rate,'Mass infall rate for ULRICH env (M_O/yr)')
  call get_parameter_value('RC',rc,'centrifugal radius for ULRICH env')
  call get_parameter_value('RHODENS1',rhodens1,'fiducial density at 1 AU for POWLAW env')
  call get_parameter_value('ENVEXP',envexp,'exponent for POWLAW env')

  if (ENVTYPE.eq.'ULRICH')  then
    ienvtype=1
  else if (ENVTYPE.eq.'POWLAW') then
    ienvtype=0
  else
    print*,'ERROR:  ENVTYPE should be set to ULRICH or POWLAW; stopping program'
    stop
  endif

  ! gap in envelope.  
  call get_parameter_value('CGAPE',cgape,'gap in envelope?')
  call get_parameter_value('RGAPE1',rgape1,'Inner gap radius (AU)')
  call get_parameter_value('RGAPE2',rgape2,'Outer gap radius (AU)')
  call get_parameter_value('CGAPEDENS',cgapedens,' scale density in gap (SCALE) or constant density (CONST)')
  call get_parameter_value('FRACTE',fracte,'gap density scale factor (if CGAPEDENS=SCALE)')
  call get_parameter_value('RHOGAPE',rhogape,'density in envelope gap if CGAPEDENS=CONST')

  if (cgape.eq.'YES') then
    igape=1
    if (cgapedens.eq.'SCALE') then
      igapedens=1
    else if (cgapedens.eq.'CONST') then
      igapedens=0
    else
	  print*,'ERROR:  CGAPEDENS should be set to SCALE or CONST; stopping program'
	  stop
    endif
  else
    igape=0
    igapedens=0
    fracte=1.d0  
  endif

  ! Fractal density variations
  call get_parameter_value('CFRACTAL',ifractal,'YES for fractal density variations, or NO')
  call get_parameter_value('DENSRATIO',densratio,'ratio of background to fractal density')
  call get_parameter_value('FRACTMASK',fractmasktmp,'apply to [disk1,disk2,enve,outflow], 0=no 1=yes')
  call get_parameter_value('IFSEED',ifseed,'random number seed, negative integer')
!  call get_parameter_value('FRACTL',fractl,'1/length scale.') 
  fractl=3.792

  fractmask(1)=fractmasktmp(1)
  fractmask(2)=fractmasktmp(2)
  fractmask(3)=fractmasktmp(3)
  fractmask(4)=fractmasktmp(4)
  fractmask(5)=fractmasktmp(1)
  fractmask(6)=fractmasktmp(2)
  fractmask(7)=fractmasktmp(3)
  fractmask(8)=fractmasktmp(4)
  
  print*,'fractmask',fractmask

  ! Cavity properties

  call get_parameter_value('CHOLE',ihole,'Bipolar cavity')
  call get_parameter_value('DUSTNAME(4)',dustname(4),'dust file, thermal grains [cavity]')
  call get_parameter_value('FSG(4)',fsg(4),'fraction of the mass in <200A grains [cavity]')
  call get_parameter_value('DUSTNAME(8)',dustname(8),'dust file, <200A grains [cavity]')

  fmass(4) = 1. - fsg(4)
  fmass(8) = fsg(4)

  if (ihole==1) then

     call get_parameter_value('CSHAPE',cshape,'')
     if (cshape.eq.'STREAM') then
        ipoly=0
        istream=1
        call get_parameter_value('RCHOLE',rchole,'Streamline hole size in AU')
 !       call get_parameter_value('THETMU0',thetmu0,'Opening angle of streamline wall (deg)')
     else
        ipoly=1
        istream=0
        call get_parameter_value('EX1',ex1,'inner Cavity wall exponent')
        call get_parameter_value('EX2',ex2,'outer Cavity wall exponent')
     end if

     if (cshape.eq.'STREAM'.and.ENVTYPE.eq.'POWLAW') then
        print*,'ERROR: CSHAPE=STREAM and ENVTYPE=POWLAW not allowed'
        print*,'For ENVTYPE=POWLAW, set CSHAPE=POLYN'
        print*,'stopping program'
        stop
     endif

     call get_parameter_value('THET1',thet1,'Opening angle of inner cavity wall')
     call get_parameter_value('Z01',z01,'z-intercept inner cavity at w=0')
     call get_parameter_value('THET2',thet2,'Opening angle of outer cavity wall')
     call get_parameter_value('Z02',z02,'z intercept outer cavity at w=0')

    if (cshape.eq.'STREAM') thetmu0=thet2   !streamline is outer surface, polyn is inner

     call get_parameter_value('EXF',exf,'exponent for cavity density power-law')
     call get_parameter_value('RHOCONST1',rhoconst1,'Coefficient for cavity density')
     call get_parameter_value('RHOCONST2',rhoconst2,'Coefficient for cavity density')

 !    thet2=thet1
 !    ex2=ex1
 !    z02=z01
 !    rhoconst2=rhoconst1

  end if

   call get_parameter_value('RHOAMB',rhoamb,'Ambient density')

  nbub=4.d0
  zbub1=250.d0
  zbub2=300.d0
  buboa=0.0d0

  ! Output parameters

  call get_parameter_value('IMFILT',output_imfilt,'Outputing images convolved with filters')
  call get_parameter_value('IMCUBE',output_imcube,'Outputing multi-wavelength data cubes')
  call get_parameter_value('NXIMG',nxhst,'number of pixels in one dimension of square image ')
 ! call get_parameter_value('IMSPARSE',output_imsparse,'Outputing multi-wavelength sparse matrice')
  output_imsparse=.false.

  ! NXHST=149

  call get_parameter_value('NPEEL',npeel,'Number of peel-off angles for output')
  call get_parameter_value('NAP',nap,'Number of apertures for output')

  allocate(thete_arr(npeel),coste_arr(npeel),sinte_arr(npeel))
  allocate(phie_arr(npeel),cospe_arr(npeel),sinpe_arr(npeel))

  allocate(aperture2(nap))

  call get_parameter_value('RMAXI',rmaxi,'Image half-size in AU')
  call get_parameter_value('APMIN',apmin,'radius of smallest aperture in AU')
  call get_parameter_value('APMAX',apmax,'radius of largest aperture in AU')

  call get_parameter_value('THETE',thete_arr,'Theta angle(s) (deg) of high S/N image(s)/SED(s)')
  call get_parameter_value('PHIE',phie_arr,'Phi angle(s) (deg) of high S/N image(s)/SED(s)')

  phie=0.d0

  if(any(thete_arr == 0.)) stop "ERROR: THETE has to be <> 0"
    
  ! Advanced parameters

!  call get_parameter_value('PARTIALPEEL',partial_peeloff,"Do partial peeling-off (only accretion and scattered photons)")
  call get_parameter_value('DIFFUSION',diffusion,"Whether to use the diffusion approximation in very optically thick regions")

!  call get_parameter_value('SPARSEF',sparse_factor,"Sparse image factor")
!  call get_parameter_value('SPARSED',sparse_dim,"Sparse image size")
!  call get_parameter_value('SPARSERMIN',sparse_rmin,"Sparse image inner image size")
!  call get_parameter_value('SPARSERMAX',sparse_rmax,"Sparse image outer image size")
  sparse_factor=3
  sparse_dim=96
  sparse_rmin=10.
  sparse_rmax=1200.  

  ! GRID

  call get_parameter_value('NRG',nrg,'Number of radial cells')
  call get_parameter_value('NTG',ntg,'Number of theta cells')
  call get_parameter_value('NPG',npg,'Number of phi cells')
  
  ! IMAGES/SEDS
  
  call get_parameter_value('NMU',nmu,'Number of theta bins')
  call get_parameter_value('NPH',nph,'Number of phi bins')
  call get_parameter_value('NFREQ',nfreq,'Number of frequencies')

  do id=1,ndg
     dustname(id) = trim(par_dir)//'/'//trim(dustname(id))
  end do

  print *,'-----------------'
  print *,'Checking fmass'
  print *,'Total disk: ',sum(fmass,mask=is_disk)
  print *,'Total envelope: ',sum(fmass,mask=is_envelope)
  print *,'Total cavity: ',sum(fmass,mask=is_cavity)
  print *,'-----------------'

end subroutine reapar




