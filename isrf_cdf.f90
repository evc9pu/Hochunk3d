subroutine isrf_cdf(isrf_int,isrf_Av)

  ! BAW 20070816, read in interstellar radiation field, (ISRF)
  ! make cdf array (cumulative probability distributino).

  use tts_mod, only : par_dir,rstar
  use isrf_mod
  use dust_mod
  use constants, only : rsol
  use opacin_mod, only : kapd
  implicit none

  !integer,parameter :: ncldy = 23
  integer,parameter :: ncldy = 44 
  integer,parameter :: nyancy = 85
  integer,parameter :: next = 126

  real*8 :: lnu(ncldy),nufnu(ncldy),freqin(nyancy),pdfin(nyancy)
  real*8 :: logfreqin(nyancy),logpdfin(nyancy),freqext(next)
  real*8 :: kapext(next)
  real*8 :: lmin,lmax,lam,isrf_int,freq,logpdf,kaptmp,rjnk
  real*8 :: kapvext,ext,isrf_Av,testnu,testwave

  integer :: i,j,nin,id

  real*8,parameter :: c = 2.99792458d14


  ! use ISRF from CLOUDY code (thanks remy!)
  ! tnuism must be log nu, and fnuism log nu fnu
  ! /* interstellar radiation field from Black 1987, "Interstellar Processes"
  ! * table of log nu, log nu*fnu taken from his figure 2 */
  ! /* >>chng 99 jun 14 energy range lowered to 1e-8 ryd */

!   lnu = (/6.00,10.72,11.00,11.23,11.47,11.55,11.85,12.26,12.54 &
!        & ,12.71,13.10,13.64,14.14,14.38,14.63,14.93,15.08,15.36,15.43 &
!        & ,16.25,17.09,18.00,23.00 /)
! 
!   nufnu = (/-16.708,-2.96,-2.47,-2.09,-2.11,-2.34,-3.66,-2.72 &
!        & ,-2.45,-2.57,-3.85,-3.34,-2.30,-1.79,-1.79,-2.34,-2.72 &
!        & ,-2.55,-2.62,-5.68,-6.45,-6.30,-11.3 /)

     ! draine's version from tom
	lnu = (/ &
	      & 11.4771,11.5740,11.6990,11.8751,12.0000, &
	      & 12.1761,12.3010,12.4771,12.5229,12.5740, &
	      & 12.6320,12.6990,12.7782,12.8751,13.0000, &
	      & 13.0792,13.1761,13.3979,13.4771,13.5740, &
	      & 13.6990,13.7782,13.8751,13.9456,14.0969, &
	      & 14.1347,14.2218,14.3979,14.5229,14.6320, &
	      & 14.7368,14.8386,14.9380,15.0792,15.1154, &
	      & 15.1427,15.1549,15.1761,15.2218,15.3218, &
	      & 15.3632,15.4357,15.4771,15.5181 &
	      &/)

	nufnu = (/ &
	         & -4.6696,-4.4737,-4.1013,-3.6490,-3.3799, &
	         & -2.9500,-2.6904,-2.5003,-2.4870,-2.4737, &
	         & -2.5007,-2.5442,-2.6904,-2.9122,-3.2924, &
	         & -3.5452,-3.8459,-3.5764,-3.4868,-3.4604, &
	         & -3.2576,-3.0555,-2.8193,-2.6563,-2.3433, &
	         & -2.2756,-2.1362,-1.9542,-1.9252,-1.9702, &
	         & -2.0637,-2.1854,-2.3470,-2.7405,-2.7218, &
	         & -2.7032,-2.6951,-2.6819,-2.6513,-2.5846, &
	         & -2.5743,-2.6490,-2.8327,-3.0116 &
	         &/)




  ! modified their numbers to get them to agree with their Figure 5

  ! data lnu /10.5,11.85,12.26,12.54
  ! $    ,12.71,13.10,13.64,14.14,14.38,14.63,14.93,15.08,15.36,15.43
  ! $    ,16.25,17.09,18.00,23.00 /

  ! data nufnu /-7.4,-3.66,-2.72
  ! $     ,-2.45,-2.57,-3.85,-3.34,-2.30,-1.79,-1.79,-2.34,-2.72
  ! $     ,-2.55,-2.62,-5.68,-6.45,-6.30,-11.3 /

  ! CLOUDY isrf.  comment out next x lines if using Yancy's
  do i=1,ncldy
     freqin(i)=10.d0**lnu(i)
     pdfin(i)=10.d0**nufnu(i)/freqin(i)
     logfreqin(i)=log10(freqin(i))
     logpdfin(i)=log10(pdfin(i))
     ! print*,'freq,pdf',freqin(i),pdfin(i)
     ! print*,'nu,nufnu',freqin(i),10.d0**nufnu(i)
  end do
  nin=ncldy

  ! c     Yancy's file.  comment out next 6 lines if using Cloudy
  ! note, this file doesn't seem right.doesn't agree with cloudy anyway
  ! open(unit=10,file='../parfiles/isrf.dat',status='unknown')
  ! do i=1,nyancy
  ! read(10,*)  freqin(i),pdfin(i)
  ! logfreqin(i)=log10(freqin(i))
  ! logpdfin(i)=log10(pdfin(i))
  ! end do
  ! close(10)
  ! nin=nyancy

  ! set up frequency array.
  ! lmin=0.01   !min wavelength (microns)
  ! lmax=3.6e4  !max  "
  lmin=c/(freqin(nin)*0.99999999d0)
  lmax=c/(freqin(1)*1.0000001d0)
  print*,'lammin,lammax of input file ',lmin,lmax
  if (lmin.lt.0.05d0) lmin=0.05d0
  if (lmax.gt.5000.d0) lmax=5000.d0
  ! if (lmin.lt.1.2398d0/numax) lmin=1.2398d0/numax
  ! if (lmax.gt.1.2398d0/numin) lmax=1.2398d0/numin
  print*,'lammin,lammax used ',lmin,lmax
  do i=1,nisrf
     lam=lmin*(lmax/lmin)**(dble(i)/dble(nisrf))
     freq_i(nisrf-i+1)=c/lam
     ! print*,'lam,freq',lam,freq_i(nisrf-i+1)
  end do

  print*,freq_i(1),freq_i(nisrf)
  print*,freqin(1),freqin(nin)

  ! read in extinction law
  open(unit=10,file=trim(par_dir)//'/r400_ice095.par' &
       & ,status='unknown')
  do i=1,next
     read(10,*) lam,rjnk,rjnk,kaptmp
     freqext(next-i+1)=c/lam
     kapext(next-i+1)=kaptmp
  end do
  close(10)
  freq=c/0.55d0
  call locate(freqext,next,freq,j)
  kapvext=kapext(j)+(freq-freqext(j))/ &
       & (freqext(j+1)-freqext(j))*(kapext(j+1)-kapext(j))
  print*,'isrf_cdf, kapvext ',kapvext

  print*,'isrf_Av',isrf_Av

  ! interpolate cloudy or yancy's isrf onto our grid
  ! At the same time, interpolate kapext and apply extinction law to isrf
  do i=1,nisrf
     freq=freq_i(i)

     ext=1.d0
     ! calculate extinction factor to apply to pdf
     call locate(freqext,next,freq,j)
     if (j.eq.next) then
        kaptmp=kapext(next)
     else if (j.eq.0) then
        kaptmp=kapext(1)
     else
        kaptmp=kapext(j)+(freq-freqext(j))/ &
             & (freqext(j+1)-freqext(j))*(kapext(j+1)-kapext(j))
     end if
     ext=exp(-kaptmp/kapvext/1.086d0*isrf_Av)
     if (ext.gt.1.d0) then
        print*,'oops! ext, kaptmp,kapvext',ext,kaptmp,kapvext
     end if
     ! print*,'freq,j,freqext(j)',freq,j,freqext(j)
     ! print*,'kaptmp,kapvext,ext',kaptmp,kapvext,ext

     call locate2(freqin,nyancy,nin,freq,j)
     ! print*,freq,freqin(1)
     ! print*,'j',j,nyancy
     ! print*,freq_i(i),freqin(j),freqin(j+1)
     ! pdf_i(i)=pdfin(j)+
     ! $        (freq-freqin(j))/
     ! $        (freqin(j+1)-freqin(j))*
     ! $        (pdfin(j+1)-pdfin(j))
     ! see if log-interpolation works better; yes, much better!
     logpdf=logpdfin(j)+ &
          & (log10(freq)-logfreqin(j))/ &
          & (logfreqin(j+1)-logfreqin(j))* &
          & (logpdfin(j+1)-logpdfin(j))
     pdf_i(i)=10.d0**logpdf*ext
     ! try log sampling in freq.  linear is giving ugly jags. so use nu*Inu
     ! note:  Inu*dnu equals (Inu*nu)log(nu)
     ! this doesn't really improve things so comment back out.
     ! logfreq_i(i)=log10(freq_i(i))
     ! nupdf_i(i)=pdf_i(i)*freq_i(i)
     ! print*,'freq,pdf,nupdf',logfreq_i(i),pdf_i(i),nupdf_i(i)
     ! print*,'freq,pdf',freq_i(i),pdf_i(i)
  end do

  ! Make CDF and log_10(CDF) for emissivity.  Also take logs of wavelengths

  ! integrate pdf to get cdf
  cdf_i(1)=0.d0
  logcdf_i(1)=0.d0
  isrfkap_int=0.d0
  do i=2,nisrf
     cdf_i(i)=cdf_i(i-1)+0.5d0*(pdf_i(i)+pdf_i(i-1)) &
          & *(freq_i(i)-freq_i(i-1))
     ! try log sampling in freq.  linear is giving ugly jags
     ! note:  Inu*dnu equals log(Inu*nu)log(nu)
     ! logcdf_i(i)=logcdf_i(i-1)+0.5d0*(nupdf_i(i)+nupdf_i(i-1))
     ! +        *(logfreq_i(i)-logfreq_i(i-1))
     testwave=2.9979e14/freq_i(i)
     testnu=1.2398d0/testwave
     call opacset(testnu)
     !do id=1,ndg
     !   kapd(id)=kappa(id)*rsol*rstar*kappav(id)
     !end do
 !    print*,'wave,kappa(1),kappa(5)',testwave,kappa(1)*kappav(1),kappa(5)*kappav(5)
   ! if (testwave.ge.0.09.and.testwave.le.1000.0) testint=testint+0.5d0*(pdf_i(i)+pdf_i(i-1)) &
     !     & *(freq_i(i)-freq_i(i-1))*kappa(5)   !note:  this requires you read in draine_opac_new.dat
       isrfkap_int=isrfkap_int+0.5d0*(pdf_i(i)+pdf_i(i-1)) &
           & *(freq_i(i)-freq_i(i-1))*kappa(5)*kappav(5)*(1.d0-albedo(5))
       
  end do
  ! print*,'logcdf_i(nisrf)',logcdf_i(nisrf)
  print*,'cdf_i(nisrf)',cdf_i(nisrf)
  print*,'integral of isrf*kappa(smallgrains)',isrfkap_int
  ! stop

  ! hopefully these will agree better after interpolating to finer
  ! grid.
  !TESTING:  put back old way for comparison
  !isrfkap_int=0.02

  ! save integral for later scaling to stellar flux
  isrf_int=cdf_i(nisrf)
  ! isrf_int=logcdf_i(nisrf)
  ! this is in units of ergs/s/sr

  ! normalize cdf
  cdf_i(1)=0.d0
  logcdf_i(1)=-99.d0
  logfreq_i(1)=0.d0
  do i=2,nisrf
     cdf_i(i)=cdf_i(i)/cdf_i(nisrf)
     logcdf_i(i)=log10(cdf_i(i))
     ! logcdf_i(i)=logcdf_i(i)/logcdf_i(nisrf)   !if using nu_Inu
     ! put freq into jon's crazy unit  (I know, it's probably brilliant)
     freq_i(i)=freq_i(i)*1.2398d0/c
     logfreq_i(i)=log10(freq_i(i))
     ! print*,'isrf freq',freq_i(i)
     ! print*,'logcdf',logcdf_i(i)
  end do

  return
end subroutine isrf_cdf
