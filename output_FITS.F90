! history:
! 2012/03/09 (mjw):  MPI stuff.  return unless myid == 0 (otherwise get I/O collisions)
!
! hi there

module output_mod

  use grid_mod
  use tts_mod
  use dust_mod
  use opacin_mod
  use out_mod
  use filt_mod
  use log_mod
  use lib_cfitsio
  use lib_array
  use messages
  use constants
#ifdef MPI
  use ttsre_mpi_mod
#endif

  implicit none
  save

  private
  public :: output_grid
  public :: output
  public :: output_accretion

  interface output_grid
     module procedure output_grid_integer
     module procedure output_grid_logical
     module procedure output_grid_real8
     module procedure output_grid_integer_4d
     module procedure output_grid_logical_4d
     module procedure output_grid_real8_4d
  end interface

  interface fits_write_list
     module procedure fits_write_list4
     module procedure fits_write_list8
  end interface

  integer,parameter :: sp = selected_real_kind(p=6,r=37)
  integer,parameter :: dp = selected_real_kind(p=15,r=307)

contains

  subroutine output_grid_init(unit,prefix,bitpix)
    implicit none
    character(len=*),intent(in) :: prefix
    integer,intent(out) :: unit
    integer,intent(in) :: bitpix
    integer :: naxes(3)
#ifdef MPI
    if (myid > 0) return
#endif

    naxes = (/nrg-1,ntg-1,npg-1/)
    call fits_open_new(unit,trim(prefix)//'.fits.gz',confirm=.false.)
    call fits_write_primary_header(unit,bitpix,naxes,.true.)
  end subroutine output_grid_init

  subroutine output_grid_logical(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    logical,intent(in) :: array(nrg,ntg,npg)
    integer :: array_int(nrg,ntg,npg)
    integer :: unit
#ifdef MPI
    if (myid > 0) return
#endif
    where(array)
       array_int = 1
    elsewhere
       array_int = 0
    end where
    call output_grid_init(unit,prefix,32)
    call fits_write_array(unit,array_int(1:nrg-1,1:ntg-1,1:npg-1))
    call output_grid_close(unit)
  end subroutine output_grid_logical

  subroutine output_grid_integer(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    integer,intent(in) :: array(nrg,ntg,npg)
    integer :: unit
#ifdef MPI
    if (myid > 0) return
#endif
    call output_grid_init(unit,prefix,32)
    call fits_write_array(unit,array(1:nrg-1,1:ntg-1,1:npg-1))
    call output_grid_close(unit)
  end subroutine output_grid_integer

  subroutine output_grid_real8(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    real*8,intent(in) :: array(nrg,ntg,npg)
    integer :: unit
#ifdef MPI
    if (myid > 0) return
#endif
    call output_grid_init(unit,prefix,-32)
    call fits_write_array(unit,array(1:nrg-1,1:ntg-1,1:npg-1))
    call output_grid_close(unit)
  end subroutine output_grid_real8

  subroutine output_grid_init_4d(unit,prefix,bitpix)
    implicit none
    character(len=*),intent(in) :: prefix
    integer,intent(out) :: unit
    integer,intent(in) :: bitpix
    integer :: naxes(4)
#ifdef MPI
    if (myid > 0) return
#endif
    naxes = (/nrg-1,ntg-1,npg-1,ndg/)
    call fits_open_new(unit,trim(prefix)//'.fits.gz',confirm=.false.)
    call fits_write_primary_header(unit,bitpix,naxes,.true.)
  end subroutine output_grid_init_4d

  subroutine output_grid_logical_4d(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    logical,intent(in) :: array(nrg,ntg,npg,ndg)
    integer :: array_int(nrg,ntg,npg,ndg)
    integer :: unit
#ifdef MPI
    if (myid > 0) return
#endif
    where(array)
       array_int = 1
    elsewhere
       array_int = 0
    end where
    call output_grid_init_4d(unit,prefix,32)
    call fits_write_array(unit,array_int(1:nrg-1,1:ntg-1,1:npg-1,1:ndg))
    call output_grid_close(unit)
  end subroutine output_grid_logical_4d

  subroutine output_grid_integer_4d(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    integer,intent(in) :: array(nrg,ntg,npg,ndg)
    integer :: unit
#ifdef MPI
    if (myid > 0) return
#endif
    call output_grid_init_4d(unit,prefix,32)
    call fits_write_array(unit,array(1:nrg-1,1:ntg-1,1:npg-1,1:ndg))
    call output_grid_close(unit)
  end subroutine output_grid_integer_4d

  subroutine output_grid_real8_4d(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    real*8,intent(in) :: array(nrg,ntg,npg,ndg)
    integer :: unit
#ifdef MPI
    if (myid > 0) return
#endif
    call output_grid_init_4d(unit,prefix,-32)
    call fits_write_array(unit,array(1:nrg-1,1:ntg-1,1:npg-1,1:ndg))
    call output_grid_close(unit)
  end subroutine output_grid_real8_4d

  subroutine output_grid_close(unit)
    implicit none
    integer,intent(in) :: unit
#ifdef MPI
    if (myid > 0) return
#endif
    call fits_create_hdu(unit,2)
    call fits_write_list(unit,'r walls','cm',rarr*rstar*rsol)
    call fits_create_hdu(unit,3)
    call fits_write_list(unit,'theta walls','rad',thetarr)
    call fits_create_hdu(unit,4)
    call fits_write_list(unit,'phi walls','rad',phiarr)
    call fits_close(unit)
  end subroutine output_grid_close

  subroutine fits_write_list4(unit,name,units,values)
    use lib_cfitsio
    implicit none
    integer,intent(in) :: unit
    character(len=*),intent(in) :: name,units
    real(sp),intent(in) :: values(:)
#ifdef MPI
    if (myid > 0) return
#endif
    call fits_table_write_header(unit,0,1,(/name/),(/'1E'/),(/units/),name)
    call fits_table_write_column(unit,name,values)
  end subroutine fits_write_list4

  subroutine fits_write_list8(unit,name,units,values)
    use lib_cfitsio
    implicit none
    integer,intent(in) :: unit
    character(len=*),intent(in) :: name,units
    real(dp),intent(in) :: values(:)
#ifdef MPI
    if (myid > 0) return
#endif
    call fits_table_write_header(unit,0,1,(/name/),(/'1E'/),(/units/),name)
    call fits_table_write_column(unit,name,values)
  end subroutine fits_write_list8

  subroutine output_atmos(wav, f_nu, wav_new)
    implicit none
    real(dp),intent(in) :: wav(:), f_nu(:), wav_new(:)
    real(dp) :: nu(size(wav)), f_lambda(size(f_nu))
    real(dp) :: f_lambda_new(size(wav_new)), nu_new(size(wav_new)), f_nu_new(size(wav_new))
    integer :: iw,unit
    real(dp) :: lstar

#ifdef MPI
    if (myid > 0) return
#endif

    nu = c / wav
    f_lambda = f_nu * nu / wav

    do iw = 1, size(wav_new)-1
       f_lambda_new(iw) = integral_log10(wav(size(wav):1:-1), f_lambda(size(wav):1:-1), wav_new(iw+1), &
            & wav_new(iw)) / abs(wav_new(iw+1) - wav_new(iw))
    end do

    f_lambda_new = f_lambda_new / integral_log10(wav_new(size(wav_new):1:-1), f_lambda_new(size(wav_new):1:-1))

    nu_new = c / wav_new
    f_nu_new = f_lambda_new * wav_new / nu_new

    lstar = 4.d0*pi*sigt*(tstar*11605.d0)**4.d0*(rstar*rsol)**2.d0
    f_nu_new = f_nu_new * lstar

    call fits_open_new(unit,'stellar_photosphere_binned.fits.gz',confirm=.false.)
    call fits_write_primary_header(unit,-32,(/0/),.true.)
    call fits_create_hdu(unit,2)
    call fits_table_write_header(unit,0,2,(/"wavelength  ","flux_stellar"/), (/"1D","1D"/), (/"microns","nuFnu  "/), 'stellar')
    call fits_table_write_column(unit,'wavelength',wav_new*1.d6)
    call fits_table_write_column(unit,'flux_stellar',nu_new * f_nu_new)
    call fits_close(unit)

  end subroutine output_atmos

  subroutine output_dust(wavarr)
    ! Output the dust properties to dust_properties.fits.gz, interpolating over wavelength

    implicit none

    ! Input
    real(dp),intent(in) :: wavarr(:)

    ! Interpolation
    integer :: id, iw
    real(dp),allocatable,dimension(:,:) :: opacity, albedo

    ! FITS output
    integer :: unit
    character(len=10) :: names(3), units(3)
    character(len=3) :: formats(3)

#ifdef MPI
    if (myid > 0) return
#endif

    ! Allocate arrays
    allocate(opacity(ndg, nfreq))
    allocate(albedo(ndg, nfreq))

    ! Interpolate in Log10 space
    do id=1,ndg
       do iw =1,nfreq
          opacity(id, iw) = interp1d_log10(lamdust(1:nlambda(id),id),kappad(1:nlambda(id),id),&
               &wavarr(iw)*10000.,bounds_error=.false.) * kappav(id)
          albedo(id, iw) = interp1d(log10(lamdust(1:nlambda(id),id)),adust(1:nlambda(id),id),&
               &log10(wavarr(iw)*10000.),bounds_error=.false.)
       end do
    end do

    ! Output to FITS file

    names(1) = "wavelength"
    names(2) = "opacity"
    names(3) = "albedo"

    formats(1) = "1E"
    write(formats(2),'(I0,"E")') ndg
    write(formats(3),'(I0,"E")') ndg

    units(1) = "microns"
    units(2) = "cm^2/g"
    units(3) = ""

    call fits_open_new(unit,'dust_properties.fits.gz',confirm=.false.)
    call fits_write_primary_header(unit,-32,[0],.true.)
    call fits_create_hdu(unit,2)
    call fits_table_write_header(unit,0,3,names, formats, units, 'DUST')
    call fits_table_write_column(unit,'wavelength',wavarr)
    call fits_table_write_column(unit,'opacity',opacity)
    call fits_table_write_column(unit,'albedo',albedo)
    call fits_close(unit)

  end subroutine output_dust

  subroutine output_accretion(l_spot, l_disk)

    use random
    use spot_mod
    implicit none

    real(dp), intent(in) :: l_spot, l_disk
    ! Disk and spot accretion luminosity

    integer :: ir, unit

    real(dp) :: t0, mdot

    ! Radial functions
    integer,parameter :: n_r = 1000
    real(dp) :: r(n_r), t(n_r), p(n_r), h(n_r)

    ! Compute accretion rate in cgs
    mdot=mdotdisk/3600.d0/24.d0/365.25d0*msol

    do ir=1,n_r

       ! Radius
       r(ir) = 10.**(real(ir-1)/real(n_r-1) * (log10(rmaxdmax)-log10(rtrunc)) + log10(rtrunc))

       ! Probability
       p(ir) = (1.d0-sqrt(1.d0/r(ir)))/r(ir)**2

       ! Scaleheight
       h(ir) = sum(fmass*z1*r(ir)**b, is_disk)

       ! Temperature
       t0 = 3.d0 * gn * (msol*massc) * mdot / &
            &    (8.d0 * pi * sigt * (r(ir)*rstar*rsol)**3.d0)
       t(ir) = ( t0 * (1.d0 - sqrt(1.d0/r(ir))) ) **0.25d0

    end do

! i cannot tell if any of the above is actually returned to the calling routine via the modules, so
! just in case, do it and then return if myid > 0
#ifdef MPI
    if (myid > 0) return
#endif


    ! Output radial functions
    call fits_open_new(unit,'accretion_info.fits.gz',confirm=.false.)
    call fits_write_primary_header(unit,-32,(/0/),.true.)
    call fits_write_keyword(unit, 'LSPOT', l_spot)
    call fits_write_keyword(unit, 'LDISK', l_disk)
    call fits_write_keyword(unit, 'SPOTLON', spotlon)
    call fits_write_keyword(unit, 'SPOTLAT', spotlat)
    call fits_write_keyword(unit, 'SPOTSIZE', thspot)
    call fits_write_keyword(unit, 'SPOTTEMP', tspot)
    call fits_create_hdu(unit,2)
    call fits_table_write_header(unit,0,4,(/"r","t","p","h"/), (/"1D","1D","1D","1D"/), (/"cm","K ","  ","cm"/), 'accretion')
    call fits_table_write_column(unit,'r',r*rstar*rsol)
    call fits_table_write_column(unit,'t',t)
    call fits_table_write_column(unit,'p',p)
    call fits_table_write_column(unit,'h',h)
    call fits_close(unit)

  end subroutine output_accretion

  subroutine output()

    use atmos_interp, only : nhnu, lognu, f

    implicit none

    character(len=20) :: beg

    real*8 :: wavarr(nfreq)

    real*8 :: dflux,scatave,cpusec,cpuhrs,cpumin,sec
    real*8 :: tave1,tave2,count,fsgtot

    real :: fluxcube(nfreq,nmu,nph,nap,no,8)
    real :: peelcube(nfreq,npeel,nap,no,8)

    integer :: peelid
    character(len=10) :: cphi,ctheta

    integer :: unit,naxesf(6),naxesp(5),naxesi(2),naxesc(3),naxes4(4),naxesl(2)

    integer :: imin,ihrs,it,ip,inub,ir,ift

    character(len=10) :: prefix(no)

    real*8 :: masstot,massdisk,masssgreg
    integer :: id
    character(len=100) :: filename,peel_suffix

    if(partial_peeloff) then
       peel_suffix = '_partial'
    else
       peel_suffix = ''
    end if

    prefix(1) = 'TOTAL'
    prefix(2) = 'STELLAR'
    prefix(3) = 'DISK'
    prefix(4) = 'ENVELOPE'
    prefix(5) = 'OUTSIDE'
    prefix(6) = 'DIRECT'
    prefix(7) = 'SCATTERED'
    prefix(8) = 'THERMAL'

    ! construct wavelength array
    do inub=1,nfreq
       wavarr(inub)=1.2398d0/(numin*nurat**((dble(inub)-0.5d0)/dble(nfreq)))
     end do


   ! call output_dust(wavarr)

   ! call output_atmos(1.2398d0/10.**lognu(1:nhnu)*1.e-6, f(1:nhnu), wavarr(1:nfreq)*1.e-6)

! i cannot tell if any of the above is actually returned to the calling routine via the modules, so
! just in case, do it and then return if myid > 0
#ifdef MPI
    if (myid > 0) return
#endif

    ! STANDARD SEDs

    fluxcube = 0.
    where(si > 0.)
       fluxcube(:,:,:,:,:,1) = si/real(nunorm)
       fluxcube(:,:,:,:,:,2) = si2*si/real(nunorm)
       fluxcube(:,:,:,:,:,3) = sq/si
       fluxcube(:,:,:,:,:,4) = sq2
       fluxcube(:,:,:,:,:,5) = su/si
       fluxcube(:,:,:,:,:,6) = su2
       fluxcube(:,:,:,:,:,7) = sv/si
       fluxcube(:,:,:,:,:,8) = sv2
    end where

    naxesf = (/nfreq,nmu,nph,nap,no,8/)

    call fits_open_new(unit,'flux_hypercube.fits.gz',confirm=.false.)
    call fits_write_primary_header(unit,-32,naxesf,.true.)
    call fits_write_array(unit,fluxcube)
    call fits_create_hdu(unit,2)
    call fits_write_list(unit,'WAVELENGTH','microns',wavarr)
    call fits_create_hdu(unit,3)
    call fits_write_list(unit,'APERTURE','AU',sqrt(aperture2)/autors)
    call fits_close(unit)

    ! luminosity sources
    naxesl = (/nfreq,no_init/)
    si_init=si_init/real(nunorm)/2.   !why do I need to divide by 2? it gives the right answer.  shoot.
    call fits_open_new(unit,'lum_sources.fits.gz',confirm=.false.)
    call fits_write_primary_header(unit,-32,naxesl,.true.)
    call fits_write_array(unit,si_init)
    call fits_create_hdu(unit,2)
    call fits_write_list(unit,'WAVELENGTH','microns',wavarr)
    call fits_close(unit)

    ! PEELOFF SEDs

    if(ipeel.eq.1) then

       peelcube=0.
       where(ti>0.)
          peelcube(:,:,:,:,1) = ti/real(nunorm)
          peelcube(:,:,:,:,2) = ti2*ti/real(nunorm)
          peelcube(:,:,:,:,3) = tq/ti
          peelcube(:,:,:,:,4) = tq2
          peelcube(:,:,:,:,5) = tu/ti
          peelcube(:,:,:,:,6) = tu2
          peelcube(:,:,:,:,7) = tv/ti
          peelcube(:,:,:,:,8) = tv2
       end where

       naxesp = (/nfreq,npeel,nap,no,8/)

!beg = 'e_'//trim(filtnames(ift))//'_'//trim(adjustl(ctheta))//'_'//trim(adjustl(cphi))

! call fits_open_new(unit,trim(beg)//'_I_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)

       call fits_open_new(unit,'peel_hypercube'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
       call fits_write_primary_header(unit,-32,naxesp,.true.)
       call fits_write_array(unit,peelcube)
       call fits_create_hdu(unit,2)
       call fits_write_list(unit,'WAVELENGTH','microns',wavarr)
       call fits_create_hdu(unit,3)
       call fits_table_write_header(unit,0,2,(/'THETA','PHI  '/),(/'1E','1E'/),(/'deg','deg'/),'ANGLE')
       call fits_table_write_column(unit,'THETA',thete_arr*rad2deg)
       call fits_table_write_column(unit,'PHI',phie_arr*rad2deg)
       call fits_create_hdu(unit,4)
       call fits_write_list(unit,'APERTURE','AU',sqrt(aperture2)/autors)
       call fits_close(unit)

    end if

    ! flux
    ! flux=flux+aflux
    ! flux which hits disk
    dflux=(tot/flux)
    ! scattered flux
    if(nscat.gt.0.d0) then
       scatave=nscat/sflux
    else
       scatave=0.d0
    end if
    sflux=sflux/flux
    ! absorbed flux
    ! aflux=aflux/np
    ! abflux=abflux/np
    ! cpu time
    call cpu_time(cpusec)
    cpuhrs=cpusec/3600.d0
    ihrs=int(cpuhrs)
    cpumin=(cpuhrs-dble(ihrs))*60.d0
    imin=int(cpumin)
    sec=(cpumin-dble(imin))*60.d0

    call diskdat_section('output')
    call diskdat_write('rmine',rmine,'envelope inner radius')
    call diskdat_write('rmax',rmax,'envelope outer radius')
    call diskdat_write('rmaxd',rmaxd(1),'disk outer radius')
    call diskdat_write('rmaxd',rmaxd(2),'disk outer radius')
    call diskdat_write('zmax',zmax,'')
    call diskdat_write('rmaxi',rmaxi,'')
    call diskdat_write('rmind',rmind(1),'')
    call diskdat_write('rmind',rmind(2),'')
    call diskdat_write('photons',np,'')
    call diskdat_write('taur(1)',taur(1),'')
    call diskdat_write('taur(2)',taur(2),'')
    call diskdat_write('taur(3)',taur(3),'')
    call diskdat_write('taur(4)',taur(4),'')
    call diskdat_write('taur(5)',taur(5),'')
    call diskdat_write('massenv',massenv,'')
    call diskdat_write('massdisk',massdisk,'')
    call diskdat_write('fluxtot',flux,'total flux')
    call diskdat_write('fluxdefr',dflux,'flux which hits disk+envelope (fractional)')
    call diskdat_write('fluxde',dflux*flux,'flux which hits disk+envelope')
    call diskdat_write('fluxsfr',sflux,'scattered flux (fractional)')
    call diskdat_write('fluxs',sflux*flux,'scattered flux')
    call diskdat_write('avescat',scatave,'ave number of scatters in this flux')
    call diskdat_write('starabfr',aflux,'flux which gets absorbed (and reemitted) by star')
    call diskdat_write('starabfr',abflux,'killed photons')
    call diskdat_write('killflux',abflux/dble(np),'killed photons (flux)')
    call diskdat_write('ltot/lsol*sumsub/sumall',ltot/lsol*sumsub/sumall,'output luminosity')
    call diskdat_write('warn',count_warnings,'# of warning messages')
    call diskdat_write('error',count_errors,'# of error messages')
    write(12,*) 'cputime: ',ihrs,' hrs ',imin,' min ',sec,' sec'

    if(ipeel.eq.1) then

       do peelid=1,npeel

          write(cphi,'(F5.1)') phie_arr(peelid)*rad2deg
          write(ctheta,'(F5.1)') thete_arr(peelid)*rad2deg

          ! Broadband images

          if(output_imfilt) then

             naxesi = (/nxhst,nxhst/)

             do ift=1,nbnd

                beg = 'e_'//trim(filtnames(ift))//'_'//trim(adjustl(ctheta))//'_'//trim(adjustl(cphi))

                call fits_open_new(unit,trim(beg)//'_I_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
                call fits_write_primary_header(unit,-32,naxesi,.false.)
                call fits_write_keyword(unit,'RMAXI',rmaxi*rstar*rsol)
                call fits_write_keyword(unit,'NXHST',nxhst)
                call fits_write_keyword(unit,'PHIE',phie_arr(peelid)*rad2deg)
                call fits_write_keyword(unit,'THETE',thete_arr(peelid)*rad2deg)
                call fits_write_array(unit,image_b_i(peelid,:,:,ift))
                call fits_create_hdu(unit,2)
                call fits_write_list(unit,'WAVELENGTH','microns',wavarr)
                call fits_close(unit)

                call fits_open_new(unit,trim(beg)//'_Q_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
                call fits_write_primary_header(unit,-32,naxesi,.false.)
                call fits_write_keyword(unit,'RMAXI',rmaxi*rstar*rsol)
                call fits_write_keyword(unit,'NXHST',nxhst)
                call fits_write_array(unit,image_b_q(peelid,:,:,ift))
                call fits_close(unit)

                call fits_open_new(unit,trim(beg)//'_U_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
                call fits_write_primary_header(unit,-32,naxesi,.false.)
                call fits_write_keyword(unit,'RMAXI',rmaxi*rstar*rsol)
                call fits_write_keyword(unit,'NXHST',nxhst)
                call fits_write_array(unit,image_b_u(peelid,:,:,ift))
                call fits_close(unit)

                call fits_open_new(unit,trim(beg)//'_V_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
                call fits_write_primary_header(unit,-32,naxesi,.false.)
                call fits_write_keyword(unit,'RMAXI',rmaxi*rstar*rsol)
                call fits_write_keyword(unit,'NXHST',nxhst)
                call fits_write_array(unit,image_b_v(peelid,:,:,ift))
                call fits_close(unit)

                call fits_open_new(unit,trim(beg)//'_PF_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
                call fits_write_primary_header(unit,-32,naxesi,.false.)
                call fits_write_keyword(unit,'RMAXI',rmaxi*rstar*rsol)
                call fits_write_keyword(unit,'NXHST',nxhst)
                call fits_write_array(unit,sqrt(image_b_q(peelid,:,:,ift)**2+image_b_u(peelid,:,:,ift)**2+image_b_v(peelid,:,:,ift)**2)*image_b_i(peelid,:,:,ift))
                call fits_close(unit)

             end do

          end if

          ! Multi-wavelength images

          if(output_imcube) then

             naxes4 = (/nxhst,nxhst,nfreq,no/)

             beg = 'e_cube'//'_'//trim(adjustl(ctheta))//'_'//trim(adjustl(cphi))

             call fits_open_new(unit,trim(beg)//'_I_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
             call fits_write_primary_header(unit,-32,naxes4,.false.)
             call fits_write_keyword(unit,'RMAXI',rmaxi*rstar*rsol)
             call fits_write_keyword(unit,'NXHST',nxhst)
             call fits_write_keyword(unit,'PHIE',phie_arr(peelid)*rad2deg)
             call fits_write_keyword(unit,'THETE',thete_arr(peelid)*rad2deg)
             call fits_write_array(unit,image_m_i(peelid,:,:,:,:))
             call fits_create_hdu(unit,2)
             call fits_write_list(unit,'WAVELENGTH','microns',wavarr)
             call fits_close(unit)

             naxesc = (/nxhst,nxhst,nfreq/)

             call fits_open_new(unit,trim(beg)//'_Q_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
             call fits_write_primary_header(unit,-32,naxesc,.false.)
             call fits_write_keyword(unit,'RMAXI',rmaxi*rstar*rsol)
             call fits_write_keyword(unit,'NXHST',nxhst)
             call fits_write_array(unit,image_m_q(peelid,:,:,:))
             call fits_close(unit)

             call fits_open_new(unit,trim(beg)//'_U_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
             call fits_write_primary_header(unit,-32,naxesc,.false.)
             call fits_write_keyword(unit,'RMAXI',rmaxi*rstar*rsol)
             call fits_write_keyword(unit,'NXHST',nxhst)
             call fits_write_array(unit,image_m_u(peelid,:,:,:))
             call fits_close(unit)

             call fits_open_new(unit,trim(beg)//'_V_img'//trim(peel_suffix)//'.fits.gz',confirm=.false.)
             call fits_write_primary_header(unit,-32,naxesc,.false.)
             call fits_write_keyword(unit,'RMAXI',rmaxi*rstar*rsol)
             call fits_write_keyword(unit,'NXHST',nxhst)
             call fits_write_array(unit,image_m_v(peelid,:,:,:))
             call fits_close(unit)

          end if

          ! Sparse images

          if(output_imsparse) then

             beg = 'e_sparse'//'_'//trim(adjustl(ctheta))//'_'//trim(adjustl(cphi))

             call nested_image_write(image_ms_i(peelid),trim(beg)//'_I_img'//trim(peel_suffix)//'.fits.gz')
             call nested_image_write(image_ms_q(peelid),trim(beg)//'_Q_img'//trim(peel_suffix)//'.fits.gz')
             call nested_image_write(image_ms_u(peelid),trim(beg)//'_U_img'//trim(peel_suffix)//'.fits.gz')
             call nested_image_write(image_ms_v(peelid),trim(beg)//'_V_img'//trim(peel_suffix)//'.fits.gz')

          end if

       end do

    end if

!!$    if (iveeg.eq.1) then
!!$       open(unit=21,file=filnam(10),status='unknown')
!!$       write(21,*) '      F         Q/F         U/F           V/F'
!!$       totvg=ivgflux
!!$       write(21,*) (totvg),(qvgflux/totvg),(uvgflux/totvg),(vvgflux/totvg)
!!$       write(21,*)'    errF       errQ/F      err(U/F)      err(V/F)'
!!$       write(21,*) (ivgerr)*(totvg),(qvgerr),(uvgerr),(vvgerr)
!!$       close(21)
!!$       ! if you really want vger images, have to assign them names
!!$       ! in namer.  for now, just assume doing one wavelength.
!!$       open(unit=21,file=filnam(11),status='unknown')
!!$       open(unit=22,file=filnam(12),status='unknown')
!!$       open(unit=23,file=filnam(13),status='unknown')
!!$       open(unit=24,file=filnam(14),status='unknown')
!!$       do it=1,ncvg
!!$          write(21,*) (ivg(it,ip),ip=1,npvg)
!!$          write(22,*) (qvg(it,ip),ip=1,npvg)
!!$          write(23,*) (uvg(it,ip),ip=1,npvg)
!!$          write(24,*) (vvg(it,ip),ip=1,npvg)
!!$       end do
!!$       close(21)
!!$       close(22)
!!$       close(23)
!!$       close(24)
!!$    end if

    call output_grid('nabs',nabs)

	call output_grid('imfract',imfract)

    where(dbig+dsg==0)
       jmean = 0.d0
    elsewhere
       jmean = dble(dsg) / dble(dbig + dsg)
    end where

    call output_grid('fvsg_em',jmean)

  !  masstot = sum(densarr) / msun
  !  massdisk = sum(densarr(:,:,:,1)*vcell + densarr(:,:,:,2)*vcell) / msun
  !  massenv = sum(densarr(:,:,:,3)*vcell) / msun
  !  masssgreg = sum(densarr(:,:,:,4)*vcell) / msun

   massdisk = sum(densarr(:,:,:,1)*vcell + densarr(:,:,:,2)*vcell) / msun
   massenv = sum(densarr(:,:,:,3)*vcell + densarr(:,:,:,4)*vcell) / msun
   masssgreg = sum(densarr(:,:,:,5)*vcell + densarr(:,:,:,6)*vcell + densarr(:,:,:,7)*vcell + densarr(:,:,:,8)*vcell) / msun
   masstot=massdisk+massenv+masssgreg
   fsgtot=masssgreg/masstot

   print*,'masstot,massenv (big grains),massdisk (big grains),masssgreg'
   print*,masstot,massenv,massdisk,masssgreg
   print*, 'fraction of sg/vsg mass ',masssgreg/masstot

   call diskdat_write('masstot',masstot,'total circumstellar mass')
   call diskdat_write('massenv',massenv,'envelope mass without sgs/vsgs')
   call diskdat_write('massdisk',massdisk,'disk mass without sgs/vsgs')
   call diskdat_write('masssgs',masssgreg,'total mass in sgs/vsgs')
   call diskdat_write('fsgtot',fsgtot,'ratio of small grains to total')

    do id=1,ndg
       write(filename,'("tmidplane",I0,".dat")') id
       open(unit=15,file=filename,status='unknown')
       it=(ntg+1)/2
       write(15,'(F9.5)',advance='no') thetarr(it)*rad2deg
       write(15,'(" / theta angle at which these values are quoted")')
       write(15,'(15X,"r(rstar)",9X,"r(au)",11X,"tau",8X,"Tdust",5X,"Tdust2")')
       do ir=1,nrg-1
          write(15,900) ir,ravearr(ir),ravearr(ir)/autors,tauarr(ir),tdust(ir,it,1,id)*11605.d0,tdust2(ir,it,1,id)*11605.d0
       end do
       close(15)
900    format(i6,3(1x,f15.5),2(1x,f10.5))
    end do

    do id=1,ndg
       write(filename,'("tave",I0,".dat")') id
       open(unit=15,file=filename,status='unknown')
       write(15,*) 'r(rstar),r(au),  Tave, Tave2'
       do ir=1,nrg-1
          tave1=0.d0
          tave2=0.d0
          count=0.d0
          do it=1,ntg-1
             do ip=1,npg-1
                tave1=tdust(ir,it,ip,id)+tave1
                tave2=tdust2(ir,it,ip,id)+tave2
                count=count+1.d0
             end do
          end do
          tave1=tave1/count
          tave2=tave2/count
          write(15,901) ir,ravearr(ir),ravearr(ir)/autors,tave1*11605.d0,tave2*11605.d0
       end do
       close(15)
901    format(i10,2(1x,f15.5),2(1x,f10.5))
    end do

  end subroutine output

end module output_mod
