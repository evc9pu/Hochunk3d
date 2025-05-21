! history:
! 2012/03/30 (mjw):  MPI stuff.  return unless myid == 0 (otherwise get I/O collisions)
!

module output_mod

  use grid_mod
  use tts_mod
  use opacin_mod
  use out_mod
  use filt_mod
  use log_mod
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

contains

  subroutine output_grid_init(prefix,lastdim)
    implicit none
    character(len=*),intent(in) :: prefix
    integer,intent(in) :: lastdim
    integer :: ir,it,ip

#ifdef MPI
    if (myid > 0) return
#endif

    open(unit=15,file=trim(prefix)//'.unf',form='unformatted')
    write(15) nrg-1,ntg-1,npg-1,lastdim
    write(15) (ravearr(ir)/autors,ir=1,size(ravearr))
    write(15) (thetavearr(it),it=1,size(thetavearr))
    write(15) (phiavearr(ip),ip=1,size(phiavearr))
  end subroutine output_grid_init

  subroutine output_grid_logical(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    logical,intent(in) :: array(nrg,ntg,npg)
    integer :: array_int(nrg,ntg,npg)
    integer :: ir,it,ip

#ifdef MPI
    if (myid > 0) return
#endif

    where(array)
       array_int = 1
    elsewhere
       array_int = 0
    end where
    call output_grid_init(prefix,1)
    write(15) (((array_int(ir,it,ip),ir=1,nrg-1),it=1,ntg-1),ip=1,npg-1)
    call output_grid_close()
  end subroutine output_grid_logical

  subroutine output_grid_integer(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    integer,intent(in) :: array(nrg,ntg,npg)
    integer :: ir,it,ip

#ifdef MPI
    if (myid > 0) return
#endif
    call output_grid_init(prefix,1)
    write(15) (((array(ir,it,ip),ir=1,nrg-1),it=1,ntg-1),ip=1,npg-1)
    call output_grid_close()
  end subroutine output_grid_integer

  subroutine output_grid_real8(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    real*8,intent(in) :: array(nrg,ntg,npg)
    integer :: ir,it,ip

#ifdef MPI
    if (myid > 0) return
#endif
    call output_grid_init(prefix,1)
    write(15) (((array(ir,it,ip),ir=1,nrg-1),it=1,ntg-1),ip=1,npg-1)
    call output_grid_close()
  end subroutine output_grid_real8

  subroutine output_grid_logical_4d(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    logical,intent(in) :: array(nrg,ntg,npg,ndg)
    integer :: array_int(nrg,ntg,npg,ndg)
    integer :: ir,it,ip,id

#ifdef MPI
    if (myid > 0) return
#endif
    where(array)
       array_int = 1
    elsewhere
       array_int = 0
    end where
    call output_grid_init(prefix,ndg)
    write(15) ((((array_int(ir,it,ip,id),ir=1,nrg-1),it=1,ntg-1),ip=1,npg-1),id=1,ndg)
    call output_grid_close()
  end subroutine output_grid_logical_4d

  subroutine output_grid_integer_4d(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    integer,intent(in) :: array(nrg,ntg,npg,ndg)
    integer :: ir,it,ip,id

#ifdef MPI
    if (myid > 0) return
#endif
    call output_grid_init(prefix,ndg)
    write(15) ((((array(ir,it,ip,id),ir=1,nrg-1),it=1,ntg-1),ip=1,npg-1),id=1,ndg)
    call output_grid_close()
  end subroutine output_grid_integer_4d

  subroutine output_grid_real8_4d(prefix,array)
    implicit none
    character(len=*),intent(in) :: prefix
    real*8,intent(in) :: array(nrg,ntg,npg,ndg)
    integer :: ir,it,ip,id

#ifdef MPI
    if (myid > 0) return
#endif
    call output_grid_init(prefix,ndg)
    write(15) ((((array(ir,it,ip,id),ir=1,nrg-1),it=1,ntg-1),ip=1,npg-1),id=1,ndg)
    call output_grid_close()
  end subroutine output_grid_real8_4d

  subroutine output_grid_close()
    implicit none
#ifdef MPI
    if (myid > 0) return
#endif
    close(15)
  end subroutine output_grid_close

  subroutine output_accretion(l_spot, l_disk)
    implicit none
    real*8, intent(in) :: l_spot, l_disk
    print *,'Accretion information is only output in FITS format'
  end subroutine output_accretion

  subroutine output()

    implicit none

    character(len=20) :: beg

    real*8 :: wavarr(nfreq)

    real*8 :: dflux,scatave,cpusec,cpuhrs,cpumin,sec
    real*8 :: totvg
    real*8 :: tave1,tave2,count,fsgtot

    real*8 :: fi(nfreq,nmu,nph,nap,no)
    real*8 :: fq(nfreq,nmu,nph,nap,no)
    real*8 :: fu(nfreq,nmu,nph,nap,no)
    real*8 :: fv(nfreq,nmu,nph,nap,no)

    real*8 :: fie(nfreq,nmu,nph,nap,no)
    real*8 :: fqe(nfreq,nmu,nph,nap,no)
    real*8 :: fue(nfreq,nmu,nph,nap,no)
    real*8 :: fve(nfreq,nmu,nph,nap,no)

    real*8 :: pfi(nfreq,npeel,nap,no)
    real*8 :: pfq(nfreq,npeel,nap,no)
    real*8 :: pfu(nfreq,npeel,nap,no)
    real*8 :: pfv(nfreq,npeel,nap,no)

    real*8 :: pfie(nfreq,npeel,nap,no)
    real*8 :: pfqe(nfreq,npeel,nap,no)
    real*8 :: pfue(nfreq,npeel,nap,no)
    real*8 :: pfve(nfreq,npeel,nap,no)

    integer :: peelid
    character(len=10) :: cphi,ctheta

    integer :: i,ia,io,imin,ihrs,ix,iy,it,ip,id,inub,ir,ift

    integer,parameter :: MXNAM=20
    character(len=25) :: filnam(MXNAM)
    character(len=50) :: filename

    character(len=10) :: suffix(no)

    character(len=100) :: fmt,peel_suffix

    real*8 :: masstot, masssgreg

#ifdef MPI
    if (myid > 0) return
#endif

    if(partial_peeloff) then
       peel_suffix = '_partial'
    else
       peel_suffix = ''
    end if

    suffix(1) = ''
    suffix(2) = '_star'
    suffix(3) = '_disk'
    suffix(4) = '_enve'
    suffix(5) = '_outi'
    suffix(6) = '_dire'
    suffix(7) = '_scat'
    suffix(8) = '_ther'

    ! construct wavelength array
    do inub=1,nfreq
       wavarr(inub)=1.2398d0/(numin*nurat**((real(inub)-0.5)/nfreq))
    end do

    ! STANDARD SEDs

    where(si.eq.0.d0)
       fi=0.  ; fq=0.  ; fu=0.  ; fv=0.
       fie=0. ; fqe=0. ; fue=0. ; fve=0.
    elsewhere
       fi = si/nunorm
       fq = sq/si
       fu = su/si
       fv = sv/si
       fie = si2*si/nunorm
       fqe = sq2
       fue = su2
       fve = sv2
    end where

    do io=1,no

       filename = 'flux'//trim(suffix(io))//'.dat'
       open(unit=11,file=filename,status='replace')

       write(11,*) nfreq,nph,nmu,nap,'   nfreq,nph,nmu,nap'

       write(11,'("#",1X,"Lambda",9X,"I",10X,"I_err",9X,"Q/I",7X,"(Q/I)_err",7X,"U/I",7X,"(U/I)_err",7X,"V/I",7X,"(V/I)_err")')

       do ia=1,nap
          do i=1,nmu
             do ip=1,nph
                do inub=1,nfreq
                   write(11,'(F10.4,8(2X,ES11.4))') wavarr(inub),fi(inub,i,ip,ia,io),fie(inub,i,ip,ia,io),&
                        & fq(inub,i,ip,ia,io),fqe(inub,i,ip,ia,io),fu(inub,i,ip,ia,io),&
                        & fue(inub,i,ip,ia,io),fv(inub,i,ip,ia,io),fve(inub,i,ip,ia,io)
                end do
             end do
          end do
       end do

       close(11)

    end do

    ! PEELOFF SEDs

    if(ipeel.eq.1) then

       where(ti.eq.0.d0)
          pfi=0.
          pfq=0.
          pfu=0.
          pfv=0.
          pfie=0.
          pfqe=0.
          pfue=0.
          pfve=0.
       elsewhere
          pfi = ti/nunorm
          pfq = tq/ti
          pfu = tu/ti
          pfv = tv/ti
          pfie = ti2*ti/nunorm
          pfqe = tq2
          pfue = tu2
          pfve = tv2
       end where

       do io=1,no
       do peelid=1,npeel

       write(cphi,'(F5.1)') phie_arr(peelid)*rad2deg
       write(ctheta,'(F5.1)') thete_arr(peelid)*rad2deg

!          filename='peel_flux'//trim(suffix(io))//trim(peel_suffix)//'.dat'
filename='peel_flux'//trim(suffix(io))//trim(peel_suffix)//'_'//trim(adjustl(ctheta))//'_'//trim(adjustl(cphi))//'.dat'
          open(unit=11,file=filename,status='unknown')

          write(11,*) nfreq,npeel,nap,'   nfreq,npeel,nap'

          write(11,'("#",1X,"Lambda",9X,"I",10X,"I_err",9X,"Q/I",7X,"(Q/I)_err",7X,"U/I",7X,"(U/I)_err",7X,"V/I",7X,"(V/I)_err")')

             do ia=1,nap
                do inub=1,nfreq
                   write(11,'(F10.4,8(2X,ES11.4))') wavarr(inub),pfi(inub,peelid,ia,io),pfie(inub,peelid,ia,io),&
                        & pfq(inub,peelid,ia,io),pfqe(inub,peelid,ia,io),pfu(inub,peelid,ia,io),&
                        & pfue(inub,peelid,ia,io),pfv(inub,peelid,ia,io),pfve(inub,peelid,ia,io)
                end do
             end do
          end do

          close(11)

       end do

    end if



    ! flux
    ! flux=flux+aflux
    ! flux which hits disk
    dflux=(tot/flux)
    ! scattered flux
    if(nscat.gt.0) then
       scatave=nscat/sflux
    else
       scatave=0.
    end if
    sflux=sflux/flux
    ! absorbed flux
    ! aflux=aflux/np
    ! abflux=abflux/np
    ! cpu time
    call cpu_time(cpusec)
    cpuhrs=cpusec/3600.
    ihrs=int(cpuhrs)
    cpumin=(cpuhrs-ihrs)*60.
    imin=int(cpumin)
    sec=(cpumin-imin)*60.

    call diskdat_section('output')
    call diskdat_write('rmine',rmine,'envelope inner radius')
    call diskdat_write('rmind',rmind,'dissk inner radius')
    call diskdat_write('rminsub',rddust,'dust destruction radius')
    call diskdat_write('rmax',rmax,'envelope outer radius')
    call diskdat_write('rmaxd',rmaxd,'disk outer radius')
    call diskdat_write('zmax',zmax,'disk scale height at rmaxd')
    call diskdat_write('rmaxi',rmaxi,'image size')
    call diskdat_write('rmind',rmind,'disk inner radius')
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
    call diskdat_write('killflux',abflux/np,'killed photons (flux)')
    call diskdat_write('ltot/lsol*sumsub/sumall',ltot/lsol*sumsub/sumall,'output luminosity')
    call diskdat_write('warn',count_warnings,'# of warning messages')
    call diskdat_write('error',count_errors,'# of error messages')
    write(12,*) 'cputime: ',ihrs,' hrs ',imin,' min ',sec,' sec'

    if(ipeel.eq.1.and.output_imfilt) then
       write(fmt,'("(",I5.5,"(ES11.4,1X))")') nxhst
       do peelid=1,npeel
          write(cphi,'(F5.1)') phie_arr(peelid)*rad2deg
          write(ctheta,'(F5.1)') thete_arr(peelid)*rad2deg
          do ift=1,nbnd
             beg = 'e_'//trim(filtnames(ift))
             beg = trim(beg)//'_'//trim(adjustl(ctheta))
             beg = trim(beg)//'_'//trim(adjustl(cphi))
             open(unit=21,file=trim(beg)//'_I_img'//trim(peel_suffix)//'.dat')
             open(unit=22,file=trim(beg)//'_Q_img'//trim(peel_suffix)//'.dat')
             open(unit=23,file=trim(beg)//'_U_img'//trim(peel_suffix)//'.dat')
             open(unit=24,file=trim(beg)//'_V_img'//trim(peel_suffix)//'.dat')
             do iy=1,nxhst
                write(21,fmt) ((image_b_i(peelid,ix,iy,ift)),ix=1,nxhst)
                write(22,fmt) ((image_b_q(peelid,ix,iy,ift)),ix=1,nxhst)
                write(23,fmt) ((image_b_u(peelid,ix,iy,ift)),ix=1,nxhst)
                write(24,fmt) ((image_b_v(peelid,ix,iy,ift)),ix=1,nxhst)
             end do
             close(21)
             close(22)
             close(23)
             close(24)
          end do
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

    where(dbig+dsg==0)
       jmean = 0
    elsewhere
       jmean = dble(dsg) / dble(dbig + dsg)
    end where

    call output_grid('dsgratio',jmean)

    !    masstot = sum(densarr*vcell) / msun
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
