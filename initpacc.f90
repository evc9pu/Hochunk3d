! ********************************************************

subroutine initpacc(ii,nphot,nstar,iter)

  ! initialize variables for new accretion photon

  use tts_mod
  use grid_mod
  use stokesvar_mod
  use taunum_mod
  use dust_mod
  use random
  use opacin_mod
  use output_mod
  use constants

  implicit none

  real*8 :: winv
  real*8 :: c2p,s2p,cp,sp,opac
  integer :: ii(3),ir,it,ip,id,nphot,idust2,iter,nstar
  integer :: try

  ! find photon position

  do

     ! only emitting photons in disk.  adjust accretion
     ! luminosity accordingly.
     rp=0.d0
     do while ((rp.ge.rmaxdmax).or.(rp.le.rddustmin))
        winv=ran()
        do while(ran().gt.(1.d0-sqrt(winv)))
           winv=ran()
        end do
        rp=1.d0/winv
     end do

     call rantheta(sp,cp,s2p,c2p)
     xp=rp*cp
     yp=rp*sp
     !20090731, weight zp by different dust-type masses in the disk
     zp = sum(fmass*z1*rp**b, mask=is_disk) * gasdev()

     ! 20090321 BAW, if RMAXD=RMAX, can emit photon at z-height outside of outer radius
     ! so resample
     if (zp**2+rp**2.lt.rmax*rmax) exit

  end do

  sip=1.d0
  sqp=0.d0
  sup=0.d0
  svp=0.d0

  rsq=xp*xp+yp*yp+zp*zp
  rtot=sqrt(rsq)
  call locate(r2arr,nrg,rsq,ir)

  ! if (ntg.gt.1) then
  cosb=zp/rtot
  ! sinb=sqrt(1.d0-cosb**2)
  call locate(costarr,ntg,cosb,it)
  ! else
  ! it=1
  ! end if
  ! if (npg.gt.1) then
  lp=atan2(yp,xp)
  if (lp.lt.0.d0) lp=lp+r2p
  call locate(phiarr,npg,lp,ip)
  ! if (ip.gt.1) then
  ! print*,'error in 2-D code; ip should be 1',ip
  ! end if
  ! else
  ! ip=1
  ! end if
  ii(1)=ir
  ii(2)=it
  ii(3)=ip

  ! randomly sample thermal dust type.
!! testing!!  20120213
 ! if(sum(densarr(ir,it,ip,:), mask=is_thermal) > 0.) then
	  if(sum(densarr(ir,it,ip,:), mask=is_thermal_disk) > 0.) then
  !   call select_dust(ir,it,ip,is_thermal,idust2)
     call select_dust(ir,it,ip,is_thermal_disk,idust2)
! this doesn't make any difference.  I don't get this!
  else
     idust2 = -1 ! indicates emit from gas
  end if

  ! SHOULD WE ALLOW ENVELOPE TYPES TO BE SELECTED?  I don't think so.  it's an accretion photon in the disk

  !     tabulate which emission process is chosen as a function of radius
  if(idust2 > 0) then
     if (is_sg(idust2)) then
        dsg(ir,it,ip)=dsg(ir,it,ip)+1
     else
        dbig(ir,it,ip)=dbig(ir,it,ip)+1
     end if
  end if

if (diffus(ir,it,ip)) then
	!testing!!!!   20120213
!    nabs(ir,it,ip,:)=nabs(ir,it,ip,:)+1
 !okay, this is crazy, but I'm going for it
 opac=kapd(1)*densarr(ir,it,ip,1)+kapd(2)*densarr(ir,it,ip,2)
 nfrac(1)=kapd(1)*densarr(ir,it,ip,1)/opac
 nfrac(2)=kapd(2)*densarr(ir,it,ip,2)/opac
 nabs(ir,it,ip,1)=nabs(ir,it,ip,1)+nfrac(1)
 nabs(ir,it,ip,2)=nabs(ir,it,ip,2)+nfrac(2)
!this doesn't make any difference.  ?????
 else if(idust2 > 0) then
    nabs(ir,it,ip,idust2)=nabs(ir,it,ip,idust2)+1
 end if



  call emit_common(ir,it,ip,idust2,iter,ii,nphot,nstar)

  return
end subroutine initpacc


! **********************************************************
