! ********************************************************

subroutine initpacc(ii,nphot,iter)

  ! initialize variables for new accretion photon

  use tts_mod
  use grid_mod
  use stokesvar_mod
  use taunum_mod
  use dust_mod
  use random
  implicit none

  real*8 :: xran,costi,sinti,cospn,phinew,winv
  real*8 :: T,T4eff,cosbnew,sinbnew,sth,mu,T4d,T4s,T4
  real*8 :: dA,suz,suy,sux,rnew

  real*8 :: c2p,s2p,cp,sp,tnew,told,tave
  real*8 :: dBBfreq,chiR,dusttemp,lucydustfreq

  integer :: ii(3),ir,it,ip,i,nphot,nacc,idust2,iter

  ! find photon position

  do

     ! only emitting photons in disk.  adjust accretion
     ! luminosity accordingly.
     rp=0.d0
     do while ((rp.ge.rmaxd).or.(rp.le.rddust))
        winv=ran()
        do while(ran().gt.(1.d0-sqrt(winv)))
           winv=ran()
        end do
        rp=1.d0/winv
     end do

     call rantheta(sp,cp,s2p,c2p)
     xp=rp*cp
     yp=rp*sp
     zp=z1(id)*rp**b*gasdev(i1)
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
  idust2=dustarr(ir,it,ip)
  if (idust2.eq.0) then
     print*,'idust2=0,ir,it,ip',ir,it,ip
  end if
  ! print*,''
  ! print*,'rsq.rtot,zp',rsq,rtot,zp,ir,it,ip,diffus(ir,it,ip)

  if (diffus(ir,it,ip)) then

  else
     nabs(ir,it,ip)=nabs(ir,it,ip)+1

     if (ilucy.eq.0) then
        told=tdust(ir,it,ip)
        tdust(ir,it,ip)=dusttemp(ir,it,ip,dble(nabs(ir,it,ip))&
             &,tdust(ir,it,ip),tscale,ltot,nphot,idust2)
        tnew=tdust(ir,it,ip)
        if (nabs(ir,it,ip).lt.2) then
           tave=0.25*told+0.75*tnew
        else
           tave=tnew
        end if
     else
        tave=tdust(ir,it,ip)
     end if
     call reemit(tave,idust2)
     uy=sint*sinp
     ux=sint*cosp
     uz=cost

  end if

  on_wall=.false.

  return
end subroutine initpacc


! **********************************************************
