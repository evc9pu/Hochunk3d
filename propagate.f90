subroutine propagate(iphot,sini,dum,nphot,nstar,idust2,iter)

  use tts_mod
  use grid_mod
  use stokesvar_mod
  use taunum_mod
  use opacin_mod
  use tauint_mod
  use peeloff_mod
  use dust_mod
  use random
  use constants

  implicit none

  real*8 :: xran,xpold,ypold,rho2,zpold,x,y,opac
  real*8 :: tsum,pol2,sipnew
  real*8 :: sini(nmu)
  integer :: ii(3),dum(3),iphot,iscat,ia
  integer :: ir,it,ip,inub,iint,iintmax,nphot,nstar
  integer :: iabs,idust2,id,iter,i

  integer :: indx(3),i_inter, ii_prev(3)

  integer :: peelid

  iintmax = 10000
  
  xpold=xp
  ypold=yp
  zpold=zp

  do i=1,3
     ii(i)=dum(i)
  end do

  ! iflag is set to 1 if photon scatters
  iflag=0
  exitflag=0
  aflag=0
  iscat=0
  iabs=0
  iint=0

  ! sample tau
  xran=ran()
  tau=-log(xran)

  ! integrate over distance until the optical depth equals tau.
  call tauint(iphot,tsum,ii(1),ii(2),ii(3))
  ir=ii(1)
  it=ii(2)
  ip=ii(3)

  if(exitflag.eq.1) go to 300
  if(aflag.eq.1) then
     write(6,*) 'shouldnt be here, aflag=1'
     return
  end if

  !  randomly sample dust type.
  call select_dust(ir,it,ip,is_any,idust2)
  
 ! print*,'is_any,idust2',is_any,idust2

  tot=tot+1.d0

  ! photon scatters/absorbs/emits until exit exits disk

  do

     if(ran().le.albedo(idust2)) then ! SCATTERED
        iflag=1
        iscat=iscat+1
        iint=iint+1
        sipnew=sip
        if (ipeel.eq.1) then
           do peelid=1,npeel
              call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp &
                   &,cost,sint,cosp,sinp,phi&
                   &,hit,htot,rsq,rtot&
                   &,tsum,ii(1),ii(2),ii(3),idust,idust2,iflag,iphot,peelid)
           end do
        end if

        pol2=sqp**2+sup**2+svp**2
        if (pol2 .gt. sip**2) then
           print*,'error, P^2, I^2 ', pol2,sip**2
           continue
        end if

        if (idust.lt.2) then
          call stokes(idust2)
        else
           print*,'oops, error, idust is wrong'
           stop
        end if

        pol2=sqp**2+sup**2+svp**2
        if (pol2 .gt. sip**2) then
           print*,'error, P^2, I^2 ', pol2,sip**2
           continue
        end if

     else ! ABSORBED + RE-EMITTED

        ! tabulate which emission process is chosen as a function of radius
        if (is_sg(idust2)) then
           dsg(ir,it,ip)=dsg(ir,it,ip)+1
        else
           dbig(ir,it,ip)=dbig(ir,it,ip)+1
        end if

        ! set photon origin to disk or envelope
        if (idust2.eq.1.or.idust2.eq.2.or.idust2.eq.5.or.idust2.eq.6) then
	       i_orig=2
        else
           i_orig=3
        end if

        if(albedo(idust2).eq.1.d0) stop 'error! albedo=1, yet photon absorbed'

        if(.not.diffus(ir,it,ip)) then
	       nabs(ir,it,ip,idust2)=nabs(ir,it,ip,idust2)+1.d0
	    else
		 ! none of this was in the code before 20120217.  obviously I'm testing
		 ! we don't update nabs for diffusion (at least we didn't in previous versions of code)
		 ! want to calculate nfrac for testing.  20120217
 		   opac=kapd(1)*densarr(ir,it,ip,1)+kapd(2)*densarr(ir,it,ip,2)
		   nfrac(1)=kapd(1)*densarr(ir,it,ip,1)/opac
		   nfrac(2)=kapd(2)*densarr(ir,it,ip,2)/opac
		 ! jon said update nabs by 2 for external photon entering diffusion layer
		 ! however, it's already too hot and this makes it hotter!
		  ! nabs(ir,it,ip,1)=nabs(ir,it,ip,1)+nfrac(1)*2.d0
		  ! nabs(ir,it,ip,2)=nabs(ir,it,ip,2)+nfrac(2)*2.d0
		! okay, how about just updating by 1
		  ! nabs(ir,it,ip,1)=nabs(ir,it,ip,1)+nfrac(1)
		  ! nabs(ir,it,ip,2)=nabs(ir,it,ip,2)+nfrac(2)
		! seems best with none.  
        endif

        call emit_common(ir,it,ip,idust2,iter,ii,nphot,nstar)

        iflag=0

        call opacset(nub)
        wave=1.2398d0/nub

        do id=1,ndg
           kapd(id)=kappa(id)*rstar*rsol*kappav(id)
        end do

        sipnew=sip
        if (ipeel.eq.1.and..not.partial_peeloff) then
           ! don't need to know idust because it's not a scattered photon
           do peelid=1,npeel
              call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp&
                   &,cost,sint,cosp,sinp,phi,&
                   &hit,htot,rsq,rtot,&
                   &tsum,ii(1),ii(2),ii(3),idust,idust2,iflag,iphot,peelid)
           end do
        end if

        iint=iint+1

     end if

     ! sample tau
     tau=-log(ran())

     ! integrate distance until optical depth equals tau.
     xpold=xp
     ypold=yp
     zpold=zp

     ! not sure if we need to do this here, might already have this
     ! rsq=xp**2+yp**2+zp**2
     ! rtot=sqrt(rsq)
     ! rpold=rp

     ux=sint*cosp
     uy=sint*sinp
     uz=cost
     ii_prev = ii
     call tauint(iphot,tsum,ii(1),ii(2),ii(3))
     ir=ii(1)
     it=ii(2)
     ip=ii(3)

     if(exitflag.eq.1) exit

      if(any(ii.ne.ii_prev)) iint = 0.

     if (iint.gt.iintmax) then
        killed(ii(1),ii(2),ii(3)) = killed(ii(1),ii(2),ii(3)) + 1
        abflux=abflux+1
        ! print*,'killing photon (too many interactions)'
        go to 400
     end if

     !  randomly sample dust type.
     call select_dust(ir,it,ip,is_any,idust2)

  end do

300 continue

  ! Bin angle
  ! NOTE:  assuming axisymmetric, and z=-z.
  ! k=int(real(nmu-1)*abs(cost)+1.5d0)
  ! for comparison to jon
  ! k=min(int(real(nmu)*abs(cost)+1),nmu)
  ! 20080826 new binning in cost, not combining about midplane

  it=int(dble(nmu)*(1.d0-cost)/2.d0)+1
  if(it.lt.0.or.it.gt.nmu) then
     print*,'cost binning error, it, nmu',it,nmu
     print*, 'cost,sint',cost,sint
     print*, 'iphot',iphot
  end if

  if (phi.lt.0.d0) phi=phi+r2p
  ip=int(dble(nph)*phi/r2p)+1
  if(ip.lt.0.or.ip.gt.nph) then
     print*, 'phi binning error, m, nph ',ip,nph
     print*, 'phi ', phi
     print*, 'iphot ',iphot
  end if

  ! Old imaging.  (not even used anymore but don't ever want to have
  ! to figure these out again)
  ! if (cost.lt.0.d0) then
  !    x = -ypold*cosp+xpold*sinp
  !    y = -zpold*sini(k)-ypold*u(k)*sinp
  !     &        -xpold*u(k)*cosp
  ! else
  !    x = ypold*cosp-xpold*sinp
  !    y = zpold*sini(k)-ypold*u(k)*sinp
  !     &        -xpold*u(k)*cosp
  ! end if

  ! 20080826 imaging with no mirror of cost (not used anymore)
  x = ypold*cosp-xpold*sinp
  y = zpold*sini(it)-ypold*u(it)*sinp-xpold*u(it)*cosp

  ! first, sum fluxes
  rho2=(x**2+y**2)
  if (rho2.lt.aperture2(nap)*1.0001d0) then
     flux=flux+1.d0
     if (iflag.eq.1) then
        sflux=sflux+1
        nscat=nscat+iscat
     end if
  end if

  ! Find frequency bin, and make sure it is inside the limits
  inub=min(max(int(dble(nfreq)*log(nub/numin)/lnurat)+1,1),nfreq)
  if (inub.lt.1.or.inub.gt.nfreq) then
     print*,'inub out of limits!,inub,nub',inub,nub
  end if

  if (i_orig.lt.1.or.i_orig.gt.5) print*,'error, IOR wrong',i_orig

  ! Loop over apertures, and bin photon into SEDs

  if(iflag==1) then
     i_inter=2
  else
     if(i_orig==1) then
        i_inter=1
     else
        i_inter=3
     end if
  end if

    indx = (/1,1+i_orig,5+i_inter/)
!    this is an array with 3 indices,  the first is 1, the second is the origin of the photon, the third is the interaction type
!     1 = all
!     2 = star origin 
!     3 = disk origin (either accretion or re-emission)
!     4 = envelope origin (re-emission)
!     5 = external illumination origin
!     6 = direct star
!     7 = scattered
!     8 = thermal emission  

  do ia=1,nap
     if(rho2 < aperture2(ia)) then

!  comment out this statement if you want background ISFR radiation that doesn't interact
       if (.not.(iflag.eq.0.and.i_orig.eq.4)) then    !ignore externally illuminated photons that don't interact

           !    add photon to flux arrays
           si(inub,it,ip,ia,indx) = si(inub,it,ip,ia,indx) + real(sip)
           sq(inub,it,ip,ia,indx) = sq(inub,it,ip,ia,indx) + real(sqp)
           su(inub,it,ip,ia,indx) = su(inub,it,ip,ia,indx) + real(sup)
           sv(inub,it,ip,ia,indx) = sv(inub,it,ip,ia,indx) + real(svp)

           !    add photon to flux error arrays
           si2(inub,it,ip,ia,indx) = si2(inub,it,ip,ia,indx) + real(sip*sip)
           sq2(inub,it,ip,ia,indx) = sq2(inub,it,ip,ia,indx) + real(sqp*sqp)
           su2(inub,it,ip,ia,indx) = su2(inub,it,ip,ia,indx) + real(sup*sup)
           sv2(inub,it,ip,ia,indx) = sv2(inub,it,ip,ia,indx) + real(svp*svp)

           nums(inub,it,ip,ia,indx) = nums(inub,it,ip,ia,indx) + 1.

           ! aveinc(inub,it,ip,ia,indx) = aveinc(inub,it,ip,ia,indx) + (cost)

       end if

     end if
  end do

! sum all photons except noninteracting outside illumination photons
if (.not.(iflag.eq.0.and.i_orig.eq.4)) sumsub=sumsub+sip
! sum all photons
  sumall=sumall+sip

  if (iflag.eq.1) scount=scount+1.d0

400 continue

  return

end subroutine propagate


