!     *********************************************************

subroutine propagate(iphot,sini,xmax,dum,nphot,idust2,icount,iter)

  use tts_mod
  use grid_mod
  use stokesvar_mod
  use taunum_mod
  use opacin_mod
  use vger_mod
  use tauint_mod
  use peeloff_mod
  use dust_mod
  use random
  implicit none

  real*8 :: xran,xpold,ypold,rpold,rho2,zpold,x,y,xmax
  real*8 :: sipold,tsum,pol2,sipnew,mu,sth,sux,suy,suz,t4d,t4s
  real*8 :: dA,t4,t,t4eff,cosbnew,sinbnew,rnew
  real*8 :: sini(nmu)
  real   :: ran2,dusttemp,dbbfreq,chiR,tnew,told,tave,bbfreq,lucydustfreq
  integer :: ii(3),dum(3),iphot,iscat,k,ia,ix,iy
  integer :: icount,i,ir,it,ip,inub,iint,iintmax,nphot,ntgh,nacc
  integer :: iabs,idust2,id,iter

  integer :: indx(3),i_inter

  integer ::peelid

  logical :: taudone

  iintmax=10000

  ! iintmax=100000000
  xpold=xp
  ypold=yp
  !     rpold=rp
  zpold=zp
  !     already calculated rsq and rtot before call to this routine.

  do i=1,3
     ii(i)=dum(i)
     !     print*,'ii ',ii(i)
  end do


  !     iflag is set to 1 if photon scatters
  iflag=0
  exitflag=0
  aflag=0
  iscat=0
  iabs=0
  iint=0

  !     sample tau
  ! xran=ran()
  ! tau=-log(xran)
  !     integrate over distance until the optical depth equals tau.
  ! print*,'hi,rsq.rtot,zp',rsq,rtot,zp,ii(1),ii(2),ii(3)
  !  call tauint(iphot,tsum,ii(1),ii(2),ii(3),idust2)
  ! print*,'hitau,rsq.rtot,zp',rsq,rtot,zp,ii(1),ii(2),ii(3)
  !     this should be equivalent
  !       print*,'exiting loop'
  ir=ii(1)
  it=ii(2)
  ip=ii(3)
  idust2=dustarr(ir,it,ip)
  if (.not.diffus(ir,it,ip)) then
     tsum=0.d0
     xran=ran()
     tau=-log(xran)
  else
     if (.not.on_wall)   !need to initialize on_wall...
     nabs(ir,it,ip)=nabs(ir,it,ip)-1
  end if

  if(aflag.eq.1) then
     write(6,*) 'shouldnt be here, aflag=1'
     return
  end if

  tot=tot+1.d0

  !     photon scatters until exit exits disk
  !  do while(iint.lt.iintmax)

  exitflag=0

  !do loop over cells
  do while (exitflag.eq.0)

     if (diffuse(ir,it,ip)) then

        nabs(ir,it,ip)=nabs(ir,it,ip)+2 ! photon crosses two walls

        ! find outward fluxes for each wall

        f =                  dtauRoss(ir,it,ip,1)               / &
             (dtauRoss(ir+1,it,ip,1)+dtauRoss(ir,it,ip,1))
        bs(1) = (1.d0-f)*t4dust(ir,it,ip) + f*t4dust(ir+1,it,ip)
        db(1) = 2.*(t4dust(ir+1,it,ip)-t4dust(ir,it,ip))    / &
             (dtauRoss(ir+1,it,ip,1)+dtauRoss(ir,it,ip,1))

        f =                  dtauRoss(ir,it,ip,1)               / &
             (dtauRoss(ir-1,it,ip,1)+dtauRoss(ir,it,ip,1))
        bs(-1) = (1.d0-f)*t4dust(ir,it,ip) + f*t4dust(ir-1,it,ip)
        db(-1) = 2.*(t4dust(ir-1,it,ip)-t4dust(ir,it,ip))   / &
             (dtauRoss(ir-1,it,ip,1)+dtauRoss(ir,it,ip,1))

        f =                  dtauRoss(ir,it,ip,2)               / &
             (dtauRoss(ir,it+1,ip,2)+dtauRoss(ir,it,ip,2))
        bs(2) = (1.d0-f)*t4dust(ir,it,ip) + f*t4dust(ir,it+1,ip)
        db(2) = 2.*(t4dust(ir,it+1,ip)-t4dust(ir,it,ip))    / &
             (dtauRoss(ir,it+1,ip,2)+dtauRoss(ir,it,ip,2))

        f =                  dtauRoss(ir,it,ip,2)               / &
             (dtauRoss(ir,it-1,ip,2)+dtauRoss(ir,it,ip,2))
        bs(-2) = (1.d0-f)*t4dust(ir,it,ip) + f*t4dust(ir,it-1,ip)
        db(-2) = 2.*(t4dust(ir,it-1,ip)-t4dust(ir,it,ip))   / &
             (dtauRoss(ir,it-1,ip,2)+dtauRoss(ir,it,ip,2))

        f =                  dtauRoss(ir,it,ip,3)               / &
             (dtauRoss(ir,it,ip+1,3)+dtauRoss(ir,it,ip,3))
        bs(3) = (1.d0-f)*t4dust(ir,it,ip) + f*t4dust(ir,it,ip+1)
        db(3) = 2.*(t4dust(ir,it,ip+1)-t4dust(ir,it,ip))    / &
             (dtauRoss(ir,it,ip+1,3)+dtauRoss(ir,it,ip,3))

        f=                   dtauRoss(ir,it,ip,3)               / &
             (dtauRoss(ir,it,ip-1,3)+dtauRoss(ir,it,ip,3))
        bs(-3) = (1.d0-f)*t4dust(ir,it,ip) + f*t4dust(ir,it,ip-1)
        db(-3) = 2.*(t4dust(ir,it,ip-1)-t4dust(ir,it,ip))   / &
             (dtauRoss(ir,it,ip-1,3)+dtauRoss(ir,it,ip,3))

        !                db = 0.

        fplus(-1) =             (bs(-1)-2./3.*db(-1)) * env%area(i,  j  ,k  ,-1)
        fplus( 1) = fplus(-1) + (bs( 1)-2./3.*db( 1)) * env%area(i+1,j  ,k  ,1)
        fplus(-2) = fplus( 1) + (bs(-2)-2./3.*db(-2)) * env%area(i,  j  ,k  ,-2)
        fplus( 2) = fplus(-2) + (bs( 2)-2./3.*db( 2)) * env%area(i,  j+1,k  ,2)
        fplus(-3) = fplus( 2) + (bs(-3)-2./3.*db(-3)) * env%area(i,  j  ,k  ,-3)
        fplus( 3) = fplus(-3) + (bs( 3)-2./3.*db( 3)) * env%area(i,  j  ,k+1,3)

        !                fplus(-1) =             env%area(i,  j  ,k  ,1)
        !                fplus( 1) = fplus(-1) + env%area(i+1,j  ,k  ,1)
        !                fplus(-2) = fplus( 1) + env%area(i,  j  ,k  ,2)
        !                fplus( 2) = fplus(-2) + env%area(i,  j+1,k  ,2)
        !                fplus(-3) = fplus( 2) + env%area(i,  j  ,k  ,3)
        !                fplus( 3) = fplus(-3) + env%area(i,  j  ,k+1,3)

        !                ftotal = fplus(-3) + fplus(-2) + fplus(-1) + &
        !                         fplus( 3) + fplus( 2) + fplus( 1)
        ftotal = fplus(3)

        ! find exit wall

        xi = ran_num()*ftotal

        if (xi .lt. fplus(-1)) then ! -r direction
           ncross_minus(ir,it,ip,1)=ncross_minus(ir,it,ip,1)+1
           ir=ir-1
           if (ir==0) then     ! hits the star
              ncross_plus(ir+1,it,ip,1)=ncross_plus(ir+1,it,ip,1)+1
              ir=1
           end if
           if(.not.diffus(ir,it,ip)) then ! moving into normal cell?
              current_wall=1
              p%pos%x = env%x%array(i+1)     ! release photon from surface
              p%pos%y = ran_num(env%y%array(j),env%y%array(j+1))
              p%pos%z = ran_num(env%z%array(k),env%z%array(k+1))
              p%n = surface_direction(-xhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i+1,j,k))
              end if
              p%w = 1.
              tau = random_optical_depth() ! find tau to next interaction

              rnew=rarr(ir+1)
              on_wall = .true.
              !          if (rnew.gt.rtot) print*,'rnew>rtot',rnew,rtot
              !          print*,'eps,rnew,rtot',eps,rnew,rtot
              !          rnew=0.9999999d0*rarr(ir) !move photon into cell above diff layer
              xp=xp*rnew/rtot
              yp=yp*rnew/rtot
              rp=sqrt(xp**2+yp**2)
              zp=zp*rnew/rtot
              !          rsq=zp**2+xp**2+yp**2
              !          rtot=sqrt(rsq)
              rtot=rnew
              rsq=rnew*rnew

              sux=-xp/rtot  !unit vector perp to diffusion surf
              suy=-yp/rtot
              suz=-zp/rtot

              call isotrp(sux,suy,suz,ux,uy,uz) !emit isotropically from surf
              cost=uz
              sint=sqrt(1.d0-cost*cost)
              cosp=ux/sint
              sinp=uy/sint
              phi=atan2(sinp,cosp)


           end if

        else if (xi .lt. fplus(1)) then ! +x direction
           ncross_plus(ir+1,j,k,1)=ncross_plus(ir+1,j,k,1)+1
           i=i+1
           if (i==size(env%x%array)) then
              exit_wall=1
              p%pos%x = env%x%array(i)
              p%pos%y = ran_num(env%y%array(j),env%y%array(j+1))
              p%pos%z = ran_num(env%z%array(k),env%z%array(k+1))
              p%n = surface_direction(xhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i-1,j,k))
              end if
              p%w = 1.
              return
           end if
           if(.not.diffus(ir,it,ip)) then
              current_wall=-1
              p%pos%x = env%x%array(i)
              p%pos%y = ran_num(env%y%array(j),env%y%array(j+1))
              p%pos%z = ran_num(env%z%array(k),env%z%array(k+1))
              p%n = surface_direction(xhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i-1,j,k))
              end if
              p%w = 1.
              tau = random_optical_depth()
           end if

        else if (xi .lt. fplus(-2)) then ! -y direction
           ncross_minus(ir,it,ip,2)=ncross_minus(ir,it,ip,2)+1
           j=j-1
           if (j==0) then
              exit_wall=2
              p%pos%x = ran_num(env%x%array(i),env%x%array(i+1))
              p%pos%y = env%y%array(j+1)
              p%pos%z = ran_num(env%z%array(k),env%z%array(k+1))
              p%n = surface_direction(-yhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i,j+1,k))
              end if
              p%w = 1.
              return
           end if
           if(.not.diffus(ir,it,ip)) then
              current_wall=2
              p%pos%x = ran_num(env%x%array(i),env%x%array(i+1))
              p%pos%y = env%y%array(j+1)
              p%pos%z = ran_num(env%z%array(k),env%z%array(k+1))
              p%n = surface_direction(-yhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i,j+1,k))
              end if
              p%w = 1.
              tau = random_optical_depth()
           end if

        else if (xi .lt. fplus(2)) then ! +y direction
           ncross_plus(ir,it+1,ip,2)=ncross_plus(ir,it+1,ip,2)+1
           j=j+1
           if (j==size(env%y%array)) then
              exit_wall=3
              p%pos%x = ran_num(env%x%array(i),env%x%array(i+1))
              p%pos%y = env%y%array(j)
              p%pos%z = ran_num(env%z%array(k),env%z%array(k+1))
              p%n = surface_direction(yhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i,j-1,k))
              end if
              p%w = 1.
              return
           end if
           if(.not.diffus(ir,it,ip)) then
              current_wall=-2
              p%pos%x = ran_num(env%x%array(i),env%x%array(i+1))
              p%pos%y = env%y%array(j)
              p%pos%z = ran_num(env%z%array(k),env%z%array(k+1))
              p%n = surface_direction(yhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i,j-1,k))
              end if
              p%w = 1.
              tau = random_optical_depth()
           end if

        else if (xi .lt. fplus(-3)) then ! -z direction
           ncross_minus(ir,it,ip,3)=ncross_minus(ir,it,ip,3)+1
           k=k-1
           if (k==0) then
              exit_wall=4
              p%pos%x = ran_num(env%x%array(i),env%x%array(i+1))
              p%pos%y = ran_num(env%y%array(j),env%y%array(j+1))
              p%pos%z = env%z%array(k+1)
              p%n = surface_direction(-zhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i,j,k+1))
              end if
              p%w = 1.
              return
           end if
           if(.not.diffus(ir,it,ip)) then
              current_wall=3
              p%pos%x = ran_num(env%x%array(i),env%x%array(i+1))
              p%pos%y = ran_num(env%y%array(j),env%y%array(j+1))
              p%pos%z = env%z%array(k+1)
              p%n = surface_direction(-zhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i,j,k+1))
              end if
              p%w = 1.
              tau = random_optical_depth()
           end if

        else ! +z direction
           ncross_plus(ir,it,ip+1,3)=ncross_plus(ir,it,ip+1,3)+1
           k=k+1
           if (k==0) then
              exit_wall=5
              p%pos%x = ran_num(env%x%array(i),env%x%array(i+1))
              p%pos%y = ran_num(env%y%array(j),env%y%array(j+1))
              p%pos%z = env%z%array(k)
              p%n = surface_direction(zhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i,j,k-1))
              end if
              p%w = 1.
              return
           end if
           if(.not.diffus(ir,it,ip)) then
              current_wall=-3
              p%pos%x = ran_num(env%x%array(i),env%x%array(i+1))
              p%pos%y = ran_num(env%y%array(j),env%y%array(j+1))
              p%pos%z = env%z%array(k)
              p%n = surface_direction(zhat)
              p%k = unit_vector(p%n.X.p%pos)
              p%I = stokes_vector(1.,0.,0.,0.d0)
              if (single_frequency) then
                 p%nu = nu_0
              else
                 p%nu=bbfreq(env%T_array(i,j,k-1))
              end if
              p%w = 1.
              tau = random_optical_depth()
           end if

        end if

     else

        !do loop over interactions within a cell
        taudone=.true.
        do while (taudone)

           call tauint_cell(iphot,tsum,ii(1),ii(2),ii(3))
           ! if (exitflag.eq.1) exit
           if (tsum.ge.tau) then
              taudone=.true.
              ir=ii(1)
              it=ii(2)
              ip=ii(3)
              idust2=dustarr(ir,it,ip)

              xran=ran()
              if(xran.le.albedo(idust2)) then
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
              else
                 !     photon absorbed, reemitted at longer wavelengths
                 !     set photon origin to disk or envelope
                 if (idust2.eq.1.or.idust2.eq.2) then
                    i_orig=2
                 else
                    i_orig=3
                 end if
                 !    i_orig=2
                 if(albedo(idust2).eq.1.d0)print*,'error! albedo=1, yet photon absorbed'
                 !     print*,'freq before abs',nub

                 nabs(ir,it,ip)=nabs(ir,it,ip)+1
                 if (ilucy.eq.0) then
                    told=tdust(ir,it,ip)
                    tdust(ir,it,ip)=dusttemp(ir,it,ip,dble(nabs(ir,it,ip))&
                         &,tdust(ir,it,ip),rstar,tstarave,nphot,idust2)
                    tnew=tdust(ir,it,ip)
                    if (nabs(ir,it,ip).lt.2) then
                       tave=0.25*told+0.75*tnew
                    else
                       tave=tnew
                    end if
                 else
                    tave=tdust(ir,it,ip)
                    !          print*,tave*11605.d0
                 end if
                 call reemit(tave,idust2)
                 uy=sint*sinp
                 ux=sint*cosp
                 uz=cost

                 iflag=0
                 !     print*,'freq after abs',nub
                 call opacset(nub)
                 wave=1.2398d0/nub

                 idust2=dustarr(ir,it,ip)
                 do id=1,4
                    kapd(id)=kappa(id)*rstar*rsol*kappav(id)
                 end do
                 sipnew=sip
                 if (ipeel.eq.1) then
                    do peelid=1,npeel
                       call peeloff(on_wall,xp,yp,zp,sipnew,sqp,sup,svp&
                            &,cost,sint,cosp,sinp,phi,&
                            &hit,htot,rsq,rtot,&
                            &tsum,ii(1),ii(2),ii(3),idust,idust2,iflag,iphot,peelid)
                    end do
                 end if

                 iint=iint+1
              end if
              xpold=xp
              ypold=yp
              zpold=zp
              !     not sure if we need to do this here, might already have this
              !     rsq=xp**2+yp**2+zp**2
              !     rtot=sqrt(rsq)
              !     rpold=rp
              ux=sint*cosp
              uy=sint*sinp
              uz=cost

              !     testing!!!!!!

              if (iint.gt.iintmax) then
                 abflux=abflux+1
                 print*,'large number of interactions',iint
                 go to 400
              end if

              tsum=0.d0
              xran=ran()
              tau=-log(xran)

           else
              taudone=.false.
           end if
        end do

     end if

  end do

  !  end do

  !     should never be here
  !     print*,' oops, should not be here, propagate, L145'

300 continue
  !     photon exits

  !     bin angle
  !     NOTE:  assuming axisymmetric, and z=-z.
  ! k=int(real(nmu-1)*abs(cost)+1.5d0)
  !     for comparison to jon
  ! k=min(int(real(nmu)*abs(cost)+1),nmu)
  !     20080826 new binning in cost, not combining about midplane
  ! print*,k
  k=int(dble(nmu)*(1-cost)/2.)+1
  if(k.lt.0.or.k.gt.nmu) then
     print*,'cost binning error, k, nmu',k,nmu
     print*, 'cost,sint',cost,sint
     print*, 'iphot',iphot
  end if

  if (phi.lt.0.d0) phi=phi+r2p
  ip=int(nph*phi/r2p)+1
  if(ip.lt.0.or.ip.gt.nph) then
     print*, 'phi binning error, m, nph ',ip,nph
     print*, 'phi ', phi
     print*, 'iphot ',iphot
  end if

  !     old imaging.  (not even used anymore but don't ever want to have
  !     to figure these out again)
  ! if (cost.lt.0.d0) then
  !    x = -ypold*cosp+xpold*sinp
  !    y = -zpold*sini(k)-ypold*u(k)*sinp
  !     &        -xpold*u(k)*cosp
  ! else
  !    x = ypold*cosp-xpold*sinp
  !    y = zpold*sini(k)-ypold*u(k)*sinp
  !     &        -xpold*u(k)*cosp
  ! end if

  !     20080826 imaging with no mirror of cost (not used anymore)
  x = ypold*cosp-xpold*sinp
  y = zpold*sini(k)-ypold*u(k)*sinp-xpold*u(k)*cosp

  !     first, sum fluxes
  rho2=(x**2+y**2)
  if (rho2.lt.aperture2(nap)*1.0001d0) then
     flux=flux+1.
     if (iflag.eq.1) then
        sflux=sflux+1
        nscat=nscat+iscat
     end if
     ! else
     !    print*,'who am I?',iphot
     !     print*,'rho2 bigger than aperture2(nap)',rho2,aperture2(nap)
     !     well, maybe the user wants rho2 bigger than aperture2(nap)
  end if

  !     find frequency bin, and make sure it is inside the limits
  inub=min(max(int(nfreq*log(nub/numin)/lnurat)+1,1),nfreq)
  if (inub.lt.1.or.inub.gt.nfreq) then
     print*,'inub out of limits!,inub,nub',inub,nub
  end if

  if (i_orig.lt.1.or.i_orig.gt.3) print*,'error, IOR wrong',i_orig

  !     loop over apertures, and bin photon into SEDs

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

  do ia=1,nap
     if(rho2 < aperture2(ia)) then

        !    add photon to flux arrays
        si(inub,k,ip,ia,indx)=si(inub,k,ip,ia,indx)+sip
        sq(inub,k,ip,ia,indx)=sq(inub,k,ip,ia,indx)+sqp
        su(inub,k,ip,ia,indx)=su(inub,k,ip,ia,indx)+sup
        sv(inub,k,ip,ia,indx)=sv(inub,k,ip,ia,indx)+svp

        !    add photon to flux error arrays
        si2(inub,k,ip,ia,indx)=si2(inub,k,ip,ia,indx)+sip*sip
        sq2(inub,k,ip,ia,indx)=sq2(inub,k,ip,ia,indx)+sqp*sqp
        su2(inub,k,ip,ia,indx)=su2(inub,k,ip,ia,indx)+sup*sup
        sv2(inub,k,ip,ia,indx)=sv2(inub,k,ip,ia,indx)+svp*svp

        nums(inub,k,ip,ia,indx)=nums(inub,k,ip,ia,indx)+1.d0

       ! aveinc(inub,k,ip,ia,indx)=aveinc(inub,k,ip,ia,indx)+abs(cost)

     end if
  end do

  ! if (iflag.eq.1.and.i_orig.eq.1) scount=scount+1.d0
  if (iflag.eq.1) scount=scount+1.d0

  !     testing**************
  !    removed images except for peeling off
  !     if (cost.ge.0) then
  ! ix=int((x+xmax)/fractx)+1
  ! iy=int((y+xmax)/fractx)+1
  ! if (ix.le.nx.and.iy.le.nx.and.ix.gt.0.d0.and.iy.gt.0.d0) then
  !    image(ix,iy,k)=image(ix,iy,k)+(sip)
  !    imagei(ix,iy,k)=imagei(ix,iy,k)+(sip)
  !    imageq(ix,iy,k)=imageq(ix,iy,k)+(sqp)
  !    imageu(ix,iy,k)=imageu(ix,iy,k)+(sup)
  !    imagev(ix,iy,k)=imagev(ix,iy,k)+(svp)
  !    image2(ix,iy,k)=image2(ix,iy,k)+sip*sip
  !    numi(ix,iy,k)=numi(ix,iy,k)+1
  ! else
  !    icount=icount+1
  ! end if

400 continue

  return
end subroutine propagate


!     *********************************************************



