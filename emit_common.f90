subroutine emit_common(ir,it,ip,idust2,iter,ii,nphot,nstar)

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
  use output_mod

  implicit none

  integer,intent(inout) :: ir,it,ip,idust2,iter,nphot,nstar,ii(3)
  integer :: im,i,id,ntest
  real*8 :: told,tnew,tave,rnew,sux,suy,suz,dA,T4d,T4s,T4,T4eff,mu,T,cosbnew,sinbnew,sth,imf,nacc,ftest
  real*8,external :: dusttemp,chiR,dBBfreq

  ntest=0
  ftest=1.d0

  if (.not.diffus(ir,it,ip)) then

     if (idust2 == -1) then ! Gas accretion emission
       
        call emit_gasacc()

     else if (is_sg(idust2)) then ! VSG+PAH emission

        ! determine which array to sample from based on mean intensity
        if (iter.eq.1) then
           im=2 ! uratio=1
		   imf=0.d0
        else
!			print*,'iter',iter
           im=imean(ir,it,ip)
           imf=imfract(ir,it,ip)
        end if
        ! im=2
        imave=imave+im
        imct=imct+1

        ! testing
        ! im=0
!		print*,'emit_common imf,imfract',imf,imfract(ir,it,ip)
        call rep_draine(im,nub,imf)

        call opacset(nub)

        call emit_sg() ! emit unpolarized photon in random direction

     else
	
        if (ilucy.eq.0) then
           told=tdust(ir,it,ip,idust2)
           ! tdust(ir,it,ip,idust2)=dusttemp(ir,it,ip,dble(nabs(ir,it,ip,idust2))&
            !    &,tdust(ir,it,ip,idust2),rstar,tstarave,nphot,idust2)
           tdust(ir,it,ip,idust2)=dusttemp(ir,it,ip,(nabs(ir,it,ip,idust2))&
                &,tdust(ir,it,ip,idust2),ltot,nphot,idust2)
           tnew=tdust(ir,it,ip,idust2)
           if (nabs(ir,it,ip,idust2).lt.2) then
              tave=0.25d0*told+0.75d0*tnew
           else
              tave=tnew
           end if
        else
           tave=tdust(ir,it,ip,idust2)
        end if
        call reemit(tave,idust2)

     end if
     
     uy=sint*sinp
     ux=sint*cosp
     uz=cost
     
  else
 
!print*,'ir,it,ip,diffus',ir,it,ip,diffus(ir,it,ip)

     ! nabs in diffusion layer only consists of accretion photons
     ! march through theta array to find theta of first non-diffusion
     ! cell
     ! call emitdiff(ir,it,ip,xp,yp,zp)

     ! find diffusion direction
     if (abs(diffdir(ir,it,ip)) .eq. 1) then
	   !radial direction

        if (diffdir(ir,it,ip).eq. 1) then
           write(*,*) 'error, positive r diffusion not supported'
           stop
        end if

        ! diffuse to front of disk
        do while (diffdir(ir,it,ip) .eq. -1 .and. diffus(ir,it,ip)) ! find diff surface
           ir=ir-1
        end do

        ir=ir+1 ! point to lower wall of 1st diffusion cell
        rnew=rarr(ir)
        on_wall = .true.

        xp=xp*rnew/rtot
        yp=yp*rnew/rtot
        rp=sqrt(xp*xp+yp*yp)
        zp=zp*rnew/rtot

        rtot=rnew
        rsq=rnew*rnew

        sux=-xp/rtot ! unit vector perp to diffusion surf
        suy=-yp/rtot
        suz=-zp/rtot

        call isotrp(sux,suy,suz,ux,uy,uz) ! emit isotropically from surf
        cost=uz
        sint=sqrt(1.d0-cost*cost)
        cosp=ux/sint
        sinp=uy/sint
        phi=atan2(sinp,cosp)

        ! find new disk surface temperature

        dA=2.d0*pi*rarr(ir)**2*(costarr(it)-costarr(it+1)) !assumes 2-D

        T4s=Tstar*Tstar
        T4s=T4s*T4s

 
!!  testing!!  20120213
    !    do id=1,ndg
    !    do id=idust2,idust2
	     do id=1,2

           if(.not.is_sg(id)) then

! 20120213 testing
    !          T4d=tdust(ir,it,ip,id)*tdust(ir,it,ip,id)
              T4d=tdust(ir-ntest,it,ip,id)*tdust(ir-ntest,it,ip,id)
              T4d=T4d*T4d

!	          print*,'tstar,nstar,dA',tstar*11605.,nstar,dA
!	          stop
       !     T4=T4d+4*pi*Rstar**2*T4s/(Lsfrac*Nphot*dA)   Rstar**2 cancels units of dA.  nstar=Lsfrac*nphot
       !   actually it probably should be ntot instead of nstar, right?
              T4=T4d+4.d0*pi*T4s/(dble(nstar)*dA*ftest) !surf temp
              tdust(ir,it,ip,id)=min(1600.d0/11605.d0,T4**0.25d0)

              ! determine temperature in disk interior
              nacc=0
              i=ir
              do while (diffdir(i,it,ip) .eq. -1.and. diffus(ir,it,ip))
                 nacc=nacc+(nabs(i,it,ip,id)) !accumulate accretion photon counter
                 i=i+1
              end do

	!		print*,'nacc',nacc

              T4d=T4 !initialize integration variable to surf temp
              i=ir+1
              do while (diffdir(i,it,ip) .eq. -1.and. diffus(ir,it,ip))
                 T4d=T4d+(3.d0*pi*T4s/(dble(nstar)*dA*ftest))*(nacc)&
                      &*kappav(id)*rsol*rstar*densarr(i-1,it,ip,id)&
                      &*chiR(tdust(i-1,it,ip,id),id)*(rarr(i)-&
                      &rarr(i-1))
    !             print*,	(3.d0*pi*T4s/(dble(nstar)*dA))*(nacc)*kappav(id)*rsol*rstar*densarr(i-1,it,ip,id)*chiR(tdust(i-1,it,ip,id),id)*(rarr(i)-rarr(i-1))**0.25d0*11605.d0
	!			 print*,'t4d',t4d**0.25d0*11605.

    
                 told=tdust(i,it,ip,id)
                 tnew=min(1600.d0/11605.d0,T4d**0.25d0)
 			     tdust(i,it,ip,id)=max(tnew,told)
                 nacc=nacc-(nabs(i-1,it,ip,id))
                 i=i+1
              end do
! stop
           end if

        end do

        ! We are in the cell just inside the diffusion layer. Sample
        ! a dust type to emit from and retrieve T4
        
! testing 20120213 !!!!
   !     call select_dust(ir,it,ip,is_thermal,idust2)
        call select_dust(ir,it,ip,is_thermal_disk,idust2)
        T4=Tdust(ir,it,ip,idust2)**4

        ! get photon temperature using grey atm approximation

        T4eff=4.d0*pi*T4s*(nacc)/(ftest*dA*dble(nstar))

        mu=ux*sux+uy*suy+uz*suz
        T=min((T4+0.75d0*mu*T4eff)**0.25d0,1600.d0/11605.d0)

        if (ilucy.eq.0) then
           nub=dBBfreq(T,idust2)
        else
           nub=BBfreq(T)
           ! nub=lucydustfreq(T,idust2)
        end if

        ir=ir-1 !move pointer to first cell below diffusion layer

        ii(1)=ir

     else if (abs(diffdir(ir,it,ip)) .eq. 2) then

        ! diffuse to top or bottom of disk

        if (diffdir(ir,it,ip) .eq. -2) then
           do while (diffus(ir,it,ip).and. diffus(ir,it,ip))
              it=it-1
           end do
           it=it+1
        else
           do while (diffus(ir,it,ip).and. diffus(ir,it,ip))
              it=it+1
           end do
        end if

        ! get photon position at diffusion layer
        cosb=zp/rtot
        sinb=sqrt(1.d0-cosb*cosb)
        cosbnew=1.00001d0*costarr(it)
        sinbnew=sqrt(1.d0-cosbnew*cosbnew)
        xp=xp*sinbnew/sinb !move photon
        yp=yp*sinbnew/sinb
        rp=sqrt(xp**2+yp**2)
        zp=rtot*cosbnew

        ! calculate mu of disk at position of photon
        mu=costarr(it)
        sth=sintarr(it)

        ! unit vector perpendicular to surface
        if (mu .gt. 0.d0) then
           sux=-mu*xp/sth/rtot
           suy=-mu*yp/sth/rtot
           suz=sth
        else
           sux=mu*xp/sth/rtot
           suy=mu*yp/sth/rtot
           suz=-sth
        end if

        ! get isotropic emission for photon
        call isotrp(sux,suy,suz,ux,uy,uz)
        cost=uz
        sint=sqrt(1.d0-cost*cost)
        cosp=ux/sint
        sinp=uy/sint
        phi=atan2(sinp,cosp)

        ! reassign cost,cosp, etc
        if (costarr(it).lt.0.d0) then
           it=it-1
        end if

        ! find new disk surface temperature
        dA=pi*sth*(rarr(ir+1)**2-rarr(ir)**2) !assumes 2-D

        T4s=Tstar*Tstar
        T4s=T4s*T4s

!!  testing!!  20120213
    !    do id=1,ndg
      !  do id=idust2,idust2
	    do id=1,2

           if(.not.is_sg(id)) then

!!  testing!!  20120213
     !         T4d=tdust(ir,it,ip,id)*tdust(ir,it,ip,id)
               if (diffdir(ir,it,ip) .eq. -2) then
                 T4d=tdust(ir,it-ntest,ip,id)*tdust(ir,it-ntest,ip,id)
               else
                 T4d=tdust(ir,it+ntest,ip,id)*tdust(ir,it+ntest,ip,id)
               end if

               T4d=T4d*T4d

              ! T4=T4d+4*pi*Rstar**2*T4s/(Lsfrac*Nphot*dA)
              T4=T4d+4.d0*pi*T4s/(dble(nstar)*dA*ftest) !surf temp
              Tdust(ir,it,ip,id)=min(1600.d0/11605.d0,T4**0.25d0)

              ! determine temperature in disk interior

              T4d=T4 !initialize integration variable to surf temp

              if (diffdir(ir,it,ip) .eq. -2) then

                 nacc=0
                 i=it
                 do while (diffdir(ir,i,ip) .eq. -2.and. diffus(ir,it,ip))
                    nacc=nacc+(nabs(ir,i,ip,id)) !accumulate accretion photon counter
                    i=i+1
                 end do

                 i=it+1
                 do while (diffdir(ir,i,ip) .eq. -2.and. diffus(ir,it,ip))
                    T4d=T4d+(3.d0*pi*T4s/(dble(nstar)*dA*ftest))*(nacc)&
                         &*kappav(id)*rsol*rstar*densarr(ir,i-1,ip,id)&
                         &*chiR(Tdust(ir,i-1,ip,id),id)&
                         &*0.5d0*(rarr(ir)+rarr(ir+1))&
                         &*(thetarr(i)-thetarr(i-1))
					told=tdust(ir,i,ip,id)
                    tnew=min(1600.d0/11605.d0,T4d**0.25d0)
                    tdust(ir,i,ip,id)=max(told,tnew)
                    nacc=nacc-(nabs(ir,i-1,ip,id))
                    i=i+1
                 end do

              else

                 nacc=0
                 i=it
                 do while (diffdir(ir,i,ip) .eq. 2.and. diffus(ir,it,ip))
                    nacc=nacc+(nabs(ir,i,ip,id)) !accumulate accretion photon counter
                    i=i-1
                 end do

                 i=it-1
                 do while (diffdir(ir,i,ip) .eq. 2.and. diffus(ir,it,ip))
                    T4d=T4d+(3.d0*pi*T4s/(dble(nstar)*dA))*(nacc)&
                         & *kappav(id)*rsol*rstar*densarr(ir,i+1,ip,id)&
                         & *chiR(tdust(ir,i+1,ip,id),id)&
                         & *0.5d0*(rarr(ir)+rarr(ir+1))&
                         & *(thetarr(i+2)-thetarr(i+1))
				   told=tdust(ir,i,ip,id)
	               tnew=min(1600.d0/11605.d0,T4d**0.25d0)
	               tdust(ir,i,ip,id)=max(told,tnew)
                     nacc=nacc-(nabs(ir,i+1,ip,id))
                    i=i-1
                 end do


              end if

           end if

        end do

        ! We are in the cell just inside the diffusion layer. Sample
        ! a dust type to emit from and retrieve T4

  !  TESTING!  20120213
 !       call select_dust(ir,it,ip,is_thermal,idust2)
        call select_dust(ir,it,ip,is_thermal_disk,idust2)
        T4=Tdust(ir,it,ip,idust2)**4

        ! get photon temperature using grey atm approximation

        ! T4eff=3.*Lacc*(Rstar/r)**3*(1.d0-sqrt(Rstar/r))*T4s
        T4eff=4.d0*pi*T4s*(nacc)/(dA*dble(nstar)*ftest)

        mu=ux*sux+uy*suy+uz*suz
        T=min((T4+0.75d0*mu*T4eff)**0.25d0,1600.d0/11605.d0)

        if (ilucy.eq.0) then
           nub=dBBfreq(T,idust2)
        else
           nub=BBfreq(T)
           ! nub=lucydustfreq(T,idust2)
        end if

        if (diffdir(ir,it,ip) .eq. -2) then
           it=it-1
        else
           it=it+1
        end if

        ii(2)=it

     else if (abs(diffdir(ir,it,ip)).eq.3) then
        write(*,*) 'error, phi diffusion not supported'
        stop
     else
        write(*,*) 'error, unknown diffusion direction'
		print*,'diffus,diffdir',diffus(ir,it,ip),diffdir(ir,it,ip)
        call output_grid('diffdir',diffdir)
		call output_grid('diffus',diffus)
		call output_grid('darr',densarr)
		call output_grid('darrtot',sum(densarr,dim=4))
		print*,'ip,ux,uy,uz,ir,it,ip,rtot',ip,ux,uy,uz,it,it,ip,rtot
        stop
     end if

  end if

end subroutine emit_common
