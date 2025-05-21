subroutine hseq(iter)
!logical function hseq()

  use tts_mod
  use constants
  use opacin_mod
  use grid_mod
  use output_mod
  use out_mod

  implicit none

  integer,intent(in) :: iter

  character(len=4) :: suffix

  real*8 :: 	rad,dr,thet,cost,sint,dt,phi,dp,densd, &
       mugas,sigma2,rmincgs,dsigma,sigma,rho1,const, &
       h2,cost2im1,cost2,dthet,Tbar,sigpw,adens,bdens,zscale,rhodens, &
       sigpwtot,sigmatot,sigvis,sigtest,dcost,testmassd(8)


  integer ::  ir,it,ip,idd,nmid,idust2,i,id

  write(suffix,'("_",I3.3)') iter

  mugas=2.3
  rmincgs=rstar*rsol
  const=k/gn/mH/mugas/massc/msol

 ! print*,'k,gn,mH,mugas,massc,msol',k,gn,mH,mugas,massc,msol

!!  set these elsewhere!
  iviscous=0

   print*,'in hseq'
   sigpwtot=0.d0
   sigmatot=0.d0

   do id=nhseq1,nhseq2    
   if(is_disk(id)) then    
      open(unit=15,file='sig'//char(48+id)//'.dat',status='unknown')
      write(15,*)' r/rstar  sigpw   hold/rstar   sigma  h/rstar tbar'
   !   do ir=2,nrg-1
      ir=1   
      rad=0.5d0*(rarr(ir+1)+rarr(ir+2))  
      do while (ir.lt.nrg-1)
         ir=ir+1   !want to start at ir=2
         rad=0.5d0*(rarr(ir)+rarr(ir+1))
         dr=rarr(ir+1)-rarr(ir)
         h2=const*rad*rmincgs
         if (igapd.eq.1.and.rad.lt.rgapd2.and.rad.gt.rgapd1) then
            adens=agap(id)
            bdens=bgap(id)
            zscale=z1gap(id)
            rhodens=rhogap0(id)
         else
            adens=a(id)
            bdens=b(id)
            zscale=z1(id)
            rhodens=rho0(id)
         endif
         do ip=1,npg-1
	       phi=0.5*(phiarr(ip)+phiarr(ip+1))
	       dp=phiarr(ip+1)-phiarr(ip)	
		!     sigpw=dsqrt(r2p)*zscale*rmincgs*rhodens*rad**(bdens-adens)*(1.d0-dsqrt(1./rad))
		     sigpw=sigarr(ir,ip,id)
         sigpwtot=sigpwtot+sigpw*rad*dr*rmincgs**2*r2p

         if (ntg.gt.1) then     !2- or 3-D atmosphere
            it=(ntg+1)/2
       !     print*,'ir,it,ip,id',ir,it,ip,id
        !    print*,'dens',densarr(ir,it,ip,id)
        !    print*,'tdust2',tdust2(ir,it,ip,id)*11605.d0
            densarr(ir,it,ip,id)=1.d0    !setting midplane density to 1.  okay
            thet=0.5d0*(thetarr(it)+thetarr(it+1))
            cost=cos(thet)
            cost2im1=cost*cost
            nmid=(ntg+1)/2
            do it=(ntg+3)/2,ntg-1   !what is this loop over, lower half?  yes
    !          print*,'densarr before resetting',densarr(ir,it,ip,id)/densarr(ir,nmid,ip,id)
               thet=0.5d0*(thetarr(it)+thetarr(it+1))
               cost=cos(thet)
               sint=sin(thet)
               cost2=cost*cost
               dt=thetarr(it+1)-thetarr(it)
		           sigma2=h2*2.d0*(tdust2(ir,it,ip,id)*tdust2(ir,it-1,ip,id)/(tdust2(ir,it,ip,id)+tdust2(ir,it-1,ip,id)))*11605.d0    
		           ! 2.d0*(tdust2(ir,it,ip,id)*tdust2(ir,it-1,ip,id)/(tdust2(ir,it,ip,id)+tdust2(ir,it-1,ip,id)))*11605.d0  is harmonic mean of T
		           densarr(ir,it,ip,id)=densarr(ir,it-1,ip,id)*tdust2(ir,it-1,ip,id)/tdust2(ir,it,ip,id)    &
		                         *exp(-0.5d0*(cost2-cost2im1)/sigma2)
			  !       print*,'cost2-cost2im1',cost2-cost2im1
			         cost2im1=cost2
        !       print*,'densarr after resetting',densarr(ir,it,ip,id)/densarr(ir,nmid,ip,id)
        		 enddo
        		 dsigma=0.d0
        		 Tbar=0.d0
        		 do it=(ntg+1)/2,ntg-1
        		    dthet=thetarr(it+1)-thetarr(it)
        		    dsigma=dsigma+2.d0*densarr(ir,it,ip,id)*rad*dthet
        		    Tbar=Tbar+2.d0*densarr(ir,it,ip,id)*tdust2(ir,it,ip,id)*rad*dthet  !density weighted ave T
        	   enddo
             Tbar=Tbar/dsigma*11605.d0
      	     sigvis=mdotcgs(id)*mugas*mH*sqrt(gn*massc*msol/(rad*rmincgs))/(3.*pi*alphad*k*Tbar*rad*rmincgs)*  &
      	     (1.-sqrt(1./rad))    !bjorkman 97, this assumes azimuthal symmetry so may not be 
      			                         !appropriate in 3-D
 !     			 print*,'rad,alphad,tbar,alphad*tbar',rad,alphad,Tbar,alphad*Tbar
  !           print*,'sigvis',sigvis
      			 sigmatot=sigmatot+sigvis*rad*dr*rmincgs**2*r2p
       !  		 print*,'sigma viscous, sigmapw',sigvis,sigpw
             if (iviscous.eq.1) then
        	      sigma=sigvis
        	   else
        	      sigma=sigpw
      		   endif
         		 rho1=sigma/(dsigma*rmincgs)
         		          !		 print*,'rho1,',rho1
           	 nmid=(ntg+1)/2
        	   do it=nmid,ntg-1
      		     densarr(ir,it,ip,id)=rho1*densarr(ir,it,ip,id)
      		     densarr(ir,ntg-it,ip,id)=densarr(ir,it,ip,id)
      		   enddo
      		   do it=1,ntg-1
!      		     if (npg.gt.1) then
!     		       do ip=1,npg-1   
!        	        phi=0.5d0*(phiarr(ip)+phiarr(ip+1))
!        		      dp=phiarr(ip+1)-phiarr(ip)
      			      densarr(ir,it,ip,id)=densarr(ir,it,1,id)
 !      	           end do
  !    		     endif
      	       enddo
!  test, rho integrated along arc should equal sigma;  works!
         !   sigtest=0.d0
         !   do it=1,ntg-1
		 !     dcost=costarr(it+1)-costarr(it)	    
		 !     print*,'dcost,rad,densarr(ir,it,1,id)',dcost,rad,densarr(ir,it,1,id)
         !     sigtest=sigtest+abs(densarr(ir,it,1,id)*rad*dcost*rmincgs)   !rad*costarr(it)*rmincgs is ds along theta
         !     print*,'sigtest = sigarr?',sigtest,sigarr(ir,id)
         !   enddo

        	else              !1-D atmosphere
        	  print*,'error, cannot do HSEQ in 1-D'
        	  stop
          endif
        
        	write(15,*) rad,rad/autors,sigpw,z1*rad**b,sigvis,sqrt(h2*tdust2(ir,nmid,1,id))*rad,tbar*11605.	
        	
          !end of radial loop?
          ! print*,sigpwtot/msun,sigmatot/msun
        end do
        end do
            close(15)
         !end of dust loop
   endif
   enddo

 ! normalize densarr to mass of disk
!  testmassd = sum(densarr(:,:,:,1)*vcell + densarr(:,:,:,2)*vcell) / msun
! recalculate mass array
! do id=1,2
! do ir=1,nrg-1
!   do it=1,ntg-1
!     do ip=1,npg-1
!	    densarr(ir,it,ip,id)=densarr(ir,it,ip,id)*massd/testmassd
!		massarr(ir,it,ip,id)=(densarr(ir,it,ip,id)*vcell(ir,it,ip))**0.25d0 !massarr is really mass**0.25
!     enddo
!   enddo
! enddo
! enddo

  do id=1,2
	testmassd(id)=sum(densarr(:,:,:,id)*vcell)/msun
	densarr(:,:,:,id)=densarr(:,:,:,id)*massd*fmass(id)/testmassd(id)
  enddo
    densarr(:,:,:,5)=densarr(:,:,:,1)*fmass(5)/fmass(1)
    densarr(:,:,:,6)=densarr(:,:,:,2)*fmass(6)/fmass(1)

  
  do id=1,ndg
	if (is_disk(id)) massarr(:,:,:,id)=(densarr(:,:,:,id)*vcell)**0.25d0
  enddo


   
   print*,'integral of sigma',sigpwtot/msun
  ! print*,'mdot used in viscous disk,alphad',mdot
  ! print*,'alphad',alphad


  call output_grid('darr' //suffix,densarr)
  call output_grid('darr',densarr)
  call output_grid('darrtot',sum(densarr,dim=4))

  !stop


  !    hseq=.true.
   !  hseq=.false.

   !check for convergence.

! stop
  return

end subroutine hseq
! end function hseq

