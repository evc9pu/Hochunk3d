subroutine diskflux()

  ! calculates flux from disk
  !
  ! history:
  ! 2000/07/28 (mjw):  change filter on rtmp to .gt.1.d0

  use tts_mod
  use stokesvar_mod
  use constants

  implicit none

  real*8 :: wavarr(nfreq), freqarr(nfreq)
  real*8 :: rnorm,rnp,image_norm,image_norm_intensity
  
  integer :: inub

  print*,'in diskflux'


  ! total flux from disk
  rnp=dble(np)

! trying to sum up output luminosity but not getting it right.
!  lumout=sum(si(:,:,:,nap,1))/rnp/nunorm*ltot/r2p/lsun
!  print*,'output luminosity',lumout

  rnorm=rnp/dble(nmu)/dble(nph)/(ltot/lsol)
  image_norm=rnp/(ltot/lsol)
  ! image_norm=rnp
  print*,'ltot check',ltot/lsol

  print*,'luminosity minus noninteracting outside illumination photons', ltot/lsol*sumsub/sumall

 ! intensity units for images
 image_norm_intensity=rnp/ltot*4.d0*pi/1.d17*(2.*rmaxi/autors/dble(nxhst))**2*pc2cm**2/(206265.d0)**2
! =4*pi*d(pc)^2 / pc2cm^2 (pc to cm) / 1.d17 (ergs to MJy) * (2*rmaxi(au)/nxhst/d(pc))^2 (pixel size in arcsec) /206265^2 (arcsec to steradians)

! construct wavelength array
  do inub=1,nfreq
     wavarr(inub)=1.2398d0/(numin*nurat**((dble(inub)-0.5d0)/dble(nfreq)))
     freqarr(inub)=c*1.d4/wavarr(inub)
  end do

                                     
  where(nums > 1.)
     si2 = sqrt((si2-si**2/nums)/(nums-1.)*nums)/si+1./nums
     sq2 = sqrt((sq2-sq**2/nums)/(nums-1.)*nums)/si
     su2 = sqrt((su2-su**2/nums)/(nums-1.)*nums)/si
     sv2 = sqrt((sv2-sv**2/nums)/(nums-1.)*nums)/si
     !  = aveinc / nums
  elsewhere
     si2 = 0.
     sq2 = 0.
     su2 = 0.
     sv2 = 0.
  end where

  si = si / real(rnorm)
  sq = sq / real(rnorm)
  su = su / real(rnorm)
  sv = sv / real(rnorm)


  si_init = si_init / real(image_norm) !initial spectrum, no angle dependence

  if(ipeel==1) then

     where(numt > 1.)
        ti2 = sqrt((ti2*numt/ti**2-1.)/(numt-1.))
        tq2 = sqrt((tq2*numt/tq**2-1.)/(numt-1.))
        tu2 = sqrt((tu2*numt/tu**2-1.)/(numt-1.))
        tv2 = sqrt((tv2*numt/tv**2-1.)/(numt-1.))
     elsewhere
        ti2 = 0.
        tq2 = 0.
        tu2 = 0.
        tv2 = 0.
     end where
  
     ti = ti / real(image_norm)
     tq = tq / real(image_norm)
     tu = tu / real(image_norm)
     tv = tv / real(image_norm)

     if(output_imfilt) then
       
       !first divide q,u,v by i, then normalize i       
       where(image_b_i > 0.)
         image_b_q = image_b_q / image_b_i
         image_b_u = image_b_u / image_b_i
         image_b_v = image_b_v / image_b_i
       elsewhere
         image_b_q = 0.
         image_b_u = 0. 
         image_b_v = 0.
       endwhere
       
       ! now normalize I to intensity units      
        image_b_i = image_b_i / real(image_norm_intensity)
     end if

     if(output_imcube) then
       
       !first divide q,u,v by i, then normalize i       
       where(image_m_i(:,:,:,:,1) > 0.)
         image_m_q = image_m_q / image_m_i(:,:,:,:,1)
         image_m_u = image_m_u / image_m_i(:,:,:,:,1)
         image_m_v = image_m_v / image_m_i(:,:,:,:,1)
       elsewhere
         image_m_q = 0.
         image_m_u = 0. 
         image_m_v = 0.
       endwhere
       
 ! now normalize I to intensity units
         do inub=1,nfreq
           image_m_i(:,:,:,inub,:) = image_m_i(:,:,:,inub,:) / real((freqarr(inub+1)-freqarr(inub))) / real(image_norm_intensity)
        enddo
     end if

     ! need to add nested image normalization

  end if

  ! check constants
  ! rnorm has factor of 2*pi because you made all the photons
  ! come out at one phi.
  ! rnorm=r2p*delnt*real(np)
  ! that was the value for the correct rnorm.  using that
  ! value, the flux at a given angle from the star, if there
  ! is not absorbing disk, is 1/4pi.  We will multiply the
  ! flux by 4pi so that it is normalized to the stellar flux.
  ! rnorm=0.5d0*dmu*real(np)/normstar
  ! rnorm=0.5d0*dmu*dble(np)/(ltot/lsol)

  ! multiply rnorm by 2 if you are doubling photons by
  ! symmetry through z axis (in diskim)
  ! rnorm=rnorm*2.d0
  ! fimage=0.d0
  ! print*,'normalizing flux image (not pol images) by ',rnorm
  ! whenever we deal with bins 1 and nmu, multiply by 2 because
  ! bin is half size of others.

  ! icenter=nx/2+1
  ! do it=1,nmu
  ! do ix=1,nx
  ! do iy=1,nx
  ! standard deviation of imagei array
  ! rtmp=dble(numi(ix,iy,it))
  ! if (image2(ix,iy,it).gt.1.d0) then
  ! print*,'rtmp',rtmp
  ! image2(ix,iy,it)=(image2(ix,iy,it)-imagei(ix,iy,it)**2/rtmp)
  ! 1           /rtmp/(rtmp-1.d0)
  ! end if
  ! if(ix.eq.icenter.and.iy.eq.icenter) then
  ! write(6,*) 'image before norm'
  ! 1           ,image(icenter,icenter,it),icenter,it
  ! end if
  ! if(image(ix,iy,it).lt.0.d0) then
  ! write(6,*) 'before normalizing'
  ! write(6,*) 'image lt 0',image(ix,iy,it),ix,iy,it
  ! end if
  ! if(it.eq.1.or.it.eq.nmu) then
  ! image(ix,iy,it)=image(ix,iy,it)/(rnorm)*2.d0
  ! else
  ! image(ix,iy,it)=image(ix,iy,it)/(rnorm)
  ! end if
  ! fimage=fimage+(image(ix,iy,it))/dble(nmu)
  ! if(ix.eq.icenter.and.iy.eq.icenter) then
  ! write(6,*) 'image(icenter,icenter,it),icenter,it'
  ! 1           ,image(icenter,icenter,it),icenter,it
  ! end if
  ! if(image(ix,iy,it).lt.0.d0) then
  ! write(6,*) 'image lt 0',image(ix,iy,it),ix,iy,it
  ! end if
  ! end do
  ! end do
  ! end do


  return

end subroutine diskflux
