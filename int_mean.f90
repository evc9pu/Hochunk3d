subroutine int_mean(nphot)

  ! BW & KW 2007.

  ! calculate mean intensity.

  use grid_mod
  use tts_mod
  use output_mod
  use draine_mod
  use isrf_mod
  use constants, only : pi

  implicit none

  integer :: i,j,k,nphot,im
  real*8 :: jave,count,ljrat_ave,jratio,vtot,jmeanave,jr,jr1,jr2

  real*8,parameter :: lsol = 3.85d33

  call output_grid('jmeanbef',jmean)
  call output_grid('volume',vcell)

  print*,'ltot',ltot/lsol

  jmeanave=0.d0
  do i=1,nrg-1
     do j=1,ntg-1
        do k=1,npg-1
           jmean(i,j,k)=jmean(i,j,k)/vcell(i,j,k)/dble(nphot)*ltot
           jmeanave=jmeanave+jmean(i,j,k)
        end do
     end do
  end do
  print*,'jmeanave',jmeanave

  call output_grid('jmean',jmean)

  ! if (ialign.eq.1) then
  ! do i=1,nrg-1
  ! do j=1,ntg-1
  ! do k=1,npg-1
  ! do l=1,nl
  ! jml(i,j,k,l)=jml(i,j,k,l)/vcell(i,j,k)*lsun/nphot
  ! +        *ltot/pi/4.d0
  ! jmlx(i,j,k,l)=jmlx(i,j,k,l)/vcell(i,j,k)*lsun/nphot
  ! +        *ltot/pi/4.d0
  ! jmly(i,j,k,l)=jmly(i,j,k,l)/vcell(i,j,k)*lsun/nphot
  ! +        *ltot/pi/4.d0
  ! jmlz(i,j,k,l)=jmlz(i,j,k,l)/vcell(i,j,k)*lsun/nphot
  ! +        *ltot/pi/4.d0
  ! end do
  ! end do
  ! end do
  ! end do
  ! end if

  open(unit=15,file='jmean_radave.dat',status='unknown')
  write(15,*) &
       & 'i,r(rstar),r(au),  jmean (rad ave) log(jmean/0.02)'
  do i=1,nrg-1
     count=0.d0
     vtot=0.d0
     jave=0.d0
     do j=1,ntg-1
        do k=1,npg-1
			imfract(i,j,k)=0.d0
!			print*,'i,j,k',i,j,k
           ! jratio(i,j,k)=0. ! zero out for main iteration loop
           jratio=jmean(i,j,k)/isrfkap_int ! Compare with average ISRF weighted by small grain opacity
 !          jratio=jmean(i,j,k)/2.d-2 ! Compare with average ISRF
           call locate(uratio_arr,nemit,jratio,im)
	           ! if (i.eq.5)
           ! $            print*,'im,log10(jratio)',im,log10(jratio)
           if (im.gt.nemit) then
!			print*,'im, jratio',im,jratio
              print*,'jratio higher than 1.e7, ',log10(jratio)
              im=nemit
           end if
           if (im.lt.1) then
              im=1
			! print*,'im,jratio',im,jratio
           end if
		if (jratio.eq.0.) then
!       print*,'warning, jratio should not be zero'
!       print*,'i,j,k',i,j,k
			imfract(i,j,k)=0.d0
		else
!			print*,'im,jratio',im,jratio
		   if (im.eq.nemit) then
			 imfract(i,j,k)=1.d0
		   else
 !			print*,uratio_arr(im),uratio_arr(im+1),jratio
            jr1=log10(uratio_arr(im))
             jr2=log10(uratio_arr(im+1))
             jr=log10(jratio)
!		     print*,'jr,jr1,jr2',jr,jr1,jr2
             imfract(i,j,k)=(jr-jr1)/(jr2-jr1)
           endif
         endif
           imean(i,j,k)=im
           vtot=vtot+vcell(i,j,k)
           jave=jratio*vcell(i,j,k)+jave
           count=count+1.d0
!		   print*,'in int_mean,imfract',imfract(i,j,k)
        end do
     end do
     jave=jave/vtot
	 ljrat_ave=0.d0
     if (jave.gt.0.d0) ljrat_ave=log10(jave)
     write(15,900) i,ravearr(i),ravearr(i)/autors &
          & ,jave,ljrat_ave
  end do
  close(15)
900 format(i10,1x,f12.0,1x,f13.5,1x,e12.5,1x,f10.5)

  return
end subroutine int_mean
