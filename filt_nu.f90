subroutine filt()

  ! read in filter functions for the images

  use filt_mod
  implicit none

  integer :: i,nf,ifilt,ich,ich2
  character(len=5) :: junk
  character(len=100) :: filnam
  real*8 :: fact,w,f,fint,w_ave,f_peak,w_peak,w_half1,w_half2
  real*8 :: w_avew,fintw

  fact=2.9979d14
  ! fact=1.2398d0

  open(unit=16,file='filters.log')

  write(16,'("    w_peak    f_peak       norm       w_half1'// &
       & '    w_half2      fwhm      w_freq     w_wave")')

  ! nfilt is number of filter files
  ! loop over filter files
  do ifilt=1,nfilt
     fint=0.d0
     filnam=trim(filter_dir)//'/'//trim(filtnames(ifilt))//'.txt'
     print*,"Reading "//trim(filnam)
     ! nf=number of lines in this filter file. subtract 1 for the header
     nf=nwfilt(ifilt)-1

     ! read in filter file
     f_peak=0.d0
     open(unit=14,file=filnam,status='old')
     read(14,*) junk
     do i=1,nf
        read(14,*) w,f
        ! convert wavelength to frequency
        filtwave(i,ifilt)=fact/w
        filtphi(i,ifilt)=f
        if (f.gt.f_peak) then
           f_peak=f
           w_peak=filtwave(i,ifilt)
        end if
     end do
     f_peak=f_peak/2.d0   !for FWHM calculation later on
     close(14)

     ! integrate filter to determine normalization factor
     ! tom says they are already normalized but let's check
     fint=0.d0
     w_ave=0.d0
     w_avew=0.d0
     fintw=0.d0
     do i=1,nf-1
        fint=0.5d0*(filtphi(i+1,ifilt)+filtphi(i,ifilt))* &
             & (filtwave(i+1,ifilt)-filtwave(i,ifilt))+fint
        fintw=0.5d0*(filtphi(i+1,ifilt)+filtphi(i,ifilt))* &
             & (fact/filtwave(i,ifilt)-fact/filtwave(i+1,ifilt))+fintw
        w_ave=w_ave+(filtwave(i+1,ifilt)-filtwave(i,ifilt))* &
             & 0.5d0*(filtphi(i+1,ifilt)+filtphi(i,ifilt))* &
             & 0.5d0*(filtwave(i+1,ifilt)+filtwave(i,ifilt))
        w_avew=w_avew+(fact/filtwave(i,ifilt)- &
             & fact/filtwave(i+1,ifilt))* &
             & 0.5d0*(filtphi(i+1,ifilt)+filtphi(i,ifilt))* &
             & 0.5d0*(fact/filtwave(i+1,ifilt) &
             & +fact/filtwave(i,ifilt))
     end do
     w_ave=w_ave/fint
     w_avew=w_avew/fintw


     ! calculate FWHM
     ich=1
     ich2=1
     do i=1,nf
        ! print*,'calculating w_half1'
        if (ich.eq.1) then
           ! print*,'ich=1',fact/filtwave(i,ifilt)
           ! $    ,filtphi(i,ifilt),f_peak
           if (filtphi(i,ifilt).gt.f_peak) then
              w_half1=filtwave(i,ifilt)
              ! print*,'wave,filtphi',fact/filtwave(i,ifilt),filtphi(i,ifilt)
              ich=0
           end if
        end if
        ! print*,'filtwave,wpeak',fact/filtwave(i,ifilt),fact/w_peak
        if (filtwave(i,ifilt).gt.w_peak) then
           ! print*,'calculating w_half2'
           ! print*,fact/filtwave(i,ifilt),f_peak
           if (ich2.eq.1) then
              ! print*,'ich2 eq 1',filtphi(i,ifilt),f_peak
              if (filtphi(i,ifilt).lt.f_peak) then
                 w_half2=filtwave(i,ifilt)
                 ich2=0
              end if
           end if
        end if
     end do

     write(16,'(F10.4,2X,ES10.3,2X,ES10.3,1X,5(F10.3,1X))') &
          & fact/w_peak,f_peak,fint,fact/w_half2,fact/w_half1, &
          & fact/w_half1-fact/w_half2,fact/w_ave,w_avew

     ! rescale filter
     do i=1,nf
        filtphi(i,ifilt)=filtphi(i,ifilt)/fint
     end do
  end do

  close(unit=16)

end subroutine filt

