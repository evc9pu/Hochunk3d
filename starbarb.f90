! **********************************************************************

subroutine atmosinit(atmosfile)

  ! ----------------------------------------------------------------------
  ! Reads in a two-column array containing wavelengths and fluxes
  ! reorder frequency array, linearly interpolate cdf
  ! ======================================================================

  use dust_mod
  use atmos_interp
  use constants
  use tts_mod
  implicit none

  character(len=80) :: atmosfile
  integer :: i,id,ioerr
  real*8 :: ftmp(nhnumax),chif(nhnumax)
  real*8 :: winvtmp(nhnumax),fbol,wtmp,tmp2
  real*8,external :: Bnu
  ! -----------------------------------------------------------------------

  open(21,file=atmosfile,status='old')
  read(21,*)

  ! Read in stellar photosphere file

  do i=1,nhnumax

     read(21,*,iostat=ioerr) wtmp,ftmp(i)
     if(ioerr.ne.0) exit

     winvtmp(i)=1.d0/wtmp
     atmoswave(i)=wtmp

     ! still need work on this.  Bnu and Fnu should have similar value 
     ! in some parts of spectrum, longwave anyway.

     if (iplanckst.eq.1) then 
        tmp2=Bnu(1.2398d0/wtmp,tstar/11605.d0)
        if (wtmp.gt.1.and.wtmp.lt.1.1) then
           print*,'ftmp(i),bnu',ftmp(i),tmp2
        endif
        ftmp(i)=tmp2
     endif

  end do

  if(ioerr.eq.0) then
     write(*,*) 'atmosinit: too many frequencies in model atmosphere'
     stop
  end if

  close(21)

  do id=1,ndg

     nhnu=i-1
     write(*,*) 'nhnu=',nhnu

     ! reorder arrays with increasing frequency
     do i=1,nhnu
        winv(nhnu-i+1)=winvtmp(i)
        f(nhnu-i+1)=ftmp(i)
     end do

     ! integrate pdf
     do i=1,nhnu
        atmosnu(i)=1.2398d0*winv(i)
        lognu(i)=log10(atmosnu(i))
        call opacset(atmosnu(i))
        chif(i)=kappa(id)*f(i)

        if(i.eq.1) then
           hnuint(i)=0.d0
        else
           hnuint(i)=hnuint(i-1)+0.5d0*(f(i)+f(i-1))* &
                & (winv(i)-winv(i-1))
	     !   test.  wavelength instead of frequency
   !        hnuint(i)=hnuint(i-1)+0.5d0*(ftmp(i)+ftmp(i-1))* &
   !            & (atmoswave(i)-atmoswave(i-1))
        end if
        ! print*,winv(i),atmosnu(i),chif(i)
     end do

     fbol=hnuint(nhnu)

     if (fbol.eq.0.d0) then
        print*,''
        print*,'***************'
        print*,'error, atmosphere file is all zeros. use another.'
        print*,'this is a problem with the downloaded kurucz'
        print*,'models.   program stopping'
        stop
     end if

     kappaf(id)=0.d0
     do i=1,nhnu
        hnuint(i)=hnuint(i)/fbol
        if (i.gt.1) loghnu(i)=log10(hnuint(i))
        if (i.eq.1) then
           kappaf(id)=0.d0
        else
           kappaf(id)=kappaf(id)+0.5d0*(chif(i)+chif(i-1))* &
                & (winv(i)-winv(i-1))
        end if
        ! print*,winv(i),chif(i),kappaf
     end do
     kappaf(id)=kappaf(id)/fbol

     fbol=fbol*2.997925d14
     write(*,*) 'teff,kappaf=',(4.d0*pi*fbol/sigt)**0.25d0, &
          & kappaf(id)

  end do

  return
end subroutine atmosinit

real*8 function atmosfreq()

  ! baw, july 2002

  use random
  use atmos_interp
  implicit none

  integer :: i
  real*8 :: xi

  real*8 :: logatm

  xi=ran()

  call locate2(hnuint,nhnumax,nhnu,dble(xi),i)

  if (i.gt.1) then
  ! log interpolation
     logatm=lognu(i)+(lognu(i+1)-lognu(i)) &
        & *(log10(xi)-loghnu(i))/(loghnu(i+1)-loghnu(i))
     atmosfreq=10.d0**logatm
  ! linear interpolation
  !   atmosfreq=atmosnu(i)+(atmosnu(i+1)-atmosnu(i)) &
  !       & *(xi-hnuint(i))/(hnuint(i+1)-hnuint(i)) 
   ! wavelength instead of frequency (didn't work)
 !  atmosfreq=atmoswave(i)+(atmoswave(i+1)-atmoswave(i)) &
 !     & * (xi-hnuint(i))/(hnuint(i+1)-hnuint(i))
 !  atmosfreq=1.2398d0/atmosfreq
  else
     ! hnuint(1)=0. so log won't work
     atmosfreq=atmosnu(i)+(atmosnu(i+1)-atmosnu(i)) &
         & *(xi-hnuint(i))/(hnuint(i+1)-hnuint(i))
   end if

  return
end function atmosfreq
