module random

!history:
! 2012/03/22 (mjw):  replace code in function RAN with a call to
!			RANDOM_STANDARD_UNIFORM (include "use" here)
!

  use random_standard_uniform_mod
  use ecuyer_cote_mod

  implicit none
  save

  integer,private :: idum=-12421412
  integer :: fractal_gen=32,current_gen

  integer,private,parameter :: izmax = 145
  real*8,private :: inczeta(izmax)

contains

  subroutine set_random_seed(iseed)
    implicit none
    integer,intent(in) :: iseed
    idum = iseed
  end subroutine set_random_seed

  subroutine rantheta(st,ct,s2t,c2t)
    implicit none
    real*8 :: st,ct,s2t,c2t
    real*8 :: x,x2,y,y2,r,r2
    r2=2.d0
    do while(r2.gt.1.d0)
       x=ran()
       x=x+x-1.d0
       y=ran()
       y=y+y-1.d0
       x2=x*x
       y2=y*y
       r2=x2+y2
    end do
    r=sqrt(r2)
    st=y/r
    ct=x/r
    s2t=y2/r2
    c2t=x2/r2
  end subroutine rantheta

  real*8 function bbfreq(teff)

    implicit none

    integer :: i      
    real*8 :: teff
    real*8 :: z0,z1,z2,z3,z4
    logical,save :: first = .true.

    if (first) then
       call zetainit()
       write(*,*) 'initializing inczeta'
       first=.false.
    end if

    z0=inczeta(izmax)*ran()
    z1=ran()
    z2=ran()
    z3=ran()
    z4=ran()

    i=1
    do while(z0 .gt. inczeta(i))
       i=i+1
    end do

    bbfreq=-log(z1*z2*z3*z4)*teff/dble(i)

  end function bbfreq

  subroutine zetainit()

    implicit none

    integer i,j
    real*8 temp

    do i=1,izmax
       inczeta(i)=0.d0
    end do

    do j=izmax,1,-1
       temp=1.d0/dble(j*j*j*j)
       do i=izmax,j,-1
          inczeta(i)=inczeta(i)+temp
       end do
    end do

  end subroutine zetainit

  real*8 function ran() result(xi)
    implicit none
!    real*8,save :: am
!    integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
!    integer, save :: ix=-1,iy=-1,k
!    if (idum <= 0 .or. iy < 0) then
!       am=nearest(1.0d0,-1.0d0)/im
!       iy=ior(ieor(888889999,abs(idum)),1)
!       ix=ieor(777755555,abs(idum))
!       idum=abs(idum)+1
!    end if
!    ix=ieor(ix,ishft(ix,13))
!    ix=ieor(ix,ishft(ix,-17))
!    ix=ieor(ix,ishft(ix,5))
!    k=iy/iq
!    iy=ia*(iy-k*iq)-ir*k
!    if (iy < 0) iy=iy+im
!    xi=am*ior(iand(im,ieor(ix,iy)),1)
     xi = random_standard_uniform()
  end function ran

real*8 function ran_old() result(xi)
  implicit none
  real*8,save :: am
  integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
  integer, save :: ix=-1,iy=-1,k
  if (idum <= 0 .or. iy < 0) then
     am=nearest(1.0d0,-1.0d0)/im
     iy=ior(ieor(888889999,abs(idum)),1)
     ix=ieor(777755555,abs(idum))
     idum=abs(idum)+1
  end if
  ix=ieor(ix,ishft(ix,13))
  ix=ieor(ix,ishft(ix,-17))
  ix=ieor(ix,ishft(ix,5))
  k=iy/iq
  iy=ia*(iy-k*iq)-ir*k
  if (iy < 0) iy=iy+im
  xi=am*ior(iand(im,ieor(ix,iy)),1)
end function ran_old

  real*8 function gasdev()
    integer,save :: iset = 0
    real*8,save :: gset
    real*8 :: fac,rsq,v1,v2
    if (iset.eq.0) then
       do
          v1=2.d0*ran()-1.d0
          v2=2.d0*ran()-1.d0
          rsq=v1*v1+v2*v2
          if(rsq.lt.1.d0.and.rsq.gt.0.d0) exit
       end do
       fac=sqrt(-2.d0*log(rsq)/rsq)
       gset=v1*fac
       gasdev=v2*fac
       iset=1
    else
       gasdev=gset
       iset=0
    end if
  end function gasdev

end module random



