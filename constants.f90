module constants

  implicit none

  real*8,parameter :: pi = 3.1415926535897932384626433d0
  real*8,parameter :: pihalf = pi * 0.5d0
  real*8,parameter :: r2p = 2.d0 * pi
  real*8,parameter :: deg2rad = pi / 180.d0
  real*8,parameter :: rad2deg = 180.d0 / pi

  real*8,parameter :: rsun = 6.9598d10
  real*8,parameter :: rsol = rsun
  real*8,parameter :: msun = 1.989d33
  real*8,parameter :: msol = msun
  real*8,parameter :: lsun = 3.839d33
  real*8,parameter :: lsol = lsun


  real*8,parameter :: gn = 6.67428d-8     ! Gravitational constant [cgs]
  real*8,parameter :: sigt = 5.670400d-5  ! Stefan-Boltzmann constant [cgs]
  real*8,parameter :: k = 1.3806503d-16   ! Plan constant [cgs]
  real*8,parameter :: mH = 1.67262158d-24 ! Hydrogen mass [cgs]
  real*8,parameter :: c = 2.99792458d10   ! Speed of light [cgs]
  real*8,parameter :: h = 6.6260755d-27   ! Boltzmann constant [cgs]
  real*8,parameter :: au2cm = 1.49598d13  ! Astronomical unit [cgs]
  real*8,parameter :: pc2cm = 3.0857d18   ! parsec (cgs)
  
  ! real,parameter :: hertzfact = 2.41805e14  ! convert frequency from jon's units to hertz

end module constants

