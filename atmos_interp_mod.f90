module atmos_interp

  implicit none
  save

  integer :: nhnu
  integer,parameter :: nhnumax=60000
  real*8,dimension(nhnumax) :: atmosnu,hnuint,lognu,loghnu,atmoswave
  real*8 :: winv(nhnumax),f(nhnumax)

end module atmos_interp
