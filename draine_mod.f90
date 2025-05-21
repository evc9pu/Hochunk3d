module draine_mod

  implicit none
  save

  integer,parameter :: ndraine = 1000
  integer,parameter :: nemit = 18

  real*8 :: pdf_d(nemit,ndraine),cdf_d(nemit,ndraine)
  real*8 :: logcdf_d(nemit,ndraine)
  real*8 :: freq_d(ndraine),logfreq_d(ndraine)
  real*8 :: uratio_arr(nemit)

end module draine_mod
