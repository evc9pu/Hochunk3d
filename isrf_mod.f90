module isrf_mod

  implicit none
  save

  integer,parameter :: nisrf = 200

  real*8 :: pdf_i(nisrf),cdf_i(nisrf),freq_i(nisrf)
  real*8 :: logfreq_i(nisrf),logcdf_i(nisrf),nupdf_i(nisrf)
  real*8 :: isrfkap_int

end module isrf_mod
