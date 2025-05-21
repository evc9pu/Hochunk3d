subroutine rep_isrf(freq)

  ! BAW, 20070815

  use isrf_mod
  use random
  
  implicit none

  integer :: ifreq
  real*8 :: xran,logfreq,freq
  
  xran=ran()

  call locate(cdf_i,nisrf,xran,ifreq)
  ! call locate(logcdf_i,nisrf,xran,ifreq)     !if using nu_Inu

  if(ifreq.eq.nisrf)ifreq=ifreq-1
  if(ifreq.eq.0) then
     print*,'ifreq=0, emit_isrf; problem?'
     ifreq=1
  end if

  logfreq=logfreq_i(ifreq)+ &
       & (log10(xran)-logcdf_i(ifreq))/ &
       ! & ((xran)-logcdf_i(ifreq))/ &   !if using nu_Inu &
       & (logcdf_i(ifreq+1)-logcdf_i(ifreq))* &
       & (logfreq_i(ifreq+1)-logfreq_i(ifreq))
  freq=10.d0**logfreq

  ! don't do logfreq interpolation
  ! freq=(xran-cdf_i(ifreq))/
  ! $     (cdf_i(ifreq+1)-cdf_i(ifreq))*
  ! $     (freq_i(ifreq+1)-freq_i(ifreq))
  ! $     + freq_i(ifreq)
  ! this looks better when spectrum is decreasing with wavelength and
  ! worse as it increases.

  return
end subroutine rep_isrf

