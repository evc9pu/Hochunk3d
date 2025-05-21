subroutine rep_draine(jmean,freq,jmf)

  ! 2007, KW. modified for my code, BAW, 20070815

  use draine_mod
  use random

  implicit none

  integer :: ifreq,jmean,jmeanp1
  real*8 :: xran,logfreq,freq,jmf,logfreq1,logfreq2

  ! need to find out if PAHs are destroyed above a certain energy density

  xran=ran()

  ! call locate(cdf_d,ndraine,xran,ifreq)
  call locate3(cdf_d,nemit,ndraine,jmean,xran,ifreq)

  if(ifreq.eq.ndraine)ifreq=ifreq-1

  if(ifreq.eq.0) then
     print*,'ifreq=0, rep_draine; problem?'
     ifreq=1
  end if

  logfreq1=(log10(xran)-logcdf_d(jmean,ifreq))/ &
       & (logcdf_d(jmean,ifreq+1)-logcdf_d(jmean,ifreq))* &
       & (logfreq_d(ifreq+1)-logfreq_d(ifreq)) &
       & + logfreq_d(ifreq)

!repeat this for jmean+1
  jmeanp1=jmean+1
 if (jmeanp1.le.nemit) then
   call locate3(cdf_d,nemit,ndraine,jmeanp1,xran,ifreq)
 
    if(ifreq.eq.ndraine)ifreq=ifreq-1
 
    if(ifreq.eq.0) then
      print*,'ifreq=0, rep_draine; problem?'
      ifreq=1
    end if
 
   logfreq2=(log10(xran)-logcdf_d(jmeanp1,ifreq))/ &
      & (logcdf_d(jmeanp1,ifreq+1)-logcdf_d(jmeanp1,ifreq))* &
      & (logfreq_d(ifreq+1)-logfreq_d(ifreq)) &
      & + logfreq_d(ifreq)

   logfreq = logfreq1 + jmf * (logfreq2 - logfreq1)
   
!	print*,'logfreq1,logfreq2,jmf,logfreq',logfreq1,logfreq2,jmf,logfreq
   else
     
     logfreq=logfreq1
   
 endif


  freq=10.d0**logfreq

  ! lamran=10.d0**loglam_c(ifreq)

  return
end subroutine rep_draine

