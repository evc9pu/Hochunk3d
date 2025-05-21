subroutine draine_cdf

  ! KW, 2007, adapted for my code, BAW, 20070815
  ! 20070817, BAW, use j_nu instead of nu*Pnu/dnu, use freq instead of lam

  use tts_mod, only : par_dir
  use draine_mod
  implicit none

  real*8 :: lam,pdftmp

  integer :: i,j

  character(len=25) :: filename(18)

  filename(1)='emit_draine_0.50.dat'
  filename(2)='emit_draine_1.00.dat'
  filename(3)='emit_draine_2.00.dat'
  filename(4)='emit_draine_5.00.dat'
  filename(5)='emit_draine_10.0.dat'
  filename(6)='emit_draine_20.0.dat'
  filename(7)='emit_draine_30.0.dat'
  filename(8)='emit_draine_1e2.dat'
  filename(9)='emit_draine_3e2.dat'
  filename(10)='emit_draine_1e3.dat'
  filename(11)='emit_draine_3e3.dat'
  filename(12)='emit_draine_1e4.dat'
  filename(13)='emit_draine_3e4.dat'
  filename(14)='emit_draine_1e5.dat'
  filename(15)='emit_draine_3e5.dat'
  filename(16)='emit_draine_1e6.dat'
  filename(17)='emit_draine_3e6.dat'
  filename(18)='emit_draine_1e7.dat'

  uratio_arr = (/0.5,1.0,2.0,5.0,10.0,20.0,30.0,1.e2,3.e2,1.e3,&
       & 3.e3,1.e4,3.e4,1.e5,3.e5,1.e6,3.e6,1.e7/)

  do j=1,nemit

     open(10,file=trim(par_dir)//'/'//filename(j),status='old')
     do i=1,ndraine

        read(10,*) lam,pdftmp
        pdf_d(j,ndraine-i+1)=pdftmp

        ! weird j&k units
        freq_d(ndraine-i+1)=1.2398d0/lam

        ! pdf_d(j,ndraine-i+1)=pdf_d(j,ndraine-i+1)*freq_d(ndraine-i+1)
        ! arrays are resorted so pdf will increase with increasing freq.

     end do
     close(10)

  end do

  ! Make CDF and log_10(CDF) for emissivity.  Also take logs of frequency

  do j=1,nemit

     cdf_d(j,1)=0.0d0
     do i=2,ndraine
        cdf_d(j,i)=cdf_d(j,i-1)+0.5d0*(pdf_d(j,i)+pdf_d(j,i-1))*(freq_d(i)-freq_d(i-1))
     end do

     cdf_d(j,1)=0.d0
     logcdf_d(j,1)=-10.d0
     logfreq_d(1)=0.d0
     do i=2,ndraine
        cdf_d(j,i)=cdf_d(j,i)/cdf_d(j,ndraine) ! normalize
        logcdf_d(j,i)=log10(cdf_d(j,i))
        logfreq_d(i)=log10(freq_d(i))
     end do

  end do

  return

end subroutine draine_cdf
