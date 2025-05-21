module filt_mod

  implicit none
  save

  character(len=50) :: filter_dir

  integer,parameter :: nwfiltmax=1200
  integer,parameter :: nfilt=21

  real*8 :: filtwave(nwfiltmax,nfilt)
  real*8 :: filtphi(nwfiltmax,nfilt)

 ! old version with SIRTF and MSX and not Herschel
 !integer,parameter :: nwfilt(nfilt) =  (/25,25,24,108,59,76,506,428,371,422,129,112,401,19,33,22,17,632,319,403,1097,7,7/)
 !new version with Herschel, no SIRTF or MSX
 integer,parameter :: nwfilt(nfilt) =  (/25,25,24,108,59,76,506,428,371,422,129,112,401,69,98,123,125,179,135,7,7/)

  character(len=6),parameter :: filtnames(nfilt) = (/'VV','RR','II',&
                                                    '2J','2H','2K',&
                                                    'I1','I2','I3','I4',&
                                                    'M1','M2','M3',&
                                                    'H1','H2','H3','H4','H5','H6',&
                                                    'Z1','Z2'/)

                                               !     'S1','S2','S3','S4','XA',&
                                               !     'XC','XD','XE','Z1','Z2'/)

end module filt_mod
