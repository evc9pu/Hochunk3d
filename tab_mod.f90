module tab_mod

  implicit none

  integer,parameter :: narr=720

  real*8,dimension(narr) :: tharr,cosarr

  real*8,dimension(narr,4) :: &
       & s11arr,s12arr,s13arr,s14arr,&
       & s21arr,s22arr,s23arr,s24arr,&
       & s31arr,s32arr,s33arr,s34arr,&
       & s41arr,s42arr,s43arr,s44arr

  real*8 :: dct

  integer :: nelems

end module tab_mod

