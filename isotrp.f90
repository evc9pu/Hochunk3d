subroutine isotrp(sux,suy,suz,ux,uy,uz)

  use random
  implicit none

  real*8 :: x,y,z,r2,r,mur,sux,suy,suz,ux,uy,uz

  r2=2.d0
  mur=0.d0

  do while(r2*r2.gt.abs(mur))

     x=ran()
     x=x+x-1.d0

     y=ran()
     y=y+y-1.d0

     z=ran()
     z=z+z-1.d0

     mur=x*sux+y*suy+z*suz

     r2=x*x+y*y+z*z

  end do

  r=sqrt(r2)
  if (mur.lt.0.d0) r=-r  !reverse the direction of inward going photons

  ux=x/r
  uy=y/r
  uz=z/r

  return
end subroutine isotrp
