SUBROUTINE locate3(xx,n1,n2,n3,x,j)
  implicit none
  ! from numerical recipes
  ! searches an ordered table, using bisection
  INTEGER j,n1,n2,n3
  REAL*8 x,xx(n1,n2)
  INTEGER jl,jm,ju
  jl=0
  ju=n2+1
  do
     if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n3,n2).gt.xx(n3,1)).eqv.(x.gt.xx(n3,jm)))then
           jl=jm
        else
           ju=jm
        end if
     else
        exit
     end if
  end do
  j=jl
  return
END SUBROUTINE locate3
