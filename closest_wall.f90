module closest_wall

  implicit none
  save

  real*8 :: tn1,tp1,tp2
  integer :: in1,ip1,ip2

  logical :: debug = .false.

contains

  subroutine reset_t()
    implicit none
    tn1 = -huge(tn1)
    tp1 = +huge(tp1)
    tp2 = +huge(tp2)
    in1 = 0
    ip1 = 0
    ip2 = 0
    if(debug) print *,'-----'
  end subroutine reset_t

  subroutine insert_t(t,i)
    implicit none
    real*8,intent(in)    :: t
    integer,intent(in)    :: i
    if(debug) print *,'inserting ',t,i
    if(t < 0.d0) then
       if(t > tn1) then
          tn1 = t
          in1 = i
       end if
    else
       if(t < tp1) then
          tp2 = tp1
          ip2 = ip1
          tp1 = t
          ip1 = i
       else if(t < tp2) then
          tp2 = t
          ip2 = i
       end if
    end if
  end subroutine insert_t

  subroutine find_next_wall(on_wall,t,i)
    implicit none
    logical(1),intent(in)  :: on_wall
    real*8,intent(out)    :: t
    integer,intent(out)    :: i
    if(on_wall) then
       if(abs(tn1) < abs(tp1)) then
          t = tp1
          i = ip1
       else
          t = tp2
          i = ip2
       end if
    else
       t = tp1
       i = ip1
    end if
    if(debug) print *,'choosing ',t,i
    if(t==huge(t)) stop "ERROR - no valid solution"
  end subroutine find_next_wall

end module closest_wall
