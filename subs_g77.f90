module trigonometry

  implicit none
  save

contains

  real*8 function cosd(angle)
    use constants
    implicit none
    real*8,intent(in) :: angle
    cosd = cos(angle * deg2rad)
  end function cosd

  real*8 function sind(angle)
    use constants
    implicit none
    real*8,intent(in) :: angle
    sind = sin(angle * deg2rad)
  end function sind

  real*8 function tand(angle)
    use constants
    implicit none
    real*8,intent(in) :: angle
    tand = tan(angle * deg2rad)
  end function tand

end module trigonometry
