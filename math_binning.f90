module math_binning

  implicit none
  save

  private
  public :: locate
  public :: ipos
  public :: xval

  interface locate
     module procedure locate_sp
     module procedure locate_dp
  end interface

  interface ipos
     module procedure ipos_sp
     module procedure ipos_dp
  end interface

  interface xval
     module procedure xval_sp
     module procedure xval_dp
  end interface

  integer,parameter :: sp = selected_real_kind(p=6,r=37)
  integer,parameter :: dp = selected_real_kind(p=15,r=307)

contains


  integer function locate_dp(xx,x)

    !**********************************************************************!
    ! This function is used to locate the position in an array
    !
    ! Versions:
    !
    ! Original code : Numerical Recipes
    ! April 2007    : TR - changed variable types
    !                      added out of bounds clause
    !**********************************************************************!

    implicit none
    real(dp), dimension(:), intent(in) :: xx
    real(dp), intent(in) :: x
    integer :: n,jl,jm,ju
    logical :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do

    if (x == xx(1)) then
       locate_dp = 1
    else if (x == xx(n)) then
       locate_dp = n-1
    else if(x > xx(n) .or. x < xx(1)) then
       locate_dp = -1
    else
       locate_dp = jl
    end if

  end function locate_dp

  integer pure function ipos_dp(xmin,xmax,x,nbin)

    !**********************************************************************!
    ! ipos : find bin corresponding to value x
    !**********************************************************************!

    implicit none

    real(dp),intent(in) :: xmin,xmax
    ! range of values

    real(dp),intent(in) :: x
    ! the value to bin

    integer,intent(in) :: nbin
    ! number of bins

    real(dp) :: frac

    if(xmax > xmin) then

    if(x < xmin) then
       ipos_dp = 0
    else if(x > xmax) then
       ipos_dp = nbin+1
    else
       frac=(x-xmin)/(xmax-xmin)
       ipos_dp=int(frac*real(nbin))+1
    end if

    else

      if(x < xmax) then
         ipos_dp = 0
      else if(x > xmin) then
         ipos_dp = nbin+1
      else
         frac=(x-xmin)/(xmax-xmin)
         ipos_dp=int(frac*real(nbin))+1
      end if

    end if

  end function ipos_dp

  real(dp) pure function xval_dp(xmin,xmax,i,nbin)

    !**********************************************************************!
    ! Find central value of bin i
    !**********************************************************************!

    implicit none

    real(dp),intent(in) :: xmin,xmax
    ! range of values

    integer,intent(in) :: i
    ! the bin number

    integer,intent(in) :: nbin
    ! number of bins

    real(dp) :: frac

    frac=(real(i-1)+0.5d0)/dble(nbin)

    xval_dp=frac*(xmax-xmin)+xmin

  end function xval_dp


  integer function locate_sp(xx,x)

    !**********************************************************************!
    ! This function is used to locate the position in an array
    !
    ! Versions:
    !
    ! Original code : Numerical Recipes
    ! April 2007    : TR - changed variable types
    !                      added out of bounds clause
    !**********************************************************************!

    implicit none
    real(sp), dimension(:), intent(in) :: xx
    real(sp), intent(in) :: x
    integer :: n,jl,jm,ju
    logical :: ascnd
    n=size(xx)
    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do

    if (x == xx(1)) then
       locate_sp = 1
    else if (x == xx(n)) then
       locate_sp = n-1
    else if(x > xx(n) .or. x < xx(1)) then
       locate_sp = -1
    else
       locate_sp = jl
    end if

  end function locate_sp

  integer pure function ipos_sp(xmin,xmax,x,nbin)

    !**********************************************************************!
    ! ipos : find bin corresponding to value x
    !**********************************************************************!

    implicit none

    real(sp),intent(in) :: xmin,xmax
    ! range of values

    real(sp),intent(in) :: x
    ! the value to bin

    integer,intent(in) :: nbin
    ! number of bins

    real(sp) :: frac

    if(xmax > xmin) then

    if(x < xmin) then
       ipos_sp = 0
    else if(x > xmax) then
       ipos_sp = nbin+1
    else
       frac=(x-xmin)/(xmax-xmin)
       ipos_sp=int(frac*real(nbin))+1
    end if

    else

      if(x < xmax) then
         ipos_sp = 0
      else if(x > xmin) then
         ipos_sp = nbin+1
      else
         frac=(x-xmin)/(xmax-xmin)
         ipos_sp=int(frac*real(nbin))+1
      end if

    end if

  end function ipos_sp

  real(sp) pure function xval_sp(xmin,xmax,i,nbin)

    !**********************************************************************!
    ! Find central value of bin i
    !**********************************************************************!

    implicit none

    real(sp),intent(in) :: xmin,xmax
    ! range of values

    integer,intent(in) :: i
    ! the bin number

    integer,intent(in) :: nbin
    ! number of bins

    real(sp) :: frac

    frac=(real(i-1)+0.5)/real(nbin)

    xval_sp=frac*(xmax-xmin)+xmin

  end function xval_sp


end module math_binning
