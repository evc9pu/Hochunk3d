subroutine Rdist(ifound,t,ux,uy,uz,x1,y1,z1,R)
  !
  ! Computes the shortest distance (t) to constant radius surface R
  ! from the point (x1,y1,z1) along the line defined by the
  ! directional cosines (ux,uy,uz).  Intersections with negative "t" are
  ! NOT selected (i.e., are OPPOSITE of desired direction of propagation).
  !
  ! NOTE:  spherical surfaces are assumed to be centered at origin --
  ! x**2 + y**2 + z**2 = R**2
  !
  ! Output values:
  ! ifound = flag for non-negative
  ! = .T. for positive "t" found, .F. NOT.
  ! t = smallest distance for surface intersection
  ! (negative value indicated desired interesction does NOT exist...
  ! but user should test for this case using "ifound".
  !
  ! Input values:
  ! ux,uy,uz = directional cosines along x,y,z directions
  ! x1,y1,z1 = point from which distance is calculated
  ! R = radius of constant surface
  !
  ! m.j. wolff/b.a. whitney, 2000/02/04
  ! history:
  !
  !

  implicit none

  ! ... scalar arguments
  real*8 :: ux,uy,uz,x1,y1,z1,t,R
  logical :: ifound

  ! ... local arrays
  real*8 :: root(2)

  ! ... local scalars
  real*8 :: bb,cc,descr,rootmax,rootmin


  bb = 2.d0*(x1*ux + y1*uy + uz*z1)
  ! simplify since aa = 1.
  cc = x1*x1 + y1*y1 + z1*z1 - R*R
  descr = bb*bb - 4.d0*cc

  ! print*,'descr = ',descr,bb,cc
  if (descr.lt.0.d0) then
     root(1) = -999.d0
     root(2) = -999.d0
     ! print*,'shit...no intersection'
  else
     descr = sqrt(descr)
     root(1) = (-bb + descr) / 2.d0
     root(2) = (-bb - descr) / 2.d0
  end if

  rootmax = max(root(1),root(2))
  rootmin = min(root(1),root(2))

  ! both roots .lt. 0
  if (rootmax.lt.0.d0) then
     ifound = .false.
     t = -1.d0
     ! one root .lt. 0
  else if (rootmin.lt.0.d0) then
     ifound = .true.
     t = rootmax
     ! both roots ge 0
  else
     ifound = .true.
     t = rootmin
  end if

  return
end subroutine Rdist




