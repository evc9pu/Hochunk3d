subroutine radface(ux,uy,uz,x1,y1,z1,R2_1,R2_2)

  use closest_wall
  !
  ! Computes the shortest distance (t) to constant radius surface from set of
  ! radii=(R1,R2) from the point (x1,y1,z1) along the line defined by the
  ! directional cosines (ux,uy,uz).  Intersections with negative "t" are
  ! NOT selected (i.e., are OPPOSITE of desired direction of propagation).
  !
  ! NOTE:  spherical surfaces are assumed to be centered at origin --
  ! x**2 + y**2 + z**2 = R**2
  !
  ! Output values:
  ! iface = which constant radius surface is closest
  ! = 0 for R1, 1 for R2, and -1 for NEITHER
  ! ifound = flag for non-negative
  ! = .T. for positive "t" found, .F. NOT.
  ! t = smallest distance for surface intersection
  ! (negative value indicated desired interesction does NOT exist...
  ! but user should test for this case using "ifound".
  !
  ! Input values:
  ! ux,uy,uz = directional cosines along x,y,z directions
  ! x1,y1,z1 = point from which distance is calculated
  ! R1,R2 = radii of two constant surfaces
  !
  ! m.j. wolff/b.a. whitney, 2000/02/04
  ! history:
  ! 00/02/28 (mjw):  disallow "0" as a valid distance
  ! 00/03/16 (mjw):  now passing R**2 instead of R
  ! using double precision (real*8)
  !

  implicit none

  ! ... scalar arguments
  real*8 :: ux,uy,uz,x1,y1,z1,R2_1,R2_2

  ! ... local scalars
  real*8 :: bb,cc,descr

  ! ... IND array indcates which "face" is associated with which roots

  bb = 2.d0*(x1*ux + y1*uy + uz*z1)

  cc = x1*x1 + y1*y1 + z1*z1 - R2_1
  descr = bb*bb - 4.d0*cc

  if (descr.gt.0.d0) then
     descr = sqrt(descr)
     call insert_t((-bb + descr)*0.5d0,1)
     call insert_t((-bb - descr)*0.5d0,1)
  end if

  cc = x1*x1 + y1*y1 + z1*z1 - R2_2
  descr = bb*bb - 4.d0*cc

  if (descr.gt.0.d0) then
     descr = sqrt(descr)
     call insert_t((-bb + descr)*0.5d0,2)
     call insert_t((-bb - descr)*0.5d0,2)
  end if

  return
end subroutine radface




