subroutine phiface(A1,B1,C1,D1,A2,B2,C2,D2,ux,uy,uz,x1,y1,z1)

  use closest_wall

  !
  ! Computes the shortest distance (t) to constant phi (azimuthal angle)
  ! surface set of (phi1,phi2) from the point (x1,y1,z1)
  ! along the line defined by the directional cosines (ux,uy,uz).
  ! Intersections with negative "t" are NOT selected (i.e., are
  ! OPPOSITE of desired direction of propagation).
  !
  ! NOTE:  constant phi surfaces are planes that contain the z-axis.
  ! Ax + By +Cz + D = 0, where the coefficients are
  ! pre-computed.
  !
  ! Output values:
  ! iface = which constant radius surface is closest
  ! = 0 for phi1, 1 for phi2, and -1 for NEITHER
  ! ifound = flag for non-negative
  ! = .T. for positive "t" found, .F. NOT.
  ! t = smallest distance for surface intersection
  ! (negative value indicated desired interesction does NOT exist...
  ! but user should test for this case using "ifound".
  !
  ! Input values:
  ! ux,uy,uz = directional cosines along x,y,z directions
  ! x1,y1,z1 = point from which distance is calculated
  ! A1...D1,A2...D2 = plane equation coefficients for two constant surfaces
  !
  ! m.j. wolff/b.a. whitney, 2000/02/04
  ! history:
  !
  ! 00/03/16 (mjw):  Using double precision (real*8).
  ! Redefine array IND as integer :: (was real).
  ! Identify only roots greater than 0.
  !
  !
  !

  implicit none

  ! ... scalar arguments
  real*8 :: A1,A2,B1,B2,C1,C2,D1,D2,ux,uy,uz,x1,y1,z1

  ! ... local scalars
  real*8 :: denom1,denom2

  denom1 = (A1*ux + B1*uy + C1*uz)
  if (denom1.ne.0.d0) then
     call insert_t(-(A1*x1 + B1*y1 + C1*z1 + D1) / denom1,5)
  end if

  denom2 = (A2*ux + B2*uy + C2*uz)
  if (denom2.ne.0.d0) then
     call insert_t(-(A2*x1 + B2*y1 + C2*z1 + D2) / denom2,6)
  end if

  return
end subroutine phiface




