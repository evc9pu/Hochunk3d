subroutine envset

  ! 2001/02/17 (baw) constants for envelope assuming TSC
  ! sets cavity shape parameters also.
  ! 2003/12/26 set c1e,c2e when hole opening angle=90.

  use tts_mod
  use opacin_mod
  use constants

  implicit none

  real*8 :: const,sphmass,rd,rdcgs,rmincgs,z1max,r1max,rtmp,r2max

  ! infall;  Rd is disk radius, units of rmin
  ! BAW 7/5/98 make rd = rchole = rc
  ! rd = rchole
  ! BAW 12/15/99 make rd=rc
  rd=rc
  rdcgs=rd*rstar*rsol

  ! rhoe0 is factor in front of density distribution, given by
  ! infall calculation.
  ! rhoe0=(infallrate)/(4*pi)/sqrt(G*Mcore)/rd**1.5

  const=1.d0/3.15576d7*sqrt(msol)/(4.d0*pi)/sqrt(gn)
  write(6,*)'const',const

  rhoe0=const*rate/sqrt(massc)/rdcgs**1.5d0
  sphmass=2.d0/3.d0/3.15576d7/sqrt(msol)*rate/sqrt(massc)/sqrt(2.d0*gn)
  sphmass=sphmass*(rmax*rstar*rsol)**1.5d0
  write(6,*) 'rhoe0 of envelope',rhoe0

  rmincgs=rstar*rsol
  write(6,*) 'rmincgs',rmincgs

  windmu0=cos(thetmu0*deg2rad)

  ! ambient density set in input file now, 10/21/02
  ! rhoamb=3.34e-20

  ! hole in bubble
  ! roa=rmax*tand(buboa)
  cosbuboa=cos(buboa*deg2rad)

  ! 1999
  ! outflow boundaries
  ! z1max is height where opening angle is measured.  take as outer
  ! bound.
  z1max=rmax
  if (thet1.lt.89.999d0) then
     r1max=z1max*tan(thet1*deg2rad)
     ! z=a+b*x**beta
     ! c1=b
     ! z01=a
     ! ex1=beta
     ! r=x  (r is cylindrical radius)
     c1e=(z1max-z01)/r1max**ex1
     print*,'z1max,z01,r1max,ex1,c1e'
     print*,z1max,z01,r1max,ex1,c1e
     ! z01 is input
     ! check
     if(z01.lt.0.d0) then
        rtmp=(-z01/c1e)**(1.d0/ex1)
        print*,'hole 1 intersects disk at ',rtmp/autors, ' AU'
     end if
     if(ipoly.eq.1) then
        r2max=z1max*tan(thet2*deg2rad)
        c2e=(z1max-z01)/r2max**ex2   !bug??????  2011 nov 2
        c2e=(z1max-z02)/r2max**ex2
        print*,'r2max,ex2,c2e',r2max,ex2,c2e
        ! z02 is input
        if(z02.lt.0.d0) then
           rtmp=(-z02/c2e)**(1.d0/ex2)
           print*,'hole 2 intersects disk at ',rtmp/autors, ' AU'
        end if
     end if
  else
     c1e=0.d0
     c2e=0.d0
  end if

  return

end subroutine envset
