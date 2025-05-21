subroutine dustmat(p1,p2,p3,p4,cost,cost2,idust2)

  ! **********************************************************************
  !
  ! this program calculates the elements of the phase matrix for a
  ! simple representation of the mrn dust mixture using the algorithms
  ! for the ultraviolet region due to richard l. white ap.j. 229, 954,
  ! 1979.
  !
  ! ***********************************************************************
  !
  ! cost = cos(angle) of scattering (i.e. angle between incident
  ! photon and scattered photon)
  ! g = mean value of cosine of scattering angle (henyey-greenstein)
  ! pl = peak linear polarization
  ! pc = peak value of linear to circular conversion
  ! sc = asymmetry of the circular polarization.
  ! p1 = intensity phase function
  ! p2 = polarization function
  ! p3 = skew polarization
  ! p4 = circular polarization
  !
  ! the scattering matrix for (i,q,u,v) is of the form
  !
  ! p1    p2    0   0
  ! p2    p1    0   0
  ! 0     0     p3  -p4
  ! 0     0     p4   p3
  !
  ! **********************************************************************
  use dust_mod
  use constants
  implicit none

  real*8 :: cost2,cost,p4,p3,p2,p1
  real*8 :: temp
  integer :: idust2
  ! data g, pl, pc, sc / 0.56, 0.51, 0.39, 1.0 /
  !
  temp=g2p1(idust2)-twog(idust2)*cost
  ! H-G function
  ! p1=onemg2(idust2)/(temp*sqrt(temp))
  ! modified Cornette & Shanks function---make sure peak is known
  p1=1.5d0*onemg2(idust2)/(2.d0+g2(idust2))*(1.d0+cost2)/(temp*sqrt(temp))

  p2 = -pl(idust2)*p1*(1.0d0-cost2)/(1.0d0+cost2)
  p3 = p1*2.0d0*cost/(1.0d0+cost2)
  if(abs(cost).gt.1.d0) then
     write(6,*) 'in dustmat, cost.gt.1',cost
     if(cost.gt.1.d0) then
        cost=1.d0
     else
        cost=-1.d0
     end if
  end if
  ! circular polarization, comment out if you don't care
  ! angle in degrees!
  ! phi=acosd(cost)
  ! f=3.13*phi*exp(-7.0*phi/180.d0)
  ! now convert to radians
  ! f2=(phi+sc*f)*deg2rad
  ! fphi= (1+3.13*sc*exp(-7.0*phi/pi))*phi
  ! c=(cos(f2))**2
  ! p4 = -pc*p1*(1-c)/(1+c)
  p4=0.d0

  return
end subroutine dustmat
