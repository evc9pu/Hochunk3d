subroutine initp_sample()

  use tts_mod
  use stokesvar_mod
  use random
  use constants

  implicit none

  real*8 :: xran

  ! sample angle of exit using appropriate limb darkening law.
  if(limb.eq.0) then
     ! isotropic intensity, sample from n=mu
     xran=ran()
     cost=sqrt(xran)
  else
     xran=ran()
     call darkening(xran,cost)
  end if
  sint=sqrt(1.d0-cost*cost)

  ! sample phi
  xran=ran()
  phi=r2p*xran
  cosp=cos(phi)
  sinp=sin(phi)

end subroutine initp_sample

subroutine initp_transform()

  use tts_mod
  use stokesvar_mod
  use random
  use constants

  implicit none

  real*8 :: costi,sinti,cospn,phinew

  ! transform to coordinate system of star
  if(abs(cosb).lt.0.9999999d0) then
     costi=cost*cosb-sint*sinb*cosp
     if(abs(costi).gt.1.d0) then
        write(6,*) 'initp: costi gt 1',costi
        if(costi.gt.1.d0) costi=1.d0
        if(cost.lt.-1.d0) costi=-1.d0
     end if
     sinti=sqrt(1.d0-costi**2)
     if(sinti.gt.0.0000001d0)then
        cospn=(cost-costi*cosb)/(sinti*sinb)
        if(abs(cospn).gt.1.d0) then
           if (abs(cospn).gt.1.01d0) &
                & write(6,*)'cospn gt 1 in initp',cospn
           if (cospn.lt.0.d0) then
              cospn=-1.d0
           else
              cospn=1.d0
           end if
        end if
        if(phi.lt.pi) then
           phinew=acos(cospn)
           if(phinew.gt.pi) then
              write(6,*) 'phinew wrong in initp'
              write(6,*) 'phinew,pi',phinew,pi
              phinew=pi
           end if
        else
           phinew=r2p-acos(cospn)
           if(phinew.lt.pi) then
              write(6,*) 'phinew wrong in initp'
              write(6,*) 'phinew,pi',phinew,pi
              phinew=pi
           end if
        end if
     else
        phinew=0.d0
     end if
     phi=phinew
     cost=costi
     sint=sinti
  else
     if(cosb.lt.0.d0) cost=-cost
  end if
  phi=phi+lp
  if(phi.gt.r2p) phi=phi-r2p
  if (phi.lt.0.d0) phi=phi+r2p
  if(phi.lt.0.d0) write(6,*)'phi lt 0',phi
  if(phi.gt.r2p) write(6,*) 'phi gt 2pi',phi
  sinp=sin(phi)
  cosp=cos(phi)

end subroutine initp_transform
