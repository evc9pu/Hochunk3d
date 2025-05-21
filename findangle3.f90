subroutine findangle3(cost1,sint1,phi1,cosp1,sinp1 &
     & ,cost2,sint2,phi2,cosp2,sinp2,costnew,sintnew,phinew,ri1)

  ! HISTORY
  ! 9/7/99  BAW, add i1 angle for stokes scattering diagram.

  ! input:
  ! cost1,sint1,phi1    --incident photon direction (t1,phi1)
  ! (or normal to surface)
  ! cost2,sint2,phi2    --scattered photon direction (t2,phi2)
  !
  ! output:
  ! costnew,sintnew,phinew   ---tnew is scattering angle in scattering frame
  ! ---phinew is azimuthal angle in frame of scattering
  ! see MARS notebook p. 4-7 for drawings of coordinate systems

  ! note:  when calling for pathfinder and MGS coords, the positions are
  ! theta1 and phi1, and the photon angles are theta2 and phi2.

  use constants

  implicit none
  real*8 :: cost1,sint1,cost2,sint2,phi1,phi2,costnew
  real*8 :: sintnew,phinew,diff,cosx,cosp1,cosp2,sinp1,sinp2,ri1,x
  ! integer :: flag
  character(len=20) :: routine
  real*8 :: check

  routine='findangle'

  costnew=cost1*cost2+(sint1*sint2*(cosp1*cosp2+sinp1*sinp2))
  costnew=check(costnew,routine)
  sintnew=sqrt(1.d0-costnew**2)

  ! sint1 and sintnew should both be greater than 0.  but they could be 0.
  if (sint1.gt.1.d-8) then
     if (sintnew.gt.1.d-8) then
        cosx=(cost2-(cost1*costnew))/(sint1*sintnew)
        cosx=check(cosx,routine)
        x=acos(cosx)
        if (phi2.gt.phi1) then
           diff=phi2-phi1
           if (diff.gt.pi) then
              diff=r2p-diff
              ri1=x
              phinew=pi+x
           else
              phinew=pi-x
              ri1=r2p-x
           end if
        else
           diff=phi1-phi2
           if (diff.gt.pi) then
              diff=r2p-diff
              phinew=pi-x
              ri1=r2p-x
           else
              phinew=pi+x
              ri1=x
           end if
        end if
     else
        phinew=0.d0        !should be random, doesn't matter, means at pole
     end if
  else
     phinew=phi2
  end if

  return
end subroutine findangle3
