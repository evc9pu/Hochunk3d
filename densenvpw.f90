subroutine densenvpw(rad,thetin,costin,sina,phi,dens,radamb,id,inhole,rmincgs)

  ! calculates density in envelope. called during grid setup
  ! so doesn't need to be optimized.

  ! history:
  ! 00/09/06 (baw): write subroutine, modified from opacenv.f
  ! 01/02/17 (baw): combine opacitybub subroutine into here.

  use tts_mod
  use grid_mod, only : is_envelope, is_cavity
  use opacin_mod
  use constants

  implicit none

  real*8 :: rad,dens,cosa,phi,thet,thetin,rmincgs
  real*8 :: costin,sina,xmu,rx,rp,zup,xmu0,factor
  real*8 :: rx2,xmu0new,zlo,zp,radamb

  integer :: iflag,id,inhole

  if(.not.is_envelope(id).and..not.is_cavity(id)) then
     write(*,*) "ERROR: dust index",id,"is not for an envelope or cavity"
     stop
  end if

  inhole=0

  ! for now, assuming +z = -z
  thet=thetin
  cosa=abs(costin)
  if (thet.gt.pihalf) thet=pi-thet
  rp=rad*sina
  zp=rad*cosa

  zup=c1e*rp**ex1+z01

  if (rad.gt.rmine) then
	dens=rhodens1*(rad)**(-envexp)
	!rmincgs converts rad (in rstar) to cm
  else
    dens=0.d0
  endif

! constant density in gap--will be overwritten in hole and by ambient density (I want this)
  if (igape.eq.1.and.igapedens.eq.0) then
! igapedens=1 means scaling density in gap; igapedens=0 means density = constant in gap
  if (rad.lt.rgape2.and.rad.gt.rgape1) then
    dens=rhogape
  endif
  endif

  ! compare to ambient density; if rho lt rhoamb, set rho to rhoamb.
  if(rad.gt.rmine.and.dens.lt.rhoamb) then
     dens=rhoamb
     if (rad.lt.radamb) radamb=rad
  end if

  if(ihole.eq.1) then
        zlo=c2e*rp**ex2+z02
        ! baw: 2/5/99 new.  doing this so holes can intersect and create
        ! different shapes!
        if (zp.gt.zlo) then
           inhole=1
           dens=rhoconst2*rad**(-exf)
           ! if(zp.lt.zflowmin) dens=0.d0
        end if
        if(zp.gt.zup) then
           inhole=1
           dens=rhoconst1*rad**(-exf)
           ! if(zp.lt.zflowmin) dens=0.d0
        end if
     if(ibub.eq.1) then
        if (rad.lt.zbub2*xmu**nbub) then
           ! idust2=4
           inhole=1
           if (rad.lt.zbub1*xmu**nbub) then
              dens=rhoconst1
           else
              dens=rhoconst2
           end if
        end if
        ! if(rp.lt.roa.and.rad.gt..8d0*zbub1) dens=rhoconst1
        if(xmu.gt.cosbuboa) dens=rhoconst1
     end if
  end if

!  scaled density in gap
  if (igape.eq.1.and.igapedens.eq.1) then
  if (rad.lt.rgape2.and.rad.gt.rgape1) then
     dens=dens*fracte 
  endif
  endif

  ! make density stop at r>rmax and r<rmine
  if (rad.ge.rmax) dens=0.d0
  if (rad.le.rmine) dens=0.d0

  ! set cavity or envelope density based on dust type
  if (is_cavity(id).and.inhole.eq.0) dens=0.d0
  if (is_envelope(id).and.inhole.eq.1) dens=0.d0

  return

end subroutine densenvpw

