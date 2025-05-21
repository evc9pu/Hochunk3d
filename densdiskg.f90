subroutine densdisk(rad,sint,cost,phi,dens,id)

  ! calculates density in disk with Gaussian scale height.
  ! called during grid setup so doesn't need to be optimized.

  ! history:
  ! 00/09/06 (baw): write subroutine, modified from opacdiskg.f

  use tts_mod
  use grid_mod, only : is_disk
  use spot_mod
  use constants

  implicit none

  real*8 :: zp,rp,rad,sint,phi,dens,cost,gexp,zmr,warp,warp2,rimcurve,gapcurve
  real*8 :: spiral,adens,bdens,rhodens,sfact,rimpuff,gappuff,densfact,densfact1,zscale,rcoord
  integer :: id,is

  if(.not.is_disk(id)) then
     write(*,*) "ERROR: dust index",id,"is not for a disk"
     stop
  end if

  rp=rad*sint
  zp=rad*cost

  if (idiskcurve.eq.1) then
     rcoord=rad   !surfaces and density follow theta instead of z.  
  else
	 rcoord=rp
  endif

  if (igapd.eq.1.and.rcoord.lt.rgapd2.and.rcoord.gt.rgapd1) then
	adens=agap(id)
	bdens=bgap(id)
	rhodens=rhogap0(id)
	zscale=z1gap(id)
  else
	adens=a(id)
	bdens=b(id)
	rhodens=rho0(id)
	zscale=z1(id)
  endif
	
  ! set up warp structure of disk--set warpheight=1 for no warp
  ! note:  assumes longitude of spots is 0 and 180, because
  ! of the way cos function behaves (+ and - signs)
  warp=0.d0
  warp2=0.d0
  if (idiskwarp.eq.1) then
     if (phi.lt.pi) then
  !      warp=warpheight*(cos(phi/2.d0))**wexp*exp(1.0d0-rp/rddust)
        warp=warpheight*(cos(phi/2.d0))**wexp*exp(-((rcoord-rddust(id))/warplength)**2)
     else
  !      warp=-warpheight*(cos(phi/2.d0))**wexp*exp(1.0d0-rp/rddust)
        warp=-warpheight*(cos(phi/2.d0))**wexp*exp(-((rcoord-rddust(id))/warplength)**2)
     end if
     if (nspot.eq.2) then
 !      warp2=warpheight*(cos(phi/2.d0-pi/2.d0))**wexp*exp(1.0d0-rp/rddust)
       warp2=warpheight*(cos(phi/2.d0-pi/2.d0))**wexp*exp(-((rcoord-rddust(id))/warplength)**2)
     else
        warp2=0.d0
     end if
  else
     warp=0.d0
     warp2=0.d0
  end if

  rimcurve=0.d0
  if (irimcurve.eq.1) then
	rimcurve=rimcurveheight*exp(-((rcoord-rddust(id))/rimcurvelength)**rimcurveexp)   !play around with exponent if you want.  0.5-1.5?
  else
	rimcurve=0
  endif

  rimpuff=0.d0
  if (irimpuff.eq.1) then
	rimpuff=rimheight*exp(-((rcoord-rddust(id))/rimlength)**2)
  else	
	rimpuff=0.d0
  endif	
  if (rcoord.lt.rddust(id)) rimpuff=0.d0

  gappuff=0.d0
  if (igappuff.eq.1.and.igapd.eq.1) then
	gappuff=gapheight*exp(-((rcoord-rgapd2)/gaplength)**2)
  else	
	gappuff=0.d0
  endif	
  if (rcoord.lt.rgapd2) gappuff=0.d0

  gapcurve=0.d0
  if (igapcurve.eq.1.and.igapd.eq.1) then
	gapcurve=gapcurveheight*exp(-((rcoord-rgapd2)/gapcurvelength)**gapcurveexp)
  else	
	gapcurve=0.d0
  endif	
  if (rcoord.lt.rgapd2) gapcurve=0.d0

  spiral=1.d0
  if (ispiral.eq.1) then
  !    spiral=ssc*abs(cos(phi*nspiral/2.d0+(lexp*log(rp/rddust)))**sexp)
  sfact=1.d0
      do is=2,SN,2
		sfact=sfact*(real(is)/(real(is)-1))
  	  end do
	  spiral=1.d0-sw+sfact*sw*sin(log(rcoord)/tan(pitch)-phi+0.7854)**SN
	  if (rcoord.lt.rspiral1.or.rcoord.gt.rspiral2) spiral=1.d0
  else
      spiral=1.d0
  endif

  zmr=zscale*(rcoord/rmin)**bdens

  densfact1=1.d0
  if (phi.lt.(pi/2.d0).or.phi.gt.(3.d0*pi/2.d0)) then
     ! warp on +z side only of spot #1 (i.e., not on bottom)
     if (zp.gt.0.d0) zmr=zmr*(1.d0+warp)
     densfact1=1.d0/(1.d0+warp)
  else
     ! warp on -z side only of spot #2 (not on top)
     if (zp.lt.0.d0) zmr=zmr*(1.d0+warp2)
     densfact1=1.d0/(1.d0+warp2)
  end if

!  zmr=zmr*(rimpuff)
!  zmr=zmr*(gappuff)

  densfact=1.d0/(1.d0+rimpuff+gappuff-rimcurve-gapcurve)*densfact1
  zmr=zmr*(1.d0+rimpuff+gappuff-rimcurve-gapcurve)

 if(rcoord.gt.rmaxd(id).or.rhodens.eq.0.d0.or.rcoord.lt.rddust(id).or.zmr.lt.1.d-10) then
     dens=0.d0
     return
  end if

  gexp=(zp/zmr)**2
  if(gexp.gt.1.d-15) then
     if(gexp.lt.120.d0) then
!        dens=densfact*spiral*rhodens*(rcoord/rmin)**(-adens)*exp(-0.5d0*gexp)* &
!             & (1.d0-sqrt(rmin/rcoord))
        dens=densfact*spiral*rhodens*(rcoord/rmin)**(-adens)*exp(-0.5d0*gexp)
     else
        dens=0.d0
     end if
  else
 !    dens=densfact*spiral*rhodens*(rcoord/rmin)**(-adens)*(1.d0-sqrt(rmin/rcoord))
     dens=densfact*spiral*rhodens*(rcoord/rmin)**(-adens)
  end if

  if (iradexp.eq.1) then
    dens=dens*exp(-rcoord/radexp)
!	print*,'rad,dens',rad,dens
   endif

 
  return

end subroutine densdisk
