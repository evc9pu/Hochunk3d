subroutine stokes6()

  ! history:
  ! 99/03/22 (baw) put in 6-element matrix, use rejection method only
  ! 95/01/17 (mjw): add calls to ERRMSG,WRIMSG
  ! 95/12/08 (baw): sample from tables instead of analytic functions

  ! This subroutine
  ! calculates scattering using tabular phase functions.
  ! this is for randomly oriented spheres or spheroids, which require
  ! four to six unique matrix elements
  ! placed the following way:
  ! |  S11 S12 0   0  |
  ! |  S12 S22 0   0  |
  ! |  0  0   S33 S34 |
  ! |  0  0  -S34 S44 |
  ! Bohren&Huffman eqn (13.21).
  ! note:  white and bh have sign of S34 switched.  I'm using B-H's
  ! convention here, but in the original Stokes routines, I used
  ! White's.

  ! This matrix is rotated into and out of the scattering frame to
  ! get stokes parameters in observer's frame

  use stokesvar_mod
  use dust_mod
  use random
  use tab_mod
  implicit none

  character(len=70) :: cmsger

  integer :: itab

  real*8 :: bmu,xran,ri1,ri3,cosi3,sini3,costp,sintp,phip,cosb2
  real*8 :: b,sinbt,sini2,bott,cosdph,cosi2,sin2i3,sin2i2,cos2i3,cos2i2
  real*8 :: sin2,cos2,sin2cos1,cos2sin1,s11,s12,s22,s33,s34,s44
  real*8 :: a11,a12,a13,a21,a22
  real*8 :: a23,a24,a31,a32,a33,a34,a42,a43,a44,cosi1,sini1,sin2i1
  real*8 :: cos2i1,a,si,sq,su,sv,rprob

  
  costp=cost
  sintp=sint
  phip=phi

  ! uses rejection method for sampling from
  ! exact phase function; can be inefficient for
  ! forward-peaked dust.

  ! *************
  ! WARNING; we aren't using this subroutine in radiative eq.
  ! if we do, we need to compute peak every time a wavelength
  ! is sampled---it is much faster than using the maximum
  ! peak computed for all wavelengths.  just try using www003 dust
  ! if you don't believe me.
  ! ************

5 continue
  htot=htot+1.d0
  xran=ran()
  bmu=1.d0-2.d0*xran
  itab=int((1.d0+bmu)/dct)+1
  if (itab.lt.1.or.itab.gt.narr) print*,'whoops,itab',itab
  s11=s11arr(itab)
  s12=s12arr(itab)
  cosb2=bmu**2
  b=cosb2-1.d0

  if(abs(bmu).gt.1.d0) then
     write(cmsger,'(a,f)')'bmu eq.1',bmu
     call ERRMSG('WARNING','STOKES',cmsger)
     if(bmu.gt.1.d0) then
        bmu=1.d0
        cosb2=1.d0
        b=0.d0
     else
        bmu=-1.d0
        cosb2=1.d0
        b=0.d0
     end if
  end if

  sinbt=sqrt(1.d0-cosb2)

  xran=ran()
  ri1=r2p*xran

  ! **** ri1 gt pi ****
  ! don't try to combine this with ri1 lt pi.  there are several
  ! slight differences and it will require more than one if statement

  if (ri1.gt.pi) then
     ri3=r2p-ri1
     cosi3=cos(ri3)
     sini3=sin(ri3)
     sin2i3=2.d0*sini3*cosi3
     cos2i3=2.d0*cosi3**2-1.d0
     a11=s11
     a12=s12*cos2i3
     a13=s12*sin2i3
     ! rejection test
     rprob=(a11*sip+a12*sqp+a13*sup)/sip
     if(rprob.gt.peak) then
        ! this shouldn't be the case since we measured peak of p1
        write(cmsger,'(a,f14.11,f14.11)')'rprob gt peak!',
        $		rprob,peak
        call ERRMSG('WARNING','STOKES',cmsger)
        peak=rprob
     end if
     xran=ran()
     if(peak*xran.gt.rprob) go to 5
     hit=hit+1
     a=rprob

     if(bmu.eq.1.d0) then
        write(cmsger,'(a,f)')'bmu eq.1',bmu
        call ERRMSG('WARNING','STOKES',cmsger)
        go to 10
     else if(bmu.eq.-1.d0) then
        write(cmsger,'(a,f)')'bmu eq.-1',bmu
        call ERRMSG('WARNING','STOKES',cmsger)
        sini3=1.d0
        sini2=1.d0
        cosi3=0.d0
        cosi2=0.d0
        cost=-costp
        sint=sintp
     else
        ! anything with p on the end is the primed angle, which is
        ! the angle from whence it came. oops--except cosbwp.
        cost=costp*bmu+sintp*sinbt*cosi3
        ! write(6,*)'costp,bmu,sintp,sinbt,cosi3,cost'
        ! write(6,*)costp,bmu,sintp,sinbt,cosi3,cost
        ! everytime we do an arcsin or arccos we have to test to see
        ! if roundoff error has caused cos or sin to be slightly
        ! greater than 1.  then we have to set them to 1.
        if (abs(cost).lt.1.d0) then
           sint=abs(sqrt(1.d0-cost**2))
           sini2=sini3*sintp/sint
           ! if (sini2.lt.0) write(6,*)'sini2 lt 0',sini2
           bott=sint*sinbt
           cosi2=costp/bott-cost*bmu/bott
           if (abs(cosi2).gt.1.0001d0) then
              write(cmsger,'(a,f)')'cosi2 big',cosi2
              call ERRMSG('WARNING','STOKES',cmsger)
           end if
        else
           sint=0.d0
           call ERRMSG('WARNING','STOKES','sint,sini2=0')
           sini2=0.0d0
           if (cost.ge.1.d0) then
              cosi2=-1.d0
              cost=1.d0
           end if
           if (cost.le.-1.d0) then
              cosi2=1.d0
              cost=-1.d0
           end if
        end if
     end if

     cosdph=-cosi2*cosi3+sini2*sini3*bmu
     if (abs(cosdph).gt.1.d0) then
        if(abs(cosdph).gt.1.0001d0) then
           write(cmsger,'(a,f)')'abs(cosdph).gt.1.0001',cosdph
           call ERRMSG('WARNING','STOKES',cmsger)
        end if
        if(cosdph.gt.1.d0) then
           cosdph=1.0d0
        else
           cosdph=-1.0d0
        end if
     end if
     phi=phip+acos(cosdph)
     if(phi.gt.r2p) phi=phi-r2p
     if(phi.lt.0.d0) phi=phi+r2p

     sin2i2=2.d0*sini2*cosi2
     cos2i2=2.d0*cosi2**2-1.d0
     sin2=sin2i2*sin2i3
     cos2=cos2i2*cos2i3
     sin2cos1=sin2i2*cos2i3
     cos2sin1=cos2i2*sin2i3

     s22=s22arr(itab)
     s33=s33arr(itab)
     s34=s34arr(itab)
     s44=s44arr(itab)

     a21=s12*cos2i2
     a22=s22*cos2-s33*sin2
     a23=s22*cos2sin1+s33*sin2cos1
     a24=s34*sin2i2
     a31=-s12*sin2i2
     a32=-s22*sin2cos1-s33*cos2sin1
     a33=-s22*sin2+s33*cos2
     a34=s34*cos2i2
     a42=s34*sin2i3
     a43=-s34*cos2i3
     a44=s44

     ! **** ri1 lt pi ****
  else

     cosi1=cos(ri1)
     sini1=sin(ri1)
     sin2i1=2.d0*sini1*cosi1
     cos2i1=2.d0*cosi1**2-1.d0
     a11=s11
     a12=s12*cos2i1
     a13=-s12*sin2i1
     ! rejection test
     rprob=(a11*sip+a12*sqp+a13*sup)/sip
     if(rprob.gt.peak) then
        ! this shouldn't be the case since we measured peak of p1
        write(cmsger,'(a,f14.11,f14.11)')'rprob gt peak!',
        $		rprob,peak
        call ERRMSG('WARNING','STOKES',cmsger)
        peak=rprob
     end if
     xran=ran()
     if(peak*xran.gt.rprob) go to 5
     hit=hit+1
     a=rprob

     if(bmu.eq.1.d0) then
        write(cmsger,'(a,f)')'bmu eq.1',bmu
        call ERRMSG('WARNING','STOKES',cmsger)
        go to 10
     else if(bmu.eq.-1.d0) then
        write(cmsger,'(a,f)')'bmu eq.-1',bmu
        call ERRMSG('WARNING','STOKES',cmsger)
        sini1=1.d0
        sini2=1.d0
        cosi3=0.d0
        cosi2=0.d0
        cost=-costp
        sint=sintp
     else

        ! anything with p on the end is the primed angle, which is
        ! the angle from whence it came. oops--except cosbwp.

        cost=costp*bmu+sintp*sinbt*cosi1
        ! write(6,*)'costp,bmu,sintp,sinbt,cosi1,cost'
        ! write(6,*)costp,bmu,sintp,sinbt,cosi1,cost
        ! everytime we do an arcsin or arccos we have to test to see
        ! if roundoff error has caused cos or sin to be slightly
        ! greater than 1.  then we have to set them to 1.
        if (abs(cost).lt.1.d0) then
           sint=abs(sqrt(1.d0-cost**2))
           sini2=sini1*sintp/sint
           ! if (sini2.lt.0) write(6,*)'sini2 lt 0',sini2
           bott=sint*sinbt
           cosi2=costp/bott-cost*bmu/bott
           if (abs(cosi2).gt.1.0001d0) then
              write(cmsger,'(a,f)')'cosi2 big',cosi2
              call ERRMSG('WARNING','STOKES',cmsger)
           end if
        else
           sint=0.d0
           call ERRMSG('WARNING','STOKES','sint,sini2=0')
           sini2=0.0d0
           if (cost.ge.1.d0) then
              cosi2=-1.d0
              cost=1.d0
           end if
           if (cost.le.-1.d0) then
              cosi2=1.d0
              cost=-1.d0
           end if
        end if
     end if
     cosdph=-cosi1*cosi2+sini1*sini2*bmu
     if (abs(cosdph).gt.1.d0) then
        if(abs(cosdph).gt.1.0001d0)then
           write(cmsger,'(a,f)')'abs(cosdph).gt.1.0001',cosdph
           call ERRMSG('WARNING','STOKES',cmsger)
        end if
        if(cosdph.gt.1.d0) then
           cosdph=1.0d0
        else
           cosdph=-1.0d0
        end if
     end if

     phi=phip-acos(cosdph)
     if(phi.gt.r2p) phi=phi-r2p
     if(phi.lt.0.d0) phi=phi+r2p

     sin2i2=2.d0*sini2*cosi2
     cos2i2=2.d0*cosi2**2-1.d0
     sin2=sin2i2*sin2i1
     cos2=cos2i2*cos2i1
     sin2cos1=sin2i2*cos2i1
     cos2sin1=cos2i2*sin2i1

     s22=s22arr(itab)
     s33=s33arr(itab)
     s34=s34arr(itab)
     s44=s44arr(itab)

     a21=s12*cos2i2
     a22=s22*cos2-s33*sin2
     a23=-s22*cos2sin1-s33*sin2cos1
     a24=-s34*sin2i2
     a31=s12*sin2i2
     a32=s22*sin2cos1+s33*cos2sin1
     a33=-s22*sin2+s33*cos2
     a34=s34*cos2i2
     a42=-s34*sin2i1
     a43=-s34*cos2i1
     a44=s44

  end if

  ! sampled from a so divide by a.
  si=(a11*sip+a12*sqp+a13*sup)/a
  sq=(a21*sip+a22*sqp+a23*sup+a24*svp)/a
  su=(a31*sip+a32*sqp+a33*sup+a34*svp)/a
  sv=(a42*sqp+a43*sup+a44*svp)/a

  ! write(6,*)'si,sq,su',si,sq,su,nscat

  sip=si
  sqp=sq
  sup=su
  svp=sv
  cosp=cos(phi)
  sinp=sin(phi)


10 continue
  return
end subroutine stokes6

! ***********************************************************************
