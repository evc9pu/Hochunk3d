subroutine fractal_sph(n,fl,const,fm1,fm2,fm3,fm4,icent)

  ! calculates density in envelope. called during grid setup
  ! so doesn't need to be optimized.

  ! history:
  ! 00/09/06 (baw): write subroutine, modified from opacenv.f
  ! 01/02/17 (baw): combine opacitybub subroutine into here.
! 2012/03/28 (mjw): replace calls to RAN() with call to RAN_OLD()


 !  use tts_mod
  use grid_mod
  use random
 ! use opacin_mod
  use constants, only : pi

  implicit none

  real*8 :: r,rmin,th,phi,rmax2,rl,const
  real*8 :: fl,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,fm1,fm2,fm3,fm4
  real*8 :: x5,y5,z5,x6,y6,z6,rmax1,xmax,xmin,ymax,ymin,zmax,zmin

  integer :: ir,it,ip,i1,i2,i3,i4,n,i5,i6,icent,id


      rl=1./fl
      
   !  this only works if the array is a density array.  
   !   do ir=1,nrg
   !      do it=1,ntg
   !         do ip=1,npg
   !            densarr(ir,it,ip,id)=const
   !         end do
   !      end do
   !   end do

      rmax1=(0.5+0.5*rl+(0.5*rl)**2+(0.5*rl)**3+(0.5*rl)**4)*sqrt(2.)
      xmin=0.
      xmax=0.
      ymin=0.
      ymax=0.
      zmin=0.
      zmax=0.
      rmax2=0.
      
      print*,'rmax1',rmax1

      do i1=1,n
        
   !     print*,'i1,n',i1,n

   !  x1,y1,z1 range from -.5 to .5
         if (icent.eq.0) then
            x1=ran_old()-0.5
            y1=ran_old()-0.5
            z1=ran_old()-0.5
         else
            x1=0.
            y1=0.
            z1=0.
         endif
         do i2=1,n
            x2=x1+(ran_old()-.5)*rl*fm1
            y2=y1+(ran_old()-.5)*rl*fm1
            z2=z1+(ran_old()-.5)*rl*fm1
            do i3=1,n
               x3=x2+(ran_old()-.5)*rl**2*fm2
               y3=y2+(ran_old()-.5)*rl**2*fm2
               z3=z2+(ran_old()-.5)*rl**2*fm2
               do i4=1,n
                  x4=x3+(ran_old()-.5)*rl**3*fm3
                  y4=y3+(ran_old()-.5)*rl**3*fm3
                  z4=z3+(ran_old()-.5)*rl**3*fm3
                  do i5=1,n
                     x5=x4+(ran_old()-.5)*rl**4*fm4
                     y5=y4+(ran_old()-.5)*rl**4*fm4
                     z5=z4+(ran_old()-.5)*rl**4*fm4
   !                  do i6=1,n
   !  fm4=1.
   !                     x6=x5+(ran_old()-.5)*rl**5*fm4
   !                     y6=y5+(ran_old()-.5)*rl**5*fm4
   !                     z6=z5+(ran_old()-.5)*rl**5*fm4

                     
   !  07/24/03 baw took this out.  can put it back in if I ever think I
   !  understand it but I don't think it hurts to take it out.  this makes
   !  it into 1 cube which repeats.  I worry it adds more density to
   !  the edges.
   !               if(x4.lt.3.e-5)then
   !                  x4=x4+.99999
   !               end if
   !               if(y4.lt.3.e-5)then
   !                  y4=y4+.99999
   !               endif
   !               if(z4.lt.3.e-5)then  !change only if 1st box
   !                  z4=z4+.99999
   !               end if
   !               if(x4.gt.0.99997)x4=x4-.9999 
   !               if(y4.gt.0.99997)y4=y4-.9999   
   !               if(z4.gt.0.99997)z4=z4-.9999 !last box

   !               ix=nxg*x4+1
   !               iy=nyg*y4+1
   !               iz=nzg*z4+1
   !  now that we've found x4,y4,z4, calculate r,theta,phi, and then
   !  ir,it,ip
   !  sphere center is at (x,y,z)=(0.5,0.5,0.5)
   !               r=(x4-0.5)**2+(y4-0.5)**2+(z4-0.5)**2
   !               r=sqrt(r)
   !               th=acos((z4-0.5)/r)
   !               phi=atan2((x4-0.5),(y4-0.5))

                     if (x5.lt.xmin) xmin=x5
                     if (x5.gt.xmax) xmax=x5
                     if (y5.lt.ymin) ymin=y5
                     if (y5.gt.ymax) ymax=y5
                     if (z5.lt.zmin) zmin=z5
                     if (z5.gt.zmax) zmax=z5

                     r=x5**2+y5**2+z5**2
                     r=sqrt(r)

                     if (r.gt.rmax2) rmax2=r
                     th=acos(z5/r)
                     phi=atan2(y5,x5)
 
                     if (phi.lt.0.) phi=phi+2.0*pi

                     if (icent.eq.0) then
                        r=r*rarr(nrg)*2.*1./(1.+rl)
                     else
                        r=r*rarr(nrg)*4.
                     endif
                     call locate(rarr,nrg,r,ir)
                     if (ir.eq.0) ir=1


                     call locate(thetarr,ntg,th,it)
   !  it=int(float(ntg)*th/pi)+1
                     call locate(phiarr,npg,phi,ip)
                     if (ip.eq.0) then
                        print*,'warning, fractal_sph, ip should be > 0'
                        ip=1
                     endif
                     if (it.eq.0) then
                        print*,'warning, fractal_sph, it should be > 0'
                        it=1
                     endif    
                     
    !                 print*,'ir,it,ip',ir,it,ip,nrg,ntg,npg              
   !  ip=int(float(npg)*phi/pi/2.0)+1
                     if (ir.lt.1.or.it.lt.1.or.ip.lt.1.or.it.gt.ntg.or.ip.gt.npg) then
   !  okay if ir > nrg since it's filling a cube
                        print*,'indices wrong in fractal_sph',ir,it,ip
                        stop
                      endif
                     if (ir.le.nrg) densvararr(ir,it,ip)=densvararr(ir,it,ip)+1.
   !  scale density to be small to keep from overflows later on...
   !  enddo !do i6
                  end do        !do i5
               end do           !do i4
            end do              !do i3
         end do                 !do i2
      end do                    !do i1
      
      print*,'range in x,y,z'
      print*,xmin,xmax,ymin,ymax,zmin,zmax
      print*,'rmax vs calculated',rmax2,rmax1
      

return

end subroutine fractal_sph

