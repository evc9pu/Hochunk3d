set_plot,'PS'

if !d.name eq 'PS' then device,xsize=7.,ysize=7.,/inches,$
/portrait,xoffset=1.,yoffset=1.,/color,bits=24,filename='diffus.ps'

nr=1L
nt=1L
np=1L

dir='../models/test2/'

close,1
openr,1,dir+'diffdir.unf',/f77_unformatted
readu,1,nr,nt,np
r=dblarr(nr)
theta=dblarr(nt)
phi=dblarr(np)
temp=lonarr(nr,nt,np)
readu,1,r
readu,1,theta
readu,1,phi
readu,1,temp
close,1
;temp=convol(temp,k,/edge_truncate)

close,1
openr,1,dir+'diffuse.unf',/f77_unformatted
readu,1,nr,nt,np
r=dblarr(nr)
theta=dblarr(nt)
phi=dblarr(np)
rho=lonarr(nr,nt,np)
readu,1,r
readu,1,theta
readu,1,phi
readu,1,rho
close,1


rmin=.05
rmax=30.
pi=3.14159265

theta=theta-pi/2.

r1=fltarr(nr,nt)
theta1=fltarr(nr,nt)

x=fltarr(nr,nt)
y=fltarr(nr,nt)
for i=0,nr-1 do begin
  for j=0,nt-1 do begin
    x(i,j)=r(i)*sin(theta(j))
    y(i,j)=r(i)*cos(theta(j))
    r1(i,j)=r(i)
    theta1(i,j)=theta(j)
  endfor
endfor

loadct,5

rim1=.25
sx=rim1/300.
sy=rim1/300.
t_int1=polar_surface(temp(*,*),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim1/2,rim1,rim1/2])
r_int1=polar_surface(rho(*,*),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim1/2,rim1,rim1/2])

rim2=2.5
sx=rim2/300.
sy=rim2/300.
t_int2=polar_surface(temp(*,*),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim2/2,rim2,rim2/2])
r_int2=polar_surface(rho(*,*),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim2/2,rim2,rim2/2])

rim3=25.
sx=rim3/300.
sy=rim3/300.
t_int3=polar_surface(temp(*,*),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim3/2,rim3,rim3/2])
r_int3=polar_surface(rho(*,*),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim3/2,rim3,rim3/2])

tvscl,t_int1,0,5.,xsize=2.,ysize=2.,/inches
tvscl,t_int2,2.5,5.,xsize=2.,ysize=2.,/inches
tvscl,t_int3,5.,5.,xsize=2.,ysize=2.,/inches

tvscl,r_int1,0,2.5,xsize=2.,ysize=2.,/inches
tvscl,r_int2,2.5,2.5,xsize=2.,ysize=2.,/inches
tvscl,r_int3,5.,2.5,xsize=2.,ysize=2.,/inches

x1=findgen(301)*rim1/301
y1=findgen(301)*rim1/301-rim1/2

x2=findgen(301)*rim2/201
y2=findgen(301)*rim2/201-rim2/2

x3=findgen(301)*rim3/301
y3=findgen(301)*rim3/301-rim3/2

contour,t_int1,x1,y1, POSITION=[0., 0.714, 0.286, 1.],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim1],$
yrange=[-rim1/2,rim1/2],$
xminor=2,yminor=2,xticks=3,yticks=3

contour,t_int1,x1,y1, POSITION=[0., 0.714, 0.286, 1.],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim1],$
yrange=[-rim1/2,rim1/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2,xticks=3,yticks=3

contour,t_int2,x2,y2, POSITION=[0.357, 0.714, 0.643, 1.],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim2],$
yrange=[-rim2/2,rim2/2],$
xminor=2,yminor=2

contour,t_int2,x2,y2, POSITION=[0.357, 0.714, 0.643, 1.],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim2],$
yrange=[-rim2/2,rim2/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2

contour,t_int3,x3,y3, POSITION=[0.714, 0.714, 1., 1.],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim3],$
yrange=[-rim3/2,rim3/2],$
xminor=2,yminor=2,xticks=3,yticks=3

contour,t_int3,x3,y3, POSITION=[0.714, 0.714, 1., 1.],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim3],$
yrange=[-rim3/2,rim3/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2,xticks=3,yticks=3


contour,r_int1,x1,y1, POSITION=[0., 0.357, 0.286, 0.643],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim1],$
yrange=[-rim1/2,rim1/2],$
xminor=2,yminor=2,xticks=3,yticks=3

contour,r_int1,x1,y1, POSITION=[0., 0.357, 0.286, 0.643],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim1],$
yrange=[-rim1/2,rim1/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2,xticks=3,yticks=3

contour,r_int2,x2,y2, POSITION=[0.357, 0.357, 0.643, 0.643],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim2],$
yrange=[-rim2/2,rim2/2],$
xminor=2,yminor=2

contour,r_int2,x2,y2, POSITION=[0.357, 0.357, 0.643, 0.643],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim2],$
yrange=[-rim2/2,rim2/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2

contour,r_int3,x3,y3, POSITION=[0.714, 0.357, 1., 0.643],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim3],$
yrange=[-rim3/2,rim3/2],$
xminor=2,yminor=2,xticks=3,yticks=3

contour,r_int3,x3,y3, POSITION=[0.714, 0.357, 1., 0.643],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim3],$
yrange=[-rim3/2,rim3/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2,xticks=3,yticks=3


;colorbar, NCOLORS=255, POSITION=[0.05, 0.2, 0.45, 0.25], range=[0,300],divisions=3,$
;title='T (Kelvin)'

;colorbar, NCOLORS=255, POSITION=[0.55, 0.2, 0.95, 0.25], range=[0,300],divisions=3,$
;title='!4q!3'




;device,xoffset=3.5,yoffset=3.5,xsize=3,ysize=3,/inches
;tvscl,t_int1


;image_cont,t_int1

;colorbar, NCOLORS=100, POSITION=[0.15, 0.85, 0.95, 0.90]

;colorbar, NCOLORS=255, POSITION=[0.1, 0.8, 0.9, 0.90], range=[0,300],divisions=3

;tvscl,congrid(r_int1^0.1,400,400),0,2.7,/inches,/nan
;tvscl,congrid(r_int2^0.1,400,400),2.6,2.7,/inches,/nan
;tvscl,congrid(r_int3^0.1,400,400),5.2,2.7,/inches,/nan


;if !d.name eq 'PS' then device,xsize=2.5,ysize=2.,/inches
;r=findgen(300)*255/299
;cb=r#replicate(1,30)
;;tvscl,cb,1.75,0.4,/inches
;tvscl,cb,0.25,-0.2,/inches

;if !d.name eq 'PS' then device,xsize=11.,ysize=8.5,/inches,$
;/portrait,xoffset=1.,yoffset=1.
;pos=[0.02,-0.01,0.5,.5]
;empty=replicate(' ',30)
;plot,/nodata,[0],yr=[0,1],xr=[0,255],xstyle=1,ystyle=1,/noerase,$
;pos=pos,yticks=1,ytickname=empty
 
;if !d.name eq 'PS' then device,xsize=11.,ysize=8.5,/inches,$
;/portrait,xoffset=0.,yoffset=0,/color,bits=24
;;pos=[0.16,0.048,0.545,0.095]
;pos=[0,0,5,5]
;empty=replicate(' ',30)
;plot,/nodata,[0],yr=[0,200],xr=[0,1],/xstyle,/ystyle,/noerase,$
;pos=pos,yticks=1,ytickname=empty


if !d.name eq 'PS' then device,/close
set_plot,'x'

end




