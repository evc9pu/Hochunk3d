;pro trho3d,inphi

; set these variables

;iave=1   ;average density and temperature over all dust types
iave=0 ; don't average density and temperature; set dust type to plot (idust)

if iave eq 1 then begin
  ibeg=0  ;  don't change this if iave=1
  iend=0  ;  don't change this if iave=1
endif else begin
; can set these however you want from 0-7
      ;0 is large thermal grain, flattened disk
      ;1 is smaller thermal grain, thicker disk
      ;2 is envelope
      ;3 is cavity
      ;4-7 is same as 0-3 but for PAH/vsg nonthermal grains
  ibeg=0   
  iend=1  
endelse

;model directory
;dir='../models/modcII/'
dir='../models/john/'

inphi=0  ; azimuthal viewing index, set to 0 if 2-D; 0-np if 3-D.

;max radius to plot in 3 windows (left middle right)
rmax1=55.  ;makes 3 plots with different outer radii
rmax2=1055.
rmax3=9511.

;disk model, e.g., modcII
rmax1=0.1  ;makes 3 plots with different outer radii
rmax2=55.
rmax3=511.

;ambient density;  will set T and rho to 0 below this density (gm/cm^3)
;rhoamb=1.e-25 
rhoamb=0. 
;note:  n_H2 = rhoamb/2/1.67e-24.  

;min max temperature to plot
mintin=10.d0   ;min temperature plotted
maxtin=1600.d0 ; max  ' '

;okay that's it for user input

set_plot,'PS'

if !d.name eq 'PS' then device,xsize=8.,ysize=5.5,/inches,$
/portrait,xoffset=1.,yoffset=1.,/color,bits=24,filename='trho_rz.ps'

for idust=ibeg,iend do begin

nr=1L
nt=1L
np=1L
nd=1l

close,1
openr,1,dir+'tarr.unf',/f77_unformatted
readu,1,nr,nt,np,nd
r=dblarr(nr)
theta=dblarr(nt)
phi=dblarr(np)
tmp=dblarr(nr,nt,np,nd)
readu,1,r
readu,1,theta
readu,1,phi
readu,1,tmp
close,1
;temp=dblarr(nr,nt,np)
;temp(0:nr-1,0:nt-1,0:np-1)=tmp(0:nr-1,0:nt-1,0:np-1,id)
tmptemp=tmp
tempd=tmp(*,*,*,idust)

print,'rmax',max(r)

a=where(phi*180./!pi ge inphi)
iphi=a(0)
print,inphi,iphi,np

close,1
openr,1,dir+'darr.unf',/f77_unformatted
readu,1,nr,nt,np,nd
r=dblarr(nr)
theta=dblarr(nt)
phi=dblarr(np)
tmp=dblarr(nr,nt,np,nd)
readu,1,r
readu,1,theta
readu,1,phi
readu,1,tmp
close,1
;rho=dblarr(nr,nt,np)
;rho(0:nr-1,0:nt-1,0:np-1)=tmp(0:nr-1,0:nt-1,0:np-1,id)
rhoave=total(tmp,4)
rhod=tmp(*,*,*,idust)
tmprho=tmp

;do a density-averaged temperature
tmp=0.
a=where(tmptemp gt 4.)
if n_elements(a) gt 1 then begin
   tmp=tmptemp*tmprho
   n=fltarr(nr,nt,np,nd)
   n(a)=1.
   top=total(tmp,4)
   bot=total(n,4)
   a=where(rhoave gt 0.)
   tempave=fltarr(nr,nt,np)
   tempave(a)=top(a)/bot(a)/rhoave(a)
 endif else begin
   tempave=tmptemp
endelse

tmp=1*0
tmptemp=0.
tmprho=0.
top=0.
bot=0.

;set rho, temp to be average or single disk temp
if (iave eq 1) then begin
  rho=rhoave
  temp=tempave
endif else begin
  rho=rhod
  temp=tempd
endelse

print,'min,max temp',min(temp),max(temp)
print,'min,max rho',min(rho),max(rho)
;if you don't want to see temperature of ambient material, do this
a=where(rho lt 1.001*rhoamb)
if a gt 0 then begin
 temp(a)=0.d0
 rho(a)=0.d0
endif

r(0)=r(1)-0.01*r(1)

mint=alog10(mintin)
maxt=alog10(maxtin)
;minr=alog10(minr)
;maxr=alog10(maxr)
maxr=alog10(max(rho))
;minr=min(rho(where(rho gt 0)))*0.9
minr=maxr-12
print,maxr,minr
temp=bytscl(alog10(temp),min=mint,max=maxt)
rho=bytscl(alog10(rho),min=minr,max=maxr)
print,min(rho),max(rho)
pi=3.14159265

;this is technically incorrect since the grid cell is centered
;slightly off the poles, but it gives a prettier picture...
theta(0)=0.
theta(nt-1)=pi

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

;iphi=0

rim1=rmax1
sx=rim1/300.
sy=rim1/300.
t_int1=polar_surface(temp(*,*,iphi),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim1/2,rim1,rim1/2])
r_int1=polar_surface(rho(*,*,iphi),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim1/2,rim1,rim1/2])

rim2=rmax2
sx=rim2/300.
sy=rim2/300.
t_int2=polar_surface(temp(*,*,iphi),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim2/2,rim2,rim2/2])
r_int2=polar_surface(rho(*,*,iphi),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim2/2,rim2,rim2/2])

rim3=rmax3
sx=rim3/300.
sy=rim3/300.
t_int3=polar_surface(temp(*,*,iphi),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim3/2,rim3,rim3/2])
r_int3=polar_surface(rho(*,*,iphi),r1,theta1,spacing=[sx,sy],$
bounds=[0,-rim3/2,rim3,rim3/2])

;print,'min(t_int1),max(t_int1)',min(t_int1),max(t_int1)

tv,t_int1,0.5,3.,xsize=2.,ysize=2.,/inches
tv,t_int2,3.,3.,xsize=2.,ysize=2.,/inches
tv,t_int3,5.5,3.,xsize=2.,ysize=2.,/inches

tv,r_int1,0.5,0.5,xsize=2.,ysize=2.,/inches
tv,r_int2,3.,0.5,xsize=2.,ysize=2.,/inches
tv,r_int3,5.5,0.5,xsize=2.,ysize=2.,/inches


x1=findgen(301)*rim1/301
y1=findgen(301)*rim1/301-rim1/2

x2=findgen(301)*rim2/301
y2=findgen(301)*rim2/301-rim2/2

x3=findgen(301)*rim3/301
y3=findgen(301)*rim3/301-rim3/2


contour,t_int1,x1,y1, POSITION=[0.0625, 0.5454, 0.3125, .9090],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim1],$
yrange=[-rim1/2,rim1/2],$
xminor=2,yminor=2,xticks=3,yticks=3,charsize=0.7,charthick=2.

contour,t_int1,x1,y1, POSITION=[0.0625, 0.5454, 0.3125, .9090],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim1],$
yrange=[-rim1/2,rim1/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2,xticks=3,yticks=3

contour,t_int2,x2,y2, POSITION=[0.375, 0.5454, 0.625, .9090],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim2],$
yrange=[-rim2/2,rim2/2],$
xminor=2,yminor=2,title='temperature',charsize=0.7,charthick=2.

contour,t_int2,x2,y2, POSITION=[0.375, 0.5454, 0.625, .9090],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim2],$
yrange=[-rim2/2,rim2/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2

contour,t_int3,x3,y3, POSITION=[0.6875, 0.5454, .9375, .9090],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim3],$
yrange=[-rim3/2,rim3/2],$
xminor=2,yminor=2,xticks=3,yticks=3,charsize=0.7,charthick=2.


contour,t_int3,x3,y3, POSITION=[0.6875, 0.5454, .9375, .9090],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim3],$
yrange=[-rim3/2,rim3/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2,xticks=3,yticks=3

barposition=[0.945,0.55,0.955,0.9]
colorbarb,position=barposition,min=mint,max=maxt,$
/vertical,format='(F5.1)',/right,CHARSIZE=0.7,ncolors=256


contour,r_int1,x1,y1, POSITION=[0.0625, 0.0909, 0.3125, 0.4545],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim1],$
yrange=[-rim1/2,rim1/2],$
xminor=2,yminor=2,xticks=3,yticks=3,charsize=0.7,charthick=2.

contour,r_int1,x1,y1, POSITION=[0.0625, 0.0909, 0.3125, 0.4545],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim1],$
yrange=[-rim1/2,rim1/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2,xticks=3,yticks=3

contour,r_int2,x2,y2, POSITION=[0.375, 0.0909, 0.625, 0.4545],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim2],$
yrange=[-rim2/2,rim2/2],$
xminor=2,yminor=2,title='density',charsize=0.7,charthick=2.

contour,r_int2,x2,y2, POSITION=[0.375, 0.0909, 0.625, 0.4545],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim2],$
yrange=[-rim2/2,rim2/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2

contour,r_int3,x3,y3, POSITION=[0.6875, 0.0909, 0.9375, 0.4545],$
color=0,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim3],$
yrange=[-rim3/2,rim3/2],$
xminor=2,yminor=2,xticks=3,yticks=3,charsize=0.7,charthick=2.

contour,r_int3,x3,y3, POSITION=[0.6875, 0.0909, 0.9375, 0.4545],$
color=250,/noerase,/nodata,xstyle=1,ystyle=1,xrange=[0,rim3],$
yrange=[-rim3/2,rim3/2],$
xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],$
ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' '],xminor=2,yminor=2,xticks=3,yticks=3

barposition=[0.945,0.1,0.955,0.45]
colorbarb,position=barposition,min=minr,max=maxr,$
/vertical,format='(F5.1)',/right,CHARSIZE=0.7,ncolors=256

erase
endfor

if !d.name eq 'PS' then device,/close
set_plot,'x'

end




