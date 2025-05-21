set_plot,'ps'
device,xsize=5.35,ysize=5.,/inches,xoffset=0.4,yoffset=2.,$
        filename='polmap.ps'

;for models without fitsio
;make fits files first (@runfits).  this has pixscale info in it.
file1='e_2H_80.0_0.0_I_img.fits'
file2='e_2H_80.0_0.0_Q_img.fits'
file3='e_2H_80.0_0.0_U_img.fits'

!p.multi=[0,1,1]

;number of contour levels to plot
nlev=8   ;plots contours in half-mag intervals
outerfrac=1.0   ;zoom of image.  fraction of outer radius. 1=rmaxi

tmpi=readfits(file1,headi)
tmpq=readfits(file2,headq)
tmpu=readfits(file3,headu)

nbin=5
nx=155  ;=5x31   assumes nx=149, need a slightly bigger image divisible by 5
;nbin=3
;nx=153  ;=3x51   assumes nx=149, image divisible by 3
datai=fltarr(nx,nx)
dataq=fltarr(nx,nx)
datau=fltarr(nx,nx)

nx2=fxpar(headi,'NAXIS1')  
;model images are square.

n1=(nx-nx2)/2
n2=nx-n1-1
print,n1,n2

datai(n1:n2,n1:n2)=tmpi
dataq(n1:n2,n1:n2)=tmpq
datau(n1:n2,n1:n2)=tmpu

ang=0.
;ang=135.
;datai=rot(datai,ang,/interp,missing=0)
;dataq=rot(dataq,ang,/interp,missing=0)
;datau=rot(datau,ang,/interp,missing=0)

pixsize=sxpar(headi,'PXSCAL1')
;for models, it's square.

;did convolution already

data1=datai

npol=nx/nbin
data2=rebin(data1,npol,npol)
datai2=rebin(datai,npol,npol)
dataq2=rebin(dataq,npol,npol)
datau2=rebin(datau,npol,npol)

maxv=max(data1)
;half-mag intervals
levarr=maxv*10^(-.2*findgen(nlev))
levarr=reverse(levarr)
psig=levarr(0)
p=fltarr(npol,npol)
thet=fltarr(npol,npol)

;pi=3.14159265

itot=0.
qtot=0.
utot=0.

for i=0,npol-1 do begin
for j=0,npol-1 do begin
   if (datai2(i,j) gt 1.e-15) then p(i,j)=sqrt(dataq2(i,j)^2+datau2(i,j)^2)/datai2(i,j)
   if (abs(dataq2(i,j)) gt 1.e-20) then thet(i,j)=0.5*atan(datau2(i,j),dataq2(i,j))
   thet(i,j) = thet(i,j)+ang*!pi/180.
   if (thet(i,j) gt !pi) then thet(i,j) = thet(i,j) - !pi
endfor
endfor

half=fix(nx/2)
x=pixsize*(findgen(nx)-half)
y=pixsize*(findgen(nx)-half)

prmax=max(x)*outerfrac
contour,data1,x,y,levels=levarr,xrange=[-prmax,prmax],yrange=[-prmax,prmax],$
xticks=4,yticks=4,$
xstyle=1,ystyle=1,$
title='',xtitle='Arcsec',ytitle='Arcsec',$
charsize=1.3,charthick=2,$
xminor=4,yminor=4,xthick=1,ythick=1,/nodata
contour,data1,x,y,levels=levarr,/overplot,$
color=100

half=fix(npol/2)
pixsizepol=float(nx)/float(npol)*pixsize
xpol=pixsizepol*(findgen(npol)-half)
ypol=pixsizepol*(findgen(npol)-half)

const2=3.
pcount=0
tcount=0
pmax=0
for i=0,npol-1 do begin
for j=0,npol-1 do begin
  j2=j
  i2=i
  if(data2(i2,j2) gt psig and p(i2,j2) ge 1.) then pcount=pcount+1
  if(data2(i2,j2) gt psig) then tcount=tcount+1
  if(data2(i2,j2) gt psig and p(i2,j2) le 1.) then begin
  if p(i2,j2) gt pmax then pmax=p(i2,j2)
  c1=sin(1.*(thet(i2,j2)))
  c2=cos(1.*(thet(i2,j2)))
  xmin=-p(i2,j2)*c1*const2
  xmax=p(i2,j2)*c1*const2
  ymin=-p(i2,j2)*c2*const2
  ymax=p(i2,j2)*c2*const2
  usersym,[xmin,xmax],[ymin,ymax],thick=2
;  usersym,[-.707,.707],[-.707,.707]
  plots,xpol(i2),ypol(j2),psym=8,$
;  color=31,$
  thick=10
  endif
endfor
endfor
print,'fraction of data points with P>100%',float(pcount)/float(tcount)

ilab=1
if ilab eq 1 then begin
   phey=const2*.5
   usersym,[-phey,phey],[2,2]
   plots,.63*prmax,0.76*prmax,psym=8,$
;   plots,1.2*prmax,0.76*prmax,psym=8,$
   thick=5
   xyouts,.68*prmax,.8*prmax,'50%',charthick=2,charsize=1.
;   xyouts,1.25*prmax,.8*prmax,'50%',charthick=2,charsize=1.
   pmax=fix(100*pmax)
;   labpm=strcompress('Pmax='+string(pmax)+'%',/remove_all)
;   xyouts,-.92*prmax,-.92*prmax,labpm,charthick=2,charsize=.75
endif

device,/close

end

