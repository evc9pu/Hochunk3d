pro im3col_fits,fileb,fileg,filer,minflux,lrange

imb=readfits(fileb,headb)
img=readfits(fileg,headg)
imr=readfits(filer,headr)

print,headb

nx=fxpar(headb,'NAXIS1')
ny=fxpar(headb,'NAXIS2')

rmaxi=fxpar(headb,'RMAXI')
rmaxi=rmaxi/1.49598e13
print,'rmaxi',rmaxi
pixsizx=2.*rmaxi/float(nx)
pixsizy=pixsizx
print,'pixsizx',pixsizx

unit='intensity'  ;new version of code has already converted to intensity so this will be skipped
if (unit eq 'flux') then begin
; convert from mJy to Jy
    imb=imb/1000.
    img=img/1000.
    imr=imr/1000.
; convert from Jy/pixel to MJy/sr
    fact=1./pixsizx/pixsizy/1.e6*(206265.)^2
    print,'fact,pixsize',fact,pixsizx,pixsizy
    imb=imb*fact
    img=img*fact
    imr=imr*fact
endif

;try to compute min/max from median
medb=median(imb)
medg=median(img)
medr=median(imr)
valmin=min([medb,medg,medr])
print,'estimated valmin,valmax = ',valmin,10^lrange*valmin
valmin=minflux
print,'using valmin,valmax = ',valmin,10^(lrange)*valmin
valmin=alog10(valmin)
valmax=valmin+lrange

;valmin=alog10(valmin)
;valmax=alog10(valmax)

imnew1=bytscl(alog10(imb),min=valmin,max=valmax)
imnew2=bytscl(alog10(img),min=valmin,max=valmax)
imnew3=bytscl(alog10(imr),min=valmin,max=valmax)

a=where(imnew1 gt 255,count)
if count gt 0 then imnew1(a)=255
a=where(imnew2 gt 255,count)
if count gt 0 then imnew2(a)=255
a=where(imnew3 gt 255,count)
if count gt 0 then imnew3(a)=255

;assume images are square!
print, 'assuming images are square!'

half=fix(nx/2)
x=pixsizx*(findgen(nx)-half)
y=pixsizy*(findgen(ny)-half)

xmax=max(x)+pixsizx*0.5
ymax=max(y)+pixsizy*0.5

print,'xmax,ymax',xmax,ymax

;make ps plot

set_plot,'ps'

device,/inches,/color,xsize=5.1,ysize=5.,bits=24,filename='3col.eps',$
	/encapsulated,/CMYK
;device,/inches,/color,xsize=3.1,ysize=3.,bits=24,filename='3col.ps',$
;	xoffset=1.5,yoffset=1.5


contour,imnew1,x,y,xrange=[-xmax,xmax],yrange=[-ymax,ymax],$
xstyle=1+4,ystyle=1+4,$
charsize=1.3,charthick=2,xmargin=[8,3],ymargin=[4,2],$
/nodata
;xticks=4,yticks=4,$
;xstyle=1,ystyle=1,$
;title=str1,$
;charsize=1.,charthick=2,xmargin=[7,3],ymargin=[3,2],$
;xminor=4,yminor=4,xthick=1,ythick=1,/nodata

; Get size of plot window in device pixels.  (got this from an IDL help file)
PX = !X.WINDOW * !D.X_VSIZE
PY = !Y.WINDOW * !D.Y_VSIZE
; Desired size of image in pixels.
SX = PX[1] - PX[0] + 1
SY = PY[1] - PY[0] + 1
;print,sx,sy

sxnew=sx
synew=sy

tv,[[[imnew3]],[[imnew2]],[[imnew1]]],px[0],py[0],xsize=sxnew,ysize=synew,$
	true=3

;levarr=maxim*10^(-3.*findgen(2))
;levarr=reverse(levarr)
;print,levarr

;add axes
axis, xaxis=0,xticks=4,xminor=4,xthick=2,charsize=1.3,xstyle=1,xtitle='AU',charthick=2
axis, xaxis=1,xticks=4,xminor=4,xthick=2,charsize=0.001,xstyle=1
axis, yaxis=0,yticks=4,yminor=4,ythick=2,charsize=1.3,ystyle=1,ytitle='AU',charthick=2
axis, yaxis=1,yticks=4,yminor=4,ythick=2,charsize=0.001,ystyle=1

device,/close

end
