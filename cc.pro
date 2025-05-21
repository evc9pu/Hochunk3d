;cc.pro

;first run maketab to calculate filter fluxes from SED models.

!p.multi=[0,1,1]
loadct,13
set_plot,'ps'
device,xsize=5.,ysize=4.,/inches,$
/portrait,xoffset=1,yoffset=1,/color,bits=24,$
filename='cc.ps'

;for fitsio version of code
readcol,'filterfluxes.txt',cosi,phi,j,h,k,b1,b2,b3,b4,b5,b6,b7,skipline=2
;for non-fitsio version of code
;readcol,'filterfluxes.txt',cosi,j,h,k,b1,b2,b3,b4,b5,b6,b7,skipline=2
ninc=n_elements(j)
nfilt=10

;compute reddening line from star fluxes
filtlam=[1.2,1.65,2.16,3.6,4.5,5.8,8.0,24.0,70.0,160.0]
Av=findgen(100)*3.
readcol,'../models/parfiles/r400_ice095.par',klam4,junk1,junk2,kkap4
;readcol,'r400.par',klam,junk1,junk2,kkap
readcol,'../models/parfiles/kmh.par',klam,junk1,junk2,kkap
filtkap=interpol(kkap,klam,filtlam)
filtkap4=interpol(kkap4,klam4,filtlam)
kapv=interpol(kkap,klam,0.55)
kapv4=interpol(kkap4,klam4,0.55)
red=fltarr(100,nfilt,2)

;Jys at 0 mag at each wavelength (from Cohen)
filtfb=[1594.,1024.,666.7,280.9,179.7,115.0,64.13,7.160,0.7972,0.160]
;when ready, pacs/spires are ,0.6751077,0.3107111,0.1281016,5.8182828E-02,2.7805042E-02,1.3336738E-02]

;convert everything to magnitudes
;used d=10pc so these are absolute mags
j=-2.5*alog10(j/filtfb(0))
h=-2.5*alog10(h/filtfb(1))
k=-2.5*alog10(k/filtfb(2))
b1=-2.5*alog10(b1/filtfb(3))
b2=-2.5*alog10(b2/filtfb(4))
b3=-2.5*alog10(b3/filtfb(5))
b4=-2.5*alog10(b4/filtfb(6))
b5=-2.5*alog10(b5/filtfb(7))
b6=-2.5*alog10(B6/filtfb(8))
b7=-2.5*alog10(b7/filtfb(9))

for i=0,nfilt-1 do begin
   red(*,i,0)=Av*filtkap(i)/kapv
   red(*,i,1)=Av*filtkap4(i)/kapv4
endfor

;pick a distance
d=140.
j=j-5+5*alog10(d)
h=h-5+5*alog10(d)
k=k-5+5*alog10(d)
b1=b1-5+5*alog10(d)
b2=b2-5+5*alog10(d)
b3=b3-5+5*alog10(d)
b4=b4-5+5*alog10(d)
b5=b5-5+5*alog10(d)
b6=b6-5+5*alog10(d)
b7=b7-5+5*alog10(d)

;1 microJy sensitivity limit for 100 Lsun source 
;i.e., assume our source is 100 times brighter and
;compare to 1 uJy.  or compare our 1Lsun source to 1.e-8 Jy.
fmin=-2.5*alog10(1.e-8/filtfb)

;don't plot fluxes below fmin
a=where(j gt fmin(0),count)
if count ne 0 then j(a)=-30000
a=where(h gt fmin(1),count)
if count ne 0 then h(a)=-29000
a=where(k gt fmin(2),count)
if count ne 0 then k(a)=-28000
a=where(b1 gt fmin(3),count)
if count ne 0 then b1(a)=-27000
a=where(b2 gt fmin(4),count)
if count ne 0 then b2(a)=-26000
a=where(b3 gt fmin(5),count)
if count ne 0 then b3(a)=-25000
a=where(b4 gt fmin(6),count)
if count ne 0 then b4(a)=-24000
a=where(b5 gt fmin(7),count)
if count ne 0 then b5(a)=-23000
a=where(b6 gt fmin(8),count)
if count ne 0 then b6(a)=-22000
a=where(b7 gt fmin(9),count)
if count ne 0 then b7(a)=-21000

;so here's a bunch of different combos.  pick one or make your own

redx=red(*,1,0)-red(*,2,0)
redy=red(*,0,0)-red(*,1,0)
redx4=red(*,1,1)-red(*,2,1)
redy4=red(*,0,1)-red(*,1,1)
;plotcol,d,redx,redy,redx4,redy4,h-k,j-h,'H-K','J-H',-.3,2.,-.5,4.2,1,ninc,-.1,-.25,5,Av

redx=red(*,2,0)-red(*,3,0)
redy=red(*,0,0)-red(*,1,0)
redx4=red(*,2,1)-red(*,3,1)
redy4=red(*,0,1)-red(*,1,1)
;plotcol,d,redx,redy,redx4,redy4,k-b1,j-h,'K-[3.6]','J-H',-.5,3.5,0,3.5,0,ninc,-0.,.6,5,Av

redy=red(*,3,0)-red(*,4,0)
redx=red(*,5,0)-red(*,6,0)
redy4=red(*,3,1)-red(*,4,1)
redx4=red(*,5,1)-red(*,6,1)
plotcol,d,redx,redy,redx4,redy4,b3-b4,b1-b2,'[5.8]-[8.0]','[3.6]-[4.5]',-.5,2.0,-.5,4,0,ninc,0,0,60,Av
;note, [3.6]-[4.5] colors go up to 15

redy=red(*,3,0)-red(*,5,0)
redx=red(*,6,0)-red(*,7,0)
redy4=red(*,3,1)-red(*,5,1)
redx4=red(*,6,1)-red(*,7,1)
;plotcol,d,redx,redy,redx4,redy4,b4-b5,b1-b3,'[8.0]-[24]','[3.6]-[5.8]',-1.,10,-1,7,0,ninc,0,0,60,Av

redy=red(*,6,0)-red(*,7,0)
redx=red(*,7,0)-red(*,8,0)
redy4=red(*,6,1)-red(*,7,1)
redx4=red(*,7,1)-red(*,8,1)
;plotcol,d,redx,redy,redx4,redy4,b5-b6,b4-b5,'[24]-[70]','[8.0]-[24]',0.0,15.0,0,10,0,ninc,0,0,180,Av

device,/close
end
