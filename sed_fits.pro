pro sed_fits,dirm,rstar,tstar,iplotstar,lumscale,intit,d,napin,$
  Jy,atfilein,planckin,Av,flux,cs,op,i1,i2,ip1,ip2,ls,ifilt,isphav

;print,dirm,rstar,tstar

print,''

if Jy eq 0 and ifilt eq 1 then begin
  print,'ERROR, you specified Jy=0 and ifilt=1'
  print,'if ifilt=1, you need Jy=1; if ifilt=0, Jy can be 0 or 1'
  print,'change this and then run again!'
  stop
endif

loadct,10
ibar=0

nap=0l
nap=napin-1l
print,'nap',nap

atfile=atfilein
cplanck=planckin

;calculate fnorm
c=2.9979d14
sigma=5.67d-5
pc=3.0857d18
lsun=3.845d33
lumstar=4.d0*!pi*(rstar*6.955d10)^2*sigma*tstar^2*tstar^2
rstar2=lumstar/(4.d0*!pi*sigma)/tstar^2/tstar^2
print,'stellar luminosity calculated from rstar,tstar ',lumstar/lsun
;print,'lsun',lsun
;lumscale is scale factor to apply source luminosity to, default = 1
fnorm=lumscale/(4.d0*!pi)/d^2/pc/pc*lsun
;fnorm=9.14706d0*lsun/(4.d0*!pi)/d^2/pc/pc
print,'fnorm',fnorm

;atmnorm=4.d0*!pi*rstar2/d/d/pc/pc  ;old kurucz models
atmnorm=rstar2/d/d/pc/pc     ;tom's atmosphere models and idlastro planck 
if cplanck eq 'no' or cplanck eq 'NO' then begin
  ;read in stellar spectrum used as input for comparison
; old format
;  readcol,atfile,lama,junk,hnu,skipline=3
  ; convert from nm
;  lama=lama*1.d-3
; new format
  readcol,atfile,lama,hnu,skipline=3
  nuhnu=hnu/lama*2.9979d14
  ;f=hnu*rstar^2/(d^2)
  ;solve for rstar using lum (in case we want to use lum=1lsun)
  nuhnu=nuhnu*atmnorm
endif else begin
  lama=(findgen(250))/float(250)*4.5
  lama=10^lama/10.
  nuhnu=planck(lama*10000.,tstar)*lama*10000.   ;wavelength units are A
  nuhnu=nuhnu*atmnorm
endelse
;print,'stellar atmosphere peak',max(nuhnu)
;print,size(nuhnu)
;print,'calculated rstar ',sqrt(rstar2)/6.96e10
a=where(nuhnu gt 0)
lama=lama(a)
nuhnu=nuhnu(a)

;read in model
filename=dirm+'/flux_hypercube.fits.gz'
fluxarr=mrdfits(filename,0,header)
wavearr=mrdfits(filename,1,headwave)

nfreq=fxpar(header,'NAXIS1')
nthet=fxpar(header,'NAXIS2')
nphi=fxpar(header,'NAXIS3')
napmax=fxpar(header,'NAXIS4')
no=fxpar(header,'NAXIS5')

print,'nfreq,nthet,nphi',nfreq,nthet,nphi

x4=fltarr(nfreq)
x4=wavearr.wavelength
x4s=x4
x4t=x4
f4=fltarr(nfreq,nthet,nphi,napmax)
f4s=f4
f4d=f4
f4t=f4
f4=fluxarr(*,*,*,*,0,0)
f4d=fluxarr(*,*,*,*,5,0)
f4s=fluxarr(*,*,*,*,6,0)
f4t=fluxarr(*,*,*,*,7,0)

f4=f4*fnorm
f4s=f4s*fnorm
f4t=f4t*fnorm
f4d=f4d*fnorm

;f4=f4/lsun

;integrate fnu*dnu=flam*dlam=(lam*flam)*d(ln(lam))
sum=0.
for it=0,nthet-1 do begin
for ip=0,nphi-1 do begin
   if flux eq 'a' then f4hi=f4(0:nfreq-1,it,ip,nap)
   if flux eq 't' then f4hi=f4t(0:nfreq-1,it,ip,nap)
   if flux eq 's' then f4hi=f4s(0:nfreq-1,it,ip,nap)
   if flux eq 'd' then f4hi=f4d(0:nfreq-1,it,ip,nap)
   x4hi=alog(x4*1.e-4)
  result=int_tabulated(x4hi,f4hi,/sort)
  sum=sum+result
endfor
endfor
print,'integrated luminosity, no foreground extinction applied'
  if flux eq 'a' then print,'integrated total luminosity (lsun)',sum/float(nthet*nphi)/fnorm
  if flux eq 't' then print,'integrated thermal luminosity (lsun)',sum/float(nthet*nphi)/fnorm
  if flux eq 'd' then print,'integrated direct stellar luminosity (lsun)',sum/float(nthet*nphi)/fnorm
  if flux eq 's' then print,'integrated scattered luminosity (lsun)',sum/float(nthet*nphi)/fnorm

;integrate stellar atmosphere file
;result2=int_tabulated(alog(lama*1.e-4),nuhnu,/sort)
;print,'integrated stellar atmosphere file (lsun)',result2/fnorm  

;print,'MC output/stellaratm',result/result2

;add in reddening
xnunew=x4(0:nfreq-1)
readcol,'../models/parfiles/kmhnew_extrap.par',klam,junk1,junk2,kkap
;readcol,'../models/parfiles/r400_ice095.par',klam,junk1,junk2,kkap
;interpolate onto flux grid
kapnew=interpol(kkap,klam,xnunew)
kapv=interpol(kkap,klam,0.55)
taul=kapnew/kapv/1.086*Av
extinct=exp(-taul)
for it=0,nthet-1 do begin
for ip=0,nphi-1 do begin
  f4(*,it,ip,nap)=f4(*,it,ip,nap)*extinct
  f4t(*,it,ip,nap)=f4t(*,it,ip,nap)*extinct
  f4s(*,it,ip,nap)=f4s(*,it,ip,nap)*extinct
  f4d(*,it,ip,nap)=f4d(*,it,ip,nap)*extinct
endfor
endfor

f4all=f4
if flux eq 't' then f4=f4t
if flux eq 'd' then f4=f4d
if flux eq 's' then f4=f4s

;f4all=f4all*fnorm
;f5=f5*fnorm

print,'max',max(f4all,/nan)

if 0 then begin
;testing for tottle. integrate flux and stellar atmosphere file over specific frequency range
;only for 1-D model
xmin=1.
xmax=30.
f1=interpol(f4,x4,xmin)
f2=interpol(f4,x4,xmax)
a=where(x4 lt xmin)  ;backwards array
a1=a(0)-1
;print,'yo',a1,x4(a1)
a=where(x4 lt xmax)
a2=a(0)
;print,'yo',a2,x4(a2)
n=a1-a2+3
fnew=fltarr(n)
fnew(0)=f2
fnew(n-1)=f1
xnew=fltarr(n)
xnew(0)=xmax
xnew(n-1)=xmin
fnew(1:n-2)=f4(a2:a1)
xnew(1:n-2)=x4(a2:a1)
;print,xnew
;print,fnew
result1=int_tabulated(alog10(xnew*1.e-4),fnew,/sort)
print,'flux integral from',xmin,' to',xmax,result1/fnorm

;not stellar atmosphere file
f1=interpol(nuhnu,lama,xmin)
f2=interpol(nuhnu,lama,xmax)
a=where(lama gt xmin)  
a1=a(0)
;print,'yo',a1,lama(a1)
a=where(lama gt xmax)
a2=a(0)-1
;print,'yo',a2,lama(a2)
n=a2-a1+3
fnew=fltarr(n)
fnew(0)=f1
fnew(n-1)=f2
xnew=fltarr(n)
xnew(0)=xmin
xnew(n-1)=xmax
fnew(1:n-2)=nuhnu(a1:a2)
xnew(1:n-2)=lama(a1:a2)
;print,xnew
;print,fnew
result2=int_tabulated(alog10(xnew*1.e-4),fnew,/sort)
print,'atm integral from',xmin,' to',xmax,result2/fnorm

print,'flux/atm ',result1/result2
endif

;convert nuhnu and f4 to Jy; units are currently in lamflam=nufnu (erg...)
;first, divide by frequency to get erg/cm2/s/hz
;fnu=lamflam/nu
;then multiply by 1.d23 to get Jy
;c is in microns/s
;x4 is in microns
;help,f4
;help,x4
;help,c
if Jy eq 1 then begin
   for it=0,nthet-1 do begin
  for ip=0,nphi-1 do begin
    f4(*,it,ip,nap)=f4(*,it,ip,nap)*1.d23/c*x4
    f4all(*,it,ip,nap)=f4all(*,it,ip,nap)*1.d23/c*x4
  endfor
  endfor
  nuhnu=nuhnu*1.d23/c*lama
endif

print,'max x4',max(x4)

;maxm=max(nuhnu,/nan)*100.
maxm=max(f4all,/nan)*2.
;maxm=100.
;maxm=1.e-10
;maxm=1.e-6
;minm=maxm*.5e-1
 minm=maxm*1.e-5
print,'max',max(f4all,/nan)
;maxm=maxm*.1
;minm=1.e-1*maxm

if jy eq 1 then begin
  ytit='Flux Density (Jy)'
endif else begin
  ytit='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)'
endelse

!psym=0
cmax=240
cmin=10
if op eq 0 then begin
   plot,x4,0.*f4(*,0,0,0),$
        /xlog,/ylog,$
        yrange=[minm,maxm],$
        xrange=[0.1,2000.],$
  ;      xrange=[.8,2],$
        xmargin=[10,3],charthick=2,$
        xstyle=1,ystyle=1,title=intit,$
        linestyle=0,ytitle=ytit,$
        xtitle='!4k!3(!4l!3m)',charsize=cs,thick=2
endif
;nuhnu=smooth(nuhnu,3)
if iplotstar eq 1 then oplot,lama,nuhnu,linestyle=1
;print,'stellar atm peak',max(nuhnu)

if isphav eq 0 then begin
;plot all angles
print,'i1,i2',i1,i2
for it=i1,i2 do begin
for ip=ip1,ip2 do begin
   cosi=1.-(float(it)+0.5)*2./float(nthet)
    theti=acos(cosi)*180./!pi
    print,'plotting viewing angle (cosi,theti)',cosi,theti
    if cosi gt 0 then begin
       ccol=(cmax-cmin)/(nthet/2+1)*it+cmin
    endif else begin
       ccol=(cmax-cmin)/(nthet/2+1)*(nthet-it)+cmin
    endelse
  if i1 eq 0 and i2 eq 0 then ccol=0
  ;ccol=200   ;for figure 1; comment this out when done
  ;ccol=0 ; for figs 25 and 30
  ;don't plot first frequency since it includes all flux in this and shortward.
  oplot,x4(0:nfreq-2),f4(0:nfreq-2,it,ip,nap),linestyle=ls,color=ccol,thick=2
endfor
endfor

endif else begin
;spherically average SED
fphi=total(f4,3,/nan)/float(nphi)
fth=total(fphi,2,/nan)/float(nthet)
oplot,x4(0:nfreq-2),fth(0:nfreq-1,nap),linestyle=ls
endelse

if ibar eq 1 then begin
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A), /FILL
x=10^(findgen(100)*.007)-.5
y=.6
for i=0,99 do begin
  ccol=(cmax-cmin)/100.*i+cmin
  plots,x(i),y,color=ccol,psym=8,symsize=.5
endfor
yc2=.1
xyouts,x(0)-.2,yc2,'i=90',charsize=.6
xyouts,x(99)-.5,yc2,'0',charsize=.6
endif

; convolve with filter functions
if ifilt eq 1 then begin
  close,1
  openw,1,'filterfluxes.txt'
  printf,1,'filter fluxes in Jy at 2mass, Spitzer, and herschel bands'
  printf,1,'cosi, phi, J, H, K, 3.6, 4.5, 5.8, 8, 24, 70, 160'
  ;loop over inclination
  for it=i1,i2 do begin
  cosi=1.-(float(it)+0.5)*2./float(nthet)
  for ip=ip1,ip2 do begin
    phi=360.*(float(ip)+0.5)/float(nphi)
    f=f4(0:nfreq-2,it,ip,nap)
    lam=x4(0:nfreq-2)
  ; integrate over filter functions
    convfilt_jy,'../models/parfiles/2J.txt',1,lam,f,flux1,0
    convfilt_jy,'../models/parfiles/2H.txt',1,lam,f,flux2,0
    convfilt_jy,'../models/parfiles/2K.txt',1,lam,f,flux3,0
    convfilt_jy,'../models/parfiles/I1.txt',1,lam,f,flux4,0
    convfilt_jy,'../models/parfiles/I2.txt',1,lam,f,flux5,0
    convfilt_jy,'../models/parfiles/I3.txt',1,lam,f,flux6,0
    convfilt_jy,'../models/parfiles/I4.txt',1,lam,f,flux7,0
    convfilt_jy,'../models/parfiles/M1.txt',1,lam,f,flux8,0
    convfilt_jy,'../models/parfiles/M2.txt',1,lam,f,flux9,0
    convfilt_jy,'../models/parfiles/M3.txt',1,lam,f,flux10,0
    printf,1,format='(12e12.4)',cosi,phi,flux1,flux2,flux3,flux4,flux5,flux6,flux7,flux8,flux9,flux10
  endfor
  endfor
  close,1
endif

;oplot,[3.4,3.4],[1.e-30,1.e30],linestyle=2
;oplot,[4.5,4.5],[1.e-30,1.e30],linestyle=2
;oplot,[5.8,5.8],[1.e-30,1.e30],linestyle=2

;test of filter functions on modcII;  works!
;oplot,[1.2,1.6,2.2,3.6,4.5,5.8,8.0,24.0,70.0,160.],[4.7927e+01,  6.5302e+01,  5.0406e+01,  3.2604e+01,  2.7353e+01,  2.5276e+01,  3.3884e+01,  1.3471e+02,  5.0827e+02,  4.2945e+02],symsize=1,psym=5

;some stuff for marta.  
imarta=0
if imarta eq 1 then begin
;flux at 70, 160 um
for it=i1,i2 do begin
for ip=ip1,ip2 do begin
fhi=f4(0:nfreq-1,it,ip,nap)
;fhism=smooth(fhi,3)
;fhi=fhism
;help,fhi
;help,x4
;for i=1,nfreq-1 do begin
;print,x4(i),fhi(i)
;endfor
;f3_6=interpol(fhi,x4,3.33)
f3_6=interpol(fhi,x4,3.6)
f4_5=interpol(fhi,x4,4.5)
f5_8=interpol(fhi,x4,5.8)
f70=interpol(fhi,x4,70.)
f160=interpol(fhi,x4,160.)
;print,'f4.5,f70,f160,f4.5/f70',f4_5,f70,f160,f4_5/f70
;print,f4_5,f70,f160,f4_5/f70
zerop=[280.9,179.7,115.0]
flux=[f3_6,f4_5,f5_8]
mags=-2.5*alog10(flux/zerop/1000.)
;print,'3.6-4.5,4.5-5.8',mags(0)-mags(1),mags(1)-mags(2)
print,mags(0)-mags(1),mags(1)-mags(2)
endfor
endfor

print,'phi average'
;average over phi
fphi=total(f4,3,/nan)
for it=i1,i2 do begin
fhi=fphi(0:nfreq-1,it,nap)
f3_6=interpol(fhi,x4,3.6)
f4_5=interpol(fhi,x4,4.5)
f5_8=interpol(fhi,x4,5.8)
zerop=[280.9,179.7,115.0]
flux=[f3_6,f4_5,f5_8]
mags=-2.5*alog10(flux/zerop/1000.)
;print,'3.6-4.5,4.5-5.8',mags(0)-mags(1),mags(1)-mags(2)
print,mags(0)-mags(1),mags(1)-mags(2)
endfor

print,'theta average'
fphi=total(f4,3,/nan)
fth=total(fphi,2,/nan)
fhi=fth(0:nfreq-1,nap)
f3_6=interpol(fhi,x4,3.6)
f4_5=interpol(fhi,x4,4.5)
f5_8=interpol(fhi,x4,5.8)
zerop=[280.9,179.7,115.0]
flux=[f3_6,f4_5,f5_8]
mags=-2.5*alog10(flux/zerop/1000.)
;print,'3.6-4.5,4.5-5.8',mags(0)-mags(1),mags(1)-mags(2)
print,mags(0)-mags(1),mags(1)-mags(2)


print,''

endif

;wave=[0.3543,0.477,0.6231,0.7625,0.9134,1.235,1.662,2.159,3.35,3.6,4.5,4.6,5.8,8,9,11.6,18,22.1,24,70,170]  ;;um
;flux_density=[0.005,0.065,1.393,1.832,5.054,11.379,25.984,27.185,24.139,22.938,44.315,70.500,93.907,287.780,451.500,900.125,2117.000,3375.686,2665.006,1911.180,1664.934]  ;;mJy
;flux_density=flux_density*1.e-3  ;;Jy
;oplot,wave,flux_density,psym=4,thick=5 ;;um & Jy
;print,max(flux_density)

;wave_luciana=[.55,1.2,1.7,2.0,3.7,4.4,5.8,8.0,11.0,22.0,22.0,70.0,430.0,830.0,1300.0] ;;um, Luciana                                                                                                                                                                 
;freq_luciana=c/wave_luciana
;flux_density_luciana1=[3.7e-13,3e-11,5e-11,4e-11,2e-11,3e-11,5e-11,1.2e-10,3e-10,5e-10,3e-10,3e-10,4e-12,1.2e-13,8e-14] ;;ergs s^-1, cm^-2, Luciana                                                                                                                 
;flux_density_luciana=(flux_density_luciana1/freq_luciana)*1.e23 ;;Jy                                                                                                                                                                                          

;oplot,wave_luciana,flux_density_luciana,psym=4,thick=5

;print,'fnorm',fnorm

end

