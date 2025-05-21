pro sed,dirm,rstar,tstar,iplotstar,lumscale,intit,d,napin,$
  Jy,atfilein,planckin,Av,flux,cs,op,i1,i2,ls,ifilt

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

;read in model, skip to aperture
print,'model directory',dirm
readcol,dirm+'/flux.dat',nfreqin,nphiin,nincin,napmaxin,numline=1,format='(i,i,i,i)'
napmax=long(napmaxin(0))
nfreq=long(nfreqin(0))
nphi=long(nphiin(0))
ninc=long(nincin(0))

print,'nfreq,nphi,ninc,napmax',nfreq,nphi,ninc,napmax
;end
if nap gt napmax+1 then begin
  print,'error  nap > napmax; stopping program'
  return
endif
x4=dblarr(ninc*nfreq*nphi)  ; shouldn't that be ninc*nfreq*nphi?????  how did this work before???
f4=x4
x4s=x4
x4d=x4
f4s=f4
f4d=f4
f4t=f4
skipline=0l
;if nphi gt 2 then begin
;   print,'error, this is only for 2-D models, nphi=1; stopping program'
;   return
;endif
readcol,dirm+'/flux.dat',x4,f4,skipline=nap*nfreq*ninc*nphi+2l,numline=nfreq*ninc*nphi,format='(d,d)'
readcol,dirm+'/flux_ther.dat',x4t,f4t,skipline=nap*nfreq*ninc*nphi+2l,numline=nfreq*ninc*nphi,format='(d,d)'
readcol,dirm+'/flux_scat.dat',x4s,f4s,skipline=nap*nfreq*ninc*nphi+2l,numline=nfreq*ninc*nphi,format='(d,d)'
readcol,dirm+'/flux_dire.dat',x4d,f4d,skipline=nap*nfreq*ninc*nphi+2l,numline=nfreq*ninc*nphi,format='(d,d)'
;readcol,dirm+'/flux_dire.dat',x4d,f4d,skipline=0,numline=nfreq*ninc*nphi,format='(d,d)'
;print,'max',max(f4,/nan)
x4t=x4
f4t=f4
;read in peeled spectrum
;readcol,dirm+'/eflux.dat',x5,f5,skipline=3+nfreq*(nap),numline=nfreq

f4=f4*fnorm
f4s=f4s*fnorm
f4t=f4t*fnorm
f4d=f4d*fnorm

;f4=f4/lsun

;integrate fnu*dnu=flam*dlam=(lam*flam)*d(ln(lam))
sum=0.
for ii=0,ninc*nphi-1 do begin
  n1=ii*nfreq
  n2=n1+nfreq-1
  x4hi=x4(n1:n2)
  f4hi=f4(n1:n2)
  if flux eq 't' then f4hi=(f4(n1:n2)-f4s(n1:n2)-f4d(n1:n2))
;  if flux eq 't' then f4hi=f4t(n1:n2)
  if flux eq 'd' then f4hi=f4d(n1:n2)
  if flux eq 's' then f4hi=f4s(n1:n2)
  x4hi=alog(x4hi*1.e-4)
  result=int_tabulated(x4hi,f4hi,/sort)
  sum=sum+result
endfor
;print,'integrated luminosity (lsun)',sum/(ninc)
;print,'integrated luminosity (lsun)',sum/(ninc)/fnorm
print,'integrated luminosity, no foreground extinction applied'
  if flux eq 'a' then print,'integrated total luminosity (lsun)',sum/float(ninc*nphi)/fnorm
  if flux eq 't' then print,'integrated thermal luminosity (lsun)',sum/float(ninc*nphi)/fnorm
  if flux eq 'd' then print,'integrated direct stellar luminosity (lsun)',sum/float(ninc*nphi)/fnorm
  if flux eq 's' then print,'integrated scattered luminosity (lsun)',sum/float(ninc*nphi)/fnorm
;print,'compare to value in disk.dat file'
;x4hi=x4(0:nfreq-1)
;f4hi=f4(0:nfreq-1)
;x4hi=alog(x4hi*1.e-4)
;result=int_tabulated(x4hi,f4hi,/sort)
;print,'edge-on result, integrating lam*flam*dloglam',result/fnorm
;x4hi=x4(0:nfreq-1)
;f4hi=f4(0:nfreq-1)/x4hi
;result=int_tabulated(x4hi,f4hi,/sort)
;print,'edge-on result, integrating flam*dlam',result/fnorm

;f4=f4*lsun
  
;add in reddening
xnunew=x4(0:nfreq-1)
readcol,'../models/parfiles/kmhnew_extrap.par',klam,junk1,junk2,kkap
;readcol,'../models/parfiles/r400_ice095.par',klam,junk1,junk2,kkap
;interpolate onto flux grid
kapnew=interpol(kkap,klam,xnunew)
kapv=interpol(kkap,klam,0.55)
taul=kapnew/kapv/1.086*Av
extinct=exp(-taul)
for ii=0,ninc-1 do begin
; subtract direct and scattered spectrum (to get thermal only)
  f4t(nfreq*ii:nfreq*ii+nfreq-1)=f4(nfreq*ii:nfreq*ii+nfreq-1)-$
     f4s(nfreq*ii:nfreq*ii+nfreq-1)-f4d(nfreq*ii:nfreq*ii+nfreq-1)
;  f4t(nfreq*ii:nfreq*ii+nfreq-1)=f4t(nfreq*ii:nfreq*ii+nfreq-1)*extinct
  f4(nfreq*ii:nfreq*ii+nfreq-1)=f4(nfreq*ii:nfreq*ii+nfreq-1)*extinct
  f4s(nfreq*ii:nfreq*ii+nfreq-1)=f4s(nfreq*ii:nfreq*ii+nfreq-1)*extinct
  f4d(nfreq*ii:nfreq*ii+nfreq-1)=f4d(nfreq*ii:nfreq*ii+nfreq-1)*extinct
endfor

f4all=f4
if flux eq 't' then f4=f4t
if flux eq 'd' then f4=f4d
if flux eq 's' then f4=f4s

;f4all=f4all*fnorm
;f5=f5*fnorm

;convert nuhnu and f4 to Jy; units are currently in lamflam=nufnu (erg...)
;first, divide by frequency to get erg/cm2/s/hz
;fnu=lamflam/nu
;c is in microns/s
;x4 is in microns
if Jy eq 1 then begin
  f4=f4/c*x4
  f4all=f4all/c*x4
;  f5=f5/c*x5
  nuhnu=nuhnu/c*lama
  ;now convert to Jy
  f4=f4*1.e23
  f4all=f4all*1.e23
;  f5=f5*1.e23
  nuhnu=nuhnu*1.e23
endif

maxm=max(f4all,/nan)*2.
;maxm=1.e-6
;minm=maxm*2.e-2
minm=maxm*1.e-5
;print,'max',max(f4all,/nan)

if jy eq 1 then begin
  ytit='Flux Density (Jy)'
endif else begin
  ytit='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)'
endelse


ii=ninc-1
!psym=0
cmax=240
cmin=10
if op eq 0 then begin
   plot,(x4(0:nfreq-1)),0.*(f4(nfreq*ii:nfreq*ii+nfreq-1)),$
        /xlog,/ylog,$
        yrange=[minm,maxm],$
;xrange=[2.,100.],$
        xrange=[0.1,2000.],$
        xmargin=[10,3],charthick=2,$
        xstyle=1,ystyle=1,title=intit,$
        linestyle=0,ytitle=ytit,$
        xtitle='!4k!3(!4l!3m)',charsize=cs,thick=2
endif
if iplotstar eq 1 then oplot,lama,nuhnu,linestyle=1
;print,'stellar atm peak',max(nuhnu)
for ii=i1,i2 do begin
  if nphi eq 1 then begin
    cosi=1.-(float(ii)+0.5)*2./float(ninc)
    theti=acos(cosi)*180./!pi
    print,'plotting viewing angle (cosi,theti)',cosi,theti
    if cosi gt 0 then begin
       ccol=(cmax-cmin)/(ninc/2+1)*ii+cmin
    endif else begin
       ccol=(cmax-cmin)/(ninc/2+1)*(ninc-ii)+cmin
    endelse
  endif else begin
     ccol=(cmax-cmin)/(ninc*nphi+1)+cmin
  endelse
  ;ccol=0   ;for figure 1; comment this out when done.
  !psym=0
  ;don't plot first frequency since it includes all flux in this and
  ;shortward.
  oplot,(x4(0:nfreq-2)),(f4(nfreq*ii:nfreq*ii+nfreq-2)),linestyle=ls,color=ccol,$
   thick=2
  ;calculate herschel fluxes real quick
  xarr=x4(0:nfreq-2)
  yarr=f4(nfreq*ii:nfreq*ii+nfreq-2)
  a=where(xarr le 102.)
  i100=a(0)
;  print,'100',yarr(i100)
  a=where(xarr le 165.)
  i160=a(0)
;  print,'160',yarr(i160)
  a=where(xarr le 255.)
  i250=a(0)
;  print,'250',yarr(i250)  
  a=where(xarr le 370.)
  i350=a(0)
;  print,'350',yarr(i350)  
  a=where(xarr le 501.)
  i500=a(0)
;  print,'500',yarr(i500)  
;  print,'log10(SPIRE250/PACS160),log10(PACS160)',alog10(yarr(i250)/yarr(i160)),alog10(yarr(i160))
endfor
;print,xarr
;oplot,x5,f5,thick=5
;oplot,[160,160],[0.001,10000.],col=155
;oplot,[250,250],[0.001,10000.],col=155

maserw=[1.2, 1.6, 2.2, 3.6, 4.5, 5.8, 8.0, 24, 70, 100, 160, 160, 250, 350, 500]
maserf=[9.2409510e-04,  4.8180381e-04,  7.6769583e-04,  3.5231863e-03,  3.6814220e-03,  1.4992598e-02, $
3.6091201e-02,  1.4586538e-01,  1.4032848e+00,  2.1700363e+00, -3.1171774e+01,  1.6200641e+00, $
6.0855746e-01,  2.3864444e-01,  5.3819758e-01]*1000.

;oplot,maserw,maserf,psym=5

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
;oplotdata,'3'

;print,'fnorm',fnorm

; convolve with filter functions
if ifilt eq 1 then begin
  close,1
  openw,1,'filterfluxes.txt'
  printf,1,'filter fluxes in Jy at 2mass, Spitzer, and herschel bands'
  printf,1,'cosi,J, H, K, 3.6, 4.5, 5.8, 8, 24, 70, 160'
  ;loop over inclination
  for ii=0,ninc-1 do begin
     n1=ii*nfreq
     n2=n1+nfreq-1
     f=f4(n1:n2)
     lam=x4(n1:n2)
     cosi=1.-(float(ii)+0.5)*2./float(ninc)
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
    printf,1,format='(12e12.4)',cosi,flux1,flux2,flux3,flux4,flux5,flux6,flux7,flux8,flux9,flux10
  endfor
  close,1
endif

;test of filter functions on modcII;  works!
;oplot,[1.2,1.6,2.2,3.6,4.5,5.8,8.0,24.0,70.0,160.],[4.7927e+01,  6.5302e+01,  5.0406e+01,  3.2604e+01,  2.7353e+01,  2.5276e+01,  3.3884e+01,  1.3471e+02,  5.0827e+02,  4.2945e+02],symsize=1,psym=5


end

