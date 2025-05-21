pro sed,dirm,nfreq,nincin,rstar,tstar,lumscale,intit,d,napin,$
        Jy,atfilein,planckin,Av,flux,cs

loadct,10
ibar=0

nap=0l
nap=napin-1l

ninc=nincin
atfile=atfilein
cplanck=planckin
;mstar=mstarin
;mdot=mdotin
;lum=lumin

;calculate fnorm
c=2.9979d14
sigma=5.67d-5
pc=3.0857d18
lsun=3.845d33
lumstar=4.d0*!pi*(rstar*6.955d10)^2*sigma*tstar^2*tstar^2
rstar2=lumstar/(4.d0*!pi*sigma)/tstar^2/tstar^2
print,'stellar luminosity ',lumstar/lsun
print,'lsun',lsun
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
print,'stellar atmosphere peak',max(nuhnu)
print,size(nuhnu)
;print,'calculated rstar ',sqrt(rstar2)/6.96e10

;read in model, skip to aperture
x4=dblarr(ninc*nfreq)
f4=x4
x4s=x4
x4d=x4
f4s=f4
f4d=f4
readcol,dirm+'/flux.dat',x4,f4,skipline=nap*nfreq*ninc+1l,numline=nfreq*ninc,format='(d,d)'
readcol,dirm+'/sflux.dat',x4s,f4s,skipline=nap*nfreq*ninc+1l,numline=nfreq*ninc,format='(d,d)'
readcol,dirm+'/dflux.dat',x4d,f4d,skipline=0,numline=nfreq*ninc,format='(d,d)'
print,'max',max(f4,/nan)
x4t=x4
f4t=f4
;read in peeled spectrum
;readcol,dirm+'/eflux.dat',x5,f5,skipline=3+nfreq*(nap),numline=nfreq

f4=f4*fnorm
f4s=f4s*fnorm
f4d=f4d*fnorm

;f4=f4/lsun

;integrate fnu*dnu=flam*dlam=(lam*flam)*d(ln(lam))
sum=0.
for ii=0,ninc-1 do begin
  n1=ii*nfreq
  n2=n1+nfreq-1
  x4hi=x4(n1:n2)
  f4hi=f4(n1:n2)
  x4hi=alog(x4hi*1.e-4)
  result=int_tabulated(x4hi,f4hi,/sort)
  sum=sum+result
endfor
;print,'integrated luminosity (lsun)',sum/(ninc)
print,'integrated luminosity (lsun)',sum/(ninc)/fnorm
print,'compare to value in disk.dat file'
x4hi=x4(0:nfreq-1)
f4hi=f4(0:nfreq-1)
x4hi=alog(x4hi*1.e-4)
result=int_tabulated(x4hi,f4hi,/sort)
print,'edge-on result, integrating lam*flam*dloglam',result/fnorm
x4hi=x4(0:nfreq-1)
f4hi=f4(0:nfreq-1)/x4hi
result=int_tabulated(x4hi,f4hi,/sort)
print,'edge-on result, integrating flam*dlam',result/fnorm

;f4=f4*lsun
  
;add in reddening
xnunew=x4(0:nfreq-1)
readcol,'../models/parfiles/kmh.par',klam,junk1,junk2,kkap
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
  f4(nfreq*ii:nfreq*ii+nfreq-1)=f4(nfreq*ii:nfreq*ii+nfreq-1)*extinct
  f4s(nfreq*ii:nfreq*ii+nfreq-1)=f4s(nfreq*ii:nfreq*ii+nfreq-1)*extinct
  f4d(nfreq*ii:nfreq*ii+nfreq-1)=f4d(nfreq*ii:nfreq*ii+nfreq-1)*extinct
  f4t(nfreq*ii:nfreq*ii+nfreq-1)=f4t(nfreq*ii:nfreq*ii+nfreq-1)*extinct
endfor

;barb test!  delete this!
;kapnew=interpol(kkap,klam,lama)
;kapv=interpol(kkap,klam,0.55)
;taul=kapnew/kapv/1.086*Av
;extinct=exp(-taul)
;nuhnu=nuhnu*extinct

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
  ;now convert to mJy
  f4=f4*1.e26
  f4all=f4all*1.e26
;  f5=f5*1.e26
  nuhnu=nuhnu*1.e26
endif

maxm=max(f4all,/nan)
;maxm=1.e-9
;minm=maxm*2.e-2
minm=maxm*1.e-5
print,'max',max(f4all,/nan)

if jy eq 1 then begin
  ytit='Flux Density (mJy)'
endif else begin
  ytit='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)'
endelse

ii=ninc-1
!psym=0
cmax=240
cmin=10
plot,(x4(0:nfreq-1)),0.*(f4(nfreq*ii:nfreq*ii+nfreq-1)),$
/xlog,/ylog,$
yrange=[minm,maxm],$
;xrange=[2.,100.],$
xrange=[0.1,5000.],$
xmargin=[10,3],charthick=2,$
xstyle=1,ystyle=1,title=intit,$
linestyle=0,ytitle=ytit,$
xtitle='!4k!3(!4l!3m)',charsize=cs,thick=2
oplot,lama,nuhnu,linestyle=1
print,'stellar atm peak',max(nuhnu)
for ii=0,ninc-1,1 do begin
  cosi=(float(ii)+0.5)/10.
  theti=acos(cosi)*180./!pi
  print,'plotting viewing angle',cosi,theti
  ccol=(cmax-cmin)/11*ii+cmin
  !psym=0
  ;don't plot first frequency since it includes all flux in this and
  ;shortward.
  oplot,(x4(0:nfreq-2)),(f4(nfreq*ii:nfreq*ii+nfreq-2)),linestyle=0,color=ccol,$
   thick=2
endfor
;oplot,x5,f5,thick=5

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

end

