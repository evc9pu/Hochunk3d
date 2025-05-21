pro sedo,dirm,rstar,tstar,iplotstar,lumscale,intit,d,napin,$
        Jy,atfilein,planckin,Av,flux,cs

print,''

loadct,10
ibar=0

nap=0l
nap=napin-1l

;ninc=nincin
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
lumstar=4.*!pi*(rstar*6.955d10)^2*sigma*tstar^2*tstar^2
rstar2=lumstar/(4.d0*!pi*sigma)/tstar^2/tstar^2
print,'stellar luminosity ',lumstar/lsun
print,'lsun',lsun
;lumscale is scale factor to apply source luminosity to, default = 1
fnorm=lumscale/(4.d0*!pi)/d^2/pc/pc*lsun
;fnorm=9.14706d0*lsun/(4.d0*!pi)/d^2/pc/pc
;print,'fnorm',fnorm

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

x4=dblarr(ninc*nfreq)
f4=x4
x4s=x4
x4d=x4
f4s=f4
f4d=f4
if nphi gt 2 then begin
   print,'error, this is only for 2-D models, nphi=2; stopping program'
   return
endif
readcol,dirm+'/flux.dat',x4,f4,skipline=nap*nfreq*ninc+2l,numline=nfreq*ninc,format='(d,d)'
readcol,dirm+'/flux_star.dat',x4s,f4s,skipline=nap*nfreq*ninc+2l,numline=nfreq*ninc,format='(d,d)'
readcol,dirm+'/flux_disk.dat',x4d,f4d,skipline=nap*nfreq*ninc+2l,numline=nfreq*ninc,format='(d,d)'
readcol,dirm+'/flux_enve.dat',x4e,f4e,skipline=nap*nfreq*ninc+2l,numline=nfreq*ninc,format='(d,d)'
print,'max',max(f4)
;read in peeled spectrum
;readcol,dirm+'/eflux.dat',x5,f5,skipline=3+nfreq*(nap),numline=nfreq

f4=f4*fnorm
f4s=f4s*fnorm
f4d=f4d*fnorm
f4e=f4e*fnorm

;f4=f4/lsun

;integrate fnu*dnu=flam*dlam=(lam*flam)*d(ln(lam))
sum=0.
for ii=0,ninc-1 do begin
  n1=ii*nfreq
  n2=n1+nfreq-1
  x4hi=x4(n1:n2)
  if flux eq 'a' then f4hi=f4(n1:n2)
  if flux eq 's' then f4hi=f4s(n1:n2)
  if flux eq 'd' then f4hi=f4d(n1:n2)
  if flux eq 'e' then f4hi=f4e(n1:n2)
  x4hi=alog(x4hi*1.e-4)
  result=int_tabulated(x4hi,f4hi,/sort)
  sum=sum+result
endfor
;print,'integrated luminosity (lsun)',sum/(ninc)
;print,'integrated luminosity (lsun)',sum/(ninc)/fnorm
;print,'compare to value in disk.dat file'
print,'integrated luminosity, no foreground extinction applied'
  if flux eq 'a' then print,'integrated total luminosity (lsun)',sum/(ninc)/fnorm
  if flux eq 's' then print,'integrated stellar luminosity (lsun)',sum/(ninc)/fnorm
  if flux eq 'd' then print,'integrated disk luminosity (lsun)',sum/(ninc)/fnorm
  if flux eq 'e' then print,'integrated envelope luminosity (lsun)',sum/(ninc)/fnorm
x4hi=x4(0:nfreq-1)
f4hi=f4(0:nfreq-1)
x4hi=alog(x4hi*1.e-4)
result=int_tabulated(x4hi,f4hi,/sort)
;print,'edge-on result, integrating lam*flam*dloglam',result/fnorm
x4hi=x4(0:nfreq-1)
f4hi=f4(0:nfreq-1)/x4hi
result=int_tabulated(x4hi,f4hi,/sort)
;print,'edge-on result, integrating flam*dlam',result/fnorm

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
  f4(nfreq*ii:nfreq*ii+nfreq-1)=f4(nfreq*ii:nfreq*ii+nfreq-1)*extinct
  f4s(nfreq*ii:nfreq*ii+nfreq-1)=f4s(nfreq*ii:nfreq*ii+nfreq-1)*extinct
  f4d(nfreq*ii:nfreq*ii+nfreq-1)=f4d(nfreq*ii:nfreq*ii+nfreq-1)*extinct
  f4e(nfreq*ii:nfreq*ii+nfreq-1)=f4e(nfreq*ii:nfreq*ii+nfreq-1)*extinct
endfor

f4all=f4
if flux eq 'e' then f4=f4e
if flux eq 'd' then f4=f4d
if flux eq 's' then f4=f4s

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

maxm=max(f4all,/nan)*2.
;maxm=1.e-9
;minm=maxm*2.e-2
minm=maxm*1.e-5
;print,'max',max(f4all)

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
if iplotstar eq 1 then oplot,lama,nuhnu,linestyle=1
print,'stellar atm peak',max(nuhnu)
for ii=0,ninc-1,1 do begin
  cosi=1.-(float(ii)+0.5)*2./float(ninc)
  theti=acos(cosi)*180./!pi
  print,'plotting viewing angle (cosi,theti)',cosi,theti
  if cosi gt 0 then begin
     ccol=(cmax-cmin)/(ninc/2+1)*ii+cmin
  endif else begin
     ccol=(cmax-cmin)/(ninc/2+1)*(ninc-ii)+cmin
  endelse
  !psym=0
  ;don't plot first frequency since it includes all flux in this and
  ;shortward.
  oplot,(x4(0:nfreq-2)),(f4(nfreq*ii:nfreq*ii+nfreq-2)),linestyle=0,$
    color=ccol,thick=2
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

