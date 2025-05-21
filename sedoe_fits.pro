pro sedoe_fits,dirm,rstar,tstar,iplotstar,lumscale,intit,d,napin,$
         Jy,atfilein,planckin,Av,flux,cs,op,ls,inpcol

print,''

loadct,10
ibar=0

;ninc=nincin
atfile=atfilein
cplanck=planckin
;mstar=mstarin
;mdot=mdotin
;lum=lumin
;nap=0l
;iap=0l
nap=napin-1l
;iap=iapin-1l

;calculate fnorm
print,'rstar,tstar',rstar,tstar
c=2.9979d14
sigma=5.67d-5
pc=3.0857d18
lsun=3.845d33
lumstar=4.d0*!pi*(rstar*6.955d10)^2*sigma*tstar^2*tstar^2
rstar2=lumstar/(4.d0*!pi*sigma)/tstar^2/tstar^2
print,'stellar luminosity calculated from rstar,tstar (Lsun)',lumstar/lsun
;print,'lsun',lsun
fnorm=lumscale/(4.d0*!pi)/d^2/pc/pc*lsun
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

;read in model
filename=dirm+'/peel_hypercube.fits.gz'
fluxarr=mrdfits(filename,0,header)
wavearr=mrdfits(filename,1,headwave)
peelarr=mrdfits(filename,2,headpeel)
aparr=mrdfits(filename,3,headap)

nfreq=fxpar(header,'NAXIS1')
npeel=fxpar(header,'NAXIS2')
napmax=fxpar(header,'NAXIS3')
no=fxpar(header,'NAXIS4')

x5=fltarr(nfreq)
x5=wavearr.wavelength
x5s=x5
x5d=x5
x5e=x5
peeltheta=peelarr.theta
peelphi=peelarr.phi
aps=aparr.aperture
print,'npeel',npeel
print,'peeled theta', peelarr.theta
print,'peeled phi',peelarr.phi
print,'apertures ',aps
f5=fltarr(nfreq,npeel,napmax)
;help,f5
f5s=f5
f5d=f5
f5t=f5
f5=fluxarr(*,*,*,0,0)
f5d=fluxarr(*,*,*,2,0)
f5s=fluxarr(*,*,*,1,0)
f5e=fluxarr(*,*,*,3,0)

;print,'max f5',max(f5,/nan)
;print,x5
f5=f5*fnorm
f5s=f5s*fnorm
f5d=f5d*fnorm
f5e=f5e*fnorm

;integrate fnu*dnu=flam*dlam=(lam*flam)*d(ln(lam))
print,'integrated flux along this direction'
print,'note: not the same as source luminosity since flux varies with viewing angle (so pole-on results are usually higher, edge-on lower)'
sum=0.d0
x5hi=alog(x5*1.e-4)
print,'nap',nap
for it=0,npeel-1 do begin
 if flux eq 'a' then f5hi=f5(0:nfreq-1,it,nap)
 if flux eq 's' then f5hi=f5s(0:nfreq-1,it,nap)
 if flux eq 'd' then f5hi=f5d(0:nfreq-1,it,nap)
 if flux eq 'e' then f5hi=f5e(0:nfreq-1,it,nap)
 result=0
 help,f5hi
 a=where(f5hi gt 0)
 print,'n_elements(a)',n_elements(a)
 if n_elements(a) gt 0 then result=int_tabulated(x5hi,f5hi,/sort)
  sum=sum+result
endfor

print,'integrated luminosity, no foreground extinction applied'
  if flux eq 'a' then print,'integrated total luminosity (lsun)',sum/float(npeel)/fnorm
  if flux eq 's' then print,'integrated star origin luminosity (lsun)',sum/float(npeel)/fnorm
  if flux eq 'd' then print,'integrated disk origin luminosity (lsun)',sum/float(npeel)/fnorm
  if flux eq 'e' then print,'integrated envelope origin luminosity (lsun)',sum/float(npeel)/fnorm


;add in reddening
xnunew=x5
readcol,'../models/parfiles/kmhnew_extrap.par',klam,junk1,junk2,kkap
;readcol,'../models/parfiles/r400_ice095_extrap.par',klam,junk1,junk2,kkap
;interpolate onto flux grid
kapnew=interpol(kkap,klam,xnunew)
kapv=interpol(kkap,klam,0.55)
taul=kapnew/kapv/1.086*Av
extinct=exp(-taul)

;.r sede_fi,f5
for it=0,npeel-1 do begin
  f5(*,it,nap)=f5(*,it,nap)*extinct
  f5s(*,it,nap)=f5s(*,it,nap)*extinct
  f5d(*,it,nap)=f5d(*,it,nap)*extinct
  f5e(*,it,nap)=f5e(*,it,nap)*extinct
endfor

;convert nuhnu and f5 to Jy; units are currently in lamflam=nufnu (erg...)
;first, divide by frequency to get erg/cm2/s/hz
;fnu=lamflam/nu
;c is in microns/s
;x5 is in microns
if Jy eq 1 then begin
  for it=0,npeel-1 do begin
    f5(*,it,nap)=f5(*,it,nap)/c*x5
    f5s(*,it,nap)=f5s(*,it,nap)/c*x5
    f5d(*,it,nap)=f5d(*,it,nap)/c*x5
    f5e(*,it,nap)=f5e(*,it,nap)/c*x5
  endfor
  nuhnu=nuhnu/c*lama
 ; print,'max f5',max(f5)
  ;now convert to mJy
  f5=f5*1.e26
  f5s=f5s*1.e26
  f5d=f5d*1.e26
  f5e=f5e*1.e26
  nuhnu=nuhnu*1.e26
endif

print,'max f5',max(f5)

maxm=max(f5,/nan)*2.
;maxm=1.e-9
;minm=maxm*2.e-2
minm=maxm*1.e-5

if flux eq 'e' then f5=f5e
if flux eq 'd' then f5=f5d
if flux eq 's' then f5=f5s

if jy eq 1 then begin
  ytit='Flux Density (mJy)'
endif else begin
  ytit='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)'
endelse

;write out wavelength and flux to a file
;for it=0,npeel-1 do begin
;  forprint,x5,f5(*,it,nap),textout=3
;endfor

;ii=ninc-1
!psym=0
cmax=240
cmin=10
if op eq 0 then begin
   plot,x5,0.*f5(*,0,0),$
   /xlog,/ylog,$
   yrange=[minm,maxm],$
   xrange=[0.1,5000.],$
   xmargin=[10,3],charthick=2,$
   ;xrange=[alog10(0.2),alog10(5000.)],$
   xstyle=1,ystyle=1,title=intit,$
   linestyle=ls,ytitle=ytit,$
   xtitle='!4k!3(!4l!3m)',charsize=cs,thick=3
endif
if iplotstar eq 1 then oplot,lama,nuhnu,linestyle=1

for it=0,npeel-1 do begin
  ccol=(cmax-cmin)/(npeel/2+1)*it+cmin
  if npeel eq 1 then ccol=inpcol
  print,ccol
  oplot,x5(1:nfreq-2),f5(1:nfreq-2,it,nap),thick=3,color=ccol,linestyle=ls
endfor


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


end

