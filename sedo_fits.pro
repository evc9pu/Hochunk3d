pro sedo_fits,dirm,rstar,tstar,iplotstar,lumscale,intit,d,napin,$
        Jy,atfilein,planckin,Av,flux,cs,it1,it2,ip1,ip2,isphav

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

;read in model
filename=dirm+'/flux_hypercube.fits.gz'
fluxarr=mrdfits(filename,0,header)
wavearr=mrdfits(filename,1,headwave)

nfreq=fxpar(header,'NAXIS1')
ninc=fxpar(header,'NAXIS2')
nphi=fxpar(header,'NAXIS3')
napmax=fxpar(header,'NAXIS4')
no=fxpar(header,'NAXIS5')

x4=fltarr(nfreq)
x4=wavearr.wavelength
x4s=x4
x4d=x4
x4e=x4
f4=fltarr(nfreq,ninc,nphi,napmax)
f4s=f4
f4d=f4
f4e=f4
f4=fluxarr(*,*,*,*,0,0)
f4d=fluxarr(*,*,*,*,2,0)
f4s=fluxarr(*,*,*,*,1,0)
f4e=fluxarr(*,*,*,*,3,0)

f4=f4*fnorm
f4s=f4s*fnorm
f4d=f4d*fnorm
f4e=f4e*fnorm

;integrate fnu*dnu=flam*dlam=(lam*flam)*d(ln(lam))
sum=0.
for it=0,ninc-1 do begin
for ip=0,nphi-1 do begin
   if flux eq 'a' then f4hi=f4(0:nfreq-1,it,ip,nap)
   if flux eq 's' then f4hi=f4s(0:nfreq-1,it,ip,nap)
   if flux eq 'd' then f4hi=f4d(0:nfreq-1,it,ip,nap)
   if flux eq 'e' then f4hi=f4e(0:nfreq-1,it,ip,nap)
   x4hi=alog(x4*1.e-4)
  result=int_tabulated(x4hi,f4hi,/sort)
  sum=sum+result
endfor
endfor
print,'integrated luminosity, no foreground extinction applied'
  if flux eq 'a' then print,'integrated total luminosity (lsun)',sum/float(ninc*nphi)/fnorm
  if flux eq 's' then print,'integrated star origin luminosity (lsun)',sum/float(ninc*nphi)/fnorm
  if flux eq 'd' then print,'integrated disk origin luminosity (lsun)',sum/float(ninc*nphi)/fnorm
  if flux eq 'e' then print,'integrated envelope origin luminosity (lsun)',sum/float(ninc*nphi)/fnorm

;add in reddening
xnunew=x4(0:nfreq-1)
readcol,'../models/parfiles/kmhnew_extrap.par',klam,junk1,junk2,kkap
;readcol,'../models/parfiles/r400_ice095.par',klam,junk1,junk2,kkap
;interpolate onto flux grid
kapnew=interpol(kkap,klam,xnunew)
kapv=interpol(kkap,klam,0.55)
taul=kapnew/kapv/1.086*Av
extinct=exp(-taul)
for it=0,ninc-1 do begin
for ip=0,nphi-1 do begin
  f4(*,it,ip,nap)=f4(*,it,ip,nap)*extinct
  f4s(*,it,ip,nap)=f4s(*,it,ip,nap)*extinct
  f4d(*,it,ip,nap)=f4d(*,it,ip,nap)*extinct
  f4e(*,it,ip,nap)=f4e(*,it,ip,nap)*extinct
endfor
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

!psym=0
cmax=240
cmin=10
plot,x4,0.*f4(*,0,0,0),$
        /xlog,/ylog,$
        yrange=[minm,maxm],$
        xrange=[0.1,2000.],$
        xmargin=[10,3],charthick=2,$
        xstyle=1,ystyle=1,title=intit,$
        linestyle=0,ytitle=ytit,$
        xtitle='!4k!3(!4l!3m)',charsize=cs,thick=2
if iplotstar eq 1 then oplot,lama,nuhnu,linestyle=1
;print,'stellar atm peak',max(nuhnu)

print,intit
if isphav eq 0 then begin
;plot all angles
for it=it1,it2 do begin
for ip=ip1,ip2 do begin
    cosi=1.-(float(it)+0.5)*2./float(ninc)
    theti=acos(cosi)*180./!pi
    ; print,'plotting viewing angle (cosi,theti)',cosi,theti
    if cosi gt 0 then begin
       ccol=(cmax-cmin)/(ninc/2+1)*it+cmin
    endif else begin
       ccol=(cmax-cmin)/(ninc/2+1)*(ninc-it)+cmin
    endelse
  ;don't plot first frequency since it includes all flux in this and shortward.
  oplot,x4(0:nfreq-2),f4(0:nfreq-2,it,ip,nap),linestyle=ls,color=ccol,thick=2
  ;calculate 4.5 micron flux
  f4sm=f4(0:nfreq-2,it,ip,nap)
  x4hi=x4(0:nfreq-2)
  f45=interpol(f4sm,x4hi,4.5)
  print,'inclination, 4.5 um flux',theti,f45
endfor
endfor

endif else begin
;spherically average SED
fphi=total(f4,3,/nan)/float(nphi)
fth=total(fphi,2,/nan)/float(ninc)
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


end

