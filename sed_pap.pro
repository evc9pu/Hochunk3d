pro sed_pap,dirm,titl,fluxfilt,napin

;loadct,10

;nfreq=250
;ninc=10

nap=0l
nap=napin-1l

;calculate fnorm
;fnorm=L/(4*pi*d^2)=(4*pi*R^2*sigma*T^4)/(4*pi*d^2)=(R^2*sigma*T^4)/d^2
;using model or specific input
c=2.9979e14
sigma=5.67e-5
pc=3.0857e18
lsun=3.845e33
rstar=2.09
tstar=4000.
d=10.   ;for conversion to absolute Mag in later program
lumstar=4.*!pi*(rstar*6.955e10)^2*sigma*tstar^2*tstar^2
print,'stellar luminosity ',lumstar/lsun
print,'...not used...'
fnorm=1.d0/(4.*!pi)/d^2/pc/pc*lsun
print,'fnorm',fnorm

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
f4err=x4
q4=x4
q4err=f4
u4=f4
u4err=f4
skipline=0l
if nphi gt 2 then begin
   print,'error, this is only for 2-D models, nphi=2; stopping program'
   return
endif

readcol,dirm+'/flux.dat',x4,f4,f4err,q4,q4err,u4,u4err,skipline=nap*nfreq*ninc+2l,numline=nfreq*ninc,format='(d,d)'

;readcol,dirm+'/flux.dat',x4,f4,p4,skipline=nap*nfreq*ninc+1l,numline=nfreq*ninc,format='(d,d)'

;read in peeled image
;readcol,dirm+'/eflux.dat',px4,pf4,skipline=3+nfreq*2

f4=f4*fnorm
p4=sqrt(q4^2+u4^2)

;convert nuhnu and f4 to Jy; units are currently in lamflam=nufnu (erg...)
;first, divide by frequency to get erg/cm2/s/hz
;fnu=lamflam/nu
;c is in microns/s
;x4 is in microns
f4=f4/c*x4
;pf4=pf4/c*px4
f4=f4*1.e23

maxm=5.e-8
minm=1.e-14
maxm=max(f4)
;maxm=1000.
minm=maxm*1.e-8
print,'max',max(f4)

ii=9
!psym=0
cmax=240
cmin=10
plot,(x4(0:nfreq-1)),0.*(f4(nfreq*ii:nfreq*ii+nfreq-1)),$
/xlog,/ylog,$
yrange=[minm,1.05*maxm],$
xrange=[0.2,5000.],$
;xrange=[alog10(0.2),alog10(5000.)],$
xstyle=1,ystyle=1,title=titl,$
;linestyle=0,ytitle='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)',$
linestyle=0,ytitle='Flux Density (Jy)',$
;xtitle='!4k!3(!4l!3m)',charsize=1.5,thick=5
xtitle='!4k!3(!4l!3m)',charsize=1.5,thick=5
for ii=0,9,1 do begin
ccol=(cmax-cmin)/11*ii+cmin
;oplot,(x(0:nfreq-1)),(f1(nfreq*ii:nfreq*ii+nfreq-1)),linestyle=2
;oplot,(x(0:nfreq-1)),(f2(nfreq*ii:nfreq*ii+nfreq-1)),linestyle=3
!psym=0
oplot,(x4(0:nfreq-1)),(f4(nfreq*ii:nfreq*ii+nfreq-1)),linestyle=0
!psym=0
;oplot,(x(0:nfreq-1)),(f(nfreq*ii:nfreq*ii+nfreq-1)),linestyle=1
endfor

;now integrate over filters to get 2mass & IRAC magnitudes

;used just for testing
filtlam=[1.2,1.65,2.16,3.6,4.5,5.8,8.0]
;print,interpol([655.,276.,248.,160,50],[2.179,3.547,3.761,4.769,8.756],5.8)
filtf=[1630.,1050.,667.,289.,183.,131.,71.]
mag=fltarr(ninc,10)

fluxfilt=fltarr(ninc,11)
;loop over inclination
for ii=0,ninc-1 do begin
;interpolate x4,f4 to filters
  n1=ii*nfreq
  n2=n1+nfreq-1
  f=f4(n1:n2)
  pol=p4(n1:n2)
  lam=x4(n1:n2)
; integrate over filter functions
  convfilt_jy,'../models/parfiles/2J.txt',1,lam,f,flux,0
  fluxfilt(ii,0)=flux
  convfilt_jy,'../models/parfiles/2H.txt',1,lam,f,flux,0
  fluxfilt(ii,1)=flux
  convfilt_jy,'../models/parfiles/2K.txt',1,lam,f,flux,0
  fluxfilt(ii,2)=flux
  convfilt_jy,'../models/parfiles/I1.txt',1,lam,f,flux,0
  fluxfilt(ii,3)=flux
  convfilt_jy,'../models/parfiles/I2.txt',1,lam,f,flux,0
  fluxfilt(ii,4)=flux
  convfilt_jy,'../models/parfiles/I3.txt',1,lam,f,flux,0
  fluxfilt(ii,5)=flux
  convfilt_jy,'../models/parfiles/I4.txt',1,lam,f,flux,0
  fluxfilt(ii,6)=flux
  convfilt_jy,'../models/parfiles/M1.txt',1,lam,f,flux,0
  fluxfilt(ii,7)=flux
  convfilt_jy,'../models/parfiles/M2.txt',1,lam,f,flux,0
  fluxfilt(ii,8)=flux
  convfilt_jy,'../models/parfiles/M3.txt',1,lam,f,flux,0
  fluxfilt(ii,9)=flux
; k-band pol
  convfilt_jy,'../models/parfiles/2K.txt',1,lam,pol,flux,0
  fluxfilt(ii,10)=flux

;test, interpolate x4,f4 to filters
   fnew=interpol(f,lam,filtlam)
   pnew=interpol(pol,lam,filtlam)
;   mag(ii,0:6)=-2.5*alog10(fnew/filtf)
   mag(ii,0:6)=fnew
print,'fluxes integrated over filters'
print,fluxfilt(ii,*)
print,'fluxes interpolated to filter wavelengths'
print,fnew
print,'polarization'
print,pnew
endfor

end

