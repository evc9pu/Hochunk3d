pro sede,dirm,rstar,tstar,iplotstar,lumscale,intit,d,napin,$
         Jy,atfilein,planckin,Av,flux,cs,op,ls,col,thetst,phist

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
;  convert from nm
;  lama=lama*1.d-3
; new format, already in microns
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

;read in model, skip to aperture   ;**!!!uncomment these on june 6
readcol,dirm+'/peel_flux'+'_'+thetst+'_'+phist+'.dat',nfreqin,junk,napmaxin,numline=1,format='(i,i,i)'
napmax=long(napmaxin(0))
nfreq=long(nfreqin(0))

;napmax=1  ;**!!!delete on june 6
;nfreq=150   ;**!!!delete on june 6

x5=dblarr(nfreq)
f5=x5
x5s=x5
x5d=x5
f5s=x5
f5d=x5

;read in peeled spectra  ;**!!!uncomment these on june 6
readcol,dirm+'/peel_flux'+'_'+thetst+'_'+phist+'.dat',x5,f5,skipline=1+nfreq*(nap),numline=nfreq,format='(d,d)'
readcol,dirm+'/peel_flux_scat'+'_'+thetst+'_'+phist+'.dat',x5s,f5s,skipline=1+nfreq*(nap),$
  numline=nfreq,format='(d,d)'
readcol,dirm+'/peel_flux_dire'+'_'+thetst+'_'+phist+'.dat',x5d,f5d,skipline=1+nfreq*(nap),numline=nfreq,format='(d,d)'
;readcol,dirm+'/peel_flux_dire'+'_'+thetst+'_'+phist+'.dat',x5d,f5d,skipline=1,numline=nfreq,format='(d,d)'
readcol,dirm+'/peel_flux_ther'+'_'+thetst+'_'+phist+'.dat',x5t,f5t,skipline=1+nfreq*(nap),numline=nfreq,format='(d,d)'

;;read in peeled spectra  ;**!!!delete on june 6
;readcol,dirm+'eflux.dat',x5,f5,skipline=1+nfreq*(nap),numline=nfreq,format='(d,d)'
;readcol,dirm+'eflux.dat',x5s,f5s,skipline=1+nfreq*(nap),$
;  numline=nfreq,format='(d,d)'
;readcol,dirm+'eflux.dat',x5d,f5d,skipline=1+nfreq*(nap),numline=nfreq,format='(d,d)'
;;readcol,dirm+'eflux.dat',x5d,f5d,skipline=1,numline=nfreq,format='(d,d)'
;readcol,dirm+'eflux.dat',x5t,f5t,skipline=1+nfreq*(nap),numline=nfreq,format='(d,d)'

;print,'max f5',max(f5,/nan)
f5=f5*fnorm
f5s=f5s*fnorm
f5d=f5d*fnorm
f5t=f5-f5s-f5d
;f5t=f5*fnorm

;integrate fnu*dnu=flam*dlam=(lam*flam)*d(ln(lam))
print,'integrated flux along this direction'
print,'note: not the same as source luminosity since flux varies with viewing angle (so pole-on results are usually higher, edge-on lower)'
;print,x5
;print,f5
sum=0.d0
;print,x5(0),f5(0)
x5hi=alog(x5*1.e-4)
if flux eq 'a' then begin
   f5hi=f5
   sum=int_tabulated(x5hi,f5hi,/sort)
   print,'integrated system flux (lsun)',sum/fnorm
;   print,f5
endif
if flux eq 't' then begin
   f5hi=f5t
   sum=int_tabulated(x5hi,f5hi,/sort)
   print,'integrated thermal flux (lsun)',sum/fnorm
   print,f5t
endif
if flux eq 'd' then begin
   f5hi=f5d
   sum=int_tabulated(x5hi,f5hi,/sort)
   print,'integrated direct stellar flux (lsun)',sum/fnorm
;   print,f5d
endif
if flux eq 's' then begin
   f5hi=f5s
   sum=int_tabulated(x5hi,f5hi,/sort)
   print,'integrated scattered flux (lsun)',sum/fnorm
;   print,f5s
endif

;convert nuhnu and f4 to Jy; units are currently in lamflam=nufnu (erg...)
;first, divide by frequency to get erg/cm2/s/hz
;fnu=lamflam/nu
;c is in microns/s
;x4 is in microns
if Jy eq 1 then begin
  f5=f5/c*x5
  f5s=f5s/c*x5
  f5d=f5d/c*x5
  f5t=f5t/c*x5
  nuhnu=nuhnu/c*lama
  ;now convert to mJy
  f5=f5*1.e26
  f5s=f5s*1.e26
  f5d=f5d*1.e26
  f5t=f5t*1.e26
  nuhnu=nuhnu*1.e26
endif

maxm=max(f5,/nan)*2.
;maxm=1.e-6
minm=maxm*1.e-5
;print,'max',max(f5,/nan)

if jy eq 1 then begin
  ytit='Flux Density (mJy)'
endif else begin
  ytit='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)'
endelse

;add in reddening
xnunew=x5
readcol,'../models/parfiles/kmhnew_extrap.par',klam,junk1,junk2,kkap
;readcol,'../models/parfiles/r400_ice095_extrap.par',klam,junk1,junk2,kkap
;interpolate onto flux grid
kapnew=interpol(kkap,klam,xnunew)
kapv=interpol(kkap,klam,0.55)
taul=kapnew/kapv/1.086*Av
extinct=exp(-taul)

; subtract direct and scattered spectrum (to get thermal only)

f5=f5*extinct
f5s=f5s*extinct
f5d=f5d*extinct
f5t=f5t*extinct

;write out wavelength and flux to a file
forprint,x5,f5,textout=3

if flux eq 't' then f5=f5t
if flux eq 'd' then f5=f5d
if flux eq 's' then f5=f5s

;ii=ninc-1
!psym=0
cmax=240
cmin=10
if op eq 0 then begin
   plot,x5,0.*f5,$
   /xlog,/ylog,$
   yrange=[minm,maxm],$
   xrange=[0.1,5000.],$
   xmargin=[10,3],charthick=2,$
   ;xrange=[alog10(0.2),alog10(5000.)],$
   xstyle=1,ystyle=1,title=intit,$
   linestyle=ls,ytitle=ytit,$
   xtitle='!4k!3(!4l!3m)',charsize=cs,thick=3
endif
;for ii=0,ninc-1,1 do begin
;ccol=(cmax-cmin)/(ninc/2+1)*ii+cmin
;!psym=0
;oplot,(x4(0:nfreq-1)),(f4(nfreq*ii:nfreq*ii+nfreq-1)),linestyle=0,color=ccol,$
; thick=2
;endfor
oplot,lama,nuhnu,linestyle=1
atom=where(lama gt 10.)
lam1=lama(atom(0)-1)
lam2=lama(atom(0))
fnu1=nuhnu(atom(0)-1)
fnu2=nuhnu(atom(0))
nutom=fnu1+(10.-lam1)/(lam2-lam1)*(fnu2-fnu1)
;print,'stellar flux at 10 microns',nutom
;lama(atom(0)-1),lama(atom(0)),nuhnu(atom(0)-1),nuhnu(atom(0))
oplot,x5,f5,thick=3,color=col,linestyle=ls
;oplot,x5,f5-f53,linestyle=2
;oplot,x5,f5+f5e,linestyle=2
;print,'max f5e',max(f5e,/nan)

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

;readcol,'IRAS04325AB.dat',lamdata,fluxdata
;oplot,lamdata,fluxdata,psym=5,symsize=.2

;oplot,[1.1,1.6,2.2,3.6,4.5,5.8,8.0,24.0,70.0],[4.,10.,30.,45.,40.,200.,450.,580.,2500.],psym=5
;oploterr,[1.1,1.6,2.2,3.6,4.5,5.8,8.0,24.0,70.0],[4.,10.,30.,45.,40.,200.,450.,580.,2500.],[1.,2.,5.,7.,6.,25.,35.,50.,1500.]

print,'fnorm',fnorm

end

