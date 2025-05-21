pro sedi_fits,dirm,rstar,tstar,iplotstar,lumscale,intit,d,$
  Jy,atfilein,planckin,cs,op,ls

;print,dirm,rstar,tstar

print,''

loadct,10
ibar=0

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

;read in model
filename=dirm+'/lum_sources.fits.gz'
fluxarr=mrdfits(filename,0,header)
wavearr=mrdfits(filename,1,headwave)

nfreq=fxpar(header,'NAXIS1')
no_init=fxpar(header,'NAXIS2')

print,'nfreq,no_init',nfreq,no_init

x4=fltarr(nfreq)
x4=wavearr.wavelength
f4=fluxarr

f4=f4*fnorm

;f4=f4/lsun

;convert nuhnu and f4 to Jy; units are currently in lamflam=nufnu (erg...)
;first, divide by frequency to get erg/cm2/s/hz
;fnu=lamflam/nu
;c is in microns/s
;x4 is in microns
if Jy eq 1 then begin
  for io=0,no_init-1 do begin
    f4(*,io)=f4(*,io)/c*x4
  endfor
  nuhnu=nuhnu/c*lama
  ;now convert to Jy
  f4=f4*1.e23
  nuhnu=nuhnu*1.e23
endif

;maxm=max(nuhnu,/nan)*100.
maxm=max(f4,/nan)*2.
;maxm=1.e-6
;minm=maxm*2.e-2
minm=maxm*1.e-5
print,'max',max(f4,/nan)

if jy eq 1 then begin
  ytit='Flux Density (Jy)'
endif else begin
  ytit='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)'
endelse

!psym=0
cmax=240
cmin=10
if op eq 0 then begin
   plot,x4,0.*f4(*,0),$
        /xlog,/ylog,$
        yrange=[minm,maxm],$
        xrange=[0.1,2000.],$
 ;       xrange=[1.,10.],$
        xmargin=[10,3],charthick=2,$
        xstyle=1,ystyle=1,title=intit,$
        linestyle=0,ytitle=ytit,$
        xtitle='!4k!3(!4l!3m)',charsize=cs,thick=2
endif
if iplotstar eq 1 then oplot,lama,nuhnu,linestyle=1,color=200
;print,'stellar atm peak',max(nuhnu)

oplot,x4(0:nfreq-2),f4(0:nfreq-2,0),thick=2,linestyle=0
oplot,x4(0:nfreq-2),f4(0:nfreq-2,1),thick=2,linestyle=1
oplot,x4(0:nfreq-2),f4(0:nfreq-2,2),thick=2,linestyle=2
oplot,x4(0:nfreq-2),f4(0:nfreq-2,3),thick=2,linestyle=3
oplot,x4(0:nfreq-2),f4(0:nfreq-2,4),thick=2,linestyle=4

end

