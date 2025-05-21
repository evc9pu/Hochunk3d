pro sed3d720,dirm,rstar,tstar,lumscale,titlin,inctitl,d,napin,$
        Jy,atfilein,planckin,Av,cs,liml,limu,ismooth,pl,thetarrst,phiarrst,phiarr,nfreq,pmin,pmax,override,npeelov

if 0 then begin
dirm='../models/mod1/'
rstar=2.09
tstar=4000.
lumscale=1.d0
titlin=''
inctitl=''
d=430.
napin=1
Jy=0   ;0 for plotting lam*Flam, 1 for plotting Jy
atfilein='../models/parfiles/kt04000g+3.5z-1.5.ascii'
planckin='no'  ;'yes' for plotting planck file, 'no' for input atmosphere file
Av=0
cs=0.9
liml=0.6
limu=-0.1
ismooth=15
pl='sed'
thetarrst=['60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0','60.0']
phiarrst=['0.0','20.0','40.0','60.0','80.0','100.0','120.0','140.0','160.0','180.0','200.0','220.0','240.0','260.0','280.0','300.0','320.0','340.0']
theta=60.0
phiarr=[0.0,20.0,40.0,60.0,80.0,100.0,120.0,140.0,160.0,180.0,200.0,220.0,240.0,260.0,280.0,300.0,320.0,340.0]
nfreq=150
endif


print,''

ibar=0

nap=0l
nap=napin-1l
titl=titlin


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
print,'stellar luminosity calculated from rstar,tstar ',lumstar/lsun
;print,'lsun',lsun
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
;print,'stellar atmosphere peak',max(nuhnu)
;print,size(nuhnu)
;print,'calculated rstar ',sqrt(rstar2)/6.96e10


;read in model, skip to aperture
nfiles=n_elements(thetarrst)
x4=dblarr(nfiles,nfreq)
f4=x4
f4err=x4
q4=x4
q4err=x4
u4=x4
u4err=x4
p4=x4
for i=0,nfiles-1 do begin
  filein='peel_flux_'+thetarrst(i)+'_'+phiarrst(i)+'.dat'
  print,'reading in',dirm+filein
  readcol,dirm+filein,nfreqin,npeelin,napmaxin,numline=1,format='(i,i,i)'
  napmax=long(napmaxin(0))
  nfreq=long(nfreqin(0))
  npeel=long(npeelin(0))
  if (override eq  1) then npeel=npeelov
  print,'nfreq,npeel,napmax'
  ;npeel should equal nfiles
  if nap gt napmax+1 then begin
    print,'error  nap > napmax; stopping program'
    return
  endif
  if napmax gt 1 then begin
     print,'error program only works with 1 aperture right now; stopping'
     return
  endif
  readcol,dirm+'/'+filein,x1,f1,f1err,q1,q1err,u1,u1err,skipline=2
  x4(i,*)=x1
  f4(i,*)=f1*fnorm
  f4err(i,*)=f1err*fnorm
  q4(i,*)=q1
  q4err(i,*)=q1err
  u4(i,*)=u1
  u4err(i,*)=u1err
endfor

a=where(f4 gt 0)
p4=p4*0.d0
p4(a)=sqrt(u4(a)^2+q4(a)^2)*100.
 
;add in reddening
xnunew=x4(0,0:nfreq-1)
readcol,'../models/parfiles/kmh.par',klam,junk1,junk2,kkap
;readcol,'../models/parfiles/r400_ice095.par',klam,junk1,junk2,kkap
;interpolate onto flux grid
kapnew=interpol(kkap,klam,xnunew)
kapv=interpol(kkap,klam,0.55)
taul=kapnew/kapv/1.086*Av
extinct=exp(-taul)

if 1 then begin
for i=0,npeel-1 do begin
  f4(i,*)=f4(i,*)*extinct
  q4(i,*)=q4(i,*)*extinct
  u4(i,*)=u4(i,*)*extinct
endfor
endif

;convert nuhnu and f4 to Jy; units are currently in lamflam=nufnu (erg...)
;first, divide by frequency to get erg/cm2/s/hz
;fnu=lamflam/nu
;c is in microns/s
;x4 is in microns
if Jy eq 1 then begin
  f4=f4/c*x4
  u4=u4/c*x4
  q4=q4/c*x4
  nuhnu=nuhnu/c*lama
  ;now convert to mJy
  f4=f4*1.e26
  q4=q4*1.e26
  u4=u4*1.e26
  nuhnu=nuhnu*1.e26
endif

maxm=max(f4,/nan)*2.
;spot19
;maxm=5.e-10
minm=maxm*2.e-4
;minm=maxm*6.e-2
;spot42
;maxm=1.e-8
;minm=maxm*1.e-4
;minm=maxm*1.e-4
;minm=maxm*2.e-2
;print,'max',max(f4all,/nan)

if jy eq 1 then begin
  ytit='Flux Density (mJy)'
endif else begin
  ytit='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)'
endelse

;title
if inctitl eq 'yes' then begin
  titl=titl+', i='+string(strcompress(thetarrst(0),/remove_all))
endif

nphi=npeel

;lightcurve data
fv=fltarr(nphi)
fr=fltarr(nphi)
fi=fltarr(nphi)
fj=fltarr(nphi)
f36=fltarr(nphi)
f45=fltarr(nphi)
f58=fltarr(nphi)
f80=fltarr(nphi)
f24=fltarr(nphi)
pv=fltarr(nphi)
pi=fltarr(nphi)
pk=fltarr(nphi)
p36=fltarr(nphi)
p45=fltarr(nphi)

ii=nphi-1
print,'nphi ',nphi
!psym=0
cmax=240
cmin=10
if pl eq 'sed' then begin
plot,x4(ii,*),0.*f4(ii,*),$
/xlog,/ylog,$
yrange=[minm,maxm],$
;xrange=[2.,100.],$
xrange=[0.1,1000.],$
xmargin=[10,1],charthick=3,$
xstyle=1,ystyle=1,title=titl,$
linestyle=0,ytitle=ytit,$
xtitle='!4k!3(!4l!3m)',charsize=cs,thick=2
endif
;loadct,0
;oplot,lama,nuhnu,linestyle=1,color=200
loadct,10
;print,'stellar atm peak',max(nuhnu)
nphi2=nphi/2   ;only plot seds for half the phase
for ii=0,nphi-1,1 do begin
; this is phi in this program.  redo some day
;  cosi=(float(ii)+0.5)/float(nphi)
;  theti=acos(cosi)*180./!pi
;  print,'plotting viewing angle (cosi,theti)',cosi,theti
;ccol=(cmax-cmin)/(nphi2+1)*ii+cmin
  ccol=(cmax-cmin)/(nphi+1)*ii+cmin
  !psym=0
  ;don't plot first frequency since it includes all flux in this and
  ;shortward.
  x4hi=x4(ii,0:nfreq-2)
  f4hi=f4(ii,0:nfreq-2)
  f4hi=f4(ii,0:nfreq-2)
  p4hi=p4(ii,0:nfreq-2)
  ismooth2=ismooth
;spot19
;  ismooth2=ismooth-1
;spot42
;  ismooth2=ismooth-3
  f4sm=smooth(f4hi,ismooth2)
  p4sm=smooth(p4hi,ismooth2)
  print,ccol
  if ii lt nphi and pl eq 'sed' then oplot,x4hi,f4sm,linestyle=0,color=ccol,$
   thick=2
;collect data for light curve
  fv(ii)=interpol(f4sm,x4hi,0.55)
  fr(ii)=interpol(f4sm,x4hi,0.67)
  fi(ii)=interpol(f4sm,x4hi,0.8)
  fj(ii)=interpol(f4sm,x4hi,1.2)
  f36(ii)=interpol(f4sm,x4hi,3.6)
  f45(ii)=interpol(f4sm,x4hi,4.5)
  f58(ii)=interpol(f4sm,x4hi,5.8)
  f80(ii)=interpol(f4sm,x4hi,8.0)
  f24(ii)=interpol(f4sm,x4hi,24.0)
  pv(ii)=interpol(p4sm,x4hi,0.55)
  pi(ii)=interpol(p4sm,x4hi,0.8)
  pk(ii)=interpol(p4sm,x4hi,2.2)
  p36(ii)=interpol(p4sm,x4hi,3.6)
  p45(ii)=interpol(p4sm,x4hi,4.5)
endfor
;oplot,x5,f5,thick=5

;plot light curve

if pl eq 'lc' then begin

;plot 2 periods worth
;nphilc=nphi*2
nphilc=nphi
phiarr2=fltarr(nphilc)
fvlc=fltarr(nphilc)
frlc=fltarr(nphilc)
filc=fltarr(nphilc)
fjlc=fltarr(nphilc)
f36lc=fltarr(nphilc)
f45lc=fltarr(nphilc)
f58lc=fltarr(nphilc)
f80lc=fltarr(nphilc)
f24lc=fltarr(nphilc)

fvlc(0:nphi-1)=fv
frlc(0:nphi-1)=fr
filc(0:nphi-1)=fi
fjlc(0:nphi-1)=fj
f36lc(0:nphi-1)=f36
f45lc(0:nphi-1)=f45
f58lc(0:nphi-1)=f58
f80lc(0:nphi-1)=f80
f24lc(0:nphi-1)=f24
phiarr2(0:nphi-1)=phiarr


;fvlc(nphi:nphilc-1)=fv
;frlc(nphi:nphilc-1)=fr
;filc(nphi:nphilc-1)=fi
;fjlc(nphi:nphilc-1)=fj
;f36lc(nphi:nphilc-1)=f36
;f45lc(nphi:nphilc-1)=f45
;f58lc(nphi:nphilc-1)=f58
;f80lc(nphi:nphilc-1)=f80
;f24lc(nphi:nphilc-1)=f24
;phiarr2(nphi:nphilc-1)=phiarr+360.


loadct,10
  delv=-2.5*alog10(fvlc/fvlc(0))
  delr=-2.5*alog10(frlc/frlc(0))
  deli=-2.5*alog10(filc/filc(0))
  delj=-2.5*alog10(fjlc/fjlc(0))
  del36=-2.5*alog10(f36lc/f36lc(0))
  del45=-2.5*alog10(f45lc/f45lc(0))
  del58=-2.5*alog10(f58lc/f58lc(0))
  del80=-2.5*alog10(f80lc/f80lc(0))
  del24=-2.5*alog10(f24lc/f24lc(0))
  
plot,phiarr2,delv,psym=2,yrange=[liml,limu],ystyle=1,xmargin=[8,3],$
        xstyle=1,charsize=cs,charthick=3,title=titl,$
	xrange=[0,720],xtitle='!4u!3',ytitle='!4D!3mag',symsize=.5,/nodata
;oplot,phiarr2,delv,color=100
oplot,phiarr2,delv
oplot,phiarr2,delv,psym=5,symsize=0.7
oplot,phiarr2,deli,color=186
oplot,phiarr2,deli,psym=1,symsize=0.7,color=186
oplot,phiarr2,delj,color=26
oplot,phiarr2,delj,psym=2,symsize=0.7,color=26
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A),/FILL
;;oplot,phiarr2,del58,psym=8,symsize=.4
oplot,phiarr2,del36,color=118
oplot,phiarr2,del36,psym=8,symsize=.4,color=118
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A)
;;oplot,phiarr2,del80,psym=8,symsize=.7
oplot,phiarr2,del45
oplot,phiarr2,del45,psym=8,symsize=.7
;oplot,phiarr2,del80,psym=5,symsize=0.7
;oplot,phiarr2,del24,psym=6,symsize=0.7

endif

if pl eq 'pol' then begin

;plot 2 periods worth
;nphilc=nphi*2
nphilc=nphi
phiarr2=fltarr(nphilc)
pvlc=fltarr(nphilc)
pilc=fltarr(nphilc)
pklc=fltarr(nphilc)
p36lc=fltarr(nphilc)
p45lc=fltarr(nphilc)

pvlc(0:nphi-1)=pv
pilc(0:nphi-1)=pi
pklc(0:nphi-1)=pk
p36lc(0:nphi-1)=p36
p45lc(0:nphi-1)=p45
phiarr2(0:nphi-1)=phiarr

;pvlc(nphi:nphilc-1)=pv
;pilc(nphi:nphilc-1)=pi
;pklc(nphi:nphilc-1)=pk
;p36lc(nphi:nphilc-1)=p36
;p45lc(nphi:nphilc-1)=p45
;phiarr2(nphi:nphilc-1)=phiarr+360.

loadct,10

print,pvlc
  
plot,phiarr2,pvlc,psym=2,yrange=[pmin,pmax],ystyle=1,xmargin=[8,3],$
        xstyle=1,charsize=cs,charthick=3,title=titl,$
	xrange=[0,720],xtitle='!4u!3',ytitle='Polarization (%)',symsize=.5,/nodata
;oplot,phiarr,delv,color=100
oplot,phiarr2,pvlc
oplot,phiarr2,pvlc,psym=5,symsize=0.7
oplot,phiarr2,pilc,color=186
oplot,phiarr2,pilc,psym=1,symsize=0.7,color=186
oplot,phiarr2,pklc,color=90
oplot,phiarr2,pklc,psym=6,symsize=0.7,color=90
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A),/FILL
;oplot,phiarr2,p36lc,psym=8,symsize=.4,color=118
A = FINDGEN(17) * (!PI*2/16.)
USERSYM, COS(A), SIN(A)
;oplot,phiarr2,p45lc,psym=5,symsize=0.7

endif




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

