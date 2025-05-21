pro pol,dirm,tit,xr,yr,napin,op,i1,i2,nsmooth,linestyl

lam2=0.3
nap=0l
nap=napin-1l


;read in model, skip to aperture
readcol,dirm+'/flux.dat',nfreqin,nphiin,nincin,napmaxin,numline=1,format='(i,i,i,i)'
napmax=long(napmaxin(0))
nfreq=long(nfreqin(0))
nphi=long(nphiin(0))
ninc=long(nincin(0))
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

readcol,dirm+'/flux.dat',x4,f4,f4err,q4,q4err,u4,u4err,$
  skipline=nap*nfreq*ninc+2l,numline=nfreq*ninc,format='(d,d)'
print,'nfreq*ninc',nfreq*ninc
print,'nfreq,nphi,ninc,napmax',nfreq,nphi,ninc,napmax
;end
if nap gt napmax+1 then begin
  print,'error  nap > napmax; stopping program'
  return
endif

;a=where(abs(q4err/q4) gt 0.5, count)
;if count gt 0 then q4(a)=0.
;a=where(abs(u4err/u4) gt 0.5, count)
;if count gt 0 then u4(a)=0.

;help,x4
lamarr=fltarr(nfreq)
lamarr[0:nfreq-1]=x4[0:nfreq-1]
a=where(lamarr lt lam2)
nlam2=a(0)-1
if nlam2 lt 0 then nlam2=nfreq
;try not using this if specifying errors correctly
;print,'nlam2',nlam2

loadct,10
;loadct,0
cmax=240
cmin=10

;nfreq=50
;x4=rebin(x4,ninc*nfreq)
;q4=rebin(q4,ninc*nfreq)

ii=9

if op eq 0 then begin
   plot,(x4(0:nfreq-1)),0.*(q4(nfreq*ii:nfreq*ii+nfreq-1)),$
   /xlog,$
   yrange=yr,charthick=2,$
   xrange=xr,$
   ;xrange=[2,5],$
   xstyle=1,ystyle=1,$
   ;ymargin=ymarg,xmargin=xmarg,$
   linestyle=linestyl,$
   ;ytitle='!4k!3F!d!4k!3!n/F',$
   ytitle='Q/I',$
   ;linestyle=0,ytitle='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)',$
   xtitle='wavelength (um)',charsize=1.3,thick=5,title=tit,/nodata
endif

for ii=i1,i2,1 do begin
  cosi=1.-(float(ii)+0.5)*2./float(ninc)
  theti=acos(cosi)*180./!pi
  print,'plotting viewing angle (cosi,theti)',cosi,theti
  if cosi gt 0 then begin
     ccol=(cmax-cmin)/(ninc/2+1)*ii+cmin
  endif else begin
     ccol=(cmax-cmin)/(ninc/2+1)*(ninc-ii)+cmin
  endelse
;plot only down to lam2 min wavelength
  x4plt=x4(0:nlam2-1)
  tmp=q4(nfreq*ii:nfreq*ii+nfreq-1)
  q4plt=tmp(0:nlam2-1)
;plot all wavelengths
;  x4plt=x4(0:nfreq-1)
;  q4plt=q4(nfreq*ii:nfreq*ii+nfreq-1)
  q4plt=smooth(q4plt,nsmooth)
  a=where(q4plt eq 0.)
  q4plt(a)=99.
  oplot,x4plt,q4plt,color=ccol,thick=2,max_value=98.,linestyle=linestyl
;  a=where(x4plt lt 1.29)
;  ij=a(0)
;  cost=float(ii)/10.+0.05
;  thet=acos(cost)*180./!pi
;  print,'theta, Jwave, %Jpol',thet,x4plt(ij),q4plt(ij)*100.
endfor
oplot,[-100,1000],[0,0],linestyle=1
print,'lambda min plotted',x4plt(nlam2-1)

end
