pro pole,dirm,tit,xr,yr,napin,iapin,nsmooth,op,thetst,phist,ccol,linestyl

lam2=0.3
nap=0l
iap=0l
nap=napin
iap=iapin-1l

;read in model, skip to aperture
readcol,dirm+'/peel_flux'+'_'+thetst+'_'+phist+'.dat',nfreqin,junk,napmaxin,numline=1,format='(i,i,i)'
napmax=long(napmaxin(0))
nfreq=long(nfreqin(0))

x4=dblarr(nfreq)
f4=x4
f4err=x4
q4=x4
q4err=x4
u4=x4
u4err=x4

;read in model, skip to largest aperture
readcol,dirm+'/peel_flux_'+thetst+'_'+phist+'.dat',x4,f4,f4err,q4,q4err,u4,u4err,skipline=3+nfreq*(iap),$
numline=nfreq,format='(d,d)'
print,'nfreq',nfreq

;undelete this sometime....
;a=where(abs(q4err/q4) gt .1, count)
;if count gt 0 then q4(a)=0.
;a=where(abs(u4err/u4) gt .1, count)
;if count gt 0 then u4(a)=0.

;help,x4
;lamarr=fltarr(nfreq)
;lamarr[0:nfreq-1]=x4[0:nfreq-1]
lamarr=x4
a=where(lamarr lt lam2)
nlam2=a(0)-1
if nlam2 lt 0 then nlam2=nfreq
;try not using this since not plotting high-error points
print,'nlam2',nlam2

;loadct,10
loadct,0

if op eq 0 then begin
   plot,x4,q4,$
   ;plot,(x4(0:nfreq-1)),0.*(q4(nfreq*ii:nfreq*ii+nfreq-1)),$
   /xlog,$
   yrange=yr,charthick=2,$
   xrange=xr,$
   ;xrange=[2,5],$
   xstyle=1,ystyle=1,$
   ;ymargin=ymarg,xmargin=xmarg,$
   linestyle=0,$
   ytitle='Q/I',$
   xtitle='wavelength (um)',charsize=1.3,thick=5,title=tit,/nodata
endif

x4plt=x4(0:nlam2-1)
q4plt=q4(0:nlam2-1)
q4plt=smooth(q4plt,nsmooth)
;x4plt=x4
;q4plt=q4
  a=where(q4plt eq 0.,ct)
  if ct gt 0 then q4plt(a)=99.
  oplot,x4plt,q4plt,thick=4,max_value=98.,linestyle=linestyl,col=ccol
oplot,[-100,1000],[0,0],linestyle=1
print,'lambda min plotted',x4plt(nlam2-1)

end
