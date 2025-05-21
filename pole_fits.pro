pro pole_fits,dirm,tit,xr,yr,napin,iapin,nsmooth,op,thetst,phist,ccol,linestyl

lam2=0.3
nap=0l
iap=0l
nap=napin
iap=iapin-1l

;read in model
filename=dirm+'/peel_hypercube.fits.gz'
fluxarr=mrdfits(filename,0,header)
wavearr=mrdfits(filename,1,headwave)

nfreq=fxpar(header,'NAXIS1')
npeel=fxpar(header,'NAXIS2')
napmax=fxpar(header,'NAXIS3')
no=fxpar(header,'NAXIS4')

x=fltarr(nfreq)
x=wavearr.wavelength
f=fltarr(nfreq,npeel,napmax)
f=fluxarr(*,*,*,0,0)
f2=fluxarr(*,*,*,0,1)
q=fluxarr(*,*,*,0,2)
q2=fluxarr(*,*,*,0,3)
u=fluxarr(*,*,*,0,4)
u2=fluxarr(*,*,*,0,5)
v=fluxarr(*,*,*,0,6)
v2=fluxarr(*,*,*,0,7)

p=sqrt(q^2+u^2+u^2)
pa=0.5*atan(u,q)

; calculate errors later.  TBD!

;a=where(abs(q4err/q4) gt 0.5, count)
;if count gt 0 then q4(a)=0.
;a=where(abs(u4err/u4) gt 0.5, count)
;if count gt 0 then u4(a)=0.

a=where(x lt lam2)
nlam2=a(0)-1
if nlam2 lt 0 then nlam2=nfreq
;try not using this if specifying errors correctly
;print,'nlam2',nlam2

loadct,10
;loadct,0

if op eq 0 then begin
   plot,x,p,$
   /xlog,$
   yrange=yr,charthick=2,$
   xrange=xr,$
   xstyle=1,ystyle=1,$
   linestyle=0,$
   ytitle='% Polarization',$
   xtitle='wavelength (um)',charsize=1.3,thick=5,title=tit,/nodata
endif

for it=0,npeel-1 do begin

  xplt=x(0:nlam2-1)
  pplt=100.*p(0:nlam2-1,it,iap)
  pplt=smooth(pplt,nsmooth)
  a=where(pplt eq 0.,ct)
  if ct gt 0 then pplt(a)=199.
  oplot,xplt,pplt,color=ccol,thick=4,max_value=198.,linestyle=linestyl

endfor

oplot,[-100,1000],[0,0],linestyle=1
print,'lambda min plotted',xplt(nlam2-1)

end
