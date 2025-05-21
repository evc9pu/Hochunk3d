pro pol_fits,dirm,tit,xr,yr,napin,op,i1,i2,nsmooth,linestyl

lam2=0.3
nap=0l
nap=napin-1l

;read in model
filename=dirm+'/flux_hypercube.fits.gz'
fluxarr=mrdfits(filename,0,header)
wavearr=mrdfits(filename,1,headwave)

nfreq=fxpar(header,'NAXIS1')
ninc=fxpar(header,'NAXIS2')
nphi=fxpar(header,'NAXIS3')
napmax=fxpar(header,'NAXIS4')
no=fxpar(header,'NAXIS5')

x=fltarr(nfreq)
x=wavearr.wavelength
f=fltarr(nfreq,ninc,nphi,napmax)
f=fluxarr(*,*,*,*,0,0)
f2=fluxarr(*,*,*,*,0,1)
q=fluxarr(*,*,*,*,0,2)
q2=fluxarr(*,*,*,*,0,3)
u=fluxarr(*,*,*,*,0,4)
u2=fluxarr(*,*,*,*,0,5)
v=fluxarr(*,*,*,*,0,6)
v2=fluxarr(*,*,*,*,0,7)

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
cmax=240
cmin=10

if op eq 0 then begin
   plot,x,0.*p(*,0,0,0),$
   /xlog,$
   yrange=yr,charthick=2,$
   xrange=xr,$
   xstyle=1,ystyle=1,$
   ;ymargin=ymarg,xmargin=xmarg,$
   linestyle=linestyl,$
   ;ytitle='!4k!3F!d!4k!3!n/F',$
   ytitle='% Polarization',$
   ;linestyle=0,ytitle='!4k!3F!d!4k!3!n (ergs s!u-1!n cm!u-2!n)',$
   xtitle='wavelength (um)',charsize=1.3,thick=5,title=tit,/nodata
endif

for it=0,ninc-1 do begin
for ip=0,nphi-1 do begin
  cosi=1.-(float(it)+0.5)*2./float(ninc)
  theti=acos(cosi)*180./!pi
  print,'plotting viewing angle (cosi,theti)',cosi,theti
  if cosi gt 0 then begin
     ccol=(cmax-cmin)/(ninc/2+1)*it+cmin
  endif else begin
     ccol=(cmax-cmin)/(ninc/2+1)*(ninc-it)+cmin
  endelse
;plot only down to lam2 min wavelength
  xplt=x(0:nlam2-1)
  pplt=100.*p(0:nlam2-1,it,ip,nap)
  pplt=smooth(pplt,nsmooth)
  a=where(pplt eq 0., ct)
  if ct gt 0 then pplt(a)=199.
  oplot,xplt,pplt,color=ccol,thick=2,max_value=198.,linestyle=linestyl
endfor
endfor
oplot,[-100,1000],[0,0],linestyle=1
print,'lambda min plotted',xplt(nlam2-1)

end