pro tinc_a,dir,tit,xtit,ytit,xmarg,ymarg

file=dir+'/tarr.unf'

loadct,10
cmax=230
cmin=10

nrg=1L
ntg=1L
npg=1L

close,1
openr,1,file,/f77_unformatted
readu,1,nrg,ntg,npg
print,'nrg,ntg,npg',nrg,ntg,npg
rarr=dblarr(nrg)
thetarr=dblarr(ntg)
phiarr=dblarr(npg)
tdust=dblarr(nrg,ntg,npg)
readu,1,rarr
readu,1,thetarr
readu,1,phiarr
readu,1,tdust
close,1

print,'max T',max(Tdust)

iave=1

;window,ip+1,xsize=500,ysize=400,retain=2
tpol=fltarr(nrg)
for j=0,nrg-1 do begin
icount=0.
tpol(j)=0.
ihi=iave
for i=0,ihi do begin
  if tdust(j,i,0) gt 0 then begin
	tpol(j)=tpol(j)+tdust(j,i,0)
	icount=icount+1.
  endif
endfor
tpol(j)=tpol(j)/icount
endfor
plot,rarr,tpol,xrange=[3e-2,9000],/ylog,/xlog,yrange=[1,5000],$
xtitle=xtit,ytitle=ytit,title=tit,charthick=2,$
xstyle=1,ystyle=1,xmargin=xmarg,ymargin=ymarg,charsize=.9,/nodata
oplot,rarr,180.*rarr^(-.33),linestyle=0,color=35,thick=2
oplot,rarr,tpol,linestyle=3,thick=2
print,'thetapol_max',thetarr(ihi)*180./!pi


;30 degrees
for j=0,nrg-1 do begin
icount=0.
tpol(j)=0.
a=where(thetarr gt 30.*!pi/180.)
ihi=a(0)
for i=ihi,ihi+iave do begin
  if tdust(j,i,0) gt 0 then begin
	tpol(j)=tpol(j)+tdust(j,i,0)
	icount=icount+1.
  endif
endfor
tpol(j)=tpol(j)/icount
endfor
;oplot,rarr,tpol,linestyle=0,color=cmin+.67*(cmax-cmin)
oplot,rarr,tpol,linestyle=2,thick=2
print,'theta',thetarr(ihi+5)*180./!pi
a=where(tpol eq 1600.00,count)
if count gt 0 then print,'T=1600 at pole at r = ',rarr(a)

;60 degrees
for j=0,nrg-1 do begin
icount=0.
tpol(j)=0.
a=where(thetarr gt 60.*!pi/180.)
ihi=a(0)
for i=ihi,ihi+iave do begin
  if tdust(j,i,0) gt 0 then begin
	tpol(j)=tpol(j)+tdust(j,i,0)
	icount=icount+1.
  endif
endfor
tpol(j)=tpol(j)/icount
endfor
;oplot,rarr,tpol,linestyle=0,color=.33*(cmax-cmin)
oplot,rarr,tpol,linestyle=1,thick=2
print,'theta',thetarr(ihi+5)*180./!pi

;90 degrees
for j=0,nrg-1 do begin
icount=0.
tpol(j)=0.
ihi=90
a=where(thetarr lt 90.*!pi/180.)
nhi=n_elements(a)
ihi=a(nhi-1)
for i=ihi,ihi+iave do begin
  if tdust(j,i,0) gt 0 then begin
	tpol(j)=tpol(j)+tdust(j,i,0)
	icount=icount+1.
  endif
endfor
tpol(j)=tpol(j)/icount
endfor
;oplot,rarr,tpol,linestyle=0,color=cmin
oplot,rarr,tpol,thick=2
print,'theta',thetarr(ihi+5)*180./!pi
;print,rarr(0:40),tpol(0:40)
a=where(tpol eq 1600.00,count)
if count gt 0 then print,'t=1600 at equator at r = ',rarr(a)
print,rarr(0:40),tpol(0:40)

end
