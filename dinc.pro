pro dinc,dir,tit,xtit,ytit,xmarg,ymarg,rin,rout

file=dir+'/darr.unf'

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

print,'max density',max(Tdust)
peak=max(Tdust)
rhomin=1.e-10*peak
rhomax=3.*peak

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
plot,rarr,tpol,xrange=[rin,rout],/ylog,/xlog,yrange=[rhomin,rhomax],$
xtitle=xtit,ytitle=ytit,title=tit,charthick=2,$
xstyle=1,ystyle=1,xmargin=xmarg,ymargin=ymarg,charsize=.9,/nodata
oplot,rarr,0.01*peak*rarr^(-1.5),linestyle=0,color=35,thick=2
;oplot,rarr,peak*rarr^(-2.0),linestyle=0,color=35,thick=2
;oplot,rarr,tpol,linestyle=3,thick=2
print,'thetapol_max',thetarr(ihi)*180./!pi

for iinc=0,9 do begin
   theta=float(iinc)/(9.)*90.
   for j=0,nrg-1 do begin
      icount=0.
      tpol(j)=0.
      a=where(thetarr gt theta*!pi/180.)
      ihi=a(0)
      tpol(j)=tdust(j,ihi,0)
   endfor
   oplot,rarr,tpol,thick=1
   print,'theta',thetarr(ihi)*180./!pi,max(tpol)
endfor

;spherically averaged
for j=0,nrg-1 do begin
   tpol(j)=total(tdust(j,0:ntg-1,0))/float(ntg)
endfor
oplot,rarr,tpol,thick=5,color=155
print,max(tpol)

end
