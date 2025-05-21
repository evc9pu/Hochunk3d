pro plotcol,d,redx,redy,redx4,redy4,x,y,xlab,ylab,xmin,xmax,ymin,ymax,il,ninc,$
	av1,av2,avmax,av

nmod=1
A = FINDGEN(17) * (!PI*2/16.)

plot,x(0:ninc-1),y(0:ninc-1),/nodata,$
xmargin=[8,1],ymargin=[4,1],$
xtitle=xlab,$
ytitle=ylab,$
xstyle=1,ystyle=1,$
xrange=[xmin,xmax],yrange=[ymin,ymax],charthick=2,charsize=1.8

m2=ninc*nmod-1
n1=m2+1
;n2=m2+1+nms-1
;MS/giant/SG branch
;n2=n1+33-1
;oplot,x(n1:n2),y(n1:n2),color=0,psym=1,symsize=.6
;n1=n2+1
;n2=n1+42-1
;oplot,x(n1:n2),y(n1:n2),color=0,psym=4,symsize=.6
;n1=n2+1
;n2=n1+6-1
;;oplot,x(n1:n2),y(n1:n2),color=0,psym=5,symsize=.6
;n1=n2+1
;n2=n1+2-1
;oplot,x(n1:n2),y(n1:n2),color=0,psym=6,symsize=.6
;n1=n2+1
;n2=n1+2-1
;oplot,x(n1:n2),y(n1:n2),color=0,psym=2,symsize=.6
;n1=n2+1
;n2=n1+1-1
;oplot,x(n1:n2),y(n1:n2),color=0,psym=7,symsize=.6
;n1=n2+1
;n2=n1+1-1
;oplot,x(n1:n2),y(n1:n2),color=0,psym=7,symsize=1.

cmax=190.
cmin=18.
USERSYM, COS(A), SIN(A),thick=3
;USERSYM, COS(A), SIN(A),/fill
numsym=8
print,''
;
;for color table 5
;ccol=[5,55,90,116,162,192]
;for color table 13
ccol=[22,74,149,255,226,203]
for ii=0,nmod-1 do begin
  print,ccol(ii)
  m1=ii*ninc
  m2=m1+ninc-1
  for iinc=0,ninc-1,1 do begin
    ssze=0.5+0.6*iinc/(ninc-1)
    plots,x(iinc+m1),y(iinc+m1),color=ccol(ii),psym=numsym,symsize=ssze,thick=2
  endfor
endfor
USERSYM, COS(A), SIN(A),/fill

ihi=where(Av gt Avmax)
ihi2=ihi(0)
av3=redx(ihi2)+av1
av4=redy(ihi2)+av2
arrow,av1,av2,av3,av4,/data
;print,av1,av2,av3,av4
;xyouts,av3,av4,strcompress(string(avmax)),charsize=.8

av3=redx4(ihi2)+av1
av4=redy4(ihi2)+av2
arrow,av1,av2,av3,av4,/data,thick=2
print,av1,av2,av3,av4
xyouts,av3,av4,strcompress(string(avmax)),charsize=.8


USERSYM, COS(A), SIN(A),thick=3
if il eq 1 then begin
  xc=0.08*(xmax-xmin)+xmin
  yc=ymax-0.11*(ymax-ymin)-findgen(nmod)/float(nmod)*.4*(ymax-ymin)
  clab=['0','Late 0','I','Late I','II','III']
  xyouts,xc+.04*(xmax-xmin),yc(0)+.04*(ymax-ymin),'Class',charsize=.8,charthick=2
  for ii=0,nmod-1 do begin
    xyouts,xc+.07*(xmax-xmin),yc(ii)-.02*(ymax-ymin),clab(ii),charsize=.8,charthick=2
    plots,xc,yc(ii),color=ccol(ii),psym=8,symsize=.6
  endfor
endif

end
