pro plotcm,d,redx,redy,x,y,xlab,ylab,xr,yr,il,nms,$
	x1,y1,x2,y2,av

ninc=10
nmod=6
A = FINDGEN(17) * (!PI*2/16.)
xmin=xr(0)
xmax=xr(1)
ymin=yr(0)
ymax=yr(1)

plot,x(0:ninc-1),y(0:ninc-1),/nodata,$
xmargin=[8,1],ymargin=[4,1],$
xtitle=xlab,$
ytitle=ylab,$
xstyle=1,ystyle=1,$
xrange=xr,$
yrange=yr,$
charthick=2,charsize=1.8

m2=ninc*nmod-1
n1=m2+1
;MS/giant/SG branch
;n2=n1+33-1
;;n2=n1+10-1
;;oplot,x(n1:n2),y(n1:n2),color=0,psym=1,symsize=.6
;n1=n2+1
;n2=n1+42-1
;;oplot,x(n1:n2),y(n1:n2),color=0,psym=4,symsize=.6
;n1=n2+1
;n2=n1+6-1
;;oplot,x(n1:n2),y(n1:n2),color=0,psym=5,symsize=.6
;n1=n2+1
;n2=n1+4-1
;;oplot,x(n1:n2),y(n1:n2),color=0,psym=6,symsize=.6
;n1=n2+1
;n2=n1+1-1
;;oplot,x(n1:n2),y(n1:n2),color=0,psym=7,symsize=.6
;n1=n2+1
;n2=n1+1-1
;oplot,x(n1:n2),y(n1:n2),color=0,psym=7,symsize=1.

cmax=190.
cmin=18.
USERSYM, COS(A), SIN(A),thick=3
;USERSYM, COS(A), SIN(A),/fill
numsym=8
print,''
;for color table 5
;ccol=[15,49,90,116,162,192]
;for color table 13
;ccol=[30,84,149,255,223,207]
ccol=[22,74,149,255,226,203]
for ii=0,nmod-1 do begin
;for ii=0,7 do begin 
;  ccol=fix((cmax-cmin)/float(nmod-1)*ii+cmin)
	print,ccol(ii)
  m1=ii*ninc
  m2=m1+ninc-1
;  print,m1,m2
  for iinc=0,ninc-1,1 do begin
;    ccol=(cmax-cmin)/ninc*iinc+cmin
    ssze=0.5+0.6*iinc/(ninc-1)
    plots,x(iinc+m1),y(iinc+m1),color=ccol(ii),psym=numsym,symsize=ssze,thick=2
  endfor
endfor
USERSYM, COS(A), SIN(A),/fill
;oplot,x(n1:n2),y(n1:n2),color=0
;for i=n1,n2 do begin 
;  oplot,redx+x(i),redy+y(i),linestyle=1
;endfor

stupid=0

;oplot,redx+x1,redy+y1,linestyle=1
;oplot,redx+x2,redy+y2,linestyle=1
;ihi=where(Av eq 30)
;if redy(ihi(0))+y2 lt ymax and redx(ihi(0))+x2 lt xmax then plots,redx(ihi)+x2,redy(ihi)+y2,psym=8,symsize=.4
;ihi=where(Av eq 60)
;if redy(ihi(0))+y2 lt ymax and redx(ihi(0))+x2 lt xmax then plots,redx(ihi)+x2,redy(ihi)+y2,psym=8,symsize=.4
;ihi=where(Av eq 90)
;if redy(ihi(0))+y2 lt ymax and redx(ihi(0))+x2 then plots,redx(ihi)+x2,redy(ihi)+y2,psym=8,symsize=.4
;ihi=where(Av eq 120)
;if redy(ihi(0))+y2 lt ymax and redx(ihi(0))+x2 then plots,redx(ihi)+x2,redy(ihi)+y2,psym=8,symsize=.4
;ihi=where(Av eq 150)
;if redy(ihi(0))+y2 lt ymax and redx(ihi(0))+x2 then plots,redx(ihi)+x2,redy(ihi)+y2,psym=8,$
;symsize=.4
;ihi=where(Av eq 180)
;if redy(ihi(0))+y2 lt ymax and redx(ihi(0))+x2 then plots,redx(ihi)+x2,redy(ihi)+y2,psym=8,$
;symsize=.4
;ihi=where(Av eq 210)
;if redy(ihi(0))+y2 lt ymax and redx(ihi(0))+x2 then plots,redx(ihi)+x2,redy(ihi)+y2,psym=8,$
;symsize=.4
;ihi=where(Av eq 240)
;if redy(ihi(0))+y2 lt ymax and redx(ihi(0))+x2 then plots,redx(ihi)+x2,redy(ihi)+y2,psym=8,$


symsize=.4
;some labels
cs=.7
if il eq 1 then begin
;  xyouts,-1,0,'MS',charsize=cs
;  xyouts,1,.6,'1 Lsun',charsize=cs
;  xyouts,2.5,1.2,'Pre-UCHII',charsize=cs
;  xyouts,4,2.3,'UCHII',charsize=cs
endif

USERSYM, COS(A), SIN(A),thick=3
;xc=findgen(150)/150.*0.45*(xmax-xmin)+xmin+0.19*(xmax-xmin)
;yc=.92*(ymax-ymin)+ymin
if il eq 1 then begin
  xc=0.08*(xmax-xmin)+xmin
  yc=ymax-0.11*(ymax-ymin)-findgen(nmod)/float(nmod)*.4*(ymax-ymin)
;  print,ymin,ymax,yc
;  yc=findgen(nmod)/float(nmod)*.4*(ymax-ymin)+ymin+0.6*(ymax-ymin)
  clab=['0','Late 0','I','Late I','II','III']
  xyouts,xc+.04*(xmax-xmin),yc(0)+.04*(ymax-ymin),'Class',charsize=.8,charthick=2
  for ii=0,nmod-1 do begin
;    ccol=(cmax-cmin)/(nmod-1)*ii+cmin
;    ccol=fix((cmax-cmin)/float(nmod-1)*ii+cmin)
;    print,ccol
    xyouts,xc+.07*(xmax-xmin),yc(ii)-.02*(ymax-ymin),clab(ii),charsize=.8,charthick=2
    plots,xc,yc(ii),color=ccol(ii),psym=8,symsize=.6
  endfor
;for i=0,149 do begin
;  ccol=(cmax-cmin)/150.*i+cmin
;  print,ccol
;  plots,xc(i),yc,color=ccol,psym=8,symsize=.6
;endfor
;yc2=.85*(ymax-ymin)+ymin
;xyouts,xmin+.08,yc2,'Class:',charsize=.8,charthick=2
;i=fix((cmin+1-cmin)/(cmax-cmin)*150.)
;xyouts,xc(i),yc2,'0',charsize=.8,charthick=2
;i=fix((56.-10.-cmin)/(cmax-cmin)*150.)
;xyouts,xc(i),yc2,'0-I',charsize=.8,charthick=2
;i=fix((92.-cmin)/(cmax-cmin)*150.)
;xyouts,xc(i),yc2,'I',charsize=.8,charthick=2
;i=fix((128.-10.-cmin)/(cmax-cmin)*150.)
;xyouts,xc(i),yc2,'I-II',charsize=.8,charthick=2
;i=fix((164.-cmin)/(cmax-cmin)*150.)
;xyouts,xc(i),yc2,'II',charsize=.8,charthick=2
;i=fix((199.-cmin)/(cmax-cmin)*150.)
;xyouts,xc(i),yc2,'III',charsize=.8,charthick=2
endif

;errorbar
;say, 10% photometric error.
;error of ratio of fluxes is sqrt(2*.1^2)=14%
;in magnitude space, error ranges between -2.5log(F2/F1*(1+.14)
; and -2.5log(F2/F1*(1-.14)) 
; -2.5*alog10(1-.14)=.16, 2.5*alog10(1.14)=.14, so use 0.15 as error
; in mag. space
;x1=xc(99)+.75
;y1=yc-.25
;err=0.15
;plots,[x1,x1+2*err],[y1+err,y1+err]
;plots,[x1+err,x1+err],[y1,y1+2*err]
end
