dir='../models/modcII'
dira='../models/modcII_spiral_nodiff'

;along phi=0.  doesn't average over phi in 3-D.
readcol,dir+'/tmidplane2.dat',i,rstar,rarr,tau,tdust,tdust2,skipline=2
readcol,dira+'/tmidplane2.dat',ia,rstara,rarra,taua,tdusta,tdust2a,skipline=2

xtit='r(AU)'
ytit='T(K)'
xmar=[8,0]
ymar=[3,2]
rin=.1
rout=100.

set_plot,'ps'
device,/inches,xsize=6,ysize=5.,xoffset=1,yoffset=1,$
	filename='tmidplane.ps',/color

!p.multi=[0,1,1]

loadct,10

plot,rarr,tdust,xrange=[rin,rout],/ylog,/xlog,yrange=[10,2000],$
;plot,rarr,tdust2,xrange=[rin,rout],yrange=[1,2000],$
xtitle=xtit,ytitle=ytit,title=tit,charthick=2,$
xstyle=1,ystyle=1,xmargin=xmarg,ymargin=ymarg,charsize=.9,/nodata
;oplot,rarr,tpol,linestyle=3,thick=2

;oplot,rarr,tdust
oplot,rarr,tdust2,color=185
oplot,rarra,tdust2a,color=95

t33=rarr^(-.33)
t40=rarr^(-.40)
ir=50
scale33=tdusta(ir)/t33(ir)
scale40=tdusta(ir)/t40(ir)
;oplot,rarr,scale33*t33,linestyle=0,color=35,thick=2
;oplot,rarr,scale40*t40,linestyle=0,color=135,thick=2
print,max(tdust)

device,/close


end
