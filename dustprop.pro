!p.multi=[0,2,2]

set_plot,'ps'
device,/inches,xsize=6.5,ysize=5.5,xoffset=1,yoffset=1,$
	filename='dustprop.ps'

;readcol,'../models/parfiles/www003.par',lam,cext,csca,kap,hgg,pmax
;readcol,'../models/parfiles/ww04.par',lam2,cext2,csca2,kap2,hgg2,pmax2
;readcol,'../models/parfiles/r400_ice095.par',lam3,cext3,csca3,kap3,hgg3,pmax3
;readcol,'../models/parfiles/kmh.par',lam4,cext4,csca4,kap4,hgg4,pmax4

readcol,'../models/parfiles/kmh.par',lam,cext,csca,kap,hgg,pmax
readcol,'../models/parfiles/kmh_no200_rescale.par',lam2,cext2,csca2,kap2,hgg2,pmax2
readcol,'../models/parfiles/draine_opac_new.dat',lam3,cext3,csca3,kap3,hgg3,pmax3
;readcol,'../models/parfiles/r400_ice095_no200_rescale.par',lam4,cext4,csca4,kap4,hgg4,pmax4

x1=.01
x2=2000.

plot,lam2,kap2,/ylog,/xlog,xrange=[x1,x2],yrange=[.003,10000],$
xstyle=1,ystyle=1,$
ymargin=[3,2],xmargin=[10,1],$
linestyle=1,ytitle='!4j!3',$
;xtitle='!4k!3(!4l!3m)',$
charsize=1.,thick=5,/nodata,charthick=2
;oplot,lam,kap,linestyle=2
oplot,lam2,kap2,linestyle=1
oplot,lam3,kap3,linestyle=0
;oplot,lam4,kap4,linestyle=3

plot,lam2,csca2/cext2,/xlog,xrange=[x1,x2],yrange=[0,.6],$
xstyle=1,ystyle=1,charthick=2,$
ymargin=[3,2],xmargin=[10,1],$
linestyle=1,ytitle='!4x!3',$
;xtitle='!4k!3(!4l!3m)',$
charsize=1.,thick=5,/nodata
;oplot,lam,csca/cext,linestyle=2
oplot,lam2,csca2/cext2,linestyle=1
oplot,lam3,csca3/cext3,linestyle=0
;oplot,lam4,csca4/cext4,linestyle=3

plot,lam2,hgg2,/xlog,xrange=[x1,x2],$
xstyle=1,ystyle=1,charthick=2,$
ymargin=[4,1],xmargin=[10,1],$
linestyle=1,ytitle='g',$
xtitle='!4k!3(!4l!3m)',charsize=1.,thick=5,/nodata
;oplot,lam,hgg,linestyle=2
oplot,lam2,hgg2,linestyle=1
oplot,lam3,hgg3,linestyle=0
;oplot,lam4,hgg4,linestyle=3

plot,lam2,pmax2,/xlog,xrange=[x1,x2],yrange=[0,1],$
xstyle=1,ystyle=1,charthick=2,$
ymargin=[4,1],xmargin=[10,1],$
linestyle=1,ytitle='P!dl!n',$
xtitle='!4k!3(!4l!3m)',charsize=1.,thick=5,/nodata
;oplot,lam,pmax,linestyle=2
oplot,lam2,pmax2,linestyle=1
oplot,lam3,pmax3,linestyle=0
;oplot,lam4,pmax4,linestyle=3

device,/close
end
