!p.multi=[0,1,1]

set_plot,'ps'
device,/inches,xsize=6.5,ysize=5.5,xoffset=1,yoffset=1,$
	filename='sigma.eps',/encapsulated

;you need to know nrg, sorry (I'm lazy)

nrg=400

dir='../models/modcII_gap_puff'
dir='../models/kate'
readcol,dir+'/sig_theta1.dat',rstar,rau,sig
readcol,dir+'/sig_theta2.dat',rstar2,rau2,sig2

a=where(sig gt 0)
plot,rau(a),sig(a),xtitle='r(AU)',ytitle='Surface Density',$
/ylog,/xlog,yrange=[0.00001,1000.],charsize=1.5,charthick=2,thick=2
oplot,rau2(a),sig2(a),linestyle=1,thick=2

;dir='../models/kate2'
;readcol,dir+'/sig_theta1.dat',rstar,rau,sig
;readcol,dir+'/sig_theta2.dat',rstar2,rau2,sig2

;a=where(sig gt 0)
;oplot,rau(a),sig(a)
;oplot,rau2(a),sig2(a),linestyle=1


device,/close

end