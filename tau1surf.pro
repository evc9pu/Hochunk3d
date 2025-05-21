pro tau1surf,dirm,rmax

set_plot,'ps'
device,filename='disksurf.ps',xsize=5.,ysize=4.,$
    xoffset=1.5,yoffset=2,/inches

;rmax=100.
!p.multi=[0,1,1]
;dirm='../models/modcII/'

readcol,dirm+'/taudsurfin.dat',rcyl,z

rad=sqrt(rcyl^2+z^2)
;rmax=max(rad)
idisk=where(rad lt rmax*.9999,ct)

if ct gt 0 then begin
  plot,rcyl(idisk),z(idisk),yrange=[0,rmax],xrange=[0,rmax]
endif else begin
  print,'no Av=1 surface in disk, optically thin?'
endelse

device,/close

end
