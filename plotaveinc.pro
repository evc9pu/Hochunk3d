!p.multi=[0,2,3]
set_plot,'ps' 
if !d.name eq 'PS' then device,xsize=6.5,ysize=9.,/inches,$
/portrait,xoffset=1.0,yoffset=1.0,filename='aveinc.ps'

dirm='../models/modcII'

;choose your aperture (from 1 to NAP in mctherm.par)
iap=1l  ; l stands for long integer, I have 5 long here.

;boxcar smoothing if you want, set to 0 if you don't
ns=5

for iinc=1,10 do begin

   aveinc,dirm,iap,iinc,ns

endfor

device,/close

end

