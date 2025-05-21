close,1
openw,1,'pmsfluxes.txt'
printf,1,'PMS fluxes at J,H,K,3.6,4.5,5.8,8,24,70,160,N'
printf,1,'each row plots the 11 wavelengths in order as above'
printf,1,'10 inclinations per model'
printf,1,'units are Jy'

;large aperture
nap=50
;middle aperture
;nap=2
;small aperture
;nap=1

ninc=10

sed_pap,'/d/glimpse142/tom/verif/3000001','3000001',f1,nap
f=f1
for i=0,ninc-1 do begin
 printf,1,format='(12e12.4)',f(i,0),f(i,1),f(i,2),f(i,3),f(i,4),f(i,5),f(i,6),$
	f(i,7),f(i,8),f(i,9),f(i,10),f(i,11)
endfor

end
