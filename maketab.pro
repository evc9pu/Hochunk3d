close,1
openw,1,'pmsfluxes.txt'
printf,1,'PMS fluxes at J,H,K,3.6,4.5,5.8,8,24,70,160'
printf,1,'each row plots the 10 wavelengths in order as above and K-band polarization in the last column'
printf,1,'10 inclinations per model'
printf,1,'units are Jy'

;pick an aperture, from 1-nap (currently code uses 5 apertures but you can
; change this in srcdust/tts.txt.  see aperture sizes in disk.dat file)
; iap=2   ;406 AU aperture (this will be redder than Figs 7a,8a in paper II)
 iap=3   ;1325 AU aperture  (slightly bluer than Figs 7a,8a)
; iap=5   ;5000 AU aperture  (should agree with figs 7b,8b; however, slight
; differences due to accretion luminosity changes in code).

ninc=10

sed_pap,'../models/modc0','Stage 0',f1,iap
f=f1
for i=0,ninc-1 do begin
 printf,1,format='(12e12.4)',f(i,0),f(i,1),f(i,2),f(i,3),f(i,4),f(i,5),f(i,6),$
	f(i,7),f(i,8),f(i,9),f(i,10)
endfor

sed_pap,'../models/modc0-I','Stage 0-I',f1,iap
f=f1
for i=0,ninc-1 do begin
 printf,1,format='(12e12.4)',f(i,0),f(i,1),f(i,2),f(i,3),f(i,4),f(i,5),f(i,6),$
	f(i,7),f(i,8),f(i,9),f(i,10)
endfor

sed_pap,'../models/modcI','Stage I',f1,iap
f=f1
for i=0,ninc-1 do begin
 printf,1,format='(12e12.4)',f(i,0),f(i,1),f(i,2),f(i,3),f(i,4),f(i,5),f(i,6),$
	f(i,7),f(i,8),f(i,9),f(i,10)
endfor

sed_pap,'../models/modcI-II','Stage I-II',f1,iap
f=f1
for i=0,ninc-1 do begin
 printf,1,format='(12e12.4)',f(i,0),f(i,1),f(i,2),f(i,3),f(i,4),f(i,5),f(i,6),$
	f(i,7),f(i,8),f(i,9),f(i,10)
endfor

sed_pap,'../models/modcII','Stage II',f1,iap
f=f1
for i=0,ninc-1 do begin
 printf,1,format='(12e12.4)',f(i,0),f(i,1),f(i,2),f(i,3),f(i,4),f(i,5),f(i,6),$
	f(i,7),f(i,8),f(i,9),f(i,10)
endfor

sed_pap,'../models/modcIII','Stage III',f1,iap
f=f1
for i=0,ninc-1 do begin
 printf,1,format='(12e12.4)',f(i,0),f(i,1),f(i,2),f(i,3),f(i,4),f(i,5),f(i,6),$
	f(i,7),f(i,8),f(i,9),f(i,10)
endfor

close,1

end
