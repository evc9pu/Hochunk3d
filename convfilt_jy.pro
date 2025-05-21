pro convfilt_jy,filtfile,skip,lam,f,fluxout,dbnu

readcol,filtfile,filtlam,filtf,skipline=skip
;convert phi_lam to phi_nu
nnu=n_elements(filtlam)
phi=fltarr(nnu)
nufilt=fltarr(nnu)
for i=0,nnu-1 do begin
  nufilt(nnu-i-1)=2.9979e14/filtlam(i)
;  phi(nnu-i-1)=filtf(i)
  if (dbnu eq 0) then phi(nnu-i-1)=filtf(i)
; according to Reach memo, should divide R by nu to get R_nu
; however, using Tom's files which took this into account, so setting
; dbnu=0 in input.
  if (dbnu eq 1) then phi(nnu-i-1)=filtf(i)/nufilt(nnu-i-1)
endfor

nu=2.9979e14/lam
fnew=interpol(f,nu,nufilt)
ftop=int_tabulated(nufilt,fnew*phi,/sort)
fbot=int_tabulated(nufilt,phi,/sort)
fluxout=ftop/fbot

;old version

; assumes lam is in microns, f in Jy
; converts f to ergs/s/cm^2/micron,
; convolves with filter,
; converts back to Jy

;f2=f/1.e23*2.9979e14/lam^2
;readcol,filtfile,filtlam,filtf,skipline=skip
;
;fnew=interpol(f2,lam,filtlam)
;ftop=int_tabulated(filtlam,fnew*filtf,/sort)
;fbot=int_tabulated(filtlam,filtf,/sort)
;fluxout=ftop/fbot

;ftop=int_tabulated(filtlam,filtlam*filtf,/sort)
;fbot=int_tabulated(filtlam,filtf,/sort)
;lameff=ftop/fbot
;print,'using filter effective lambda of ',lameff
;fluxout=fluxout*1.e23/2.9979e14*lameff^2

end
