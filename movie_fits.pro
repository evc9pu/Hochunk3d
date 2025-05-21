pro movie_fits,dir,d,rmaxi,nx,gsig,Av,thetstr,phiarrst,filter,wavelength,minval,TEST=test

;minval=1000000.
     
if keyword_set(TEST) then begin
   n=1
   print, 'in testing mode, making 1 image only'
endif else begin
   n=n_elements(phiarrst)
endelse

print,n

for i=0,n-1 do begin

  fileroot='e_'+filter+'_'+thetstr+'_'+phiarrst(i)+'_I_img'
    
  filename=dir+fileroot+'.fits.gz'

  im=readfits(filename,head)

  nx=fxpar(head,'NAXIS1')
  ny=fxpar(head,'NAXIS2')
  pixsizx=sxpar(head,'PXSCAL1')
  pixsizy=sxpar(head,'PXSCAL2')
  
  print,min(im),max(im)
  lrange=3
  valmin=alog10(minval)
  valmax=valmin+lrange

  imnew=bytscl(alog10(im),min=valmin,max=valmax)

  imnew2=congrid(imnew,2*nx,2*ny,/interp)
  tv,imnew2
  wait,0.2

  write_png,fileroot+'.png',imnew2

endfor

end
