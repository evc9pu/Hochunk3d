pro movie3col,dir,d,rmaxi,nx,gsig,Av,thetstr,phiarrst,filters,wavelengths,minvals,TEST=test
     
if keyword_set(TEST) then begin
   n=1
   print, 'in testing mode, making 1 image only'
endif else begin
   n=n_elements(phiarrst)
endelse

print,n

for i=0,n-1 do begin

  filerootr='e_'+filters(2)+'_'+thetstr+'_'+phiarrst(i)+'_I_img'
  filerootg='e_'+filters(1)+'_'+thetstr+'_'+phiarrst(i)+'_I_img'
  filerootb='e_'+filters(0)+'_'+thetstr+'_'+phiarrst(i)+'_I_img'
  
  print,filerootr
  print,filerootg
  print,filerootb
  
  makefits,dir,d,rmaxi,nx,gsig,filerootr,wavelengths(2),Av
  makefits,dir,d,rmaxi,nx,gsig,filerootg,wavelengths(1),Av
  makefits,dir,d,rmaxi,nx,gsig,filerootb,wavelengths(0),Av
  
  filenamer=filerootr+'.fits'
  filenameg=filerootg+'.fits'
  filenameb=filerootb+'.fits'

  imr=readfits(filenamer,head)
  img=readfits(filenameg,head)
  imb=readfits(filenameb,head)

  nx=fxpar(head,'NAXIS1')
  ny=fxpar(head,'NAXIS2')

  lrange=3
  valminr=alog10(minvals(2))
  valmaxr=valminr+lrange
  valming=alog10(minvals(1))
  valmaxg=valming+lrange
  valminb=alog10(minvals(0))
  valmaxb=valminb+lrange

  imnewr=bytscl(alog10(imr),min=valminr,max=valmaxr)
  imnewg=bytscl(alog10(img),min=valming,max=valmaxg)
  imnewb=bytscl(alog10(imb),min=valminb,max=valmaxb)

  imnew2r=congrid(imnewr,2*nx,2*ny,/interp)
  imnew2g=congrid(imnewg,2*nx,2*ny,/interp)
  imnew2b=congrid(imnewb,2*nx,2*ny,/interp)
  imnew3col=fltarr(3,2*nx,2*ny)
  imnew3col(2,*,*)=imnew2b
  imnew3col(1,*,*)=imnew2g
  imnew3col(0,*,*)=imnew2r
  window,0,xsize=6*nx,ysize=2*ny
  
  tv,imnew2b,0
  tv,imnew2g,1
  tv,imnew2r,2
  
  window,1,xsize=2*nx,ysize=2*ny
  tv,imnew3col,true=1
  wait,0.1

  fileroot='img_'+filters(0)+'_'+filters(1)+'_'+filters(2)
  write_png,fileroot+'.png',imnew3col

endfor

end
