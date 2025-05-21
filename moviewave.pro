pro moviewave,dir,thetstr,phistr,minval,TEST=test

;minval=1000000.
 
fileroot=dir+'e_cube_'+thetstr+'_'+phistr+'_I_img'
filename=dir+'e_cube_'+thetstr+'_'+phistr+'_I_img.fits.gz'
img=mrdfits(filename,0,header)
wave=mrdfits(filename,1,headwave)
n=n_elements(wave)
nx=fxpar(header,'NAXIS1')
ny=fxpar(header,'NAXIS2')   
im=fltarr(nx,ny)     
        
if keyword_set(TEST) then begin
   i1=96
   i2=96
   print, 'in testing mode, making 1 image only (1 micron)'
endif else begin
   i1=0
   i2=n-1-40
endelse

print,n,nx,ny

lrange=3
valmin=alog10(minval)
valmax=valmin+lrange

window,1,xsize=2*nx,ysize=2*ny

for i=i2,i1,-1 do begin
  
  print,n-i,wave(i).wavelength
  
  im=img(*,*,i)
  imnew=bytscl(alog10(im),min=valmin,max=valmax)

  imnew2=congrid(imnew,2*nx,2*ny,/interp)
  tv,imnew2
  wait,0.1
 ;  stop

  write_png,'wave'+strcompress(string(n-i),/remove_all)+'.png',imnew2

endfor

end
