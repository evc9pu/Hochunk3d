pro makefitspol,dirm,d,rmax,nx,gsig,fileroot1,fileroot2,fileroot3

ilum=0  ;new version of code.  luminosity is already folded in
lum=1.
print,dirm,fileroot1

filename=dirm+fileroot1+'.dat'
readdat,filename,nx,im1
filename=dirm+fileroot2+'.dat'
readdat,filename,nx,im2
filename=dirm+fileroot3+'.dat'
readdat,filename,nx,im3
;im1=I, im2=Q, im3=U
im=(im2^2+im3^2)
im=sqrt(im)*im1  ;should be polarized flux.  the rest should be same as makefits...
;im=im1

;calculate ng, gaussian FWHM in pixels from gsig
;first compute pixel size 
rmaxi=rmax
pixsize=2.*rmaxi/float(nx)  ; in AU
ng=fix(gsig/pixsize)
if ng eq 0 then ng=1
print,'number of pixels for gaussian FWHM',ng

wvoigt=0.0  ;too tricky to explain how to use
ngwidth=30
if wvoigt lt 0.0001 then begin
print,'making gauss psf'
  makegauss,ng,img
endif else begin
  ;play with voigt profile
  img=psf_voigt_2d(ngwidth,wvoigt,ng)
endelse
print,'done making gauss psf'
;im=convol(im,img,/edge_truncate)
;this is much faster!
print,'convolving image'
help,im1,im2,im3,im,img
im=convolve(im,img,FT_PSF=psf_FT)
print,'done convolving image'

; don't convert anymore, already have units of intensity

mkhdr,head,im
sxaddpar,head,'PXSCAL1',pixsize,' AU','DATE'
sxaddpar,head,'PXSCAL2',pixsize,' AU','DATE'

writefits,fileroot1+'PF.fits', im, head

end
