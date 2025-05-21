pro makefits,dirm,d,rmax,nx,gsig,fileroot,wave,Av

ilum=0  ;new version of code.  luminosity is already folded in
lum=1.
print,dirm,fileroot

filename=dirm+fileroot+'.dat'

;calculate ng, gaussian FWHM in pixels from gsig
;first compute pixel size 
rmaxi=rmax
pixsize=2.*rmaxi/float(nx)   ; in AU
ng=fix(gsig/pixsize)
if ng eq 0 then ng=1
print,'number of pixels for gaussian FWHM',ng

readdat,filename,nx,im
wvoigt=0.0  ;too tricky to explain how to use
ngwidth=30
if ng gt 1 then begin
  if wvoigt lt 0.0001 then begin
  print,'making gauss psf'
    makegauss,ng,img
  endif else begin
    ;play with voigt profile
    img=psf_voigt_2d(ngwidth,wvoigt,ng)
  endelse
  ;im=convol(im,img,/edge_truncate)
  ;this is much faster!
  im=convolve(im,img,FT_PSF=psf_FT)
endif

;add in foreground reddening
if Av gt 0.0 then begin
  readcol,'../models/parfiles/kmh.par',klam,junk1,junk2,kkap
  kapnew=interpol(kkap,klam,wave)
  kapv=interpol(kkap,klam,0.55)
  taul=kapnew/kapv/1.086*Av
  extinct=exp(-taul)
  im=im*extinct
endif

mkhdr,head,im
sxaddpar,head,'PXSCAL1',pixsize,' AU','DATE'
sxaddpar,head,'PXSCAL2',pixsize,' AU','DATE'

writefits,fileroot+'.fits', im, head

print, 'max image ',max(im)

end
