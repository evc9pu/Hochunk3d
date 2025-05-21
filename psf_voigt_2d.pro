function psf_voigt_2d, npix, alpha, u, cntrd=cntrd

;+
; Crude routine to generate a 2D (normalized) psf using the voigt line profile.
; PSF is assume to be aziumthally symmetric.
;
; INPUTS:
; npix = number of pixels in each of the x and y directions
; alpha = value of the "a" or alpha parameters for call to VOIGT (see
;		IDL documentation).
; OUTPUT:
; psf = 2D array with normalized psf
;
; EXAMPLE:
; generate a 101x101 psf with alpha=5.
;
; IDL> psf = psf_voigt_2d(101,5.)
;
; history:
; 2003/11/02 (mjw):  created
;-

if N_ELEMENTS(npix) ne 1 then begin
     message,'usage:  psf = psf_voigt_2d(npix,alpha,u)'
     retall
endif


; set up psf array
psf = fltarr(npix,npix)

; deal with centroid
if N_elements( cntrd ) LE 0 then begin
   cntrd=(npix-1)/2. 
endif else begin
   cntrd=cntrd-0.5
endelse
if N_elements( alpha ) LE 0 then alpha=0.1*npix

; construct "centered" x and y arrays
x = findgen(npix) - cntrd
y = findgen(npix) - cntrd

; generate PSF function for each pixel, using radius as the dimensional
; input to VOIGT
print,u
for i=0,npix-1 do begin
   for j=0,npix-1 do begin
      r = (sqrt(x[i]^2 + y[j]^2))/u*1.67   ;1.67 gives the right FWHM factor
                                          ;don't understand. ????
      psf[i,j] = voigt(alpha,r)
   endfor
endfor

; normalize integral of PSF to 1.
norm = total(psf)
psf = psf / norm

; Fichons!
return,psf
end
