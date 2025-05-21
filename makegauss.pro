pro makegauss,FWHM,k

n1=fix(6.*FWHM)
if n1 gt 145 then n1=145
sig=FWHM/2.354
sigsq=sig*sig
xc = rebin( findgen(n1) - n1/2., n1, n1)
yc = transpose(xc)
k = exp( -(xc^2 + yc^2) / (2.*sigsq) )

;im = convol(im,k,/edge_truncate)

tot=total(k)
k=k/tot

end
