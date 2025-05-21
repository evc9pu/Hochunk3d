pro lum,rstar,tstar
;this may bomb if you don't use real numbers, e.g.,
; lum,2.0,4000.0

c=2.9979d14
;d=140.
sigma=5.67d-5
;pc=3.0857e18
lsun=3.845d33
l=4.d0*!pi*(rstar*6.955d10)^2*sigma*tstar^2*tstar^2
l=l/lsun
print,l

end