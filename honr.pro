;calculate disk scale height at various radii

mstar=0.5  ;solar units
rstar=2. ; solar units
;rad=6.72   ;dust sublimation radius, units of rstar
;rad=1.
;T=1600.
Tstar=4000.
Tsub=1600.

k=1.38e-16
mu=2.3
mH=1.67e-24
msol=1.989d33
rsol=6.9598d10
G=6.6732d-8

;if you choose your favorite rad and T
;r=rad*rstar*rsol
;M=mstar*msol
;honr=sqrt(k*T/G/M*r/mu/mH)
;h=honr*rad
;print,'h on r ',honr
;print,'h (r_sub) ',h
;h0=honr*rad/rad^1.25
;print,'h_0 ',h0

;empirical fit to dust destruction radius
rsub=(Tsub/Tstar)^(-2.1)
print,''
print,'r_sub based on empirical fit ',rsub
T=Tsub
rad=rsub
r=rad*rstar*rsol
M=mstar*msol
honr=sqrt(k*T/G/M*r/mu/mH)
h=honr*rad
print,'h on r at r_sub',honr
print,'h (r_sub) ',h
b=1.25
h0=honr*rad/rad^b
print,'h_0 at rstar based on h(r_sub)',h0
print,'assumes Tsub=',Tsub,'    beta =',b

end
