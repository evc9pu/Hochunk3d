;calculate disk accretion rates based on star/disk parameters

;modify these
mstar=0.5
Rstar=2.09
rdisk=7.  ;inner disk radius in units of rstar
;rho0=3.474340464752512E-06
rho0=1.345368785627549E-06  ;get this from log file
h0=0.01
Tstar=4000.
alphad=0.01

;constants
msol=1.989d33
rsol=6.9598d10
G=6.6732d-8
Lsun=3.845e33
sigma=5.67e-5

mstar=mstar*msol
rstar=rstar*rsol
rdisk=rdisk*rstar
h0=h0*rstar
Lstar=4.d0*!pi*sigma*Tstar^4*Rstar^2

Vc=sqrt(G*Mstar/Rstar)

mdot=sqrt(18*!pi^3)*alphad*Vc*rho0*h0^3/Rstar

;Lacc=G*Mstar*Mdot/2.d0/Rstar
Lacc=G*Mstar*Mdot/2.d0/Rdisk
Lacc=Lacc/Lsun

mdot=mdot*3600.*24.*365.25/msol

print,'Lstar/Lsun',Lstar/Lsun
print,'Lacc/Lstar',Lacc
print,'mdot',mdot

end
