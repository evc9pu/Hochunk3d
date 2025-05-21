

;uncomment filename and extrapolation exponent and comment the others
filnam='kmh'
;filnam='r400_ice095'
;filnam='ww04'   ;note, original ww04 seems fine, goes out to 4.2mm, why, I don't know.
;filnam='www003'
;filnam='www006'
exp=-2.  ;kmh and r400_ice095
;exp=-1.57  ; ww04
;exp=-1.3   ;www003
;exp=-0.6  ;www006

readcol,filnam+'.par',lam,cext,csca,kap,hgg,pmax,thet
;thet is garbage

window,0,retain=2
plot,lam,kap,/xlog,/ylog,xrange=[1.,50000],yrange=[0.00001,100]

;extrapolate to 4.2 mm
a=where(lam gt 500)
i=a(0)
kap42mm=kap(i)*(4.2e4/lam(i))^exp

lam=[lam,4.2e4]
kap=[kap,kap42mm]
n=n_elements(cext)
cext=[cext,cext(n-1)]
csca=[csca,csca(n-1)]
hgg=[hgg,hgg(n-1)]
pmax=[pmax,pmax(n-1)]
thet=[thet,thet(n-1)]

oplot,lam,kap,linestyle=1

close,1
openw,1,filnam+'_extrap.par'
n=n_elements(lam)
for i=0,n-1 do begin
printf,1,lam(i),cext(i),csca(i),kap(i),hgg(i),pmax(i),thet(i),$
	format='(7(1x,e11.4))'
endfor
close,1

end