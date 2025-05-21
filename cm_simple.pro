;cc.pro

;plotsed has an option to calculate filter fluxes from SED models.

!p.multi=[0,1,1]
loadct,13
set_plot,'ps'
device,xsize=5.,ysize=4.,/inches,$
/portrait,xoffset=1,yoffset=1,/color,bits=24,$
;filename='cc.ps'
filename='cc.eps',/encapsulated


;for fitsio version of code
;readcol,'filterfluxes.txt',cosi,phi,j,h,k,b1,b2,b3,b4,b5,b6,b7,skipline=2
readcol,'filter_AeBe_disk_nopah.txt',cosi,phi,j,h,k,b1,b2,b3,b4,b5,b6,b7,skipline=2
;readcol,'filter_AeBe_env.txt',cosi,phi,j,h,k,b1,b2,b3,b4,b5,b6,b7,skipline=2
;for non-fitsio version of code
;readcol,'filterfluxes.txt',cosi,j,h,k,b1,b2,b3,b4,b5,b6,b7,skipline=2
ninc=n_elements(j)
nfilt=10


;Jys at 0 mag at each wavelength
filtfb=[1594.,1024.,666.7,280.9,179.7,115.0,64.13,7.160,0.7972,0.160]
;when ready, pacs/spires are ,0.6751077,0.3107111,0.1281016,5.8182828E-02,2.7805042E-02,1.3336738E-02]

;convert everything to magnitudes
;used d=10pc so these are absolute mags
j=-2.5*alog10(j/filtfb(0))
h=-2.5*alog10(h/filtfb(1))
k=-2.5*alog10(k/filtfb(2))
b1=-2.5*alog10(b1/filtfb(3))
b2=-2.5*alog10(b2/filtfb(4))
b3=-2.5*alog10(b3/filtfb(5))
b4=-2.5*alog10(b4/filtfb(6))
b5=-2.5*alog10(b5/filtfb(7))
b6=-2.5*alog10(B6/filtfb(8))
b7=-2.5*alog10(b7/filtfb(9))

;pick a distance
d=140.
j=j-5+5*alog10(d)
h=h-5+5*alog10(d)
k=k-5+5*alog10(d)
b1=b1-5+5*alog10(d)
b2=b2-5+5*alog10(d)
b3=b3-5+5*alog10(d)
b4=b4-5+5*alog10(d)
b5=b5-5+5*alog10(d)
b6=b6-5+5*alog10(d)
b7=b7-5+5*alog10(d)

;here's a bunch of sample plots, pick one or make your own.

;format is:  plotcol_simp,d,x,y,xlab,ylab,xmin,xmax,ymin,ymax,ninc,op,ccol   
; x, y are the colors you want, e.g., h-k, j-h
; xlab, ylab are axis labels
; xmin, xmax, ymin, ymax are plotting ranges
; ninc is number of inclinations
; op is 0 if you are plotting one set for the first time;  if overplotting, set to 1
; ccol is color number
; circ = 'fill' or 'open', symbol size

ccol=55
circ='open'

;plotcol_simp,d,h-k,j-h,cosi,'H-K','J-H',-.3,2.,-.5,4.2,ninc,0,ccol,circ

;plotcol_simp,d,k-b1,j-h,cosi,'K-[3.6]','J-H',-.5,3.5,0,3.5,ninc,0,ccol,circ

plotcol_simp,d,b3-b4,b1-b2,cosi,'[5.8]-[8.0]','[3.6]-[4.5]',0.,3.5,0.,1.5,ninc,0,ccol,circ

;plotcol_simp,d,b4-b5,b1-b3,cosi,'[8.0]-[24]','[3.6]-[5.8]',-1.,10,-1,7,ninc,0,ccol,circ

;plotcol_simp,d,b5-b6,b4-b5,cosi,'[24]-[70]','[8.0]-[24]',0.0,15.0,0,10,ninc,0,ccol,circ

; I'm sorry, you're just going to have to repeat most of these lines to make overplots:


readcol,'filter_AeBe_disk_pah.txt',cosi,phi,j,h,k,b1,b2,b3,b4,b5,b6,b7,skipline=2
;readcol,'filter_AeBe_env_pah.txt',cosi,phi,j,h,k,b1,b2,b3,b4,b5,b6,b7,skipline=2

ninc=n_elements(j)

;convert everything to magnitudes
;used d=10pc so these are absolute mags
j=-2.5*alog10(j/filtfb(0))
h=-2.5*alog10(h/filtfb(1))
k=-2.5*alog10(k/filtfb(2))
b1=-2.5*alog10(b1/filtfb(3))
b2=-2.5*alog10(b2/filtfb(4))
b3=-2.5*alog10(b3/filtfb(5))
b4=-2.5*alog10(b4/filtfb(6))
b5=-2.5*alog10(b5/filtfb(7))
b6=-2.5*alog10(B6/filtfb(8))
b7=-2.5*alog10(b7/filtfb(9))

j=j-5+5*alog10(d)
h=h-5+5*alog10(d)
k=k-5+5*alog10(d)
b1=b1-5+5*alog10(d)
b2=b2-5+5*alog10(d)
b3=b3-5+5*alog10(d)
b4=b4-5+5*alog10(d)
b5=b5-5+5*alog10(d)
b6=b6-5+5*alog10(d)
b7=b7-5+5*alog10(d)

ccol=155
circ='open'
plotcol_simp,d,b3-b4,b1-b2,cosi,'[5.8]-[8.0]','[3.6]-[4.5]',-.5,2.0,-.5,4,ninc,1,ccol,circ


device,/close
end
