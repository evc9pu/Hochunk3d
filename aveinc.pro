pro aveinc,dirm,iapin,iincin,ns

;iinc ranges from 1-10 unless you change ttsre program
;iapin ranges from 1-5

ninc=10l  ;change this if you change ttsre program
nap=50l  ; "
nfreq=250l   ; "

iap=iapin-1l 
iinc=iincin-1l

print,'skipping',iap*nfreq*ninc+iinc*nfreq+1l,' lines'
readcol,dirm+'/aveinc.dat',x4,f4,skipline=iap*nfreq*ninc+iinc*nfreq+1l,numline=nfreq,format='(d,d)'

;readcol,dirm+'/flux.dat',x4,f4,skipline=nap*nfreq*ninc+1l,numline=nfreq*ninc,format='(d,d)'

mu=float(iinc)/float(ninc)+0.05
;if iinc eq 0 then mu=0.
print,'mu',mu,'theta',acos(mu)*180./!pi
thet=acos(mu)*180./!pi

;!p.multi=[0,1,1]
;set_plot,'ps'
;if !d.name eq 'PS' then device,xsize=6.,ysize=4,/inches,$
;/portrait,xoffset=0.2,yoffset=2,bits=24,filename='aveinc.ps',/color

yr1=acos(mu+0.05)*180./!pi
yr2=acos(mu-0.05)*180./!pi

;if iinc eq 0 then begin
;  yr2=acos(0.05)*180./!pi
;  yr1=90.
;endif

;print,f4

if ns gt 0 then begin
   a=where(f4 lt 89.9)
   result=moment(f4(a))
   print,'ave viewing angle ',result(0)
   f4smooth=f4
   f4smooth(a)=smooth(f4(a),ns)
   print,'smoothing by ',ns
endif else begin
   f4smooth=f4
endelse


plot,x4,f4smooth,title='theta='+strcompress(thet),xrange=[0.5,1000],/xlog,$
   yrange=[yr1-.05*yr1,yr2+0.05*yr2],charsize=1.5,max_value=89.99
oplot,[0.01,10000.],[thet,thet],linestyle=2
oplot,[0.01,10000.],[yr1,yr1],linestyle=1
oplot,[0.01,10000.],[yr2,yr2],linestyle=1

;device,/close

end
