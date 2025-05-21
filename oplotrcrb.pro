pro oplotrcrb

readcol,'rcrb-sed-barb.txt',lam1,f1,ferr1

readcol,'rcrb3.sws',lam2,f2,f2err


;convert to mJy
f1=f1*1000.
f2=f2*1000.

oplot,lam1,f1,psym=5,symsize=0.2
oplot,lam2,f2

end