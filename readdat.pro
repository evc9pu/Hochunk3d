pro readdat,filename,nx,out

out=fltarr(nx,nx)

close,1
openr,1,filename
readf,1,out

end
