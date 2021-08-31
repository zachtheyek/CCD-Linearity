pro linearity,filelist

  nfiles=n_elements(filelist)
  nx=3358
  ny=2536
  bp=28000

  for i=0,nfiles-1 do begin
     fits_read,filelist[i],image,header
     
     w1=where(image lt bp)
     image[w1]*=1.0014226
     
     w2=where(image ge bp)
     image[w2]*=1.0180221

     writefits,filelist[i]+'_linearized.fits',image,header
  endfor

end
