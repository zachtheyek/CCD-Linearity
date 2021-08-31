pro make_master,name

  nx=3358
  ny=2536
  readcol,'filelist_'+name+'.txt',file,format='a'
  nfiles=n_elements(file)

  all_images=dblarr(nx,ny,nfiles)

  for i=0,nfiles-1 do begin
     fits_read,file[i],image,header
     all_images[*,*,i]=image
  endfor

  image_mean=dblarr(nx,ny)
  image_unc=dblarr(nx,ny)
  for i=0,nx-1 do begin
     for j=0,ny-1 do begin
        image_mean[i,j]=mean(all_images[i,j,*])
        image_unc[i,j]=stddev(all_images[i,j,*])/sqrt(1.0*nfiles)     ;stddev of the mean
     endfor
  endfor

  writefits,'linearity.'+name+'_master.fits',image_mean
  writefits,'linearity.'+name+'_uncertainty.fits',image_unc

end
