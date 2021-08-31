pro plot_fits,makeplot_alldata,pass  
     
  nx=3358
  ny=2536
  bin_step=1000.
  nbins=floor(65000./bin_step)
  exptime=['0.500','0.750','1.000','1.250','1.500','2.000','4.000','7.000','10.000','13.000','15.000','16.000','17.000','18.000','19.000']
  nexp=n_elements(exptime)
     
  fits_read,'bias_master.fits',bias,header_bias
  fits_read,'bias_uncertainty.fits',unc_bias,header_unc_bias
     
  lin_sub_bias=dblarr(nx,ny,nexp)
  lin_unc=dblarr(nx,ny,nexp)
  for i=0,nexp-1 do begin
     fits_read,'linearity.'+exptime[i]+'_master.fits',linearity,header_linearity
     fits_read,'linearity.'+exptime[i]+'_uncertainty.fits',unc_linearity,header_unc_linearity
     lin_sub_bias[*,*,i]=linearity-bias
     lin_unc[*,*,i]=unc_linearity
  endfor
     
     
  if(pass eq 1) then begin
     bin_min=dindgen(nbins)*bin_step
     bin_max=bin_min+bin_step
     bin_sum=dblarr(nbins)
     bin_num=lonarr(nbins)
  endif
  if(pass eq 2) then begin
     readcol,'crf_bin_pass1.txt',bin_min,bin_max,bin_sum,bin_num,bin_mean
     nbins=n_elements(bin_min)
     bin_std=dblarr(nbins)
 endif
 
  if(makeplot_alldata eq 'yes') then begin
     window,0,xsize=750,ysize=750
     !p.thick=2
     !x.thick=2
     !y.thick=2
     !p.charsize=1.5
     !p.charthick=1.5
  endif

  for j=0,nx-1 do begin

     print,'Starting ',j,' of ',nx-1
     
     for k=0,ny-1 do begin

        lincounts=reform(lin_sub_bias[j,k,*])
        unc_lin=sqrt((reform(lin_unc[j,k,*])^2.)+(unc_bias[j,k]^2.)) ;add bias and linearity uncertainties in quadrature
        linefit=fitaline(exptime[0:4],lincounts[0:4],unc_lin[0:4])

        if((j eq 0) and (k eq 0)) then begin
           set_plot,'ps'
           device,filename='exptime_counts_plot.eps',/encaps,xsize=30,ysize=30
           !p.thick=4
           !x.thick=4
           !y.thick=4
           !p.charsize=3
           !p.charthick=3
           
           plot,exptime,lincounts,psym=7,/ystyle,xtitle='exposure time (s)',ytitle='counts'
           errorbars,exptime,lincounts,unc_lin
           oplot,[0,60],linefit[0]*[0,60]+linefit[2]
           device,/close
           set_plot,'x'
        endif 

        expcounts=linefit[0]*exptime+linefit[2]  ;use the best-fit line to obtain expected counts at each exposure time
        unc_exp=linefit[1]*exptime+linefit[3]    ;use uncertainty in best-fit line to get uncertainties in expected counts
  
        crf=expcounts/lincounts                              ;calculate correction factor = expected counts / measured counts
        unc_crf=crf*(unc_exp/expcounts+unc_lin/lincounts)    ;use propagation of errors to calculate correction factor uncertainties

        for z=0,nbins-1 do begin

           w=where((lincounts ge bin_min[z]) and (lincounts lt bin_max[z]))

           if(w[0] ne -1) then begin

              if(pass eq 1) then begin
                 bin_sum[z]=bin_sum[z]+total(crf[w])
                 bin_num[z]=bin_num[z]+n_elements(w)
              endif
             if(pass eq 2) then begin
                for index=0,n_elements(w)-1 do begin
                   bin_std[z]=bin_std[z]+((crf[w[index]]-bin_mean[z])^2.0)
                endfor
             endif
          endif
        endfor
       
        if(makeplot_alldata eq 'yes') then begin
           if((j eq 0) and (k eq 0)) then begin
              plot,lincounts,crf,psym=7,yrange=[0.9,1.1],/ystyle,xtitle='counts',ytitle='correction factor'
           endif else begin
              oplot,lincounts,crf,psym=7
              errorbars,lincounts,crf,unc_crf,color=150
           endelse
        endif       
     endfor
  endfor

  if(pass eq 1) then begin
     openw,1,'crf_bin_pass1.txt'
     for z=0,nbins-1 do begin
        printf,1,bin_min[z],bin_max[z],bin_sum[z],bin_num[z],bin_sum[z]/bin_num[z]
     endfor
     close,1
  endif

  if(pass eq 2) then begin
     openw,1,'crf_binned.txt'
     for z=0,nbins-1 do begin
        printf,1,bin_min[z],bin_max[z],bin_sum[z]/bin_num[z],sqrt(bin_std[z]/bin_num[z]),sqrt(bin_std[z]/bin_num[z])/sqrt(bin_num[z])
     endfor
     close,1
     
     set_plot,'ps'
     device,filename='counts_crf_plot.eps',/encaps,xsize=30,ysize=30
     !p.thick=4
     !x.thick=4
     !y.thick=4
     !p.charsize=3
     !p.charthick=3
     
     readcol,'crf_binned.txt',min,max,mean,std,std_mean
     plot,(min+max)/2.,mean,psym=7,xrange=[0,65000],yrange=[0.95,1.05],/ystyle,xtitle='counts',ytitle='correction factors'
     errorbars,(min+max)/2.,mean,std
     device,/close
     set_plot,'x'     
  endif
  
  stop
end
