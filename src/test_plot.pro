pro test_plot

  
  ;;;y=mx+c
  
  set_plot,'ps'
  device,filename='mx_c.eps',/encaps,xsize=30,ysize=30
  !p.thick=4
  !x.thick=4
  !y.thick=4
  !p.charsize=3
  !p.charthick=3
  
  readcol,'crf_binned.txt',min,max,mean,std,std_mean,skipline=1
  plot,(min+max)/2.,mean,psym=7,xrange=[0,65000],yrange=[0.95,1.05],/ystyle,xtitle='counts',ytitle='correction factors'
  errorbars,(min+max)/2.,mean,std
  linefit=fitaline((min+max)/2.,mean,std)
  oplot,[0,65000],linefit[0]*[0,65000]+linefit[2]
  device,/close
  set_plot,'x'

  cs=0                          ;chi-squared
  for i=0,n_elements(mean)-1 do begin
     cs+=((mean[i]-(linefit[0]*(min[i]+max[i])/2.+linefit[2]))/std[i])^2
  endfor
  dof=n_elements(mean)-2        ;degrees of freedom
  css=cs/dof                    ;chi-squared statistic
  print,css

  print,linefit[0],linefit[2]   ;y=linefit[0]*x+linefit[2]

  
  ;;;y=a*sqrt(x)
  
  set_plot,'ps'
  device,filename='sqrt_x.eps',/encaps,xsize=30,ysize=30
  !p.thick=4
  !x.thick=4
  !y.thick=4
  !p.charsize=3
  !p.charthick=3

  readcol,'crf_binned.txt',min,max,mean,std,std_mean,skipline=1
  plot,(min+max)/2.,mean,psym=7,xrange=[0,65000],yrange=[0.95,1.05],/ystyle,xtitle='counts',ytitle='correction factors'
  errorbars,(min+max)/2.,mean,std
  linefit=fitaline((min+max)/2.,(mean-1)^2,2*std*mean)
  xplot=dindgen(65000)
  oplot,xplot,sqrt(linefit[0]*xplot)+1
  device,/close
  set_plot,'x'

  cs=0
  for j=0,n_elements(mean)-1 do begin
     cs+=((mean[j]-(linefit[0]^2*(min[j]+max[j])/2+1))/std[j])^2
  endfor
  dof=n_elements(mean)-2
  css=cs/dof
  print,css

  print,linefit[0]              ;y=linefit[0]*sqrt(x)+1


  ;;;piecewise

  readcol,'crf_binned.txt',min,max,mean,std,std_mean,skipline=1

  x=(min+max)/2
  y=mean
  w_mean1=dblarr(65)
  w_mean2=dblarr(65)
  css=dblarr(65)
  index=0
  
  for bp=0,64000,1000 do begin
     w1=where(x lt bp)
     num1=0
     dem1=0
     for a=0,n_elements(w1)-1 do begin
        weight1=1/std[w1[a]]^2
        num1+=y[w1[a]]*weight1
        dem1+=weight1
     endfor
     w_mean1[index]=num1/dem1

     w2=where(x ge bp)
     num2=0
     dem2=0
     for b=0,n_elements(w2)-1 do begin
        weight2=1/std[w2[b]]^2
        num2+=y[w2[b]]*weight2
        dem2+=weight2
     endfor
     w_mean2[index]=num2/dem2

     cs=0
     dof=n_elements(mean)-2          ;we have 2 free parameters: the 2 constants that we're overplotting
     for c=0,n_elements(x)-1 do begin
        if (x[c] lt bp) then model=w_mean1[index]
        if (x[c] ge bp) then model=w_mean2[index]
        cs+=((y[c]-model)/std[c])^2
     endfor
     css[index]=cs/dof
     
     index+=1
  endfor

  set_plot,'ps'
  device,filename='min_breakpoint.eps',/encaps,xsize=30,ysize=30
  !p.thick=4
  !x.thick=4
  !y.thick=4
  !p.charsize=3
  !p.charthick=3

  break_point=(dindgen(65)+1)*1000
  plot,break_point,css,psym=7,xrange=[0,65000],yrange=[0,2],/ystyle,xtitle='break points',ytitle='reduced-chi squared'
  device,/close
  set_plot,'x'

  w_min=where(css eq min(css))       ;min(css) gives us the min reduced-chi squared, and the where() function saves the index of the min(css) in w_min
  bp_min=break_point[w_min]          ;saves the break point where min(css) occurs in bp_min
  mean1=w_mean1[w_min]
  mean2=w_mean2[w_min]
  
  set_plot,'ps'
  device,filename='piecewise.eps',/encaps,xsize=30,ysize=30
  !p.thick=4
  !x.thick=4
  !y.thick=4
  !p.charsize=3
  !p.charthick=3

  plot,(min+max)/2.,mean,psym=7,xrange=[0,65000],yrange=[0.95,1.05],/ystyle,xtitle='counts',ytitle='correction factors'
  errorbars,(min+max)/2.,mean,std
  oplot,[0,bp_min],[mean1,mean1]
  oplot,[bp_min,65000],[mean2,mean2]
  device,/close
  set_plot,'x'

  print,min(css)

  print,bp_min,mean1,mean2
  
end
