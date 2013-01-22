PRO source_counts_6C_catalog,hist_norm,flux_bin,binsize=binsize,min_use=min_use

filename='6c_full_radio_catalog.fit'
filename2=rootdir('mwa') + filename

data_arr=Ptrarr(4,/allocate)
*data_arr[0]=mrdfits(filename2,1,header0)
*data_arr[1]=mrdfits(filename2,2,header1)
*data_arr[2]=mrdfits(filename2,3,header2)
*data_arr[3]=mrdfits(filename2,4,header3)

ra_tag='_RAJ2000'
dec_tag='_DEJ2000'
flux_tag='F151P'

IF N_Elements(binsize) EQ 0 THEN binsize=0.15

flux_arr=Ptrarr(4,/allocate)
weight_arr=Ptrarr(4,/allocate)
area_arr=fltarr(4)
plots=ptrarr(4+1,/allocate)
names=strarr(4+1)
min_arr=fltarr(4)
max_arr=fltarr(4)

FOR i=0,3 DO BEGIN
    ra_use=(*data_arr[i]).(where(tag_names(*data_arr[i]) EQ ra_tag))
    dec_use=(*data_arr[i]).(where(tag_names(*data_arr[i]) EQ dec_tag))
    flux_use=(*data_arr[i]).(where(tag_names(*data_arr[i]) EQ flux_tag))
    
    ;the following is all for computing the effective area of the survey.
    ;as a check, should also look up the reported area!
    ra_center=Median(ra_use)
    dec_center=Median(dec_use)
    
    radec2xy,ra_use,dec_use,x_use,y_use,1
    dimension_use=Max(x_use)-Min(x_use)+1
    elements_use=Max(x_use)-Min(x_use)+1
    
    image_use=fltarr(dimension_use,elements_use)
    image_use[x_use,y_use]=flux_use
    weightgenerate,image_use,weights_use
    
    area_arr[i]=Total(weights_use)
;    *map_arr[i]=image_use
    *flux_arr[i]=flux_use
    *weight_arr[i]=weights_use
    
    IF Keyword_Set(min_use) THEN hist_I=histogram(Alog10(flux_use),locations=loc_I,binsize=binsize,min=min_use,omin=omin,omax=omax) $
        ELSE hist_I=histogram(Alog10(flux_use),locations=loc_I,binsize=binsize,omin=omin,omax=omax)
    min_arr[i]=omin
    max_arr[i]=omax
ENDFOR

print,area_arr
print,'Total area: ',Total(area_arr)
min_use=min(min_arr)
max_use=Max(max_arr)

FOR i=0,3 DO BEGIN
    flux_use=*flux_arr[i]
    hist_I=histogram(Alog10(flux_use),locations=loc_I,binsize=binsize,min=min_use,max=max_use)
    
    IF i EQ 0 THEN BEGIN
        nbins=N_Elements(hist_I)
        hist_tot=fltarr(nbins)
        hist_n=fltarr(nbins)
    ENDIF
    
    loc_I+=binsize/2.
    flux_I_bin=10.^loc_I
    binsize_use=(10.^(loc_I+binsize/2.)-10.^(loc_I-binsize/2.))
    
    hist_I_norm=hist_I/binsize_use/area_arr[i]
    
    ind_use=where(hist_I_norm,n_use2)
    
    hist_tot[ind_use]+=hist_I_norm[ind_use]
    hist_n[ind_use]+=1.
ENDFOR

ind_use=where(hist_n GT 0,n_use)
hist_norm=hist_tot[ind_use]/hist_n[ind_use]
flux_bin=flux_I_bin[ind_use]


END

;    
;    plot_single=fltarr(2,n_use2)
;    plot_single[0,*]=flux_I_bin[ind_use]
;    plot_single[1,*]=hist_I_norm[ind_use]
;    *plots[i]=plot_single
;    names[i]='6C field '+Strtrim(String(i+1),2)
;    min_arr[i]=omin
;    max_arr[i]=omax
;    
;    IF i EQ 0 THEN BEGIN
;        x_low=min(plot_single[0,*])
;        x_high=max(plot_single[0,*])
;        y_low=min(plot_single[1,*])
;        y_high=max(plot_single[1,*])
;    ENDIF ELSE BEGIN
;        x_low=x_low<min(plot_single[0,*])
;        x_high=x_high>max(plot_single[0,*])
;        y_low=y_low<min(plot_single[1,*])
;        y_high=y_high>max(plot_single[1,*])
;    ENDELSE
;ENDFOR
;
;
;
;
;hist_theory1=source_counts_theory(flux_I_bin,g1=1.75)
;hist_theory2=source_counts_theory(flux_I_bin,g1=1.18)
;
;loadct,39
;lstyle_list=[2,1,3,5,4]
;leg_x=400 ;360
;leg_y=460 ;480
;leg_y2=leg_y+4
;plot,flux_I_bin,hist_theory1,/xlog,/ylog,xtitle='Source flux (Jy)',ytitle='dN(S)/dS',$
;    xrange=xrange,yrange=yrange,background=!white,color=!black,xmargin=[12,1],ymargin=[6,1],thick=2,xcharsize=1.5,ycharsize=1.5
;plots,[leg_x-30,leg_x],[leg_y2,leg_y2],/device,color=!black,thick=2
;XYoutS,leg_x,leg_y,"  6C power law fit.",color=!black,/device,charsize=1.5;,charthick=2
;
;
;color_list=[!black,!red,!blue,!green,!magenta,!cyan,!yellow]
;FOR i=0,3 DO BEGIN
;
;    plot_x=Reform((*plots[i])[0,*])
;    plot_y=Reform((*plots[i])[1,*])
;    oplot,plot_x,plot_y,color=color_list[i+1],thick=2,linestyle=lstyle_list[i],psym=10
;    plots,[leg_x-30,leg_x],[leg_y2,leg_y2]-20*(i+1),/device,color=color_list[i+1],line=lstyle_list[i],thick=2;,psym=psym_style
;    XYoutS,leg_x,leg_y-20*(i+1),names[i],color=!black,/device,charsize=1.5
;ENDFOR
;
;plot_image=tvread(true=1)
;
;result_filepath=filepath('6C catalog source counts.png',Root_dir=rootdir('mwa'),subdir='temporary')
;write_png,result_filepath,plot_image,/verbose
;
;END
