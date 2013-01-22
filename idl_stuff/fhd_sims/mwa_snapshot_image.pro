PRO mwa_snapshot_image,polarization=polarization,time_snap=time_snap,$
    data_directory=data_directory,filename=filename,restore_clean=restore_clean,$
    image_create=image_create,zoom_image=zoom_image,beam_threshold=beam_threshold,$
    max_image=max_image
;snapshot data must have been gridded previously, and the transfer functions generated
heap_gc
IF N_Elements(polarization) EQ 0 THEN polarization=0
IF N_Elements(time_snap) EQ 0 THEN time_snap=0

COMMON obs_params,dimension,elements,degpix,kbinsize,azimuth_start,elevation_start,rotation,JD0
IF N_Elements(dimension) EQ 0 THEN dimension=1024.
IF N_Elements(elements) EQ 0 THEN elements=1024.
IF N_Elements(baseline_threshold) EQ 0 THEN baseline_threshold=50.
IF N_Elements(beam_threshold) EQ 0 THEN beam_threshold=0.1

COMMON data_params,flag_arr,vis_xx,vis_yy,vis_xy,vis_yx,file_path
COMMON visibility_params,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,freq_bin_i,file_path2
COMMON psf_params,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals

IF N_Elements(data_directory) EQ 0 THEN BEGIN
    ;data_directory='DATA\r4\clmw\X13\HydA_121_20100322135442'
    ;data_directory='DATA\r4\clmw\X14\CenA_20100921082132'
    ;;data_directory='DATA\r4\clmw\X14\PicA_121_20100924211938'
   data_directory=''
ENDIF
IF N_Elements(filename) EQ 0 THEN BEGIN
    ;filename='HydA_121_20100322135442'
    ;filename='CenA'
    filename='PicA_121_20100924211938.cal'
    ;filename='PicA_121_20100924211938'
ENDIF
ext='.UVFITS'
file_path=filepath(filename,root_dir=rootdir('mwa'),subdir=data_directory)
file_path2=file_path

pol_names=['xx','yy','xy','yx']
time_name=(Sindgen2(100,char=2))[time_snap] ;inelegant, but it works
;different frequency intervals can be supported later
freq_name='00'
file_identifier=String(format='(A2,A2,"_",A2)',pol_names[polarization],time_name,freq_name)
;file_name_base=String(format='("_xferfn_",A2,A2,"_",A2)',pol_names[pol_i],time_name,freq_names[freq_i])

IF Keyword_Set(restore_clean) THEN restore,file_path+'_clean_'+file_identifier+'.sav' ELSE BEGIN
    beam_setup,data_directory=data_directory,filename=filename,/restore_last
    restore,filename=file_path+'_uv_'+file_identifier+'.sav' ; dirty_uv,weights_grid
    ;restore,filename=file_path+'_dirty_'+file_identifier+'.sav' ; dirty_image
    restore,file_path+'_xferfn_'+file_identifier+'.sav' ;xfer_fn
    max_iter=5000.
    gain_factor=0.1
    center_gain=0.95
    transfer_interval=20.
    transfer_threshold=0.85
    visibility_image_clean,dirty_uv,clean_image,xfer_fn,source_array=source_array,beam_threshold=beam_threshold,$
        model_uv=model_uv,max_iter=max_iter,gain_factor=gain_factor,center_gain=center_gain,$
        transfer_interval=transfer_interval,transfer_threshold=transfer_threshold,timing=t_clean,$
        beam_base=beam_base,polarization=polarization,baseline_threshold=baseline_threshold
    Ptr_free,psf_base,psf_residuals_i,psf_residuals_val,psf_xvals,psf_yvals
    psf_residuals_n=0
    save,clean_image,model_uv,source_array,beam_base,filename=file_path+'_clean_'+file_identifier+'.sav'
ENDELSE

IF restore_clean EQ 2 THEN restore,filename=file_path+'_clean2_'+file_identifier+'.sav' ELSE BEGIN
    IF restore_clean EQ 1 THEN BEGIN
        restore,filename=file_path+'_uv_'+file_identifier+'.sav' ; dirty_uv,weights_grid
        ;restore,filename=file_path+'_dirty_'+file_identifier+'.sav' ; dirty_image
        restore,file_path+'_xferfn_'+file_identifier+'.sav' ;xfer_fn    
    ENDIF
    weights=abs(visibility_transfer_apply(complexarr(dimension,elements)+complex(1,0),xfer_fn,/new))
    normalization=1./Mean(weights[where(weights)])
    n_source0=(size(source_array,/dimension))[1]
    source_array1=visibility_source_extract(clean_image,beam_use=beam_base,max_sources=1000.,timing=t_extract,beam_threshold=beam_threshold)
    n_source1=(size(source_array1,/dimension))[1]
    source_array2=fltarr(4,n_source0+n_source1)
    source_array2[*,0:n_source0-1]=source_array
    source_array2[*,n_source0:*]=source_array1
    source_array2=clean_source_condense(source_array2,radius=0.5)
    normalization2=source_array2[2,0]/source_array2[3,0]
    n_sources=(size(source_array2,/dimension))[1]
    
    ;columns of source_array are: 0:x, 1:y, 2:RA, 3:Dec, 4:co-added flux, 5:pixel index, 6:varies
    source_array=fltarr(7,n_sources)
    i1=[0,1,4] & i2=[0,1,3]
    FOR i=0,2 DO source_array[i1[i],*]=source_array2[i2[i],*]
    source_array[0,*]+=dimension/2.
    source_array[1,*]+=elements/2.
    source_array[5,*]=Round(source_array[0,*])+Round(source_array[1,*])*dimension
    
    model_uv2=visibility_source_uv_grid(source_array)*normalization2 ; off by a normalization factor!
    dirty_image=dirty_image_generate(dirty_uv,baseline_threshold=baseline_threshold,normalization=normalization)
    model_image=dirty_image_generate(visibility_transfer_apply(model_uv2,xfer_fn,/new),baseline_threshold=baseline_threshold,normalization=normalization)
    clean_image2=dirty_image-model_image
    
    source_image=clean_source_restore(fltarr(dimension,elements),source_array)
    
    save,clean_image2,model_uv2,source_array2,source_image,beam_base,filename=file_path+'_clean2_'+file_identifier+'.sav'
ENDELSE
stop
IF Keyword_Set(image_create) THEN BEGIN
    image_region=where(beam_base GT 0.1,n_image,complement=image_cut,ncomplement=n_cut)
    clean_residual=clean_image/beam_base
    clean_residual2=clean_image2/beam_base
    source_imageB=source_image/beam_base
    source_imageB[where(source_imageB)]=Alog10(source_imageB[where(source_imageB)])

    recon_image=clean_residual2+source_imageB
    IF n_cut GT 0 THEN clean_residual[image_cut]=0
    
    IF Keyword_Set(max_image) THEN high=max_image ELSE high=Max(clean_residual[image_region])
    low=-high    
    IF Keyword_Set(max_image) THEN high2=max_image ELSE high2=Max(clean_residual2[image_region])
    low2=-high2
    high3=Max(source_imageB[image_region])
    low3=0
    region_xy=array_indices(clean_residual,image_region)
    zoom_x_low=min(region_xy[0,*])
    zoom_x_high=max(region_xy[0,*])
    zoom_y_low=min(region_xy[1,*])
    zoom_y_high=max(region_xy[1,*])
    IF Keyword_Set(zoom_image) THEN BEGIN
        image_out=clean_residual[zoom_x_low:zoom_x_high,zoom_y_low:zoom_y_high]
        Imagefast,image_out,data_dir=data_directory+'\Images',filename=filename+'_cleanresidual_'+file_identifier,$
            /right_colorbar,/zero_black,project='mwa',high=high,low=low
        image_out=clean_residual2[zoom_x_low:zoom_x_high,zoom_y_low:zoom_y_high]
        Imagefast,image_out,data_dir=data_directory+'\Images',filename=filename+'_cleanresidual2_'+file_identifier,$
            /right_colorbar,/zero_black,project='mwa',high=high2,low=low2
        image_out=source_imageB[zoom_x_low:zoom_x_high,zoom_y_low:zoom_y_high]
        Imagefast,image_out,data_dir=data_directory+'\Images',filename=filename+'_sources_'+file_identifier,$
            /right_colorbar,/zero_black,project='mwa',high=high3,low=low3
        
    ENDIF ELSE BEGIN
        image_out=clean_residual
        Imagefast,image_out,data_dir=data_directory+'\Images',filename=filename+'_cleanresidual_'+file_identifier,$
            /right_colorbar,/zero_black,project='mwa',high=high,low=low
        image_out=clean_residual2
        Imagefast,image_out,data_dir=data_directory+'\Images',filename=filename+'_cleanresidual2_'+file_identifier,$
            /right_colorbar,/zero_black,project='mwa',high=high2,low=low2
        image_out=source_imageB
        Imagefast,image_out,data_dir=data_directory+'\Images',filename=filename+'_sources_'+file_identifier,$
            /right_colorbar,/zero_black,project='mwa',high=high3,low=low3
    ENDELSE
    
ENDIF

;field_name=Strmid(filename,0,4)
;parkes_catalog_read,source_array_anc,field_name=field_name,x_mirror=1,y_mirror=1,degpix=2.*!RaDeg/1024./0.9167
;base_color=source_array[4,0]/source_array_anc[4,0]
;source_array_anc[4,*]*=base_color
;source_array_anc_test=source_array_anc
;source_array_anc_test[0,*]+=dimension/2.-(source_array_anc_test[0,0])
;source_array_anc_test[1,*]+=elements/2.-(source_array_anc_test[1,0])
;anc_image=clean_source_restore(fltarr(dimension,elements),source_array_anc_test)
;
;mwa_align_pipe,header,clean_image,match_test,Ancillary_image_use,$
;    error=error,source_array=source_array,ancillary_source_array=source_array_anc
END
