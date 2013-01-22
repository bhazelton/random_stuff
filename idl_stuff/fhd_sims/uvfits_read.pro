PRO uvfits_read,dirty_image,dirty_image_UV,weights_UV,test=test,$
    calibrate=calibrate,baseline_threshold=baseline_threshold,$
    restore_last=restore_last,no_save=no_save, restore_beam = restore_beam
IF N_Elements(test) EQ 0 THEN test=0
;IF N_Elements(calibrate) EQ 0 THEN calibrate=1
IF N_Elements(baseline_threshold) EQ 0 THEN baseline_threshold=50.
heap_gc

COMMON obs_params,dimension,elements,degpix,kbinsize,azimuth_start,elevation_start,rotation,JD0
IF N_Elements(dimension) EQ 0 THEN dimension=1024.
IF N_Elements(elements) EQ 0 THEN elements=1024.

COMMON data_params,flag_arr,vis_xx,vis_yy,vis_xy,vis_yx,file_path
COMMON visibility_params,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,freq_bin_i,file_path2
COMMON psf_params,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals
;leave both hard-coded for now
;data_directory='DATA\r4\clmw\X13\HydA_121_20100322135442'
;filename='HydA_121_20100322135442'

;data_directory='DATA\r4\clmw\X14\CenA_20100921082132'
;filename='CenA'
;data_directory='DATA\r4\clmw\X14\PicA_121_20100924211938'
data_directory=''
filename='PicA_121_20100924211938.cal'
;filename='PicA_121_20100924211938'
ext='.UVFITS'
file_path=filepath(filename,root_dir=rootdir('mwa'),subdir=data_directory)
file_path2=file_path
;temp_file=filepath('default',root=rootdir('mwa'),subdir='temporary')
temp_file=filepath('default',root=rootdir('mwa'),subdir='')

pol_names=['xx','yy','xy','yx']
time_name=(Sindgen2(100,char=2))[0] ;generalize later to allow multiple time files
;different frequency intervals can be supported later
freq_name='00'
file_identifier=Strarr(4)
FOR pol=0,3 DO file_identifier[pol]=String(format='(A2,A2,"_",A2)',pol_names[pol],time_name,freq_name)

IF Keyword_Set(restore_last) THEN BEGIN 
    beam_setup,data_directory=data_directory,filename=filename,/restore_last,file_identifier=file_identifier[0]
    ;;npol=2
    npol=1
    FOR pol=0,npol-1 DO BEGIN
        restore,filename=file_path+'_uv_'+file_identifier[pol]+'.sav'
        restore,filename=file_path+'_dirty_'+file_identifier[pol]+'.sav'
        restore,filename=file_path+'_xferfn_'+file_identifier[pol]+'.sav'            
        CASE pol OF 
            0:BEGIN
                dirty_xx_UV=temporary(dirty_uv)
                xfer_fn_xx=temporary(xfer_fn)
            END
            1:BEGIN
                dirty_yy_UV=temporary(dirty_uv)
                xfer_fn_yy=temporary(xfer_fn)
            END
            2:BEGIN
                dirty_xy_UV=temporary(dirty_uv)
                xfer_fn_xy=temporary(xfer_fn)            
            END
            3:BEGIN
                dirty_yx_UV=temporary(dirty_uv)
                xfer_fn_yx=temporary(xfer_fn)            
            END
        ENDCASE    
    ENDFOR
;    restore,filename=file_path+'_xx_uv.sav'
;    restore,filename=file_path+'_yy_uv.sav'
;;    restore,filename=file_path+'_xy_uv.sav'
;;    restore,filename=file_path+'_yx_uv.sav'
;;    restore,filename=file_path+'_data.sav'
;;    restore,filename=file_path+'_dirty_uv.sav'
;    restore,filename=file_path+'_dirty_xx.sav'
;    restore,filename=file_path+'_dirty_yy.sav'
;;    restore,filename=file_path+'_dirty_xy.sav'
;;    restore,filename=file_path+'_dirty_yx.sav'
;    visibility_transfer_generate,xfer_fn_xx,/restore_last,timing=t_xfer_gen0,/new,polarization=0
;    visibility_transfer_generate,xfer_fn_yy,/restore_last,timing=t_xfer_gen1,/new,polarization=1
;    visibility_transfer_generate,xfer_fn_xy,/restore_last,timing=t_xfer_gen2,/new,polarization=2
;    visibility_transfer_generate,xfer_fn_yx,/restore_last,timing=t_xfer_gen3,/new,polarization=3
;    weights_uv=Real_part(visibility_transfer_apply(complexarr(512,512)+1,xfer_fn))
ENDIF ELSE BEGIN
;    observation_setup,filename=filename,data_directory=data_directory
   beam_setup,data_directory=data_directory,filename=filename,restore_last=restore_beam,file_identifier=file_identifier[0]

   exten_file = file_path+'_extensions'+ext
   ftest = file_test(exten_file) *  (1 - file_test(exten_file, /zero_length))
   if ftest ne 0 then info_struct=mrdfits(exten_file,1,info_header,/silent) $
   else info_struct=mrdfits(file_path+ext,1,info_header,/silent)

   data_struct=mrdfits(file_path+ext,0,data_header,/silent)
    
    ;fits_open,file_path,fcb,no_abort=1,message=message
    n_extensions=sxpar(data_header,'nextend')
    naxis=sxpar(data_header,'naxis') ;6
    n_grp_params=sxpar(data_header,'pcount') ;5
    gcount=sxpar(data_header,'gcount') ;variable, based on length of observation
    n_complex=sxpar(data_header,'naxis2') ;3 columns are amplitude, phase (degrees), weights
    n_polarizations=sxpar(data_header,'naxis3') ;4 columns are xx, yy, xy, yx
    n_frequencies=sxpar(data_header,'naxis4') ;768
    freq_ref=sxpar(data_header,'crval4') ;1.5424E8
    freq_bandwidth=sxpar(data_header,'cdelt4') ;40000 
    freq_ref_i=sxpar(data_header,'crpix4') ;368
    frequency_array=(Dindgen(n_frequencies)-freq_ref_i)*freq_bandwidth+freq_ref
    wavelength_array=299792458./frequency_array
     
    n_fields=sxpar(info_header,'tfields') ;12
    n_tiles=sxpar(info_header,'naxis2') ;32
    
    grp_row_size=n_complex*n_polarizations*n_frequencies*gcount
    
    pol_dim=2
    freq_dim=4
    amp_index=0
    phase_index=1
    flag_index=2
    ;fits_read,fcb,0,header,/no_abort,exten_no=0,/header_only
;    commented out, since not actually used
;    ra_refval=sxpar(data_header,'crval5')
;    dec_refval=sxpar(data_header,'crval6')
;    p_uu0=sxpar(data_header,'pzero1')
;    p_uus=sxpar(data_header,'pscal1')
;    p_vv0=sxpar(data_header,'pzero2')
;    p_vvs=sxpar(data_header,'pscal2')
;    p_ww0=sxpar(data_header,'pzero3')
;    p_wws=sxpar(data_header,'pscal3')
;    p_baseline0=sxpar(data_header,'pzero4')
;    p_baselines=sxpar(data_header,'pscal4')
;    p_date0=sxpar(data_header,'pzero5')
;    p_dates=sxpar(data_header,'pscal5')
    
    uu_arr=reform(data_struct.params[0])
    vv_arr=reform(data_struct.params[1])
    ww_arr=reform(data_struct.params[2])
    baseline_arr=reform(data_struct.params[3])
    date_arr=reform(data_struct.params[4])
    data_array=data_struct.array
    data_struct=0. ;free memory
    
    flag_center=1
    flag_edge=3
    flag_arr=Reform(data_array[flag_index,*,*,*])
    coarse_channel_id=Floor(indgen(n_frequencies)/24.)
    coarse_channel_pos=indgen(n_frequencies) mod 32
    flag_arr[*,where(coarse_channel_pos EQ 16),*]=0
    flag_arr[*,where((coarse_channel_pos LT flag_edge) OR (coarse_channel_pos GT 32-flag_edge-1)),*]=0
    
;    amp_hist0=simple_histogram(data_array[amp_index,0,*,*],center=amp_mean,sigma=amp_sigma)
;    phase_hist0=simple_histogram(data_array[phase_index,0,*,*],center=phase_mean,sigma=phase_sigma)
    ;flag_amp0=where(Abs(data_array[amp_index,0,*,*]-amp_mean) GT 3.*amp_sigma,n_ampflag)
    ;flag_phase0=where(Abs(data_array[phase_index,0,*,*]-phase_mean) GT 3.*phase_sigma,n_phaseflag)
    ;IF n_ampflag GT 0 THEN flag_arr[flag_amp0]=-2
    ;IF n_phaseflag GT 0 THEN flag_arr[flag_phase0]=-2
    
    ;;tile locs are Eastings (meters), Northings (meters), elevation (meters)
    ;textfast,tile_locs,/read,filename='32T_tile_locations',filepathfull='Data',rootdir=rootdir('mwa'),column_list=[2,3,4]
    time_order_i=Sort(date_arr)
    date_arr=date_arr[time_order_i]
    times=date_arr[Uniq(date_arr)]
    times=(times-min(times))*3600.*24.
    baseline_arr=baseline_arr[time_order_i]
    data_array=data_array[*,*,*,time_order_i]
    ;256 tile upper limit is hard-coded in CASA format
    ;these tile numbers have been verified to be correct
    tile_A=Floor(baseline_arr/256) 
    tile_B=Fix(baseline_arr mod 256)
    uu_arr=uu_arr[time_order_i]
    vv_arr=vv_arr[time_order_i]
    ww_arr=ww_arr[time_order_i]
    b0i=Uniq(date_arr)
    time_step=(date_arr[b0i[1]]-date_arr[b0i[0]])*24.*3600.
    time_total=(Max(date_arr)-Min(date_arr))*24.*3600.
    nb=N_Elements(b0i)
    bin_start=fltarr(nb) & bin_start[1:*]=b0i[0:nb-2]+1
    bin_end=b0i
    time_bin=fltarr(2,nb) & time_bin[0,*]=bin_start & time_bin[1,*]=bin_end
    bin_width=fltarr(nb)
    bin_width[0]=b0i[0]+1
    FOR i=1,nb-1 DO bin_width[i]=b0i[i]-b0i[i-1]
    bin_width_c=total(bin_width,/cumulative)
    bin_offset=fltarr(nb) & bin_offset[1:*]=total(bin_width[0:nb-2],/cumulative)
    nb2=b0i[0]+1
    
    ;save configuration data for use with simulations
    file_id=Strmid(file_identifier[0],2,5)
    IF not Keyword_Set(no_save) THEN save,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,flag_arr,file=file_path+'_uv'+file_id+'.sav'
    save,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,flag_arr,file=temp_file+'_uv'+file_id+'.sav'
    
;    ;Calibration
;    ;NOT COMPATIBLE WITH MULTIPLE POLARIZATIONS YET!!
;    freq_bin=32.*freq_bandwidth ;Hz
;    freq_hist=histogram(frequency_array,locations=freq_center,binsize=freq_bin,reverse_ind=freq_ri)
;    nf=N_Elements(freq_hist)
;    baseline_matrix=0
;    tile_matrix=0
;    baseline_dist=Sqrt(uu_arr^2.+vv_arr^2.)
;    calibrated_phase=Reform(data_array[1,0,*,*])
;    IF Keyword_Set(calibrate) THEN BEGIN
;        FOR f=0,nf-1 DO BEGIN
;            IF freq_ri[f] GE freq_ri[f+1] THEN CONTINUE
;            freq_inds=freq_ri[freq_ri[f]:freq_ri[f+1]-1]
;            phase_visibilities=Reform(data_array[1,0,freq_inds,*])
;            visibility_calibrate,phase_visibilities,tile_calibration_phase,baseline_arr,time_arr=date_arr,$
;                tile_matrix=tile_matrix,baseline_matrix=baseline_matrix,flag_arr=flag_arr[freq_inds,*],$
;                baseline_threshold=baseline_threshold,baseline_dist=baseline_dist,frequency=freq_center[f]
;            calibrated_phase[freq_ri[freq_ri[f]:freq_ri[f+1]-1],*]=phase_visibilities
;    ;        FOR ff=freq_ri[freq_ri[f]],freq_ri[freq_ri[f+1]-1] DO data_array[1,0,ff,*]-=calibration_phase
;        ENDFOR
;    ENDIF

;If the visibilities have been calibrated in CASA, the filename will end with .cal and phase will be in RADIANS instead of DEGREES
    IF Strmatch(filename,'*.cal') EQ 0 THEN BEGIN
        vis_xx=reform(data_array[amp_index,0,*,*]*Exp(Complex(0,1)*data_array[phase_index,0,*,*]*!DtoR))
        vis_yy=reform(data_array[amp_index,1,*,*]*Exp(Complex(0,1)*data_array[phase_index,1,*,*]*!DtoR))
        vis_xy=reform(data_array[amp_index,2,*,*]*Exp(Complex(0,1)*data_array[phase_index,2,*,*]*!DtoR))
        vis_yx=reform(data_array[amp_index,3,*,*]*Exp(Complex(0,1)*data_array[phase_index,3,*,*]*!DtoR))
    ENDIF ELSE BEGIN
        vis_xx=reform(data_array[amp_index,0,*,*]*Exp(Complex(0,1)*data_array[phase_index,0,*,*]))
        vis_yy=reform(data_array[amp_index,1,*,*]*Exp(Complex(0,1)*data_array[phase_index,1,*,*]))
        vis_xy=reform(data_array[amp_index,2,*,*]*Exp(Complex(0,1)*data_array[phase_index,2,*,*]))
        vis_yx=reform(data_array[amp_index,3,*,*]*Exp(Complex(0,1)*data_array[phase_index,3,*,*]))
    ENDELSE
    data_array=0 ;free memory
    IF not Keyword_Set(no_save) THEN save,flag_arr,vis_xx,vis_yy,vis_xy,vis_yx,filename=file_path+'_data'+file_id+'.sav'
    save,flag_arr,vis_xx,vis_yy,vis_xy,vis_yx,filename=temp_file+'_data'+file_id+'.sav'
    
;    sim_vis0=visibility_source_setup([0,0,1],timing=t_sim_vis0)
    t_grid0=(t_grid1=(t_grid2=(t_grid3=0)))
    t_xfer_gen0=(t_xfer_gen1=(t_xfer_gen2=(t_xfer_gen3=0)))
    test_xx=where(flag_arr[0,*,*] GT 0,n_xx) & test_xx=0
    test_yy=where(flag_arr[1,*,*] GT 0,n_yy) & test_yy=0
    test_xy=where(flag_arr[2,*,*] GT 0,n_xy) & test_xy=0
    test_yx=where(flag_arr[3,*,*] GT 0,n_yx) & test_yx=0

   IF n_xx GT 0 THEN BEGIN
        dirty_UV=visibility_grid(vis_xx,reform(flag_arr[0,*,*]),weights=weights_grid,timing=t_grid0,test=test,polarization=0)
        IF not Keyword_Set(no_save) THEN save,dirty_UV,weights_grid,filename=file_path+'_uv_'+file_identifier[0]+'.sav'
        dirty_xx_UV=dirty_UV
        weights_xx_grid=weights_grid
    ENDIF ELSE dirty_xx_UV=complexarr(dimension,elements)
    IF n_yy GT 0 THEN BEGIN
        dirty_UV=visibility_grid(vis_yy,reform(flag_arr[1,*,*]),weights=weights_grid,timing=t_grid1,test=test,polarization=1)
        IF not Keyword_Set(no_save) THEN save,dirty_uv,weights_grid,filename=file_path+'_uv_'+file_identifier[1]+'.sav'
        dirty_yy_UV=dirty_UV
        weights_yy_grid=weights_grid
    ENDIF ELSE dirty_yy_UV=complexarr(dimension,elements)
    IF n_xy GT 0 THEN BEGIN
        dirty_UV=visibility_grid(vis_xy,reform(flag_arr[2,*,*]),weights=weights_grid,timing=t_grid2,test=test,polarization=2)
        IF not Keyword_Set(no_save) THEN save,dirty_uv,weights_grid,filename=file_path+'_uv_'+file_identifier[2]+'.sav'
        dirty_xy_UV=dirty_UV
        weights_xy_grid=weights_grid
    ENDIF ELSE dirty_xy_UV=complexarr(dimension,elements)
    IF n_yx GT 0 THEN BEGIN
        dirty_UV=visibility_grid(vis_yx,reform(flag_arr[3,*,*]),weights=weights_grid,timing=t_grid3,test=test,polarization=3)
        IF not Keyword_Set(no_save) THEN save,dirty_uv,weights_grid,filename=file_path+'_uv_'+file_identifier[3]+'.sav'
        dirty_yx_UV=dirty_UV
        weights_yx_grid=weights_grid        
    ENDIF ELSE dirty_yx_UV=complexarr(dimension,elements)
    
    IF n_xx GT 0 THEN BEGIN
        visibility_transfer_generate,xfer_fn,restore_last=0,timing=t_xfer_gen0,/new,polarization=0,flag_arr=reform(flag_arr[0,*,*]),file_identifier=file_identifier[0]
        xfer_fn=0
    ENDIF
    IF n_yy GT 0 THEN BEGIN
        visibility_transfer_generate,xfer_fn,restore_last=0,timing=t_xfer_gen1,/new,polarization=1,flag_arr=reform(flag_arr[1,*,*]),file_identifier=file_identifier[1]
        xfer_fn=0
    ENDIF
    IF n_xy GT 0 THEN BEGIN
        visibility_transfer_generate,xfer_fn,restore_last=0,timing=t_xfer_gen2,/new,polarization=2,flag_arr=reform(flag_arr[2,*,*]),file_identifier=file_identifier[2]
        xfer_fn=0
    ENDIF
    IF n_yx GT 0 THEN BEGIN
        visibility_transfer_generate,xfer_fn,restore_last=0,timing=t_xfer_gen3,/new,polarization=3,flag_arr=reform(flag_arr[3,*,*]),file_identifier=file_identifier[3]
        xfer_fn=0
    ENDIF
    dirty_xx_img=dirty_image_generate(dirty_xx_UV,baseline_threshold=baseline_threshold)
    dirty_yy_img=dirty_image_generate(dirty_yy_uv,baseline_threshold=baseline_threshold)
    dirty_xy_img=dirty_image_generate(dirty_xy_uv,baseline_threshold=baseline_threshold)
    dirty_yx_img=dirty_image_generate(dirty_yx_uv,baseline_threshold=baseline_threshold)
    IF not Keyword_Set(no_save) THEN BEGIN
        save,dirty_xx_img,filename=file_path+'_dirty_'+file_identifier[0]+'.sav'
        save,dirty_yy_img,filename=file_path+'_dirty_'+file_identifier[1]+'.sav'
        save,dirty_xy_img,filename=file_path+'_dirty_'+file_identifier[2]+'.sav'
        save,dirty_yx_img,filename=file_path+'_dirty_'+file_identifier[3]+'.sav'
    ENDIF
    print,'Gridding time:',t_grid0,t_grid1,t_grid2,t_grid3
    print,'Time spent building the transfer function:',t_xfer_gen0,t_xfer_gen1,t_xfer_gen2,t_xfer_gen3
ENDELSE

beam_base_uv=Ptrarr(4,/allocate)
beam_base=Ptrarr(4,/allocate)
beam_corr=Ptrarr(4,/allocate)
beam_i=Ptrarr(4,/allocate)
beam_i2=Ptrarr(4,/allocate)
beam_i3=Ptrarr(4,/allocate)
weights=Ptrarr(4,/allocate)
dirty_image=Ptrarr(4,/allocate)
clean0_image=Ptrarr(4,/allocate)
clean0b_image=Ptrarr(4,/allocate)
clean_image_arr=Ptrarr(4,/allocate)
dirty_uv=Ptrarr(4,/allocate)
xfer_fn_arr=Ptrarr(4,/allocate)
source_array2=Ptrarr(4,/allocate)
model_uv_arr=Ptrarr(4,/allocate)

test_beam=Ptrarr(4,/allocate)
test_corr=Ptrarr(4,/allocate)
;*weights[0]=weights_xx_grid & *weights[1]=weights_yy_grid ;& *weights[2]=weights_xy_grid & *weights[3]=weights_yx_grid 
*dirty_image[0]=dirty_xx_img & *dirty_image[1]=dirty_yy_img ;& *dirty_image[2]=dirty_xy_img & *dirty_image[3]=dirty_yx_img
*dirty_uv[0]=dirty_xx_uv & *dirty_uv[1]=dirty_yy_uv ;& *dirty_uv[2]=dirty_xy_uv & *dirty_uv[3]=dirty_yx_uv
*xfer_fn_arr[0]=xfer_fn_xx & *xfer_fn_arr[1]=xfer_fn_yy 
scale=fltarr(4)
beam_threshold=2./3.
FOR pol=0,1 DO BEGIN
    *beam_base_uv[pol]=fltarr(dimension,elements)
    (*beam_base_uv[pol])[dimension/2.-Floor(psf_dim/2.):dimension/2.-Floor(psf_dim/2.)+psf_dim-1,elements/2.-Floor(psf_dim/2.):elements/2.-Floor(psf_dim/2.)+psf_dim-1]$
        =*psf_base[pol,max(freq_bin_i)/2,0,0]
    test_fit=gauss2dfit(*beam_base_uv[pol])
    test_fit-=Min(test_fit)
    test_fit[dimension/2.-Floor(psf_dim/2.):dimension/2.-Floor(psf_dim/2.)+psf_dim-1,elements/2.-Floor(psf_dim/2.):elements/2.-Floor(psf_dim/2.)+psf_dim-1]$
        =*psf_base[pol,max(freq_bin_i)/2,0,0]
    *test_beam[pol]=fft_shift(real_part(FFT(fft_shift(test_fit),/inverse))) 
    *test_beam[pol]/=Max(*test_beam[pol])   
    *beam_base[pol]=fft_shift(real_part(FFT(fft_shift(*beam_base_uv[pol]),/inverse)))
    *beam_base[pol]/=max(*beam_base[pol])
    *beam_i[pol]=region_grow(*beam_base[pol],dimension/2.+dimension*elements/2.,threshold=[beam_threshold,1.])
    *beam_i2[pol]=region_grow(*beam_base[pol],dimension/2.+dimension*elements/2.,threshold=[beam_threshold/2.,1.])
    *beam_corr[pol]=fltarr(dimension,elements)
    (*beam_corr[pol])[*beam_i2[pol]]=1./(*beam_base[pol])[*beam_i2[pol]]
    *beam_i3[pol]=region_grow(*test_beam[pol],dimension/2.+dimension*elements/2.,threshold=[beam_threshold/2.,1.])
    *test_corr[pol]=fltarr(dimension,elements)
    (*test_corr[pol])[*beam_i3[pol]]=1./(*test_beam[pol])[*beam_i3[pol]]
    
    *weights[pol]=real_part(visibility_transfer_apply(complexarr(dimension,elements)+1,*xfer_fn_arr[pi],/new))
    scale[pol]=(linfit((dirty_image_generate(*weights[pol],base=baseline_threshold))[*beam_i[pol]],(*dirty_image[pol])[*beam_i[pol]]))[1]
    *clean0_image[pol]=(*dirty_image[pol]-dirty_image_generate(*weights[pol]*scale[pol],base=baseline_threshold))*(*beam_corr[pol])
    *clean0b_image[pol]=(*dirty_image[pol]-dirty_image_generate(*weights[pol]*scale[pol],base=baseline_threshold))*(*test_corr[pol])
    
    max_iter=5000.
    visibility_image_clean,*dirty_uv[pol],clean_image_single,*xfer_fn_arr[pol],max_iter=max_iter,$
        source_array=source_array_single,model_uv=model_uv_single,gain_factor=0.07,timing=t_clean
    *clean_image_arr[pol]=clean_image_single
    *source_array2[pol]=clean_source_condense(source_array_single,radius=.25)
    *model_uv_arr[pol]=visibility_source_uv_grid(*source_array2[pol])
ENDFOR
IF not Keyword_Set(no_save) THEN save,beam_base_uv,beam_base,beam_corr,weights,dirty_image,clean_image_arr,dirty_uv,xfer_fn_arr,source_array2,model_uv_arr,filename=temp_file+'_clean.sav'

END
