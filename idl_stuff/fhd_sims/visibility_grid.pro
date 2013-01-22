FUNCTION visibility_grid,visibility_array,flag_arr,weights=weights,timing=timing,test=test,polarization=polarization, $
                         freq_inds_to_use = freq_inds_to_use
;This function grids any visibilities to a UV plane.
;visibility_analysis_setup must be run first to set up the psf

t0=Systime(1)
COMMON obs_params,dimension,elements,degpix,kbinsize,azimuth_start,elevation_start,rotation,JD0
;COMMON psf_params,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals
COMMON visibility_params,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,freq_bin_i,file_path

dirty_image_UV=Dcomplexarr(dimension,elements)
weights=fltarr(dimension,elements)
weights_comp=fltarr(dimension,elements)

kx_arr=uu_arr/kbinsize
ky_arr=vv_arr/kbinsize
;;nbaselines=bin_offset[1] -- changed to support n_samples=1
n_samples=N_Elements(bin_offset)
nbaselines=n_elements(baseline_arr)/n_samples

if n_elements(freq_inds_to_use) eq 0 then freq_inds_to_use = indgen(n_elements(frequency_array))
n_frequencies=N_Elements(frequency_array[freq_inds_to_use])
n_freq_bin=N_Elements(freq_bin_i[freq_inds_to_use])
t1=0
t2=0
t3=0
time_check_interval=Ceil(nbaselines/20.)
T0a=Systime(1)

;;IF N_Elements(flag_arr) EQ 0 THEN flag_arr=fltarr(n_frequencies,nbaselines*bin_offset[1])+1.
IF N_Elements(flag_arr) EQ 0 THEN flag_arr=fltarr(n_frequencies,nbaselines*n_samples)+1.
psf_i=1
FOR bi=0.,nbaselines-1 DO BEGIN
;    bi=vis_i mod nbaselines
;    fi=Floor(vis_i/nbaselines)
;    baseline_i=baseline_arr[bi]
    baseline_i=bi mod nbaselines
    
    IF Keyword_Set(time_check_interval) AND (bi mod time_check_interval EQ 0) THEN BEGIN
        t0b=Systime(1)-t0a
        print,Strcompress(String(format='("Time elapsed:",I," Estimated time remaining:",I)',t0b,(t0b/bi)*(nbaselines-bi)))
    ENDIF 
    IF Keyword_set(test) THEN test_img=dblarr(dimension,elements)
    IF kx_arr[bi] EQ 0 THEN CONTINUE ;skip the auto-correlations
    FOR fi=0,n_frequencies-1 DO BEGIN
        freq_i=freq_bin_i[freq_inds_to_use[fi]]
        freq=frequency_array[freq_inds_to_use[fi]]
        FOR ti=0.,n_samples-1 DO BEGIN
            bi_use=bi+bin_offset[ti]
                        
            IF flag_arr[fi,bi_use] LE 0 THEN CONTINUE
            vis_value=visibility_array[fi,bi_use]
            IF vis_value EQ 0 THEN CONTINUE
            t1_0=Systime(1)
            psf_val=visibility_psf(freq_i,baseline_i,polarization=polarization,xcenter=kx_arr[bi_use]*freq,ycenter=ky_arr[bi_use]*freq,psf_i=psf_i)
;            psf_val=visibility_psf(freq_i,baseline_i,polarization=polarization,xcenter=kx_arr[bi_use]*freq-.5,ycenter=ky_arr[bi_use]*freq-.5,psf_i=psf_i)
            t2_0=Systime(1)
            t1+=t2_0-t1_0
            psf_x=psf_i mod dimension
            psf_y=Floor(psf_i/dimension)
            psf_x2=dimension-0.-psf_x
            psf_y2=elements-0.-psf_y
            psf_i2=psf_x2+psf_y2*dimension
    ;        dirty_image_UV[psf_x2,psf_y2]+=psf_val*Conj(Mean(visibility_array[fi,baseline_i+bin_offset]))
    ;        dirty_image_UV[psf_i]+=psf_val*Mean(visibility_array[fi,baseline_i+bin_offset])
    
            vis_psf_vals=psf_val*vis_value
            t3_0=Systime(1)
            t2+=t3_0-t2_0
;            IF Median(psf_y) GT elements/2. THEN BEGIN
;                dirty_image_UV[psf_i]+=Conj(vis_psf_vals)
;                dirty_image_UV[psf_i2]+=vis_psf_vals
;                IF Keyword_set(test) THEN BEGIN
;                    weights_comp[psf_i]-=psf_val
;                    weights_comp[psf_i2]+=psf_val
;                ENDIF            
;            ENDIF ELSE BEGIN
                dirty_image_UV[psf_i2]+=Conj(vis_psf_vals)
                dirty_image_UV[psf_i]+=vis_psf_vals
                IF Keyword_set(test) THEN BEGIN
                    weights_comp[psf_i2]-=psf_val
                    weights_comp[psf_i]+=psf_val
                ENDIF
;            ENDELSE
            weights[psf_i]+=psf_val
            weights[psf_i2]+=psf_val
            t3+=Systime(1)-t3_0
            break_loc=1
        ENDFOR
    ENDFOR
ENDFOR
IF Keyword_set(test) THEN test=weights_comp
timing=Systime(1)-t0
RETURN,dirty_image_UV
END
