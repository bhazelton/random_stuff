FUNCTION visibility_source_setup,source_array,fft_source_image=fft_source_image,timing=timing
;This function is only used when generating simulations of visibilities
t0=Systime(1)
COMMON obs_params,dimension,elements,degpix,kbinsize
COMMON psf_params,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n;,psf_xvals,psf_yvals
COMMON visibility_params,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,freq_bin_i,file_path
COMMON switches,source_switch,psf_switch

kx_arr=uu_arr/kbinsize
ky_arr=vv_arr/kbinsize
nbaselines=bin_offset[1]
n_samples=N_Elements(bin_offset)
n_frequencies=N_Elements(frequency_array)
n_freq_bin=N_Elements(freq_bin_i)
n_uniq_vis=nbaselines*n_freq_bin
visibility_array=Dcomplexarr(n_frequencies,nbaselines*n_samples)

i=complex(0,1)

;easy to write way (slightly less accurate, probably very slow)
fft_source_image=dcomplexarr(dimension,elements)
xvals=Reform(source_array[0,*])
yvals=Reform(source_array[1,*])
flux=Reform(source_array[2,*])
n_sources=N_Elements(flux)
IF Arg_present(fft_source_image) OR (source_switch EQ 1) THEN $
    fft_source_image=visibility_source_uv_grid(source_array)

psf_xvals0=meshgrid(psf_dim,psf_dim,1)-psf_dim/2.
psf_yvals0=meshgrid(psf_dim,psf_dim,2)-psf_dim/2.

FOR bi=0.,nbaselines-1 DO BEGIN
;    bi=vis_i mod nbaselines
;    fi=Floor(vis_i/nbaselines)
;    baseline_i=baseline_arr[bi]
    baseline_i=bi mod nbaselines
    IF kx_arr[bi] EQ 0 THEN CONTINUE ;skip the auto-correlations
    FOR fi=0,n_frequencies-1 DO BEGIN
        freq_i=freq_bin_i[fi]
        freq=frequency_array[fi]
        
        CASE source_switch OF
            0:FOR ti=0,n_samples-1 DO BEGIN
                bi_use=bi+bin_offset[ti]
                visibility_array[fi,bi_use]=$
                    Total(flux*Exp(i*(-2d*!DPi/dimension)*(kx_arr[bi_use]*freq*xvals+ky_arr[bi_use]*freq*yvals)))
            ENDFOR
            1:FOR ti=0,n_samples-1 DO BEGIN
                bi_use=bi+bin_offset[ti]
                visibility_array[fi,bi_use]=$
                    visibility_psf(freq_i,baseline_i,xcenter=kx_arr[bi_use]*freq,ycenter=ky_arr[bi_use]*freq,image_sample=fft_source_image)
            ENDFOR
            2:BEGIN
                psf_use=visibility_psf(freq_i,baseline_i)
                FOR ti=0,n_samples-1 DO BEGIN
                    bi_use=bi+bin_offset[ti]
                    FOR si=0,n_sources-1 DO visibility_array[fi,bi_use]+=$
                        Total(psf_use*flux[si]*Exp(i*(-2d*!DPi/dimension)*((kx_arr[bi_use]*freq+psf_xvals0)*xvals[si]+(ky_arr[bi_use]*freq+psf_yvals0)*yvals[si])))
                ENDFOR
            END
        ENDCASE
        
    ENDFOR
ENDFOR

timing=Systime(1)-t0
RETURN,visibility_array
END