PRO visibility_image_clean,image_uv,clean_image,xfer_fn,source_array=source_array,beam_threshold=beam_threshold,$
    model_uv=model_uv,max_iter=max_iter,gain_factor=gain_factor,center_gain=center_gain,$
    transfer_interval=transfer_interval,transfer_threshold=transfer_threshold,timing=timing,$
    beam_base=beam_base,polarization=polarization,baseline_threshold=baseline_threshold
    
t0=Systime(1)    
heap_gc
COMMON obs_params,dimension,elements,degpix,kbinsize,azimuth_start,elevation_start,rotation,JD0
COMMON psf_params,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals
COMMON visibility_params,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,freq_bin_i,file_path2
IF N_Elements(baseline_threshold) EQ 0 THEN baseline_threshold=50.
IF N_Elements(gain_factor) EQ 0 THEN gain_factor=0.2
IF N_Elements(center_gain) EQ 0 THEN center_gain=0.95 ;set to 0 to NOT remove a source at the exact center
IF N_Elements(max_iter) EQ 0 THEN max_iter=100./gain_factor
IF N_Elements(check_iter) EQ 0 THEN check_iter=Round(10./gain_factor)
;transfer_interval is max # of iter. before running full transfer fn
IF N_Elements(transfer_interval) EQ 0 THEN transfer_interval=20. 
;if source flux has dropped by transfer_threshold, then run full transfer fn even if not yet at transfer_interval 
IF N_Elements(transfer_threshold) EQ 0 THEN transfer_threshold=0.85
IF N_Elements(polarization) EQ 0 THEN polarization=0.
radius=3.
;beam_threshold=2./3. ;fraction of beam center to include in finding sources
IF N_Elements(beam_threshold) EQ 0 THEN beam_threshold=0.1 ;fraction of beam center to include in finding sources

weights=abs(visibility_transfer_apply(complexarr(dimension,elements)+complex(1,0),xfer_fn,/new))
weights_inv=weight_invert(weights,Max(weights)/100.)
mask_uv=fltarr(dimension,elements)
mask_uv[where(weights_inv)]=1.

IF N_Elements(beam_base) EQ 0 THEN BEGIN
    beam_base_uv=fltarr(dimension,elements)
    beam_base_uv[dimension/2.-Floor(psf_dim/2.):dimension/2.-Floor(psf_dim/2.)+psf_dim-1,elements/2.-Floor(psf_dim/2.):elements/2.-Floor(psf_dim/2.)+psf_dim-1]$
        =*psf_base[polarization,max(freq_bin_i)/2,0,0]
    beam_base=fft_shift(real_part(FFT(fft_shift(beam_base_uv),/inverse)))
ENDIF
beam_base/=max(beam_base)
beam_i=region_grow(beam_base,dimension/2.+dimension*elements/2.,threshold=[beam_threshold,1.])
beam_correction=fltarr(dimension,elements)
beam_correction[beam_i]=1./beam_base[beam_i]

;gain_check=Dirty_image_generate(weights,baseline_threshold=baseline_threshold)
;gain_check2=Smooth(gain_check^2.,2.*radius,/edge)
;gain_factor2=Sqrt(Max(gain_check2))*radius^2.  
;gain_factor_use=gain_factor/gain_factor2;*Total(abs(gain_check))

xvals=meshgrid(dimension,elements,1)-dimension/2
yvals=meshgrid(dimension,elements,2)-elements/2
icomp=Complex(0,1)

normalization=1./Mean(weights[where(weights)])
gain_check=Dirty_image_generate(weights,baseline_threshold=baseline_threshold,normalization=normalization)
gain_factor2=1./gain_check[dimension/2.,elements/2.]
gain_factor_use=gain_factor*gain_factor2

image_uv_use=image_uv*weights_inv
dirty_image=dirty_image_generate(image_uv,baseline_threshold=baseline_threshold,mask=mask_uv2,normalization=normalization)
center_source_image=Dirty_image_generate(weights,baseline_threshold=baseline_threshold,normalization=normalization)
scale_center=(linfit(center_source_image[beam_i],dirty_image[beam_i]))[1]
source_array1=[0,0,scale_center*center_gain,scale_center/gain_factor2]
model_uv_full=source_array1[2]*Exp(icomp*(-2d*!DPi/dimension)*(source_array1[0]*xvals+source_array1[1]*yvals))

;weights is by definition the transfer function applied to a source at the center, but run transfer function for general source location case
model_uv_xfer=visibility_transfer_apply(model_uv_full,xfer_fn,/new)
model_uv=dcomplexarr(dimension,elements)
source_array=Fltarr(4,max_iter)
source_array[*,0]=source_array1
converge_check=fltarr(Ceil(max_iter*gain_factor))
converge_check2=fltarr(max_iter)
converge_check[0]=Stddev(dirty_image,/nan)
converge_check2[0]=Stddev(dirty_image)

t1=0 & t2=0 & t3=0 & t4=0
i2=0 & i3=0
FOR i=1L,max_iter-1 DO BEGIN
    t1_0=Systime(1)
    model_image_xfer=dirty_image_generate(model_uv_xfer,baseline_threshold=baseline_threshold,normalization=normalization)
    model_image_est=dirty_image_generate(model_uv*weights,baseline_threshold=baseline_threshold,normalization=normalization)*beam_base^2.
    model_image=model_image_xfer+model_image_est
    
    image_use=dirty_image-model_image
    t2_0=Systime(1)
    t1+=t2_0-t1_0    
    source_array1=single_source_extract(image_use,beam_correction,radius=radius,gain_factor=gain_factor_use,error=error);,/median_remove)
    IF Keyword_Set(error) THEN BEGIN
        print,String(format='("Failure to centroid after",I," iterations, at (",I,",",I,")")',i,source_array1[0],source_array1[1])
        BREAK
    ENDIF
    converge_check2[i]=Stddev(image_use,/nan)
    source_array[*,i]=source_array1
    t3_0=Systime(1)
    t2+=t3_0-t2_0
    source_uv=source_array1[2]*Exp(icomp*(-2d*!DPi/dimension)*(source_array1[0]*xvals+source_array1[1]*yvals));/(dimension*elements)
    model_uv+=source_uv
    model_uv_full+=source_uv
    t4_0=Systime(1)
    t3+=t4_0-t3_0
    
    IF i3 EQ 0 THEN flux_ref=source_array1[2]
    IF (i3 GE transfer_interval) OR (source_array1[2] LT flux_ref*transfer_threshold) OR ((i+1) mod check_iter EQ 0) THEN BEGIN
        model_uv_xfer=visibility_transfer_apply(model_uv_full,xfer_fn,/new)
        model_uv=Dcomplexarr(dimension,elements)
        i3=0
    ENDIF ELSE i3+=1
    t4+=Systime(1)-t4_0
    
    IF (i mod check_iter EQ 0) AND (i GT 0) THEN BEGIN
        i2+=1
        t10=Systime(1)-t0
        print,StrCompress(String(format='("At iteration ",I," after ",I," seconds (convergence:",F,")")',i,t10,Stddev(image_use,/nan)))
        converge_check[i2]=Stddev(image_use,/nan)
        IF converge_check[i2] GT converge_check[i2-1] THEN BEGIN
            print,'Break after iteration',i,' from lack of convergence'
            converge_check2=converge_check2[0:i]
            BREAK
        ENDIF
    ENDIF
ENDFOR

source_array=source_array[*,0:(i<(max_iter-1))]

model_uv=model_uv_full
model_uv_xfer=visibility_transfer_apply(model_uv,xfer_fn,/new)
dirty_source_image=dirty_image_generate(model_uv_xfer,baseline_threshold=baseline_threshold,normalization=normalization)
clean_image=dirty_image_generate((image_uv-model_uv_xfer),baseline_threshold=baseline_threshold,normalization=normalization)

timing=Systime(1)-t0
END

