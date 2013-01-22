FUNCTION visibility_transfer_apply,image,xfer_fn,mask,weights=weights,timing=timing,new=new,transpose=transpose
t0=Systime(1)
IF N_Elements(new) EQ 0 THEN new=1

COMMON obs_params,dimension,elements,degpix,kbinsize,azimuth_center,elevation_center,rotation
COMMON psf_params,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals

if n_elements(dimension) eq 0 or n_elements(elements) eq 0 then begin
   dims = size(image,/dimension)
   dimension = dims[0]
   elements = dims[1]
endif


IF Keyword_Set(new) THEN BEGIN
    image_real_vector=reform(Real_part(image),Float(dimension)*elements)
    image_comp_vector=reform(Imaginary(image),Float(dimension)*elements)
    
;    result_image_real = SPRSAX2(xfer_fn,image_real_vector,transpose=transpose,mask=mask)
;    result_image_comp = SPRSAX2(xfer_fn,image_comp_vector,transpose=transpose,mask=mask)
    SPRSAX2,xfer_fn,image_real_vector,result_image_real,image_comp_vector,result_image_comp,transpose=0,mask=0
    result_image=Complex(result_image_real,result_image_comp)
    result_image=reform(result_image,dimension,elements)
    result_image_conj=Shift(Reverse(reverse(Conj(result_image),1),2),1,1)
    result_image+=result_image_conj
ENDIF ELSE BEGIN
    conj_image=Conj(image)
    result_image=complexarr(dimension,elements)
    weights=fltarr(dimension,elements)
    psf_dim2=2*psf_dim
    sub_xv=meshgrid(psf_dim2,1)-psf_dim-1
    sub_yv=meshgrid(psf_dim2,2)-psf_dim-1
    
    FOR xi=0,dimension-1 DO BEGIN
        FOR yi=0,elements-1 DO BEGIN
            xfer_fn_sub=*xfer_fn[xi,yi];,*,*]
            wt=Total(xfer_fn_sub)
            IF wt EQ 0 THEN CONTINUE
            xi2=dimension-xi
            yi2=elements-yi
            weights[xi,yi]+=wt
            weights[xi2,yi2]+=wt
            xii_arr=sub_xv+xi
            yii_arr=sub_yv+yi
            result_image[xi,yi]+=Total(xfer_fn_sub*image[xii_arr,yii_arr])
            result_image[xi2,yi2]+=Total(xfer_fn_sub*conj_image[xii_arr,yii_arr])
        ENDFOR
    ENDFOR
ENDELSE

;;normalization=2.*!Pi*!Radeg^2.
;normalization=((dimension*elements)/(2.*!Pi));^2.
;result_image/=normalization

timing=Systime(1)-t0
RETURN,result_image
END
