FUNCTION dirty_image_generate,dirty_image_uv,baseline_threshold=baseline_threshold,mask=mask,$
    normalization=normalization,pad=pad,resize=resize
IF N_Elements(baseline_threshold) EQ 0 THEN baseline_threshold=0.

IF N_Elements(normalization) EQ 0 THEN normalization=1
dimension=(size(dirty_image_uv,/dimension))[0]
elements=(size(dirty_image_uv,/dimension))[1]
rarray=Sqrt((meshgrid(dimension,1)-dimension/2)^2.+(meshgrid(elements,2)-elements/2.)^2.)
cut_i=where(rarray LT baseline_threshold,n_cut)
mask=fltarr(dimension,elements)+1.
IF n_cut GT 0 THEN mask[cut_i]=0
cut_i=where(dirty_image_uv*mask EQ 0,n_cut,comp=keep_i,ncomp=n_keep)
di_uv_use=dirty_image_uv*mask
;IF n_keep GT 0 THEN di_uv_use[keep_i]-=Mean(di_uv_use[keep_i])
;noise_pad=RandomN(seed,n_cut)*stddev(di_uv_use[keep_i])/100.
;di_uv_use[cut_i]=noise_pad
;IF Keyword_Set(pad) THEN BEGIN
;    di_uv_temp=dcomplexarr(dimension+2*pad,elements+2*pad)
;    di_uv_temp[pad:dimension+pad-1,pad:elements+pad-1]=di_uv_use
;    di_uv_use=di_uv_temp
;ENDIF
IF Keyword_Set(resize) THEN BEGIN
    dimension2=dimension*resize
    elements2=elements*resize
    di_uv_real=Real_part(di_uv_use)
    di_uv_img=Imaginary(di_uv_use)
    di_uv_real=REBIN(di_uv_real,dimension2,elements2)
    di_uv_img=REBIN(di_uv_img,dimension2,elements2)
    di_uv_use=complex(di_uv_real,di_uv_img)    
ENDIF

dirty_image=Real_part(fft_shift(FFT(fft_shift(di_uv_use),/inverse)))*normalization;*(2*!Pi)/(dimension*elements)

RETURN,dirty_image
END