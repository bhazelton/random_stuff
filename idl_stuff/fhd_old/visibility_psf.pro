FUNCTION visibility_psf,freq_i,baseline_i,xcenter=xcenter,ycenter=ycenter,$
    image_sample=image_sample,psf_i=psf_i,short_circuit_flag=short_circuit_flag,$
    xmin=xmin,ymin=ymin,xcen0=xcen0,ycen0=ycen0,polarization=polarization;,conjugate=conjugate
COMMON obs_params,dimension,elements,degpix,kbinsize
COMMON psf_params,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals
COMMON switches,source_switch,psf_switch
IF N_Elements(xcenter) EQ 0 THEN BEGIN short_circuit_flag=1 & xcenter=0 & ycenter=0 & ENDIF

;If keyword image_sample is set, then this function will instead return Total(psf*image_sample) (this saves computations)
;tile_A=Floor(baseline_i/256) 
;tile_B=Fix(baseline_i mod 256)

CASE psf_switch OF
    1: BEGIN
        ;high resolution nearest neighbor method
;        psf_use=psf_base[polarization,freq_i,*,*]
;;        <<<PSF RESIDUALS ARE NOT YET SUPPORTED IN THIS FORMAT!!!>>>
;        IF psf_residuals_n[freq_i,baseline_i] GT 0 THEN $
;            psf_use[*psf_residuals_i[freq_i,baseline_i]]+=*psf_residuals_val[freq_i,baseline_i]
        x_offset=Round((Ceil(xcenter)-xcenter)*psf_resolution) mod psf_resolution    
        y_offset=Round((Ceil(ycenter)-ycenter)*psf_resolution) mod psf_resolution
        xcen0=Round(xcenter+x_offset/psf_resolution+dimension/2.) ;do this after offset, in case it has rounded to the next grid point
        ycen0=Round(ycenter+y_offset/psf_resolution+elements/2.)
        psf_use=*psf_base[polarization,freq_i,x_offset,y_offset]
    END
    2: BEGIN
        ;Interpolation method        
        psf_use=psf_base[polarization,freq_i,*,*]
        IF psf_residuals_n[polarization,freq_i,baseline_i] GT 0 THEN $
            psf_use[*psf_residuals_i[polarization,freq_i,baseline_i]]+=*psf_residuals_val[polarization,freq_i,baseline_i]
        x_offset=Ceil(xcenter)-xcenter   
        y_offset=Ceil(ycenter)-ycenter
        xcen0=Round(xcenter+x_offset+dimension/2.) 
        ycen0=Round(ycenter+y_offset+elements/2.)
        psf_use=Interpolate(psf_use,(findgen(psf_dim)+x_offset)*psf_resolution,(findgen(psf_dim)+y_offset)*psf_resolution,/grid)
    END    
ENDCASE
xmin=Floor(xcen0-psf_dim/2.) & xmax=xmin+psf_dim-1
ymin=Floor(ycen0-psf_dim/2.) & ymax=ymin+psf_dim-1

IF Keyword_Set(short_circuit_flag) THEN RETURN, psf_use

CASE 1 OF
    Keyword_Set(image_sample) : result=Total(psf_use*image_sample[xmin:xmax,ymin:ymax],/double)
    Keyword_Set(psf_i) : BEGIN
        psf_use_i=where(psf_use)
        result=psf_use[psf_use_i]
        psf_i=(meshgrid(psf_dim,psf_dim,1)+xmin)+(meshgrid(psf_dim,psf_dim,2)+ymin)*dimension
        psf_i=psf_i[psf_use_i]
    END
    ELSE : BEGIN
        result=fltarr(dimension,elements)
        result[xmin:xmax,ymin:ymax]=psf_use
    ENDELSE
ENDCASE

RETURN, result
END
