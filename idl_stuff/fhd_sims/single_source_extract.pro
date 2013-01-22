FUNCTION single_source_extract,image,beam_correction,radius=radius,gain_factor=gain_factor,$
    image_correct=image_correct,absolute=absolute,error=error,median_remove=median_remove,source_i=source_i

COMMON obs_params,dimension,elements,degpix,kbinsize
IF N_Elements(radius) EQ 0 THEN radius=3.
IF N_Elements(gain_factor) EQ 0 THEN gain_factor=1.

image_use=image*beam_correction;^2.
IF Keyword_Set(median_remove) THEN image_use-=Median(image_use,radius*5.)
IF Keyword_Set(absolute) THEN BEGIN
    dummy_flux=Max(Abs(image_use),source_i)
    source_flux=image_use[source_i]
ENDIF ELSE source_flux=Max(image_use,source_i)
sx=(source_i mod dimension)
sy=Floor(source_i/dimension)
gcntrd,image_use,sx,sy,xcen,ycen,radius,/silent 

IF xcen EQ -1 THEN BEGIN
    error=1 
    xcen=sx
    ycen=sy
ENDIF ELSE error=0

xcen2=xcen-dimension/2.
ycen2=ycen-elements/2.
source_flux_corr2=source_flux*beam_correction[source_i]
source_flux_corr1=source_flux_corr2*gain_factor

single_source_array=[xcen2,ycen2,source_flux_corr1,source_flux_corr2]

RETURN,single_source_array
END