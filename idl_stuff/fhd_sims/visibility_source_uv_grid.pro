FUNCTION visibility_source_uv_grid,source_array, max_sources=max_sources, timing=timing, u_dim = u_dim, v_dim = v_dim
t0=Systime(1)
    
if n_elements(u_dim) eq 0 or n_elements(v_dim) eq 0 then begin
   COMMON obs_params,dimension,elements,degpix,kbinsize,azimuth_center,elevation_center,rotation
endif else begin
   dimension = u_dim
   elements = v_dim
endelse


icomp=complex(0,1)
fft_source_image=dcomplexarr(dimension,elements)
n_sources=N_Elements(source_array[0,*])

IF Keyword_Set(max_sources) THEN n_sources=max_sources<n_sources
;xvals=Reform(source_array[0,0:n_sources-1])
;yvals=Reform(source_array[1,0:n_sources-1])
;flux=Reform(source_array[2,0:n_sources-1])
;
;FOR ii=0,dimension-1 DO BEGIN
;    FOR jj=0,elements-1 DO BEGIN
;        fft_source_image[ii,jj]=Total(flux*Exp(icomp*(-2d*!DPi/dimension)*((ii-dimension/2.)*xvals+(jj-elements/2.)*yvals)))
;    ENDFOR
;ENDFOR

sdim=(size(source_array,/dimension))[0]
IF sdim EQ 4 THEN BEGIN
    ;source array format compatible with visibility_image_clean.pro
    sx=Reform(source_array[0,*])
    sy=Reform(source_array[1,*])
    sf=Reform(source_array[3,*])
ENDIF ELSE BEGIN
    ;source array format compatible with general programs
    sx=Reform(source_array[0,*])-dimension/2.
    sy=Reform(source_array[1,*])-elements/2.
    sf=Reform(source_array[4,*])
ENDELSE

xvals=meshgrid(dimension,elements,1)-dimension/2
yvals=meshgrid(dimension,elements,2)-elements/2
FOR si=0L,n_sources-1 DO fft_source_image+=sf[si]*Exp(icomp*(-2d*!DPi/dimension)*(sx[si]*xvals+sy[si]*yvals))

timing=Systime(1)-t0
RETURN,fft_source_image

END
