FUNCTION clean_source_restore,clean_residual,source_array,width=width,source_image=source_image,ring_radius=ring_radius

IF N_Elements(width) EQ 0 THEN width=1.
IF N_Elements(ring_radius) EQ 0 THEN ring_radius=0

dimension=(size(clean_residual,/dimension))[0]
elements=(size(clean_residual,/dimension))[1]
x_vals=meshgrid(dimension,elements,1)-dimension/2.
y_vals=meshgrid(dimension,elements,2)-elements/2.

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
ns=N_Elements(sx)

source_image=fltarr(dimension,elements)
IF Keyword_Set(ring_radius) THEN FOR si=0,ns-1 DO $
    source_image+=sf[si]*Exp(-Abs(((x_vals-sx[si])/width)^2.+((y_vals-sy[si])/width)^2.-2.*ring_radius)/2.) $
ELSE FOR si=0.,ns-1 DO $
    source_image+=sf[si]*Exp(-(((x_vals-sx[si])/width)^2.+((y_vals-sy[si])/width)^2.)/2.) 
recon_image=clean_residual+source_image
RETURN,recon_image
END

