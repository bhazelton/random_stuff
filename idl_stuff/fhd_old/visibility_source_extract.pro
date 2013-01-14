FUNCTION visibility_source_extract,image,beam_use=beam_use,max_sources=max_sources,timing=timing,beam_threshold=beam_threshold
t0=Systime(1)
IF N_Elements(beam_threshold) EQ 0 THEN beam_threshold=0.1
LocalMaxRadius=3
fit_radius=Floor(LocalMaxRadius*3./2.)
dist_tol=1.
ellipticity_threshold=2. ;Max or 1/Min ellipticity of local maxima to consider sources

dimension=(size(image,/dimension))[0]
elements=(size(image,/dimension))[1]
IF N_Elements(beam_use) EQ 0 THEN beam_use=fltarr(dimension,elements)+1.
beam_i=region_grow(beam_use,dimension/2.+dimension*elements/2.,threshold=[beam_threshold,1.])
beam_corr=fltarr(dimension,elements) & beam_corr[beam_i]=1./beam_use[beam_i]
mask_use=fltarr(dimension,elements) & mask_use[beam_i]=1.
WorkImage=Image*beam_corr^2.
image_max=max_filter(WorkImage,radius=LocalMaxRadius,/circle,/edge,/double,mask=mask_use)
;image_min=max_filter(WorkImage,2*LocalMaxRadius,/circle,/edge,/min,/double)
smooth_image=max_filter(WorkImage,radius=2*LocalMaxRadius,/circle,/edge,/median,/double,mask=mask_use)
source_find_image=WorkImage-smooth_image

source_find_image*=mask_use
sigma=Stddev(source_find_image[beam_i])
nsigma=1.

candidates=where((WorkImage EQ image_max) AND (source_find_image GT sigma*nsigma),n_candidates)
order=Reverse(Sort(source_find_image[candidates]))
candidates=candidates[order]
xvals=Float(candidates mod dimension)
yvals=Float(Floor(candidates/dimension))
flux=source_find_image[candidates]

dist_array=fltarr(n_candidates)
ellipticity_array=fltarr(n_candidates)
amplitude_array=fltarr(n_candidates)
flag=fltarr(n_candidates)
;edge_dist=4.

!Quiet=1
FOR c_i=0.,n_candidates-1 DO BEGIN
;    IF Abs(xvals[c_i]-dimension)<xvals[c_i] LE edge_dist THEN CONTINUE
;    IF Abs(yvals[c_i]-elements)<yvals[c_i] LE edge_dist THEN CONTINUE
    x_low=Round(xvals[c_i]-fit_radius)>0 & x_high=Round(xvals[c_i]+fit_radius)<(dimension-1)
    y_low=Round(yvals[c_i]-fit_radius)>0 & y_high=Round(yvals[c_i]+fit_radius)<(elements-1)
    sub_image=source_find_image[x_low:x_high,y_low:y_high]
;    sub_image_smooth=smooth_image[x_low:x_high,y_low:y_high];+meanval
    sub_mask=mask_use[x_low:x_high,y_low:y_high]
    sub_image=sub_image*sub_mask;+sub_image_smooth*(1-sub_mask)
    sub_dimension=x_high-x_low+1
    sub_elements=y_high-y_low+1
    center_x=xvals[c_i]-x_low
    center_y=yvals[c_i]-y_low
    fwhm_use=localmaxradius
    gcntrd, sub_image, center_x, center_y, xcen, ycen, fwhm_use, /SILENT
    IF xcen EQ -1 THEN BEGIN
        flag[c_i]=1
        CONTINUE
    ENDIF
    single_fit=gauss2Dfit(sub_image,single_params,/tilt)
    amplitude=single_params[1]
    width_x=single_params[2]
    width_y=single_params[3]
    single_ellipticity=(width_x/width_y)>(width_y/width_x)

    d=Sqrt((center_x-xcen)^2.+(center_y-ycen)^2.)
    dist_array[c_i]=d
    ellipticity_array[c_i]=single_ellipticity
    amplitude_array[c_i]=amplitude
    IF d GE dist_tol THEN flag[c_i]=+2
    IF (single_ellipticity GE ellipticity_threshold) OR (single_ellipticity LE 0) THEN flag[c_i]+=4
    IF amplitude GE 10.*sub_image[center_x,center_y] THEN flag[c_i]+=8
    IF amplitude LE 0 THEN flag[c_i]+=8
    xvals[c_i]=xvals[c_i]-center_x+xcen
    yvals[c_i]=yvals[c_i]-center_y+ycen
    flux[c_i]=amplitude
ENDFOR
!Quiet=0

IF Total(flag) NE 0 THEN BEGIN
    keep_i=where(flag EQ 0,n_keep)
    order_keep=Reverse(Sort(flux[keep_i]))
    IF Keyword_Set(max_sources) THEN n_keep=n_keep<max_sources
    order=(keep_i[order_keep])[0:n_keep-1]
ENDIF ELSE BEGIN
    order=Reverse(Sort(flux))
    IF Keyword_Set(max_sources) THEN n_keep=n_candidates<max_sources
    order=order[0:n_keep-1]
ENDELSE

xvals=xvals[order]
yvals=yvals[order]
flux=flux[order]
candidates=candidates[order]
dist_array=dist_array[order]
ellipticity_array=ellipticity_array[order]
;flag=flag[order]

;columns of source_array are: 0 x, 1 y, 2 UV flux, 3 image flux
source_array=fltarr(4,n_keep)
source_array[0,*]=xvals-dimension/2.
source_array[1,*]=yvals-elements/2.
source_array[2,*]=flux/(dimension*elements)
source_array[3,*]=flux
;source_array[5,*]=candidates
;source_array[6,*]=

timing=Systime(1)-t0
RETURN,source_array
END