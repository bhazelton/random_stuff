FUNCTION mwa_tile_beam_generate,antenna_gain_arr,antenna_beam_arr,$
    frequency=frequency,angle_offset=angle_offset,polarization=polarization,$
    psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
    foreshorten_U=foreshorten_U,foreshorten_V=foreshorten_V,normalization=normalization,$
    xvals=xvals,yvals=yvals,zenith_angle_offset=zenith_angle_offset,window=window,azimuth_angle_offset=azimuth_angle_offset,$
    dimension=dimension,elements=elements,za_arr=za_arr,az_arr=az_arr

compile_opt idl2,strictarrsubs  
;indices of antenna_gain_arr correspond to these antenna locations:
;12 13 14 15
;8  9  10 11
;4  5  6  7
;0  1  2  3 
;
;polarization 1: x, 2: y
;angle offset is the rotation of the entire tile in current coordinates in DEGREES
; (this should be the rotation between E-W or N-S and Ra-Dec)
kbinsize_use=kbinsize/psf_resolution
;dimension=psf_dim
;IF N_Elements(frequency) EQ 0 THEN frequency=1.36E8
;IF N_Elements(polarization) EQ 0 THEN polarization=1

;IF N_Elements(foreshorten_U) EQ 0 THEN foreshorten_U=1.
;IF N_Elements(foreshorten_V) EQ 0 THEN foreshorten_V=1.
IF N_Elements(normalization) EQ 0 THEN normalization=1.
;IF N_Elements(zenith_angle_offset) EQ 0 THEN za=0. ELSE za=zenith_angle_offset
;IF N_Elements(azimuth_angle_offset) EQ 0 THEN az=0. ELSE az=azimuth_angle_offset
psf_dim=Float(psf_dim)
psf_dim2=psf_dim*psf_resolution
za=za_arr[psf_dim2/2.,psf_dim2/2.]
az=az_arr[psf_dim2/2.,psf_dim2/2.]
;IF N_Elements(elements) EQ 0 THEN elements=dimension ELSE elements=Float(elements)

antenna_spacing=1.1 ;meters (design) ;1.071
antenna_length=29.125*2.54/100. ;meters (measured)
;groundscreen_width=5. ;meters (guess!)
antenna_height=0.35 ;meters (rumor)

;IF polarization EQ 2 THEN az-=90.

;kbinsize_use=kbinsize/psf_resolution
degpix_use=!RaDeg/(kbinsize_use*psf_dim2) 
;kconv= (2*!Pi)*(frequency/299792458.)/kbinsize_use;factor to convert meters to pixels
Kconv=(2*!Pi)*(frequency/299792458.) ;wavenumber (radians/meter)
wavelength=299792458./frequency

IF Keyword_Set(antenna_beam_arr) THEN IF Keyword_Set(*antenna_beam_arr[0]) THEN BEGIN
    tile_beam=fltarr(psf_dim2,psf_dim2)
    FOR i=0,15 DO tile_beam+=*antenna_beam_arr[i]*antenna_gain_arr[i]
    
;    uv_mask=fltarr(psf_dim2,psf_dim2)
;    threshold=Max(tile_beam)/1E6
;    beam_i=region_grow(tile_beam,psf_dim2*(1.+psf_dim2)/2.,thresh=[Max(tile_beam)/1E6,Max(tile_beam)])
;    uv_mask[beam_i]=1.
    
    tile_beam*=normalization;*uv_mask
    tile_beam=tile_beam
    RETURN,tile_beam
ENDIF

xc_arr=Reform((meshgrid(4,4,1)-1.5)*antenna_spacing,16) ;dipole east position (meters)
yc_arr=Reform((meshgrid(4,4,2)-1.5)*antenna_spacing,16) ;dipole north position (meters)
zc_arr=fltarr(16)

xvals=(meshgrid(psf_dim2,psf_dim2,1)-psf_dim2/2.)*degpix_use
yvals=(meshgrid(psf_dim2,psf_dim2,2)-psf_dim2/2.)*degpix_use
;za_arr=za+Reform((xvals*Sin(az*!DtoR)+yvals*Cos(az*!DtoR)),dimension*elements)
;az_arr=az+Reform((xvals*Cos(az*!DtoR)-yvals*Sin(az*!DtoR)),dimension*elements)
term_A=Tan(az*!DtoR)
term_B=za*!DtoR
xc=Sqrt((term_B^2.)/(1+term_A^2.))
yc=term_A*xc
za_arr_use=Reform(za_arr,(psf_dim2)^2.)
az_arr_use=Reform(az_arr,(psf_dim2)^2.)
;za_arr=Reform(Sqrt((xvals*!DtoR-xc)^2.+(yvals*!DtoR-yc)^2.),(psf_dim2)^2.)*!Radeg
;az_arr=Reform(Atan(yvals*!DtoR-yc,xvals*!DtoR-xc),(psf_dim2)^2.)*!Radeg

;!!!THIS SHOULD REALLY BE READ IN FROM A FILE!!!
;beamformer phase setting (meters) 
D0_d=xc_arr*sin(za*!DtoR)*Sin(az*!DtoR)+yc_arr*Sin(za*!DtoR)*Cos(az*!DtoR)

proj_east=Reform(xvals,(psf_dim2)^2.);Sin(za_arr*!DtoR)*Sin(az_arr*!DtoR)
proj_north=Reform(yvals,(psf_dim2)^2.);Sin(za_arr*!DtoR)*Cos(az_arr*!DtoR)
proj_z=Cos(za_arr_use*!DtoR)

;phase of each dipole for the source (relative to the beamformer settings)
D_d=proj_east#xc_arr+proj_north#yc_arr+proj_z#zc_arr-replicate(1,(psf_dim2)^2.)#D0_d
D_d=Reform(D_d,psf_dim2,psf_dim2,16)

groundplane=2*Sin(Cos(za_arr_use*!DtoR)#(Kconv*(antenna_height+zc_arr))) ;looks correct
groundplane=Reform(groundplane,psf_dim2,psf_dim2,16)

IF polarization EQ 1 THEN projection=Sqrt(1-(Sin(za_arr_use*!DtoR)*Sin(az_arr_use*!DtoR))^2.)#replicate(1,16) $
    ELSE projection=Sqrt(1-(Sin(za_arr_use*!DtoR)*Cos(az_arr_use*!DtoR))^2.)#replicate(1,16) 
projection=Reform(projection,psf_dim2,psf_dim2,16)

;leakage_xtoy=0.
;leakage_ytox=0.

ii=Complex(0,1)
;IF polarization EQ 2 THEN pol=Cos(xvals*!DtoR-xc)^2. ELSE pol=Cos(yvals*!DtoR-yc)^2.
IF polarization EQ 1 THEN pol=Cos(xvals*!DtoR-xc) ELSE pol=Cos(yvals*!DtoR-yc)
;IF polarization EQ 1 THEN pol=(1.-((xvals*!DtoR-xc)^2.)/2.)>0. ELSE pol=(1.-((yvals*!DtoR-yc)^2.)/2.)>0.
dipole_gain_arr=groundplane*projection*Exp(-ii*Kconv*D_d*!DtoR)
horizon_test=where(abs(za_arr_use) GE 90.,n_horizon_test)

horizon_mask=fltarr(psf_dim2,psf_dim2)+1
IF n_horizon_test GT 0 THEN horizon_mask[horizon_test]=0
    
;    FOR i=0,15 DO dipole_gain_arr[*,*,i]*=horizon_mask
;ENDIF

IF not Keyword_Set(antenna_beam_arr) THEN antenna_beam_arr=Ptrarr(16,/allocate)
FOR i=0,15 DO BEGIN
;    gain_pad=fltarr(psf_dim2,psf_dim2)
;    gain_pad[psf_dim2/2-psf_dim/2:psf_dim2/2+psf_dim/2-1,psf_dim2/2-psf_dim/2:psf_dim2/2+psf_dim/2-1]$
;        =dipole_gain_arr[*,*,i]*pol
;    gain_pad=dipole_gain_arr[*,*,i]*pol
;    kgain=fft_shift(FFT(fft_shift(gain_pad)))
;    kgain=real_part(kgain)
;    *antenna_beam_arr[i]=kgain
    *antenna_beam_arr[i]=dipole_gain_arr[*,*,i]*horizon_mask*pol
ENDFOR

tile_beam=fltarr(psf_dim2,psf_dim2)
FOR i=0,15 DO tile_beam+=*antenna_beam_arr[i]*antenna_gain_arr[i]

;uv_mask=fltarr(psf_dim2,psf_dim2)
;threshold=Max(tile_beam)/1E6
;beam_i=region_grow(tile_beam,psf_dim2*(1.+psf_dim2)/2.,thresh=[Max(tile_beam)/1E6,Max(tile_beam)])
;uv_mask[beam_i]=1.

tile_beam*=normalization;*uv_mask

;tile_beam=tile_beam>0.
RETURN,tile_beam

END