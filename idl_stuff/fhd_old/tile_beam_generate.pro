FUNCTION tile_beam_generate,antenna_gain_arr,antenna_beam_arr, $
    frequency=frequency,angle_offset=angle_offset,polarization=polarization,$
    psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
    foreshorten_U=foreshorten_U,foreshorten_V=foreshorten_V,normalization=normalization,$
    xvals=xvals,yvals=yvals,zenith_angle_offset=zenith_angle_offset,window=window,azimuth_angle_offset=azimuth_angle_offset,$
    dimension=dimension,elements=elements

;indices of antenna_gain_arr correspond to these antenna locations:
;12 13 14 15
;8  9  10 11
;4  5  6  7
;0  1  2  3 
;
;polarization 1: x, 2: y
;angle offset is the rotation of the entire tile in current coordinates in DEGREES
; (this should be the rotation between E-W or N-S and Ra-Dec)
;IF N_Elements(psf_dim) EQ 0 THEN psf_dim=16
;IF N_Elements(psf_resolution) EQ 0 THEN psf_resolution=32.
;IF N_Elements(kbinsize) EQ 0 THEN kbinsize=0.5
kbinsize_use=kbinsize/psf_resolution
;;dimension=psf_dim
dimension = psf_dim * psf_resolution ;; Bryna -- prevent chopping in image space
;IF N_Elements(frequency) EQ 0 THEN frequency=1.36E8
;IF N_Elements(polarization) EQ 0 THEN polarization=1

IF N_Elements(foreshorten_U) EQ 0 THEN foreshorten_U=1.
IF N_Elements(foreshorten_V) EQ 0 THEN foreshorten_V=1.
IF N_Elements(normalization) EQ 0 THEN normalization=1.
IF N_Elements(zenith_angle_offset) EQ 0 THEN za=0. ELSE za=zenith_angle_offset
IF N_Elements(azimuth_angle_offset) EQ 0 THEN az=0. ELSE az=azimuth_angle_offset
dimension=Float(dimension)
IF N_Elements(elements) EQ 0 THEN elements=dimension ELSE elements=Float(elements)

antenna_spacing=1.1 ;meters (design) ;1.071
antenna_length=29.125*2.54/100. ;meters (measured)
antenna_height=0.35 ;meters (rumor)
dipole_phase=0.
;IF polarization EQ 1 THEN dipole_phase=0. ELSE dipole_phase=1.
;IF polarization EQ 2 THEN az-=90.

;kbinsize_use=kbinsize/psf_resolution
;degpix=!RaDeg/(kbinsize_use*dimension*psf_resolution)
degpix=!RaDeg/(kbinsize_use*dimension) ;; Bryna
;kconv= (2*!Pi)*(frequency/299792458.)/kbinsize_use;factor to convert meters to pixels
Kconv=(2*!Pi)*(frequency/299792458.) ;wavenumber (radians/meter)
wavelength=299792458./frequency

IF Keyword_Set(antenna_beam_arr) THEN IF Keyword_Set(*antenna_beam_arr[0]) THEN BEGIN
    ;;tile_beam=fltarr(dimension*psf_resolution,dimension*psf_resolution)
    tile_beam=fltarr(dimension,elements) ;; Bryna
    FOR i=0,15 DO tile_beam+=*antenna_beam_arr[i]*antenna_gain_arr[i]
    tile_beam*=normalization
    RETURN,tile_beam
ENDIF

xc_arr=Reform((meshgrid(4,4,1)-1.5)*antenna_spacing,16) ;dipole east position (meters)
yc_arr=Reform((meshgrid(4,4,2)-1.5)*antenna_spacing,16) ;dipole north position (meters)
zc_arr=fltarr(16)

;;xvals=(meshgrid(dimension,elements,1)-dimension/2.)*degpix
;;yvals=(meshgrid(dimension,elements,2)-elements/2.)*degpix
xvals=(meshgrid(dimension,elements,1)-dimension/2.)*degpix_use ;; Bryna
yvals=(meshgrid(dimension,elements,2)-elements/2.)*degpix_use ;; Bryna
;za_arr=za+Reform((xvals*Sin(az*!DtoR)+yvals*Cos(az*!DtoR)),dimension*elements)
;az_arr=az+Reform((xvals*Cos(az*!DtoR)-yvals*Sin(az*!DtoR)),dimension*elements)
term_A=Tan(az*!DtoR)
term_B=za
xc=Sqrt((term_B^2.)/(1+term_A^2.))
yc=term_A*xc
za_arr=Reform(Sqrt((xvals-xc)^2.+(yvals-yc)^2.),dimension*elements)
az_arr=Reform(Atan(yvals-yc,xvals-xc),dimension*elements)*!Radeg

;beamformer phase setting (meters)
D0_d=xc_arr*sin(za*!DtoR)*Sin(az*!DtoR)+yc_arr*Sin(za*!DtoR)*Cos(az*!DtoR)

proj_east=Reform(xvals,dimension*elements);Sin(za_arr*!DtoR)*Sin(az_arr*!DtoR)
proj_north=Reform(yvals,dimension*elements);Sin(za_arr*!DtoR)*Cos(az_arr*!DtoR)
proj_z=Cos(za_arr*!DtoR)

;phase of each dipole for the source (relative to the beamformer settings)
D_d=proj_east#xc_arr+proj_north#yc_arr+proj_z#zc_arr-replicate(1,dimension*elements)#D0_d
D_d=Reform(D_d,dimension,elements,16)

groundplane=2*Sin(Cos(za_arr*!DtoR)#(Kconv*(antenna_height+zc_arr))) ;looks correct
groundplane=Reform(groundplane,dimension,elements,16)

projection=Sqrt(1-(Sin(za_arr*!DtoR)*Sin(az_arr*!DtoR))^2.)#replicate(1,16)
projection=Reform(projection,dimension,elements,16)

;leakage_xtoy=0.
;leakage_ytox=0.

ii=Complex(0,1)
IF polarization EQ 1 THEN pol=Cos(xvals*!DtoR)^2. ELSE pol=Cos(yvals*!DtoR)^2.
dipole_gain_arr=Exp(-ii*dipole_phase)*groundplane*projection*Exp(-ii*Kconv*D_d*!DtoR)
horizon_test=where(abs(za_arr) GE 90.,n_horizon_test)
IF n_horizon_test GT 0 THEN BEGIN
    horizon_mask=fltarr(dimension,elements)+1
    horizon_mask[horizon_test]=0
    FOR i=0,15 DO dipole_gain_arr[*,*,i]*=horizon_mask
ENDIF

;; make a uv mask to leave only the part that is detectable by the tile
tile_length = 3.*antenna_spacing + antenna_length
tile_length_uvpix = ceil(tile_length * Kconv / kbinsize_use * 2)/2
uv_mask = fltarr(dimension, elements)
uv_mask[(dimension - tile_length_uvpix)/2+1:(dimension + tile_length_uvpix)/2, $
        (dimension - tile_length_uvpix)/2+1:(dimension + tile_length_uvpix)/2] = 1

IF not Keyword_Set(antenna_beam_arr) THEN antenna_beam_arr=Ptrarr(16,/allocate)
FOR i=0,15 DO BEGIN
    gain_pad=fltarr(psf_dim*psf_resolution,psf_dim*psf_resolution)
    ;;gain_pad[psf_dim*psf_resolution/2-psf_dim/2:psf_dim*psf_resolution/2+psf_dim/2-1,psf_dim*psf_resolution/2-psf_dim/2:psf_dim*psf_resolution/2+psf_dim/2-1]$
    gain_pad=dipole_gain_arr[*,*,i]*pol ;; Bryna
    kgain=fft_shift(FFT(fft_shift(gain_pad)))
    kgain=real_part(kgain) * uv_mask
    *antenna_beam_arr[i]=kgain
ENDFOR

;;tile_beam=fltarr(dimension*psf_resolution,dimension*psf_resolution)
tile_beam=fltarr(dimension,elements) ;; Bryna
FOR i=0,15 DO tile_beam+=*antenna_beam_arr[i]*antenna_gain_arr[i]
tile_beam*=normalization

RETURN,tile_beam

;gain=Total(real_part(dipole_gain_arr),3)*pol
;gain_pad=fltarr(psf_dim*psf_resolution,psf_dim*psf_resolution)
;gain_pad[psf_dim*psf_resolution/2-psf_dim/2:psf_dim*psf_resolution/2+psf_dim/2-1,psf_dim*psf_resolution/2-psf_dim/2:psf_dim*psf_resolution/2+psf_dim/2-1]=gain
;kgain=fft_shift(FFT(fft_shift(gain_pad)))
;kgain=real_part(kgain)
;kgain/=max(kgain)
;gain/=max(gain)
;filename=String(format='("Beam at f= ",I3,"MHz")',Round(frequency/1E6))
;imagefast,gain,data_dir='temporary',proj='mwa',/right,low=min(gain),high=max(gain),filename=filename
;filename=String(format='("Beam (K) at f= ",I3,"MHz")',Round(frequency/1E6))
;imagefast,kgain,data_dir='temporary',proj='mwa',/right,low=min(kgain),high=max(kgain),filename=filename
END




;kconv= 1./(2*!Pi*wavelength*kbinsize_use);factor to convert meters to pixels
;;kconv= 1./(wavelength*kbinsize_use);factor to convert meters to pixels
;ant_spc=antenna_spacing*kconv
;ant_len=antenna_length*kconv
;;ant_ht=antenna_height*kconv
;;ant_ratio=(299792458./frequency)*antenna_length
;ant_ratio=ant_len*kbinsize_use
;psf_dim2=1.5*psf_dim ;use a larger size box to allow cropping after rotation
;
;;generate the beam in unrotated coordinates
;psf_cen=Floor(psf_dim2*psf_resolution/2.)
;xc_arr=(meshgrid(4,4,1)-1.5)*ant_spc+psf_cen
;yc_arr=(meshgrid(4,4,2)-1.5)*ant_spc+psf_cen
;;uc_arr=(xc_arr*Cos(angle_offset*!DtoR)-yc_arr*Sin(angle_offset*!DtoR))*foreshorten_U
;;vc_arr=(xc_arr*Sin(angle_offset*!DtoR)+yc_arr*Cos(angle_offset*!DtoR))*foreshorten_V
;
;xvals=meshgrid(psf_dim2*psf_resolution,psf_dim2*psf_resolution,1)
;yvals=meshgrid(psf_dim2*psf_resolution,psf_dim2*psf_resolution,2)
;;rvals=Sqrt((xvals-psf_cen)^2.+(yvals-psf_cen)^2.)
;;uvals=(xvals*Cos(angle_offset*!DtoR)-yvals*Sin(angle_offset*!DtoR))*foreshorten_U
;;vvals=(xvals*Sin(angle_offset*!DtoR)+yvals*Cos(angle_offset*!DtoR))*foreshorten_V
;;border_pix=where((xvals EQ 0) OR (yvals EQ 0) OR (xvals EQ psf_dim2*psf_resolution-1) OR (yvals EQ psf_dim2*psf_resolution-1))
;
;tile_beam_base=fltarr(psf_dim2*psf_resolution,psf_dim2*psf_resolution)
;FOR antenna=0,15 DO BEGIN
;    IF Arg_present(antenna_beam_arr) THEN IF Keyword_Set(*antenna_beam_arr[antenna]) THEN BEGIN
;        ant_beam=*antenna_beam_arr[antenna] 
;        ant_calc=0
;    ENDIF ELSE ant_calc=1 ELSE ant_calc=1
;    IF Keyword_Set(ant_calc) THEN BEGIN
;        xvals1=xvals-xc_arr[antenna]
;        yvals1=yvals-yc_arr[antenna]
;;        uvals1=xvals-uc_arr[antenna]
;;        vvals1=yvals-vc_arr[antenna]
;;        rvals1=Sqrt(uvals1^2.+vvals1^2.)
;;        angle1=Atan(vvals1,uvals1)
;        
;        rvals1=Sqrt(xvals1^2.+yvals1^2.)>1.;/ant_len
;        angle1=Atan(yvals1,xvals1)
;        IF polarization EQ 2 THEN angle1-=!Pi/2.
;;        angle_term=Sin(angle1)^2.
;;        angle_term2=Abs(Cos(2.*Cos(angle1)*(2*antenna_height/(wavelength*Cos(zenith_angle_offset*!DtoR)))))
;;        ant_beam=real_part(fft_shift(FFT(fft_shift(angle_term*angle_term2))))
;        angle_term=(Cos((!Pi/2.)*Cos(angle1))/Sin(angle1))^2.
;
;        ant_beam=angle_term/rvals1^2.
;;        ant_beam=2.*Abs(ant_beam)*Abs(Cos((!Pi*ant_ratio)*Sin(angle1))) ;this term is wrong!
;;        ant_beam=2.*Abs(ant_beam)*Abs(Sin((!Pi*ant_ratio)*Sin(angle1))) ;this term is wrong!
;;        ant_beam=2.*Abs(ant_beam)*Abs(Sin((!Pi*antenna_height/wavelength/Sin(zenith_angle_offset*!DtoR))*Sin(angle1)))
;;        ant_beam=2.*Abs(ant_beam)*Abs(Cos(2.*Sin(angle1)*(2*antenna_height/(wavelength*Cos(zenith_angle_offset*!DtoR)))))
;        ant_beam=2.*Abs(ant_beam)*Abs(Cos(2.*Cos(angle1)*(2*antenna_height/(wavelength*Cos(zenith_angle_offset*!DtoR)))))
;;        ant_beam-=Max(ant_beam[border_pix])
;;        ant_beam=ant_beam>0
;        IF Arg_present(antenna_beam_arr) THEN *antenna_beam_arr[antenna]=ant_beam
;    ENDIF
;;    tile_beam_base+=antenna_gain_arr[antenna]*ant_beam2
;    tile_beam_base+=antenna_gain_arr[antenna]*ant_beam
;ENDFOR
;;tile_beam=tile_beam_base
;tile_beam=rot(tile_beam_base,-angle_offset,missing=0,/interp)
;psf_cen0=Floor(psf_dim*psf_resolution/2.)
;tile_beam=congrid(tile_beam,psf_dim2*psf_resolution*foreshorten_U,psf_dim2*psf_resolution*foreshorten_V,/center)
;tile_beam=tile_beam[Round(psf_cen*foreshorten_U)-psf_cen0:Round(psf_cen*foreshorten_U)-psf_cen0+psf_dim*psf_resolution-1,Round(psf_cen*foreshorten_V)-psf_cen0:Round(psf_cen*foreshorten_V)-psf_cen0+psf_dim*psf_resolution-1]
;;tile_beam=tile_beam[psf_cen-psf_cen0:psf_cen-psf_cen0+psf_dim*psf_resolution-1,psf_cen-psf_cen0:psf_cen-psf_cen0+psf_dim*psf_resolution-1]
;IF Keyword_Set(window) THEN tile_beam*=Hanning(psf_dim*psf_resolution,psf_dim*psf_resolution)
;tile_beam*=normalization
;
;RETURN,tile_beam 
;END
