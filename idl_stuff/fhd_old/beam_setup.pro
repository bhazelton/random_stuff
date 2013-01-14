PRO beam_setup,data_directory=data_directory,filename=filename,restore_last=restore_last,$
    residual_tolerance=residual_tolerance,residual_threshold=residual_threshold,convolve_pad=convolve_pad,$
    file_identifier=file_identifier
;restore_last skips re-calculating beam parameters, 
; but does not skip loading visibility parameters at the start

COMMON obs_params,dimension,elements,degpix,kbinsize,azimuth_center,elevation_center,rotation
COMMON psf_params,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals
COMMON visibility_params,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,freq_bin_i,file_path
COMMON switches,source_switch,psf_switch

;Fixed parameters (though they will be over-written if /restore_last and the actual parameters are different)
psf_switch=1
dimension=(elements=1024.) ;=2048?
kx_span=512d ;Units are # of wavelengths
ky_span=kx_span
kbinsize=kx_span/dimension
degpix=!RaDeg/(kbinsize*dimension)
psf_resolution=32. ;=32?
psf_dim=16. ;=16?

;residual_tolerance is residual as fraction of psf_base above which to include 
IF N_Elements(residual_tolerance) EQ 0 THEN residual_tolerance=1./100.  
;residual_threshold is minimum residual above which to include
IF N_Elements(residual_threshold) EQ 0 THEN residual_threshold=0.
IF N_Elements(convolve_pad) EQ 0 THEN convolve_pad=psf_resolution*Round(psf_dim/8.)
IF N_Elements(data_directory) EQ 0 THEN data_directory='DATA\r4\clmw\X14\PicA_121_20100924211938'
IF N_Elements(filename) EQ 0 THEN filename='PicA_121_20100924211938.cal'
IF N_Elements(file_identifier) EQ 0 THEN file_id='00_00' ELSE file_id=Strmid(file_identifier,2,5)
;filename='PicA_121_20100924211938'
IF N_Elements(file_path) EQ 0 THEN file_path=filepath(filename,root_dir=rootdir('mwa'),subdir=data_directory)
n_tiles=32
n_frequencies=768.
n_pol=4
c_light=299792458.

IF not Keyword_Set(restore_last) THEN observation_setup,filename=filename,data_directory=data_directory

header_filename=filename+'_header'
tile_gain_x_filename=filename+'_tile_gains_x'
tile_gain_y_filename=filename+'_tile_gains_y'
textfast,gain_array_X,/read,filename=tile_gain_x_filename,root=rootdir('mwa'),filepathfull=data_directory
textfast,gain_array_Y,/read,filename=tile_gain_y_filename,root=rootdir('mwa'),filepathfull=data_directory
FitsFast,uvw_baseline,data_header,/read,filename=header_filename,root=rootdir('mwa'),filepathfull=data_directory
;;textfast,tile_locs,/read,filename='32T_tile_locations',filepathfull='guides',root=rootdir('mwa'),first_line=1,column=indgen(3)+2
textfast,tile_locs,/read,filename='32T_tile_locations',root=rootdir('mwa'),first_line=1,column=indgen(3)+2
xtile=reform(tile_locs[0,*])
ytile=reform(tile_locs[1,*])
ztile=reform(tile_locs[2,*])
tile_dist0=fltarr(32,32)
FOR i=0,31 DO FOR j=i,31 DO tile_dist0[i,j]=Sqrt((xtile[i]-xtile[j])^2.+(ytile[i]-ytile[j])^2.)

uu_arr=reform(uvw_baseline[0,*])
vv_arr=reform(uvw_baseline[1,*])
ww_arr=reform(uvw_baseline[2,*])
baseline_arr=reform(uvw_baseline[3,*])
date_arr=reform(uvw_baseline[4,*])

tile_A=Floor(baseline_arr/256) 
tile_B=Fix(baseline_arr mod 256)
dx_arr=(xtile[tile_A-1]-xtile[tile_B-1])/c_light
dy_arr=(ytile[tile_A-1]-ytile[tile_B-1])/c_light

freq_ref=sxpar(data_header,'crval4') ;1.5424E8
freq_bandwidth=sxpar(data_header,'cdelt4') ;40000 
freq_ref_i=sxpar(data_header,'crpix4') ;368
frequency_array=(Dindgen(n_frequencies)-freq_ref_i)*freq_bandwidth+freq_ref
freq_bin=32.*freq_bandwidth  ;Hz
freq_hist=histogram(frequency_array,locations=freq_center,binsize=freq_bin,reverse_ind=freq_ri)
nfreq_bin=N_Elements(freq_hist)
freq_bin_i=fltarr(n_frequencies)
FOR bin=0,nfreq_bin-1 DO IF freq_ri[bin] LT freq_ri[bin+1] THEN freq_bin_i[freq_ri[freq_ri[bin]:freq_ri[bin+1]-1]]=bin


  kx_arr0=uu_arr#frequency_array
  ky_arr0=vv_arr#frequency_Array
  print,minmax(kx_Arr0)
  print, minmax(ky_arr0)
 stop

b0i=Uniq(date_arr)
nb=N_Elements(b0i)
bin_width=fltarr(nb)
bin_width[0]=b0i[0]+1
FOR i=1,nb-1 DO bin_width[i]=b0i[i]-b0i[i-1]
bin_offset=fltarr(nb) & bin_offset[1:*]=total(bin_width[0:nb-2],/cumulative)
nbaselines=bin_offset[1]

COMMON SITE, lat, lon, tzone
lat=sxpar(data_header,'lat') ;-26.703319 ;latitude (degrees)
lon=sxpar(data_header,'lon') ;116.67081 ;longitude (degrees)
alt=sxpar(data_header,'alt') ;377.83 ;altitude (meters)
tzone=0 ;assume for now that we use UTC

obsra=sxpar(data_header,'obsra')
obsdec=sxpar(data_header,'obsdec')
date_obs=sxpar(data_header,'date-obs')
year=Float(Strmid(date_obs,0,4))
month=Float(Strmid(date_obs,5,2))
day=Float(Strmid(date_obs,8,2))
Juldate,[year,month,day],jdate & jdate+=2400000.5 ;conversion between two forms of Julian dates
Julian_time=Jdate+date_arr-0.5 ;Julian dates are zeroed at noon, not midnight
zenpos,Julian_time,zen_ra,zen_dec
;angle1=angle_difference(zen_dec*!Radeg,zen_ra*!Radeg,obsdec,obsra,/degree)

;foreshorten_U=Cos(zen_ra-obsra*!DtoR)
;foreshorten_V=Cos(zen_dec-obsdec*!DtoR)
;angle_offset=obsra*!DtoR-zen_ra
zenith_angle_offset=Angle_difference(Median(zen_dec)*!Radeg,Median(zen_ra)*!Radeg,obsdec,obsra,/degree,/nearest)
azimuth_angle_offset=Median(zen_ra)*!Radeg-obsra
azimuth_center=180.+azimuth_angle_offset
elevation_center=90.-zenith_angle_offset
rotation=azimuth_angle_offset

;test_u=(dx_arr*Cos(angle_offset)-dy_arr*Sin(angle_offset))*foreshorten_U
;test_v=(dy_arr*Cos(angle_offset)+dx_arr*Sin(angle_offset))*foreshorten_V

IF Keyword_Set(restore_last) THEN BEGIN
    restore,file_path+'_beams'+file_id+'.sav'
    RETURN
ENDIF

;begin forming psf
psf_residuals_i=Ptrarr(n_pol,nfreq_bin,nbaselines,/allocate) ;contains arrays of pixel indices of pixels with modified psf for a given baseline id
psf_residuals_val=Ptrarr(n_pol,nfreq_bin,nbaselines,/allocate) ;contains arrays of values corresponding to the pixel indices above
psf_residuals_n=fltarr(n_pol,nfreq_bin,nbaselines) ;contains the total number of modified pixels for each baseline id

psf_base=Ptrarr(n_pol,nfreq_bin,psf_resolution,psf_resolution,/allocate)
psf_xvals=Ptrarr(psf_resolution,psf_resolution,/allocate)
psf_yvals=Ptrarr(psf_resolution,psf_resolution,/allocate)
xvals_i=meshgrid(psf_dim,psf_dim,1)*psf_resolution
yvals_i=meshgrid(psf_dim,psf_dim,2)*psf_resolution
xvals=meshgrid(psf_dim*psf_resolution,psf_dim*psf_resolution,1)/psf_resolution-psf_dim/2.
yvals=meshgrid(psf_dim*psf_resolution,psf_dim*psf_resolution,2)/psf_resolution-psf_dim/2.

FOR i=0,psf_resolution-1 DO FOR j=0,psf_resolution-1 DO BEGIN 
    *psf_xvals[i,j]=xvals[xvals_i+i,yvals_i+j]
    *psf_yvals[i,j]=yvals[xvals_i+i,yvals_i+j]
ENDFOR

;polarization ids are 0:XX, 1:YY, 2:XY, 3:YX
;fshort_U=median(foreshorten_U)
;fshort_V=median(foreshorten_V)
;ang_off=median(angle_offset)*!Radeg
gain_tile_i=reform(gain_array_X[0,*])
gain_freq_bin_i=findgen(N_Elements(gain_tile_i)) mod nfreq_bin
pol_arr=[[1,1],[2,2],[1,2],[2,1]]
t1=Systime(1)
FOR pol_i=0,n_pol-1 DO BEGIN
    pol1=pol_arr[0,pol_i]
    pol2=pol_arr[1,pol_i]
    gain1_full=(pol1 EQ 1) ? gain_array_X:gain_array_Y
    gain2_full=(pol2 EQ 1) ? gain_array_X:gain_array_Y
    
    FOR freq_i=0,nfreq_bin-1 DO BEGIN
        IF freq_i mod 4 EQ 3 THEN BEGIN
            t1b=Systime(1)-t1
            iter_past=pol_i*nfreq_bin*nbaselines+freq_i*nbaselines
            iter_total=n_pol*nfreq_bin*nbaselines
            print,Strcompress(String(format='("Time elapsed:",I," estimated time remaining:",I)',t1b,t1b*(iter_total-iter_past)/iter_past))
        ENDIF
        antenna_beam_arr1=Ptrarr(16,/allocate)
        antenna_beam_arr2=Ptrarr(16,/allocate)
        beam1_arr=Ptrarr(n_tiles,/allocate)
        beam2_arr=Ptrarr(n_tiles,/allocate)
        
        gain1=gain1_full[1:*,where(gain_freq_bin_i EQ freq_i)]
        gain2=gain2_full[1:*,where(gain_freq_bin_i EQ freq_i)]
        gain1_avg=Median(gain1,dimension=1)
        gain2_avg=Median(gain2,dimension=1)
        beam_test1=tile_beam_generate(fltarr(16)+1.,antenna_beam_arr1,$
            frequency=freq_center[freq_i],polarization=pol1,$
            psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
            zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
            normalization=1)
        IF pol2 EQ pol1 THEN antenna_beam_arr2=antenna_beam_arr1
        beam_test2=tile_beam_generate(fltarr(16)+1.,antenna_beam_arr2,$
            frequency=freq_center[freq_i],polarization=pol2,$
            psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
            zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
            normalization=1)
        gain1_normalization=1./(Total(beam_test1)/psf_resolution^2.)
        gain2_normalization=1./(Total(beam_test2)/psf_resolution^2.)
        beam1_0=tile_beam_generate(gain1_avg,antenna_beam_arr1,$
            frequency=freq_center[freq_i],polarization=pol1,$
            psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
            zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
            normalization=gain1_normalization)
        beam2_0=tile_beam_generate(gain2_avg,antenna_beam_arr2,$
            frequency=freq_center[freq_i],polarization=pol2,$
            psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
            zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
            normalization=gain2_normalization)
        psf_base0=Convolve2(beam1_0,beam2_0,pad=convolve_pad,/absolute)
        
        FOR tile_i=0,n_tiles-1 DO BEGIN
            *beam1_arr[tile_i]=tile_beam_generate(gain1[*,tile_i],antenna_beam_arr1,$
                frequency=freq_center[freq_i],polarization=pol1,$
                psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
                zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
                normalization=gain1_normalization)
            *beam2_arr[tile_i]=tile_beam_generate(gain2[*,tile_i],antenna_beam_arr2,$
                frequency=freq_center[freq_i],polarization=pol2,$
                psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
                zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
                normalization=gain2_normalization)
        ENDFOR
        
        FOR bi=0,nbaselines-1 DO BEGIN
            IF Min((gain1[*,tile_A[bi]-1]-gain1_avg EQ fltarr(16)) AND (gain2[*,tile_B[bi]-1]-gain2_avg EQ fltarr(16))) THEN BEGIN
                psf_residuals_n[pol_i,freq_i,bi]=0
                CONTINUE
            ENDIF
;;            beam1=tile_beam_generate(gain1[*,tile_A[bi]],antenna_beam_arr1,$
;;                frequency=freq_center[freq_i],angle_offset=ang_off,polarization=pol1,$
;;                psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
;;                foreshorten_U=fshort_U,foreshorten_V=fshort_V)
;;            beam2=tile_beam_generate(gain2[*,tile_B[bi]],antenna_beam_arr2,$
;;                frequency=freq_center[freq_i],angle_offset=ang_off,polarization=pol2,$
;;                psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
;;                foreshorten_U=fshort_U,foreshorten_V=fshort_V)
;            beam1=*beam1_arr[tile_A[bi]-1]
;            beam2=*beam2_arr[tile_B[bi]-1]
            psf_single=Convolve2(*beam1_arr[tile_A[bi]-1],*beam2_arr[tile_B[bi]-1],pad=convolve_pad,/absolute)
            residual_single=psf_single-psf_base0
            i_res=where(residual_single GE ((psf_base0*residual_tolerance)>residual_threshold),nres)
            psf_residuals_n[pol_i,freq_i,bi]=nres
            IF nres GT 0 THEN BEGIN
                *psf_residuals_i[pol_i,freq_i,bi]=i_res
                *psf_residuals_val[pol_i,freq_i,bi]=residual_single[i_res]
            ENDIF
        ENDFOR
        Ptr_free,antenna_beam_arr1,antenna_beam_arr2,beam1_arr,beam2_arr
        FOR i=0,psf_resolution-1 DO FOR j=0,psf_resolution-1 DO $
            *psf_base[pol_i,freq_i,i,j]=psf_base0[xvals_i+i,yvals_i+j];/Total(psf_base0[xvals_i+i,yvals_i+j]) 
    ENDFOR
ENDFOR

save,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals,$
    filename=file_path+'_beams'+file_id+'.sav'

END
