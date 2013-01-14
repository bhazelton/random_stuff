

pro sim_beam_setup, data_directory = data_directory, filename = filename, restore_last = restore_last, t32 = t32, $
                    define_baselines = define_baselines, baseline_spacing = baseline_spacing, baseline_layout = baseline_layout, $
                    residual_tolerance = residual_tolerance, residual_threshold = residual_threshold, convolve_pad = convolve_pad, $
                    file_identifier = file_identifier, fine_res = fine_res, use_outliers = use_outliers, $
                    beam_shape_file = beam_shape_file, nfreq_per_image = nfreq_per_image

  COMMON obs_params, dimension, elements, degpix, kbinsize, azimuth_center, elevation_center, rotation
  COMMON psf_params, psf_base, psf_dim, psf_resolution, psf_residuals_i, psf_residuals_val, psf_residuals_n, psf_xvals,psf_yvals
  COMMON visibility_params, uu_arr, vv_arr, frequency_array, date_arr, baseline_arr, bin_offset, freq_bin_i, file_path
  COMMON switches, source_switch, psf_switch

;Fixed parameters (though they will be over-written if /restore_last and the actual parameters are different)
  psf_switch=1
  if keyword_set(t32) or keyword_set(define_baselines) then begin
     dimension=(elements=1024.) ;=2048?
     kx_span=512d              ;Units are # of wavelengths
     ky_span=kx_span
     kbinsize=kx_span/dimension
     degpix=!RaDeg/(kbinsize*dimension)
     psf_resolution=32.         ;=32?
     psf_dim=16.                 ;=16?
  endif else begin
     dimension=(elements=2048.) ;=2048?
     if keyword_set(fine_res) then kx_span = 2048d else kx_span=4096d              ;Units are # of wavelengths
     ky_span=kx_span
     kbinsize=kx_span/dimension
     degpix=!RaDeg/(kbinsize*dimension)
     psf_resolution=64.         ;=32?
     psf_dim=8.                 ;=16?
  endelse
  psf_dim2=psf_dim*psf_resolution
  degpix_use=!RaDeg/(kbinsize*psf_dim)

;residual_tolerance is residual as fraction of psf_base above which to include 
  IF N_Elements(residual_tolerance) EQ 0 THEN residual_tolerance=1./100.  
;residual_threshold is minimum residual above which to include
  IF N_Elements(residual_threshold) EQ 0 THEN residual_threshold=0.
  IF N_Elements(convolve_pad) EQ 0 THEN convolve_pad=psf_resolution*Round(psf_dim/8.)
  IF N_Elements(data_directory) EQ 0 THEN data_directory='data'
  IF N_Elements(filename) EQ 0 THEN filename='Simulation'
  IF N_Elements(file_identifier) EQ 0 THEN file_id='00_00' ELSE file_id=Strmid(file_identifier,2,5)

  IF N_Elements(file_path) EQ 0 THEN file_path=filepath(filename,root_dir=rootdir('mwa'),subdir=data_directory) 
  
  n_frequencies=768.
  n_pol=4
  c_light=299792458.

  if keyword_set(define_baselines) then begin
     case baseline_layout of
        'simple': begin
           spacing = 5*baseline_spacing
           radius = 300

           locs = dindgen(2*radius/spacing+1)*spacing - radius
           wh = where(locs eq 0, count, complement = wh_n0)
           if count gt 0 then locs = locs[wh_n0]
           nlocs = n_elements(locs)
           xlocs = [locs, dblarr(nlocs)]
           ylocs = [dblarr(nlocs), locs]
        end
        'flat': begin
           spacing = 5*baseline_spacing
           radius = 300
           locs = dindgen(2*radius/spacing+1)*spacing - radius
           nlocs = n_elements(locs)
           xlocs = reform(rebin(locs, nlocs, nlocs), nlocs*nlocs)
           ylocs = reform(rebin(reform(locs, 1, nlocs), nlocs, nlocs), nlocs*nlocs)
        end
        'single': begin
           spacing = 5d*baseline_spacing
           xlocs = [spacing]
           ylocs = [0d]
        end
        else: stop
     endcase

     nbaselines = n_elements(xlocs)
          
     uu_arr=(xlocs)/c_light
     vv_arr=(ylocs)/c_light
 
     tile_A = intarr(nbaselines)+1
     tile_B = indgen(nbaselines)+2
     baseline_arr = tile_A*256 + tile_B
     n_tiles = nbaselines+1

  endif else begin
     ;;n_tiles=32
     if keyword_set(t32) then $
        textfast,tile_locs,/read,filename='32T_tile_locations',root=rootdir('mwa'),first_line=1,column=indgen(3)+2 $
     else if not keyword_set(use_outliers) then $
        textfast,tile_locs,/read,filename='496T_tile_locations',root=rootdir('mwa'),first_line=1,column=indgen(2)+1 $
     else textfast,tile_locs,/read,filename='512T_tile_locations',root=rootdir('mwa'),first_line=1,column=indgen(2)+1
     xtile=reform(tile_locs[0,*])
     ytile=reform(tile_locs[1,*])
     ;;ztile=reform(tile_locs[2,*])

     if keyword_set(t32) then n_tiles = 32L else if not keyword_set(use_outliers) then n_tiles = 496L else n_tiles=512L
     if n_elements(xtile) ne n_tiles then message, 'Number of tiles does not match number of given locations.'

     tile_dist0=fltarr(n_tiles,n_tiles)
     FOR i=0,n_tiles-1 DO FOR j=i,n_tiles-1 DO tile_dist0[i,j]=Sqrt((xtile[i]-xtile[j])^2.+(ytile[i]-ytile[j])^2.)
     
     tile_A = reform(rebin(indgen(n_tiles)+1, n_tiles, n_tiles), n_tiles*n_tiles)
     tile_B = reform(rebin(reform(indgen(n_tiles)+1, 1, n_tiles), n_tiles, n_tiles), n_tiles*n_tiles)

     ;; get unique baselines
     temp = double(tile_A)^3/double(tile_B) + double(tile_B)^3/double(tile_A)
     uniq_inds = uniq(temp, sort(temp))
     nbaselines = n_elements(uniq_inds)
     
     if nbaselines ne ((n_tiles^2-n_tiles)/2 + n_tiles) then stop
     
     tile_A = tile_A[uniq_inds]
     tile_B = tile_B[uniq_inds]
     baseline_arr = tile_A*256 + tile_B
     
     uu_arr=(xtile[tile_A-1]-xtile[tile_B-1])/c_light
     vv_arr=(ytile[tile_A-1]-ytile[tile_B-1])/c_light
     ;;ww_arr=(ztile[tile_A-1]-ztile[tile_B-1])/c_light
  endelse
 
  tile_gain_x_filename=filename+'_tile_gains_x'
  tile_gain_y_filename=filename+'_tile_gains_y'

  lat=-26.703319                ;degrees
  lon=116.67081                 ;degrees
  alt=377.83                    ;meters
  tzone=0                       ;assume for now that we use UTC

  if n_elements(nfreq_per_image) eq 0 then begin
     nfreq_bin=24.
     nfreq_per_image = n_frequencies / nfreq_bin
  endif else begin
     nfreq_bin = n_frequencies / nfreq_per_image
  endelse
  gain_array_x=fltarr(17,nfreq_bin*n_tiles)+1.
  gain_array_x[0,*]=Floor(indgen(nfreq_bin*n_tiles)/nfreq_bin)+1
  gain_array_y = gain_array_x

  tile_gain_x_filepath=filepath(tile_gain_x_filename,root_dir=rootdir('mwa'),subdir=data_directory)
  tile_gain_y_filepath=filepath(tile_gain_y_filename,root_dir=rootdir('mwa'),subdir=data_directory)

;do not overwrite a gain_array if one already exists (it's either real data, or the same default data as this!)
  IF file_test(tile_gain_x_filepath) EQ 0 THEN $
     textfast,gain_array_x,/write,filename=tile_gain_x_filename,root=rootdir('mwa'),filepathfull=data_directory
  IF file_test(tile_gain_y_filepath) EQ 0 THEN $
     textfast,gain_array_y,/write,filename=tile_gain_y_filename,root=rootdir('mwa'),filepathfull=data_directory



  freq_ref = 1.5E8
  freq_bandwidth = 40000d
  freq_ref_i = 384
  frequency_array = (Dindgen(n_frequencies)-freq_ref_i)*freq_bandwidth+freq_ref
  freq_bin=nfreq_per_image*freq_bandwidth   ;Hz
  freq_hist=histogram(frequency_array,locations=freq_center,binsize=freq_bin,reverse_ind=freq_ri)
  nfreq_bin=N_Elements(freq_hist)
  freq_bin_i=fltarr(n_frequencies)
  FOR bin=0,nfreq_bin-1 DO IF freq_ri[bin] LT freq_ri[bin+1] THEN freq_bin_i[freq_ri[freq_ri[bin]:freq_ri[bin+1]-1]]=bin

  bin_offset=0

  ;; kx_arr0=uu_arr#frequency_array
  ;; ky_arr0=vv_arr#frequency_Array
  ;; print,minmax(kx_Arr0)
  ;; print, minmax(ky_arr0)
  ;; stop


;; for simulation, point at zenith
  obsra = 0
  obsdec = 0
  zenra = 0
  zendec = 0
  zenith_angle_offset=Angle_difference(zendec,zenra,obsdec,obsra,/degree,/nearest)
  azimuth_angle_offset=zenra-obsra

  ;; these are used in common blocks
  azimuth_center=180.+azimuth_angle_offset
  elevation_center=90.-zenith_angle_offset
  rotation=azimuth_angle_offset


  beam_file = file_path+'_beams'+file_id+'.sav'
  if file_test(beam_file) eq 0 then restore_last = 0

  IF Keyword_Set(restore_last) THEN BEGIN
     restore,file_path+'_beams'+file_id+'.sav'
     RETURN
  ENDIF

;begin forming psf
  psf_residuals_i=Ptrarr(n_pol,nfreq_bin,nbaselines,/allocate) ;contains arrays of pixel indices of pixels with modified psf for a given baseline id
  psf_residuals_val=Ptrarr(n_pol,nfreq_bin,nbaselines,/allocate) ;contains arrays of values corresponding to the pixel indices above
  psf_residuals_n=fltarr(n_pol,nfreq_bin,nbaselines)             ;contains the total number of modified pixels for each baseline id

  psf_base=Ptrarr(n_pol,nfreq_bin,psf_resolution,psf_resolution,/allocate)
  psf_xvals=Ptrarr(psf_resolution,psf_resolution,/allocate)
  psf_yvals=Ptrarr(psf_resolution,psf_resolution,/allocate)
  xvals_i=meshgrid(psf_dim,psf_dim,1)*psf_resolution
  yvals_i=meshgrid(psf_dim,psf_dim,2)*psf_resolution
  xvals=meshgrid(psf_dim2,psf_dim2,1)/psf_resolution-psf_dim/2.
  yvals=meshgrid(psf_dim2,psf_dim2,2)/psf_resolution-psf_dim/2.

  xvals2=(meshgrid(psf_dim2,psf_dim2,1)-psf_dim2/2.)*degpix_use+obsra
  yvals2=(meshgrid(psf_dim2,psf_dim2,2)-psf_dim2/2.)*degpix_use+obsdec
  za_arr=Sqrt((xvals2-zenra)^2.+(yvals2-zendec)^2.)
  az_arr=Atan((yvals2-zendec),(xvals2-zenra))*!RaDeg
  az_arr=360-((az_arr+270) mod 360)
  za_arr0=Sqrt((xvals2-obsra)^2.+(yvals2-obsdec)^2.)
  az_arr0=Atan((yvals2-obsdec),(xvals2-obsra))*!RaDeg
  az_arr0=360-((az_arr0+270) mod 360)

  FOR i=0,psf_resolution-1 DO FOR j=0,psf_resolution-1 DO BEGIN 
     *psf_xvals[i,j]=xvals[xvals_i+i,yvals_i+j]
     *psf_yvals[i,j]=yvals[xvals_i+i,yvals_i+j]
  ENDFOR

  fine_beam = fltarr(psf_dim*psf_resolution, psf_dim*psf_resolution, nfreq_bin, n_pol)

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
        ;;gain1_avg=Median(gain1,dimension=1)
        ;;gain2_avg=Median(gain2,dimension=1)
        gain1_avg=Median(gain1,dimension=2)
        gain2_avg=Median(gain2,dimension=2)
        ;; beam_test1=tile_beam_generate(fltarr(16)+1.,antenna_beam_arr1, $ 
        ;;                               frequency=freq_center[freq_i],polarization=pol1,$
        ;;                               psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
        ;;                               zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
        ;;                               normalization=1)
        ;;IF pol2 EQ pol1 THEN antenna_beam_arr2=antenna_beam_arr1
        ;; beam_test2=tile_beam_generate(fltarr(16)+1.,antenna_beam_arr2,$
        ;;                               frequency=freq_center[freq_i],polarization=pol2,$
        ;;                               psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
        ;;                               zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
        ;;                               normalization=1)
        ;;gain1_normalization=1./(Total(beam_test1)/psf_resolution^2.)
        ;;gain2_normalization=1./(Total(beam_test2)/psf_resolution^2.)
        ;; beam1_0=tile_beam_generate(gain1_avg,antenna_beam_arr1,$
        ;;                            frequency=freq_center[freq_i],polarization=pol1,$
        ;;                            psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
        ;;                            zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
        ;;                            normalization=gain1_normalization)
        beam1_0=mwa_tile_beam_generate(gain1_avg,antenna_beam_arr1,$
                                   frequency=freq_center[freq_i],polarization=pol1,za_arr=za_arr,az_arr=az_arr,$
                                   psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize)
        ;; beam2_0=tile_beam_generate(gain2_avg,antenna_beam_arr2,$
        ;;                            frequency=freq_center[freq_i],polarization=pol2,$
        ;;                            psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
        ;;                            zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
        ;;                            normalization=gain2_normalization)
        IF pol2 EQ pol1 THEN antenna_beam_arr2=antenna_beam_arr1
        beam2_0=mwa_tile_beam_generate(gain2_avg,antenna_beam_arr2,$
                                   frequency=freq_center[freq_i],polarization=pol2,za_arr=za_arr,az_arr=az_arr,$,$
                                   psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize)
        ;;psf_base0=Convolve2(beam1_0,beam2_0,pad=convolve_pad,/absolute)
        psf_base0=dirty_image_generate(beam1_0*beam2_0)
        uv_mask=fltarr(psf_dim2,psf_dim2)
        threshold=Max(psf_base0)/1E6
        beam_i=region_grow(psf_base0,psf_dim2*(1.+psf_dim2)/2.,thresh=[Max(psf_base0)/1E6,Max(psf_base0)])
        uv_mask[beam_i]=1.
        psf_base0*=uv_mask

        fine_beam[*, *, freq_i, pol_i] = psf_base0
        
        beam_test1=mwa_tile_beam_generate(fltarr(16)+1.,antenna_beam_arr1, $ 
                                          frequency=freq_center[freq_i],polarization=pol1,za_arr=za_arr,az_arr=az_arr,$
                                          psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize)
        beam_test2=mwa_tile_beam_generate(fltarr(16)+1.,antenna_beam_arr2,$
                                      frequency=freq_center[freq_i],polarization=pol2,za_arr=za_arr,az_arr=az_arr,$
                                      psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize)
        psf_gain_test=dirty_image_generate(beam_test1*beam_test2)*uv_mask
        gain_normalization=1./(Total(psf_gain_test)/psf_resolution^2.)
        psf_base0*=gain_normalization

        FOR tile_i=0,n_tiles-1 DO BEGIN
           ;; *beam1_arr[tile_i]=tile_beam_generate(gain1[*,tile_i],antenna_beam_arr1,$
           ;;                                       frequency=freq_center[freq_i],polarization=pol1,$
           ;;                                       psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
           ;;                                       zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
           ;;                                       normalization=gain1_normalization)
           *beam1_arr[tile_i]=mwa_tile_beam_generate(gain1[*,tile_i],antenna_beam_arr1,$
                                                 frequency=freq_center[freq_i],polarization=pol1,za_arr=za_arr,az_arr=az_arr,$
                                                 psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize)
           ;; *beam2_arr[tile_i]=tile_beam_generate(gain2[*,tile_i],antenna_beam_arr2,$
           ;;                                       frequency=freq_center[freq_i],polarization=pol2,$
           ;;                                       psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize,$
           ;;                                       zenith_angle_offset=zenith_angle_offset,azimuth_angle_offset=azimuth_angle_offset,$
           ;;                                       normalization=gain2_normalization)
           *beam2_arr[tile_i]=mwa_tile_beam_generate(gain2[*,tile_i],antenna_beam_arr2,$
                                                 frequency=freq_center[freq_i],polarization=pol2,za_arr=za_arr,az_arr=az_arr,$
                                                 psf_dim=psf_dim,psf_resolution=psf_resolution,kbinsize=kbinsize)
        ENDFOR
        
        FOR bi=0L,nbaselines-1 DO BEGIN
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
           ;;psf_single=Convolve2(*beam1_arr[tile_A[bi]-1],*beam2_arr[tile_B[bi]-1],pad=convolve_pad,/absolute)
           psf_single=dirty_image_generate(*beam1_arr[tile_A[bi]-1],*beam2_arr[tile_B[bi]-1])*uv_mask
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
           *psf_base[pol_i,freq_i,i,j]=psf_base0[xvals_i+i,yvals_i+j] ;/Total(psf_base0[xvals_i+i,yvals_i+j]) 
     ENDFOR
  ENDFOR
  
  if n_elements(beam_shape_file) ne 0 then begin
     beam_uv = fine_beam
     save, file=beam_shape_file, beam_uv, psf_resolution, psf_dim
  endif

  save,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals,$
       filename=file_path+'_beams'+file_id+'.sav'

end
