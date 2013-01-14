
pro sim_3dps, simfile, psf = psf, refresh = refresh, no_kzero = no_kzero, beam_exp = beam_exp, std_power = std_power, $
              no_weighting = no_weighting, add_noise = add_noise, pub = pub, fill_holes = fill_holes, quiet = quiet, $
              eor_only = eor_only, test_power_shape = test_power_shape, eor_test = eor_test, clean_type = clean_type

  if n_elements(fill_holes) eq 0 then fill_holes = 0
  if n_elements(beam_exp) eq 0 then beam_exp = 0
  if keyword_set(eor_test) then eor_only = 1

  clean_type_enum = ['hmf', 'iterate', 'fit']
  if n_elements(clean_type) ne 0 then begin
     wh = where(clean_type_enum eq clean_type, count)
     if count eq 0 then message, 'Clean type not recognized'

     if beam_exp ne 0 then message, 'Cleaning is not compatible with primary beam removal'
     if keyword_set(eor_test) then message, 'Cleaning is not compatible with eor_test'
  endif

  if n_elements(simfile) eq 0 then begin 
     froot = base_path() + 'fhd/simulations/'
     simfile = froot + 'sim_496t_0025_uvf.idlsave'
  endif

  temp = strpos(simfile, '/', /reverse_search)
  froot = strmid(simfile, 0, temp+1)
  infilebase = strmid(simfile, temp+1)

  tests = [strmatch(infilebase, '*32t*', /fold_case), strmatch(infilebase, '*512t*', /fold_case), $
           strmatch(infilebase, '*496t*', /fold_case)]
  arrays = ['32t', '512t', '496t']
  if round(total(tests)) eq 1 then begin
     array = arrays[where(tests eq 1)]
     
     if n_elements(clean_type) ne 0 then begin
        if tests[0] eq 1 then t32=1 else begin
           if tests[1] eq 1 then use_outliers = 1
           if strmatch(infilebase, '*fine*', /fold_case) then fine_res = 1
        endelse
     endif
  endif else begin
     if n_elements(clean_type) ne 0 then define_baselines = 1
     test = strmatch(infilebase, '*simple*', /fold_case)
     if test eq 1 then begin
        pos1 = strpos(infilebase, 'simple')
        pos2 = strpos(infilebase, '_', pos1)
        array = strmid(infilebase, pos1, pos2-pos1)

        if n_elements(clean_type) ne 0 then begin
           baseline_layout = 'simple'
           baseline_spacing = strmid(array, 5)
        endif
     endif else begin
        test = strmatch(infilebase, '*single*', /fold_case)
        if test eq 1 then begin
           pos1 = strpos(infilebase, 'single')
           pos2 = strpos(infilebase, '_', pos1)
           array = strmid(infilebase, pos1, pos2-pos1)

           if n_elements(clean_type) ne 0 then begin
              baseline_layout = 'single'
              baseline_spacing = strmid(array, 5)
           endif
        endif else message, 'array not recognized'
     endelse
  endelse

  info_file = froot + 'sim_' + array + '_info.idlsave'
  weights_file = froot + 'sim_' + array + '_weights.idlsave'

  if keyword_set(psf) then begin
     filebase = 'psf_' + array + '_uvf'
  endif else begin
     if keyword_set(eor_only) then begin
        if keyword_set(eor_test) then eor_tag = 'test' else eor_tag = ''
        filebase = 'eor' + eor_tag + '_' + array + '_uvf'
     endif else if n_elements(test_power_shape) ne 0 then $
        filebase = test_power_shape + 'power_' + array + '_uvf' $
     else filebase = (strsplit(infilebase, '.idlsave', /regex, /extract))[0]
  endelse

  initial_uv_savefile = froot + filebase + '_initial.idlsave'

  if keyword_set(add_noise) then filebase = filebase + '_noise'
  if n_elements(clean_type) ne 0 then begin
     case clean_type of
        'hmf': filebase = filebase + '_clean'
        'iterate': filebase = filebase + '_cleanplus'
        'fit': filebase = filebase + '_cleanfit'
     endcase
  endif
  if not keyword_set(eor_test) then $
     if beam_exp eq 1 then filebase = filebase + '_beam' else if beam_exp eq 0 then filebase = filebase + '_holo'
  if keyword_set(std_power) then filebase = filebase + '_sp'
  save_file = froot + filebase + '_power.idlsave'

  if n_elements(clean_type) ne 0 then if clean_type ne 'fit' then model_save = froot + filebase + '_modeluv.idlsave'

  test_save = file_test(save_file) *  (1 - file_test(save_file, /zero_length))

  if test_save eq 0 or keyword_set(refresh) then begin

     restore, info_file
     restore, weights_file

     if keyword_set(add_noise) then begin
        noise_filebase = 'sim_' + array + '_noise_uvf'
        noise_file = froot + noise_filebase + '.idlsave'
        restore, noise_file
     endif
     
     if keyword_set(eor_only) then begin
        eor_filebase = 'sim_' + array + '_eor_uvf'
        if keyword_set(eor_test) then eor_filebase = eor_filebase + '_initial'
        eor_file = froot + eor_filebase + '.idlsave'
        restore, eor_file
        if n_elements(eor_cube) eq 0 then eor_cube = temporary(eor_uvf)
        dims = size(eor_cube, /dimensions)
     endif else if n_elements(test_power_shape) ne 0 then begin
        test_power_filebase = 'sim_' + array + '_' + test_power_shape + '_uvf'
      
        test_power_file = froot + test_power_filebase + '.idlsave'
        restore, test_power_file
        dims = size(test_power_cube, /dimensions)

     endif else begin
        restore, simfile
        dims = size(uvf_cube, /dimensions)
     endelse

     ;; check whether or not the frequencies are evenly spaced.
     n_freq = n_elements(frequencies)
     freq_diff = frequencies - shift(frequencies, 1)
     freq_diff = freq_diff[1:*]
     
     z0_freq = 1420.40 ;; MHz
     redshifts = z0_freq/frequencies - 1
     cosmology_measures, redshifts, comoving_dist_los = comov_dist_los
     
     comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
     comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
     
     if max(freq_diff-freq_diff[0]) gt 1e-12 then begin
        ;; frequencies are not evenly spaced, need to be careful about z_mpc_delta/mean
        
        nominal_freqs = dindgen(floor(((max(frequencies)-min(frequencies))/freq_resolution))+1)*freq_resolution + min(frequencies)
        nominal_z = z0_freq/nominal_freqs - 1
        comoving_distance_los, nominal_z, nominal_comov_dist_los
        nominal_comov_diffs = nominal_comov_dist_los - shift(nominal_comov_dist_los, -1)
        nominal_comov_diffs = nominal_comov_diffs[0:n_elements(nominal_comov_diffs)-2]

        z_mpc_delta = mean(nominal_comov_diffs)
        z_mpc_mean = mean(nominal_comov_dist_los)
        
     endif else begin
        
        z_mpc_delta = mean(comov_los_diff)
        z_mpc_mean = mean(comov_dist_los)
        n_kz = n_freq
        
     endelse
     
     z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
     kz_mpc_range =  (2d*!pi) / (z_mpc_delta)
     kz_mpc_delta = (2d*!pi) / z_mpc_length
     kz_mpc = dindgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2
     if n_elements(n_kz) ne 0 then begin
        if n_elements(kz_mpc) ne n_kz then stop
     endif else n_kz = n_elements(kz_mpc)
    
     print, 'z delta: ', z_mpc_delta
     print,  'kz delta: ', kz_mpc_delta

     x_rad_delta = abs(degpix) * !pi / 180d
     n_kx = dims[0]
     x_rad_length = dims[0] * x_rad_delta
     x_mpc_delta = x_rad_delta * z_mpc_mean
     x_mpc_length = x_rad_length * z_mpc_mean
     kx_mpc_range = (2d*!pi) / x_mpc_delta
     kx_mpc_delta = (2d*!pi) / x_mpc_length
     kx_mpc = (dindgen(n_kx)-n_kx/2) * kx_mpc_delta

     y_rad_delta = abs(degpix) * !pi / 180d
     n_ky = dims[1]
     y_rad_length = dims[1] * y_rad_delta
     y_mpc_delta = y_rad_delta * z_mpc_mean
     y_mpc_length = y_rad_length * z_mpc_mean
     ky_mpc_range = (2d*!pi) / y_mpc_delta
     ky_mpc_delta = (2d*!pi) / y_mpc_length
     ky_mpc = (dindgen(n_ky)-n_ky/2) * ky_mpc_delta
     
     ;; beam_diameter_rad = (3d * 10^8d) / (frequencies * 10^6d * max_baseline)
     ;; beam_area_str = !pi * beam_diameter_rad^2d /4d

     ;; conv_factor = (10^(double(-26+16+3-12+23)) * 9d) / (beam_area_str * 2d * frequencies^2d * 1.38)
     ;; if max(conv_factor-conv_factor[0]) gt 1e-8 then stop else conv_factor = conv_factor[0]
     conv_factor = 2d * max_baseline^2d / (!dpi * 1.38)

     kperp_lambda_conv = z_mpc_mean / (2d*!pi)

     if max(abs(imaginary(weights_cube))) eq 0 then weights_cube = real_part(weights_cube) $
     else stop

     if keyword_set(psf) then begin
        sim_cube = weights_cube * conv_factor
     endif else begin
        if keyword_set(eor_only) then sim_cube = temporary(eor_cube) $
        else if n_elements(test_power_shape) ne 0 then sim_cube = test_power_cube $
        else begin
           sim_cube = uvf_cube * conv_factor
           undefine, my_sources
           ;; save initial uv slice   
           uv_slice = uvf_slice(reform(my_model_uv, n_kx, n_ky, 1), kx_mpc, ky_mpc, [mean(frequencies)], kperp_lambda_conv, $
                                slice_axis = 2, slice_inds = 0, slice_savefile = initial_uv_savefile)
           undefine, my_model_uv
        endelse
     endelse
     undefine, uvf_cube

     if keyword_set(eor_test) then weights_cube = temporary(weights_cube)*0d + 1d
     sigma2_cube = 1d/(weights_cube)
     wh_wt0 = where(weights_cube eq 0, count_wt0)
     ;; wh = where(weights_cube le 1e-10, count)
     if count_wt0 ne 0 then sigma2_cube[wh_wt0] = 0
     undefine, weights_cube

     mask = dblarr(n_kx, n_ky, n_kz) + 1
     if count_wt0 gt 0 then mask[wh_wt0] = 0
     ;; n_pix_contrib = total(total(mask, 2), 1)
     n_freq_contrib = total(mask, 3)
     wh_nofreq = where(n_freq_contrib eq 0, count_nofreq)
     undefine, mask

      print, 'pre-weighting sum(sim_cube^2)*z_delta:', total(abs(sim_cube)^2d)*z_mpc_delta
;;if keyword_set(eor_only) then stop

     ;; divide by weights (~array beam) to estimate true sky
     sim_cube = sim_cube * sigma2_cube
     if count_wt0 ne 0 then sim_cube[wh_wt0] = 0
     
     print, 'sum(sim_cube^2)*z_delta (after weighting):', total(abs(sim_cube)^2d)*z_mpc_delta

     ;; for eor_test, skip beam corrections
     if not keyword_set(eor_test) then begin


        ;; Remove 1 or 2 factors of tile (primary) beam (in image space!) 
        ;; to go to detected frame or true sky frame
        if beam_exp gt 0 then begin
           for i = 0, n_freq - 1 do begin
 

              beam_base=fft_shift(real_part(FFT(fft_shift(tile_beam_uv[*,*,i]),/inverse)))
              beam_norm = max(beam_base)
              beam_base/=beam_norm
              beam_factor = (beam_base) ^ beam_exp
              wh_beam0 = where(beam_factor eq 0, count_beam0)
          
              temp = fft_shift(FFT(fft_shift(sim_cube[*,*,i]),/inverse))
             
              sim_image = temp / beam_factor
              if count_beam0 gt 0 then sim_image[wh_beam0] = 0

              temp2 = fft_shift(FFT(fft_shift(sim_image)))

              sim_cube[*,*,i] = temp2
              temp2=0
              temp=0
              
              ;;test normalization
              ;; temp = fft_shift(FFT(fft_shift(weights_cube[*,*,i]),/inverse))
              ;; weights_image = temp / (tile_beam_image[*,*,i]) ^ beam_exp
              
              ;; if keyword_set(add_noise) then begin
              ;;    ;; noise is defined in detected frame not holographic frame
              ;;    sky_noise_image = fft_shift(FFT(fft_shift(noise_cube[*,*,i]))) / (tile_beam_image) ^ (beam_exp - 1d)
              ;;    noise_cube[*,*,i] = fft_shift(FFT(fft_shift(sky_noise_image),/inverse))
              ;; endif
           endfor

           ;; save some slices of the sim cube (after beam correction)
           uf_savefile = froot + filebase + '_uf_plane_conv.idlsave'
           uf_slice = uvf_slice(sim_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 1, slice_inds = n_ky/2, $
                                slice_savefile = uf_savefile)
           
           vf_savefile = froot + filebase + '_vf_plane_conv.idlsave'
           vf_slice = uvf_slice(sim_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 0, slice_inds = n_kx/2, $
                                slice_savefile = vf_savefile)
           
           if max(abs(vf_slice)) eq 0 then begin
              nloop = 0
              while max(abs(vf_slice)) eq 0 do begin
                 nloop = nloop+1
                 vf_slice = uvf_slice(sim_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 0, $
                                      slice_inds = n_kx/2+nloop, slice_savefile = vf_savefile)
              endwhile
           endif

           uv_savefile = froot + filebase + '_uv_plane_conv.idlsave'
           uv_slice = uvf_slice(sim_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 2, slice_inds = 0, $
                                slice_savefile = uv_savefile)
           
        endif
     endif
     undefine, tile_beam_image

     if count_wt0 ne 0 then sim_cube[wh_wt0] = 0

     ;; if keyword_set(add_noise) then begin
     ;;    norm_max = max(abs(sim_cube))/max(abs(noise_cube))

     ;;   ;; hist = histogram(abs(sky_noise), binsize = 0.01, min = 0.01, locations = locs)
     ;;    ;; sim_hist = histogram(abs(sim_cube), binsize = 1e-8, min = 1e-8, locations = sim_locs)
     ;;    ;; norm_peak = sim_locs[(where(sim_hist eq max(sim_hist)))[0]] / locs[(where(hist eq max(hist)))[0]]

     ;;    sim_cube = sim_cube + temporary(noise_cube) * norm_max/10d
     ;; endif

     if n_elements(clean_type) ne 0 then begin

        if clean_type eq 'fit' then begin
           print, max(abs(sim_cube))

           ;; first subtract off k0
           k0_real = total(real_part(sim_cube), 3, /nan) / n_freq_contrib
           if count_nofreq gt 0 then k0_real[wh_nofreq]=0      
           k0_imag = total(imaginary(sim_cube), 3, /nan) / n_freq_contrib
           if count_nofreq gt 0 then k0_imag[wh_nofreq]=0

           pred_cube = rebin(k0_real, n_kx, n_ky, n_kz) + complex(0,1) * rebin(k0_imag, n_kx, n_ky, n_kz)
           if count_wt0 gt 0 then pred_cube[wh_wt0] = 0
           undefine, k0_real
           undefine, k0_imaginary

           new_cube = sim_cube - temporary(pred_cube)
           print, max(abs(new_cube))

           ;; then fit linear term and subtract it off
           ;; restore response to linear ramp
           x_savefile = froot + 'sim_' + array + '_' + 'xcom.idlsave'
           test_com_file = file_test(x_savefile)  *  (1 - file_test(x_savefile, /zero_length))
           if test_com_file eq 0 then sim_uv_com
           restore, x_savefile
           undefine, uv_x_in

           com_mean = total(uv_x_com, 3) / n_freq_contrib
           if count_nofreq gt 0 then com_mean[wh_nofreq] = 0
           x_com_norm = temporary(uv_x_com) - rebin(com_mean, n_kx, n_ky, n_kz)
           undefine, com_mean
           if count_wt0 gt 0 then x_com_norm[wh_wt0]=0

           ;; fit amplitude of linear term
           lin_amp_weights = total(x_com_norm^2d, 3)
           if count_nofreq gt 0 then lin_amp_weights[wh_nofreq] = 0
           real_lin_amp = total(real_part(new_cube) * x_com_norm,3) / lin_amp_weights
           imag_lin_amp = total(imaginary(new_cube) * x_com_norm,3) / lin_amp_weights
           if count_nofreq gt 0 then begin
              real_lin_amp[wh_nofreq] = 0
              imag_lin_amp[wh_nofreq] = 0
           endif
           pred_cube = x_com_norm * (rebin(real_lin_amp, n_kx, n_ky, n_kz) + complex(0,1) * rebin(imag_lin_amp, n_kx, n_ky, n_kz))
           if count_wt0 gt 0 then pred_cube[wh_wt0] = 0
           undefine, real_lin_amp
           undefine, imag_lin_amp

           new_cube = temporary(new_cube) - temporary(pred_cube)
           print, max(abs(new_cube))

           
           ;; now fit & subtract off 2nd derivative term
           ;; restore response to x^2 ramp
           x_savefile = froot + 'sim_' + array + '_' + 'xcurv.idlsave'
           test_com_file = file_test(x_savefile)  *  (1 - file_test(x_savefile, /zero_length))
           if test_com_file eq 0 then sim_uv_com
           restore, x_savefile
           undefine, uv_x_in
           
           com_curv_mean = total(uv_x_com, 3) / n_freq_contrib
           if count_nofreq gt 0 then com_curv_mean[wh_nofreq] = 0
           x_com_curv_norm = temporary(uv_x_com) - rebin(com_curv_mean, n_kx, n_ky, n_kz)
           undefine, com_curv_mean
           if count_wt0 gt 0 then x_com_curv_norm[wh_wt0]=0

           ;; subtract off linear term from x^2 response
           curv_lin_amp = total(x_com_curv_norm * x_com_norm,3) / lin_amp_weights
           if count_nofreq gt 0 then curv_lin_amp[wh_nofreq] = 0
           x_com_curv_norm = x_com_curv_norm - x_com_norm * rebin(curv_lin_amp, n_kx, n_ky, n_kz)
           undefine, curv_lin_amp

           ;; fit amplitude of curvature term
           curv_amp_weights = total(x_com_curv_norm^2d, 3)
           if count_nofreq gt 0 then curv_amp_weights[wh_nofreq] = 0
           real_curv_amp = total(real_part(new_cube) * x_com_curv_norm,3) / curv_amp_weights
           imag_curv_amp = total(imaginary(new_cube) * x_com_curv_norm,3) / curv_amp_weights
           if count_nofreq gt 0 then begin
              real_curv_amp[wh_nofreq] = 0
              imag_curv_amp[wh_nofreq] = 0
           endif
           pred_cube = x_com_curv_norm * (rebin(real_curv_amp, n_kx, n_ky, n_kz) + $
                                          complex(0,1) * rebin(imag_curv_amp, n_kx, n_ky, n_kz))
           if count_wt0 gt 0 then pred_cube[wh_wt0] = 0
           undefine, real_curv_amp
           undefine, imag_curv_amp

           new_cube = temporary(new_cube) - temporary(pred_cube)
           print, max(abs(new_cube))

           undefine, x_com_norm
           undefine, lin_amp_weights
           undefine, x_com_curv_norm
           undefine, curv_amp_weights

           ;; ;; restore response to linear ramp
           ;; x_savefile = froot + 'sim_' + array + '_' + 'xcom.idlsave'
           ;; test_com_file = file_test(x_savefile)  *  (1 - file_test(x_savefile, /zero_length))
           ;; if test_com_file eq 0 then sim_uv_com
           ;; restore, x_savefile
           ;; undefine, uv_x_in

           ;; com_mean = total(uv_x_com, 3) / n_freq_contrib
           ;; if count_nofreq gt 0 then com_mean[wh_nofreq] = 0
           ;; x_com_norm = temporary(uv_x_com) - rebin(com_mean, n_kx, n_ky, n_kz)
           ;; if count_wt0 gt 0 then x_com_norm[wh_wt0]=0
           
           ;; ;; figure out local phase slope in uv_x_com and compare with
           ;; ;; slope in k0 to determine amplitude of x_com_norm
           ;; meas_len = 3 ;; optimized, > 1 b/c of antenna size
           ;; com_slope = (shift(com_mean, [-1*meas_len,0]) - shift(com_mean, [meas_len,0])) / (2d*meas_len)
           ;; ;;undefine, com_mean

           ;; ;; restore response to x^2 ramp
           ;; x_savefile = froot + 'sim_' + array + '_' + 'xcurv.idlsave'
           ;; test_com_file = file_test(x_savefile)  *  (1 - file_test(x_savefile, /zero_length))
           ;; if test_com_file eq 0 then sim_uv_com
           ;; restore, x_savefile
           ;; undefine, uv_x_in
           
           ;; com_curv_mean = total(uv_x_com, 3) / n_freq_contrib
           ;; if count_nofreq gt 0 then com_curv_mean[wh_nofreq] = 0
           ;; x_com_curv_norm = temporary(uv_x_com) - rebin(com_curv_mean, n_kx, n_ky, n_kz)
           ;; if count_wt0 gt 0 then x_com_curv_norm[wh_wt0]=0

           ;; ;; subtract off linear term from x^2 response
           ;; com_curv_slope = (shift(com_curv_mean, [-1*meas_len,0]) - shift(com_curv_mean, [meas_len,0])) / (2d*meas_len)
           ;; com_curv = (shift(com_curv_mean, [-1*meas_len,0]) + shift(com_curv_mean, [meas_len,0]) - 2d * com_curv_mean) /meas_len^2d
           ;; undefine, com_curv_mean
           ;; lin_amp = temporary(com_curv_slope) / com_slope
           ;; if count_nofreq gt 0 then lin_amp[wh_nofreq] = 0
           ;; x_com_curv_norm = x_com_curv_norm - x_com_norm * rebin(lin_amp, n_kx,n_ky,n_kz)

           ;; undefine, lin_amp
           ;; if count_wt0 gt 0 then x_com_curv_norm[wh_wt0]=0

           ;; k0_real = total(real_part(sim_cube), 3, /nan) / n_freq_contrib
           ;; if count_nofreq gt 0 then k0_real[wh_nofreq]=0      
           ;; k0_real_slope = (shift(k0_real, [-1*meas_len,0]) - shift(k0_real, [meas_len,0])) / (2d*meas_len)
           ;; k0_real_curv = (shift(k0_real, [-1*meas_len,0]) + shift(k0_real, [meas_len,0]) + 2d * k0_real) / meas_len^2d

           ;; ;;real_lin_amp = temporary(k0_real_slope) / com_slope
           ;; real_lin_amp = k0_real_slope / com_slope
           ;; if count_nofreq gt 0 then real_lin_amp[wh_nofreq] = 0
           ;; ;;real_curv_amp = temporary(k0_real_curv) / com_curv
           ;; real_curv_amp = k0_real_curv / com_curv
           ;; if count_nofreq gt 0 then real_curv_amp[wh_nofreq] = 0

           ;; real_cube = rebin(k0_real, n_kx,n_ky,n_kz) + x_com_norm * rebin(real_lin_amp, n_kx,n_ky,n_kz) ;;+ $
           ;;             ;;x_com_curv_norm * rebin(real_curv_amp, n_kx,n_ky,n_kz)
           ;; if count_wt0 gt 0 then real_cube[wh_wt0]=0

           ;; undefine, real_lin_amp
           ;; undefine, real_curve_amp

           ;; k0_imag = total(imaginary(sim_cube), 3, /nan) / n_freq_contrib
           ;; if count_nofreq gt 0 then k0_imag[wh_nofreq]=0
        
           ;; k0_imag_slope = (shift(k0_imag, [-1*meas_len,0]) - shift(k0_imag, [meas_len,0])) / (2d*meas_len)
           ;; k0_imag_curv = (shift(k0_imag, [-1*meas_len,0]) + shift(k0_imag, [meas_len,0]) + 2d * k0_imag) / meas_len^2d

           ;; imag_lin_amp = k0_imag_slope / com_slope
           ;; if count_nofreq gt 0 then imag_lin_amp[wh_nofreq] = 0
           ;; imag_curv_amp = k0_imag_curv / com_curv
           ;; if count_nofreq gt 0 then imag_curv_amp[wh_nofreq] = 0

           ;; imag_cube = rebin(k0_imag, n_kx,n_ky,n_kz) + x_com_norm * rebin(imag_lin_amp, n_kx,n_ky,n_kz) ;;+ $
           ;;             ;;x_com_curv_norm * rebin(imag_curv_amp, n_kx,n_ky,n_kz)
           ;; if count_wt0 gt 0 then imag_cube[wh_wt0]=0

           ;; pred_cube = real_cube + complex(0,1) * imag_cube
       
           ;; new_cube = sim_cube - temporary(pred_cube)
           ;; new_cube = sim_cube * exp(dcomplex(0, -1d) * pred_phase)
           
        endif else begin
           ;; remove the contribution to the power from the edges of the antennas
           ;; by calculating the response (xfer fn) based on the k0 mode
           
           ;; sorted_mask = sort_nd(temporary(mask) * dindgen(n_kx, n_ky, n_kz), 3)
           
           ;; ;; unwrap phases: blantantly cribbed from JD Smith's phunwrap
           ;; ;; have to sort by mask so that unmeasured values between real values don't screw up wrapping detection
           ;; phases = atan(sim_cube[sorted_mask],/phase)
           ;; delta_phase = phases - shift(phases, [0,0,1])
           ;; delta_phase[*,*,0] = 0
           ;; p = 2d*!dpi * (fix((delta_phase GT !dpi) EQ 1) - fix((delta_phase LT (-1)*!dpi) EQ 1))
           ;; undefine, delta_phase
           ;; r = total(temporary(p), /cumulative, 3)
           ;; new_phase = temporary(phases) - temporary(r)
           
           ;; ;; undo sorting to get back to actual ordering
           ;; unwrapped_phases = new_phase[sort_nd(sorted_mask, 3)]
           ;; undefine, new_phase
           ;; if count_wt0 gt 0 then unwrapped_phases[wh_wt0]=0
           
           ;; k0_phase_est = total(unwrapped_phases, 3, /nan) / n_freq_contrib ;; this is the k0 mode, faster than taking fft
           ;; undefine, unwrapped_phases
           ;; if count_nofreq gt 0 then k0_phase_est[wh_nofreq]=0
           
           k0_est = total(sim_cube, 3, /nan) / n_freq_contrib
           if count_nofreq gt 0 then k0_est[wh_nofreq]=0      


           nloop = 0
           sig = 0
           if clean_type eq 'iterate' then loop_max = 5 else loop_max = 1
           while sig eq 0 do begin
              
              if nloop eq 0 then begin
                 ;; model_uv = exp(complex(0,1) * k0_phase_est)
                 model_uv = k0_est / conv_factor

                 print, max(abs(model_uv) * conv_factor)
                 model_elems = model_uv
              endif else begin
                 print, max(abs(uv_est) * conv_factor)
                 model_uv = model_uv + uv_est
                 model_elems = [[[model_elems]],[[uv_est]]]
              endelse
              
              model_response = sim_xfer_apply(model_uv, t32 = t32, define_baselines = define_baselines, $
                                              baseline_spacing = baseline_spacing, baseline_layout = baseline_layout, $
                                              fine_res = fine_res, use_outliers = use_outliers)
              
              ;; divide by weights & apply factors
              model_response = model_response * conv_factor * sigma2_cube
              if count_wt0 ne 0 then model_response[wh_wt0] = 0
              
              ;; subtract off the contribution
              new_cube = sim_cube - model_response
              
              nloop = nloop+1
              if nloop eq loop_max then sig = 1 else begin

                 ;; ;; unwrap phases: blantantly cribbed from JD Smith's phunwrap
                 ;; ;; have to sort by mask so that unmeasured values between real values don't screw up wrapping detection
                 ;; phases = atan(new_cube[sorted_mask],/phase)
                 ;; delta_phase = phases - shift(phases, [0,0,1])
                 ;; delta_phase[*,*,0] = 0
                 ;; p = 2d*!dpi * (fix((delta_phase GT !dpi) EQ 1) - fix((delta_phase LT (-1)*!dpi) EQ 1))
                 ;; undefine, delta_phase
                 ;; r = total(temporary(p), /cumulative, 3)
                 ;; new_phase = temporary(phases) - temporary(r)
                 
                 ;; ;; undo sorting to get back to actual ordering
                 ;; unwrapped_phases = new_phase[sort_nd(sorted_mask, 3)]
                 ;; undefine, new_phase
                 ;; if count_wt0 gt 0 then unwrapped_phases[wh_wt0]=0

                 ;; res_k0_phase = total(temporary(unwrapped_phases), 3, /nan) / n_freq_contrib
                 ;; if count_nofreq gt 0 then res_k0_phase[wh_nofreq]=0

                 ;; uv_est = (abs(new_cube)/(conv_factor*fudge_factor)) * exp(complex(0,1) * res_k0_phase)
                 ;; if count_nofreq gt 0 then uv_est[wh_nofreq] = 0

                 uv_est = total(new_cube, 3, /nan) / n_freq_contrib
                 if count_nofreq gt 0 then uv_est[wh_nofreq]=0      

                 uv_est = uv_est / conv_factor
              endelse
           endwhile
           save, file = model_save, model_elems
           
           uv_est = 0
           model_uv = 0
           model_elems = 0
        endelse

        print, minmax(abs(sim_cube))
        print, minmax(abs(new_cube))
        sim_cube = temporary(new_cube)
     endif else undefine, mask

     ;; save some slices of the sim cube (post cleaning)
     uf_savefile = froot + filebase + '_uf_plane.idlsave'
     uf_slice = uvf_slice(sim_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 1, slice_inds = n_ky/2, $
                slice_savefile = uf_savefile)

     vf_savefile = froot + filebase + '_vf_plane.idlsave'
     vf_slice = uvf_slice(sim_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 0, slice_inds = n_kx/2, $
                slice_savefile = vf_savefile)

     if max(abs(vf_slice)) eq 0 then begin
        nloop = 0
        while max(abs(vf_slice)) eq 0 do begin
           nloop = nloop+1
           vf_slice = uvf_slice(sim_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 0, slice_inds = n_kx/2+nloop, $
                                slice_savefile = vf_savefile)
        endwhile
     endif

     uv_savefile = froot + filebase + '_uv_plane.idlsave'
     uv_slice = uvf_slice(sim_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 2, slice_inds = 0, $
                slice_savefile = uv_savefile)

     ;; need to cut uvf cubes in half because image is real -- we'll cut in v
     sim_cube = sim_cube[*, n_ky/2:n_ky-1,*]
     sigma2_cube = sigma2_cube[*, n_ky/2:n_ky-1,*]
     n_freq_contrib = n_freq_contrib[*, n_ky/2:n_ky-1]

     ky_mpc = ky_mpc[n_ky/2:n_ky-1]
     n_ky = n_elements(ky_mpc)

     print, 'sum(sim_cube^2)*z_delta (after cut):', total(abs(sim_cube)^2d)*z_mpc_delta

     ;; now take FFT
     sim_ft = fft(sim_cube, dimension=3) * n_freq * z_mpc_delta / (2d*!dpi)
     ;; put k0 in middle of cube
     sim_ft = shift(sim_ft, [0,0,n_kz/2])
     ;;undefine, sim_cube

     print, 'full_ft^2d integral:', total(abs(sim_ft)^2d)
     print, 'full_ft^2d integral * 2pi*delta_k^2d:', total(abs(sim_ft)^2d) * kz_mpc_delta * 2d * !dpi
 
     ;; factor to go to eor theory FT convention
     ;; Only 1 factor of 2pi b/ only transforming along z
     ;; factor = (2d*!pi)
     ;; if not keyword_set(eor_test) then sim_ft = factor * temporary(sim_ft)

     ;; print, 'full_ft^2d integral (after theory factor):', total(abs(sim_ft)^2d)

     n_val = kz_mpc / kz_mpc_delta
     a_0 = 2d * sim_ft[*,*,where(n_val eq 0)]
     a_n = sim_ft[*,*, where(n_val gt 0)] + sim_ft[*,*, reverse(where(n_val lt 0))]
     b_n = complex(0,1) * (sim_ft[*,*, where(n_val gt 0)] - sim_ft[*,*, reverse(where(n_val lt 0))])
     undefine, sim_ft
   
     kz_mpc = kz_mpc[where(n_val ge 0)]
     n_kz = n_elements(kz_mpc)

     if keyword_set(std_power) then begin
        ;; for standard power calc. just need ft of sigma2
        sigma2_ft = fft(sigma2_cube, dimension=3) * n_freq * z_mpc_delta / (2d*!dpi)
        sigma2_ft = shift(sigma2_ft, [0,0,n_kz/2])

        undefine, sigma2_cube

        ;; factor to go to eor theory FT convention
        ;;sigma2_ft = factor * temporary(sigma2_ft)

        power_3d = dblarr(n_kx, n_ky, n_kz)
        power_3d[*,*,0] = (a_0 * conj(a_0))/4d
        power_3d[*,*,1:n_kz-1] = ((a_n * conj(a_n)) + (b_n * conj(b_n)))/2d
    
        ;; power_3d[*,*,0] = real_part(sim_ft[*,*,where(n_val eq 0)] * conj(sim_ft[*,*,where(n_val eq 0)]))
        ;; power_3d[*,*,where(n_val gt 0)] = real_part(sim_ft[*,*,where(n_val gt 0)] * conj(sim_ft[*,*,where(n_val gt 0)])) + $
        ;;                   reverse(real_part(sim_ft[*,*,where(n_val lt 0)] * conj(sim_ft[*,*,where(n_val lt 0)])), 3)
        ;;undefine, sim_ft

        sigma_a0 = 2d * abs(sigma2_ft[*,*,where(n_val eq 0)])
        sigma_an_bn = sqrt(abs(sigma2_ft[*,*, where(n_val gt 0)])^2d + abs(sigma2_ft[*,*, reverse(where(n_val lt 0))])^2d)

        sigma2_3d = dblarr(n_kx, n_ky, n_kz)
        sigma2_3d[*,*,0] = abs(a_0) * sigma_a0 / 2d
        sigma2_3d[*,*,1:n_kz-1] = sqrt(abs(a_n)^2d + abs(b_n)^2d) * sigma_an_bn
        undefine, sigma2_ft

        weights_3d = 1d/sigma2_3d
        wh_sig0 = where(sigma2_3d eq 0, count_sig0)
        if count_sig0 gt 0 then weights_3d[wh_sig0] = 0
        sigma2_3d=0

     endif else begin  
 
        ;; force pixels with only 1 measurement to only have k0
        wh_1freq = where(n_freq_contrib eq 1, count_1freq)
        if count_1freq gt 0 then begin
           mask = n_freq_contrib * 0 + 1
           mask[wh_1freq] = 0
           mask = rebin(temporary(mask), n_kx, n_ky, n_kz-1)

           a_n = temporary(a_n) * mask
           b_n = temporary(b_n) * mask
        endif

        sim_cos = complex(dblarr(n_kx, n_ky, n_kz))
        sim_sin = complex(dblarr(n_kx, n_ky, n_kz))
        sim_cos[*, *, 0] = a_0 /2d
        sim_cos[*, *, 1:n_kz-1] = a_n
        sim_sin[*, *, 1:n_kz-1] = b_n
 
        ;; for new power calc, need cos2, sin2, cos*sin transforms
        ;; have to do this in a for loop for memory's sake
        covar_cos = dblarr(n_kx, n_ky, n_kz)
        covar_sin = dblarr(n_kx, n_ky, n_kz)
        covar_cross = dblarr(n_kx, n_ky, n_kz)
        
        ;; comov_dist_los goes from large to small z
        z_relative = dindgen(n_freq)*z_mpc_delta
        freq_kz_arr = rebin(reform(rebin(reform(kz_mpc, 1, n_kz), n_freq, n_kz) * $
                                   rebin(z_relative, n_freq, n_kz), 1, n_freq, n_kz), n_ky, n_freq, n_kz)
        
        cos_arr = cos(freq_kz_arr)
        sin_arr = sin(freq_kz_arr)
        ;;z_exp_arr = exp(-1*dcomplex(0,1)*freq_kz_arr)
        
        for i=0, n_kx-1 do begin
           if max(n_freq_contrib[i,*]) eq 0 then continue
           wh_divide = where(n_freq_contrib[i,*] gt 0, count_divide)
           n_divide = double(rebin(reform(n_freq_contrib[i,wh_divide]),count_divide, n_kz))

           sigma2_arr = rebin(reform(sigma2_cube[i,wh_divide,*]), count_divide, n_freq, n_kz) 

           covar_cos[i,wh_divide,*] = total(sigma2_arr*cos_arr[wh_divide, *, *]^2d, 2)/n_divide
           covar_sin[i,wh_divide,*] = total(sigma2_arr*sin_arr[wh_divide, *, *]^2d, 2)/n_divide
           covar_cross[i,wh_divide,*] = total(sigma2_arr*cos_arr[wh_divide, *, *]*sin_arr[wh_divide, *, *], 2)/n_divide
        endfor
        
        ;; force pixels with only 1 measurement to only have k0
        wh_1freq = where(n_freq_contrib eq 1, count_1freq)
        if count_1freq gt 0 then begin
           mask = [[[dblarr(n_kx, n_ky) + 1]], [[temporary(mask)]]]
           covar_cos = temporary(covar_cos) * mask
           covar_sin = temporary(covar_sin) * mask
           covar_cross = temporary(covar_cross) * mask
           undefine, mask
        endif
        
        ;; cos 0 term has different normalization
        covar_cos[*,*,0] = covar_cos[*,*,0]/4d
       
        undefine, sigma2_cube
        undefine, freq_kz_arr
        undefine, cos_arr
        undefine, sin_arr
        undefine, sigma2_arr

        ;; factor to go to eor theory FT convention
        ;; I don't think I need these factors in the covariance
        ;; matrix because I've use the FT & inv FT -- should cancel
        ;; covar_cos = factor * temporary(covar_cos2)
        ;; covar_sin = factor * temporary(covar_sin2)
        ;; covar_cross = factor * temporary(covar_cross)

        ;; get rotation angle to diagonalize covariance block
        theta = atan(2*covar_cross, covar_cos - covar_sin)/2d
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        undefine, theta

        sigma1_2 = covar_cos*cos_theta^2 + 2d*covar_cross*cos_theta*sin_theta + covar_sin*sin_theta^2d
        sigma2_2 = covar_cos*sin_theta^2d - 2d*covar_cross*cos_theta*sin_theta + covar_sin*cos_theta^2d

        undefine, covar_cos
        undefine, covar_sin
        undefine, covar_cross
       
        sim1 = sim_cos*cos_theta + sim_sin*sin_theta
        sim2 = (-1d)*sim_cos*sin_theta + sim_sin*cos_theta
        undefine, sim_cos
        undefine, sim_sin

        weights_1 = (1/sigma1_2)^2d
        term1 = real_part(sim1 * conj(sim1))*weights_1
        wh_sig0 = where(sigma1_2^2d eq 0, count_sig0)
        if count_sig0 ne 0 then begin
           weights_1[wh_sig0] = 0
           term1[wh_sig0] = 0
        endif
        
        weights_2 = (1/sigma2_2)^2d
        term2 = real_part(sim2 * conj(sim2))*weights_2
        wh_sig0 = where(sigma2_2^2d eq 0, count_sig0)
        if count_sig0 ne 0 then begin
           weights_2[wh_sig0] = 0
           term2[wh_sig0] = 0
        endif
        undefine, sim1
        undefine, sim2
        undefine, sigma1_2
        undefine, sigma2_2   

        weights_3d = weights_1 + weights_2
        power_3d = (term1 + term2) / weights_3d
        ;;power_error = 1/weights
        wh_wt0 = where(weights_3d eq 0, count_wt0)
        if count_wt0 ne 0 then begin
           power_3d[wh_wt0] = 0
        endif

        undefine, term1
        undefine, term2
        undefine, weights1
        undefine, weights2

     endelse

     save, file = save_file, power_3d, weights_3d, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv, n_freq_contrib

  endif else restore, save_file

  print, 'power integral:', total(power_3d)

  n_kx = n_elements(kx_mpc)
  n_ky = n_elements(ky_mpc)
  n_kz = n_elements(kz_mpc)

  if keyword_set (no_kzero) then begin 
     ;; leave out kz=0 -- full of foregrounds
     kz_mpc = kz_mpc[1:*]
     power_3d = temporary(power_3d[*, *, 1:*])
     weights_3d = temporary(weights_3d[*,*,1:*])
     n_kz = n_elements(kz_mpc)
  endif

  fadd = ''
  if keyword_set(no_weighting) then fadd = fadd + '_nowt'
  if keyword_set(no_kzero) then fadd = fadd + '_nok0'
  if keyword_set(fill_holes) then fadd = fadd + '_nohole'

  savefile = froot + filebase + fadd + '_2dkpower.idlsave'
  savefile_lin = froot + filebase + fadd + '_linkpar_2dkpower.idlsave'

  print, 'Binning to 2D power spectrum'

  bins_per_decade = 8
  if keyword_set(no_weighting) then $
     power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, $
                                       binned_weights = binned_weights, bins = bins_per_decade, fill_holes = fill_holes) $
  else power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, $
                                         weights = weights_3d, binned_weights = binned_weights, bins = bins_per_decade, $
                                         fill_holes = fill_holes)
  power = power_rebin
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  weights = binned_weights
  
  save, file = savefile, power, weights, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv

  ;; also bin with linear k_par

  if keyword_set(no_weighting) then $
     power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, /linear_kpar, $
                                       binned_weights = binned_weights, bins = bins_per_decade, fill_holes = fill_holes) $
  else power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, /linear_kpar, $
                                         weights = weights_3d, binned_weights = binned_weights, bins = bins_per_decade, $
                                         fill_holes = fill_holes)
  power = power_rebin
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  weights = binned_weights
  
  save, file = savefile_lin, power, weights, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv

 
  plotfile_path = base_path() + 'power_spectrum/plots/'
  plotfile = plotfile_path + filebase + fadd + '_kspace_power'
  weight_plotfile = plotfile_path + filebase + fadd + '_kspace_weights'

  kperp_plot_range = [min(kperp_edges), 0.3]
  if not keyword_set(quiet) then begin
     kpower_2d_plots, savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = data_range, pub = pub, plotfile = plotfile
     kpower_2d_plots, savefile, /plot_weights, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      pub = pub, plotfile = weight_plotfile, window_num = 2, title = 'Weights'
  endif

  ;; now do slices    
  yslice_savefile = froot + filebase + '_xz_plane.idlsave'
  yslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, weights_3d, kperp_lambda_conv, slice_axis = 1, slice_inds = 0, $
                       slice_weights = yslice_weights, slice_savefile = yslice_savefile)

  xslice_savefile = froot + filebase + '_yz_plane.idlsave'
  xslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, weights_3d, kperp_lambda_conv, slice_axis = 0, slice_inds = n_kx/2, $
                              slice_weights = xslice_weights, slice_savefile = xslice_savefile)
  if max(xslice_power) eq 0 then begin
     nloop = 0
     while max(xslice_power) eq 0 do begin
         nloop = nloop+1
         xslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, weights_3d, kperp_lambda_conv, slice_axis = 0, $
                                     slice_inds = n_kx/2+nloop, slice_weights = xslice_weights, slice_savefile = xslice_savefile)
     endwhile
  endif

  zslice_savefile = froot + filebase + '_xy_plane.idlsave'
  zslice_power = kpower_slice(power_3d, kx_mpc, ky_mpc, kz_mpc, weights_3d, kperp_lambda_conv, slice_axis = 2, slice_inds = 1, $
                       slice_weights = zslice_weights, slice_savefile = zslice_savefile)

  ;; also make binned versions of x & y slices
  yslice_binned_savefile = froot + filebase + fadd + '_xz_plane_binned.idlsave'
  if keyword_set(no_weighting) then $
    yslice_power_rebin = kspace_rebinning_2d(yslice_power, kx_mpc, [ky_mpc[0]], kz_mpc, kperp_edges_mpc, kpar_edges_mpc, $
                                             binned_weights = binned_weights, bins = bins_per_decade, fill_holes = fill_holes) $
  else yslice_power_rebin = kspace_rebinning_2d(yslice_power, kx_mpc, [ky_mpc[0]], kz_mpc, kperp_edges_mpc, kpar_edges_mpc, $
                                               weights = yslice_weights, binned_weights = binned_weights, bins = bins_per_decade, $
                                               fill_holes = fill_holes)
     
  power = yslice_power_rebin
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  weights = binned_weights
  
  save, file = yslice_binned_savefile, power, weights, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv
  
  if not keyword_set(quiet) then begin
     kpower_2d_plots, yslice_binned_savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = data_range, pub = pub, window_num = 3, title = 'XZ plane'
  endif


  xslice_binned_savefile = froot + filebase + fadd + '_yz_plane_binned.idlsave'
  if keyword_set(no_weighting) then $
    xslice_power_rebin = kspace_rebinning_2d(xslice_power, [kx_mpc[n_kx/2]], ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, $
                                             binned_weights = binned_weights, bins = bins_per_decade, fill_holes = fill_holes) $
  else xslice_power_rebin = kspace_rebinning_2d(xslice_power, [kx_mpc[n_kx/2]], ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, $
                                               weights = xslice_weights, binned_weights = binned_weights, bins = bins_per_decade, $
                                               fill_holes = fill_holes)
     
  power = xslice_power_rebin
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  weights = binned_weights
  
  save, file = xslice_binned_savefile, power, weights, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv
     
  if not keyword_set(quiet) then begin
     kpower_2d_plots, xslice_binned_savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = data_range, pub = pub, window_num = 4, title = 'YZ plane'
  endif




  print, 'Binning to 1D power spectrum'
  plotfile = plotfile_path + filebase + fadd + '_1d_power'

  bins_per_decade = 10d
  if keyword_set(no_weighting) then $
     power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, bins = bins_per_decade, $
                                    binned_weights = weights_1d, mask = mask, pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, $
                                    k2_mask = k2_mask,  k3_mask = k3_mask, edge_on_grid = edge_on_grid, match_datta = match_datta)$
  else power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, bins = bins_per_decade, weights = weights_3d, $
                                      binned_weights = weights_1d, mask = mask, pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, $
                                      k2_mask = k2_mask,  k3_mask = k3_mask, edge_on_grid = edge_on_grid, match_datta = match_datta)

  power = power_1d
  weights = weights_1d
  k_edges = k_edges_mpc

  savefile = froot + filebase + fadd + '_1dkpower.idlsave'
  save, file = savefile, power, weights, k_edges, bins_per_decade

  eor_file_1d = base_path() + 'power_spectrum/eor_data/eor_power_1d.idlsave'
  file_arr = [savefile, eor_file_1d]
  if keyword_set(eor_only) then begin
     if keyword_set(eor_test) then names_arr = 'Input EoR' else names_arr = 'Simulated EoR'
  endif else names_arr = 'Simulation PS'
  names_arr = [names_arr, 'EoR signal']
  colors_arr = [0, 254]

    if not keyword_set(quiet) then begin
       kpower_1d_plots, file_arr, window_num = 5, names = names_arr, colors = colors_arr
    endif
end
