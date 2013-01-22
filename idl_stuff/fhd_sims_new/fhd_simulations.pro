
pro fhd_simulations, nfreq_per_image = nfreq_per_image, x_offset = x_offset, y_offset = y_offset, flux = flux, $
                     sim_file = sim_file, info_file = info_file, weights_file = weights_file, no_sim = no_sim, noise_gen = noise_gen, $
                     noise_file = noise_file, eor_gen = eor_gen, eor_file = eor_file, test_power_shape = test_power_shape, $
                     test_power_file = test_power_file, no_save = no_save, t32 = t32, define_baselines = define_baselines, $
                     baseline_spacing = baseline_spacing, baseline_layout = baseline_layout, fine_res = fine_res, $
                     use_outliers = use_outliers, refresh = refresh
                       

  if keyword_set(define_baselines) then begin
     if n_elements(baseline_spacing) eq 0 then baseline_spacing = 4
     if n_elements(baseline_layout) eq 0 then baseline_layout = 'simple'

     spacing_name = number_formatter(string(baseline_spacing))
     filename = 'Simulation_' + baseline_layout + spacing_name

  endif else if keyword_set(t32) then filename = 'Simulation_32t' else begin
     filename = 'Simulation'
     if keyword_set(use_outliers) then filename = filename + '_512t' else filename = filename + '_496t'
     if keyword_set(fine_res) then filename = filename + '_fine'
  endelse
       
  file_path = base_path('data') + 'fhd_simulations/fhd_files/'
  sim_path = base_path('data') + 'fhd_simulations/'
  pol_names = ['xx','yy','xy','yx']

  if keyword_set(simple_baselines) then tag =  baseline_layout + spacing_name + '_' $
  else if keyword_set(t32) then tag = '32t_' else begin
     if keyword_set(use_outliers) then tag = '512t_' else tag = '496t_'
     if keyword_set(fine_res) then tag = tag + 'fine_'
  endelse
  if n_elements(sim_file) ne 1 then sim_file = sim_path + 'sim_' + tag + 'uvf.idlsave'
  if n_elements(info_file) ne 1 then info_file = sim_path + 'sim_' + tag + 'info.idlsave'
  if n_elements(weights_file) ne 1 then weights_file = sim_path + 'sim_' + tag + 'weights.idlsave'
  if n_elements(beam_shape_file) ne 1 then beam_shape_file = sim_path + 'sim_' + tag + 'beam.idlsave'
  if keyword_set(noise_gen) and n_elements(noise_file) ne 1 then noise_file = sim_path + 'sim_' + tag + 'noise_uvf.idlsave'

  if keyword_set(refresh) then restore_last=0 else $
     restore_last = 1

  if n_elements(nfreq_per_image) eq 0 then nfreq_per_image = 24
  sim_beam_setup,data_directory=data_directory,filename=filename, t32 = t32, define_baselines = define_baselines, $
                 baseline_spacing = baseline_spacing, baseline_layout = baseline_layout, fine_res = fine_res, $
                 use_outliers = use_outliers, beam_shape_file = beam_shape_file, nfreq_per_image = nfreq_per_image, $
                 restore_last = restore_last
;;stop
  ;;beam_setup,data_directory=data_directory,filename=filename,/restore_last

  c_light=299792458.
  max_baseline = max(sqrt((uu_arr*c_light)^2d + (vv_arr*c_light)^2d))

  n_samples=N_Elements(bin_offset)
  nbaselines=n_elements(baseline_arr)/n_samples

  ntiles = floor(sqrt(2d*nbaselines))
  nbaselines_per_xferfn = long(ntiles) * (long(ntiles) - 1L) * long(nfreq_per_image)

  delta_freq = frequency_array[1] - frequency_array[0]
  nfreq = n_elements(frequency_array) / nfreq_per_image
  freq_resolution = delta_freq * nfreq_per_image / 1e6 ;; in MHz

  norm_source = reform([0, 0, 1], 3, 1) ;; columns are: x location (zero centered), y location (zero centered), flux
  n_sources=1
  norm_source_array = reform(fltarr(7, n_sources), 7, n_sources)
  norm_source_array[0,*]=norm_source[0,*] + dimension/2.
  norm_source_array[1,*]=norm_source[1,*] + elements/2.
  norm_source_array[4,*]=norm_source[2,*]
  norm_source_array[5,*]=Round(norm_source_array[0,*])+Round(norm_source_array[1,*])*dimension

  n_sources = max([n_elements(x_offset), n_elements(y_offset), n_elements(flux)]) 
  if n_sources eq 1 then begin
     if n_elements(x_offset) eq 0 then x_offset = 25.3
     if n_elements(y_offset) eq 0 then y_offset = 0
     if n_elements(flux) eq 0 then flux = 0.5
  endif else begin
     if n_elements(x_offset) ne n_sources or n_elements(y_offset) ne n_sources or n_elements(flux) ne n_sources then $
        message, 'number of elements in x_offset, y_offset, flux must be equal'
  endelse
  ;; columns are: x location (zero centered), y location (zero centered), flux
  my_sources = fltarr(3, n_sources)
  my_sources[0,*] = x_offset
  my_sources[1,*] = y_offset
  my_sources[2,*] = flux

  my_source_array = fltarr(7, n_sources)
  my_source_array[0,*]=my_sources[0,*] + dimension/2.
  my_source_array[1,*]=my_sources[1,*] + elements/2.
  my_source_array[4,*]=my_sources[2,*]
  my_source_array[5,*]=Round(my_source_array[0,*])+Round(my_source_array[1,*])*dimension

  if keyword_set(eor_gen) then begin
     if n_elements(eor_file) ne 1 then eor_file = sim_path + 'sim_' + tag + 'eor_uvf.idlsave'
     kx_arr = (indgen(dimension) - dimension/2)*kbinsize
     ky_arr = (indgen(elements) - elements/2)*kbinsize
     freq_arr = indgen(nfreq) * freq_resolution + mean(frequency_array[indgen(nfreq_per_image)] / 1e6) ;; in MHz

     eor_uvf = eor_sim(kx_arr, ky_arr, freq_arr)
     save, file = sim_path + 'sim_' + tag + 'eor_uvf_initial.idlsave', eor_uvf, kx_arr, ky_arr, freq_arr
  endif

  if n_elements(test_power_shape) gt 1 then message, 'Only one test power shape can be run at a time'
  if n_elements(test_power_shape) ne 0 then begin
     case test_power_shape of
        'flat': begin
           if n_elements(test_power_file) ne 1 then $
              test_power_file = sim_path + 'sim_' + tag + test_power_shape + '_uvf.idlsave'
           sigma = 100
    
           kx_arr = (indgen(dimension) - dimension/2)*kbinsize
           ky_arr = (indgen(elements) - elements/2)*kbinsize
           freq_arr = indgen(nfreq) * freq_resolution + mean(frequency_array[indgen(nfreq_per_image)] / 1e6) ;; in MHz

           test_power_uvf = eor_sim(kx_arr, ky_arr, freq_arr, flat_sigma = sigma)       
           save, file = sim_path + 'sim_' + tag + test_power_shape +'_uvf_initial.idlsave', $
                 test_power_uvf, kx_arr, ky_arr, freq_arr
        end
        else: message, 'Test power shape not recognized'
     endcase
  endif

  weights_cube = make_array(dimension, elements, nfreq, /float)
  tile_beam_image = make_array(dimension, elements, nfreq, /float)
  if not keyword_set(no_sim) then uvf_cube = make_array(dimension, elements, nfreq, /complex)
  if keyword_set(noise_gen) then noise_cube = make_array(dimension, elements, nfreq, /complex)
  if keyword_set(eor_gen) then eor_cube = make_array(dimension, elements, nfreq, /complex)
  if n_elements(test_power_shape) ne 0 then test_power_cube = make_array(dimension, elements, nfreq, /complex)
  frequencies = dblarr(nfreq)
  normalizations = dblarr(nfreq)
 

  ;; need to make image to get normalization (to get to Jy)
  ;; for now force 1 Jy to give all 1's on the uv plane
  ;; (should maybe have factors of 2 pi...)
  test_gridding = visibility_source_uv_grid(norm_source_array)
  norm_gridding = 1/max(abs(test_gridding))

  if not keyword_set(no_sim) then begin
     ;; make model uv sky (same for all freq)
     my_model_uv=visibility_source_uv_grid(my_source_array, timing=timing) * norm_gridding
     print, 'source gridding (s): ',timing
  endif
stop
  for i=0, nfreq-1 do begin
     inds_to_use = indgen(nfreq_per_image) + nfreq_per_image * i
     print, 'frequency_bins: ' + string(min(inds_to_use), format = '(i3.3)') + ' to ' + string(max(inds_to_use), format = '(i3.3)')

     file_identifier=String(format='(A2,"_f",A3,"_f",A3)',pol_names[0], string(min(inds_to_use), format = '(i3.3)'), $
                            string(max(inds_to_use), format = '(i3.3)'))
     xfer_file = file_path+'_xferfn_'+file_identifier+'.sav'
     ftest_new = file_test(xfer_file) *  (1 - file_test(xfer_file, /zero_length))

     if ftest_new eq 0 or keyword_set(refresh) then begin
        xfer_file_old = file_path+'_xferfn_'+file_identifier+'_old.sav'
        ftest_old = file_test(xfer_file_old) *  (1 - file_test(xfer_file_old, /zero_length))

        if ftest_old eq 0 or keyword_set(refresh) then $
           visibility_transfer_generate,xfer_fn,/new,restore_last=0,timing=t_xfer_gen0,polarization=0, $
                                        file_identifier=file_identifier, freq_inds_to_use = inds_to_use  $
        else begin
           restore, xfer_file_old
           file_name_base='_xferfn_'+file_identifier
           xfer_fn=visibility_transfer_convert(xfer_fn,restore_last=0,file_name_base=file_name_base,psf_dim=psf_dim)
        endelse
        heap_free, xfer_fn
        xfer_fn=0
        file_delete, xfer_file_old, /allow_nonexistent
     endif

     restore, xfer_file

     ;; get tile beam in image frame -- needed for going to true sky
     beam_base_uv=fltarr(dimension,elements)
     beam_base_uv[dimension/2.-Floor(psf_dim/2.):dimension/2.-Floor(psf_dim/2.)+psf_dim-1, $
                  elements/2.-Floor(psf_dim/2.):elements/2.-Floor(psf_dim/2.)+psf_dim-1] = *psf_base[0,Floor(nfreq_per_image/2.),0,0]
     beam_base=fft_shift(real_part(FFT(fft_shift(beam_base_uv),/inverse)))
     beam_norm = max(beam_base)
     beam_base/=beam_norm

     ;; weights: all 1's on the uv plane
     norm_model_uv = make_array(dimension, elements, /dcomplex, value=1d)
     weights=visibility_transfer_apply(norm_model_uv, xfer_fn, /new)
     if max(abs(imaginary(weights))) eq 0 then weights = real_part(weights) $
     else stop

     ;; normalize such that the integral of the weights is equal to the number of gridded baselines
     ;; should eventually include the effective area of the tiles...
     normalization = double(nbaselines_per_xferfn) / (total(abs(weights), /double) * kbinsize^2d)
     weights = weights * normalization

     if not keyword_set(no_sim) then begin
        sim_uv = visibility_transfer_apply(my_model_uv, xfer_fn, /new)
        uvf_cube[*, *, i] = sim_uv * normalization
     endif

     if keyword_set(eor_gen) then begin
        eor_uv = visibility_transfer_apply(eor_uvf[*,*,i], xfer_fn, /new)
        eor_cube[*, *, i] = eor_uv * normalization
     endif

     if n_elements(test_power_shape) ne 0 then begin
        test_power_uv = visibility_transfer_apply(test_power_uvf[*,*,i], xfer_fn, /new)
        test_power_cube[*, *, i] = test_power_uv * normalization
     endif

     heap_free, xfer_fn
     xfer_fn=0

     if keyword_set(noise_gen) then begin
        if i eq 0 then seed = systime(1)
        noise_vis_arr = randomn(seed, nfreq_per_image, nbaselines) + randomn(seed, nfreq_per_image, nbaselines) * complex(0,1)
        noise_uv = visibility_grid(noise_vis_arr, polarization=0, freq_inds_to_use = inds_to_use)
        noise_cube[*, *, i] = noise_uv
     endif

     this_freq = mean(frequency_array[inds_to_use]) / 1e6 ;; in MHz

     weights_cube[*, *, i] = weights
     frequencies[i] = this_freq
     tile_beam_image[*,*,i] = beam_base
     normalizations[i] = normalization

  endfor
  if not keyword_set(no_save) then begin
     save, file = info_file, frequencies, degpix, freq_resolution, max_baseline, normalizations
     save, file = weights_file, weights_cube, tile_beam_image

     if not keyword_set(no_sim) then save, file=sim_file, uvf_cube, my_sources, my_model_uv, norm_gridding
     if keyword_set(noise_gen) then save, file = noise_file, noise_cube
     if keyword_set(eor_gen) then save, file = eor_file, eor_cube
     if n_elements(test_power_shape) ne 0 then save, file = test_power_file, test_power_cube

  endif

end
