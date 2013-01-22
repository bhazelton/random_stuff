
function sim_xfer_apply, uv_in, t32 = t32, define_baselines = define_baselines, baseline_spacing = baseline_spacing, $
                         baseline_layout = baseline_layout, fine_res = fine_res, use_outliers = use_outliers

  if keyword_set(define_baselines) then begin
     if n_elements(baseline_spacing) eq 0 then baseline_spacing = 4
     if n_elements(baseline_layout) eq 0 then baseline_layout = 'simple'

     spacing_name = number_formatter(string(baseline_spacing), flag = 2)
     tag = baseline_layout + spacing_name + '_'
     filename = 'Simulation_' + baseline_layout + spacing_name

  endif else if keyword_set(t32) then tag = '32t_' else begin
     filename = 'Simulation'

     if keyword_set(use_outliers) then begin
        filename = filename + '_512t'
        tag = '512t_' 
     endif else begin
        tag = '496t_'
        filename = filename + '_496t'
     endelse
     if keyword_set(fine_res) then begin
        tag = tag + 'fine_'
        filename = filename + '_fine'
     endif
  endelse
  info_file = rootdir('mwa') + 'simulations/' + 'sim_' + tag + 'info.idlsave'

  restore, info_file
  delta_freq = 0.04d
  nfreq_per_image = round(freq_resolution / delta_freq)
  nfreq = n_elements(frequencies)
  dims = [size(uv_in, /dimension), nfreq]

  if n_elements(normalizations) eq 0 then norm_gen = 1 else norm_gen = 0

  data_directory='data'
  file_path=filepath(filename,root_dir=rootdir('mwa'),subdir=data_directory)
  pol_names=['xx','yy','xy','yx']

  if norm_gen eq 1 then begin
     print, 'No normalizations, have to regenerate'

     normalizations = dblarr(nfreq)

     ;; source array for weights
     norm_source = reform([0, 0, 1], 3, 1) ;; columns are: x location (zero centered), y location (zero centered), flux
     n_sources=1
     norm_source_array = reform(fltarr(7, n_sources), 7, n_sources)
     norm_source_array[0,*]=norm_source[0,*] + dims[0]/2.
     norm_source_array[1,*]=norm_source[1,*] + dims[1]/2.
     norm_source_array[4,*]=norm_source[2,*] + dims[1]/2.
     norm_source_array[5,*]=Round(norm_source_array[0,*])+Round(norm_source_array[1,*])*dims[0]
  endif
  
  uvf_cube = make_array(dims, /complex)
  for i=0, nfreq-1 do begin
     inds_to_use = indgen(nfreq_per_image) + nfreq_per_image * i
     ;;print, 'frequency_bins: ' + string(min(inds_to_use), format = '(i3.3)') + ' to ' + string(max(inds_to_use), format = '(i3.3)')

     file_identifier=String(format='(A2,"_f",A3,"_f",A3)',pol_names[0], string(min(inds_to_use), format = '(i3.3)'), $
                            string(max(inds_to_use), format = '(i3.3)'))
     xfer_file = file_path+'_xferfn_'+file_identifier+'.sav'
     ftest_new = file_test(xfer_file) *  (1 - file_test(xfer_file, /zero_length))
     if ftest_new eq 0 then message, 'xfer file does not exist: ' + xfer_file

     restore, xfer_file

     if norm_gen eq 1 then begin
        ;; weights: 1 jy source on center
        norm_model_uv = visibility_source_uv_grid(norm_source_array, u_dim = dims[0], v_dim = dims[1]) 
        weights=visibility_transfer_apply(norm_model_uv, xfer_fn, /new)
        
        ;; need to make image to get normalization (to get to Jy)
        ;; divide by (tile beam)^2 to get to true sky
        norm_image = dirty_image_generate(weights)
        beam2 = (tile_beam_image[*,*,i])^2d
        wh = where(beam2 eq 0, count)
        true_norm_image = norm_image / beam2
        if count gt 0 then true_norm_image[wh] = 0
        
        normalizations[i] = 1./true_norm_image[dims[0]/2, dims[1]/2]
     endif

     sim_uv = visibility_transfer_apply(uv_in, xfer_fn, /new)
     uvf_cube[*, *, i] = sim_uv * normalizations[i]
  endfor

  if norm_gen eq 1 then save, file = info_file, frequencies, degpix, freq_resolution, max_baseline, normalizations

  return, uvf_cube

end
