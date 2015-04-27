pro simple_ps, uv_data, uv_arr, freq_arr, uv_weights = uv_weights, spec_window_type=spec_window_type, $
    beam2_int = beam2_int
  
  if n_elements(uv_weights) eq 0 then uv_weights = uv_data*0. + 1.
    
  delta_u = uv_arr[1] - uv_arr[0]
  delta_v = uv_arr[1] - uv_arr[0]
  
    if min(freq_arr) gt 1e8 then frequencies = freq_arr / 1e6 else frequencies = freq_arr
  f_delta = frequencies[1]-frequencies[0]
  n_kz = n_elements(frequencies)
  
  z0_freq = 1420.40 ;; MHz
  ;; make sure frequencies are in MHz
  
  redshifts = z0_freq/frequencies - 1
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los
  
  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  z_mpc_delta = float(mean(comov_los_diff))
  z_mpc_mean = float(mean(comov_dist_los))
  
  ;; converting from Jy (in u,v,f) to mK*str (10^-26 * c^2 * 10^3/ (2*f^2*kb))
  conv_factor = float((3e8)^2 / (2. * (frequencies*1e6)^2. * 1.38065))
  
  ;; convert from Jy -> mk*str -> mK*Mpc^2
  conv_factor = conv_factor * z_mpc_mean^2.
  
  ;; convert from uv (in wavelengths) to kx/ky in inverse comoving Mpc
  kx_mpc = uv_arr * (2.*!pi) / z_mpc_mean
  kx_mpc_delta = delta_u * (2.*!pi) / z_mpc_mean
  n_kx = n_elements(kx_mpc)
  
  ky_mpc = uv_arr * (2.*!pi) / z_mpc_mean
  ky_mpc_delta = delta_v * (2.*!pi) / z_mpc_mean
  n_ky = n_elements(ky_mpc)
  
  z_mpc_length = float(max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta)
  kz_mpc_range =  (2.*!pi) / (z_mpc_delta)
  kz_mpc_delta = (2.*!pi) / z_mpc_length
  kz_mpc = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2.
  if n_elements(kz_mpc) ne n_kz then stop
  
  
  if keyword_set(calc_al_weights) then begin
    signal_mk = uv_data/al_weights
    wh_wt0 = where(al_weights eq 0, count_wt0, complement=wh_wt_gt0)
    if count_wt0 gt 0 then signal_mk[wh_wt0] = 0
  endif else begin
    signal_mk = uv_data/uv_weights
    wh_wt0 = where(uv_weights eq 0, count_wt0, complement=wh_wt_gt0)
    if count_wt0 gt 0 then signal_mk[wh_wt0] = 0
  endelse
  ;; convert from Jy -> mK Mpc^2
  for i=0, n_kz-1 do signal_mk[*,*,i] = signal_mk[*,*,i] * conv_factor[i]
  
  print, 'stddev of uvf_cube in mK: ' + number_formatter(stddev(abs(signal_mk)), format='(e10.2)')
  
  ;; apply spectral windowing function if desired
  if n_elements(spec_window_type) ne 0 then begin
    window = spectral_window(n_freq, type = spec_window_type, /periodic)
    
    norm_factor = sqrt(n_freq/total(window^2.))
    
    window = window * norm_factor
    
    window_expand = rebin(reform(window, 1, 1, n_freq), n_kx, n_ky, n_freq, /sample)
    
    signal_mk = signal_mk * window_expand
  endif
  
  signal_mk_k = fft(temporary(signal_mk), dimension = 3) * z_mpc_delta * n_kz
  ;; shift
  signal_mk_k = shift(temporary(signal_mk_k), [0,0,n_kz/2])
  
  print, 'stddev of k_cube in mK: ' + number_formatter(stddev(abs(signal_mk_k)), format='(e10.2)')
  
  power_3d = abs(signal_mk_k)^2.
  
  power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, k_bin = k1_bin, log_k = log_k1d, $
    binned_weights = weights_1d)
  k_centers = (k_edges[1:*] - k_edges[0:n_elements(power)-12])/2.
  
  ;; adjust for beam window function
  if n_elements(beam2_int) gt 0 then begin
    window_int = beam2_int * ((2.*!pi)/kz_mpc_delta)
    
    volume = ((2.*!pi)^3./(kx_mpc_delta*ky_mpc_delta*kz_mpc_delta))
  endif else begin
    window_int = ((2.*!pi)^3./(kx_mpc_delta*ky_mpc_delta*kz_mpc_delta))
  endelse
  
  power_1d = power_1d / window_int
  print, 'window_int', window_int
  
  print, 'mean of power in mK (final): ' + number_formatter(mean(power_1d), format='(e10.2)')
  
  
  case strlowcase(!version.os_family) OF
    'windows': split_delim = ';'
    'unix':    split_delim = ':'
  endcase
  path_dirs = strsplit(!path, split_delim, /extract)
  
  fhd_catalog_loc = strpos(path_dirs, 'catalog_data')
  wh_catalog = where(fhd_catalog_loc gt 0, count_catalog)
  if count_catalog gt 0 then begin
    file_path = path_dirs[wh_catalog[0]]
    ;; make sure file_path has a path separator at the end
    pos = strpos(file_path, path_sep(), /reverse_search)
    if pos+1-strlen(file_path) lt 0 then file_path = file_path + path_sep()
    
    eor_file_1d = file_path + 'eor_power_1d.idlsave'
    flat_file_1d = file_path + 'flat_power_1d.idlsave'
    
    flat_kcenters = getvar_savefile(flat_file_1d, 'k_centers')
    flat_power = getvar_savefile(flat_file_1d, 'power')
    
    eor_kcenters = getvar_savefile(eor_file_1d, 'k_centers')
    eor_power = getvar_savefile(eor_file_1d, 'power')
  endif
  
  
  if not keyword_set(no_plots) then begin
    if windowavailable(1) then begin
      wset, 1
      if !d.y_size/!d.x_size gt 2 or !d.x_size/!d.y_size gt 2 then make_win = 1 else make_win = 0
    endif else make_win = 1
    if make_win eq 1 then window, 1
    
    if keyword_set(flat_sigma) then yrange = minmax([max(power), new_power2[where(new_power2 gt 0)], 1e5, 1e7]) else yrange = minmax([power, new_power2[where(new_power2 gt 0)]])
    
    cgplot, eor_kcenters, eor_power, psym=-3, /ylog, /xlog, title = title, yrange = yrange
    cgplot, k_centers, power_1d, psym=6, /over, color='red'
    cgplot, flat_kcenters, flat_power, psym=-3, /over, color='blue'
    
    kpower_2d_plots, power_savefile, power = power_2d, kperp_edges = [locs_kpar-kx_mpc_delta/2., max(locs_kpar) + kx_mpc_delta/2.], $
      kpar_edges = [locs_kz-kz_mpc_delta/2., max(locs_kz) + kz_mpc_delta/2.], window_num=2
  endif
    
   
  
end
