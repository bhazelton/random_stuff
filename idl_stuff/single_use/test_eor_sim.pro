pro test_eor_sim, delta_uv = delta_uv, uv_max = uv_max, f_avg = f_avg, uv_avg = uv_avg, spec_window_type = spec_window_type, $
    apply_beam = apply_beam, use_sim = use_sim, flat_sigma = flat_sigma, save_cubefile = save_cubefile
    
    
  if keyword_set(use_sim) then begin
    restore, base_path('data') +'fhd_sim_data/snap_highf_eor_nomu_newconv/1061316176_input_model.sav' ; model_uvf, uv_arr, freq_arr
    eor_uvf_cube = temporary(model_uvf)*2.
    n_freq = n_elements(freq_arr)
    n_uv = n_elements(uv_arr)
    
    delta_uv = uv_arr[1] - uv_arr[0]
    uv_max = max(uv_arr)
    title = 'delta uv: ' + number_formatter(delta_uv) + ', max uv: ' + number_formatter(uv_max)
    
    
    if keyword_set(apply_beam) then begin
      ;      beam2_image = getvar_savefile(base_path('data') +'fhd_sim_data/snap_highf_eor_nomu_newconv/1061316176_initial_beam2_image.sav', 'beam2_xx_image')
      stop
      beam2_image = getvar_savefile(base_path('data') +'fhd_sim_data/snap_highf_eor_nomu_newconv/1061316176_initial_beam2_image.sav', beam2_xx_image)
      
      temp = shift(fft(fft(shift(eor_uvf_cube,n_uv/2,n_uv/2,0), dimension=1),dimension=2),n_uv/2,n_uv/2,0)
      stop
      temp = temp * sqrt(beam2_image)
      stop
      undefine, beam2_image
      temp = shift(fft(fft(shift(temp,n_uv/2,n_uv/2,0), dimension=1, /inverse), dimension=2, /inverse),n_uv/2,n_uv/2,0)
      stop
      eor_uvf_cube = temp
    endif
    
  endif else begin
    freq_arr = findgen(384)*0.08+120
    n_freq = n_elements(freq_arr)
    
    if n_elements(delta_uv) eq 0 then delta_uv = 4.
    if n_elements(uv_max) eq 0 then uv_max = 200.
    uv_arr = findgen(uv_max*2./delta_uv + 1)*delta_uv-uv_max
    n_uv = n_elements(uv_arr)
    
    title = 'delta uv: ' + number_formatter(delta_uv) + ', max uv: ' + number_formatter(uv_max)
    eor_uvf_cube = eor_sim(uv_arr, uv_arr, freq_arr, flat_sigma = flat_sigma)
    
    print, 'stddev of uvf_cube: ' + number_formatter(stddev(abs(eor_uvf_cube)), format='(e10.2)')
    
    if keyword_set(apply_beam) then begin
    
      degpix = !radeg / uv_max
      theta_vals = (findgen(n_uv)-n_uv/2) * degpix
      sigma_beam = max(theta_vals)/(sqrt(2.*alog(2)))
      gaussian_theta = exp(-1 * theta_vals^2./(2.*sigma_beam^2.))
      gaussian_theta = gaussian_theta / max(gaussian_theta)
      beam_tophat = fltarr(n_uv,n_uv,n_freq)+1.
      beam = rebin(gaussian_theta, n_uv, n_uv, n_freq, /sample) * rebin(reform(gaussian_theta, 1, n_uv), n_uv, n_uv, n_freq, /sample)
      
      temp = shift(fft(fft(shift(eor_uvf_cube,n_uv/2,n_uv/2,0), dimension=1),dimension=2),n_uv/2,n_uv/2,0)
      temp = temp * beam
      temp = shift(fft(fft(shift(temp,n_uv/2,n_uv/2,0), dimension=1, /inverse), dimension=2, /inverse),n_uv/2,n_uv/2,0)
      
      eor_uvf_cube = temp
    endif
    
  endelse
  
  if n_elements(f_avg) gt 0 then if f_avg gt 1 then begin
    nf_new = floor(n_freq / f_avg)
    temp = complex(fltarr(n_uv, n_uv, nf_new))
    temp_f = fltarr(nf_new)
    for i=0, nf_new-1 do begin
      temp[*,*,i] = total(eor_uvf_cube[*,*,i*f_avg:(i+1)*f_avg-1], 3) / f_avg
      temp_f[i] = total(freq_arr[i*f_avg:(i+1)*f_avg-1]) / f_avg
    endfor
    
    title = title + ', f avg: ' + number_formatter(f_avg)
    
    eor_uvf_cube = temporary(temp)
    freq_arr = temporary(temp_f)
    n_freq = nf_new
  endif
  
  
  if n_elements(uv_avg) gt 0 then if uv_avg gt 1 then begin
    nuv_new = floor(n_uv / uv_avg)
    temp = complex(fltarr(nuv_new, n_uv, n_freq))
    temp_uv = fltarr(nuv_new)
    for i=0, nuv_new-1 do begin
      temp[i,*,*] = total(eor_uvf_cube[i*uv_avg:(i+1)*uv_avg-1,*,*], 1) / uv_avg
      temp_uv[i] = total(uv_arr[i*uv_avg:(i+1)*uv_avg-1]) / uv_avg
    endfor
    
    temp2 = complex(fltarr(nuv_new, nuv_new, n_freq))
    for i=0, nuv_new-1 do temp2[*,i,*] = total(temp[*,i*uv_avg:(i+1)*uv_avg-1,*], 2) / uv_avg
    undefine, temp
    
    title = title + ', uv avg: ' + number_formatter(uv_avg)
    
    eor_uvf_cube = temporary(temp2)
    uv_arr = temporary(temp_uv)
    n_uv = nuv_new
  endif
  
  restore, filepath('eor_power_1d.idlsave',root=rootdir('FHD'),subdir='catalog_data')
  npts_log = n_elements(k_centers)
  
  log_diff =  alog10(k_centers) - shift(alog10(k_centers), 1)
  log_diff = log_diff[1:*]
  log_binsize = log_diff[0]
  
  
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
  
  
  signal_mk = eor_uvf_cube
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
  
  k_arr_2d = sqrt(reform(rebin(kx_mpc, n_kx, n_ky, /sample)^2. + rebin(reform(ky_mpc, 1, n_ky), n_kx, n_ky, /sample)^2.,n_kx*n_ky))
  hist_kpar = histogram(k_arr_2d, binsize = kx_mpc_delta, min = 0, locations = locs_kpar, reverse_indices=ri_kpar)
  temp_power = fltarr(n_elements(hist_kpar), n_kz)
  temp_signal = reform(signal_mk_k, n_kx*n_ky, n_freq)
  for i=0, n_elements(hist_kpar)-1 do begin
    if hist_kpar[i] eq 0 then continue
    if hist_kpar[i] gt 1 then temp_power[i, *] = total(abs(temp_signal[ri_kpar[ri_kpar[i]:ri_kpar[i+1]-1], *])^2.,1)/hist_kpar[i] $
    else temp_power[i, *] = abs(temp_signal[ri_kpar[ri_kpar[i]:ri_kpar[i+1]-1], *])^2.
  endfor
  undefine, temp_signal
  kz_bin = abs(round(kz_mpc / kz_mpc_delta))
  hist_kz = histogram(kz_bin, binsize = 1, min=0, locations = locs_kz, reverse_indices=ri_kz)
  locs_kz = locs_kz*kz_mpc_delta
  power_2d = fltarr(n_elements(hist_kpar), n_elements(hist_kz))
  for i=0, n_elements(hist_kz)-1 do begin
    if hist_kz[i] eq 0 then continue
    if hist_kz[i] gt 1 then power_2d[*,i] = total(temp_power[*, ri_kz[ri_kz[i]:ri_kz[i+1]-1]],2)/hist_kz[i] else power_2d[*,i] = temp_power[*, ri_kz[ri_kz[i]:ri_kz[i+1]-1]]
  endfor
  
  undefine, temp_power
  
  k_arr = sqrt(rebin(kx_mpc, n_kx, n_ky, n_kz, /sample)^2. + rebin(reform(ky_mpc, 1, n_ky), n_kx, n_ky, n_kz, /sample)^2. + $
    rebin(reform(kz_mpc, 1, 1, n_kz), n_kx, n_ky, n_kz, /sample)^2.)
    
  hist = histogram(alog10(k_arr), binsize = log_binsize, min = min(alog10(k_centers))-log_binsize/2., locations = locs, reverse_indices=ri)
  new_power2 = fltarr(n_elements(hist))
  for i=0, n_elements(hist)-1 do if hist[i] gt 1 then new_power2[i] = mean(abs(signal_mk_k[ri[ri[i]:ri[i+1]-1]])^2.)
  
  print, 'mean of power in mK (before window divided out): ' + number_formatter(mean(new_power2), format='(e10.2)')
  
  ;; adjust for window function
  if keyword_set(apply_beam) then begin
    xy_mpc_delta = (2.*!pi) / (n_uv * kx_mpc_delta)
    
    window_int = total(beam[*,*,0]^2.*xy_mpc_delta^2.) * ((2.*!pi)/kz_mpc_delta)
    window_int_tophat = total(beam_tophat[*,*,0]^2.*xy_mpc_delta^2.) * ((2.*!pi)/kz_mpc_delta)
    volume = ((2.*!pi)^3./(kx_mpc_delta*ky_mpc_delta*kz_mpc_delta))
  endif else begin
    window_int = ((2.*!pi)^3./(kx_mpc_delta*ky_mpc_delta*kz_mpc_delta))
  endelse
  new_power2 = new_power2 / window_int
  print, 'window_int', window_int
  
  print, 'mean of power in mK (final): ' + number_formatter(mean(new_power2), format='(e10.2)')
  
  if windowavailable(1) then wset, 1 else window, 1
  
  cgplot, k_centers, power, psym=4, /ylog, /xlog, title = title
  cgplot, 10.^(locs+log_binsize/2.), new_power2, psym=6, /over, color='red'
  if keyword_set(flat_sigma) then cgplot, k_centers, (power*0d + max(power)), psym=4, /over, color='blue'
  
  kpower_2d_plots, power_savefile, power = power_2d, kperp_edges = [locs_kpar-kx_mpc_delta/2., max(locs_kpar) + kx_mpc_delta/2.], $
    kpar_edges = [locs_kz-kz_mpc_delta/2., max(locs_kz) + kz_mpc_delta/2.], window_num=2
    
  ratio2 = new_power2/power
  
  if n_elements(save_cubefile) ne 0 then begin
  
    model_uv_arr = temporary(eor_uvf_cube)/2. ;; divide by 2 to get to xx polarization
    weights_uv_arr = fltarr(n_kx, n_ky, n_kz) +1.
    variance_uv_arr = weights_uv_arr
    beam2_image = weights_uv_arr
    
    if keyword_set(apply_beam) then begin
      temp = shift(fft(fft(shift(weights_uv_arr,n_uv/2,n_uv/2,0), dimension=1),dimension=2),n_uv/2,n_uv/2,0)
      temp = temp * beam
      temp = shift(fft(fft(shift(temp,n_uv/2,n_uv/2,0), dimension=1, /inverse), dimension=2, /inverse),n_uv/2,n_uv/2,0)
      weights_uv_arr = temp
      
      temp = shift(fft(fft(shift(variance_uv_arr,n_uv/2,n_uv/2,0), dimension=1),dimension=2),n_uv/2,n_uv/2,0)
      temp = temp * beam^2.
      temp = shift(fft(fft(shift(temp,n_uv/2,n_uv/2,0), dimension=1, /inverse), dimension=2, /inverse),n_uv/2,n_uv/2,0)
      variance_uv_arr = temp
      
      beam2_image = beam^2.
    endif
    
    vis_noise = ptr_new(fltarr(1, n_freq) + 1.)
    obs = {max_baseline:uv_max, obsra:0, obsdec:0, zenra:0, zendec:0, n_freq:n_freq, degpix:(180./!pi)/(n_kx*delta_uv), $
      kpix:delta_uv, dimension:n_uv, elements:n_uv, freq:freq_arr*1e6, time_res:2, $
      n_vis:(128*127/2.)*60.*n_freq, nf_vis:fltarr(n_freq)+(128*127/2.)*60., vis_noise:vis_noise}
      
    for i=0, n_elements(save_cubefile)-1 do begin
      save, file=save_cubefile[i], weights_uv_arr, model_uv_arr, variance_uv_arr, beam2_image, obs
    endfor
  endif
  
end
