pro test_eor_sim, delta_uv = delta_uv, uv_max = uv_max, f_avg = f_avg, uv_avg = uv_avg

  freq_arr = findgen(384)*0.08+120
  n_freq = n_elements(freq_arr)
  
  if n_elements(delta_uv) eq 0 then delta_uv = 4.
  if n_elements(uv_max) eq 0 then uv_max = 200.
  uv_arr = findgen(uv_max*2./delta_uv + 1)*delta_uv-uv_max
  n_uv = n_elements(uv_arr)
  
  title = 'delta uv: ' + number_formatter(delta_uv) + ', max uv: ' + number_formatter(uv_max)
  eor_uvf_cube = eor_sim(uv_arr, uv_arr, freq_arr)
  
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
  
  k_arr = sqrt(rebin(kx_mpc, n_kx, n_ky, n_kz, /sample)^2. + rebin(reform(ky_mpc, 1, n_ky), n_kx, n_ky, n_kz, /sample)^2. + $
    rebin(reform(kz_mpc, 1, 1, n_kz), n_kx, n_ky, n_kz, /sample)^2.)
    
  signal_mk = eor_uvf_cube
  ;; convert from Jy -> mK Mpc^2
  for i=0, n_kz-1 do signal_mk[*,*,i] = signal_mk[*,*,i] * conv_factor[i]
  
  signal_mk_k = fft(temporary(signal_mk), dimension = 3) * z_mpc_delta * n_kz
  ;; shift
  signal_mk_k = shift(temporary(signal_mk_k), [0,0,n_kz/2])
  
  hist = histogram(alog10(k_arr), binsize = log_binsize, min = min(alog10(k_centers))-log_binsize/2., locations = locs, reverse_indices=ri)
  new_power2 = fltarr(n_elements(hist))
  for i=0, n_elements(hist)-1 do if hist[i] gt 1 then new_power2[i] = mean(abs(signal_mk_k[ri[ri[i]:ri[i+1]-1]])^2.)
  
  ;; adjust for window function
  window_int = ((2.*!pi)^3./(kx_mpc_delta*ky_mpc_delta*kz_mpc_delta))
  new_power2 = new_power2 / window_int
  print, 'window_int', window_int
  
  cgplot, k_centers, power, psym=4, /ylog, /xlog, title = title
  cgplot, 10.^(locs+log_binsize/2.), new_power2, psym=6, /over, color='red'
  
  
  ratio2 = new_power2/power
  
end
