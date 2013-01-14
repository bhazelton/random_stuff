

function eor_sim, kx_arr, ky_arr, freq_arr, seed = seed, flat_sigma = flat_sigma

  if n_elements(seed) eq 0 then seed = systime(1)

  kx_diff = kx_arr - shift(kx_arr, 1)
  kx_diff = kx_diff[0:*]

  ky_diff = ky_arr - shift(ky_arr, 1)
  ky_diff = ky_diff[0:*]

  z0_freq = 1420.40 ;; MHz
  redshifts = z0_freq/freq_arr - 1
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los
     
  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  z_mpc_delta = mean(comov_los_diff)
  z_mpc_mean = mean(comov_dist_los)

  z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
  kz_mpc_range =  (2d*!pi) / z_mpc_delta ;; no factor of 2 -- need as many kz as freq
  kz_mpc_delta = (2d*!pi) / z_mpc_length
  kz_mpc_final = dindgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - (kz_mpc_range/2-kz_mpc_delta)
  ;;kz_mpc = kz_mpc_final[n_elements(kz_mpc_final)/2-1:*]
  kz_mpc = temporary(kz_mpc_final)

  print, 'z delta: ', z_mpc_delta
  print,  'kz delta: ', kz_mpc_delta

  ;; convert kx_arr & ky_arr (in wavelengths) to inverse comoving Mpc
  kx_mpc = kx_arr * (2d*!pi)/ z_mpc_mean
  kx_mpc_delta = kx_diff[0] * (2d*!pi)/ z_mpc_mean

  ky_mpc = ky_arr * (2d*!pi) / z_mpc_mean
  ky_mpc_delta = ky_diff[0] * (2d*!pi) / z_mpc_mean

  n_kx = n_elements(kx_mpc)
  n_ky = n_elements(ky_mpc)
  n_kz = n_elements(kz_mpc)

  restore, base_path() + 'power_spectrum/eor_data/eor_power_1d.idlsave' ;;k_centers, power

  npts_log = n_elements(k_centers)

  log_diff =  alog10(k_centers) - shift(alog10(k_centers), 1)
  log_diff = log_diff[1:*]
  log_binsize = log_diff[0]

  if n_elements(flat_sigma) ne 0 then begin
     power_3d = dblarr(n_kx, n_ky, n_kz) + flat_sigma

  endif else begin

     k_arr = sqrt(rebin(kx_mpc, n_kx, n_ky, n_kz)^2d + rebin(reform(ky_mpc, 1, n_ky), n_kx, n_ky, n_kz)^2d + $
                  rebin(reform(kz_mpc, 1, 1, n_kz), n_kx, n_ky, n_kz)^2d)
     wh0 = where(k_arr eq 0, count)
     if count ne 0 then k_arr[wh0] = min(k_centers)
     
     result = 10^(interpol(alog10(power), alog10(k_centers), alog10(k_arr)))
     
     power_3d = reform(temporary(result), n_kx, n_ky, n_kz)
     
     mu = rebin(reform(abs(kz_mpc), 1, 1, n_kz), n_kx, n_ky, n_kz) / temporary(k_arr)
     power_3d = power_3d * (1 + 2 * mu^2d + mu^4d)
     
     mu=0
  endelse


  bins_per_decade = 8
  power_2d = kspace_rebinning_2d(power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, bins = bins_per_decade, $
                                 binned_weights = weights_2d, fill_holes=0)
  
  temp = power
  power = power_2d
  weights = weights_2d  
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc

  filename = base_path() + 'power_spectrum/eor_data/sim_2d.idlsave'
  save, file = filename, power, weights, kperp_edges, kpar_edges
  power = temp

  kperp_plot_range = [min(kperp_edges), 0.3]
  plotfile = base_path() + 'power_spectrum/plots/eor_sim_initial_kspace_power'
  kpower_2d_plots, filename, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = [1e4, 1e8], $
                   title = 'EoR Power', pub = pub, plotfile = plotfile


  bins_per_decade = 10d
  power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, bins = bins_per_decade, $
                                 binned_weights = weights_1d)

  temp = power
  power = power_1d
  weights = weights_1d  
  k_edges = k_edges_mpc

  filename = base_path() + 'power_spectrum/eor_data/sim_1d.idlsave'
  save, file = filename, power, weights, k_edges, bins_per_decade
  power = temp

  eor_file_1d = base_path() + 'power_spectrum/eor_data/eor_power_1d.idlsave'
  file_arr = [filename, eor_file_1d]
  names_arr = ['Constructed EoR cube', 'EoR signal']
  colors_arr = [0, 254]
  kpower_1d_plots, file_arr, window_num = 2, names = names_arr, colors = colors_arr


  signal_amp = sqrt(power_3d)
  signal_phase = randomu(seed, n_kx, n_ky, n_kz) * 2d * !pi

  signal = temporary(signal_amp) * exp(dcomplex(0,1) * temporary(signal_phase))

  ;; shift it so that it's as expected when we take the fft
  signal = shift(temporary(signal), [0,0,n_kz/2+1])

  print, 'signal^2d integral:', total(abs(signal)^2d)
  print, 'signal^2d integral * 2pi*delta_k^2d:', total(abs(signal)^2d) * kz_mpc_delta * 2d * !dpi

  ;;temp = conj(reverse(signal[*,*,1:n_kz-2],3))
  ;;signal = [[[signal]], [[temp]]]

  ;; fourier transform along z direction to get to uvf space
  temp = fft(temporary(signal), dimension = 3, /inverse) * kz_mpc_delta

  print, 'sum(uvf signal^2)*z_delta:', total(abs(temp)^2d)*z_mpc_delta

  signal = temp
  return, signal
end
