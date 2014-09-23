pro test_spectral_window, spec_window_type = spec_window_type, n_sets = n_sets, no_noise = no_noise

  if n_elements(n_sets) eq 0 then n_sets = 2
  if n_sets lt 1 or n_sets gt 2 then message, 'n_sets must be set to 1 or 2'
  
  n_trials = 1000
  n_vals = 100
  
  x_delta = 2
  x_vals = findgen(n_vals)*x_delta - n_vals*x_delta/2.
  
  k_delta = 2.*!pi/(n_vals*x_delta)
  k_vals_orig = findgen(n_vals)*k_delta - n_vals*k_delta/2.
  k_ind = round(k_vals_orig / k_delta)
  k_vals_orig[where(k_ind eq 0)] = 0
  
  power_in = fltarr(n_vals) + 1.
  ;power_in = abs(k_ind)/float(n_vals) +1.
  ;; adjust for measurement volume
  signal2 = power_in * (2.*!pi)/(k_delta)
  signal_amp = sqrt(signal2/2.)
  
  signal_real = randomn(seed, n_vals, n_trials) * rebin(signal_amp, n_vals, n_trials, /sample)
  signal_imaginary = randomn(seed, n_vals, n_trials) * rebin(signal_amp, n_vals, n_trials, /sample)
  signal_k = temporary(signal_real) + complex(0,1) * temporary(signal_imaginary)
  
  signal_x = fft(shift(signal_k, n_vals/2), dimension=1, /inverse) * k_delta / (2.*!pi)
  
  if not keyword_set(no_noise) then begin
    sigma_x = fltarr(n_vals) + 0.5
    noise_real = randomn(seed, n_vals, n_trials) * rebin(sigma_x, n_vals, n_trials, /sample)
    noise_imaginary = randomn(seed, n_vals, n_trials) * rebin(sigma_x, n_vals, n_trials, /sample)
    noise1_x = temporary(noise_real) + complex(0,1) * temporary(noise_imaginary)
    
    if n_sets eq 2 then begin
      noise_real = randomn(seed, n_vals, n_trials) * rebin(sigma_x, n_vals, n_trials, /sample)
      noise_imaginary = randomn(seed, n_vals, n_trials) * rebin(sigma_x, n_vals, n_trials, /sample)
      noise2_x = temporary(noise_real) + complex(0,1) * temporary(noise_imaginary)
    endif
  endif else begin
    sigma_x = fltarr(n_vals) + 1
    noise1_x = fltarr(n_vals, n_trials)
    if n_sets eq 2 then noise2_x = noise1_x
  endelse
  
  data1_x = signal_x + noise1_x
  if n_sets eq 2 then data2_x = signal_x + noise2_x
  ;; adjust for measurement volume
  data1_x = data1_x / sqrt((2.*!pi)/(k_delta))
  if n_sets eq 2 then data2_x = data2_x / sqrt((2.*!pi)/(k_delta))
  var_x = sigma_x^2. / ((2.*!pi)/(k_delta))
  
  if n_sets eq 1 then begin
    data_sum = data1_x
    var_sum = var_x
  endif else begin
    data_sum = (data1_x + data2_x)/2.
    data_diff = (data1_x - data2_x)/2.
    var_sum = var_x/2.
  endelse
  
  ;; apply spectral windowing function if desired
  if n_elements(spec_window_type) ne 0 then begin
    window = spectral_window(n_vals, type = spec_window_type, /periodic)
    
    norm_factor = sqrt(n_vals/total(window^2.))
    
    window = window * norm_factor
    
    window_expand = rebin(window, n_vals, n_trials, /sample)
    
    data_sum = data_sum * window_expand
    if n_sets eq 2 then data_diff = data_diff * window_expand
    
    var_sum = var_sum * temporary(window_expand^2.)
  endif
  
  
  data_sum_ft = shift(fft(data_sum, dimension=1), n_vals/2) * n_vals * x_delta
  if n_sets eq 2 then data_diff_ft = shift(fft(data_diff, dimension=1), n_vals/2) * n_vals * x_delta
  
  
  k_vals = k_vals_orig[where(k_ind ge 0)]
  n_k = n_elements(k_vals)
  
  zero_ind = where(k_ind eq 0)
  pos_inds = where(k_ind gt 0)
  neg_inds = reverse(where(k_ind lt 0 and abs(k_ind) le max(k_ind)))
  a1_0 = data_sum_ft[zero_ind,*]
  a1_n = (data_sum_ft[pos_inds,*] + data_sum_ft[neg_inds,*])/2.
  b1_n = complex(0,1) * (data_sum_ft[pos_inds,*] - data_sum_ft[neg_inds,*])/2.
  
  if n_sets eq 2 then begin
    a2_0 = data_diff_ft[zero_ind,*]
    a2_n = (data_diff_ft[pos_inds,*] + data_diff_ft[neg_inds,*])/2.
    b2_n = complex(0,1) * (data_diff_ft[pos_inds,*] - data_diff_ft[neg_inds,*])/2.
  endif
  
  data_sum_cos = complex(fltarr(n_k, n_trials))
  data_sum_sin = complex(fltarr(n_k, n_trials))
  data_sum_cos[0,*] = a1_0
  data_sum_cos[1:n_k-1,*] = a1_n
  data_sum_sin[1:n_k-1,*] = b1_n
  
  if n_sets eq 2 then begin
    data_diff_cos = complex(fltarr(n_k, n_trials))
    data_diff_sin = complex(fltarr(n_k, n_trials))
    data_diff_cos[0,*] = a2_0
    data_diff_cos[1:n_k-1,*] = a2_n
    data_diff_sin[1:n_k-1,*] = b2_n
  endif
  
  if max(abs(var_sum)) gt 0 then begin
    covar_cos = fltarr(n_k, n_trials)
    covar_sin = fltarr(n_k, n_trials)
    covar_cross = fltarr(n_k, n_trials)
    
    x_k_arr = rebin(reform(k_vals, 1, n_k), n_vals, n_k) * rebin(x_vals, n_vals, n_k)
    
    cos_arr = cos(x_k_arr)
    sin_arr = sin(x_k_arr)
    
    sigma2 = rebin(var_sum, n_vals, n_trials)
    
    covar_cos = transpose(matrix_multiply(sigma2, cos_arr^2d, /atranspose)) * (x_delta)^2.
    covar_sin = transpose(matrix_multiply(sigma2, sin_arr^2d, /atranspose)) * (x_delta)^2.
    covar_cross = transpose(matrix_multiply(sigma2, cos_arr*sin_arr, /atranspose)) * (x_delta)^2.
    
    theta = atan(2.*covar_cross, covar_cos - covar_sin)/2.
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    
    sigma2_1 = covar_cos*cos_theta^2. + 2.*covar_cross*cos_theta*sin_theta + covar_sin*sin_theta^2.
    sigma2_2 = covar_cos*sin_theta^2. - 2.*covar_cross*cos_theta*sin_theta + covar_sin*cos_theta^2.
    
    data_sum_1 = data_sum_cos*cos_theta + data_sum_sin*sin_theta
    data_sum_2 = (-1d)*data_sum_cos*sin_theta + data_sum_sin*cos_theta
    
    if n_sets eq 2 then begin
      data_diff_1 = data_diff_cos*cos_theta + data_diff_sin*sin_theta
      data_diff_2 = (-1d)*data_diff_cos*sin_theta + data_diff_sin*cos_theta
    endif
  endif else begin
    data_sum_1 = data_sum_cos
    data_sum_2 = data_sum_sin
    
    if n_sets eq 2 then begin
      data_diff_1 = data_diff_cos
      data_diff_2 = data_diff_sin
    endif
  endelse
  
  if n_elements(sigma2_1) gt 0 then begin
    power_weights1 = 1d/(4*(sigma2_1)^2d)
    wh_sig1_0 = where(sigma2_1^2d eq 0, count_sig1_0)
    if count_sig1_0 ne 0 then power_weights1[wh_sig1_0] = 0
    
    power_weights2 = 1d/(4*(sigma2_2)^2d) ;; inverse variance
    wh_sig2_0 = where(sigma2_2^2d eq 0, count_sig2_0)
    if count_sig2_0 ne 0 then power_weights2[wh_sig2_0] = 0
  endif else begin
    power_weights1 = fltarr(n_k, n_trials)+1
    power_weights2 = fltarr(n_k, n_trials)+1
  endelse
  
  
  if n_sets eq 2 then begin
    term1 = (abs(data_sum_1)^2. - abs(data_diff_1)^2.) * power_weights1
    term2 = (abs(data_sum_2)^2. - abs(data_diff_2)^2.) * power_weights2
    
    noise_t1 = abs(data_diff_1)^2. * power_weights1
    noise_t2 = abs(data_diff_2)^2. * power_weights2
  endif else begin
    term1 = abs(data_sum_1)^2. * power_weights1
    term2 = abs(data_sum_2)^2. * power_weights2
  endelse
  
  
  ;; Factor of 2 because we're adding the cosine & sine terms
  noise_expval = (sqrt(power_weights1 + power_weights2))*2.
  ;; except for kparallel=0 b/c there's only one term
  noise_expval[0,*] = noise_expval[0,*]/2.
  
  weights = (power_weights1 + power_weights2)
  
  ;; multiply by 2 because power is generally the SUM of the cosine & sine powers
  power = (term1 + term2)*2.
  if n_sets eq 2 then noise = (noise_t1 + noise_t2)*2
  ;; except for kparallel=0 b/c there's only one term
  power[0,*] = power[0,*]/2.
  if n_sets eq 2 then noise[0,*] = noise[0,*]/2.
  
  power = power / weights
  if n_sets eq 2 then noise = noise / weights
  noise_expval = noise_expval / weights
  
  wh_wt0 = where(weights eq 0, count_wt0)
  if count_wt0 ne 0 then begin
    power[wh_wt0] = 0
    if n_sets eq 2 then noise[wh_wt0] = 0
    noise_expval[wh_wt0] = 0
  endif
  
  ;; variance = 4/weights b/c of factors of 2 in power
  var = 4./weights[*,0]
  ;; except for kparallel=0 b/c there's only one term
  var[0] = var[0]/4.
  
  power_1d = total(power, 2)/n_trials
  if n_sets eq 2 then begin
    noise_1d = total(noise, 2)/n_trials
    var_1d = total((noise - rebin(noise_1d, n_k, n_trials))^2.,2)/(n_trials-1)
  endif
  noise_exp_1d = noise_expval[*,0]
  
  if windowavailable(1) then wset, 1 else window, 1
  cgplot, k_vals_orig, power_in, yrange = minmax([0,power_in, power_1d]), title='power', color='blue'
  cgplot, k_vals, power_1d, /over, color='red', psym=-4
  
  if n_sets eq 2 then begin
    if windowavailable(2) then wset, 2 else window, 2
    cgplot, var, yrange=minmax([0, var, var_1d]), title = 'Variance', color='blue'
    cgplot, var_1d, /over, color='red', psym=-4
  endif
  
  
end
