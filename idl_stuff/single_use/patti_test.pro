
function model, a, b, freq, times
  return, a*cos(times*freq) + b*sin(times*freq)
end

function chi2, a, b, freq, times, data, sigma2
  model = model(a, b, freq, times)
  return, total((data-model)^2d/sigma2)
end

function second_deriv_matrix, a0, b0, freq,  times, data, sigma2

  delta = 0.01
  
  ;; first calculate second deriv. in each primary direction
  central_chi2 = chi2(a0, b0, freq, times, data, sigma2)
  xplus_chi2 = chi2(a0+delta, b0, freq, times, data, sigma2)
  xminus_chi2 = chi2(a0-delta, b0, freq, times, data, sigma2)
  xxderiv = (xplus_chi2 - 2*central_chi2 + xminus_chi2) / delta^2d

  yplus_chi2 = chi2(a0, b0+delta, freq, times, data, sigma2)
  yminus_chi2 = chi2(a0, b0-delta, freq, times, data, sigma2)
  yyderiv = (yplus_chi2 - 2*central_chi2 + yminus_chi2) / delta^2d

  ;; now do off diagonal derivative
  xplus_yplus_chi2 = chi2(a0+delta, b0+delta, freq, times, data, sigma2)
  xplus_yminus_chi2 = chi2(a0+delta, b0-delta, freq, times, data, sigma2)
  xminus_yplus_chi2 = chi2(a0-delta, b0+delta, freq, times, data, sigma2)
  xminus_yminus_chi2 = chi2(a0-delta, b0-delta, freq, times, data, sigma2)


  xyderiv = (xplus_yplus_chi2 - xplus_yminus_chi2 - xminus_yplus_chi2 + xminus_yminus_chi2) / (4*delta^2d)

  return, [[xxderiv, xyderiv],[xyderiv,yyderiv]]

end



pro patti_test, new_data = new_data, noise_phase = noise_phase, noise_amp = noise_amp, noise_freq = noise_freq, $
                use_fisher_cut = use_fisher_cut, std_power = std_power

  textfast,patti_data,/read,filename='sprouts.dat',root=base_path()+'/single_use',column=indgen(3),/use_exact
  
  times = reform(patti_data[0,*])
  data = reform(patti_data[1,*])
  sigma2 = reform(patti_data[2,*]^2.)

  npts = n_elements(times)

  t_delta = times[1] - times[0]
  t_length = max(times) - min(times) + t_delta
  freq_delta = 2d*!pi/t_length
  freq_range = 2d*!pi/(2*t_delta) ;; factor of 2 b/c only keep half
  frequencies = dindgen(round(freq_range / freq_delta)) * freq_delta
  n_freq = n_elements(frequencies)
  if n_freq*2 ne npts then stop

  freq = !pi
  wh = where(frequencies eq !pi, count)
  if count eq 1 then freq_ind = wh[0] else stop

  if n_elements(noise_phase) ne 0 or n_elements(noise_amp) ne 0 or n_elements(noise_freq) ne 0 then new_noise = 1 else new_noise = 0

  if keyword_set(new_data) or new_noise eq 1 then begin
     a0 = 1
     b0 = 0
     signal = model(a0, b0, freq, times)
 
     rands = randomn(seed, npts)

     if new_noise eq 1 then begin
        if n_elements(noise_phase) eq 0 then noise_phase = 0
        if n_elements(noise_amp) eq 0 then noise_amp = 0.5
        if n_elements(noise_freq) eq 0 then noise_freq = freq

        noise = abs(noise_amp * cos(noise_freq*times + noise_phase))
        sigma2 = noise^2d
     endif else noise = sqrt(sigma2)

 
     data = signal + noise*rands
  endif

  data_ft = fft(data)
  ;; drop negative frequencies
  data_ft = data_ft[0:n_freq-1]

  sigma2_inv_arr = rebin(1d/sigma2, npts, n_freq)
  freq_time_arr = rebin(reform(frequencies, 1, n_freq), npts, n_freq) * rebin(times, npts, n_freq)
  cos_arr = cos(freq_time_arr)
  sin_arr = sin(freq_time_arr)
  
  ;; norm = 1/double(npts)
  ;; norm = 2/double(npts)^2d
  norm = 1d
  f_cos2 = norm * total(sigma2_inv_arr*cos_arr^2d, 1)
  f_sin2 = norm * total(sigma2_inv_arr*sin_arr^2d, 1)
  f_cross = norm * total(sigma2_inv_arr*cos_arr*sin_arr, 1)

  fisher_rot = [[f_cos2[freq_ind], f_cross[freq_ind]], [f_cross[freq_ind], f_sin2[freq_ind]]]
  covar_rot = invert(fisher_rot, status)
  if status ne 0 then begin
     print, 'Inversion was not successful, using pseudo inverse'
     covar_rot = idl_pseudo_inverse(fisher_rot)
  endif

  fisher_ff_cut = dblarr(n_freq*2, n_freq*2)
  covar_ff_cut = dblarr(n_freq*2, n_freq*2)
  
  fisher_diag = dblarr(n_freq*2)
  fisher_diag[2*indgen(n_freq)] = f_cos2
  fisher_diag[2*indgen(n_freq)+1] = f_sin2
  fisher_off_diag = dblarr(n_freq*2-1)
  fisher_off_diag[2*indgen(n_freq)] = f_cross
  fisher_ff_cut = diag_matrix(fisher_diag) + diag_matrix(fisher_off_diag, -1) + diag_matrix(fisher_off_diag, 1)

  for i = 0, n_freq-1 do covar_ff_cut[2*i:2*i+1,2*i:2*i+1] = invert(fisher_ff_cut[2*i:2*i+1,2*i:2*i+1])

  ;; ;; covar_ff2 = invert(fisher_ff_cut, status)
  ;; ;; if status ne 0 then begin
  ;; ;;    print, 'Inversion was not successful, using pseudo inverse'
  ;; ;;    covar_ff2 = idl_pseudo_inverse(fisher_ff_cut)
  ;; ;; endif

  ;; ft_matrix = dblarr(npts, n_freq*2)
  ;; ft_matrix[*, indgen(n_freq)*2] = cos_arr
  ;; ft_matrix[*, indgen(n_freq)*2+1] = sin_arr
  ;; ft_matrix = ft_matrix / sqrt(n_freq)

  ;; fisher_tt = diag_matrix(1/sigma2)
  ;; fisher_ff = ft_matrix ## fisher_tt ## transpose(ft_matrix)
  ;; ;; covar_ff3 = invert(fisher_ff, status, /double)
  ;; ;; if status ne 0 then begin
  ;; ;;    print, 'Inversion was not successful, using pseudo inverse'
  ;; ;;    covar_ff3 = idl_pseudo_inverse(fisher_ff)
  ;; ;; endif

  ;; covar_tt = diag_matrix(sigma2)
  ;; covar_ff = ft_matrix ## covar_tt ## transpose(ft_matrix)
  ;; covar_ff[0,0] = covar_ff[0,0]/4d

  ;; ;;temp = ft_matrix ## transpose(ft_matrix)
  ;; if keyword_set(use_fisher_cut) then covar_use = covar_ff_cut else covar_use = covar_ff

  ;; covar_diag = diag_matrix(covar_use)
  ;; covar_cos1 = covar_diag[2*indgen(n_freq)]
  ;; covar_sin1 = covar_diag[2*indgen(n_freq)+1]

  ;; covar_off_diag = diag_matrix(covar_use, 1)
  ;; covar_cross1 = covar_off_diag[2*indgen(n_freq)]

  sigma2_arr = rebin(sigma2, npts, n_freq)
  covar_cos = total(sigma2_arr*cos_arr^2d, 1)/double(n_freq)
  covar_sin = total(sigma2_arr*sin_arr^2d, 1)/double(n_freq)
  covar_cross = total(sigma2_arr*cos_arr*sin_arr, 1)/double(n_freq)

  covar_cos[0] = covar_cos[0]/4d


  theta = atan(2*covar_cross, covar_cos - covar_sin)/2d
  cos_theta = cos(theta)
  sin_theta = sin(theta)

  sigma1_2 = covar_cos*cos_theta^2 + 2d*covar_cross*cos_theta*sin_theta + covar_sin*sin_theta^2d
  sigma2_2 = covar_cos*sin_theta^2d - 2d*covar_cross*cos_theta*sin_theta + covar_sin*cos_theta^2d

  print, 'theta (deg), sigma1^2, sigma2^2'
  print, theta[freq_ind]*180/!pi, sigma1_2[freq_ind], sigma2_2[freq_ind]

  data_cos = real_part(data_ft)
  data_sin = imaginary(data_ft)
  data1 = data_cos*cos_theta + data_sin*sin_theta
  data2 = (-1d)*data_cos*sin_theta + data_sin*cos_theta

  if keyword_set(std_power) then begin
     ;; calculate the error on data1^2 & data2^2 and propagate to get error on power
     error1 = 2d * abs(data1) * sigma1_2
     error2 = 2d * abs(data2) * sigma2_2
     
     power = data1^2d + data2^2d
     power_error = sqrt(error1^2d + error2^2d)
  endif else begin
     ;; calculate power as weighted sum of data1^2 & data2^2
     
     weights_1 = (1/sigma1_2)^2d
     term1 = data1^2d*weights_1
     wh = where(sigma1_2 eq 0, count)
     if count ne 0 then begin
        weights_1[wh] = 0
        term1[wh] = 0
     endif
     
     weights_2 = (1/sigma2_2)^2d
     term2 = data2^2d*weights_2
     wh = where(sigma2_2 eq 0, count)
     if count ne 0 then begin
        weights_2[wh] = 0
        term2[wh] = 0
     endif
     
     weights = weights_1 + weights_2
     power = (term1 + term2) / weights
     power_error = 1/weights
     wh = where(weights eq 0, count)
     if count ne 0 then begin
        power[wh] = 0
        power_error[wh] = 0
     endif

  endelse
  
  a0 = 1
  b0 = 0
  freq = frequencies[freq_ind]
  model = model(a0, b0, freq, times)
  envelope1 = model - sigma2
  envelope2 = model + sigma2
  xxx=[times,reverse(times)]
  yyy=[envelope1,reverse(envelope2)]

  yrange = [-2.5, 2.5]
  xrange = minmax(times)

  if windowavailable(1) then wset, 1 else window, 1, xsize = 700, ysize = 700
  !p.multi = [0, 1, 2]
  plot, times, data, psym = 3, yrange = yrange, xrange = xrange, ystyle = 1, xtitle = 'times', ytitle = 'amplitude', $
        title = 'Data with model (red) and 1 sigma error region (blue)'
  polyfill,xxx,yyy,color=75
  oplot, times, data, color =0, psym = 1
  oplot, times, model, color = 254

  yrange = [-1.5, 1.5]
  xxx=[times,reverse(times)]
  yyy=[-1*sigma2,reverse(sigma2)]
  plot, times, data-model, psym = 3, yrange = yrange, xrange = xrange, ystyle = 1, xtitle = 'times', ytitle = 'residual amplitude', $
        title = 'Residuals with 1 sigma error region (blue)'
  polyfill,xxx,yyy,color=75
  oplot, times, data-model, color =0, psym = 1
  oplot, times, model-model, color = 254

  !p.multi = 0

  model_power = abs((fft(model))[0:n_freq-1])^2d
  ;; yrange = minmax([power, power_error, model_power])
  yrange = [0,0.3]
  xrange = minmax(frequencies)
  if windowavailable(2) then wset, 2 else window, 2, xsize = 700, ysize = 350
  plot, frequencies, power, yrange = yrange, xrange = xrange, xstyle = 1, xtitle = 'frequencies', $
           ytitle = 'power', title = 'Power with power errors (blue)', psym=4
  oplot, frequencies, model_power, color = 254,psym=4
  oploterror, frequencies, power, power_error,/nohat, errcolor = 75, psym=4


  if windowavailable(3) then wset, 3 else window, 3, xsize = 700, ysize = 350
  binsize = 1e-5
  power_hist = histogram(power, locations = locs, binsize = binsize, omin=min_val, omax = max_val)
  power_err_hist = histogram(power_error, binsize = binsize, min=min_val, max = max_val)

  x_range = [0, max([locs[where(power_err_hist eq max(power_err_hist))] + binsize/2d, 5e-4])]
  plot, locs+binsize/2., power_hist, xrange = x_range, title = 'Power and power_error (blue) histogram'
  oplot, locs+binsize/2., power_err_hist, color = 75

  hessian = second_deriv_matrix(a0, b0, freq,  times, data, sigma2)
  covar_chi2 = invert (0.5 * hessian)
  
  print, 'Covariance based on chi2'
  print, covar_chi2

  print, 'Covariance based on rotation'
  print, covar_rot

end
