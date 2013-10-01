pro fhd_caltest

  filename = base_path('data')+'single_use/fhd_caltest.idlsave'
  restore, filename
  
  
  dist_hist_avg = histogram(dist_arr_avg, min=0, binsize=2, locations = locs, reverse_indices = ri_avg)
  dist_hist_full = histogram(dist_arr_full, min=0, binsize=2, reverse_indices = ri_full)
  
  mean_amp_avg = fltarr(n_elements(locs))
  max_amp_avg = fltarr(n_elements(locs))
  var_amp_avg = fltarr(n_elements(locs))
  mean_amp_model_avg = fltarr(n_elements(locs))
  max_amp_model_avg = fltarr(n_elements(locs))
  var_amp_model_avg = fltarr(n_elements(locs))
  
  mean_amp_full = fltarr(n_elements(locs))
  max_amp_full = fltarr(n_elements(locs))
  var_amp_full = fltarr(n_elements(locs))
  mean_amp_model_full = fltarr(n_elements(locs))
  max_amp_model_full = fltarr(n_elements(locs))
  var_amp_model_full = fltarr(n_elements(locs))
  for i=0, n_elements(locs)-1 do begin
    if dist_hist_avg[i] gt 0 then begin
      ave_set = ri_avg[ri_avg[i]:ri_avg[i+1]-1]
      ave_data = vis_avg[ave_set]
      ave_model = model_avg[ave_set]
      mean_amp_avg[i] = mean(abs(ave_data))
      max_amp_avg[i] = max(abs(ave_data))
      var_amp_avg[i] = variance(abs(ave_data))
      mean_amp_model_avg[i] = mean(abs(ave_model))
      max_amp_model_avg[i] = max(abs(ave_model))
      var_amp_model_avg[i] = variance(abs(ave_model))
    endif
    
    if dist_hist_full[i] gt 0 then begin
      full_set = ri_full[ri_full[i]:ri_full[i+1]-1]
      full_data = vis_full[full_set]
      full_model = model_full[full_set]
      mean_amp_full[i] = mean(abs(full_data))
      max_amp_full[i] = max(abs(full_data))
      var_amp_full[i] = variance(abs(full_data))
      mean_amp_model_full[i] = mean(abs(full_model))
      max_amp_model_full[i] = max(abs(full_model))
      var_amp_model_full[i] = variance(abs(full_model))
    endif
  endfor
  
  cgplot, locs, mean_amp_full, /ylog, /xlog, xrange = [1, 2e3], yrange = [1e-1, 1e6], psym=10
  cgplot, locs, max_amp_full, /over, color = 'dark grey', psym=10
  cgplot, locs, sqrt(var_amp_full), /over, color = 'medium grey', psym=10
  cgplot, locs, mean_amp_avg, /over, color = 'red', psym=10
  cgplot, locs, max_amp_avg, /over, color = 'salmon', psym=10
  cgplot, locs, sqrt(var_amp_avg), /over, color = 'pink', psym=10
  al_legend, ['avg mean', 'avg sigma', 'avg max', 'full mean', 'full sigma', 'full max'], $
    textcolors = ['red', 'pink', 'salmon', 'black', 'medium grey', 'dark grey'], /right, box=0
    
  window, /free
  cgplot, locs, mean_amp_model_full, /ylog, /xlog, yrange = [1, 1e4], xrange = [1, 2e3], psym=10
  cgplot, locs, max_amp_model_full, /over, color = 'dark grey', psym=10
  cgplot, locs, sqrt(var_amp_model_full), /over, color = 'medium grey', psym=10
  cgplot, locs, mean_amp_model_avg, /over, color = 'red', psym=10
  cgplot, locs, max_amp_model_avg, /over, color = 'salmon', psym=10
  cgplot, locs, sqrt(var_amp_model_avg), /over, color = 'pink', psym=10
  al_legend, ['avg mean', 'avg sigma', 'avg max', 'full mean', 'full sigma', 'full max'], $
    textcolors = ['red', 'pink', 'salmon', 'black', 'medium grey', 'dark grey'], /right, box=0
    
  window, /free
  cgplot, locs, mean_amp_full, xrange = [0, 200], yrange = [0, 200], psym=10
  cgplot, locs, mean_amp_model_full*5, /over, color='blue', psym=10
  
end
