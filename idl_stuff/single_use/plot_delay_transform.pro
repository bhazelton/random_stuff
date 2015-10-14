pro plot_delay_transform, save_file, window_num = window_num

  filebase = cgrootname(save_file, directory = dir, extension = exten)
  flag_file = dir + (strsplit(filebase, 'vis', /extract))[0] + 'flags.sav'
  
  metadata_dir = (strsplit(dir, 'vis_data',/extract, /regex))[0] + 'metadata/'
  params_file = metadata_dir + (strsplit(filebase, 'vis', /extract))[0] + 'params.sav'
  
  restore, save_file
  if n_elements(vis_ptr) gt 0 then vis = *vis_ptr else if n_elements(vis_model_ptr) then vis = *vis_model_ptr else $
    message, 'save_file does not contain vis_ptr or vis_model_ptr'
  undefine_fhd, vis_ptr
  
  freq = (*obs.baseline_info).freq
  n_freq = n_elements(freq)
  delta_freq = freq[1] - freq[0]
  
  delta_delay = 1./(n_freq * delta_freq)
  delays = (findgen(n_freq)-n_freq/2)*delta_freq
  
  
  restore, params_file
  ;; convert u,v,w from light travel time in seconds to meters
  uu = params.uu * 3e8
  vv = params.vv * 3e8
  ww = params.ww * 3e8
  undefine_fhd, params
  
  baseline_length = sqrt(uu^2. + vv^2.)
  wh_cross = where(baseline_length gt 0, count_cross)
  if count_cross eq 0 then stop
  vis = vis[*, wh_cross]
  baseline_length = baseline_length[wh_cross]
  
  vis_ft = fft(vis, dimension=1)
  vis_ft = shift(vis_ft, n_freq/2)
  
  len_sort = sort(baseline_length)
  baseline_length = baseline_length[len_sort]
  horizon_delay = baseline_length / 3e8
  horizon_delay_bin_upper = (horizon_delay / delta_delay) + n_freq/2
  horizon_delay_bin_lower = (horizon_delay / delta_delay)*(-1) + n_freq/2
  
  vis_ft_sort = temporary(vis_ft[*, len_sort])
  vis_sort = temporary(vis[*, len_sort])
  
  
  quick_image, transpose(abs(vis_sort)), window_num = window_num+1, title = filebase, ytitle = 'freq channel', xtitle = 'baseline number (ordered by length)'
  
  quick_image, transpose(abs(vis_ft_sort)), window_num = window_num, title = filebase, ytitle = 'delay channel', xtitle = 'baseline number (ordered by length)'
  cgplot, horizon_delay_bin_upper, /over
  cgplot, horizon_delay_bin_lower, /over
    
end