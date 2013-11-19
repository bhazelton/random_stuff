pro plot_gain_track, gain_track, tile, gain_fit_mode=gain_fit_mode

  track_dims = size(gain_track, /dimension)
  case n_elements(track_dims) of
    3: gain_track_use = reform(gain_track[*,tile,*])
    2: gain_track_use = gain_track
    else: stop
  endcase
  
  n_mode = track_dims[0]
  n_iter = n_elements(gain_track[0,*])
  
  fit_dims = size(gain_fit_mode, /dimension)
  if fit_dims[0] gt 0 then begin
    case n_elements(fit_dims) of
      2: gain_fit_use = reform(gain_fit_mode[*, tile])
      1: gain_fit_use = gain_fit_mode
      else: stop
    endcase
  endif
  
  multi_setting = !p.multi
  !p.multi=[0, 1, n_mode]
  for mode_i=0, n_mode-1 do begin
  if n_elements(gain_fit_use) gt 0 then yrange = minmax([reform(gain_track_use[mode_i, *]), gain_fit_use[mode_i], 0.1, -0.1]) $
  else yrange = minmax([reform(gain_track_use[mode_i, *]), 0.1, -0.1])
    cgplot, gain_track_use[mode_i, *], yrange=yrange;*[.9,1.1]
    if n_elements(gain_fit_use) gt 0 then cgplot, [0, n_iter], gain_fit_use[mode_i]*[1,1],/over, color='blue'
  end
  !p.multi=multi_setting
end
