pro test_hpx_weights, pol, log=log, variance=variance, sum_freq = sum_freq, all_freq = all_freq

  datafile = base_path('data') + 'fhd_ps_data/Combined_obs_EOR1_P00_145_20110926193959-EOR1_P00_145_20110926200503_even_cube.sav' 
  ;;datafile = base_path('data') + 'fhd_ps_data/Combined_obs_EOR1_P00_145_20110926193959-EOR1_P00_145_20110926200503_odd_cube.sav' 
  ;;datafile = base_path('data') + 'fhd_ps_data/multi_freq_residuals_cube_healpix.sav'

  if n_elements(pol) eq 0 then pol = 'xx'
  wh_pol = where(['xx', 'yy'] eq pol, count_pol)
  if count_pol eq 0 then message, 'pol not recognized'

  pixel_nums = getvar_savefile(datafile, 'hpx_inds')
  if keyword_set(variance) then hpx_weight = getvar_savefile(datafile, 'variance_' + pol + '_cube') $
  else hpx_weight = getvar_savefile(datafile, 'weights_' + pol + '_cube')
  nside = getvar_savefile(datafile, 'nside')
  
  ;; Angular resolution is given in Healpix paper in units of arcminutes, need to convert to radians
  ang_resolution = sqrt(3d/!pi) * 3600d/nside * (1d/60d) * (!pi/180d)
  degpix = ang_resolution * 180d / !dpi

  pix2ang_ring, nside, pixel_nums, pix_center_theta, pix_center_phi

  theta_deg = pix_center_theta * 180. / !pi
  phi_deg = pix_center_phi * 180. / !pi
  xrange = [floor(min(theta_deg)/10.)*10, ceil(max(theta_deg)/10.)*10]
  yrange = [floor(min(phi_deg)/10.)*10, ceil(max(phi_deg)/10.)*10]
  ;;resolution = intarr(2) + 1000
  resolution = [xrange[1]-xrange[0], yrange[1]-yrange[0]]*1.5/degpix
  
  temp = total(abs(hpx_weight), 1)
  wh_freq0 = where(temp eq 0, count_freq0, complement=wh_freq_n0, ncomplement = count_freq_n0)
  if count_freq_n0 eq 0 then message, 'weights are all 0'
  if count_freq0 gt 0 then print, number_formatter(count_freq0) + ' of ' + number_formatter(n_elements(temp)) + $
                                  ' frequencies have all 0 weights'

  if keyword_set(all_freq) then begin
     data = hpx_weight[*,wh_freq_n0]
     nsteps = count_freq_n0
     freq_str = 'Channel ' + number_formatter(wh_freq_n0)
  endif else if keyword_set(sum_freq) then begin
     data = total(hpx_weight[*,wh_freq_n0], 2)
     freq_str = 'All non-zero Channels'
  endif else begin
     freq_ind = wh_freq_n0[0]
     data = hpx_weight[*,freq_ind]
     freq_str = 'Channel ' + number_formatter(freq_ind)
  endelse

  data_range = minmax(data)

  tvlct, r, g, b, /get

  round_log = floor(alog10(max(abs(data_range))))-1
  plot_range = [floor(data_range[0]/10.^round_log)*10.^round_log, ceil(data_range[1]/10.^round_log)*10.^round_log]

  nlevels = round((plot_range[1]-plot_range[0])/10.^(round_log-1) / 10.)
  while nlevels gt 10 or nlevels lt 5 do begin
     if nlevels gt 10 then nlevels = nlevels / 2.
     if nlevels lt 5 then nlevels = nlevels * 2.
  endwhile

  cgloadct, 25, /brewer, /reverse, ncolors=nlevels, bottom=1
  levels = indgen(nlevels) * (plot_range[1]-plot_range[0]) / (nlevels-1) + plot_range[0]

  window_num = 1
  if keyword_set(all_freq) then begin
 
     ncol = round(sqrt(nsteps) * 1.5)

     data_aspect = (yrange[1]-yrange[0]) / float((xrange[1]-xrange[0]))
     position_struct = calculate_plot_positions(data_aspect, nplots = nsteps, /no_colorbar, nrow = nrow)

     if windowavailable(window_num) then begin 
        wset, window_num
        if !d.x_size ne position_struct.wsize[0] or !d.y_size ne position_struct.wsize[1] then make_win = 1 else make_win = 0
     endif else make_win = 1
     if make_win eq 1 then window, window_num, xsize = position_struct.wsize[0], ysize = position_struct.wsize[1]

     for i=0, nsteps-1 do begin
        cgcontour, data[*,i], /noerase, theta_deg, phi_deg, position = position_struct.plot_pos[i,*], /irregular, $
                   resolution=resolution, levels = levels, c_colors = indgen(nlevels)+1, /fill, xrange = xrange, yrange = yrange, $
                   xtitle = 'theta (degrees)', ytitle = 'phi (degrees)', title = 'Weights ' + strupcase(pol) + ' ' + freq_str[i]

        ;;cgcolorbar, range=plot_range, position = position_struct.cb_pos[i,*], divisions= nlevels, ncolors = nlevels, $
        ;;            bottom=1, /discrete,/vertical

     endfor
  endif else begin
     cgcontour, data, theta_deg, phi_deg, position = [.1,.1,.75,.9], /irregular, resolution=resolution, levels = levels, $
                c_colors = indgen(nlevels)+1, /fill, xrange = xrange, yrange = yrange, $
                xtitle = 'theta (degrees)', ytitle = 'phi (degrees)', title = 'Weights ' + strupcase(pol) + ' ' + freq_str

  ;;contour, data, theta_deg, phi_deg, position = [.1,.1,.75,.9], /irregular, levels = levels, $
  ;;         c_colors = indgen(nlevels)+1, /fill, xrange = xrange, yrange = yrange, $
  ;;         xtitle = 'theta (degrees)', ytitle = 'phi (degrees)', title = 'Weights ' + strupcase(pol)

  ;;contour, data, theta_deg, phi_deg, /overplot, /irregular, levels = levels, c_colors = 0

     cgcolorbar, range=plot_range, position = [.87, .1,.9,.9], divisions= nlevels, ncolors = nlevels, $
                 bottom=1, /discrete,/vertical

  endelse

  tvlct, r, g, b

end
