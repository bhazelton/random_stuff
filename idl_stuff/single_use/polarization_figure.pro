pro polarization_figure, pub = pub, grey_scale = grey_scale, hinv = hinv

  if n_elements(hinv) eq 0 then hinv = 1

  froot = base_path('data') + 'fhd_simulations_old/'
  array = '496t'

  info_file = froot + 'sim_496t_info.idlsave'
  restore, info_file   

  sim_num = 0
  names = string([25, 50, 100, 200, 400, 600, 800, 1000], format = '(i04)')
  ;;names = 'offaxis'
  ;;names = 'fullsky'

  filebase = 'sim_496t_' + names[sim_num] + '_uvf'
  initial_uv_file = froot + filebase + '_initial.idlsave'

  truesky_uv_beamcorr_file = froot + filebase + '_truesky_uv_plane_beamcorr.idlsave'
  truesky_uv_file = froot + filebase + '_truesky_uv_plane.idlsave'

  holo_uv_file = froot + filebase + '_holo_uv_plane.idlsave'

  weights_uv_file = froot + 'sim_496t_weights_uv_plane.idlsave'

  if keyword_set(hinv) then xrange = [.17,.31] else xrange = [.12,.22]
  yrange = xrange
  

  plotfile = base_path('plots') + 'single_use/polarization_uv_figure.eps'

  window_num = 1
  start_multi_params = {ncol:2, nrow:2, ordering:'row'}


  uvf_slice_plot, initial_uv_file, multi_pos = positions, start_multi_params = start_multi_params, plot_xrange = xrange, $
                  plot_yrange = yrange, title = 'Initial', /baseline_axis, grey_scale = grey_scale, $
                  plotfile = plotfile, window_num = window_num, pub = pub, hinv = hinv

  uvf_slice_plot, holo_uv_file, multi_pos = positions[*,1], plot_xrange = xrange, plot_yrange = yrange, $
                  title = 'Gridded', /baseline_axis, /color_0amp, grey_scale = grey_scale, pub = pub, hinv = hinv
  uvf_slice_plot, truesky_uv_beamcorr_file, multi_pos = positions[*,2], plot_xrange = xrange, plot_yrange = yrange, $
                  title = 'Beam corrected!C(before 0 weight areas removed)', /baseline_axis, /color_0amp, $
                  grey_scale = grey_scale, pub = pub, hinv = hinv
  uvf_slice_plot, truesky_uv_file, multi_pos = positions[*,3], plot_xrange = xrange, plot_yrange = yrange, $
                  title = 'Beam corrected', /baseline_axis, /color_0amp, grey_scale = grey_scale, pub = pub, hinv = hinv
 
  if keyword_set(pub) then begin
     psoff
     wdelete, window_num
  endif
  undefine, positions

  ;; make some histograms of differences from input uv toward center of uv plane.
  
  restore, initial_uv_file
  initial_uv = uvf_slice

  restore, holo_uv_file
  holo_uv = uvf_slice

  restore, truesky_uv_beamcorr_file
  truesky_nowt_uv = uvf_slice

  restore, truesky_uv_file
  truesky_uv = uvf_slice

  restore, weights_uv_file
  weights_uv = uvf_slice

  holo_diff = !dpi - abs(abs(atan(holo_uv,/phase) - atan(initial_uv, /phase)) - !dpi)
  truesky_nowt_diff = !dpi - abs(abs(atan(truesky_nowt_uv,/phase) - atan(initial_uv, /phase)) - !dpi)
  truesky_diff = !dpi - abs(abs(atan(truesky_uv,/phase) - atan(initial_uv, /phase)) - !dpi)

  radius = xarr[1024+64]
  pix_radii = sqrt(rebin(xarr, n_elements(xarr), n_elements(yarr))^2 + $
                   rebin(reform(yarr, 1, n_elements(yarr)), n_elements(xarr), n_elements(yarr))^2)
  wh_in = where(pix_radii le radius, count_in)
  if count_in eq 0 then message, 'radius excludes all pixels'
  npix = float(count_in)
  npix_total = n_elements(xarr) * n_elements(yarr)

  binsize = 0.001
  hist_max = ceil(max(truesky_diff[wh_in])/binsize) * binsize
  holo_hist = histogram(holo_diff[wh_in], binsize=binsize, min = 0, max = hist_max, locations = hist_locs)
  truesky_hist = histogram(truesky_diff[wh_in], binsize=binsize, min = 0, max = hist_max)
  
  if keyword_set(pub) then begin
     charthick = 3
     thick = 3
     xthick = 3
     ythick = 3
     charsize = 2
     font = 1
     
  endif else begin
     charthick = 1
     thick = 1
     xthick = 1
     ythick = 1
     font = -1
     charsize = 1
   endelse

  plotfile = base_path('plots') + 'single_use/polarization_hist_figure.eps'
  window_num=2
  start_multi_params = {ncol:2, nrow:1, ordering:'row'}

  plot_range = 1024+[-1,1]*256
  if keyword_set(hinv) then plot_xrange = xarr[plot_range] / hubble_param else plot_xrange = xarr[plot_range]
  if keyword_set(hinv) then plot_yrange = yarr[plot_range] / hubble_param else plot_yrange = yarr[plot_range]
  uvf_slice_plot, weights_uv_file, multi_pos = positions, start_multi_params = start_multi_params, type='weights', $
                  plot_xrange = plot_xrange,  plot_yrange = plot_yrange, /no_title, /baseline_axis, $
                  grey_scale = grey_scale, plotfile = plotfile, window_num = window_num, pub = pub, hinv = hinv
  tvcircle, radius, 0, 0, color = 'black', thick = thick
  cgplots, [xrange, reverse(xrange), xrange[0]], [yrange[0], yrange, reverse(yrange)], color = 'black', thick = thick/2


  pos_len = [positions[2,1] - positions[0,1], positions[3,1] - positions[1,1]]
  pos_center = [pos_len[0]/2d + positions[0,1], pos_len[1]/2d + positions[1,1]]
  margin = [0.1, 0.08, 0.01, 0.07]
  
  new_len = pos_len * [1 - margin[2] - margin[0], 1 - margin[3] - margin[1]]
  new_pos = [pos_center[0] - new_len[0]/2. + margin[0]*pos_len[0], $
             pos_center[1] - new_len[1]/2. + margin[1]*pos_len[1], $
             pos_center[0] + new_len[0]/2. - margin[2]*pos_len[0], $
             pos_center[1] + new_len[1]/2. - margin[3]*pos_len[1]]

  cgplot, hist_locs*180d/!dpi, holo_hist/npix, position = new_pos, psym=10, color = 'black', xstyle=1, ytitle = 'fraction of pixels', $
          xtitle = 'Phase difference from initial uv (degrees)', /noerase, thick = thick, charthick = charthick, $
          xthick = xthick, ythick = ythick, charsize = charsize, font = font
  cgplot, /overplot, hist_locs*180d/!dpi, truesky_hist/npix, psym=10, color = 'red', thick = thick
  al_legend, ['Gridded', 'Beam corrected'], textcolor = ['black', 'red'], box = 0, /right, charthick = charthick, $
             charsize = charsize, font = font


  if keyword_set(pub) then begin
     psoff
     wdelete, window_num
  endif


end
