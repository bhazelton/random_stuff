

pro multibaseline_figures, refresh = refresh, pub = pub, grey_scale = grey_scale

  savefile = base_path('data') + 'single_use/multibaseline_data.idlsave'
  ftest = file_test(savefile) *  (1 - file_test(savefile, /zero_length))

  if ftest eq 0 or keyword_set(refresh) then begin
     froot = base_path('data') + 'fhd_simulations_old/'
     simfile = froot + 'sim_496t_offaxis_uvf.idlsave'
     info_file = froot + 'sim_496t_info.idlsave'
     weights_file = froot + 'sim_496t_weights.idlsave'
     beam_file = froot + 'sim_496t_beam.idlsave'
     
     tile_locs_file = '496T_tile_locations.txt'
     
     textfast,tile_locs,/read,filename=tile_locs_file,root=rootdir('mwa'),first_line=1,column=indgen(2)+1
     restore, info_file
     restore, simfile
     restore, weights_file
     undefine, tile_beam_image
     
     restore, beam_file
     beam_uv = beam_uv[*,*,*,0] ;; just go to XX pol for now
     
     sim_cube = temporary(uvf_cube) / weights_cube
     wh_wt0 = where(weights_cube eq 0, count_wt0)
     if count_wt0 gt 0 then sim_cube[wh_wt0] = 0
     
     n_tiles = 496L
     
     dims = size(sim_cube, /dimension)
     n_u = dims[0]
     n_v = dims[1]
     n_freq = n_elements(frequencies)
     
     uv_binsize = 1/ ((abs(degpix) * !pi / 180d) * n_u)
     if n_u eq n_v then uv_arr = (dindgen(n_u) - n_u/2d) * uv_binsize else stop
     
     xtile=reform(tile_locs[0,*])
     ytile=reform(tile_locs[1,*])
     
     tile_A = reform(rebin(indgen(n_tiles)+1, n_tiles, n_tiles), n_tiles*n_tiles)
     tile_B = reform(rebin(reform(indgen(n_tiles)+1, 1, n_tiles), n_tiles, n_tiles), n_tiles*n_tiles)
     
     ;; get unique baselines
     temp = double(tile_A)^3/double(tile_B) + double(tile_B)^3/double(tile_A)
     uniq_inds = uniq(temp, sort(temp))
     nbaselines = n_elements(uniq_inds)
     
     if nbaselines ne ((n_tiles^2-n_tiles)/2 + n_tiles) then stop
     
     tile_A = tile_A[uniq_inds]
     tile_B = tile_B[uniq_inds]
     
     c_light=299792458d
     baseline_u = rebin((xtile[tile_A-1]-xtile[tile_B-1])/c_light, nbaselines, n_freq) * $
                  rebin(reform(frequencies * 1e6, 1, n_freq), nbaselines, n_freq)
     baseline_v = rebin((ytile[tile_A-1]-ytile[tile_B-1])/c_light, nbaselines, n_freq) * $
                  rebin(reform(frequencies * 1e6, 1, n_freq), nbaselines, n_freq)
     
     ;;beam_radii = (dblarr(n_freq) + 10d/c_light) * frequencies * 1e6
     
     ;; make beam shape based on contour of beam
     n_uv_beam =  psf_resolution*psf_dim
     beam_uv_binsize = uv_binsize/psf_resolution
     beam_uv_arr  = (dindgen(n_uv_beam) - n_uv_beam/2d) * beam_uv_binsize
     ;; normalize to make integral = 1
     beam_uv = beam_uv / rebin(reform(total(total(beam_uv, 2), 1), 1, 1, n_freq), n_uv_beam, n_uv_beam, n_freq)
     
     ;;beam_cutoff_factor = 0.01
     beam_radii = dblarr(n_freq)
     window, 3
     for i = 0, n_freq-1 do begin
        level = max(beam_uv[*,*,i])*0.002d

        ;;plot, beam_uv_arr, beam_uv_arr, /xstyle,/ystyle
        ;;cgimage, beam_uv[*,*,i], /noerase, /overplot,/scale
        ;;cgcontour, beam_uv[*,*,i], beam_uv_arr, beam_uv_arr, /overplot, levels=level, c_colors=255

        cgcontour, beam_uv[*,*,i], beam_uv_arr, beam_uv_arr, levels=level, $
                   path_info = beam_contour_info, path_xy = beam_contour_xy, /path_data_coords
        beam_radii[i] = max(sqrt(beam_contour_xy[0,*]^2d + beam_contour_xy[1,*]^2d))

        plot, beam_uv_arr, beam_uv_arr, /xstyle,/ystyle
        cgimage, alog10(beam_uv[*,*,i]), /noerase, /overplot,/scale
        oplot, beam_contour_xy[0,*], beam_contour_xy[1,*], color=255
             
        tag = 'freq_' + number_formatter(i)
        if i eq 0 then beam_shape = create_struct(tag, beam_contour_xy) else $
           beam_shape = create_struct(beam_shape, tag, beam_contour_xy)
     endfor
     wdelete, 3
    
     ;; get input cube
     ;; columns are: x location (zero centered), y location (zero centered), flux
     my_sources = reform(my_sources, 3, n_elements(my_sources)/3)
     n_sources = (size(my_sources, /dimension))[1]
     my_source_array = reform(fltarr(7, n_sources), 7, n_sources)
     my_source_array[0,*]=my_sources[0,*] + n_u/2.
     my_source_array[1,*]=my_sources[1,*] + n_v/2.
     my_source_array[4,*]=my_sources[2,*] + n_v/2.
     my_source_array[5,*]=Round(my_source_array[0,*])+Round(my_source_array[1,*])*n_u
     
     model_uv = visibility_source_uv_grid(my_source_array, u_dim = n_u, v_dim = n_v)
     
     ;; get gridding normalizations
     norm_sources = reform([0,0,1], 3, 1)
     n_sources = (size(norm_sources, /dimension))[1]
     norm_source_array = reform(fltarr(7, n_sources), 7, n_sources)
     norm_source_array[0,*]=norm_sources[0,*] + n_u/2.
     norm_source_array[1,*]=norm_sources[1,*] + n_v/2.
     norm_source_array[4,*]=norm_sources[2,*] + n_v/2.
     norm_source_array[5,*]=Round(norm_source_array[0,*])+Round(norm_source_array[1,*])*n_u
     test_gridding = visibility_source_uv_grid(norm_source_array)
     norm_gridding = 1/max(abs(test_gridding))
     
     model_uv = model_uv * norm_gridding

     save, file = savefile, sim_cube, weights_cube, model_uv, uv_arr, baseline_u, baseline_v, frequencies, beam_uv, beam_shape, $
           beam_radii, nbaselines, beam_uv_binsize, beam_uv_arr, psf_resolution, n_uv_beam, n_freq, n_u, n_v, uv_binsize
  endif else restore, savefile


  ;; Beam figure
  if keyword_set(grey_scale) then plotfile = base_path('plots') + 'single_use/multibaseline_beam_grey.eps' $
  else plotfile = base_path('plots') + 'single_use/multibaseline_beam.eps'

  ;;fractions = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, .1, .2, .5]
  ;;fractions = [0.002, 0.1*(indgen(9)+1)]

  beam_freq_ind = 16

  ;; set up beam uv plot
  fractions = [0.0002, 0.002, 0.02, .2, .5]
  levels = max(beam_uv[*,*,0])*fractions
  colors = strarr(n_elements(levels)) + 'darkgrey'
  colors[where(fractions eq 0.002)] = 'black'
  if keyword_set(grey_scale) then colors[n_elements(levels)-1] = 'grey'
  annotations = number_formatter(fractions*100) + ' %'
  linestyles = intarr(n_elements(levels))+2
  linestyles[where(fractions eq 0.002)] = 0
  linestyles[where(fractions eq 0.5)] = 1

  range = [-1,1]*round(beam_radii[beam_freq_ind])
  ind_range = where(beam_uv_arr gt range[0] and beam_uv_arr lt range[1], count_inds)
  if count_inds lt 1 then stop
  
  beam_use = beam_uv[ind_range,ind_range,beam_freq_ind]
  uv_arr_use = beam_uv_arr[ind_range]

  beam_norm = beam_use/max(beam_use) * 255

  ;; now work on beam theta plot
  beam_theta = abs(fft_shift(fft(beam_uv[*,*,beam_freq_ind])) * n_uv_beam^2. * beam_uv_binsize^2.)
  beam_degpix = (1/(n_uv_beam * beam_uv_binsize))*180/!pi
  beam_deg_arr = (dindgen(n_uv_beam) - n_uv_beam/2d) * beam_degpix

  deg_range = [-1,1] * 90
  deg_ind_range = where(beam_deg_arr gt deg_range[0] and beam_deg_arr lt deg_range[1], count_deg_inds)
  if count_deg_inds lt 1 then stop

  beam_deg_use = beam_theta[deg_ind_range, deg_ind_range, 0]
  deg_arr_use = beam_deg_arr[deg_ind_range]

  beam_deg_norm = beam_deg_use/max(beam_deg_use)*255

  cb_size = 0.025
  margin = [0.12, 0.14, 0.02, 0.04]
  cb_margin = [0.14, 0.02] 
  plot_pos = [margin[0], margin[1], (1-cb_margin[1]-cb_size-cb_margin[0]-margin[2]), (1-margin[3])]
  cb_pos = [(1-cb_margin[1]-cb_size), margin[1], (1-cb_margin[1]), (1-margin[3])]
  
  plot_len = [plot_pos[2]-plot_pos[0], plot_pos[3] - plot_pos[1]]
  plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])

  data_aspect=1
  aspect_ratio =  data_aspect /plot_aspect
  if aspect_ratio gt 1 then begin
     y_factor = aspect_ratio
     x_factor = 1.
  endif else begin   
     y_factor = 1.
     x_factor = 1./aspect_ratio
  endelse
  
  ncol = 1
  nrow = 2
  col_val = [0,0]
  row_val = [0,1]

  ;; taken from kpower_2d_plots
  multi_pos = fltarr(4,2)
  multi_pos[0,*] = col_val/double(ncol)
  multi_pos[1,*] = row_val/double(nrow)
  multi_pos[2,*] = (col_val+1)/double(ncol)
  multi_pos[3,*] = (row_val+1)/double(nrow)

  max_ysize = 1000
  max_xsize = 1200
  base_size = 600

  ;; define window size based on aspect ratio
  base_size_use = base_size
  xsize = round(base_size * x_factor * ncol)
  ysize = round(base_size * y_factor * nrow)
  while (ysize gt max_ysize) or (xsize gt max_xsize) do begin
     base_size_use = base_size_use - 100
     xsize = round(base_size_use * x_factor * ncol)
     ysize = round(base_size_use * y_factor * nrow)
  endwhile

  new_pos = fltarr(4,2)
  new_cb_pos = fltarr(4,2)
  for i=0, 1 do begin
     multi_xlen = (multi_pos[2,i]-multi_pos[0,i])
     multi_ylen = (multi_pos[3,i]-multi_pos[1,i])
     multi_center = [multi_pos[0,i] + multi_xlen/2d, multi_pos[1,i] + multi_ylen/2d]

     multi_size = [xsize*multi_xlen, ysize*multi_ylen]
     multi_aspect = multi_size[1]/float(multi_size[0])
     new_aspect = aspect_ratio/multi_aspect
     if new_aspect gt 1 then begin
        y_factor = 1.
        x_factor = 1/new_aspect
     endif else begin          
        y_factor = new_aspect
        x_factor = 1.
     endelse
 
     new_xlen = multi_xlen*x_factor
     new_ylen = multi_ylen*y_factor
     new_multi = [multi_center[0] - new_xlen/2d, multi_center[1] - new_ylen*y_factor/2d, $
                  multi_center[0] + new_xlen/2d, multi_center[1] + new_ylen*y_factor/2d]
  
     new_pos[*,i] = [new_xlen * plot_pos[0] + new_multi[0], new_ylen * plot_pos[1] + new_multi[1], $
                new_xlen * plot_pos[2] + new_multi[0], new_ylen * plot_pos[3] + new_multi[1]]
     
     new_cb_pos[*,i] = [new_xlen * cb_pos[0] + new_multi[0], new_ylen * cb_pos[1] + new_multi[1], $
                   new_xlen * cb_pos[2] + new_multi[0], new_ylen * cb_pos[3] + new_multi[1]]
  endfor
  plot_pos = new_pos
  cb_pos = new_cb_pos

  window_num = 1
  if windowavailable(window_num) then begin 
     wset, window_num
     if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
  endif else make_win = 1
  if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize

  if keyword_set(pub) then begin
     charthick = 3
     thick = 3
     xthick = 3
     ythick = 3
     font = 1
     charsize = 2

     pson, file = plotfile, /eps 
  endif else begin
     charthick = 1
     thick = 1
     xthick = 1
     ythick = 1
     font = -1
     charsize = 1

  endelse

  if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 13, /brewer
  axkeywords = {xstyle: 1, ystyle: 1, thick: thick, charthick: charthick,  $
                charsize: charsize, font: font, xthick: xthick, ythick: ythick} 

  cgimage, beam_norm, position = plot_pos[*,1], xrange = range, yrange = range, /nointerp, $
           axkeywords = axkeywords, /axes, xtitle = textoidl('u (\lambda)', font = font), $
           ytitle = textoidl('v (\lambda)', font = font)

  cgcontour, beam_use, uv_arr_use, uv_arr_use, /onimage, levels=levels, c_linestyle=linestyles, $
             c_colors=colors, c_annotation=annotations, charsize = charsize, font = font
             
  cgcolorbar, position = cb_pos[*,1], /vertical, annotatecolor = 'black', minrange=0, maxrange = 1, $
              title = 'Baseline Power Response', charsize = charsize, charthick = charthick, xthick = xthick, $
              ythick = ythick, font = font


  ;;cgimage, beam_deg_norm, /noerase, position = plot_pos[*,0], xrange = deg_range, yrange = deg_range, $
  ;;         axkeywords = axkeywords, /axes, xtitle = textoidl('\theta_x (degrees)', font = font), $
  ;;         ytitle = textoidl('\theta_x (degrees)', font = font), title = number_formatter(frequencies[0]) + ' MHz'

  ;;cgcontour, beam_deg_use, deg_arr_use, deg_arr_use, /onimage, levels=levels, $
  ;;           c_colors='dark_grey', c_annotation=annotations, charsize = charsize, font = font

  ;;cgcolorbar, position = cb_pos[*,0], /vertical, annotatecolor = 'black', minrange=0, maxrange = 1, $
  ;;            title = 'Baseline Power Response', charsize = charsize, charthick = charthick, xthick = xthick, $
  ;;            ythick = ythick, font = font

  nlevels = 7
  deg_fractions = 10^(findgen(nlevels)*3/(nlevels-1)-3)

  deg_levels = max(beam_deg_use)*deg_fractions

  nlevels = n_elements(deg_levels)
  deg_annotations = strarr(n_elements(nlevels)) + ' '
  ticknames = strarr(nlevels)+ ' '
  wh_decade = where(round(((alog10(deg_fractions) mod 1)*10))/10. eq 0, count_decade)
  if count_decade eq 0 then stop
  ticknames[wh_decade] = number_formatter(deg_fractions[wh_decade], format = '(e0)',/print_exp)

  cgloadct, 0, /reverse, ncolors=nlevels, bottom=1

  cgcontour, beam_deg_use, deg_arr_use, deg_arr_use, levels=deg_levels, /noerase, /fill, position = plot_pos[*,0], $
             xrange = deg_range, yrange = deg_range, xtitle = textoidl('\theta_x (degrees)', font = font), $
             ytitle = textoidl('\theta_x (degrees)', font = font), $
             c_colors = indgen(nlevels)+1, c_annotation=deg_annotations, xstyle=1, ystyle=1, thick=thick, charthick=charthick,  $
             charsize=charsize, font=font, xthick = xthick, ythick = ythick

  cgcolorbar, position = cb_pos[*,0], /vertical, annotatecolor = 'black', range=minmax(deg_fractions), divisions= nlevels-1, $
              ncolors = nlevels-1, bottom=2, /discrete, /ylog, ticknames = ticknames, title = 'Baseline Power Response', $
              charsize = charsize, charthick = charthick, xthick = xthick, ythick = ythick, font = font

  if keyword_set(pub) then begin
     psoff
     wdelete, window_num
  endif


  ;; Main figure

  if keyword_set(grey_scale) then plotfile = base_path('plots') + 'single_use/multibaseline_fig_grey.eps' $
  else plotfile = base_path('plots') + 'single_use/multibaseline_fig.eps'

  u_pix = 945
  v_pix = 945

  pix_vals = reform(sim_cube[u_pix, v_pix, *])
  undefine, sim_cube
  weights = reform(weights_cube[u_pix, v_pix, *])
  undefine, weights_cube
  sigma = 1/sqrt(weights)
  if min(weights) eq 0 then stop
  
  input_val = model_uv[u_pix, v_pix]
  amp_residual = abs(pix_vals) - abs(input_val)
  phase_residual = (atan(pix_vals,/phase) - atan(input_val,/phase)) * 180d/ !dpi
  
  error_scale_factor = 0.05
  sigma = sigma * error_scale_factor
  amp_error = sqrt(2d) * sigma
  phase_error = (sigma / abs(pix_vals)) * 180d/ !dpi

  if max(abs(phase_residual)) ge 180 then begin
     ;; unwind phases
     wh1 = where(phase_residual lt -180, count1)
     if count1 gt 0 then phase_residual[wh1] = phase_residual[wh1] + 360d

     wh2 = where(phase_residual gt 180, count2)
     if count2 gt 0 then phase_residual[wh2] = phase_residual[wh2] - 360d
  endif

  amp_range = [-1,1] * max(abs(amp_residual))
  res_phase_range = [-1,1] * max(abs(phase_residual))

  uv_loc = [uv_arr[u_pix], uv_arr[v_pix]]

  baseline_dist = sqrt((baseline_u - uv_loc[0])^2d + (baseline_v - uv_loc[1])^2d)
  wh_contrib = where(baseline_dist / rebin(reform(beam_radii, 1, n_freq), nbaselines, n_freq) le 1d, count_wh)
  if count_wh eq 0 then stop

  baseline_inds = wh_contrib mod nbaselines
  freq_inds = wh_contrib / nbaselines
  ncontrib_freq = histogram(freq_inds, min=0, max=n_freq-1)
  ncontrib_baseline = histogram(baseline_inds, locations = baseline_nums, reverse_indices = ri)
  wh_hist_n0 = where(ncontrib_baseline gt 0)
  baseline_contrib = baseline_nums[wh_hist_n0]
  nb_contrib = n_elements(baseline_contrib)

  bcontrib_locs = [[[baseline_u[baseline_contrib, *]]], [[baseline_v[baseline_contrib, *]]]]
  mask = intarr(nb_contrib, n_freq)
  for i=0L, nb_contrib-1 do mask[i, freq_inds[ri[ri[wh_hist_n0[i]]:ri[wh_hist_n0[i]+1]-1]]] = 1

  ;;xrange = minmax(bcontrib_locs[*,*,0]) + [-1d, 1d]*beam_radii
  ;;yrange = minmax(bcontrib_locs[*,*,1]) + [-1d, 1d]*beam_radii

  xrange = uv_loc[0] + [-1d, 1d]*max(beam_radii)*4
  yrange = uv_loc[1] + [-1d, 1d]*max(beam_radii)*4

  wh_u = where(uv_arr gt xrange[0] and uv_arr lt xrange[1], count_u)
  if count_u eq 0 then stop
  u_inds_range = [max([min(wh_u)-1, 0]), min([max(wh_u)+1, n_u])]
  
  wh_v = where(uv_arr gt yrange[0] and uv_arr lt yrange[1], count_v)
  if count_v eq 0 then stop
  v_inds_range = [max([min(wh_v)-1, 0]), min([max(wh_v)+1, n_u])]

  xrange = uv_arr[u_inds_range]
  yrange = uv_arr[v_inds_range]
 
  map_size = [(u_inds_range[1]-u_inds_range[0]), (v_inds_range[1]-v_inds_range[0])] * psf_resolution

  uv_plot_arr = atan(model_uv[wh_u,*], /phase)
  ;; uv_plot_arr = congrid(uv_plot_arr[*, wh_v], count_u*10, count_v*10)
  uv_plot_arr = congrid(uv_plot_arr[*, wh_v], map_size[0], map_size[1], cubic=-0.5)
  phase_range = [min(uv_plot_arr), max(uv_plot_arr)]

  ;; normalize to +/- pi
  norm_plot_arr = (uv_plot_arr - phase_range[0])*255/(phase_range[1] - phase_range[0])
  
  ;; make array for heat map & for outline
  weights_map = dblarr(map_size[0], map_size[1], n_freq)
  influence_map = dblarr(map_size[0], map_size[1], n_freq)
  map_uv_pix = [u_pix-u_inds_range[0], v_pix-v_inds_range[0]]*psf_resolution
  map_uv_arr =[[[indgen(map_size[0]) * uv_binsize / psf_resolution + xrange[0]]], $
               [[indgen(map_size[1]) * uv_binsize / psf_resolution] + yrange[0]]]

  window, 3
  cgloadct, 13, /brewer
  val_add = dblarr(nb_contrib, n_freq) 
  for j=0, n_freq-1 do begin
     for i=0, nb_contrib-1 do begin
        if mask[i,j] eq 0 then continue
        pix_center = round((bcontrib_locs[i,j,*] - uv_loc) / beam_uv_binsize) + map_uv_pix
        pix_u_range = pix_center[0]-n_uv_beam/2 + [0, n_uv_beam-1]
        pix_v_range = pix_center[1]-n_uv_beam/2 + [0, n_uv_beam-1]
        beam_pix_uvloc = pix_center - map_uv_pix + n_uv_beam/2

        beam_u_range = [0, n_uv_beam-1]
        beam_v_range = [0, n_uv_beam-1]
        if pix_u_range[0] lt 0 then begin
           beam_u_range[0] = (-1)*pix_u_range[0]
           pix_u_range[0] = 0
        endif
        if pix_v_range[0] lt 0 then begin
           beam_v_range[0] = (-1)*pix_v_range[0]
           pix_v_range[0] = 0
        endif
        if pix_u_range[1] gt map_size[0]-1 then begin
           beam_u_range[1] = (n_uv_beam-1) + (map_size[0]-1) - pix_u_range[1]
           pix_u_range[1] = map_size[0] - 1
        endif
        if pix_v_range[1] gt map_size[1]-1 then begin
           beam_v_range[1] = (n_uv_beam-1) + (map_size[1]-1) - pix_v_range[1]
           pix_v_range[1] = map_size[1]-1
        endif

        beam_to_add = beam_uv[beam_u_range[0]:beam_u_range[1],beam_v_range[0]:beam_v_range[1],j]
        weight_at_uvpix = beam_uv[beam_pix_uvloc[0], beam_pix_uvloc[1], j]

        weights_map[pix_u_range[0]:pix_u_range[1], pix_v_range[0]:pix_v_range[1],j] += beam_to_add
        influence_map[pix_u_range[0]:pix_u_range[1], pix_v_range[0]:pix_v_range[1],j] += weight_at_uvpix * beam_to_add

        cube_inds = reform(round((bcontrib_locs[i,j,*] - uv_loc) / uv_binsize)+[u_pix,v_pix])
        model_val = model_uv[cube_inds[0],cube_inds[1]]
        val_add[i,j] = model_val*weight_at_uvpix

     endfor
;;stop
     ;; make contours to ouline all contributing baselines
     level = max(beam_uv[*,*,j])*0.002d
     cgcontour, weights_map[*,*,j], map_uv_arr[*,0], map_uv_arr[*,0], levels=level, $
                path_info = outline_contour_info, path_xy = outline_contour_xy, /path_data_coords

     plot, map_uv_arr[*,0], map_uv_arr[*,1], /xstyle,/ystyle
     cgimage, weights_map[*,*,j], /noerase, /overplot,/scale
     oplot, outline_contour_xy[0,*], outline_contour_xy[1,*], color=255
  
     tag = 'freq_' + number_formatter(j)
     if j eq 0 then outline_shape = create_struct(tag, outline_contour_xy) else $
        outline_shape = create_struct(outline_shape, tag, outline_contour_xy)
  endfor
  wdelete, 3

  ;; weights_color_range = [0, 255]
  ;; weights_n_colors = weights_color_range[1] - weights_color_range[0]
  ;; norm_weight_map = weights_map * weights_n_colors / max(weights_map) + weights_color_range[0]
  influence_color_range = [0, 255]
  influence_n_colors = influence_color_range[1] - influence_color_range[0]

  norm_influence_map = influence_map * 0d
  for j=0, n_freq -1 do norm_influence_map[*,*,j] = (influence_map[*,*,j]/max(influence_map[*,*,j])) * influence_n_colors + $
     influence_color_range[0]

  ;;freq_inds_plot = [8, 18, 24]
  freq_inds_plot = [8, 18, 24]
  n_freq_plot = n_elements(freq_inds_plot)
 

  tvlct, r, g, b, /get

  positions = fltarr(4, 3*n_freq_plot + 1)
  cb_positions = fltarr(4, 3)

  cb_size_punits = 0.1
  margin1_punits = [0.2, 0.2]
  margin2_punits = [0.05, 0.1]

  xlen_punits = 4.*(margin1_punits[0] + margin2_punits[0]) + 3. + cb_size_punits
  ylen_punits = 4.*(margin1_punits[1] + margin2_punits[1] + 1.)

  plot_size = 1./[xlen_punits, ylen_punits]
  cb_size = [cb_size_punits, 1.] / [xlen_punits, ylen_punits]
  margin1 = margin1_punits / [xlen_punits, ylen_punits]
  margin2 = margin2_punits / [xlen_punits, ylen_punits]

  for i=0, n_freq_plot-1 do begin
     positions[0, 3*i:3*(i+1)-1] = (i+1.)*margin1[0] + i*(plot_size[0] + margin2[0])
     positions[2, 3*i:3*(i+1)-1] = (i+1.)*(margin1[0] + plot_size[0]) + i*margin2[0]
  endfor
  
  for i=0, 2 do begin
     positions[1, indgen(n_freq_plot)*3+i] = (4-i)*margin1[1] + (3-i)*(plot_size[1] + margin2[1])
     positions[3, indgen(n_freq_plot)*3+i] = (4-i)*(margin1[1] + plot_size[1]) + (3-i)*(margin2[1])
     cb_positions[1, i] = (4-i)*margin1[1] + (3-i)*(plot_size[1] + margin2[1])
     cb_positions[3, i] = (4-i)*(margin1[1] + plot_size[1]) + (3-i)*(margin2[1])
  endfor

  positions[*, n_freq_plot*3] = [margin1[0], margin1[1], 1-margin2[0], margin1[1] + plot_size[1]]
  cb_positions[0, *] = 4.*margin1[0] + 3*(plot_size[0] + margin2[0])
  cb_positions[2, *] = 1.-margin2[0]

  size_factor = 200
  xsize = xlen_punits * size_factor
  ysize = ylen_punits * size_factor
  window_num = 2

  if windowavailable(window_num) then begin 
     wset, window_num
     if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
  endif else make_win = 1
  if make_win then cgdisplay, wid=window_num, xsize = xsize, ysize = ysize, color='white'
  cgerase, 'white'

   if keyword_set(pub) then begin
     charthick = 2
     thick = 2
     xthick = 2
     ythick = 2
     font = 1
     charsize = 1

     pson, file = plotfile, /eps 
  endif else begin
     charthick = 1
     thick = 1
     xthick = 1
     ythick = 1
     font = -1
     charsize = 0.8

  endelse

  if keyword_set(grey_scale) then begin
     residual_plot_color = 'dark grey'
     uv_mark_color = 'black'
     arrow_color = 'black'
     freq_line_color = 'black'
     baseline_loc_color = 'charcoal'
  endif else begin
     residual_plot_color = 'red6'
     uv_mark_color = 'red6'
     ;;arrow_color = 'grn5'
     ;;arrow_color = 'red6'
     arrow_color = 'dark grey'
     freq_line_color = 'pbg5'
     baseline_loc_color = 'black'
  endelse

  phase_plot_range = [-1,1]*max(abs(phase_residual))
  delta_freq = frequencies[1]-frequencies[0]
  freq_plot = [frequencies[0]-delta_freq/2, frequencies, frequencies[n_freq-1]+delta_freq/2]
  res_plot = [phase_residual[0], phase_residual, phase_residual[n_freq-1]]
  cgplot, freq_plot, res_plot*0, yrange = phase_plot_range, xstyle=1, xtitle = 'f (MHz)', ytick_get = yticks, ytitle = 'degrees', $
          title = 'Residual Phase', position = positions[*, n_freq_plot*3], psym=-3, axiscolor='black', color = 'black', $
          charsize = charsize, thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, font = font
  for i=0, n_freq_plot-1 do cgplot, /overplot, rebin([frequencies[freq_inds_plot[i]]]-0.2*delta_freq, 2), minmax(yticks), $
                                    color = freq_line_color, linestyle=2, thick=thick

  temp = residual_plot_color
  oploterror, frequencies, phase_residual, phase_error, psym=10, errcolor = temp, $
              /nohat, thick = thick, errthick = thick
  cgplot, /overplot, freq_plot, res_plot, psym=10, color = residual_plot_color, thick = thick

  contrib_inds_set = where(total(mask[*,freq_inds_plot], 2) gt 0, count_contrib_set)

  xrange = minmax(map_uv_arr[*,0])
  yrange = minmax(map_uv_arr[*,1])
  for k=0, n_freq_plot-1 do begin
     j = freq_inds_plot[k]

     contrib_inds = where(mask[*,j] eq 1, count_contrib, complement = noncontrib_inds, ncomplement = count_noncontrib)
     contrib_dists = sqrt((reform(bcontrib_locs[contrib_inds[*],j,0]) - uv_loc[0])^2d + $
                          (reform(bcontrib_locs[contrib_inds[*],j,1]) - uv_loc[1])^2d)
     dist_order = sort(contrib_dists)

     cgplot, [uv_loc[0], uv_loc[1]], color = 'black', xrange = xrange, yrange = yrange, xstyle=5, ystyle=5, /nodata, $
           title = number_formatter(frequencies[j]) + ' MHz', position = positions[*, 3*k], /noerase, charsize = charsize, $
           charthick = charthick, font = font 

     if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 13, /brewer
     cgimage, norm_influence_map[*,*,j], /noerase, /overplot, /nointerp


     for i=0, nb_contrib-1 do cgplot, /overplot, bcontrib_locs[i,j,0]+beam_shape.(j)[0,*], $
                                      bcontrib_locs[i,j,1]+beam_shape.(j)[1,*], psym=-3, color = 'light grey', $
                                      thick = thick/4d
     ;; for i=0, count_noncontrib-1 do cgplot, /overplot, bcontrib_locs[noncontrib_inds[i],j,0]+beam_shape.(j)[0,*], $
     ;;                                     bcontrib_locs[noncontrib_inds[i],j,1]+beam_shape.(j)[1,*], psym=-3, color = 'light grey', $
     ;;                                       thick = thick
     ;; for i=0, count_contrib_set-1 do cgplot, /overplot, bcontrib_locs[contrib_inds_set[i],j,0]+beam_shape.(j)[0,*], $
     ;;                                        bcontrib_locs[contrib_inds_set[i],j,1]+beam_shape.(j)[1,*], psym=-3, $
     ;;                                        color = 'light grey', thick = thick
     ;; for i=0, count_contrib-1 do cgplot, /overplot, bcontrib_locs[contrib_inds[i],j,0]+beam_shape.(j)[0,*], $
     ;;                                     bcontrib_locs[contrib_inds[i],j,1]+beam_shape.(j)[1,*], psym=-3, color = 'grey', $
     ;;                                     thick = thick

     for i=0, n_freq_plot-1 do begin
        j2 = freq_inds_plot[i]
        if j2 eq j then continue
        freq_ratio = frequencies[j]/frequencies[j2]

        cgplot, /overplot,outline_shape.(j2)[0,*]*freq_ratio, outline_shape.(j2)[1,*]*freq_ratio, psym=-3, color = 'grey', $
                thick = thick
     endfor

     cgplot, /overplot,outline_shape.(j)[0,*], outline_shape.(j)[1,*], psym=-3, color = 'black', thick = thick

     ;; arrow_slope = 1
     ;; arrow_end_x = [[0], [-161],[-158]]
     ;; arrow_end_y = [[0], [-169],[-168]]
     ;; arrow_start_x = [[0], [-154],[-152]]
     ;; arrow_start_y = (-1)*arrow_slope*(arrow_end_x - arrow_start_x) + arrow_end_y
     ;; cgarrow, arrow_start_x[k], arrow_start_y[k], arrow_end_x[k], arrow_end_y[k], color =  arrow_color, /solid, /data, $
     ;;          hsize = !D.X_SIZE / 96., thick=thick 
     

     cgplot, /overplot, [uv_loc[0]], [uv_loc[1]], psym=1, thick=2, color = uv_mark_color

     cgaxis, xaxis=0, xtick_get = xticks, xtitle = textoidl('u (\lambda)', font = font), xrange = xrange, xstyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, xthick = xthick, font = font
     cgaxis, yaxis=0, ytick_get = yticks, ytitle = textoidl('v (\lambda)', font = font), yrange = yrange, ystyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, ythick = ythick, font = font
     cgaxis, xaxis=1, xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), xrange = xrange, xstyle = 1, $
             color = 'black', xthick = xthick
     cgaxis, yaxis=1, ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), yrange = yrange, ystyle = 1, $
             color = 'black', ythick = ythick

     

     cgplot, [uv_loc[0]], [uv_loc[1]], color = 'black', xrange = xrange, yrange = yrange, xstyle=5, ystyle=5, /nodata, $
           title = number_formatter(frequencies[j]) + ' MHz', position = positions[*, 3*k+1], /noerase, charsize = charsize, $
           charthick = charthick, font = font

     if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 25, /brewer
     cgimage, norm_plot_arr, /overplot,/noerase

     for i=0, count_contrib-1 do cgplot, /overplot, [bcontrib_locs[contrib_inds[i],j,0]], [bcontrib_locs[contrib_inds[i],j,1]], $
                                         psym = 16, color = baseline_loc_color, symsize=0.6
     cgplot, /overplot,outline_shape.(j)[0,*], outline_shape.(j)[1,*], psym=-3, color = 'black', thick = thick

     cgplot, /overplot, [uv_loc[0]], [uv_loc[1]], psym=1, thick=2, color = uv_mark_color

     cgaxis, xaxis=0, xtick_get = xticks, xtitle = textoidl('u (\lambda)', font = font), xrange = xrange, xstyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, xthick = xthick, font = font
     cgaxis, yaxis=0, ytick_get = yticks, ytitle = textoidl('v (\lambda)', font = font), yrange = yrange, ystyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, ythick = ythick, font = font
     cgaxis, xaxis=1, xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), xrange = xrange, xstyle = 1, $
             color = 'black', xthick = xthick
     cgaxis, yaxis=1, ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), yrange = yrange, ystyle = 1, $
             color = 'black', ythick = ythick


     cgplot, [uv_loc[0]], [uv_loc[1]], color = 'black', xrange = xrange, yrange = yrange, xstyle=5, ystyle=5, /nodata, $
             title = number_formatter(frequencies[j]) + ' MHz', position = positions[*, 3*k+2], /noerase, charsize = charsize, $
             charthick = charthick, font = font

     if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 25, /brewer
     tvlct, r1, g1, b1, /get
     dims = size(norm_plot_arr, /dimension)
     norm_plot_arr_rgb = dblarr(3, dims[0], dims[1])
     norm_plot_arr_rgb[0,*,*] = r1[norm_plot_arr]
     norm_plot_arr_rgb[1,*,*] = g1[norm_plot_arr]
     norm_plot_arr_rgb[2,*,*] = b1[norm_plot_arr]
     
     cgloadct, 0
     tvlct, r2, g2, b2, /get
     dims = size(norm_influence_map[*,*,j], /dimension)
     norm_influence_map_rgb = bytarr(3, dims[0], dims[1])
     norm_influence_map_rgb[0,*,*] = r2[norm_influence_map[*,*,j]]
     norm_influence_map_rgb[1,*,*] = g2[norm_influence_map[*,*,j]]
     norm_influence_map_rgb[2,*,*] = b2[norm_influence_map[*,*,j]]
     
     ;; alpha = dblarr(3, dims[0], dims[1])
     ;; for i=0,2 do alpha[i,*,*] = (norm_influence_map[*,*,j] / max(norm_influence_map[*,*,j]))^(1/6d)
     ;; blend_image = norm_plot_arr_rgb * alpha + (255d) * (1-alpha)

     ;;blend_image = norm_plot_arr_rgb*0.75 + norm_influence_map_rgb * .25
     ;multiplier = (1d - norm_influence_map_rgb/255d)*.5+.5
     multiplier = (norm_influence_map_rgb/255d)*.3+.6
     blend_image = norm_plot_arr_rgb * multiplier

     cgimage, blend_image, /overplot,/noerase
     ;;cgblendimage, norm_plot_arr_rgb, norm_influence_map_rgb, alpha=.75, /overplot,/noerase
  
     cgplot, /overplot,outline_shape.(j)[0,*], outline_shape.(j)[1,*], psym=-3, color = 'black', thick = thick

     cgplot, /overplot, [uv_loc[0]], [uv_loc[1]], psym=1, thick=2, color = uv_mark_color

     cgaxis, xaxis=0, xtick_get = xticks, xtitle = textoidl('u (\lambda)', font = font), xrange = xrange, xstyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, xthick = xthick, font = font
     cgaxis, yaxis=0, ytick_get = yticks, ytitle = textoidl('v (\lambda)', font = font), yrange = yrange, ystyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, ythick = ythick, font = font
     cgaxis, xaxis=1, xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), xrange = xrange, xstyle = 1, $
             color = 'black', xthick = xthick
     cgaxis, yaxis=1, ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), yrange = yrange, ystyle = 1, $
             color = 'black', ythick = ythick


     ;; cgplot, [uv_loc[0]], [uv_loc[1]], color = 'black', xrange = xrange, yrange = yrange, xstyle=1, ystyle=1, /nodata, $
     ;;       title = number_formatter(frequencies[j]) + ' MHz', xtitle = textoidl('u (\lambda)', font = font), $
     ;;       ytitle = textoidl('v (\lambda)', font = font), position = positions[*, 3*k+2], /noerase, charsize = charsize, $
     ;;       thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, font = font

     ;; for i=0, count_contrib-1 do begin
     ;;    ind = (reverse(dist_order))[i]

     ;;    bcontrib_pix = reform(round(bcontrib_locs[contrib_inds[ind],j,*] / uv_binsize + [n_u/2, n_v/2]))
     ;;    bcontrib_phase = atan(model_uv[bcontrib_pix[0], bcontrib_pix[1]], /phase)
     ;;    color = round((bcontrib_phase - phase_range[0])*255/(phase_range[1] - phase_range[0]))
 
     ;;    if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 25, /brewer

     ;;    cgcolorfill, bcontrib_locs[contrib_inds[ind],j,0]+reform(beam_shape.(j)[0,*]), $
     ;;                 bcontrib_locs[contrib_inds[ind],j,1]+reform(beam_shape.(j)[1,*]), color = color

     ;;    cgplot, /overplot, bcontrib_locs[contrib_inds[ind],j,0]+beam_shape.(j)[0,*], $
     ;;           bcontrib_locs[contrib_inds[ind],j,1]+beam_shape.(j)[1,*], psym=-3, color = 'black', thick=thick

     ;; endfor

     ;; cgplot, /overplot, [uv_loc[0]], [uv_loc[1]], psym=1, thick=2, color = uv_mark_color

  endfor

  if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 13, /brewer
  cgcolorbar, position = cb_positions[*, 0], /vertical, annotatecolor = 'black', minrange=0, maxrange = 1, $
              title = 'Influence', charsize = charsize, charthick = charthick, xthick = xthick, $
              ythick = ythick, font = font

  if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 25, /brewer
  cgcolorbar, position = cb_positions[*, 1], /vertical, annotatecolor = 'black', minor=5, minrange=phase_range[0]*180d/!dpi, $
              maxrange = phase_range[1]*180d/!dpi, title = 'Phase (degrees)', charsize = charsize, $
              charthick = charthick, xthick = xthick, ythick = ythick, font = font

  if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 25, /brewer
  cgcolorbar, position = cb_positions[*, 2], /vertical, annotatecolor = 'black', minor=5, minrange=phase_range[0]*180d/!dpi, $
              maxrange = phase_range[1]*180d/!dpi, title = 'Phase (degrees)', charsize = charsize, $
              charthick = charthick, xthick = xthick, ythick = ythick, font = font

  if keyword_set(pub) then begin
     psoff
     wdelete, window_num
  endif


  ;; Main figure movie frames
  positions = fltarr(4, 4)
  cb_positions = fltarr(4, 3)

  cb_size_punits = 0.1
  margin1_punits = [0.4, 0.35]
  margin2_punits = [0.1, 0.25]

  xlen_punits = 6.*(margin1_punits[0] + margin2_punits[0]) + 3. + 3.*cb_size_punits
  ylen_punits = 2.*(margin1_punits[1] + margin2_punits[1] + 1.)

  plot_size = 1./[xlen_punits, ylen_punits]
  cb_size = [cb_size_punits, 1.] / [xlen_punits, ylen_punits]
  margin1 = margin1_punits / [xlen_punits, ylen_punits]
  margin2 = margin2_punits / [xlen_punits, ylen_punits]


  for i=0, 2 do begin
     positions[0, i] = (i+1.)*margin1[0] + i*(plot_size[0] + cb_size[0] + 2.*margin2[0] + margin1[0])
     positions[2, i] = (i+1.)*(margin1[0] + plot_size[0]) + i*(cb_size[0] + 2.*margin2[0] + margin1[0])
     cb_positions[0, i] = (i+1.)*(2.*margin1[0] + plot_size[0] + margin2[0]) + i*(cb_size[0] + margin2[0])
     cb_positions[2, i] = (i+1.)*(2.*margin1[0] + plot_size[0] + margin2[0] + cb_size[0]) + i*margin2[0]
  endfor

  positions[1, 0:2] = 2.*margin1[1] + (plot_size[1] + margin2[1])
  positions[3, 0:2] = 2.*(margin1[1] + plot_size[1]) + margin2[1]
  cb_positions[1, *] = 2.*margin1[1] + (plot_size[1] + margin2[1])
  cb_positions[3, *] = 2.*(margin1[1] + plot_size[1]) + margin2[1]

 
  positions[*, 3] = [margin1[0], margin1[1], 1-margin2[0], margin1[1] + plot_size[1]]

  size_factor = 200
  xsize = xlen_punits * size_factor
  ysize = ylen_punits * size_factor
  window_num = 3

  if windowavailable(window_num) then begin 
     wset, window_num
     if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
  endif else make_win = 1
  if make_win then cgdisplay, wid=window_num, xsize = xsize, ysize = ysize, color='white'

   if keyword_set(pub) then begin
     charthick = 3
     thick = 3
     xthick = 3
     ythick = 3
     font = 1
     charsize = 2
  endif 

  for j=0, n_freq-1 do begin
     if keyword_set(grey_scale) then plotfile = base_path('plots') + 'single_use/multibaseline_movie_frames/multibaseline_freq' + $
        number_formatter(j) + '_grey.eps' $
     else plotfile = base_path('plots') + 'single_use/multibaseline_movie_frames/multibaseline_freq' + number_formatter(j) + '.eps'

     cgerase, 'white'
     if keyword_set(pub) then pson, file = plotfile, /eps 

     ;; bottom plot
     cgplot, freq_plot, res_plot*0, yrange = phase_plot_range, xstyle=1, xtitle = 'f (MHz)', ytick_get = yticks, ytitle = 'degrees', $
             title = 'Residual Phase', position = positions[*, 3], psym=-3, axiscolor='black', color = 'black', $
             charsize = charsize, thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, font = font
     cgplot, /overplot, rebin([frequencies[j]]-0.1*delta_freq, 2), minmax(yticks), color = freq_line_color, linestyle=2, thick=thick

     temp = residual_plot_color
     oploterror, frequencies, phase_residual, phase_error, psym=10, errcolor = temp, $
                 /nohat, thick = thick, errthick = thick
     cgplot, /overplot, freq_plot, res_plot, psym=10, color = residual_plot_color, thick = thick

     ;; influence plot
     contrib_inds = where(mask[*,j] eq 1, count_contrib, complement = noncontrib_inds, ncomplement = count_noncontrib)
     contrib_dists = sqrt((reform(bcontrib_locs[contrib_inds[*],j,0]) - uv_loc[0])^2d + $
                          (reform(bcontrib_locs[contrib_inds[*],j,1]) - uv_loc[1])^2d)
     dist_order = sort(contrib_dists)

     cgplot, [uv_loc[0], uv_loc[1]], color = 'black', xrange = xrange, yrange = yrange, xstyle=5, ystyle=5, /nodata, $
           title = number_formatter(frequencies[j]) + ' MHz', position = positions[*, 0], /noerase, charsize = charsize, $
           charthick = charthick, font = font 

     if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 13, /brewer
     cgimage, norm_influence_map[*,*,j], /noerase, /overplot, /nointerp


     for i=0, nb_contrib-1 do cgplot, /overplot, bcontrib_locs[i,j,0]+beam_shape.(j)[0,*], $
                                      bcontrib_locs[i,j,1]+beam_shape.(j)[1,*], psym=-3, color = 'light grey', $
                                      thick = thick/4d

     cgplot, /overplot,outline_shape.(j)[0,*], outline_shape.(j)[1,*], psym=-3, color = 'black', thick = thick
     

     cgplot, /overplot, [uv_loc[0]], [uv_loc[1]], psym=1, thick=2, color = uv_mark_color

     cgaxis, xaxis=0, xtick_get = xticks, xtitle = textoidl('u (\lambda)', font = font), xrange = xrange, xstyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, xthick = xthick, font = font
     cgaxis, yaxis=0, ytick_get = yticks, ytitle = textoidl('v (\lambda)', font = font), yrange = yrange, ystyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, ythick = ythick, font = font
     cgaxis, xaxis=1, xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), xrange = xrange, xstyle = 1, $
             color = 'black', xthick = xthick
     cgaxis, yaxis=1, ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), yrange = yrange, ystyle = 1, $
             color = 'black', ythick = ythick

     cgcolorbar, position = cb_positions[*, 0], /vertical, annotatecolor = 'black', minrange=0, maxrange = 1, $
                 title = 'Influence', charsize = charsize, charthick = charthick, xthick = xthick, $
                 ythick = ythick, font = font


     ;; phase plot with integration region & baselines
     cgplot, [uv_loc[0]], [uv_loc[1]], color = 'black', xrange = xrange, yrange = yrange, xstyle=5, ystyle=5, /nodata, $
           title = number_formatter(frequencies[j]) + ' MHz', position = positions[*, 1], /noerase, charsize = charsize, $
           charthick = charthick, font = font

     if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 25, /brewer
     cgimage, norm_plot_arr, /overplot,/noerase

     for i=0, count_contrib-1 do cgplot, /overplot, [bcontrib_locs[contrib_inds[i],j,0]], [bcontrib_locs[contrib_inds[i],j,1]], $
                                         psym = 16, color = baseline_loc_color, symsize=0.6
     cgplot, /overplot,outline_shape.(j)[0,*], outline_shape.(j)[1,*], psym=-3, color = 'black', thick = thick

     cgplot, /overplot, [uv_loc[0]], [uv_loc[1]], psym=1, thick=2, color = uv_mark_color

     cgaxis, xaxis=0, xtick_get = xticks, xtitle = textoidl('u (\lambda)', font = font), xrange = xrange, xstyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, xthick = xthick, font = font
     cgaxis, yaxis=0, ytick_get = yticks, ytitle = textoidl('v (\lambda)', font = font), yrange = yrange, ystyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, ythick = ythick, font = font
     cgaxis, xaxis=1, xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), xrange = xrange, xstyle = 1, $
             color = 'black', xthick = xthick
     cgaxis, yaxis=1, ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), yrange = yrange, ystyle = 1, $
             color = 'black', ythick = ythick

     cgcolorbar, position = cb_positions[*, 1], /vertical, annotatecolor = 'black', minor=5, minrange=phase_range[0]*180d/!dpi, $
                 maxrange = phase_range[1]*180d/!dpi, title = 'Phase (degrees)', charsize = charsize, $
                 charthick = charthick, xthick = xthick, ythick = ythick, font = font

     ;; phase plot with influence transparency
    cgplot, [uv_loc[0]], [uv_loc[1]], color = 'black', xrange = xrange, yrange = yrange, xstyle=5, ystyle=5, /nodata, $
             title = number_formatter(frequencies[j]) + ' MHz', position = positions[*, 2], /noerase, charsize = charsize, $
             charthick = charthick, font = font

     if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 25, /brewer
     tvlct, r1, g1, b1, /get
     dims = size(norm_plot_arr, /dimension)
     norm_plot_arr_rgb = dblarr(3, dims[0], dims[1])
     norm_plot_arr_rgb[0,*,*] = r1[norm_plot_arr]
     norm_plot_arr_rgb[1,*,*] = g1[norm_plot_arr]
     norm_plot_arr_rgb[2,*,*] = b1[norm_plot_arr]
     
     cgloadct, 0
     tvlct, r2, g2, b2, /get
     dims = size(norm_influence_map[*,*,j], /dimension)
     norm_influence_map_rgb = bytarr(3, dims[0], dims[1])
     norm_influence_map_rgb[0,*,*] = r2[norm_influence_map[*,*,j]]
     norm_influence_map_rgb[1,*,*] = g2[norm_influence_map[*,*,j]]
     norm_influence_map_rgb[2,*,*] = b2[norm_influence_map[*,*,j]]
     
     multiplier = (norm_influence_map_rgb/255d)*.3+.6
     blend_image = norm_plot_arr_rgb * multiplier

     cgimage, blend_image, /overplot,/noerase
  
     cgplot, /overplot,outline_shape.(j)[0,*], outline_shape.(j)[1,*], psym=-3, color = 'black', thick = thick

     cgplot, /overplot, [uv_loc[0]], [uv_loc[1]], psym=1, thick=2, color = uv_mark_color

     cgaxis, xaxis=0, xtick_get = xticks, xtitle = textoidl('u (\lambda)', font = font), xrange = xrange, xstyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, xthick = xthick, font = font
     cgaxis, yaxis=0, ytick_get = yticks, ytitle = textoidl('v (\lambda)', font = font), yrange = yrange, ystyle = 1, $
             color = 'black', charsize = charsize, charthick = charthick, ythick = ythick, font = font
     cgaxis, xaxis=1, xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), xrange = xrange, xstyle = 1, $
             color = 'black', xthick = xthick
     cgaxis, yaxis=1, ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), yrange = yrange, ystyle = 1, $
             color = 'black', ythick = ythick

     if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 25, /brewer
     cgcolorbar, position = cb_positions[*, 2], /vertical, annotatecolor = 'black', minor=5, minrange=phase_range[0]*180d/!dpi, $
                 maxrange = phase_range[1]*180d/!dpi, title = 'Phase (degrees)', charsize = charsize, $
                 charthick = charthick, xthick = xthick, ythick = ythick, font = font

     if keyword_set(pub) then psoff

  endfor
  tvlct, r, g, b

  if keyword_set(pub) then wdelete, window_num
     


  froot = base_path('data') + 'fhd_simulations_old/'
  simfile = froot + 'sim_496t_fullsky_uvf_holo_2dkpower.idlsave'

  if keyword_set(grey_scale) then plotfile = base_path('plots') + 'single_use/multibaseline_ps_grey.eps' $
  else plotfile = base_path('plots') + 'single_use/multibaseline_ps.eps'

  info_file = froot + 'sim_496t_info.idlsave'
  restore, info_file   
  xy_length = 2048
  deg_offset = xy_length * degpix / sqrt(2d)
  rad_offset = deg_offset * !pi/180d

  redshift = 8
  cosmology_measures, redshift, wedge_factor = wedge_factor
  source_dists = rad_offset
  wedge_amp = wedge_factor * source_dists

  kperp_plot_range = [6e-3, 0.2]
  data_range = [1e13, 1e18]

  kpower_2d_plots, simfile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = data_range, $
                   /no_title, pub = pub, grey_scale = grey_scale, /baseline_axis, /plot_wedge_line, wedge_amp = wedge_amp, $
                   plotfile = plotfile, window_num = 4


end
