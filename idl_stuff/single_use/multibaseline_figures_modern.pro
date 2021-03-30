pro multibaseline_figures_modern, use_sim_residual = use_sim_residual, $
   make_movie_frames = make_movie_frames, plot_2d_ps = plot_2d_ps, $
   refresh = refresh, png = png, pdf=pdf, eps=eps, grey_scale = grey_scale

  if n_elements(use_sim_residual) eq 0 then use_sim_residual = 1
  if n_elements(make_movie_frames) eq 0 then make_movie_frames = 0
  if n_elements(plot_2d_ps) eq 0 then plot_2d_ps = 0

  folder_name = base_path('data') + 'fhd_kernel_window_control_run_leave_flags_write_uvf_single_source'
  obs_name = '1094485088_first_half_flagged'
  savefile = folder_name + path_sep()+ 'multibaseline_data_modern.idlsave'
  ftest = file_test(savefile) *  (1 - file_test(savefile, /zero_length))

   model_obs_name = '1094485088_whole_obs_no_flags'

   plot_exten = ''
   if keyword_set(png) then begin
      plot_exten = '.png'
      delete_ps = 1
   endif else if keyword_set(pdf) then begin
      plot_exten = '.pdf'
      delete_ps = 1
   endif else if keyword_set(eps) then begin
      plot_exten = '.eps'
      delete_ps = 0
   endif

   if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then begin
      pub = 1
      movie_frame_path = base_path('plots') + 'single_use/multibaseline_modern_movie_frames
      if not file_test(movie_frame_path, /directory) then begin
         file_mkdir, movie_frame_path
      endif
   endif else begin
      pub = 0
   endelse

   if ftest eq 0 or keyword_set(refresh) then begin

      refresh_info = 0
      uvf_input = 1
      exact_obsnames = 1

      ; sky_model_file = filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')
      ; sky_model_file = folder_name  + path_sep() + 'test_source_1Jy_RA_350_dec_-26.8.sav'
      sky_model_file = filepath('test_source_1Jy_RA_350_dec_-26.8.sav',root=rootdir('FHD'),subdir='catalog_data/simulation')

      obs_info = ps_filenames(folder_name, obs_name, $
         exact_obsnames = exact_obsnames, uvf_input = uvf_input, $
         data_subdirs = data_subdirs, ps_foldernames = ps_foldername, $
         save_paths = save_path, plot_paths = plot_path, refresh_info = refresh_info)

      if obs_info.info_files[0] ne '' then begin
         datafile = obs_info.info_files[0]
      endif else begin
         datafile = obs_info.cube_files.(0)
      endelse

      if not file_test(save_path, /directory) then file_mkdir, save_path

      uvf_options = create_uvf_options(delta_uv_lambda = delta_uv_lambda, $
         max_uv_lambda = max_uv_lambda, full_image = full_image, image_clip = image_clip, $
         uv_avg = uv_avg, uv_img_clip = uv_img_clip, require_radec = require_radec, $
         dft_fchunk = dft_fchunk, no_dft_progress = no_dft_progress)

      ps_options = create_ps_options(ave_removal = ave_removal, wt_cutoffs = wt_cutoffs, $
         wt_measures = wt_measures, spec_window_type = spec_window_type, $
         no_spec_window = no_spec_window, allow_beam_approx = allow_beam_approx, $
         std_power = std_power, no_wtd_avg = no_wtd_avg, $
         inverse_covar_weight = inverse_covar_weight)

      refresh_options = create_refresh_options(refresh_dft = refresh_dft, $
      refresh_beam = refresh_beam, refresh_ps = refresh_ps, refresh_kcube = refresh_ps, $
      refresh_binning = refresh_binning, refresh_info = refresh_info)

      file_struct_arr = fhd_file_setup(datafile, beamfile = beamfiles, /sim, $
         uvf_input = uvf_input, savefilebase = savefilebase, save_path = save_path, $
         freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
         refresh_info = refresh_options.refresh_info, uvf_options = uvf_options, $
         ps_options = ps_options)


      model_obs_info = ps_filenames(folder_name, model_obs_name, $
         exact_obsnames = exact_obsnames, uvf_input = uvf_input, $
         data_subdirs = data_subdirs, ps_foldernames = ps_foldername, $
         save_paths = save_path, plot_paths = plot_path, refresh_info = refresh_info)

      if model_obs_info.info_files[0] ne '' then begin
         model_datafile = model_obs_info.info_files[0]
      endif else begin
         model_datafile = model_obs_info.cube_files.(0)
      endelse

      model_file_struct_arr = fhd_file_setup(model_datafile, beamfile = beamfiles, /sim, $
         uvf_input = uvf_input, savefilebase = savefilebase, save_path = save_path, $
         freq_ch_range = freq_ch_range, freq_flags = freq_flags, freq_flag_name = freq_flag_name, $
         refresh_info = refresh_options.refresh_info, uvf_options = uvf_options, $
         ps_options = ps_options)

      ; just use xx pol for now:
      file_struct = file_struct_arr[0]
      model_file_struct = model_file_struct_arr[0]

      filebase = cgRootName(file_struct.datafile[0], directory=froot)
      end_pos_even = strpos(strmid(filebase, 0), '_even')
      end_pos_odd = strpos(strmid(filebase, 0), '_odd')
      end_pos_cube = strpos(strmid(filebase, 0), '_gridded') ;; always > -1
      end_pos = end_pos_even > end_pos_odd
      if end_pos eq -1 then end_pos = end_pos_cube
      obs_name = strmid(filebase, 0, end_pos)

      metadata_path = folder_name + path_sep() + 'metadata'
      obs_file = metadata_path + path_sep() + obs_name + '_obs.sav'
      layout_file = metadata_path + path_sep() + obs_name + '_layout.sav'
      params_file = metadata_path + path_sep() + obs_name + '_params.sav'
      ; the psf structure is in the beams file:
      beams_file = folder_name + path_sep() + 'beams' + path_sep() + obs_name + '_beams.sav'
      jones_file = folder_name + path_sep() + 'beams' + path_sep() + obs_name + '_jones.sav'

      obs = getvar_savefile(obs_file, 'obs')
      n_tiles = obs.n_tile

      uv_binsize = obs.kpix

      layout = getvar_savefile(layout_file, 'layout')
      array_center_ecef = layout.array_center

      ; see wikipedia geodetic_datum and Datum transformations of
      ; GPS positions PDF in pyuvdata's docs/references folder
      gps_b = 6356752.31424518
      gps_a = 6378137
      e_squared = 6.69437999014e-3
      e_prime_squared = 6.73949674228e-3
      gps_p = sqrt(array_center_ecef[0]^2 + array_center_ecef[1]^2)
      gps_theta = atan(array_center_ecef[2] * gps_a, gps_p * gps_b)
      latitude = atan(array_center_ecef[2] + e_prime_squared * gps_b * sin(gps_theta)^3, $
                     gps_p - e_squared * gps_a * cos(gps_theta)^3)

      longitude = atan(array_center_ecef[1], array_center_ecef[0])
      gps_n = gps_a / sqrt(1 - e_squared * sin(latitude)^2)
      altitude = (gps_p / cos(latitude)) - gps_n

      antenna_enu = dblarr(3, n_tiles)
      antenna_enu[*, 0] = (-sin(longitude) * layout.antenna_coords[*, 0] + $
                           cos(longitude) * layout.antenna_coords[*, 1])
      antenna_enu[*, 1] = (-sin(latitude) * cos(longitude) * layout.antenna_coords[*, 0] $
                           - sin(latitude) * sin(longitude) * layout.antenna_coords[*, 1] $
                           + cos(latitude) * layout.antenna_coords[*, 2])
      antenna_enu[*, 2] = (cos(latitude) * cos(longitude) * layout.antenna_coords[*, 0] $
                           + cos(latitude) * sin(longitude) * layout.antenna_coords[*, 1] $
                           + sin(latitude) * layout.antenna_coords[*, 2])

      tile_A = (*obs.baseline_info).tile_a
      tile_B = (*obs.baseline_info).tile_b

      ; if obs.nbaseline_times ne ((n_tiles^2-n_tiles)/2 + n_tiles) then stop
      nbaseline_times = n_elements(tile_A)
      params = getvar_savefile(params_file, 'params')

      frequencies = file_struct.frequencies
      n_freq = n_elements(frequencies)

      baseline_u = rebin(params.uu, nbaseline_times, n_freq, /sample) * $
                  rebin(reform(frequencies * 1e6, 1, n_freq), nbaseline_times, n_freq, /sample)
      baseline_v = rebin(params.vv, nbaseline_times, n_freq, /sample) * $
                  rebin(reform(frequencies * 1e6, 1, n_freq), nbaseline_times, n_freq, /sample)

      ;; get beams
      psf = getvar_savefile(beams_file, 'psf')
      bl_ind_use = 0
      n_beam_freq = max(psf.fbin_i) + 1
      beam_uv = dblarr(psf.dim, psf.dim, n_beam_freq)
      for fi=0, n_beam_freq-1 do begin
         beam_uv[*, *, fi] = $
            reform(*(*(*psf.beam_ptr)[file_struct.pol_index, fi, bl_ind_use])[0,0], $
                  psf.dim, psf.dim)
      endfor
      beam_fbin = psf.fbin_i

      ;; normalize to make integral = 1
      beam_uv = beam_uv / rebin(reform($
      total(total(beam_uv, 2), 1), 1, 1, n_beam_freq), psf.dim, psf.dim, n_beam_freq, /sample)
      beam_uv_arr  = (dindgen(psf.dim) - psf.dim/2d) * uv_binsize

      beam_cutoff_factor = 0.001d
      beam_radii = dblarr(n_beam_freq)
      window, 3
      for i = 0, n_beam_freq-1 do begin
         level = max(beam_uv[*,*,i])*0.001d

         cgcontour, beam_uv[*,*,i], beam_uv_arr, beam_uv_arr, levels=level, $
                     path_info = beam_contour_info, path_xy = beam_contour_xy, /path_data_coords
         beam_radii[i] = max(sqrt(beam_contour_xy[0,*]^2d + beam_contour_xy[1,*]^2d))

         quick_image, beam_uv[*,*,i], beam_uv_arr, beam_uv_arr, /log, window=3, $
            title='Frequency bin ' + number_formatter(i), xtitle = 'u (pixel number)', $
            ytitle = 'v (pixel number)', cb_title = 'beam value'
         cgplot, beam_contour_xy[0,*], beam_contour_xy[1,*], color=255, /over

         tag = 'freq_bin_' + number_formatter(i)
         if i eq 0 then beam_shape = create_struct(tag, beam_contour_xy) else $
            beam_shape = create_struct(beam_shape, tag, beam_contour_xy)
         
      endfor
      wdelete, 3

      nfiles = n_elements(file_struct.datafile)

      weights_cube1 = get_cube_uvf_input(file_struct.weightfile[0], $
         file_struct.weightvar, n_freq, file_struct.pol_index)
      if nfiles eq 2 then begin
         weights_cube2 = get_cube_uvf_input(file_struct.weightfile[1], $
            file_struct.weightvar, n_freq, file_struct.pol_index)
         weights_cube = temporary(weights_cube1) + temporary(weights_cube2)
      endif else weights_cube = temporary(weights_cube1)

      data_cube1 = get_cube_uvf_input(file_struct.datafile[0], $
         file_struct.datavar, n_freq, file_struct.pol_index)
      if nfiles eq 2 then begin
         data_cube2 = get_cube_uvf_input(file_struct.datafile[1], $
            file_struct.datavar, n_freq, file_struct.pol_index)
         uvf_cube = temporary(data_cube1) + temporary(data_cube2)
      endif else uvf_cube = temporary(data_cube1)

      model_wt_cube1 = get_cube_uvf_input(model_file_struct.datafile[0], $
         model_file_struct.weightvar, n_freq, file_struct.pol_index)
      if nfiles eq 2 then begin
         model_wt_cube2 = get_cube_uvf_input(model_file_struct.datafile[1], $
            model_file_struct.weightvar, n_freq, file_struct.pol_index)
         model_weights_cube = temporary(model_wt_cube1) + temporary(model_wt_cube2)
      endif else model_weights_cube = temporary(model_wt_cube1)

      model_cube1 = get_cube_uvf_input(model_file_struct.datafile[0], $
         model_file_struct.datavar, n_freq, file_struct.pol_index)
      if nfiles eq 2 then begin
         model_cube2 = get_cube_uvf_input(model_file_struct.datafile[1], $
            model_file_struct.datavar, n_freq, file_struct.pol_index)
         model_uvf_cube = temporary(model_cube1) + temporary(model_cube2)
      endif else model_uvf_cube = temporary(model_cube1)


;       quick_image, (abs(uvf_cube))[*,200,*], title='|uvf_cube|'
;       quick_image, (abs(model_uvf_cube))[*,200,*], window=2, title='|model_cube|'
;       uvf_abs_diff = abs(model_uvf_cube) - abs(uvf_cube)
;       quick_image, uvf_abs_diff[*,200,*], window=3, title='|model_cube|-|uvf_cube|'
;       uvf_abs_diff_ratio = uvf_abs_diff/abs(model_uvf_cube)
;       uvf_abs_diff_ratio[where(abs(model_uvf_cube) eq 0)] = 0
;       quick_image, uvf_abs_diff_ratio[*,200,*], window=7, title='(|model_cube|-|uvf_cube|)/|model_cube|', /log

;       quick_image, weights_cube[*,200,*], window=4, title='weights_cube'
;       quick_image, model_weights_cube[*,200,*], window=5, title='model_weights_cube'
;       weights_diff = model_weights_cube - weights_cube
;       quick_image, weights_diff[*,200,*], window=6, title='model_weights_cube-weights_cube'
;       weights_diff_ratio = weights_diff/model_weights_cube
;       weights_diff_ratio[where(abs(model_weights_cube) eq 0)] = 0
;       quick_image, weights_diff_ratio[*,200,*], window=8, title='(model_weights_cube-weights_cube)/model_weights_cube', /log
; stop

      sim_cube = temporary(uvf_cube) / weights_cube
      wh_wt0 = where(weights_cube eq 0, count_wt0)
      if count_wt0 gt 0 then sim_cube[wh_wt0] = 0

      sim_model_cube = temporary(model_uvf_cube) / model_weights_cube
      wh_wt0 = where(model_weights_cube eq 0, count_wt0)
      if count_wt0 gt 0 then sim_model_cube[wh_wt0] = 0

      dims = size(sim_cube, /dimension)
      n_u = dims[0]
      n_v = dims[1]
      uv_arr = (dindgen(n_u) - n_u/2d) * uv_binsize

      catalog = getvar_savefile(sky_model_file, 'catalog')
      apply_astrometry, obs, ra_arr=catalog.ra, dec_arr=catalog.dec, x_arr=x_arr, y_arr=y_arr, /ad2xy
      catalog.x=x_arr
      catalog.y=y_arr

      jones = getvar_savefile(jones_file, 'jones')

      source_dft_multi,obs,jones,catalog,model_uv_ptr,spectral_uv_full,xvals=xvals,yvals=yvals,uv_i_use=uv_i_use,$
         conserve_memory=conserve_memory,frequency=frequency,dft_threshold=dft_threshold,silent=silent,$
         dimension=dimension,elements=elements,n_pol=n_pol,spectral_model_uv_arr=spectral_model_uv_arr,$
         n_spectral=n_spectral,flatten_spectrum=flatten_spectrum,double_precision=double_precision,$
         gaussian_source_models = gaussian_source_models,_Extra=extra

      model_uv = *model_uv_ptr[file_struct.pol_index]

      psf_dim = psf.dim
      ntimes = obs.n_time
      degpix = obs.degpix

      save, file = savefile, sim_cube, weights_cube, sim_model_cube, model_weights_cube, $
            model_uv, uv_arr, $
            baseline_u, baseline_v, frequencies, beam_uv, beam_shape, $
            beam_radii, nbaseline_times, ntimes, degpix, beam_uv_arr, n_freq, $
            psf_dim, beam_fbin, n_u, n_v, uv_binsize, beam_cutoff_factor
   endif else restore, savefile


   ;; Beam figure
   if keyword_set(grey_scale) then begin
      plotfile = base_path('plots') + 'single_use/multibaseline_modern_beam_grey' + plot_exten
   endif else begin
      plotfile = base_path('plots') + 'single_use/multibaseline_modern_beam' + plot_exten
   endelse

   ; beam_freq_ind = 16
   beam_freq_ind = 0

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
   beam_theta = abs(fft_shift(fft(beam_uv[*,*,beam_freq_ind])) * psf_dim^2. * uv_binsize^2.)
   beam_degpix = (1/(psf_dim * uv_binsize))*180/!pi
   beam_deg_arr = (dindgen(psf_dim) - psf_dim/2d) * beam_degpix

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

      cgps_open, plotfile, /font, encapsulated=eps;, /nomatch, inches=sizes.inches, $
         ;xsize=sizes.xsize, ysize=sizes.ysize, xoffset=sizes.xoffset, $
         ;yoffset=sizes.yoffset, landscape = landscape
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
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
      wdelete, window_num
   endif


   ;; Main figure

   if keyword_set(grey_scale) then begin
      plotfile = base_path('plots') + 'single_use/multibaseline_modern_fig_grey' + plot_exten
   endif else begin
      plotfile = base_path('plots') + 'single_use/multibaseline_modern_fig' + plot_exten
   endelse

   bl_keep_fraction = 1 ; 0.1
   if bl_keep_fraction ne 1 then begin
      nblts_to_keep = round(nbaseline_times * bl_keep_fraction)
      baseline_u = baseline_u[indgen(nblts_to_keep), *]
      baseline_v = baseline_v[indgen(nblts_to_keep), *]
   endif

   shade_freqs_range = [181, 188]
   sim_uv_extent = (size(sim_cube, /dim))[0]
   ; model_uv has a much larger extent, need to trim it down
   model_uv_extent = (size(model_uv, /dim))[0]

   model_uv_range_use = [model_uv_extent/2 - sim_uv_extent/2 + 1, model_uv_extent/2 + sim_uv_extent/2]
   model_uv_use = model_uv[model_uv_range_use[0]:model_uv_range_use[1], model_uv_range_use[0]:model_uv_range_use[1]]

   u_pix = 200
   v_pix = 110
   uv_loc = [uv_arr[u_pix], uv_arr[v_pix]]

   ; quick_image, atan(sim_cube[u_pix,*,*],/phase), uv_arr, frequencies, window = 4
   ; cgplot, /over, [uv_loc[1], uv_loc[1]], [0, max(frequencies)], color="black"
   ; cgplot, /over, [min(uv_arr), max(uv_arr)], [shade_freqs_range[0], shade_freqs_range[0]], color="black"
   ; cgplot, /over, [min(uv_arr), max(uv_arr)], [shade_freqs_range[1], shade_freqs_range[1]], color="black"

   quick_image, weights_cube[u_pix,*,*], uv_arr, frequencies, window = 5, $
      ytitle='frequency (MHz)', xtitle='v', title='Weights'
   cgplot, /over, [uv_loc[1], uv_loc[1]], [0, max(frequencies)], color="black"
   cgplot, /over, [min(uv_arr), max(uv_arr)], [shade_freqs_range[0], shade_freqs_range[0]], color="black"
   cgplot, /over, [min(uv_arr), max(uv_arr)], [shade_freqs_range[1], shade_freqs_range[1]], color="black"

   model_uv_cube = rebin(reform(real_part(model_uv_use), sim_uv_extent, sim_uv_extent, 1), $
      sim_uv_extent, sim_uv_extent, n_freq, /sample) $
      + complex(0,1) * rebin(reform(imaginary(model_uv_use), sim_uv_extent, sim_uv_extent, 1), $
         sim_uv_extent, sim_uv_extent, n_freq, /sample)

   sim_model_diff = sim_cube - model_uv_cube
   quick_image, atan(sim_model_diff[u_pix,*,*],/phase), uv_arr, frequencies, window = 4, $
      ytitle='frequency (MHz)', xtitle='v', title='Phase (sim - model)'
   cgplot, /over, [uv_loc[1], uv_loc[1]], [0, max(frequencies)], color="black"
   cgplot, /over, [min(uv_arr), max(uv_arr)], [shade_freqs_range[0], shade_freqs_range[0]], color="black"
   cgplot, /over, [min(uv_arr), max(uv_arr)], [shade_freqs_range[1], shade_freqs_range[1]], color="black"

   sim_model_diff2 = sim_cube - sim_model_cube
   quick_image, atan(sim_model_diff2[u_pix,*,*],/phase), uv_arr, frequencies, window = 6, $
      ytitle='frequency (MHz)', xtitle='v', title='Phase (sim - sim_model)'
   cgplot, /over, [uv_loc[1], uv_loc[1]], [0, max(frequencies)], color="black"
   cgplot, /over, [min(uv_arr), max(uv_arr)], [shade_freqs_range[0], shade_freqs_range[0]], color="black"
   cgplot, /over, [min(uv_arr), max(uv_arr)], [shade_freqs_range[1], shade_freqs_range[1]], color="black"

; stop
   undefine_fhd, sim_model_diff, sim_model_diff2

   ; window, 3
   ; u_range = round(5 * 2 / uv_binsize) * [-1, 1] + uv_loc[0]
   ; v_range = round(5 * 2 / uv_binsize) * [-1, 1] + uv_loc[1]
   ; cgplot, reform(baseline_u, n_elements(baseline_u)), reform(baseline_v, n_elements(baseline_u)), $
   ;     xrange = u_range, yrange = v_range
   ; cgplot, /overplot, [uv_loc[0]], [uv_loc[1]], psym=1, thick=2, color = "red"
   ; stop

   pix_vals = reform(sim_cube[u_pix, v_pix, *])
   ; undefine, sim_cube
   weights = reform(weights_cube[u_pix, v_pix, *])
   ; undefine, weights_cube
   sigma = 1/sqrt(weights)
   if min(weights) eq 0 then stop

   ; model_uv has a much larger extent, need to trim it down
   model_uv_extent = (size(model_uv, /dim))[0]

   model_uv_range_use = [model_uv_extent/2 - sim_uv_extent/2 + 1, model_uv_extent/2 + sim_uv_extent/2]
   model_uv_use = model_uv[model_uv_range_use[0]:model_uv_range_use[1], model_uv_range_use[0]:model_uv_range_use[1]]
   input_val = model_uv_use[u_pix, v_pix]
   amp_residual = abs(pix_vals) - abs(input_val)
   phase_residual = (atan(pix_vals,/phase) - atan(input_val,/phase)) * 180d/ !dpi

   ; same calc but with sim_model
   sim_input_val = reform(sim_model_cube[u_pix, v_pix, *])
   sim_amp_residual = abs(pix_vals) - abs(sim_input_val)
   sim_phase_residual = (atan(pix_vals,/phase) - atan(sim_input_val,/phase)) * 180d/ !dpi

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

   sim_amp_range = [-1,1] * max(abs(sim_amp_residual))
   sim_res_phase_range = [-1,1] * max(abs(sim_phase_residual))

   baseline_dist = sqrt((baseline_u - uv_loc[0])^2d + (baseline_v - uv_loc[1])^2d)
   beam_radii_freq = dblarr(n_freq)
   for fbin_i=0, max(beam_fbin) do begin
      fbin_inds = where(beam_fbin eq fbin_i, count_fbin)
      if count_fbin gt 0 then beam_radii_freq[fbin_inds] = beam_radii[fbin_i]
   endfor
   if bl_keep_fraction ne 1 then begin
      wh_contrib = where(baseline_dist / rebin(reform(beam_radii_freq, 1, n_freq), nblts_to_keep, n_freq, /sample) le 1d, count_wh)
   endif else begin
      wh_contrib = where(baseline_dist / rebin(reform(beam_radii_freq, 1, n_freq), nbaseline_times, n_freq, /sample) le 1d, count_wh)
   endelse
   if count_wh eq 0 then begin
      print, "No contributing baselines at any frequencies"
      stop
   endif

   if bl_keep_fraction ne 1 then begin
      baseline_inds = wh_contrib mod nblts_to_keep
      freq_inds = wh_contrib / nblts_to_keep
   endif else begin
      baseline_inds = wh_contrib mod nbaseline_times
      freq_inds = wh_contrib / nbaseline_times
   endelse
   ncontrib_freq = histogram(freq_inds, min=0, max=n_freq-1)
   if min(ncontrib_freq) eq 0 then begin
      print, "No contributing baselines at some frequencies"
      stop
   endif
   ncontrib_baseline = histogram(baseline_inds, locations = baseline_nums, reverse_indices = ri)
   wh_hist_n0 = where(ncontrib_baseline gt 0)
   baseline_contrib = baseline_nums[wh_hist_n0]
   nb_contrib = n_elements(baseline_contrib)

   bcontrib_locs = [[[baseline_u[baseline_contrib, *]]], [[baseline_v[baseline_contrib, *]]]]
   mask = intarr(nb_contrib, n_freq)
   for i=0L, nb_contrib-1 do mask[i, freq_inds[ri[ri[wh_hist_n0[i]]:ri[wh_hist_n0[i]+1]-1]]] = 1

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

   map_size = [(u_inds_range[1]-u_inds_range[0]), (v_inds_range[1]-v_inds_range[0])]

   uv_plot_arr = atan(model_uv_use[wh_u,*], /phase)
   ;; uv_plot_arr = congrid(uv_plot_arr[*, wh_v], count_u*10, count_v*10)
   uv_plot_arr = congrid(uv_plot_arr[*, wh_v], map_size[0], map_size[1], cubic=-0.5)
   phase_range = [min(uv_plot_arr), max(uv_plot_arr)]

   ;; normalize to +/- pi
   norm_plot_arr = (uv_plot_arr - phase_range[0])*255/(phase_range[1] - phase_range[0])

   ;; make array for heat map & for outline
   weights_map = dblarr(map_size[0], map_size[1], n_freq)
   influence_map = dblarr(map_size[0], map_size[1], n_freq)
   map_uv_pix = [u_pix-u_inds_range[0], v_pix-v_inds_range[0]]
   map_uv_arr =[[[indgen(map_size[0]) * uv_binsize + xrange[0]]], $
                  [[indgen(map_size[1]) * uv_binsize] + yrange[0]]]

   window, 3
   cgloadct, 13, /brewer
   val_add = dblarr(nb_contrib, n_freq)
   influence_peak_uv_loc = dblarr(2, n_freq)
   influence_cog_uv_loc = dblarr(2, n_freq)
   for j=0, n_freq-1 do begin
      for i=0, nb_contrib-1 do begin
         if mask[i,j] eq 0 then continue
         beam_find = beam_fbin[j]

         pix_center = round((bcontrib_locs[i,j,*] - uv_loc)) + map_uv_pix
         pix_u_range = pix_center[0]-psf_dim/2 + [0, psf_dim-1]
         pix_v_range = pix_center[1]-psf_dim/2 + [0, psf_dim-1]
         beam_pix_uvloc = pix_center - map_uv_pix + psf_dim/2

         beam_u_range = [0, psf_dim-1]
         beam_v_range = [0, psf_dim-1]
         if pix_u_range[0] lt 0 then begin
            beam_u_range[0] = (-1)*pix_u_range[0]
            pix_u_range[0] = 0
         endif
         if pix_v_range[0] lt 0 then begin
            beam_v_range[0] = (-1)*pix_v_range[0]
            pix_v_range[0] = 0
         endif
         if pix_u_range[1] gt map_size[0]-1 then begin
            beam_u_range[1] = (psf_dim-1) + (map_size[0]-1) - pix_u_range[1]
            pix_u_range[1] = map_size[0] - 1
         endif
         if pix_v_range[1] gt map_size[1]-1 then begin
            beam_v_range[1] = (psf_dim-1) + (map_size[1]-1) - pix_v_range[1]
            pix_v_range[1] = map_size[1]-1
         endif

         beam_to_add = beam_uv[beam_u_range[0]:beam_u_range[1],beam_v_range[0]:beam_v_range[1],beam_find]
         weight_at_uvpix = beam_uv[beam_pix_uvloc[0], beam_pix_uvloc[1], beam_find]

         weights_map[pix_u_range[0]:pix_u_range[1], pix_v_range[0]:pix_v_range[1],j] += beam_to_add
         influence_map[pix_u_range[0]:pix_u_range[1], pix_v_range[0]:pix_v_range[1],j] += weight_at_uvpix * beam_to_add

         cube_inds = reform(round((bcontrib_locs[i,j,*] - uv_loc) / uv_binsize)+[u_pix,v_pix])
         model_val = model_uv_use[cube_inds[0],cube_inds[1]]
         val_add[i,j] = model_val*weight_at_uvpix

      endfor
      if max(weights_map[*,*,j]) eq 0 then begin
         print, "no contributing baselines at freq ind " + string(j) + " (" + $
            number_formatter(frequencies[j]) + " MHz)"
         continue
      endif

      ;; make contours to ouline all contributing baselines
      level = max(beam_uv[*,*,beam_find])*beam_cutoff_factor
      cgcontour, weights_map[*,*,j], map_uv_arr[*,0], map_uv_arr[*,1], levels=level, $
                  path_info = outline_contour_info, path_xy = outline_contour_xy, /path_data_coords

      wh_max_influence = where(influence_map[*, *, j] eq max(influence_map[*, *, j]), count_infl_max)
      if count_infl_max eq 1 then begin
         influence_peak_uv_loc[0, j] = map_uv_arr[wh_max_influence mod map_size[0], 0]
         influence_peak_uv_loc[1, j] = map_uv_arr[wh_max_influence / map_size[0], 1]
      endif else begin
         stop
      endelse
     
      influence_cog_uv_loc[0, j] = total(influence_map[*, *, j] * $
         rebin(map_uv_arr[*,0], map_size[0], map_size[1]))/total(influence_map[*, *, j])
      influence_cog_uv_loc[1, j] = total(influence_map[*, *, j] * $
         rebin(reform(map_uv_arr[*,1], 1, map_size[1]), map_size[0], map_size[1]))/total(influence_map[*, *, j])

      plot, map_uv_arr[*,0], map_uv_arr[*,1], /xstyle,/ystyle
      cgimage, weights_map[*,*,j], /noerase, /overplot,/scale
      oplot, outline_contour_xy[0,*], outline_contour_xy[1,*], color=255

      tag = 'freq_' + number_formatter(j)
      if j eq 0 then outline_shape = create_struct(tag, outline_contour_xy) else $
         outline_shape = create_struct(outline_shape, tag, outline_contour_xy)
   endfor
   wdelete, 3


   ;; make center of mass track plot
   
   window_num=3
   track_xlen = max(influence_cog_uv_loc[0,*]) - min(influence_cog_uv_loc[0,*])
   track_ylen = max(influence_cog_uv_loc[1,*]) - min(influence_cog_uv_loc[1,*])
   track_maxlen = max([track_xlen, track_ylen])
   track_xcenter = min(influence_cog_uv_loc[0,*]) + track_xlen/2
   track_ycenter = min(influence_cog_uv_loc[1,*]) + track_ylen/2
   plot_xrange = track_xcenter + [-1,1]*track_maxlen/2 +[-0.5, 0.5]
   plot_yrange = track_ycenter + [-1,1]*track_maxlen/2 +[-0.5, 0.5]

   ;; Work out plot & colorbar positions
   ;; in units of plot area (incl. margins)
   cb_size = 0.025
   margin = [0.2, 0.2, 0.02, 0.1]
  
   cb_margin = [0.2, 0.02]
  
   plot_pos = [margin[0], margin[1], (1-cb_margin[1]-cb_size-cb_margin[0]-margin[2]), (1-margin[3])]
   cb_pos = [(1-cb_margin[1]-cb_size), margin[1], (1-cb_margin[1]), (1-margin[3])]
  
   plot_len = [plot_pos[2]-plot_pos[0], plot_pos[3] - plot_pos[1]]
   plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])
  
   data_aspect = 1
   aspect_ratio =  data_aspect /plot_aspect
   if aspect_ratio gt 1 then begin
      y_factor = aspect_ratio
      x_factor = 1.
   endif else begin
      y_factor = 1.
      x_factor = 1./aspect_ratio
   endelse
  
   screen_size = get_screen_size()
   max_xsize = screen_size[0]
   max_ysize = screen_size[1]
   base_size = 600

   xsize = round(base_size * x_factor)
   ysize = round(base_size * y_factor)

   if windowavailable(window_num) then begin
      wset, window_num
      if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
   endif else make_win = 1
   if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
         
   cgerase

   if keyword_set(pub) then begin
      ps_aspect = y_factor / x_factor
      
      if ps_aspect lt 1 then landscape = 1 else landscape = 0
      IF Keyword_Set(eps) THEN landscape = 0
      sizes = cgpswindow(LANDSCAPE=landscape, aspectRatio = ps_aspect, /sane_offsets)
            
      cgps_open, savefile, /font, encapsulated=eps, /nomatch, inches=sizes.inches, $
         xsize=sizes.xsize, ysize=sizes.ysize, $
         xoffset=sizes.xoffset, yoffset=sizes.yoffset, landscape = landscape

      font = 1

   endif else begin
      font = -1
   endelse
   charsize = 1.5

   colors = round(cgScaleVector(findgen(n_freq), 0, 255))

   infl_track_plotfile = base_path('plots') + 'single_use/multibaseline_modern_track' + plot_exten

   if keyword_set(pub) then begin
      cgps_open, infl_track_plotfile, /font, encapsulated=eps
   endif

   cgplot, influence_cog_uv_loc[0,*], influence_cog_uv_loc[1,*], /nodata, position = plot_pos, $
      charsize = charsize, font = font, $
      xrange=plot_xrange, yrange=plot_yrange, title='Influence Center of Mass', $
      xtitle=textoidl('u (\lambda)', font = font), ytitle=textoidl('v (\lambda)', font = font)
   cgloadct, 25, /brewer
   for j=0, n_freq-2 do begin
      if frequencies[j] ge shade_freqs_range[0] and frequencies[j] le shade_freqs_range[1] then begin
         thickness=4
      endif else begin
         thickness=1
      endelse
      cgplot, /over, [influence_cog_uv_loc[0, j], influence_cog_uv_loc[0, j+1]], $
         [influence_cog_uv_loc[1, j], influence_cog_uv_loc[1, j+1]], color=strtrim(colors[j],2), thick=thickness
   endfor
   for j=0, n_freq-1 do begin
      if frequencies[j] ge shade_freqs_range[0] and frequencies[j] le shade_freqs_range[1] then begin
         psym=16
         symsize=1
      endif else begin
         psym=9
         symsize=0.8
      endelse
      cgplot, /over, influence_cog_uv_loc[0, j], influence_cog_uv_loc[1, j], $
         psym=psym, symsize=symsize, color=strtrim(colors[j],2)
   endfor

   cgcolorbar, range=minmax(frequencies), position = cb_pos, /vertical, $
      charsize = charsize, font = font, title = "Frequency (MHz)"

   if keyword_set(pub) then begin
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
      wdelete, window_num
   endif

   ;; weights_color_range = [0, 255]
   ;; weights_n_colors = weights_color_range[1] - weights_color_range[0]
   ;; norm_weight_map = weights_map * weights_n_colors / max(weights_map) + weights_color_range[0]
   influence_color_range = [0, 255]
   influence_n_colors = influence_color_range[1] - influence_color_range[0]

   norm_influence_map = influence_map * 0d
   for j=0, n_freq -1 do norm_influence_map[*,*,j] = (influence_map[*,*,j]/max(influence_map[*,*,j])) * influence_n_colors + $
      influence_color_range[0]

   ;;freq_inds_plot = [8, 18, 24]
   freq_inds_plot = [80, 110, 150]
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

      cgps_open, plotfile, /font, encapsulated=eps;, /nomatch, inches=sizes.inches, $
         ;xsize=sizes.xsize, ysize=sizes.ysize, xoffset=sizes.xoffset, $
         ;yoffset=sizes.yoffset, landscape = landscape
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

   if use_sim_residual then begin
      phase_plot_range = sim_res_phase_range
      res_plot = [sim_phase_residual[0], sim_phase_residual, sim_phase_residual[n_freq-1]]
   endif else begin
      phase_plot_range = res_phase_range
      res_plot = [phase_residual[0], phase_residual, phase_residual[n_freq-1]]
   endelse

   delta_freq = frequencies[1]-frequencies[0]
   freq_plot = [frequencies[0]-delta_freq/2, frequencies, frequencies[n_freq-1]+delta_freq/2]
   cgplot, freq_plot, res_plot*0, yrange = phase_plot_range, xstyle=1, xtitle = 'f (MHz)', ytick_get = yticks, ytitle = 'degrees', $
            title = 'Residual Phase', position = positions[*, n_freq_plot*3], psym=-3, axiscolor='black', color = 'black', $
            charsize = charsize, thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, font = font
   cgcolorfill, [shade_freqs_range, reverse(shade_freqs_range)], $
      [phase_plot_range[0], phase_plot_range[0], phase_plot_range[1], phase_plot_range[1]]*.99, $
      color='light grey'
   cgplot, /overplot, [frequencies[0], frequencies[-1]], [0, 0], color='black'
   for i=0, n_freq_plot-1 do cgplot, /overplot, rebin([frequencies[freq_inds_plot[i]]]-0.2*delta_freq, 2, /sample), minmax(phase_plot_range), $
      color = freq_line_color, linestyle=2, thick=thick

   temp = residual_plot_color
   if use_sim_residual eq 0 then begin
      ; the error bars are based on the weights -- they only really make sense for the original
      ; comparison to the analytic uvplane, not to the comparison to unflagged gridded cube
      oploterror, frequencies, phase_residual, phase_error, psym=10, errcolor = temp, $
                  /nohat, thick = thick, errthick = thick
   endif
   cgplot, /overplot, freq_plot, res_plot, psym=10, color = residual_plot_color, thick = thick

   contrib_inds_set = where(total(mask[*,freq_inds_plot], 2) gt 0, count_contrib_set)

   xrange = minmax(map_uv_arr[*,0])
   yrange = minmax(map_uv_arr[*,1])
   for k=0, n_freq_plot-1 do begin
      j = freq_inds_plot[k]
      beam_find = beam_fbin[j]

      contrib_inds = where(mask[*,j] eq 1, count_contrib, complement = noncontrib_inds, ncomplement = count_noncontrib)
      contrib_dists = sqrt((reform(bcontrib_locs[contrib_inds[*],j,0]) - uv_loc[0])^2d + $
                           (reform(bcontrib_locs[contrib_inds[*],j,1]) - uv_loc[1])^2d)
      dist_order = sort(contrib_dists)

      cgplot, [uv_loc[0], uv_loc[1]], color = 'black', xrange = xrange, yrange = yrange, xstyle=5, ystyle=5, /nodata, $
            title = number_formatter(frequencies[j]) + ' MHz', position = positions[*, 3*k], /noerase, charsize = charsize, $
            charthick = charthick, font = font

      if keyword_set(grey_scale) then cgloadct, 0, /reverse else cgloadct, 13, /brewer
      cgimage, norm_influence_map[*,*,j], /noerase, /overplot, /nointerp


      for i=0, nb_contrib-1 do cgplot, /overplot, bcontrib_locs[i,j,0]+beam_shape.(beam_find)[0,*], $
                                       bcontrib_locs[i,j,1]+beam_shape.(beam_find)[1,*], psym=-3, color = 'light grey', $
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
      ;;    bcontrib_phase = atan(model_uv_use[bcontrib_pix[0], bcontrib_pix[1]], /phase)
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
      cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
      wdelete, window_num
   endif

   if make_movie_frames gt 0 then begin
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
         beam_find = beam_fbin[j]

         if keyword_set(grey_scale) then begin
            plotfile = base_path('plots') + 'single_use/multibaseline_modern_movie_frames/multibaseline_modern_freq' $
               + number_formatter(j) + '_grey' + plot_exten
         endif else begin
            plotfile = base_path('plots') + 'single_use/multibaseline_modern_movie_frames/multibaseline_modern_freq' $
               + number_formatter(j) + plot_exten
         endelse

         cgerase, 'white'
         if keyword_set(pub) then begin
            cgps_open, plotfile, /font, encapsulated=eps;, /nomatch, inches=sizes.inches, $
               ;xsize=sizes.xsize, ysize=sizes.ysize, xoffset=sizes.xoffset, $
               ;yoffset=sizes.yoffset, landscape = landscape
         endif

         ;; bottom plot
         cgplot, freq_plot, res_plot*0, yrange = phase_plot_range, xstyle=1, xtitle = 'f (MHz)', ytick_get = yticks, ytitle = 'degrees', $
                  title = 'Residual Phase', position = positions[*, 3], psym=-3, axiscolor='black', color = 'black', $
                  charsize = charsize, thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, font = font
         cgplot, /overplot, rebin([frequencies[j]]-0.1*delta_freq, 2, /sample), minmax(yticks), color = freq_line_color, linestyle=2, thick=thick

         temp = residual_plot_color
         if use_sim_residual eq 0 then begin
            ; the error bars are based on the weights -- they only really make sense for the original
            ; comparison to the analytic uvplane, not to the comparison to unflagged gridded cube
            oploterror, frequencies, phase_residual, phase_error, psym=10, errcolor = temp, $
                        /nohat, thick = thick, errthick = thick
         endif
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


         for i=0, nb_contrib-1 do cgplot, /overplot, bcontrib_locs[i,j,0]+beam_shape.(beam_find)[0,*], $
                                          bcontrib_locs[i,j,1]+beam_shape.(beam_find)[1,*], psym=-3, color = 'light grey', $
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

         if keyword_set(pub) then cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600


      endfor
      tvlct, r, g, b

      if keyword_set(pub) then wdelete, window_num
   endif

   if plot_2d_ps gt 0 then begin

      if keyword_set(grey_scale) then begin
         plotfile = base_path('plots') + 'single_use/multibaseline_modern_ps_grey' + plot_exten
      endif else begin
         plotfile = base_path('plots') + 'single_use/multibaseline_modern_ps' + plot_exten
      endelse

      info_file = folder_name + path_sep() + 'multibaseline_data_modern.idlsave'
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

      kpower_2d_plots, simfile, kperp_plot_range = kperp_plot_range, $
         kpar_plot_range = kpar_plot_range, data_range = data_range,  /no_title, $
         png = png, pdf=pdf, eps=eps, /baseline_axis, /plot_wedge_line, $
         wedge_amp = wedge_amp, plotfile = plotfile, window_num = 7

   endif

end
