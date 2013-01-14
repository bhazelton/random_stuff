
pro healpix_to_2dps, info_file, file_num = file_num, refresh = refresh, lsp = lsp, quiet = quiet, pub = pub, $
                     plot_unweighted = plot_unweighted, beam2 = beam2, plot_map = plot_map

  froot = base_path() + 'power_spectrum/healpix_maps/'
  if n_elements(info_file) eq 0 then info_file = froot + 'int_102.148MHz_info.txt'

  if n_elements(file_num) eq 0 then file_num = 0
  if n_elements(file_num) gt 1 then message, 'file_num must be a scalar.'
  
  info = read_info_file(info_file)
  nfiles = n_elements(info.files)
  if file_num ge n_elements(info.files) then message,  'file_num must be less than the number of files listed -- ' + nfiles
  

  ;; restore metadata
  test_file = file_test(info.metadata_file) *  (1 - file_test(info.metadata_file, /zero_length))  
  if test_file eq 0 then begin
     run_frame = 0
     answer=''
     read, answer, prompt = 'Metadata file does not exist. Run healpix_get_ft_frame now? (y/n)'
     case answer of
        'y':  run_frame = 1
        'yes': run_frame = 1
        else: run_frame = 0
     endcase
     if run_frame eq 1 then healpix_get_ft_frame, info_file else $
        message, 'healpix_get_ft_frame must be run to generate metadata before a power spectrum can be generated.'
  endif
  restore, info.metadata_file
  frequency = frequencies[file_num]

  map_file = info.files[file_num]
  read_fits_map, map_file, map_out, hdr, ehdr, nside=nside, ordering=ordering, coordsys=coordsys
  
  pix_nums = long(map_out[*,pix_ind])
  data = double(map_out[*,data_ind])
  if weight_ind ne -1 then weight = double(map_out[*,weight_ind])
  map_out = 0

  ;; find the indices for the included pixels
  min_val = min([pix_nums, pixels])
  wh = where(histogram(pix_nums, min=min_val, reverse_indices = ri) gt 0 AND histogram(pixels, min=min_val) gt 0, count)
  pix_inds = lonarr(count)
  for j = 0L, count-1 do pix_inds[j] = ri[ri[wh[j]]:ri[wh[j]+1]-1]
  wh = 0

  image_vals = data[pix_inds] * mk_conv_factors[file_num]
  if n_elements(weight) ne 0 then weight_vals = weight[pix_inds] * mk_conv_factors[file_num]

  if keyword_set(plot_map) then begin
     orthview, map_file, data_ind+1, /graticule, glsize = 1, rot = [0, -45], /half_sky, window=3
     if weight_ind ne -1 then orthview, map_file, weight_ind+1, /graticule, glsize = 1, rot = [0, -45], /half_sky, window=4

     plot_irreg_map, image_vals, x_rad_vals*180/!pi, y_rad_vals*180/!pi, ang_resolution*180/!pi, window_num = 5, $
                     title = 'Image (Jy/beam)', /interp

     if n_elements(weight_vals) ne 0 then $
        plot_irreg_map, weight_vals, x_rad_vals*180/!pi, y_rad_vals*180/!pi, ang_resolution*180/!pi, window_num = 6, $
                        title = 'Weights (Jy/beam)', /interp
  endif

  n_kx = n_elements(kx_rad_vals)
  n_ky = n_elements(ky_rad_vals)
  n_pix = n_elements(pix_inds)

  ;; do FT.
  if keyword_set(lsp) then begin
     savefile =  strsplit(map_file, '.fits', /regex, /extract) + '_LSP.idlsave'
     test_save = file_test(savefile) *  (1 - file_test(savefile, /zero_length))
     if test_save ne 0 and not keyword_set(refresh) then begin
        restore, savefile
        refresh = 0
        if n_elements(save_data) ne n_elements(data_vals_r) or n_elements(save_kx) ne n_elements(kx_rad_vals) or $
           n_elements(save_ky) ne n_elements(ky_rad_vals) or n_elements(transform) lt 1 then refresh = 1 $
        else if total(save_x-x_rad_vals) ne 0 or total(save_y-y_rad_vals) ne 0 or total(save_kx-kx_rad_vals) ne 0 or $
           total(save_ky-ky_rad_vals) ne 0 or total(save_data-image_vals) ne 0 then refresh = 1
     endif else refresh = 1
     
     if refresh eq 1 then begin
        ;; convention in lomb_scargle is for k/2pi
        ;;power = lomb_scargle(x_vals_r, y_vals_r, data_vals_r, kx_vals/(2d*!pi), ky_vals/(2d*!pi), timing = ls_time)
        data_ft = lomb_scargle_2D_fast(x_rad_vals, y_rad_vals, image_vals, kx_rad_vals/(2d*!pi), ky_rad_vals/(2d*!pi), $
                                       timing = ls_time, /transform) 
        if n_elements(weight_vals) ne 0 then $
           weight_ft = lomb_scargle_2D_fast(x_rad_vals, y_rad_vals, weight_vals, kx_rad_vals/(2d*!pi), ky_rad_vals/(2d*!pi), $
                                            timing = ls_time, /transform) $
        else weight_ft = data_ft*0 + 1d
        print, "Lomb-Scargle time: " + string(ls_time)

        save_x = x_rad_vals
        save_y = y_rad_vals
        save_data = image_vals
        if n_elements(weight_vals) ne 0 then save_weight = weight_vals
        save_kx = kx_rad_vals
        save_ky = ky_rad_vals
        save, file = savefile, data_ft, weight_ft, save_x, save_y, save_data, save_weight, save_kx, save_ky
     endif
     
  endif else begin
     savefile =  strsplit(map_file, '.fits', /regex, /extract) + '_DFT.idlsave'
     test_save = file_test(savefile) *  (1 - file_test(savefile, /zero_length))
     if test_save ne 0 and not keyword_set(refresh)then begin
        restore, savefile
        refresh = 0
        if n_elements(save_data) ne n_pix or n_elements(save_kx) ne n_kx or n_elements(save_ky) ne n_ky or $
           n_elements(data_ft) ne n_kx*n_ky or n_elements(weight_ft) ne n_kx*n_ky then refresh = 1 $
        else if total(save_x-x_rad_vals) ne 0 or total(save_y-y_rad_vals) ne 0 or total(save_kx-kx_rad_vals) ne 0 or $
           total(save_ky-ky_rad_vals) ne 0 or total(save_data-image_vals) ne 0 then refresh = 1
     endif else refresh = 1
     
     if refresh eq 1 then begin
        ;; transform_vec = discrete_ft_2D_vec(x_rad_vals, y_rad_vals, image_vals, kx_rad_vals, ky_rad_vals, $
        ;;                                    timing = ft_time)
        ;; print, "vectorized Discrete FT time: " + string(ft_time)
       
        if n_elements(weight_vals) ne 0 then begin
           arr = [[image_vals],[weight_vals]]
           transform = discrete_ft_2D_fast(x_rad_vals, y_rad_vals, arr, kx_rad_vals, ky_rad_vals, timing = ft_time)
           print, "partially vectorized Discrete FT time: " + string(ft_time)
           
           data_ft = reform(transform[*,*,0])
           weight_ft = reform(transform[*,*,1])
        endif else begin
           data_ft = discrete_ft_2D_fast(x_rad_vals, y_rad_vals, image_vals, kx_rad_vals, ky_rad_vals, $
                                         timing = ft_time)
           print, "partially vectorized Discrete FT time: " + string(ft_time)

           weight_ft = data_ft*0 + 1d
        endelse

        save_x = x_rad_vals
        save_y = y_rad_vals
        save_data = image_vals
        if n_elements(weight_vals) ne 0 then save_weight = weight_vals
        save_kx = kx_rad_vals
        save_ky = ky_rad_vals
        save, file = savefile, data_ft, weight_ft, save_x, save_y, save_data, save_weight, save_kx, save_ky
     endif
  endelse
  
  if keyword_set(beam2) then sigma = 1/abs(weight_ft)^2d $
  else sigma = 1/weight_ft
  wh = where(weight_ft eq 0, count)

  sigma2 = abs(sigma)^2
  if count gt 0 then sigma2[wh] = 0

  weighted_image = data_ft * sigma
  if count gt 0 then weighted_image[wh] = 0

  if keyword_set(plot_unweighted) then power = real_part(data_ft * conj(data_ft)) / (n_kx * n_ky) $
  else power = real_part(weighted_image * conj(weighted_image)) / (n_kx * n_ky)


  if not keyword_set(quiet) then begin
     ;; put kx & ky in 1/Mpc so I can plot them nicely
     z0_freq = 1420.40 ;; MHz
     ;; convert to MHz
     redshift = z0_freq/frequency - 1
     cosmology_measures, redshift, comoving_dist_los = comov_dist_los
     kx_mpc = kx_rad_vals / comov_dist_los
     ky_mpc = ky_rad_vals / comov_dist_los

    temp = strsplit(strmid(map_file, strpos(map_file, '/', /reverse_search)+1), '.fits', /regex, /extract)

     plotfile_path = base_path() + 'power_spectrum/plots/rts_healpix'
     plotfile = plotfile_path + temp + '_power'
     wt_plotfile = plotfile_path + temp + '_weight2'
     if keyword_set(plot_unweighted) then begin
        plotfile = plotfile + '_nowt'
        title = 'Unweighted P!Dk!N (mK!U2!N Mpc!U-2!N)'
     endif else begin
        if keyword_set(beam2) then plotfile = plotfile + '_beam2'
        title = 'P!Dk!N (mK!U2!N Mpc!U-2!N)'
     endelse
     if keyword_set(lsp) then begin
        plotfile = plotfile + '_lsp' 
        wt_plotfile = wt_plotfile + '_lsp'
     endif

     wt_title = 'Weight!U2!N'

     kpower_slice_plot, reform(power, n_kx, n_ky, 1), kx_mpc, ky_mpc, [0], slice_axis = 2, slice_ind = 0, $
     title = title, pub = pub, plotfile = plotfile

     kpower_slice_plot, reform(abs(weight_ft)^2, n_kx, n_ky, 1), kx_mpc, ky_mpc, [0], slice_axis = 2, slice_ind = 0, $
     title = wt_title, pub = pub, plotfile = wt_plotfile, window_num = 2

  endif

 
end
