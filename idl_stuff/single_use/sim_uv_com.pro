pro sim_uv_com, t32 = t32, define_baselines = define_baselines, baseline_spacing = baseline_spacing, $
                baseline_layout = baseline_layout, fine_res = fine_res, use_outliers = use_outliers, curvature = curvature

  if keyword_set(define_baselines) then begin
     if n_elements(baseline_spacing) eq 0 then baseline_spacing = 4
     if n_elements(baseline_layout) eq 0 then baseline_layout = 'simple'

     spacing_name = strtrim_plus(string(baseline_spacing), flag = 2)
     tag = baseline_layout + spacing_name + '_'
     filename = 'Simulation_' + baseline_layout + spacing_name

  endif else if keyword_set(t32) then tag = '32t_' else begin
     filename = 'Simulation'

     if keyword_set(use_outliers) then begin
        filename = filename + '_512t'
        tag = '512t_' 
     endif else begin
        tag = '496t_'
        filename = filename + '_496t'
     endelse
     if keyword_set(fine_res) then begin
        tag = tag + 'fine_'
        filename = filename + '_fine'
     endif
  endelse

  froot = base_path('data') + 'fhd_simulations_old/'
  info_file = froot + 'sim_' + tag + 'info.idlsave'
  weights_file = froot + 'sim_' + tag + 'weights.idlsave'

  ;;restore, info_file
  restore, weights_file
  if max(abs(imaginary(weights_cube))) eq 0 then weights_cube = real_part(weights_cube) else stop

  dims = size(weights_cube, /dimension)
  temp = n_elements(temporary(tile_beam_image))

  if keyword_set(curvature) then begin
     type_name = 'curv'

     x_arr = rebin(((dindgen(dims[0])-dims[0]/2)*2d/dims[0])^2d, dims[0], dims[1])
     y_arr = rebin(reform(((dindgen(dims[1])-dims[1]/2)*2d/dims[1])^2d, 1, dims[1]), dims[0], dims[1])
     uv_x_in = x_arr ;; real ramp -- has same symmetry as real (not imaginary b/ complex conjugation)
     uv_y_in = y_arr

  endif else begin
     type_name = 'com'

     x_arr = rebin((dindgen(dims[0])-dims[0]/2)*2d/dims[0], dims[0], dims[1])
     y_arr = rebin(reform((dindgen(dims[1])-dims[1]/2)*2d/dims[1], 1, dims[1]), dims[0], dims[1])
     uv_x_in = complex(0,1) * x_arr ;; imaginary ramp -- same symmetry with avoid complex conjugation
     uv_y_in = complex(0,1) * y_arr
 endelse

  sigma2_cube = 1d/(weights_cube)
  wh = where(weights_cube eq 0, count)
  if count ne 0 then sigma2_cube[wh] = 0

  uv_x_com = sim_xfer_apply(uv_x_in, t32 = t32, define_baselines = define_baselines, baseline_spacing = baseline_spacing, $
                           baseline_layout = baseline_layout, fine_res = fine_res, use_outliers = use_outliers)
  if max(abs(imaginary(uv_x_com))) eq 0 then uv_x_com = real_part(temporary(uv_x_com))
  if max(abs(imaginary(uv_x_in))) eq 0 then uv_x_in = real_part(temporary(uv_x_in))
  if max(abs(real_part(uv_x_com))) eq 0 then uv_x_com = imaginary(temporary(uv_x_com))
  if max(abs(real_part(uv_x_in))) eq 0 then uv_x_in = imaginary(temporary(uv_x_in))
  uv_x_com = uv_x_com * sigma2_cube
  if count ne 0 then uv_x_com[wh] = 0

  x_savefile = froot + 'sim_' + tag + 'x'+ type_name + '.idlsave'
  save, file = x_savefile, uv_x_in, uv_x_com
  undefine, uv_x_com
  undefine, uv_x_in

  uv_y_com = sim_xfer_apply(uv_y_in, t32 = t32, define_baselines = define_baselines, baseline_spacing = baseline_spacing, $
                            baseline_layout = baseline_layout, fine_res = fine_res, use_outliers = use_outliers)
  if max(abs(imaginary(uv_y_com))) eq 0 then uv_y_com = real_part(uv_y_com)
  if max(abs(imaginary(uv_y_in))) eq 0 then uv_y_in = real_part(uv_y_in)
  if max(abs(real_part(uv_y_com))) eq 0 then uv_y_com = imaginary(temporary(uv_y_com))
  if max(abs(real_part(uv_y_in))) eq 0 then uv_y_in = imaginary(temporary(uv_y_in))
  uv_y_com = uv_y_com * sigma2_cube
  if count ne 0 then uv_y_com[wh] = 0

  y_savefile = froot + 'sim_' + tag + 'y'+ type_name + '.idlsave'
  save, file = y_savefile, uv_y_in, uv_y_com

end
