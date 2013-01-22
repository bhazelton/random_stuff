pro plot_freq, xcoord, ycoord, model_cube, sim_cube, type = type
  
  if n_elements(type) eq 0 then type = 'phase'
  type_enum = ['amp', 'phase', 'real', 'imaginary']
  wh_type = where(type_enum eq type, count)
  if count eq 0 then message, 'type not recognized. Type must be one of: ' + type_enum

  dims = size(model_cube,/dimension)

  case type of
     'amp' : begin
        model_line = abs(model_cube[xcoord, ycoord, *])
        sim_line = abs(sim_cube[xcoord, ycoord,*])
     end
     'phase': begin
        model_line = atan(model_cube[xcoord, ycoord, *],/phase)
        sim_line = atan(sim_cube[xcoord, ycoord,*],/phase)
     end
     'real': begin
        model_line = real_part(model_cube[xcoord, ycoord, *])
        sim_line = real_part(sim_cube[xcoord, ycoord,*])
     end
     'imaginary': begin
        model_line = imaginary(model_cube[xcoord, ycoord, *])
        sim_line = imaginary(sim_cube[xcoord, ycoord,*])
     end
  endcase



  yrange = minmax([model_line, sim_line])
  xrange = [0, dims[2]-1]

  plot, model_line, xrange = xrange, yrange = yrange, xstyle=1, ystyle=1
  oplot, sim_line, color = 254

end


pro sim_test, sim_num = sim_num, t32 = t32, baseline_layout = baseline_layout, baseline_spacing = baseline_spacing, $
              use_outliers = use_outliers, plot_diff = plot_diff, clean_type = clean_type, off_axis = off_axis, $
              full_sky = full_sky, source_radius = source_radius

  froot = base_path('data') + 'fhd_simulations_old/'

  n_space = n_elements(baseline_spacing)
  n_layout = n_elements(baseline_layout)

  if n_elements(sim_num) eq 0 then sim_num = 0

  if max([n_space, n_layout, n_elements(sim_num)]) le 1 then n_sims = 1 $
  else if max([n_space, n_layout, n_elements(sim_num)]) gt 2 then message, 'sim_test can only compare 2 simulations' $
  else begin
     n_sims = 2
     if n_space eq 1 then baseline_spacing = [baseline_spacing, baseline_spacing]
     if n_layout eq 1 then baseline_layout = [baseline_layout, baseline_layout]
  endelse     
  

  if n_elements(baseline_layout) gt 0 or n_elements(baseline_spacing) gt 0 then begin
     if n_elements(baseline_spacing) eq 0 then if n_sims eq 1 then baseline_spacing = 4 else baseline_spacing = fltarr(2)+4
     if n_elements(baseline_layout) eq 0 then if n_sims eq 1 then baseline_layout = 'simple' $
     else baseline_layout = strarr(2) + 'simple'
     
     names = string([3, 6, 13, 25, 50, 75, 100, 125, 250, 500], format = '(i04)')
     ;; offsets = [3.16, 6.45, 12.65, 25.09, 49.91, 75.3, 100.04, 124.58, 250.26, 499.18]
     
     array = strarr(n_sims)
     for i=0, n_sims-1 do begin
        
        spacing_name = number_formatter(string(baseline_spacing[i]), flag=2)
        array[i] = baseline_layout[i] + spacing_name
     endfor
  endif else if keyword_set(t32) then begin
     names = string([3, 6, 13, 25, 50, 75, 100, 125, 250, 500], format = '(i04)')
     ;; offsets = [3.16, 6.45, 12.65, 25.09, 49.91, 75.3, 100.04, 124.58, 250.26, 499.18]
     if n_sims eq 1 then array = '32t' else array = strarr(2) + '32t'
  endif else begin
     ;; offsets = [25.3, 51.6, 101.2, 200.7, 399.3, 602.4, 800.3, 998.8]
     names = string([25, 50, 100, 200, 400, 600, 800, 1000], format = '(i04)')
     
     if keyword_set(use_outliers) then array = '512t' else array = '496t'
     if n_sims eq 2 then array = starr(2) + array
  endelse
  
  if keyword_set(full_sky) then begin
     if n_elements(source_radius) eq 0 then source_radius = 3
     names = 'fullsky_' + number_formatter(source_radius) + 'deg'
   endif else begin
     if keyword_set(off_axis) then names = 'offaxis'
  endelse


  fbase_arr = 'sim_' + array + '_' + names[sim_num] + '_uvf'
  sim_file = froot + fbase_arr + '.idlsave'
  weights_file = froot + 'sim_' + array + '_weights.idlsave'
  
  ;; setup first simulation
  restore, sim_file[0]
  uvf_cube1 = temporary(uvf_cube)
  
  dims = size(uvf_cube1, /dimensions)
  n_kx = dims[0]
  n_ky = dims[1]
  n_kz = dims[2]

  restore, weights_file[0]
  weights_cube1 = temporary(weights_cube)
  ;;temp = n_elements(temporary(tile_beam_image))

  sigma2_cube1 = 1d/(weights_cube1)
  wh1 = where(weights_cube1 eq 0, count1)
  if count1 ne 0 then sigma2_cube1[wh1] = 0

  sim_cube1 = uvf_cube1 * sigma2_cube1
  if count1 ne 0 then sim_cube1[wh1] = 0
  
  undefine, uvf_cube1
  undefine, weights_cube1
  undefine, sigma2_cube1

  ;; get input cube
  ;; columns are: x location (zero centered), y location (zero centered), flux
  my_sources1 = reform(my_sources, 3, n_elements(my_sources)/3)
  n_sources = (size(my_sources1, /dimension))[1]
  my_source_array1 = reform(fltarr(7, n_sources), 7, n_sources)
  my_source_array1[0,*]=my_sources1[0,*] + n_kx/2.
  my_source_array1[1,*]=my_sources1[1,*] + n_ky/2.
  my_source_array1[4,*]=my_sources1[2,*] + n_ky/2.
  my_source_array1[5,*]=Round(my_source_array1[0,*])+Round(my_source_array1[1,*])*n_kx

  model_uv1 = visibility_source_uv_grid(my_source_array1, u_dim = n_kx, v_dim = n_ky)
  model_uvf1 = make_array(n_kx, n_ky,  n_kz, /dcomplex)
  model_uvf1[*] = model_uv1[*] # (1.0 + dblarr(n_kz))


  ;; need to cut uvf cubes in half (in u or v -- we'll cut in v) before proceeding
  sim_cube1_clip = sim_cube1[*, n_ky/2:n_ky-1,*]

  sim_ft1 = fft(sim_cube1_clip, dimension=3)
  sim_ft1 = shift(sim_ft1, [0,0,n_kz/2])

  wset, 1
  cgimage, alog(abs(sim_ft1[*,*,0])),/scale
  wset, 2
  cgimage, alog(abs(sim_ft1[*,0,*])),/scale

  stop

  if windowavailable(1) then wset, 1 else window, 1
  if keyword_set(plot_diff) then begin
     model_uvf1 = model_uvf1/ max(abs(model_uvf1)) * max(abs(sim_cube1))
     temp1 = atan(model_uvf1[*,512,*],/phase) - atan(sim_cube1[*,512,*],/phase)
     thresh1 = stdev(temp1)
     wh_low1 = where(temp1 lt -1 * thresh1)
     wh_high1 = where(temp1 gt thresh1)
     temp1[wh_low1] = -1*thresh1
     temp1[wh_high1] = thresh1
     
     cgimage, temp1
  endif else cgimage, atan(sim_cube1[*,512,*],/phase)

  if n_sims eq 2 then begin
     restore, sim_file[1]
     uvf_cube2 = temporary(uvf_cube)
     
     restore, weights_file[1]
     weights_cube2 = temporary(weights_cube)
     
     sigma2_cube2 = 1d/(weights_cube2)
     wh2 = where(weights_cube2 eq 0, count2)
     if count2 ne 0 then sigma2_cube2[wh2] = 0
     
     sim_cube2 = uvf_cube2 * sigma2_cube2
     if count2 ne 0 then sim_cube2[wh2] = 0
     
     undefine, uvf_cube1
     undefine, weights_cube1
     undefine, sigma2_cube1

     dims = size(uvf_cube2, /dimensions)
     n_kx = dims[0]
     n_ky = dims[1]
     n_kz = dims[2]
     
     ;; get input cube
     ;; columns are: x location (zero centered), y location (zero centered), flux
     my_sources2 = reform(my_sources, 3, n_elements(my_sources)/3)
     n_sources = (size(my_sources2, /dimension))[1]
     my_source_array2 = reform(fltarr(7, n_sources), 7, n_sources)
     my_source_array2[0,*]=my_sources2[0,*] + n_kx/2.
     my_source_array2[1,*]=my_sources2[1,*] + n_ky/2.
     my_source_array2[4,*]=my_sources2[2,*] + n_ky/2.
     my_source_array2[5,*]=Round(my_source_array2[0,*])+Round(my_source_array2[1,*])*n_kx
     
     model_uv2 = visibility_source_uv_grid(my_source_array2, u_dim = n_kx, v_dim = n_ky)
     model_uvf2 = make_array(n_kx, n_ky,  n_kz, /dcomplex)
     model_uvf2[*] = model_uv2[*] # (1.0 + dblarr(n_kz))
     
     ;; need to cut uvf cubes in half (in u or v -- we'll cut in v) before proceeding
     sim_cube2_clip = sim_cube2[*, n_ky/2:n_ky-1,*]
     
     sim_ft2 = fft(sim_cube2_clip, dimension=3)
     sim_ft2 = shift(sim_ft2, [0,0,n_kz/2])
     
     if windowavailable(2) then wset, 2 else window, 2
     if keyword_set(plot_diff) then begin
        model_uvf2 = model_uvf2/ max(abs(model_uvf2)) * max(abs(uvf_cube2))
        temp2 = atan(model_uvf2[*,512,*],/phase) - atan(sim_cube2[*,512,*],/phase)
        thresh2 = stdev(temp2)
        wh_low2 = where(temp2 lt -1 * thresh2)
        wh_high2 = where(temp2 gt thresh2)
        temp2[wh_low2] = -1*thresh2
        temp2[wh_high2] = thresh2
        
        cgimage, temp2
     endif else cgimage, atan(sim_cube2[*,512,*],/phase)
 endif
  
end
