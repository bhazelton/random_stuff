

pro sim_plots, sim_num = sim_num, no_kzero = no_kzero, pub = pub, refresh_ps = refresh_ps, refresh_binning = refresh_binning, $
               plot_slice = plot_slice, std_power = std_power, beam_exp = beam_exp, no_weighting = no_weighting, $
               add_noise = add_noise, t32 = t32, baseline_layout = baseline_layout, baseline_spacing = baseline_spacing, $
               fill_holes = fill_holes, data_range = data_range, compare_1d = compare_1d, psf = psf, plot_weights = plot_weights, $
               cont_data_lims = cont_data_lims, eor_only = eor_only, test_power_shape = test_power_shape, eor_test = eor_test, $
               norm_2d = norm_2d, grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, slice_nobin = slice_nobin, $
               undefined_data_range = undefined_data_range, slice_linear_axes = slice_linear_axes, use_outliers = use_outliers, $
               baseline_axis = baseline_axis, delay_axis = delay_axis, plot_uvf = plot_uvf, plot_urange = plot_urange, $
               plot_vrange = plot_vrange, uvf_data_range = uvf_data_range, uvf_conv = uvf_conv, uvf_type = uvf_type, $
               clean_type = clean_type, clean_ratio = clean_ratio, delta = delta, eor_ratio = eor_ratio, log_kpar = log_kpar, $
               log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, $
               off_axis = off_axis, full_sky = full_sky, source_radius = source_radius

  if n_elements(beam_exp) eq 0 then beam_exp = 0
  if keyword_set(eor_test) then eor_only = 1
  if keyword_set(eor_ratio) then begin
     if keyword_set(clean_ratio) then message, 'both eor_ratio and clean_ratio cannot be set.'
     if keyword_set(eor_test) then message, 'Cannot set both eor_ratio and eor_test.'
     eor_only = 1
  endif

  ;; default to absolute value for uvf plots
  if keyword_set(plot_uvf) and n_elements(uvf_type) eq 0 then begin
     if keyword_set(full_sky) then uvf_type = 'abs' else uvf_type = 'phase'
  endif

  ;; default to including baseline axis
  if n_elements(baseline_axis) eq 0 then baseline_axis = 1
  ;; default to including delay axis
  if n_elements(delay_axis) eq 0 then delay_axis = 1

  if keyword_set(clean_ratio) and n_elements(clean_type) eq 0 then clean_type = 'hmf'
  clean_type_enum = ['hmf', 'iterate', 'fit']
  if n_elements(clean_type) ne 0 then begin
     wh = where(clean_type_enum eq clean_type, count)
     if count eq 0 then message, 'Clean type not recognized'

     if beam_exp ne 0 then message, 'Cleaning is not compatible with primary beam removal'
     if keyword_set(eor_test) then message, 'Cleaning is not compatible with eor_test'
  endif

  froot = base_path('data') + 'fhd_simulations_old/'

  if keyword_set(baseline_layout) or keyword_set(baseline_spacing) then begin
     if n_elements(baseline_spacing) eq 0 then baseline_spacing = 4
     if n_elements(baseline_layout) eq 0 then baseline_layout = 'simple'

     names = string([3, 6, 13, 25, 50, 75, 100, 125, 250, 500], format = '(i04)')
     offsets = [3.16, 6.45, 12.65, 25.09, 49.91, 75.3, 100.04, 124.58, 250.26, 499.18]
  
     spacing_name = number_formatter(string(baseline_spacing), flag = 2)
     array = baseline_layout + spacing_name
     tile_tag = '_' + array

     plot_folder = 'sim_simple_baselines/'
  endif else if keyword_set(t32) then begin
     names = string([3, 6, 13, 25, 50, 75, 100, 125, 250, 500], format = '(i04)')
     offsets = [3.16, 6.45, 12.65, 25.09, 49.91, 75.3, 100.04, 124.58, 250.26, 499.18]
     array = '32t'
     tile_tag = '_32t'

     plot_folder = 'sim_32t/'
  endif else begin
     offsets = [25.3, 51.6, 101.2, 200.7, 399.3, 602.4, 800.3, 998.8]
     names = string([25, 50, 100, 200, 400, 600, 800, 1000], format = '(i04)')
 
     ;; max_good = 3
     ;; offsets = offsets[0:max_good]
     ;; names = names[0:max_good]

     if keyword_set(use_outliers) then array = '512t' else array = '496t'
     tile_tag = '_' + array

     plot_folder = 'sim_' + array + '/'
  endelse

  info_file = froot + 'sim_' + array + '_info.idlsave'
  restore, info_file   

  sim_plot_folder = ''
  if keyword_set(full_sky) then begin
     sim_plot_folder = 'fullsky/'

     if n_elements(source_radius) eq 0 then begin
        names = 'fullsky'
        if keyword_set(t32) then xy_length = 1024 else xy_length = 2048
        deg_offsets = xy_length * degpix / sqrt(2d)
        deg_offset_str = '<' + number_formatter(deg_offsets, format = '(g4.2)')
     endif else begin
        names = 'fullsky_' + number_formatter(source_radius) + 'deg'
        deg_offsets = source_radius
        deg_offset_str = '<' + number_formatter(deg_offsets, format = '(g4.2)')
     endelse
  endif else begin
     if keyword_set(off_axis) then begin
        y_offsets = offsets[1]
        x_offsets = offsets[0]
        names = 'offaxis'
     endif else begin
        x_offsets = offsets
        y_offsets = fltarr(n_elements(offsets))
     endelse

     deg_offsets = sqrt(x_offsets^2d + y_offsets^2d) * degpix
     deg_offset_str = number_formatter(deg_offsets, format = '(g4.2)')
  endelse

  rad_offsets =  deg_offsets *!pi/180d

  if n_elements(sim_num) eq 0 then sim_num = 0
  if min(sim_num) lt 0 or max(sim_num) gt n_elements(names)-1 then $
     message, 'sim num must be between 0 and ' + number_formatter(n_elements(names)-1)
  n_sims = n_elements(sim_num)  
  
  fadd = ''
  eor_test_fadd = ''
  if keyword_set(add_noise) then fadd = fadd + '_noise'
  if n_elements(clean_type) ne 0 then begin
     case clean_type of
        'hmf': clean_tag = '_clean'
        'iterate': clean_tag = '_cleanplus'
        'fit': clean_tag = '_cleanfit'
     endcase
     fadd = fadd + clean_tag
  endif
  case beam_exp of
     0: fadd = fadd + '_holo'
     1: fadd = fadd + '_beam'
     2: fadd = fadd + '_truesky'
     else: message, 'division by more than 2 factors of the primary beam is unsupported.'
  endcase
  if keyword_set(std_power) then begin
     fadd = fadd + '_sp'
     eor_test_fadd = eor_test_fadd  + '_sp'
  endif
  fadd_3d = fadd
  if keyword_set(no_weighting) then begin
     fadd = fadd + '_nowt'
     eor_test_fadd = eor_test_fadd  + '_nowt'
  endif
  if keyword_set(no_kzero) then begin
     fadd = fadd + '_nok0'
     eor_test_fadd = eor_test_fadd  + '_nok0'
  endif

  fadd_2dbin = ''
  if keyword_set(fill_holes) then fadd_2dbin = fadd_2dbin + '_nohole'
  if keyword_set(log_kpar) then fadd_2d = fadd_2d + '_logkpar'
  if keyword_set(log_kperp) then fadd_2d = fadd_2d + '_logkperp'

  fbase_arr = 'sim' + tile_tag + '_' + names[sim_num] + '_uvf'
  sim_file = froot + fbase_arr + '.idlsave'

  savefile_3d = froot + fbase_arr + fadd_3d + '_power.idlsave'
  savefiles_2d = froot + fbase_arr + fadd + fadd_2dbin + '_2dkpower.idlsave'
  savefiles_1d = froot + fbase_arr + fadd +'_1dkpower.idlsave'
  eor_file_1d = base_path('data') + 'eor_data/eor_power_1d.idlsave'
  eor_file_1d_input = froot + 'eortest' + tile_tag + '_uvf' + eor_test_fadd +'_1dkpower.idlsave'

  if keyword_set(clean_ratio) then begin
     
    case clean_type of
        'hmf': title_str = 'Cleaned'
        'iterate': title_str = 'Cleaned+'
        'fit': title_str = 'Cleaned Fit'
     endcase

     reg_fadd = ''
     if keyword_set(add_noise) then reg_fadd = reg_fadd + '_noise'
     reg_fadd = reg_fadd + '_holo'
     if keyword_set(std_power) then reg_fadd = reg_fadd + '_sp'
     reg_fadd_3d = reg_fadd
     if keyword_set(no_weighting) then reg_fadd = reg_fadd + '_nowt'
     if keyword_set(no_kzero) then reg_fadd = reg_fadd + '_nok0'
 
     reg_fadd_2dbin = ''
     if keyword_set(fill_holes) then reg_fadd_2dbin = reg_fadd_2dbin + '_nohole'
     if keyword_set(log_kpar) then fadd_2d = fadd_2d + '_logkpar'
     if keyword_set(log_kperp) then fadd_2d = fadd_2d + '_logkperp'

     reg_savefile_3d = froot + fbase_arr + reg_fadd_3d + '_power.idlsave'
     reg_savefiles_2d = froot + fbase_arr + reg_fadd + reg_fadd_2dbin + '_2dkpower.idlsave'

     ftests_reg = file_test(reg_savefiles_2d) *  (1 - file_test(reg_savefiles_2d, /zero_length))

     compare_1d = 1
  endif 

  if keyword_set(compare_1d) then begin
     if keyword_set(clean_ratio) then begin
        ;; if doing clean_ratio, then compare clean & regular
        comp_add = reg_fadd
     endif else begin
        ;; otherwise compare standard, no weight & regular
        comp_add = ''
        if keyword_set(add_noise) then comp_fadd = comp_fadd + '_noise'
        if n_elements(clean_type) ne 0 then comp_fadd = comp_fadd + clean_tag
        if beam_exp eq 1 then comp_fadd = comp_fadd + '_beam' else if beam_exp eq 0 then comp_fadd = comp_fadd + '_holo'
        case beam_exp of
           0: comp_fadd = comp_fadd + '_holo'
           1: comp_fadd = comp_fadd + '_beam'
           2: comp_fadd = comp_fadd + '_truesky'
           else: message, 'division by more than 2 factors of the primary beam is unsupported.'
        endcase
        comp_add = comp_add + '_sp_nowt'
        if keyword_set(no_kzero) then comp_add = comp_add + '_nok0'
     endelse
     savefiles_1d_compare = froot + fbase_arr + comp_add + '_1dkpower.idlsave'
     ftests_1d_comp = file_test(savefiles_1d_compare) *  (1 - file_test(savefiles_1d_compare, /zero_length))
  endif

  ftests_2d = file_test(savefiles_2d) *  (1 - file_test(savefiles_2d, /zero_length))
  ftests_1d = file_test(savefiles_1d) *  (1 - file_test(savefiles_1d, /zero_length))
  ftests = ftests_1d * ftests_2d

  for i=0, n_sims-1 do begin
     if ftests[i] eq 0 or keyword_set(refresh_ps) or keyword_set(refresh_binning) then sim_3dps, sim_file[i], $
        refresh = refresh_ps, no_kzero = no_kzero, beam_exp = beam_exp, std_power = std_power, no_weighting = no_weighting, $
        add_noise = add_noise, clean_type = clean_type, fill_holes = fill_holes, log_kpar = log_kpar, log_kperp = log_kperp, $
        kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, /quiet

     if keyword_set(clean_ratio) then begin
        ;; make regular file if it's missing
        if ftests_reg[i] eq 0 or keyword_set(refresh_ps) or keyword_set(refresh_binning) then sim_3dps, sim_file[i], $
           refresh = refresh_ps, no_kzero = no_kzero, beam_exp = beam_exp, std_power = std_power, no_weighting = no_weighting, $
           add_noise = add_noise, fill_holes = fill_holes, log_kpar = log_kpar, log_kperp = log_kperp, $
           kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, /quiet
     endif else if keyword_set(compare_1d) then begin
        ;; make comparision file if it's missing
        if ftests_1d_comp[i] eq 0 or keyword_set(refresh_ps) or keyword_set(refresh_binning) then sim_3dps, sim_file[i], $
           refresh = refresh_ps, beam_exp = beam_exp, /std_power, /no_weighting, no_kzero = no_kzero, add_noise = add_noise, $
           clean_type = clean_type, fill_holes = fill_holes, log_kpar = log_kpar, log_kperp = log_kperp, $
           kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, /quiet
     endif
  endfor

  ;; test for input 1D EOR file
  ftests_input_eor = file_test(eor_file_1d_input) *  (1 - file_test(eor_file_1d_input, /zero_length))
  if ftests_input_eor eq 0 then $
     sim_3dps, sim_file, /eor_test, no_kzero = no_kzero, std_power = std_power, no_weighting = no_weighting, fill_holes = fill_holes, $
               log_kpar = log_kpar, log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, $
               k1d_bin = k1d_bin, /quiet


  if keyword_set(psf) then begin
     psf_fbase = 'psf' + tile_tag + '_uvf'
     savefile_psf_base = froot + psf_fbase + fadd
     savefile_psf = savefile_psf_base +'_2dkpower.idlsave'
     psf_ftest = file_test(savefile_psf) *  (1 - file_test(savefile_psf, /zero_length))

     if psf_ftest eq 0 or keyword_set(refresh_ps) or keyword_set(refresh_binning) then sim_3dps, sim_file, /psf, $
        refresh = refresh_ps,  no_kzero = no_kzero, beam_exp = beam_exp, std_power = std_power, no_weighting = no_weighting, $
        add_noise = add_noise, clean_type = clean_type, fill_holes = fill_holes, log_kpar = log_kpar, log_kperp = log_kperp, $
        kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, /quiet
  endif

  if keyword_set(eor_only) then begin
     if keyword_set(eor_test) then eor_tag = 'test' else eor_tag = ''
     eor_fbase = 'eor' + eor_tag + tile_tag + '_uvf'
     
     savefile_eor_base = froot + eor_fbase
     if keyword_set(eor_test) then savefile_eor_base = savefile_eor_base + eor_test_fadd $
     else savefile_eor_base = savefile_eor_base + fadd
     savefile_eor = savefile_eor_base + fadd_2dbin + '_2dkpower.idlsave'
     savefile_eor_1d = savefile_eor_base +'_1dkpower.idlsave'
     
     eor_ftest_2d = file_test(savefile_eor) *  (1 - file_test(savefile_eor, /zero_length))
     eor_ftest_1d = file_test(savefile_eor_1d) *  (1 - file_test(savefile_eor_1d, /zero_length))
     eor_ftest = eor_ftest_2d * eor_ftest_1d
    
     if eor_ftest eq 0 or keyword_set(refresh_ps) or keyword_set(refresh_binning) then sim_3dps, sim_file, /eor_only, $
        refresh = refresh_ps,  no_kzero = no_kzero, beam_exp = beam_exp, std_power = std_power, no_weighting = no_weighting, $
        add_noise = add_noise, fill_holes = fill_holes, clean_type = clean_type, log_kpar = log_kpar, log_kperp = log_kperp, $
        kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, /quiet, eor_test = eor_test
  endif

  if n_elements(test_power_shape) ne 0 then begin
      test_power_fbase = test_power_shape + 'power' + tile_tag + '_uvf'
 
      savefile_test_power_base = froot + test_power_fbase + fadd
      savefile_test_power = savefile_test_power_base + fadd_2dbin +'_2dkpower.idlsave'
      savefile_test_power_1d = savefile_test_power_base +'_1dkpower.idlsave'
    
      test_power_ftest_2d = file_test(savefile_test_power) *  (1 - file_test(savefile_test_power, /zero_length))
      test_power_ftest_1d = file_test(savefile_test_power_1d) *  (1 - file_test(savefile_test_power_1d, /zero_length))
      test_power_ftest = test_power_ftest_2d * test_power_ftest_1d
     
     if test_power_ftest eq 0 or keyword_set(refresh_ps) or keyword_set(refresh_binning) then sim_3dps, sim_file, $
        refresh = refresh_ps,  no_kzero = no_kzero, beam_exp = beam_exp, std_power = std_power, no_weighting = no_weighting, $
        add_noise = add_noise, fill_holes = fill_holes, clean_type = clean_type, log_kpar = log_kpar, log_kperp = log_kperp, $
        kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, /quiet, test_power_shape = test_power_shape
  endif

  restore, savefiles_2d[0]
  if keyword_set(clean_ratio) and n_elements(data_range) eq 0 and not keyword_set(undefined_data_range) then begin
     if keyword_set(full_sky) then data_range = [1e-3, 1e0] else data_range = [1e-6, 1e0]
  endif else if keyword_set(eor_ratio) and n_elements(data_range) eq 0 and not keyword_set(undefined_data_range) then $
     data_range = [1e-14, 1.001e-8]
  if keyword_set(baseline_layout) then begin
     case baseline_layout of
        'simple': begin
           if n_elements(kperp_plot_range) eq 0 then kperp_plot_range = [2.5e-3, 0.08] $
           else kperp_plot_range = [max([kperp_plot_range[0], min(kperp_edges)]), min([kperp_plot_range[1], max(kperp_edges)])]
           if n_elements(data_range) eq 0 and not keyword_set(undefined_data_range) then begin
              if keyword_set(norm_2d) then data_range = [1e-10, 1e0] else $
                 if beam_exp lt 2 then data_range = [1e18, 1e23] $ 
                 else data_range = [1e50, 1e56]
           endif
              if not keyword_set(undefined_data_range) then data_range_1d = [1e3, 1e7]
           end
        'single': begin
           if n_elements(kperp_plot_range) eq 0 then kperp_plot_range = [2.5e-3, 0.08] $
           else kperp_plot_range = [max([kperp_plot_range[0], min(kperp_edges)]), min([kperp_plot_range[1], max(kperp_edges)])]
           if n_elements(data_range) eq 0 and not keyword_set(undefined_data_range) then begin
              if keyword_set(norm_2d) then data_range = [1e-10, 1e0] else $
                 if beam_exp lt 2 then data_range = [1e18, 1e23] $ 
                 else data_range = [1e50, 1e56]
           endif
           if not keyword_set(undefined_data_range) then data_range_1d = [1e3, 1e7]

        end
     endcase
  endif else begin
     case array of
        '32t': begin
           if n_elements(kperp_plot_range) eq 0 then kperp_plot_range = [min(kperp_edges), 0.005] $
           else kperp_plot_range = [max([kperp_plot_range[0], min(kperp_edges)]), min([kperp_plot_range[1], max(kperp_edges)])]
           if n_elements(data_range) eq 0 and not keyword_set(undefined_data_range) then begin
              if keyword_set(norm_2d) then data_range = [1e-10, 1e0] else $
                 data_range = [1e8, 1e18]
           endif
           if not keyword_set(undefined_data_range) then data_range_1d = [1e2, 1e6]
        end
        '512t': begin
           if keyword_set(cont_data_lims) then begin
              kperp_plot_range = [3e-3, 0.15]
              kpar_plot_range = [2e-2, max(kpar_edges)]
           endif
           
           if n_elements(kperp_plot_range) eq 0 then kperp_plot_range = [min(kperp_edges), 0.3] $
           else kperp_plot_range = [max([kperp_plot_range[0], min(kperp_edges)]), min([kperp_plot_range[1], max(kperp_edges)])]
           if n_elements(data_range) eq 0 and not keyword_set(undefined_data_range) then begin
              if keyword_set(norm_2d) then data_range = [1e-7, 1e0] else $
                 data_range = [1e21, 1e25]
           endif
           if not keyword_set(undefined_data_range) then if beam_exp lt 2 then data_range_1d = [1e-2, 1e11] $
           else data_range_1d = [1e-3, 1e16]
        end
        '496t': begin
          if keyword_set(cont_data_lims) then begin
              kperp_plot_range = [3e-3, 0.15]
              kpar_plot_range = [2e-2, max(kpar_edges)]
           endif
           
           if n_elements(kperp_plot_range) eq 0 then kperp_plot_range = [5d/kperp_lambda_conv, 0.3] $
           else kperp_plot_range = [max([kperp_plot_range[0], min(kperp_edges)]), min([kperp_plot_range[1], max(kperp_edges)])]
           if n_elements(data_range) eq 0 and not keyword_set(undefined_data_range) then begin
              if keyword_set(norm_2d) then data_range = [1e-7, 1e0] else $
                 if keyword_set(full_sky) then data_range = [1e11, 1e18] else data_range = [1e7, 1e15]
           endif
           if not keyword_set(undefined_data_range) then if keyword_set(delta) then data_range_1d = [1e-3, 1e6] $
           else data_range_1d = [1e-3, 1e16]
 
           if keyword_set(plot_uvf) and n_elements(uvf_data_range) eq 0 and not keyword_set(undefined_data_range) then begin
              if uvf_type eq 'abs' then begin
                 if keyword_set(full_sky) then uvf_data_range = [0, 5.5e7]
              endif
           endif
        end
     endcase
  endelse
  if not keyword_set(undefined_data_range) then if keyword_set(slice_nobin) then slice_data_range = data_range * [1e-1,1e1] $
  else slice_data_range = data_range

  plot_fadd = ''
  if keyword_set(cont_data_lims) then plot_fadd = plot_fadd + '_contig'
  if keyword_set(norm_2d) then plot_fadd = plot_fadd + '_norm'
  if keyword_set(grey_scale) then plot_fadd = plot_fadd + '_grey'

  plotfile_path = base_path('plots') + 'power_spectrum/' + plot_folder
  if n_sims eq 1 then plotfile_base = plotfile_path + sim_plot_folder + 'sim' + tile_tag + '_' + names[sim_num] + fadd $
  else plotfile_base = plotfile_path + sim_plot_folder + 'sim' + tile_tag + '_' + strjoin(string(sim_num, format='(i1)')) + fadd
  if keyword_set(clean_ratio) then plotfile_2d = plotfile_base + fadd_2dbin + clean_tag + '_ratio' + plot_fadd + '.eps' $
  else if keyword_set(eor_ratio) then plotfile_2d = plotfile_base + fadd_2dbin + '_eor' + '_ratio' + plot_fadd + '.eps' $
  else plotfile_2d = plotfile_base + fadd_2dbin + '_kspace_power' + plot_fadd + '.eps'
  
  if keyword_set(delta) then plotfile_1d = plotfile_base + '_1d_delta.eps' else plotfile_1d = plotfile_base + '_1d_power.eps'

  if not keyword_set(slice_nobin) then slice_fadd = '_binned' else slice_fadd = ''
  yslice_plotfile = plotfile_base + '_xz_plane' + plot_fadd + slice_fadd + '.eps'
  xslice_plotfile = plotfile_base + '_yz_plane' + plot_fadd + slice_fadd + '.eps'
  zslice_plotfile = plotfile_base + '_xy_plane' + plot_fadd + '.eps' ;; zslice is never rebinnined

  if keyword_set(plot_wedge_line) then begin
     redshift = 8
     cosmology_measures, redshift, wedge_factor = wedge_factor
     source_dists = rad_offsets[sim_num]
     wedge_amp = wedge_factor * source_dists
  endif else wedge_amp = rad_offsets*0.


  if keyword_set(psf) then begin
     psf_plotfile = plotfile_path + 'psf/' + 'psf' + tile_tag + fadd + fadd_2dbin + '_kspace_power'
     if not keyword_set(undefined_data_range) then if keyword_set(norm_2d) then psf_data_range = [1e-10, 1e0] $
     else psf_data_range = [1e-17, 1e3]
     kpower_2d_plots, savefile_psf, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = psf_data_range, pub = pub, plotfile = psf_plotfile, title = 'PSF', window_num = 9, $
                      norm_2d = norm_2d, grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, $
                      baseline_axis = baseline_axis, delay_axis = delay_axis
                     
  endif

  if keyword_set(eor_only) then begin
     eor_plotfile = plotfile_path + 'eor/' + 'eor' + eor_tag + '_sim' + tile_tag + fadd + fadd_2dbin + '_kspace_power'
     if not keyword_set(undefined_data_range) then if keyword_set(norm_2d) then eor_data_range = [1e-10, 1e0] $
     else if keyword_set(eor_test) then eor_data_range = [1e4, 5e7] else eor_data_range = [1e4, 5e6]
     kpower_2d_plots, savefile_eor, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = eor_data_range, pub = pub, plotfile = eor_plotfile, title = 'EoR ' + eor_tag, $
                      window_num = 10, norm_2d = norm_2d, grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, $
                      wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis
  endif

  if n_elements(test_power_shape) ne 0 then begin
     test_power_plotfile = plotfile_path + 'test_power/' + test_power_shape + 'power' + tile_tag + fadd + fadd_2dbin + '_kspace_power'
     if not keyword_set(undefined_data_range) then if keyword_set(norm_2d) then test_data_range = [1e-10, 1e0] $
     else begin
        case beam_exp of
           0: test_data_range = [1e8, 1e10]
           1: test_data_range = [1e9, 2e11]
           2: test_data_range = [1e11, 1e13]
        endcase
     endelse

     kpower_2d_plots, savefile_test_power, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = test_data_range, pub = pub, plotfile = test_power_plotfile, $
                      title = test_power_shape + ' Power test', window_num = 10, norm_2d = norm_2d, grey_scale = grey_scale, $
                      plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis
  endif
  
  if keyword_set(full_sky) then begin
     initial_uv_savefile = froot + fbase_arr + '_initial.idlsave'
     
     ;; plot input sky
     initial_image_plotfile = plotfile_path + sim_plot_folder + 'sim' + tile_tag + '_' + names[sim_num] + '_image_initial.eps'
     uvf_slice_plot, initial_uv_savefile, pub = pub, plotfile = initial_image_plotfile, window_num = 11, $
                     title = array + ' initial image', grey_scale = grey_scale, /mark_0, /image_space

     
     if keyword_set(plot_uvf) then begin
        ;; plot input uv
        initial_uv_plotfile = froot + fbase_arr + '_initial.eps'
        
        uvf_slice_plot, initial_uv_savefile, pub = pub, plotfile = initial_uv_plotfile, window_num = 10, plot_xrange = plot_urange, $
                        plot_yrange = plot_vrange, title = array + ' full sky initial uv plane (' + uvf_type + ')', $
                        grey_scale = grey_scale, baseline_axis = baseline_axis, /mark_0, type = uvf_type
     endif
  endif

  if keyword_set(plot_weights) then begin
     weight_plotfile = plotfile_path + 'weights/' + 'weights' + tile_tag + fadd + fadd_2dbin + '_kspace'
     if not keyword_set(undefined_data_range) then if keyword_set(norm_2d) then weights_range = [1e-10, 1e0] $
     else weights_range = [1e3, 1e8]

     kpower_2d_plots, savefiles_2d[0], /plot_weights, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      pub = pub, plotfile = weight_plotfile, window_num = 11, norm_2d = norm_2d, grey_scale = grey_scale, $
                      plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, baseline_axis = baseline_axis, delay_axis = delay_axis
  endif

  if n_sims eq 1 then begin
     if keyword_set(clean_ratio) then begin
        file_arr_2d = [savefiles_2d, reg_savefiles_2d]
        ratio = 1
     endif else if keyword_set(eor_ratio) then begin
        file_arr_2d = [savefile_eor, savefiles_2d] 
        ratio = 1
     endif else file_arr_2d = savefiles_2d
     
     if keyword_set(full_sky) then title_note = deg_offset_str[sim_num] + '!Uo!N' $
     else title_note =  ' [' + deg_offset_str[sim_num] + '!Uo!N, 0!Uo!N]'

     kpower_2d_plots, file_arr_2d, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = data_range, pub = pub, plotfile = plotfile_2d, title = array + ' ' + title_note, $
                      norm_2d = norm_2d, grey_scale = grey_scale,  plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, $
                      baseline_axis = baseline_axis, delay_axis = delay_axis, ratio = ratio
     
     file_arr = [savefiles_1d, eor_file_1d_input, eor_file_1d]
     names_arr = ['Simulation PS', 'input EoR', 'general EoR']
     colors_arr = [0, 210, 254]
     if keyword_set(eor_only) then begin
        file_arr = [file_arr, savefile_eor_1d]
        names_arr = [names_arr,'Simulated EoR']
        colors_arr = [colors_arr, 75]
     endif 
     if n_elements(test_power_shape) ne 0 then begin
       file_arr = [file_arr, savefile_test_power_1d]
       names_arr = [names_arr,test_power_shape + ' Power test']
       colors_arr = [colors_arr, 100]
    endif 


     if keyword_set(compare_1d) then begin
        file_arr = [savefiles_1d_compare, file_arr]
        if keyword_set(clean_ratio) then names_arr = [title_str + ' PS', names_arr] else names_arr = ['Naive PS', names_arr]
        colors_arr = [120, colors_arr]
     endif

     kpower_1d_plots, file_arr, pub = pub, plotfile = plotfile_1d, window_num = 2, names = names_arr, data_range = data_range_1d, $
                      colors = colors_arr, delta = delta

     if keyword_set(plot_slice) then begin
        if keyword_set(full_sky) then title_note = deg_offset_str[sim_num] + '!Uo!N' $
        else title_note =  ' [' + deg_offset_str[sim_num] + '!Uo!N, 0!Uo!N]'

         if keyword_set(slice_nobin) then begin
           yslice_savefile = froot + fbase_arr + fadd_3d + '_xz_plane.idlsave'
           xslice_savefile = froot + fbase_arr + fadd_3d + '_yz_plane.idlsave'

           kpower_slice_plot, yslice_savefile, data_range = slice_data_range, pub = pub, plotfile = yslice_plotfile, $
                              window_num = 3, title = array + ' XZ plane ' + title_note, $
                              grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, $
                              linear_axes = slice_linear_axes, baseline_axis = baseline_axis, delay_axis = delay_axis
           
           kpower_slice_plot, xslice_savefile, data_range = slice_data_range, pub = pub, plotfile = xslice_plotfile, $
                              window_num = 4, title = array + ' YZ plane ' + title_note, $
                              grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, $
                              linear_axes = slice_linear_axes, baseline_axis = baseline_axis, delay_axis = delay_axis
        endif else begin
           yslice_savefile = froot + fbase_arr + fadd_3d + '_xz_plane_binned.idlsave'

           kpower_2d_plots, yslice_savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                            data_range = slice_data_range, pub = pub, plotfile = yslice_plotfile, window_num = 3, $
                            title = array + ' XZ plane ' + title_note, norm_2d = norm_2d, $
                            grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp, $
                            baseline_axis = baseline_axis, delay_axis = delay_axis
           
           xslice_savefile = froot + fbase_arr + fadd_3d + '_yz_plane_binned.idlsave'
           
           kpower_2d_plots, xslice_savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                            data_range = slice_data_range, pub = pub, plotfile = xslice_plotfile, window_num = 4, $
                            title = array + ' YZ plane ' + title_note, norm_2d = norm_2d, $
                            grey_scale = grey_scale, baseline_axis = baseline_axis, delay_axis = delay_axis
        endelse

        zslice_savefile = froot + fbase_arr + fadd_3d + '_xy_plane.idlsave'

        kpower_slice_plot, zslice_savefile, data_range = slice_data_range, pub = pub, plotfile = zslice_plotfile, $
                           window_num = 5, title = array + ' XY plane ' + title_note, grey_scale = grey_scale, $
                           linear_axes = slice_linear_axes, baseline_axis = baseline_axis, delay_axis = delay_axis
     endif

     if keyword_set(plot_uvf) then begin
        if keyword_set(uvf_conv) then uvf_fadd = '_conv' else uvf_fadd = ''

        uvf_plot_fadd = uvf_fadd
        if keyword_set(grey_scale) then uvf_plot_fadd = uvf_plot_fadd + '_grey'
        uf_plotfile = plotfile_base + '_uf_plane' + uvf_plot_fadd + '.eps'
        vf_plotfile = plotfile_base + '_vf_plane' + uvf_plot_fadd + '.eps'
        uv_plotfile = plotfile_base + '_uv_plane' + uvf_plot_fadd + '.eps'

        uf_savefile = froot + fbase_arr + fadd_3d + '_uf_plane' + uvf_fadd + '.idlsave'
        vf_savefile = froot + fbase_arr + fadd_3d + '_vf_plane' + uvf_fadd + '.idlsave'
        uv_savefile = froot + fbase_arr + fadd_3d + '_uv_plane' + uvf_fadd + '.idlsave'

        if keyword_set(full_sky) then begin
           title_note = deg_offset_str[sim_num] + '!Uo!N'
        endif else begin
           title_note =  ' [' + deg_offset_str[sim_num] + '!Uo!N, 0!Uo!N]'
        endelse


        uvf_slice_plot, uf_savefile, pub = pub, plotfile = uf_plotfile, window_num = 6, plot_xrange = plot_urange, $
                        data_range = uvf_data_range, title = array + ' uf plane (' + uvf_type + ') ' + title_note, $
                        grey_scale = grey_scale, baseline_axis = baseline_axis, /mark_0, type = uvf_type
        
        uvf_slice_plot, vf_savefile, pub = pub, plotfile = vf_plotfile, window_num = 7, plot_xrange = plot_vrange, $
                        data_range = uvf_data_range, title = array + ' vf plane (' + uvf_type + ') ' + title_note, $
                        grey_scale = grey_scale, baseline_axis = baseline_axis, /mark_0, type = uvf_type
        
        uvf_slice_plot, uv_savefile, pub = pub, plotfile = uv_plotfile, window_num = 8, plot_xrange = plot_urange, $
                        plot_yrange = plot_vrange, title = array + ' uv plane (' + uvf_type + ') ' + title_note, $
                        data_range = uvf_data_range, grey_scale = grey_scale, baseline_axis = baseline_axis, /mark_0, type = uvf_type

     endif

  endif else begin

     ncol = ceil(sqrt(n_sims))
     nrow = ceil(n_sims / double(ncol))
     multi_params = {ncol:ncol, nrow:nrow, ordering:'col'}

     window_num = 1

     for i=0, n_sims-1 do begin
        if keyword_set(clean_ratio) then begin
           file_arr_2d = [savefiles_2d[i], reg_savefiles_2d[i]]
           ratio = 1
        endif else if keyword_set(eor_ratio) then begin
           file_arr_2d = [savefile_eor, savefiles_2d[i]]
           ratio = 1
        endif else file_arr_2d = savefiles_2d[i]
 
        if i gt 0 then pos_use = positions[*,i] else start_multi_params = multi_params

        kpower_2d_plots, file_arr_2d, multi_pos = pos_use, start_multi_params = start_multi_params, $
                         kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = data_range, pub = pub, $
                         title = array + ' [' + deg_offset_str[i] + '!Uo!N, 0!Uo!N]', norm_2d = norm_2d, grey_scale = grey_scale, $
                         plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp[i], baseline_axis = baseline_axis, $
                         delay_axis = delay_axis, ratio = ratio, plotfile = plotfile_2d, window_num = window_num
        if keyword_set(undefined_data_range) then temp = n_elements(temporary(data_range))

        if i eq 0 then begin
           positions = pos_use
           undefine, start_multi_params
        endif
     endfor

     if keyword_set(pub) then begin
        psoff
        wdelete, window_num
     endif

     xsize = !d.x_vsize
     ysize = !d.y_vsize
     window_num = 2
     if windowavailable(window_num) then begin 
        wset, window_num
        if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
     endif else make_win = 1
     if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
     erase

     if keyword_set(pub) then begin
        pson, file = plotfile_1d, /eps
     endif

     for i=0, n_sims-1 do begin
        if keyword_set(eor_only) then begin
           file_arr = [savefiles_1d[i], savefile_eor_1d, eor_file_1d]
           names_arr = ['Simulation PS','Simulated EoR', 'EoR signal']
        endif else begin
           file_arr = [savefiles_1d[i], eor_file_1d]
           names_arr = ['Simulation PS', 'EoR signal']
        endelse
        
        if keyword_set(compare_1d) then begin
           file_arr = [savefiles_1d_compare[i], file_arr]
           if keyword_set(clean_ratio) then names_arr = [title_str + ' PS', names_arr] else names_arr = ['Naive PS', names_arr]
        endif

        kpower_1d_plots, file_arr, multi_pos = positions[*,i], pub = pub, names = names_arr, data_range = data_range_1d 
     endfor

     if keyword_set(pub) then begin
        psoff
        wdelete, window_num
     endif


     if keyword_set(plot_slice) then begin

        window_num = 3

        if keyword_set(slice_nobin) then yslice_savefiles = froot + fbase_arr + fadd_3d + '_xz_plane.idlsave' $
        else yslice_savefiles = froot + fbase_arr + fadd_3d + '_xz_plane_binned.idlsave'

        undefine, pos_use
        for i=0, n_sims-1 do begin
           if i gt 0 then pos_use = positions[*,i] else start_multi_params = multi_params
        
           if keyword_set(slice_nobin) then begin

              kpower_slice_plot, yslice_savefiles[i], multi_pos = pos_use, start_multi_params = start_multi_params, $
                                 data_range = slice_data_range, pub = pub, plotfile = yslice_plotfile, window_num = window_num, $
                                 title = array + ' XZ plane [' + deg_offset_str[i] + '!Uo!N, 0!Uo!N]', $
                                 grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp[i], $
                                 linear_axes = slice_linear_axes, baseline_axis = baseline_axis, delay_axis = delay_axis
           endif else begin
              kpower_2d_plots, yslice_savefiles[i], multi_pos = pos_use, start_multi_params = start_multi_params, $
                               kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = slice_data_range, $
                               title = array + ' XZ plane [' + deg_offset_str[i] + '!Uo!N, 0!Uo!N]', pub = pub, $
                               plotfile = yslice_plotfile, window_num = window_num, $
                               norm_2d = norm_2d, grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, $
                               wedge_amp = wedge_amp[i], baseline_axis = baseline_axis, delay_axis = delay_axis
           endelse

           if keyword_set(undefined_data_range) then temp = n_elements(temporary(data_range))
           if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
           endif
        endfor

        if keyword_set(pub) then begin
           psoff
           wdelete, window_num
        endif
        
        
        if keyword_set(slice_nobin) then xslice_savefiles = froot + fbase_arr + fadd_3d + '_yz_plane.idlsave' $
        else xslice_savefiles = froot + fbase_arr + fadd_3d + '_yz_plane_binned.idlsave'

        window_num = 4
        
        undefine, pos_use
        for i=0, n_sims-1 do begin
           if i gt 0 then pos_use = positions[*,i] else start_multi_params = multi_params

           if keyword_set(slice_nobin) then begin
              kpower_slice_plot, xslice_savefiles[i], multi_pos = pos_use, start_multi_params = start_multi_params, $
                                 data_range = slice_data_range, pub = pub, plotfile = xslice_plotfile, window_num = window_num, $
                                 title = array + ' YZ plane [' + deg_offset_str[i] + '!Uo!N, 0!Uo!N]', $
                                 grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp[i], $
                                 linear_axes = slice_linear_axes, baseline_axis = baseline_axis, delay_axis = delay_axis
              
           endif else begin
              kpower_2d_plots, xslice_savefiles[i], multi_pos = pos_use, start_multi_params = start_multi_params, $
                               kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = slice_data_range, $
                               title = array + ' YZ plane [' + deg_offset_str[i] + '!Uo!N, 0!Uo!N]', pub = pub, $
                               plotfile = xslice_plotfile, window_num = window_num, $
                               norm_2d = norm_2d, grey_scale = grey_scale, plot_wedge_line = plot_wedge_line, $
                               wedge_amp = wedge_amp[i], baseline_axis = baseline_axis, delay_axis = delay_axis
           endelse

           if keyword_set(undefined_data_range) then temp = n_elements(temporary(data_range))
           if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
           endif
        endfor

        if keyword_set(pub) then begin
           psoff
           wdelete, window_num
        endif
        

        ;; zslice is never binned
        zslice_savefiles = froot + fbase_arr + fadd_3d + '_xy_plane.idlsave'

        window_num = 5
       
        undefine, pos_use
         for i=0, n_sims-1 do begin     
           if i gt 0 then pos_use = positions[*,i] else start_multi_params = multi_params

           kpower_slice_plot, zslice_savefiles[i], multi_pos = pos_use, start_multi_params = start_multi_params, $
                              data_range = slice_data_range, pub = pub, plotfile = zslice_plotfile, window_num = window_num, $
                              title = array + ' XY plane [' + deg_offset_str[i] + '!Uo!N, 0!Uo!N]', $
                              grey_scale = grey_scale, linear_axes = slice_linear_axes, baseline_axis = baseline_axis, $
                              delay_axis = delay_axis
           
           if keyword_set(undefined_data_range) then temp = n_elements(temporary(data_range))
           if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
           endif
        endfor

        if keyword_set(pub) then begin
           psoff
           wdelete, window_num
        endif

     endif


     if keyword_set(plot_uvf) then begin
        if keyword_set(grey_scale) then uvf_fadd = '_grey' else uvf_fadd = ''
        uf_plotfile = plotfile_base + '_uf_plane' + uvf_fadd + '.eps'
        vf_plotfile = plotfile_base + '_vf_plane' + uvf_fadd + '.eps'
        uv_plotfile = plotfile_base + '_uv_plane' + uvf_fadd + '.eps'

        uf_savefile = froot + fbase_arr + fadd_3d + '_uf_plane.idlsave'
        vf_savefile = froot + fbase_arr + fadd_3d + '_vf_plane.idlsave'
        uv_savefile = froot + fbase_arr + fadd_3d + '_uv_plane.idlsave'

        window_num = 6

        undefine, pos_use
        for i=0, n_sims-1 do begin
           if i gt 0 then pos_use = positions[*,i] else start_multi_params = multi_params

          uvf_slice_plot, uf_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, pub = pub, $
                          grey_scale = grey_scale, baseline_axis = baseline_axis, data_range = uvf_data_range, $
                          title = array + ' uf plane phase [' + deg_offset_str[sim_num[i]] + '!Uo!N, 0!Uo!N]', $
                          plotfile = uf_plotfile, window_num = window_num
                        
          if i eq 0 then begin
             positions = pos_use
             undefine, start_multi_params
          endif
       endfor

        if keyword_set(pub) then begin
           psoff
           wdelete, window_num
        endif

        
        window_num = 7

        undefine, pos_use
        for i=0, n_sims-1 do begin
           if i gt 0 then pos_use = positions[*,i] else start_multi_params = multi_params

           uvf_slice_plot, vf_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, pub = pub, $
                           grey_scale = grey_scale, baseline_axis = baseline_axis, data_range = uvf_data_range, $
                           title = array + ' vf plane phase [' + deg_offset_str[sim_num[i]] + '!Uo!N, 0!Uo!N]', $
                           plotfile = vf_plotfile, window_num = window_num
                        
          if i eq 0 then begin
             positions = pos_use
             undefine, start_multi_params
          endif
       endfor
        
        if keyword_set(pub) then begin
           psoff
           wdelete, window_num
        endif


        window_num = 8

        undefine, pos_use
         for i=0, n_sims-1 do begin
           if i gt 0 then pos_use = positions[*,i] else start_multi_params = multi_params
           
           uvf_slice_plot, uv_savefile[i], multi_pos = pos_use, start_multi_params = start_multi_params, pub = pub, $
                           grey_scale = grey_scale, baseline_axis = baseline_axis, data_range = uvf_data_range, $
                           title = array + ' uv plane phase [' + deg_offset_str[sim_num[i]] + '!Uo!N, 0!Uo!N]', $
                           plotfile = uv_plotfile, window_num = window_num
                        
           if i eq 0 then begin
              positions = pos_use
              undefine, start_multi_params
           endif
        endfor
        
        
        if keyword_set(pub) then begin
           psoff
           wdelete, window_num
        endif

     endif


  endelse



end
