
pro ps_crosspower, info_files, comp_type = comp_type, baseline_axis = baseline_axis, linear_kpar = linear_kpar, pub = pub

  ;; default to including baseline axis
  if n_elements(baseline_axis) eq 0 then baseline_axis = 1

  types_enum = ['interleave', 'stokes']
  if n_elements(comp_type) eq 0 then comp_type = 'stokes'
  wh_type = where(types_enum eq comp_type, test_type)
  if test_type eq 0 then message, 'unrecognized comparison type'

  if n_elements(info_files) eq 0 then begin
     froot = base_path() + 'power_spectrum/healpix_maps/'
     ;;info_files = froot + ['even', 'odd'] + '_info.txt'
     fname = 'lambda_15_unsubtracted_0.6MHz'
     info_files = froot + fname + '/' + fname + ['', '_stokesq'] +  '_info.txt'
  endif else if n_elements(info_files) ne 2 then message, 'Exactly two files must be listed in info_files'

  plotfile_path = base_path() + 'power_spectrum/plots/'

  info1 = read_info_file(info_files[0])
  info2 = read_info_file(info_files[1])

  power_file1 = info1.power_file
  power_file2 = info2.power_file

  ftest1 = file_test(power_file1) *  (1 - file_test(power_file1, /zero_length))  
  if ftest1 eq 0 then message, 'power file does not exist for info file: ' + info_files[0]

  ftest2 = file_test(power_file2) *  (1 - file_test(power_file2, /zero_length))  
  if ftest2 eq 0 then message, 'power file does not exist for info file: ' + info_files[1]

  
  restore, power_file1
  ft1 = image_ft
  kx1 = kx_mpc
  ky1 = ky_mpc
  kz1 = kz_mpc

  restore, power_file2
  ft2 = image_ft
  kx2 = kx_mpc
  ky2 = ky_mpc
  kz2 = kz_mpc

  if not array_equal(kx1, kx1) or not array_equal(ky1, ky2) or not array_equal(kz1, kz2) then $
     message, 'The k values for each power spectrum are not equal.'


  cross_power = ft1 * conj(ft2)
  power1 = real_part(ft1 * conj(ft1))
  power2 = real_part(ft2 * conj(ft2))

  froot1 = strmid(info_files[0], 0, strpos(info_files[0], '/', /reverse_search)) + '/'
  froot2 = strmid(info_files[0], 0, strpos(info_files[0], '/', /reverse_search)) + '/'
  if not strcmp(froot1, froot2, /fold_case) then print, 'info files have different paths, saving to ' + froot1
  froot = froot1
  
  min_len = min([strlen(info1.id), strlen(info2.id)])
  for i=1, min_len do begin
     if strcmp(info1.id, info2.id, i) eq 0 then begin
        match_len=i-1
        break
     endif
  endfor
  if n_elements(match_len) eq 0 then match_len = min_len
  base_name = strmid(info1.id, 0, match_len)
  sub_names = [strmid(info1.id, match_len), strmid(info2.id, match_len)]
  for i=0, 1 do begin
     if sub_names[i] eq '' then sub_names[i] = '_reg'
     if strmid(sub_names[i], 0, 1) ne '_' then sub_names[i] = '_' + subnames[i]
  endfor
  savefile_base = froot + base_name + sub_names[0] + sub_names[1]
  plotfile_base = plotfile_path + base_name + sub_names[0] + sub_names[1]

  if keyword_set(linear_kpar) then begin
     savefile_base = savefile_base + '_linkpar'
     plotfile_base = plotfile_base + '_linkpar'
  endif     

  case comp_type of
     'interleave': begin
        real_cross = real_part(cross_power)
        im_cross = imaginary(cross_power)
        undefine, cross_power
        
        cross_ratio = abs(real_cross) / abs(im_cross)
        wh_im_zero = where(im_cross eq 0, count_im_zero, complement = wh_comp, ncomplement = count_comp)
        if count_im_zero gt 0 and count_im_zero lt n_elements(im_cross) then cross_ratio[wh_im_zero] = max(cross_ratio[wh_comp]) $
        else if count_im_zero eq n_elements(im_cross) then cross_ratio = cross_ratio * 0d
        
        ;; power_diff = power1 - power2
        power_diff = power1 + power2 - 2d*real_cross
        ;;power_diff = power1/real_cross
        
        bins_per_decade = 10
        real_cross_rebin = kspace_rebinning_2d(real_cross, kx1, ky1, kz1, kperp_edges_mpc, kpar_edges_mpc, bins = bins_per_decade, $
                                               linear_kpar = linear_kpar)
        
        power = real_cross_rebin
        kperp_edges = kperp_edges_mpc
        kpar_edges = kpar_edges_mpc
        
        savefile = savefile_base + '_2dcrosspower.idlsave'
        save, file = savefile, power, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv
        undefine, power
        
        plotfile = plotfile_base + '_2dcrosspower.eps'
        
        kpower_2d_plots, savefile, plotfile = plotfile, title = 'Real Crosspower', baseline_axis = baseline_axis, $
                         pub = pub
        
        
        im_cross_rebin = kspace_rebinning_2d(im_cross, kx1, ky1, kz1, kperp_edges_mpc, kpar_edges_mpc, bins = bins_per_decade, $
                                             linear_kpar = linear_kpar)
        
        power = im_cross_rebin
        kperp_edges = kperp_edges_mpc
        kpar_edges = kpar_edges_mpc
        savefile = savefile_base + '_2dimagcrosspower.idlsave'
        save, file = savefile, power, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv
        undefine, power

        
        plotfile = plotfile_base + '_2dimagcrosspower.eps'
        
        kpower_2d_plots, savefile, plotfile = plotfile, title = '|Imaginary Crosspower|', window_num = 2, $
                         baseline_axis = baseline_axis, pub = pub, color_profile = 'abs'
        
        
        cross_ratio_rebin = kspace_rebinning_2d(cross_ratio, kx1, ky1, kz1, kperp_edges_mpc, kpar_edges_mpc, bins = bins_per_decade, $
                                                linear_kpar = linear_kpar)
        
        if count_im_zero eq 1 and wh_im_zero[0] eq where(kx1 eq 0) + n_elements(kx1)*where(ky1 eq 0) then $
           cross_ratio_rebin[0] = max(cross_ratio_rebin[1:*])*10d
        
        power = cross_ratio_rebin
        kperp_edges = kperp_edges_mpc
        kpar_edges = kpar_edges_mpc
        savefile = savefile_base + '_2dcrosspower_ratio.idlsave'
        save, file = savefile, power, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv
        undefine, power
        
        plotfile = plotfile_base + '_2dcrosspower_ratio.eps'
        
        kpower_2d_plots, savefile, plotfile = plotfile, title = '|Real Crosspower|/|Imaginary Crosspower|', window_num = 3, $
                         baseline_axis = baseline_axis, pub = pub, /no_units
        
        
        power_diff_rebin = kspace_rebinning_2d(power_diff, kx1, ky1, kz1, kperp_edges_mpc, kpar_edges_mpc, bins = bins_per_decade, $
                                               linear_kpar = linear_kpar)
        
        power = power_diff_rebin
        kperp_edges = kperp_edges_mpc
        kpar_edges = kpar_edges_mpc
        savefile = savefile_base + '_2dpowerdiff.idlsave'
        save, file = savefile, power, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv
        undefine, power
        
        plotfile = plotfile_base + '_2dpowerdiff.eps'
        
        kpower_2d_plots, savefile, plotfile = plotfile, title = 'Sum of Autopowers - 2*Real Crosspower', $
                         window_num = 4, baseline_axis = baseline_axis, pub = pub, color_profile = 'abs'
     end
     
     'stokes': begin
        ;;undefine, power1
        ;;undefine, power2

        abs_cross = abs(cross_power)

        bins_per_decade = 10
        power_abs_rebin = kspace_rebinning_2d(abs_cross/(power1), kx1, ky1, kz1, kperp_edges_mpc, kpar_edges_mpc, $
                                              bins = bins_per_decade, linear_kpar = linear_kpar)

        power = power_abs_rebin
        kperp_edges = kperp_edges_mpc
        kpar_edges = kpar_edges_mpc
        savefile = savefile_base + '_2dabscrosspower.idlsave'
        save, file = savefile, power, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv
        ;;undefine, power
        ;;undefine, abs_cross

        min_baseline = 15d
        min_kperp = min_baseline / kperp_lambda_conv

        mean_freq = 188.8d
        max_baseline = 342.497 / (300d/mean_freq)
        max_kperp = max_baseline / kperp_lambda_conv

        kperp_plot_range = [min_kperp, max_kperp]
  
        
        plotfile = plotfile_base + '_2dabscrosspower.eps'
        
        kpower_2d_plots, savefile, plotfile = plotfile, title = '<|crosspower|/Stokes I power>', baseline_axis = baseline_axis, $
                         pub = pub, kperp_plot_range = kperp_plot_range, /no_units


        crosspower_rebin = kspace_rebinning_2d(cross_power/(power1), kx1, ky1, kz1, kperp_edges_mpc, kpar_edges_mpc, $
                                               bins = bins_per_decade, linear_kpar = linear_kpar)

        power = abs(crosspower_rebin)
        kperp_edges = kperp_edges_mpc
        kpar_edges = kpar_edges_mpc
        savefile = savefile_base + '_2dcrosspower.idlsave'
        save, file = savefile, power, kperp_edges, kpar_edges, bins_per_decade, kperp_lambda_conv
        undefine, power
        
        plotfile = plotfile_base + '_2dcrosspower.eps'
        
        kpower_2d_plots, savefile, plotfile = plotfile, title = '|<crosspower/Stokes I power>|', window_num=2, $
                         baseline_axis = baseline_axis, pub = pub, kperp_plot_range = kperp_plot_range, /no_units

;;stop
     end
  endcase
  
end
