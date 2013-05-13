pro gianni_paper_plots, pub = pub


  froot = base_path('data') + 'healpix_maps/'
  fname_arr = 'lambda_15_' + ['unsubtracted', '5Jy_subtracted', '5Jy_subtracted_PCA_5comp'] + '_0.6MHz'
  title = ['Unsubtracted', '5Jy Subtracted', '5Jy & 5 comp. PCA Subtracted']
  fbase_arr = froot + fname_arr + '/' + fname_arr
  stokes_name = ''

  nfiles = n_elements(fbase_arr)
  savefiles_2d = fbase_arr + '_nohole_logkperp_2dkpower.idlsave'
  ftests_2d = file_test(savefiles_2d) *  (1 - file_test(savefiles_2d, /zero_length))

  info_files = fbase_arr + '_info.txt'

  ;;refresh_3dbin=1

  for i = 0, nfiles-1 do $
     if ftests_2d[i] eq 0 or keyword_set(refresh_2d) or keyword_set(refresh_3dbin) then $
        healpix_to_3dps, info_files[i], refresh_2d = refresh_2d, refresh_3dbin = refresh_3dbin, /quiet, $
                         /fill_holes, /log_kperp, kperp_bin=.1, /log_k1d, k1d_bin = .1

  data_range = [1e-4, 1e0]
  ratio_data_range = [1e-3, 1e0]
  ;;diff_data_range = [1-4, 1e0]

  restore, savefiles_2d[0]
  undefine, power
  undefine, weights
  min_baseline = 15d
  min_kperp = min_baseline / kperp_lambda_conv

  mean_freq = 188.8d
  max_baseline = 342.497 / (300d/mean_freq)
  max_kperp = max_baseline / kperp_lambda_conv

  kperp_plot_range = [min_kperp, max_kperp]

  plotfiles = base_path('plots') + 'single_use/gianni_paper_' + $
              ['initial_fig.eps', 'ps_sub_fig.eps', 'ps_sub_ratio_fig.eps', 'ps_sub_diff_fig.eps', 'pca5_sub_fig.eps', $
               'pca5_sub_ratio_fig.eps', 'pca5_sub_diff_fig.eps']

  ;; make initial power figure
  kpower_2d_plots, savefiles_2d[0], kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range=data_range, $
                   /no_title, pub = pub, /baseline_axis, charsize = charsize, plotfile=plotfiles[0], $
                   cb_margin = [0.14, 0.04], margin = [0.15, 0.1, 0.02, 0.11], /norm_2d, norm_factor = norm_factor

  ;; make point-source subtracted figures  
  kpower_2d_plots, savefiles_2d[1], kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range=data_range, $
                   /no_title, pub = pub, /baseline_axis, charsize = charsize, plotfile=plotfiles[1], window_num = 2, $
                   cb_margin = [0.14, 0.04], margin = [0.15, 0.1, 0.02, 0.11], /norm_2d, norm_factor = norm_factor

  kpower_2d_plots, [savefiles_2d[1], savefiles_2d[0]], /ratio, kperp_plot_range = kperp_plot_range, $
                   kpar_plot_range = kpar_plot_range, data_range = ratio_data_range, /no_title, pub = pub, /baseline_axis, $
                   charsize = charsize, plotfile=plotfiles[2], window_num = 3, $
                   cb_margin = [0.14, 0.04], margin = [0.15, 0.1, 0.02, 0.11]
  
  kpower_2d_plots, [savefiles_2d[0], savefiles_2d[1]], /diff, kperp_plot_range = kperp_plot_range, $
                   kpar_plot_range = kpar_plot_range, data_range = diff_data_range, /no_title, pub = pub, /baseline_axis, $
                   charsize = charsize, plotfile=plotfiles[3], window_num = 4, /norm_2d, norm_factor = norm_factor, $
                   cb_margin = [0.14, 0.04], margin = [0.15, 0.1, 0.02, 0.11]

  ;; make PCA subtracted figure

  kpower_2d_plots, savefiles_2d[2], kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range=data_range, $
                   /no_title, pub = pub, /baseline_axis, charsize = charsize, plotfile=plotfiles[4], window_num = 5, $
                   cb_margin = [0.14, 0.04], margin = [0.15, 0.1, 0.02, 0.11], /norm_2d, norm_factor = norm_factor

  kpower_2d_plots, [savefiles_2d[2], savefiles_2d[0]], /ratio, kperp_plot_range = kperp_plot_range, $
                   kpar_plot_range = kpar_plot_range, data_range = ratio_data_range, /no_title, pub = pub, /baseline_axis, $
                   charsize = charsize, plotfile=plotfiles[5], window_num = 6, $
                   cb_margin = [0.14, 0.04], margin = [0.15, 0.1, 0.02, 0.11]
  
  kpower_2d_plots, [savefiles_2d[0], savefiles_2d[2]], /diff, kperp_plot_range = kperp_plot_range, $
                   kpar_plot_range = kpar_plot_range, data_range = diff_data_range, /no_title, pub = pub, /baseline_axis, $
                   charsize = charsize, plotfile=plotfiles[6], window_num = 7, /norm_2d, norm_factor = norm_factor, $
                   cb_margin = [0.14, 0.04], margin = [0.15, 0.1, 0.02, 0.11]

 
end
