pro gianni_paper_plots, pub = pub


  froot = base_path() + 'power_spectrum/healpix_maps/'
  fname_arr = 'lambda_15_' + ['unsubtracted', '5Jy_subtracted', '5Jy_subtracted_PCA_5comp'] + '_0.6MHz'
  title = ['Unsubtracted', '5Jy Subtracted', '5Jy & 5 comp. PCA Subtracted']
  fbase_arr = froot + fname_arr + '/' + fname_arr
  stokes_name = ''

  nfiles = n_elements(fbase_arr)
  savefiles_2d = fbase_arr + '_nohole_linkpar_2dkpower.idlsave'
  ftests_2d = file_test(savefiles_2d) *  (1 - file_test(savefiles_2d, /zero_length))

  for i = 0, nfiles-1 do $
     if ftests_2d[i] eq 0 or keyword_set(refresh_2d) or keyword_set(refresh_3dbin) then $
        healpix_to_3dps, info_files[i], refresh_2d = refresh_2d, refresh_3dbin = refresh_3dbin, /quiet, $
                         /fill_holes

  data_range = [1e19, 1e25]
  ratio_data_range = [1e-3, 1e0]

  restore, savefiles_2d[0]
  undefine, power
  undefine, weights
  min_baseline = 15d
  min_kperp = min_baseline / kperp_lambda_conv

  mean_freq = 188.8d
  max_baseline = 342.497 / (300d/mean_freq)
  max_kperp = max_baseline / kperp_lambda_conv

  kperp_plot_range = [min_kperp, max_kperp]

  plotfiles = base_path() + 'single_use/gianni_paper_' + ['initial_fig.eps', 'ps_sub_fig.eps', 'pca5_sub_fig.eps']

  ;; make initial power figure
  kpower_2d_plots, savefiles_2d[0], kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = data_range, $
                   /no_title, pub = pub, /baseline_axis, charsize = charsize, plotfile=plotfiles[0], $
                   cb_margin = [0.14, 0.04], margin = [0.15, 0.1, 0.02, 0.11] 

  ;; make point-source subtracted figure
  ncol = 2
  nrow = 1
  positions = fltarr(4, ncol*nrow)
  
  row_val = reverse(reform(rebin(reform(indgen(nrow), 1, nrow), ncol, nrow), ncol*nrow))
  col_val = reform(rebin(indgen(ncol), ncol, nrow), ncol*nrow)
  
  xmargin = 0.
  ymargin = 0.
  
  positions[0,*] = col_val/double(ncol)+xmargin
  positions[1,*] = row_val/double(nrow)+ymargin
  positions[2,*] = (col_val+1)/double(ncol)-xmargin
  positions[3,*] = (row_val+1)/double(nrow)-ymargin
  
  max_ysize = 1100
  multi_aspect = 1.5
  xsize = round((max_ysize/nrow) * ncol/multi_aspect)
  ysize = max_ysize
  
  window_num = 2
  if windowavailable(window_num) then begin 
     wset, window_num
     if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
  endif else make_win = 1
  if make_win eq 1 then cgdisplay, xsize, ysize, wid=window_num, xsize = xsize, ysize = ysize, color = 'white'
  cgerase, 'white'
  
  if keyword_set(pub) then begin
     pson, file = plotfiles[1], /eps
  endif

  kpower_2d_plots, savefiles_2d[1],  multi_pos = positions[*,0], multi_aspect = multi_aspect, kperp_plot_range = kperp_plot_range, $
                   kpar_plot_range = kpar_plot_range, data_range = data_range, /no_title, pub = pub, /baseline_axis

  kpower_2d_plots, [savefiles_2d[1], savefiles_2d[0]], /ratio, multi_pos = positions[*,1], multi_aspect = multi_aspect, $
                   kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = ratio_data_range, $
                   /no_title, pub = pub, /baseline_axis

  if keyword_set(pub) then begin
     psoff
     wdelete, window_num
  endif

  ;; make PCA subtracted figure
  window_num = 3
  if windowavailable(window_num) then begin 
     wset, window_num
     if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
  endif else make_win = 1
  if make_win eq 1 then cgdisplay, xsize, ysize, wid=window_num, xsize = xsize, ysize = ysize, color = 'white'
  cgerase, 'white'
  
  if keyword_set(pub) then begin
     pson, file = plotfiles[2], /eps
  endif

  kpower_2d_plots, savefiles_2d[2],  multi_pos = positions[*,0], multi_aspect = multi_aspect, kperp_plot_range = kperp_plot_range, $
                   kpar_plot_range = kpar_plot_range, data_range = data_range, /no_title, pub = pub, /baseline_axis

  kpower_2d_plots, [savefiles_2d[2], savefiles_2d[0]], /ratio, multi_pos = positions[*,1], multi_aspect = multi_aspect, $
                   kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = ratio_data_range, $
                   /no_title, pub = pub, /baseline_axis

  if keyword_set(pub) then begin
     psoff
     wdelete, window_num
  endif

end
