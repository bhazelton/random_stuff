


pro survey_plots, pub = pub, no_kzero = no_kzero, refresh_2d = refresh_2d, refresh_3dbin = refresh_3dbin, grey_scale = grey_scale, $
                  baseline_axis = baseline_axis, compare_1d = compare_1d, field = field, fill_holes = fill_holes, norm_2d = norm_2d, $
                  norm_factor = norm_factor, no_title = no_title, coarse = coarse, save_1d_text = save_1d_text, $
                  freq_band = freq_band, delta = delta, mid = mid, orig = orig, stokes = stokes, linear_kpar = linear_kpar, $
                  plot_wedge_line = plot_wedge_line

  ;; default to including baseline axis
  if n_elements(baseline_axis) eq 0 then baseline_axis = 1

  ;; default to using delta for 1d plots
  if n_elements(delta) eq 0 then delta=1

  if n_elements(norm_factor) ne 0 then norm_2d=1

  froot = base_path() + 'power_spectrum/healpix_maps/'

  if keyword_set(mid) or keyword_set(orig) then begin
     fname_arr = ['dirty_survey', 'clean_survey']
     if keyword_set(mid) then fname_arr = 'new_' + fname_arr
     if keyword_set(coarse) then fname_arr = fname_arr + '_coarse'
     fbase_arr = froot + fname_arr + '/' + fname_arr
     
     if n_elements(field) ne 0 then begin
        case field of
           '2122': begin
              field_tag = '_2122'
              field_title = ' RA 21-22'
           end
           '0405': begin
              field_tag = '_0405'
              field_title = ' RA 4-5'
           end
           else: message, 'unknown field'
        endcase
        fbase_arr = fbase_arr + field_tag
        fname_arr = fname_arr + field_tag
     endif else begin
        field_title = ''
        field_tag = ''
     endelse
     title = ['Dirty Survey', 'Clean Survey'] + field_title

     if keyword_set(freq_band) then begin
        fband_names = 'fband' + ['1', '2', '3', '4']
        fbase_arr = [fbase_arr[0] + '_' + fband_names, fbase_arr[1] + '_' + fband_names]
        
        fname_arr = fname_arr + '_fband'
     endif

  endif else begin
     fname_arr = 'lambda_15_' + ['unsubtracted', '5Jy_subtracted', '5Jy_subtracted_PCA_5comp', '5Jy_subtracted_PCA_10comp'] + '_0.6MHz'
     title = ['Unsubtracted', '5Jy Subtracted', '5Jy & 5 comp. PCA Subtracted', '5Jy & 10 comp. PCA Subtracted']
     fbase_arr = froot + fname_arr + '/' + fname_arr

     if n_elements(stokes) eq 0 then stokes = 0 else if stokes gt 2 or stokes lt 0 then message, $
        'Stokes must be an integer between 0 & 2 (I,Q,U)'
     if stokes gt 0 then begin
        stokes_name = '_stokes' + (['q','u'])[stokes-1]
        fbase_arr = fbase_arr[0:1] + stokes_name
        fname_arr = fname_arr[0:1] + stokes_name
        title = title[0:1] + ' Stokes ' + (['Q','U'])[stokes-1]
     endif else stokes_name = ''

  endelse
  
  nfiles = n_elements(fbase_arr)
  
  info_files = fbase_arr + '_info.txt'
  
  for i=0, nfiles-1 do begin
     info = read_info_file(info_files[i])
     test_file = file_test(info.metadata_file) *  (1 - file_test(info.metadata_file, /zero_length))  
     if test_file eq 0 then begin
        print, 'Metadata file does not exist. Running healpix_get_ft_frame now'
        healpix_get_ft_frame, info_files[i]
     endif
     undefine, info
  endfor
  
  fadd = ''
  if keyword_set(no_kzero) then fadd = fadd + '_nok0'
  fadd_1d = fadd
  if keyword_set(fill_holes) then fadd = fadd + '_nohole'
  if keyword_set(linear_kpar) then fadd = fadd + '_linkpar'

  savefiles_2d = fbase_arr + fadd + '_2dkpower.idlsave'
  savefiles_1d = fbase_arr + fadd_1d + '_1dkpower.idlsave'
  eor_file_1d = base_path() + 'power_spectrum/eor_data/eor_power_1d.idlsave'

  ftests_2d = file_test(savefiles_2d) *  (1 - file_test(savefiles_2d, /zero_length))
  ftests_1d = file_test(savefiles_1d) *  (1 - file_test(savefiles_1d, /zero_length))
  ftests = ftests_2d * ftests_1d

  if keyword_set(no_kzero) and keyword_set(compare_1d) then compare_1d = 0

  if keyword_set(compare_1d) then begin
     savefiles_1d_compare = fbase_arr + '_nok0_1dkpower.idlsave'
     ftests_1d_comp = file_test(savefiles_1d_compare) *  (1 - file_test(savefiles_1d_compare, /zero_length))
  endif

  for i = 0, nfiles-1 do begin
     if ftests[i] eq 0 or keyword_set(refresh_2d) or keyword_set(refresh_3dbin) then $
        healpix_to_3dps, info_files[i], no_kzero = no_kzero, refresh_2d = refresh_2d, refresh_3dbin = refresh_3dbin, /quiet, $
                         fill_holes = fill_holes

     if keyword_set(compare_1d) then $
        if ftests_1d_comp[i] eq 0 or keyword_set(refresh_2d) or keyword_set(refresh_3dbin) then $
           healpix_to_3dps, info_files[i], /no_kzero, refresh_2d = refresh_2d, refresh_3dbin = refresh_3dbin, /quiet, $
                            fill_holes = fill_holes
  endfor

  plot_fadd = ''
  if keyword_set(norm_2d) then plot_fadd = plot_fadd + '_norm'
  if keyword_set(grey_scale) then plot_fadd = plot_fadd + '_grey'

  plotfiles_2d =  base_path() + 'power_spectrum/plots/survey/' + fname_arr + fadd + $
                  '_kspace_power' + plot_fadd + '.eps'
  plotfiles_1d =  base_path() + 'power_spectrum/plots/survey/'

  if keyword_set(mid) or keyword_set(orig) then begin
     if keyword_set(mid) then plotfiles_1d = plotfiles_1d + 'new_'
     plotfiles_1d = plotfiles_1d +'survey'
     if keyword_set(coarse) then plotfiles_1d = plotfiles_1d + '_coarse'
     if keyword_set(freq_band) then plotfiles_1d = plotfiles_1d + '_fband'
     plotfiles_1d = plotfiles_1d + field_tag
  endif else begin
     plotfiles_1d = plotfiles_1d + 'lambda_15' + stokes_name
     if keyword_set(freq_band) then plotfiles_1d = plotfiles_1d + '_fband'
  endelse
  plotfiles_1d = plotfiles_1d + fadd_1d 
  if keyword_set(delta) then plotfiles_1d = plotfiles_1d + '_1d_delta' else plotfiles_1d = plotfiles_1d + '_1d_power'
  
  if keyword_set(norm_2d) then data_range = [1e-6, 1e0] $
  else if keyword_set(mid) or (keyword_set(orig) and n_elements(field) ne 0) then data_range = [1e17, 5e21] $
  else if keyword_set(orig) and n_elements(field) eq 0 then data_range = [1e18, 1e23] $
  else data_range = [1e19, 1e25]

  ;; Gianni removed baselines < 40 lambda (15 lambda in new surveys) Don't plot kperp (much) below this.
  restore, savefiles_2d[0]
  undefine, power
  undefine, weights
  if keyword_set(orig) then min_baseline = 38d else min_baseline = 15d
  min_kperp = min_baseline / kperp_lambda_conv

  mean_freq = 188.8d
  max_baseline = 342.497 / (300d/mean_freq)
  max_kperp = max_baseline / kperp_lambda_conv
  ;; print, min_kperp
  ;;min_kperp=2e-3

  kperp_plot_range = [min_kperp, max_kperp]
  
  if keyword_set(plot_wedge_line) then begin
     z0_freq = 1420.40 ;; MHz
     redshift = z0_freq/mean_freq - 1
     cosmology_measures, redshift, wedge_factor = wedge_factor
     wedge_amp = wedge_factor * (20d*!dpi/180d)
  endif


  if not keyword_set(freq_band) then begin
     
     for i=0, nfiles-1 do kpower_2d_plots, savefiles_2d[i], kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                                           data_range = data_range, pub = pub, plotfile = plotfiles_2d[i], no_title = no_title, $
                                           window_num = 1+i, grey_scale = grey_scale, baseline_axis = baseline_axis, title = title[i],$
                                           norm_2d = norm_2d, norm_factor = norm_factor, plot_wedge_line = plot_wedge_line, $
                                           wedge_amp = wedge_amp
  
  endif else begin
     n_plots = n_elements(fband_names)

     ncol = ceil(sqrt(n_plots))
     nrow = ceil(n_plots / double(ncol))
     positions = fltarr(4, ncol*nrow)

     row_val = reverse(reform(rebin(indgen(nrow), nrow, ncol), ncol*nrow))
     col_val = reform(rebin(reform(indgen(ncol), 1, ncol), nrow, ncol), ncol*nrow)

     if keyword_set(pub) then begin
        xmargin = 0.025
        ymargin = 0.025
     endif else begin
        xmargin = 0.0125
        ymargin = 0.0125
     endelse

     positions[0,*] = col_val/double(ncol)+xmargin
     positions[1,*] = row_val/double(nrow)+ymargin
     positions[2,*] = (col_val+1)/double(ncol)-xmargin
     positions[3,*] = (row_val+1)/double(nrow)-ymargin

     max_ysize = 800
     if n_plots gt 9 then begin
        multi_aspect =0.9 
        xsize = round((max_ysize/nrow) * ncol/multi_aspect)
        ysize = max_ysize
     endif else begin
        if keyword_set(baseline_axis) then multi_aspect = 0.75 else multi_aspect =0.5
        xsize = round((max_ysize/nrow) * ncol/multi_aspect)
        ysize = max_ysize
     endelse

     for j=0, nfiles/n_plots do begin

        window_num = j+1
        if windowavailable(window_num) then begin 
           wset, window_num
           if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
        endif else make_win = 1
        if make_win eq 1 then window, window_num, xsize = xsize, ysize = ysize
        erase

        if keyword_set(pub) then pson, file = plotfiles_2d[j], /eps
  
        for i=0, n_plots-1 do kpower_2d_plots, savefiles_2d[i+j*n_plots], multi_pos = positions[*,i], multi_aspect = multi_aspect, $
                                               kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                                               data_range = data_range, pub = pub, title = title[j] + ' ' + fband_names[i], $
                                               norm_2d = norm_2d,norm_factor = norm_factor, grey_scale = grey_scale, $
                                               baseline_axis = baseline_axis, plot_wedge_line = plot_wedge_line, wedge_amp = wedge_amp
                         
        if keyword_set(pub) then begin
           psoff
           wdelete, window_num
        endif
     endfor

  endelse

  file_arr = [savefiles_1d, eor_file_1d]
  if keyword_set(freq_band) then begin
     names_arr = title[0] + fband_names
     for i=1, n_elements(title) do names_arr = [names_arr, title[i] + fband_names]
     names_arr = [names_arr, 'EOR signal']

     if n_elements(title) eq 2 then colors_arr = [indgen(4)*15+30, indgen(4)*15+200, 0] else stop
  endif else begin
     names_arr = [title, 'EOR signal']
     colors_arr = [indgen(nfiles) * (254-45)/(nfiles-1) + 45, 0]

     if keyword_set(compare_1d) then begin
        file_arr = [file_arr, savefiles_1d_compare]
        names_arr = [names_arr, title + ' w/o k0']
        colors_arr = [colors_arr, indgen(nfiles) * (210-75)/(nfiles-1)+75]
     endif
  endelse

  if keyword_set(no_kzero) then begin
     if keyword_set(delta) then data_range_1d = [1e-1, 1e8] else data_range_1d = [1e3, 1e22]
     k_range = [1e-2, 1e0]
  endif else begin
     if keyword_set(delta) then data_range_1d = [1e-2, 1e10] else data_range_1d = [1e3, 1e25]
     k_range = [1e-3, 1e0]
  endelse
  kpower_1d_plots, file_arr, window_num=5, data_range = data_range_1d, pub = pub, k_range = k_range, $
                   names = names_arr, colors = colors_arr, plotfile = plotfiles_1d, save_text = save_1d_text, delta = delta
  ;; kpower_1d_plots, file_arr, window_num=3, pub = pub, k_range = k_range, names = names_arr, $
  ;;                  colors = colors_arr, plotfile = plotfiles_1d, save_text = save_1d_text
 
end
