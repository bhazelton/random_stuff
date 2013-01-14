

pro datta_sensitivity_plots, pub = pub, match = match, contours = contours, shading = shading


  restore, base_path() + 'single_use/datta_data/' + 'MWA_GSM_0.1err_2dkpower_edgegrid_match.idlsave'
  
  res_power = power
  res_kperp_edges = kperp_edges
  res_kpar_edges = kpar_edges
  res_log_binsize = 1d/bins_per_decade


  restore, base_path() + 'single_use/datta_data/' + 'eor_signal.idlsave'
  eor_log_binsize = 1d/10d ;; 10 bins per decade
  eor_kperp_edges = 10^([alog10(eor_kperp[0]) - eor_log_binsize, alog10(eor_kperp)])
  eor_kpar_edges = 10^([alog10(eor_kpar[0]) - eor_log_binsize, alog10(eor_kpar)])
  eor_power = 10^eor_power ;; saved as log power

  kperp_range = [max([min(eor_kperp_edges), min(res_kperp_edges)]), min([max(eor_kperp_edges), max(res_kperp_edges)])]
  kpar_range = [max([min(eor_kpar_edges), min(res_kpar_edges)]), min([max(eor_kpar_edges), max(res_kpar_edges)])]

  ;; limit residuals to calculated ranges
  wh_kperp_under = where(res_kperp_edges lt kperp_range[0], count_kperp_under)
  wh_kperp_over = where(res_kperp_edges gt kperp_range[1], count_kperp_over)
  if count_kperp_under ne 0 or count_kperp_over ne 0 then begin
     wh_inrange = where(res_kperp_edges ge kperp_range[0] and res_kperp_edges[1:*] le kperp_range[1])
     res_power = res_power[wh_inrange, *]
     res_kperp_edges = res_kperp_edges[wh_inrange, max(wh_inrange)+1]
  endif

  wh_kpar_under = where(res_kpar_edges lt kpar_range[0], count_kpar_under)
  wh_kpar_over = where(res_kpar_edges gt kpar_range[1], count_kpar_over)
  if count_kpar_under ne 0 or count_kpar_over ne 0 then begin
     wh_inrange = where(res_kpar_edges ge kpar_range[0] and res_kpar_edges[1:*] le kpar_range[1])
     res_power = res_power[*, wh_inrange]
     res_kpar_edges = res_kpar_edges[wh_inrange, max(wh_inrange)+1]
  endif

  dims = size(res_power, /dimension)
  n_kperp = dims[0]
  n_kpar = dims[1]

  ;; rebin eor signal to have same bins as residual
  ;; find equivalent bin locations, if they exist
  test_eq_val = 10^(-4d)
  eor_kperp_bin = intarr(n_kperp)
  eor_kpar_bin = intarr(n_kpar)
  for i = 0, max([n_kperp, n_kpar])-1 do begin
     if i lt n_kperp then begin
        wh = where(abs(eor_kperp_edges - res_kperp_edges[i]) lt test_eq_val and $
                   abs(eor_kperp_edges[1:*] - res_kperp_edges[i+1]) lt test_eq_val, count)
        if count eq 1 then eor_kperp_bin[i] = wh[0] else eor_kperp_bin[i] = -1
     endif
     
     if i lt n_kpar then begin
        wh = where(abs(eor_kpar_edges - res_kpar_edges[i]) lt test_eq_val and $
                   abs(eor_kpar_edges[1:*] - res_kpar_edges[i+1]) lt test_eq_val, count)
        if count eq 1 then eor_kpar_bin[i] = wh[0] else eor_kpar_bin[i] = -1
     endif
  endfor
        
  if min(eor_kperp_bin) ge 0 and min(eor_kpar_bin) ge 0 then begin
     ;; all the bins line up with residual bins
     eor_rebin = eor_power[eor_kperp_bin, *]
     eor_rebin = eor_rebin[*, eor_kpar_bin]
  endif else begin
     ;; the bins don't all line up, have to do some interpolation
     stop
     eor_rebin = dblarr(n_kperp, n_kpar)
     for i = 0, n_kperp do begin
        for j=0, n_kpar do begin
           if eor_kperp_edges[i] gt -1 and eor_kpar_edges[j] gt -1 then $
              eor_rebin[i,j] = eor_power[eor_kperp_edges[i], eor_kpar_edges[j]] $
           else begin
              stop
              wh_kperp_inrange = where((eor_kperp_edges - res_kperp_edges[i]) ge -1*test_eq_val $
                                       and (eor_kperp_edges - res_kperp_edges[j]) le test_eq_val, count_kperp)
              wh_kpar_inrange = where((eor_kpar_edges-res_kpar_edges[i]) ge -1*test_eq_val $
                                      and (eor_kpar_edges-res_kpar_edges[j]) le test_eq_val, count_kpar)
   
              kperp_interp_loc = (alog10(res_kperp_edges[i]) + res_log_binsize/2d - min(alog10(eor_kperp))) / eor_log_binsize
              kpar_interp_loc = (alog10(res_kpar_edges[i]) + res_log_binsize/2d - min(alog10(eor_kpar))) / eor_log_binsize
          
              if eor_kperp_edges[i] gt -1 then eor_rebin[i,j] = interpolate(eor_power[eor_kperp_edges[i], *], kpar_interp_loc) $
              else if eor_kpar_edges[j] gt -1 then eor_rebin[i,j] = interpolate(eor_power[*, eor_kpar_edges[i]], kperp_interp_loc) $
              else eor_rebin[i,j] = interpolate(eor_power, kperp_interp_loc, kpar_interp_loc, /grid)
                 
           endelse
        endfor
     endfor

  endelse

  sig2res = eor_rebin / res_power
  sig2res_log = alog10(sig2res)

  tvlct, r, g, b, /get

  ;; expand image array to prevent interpolation in postscript, normalize
  color_range = [45, 254]
  n_colors = color_range[1] - color_range[0]

  
  nplots = 3
  power_plot = dblarr(nplots, n_kperp*10, n_kpar*10)
  power_plot[0, *, *] = congrid(res_power, n_kperp*10, n_kpar*10)
  power_plot[1, *, *] = congrid(eor_rebin, n_kperp*10, n_kpar*10)
  power_plot[2, *, *] = congrid(sig2res, n_kperp*10, n_kpar*10)
 
  power_log = dblarr(nplots, n_kperp*10, n_kpar*10)
  power_log_norm = dblarr(nplots, n_kperp*10, n_kpar*10)
  data_range = dblarr(nplots, 2)

  data_range[0,*] = [10^(-2d), 10^8d]
  data_range[1,*] = data_range[0,*]
  data_range[2,*] = minmax(power_plot[2,*, *])

  log_data_range =  alog10(data_range)

  for i=0, nplots-1 do begin

     power_log[i,*, *] = alog10(power_plot[i,*, *])
     temp = power_log[i,*, *]
     wh_under = where(temp lt log_data_range[i, 0], count)
     if count ne 0 then temp[wh_under] = log_data_range[i, 0]
     wh_over = where(temp gt log_data_range[i, 1], count)
     if count ne 0 then temp[wh_over] = log_data_range[i, 1]
     power_log[i,*, *] = temp
     power_log_norm[i,*, *] = (power_log[i,*, *]  - log_data_range[i, 0]) * n_colors / $
                              (log_data_range[i, 1] - log_data_range[i, 0]) + color_range[0]
  endfor

  cb_size = 0.025               ;; in units of plot area (incl. margins)
  big_margin = [0.15, 0.2]      ;; in units of plot area (incl. margins)
  small_margin = [0.05, 0.1]    ;; in units of plot area (incl. margins)
  
  plot_pos = [big_margin[0], big_margin[1], (1-small_margin[0]*3-cb_size), (1-small_margin[1])]
  cb_pos = [(1-small_margin[0]-cb_size), big_margin[1], (1-small_margin[0]), (1-small_margin[1])]
  cb_title_pos = [1-small_margin[0]*2-cb_size, 0.5]

  plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])
  
  kpar_length_log = max(alog10(res_kpar_edges)) - min(alog10(res_kpar_edges))
  kperp_length_log = max(alog10(res_kperp_edges)) - min(alog10(res_kperp_edges))
  data_aspect = (kpar_length_log / kperp_length_log)
  aspect_ratio =  data_aspect /plot_aspect
  
  if aspect_ratio gt 1 then begin
     y_factor = 1.
     x_factor = 1/aspect_ratio
  endif else begin
     y_factor = aspect_ratio
     x_factor = 1.
  endelse
  
  if keyword_set(pub) then begin
     charthick = 3
     thick = 3
     xthick = 3
     ythick = 3
     charsize = 2
     font=1
     
     window, 4, xsize = 700 * x_factor, ysize = 700 * y_factor
     
     plotfiles = ['residual', 'eor', 'sensitivity']
     if keyword_set(contours) then plotfiles = plotfiles + ['_sens_contours', '_sens_contours', '_contours']
     if keyword_set(shading) then plotfiles = plotfiles + '_shaded'
     plotfiles = base_path() + 'single_use/datta_plots/datta_GSM_0.1err_' + plotfiles + '.eps'
     
  endif else begin
     charthick = 1
     thick = 1
     xthick = 1
     ythick = 1
     charsize = 1.25
     font = -1
  endelse
  
  loadct, 39
  
  plot_kperp = res_kperp_edges
  plot_kpar = res_kpar_edges

  ;; get contours and mask
  contour_sig2res = [sig2res[0,*], sig2res, sig2res[n_kperp-1, *]]
  contour_sig2res = transpose([transpose(contour_sig2res[*,0]), transpose(contour_sig2res), $
                               transpose(contour_sig2res[*, n_kpar-1])])
  contour_kperp = [10^(alog10(res_kperp_edges[0]) - res_log_binsize/2d),10^(alog10(res_kperp_edges) + res_log_binsize/2d)]
  contour_kpar = [10^(alog10(res_kpar_edges[0]) - res_log_binsize/2d), 10^(alog10(res_kpar_edges) + res_log_binsize/2d)]
  
  wh_ge1 = where(contour_sig2res ge 1d, count_ge1)
  if count_ge1 eq 0 then message, 'No area with sig/noise gt 1' $
  else begin
     kperp_inds = wh_ge1 mod (n_kperp + 2)
     kpar_inds = wh_ge1 / (n_kperp + 2)
     ;; limit wh_ge1 to not be on edge
     wh_notedge = where(kperp_inds gt 0 and kpar_inds gt 0 and kperp_inds lt n_kperp+1 and kpar_inds lt n_kpar+1, count_notedge)
     if count_notedge eq 0 then message, 'Serious error in expanding sig2res array for contour' $
     else begin
        new_wh_ge1 = wh_ge1[wh_notedge]
        kperps = contour_kperp[kperp_inds[wh_notedge]]
        wh_start = where(kperps eq min(kperps))
        region_start = new_wh_ge1[where(kperps eq min(kperps))]
        
        region = region_grow(contour_sig2res, region_start, /all_neighbors, threshold = [1d, max(contour_sig2res)])
        mask = intarr(size(contour_sig2res, /dimensions)) + 1
        mask[region] = 0
        
     endelse
  endelse
  
  if keyword_set(shading) then titles = ['Residual with shaded exclusion zone', 'EOR signal with shaded exclusion zone', $
                  'EOR signal/residual with shaded exclusion zone'] $
  else if keyword_set(contours) then titles = ['Residual with signal/residual contours', 'EOR signal with signal/residual contours', $
                                               'EOR signal/residual'] $
  else titles = ['Residual', 'EOR signal', 'EOR signal/residual']
  
  for i=0, nplots-1 do begin
     if keyword_set(pub) then pson, file = plotfiles[i], /eps $
     else if windowavailable(i+1) then wset, i+1 else window, i+1, xsize = 700 * x_factor, ysize = 700 * y_factor

     if keyword_set(shading) then $
        contour, mask, contour_kperp, contour_kpar, levels = [1], c_labels = [0], path_info = path_info, path_xy = path_xy, $
                 /path_data_coords, /path_double

     plot, plot_kperp, plot_kpar, /xlog, /ylog, /nodata, xstyle=5, ystyle=5, xrange = minmax(plot_kperp), $
           yrange = minmax(plot_kpar), title = titles[i], position = plot_pos, thick = thick, charthick = charthick, $
           xthick = xthick, ythick = ythick, charsize = charsize, font = font
     tvimage, power_log_norm[i,*, *], /nointerp,/overplot
     if keyword_set(contours) then $
        contour, contour_sig2res, contour_kperp, contour_kpar, levels = [.1, 1, 10, 100], c_annotation = ['0.1', '1', '10', '100'], $
                 c_colors = [0], c_labels = intarr(4)+1, /overplot, c_thick = thick, c_charthick = charthick, c_charsize = charsize, $
                 font = font
     if keyword_set(shading) then begin
           x_vals = reform(path_xy[0,*])
           y_vals = reform(path_xy[1,*])

           polyfill, x_vals, y_vals, /line_fill, orientation = 45, /data, noclip=0, spacing = 0.1
           polyfill, x_vals, y_vals, /line_fill, orientation = -45, /data, noclip=0, spacing = 0.1
     endif

     axis, xaxis=0, xtick_get = xticks, xtitle = 'k!Dperp!N (Mpc!E-1!N)', xrange = minmax(plot_kperp), $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
           xtickformat = 'exponent', xstyle = 1
     axis, yaxis=0, ytick_get = yticks, ytitle = 'k!Dpar!N (Mpc!E-1!N)', yrange = minmax(plot_kpar), $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, $
        ytickformat = 'exponent', ystyle = 1
     axis, xaxis=1, xrange = minmax(plot_kperp), xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, xstyle = 1
     axis, yaxis=1, yrange = minmax(plot_kpar), ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font, ystyle = 1

     temp = [ceil(log_data_range[i,0]), floor(log_data_range[i,1])]
     tick_vals = 10^(dindgen(temp[1]-temp[0] + 1) + temp[0])
     if (alog10(tick_vals[0]) - log_data_range[i,0]) gt 10^(-3d) then begin
        cb_ticknames = [' ', string(round(alog10(tick_vals)))]
        cb_ticks = [color_range[0], (alog10(tick_vals) - log_data_range[i,0]) * n_colors / $
                    (log_data_range[i,1] - log_data_range[i,0]) + color_range[0]] - color_range[0]
     endif else begin
        cb_ticknames = string(round(alog10(tick_vals)))
        cb_ticks = ((alog10(tick_vals) - log_data_range[i,0]) * n_colors / $
                    (log_data_range[i,1] - log_data_range[i,0]) + color_range[0]) - color_range[0]
     endelse

     cgcolorbar, color = 0, /vertical, position = cb_pos, bottom = color_range[0], ncolors = n_colors, yminor = 0, $
                 ticknames = cb_ticknames, ytickv = cb_ticks, divisions = n_elements(cb_ticks) -1, charthick = charthick, $
                 xthick = xthick, ythick = ythick, charsize = charsize, font = font

     xyouts, cb_title_pos[0], cb_title_pos[1], 'log!D10!N P!Dk!N (mK!E2!N Mpc!E3!N)', orientation = 90, /normal, alignment=0.5, $
             charthick = charthick, charsize = charsize, font = font

     if keyword_set(pub) then psoff
  endfor 

  if keyword_set(pub) then wdelete, 4
  
  tvlct, r, g, b


  ;; Now do 1D plotting
  ;; get 1d spectrum with kpar=0 mode
  if keyword_set(match) then restore, base_path() + 'single_use/datta_data/' + 'MWA_GSM_0.1err_1dkpower_edgegrid_match.idlsave' $
  else restore, base_path() + 'single_use/datta_data/' + 'MWA_GSM_0.1err_1dkpower.idlsave'
  res_1d = power
  res_k_edges = k_edges
  res_bpd = bins_per_decade
  res_nk = n_elements(res_1d)
  res_k_centers = 10^(alog10(res_k_edges[1:*]) - (1d/res_bpd)/2.)
  
  ;; without kpar=0 mode
  if keyword_set(match) then restore, base_path() + 'single_use/datta_data/' + $
                                      'MWA_GSM_0.1err_1dkpower_nok0_edgegrid_match.idlsave' $
  else restore, base_path() + 'single_use/datta_data/' + 'MWA_GSM_0.1err_1dkpower_nok0.idlsave'
  res_1d_nok0 = power
  res_k_edges_nok0 = k_edges
  res_bpd_nok0 = bins_per_decade
  res_nk_nok0 = n_elements(res_1d_nok0)
  res_k_centers_nok0 = 10^(alog10(res_k_edges_nok0[1:*]) - (1d/res_bpd_nok0)/2.)
  
  ;;for comparison, rebin from 2d spectrum
  binning_kperp = 10^(alog10(res_kperp_edges[1:*]) - res_log_binsize/2d)
  binning_kperp[0] = 0
  binning_kpar = 10^(alog10(res_kpar_edges[1:*]) - res_log_binsize/2d)
  binning_kpar[0] = 0
  res_bpd_from2d = 5
  res_1d_from2d = kspace_rebinning_1d(res_power, binning_kperp, binning_kpar, 0, res_k_edges_from2d, $
                                      bins_per_decade = res_bpd_from2d, /quiet)
  res_nk_from2d = n_elements(res_1d_from2d)
  res_k_centers_from2d = 10^(alog10(res_k_edges_from2d[1:*]) - (1d/res_bpd_from2d)/2.)

  ;;for comparison, rebin from 2d spectrum with mask
  binning_mask = 1- mask[1:n_kperp, 1:n_kpar]
  res_1d_from2d_masked = kspace_rebinning_1d(res_power, binning_kperp, binning_kpar, bins_per_decade = res_bpd_from2d, $
                                             mask = binning_mask, /pixelwise_mask, /quiet)

  ;;for comparison, rebin from 3d spectrum with mask
  savefile = base_path() + 'single_use/datta_data/' + 'MWA_GSM_0.1err_1dkpower_masked.idlsave'
  test_save = file_test(savefile) *  (1 - file_test(savefile, /zero_length))
  refresh = 0
  if test_save eq 1 then begin
     restore, savefile
     if total(abs(size(binning_mask, /dimensions) - size(mask, /dimensions))) ne 0 then refresh = 1
  endif else refresh = 1
  if refresh eq 1 then begin
     datta_3d, mask = binning_mask, k1_mask = res_kperp_edges, k2_mask = res_kpar_edges
     restore, savefile
  endif
  res_1d_masked = power
  res_k_edges_masked = k_edges
  res_bpd_masked = bins_per_decade
  res_nk_masked = n_elements(res_1d_masked)
  res_k_centers_masked = 10^(alog10(res_k_edges_masked[1:*]) - (1d/res_bpd_masked)/2.)

  restore, base_path() + 'single_use/datta_data/eor_power_1d.idlsave'
  eor_1d = power
  eor_k_centers = k_centers
  
  res_theory = sqrt(res_1d * res_k_centers^3d / (2d*!pi^2d))
  res_theory_nok0 = sqrt(res_1d_nok0 * res_k_centers_nok0^3d / (2d*!pi^2d))
  res_theory_masked = sqrt(res_1d_masked * res_k_centers_masked^3d / (2d*!pi^2d))
  res_theory_from2d = sqrt(res_1d_from2d * res_k_centers_from2d^3d / (2d*!pi^2d))
  res_theory_from2d_masked = sqrt(res_1d_from2d_masked * res_k_centers_from2d^3d / (2d*!pi^2d))
  eor_theory = sqrt(eor_1d * eor_k_centers^3d / (2d*!pi^2d))
  
  if windowavailable(nplots+1) then wset, nplots+1 else window, nplots+1
  
  xrange = minmax(res_k_centers)
  yrange = 10^[-5d, 3d]
  plot, res_k_centers, res_theory, /ylog, /xlog, xrange = xrange, yrange = yrange, xstyle=1, ystyle=1, xtitle = 'k (Mpc!E-1!N)', $
        ytitle = '(K!E3!N P!Dk!N(2pi!E2!N))!E1/2!N (mK)', psym=-3
  oplot, res_k_centers, res_theory, psym=-3, color = 254

  oplot, eor_k_centers, eor_theory, psym=-3, color = 0
  oplot, res_k_centers_nok0, res_theory_nok0, psym=-3, color = 45
  oplot, res_k_centers_masked, res_theory_masked, psym=-3, color = 85, linestyle = 2
  ;; oplot, res_k_centers_from2d, res_theory_from2d, psym=-3, color = 254
  ;; oplot, res_k_centers_from2d, res_theory_from2d_masked, psym=-3, color = 85
  
  
  
end
