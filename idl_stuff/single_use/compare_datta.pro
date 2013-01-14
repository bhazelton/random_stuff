

pro compare_datta

  savefile = base_path() + 'single_use/datta_data/' + 'MWA_GSM_0.1err'
  savefile = savefile + '_2dkpower_edgegrid_match.idlsave'
  restore, savefile

  ;; power, kperp_edges, kpar_edges, bins_per_decade
  if bins_per_decade ne 10d then message, 'Not sampled with 10 bins/decade'
  my_power = power
  log_binsize = 1d/bins_per_decade
  my_kperp_edges = kperp_edges
  my_kpar_edges = kpar_edges
  ;; my_kperp = 10^(alog10(kperp_edges[1:*]) - log_binsize/2d)
  ;; my_kpar = 10^(alog10(kpar_edges[1:*]) - log_binsize/2d)
  my_kperp = kperp_edges[1:*]
  my_kpar = kpar_edges[1:*]

  restore, base_path() + 'single_use/datta_data/' + 'figure_10a.idlsave'
  ;; res_power, res_kperp, res_kpar
  datta_bins_per_decade = 10d

  res_kpar[1] = 10^(alog10(res_kpar[1]) + log_binsize)

  ;; the datta array has some bins missing -- copy from lower bins
  datta_kperp = 10^(dindgen(round((max(alog10(res_kperp)) - min(alog10(res_kperp)))*datta_bins_per_decade) + 1) / $
                    datta_bins_per_decade + min(alog10(res_kperp)))
  datta_kpar = 10^(dindgen(round((max(alog10(res_kpar)) - min(alog10(res_kpar)))*datta_bins_per_decade) + 1) / $
                    datta_bins_per_decade + min(alog10(res_kpar)))

  test_eq_val = 10^(-4d)
  kperp_inds_to_use = dblarr(n_elements(datta_kperp))
  kpar_inds_to_use = dblarr(n_elements(datta_kpar))
  for i=0, max([n_elements(datta_kperp), n_elements(datta_kpar)]) -1 do begin
     if i lt n_elements(datta_kperp) then begin
        wh_kperp = where(abs(res_kperp - datta_kperp[i]) lt test_eq_val, count)
        if count eq 1 then kperp_inds_to_use[i] = wh_kperp[0] $
        else kperp_inds_to_use[i] = kperp_inds_to_use[i-1]
     endif
     
     if i lt n_elements(datta_kpar) then begin
        wh_kpar= where(abs(res_kpar - datta_kpar[i]) lt test_eq_val, count)
        if count eq 1 then kpar_inds_to_use[i] = wh_kpar[0] $
        else kpar_inds_to_use[i] = kpar_inds_to_use[i-1]
     endif
  endfor
  
  datta_power = dblarr(n_elements(datta_kperp), n_elements(datta_kpar))
  for i=0, n_elements(datta_kperp)-1 do begin
     for j=0, n_elements(datta_kpar) -1 do datta_power[i,j] = res_power[kperp_inds_to_use[i], kpar_inds_to_use[j]]
  endfor

  if keyword_set(edge_on_grid) then begin
     datta_kperp = 10^(alog10(datta_kperp) - log_binsize/2d)
     datta_kpar = 10^(alog10(datta_kpar) - log_binsize/2d)
  endif

  ;;res_power is log -- convert to linear
  datta_power = 10^datta_power

  kperp_range = [max([min(my_kperp), min(datta_kperp)]), min([max(my_kperp), max(datta_kperp)])]
  kpar_range = [max([min(my_kpar), min(datta_kpar)]), min([max(my_kpar), max(datta_kpar)])]

  ;; limit my data to calculated ranges
  wh_kperp_under = where(my_kperp - kperp_range[0] lt -1*test_eq_val, count_kperp_under)
  wh_kperp_over = where(my_kperp - kperp_range[1] gt test_eq_val, count_kperp_over)
  if count_kperp_under ne 0 or count_kperp_over ne 0 then begin
     wh_inrange = where(my_kperp - kperp_range[0] gt -1*test_eq_val and my_kperp - kperp_range[1] lt test_eq_val)
     my_power = my_power[wh_inrange, *]
     my_kperp = my_kperp[wh_inrange]
     my_kperp_edges = my_kperp_edges[[wh_inrange, max(wh_inrange) + 1]]
  endif

  wh_kpar_under = where(my_kpar - kpar_range[0] lt -1*test_eq_val, count_kpar_under)
  wh_kpar_over = where(my_kpar - kpar_range[1] gt test_eq_val, count_kpar_over)
  if count_kpar_under ne 0 or count_kpar_over ne 0 then begin
     wh_inrange = where(my_kpar - kpar_range[0] gt -1*test_eq_val and my_kpar - kpar_range[1] lt test_eq_val)
     my_power = my_power[*, wh_inrange]
     my_kpar = my_kpar[wh_inrange]
     my_kpar_edges = my_kpar_edges[[wh_inrange, max(wh_inrange) + 1]]
  endif

  ;; limit datta data to calculated ranges
  wh_kperp_under = where(datta_kperp - kperp_range[0] lt -1*test_eq_val, count_kperp_under)
  wh_kperp_over = where(datta_kperp - kperp_range[1] gt test_eq_val, count_kperp_over)
  if count_kperp_under ne 0 or count_kperp_over ne 0 then begin
     wh_inrange = where(datta_kperp - kperp_range[0] gt -1*test_eq_val and datta_kperp - kperp_range[1] lt test_eq_val)
    datta_power = datta_power[wh_inrange, *]
    datta_kperp = datta_kperp[wh_inrange]
  endif

  wh_kpar_under = where(datta_kpar - kpar_range[0] lt -1*test_eq_val, count_kpar_under)
  wh_kpar_over = where(datta_kpar - kpar_range[1] gt test_eq_val, count_kpar_over)
  if count_kpar_under ne 0 or count_kpar_over ne 0 then begin
     wh_inrange = where(datta_kpar - kpar_range[0] gt -1*test_eq_val and datta_kpar - kpar_range[1] lt test_eq_val)
     datta_power = datta_power[*, wh_inrange]
     datta_kpar = datta_kpar[wh_inrange]
  endif

  n_kperp = n_elements(my_kperp)
  n_kpar = n_elements(my_kpar)

  power_ratio = my_power / datta_power
  power_diff = my_power - datta_power
 
  color_range = [45, 254]
  n_colors = color_range[1] - color_range[0]
  
  
  power_plot = dblarr(4, n_kperp*10, n_kpar*10)
  power_plot[0, *, *] = congrid(datta_power, n_kperp*10, n_kpar*10)
  power_plot[1, *, *] = congrid(my_power, n_kperp*10, n_kpar*10)
  power_plot[2, *, *] = congrid(power_ratio, n_kperp*10, n_kpar*10)
  power_plot[3, *, *] = congrid(power_diff, n_kperp*10, n_kpar*10)

  power_log = dblarr(4, n_kperp*10, n_kpar*10)
  power_log_norm = dblarr(4, n_kperp*10, n_kpar*10)
  data_range = dblarr(4, 2)

  ;; data_range[0,*] = alog10(minmax(power_plot[0:1,*, *]))
  ;; data_range[1,*] = data_range[0,*]
  ;; data_range[2,*] = alog10(minmax(power_plot[2,*, *]))

  data_range[0,*] = [10^(-2d), 10^8d]
  data_range[1,*] = data_range[0,*]
  data_range[2,*] = minmax(power_plot[2,*, *])
  data_range[3,*] = minmax(power_plot[3,*, *])

  log_data_range = alog10(data_range)

  for i=0, 3 do begin
     ;; data_range[i,*] = alog10(minmax(power_plot[i,*, *]))
 
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

  tvlct, r, g, b, /get
   

  cb_size = 0.025               ;; in units of plot area (incl. margins)
  big_margin = [0.15, 0.2]      ;; in units of plot area (incl. margins)
  small_margin = [0.05, 0.1]    ;; in units of plot area (incl. margins)
  
  plot_pos = [big_margin[0], big_margin[1], (1-2*small_margin[0]-cb_size), (1-small_margin[1])]
  cb_pos = [(1-small_margin[0]-cb_size), big_margin[1], (1-small_margin[0]), (1-small_margin[1])]
  
  plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])
  
  kpar_length_log = max(alog10(my_kpar_edges)) - min(alog10(my_kpar_edges))
  kperp_length_log = max(alog10(my_kperp_edges)) - min(alog10(my_kperp_edges))
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
     font = 1
     
     window, 4, xsize = 700 * x_factor, ysize = 700 * y_factor
     
     
     plotfile = 'datta_kspace_comp'
     if keyword_set(bw) then plotfile = plotfile + '_bw'
     plotfile = base_path() + 'single_use/datta_plots/' + plotfile
     
     pson, file = plotfile + '.eps', /eps 
  endif

  loadct, 39
  
  plot_kperp = my_kperp_edges
  plot_kpar = my_kpar_edges

  if keyword_set(pub) then nplots = 1 $
  else nplots = 4
  
  titles = ['Datta Result ', 'My Result ', 'My Result/Datta Result ', 'My Result-Datta Result '] + 'P!Dk!N (mK!E2!N Mpc!E3!N)'
  for i=0, nplots-1 do begin
     if keyword_set(pub) then ind = 2 $
     else begin 
        ind = i
        if windowavailable(ind+1) then wset, ind+1 else window, ind+1, xsize = 700 * x_factor, ysize = 700 * y_factor
     endelse

     image = power_log_norm[ind,*,*]
  
     plot, plot_kperp, plot_kpar, /xlog, /ylog, /nodata, xstyle=4, ystyle=4,   title = titles[ind], $
           position = plot_pos, thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, $
           font = font
     cgimage, image, /nointerp,/overplot
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
     
     temp = [ceil(log_data_range[ind, 0]), floor(log_data_range[ind, 1])]
     tick_vals = 10^(dindgen(temp[1]-temp[0] + 1) + temp[0])
     if (alog10(tick_vals[0]) -log_data_range[ind, 0]) gt 10^(-3d) then begin
        cb_ticknames = [' ', string(round(alog10(tick_vals)))]
        cb_ticks = [color_range[0], (alog10(tick_vals) - log_data_range[ind, 0]) * n_colors / $
                    (log_data_range[ind, 1] - log_data_range[ind, 0]) + color_range[0]] - color_range[0]
     endif else begin
        cb_ticknames = string(round(alog10(tick_vals)))
        cb_ticks = ((alog10(tick_vals) - log_data_range[ind, 0]) * n_colors / $
                    (log_data_range[ind, 1] - log_data_range[ind, 0]) + color_range[0]) - color_range[0]
     endelse
     
     cgcolorbar, color = 0, /vertical, position = cb_pos, bottom = color_range[0], ncolors = n_colors, yminor = 0, $
                   ticknames = cb_ticknames, ytickv = cb_ticks, divisions = n_elements(cb_ticks) -1, charsize = charsize, font = font
  endfor

  if keyword_set(pub) then begin
     psoff
     wdelete, 4
  endif
  
  tvlct, r, g, b




end
