pro sim_clean_fig, pub = pub, kpar_log = kpar_log, grey_scale = grey_scale

  if keyword_set(grey_scale) then   plotfile = base_path() + 'single_use/sim_clean_fig_grey.eps' $
  else plotfile = base_path() + 'single_use/sim_clean_fig.eps'

  froot = base_path() + 'fhd/simulations/'
  simfile = froot + 'sim_496t_fullsky_uvf_holo_linkpar_2dkpower.idlsave'
  clean_simfile = froot + 'sim_496t_fullsky_uvf_cleanplus_holo_linkpar_2dkpower.idlsave'

  info_file = froot + 'sim_496t_info.idlsave'
  restore, info_file   
  xy_length = 2048
  deg_offset = xy_length * degpix / sqrt(2d)
  rad_offset = deg_offset * !pi/180d

  redshift = 8
  cosmology_measures, redshift, wedge_factor = wedge_factor
  source_dists = rad_offset
  wedge_amp = wedge_factor * source_dists

  restore, simfile 
  kperp_centers = (kperp_edges[1:*] - kperp_edges[0:n_elements(kperp_edges)-2])/2d + kperp_edges[0:n_elements(kperp_edges)-2]
  
  kpar_plot_edges = kpar_edges
  if keyword_set(kpar_log) then kpar_plot_edges[0] = kpar_edges[1]^2d/kpar_edges[2]
  kpar_plot_centers = (kpar_plot_edges[1:*] - kpar_plot_edges[0:n_elements(kpar_edges)-2])/2d + $
                      kpar_plot_edges[0:n_elements(kpar_plot_edges)-2]
  

  ;;kperp_locs = [1.3e-2, 1e-1]
 ;; kperp_locs = [3e-2, 5e-2]
  kperp_locs = [4e-2]
  nlocs = n_elements(kperp_locs)
  bin_num = intarr(nlocs)  
  for i=0, nlocs-1 do bin_num[i] = max(where(kperp_edges lt kperp_locs[i]))
       
  sim_power_trace = power[bin_num, *]

  restore, clean_simfile
  clean_power_trace = power[bin_num, *]

  if keyword_set(grey_scale) then colors = ['dark grey', 'grey']  else colors = ['PBG5', 'red6'] 

  ;; divide plot area in quarters
  ncol = 2
  nrow = 2
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
  multi_aspect = .95
  xsize = round((max_ysize/nrow) * ncol/multi_aspect)
  ysize = max_ysize
 
  margin1 = [0.2, 0.15]
  margin2 = [0.1, 0.1]
  plot_pos = [margin1[0], margin1[1], (1-margin2[0]), (1-margin2[1])]

  multi_xlen = (positions[2,2]-positions[0,2])
  multi_ylen = (positions[3,2]-positions[1,2])
  
  new_pos = [multi_xlen * plot_pos[0] + positions[0,2], multi_ylen * plot_pos[1] + positions[1,2], $
             multi_xlen * plot_pos[2] + positions[0,2], multi_ylen * plot_pos[3] + positions[1,2]]
 

  window_num = 1
  if windowavailable(window_num) then begin 
     wset, window_num
     if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
  endif else make_win = 1
  if make_win eq 1 then cgdisplay, xsize, ysize, wid=window_num, xsize = xsize, ysize = ysize, color = 'white'
  cgerase, 'white'
  
  if keyword_set(pub) then begin
     charthick = 3
     thick = 3
     xthick = 3
     ythick = 3

     min_len = min([multi_xlen, multi_ylen])
     charsize = 5d * min_len

     font = 1

     pson, file = plotfile, /eps
  endif else begin
     thick = 1
     charsize=1
  endelse

  kperp_plot_range = [6e-3, 0.3]
  data_range = [1e11, 1e18]
  ratio_data_range = [1e-3, 1e0]
  ;;yrange = minmax([sim_power_trace, clean_power_trace])
  yrange = data_range


  kpower_2d_plots, simfile, multi_pos = positions[*,0], multi_aspect = multi_aspect, kperp_plot_range = kperp_plot_range, $
                   kpar_plot_range = kpar_plot_range, data_range = data_range, title='Initial Power', pub = pub, $
                   grey_scale = grey_scale, /baseline_axis, /plot_wedge_line, wedge_amp = wedge_amp
  for i=0, nlocs-1 do cgplot, /overplot, [kperp_centers[bin_num[i]], kperp_centers[bin_num[i]]], $
                              [kpar_edges[1]^2d/kpar_edges[2], max(kpar_edges)], color = colors[i], thick = thick+1, linestyle = 0

  kpower_2d_plots, clean_simfile, multi_pos = positions[*,1], multi_aspect = multi_aspect, kperp_plot_range = kperp_plot_range, $
                   kpar_plot_range = kpar_plot_range, data_range = data_range, title='Post-subtraction Power', pub = pub, $
                   grey_scale = grey_scale, /baseline_axis, /plot_wedge_line, wedge_amp = wedge_amp
  for i=0, nlocs-1 do cgplot, /overplot, [kperp_centers[bin_num[i]], kperp_centers[bin_num[i]]], $
                              [kpar_edges[1]^2d/kpar_edges[2], max(kpar_edges)], color = colors[i], thick = thick-1, linestyle = 0

  kpower_2d_plots, [clean_simfile, simfile], multi_pos = positions[*,3], multi_aspect = multi_aspect, $
                   kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, data_range = ratio_data_range, $
                   title='Ratio of Post-subtraction to Inital Power', pub = pub, grey_scale = grey_scale, /baseline_axis, /ratio, $
                   /plot_wedge_line, wedge_amp = wedge_amp


  cgplot, kpar_plot_centers, sim_power_trace[0,*], /noerase, /ylog, xlog=kpar_log, color=colors[0], axiscolor='black', $
          position = new_pos,  ytickformat = 'exponent', xtickformat='exponent', ytitle = textoidl(' (mK^2 Mpc^3)', font = font), $
          xtitle = textoidl('k_{||} (Mpc^{-1})', font = font), xstyle=1, yrange = yrange, ystyle=1, psym=10, thick = thick+1, $
          charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font

  cgplot, /overplot, kpar_plot_centers, clean_power_trace[0,*], color = colors[0], psym=10, thick = thick-1

  line_names = textoidl('k_{perp}=' + number_formatter(kperp_centers[bin_num[0]], format = '(e10.1)', /print_exp) + ' Mpc^{-1}', $
                        font = font)
  for i=1, nlocs-1 do begin
     cgplot, /overplot, kpar_plot_centers, sim_power_trace[i,*], color = colors[i], psym=10, thick = thick+1
     cgplot, /overplot, kpar_plot_centers, clean_power_trace[i,*], color = colors[i], psym=10, thick = thick-1
     line_names = [line_names, textoidl('k_{perp}=' + number_formatter(kperp_centers[bin_num[i]], format = '(e10.1)', /print_exp) + $
                                        ' Mpc^{-1}', font = font)]
  endfor

  al_legend, line_names, box=0, textcolors=colors, /right, charthick = charthick, charsize = charsize, font = font

  if keyword_set(pub) then begin
     psoff
     wdelete, window_num
  endif

end
