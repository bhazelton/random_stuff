pro test_fhd_sim_suite, uvf_input = uvf_input, flat_sigma = flat_sigma, recalculate_ps = recalculate_ps, $
    no_ps_plots = no_ps_plots
    
  sample_factors = [.0001,.001,.01]
  
  if n_elements(flat_sigma) eq 0 then flat_sigma = 1
  if keyword_set(flat_sigma) then sim_in = 'flat' else sim_in = 'eor'
  folder_names = 'fhd_bjh_arrsim_' + sim_in + '_density' + number_formatter(sample_factors)
  
  if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then pub = 1 else pub = 0
  if pub eq 1 then begin
    if not (keyword_set(png) or keyword_set(eps) or keyword_set(pdf)) then begin
      basename = cgRootName(plotfile, directory=directory, extension=extension)
      
      case extension of
        'eps': eps=1
        'png': png=1
        'pdf': pdf=1
        '': png = 1
        else: begin
          print, 'Unrecognized extension, using png'
          png = 1
        end
      endcase
      
    endif
    
    if keyword_set(png) and keyword_set(eps) and keyword_set(pdf) then begin
      print, 'only one of eps, pdf and png can be set, using png'
      eps = 0
    endif
    
    if keyword_set(png) then begin
      plot_exten = '.png'
      delete_ps = 1
    endif else if keyword_set(pdf) then begin
      plot_exten = '.pdf'
      delete_ps = 1
    endif else if keyword_set(eps) then begin
      plot_exten = '.eps'
      delete_ps = 0
    endif
  endif
  
  if n_elements(window_num) eq 0 then window_num = 2
  
  if n_elements(sample_inds) gt 0 then sample_factors = sample_factors[sample_inds]
  
  sim_powers = fltarr(nbeams, n_elements(sample_factors))
  sim_ave_powers = fltarr(nbeams, n_elements(sample_factors))
  sim_wt_ave_powers = fltarr(nbeams, n_elements(sample_factors))
  sim_ave_weights = fltarr(nbeams, n_elements(sample_factors))
  
  if keyword_set(pub) then no_plots=1
  if keyword_set(no_ps_plots) then begin
    no_plots = 1
    plot_stdset = 0
  endif
  
  plotfile = base_path('plots') + 'power_spectrum/fhd_sim/power_ratio_vs_weight_'
  
  plotfile_freq = base_path('plots') + 'power_spectrum/fhd_sim/power_per_freq_vs_weight_'
  
  t0 = systime(1)
  
  nbeams=1
  for i=0, nbeams-1 do begin
    for j=0, n_elements(sample_factors)-1 do begin
    
      print, 'calculating ps for ' + folder_name
      hellebore_wrapper, folder_names[i],/sim, /no_spec_window, png = png, $
        refresh_ps=recalculate_ps, refresh_beam=recalculate_ps, refresh_info=recalculate_ps, $
        cube_power_info = cube_power_info, plot_stdset = plot_stdset, uvf_input = uvf_input
        
      sim_ave_powers[i,j] = cube_power_info.ave_power[1]
      sim_wt_ave_powers[i,j] = cube_power_info.wt_ave_power[1]
      sim_ave_weights[i,j] = cube_power_info.ave_weights[1]
      
      if i eq 0 and j eq 0 then begin
        dims = size(cube_power_info.ave_power_freq, /dimension)
        sim_ave_power_freq = fltarr(nbeams, n_elements(sample_factors), dims[1])
        sim_wt_ave_power_freq = fltarr(nbeams, n_elements(sample_factors), dims[1])
        sim_ave_weights_freq = fltarr(nbeams, n_elements(sample_factors), dims[1])
      endif
      sim_ave_power_freq[i,j,*] = cube_power_info.ave_power_freq[1,*]
      sim_wt_ave_power_freq[i,j,*] = cube_power_info.wt_ave_power_freq[1,*]
      sim_ave_weights_freq[i,j,*] = cube_power_info.ave_weights_freq[1,*]
      
      if i+j eq 0 then flat_power = cube_power_info.flat_power else if flat_power ne cube_power_info.flat_power then print, 'flat powers do not agree'
      
    endfor
  endfor
  t1=systime(1)
  
  run_time = t1-t0
  if run_time lt 60 then time_str = number_formatter(run_time) + ' s' $
  else if run_time lt 3600 then time_str = number_formatter(run_time/60.) + ' m' else time_str = number_formatter(run_time/3600.) + ' h'
  print, 'run time: ' + time_str
  
  if keyword_set(recalculate_all) then print, sim_powers
  
  yrange=[0,1]
  xrange=[-1,9]
  
  colors1 = ['cyan', 'light salmon', 'light sea green', 'peru']
  colors2 = ['blue', 'red', 'sea green', 'chocolate']
  
  tvlct, r, g, b, /get
  
  if keyword_set(pub) then begin
    plotfile = plotfile + plot_exten
    
    charthick = 3
    thick = 3
    xthick = 3
    ythick = 3
    charsize = 2
    font = 1
    legend_charsize = 2
    
    if n_elements(multi_pos) eq 0 then begin
      cgps_open, plotfile, /font, encapsulated=eps, landscape=1, pagetype='letter'
    endif
    
  endif else if windowavailable(window_num) then wset, window_num else window, window_num
  
  
  cgplot, sim_ave_weights[0,*], sim_ave_powers[0,*]*0+0.52, color='black', yrange = yrange, xtitle='ave weight', ytitle = 'power ratio', xrange=xrange, $
    thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font
  for i=0, nbeams-1 do cgplot, sim_ave_weights[i,*], sim_ave_powers[i,*]/flat_power, color=colors1[i], psym=-4, /over, thick = thick
  for i=0, nbeams-1 do cgplot, sim_ave_weights[i,*], sim_wt_ave_powers[i,*]/flat_power, color=colors2[i], /over, psym=-4, thick = thick
  
  beam_str = strarr(nbeams)
  for i=0, nbeams-1 do begin
    case sim_beam[i] of
      1: beam_str[i] = 'zenith'
      2: beam_str[i] = 'off zenith'
    endcase
  endfor
  
  al_legend, ['weighted power ave ' + beam_str, 'straight power ave ' + beam_str,'0.52'],$
    textcolor = [colors2[0:nbeams-1], colors1[0:nbeams-1], 'black'], box = 0, /right, charsize = legend_charsize, charthick = charthick
    
  if keyword_set(pub) and n_elements(multi_pos) eq 0 then begin
    cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
  endif
  
  nrow = 2
  ncol = nbeams
  
  if keyword_set(pub) then plotfiles_use = plotfile_freq
  
  for i=0, nbeams-1 do begin
    if i eq 0 then begin
      start_multi_params = {ncol:ncol, nrow:nrow, ordering:'col'}
      undefine, positions, pos_use
      noerase = 0
    endif else pos_use = positions[*,2*i]
    
    quick_image, sim_ave_power_freq[i,*,*], total(reform(sim_ave_weights_freq[i,*,*]),2)/dims[1], $
      start_multi_params = start_multi_params, multi_pos = pos_use, $
      xtitle = 'ave weight', ytitle = 'frequency channel', title = 'straight power ave ' + beam_str[i], $
      png = png, eps = eps, pdf = pdf, alphabackgroundimage = alphabackgroundimage, savefile = plotfiles_use, noerase = noerase
      
    if i eq 0 then begin
      positions = pos_use
      undefine, start_multi_params
      noerase = 1
    endif
    
    pos_use = positions[*,2*i+1]
    quick_image, sim_wt_ave_power_freq[i,*,*], total(reform(sim_ave_weights_freq[i,*,*]),2)/dims[1], $
      start_multi_params = start_multi_params, multi_pos = pos_use, $
      xtitle = 'ave weight', ytitle = 'frequency channel', title = 'weighted power ave ' + beam_str[i], $
      png = png, eps = eps, pdf = pdf, alphabackgroundimage = alphabackgroundimage, savefile = plotfiles_use, noerase = noerase
      
  endfor
  
  if keyword_set(pub) then begin
    cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
  endif
  
  tvlct, r, g, b
end
