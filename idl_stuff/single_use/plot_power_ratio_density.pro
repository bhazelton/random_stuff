pro plot_power_ratio_density, png = png, eps = eps, pdf = pdf, xlog=xlog, aspect=aspect

  if n_elements(aspect) eq 0 then aspect = 0.5

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

  plotfile = '/Users/bryna/Projects/Physics/data_files/fhd_ps_data/sim_suite_plots/sim_suite_power_ratio_paper'

  endif

  data_file = '/Users/bryna/Projects/Physics/data_files/fhd_ps_data/sim_suite_plots/arrsim_flat_realsky_cube_power_info.idlsave'

  restore, data_file

  min_ind_use = 7
  inds_use = indgen(12-min_ind_use) + min_ind_use

  if keyword_set(pub) then begin
    if keyword_set(xlog) then plotfile = plotfile + '_xlog'
    plotfile = plotfile + plot_exten

    charthick = 3
    thick = 3
    xthick = 3
    ythick = 3
    charsize = 2
    font = 1
    legend_charsize = 2
    symsize = 2

  cgps_open, plotfile, /font, encapsulated=eps, landscape=1, pagetype='letter'
  endif

  if keyword_set(xlog) then begin
    weight_plot_range = [0.0008, 0.03]
  endif else begin
    weight_plot_range = [0, 0.025]
  endelse
  asymptote_line = 0.5+intarr(2)

  cgplot, weight_plot_range, asymptote_line, yrange=[0,1], $
    ytitle = 'Power Ratio', xtitle = 'Weights', xlog=xlog, xrange=weight_plot_range, $
    thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, $
    charsize = charsize, font = font, position = position, aspect = aspect, linestyle=2
  cgplot, sim_ave_weights[inds_use], sim_wt_ave_powers[inds_use]/flat_power, $
    /over, color='blue', psym=-4, thick = thick, symsize=symsize

  if keyword_set(pub) and n_elements(multi_pos) eq 0 then begin
    cgps_close, png = png, pdf = pdf, delete_ps = delete_ps, density=600
  endif

end
