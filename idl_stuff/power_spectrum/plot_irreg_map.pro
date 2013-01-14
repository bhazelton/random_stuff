

pro plot_irreg_map, data_values, x_values, y_values, resolution, plot_sampling = plot_sampling, window_num = window_num, $
                    plotfile = plotfile, pub = pub, title = title, xtitle = xtitle, ytitle = ytitle, interp = interp

  if n_elements(window_num) eq 0 then window_num = 1

  dims = size(data_values, /dimension)
  if n_elements(dims) ne 1 then message, 'data_values must be a vector'

  npts = n_elements(data_values)
  if n_elements(x_values) ne npts or n_elements(y_values) ne npts then $
     message, 'x_values and y_values must have the same number of elements as data_values'

  case n_elements(resolution) of
     1: begin
        x_res = resolution[0]
        y_res = resolution[0]
     end
     2: begin
         x_res = resolution[0]
        y_res = resolution[1]
     end
     else: message, 'resolution must have 1 or 2 elements'
  endcase


  ;; get current colortable so it can be restored later
  tvlct, r, g, b, /get
  loadct,39

  
  ;; Work out plot & colorbar positions
  cb_size = 0.025               ;; in units of plot area (incl. margins)
  big_margin = [0.15, 0.2]      ;; in units of plot area (incl. margins)
  small_margin = [0.05, 0.1]    ;; in units of plot area (incl. margins)
  cb_big_margin = 0.08
  
  plot_pos = [big_margin[0], big_margin[1], (1-small_margin[0]-cb_big_margin[0]-cb_size), (1-small_margin[1])]
  cb_pos = [(1-small_margin[0]-cb_size), big_margin[1], (1-small_margin[0]), (1-small_margin[1])]
  
  plot_aspect = (plot_pos[3] - plot_pos[1]) / (plot_pos[2] - plot_pos[0])

  x_length = max(x_values) - min(x_values) + x_res
  y_length = max(y_values) - min(y_values) + y_res

  data_aspect = (y_length / x_length)
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
     
     window, window_num, xsize = 700 * x_factor, ysize = 700 * y_factor
    
     pson, file = plotfile, /eps 
  endif else begin
     if windowavailable(window_num) then wset, window_num else window, window_num, xsize = 700 * x_factor, ysize = 700 * y_factor
     charsize = 1
  endelse

  if n_elements(data_range) lt 2 then data_range = minmax(data_values)
  color_range = [0, 254]
  n_colors = color_range[1] - color_range[0]

  if keyword_set(interp) then begin
     ;; interpolate data to get a nice background image
     nx_grid = floor(x_length/x_res)*2
     ny_grid = floor(y_length/y_res)*2
     ;; input x&y values need to be in units of nx/ny_grid and
     ;; strictly less than nx/ny_grid
     x_vals_adj = (x_values-min(x_values))*(nx_grid-.01)/(max(x_values)-min(x_values))
     y_vals_adj = (y_values-min(y_values))*(ny_grid-.01)/(max(y_values)-min(y_values))
     
     interp_image = tsc(data_values, x_vals_adj, nx_grid, y_vals_adj, ny_grid, /average, /isolated)

     interp_image_norm = (interp_image-data_range[0])*n_colors/(data_range[1]-data_range[0]) + color_range[0]
     wh = where(interp_image_norm lt color_range[0], count)
     if count ne 0 then interp_image_norm[wh] = color_range[0]
     wh = where(interp_image_norm gt color_range[1], count)
     if count ne 0 then interp_image_norm[wh] = color_range[1]
  endif

  if n_elements(title) ne 0 then plot_title = title $
  else begin
     if keyword_set(plot_sampling) then plot_title = 'Image (sampled at white points) (Jy/beam)' $
     else plot_title = 'Image (Jy/beam)'
  endelse

  theta_char = greek('theta')
  if n_elements(xtitle) ne 0 then plot_xtitle = xtitle else plot_xtitle = theta_char+'x (deg)'
  if n_elements(ytitle) ne 0 then plot_ytitle = ytitle else plot_ytitle = theta_char+'y (deg)'


  loadct,39
  plot, x_values, y_values, /nodata, title = plot_title, position = plot_pos, $
        xrange = minmax(x_values), yrange = minmax(y_values), xstyle = 5, ystyle = 5, $
        thick = thick, charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font

  if keyword_set(interp) then begin
     cgimage, interp_image_norm, /nointerp,/overplot,/noerase
     axis, xaxis=0, xtick_get = xticks, xrange = minmax(x_values), xstyle = 1, xtitle = plot_xtitle, $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font
     axis, yaxis=0, ytick_get = yticks, yrange = minmax(y_values), ystyle = 1, ytitle = plot_ytitle, $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font
     axis, xaxis=1, xrange = minmax(x_values), xstyle = 1, xtickv = xticks, $
           xtickname = replicate(' ', n_elements(xticks)), charthick = charthick, xthick = xthick, ythick = ythick, $
           charsize = charsize, font = font
     axis, yaxis=1, yrange = minmax(y_values), ystyle = 1, ytickv = yticks, $
           ytickname = replicate(' ', n_elements(yticks)), charthick = charthick, xthick = xthick, ythick = ythick, $
           charsize = charsize, font = font
  endif else begin
     data_colors = round((data_values-data_range[0])*n_colors/(data_range[1]-data_range[0])) + color_range[0]

     plotsym, 0, 0.5, /fill
     for i=0L, npts-1 do oplot, [x_values[i]], [y_values[i]], psym = 8, color = data_colors[i]
  endelse


if keyword_set(plot_sampling) then oplot, x_values, y_values, psym = 3, color = 255
  

  exps = floor(alog10(abs(data_range)))

  temp = [ceil(data_range[0]/(10d^exps[0])*10)/10d*(10d^exps[0]), $
          floor(data_range[1]/(10d^exps[1])*10)/10d*(10d^exps[1])]

  n_div = 6
  div_size = (temp[1]-temp[0]) / n_div
  tick_vals = dindgen(n_div+1)*div_size + temp[0]

  names = string(tick_vals, format = '(e0.1)')

  if alog10(abs(tick_vals[0] - data_range[0])) gt 10^(-3d) then begin
     cb_ticknames = [' ', names]
     cb_ticks = [color_range[0], (tick_vals - data_range[0]) * n_colors / $
                 (data_range[1] - data_range[0]) + color_range[0]] - color_range[0]
  endif else begin
     cb_ticknames = names
     cb_ticks = ((tick_vals - data_range[0]) * n_colors / $
                 (data_range[1] - data_range[0]) + color_range[0]) - color_range[0]
  endelse

  if alog10(abs(data_range[1] - max(tick_vals))) gt 10^(-3d) then begin
     cb_ticknames = [cb_ticknames, ' ']
     cb_ticks = [cb_ticks, color_range[1]-color_range[0]+1]
  endif

  cgcolorbar, color = 0, /vertical, position = cb_pos, bottom = color_range[0], ncolors = n_colors+1, yminor = 0, $
              ticknames = cb_ticknames, ytickv = cb_ticks, yticks = n_elements(cb_ticks) -1, charsize = charsize, font = font

  if keyword_set(pub) then begin
     psoff
     wdelete, window_num
  endif
  
  tvlct, r, g, b

end
