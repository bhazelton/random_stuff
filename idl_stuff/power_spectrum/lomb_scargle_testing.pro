

pro lomb_scargle_testing, grid_type = grid_type, dbl_pi = dbl_pi

  if not keyword_set(grid_type) then grid_type = 'even'
  if keyword_set(dbl_pi) then pi = !dpi else pi = double(!pi)  


  case grid_type of
     'even': begin
        num=30
        spacing1 = 1d/num
        spacing2 = 1d/num
        
        locs1 = rebin(dindgen(num)*spacing1, num, num)
        locs2 = rebin(reform(dindgen(num)*spacing2, 1, num), num, num)
        n_pts = n_elements(locs1)
        
        locs1 = reform(locs1, n_pts)
        locs2 = reform(locs2, n_pts)

     end
     
     'squished': begin
        num=30
        spacing1 = dindgen(num)*0.1d
        spacing2 = 1
        
        locs1 = rebin(dindgen(num)*spacing1, num, num)
        locs2 = rebin(reform(dindgen(num)*spacing2, 1, num), num, num)
        n_pts = n_elements(locs1)

        locs1 = reform(locs1, n_pts)
        locs2 = reform(locs2, n_pts)
     end

     'fan': begin
        num = 30
        spacing1 = 1
        spacing2 = 1
        
        locs1 = rebin(dindgen(num)*spacing1, num, num)
        locs2 = rebin(reform(dindgen(num)*spacing2, 1, num), num, num)

        spacing1 = dblarr(num-1, num)
        for i=0, num -1 do begin
           locs1[*,i] = locs1[*,i]*(1d + 0.1d*i)
           temp = locs1[*,i] - shift(locs1[*,i], 1)
           temp = temp[1:*]
           spacing1[*,i] = temp
        endfor

        n_pts = n_elements(locs1)
        locs1 = reform(locs1, n_pts)
        locs2 = reform(locs2, n_pts)

     end
     
     'random' : begin
        size=30
        num = 100

        seed = systime(1)
        locs1 = randomu(seed, num^2)*(size-1)
        locs2 = randomu(seed, num^2)*(size-1)
        n_pts = n_elements(locs1)

        spacing1 = (spacing2 = double(size)/double(num))

     end

     else: message, 'invalid grid type'
  endcase

  vals = dblarr(n_pts)
  ;; for i=0, n_pts-1 do vals[i] = sin(5*pi*locs1[i]/max(locs1))
  for i=0, n_pts-1 do vals[i] = sin(5*pi*locs1[i]/max(locs1)) * cos(3*pi*locs2[i]/max(locs2))

  !p.multi = [0,2,3]
  ;; get current color table so we can restore it.
  tvlct, r, g, b, /get
  loadct,3
  if windowavailable(1) then wset, 1 else window, 1,  xsize = 300*2, ysize = 300*3

  data_colors = round(254*(vals-min(vals))/(max(vals) - min(vals)))
  plot, locs1, locs2, psym = 3, color = 0, xtitle = 'x', ytitle = 'y', title = 'Image (sampled at white points)', $
        xrange = minmax(locs1), yrange = minmax(locs2), xstyle = 1, ystyle = 1
  plotsym, 0, 0.8, /fill
  for i=0L, n_pts-1 do oplot, [locs1[i]], [locs2[i]], psym = 8, color = data_colors[i]
  oplot, locs1, locs2, psym = 3, color = 255

  if grid_type eq 'even' then begin
     ft = fft(reform(vals, 30,30))
     ft = ft[0:15, 0:15]
     ft_power = real_part(ft * conj(ft)) * n_elements(locs1)
     
     ian_ls = ian_lomb_scargle_grid(reform(vals, 30,30), delta1 = ian_ls_delta1, delta2 = ian_ls_delta2)
     ian_ls = ian_ls[0:15, 0:15]
     ian_ls_power = real_part(ian_ls * conj(ian_ls))* n_elements(locs1)
  endif

  ;; work out k values
  k1_max = (2d*pi)/(2*spacing1)
  k2_max = (2d*pi)/(2*spacing2)
 
  delta_k1 = (2d*pi)/(max(locs1) - min(locs1)+spacing1)
  delta_k2 = (2d*pi)/(max(locs2) - min(locs2)+spacing2)

  ;; define locations (in k) to take FT
  n_k1 = round(k1_max/delta_k1) + 1
  n_k2 = round(k2_max/delta_k2) + 1
  k1_locs = indgen(n_k1) * delta_k1
  k2_locs = indgen(n_k2) * delta_k2

  if grid_type eq 'even' then begin
     ls_delta1 = dblarr(n_k1)
     ls_delta2 = dblarr(n_k2)
  endif


  lst = lomb_scargle_2d_fast(locs1, locs2, vals, k1_locs/(2d*pi), k2_locs/(2d*pi), /transform, delta1 = ls_delta1, $
                             delta2 = ls_delta2, dbl_pi = dbl_pi)
  lsp = lomb_scargle_2d_fast(locs1, locs2, vals, k1_locs/(2d*pi), k2_locs/(2d*pi), dbl_pi = dbl_pi)

  dft = discrete_ft_2d(locs1, locs2, vals, k1_locs, k2_locs)
  dft_power = real_part(dft * conj(dft))/n_elements(locs1)

  lst_power = real_part(lst * conj(lst)) / n_elements(locs1)
  diff = lsp - lst_power

  diff_dft = dft_power - lsp;;lst_power

  temp = [lsp, lst_power, abs(diff), dft_power, abs(diff_dft)]
  temp = temp(where(temp ne 0d))
  range = minmax(alog10(temp))
  lsp_norm =  (alog10(lsp)-range[0])*255/(range[1]-range[0])
  lst_power_norm =  (alog10(lst_power)-range[0])*255/(range[1]-range[0])
  diff_norm =  (alog10(abs(diff))-range[0])*255/(range[1]-range[0])
  dft_power_norm = (alog10(dft_power)-range[0])*255/(range[1]-range[0])
  dft_diff_norm = (alog10(abs(diff_dft))-range[0])*255/(range[1]-range[0])

  ;; range = minmax(temp)
  ;; lsp_norm =  (lsp-range[0])*255/(range[1]-range[0])
  ;; lst_power_norm =  (lst_power-range[0])*255/(range[1]-range[0])
  ;; diff_norm =  (diff-range[0])*255/(range[1]-range[0])
  ;; dft_power_norm = (dft_power-range[0])*255/(range[1]-range[0])
  ;; dft_diff_norm = (diff_dft-range[0])*255/(range[1]-range[0])

  inc_factor = 100
  images = {lsp: congrid(lsp_norm, n_elements(k1_locs)*inc_factor, n_elements(k2_locs)*inc_factor), $
            lst_power: congrid(lst_power_norm, n_elements(k1_locs)*inc_factor, n_elements(k2_locs)*inc_factor), $
            diff: congrid(diff_norm, n_elements(k1_locs)*inc_factor, n_elements(k2_locs)*inc_factor), $
            dft_power: congrid(dft_power_norm, n_elements(k1_locs)*inc_factor, n_elements(k2_locs)*inc_factor), $
            dft_diff: congrid(dft_diff_norm, n_elements(k1_locs)*inc_factor, n_elements(k2_locs)*inc_factor)}
  titles = ['log LSP', 'log LST Squared', 'log abs(LSP-LST^2)', 'log DFT power', 'log abs(DFT^2) - LSP']


  for i=0, 4 do begin
     plot, k1_locs, k2_locs, /nodata, xstyle=4, ystyle=4, title = titles[i], thick = thick, charthick = charthick, $
           xthick = xthick, ythick = ythick, charsize = charsize, font = font
     tvimage, images.(i), /nointerp,/overplot
     axis, xaxis=0, xtick_get = xticks, xtitle = 'kx', xrange = minmax(k1_locs), xstyle = 1, $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font
     axis, yaxis=0, ytick_get = yticks, ytitle = 'ky', yrange = minmax(k2_locs), ystyle = 1, $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font
     axis, xaxis=1, xrange = minmax(k1_locs), xstyle = 1, xtickv = xticks, xtickname = replicate(' ', n_elements(xticks)), $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font
     axis, yaxis=1, yrange = minmax(k2_locs), ystyle = 1, ytickv = yticks, ytickname = replicate(' ', n_elements(yticks)), $
           charthick = charthick, xthick = xthick, ythick = ythick, charsize = charsize, font = font
  endfor

  !p.multi=0
  tvlct, r, g, b

stop








end
