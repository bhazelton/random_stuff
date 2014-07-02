; Vectorized 2D discrete fourier transform of unevenly spaced data
; locations1 & 2 are x/y values of data points
; fourier_locs1 & 2 are kx/ky values to test at


function discrete_ft_2d_vec, locations1, locations2, data, k1, k2, timing = timing

  time0 = systime(1)
  if n_elements(locations1) ne n_elements(data) or n_elements(locations2) ne n_elements(data) then message, $
     'locations1 & 2 must have same number of elements as data.'

  n_k1 = n_elements(k1)
  n_k2 = n_elements(k2)
  n_pts = n_elements(data)

  ;; x_exp = exp(-1 * 2 * !pi * dcomplex(0,1) * rebin(reform(fourier_locs1, 1, n_k1), n_pts, n_k1, n_k2) * $
  ;;                                                  rebin(locations1, n_pts, n_k1, n_k2))
  ;; y_exp = exp(-1 * 2 * !pi * dcomplex(0,1) * rebin(reform(fourier_locs2, 1, 1, n_k2), n_pts, n_k1, n_k2) * $
  ;;             rebin(locations2, n_pts, n_k1, n_k2))
  print, 'generating x exponential array' 
  x_exp = exp(-1*dcomplex(0,1) * rebin(rebin(reform(k1, 1, n_k1), n_pts, n_k1) * $
                                       rebin(locations1, n_pts, n_k1), n_pts, n_k1, n_k2))
  print, 'generating y exponential array' 
  y_exp = exp(-1*dcomplex(0,1) * rebin(reform(rebin(reform(k2, 1, n_k2), n_pts, n_k2) * $
                                              rebin(locations2, n_pts, n_k2), npts, 1, n_k2), n_pts, n_k1, n_k2))
  print, 'multiplying and summing' 
  ft = total(rebin(data, n_pts, n_k1, n_k2)*temporary(x_exp)*temporary(y_exp), 1)

  time1= systime(1)
  timing = time1-time0
  
  return, ft
end
