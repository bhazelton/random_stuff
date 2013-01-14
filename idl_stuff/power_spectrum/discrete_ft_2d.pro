; Straight forward 2D discrete fourier transform of unevenly spaced data
; locations1 & 2 are x/y values of data points
; fourier_locs1 & 2 are kx/ky values to test at

function discrete_ft_2d, locations1, locations2, data, k1, k2, timing = timing

  time0 = systime(1)
  if n_elements(locations1) ne n_elements(data) or n_elements(locations2) ne n_elements(data) then message, $
     'locations1 & 2 must have same number of elements as data.'

  n_k1 = n_elements(k1)
  n_k2 = n_elements(k2)
  n_pts = n_elements(data)

  ft = dcomplex(dblarr(n_k1, n_k2))
  for i=0, n_k1-1 do begin
     x_exp = exp(-1*dcomplex(0,1)*k1[i]*locations1)
     for j=0, n_k2-1 do begin
        y_exp = exp(-1*dcomplex(0,1)*k2[j]*locations2)
        ft[i,j] = total(data*x_exp*y_exp)
     endfor
  endfor
  time1= systime(1)
  timing = time1-time0
  
  return, ft
end
