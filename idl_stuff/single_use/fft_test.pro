pro fft_test, fn_type = fn_type


  n_pts = 10000
  x_vals = dindgen(n_pts)/n_pts * 40-20
  x_delta = x_vals[1]-x_vals[0]
  k_vals = dindgen(n_pts)/n_pts * 2d*!dpi/x_delta -1d*!dpi/x_delta
  
  if n_elements(fn_type) eq 0 then fn_type = 'gaussian'

  case fn_type of 
     'gaussian': begin
        sigma = 2
        fx = exp(-1d*x_vals^2d/(2d*sigma^2d))
        integral = sqrt(2d*!dpi) * sigma * exp(-1d*k_vals^2d*sigma^2d/2d)
     end
     'exp_absx': begin
        a=3d
        fx = exp(-1d*a*abs(x_vals))
        integral = (2d * a)/ (a^2d + k_vals^2d)
     end
     'sinx_x': begin
        a=3d
        fx = sin(a*x_vals)/x_vals
        wh_x0 = where(x_vals eq 0, count_x0)
        if count_x0 gt 0 then fx[wh_x0] = 0

        integral = dblarr(n_pts)+ !dpi
        wh_abs_gta = where(abs(k_vals) gt a, count_abs_gta)
        if count_abs_gta gt 0 then integral[wh_abs_gta] = 0
     end
  endcase
  analytic = integral / (2d*!dpi)
  
  basic_fft = shift(fft(fx), n_pts/2)
  fft = (1/(2d*!dpi)) * (basic_fft * n_pts) * x_delta

  exp_arr = exp(-1d * dcomplex(0,1) * rebin(x_vals, n_pts, n_pts) * rebin(reform(k_vals, 1, n_pts), n_pts, n_pts))

  dft = (1/(2d*!dpi)) * total(rebin(fx, n_pts, n_pts) * exp_arr, 1) * x_delta


  plot, k_vals, analytic, xrange = [-10,10], yrange = [0, max(analytic)* 2]
  oplot, k_vals, abs(fft), color=254, linestyle=3
    
  oplot, k_vals, abs(dft), color=75, linestyle=2
  
  ;;stop
  
end
