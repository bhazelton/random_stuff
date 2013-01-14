pro chi2_gen, fit, data, upper_err, lower_err, chi2, dof = dof, rchi2 = rchi2, $
              error_statistic = error_statistic, no_err = no_err

n_pts = n_elements(data)

if n_elements(fit) ne n_pts then message, $
  'data & fit arrays must have equal length.'

if not keyword_set(no_err) and (n_elements(upper_err) ne n_pts or $
                                n_elements(lower_err) ne n_pts) then begin
   default, error_statistic, 'poisson'
   allowed_statistics = ['poisson', 'gaussian']
   
   stat_ind = where(allowed_statistics eq error_statistic)
   if stat_ind[0] eq -1 then message, 'error_statistic not recognized'
   
   if error_statistic eq 'poisson' then begin
      poisson_errors, data, limits
      upper_err = limits.upper
      lower_err = limits.lower
   endif else begin
      upper_err = data + sqrt(data)
      lower_err = data - sqrt(data)
   endelse
endif


diff = fit - data

if not keyword_set(no_err) then begin
   mask = fltarr(n_pts)
   wh = where(diff ge 0, nwh)
   if nwh ge 1 then mask[wh] = 1
   
   sigma = [upper_err - data] * mask + [data - lower_err] * (1 - mask)
   
   chi2 = total((diff / sigma) ^ 2)
endif else begin
   chi2 = total(diff ^ 2 / fit)
endelse   

if n_elements(dof) gt 0 then begin
   if dof gt 0 then rchi2 = chi2 / dof $
   else rchi2 = -1
endif

end
