pro fitting, x, y, result, chisq, rchisq, y_calc, status = status, y_err = y_err, no_err = no_err, linear = linear, poly = poly, $
             degree = degree, init_exp = init_exp, nterms_exp = nterms_exp, gauss = gauss, nterms_gauss = nterms_gauss, $
             init_gauss = init_gauss, init_lor = init_lor, nterms_lor = nterms_lor, init_hybrid = init_hybrid, $
             init_power = init_power, itmax = itmax, outlier_rej = outlier_rej, out_rej_pc = out_rej_pc, xy_mod = xy_mod, rej = rej, $
             quiet = quiet

;y_err: structure with upper and lower errors for y array (form:
;y_err.upper, yerr.lower). If not set, y errors are set to poisson errors.

;outlier_rej: sets oulier rejection
;out_rej_pc: sets % of database to be rejected (default=5%)
;xy_mod: x & y arrays after outlier rejection
;rej: rejected data points

n = n_elements(y)
out_rej_pc = 5
rej_num = floor(n*out_rej_pc/100)

status = 0

if keyword_set(y_err) and keyword_set(no_err) then $
   message, 'no_err keyword set but error values provided'

if not keyword_set(y_err) then begin
   if min(y) lt 0 then no_err = 1

   if not keyword_set(no_err) then begin
      poisson_errors, y, y_err

      sym_errs = dblarr(n)
      for i = 0, n - 1 do sym_errs[i] = mean([y_err.upper[i] - y[i], $
                                              y[i] - y_err.lower[i]])
   endif else begin
      y_err = {upper:-1, lower:-1}
   endelse
endif

if keyword_set(linear) then begin
    test=0
    n_rej=0
    y_mod=y
    ymod_err=y_err
    x_mod=x
    xymod_lin=[[x_mod],[y_mod]]
    if keyword_set(outlier_rej) then lin_rej=dblarr(rej_num, 2) else lin_rej=0
    while test eq 0 do begin
        lin_res=linfit(x_mod, y_mod, /double, measure_errors=sym_errs)
        lin_calc=lin_res[0]+lin_res[1]*x_mod
        
        
        chi2_gen, lin_calc, y_mod, ymod_err.upper, ymod_err.lower, lin_chisq, $
                  dof = n_elements(lin_calc) - n_elements(lin_res), $
                  rchi2 = lin_rchisq, no_err = no_err

        if keyword_set(outlier_rej) and n_rej lt rej_num then begin
            index=where(lin_chisq eq max(lin_chisq))
            mask=dblarr(n_elements(y_mod))+1
            mask[index]=0
            lin_rej[n_rej,0]=x_mod[index]
            lin_rej[n_rej,1]=y_mod[index]
            y_mod=y_mod[where(mask eq 1)]
            ymod_err={upper:ymod_err.upper[where(mask eq 1)], $
                      lower:ymod_err.lower[where(mask eq 1)]}
            x_mod=x_mod[where(mask eq 1)]
            xymod_lin=[[x_mod],[y_mod]]
            n_rej=n_rej+1
        endif else test=1
    endwhile
endif else begin
    lin_res=0
    lin_chisq=0
    lin_rchisq=0
    lin_calc=0
    xymod_lin=0
    lin_rej=0
endelse

if keyword_set(poly) then begin
    default, degree, 2
    test=0
    n_rej=0
    y_mod=y
    ymod_err=y_err
    x_mod=x
    xymod_poly=[[x_mod],[y_mod]]
    if keyword_set(outlier_rej) then poly_rej=dblarr(rej_num, 2) $
    else poly_rej=0
    while test eq 0 do begin
        poly_res=svdfit(x_mod, y_mod, degree+1, /double, $
                       measure_errors=sym_errs)
        poly_calc=poly_res[0]+x_mod*0
        for i=1, degree do poly_calc=poly_calc+poly_res[i]*x_mod^i

        chi2_gen, poly_calc, y_mod, ymod_err.upper, ymod_err.lower, $
                  poly_chisq, dof=n_elements(poly_calc)-n_elements(poly_res), $
                  rchi2=poly_rchisq, no_err = no_err

        if keyword_set(outlier_rej) and n_rej lt rej_num then begin
            index=where(poly_chisq eq max(poly_chisq))
            mask=dblarr(n_elements(y_mod))+1
            mask[index]=0
            poly_rej[n_rej,0]=x_mod[index]
            poly_rej[n_rej,1]=y_mod[index]
            y_mod=y_mod[where(mask eq 1)]
            ymod_err={upper:ymod_err.upper[where(mask eq 1)], $
                      lower:ymod_err.lower[where(mask eq 1)]}
            x_mod=x_mod[where(mask eq 1)]
            xymod_poly=[[x_mod],[y_mod]]
            n_rej=n_rej+1
        endif else test=1
    endwhile
endif else begin
    poly_res=0
    poly_chisq=0
    poly_rchisq=0
    poly_calc=0
    xymod_poly=0
    poly_rej=0
endelse


if keyword_set(init_exp) then begin
    default, nterms_exp, 2
    if nterms_exp lt 2 or nterms_exp gt 5 then message, $
      'nterms_exp must be between 2 & 5'
    if n_elements(init_exp) ne 2 then $
      message, 'exp needs 2 parameters: y = a0 * exp(a1 * x)'
    a = init_exp
    default, itmax, 50
    test = 0
    n_rej = 0
    y_mod = y
    ymod_err = y_err
    x_mod = x
    xymod_exp = [[x_mod],[y_mod]]
    if keyword_set(outlier_rej) then exp_rej = dblarr(rej_num, 2) $
    else exp_rej=0
    while test eq 0 do begin
        weights=1.0/(sym_errs^2)
        case nterms_exp of
            2: exp_calc = curvefit(x, y, weights, a, itmax = itmax, $
                                   function_name = 'exponential', $
                                   status = status)
            3: begin
                a = [a, 0]
                exp_calc = curvefit(x, y, weights, a, itmax = itmax, $
                                    function_name = 'exponential3', $
                                    status = status)
            end
            4: begin
                a = [a, 0, 0]
                exp_calc = curvefit(x, y, weights, a, itmax = itmax, $
                                    function_name = 'exponential4', $
                                    status = status)
            end
            5: begin
                a = [a, 0, 0, 0]
                exp_calc = curvefit(x, y, weights, a, itmax = itmax, $
                                    function_name = 'exponential5', $
                                    status = status)
            end
        endcase
        exp_res = a

        chi2_gen, exp_calc, y_mod, ymod_err.upper, ymod_err.lower, exp_chisq, $
                  dof = n_elements(exp_calc)-n_elements(exp_res), $
                  rchi2 = exp_rchisq, no_err = no_err

        if keyword_set(outlier_rej) and n_rej lt rej_num then begin
            index=where(exp_chisq eq max(exp_chisq))
            mask=dblarr(n_elements(y_mod))+1
            mask[index]=0
            exp_rej[n_rej,0]=x_mod[index]
            exp_rej[n_rej,1]=y_mod[index]
            y_mod=y_mod[where(mask eq 1)]
            ymod_err={upper:ymod_err.upper[where(mask eq 1)], $
                      lower:ymod_err.lower[where(mask eq 1)]}
            x_mod=x_mod[where(mask eq 1)]
            xymod_exp=[[x_mod],[y_mod]]
            n_rej=n_rej+1
        endif else test=1
    endwhile
endif else begin
    exp_res = 0
    exp_chisq = 0
    exp_rchisq = 0
    exp_calc = 0
    xymod_exp = 0
    exp_rej = 0
endelse

if keyword_set(gauss) then begin
    if keyword_set(init_gauss) then begin
        if not keyword_set(nterms_gauss) then $
          nterms_gauss = n_elements(init_gauss) $
        else if nterms_gauss ne n_elements(init_gauss) then message, $
          'nterms_gauss must equal n_elements(init_gauss)'
        estimates=init_gauss
    endif
    default, nterms_gauss, 3
    if nterms_gauss lt 3 or nterms_gauss gt 6 then message, $
      'nterms_gauss must be between 3 & 6'
    test=0
    n_rej=0
    y_mod=y
    ymod_err=y_err
    x_mod=x
    xymod_gauss=[[x_mod],[y_mod]]
    if keyword_set(outlier_rej) then gauss_rej=dblarr(rej_num, 2) $
    else gauss_rej=0
    while test eq 0 do begin
        gauss_calc=gaussfit(x_mod, y_mod, gauss_res, measure_errors=sym_errs, $
                           nterms=nterms_gauss, estimates=estimates)

        chi2_gen, gauss_calc, y_mod, ymod_err.upper, ymod_err.lower, $
                  gauss_chisq, dof = n_elements(gauss_calc) - $
                  n_elements(gauss_res), rchi2=gauss_rchisq, no_err = no_err

        if keyword_set(outlier_rej) and n_rej lt rej_num then begin
            index=where(gauss_chisq eq max(gauss_chisq))
            mask=dblarr(n_elements(y_mod))+1
            mask[index]=0
            gauss_rej[n_rej,0]=x_mod[index]
            gauss_rej[n_rej,1]=y_mod[index]
            y_mod=y_mod[where(mask eq 1)]
            ymod_err={upper:ymod_err.upper[where(mask eq 1)], $
                      lower:ymod_err.lower[where(mask eq 1)]}
            x_mod=x_mod[where(mask eq 1)]
            xymod_gauss=[[x_mod],[y_mod]]
            n_rej=n_rej+1
        endif else test=1
    endwhile
endif else begin
    gauss_res=0
    gauss_chisq=0
    gauss_rchisq=0
    gauss_calc=0
    xymod_gauss=0
    gauss_rej=0
endelse

if keyword_set(init_lor) then begin
    default, nterms_lor, 3
    if nterms_lor lt 3 or nterms_lor gt 6 then message, $
      'nterms_lor must be between 3 & 6'
    if n_elements(init_lor) ne 3 then $
      message, 'lorentzian needs 3 parameters, a0=peak height,' + $
      ' a1=peak location, a2=FWHM/2'
    a = init_lor
    default, itmax, 50
    test = 0
    n_rej = 0
    y_mod = y
    ymod_err = y_err
    x_mod = x
    xymod_lor = [[x_mod],[y_mod]]
    if keyword_set(outlier_rej) then lor_rej=dblarr(rej_num, 2) $
    else lor_rej=0
    while test eq 0 do begin
        weights=1.0/(sym_errs^2)

        case nterms_lor of
            3: lor_calc = curvefit(x, y, weights, a, itmax = itmax, $
                                   function_name = 'lorentzian', $
                                   status = status)
            4: begin
                a = [a, 0]
                lor_calc = curvefit(x, y, weights, a, itmax = itmax, $
                                    function_name = 'lorentzian4', $
                                    status = status)
            end
            5: begin
                a = [a, 0, 0]
                lor_calc = curvefit(x, y, weights, a, itmax = itmax, $
                                    function_name = 'lorentzian5', $
                                    status = status)
            end
            6: begin
                a = [a, 0, 0, 0]
                lor_calc = curvefit(x, y, weights, a, itmax = itmax, $
                                    function_name = 'lorentzian6', $
                                    status = status)
            end
        endcase
        lor_res=a

        chi2_gen, lor_calc, y_mod, ymod_err.upper, ymod_err.lower, lor_chisq, $
                  dof = n_elements(lor_calc) - n_elements(lor_res), $
                  rchi2 = lor_rchisq, no_err = no_err

        if keyword_set(outlier_rej) and n_rej lt rej_num then begin
            index=where(lor_chisq eq max(lor_chisq))
            mask=dblarr(n_elements(y_mod))+1
            mask[index]=0
            lor_rej[n_rej,0]=x_mod[index]
            lor_rej[n_rej,1]=y_mod[index]
            y_mod=y_mod[where(mask eq 1)]
            ymod_err={upper:ymod_err.upper[where(mask eq 1)], $
                      lower:ymod_err.lower[where(mask eq 1)]}
            x_mod=x_mod[where(mask eq 1)]
            xymod_lor=[[x_mod],[y_mod]]
            n_rej=n_rej+1
        endif else test=1
    endwhile
endif else begin
    lor_res=0
    lor_chisq=0
    lor_rchisq=0
    lor_calc=0
    xymod_lor=0
    lor_rej=0
endelse


if keyword_set(init_hybrid) then begin
    if n_elements(init_hybrid) ne 6 then $
      message, 'hybrid lorentzian/gauss needs 6 parameters,' + $
      'a0=peak height/2, a1=peak location, a2=sigma, a3=peak height/2,' + $
      'a4=peak location, a5=FWHM/2'
    a = [init_hybrid, 0]
    default, itmax, 50
    test = 0
    n_rej = 0
    y_mod = y
    ymod_err = y_err
    x_mod = x
    xymod_hybrid = [[x_mod],[y_mod]]
    if keyword_set(outlier_rej) then hybrid_rej=dblarr(rej_num, 2) $
    else hybrid_rej=0
    while test eq 0 do begin
        weights=1.0/(sym_errs^2)

        hybrid_calc = curvefit(x, y, weights, a, itmax = itmax, $
                               function_name = 'lor_gauss_hybrid', $
                               status = status)
        hybrid_res = a

        chi2_gen, hybrid_calc, y_mod, ymod_err.upper, ymod_err.lower, $
                  hybrid_chisq, dof = n_elements(hybrid_calc) - $
                  n_elements(hybrid_res), rchi2 = hybrid_rchisq, no_err = no_err

        if keyword_set(outlier_rej) and n_rej lt rej_num then begin
            index=where(hybrid_chisq eq max(hybrid_chisq))
            mask=dblarr(n_elements(y_mod))+1
            mask[index]=0
            hybrid_rej[n_rej,0]=x_mod[index]
            hybrid_rej[n_rej,1]=y_mod[index]
            y_mod=y_mod[where(mask eq 1)]
            ymod_err={upper:ymod_err.upper[where(mask eq 1)], $
                      lower:ymod_err.lower[where(mask eq 1)]}
            x_mod=x_mod[where(mask eq 1)]
            xymod_hybrid=[[x_mod],[y_mod]]
            n_rej=n_rej+1
        endif else test=1
    endwhile
endif else begin
    hybrid_res=0
    hybrid_chisq=0
    hybrid_rchisq=0
    hybrid_calc=0
    xymod_hybrid=0
    hybrid_rej=0
endelse


if keyword_set(init_power) then begin
    if n_elements(init_power) ne 2 then $
      message, 'power law needs 2 parameters: y = a0 * x ^ a1)'
    a = init_power
    default, itmax, 50
    test = 0
    n_rej = 0
    y_mod = y
    ymod_err = y_err
    x_mod = x
    xymod_power = [[x_mod],[y_mod]]
    if keyword_set(outlier_rej) then power_rej = dblarr(rej_num, 2) $
    else power_rej=0
    while test eq 0 do begin
        weights = 1.0 / (sym_errs ^ 2)

        power_calc = curvefit(x, y, weights, a, itmax = itmax, $
                                   function_name = 'power_law', $
                                   status = status)
        power_res = a

        chi2_gen, power_calc, y_mod, ymod_err.upper, ymod_err.lower, $
                  power_chisq, dof = n_elements(power_calc) - $
                  n_elements(power_res), rchi2 = power_rchisq, no_err = no_err

        if keyword_set(outlier_rej) and n_rej lt rej_num then begin
            index = where(power_chisq eq max(power_chisq))
            mask = dblarr(n_elements(y_mod)) + 1
            mask[index] = 0
            power_rej[n_rej,0] = x_mod[index]
            power_rej[n_rej,1] = y_mod[index]
            y_mod = y_mod[where(mask eq 1)]
            ymod_err = {upper:ymod_err.upper[where(mask eq 1)], $
                        lower:ymod_err.lower[where(mask eq 1)]}
            x_mod = x_mod[where(mask eq 1)]
            xymod_power = [[x_mod], [y_mod]]
            n_rej = n_rej + 1
        endif else test = 1
    endwhile
endif else begin
    power_res = 0
    power_chisq = 0
    power_rchisq = 0
    power_calc = 0
    xymod_power = 0
    power_rej = 0
endelse



result={lin:lin_res, poly:poly_res, exp:exp_res, gauss:gauss_res, lor:lor_res,$
       hybrid:hybrid_res, power:power_res}

chisq={lin:lin_chisq, poly:poly_chisq, exp:exp_chisq, gauss:gauss_chisq, $
       lor:lor_chisq, hybrid:hybrid_chisq, power:power_chisq}

;reduced chisq:
rchisq={lin:lin_rchisq, poly:poly_rchisq, exp:exp_rchisq, gauss:gauss_rchisq, $
       lor:lor_rchisq, hybrid:hybrid_rchisq, power:power_rchisq}

y_calc={lin:lin_calc, poly:poly_calc, exp:exp_calc, gauss:gauss_calc, $
        lor:lor_calc, hybrid:hybrid_calc, power:power_calc}

xy_mod={lin:xymod_lin, poly:xymod_poly, exp:xymod_exp, gauss:xymod_gauss, $
       lor:xymod_lor, hybrid:xymod_hybrid, power:xymod_power}

rej={lin:lin_rej, poly:poly_rej, exp:exp_rej, gauss:gauss_rej, lor:lor_rej, $
     hybrid:hybrid_rej, power:power_rej}

if not keyword_set(quiet) then begin
   if chisq.lin ne 0 then print, 'linear params      = ',result.lin
   if chisq.poly ne 0 then print, 'polynomial params = ',result.poly
   if chisq.exp ne 0 then print, 'exponential params = ',result.exp
   if chisq.gauss ne 0 then print, 'gaussian params = ',result.gauss
   if chisq.lor ne 0 then print, 'lorentzian params = ',result.lor
   if chisq.hybrid ne 0 then print, 'gauss/lorentz hybrid params = ',$
                                    result.hybrid
   if chisq.power ne 0 then print, 'power law params = ',result.power

   print,''
   if chisq.lin ne 0 then print, 'linear chisq, reduced chisq      = ', $
                                 chisq.lin, ',', rchisq.lin
   if chisq.poly ne 0 then print,'polynomial chisq, reduced chisq  = ', $
                                 chisq.poly, ',', rchisq.poly
   if chisq.exp ne 0 then print, 'exponential chisq, reduced chisq = ', $
                                 chisq.exp, ',', rchisq.exp
   if chisq.gauss ne 0 then print, 'gaussian chisq, reduced chisq = ', $
                                   chisq.gauss, ',', rchisq.gauss
   if chisq.lor ne 0 then print, 'lorentzian chisq, reduced chisq = ', $
                                 chisq.lor, ',', rchisq.lor
   if chisq.hybrid ne 0 then print, $
      'gauss/lorentz hybrid chisq, reduced chisq = ', chisq.hybrid, ',', $
      rchisq.hybrid
   if chisq.power ne 0 then print, 'power law chisq, reduced chisq = ', $
                                   chisq.power, ',', rchisq.power
endif
end
