pro lnsrch2, xold, fold, g, p, x, f, stpmax, check, func

  tol = 1e-7
  alf = 1e-4
  
  check=0
  
  p_use=p
  if sqrt(total(p^2.)) gt stpmax then p_use = p*stpmax/sqrt(total(p^2.))
  
  slope = total(g*p_use)
  if slope ge 0 then stop;message, 'Roundoff problem in lnsrch'
  
  p=p_use
  
  temp = abs(p)/(abs(xold) > 1)
  test = max(temp)
  
  alamin = tol/test
  alam = 1.
  
  while(1) do begin
    x = xold + alam*p
    f = call_function(func, x)
    if alam lt alamin then begin
      x = xold
      check=1
      return
    endif else if f le (fold+alf*alam*slope) then begin
      return
    endif else begin
      if alam eq 1 then tmplam = (-1. * slope)/(2.*(f-fold-slope)) else begin
        rhs1 = f-fold-alam*slope
        rhs2 = f2-fold-alam2*slope
        a = (rhs1/alam^2. - rhs2/alam2^2.)/(alam-alam2)
        b = (-1.*alam2*rhs1/alam^2. + alam*rhs2/alam2^2.)/(alam-alam2)
        
        if a eq 0 then tmplam = -1*slope/(2*b) else begin
          disc = b^2.-3.*a*slope
          if disc lt 0 then tmplam = 0.5*alam else if b le 0 then tmplam = (-1.*b+sqrt(disc))/(3.*a) else tmplam = -1*slope/(b+sqrt(disc))
        endelse
        if tmplam gt 0.5*alam then tmp = 0.5 * alam
      endelse
    endelse
    alam2 = alam
    f2 = f
    alam = tmplam > 0.1*alam
  endwhile
  
end

pro dfpmin2, p, gtol, fmin, func, dfunc, double=double, eps=eps, iter=iter, itmax=itmax, stepmax = stepmax

  if n_elements(itmax) eq 0 then itmax=200
  if n_elements(stepmax) eq 0 then stepmax=100
  if n_elements(eps) eq 0 then eps=3e-8
  tolx = 4 * eps
  
  nparams = n_elements(p)
  
  fp = call_function(func, p)
  g = call_function(dfunc, p)
  hessian = fltarr(nparams, nparams)
  hessian[lindgen(nparams) * (nparams+1)] = 1.0
  xi = -1.*g
  
  stpmax = stepmax*(sqrt(total(p^2)) > float(nparams))
  
  for i=0, itmax-1 do begin
    iter=i
    
    lnsrch2, p, fp, g, xi, pnew, fmin, stpmax, check, func
    
    fp = fmin
    xi = pnew - p
    p = pnew
    
    temp = abs(xi)/(abs(p) > 1)
    test = max(temp)
    if test lt tolx then return
    
    dg = g
    g = call_function(dfunc, p)
    
    den = fmin > 1
    temp = abs(g)*(abs(p) > 1)/den
    test = max(temp)
    if test lt gtol then return
    
    dg = g-dg
    hdg = matrix_multiply(hessian, dg)
    fac = total(dg*xi)
    fae = total(dg*hdg)
    sumdg = total(dg^2.)
    sumxi = total(xi^2.)
    
    if fac gt sqrt(eps*sumdg*sumxi) then begin
      fac = 1./fac
      fad = 1./fae
      dg = fac*xi-fad*hdg
      hessian += fac*matrix_multiply(xi, xi) - fad*matrix_multiply(hdg, hdg) + fae*matrix_multiply(dg, dg)
    endif
    xi = -1*matrix_multiply(hessian, g)
    
  endfor
;print, 'max iterations reached, convergence not complete'
  
end