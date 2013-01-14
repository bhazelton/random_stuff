FUNCTION weight_invert,weights,threshold
result=fltarr(size(weights,/dimension))

IF Keyword_Set(threshold) THEN i_use=where(weights GE threshold,n_use) $
    ELSE i_use=where(weights,n_use)
IF n_use GT 0 THEN result[i_use]=1./weights[i_use]

RETURN,result
END