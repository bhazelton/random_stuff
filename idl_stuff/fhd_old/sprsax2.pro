PRO sprsax2,A,X,B,X2,B2,double=double,transpose=transpose,mask=mask
;slight modification to sprsax to allow much larger arrays
;also modified to more efficiently use sparse vectors if mask is supplied 
; <<< MASK ONLY WORKS WITH TRANSPOSE!>>>
;set keyword "transpose" to use the transpose of A instead, without having to do extra calculations

sa=A.sa
ija=A.ija-1

n=A.ija[0]-2. ;DO NOT include an extra -1 here. It MUST be ija[0]-2.
IF Keyword_Set(double) THEN B=dblarr(n) ELSE B=Fltarr(n)

IF N_Params() LE 3 THEN BEGIN
    FOR i=0.,n-1 DO BEGIN
        B[i]=sa[i]*X[i]
        i2=ija[i+1]-1
        i1=ija[i]
        IF i2 LT i1 THEN CONTINUE
        sa_sub=sa[i1:i2]
        ija_sub=ija[i1:i2]
        B[i]+=Total(sa_sub*X[ija_sub],/double)
    ENDFOR
ENDIF ELSE BEGIN
    B2=B
    FOR i=0.,n-1 DO BEGIN
        B[i]=sa[i]*X[i]
        B2[i]=sa[i]*X2[i]
        i2=ija[i+1]-1
        i1=ija[i]
        IF i2 LT i1 THEN CONTINUE
        sa_sub=sa[i1:i2]
        ija_sub=ija[i1:i2]
        B[i]+=Total(sa_sub*X[ija_sub],/double)
        B2[i]+=Total(sa_sub*X2[ija_sub],/double)
    ENDFOR
ENDELSE
END

;;FOR i=0.,n-1 DO BEGIN
;;    B[i]=A.sa[i]*x[i]
;;    FOR k=A.ija[i]-1,A.ija[i+1]-1-1 DO B[i]+=A.sa[k]*x[A.ija[k]-1]
;;ENDFOR
;
;IF N_Elements(mask) EQ N_Elements(x) THEN BEGIN
;    i_use=where(mask,n_use)
;    IF Keyword_Set(transpose) THEN BEGIN
;        FOR i1=0.,n_use-1 DO BEGIN
;            i=i_use[i1]
;            B[i]=sa[i]*x[i]
;            FOR k=ija[i],ija[i+1]-1 DO B[ija[k]]+=sa[k]*x[i]
;        ENDFOR
;    ENDIF ELSE BEGIN
;        FOR i=0.,n-1 DO BEGIN
;            B[i]=sa[i]*x[i]
;            FOR k=ija[i],ija[i+1]-1 DO B[i]+=sa[k]*x[ija[k]]
;        ENDFOR
;    ENDELSE
;ENDIF ELSE BEGIN
;    IF Keyword_Set(transpose) THEN BEGIN
;        FOR i=0.,n-1 DO BEGIN
;            B[i]=sa[i]*x[i]
;            FOR k=ija[i],ija[i+1]-1 DO B[ija[k]]+=sa[k]*x[i]
;        ENDFOR
;    ENDIF ELSE BEGIN
;        FOR i=0.,n-1 DO BEGIN
;            B[i]=sa[i]*x[i]
;            FOR k=ija[i],ija[i+1]-1 DO B[i]+=sa[k]*x[ija[k]]
;        ENDFOR
;    ENDELSE
;ENDELSE
;IF Keyword_Set(transpose) THEN BEGIN
;    FOR i=0.,n-1 DO BEGIN
;        B[i]=sa[i]*x[i]
;        FOR k=ija[i],ija[i+1]-1 DO B[ija[k]]+=sa[k]*x[i]
;    ENDFOR
;ENDIF ELSE BEGIN
;    FOR i=0.,n-1 DO BEGIN
;        B[i]=sa[i]*x[i]
;        FOR k=ija[i],ija[i+1]-1 DO B[i]+=sa[k]*x[ija[k]]
;    ENDFOR
;ENDELSE
;RETURN,B
;END