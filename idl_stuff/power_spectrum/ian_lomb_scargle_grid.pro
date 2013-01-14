FUNCTION ian_lomb_scargle_grid,image,degpix=degpix,mask=mask,timing=timing,inverse=inverse,$
    delta1=delta1,delta2=delta2
;wrapper routine to make lomb-scargle algorithm mimic IDL FFT for REGULARLY GRIDDED DATA (holes are OK)
;set keyword inverse to transform back to image space from Fourier space. NOT YET WORKING!

dimension=(size(image,/dimension))[0]
elements=(size(image,/dimension))[1]

IF N_Elements(degpix) EQ 0 THEN degpix=1d/dimension
sterad=(degpix*!DtoR)^2.
n_k=dimension<elements
kbinsize=1d/(degpix*n_k)

IF Keyword_Set(mask) THEN $
    IF N_Elements(mask) EQ N_Elements(image) THEN i_use=where(mask,n_use,comp=mask_i) ELSE i_use=where(image,n_use,comp=mask_i) $
ELSE i_use=findgen(dimension*elements)

delta1=dblarr(n_k)
delta2=dblarr(n_k)

IF Keyword_Set(inverse) THEN BEGIN
    xvals=dindgen(dimension)*degpix
    yvals=dindgen(elements)*degpix
    
    kvals=dist(n_k,n_k)*kbinsize
    kxvals=kvals[*,0]#Replicate(1,elements)
    kyvals=kvals[0,*]##Replicate(1,dimension)
    
    kxvals_use=kxvals[i_use]
    kyvals_use=kyvals[i_use]
    data_use=image[i_use]
    
    IF Keyword_Set(mask) THEN BEGIN
        IF n_use GE 0.5*dimension*elements THEN i_use2=mask_i ELSE i_use2=i_use 
        FOR k=1,n_k-1 DO BEGIN
            delta1[k]=atan(total(sin(4d*!Dpi*xvals[k]*kxvals[i_use2])), total(cos(4d*!Dpi*xvals[k]*kxvals[i_use2]))) / 2d
            delta2[k]=atan(total(sin(4d*!Dpi*yvals[k]*kyvals[i_use2])), total(cos(4d*!Dpi*yvals[k]*kyvals[i_use2]))) / 2d
        ENDFOR
    ENDIF 
    
    image_ft=ian_lomb_scargle_test(kxvals_use,kyvals_use,data_use,xvals,yvals,delta1=delta1,delta2=delta2,timing=timing)
ENDIF ELSE BEGIN
;    kvals=dist(n_k,n_k)*kbinsize
;    kxvals=Reform(kvals[*,0])
;    kyvals=Reform(kvals[0,*])
    kxvals=dindgen(n_k)*kbinsize
    kyvals=dindgen(n_k)*kbinsize
    
;;    xvals=meshgrid(dimension,elements,1)*degpix
;;    yvals=meshgrid(dimension,elements,2)*degpix
    xvals = rebin(dindgen(dimension)*degpix, dimension, elements)
    yvals = rebin(reform(dindgen(elements)*degpix, 1, elements), dimension, elements)

    xvals_use=xvals[i_use]
    yvals_use=yvals[i_use]
    data_use=image[i_use]
    
    IF Keyword_Set(mask) THEN BEGIN
        IF n_use GE 0.5*dimension*elements THEN i_use2=mask_i ELSE i_use2=i_use 
        FOR k=1,n_k-1 DO BEGIN
            delta1[k]=atan(total(sin(4d*!Dpi*kxvals[k]*xvals[i_use2])), total(cos(4d*!Dpi*kxvals[k]*xvals[i_use2]))) / 2d
            delta2[k]=atan(total(sin(4d*!Dpi*kyvals[k]*yvals[i_use2])), total(cos(4d*!Dpi*kyvals[k]*yvals[i_use2]))) / 2d
        ENDFOR
    ENDIF
    
    image_ft=ian_lomb_scargle_test(xvals_use,yvals_use,data_use,kxvals,kyvals,delta1=delta1,delta2=delta2,timing=timing)
ENDELSE

RETURN,image_ft
END
