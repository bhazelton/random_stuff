function ian_lomb_scargle_test, locations1, locations2, data, fourier_locs1, fourier_locs2, timing = timing, $
    delta1 = delta1, delta2 = delta2
;; find (2D) lomb-scargle periodogram for unevenly spaced data
;; locations1 & 2 are x/y values of points
;; fourier_locs1 & 2 are kx/ky values to test at
;; modified version returns Fourier transform instead of power

time0 = systime(1)
if n_elements(locations1) ne n_elements(data) or n_elements(locations2) ne n_elements(data) then message, $
    'locations1 & 2 must have same number of elements as data.'

n_k1 = n_elements(fourier_locs1)
n_k2 = n_elements(fourier_locs2)
n_pts = n_elements(data)
ii=Dcomplex(0,1)

;allow delta 1 and 2 to be supplied (Useful for any form of regular grid)
IF N_Elements(delta1) NE n_k1 THEN BEGIN
    delta1=dblarr(n_k1)
    FOR i=1,n_k1-1 DO delta1[i] = $
        atan(total(sin(4d*!pi*fourier_locs1[i]*locations1)), total(cos(4d*!pi*fourier_locs1[i]*locations1))) / 2d
ENDIF

IF N_Elements(delta2) NE n_k2 THEN BEGIN    
    delta2=dblarr(n_k2)
    FOR j=1,n_k2-1 DO delta2[j] = $
        atan(total(sin(4d*!pi*fourier_locs2[j]*locations2)), total(cos(4d*!pi*fourier_locs2[j]*locations2))) / 2d
ENDIF

term1a=(term1b=(term2a=(term2b=(term3a=(term3b=(term4a=(term4b=dblarr(n_k1,n_k2)))))))) 
FOR i=0.,n_pts-1 DO BEGIN
    arg1=(2d*!Dpi*fourier_locs1 * locations1[i]) - delta1
    arg2=(2d*!Dpi*fourier_locs2 * locations2[i]) - delta2
    sinx = sin(arg1)
    cosx = cos(arg1)
    siny = sin(arg2)
    cosy = cos(arg2)
    
    di=data[i]
    term1a+=di*sinx#siny 
    term2a+=di*sinx#cosy 
    term3a+=di*cosx#siny 
    term4a+=di*cosx#cosy 
    
    ;I am getting floating point errors from these terms -> supply them ahead of time for gridded data
    term1b+=(sinx^2d)#(siny^2d)
    term2b+=(sinx^2d)#(cosy^2d)
    term3b+=(cosx^2d)#(siny^2d)
    term4b+=(cosx^2d)#(cosy^2d)  
ENDFOR

;some function is limiting results to floating-point accuracy, so terms are not going to zero properly
float_err=1D/2^31D 
;float_err=Abs(!Pi-!Dpi)/(!Pi+!DPi)
izero=where(Abs(term1a)/max(abs(term1a)) LE float_err,n0) & IF n0 GT 0 THEN term1a[izero]=0
izero=where(Abs(term1b)/max(abs(term1b)) LE float_err,n0) & IF n0 GT 0 THEN term1b[izero]=0
izero=where(Abs(term2a)/max(abs(term2a)) LE float_err,n0) & IF n0 GT 0 THEN term2a[izero]=0
izero=where(Abs(term2b)/max(abs(term2b)) LE float_err,n0) & IF n0 GT 0 THEN term2b[izero]=0
izero=where(Abs(term3a)/max(abs(term3a)) LE float_err,n0) & IF n0 GT 0 THEN term3a[izero]=0
izero=where(Abs(term3b)/max(abs(term3b)) LE float_err,n0) & IF n0 GT 0 THEN term3b[izero]=0
izero=where(Abs(term4a)/max(abs(term4a)) LE float_err,n0) & IF n0 GT 0 THEN term4a[izero]=0
izero=where(Abs(term4b)/max(abs(term4b)) LE float_err,n0) & IF n0 GT 0 THEN term4b[izero]=0
term1=(term2=(term3=(term4=dblarr(n_k1,n_k2))))

normalization=dblarr(n_k1,n_k2)
; +/- signs chosen to match output of IDL FFT for regularly gridded data
i1use=where(term1b,n1use) & IF n1use GT 0 THEN BEGIN term1[i1use]=1D*term1a[i1use]/Sqrt(term1b[i1use]) & normalization[i1use]+=1d & ENDIF
i2use=where(term2b,n2use) & IF n2use GT 0 THEN BEGIN term2[i2use]=-term2a[i2use]/Sqrt(term2b[i2use]) & normalization[i2use]+=1d & ENDIF
i3use=where(term3b,n3use) & IF n3use GT 0 THEN BEGIN term3[i3use]=-term3a[i3use]/Sqrt(term3b[i3use]) & normalization[i3use]+=1d & ENDIF
i4use=where(term4b,n4use) & IF n4use GT 0 THEN BEGIN term4[i4use]=+1D*term4a[i4use]/Sqrt(term4b[i4use]) & normalization[i4use]+=1d & ENDIF

;ft=(Sqrt(term1^2.+term4^2.)+ii*Sqrt(term2^2.+term3^2.))
;normalization=sqrt(normalization)
;normalization*=n_pts
;ft/=normalization

ft=(Sqrt(term1^2.+term4^2.)+ii*Sqrt(term2^2.+term3^2.))/Sqrt(2d)   
normalization=1d/Sqrt(n_pts)
ft*=normalization

time1 = systime(1)
timing = time1-time0

return, ft
end
