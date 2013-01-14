FUNCTION visibility_transfer_convert,xfer_fn,restore_last=restore_last,timing=timing,file_name_base=file_name_base,$
    psf_dim=psf_dim,threshold=threshold
t0=Systime(1)
IF N_Elements(restore_last) EQ 0 THEN restore_last=1

COMMON obs_params,dimension,elements,degpix,kbinsize,azimuth_center,elevation_center,rotation
;COMMON psf_params,psf_base,psf_dim
COMMON visibility_params,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,freq_bin_i,file_path
IF Keyword_Set(restore_last) THEN BEGIN
    restore,file_path+file_name_base+'.sav'
;    xfer_fn=read_spr(file_path+file_name_base+'.as')
    timing=Systime(1)-t0
    RETURN,xfer_fn
ENDIF
IF N_Elements(threshold) EQ 0 THEN threshold=0.

;convert pointer array transfer function to a sparse matrix

;Result = SPRSIN(Columns, Rows, Values, N [, /DOUBLE] [, THRESHOLD=value])
;Result = LINBCG( A, B, X [, /DOUBLE] [, ITOL={4 | 5 | 6 | 7}] [, TOL=value] [, ITER=variable] [, ITMAX=value] )
;A: A row-indexed sparse array created by the SPRSIN function.
;B: An n-element vector containing the right-hand side of the linear system Ax=b.
;X: An n-element vector containing the initial solution of the linear system.

psf_dim2=2*psf_dim
psf_n=psf_dim2^2.
psf_i=Indgen(psf_n)
sub_xv=meshgrid(psf_dim2,1)-psf_dim
sub_yv=meshgrid(psf_dim2,2)-psf_dim

n=dimension*elements

print, 'Starting diagonal terms'
t1 = systime(1)

n_arr=fltarr(dimension,elements)
diag_vals=fltarr(dimension,elements)
FOR xi=1,dimension-2 DO BEGIN
    FOR yi=1,elements-2 DO BEGIN
        temp_arr=*xfer_fn[xi,yi]
        diag_vals[xi,yi]=temp_arr[psf_dim,psf_dim]
        temp_arr[psf_dim,psf_dim]=0
        i1=where(Abs(temp_arr) GT threshold,n1)
;        IF temp_arr[psf_dim,psf_dim] NE 0 THEN n1-=1
        n_arr[xi,yi]=n1
    ENDFOR
ENDFOR

t2 = systime(1)
print, 'Finished diagonal terms, time = ' + string(t2-t1)

;i_use=where(n_arr,n_use)
;offset_i=Total(n_arr,/cumulative)

;length of the vectors containing the sparse array indices and values
n_sparse=dimension*elements+Total(n_arr,/double)+1.
ija=Lon64arr(n_sparse) ;sparse indices, following the convention of sprsin (NumRec in C, ch 2.7)
sa=Fltarr(n_sparse)
sa[0:n-1]=reform(diag_vals,n)

ki=Long64(n)+1
ija[0]=n+2 ;index to first off-diagonal term

print, 'Starting off diagonal terms'
t3 = systime(1)

FOR i=0.,n-1 DO BEGIN
;FOR ii=0.,n_use-1 DO BEGIN
    IF n_arr[i] GT 0 THEN BEGIN
        xi=Float(i mod dimension)
        yi=Float(Floor(i/dimension))
        xfer_fn_sub=*xfer_fn[xi,yi]
        *xfer_fn[xi,yi]=0 ;free memory as soon as it's read
        xfer_fn_sub[psf_dim,psf_dim]=0 ;diagonal terms are treated previously
        j_use=where(Abs(xfer_fn_sub) GT threshold,n_use)
        
;        wt=Total(xfer_fn_sub)
;        IF wt EQ 0 THEN CONTINUE ;these SHOULD all be nonzero already, but check anyway
        IF n_use GT 0 THEN BEGIN ;use 'IF' instead of a BREAK, though that would work too.
            xii_arr=sub_xv[j_use]+xi
            yii_arr=sub_yv[j_use]+yi    
;            ind2=xii_arr+yii_arr*dimension
;            vals=xfer_fn_sub[j_use]
            sa[ki:ki+n_use-1]=xfer_fn_sub[j_use]
            ija[ki:ki+n_use-1]=xii_arr+yii_arr*dimension
            ki+=n_use
        ENDIF
;        FOR j=0,n_use-1 DO BEGIN
;            ki+=1
;            sa[ki]=vals[j]
;            ija[ki]=ind2[j]
;        ENDFOR
    ENDIF
    ija[i+1]=ki+1 ;store index to next row
ENDFOR

t4 = systime(1)
print, 'Finished off diagonal terms, time = ' + string(t4-t3)


;;heap_free,xfer_fn
 Ptr_free,xfer_fn
 heap_gc
;; xfer_fn={xfer_fn,ija:Temporary(ija),sa:Temporary(sa)}

t5 = systime(1)
print, 'Creating structure. Pointer free time = ' + string(t5-t4)

xfer_fn={ija:Temporary(ija),sa:Temporary(sa)}
t6 = systime(1)
print, 'Created structure, time = ' + string(t6-t5)

timing=Systime(1)-t0
save,xfer_fn,filename=file_path+file_name_base+'.sav'

t5 = systime(1)
print, 'Finished xfer function conversion, total time = ' + string(t5-t0)

RETURN,xfer_fn
END

;Result = LINBCG( A, B, X [, /DOUBLE] [, ITOL={4 | 5 | 6 | 7}] [, TOL=value] [, ITER=variable] [, ITMAX=value] )
;A: A row-indexed sparse array created by the SPRSIN function.
;B: An n-element vector containing the right-hand side of the linear system Ax=b.
;X: An n-element vector containing the initial solution of the linear system.

;image_uv_real=Reform(Real_Part(image_uv),n)
;image_uv_comp=Reform(Imaginary(image_uv),n)
;image_uv_init=Dcomplexarr(dimension,elements)
;image_uv_init[where(weights)]=image_uv[where(weights)]/weights[where(weights)]
;image_uv_real_init=Reform(Real_Part(image_uv_init),n)
;image_uv_comp_init=Reform(Imaginary(image_uv_init),n)
;;image_uv_real_init=image_uv_real & image_uv_real_init[where(weights)]/=weights[where(weights)]
;;image_uv_comp_init=image_uv_comp & image_uv_comp_init[where(weights)]/=weights[where(weights)]
;
;;included to illustrate use
;;xfer_fn2_real_Transpose = SPRSTP(xfer_fn2_real)
;;xfer_fn2_comp_Transpose = SPRSTP(xfer_fn2_comp)
;;dirty_image_uv = SPRSAX(xfer_fn2_real,true_image_uv)
;
;image_uv_real_estimate=Linbcg(xfer_fn2_real,image_uv_real,image_uv_real_init,/double)
;image_uv_comp_estimate=Linbcg(xfer_fn2_comp,image_uv_comp,image_uv_comp_init,/double)
;image_uv_estimate=complex(image_uv_real_estimate,image_uv_comp_estimate)
;image_uv_estimate=Reform(image_uv_estimate,dimension,elements)
;timing=Systime(1)-t0
;END
