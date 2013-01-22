PRO visibility_transfer_generate,xfer_fn,restore_last=restore_last,timing=timing,test=test,new=new,$
    polarization=polarization,flag_arr=flag_arr,file_identifier=file_identifier, freq_inds_to_use = freq_inds_to_use
t0=Systime(1)

COMMON obs_params,dimension,elements,degpix,kbinsize,azimuth_center,elevation_center,rotation
COMMON psf_params,psf_base,psf_dim,psf_resolution,psf_residuals_i,psf_residuals_val,psf_residuals_n,psf_xvals,psf_yvals
COMMON visibility_params,uu_arr,vv_arr,frequency_array,date_arr,baseline_arr,bin_offset,freq_bin_i,file_path
COMMON switches,source_switch,psf_switch
;IF N_Elements(polarization) EQ 0 THEN polarization=0 ;could allow default, but disallow to make sure it is specified 
file_name_base='_xferfn_'+file_identifier

IF Keyword_Set(restore_last) THEN BEGIN
    IF Keyword_Set(new) THEN restore,file_path+file_name_base+'.sav' ELSE restore,file_path+file_name_base+'_old.sav'
    timing=Systime(1)-t0
    RETURN
ENDIF

source_switch=1 ;only works with this option selected
kx_arr=uu_arr/kbinsize
ky_arr=vv_arr/kbinsize
;;nbaselines=bin_offset[1] -- changed to support n_samples=1
n_samples=N_Elements(bin_offset)
nbaselines=n_elements(baseline_arr)/n_samples

if n_elements(freq_inds_to_use) eq 0 then freq_inds_to_use = indgen(n_elements(frequency_array))
n_frequencies=N_Elements(frequency_array[freq_inds_to_use])
n_freq_bin=N_Elements(freq_bin_i[freq_inds_to_use])
psf_dim2=2*psf_dim

xfer_fn=Ptrarr(dimension,elements,/allocate)
FOR i=0,dimension-1 DO FOR j=0,elements-1 DO *xfer_fn[i,j]=fltarr(psf_dim2,psf_dim2)

time_check_interval=Ceil(nbaselines/20.)
t1=Systime(1)
t2=0
t3=0

xmin_old=-1
ymin_old=-1
box_arr=Ptrarr(psf_dim,psf_dim,/allocate)   
FOR bi=0.,nbaselines-1 DO BEGIN

   IF Keyword_Set(flag_arr) THEN IF total(flag_arr) LE 0 THEN CONTINUE

    IF Keyword_Set(time_check_interval) AND (bi mod time_check_interval EQ 0) THEN BEGIN
        t1b=Systime(1)-t1
        print,Strcompress(String(format='("Time elapsed:",I," Estimated time remaining:",I)',t1b,(t1b/bi)*(nbaselines-bi)))
    ENDIF 
    
    IF kx_arr[bi] EQ 0 and ky_arr[bi] eq 0 THEN CONTINUE ;skip the auto-correlations

    FOR fi=0.,n_frequencies-1 DO BEGIN
        freq_i=freq_bin_i[freq_inds_to_use[fi]]
        freq=frequency_array[freq_inds_to_use[fi]]
        FOR ti=0.,n_samples-1 DO BEGIN
;            bi_use=bi+bin_offset[ti]                        
;            IF Keyword_Set(flag_arr) THEN IF flag_arr[fi,bi_use] LE 0 THEN CONTINUE
;            xc=kx_arr[bi_use]*freq
;            yc=ky_arr[bi_use]*freq
;            t2a=Systime(1)
;            psf_use=visibility_psf(freq_i,bi,polarization=polarization,xcenter=xc,ycenter=yc,xmin=xmin,ymin=ymin,xcen0=xcen0,ycen0=ycen0,/short)
;            psf_i=where(psf_use,n_use)
;            psf_vals=psf_use[psf_i]
;            psf_x=psf_i mod psf_dim
;            psf_y=Floor(psf_i/psf_dim)
;            t2+=Systime(1)-t2a
;            t3a=Systime(1)
;            FOR ni=0,n_use-1 DO BEGIN
;                i=psf_x[ni]
;                j=psf_y[ni]
;                (*xfer_fn[xmin+i,ymin+j])[psf_dim-i:2*psf_dim-i-1,psf_dim-j:2*psf_dim-j-1]+=psf_vals[ni]*psf_use
;            ENDFOR
;            t3+=Systime(1)-t3a

            bi_use=bi+bin_offset[ti]                        
            IF Keyword_Set(flag_arr) THEN IF flag_arr[fi,bi_use] LE 0 THEN CONTINUE
            xc=kx_arr[bi_use]*freq
            yc=ky_arr[bi_use]*freq
            t2a=Systime(1)
            psf_use=visibility_psf(freq_i,bi,polarization=polarization,xcenter=xc,ycenter=yc,xmin=xmin,ymin=ymin,xcen0=xcen0,ycen0=ycen0,/short)

            IF (xmin_old NE xmin) OR (ymin_old NE ymin) THEN BEGIN
                FOR i=0,psf_dim-1 DO FOR j=0,psf_dim-1 DO BEGIN
                    IF xmin_old GE 0 THEN (*xfer_fn[xmin_old+i,ymin_old+j])[psf_dim-i:2*psf_dim-i-1,psf_dim-j:2*psf_dim-j-1]+=*box_arr[i,j]
                    *box_arr[i,j]=fltarr(psf_dim,psf_dim)
                ENDFOR
                xmin_old=xmin
                ymin_old=ymin
            ENDIF
            psf_i=where(psf_use,n_use)
            if n_use eq 0 then stop else psf_vals=psf_use[psf_i]
            psf_x=psf_i mod psf_dim
            psf_y=Floor(psf_i/psf_dim)
            t2+=Systime(1)-t2a
            t3a=Systime(1)
            FOR ni=0,n_use-1 DO BEGIN
                i=psf_x[ni]
                j=psf_y[ni]
;                (*xfer_fn[xmin+i,ymin+j])[psf_dim-i:2*psf_dim-i-1,psf_dim-j:2*psf_dim-j-1]+=psf_vals[ni]*psf_use
                *box_arr[i,j]+=psf_vals[ni]*psf_use
            ENDFOR
            t3+=Systime(1)-t3a
;            IF Keyword_Set(test) THEN test1=fltarr(psf_dim*psf_dim2,psf_dim*psf_dim2)
        ENDFOR
    ENDFOR
ENDFOR
;Process LAST box
IF xmin_old GE 0 then FOR i=0,psf_dim-1 DO FOR j=0,psf_dim-1 DO $
   (*xfer_fn[xmin_old+i,ymin_old+j])[psf_dim-i:2*psf_dim-i-1,psf_dim-j:2*psf_dim-j-1]+=*box_arr[i,j]
                    
t1=Systime(1)-t1

IF Keyword_Set(new) THEN xfer_fn=visibility_transfer_convert(xfer_fn,restore_last=0,file_name_base=file_name_base,psf_dim=psf_dim) $
else save,xfer_fn,file=file_path+file_name_base+'_old.sav'
timing=Systime(1)-t0

END
