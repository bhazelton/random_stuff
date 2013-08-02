pro fhd_sim, version = version

except=!except
!except=0 
heap_gc

IF N_Elements(version) EQ 0 THEN version=0

IF N_Elements(recalculate_all) EQ 0 THEN recalculate_all=1
IF N_Elements(cleanup) EQ 0 THEN cleanup=0
IF N_Elements(ps_export) EQ 0 THEN ps_export=1
IF N_Elements(version) EQ 0 THEN version=0
IF N_Elements(channel) EQ 0 THEN channel=145
image_filter_fn='filter_uv_radial' ;applied ONLY to output images

data_directory=rootdir('mwa')+filepath('',root='DATA',subdir=['X16','EOR1',Strn(Floor(channel))])
vis_file_list=file_search(data_directory,'*_cal.uvfits',count=n_files)
fhd_file_list=fhd_path_setup(vis_file_list,version=version)
healpix_path=fhd_path_setup(output_dir=data_directory,subdir='Healpix',output_filename='Combined_obs',version=version)
catalog_file_path=filepath('MRC_full_radio_catalog.fits',root=rootdir('FHD'),subdir='catalog_data')


dimension=1024.
max_sources=10000.
pad_uv_image=2.
complex_beam=1
rephase_to_zenith=1
FoV=90.
no_ps=1 ;don't save postscript copy of images
min_baseline=12.

n_files=N_Elements(vis_file_list)

;Set which files to restore or recalculate (if the file is not found and needed, it will be recalculated
IF N_Elements(double_precison_beam) EQ 0 THEN double_precison_beam=0
IF N_Elements(beam_recalculate) EQ 0 THEN beam_recalculate=recalculate_all
IF N_Elements(healpix_recalculate) EQ 0 THEN healpix_recalculate=0
IF N_Elements(flag) EQ 0 THEN flag=0
IF N_Elements(grid) EQ 0 THEN grid=recalculate_all

;Set up gridding and deconvolution parameters
IF N_Elements(n_pol) EQ 0 THEN n_pol=2


IF N_Elements(beam_recalculate) EQ 0 THEN beam_recalculate=1
IF N_Elements(grid_recalculate) EQ 0 THEN grid_recalculate=1
IF N_Elements(healpix_recalculate) EQ 0 THEN healpix_recalculate=0

ext='.uvfits'
fhd_dir=file_dirname(file_path_fhd)
basename=file_basename(file_path_fhd)
header_filepath=file_path_fhd+'_header.sav'
flags_filepath=file_path_fhd+'_flags.sav'
;vis_filepath=file_path_fhd+'_vis.sav'
obs_filepath=file_path_fhd+'_obs.sav'
params_filepath=file_path_fhd+'_params.sav'
hdr_filepath=file_path_fhd+'_hdr.sav'
fhd_filepath=file_path_fhd+'_fhd.sav'
autocorr_filepath=file_path_fhd+'_autos.sav'
cal_filepath=file_path_fhd+'_cal.sav'

pol_names=['xx','yy','xy','yx','I','Q','U','V']

IF Keyword_Set(n_pol) THEN n_pol1=n_pol ELSE n_pol1=1

data_flag=file_test(hdr_filepath) AND file_test(flags_filepath) AND file_test(obs_filepath) AND file_test(params_filepath)

IF Keyword_Set(beam_recalculate) OR Keyword_Set(grid_recalculate) OR $
   ~data_flag THEN data_flag=1 ELSE data_flag=0

IF Keyword_Set(force_data) THEN data_flag=1
IF Keyword_Set(force_no_data) THEN data_flag=0


IF Keyword_Set(data_flag) THEN BEGIN
    ;info_struct=mrdfits(filepath(filename+ext,root_dir=rootdir('mwa'),subdir=data_directory),2,info_header,/silent)
    IF file_test(file_path_vis) EQ 0 THEN BEGIN
        message,"File: "+file_path_vis+" not found!"
    ENDIF
    data_struct=mrdfits(file_path_vis,0,data_header0,/silent)
    ;; testing for export type. If multibeam, then read original header
    casa_type = strlowcase(strtrim(sxpar(data_header0, 'OBJECT'), 2))
    if casa_type ne 'multi' then use_calheader = 1 else use_calheader = 0
    if use_calheader eq 0 then begin
      IF file_test(header_filepath) EQ 0 THEN uvfits_header_casafix,file_path_vis,file_path_fhd=file_path_fhd
      RESTORE,header_filepath
      hdr=vis_header_extract(data_header,header2=data_header0, params = data_struct.params)
    endif else begin
      hdr=vis_header_extract(data_header0, params = data_struct.params)
    
    endelse
    
    params=vis_param_extract(data_struct.params,hdr)
    obs=vis_struct_init_obs(hdr,params,n_pol=n_pol,_Extra=extra)
    kbinsize=obs.kpix
    degpix=obs.degpix
    dimension=obs.dimension
    
    IF Keyword_Set(rephase_to_zenith) THEN BEGIN
        obs0=obs
        hdr0=hdr
        phasera=obs0.obsra
        phasedec=obs0.obsdec        
        hdr.obsra=obs0.zenra
        hdr.obsdec=obs0.zendec
        obs.obsx=obs0.zenx
        obs.obsy=obs0.zeny
        obs=vis_struct_init_obs(hdr,params,n_pol=n_pol,phasera=phasera,phasedec=phasedec,_Extra=extra)
;        obs=vis_struct_init_obs(hdr,params,n_pol=n_pol,obsx=obs0.zenx,obsy=obs0.zeny,$
;            phasera=phasera,phasedec=phasedec,_Extra=extra)
;        obs=vis_struct_init_obs(hdr,params,n_pol=n_pol,obsx=obs.zenx,obsy=obs.zeny,_Extra=extra)
        phase_shift=1.
;        obs1=vis_struct_init_obs(hdr,params,n_pol=n_pol,phasera=phasera,phasedec=phasedec,_Extra=extra)
;        kx_arr=params.uu/kbinsize
;        ky_arr=params.vv/kbinsize
;        xcen=Float(obs.freq#kx_arr)
;        ycen=Float(obs.freq#ky_arr)
;        phase_shift=Exp((2.*!Pi*Complex(0,1)/dimension)*((obs.obsx-obs1.obsx)*xcen+(obs.obsy-obs1.obsy)*ycen))
;        obs=vis_struct_init_obs(hdr,params,n_pol=n_pol,_Extra=extra)
    ENDIF ELSE phase_shift=1.
    
    IF Tag_exist(obs,'freq') THEN freq_arr=obs.freq ELSE freq_arr=(*obs.baseline_info).freq
    
    IF Keyword_Set(freq_start) THEN bw_start=(freq_start*1E6)>Min(freq_arr) ELSE bw_start=Min(freq_arr)
    IF Keyword_Set(freq_end) THEN bw_end=(freq_end*1E6)<Max(freq_arr) ELSE bw_end=Max(freq_arr)
    bandwidth=Round((bw_end-bw_start)/1E5)/10.
    fov=dimension*degpix
    k_span=kbinsize*dimension
    
    
    print,String(format='("Image size used: ",A," pixels")',Strn(dimension))
    print,String(format='("Image resolution used: ",A," degrees/pixel")',Strn(degpix))
    print,String(format='("Approx. beam area: ",A," pixels")',Strn((!RaDeg/(obs.MAX_BASELINE/obs.KPIX)/obs.degpix)))
    print,String(format='("Field of view used: ",A," degrees")',Strn(fov))
    print,String(format='("Frequency range: ",A,"-",A," MHz")',Strn(Round((bw_start)/1E5)/10.),Strn(Round((bw_end)/1E5)/10.))
    print,String(format='("UV resolution used: ",A," wavelengths/pixel")',Strn(kbinsize))
    print,String(format='("UV image size used: ",A," wavelengths")',Strn(k_span))
    print,String(format='("Min baseline: ",A," wavelengths")',Strn(obs.min_baseline))
    print,String(format='("Max baseline: ",A," wavelengths")',Strn(obs.max_baseline))
    print,String(format='("Observation coordinates: ",A," ",A,A)',$
        Strn(obs.obsra,length=6),(obs.obsdec GE 0) ? '+':'-',Strn(Abs(obs.obsdec),length=5))
    print,String(format='("Zenith coordinates: ",A," ",A,A)',$
        Strn(obs.zenra,length=6),(obs.zendec GE 0) ? '+':'-',Strn(Abs(obs.zendec),length=5))
    IF (obs.phasera NE obs.obsra) OR (obs.phasedec NE obs.obsdec) THEN $
        print,String(format='("Image phased to coordinates: ",A," ",A,A)',$
            Strn(obs.phasera,length=6),(obs.phasedec GE 0) ? '+':'-',Strn(Abs(obs.phasedec),length=5))
    
    IF Tag_exist(obs,'alpha') THEN alpha=obs.alpha ELSE alpha=0.
    print,String(format='("Spectral index fit: ",A)',Strn(alpha))
    pol_dim=hdr.pol_dim
    freq_dim=hdr.freq_dim
    real_index=hdr.real_index
    imaginary_index=hdr.imaginary_index
    flag_index=hdr.flag_index
    n_pol=obs.n_pol
    n_freq=obs.n_freq
    
    data_array=data_struct.array[*,0:n_pol-1,*]
    data_struct=0. ;free memory
    
    ;Read in or construct a new beam model. Also sets up the structure PSF
    print,'Calculating beam model'
    psf=beam_setup(obs,file_path_fhd,restore_last=(Keyword_Set(beam_recalculate) ? 0:1),silent=silent,timing=t_beam,_Extra=extra)
    IF Keyword_Set(t_beam) THEN print,'Beam modeling time: ',t_beam
    
    vis_arr=Ptrarr(n_pol,/allocate)
    flag_arr=Ptrarr(n_pol,/allocate)
    FOR pol_i=0,n_pol-1 DO BEGIN
        *vis_arr[pol_i]=Complex(reform(data_array[real_index,pol_i,*,*]),$
            Reform(data_array[imaginary_index,pol_i,*,*]))*phase_shift
        *flag_arr[pol_i]=reform(data_array[flag_index,pol_i,*,*])
    ENDFOR
    ;free memory
    data_array=0 
    flag_arr0=0
 endif

;; make model uv image
CASE 1 OF
    Keyword_Set(complex_beam) AND Keyword_Set(double_precison_beam): model_uv=Dcomplexarr(dimension,dimension)
    Keyword_Set(double_precison_beam): model_uv=Dblarr(dimension,dimension)
    Keyword_Set(complex_beam): model_uv=Complexarr(dimension,dimension)
    ELSE: model_uv=Fltarr(dimension,dimension)
ENDCASE

;; just one non-zero mode for now
pix_inds = [550, 620]
model_uv[pix_inds] = 1.

model_vis = visibility_degrid(model_uv, flag_arr,obs,psf,params,$
    timing=timing,polarization=polarization,silent=silent,$
    complex=complex_beam,double=double_precison_beam)

stop

!except=except

END
