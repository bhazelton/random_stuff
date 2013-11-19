PRO testcal,cleanup=cleanup,ps_export=ps_export,recalculate_all=recalculate_all,export_images=export_images,version=version,$
    beam_recalculate=beam_recalculate,healpix_recalculate=healpix_recalculate,mapfn_recalculate=mapfn_recalculate,$
    grid=grid,deconvolve=deconvolve,channel=channel,_Extra=extra
except=!except
!except=0 
heap_gc

calibrate_visibilities=1
IF N_Elements(recalculate_all) EQ 0 THEN recalculate_all=0
IF N_Elements(export_images) EQ 0 THEN export_images=1
IF N_Elements(cleanup) EQ 0 THEN cleanup=1
IF N_Elements(ps_export) EQ 0 THEN ps_export=0
IF N_Elements(deconvolve) EQ 0 THEN deconvolve=0
IF N_Elements(mapfn_recalculate) THEN mapfn_recalculate=0
IF N_Elements(healpix_recalculate) EQ 0 THEN healpix_recalculate=0
IF N_Elements(beam_recalculate) EQ 0 THEN beam_recalculate=0
IF N_Elements(vis_baseline_hist) EQ 0 THEN vis_baseline_hist=1
IF N_Elements(version) EQ 0 THEN version='bjh_1'
image_filter_fn='' ;applied ONLY to output images

;;IF StrLowCase(!version.os_family) EQ 'unix' THEN data_directory=rootdir('mwa')+filepath('',root='DATA',subdir=['128T','test']) $
;;    ELSE data_directory=rootdir('mwa')+filepath('',root='DATA3',subdir=['128T','testcal'])
;;data_directory = base_path('data') + 'fhd_ps_data/EOR1_P00_145_20110926193959/'
data_directory = base_path('data') + 'fhd_ps_data/1061316296/'
vis_file_list=file_search(data_directory,'*.uvfits',count=n_files)
fhd_file_list=fhd_path_setup(vis_file_list,version=version,_Extra=extra)
healpix_path=fhd_path_setup(output_dir=data_directory,subdir='Healpix',output_filename='Combined_obs',version=version,_Extra=extra)
;catalog_file_path=filepath('MRC_full_radio_catalog.fits',root=rootdir('FHD'),subdir='catalog_data')
calibration_catalog_file_path=filepath('mwa_commissioning_source_list.sav',root=rootdir('FHD'),subdir='catalog_data')

;noise_calibrate=0
;align=0
dimension=2048.
max_sources=0
pad_uv_image=2.
FoV=100.
no_ps=1 ;don't save postscript copy of images
psf_dim=10
min_baseline=1.
min_cal_baseline=50.
flag_nsigma=20.
ring_radius=6.*pad_uv_image
nfreq_avg=16
max_calibration_sources=300.
weights_grid=0
psf_resolution=8.

calib_freq_func = 1
max_cal_iter=30L
initial_calibration = data_directory + 'fhd_bjh_1/1061316296_cal.sav'
;calibration_visibilities_subtract=1

general_obs,cleanup=cleanup,ps_export=ps_export,recalculate_all=recalculate_all,export_images=export_images,version=version,$
    beam_recalculate=beam_recalculate,healpix_recalculate=healpix_recalculate,mapfn_recalculate=mapfn_recalculate,$
    grid=grid,deconvolve=deconvolve,image_filter_fn=image_filter_fn,data_directory=data_directory,output_directory=output_directory,$
    vis_file_list=vis_file_list,fhd_file_list=fhd_file_list,healpix_path=healpix_path,catalog_file_path=catalog_file_path,$
    dimension=dimension,max_sources=max_sources,pad_uv_image=pad_uv_image,psf_dim=psf_dim,$
    FoV=FoV,no_ps=no_ps,min_baseline=min_baseline,flag_visibilities=flag_visibilities,flag_calibration=flag_calibration,$
    calibrate_visibilities=calibrate_visibilities,calibration_catalog_file_path=calibration_catalog_file_path,$
    ring_radius=ring_radius,flag_nsigma=flag_nsigma,nfreq_avg=nfreq_avg,max_calibration_sources=max_calibration_sources,$
    calibration_visibilities_subtract=calibration_visibilities_subtract,weights_grid=weights_grid,/mark_zenith,$
    vis_baseline_hist=vis_baseline_hist,psf_resolution=psf_resolution,no_rephase=no_rephase,show_obsname=show_obsname,$
    silent=silent,smooth_width=smooth_width,gain_factor=gain_factor,combine_obs=combine_obs,no_fits=no_fits,$
    min_cal_baseline=min_cal_baseline,calib_freq_func=calib_freq_func,$
    max_cal_iter=max_cal_iter, initial_calibration=initial_calibration, _Extra=extra

!except=except
END
