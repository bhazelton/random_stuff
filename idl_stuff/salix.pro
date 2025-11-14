pro salix,obs_id=obs_id,version=version,outdir=output_directory,indir=input_folder
  except=!except
  !except=0
  heap_gc
  if ~keyword_set(output_directory) then print, 'outdir not set' 
  if ~keyword_set(version) then print, 'version not set'
  if ~keyword_set(obs_id) then print, 'obsid not set'
  ; parse command line args
  compile_opt strictarr
  args = Command_Line_Args(count=nargs)
  if nargs ne 0 then begin
    obs_id = args[0]
    output_directory = args[1]
    version = args[2]
  
    if nargs gt 3 then input_folder = args[3]
    ; if nargs gt 3 then platform = args[3] else platform = '' ;indicates if running on AWS
  endif

  cmd_args={version:version}

  case version of

    'freq_cal_only': begin
    calibration_catalog_file_path=filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')
    restrict_hpx_inds='EoR0_high_healpix_inds_3x.idlsave'
    
    ; these keywords are not in the dictionary
    max_cal_iter=1000L
    model_delay_filter=1

    ; turn off beam averaging
    beam_nfreq_avg=1

    ; recalculate all products. I'm not sure why this is set
    ; recalculate_all=1
    
    ; set flux cutoff for calibration sources
    calibration_flux_threshold=0.1

    ; turn off time averaging
    cal_time_average=0

    ; sidelobe subtraction catalogs
    calibration_subtract_sidelobe_catalog=filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')
    model_subtract_sidelobe_catalog=filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')
    subtract_sidelobe_catalog=filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')
    
    ; prevent power leakage by not gridding longer baselines
    ps_kspan=200.
    
    ; this is only a cal run so stop before gridding
    cal_stop=1
    
    ; use Ian's speedup
    use_adaptive_calibration_gain=1
    calibration_base_gain=0.5

    ; turn off bandpass, polynomial, and cable reflection fitting
    bandpass_calibrate=0 
    calibration_polyfit=0
    end

    'Wilensky_thesis_cal': begin
    calibration_catalog_file_path=filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')
    restrict_hpx_inds='EoR0_high_healpix_inds_3x.idlsave'

    ; these keywords are not in the dictionary
    max_cal_iter=1000L
    model_delay_filter=1

    ; turn off beam averaging (this is set as 16 in eor_defaults)
    ; beam_nfreq_avg=1

    ; recalculate all products. I'm not sure why this is set
    ; recalculate_all=1

    ; set flux cutoff for calibration sources
    calibration_flux_threshold=0.1

    ; turn off time averaging
    cal_time_average=0

    ; sidelobe subtraction
    calibration_subtract_sidelobe_catalog=filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')
    model_subtract_sidelobe_catalog=filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')
    subtract_sidelobe_catalog=filepath('GLEAM_v2_plus_rlb2019.sav',root=rootdir('FHD'),subdir='catalog_data')

    ; prevent power leakage by not gridding longer baselines (this is set as 600 in eor_defaults)
    ps_kspan=200.

    ; this is only a cal run so stop before gridding
    cal_stop=1

    ; use Ian's speedup
    use_adaptive_calibration_gain=1
    calibration_base_gain=0.5

    ; other calibration settings
    cal_reflection_mode_theory=1
    cal_mode_fit=[90,150,230,320] 
    auto_ratio_calibration=1 ; use wengyang's auto calibration
    calibration_auto_fit=0 ; this is default
    digital_gain_jump_polyfit=0 ; this is default
    end

    'Wilensky_thesis_image': begin

    ; set file pathways
    model_uv_transfer='uvfits/transfer/' + obs_id + '_model_uv_arr.sav'
    transfer_calibration = 'uvfits/transfer/' + obs_id + '_cal.sav'
    restrict_hpx_inds='EoR0_high_healpix_inds_3x.idlsave'
    
    ; apply 'Blackman-Harris^2' window function to gridding kernel
    kernel_window=1

    ; not in dictionary
    debug_dim=1

    ; set the factor at which to clip the beam model
    beam_mask_threshold=1e3

    ; turn off beam averaging (this is set at 16 in eor_defaults)
    beam_nfreq_avg=1
    
    ; don't save out calibrated visibilities
    return_cal_visibilities=0

    ; make visibilities for the subtraction model separately from the model used in calibration    
    model_visibilities=1

    ; prevent power leakage by not gridding longer baselines (this is set at 600 in eor_defaults)
    ps_kspan=200.

    ; not sure why calibration keywords are needed for gridding
    cal_time_average=0
    auto_ratio_calibration=1
    cal_bp_transfer=0

    end

  endcase

  if ~keyword_set(vis_file_list) then begin
    if keyword_set(input_folder) then begin
      vis_file_list = STRING(input_folder) + '/' + STRING(obs_id) + '.uvfits'
    endif else begin
      SPAWN, 'read_uvfits_loc.py -v ' + STRING(uvfits_version) + ' -s ' + $
        STRING(uvfits_subversion) + ' -o ' + STRING(obs_id), vis_file_list
    endelse
  endif

  undefine, uvfits_subversion, uvfits_version, input_folder

  fhd_file_list=fhd_path_setup(vis_file_list,version=version,output_directory=output_directory)
  healpix_path=fhd_path_setup(output_dir=output_directory,subdir='Healpix',output_filename='Combined_obs',version=version)


  ; Set global defaults and bundle all the variables into a structure.
  ; Any keywords set on the command line or in the top-level wrapper will supercede these defaults
  eor_wrapper_defaults,extra
  fhd_depreciation_test, _Extra=extra

  print,""
  print,"Keywords set in wrapper:"
  print,structure_to_text(extra)
  print,""

  general_obs,_Extra=extra

end
