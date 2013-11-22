pro healpix_wrapper, rts = rts, version = version, refresh_dft = refresh_dft, refresh_ps = refresh_ps, dft_ian = dft_ian, $
    refresh_binning = refresh_binning, pol_inc = pol_inc, sim = sim, $
    no_spec_window = no_spec_window, spec_window_type = spec_window_type, noise_sim = noise_sim, $
    cut_image = cut_image, individual_plots = individual_plots, plot_filebase = plot_filebase, pub = pub, $
    kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis
    
  ;; The only required input is the datafile name (including the full path)
    
  if keyword_set(rts) then begin
    froot = base_path('data') + 'rts_data/test2/'
    
    data_dir = froot + 'BdaggerV/'
    datafiles = file_search(data_dir + '*.fits')
    
    weights_dir = froot + 'Bdagger1/'
    weightfiles = file_search(weights_dir + '*.fits')
    
    variance_dir = froot + 'BdaggerB/'
    variancefiles = file_search(variance_dir + '*.fits')
    
    datafile =  rts_fits2imagecube(datafiles, weightfiles, variancefiles, pol_inc, save_path = froot)
    
  endif else if keyword_set(sim) then begin
    datafile = base_path('data') + 'fhd_sim_data/fhd_v300/Healpix/Sim_obs_' + ['even','odd']+ '_cube.sav'
    
  endif else begin
    ;;datafile = base_path('data') + 'fhd_ps_data/X16_EOR1/multi_freq_residuals_cube_healpix.sav'
    ;;datafile = base_path('data') + 'fhd_ps_data/X16_EOR1/fhd_v101/Combined_obs_EOR1_P00_145_20110926193959-' + $
    ;;           'EOR1_P00_145_20110926193959_' + $
    ;;datafile = base_path('data') + 'fhd_ps_data/X16_EOR1/fhd_v15/Combined_obs_EOR1_P00_145_20110926193959-' + $
    ;;           'EOR1_P00_145_20110926200503_' + $
    ;;datafile = base_path('data') + 'fhd_ps_data/128T_cubes/firstpass/Combined_obs_1061316296-1061316296_' + $
    datafile = base_path('data') + 'fhd_ps_data/128T_cubes/deconvolved/Combined_obs_1061316296-1061316296_' + $
      ['even','odd']+ '_cube.sav'
  endelse
  
  
  
  ;; dft_fchunk applies only to Healpix datasets (it's ignored otherwise) and it specifies how many frequencies to process
  ;;   through the dft at once. This keyword allows for trade-offs between memory use and speed.
  ;; The optimum setting varies by computer and the speed is NOT a linear function of this parameter
  ;;   (it's not even monotonic) so some experimenting is required. The memory required is approximately linear --
  ;;   the higher this value the more memory required.
  ;; The maximum value of this parameter is the number of frequency slices in the cube
  ;;   (if it's set too large it will be reduced to the maximum)
  
  dft_fchunk = 1
  
  ;; save_path specifies a location to save the power spectrum files.
  ;; This is also where the code looks for intermediate save files to avoid re-running code.
  ;; If this is parameter is not set, the files will be saved in the same directory as the datafile.
  
  ;; the following sets the save_path to a 'psv'+version directory inside the datafile[0] directory and
  ;; creates the directory if it doesn't exist
  if n_elements(version) gt 0 then begin
    save_path = file_dirname(datafile[0], /mark_directory) + 'psv' + number_formatter(version) + path_sep()
    if not file_test(save_path, /directory) then file_mkdir, save_path
  endif
  
  if keyword_set(rts) then std_savepath = base_path('data') + 'rts_data/' $
  else if keyword_set(sim) then std_savepath = base_path('data') + 'fhd_sim_data/' else std_savepath = base_path('data') + 'fhd_ps_data/'
  
  if n_elements(save_path) gt 0 then begin
    pos = strpos(save_path, std_savepath)
    if pos ne -1 then save_path_ext = strmid(save_path, pos + strlen(std_savepath)) else save_path_ext = ''
  endif else begin
    pos = strpos(file_dirname(datafile[0], /mark_directory), std_savepath)
    if pos ne -1 then save_path_ext = strmid(file_dirname(datafile[0], /mark_directory), pos + strlen(std_savepath)) $
    else save_path_ext = ''
  endelse
  
  ;; savefilebase specifies a base name to use for the save files
  
  
  ;; plot_path specifies a location to save plot files.
  ;; If this parameter is not set, the plots will be saved in the same directory as the datafile.
  if keyword_set(rts) then plot_path = base_path('plots') + 'power_spectrum/rts_data/' + save_path_ext $
  else if keyword_set(sim) then plot_path = base_path('plots') + 'power_spectrum/fhd_sim/' + save_path_ext $
else plot_path = base_path('plots') + 'power_spectrum/fhd_data/' + save_path_ext
if not file_test(plot_path, /directory) then file_mkdir, plot_path

  ;; plot_filebase specifies a base name to use for the plot files


;; freq_ch_range specifies which frequency channels to include in the power spectrum.
;; Fewer number of channels makes the dfts faster
;;freq_ch_range = [0, 191]
;; freq_ch_range = [288, 480]
;;freq_ch_range = [575, 767]


;; pol_inc specifies which polarizations to generate the power spectra for.
;; The default is ['xx,'yy']
;;pol_inc = ['xx']


;; cut_image keyword only applies to Healpix datasets. It allows for limiting the field of view in the
;; image plane to match calculated k-modes (centered on image center).
;; Currently defaults to on. Set equal to 0 to turn it off, 1 to turn it on

;; There are 3 refresh flags to indicate that various stages should be re-calculated
;;   (rather than using previous save files if they exist).
;; If an early stage is recalculated, all subsequent stages will also be recalculated
;; The earliest stage is refresh_dft, which is only used for Healpix datasets (it's ignored otherwise)
;; The next stage is refresh_ps and the last stage is refresh_binning.
;; To set any of these flags, set them equal to 1 (true)

;;refresh_dft=1
;;refresh_ps = 1
;;refresh_binning=1

;; options for spectral windowing:
;; available window funtions are: ['Hann', 'Hamming', 'Blackman', 'Nutall', 'Blackman-Nutall', 'Blackman-Harris']
;; Default is to use Blackman-Harris, for no spectral windowing set no_spec_window = 1
;; To use another window type use the spec_window_type keyword, eg spec_window_type = 'hann'

;; options for binning:
;; log_kperp, log_kpar and log_k1d are flags: set to 1 (true) for logarithmic bins
;; kperp_bin, kpar_bin and k1d_bin take scalar values to control bin sizes.
;;   (The actual binsize for linear binning and the log binsize for log binning -- bins per decade = 1/binsize)

;; log_kperp = 1
;; log_kpar = 1


;; options for plotting:
;; kperp_linear_axis is a flag, set to 1 to use a linear kperp axis (default is log axis)
;; kpar_linear_axis is a flag, set to 1 to use a linear kpar axis (default is log axis)
;; data_range specifies the min & max value of the signal colorbar (values outside that range are clipped to those values)
;; sigma_range, nev_range, snr_range, noise_range, nnr_range control the other colorbar ranges
;; baseline_axis is a flag (defaulted to true) to mark baseline; length along top axis of 2d plots (set to 0 to turn off)
;; delay_axis is a flag (defaulted to true) to mark delay time along right axis of 2d plots (set to 0 to turn off)
;; hinv is a flag (defaulted to true) to use h^-1 Mpc rather than physical Mpc in plot units (set to 0 to turn off)
;; plot_wedge_line is a flag (defaulted to true) to plot a line marking the wedge (both horizon & FoV) (set to 0 to turn off)
;; grey_scale is a flag to use a black/white color scale rather than the default color scale
;; pub is a flag to make save plots as eps files rather than displaying to the screen

;; kperp_linear_axis = 1
;; kpar_linear_axis = 1

data_range = [5e-2, 5e7]
noise_range = [5e3, 5e5]
nnr_range = [1e-2, 2]
snr_range = [5e-5, 5e3]


fhd_data_plots, datafile, dft_fchunk=dft_fchunk, plot_path = plot_path, plot_filebase = plot_filebase, save_path = save_path, savefilebase = savefilebase, $
  pol_inc = pol_inc, rts = rts, $
  refresh_dft = refresh_dft, refresh_ps = refresh_ps, refresh_binning = refresh_binning, $
  freq_ch_range = freq_ch_range, no_spec_window = no_spec_window, spec_window_type = spec_window_type, $
  noise_sim = noise_sim, cut_image = cut_image, dft_ian = dft_ian, $
  log_kpar = log_kpar, log_kperp = log_kperp, kpar_bin = kpar_bin, kperp_bin = kperp_bin, $
  log_k1d = log_k1d, k1d_bin = k1d_bin, kperp_linear_axis = kperp_linear_axis, kpar_linear_axis = kpar_linear_axis, $
  data_range = data_range, sigma_range = sigma_range, nev_range = nev_range, snr_range = snr_range, noise_range = noise_range, nnr_range = nnr_range, $
  baseline_axis = baseline_axis, delay_axis = delay_axis, hinv = hinv, $
  plot_wedge_line = plot_wedge_line, grey_scale = grey_scale, individual_plots = individual_plots, pub = pub
  
end
