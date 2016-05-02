pro plot_josh_ps, window_num = window_num,  png = png, pdf = pdf, eps = eps

  file_dir = '/Users/bryna/Dropbox/Documents/Physics roaming/MWAPipelinepaper_pspec_plots/Dillon_Files/'
  power_file = 'P_k_.txt'
  kperp_edges_file = 'kPerpBinEdges.txt'
  kpar_edges_file = 'kParaBinEdges.txt'

  textfast, power, file_path = file_dir + power_file, /read
  textfast, kperp_edges, file_path = file_dir + kperp_edges_file, /read
  textfast, kpar_edges, file_path = file_dir + kpar_edges_file, /read

  kperp_edges = [kperp_edges[0,0], reform(kperp_edges[1, *])]
  kpar_edges = [kpar_edges[0,0], reform(kpar_edges[1, *])]

  ; ------   assumptions:
  n_freq = 192
  f_delta = .16 ; in MHz
  frequencies = findgen(n_freq)*f_delta + 167.155 ; in MHz
  ; ------


  z0_freq = 1420.40d ;; MHz
  redshifts = z0_freq/frequencies - 1d
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los, hubble_param = hubble_param

  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  z_mpc_delta = float(mean(comov_los_diff))
  z_mpc_mean = float(mean(comov_dist_los))

  kperp_lambda_conv = z_mpc_mean / (2.*!pi)

  delay_delta = 1e9/(n_freq*f_delta*1e6) ;; equivilent delay bin size for kparallel
  delay_max = delay_delta * n_freq/2.    ;; factor of 2 b/c of neg/positive
  delay_params = [delay_delta, delay_max]

  ; convert into actual Mpc because that's what eppsilon expects
  power = power / (hubble_param)^3d
  kperp_edges = kperp_edges * hubble_param
  kpar_edges = kpar_edges * hubble_param


  ;; assume 20 degrees from pointing center to first null
  mean_redshift = mean(redshifts)

  cosmology_measures, mean_redshift, wedge_factor = wedge_factor
  source_dist = 20d * !dpi / 180d
  fov_amp = wedge_factor * source_dist

  ;; use 90 degrees from pointing center to horizon
  horizon_amp = wedge_factor * (90d * !dpi / 180d)

  wedge_amp = [fov_amp, horizon_amp]

  kperp_plot_range = ([5.5, 65.]/kperp_lambda_conv)/hubble_param
  force_kperp_axis_range = ([5.5, 300.]/kperp_lambda_conv)/hubble_param
  kpar_plot_range=[-1,1.74]
  ;force_kpar_axis_range = [1e-3,1.74*2]

  if keyword_set(png) or keyword_set(eps) or keyword_set(pdf) then plotfile = file_dir + 'Dillon_thru_eppsilon_2dkpower'
  data_range = [1e6, 1e14]
  cb_margin = [0.3, 0.02]

  kpower_2d_plots, power = power, kperp_edges = kperp_edges, kpar_edges = kpar_edges, $
    delay_params = delay_params, kperp_lambda_conv = kperp_lambda_conv, wedge_amp = wedge_amp, $
    hubble_param = hubble_param, kperp_plot_range = kperp_plot_range, kpar_plot_range=kpar_plot_range, $
    force_kperp_axis_range = force_kperp_axis_range, force_kpar_axis_range = force_kpar_axis_range, data_range = data_range, $
    /baseline_axis, /delay_axis, /hinv, /plot_wedge_line, window_num = window_num, $
    png = png, pdf = pdf, eps = eps, plotfile = plotfile, cb_margin=cb_margin

end