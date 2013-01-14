

pro datta_2d_plots, refresh_binning = refresh_binning, refresh_fft = refresh_fft, plot_weights = plot_weights, psf = psf, $
                    pub = pub, no_kzero = no_kzero, edge_on_grid = edge_on_grid, match_datta = match_datta, $
                    std_power = std_power, no_weighting = no_weighting, window_num = window_num, from_3d = from_3d

  if keyword_set(match_datta) then edge_on_grid = 1

  savefile = base_path() + 'single_use/datta_data/'
  if keyword_set(psf) then savefile = savefile + 'MWA_PSF' else savefile = savefile + 'MWA_GSM_0.1err'
  
  fadd = ''
  if keyword_set(std_power) then fadd = fadd + '_sp'
  if keyword_set(no_weighting) then fadd = fadd + '_nowt'
  if keyword_set(no_kzero) then fadd = fadd + '_nok0'
  if keyword_set(edge_on_grid) then fadd = fadd + '_edgegrid'
  if keyword_set(match_datta) then fadd = fadd + '_match'
  savefile = savefile + fadd + '_2dkpower.idlsave'

  if keyword_set(refresh_fft) then datta_3d, /refresh, psf = psf, no_kzero = no_kzero, edge_on_grid = edge_on_grid, $
     match_datta = match_datta, std_power = std_power, no_weighting = no_weighting, /quiet $
  else if not keyword_set(from_3d) then begin
     test_save = file_test(savefile) *  (1 - file_test(savefile, /zero_length))
     if test_save eq 0 or keyword_set(refresh_binning) then datta_3d, psf = psf, no_kzero = no_kzero, edge_on_grid = edge_on_grid, $
        match_datta = match_datta, std_power = std_power, no_weighting = no_weighting, /quiet
  endif

  if keyword_set(match_datta) then begin
     restore, savefile
     kperp_plot_range = minmax(kperp_edges)
     kperp_plot_range[1] = 2d
     ;;kperp_plot_range[1] = 2d
  
     ;;kpar_plot_range = minmax(kpar_edges)
     ;;kpar_plot_range[1] = 2d

     data_range = [10^(-2d), 10^8d]
  endif
  
  plotfile = base_path() + 'single_use/datta_plots/datta'
  if keyword_set(psf) then plotfile = plotfile + '_psf'
  plotfile = plotfile +'_kspace_power' + fadd

  kpower_2d_plots, savefile, plot_weights = plot_weights, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                   data_range = data_range, pub = pub, plotfile = plotfile, window_num = window_num
end
