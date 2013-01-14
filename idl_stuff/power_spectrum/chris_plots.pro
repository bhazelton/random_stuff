

pro chris_plots, no_kzero = no_kzero, pub = pub, refresh_ps = refresh_ps, refresh_binning = refresh_binning, $
               std_power = std_power, no_weighting = no_weighting, kperp_plot_range = kperp_plot_range, data_range = data_range

 froot = base_path() + 'power_spectrum/chris_maps/'
 fadd = ''
 if keyword_set(std_power) then fadd = fadd + '_sp'
 if keyword_set(no_weighting) then fadd = fadd + '_nowt'
 if keyword_set(flag_weights) then fadd = fadd + '_flagwt'

 info_file = froot + 'chris_info.txt'
 filebase = repstr(info_file, '_info.txt','') + fadd
 savefile = filebase + '_2dkpower' + '.idlsave'

 ftest = file_test(savefile) *  (1 - file_test(savefile, /zero_length))

 if ftest eq 0 or keyword_set(refresh_ps) or keyword_set(refresh_binning) then $
    reg_grid_3dps, info_file, std_power = std_power, no_weighting = no_weighting, refresh = refresh, no_kzero = no_kzero, $
                   flag_weights = flag_weights

 restore, savefile

 plotfile_path = base_path() + 'power_spectrum/plots/'
 plotfile_base = plotfile_path + 'chris' + fadd
 if keyword_set(no_kzero) then plotfile_base = plotfile_base + '_nok0'
 plotfile = plotfile_base + '_kspace_power'
 weight_plotfile = plotfile_base + '_kspace_weights'
     
 data_range = [1e9, 1e25]
 ;;weights_range = [1e3, 1e8]
 ;;kperp_plot_range = [min(kperp_edges), 0.005]
 kpower_2d_plots, savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                  data_range = data_range, pub = pub, plotfile = plotfile
 kpower_2d_plots, savefile, /plot_weights, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                  pub = pub, plotfile = weight_plotfile, window_num = 2, title = 'Weights'
 


end
