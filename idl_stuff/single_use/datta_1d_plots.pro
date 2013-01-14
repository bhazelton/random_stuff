

pro datta_1d_plots, refresh_binning = refresh_binning, refresh_fft = refresh_fft, plot_weights = plot_weights, psf = psf, pub = pub, $
                    no_kzero = no_kzero, edge_on_grid = edge_on_grid, match_datta = match_datta, std_power = std_power, $
                    no_weighting = no_weighting, window_num = window_num, from_3d = from_3d, compare = compare

  if keyword_set(match_datta) then edge_on_grid = 1

  savefile = base_path() + 'single_use/datta_data/'
  if keyword_set(psf) then savefile = savefile + 'MWA_PSF' else savefile = savefile + 'MWA_GSM_0.1err'
  
  if keyword_set(compare) then begin
     savefile = [savefile + '_sp_nowt_1dkpower.idlsave', savefile + '_nok0_1dkpower.idlsave', $
                 base_path()+'single_use/datta_data/eor_power_1d.idlsave']

     if keyword_set(refresh_fft) then begin
        datta_3d, /refresh, psf = psf, /no_kzero, edge_on_grid = edge_on_grid, match_datta = match_datta, /quiet
        datta_3d, /refresh, psf = psf, edge_on_grid = edge_on_grid, match_datta = match_datta, /std_power, $
                  /no_weighting, /quiet
     endif else if not keyword_set(from_3d) then begin
        test_save = file_test(savefile[0]) *  (1 - file_test(savefile[0], /zero_length))
        if test_save eq 0 or keyword_set(refresh_binning) then datta_3d, psf = psf, /no_kzero, edge_on_grid = edge_on_grid, $
           match_datta = match_datta, /quiet

        test_save = file_test(savefile[1]) *  (1 - file_test(savefile[1], /zero_length))
        if test_save eq 0 or keyword_set(refresh_binning) then datta_3d, psf = psf, edge_on_grid = edge_on_grid, $
           match_datta = match_datta, /std_power, /no_weighting, /quiet
     endif

  endif else begin
     fadd = ''
     if keyword_set(std_power) then fadd = fadd + '_sp'
     if keyword_set(no_weighting) then fadd = fadd + '_nowt'
     if keyword_set(no_kzero) then fadd = fadd + '_nok0'
     if keyword_set(edge_on_grid) then fadd = fadd + '_edgegrid'
     if keyword_set(match_datta) then fadd = fadd + '_match'
     savefile = savefile + fadd + '_1dkpower.idlsave'

     if keyword_set(refresh_fft) then datta_3d, /refresh, psf = psf, no_kzero = no_kzero, edge_on_grid = edge_on_grid, $
        match_datta = match_datta, std_power = std_power, no_weighting = no_weighting, /quiet $
     else if not keyword_set(from_3d) then begin
        test_save = file_test(savefile) *  (1 - file_test(savefile, /zero_length))
        if test_save eq 0 or keyword_set(refresh_binning) then datta_3d, psf = psf, no_kzero = no_kzero, edge_on_grid = edge_on_grid, $
           match_datta = match_datta, std_power = std_power, no_weighting = no_weighting, /quiet
     endif
  endelse


  plotfile = base_path() + 'single_use/datta_plots/datta'
  if keyword_set(psf) then plotfile = plotfile + '_psf'
  if keyword_set(compare) then plotfile = plotfile +'_1d_power_compare' else plotfile = plotfile +'_1d_power' + fadd
 
  data_range = 10^[-3d, 5d]

  if keyword_set(compare) then begin
     kpower_1d_plots, savefile, plot_weights = plot_weights, data_range = data_range, pub = pub, plotfile = plotfile, $
                      window_num = window_num, names = ['Naive PS', 'Weighted PS', 'EoR Signal'], colors = [254, 0, 75]
  endif else kpower_1d_plots, savefile, plot_weights = plot_weights, data_range = data_range, pub = pub, plotfile = plotfile, $
                              window_num = window_num

end
