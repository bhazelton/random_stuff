

pro simulation_wrapper, t32 = t32, define_baselines = define_baselines, baseline_spacing = baseline_spacing, $
                        baseline_layout = baseline_layout, fine_res = fine_res, use_outliers = use_outliers, $
                        refresh_sim = refresh_sim, eor_refresh = eor_refresh, off_axis = off_axis, refresh_hmf = refresh_hmf, $
                        full_sky = full_sky, source_radius = source_radius

  if keyword_set(refresh_hmf) then refresh_sim = 1
  if keyword_set(refresh_sim) then eor_refresh = 1


  if keyword_set(define_baselines) then begin
     if n_elements(baseline_spacing) eq 0 then baseline_spacing = 4
     if n_elements(baseline_layout) eq 0 then baseline_layout = 'simple'

     x_offsets = [3.16, 6.45, 12.65, 25.09, 49.91, 75.3, 100.04, 124.58, 250.26, 499.18]
     names = string([3, 6, 13, 25, 50, 75, 100, 125, 250, 500], format = '(i04)')

     spacing_name = number_formatter(string(baseline_spacing), flag = 2)
     tag = baseline_layout + spacing_name + '_'

  endif else if keyword_set(t32) then begin
     x_offsets = [3.16, 6.45, 12.65, 25.09, 49.91, 75.3, 100.04, 124.58, 250.26, 499.18]
     names = string([3, 6, 13, 25, 50, 75, 100, 125, 250, 500], format = '(i04)')
     tag = '32t_'
  endif else begin
     x_offsets = [25.3, 51.6, 101.2, 200.7, 399.3, 602.4, 800.3, 998.8]
     names = string([25, 50, 100, 200, 400, 600, 800, 1000], format = '(i04)')
     
     if keyword_set(use_outliers) then tag = '512t_' else tag = '496t_'
     if keyword_set(fine_res) then tag = tag + 'fine_'
  endelse

   if keyword_set(off_axis) then begin
     y_offsets = -1. * x_offsets[2]
     x_offsets = x_offsets[0]
     names = 'offaxis'
  endif

   sim_froot = base_path('data') + 'fhd_simulations_old/'
   noise_file = sim_froot + 'sim_' + tag + 'noise_uvf.idlsave'
   test_noise = file_test(noise_file)  *  (1 - file_test(noise_file, /zero_length))
   
   eor_file = sim_froot + 'sim_' + tag + 'eor_uvf.idlsave'
   test_eor = file_test(eor_file)  *  (1 - file_test(eor_file, /zero_length))
   
   info_file = sim_froot + 'sim_' + tag + 'info.idlsave'
   weights_file = sim_froot + 'sim_' + tag + 'weights.idlsave'
   
   if keyword_set(full_sky) then begin

     restore, info_file
     if keyword_set(t32) then xy_length = 1024 else xy_length = 2048

     if n_elements(source_radius) eq 0 then begin
        good_radius_deg = ceil(xy_length * degpix / sqrt(2d))
        names = 'fullsky'
     endif else begin
        good_radius_deg = source_radius
        names = 'fullsky_' + number_formatter(good_radius_deg) + 'deg'
     endelse

     good_radius_pix = good_radius_deg / degpix
     npix_good = !dpi * good_radius_pix^2d
     good_area = !dpi * good_radius_deg^2d
     
     source_counts_6C_catalog,hist_norm,flux_bin, binsize=log_binsize
     log_flux = alog10(flux_bin)
     log_hist = alog10(hist_norm)
     
     ;; get good range for fitting
     wh_fit = where(log_flux lt 1. and log_flux gt (-0.5), count_fit)
     if count_fit lt 3 then stop
     lin_params = linfit(log_flux[wh_fit], log_hist[wh_fit])
     slope = lin_params[1]
     min_log_flux = alog10(0.1d)
     max_log_flux = alog10(1d)
     prob_range = minmax([slope*min_log_flux, slope*max_log_flux])
     
     npts = round(npix_good)
     rand_vals = randomu(seed, 4, npts)
     log_flux_vals = rand_vals[0,*]*(max_log_flux-min_log_flux)+min_log_flux
     test_vals = rand_vals[1,*]*(prob_range[1]-prob_range[0])+prob_range[0]
     x_vals = rand_vals[2,*]*(2d*good_radius_pix)-good_radius_pix
     y_vals = rand_vals[3,*]*(2d*good_radius_pix)-good_radius_pix
     radii = sqrt(x_vals^2d + y_vals^2d)

     prob_vals = slope * log_flux_vals
     
     good_pts = where(test_vals le prob_vals and radii le good_radius_pix $
                      and abs(x_vals) le xy_length/2d and abs(y_vals) le xy_length/2d, count_good)
     ;;cghistoplot, flux_vals[good_pts], binsize=0.1
     if count_good gt 5000 then good_pts = good_pts[0:5000]
     
     x_offsets = x_vals[good_pts]
     y_offsets = y_vals[good_pts]
     fluxes = 10d^log_flux_vals[good_pts]

  endif
  
  sim_files = sim_froot + 'sim_' + tag + names + '_uvf.idlsave'

  nsims = n_elements(sim_files)
  for i = 0, nsims-1 do begin
     test = file_test(sim_files[i])  *  (1 - file_test(sim_files[i], /zero_length))

     if i eq 0 and test_noise eq 0 then begin
        if test eq 1 and not keyword_set(refresh_sim) then no_sim = 1

        test = 0
        noise = 1
     endif else noise = 0

     if i eq 1 and (test_eor eq 0 or keyword_set(eor_refresh)) then begin
        if test eq 1 and not keyword_set(refresh_sim) then no_sim = 1

        test = 0
        eor = 1
     endif else eor = 0

     if keyword_set(refresh_hmf) and i eq 0 then new_hmf = 1 else new_hmf = 0

     if test eq 0 or keyword_set(refresh_sim) then begin
        if keyword_set(full_sky) then begin
           x_offset = x_offsets
           y_offset = y_offsets
           flux = fluxes
        endif else begin
           x_offset = x_offset[i]
           if keyword_set(offaxis) then y_offset = y_offset[i]
        endelse

        fhd_simulations, x_offset = x_offset, y_offset = y_offset, flux = flux, sim_file = sim_files[i], info_file = info_file, $
                         weights_file = weights_file, t32 = t32, noise_file = noise_file, noise_gen = noise, eor_gen = eor, $
                         eor_file = eor_file, define_baselines = define_baselines, baseline_spacing = baseline_spacing, $
                         baseline_layout = baseline_layout, no_sim = no_sim, fine_res = fine_res, use_outliers = use_outliers, $
                         refresh = new_hmf
     endif
 endfor


end
