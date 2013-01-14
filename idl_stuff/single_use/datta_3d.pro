


pro datta_3d, refresh = refresh, psf = psf, no_kzero = no_kzero, mask = mask, pixelwise_mask = pixelwise_mask, $
              k1_mask = k1_mask, k2_mask = k2_mask, k3_mask = k3_mask, edge_on_grid = edge_on_grid, match_datta = match_datta, $
              std_power = std_power, no_weighting = no_weighting, quiet = quiet

  froot =  base_path() + 'single_use/datta_data/'

  ;; first work out psf parameters
  psf_fitsfile = froot + 'MWA_PSF.fits'
  psf_header = headfits(psf_fitsfile)
  psf_dims = sxpar(psf_header, 'NAXIS*')

  ;; get number of dimensions (larger than 1)
  wh_dims = where(psf_dims gt 1, n_psfdims)

  if (n_psfdims ne 3) then message, 'psf data is not 3 dimensional'

  psf_dim_type = sxpar(psf_header,'CTYPE*')
  psf_dim_units = sxpar(psf_header, 'CUNIT*')
  psf_ref_pixels = sxpar(psf_header,'CRPIX*') ;; these are 1-based NOT 0-based
  psf_ref_pix_vals = sxpar(psf_header,'CRVAL*')
  psf_ref_pix_deltas = sxpar(psf_header,'CDELT*')

  ;; drop any shallow dimensions
  if (n_psfdims ne n_elements(psf_dims)) then begin
     psf_dims = psf_dims[wh_dims]
     psf_dim_type = psf_dim_type[wh_dims]
     psf_dim_units = psf_dim_units[wh_dims]
     psf_ref_pixels = psf_ref_pixels[wh_dims]
     psf_ref_pix_vals = psf_ref_pix_vals[wh_dims]
     psf_ref_pix_deltas = psf_ref_pix_deltas[wh_dims]
  endif

  ;; then work out image parameters
  fitsfile = froot + 'MWA_GSM_0.1err.fits'

  header = headfits(fitsfile)
  dims = sxpar(header, 'NAXIS*')

  ;; get number of dimensions (larger than 1)
  wh_dims = where(dims gt 1, n_dims)

  if (n_dims ne 3) then message, 'Data is not 3 dimensional'

  dim_type = sxpar(header,'CTYPE*')
  dim_units = sxpar(header, 'CUNIT*')
  ref_pixels = sxpar(header,'CRPIX*') ;; these are 1-based NOT 0-based
  ref_pix_vals = sxpar(header,'CRVAL*')
  ref_pix_deltas = sxpar(header,'CDELT*')

  ;; drop any shallow dimensions
  if (n_dims ne n_elements(dims)) then begin
     dims = dims[wh_dims]
     dim_type = dim_type[wh_dims]
     dim_units = dim_units[wh_dims]
     ref_pixels = ref_pixels[wh_dims]
     ref_pix_vals = ref_pix_vals[wh_dims]
     ref_pix_deltas = ref_pix_deltas[wh_dims]
  endif

  ;; check the dimensional ordering
  x_dim = where(strpos(strlowcase(dim_type), 'ra') ne -1, count)
  if count ne 1 or x_dim[0] ne 0 then message, 'First dimension must be contain "ra".'
  if strpos(strlowcase(dim_units[x_dim]), 'deg') eq -1 then message, 'RA dimension does not have units of degrees'

  y_dim = where(strpos(strlowcase(dim_type), 'dec') ne -1, count)
  if count ne 1 or y_dim[0] ne 1 then message, 'Second dimension must be contain "dec".'
  if strpos(strlowcase(dim_units[y_dim]), 'deg') eq -1 then message, 'DEC dimension does not have units of degrees'

  z_dim = where(strpos(strlowcase(dim_units), 'hz') ne -1, count)
  if count ne 1 or z_dim[0] ne 2 then message, 'Third dimension must be in units of hz.'

  ;;check that coordinates are the same on the image and psf
  if total(dims - psf_dims) ne 0 then stop
  if total(ref_pix_vals - psf_ref_pix_vals) ne 0 then stop
  if total(ref_pix_deltas - psf_ref_pix_deltas) ne 0 then stop

  ;; conversion to comoving Mpc
  z0_freq = 1420.40 ;; MHz
  freq_vals = dindgen(dims[2]) * ref_pix_deltas[2] + ref_pix_vals[2]

  ;; convert to MHz
  if (strtrim(strlowcase(dim_units[2])) eq 'hz') then freq_vals = freq_vals / (1d*1e6) $
  else if (strtrim(strlowcase(dim_units[2])) eq 'khz') then freq_vals = freq_vals / (1d*1e3) $
  else if (strtrim(strlowcase(dim_units[2])) ne 'mhz') then message, 'Unknown units on frequency dimension (not Hz, kHz or MHz)'
  redshifts = z0_freq/freq_vals - 1
  comoving_distance_los, redshifts, comov_dist_los


  ;; now calculate k values for FFT (in inverse Mpc) 
  ;; factors of 2 b/c image was oversampled by a factor of 2
  x_rad_delta = abs(ref_pix_deltas[0]) * !pi / 180d
  n_kx = dims[0] / 2
  x_rad_length = dims[0] * abs(ref_pix_deltas[0]) * !pi / 180d
  x_mpc_delta = x_rad_delta * mean(comov_dist_los)
  x_mpc_length = x_rad_length * mean(comov_dist_los)
  kx_mpc_range = (2d*!pi) / (2*x_mpc_delta)
  kx_mpc_delta = (2d*!pi) / x_mpc_length
  kx_mpc_fft = (dindgen(n_kx)-n_kx/2+1) * kx_mpc_delta

  y_rad_delta = abs(ref_pix_deltas[1]) * !pi / 180d
  n_ky = dims[1] / 2
  y_rad_length = dims[1] * abs(ref_pix_deltas[1]) * !pi / 180d
  y_mpc_delta = y_rad_delta * mean(comov_dist_los)
  y_mpc_length = y_rad_length * mean(comov_dist_los)
  ky_mpc_range = (2d*!pi) / (2*y_mpc_delta)
  ky_mpc_delta = (2d*!pi) / y_mpc_length
  ky_mpc_fft = (dindgen(n_ky)-n_ky/2+1) * ky_mpc_delta

  kx_ind_range = [n_kx/2, 3*n_kx/2-1]
  ky_ind_range = [n_ky/2, 3*n_ky/2-1]

  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  
  z_mpc_delta = mean(comov_los_diff)
  z_mpc_length = max(comov_dist_los) - min(comov_dist_los)+mean(comov_los_diff)
  kz_mpc_range =  (2d*!pi) / z_mpc_delta
  kz_mpc_delta = (2d*!pi) / z_mpc_length
  n_freq = dims[2]
  ;; we will drop the negative kz elements, so only half as many as frequencies.
  if (n_freq mod 2) eq 0 then n_kz = n_freq / 2 else n_kz = (n_freq + 1)/2
  kz_mpc = dindgen(n_kz) * kz_mpc_delta


  temp = strpos(fitsfile, '/', /reverse_search)
  img_filebase = strsplit(strmid(fitsfile, temp+1), '.fits', /extract, /regex)
  ;;if n_elements(img_filebase) gt 1 then stop

  temp = strpos(psf_fitsfile, '/', /reverse_search)
  psf_filebase = strsplit(strmid(psf_fitsfile, temp+1), '.fits', /extract, /regex)
  ;;if n_elements(psf_filebase) gt 1 then stop

  if keyword_set(psf) then filebase = psf_filebase[0] $
  else filebase = img_filebase[0]

  if keyword_set(std_power) then filebase = filebase + '_sp'

  savefile = froot + filebase + '_power.idlsave'
  test_save = file_test(savefile) *  (1 - file_test(savefile, /zero_length))
  if test_save ne 0 and not keyword_set(refresh) then begin
     print, 'Save file exists, restoring now.'
     restore, savefile
     refresh = 0

     if n_elements(power_3d) ne 0 then save_dims = size(power_3d, /dimensions) $
     else save_dims = size(save_power,/dimensions)
     dims_expected = [n_kx, n_ky, n_kz]
     if total(save_dims-dims_expected) ne 0 then refresh = 1
  endif else refresh = 1

  if refresh eq 1 then begin

     ;; convert to mK (from Jy/beam)
     ;; first get pixel size in strad.
     pix_size_str = abs(ref_pix_deltas[y_dim]) * abs(ref_pix_deltas[x_dim]) * (!pi / 180d)^2
     pix_size_str = pix_size_str[0]
     ;; powers of ten are from: Jy, c^2, mK, MHz^2, kB
     conv_factor = (10^(double(-26+16+3-12+23)) * 9d) / (pix_size_str * 2d * freq_vals^2 * 1.38)

     times = dblarr(n_freq+1)
     times[0] = systime(1)

     uvf_cube = make_array(n_kx, n_ky, n_freq,/dcomplex)
     sigma2_cube = make_array(n_kx, n_ky, n_freq,/dcomplex)

     for i=0, n_freq-1 do begin      

        if i eq 0 then print, 'Array allocation time: ' + string(systime(1) - times[0], format = '(f7.2)') $
        else begin
           times[i] = systime(1)
           time_per_bin = (times[i] - times[0])/i
           time_last_bin = times[i]-times[i-1]
           if (time_last_bin gt 3*time_per_bin) or i mod 16 eq 0 then $
           print, 'On freq bin ' + string(i, format = '(i3.3)') + ', time per bin: ' + string(time_per_bin, format = '(f7.2)') + $
                  ', time last bin ' + string(time_last_bin, format = '(f7.2)')
        endelse

        psf_image = double(readfits(psf_fitsfile, nslice = i, /silent))
        if psf_ref_pix_deltas[0] lt 0 then psf_image = reverse(psf_image)
        if psf_ref_pix_deltas[1] lt 0 then psf_image = reverse(psf_image, 2)
        psf_image = psf_image * conv_factor[i]

        psf_ft = FFT(psf_image, /double)
        psf_image=0
  
        ;; now need to shift the FFT to put the zero at the center
        center = size(psf_ft, /dimension)/2.-1
        if (total(ceil(center)-floor(center)) ne 0) then stop $
        else center = floor(center)
        psf_ft = shift(temporary(psf_ft), center)
        psf_ft = psf_ft[kx_ind_range[0]:kx_ind_range[1], ky_ind_range[0]:ky_ind_range[1]]

        weight = psf_ft
        min_val = 0
        wh = where(abs(weight) le min_val, count)
 
        sigma2 = abs(1/weight)^2d
        if count gt 0 then sigma2[wh] = 0

        if keyword_set(psf) then begin
           uvf = psf_ft / weight
           if count gt 0 then uvf[wh] = 0
        endif else begin
           image = double(readfits(fitsfile, nslice = i, /silent))     
           if ref_pix_deltas[0] lt 0 then image = reverse(image)
           if ref_pix_deltas[1] lt 0 then image = reverse(image, 2)
           image = image * conv_factor[i]

           ;; do ft in 2 spacial directions
           uv = FFT(image, /double)
           image = 0
           uv = shift(temporary(uv), center)
           uv = uv[kx_ind_range[0]:kx_ind_range[1], ky_ind_range[0]:ky_ind_range[1]]

           uvf =  uv / weight
           if count gt 0 then uvf[wh] = 0

           ;; threshold = 100
           ;; wh = where(uvf gt threshold, count)
           ;; if count gt 0 then uvf[wh] = 0
        endelse
;;stop
        psf_ft = 0
        sigma=0
        wh=0
        uv=0

        uvf_cube[*,*,i] = temporary(uvf)
        sigma2_cube[*,*,i] = temporary(sigma2)
     endfor
     times[n_freq] = systime(1)
     time_per_bin = (times[n_freq] - times[0])/n_freq
     time_last_bin = times[n_freq]-times[n_freq-1]
     print, 'finished freq loop ' + string(i, format = '(i3.3)') + ', time per bin: ' + string(time_per_bin, format = '(f7.2)') + $
                  ', time last bin ' + string(time_last_bin, format = '(f7.2)')

     ;; do fft in third dimension and drop the negative kparallels because the images are real

     ft = FFT(temporary(uvf_cube), dimension=3)
     ft = ft[*,*,0:n_kz-1]
     
     if keyword_set(std_power) then begin
        ;; for standard power calc. just need ft of sigma2
        sigma2_ft = FFT(temporary(sigma2_cube), dimension=3)
        sigma2_ft = sigma2_ft[*,*,0:n_kz-1]
     endif else begin
        ;; for new power calc, need cos2, sin2, cos*sin transforms
        
        ;; have to do this in a for loop for memory's sake
        covar_cos = dblarr(n_kx, n_ky, n_kz)
        covar_sin = dblarr(n_kx, n_ky, n_kz)
        covar_cross = dblarr(n_kx, n_ky, n_kz)
        
        ;; comov_dist_los goes from large to small z
        z_relative = dindgen(n_freq)*z_mpc_delta
        freq_kz_arr = rebin(reform(rebin(reform(kz_mpc, 1, n_kz), n_freq, n_kz) * $
                                   rebin(z_relative, n_freq, n_kz), 1, n_freq, n_kz), n_ky, n_freq, n_kz)
        
        cos_arr = cos(freq_kz_arr)
        sin_arr = sin(freq_kz_arr)
        z_exp_arr = exp(-1*dcomplex(0,1)*freq_kz_arr)
        
        for i=0, n_kx-1 do begin
           sigma2_arr = rebin(reform(sigma2_cube[i,*,*]), n_ky, n_freq, n_kz) 
           
           covar_cos[i,*,*] = total(sigma2_arr*cos_arr^2d, 2)/double(n_kz)
           covar_sin[i,*,*] = total(sigma2_arr*sin_arr^2d, 2)/double(n_kz)
           covar_cross[i,*,*] = total(sigma2_arr*cos_arr*sin_arr, 2)/double(n_kz)
        endfor
        
        ;; cos 0 term has different normalization
        covar_cos[*,*,0] = covar_cos[*,*,0]/4d
        
        freq_kz_arr = 0
        cos_arr = 0
        sin_arr = 0
        sigma2_arr = 0

     endelse

     sim_cube = 0
     sigma2_cube = 0


    ;; factor to go to eor theory FT convention
     deltas = [x_mpc_delta, y_mpc_delta, z_mpc_delta]
     factor = (1d/(2d*!pi))^3d * product(deltas) * product([n_kx, n_ky, n_freq])
     ft = factor * temporary(ft)
     if keyword_set(std_power) then sigma2_ft = factor * temporary(sigma2_ft) $
     else begin
        ;; I don't think I need these factors in the covariance
        ;; matrix because I've use the FT & inv FT -- should cancel
        ;; covar_cos = factor * temporary(covar_cos2)
        ;; covar_sin = factor * temporary(covar_sin2)
        ;; covar_cross = factor * temporary(covar_cross)
     endelse
    
     if keyword_set(std_power) then begin
        power_3d = real_part(ft * conj(ft))
        sigma2_3d = real_part(sigma2_ft * conj(sigma2_ft))
        ft = 0
        sigma2_ft = 0

        weights_3d = 1d/sigma2_3d
        wh = where(sigma2_3d eq 0, count)
        if count gt 0 then weights_3d[wh] = 0
        sigma2_3d=0

     endif else begin
        
        ;; get rotation angle to diagonalize covariance block
        theta = atan(2*covar_cross, covar_cos - covar_sin)/2d
        cos_theta = cos(theta)
        sin_theta = sin(theta)

        theta = 0

        sigma1_2 = covar_cos*cos_theta^2 + 2d*covar_cross*cos_theta*sin_theta + covar_sin*sin_theta^2d
        sigma2_2 = covar_cos*sin_theta^2d - 2d*covar_cross*cos_theta*sin_theta + covar_sin*cos_theta^2d

        covar_cos = 0
        covar_sin = 0
        covar_cross = 0

        data_cos = real_part(ft)
        data_sin = imaginary(ft)
        ft = 0
        
        data1 = data_cos*cos_theta + data_sin*sin_theta
        data2 = (-1d)*data_cos*sin_theta + data_sin*cos_theta
        data_cos = 0
        data_sin = 0

        weights_1 = (1/sigma1_2)^2d
        term1 = data1^2d*weights_1
        wh = where(sigma1_2^2d eq 0, count)
        if count ne 0 then begin
           weights_1[wh] = 0
           term1[wh] = 0
        endif
        
        weights_2 = (1/sigma2_2)^2d
        term2 = data2^2d*weights_2
        wh = where(sigma2_2^2d eq 0, count)
        if count ne 0 then begin
           weights_2[wh] = 0
           term2[wh] = 0
        endif
        data1=0
        data2=0
        sigma1_2=0
        sigma2_2=0       

        weights_3d = weights_1 + weights_2
        power_3d = (term1 + term2) / weights_3d
        ;;power_error = 1/weights
        wh = where(weights_3d eq 0, count)
        if count ne 0 then begin
           power_3d[wh] = 0
           ;;power_error[wh] = 0
        endif
        term1=0
        term2=0
        weights1 = 0
        weights2 = 0

     endelse

     ;; save power and weights
     save, file = savefile, power_3d, weights_3d

  endif

  if keyword_set(match_datta) then edge_on_grid = 1

  kx_mpc = kx_mpc_fft
  ky_mpc = ky_mpc_fft
  n_kz = n_elements(kz_mpc)

  if keyword_set (no_kzero) then begin 
     ;; leave out kz=0 -- full of foregrounds
     kz_mpc = kz_mpc[1:*]
     power_3d = temporary(power_3d[*, *, 1:*])
     weights_3d = temporary(weights_3d[*, *, 1:*])
     n_kz = n_elements(kz_mpc)
  endif

  print, 'Binning to 2D power spectrum'

  bins_per_decade = 10d
  if keyword_set(no_weighting) then $
     power_2d = kspace_rebinning_2d(power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, binned_weights = weights_2d, $
                                    bins = bins_per_decade, edge_on_grid = edge_on_grid, match_datta = match_datta) $
  else power_2d = kspace_rebinning_2d(power_3d, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, weights = weights_3d, $
                                      binned_weights = weights_2d, bins = bins_per_decade, edge_on_grid = edge_on_grid, $
                                      match_datta = match_datta)

  power = power_2d
  weights = weights_2d
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
 
  fadd = ''
  if keyword_set(no_weighting) then fadd = fadd + '_nowt'
  if keyword_set(no_kzero) then fadd = fadd + '_nok0'
  if keyword_set(edge_on_grid) then fadd = fadd + '_edgegrid'
  if keyword_set(match_datta) then fadd = fadd + '_match'

  savefile = froot + filebase + fadd+ '_2dkpower.idlsave'
  save, file = savefile, power, weights, kperp_edges, kpar_edges, bins_per_decade
  power = 0

  if not keyword_set(quiet) then begin
     datta_2d_plots, no_kzero = no_kzero, edge_on_grid = edge_on_grid, match_datta = match_datta, psf = psf, $
                     std_power = std_power, no_weighting = no_weighting, /from_3d

     datta_2d_plots, /plot_weights, window_num = 2, no_kzero = no_kzero, edge_on_grid = edge_on_grid, $
                     match_datta = match_datta, psf = psf, std_power = std_power, no_weighting = no_weighting, /from_3d
  endif  

  print, 'Binning to 1D power spectrum'
  
  bins_per_decade = 10d
  if keyword_set(no_weighting) then $
     power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, bins = bins_per_decade, $
                                    binned_weights = weights_1d, mask = mask, pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, $
                                    k2_mask = k2_mask,  k3_mask = k3_mask, edge_on_grid = edge_on_grid, match_datta = match_datta)$
  else power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, bins = bins_per_decade, weights = weights_3d, $
                                      binned_weights = weights_1d, mask = mask, pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, $
                                      k2_mask = k2_mask,  k3_mask = k3_mask, edge_on_grid = edge_on_grid, match_datta = match_datta)

  power = power_1d
  weights = weights_1d
  k_edges = k_edges_mpc

  savefile = froot + filebase + fadd + '_1dkpower.idlsave'
  if n_elements(mask) ne 0 then save, file = savefile, power, weights, k_edges, bins_per_decade, mask $ 
  else save, file = savefile, power, weights, k_edges, bins_per_decade

  if not keyword_set(quiet) then $
     datta_1d_plots, no_kzero = no_kzero, edge_on_grid = edge_on_grid, match_datta = match_datta, window_num=3

end
