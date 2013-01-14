

pro datta_psf, smooth = smooth, refresh = refresh, edge_on_grid = edge_on_grid, match_datta = match_datta

  ;; first work out uvdensity parameters
  uv_fitsfile = base_path() + 'single_use/datta_data/' + 'uv_density_MWA.fits'
  uv_header = headfits(uv_fitsfile)
  uv_dims = sxpar(uv_header, 'NAXIS*')

  ;; get number of dimensions (larger than 1)
  wh_dims = where(uv_dims gt 1, n_uvdims)

  if (n_uvdims ne 3) then message, 'uvdensity data is not 3 dimensional'

  uv_dim_type = sxpar(uv_header,'CTYPE*')
  uv_dim_units = sxpar(uv_header, 'CUNIT*')
  uv_ref_pixels = sxpar(uv_header,'CRPIX*') ;; these are 1-based NOT 0-based
  uv_ref_pix_vals = sxpar(uv_header,'CRVAL*')
  uv_ref_pix_deltas = sxpar(uv_header,'CDELT*')

  ;; drop any shallow dimensions
  if (n_uvdims ne n_elements(uv_dims)) then begin
     uv_dims = uv_dims[wh_dims]
     uv_dim_type = uv_dim_type[wh_dims]
     uv_dim_units = uv_dim_units[wh_dims]
     uv_ref_pixels = uv_ref_pixels[wh_dims]
     uv_ref_pix_vals = uv_ref_pix_vals[wh_dims]
     uv_ref_pix_deltas = uv_ref_pix_deltas[wh_dims]
  endif


  ;; check the dimensional ordering
  u_dim = where(strpos(strlowcase(uv_dim_type), 'u') ne -1, count)
  if count ne 1 or u_dim[0] ne 0 then message, 'First dimension must be contain "u".'
  if strpos(strlowcase(uv_dim_units[u_dim]), 'lambda') eq -1 then message, 'u dimension does not have units of lambda'

  v_dim = where(strpos(strlowcase(uv_dim_type), 'v') ne -1, count)
  if count ne 1 or v_dim[0] ne 1 then message, 'Second dimension must be contain "v".'
  if strpos(strlowcase(uv_dim_units[v_dim]), 'lambda') eq -1 then message, 'v dimension does not have units of lambda'

  f_dim = where(strpos(strlowcase(uv_dim_units), 'hz') ne -1, count)
  if count ne 1 or f_dim[0] ne 2 then message, 'Third dimension must be in units of hz.'

  ;; conversion to comoving Mpc
  z0_freq = 1420.40 ;; MHz
  ;; use min on one element arrays to turn them into non-array values
  uv_freq_vals = dindgen(uv_dims[2])*uv_ref_pix_deltas[2] + (uv_ref_pix_vals[2] + (uv_ref_pixels[2] -1) * uv_ref_pix_deltas[2])

  ;; convert to MHz
  if (strtrim(strlowcase(uv_dim_units[2])) eq 'hz') then uv_freq_vals = uv_freq_vals / (1d*1e6) $
  else if (strtrim(strlowcase(uv_dim_units[2])) eq 'khz') then uv_freq_vals = uv_freq_vals / (1d*1e3) $
  else if (strtrim(strlowcase(uv_dim_units[2])) ne 'mhz') then message, 'Unknown units on frequency dimension (not Hz, kHz or MHz)'
  redshifts = z0_freq/uv_freq_vals - 1
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los

  ;; get uv coordinates in 1/Mpc
  u_lambda_vals = dindgen(uv_dims[0]) * abs(uv_ref_pix_deltas[0]) + $
                  (uv_ref_pix_vals[0] + (uv_ref_pixels[0] -1) * uv_ref_pix_deltas[0])
  v_lambda_vals = dindgen(uv_dims[0]) * abs(uv_ref_pix_deltas[0]) + $
                  (uv_ref_pix_vals[0] + (uv_ref_pixels[0] -1) * uv_ref_pix_deltas[0])
  u_mpc = (2d * !pi) * u_lambda_vals / mean(comov_dist_los)
  v_mpc = (2d * !pi) * v_lambda_vals / mean(comov_dist_los)
  n_u = uv_dims[0]
  n_v = uv_dims[1]

  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  
  z_mpc_delta = mean(comov_los_diff)
  z_mpc_length = max(comov_dist_los) - min(comov_dist_los)+mean(comov_los_diff)
  kz_mpc_range =  (2d*!pi) / z_mpc_delta
  kz_mpc_delta = (2d*!pi) / z_mpc_length
  n_kz = uv_dims[2]
  kz_mpc_fft = (dindgen(n_kz)-n_kz/2+1) * kz_mpc_delta
  kz_min_ind = n_kz/2-1
  kz_max_ind = kz_min_ind + n_kz/2 - 1

  savefile = strsplit(uv_fitsfile, '.fits', /regex, /extract) + '_power'
  if keyword_set(smooth) then savefile = savefile + '_smooth' + strsplit(string(smooth), /extract)
  savefile = savefile + '.idlsave'
stop 
  test_save = file_test(savefile) *  (1 - file_test(savefile, /zero_length))
  if test_save ne 0 and not keyword_set(refresh) then begin
     print, 'Save file exists, restoring now.'
     restore, savefile
     refresh = 0

     save_dims = size(power_3d, /dimensions)
     dims_expected = [n_u, n_v, n_kz/2]
     if total(save_dims-dims_expected) ne 0 then refresh = 1

  endif else begin

     ;; ;; convert to mK (from Jy/beam)
     ;; ;; first get pixel size in strad.
     ;; pix_size_str = abs(ref_pix_deltas[y_dim]) * abs(ref_pix_deltas[x_dim]) * (!pi / 180d)^2
     ;; pix_size_str = pix_size_str[0]
     ;; ;; powers of ten are from: Jy, c^2, mK, MHz^2, kB
     ;; conv_factor = (10^(double(-26+16+3-12+23)) * 9d) / (pix_size_str * 2d * uv_freq_vals^2 * 1.38)

     uvdensity = double(readfits(uv_fitsfile, /silent))
     if uv_ref_pix_deltas[0] lt 0 then uvdensity = reverse(uvdensity)
     if uv_ref_pix_deltas[1] lt 0 then uvdensity = reverse(uvdensity, 2)

     if keyword_set(smooth) then begin
        kernel_sigma = smooth
        kernel_size = ceil(5*kernel_sigma) * 2 + 1
        vals = dindgen(kernel_size) - (kernel_size-1)/2
        x_vals = rebin(vals, kernel_size, kernel_size)
        y_vals = rebin(reform(vals, 1, kernel_size), kernel_size, kernel_size)
        kernel = exp(-1 * x_vals^2d / (2d*kernel_sigma^2d) - y_vals^2d / (2d*kernel_sigma^2d))
        kernel = kernel / total(kernel)
        kernel = reform(kernel, kernel_size, kernel_size, 1)
     
        uvdensity = convol(temporary(uvdensity), kernel)
     endif

     ft = fft(temporary(uvdensity), dimension=3)

     ;; drop the negative kparallels because the images are real
     limitedfft = ft[*,*,0:n_kz/2 - 1]
     ft=0

     ;; factor to go to eor theory FT convention
     ;; deltas = [x_mpc_delta, y_mpc_delta, z_mpc_delta]
     ;; factor = (1d/(2d*!pi))^3d * product(deltas) * product([n_u, n_v, n_kz])
     ;; limitedfft = factor * temporary(limitedfft)

     power_3d = abs(temporary(limitedfft))^2d
     save, file = savefile, power_3d

  endelse

  if keyword_set(match_datta) then edge_on_grid = 1

  kz_mpc = kz_mpc_fft[kz_min_ind:kz_max_ind]
  n_kz = n_elements(kz_mpc)

  if keyword_set (no_kzero) then begin 
     ;; leave out kz=0 -- full of foregrounds
     kz_mpc = kz_mpc[1:*]
     power_3d = temporary(power_3d[*, *, 1:*])
     n_kz = n_elements(kz_mpc)
  endif

  print, 'Binning to 2D power spectrum'

  bins_per_decade = 10d
  power_2d = kspace_rebinning_2d(power_3d, u_mpc, v_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, $
                                 bins = bins_per_decade, edge_on_grid = edge_on_grid, $
                                 match_datta = match_datta, /quiet)

  power = power_2d
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  savefile = strsplit(uv_fitsfile, '.fits', /regex, /extract)
  savefile = savefile + '_2dkpower'
  if keyword_set(smooth) then savefile = savefile + '_smooth' + strsplit(string(smooth), /extract)
  if keyword_set(no_kzero) then savefile = savefile + '_nok0'
  if keyword_set(edge_on_grid) then savefile = savefile + '_edgegrid'
  if keyword_set(match_datta) then savefile = savefile + '_match'
  savefile = savefile + '.idlsave'
  save, file = savefile, power, kperp_edges, kpar_edges, bins_per_decade
  power = 0

  datta_2d_plots, /psf, smooth = smooth, no_kzero = no_kzero, edge_on_grid = edge_on_grid, match_datta = match_datta


  print, 'Binning to 1D power spectrum'
  
  bins_per_decade = 10d
  power_1d = kspace_rebinning_1d(power_3d, u_mpc, v_mpc, kz_mpc, k_edges_mpc, bins = bins_per_decade, $
                                 mask = mask, pixelwise_mask = pixelwise_mask, k1_mask = k1_mask, k2_mask = k2_mask, $
                                 k3_mask = k3_mask, edge_on_grid = edge_on_grid, match_datta = match_datta, /quiet)
  power_3d=0

  power = power_1d
  k_edges = k_edges_mpc
  savefile = strsplit(uv_fitsfile, '.fits', /regex, /extract)
  savefile = savefile + '_1dkpower'
  if keyword_set(smooth) then savefile = savefile + '_smooth' + strsplit(string(smooth), /extract)
  if n_elements(mask) ne 0 then savefile = savefile + '_masked'
  if keyword_set(no_kzero) then savefile = savefile + '_nok0'
  if keyword_set(edge_on_grid) then savefile = savefile + '_edgegrid'
  if keyword_set(match_datta) then savefile = savefile + '_match'
  savefile = savefile + '.idlsave'
  if n_elements(mask) ne 0 then save, file = savefile, power, k_edges, bins_per_decade, mask $ 
  else save, file = savefile, power, k_edges, bins_per_decade
  power = 0

end
