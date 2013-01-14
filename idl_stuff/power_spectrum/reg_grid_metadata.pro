

pro reg_grid_metadata, info_file


  froot = base_path() + 'power_spectrum/chris_maps/'
  if n_elements(info_file) eq 0 then info_file = froot + 'chris_info.txt'
  
  info = read_info_file(info_file)
  
  nfiles = n_elements(info.files)
  for i=0, nfiles-1 do begin
     map_file = info.files[i]
     header = headfits(map_file)
     
     dims = sxpar(header, 'NAXIS*')
     dim_type = sxpar(header,'CTYPE*')
     dim_units = sxpar(header, 'CUNIT*')
     cr_pixels = sxpar(header,'CRPIX*') ;; these are 1-based NOT 0-based
     cr_vals = sxpar(header,'CRVAL*')
     pix_deltas = sxpar(header,'CDELT*')
     
     wh_freq = where(strtrim(strlowcase(dim_units)) eq 'hz', count)
     if count eq 1 then begin
        freq_dim = wh_freq[0]
        n_freq = dims[freq_dim]
        if n_freq eq 1 then this_freq = cr_vals[freq_dim] $
        else begin
           this_freq = (dindgen(n_freq) - (cr_pixels[freq_dim]-1)) * pix_deltas[freq_dim] + cr_vals[freq_dim]
           ;; for j=0, n_freq-1 do begin
           ;;    image = (readfits(info.files[i], nslice = j, /silent))
           ;;    if info.weights_file ne '' then weights = (readfits(info.weights_file[i], nslice = j, /silent))

           ;;    if n_elements(size(image, /dimension)) gt 2 then image = reform(image[*,*,0])
           ;;    if n_elements(size(weights, /dimension)) gt 2 then weights = reform(weights[*,*,0])

           ;;    if max(image) gt 0 and max(weights) gt 0 then if n_elements(good) eq 0 then good = j else good = [good, j]

           ;; endfor
        endelse

        if i eq 0 then begin
           frequencies = this_freq
           ;; wh_freq_good = good
        endif else begin
           ;; wh_freq_good = [wh_freq_good, good + n_elements(frequencies)]
           frequencies = [frequencies, this_freq]
        endelse
     endif else if count eq 0 then message, 'Could not find a ctype with frequency units' $
     else message, 'Found multiple ctypes with frequency units'
     
     if i eq 0 then begin
        ref_dims = dims
        ref_dim_type = dim_type
        ref_dim_units = dim_units
        ref_cr_pixels = cr_pixels
        ref_cr_vals = cr_vals
        ref_pix_deltas = pix_deltas
     endif else begin
        if total(dims ne ref_dims) ne 0 then $
           message, 'file ' + map_file + ' does not have the same dimensions as the first listed file.'
        
        if total(dim_type ne ref_dim_type) ne 0 then $
           message, 'file ' + map_file + ' does not have the same ctypes as the first listed file.'
        
        if total(dim_units ne ref_dim_units) ne 0 then $
           message, 'file ' + map_file + ' does not have the same cunits as the first listed file.'
        
        if total(abs(cr_pixels - ref_cr_pixels)) ne 0 then $
           message, 'file ' + map_file + ' does not have the same crpixs as the first listed file.'
        
        if total(abs(cr_vals[0:1] - ref_cr_vals[0:1])) ne 0 then $
           message, 'file ' + map_file + ' does not have the same crvals as the first listed file.'
        
        if total(abs(pix_deltas - ref_pix_deltas)) ne 0 then $
           message, 'file ' + map_file + ' does not have the same cdelts as the first listed file.'
     endelse
  endfor
  frequencies = frequencies[0:192]
  ;; wh = where(wh_freq_good le 99, count)
  ;; wh_freq_good = wh_freq_good[0:count] - 3

  n_freq = n_elements(frequencies)


  ;; Calculate k step size and range, which are given by spatial resolution & size of field
  
  ;; conversion to comoving Mpc
  z0_freq = 1420.40 ;; MHz
  
  ;; convert to MHz
  if (strtrim(strlowcase(ref_dim_units[wh_freq])) eq 'hz') then frequencies = frequencies / (1d*1e6) $
  else if (strtrim(strlowcase(ref_dim_units[wh_freq])) eq 'khz') then frequencies = frequencies / (1d*1e3) $
  else if (strtrim(strlowcase(ref_dim_units[wh_freq])) ne 'mhz') then $
     message, 'Unknown units on frequency dimension (not Hz, kHz or MHz)'
  
  redshifts = z0_freq / frequencies - 1
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los

  x_rad_delta = abs(ref_pix_deltas[0]) * !pi / 180d
  n_kx = ref_dims[0]
  x_rad_length = ref_dims[0] * abs(ref_pix_deltas[0]) * !pi / 180d
  x_mpc_delta = x_rad_delta * mean(comov_dist_los)
  x_mpc_length = x_rad_length * mean(comov_dist_los)
  kx_mpc_range = (2d*!pi) / x_mpc_delta
  kx_mpc_delta = (2d*!pi) / x_mpc_length
  kx_mpc = (dindgen(n_kx)-n_kx/2+1) * kx_mpc_delta

  y_rad_delta = abs(ref_pix_deltas[1]) * !pi / 180d
  n_ky = ref_dims[1]
  y_rad_length = ref_dims[1] * abs(ref_pix_deltas[1]) * !pi / 180d
  y_mpc_delta = y_rad_delta * mean(comov_dist_los)
  y_mpc_length = y_rad_length * mean(comov_dist_los)
  ky_mpc_range = (2d*!pi) / y_mpc_delta
  ky_mpc_delta = (2d*!pi) / y_mpc_length
  ky_mpc = (dindgen(n_ky)-n_ky/2+1) * ky_mpc_delta

  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  
  z_mpc_delta = mean(comov_los_diff)
  z_mpc_length = max(comov_dist_los) - min(comov_dist_los)+mean(comov_los_diff)
  kz_mpc_range =  (2d*!pi) / z_mpc_delta
  kz_mpc_delta = (2d*!pi) / z_mpc_length
  ;; we will drop the negative kz elements, so only half as many as frequencies.
  if n_freq le 2 then n_kz = n_freq else if (n_freq mod 2) eq 0 then n_kz = n_freq / 2 else n_kz = (n_freq + 1)/2
  kz_mpc = dindgen(n_kz) * kz_mpc_delta

  ;; convert to mK (from Jy/beam)
  
  ;; calculate beam diameter in radians = c/freq*max_baseline
  max_baseline = 342.497
  beam_diameter_rad = (3d * 10^8d) / (frequencies * 10^6d * max_baseline)
  beam_area_str = !pi * beam_diameter_rad^2d /4d
  
  ;; powers of ten are from: Jy, c^2, mK, MHz, kB
  mk_conv_factors = (10^(double(-26+16+3-12+23)) * 9d) / (beam_area_str * 2d * frequencies^2 * 1.38)

  save, file = info.metadata_file, kx_mpc, ky_mpc, kz_mpc, x_mpc_delta, y_mpc_delta, z_mpc_delta, mk_conv_factors, n_freq
end
