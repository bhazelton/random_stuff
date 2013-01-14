

pro reg_grid_3dps, info_file, std_power = std_power, no_weighting = no_weighting, refresh = refresh, no_kzero = no_kzero, $
                   flag_weights = flag_weights

  froot = base_path() + 'power_spectrum/chris_maps/'
  if n_elements(info_file) eq 0 then info_file = froot + 'chris_info.txt'

  info = read_info_file(info_file)
  nfiles = n_elements(info.files)

  ;; restore metadata
  test_file = file_test(info.metadata_file) *  (1 - file_test(info.metadata_file, /zero_length))  
  if test_file eq 0 then begin
     reg_grid_metadata, info_file 
  endif
  restore, info.metadata_file

  n_kx = n_elements(kx_mpc)
  n_ky = n_elements(ky_mpc)
  n_kz = n_elements(kz_mpc)

  fadd = ''
  if keyword_set(std_power) then fadd = fadd + '_sp'
  if keyword_set(no_weighting) then fadd = fadd + '_nowt'
  if keyword_set(flag_weights) then fadd = fadd + '_flagwt'

  save_file = strsplit(info.power_file, '.idlsave', /regex, /extract) + fadd + '.idlsave'
  
  test_save = file_test(save_file) *  (1 - file_test(save_file, /zero_length))

  if test_save eq 0 or keyword_set(refresh) then begin
     uvf_cube = complex(fltarr([n_kx, n_ky, n_freq]))
     sigma2_cube = complex(fltarr([n_kx, n_ky, n_freq]))

     print, 'Begin reading in images'
     for i=0l, n_freq-1 do begin
       
        if nfiles eq 1 then begin
           image = double(readfits(info.files[0], nslice = i, /silent))
           if keyword_set(flag_weights) then if max(image) eq 0 then weights = image * 0d else weights = image * 0d + 1d $
           else if info.weights_file ne '' then weights = double(readfits(info.weights_file[0], nslice = i, /silent))
        endif else if nfiles eq n_freq then begin
           image = double(readfits(info.files[i], /silent))
           if keyword_set(flag_weights) then if max(image) eq 0 then weights = image * 0d else weights = image * 0d + 1d $
           else if info.weights_file ne '' then weights = double(readfits(info.weights_file[i], /silent))
        endif
        if n_elements(size(image, /dimension)) gt 2 then image = reform(image[*,*,0])
        if n_elements(size(weights, /dimension)) gt 2 then weights = reform(weights[*,*,0])
    
        if n_elements(weights) eq 0 then weights = fltarr(n_kx, n_ky) + 1d
        image = image * mk_conv_factors[i]
        weights = weights * mk_conv_factors[i]
    
        ;; do ft in 2 spacial directions
        ft = FFT(image)           
        weight_ft = FFT(weights)           
        image = 0
        weights = 0

        ;; now need to shift the FFT to put the zero at the center (align with uvdensity)
        center = [n_kx, n_ky]/2.-1
        if (total(ceil(center)-floor(center)) ne 0) then stop $
        else center = floor(center)
        
        uvf_cube[*,*,i] = shift(temporary(ft), center)
        
        wh = where(weight_ft eq 0, count)
        sigma2 = 1/shift(temporary(weight_ft), center)
        if count gt 0 then sigma2[wh]=0

        sigma2_cube[*,*,i] = temporary(sigma2)

     endfor
     print, 'starting fft in freq. dimension'

     ;; do fft in third dimension
     ft = FFT(temporary(uvf_cube), dimension=3)
     ft = ft[*,*,0:n_kz-1]
     
     if keyword_set(std_power) then begin
        ;; for standard power calc. just need ft of sigma2
        sigma2_ft = fft(sigma2_cube, dimension=3)
        sigma2_ft = sigma2_ft[*,*,0:n_kz-1]
     endif else begin
        ;; for new power calc, need cos2, sin2, cos*sin transforms
        
        ;; have to do this in a for loop for memory's sake
        covar_cos = fltarr(n_kx, n_ky, n_kz)
        covar_sin = fltarr(n_kx, n_ky, n_kz)
        covar_cross = fltarr(n_kx, n_ky, n_kz)
        
        ;; comov_dist_los goes from large to small z
        z_relative = findgen(n_freq)*z_mpc_delta
        freq_kz_arr = rebin(reform(rebin(reform(kz_mpc, 1, n_kz), n_freq, n_kz) * $
                                   rebin(z_relative, n_freq, n_kz), 1, n_freq, n_kz), n_ky, n_freq, n_kz)
        
        cos_arr = cos(freq_kz_arr)
        sin_arr = sin(freq_kz_arr)
        z_exp_arr = exp(-1*complex(0,1)*freq_kz_arr)
        
        print, 'starting covariance calc.'
        for i=0, n_kx-1 do begin
           ;;sigma2_arr = rebin(reform(sigma2_cube[i,*,*]), n_ky, n_freq, n_kz) 
           sigma2_arr = make_array(n_ky, n_freq, n_kz,/dcomplex)
           u = 1 + dblarr(n_kz)
           image = sigma2_cube[i,*,*]
           sigma2_arr[*] = image[*] # u

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

     uvf_cube = 0
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

     print, 'starting power calculation'
     if keyword_set(std_power) then begin
        power_3d = real_part(ft * conj(ft))
        sigma2_3d = real_part(sigma2_ft * conj(sigma2_ft))
        ft = 0
        sigma2_ft = 0

        weights = 1d/sigma2_3d
        wh = where(sigma2_3d eq 0, count)
        if count gt 0 then weights[wh] = 0
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
        sim1=0
        sim2=0
        sigma1_2=0
        sigma2_2=0       

        weights = weights_1 + weights_2
        power_3d = (term1 + term2) / weights
        ;;power_error = 1/weights
        wh = where(weights eq 0, count)
        if count ne 0 then begin
           power_3d[wh] = 0
           ;;power_error[wh] = 0
        endif
        term1=0
        term2=0
        weights1 = 0
        weights2 = 0
     endelse
 
     print, 'saving power file'
     save, file = save_file, power_3d, weights

  endif else restore, save_file

  if keyword_set (no_kzero) then begin 
     ;; leave out kz=0 -- full of foregrounds
     kz_mpc = kz_mpc[1:*]
     power_3d = temporary(power_3d[*, *, 1:*])
     n_kz = n_elements(kz_mpc)
  endif

  print, 'Rebinning to 2D'
  if keyword_set(no_weighting) then power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, $
                                                                       kpar_edges_mpc, binned_weights = binned_weights, bins = 10) $
  else power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, $
                                         kpar_edges_mpc, weights = weights, binned_weights = binned_weights, bins = 10)

  dims = size(power_rebin, /dimension)
  n_kperp = dims[0]
  n_kpar = dims[1]

  savefile = repstr(info_file, '_info.txt','') + fadd + '_2dkpower' + '.idlsave'

  power = power_rebin
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  temp = weights
  weights = binned_weights
  bins_per_decade = 10

  save, file = savefile, power, weights, kperp_edges, kpar_edges, bins_per_decade
  weights = temp
stop
  kpower_2d_plots, savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                   data_range = data_range, pub = pub, plotfile = plotfile

end
