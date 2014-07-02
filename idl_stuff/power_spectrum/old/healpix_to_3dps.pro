
; Keywords:
;    gamma: power law index for color table


pro healpix_to_3dps, info_file, refresh_2d = refresh_2d, refresh_3dbin = refresh_3dbin, no_kzero = no_kzero, log_kpar = log_kpar, $
                     log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, log_k1d = log_k1d, k1d_bin = k1d_bin, $
                     data_range = data_range, pub = pub, fill_holes = fill_holes, quiet = quiet, plot_uvf = plot_uvf

  if keyword_set(refresh_2d) then refresh_3dbin = 1

  if n_elements(data_range) gt 0 and n_elements(data_range) ne 2 then $s
     message, 'data_range must be a 2 element vector with the min and max values for clipping the data.'

  if n_elements(fill_holes) eq 0 then fill_holes = 0

  froot = base_path() + 'power_spectrum/healpix_maps/'
  if n_elements(info_file) eq 0 then info_file = froot + '2455280.859815_info.txt' else begin
     pos = strpos(strlowcase(info_file), '/users')
     if pos eq -1 then info_file = froot + info_file
  endelse


  plotfile_path = base_path() + 'power_spectrum/plots/rts_healpix'


  info = read_info_file(info_file)
  nfiles = n_elements(info.files)

  ;; freq_chan_width = info.freq_res
  ;; if info.freq_res le 0 then message, 'frequency resolution must be included in the info file.'

  ;; restore metadata
  test_file = file_test(info.metadata_file) *  (1 - file_test(info.metadata_file, /zero_length))  
  if test_file eq 0 then begin
     print, 'Metadata file does not exist. Running healpix_get_ft_frame now'
stop
     healpix_get_ft_frame, info_file 
  endif
  restore, info.metadata_file

  n_kx = n_elements(kx_rad_vals)
  n_ky = n_elements(ky_rad_vals)

  ;; conversion to comoving Mpc
  z0_freq = 1420.40 ;; MHz
  ;; convert to MHz
  redshifts = z0_freq/frequencies - 1
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los

  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  
  ;; check whether or not the frequencies are evenly spaced.
  n_freq = n_elements(frequencies)
  freq_diff = frequencies - shift(frequencies, 1)
  freq_diff = freq_diff[1:*]
  if total(freq_diff-freq_diff[0]) ne 0 then begin
     print, 'Frequencies are not evenly spaced.  Using a DFT.'
     freq_dft = 1

     nominal_freqs = dindgen(floor(((max(frequencies)-min(frequencies))/freq_resolution))+1)*freq_resolution + min(frequencies)
     nominal_z = z0_freq/nominal_freqs - 1
     cosmology_measures, nominal_z, comoving_dist_los = nominal_comov_dist_los
     nominal_comov_diffs = nominal_comov_dist_los - shift(nominal_comov_dist_los, -1)
     nominal_comov_diffs = nominal_comov_diffs[0:n_elements(nominal_comov_diffs)-2]

     z_mpc_delta = mean(nominal_comov_diffs)
     z_mpc_mean = mean(nominal_comov_dist_los)

  endif else begin
     freq_dft = 0

     z_mpc_delta = mean(comov_los_diff)
     z_mpc_mean = mean(comov_dist_los)

     ;; we will drop the negative kz elements, so only half as many as frequencies.
     if n_freq le 2 then n_kz = n_freq else if (n_freq mod 2) eq 0 then n_kz = n_freq / 2 else n_kz = (n_freq + 1)/2

  endelse

  z_mpc_length = max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta
  kz_mpc_range =  (2d*!pi) / (2*z_mpc_delta) ;; factor of 2 b/c only keep half
  kz_mpc_delta = (2d*!pi) / z_mpc_length
  kz_mpc = dindgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta
  if n_elements(n_kz) ne 0 then begin
     if n_elements(kz_mpc) ne n_kz then stop
  endif else n_kz = n_elements(kz_mpc)

  kx_mpc = kx_rad_vals / z_mpc_mean
  kx_mpc_delta = kx_mpc[1] - kx_mpc[0]
  kx_mpc_range = n_kx * kx_mpc_delta
  x_mpc_delta = (2d*!pi)/kx_mpc_range
  x_mpc_length = (2d*!pi)/kx_mpc_delta
 
  ky_mpc = ky_rad_vals / z_mpc_mean
  ky_mpc_delta = ky_mpc[1] - ky_mpc[0]
  ky_mpc_range = n_ky * ky_mpc_delta
  y_mpc_delta = (2d*!pi)/ky_mpc_range
  y_mpc_length = (2d*!pi)/ky_mpc_delta
  
  voxel_vol_x = x_mpc_delta * y_mpc_delta* z_mpc_delta
  voxel_vol_k = kx_mpc_delta*ky_mpc_delta*kz_mpc_delta

  kperp_lambda_conv = z_mpc_mean / (2d*!pi)

  save_file = strsplit(info.power_file, '.idlsave', /regex, /extract)
  save_file = save_file + '.idlsave'
  
  uf_savefile = strsplit(info_file, '_info.txt', /regex, /extract) + '_uf_plane.idlsave'
  vf_savefile = strsplit(info_file, '_info.txt', /regex, /extract) + '_vf_plane.idlsave'
  uv_savefile = strsplit(info_file, '_info.txt', /regex, /extract) + '_uv_plane.idlsave'

  test_save = file_test(save_file) *  (1 - file_test(save_file, /zero_length))
  if test_save eq 0 or keyword_set(refresh_3dbin) then file_match = 0 $
  else begin
     new_kx_mpc = kx_mpc
     new_ky_mpc = ky_mpc
     new_kz_mpc = kz_mpc

     ;; want to restore ONLY the kx/ky/kz_mpc to check that they're the same
     sObj = obj_new('IDL_Savefile', save_file) 
     sObj -> restore, ['kx_mpc', 'ky_mpc', 'kz_mpc']

     if n_elements(kx_mpc) ne n_kx or n_elements(ky_mpc) ne n_ky or n_elements(kz_mpc) ne n_kz or total(abs(kx_mpc - new_kx_mpc)) $
        or total(abs(ky_mpc - new_ky_mpc)) or total(abs(kz_mpc - new_kz_mpc)) then begin
       file_match = 0
        kx_mpc = new_kx_mpc
        ky_mpc = new_ky_mpc
        kz_mpc = new_kz_mpc
     endif else file_match = 1

     new_kx_mpc = 0
     new_ky_mpc = 0
     new_kz_mpc = 0
  endelse 

  if keyword_set(refresh_3dbin) or file_match eq 0 then begin

     ;; First do discrete transform over x & y
     ;; check if uf_file exists
     uvf_file = strsplit(info.uvf_file, '.idlsave', /regex, /extract)
     uvf_file = uvf_file + '.idlsave'
  
     test_save = file_test(uvf_file) *  (1 - file_test(uvf_file, /zero_length))
     if test_save eq 0 or keyword_set(refresh_2d) then uvf_file_match = 0 $
     else begin
        new_kx_rad_vals = kx_rad_vals
        new_ky_rad_vals = ky_rad_vals
        new_x_rad_vals = x_rad_vals
        new_y_rad_vals = y_rad_vals      
        new_mk_conv_factors = mk_conv_factors

        ;; want to restore ONLY the kx/ky/kz_mpc to check that they're the same
        sObj = obj_new('IDL_Savefile', uvf_file) 
        sObj -> restore, ['kx_rad_vals', 'ky_rad_vals', 'x_rad_vals', 'y_rad_vals', 'mk_conv_factors']
        
        if total(abs(kx_rad_vals - new_kx_rad_vals)) or total(abs(ky_rad_vals - new_ky_rad_vals)) or $
           total(abs(x_rad_vals - new_x_rad_vals)) or total(abs(y_rad_vals - new_y_rad_vals)) or $
           total(abs(mk_conv_factors - new_mk_conv_factors)) then begin

           uvf_file_match = 0
           kx_rad_vals = new_kx_rad_vals
           ky_rad_vals = new_ky_rad_vals
           x_rad_vals = new_x_rad_vals
           y_rad_vals = new_y_rad_vals      
           mk_conv_factors = new_mk_conv_factors
        endif else uvf_file_match = 1
       
        new_kx_rad_vals = 0
        new_ky_rad_vals = 0
        new_x_rad_vals = 0
        new_y_rad_vals = 0
        new_mk_conv_factors = 0
     endelse
   
     if keyword_set(refresh_2d) or uvf_file_match eq 0 then begin
       
        ;; get images
        image_arr = dblarr(n_elements(pixels), nfiles)
        for i = 0, nfiles -1 do begin
           ;;print, 'file ' + number_formatter(i) + ' of ' + number_formatter(nfiles)
           map_file = info.files[i]

           if info.format eq 'fits' then begin
              read_fits_map, map_file, map_out, hdr, ehdr, nside=nside, ordering=ordering, coordsys=coordsys
           
              pix_nums = long(map_out[*,pix_ind])
              data = double(map_out[*,data_ind])
           endif else begin
              ncol = max([pix_ind, data_ind, weight_ind]) + 1
              fmt = ''
              for j=0, ncol-1 do begin
                 if pix_ind eq j then fmt = fmt + 'L' else if data_ind eq j then fmt = fmt + 'F' else fmt = fmt + 'X'
                 if j lt ncol-1 then fmt = fmt + ','
              endfor

              if weight_ind ge 0 then begin
                 ind_arr = [pix_ind, data_ind, weight_ind] 
                 tags = ['pix_nums', 'data', 'weights']
                 tags = tags[sort(ind_arr)]

                 readcol, map_file, format = fmt, v1, v2, v3, count = nread, nlines = nlines, /silent
                 if nread ne nlines then stop
                 temp_struct = create_struct(tags[0], v1, tags[1], v2, tags[2], v3)
                 
                 weights = temp_struct.weights
              endif else begin
                 ind_arr = [pix_ind, data_ind]
                 tags = ['pix_nums', 'data']
                 tags = tags[sort(ind_arr)]

                 readcol, map_file, format = fmt, v1, v2, count = nread, nlines = nlines, /silent
                 if nread ne nlines then stop
                 temp_struct = create_struct(tags[0], v1, tags[1], v2)
              endelse

              data = temp_struct.data
              pix_nums = temp_struct.pix_nums
              undefine, temp_struct
           endelse
           
           ;; find the indices for the included pixels
           min_val = min([pix_nums, pixels])
           wh = where(histogram(pix_nums, min=min_val, reverse_indices = ri) gt 0 AND histogram(pixels, min=min_val) gt 0, count)
           pix_inds = lonarr(count)
           for j = 0L, count-1 do pix_inds[j] = ri[ri[wh[j]]:ri[wh[j]+1]-1]
           wh = 0
           
           image_arr[*,i] = data[pix_inds] * mk_conv_factors[i]
           pix_inds = 0
        endfor

        print, 'image variance:', variance(image_arr)

        uvf_cube = discrete_ft_2D_fast(x_rad_vals, y_rad_vals, image_arr, kx_rad_vals, ky_rad_vals, $
                                          timing = ft_time)           
        print, "partially vectorized Discrete FT time: " + string(ft_time)
        print, "time per image: " + string(ft_time/double(nfiles*2))
                
        save, file = uvf_file, uvf_cube, x_rad_vals, y_rad_vals, kx_rad_vals, ky_rad_vals, mk_conv_factors, $
              kperp_lambda_conv
   
     endif else restore, uvf_file
     if n_elements(perp_ft_arr) ne 0 then begin
        uvf_cube = perp_ft_arr
        save, file = uvf_file, uvf_cube, x_rad_vals, y_rad_vals, kx_rad_vals, ky_rad_vals, mk_conv_factors, $
              kperp_lambda_conv
     endif

     ;; save some uvf slices
     uf_slice = uvf_slice(uvf_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 1, slice_inds = n_ky/2, $
                slice_savefile = uf_savefile)
     vf_slice = uvf_slice(uvf_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 0, slice_inds = n_kx/2, $
                slice_savefile = vf_savefile)
     uv_slice = uvf_slice(uvf_cube, kx_mpc, ky_mpc, frequencies, kperp_lambda_conv, slice_axis = 2, slice_inds = 0, $
                slice_savefile = uv_savefile)

     weights = dblarr(size(uvf_cube,/dimension)) + 1d

     if freq_dft eq 1 then begin
        real_z_mpc = rebin(reform(comov_dist_los, 1, 1, nfiles), n_kx, n_ky, nfiles)

        image_ft = dcomplex(dblarr(n_kx, n_ky, n_kz))
        for i=0, n_kz-1 do begin
           z_exp = exp(-1*dcomplex(0,1)*kz_mpc[i]*real_z_mpc)
           image_ft[*,*,i] = total(uvf_cube*z_exp, 3)
        endfor
        
        ;; fix normalization -- divide by # points in ft & multiply by 2pi
        ;; (factor of 2 b/c same as taking full cube & dropping negative freqs)
        ;;image_ft = temporary(image_ft) / (n_kz*2d)
        image_ft = temporary(image_ft) * voxel_vol_x / (2*!dpi)^3d

     endif else begin
        ;; Now do FFT over frequency dimension
        image_ft = FFT(uvf_cube, dimension=3) * nfiles * voxel_vol_x / (2*!dpi)^3d

        ;; drop the negative kparallels because the images are real
        image_ft = image_ft[*,*,0:n_kz-1]
     endelse

     print, 'full_ft^2d integral:', total(abs(image_ft)^2d)
     print, 'full_ft^2d integral * delta_k^6:', total(abs(image_ft)^2d)*voxel_vol_k^2d

     ;; factor to go to eor theory FT convention
     ;; deltas = [x_mpc_delta, y_mpc_delta, z_mpc_delta]
     ;;factor = (1d/(2d*!pi))^3d * product(deltas) * product([n_kx, n_ky, n_freq])
  
     ;;print, 'full_ft^2d integral (after theory factor):', total(abs(image_ft)^2d)

     save, file = save_file, image_ft, kx_mpc, ky_mpc, kz_mpc, kperp_lambda_conv

  endif else restore, save_file

  ;; square to get power
  power_3d = real_part(image_ft * conj(image_ft))
  image_ft = 0

  if keyword_set (no_kzero) then begin 
     ;; leave out kz=0 -- full of foregrounds
     kz_mpc = kz_mpc[1:*]
     power_3d = temporary(power_3d[*, *, 1:*])
     n_kz = n_elements(kz_mpc)
  endif

  print, 'power integral:', total(power_3d)
  print, 'power integral * delta_k^6:', total(power_3d)*voxel_vol_k^2d

  print, 'Binning to 2D power spectrum'

  savefile_base = strsplit(info_file, '_info.txt', /regex, /extract)
  fadd = ''
  if keyword_set(no_kzero) then fadd = fadd + '_nok0'

  fadd_2d = ''
  if keyword_set(fill_holes) then fadd_2d = fadd_2d + '_nohole'
  if keyword_set(log_kpar) then fadd_2d = fadd_2d + '_logkpar'
  if keyword_set(log_kperp) then fadd_2d = fadd_2d + '_logkperp'


  savefile = savefile_base + fadd + fadd_2d + '_2dkpower.idlsave'
 
  power_rebin = kspace_rebinning_2d(power_3D, kx_mpc, ky_mpc, kz_mpc, kperp_edges_mpc, kpar_edges_mpc, log_kpar = log_kpar, $
                                    log_kperp = log_kperp, kperp_bin = kperp_bin, kpar_bin = kpar_bin, $
                                    binned_weights = weights_2d, fill_holes = fill_holes)

  print, 'power_2d*weights integral:', total(weights_2d*power_rebin)

  power = power_rebin
  weights = weights_2d
  kperp_edges = kperp_edges_mpc
  kpar_edges = kpar_edges_mpc
  save, file = savefile, power, weights, kperp_edges, kpar_edges, kperp_bin, kpar_bin, kperp_lambda_conv
  power = 0


 
  plotfile_base = plotfile_path + info.id
  if keyword_set(no_kzero) then plotfile_base = plotfile_base + '_nok0'
  plotfile = plotfile_base + '_kspace_power'
 ;;plotfile = plotfile + '.eps'

 ;; weight_plotfile = plotfile_path + info.id + '_weight'

  if not keyword_set(quiet) then begin
     kpower_2d_plots, savefile, kperp_plot_range = kperp_plot_range, kpar_plot_range = kpar_plot_range, $
                      data_range = data_range, pub = pub, plotfile = plotfile
     
  endif

  print, 'Binning to 1D power spectrum'
  plotfile = plotfile_base + '_1d_power'

  fadd_1d = ''
  if keyword_set(log_k) then fadd_1d = fadd_1d + '_logk'
  savefile = savefile_base + fadd + fadd_1d + '_1dkpower.idlsave'

  power_1d = kspace_rebinning_1d(power_3d, kx_mpc, ky_mpc, kz_mpc, k_edges_mpc, k_bin = k1d_bin, log_k = log_k1d, $
                                 binned_weights = weights_1d)
             
  print, 'power_1d*weights integral:', total(weights_1d*power_1d)
  
  if keyword_set(log_k1d) then k_mid = 10^(alog10(k_edges_mpc[1:*]) - (k1d_bin)/2.) $
  else k_mid = k_edges_mpc[1:*] + k1d_bin/2.
  print, 'k^2d*power_1d*dk integral * 4pi *delta_k^3:', $
         total(power_1d * k_mid^2d * (k_edges_mpc - shift(k_edges_mpc, 1))[1:*])*voxel_vol_k*4d*!dpi

  power = power_1d
  weights = weights_1d
  k_edges = k_edges_mpc

  savefile = savefile_base + fadd_1d + '_1dkpower.idlsave'
  save, file = savefile, power, weights, k_edges, k_bin

  if not keyword_set(quiet) then kpower_1d_plots, savefile, window_num=2
 
  if keyword_set(plot_uvf) then begin
     uvf_slice_plot, uf_savefile, window = 3
     uvf_slice_plot, vf_savefile, window = 4
     uvf_slice_plot, uv_savefile, window = 5
  endif

end
