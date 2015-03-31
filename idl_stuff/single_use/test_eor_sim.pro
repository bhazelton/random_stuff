pro test_eor_sim, delta_uv = delta_uv, uv_max = uv_max, n_freq = n_freq, f_avg = f_avg, uv_avg = uv_avg, spec_window_type = spec_window_type, $
    flat_sigma = flat_sigma, no_distrib = no_distrib, use_sim = use_sim, apply_beam = apply_beam, sim_beam = sim_beam, $
    sample_factor = sample_factor, uniform_sampling = uniform_sampling, beam_size_factor = beam_size_factor, calc_al_weights = calc_al_weights, $
    save_cubefile = save_cubefile, no_plots = no_plots, sim_power=sim_power, ver_name = ver_name, al_norm2 = al_norm2
    
    
  if keyword_set(apply_beam) then begin
    if n_elements(sim_beam) eq 0 then $
      if n_elements(beam_size_factor) eq 0 then beam_size_factor=10. $
    else print, 'beam_size_factor has no effect if sim_beam is set. Ignoring.'
  endif
  
  if n_elements(save_cubefile) gt 0 and size(save_cubefile,/type) ne 7 then begin
    save_path = base_path('data') +'fhd_sim_data/snap_highf_noinst_'
    if n_elements(ver_name) gt 0 then save_path = save_path + ver_name + '_'
    if keyword_set(al_norm2) then save_path = save_path + 'norm2_'
    if keyword_set(flat_sigma) then begin
      if keyword_set(no_distrib) then save_path = save_path + 'delta_image_' $
      else save_path = save_path + 'flat_'
    endif else save_path = save_path + 'eor_nomu_'
    
    if keyword_set(apply_beam) then begin
      if n_elements(sim_beam) gt 0 then begin
        save_path = save_path + 'beaminst'
        if sim_beam gt 1 then save_path = save_path + number_formatter(sim_beam)
        save_path = save_path + '_'
      endif else save_path = save_path + 'beam' + number_formatter(beam_size_factor) + '_'
    endif
    
    if keyword_set(uniform_sampling) then save_path = save_path + 'sampleuniform_' $
    else if n_elements(sample_factor) gt 0 then save_path = save_path + 'samplefactor' + number_formatter(sample_factor) + '_'
    
    save_path = save_path + 'uvin/'
    save_cubefile = save_path + 'sim_noinst_' + ['even_xx', 'odd_xx', 'even_yy', 'odd_yy'] + '_uvf.sav'
    
    if not file_test(save_path, /directory) then file_mkdir, save_path
  endif
  
  if keyword_set(use_sim) then begin
    restore, base_path('data') +'fhd_sim_data/snap_highf_eor_nomu_beamfix/1061316176_input_model.sav' ; model_uvf, uv_arr, freq_arr
    eor_uvf_cube = temporary(model_uvf)*2.
    n_freq = n_elements(freq_arr)
    n_uv = n_elements(uv_arr)
    
    delta_uv = uv_arr[1] - uv_arr[0]
    uv_max = max(uv_arr)
    title = 'delta uv: ' + number_formatter(delta_uv) + ', max uv: ' + number_formatter(uv_max)
    
    
    if keyword_set(apply_beam) then begin
      ;      beam2_image = getvar_savefile(base_path('data') +'fhd_sim_data/snap_highf_eor_nomu_newconv/1061316176_initial_beam2_image.sav', 'beam2_xx_image')
    
      beam2_image = getvar_savefile(base_path('data') +'fhd_sim_data/snap_highf_eor_nomu_beamfix/1061316176_initial_beam2_image.sav', beam2_xx_image)
      beam=sqrt(beam2_image)
      
      temp = shift(fft(fft(shift(eor_uvf_cube,n_uv/2,n_uv/2,0), dimension=1),dimension=2),n_uv/2,n_uv/2,0)
      temp = temp * beam^2.
      temp = shift(fft(fft(shift(temp,n_uv/2,n_uv/2,0), dimension=1, /inverse), dimension=2, /inverse),n_uv/2,n_uv/2,0)
      
      eor_uvf_cube = temp
    endif
    
  endif else begin
    if keyword_set(apply_beam) and n_elements(sim_beam) gt 0 then begin
    
      case sim_beam of
        1: beam_path = base_path('data') +'fhd_sim_data/snap_highf_eor_nomu_beamfix_uvin/1061316176_gridded_beam2_image.sav'
        2: beam_path = base_path('data') +'fhd_sim_data/offzen_high_flat_uvin/1061323008_gridded_beam2_image.sav'
        else: stop
      endcase
      obs = getvar_savefile(beam_path, 'obs_out')
      
      beam_uv_pix = obs.kpix
      if n_elements(delta_uv) eq 0 then delta_uv = beam_uv_pix
      
      if n_elements(uv_max) eq 0 then uv_max = 200.
      uv_arr = findgen(uv_max*2./delta_uv)*delta_uv-uv_max
      n_uv = n_elements(uv_arr)
      
      delta_theta = 1./(n_uv * delta_uv)
      
      n_uv_beam = obs.elements
      delta_theta_beam = 1./(n_uv_beam * beam_uv_pix)
      
      n_freq = obs.n_freq
      freq_arr = (*obs.BASELINE_INFO).freq/1e6
      
      nf_vis = obs.nf_vis
      wh_no_vis = where(nf_vis eq 0, count_no_vis)
      
      beam2_image = getvar_savefile(beam_path, 'beam2_xx_image')
      beam2_image = beam2_image / rebin(reform(nf_vis, 1, 1, n_freq), n_uv_beam, n_uv_beam, n_freq, /sample)
      if count_no_vis gt 0 then beam2_image[*,*, wh_no_vis] = 0
      beam=sqrt(beam2_image)
      
      beam_total_freq = total(total(beam, 2), 1)
      wh_beam0 = where(beam_total_freq eq 0, count_beam0, complement = wh_beam_n0, ncomplement = count_beam_n0)
      if count_beam0 gt 0 then begin
        if count_beam_n0 eq 0 then message, 'beam cube is entirely 0'
        
        temp = abs(rebin(wh_beam_n0, count_beam_n0, count_beam0, /sample) - rebin(reform(wh_beam0, 1, count_beam0), count_beam_n0, count_beam0, /sample))
        temp2 = min(temp, min_locs, dimension=1)
        beam_inds_sub = wh_beam_n0[min_locs mod count_beam_n0]
        
        beam[*,*,wh_beam0] = beam[*,*,beam_inds_sub]
        
      endif
      
      if delta_uv ne beam_uv_pix then begin
        if delta_uv gt beam_uv_pix then begin
          ;; clip image beam to get correct delta_uv
          img_extent = 1./delta_uv
          n_img = img_extent/delta_theta_beam
          
          beam = beam[n_uv_beam/2 - n_img/2:n_uv_beam/2+n_img/2-1, n_uv_beam/2 - n_img/2:n_uv_beam/2+n_img/2-1, *]
          n_uv_beam = n_img
        endif else stop
        
        
      endif
      
      
      undefine_fhd, obs
    endif else begin
    
      if n_elements(n_freq) eq 0 then n_freq = 384
      freq_arr = findgen(n_freq)*0.08+120
      
      if n_elements(delta_uv) eq 0 then delta_uv = 4.
      if n_elements(uv_max) eq 0 then uv_max = 200.
      uv_arr = findgen(uv_max*2./delta_uv)*delta_uv-uv_max
      n_uv = n_elements(uv_arr)
      
      delta_theta = 1./(n_uv * delta_uv)
    endelse
    
    title = 'delta uv: ' + number_formatter(delta_uv) + ', max uv: ' + number_formatter(uv_max)
    eor_uvf_cube = eor_sim(uv_arr, uv_arr, freq_arr, flat_sigma = flat_sigma, no_distrib = no_distrib)
    
    print, 'stddev of uvf_cube: ' + number_formatter(stddev(abs(eor_uvf_cube)), format='(e10.2)')
    
    if keyword_set(apply_beam) then begin
    
      if n_elements(sim_beam) gt 0 then begin
        beam_uv = shift(fft(fft(shift(beam,n_uv_beam/2,n_uv_beam/2,0), dimension=1, /inverse), dimension=2, /inverse),n_uv_beam/2,n_uv_beam/2,0) * delta_theta_beam^2.
        
        wh_small = where(beam_uv lt max(beam_uv)*0.01, count_small)
        if count_small gt 0 then beam_uv[wh_small] = 0. else stop
        
        beam_uv = beam_uv[n_uv_beam/2-n_uv/2:n_uv_beam/2+n_uv/2-1, n_uv_beam/2-n_uv/2:n_uv_beam/2+n_uv/2-1,*]
        beam = shift(fft(fft(shift(beam_uv,n_uv/2,n_uv/2,0), dimension=1),dimension=2),n_uv/2,n_uv/2,0) * n_uv^2. * delta_uv^2.
        
      endif else begin
        theta_vals1d = (findgen(n_uv)-n_uv/2) * delta_theta
        theta_vals2d = sqrt(rebin(theta_vals1d, n_uv, n_uv, /sample)^2. + rebin(reform(theta_vals1d, 1, n_uv), n_uv, n_uv, /sample)^2.)
        sigma_beam = (delta_theta*n_uv/2.)/beam_size_factor
        gaussian_theta = exp(-1 * theta_vals2d^2./(2.*sigma_beam^2.))
        gaussian_theta = gaussian_theta / max(gaussian_theta)
        beam_tophat = fltarr(n_uv,n_uv,n_freq)+1.
        beam = rebin(gaussian_theta, n_uv, n_uv, n_freq, /sample)
        
        beam_uv = shift(fft(fft(shift(beam,n_uv/2,n_uv/2,0), dimension=1, /inverse), dimension=2, /inverse),n_uv/2,n_uv/2,0) * delta_theta^2.
        uv_dist = sqrt(rebin(uv_arr, n_uv, n_uv,/sample)^2. + rebin(reform(uv_arr, 1, n_uv), n_uv, n_uv, /sample)^2.)
        hist = histogram(uv_dist, binsize = delta_uv, locations = locs, reverse_indices=ri)
        min_per_dist = fltarr(n_elements(hist))
        beam_0 = beam_uv[*,*,0]
        for i=0, n_elements(hist)-1 do min_per_dist[i] = min(beam_0[ri[ri[i]:ri[i+1]-1]])
        
        wh_small = where(min_per_dist lt 1e-6, count_small)
        if count_small gt 0 then begin
          beam_uv = reform(beam_uv, n_uv^2., n_freq)
          for i=0, count_small-1 do beam_uv[ri[ri[wh_small[i]]:ri[wh_small[i]+1]-1],*]=0.
          beam_uv = reform(beam_uv, n_uv, n_uv, n_freq)
          beam = shift(fft(fft(shift(beam_uv,n_uv/2,n_uv/2,0), dimension=1),dimension=2),n_uv/2,n_uv/2,0) * n_uv^2. * delta_uv^2.
        endif
      endelse
      
      model_image = shift(fft(fft(shift(eor_uvf_cube,n_uv/2,n_uv/2,0), dimension=1),dimension=2),n_uv/2,n_uv/2,0) * n_uv^2. * delta_uv^2.
      
      print, 'beam fraction 1 pixel away: ', abs(beam_uv[n_uv/2,n_uv/2-1,0])/abs(beam_uv[n_uv/2,n_uv/2,0])
      
      model_image = model_image * beam
      
      if not keyword_set(no_plots) then quick_image, beam_uv[*,*,0], window=3,/log, title = 'uv beam'
      
      if n_elements(sample_factor) gt 0 or keyword_set(uniform_sampling) then begin
      
        if keyword_set(uniform_sampling) then begin
          nsample = n_uv^2.
          upix_samples = rebin(indgen(n_uv), n_uv, n_uv, /sample)
          vpix_samples = rebin(reform(indgen(n_uv), 1, n_uv), n_uv, n_uv, /sample)
        endif else begin
          nsample = round(float(n_uv^2.) * sample_factor, /L64)
          upix_samples = round(randomu(seed, nsample)*(n_uv-1), /L64)
          vpix_samples = round(randomu(seed, nsample)*(n_uv-1), /L64)
        endelse
        
        if not keyword_set(no_plots) then quick_image, eor_uvf_cube[*,*,1], window=4, title = 'uv model'
        input_model_power = total(abs(eor_uvf_cube[*,*,1])^2.) / (n_uv^2.)
        
        model_uv = shift(fft(fft(shift(model_image,n_uv/2,n_uv/2,0), dimension=1, /inverse), dimension=2, /inverse),n_uv/2,n_uv/2,0) * delta_theta^2.
        convolved_model = model_uv
        if not keyword_set(no_plots) then quick_image, convolved_model[*,*,1], window=5, title = 'uv convolved model', data_range = data_range_convol
        convolved_model_power = total(abs(convolved_model[*,*,1])^2.) / (n_uv^2.)
        
        u_kernel_range = minmax(where(total(total(beam_uv, 3),2) gt 0))
        v_kernel_range = minmax(where(total(total(beam_uv, 3),1) gt 0))
        if (u_kernel_range[1]-u_kernel_range[0]+1) ne (v_kernel_range[1]-v_kernel_range[0]+1) then begin
          n_uv_kernel = max([(u_kernel_range[1]-u_kernel_range[0]+1), (v_kernel_range[1]-v_kernel_range[0]+1)], max_loc)
          n_diff = abs((u_kernel_range[1]-u_kernel_range[0]+1) - (v_kernel_range[1]-v_kernel_range[0]+1))
          case max_loc of
            0: begin
              v_kernel_range = v_kernel_range + [-1,1]*n_diff/2
              if (u_kernel_range[1]-u_kernel_range[0]+1) ne (v_kernel_range[1]-v_kernel_range[0]+1) then stop
            end
            1: begin
              u_kernel_range = u_kernel_range + [-1,1]*n_diff/2
              if (u_kernel_range[1]-u_kernel_range[0]+1) ne (v_kernel_range[1]-v_kernel_range[0]+1) then stop
            end
          endcase
          beam_kernel = beam_uv[u_kernel_range[0]:u_kernel_range[1],v_kernel_range[0]:v_kernel_range[1],*]
        endif else begin
          beam_kernel = beam_uv[u_kernel_range[0]:u_kernel_range[1],v_kernel_range[0]:v_kernel_range[1],*]
          n_uv_kernel = n_elements(beam_kernel[*,0,0])
        endelse
        
        print, 'uv kernel width (pixels): ', n_uv_kernel
        
        model_uv_sampled=reform(model_uv*0., n_uv^2, n_freq)
        weights_uv_arr=reform(model_uv*0., n_uv^2, n_freq)
        variance_uv_arr=reform(abs(model_uv)*0., n_uv^2, n_freq)
        if keyword_set(calc_al_weights) then hmf = complex(fltarr(n_uv^2, n_uv^2, n_freq))
        beam_kernel = reform(beam_kernel, n_uv_kernel^2, n_freq)
        
        for i=0L, nsample-1 do begin
          u_ind_range = upix_samples[i] - n_uv_kernel/2 + [0, n_uv_kernel-1]
          v_ind_range = vpix_samples[i] - n_uv_kernel/2 + [0, n_uv_kernel-1]
          
          u_lim = [0, n_uv_kernel-1]
          v_lim = [0, n_uv_kernel-1]
          if u_ind_range[0] lt 0 then begin
            u_lim[0] = 0 - u_ind_range[0]
            u_ind_range[0] = 0
          endif
          if v_ind_range[0] lt 0 then begin
            v_lim[0] = 0 - v_ind_range[0]
            v_ind_range[0] = 0
          endif
          if u_ind_range[1] ge n_uv then begin
            u_lim[1] =  (n_uv_kernel-1) - (u_ind_range[1] - (n_uv - 1))
            u_ind_range[1] = n_uv-1
          endif
          if v_ind_range[1] ge n_uv then begin
            v_lim[1] =  (n_uv_kernel-1) - (v_ind_range[1] - (n_uv - 1))
            v_ind_range[1] = n_uv-1
          endif
          
          n_uinds = u_ind_range[1]-u_ind_range[0]+1
          n_vinds = v_ind_range[1]-v_ind_range[0]+1
          inds_use = rebin(lindgen(n_uinds)+u_ind_range[0], n_uinds, n_vinds, /sample) + $
            rebin(reform(lindgen(n_vinds)+v_ind_range[0], 1, n_vinds), n_uinds, n_vinds, /sample)*n_uv
          kernel_inds_use = rebin(lindgen(n_uinds)+u_lim[0], n_uinds, n_vinds, /sample) + $
            rebin(reform(lindgen(n_vinds)+v_lim[0], 1, n_vinds), n_uinds, n_vinds, /sample)*n_uv_kernel
            
          model_uv_sampled[inds_use,*] += beam_kernel[kernel_inds_use, *] * rebin_complex(model_uv[upix_samples[i], vpix_samples[i],*], n_uinds, n_vinds, n_freq)
          
          weights_uv_arr[inds_use,*] += beam_kernel[kernel_inds_use, *]
          
          variance_uv_arr[inds_use,*] += abs(beam_kernel[kernel_inds_use, *])^2.
          
          if keyword_set(calc_al_weights) then $
            for j=0, n_freq-1 do hmf[inds_use, inds_use, j] += matrix_multiply(beam_kernel[kernel_inds_use, j], conj(beam_kernel[kernel_inds_use, j]))
            
        endfor
        model_uv = reform(temporary(model_uv_sampled), n_uv, n_uv, n_freq)
        weights_uv_arr=reform(weights_uv_arr, n_uv, n_uv, n_freq)
        variance_uv_arr=reform(variance_uv_arr, n_uv, n_uv, n_freq)
        
        if keyword_set(calc_al_weights) then begin
          hmf2 = temporary(abs(hmf)^2.)
          if keyword_set(al_norm2) then al_weights = sqrt(reform(total(hmf2,1), n_uv, n_uv, n_freq))*delta_uv^2. $
          else al_weights = sqrt(reform(total(hmf2,1)*delta_uv^2., n_uv, n_uv, n_freq))
          undefine, hmf2
        endif
        
      endif else begin
        weights_uv_arr = fltarr(n_uv, n_uv, n_freq) +1.
        variance_uv_arr = weights_uv_arr
        
        ;model_image = model_image * beam
        model_uv = shift(fft(fft(shift(model_image,n_uv/2,n_uv/2,0), dimension=1, /inverse), dimension=2, /inverse),n_uv/2,n_uv/2,0) * delta_theta^2.
      endelse
      
      eor_uvf_cube = temporary(model_uv)
      beam2_image = abs(beam)^2.
      
    endif else begin
      weights_uv_arr = fltarr(n_uv, n_uv, n_freq) +1.
      variance_uv_arr = weights_uv_arr
      beam2_image =weights_uv_arr
      
    endelse
    
  endelse
  
  if n_elements(f_avg) gt 0 then if f_avg gt 1 then begin
    nf_new = floor(n_freq / f_avg)
    temp = complex(fltarr(n_uv, n_uv, nf_new))
    temp_f = fltarr(nf_new)
    for i=0, nf_new-1 do begin
      temp[*,*,i] = total(eor_uvf_cube[*,*,i*f_avg:(i+1)*f_avg-1], 3) / f_avg
      temp_f[i] = total(freq_arr[i*f_avg:(i+1)*f_avg-1]) / f_avg
    endfor
    
    title = title + ', f avg: ' + number_formatter(f_avg)
    
    eor_uvf_cube = temporary(temp)
    freq_arr = temporary(temp_f)
    n_freq = nf_new
  endif
  
  
  if n_elements(uv_avg) gt 0 then if uv_avg gt 1 then begin
    nuv_new = floor(n_uv / uv_avg)
    temp = complex(fltarr(nuv_new, n_uv, n_freq))
    temp_uv = fltarr(nuv_new)
    for i=0, nuv_new-1 do begin
      temp[i,*,*] = total(eor_uvf_cube[i*uv_avg:(i+1)*uv_avg-1,*,*], 1) / uv_avg
      temp_uv[i] = total(uv_arr[i*uv_avg:(i+1)*uv_avg-1]) / uv_avg
    endfor
    
    temp2 = complex(fltarr(nuv_new, nuv_new, n_freq))
    for i=0, nuv_new-1 do temp2[*,i,*] = total(temp[*,i*uv_avg:(i+1)*uv_avg-1,*], 2) / uv_avg
    undefine, temp
    
    title = title + ', uv avg: ' + number_formatter(uv_avg)
    
    eor_uvf_cube = temporary(temp2)
    uv_arr = temporary(temp_uv)
    n_uv = nuv_new
  endif
  
  restore, filepath('eor_power_1d.idlsave',root=rootdir('FHD'),subdir='catalog_data')
  npts_log = n_elements(k_centers)
  
  log_diff =  alog10(k_centers) - shift(alog10(k_centers), 1)
  log_diff = log_diff[1:*]
  log_binsize = log_diff[0]
  
  
  delta_u = uv_arr[1] - uv_arr[0]
  delta_v = uv_arr[1] - uv_arr[0]
  
  
  if min(freq_arr) gt 1e8 then frequencies = freq_arr / 1e6 else frequencies = freq_arr
  f_delta = frequencies[1]-frequencies[0]
  n_kz = n_elements(frequencies)
  
  z0_freq = 1420.40 ;; MHz
  ;; make sure frequencies are in MHz
  
  redshifts = z0_freq/frequencies - 1
  cosmology_measures, redshifts, comoving_dist_los = comov_dist_los
  
  comov_los_diff = comov_dist_los - shift(comov_dist_los, -1)
  comov_los_diff = comov_los_diff[0:n_elements(comov_dist_los)-2]
  z_mpc_delta = float(mean(comov_los_diff))
  z_mpc_mean = float(mean(comov_dist_los))
  
  ;; converting from Jy (in u,v,f) to mK*str (10^-26 * c^2 * 10^3/ (2*f^2*kb))
  conv_factor = float((3e8)^2 / (2. * (frequencies*1e6)^2. * 1.38065))
  
  ;; convert from Jy -> mk*str -> mK*Mpc^2
  conv_factor = conv_factor * z_mpc_mean^2.
  
  ;; convert from uv (in wavelengths) to kx/ky in inverse comoving Mpc
  kx_mpc = uv_arr * (2.*!pi) / z_mpc_mean
  kx_mpc_delta = delta_u * (2.*!pi) / z_mpc_mean
  n_kx = n_elements(kx_mpc)
  
  ky_mpc = uv_arr * (2.*!pi) / z_mpc_mean
  ky_mpc_delta = delta_v * (2.*!pi) / z_mpc_mean
  n_ky = n_elements(ky_mpc)
  
  z_mpc_length = float(max(comov_dist_los) - min(comov_dist_los) + z_mpc_delta)
  kz_mpc_range =  (2.*!pi) / (z_mpc_delta)
  kz_mpc_delta = (2.*!pi) / z_mpc_length
  kz_mpc = findgen(round(kz_mpc_range / kz_mpc_delta)) * kz_mpc_delta - kz_mpc_range/2.
  if n_elements(kz_mpc) ne n_kz then stop
  
  
  if keyword_set(calc_al_weights) then begin
    signal_mk = eor_uvf_cube/al_weights
    wh_wt0 = where(al_weights eq 0, count_wt0, complement=wh_wt_gt0)
    if count_wt0 gt 0 then signal_mk[wh_wt0] = 0
  endif else begin
    signal_mk = eor_uvf_cube/weights_uv_arr
    wh_wt0 = where(weights_uv_arr eq 0, count_wt0, complement=wh_wt_gt0)
    if count_wt0 gt 0 then signal_mk[wh_wt0] = 0
  endelse
  ;; convert from Jy -> mK Mpc^2
  for i=0, n_kz-1 do signal_mk[*,*,i] = signal_mk[*,*,i] * conv_factor[i]
  
  print, 'stddev of uvf_cube in mK: ' + number_formatter(stddev(abs(signal_mk)), format='(e10.2)')
  
  ;; apply spectral windowing function if desired
  if n_elements(spec_window_type) ne 0 then begin
    window = spectral_window(n_freq, type = spec_window_type, /periodic)
    
    norm_factor = sqrt(n_freq/total(window^2.))
    
    window = window * norm_factor
    
    window_expand = rebin(reform(window, 1, 1, n_freq), n_kx, n_ky, n_freq, /sample)
    
    signal_mk = signal_mk * window_expand
  endif
  
  signal_mk_k = fft(temporary(signal_mk), dimension = 3) * z_mpc_delta * n_kz
  ;; shift
  signal_mk_k = shift(temporary(signal_mk_k), [0,0,n_kz/2])
  
  print, 'stddev of k_cube in mK: ' + number_formatter(stddev(abs(signal_mk_k)), format='(e10.2)')
  
  k_arr_2d = sqrt(reform(rebin(kx_mpc, n_kx, n_ky, /sample)^2. + rebin(reform(ky_mpc, 1, n_ky), n_kx, n_ky, /sample)^2.,n_kx*n_ky))
  hist_kpar = histogram(k_arr_2d, binsize = kx_mpc_delta, min = 0, locations = locs_kpar, reverse_indices=ri_kpar)
  temp_power = fltarr(n_elements(hist_kpar), n_kz)
  temp_signal = reform(signal_mk_k, n_kx*n_ky, n_freq)
  for i=0, n_elements(hist_kpar)-1 do begin
    if hist_kpar[i] eq 0 then continue
    if hist_kpar[i] gt 1 then temp_power[i, *] = total(abs(temp_signal[ri_kpar[ri_kpar[i]:ri_kpar[i+1]-1], *])^2.,1)/hist_kpar[i] $
    else temp_power[i, *] = abs(temp_signal[ri_kpar[ri_kpar[i]:ri_kpar[i+1]-1], *])^2.
  endfor
  undefine, temp_signal
  kz_bin = abs(round(kz_mpc / kz_mpc_delta))
  hist_kz = histogram(kz_bin, binsize = 1, min=0, locations = locs_kz, reverse_indices=ri_kz)
  locs_kz = locs_kz*kz_mpc_delta
  power_2d = fltarr(n_elements(hist_kpar), n_elements(hist_kz))
  for i=0, n_elements(hist_kz)-1 do begin
    if hist_kz[i] eq 0 then continue
    if hist_kz[i] gt 1 then power_2d[*,i] = total(temp_power[*, ri_kz[ri_kz[i]:ri_kz[i+1]-1]],2)/hist_kz[i] else power_2d[*,i] = temp_power[*, ri_kz[ri_kz[i]:ri_kz[i+1]-1]]
  endfor
  
  undefine, temp_power
  
  k_arr = sqrt(rebin(kx_mpc, n_kx, n_ky, n_kz, /sample)^2. + rebin(reform(ky_mpc, 1, n_ky), n_kx, n_ky, n_kz, /sample)^2. + $
    rebin(reform(kz_mpc, 1, 1, n_kz), n_kx, n_ky, n_kz, /sample)^2.)
    
  ;; only include non-zero locations
  if count_wt0 then k_arr = k_arr[wh_wt_gt0]
  
  hist = histogram(alog10(k_arr), binsize = log_binsize, min = min(alog10(k_centers))-log_binsize/2., locations = locs, reverse_indices=ri)
  new_power2 = fltarr(n_elements(hist))
  for i=0, n_elements(hist)-1 do if hist[i] gt 1 then new_power2[i] = mean(abs(signal_mk_k[ri[ri[i]:ri[i+1]-1]])^2.)
  
  print, 'mean of power in mK (before window divided out): ' + number_formatter(mean(new_power2), format='(e10.2)')
  
  ;; adjust for window function
  if keyword_set(apply_beam) then begin
    xy_mpc_delta = (2.*!pi) / (n_uv * kx_mpc_delta)
    
    window_int = total(abs(beam[*,*,0])^2.*xy_mpc_delta^2.) * ((2.*!pi)/kz_mpc_delta)
    
    volume = ((2.*!pi)^3./(kx_mpc_delta*ky_mpc_delta*kz_mpc_delta))
  endif else begin
    window_int = ((2.*!pi)^3./(kx_mpc_delta*ky_mpc_delta*kz_mpc_delta))
  endelse
  new_power2 = new_power2 / window_int
  print, 'window_int', window_int
  
  print, 'mean of power in mK (final): ' + number_formatter(mean(new_power2), format='(e10.2)')
  
  if not keyword_set(no_plots) then begin
    if windowavailable(1) then begin
      wset, 1
      if !d.y_size/!d.x_size gt 2 or !d.x_size/!d.y_size gt 2 then make_win = 1 else make_win = 0
    endif else make_win = 1
    if make_win eq 1 then window, 1
    
    if keyword_set(flat_sigma) then yrange = minmax([max(power), new_power2[where(new_power2 gt 0)], 1e5, 1e7]) else yrange = minmax([power, new_power2[where(new_power2 gt 0)]])
    
    cgplot, k_centers, power, psym=-3, /ylog, /xlog, title = title, yrange = yrange
    cgplot, 10.^(locs+log_binsize/2.), new_power2, psym=6, /over, color='red'
    if keyword_set(flat_sigma) then cgplot, k_centers, (power*0d + max(power)), psym=-3, /over, color='blue'
    
    kpower_2d_plots, power_savefile, power = power_2d, kperp_edges = [locs_kpar-kx_mpc_delta/2., max(locs_kpar) + kx_mpc_delta/2.], $
      kpar_edges = [locs_kz-kz_mpc_delta/2., max(locs_kz) + kz_mpc_delta/2.], window_num=2
  endif
  
  ratio2 = new_power2/power
  
  if not keyword_set(no_plots) then quick_image, abs(weights_uv_arr[*,*,1]), window=6, title = 'uv weights'
  
  model_uv_arr = eor_uvf_cube/2.
  temp = eor_uvf_cube[*,*,1] / weights_uv_arr[*,*,1]
  wh_wt0 = where(weights_uv_arr[*,*,1] eq 0, count_wt0, ncomplement=count_gt0)
  if count_wt0 gt 0 then temp[wh_wt0] = 0
  if not keyword_set(no_plots) then quick_image, temp, window=7, title = 'uv gridded model/weights';,/log, color_profile='sym_log'
  
  gridded_divided_model_power = total(abs(temp)^2.) / (count_gt0)
  
  beam2_factor = total(beam2_image[*,*,1]) / n_uv^2.
  sim_power = gridded_divided_model_power/beam2_factor
  
  if n_elements(input_model_power) gt 0 then begin
    print, 'input model power per pixel: ' + number_formatter(input_model_power, format='(e8.2)')
    print, 'convolved model power per pixel: ' + number_formatter(convolved_model_power/beam2_factor, format='(e8.2)')
    print, 'gridded model/weights power per pixel: '+ number_formatter(gridded_divided_model_power/beam2_factor, format='(e8.2)')
  endif
  
  if not keyword_set(no_plots) then quick_image, model_uv_arr[*,*,1], window=8, title = 'uv gridded model'
  
  if keyword_set(calc_al_weights) then begin
    if not keyword_set(no_plots) then quick_image, abs(al_weights[*,*,1]), window=9, title = 'AL uv weights'
    temp_al = eor_uvf_cube[*,*,1] / al_weights[*,*,1]
    wh_alwt0 = where(al_weights[*,*,1] eq 0, count_alwt0, ncomplement=count_al_gt0)
    if count_alwt0 gt 0 then temp_al[wh_alwt0] = 0
    if not keyword_set(no_plots) then quick_image, temp_al, window=10, title = 'uv gridded model/(AL weights)';,/log, color_profile='sym_log'
  endif
  
  if n_elements(save_cubefile) ne 0 then begin
  
    vis_noise = ptr_new(fltarr(1, n_freq) + 1.)
    obs = {max_baseline:uv_max, obsra:0, obsdec:0, zenra:0, zendec:0, n_freq:n_freq, degpix:(180./!pi)*delta_theta, $
      kpix:delta_uv, dimension:n_uv, elements:n_uv, freq:freq_arr*1e6, time_res:2, $
      n_vis:n_sample*n_freq, nf_vis:fltarr(n_freq)+nsample, vis_noise:vis_noise}
      
    for i=0, n_elements(save_cubefile)-1 do begin
      if keyword_set(calc_al_weights) then if i ge 2 then begin
        weights_uv_arr = al_weights
      endif
      
      save, file=save_cubefile[i], weights_uv_arr, model_uv_arr, variance_uv_arr, beam2_image, obs
    endfor
  endif
  
end
