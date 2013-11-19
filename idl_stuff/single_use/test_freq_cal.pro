function calc_freq_mode2, mode, val, mask

  dims = size(mask, /dim)
  
  case n_elements(dims) of
    1: begin
      n_freq = dims[0]
      n_val = 1
    end
    
    2: begin
      n_freq = dims[0]
      n_val = dims[1]
    end
    else: stop
  end
  
  if n_elements(val) ne n_val then stop
  
  val_arr = matrix_multiply(dblarr(n_freq)+1, val)
  x_arr = (dindgen(n_freq)/(n_freq-1))*2-1
  if n_val gt 1 then x_arr = rebin(x_arr, n_freq, n_val,/sample)
  
  case mode of
    0: freq_arr = val_arr
    1: freq_arr = val_arr * x_arr
    2: freq_arr = val_arr * (3*x_arr^2. - 1)/2.
  endcase
  
  ;freq_arr = sqrt((2*mode+1)/2.)*freq_arr
  
  return, freq_arr*mask
end

function calc_gain2, gain_modes, mode_types, mode_num, mask, amp=amp, phase=phase

  n_modes = n_elements(mode_types)
  dims_gain = size(gain_modes, /dimension)
  dims_mask = size(mask, /dimension)
  
  if n_elements(dims_gain) gt 1 then if dims_gain[1] ne dims_mask[1] then message, 'second dimension of gain_modes and mask must match'
  
  if dims_gain[0] ne n_modes then message, 'first dimension of gain_modes must match length of mode_types'
  
  
  amp = dblarr(dims_mask)
  phase = dblarr(dims_mask)
  for fi=0L,n_modes-1 do begin
    if mode_types[fi] eq 'amp' then amp += calc_freq_mode2(mode_num[fi], reform(gain_modes[fi,*]), mask) $
    else phase += calc_freq_mode2(mode_num[fi], reform(gain_modes[fi,*]), mask)
  endfor
  
  if total(abs(amp)) eq 0 then amp +=1. ;; for phase-only calibration
  if total(abs(phase)) gt 0 then gain = amp * exp(complex(0,1)*phase) else gain = amp
  
  
  return, gain
end

function amoeba_single, params

  common single_data, amoeba_model, amoeba_data, amoeba_gain_B, f_arr, mode_types, mode_num, mask
  
  dims = size(mask, /dimension)
  if n_elements(dims) gt 1 then n_vis = dims[1] else n_vis = 1
  
  if n_vis gt 1 then p = rebin(params, n_elements(params), n_vis) else p=params
  gain_A = calc_gain2(p, mode_types, mode_num, mask)
  
  data_use = amoeba_data*mask
  model_use = amoeba_model*mask
  
  diff = model_use*Conj(amoeba_gain_B)*gain_A - data_use
  
  return, total(abs(diff)^2.)
  
end

function amoeba_set, params

  common set_block, model, data, gain_B, f_arr, mode_types, mode_num, mask, inter_set_vis_inds
  
  n_tiles = n_tags(mask)
  p = reform(params, n_elements(mode_num), n_tiles)
  
  chi2_set = dblarr(n_tiles)
  
  for i=0, n_tiles-1 do begin
    dims = size(mask.(i), /dimension)
    if n_elements(dims) gt 1 then n_vis = dims[1] else n_vis = 1
    
    if n_vis gt 1 then p2 = rebin(p[*,i], n_elements(mode_num), n_vis) else p2=p[*,i]
    gain_A = calc_gain2(p2, mode_types, mode_num, mask.(i))
    
    data_use = data.(i)*mask.(i)
    model_use = model.(i)*mask.(i)
    
    diff = model_use*Conj(gain_B.(i))*gain_A - data_use
    chi2_set[i] = total(abs(diff)^2.)
  endfor
  return, total(chi2_set)
  
end


function amoeba_multi, params

  common multi_data, amoeba_model_multi, amoeba_data_multi, tile_B, tile_A, vis_tile_inds, f_arr_multi, n_tiles, n_modes, n_freqs, mode_types, mode_num, mask
  
  p = reform(params, n_modes, n_tiles)
  
  gains = calc_gain2(p, mode_types, mode_num, mask)
  
  gain_A = gains[*, tile_A]
  
  gain_B = gains[*, tile_B]
  
  diff = amoeba_model_multi*gain_A*Conj(gain_B) - amoeba_data_multi*mask[*,tile_A]*mask[*,tile_B]
  
  return, total(abs(diff)^2.)
  
end

function dfpmin_grad, params

  common single_data, amoeba_model, amoeba_data, amoeba_gain_B, f_arr, mode_types, mode_num, mask
  
  dims = size(mask, /dimension)
  if n_elements(dims) gt 1 then n_vis = dims[1] else n_vis = 1
  
  if n_vis gt 1 then p = rebin(params, n_elements(params), n_vis) else p=params
  gain_A = calc_gain2(p, mode_types, mode_num, mask, phase=phase_A)
  
  data_use = amoeba_data*mask
  model_use = amoeba_model*mask
  
  amp_grad_factor = 2.*abs(model_use)^2.*abs(amoeba_gain_B)^2.*abs(gain_A) - 2*Real_part(model_use*Conj(amoeba_gain_B)*Conj(data_use)*exp(complex(0,1)*phase_A))
  phase_grad_factor = 2*Imaginary(model_use*gain_A*Conj(amoeba_gain_B)*Conj(data_use))
  
  grad = params*0
  for i=0, n_elements(mode_types)-1 do begin
    if mode_types[i] eq 'amp' then grad[i] = total(amp_grad_factor*f_arr^(mode_num[i])) $
    else grad[i] = total(phase_grad_factor*f_arr^(mode_num[i]))
    
  endfor
  
  return, grad
  
end

function dfpmin_grad_set, params

  common set_block, model, data, gain_B, f_arr, mode_types, mode_num, mask, inter_set_vis_inds
  
  n_tiles = n_tags(mask)
  p = reform(params, n_elements(mode_num), n_tiles)
  
  grad_set = dblarr(n_elements(mode_num), n_tiles)
  
  for i=0, n_tiles-1 do begin
    dims = size(mask.(i), /dimension)
    if n_elements(dims) gt 1 then n_vis = dims[1] else n_vis = 1
    
    if n_vis gt 1 then p2 = rebin(p[*,i], n_elements(mode_num), n_vis) else p2=p[*,i]
    gain_A = calc_gain2(p2, mode_types, mode_num, mask.(i), phase=phase_A)
    
    data_use = data.(i)*mask.(i)
    model_use = model.(i)*mask.(i)
    
    amp_grad_factor = 2.*abs(model_use)^2.*abs(gain_B.(i))^2.*abs(gain_A) - 2*Real_part(model_use*Conj(gain_B.(i))*Conj(data_use)*exp(complex(0,1)*phase_A))
    phase_grad_factor = 2*Imaginary(model_use*gain_A*Conj(gain_B.(i))*Conj(data_use))
    
    if n_tiles gt 1 then begin
      if tag_exist(inter_set_vis_inds, 'ts'+number_formatter(i)) then begin
        amp_grad_factor[*, inter_set_vis_inds.(i)] = 2*amp_grad_factor[*, inter_set_vis_inds.(i)]
        phase_grad_factor[*, inter_set_vis_inds.(i)] = 2*phase_grad_factor[*, inter_set_vis_inds.(i)]
      endif
    endif
    
    for j=0, n_elements(mode_types)-1 do begin
      if mode_types[j] eq 'amp' then grad_set[j, i] = total(amp_grad_factor*f_arr.(i)^(mode_num[j])) $
      else grad_set[j, i] = total(phase_grad_factor*f_arr.(i)^(mode_num[j]))
    endfor
  endfor
  
  
  return, reform(grad_set, n_elements(mode_num)*n_tiles)
  
end


function dfpmin_grad_multi, params

  common multi_data, amoeba_model_multi, amoeba_data_multi, tile_B, tile_A, vis_tile_inds, f_arr_multi, n_tiles, n_modes, n_freqs, mode_types, mode_num, mask
  
  
  p = reform(params, n_modes, n_tiles)
  gains = calc_gain2(p, mode_types, mode_num, mask, phase=phase)
  
  gain_A = gains[*, tile_A]
  gain_B = gains[*, tile_B]
  phase_A = phase[*, tile_A]
  phase_B = phase[*, tile_B]
  
  data_use = amoeba_data_multi*mask[*,tile_A]*mask[*,tile_B]
  model_use = amoeba_model_multi*mask[*,tile_A]*mask[*,tile_B]
  
  grad = p*0
  
  for i=0, n_tiles-1 do begin
    inds = vis_tile_inds.(i)
    
    amp_grad_factor = 2.*abs(model_use[*, inds])^2.*abs(gain_B[*, inds])^2.*abs(gain_A[*, inds]) - 2*Real_part(model_use[*, inds]*Conj(gain_B[*, inds])*Conj(data_use[*, inds])*exp(complex(0,1)*phase_A[*, inds]))
    phase_grad_factor = 2*Imaginary(model_use[*, inds]*gain_A[*, inds]*Conj(gain_B[*, inds])*Conj(data_use[*, inds]))
    
    for j=0, n_elements(mode_types)-1 do begin
      if mode_types[j] eq 'amp' then grad[j,i] = total(amp_grad_factor*f_arr_multi[*, inds]^(mode_num[j])) $
      else grad[j,i] = total(phase_grad_factor*f_arr_multi[*, inds]^(mode_num[j]))
    endfor
    
  endfor
  
  return, grad
  
end

pro test_freq_cal, amoeba = amoeba, multi=multi, tile_sets=tile_sets

  n_tile = 25
  n_freq = 384
  
  tile_A_i = rebin(indgen(n_tile), n_tile, n_tile, /sample)
  tile_B_i = rebin(reform(indgen(n_tile), 1, n_tile), n_tile, n_tile, /sample)
  
  wh_keep = where(tile_A_i gt tile_B_i, n_vis)
  tile_A_i = tile_A_i[wh_keep]
  tile_B_i = tile_B_i[wh_keep]
  A_ind=[tile_A_i, tile_B_i]
  B_ind=[tile_B_i, tile_A_i]
  
  ref_tile = 1
  
  phase_modes = 2
  amp_modes = 3
  n_mode = phase_modes + amp_modes
  if phase_modes gt 0 then begin
    if amp_modes gt 0 then begin
      mode_type = [strarr(phase_modes) + 'phase', strarr(amp_modes) + 'amp']
      mode_num = [indgen(phase_modes), indgen(amp_modes)]
    endif else begin
      mode_type = strarr(phase_modes) + 'phase'
      mode_num = indgen(phase_modes)
    endelse
  endif else begin
    mode_type = strarr(amp_modes) + 'amp'
    mode_num = indgen(amp_modes)
  endelse
  
  mask = fltarr(n_freq, n_tile) + 1
  
  ;seed=100
  
  ;  true_gain_mode = dblarr(n_mode, n_tile)
  ;  if amp_modes gt 0 then begin
  ;    true_gain_mode[where(mode_type eq 'amp' and mode_num eq 0), *] = abs(randomn(seed, 1, n_tile)*.25 + 1)
  ;    if amp_modes gt 1 then true_gain_mode[where(mode_type eq 'amp' and mode_num gt 0), *] = randomn(seed, amp_modes-1, n_tile)*.1
  ;  endif
  ;
  ;  if phase_modes gt 0 then begin
  ;    true_gain_mode[where(mode_type eq 'phase'), *] = randomn(seed, phase_modes, n_tile)*!pi/8
  ;    true_gain_mode[where(mode_type eq 'phase'), ref_tile] = 0
  ;  endif
  ;  noise_amp = .05
  ;  noise = randomn(seed, n_freq, n_tile)*noise_amp + complex(0,1)*randomn(seed, n_freq, n_tile)*noise_amp
  ;  true_gain0 = calc_gain2(true_gain_mode, mode_type, mode_num, mask)
  ;  true_gain = true_gain0 + noise*mask
  
  file_path = base_path('data') + 'fhd_ps_data/1061316296/fhd_bjh_1/1061316296'
  cal_file = file_path + '_cal.sav'
  cal_init=getvar_savefile(cal_file, 'cal')
  
  obs_file = file_path + '_obs.sav'
  obs_init = getvar_savefile(obs_file, 'obs')
  freq_use = (*obs_init.baseline_info).freq_use
  
  true_gain=(*cal_init.gain[0])[*, 0:n_tile]
  
  mask = rebin(freq_use, n_freq, n_tile)
  
  temp = dblarr(n_mode, n_tile)
  f_arr1 = (dindgen(n_freq)/(n_freq-1))*2-1
  for tile_i=0L, n_tile-1 do begin
    wh_f_use = where(mask[*,tile_i] gt 0, count_f_use)
    if count_f_use eq 0 then continue
    temp[where(mode_type eq 'amp'),tile_i] = poly_fit(f_arr1[wh_f_use], abs(true_gain[[wh_f_use],tile_i]), amp_modes-1)
    temp[where(mode_type eq 'phase'),tile_i] = poly_fit(f_arr1[wh_f_use], atan(true_gain[[wh_f_use],tile_i],/phase), phase_modes-1)
  end
  fit_gain_mode = temp
  
  ;; get legendre polynomial coefficients
  if amp_modes gt 1 then begin
    fit_gain_mode[where(mode_type eq 'amp' and mode_num eq 2),*] = temp[where(mode_type eq 'amp' and mode_num eq 2),*] * 2/3
    fit_gain_mode[where(mode_type eq 'amp' and mode_num eq 0),*] = temp[where(mode_type eq 'amp' and mode_num eq 2),*] * 1/3 + $
      temp[where(mode_type eq 'amp' and mode_num eq 0),*]
  endif
  fit_gain = calc_gain2(fit_gain_mode, mode_type, mode_num, mask)
  
  ;print, true_gain_mode
  ;print, max(abs(true_gain_mode-fit_gain_mode))
  
  ;vis_model1 = dblarr(n_freq, n_vis) + 1
  vis_model1 = exp(complex(0,1)*randomn(seed, n_freq, n_vis)*!pi)
  vis_model2=[[vis_model1], [Conj(vis_model1)]]
  
  noise_amp1 = .05
  noise1 = randomn(seed, n_freq, n_vis)*noise_amp1 + complex(0,1)*randomn(seed, n_freq, n_vis)*noise_amp1
  if phase_modes gt 0 then vis_data = complex(dblarr(n_freq, n_vis)) else vis_data = dblarr(n_freq, n_vis)
  for i=0, n_vis-1 do vis_data[*,i] = true_gain[*, tile_A_i[i]] * Conj(true_gain[*, tile_B_i[i]]) * vis_model1[*,i]
  vis_data1 = vis_data;+noise1
  vis_data2 = [[vis_data1], [Conj(vis_data1)]]
  
  gain_arr_mode = dblarr(n_mode, n_tile)
  if amp_modes gt 0 then gain_arr_mode[where(mode_type eq 'amp' and mode_num eq 0),*] = 1.
  
  gain_curr_mode=gain_arr_mode
  ;gain_curr_mode = fit_gain_mode
  gain_curr = calc_gain2(gain_curr_mode, mode_type, mode_num, mask)
  
  
  if keyword_set(multi) then scale = fltarr(n_mode, n_tile) else scale = fltarr(n_mode)
  if amp_modes gt 0 then begin
    scale[where(mode_type eq 'amp' and mode_num eq 0), *] = 1
    if amp_modes gt 1 then scale[where(mode_type eq 'amp' and mode_num gt 0), *] = .5
  endif
  if phase_modes gt 0 then scale[where(mode_type eq 'phase'), *] = !pi
  
  
  max_cal_iter=30
  gain_track = dblarr(n_mode, n_tile, max_cal_iter)
  n_vis_use = lonarr(n_tile)
  
  hist_A = histogram(A_ind, min=0, max = n_tile-1, reverse_indices = ri_A)
  
  for i=0, n_tile-1 do begin
    inds = ri_A[ri_A[i]:ri_A[i+1]-1]
    
    if i eq 0 then vis_tile_inds = create_struct('t'+number_formatter(i), inds) $
    else vis_tile_inds = create_struct(vis_tile_inds, 't'+number_formatter(i), inds)
  endfor
  
  
  if keyword_set(multi) then begin
  
    common multi_data, amoeba_model_multi, amoeba_data_multi, am_tile_B, am_tile_A, am_vis_tile_inds, f_arr_multi, am_n_tile, am_n_mode, am_n_freq, am_mode_types, am_mode_num, am_mask
    f_arr_multi = rebin((dindgen(n_freq)/(n_freq-1))*2-1, n_freq, n_vis*2)
    am_tile_B = B_ind
    am_tile_A = A_ind
    am_n_tile = n_tile
    am_n_mode = n_mode
    am_n_freq = n_freq
    am_mode_types = mode_type
    am_mode_num = mode_num
    am_vis_tile_inds = vis_tile_inds
    
    am_mask = mask
    
    amoeba_model_multi = vis_model2
    amoeba_data_multi = vis_data2
    
    if not keyword_set(amoeba) then begin
    
      p=reform(gain_arr_mode, n_mode*n_tile)
      time0 = systime(1)
      dfpmin, p, 1.0e-5, fval, 'amoeba_multi', 'dfpmin_grad_multi', iter=ncalls, stepmax = max(scale)
      time1 = systime(1)
      gain_new_mode=p
      
      print, 'dfp time:', time1-time0
      
    endif else begin
    
      time0 = systime(1)
      gain_new_mode = amoeba(1.0e-6, function_name = 'amoeba_multi', P0 = reform(gain_arr_mode, n_mode*n_tile), scale = scale, ncalls=ncalls, function_value=fval)
      time1 = systime(1)
      
      print, 'amoeba time:', time1-time0
    endelse
    
    if n_elements(gain_new_mode) eq 1 then stop
    
    gain_curr_mode = reform(gain_new_mode, n_mode, n_tile)
    if phase_modes gt 0 then begin
      gain_curr_mode[where(mode_type eq 'phase'), *] -= rebin(gain_curr_mode[where(mode_type eq 'phase'), ref_tile], phase_modes, n_tile)
      
      gain_curr_mode[where(mode_type eq 'phase' and mode_num eq 0),*] = gain_curr_mode[where(mode_type eq 'phase' and mode_num eq 0),*] mod !pi
    endif
    
    gain_curr = calc_gain2(gain_curr_mode, mode_type, mode_num, mask)
    stop
    gain_arr_mode=gain_curr_mode
    gain_arr = gain_curr
    
  endif else if keyword_set(tile_sets) then begin
  
    common set_block, set_model, set_data, set_gain_B, f_arr_set, set_mode_types, set_mode_num, set_mask, inter_set_vis_inds
    
    tiles_per_set = 2
    n_sets = ceil(n_tile/float(tiles_per_set))
    set_sizes = intarr(n_sets) + tiles_per_set
    if n_sets gt 1 then set_sizes[n_sets-1] = n_tile - total(set_sizes[0:n_sets-2])
    set_edges = [0, total(set_sizes, /cumulative)]
    
    ncalls_track = lonarr(n_sets, max_cal_iter)
    
    time0 = systime(1)
    FOR i=0L,(max_cal_iter-1)>1 DO BEGIN
      phase_fit_iter=Floor(max_cal_iter/3.)
      
      gain_new_mode=dblarr(n_mode, n_tile)
      gain_track[*, *, i] = gain_curr_mode
      
      if phase_fit_iter-i GT 0 then begin
        set_mode_types = mode_type[where(mode_type eq 'phase'), *]
        set_mode_num = mode_num[where(mode_type eq 'phase'), *]
      endif else begin
        set_mode_types = mode_type
        set_mode_num = mode_num
      endelse
      
      FOR set_i=0L,n_sets-1 DO begin
        n_vis_set = intarr(set_sizes[set_i])
        
        for ti=0, set_sizes[set_i]-1 do begin
        
          if ti eq 0 then begin
            inds_set = create_struct('ts'+number_formatter(ti), vis_tile_inds.(set_edges[set_i]))
            set_model = create_struct('ts'+number_formatter(ti), vis_model2[*,inds_set.(ti)])
            set_data = create_struct('ts'+number_formatter(ti), vis_data2[*, inds_set.(ti)])
            set_gain_B = create_struct('ts'+number_formatter(ti), gain_curr[*,B_ind[inds_set.(ti)]])
            set_mask = create_struct('ts'+number_formatter(ti), mask[*,A_ind[inds_set.(ti)]])
            f_arr_set = create_struct('ts'+number_formatter(ti), rebin((dindgen(n_freq)/(n_freq-1))*2-1, n_freq, n_elements(inds_set.(ti))))
          endif else begin
            inds_set = create_struct(inds_set, 'ts'+number_formatter(ti), vis_tile_inds.(set_edges[set_i]))
            set_model = create_struct(set_model, 'ts'+number_formatter(ti), vis_model2[*,inds_set.(ti)])
            set_data = create_struct(set_data, 'ts'+number_formatter(ti), vis_data2[*, inds_set.(ti)])
            set_gain_B = create_struct(set_gain_B, 'ts'+number_formatter(ti), gain_curr[*,B_ind[inds_set.(ti)]])
            set_mask = create_struct(set_mask, 'ts'+number_formatter(ti), mask[*,A_ind[inds_set.(ti)]])
            f_arr_set = create_struct(f_arr_set, 'ts'+number_formatter(ti), rebin((dindgen(n_freq)/(n_freq-1))*2-1, n_freq, n_elements(inds_set.(ti))))
          endelse
          n_vis_set[ti] = n_elements(inds_set.(ti))
        endfor
        
        if i eq 0 then n_vis_use[set_edges[set_i]:set_edges[set_i+1]-1] = n_vis_set
        
        ;f_arr_set = rebin((dindgen(n_freq)/(n_freq-1))*2-1, [n_freq, n_vis_set])
        
        temp = lonarr(set_sizes[set_i], set_sizes[set_i])
        for ti=0, set_sizes[set_i]-1 do begin
          temp[ti, ti]=-1
          if set_sizes[set_i] eq 1 then continue
          for ti2=ti+1, set_sizes[set_i]-1 do begin
            temp[ti, ti2] = where(B_ind[inds_set.(ti)] eq set_edges[set_i]+ti2)
            
            temp[ti2, ti]=temp[ti, ti2]
          endfor
          wh_inter = where(temp[ti,*] gt -1, count_inter)
          if count_inter eq 0 then continue
          if ti eq 0 then inter_set_vis_inds = create_struct('ts'+number_formatter(ti), temp[ti, wh_inter]) $
          else inter_set_vis_inds = create_struct(inter_set_vis_inds, 'ts'+number_formatter(ti), temp[ti, wh_inter])
        endfor
        
        
        if not keyword_set(amoeba) then begin
        
          if phase_fit_iter-i GT 0 then p=reform(gain_curr_mode[where(mode_type eq 'phase'), set_edges[set_i]:set_edges[set_i+1]-1], phase_modes*set_sizes[set_i]) $
          else p=reform(gain_curr_mode[*, set_edges[set_i]:set_edges[set_i+1]-1], n_mode*set_sizes[set_i])
          
          dfpmin, p, 1.0e-3, fval, 'amoeba_set', 'dfpmin_grad_set', iter=ncalls, itmax=1000, stepmax = !pi/180.
          if phase_fit_iter-i GT 0 then begin
            gain_new_mode[where(mode_type eq 'phase'), set_edges[set_i]:set_edges[set_i+1]-1]=p
            gain_new_mode[where(mode_type eq 'amp'), set_edges[set_i]:set_edges[set_i+1]-1] = gain_curr_mode[where(mode_type eq 'amp'), set_edges[set_i]:set_edges[set_i+1]-1]
          endif else gain_new_mode[*, set_edges[set_i]:set_edges[set_i+1]-1]=p
          
        endif else begin
          stop
          gain_new_mode[*, tile_i] = amoeba(1.0e-5, function_name = 'amoeba_single', P0 = gain_curr_mode[*, tile_i], $
            scale = scale, ncalls=ncalls, function_value=fval)
        endelse
        ncalls_track[set_i, i] = ncalls
      endfor
      
      gain_old_mode=gain_curr_mode
      gain_curr_mode=(gain_new_mode+gain_old_mode)/2.
      ;gain_curr_mode=gain_new_mode
      
      if amp_modes gt 0 then begin
        IF phase_fit_iter-i GT 0  then gain_curr_mode[where(mode_type eq 'amp'),*] = gain_old_mode[where(mode_type eq 'amp'),*]
        
        if min(gain_curr_mode[where(mode_type eq 'amp' and mode_num eq 0),*]) lt 0 then $
          gain_curr_mode[where(mode_type eq 'amp' and mode_num eq 0),*] = abs(gain_curr_mode[where(mode_type eq 'amp' and mode_num eq 0),*])
      endif
      
      if phase_modes gt 0 then begin
        gain_curr_mode[where(mode_type eq 'phase'), *] -= rebin(gain_curr_mode[where(mode_type eq 'phase'), ref_tile], phase_modes, n_tile)
        gain_curr_mode[where(mode_type eq 'phase' and mode_num eq 0),*] = gain_curr_mode[where(mode_type eq 'phase' and mode_num eq 0),*] mod !pi
      endif
      
      gain_old = gain_curr
      gain_curr = calc_gain2(gain_curr_mode, mode_type, mode_num, mask)
      
    ;stop
    endfor
    time1 = systime(1)
    if not keyword_set(amoeba) then print, 'dfp loop time: ', time1-time0 else print, 'amoeba loop time: ', time1-time0
    
    if phase_modes gt 0 then gain_curr_mode[where(mode_type eq 'phase'), *] -= rebin(gain_curr_mode[where(mode_type eq 'phase'), ref_tile], phase_modes, n_tile)
    gain_curr = calc_gain2(gain_curr_mode, mode_type, mode_num, mask)
    
    ;print, max(abs(gain_curr_mode-true_gain_mode))
    print, max(abs(gain_curr_mode-fit_gain_mode))
    stop
    gain_arr_mode=gain_curr_mode
    
    gain_arr = gain_curr
    
    
    
    
  endif else begin
  
    common single_data, amoeba_model, amoeba_data, amoeba_gain_B, f_arr, as_mode_types, as_mode_num, as_mask
    
    time0 = systime(1)
    FOR i=0L,(max_cal_iter-1)>1 DO BEGIN
      phase_fit_iter=Floor(max_cal_iter/3.)
      
      gain_new_mode=dblarr(n_mode, n_tile)
      gain_track[*, *, i] = gain_curr_mode
      
      if phase_fit_iter-i GT 0 then begin
        as_mode_types = mode_type[where(mode_type eq 'phase'), *]
        as_mode_num = mode_num[where(mode_type eq 'phase'), *]
      endif else begin
        as_mode_types = mode_type
        as_mode_num = mode_num
      endelse
      
      ncalls_track = lonarr(n_tile, max_cal_iter)
      
      FOR tile_i=0L,n_tile-1 DO begin
        inds_A = vis_tile_inds.(tile_i)
        n_vis_A = n_elements(inds_A)
        
        if i eq 0 then n_vis_use[tile_i] = n_vis_A
        
        f_arr = rebin((dindgen(n_freq)/(n_freq-1))*2-1, n_freq, n_vis_A)
        amoeba_model = vis_model2[*,inds_A]
        amoeba_data = vis_data2[*, inds_A]
        amoeba_gain_B = gain_curr[*,B_ind[inds_A]]
        as_mask = mask[*,A_ind[inds_A]]
        
        if not keyword_set(amoeba) then begin
        
          if phase_fit_iter-i GT 0 then p=gain_curr_mode[where(mode_type eq 'phase'), tile_i] else p=gain_curr_mode[*, tile_i]
          
          dfpmin, p, 1.0e-3, fval, 'amoeba_single', 'dfpmin_grad', iter=ncalls, itmax=1000, stepmax = !pi/180.
          if phase_fit_iter-i GT 0 then begin
            gain_new_mode[where(mode_type eq 'phase'), tile_i]=p
            gain_new_mode[where(mode_type eq 'amp'),*] = gain_curr_mode[where(mode_type eq 'amp'),*]
          endif else gain_new_mode[*, tile_i]=p
          
        endif else begin
          gain_new_mode[*, tile_i] = amoeba(1.0e-5, function_name = 'amoeba_single', P0 = gain_curr_mode[*, tile_i], $
            scale = scale, ncalls=ncalls, function_value=fval)
        endelse
        ncalls_track[tile_i, i] = ncalls
      endfor
      
      gain_old_mode=gain_curr_mode
      gain_curr_mode=(gain_new_mode+gain_old_mode)/2.
      ;gain_curr_mode=gain_new_mode
      
      if amp_modes gt 0 then begin
        IF phase_fit_iter-i GT 0  then gain_curr_mode[where(mode_type eq 'amp'),*] = gain_old_mode[where(mode_type eq 'amp'),*]
        
        if min(gain_curr_mode[where(mode_type eq 'amp' and mode_num eq 0),*]) lt 0 then $
          gain_curr_mode[where(mode_type eq 'amp' and mode_num eq 0),*] = abs(gain_curr_mode[where(mode_type eq 'amp' and mode_num eq 0),*])
      endif
      
      if phase_modes gt 0 then begin
        gain_curr_mode[where(mode_type eq 'phase'), *] -= rebin(gain_curr_mode[where(mode_type eq 'phase'), ref_tile], phase_modes, n_tile)
        gain_curr_mode[where(mode_type eq 'phase' and mode_num eq 0),*] = gain_curr_mode[where(mode_type eq 'phase' and mode_num eq 0),*] mod !pi
      endif
      
      gain_old = gain_curr
      gain_curr = calc_gain2(gain_curr_mode, mode_type, mode_num, mask)
      
    ;stop
    endfor
    time1 = systime(1)
    if not keyword_set(amoeba) then print, 'dfp loop time: ', time1-time0 else print, 'amoeba loop time: ', time1-time0
    
    if phase_modes gt 0 then gain_curr_mode[where(mode_type eq 'phase'), *] -= rebin(gain_curr_mode[where(mode_type eq 'phase'), ref_tile], phase_modes, n_tile)
    gain_curr = calc_gain2(gain_curr_mode, mode_type, mode_num, mask)
    
    ;print, max(abs(gain_curr_mode-true_gain_mode))
    print, max(abs(gain_curr_mode-fit_gain_mode))
    stop
    gain_arr_mode=gain_curr_mode
    
    gain_arr = gain_curr
  endelse
  
  
end