pro test_ps_suite, recalculate_all = recalculate_all, uv_max=uv_max

  beam_size_factors = [1,3,5,10]
  sample_factors = [.01,.1,.25,.5,1,2]
  
  sim_powers = fltarr(n_elements(beam_size_factors), n_elements(sample_factors))
  if n_elements(uv_max) eq 0 then uv_max=100
  
  t0 = systime(1)
  for i=0, n_elements(beam_size_factors)-1 do begin
    for j=0, n_elements(sample_factors)-1 do begin
      if keyword_set(recalculate_all) then begin
        print, 'generating cubes for beam' + number_formatter(beam_size_factors[i]) + ' samplefactor' + number_formatter(sample_factors[j])
        
        
        test_eor_sim, /flat_sigma, /apply_beam, sample_factor=sample_factors[j], beam_size_factor=beam_size_factors[i], $
          uv_max=uv_max, /save_cubefile, /no_plots, /calc_al_weights, sim_power = temp
          
        sim_powers[i,j] = temp
        
        folder_name = 'snap_highf_noinst_flat_beam' + number_formatter(beam_size_factors[i]) + '_samplefactor' + number_formatter(sample_factors[j]) + '_uvin'
        print, 'calculating ps for ' + folder_name
      endif
      
      hellebore_wrapper, folder_name,/sim, /no_spec_window,/png, refresh_beam=recalculate_all, refresh_info=recalculate_all
    endfor
  endfor
  t1=systime(1)
  
  run_time = t1-t0
  if run_time lt 60 then time_str = number_formatter(run_time) + ' s' $
  else if run_time lt 3600 then time_str = number_formatter(run_time/60.) + ' m' else time_str = number_formatter(run_time/3600.) + ' h'
  print, 'run time: ' + time_str
  
  if keyword_set(recalculate_all) then print, sim_powers
  
end
