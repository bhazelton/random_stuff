pro read_datta_data

  froot = base_path() + 'single_use/datta_data/'

  axis_record = {value:0d}
  max_records = 100

  fig_10a_kperp = dblarr(max_records)
  infile = froot + '10a_c1.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     fig_10a_kperp[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  fig_10a_kperp = fig_10a_kperp[0:nrecords - 1]
  
  fig_10a_kpar = dblarr(max_records)
  infile = froot + '10a_c2.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     fig_10a_kpar[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  fig_10a_kpar = fig_10a_kpar[0:nrecords - 1]
  
  power_record = {values:dblarr(n_elements(fig_10a_kpar))}
  fig_10a_power = dblarr(n_elements(fig_10a_kperp), n_elements(fig_10a_kpar))
  infile = froot + '10a_c3.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, power_record, format=fmt
     fig_10a_power[nrecords, *] = power_record.values
     nrecords = nrecords + 1L
  endwhile
  free_lun, lun
  if nrecords ne n_elements(fig_10a_kperp) then message, 'incorrect number of elements in figure 10a col. 3'
  
  savefile = froot + 'figure_10a.idlsave'
  res_kperp = fig_10a_kperp
  res_kpar = fig_10a_kpar
  res_power = fig_10a_power
  save, file = savefile, res_kperp, res_kpar, res_power



  fig_10b_kperp = dblarr(max_records)
  infile = froot + '10b_c1.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     fig_10b_kperp[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  fig_10b_kperp = fig_10b_kperp[0:nrecords - 1]
  
  fig_10b_kpar = dblarr(max_records)
  infile = froot + '10b_c2.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     fig_10b_kpar[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  fig_10b_kpar = fig_10b_kpar[0:nrecords - 1]
  
  power_record = {values:dblarr(n_elements(fig_10b_kpar))}
  fig_10b_power = dblarr(n_elements(fig_10b_kperp), n_elements(fig_10b_kpar))
  infile = froot + '10b_c3.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, power_record, format=fmt
     fig_10b_power[nrecords, *] = power_record.values
     nrecords = nrecords + 1L
  endwhile
  free_lun, lun
  if nrecords ne n_elements(fig_10b_kperp) then message, 'incorrect number of elements in figure 10b col. 3'
  
  savefile = froot + 'figure_10b.idlsave'
  res_kperp = fig_10b_kperp
  res_kpar = fig_10b_kpar
  res_power = fig_10b_power
  save, file = savefile, res_kperp, res_kpar, res_power

  
  

  fig_11a_kperp = dblarr(max_records)
  infile = froot + '11a_c1.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     fig_11a_kperp[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  fig_11a_kperp = fig_11a_kperp[0:nrecords - 1]
  
  fig_11a_kpar = dblarr(max_records)
  infile = froot + '11a_c2.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     fig_11a_kpar[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  fig_11a_kpar = fig_11a_kpar[0:nrecords - 1]
  
  power_record = {values:dblarr(n_elements(fig_11a_kpar))}
  fig_11a_power = dblarr(n_elements(fig_11a_kperp), n_elements(fig_11a_kpar))
  infile = froot + '11a_c3.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, power_record, format=fmt
     fig_11a_power[nrecords, *] = power_record.values
     nrecords = nrecords + 1L
  endwhile
  free_lun, lun
  if nrecords ne n_elements(fig_11a_kperp) then message, 'incorrect number of elements in figure 11a col. 3'
  
  savefile = froot + 'figure_11a.idlsave'
  res_kperp = fig_11a_kperp
  res_kpar = fig_11a_kpar
  res_power = fig_11a_power
  save, file = savefile, res_kperp, res_kpar, res_power
  


  fig_11b_kperp = dblarr(max_records)
  infile = froot + '11b_c1.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     fig_11b_kperp[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  fig_11b_kperp = fig_11b_kperp[0:nrecords - 1]

  fig_11b_kpar = dblarr(max_records)
  infile = froot + '11b_c2.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     fig_11b_kpar[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  fig_11b_kpar = fig_11b_kpar[0:nrecords - 1]
  
  power_record = {values:dblarr(n_elements(fig_11b_kpar))}
  fig_11b_power = dblarr(n_elements(fig_11b_kperp), n_elements(fig_11b_kpar))
  infile = froot + '11b_c3.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, power_record, format=fmt
     fig_11b_power[nrecords, *] = power_record.values
     nrecords = nrecords + 1L
  endwhile
  free_lun, lun
  if nrecords ne n_elements(fig_11b_kperp) then message, 'incorrect number of elements in figure 11b col. 3'
  
  savefile = froot + 'figure_11b.idlsave'
  res_kperp = fig_11b_kperp
  res_kpar = fig_11b_kpar
  res_power = fig_11b_power
  save, file = savefile, res_kperp, res_kpar, res_power
   




  eor_arr_size = [51, 51]
  eor_record = {values:dblarr(eor_arr_size[0])}
  eor_kx = dblarr(eor_arr_size)
  infile = froot + 'p21_kx.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, eor_record, format=fmt
     eor_kx[nrecords, *] = eor_record.values
     nrecords = nrecords + 1L
  endwhile
  free_lun, lun
  if nrecords ne eor_arr_size[1] then message, 'incorrect number of elements in eor kx array'
 
  eor_ky = dblarr(eor_arr_size)
  infile = froot + 'p21_ky.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, eor_record, format=fmt
     eor_ky[nrecords, *] = eor_record.values
     nrecords = nrecords + 1L
  endwhile
  free_lun, lun
  if nrecords ne eor_arr_size[1] then message, 'incorrect number of elements in eor ky array'

  eor_power_linear = dblarr(eor_arr_size)
  infile = froot + 'p21_2d.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, eor_record, format=fmt
     eor_power_linear[nrecords, *] = eor_record.values
     nrecords = nrecords + 1L
  endwhile
  free_lun, lun
  if nrecords ne eor_arr_size[1] then message, 'incorrect number of elements in eor power array'

  savefile = froot + 'eor_signal.idlsave'
  eor_kperp = reform(eor_kx[0,*])
  eor_kpar = reform(eor_ky[*, 0])
  eor_power = alog10(eor_power_linear)
  save, file = savefile, eor_kperp, eor_kpar, eor_power
 



  noise_kpar = dblarr(max_records)
  infile = froot + 'noise_kpar.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     noise_kpar[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  noise_kpar = noise_kpar[0:nrecords - 1]

  noise_kperp = dblarr(max_records)
  infile = froot + 'noise_kperp.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, axis_record, format=fmt
     noise_kperp[nrecords] = axis_record.value
     nrecords = nrecords + 1L
     if nrecords eq max_records then message, 'maximum record reached'
  endwhile
  free_lun, lun
  noise_kperp = noise_kperp[0:nrecords - 1]

  power_record = {values:dblarr(n_elements(noise_kpar))}
  noise_power_linear = dblarr(n_elements(noise_kperp), n_elements(noise_kpar))
  infile = froot + 'thermal_noise.txt'
  openr, lun, infile, /get_lun
  nrecords = 0
  while eof(lun) ne 1 do begin 
     readf, lun, power_record, format=fmt
     noise_power_linear[nrecords, *] = power_record.values
     nrecords = nrecords + 1L
  endwhile
  free_lun, lun
  if nrecords ne n_elements(noise_kperp) then message, 'incorrect number of elements in thermal_noise.txt'
  
  savefile = froot + 'thermal_noise.idlsave'
  noise_power = alog10(noise_power_linear)
  save, file = savefile, noise_kperp, noise_kpar, noise_power

end
