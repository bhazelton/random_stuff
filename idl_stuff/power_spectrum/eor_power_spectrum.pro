pro eor_power_spectrum

  froot = base_path('data') + 'eor_data/'
  filename = 'z8_mwa_t300_f2.dat'
  n_header_lines = 7
  record = {k21:0d, p21:0d, error:0d, upper:0d, lower:0d, del_0sq:0d, d_del_0sq:0d, del_2sq:0d, d_del_2sq:0d, $
            del_4sq:0d, d_del_4sq:0d}
  max_record = 40
  data = replicate(record, max_record)
  line = 0
  nrecords = 0
  openr, lun, filename, /get_lun
  while not eof(lun) do begin
     if (line lt n_header_lines) then begin
        temp =''
        readf, lun, temp
     endif else begin
        if nrecords ge max_record then message, "Maximum number of records reached"
        readf, lun, record
        data[nrecords] = record
        nrecords = nrecords+1
     endelse
     
     line = line + 1
  endwhile
  data = data[0:nrecords-1]
  free_lun, lun

  k_centers = data.k21

  ;; data.p21 is delta^2  (delta = sqrt(power *k^3/(2pi^2)))
  ;; so power = delta^2 * (2pi^2) / k^3
  power = data.p21 * (2d * !dpi^2d) / k_centers^3
 
  savefile = froot + 'eor_power_1d.idlsave'
  save, file = savefile, power, k_centers

end
