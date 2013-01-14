;; Little program to read in temperature data output from psql on the new receivers

pro read_temp_data

  froot = base_path() + 'MWA/recv_temp_data/'

  filename = froot + 'temp_out.txt'

  temp_record = {begintime:0, atim_temp:0., asc1_temp1:0., asc2_temp1:0., adfb1_temps:fltarr(8), $
                 adfb2_temps:fltarr(8), compressor:0., evaporator:0., outside:0., inside:0., condensor:0.}

  max_records = 100
  data = replicate(temp_record, 

  openr, lun, infile, /get_lun
  nrecords = 0
  temp = ''
  while eof(lun) ne 1 do begin 
     readf, lun, temp

     if strpos(temp, 'row') ne -1 then test_end = 1 else test_end = 0

     if n_records gt 0 and test_end eq 0 then begin
        col_res = strsplit(temp, ' |', /extract)

        if n_elements (col_result) ne 11 then message, 'incorrect number of columns'

        for i=0, 10 do begin
           if i eq 4 or i eq 5 then begin
              res = strsplit(col_result[i], ' {},', /extract)
              temp_record.(i) = res
stop
           endif else temp_record.(i) = col_result[i]
        endfor
stop
        data[nrecords] = temp_record
        nrecords = nrecords + 1
        if nrecords ge max_records-1 then message, 'maximum number of records reached'
     endif

  endwhile
  free_lun, lun
  data = data[0:nrecords-1]






end
