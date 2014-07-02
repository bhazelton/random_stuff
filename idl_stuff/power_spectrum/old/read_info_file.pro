

function read_info_file, info_file

 test_info = file_test(info_file) *  (1 - file_test(info_file, /zero_length))  
 if test_info eq 0 then message, 'Info file does not exist.'

 openr, lun, info_file, /get_lun
 while not eof(lun) do begin
       str=''
       readf, lun, str, format = '(a)'

       temp = strsplit(str, ' ,;'+STRING(9b), /extract, count=count)
       
       if strmid(str,0,1) eq ';' then begin
          case temp[0] of
             'id':  begin
                if n_elements(id) ne 0 then message, 'More than 1 id entry'
                if count ne 2 then message, 'Info file formatting error: there must be only 1 argument for id' else id = temp[1]
             end
             
             'metadata_file': begin
                if n_elements(metadata_file) then message, 'More than 1 metadata_file entry'
                if count ne 2 then message, 'Info file formatting error: there must be only 1 argument for metadata_file' $
                else metadata_file = temp[1]
             end

             'power_file': begin
                if n_elements(power_file) then message, 'More than 1 power_file entry'
                if count ne 2 then message, 'Info file formatting error: there must be only 1 argument for power_file' $
                else power_file = temp[1]
             end


             'weights_file': begin
                if n_elements(weights_file) then message, 'More than 1 weights_file entry'                
                if count ne 2 then message, 'Info file formatting error: there must be only 1 argument for weights_file' $
                else begin
                   if n_elements(weights_file) eq 0 then weights_file = temp[1] else weights_file = [weights_file, temp[1]]
                endelse
             end

             'uvf_file': begin
                if n_elements(uvf_file) then message, 'More than 1 uvf_file entry'                               
                if count ne 2 then message, 'Info file formatting error: there must be only 1 argument for uvf_file' $
                else uvf_file = temp[1]
             end

             'format': begin
                if n_elements(format) then message, 'More than 1 format entry'                               
                temp = strsplit(str, ' ,;[]'+STRING(9b), /extract, count=fmt_count)
                format = temp[1]
                wh_match = where(['fits', 'text'] eq format, count_match)
                if count_match eq 0 then message, 'unrecognized data format type'

                if format eq 'text' then begin
                   if fmt_count lt 4 then $
                      message, 'for text files, format entry should list columns (required: pixnum, data; also allowed: weight, skip)'
                
                   for i=2, fmt_count -1 do begin
                      case temp[i] of
                         'pixnum': begin
                            if n_elements(pix_ind) ne 0 then message, 'more than one pixnum column'
                            pix_ind = i-2
                         end
                         'data':begin
                            if n_elements(data_ind) ne 0 then message, 'more than one data column'
                            data_ind = i-2
                         end
                         'weight':begin
                            if n_elements(weight_ind) ne 0 then message, 'more than one weight column'
                            weight_ind = i-2
                         end
                         'skip': if n_elements(skip_ind) eq 0 then skip_ind = i-2 else skip_ind = [skip_ind, i-2]
                         else: message, 'unrecognized column type (allowed: pixnum, data, weight, skip)'
                      endcase
                   endfor
                   if n_elements(pix_ind) eq 0 then message, 'no pixnum column specified'
                   if n_elements(data_ind) eq 0 then message, 'no data column specified'
                   if n_elements(weight_ind) eq 0 then weight_ind = -1
                   if n_elements(skip_ind) eq 0 then skip_ind = -1
                   text_cols = {pix_ind:pix_ind, data_ind:data_ind, weight_ind:weight_ind, skip_ind:skip_ind}
                endif
             end

             'theta_range': begin
                if n_elements(theta_range) then message, 'More than 1 theta_range entry'                               
                temp = strsplit(str, ' ,;[]'+STRING(9b), /extract, count=theta_count)
                n_theta = theta_count - 1
                if n_theta ne 2 then message, 'Theta range must have 2 elements'
                theta_range = temp[0:1]
             end

             'phi_range': begin
                if n_elements(phi_range) then message, 'More than 1 phi_range entry'                               
                temp = strsplit(str, ' ,;[]'+STRING(9b), /extract, count=phi_count)
                n_phi = phi_count - 1
                if n_phi ne 2 then message, 'Phi range must have 2 elements'
                phi_range = temp[1:2]
             end

             'frequencies': begin
                if n_elements(frequencies) then message, 'More than 1 frequencies entry'                               
                temp = strsplit(str, ' ,;[]'+STRING(9b), /extract, count=freq_count)
                n_freq = freq_count - 1
                frequencies = dblarr(n_freq)
                for i=0, n_freq-1 do frequencies[i] = double(temp[i+1])
             end

             else: if count ne 0 then print, 'Unrecognized entry in info file: ' + str
          endcase
       endif else begin
          if count gt 1 then message, 'Info file formatting error: non-header lines must contain only 1 file name'
          if count eq 1 then if n_elements(files) eq 0 then files = temp[0] else files = [files, temp[0]]
       endelse
    endwhile
 free_lun, lun
 nfiles = n_elements(files)
 n_wt_files = n_elements(weight_file)
 if n_wt_files gt 0 and n_wt_files ne nfiles then message, 'Number of weight files must match number of data files.'

 rewrite = 0
 if n_elements(id) eq 0 then begin
    print, 'no id string given, generating one.'
    temp = strsplit(info_file, '/', /extract)
    id = strsplit(temp[n_elements(temp)-1], '_info.txt', /regex, /extract)
    rewrite = 1
 endif
 
 if n_elements(metadata_file) eq 0 then begin
    print, 'no metadata file name given, generating one.'
    metadata_file = strsplit(info_file, '_info.txt', /regex, /extract) + '_meta.idlsave'
    rewrite = 1
 endif
   
 if n_elements(power_file) eq 0 then begin
    print, 'no 3D power file name given, generating one.'
    power_file = strsplit(info_file, '_info.txt', /regex, /extract) + '_power.idlsave'
    rewrite = 1
 endif

 if n_elements(uvf_file) eq 0 then begin
    print, 'no uvf file name given, generating one.'
    uvf_file = strsplit(info_file, '_info.txt', /regex, /extract) + '_uvf.idlsave'
    rewrite = 1
 endif

 if n_elements(format) eq 0 then begin
    print, 'no format specified, assuming fits'
    format = 'fits'
    rewrite = 1
 endif

  if n_elements(frequencies) eq 0 then begin
    stop_pos = strpos(strlowcase(files), 'mhz') - 1
    if min(stop_pos) gt -1 then begin
       if total(stop_pos - min(stop_pos)) ne 0 then begin
          start_pos = intarr(nfiles)
          for i=0, nfiles -1 do start_pos[i] = strpos(files[i], '_', stop_pos[i], /reverse_search) + 1
       endif else start_pos = strpos(files, '_', stop_pos[0], /reverse_search) + 1
       frequencies = double(strmid(files, reform(start_pos, 1, nfiles), reform(stop_pos-start_pos, 1, nfiles)))
       if min(frequencies) gt 0 then if min((frequencies - shift(frequencies,1))[1:*]) gt 0 then begin
          print, 'no frequencies given but could be extracted from filename, using the extracted frequencies.'

          rewrite = 1 
       endif else frequencies = -1d
    endif else frequencies = -1d
 endif

 if n_elements(weights_file) eq 0 then weights_file = ''
 if n_elements(theta_range) eq 0 then theta_range = [0, 180]
 if n_elements(phi_range) eq 0 then phi_range = [0, 360] 

 info = {id:id, metadata_file:metadata_file, uvf_file:uvf_file, power_file:power_file, format:format, files:files, $
         frequencies:frequencies, weights_file:weights_file, theta_range:theta_range, phi_range:phi_range}
 if format eq 'text' then info = create_struct(info, 'text_cols', text_cols)

 if rewrite then begin
    
    openw, lun, info_file, /get_lun

    printf, lun, ';; id' + string(9B) + id
    printf, lun, ';; metadata_file' + string(9B) + metadata_file
    printf, lun, ';; uvf_file' + string(9B) + uvf_file
    printf, lun, ';; power_file' + string(9B) + power_file
    printf, lun, ';; format' + string(9B) + format
    if frequencies[0] gt 0 then begin
       printstr = ';; frequencies' + string(9B) + '[' 
       for i=0, nfiles-2 do printstr = printstr + strtrim(string(frequencies[i]), 2) + ', '
       printstr = printstr + strtrim(string(frequencies[nfiles-1]), 2) + ']'
       printf, lun, printstr
    endif
       printf, lun, ''

    for i=0, n_elements(files)-1 do printf, lun, files[i]

    free_lun, lun

 endif

 return, info
end
