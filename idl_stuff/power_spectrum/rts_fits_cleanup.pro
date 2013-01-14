

pro rts_fits_cleanup, fits_files, quiet = quiet

  nfiles = n_elements(fits_files)

  for i=0, nfiles-1 do begin
     fits_file = fits_files[i]

     hdr = headfits(fits_file, exten=1)
     
     
     ttypes = sxpar(hdr, 'TTYPE*')
     tforms = sxpar(hdr, 'TFORM*')
     
     res = strmatch(ttypes, '*weight*', /fold_case)
     
     first_chars = strmid(tforms, 0, 1)
     bad_first = strmatch(first_chars, '[!a-z]', /fold_case)
    
     ;; only fix weight columns that don't start with a letter
     wh_fix = where(res*bad_first eq 1, count)
     
     if count gt 0 then begin
        weight_par_names = 'TFORM' + string(wh_fix+1, format='(i1)')
        
        for j=0, count-1 do sxaddpar, hdr, weight_par_names[j], 'E'
        
        modfits, fits_file, 0, hdr, exten=1
        
        if not keyword_set(quiet) then print, 'File ' + fits_file + ' cleaned up successfully'
        
     endif else if not keyword_set(quiet) then print, 'No errors in file ' + fits_file
  endfor

end
