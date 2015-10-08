pro fhd_source_plots, filename

  if file_test(filename) eq  0 then message, 'file not found'
  
  ;#id  x_loc y_loc RA  Dec S/N radius  avg_beam  XX_apparent YY_apparent Stokes_I_fit  Stokes_Q_fit  Stokes_I_res  Stokes_Q_res  Extended
  
  textfast,data,header,file_path=filename,/read, first_line=1, column_list=[0, 3, 4] ;; id, RA, DEC
  
  ra_vals = data[1,*]
  dec_vals = data[2,*]
  
  wh_gt180 = where(ra_vals gt 180, count_gt180)
  if count_gt180 gt 0 then ra_vals[wh_gt180] = ra_vals[wh_gt180] - 360.
  
  cgplot, ra_vals, dec_vals, psym=1
  
end