pro healpix_idlsave2fits

  filename = base_path('data') + '1061316176_even_cube.sav'
  
  file_basename = cgrootname(filename, directory=froot, extension=extension)
  fits_filename = froot + file_basename + '.fits'
  
  ;restore, filename
  
  void = getvar_savefile(filename, names = varnames)
  ;wh_obs = where(strmatch(varnames, 'obs*',/fold_case) gt 0, count_obs)
  ;obs_arr = getvar_savefile(filename, varnames[wh_obs[0]])
  
  print, 'reading data'
  ordering = 'ring'
  t0 = systime(1)
  nside = getvar_savefile(filename, 'nside')
  pixel = getvar_savefile(filename, 'hpx_inds')
  signal = (getvar_savefile(filename, 'model_xx_cube'))[*,0]
  weights = (getvar_savefile(filename, 'weights_xx_cube'))[*,0]
  variances = (getvar_savefile(filename, 'variance_xx_cube'))[*,0]
  t1 = systime(1)
  
  ;ring2nest, nside, hpx_inds, hpx_inds_nest
  
  ;write_fits_cut4,fits_filename,pixels,signal,weights,variances,nside=nside, ordering=ordering,coord='C'
  
  local_header = [' ']
  ; insert units
  add_units_fits, local_header, units='      ', col=1, err=err
  add_units_fits, local_header, units='Jy    ', col=2, err=err
  add_units_fits, local_header, units='      ', col=3, err=err
  add_units_fits, local_header, units='      ', col=4, err=err
  if (err ne 0) then message,'Error while writing header'
  
  add_ordering_fits, local_header, nested=nested, ring=ring, ordering=ordering,error=error
  if (error ne 0) then message,'Error while writing header'
  
  ;; add 3rd dimension info
  
  
  ; create structures
  prim_st = 0
  exten_st = create_struct('HDR',local_header,'PIXEL',round(pixel),'SIGNAL',signal*1.0,'N_OBS',round(weights),'SERROR',variances*1.0)
  
  print, 'writing data'
  t2 = systime(1)
  write_fits_sb, fits_filename, prim_st, exten_st,  $
    Coordsys='C', Nside = nside, /partial, extension=0
  t3 = systime(1)
  
  print, 'reading time', t1-t0
  print, 'writing time', t3-t2
  
end