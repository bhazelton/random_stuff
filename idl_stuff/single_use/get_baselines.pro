function get_baselines, freq_mhz = freq_mhz, quiet = quiet

  filepath = base_path() + 'fhd_sims/'
  
  textfast,tile_locs,/read,file_path=filepath + '32T_tile_locations', first_line=1,column=indgen(3)+2
 
  ;;textfast,tile_locs,/read,filename= filepath +'496T_tile_locations', first_line=1,column=indgen(2)+1
  ;;textfast,tile_locs,/read,filename= filepath +'512T_tile_locations', first_line=1,column=indgen(2)+1

  xtile=reform(tile_locs[0,*])
  ytile=reform(tile_locs[1,*])

  n_tiles = n_elements(xtile)

  tile_A = reform(rebin(indgen(n_tiles)+1, n_tiles, n_tiles), n_tiles*n_tiles)
  tile_B = reform(rebin(reform(indgen(n_tiles)+1, 1, n_tiles), n_tiles, n_tiles), n_tiles*n_tiles)
  
  ;; get unique baselines
  temp = double(tile_A)^3/double(tile_B) + double(tile_B)^3/double(tile_A)
  uniq_inds = uniq(temp, sort(temp))
  nbaselines = n_elements(uniq_inds)
  
  if nbaselines ne ((n_tiles^2-n_tiles)/2 + n_tiles) then stop
  
  tile_A = tile_A[uniq_inds]
  tile_B = tile_B[uniq_inds]
  baseline_arr = tile_A*256 + tile_B
  
  uu_arr=(xtile[tile_A-1]-xtile[tile_B-1])
  vv_arr=(ytile[tile_A-1]-ytile[tile_B-1])

  if n_elements(freq_mhz) eq 1 then begin
     c_light=299792458.
     uu_arr = uu_arr * freq_mhz * 1e6 / c_light
     vv_arr = vv_arr * freq_mhz * 1e6 / c_light
     xtitle = 'baseline length ' + textoidl('(\lambda)', font = font) + ' (' + number_formatter (freq_mhz) + ' Mhz)'
  endif else xtitle = 'baseline length (m)'

  baseline_dist = sqrt(uu_arr^2. + vv_arr^2.)
  
  ;; remove autos
  baseline_dist = baseline_dist[where(baseline_dist gt 0)]
  
  if not keyword_set(quiet) then quick_histplot, baseline_dist, /logdata, xtitle = xtitle, binsize=0.1

  return, baseline_dist

end
