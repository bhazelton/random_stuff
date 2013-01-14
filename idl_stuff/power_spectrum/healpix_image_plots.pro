

pro healpix_image_plots, fitsfile, projection = projection

  read_fits_map, fitsfile, map_out, hdr, ehdr, nside=nside, ordering=ordering, coordsys=coordsys

  allowed_projections = ['orthographic', 'mollweide', 'gnomonic', 'cartesian']
  if n_elements(projection) eq 0 then projection = 'orthographic'

  wh_proj = where(allowed_projections eq projection, count)
  if count eq 0 then message, 'Projection must be one of: ' + allowed_projections

  ;; get current colortable so it can be restored later
  tvlct, r, g, b, /get
  loadct,39

  case projection of
     'orthographic': begin
          orthview, map_file, 'Stokes I', /graticule, glsize = 1, rot = [90, -45], /half_sky, window=1
          orthview, map_file, 'Stokes U', /graticule, glsize = 1, rot = [90, -45], /half_sky, window=2
       end

     'mollweide': begin
        mollview, map_file, 'Stokes I', /graticule, glsize = 1, rot = [90, -45]
        mollview, map_file, 'Stokes U', /graticule, glsize = 1, rot = [90, -45]
     end

     'gnomonic': begin
        gnomview, map_file, 'Stokes I'
        gnomview, map_file, 'Stokes U'
     end

     'cartesian': begin
        cartview, map_file, 'Stokes I'
        cartview, map_file, 'Stokes U'
     end
  endcase

  pix_nums = long(map_out[*,0])
  pix2ang_ring, nside, pix_nums, pix_center_theta, pix_center_phi

  data = double(map_out[*,1])
  
  ;; window, 3, xsize = 600, ysize = 600
  ;; map_set, -45, -90, /orthographic, /grid
     
  ;;convert to degrees E lon & lat [-90 -> 90]
  lons = -1 * pix_center_phi * 180. / !pi
  lats = 90 - pix_center_theta * 180. / !pi
     
  ;; for i=0L, obs_npix-1 do oplot, [lons[i]], [lats[i]], psym = 3, color = round(254*(obs_npix-i)/obs_npix)
     
  data_minmax = minmax(data)
  data_colors = round(254*(data-data_minmax[0])/(data_minmax[1] - data_minmax[0]))
  ;;for i=0L, obs_npix-1 do oplot, [lons[i]], [lats[i]], psym = 3, color = data_colors[i]
     
     
  window, 4, xsize = 600, ysize = 600
  plot, lons, lats, psym =3
  for i=0L, obs_npix-1 do oplot, [lons[i]], [lats[i]], psym = 3, color = round(254*(obs_npix-i)/obs_npix)
  for i=0L, obs_npix-1 do oplot, [lons[i]], [lats[i]], psym = 3, color = data_colors[i]
  





end
