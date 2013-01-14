

pro survey_image, pub = pub

  froot = base_path() + 'power_spectrum/healpix_maps/'
  filenames = 'mosaic_restored_' + ['0', '1'] + '.fits'
  fitsfiles = froot + filenames

  map_header1 = headfits(fitsfiles[0])
  dims1 = sxpar(map_header1, 'NAXIS*')
  ref_pixels1 = sxpar(map_header1,'CRPIX*') ;; these are 1-based NOT 0-based
  ref_pix_vals1 = sxpar(map_header1,'CRVAL*')
  ref_pix_deltas1 = sxpar(map_header1,'CDELT*')

  x1_coords = (findgen(dims1[0])-ref_pixels1[0])*ref_pix_deltas1[0] + ref_pix_vals1[0]
  y1_coords = (findgen(dims1[1])-ref_pixels1[1])*ref_pix_deltas1[1] + ref_pix_vals1[1]

  data1 = readfits(fitsfiles[0],/silent)
  mask1 = data1[*,*,4]
  map1 = data1[*,*,0]

  map_header2 = headfits(fitsfiles[1])
  dims2 = sxpar(map_header2, 'NAXIS*')
  ref_pixels2 = sxpar(map_header2,'CRPIX*') ;; these are 1-based NOT 0-based
  ref_pix_vals2 = sxpar(map_header2,'CRVAL*')
  ref_pix_deltas2 = sxpar(map_header2,'CDELT*')

  x2_coords = (findgen(dims2[0])-ref_pixels2[0])*ref_pix_deltas2[0] + ref_pix_vals2[0]
  y2_coords = (findgen(dims2[1])-ref_pixels2[1])*ref_pix_deltas2[1] + ref_pix_vals2[1]

  data2 = readfits(fitsfiles[1],/silent)
  mask2 = data2[*,*,4]
  map2 = data2[*,*,0]


  wh1 = where(mask1 gt 0, count1)
  x1_inds = wh1 mod dims1[0]
  y1_inds = wh1 / dims1[0]

  xyad, map_header1, x1_inds, y1_inds, ra_vals1, dec_vals1, /celestial, print=0

  wh2 = where(mask2 gt 0, count2)
  x2_inds = wh2 mod dims2[0]
  y2_inds = wh2 / dims2[0]

  xyad, map_header2, x2_inds, y2_inds, ra_vals2, dec_vals2, /celestial, print=0

  ra_vals = [ra_vals2-360d, ra_vals1]
  dec_vals = [dec_vals2, dec_vals1]
  data = [data2[wh2], data1[wh1]]

  tvlct, r, g, b, /get
  loadct, 39

  minmax_data = minmax(data)
  colors = (data - minmax_data[0]) * 254/(minmax_data[1] - minmax_data[0])

  if keyword_set(pub) then begin
     charthick = 3
     thick = 3
     xthick = 3
     ythick = 3
     charsize = 2
     font = 1
 
     plotfile = base_path() + 'power_spectrum/plots/survey/survey_image.eps'

     window, 1, xsize = 600, ysize = 300
     pson, file = plotfile, /eps 
  endif else window, 1, xsize = 600, ysize = 300


  map_set, position =[0.1,0.1,0.9,0.9], limit=[-45, -60, 0, 110]
  PlotS, ra_vals, dec_vals, PSym=3, Color=colors, SymSize=0.5
  map_grid, /box, color = 0, label=1

  if keyword_set(pub) then begin
     psoff
     wdelete, 1
  endif

  tvlct, r, g, b

end
