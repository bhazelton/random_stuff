pro polarization_figure, pub = pub, grey_scale = grey_scale


  froot = base_path('data') + 'fhd_simulations_old/'
  array = '496t'

  info_file = froot + 'sim_496t_info.idlsave'
  restore, info_file   

  sim_num = 0
  names = string([25, 50, 100, 200, 400, 600, 800, 1000], format = '(i04)')
  ;;names = 'offaxis'
  ;;names = 'fullsky'

  filebase = 'sim_496t_' + names[sim_num] + '_uvf'
  initial_uv_file = froot + filebase + '_initial.idlsave'

  truesky_uv_beamcorr_file = froot + filebase + '_truesky_uv_plane_beamcorr.idlsave'
  truesky_uv_file = froot + filebase + '_truesky_uv_plane.idlsave'

  holo_uv_file = froot + filebase + '_holo_uv_plane.idlsave'

  xrange = [.13,.23]
  yrange = xrange
  

  plotfile = base_path('plots') + 'single_use/polarization_figure.eps'

  ;; work out plot positions
  ;; positions = fltarr(4, 4)
  ;; cb_positions = fltarr(4, 1)

  ;; cb_size_punits = 0.1
  ;; margin_punits = [0.2, 0.2, 0.05, 0.1]

  ;; xlen_punits = 2.*(margin_punits[0] + margin_punits[2] + 1.) 
  ;; ylen_punits = 2.*(margin_punits[1] + margin_punits[3] + 1.) + cb_size_punits

  ;; plot_size = 1./[xlen_punits, ylen_punits]
  ;; cb_size = [1., cb_size_punits] / [xlen_punits, ylen_punits]
  ;; margin = [margin_punits[0:1] / [xlen_punits, ylen_punits], margin_punits[2:3] / [xlen_punits, ylen_punits]]

  ;; for i=0, 1 do begin
  ;;    positions[0, 2*i:2*(i+1)-1] = (i+1.)*margin[0] + i*(plot_size[0] + margin[2])
  ;;    positions[2, 2*i:2*(i+1)-1] = (i+1.)*(margin[0] + plot_size[0]) + i*margin[2]
  ;;    positions[1, [0,1]*2+i] = 1-[(i+1)*(margin[3]+plot_size[1]) + i*margin[1]]
  ;;    positions[3, [0,1]*2+i] = 1-[(i+1)*margin[3]+ i*(plot_size[1] + margin[1])]
  ;; endfor
  
  ;; cb_positions = [margin[0], margin[1], 1-margin[2], margin[1] + cb_size[1]]

  ;; size_factor = 200
  ;; xsize = xlen_punits * size_factor
  ;; ysize = ylen_punits * size_factor


  positions = fltarr(4, 4)
  for i=0, 1 do begin
     positions[0, 2*i:2*(i+1)-1] = i * 0.5
     positions[2, 2*i:2*(i+1)-1] = (i+1) * 0.5
     positions[1, [0,1]*2+i] = 1-((i+1)*0.5)
     positions[3, [0,1]*2+i] = 1-(i*0.5)
  endfor
  xsize = 700
  ysize = xsize

  window_num = 1

  if windowavailable(window_num) then begin 
     wset, window_num
     if !d.x_size ne xsize or !d.y_size ne ysize then make_win = 1 else make_win = 0
  endif else make_win = 1
  if make_win then cgdisplay, wid=window_num, xsize = xsize, ysize = ysize, color='white'
  cgerase, 'white'

  if keyword_set(pub) then begin
     charthick = 2
     thick = 2
     xthick = 2
     ythick = 2
     font = 1
     charsize = 1

     pson, file = plotfile, /eps 
  endif else begin
     charsize = 0.8

  endelse


  uvf_slice_plot, initial_uv_file, multi_pos = positions[*,0], multi_aspect = 1, plot_xrange = xrange, plot_yrange = yrange, $
                  title = 'Initial uv plane', grey_scale = grey_scale
  uvf_slice_plot, holo_uv_file, multi_pos = positions[*,2], multi_aspect = 1, plot_xrange = xrange, plot_yrange = yrange, $
                  title = 'Holographic uv plane', /color_0amp, grey_scale = grey_scale
  uvf_slice_plot, truesky_uv_beamcorr_file, multi_pos = positions[*,1], multi_aspect = 1, plot_xrange = xrange, plot_yrange = yrange, $
                  title = 'True Sky uv plane!C(before 0 weight areas removed)', /color_0amp, grey_scale = grey_scale
  uvf_slice_plot, truesky_uv_file, multi_pos = positions[*,3], multi_aspect = 1, plot_xrange = xrange, plot_yrange = yrange, $
                  title = 'True Sky uv plane', /color_0amp, grey_scale = grey_scale
 
  if keyword_set(pub) then begin
     psoff
     wdelete, window_num
  endif




  ;; uvf_slice_plot, initial_uv_file, window_num = 1, plot_xrange = xrange, plot_yrange = yrange, $
  ;;                 title = 'Initial uv plane', grey_scale = grey_scale
  ;; uvf_slice_plot, holo_uv_file, window_num = 2, plot_xrange = xrange, plot_yrange = yrange, $
  ;;                 title = 'Holographic uv plane', /mark_0amp, grey_scale = grey_scale
  ;; uvf_slice_plot, truesky_uv_beamcorr_file, window_num = 3, plot_xrange = xrange, plot_yrange = yrange, $
  ;;                 title = 'True Sky uv plane (before 0 weights removed)', /mark_0amp, grey_scale = grey_scale
  ;; uvf_slice_plot, truesky_uv_file, window_num = 4, plot_xrange = xrange, plot_yrange = yrange, $
  ;;                 title = 'True Sky uv plane', /mark_0amp, grey_scale = grey_scale


end
