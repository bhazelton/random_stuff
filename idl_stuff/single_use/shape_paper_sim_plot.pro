pro shape_paper_sim_plot, baseline_axis = baseline_axis, pub = pub, grey_scale = grey_scale


  froot = base_path() + 'fhd/simulations/'

  sim_num = 7
  baseline_layout = ['simple', 'single']
  baseline_spacing = [4,30]

  names = string([3, 6, 13, 25, 50, 75, 100, 125, 250, 500], format = '(i04)')
  offsets = [3.16, 6.45, 12.65, 25.09, 49.91, 75.3, 100.04, 124.58, 250.26, 499.18]
  
  spacing_name = strtrim_plus(string(baseline_spacing), flag = 2)
  array = baseline_layout + spacing_name
  tile_tag = '_' + array
  

  info_file = froot + 'sim_' + array + '_info.idlsave'
  restore, info_file[0]
  deg_offsets = offsets * degpix
  rad_offsets =  deg_offsets *!pi/180d
  deg_offset_str = string(deg_offsets, format = '(g4.2)')

  fadd = '_holo'
  fbase_arr = 'sim' + tile_tag + '_' + names[sim_num] + '_uvf'
  uf_savefiles = froot + fbase_arr + fadd + '_uf_plane.idlsave'

  new_uf_savefile = froot + 'sim' + tile_tag[0] + tile_tag[1] + '_' + names[sim_num] + '_spliced_uf.idlsave'

  plotfile_path = base_path() + 'single_use/'
  plotfile = plotfile_path + 'sim' + tile_tag[0] + tile_tag[1] + '_' + names[sim_num] + '_spliced_uf.eps'

  restore, uf_savefiles[0]
  spliced_uf = reform(uvf_slice, n_elements(xarr), n_elements(yarr))

  xdelta = xarr[1] - xarr[0]
  xarr_edges = [xarr - xdelta/2, max(xarr) + xdelta/2]
  temp = where(total(abs(reform(uvf_slice, n_elements(xarr), n_elements(yarr))),2) gt 0, count)
  if count gt 1 then plot_xrange =  minmax(xarr_edges[[temp, max(temp)+1]])

  restore, uf_savefiles[1]
  wh_pos = where(xarr gt 0)
  spliced_uf[wh_pos, *] = (reform(uvf_slice, n_elements(xarr), n_elements(yarr)))[wh_pos, *]

  uvf_slice = spliced_uf
  save, file = new_uf_savefile, uvf_slice, kperp_lambda_conv, slice_axis, slice_inds, xarr, yarr, slice_name, plane_name, $
        plot_xname, plot_yname

  uvf_slice_plot, new_uf_savefile, pub = pub, plotfile = plotfile, window_num = 1, $
                  title = ' phase [' + deg_offset_str[sim_num] + '!Uo!N, 0!Uo!N]', grey_scale = grey_scale, $
                  baseline_axis = baseline_axis, plot_xrange = plot_xrange, /mark_0

end
