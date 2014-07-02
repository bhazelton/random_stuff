pro cal_4th_line_test

  filename = base_path('data') + 'fhd_ps_data/128T_cubes/bandpass_txt/1061316296_bandpass.txt'
  cube_filename = base_path('data') + 'fhd_ps_data/128T_cubes/aug23_3hr_decon_bp/Combined_obs_1061311664-1061323008_cube__even_odd_joint_dirty_yy_bh_power.idlsave'
  power_cube = getvar_savefile(cube_filename, 'power_3d')
  textfast, data, file_path = filename,/read
  mask = float(reform(data[1,*]) gt 0)
  mask = float(data[1,*] gt 0)
  maskstep = mask
  maskstep[256:*] = maskstep[256:*]*.9
  ft_data = fft_shift(fft(data[1,*]))
  ft_mask = fft_shift(fft(mask))
  ft_maskstep = fft_shift(fft(maskstep))
  power_line = power_cube[99, 0, *]
  ft_4th = ft_mask*0
  ft_4th[192+38] = 0.1
  mask_4th = (sin(findgen(1,384)*2.*!pi/10.1)*.05+mask)*mask
  
  cgplot, data[0,*], data[1,*], psym=3, yrange = [0.8,1.2]
  cgplot, data[0,*], mask, color ='red',/over
  cgplot, data[0,*], maskstep, color ='blue',/over
  window,/free
  cgplot, abs(ft_data[192:288]), yrange = [-0.1, 0.8]
  cgplot, abs(ft_mask[192:288]), /over, color='red'
  cgplot, abs(ft_maskstep[192:288]), /over, color='blue'
  cgplot, abs(ft_4th[192:288]), /over, color='green'
  cgplot, alog10(power_line)/8., /over, color='orange'
  window,/free
  cgplot, sin(findgen(1,384)*2.*!pi/10.1)*.1+mask
  cgplot, mask_4th
  
end