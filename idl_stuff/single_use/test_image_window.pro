uvf_folder = '/Volumes/Data1/fhd_eppsilon_test_data/fhd_nb_Aug2017_savedbp_w_cable_w_digjump/ps_master/data/uvf_cubes/'
file_start = 'Combined_obs_Aug23_longrunstyle_'
even_label = 'even_cubeXX_'
odd_label = 'odd_cubeXX_'
file_end = 'uvf.idlsave'

big_bh_label = 'bh_'
big_label = 'fullimg_'
small_label = ''

dirty_label = 'dirty_'
model_label = 'model_'
res_label = 'res_'
wt_label = 'weights_'

big_bh_res_file = uvf_folder + file_start + even_label + big_bh_label + res_label + file_end
big_res_file = uvf_folder + file_start + even_label + big_label + res_label + file_end
small_res_file = uvf_folder + file_start + even_label + small_label + res_label + file_end

big_bh_res_odd_file = uvf_folder + file_start + odd_label + big_bh_label + res_label + file_end
big_res_odd_file = uvf_folder + file_start + odd_label + big_label + res_label + file_end
small_res_odd_file = uvf_folder + file_start + odd_label + small_label + res_label + file_end

big_bh_dirty_file = uvf_folder + file_start + even_label + big_bh_label + dirty_label + file_end
big_dirty_file = uvf_folder + file_start + even_label + big_label + dirty_label + file_end
small_dirty_file = uvf_folder + file_start + even_label + small_label + dirty_label + file_end

big_bh_model_file = uvf_folder + file_start + even_label + big_bh_label + model_label + file_end
big_model_file = uvf_folder + file_start + even_label + big_label + model_label + file_end
small_model_file = uvf_folder + file_start + even_label + small_label + model_label + file_end

big_bh_wt_file = uvf_folder + file_start + even_label + big_bh_label + wt_label + file_end
big_wt_file = uvf_folder + file_start + even_label + big_label + wt_label + file_end
small_wt_file = uvf_folder + file_start + even_label + small_label + wt_label + file_end


big_bh_res = getvar_savefile(big_bh_res_file, 'data_cube')
big_res = getvar_savefile(big_res_file, 'data_cube')
small_res = getvar_savefile(small_res_file, 'data_cube')

big_bh_res_odd = getvar_savefile(big_bh_res_odd_file, 'data_cube')
big_res_odd = getvar_savefile(big_res_odd_file, 'data_cube')
small_res_odd = getvar_savefile(small_res_odd_file, 'data_cube')

big_bh_res_diff = big_bh_res - big_bh_res_odd
big_res_diff = big_res - big_res_odd
small_res_diff = small_res - small_res_odd

big_bh_dirty = getvar_savefile(big_bh_dirty_file, 'data_cube')
big_dirty = getvar_savefile(big_dirty_file, 'data_cube')
small_dirty = getvar_savefile(small_dirty_file, 'data_cube')

big_bh_model = getvar_savefile(big_bh_model_file, 'data_cube')
big_model = getvar_savefile(big_model_file, 'data_cube')
small_model = getvar_savefile(small_model_file, 'data_cube')


big_bh_wt = getvar_savefile(big_bh_wt_file, 'weights_cube')
big_wt = getvar_savefile(big_wt_file, 'weights_cube')
small_wt = getvar_savefile(small_wt_file, 'weights_cube')

big_bh_var = getvar_savefile(big_bh_wt_file, 'variance_cube')
big_var = getvar_savefile(big_wt_file, 'variance_cube')
small_var = getvar_savefile(small_wt_file, 'variance_cube')

big_bh_var_ave = mean(big_bh_var, dimension=3)
big_var_ave = mean(big_var, dimension=3)
small_var_ave = mean(small_var, dimension=3)

big_bh_var_calc = stddev(real_part(big_bh_res_diff), dimension=3)^2
big_var_calc = stddev(real_part(big_res_diff), dimension=3)^2
small_var_calc = stddev(real_part(small_res_diff), dimension=3)^2


dims = size(big_bh_wt, /dim)
x_inds = indgen(ceil(dims[0]/3.))*3
y_inds = indgen(ceil(dims[1]/3.))*3

big_bh_dec_res = big_bh_res[x_inds, *, *]
big_bh_dec_res = big_bh_dec_res[*, y_inds, *]

big_bh_dec_wt = big_bh_wt[x_inds, *, *]
big_bh_dec_wt = big_bh_dec_wt[*, y_inds, *]

big_bh_dec_var = big_bh_var[x_inds, *, *]
big_bh_dec_var = big_bh_dec_var[*, y_inds, *]

quick_image, abs(big_bh_wt[*,*,0]),/log, title='big BH weight', window=1
quick_image, abs(big_wt[*,*,0]),/log, title='big weight', window=2
quick_image, abs(small_wt[*,*,0]),/log, title='small weight', window=3
quick_image, abs(big_bh_dec_wt[*,*,0]),/log, title='big BH weight decimated', window=4

quick_image, abs(big_bh_res[*,*,0]),/log, title='big BH residual', window=1
quick_image, abs(big_res[*,*,0]),/log, title='big residual', window=2
quick_image, abs(small_res[*,*,0]),/log, title='small residual', window=3
quick_image, abs(big_bh_dec_res[*,*,0]),/log, title='big BH residual decimated', window=4

quick_image, abs(big_bh_dirty[*,*,0]),/log, title='big BH dirty', window=1
quick_image, abs(big_dirty[*,*,0]),/log, title='big dirty', window=2
quick_image, abs(small_dirty[*,*,0]),/log, title='small dirty', window=3

quick_image, abs(big_bh_model[*,*,0]),/log, title='big BH model', window=1
quick_image, abs(big_model[*,*,0]),/log, title='big model', window=2
quick_image, abs(small_model[*,*,0]),/log, title='small model', window=3

quick_image, abs(big_bh_res[*,*,0]/big_bh_wt[*,*,0]),/log, title='big BH residual/weight', window=1
quick_image, abs(big_res[*,*,0]/big_wt[*,*,0]),/log, title='big residual/weight', window=2
quick_image, abs(small_res[*,*,0]/small_wt[*,*,0]),/log, title='small residual/weight', window=3
quick_image, abs(big_bh_dec_res[*,*,0]/big_bh_dec_wt[*,*,0]),/log, title='big BH residual/weight decimated', window=4

quick_image, abs(big_bh_dirty[*,*,0]/big_bh_wt[*,*,0]),/log, title='big BH dirty/weight', window=1
quick_image, abs(big_dirty[*,*,0]/big_wt[*,*,0]),/log, title='big dirty/weight', window=2
quick_image, abs(small_dirty[*,*,0]/small_wt[*,*,0]),/log, title='small dirty/weight', window=3

quick_image, abs(big_bh_model[*,*,0]/big_bh_wt[*,*,0]),/log, title='big BH model/weight', window=1
quick_image, abs(big_model[*,*,0]/big_wt[*,*,0]),/log, title='big model/weight', window=2
quick_image, abs(small_model[*,*,0]/small_wt[*,*,0]),/log, title='small model/weight', window=3

quick_image, abs(big_bh_var[*,*,0]/(big_bh_wt[*,*,0])^2),/log, title='big BH var/weight^2', window=1
quick_image, abs(big_var[*,*,0]/(big_wt[*,*,0])^2),/log, title='big var/weight^2', window=2
quick_image, abs(small_var[*,*,0]/(small_wt[*,*,0])^2),/log, title='small var/weight^2', window=3
quick_image, abs(big_bh_dec_var[*,*,0]/(big_bh_dec_wt[*,*,0])^2),/log, title='big BH var/weight^2 decimated', window=4

quick_image, big_bh_var_calc,/log, title='big BH var calc', window=1
quick_image, big_var_calc,/log, title='big var calc', window=2
quick_image, small_var_calc,/log, title='small var calc', window=3

quick_image, big_bh_var_calc/big_bh_var_ave,/log, title='big BH var calc', window=1
quick_image, big_var_calc/big_var_ave,/log, title='big var calc', window=2
quick_image, small_var_calc/small_var_ave,/log, title='small var calc', window=3
