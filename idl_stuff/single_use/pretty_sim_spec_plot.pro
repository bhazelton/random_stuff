pro pretty_sim_spec_plot

input_file = '~/FHD/catalog_data/eor_power_1d.idlsave'
file = '/data4/MWA/FHD_Aug23/fhd_bjh_arrsim_eor_realsky/ps/1061316176_gridded_uvf__even_odd_joint_model_xx_bh_kperp_density_min_gt1_1dkpower.idlsave'

file_arr = [file, input_file]
psyms = [10, -3]
titles_use = ['simulation', 'EoR signal']
colors = ['red', 'black']

plotfile_use = '/data4/MWA/FHD_Aug23/fhd_bjh_arrsim_eor_realsky/ps/plots/1061316176_gridded_uvf__even_odd_joint_model_xx_bh_kperp_density_min_gt1_1dkpower_pretty.pdf'


kpower_1d_plots, file_arr, window_num = 1, colors = colors, names = titles_use, psyms = psyms, /hinv, $
        png = png, eps = eps, pdf = pdf, plotfile = plotfile_use, k_range = [2e-2, 3], data_range = [1e2,1e6]

end