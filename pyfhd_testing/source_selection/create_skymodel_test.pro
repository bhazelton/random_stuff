pyfhd_test_path = "/Users/bryna/Projects/Physics/data_files/pyfhd_test_data/"

data_dir = pyfhd_test_path + "gridding/visibility_degrid/"

test_input1 = data_dir+"input_1.sav"

obs = getvar_savefile(test_input1, "obs")
psf = getvar_savefile(test_input1, "psf")

catalog_path = "/Users/bryna/Projects/Physics/FHD/catalog_data/GLEAM_v2_plus_rlb2019.sav"

source_list = generate_source_cal_list(obs, psf, catalog_path=catalog_path)



; setup unit test
catalog_path = "/Users/bryna/Projects/Physics/FHD/catalog_data/GLEAM_v2_plus_rlb2019.sav"

catalog = load_source_catalog(catalog_path, varname='catalog')

test_file = "gleam_v2_rlb2019_cut.sav"
n_src = n_elements(catalog)

; randomly select 100th of the sources for the test catalog
inds_use = cgRandomIndices(n_src, n_src/100, SEED=100)
catalog = catalog[inds_use]
save, catalog, filename = test_file


catalog_path="gleam_v2_rlb2019_cut.sav"
obs_file = "/Users/bryna/Projects/Physics/pyfhd-datasets/fhd_extracts/MWA/std_1061316296/1061316296_obs.sav"
psf_file = "/Users/bryna/Projects/Physics/pyfhd-datasets/fhd_extracts/MWA/std_1061316296/cut_down_psf.sav"

obs = getvar_savefile(obs_file, "obs")
psf = getvar_savefile(psf_file, "psf")

; set obs freq_array to match psf freq_array
freqs = psf.freq
n_freq = n_elements(freqs)
obs.n_freq = n_freq
obs.nf_vis = obs.nf_vis[*, 192:193]
obs.freq_center = mean(freqs)

orig_bi = (*obs.baseline_info)
new_bi = {tile_A:orig_bi.tile_A,tile_B:orig_bi.tile_B,bin_offset:orig_bi.bin_offset,Jdate:orig_bi.Jdate,freq:freqs,fbin_i:lindgen(n_freq),$
freq_use:fltarr(n_freq) + 1,tile_use:orig_bi.tile_use,time_use:orig_bi.time_use,tile_names:orig_bi.tile_names,tile_height:orig_bi.tile_height,tile_flag:orig_bi.tile_flag}
obs.baseline_info = ptr_new(new_bi)

source_array = generate_source_cal_list(obs, psf, catalog_path=catalog_path)

cal_src_list_file = "gleam_v2_rlb2019_cut_cal_src_list_end.sav"
save, filename = cal_src_list_file, source_array