pyfhd_test_path = "/Users/bryna/Projects/Physics/data_files/pyfhd_test_data/"

data_dir = pyfhd_test_path + "gridding/visibility_degrid/"

test_input1 = data_dir+"input_1.sav"
test_input2 = data_dir+"input_2.sav"
test_input3 = data_dir+"input_3.sav"

test_output1 = data_dir+"output_1_new.sav"
test_output2 = data_dir+"output_2_new.sav"
test_output3 = data_dir+"output_3_new.sav"

; image_uv1 is the same as image_uv3. image_uv2 is different
image_uv1 = getvar_savefile(test_input1, "image_uv")

; vis_weights are all identical, just use the first one
vis_weights = getvar_savefile(test_input1, "vis_weight_ptr")
; obs are all identical, just use the first one
obs = getvar_savefile(test_input1, "obs")

; psf 1 & 2 are the same, psf3 has a psf.image_info.image_power_beam_arr defined.
; not clear that matters...
psf1 = getvar_savefile(test_input1, "psf")

; params are all identical, just use the first one
params = getvar_savefile(test_input1, "params")

; polarization1=0, polarization2=2, polarization3=0
polarization1 = getvar_savefile(test_input1, "polarization")

; fill_model_visibilities are all 1. use the first one
fill_model_visibilities = getvar_savefile(test_input1, "fill_model_visibilities")

; spectral_model_uv_arr are all 1. use the first one
spectral_model_uv_arr = getvar_savefile(test_input1, "spectral_model_uv_arr")

; conserve_memory are all the same. use the first one
conserve_memory = getvar_savefile(test_input1, "conserve_memory")

beam_per_baseline = 0
uv_grid_phase_only = 1

vis_return_orig1 = getvar_savefile(test_output1, "vis_return")
vis_return1 = visibility_degrid( $
    image_uv1, $
    vis_weights, $
    obs, $
    psf1, $
    params, $
    polarization=polarization1, $
    fill_model_visibilities=fill_model_visibilities, $
    spectral_model_uv_arr=spectral_model_uv_arr, $
    beam_per_baseline=beam_per_baseline, $
    uv_grid_phase_only=uv_grid_phase_only, $
    conserve_memory=conserve_memory)

max(abs(*vis_return1 - *vis_return_orig1))
; test_output1_new = data_dir+"output_1_new.sav"
; vis_return = vis_return1
; save, filename=test_output1_new, vis_return



image_uv2 = getvar_savefile(test_input2, "image_uv")
psf2 = getvar_savefile(test_input2, "psf")
polarization2 = getvar_savefile(test_input2, "polarization")

vis_return_orig2 = getvar_savefile(test_output2, "vis_return")
vis_return2 = visibility_degrid( $
    image_uv2, $
    vis_weights, $
    obs, $
    psf2, $
    params, $
    polarization=polarization2, $
    fill_model_visibilities=fill_model_visibilities, $
    spectral_model_uv_arr=spectral_model_uv_arr, $
    beam_per_baseline=beam_per_baseline, $
    uv_grid_phase_only=uv_grid_phase_only, $
    conserve_memory=conserve_memory)

max(abs(*vis_return2 - *vis_return_orig2))
; test_output2_new = data_dir+"output_2_new.sav"
; vis_return = vis_return2
; save, filename=test_output2_new, vis_return


image_uv3 = getvar_savefile(test_input3, "image_uv")
psf3 = getvar_savefile(test_input3, "psf")
polarization3 = getvar_savefile(test_input3, "polarization")
beam_per_baseline = 1
psf3.interpolate_kernel = 0

vis_return_orig3 = getvar_savefile(test_output3, "vis_return")
vis_return3 = visibility_degrid( $
    image_uv3, $
    vis_weights, $
    obs, $
    psf3, $
    params, $
    polarization=polarization3, $
    fill_model_visibilities=fill_model_visibilities, $
    spectral_model_uv_arr=spectral_model_uv_arr, $
    beam_per_baseline=beam_per_baseline, $
    uv_grid_phase_only=uv_grid_phase_only, $
    conserve_memory=conserve_memory)

max(abs(*vis_return3 - *vis_return_orig3))
; test_output3_new = data_dir+"output_3_new.sav"
; vis_return = vis_return3
; save, filename=test_output3_new, vis_return
