pyfhd_test_path = "/Users/bryna/Projects/Physics/data_files/pyfhd_test_data/"

data_dir = pyfhd_test_path + "gridding/visibility_grid/"

test_input = data_dir+"input_1.sav"

test_output = data_dir+"output_1.sav"

visibility_ptr = getvar_savefile(test_input, "visibility_ptr")

vis_weights = getvar_savefile(test_input, "vis_weight_ptr")

obs = getvar_savefile(test_input, "obs")

psf = getvar_savefile(test_input, "psf")

params = getvar_savefile(test_input, "params")

polarization = getvar_savefile(test_input, "polarization")

_ = getvar_savefile(test_input, names=names)

if max(strmatch(names, "uniform_flat", /fold_case)) > 0 then uniform_flat = getvar_savefile(test_input, "uniform_flag")

if max(strmatch(names, "no_conjugate", /fold_case)) > 0 then no_conjugate = getvar_savefile(test_input, "no_conjugate")

if max(strmatch(names, "model_ptr", /fold_case)) > 0 then model_ptr = getvar_savefile(test_input, "model_ptr")

if max(strmatch(names, "fi_use", /fold_case)) > 0 then fi_use = getvar_savefile(test_input, "fi_use")

if max(strmatch(names, "bi_use", /fold_case)) > 0 then bi_use = getvar_savefile(test_input, "bi_use")

uv_image_return_orig = getvar_savefile(test_output, "image_uv")

uv_image_return = visibility_grid( $
    visibility_ptr, $
    vis_weights, $
    obs, $
    status_str, $
    psf, $
    params, $
    polarization=polarization, $
    uniform_flag=uniform_flag, $
    no_conjugate=no_conjugate, $
    model_ptr=model_ptr, $
    fi_use=fi_use, $
    bi_use=bi_use)

max(abs(*uv_image_return - *uv_image_return_orig))

stop

