from logging import Logger

import numpy as np
import matplotlib.pyplot as plt
from pyradiosky import SkyModel

from pyfhd.io.pyfhd_io import recarray_to_dict
from pyfhd.source_modeling.source_utils import create_skymodel
from pyfhd.pyfhd_tools.test_utils import get_savs, sav_file_rearrange_psf


catalog_path = (
    "/Users/bryna/Projects/Physics/FHD/catalog_data/GLEAM_v2_plus_rlb2019.sav"
)

pyfhd_test_path = "/Users/bryna/Projects/Physics/data_files/pyfhd_test_data/"

data_dir = pyfhd_test_path + "gridding/visibility_degrid/"

h5_save_dict = get_savs(data_dir, "input_1.sav")

h5_save_dict["psf"] = sav_file_rearrange_psf(h5_save_dict["psf"])
h5_save_dict = recarray_to_dict(h5_save_dict)

h5_save_dict["obs"]["n_baselines"] = h5_save_dict["obs"]["nbaselines"]
h5_save_dict["obs"]["dimension"] = int(h5_save_dict["obs"]["dimension"])


# Transpose the model if it exists
h5_save_dict["pyfhd_config"] = {
    "interpolate_kernel": h5_save_dict["psf"]["interpolate_kernel"],
    "psf_dim": h5_save_dict["psf"]["dim"],
    "psf_resolution": h5_save_dict["psf"]["resolution"],
    "beam_mask_threshold": h5_save_dict["psf"]["beam_mask_threshold"],
    "beam_clip_floor": h5_save_dict["extra"]["beam_clip_floor"],
    # need this to be defined (not actually used)
    "image_filter": "filter_uv_uniform",
}

h5_save_dict["vis_weight_ptr"] = h5_save_dict["vis_weight_ptr"].T

refraction = "idl"
sky = create_skymodel(
    obs=h5_save_dict["obs"],
    psf=h5_save_dict["psf"],
    logger=Logger,
    catalog_path=catalog_path,
    refraction=refraction,
)

ref_str = ""
if refraction is not None:
    ref_str = f"_refract_{refraction}"

idl_dbl = True

idl_str = ""
if idl_dbl:
    idl_str = "_dbl"

# FHD catalog after cuts:
# fhd_cat_file = "/Users/bryna/Projects/Physics/FHD/cal_source_list_after_cuts.sav"
fhd_cat_file = (
    f"/Users/bryna/Projects/Physics/FHD/cal_source_list_after_cuts{idl_str}.sav"
)
# pyfhd_cat_file = "/Users/bryna/Projects/Physics/PyFHD/skymodel_after_cuts_fk5.skyh5"
pyfhd_cat_file = (
    f"/Users/bryna/Projects/Physics/PyFHD/skymodel_after_cuts{ref_str}_fk5.skyh5"
)

fhd_sky = SkyModel.from_file(
    fhd_cat_file, extra_columns={"x": "image_x", "y": "image_y", "beam_I": "beam_I"}
)
pyfhd_sky = SkyModel.from_file(pyfhd_cat_file)

fhd_in_pyfhd = np.isin(fhd_sky.name.astype(int), pyfhd_sky.name.astype(int))
pyfhd_in_fhd = np.isin(pyfhd_sky.name.astype(int), fhd_sky.name.astype(int))

fhd_unique_inds = np.nonzero(~fhd_in_pyfhd)[0]
pyfhd_unique_inds = np.nonzero(~pyfhd_in_fhd)[0]

plotfile = f"fhd_pyfhd_catalog_compare{idl_str}{ref_str}.png"

fig = plt.figure()
ax = plt.subplot(111)

ax.scatter(
    fhd_sky.ra.wrap_at("180d").deg,
    fhd_sky.dec.deg,
    label=f"both fhd & pyfhd ({fhd_in_pyfhd.sum()})",
)
ax.scatter(
    fhd_sky.ra[fhd_unique_inds].wrap_at("180d").deg,
    fhd_sky.dec[fhd_unique_inds].deg,
    label=f"fhd only ({fhd_unique_inds.size})",
)
ax.scatter(
    pyfhd_sky.ra[pyfhd_unique_inds].wrap_at("180d").deg,
    pyfhd_sky.dec[pyfhd_unique_inds].deg,
    label=f"pyfhd only ({pyfhd_unique_inds.size})",
)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.1), ncol=3)
ax.set_xlabel("RA (deg)")
ax.xaxis.set_inverted(True)
ax.set_ylabel("Dec (deg)")
plt.savefig(fname=plotfile)
plt.close(fig)

same_ids = np.intersect1d(fhd_sky.name.astype(int), pyfhd_sky.name.astype(int))

fhd_sort_idx = fhd_sky.name.astype(int).argsort()
fhd_same_inds = fhd_sort_idx[
    np.searchsorted(fhd_sky.name.astype(int), same_ids, sorter=fhd_sort_idx)
]

pyfhd_sort_idx = pyfhd_sky.name.astype(int).argsort()
pyfhd_same_inds = pyfhd_sort_idx[
    np.searchsorted(pyfhd_sky.name.astype(int), same_ids, sorter=pyfhd_sort_idx)
]

np.testing.assert_allclose(
    fhd_sky.ra[fhd_same_inds].deg, pyfhd_sky.ra[pyfhd_same_inds].deg
)
np.testing.assert_allclose(
    fhd_sky.dec[fhd_same_inds].deg, pyfhd_sky.dec[pyfhd_same_inds].deg
)

diff_x = (
    fhd_sky.extra_columns["image_x"][fhd_same_inds]
    - pyfhd_sky.extra_columns["image_x"][pyfhd_same_inds]
)
diff_y = (
    fhd_sky.extra_columns["image_y"][fhd_same_inds]
    - pyfhd_sky.extra_columns["image_y"][pyfhd_same_inds]
)

abs_diff = np.sqrt(diff_x**2 + diff_y**2)
beam_diff = (
    fhd_sky.extra_columns["beam_I"][fhd_same_inds]
    - pyfhd_sky.extra_columns["beam_I"][pyfhd_same_inds]
)

plotfile = f"fhd_pyfhd_catalog_compare_pixel_beam{idl_str}{ref_str}.png"

plots = [
    {
        "ax_inds": [0, 0],
        "value": diff_x,
        "title": "x pixel diff (FHD - pyfhd)",
        "colorbar_label": "fractional x pixel deviation",
        "cmap": "BrBG",
    },
    {
        "ax_inds": [0, 1],
        "value": diff_y,
        "title": "y pixel diff (FHD - pyfhd)",
        "colorbar_label": "fractional y pixel deviation",
        "cmap": "PuOr",
    },
    {
        "ax_inds": [1, 0],
        "value": abs_diff,
        "title": "pixel distance (FHD vs pyfhd)",
        "colorbar_label": "fractional pixel deviation",
        "cmap": "viridis",
    },
    {
        "ax_inds": [1, 1],
        "value": beam_diff,
        "title": "I Beam diff (FHD - pyfhd)",
        "colorbar_label": "fractional I beam deviation",
        "cmap": "viridis",
    },
]

fig, axes = plt.subplots(2, 2)

for pltdict in plots:
    ax = axes[pltdict["ax_inds"][0]][pltdict["ax_inds"][1]]
    sc = ax.scatter(
        fhd_sky.ra[fhd_same_inds].wrap_at("180d").deg,
        fhd_sky.dec[fhd_same_inds].deg,
        c=pltdict["value"],
        cmap=pltdict["cmap"],
        s=1.0,
    )
    ax.set_xlabel("RA (deg)")
    ax.xaxis.set_inverted(True)
    ax.set_ylabel("Dec (deg)")
    ax.set_title(pltdict["title"])
    plt.colorbar(sc)

plt.tight_layout()
plt.savefig(fname=plotfile)
plt.close(fig)

breakpoint()
