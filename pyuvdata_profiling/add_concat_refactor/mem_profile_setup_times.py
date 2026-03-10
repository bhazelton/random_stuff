import numpy as np

from pyuvdata import UVData

filename = "/Users/bryna/Projects/Physics/data_files/mwax_uvfits/ssins_phased_gains_2.0sec_80kHz/1321445344.uvfits"

uvd_full = UVData.from_file(filename)

times = np.unique(uvd_full.time_array)
uvd1 = uvd_full.select(times=times[0 : len(times) // 2], inplace=False)
uvd2 = uvd_full.select(times=times[len(times) // 2 :], inplace=False)

file1 = "split_times_1.uvh5"
file2 = "split_times_2.uvh5"

uvd1.write_uvh5(file1)
uvd2.write_uvh5(file2)
