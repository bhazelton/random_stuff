import numpy as np

from pyuvdata import UVData

filename = "/Users/bryna/Projects/Physics/data_files/mwax_uvfits/ssins_phased_gains_2.0sec_80kHz/1321445344.uvfits"

uvd_full = UVData.from_file(filename)

uvd1 = uvd_full.select(freq_chans=np.arange(uvd_full.Nfreqs//2), inplace=False)
uvd2 = uvd_full.select(freq_chans=np.arange(uvd_full.Nfreqs//2, uvd_full.Nfreqs), inplace=False)

file1 = "split_freq_1.uvh5"
file2 = "split_freq_2.uvh5"

uvd1.write_uvh5(file1)
uvd2.write_uvh5(file2)
