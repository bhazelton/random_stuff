from pyuvdata import UVData

filename = "/Users/bryna/Projects/Physics/data_files/mwax_uvfits/ssins_phased_gains_2.0sec_80kHz/1321445344.uvfits"

uvd = UVData.from_file(filename)
