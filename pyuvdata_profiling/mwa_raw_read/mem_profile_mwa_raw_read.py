from glob import glob

from pyuvdata import UVData

folder = "/Users/bryna/Projects/Physics/data_files/mwax_raw/1442248336_941555_vis/"

filelist = glob(folder + "*fits")

uvd = UVData.from_file(filelist)
