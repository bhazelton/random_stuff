from pyuvdata import UVData

file1 = "split_freq_1.uvh5"
file2 = "split_freq_2.uvh5"

uvd = UVData.from_file([file1, file2])
