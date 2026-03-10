from pyuvdata import UVData

file1 = "split_times_1.uvh5"
file2 = "split_times_2.uvh5"

uvd = UVData.from_file([file1, file2])
