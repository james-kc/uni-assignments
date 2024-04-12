from os import listdir
import numpy as np

data_filenames = listdir("PHY4005/data")

a = {filename:np.load(f"PHY4005/data/{filename}") for filename in data_filenames}

print(a[data_filenames[0]])