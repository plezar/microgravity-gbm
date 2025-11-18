from spatialdata_io import xenium
import os

path = '/Users/mzarodniuk/Documents/Scripts/Alice_Xenium/data'
dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

for d in dirs:
    sdata = xenium(os.path.join(path, d))
    sdata.write(os.path.join(path, f'{d}.zarr'))


