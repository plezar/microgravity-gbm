from spatialdata_io import xenium
import os

path = '/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/xenium_out'
out_path = '/Users/mzarodniuk/Documents/Scripts/microgravity-gbm/02_st-xenium/data/spatial_data_zarr'
dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]

for d in dirs:
    sdata = xenium(os.path.join(path, d))
    sdata.write(os.path.join(out_path, f'{d}.zarr'))