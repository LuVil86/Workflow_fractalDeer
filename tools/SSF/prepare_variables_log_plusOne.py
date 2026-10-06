# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import rasterio as rio
from rasterio.windows import Window
from rasterio.transform import from_origin
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

#import raster
#raster_dir = "/media/loreto/Grande/ie-ofev-24-25/variables/human"
raster_dir = "/media/loreto/Grande/ie-ofev-24-25/variables/habitat"
#raster_dir = "/media/loreto/Grande/ie-ofev-24-25/variables/human/others"

DensVars=[
    #'Density_Buildings_50_opt2',
#           'Density_Buildings_100_opt2', 
#           "Density_Buildings_200_opt2", 
#           "Density_Buildings_400_opt2",
#          "Density_Merge_RoadPrimary__50_opt",
#          "Density_Merge_RoadPrimary__100_opt",
#          "Density_Merge_RoadPrimary__200_opt",
#          "Density_Merge_RoadPrimary__400_opt",
#          "Density_Merge_RoadSecondary__50_opt",
#          "Density_Merge_RoadSecondary__100_opt",
#          "Density_Merge_RoadSecondary__200_opt",
#        "Density_Merge_RoadSecondary_400_opt"]
"Density_Forest_50_interpol_opt",
           "Density_Forest_100_interpol_opt",
           "Density_Forest_200_interpol_opt",
          "Density_Forest_400_interpol_opt"]

# DistVars=["Dist_Merge_Bati_16b",
#           "Dist_Merge_RoadPrimary_16b",
#           "Dist_Merge_RoadSecondary_16b"
#           ]

for DensVar in DensVars:
#for DistVar in DistVars:
    # raster_input = f"{raster_dir}/{DistVar}.tif"
    # raster_output=f"{raster_dir}/log_distances/{DistVar}_log+1.tif"
    raster_input = f"{raster_dir}/{DensVar}.tif"
    raster_output=f"{raster_dir}/log_density/{DensVar}_log+1.tif" 
    with rio.open(raster_input) as src:
        arr=src.read(1)
        arr_log=np.log(arr+1)
        arr_log=arr_log.astype(np.float32)
        meta=src.meta
        meta.update(driver="GTiff",
                     dtype=arr_log.dtype
                     )
    with rio.Env(CHECK_DISK_FREE_SPACE="NO"):
        with rio.open(raster_output, "w", **meta) as dest:
            dest.write(arr_log,1)     
