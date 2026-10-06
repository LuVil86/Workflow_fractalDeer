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
#raster_dir = "/media/loreto/Grande/ie-ofev-24-25/variables/human/others"
# raster_dirIn = "/media/loreto/Grande/ie-ofev-24-25/variables/human/others/log_distances"
# raster_dirOut = "/media/loreto/Grande/ie-ofev-24-25/variables/human/others/scld_log_distances"
# raster_dirIn = "/media/loreto/Grande/ie-ofev-24-25/variables/human/log_density"
# raster_dirOut = "/media/loreto/Grande/ie-ofev-24-25/variables/human/scld_log_density"
# raster_dirIn = "/media/loreto/Grande/ie-ofev-24-25/variables/habitat/log_density"
# raster_dirOut = "/media/loreto/Grande/ie-ofev-24-25/variables/habitat/scld_log_density"
raster_dirIn = "/media/loreto/Grande/ie-ofev-24-25/variables/relief"
raster_dirOut = "/media/loreto/Grande/ie-ofev-24-25/variables/relief/scaled"
# Vars=['Dist_Merge_Bati_16b_log+1', 
#           'Dist_Merge_RoadPrimary_16b_log+1',
#           'Dist_Merge_RoadSecondary_16b_log+1'
#           ]

# Vars=['Density_Buildings_50_opt2_log+1',
#       'Density_Buildings_100_opt2_log+1', 
#       "Density_Buildings_200_opt2_log+1", 
#       "Density_Buildings_400_opt2_log+1",
#       "Density_Merge_RoadPrimary__50_opt_log+1",
#       "Density_Merge_RoadPrimary__100_opt_log+1",
#       "Density_Merge_RoadPrimary__200_opt_log+1",
#       "Density_Merge_RoadPrimary__400_opt_log+1",
#       "Density_Merge_RoadSecondary__50_opt_log+1",
#       "Density_Merge_RoadSecondary__100_opt_log+1",
#       "Density_Merge_RoadSecondary__200_opt_log+1",
#       "Density_Merge_RoadSecondary_400_opt_log+1"
#       ]

# Vars = ["Density_Forest_50_interpol_opt_log+1",
#         "Density_Forest_100_interpol_opt_log+1",
#         "Density_Forest_200_interpol_opt_log+1",
#         "Density_Forest_400_interpol_opt_log+1"]

Vars = ["Altitude_5m_16b",
        "Exposition_5m_16b",
        "Slope_5m_8b"
        ]

for Var in Vars: 
    raster_input = f"{raster_dirIn}/{Var}.tif"
    raster_output=f"{raster_dirOut}/{Var}_scaled.tif"
    with rio.open(raster_input) as src:
        arr=src.read(1)
        arr_mean=np.mean(arr)
        arr_std=np.std(arr)
        arr_scaled1=(arr-arr_mean)/arr_std
        arr_scaled2=arr_scaled1.astype(np.float32)
        meta=src.meta
        meta.update(driver="GTiff",
                     dtype=arr_scaled2.dtype
                     )
        
    del arr_scaled1
    print(f"scaling for {raster_output} done !")
    with rio.open(raster_output, "w", **meta) as dest:
        dest.write(arr_scaled2,1)     
    del arr_scaled2
        