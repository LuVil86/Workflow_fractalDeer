#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 11:41:12 2025

@author: lucas villard
"""

import rasterio
from rasterio.merge import merge
import glob
import os
from rasterio.crs import CRS
from rasterio.enums import Resampling
# *** INPUT PARAMETERS ***** #

tileFolder = "/home/loreto/Documents/ebedu/Resultats_LU/sum_tuile_c1_1"
outFolder = "/home/loreto/Documents/ebedu/Resultats_LU/pinch-points"
outFile = os.path.join(outFolder, "mosaic_c1_1.tif") 

# *************************** #

if __name__=='__main__':
    print("******** CREATE MOSAIC FROM TILES *******\n")
    print(f"folder with tiles : {tileFolder}")
    print(f"output file name : {outFile}")
    print("*"*50)
    
    if not os.path.exists(outFolder):
	    print("the specified output folder does not exists :: creating it")
	    os.makedirs(outFolder)
    print("running merge function...")
    mosaic, out_transform = merge(glob.glob(tileFolder+"/"+"*.tif"), method="sum", resampling=Resampling.bilinear)
    print("--> mosaic done ! writing output...")
    #finArray = np.squeeze(mosaic, axis=0)
    # Copy the metadata

    # Update the metadata
    out_meta ={"driver": "GTiff",
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                    "transform": out_transform,
                    "count" : 1,
                    "nodata" : -9999,
                    "dtype" : 'float32',
                    "crs": CRS.from_epsg(2056)
                    }
                    

    # Write the mosaic raster to disk
    try:
        with rasterio.open(outFile, "w", **out_meta) as dest:
            dest.write(mosaic[0,:,:],1)
    except:
        print("ERROR : there was an error writing mosaic... aborting script")
    print(" >>>> SCRIPT COMPLETED")
