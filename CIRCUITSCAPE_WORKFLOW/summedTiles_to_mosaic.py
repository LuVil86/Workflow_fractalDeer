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
import numpy as np
from rasterio.crs import CRS
# *** INPUT PARAMETERS ***** 

tilePath = "/home/luvil/test_cleanPathway_fillHoles/ZN_tiles/"

outFile = os.path.join('/home/luvil/test_cleanPathway_fillHoles/test_mosaic', "mosaic_test.tif") 


if __name__=='__main__':
    print("******** CREATE MOSAIC FROM TILES *******\n")
    print(f"folder with tiles : {tilePath}")
    print(f"output file name : {outFile}")
    print("*"*50)
    print("running merge function...")
    mosaic, out_transform = merge(glob.glob(tilePath+"*.tif"), method="max")
    print("--> mosaic done")
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
