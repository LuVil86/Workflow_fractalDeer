# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import rasterio as rio
import numpy as np
import pandas as pd

####read generic raster
inputRaster = '/media/loreto/Grande/ie-ofev-24-25/variables/habitat_generique_07avril25/test_crapaudCo.tif'
outputRaster = '/media/loreto/Grande/ie-ofev-24-25/variables/habitat_generique_07avril25/test_crapaudCo_output.tif'
####read correspondace file
data = pd.read_csv('/media/loreto/Grande/ie-ofev-24-25/variables/habitat_cerf_15avr25/habitat_classification_cerf_15avril25.csv')

###create a dictionary for corrrespondances from a csv
data_dict = data.set_index('value').to_dict()['new_raster_value']

#reclass raster
with rio.open(inputRaster, 'r') as inp:
        out_meta=inp.meta
        ncol=inp.width
        nrow=inp.height
        inp_transform = inp.transform ### !! the "transform" data indicates coordinates of the "upper-left" corner !!
        cellSize=inp_transform[0]
        habRast=inp.read(1)
        out = np.empty(habRast.shape, dtype=np.dtype('uint8'))
        for key, value in data_dict.items():
            out[habRast==key] = value
        out_meta.update({"driver": "GTiff",
                         "height": out.shape[0],
                         "width": out.shape[1],
                         "dtype":"uint8",
                         #"nodata":-9999
                         })
        print("-- writing raster... ")
        try:
            with rio.open(outputRaster,"w",**out_meta) as dst:
                    dst.write(out,1)   
        except Exception as exc:
            print("writing raster generated an exception : ", exc)