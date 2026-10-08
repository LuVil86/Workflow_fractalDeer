#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13 10:10:45 2026

@author: loreto
"""

print("******* WX TO RESISTANCES v1.0 *******")
import rasterio as rio
import numpy as np


def computeResistance_keeley(expRaster, minExp, maxExp, keeleyFactor = 4):
    hsI=(expRaster-minExp)/(maxExp-minExp)    
    
    kterm=np.exp(-1*keeleyFactor*hsI, dtype=np.float32)
    del hsI
    res = 100-99*((1-kterm)/(1-np.exp(-keeleyFactor, dtype=np.float32)))
    del kterm
    return res


CHmask = rio.open('zip+file:///media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/Base data/mask_hermine_CH_2500_5m_BINARY.zip!mask_hermine_CH_2500_5m_BINARY.tif').read(1)



src = rio.open('/media/loreto/Grande/ie-ofev-24-25/sdm_hermine/wx_GLM_workshop2_AGGR_BCKG.tif')
wx = src.read(1)

truncation = True

if truncation:
    quantExp = np.nanquantile(wx[CHmask==1], q=[0.05,0.95])
    wx[wx>=quantExp[1]] = quantExp[1]
    wx[wx<=quantExp[0]] = quantExp[0]


minExp = np.nanmin(wx[CHmask==1])
maxExp = np.nanmax(wx[CHmask==1])

print(f" ## min value of wx :  {minExp} max value of wx :  {maxExp} ##")
print(f" quantile 5% :{quantExp[0]}, quantile  95%  : {quantExp[1]}")
resistance = computeResistance_keeley(wx, minExp, maxExp, 4)

out_meta = src.meta
out_meta.update({"driver": "GTiff",
                   "dtype":"int16",
                   "nodata": -9999
                   })


outFile = '/media/loreto/Grande/ie-ofev-24-25/sdm_hermine/res_raster_allS_cleaned/resistances_GLM_workshop2_AGGR_BCKG.tif'
print("--- writing output... ")
try:
    with rio.open(outFile,"w",**out_meta) as dst:
            dst.write(resistance,1)   
except Exception as exc:
    print("writing raster generated an exception : ", exc)
print(" >>>> SCRIPT COMPLETED")
