#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 18:35:16 2023

@author: loreto
"""

### This should calculate the sum of a lot of rasters as 
##long as they are all the same resolution, extent, CRS, 
##etc. I put in an assert statement to double check that.


import rasterio
import os

#cutTileFolder="/media/loreto/Linux/colo_tiles_1/cut_tiles"#projet coloplato
#cutTileFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_plateau/cut_tuile"#projet ie-ofev cerf-plateau
cutTileFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_plateau/cut_tuile_for_barriers_2_10000"#projet ie-ofev cerf-plateau #test for barrier values
#outFolder="/media/loreto/Linux/colo_tiles_1/summed_tiles"#projet coloplato
#outFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_plateau/sum_tuile"#projet ie-ofev cerf-plateau
outFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_plateau/sum_tuile_for_barriers_2_10000"#projet ie-ofev cerf-plateau #test for barrier values

if __name__=='__main__':
    print("******** SUM EAST-WEST AND NORTH-SOUTH TILES FROM CIRCUITSCAPE OUTPUTS *******\n")
    print(f"folder with tiles : {cutTileFolder}")
    print(f"output folder : {outFolder}")
    print("*"*100)
    tileList=[]
    for file in os.listdir(cutTileFolder):
        if file.endswith("ns_cum_curmap.tif"):
            tileList.append(int(file.split("_")[2]))
            
    for tileNumber in tileList:
        matchList=[]    
        for file in os.listdir(cutTileFolder):
            if file.find(f"_{tileNumber}_")!=-1:
                with rasterio.open(os.path.join(cutTileFolder,file)) as src:            
                    matchList.append(src.read())
                    result_profile=src.profile

        if len(matchList) == 2:
            try:
                with rasterio.open(os.path.join(outFolder, f't_{tileNumber}.tif'), 'w', **result_profile) as dst:
                    rasterFinal=matchList[0]+matchList[1]
                    dst.write(rasterFinal)            
                    print(f"tile {tileNumber} done")
            except:
                print(f"WARNING : there was an error writing summed tile N°{tileNumber} : skipping this tile...")
                continue
        else:
            print(f"WARHING : the number of CircuitScape tile results for the tile {tileNumber} is not 2 : skipping the tile")
            continue
        
