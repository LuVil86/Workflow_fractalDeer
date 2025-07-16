#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 12:10:17 2023

@author: lucas et loreto both autheurs are first autheurs
"""

#import libs
import rasterio
import os
import geopandas as gpd
from rasterio.mask import mask
#CScapeFolder="/media/loreto/Linux/colo_tiles_1/out_tiles" #projet coloplato
#CScapeFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_plateau/tuile"#ie-ofev:cerf-plateau
CScapeFolder="/home/loreto/Documents/tuile_0"

#shpFolder="/media/loreto/Linux/colo_tiles_1/grille_20000"#projet coloplato
#shpFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_plateau/grille_2000_plateau"#projet ie-ofev:cerf-plateau
shpFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_alps/grille_3000_alps"#projet ie-ofev:cerf-plateau #test for value for barriers

#outFolder="/media/loreto/Linux/colo_tiles_1/cut_tiles"#projet ie-ofev:cerf-plateau
#outFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_plateau/cut_tuile"#projet ie-ofev:cerf-plateau
outFolder="/media/loreto/Grande/ie-ofev-24-25/csc_cerf_alps/cut_tuile_alps_model0"#projet ie-ofev:cerf-plateau #test for value for barriers

if __name__=="__main__":
    print("******** CUT TILES FROM CIRCUITSCAPE OUTPUTS *******\n")
    print(f"folder with tiles : {CScapeFolder}")
    print(f"folder with shapefiles as extend : {shpFolder}")
    print(f"output folder : {outFolder}")
    print("*"*100)
    for file in os.listdir(CScapeFolder):
        if file.endswith(".asc"):
            tileNumber=file.split("_")[1]
            print(f"--> processing tile {tileNumber} ...")
            shapeName=f"id_{tileNumber}.shp"
            try:
                tileShape=gpd.read_file(os.path.join(shpFolder, shapeName))
            except:
                print("Corresponding shapefile not Found")

            print(f" corresponding shapefile : {shapeName}")
            coord=tileShape.geometry
            
            with rasterio.open(os.path.join(CScapeFolder,file)) as src:
                out_image, out_transform = mask(src,coord,crop=True)
                out_meta = src.meta
                out_meta.update({'driver':'GTiff', 
                                'height':out_image.shape[1],
                                'width':out_image.shape[2],
                                'transform':out_transform})
                
            try:
                with rasterio.open(os.path.join(outFolder, f"cut_{os.path.splitext(file)[0]}.tif"),'w',**out_meta) as dest:
                    dest.write(out_image)
                print(f"*** Tile {tileNumber} successfully cut and written ***")
            except:
                print(f"WARNING : there was an error writing cut tile N°{tileNumber} : skipping this tile...")
                continue
            
