
import geopandas as gpd
import pandas as pd
import denoiseFunctions
import sys
import os
workDir="/".join(sys.argv[0].split("/")[:-1])

animalName="20074"
deerYear="deerYear_2011-2012"
tmp=pd.read_csv(os.path.join(workDir,animalName, f"{animalName}_{deerYear}.csv" ))
gdf_cerf=gpd.GeoDataFrame(tmp, geometry=gpd.GeoSeries.from_wkt(tmp.geometry),crs=2056)


### parameters for denoising
sh_f1=-0.6
multiplier1_f1=4
multiplier2_f1=4
speed_f1=2000 ### > 2km de déplacement 
netSpeed_f1=True
gdf_denoised=denoiseFunctions.denoisePipeline(gdf_cerf,
netSpeed=netSpeed_f1,
angle_sharpness=sh_f1,
net2Points_Multiplier=multiplier2_f1,
filterIndex=1,
maxSpeed=speed_f1)

gdf_denoised.to_csv(os.path.join(workDir,animalName, f"{animalName}_{deerYear}_denoised.csv"), index=False)
