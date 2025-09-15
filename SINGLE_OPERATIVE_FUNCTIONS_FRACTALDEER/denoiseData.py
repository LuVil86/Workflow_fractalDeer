
import geopandas as gpd
import pandas as pd
import denoiseFunctions
import sys
import os
import matplotlib.pyplot as plt
#workDir="/".join(sys.argv[0].split("/")[:-1])


#### infile and nbPass

workDir = "/home/luvil/IE-OFEV/new_cerfs_valais_2025/"
animalName="ID225"
deerYear="deerYear_2023-2024"

nbPass=1


### denoise parameters
sh_f1=-0.6
multiplier1_f1=4
multiplier2_f1=4
speed_f1=1000 ### > 2km de déplacement à l'heure
netSpeed_f1=True



if nbPass == 0:
   tmp=pd.read_csv(os.path.join(workDir,animalName, f"{animalName}_{deerYear}.csv" ))
else:
   tmp=pd.read_csv(os.path.join(workDir,animalName, f"{animalName}_{deerYear}_denoised_{nbPass}.csv" ))


gdf_cerf=gpd.GeoDataFrame(tmp, geometry=gpd.GeoSeries.from_wkt(tmp.geometry),crs=2056)




gdf_denoised=denoiseFunctions.denoisePipeline(gdf_cerf,
netSpeed=netSpeed_f1,
angle_sharpness=sh_f1,
net2Points_Multiplier=multiplier2_f1,
filterIndex=nbPass,
maxSpeed=speed_f1)

plt.figure(figsize=(20,10))
ax1=plt.subplot(1,2,1)
ax1.ticklabel_format(style='plain')
if nbPass == 0:
    plt.title(f"Données GPS Originales : { animalName}  : {gdf_cerf.crs}",fontsize=20)
else:
   plt.title(f"filtrage n° {nbPass} : {animalName}",fontsize=20)

plt.xlabel("Coordonnée X  [m]")
plt.ylabel("Coordonnée Y  [m]")
plt.scatter(gdf_cerf.geometry.x,gdf_cerf.geometry.y, color="black")
plt.plot(gdf_cerf.geometry.x,gdf_cerf.geometry.y,linestyle="--", color="dimgrey")


ax2=plt.subplot(1,2,2, sharex=ax1, sharey=ax1)
ax2.ticklabel_format(style='plain')
plt.title(f"filtrage n° {nbPass+1} : {animalName}",fontsize=20)
plt.xlabel("Coordonnée X [m]")
plt.ylabel("Coordonnée Y [m]")
plt.scatter(gdf_denoised.geometry.x,gdf_denoised.geometry.y, color="orange",label="coordonnées GPS")
plt.plot(gdf_denoised.geometry.x,gdf_denoised.geometry.y,linestyle="--", color="dimgrey",label="Trajet")

plt.plot([], [], ' ', label=f"angle minimum :{sh_f1}")
plt.plot([], [], ' ', label=f"maxDist : {speed_f1}")
plt.plot([], [], ' ', label=f"ratio : {multiplier2_f1}")
plt.legend(loc="best")
gdf_denoised.to_csv(os.path.join(workDir,animalName, f"{animalName}_{deerYear}_denoised_{nbPass+1}.csv"), index=False)


### plt.show halt the script until the figure is closed
plt.show()
#plt.savefig(os.path.join(workDir,animalName, f"{animalName}_{deerYear}._1st_denoising.png"))
