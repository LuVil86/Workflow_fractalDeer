# -*- coding: utf-8 -*-

import numpy as np
import os
from datetime import datetime
from osgeo import gdal,osr
from scipy import ndimage as ndi

def labelPatches(habitatRasterPath,selectedMilieux, outRaster=os.path.join(os.getcwd(),'patches.tif')):
    start = datetime.now()

    
    tmp=gdal.Open(habitatRasterPath, gdal.GA_ReadOnly)
    matrice = np.array(tmp.GetRasterBand(1).ReadAsArray())
    print(f"size of the input raster : {matrice.shape[0]} , {matrice.shape[1]}")



    mat_milieu = np.zeros(matrice.shape, dtype='int32')       #I had to add dtype=int in order to have integer only

    for milieu in selectedMilieux:
        mat_milieu[matrice==milieu] = 1 

    #mat_milieu= ndi.binary_closing(mat_milieu,structure=ndi.generate_binary_structure(2,2))
    labeled_array, num_features = ndi.label(mat_milieu, structure=ndi.generate_binary_structure(2,2))  #generate_binary_structure(2,2) for queen analyse (diagonal)

    print("patches annotations done in: ", datetime.now()-start)

    #change 0 in no data
    labeled_array[labeled_array==0] = -999

    print("Numb of clusters: ",num_features)

    band=tmp.GetRasterBand(1)
    geotransform = tmp.GetGeoTransform()
    wkt = tmp.GetProjection()
    driver = gdal.GetDriverByName("GTiff")
    output_file=outRaster

    dst_ds = driver.Create(output_file,
                        band.XSize,
                        band.YSize,
                        1,gdal.GDT_Int32)   #gdal.GDT_Int16
    
    dst_ds.GetRasterBand(1).WriteArray(labeled_array)
    #setting nodata value
    dst_ds.GetRasterBand(1).SetNoDataValue(-999)
    #setting extension of output raster
    # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
    dst_ds.SetGeoTransform(geotransform)
    # setting spatial reference of output raster
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt)
    dst_ds.SetProjection( srs.ExportToWkt() )

    print("****** DONE *******")
    print (f" total elapsed time : {datetime.now()-start}")

if __name__=="__main__":

    ##### parameters
    inputRasterPath="/home/lucas/FractalDeer_project/forest_patches/test_clean_secondaryRoads_third.tif"
    outputRasterPath=""
    labelPatches(habitatRasterPath=inputRasterPath, selectedMilieux=[19,20,21], outRaster=outputRasterPath)
