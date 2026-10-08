#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  2 11:33:57 2025

@author: loreto
"""


import numpy as np
import pandas as pd
import rasterio as rio
from datetime import datetime
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
#assert np.__version__>=1.24

def computeRSF(enviCoeff, habRasterPath, habCoefs, x1,x2,y1,y2,k):
    #print(f" running resistance calculation for tile {k}")  
##### assign coefficients to habitat classes ########
    with rio.open(habRasterPath) as habPath:
        habVar = habPath.read(1,window=Window.from_slices((y1, y2+1), (x1, x2+1)))
        
        #habTerm = np.empty(habVar.shape, dtype=np.dtype('float32'))
        habTerm = np.zeros(habVar.shape, dtype=np.dtype('float32'))
        ###entry coefficients of landuses
        ##create a dictionary for corrrespondances from a csv
        data_dict = habCoefs.set_index('value').to_dict()['coef']
        #Calculate the habitat terms b*landuse
        for key, value in data_dict.items():
            habTerm[habVar==key] = value
### stock habitat coefficient raster within termList
    TermList=[]
    TermList.append(habTerm)

###### multiply each continuous variable raster with its respective coefficient 
###### --> add them to TermList    
    for covar in enviCoeff:
        with rio.open(covar) as src:
            arr=src.read(1,window=Window.from_slices((y1, y2+1), (x1, x2+1)))
            TermList.append((arr/1000)*enviCoeff[covar])
     
    #### stack everything #######   
    TermStack=np.stack(TermList, axis=0)
    sumTerms = np.sum(TermStack, axis=0, dtype=np.float32)
    expsumTerms = np.exp(sumTerms, dtype=np.float32)
    return expsumTerms

def computeResistance(expRaster, minExp, maxExp):
    
    hsI=(expRaster-minExp)/(maxExp-minExp)    ### transform betweeen 0 and 1
    kterm=np.exp(-1*4*hsI, dtype=np.float32)
    res = 100-99*((1-kterm)/(1-np.exp(-4, dtype=np.float32)))
    #return res
    return expRaster
    

if __name__=="__main__":
    start = datetime.now()

    ################################ *** input parameters *** ###########################################
        
    habRasterPath='/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Habitats/Hermine/HabitatMap_v1_50_hermine_260106_Nibble2.tif'
    
    habCoefs= pd.read_csv('/media/loreto/Grande/ie-ofev-24-25/sdm_hermine/res_raster_AggrBckg50m/unique_habitat_hermine_for_resistances_06janv26_coef_glm13_PA3.csv')
    
    enviCoeff={'/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Density_Buildings_100_opt2_scaled.tif':-0.39702,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Density_SPBherb_50_optc_scaled.tif':0.1991684,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Density_SPBlign_400_optc_scaled.tif':0.0993502,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Density_Forest_400_opt2_scaled.tif':-0.4321434,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Density_Grasslands_50_optc_scaled.tif':0.2414276,
               #'/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Density_interfaces_ligneux_cultures_50_opt_scaled.tif':-0.0127027,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Density_interfaces_ligneux_prairies_50_opt_scaled.tif':0.0807438,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Dist_Habitat_Lakes_16b_scaled.tif':-0.4596374,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Dist_Habitat_Ponds_32b_scaled.tif':-0.1890168,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Dist_Merge_Rivers_16b_scaled.tif':-0.3014708,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Dist_Merge_Bati_16b_scaled.tif':-0.505005,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Dist_Merge_Autobahn_16b_scaled.tif':-0.1406077,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Dist_Merge_RoadPrimary_16b_scaled.tif':-0.2557935,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Dist_Merge_RoadSecondary_250701_16b_scaled.tif':-0.7536786,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/LogScaled/CHmax2500__Dist_Merge_Paths_16b_scaled.tif':-0.399277,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/Scaled/CHmax2500__Altitude_5m_16b_scaled.tif':0.5758556,
               #'/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/Scaled/CHmax2500__Exposition_5m_16b_scaled.tif':0.0237715,
               '/media/loreto/NAS_DEVELOPPEMENT/IE_OFEV/Variables/Scaled/CHmax2500__Slope_5m_8b_scaled.tif':-1.0647476
               }
                     
                       
    tileSize=(6000,6000)
    nbProcessors = 10
        
        
        
    outFile = '/media/loreto/Grande/ie-ofev-24-25/sdm_hermine/res_raster_AggrBckg50m/output/wx_hermine_AggrBckg50m_glm13_PA3_06janv26.tif'
    
    ############################################################################################3
        
    finResults={}
    with rio.open(habRasterPath) as inp:
        out_meta=inp.meta
        ncol=inp.width
        nrow=inp.height
        inp_transform = inp.transform ### !! the "transform" data indicates coordinates of the "upper-left" corner !!
        cellSize=inp_transform[0]
        nbSplitY=int(np.ceil(nrow/tileSize[0]))
        nbSplitX=int(np.ceil(ncol/tileSize[1]))
               
        print("#"*10)
        print("input raster features : ")
        print(f"dtype : {inp.dtypes}")
        print(f"{nrow} row X {ncol} columns ")
        print(f"cell size : {cellSize}")
        print(f"'nodata' code : {inp.nodata}")
        
        print("number of tiles :")
        print(f"{nbSplitY} tiles per row X {nbSplitX} tiles per column")
    
    
            
        splitCoordsXStart=[ i[0] for i in np.array_split(np.arange(ncol),nbSplitX)] 
        splitCoordsXStop=[ i[-1] for i in np.array_split(np.arange(ncol),nbSplitX)]
        splitCoordsYStart=[ i[0] for i in np.array_split(np.arange(nrow),nbSplitY)]
        splitCoordsYStop=[ i[-1] for i in np.array_split(np.arange(nrow),nbSplitY)]
        tileIndex=[ (i,j) for i in range(len(splitCoordsYStart)) for j in range(len(splitCoordsXStart)) ]
        inp.close()
    
    with ProcessPoolExecutor(max_workers=nbProcessors) as executor:
            poolDF={}
            minExp=[]
            maxExp=[]
            
            for k in tileIndex:        
                    poolDF[executor.submit(computeRSF,
    									   enviCoeff,
                                           habRasterPath,
                                           habCoefs,
                                           splitCoordsXStart[k[1]],
                                           splitCoordsXStop[k[1]],
                                           splitCoordsYStart[k[0]],
                                           splitCoordsYStop[k[0]],
                                           k)]=k
                                           
            for future in as_completed(poolDF):
                try:
                   finResults[poolDF[future]]=future.result()
                   minExp.append(np.min(finResults[poolDF[future]]))
                   maxExp.append(np.max(finResults[poolDF[future]]))
                   print(f"tile {poolDF[future]} done")

                except Exception as exc:
                    print('%r generated an exception: %s' % (poolDF[future], exc))
    
    
    ######### ***GATHERING RESULTS **** ########
        
    tmpRow=[]
    smallerExp=min(minExp)
    biggerExp=max(maxExp)
    for i in range(len(splitCoordsYStart)):
            tmpRow.append(np.hstack(tuple([computeResistance(finResults[(i,j)], smallerExp,biggerExp) for j in range(len(splitCoordsXStart))])))
            
            print(f"column {i} stacking done")
            
    ### remove MP results #####
    finResults=None



    #print(f"{i.shape[0]} rows  : {i.shape[1]} columns")
   
    #finalMat=np.vstack(tuple([i for i in tmpRow])).astype(np.int16)
    finalMat=np.vstack(tuple([i for i in tmpRow]))
    
    print("row stacking done")
    print("  ##############  ARRAY DONE  ###########")
    print(f" -- final array size : {finalMat.shape[0]} rows X {finalMat.shape[1]} columns")
    print (f" total elapsed time : {datetime.now()-start}")
    print(f"smaller value of minExp : {smallerExp}, bigger value of maxExp : {biggerExp} ")
    
    # out_meta.update({"driver": "GTiff",
    #                  "height": finalMat.shape[0],
    #                  "width": finalMat.shape[1],
    #                  "dtype":"int16",
    #                  "nodata":-9999
    #                  })
   
    out_meta.update({"driver": "GTiff",
                     "height": finalMat.shape[0],
                     "width": finalMat.shape[1],
                     "dtype":"float32",
                     "nodata":np.nan
                     })

    
    print("-- writing raster... ")
    try:
        with rio.open(outFile,"w",**out_meta) as dst:
                dst.write(finalMat,1)   
    except Exception as exc:
        print("writing raster generated an exception : ", exc)
    print(" >>>> SCRIPT COMPLETED")
