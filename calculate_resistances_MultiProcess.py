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

def computeRSF(enviVars,enviCoefs ,enviRasterDir, habRasterPath, habCoefs, x1,x2,y1,y2):
        
   
    with rio.open(habRasterPath) as habPath:
        habVar = habPath.read(1,window=Window.from_slices((y1, y2+1), (x1, x2+1)))
        
        habTerm = np.empty(habVar.shape, dtype=np.dtype('float32'))
        ###entry coefficients of landuses
        ##create a dictionary for corrrespondances from a csv
        data_dict = habCoefs.set_index('value').to_dict()['coef']
        #Calculate the habitat terms b*landuse
        for key, value in data_dict.items():
            habTerm[habVar==key] = value
    a=0
    TermList=[]
    TermList.append(habTerm)
    for enviVar in enviVars: 
        enviVar_input = f"{enviRasterDir}/{enviVar}.tif"
        with rio.open(enviVar_input) as src:
            arr=src.read(1,window=Window.from_slices((y1, y2+1), (x1, x2+1)))
            TermList.append(arr*enviCoefs[a])
            a=a+1
            
    TermStack=np.stack(TermList, axis=0)

    sumTerms = np.sum(TermStack, axis=0, dtype=np.float32)
    expsumTerms = np.exp(sumTerms, dtype=np.float32)
    
    return expsumTerms

def computeResistance(expRaster, minExp, maxExp):
    hsI=(expRaster-minExp)/(maxExp-minExp)    
    kterm=np.exp(-1*4*hsI, dtype=np.float32)
    res = 100-99*((1-kterm)/(1-np.exp(-4, dtype=np.float32)))
    return res

if __name__=="__main__":
    start = datetime.now()
    ################################ *** input parameters *** ###########################################
    
    habRasterPath='/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/ssf_raster/input/HabitatMap_cerf_6mai25_for_resistances.tif'
    habCoefs= pd.read_csv('/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/ssf_raster/unique_habitat_cerf_for_resistance_coef.csv')
       
       
    enviRasterDir = '/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/ssf_raster/input'
    outFile = '/media/loreto/Grande/ie-ofev-24-25/ssf_plateau/ssf_raster/output/resistance_cerf_plateau_2juin25.tif'
       
    enviVars=['Density_Buildings_100_opt2_log+1_scaled',
                 'Density_Merge_RoadPrimary__100_opt_log+1_scaled',
                 'Density_Merge_RoadSecondary__50_opt_log+1_scaled',
                 'Dist_Merge_Bati_16b_log+1_scaled',
                 'Dist_Merge_RoadPrimary_16b_log+1_scaled',
                 'Density_Forest_200_opt2_log+1_scaled']
       
    enviCoefs=[-0.2747795, -0.1703838, -0.1507785, 0.2947973, -0.2489365, 0.8208416]#write coef for continuous variables by enviVars order
       
    tileSize=(6000,6000)
    nbProcessors = 10
    
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
        print(f"{nrow} row X {ncol} columns ")
        print(f"cell size : {cellSize}")
        print(f"'nodata' code : {inp.nodata}")

        print("number of tiles :")
        print(f"{nbSplitY} tiles per row X {nbSplitX} tiles per column")


        
        splitCoordsXStart=[ i[0] for i in np.array_split(np.arange(ncol),nbSplitX)]  ## number of splits in the x coordinates (so this is the "vertical cuts")
        splitCoordsXStop=[ i[-1] for i in np.array_split(np.arange(ncol),nbSplitX)]
        splitCoordsYStart=[ i[0] for i in np.array_split(np.arange(nrow),nbSplitY)]  ## number of splits in the x coordinates (so this is the "horizontal cuts")
        splitCoordsYStop=[ i[-1] for i in np.array_split(np.arange(nrow),nbSplitY)]
        tileIndex=[ (i,j) for i in range(len(splitCoordsYStart)) for j in range(len(splitCoordsXStart)) ]
    inp.close()
    with ProcessPoolExecutor(max_workers=nbProcessors) as executor:
        poolDF={}
        minExp=[]
        maxExp=[]

        for k in tileIndex:        
                poolDF[executor.submit(computeRSF,enviVars,
                                       enviCoefs,
                                       enviRasterDir,
                                       habRasterPath,
                                       habCoefs,
                                       splitCoordsXStart[k[1]],
                                       splitCoordsXStop[k[1]],
                                       splitCoordsYStart[k[0]],
                                       splitCoordsYStop[k[0]],
                                       )]=k
                                       
            
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



    for i in tmpRow:
        print(f"{i.shape[0]} rows  : {i.shape[1]} columns")

    finalMat=np.vstack(tuple([i for i in tmpRow])).astype(np.int16)
    print("row stacking done")
    print("  ##############  ARRAY DONE  ###########")
    print(f" -- final array size : {finalMat.shape[0]} rows X {finalMat.shape[1]} columns")
    print (f" total elapsed time : {datetime.now()-start}")
    
    print(f"smaller value of minExp : {min(minExp)}, bigger value of maxExp : {max(maxExp)} ")
    out_meta.update({"driver": "GTiff",
                                                "height": finalMat.shape[0],
                                                "width": finalMat.shape[1],
                                                "dtype":"int16",
                                                "nodata":-9999
                                                })

    print("-- writing raster... ")
    try:
        with rio.open(outFile,"w",**out_meta) as dst:
                dst.write(finalMat,1)   
    except Exception as exc:
        print("writing raster generated an exception : ", exc)
