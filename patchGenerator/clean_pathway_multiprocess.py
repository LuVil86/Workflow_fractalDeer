
import numpy as np
import scipy.ndimage as ndi
import rasterio as rio
from datetime import datetime
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
#assert np.__version__>=1.24

def clean_pathway(inputRaster,x1,x2,y1,y2,pathway, milieux, filter_size = 3) : 
    with rio.open(inputRaster) as rasterBuffer:  
        matrix=rasterBuffer.read(1, window=Window.from_slices((y1, y2+1), (x1, x2+1)))
        finMatrix=np.zeros(matrix.shape) 
        #print(f"input tile size : {finMatrix.shape[0]} rows X {finMatrix.shape[1]} columns")
        for milieu in milieux:                       #matrix (matrix) = matrix of the habitats | pathway (integer) = code for pathway in the matrix | milieu (integer) = code for the milieu (habitat) in the matrix
            fs = filter_size
            add_mat = int((fs-1) / 2)                                                           #add_mat = 1 (size of the edge added arond the matrix)
            inflatedMat = np.zeros((matrix.shape[0]+add_mat*2,matrix.shape[1]+add_mat*2))    #create empty matrix (full of 0). nb line = nb line of "matrix"+ 2. nb column = nb columns of "matrix"+ 2 (need *2 to have a  line more at left/right and up/down)
            inflatedMat[add_mat:-add_mat,add_mat:-add_mat] = matrix                          #inlude "matrix" inside matrix_expensa. The size of the buffer around the matrix depend of "add_mat"   Example  [[0. 0. 0. 0. 0.]            
            mat_pathway = inflatedMat == pathway                                             #extract the pathways from matrix_expensa (where matrix_expensa = pathway code => True)                          [0. 1. 1. 1. 0.]
            mat_milieu = inflatedMat == milieu                                                #extract the milieu from matrix_expensa (where matrix_expensa = milieu code => True)                            [0. 1. 1. 1. 0.] 
            for i in range(add_mat,inflatedMat.shape[0]-add_mat) :                           #shape[0] = row                                                                                                  [0. 1. 1. 1. 0.]
                for j in range(add_mat,inflatedMat.shape[1]-add_mat):                        #shape[1] = column                                                                                               [0. 0. 0. 0. 0.]]
                    aux = inflatedMat[i-add_mat:i+add_mat+1, j-add_mat:j+add_mat+1]          #define the area around the point (+1 to have a 3x3 matrix with element i,j in the center)                                                                         
                    if mat_pathway[i,j] and milieu in aux:                                      #if pathway next to milieu
                        mat_milieu[i,j] = True                                                       #in the milieu matrix, what was count as pathway (False) is now count as milieu (True)
            inflatedMat[mat_milieu] = milieu                                                  #in the matrix_expensa, where mat_milieu = True, replace by code milieu

            tmp=inflatedMat[add_mat:-add_mat,add_mat:-add_mat]
            los = int(tmp==milieu)
            los2 = ndi.binary_fill_holes(los)
            los2 = los2[los2==1]
            finMatrix[los2] = milieu                                                    #in the matrice, where los = True, replace by code milieu
            finMatrix[finMatrix==0] = -999
            los=None
            los2=None
        return finMatrix.astype(np.int16)

if __name__=="__main__":
    start = datetime.now()
    ######### input parameters ########
    nbProcessors=8
    toRemove= 24
    milieux=[15,16,17]
    tileSize=(676,742)
    filterSize=3

    ### input and output
    inputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_15avril25_clip.tif"
    outputRaster="/home/luvil/test_cleanPathway_fillHoles/results_raster.tif"


    finResults={}
    with rio.open(inputRaster) as inp:
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
        if nbSplitY == 1 and nbSplitX ==1 :
            print("!! ERROR only one tile would be generated with the parameters you specified : use the non-multiprocess script instead")
            raise ValueError

        
        splitCoordsXStart=[ i[0] for i in np.array_split(np.arange(ncol),nbSplitX)]  ## number of splits in the x coordinates (so this is the "vertical cuts")
        splitCoordsXStop=[ i[-1] for i in np.array_split(np.arange(ncol),nbSplitX)]
        splitCoordsYStart=[ i[0] for i in np.array_split(np.arange(nrow),nbSplitY)]  ## number of splits in the x coordinates (so this is the "horizontal cuts")
        splitCoordsYStop=[ i[-1] for i in np.array_split(np.arange(nrow),nbSplitY)]
        tileIndex=[ (i,j) for i in range(len(splitCoordsYStart)) for j in range(len(splitCoordsXStart)) ]
    inp.close()
    with ProcessPoolExecutor(max_workers=nbProcessors) as executor:
        poolDF={}


        for k in tileIndex:        
                poolDF[executor.submit(clean_pathway,inputRaster,splitCoordsXStart[k[1]],
                                       splitCoordsXStop[k[1]],
                                       splitCoordsYStart[k[0]],
                                       splitCoordsYStop[k[0]],
                                       toRemove,
                                       milieux,
                                       filterSize)]=k
                                       
            
        for future in as_completed(poolDF):
            
            try:
               
               finResults[poolDF[future]]=future.result()
               print(f"tile {poolDF[future]} done")
            except Exception as exc:
                print('%r generated an exception: %s' % (poolDF[future], exc))

 #   for i in range(len(splitCoordsYStart)):
 #       for j in range(len(splitCoordsXStart)):
 #           print(f"tile ({i}, {j}) : size = {finResults[(i,j)].shape[0]} rows X {finResults[(i,j)].shape[1]} columns" )
    
    tmpRow=[]
    for i in range(len(splitCoordsYStart)):
         #print("entering stack command..")
        # if i == 0:
        #    finalMat=np.hstack(tuple([finResults[(i,j)] for j in range(len(splitCoordsXStart))]))
        #    print("first row done")
        # else:
        #     tmp=np.hstack(tuple([finResults[(i,j)] for j in range(len(splitCoordsXStart))]))
        #     finalMat=np.vstack((finalMat,tmp ))
         #     tmp=None
        #     print(f"row {i+1} done")
        tmpRow.append(np.hstack(tuple([finResults[(i,j)] for j in range(len(splitCoordsXStart))])))
        
        print("column stacking done")
    finResults=None
    print(tmpRow[0].dtype)

    for i in tmpRow:
        print(f"{i.shape[0]} rows  : {i.shape[1]} columns")

    finalMat=np.vstack(tuple([i for i in tmpRow])).astype(np.int16)
    print("row stacking done")

    print("  ##############  ARRAY DONE  ###########")
    print (f" total elapsed time : {datetime.now()-start}")
    out_meta.update({"driver": "GTiff",
                                                "height": finalMat.shape[0],
                                                "width": finalMat.shape[1],
                                                "dtype":"int16",
                                                "nodata":-999
                                                })

    print("-- trying to write raster... ")
    try:
        with rio.open(outputRaster,"w",**out_meta) as dst:
                dst.write(finalMat,1)   
    except Exception as exc:
        print("writing raster generated an exception : ", exc)
