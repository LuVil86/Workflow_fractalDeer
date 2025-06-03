
import numpy as np
import scipy.ndimage as ndi
import rasterio as rio
from datetime import datetime
from rasterio.windows import Window
from concurrent.futures import ProcessPoolExecutor, as_completed
#assert np.__version__>=1.24

def clean_pathway(inputRaster,pathway, milieux, filter_size = 3) :  
    matrix=inputRaster.read(1)
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
        los = tmp==milieu
        los2=np.array(los, dtype=int)
        los2 = ndi.binary_fill_holes(los)
        los2 = np.where(los2==1,True,False)
        finMatrix[los2] = milieu                                                    #in the matrice, where los = True, replace by code milieu
        finMatrix[finMatrix==0] = -999
        los=None
        los2=None
    return finMatrix.astype(np.int16)


def clean_pathway_2(inputRaster,pathway, milieux, filter_size = 3) :  
    matrix=inputRaster.read(1)
    toRemove=[pathway]
    toKeep=[*milieux]
    toCross=[*milieux,pathway]
    print(toRemove)
    print(toCross)
    mask1 = np.array([[elem in toRemove for elem in row] for row in matrix]) ## mask1 avec juste les chemins
    mask1 = np.where(mask1==True, 1,0)
    #tmp = ndi.binary_erosion(mask1, iterations=1) ## erosion des chemins
    mask2 = np.array([[elem in toKeep for elem in row] for row in matrix]) ## masque avec les forêts 
    mask2=np.where(mask2==True, 1,0)
    finMask=ndi.binary_propagation(input=mask2, mask=mask1)
    finMask2=ndi.binary_erosion(finMask,iterations=3,mask=mask1)
    matrix[finMask2==1] = milieux[0] 
    #finMask = ndi.binary_fill_holes(mask2)
    #tmp2 = ndi.binary_closing(mask2)
    #finMask=np.add(tmp,mask2)
    
    
    #mask1= np.where(mask1==1, True,False)
    #print(mask1)
   
    return matrix


if __name__=="__main__":
    start = datetime.now()
    ######### input parameters ########
    toRemove= 24
    milieux=[15]
    filterSize=3

    ### input and output
    inputRaster="/home/luvil/test_cleanPathway_fillHoles/HabitatMap_cerf_15avril25_clip.tif"
    outputRaster="/home/luvil/test_cleanPathway_fillHoles/results_raster.tif"


    with rio.open(inputRaster) as inp:
        out_meta=inp.meta
        ncol=inp.width
        nrow=inp.height
        inp_transform = inp.transform ### !! the "transform" data indicates coordinates of the "upper-left" corner !!
        cellSize=inp_transform[0]

        result = clean_pathway_2(inp, toRemove, milieux, filterSize )

    out_meta.update({"driver": "GTiff",
                                                "height": result.shape[0],
                                                "width": result.shape[1],
                                                "dtype":"int16",
                                                "nodata":-999
                                                })

    print("-- trying to write raster... ")
    try:
        with rio.open(outputRaster,"w",**out_meta) as dst:
                dst.write(result,1)   
    except Exception as exc:
        print("writing raster generated an exception : ", exc)
