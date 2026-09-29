from datetime import date, datetime, time
from scipy import stats
from math import atan2, cos, log, log2, pi
from statistics import mean
import random
import adjustText
import os
import geopandas as gp
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
from IPython.display import Markdown, display
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from scipy.interpolate import interp1d
from scipy.stats import pearsonr
from shapely.geometry import LineString, Point,box

from sklearn.mixture import GaussianMixture
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics import homogeneity_score
import rasterio as rio              
from rasterio.mask import mask
import json as jsn   
import cmath
from tslearn.clustering import TimeSeriesKMeans
import seaborn as sns
from scipy.stats import linregress
class NotEnoughPointError(ValueError):
    pass

def signifStar(pValue):
    if np.isnan(pValue):
        sig="nan"
        return sig
    if pValue<0.05:
        if pValue<0.01:
            if pValue<0.001:
                sig="***"
            else:
                sig="**"
        else:
            sig="*"
    else:
        sig="ns"
    return sig


def mapSeasonIndex(seasonValue):
  corr={"Juin-Aout":2,
  "Septembre-Novembre":3,
  "Decembre-Fevrier":4,
  "Mars-Mai":1}
  return corr[seasonValue]
  

def halfHigherPointDist(gdf_cerf:gp.GeoDataFrame):
    from scipy.spatial.distance import pdist, squareform
    A=np.ndarray(shape=(gdf_cerf.shape[0],2),buffer=np.array([[pt.x,pt.y] for pt in gdf_cerf.geometry]))
    D=squareform(pdist(A))
    N = np.max(D)
    return N/2

def thirdHigherPointDist(gdf_cerf:gp.GeoDataFrame):
    from scipy.spatial.distance import pdist, squareform
    A=np.ndarray(shape=(gdf_cerf.shape[0],2),buffer=np.array([[pt.x,pt.y] for pt in gdf_cerf.geometry]))
    D=squareform(pdist(A))
    N = np.max(D)
    return N/3
def createInputNams(gdf_cerf):
    animalTest=gdf_cerf.prenom.unique()[0]
    tmp=pd.DataFrame({"x":gdf_cerf.geometry.x, "y":gdf_cerf.geometry.y, "time":range(0,gdf_cerf.shape[0])})

    resultPath=os.path.join(os.getcwd(), animalTest)
    try:
        os.mkdir(resultPath)
    except OSError as error:
        print(error)   
    tmp.to_csv(os.path.join(resultPath,f"{animalTest}_input_file_for_VFractal_nams.csv"), index=False,header=False)
from scipy.cluster.hierarchy import dendrogram
def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    distance = np.arange(model.children_.shape[0])

    linkage_matrix = np.column_stack(
        [model.children_, distance, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)


def reverseCol(col:pd.Series):
    newS=col.iloc[::-1]
    newS.reset_index(inplace=True,drop=True)
    return(newS)

def reverseTraj(gdf_cerf:gp.GeoDataFrame):
    newTraj=gdf_cerf["geometry"].iloc[::-1]
    newTraj.reset_index(inplace=True,drop=True)
    gdf_cerf["geometry"]=newTraj
    return gdf_cerf

def getRandomStart(trajData:gp.GeoDataFrame, pathProportion=0.4):
    nrow=trajData.shape[0]
    x=random.randint(0,round(nrow*pathProportion))
    newTraj=trajData.copy()
    newTraj=newTraj[x:].reset_index(drop=True)
    return newTraj

def interpolateGPSpath(GPSPath, pointMultiplier=50, fileName="test"):
    if not isinstance(GPSPath, gp.GeoDataFrame):
        print("input data are not GeoDataFrame")
        raise ValueError
    try:
        startDateTime=pd.Timestamp(GPSPath.dateTime.iloc[0])
        endDateTime=pd.Timestamp(GPSPath.dateTime.iloc[-1])
        sampleTime = pd.date_range(startDateTime.value,endDateTime.value, GPSPath.shape[0]*pointMultiplier)
    except KeyError as err:
        print(f""" {err} : some columns or column names are not found in the GeoDataFrame provided.
        Make sure you have a 'UTC_DATE'and 'UTC_TIME' column""")

    dateTimeInt = GPSPath.dateTime.values.astype(int)
    sDateTimeInt = sampleTime.astype(int)

    ## interpolation sur les données temporelles et GPS réelles
    LongInte = interp1d(dateTimeInt, GPSPath.geometry.x)
    LatInte = interp1d(dateTimeInt, GPSPath.geometry.y)

    ### génération du modèle interpolé sur un range de temps défini
    LongitudeInterpolation = LongInte(sDateTimeInt)
    LatitudeInterpolation = LatInte(sDateTimeInt)


    trajet_interp_demo=gp.GeoDataFrame({"dateTime":sDateTimeInt, "geometry":gp.points_from_xy(LongitudeInterpolation, LatitudeInterpolation),"dataType":"interpolatedData"}, crs="epsg:2056")
    #trajet_interp_demo.rename(columns={0:"dateTime", "geometry":"geometry"},inplace=True)
    obsData=gp.GeoDataFrame({"dateTime":dateTimeInt, "geometry":GPSPath.geometry, "dataType":"obsData"},crs="epsg:2056")

    finalDF=pd.concat([trajet_interp_demo[1:],obsData])
    print(finalDF.shape)
    finalDF.sort_values(by=["dateTime"], inplace=True, ignore_index=True)
    finalDF.reset_index(drop=True, inplace=True)
    trajet_interp_demo=finalDF.copy()
    del(finalDF)
    trajet_interp_demo.to_file(f"{fileName}_interpolation_{pointMultiplier}.shp")
    return trajet_interp_demo
    
def writeVFractalFile(gdf_cerf:gp.GeoDataFrame, workDir=""):
    prenom=str(gdf_cerf.prenom.unique()[0])
    annee=str(gdf_cerf.deerYear.unique()[0])
    tmp=pd.DataFrame({"x":gdf_cerf.x_coord, "y":gdf_cerf.y_coord,"time":gdf_cerf.index})

    resultPath=os.path.join(workDir,"VFRactal_from_Nams")
    csvPath=os.path.join(resultPath,f"{prenom}_{annee}_allFixes_formatted.csv")
    tmp.to_csv(os.path.join(resultPath,f"{prenom}_{annee}_input_for_VFractal.txt"), index=False)  

def angle(A, B, C, /):
    Ax, Ay = A[0]-B[0], A[1]-B[1]
    Cx, Cy = C[0]-B[0], C[1]-B[1]
    a = atan2(Ay, Ax)
    c = atan2(Cy, Cx)
    if a < 0: a += pi*2
    if c < 0: c += pi*2
    return (pi*2 + c - a) if a > c else (c - a)

def NetCos(Net, Cos,step):
    eq7denominator=log2(Cos+1)+1
    dByAngle=2/eq7denominator
    eq4denominator=log(Net/step)
    dByNet=log(2)/eq4denominator

    return mean([dByAngle, dByNet])


def computeVFractal(step, Trajectory, showPlot=False, getGeoDataFrame=False, equation_1996=False, nsims=100, correlationPlot=False, successive=False):

    Trajet=LineString([Point(i) for i in Trajectory.geometry] )
    #assert step <= Trajet.length/3, "step divider is too big"

    red = reDiscretizePoints(Trajectory,step)
    finalPoint=complex(real=round(Trajectory.iloc[-1].geometry.x), imag=round(Trajectory.iloc[-1].geometry.y))
    
    finalPointRow=pd.DataFrame({"coord_X":finalPoint.real, "coord_Y":finalPoint.imag, "polar":finalPoint, "indexJ":Trajectory.index[-1]}, index=[0])
    red=pd.concat([red,finalPointRow], ignore_index=True)

    pointList=[Point(x,y) for x,y in zip(red.coord_X, red.coord_Y)]
    if successive:
        if len(pointList)<8:
            raise NotEnoughPointError
    else:
        if len(pointList)<12:
            raise NotEnoughPointError
    indexJ=[i for i in red.indexJ]
    
    intersX=[point.x for point in pointList]
    intersY=[point.y for point in pointList]

    i=0
    pointID=[]
    cosineList=[]
    netList=[]
    fractalList=[]
    nearestPoint=[]
    nearestNextPoint=[]
    pointListStep=[]
    pointNplusOne=[]
    pointNplusTwo=[]
    while True:
        if i+2 >= len(pointList):
            break
        p1=[pointList[i].x, pointList[i].y]
        p2=[pointList[i+1].x, pointList[i+1].y]
        p3=[pointList[i+2].x, pointList[i+2].y]
        
        pointID.append(i)
        nearestPoint.append(indexJ[i])
        nearestNextPoint.append(indexJ[i+2])
        pointListStep.append(pointList[i])
        pointNplusOne.append(pointList[i+1])
        pointNplusTwo.append(pointList[i+2])
        ### quand l'animal ne change pas ou très peu de direction, Theta est très petit. Et inversément
        theta=angle(p1,p2,p3)
        if theta > pi:
            theta = theta-pi
        else:
            theta=pi-theta
        cosine=cos(theta)

        cosineList.append(cosine)
        netList.append(Point(p1).distance(Point(p3)))
        ## cosine ranges from 1 (ie turning angle of 0°) to -1 (ie turing angle of 180°). 0 is at 90°
        ### if turning angle is > 90 (i.e when the animals goes slighlty "backwards"), we consider a fractal value of 2.
        if cosine > 0:
              d_den=log2(cosine+1)
              d_den=d_den+1
              d = 2/d_den
        
        else:
            d = 2
        fractalList.append(d)
        if successive:
            i = i + 1
        else:
            i = i + 2
        
    if len(cosineList)<=2:
        mFractal = np.nan
        varFractal=np.nan
        mNet = np.nan
        varNet=np.nan
        corr = np.nan
        mCosine=np.nan
        varCosine=np.nan
        pValue=np.nan
    else:
        if equation_1996:
            ## Cos equation
            mCosine = mean(cosineList)
        
            ## Net equation
            mNet=mean(netList)
            mCosine = mean(cosineList)
            mFractal= NetCos(Net=mNet, Cos=mCosine, step=step)
            varFractal=0
            varNet = np.var(netList)
            varCosine = np.var(cosineList)
        else:
            mFractal = mean(fractalList)
            varFractal=np.var(fractalList)
            mNet = mean(netList)
            varNet = np.var(netList)
            mCosine = mean(cosineList)
            varCosine = np.var(cosineList)

        #corr = pearsonr(x=cosineList[:-1], y=cosineList[1:])[0]
        slope, intercept, r_value, p_value, std_err=linregress(x=cosineList[:-1], y=cosineList[1:])
        corr=r_value
        try:
            pValue = round(pearsonr(x=cosineList[:-1], y=cosineList[1:])[1],4)
        except:
            pValue=np.nan
    if correlationPlot:
        plt.scatter(cosineList[:-1], cosineList[1:])
        plt.xlim(-1.10,1.10)
        plt.ylim(-1.1,1.1)
        plt.axhline(y=0, color='r', linestyle='--')
        plt.axvline(x=0, color='r', linestyle='--')
    if showPlot:
        stepX=[point.x for point in pointListStep]
        stepY=[point.y for point in pointListStep]
        base=10

        fontprops = fm.FontProperties(size=base)

        cerf_bounds=[*Trajectory.total_bounds]
        xSize=cerf_bounds[2]-cerf_bounds[0] #10
        ySize=cerf_bounds[3]-cerf_bounds[1] # 30
        yRatio=xSize/ySize
        if yRatio<1:
            fig,ax=plt.subplots(figsize=(base, base/yRatio),facecolor="white")
        else:
            fig,ax = plt.subplots(figsize=(base*yRatio,base),facecolor="white")

        ax.plot(Trajectory.geometry.x, Trajectory.geometry.y,linewidth=3,linestyle="-", color="C0", label="Original path")

        ax.scatter(intersX,intersY, color="red",s=4)
        ax.plot(intersX,intersY,linestyle="-", linewidth=3,color="red", label="Step path")
        el = patches.Ellipse((2, -1), 0.5, 0.5)
        ax.annotate("START",xy=(Trajectory.iloc[0].geometry.x, Trajectory.iloc[0].geometry.y),    xytext=(-20, -10), textcoords='offset points',
            size=10,
            #bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="simple",
                            fc="black", ec="none",
                            patchB=el,
                            connectionstyle="arc3,rad=0.3")
        )
        ax.annotate("FINISH",xy=(Trajectory.iloc[-1].geometry.x, Trajectory.iloc[-1].geometry.y), xytext=(20, 10), textcoords='offset points',

            size=10,
            #bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="simple",
                            fc="black", ec="none",
                            patchB=el,
                            connectionstyle="arc3,rad=0.3")
        )
        #ax.plot(stepX, stepY, '--', color="lightgrey", label = "Net distance")
        scalebar = AnchoredSizeBar(ax.transData,
                           step,
                           f"{step} m",
                           'upper center', 
                           pad=0.5,
                           color='black',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
        plt.gca().add_artist(scalebar)
        plt.gca().add_artist(ax.legend(loc="best", fontsize=10))
        for i in range(len(intersX)):
            ax.annotate(i,(intersX[i], intersY[i]))
        display(Markdown('<p><strong>Mean Fractal Value :</strong> {}<br><strong>Mean Scale Value :</strong> {}<br><strong>Correlation :</strong> {}</p>'.format(mFractal,mNet,corr)))       
  
    if getGeoDataFrame:
       # print(f"{len(pointID)}, {len(pointList)},{len(nearestPoint)},{len(nearestNextPoint)}")
        finalDF = gp.GeoDataFrame({"pointID":pointID,
        "vertexID":[x+1 for x in pointID],
        "geometry":pointListStep,
        "pointNplusOne":pointNplusOne,
        "pointNplusTwo":pointNplusTwo,
        "nearestPoint":nearestPoint,
        "nearestNextPoint":nearestNextPoint,
        "fractalList":fractalList,
        "angleList":cosineList,
         "netValue":netList })
        return finalDF
    else:
        return step,mFractal,varFractal, mNet,varNet,mCosine, varCosine,corr,pValue



def showTrajectoryOnMap(subset_cerf:gp.GeoDataFrame, backgroundRaster=None, groupBy="prenom", savePlot=False, workDir=""):
    colorSeasons={1:"greenyellow",2:"orangered",3:"goldenrod",4:"lightskyblue"}
    if not isinstance(subset_cerf, gp.GeoDataFrame):
        print("input data are not GeoDataFrame")
        raise ValueError
    base=10 

    if backgroundRaster is not None:
        cerf_bounds=box(*subset_cerf.total_bounds)
        cerf_bounds=cerf_bounds.buffer(distance=4000,resolution=1).envelope
        geo_bounds=gp.GeoDataFrame({"geometry":cerf_bounds},index=[0], crs=subset_cerf.crs)
  
        habitat=rio.open(backgroundRaster)
        def getFeatures(gdf):
            return [jsn.loads(gdf.to_json())['features'][0]['geometry']]
            
        coords=getFeatures(geo_bounds)
        out_image,out_transform=mask(dataset=habitat,shapes=coords,crop=True)
        out_meta=habitat.meta
        out_meta.update({"driver": "GTiff",
                        "height": out_image.shape[1],
                        "width": out_image.shape[2],
                        "transform": out_transform
                        })
        with rio.open(f"cropped_{backgroundRaster}", "w", **out_meta) as dest:
            dest.write(out_image)
        clipped=rio.open(f"cropped_{backgroundRaster}")
        ### retransformation en dataframe pour simplifier la lecture du code #####
        fig,ax=plt.subplots(tight_layout=True, figsize=(base,base), facecolor="white")
        rio.plot.show(clipped, ax=ax,cmap="Greys")
    else:
        cerf_bounds=[*subset_cerf.total_bounds]
        xSize=cerf_bounds[2]-cerf_bounds[0] #10
        ySize=cerf_bounds[3]-cerf_bounds[1] # 30
        yRatio=xSize/ySize
        if yRatio<1:
            fig,ax=plt.subplots(figsize=(base, base/yRatio),facecolor="white")
        else:
            fig,ax =plt.subplots(figsize=(base*yRatio,base),facecolor="white")

    fontprops = fm.FontProperties(size=10)
    scalebar = AnchoredSizeBar(ax.transData,
                            5000,
                            '5 km',
                            'upper center', 
                            pad=0.5,
                            color='black',
                            size_vertical=1,
                            fontproperties=fontprops)
    plt.gca().add_artist(scalebar)
    ## group by seasons AND animal
    #groups = dt_plot.groupby(["prenom","saison"])
    ## group by animal only
        
    if groupBy== "saison":
        subset_cerf["saisonIndex"]=subset_cerf["saison"].apply(mapSeasonIndex)

        groups=subset_cerf.groupby(["saisonIndex"])
        for name, group in groups:
            ax.plot(np.array(group["geometry"].x), np.array(group["geometry"].y), marker="o", linestyle="-",linewidth=3, label=str(group["saison"].unique()[0]),alpha=0.4,color=colorSeasons[group.saisonIndex.unique()[0]])
        ax.set_title(str(subset_cerf.prenom.unique()[0]) +"_"+ str(subset_cerf.deerYear.unique()[0]) +" : observed path")

    else:
        groups = subset_cerf.groupby(["prenom"])
        for name, group in groups:
            ax.plot(group["geometry"].x, group["geometry"].y, linestyle="-",linewidth=3, label=name)

    el = patches.Ellipse((2, -1), 0.5, 0.5)
    ax.annotate("START",xy=(subset_cerf.iloc[0].geometry.x, subset_cerf.iloc[0].geometry.y),    xytext=(-20, -10), textcoords='offset points',
            size=10,
            #bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="simple",
                            fc="black", ec="none",
                            patchB=el,
                            connectionstyle="arc3,rad=0.3")
    )
    ax.annotate("FINISH",xy=(subset_cerf.iloc[-1].geometry.x, subset_cerf.iloc[-1].geometry.y), xytext=(20, 10), textcoords='offset points',
            size=10,
            #bbox=dict(boxstyle="round", fc="0.8"),
            arrowprops=dict(arrowstyle="simple",
                            fc="black", ec="none",
                            patchB=el,
                            connectionstyle="arc3,rad=0.3")
    )
                            
    ax.legend(loc="lower left",fontsize=8,framealpha=1,markerscale=2)

    if savePlot:
        resultPath=os.path.join(workDir, subset_cerf.prenom.unique()[0])
        try:
            os.mkdir(resultPath)
        except OSError as error:
            pass
        plt.savefig(os.path.join(resultPath, f"{subset_cerf.prenom.unique()[0]}_{subset_cerf.deerYear.unique()[0]}_observed_path.png" ), transparent=False)


def reDiscretizePoints(trajet:gp.GeoDataFrame, R):
    indexJ=[0]
    p=trajet.geometry.apply(lambda row : complex(real=row.x, imag=row.y))
    result=[p[0]]

    result.extend([complex() for i in range(127)])
    ## index : I, j, i k
    I=0 ##  enlever 1 par rapport au code R
    j=1 ## enlever 1  par rapport au code R
    while j<=len(p):
        k=None
        for i in range(j,len(p)):
            
        
            d=cmath.polar(p[i]-result[I])[0] ## Mod
            if d>=R:
                k=i
                break
        if k is None:
            break
        j=k
        indexJ.append(j)
        XI=result[I].real
        xk_1=p[k-1].real
        YI=result[I].imag
        yk_1=p[k-1].imag
        lambda1=cmath.polar(np.diff([p[k-1], p[k]])[0])[1]
        cos_l=cmath.cos(lambda1)
        sin_l=cmath.sin(lambda1)

        U=(XI-xk_1)*cos_l + (YI-yk_1)*sin_l
        V=(YI-yk_1)*cos_l - (XI-xk_1)*sin_l
        H=U+cmath.sqrt(abs(R**2 - V**2))
        XIp1=H*cos_l+xk_1
        YIp1=H*sin_l+yk_1
        if len(result)<=I+1:
            result.extend([complex() for i in range(len(result))])
        result[I+1]=complex(real=round(XIp1), imag=round(YIp1))
        I=I+1
    #fin=pd.DataFrame({"coord_X":[x.real for x in result[:I]], "coord_Y":[x.imag for x in result[:I]], "polar": result[:I], "indexJ":indexJ[:I]})
    fin=pd.DataFrame({"coord_X":[x.real for x in result[:I+1]], "coord_Y":[x.imag for x in result[:I+1]], "polar": result[:I+1], "indexJ":indexJ[:I+1]})
    

    return fin


    
def TrajFractalDimensionValues(trajet:gp.GeoDataFrame, stepSizes:list, adjustD=False):
    redLength=[]
    for stepSize in stepSizes:
        red=reDiscretizePoints(trajet, stepSize)
       
        npd=np.diff(red.iloc[:].polar).tolist()
        red_length=sum([cmath.polar(x)[0] for x in npd])
        if adjustD:
            lastRed=red.iloc[-1].polar
            lastTrj=complex(real=trajet.iloc[-1].geometry.x, imag=trajet.iloc[-1].geometry.y)
            adjLength=red_length+cmath.polar(lastRed-lastTrj)[0]
            redLength.append(adjLength)
        else:
            redLength.append(red_length)
        
    fin=pd.DataFrame({"stepSize":stepSizes, "length":redLength})
    return fin



def getBehaviourVector(selectedStep:int, 
                       trajData:gp.GeoDataFrame, 
                       showPlot=False,
                       correlationPlot=False, 
                       nClass=2, 
                       testParameter="net distance",
                         method="gaussian",
                         successive=False,
                         smooth=False,
                         cutoff=0):
    animalID=trajData["prenom"].unique()[0]
    color={-1:"black",0:"red", 1:"blue",2:"green",3:"orange",4:"pink",5:"magenta",6:"cyan",7:"yellow"}

    df_selected=computeVFractal(selectedStep, trajData, getGeoDataFrame=True,successive=successive)
    df_selected["diffNearest"]=df_selected["nearestNextPoint"]-df_selected["nearestPoint"]
    corrDF=pd.DataFrame({"angN":df_selected["angleList"][:-1].reset_index(drop=True),"angNplus1":df_selected["angleList"][1:].reset_index(drop=True),"diffNearest":stats.zscore(df_selected["diffNearest"][:-1],nan_policy="omit").reset_index(drop=True)})
    #complete=pd.DataFrame({"angN":[0], "angNplus1":[0]})
    

    #corrDF=pd.concat([corrDF,complete])
    if testParameter=="netDistance":
        X_full = np.array(df_selected["netValue"])
        #X_full=X_full.reshape(-1,1)
    elif testParameter=="angle":
        X_full = np.array(df_selected["angleList"])
        #X_full=X_full.reshape(-1,1)
    elif testParameter=="autocorrelation":
        X_full = np.array(corrDF[["angN","angNplus1"]])
    elif testParameter=="timeElapsed_correlation":
        X_full = np.array(corrDF["angN"]+corrDF["angNplus1"],corrDF["diffNearest"])        
    elif testParameter=="both":
        X_full = np.array(df_selected[["netValue","angleList"]])
    elif testParameter=="index_angle":
        X_full = np.array(df_selected[["pointID","angleList"]])
    elif testParameter=="timeElapsed_angle":
        X_full = np.array(df_selected[["nearestPoint","angleList"]])
    elif testParameter=="timeElapsed_netDistance":
        X_full = np.array(df_selected[["nearestPoint","netValue"]])
    elif testParameter=="timeElapsed_both":
        X_full = np.array(df_selected[["nearestPoint","angleList","netValue"]])
    elif testParameter=="timeDiff":
        X_full = np.array(df_selected[["diffNearest","angleList"]])
    elif testParameter=="full":

        X_full = np.array(df_selected[["nearestPoint","nearestNextPoint","angleList","netValue"]])
   
    else:
        raise ValueError("Invalid string for testParameter : choose between 'netDistance', 'angle', 'both', 'timeElapsed','timeElapsed_angle ,'diffTime'" )
    if smooth:
        
        kernel_size = round(df_selected.shape[0]/5)
        kernel = np.ones(kernel_size) / kernel_size
        data_convolved = np.convolve(X_full, kernel, mode='same')
        X_full=data_convolved.reshape(-1,1)

    if method=="gaussianMixture":
        if len(X_full.shape)==1:
            X_full=X_full.reshape(-1,1)
        y_pred = GaussianMixture(n_components=nClass, random_state=42).fit(X_full).predict(X_full)
    elif method=="timeSeriesKmeans":
        if len(X_full.shape)==1:
            X_full=X_full.reshape(-1,1)
            
        km = TimeSeriesKMeans(n_clusters=nClass, metric="softdtw",random_state=666)
        y_pred = km.fit_predict(X_full)
        corrDF["y_pred"]=y_pred
#        corrDF.to_csv(f"output_{round(selectedStep)}.csv")
    elif method=="kmeans":
        if len(X_full.shape)==1:
            X_full=X_full.reshape(-1,1)
        km = KMeans(n_clusters=nClass, random_state=666, n_init=100)
        y_pred = km.fit_predict(X_full)

    elif method=="hclust":
        from sklearn.cluster import DBSCAN
        if len(X_full.shape)==1:
            X_full=X_full.reshape(-1,1)
        km = AgglomerativeClustering(n_clusters=nClass,distance_threshold=None, linkage="ward",metric="euclidean")
        #km = DBSCAN(eps=0.15, min_samples=nClass)
        #km = AgglomerativeClustering(distance_threshold=0, n_clusters=None,linkage="complete",metric="euclidean")
        
        #plot_dendrogram(km.fit(X_full), truncate_mode="level", p=nClass)
        y_pred = km.fit_predict(X_full)        
    elif method=="deterministic":
        y_pred=[]
        for value in X_full:
            if value<cutoff:
                y_pred.append(0)
            else:
                y_pred.append(1)
    elif method=="none":
        pass
    else:
        raise ValueError("Invalid string for method : choose between 'gaussian' or  'timeSeries' ")
    
    if testParameter=="autocorrelation":
        y_pred=np.append(y_pred,-1)

    if testParameter=="timeElapsed_correlation":
        y_pred=np.append(y_pred,-1)


    df_selected["behaviour"]=y_pred

    if correlationPlot:
        plt.scatter(corrDF.angN, corrDF.angNplus1, color=[color[i] for i in y_pred[:-1]])    
        plt.xlim(-1.10,1.10)
        plt.ylim(-1.1,1.1)
        plt.axhline(y=0, color='r', linestyle='--')
        plt.axvline(x=0, color='r', linestyle='--')    
        plt.ylabel("cosine of the next turning angle")
        plt.xlabel("cosine of the current turning angle")
    if showPlot:
        #colorBehaviour=np.append(-1, df_selected.behaviour[:-1])
        colorBehaviour=np.array(df_selected.behaviour)
        # Create 4x4 Grid
        fig=plt.figure(tight_layout=True, figsize=(10,10), facecolor="white")
        fig.suptitle(f"Animal : {animalID} : testParameter = {testParameter}, method = {method}, nClass={nClass}",fontsize=14)
        gs = fig.add_gridspec(nrows=2, ncols=2)

        # Create Three Axes Objects

        ax1=fig.add_subplot(gs[0, 0])
        ax1.hist(df_selected.netValue, color="white",bins=20)
        for i in range(nClass):
            ax1.hist(df_selected.netValue[df_selected.behaviour==i],alpha=0.5, color=color[i])
            plt.xlabel("Net distance")


        ax2 = fig.add_subplot(gs[0, 1])
        ax2.scatter(df_selected["pointID"], df_selected["netValue"], color=[color[i] for i in colorBehaviour])
        ax2.plot(df_selected["pointID"], df_selected["netValue"], linestyle="--")
        plt.ylabel("Net distance")
        plt.xlabel("Successive step points")
        for i in range(len(df_selected)):
            plt.annotate(df_selected.loc[i,"pointID"], xy=(df_selected.loc[i,"pointID"],df_selected.loc[i,"netValue"]))


        ax3 = fig.add_subplot(gs[1, 0])
        ax3.hist(df_selected.angleList, color="white",bins=20)
        for i in range(nClass):
            ax3.hist(df_selected.angleList[df_selected.behaviour==i],alpha=0.5, color=color[i])
            plt.xlabel("Cosine of turning angle")

        ax4 = fig.add_subplot(gs[1, 1])
        ax4.scatter(df_selected["pointID"], df_selected["angleList"], color=[color[i] for i in colorBehaviour])
        ax4.plot(df_selected["pointID"], df_selected["angleList"], linestyle="--")
        plt.ylabel("Cosine of turning angle")
        plt.xlabel("Successive step points")
        for i in range(len(df_selected)):
            plt.annotate(df_selected.loc[i,"pointID"], xy=(df_selected.loc[i,"pointID"],df_selected.loc[i,"angleList"]))
    
    

        #ax3.scatter(trajData2.geometry.x, trajData2.geometry.y, color = [color[i] for i in trajData2.behaviour])
        #ax3.scatter(trajData.loc[trajData["dataType"]=="obsData"].geometry.x,trajData.loc[trajData["dataType"]=="obsData"].geometry.y, color="black")
        #ax3.scatter(df_selected.geometry.x, df_selected.geometry.y, color="black")
    return df_selected


################################################################################################################################################################################
##########################   OLD FUNCTIONS #####################################################################################################################################
################################################################################################################################################################################

def computeVFractal_oldMethod(step, interpTrajectory, showPlot=False, getGeoDataFrame=False, equation_1996=False, nsims=100, correlationPlot=False):
    Trajet=LineString([Point(i) for i in interpTrajectory.geometry] )
    assert step <= Trajet.length/3, "step divider is too big"
    pointList = []
    pointListIndex=[]
    startPoint = interpTrajectory.geometry.iloc[0,]

    pointList.append(startPoint)
    pointListIndex.append(0)
    n=0
    for p in interpTrajectory.geometry:
        if pointList[-1].distance(p) > step:
            pointList.append(p)
            pointListIndex.append(n)
        n=n+1
    intersX=[point.x for point in pointList]
    intersY=[point.y for point in pointList]

    i=0
    pointID=[]
    cosineList=[]
    netList=[]
    fractalList=[]
    pointListStep=[]
    pointListIndexFin=[]
    pointListNextIndexFin=[]
    while True:
        if i+2 >= len(pointList):
            break
        p1=[pointList[i].x, pointList[i].y]
        p2=[pointList[i+1].x, pointList[i+1].y]
        p3=[pointList[i+2].x, pointList[i+2].y]
        #pointListStep.append(i)
        
        pointID.append(i)
        pointListStep.append(pointList[i])
        pointListIndexFin.append(pointListIndex[i])
        pointListNextIndexFin.append(pointListIndex[i+2])

        
        ### quand l'animal ne change pas ou très peu de direction, Theta est très petit. Et inversément
        theta=angle(p1,p2,p3)
        if theta > pi:
            theta = theta-pi
        else:
            theta=pi-theta
        cosine=cos(theta)

        cosineList.append(cosine)
        netList.append(Point(p1).distance(Point(p3)))
        ## cosine ranges from 1 (ie turning angle of 0°) to -1 (ie turing angle of 180°). 0 is at 90°
        ### if turning angle is > 90 (i.e when the animals goes slighlty "backwards"), we consider a fractal value of 2.
        if cosine > 0:
              d_den=log2(cosine+1)
              d_den=d_den+1
              d = 2/d_den
        
        else:
            d = 2
        fractalList.append(d)
        i = i + 2
        
    if len(cosineList)<=2:
        mFractal = np.nan
        varFractal=np.nan
        mNet = np.nan
        varNet=np.nan
        corr = np.nan
        mCosine=np.nan
        varCosine=np.nan
    else:
        if equation_1996:
            ## Cos equation
            mCosine = mean(cosineList)
        
            ## Net equation
            mNet=mean(netList)
            mCosine = mean(cosineList)
            mFractal= NetCos(Net=mNet, Cos=mCosine, step=step)
            varFractal=0
            varNet = np.var(netList)
            varCosine = np.var(cosineList)
        else:
            mFractal = mean(fractalList)
            varFractal=np.var(fractalList)
            mNet = mean(netList)
            varNet = np.var(netList)
            mCosine = mean(cosineList)
            varCosine = np.var(cosineList)

        corr = pearsonr(x=cosineList[:-1], y=cosineList[1:])[0]    
    if correlationPlot:
        plt.scatter(cosineList[:-1], cosineList[1:])
    if showPlot:
        stepX=[point.x for point in pointListStep]
        stepY=[point.y for point in pointListStep]
        fontprops = fm.FontProperties(size=18)
        fig,ax = plt.subplots(figsize=(15,15))
        ax.plot(interpTrajectory.geometry.x, interpTrajectory.geometry.y, color="yellow", label="trajet interpolé")
        ax.scatter(intersX,intersY, color="red",s=4)
        ax.plot(intersX,intersY, color="red")

        ax.plot(stepX, stepY, '--', color="lightgrey", label = "Net")
        scalebar = AnchoredSizeBar(ax.transData,
                           step,
                           f"{step} m",
                           'upper center', 
                           pad=0.5,
                           color='black',
                           frameon=False,
                           size_vertical=1,
                           fontproperties=fontprops)
        plt.gca().add_artist(scalebar)
        plt.gca().add_artist(ax.legend(loc="best", fontsize=20))
        for i in range(len(intersX)):
            ax.annotate(i,(intersX[i], intersY[i]))
        display(Markdown('<p><strong>Mean Fractal Value :</strong> {}<br><strong>Mean Scale Value :</strong> {}<br><strong>Correlation :</strong> {}</p>'.format(mFractal,mNet,corr)))       
  
    if getGeoDataFrame:
        finalDF = gp.GeoDataFrame({"pointID":pointID,"geometry":pointListStep,"pointIndex":pointListIndexFin,"nextPointIndex":pointListNextIndexFin,"fractalList":fractalList,"angleList":cosineList, "netValue":netList })
        return finalDF
    else:
        return step,mFractal,varFractal, mNet,varNet,mCosine, varCosine,corr



def getBehaviourVector_oldMethod(selectedStep:int, trajData:gp.GeoDataFrame, showPlot=False, nClass=2, testParameter="net distance"):
    df_selected=computeVFractal_oldMethod(selectedStep, trajData, getGeoDataFrame=True)
    trajData2=trajData.copy()
    if testParameter=="net distance":
        X_full = np.array(df_selected["netValue"])
        X_full=X_full.reshape(-1,1)
    elif testParameter=="angle":
        X_full = np.array(df_selected["angleList"])
        X_full=X_full.reshape(-1,1)
    elif testParameter=="both":
        X_full = np.array(df_selected[["netValue","angleList"]])
    else:
        raise ValueError("Invalid string for testParameter : choose between 'net distance', 'angle' or 'both'" )
    
    X_train, X_test = train_test_split(X_full, test_size=0.2, random_state=42)
    y_pred = GaussianMixture(n_components=nClass, random_state=42).fit(X_full).predict(X_full)
    df_selected["behaviour_class"]=y_pred
    behaviourVector=[]
    for i in range(df_selected.shape[0]):
        segment=df_selected.nextPointIndex[i]-df_selected.pointIndex[i]
        vec=[df_selected.behaviour_class[i]]*segment
        behaviourVector.extend(vec)
        
    fillData = [-1]*int(trajData2.shape[0]-len(behaviourVector))
    behaviour=pd.Series(behaviourVector+fillData, dtype="Int64")
    trajData2["behaviour"]=behaviour
    #trajData2.dropna(subset = ["behaviour"], inplace=True)

    if showPlot:
        # Create 4x4 Grid
        fig=plt.figure(tight_layout=True, figsize=(15,15))
        gs = fig.add_gridspec(nrows=2, ncols=2)
        color={0:"red", 1:"blue",2:"green",3:"orange"}

        # Create Three Axes Objects
        ax1=fig.add_subplot(gs[0, 0])
        ax1.hist(df_selected.netValue, color="white",bins=20)
        for i in range(nClass):
            ax1.hist(df_selected.netValue[df_selected.behaviour_class==i],alpha=0.5, color=color[i])
        ax2 = fig.add_subplot(gs[0, 1])
        ax2.scatter(df_selected["pointID"], df_selected["netValue"], color=[color[i] for i in y_pred])
        ax2.plot(df_selected["pointID"], df_selected["netValue"], linestyle="--")
        for i in range(len(df_selected)):
            plt.annotate(df_selected.loc[i,"pointID"], xy=(df_selected.loc[i,"pointID"],df_selected.loc[i,"netValue"]))
        ax3 = fig.add_subplot(gs[1, 0])
        ax3.hist(df_selected.angleList, color="white",bins=20)
        for i in range(nClass):
            ax3.hist(df_selected.angleList[df_selected.behaviour_class==i],alpha=0.5, color=color[i])
        ax4 = fig.add_subplot(gs[1, 1])
        ax4.scatter(df_selected["pointID"], df_selected["angleList"], color=[color[i] for i in y_pred])
        ax4.plot(df_selected["pointID"], df_selected["angleList"], linestyle="--")
        for i in range(len(df_selected)):
            plt.annotate(df_selected.loc[i,"pointID"], xy=(df_selected.loc[i,"pointID"],df_selected.loc[i,"angleList"]))
    
    

        #ax3.scatter(trajData2.geometry.x, trajData2.geometry.y, color = [color[i] for i in trajData2.behaviour])
        #ax3.scatter(trajData.loc[trajData["dataType"]=="obsData"].geometry.x,trajData.loc[trajData["dataType"]=="obsData"].geometry.y, color="black")
        #ax3.scatter(df_selected.geometry.x, df_selected.geometry.y, color="black")
    return trajData2
def getVfractalFigures(VfractalResults, savePlot=False, printValues=False, saveValues=False,fromSimulation=False, workDir="", param=[0,0,0], scale="log"):
    fig,axes=plt.subplots(2,2,figsize=(20,12), facecolor="white",tight_layout=True)
    axes=axes.ravel()
    gs = fig.add_gridspec(nrows=2, ncols=2)
    fig.suptitle(f"Red deer '{VfractalResults.prenom.unique()[0]}' {VfractalResults.deerYear.unique()[0]} : Fractal statistics", fontsize=14)

    axes[0].plot(np.array(VfractalResults.sort_values(by="stepSize").stepSize), np.array(VfractalResults.sort_values(by="stepSize").meanFrac),"o-")
    if fromSimulation:
        axes[0].plot(np.array(VfractalResults.stepSize), np.array(VfractalResults.lowQuantFrac), linestyle="--",color="lightGray")
        axes[0].plot(np.array(VfractalResults.stepSize), np.array(VfractalResults.upQuantFrac), linestyle="--",color="lightGray")
    axes[0].set_xlabel("step divider size [m]")
    axes[0].set_ylabel("mean fractal value")
    if scale=="log":
        axes[0].set_xscale("log")
    axes[0].grid(axis="x", which="both")
    
    def appendText(figure=axes[0],x="stepSize", y="meanFrac", display="all"):
        texts=[]
        xPos=[]
        yPos=[]
        confValue=[]
        stepVectorIndex=[]
        if display == "all":
            for index,row in VfractalResults.iterrows():
                texts.append(figure.text(x=row[x], y=row[y], s=index))
                xPos.append(row[x])
                yPos.append(row[y])
            adjustText.adjust_text(texts,x=xPos, y=yPos,arrowprops=dict(arrowstyle='->', color='red'),ax=figure,lim=100)
        elif display == "confidence":
            upConfCorr95=np.quantile(VfractalResults[y],q=0.95)
            upConfCorr90=np.quantile(VfractalResults[y],q=0.90)
            for index,row in VfractalResults.iterrows():
                if row[y]>=upConfCorr95:
                    texts.append(figure.text(x=row[x], y=row[y], s=index, color="purple"))
                    xPos.append(row[x])
                    yPos.append(row[y])
                    confValue.append("95%")
                    stepVectorIndex.append(index)
                elif upConfCorr90<=row[y]<upConfCorr95:
                    texts.append(figure.text(x=row[x], y=row[y], s=index, color="green"))
                    xPos.append(row[x])
                    yPos.append(row[y])
                    confValue.append("90%")
                    stepVectorIndex.append(index)

                else:
                    pass
               
            if texts:
                finDF=VfractalResults.iloc[stepVectorIndex,[2,5,9]].reset_index(drop=True)
                finDF.insert(loc=0, column="stepVectorIndex", value=stepVectorIndex)
                finDF["concerns"]=y
                finDF["value is above"]=confValue
                finDF.sort_values(by=y, axis=0, inplace=True)
                adjustText.adjust_text(texts,x=xPos, y=yPos,arrowprops=dict(arrowstyle='->', color='red'),ax=figure,lim=100)
                figure.axhline(y=upConfCorr90, color='green', linestyle='dotted')
                figure.annotate("90%", xy=(figure.get_xlim()[1],upConfCorr90),color="green")
                figure.axhline(y=upConfCorr95, color='purple', linestyle='dotted')
                figure.annotate("95%", xy=(figure.get_xlim()[1],upConfCorr95),color="purple")
                return finDF

            

            
        else:
            print(f"cannot interpret '{display}' as a method. please provide either 'all' or 'confidence' ")
            raise ValueError



    fracDF=appendText(figure=axes[0],x="stepSize", y="meanFrac", display="confidence")



    axes[1].plot(np.array(VfractalResults.sort_values(by="stepSize").stepSize), np.array(VfractalResults.sort_values(by="stepSize").varFrac),"o-")
    if fromSimulation:
        axes[1].plot(np.array(VfractalResults.stepSize), np.array(VfractalResults.lowQuantVarFrac), linestyle="--",color="lightGray")
        axes[1].plot(np.array(VfractalResults.stepSize), np.array(VfractalResults.upQuantVarFrac), linestyle="--",color="lightGray")
    axes[1].set_xlabel("step divider size [m]")
    axes[1].set_ylabel("variance of fractal value")
    if scale=="log":
        axes[1].set_xscale("log")
    axes[1].grid(axis="x", which="both")
    #appendText(figure=axes[1],x="stepSize", y="varFrac", display="all")


    #### bottom left plot

    axes[2].set_ylim(-1,1)
    axes[2].plot(np.array(VfractalResults.sort_values(by="stepSize").stepSize), np.array(VfractalResults.sort_values(by="stepSize")["autocorrelation"]),"o-")
    if fromSimulation:
        axes[2].plot(np.array(VfractalResults.stepSize), np.array(VfractalResults.lowQuantCorr), linestyle="--",color="lightGray")
        axes[2].plot(np.array(VfractalResults.stepSize), np.array(VfractalResults.upQuantCorr), linestyle="--",color="lightGray")
    if scale=="log":
        axes[2].set_xscale("log")
    axes[2].set_xlabel("step divider size [m]")
    axes[2].set_ylabel("correlation")
    axes[2].axhline(y=0, color='r', linestyle='--')
    axes[2].grid(axis="x", which="both")


    autDF=appendText(figure=axes[2], x="stepSize", y="autocorrelation", display="confidence")
  



    #### bottom right plot
    VfractalResults["zs"]=stats.zscore(VfractalResults.varFrac,nan_policy="omit")
    axes[3].scatter(VfractalResults.zs,VfractalResults.autocorrelation)
    axes[3].set_xlabel("standardized Variance")
    axes[3].set_ylabel("correlation")
    axes[3].axhline(y=0, color='r', linestyle='--')
    axes[3].axvline(x=0, color='r', linestyle='--')

    appendText(figure=axes[3], x="zs", y="autocorrelation", display="confidence")

    resultPath=os.path.join(workDir, VfractalResults.prenom.unique()[0])
    try:
        os.mkdir(resultPath)
    except OSError as error:
        pass
    if savePlot:

        if fromSimulation:
            plt.savefig(os.path.join(resultPath,f"{VfractalResults.prenom.unique()[0]}_{VfractalResults.deerYear.unique()[0]}_Fractal plots_{param[0]}_{param[1]}_{param[2]}_from_simulation.png"), transparent=False)
        else:
            plt.savefig(os.path.join(resultPath,f"{VfractalResults.prenom.unique()[0]}_{VfractalResults.deerYear.unique()[0]}_Fractal plots_{param[0]}_{param[1]}_{param[2]}_scales_observed_path.png"), transparent=False)
    if saveValues:
        finDF=pd.concat([autDF, fracDF])
        f=open(os.path.join(resultPath,f"{VfractalResults.prenom.unique()[0]}_{VfractalResults.deerYear.unique()[0]}_zeros.txt"),"a")
        if fromSimulation:
            f.write(f"**** {VfractalResults.prenom.unique()[0]}_{VfractalResults.deerYear.unique()[0]} from simulation ****\n\n")
            f.write(f"max value for autocorrelation : {VfractalResults.iloc[np.argmax(VfractalResults.autocorrelation)].stepSize}\n")
            f.write(f"max value for fractal : {VfractalResults.iloc[np.argmax(VfractalResults.meanFrac)].stepSize}\n\n")
            f.write(f"**** up quantiles table of fractal statistics ****\n")
            f.write(finDF.to_string())
            f.write('\n\n')
            f.close()
        else:
            f.write(f"**** {VfractalResults.prenom.unique()[0]}_{VfractalResults.deerYear.unique()[0]} from observed path ****\n\n")
            f.write(f"max value for autocorrelation : {VfractalResults.iloc[np.argmax(VfractalResults.autocorrelation)].stepSize}\n")
            f.write(f"max value for fractal : {VfractalResults.iloc[np.argmax(VfractalResults.meanFrac)].stepSize}\n\n")
            f.write(f"**** up quantiles table of fractal statistics ****\n")
            f.write(finDF.to_string())
            f.write("\n\n")
            f.close()


    if printValues:
        return pd.concat([autDF, fracDF])
   

        
        
        

def fractalMultiProcess(run,gdf_cerf,stepVector, randomProportion=0.10, fixedProportion=200, fixed=False):

    meanScale=[]
    varScale=[]
    meanFrac=[]
    varFrac=[]
    stepSize=[]
    meanCosine=[]
    varCosine=[]
    accCorr=[]
    pVal=[]
    #### start from a random GPS fix that is located at the beginning of the path (select the proportion)
    if run!=0:
        if fixed:
            randomStartTraj=gdf_cerf.iloc[fixedProportion*run:].reset_index(drop=True)
        else:
            randomStartTraj=getRandomStart(gdf_cerf,pathProportion=randomProportion)
    else:
        randomStartTraj=gdf_cerf
    #print(f"starting vFractal calculation for run {run+1}...")

    #####
    for absoluteDistance in stepVector:
        try:
            step,mFrac,vFrac, mScale,vScale,mCosine,vCosine,corr,pV = computeVFractal(absoluteDistance,randomStartTraj,successive=True,equation_1996=False )
            meanFrac.append(mFrac)
            varFrac.append(vFrac)
            meanScale.append(mScale)
            stepSize.append(step)
            varScale.append(vScale)
            accCorr.append(corr)
            meanCosine.append(mCosine)
            varCosine.append(vCosine)
            pVal.append(pV)
        except NotEnoughPointError as nepe:
            break
        except Exception as exc1:
            print(exc1)
    df_stat_sim=pd.DataFrame({"prenom":gdf_cerf.prenom.unique()[0], "deerYear": gdf_cerf.deerYear.unique()[0],"stepSize":stepSize,
    "meanScale":meanScale,
    "varScale":varScale,
    "meanFrac":meanFrac,
    "varFrac":varFrac,
    "meanCosine":meanCosine, 
    "varCosine":varCosine, 
    "autocorrelation":accCorr,
    "pValue":pVal})
    df_stat_sim["runNo"]=run
    print(f"run {run+1} complete")
    return df_stat_sim
        




###############################################################################
###############3 **** SHOW MULTIPLE TRAJECTORIES IN A SINGLE PLOT **** ########
#################################################################################




def showMultipleTrajectoriesOnMap(gdf_list, backgroundRaster=None, groupBy="prenom", savePlot=False, workDir="", plotNames=None):
    """
    Affiche plusieurs trajectoires sur le même graphique.
    gdf_list : Une liste de GeoDataFrames (un par individu/trajectoire) ou un unique GeoDataFrame contenant tout le monde.
    plotNames : Une liste de noms (str) correspondante à chaque GeoDataFrame pour le titre/sauvegarde (optionnel).
    """
    # 1. Harmonisation de l'entrée : si on passe un seul GeoDataFrame, on le met dans une liste
    if isinstance(gdf_list, gp.GeoDataFrame):
        gdf_list = [gdf_list]
    elif not isinstance(gdf_list, list) or not all(isinstance(x, gp.GeoDataFrame) for x in gdf_list):
        print("input data must be a GeoDataFrame or a list of GeoDataFrames")
        raise ValueError

    colorSeasons = {1: "greenyellow", 2: "orangered", 3: "goldenrod", 4: "lightskyblue"}
    base = 10 

    # 2. Calcul des limites globales combinant tous les GeoDataFrames
    all_bounds = np.array([gdf.total_bounds for gdf in gdf_list]) # [minx, miny, maxx, maxy]
    global_bounds = [
        all_bounds[:, 0].min(),  # minx global
        all_bounds[:, 1].min(),  # miny global
        all_bounds[:, 2].max(),  # maxx global
        all_bounds[:, 3].max()   # maxy global
    ]

    # 3. Gestion du fond de carte et de la figure
    if backgroundRaster is not None:
        from shapely.geometry import box
        import rasterio as rio
        from rasterio.mask import mask
        
        # Utilisation des limites globales avec le buffer
        cerf_bounds = box(*global_bounds)
        cerf_bounds = cerf_bounds.buffer(distance=4000, resolution=1).envelope
        geo_bounds = gp.GeoDataFrame({"geometry": cerf_bounds}, index=[0], crs=gdf_list[0].crs)
  
        habitat = rio.open(backgroundRaster)
        def getFeatures(gdf):
            return [jsn.loads(gdf.to_json())['features'][0]['geometry']]
            
        coords = getFeatures(geo_bounds)
        out_image, out_transform = mask(dataset=habitat, shapes=coords, crop=True)
        out_meta = habitat.meta
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })
        
        cropped_path = f"cropped_{os.path.basename(backgroundRaster)}"
        with rio.open(cropped_path, "w", **out_meta) as dest:
            dest.write(out_image)
        clipped = rio.open(cropped_path)
        
        fig, ax = plt.subplots(tight_layout=True, figsize=(base, base), facecolor="white")
        rio.plot.show(clipped, ax=ax, cmap="Greys")
    else:
        xSize = global_bounds[2] - global_bounds[0]
        ySize = global_bounds[3] - global_bounds[1]
        yRatio = xSize / ySize
        if yRatio < 1:
            fig, ax = plt.subplots(figsize=(base, base / yRatio), facecolor="white")
        else:
            fig, ax = plt.subplots(figsize=(base * yRatio, base), facecolor="white")

    # 4. Barre d'échelle
    fontprops = fm.FontProperties(size=10)
    scalebar = AnchoredSizeBar(ax.transData,
                                5000,
                                '5 km',
                                'upper center', 
                                pad=0.5,
                                color='black',
                                size_vertical=1,
                                fontproperties=fontprops)
    ax.add_artist(scalebar)

    # 5. Boucle pour tracer chaque trajectoire de la liste
    el = patches.Ellipse((2, -1), 0.5, 0.5)
    
    for i, subset_cerf in enumerate(gdf_list):
        # On définit un label ou un préfixe pour distinguer les individus si groupBy == "saison"
        label_prefix = f"{subset_cerf.prenom.unique()[0]} - " if len(gdf_list) > 1 else ""

        if groupBy == "saison":
            # Note: Assurez-vous que la fonction mapSeasonIndex est bien définie dans votre script
            subset_cerf["saisonIndex"] = subset_cerf["saison"].apply(mapSeasonIndex)
            groups = subset_cerf.groupby(["saisonIndex"])
            for name, group in groups:
                ax.plot(np.array(group["geometry"].x), np.array(group["geometry"].y), 
                        marker="o", linestyle="-", linewidth=3, 
                        label=f"{label_prefix}{str(group['saison'].unique()[0])}", 
                        alpha=0.4, color=colorSeasons[group.saisonIndex.unique()[0]])
        else:
            groups = subset_cerf.groupby(["prenom"])
            for name, group in groups:
                ax.plot(group["geometry"].x, group["geometry"].y, linestyle="-", linewidth=3, label=name)

        # Ajout des marqueurs START et FINISH pour chaque trajectoire
        ax.annotate(f"START ({subset_cerf.prenom.unique()[0]})", 
                    xy=(subset_cerf.iloc[0].geometry.x, subset_cerf.iloc[0].geometry.y),    
                    xytext=(-20, -10), textcoords='offset points', size=8,
                    arrowprops=dict(arrowstyle="simple", fc="black", ec="none", patchB=el, connectionstyle="arc3,rad=0.3"))
        
        ax.annotate(f"FINISH ({subset_cerf.prenom.unique()[0]})", 
                    xy=(subset_cerf.iloc[-1].geometry.x, subset_cerf.iloc[-1].geometry.y), 
                    xytext=(20, 10), textcoords='offset points', size=8,
                    arrowprops=dict(arrowstyle="simple", fc="black", ec="none", patchB=el, connectionstyle="arc3,rad=0.3"))

    # 6. Titre global
    if plotNames:
        ax.set_title(f"Observed paths: {', '.join(plotNames)}")
    else:
        all_names = [gdf.prenom.unique()[0] for gdf in gdf_list]
        ax.set_title(f"Observed paths: {', '.join(list(set(all_names)))}")
                            
    ax.legend(loc="lower left", fontsize=8, framealpha=1, markerscale=2)

    # 7. Sauvegarde du graphique combiné
    if savePlot:
        save_name = f"multiple_paths_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png"
        if plotNames:
            save_name = f"{'_'.join(plotNames)}_observed_path.png"
        
        resultPath = os.path.join(workDir, "multi_trajectories")
        os.makedirs(resultPath, exist_ok=True)
        plt.savefig(os.path.join(resultPath, save_name), transparent=False)
