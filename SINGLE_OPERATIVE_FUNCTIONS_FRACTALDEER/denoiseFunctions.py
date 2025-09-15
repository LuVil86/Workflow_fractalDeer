
import os
from math import atan2, pi,cos
import pandas as pd
import geopandas as gp
from shapely.geometry import Point
import pytz
utc=pytz.UTC

workDir=os.getcwd()
# Common imports

#### FONCTION CALCUL DE L'ANGLE THETA 

def angle(A, B, C, /):
    Ax, Ay = A[0]-B[0], A[1]-B[1]
    Cx, Cy = C[0]-B[0], C[1]-B[1]
    a = atan2(Ay, Ax)
    c = atan2(Cy, Cx)
    if a < 0: a += pi*2
    if c < 0: c += pi*2
    return (pi*2 + c - a) if a > c else (c - a)



def calculateNetAngle(Trajectory:gp.GeoDataFrame):
    pointList=[Point(i) for i in Trajectory.geometry]
    #pointList=Trajectory.geometry
    i=0
    angleList=[]
    scaleList=[]
    pointListStep=[]
    pointListIndexFin=[]
    pointListNextIndexFin=[]
    distNextPoint=[]
    while True:
        if i+2 >= len(pointList):
            break
        p1=[pointList[i].x, pointList[i].y]
        p2=[pointList[i+1].x, pointList[i+1].y]
        p3=[pointList[i+2].x, pointList[i+2].y]
        
        

        pointListStep.append(pointList[i])
        pointListIndexFin.append(i)
        pointListNextIndexFin.append(i+2)
        distNextPoint.append(Point(p1).distance(Point(p2)))
        
        ### quand l'animal ne change pas ou très peu de direction, Theta est très petit. Et inversément
        theta=angle(p1,p2,p3)
        if theta > pi:
            theta = theta-pi
        else:
            theta=pi-theta
        cosinus=cos(theta)
        angleList.append(cosinus)
        scaleList.append(Point(p1).distance(Point(p3)))
        i = i + 1
        

    finalDF = gp.GeoDataFrame({"pointIndex":pointListIndexFin,
    "geometry":pointListStep,
    "nextPointIndex":pointListNextIndexFin,
    "angle":angleList, 
    "netDisplacement":scaleList,
    "distNextPoint":distNextPoint})
    return finalDF
def calculateSpeed(Trajectory:gp.GeoDataFrame):
  pointList=[Point(i) for i in Trajectory.geometry]
  i=0
  speed=[]
  while True:
    if i+1>=len(pointList):
      break
    p1=pointList[i]
    p2=pointList[i+1]
    speed.append(p1.distance(p2))
    i=i+1
  speed.append(0)
  return speed

def detectHighDist(Trajectory:gp.GeoDataFrame, maxDist:float):
    pointList=[Point(i) for i in Trajectory.geometry]
    #pointList=Trajectory.geometry
    i=0
    turnover=[0]
    while True:
        if i+1>=len(pointList):
            break
        p1=pointList[i]
        p2=pointList[i+1]
        d=p1.distance(p2)
        if(d>maxDist):
            turnover.append(i)
        i=i+1
    turnover.append(len(pointList))
    return turnover
def denoiseData(angleNetDF:gp.GeoDataFrame):
    i=0
    netDisplacement=[]
    diffNet1=[]
    netDisplacementNextPoint=[]
    diffNet2=[]
    netDisplacementTwoPointsAfter=[]
    pointIndex=[]
    angleBetween=[]
    speedNextPoint=[]
    while True:
        if i+2 >= len(angleNetDF):
            break
        pointIndex.append(int(i))
        netDisplacement.append(angleNetDF.loc[i,"netDisplacement"])
        netDisplacementNextPoint.append(angleNetDF.loc[i+1,"netDisplacement"])
        netDisplacementTwoPointsAfter.append(angleNetDF.loc[i+2,"netDisplacement"])
        angleBetween.append(angleNetDF.loc[i+1,"angle"])
        speedNextPoint.append(angleNetDF.loc[i+1, "distNextPoint"])
        if angleNetDF.loc[i+1,"netDisplacement"]!=0:
                diffNet1.append(angleNetDF.loc[i,"netDisplacement"]/angleNetDF.loc[i+1,"netDisplacement"])
        else:
                diffNet1.append(0)

        if angleNetDF.loc[i+2,"netDisplacement"]!=0:
                diffNet2.append(angleNetDF.loc[i,"netDisplacement"]/angleNetDF.loc[i+2,"netDisplacement"])
        else:
                diffNet2.append(0)
         
        i=i+1       
                
    
    finalDF = pd.DataFrame({"pointIndex":pointIndex,
    "netCurrentPoint":netDisplacement,
    "ratioNetNextPoint":diffNet1, 
    "ratioNet2PointsAfter":diffNet2,
     "angleBetween":angleBetween,
     "speedNextPoint":speedNextPoint})
    return(finalDF)

        
def inBetween(x,multiplier):
    invertMultiplier = 1/multiplier
    if x>1:
        return x<multiplier
    else:
        return x>invertMultiplier

def removeHighDist(gdf, maxSpeed):
    test=calculateNetAngle(gdf)
    ludicrousSpeed=test.loc[(test.netDisplacement>maxSpeed)].reset_index(drop=True)
    if len(ludicrousSpeed)!=0:
        stillFound=True
        indexToRemove=[i+2 for i in ludicrousSpeed.pointIndex]

        gdf_denoised=gdf.drop(gdf.index[indexToRemove]).reset_index(drop=True)
        return stillFound,gdf_denoised
    else:
        stillFound=False
        return stillFound,gdf_denoised

def denoisePipeline(gdf,netSpeed=True,angle_sharpness=-0.6, netNext_multiplier=4, net2Points_Multiplier=1.5,filterIndex=1,maxSpeed=10000):
    
    test=calculateNetAngle(gdf)
    test2=denoiseData(test)
    ### first locate the points with very sharp turning angle
    sharp=test2.loc[(test2.angleBetween<angle_sharpness)].reset_index(drop=True)
    ## check the netDist between the problematic point and the point n-2
    if netSpeed:
        sharp2=sharp.loc[(sharp.netCurrentPoint>=maxSpeed)].reset_index(drop=True)
    else:
        sharp2=sharp.loc[(sharp.speedNextPoint>=maxSpeed)].reset_index(drop=True)

        ## then check ratio of NET2 to NET1
    #sharp2=sharp.loc[(sharp.ratioNetNextPoint>netNext_multiplier) ].reset_index(drop=True)
    #sharpHigh=sharp2.loc[(sharp.ratioNetNextPoint>netNext_multiplier)  ].reset_index(drop=True)
    sharpHigh=sharp2[sharp2.ratioNet2PointsAfter.apply(lambda x : inBetween(x,net2Points_Multiplier))].reset_index(drop=True)
    ## then check if the length of the net displacement to reach the point with sharp angle is high
    #sharpHigh=sharp2.loc[(sharp2.netCurrentPoint>=maxSpeed)].reset_index(drop=True)
    #sharp3=sharp2.loc[(sharp2.ratioNetNextPoint>netNext_multiplier) ].reset_index(drop=True)
    print(f"""
    #########################################################
    number of points on the original dataset : {len(gdf)}
    number of sharp angles detected : {len(sharp)}
    number of sharp angles with  netDist bigger than {maxSpeed} from previous point :{len(sharp2)}
    number of points with ratio vertex above {net2Points_Multiplier} : {len(sharpHigh)}
    ====> NUMBER OF POINTS REMOVED : {len(sharpHigh)} <=====
     #########################################################
     
    """)
    if len(sharpHigh)>0:
        indexToRemove=[i+2 for i in sharpHigh.pointIndex]
        print("indexes of removed points : ")
        print([i for i in indexToRemove])
        gdf_denoised=gdf.drop(gdf.index[indexToRemove]).reset_index()
        gdf_denoised.rename({"index":f"indFilt_{filterIndex}"},axis=1, inplace=True)
        print(f" --> number of points after denoising : {len(gdf_denoised)}")

    else:
        print("no problematic points were detected with the provided parameters")
        return [gdf,0]

    
    if "level_0" in gdf_denoised.columns:
        return [gdf_denoised.drop(["level_0"], axis=1), len(sharpHigh)]
    else:
        return [gdf_denoised,len(sharpHigh)]

