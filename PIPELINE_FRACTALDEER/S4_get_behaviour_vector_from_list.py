
import pandas as pd
import os
import sys
import fractalFunctions
import geopandas as gp
import numpy as np
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed, ThreadPoolExecutor

def launcher(inFilePath,outDir,prenom,deerYear, stepSize, check_existing_file):
    if(isinstance(prenom, int)):
        prenom=str(prenom)
    if stepSize=="NA":
        print("no stepSize specified : aborting script")
        return
    
    outFile=os.path.join(outDir,f"{prenom}_{deerYear}_stepSize_{round(stepSize)}_autocorrelation_timeSeriesKmeans_2_classes.csv")
    stepOutFile=os.path.join(outDir,f"DIVIDER_PATH_{prenom}_{deerYear}_stepSize_{round(stepSize)}_autocorrelation_timeSeriesKmeans_2_classes.csv")

    print(f"""==> computing selected step size : {stepSize}...
          """)
    if check_existing_file:
        if os.path.exists(outFile):
            print(f"The file {outFile} already exists ! skipping the computation...")
            return
    else:
        try:
            tmp=pd.read_csv(inFilePath)
            subset_cerf=gp.GeoDataFrame(tmp, geometry=gp.GeoSeries.from_wkt(tmp.geometry),crs=2056)
        
        except FileNotFoundError as fnf:
            print(f"the input GPS dataframe {inFilePath} has not been found.. did you forget to generate it ?")
        
    
    behaviourDF = fractalFunctions.getBehaviourVector(selectedStep=stepSize,trajData=subset_cerf,showPlot=False,
    nClass=2,
    testParameter="autocorrelation", 
    method="timeSeriesKmeans",
    successive=True,
    smooth=False,
    cutoff=0)

    #### infer behaviour from stepSize path to GPS fixes
    behaviourVector=[-1]*subset_cerf.shape[0]
    colorBehaviour=np.append(-1, behaviourDF.behaviour)  

    for i in range(behaviourDF.shape[0]-1):
        indexRange=[behaviourDF.iloc[i].nearestPoint, behaviourDF.iloc[i].nearestNextPoint]
        segment=behaviourDF.iloc[i].nearestNextPoint-behaviourDF.iloc[i].nearestPoint
        vec=[colorBehaviour[i]]*segment
        behaviourVector = behaviourVector[:indexRange[0]]+vec+behaviourVector[indexRange[1]:]
        # print(f"index {indexRange[0]} to index {indexRange[1]}, {segment} points will have the behaviour : {colorBehaviour[i]}")
        color={-1:"black",0:"red", 1:"blue",2:"green",3:"orange",4:"pink",5:"magenta",6:"cyan",7:"yellow"}
    pathNo=[]
    c=0
    for i in range(len(behaviourVector)-1):
        if behaviourVector[i]!=behaviourVector[i+1]:
                pathNo.append(c)
                c=c+1

        else:
            pathNo.append(c)
    pathNo.append(c)

    #### output of classification files #####
    tmpOut=subset_cerf.copy()
    tmpOut["behaviour"]=behaviourVector
    tmpOut["path_no"]=pathNo
    tmpOut.drop(columns=["dateTime"],axis=1).iloc[:, 1:].to_csv(outFile, index=False)
    finSTR = f"---> behaviour vector for stepSize {stepSize} successfully added. The output file is {outFile}"
    #### output of step classification file #####
    #stepOut=behaviourDF.copy()
    #stepOut["behaviour"]=colorBehaviour[:-1]

    #stepOut.to_csv(stepOutFile)
    return finSTR





def run(inFilePath, stepListFile, check_existing_file):    
    

    print("**** RUN S4_get_behaviour_vector_from_list.py ****")
    try:
        print(f"candidate step list : {stepListFile}")
        stepList=pd.read_csv(os.path.abspath(stepListFile))
        outDir="/".join(stepListFile.split("/")[:-1])
        print(f"input file : {inFilePath}")
        print(f" output directory : {outDir}")
        print(f"overwriting file ? {check_existing_file}")
    except FileNotFoundError as fnf:
        print("the list you provided does not exist")
    if stepList.empty:
        print("the stepListFile is empty : aborting script..")
        sys.exit(1)
    else:   
        with ProcessPoolExecutor(max_workers=32) as executor:
                
            results=list(executor.map(launcher, [inFilePath]*len(stepList),[outDir]*len(stepList), stepList.prenom, stepList.deerYear, stepList.stepSize, [check_existing_file]*len(stepList)))
        
            for res in results:
                print(res)

    print("******* Script S4 DONE *********")
    executor.shutdown()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=''' ** run the "getBehaviourVector" function from a list of divider length associated with an animalName and deerYear :
    traditionally, the list would be the output of the "get_Zscore_CRW.R" script ''')
    parser.add_argument("inFilePath", type=str, help=''' FULL path to the GPS file in .csv''')
    parser.add_argument("stepListFile", type=str,help=''' FULL path to the csv file (table) with the candidate stepSizes. The table should have at least the three following
        columms : "prenom", "deerYear" and "stepSize"
        NOTE : the script will try to find the [animalName]_[deerYear].csv within the [animalName] folder at the basis of the folder where the script is.. 
        so be careful where your GPS files are !!   ''')
    parser.add_argument("--check_existing_file", action="store_true",help='''
    whether you want to check if the behaviour files already exists and if so, NOT overwrite them''')
    args=parser.parse_args()

    #run(args.inFilePath, args.stepListFile, args.check_existing_file)

