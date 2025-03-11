
import pandas as pd
import os
from io import StringIO
import geopandas as gp
import pytz
import astral
from astral.sun import sun
import locale
from datetime import datetime, date,time
import sys
'''
Author : Lucas Villard 
Last revised : Jan 28 2025

Usage : Scripts to format GPS data from collared red deers for the subsequent workflow.

This script works for traditional outputs of GPSplusX.

--- Inputs ---

the script needs the GPSplusX output file , the name you want to attribute to the individual as well as its sex.

-- Outputs ---


The script outputs one file with every valid fixes in the form of "[animalName] _allFixes_formatted.csv" and
several sub-files for every deerYear found in the data in the form of "[animalName]_[deerYear].csv"

all these files will be created in a folder named after the animalName, which location will be on the path where this script is

'''


############# *** INPUT PARAMETERS *** ##################################

animalName="Zima"
sexe="femelle"
inputFilePath = "/home/luvil/IE-OFEV/cerf_bern_19feb25/ZIMA_Collar37172_20250219120214.csv"
overwrite=True
workDir="/home/luvil/IE-OFEV/cerf_bern_19feb25"

#########################################################################


def mapSeasons(monthValue):
    corr={1:"Decembre-Fevrier",
          2:"Decembre-Fevrier",
          3:"Mars-Mai",
          4:"Mars-Mai",
          5:"Mars-Mai",
          6:"Juin-Aout",
          7:"Juin-Aout",
          8:"Juin-Aout",
          9:"Septembre-Novembre",
          10:"Septembre-Novembre",
          11:"Septembre-Novembre",
          12:"Decembre-Fevrier"}
    return corr[monthValue]   

def getDeerYear(gdf_cerf,yearList:list):
  deerYearVector=pd.Series(["NA" for i in range(len(gdf_cerf))])
  yearList.append(yearList[-1]+1)
  for i in yearList:
    if pd.to_datetime(i, format="%Y").is_leap_year:
        deerYear=[pd.Timestamp(i-1,3,1),pd.Timestamp(i,2,29)+pd.Timedelta(days=1)]
    else:
        deerYear=[pd.Timestamp(i-1,3,1),pd.Timestamp(i,2,28)+pd.Timedelta(days=1)]

    #deerYear=[pd.Timestamp(i-1,6,1),pd.Timestamp(i,5,31)+pd.Timedelta(days=1)]    
    
    for index, row in gdf_cerf.iterrows():
         if deerYear[0]<row.dateTime<=deerYear[1]:
              deerYearVector.iloc[index]=f"deerYear_{i-1}-{i}"
  return deerYearVector

def getSunPeriod(gdf_cerf:gp.GeoDataFrame, country:str):
    periodList=[]
    longLat=gdf_cerf.geometry.to_crs("epsg:4326")
    for index, row in gdf_cerf.iterrows():

        if country=="Switzerland":
            local_tz = pytz.timezone('Europe/Zurich')
            loc_cerf=astral.LocationInfo("Suisse",region="Switzerland",timezone="Europe/Zurich",latitude=longLat.iloc[index].y, longitude=longLat.iloc[index].x)
        elif country=="France":
            local_tz = pytz.timezone('Europe/Paris')
            loc_cerf=astral.LocationInfo("Haute-Savoie",region="France",timezone="Europe/Paris",latitude=longLat.iloc[index].y, longitude=longLat.iloc[index].x)
        else:
            print("no country was specified")
            raise Exception


        s = sun(loc_cerf.observer, date=row.dateTime)
        if local_tz.localize(row.dateTime)<=s["dawn"]:
            period="nuit"
        elif local_tz.localize(row.dateTime)<=s["sunrise"]:
            period="nuit"
        elif local_tz.localize(row.dateTime)<=s["sunset"]:
            period="jour"
        elif local_tz.localize(row.dateTime)<=s["dusk"]:
            period="jour"
        else:
            period="nuit"
        periodList.append(period)
    return periodList



if not sexe in ["male", "femelle", "unknown"]:
    print( " 'sex' variable should be either 'male', 'femelle' or 'unknown'")
    raise ValueError

#tmp=pd.read_csv(inputFilePath,sep=",", encoding='utf-8', decimal=".", skip_blank_lines=True)
tmp=pd.read_csv(inputFilePath,sep=None, engine="python", encoding="iso8859_2")
#tmp=pd.read_csv(inputFilePath, sep=";", encoding="iso8859_2", lineterminator="\r",decimal=",",skip_blank_lines=True)
neededColumns=['UTC_Date', 'UTC_Time', 'ECEF_X [m]', 'ECEF_Y [m]', 'Longitude [°]','Latitude [°]']
if not all([colNames in tmp.columns for colNames in neededColumns ]):
    print(' There are missing columns or wrong column names in your input data file : required columns : "UTC_Date", "UTC_Time", "ECEF_X [m]", "ECEF_Y [m]", "Longitude [°]","Latitude [°]" ')
    raise KeyError

tmp.dropna(axis=0, subset="UTC_Date", inplace=True)
tmp.dropna(axis=0, subset=["ECEF_X [m]","ECEF_Y [m]"], inplace=True)
tmp["prenom"]=animalName
try:
    tmp["dateTime"]=pd.to_datetime(tmp.UTC_Date+" "+tmp.UTC_Time, format="%d/%m/%Y %H:%M:%S")
except ValueError:
    tmp["dateTime"]=pd.to_datetime(tmp.UTC_Date+" "+tmp.UTC_Time, format="%d.%m.%Y %H:%M:%S")

tmp["year"]=pd.Series(tmp.dateTime).dt.year
tmp["mois"]=pd.Series(tmp.dateTime).dt.month_name(locale=locale.getlocale())
tmp["heure"]=pd.Series(tmp.dateTime).dt.hour
tmp["saison"]=tmp["dateTime"].dt.month.apply(mapSeasons)
tmp["sexe"]=sexe
nbYear=tmp["year"].unique().tolist()

subset_cerf=gp.GeoDataFrame({"prenom":tmp.prenom,
                             "sexe":tmp.sexe,
                             "annee":tmp.year,
                             "mois":tmp.mois,
                            "saison":tmp.saison, 
                            "heure":tmp.heure,
                            "dateTime":tmp.dateTime,
                            "geometry":gp.points_from_xy(tmp["Longitude [°]"], tmp["Latitude [°]"]),
                            "UTC_DATE":tmp.dateTime.dt.date,
                            "UTC_TIME":tmp.dateTime.dt.time
                            },
                            crs="EPSG:4326").reset_index(drop=True)
subset_cerf.to_crs(epsg=2056, inplace=True)
del(tmp)
print(f"********** {animalName} **********")
print(f"start : {subset_cerf.iloc[0].dateTime}, end : {subset_cerf.iloc[-1].dateTime} ")
subset_cerf["deerYear"]=getDeerYear(subset_cerf, yearList=nbYear)
subset_cerf["jourNuit"]=getSunPeriod(subset_cerf, country="Switzerland")
print(subset_cerf["deerYear"].value_counts())

gb = subset_cerf.groupby(['prenom', 'deerYear'])
agregat=gb.agg(
nbFixes=('dateTime', "count"),
startDT=('dateTime', "min"),
endDT=('dateTime', "max")
).reset_index()


## output writing ##
resultPath=os.path.join(workDir, animalName)

for yearTest in subset_cerf["deerYear"].unique():
    print(yearTest)
    subYear=subset_cerf.loc[subset_cerf["deerYear"]==yearTest].reset_index(drop=True)
    subYear.sort_values(by="dateTime", axis=0,inplace=True)
    try:
        os.mkdir(resultPath)
    except OSError as error:
       pass
    csvPath=os.path.join(resultPath,f"{animalName}_{yearTest}.csv")
    if not os.path.isfile(csvPath) or overwrite:
        subYear.to_csv(os.path.join(resultPath,f"{animalName}_{yearTest}.csv"), index=False)   
    else:
        print("file already exists : no overwriting")
try:
    os.mkdir(resultPath)
except OSError as error:
    pass
csvPath=os.path.join(resultPath,f"{animalName}_allFixes_formatted.csv")
if not os.path.isfile(csvPath) or overwrite:
    subset_cerf.to_csv(os.path.join(resultPath,f"{animalName}_allFixes_formatted.csv"), index=False)   
    subset_cerf.drop(["UTC_DATE", "UTC_TIME"],axis=1).to_file(os.path.join(resultPath,f"{animalName}_allFixes_formatted.gpkg"), driver='GPKG', layer=animalName,overwrite=overwrite)
else:
    print("file already exists : no overwriting")


agregat["Duration_days"]=(agregat["endDT"]-agregat["startDT"]).dt.days
agregat["is_full_year"]=agregat["Duration_days"].apply(lambda x : "yes" if x>=364 else "no")
agregat.to_csv(os.path.join(resultPath,f"{animalName}_fixes_counts.csv"), index=False)
