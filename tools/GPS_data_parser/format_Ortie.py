
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
import numpy as np



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

### environment parameters ###
workDir="/".join(sys.argv[0].split("/")[:-1]) ### WRITE RESULTS WHERE THE SCRIPT IS LOCATED

overwrite=True

###############################



allDat=pd.read_csv("/home/luvil/IE-OFEV/cerf_jura/ortie.csv", delimiter=",", decimal=".")
animalName = "Ortie"

tmp=allDat.loc[allDat["animals_id"]==animalName].copy()
# tmp.dropna(axis=0, subset="acquisition_time", inplace=True)
# tmp.dropna(axis=0, subset=["x_lv95","y_lv95"], inplace=True)
tmp["prenom"]=animalName
tmp["dateTime"]=pd.to_datetime(tmp.acquisition_date+" "+tmp.acquisition_time, format="%Y/%m/%d %H:%M:%S")
tmp["year"]=pd.Series(tmp.dateTime).dt.year
tmp["mois"]=pd.Series(tmp.dateTime).dt.month_name(locale=locale.getlocale())
tmp["heure"]=pd.Series(tmp.dateTime).dt.hour
tmp["saison"]=tmp["dateTime"].dt.month.apply(mapSeasons)
#tmp["sexe"]=cerfInfo.loc[cerfInfo["animals_original_id"]==animalName].sex
nbYear=tmp["year"].unique().tolist()
tmp = tmp.loc[tmp.duplicated(subset="dateTime", keep="first")==False].copy()

subset_cerf=gp.GeoDataFrame({"prenom":tmp.prenom,
                            "sexe":"male",
                            "annee":tmp.year,
                            "mois":tmp.mois,
                            "saison":tmp.saison, 
                            "heure":tmp.heure,
                            "dateTime":tmp.dateTime,
                            "geometry":gp.points_from_xy(tmp["longitude"].replace(",", ".", regex=True),tmp["latitude"].replace(",", ".", regex=True)),
                            "UTC_DATE":tmp.dateTime.dt.date,
                            "UTC_TIME":tmp.dateTime.dt.time
                            },
                            crs="EPSG:4326").reset_index(drop=True).sort_values(by="dateTime", axis=0)
del(tmp)
subset_cerf.to_crs(epsg=2056, inplace=True)
print(f"start : {subset_cerf.iloc[0].dateTime}, end : {subset_cerf.iloc[-1].dateTime} ")
subset_cerf["deerYear"]=getDeerYear(subset_cerf, yearList=nbYear)
subset_cerf["jourNuit"]=getSunPeriod(subset_cerf, country="Switzerland")
#print(subset_cerf["deerYear"].value_counts())


gb = subset_cerf.groupby(['prenom', 'deerYear'])
agregat=gb.agg(
nbFixes=('dateTime', "count"),
startDT=('dateTime', "min"),
endDT=('dateTime', "max")
).reset_index()
finalStat=agregat
print(finalStat)
resultPath=os.path.join(workDir, str(animalName))

for yearTest in subset_cerf["deerYear"].unique():
    subYear=subset_cerf.loc[subset_cerf["deerYear"]==yearTest].reset_index(drop=True)
    subYear.sort_values(by="dateTime", axis=0,inplace=True)
    try:
        os.mkdir(resultPath)
    except OSError as error:
        pass
    csvPath=os.path.join(resultPath,f"{str(animalName)}_{yearTest}.csv")
    if not os.path.isfile(csvPath) or overwrite:
        
        subYear.to_csv(os.path.join(resultPath,f"{str(animalName)}_{yearTest}.csv"), index=False)   
    else:
        print("file already exists : no overwriting")
try:
    os.mkdir(resultPath)
except OSError as error:
    pass
csvPath=os.path.join(resultPath,f"{animalName}_allFixes_formatted.csv")
if not os.path.isfile(csvPath) or overwrite:
    subset_cerf.to_csv(os.path.join(resultPath,f"{animalName}_allFixes_formatted.csv"), index=False)   
    subset_cerf.drop(["UTC_DATE", "UTC_TIME"],axis=1).to_file(os.path.join(workDir, "Ortie_allFixes_formatted.gpkg"), driver='GPKG', layer=str(animalName))
    
else:
    print("file already exists : no overwriting")

finalStat["Duration_days"]=(finalStat["endDT"]-finalStat["startDT"]).dt.days
finalStat["is_full_year"]=finalStat["Duration_days"].apply(lambda x : "yes" if x>=364 else "no")
finalStat.to_csv(os.path.join(workDir, "Ortie_fixes_counts.csv"), index=False)
