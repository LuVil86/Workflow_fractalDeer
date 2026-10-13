import os
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk
import pandas as pd
import geopandas as gp
import locale
import astral
from astral.sun import sun
import pytz
from datetime import date

# --- Vos fonctions de mapping personnalisées --- #

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


def getSunPeriod(gdf:gp.GeoDataFrame, tz):
    assert all([isinstance(i, date) for i in gdf.dateTime])
    periodList=[]
    longLat=gdf.geometry
    local_tz = pytz.timezone(tz)

    for index, row in gdf.iterrows():
        loc_cerf=astral.LocationInfo("anyCountry",region="anyRegion",timezone=tz,latitude=longLat.iloc[index].y, longitude=longLat.iloc[index].x)

        s = sun(loc_cerf.observer, date=row.dateTime, tzinfo=local_tz)
        if index == 1:
            print(s)
        if local_tz.localize(row.dateTime)<=s["dawn"]:
            period="nuit"
        elif local_tz.localize(row.dateTime)<=s["sunrise"]:
            period="aube"
        elif local_tz.localize(row.dateTime)<=s["sunset"]:
            period="jour"
        elif local_tz.localize(row.dateTime)<=s["dusk"]:
            period="crepuscule"
        else:
            period="nuit"
        periodList.append(period)
    return periodList

def selectFile():
    fichier = filedialog.askopenfilename(
        title="Sélectionner le fichier CSV", 
        filetypes=[("Fichiers CSV", "*.csv"), ("Tous les fichiers", "*.*")]
    )
    if fichier:
        entry_file_path.set(fichier)

def run():
    # On récupère les variables et on les passe à open_child
    open_child(root, entry_animal, entry_file_path, entry_encoding)

def open_child(parent, var_animal, var_file_path, var_encoding):
    # Récupération des valeurs textuelles des variables passées en argument
    animal_test = var_animal.get().strip()
    file_path = var_file_path.get()
    file_encoding = var_encoding.get()

    # Vérification que les champs ne sont pas vides
    if not animal_test:
        messagebox.showwarning("Erreur", "Veuillez entrer un nom d'animal.")
        return
    if not file_path or not os.path.exists(file_path):
        messagebox.showwarning("Erreur", "Le fichier sélectionné est invalide.")
        return

    # Chargement du fichier selon l'encodage choisi
    try:
        if file_encoding == "windows":
            tmp = pd.read_csv(file_path, sep=";", encoding="iso8859_2", lineterminator="\r", decimal=",", skip_blank_lines=True)
        else:
            tmp = pd.read_csv(file_path, skip_blank_lines=True)
    except Exception as e:
        messagebox.showerror("Erreur", f"Impossible de lire le fichier : {e}")
        return

    availableColumns = list(tmp.columns.values)

    # Fenêtre enfant
    child = tk.Toplevel(parent)
    child.geometry("400x500")
    child.title("Configuration des colonnes")
    child.grab_set()
    child.transient(parent)

    def runChild():
        global dateTime, subset_cerf
        UTC_date = date_combo.get()
        UTC_time = time_combo.get()
        longitude = long_combo.get()
        latitude = lat_combo.get()
        tz = tz_combo.get()
        try:
            # Construction de la Serie DateTime
            dateTime = pd.to_datetime(tmp[UTC_date] + " " + tmp[UTC_time])
            print(dateTime)
            year = dateTime.dt.year
            mois = dateTime.dt.month_name(locale=locale.getlocale())
            heure = dateTime.dt.hour
            saison = dateTime.dt.month.apply(mapSeasons)
            nbYear = year.unique().tolist()
            
            
            subset_cerf=gp.GeoDataFrame({"prenom":animal_test,
                                            "annee":year,
                                            "mois":mois,
                                            "heure":heure,
                                            "dateTime":dateTime,
                                            "saison":saison, 
                                             "geometry":gp.points_from_xy(tmp[longitude], tmp[latitude]),
                                            "UTC_DATE":dateTime.dt.date,
                                            "UTC_TIME":dateTime.dt.time
                                            },
                                            crs="EPSG:4326").reset_index(drop=True).sort_values(by="dateTime", axis=0)
            subset_cerf["deerYear"]=getDeerYear(subset_cerf, yearList=nbYear)
            subset_cerf["sunPeriod"]=getSunPeriod(subset_cerf, tz)

            print(subset_cerf)
            child.grab_release()
            child.destroy()
            root.destroy()
        except Exception as e:
            messagebox.showerror("Erreur de traitement", f"Une erreur est survenue : {e}")

    # Interface de la fenêtre enfant
    tk.Label(child, text="Quelle colonne correspond à UTC_date ?", font=('Arial', 11)).pack(pady=10)
    date_combo = ttk.Combobox(child, values=availableColumns, state="readonly")
    date_combo.pack(pady=5)
    if availableColumns: date_combo.set(availableColumns[0])

    tk.Label(child, text="Quelle colonne correspond à UTC_time ?", font=('Arial', 11)).pack(pady=10)
    time_combo = ttk.Combobox(child, values=availableColumns, state="readonly")
    time_combo.pack(pady=5)
    if availableColumns: time_combo.set(availableColumns[0])
    tk.Label(child, text="Quelle time zone ?", font=('Arial', 11)).pack(pady=10)
    tz_combo = ttk.Combobox(child, values=pytz.all_timezones, state="readonly")
    tz_combo.pack(pady=5)
    tk.Label(child, text="Quelle colonne correspond à la longitude ?", font=('Arial', 11)).pack(pady=10)
    long_combo = ttk.Combobox(child, values=availableColumns, state="readonly")
    long_combo.pack(pady=5)
    if availableColumns: long_combo.set(availableColumns[0])
    tk.Label(child, text="Quelle colonne correspond à la latitude ?", font=('Arial', 11)).pack(pady=10)
    lat_combo = ttk.Combobox(child, values=availableColumns, state="readonly")
    lat_combo.pack(pady=5)
    if availableColumns: long_combo.set(availableColumns[0])
    ttk.Button(child, text="Valider et envoyer", command=runChild).pack(pady=20)






# --- Création de l'interface graphique Tkinter ---
root = tk.Tk()
root.title("Configuration du traitement des données")
root.geometry("500x450")

# Variables globales Tkinter
entry_file_path = tk.StringVar()
entry_encoding = tk.StringVar(value="windows") # Valeur par défaut "windows"

# Section 1 : Sélection du fichier
lbl_file = tk.Label(root, text="Fichier de données :")
lbl_file.pack(anchor="w", padx=10, pady=5)
frame_file = tk.Frame(root)
frame_file.pack(fill="x", padx=10)

entry_file = tk.Entry(frame_file, textvariable=entry_file_path, state="readonly", width=50)
entry_file.pack(side="left", fill="x", expand=True)
btn_browse = tk.Button(frame_file, text="Parcourir...", command=selectFile)
btn_browse.pack(side="right", padx=5)

# Encodage
lbl_enc = tk.Label(root, text="Encodage du fichier :")
lbl_enc.pack(anchor="w", padx=10, pady=5)
tk.Radiobutton(root, text="Windows (iso8859_2)", variable=entry_encoding, value="windows").pack(anchor="w", padx=20)
tk.Radiobutton(root, text="Unicode / Autre", variable=entry_encoding, value="unicode").pack(anchor="w", padx=20)

# Section 2 : Saisie de l'animal
lbl_animal = tk.Label(root, text="Nom de l'animal (animalTest) :")
lbl_animal.pack(anchor="w", padx=10, pady=5)

# Utilisation d'un StringVar aussi pour l'animal pour faciliter le transfert
entry_animal = tk.StringVar(value="Cerf_Test")
entry_animal_widget = tk.Entry(root, textvariable=entry_animal, width=50)
entry_animal_widget.pack(fill="x", padx=10, pady=5)

# Section 3 : Bouton de validation
btn_submit = tk.Button(root, text="Valider et Lancer le traitement", command=run, bg="#4CAF50", fg="white")
btn_submit.pack(pady=20)

# Lancement de la boucle principale
root.mainloop()