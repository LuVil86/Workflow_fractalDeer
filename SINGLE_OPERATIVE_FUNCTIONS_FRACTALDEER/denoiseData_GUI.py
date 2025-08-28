
import geopandas as gpd
import pandas as pd
import denoiseFunctions
import sys
import os
import tkinter as tk
import pyperclip
import matplotlib.pyplot as plt
import matplotlib.pylab as p
from matplotlib.text	import Text
from matplotlib.widgets import Button, TextBox
from matplotlib.backend_bases import MouseButton
import matplotlib.gridspec as gridspec
from matplotlib import use
#use('QT5Agg')
from matplotlib.backends import qt_compat


class denoiseGUI(object):
    def __init__(self):
#### infile and nbPass

        self.workDir = None
        self.animalName= None
        self.deerYear= None

        self.gdf_cerf=None
        self.gdf_denoised=None
        self.nbPass=0


        ### denoise parameters
        self.sh_f1=-0.6
        self.multiplier1_f1=4
        self.multiplier2_f1=4
        self.speed_f1=2000 ### > 2km de déplacement à l'heure
        self.netSpeed_f1=True

        self.fig = plt.figure(figsize=(20,20))
        self.fig.subplots_adjust(bottom=0.04, right=0.95, top=0.95, wspace=0.1, hspace=0)
        gs = gridspec.GridSpec(7, 7, height_ratios=[1,1,1,1,1,1,10],width_ratios=[1,1,1,1,1,1,1])
        

        btnRunContainer=self.fig.add_subplot(gs[4,0:5])
        self.btnRun = Button(btnRunContainer, label="RUN", color="lightblue", hovercolor="blue")
        self.btnRunContainer = btnRunContainer


        ax1=self.fig.add_subplot(gs[6,0:3])
        ax1.ticklabel_format(style='plain')
        ax1.set_xlabel("Coordonnée X  [m]")
        ax1.set_ylabel("Coordonnée Y  [m]")
        self.ax1 = ax1

        ax2=self.fig.add_subplot(gs[6,4:])
        ax2.ticklabel_format(style='plain')
        ax2.set_xlabel("Coordonnée X [m]")
        ax2.set_ylabel("Coordonnée Y [m]")

        self.ax2 = ax2


        txt1=self.fig.add_subplot(gs[0,:])
        self.tb1=TextBox(txt1, label="Workdir",initial='/home/luvil/IE-OFEV/new_cerfs_valais_2025/')
        self.txt1=txt1

        txt2=self.fig.add_subplot(gs[1,:])
        self.tb2=TextBox(txt2, label="animal name", initial="ID225", color="white", hovercolor="white")
        self.txt2=txt2

        txt3=self.fig.add_subplot(gs[2,:])
        self.tb3=TextBox(txt3, label="deer year", initial="deerYear_2023-2024", color="white", hovercolor="white")
        self.txt3=txt3

        txt4=self.fig.add_subplot(gs[3,0])
        self.tb4=TextBox(txt4, label="pass n°",textalignment="center", initial=self.nbPass, color="white", hovercolor="lightgrey")
        self.txt4=txt4

        txt5=self.fig.add_subplot(gs[3,2])
        self.tb5=TextBox(txt5, label="sharpness",textalignment="center", initial=self.sh_f1, color="white", hovercolor="lightgrey")
        self.txt5=txt5

        txt6=self.fig.add_subplot(gs[3,4])
        self.tb6=TextBox(txt6, label="speed",textalignment="center",initial=self.speed_f1, color="white", hovercolor="lightgrey")
        self.txt6=txt6

        txt7=self.fig.add_subplot(gs[3,6])
        self.tb7=TextBox(txt7, label="ratio",textalignment="center",initial=self.multiplier2_f1, color="white", hovercolor="lightgrey")
        self.txt7=txt7

        results1=self.fig.add_subplot(gs[4,6])
        self.res=TextBox(results1, label="nb Points removed",textalignment="center", color="red", hovercolor="lightgrey")
        self.results1=results1

        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key_press)
    def on_key_press(self,event):
    # Check if the key pressed was Ctrl-V (paste)
        if event.key == "ctrl+v" and event.inaxes in [self.txt1]:
            self.tb1.set_val("")
            self.tb1.text_disp.set_color("black")

            self.tb1.set_val(pyperclip.paste().rstrip().strip())
        elif event.key == "ctrl+v" and event.inaxes in [self.txt2]:
            self.tb2.set_val("")
            self.tb2.set_val(pyperclip.paste().rstrip().strip())
            self.tb2.text_disp.set_color("black")

        elif event.key == "ctrl+v" and event.inaxes in [self.txt3]:
            self.tb3.set_val("")
            self.tb3.set_val(pyperclip.paste().rstrip().strip())
            self.tb3.text_disp.set_color("black")
        self.fig.canvas.draw()
        

    def on_click(self,event):
        if event.button == 1 and event.inaxes in [self.txt1]:
            self.tb1.text_disp.set_color("red")
            self.fig.canvas.draw()
        if event.button == 1 and event.inaxes in [self.txt2]:
            self.tb2.text_disp.set_color("red")
            self.fig.canvas.draw()
        if event.button == 1 and event.inaxes in [self.txt3]:
            self.tb3.text_disp.set_color("red")
            self.fig.canvas.draw()
        if event.button == 1 and event.inaxes in [self.btnRunContainer]:
            
            self.workDir = self.tb1.text
            self.animalName = self.tb2.text
            self.deerYear = self.tb3.text
            self.nbPass = int(self.tb4.text)
            self.sh_f1 = float(self.tb5.text)
            self.speed_f1 = int(self.tb6.text)
            self.multiplier2_f1 = float(self.tb7.text)
            if self.nbPass == 0:
                tmp=pd.read_csv(os.path.join(self.workDir,self.animalName, f"{self.animalName}_{self.deerYear}.csv" ))
            else:
                tmp=pd.read_csv(os.path.join(self.workDir,self.animalName, f"{self.animalName}_{self.deerYear}_denoised_{self.nbPass}.csv" ))

            self.gdf_cerf=gpd.GeoDataFrame(tmp, geometry=gpd.GeoSeries.from_wkt(tmp.geometry),crs=2056)


            self.gdf_denoised,self.nbPointRemoved=denoiseFunctions.denoisePipeline(self.gdf_cerf,
            netSpeed=self.netSpeed_f1,
            angle_sharpness=self.sh_f1,
            net2Points_Multiplier=self.multiplier2_f1,
            filterIndex=self.nbPass,
            maxSpeed=self.speed_f1)

            ## show nb points removed
            self.res.text_disp.set_color('white')
            self.res.text_disp.set_fontsize(18)
            self.res.set_val(self.nbPointRemoved)

            ### show trajectories
            self.ax1.clear()
            self.ax2.clear()
            
            if self.nbPass == 0:
                self.ax1.set_title(f"Données GPS Originales : {self.animalName}  : {self.gdf_cerf.crs}",fontsize=10)
            else:
                self.ax1.set_title(f"filtrage n° {self.nbPass} : {self.animalName}",fontsize=10)


            self.ax1.scatter(self.gdf_cerf.geometry.x,self.gdf_cerf.geometry.y, color="black")
            self.ax1.plot(self.gdf_cerf.geometry.x,self.gdf_cerf.geometry.y,linestyle="--", color="dimgrey")


            self.ax2.set_title(f"filtrage n° {self.nbPass+1} : {self.animalName}",fontsize=10)
            self.ax2.scatter(self.gdf_denoised.geometry.x,self.gdf_denoised.geometry.y, color="orange",label="coordonnées GPS")
            self.ax2.plot(self.gdf_denoised.geometry.x,self.gdf_denoised.geometry.y,linestyle="--", color="dimgrey",label="Trajet")

            self.ax2.plot([], [], ' ', label=f"angle minimum :{self.sh_f1}")
            self.ax2.plot([], [], ' ', label=f"maxDist : {self.speed_f1}")
            self.ax2.plot([], [], ' ', label=f"ratio : {self.multiplier2_f1}")
            self.ax2.legend(loc="best")
            self.fig.canvas.draw()
                
            self.gdf_denoised.to_csv(os.path.join(self.workDir,self.animalName, f"{self.animalName}_{self.deerYear}_denoised_{self.nbPass+1}.csv"), index=False)


### plt.show halt the script until the figure is closed
#plt.savefig(os.path.join(workDir,animalName, f"{animalName}_{deerYear}._1st_denoising.png"))
if __name__=="__main__":
   gui = denoiseGUI()
   p.show()