import sys
import os

print("\n"+"######### **** compute Fractal on observed and simulated trajectory ********* ########")

#### *** input parameters *** ######

workDir=os.path.dirname(sys.argv[2])
fileName=os.path.splitext(os.path.basename(sys.argv[1]))[0]
print(fileName)
animalTest=fileName.split("_")[0]

yearTest="_".join(fileName.split("_")[1:])
print(f" deer : {animalTest} -- deerYear : {yearTest}")

print("\n"+"######### **** compute Fractal on observed and simulated trajectory ********* ########")
print(f" deer : {animalTest} -- deerYear : {yearTest}")