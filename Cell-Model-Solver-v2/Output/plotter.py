import os
import sys
import numpy as np
from matplotlib import pyplot

def setNameStates (modelName):
    modelName = modelName[0:len(modelName)-1]   # Eliminate '\n'
    if (modelName == "Mitchell"):
        nameStates = ["V","m"]
    elif (modelName == "Noble"):
        nameStates = ["V","m","h","n"]
    elif (modelName == "DiFrancesco"):
        nameStates = ["V","Kc","Ki","Nai","y","x","Cai","s","m","h","d","f","f2","Ca_up","Ca_rel","p"]
    elif (modelName == "LuoRudy"):
        nameStates = ["V","m","h","j","Cai","d","f","X"]
    elif (modelName == "FitzHugh"):
        nameStates = ["V","h"]
    else:
        print("[-] Invalid model name!")
        sys.exit(-1)
    return nameStates

def plotSolution (nameStates):
    data = np.genfromtxt(open("Output/solution.dat","r"))
    for i in range(len(nameStates)):
        pyplot.clf()
        pyplot.title(nameStates[i] + " x Time")
        pyplot.xlabel("Time")
        pyplot.ylabel(nameStates[i])
        pyplot.plot(data[:,0],data[:,i+1],label=nameStates[i],linewidth=2,color="black")
        pyplot.grid()
        pyplot.legend(loc=0,fontsize=15)
        pyplot.savefig("Output/" + nameStates[i] + ".pdf")
        #pyplot.show()

def main ():
    with open("Output/logfile.log","r") as f:
        modelName = f.readline()
    nameStates = setNameStates(modelName)
    plotSolution(nameStates)

if __name__ == "__main__":
    main()