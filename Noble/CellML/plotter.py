import os
import numpy as np
from matplotlib import pyplot

data = np.genfromtxt(open("solution.dat","r"))
name_states = ["V","m"]
for i in range(len(name_states)):
    pyplot.clf()
    pyplot.title(name_states[i] + " x Time")
    pyplot.xlabel("Time")
    pyplot.ylabel(name_states[i])
    pyplot.plot(data[:,0],data[:,i+1],label=name_states[i],linewidth=2,color="black")
    pyplot.grid()
    pyplot.legend(loc=0,fontsize=15)
    pyplot.savefig("Output/" + name_states[i] + ".pdf")

