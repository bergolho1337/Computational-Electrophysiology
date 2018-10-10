import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def plot_potential (sv):
    plt.grid()
    plt.plot(sv[:,0],sv[:,1],label="Vm",c="black")
    plt.xlabel("t (ms)",fontsize=15)
    plt.ylabel("V (mV)",fontsize=15)
    plt.title("Action potential",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()    

def main():

    if (len(sys.argv) != 2):
        print("Usage:> %s <input_file>" % (sys.argv[0]))
        sys.exit(1)
    else:
        input_file = sys.argv[1]

        data = np.genfromtxt(input_file)

        plot_potential(data)        
    

if __name__ == "__main__":
    main()