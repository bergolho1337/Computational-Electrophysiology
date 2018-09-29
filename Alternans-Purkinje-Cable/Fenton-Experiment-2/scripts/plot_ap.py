import os
import sys
import numpy as np
from matplotlib import pyplot

def plot_ap (input_filename, output_filename):
    data = np.genfromtxt(open(input_filename,"r"))

    pyplot.clf()
    pyplot.title("Action Potential")
    pyplot.xlabel("Time (ms)")
    pyplot.ylabel("V (mV)")
    pyplot.plot(data,label="v",linewidth=2,color="black")
    pyplot.grid()
    pyplot.legend(loc=0,fontsize=15)
    pyplot.show()
    #pyplot.savefig(output_filename)
    print("[+] Output file %s save with sucess !" % (output_filename))

def main():
    if (len(sys.argv) != 3):
        print("------------------------------------------------------")
        print("Usage:> " + sys.argv[0] + " <input_data_file> <output_file>")
        print("------------------------------------------------------")
        sys.exit(1)
    
    plot_ap(sys.argv[1],sys.argv[2])

if __name__ == "__main__":
    main()
