import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def calc_derivatives (t,vm):
    n = len(vm)-1
    dvdt = {}
    for i in range(n):
        dvdt[t[i]] = vm[i+1]-vm[i] 
    return dvdt

def calc_apd (sv):
    t = sv[:,0]
    vm = sv[:,1]
    
    # Calculate APD resting value
    min_value = np.amin(vm)
    max_value = np.amax(vm)
    APD_90_value = abs(max_value - min_value)*0.1 + min_value

    dvdt = calc_derivatives(t,vm)
    max_dvdt = sorted(dvdt.iteritems(), key=lambda (k,v): (v,k), reverse=True)
    for i in range(3):
        print max_dvdt[i]
    #for key, value in sorted(dvdt.iteritems(), key=lambda (k,v): (v,k), reverse=True):
    #    print "%s: %s" % (key, value)
    




def main():

    if (len(sys.argv) != 2):
        print("Usage:> %s <input_file>" % (sys.argv[0]))
        sys.exit(1)
    else:
        input_file = sys.argv[1]

        data = np.genfromtxt(input_file)

        calc_apd(data)        
    

if __name__ == "__main__":
    main()