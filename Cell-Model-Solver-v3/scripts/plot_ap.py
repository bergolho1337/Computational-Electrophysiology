# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_solution (filename):
    data = np.genfromtxt(filename, delimiter=',')
    return data

def show_solution (data):
    plt.plot(data[:,0],data[:,1],label="V",c="black")

def write_solution ():
    plt.grid()
    plt.xlabel(u"t (ms)",fontsize=15)
    plt.ylabel(u"V (mV)",fontsize=15)
    plt.title(u"Action Potential",fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()
    #plt.savefig("output.pdf")

def main():

    solution_filename = sys.argv[1]

    sv = read_solution(solution_filename)
    
    show_solution(sv)
    write_solution()


if __name__ == "__main__":
    main()