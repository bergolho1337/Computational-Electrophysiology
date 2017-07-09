#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plot
import sys

for i in range(10):
    filename = "cell"+str(i)+".dat"
    print("Ploting file %s" % filename)
    data = np.genfromtxt(filename,delimiter='\t',names=['t','V','ca'])

    # Grafico de V
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("V")
    axis.set_xlabel('t')
    axis.set_ylabel('V')
    axis.plot(data['t'],data['V'],c='b',label='V')
    leg = axis.legend()
    plot.grid()
    plot.ylim([-100,60])
    plot.savefig(filename+"_v.png")
    #plot.show()

    # Grafico de m
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("Ca")
    axis.set_xlabel('t')
    axis.set_ylabel('Ca')
    axis.plot(data['t'],data['ca'],c='g',label='ca')
    leg = axis.legend()
    plot.grid()
    plot.ylim([0,0.0003])
    plot.savefig(filename+"_ca.png")
    #plot.show()
