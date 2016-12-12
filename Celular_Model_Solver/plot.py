#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plot
import sys

def plotHaq ():
    print "Haq"

def plotRudy ():
    print "Li & Rudy"

def plotFitzhugh ():
    data = np.genfromtxt("data.dat",delimiter=' ',names=['t','V','w'])
    # V
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("V")
    axis.set_xlabel('t')
    axis.set_ylabel('V')
    axis.plot(data['t'],data['V'],c='b',label='V')
    leg = axis.legend()
    plot.grid()
    plot.show()

    # w
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("w")
    axis.set_xlabel('t')
    axis.set_ylabel('w')
    axis.plot(data['t'],data['w'],c='g',label='w')
    leg = axis.legend()
    plot.grid()
    plot.show()

def plotMitchell ():
    data = np.genfromtxt("data.dat",delimiter=' ',names=['t','V','h'])
    # V
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("V")
    axis.set_xlabel('t')
    axis.set_ylabel('V')
    axis.plot(data['t'],data['V'],c='b',label='V')
    leg = axis.legend()
    plot.grid()
    plot.show()

    # h
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("h")
    axis.set_xlabel('t')
    axis.set_ylabel('h')
    axis.plot(data['t'],data['h'],c='g',label='h')
    leg = axis.legend()
    plot.grid()
    plot.show()

def plotNoble ():
    data = np.genfromtxt("data.dat",delimiter=' ',names=['t','V','m','h','n'])
    # V
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("V")
    axis.set_xlabel('t')
    axis.set_ylabel('V')
    axis.plot(data['t'],data['V'],c='b',label='V')
    leg = axis.legend()
    plot.grid()
    plot.show()

    # m
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("m")
    axis.set_xlabel('t')
    axis.set_ylabel('m')
    axis.plot(data['t'],data['m'],c='g',label='m')
    leg = axis.legend()
    plot.grid()
    plot.show()

    # h
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("h")
    axis.set_xlabel('t')
    axis.set_ylabel('h')
    axis.plot(data['t'],data['h'],c='r',label='h')
    leg = axis.legend()
    plot.grid()
    plot.show()

    # n
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("n")
    axis.set_xlabel('t')
    axis.set_ylabel('n')
    axis.plot(data['t'],data['n'],c='k',label='n')
    leg = axis.legend()
    plot.grid()
    plot.show()

def plotHodkin ():
    data = np.genfromtxt("data.dat",delimiter=' ',names=['t','V','m','n','h'])
    # V
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("V")
    axis.set_xlabel('t')
    axis.set_ylabel('V')
    axis.plot(data['t'],data['V'],c='b',label='V')
    leg = axis.legend()
    plot.grid()
    plot.show()

    # m
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("m")
    axis.set_xlabel('t')
    axis.set_ylabel('m')
    axis.plot(data['t'],data['m'],c='g',label='m')
    leg = axis.legend()
    plot.grid()
    plot.show()

    # n
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("n")
    axis.set_xlabel('t')
    axis.set_ylabel('n')
    axis.plot(data['t'],data['n'],c='r',label='n')
    leg = axis.legend()
    plot.grid()
    plot.show()

    # h
    fig = plot.figure()
    axis = fig.add_subplot(111)
    axis.set_title("h")
    axis.set_xlabel('t')
    axis.set_ylabel('h')
    axis.plot(data['t'],data['h'],c='k',label='h')
    leg = axis.legend()
    plot.grid()
    plot.show()


def main():

    if (len(sys.argv)-1 < 1 ):
        print "Usage:> python " + sys.argv[0] + " <id_model>"
        sys.exit(-1)

    id_model = int(sys.argv[1])

    if id_model == 1:
        plotMitchell()
    elif id_model == 2:
        plotNoble()
    elif id_model == 3:
        plotHodkin()
    elif id_model == 4:
        plotFitzhugh()
    elif id_model == 5:
        plotRudy()
    elif id_model == 6:
        plotHaq()

if __name__ == "__main__":
    main()
