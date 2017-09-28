#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plot
import sys

data = np.genfromtxt("solution.dat",delimiter=' ',names=['t','V','m','h','n'])

# Grafico de V
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("V")
axis.set_xlabel('t')
axis.set_ylabel('V')
axis.plot(data['t'],data['V'],c='b',label='V')
leg = axis.legend()
plot.grid()
plot.show()

# Grafico de m
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("m")
axis.set_xlabel('t')
axis.set_ylabel('m')
axis.plot(data['t'],data['m'],c='g',label='m')
leg = axis.legend()
plot.grid()
plot.show()

# Grafico de h
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("h")
axis.set_xlabel('t')
axis.set_ylabel('h')
axis.plot(data['t'],data['h'],c='r',label='h')
leg = axis.legend()
plot.grid()
plot.show()

# Grafico de n
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_title("n")
axis.set_xlabel('t')
axis.set_ylabel('n')
axis.plot(data['t'],data['n'],c='k',label='n')
leg = axis.legend()
plot.grid()
plot.show()

# Grafico dos gates
fig = plot.figure()
axis = fig.add_subplot(111)
axis.set_xlabel('t')
axis.plot(data['t'],data['m'],c='g',label='m')
axis.plot(data['t'],data['h'],c='r',label='h')
axis.plot(data['t'],data['n'],c='k',label='n')
leg = axis.legend()
plot.grid()
plot.show()
