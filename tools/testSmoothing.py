# -*- coding: utf-8 -*-
"""
Created on Wed May 10 12:00:30 2017

@author: fernandi
"""
import numpy as np
import matplotlib.pyplot as plt

def movingAvg(data, window = 10):
    result= np.cumsum(data, dtype=float)
    result[ window :] = result[window:] - result[:-window]
    result = result[ windowSize - 1:] / window
    ratioAbsc = float(len(data))/float(len(result))
    xMovAvg = [(i + 1) * ratioAbsc for i in range(len(result))]
    return result, xMovAvg

length = 10
x = [(np.exp(length)-np.exp(length-i)) + np.random.normal(0.0,np.exp(length)/10) 
     for i in np.arange(0.1,length,0.1)]

windowSize = 3

#movAvg = np.cumsum(x, dtype=float)
#movAvg[ windowSize :] = movAvg[ windowSize :] - movAvg[:-windowSize]
#movAvg = movAvg[ windowSize - 1:] / windowSize
movAvg,xMovAvg = movingAvg(x,windowSize)

deriv = [ np.absolute(movAvg[i+10]-movAvg[i]) / np.max(movAvg) for i in range(len(movAvg) - 10)]


#plt.plot(deriv)
#xMovAvg = [i + windowSize/2 for i in range(len(movAvg))]
#xRatio = len(x) / len(movAvg)
#xMovAvg = [(i + 1) * xRatio for i in range(len(movAvg))]
#plt.plot(xMovAvg,movAvg)
#plt.plot(x)


movAvg2,xMovAvg2 = movingAvg(movAvg,windowSize)
#xRatio2 = len(x) / len(movAvg2)
#xMovAvg2 = [(i + 1) * xRatio2 for i in range(len(movAvg2))]

movAvg3,xMovAvg3 = movingAvg(movAvg2,windowSize)
#xRatio3 = len(x) / len(movAvg3)
#xMovAvg3 = [(i + 1) * xRatio3 for i in range(len(movAvg3))]

movAvg4,xMovAvg4 = movingAvg(movAvg3,windowSize)
movAvg5,xMovAvg5 = movingAvg(movAvg4,windowSize)
movAvg6,xMovAvg6 = movingAvg(movAvg5,windowSize)
movAvg7,xMovAvg7 = movingAvg(movAvg6,windowSize)

plt.plot(xMovAvg7,movAvg7)
#plt.plot(xMovAvg6,movAvg6)
#plt.plot(xMovAvg5,movAvg5)
#plt.plot(xMovAvg4,movAvg4)
#plt.plot(xMovAvg3,movAvg3)
#plt.plot(xMovAvg2,movAvg2)
#plt.plot(xMovAvg,movAvg)
plt.plot(x)
