# -*- coding: utf-8 -*-
"""
Created on Fri May 12 15:22:53 2017

@author: fernandi
"""
import numpy as np
import matplotlib.pyplot as plt
import glob
import re

historyLineages = np.genfromtxt("./expRdmTopoInitNoEvoTopoStopOnError-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d10-lfalse-rtrue-ifalse-i50-i150-p0.0/lastPop.history")

historyLineages = np.genfromtxt("./expRdmTopoInitNoEvoTopoStopOnError2-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d100-lfalse-rtrue-ifalse-i50-i150-p0.0/lastPop.history")

for i in range(len(historyLineages)): plt.plot(historyLineages[i])
    
    
    
dirname = "./expRdmTopoInitNoEvoTopoStopOnError3-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d100-lfalse-rtrue-ifalse-i50-i150-p0.0/"
dirname = "./expRdmTopoInitNoEvoTopoStopOnError4-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d100-lfalse-rtrue-ifalse-i50-i150-p0.0/"

dirname = "./expRdmTopoInitNoEvoTopoStopOnError5-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d100-lfalse-rtrue-ifalse-i50-i150-p0.0/"

dirname = "./expRdmTopoInitNoEvoTopoStopOnError6-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d20-lfalse-rtrue-ifalse-i50-i150-p0.0/"

dirname = "./expRdmTopoInitNoEvoTopoStopOnError7-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d20-lfalse-rtrue-ifalse-i50-i150-p0.0/"

dirname = "./expRdmTopoInitNoEvoTopoStopOnError8-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d20-lfalse-rtrue-ifalse-i50-i150-p0.0/"

dirname = "./expRdmTopoInitNoEvoTopoStopOnError9-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d20-lfalse-rtrue-ifalse-i50-i150-p0.0/"

dirname = "./expRdmTopoInitNoEvoTopoStopOnError10-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d20-lfalse-rtrue-ifalse-i50-i150-p0.0/"
#read all pop.i files in current folder
filenames =  glob.glob(dirname + "pop*.history") #1-runFolder-

def numKey(s):
    regex = re.compile(r'\d+')
    return int(regex.findall(s)[-1]) 

filenames_sorted = sorted(filenames,key=numKey)

regex = re.compile(r'd\d+')
interv = int(regex.findall(dirname)[-1][1:])

for i,fname in enumerate(filenames_sorted):
    lineage = np.genfromtxt(fname)
    for j in range(len(lineage)):        
        plt.plot(lineage[j])
        plt.ylim([-0.05,1.05])
        plt.xlim([0,len(filenames)])    
    isT1 = True
    xcoords =  np.arange(0,len(filenames),interv)    
    for xc in xcoords:
        plt.axvline(x=xc,color='gray')
        if(isT1):
            plt.axvspan(xc, xc+interv, facecolor='grey', alpha=0.07)
        else:
            plt.axvspan(xc, xc+interv, facecolor='yellow', alpha=0.07)
        isT1 = not isT1
    print(str(i+1))
    plt.savefig(dirname + str(i+1)+'Lineage.png')
    plt.clf()
    

print("avconv -loglevel quiet -y -r 4 -start_number 1 -i " + dirname + "%dLineage.png -b:v 1000k -vcodec mpeg4 -t 100 " + dirname + "lineage.mp4")
print("avconv -loglevel quiet -y -i " + dirname + "lineage.mp4 " + dirname + "lineage.flv")
 
 
 
 
dirname = "./expRdmTopoInitNoEvoTopoStopOnError6-w0.1-c0.75-m0.6-m0.0-m0.0-m0.0-p100-d20-lfalse-rtrue-ifalse-i50-i150-p0.0/"

fname = dirname + "lastPop.history"
lineage = np.genfromtxt(fname)
for j in range(len(lineage)):        
    plt.plot(lineage[j])
    plt.ylim([-0.05,1.05])
    plt.xlim([0,len(lineage[0])])    
isT1 = True
xcoords =  np.arange(0,len(lineage[0]),interv)    
for xc in xcoords:
    plt.axvline(x=xc,color='gray')
    if(isT1):
        plt.axvspan(xc, xc+interv, facecolor='grey', alpha=0.07)
    else:
        plt.axvspan(xc, xc+interv, facecolor='yellow', alpha=0.07)
    isT1 = not isT1
    
 