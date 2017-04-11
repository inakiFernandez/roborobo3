# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 18:23:06 2017

@author: fernandi
"""


fig1 = plt.figure(1)
axis1 = fig1.add_subplot(111)
for i,c in enumerate(colorArr):
    patchC = matplotlib.patches.Rectangle((i/8.0,0.1),0.125,0.6,facecolor=c)
    axis1.add_patch(patchC)
#axis1.text(0.05,0.8, "[(0.0, 255.0, 0.0), (32.0, 223.0, 0.0), (64.0, 191.0, 0.0), (96.0, 159.0, 0.0), (128.0, 127.0, 0.0), (160.0, 95.0, 0.0), (192.0, 63.0, 0.0), (224.0, 31.0, 0.0)]")



axis1.text(0.05,0.8, 
           "[(0.0, 0.0, 1.0), (0.0, 0.5686274509803921, 1.0), (0.0, 1.0, 0.8549019607843137), (0.0, 1.0, 0.2823529411764706), (0.2823529411764706, 1.0, 0.0), (0.8549019607843137, 1.0, 0.0), (1.0, 0.5686274509803921, 0.0), (1.0, 0.0, 0.0)]")