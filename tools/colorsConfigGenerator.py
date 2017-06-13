# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 08:31:19 2017

@author: fernandi
"""
#
nItems = 200 # 30
nColors = 2 # 8
colorArr = []
if nColors != 2:
    for j in range(nColors):
        absVal = (256/nColors) * j;
        colorArr += [(absVal, 255 - absVal, 0)]
else:
    colorArr += [(0, 255, 0)]
    colorArr += [(255, 0, 0)]
#Divide by colors
outputStr = ""
for i in range(nItems):
    #outputStr += ('#%02x%02x%02x' % (colorArr[i % nColors])) + "\n"
    outputStr += ('physicalObject['+ str(i) + '].displayColorRed = %d' %  (colorArr[i % nColors][0])) + '\n'
    outputStr += ('physicalObject['+ str(i) + '].displayColorGreen = %d' %  (colorArr[i % nColors][1])) + '\n'
    outputStr += ('physicalObject['+ str(i) + '].displayColorBlue = %d' %  (colorArr[i % nColors][2])) + '\n'
    

text_file = open("/home/fernandi/git/roborobo3/logs/tmpColors-2colors.txt", "w")
#text_file.write()
#text_file.close()
print(outputStr, file=text_file)
