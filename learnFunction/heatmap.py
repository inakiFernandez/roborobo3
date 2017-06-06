# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:41:07 2017

@author: fernandi

Plot two 2D heatmap with two ordered data files
"""
import matplotlib
matplotlib.use('svg')
import matplotlib.pyplot as plt

from matplotlib import gridspec
import matplotlib.colors as clrs
import os, sys
import numpy as np
from numpy import genfromtxt
import scipy as sp
import scipy.stats as stats
import scipy.interpolate as intp
import argparse
import glob
from scipy.stats import gaussian_kde
from scipy.stats import chi2
from statsmodels.regression.quantile_regression import QuantReg
import statsmodels.formula.api as smf
import warnings
from mpl_toolkits.mplot3d import Axes3D
from pylab import MaxNLocator

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Normalize or not.')
    parser.add_argument('datapath',
                    help='folder with heatmap data and for heatmap images')
    parser.add_argument('nbGen', type=int, help='number of generations')
    parser.add_argument('--norm', action='store_true',
                    help='normalize distances in heatmap?')
    parser.add_argument('--rank', action='store_true',
                    help='rank pheno and geno distances?')
    parser.add_argument('numRuns', type=int, help='number of runs experiment')
    args = parser.parse_args()
    normalize = args.norm
    datapath = args.datapath
    doRank = args.rank
    nbGen = args.nbGen 
    numRuns = args.numRuns
    #For Chi square distribution (k)
    #nbSamples = 50
    
    doFast = True
    dirFiles = glob.glob(datapath + "/*istance*") # os.listdir(datapath)
    
    dirFiles.sort(key=lambda f: int(f.split('Distance')[1].split('.')[0]))
    speciesThreshold = 3.0 #To change if needed
    
    bgcolor="gainsboro"
    gridcolor="white"
    plt.ioff()
    if normalize:
        bigData = []
        maxBeh = []
        for fname in dirFiles:
            if fname.endswith(".csv") and ("phenotypic" in fname):
                maxBeh.append(np.max(np.concatenate(genfromtxt(fname, delimiter=',')).ravel()))
                #bigData.append(np.concatenate(genfromtxt(#datapath + "/" +
                #                                         fname, delimiter=',')).ravel())
        #bigData = np.concatenate(bigData).ravel()
        #maxBeh = np.max(bigData)
        bigData = []
        maxGen = []
        for fname in dirFiles:
            if fname.endswith(".csv") and ("genotypic" in fname):
                generationData = genfromtxt(fname, delimiter=',')
                maxGen.append(np.max(np.concatenate(generationData).ravel()))
                bigData.append(np.concatenate(generationData).ravel())

        bigData = np.concatenate(bigData).ravel()
        maxGenAll = np.max(bigData)
        #Maybe TODO draw index of best
        #speciesThreshold /= maxGen

    ranksWilcox = []
    kendalltau = []
    spearmanR = []
    genoAllDist = []
    phenoAllDist = []
    print("Start heatmap treatment")
    if not doFast:
        print("Begin if")
        for idx in range(int(nbGen)):
            fnameG = datapath + "/genotypicDistance" + str(idx + 1) + ".csv"
            fnamePh = datapath + "/phenotypicDistance" + str(idx + 1) + ".csv"
    
            dataG = genfromtxt(fnameG, delimiter=',')
            dataPh = genfromtxt(fnamePh, delimiter=',')
    
            dataGNoDiag = [[dataG[i][j] for j in range(len(dataG[i])) if j not in [i]] 
                           for i in range(len(dataG))]
    
            dataPhNoDiag = [[dataPh[i][j] for j in range(len(dataPh[i])) if j not in [i]]
                            for i in range(len(dataPh))]
    
            genoAllDist.extend(dataGNoDiag)
            phenoAllDist.extend(dataPhNoDiag)
            
            if(doRank):
                x = np.concatenate(dataGNoDiag).ravel()
                y = np.concatenate(dataPhNoDiag).ravel()
                
    
                ranksPh = stats.rankdata(dataPhNoDiag)            
                ranksG = stats.rankdata(dataGNoDiag)
                kendalltau.append(stats.kendalltau(ranksG,ranksPh))
                spearmanR.append(stats.spearmanr(x,y))
    
    
                #Cluster scatterplot geno-pheno in discretized bins         
                fig = plt.figure("Bins scatter GenoPheno") # ,facecolor=bgcolor)        
                
                fig.set_size_inches(8,6)
                gs = gridspec.GridSpec(2, 3, width_ratios=[1,3.5,0.1], height_ratios=[3.5,1])
                
                ax = plt.subplot(gs[0,1],facecolor=bgcolor)
                #xN = [d/len(x) for d in x]
                #yN = [d/len(y) for d in y]
                nbBins = 20
                vmaxLogScale = int(np.power(10,1+int(np.trunc(np.log10(len(x))))))            
                   
                hist2DGenoPheno = np.histogram2d( # plt.hist2D(
                                             x, y, bins=nbBins, #per dimension
                                             range=[[0.0, maxGenAll], [0.0, 2.0]],
                                             #cmap="viridis",
                                             #norm = clrs.LogNorm(vmin=0.1, vmax=vmaxLogScale),
                                             #cmin=-1
                                             )
                xedges = hist2DGenoPheno[1]
                yedges = hist2DGenoPheno[2]
                
                for i,r in enumerate(hist2DGenoPheno[0]):
                    for j,f in enumerate(r):
                        if f==0.0 :
                            hist2DGenoPheno[0][i][j]=1e-1
                norm = clrs.LogNorm(vmin=0.1,vmax=vmaxLogScale)   
                plt.pcolormesh(xedges,yedges,hist2DGenoPheno[0].T, #imshow
                           #interpolation='nearest', 
                           #origin='low', 
                           norm = norm,
                           #extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]]
                           )
                ax.set_xlim([0.0,maxGenAll])
                ax.set_ylim([0.0,2.0])
                #ax.set_adjustable('box-forced')
                #ax.autoscale(False)
                #ax.yaxis.set_tick_params(labelsize=5)
                plt.tick_params(axis='both',which='both',bottom='off',top='off',
                                labelbottom='off',left='off',right='off',labelleft='off')
                #norm = clrs.Normalize(vmin=0.0,vmax=1.0)            
                #sm = plt.cm.ScalarMappable(cmap="viridis", 
                #                           norm=clrs.LogNorm(vmin=0.1, vmax=vmaxLogScale) # len(x))                                       
                                               #plt.Normalize(vmin=0.0,vmax=len(x))
                #                          )
                #sm.set_array([])            
                #logTickSet = [np.power(10.0,exp) for exp in range(-1,2+int(np.trunc(np.log10(len(x)))),1)] 
                #print(logTickSet)    
                #print(len(x))    
                plt.title("It. " + str(idx + 1) +" bins for genotypic-phenotypic relationship")                         
                axClr = plt.subplot(gs[0,2],facecolor=bgcolor)
                #clrBar = plt.colorbar(extend='min')#sm,extend='min',ticks=logTickSet)                                  
                                      #range(0,len(x), np.trunc(len(x)/5))]) # , norm=norm) #hist2DGenoPheno)
                #clrBar.ax.set_yticklabels(["0"] + [r"$10^{"+str(exp)+"}$" for exp in range(0,2+int(np.trunc(np.log10(len(x)))),1)])
    
                logTickSet = [np.power(10.0,exp) for exp in range(-1,2+int(np.trunc(np.log10(len(x)))),1)] 
                cb = matplotlib.colorbar.ColorbarBase(axClr,cmap='viridis',norm=norm,ticks=logTickSet)
                cb.ax.set_yticklabels(["0"] + [r"$10^{"+str(exp)+"}$" for exp in range(0,2+int(np.trunc(np.log10(len(x)))),1)])
               
                
                axl = plt.subplot(gs[0,0], sharey=ax,facecolor=bgcolor)
                plt.grid(True,color=gridcolor)
                
                density = gaussian_kde(y)
                covarFactor = 0.20
                xs = np.linspace(0.0,2.0,1000)
                density.covariance_factor = lambda : covarFactor
                density._compute_covariance()
                maxX = np.max(density(xs))
                plt.plot(density(xs)/maxX,xs)
                axl.set_ylabel("Phenotypic distance")
                            
                plt.ylim(0.0,2.0)
                plt.xlim(0.0,1.05) # maxX)
                plt.gca().invert_xaxis()
    
                axb = plt.subplot(gs[1,1], sharex=ax,facecolor=bgcolor)
                plt.grid(True,color=gridcolor)
                
                density = gaussian_kde(x)
                ys = np.linspace(0.0,maxGenAll,1000)
                density.covariance_factor = lambda : covarFactor
                density._compute_covariance()
                maxY = np.max(density(ys))
                plt.plot(ys,density(ys)/maxY)
                axb.set_xlabel("Genotypic distance")
                plt.ylim(0.0,1.05) # maxY)
                plt.xlim(0.0,maxGenAll)
                #axb.set_adjustable('box-forced')
                #axb.autoscale(False)
                
                
                plt.savefig(datapath + "/binsScatterRelGenoPheno" + str(idx + 1) + ".png",dpi=100)
                plt.close(fig)
                
                figGenoPheno = plt.figure("Scatter GenoPheno",facecolor=bgcolor)
                #plt.grid(True,color=gridcolor)
                # definitions for the axes            
                gs = gridspec.GridSpec(2, 2, width_ratios=[1,3], height_ratios=[3,1])
                
                ax = plt.subplot(gs[0,1],facecolor=bgcolor)
                plt.grid(True,color=gridcolor)
                plt.scatter(x,y,c='blue',s=5,alpha=np.min([100/len(x),1]))
                ax.set_xlim([0.0,maxGenAll])
                ax.set_ylim([0.0,2.0])
                
                warnings.filterwarnings("ignore")
                #mod = QuantReg(y,x)
                mod = smf.quantreg('y~x',{'x': x, 'y': y})
                quantiles = np.array([0.05,0.25,0.5,0.75,0.95])# np.arange(.05, .96, .1)
                
                def fit_model(q):
                    res = mod.fit(q=q)                
                    return [q, res.params['Intercept'], res.params['x']] + res.conf_int().ix['x'].tolist()
    
                models = [fit_model(xvar) for xvar in quantiles]
                for i in range(len(quantiles)): # models.shape[0]):
                    yQ = models[i][1] + models[i][2] * x
                    quantcolor='grey'
                    quantwidth=0.2
                    if quantiles[i] == 0.5:
                        quantcolor='red'
                        quantwidth=1
                    ax.plot(x, yQ, linestyle=':', color=quantcolor,linewidth=quantwidth)
                
                warnings.filterwarnings('default')
    
                #modeTheilSlopes = stats.theilslopes(y, x, alpha=0.95)
                
                #plt.plot(x, modeTheilSlopes[1] + modeTheilSlopes[2] * x, '--',c='gray',linewidth=0.2)
                #plt.plot(x, modeTheilSlopes[1] + modeTheilSlopes[3] * x, '--',c='gray',linewidth=0.2)
                #plt.plot(x, modeTheilSlopes[1] + modeTheilSlopes[0] * x, 'r-',linewidth=0.3)
    
                plt.ylim(0.0,2.0)
                plt.xlim(0.0,maxGenAll)
                plt.title("It. " + str(idx + 1) + " genotype and phenotype distances")
                axl = plt.subplot(gs[0,0], sharey=ax,facecolor=bgcolor)
                plt.grid(True,color=gridcolor)
                
                #weights = np.ones_like(y)/len(y)
                #n =int(len(x)/100)
                #histy=axl.hist(y, bins=n, orientation='horizontal',weights=weights,
                #               histtype="stepfilled")
                density = gaussian_kde(y)
                #covarFactor = 0.33
                xs = np.linspace(0.0,2.0,10000)
                density.covariance_factor = lambda : covarFactor
                density._compute_covariance()
                maxX = np.max(density(xs))
                plt.plot(density(xs)/maxX,xs)
                axl.set_ylabel("Phenotypic distance")
    
                #chi2Dist = chi2(nbSamples,np.mean(y))
                #plt.plot(chi2Dist.pdf(xs)/np.max(chi2Dist.pdf(xs)),xs,c='red')
                
                plt.ylim(0.0,2.0)
                plt.xlim(0.0,1.05) # maxX)
                plt.gca().invert_xaxis()
                #pYD, yD = np.histogram(y, bins=n) # bin it into n = N/10 bins
                #yD = yD[:-1] + (yD[1] - yD[0])/2   # convert bin edges to centers
                #f = intp.UnivariateSpline(yD, pYD, s=n)
                #plt.plot(yD, f(yD))
                axb = plt.subplot(gs[1,1], sharex=ax,facecolor=bgcolor)
                plt.grid(True,color=gridcolor)
                #weights = np.ones_like(x)/len(x)
                #histx = axb.hist(x, bins=n,weights=weights,histtype="stepfilled")            
                density = gaussian_kde(x)
                ys = np.linspace(0.0,maxGenAll,10000)
                density.covariance_factor = lambda : covarFactor
                density._compute_covariance()
                maxY = np.max(density(ys))
                plt.plot(ys,density(ys)/maxY)
                plt.ylim(0.0,1.05) # maxY)
                plt.xlim(0.0,maxGenAll)
                axb.set_xlabel("Genotypic distance")
                 
                plt.savefig(datapath + "/scatterRelGenoPheno" + str(idx + 1) + ".png")
                plt.close(figGenoPheno)
    
            # ATTENTION: bricole pour améliorer la visualisation (coeffVis)
            coeffVis=1.0
            if normalize:
                dataG = [[np.min([(datum/maxGen[idx]) * coeffVis,1.0]) 
                          for datum in dataRow] for dataRow in dataG]
                dataPh = [[np.min([(datum/maxBeh[idx]) * coeffVis,1.0]) 
                           for datum in dataRow] for dataRow in dataPh]
    
            fig = plt.figure(fnameG)
            
            fig.set_size_inches(8, 6)        
            heatmap = plt.imshow(dataG, cmap=plt.cm.RdBu, interpolation='nearest') #, vmin=0.0, vmax=1.0)
            """        
            ax = fig.add_subplot(111,projection='3d')
            
            for i,row in enumerate(dataG):
                #for d in row:
                ys = [i] * len(row)
                xs = np.arange(len(row))
                zs = np.zeros_like(xs)
                ax.bar3d(xs, ys, zs, 1,1,row, 
                         color=plt.cm.RdBu(row/max(row)), alpha=0.3)
    
            ya = ax.get_yaxis()
            xa = ax.get_xaxis()
            ya.set_major_locator(MaxNLocator(integer=True))
            xa.set_major_locator(MaxNLocator(integer=True))
            """
            plt.title(fnameG)
            cbar = plt.colorbar(heatmap)
            #    plotThreshold = speciesThreshold
            #    if not normalize:
            #        plotThreshold /= np.max(data)
            #    cbar.ax.hlines(plotThreshold, 0, 1, colors = 'b', linewidth = 1)
            plt.savefig(os.path.splitext(fnameG)[0] + ".png",dpi=100)
            plt.close(fig)
    
            fig = plt.figure(fnamePh)
            
            fig.set_size_inches(8,6)
            heatmap = plt.imshow(dataPh, cmap=plt.cm.RdBu, interpolation='nearest') #, vmin=0.0, vmax=1.0)
            """
            ax = fig.add_subplot(111,projection='3d')
            
            for i,row in enumerate(dataPh):            
                ys = [i] * len(row)
                xs = np.arange(len(row))
                zs = np.zeros_like(xs)
                ax.bar3d(xs, ys, zs, 1,1,row, 
                         color=plt.cm.RdBu(row/max(row)), alpha=0.3)
            
            ya = ax.get_yaxis()
            xa = ax.get_xaxis()
            ya.set_major_locator(MaxNLocator(integer=True))
            xa.set_major_locator(MaxNLocator(integer=True))
            """
            plt.title(fnamePh)
            cbar = plt.colorbar(heatmap)
            plt.savefig(os.path.splitext(fnamePh)[0] + ".png",dpi=100)
            plt.close(fig)
            #TODO Done to speed up post analysis data generation
            #if doFast:
                #break               
    else:
        print("Begin else")
        for idx in range(int(nbGen)):  
            fnameG = datapath + "/genotypicDistance" + str(idx + 1) + ".csv"
            fnamePh = datapath + "/phenotypicDistance" + str(idx + 1) + ".csv"
            print("It. " + str(idx))
            #Cluster scatterplot geno-pheno in discretized bins         
            fig = plt.figure("Bins scatter GenoPheno") # ,facecolor=bgcolor)                            
            fig.set_size_inches(8,6)
            gs = gridspec.GridSpec(2, 3, width_ratios=[1,3.5,0.1], height_ratios=[3.5,1])                    
            ax = plt.subplot(gs[0,1],facecolor=bgcolor)
            ax.set_xlim([0.0,maxGenAll])
            ax.set_ylim([0.0,2.0])
            plt.title("It. " + str(idx + 1) +" bins for genotypic-phenotypic relationship")                         
            axClr = plt.subplot(gs[0,2],facecolor=bgcolor)
            axl = plt.subplot(gs[0,0], sharey=ax,facecolor=bgcolor)
            plt.grid(True,color=gridcolor)                    
            axl.set_ylabel("Phenotypic distance")                                
            plt.ylim(0.0,2.0)
            plt.xlim(0.0,1.05) # maxX)
            plt.gca().invert_xaxis()        
            axb = plt.subplot(gs[1,1], sharex=ax,facecolor=bgcolor)
            plt.grid(True,color=gridcolor)                    
            axb.set_xlabel("Genotypic distance")
            plt.ylim(0.0,1.05) 
            plt.xlim(0.0,maxGenAll)                    
            plt.savefig(datapath + "/binsScatterRelGenoPheno" + str(idx + 1) + ".png",dpi=100)
            plt.close(fig)

            
            figGenoPheno = plt.figure("Scatter GenoPheno",facecolor=bgcolor)
            gs = gridspec.GridSpec(2, 2, width_ratios=[1,3], height_ratios=[3,1])                    
            ax = plt.subplot(gs[0,1],facecolor=bgcolor)
            plt.grid(True,color=gridcolor)
            ax.set_xlim([0.0,maxGenAll])
            ax.set_ylim([0.0,2.0])                   
            #warnings.filterwarnings("ignore")

            plt.ylim(0.0,2.0)
            plt.xlim(0.0,maxGenAll)
            plt.title("It. " + str(idx + 1) + " genotype and phenotype distances")
            axl = plt.subplot(gs[0,0], sharey=ax,facecolor=bgcolor)
            plt.grid(True,color=gridcolor)                    
            axl.set_ylabel("Phenotypic distance")        
            plt.ylim(0.0,2.0)
            plt.xlim(0.0,1.05) 
            plt.gca().invert_xaxis()
            axb = plt.subplot(gs[1,1], sharex=ax,facecolor=bgcolor)
            plt.grid(True,color=gridcolor)
            plt.ylim(0.0,1.05) 
            plt.xlim(0.0,maxGenAll)
            axb.set_xlabel("Genotypic distance")                    
            plt.savefig(datapath + "/scatterRelGenoPheno" + str(idx + 1) + ".png")
            plt.close(figGenoPheno)

            fig = plt.figure(fnameG)                
            fig.set_size_inches(8, 6)   
            plt.plot(range(1))
            plt.title(fnameG)                
            plt.savefig(os.path.splitext(fnameG)[0] + ".png",dpi=100)
            plt.close(fig)
    
            fig = plt.figure(fnamePh)                
            fig.set_size_inches(8,6)
            plt.plot(range(1))
            plt.title(fnamePh)                
            plt.savefig(os.path.splitext(fnamePh)[0] + ".png",dpi=100)
            plt.close(fig)
                     

    if(doRank):
        if not doFast:
            minKendall = np.min([datum[0] for datum in kendalltau])
            minPValKendall = np.min([datum[1] for datum in kendalltau])
            minSpearman = np.min([datum[0] for datum in spearmanR])
            minPValSpearman = np.min([datum[1] for datum in spearmanR])
            minCorrRanks = np.min([minKendall,minSpearman])
            minPValRanks = np.min([minPValKendall,minPValSpearman])
    
            kendallFile = open(datapath + '/kendalltau' + str(numRuns) + '.data', 'w')
            for item in kendalltau:
                kendallFile.write("%s,%s\n" % item)
            kendallFile.close()
            
            figureKendall= plt.figure("Correlation p-values",facecolor=bgcolor)
            plt.plot([datum[1] for datum in kendalltau],label="Kendall-τ p-val.")
            plt.plot([datum[1] for datum in spearmanR],label="Spearman-ρ p-val.")
            plt.legend(loc='upper right', fontsize = 9)
            plt.ylim(ymin=minPValRanks)
            plt.grid(True,color=gridcolor)   
            figureKendall.get_axes()[0].set_facecolor(bgcolor)
            plt.savefig(datapath + "/correlationPval" + str(numRuns) + ".png")
            plt.close(figureKendall)
            
            figureKendall= plt.figure("Correlation Stats",facecolor=bgcolor)
            
            plt.plot([datum[0] for datum in kendalltau],label="Kendall-τ")
            plt.plot([datum[0] for datum in spearmanR],label="Spearman-ρ")
            plt.legend(loc='lower right', fontsize = 9)        
            plt.ylim(ymin=minCorrRanks,ymax=1.0)
            plt.grid(True,color=gridcolor)        
            figureKendall.get_axes()[0].set_facecolor(bgcolor)
            plt.savefig(datapath + "/correlationStats" +  str(numRuns) + ".png")
            plt.close(figureKendall)
    
            spearmanFile = open(datapath + '/spearman' +  str(numRuns) + '.data', 'w')
            for item in spearmanR:
                spearmanFile.write("%s,%s\n" % item)
            spearmanFile.close()
            #x = np.concatenate(genoAllDist).ravel()
            #y = np.concatenate(phenoAllDist).ravel()
            
            scatterFig =  plt.figure("scatterDist",facecolor=bgcolor)
            plt.grid(True,color=gridcolor)
            gs = gridspec.GridSpec(2, 1, width_ratios=[1], height_ratios=[1,28])
            ax = plt.subplot(gs[1,0],facecolor=bgcolor)
            plt.grid(True,color=gridcolor)        
            #print(len(genoAllDist))
            #print(genoAllDist[0])
            #plt.ion()
            
            for idx in range(nbGen):
                plt.scatter(genoAllDist[idx],phenoAllDist[idx],s=2,
                            alpha=1.0, #np.max([np.min([100.0/(nbGen),0.7]),0.1]), # *len(genoAllDist[idx])
                            c=plt.get_cmap('viridis')(idx/nbGen))
                #plt.show()
                #print(str(idx) + ',  ' + str(len(genoAllDist[idx])))
                
            ax.set_xlim([0.0,maxGenAll])
            ax.set_ylim([0.0,2.0])
            ax.set_xlabel("Genotypic distance")
            ax.set_ylabel("Phenotypic distance")
                
            viridisCMap = plt.get_cmap('viridis') # .Normalize(vmin=5, vmax=10)
            ax2 = plt.subplot(gs[0,0],facecolor=bgcolor)
            plt.grid(True,color=gridcolor)
            norm = clrs.Normalize(vmin=1,vmax=nbGen)
            cb = matplotlib.colorbar.ColorbarBase(ax2,cmap=viridisCMap,
                                             norm=norm,                                         
                                             orientation='horizontal') 
            cb.ax.xaxis.set_ticks_position('top') 
            plt.title("Genotypic vs. phenotypic distance. All generations.", y=-2.30)
            
            plt.ylim(0.0,2.0)
            plt.xlim(xmin=0.0)
            plt.savefig(datapath + "/scatterAllRelGenoPheno" +  str(numRuns) + ".png")
            plt.close(scatterFig)
        else:                       
            figureKendall= plt.figure("Correlation p-values",facecolor=bgcolor)            
            plt.grid(True,color=gridcolor)   
            figureKendall.get_axes()[0].set_facecolor(bgcolor)
            plt.savefig(datapath + "/correlationPval" + str(numRuns) + ".png")
            plt.close(figureKendall)
            
            figureKendall= plt.figure("Correlation Stats",facecolor=bgcolor)            
            
            plt.grid(True,color=gridcolor)        
            figureKendall.get_axes()[0].set_facecolor(bgcolor)
            plt.savefig(datapath + "/correlationStats" +  str(numRuns) + ".png")
            plt.close(figureKendall)
                
            scatterFig =  plt.figure("scatterDist",facecolor=bgcolor)
            plt.grid(True,color=gridcolor)
            gs = gridspec.GridSpec(2, 1, width_ratios=[1], height_ratios=[1,28])
            ax = plt.subplot(gs[1,0],facecolor=bgcolor)
            plt.grid(True,color=gridcolor)        
            
            ax.set_xlim([0.0,maxGenAll])
            ax.set_ylim([0.0,2.0])
            ax.set_xlabel("Genotypic distance")
            ax.set_ylabel("Phenotypic distance")
                
            viridisCMap = plt.get_cmap('viridis') # .Normalize(vmin=5, vmax=10)
            ax2 = plt.subplot(gs[0,0],facecolor=bgcolor)
            plt.grid(True,color=gridcolor)
            norm = clrs.Normalize(vmin=1,vmax=nbGen)
            cb = matplotlib.colorbar.ColorbarBase(ax2,cmap=viridisCMap,
                                             norm=norm,                                         
                                             orientation='horizontal') 
            cb.ax.xaxis.set_ticks_position('top') 
            plt.title("Genotypic vs. phenotypic distance. All generations.", y=-2.30)
            
            plt.ylim(0.0,2.0)
            plt.xlim(xmin=0.0)
            plt.savefig(datapath + "/scatterAllRelGenoPheno" +  str(numRuns) + ".png")
            plt.close(scatterFig)
    plt.show()

