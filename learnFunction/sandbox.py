# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 15:30:47 2017

@author: fernandi
"""
import resource

resource.setrlimit(resource.RLIMIT_DATA,(resource.getrlimit(resource.RLIMIT_DATA)[1],
                                                                    resource.getrlimit(resource.RLIMIT_DATA)[1]))
        resource.setrlimit(resource.RLIMIT_MEMLOCK,((resource.getrlimit(resource.RLIMIT_MEMLOCK)[1],
                                                                    resource.getrlimit(resource.RLIMIT_MEMLOCK)[1])))
        resource.setrlimit(resource.RLIMIT_AS,(resource.getrlimit(resource.RLIMIT_AS)[1],
                                                                    resource.getrlimit(resource.RLIMIT_AS)[1]))
        #resource.setrlimit(resource.RLIMIT_VMEM,(resource.getrlimit(resource.RLIMIT_VMEM)[1],
        #                                                            resource.getrlimit(resource.RLIMIT_VMEM)[1]))
        resource.setrlimit(resource.RLIMIT_STACK,(resource.getrlimit(resource.RLIMIT_STACK)[1],
                                                                    resource.getrlimit(resource.RLIMIT_STACK)[1]))





#import statsmodels.api as sm
#import statsmodels.formula.api as smf
#from statsmodels.regression.quantile_regression import QuantReg




        """
        quantiles = np.arange(.05, .96, .1)
        mod = smf.quantreg('genotypDist ~ phenotypDist',
                           [np.concatenate(genoAllDist).ravel(),np.concatenate(phenoAllDist).ravel()])
        res = mod.fit(q=.5)

        def fit_model(q):
            res = mod.fit(q=q)
            return [q, res.params['Intercept'], res.params['phenotypDist']] +  res.conf_int().ix['phenotypDist'].tolist()

        models = [fit_model(x) for x in quantiles]
        for i in range(models.shape[0]):
            y = get_y(models.a[i], models.b[i])
            ax.plot(np.concatenate(genoAllDist).ravel(), y, linestyle='dotted', color='grey')
        """
        
        
           #pXD, xD = np.histogram(x, bins=n) 
            #xD = xD[:-1] + (xD[1] - xD[0])/2   # convert bin edges to centers
            #f = intp.UnivariateSpline(xD, pXD, s=n)
            #plt.plot(xD, f(xD))            
            #axScatter = plt.axes(rect_scatter)
            #axHistx = plt.axes(rect_histx)
            #axHisty = plt.axes(rect_histy)
            #axHistx.xaxis.set_major_formatter(nullfmt)
            #axHisty.yaxis.set_major_formatter(nullfmt)
            #figGenoPheno.axes()
          

            binwidth = 0.25
            #xymax = np.max([np.max(np.fabs(x)), np.max(np.fabs(y))])
            #lim = (int(xymax/binwidth) + 1) * binwidth
            #axScatter.set_xlim((-lim, lim))
            #axScatter.set_ylim((-lim, lim))
            #bins = np.arange(-lim, lim + binwidth, binwidth)
            #axHistx.hist(x, bins=bins)
            #axHisty.hist(y, bins=bins, orientation='horizontal')
            #axHistx.set_xlim(axScatter.get_xlim())
            #axHisty.set_ylim(axScatter.get_ylim())
            


            
            #ranksWilcox.append(stats.wilcoxon(x,y))
            
        """        
        rankFile = open(datapath + '/wilcox.data', 'w')
        for item in ranksWilcox:
            rankFile.write("%s\n" % item[1])
        figureRanks = plt.figure("ranks")
        plt.plot([datum[1] for datum in ranksWilcox])
        plt.savefig(datapath + "/wilcox.png")
        plt.close(figureRanks)
        """
        
        
        
                #cb.set_ticks(range(nbGen))
        #cb.set_ticklabels(range(nbGen))
        #xAx2 = matplotlib.axis.XAxis(ax2) 
        #xAx2.set_ticks_position("top")                                 
        
        #plt.scatter(x,y, c='blue',s=2,alpha=0.03)
        #plt.plot(x, m*x + b, '-',c='green')
        #plt.plot(x, modeTheilSlopes[1] + modeTheilSlopes[0] * x, 'r-')
        #plt.plot(x, modeTheilSlopes[1] + modeTheilSlopes[2] * x, '--',c='gray')
        #plt.plot(x, modeTheilSlopes[1] + modeTheilSlopes[3] * x, '--',c='gray')
        
        
        
        
#mat = [[1, 2, 3, 1], [4, 5, 6, 4], [7, 7, 10, 10], [3, 7, 9, 4]]
#plt.matshow(mat,cmap = plt.cm.bwr) # ,origin='lower')

"""
import numpy.random

# Generate some test data
n = 16
x = np.random.randn(n)
y = np.random.randn(n)

data = [[x[i],y[i]] for i in range(n)]

heatmap, xedges, yedges = np.histogram2d(x, y, bins=4)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

plt.clf()
plt.imshow(heatmap.T, extent=extent, origin='lower')
plt.show()
"""

#plt.matshow(data,cmap = plt.cm.gray) #bwr) # ,origin='lower')





"""        
        figureSpearman= plt.figure("Spearman Rank Rho",facecolor=bgcolor)
        plt.grid(True,color=gridcolor)
        plt.plot([datum[0] for datum in spearmanR])
        plt.ylim(ymin=minCorrRanks)
        plt.savefig(datapath + "/spearmanRRho.png")
        plt.close(figureSpearman)
        """
        #scatterFile = open(datapath + '/scatter.data', 'w')

        #for item in [(g,p) for g in genoAllDist for p in phenoAllDist]:
        #    scatterFile.write("%s,%s\n" % item)
        #scatterFile.close()
        
        
                figureSpearman= plt.figure("Spearman Rank p-value",facecolor=bgcolor)
        plt.grid(True,color=gridcolor)
        plt.plot([datum[1] for datum in spearmanR])
        plt.ylim(ymin=minPValRanks)
        plt.savefig(datapath + "/spearmanRP-value.png")
        plt.close(figureSpearman)
        
        
        
        m, b = np.polyfit(x,y, 1)