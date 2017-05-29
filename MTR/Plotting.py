import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas import Series, DataFrame
plt.style.use('ggplot')
from copy import copy

def data_mc(df_data,df_mc,var,bins,weight=None,mclabels=[""],normalize=True):
    
    weights = df_mc[weight].values if weight else None
    
    hists_mc = map(lambda x: np.histogram(df_mc['%s%s' % (var,x)].values,density=normalize,weights=weights,bins=bins)[0], 
                   mclabels)
    hist_data,_ =  np.histogram(df_data[var].values, density=False, bins=bins)

    if normalize:
        hists_mc = map(lambda x: x*hist_data.sum() / x.sum(), hists_mc)    
    
    return hist_data,hists_mc


def draw_data_mc(bins,data,mc,figsize=(8,6),var=None,logy=False,ratio=False,
                **kwargs):
    
    mcstyle=dict(alpha=0.5,linewidth=0)
    mcstyle.update(kwargs)
    datastyle=copy(kwargs)
    datastyle["alpha"] = 1. # No transparency for data
    
    binw=bins[1]-bins[0]
    if ratio:
        fig, axes = plt.subplots(2,figsize=figsize,sharex=True,gridspec_kw = {'height_ratios':[3, 1]})
        top = axes[0]
        bottom = axes[1]
        fig.tight_layout()
    else:
        fig = plt.figure(figsize=figsize)
        axes = None
        top = plt
    
    # FIXME: assumes uniform binning
    xc = bins[1:]-binw*0.5
    
    #print mc
    for hist, style in mc:
        pkwargs = copy(mcstyle)
        pkwargs.update(style)
        top.plot(xc+binw*0.5,hist,**pkwargs)
    top.errorbar( xc+binw*0.5, data,ls='None', xerr=np.ones_like(data)*binw*0.5, yerr=np.sqrt(data), color='black', 
                 label='Data', fmt='o', **datastyle )
        
    if axes == None: axes = fig.axes
    
    if ratio:
        ratios = []
        for hist, style in mc:
            rdata = data / hist
            rdata_err = np.sqrt(data) / hist
            ratios.append((rdata,rdata_err))
            rkwargs = {}
            if len(mc) == 1: rkwargs['color'] = 'black'
            elif "color" in style: rkwargs['color'] = style['color']
            rkwargs.update(datastyle)
            bottom.errorbar( xc+binw*0.5, rdata,ls='None', xerr=np.ones_like(rdata)*binw*0.5, yerr=rdata_err, 
                        **rkwargs)
        
        bottom.plot( (bins[0],bins[-1]), (1,1), 'k-' )
        bottom.set_ylabel('Data / MC')
        bottom.set_ylim(0,2)
    
    if logy:
        axes[0].set_yscale('log')
    axes[0].set_xlim(bins[0],bins[-1])
    
    unit = None    
    if var != None:
        if type(var) != str:
            var, unit = var
        if unit: var += " (%s)" % unit
        axes[-1].set_xlabel(var)
    ylabel = 'Events / %1.3g' % binw
    if unit:
        ylabel += ' %s' % unit
    axes[0].set_ylabel(ylabel)

    top.legend(loc='best')
        

#TODO: insert weighted error: check np.histogram function for eveything. 
def plotHistErr(Y, X, Weight=1, Color="r"):
    if (type(Weight) is int): Weight= np.ones(len(Y))
    #Do histograms
    bins, _edges = np.histogram(Y, X, weights= Weight)
    edges        = _edges[:len(_edges)-1]
    bincenters   = 0.5*(_edges[1:]+_edges[:-1])
    #Calculate error normalized in height
    meanStd      = np.sqrt(bins)/float(np.sum(bins))
    #Normalize in height
    bins         = bins/float(np.sum(bins))
    #Plot
    plt.errorbar(bincenters, bins, yerr=meanStd, color=Color, capsize=0)

def getHistErr(Y, X, Weight=1):
    if (type(Weight) is int): Weight= np.ones(len(Y))
    #Do histograms
    bins, _edges = np.histogram(Y, X, weights= Weight)
    edges        = _edges[:len(_edges)-1]
    #Calculate error normalized in height
    bincenters   = 0.5*(_edges[1:]+_edges[:-1])
    meanStd      = np.sqrt(bins)/float(np.sum(bins))
    #Normalize in height
    bins         = bins/float(np.sum(bins))
    return bins, edges, meanStd


#Propagazione degli errori con i pesi???
def plotHistRatio(Y1,Y2,X,Weight1=1,Weight2=1, Color="r"):
    #Do histogram1
    bins1,edges1,meanStd1=getHistErr(Y1,X,Weight1)
    #Do histogram2
    bins2,edges2,meanStd2=getHistErr(Y2,X,Weight2)
    #do the ratio
    ratio = np.true_divide(bins1,bins2)
    #Propagate relative
    ratioErr = 0.5*np.sqrt((meanStd1/bins1)**2+(meanStd2/bins2)**2)*ratio
    #Plot
    plt.errorbar(edges1,ratio, yerr=ratioErr,capsize=0,color=Color)# uncorrected

def getHistRatio(Y1,Y2,X,Weight1=1,Weight2=1, Color="r"):
    #Do histogram1
    bins1,edges1,meanStd1=getHistErr(Y1,X,Weight1)
    #Do histogram2
    bins2,edges2,meanStd2=getHistErr(Y2,X,Weight2)
    #do the ratio
    ratio = np.true_divide(bins1,bins2)
    #Propagate relative
    ratioErr = 0.5*np.sqrt((bins1/meanStd1)**2+(bins2/meanStd2)**2)*ratio
    return ratio, edges1, ratioErr  

def WeightedMean(values, weights): 
        total_weight = weights.sum()
        #solve problem with division by zero in some bins.
        if total_weight==0:return values.mean()
        else: return (values * weights).sum() / float(total_weight)

def WeightedStd(values, weights):
        total_weight = weights.sum()
        if total_weight==0:return values.std()
        else: 
            WeightedMean = (values * weights).sum() / float(total_weight)
            return np.sqrt((weights*(values -WeightedMean)**2).sum() / total_weight)

def Profile(x,y,nbins,xmin,xmax,ymin,ymax, ax, color,weight=1):
    if (type(weight) is int): weight=np.ones(len(x))
    df = DataFrame({'x' : x , 'y' : y, 'weight': weight})   
    df=df.query("x>="+str(xmin)+" & x<="+str(xmax)).query("y>="+str(ymin)+" & y<="+str(ymax))
    binedges = xmin + ((xmax-xmin)/nbins) * np.arange(nbins+1)
    df['bin'] = np.digitize(df['x'],binedges) 
    bincenters = xmin + ((xmax-xmin)/nbins)*np.arange(nbins) + ((xmax-xmin)/(2*nbins))
    ProfileFrame = DataFrame({'bincenters' : bincenters, 'N' : df['bin'].value_counts(sort=False)},index=range(1,nbins+1))
    bins = ProfileFrame.index.values
    for bin in bins:
        y=df.ix[df['bin']==bin,'y']
        w=df.ix[df['bin']==bin,'weight']
        if (type(weight) is int):
            ProfileFrame.ix[bin,'ymean']      = y.mean()
            ProfileFrame.ix[bin,'yStandDev']  = y.std()
        else:
            ProfileFrame.ix[bin,'ymean']      = WeightedMean(y,w)
            ProfileFrame.ix[bin,'yStandDev']  = WeightedStd(y,w)
        ProfileFrame.ix[bin,'yMeanError'] = ProfileFrame.ix[bin,'yStandDev'] / np.sqrt(ProfileFrame.ix[bin,'N'])
    ProfileFrame['ymean']      = ProfileFrame['ymean'].fillna(0)
    ProfileFrame['yMeanError'] = ProfileFrame['yMeanError'].fillna(0)
    ax.errorbar(ProfileFrame['bincenters'], ProfileFrame['ymean'], yerr=np.array(ProfileFrame['yMeanError']), c=color,capsize=0)
    return ax

def exportToPdf(name='sample',section='sample',mode='a'):
    f = open('Plots/'+section+'.tex', mode)
    f.write('\\begin{figure}[H]')
    f.write('\\begin{center}')
    f.write('\\includegraphics[width=1.0\linewidth]{'+name+'.png}')
    f.write('\\end{center}')
    f.write('\\end{figure}')
    f.close()