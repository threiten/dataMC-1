import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pandas import Series, DataFrame
plt.style.use('ggplot') 

#TODO: insert weighted error
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
        total_weight = (weights * weights).sum()
        if total_weight==0: return values.std()
        else: return 1.0/ np.sqrt(float(total_weight))

def Profile(x,y,nbins,xmin,xmax,ax, color,weight=1):
    if (type(weight) is int): weight=np.ones(len(x))
    df = DataFrame({'x' : x , 'y' : y, 'weight': weight})    
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