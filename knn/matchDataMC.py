import ROOT
from ROOT import TCanvas
from ROOT import TFile
from ROOT import TH1F
from ROOT import TH2F
import math
import time
import numpy as np
from array import array 
from numpy import mean, sqrt




# Parameters
#
# number of DATA neighbors per MC point 
nNeighbors = 10




# Euclidean metrics scale factors for eta, phi, pt to give them the same weight
# (just need to roughly match the order of magnitude)
# x -> x/kx
#
keta = 5
kphi = 6.28
kpt  = 200 




# input file with toys
#
file = TFile.Open('/Users/mauro/CMS/Physics/dataMC/toy/dataMC.root', 'read')




# load events from the tree into an array
#
start = time.time()
dataPos = []  # position in pt,eta,phi
dataVars = [] # vatriables I'm interested in 
count = 0
for event in file.data:
    count += 1
    if (count > 1000): break
    X = [event.pt/kpt, event.eta/keta, event.phi/kphi]
    dataPos.append(X)
    Y = [event.R9, event.sieie]
    dataVars.append(Y)
end = time.time()
print "TIME: read the tree  = ", (end-start), "seconds"



    
# Build the KNN model on data
# data will represent the truth of a given MC event
# http://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html#sklearn.neighbors.KNeighborsClassifier.kneighbors
#
start = time.time()
from sklearn.neighbors import NearestNeighbors
neigh = NearestNeighbors(n_neighbors=nNeighbors)
neigh.fit(dataPos) # giving only one sample ito builds only the tree
end = time.time()
print "TIME: build the model = ", (end-start), "seconds"




# for each MC event get the neighbors in DATA from the previous model
# (I will use the neighbors to biuld the "truth from data")
#
count = 0
veta = []
R9Min       =   0.8; 
R9Max       =   1.0; 
sieieMin    =   0;
sieieMax    =   0.02;
h_R9MC       = TH1F('h_R9MC'      ,'MC',    100, R9Min, R9Max)
h_R9Data     = TH1F('h_R9Data'    ,'Data',  100, R9Min, R9Max)
h_R9truth    = TH1F('h_R9truth'   ,'truth', 100, R9Min, R9Max)
h_SieieMC    = TH1F('h_SieieMC'   ,'MC',    100, sieieMin, sieieMax)
h_SieieData  = TH1F('h_SieieData' ,'Data',  100, sieieMin, sieieMax)
h_Sieietruth = TH1F('h_Sieietruth','truth', 100, sieieMin, sieieMax)

h_R9etaMC       = TH2F('h_R9etaMC'      ,'MC',    100, -2.5, 2.5, 100, R9Min, R9Max)
h_R9etaData     = TH2F('h_R9etaData'    ,'Data',  100, -2.5, 2.5, 100, R9Min, R9Max)
h_R9etatruth    = TH2F('h_R9etatruth'   ,'truth', 100, -2.5, 2.5, 100, R9Min, R9Max)
h_SieieetaMC    = TH2F('h_SieieetaMC'   ,'MC',    100, -2.5, 2.5, 100, sieieMin, sieieMax)
h_SieieetaData  = TH2F('h_SieieetaData' ,'Data',  100, -2.5, 2.5, 100, sieieMin, sieieMax)
h_Sieieetatruth = TH2F('h_Sieieetatruth','truth', 100, -2.5, 2.5, 100, sieieMin, sieieMax)

for eventMC in file.mc:
    count += 1
    if (count > 1000): break
    if (count%100==0): print count
    
    MCpos  = [[eventMC.pt/kpt, eventMC.eta/keta, eventMC.phi/kphi]]
    MCvars = [[eventMC.R9, eventMC.sieie]]
    print "MC   = ", MCpos, MCvars
    h_R9etaMC   .Fill(eventMC.eta, eventMC.R9)
    h_SieieetaMC.Fill(eventMC.eta, eventMC.sieie)
    h_R9MC   .Fill(eventMC.R9)
    h_SieieMC.Fill(eventMC.sieie)
    
    distEvt = neigh.kneighbors(MCpos) # Returns indices of and distances to the N neighbors of each point.
    dists = distEvt[0]   # distances
    vidx  = distEvt[1]   # indices

    meanR9 = 0
    meanSieie = 0
    vR9    = []
    vSieie = []
    for i in vidx[0]:  # 0 = R9, 1 = sieie
        print "DATA = ", i, dataPos[i], dataVars[i]
        R9sieie = dataVars[i]
        vR9.append(R9sieie[0])
        vSieie.append(R9sieie[1])
        h_R9etaData   .Fill(eventMC.eta, eventMC.R9)
        h_SieieetaData.Fill(eventMC.eta, eventMC.sieie)
        h_R9Data   .Fill(R9sieie[0])
        h_SieieData.Fill(R9sieie[1])
        
    meanR9    = np.average(vR9)
    meanSieie = np.average(vSieie)
    print "mean R9, sieie =                                                        =  ", meanR9, meanSieie
    h_R9etatruth   .Fill(eventMC.eta, eventMC.R9)
    h_Sieieetatruth.Fill(eventMC.eta, eventMC.sieie)
    h_R9truth   .Fill(meanR9)
    h_Sieietruth.Fill(meanSieie)

    
# PLOTS
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c = TCanvas("c","c",600,600)
c.Divide(2,2)
c.cd(1)
h_R9truth   .SetLineColor(1)
h_R9truth   .SetLineStyle(2)
h_R9truth   .DrawNormalized("")
h_R9MC      .SetLineColor(2)
h_R9MC      .DrawNormalized("same")
h_R9Data    .SetLineColor(1)
h_R9Data    .DrawNormalized("same")
c.cd(2)
h_Sieietruth   .SetLineColor(1)
h_Sieietruth   .SetLineStyle(2)
h_Sieietruth   .DrawNormalized("")
h_SieieMC      .SetLineColor(2)
h_SieieMC      .DrawNormalized("same")
h_SieieData    .SetLineColor(1)
h_SieieData    .DrawNormalized("same")
c.cd(3)
h_R9etatruth   .SetLineColor(1)
h_R9etatruth   .SetLineStyle(2)
h_R9etatruth   .Draw("")
h_R9etaMC      .SetLineColor(2)
h_R9etaMC      .Draw("same")
h_R9etaData    .SetLineColor(1)
h_R9etaData    .Draw("same")
c.cd(4)
h_Sieieetatruth   .SetLineColor(1)
h_Sieieetatruth   .SetLineStyle(2)
h_Sieieetatruth   .Draw("")
h_SieieetaMC      .SetLineColor(2)
h_SieieetaMC      .Draw("same")
h_SieieetaData    .SetLineColor(1)
h_SieieetaData    .Draw("same")
c.Print("dataMC.png")
