#match quantiles
import ROOT
from ROOT import TCanvas
from ROOT import TFile
from ROOT import TH1F
import math
import time
import numpy as np
from array import array 




# for testing, set max number of events
nMax = 10000
nBins = 100



# input file with toys
#
file = TFile.Open('/Users/mauro/CMS/Physics/dataMC/toy/dataMC.root', 'read')




# output file with data/mc corrections on matched quantiles
#
nameOutputFile = '/Users/mauro/CMS/Physics/dataMC/knn/corrDataMC.root'
fout = TFile.Open(nameOutputFile, 'recreate')




# load events from the tree into an array
#
start = time.time()
dataPos = []  # position in pt,eta,phi
dataVars = [] # vatriables I'm interested in 
count = 0
R9Min       =   0.8; 
R9Max       =   1.0; 
sieieMin    =   0;
sieieMax    =   0.02;
h_R9MC        = TH1F('h_R9MC'        ,'MC',    nBins, R9Min, R9Max)
h_R9Data      = TH1F('h_R9Data'      ,'data',  nBins, R9Min, R9Max)
h_SieieMC     = TH1F('h_SieieMC'     ,'MC',    nBins, sieieMin, sieieMax)
h_SieieData   = TH1F('h_SieieData'   ,'Data',  nBins, sieieMin, sieieMax)
cdf_R9MC      = TH1F('cdf_R9MC'      ,'MC',    nBins, R9Min, R9Max)
cdf_R9Data    = TH1F('cdf_R9Data'    ,'data',  nBins, R9Min, R9Max)
cdf_SieieMC   = TH1F('cdf_SieieMC'   ,'MC',    nBins, sieieMin, sieieMax)
cdf_SieieData = TH1F('cdf_SieieData' ,'Data',  nBins, sieieMin, sieieMax)
vdata = []
print "Read data"
for event in file.data:
    count += 1
    if (count > nMax): break
    if (count%1000==0): print count
    h_R9Data   .Fill(event.R9)
    h_SieieData.Fill(event.sieie)
    X = [event.pt, event.eta, event.phi, event.R9, event.sieie]
    vdata.append(X)



    
print "Read MC"
count = 0
vmc = []
for eventMC in file.mc:
    count += 1
    if (count > nMax): break
    if (count%1000==0): print count
    h_R9MC   .Fill(eventMC.R9)
    h_SieieMC.Fill(eventMC.sieie)
    X = [eventMC.pt, eventMC.eta, eventMC.phi, eventMC.R9, eventMC.sieie]
    vmc.append(X)
end = time.time()
print "TIME: read data and mc in ", (end-start), "seconds"




# normalize the histrograms and get the cumulative of all distributions
h_R9MC     .Scale(1./h_R9MC     .Integral())
h_R9Data   .Scale(1./h_R9Data   .Integral())
h_SieieMC  .Scale(1./h_SieieMC  .Integral())
h_SieieData.Scale(1./h_SieieData.Integral())

cdf_R9MC      = h_R9MC     .GetCumulative()
cdf_R9Data    = h_R9Data   .GetCumulative()
cdf_SieieMC   = h_SieieMC  .GetCumulative()
cdf_SieieData = h_SieieData.GetCumulative()




# PLOTS
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c = TCanvas("c","c",600,600)
c.Divide(2,2)
c.cd(1)
h_R9MC      .SetLineColor(2)
h_R9MC      .Draw("")
h_R9Data    .SetLineColor(1)
h_R9Data    .Draw("same")
c.cd(2)
h_SieieMC      .SetLineColor(2)
h_SieieMC      .Draw("")
h_SieieData    .SetLineColor(1)
h_SieieData    .Draw("same")
c.cd(3)
cdf_R9MC      .SetLineColor(2)
cdf_R9MC      .Draw("")
cdf_R9Data    .SetLineColor(1)
cdf_R9Data    .Draw("same")
c.cd(4)
cdf_SieieMC      .SetLineColor(2)
cdf_SieieMC      .Draw("")
cdf_SieieData    .SetLineColor(1)
cdf_SieieData    .Draw("same")




# given x on the axis of h1, returns the corresponding value on h2
def mc2data(x, h1, h2): #g --> modify: get the histograms instead of the CDF
    c = h1.GetBinContent(h1.GetXaxis().FindBin(x))
    #print "first bin above", h2.FindFirstBinAbove(c)-1, " ", h2.GetBinCenter(h2.FindFirstBinAbove(c)-1)
    return h2.GetBinCenter(h2.FindFirstBinAbove(c)-1)
    return -999999




#shift MC onto data
#
# create the trees (copy for simplicity also the unchanged data one) to store the new values
tmc   = ROOT.TTree("mcCorr"  , "MC data-MC corrected tree")
tdata = ROOT.TTree("data", "Data tree")
# variables for the branches for mc
bPtmc     = np.zeros(1, dtype=float)
bEtamc    = np.zeros(1, dtype=float)
bPhimc    = np.zeros(1, dtype=float)
bR9mc2    = np.zeros(1, dtype=float)
bSieiemc2 = np.zeros(1, dtype=float)
# create the branches for mc
tmc.Branch('pt'    , bPtmc    , 'bPtmc/D'      )
tmc.Branch('eta'   , bEtamc   , 'bEtamc/D'     )
tmc.Branch('phi'   , bPhimc   , 'bPhimc/D'     )
tmc.Branch('R9'    , bR9mc2   , 'bR9mc2/D'    )
tmc.Branch('sieie' , bSieiemc2, 'bSieiemcm2/D')
# variables for the branches for data
bPtdata    = np.zeros(1, dtype=float)
bEtadata   = np.zeros(1, dtype=float)
bPhidata   = np.zeros(1, dtype=float)
bR9data    = np.zeros(1, dtype=float)
bSieiedata = np.zeros(1, dtype=float)
# create the branches for data
tdata.Branch('pt'   , bPtdata   , 'bPtdata/D'    )
tdata.Branch('eta'  , bEtadata  , 'bEtadata/D'   )
tdata.Branch('phi'  , bPhidata  , 'bPhidata/D'   )
tdata.Branch('R9'   , bR9data   , 'bR9data/D'    )
tdata.Branch('sieie', bSieiedata, 'bSieiedatam/D')

# book  xcheck histograms
h_R9MC2        = TH1F('h_R9MC2'        ,'MC transformed',    nBins, R9Min, R9Max)
h_SieieMC2     = TH1F('h_SieieMC2'     ,'MC transformed',    nBins, sieieMin, sieieMax)
# fill mc histograms and tree
for event in vmc:
    R9mc    = event[3] # index 3 is R9
    Sieiemc = event[4] # index 4 is Sieie
    h_R9MC2   .Fill(mc2data(R9mc   , cdf_R9MC, cdf_R9Data))
    h_SieieMC2.Fill(mc2data(Sieiemc, cdf_SieieMC, cdf_SieieData))
    bPtmc    [0] = event[0]
    bEtamc   [0] = event[1]
    bPhimc   [0] = event[2]
    bR9mc2   [0] = mc2data(R9mc   , cdf_R9MC, cdf_R9Data)
    bSieiemc2[0] = mc2data(Sieiemc, cdf_SieieMC, cdf_SieieData)
    tmc.Fill()
# fill data tree
for event in vdata:
    bPtdata   [0] = event[0]
    bEtadata  [0] = event[1]
    bPhidata  [0] = event[2]
    bR9data   [0] = event[3]
    bSieiedata[0] = event[4]
    tdata.Fill()


tmc.Write()
tdata.Write()

c.cd(1)
h_R9MC2.SetLineColor(3)
h_R9MC2.DrawNormalized("same")

c.cd(2)
h_SieieMC2.SetLineColor(3)
h_SieieMC2.DrawNormalized("okisame")

c.Print("DataMC.png")

# close files
fout.Close()
file.Close()


