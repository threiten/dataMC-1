import ROOT
from ROOT import TCanvas
from ROOT import TFile
from ROOT import TH1F, TH2F

from root_numpy import root2array, tree2array
from root_numpy import fill_hist

import numpy as np

import matplotlib.pyplot as plt

from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor

import pandas as pd

import root_pandas as rpd
from root_pandas import read_root

import time

from sklearn.externals import joblib

class mycolors:
   red = '\033[91m'
   green = '\033[92m'
   blue = '\033[94m'
   yan = '\033[96m'
   cWhite = '\033[97m'
   yellow = '\033[93m'
   magenta = '\033[95m'
   grey = '\033[90m'
   black = '\033[90m'
   default = '\033[0m'
   










   
#
# Prepare the data
# --------------------------------------------------------------------------------
#
   
# mc spigazzi's structure
# inputDir = "/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_MC_05JulyEcorr_cert/"
# fname    = inputDir+"output.root"
# dir      = "cic/trees/"


# mc
#
# inputDir = "/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_v1_mc/"
# dir      = "cicNoSigmaIetaIeta/trees/"
# trees =  [dir+"DYToEE_powheg_13TeV_EBHighR9",
#           dir+"DYToEE_powheg_13TeV_EBLowR9",
#           dir+"DYToEE_powheg_13TeV_EEHighR9",
#           dir+"DYToEE_powheg_13TeV_EELowR9" ]


# data
#
inputDir = "/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_v1_ichep/"
dir      = "cicNoSigmaIetaIeta/trees/"
trees =  [dir+"Data_13TeV_EBHighR9",
          dir+"Data_13TeV_EBLowR9",
          dir+"Data_13TeV_EEHighR9",
          dir+"Data_13TeV_EELowR9" ]


# common to data and mc
#
fname    = inputDir+"output.root"

# branches
#
evtBranches         = ["rho", "nvtx"]

trgBranches         = [ "leadHLT_Ele27_WPTight_Gsf_vMatch", "subleadHLT_Ele27_WPTight_Gsf_vMatch" ]

recoLeadBranches    = ["leadPt", "leadScEta", "leadPhi",
                       "leadR9", "leadS4", "leadSigmaIeIe", "leadEtaWidth", "leadPhiWidth", "leadCovarianceIphiIphi", "leadSigmaRR"]

recoSubleadBranches = ["subleadPt", "subleadScEta", "subleadPhi",
                       "subleadR9", "subleadS4", "subleadSigmaIeIe", "subleadEtaWidth", "subleadPhiWidth", "subleadCovarianceIphiIphi", "subleadSigmaRR"]

recoBranches        = evtBranches  + trgBranches + recoLeadBranches + recoSubleadBranches
print "NT branches: "
print recoBranches

# use common names for the traning dataset
#
uniformColumnsNames = ["rho", "nvtx" ,                       
                        "Pt", "ScEta", "Phi",
                        "R9", "S4", "SigmaIeIe", "EtaWidth", "PhiWidth", "CovarianceIphiIphi", "SigmaRR"]

# read the trees
#
print "Reading trees"
for t in trees:
   df = rpd.read_root(fname,t,columns=recoBranches)
print "Number of events  "
print df.count()

# trigger matching logic :
# get the indices of the lead matching the trigger and use them to select the sublead
# get the indices of the sublead matching the trigger and use them to select the lead
# NB: I can use the event twice if both lead and sublead trigger
#
# the trigger decision is stored as a float -> convert it to bool
# print df["leadHLT_Ele27_WPTight_Gsf_vMatch"].sum()
# print df["subleadHLT_Ele27_WPTight_Gsf_vMatch"].sum()
df.leadHLT_Ele27_WPTight_Gsf_vMatch    = df.leadHLT_Ele27_WPTight_Gsf_vMatch.astype(bool)
df.subleadHLT_Ele27_WPTight_Gsf_vMatch = df.subleadHLT_Ele27_WPTight_Gsf_vMatch.astype(bool)
idx_trgLead    = df[df["leadHLT_Ele27_WPTight_Gsf_vMatch"]].index.tolist()
idx_trgSublead = df[df["subleadHLT_Ele27_WPTight_Gsf_vMatch"]].index.tolist()
print "# lead trig   = ", len(idx_trgLead)
print "# sublead trg = ", len(idx_trgSublead)


# if the lead triggered use the sublead
#
dataSublead = df[ evtBranches + recoSubleadBranches ]
dataSublead = dataSublead.loc[idx_trgLead]
dataSublead.columns = uniformColumnsNames
#print dataSublead.count()
#
# if the sublead triggered use the lead
#
dataLead    = df[ evtBranches + recoLeadBranches ]
dataLead    = dataLead.loc[idx_trgSublead]
dataLead.columns = uniformColumnsNames
#print dataLead.count()


# concatenate leading and subleadind
frames = [dataLead, dataSublead]
data = pd.concat(frames)
#print data.count()
# reset the rows indexing
df = data.reset_index()
# print data


# reshuffle events
#
print "Reshuffle events"
rndseed = 12345
np.random.seed(rndseed)
df['random_index'] = np.random.permutation(range(df.index.size))
df.sort_values(by='random_index',inplace=True)
df.set_index('random_index',inplace=True)
# print df

maxEvents = 5

# apply basic selection
#
ptmin  =  20.
ptmax  =  60.
etamin = -2.5
etamax =  2.5
phimin = -3.14
phimax =  3.14
df = df.query('@ptmin < Pt and Pt < @ptmax and @etamin < ScEta and ScEta < @etamax and @phimin < Phi and Phi < @phimax')

print mycolors.green+"Data Frame with nEvt = "+mycolors.default, maxEvents
df = df[0:maxEvents]

# quantile regressions features
X     = df.ix[:,['Pt', 'ScEta', 'Phi', 'rho']]
# targets
R9                 = df['R9']
S4                 = df["S4"]
sieie              = df["SigmaIeIe"]
EtaWidth           = df["EtaWidth"]
PhiWidth           = df["PhiWidth"]
CovarianceIphiIphi = df["CovarianceIphiIphi"]
SigmaRR            = df["SigmaRR"]           

print X
print R9
print sieie
print S4


# train first step in SST
#
print mycolors.green+"Training X,R9 = "+mycolors.default

clf_R9 = GradientBoostingRegressor(loss='ls',
                                   n_estimators=250, max_depth=3,
                                   learning_rate=.1, min_samples_leaf=9,
                                   min_samples_split=9)
clf_R9.fit(X, R9)
#
#
clf_sieie = GradientBoostingRegressor(loss='ls', 
                                      n_estimators=250, max_depth=3,
                                      learning_rate=.1, min_samples_leaf=9,
                                      min_samples_split=9)
clf_sieie.fit(X, sieie)
#
#
clf_S4 = GradientBoostingRegressor(loss='ls', 
                                   n_estimators=250, max_depth=3,
                                   learning_rate=.1, min_samples_leaf=9,
                                   min_samples_split=9)
clf_S4.fit(X, S4)
#
#

# meta variables generation
yhat_R9    = pd.DataFrame(clf_R9.predict(X)   , columns=["yhat_R9"]    )
yhat_sieie = pd.DataFrame(clf_sieie.predict(X), columns=["yhat_sieie"] )
yhat_S4    = pd.DataFrame(clf_S4   .predict(X), columns=["yhat_S4"]    )

# second training stage
print X
print yhat_R9
print yhat_sieie
print yhat_S4

frames = [X, yhat_R9, yhat_sieie, yhat_S4 ]
metaBLock = pd.concat(frames, axis=1)
print metaBLock


clf2_R9 = GradientBoostingRegressor(loss='ls',
                                   n_estimators=250, max_depth=3,
                                   learning_rate=.1, min_samples_leaf=9,
                                   min_samples_split=9)
clf2_R9.fit(metaBLock, R9)
#
#
clf2_sieie = GradientBoostingRegressor(loss='ls', 
                                      n_estimators=250, max_depth=3,
                                      learning_rate=.1, min_samples_leaf=9,
                                      min_samples_split=9)
clf2_sieie.fit(metaBLock, sieie)
#
#
clf2_S4 = GradientBoostingRegressor(loss='ls', 
                                   n_estimators=250, max_depth=3,
                                   learning_rate=.1, min_samples_leaf=9,
                                   min_samples_split=9)
clf2_S4.fit(metaBLock, S4)
