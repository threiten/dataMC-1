import sys
from quantileRegression import quantileRegression
import pandas as pd
import numpy as np

# Nevt passing selection
# mc   4545316 : train on 2000000
# data 1441902 : train on all

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


# set the quantiles
quantiles = [ 0.11] 

Y = "R9" #  "R9", "S4", "SigmaIeIe", "EtaWidth", "PhiWidth", "CovarianceIphiIphi", "SigmaRR"]

outputDir = "/Users/mauro/CMS/Physics/dataMC/MTR/weightsTest/"

# data
#--------------------------------------------------------------------------------
#
print mycolors.green+"Training quantile regressions on Data "+mycolors.default

qr_data = quantileRegression("data")
startEvt = 1
stopEvt   = 1000
qr_data.loadDF("/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_v1_ichep/",
               "cicNoSigmaIetaIeta/trees/",
               ["Data_13TeV_EBHighR9", "Data_13TeV_EBLowR9", "Data_13TeV_EEHighR9", "Data_13TeV_EELowR9" ],
               startEvt, stopEvt, 12345)

print "Data Frame has ", qr_data.getNEntries(), "entries"

for q in quantiles:
#    qr_data.trainQuantile(Y, q, outputDir, "EB", 3, 9)
#    qr_data.trainQuantile(Y, q, outputDir, "EE", 3, 9)
    qr_data.trainQuantile(Y, q, outputDir, ""  , 3, 9)



# MC
# --------------------------------------------------------------------------------
#
print mycolors.green+"Training quantile regressions on MC "+mycolors.default

qr_mc = quantileRegression( "mc" )
startEvt = 1
stopEvt   = 1000
qr_mc.loadDF("/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_v1_mc/",
             "cicNoSigmaIetaIeta/trees/",
             ["DYToEE_powheg_13TeV_EBHighR9", "DYToEE_powheg_13TeV_EBLowR9", "DYToEE_powheg_13TeV_EEHighR9", "DYToEE_powheg_13TeV_EELowR9" ],
             startEvt, stopEvt, 12345)

print "Data Frame has ", qr_mc.getNEntries(), "entries"

for q in quantiles:
#    qr_mc.trainQuantile(Y, q, outputDir, "EB", 3, 9)
#    qr_mc.trainQuantile(Y, q, outputDir, "EE", 3, 9)
    qr_mc.trainQuantile(Y, q, outputDir, ""  , 3, 9)

