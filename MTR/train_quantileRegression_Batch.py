import sys
import os
from quantileRegression import quantileRegression

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


print mycolors.green+"Training quantile regressions on"+mycolors.default
dataMC = sys.argv[1]
print "Data/MC = ", dataMC
Y = sys.argv[2]
print "Y = ", Y
quantiles  = [ float(sys.argv[3]) ]
print "Quantile = ", quantiles
startEvt = int(sys.argv[4])
stopEvt  = int(sys.argv[5])
print "Events in [", startEvt, ", ", stopEvt, "]"
EBEE      = int(sys.argv[6])
maxDepth  = int(sys.argv[7])
minLeaf   = int(sys.argv[8])

qr = quantileRegression(sys.argv[1])

outputDir = "/mnt/t3nfs01/data01/shome/mdonega/dataMC/MTR/weights_MaxDepth_" + str(maxDepth) + "_minLeaf_" + str(minLeaf)
if (not os.path.exists(outputDir)):
   print 'Creating output dir:', outputDir
   os.mkdir(outputDir)

if dataMC == "data":
   qr.loadDF("/mnt/t3nfs01/data01/shome/mdonega/dataMC/nt/double_ele_spring16v2_sync_v1_ichep/",
             "cicNoSigmaIetaIeta/trees/",
             ["Data_13TeV_EBHighR9", "Data_13TeV_EBLowR9", "Data_13TeV_EEHighR9", "Data_13TeV_EELowR9" ],
             startEvt, stopEvt, 12345)
   for q in quantiles:
      qr.trainQuantile(Y, q, outputDir, EBEE, maxDepth, minLeaf)

elif dataMC == "mc":
   qr.loadDF("/mnt/t3nfs01/data01/shome/mdonega/dataMC/nt/double_ele_spring16v2_sync_v1_mc/",
             "cicNoSigmaIetaIeta/trees/",
             ["DYToEE_powheg_13TeV_EBHighR9", "DYToEE_powheg_13TeV_EBLowR9", "DYToEE_powheg_13TeV_EEHighR9", "DYToEE_powheg_13TeV_EELowR9" ],
             startEvt, stopEvt, 12345)
   for q in quantiles:
      qr.trainQuantile(Y, q, outputDir, EBEE, maxDepth, minLeaf)

else: print " ERROR: choose data or mc"



    
