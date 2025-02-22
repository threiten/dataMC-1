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
imaxDepth  = int(sys.argv[6])
iminLeaf   = int(sys.argv[7])




outputDir = "/mnt/t3nfs01/data01/shome/giulioisac/dataMC/MTR/weights_MaxDepth_" + str(imaxDepth) + "_minLeaf_" + str(iminLeaf)
if (not os.path.exists(outputDir)):
   print 'Creating output dir:', outputDir
   os.mkdir(outputDir)

if dataMC == "data":
   for q in quantiles:
      qr = quantileRegression(sys.argv[1])
      qr.loadDFh5("/mnt/t3nfs01/data01/shome/giulioisac/dataMC/MTR/df_data_0-2000000.h5", startEvt, stopEvt)
      qr.trainQuantile(Y, q, outputDir, maxDepth = imaxDepth, minLeaf = iminLeaf,  useWeights = False)

elif dataMC == "mc":
   for q in quantiles:
      qr = quantileRegression(sys.argv[1])
      qr.loadDFh5("/mnt/t3nfs01/data01/shome/giulioisac/dataMC/MTR/df_mc_0-2000000.h5", startEvt, stopEvt)
      qr.trainQuantile(Y, q, outputDir, maxDepth = imaxDepth, minLeaf = iminLeaf,  useWeights = False)


else: print " ERROR: choose data or mc"



    
