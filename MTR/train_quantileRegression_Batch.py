import sys
import os
from import_file import import_file

QReg = import_file("/mnt/t3nfs01/data01/shome/threiten/QReg/dataMC-1/MTR/quantileRegression")

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
sEBEE      = sys.argv[8]

# qr = quantileRegression(sys.argv[1])

outputDir = "/mnt/t3nfs01/data01/shome/threiten/QReg/ReReco17_data/camp_3_0_X/weights"
#cols=['probePt', 'probePhi', 'probeScEta', 'rho', 'weight', 'run', 'probeR9', 'probeS4', 'probeEtaWidth', 'probePhiWidth', 'probeSigmaIeIe', 'probeCovarianceIetaIphi','probeCovarianceIphiIphi']
if not os.path.exists(outputDir):
   print 'Creating output dir:', outputDir
   os.mkdir(outputDir)

if dataMC == "data":

   for q in quantiles:

      qr = QReg.quantileRegression(sys.argv[1])
      # to reduce memory consuption just load the locally pre-made h5 file
      qr.loadDFh5("/mnt/t3nfs01/data01/shome/threiten/QReg/ReReco17_data/camp_3_0_X/df_data_kinSS_All.h5", startEvt, stopEvt, 12345, False)
      qr.trainQuantile(Y, q, outputDir, EBEE = sEBEE, maxDepth = imaxDepth, minLeaf = iminLeaf,  useWeights = False)

elif dataMC == "mc":
   for q in quantiles:
      qr = QReg.quantileRegression(sys.argv[1])
      # to reduce memory consuption just load the locally pre-made h5 file
      qr.loadDFh5("/mnt/t3nfs01/data01/shome/threiten/QReg/ReReco17_data/camp_3_0_X/df_mc_kinSS_All.h5", startEvt, stopEvt, 12345, False)
      qr.trainQuantile(Y, q, outputDir, EBEE = sEBEE, maxDepth = imaxDepth, minLeaf = iminLeaf,  useWeights = False)


else: print " ERROR: choose data or mc"
