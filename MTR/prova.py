print "inside!"
from quantileRegression import quantileRegression
import sys
import os
print "inside the repository, create it"
outputDir = "/mnt/t3nfs01/data01/shome/giulioisac/dataMC/MTR/corrections"
print "create y list"
ylist = ["R9", "S4", "SigmaIeIe", "EtaWidth", "PhiWidth", "CovarianceIphiIphi",'CovarianceIetaIphi']#, "SigmaRR", 'PhoIso03', 'ChIso03', 'ChIso03worst']
# Input variables used for the regression
X = ['Pt', 'ScEta', 'Phi', 'rho']  # <<-- list
print X
# Set the quantiles 
quantiles = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ]
print quantiles
qr_mc1= quantileRegression("mc")
qr_mc1.loadDFh5("/mnt/t3nfs01/data01/shome/giulioisac/dataMC/MTR/df_mc_Final.h5", 0, 9500000)
print "loaded!"
qr_mc1.correctAll(X, ylist, quantiles,EBEE="EE", relativePath= "weights_Period_") 