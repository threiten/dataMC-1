from quantileRegression import quantileRegression

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

# data
#--------------------------------------------------------------------------------
#
print mycolors.green+"Training quantile regressions on Data "+mycolors.default

qr_data = quantileRegression("data")

qr_data.loadDF("/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_v1_ichep/",
               "cicNoSigmaIetaIeta/trees/",
               ["Data_13TeV_EBHighR9", "Data_13TeV_EBLowR9", "Data_13TeV_EEHighR9", "Data_13TeV_EELowR9" ],
               100000)

# train n-quantiles regressions
#quantiles = [ 0.25, 0.5, 0.75 ]
quantiles = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ]
for q in quantiles:
    qr_data.trainQuantile(q)

    
# MC
# --------------------------------------------------------------------------------
#

print mycolors.green+"Training quantile regressions on MC "+mycolors.default

qr_mc = quantileRegression( "mc" )

qr_mc.loadDF("/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_v1_mc/",
             "cicNoSigmaIetaIeta/trees/",
             ["DYToEE_powheg_13TeV_EBHighR9", "DYToEE_powheg_13TeV_EBLowR9", "DYToEE_powheg_13TeV_EEHighR9", "DYToEE_powheg_13TeV_EELowR9" ],
             100000)

# train n-quantiles regressions
#quantiles = [ 0.25, 0.5, 0.75 ]
quantiles = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ]
for q in quantiles:
    qr_mc.trainQuantile(q)
    
