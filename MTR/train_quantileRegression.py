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

# set the quantiles
quantiles = [ 0.11, 0.22 ] 
# quantiles = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ]
# quantiles = [ 0.05,  0.15,  0.25,  0.35,  0.45,  0.55,  0.65,  0.75,  0.85,  0.95]
# quantiles = [         0.025,  0.05 ,  0.075,  0.1  ,  0.125,  0.15 ,  0.175,
#               0.2  ,  0.225,  0.25 ,  0.275,  0.3  ,  0.325,  0.35 ,  0.375,
#               0.4  ,  0.425,  0.45 ,  0.475,  0.5  ,  0.525,  0.55 ,  0.575,
#               0.6  ,  0.625,  0.65 ,  0.675,  0.7  ,  0.725,  0.75 ,  0.775,
#               0.8  ,  0.825,  0.85 ,  0.875,  0.9  ,  0.925,  0.95 ,  0.975]


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
    qr_data.trainQuantile(q)

    
# MC
# --------------------------------------------------------------------------------
#
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
    qr_mc.trainQuantile(q)
    
