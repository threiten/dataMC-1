from quantileRegression import quantileTraining

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

dataTrees =  ["Data_13TeV_EBHighR9",
              "Data_13TeV_EBLowR9",
              "Data_13TeV_EEHighR9",
              "Data_13TeV_EELowR9" ]

qt_data = quantileTraining( "data",
                            "/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_v1_ichep/",
                            "cicNoSigmaIetaIeta/trees/",
                            dataTrees)

nEvts = 1000
qt_data.loadDF(nEvts)

quantiles = [ 0.25, 0.5, 0.75 ]
for q in quantiles:
    qt_data.trainQuantile(q)

    
# MC
# --------------------------------------------------------------------------------
#

print mycolors.green+"Training quantile regressions on MC "+mycolors.default

mcTrees =  ["DYToEE_powheg_13TeV_EBHighR9",
            "DYToEE_powheg_13TeV_EBLowR9",
            "DYToEE_powheg_13TeV_EEHighR9",
            "DYToEE_powheg_13TeV_EELowR9" ]

qt_mc = quantileTraining( "mc",
                          "/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_v1_mc/",
                          "cicNoSigmaIetaIeta/trees/",
                          mcTrees)

nEvts = 1000
qt_mc.loadDF(nEvts)

quantiles = [ 0.25, 0.5, 0.75 ]
for q in quantiles:
    qt_mc.trainQuantile(q)
    
