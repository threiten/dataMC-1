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

mcTrees =  [#"DYToEE_powheg_13TeV_EBHighR9",
            #"DYToEE_powheg_13TeV_EBLowR9",
            #"DYToEE_powheg_13TeV_EEHighR9",
            "DYToEE_powheg_13TeV_EELowR9" ]

qr_mc = quantileRegression( "mc",
                            "/Users/mauro/CMS/Physics/dataMC/nt/double_ele_spring16v2_sync_v1_mc/",
                            "cicNoSigmaIetaIeta/trees/",
                            mcTrees)


print "Load the mc"
nEvts = 3
qr_mc.loadDF(nEvts)

# get MC and data weights
#quantiles = [ 0.25, 0.5, 0.75 ]
quantiles = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 ]
filename = "./weights/mc_weights"
qr_mc.getWeights("mc", filename, quantiles)
filename = "./weights/data_weights"
qr_mc.getWeights("data", filename, quantiles)

print "Get the corrected data-mc values"
X = ['Pt', 'ScEta', 'Phi', 'rho']
y = "R9"
qr_mc.correctY(X, y, quantiles)
y_in   = qr_mc.getY(y)
y_corr = qr_mc.getCorrectedY(y)
print y_in
print y_corr
