import numpy as np
import matplotlib.pyplot as plt
from root_numpy import root2array, tree2array
from root_numpy import fill_hist
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.ensemble import GradientBoostingRegressor
import pandas as pd
import root_pandas as rpd
from   root_pandas import read_root
import time
import pickle
import gzip

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
# 
# --------------------------------------------------------------------------------
#
class quantileRegression:

   def __init__(self, dmc, iD, tD, t):
      self.dataMC   = dmc
      self.inputDir = iD
      self.treeDir  = tD
      self.trees    = t 
      

   fname    = "output.root"

   evtBranches         = ["rho", "nvtx"]

   trgBranches         = [ "leadHLT_Ele27_WPTight_Gsf_vMatch", "subleadHLT_Ele27_WPTight_Gsf_vMatch" ]
   
   recoLeadBranches    = ["leadPt", "leadScEta", "leadPhi",
                          "leadR9", "leadS4", "leadSigmaIeIe", "leadEtaWidth", "leadPhiWidth", "leadCovarianceIphiIphi", "leadSigmaRR"]
   
   recoSubleadBranches = ["subleadPt", "subleadScEta", "subleadPhi",
                          "subleadR9", "subleadS4", "subleadSigmaIeIe", "subleadEtaWidth", "subleadPhiWidth", "subleadCovarianceIphiIphi", "subleadSigmaRR"]
   
   data_recoBranches   = evtBranches  + trgBranches + recoLeadBranches + recoSubleadBranches
   mc_recoBranches     = evtBranches  + recoLeadBranches + recoSubleadBranches

   df = 0

   


   
   # load the dataframe from the input files
   # 
   # --------------------------------------------------------------------------------
   #
   def loadDF(self, maxEvents = -1):
      
      fname    = self.inputDir+self.fname
      
      print "NT branches: "
      recoBranches = self.data_recoBranches 
      if self.dataMC == "mc" :
         recoBranches = self.mc_recoBranches
      print recoBranches
            
      # use common names for the traning dataset
      #
      uniformColumnsNames = ["rho", "nvtx" ,                       
                             "Pt", "ScEta", "Phi",
                             "R9", "S4", "SigmaIeIe", "EtaWidth", "PhiWidth", "CovarianceIphiIphi", "SigmaRR"]
      
      # attach the tree structure to the tree names
      #
      trees = []
      for t in self.trees:
         trees.append(self.treeDir+t)
      self.trees = trees
      # print trees

      # read the trees into a dataframe
      #
      print "Reading trees"
      for t in self.trees:
         df = rpd.read_root(fname,t,columns=recoBranches)
      print "Number of events  "
      print df.count()

      if (self.dataMC == "data" ):
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
         dataSublead = df[ self.evtBranches + self.recoSubleadBranches ]
         dataSublead = dataSublead.loc[idx_trgLead]
         dataSublead.columns = uniformColumnsNames
         #print dataSublead.count()
         #
         # if the sublead triggered use the lead
         #
         dataLead    = df[ self.evtBranches + self.recoLeadBranches ]
         dataLead    = dataLead.loc[idx_trgSublead]
         dataLead.columns = uniformColumnsNames
         #print dataLead.count()
         
         
         # concatenate leading and subleadind
         frames = [dataLead, dataSublead]
         data = pd.concat(frames)
         #print data.count()
         # reset the rows indexing
         df = data.reset_index()
         #print data

         
      if (self.dataMC == "mc" ):
         # no trigger matching on mc
         # use both lead and sublead
         dataSublead = df[ self.evtBranches + self.recoSubleadBranches ]
         dataSublead.columns = uniformColumnsNames
         #print dataSublead.count()
         #
         # if the sublead triggered use the lead
         #
         dataLead    = df[ self.evtBranches + self.recoLeadBranches ]
         dataLead.columns = uniformColumnsNames
         #print dataLead.count()
         
         
         # concatenate leading and subleadind
         frames = [dataLead, dataSublead]
         data = pd.concat(frames)
         #print data.count()
         # reset the rows indexing
         df = data.reset_index()
         #print data

         
      
      # reshuffle events
      #
      print "Reshuffle events"
      rndseed = 12345
      np.random.seed(rndseed)
      df['random_index'] = np.random.permutation(range(df.index.size))
      df.sort_values(by='random_index',inplace=True)
      df.set_index('random_index',inplace=True)
      # print df
      
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
      if (maxEvents != -1):
         print "Running on ", maxEvents
         df = df[0:maxEvents]
      else:
         print "Running on all events"
      
      self.df = df





   # run the trainings
   # 
   # --------------------------------------------------------------------------------
   #
   def trainQuantile(self, alpha):

      # quantile regressions features
      X     = self.df.ix[:,['Pt', 'ScEta', 'Phi', 'rho']]
      # target
      R9    = self.df['R9']

      # train quantile regression
      #
      print mycolors.green+"Train q = "+mycolors.default, alpha
      clf = GradientBoostingRegressor(loss='quantile', alpha=alpha,
                                      n_estimators=250, max_depth=3,
                                      learning_rate=.1, min_samples_leaf=9,
                                      min_samples_split=9)
      t0 = time.time()
      clf.fit(X, R9)
      t1 = time.time()
      print " time = ", t1-t0
      print "Predict"
      t0 = time.time()
      y_upper = clf.predict(X)
      t1 = time.time()
      print " time = ", t1-t0
      y_upper1d=y_upper.ravel()
      print "Save weights"
      outputName = ""
      if (self.dataMC == "data"):
         outputName = "./weights/data_weights_" + str(alpha) + ".pkl"
      else :
         outputName = "./weights/mc_weights_" + str(alpha) + ".pkl"         
      pickle.dump(clf, gzip.open(outputName, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)





#MDDB #
#MDDB # 
#MDDB # --------------------------------------------------------------------------------
#MDDB #
#MDDB class quantileApply:      
#MDDB 
#MDDB    def getCorrectionta(self, y):
#MDDB 
#MDDB       # quantile regressions features
#MDDB       X     = self.df.ix[:,['Pt', 'ScEta', 'Phi', 'rho']]
#MDDB       # target
#MDDB       R9    = self.df['R9']
#MDDB 
#MDDB       # Read the weights
#MDDB       mcWeights   = "./weights/data_weights_" + str(alpha) + ".pkl"
#MDDB       dataWeights = "./weights/mc_weights_"   + str(alpha) + ".pkl"         
#MDDB 
#MDDB       mcclf   = pickle.load(gzip.open(mcWeights))
#MDDB       dataclf = pickle.load(gzip.open(dataWeights))
#MDDB 
#MDDB       t0 = time.time()
#MDDB       print "Predict"
#MDDB       y_mc = clf.predict(X)
#MDDB       t1 = time.time()
#MDDB       print " time = ", t1-t0
#MDDB       y_upper1d=y_upper.ravel()
#MDDB       print "Save weights"
#MDDB       outputName = ""
#MDDB       if (self.dataMC == "data"):
#MDDB          outputName = "./weights/data_weights_" + str(alpha) + ".pkl"
#MDDB       else :
#MDDB          outputName = "./weights/mc_weights_" + str(alpha) + ".pkl"         
#MDDB       pickle.dump(clf, gzip.open(outputName, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
#MDDB       ## # To retrieve:
#MDDB       ## clf = pickle.load(gzip.open("test.pkl"))

    


## #
## # PLOTS
## # --------------------------------------------------------------------------------
## #
## 
## # values used to produce scatter plots
## pt    = df.ix[:,['Pt']]
## eta   = df.ix[:,['ScEta']]
## phi   = df.ix[:,['Phi']]
## rho   = df.ix[:,['rho']]
## 
## print "Plotting"
## 
## # projection var vs eta
## def projection(nbins, min, max, varbins, var, y_lower1d, y_pred1d, y_upper1d):
##     bins  = np.linspace(min, max,nbins+1)
##     for i in range(0,nbins):
##         varbins.append(0.5*(bins[i]+bins[i+1]))
##     inds  = np.digitize(var, bins)    
##     varR9_upper = np.zeros(nbins)
##     varR9       = np.zeros(nbins)
##     varR9_lower = np.zeros(nbins)
##     nvar        = np.zeros(nbins)
##     var1d       = np.ravel(var)
##     for n in range(var.size):
##        # print inds[n]-1 , "   " , bins[inds[n]-1], "<=", var1d[n], "<", bins[inds[n]]
##        varR9_upper[inds[n]-1] = varR9_upper[inds[n]-1] + y_upper1d[n]
##        varR9      [inds[n]-1] = varR9      [inds[n]-1] + y_pred1d[n]
##        varR9_lower[inds[n]-1] = varR9_lower[inds[n]-1] + y_lower1d[n]
##        nvar[inds[n]-1] += 1
## 
##     meanR9var_upper = varR9_upper / nvar
##     meanR9var       = varR9       / nvar
##     meanR9var_lower = varR9_lower / nvar
##     return varbins, meanR9var_lower, meanR9var, meanR9var_upper
## 
## # R9 vs pt
## nbins = 10
## meanR9pt_lower = []
## meanR9pt       = [] 
## meanR9pt_upper = [] 
## ptbins         = []
## ptbins, meanR9pt_lower, meanR9pt, meanR9pt_upper = projection(nbins, ptmin, ptmax, ptbins, pt, y_lower1d, y_pred1d, y_upper1d)
## 
## # R9 vs eta
## nbins = 10
## meanR9eta_lower = []
## meanR9eta       = [] 
## meanR9eta_upper = [] 
## etabins         = []
## etabins, meanR9eta_lower, meanR9eta, meanR9eta_upper = projection(nbins, etamin, etamax, etabins, eta, y_lower1d, y_pred1d, y_upper1d)
## 
## # R9 vs phi
## nbins = 10
## meanR9phi_lower = []
## meanR9phi       = [] 
## meanR9phi_upper = [] 
## phibins         = []
## phibins, meanR9phi_lower, meanR9phi, meanR9phi_upper = projection(nbins, phimin, phimax, phibins, phi, y_lower1d, y_pred1d, y_upper1d)
## 
## 
## # 2D plots
## #
## 
## def plot2D(xlabel, ylabel, Xvar, Yvar, varbins, meanYvar_upper, meanYvar, meanYvar_lower, outfile ):
##     fig  = plt.figure()
##     fig.patch.set_facecolor('white')
##     axes = plt.axes()
##     plt.grid()
##     plt.xlabel(xlabel)
##     plt.ylabel(ylabel)
##     axes.xaxis.set_label_coords(1., -0.07)
##     axes.yaxis.set_label_coords(-0.07, 1.)
##     plt.scatter(Xvar   ,Yvar           , color='grey'  , marker='+')
##     plt.scatter(varbins,meanYvar_upper , color='blue'  , marker="H", linestyle='-')
##     plt.plot   (varbins,meanYvar_upper , linestyle='-')
##     plt.scatter(varbins,meanYvar       , color='green' , marker="H", linestyle='-')
##     plt.plot   (varbins,meanYvar       , linestyle='-')
##     plt.scatter(varbins,meanYvar_lower , color='red'   , marker="H", linestyle='-')
##     plt.plot   (varbins,meanYvar_lower , linestyle='-')
##     plt.plot()
##     # plt.show()
##     fig.savefig(outfile)
## 
## plot2D('$pt$'  , 'R$_9$', pt  ,R9, ptbins,  meanR9pt_upper,  meanR9pt,  meanR9pt_lower,  "./meanR9_pt.pdf")
## plot2D('$\eta$', 'R$_9$', eta ,R9, etabins, meanR9eta_upper, meanR9eta, meanR9eta_lower, "./meanR9_eta.pdf")
## plot2D('$\phi$', 'R$_9$', phi ,R9, phibins, meanR9phi_upper, meanR9phi, meanR9phi_lower, "./meanR9_phi.pdf")
## 
## 
## 
## 
## # not finished 
## # not finished #
## # not finished # Set up the traning for N-quantiles
## # not finished # --------------------------------------------------------------------------------
## # not finished #
## # not finished 
## # not finished maxEvents = 1000
## # not finished 
## # not finished # apply basic selection
## # not finished #
## # not finished ptmin  =  20.
## # not finished ptmax  =  60.
## # not finished etamin = -2.5
## # not finished etamax =  2.5
## # not finished phimin = -3.14
## # not finished phimax =  3.14
## # not finished df = df.query('@ptmin < Pt and Pt < @ptmax and @etamin < ScEta and ScEta < @etamax and @phimin < Phi and Phi < @phimax')
## # not finished 
## # not finished print mycolors.green+"Data Frame with nEvt = "+mycolors.default, maxEvents
## # not finished df = df[0:maxEvents]
## # not finished 
## # not finished # quantile regressions features
## # not finished X     = df.ix[:,['Pt', 'ScEta', 'Phi', 'rho']]
## # not finished # target
## # not finished R9    = df['R9']
## # not finished 
## # not finished # train quantile regression
## # not finished # e.g. 3 --> 0.25, 0.5, 0.75 (0 and 1 are not allowed)
## # not finished nquantiles = 9
## # not finished #
## # not finished y = []
## # not finished 
## # not finished for iq  in range(1,nquantiles+1): # 0.1,0.2, ... , 0.9
## # not finished    alpha = iq*1./(nquantiles+1)
## # not finished    
## # not finished    print mycolors.green+"Train q = "+mycolors.default, alpha
## # not finished    clf = GradientBoostingRegressor(loss='quantile', alpha=alpha,
## # not finished                                    n_estimators=250, max_depth=3,
## # not finished                                    learning_rate=.1, min_samples_leaf=9,
## # not finished                                    min_samples_split=9)
## # not finished    t0 = time.time()
## # not finished    clf.fit(X, R9)
## # not finished    t1 = time.time()
## # not finished    print " time = ", t1-t0
## # not finished    print "Predict"
## # not finished    t0 = time.time()
## # not finished    y.append(clf.predict(X))
## # not finished    t1 = time.time()
## # not finished    print " time = ", t1-t0
## # not finished 
## # not finished #for iq  in range(1,nquantiles+1):
## # not finished #   print iq, y[iq-1]
## # not finished    
## # not finished 
## # not finished #
## # not finished # PLOTS
## # not finished # --------------------------------------------------------------------------------
## # not finished #
## # not finished 
## # not finished # values used to produce scatter plots
## # not finished pt    = df.ix[:,['Pt']]
## # not finished eta   = df.ix[:,['ScEta']]
## # not finished phi   = df.ix[:,['Phi']]
## # not finished rho   = df.ix[:,['rho']]
## # not finished 
## # not finished print "Plotting"
## # not finished 
## # not finished # projection data and predictions vs var
## # not finished def projection(nbins, min, max, varbins, var, y, nquantiles):
## # not finished    bins  = np.linspace(min, max,nbins+1)
## # not finished    # bin variable var
## # not finished    for i in range(0,nbins):
## # not finished        varbins.append(0.5*(bins[i]+bins[i+1]))
## # not finished    inds  = np.digitize(var, bins)
## # not finished    varTarget = []
## # not finished    y1d = []
## # not finished    nvar = []
## # not finished    for iq in range (1,nquantiles+1):
## # not finished       varTarget   .append(np.zeros(nbins))
## # not finished       y1d         .append(y[iq-1].ravel())
## # not finished       nvar        .append(np.zeros(nbins))
## # not finished    var1d       = np.ravel(var)
## # not finished    # projections
## # not finished    meanTargetvar = []
## # not finished    for n in range(var.size):
## # not finished       print n, inds[n]-1 , "   " , bins[inds[n]-1], "<=", var1d[n], "<", bins[inds[n]]
## # not finished       for iq in range (1,nquantiles+1):
## # not finished          # print iq, nquantiles
## # not finished          varTarget[iq-1][inds[n]-1] = varTarget[iq-1][inds[n]-1] + y1d[iq-1][n]
## # not finished       nvar[iq-1][inds[n]-1] += 1
## # not finished       print nvar
## # not finished    for iq in range (1,nquantiles+1):            
## # not finished       meanTargetvar .append( varTarget[iq-1] / nvar[iq-1])
## # not finished    return varbins, meanTargetvar
## # not finished 
## # not finished # R9 vs pt
## # not finished nbins = 10
## # not finished meanR9pt       = [] 
## # not finished ptbins         = []
## # not finished ptbins, meanR9pt = projection(nbins, ptmin, ptmax, ptbins, pt, y, nquantiles)
## # not finished 
## # not finished # # R9 vs eta
## # not finished # nbins = 10
## # not finished # meanR9eta_lower = []
## # not finished # meanR9eta       = [] 
## # not finished # meanR9eta_upper = [] 
## # not finished # etabins         = []
## # not finished # etabins, meanR9eta = projection(nbins, etamin, etamax, etabins, eta, y, nquantiles)
## # not finished # 
## # not finished # # R9 vs phi
## # not finished # nbins = 10
## # not finished # meanR9phi_lower = []
## # not finished # meanR9phi       = [] 
## # not finished # meanR9phi_upper = [] 
## # not finished # phibins         = []
## # not finished # phibins, meanR9phi = projection(nbins, phimin, phimax, phibins, phi, y, nquantiles)
## # not finished # 
## # not finished 
## # not finished # 2D plots
## # not finished #
## # not finished 
## # not finished def plot2D(xlabel, ylabel, Xvar, Yvar, varbins, meanYvar, nquantiles, outfile ):
## # not finished     fig  = plt.figure()
## # not finished     fig.patch.set_facecolor('white')
## # not finished     axes = plt.axes()
## # not finished     plt.grid()
## # not finished     plt.xlabel(xlabel)
## # not finished     plt.ylabel(ylabel)
## # not finished     axes.xaxis.set_label_coords(1., -0.07)
## # not finished     axes.yaxis.set_label_coords(-0.07, 1.)
## # not finished     plt.scatter(Xvar   ,Yvar           , color='grey'  , marker='+')
## # not finished     for iq in range(1,nquantiles+1):
## # not finished        plt.scatter(varbins,meanYvar[iq-1]       , color='green' , marker="H", linestyle='-')
## # not finished        plt.plot   (varbins,meanYvar[iq-1]       , linestyle='-')
## # not finished     plt.plot()
## # not finished     # plt.show()
## # not finished     fig.savefig(outfile)
## # not finished 
## # not finished plot2D('$pt$'  , 'R$_9$', pt  ,R9, ptbins,  meanR9pt, nquantiles, "./meanR9_pt.pdf")
## # not finished #plot2D('$\eta$', 'R$_9$', eta ,R9, etabins, meanR9eta, nquantiles, "./meanR9_eta.pdf")
## # not finished #plot2D('$\phi$', 'R$_9$', phi ,R9, phibins, meanR9phi, nquantiles, "./meanR9_phi.pdf")
