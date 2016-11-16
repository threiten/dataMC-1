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

   def __init__(self, dmc):
      self.dataMC   = dmc

   inputDir = ""

   treeDir  = ""

   trees    = ""

   fname               = "output.root"

   evtBranches         = ["rho", "nvtx"]

   trgBranches         = [ "leadHLT_Ele27_WPTight_Gsf_vMatch", "subleadHLT_Ele27_WPTight_Gsf_vMatch" ]
   
   recoLeadBranches    = ["leadPt", "leadScEta", "leadPhi",
                          "leadR9", "leadS4", "leadSigmaIeIe", "leadEtaWidth", "leadPhiWidth", "leadCovarianceIphiIphi", "leadSigmaRR"]
   
   recoSubleadBranches = ["subleadPt", "subleadScEta", "subleadPhi",
                          "subleadR9", "subleadS4", "subleadSigmaIeIe", "subleadEtaWidth", "subleadPhiWidth", "subleadCovarianceIphiIphi", "subleadSigmaRR"]
   
   data_recoBranches   = evtBranches  + trgBranches + recoLeadBranches + recoSubleadBranches
   mc_recoBranches     = evtBranches  + recoLeadBranches + recoSubleadBranches

   df = 0

   y_corr = 0

   mcclf = []

   dataclf = []

   ptmin  =  20.
   ptmax  =  60.
   etamin = -2.5
   etamax =  2.5
   phimin = -3.14
   phimax =  3.14








   # load the dataframe from the input files
   # 
   # --------------------------------------------------------------------------------
   #
   def loadDF(self, iDir, tDir, t, start, stop, rndm = 12345):
      
      self.inputDir = iDir
      self.treeDir  = tDir
      self.trees = t 
      
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

               
      # print df
         
      # apply basic selection
      #
      df = df.query('@self.ptmin < Pt and Pt < @self.ptmax and @self.etamin < ScEta and ScEta < @self.etamax and @self.phimin < Phi and Phi < @self.phimax')
      
      print mycolors.green+"Apply basic selection"+mycolors.default
      print " ptmin  = ", self.ptmin ,"\n ptmax  = ", self.ptmax , " \n etamin = ", self.etamin, " \n etamax = ", self.etamax, " \n phimin = ", self.phimin, " \n phimax = ", self.phimax

      # print df  

      # reshuffle events
      #
      print mycolors.green+"Reshuffle events"+mycolors.default, "rndm seed  = ", rndm
      rndseed = rndm
      np.random.seed(rndseed)
      #      df['random_index'] = np.random.permutation(range(df.index.size))
      #      df.sort_values(by='random_index',inplace=True)
      #      df.set_index('random_index',inplace=True)
      #      print df
      index = list(df.index)
      np.random.shuffle(index)


      # Select a subset of events
      if   start == -1 :
         print "Invalid start-evt = -1 "
         return
      if stop  == -1 :
         stop = len(df.index)

      print mycolors.green+"Selecting events",mycolors.default, " [", start, ", ", stop, "]  out of ", len(df.index)
      index = index[start:stop]
      # df = df[start:stop]

      df = df.ix[index]
      df.reset_index(drop=True, inplace=True)

      # print df

      self.df = df

      print "DataFrame size = ", len(self.df.index)








      
   # geat the array of y
   # 
   # --------------------------------------------------------------------------------
   #
   def getNEntries(self):
      return self.df.count()
      










   # geat the array of y
   # 
   # --------------------------------------------------------------------------------
   #
   def getY(self, y):
      print 'Y = ', y
      return self.df[y]











   # run the trainings
   # 
   # --------------------------------------------------------------------------------
   #   
   def trainQuantile(self, var, alpha, pathWeights):

      # quantile regressions features
      X     = self.df.loc[:,['Pt', 'ScEta', 'Phi', 'rho']]
      # target
      Y     = self.df[var]
      
      # train quantile regression
      #
      print mycolors.green+"Train q = "+mycolors.default, alpha
      clf = GradientBoostingRegressor(loss='quantile', alpha=alpha,
                                      n_estimators=250, max_depth=3,
                                      learning_rate=.1, min_samples_leaf=9,
                                      min_samples_split=9)
      t0 = time.time()
      clf.fit(X, Y)
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
         outputName = pathWeights+"/data_weights_" + Y + '_' + str(alpha) + ".pkl"
      else :
         outputName = pathWeights+"/mc_weights_" + Y + '_' + str(alpha) + ".pkl"         
      pickle.dump(clf, gzip.open(outputName, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

















   # get the data regression weights
   # e.g filename = "./weights/data_weights" the quantile and .pkl are added here
   #
   # --------------------------------------------------------------------------------
   #
   def loadDataWeights(self, filename, quantiles):

      dbg = False
                  
      if dbg : print "Read weights for data"

      dataclf  = []
      for q in quantiles:
         dataWeights   = filename + "_" + str(q) + ".pkl"
         self.dataclf  .append(pickle.load(gzip.open(dataWeights)))
      if dbg : print "DATA weights : ", dataclf


      
   # get the MC regression weights
   # e.g filename = "./weights/mc_weights" the quanitle and .pkl are added here
   #
   # --------------------------------------------------------------------------------
   #
   def loadMcWeights(self, filename, quantiles):

      dbg = False
                  
      if dbg : print "Read weights for MC"

      mcclf    = []
      for q in quantiles:
         mcWeights = filename + "_" + str(q) + ".pkl"         
         self.mcclf    .append(pickle.load(gzip.open(mcWeights)))
      if dbg : print "MC   weights : ", mcclf











      
   # get the names to be used for X and y and fill the corrected vector
   # 
   # --------------------------------------------------------------------------------
   #
   def correctY(self, x, y, quantiles):
      
      dbg = False
      
      print "Get corrections for ", y, " with quantiles ", quantiles

      y_tmp = []
      
      # quantile regressions features
      X    = self.df.loc[:,x]
      # target e.g. y = "R9"
      Y    = self.df[y]
      print "Features: X = ", x, " target y = ", y

      
      if dbg : print "Predict MC and DATA for all quantiles"
      y_mc   = [] # list storing the n=q predictions on mc   for each event
      y_data = [] # list storing the n=q predictions on data for each event
      for q in range(0,len(quantiles)):
         y_mc  .append(self.mcclf[q]  .predict(X))
         y_data.append(self.dataclf[q].predict(X))
      if dbg : print " Initial value: Y = ", Y
      if dbg : print "                X = ", X
      if dbg : print " MC-regression = "  ,  y_mc
      if dbg : print " DATA-regression = ",  y_data

      # loop over the events #  <<-- this sucks... eventually should be vectorized, but I don't know how
      for ievt in range(0,len(Y)): 

         if dbg : print "#evt = ", ievt
         
         # brute force loop over quantiles predictions for MC
         # I would need anyway a loop to copy them in an array to use smarter search
         qmc_low  = 0
         qmc_high = 0
         q = 0         
         while q < len(quantiles): # while + if, to avoid bumping the range
            if dbg : print "mc   ", q, len(quantiles), y_mc[q][ievt],  Y[ievt] 
            if y_mc[q][ievt] < Y[ievt]:
               q+=1
            else:
               break
         if q == 0:
            qmc_low  = 0
            qmc_high = y_mc[0][ievt]
         if q < len(quantiles):
            qmc_low  = y_mc[q-1][ievt]
            qmc_high = y_mc[q ][ievt]
         else:
            qmc_low  = y_mc[q-1][ievt]
            qmc_high = 1            
         if dbg : print "mc-quantile    ", q, " --> [ ", qmc_low ,qmc_high, " ]"

         #
         # q is the one we find on mc
         if dbg:
            qt = 0 
            while qt < len(quantiles):
               print "data ", qt, len(quantiles), y_data[qt][ievt],  Y[ievt]
               qt+=1
         qdata_low  = 0
         qdata_high = 0
         if q == 0:
            qdata_low  = 0
            qdata_high = y_data[0][ievt]
         if q < len(quantiles):
            qdata_low  = y_data[q-1][ievt]
            qdata_high = y_data[q ][ievt]
         else:
            qdata_low  = y_data[q-1][ievt]
            qdata_high = 1
         if dbg : print "data-quantiles  --> [ ", qdata_low ,qdata_high, " ]"



         # interplopate the correction
         y_corr = (qdata_high-qdata_low)/(qmc_high-qmc_low) * (Y[ievt] - qmc_low) + qdata_low
         if dbg : print "Apply correction: Input value = ", Y[ievt], " --> corrected value = ", y_corr

         y_tmp.append(y_corr)
         
      #self.y_corr = y_tmp
      ycorr = y+"_corr"
      self.df[ycorr] = y_tmp
      # print self.df[ycorr]









     
   # get the corrected array of y
   # 
   # --------------------------------------------------------------------------------
   #
   def getCorrectedY(self, y):
      ycorr = y+"_corr"
      print 'Corrected Y = ', ycorr
      return self.df[ycorr]
      # return self.y_corr












   # Scatter plots to check the quantiles
   # 
   # --------------------------------------------------------------------------------
   #
   def plotQuantiles(self, quantiles, xVar, nbins, xMin, xMax, yVar, xLabel , yLabel, outfile ): # << R9 hardcoded !!! generalize for all vars

      xx = self.df.ix[:,[xVar]]

      # quantile regressions features
      X     = self.df.loc[:,['Pt', 'ScEta', 'Phi', 'rho']]
      # target
      yy    = self.df[yVar]

      # compute the predictions
      y = []
      y_1d = []
      for iq in range(0,len(quantiles)):
         # print 'Predict ', quantiles[iq], ' quantile'
         pred = 0
         if (self.dataMC == "mc"):
            pred = self.mcclf[iq].predict(X)
         if (self.dataMC == "data"):
            pred = self.dataclf[iq].predict(X)
         y.append(pred)
         y_1d.append(pred.ravel())


      # projection y vs x for the scatter plot
      def projection(nbins, min, max, var, y_1d):
         varbins = []
         bins  = np.linspace(min, max,nbins+1)
         for i in range(0,nbins):
            varbins.append(0.5*(bins[i]+bins[i+1]))
         varyy = []
         for iq in range(0,len(quantiles)):
            y_pred      = y_1d[iq]
            inds        = np.digitize(var, bins)    
            vyy         = np.zeros(nbins)
            nvar        = np.zeros(nbins)
            var1d       = np.ravel(var)
            for n in range(var.size):
               # print inds[n]-1 , "   " , bins[inds[n]-1], "<=", var1d[n], "<", bins[inds[n]]
               vyy[inds[n]-1] = vyy[inds[n]-1] + y_pred[n]
               nvar[inds[n]-1] += 1
            meanyyvar = vyy / nvar # the use of the mean is not great...
            varyy.append(meanyyvar)
         return varbins, varyy
     
      # yy vs xx
      meanyyxx = []
      xxbins   = []
      xxbins, meanyyxx = projection(nbins, xMin, xMax, xx, y_1d)
      
      def plot2D(xlabel, ylabel, Xvar, Yvar, varbinning, meanY, outfile ):
         fig  = plt.figure()
         # fig.patch.set_facecolor('white')
         axes = plt.axes()
         plt.grid()
         plt.xlabel(xlabel)
         plt.ylabel(ylabel)
         axes.xaxis.set_label_coords(1., -0.07)
         axes.yaxis.set_label_coords(-0.07, 1.)
         plt.scatter(Xvar   ,Yvar           , color='grey'  , marker='+')
         for iq in range(0,len(quantiles)):
            if (quantiles[iq] == 0.5):
               plt.plot   (varbinning, meanY[iq] , linestyle='-', color='red', marker="H")
            else:
               plt.plot   (varbinning, meanY[iq] , linestyle='-', color='blue', marker="H")
            plt.plot()
         fig.savefig(outfile)
         return plt

      
      plot = plot2D(xLabel, yLabel, xx  ,yy, xxbins,  meanyyxx, outfile)
      
      return plot
