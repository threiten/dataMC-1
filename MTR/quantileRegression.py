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
import bisect
# from joblib import Parallel, delayed
from sklearn.externals.joblib import Parallel, parallel_backend, register_parallel_backend
from joblib import delayed
import os
import ROOT as rt
import copy as cp

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




# -------------------------------------------------------------------------------------
def setupJoblib(ipp_profile='default'):
    
    import ipyparallel as ipp
    from ipyparallel.joblib import IPythonParallelBackend
    global joblib_rc,joblib_view,joblib_be
    joblib_rc = ipp.Client(profile=ipp_profile)
    joblib_view = joblib_rc.load_balanced_view()
    joblib_be = IPythonParallelBackend(view=joblib_view)
    
    register_parallel_backend('ipyparallel',lambda : joblib_be,make_default=True)


#
# helper class to peform quantile-based corrections
#
# --------------------------------------------------------------------------------
class Corrector:
   
   # store regressors
   def __init__(self,mcclf,dataclf,X,Y):
      self.mcqtls   = np.array([clf.predict(X) for clf in mcclf])
      self.dataqtls = np.array([clf.predict(X) for clf in dataclf])

      self.Y = Y
      
   # correction is actually done here
   def correctEvent(self,iev):

      mcqtls = self.mcqtls[:,iev]
      dataqtls = self.dataqtls[:,iev]
      Y = self.Y[iev]
      
      qmc   = bisect.bisect_right(mcqtls,Y)
      if qmc == 0:
         qmc_low,qdata_low   = 0,0                              # all shower shapes have a lower bound at 0
         qmc_high,qdata_high = mcqtls[qmc],dataqtls[qmc]
      elif qmc < len(mcqtls):
         qmc_low,qdata_low   = mcqtls[qmc-1],dataqtls[qmc-1]
         qmc_high,qdata_high = mcqtls[qmc],dataqtls[qmc]
      else:
         qmc_low,qdata_low   = mcqtls[qmc-1],mcqtls[qmc-1]
         qmc_high,qdata_high = mcqtls[-1]*1.2,dataqtls[-1]*1.2   # some variables (e.g. sigmaRR) have values above 1
         # to set the value for the highest quantile 20% higher
                                                                       
      return (qdata_high-qdata_low)/(qmc_high-qmc_low) * (Y - qmc_low) + qdata_low

   def __call__(self):
      return np.array([ self.correctEvent(iev) for iev in xrange(self.Y.size) ]).ravel()

def applyCorrection(mcclf,dataclf,X,Y):
   return Corrector(mcclf,dataclf,X,Y)()


#
# helper class to evaluate ID MVA
# 
# --------------------------------------------------------------------------------
#
class IdMvaComputer:

   def __init__(self,wd,weightsEB,weightsEE,correct=[]):
      rt.gROOT.LoadMacro(os.path.join(wd,"phoIDMVAonthefly.C"))
      
      self.rhoSubtraction = False
      if type(correct) == dict:
         self.rhoSubtraction = correct["rhoSubtraction"]
         correct = correct["correct"]

      self.X = rt.phoIDInput()
      self.readerEB = rt.bookReadersEB(weightsEB, self.X)
      self.readerEE = rt.bookReadersEE(weightsEE, self.X, self.rhoSubtraction)

      # print ("IdMvaComputer.__init__")
      
      columns = ["ScEnergy","ScEta","rho","R9","SigmaIeIe","PhiWidth","EtaWidth","CovarianceIetaIphi","S4","PhoIso03","ChIso03","ChIso03worst","SigmaRR","ScPreshowerEnergy","Pt"]


      if self.rhoSubtraction:
         self.effareas = np.array([[0.0000, 0.1210],   
                                   [1.0000, 0.1107],
                                   [1.4790, 0.0699],
                                   [2.0000, 0.1056],
                                   [2.2000, 0.1457],
                                   [2.3000, 0.1719],
                                   [2.4000, 0.1998],
         ])
         
      # make list of input columns
      self.columns = map(lambda x: x+"_corr" if x in correct else x, columns)
      
   def __call__(self,X):

      # make sure of order of the input columns and convert to a numpy array
      Xvals = X[ self.columns ].values

      return np.apply_along_axis( self.predict, 1, Xvals ).ravel()
      
   def predict(self,row):
      return self.predictEB(row) if np.abs(row[1]) < 1.5 else self.predictEE(row)
      # return self.predictEB(row)
      
   def predictEB(self,row):
      # use numeric indexes to speed up
      # print ("IdMvaComputer.predictEB")
      self.X.phoIdMva_SCRawE_          = row[0]
      self.X.phoIdMva_ScEta_           = row[1]
      self.X.phoIdMva_rho_             = row[2]
      self.X.phoIdMva_R9_              = row[3]
      self.X.phoIdMva_covIEtaIEta_     = row[4] # this is really sieie
      self.X.phoIdMva_PhiWidth_        = row[5]
      self.X.phoIdMva_EtaWidth_        = row[6]
      self.X.phoIdMva_covIEtaIPhi_     = row[7]
      self.X.phoIdMva_S4_              = row[8]
      self.X.phoIdMva_pfPhoIso03_      = row[9]
      self.X.phoIdMva_pfChgIso03_      = row[10]
      self.X.phoIdMva_pfChgIso03worst_ = row[11]
      return self.readerEB.EvaluateMVA("BDT")

   def effArea(self,eta):
      ibin = min(self.effareas.shape[0]-1,bisect.bisect_left(self.effareas[:,0],eta))
      return self.effareas[ibin,1]
   
   def predictEE(self,row):
      # print ("IdMvaComputer.predictEB")
      self.X.phoIdMva_SCRawE_          = row[0]
      self.X.phoIdMva_ScEta_           = row[1]
      self.X.phoIdMva_rho_             = row[2]
      self.X.phoIdMva_R9_              = row[3]
      self.X.phoIdMva_covIEtaIEta_     = row[4] # this is really sieie
      self.X.phoIdMva_PhiWidth_        = row[5]
      self.X.phoIdMva_EtaWidth_        = row[6]
      self.X.phoIdMva_covIEtaIPhi_     = row[7]
      self.X.phoIdMva_S4_              = row[8]
      self.X.phoIdMva_pfPhoIso03_      = row[9]
      if self.rhoSubtraction: self.X.phoIdMva_pfPhoIso03_ = max(2.5, self.X.phoIdMva_pfPhoIso03_ - self.effArea( np.abs(self.X.phoIdMva_ScEta_))*self.X.phoIdMva_rho_ - 0.0034*row[14] )
      self.X.phoIdMva_pfChgIso03_      = row[10]
      self.X.phoIdMva_pfChgIso03worst_ = row[11]
      self.X.phoIdMva_ESEffSigmaRR_    = row[12]
      esEn                             = row[13]
      ScEn                             = row[0]
      self.X.phoIdMva_esEnovSCRawEn_ = esEn/ScEn
      return self.readerEE.EvaluateMVA("BDT")


def computeIdMva(wd,weightsEB,weightsEE,correct,X):
   return IdMvaComputer(wd,weightsEB,weightsEE,correct)(X)

#
# helper class to peform stochastic isolation corrections
# 
# --------------------------------------------------------------------------------
class IsolationCorrector:

   def __init__(self,wd,corr_file):

      if os.path.exists(os.path.join(wd,"../phoIsoStoch/IsolationCorrection_C.so")):
         rt.gSystem.Load(os.path.join(wd,"../phoIsoStoch/IsolationCorrection_C.so"))
      else:
         rt.gROOT.LoadMacro(os.path.join(wd,"../phoIsoStoch/IsolationCorrection.C"))
      self.isoCorr = rt.IsolationCorrection(corr_file)

   def __call__(self,X):
      
      Xvals = X[ ['rho','ScEta','PhoIso03'] ].values

      return np.apply_along_axis( self.predict, 1, Xvals ).ravel()
   
   def predict(self,row):
      rho,eta,iso = row[0],np.abs(row[1]),row[2]
      return iso+self.isoCorr.getExtra(eta,rho)

def applyIsoCorrection(wd,corr_file,X):
   return IsolationCorrector(wd,corr_file)(X)


#
# 
# --------------------------------------------------------------------------------
#
class quantileRegression:

   def __init__(self, label):
      self.label   = label
      self.dataMC  = label.split("_",1)[0]

      self.inputDir = ""

      self.treeDir  = ""

      self.trees    = ""

      self.fname               = "output.root"

      self.evtBranches         = ["run", "rho", "nvtx", "mass", "weight","puweight"]
      
      self.trgBranches         = [ "leadHLT_Ele27_WPTight_Gsf_vMatch", "subleadHLT_Ele27_WPTight_Gsf_vMatch" ]
      
      self.eleMatchBranches    = [ "leadEleMatch", "subleadEleMatch" ]
      
      self.recoLeadBranches    = ["leadPt", "leadScEta", "leadPhi",
                                  'leadScEnergy', 'leadScPreshowerEnergy', "leadSigmaRR" ,
                                  'leadPhoIso03', 'leadChIso03', 'leadChIso03worst', 
                                  'leadPhoIDMVA' ] # , 'leadSigEOverE','leadRecoSigEOverE','leadUnsmearedSigmaEoE','leadAfterSSTrSigEOverE']
      
      self.recoSubleadBranches    = ["subleadPt", "subleadScEta", "subleadPhi",
                                  'subleadScEnergy', 'subleadScPreshowerEnergy', "subleadSigmaRR" ,
                                  'subleadPhoIso03', 'subleadChIso03', 'subleadChIso03worst', 
                                  'subLeadPhoIDMVA'] # , 'subleadSigEOverE','subleadRecoSigEOverE','subleadUnsmearedSigmaEoE','subleadAfterSSTrSigEOverE']

      self.recoLeadSSBranches  = ["leadR9", "leadS4", "leadEtaWidth", "leadPhiWidth", 
                                  "leadSigmaIeIe", 'leadCovarianceIetaIphi', "leadCovarianceIphiIphi" ]

      self.recoSubleadSSBranches  = ["subleadR9", "subleadS4", "subleadEtaWidth", "subleadPhiWidth", 
                                     "subleadSigmaIeIe", 'subleadCovarianceIetaIphi', "subleadCovarianceIphiIphi" ]
      
      # the uncorrected variables are used only if you switched on the corrections in flashgg
      # NB: there is no PhiWIdth uncorrected !
      #      self.recoLeadUncorrSSBranches  = ["leadUncorrR9", "leadUncorrS4", "leadUncorrEtaWidth", "leadPhiWidth",
      #                                        "leadUncorr_non5x5_sigmaIetaIeta", "leadCovarianceIetaIphi", "leadCovarianceIphiIphi"]
      #
      #
      # NB: there is no PhiWIdth uncorrected !
      #      self.recoSubleadUncorrSSBranches  = ["subleadUncorrR9", "subleadUncorrS4", "subleadUncorrEtaWidth", "subleadPhiWidth", 
      #                                        "subleadUncorr_non5x5_sigmaIetaIeta", "subleadCovarianceIetaIphi", "subleadCovarianceIphiIphi"]
      
      self.data_recoBranches   = self.evtBranches  + self.trgBranches + self.eleMatchBranches + self.recoLeadBranches + self.recoLeadSSBranches + self.recoSubleadBranches + self.recoSubleadSSBranches 
      self.mc_recoBranches     = self.evtBranches  +                    self.eleMatchBranches + self.recoLeadBranches + self.recoLeadSSBranches + self.recoSubleadBranches + self.recoSubleadSSBranches 

      self.df = 0
      
      self.y_corr = 0
      
      self.mcclf = []
      
      self.dataclf = []
      
      self.ptmin  =  25.
      self.ptmax  =  150.
      self.etamin = -2.5
      self.etamax =  2.5
      self.phimin = -3.14
      self.phimax =  3.14








   # load the dataframe from the input files
   # 
   # --------------------------------------------------------------------------------
   #
   def loadDF(self, iDir, tDir, t, start, stop, runStart = 0, runStop = 999999999, rndm = 12345):

      dbg = False
      
      self.inputDir = iDir
      self.treeDir  = tDir
      self.trees = t 
      
      fname    = self.inputDir+self.fname
      
      print "NT branches: "
      recoBranches = self.data_recoBranches 
      if self.dataMC == "mc" :
         recoBranches = self.mc_recoBranches
      #print recoBranches
            
      # use common names for the traning dataset
      #
      uniformColumnsNames = ["run", "rho", "nvtx" , "mass", "weight","puweight", 
                             # "SigMoM", "RecoSigMoM",
                             "Pt", "ScEta", "Phi",
                             'ScEnergy', 'ScPreshowerEnergy', "SigmaRR" ,
                             'PhoIso03', 'ChIso03', 'ChIso03worst' ,
                             'PhoIDMVA', #  'SigEOverE','RecoSigEOverE','UnsmearedSigEOverE','AfterSSTrSigEOverE',
                             "R9", "S4", "EtaWidth", "PhiWidth", 
                             "SigmaIeIe", 'CovarianceIetaIphi', "CovarianceIphiIphi"]

                                                                                
      # attach the tree structure to the tree names
      #
      trees = []
      for t in self.trees:
         trees.append(self.treeDir+t)
      self.trees = trees
      #print "trees: ",trees

      # read the trees into a dataframe
      #
      print "Adding trees into a DataFrame"
      i = 0
      for t in self.trees:
         print "  adding ", t
         if i == 0:            
            df = rpd.read_root(fname,t,columns=recoBranches)
         else:
            df = pd.concat([df, rpd.read_root(fname,t,columns=recoBranches)])
         if dbg : print df.count()
         i+=1
      
      print "number of events:", len(df.index)

      # select run period on data

      if (self.dataMC == "data" ):
         df = df.query('@runStart <= run and run <= @runStop')

      df = df.reset_index() # otherwise it will keep the index of each of the added dataframe
      # print df 

#	      # sigmaM/M
#	      df["SigMoM"]     = 0.5*np.sqrt(np.power(df.leadSigEOverE,2) + np.power(df.subleadSigEOverE,2))
#	      df["RecoSigMoM"] = 0.5*np.sqrt(np.power(df.leadRecoSigEOverE,2) + np.power(df.subleadRecoSigEOverE,2) + 
#	                                     (np.power(df.leadSigEOverE,2) - np.power(df.leadUnsmearedSigmaEoE, 2)) +
#	                                     ( np.power(df.subleadSigEOverE,2) - np.power(df.subleadUnsmearedSigmaEoE, 2)))      
#	      self.evtBranches += [ "SigMoM", "RecoSigMoM" ]

      if (self.dataMC == "data" ):
         # trigger matching logic :         
         # get the indices of the lead matching the trigger and use them to select the sublead
         # get the indices of the sublead matching the trigger and use them to select the lead
         # NB: I can use the event twice if both lead and sublead trigger

         print "Count df"
         if dbg : print df.count()
         #
         # Select only lead / sublead candidates with eleId=isHLTsafe
         if dbg : print df["leadEleMatch"].sum()
         if dbg : print df["subleadEleMatch"].sum()
         # the matching decision is stored as a float -> convert it to bool
         df.leadEleMatch   = df.leadEleMatch.astype(bool)
         df.subleadEleMatch = df.subleadEleMatch.astype(bool)
         idx_eleMatchLead    = df[df["leadEleMatch"]].index.tolist()
         idx_eleMatchSublead = df[df["subleadEleMatch"]].index.tolist()
         print "# lead eleMatch    = ", len(idx_eleMatchLead)
         print "# sublead eleMatch = ", len(idx_eleMatchSublead)
         
         #
         # print df["leadHLT_Ele27_WPTight_Gsf_vMatch"].sum()
         # print df["subleadHLT_Ele27_WPTight_Gsf_vMatch"].sum()
         # the trigger decision is stored as a float -> convert it to bool
         df.leadHLT_Ele27_WPTight_Gsf_vMatch    = df.leadHLT_Ele27_WPTight_Gsf_vMatch.astype(bool)
         df.subleadHLT_Ele27_WPTight_Gsf_vMatch = df.subleadHLT_Ele27_WPTight_Gsf_vMatch.astype(bool)
         idx_trgLead    = df[df["leadHLT_Ele27_WPTight_Gsf_vMatch"]].index.tolist()
         idx_trgSublead = df[df["subleadHLT_Ele27_WPTight_Gsf_vMatch"]].index.tolist()
         print "# lead trig   = ", len(idx_trgLead)
         print "# sublead trg = ", len(idx_trgSublead)

         #
         # select trigger AND ele safe match for lead / sublead
         idx_trgEleMatchLead    = np.intersect1d(idx_trgLead,    idx_eleMatchLead)
         idx_trgEleMatchSublead = np.intersect1d(idx_trgSublead, idx_eleMatchSublead) 
         print "# lead trig and eleMatch    = ", len(idx_trgEleMatchLead)
         print "# sublead trig and eleMatch = ", len(idx_trgEleMatchSublead)
         
         #print "evtbranches", self.evtBranches
         #print  "unif", uniformColumnsNames
         # if the lead triggered use the sublead
         #
         print "Data Sublead"
         dataSublead = df[ self.evtBranches + self.recoSubleadBranches + self.recoSubleadSSBranches]
         dataSublead = dataSublead.loc[idx_trgEleMatchLead]
         print dataSublead.columns         
         dataSublead.columns = uniformColumnsNames
         print "lead size: ", len(dataSublead.index)
         #
         # if the sublead triggered use the lead
         #
         print "Data Lead"
         dataLead    = df[ self.evtBranches + self.recoLeadBranches + self.recoLeadSSBranches]
         dataLead    = dataLead.loc[idx_trgEleMatchSublead]
         dataLead.columns = uniformColumnsNames
         print "sublead size: ", len(dataLead.index)
         
         #
         # concatenate leading and subleadind
         frames = [dataLead, dataSublead]
         data = pd.concat(frames)
         if dbg :
            print "Data count"
            print data.count()
         # reset the rows indexing
         df = data.reset_index()
         if dbg : print df


      if (self.dataMC == "mc" ):
         # no trigger matching on mc but still require eleId=isHLTsafe
         # print df["leadEleMatch"].sum()
         # print df["subleadEleMatch"].sum()
         # the matching decision is stored as a float -> convert it to bool
         df.leadEleMatch   = df.leadEleMatch.astype(bool)
         df.subleadEleMatch = df.subleadEleMatch.astype(bool)
         idx_eleMatchLead    = df[df["leadEleMatch"]].index.tolist()
         idx_eleMatchSublead = df[df["subleadEleMatch"]].index.tolist()
         print "# lead eleMatch    = ", len(idx_eleMatchLead)
         print "# sublead eleMatch = ", len(idx_eleMatchSublead)

         print "MC Sublead"
         dataSublead = df[ self.evtBranches + self.recoSubleadBranches + self.recoSubleadSSBranches]
         dataSublead = dataSublead.loc[idx_eleMatchLead]
         dataSublead.columns = uniformColumnsNames
         print "lead size: ",len(dataSublead.index)
         #
         print "MC Lead"
         dataLead    = df[ self.evtBranches + self.recoLeadBranches + self.recoLeadSSBranches]
         dataLead    = dataLead.loc[idx_eleMatchSublead]
         dataLead.columns = uniformColumnsNames
         print "sublead size: ",len(dataLead.index)
         # 

         # concatenate leading and subleadind
         frames = [dataLead, dataSublead]
         data = pd.concat(frames)
         # reset the rows indexing
         df = data.reset_index()
         if dbg : print data.count()


      # smear the Reco sigmaE/E
      # 
#       df["RecoSigEOverEsmear"] = np.sqrt(np.power(df.RecoSigEOverE,2) + (np.power(df.SigEOverE,2) - np.power(df.UnsmearedSigEOverE, 2)))
      

      print "Count final dataset"
      print len(df.index)
      #print df.columns
            



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
      ## if stop  == -1 :
      ##    stop = len(df.index)

      print mycolors.green+"Selecting events",mycolors.default, " [", start, ", ", stop, "]  out of ", len(df.index)
      index = index[start:stop]
      # df = df[start:stop]

      df = df.ix[index]
      df.reset_index(drop=True, inplace=True)

      # print df
      self.df = df

      print "DataFrame size = ", len(self.df.index)

      # save the DF locally for future use
      dfname =  'df_' + self.label + '_' + str(start) + '-' + str(stop) + '.h5'
      hdf = pd.HDFStore(dfname)
      hdf.put('df', self.df)
      hdf.close()
         







      

   # load the dataframe from an already existing h5 file
   # 
   # --------------------------------------------------------------------------------
   #
   def loadDFh5(self, h5name, start, stop):

      df = 0
      import os.path
      if os.path.exists(h5name):
         print 'Loading dataframe from : ', h5name
         df = pd.read_hdf(h5name, 'df')
      else:
         print "The h5 file ", h5name, " does not exist"
         return

      # Select a subset of events
      if   start == -1 :
         print "Invalid start-evt = -1 "
         return
      if stop  == -1 :
         stop = len(df.index)

      print mycolors.green+"Selecting events",mycolors.default, " [", start, ", ", stop, "]  out of ", len(df.index)
      # index = index[start:stop]
      df = df[start:stop]

      # df = df.ix[index]
      df.reset_index(drop=True)

      self.df = df
      print "number of events:", len(df.index)






         
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











   # returns a DF satisfying the query (the self.df remains untouched! )
   # 
   # --------------------------------------------------------------------------------
   #
   def applyCutsToDF(self, var, min, max, inout):


      df = self.df       
      
      querystr = ''
      if inout == 'inside':
         querystr = '@min < {} and {} < @max'.format(var,var)
         print min, ' < ', var, ' and ', var, ' < ', max 
      elif inout == 'outside':
         querystr = '{}<@min or @max <{}'.format(var,var)
         print var, ' < ', min, ' or ', max, ' < ', var 

      df = df.query(querystr)
      #print df

      # reset the rows indexing
      df = df.reset_index(drop=True)

      self.df = df   
   


   # run the trainings
   # 
   # --------------------------------------------------------------------------------
   #   
   def trainQuantile(self, var, alpha, pathWeights, EBEE ="", maxDepth = 3, minLeaf = 9, useWeights = False):
      #self.df = self.df.query('weight>0')
      if   EBEE == 'EB':
         self.applyCutsToDF('ScEta', -1.4442, 1.4442, 'inside')
      elif EBEE == 'EE':
         self.applyCutsToDF('ScEta', -1.57, 1.57, 'outside')
      elif EBEE == '':
         print "Training both EB and EE together"

      # quantile regressions features
      X     = self.df.loc[:,['Pt', 'ScEta', 'Phi', 'rho']]
      # target
      Y     = self.df[var]
      #event weight
      w     = self.df['weight']
      #   The isolation variables have a discrete single value at zero and then a smooth distribution
      #   To avoid degeneracies in the quantile matching I artificially smear them subracting rho
      #   (like a PU correction)
      if var == "PhoIso03":
         self.df['PhoIso03rho'] = self.df['PhoIso03'] - 0.1*self.df['rho']
         Y = self.df['PhoIso03rho']
         var = var +"rho"
      if var == "ChIso03":
         self.df['ChIso03rho'] = self.df['ChIso03'] - 0.1*self.df['rho']
         Y = self.df['ChIso03rho']
         var = var +"rho"
      if var == "ChIso03worst":
         self.df['ChIso03rhoworst'] = self.df['ChIso03worst'] - 0.1*self.df['rho']
         Y = self.df['ChIso03rhoworst']
         var = var +"rho"

      print "Data runMin= ", self.df['run'].min()," - runMax = ", self.df['run'].max()

      # train quantile regression
      #
      print mycolors.green+"Train q = "+mycolors.default, alpha
      clf = GradientBoostingRegressor(loss='quantile', alpha=alpha,
                                      n_estimators=250, max_depth=maxDepth,
                                      learning_rate=.1, min_samples_leaf=minLeaf,
                                      min_samples_split=minLeaf)
      t0 = time.time()
      if (useWeights) :
       #  w=self.df.copy()
       #  w.loc[w.query("weight<0").index,"weight"]=0
       #  w=w['weight']
       #  w.abs() for absolute values
         clf.fit(X, Y, w)
      else:
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
         if   EBEE != "":
            outputName = pathWeights+"/data_weights_" + EBEE + "_" + var + "_" + str(alpha) + ".pkl"
         else:
            outputName = pathWeights+"/data_weights_" + var + "_" + str(alpha) + ".pkl"
      else :
         if   EBEE != "":
            outputName = pathWeights+"/mc_weights_" + EBEE + "_" + var + "_" + str(alpha) + ".pkl"
         else:
            outputName = pathWeights+"/mc_weights_" + var + "_" + str(alpha) + ".pkl"         

      print outputName
      #work here for the protocol and understand to lower the  check      
      pickle.dump(clf, gzip.open(outputName, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)









   # traing a final regression on the mc Y_corr variables.
   # To use it, the user will give (X,Y) and get Y_corr.
   # 
   # --------------------------------------------------------------------------------
   #   
   def trainOnCorrections(self, x, ylist, quantiles, outputPath, maxDepth = 3, minLeaf = 9):
      
      print "Correct all variables", ylist
      # This loads the weights for data and mc, compute the correction and apply it to the variable
      self.correctAllY(x, ylist, quantiles )      
      # print self.df
      
      x4vars = x
      for y in ylist:         
         
         x = x4vars + [y]
         X     = self.df.loc[:,x]
         ycorr = y+"_corr"                                                          
         # target as the difference between the corrected and the non corrected
         Y     = self.df[ycorr]-self.df[y]
         print "Training on ", x, " for ", ycorr, " - " , y
         
         # train regression on the corrected variables
         #
         print mycolors.green+"Training the final regression for "+mycolors.default, ycorr
         clf = GradientBoostingRegressor(loss='ls',
                                         n_estimators=250, max_depth=maxDepth,
                                         learning_rate=.1, min_samples_leaf=minLeaf,
                                         min_samples_split=minLeaf)
         
         t0 = time.time()
         clf.fit(X, Y)
         t1 = time.time()
         print " time = ", t1-t0
         
         # self.mcclf  =clf  # this is for debugging only
        
         print "Save weights"
         outputName = outputPath+"/weights_corrections_" + y + ".pkl"
         pickle.dump(clf, gzip.open(outputName, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
         






         

   # Apply the corrections regression to Y
   #
   # --------------------------------------------------------------------------------
   #
   def predictY(self, x, y, filename):

      dbg = True

      # quantile regressions features
      X     = self.df.loc[:,x]
      # target
      Y     = self.df.loc[:,y]


      # get the correction regression weights
      # e.g filename = "./weightsCorrections/weights_corrections" the quanitle and .pkl are added here
      if dbg : print "Read weights for MC for variable ", y
      weights = filename + "_" + y + ".pkl"
      if dbg : print "Correction file : ", weights
      predclf = pickle.load(gzip.open(weights))
      if dbg : print "Correction weights for ", y, " : ", predclf

      return predclf.predict(X)
      
   # get the data regression weights
   # e.g filename = "./weights/data_weights" the quantile and .pkl are added here
   #
   # --------------------------------------------------------------------------------
   #
   def loadDataWeights(self, filename, var, quantiles):

      dbg = False
                  
      if dbg : print "Read weights for data for variable ", var

      self.dataclf = []
      for q in quantiles:
         dataWeights   = filename + "_" + var + '_' + str(q) + ".pkl"
         self.dataclf  .append(pickle.load(gzip.open(dataWeights)))
      if dbg : print "DATA weights : ", self.dataclf
      
   # get the MC regression weights
   # e.g filename = "./weights/mc_weights" the quanitle and .pkl are added here
   #
   # --------------------------------------------------------------------------------
   #
   def loadMcWeights(self, filename, var, quantiles):

      dbg = False
                  
      if dbg : print "Read weights for MC for variable ", var

      self.mcclf = []
      for q in quantiles:
         mcWeights = filename + "_" + var + '_' + str(q) + ".pkl"         
         self.mcclf    .append(pickle.load(gzip.open(mcWeights)))
      if dbg : print "MC   weights : ", mcclf

   # get the names to be used for X and y and fill the corrected vector
   # 
   # --------------------------------------------------------------------------------
   #
   def correctY(self, x, y, quantiles ):
      
      dbg = False
      
      print "Get corrections for ", y, " with quantiles ", quantiles

      y_tmp = []
      
      # quantile regressions features
      X    = self.df.loc[:,x]
      # target e.g. y = "R9"
      Y    = self.df[y]
      print "Features: X = ", x, " target y = ", y
      # print X, Y

      if y == 'PhoIso03' or y == 'ChIso03' or y == 'ChIso03worst':
         Y = Y - 0.1*self.df['rho']
      
      if dbg : print "Predict MC and DATA for all quantiles"
      y_mc   = [] # list storing the n=q predictions on mc   for each event
      y_data = [] # list storing the n=q predictions on data for each event

      for q in range(0, len(quantiles)):
         y_mc  .append(self.mcclf[q]  .predict(X))
         y_data.append(self.dataclf[q].predict(X))
      if dbg : print " Initial value: Y = ", Y
      if dbg : print "                X = ", X
      if dbg : print " MC-regression = "  ,  y_mc
      if dbg : print " DATA-regression = ",  y_data

         
      # loop over the events #  <<-- this sucks... eventually should be vectorized, but I don't know how
      for ievt in range(0, len(Y)): 

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
            qmc_low  = 0                               # all shower shapes have a lower bound at 0
            qmc_high = y_mc[0][ievt]
         elif q < len(quantiles):
            qmc_low  = y_mc[q-1][ievt]
            qmc_high = y_mc[q ][ievt]
         else:
            qmc_low  = y_mc[q-1][ievt]
            qmc_high = quantiles[len(quantiles)-1]*1.2 # some variables (e.g. sigmaRR) have values above 1
                                                       # to set the value for the highest quantile 20% higher
         if dbg : print "mc-quantile    ", q, " --> [ ", qmc_low ,qmc_high, " ]"
         
         #
         # q is the one we find on mc
         if dbg:
            qt = 0 
            while qt < len(quantiles):
               print "data ", qt, len(quantiles), y_data[qt][ievt] # ,  Y[ievt]
               qt+=1
         qdata_low  = 0
         qdata_high = 0
         if q == 0:
            qdata_low  = 0                              # all shower shapes have a lower bound at 0
            qdata_high = y_data[0][ievt]
         elif q < len(quantiles):
            qdata_low  = y_data[q-1][ievt]
            qdata_high = y_data[q ][ievt]
         else:
            qdata_low  = y_data[q-1][ievt]
            qdata_high = quantiles[len(quantiles)-1]*1.2 # see comment above for mc            
         if dbg : print "data-quantiles  --> [ ", qdata_low ,qdata_high, " ]"



         # interplopate the correction
         y_corr = (qdata_high-qdata_low)/(qmc_high-qmc_low) * (Y[ievt] - qmc_low) + qdata_low
         if dbg : print "Apply correction: Input value = ", Y[ievt], " --> corrected value = ", y_corr

         y_tmp.append(y_corr)
         

      #   The isolation variables have a discrete single value at zero and then a smooth distribution
      #   To avoid degeneracies in the quantile matching I artificially smear them subracting rho
      #   (like a PU correction)
      if y == "PhoIso03" or y == "ChIso03" or y == "ChIso03worst":
         y_tmp = y_tmp + 0.1*self.df['rho']

      #self.y_corr = y_tmp
      ycorr = y+"_corr"
      self.df[ycorr] = y_tmp
      # print self.df[ycorr]
      
   def correctYTime(self, x, y, quantiles):
      dbg= False  
      print "Get corrections for ", y, " with quantiles ", quantiles

      y_tmp = []
      x2=['Pt', 'ScEta', 'Phi', 'rho',"runperiod"]
      # quantile regressions features
      X1    = self.df.loc[:,x]
      X2    = self.df.loc[:,x2]
      # target e.g. y = "R9"
      Y    = self.df[y]
      print "Features: X = ", x, " target y = ", y, "for mc"
      print "Features: X = ", x2, " target y = ", y,"for data"

      if y == 'PhoIso03' or y == 'ChIso03' or y == 'ChIso03worst':
         Y = Y - 0.1*self.df['rho']
      
      if dbg : print "Predict MC and DATA for all quantiles"
      y_mc   = [] # list storing the n=q predictions on mc   for each event
      y_data = [] # list storing the n=q predictions on data for each event

      for q in range(0, len(quantiles)):
         y_mc  .append(self.mcclf[q]  .predict(X1))
         y_data.append(self.dataclf[q].predict(X2))
      if dbg : print " Initial value: Y = ", Y
      if dbg : print "                X = ", X
      if dbg : print " MC-regression = "  ,  y_mc
      if dbg : print " DATA-regression = ",  y_data

         
      # loop over the events #  <<-- this sucks... eventually should be vectorized, but I don't know how
      for ievt in range(0, len(Y)): 

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
            qmc_low  = 0                               # all shower shapes have a lower bound at 0
            qmc_high = y_mc[0][ievt]
         elif q < len(quantiles):
            qmc_low  = y_mc[q-1][ievt]
            qmc_high = y_mc[q ][ievt]
         else:
            qmc_low  = y_mc[q-1][ievt]
            qmc_high = quantiles[len(quantiles)-1]*1.2 # some variables (e.g. sigmaRR) have values above 1
                                                       # to set the value for the highest quantile 20% higher
         if dbg : print "mc-quantile    ", q, " --> [ ", qmc_low ,qmc_high, " ]"
         
         #
         # q is the one we find on mc
         if dbg:
            qt = 0 
            while qt < len(quantiles):
               print "data ", qt, len(quantiles), y_data[qt][ievt] # ,  Y[ievt]
               qt+=1
         qdata_low  = 0
         qdata_high = 0
         if q == 0:
            qdata_low  = 0                              # all shower shapes have a lower bound at 0
            qdata_high = y_data[0][ievt]
         elif q < len(quantiles):
            qdata_low  = y_data[q-1][ievt]
            qdata_high = y_data[q ][ievt]
         else:
            qdata_low  = y_data[q-1][ievt]
            qdata_high = quantiles[len(quantiles)-1]*1.2 # see comment above for mc            
         if dbg : print "data-quantiles  --> [ ", qdata_low ,qdata_high, " ]"



         # interplopate the correction
         y_corr = (qdata_high-qdata_low)/(qmc_high-qmc_low) * (Y[ievt] - qmc_low) + qdata_low
         if dbg : print "Apply correction: Input value = ", Y[ievt], " --> corrected value = ", y_corr

         y_tmp.append(y_corr)
         

      #   The isolation variables have a discrete single value at zero and then a smooth distribution
      #   To avoid degeneracies in the quantile matching I artificially smear them subracting rho
      #   (like a PU correction)
      if y == "PhoIso03" or y == "ChIso03" or y == "ChIso03worst":
         y_tmp = y_tmp + 0.1*self.df['rho']

      #self.y_corr = y_tmp
      ycorr = y+"_corr"
      self.df[ycorr] = y_tmp
      # print self.df[ycorr]

   # compute ID mvas with different sets of corrected variables
   # 
   # --------------------------------------------------------------------------------
   #
   def computeIdMvas(self,mvas,weights,n_jobs=1):
      wd = os.getcwd()

      weightsEB,weightsEE = map(lambda x: os.path.join(wd,x), weights )
      for name,correctedVariables in mvas:
         self.computeIdMva(name,wd,weightsEB,weightsEE,correctedVariables,n_jobs)

   def computeIdMva(self,name,wd,weightsEB,weightsEE,correctedVariables,n_jobs):
      stride = self.df.index.size / n_jobs
      print("Computing %s, correcting %s" % (name,correctedVariables) )
      Y = np.concatenate(Parallel(n_jobs=n_jobs,verbose=20)(
         delayed(computeIdMva)(wd,weightsEB,weightsEE,correctedVariables,self.df.loc[ch:ch+stride-1])
         for ch in xrange(0,self.df.index.size,stride) )
      )      

      ## Y = computeIdMva(wd,weightsEB,weightsEE,correctedVariables,self.df)

      self.df[name] = Y

      # compute ID mvas with different sets of corrected variables

   # correct photon isolation
   # 
   # --------------------------------------------------------------------------------
   def correctPhoIso(self,corr_file,n_jobs=1):
      wd = os.getcwd()
      stride = self.df.index.size / n_jobs
      corr_file = os.path.join(wd,corr_file)
      print("Computing corrected photon isolation using %s" % corr_file )
      rt.gROOT.LoadMacro(os.path.join(wd,"../phoIsoStoch/IsolationCorrection.C+"))
      Y = np.concatenate(Parallel(n_jobs=n_jobs,verbose=20)(
         delayed(applyIsoCorrection)(wd,corr_file,self.df.loc[ch:ch+stride-1])
         for ch in xrange(0,self.df.index.size,stride) )
      )      
      
      self.df["PhoIso03_corr"] = Y

   # get the names to be used for X and y and fill the corrected vector
   # 
   # --------------------------------------------------------------------------------
   #
   def correctYfast(self, x, y, quantiles, n_jobs=1, store=True ):
      
      dbg = False
      
      print "Get corrections for ", y, " with quantiles ", quantiles

      y_tmp = []
      
      # quantile regressions features
      X    = self.df.loc[:,x].values
      # target e.g. y = "R9"
      Y    = self.df[y]
      print "Features: X = ", x, " target y = ", y
      # print X, Y
      
      if y == 'PhoIso03' or y == 'ChIso03' or y == 'ChIso03worst':
         Y = Y - 0.1*self.df['rho']
      Y = Y.values.reshape(-1,1)
      Z = np.hstack([X,Y])
               
      Ycorr = np.concatenate(Parallel(n_jobs=n_jobs,verbose=20)(delayed(applyCorrection)(self.mcclf,self.dataclf,ch[:,:-1],ch[:,-1])
                                      for ch in np.array_split(Z,n_jobs) ) )
      
      if store:
         self.df[y+"_corr"] = Ycorr
         
      return Ycorr

     
   # get the corrected array of y
   # 
   # --------------------------------------------------------------------------------
   #
   def getCorrectedY(self, y):
      ycorr = y+"_corr"
      print 'Corrected Y = ', ycorr
      return self.df[ycorr]
      # return self.y_corr


   # produce a dataset with [X, all-Ycorr]
   # 
   # --------------------------------------------------------------------------------
   #
   def correctAllY(self, x, ylist, quantiles, n_jobs=1, forceComputeCorrections = False, EBEE="", relativePath=''):

      import os.path      
      corrTargetsName = 'correctedTargets'+relativePath
      if   EBEE == 'EB':
         corrTargetsName += '_EB'
      if   EBEE == 'EE':
         corrTargetsName += '_EE'
      corrTargetsName += '.h5'

      if ((os.path.exists(corrTargetsName)) and not( forceComputeCorrections )) :
         print 'Loading corrected targets from : ', corrTargetsName         
         self.df = pd.read_hdf(corrTargetsName, 'df')
         return
      
      else:          
         print 'Corrected variables file (e.g. ', corrTargetsName, ' ) does not exists. This will take a while...'

         # Here you are cutting out part of the DF !
         if   EBEE == 'EB':
            print "Correct EB :"
            self.applyCutsToDF('ScEta', -1.4442, 1.4442, 'inside')
         elif EBEE == 'EE':
            print "Correct EE :"
            self.applyCutsToDF('ScEta', -1.57, 1.57, 'outside')
         else:
            print "Correct both EB and EE together"
         
         mcfilename   = "./weights/"+relativePath+"/mc_weights"
         datafilename = "./weights/"+relativePath+"/data_weights"
         if   EBEE == 'EB':
            mcfilename   += "_EB"
            datafilename += "_EB"
         elif EBEE == 'EE':
            mcfilename   += "_EE"
            datafilename += "_EE"

         for Y in ylist:
            
            Yvar = Y
            if Y == 'PhoIso03' or Y == 'ChIso03' or Y == 'ChIso03worst':
               Y = Y + 'rho'
         
            print "Loading mc weights for ", Y, " : "
            print "   ", mcfilename
            self.loadMcWeights(mcfilename, Y, quantiles)      

            print "Loading data weights for ", Y
            print "   ", datafilename
            self.loadDataWeights(datafilename, Y, quantiles)      

            # print self.df
            #self.correctYfast(x, Yvar, quantiles, n_jobs=n_jobs )
            if ( (Y=="SigmaIeIe") and (EBEE=="EB")):
                self.df["SigmaIeIe_corr"]=self.df["SigmaIeIe"]
            else:
                self.correctY(x, Yvar, quantiles) #, n_jobs=n_jobs )

         if EBEE != '':
            print "Writing correctedTargets_",EBEE,".h5"
            hdf = pd.HDFStore('correctedTargets'+relativePath+'_'+EBEE+'.h5')
            hdf.put('df', self.df)
            hdf.close()
         else:
            print "Writing correctedTargets.h5"
            hdf = pd.HDFStore('correctedTargets'+relativePath+'.h5')
            hdf.put('df', self.df)
            hdf.close()
         
      #print 'Final datafame:', self.df

   # Scatter plots to check the quantiles
   # 
   # --------------------------------------------------------------------------------
   #
    #
   def correctAllYTime(self, x, ylist, quantiles, n_jobs=1, forceComputeCorrections = False, EBEE="", relativePath='',runperiod=''):

      import os.path      
      corrTargetsName = 'correctedTargetsPeriod'+str(runperiod)+relativePath
      if   EBEE == 'EB':
         corrTargetsName += '_EB'
      if   EBEE == 'EE':
         corrTargetsName += '_EE'
      corrTargetsName += '.h5'

      if ((os.path.exists(corrTargetsName)) and not( forceComputeCorrections )) :
         print 'Loading corrected targets from : ', corrTargetsName         
         self.df = pd.read_hdf(corrTargetsName, 'df')
         return
      
      else:          
         print 'Corrected variables file (e.g. ', corrTargetsName, ' ) does not exists. This will take a while...' 
      
         # Here you are cutting out part of the DF !
         if   EBEE == 'EB':
            print "Correct EB :"
            self.applyCutsToDF('ScEta', -1.4442, 1.4442, 'inside')
         elif EBEE == 'EE':
            print "Correct EE :"
            self.applyCutsToDF('ScEta', -1.57, 1.57, 'outside')
         else:
            print "Correct both EB and EE together"
         
         mcfilename   = "./weights/"+relativePath+"_"+str(runperiod)+"/mc_weights"
         datafilename = "./weights/"+relativePath+"_Tot/data_weights"
         if   EBEE == 'EB':
            mcfilename   += "_EB"
            datafilename += "_EB"
         elif EBEE == 'EE':
            mcfilename   += "_EE"
            datafilename += "_EE"

         for Y in ylist:

            Yvar = Y
            if Y == 'PhoIso03' or Y == 'ChIso03' or Y == 'ChIso03worst':
               Y = Y + 'rho'
         
            print "Loading mc weights for ", Y, " : "
            print "   ", mcfilename
            self.loadMcWeights(mcfilename, Y, quantiles)      

            print "Loading data weights for ", Y
            print "   ", datafilename
            self.loadDataWeights(datafilename, Y, quantiles)      

            self.correctYTime(x, Yvar, quantiles) #, n_jobs=n_jobs )

         if EBEE != '':
            print "Writing correctedTargets_",EBEE,".h5"
            hdf = pd.HDFStore('correctedTargets'+relativePath+'_'+EBEE+'.h5')
            hdf.put('df', self.df)
            hdf.close()
         else:
            print "Writing correctedTargets.h5"
            hdf = pd.HDFStore('correctedTargets'+relativePath+'.h5')
            hdf.put('df', self.df)
            hdf.close()

   def correctAllTime(self, x, ylist, quantiles, n_jobs=1, forceComputeCorrections = False, EBEE="", relativePath=''):
      df1=cp.deepcopy(self)
      df2=cp.deepcopy(self)
      df3=cp.deepcopy(self)
      df4=cp.deepcopy(self)
      df5=cp.deepcopy(self)
      df1.df=df1.df.query("runperiod==1")
      df2.df=df2.df.query("runperiod==2")
      df3.df=df3.df.query("runperiod==3")
      df4.df=df4.df.query("runperiod==4")
      df5.df=df5.df.query("runperiod==5")
      df1.correctAllYTime(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath,1)
      df2.correctAllYTime(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath,2)
      df3.correctAllYTime(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath,3)
      df4.correctAllYTime(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath,4)
      df5.correctAllYTime(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath,5)
      frame=[df1.df,df2.df,df3.df,df4.df,df5.df]
      dataframe=pd.concat(frame)
      self.df=dataframe.reset_index(drop=True)
   
   
   def correctAll(self, x, ylist, quantiles, n_jobs=1, forceComputeCorrections = False, EBEE="", relativePath=''):
      df1=cp.deepcopy(self)
      df2=cp.deepcopy(self)
      df3=cp.deepcopy(self)
      df4=cp.deepcopy(self)
      df5=cp.deepcopy(self)
      df1.df=df1.df.query("runperiod==1")
      df2.df=df2.df.query("runperiod==2")
      df3.df=df3.df.query("runperiod==3")
      df4.df=df4.df.query("runperiod==4")
      df5.df=df5.df.query("runperiod==5")
      df1.correctAllY(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath+str(1))
      df2.correctAllY(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath+str(2))
      df3.correctAllY(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath+str(3))
      df4.correctAllY(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath+str(4))
      df5.correctAllY(x, ylist, quantiles, n_jobs, forceComputeCorrections, EBEE, relativePath+str(5))
      frame=[df1.df,df2.df,df3.df,df4.df,df5.df]
      dataframe=pd.concat(frame)
      self.df=dataframe.reset_index(drop=True)
   
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


   # brute force access to X,Y, y_mc, y_data for debugging in notebooks
   #
   # --------------------------------------------------------------------------------
   #
   def forDebug(self, x, y, quantiles):#, X, Y, y_mc, y_data ):
      
      # quantile regressions features
      X    = self.df.loc[:,x]
      # target e.g. y = "R9"
      Y    = self.df[y]
      print "Features: X = ", x, " target y = ", y

      
      y_mc   = [] # list storing the n=q predictions on mc   for each event
      y_data = [] # list storing the n=q predictions on data for each event

      for q in range(0, len(quantiles)):
         y_mc  .append(self.mcclf[q]  .predict(X))
         y_data.append(self.dataclf[q].predict(X))

      return X, Y, y_mc, y_data 

   def getDF(self):#, X, Y, y_mc, y_data ):
      return self.df


   def getCLF(self):#, X, Y, y_mc, y_data ):
      return self.mcclf


