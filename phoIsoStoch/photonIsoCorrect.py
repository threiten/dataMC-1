from __future__ import print_function

import ROOT as RT
import numpy as np
import matplotlib.pyplot as plt
import math

import sys

# ------------------------------------------------------------------------------------------------------
class IO:
    
    # dataSrc = "root://t3dcachedb03.psi.ch//pnfs/psi.ch/cms/trivcat/store/user/musella/"
    dataSrc = "/mnt/t3nfs01/data01/shome/musella/Analysis/jupyter/dataMC/"
    objs = []
    
    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def getTree(fname,treeName,relPath=True):
        if relPath:
            fname = "%s/%s" % (IO.dataSrc,fname)
        fin=RT.TFile.Open(fname)
        if not fin:
            raise Exception("Could not open %s" % (fname))
        tree=fin.Get(treeName)
        try:
            tree.GetName()
        except:
            raise Exception("No tree %s found in %s" % (treeName,fname))
        IO.objs.append(fin)
        return tree
    
    
# ------------------------------------------------------------------------------------------------------
class Corrector:
    
    # ------------------------------------------------------------------------------------------------------
    def __init__(self,treeData=None,treeMC=None):
        self.treeData = treeData
        self.treeMC   = treeMC
        
        RT.RooMsgService.instance().setGlobalKillBelow(RT.RooFit.FATAL)
        self.reset()
        
    # ------------------------------------------------------------------------------------------------------
    def reset(self):
        self.xvar    = None
        self.multvar = None
        self.morph   = None

        self.rooData = None
        self.rooMC   = None

        self.histData = None
        self.histMC   = None
        
        self.rhobins  = None
        self.rhoData  = None
        self.rhoMC    = None

        self.mults    = []

    # ------------------------------------------------------------------------------------------------------
    def mkIsoHistos(self, rhobins, var, name, binning, weight="weight"):
        
        self.reset()

        self.getXvar(bounds=[binning[0],binning[0],binning[-1]]).setBinning(
            RT.RooBinning(binning.size-1,binning))
                
        ret = []
        rhoprofs = []
        for tree,pfx in (self.treeData,"hdata_"),(self.treeMC,"hmc_"):
            histos = []
            print(tree,pfx,file=sys.stderr)
            if tree == None:
                print("Tree %s is None" % pfx[1:].replace("_",""), file=sys.stderr)
                ret.append(histos)
            if len(rhobins.shape) > 1:
                rhopairs = rhobins
            else:
                rhopairs = zip(rhobins,rhobins[1:])
            for rhomin,rhomax in rhopairs:
                hname = pfx + (name % (rhomin, rhomax)).replace(".","p")
                print(hname,file=sys.stderr)
                hist = RT.TH1D(hname, hname, binning.size-1, binning)
                
                tree.Draw("%s>>+%s" % (var,hname), "(%s)*(rho > %f && rho <= %f)" % (weight, rhomin, rhomax), "goff")
                hist.SetDirectory(0)
                
                histos.append(hist)
            ret.append(histos)
            
            hrhoname = pfx+"rho"
            if len(rhobins.shape) == 1:
                rhoProfile = RT.TProfile(hrhoname,hrhoname,rhobins.size-1,rhobins)
                tree.Draw("rho:rho>>+%s" % hrhoname, weight, "goff")
            else:
                rhoProfile = RT.TH1D(hrhoname,hrhoname,rhobins.shape[0],np.linspace(0,rhobins.shape[0]+1,rhobins.shape[0]))
                ibin = 0
                for rhomin,rhomax in rhopairs:
                    ibin += 1
                    tmp = RT.TH1D("tmp", "tmp", 100, np.linspace(0,100,101))
                    nentries = tree.Draw("rho>>+tmp", "(%s)*(rho > %f && rho <= %f)" % (weight, rhomin, rhomax), "goff")
                    tmp.SetDirectory(0)
                    rhoProfile.SetBinContent(ibin,tmp.GetMean())
                    rhoProfile.SetBinError(ibin,tmp.GetRMS())
            rhoProfile.SetDirectory(0)
            rhoprofs.append(rhoProfile)
            

        self.histData = ret[0]
        self.histMC   = ret[1]

        self.rhobins = rhobins
        ## if len(rhobins.shape) == 1:
        self.rhoData = rhoprofs[0]
        self.rhoMC   = rhoprofs[1]
            
        return ret

    # ------------------------------------------------------------------------------------------------------
    def getXvar(self,bounds=[0.,0.,4.]):
        if not self.xvar:
            self.xvar = RT.RooRealVar("x","x",*bounds)
        return self.xvar

    # ------------------------------------------------------------------------------------------------------
    def mkDatasets(self):
        xvar = self.xvar
        
        self.rooData = map(lambda x: Corrector.buildDataset(x,xvar), self.histData )
        self.rooMC   = map(lambda x: Corrector.buildDataset(x,xvar), self.histMC   )
        
        
    # ------------------------------------------------------------------------------------------------------
    def buildPdf(self,grid,nentries,imult=0,returnHists=False):
        
        print("buildPdf",grid,nentries,imult)
        interpolated = Corrector.interpolate(self.histMC[imult],grid,nentries)
        
        xvar = self.getXvar()
        xlist = RT.RooArgList(xvar)
        
        self.interpPdfs = []
        self.interpDsets = []
        for h in interpolated:
            minval = h.Integral()*1.e-6
            Corrector.protect(h,minval)
            ds = RT.RooDataHist(h.GetName(),h.GetName(),xlist,h)
            self.interpDsets.append(ds)
            pdf = RT.RooHistPdf("%s_pdf" % h.GetName(), "%s_pdf" % h.GetName(), RT.RooArgSet(xlist), ds)
            self.interpPdfs.append(pdf)
            
        self.multvar = RT.RooRealVar("mult","mult",grid[0],grid[0],grid[-1])    
        self.morph   = Corrector.mkMomentMorph("morph",xlist,self.interpPdfs,self.multvar,grid)
        
        if returnHists:
            return self.morph,interpolated            
        return self.morph

    # ------------------------------------------------------------------------------------------------------
    def fitMultiplicity(self,imult,sameMult=False,multOffset=0,buildOpts={},plot=False):
        if not self.rooData or not self.rooMC:
            self.mkDatasets()
        if sameMult: 
            buildOpts["imult"] = max(0,imult-multOffset)
        if not self.morph or sameMult:
            self.buildPdf(**buildOpts)
            
        xvar = self.getXvar()
         
        data = self.rooData[imult]
        mc   = self.rooMC[imult]
        
        multData = Corrector.fit(data,self.morph,self.multvar)
        if plot:
            frameData = xvar.frame()
            mcpdf =  RT.RooHistPdf(mc.GetName(),mc.GetName(),RT.RooArgSet(xvar), mc)
            data.plotOn(frameData,RT.RooFit.MarkerColor(RT.kBlack),RT.RooFit.Invisible())
            self.morph.plotOn(frameData,RT.RooFit.LineColor(RT.kRed))
            mcpdf.plotOn(frameData,RT.RooFit.LineColor(RT.kBlue))
            data.plotOn(frameData,RT.RooFit.MarkerColor(RT.kBlack))
        
        multMC   = Corrector.fit(mc,self.morph,self.multvar)
        if plot:
            frameMC = xvar.frame()
            mc.plotOn(frameMC,RT.RooFit.Invisible())
            self.morph.plotOn(frameMC,RT.RooFit.LineColor(RT.kRed))
            mc.plotOn(frameMC,RT.RooFit.MarkerStyle(RT.kOpenCircle),RT.RooFit.MarkerColor(RT.kBlue))
        
        
        rhoAvgData = self.rhoData.GetBinContent(imult+1)
        rhoAvgMC   = self.rhoMC.GetBinContent(imult+1)
            
        self.mults.append( (imult,rhoAvgData,rhoAvgMC,multData,multMC) )
        
        if plot:
            return (frameData,frameMC),self.mults[-1]
            
        return self.mults[-1]
        
    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def fit(data,pdf,poi):
        
        nll = pdf.createNLL(data)
        minim = RT.RooMinimizer(nll)
        # minim.setMinimizerType("Minuit2")
        # minim.setStrategy(2)
        minim.setPrintLevel(-1)
        minim.migrad()
        minim.minos()
        
        return poi.getVal(),poi.getAsymErrorLo(),poi.getAsymErrorHi()

    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def protect(h, thr):
        for ibin in xrange(h.GetNbinsX()):
            cont = h.GetBinContent(ibin+1)
            if cont < thr: h.SetBinContent(ibin+1,thr)

    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def generate_one(h, thr):
        val = h.GetRandom()
        if val < thr: return 0.
        return val

    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def generate_one_smear(h, thr, sigma):
        val = h.GetRandom() +  RT.gRandom.Gaus(0.,sigma)
        if val < thr: return 0.
        return val

    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def generate_n(h, n, thr):
        gen = np.array(map(lambda x: Corrector.generate_one(h,thr), xrange(n)))
        # print(gen,gen.sum())
        return gen.sum()

    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def generate_n_smear(h, n, thr, sigma):
        gen = np.array(map(lambda x: Corrector.generate_one_smear(h,thr,sigma), xrange(n)))
        # print(gen,gen.sum())
        return gen.sum()
        
    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def generate(base,mult,nentries,smear=None,usePoisson=False):
        
        name = ("%s_mult_%1.2f" % (base.GetName(),mult)) if not smear else ("%s_mult_%1.2f_%1.2f" % (base.GetName(),mult,smear))
        name = name.replace(".","p")
        generated = base.Clone( name  )
        generated.Reset("ICE")
        
        thr = base.GetXaxis().GetBinLowEdge(2)
        if usePoisson:
            ngen = lambda : RT.gRandom.Poisson(mult)
            name += "_poiss"
        else:
            floor = math.floor(mult)
            remind = mult - floor 
            floor = int(floor)
            mult = int(mult)
            if remind > 0:
                ngen = lambda : floor if RT.gRandom.Uniform() < remind else floor+1
            else:
                ngen = lambda : mult    
        
        if smear:
            gen_n = lambda : Corrector.generate_n_smear(base,ngen(),thr,smear)
        else:
            gen_n = lambda : Corrector.generate_n(base,ngen(),thr)
        map(lambda x: generated.Fill(gen_n()), xrange(nentries) )
        
        return generated
        
    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def interpolate(base,grid,nentries,**kwargs):
        return map(lambda x: Corrector.generate(base,x,nentries,**kwargs), grid)
        
    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def mkMomentMorph(name,X,pdfs,theta,values):
        rooPdfs = RT.RooArgList()
        rooValues = RT.TVectorD(len(values))
        for ival,(val,pdf) in enumerate(zip(values,pdfs)):
            rooPdfs.add(pdf)
            rooValues[ival] = val
        morph = RT.RooMomentMorph(name,name,theta,X,rooPdfs,rooValues,RT.RooMomentMorph.Linear)
        morph.useHorizontalMorphing(False)
        IO.objs.append(rooPdfs)
        return morph
        

    # ------------------------------------------------------------------------------------------------------
    @staticmethod
    def buildDataset(h,X):
        return RT.RooDataHist(h.GetName(),h.GetName(),RT.RooArgList(X),h)
        

# ------------------------------------------------------------------------------------------------------
def runEtaBin(treeData,treeMC,rhobins,isobins,etaLabel,etaCut,sameMult=False,multOffset=0,drawTH1s=False,histosOnly=False):

    # print(treeData,treeMC,rhobins,isobins,etaLabel,etaCut,sameMult,multOffset,file=sys.stderr)
    if treeData:
        treeData = IO.getTree(*treeData)
    else:
        print("No tree for data. Will just make MC histograms.", file=sys.stderr)
        histosOnly = True
    treeMC = IO.getTree(*treeMC)
    ## print (treeData, file=sys.stderr)
    ## print (treeMC,file=sys.stderr)
    
    ## print("making corrector",file=sys.stderr)
    correct = Corrector(treeData,treeMC)

    ## histos = correct.mkIsoHistos(rhobins,"probe_phoiso",etaLabel+"_%1.2f_%1.2f",isobins,"weight*(%s)" % etaCut)
    ## print("making histos",file=sys.stderr)
    histos = correct.mkIsoHistos(rhobins,"PhoIso03",etaLabel+"_%1.2f_%1.2f",isobins,"1.*(%s)" % etaCut)
    ## print("done",file=sys.stderr)

    
    if histosOnly:
        return etaLabel,histos
    
    ## print(zip(rhobins,rhobins[1:]))
    ## 
    ## print(isobins)
    ## 
    ## print(histos)

    sqrtn = int(math.ceil(math.sqrt(len(histos[0]))))
    ## canvSize = min(2,sqrtn)*100
    canvSize = min(2,sqrtn)*200
    if drawTH1s:
        canv = RT.TCanvas("c","c",canvSize*7,canvSize*5)
        canv.Divide(sqrtn,sqrtn)
        for ih,hists in enumerate(zip(histos[0],histos[1])):
            data,mc = hists
            canv.cd(ih+1)            
            mc.GetXaxis().SetTitle("Iso_{#gamma}")
            mc.GetYaxis().SetTitle("A.U.")
            data.GetXaxis().SetTitle("Iso_{#gamma}")
            data.GetYaxis().SetTitle("A.U.")
            
            mc.SetLineColor(RT.kBlue)
            mc.DrawNormalized("hist")
            data.SetMarkerColor(RT.kBlack)
            data.SetLineColor(RT.kBlack)
            data.DrawNormalized("pe same")
        
        canv.Draw()
        for ih in xrange(len(histos[0])):
            canv.cd(ih+1)            
            RT.gPad.SetLogy()
        canv.Draw()

    if not sameMult:
        grid = np.linspace(0.5,10,20)
    else:
        grid = np.linspace(0.5,4,8)
        
    # frame = self.getXvar().frame()
    allframes = []
    canvData=RT.TCanvas("cd","cd",canvSize*7,canvSize*5)
    canvData.Divide(sqrtn,sqrtn)
    canvMC=RT.TCanvas("ce","ce",canvSize*7,canvSize*5)
    canvMC.Divide(sqrtn,sqrtn)
    def fit(imult):
        frames, mults = correct.fitMultiplicity(imult,sameMult=sameMult,
                                                multOffset=multOffset,
                                                buildOpts=dict(grid=grid,nentries=500000),
                                                plot=True)
        #mults = correct.fitMultiplicity(imult,plot=False)

        canvData.cd(imult+1)
        frames[0].Draw()
        frames[0].SetTitle("Data - <#rho> = %1.1f - n_{eff} = %1.1f +-%1.2f" % (mults[1],mults[3][0],0.5*(abs(mults[3][1])+abs(mults[3][2]))) )
        frames[0].GetXaxis().SetTitle("Iso_{#gamma}")
        # frames[0].GetYaxis().SetTitle("Events / bin")
        
        canvMC.cd(imult+1)
        frames[1].Draw()
        frames[1].SetTitle("MC - <#rho> = %1.1f - n_{eff} = %1.1f +- %1.0f" % (mults[2],mults[4][0],0.5*(abs(mults[4][1])+abs(mults[4][2]))) )
        frames[1].GetXaxis().SetTitle("Iso_{#gamma}")
        # frames[1].GetYaxis().SetTitle("Events / bin")

        allframes.append(frames)
        return mults

    
    mults = map( fit,  xrange(rhobins.shape[0]) )
    canvData.Draw()
    # canvMC.Draw()
    for ih in xrange(rhobins.shape[0]):
        canvData.cd(ih+1)
        RT.gPad.SetLogy()
        canvMC.cd(ih+1)
        RT.gPad.SetLogy()
    canvData.Draw()
    canvMC.Draw()
    
    npmults = np.array(map(lambda x: (x[1], x[2], x[3][0], x[4][0], (x[1]*x[4][0])/(x[2]*x[3][0])), mults))
    # print(npmults)
    mdata,qdata = np.polyfit(npmults[:,0],npmults[:,2],1)
    mmc, qmc = np.polyfit(npmults[:,1],npmults[:,3],1)
    
    plt.plot(npmults[:,0],npmults[:,2],'o',label="Data %1.2f + %1.2f<rho>" % (qdata,mdata))
    plt.plot(npmults[:,1],npmults[:,3],'o',label="MC %1.2f + %1.2f<rho>" % (qmc,mmc))
    
    
    plt.plot(npmults[:,0], mdata*npmults[:,0] + qdata, 'b-')
    plt.plot(npmults[:,1], mmc*npmults[:,1] + qmc, 'g-')
    
    plt.title(etaLabel)
    plt.xlabel("<rho>")
    plt.ylabel("effective multiplicity")
    plt.legend(loc='best')
    plt.show()

    plt.plot(npmults[:,0],npmults[:,2]/npmults[:,3],'o-')
    
    plt.title(etaLabel)
    plt.xlabel("<rho>")
    plt.ylabel("Data / MC")
    plt.legend(loc='best')
    plt.show()

    return etaLabel,rhobins.size,isobins[0],isobins[-1],sameMult,multOffset,npmults[:,2]/npmults[:,3],npmults[:,3],npmults[:,0],npmults[:,1],mdata,qdata,mmc,qmc,histos
    
