#include <iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "TRandom.h"

TRandom rndm;

float generateGaus(float min, float max, float mean = 0, float width = 1)
{
  float x = -9999;
  while (x <min || x > max){
    x = rndm.Gaus(mean, width);
  }
  return x;
}

float generateFlat(float min, float max)
{
  float x = rndm.Uniform(min, max);
  return x;
}



int toy(int nevt = 1000)
{
  rndm.SetSeed(1);

  // same for data and MC
  float etaMin   = -2.5;
  float etaMax   =  2.5;
  float etaWidth =  1.0;
  float phiMin   = -3.14;
  float phiMax   =  3.14;
  float ptMin    =  10;
  float ptMax    =  200;
  float ptWidth  =   20;

  // data
  float R9Min    =   0.8; 
  float R9Max    =   1.0; 
  float R9MeanData   =   0.95;
  float R9WidthData  =   0.01;
  float sieieMin    =   0;
  float sieieMax    =   0.02;
  float sieieMeanData   =  0.008;
  float sieieWidthData  =  0.001;

  // mc
  float R9MeanMC   =   0.96;
  float R9WidthMC  =   0.01;
  float sieieMeanMC  =  0.007;
  float sieieWidthMC =  0.001;


  // eta dependences
  TF1 *fR9etaData = new TF1("fR9etaData","[0]*x*x+[1]",-2.5,2.5);
  fR9etaData->SetParameters(-0.024,R9MeanData);
  TF1 *fsieieetaData = new TF1("fsieieetaData","[0]*x*x+[1]",-2.5,2.5);
  fsieieetaData->SetParameters(-0.00032,sieieMeanData);

  TF1 *fR9etaMC = new TF1("fR9etaMC","[0]*x*x+[1]",-2.5,2.5);
  fR9etaMC->SetParameters(-0.030,R9MeanMC);
  TF1 *fsieieetaMC = new TF1("fsieieetaMC","[0]*x*x+[1]",-2.5,2.5);
  fsieieetaMC->SetParameters(-0.00064,sieieMeanMC);

  TH1F* hdata_eta   = new TH1F("hdata_eta"  ,"hdata_eta"  , 100, etaMin   , etaMax);
  TH1F* hdata_phi   = new TH1F("hdata_phi"  ,"hdata_phi"  , 100, phiMin   , phiMax);
  TH1F* hdata_pt    = new TH1F("hdata_pt"   ,"hdata_pt"   , 100, ptMin    , ptMax );
  TH1F* hdata_R9    = new TH1F("hdata_R9"   ,"hdata_R9"   , 100, R9Min    , R9Max);
  TH1F* hdata_sieie = new TH1F("hdata_sieie","hdata_sieie", 100, sieieMin , sieieMax);

  TH1F* hMC_eta   = new TH1F("hMC_eta"  ,"hMC_eta"  , 100, etaMin, etaMax);
  TH1F* hMC_phi   = new TH1F("hMC_phi"  ,"hMC_phi"  , 100, phiMin, phiMax);
  TH1F* hMC_pt    = new TH1F("hMC_pt"   ,"hMC_pt"   , 100, ptMin,  ptMax );
  TH1F* hMC_R9    = new TH1F("hMC_R9"   ,"hMC_R9"   , 100, R9Min    , R9Max);
  TH1F* hMC_sieie = new TH1F("hMC_sieie","hMC_sieie", 100, sieieMin , sieieMax);

  TH2F *hh_etaR9data = new TH2F("hh_etaR9data","hh_etaR9data", 100, etaMin, etaMax, 100, R9Min, R9Max);
  TH2F *hh_etaR9MC   = new TH2F("hh_etaR9MC"  ,"hh_etaR9MC"  , 100, etaMin, etaMax, 100, R9Min, R9Max);
  TH2F *hh_etasieiedata = new TH2F("hh_etasieiedata","hh_etasieiedata", 100, etaMin, etaMax, 100, sieieMin, sieieMax);
  TH2F *hh_etasieieMC   = new TH2F("hh_etasieieMC"  ,"hh_etasieieMC"  , 100, etaMin, etaMax, 100, sieieMin, sieieMax);
  
  TFile *f = new TFile("dataMC.root","RECREATE","Data and MC samples");

  TTree *data  = new TTree("data","Data events");

  float eta   = -99999;
  float phi   = -99999;
  float pt    = -99999;
  float R9    = -99999;
  float sieie = -99999;

  data->Branch("eta"   , &eta   ," eta/F");
  data->Branch("phi"   , &phi   , "phi/F");
  data->Branch("pt"    , &pt    , "pt/F");
  data->Branch("R9"    , &R9    , "R9/F");
  data->Branch("sieie" , &sieie , "sieie/F");
  
  int ndata = nevt;
  for (int i = 0 ; i < ndata ; ++i){   

    pt    = generateGaus(ptMin , ptMax,  0, ptWidth);
    eta   = generateGaus(etaMin, etaMax, 0, etaWidth);
    phi   = generateFlat(phiMin, phiMax);
    R9    = generateGaus(R9Min , R9Max,  fR9etaData->Eval(eta), R9WidthData);
    sieie = generateGaus(sieieMin , sieieMax,  fsieieetaData->Eval(eta), sieieWidthData);
    data->Fill();
    hdata_eta   ->Fill(eta);
    hdata_phi   ->Fill(phi);
    hdata_pt    ->Fill(pt );
    hdata_R9    ->Fill(R9 );
    hdata_sieie ->Fill(sieie );
    hh_etaR9data->Fill(eta,R9);
    hh_etasieiedata->Fill(eta,sieie);
  }
  data->Write();
  data->Print();
  
  TTree *MC  = new TTree("mc","MC events");
  MC->Branch("eta"   , &eta   ," eta/F");
  MC->Branch("phi"   , &phi   , "phi/F");
  MC->Branch("pt"    , &pt    , "pt/F");
  MC->Branch("R9"    , &R9    , "R9/F");
  MC->Branch("sieie" , &sieie , "sieie/F");
  
  int nMC = nevt;
  for (int i = 0 ; i < nMC ; ++i){   

    eta   = generateGaus(etaMin, etaMax, 0, etaWidth);
    phi   = generateFlat(phiMin, phiMax);
    pt    = generateGaus(ptMin , ptMax,  0, ptWidth);
    R9    = generateGaus(R9Min , R9Max,  fR9etaMC->Eval(eta), R9WidthMC);
    sieie = generateGaus(sieieMin , sieieMax,  fsieieetaMC->Eval(eta), sieieWidthMC);
    MC->Fill();
    hMC_eta   ->Fill(eta);
    hMC_phi   ->Fill(phi);
    hMC_pt    ->Fill(pt );
    hMC_R9    ->Fill(R9 );
    hMC_sieie ->Fill(sieie );
    hh_etaR9MC->Fill(eta,R9);
    hh_etasieieMC->Fill(eta,sieie);    
  }
  MC->Write();
  MC->Print();

  f->Close();
  
  TCanvas *c = new TCanvas("c","c",600,600);
  c->cd();
  c->Divide(3,3);
  c->cd(1);
  hdata_eta->SetLineColor(2);
  hdata_eta->Draw();
  hMC_eta->SetLineColor(4);
  hMC_eta->Draw("same");
  c->cd(2);
  hdata_phi->SetMinimum(0);
  hdata_phi->SetLineColor(2);
  hdata_phi->Draw();
  hMC_phi->SetLineColor(4);
  hMC_phi->Draw("same");
  c->cd(3);
  hdata_pt->SetLineColor(2);
  hdata_pt->Draw();
  hMC_pt->SetLineColor(4);
  hMC_pt->Draw("same");
  c->cd(4);
  hdata_R9->SetLineColor(2);
  hdata_R9->Draw();
  hMC_R9->SetLineColor(4);
  hMC_R9->Draw("same");
  c->cd(5);
  hdata_sieie->SetLineColor(2);
  hdata_sieie->Draw();
  hMC_sieie->SetLineColor(4);
  hMC_sieie->Draw("same");

  c->cd(7);
  hh_etaR9data->SetMarkerColor(2);
  hh_etaR9data->SetMarkerStyle(24);
  hh_etaR9data->Draw();
  hh_etaR9MC->SetMarkerColor(4);
  hh_etaR9MC->SetMarkerStyle(24);
  hh_etaR9MC->Draw("same");
  c->cd(8);
  hh_etasieiedata->SetMarkerColor(2);
  hh_etasieiedata->SetMarkerStyle(24);
  hh_etasieiedata->Draw();
  hh_etasieieMC->SetMarkerColor(4);
  hh_etasieieMC->SetMarkerStyle(24);
  hh_etasieieMC->Draw("same");

  c->Update();
  getchar();

  return 0;
}
