#define MiniTree_cxx
#include "MiniEffi.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "HistogramManager.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>

#include "TTree.h"

void SetSebStyle()
{
 gStyle->SetTitleFillColor(42);
 gStyle->SetTitleFont(1);
 gStyle->SetStatColor(29);
 gStyle->SetCanvasColor(25);   
 gStyle->SetOptStat(1111111);
 gStyle->SetHistFillColor(5);
}

void SaveInFile(TH1* ahisto, TFile* afile)
{
 if (!ahisto) { std::cout << "!! no histo !!" << std::endl; return ;}
 TDirectory* current = gDirectory ;
 afile->cd();
 ahisto->Write();
 current->cd();
}

void MiniTree::Loop(int aNN, float aTagCut, float aPtMin, float aPtMax, 
                float aEtaMin, float aEtaMax, float aFreeCut, int aIntCut, 
		TString afilename, TString aweightFileMVA)
{
///////////////////////////////////////////////////////////////////
//                  T I T L E      C A R D S                     //

//$$
  int    NN           = aNN;      // 0/1 JP, 2/3 TC, 4/5 SSV, 6 CSV, 7 SL
  float  TagCut       = aTagCut;  // tag cut
  float  PtMin        = aPtMin;   // pt jet min
  float  PtMax        = aPtMax;   // pt jet max
  float  EtaMin       = aEtaMin;  // eta jet min
  float  EtaMax       = aEtaMax;  // eta jet max
  float  FreeCut      = aFreeCut; 
  int    IntCut       = aIntCut; 
  TString  filename   = afilename;
  TString  weightFile_ = aweightFileMVA;
//$$

//$$
  float HTcut =  aTagCut; 
//   float dRcut =  FreeCut; 
//$$

// to recompute the BDT output from Ntuple info
  bool ReadXml = false;
  
//////////////// Preselection cuts:
// default:
  int LepNumCut = 100;
  float LepPtCut = 15., ElEtaCut = 2.5;
  float MuPtCut1 = 25., MuPtCut2 = 10., MuEtaCut = 10.;
  float PtmumuCut = 0., MmumuCutInf = 10., MmumuCutSup = 1000.;
  int nJetCut = 0;
  float JetPtCut = 20., JetEtaCut = 2.4, HtCut = 0.;
  float JetPtCut4 = 0., JetPtCut5 = 0., JetPtCut6 = 0., JetPtCut7 = 0.; 
  float SVrCut = 73., SVzCut = 186.;

///////////////////////////////////////////////////////////////////
 
//$$
//$$ TH1F::SetDefaultSumw2(kTRUE);
//$$

//**********************************
// Histograms
//**********************************

 TH1F* hData_Hemi_SecVtx_size = new TH1F("hData_Hemi_SecVtx_size","",10,-0.5,9.5);

 TH1F* hGen_Msmu = new TH1F("hGen_Msmu","",361,159.5,520.5);
 TH1F* hGen_Mneu = new TH1F("hGen_Mneu","",361,159.5,520.5);
 TH1F* hGen_ct0  = new TH1F("hGen_ct0","",200,0.,200.);

 TH1F* hGen_ct0_smu200_neu180 = new TH1F("hGen_ct0_smu200_neu180","",200,0.,200.);
 TH1F* hGen_ct0_smu250_neu180 = new TH1F("hGen_ct0_smu250_neu180","",200,0.,200.);
 TH1F* hGen_ct0_smu250_neu200 = new TH1F("hGen_ct0_smu250_neu200","",200,0.,200.);
 TH1F* hGen_ct0_smu250_neu230 = new TH1F("hGen_ct0_smu250_neu230","",200,0.,200.);
 TH1F* hGen_ct0_smu300_neu180 = new TH1F("hGen_ct0_smu300_neu180","",200,0.,200.);
 TH1F* hGen_ct0_smu300_neu200 = new TH1F("hGen_ct0_smu300_neu200","",200,0.,200.);
 TH1F* hGen_ct0_smu300_neu250 = new TH1F("hGen_ct0_smu300_neu250","",200,0.,200.);
 TH1F* hGen_ct0_smu300_neu280 = new TH1F("hGen_ct0_smu300_neu280","",200,0.,200.);
 TH1F* hGen_ct0_smu350_neu180 = new TH1F("hGen_ct0_smu350_neu180","",200,0.,200.);
 TH1F* hGen_ct0_smu350_neu200 = new TH1F("hGen_ct0_smu350_neu200","",200,0.,200.);
 TH1F* hGen_ct0_smu350_neu250 = new TH1F("hGen_ct0_smu350_neu250","",200,0.,200.);
 TH1F* hGen_ct0_smu350_neu300 = new TH1F("hGen_ct0_smu350_neu300","",200,0.,200.);
 TH1F* hGen_ct0_smu350_neu330 = new TH1F("hGen_ct0_smu350_neu330","",200,0.,200.);
 TH1F* hGen_ct0_smu400_neu180 = new TH1F("hGen_ct0_smu400_neu180","",200,0.,200.);
 TH1F* hGen_ct0_smu400_neu200 = new TH1F("hGen_ct0_smu400_neu200","",200,0.,200.);
 TH1F* hGen_ct0_smu400_neu250 = new TH1F("hGen_ct0_smu400_neu250","",200,0.,200.);
 TH1F* hGen_ct0_smu400_neu300 = new TH1F("hGen_ct0_smu400_neu300","",200,0.,200.);
 TH1F* hGen_ct0_smu400_neu350 = new TH1F("hGen_ct0_smu400_neu350","",200,0.,200.);
 TH1F* hGen_ct0_smu400_neu380 = new TH1F("hGen_ct0_smu400_neu380","",200,0.,200.);
 TH1F* hGen_ct0_smu450_neu180 = new TH1F("hGen_ct0_smu450_neu180","",200,0.,200.);
 TH1F* hGen_ct0_smu450_neu200 = new TH1F("hGen_ct0_smu450_neu200","",200,0.,200.);
 TH1F* hGen_ct0_smu450_neu250 = new TH1F("hGen_ct0_smu450_neu250","",200,0.,200.);
 TH1F* hGen_ct0_smu450_neu300 = new TH1F("hGen_ct0_smu450_neu300","",200,0.,200.);
 TH1F* hGen_ct0_smu450_neu350 = new TH1F("hGen_ct0_smu450_neu350","",200,0.,200.);
 TH1F* hGen_ct0_smu450_neu400 = new TH1F("hGen_ct0_smu450_neu400","",200,0.,200.);
 TH1F* hGen_ct0_smu450_neu430 = new TH1F("hGen_ct0_smu450_neu430","",200,0.,200.);
 TH1F* hGen_ct0_smu500_neu180 = new TH1F("hGen_ct0_smu500_neu180","",200,0.,200.);
 TH1F* hGen_ct0_smu500_neu200 = new TH1F("hGen_ct0_smu500_neu200","",200,0.,200.);
 TH1F* hGen_ct0_smu500_neu250 = new TH1F("hGen_ct0_smu500_neu250","",200,0.,200.);
 TH1F* hGen_ct0_smu500_neu300 = new TH1F("hGen_ct0_smu500_neu300","",200,0.,200.);
 TH1F* hGen_ct0_smu500_neu350 = new TH1F("hGen_ct0_smu500_neu350","",200,0.,200.);
 TH1F* hGen_ct0_smu500_neu400 = new TH1F("hGen_ct0_smu500_neu400","",200,0.,200.);
 TH1F* hGen_ct0_smu500_neu450 = new TH1F("hGen_ct0_smu500_neu450","",200,0.,200.);
 TH1F* hGen_ct0_smu500_neu480 = new TH1F("hGen_ct0_smu500_neu480","",200,0.,200.);

 TH1F* hSim_Filter_smu200_neu180 = new TH1F("hSim_Filter_smu200_neu180","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu250_neu180 = new TH1F("hSim_Filter_smu250_neu180","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu250_neu200 = new TH1F("hSim_Filter_smu250_neu200","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu250_neu230 = new TH1F("hSim_Filter_smu250_neu230","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu300_neu180 = new TH1F("hSim_Filter_smu300_neu180","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu300_neu200 = new TH1F("hSim_Filter_smu300_neu200","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu300_neu250 = new TH1F("hSim_Filter_smu300_neu250","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu300_neu280 = new TH1F("hSim_Filter_smu300_neu280","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu350_neu180 = new TH1F("hSim_Filter_smu350_neu180","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu350_neu200 = new TH1F("hSim_Filter_smu350_neu200","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu350_neu250 = new TH1F("hSim_Filter_smu350_neu250","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu350_neu300 = new TH1F("hSim_Filter_smu350_neu300","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu350_neu330 = new TH1F("hSim_Filter_smu350_neu330","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu400_neu180 = new TH1F("hSim_Filter_smu400_neu180","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu400_neu200 = new TH1F("hSim_Filter_smu400_neu200","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu400_neu250 = new TH1F("hSim_Filter_smu400_neu250","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu400_neu300 = new TH1F("hSim_Filter_smu400_neu300","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu400_neu350 = new TH1F("hSim_Filter_smu400_neu350","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu400_neu380 = new TH1F("hSim_Filter_smu400_neu380","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu450_neu180 = new TH1F("hSim_Filter_smu450_neu180","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu450_neu200 = new TH1F("hSim_Filter_smu450_neu200","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu450_neu250 = new TH1F("hSim_Filter_smu450_neu250","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu450_neu300 = new TH1F("hSim_Filter_smu450_neu300","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu450_neu350 = new TH1F("hSim_Filter_smu450_neu350","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu450_neu400 = new TH1F("hSim_Filter_smu450_neu400","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu450_neu430 = new TH1F("hSim_Filter_smu450_neu430","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu500_neu180 = new TH1F("hSim_Filter_smu500_neu180","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu500_neu200 = new TH1F("hSim_Filter_smu500_neu200","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu500_neu250 = new TH1F("hSim_Filter_smu500_neu250","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu500_neu300 = new TH1F("hSim_Filter_smu500_neu300","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu500_neu350 = new TH1F("hSim_Filter_smu500_neu350","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu500_neu400 = new TH1F("hSim_Filter_smu500_neu400","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu500_neu450 = new TH1F("hSim_Filter_smu500_neu450","",2,-0.5,1.5);
 TH1F* hSim_Filter_smu500_neu480 = new TH1F("hSim_Filter_smu500_neu480","",2,-0.5,1.5);

 TH1F* hSim_Hemi_dR_smu200_neu180 = new TH1F("hSim_Hemi_dR_smu200_neu180","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu250_neu180 = new TH1F("hSim_Hemi_dR_smu250_neu180","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu250_neu200 = new TH1F("hSim_Hemi_dR_smu250_neu200","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu250_neu230 = new TH1F("hSim_Hemi_dR_smu250_neu230","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu300_neu180 = new TH1F("hSim_Hemi_dR_smu300_neu180","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu300_neu200 = new TH1F("hSim_Hemi_dR_smu300_neu200","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu300_neu250 = new TH1F("hSim_Hemi_dR_smu300_neu250","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu300_neu280 = new TH1F("hSim_Hemi_dR_smu300_neu280","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu350_neu180 = new TH1F("hSim_Hemi_dR_smu350_neu180","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu350_neu200 = new TH1F("hSim_Hemi_dR_smu350_neu200","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu350_neu250 = new TH1F("hSim_Hemi_dR_smu350_neu250","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu350_neu300 = new TH1F("hSim_Hemi_dR_smu350_neu300","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu350_neu330 = new TH1F("hSim_Hemi_dR_smu350_neu330","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu400_neu180 = new TH1F("hSim_Hemi_dR_smu400_neu180","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu400_neu200 = new TH1F("hSim_Hemi_dR_smu400_neu200","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu400_neu250 = new TH1F("hSim_Hemi_dR_smu400_neu250","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu400_neu300 = new TH1F("hSim_Hemi_dR_smu400_neu300","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu400_neu350 = new TH1F("hSim_Hemi_dR_smu400_neu350","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu400_neu380 = new TH1F("hSim_Hemi_dR_smu400_neu380","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu450_neu180 = new TH1F("hSim_Hemi_dR_smu450_neu180","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu450_neu200 = new TH1F("hSim_Hemi_dR_smu450_neu200","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu450_neu250 = new TH1F("hSim_Hemi_dR_smu450_neu250","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu450_neu300 = new TH1F("hSim_Hemi_dR_smu450_neu300","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu450_neu350 = new TH1F("hSim_Hemi_dR_smu450_neu350","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu450_neu400 = new TH1F("hSim_Hemi_dR_smu450_neu400","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu450_neu430 = new TH1F("hSim_Hemi_dR_smu450_neu430","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu500_neu180 = new TH1F("hSim_Hemi_dR_smu500_neu180","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu500_neu200 = new TH1F("hSim_Hemi_dR_smu500_neu200","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu500_neu250 = new TH1F("hSim_Hemi_dR_smu500_neu250","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu500_neu300 = new TH1F("hSim_Hemi_dR_smu500_neu300","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu500_neu350 = new TH1F("hSim_Hemi_dR_smu500_neu350","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu500_neu400 = new TH1F("hSim_Hemi_dR_smu500_neu400","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu500_neu450 = new TH1F("hSim_Hemi_dR_smu500_neu450","",33,0.,3.3);
 TH1F* hSim_Hemi_dR_smu500_neu480 = new TH1F("hSim_Hemi_dR_smu500_neu480","",33,0.,3.3);
 
 TH1F* hData_Hemi_Vtx_dist_smu200_neu180 = new TH1F("hData_Hemi_Vtx_dist_smu200_neu180","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu250_neu180 = new TH1F("hData_Hemi_Vtx_dist_smu250_neu180","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu250_neu200 = new TH1F("hData_Hemi_Vtx_dist_smu250_neu200","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu250_neu230 = new TH1F("hData_Hemi_Vtx_dist_smu250_neu230","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu300_neu180 = new TH1F("hData_Hemi_Vtx_dist_smu300_neu180","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu300_neu200 = new TH1F("hData_Hemi_Vtx_dist_smu300_neu200","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu300_neu250 = new TH1F("hData_Hemi_Vtx_dist_smu300_neu250","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu300_neu280 = new TH1F("hData_Hemi_Vtx_dist_smu300_neu280","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu350_neu180 = new TH1F("hData_Hemi_Vtx_dist_smu350_neu180","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu350_neu200 = new TH1F("hData_Hemi_Vtx_dist_smu350_neu200","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu350_neu250 = new TH1F("hData_Hemi_Vtx_dist_smu350_neu250","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu350_neu300 = new TH1F("hData_Hemi_Vtx_dist_smu350_neu300","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu350_neu330 = new TH1F("hData_Hemi_Vtx_dist_smu350_neu330","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu400_neu180 = new TH1F("hData_Hemi_Vtx_dist_smu400_neu180","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu400_neu200 = new TH1F("hData_Hemi_Vtx_dist_smu400_neu200","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu400_neu250 = new TH1F("hData_Hemi_Vtx_dist_smu400_neu250","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu400_neu300 = new TH1F("hData_Hemi_Vtx_dist_smu400_neu300","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu400_neu350 = new TH1F("hData_Hemi_Vtx_dist_smu400_neu350","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu400_neu380 = new TH1F("hData_Hemi_Vtx_dist_smu400_neu380","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu450_neu180 = new TH1F("hData_Hemi_Vtx_dist_smu450_neu180","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu450_neu200 = new TH1F("hData_Hemi_Vtx_dist_smu450_neu200","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu450_neu250 = new TH1F("hData_Hemi_Vtx_dist_smu450_neu250","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu450_neu300 = new TH1F("hData_Hemi_Vtx_dist_smu450_neu300","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu450_neu350 = new TH1F("hData_Hemi_Vtx_dist_smu450_neu350","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu450_neu400 = new TH1F("hData_Hemi_Vtx_dist_smu450_neu400","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu450_neu430 = new TH1F("hData_Hemi_Vtx_dist_smu450_neu430","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu500_neu180 = new TH1F("hData_Hemi_Vtx_dist_smu500_neu180","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu500_neu200 = new TH1F("hData_Hemi_Vtx_dist_smu500_neu200","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu500_neu250 = new TH1F("hData_Hemi_Vtx_dist_smu500_neu250","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu500_neu300 = new TH1F("hData_Hemi_Vtx_dist_smu500_neu300","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu500_neu350 = new TH1F("hData_Hemi_Vtx_dist_smu500_neu350","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu500_neu400 = new TH1F("hData_Hemi_Vtx_dist_smu500_neu400","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu500_neu450 = new TH1F("hData_Hemi_Vtx_dist_smu500_neu450","",200,0.,200.);
 TH1F* hData_Hemi_Vtx_dist_smu500_neu480 = new TH1F("hData_Hemi_Vtx_dist_smu500_neu480","",200,0.,200.);

 TH1F* hSim_Hemi_Vtx_dist_ping_smu200_neu180 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu200_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu250_neu180 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu250_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu250_neu200 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu250_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu250_neu230 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu250_neu230","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu300_neu180 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu300_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu300_neu200 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu300_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu300_neu250 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu300_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu300_neu280 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu300_neu280","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu350_neu180 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu350_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu350_neu200 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu350_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu350_neu250 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu350_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu350_neu300 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu350_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu350_neu330 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu350_neu330","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu400_neu180 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu400_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu400_neu200 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu400_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu400_neu250 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu400_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu400_neu300 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu400_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu400_neu350 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu400_neu350","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu400_neu380 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu400_neu380","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu450_neu180 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu450_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu450_neu200 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu450_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu450_neu250 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu450_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu450_neu300 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu450_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu450_neu350 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu450_neu350","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu450_neu400 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu450_neu400","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu450_neu430 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu450_neu430","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu500_neu180 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu500_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu500_neu200 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu500_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu500_neu250 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu500_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu500_neu300 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu500_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu500_neu350 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu500_neu350","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu500_neu400 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu500_neu400","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu500_neu450 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu500_neu450","",200,0.,200.);
 TH1F* hSim_Hemi_Vtx_dist_ping_smu500_neu480 = new TH1F("hSim_Hemi_Vtx_dist_ping_smu500_neu480","",200,0.,200.);

 TH1F* hData_Hemi_1Vtx_dist_smu200_neu180 = new TH1F("hData_Hemi_1Vtx_dist_smu200_neu180","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu250_neu180 = new TH1F("hData_Hemi_1Vtx_dist_smu250_neu180","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu250_neu200 = new TH1F("hData_Hemi_1Vtx_dist_smu250_neu200","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu250_neu230 = new TH1F("hData_Hemi_1Vtx_dist_smu250_neu230","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu300_neu180 = new TH1F("hData_Hemi_1Vtx_dist_smu300_neu180","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu300_neu200 = new TH1F("hData_Hemi_1Vtx_dist_smu300_neu200","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu300_neu250 = new TH1F("hData_Hemi_1Vtx_dist_smu300_neu250","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu300_neu280 = new TH1F("hData_Hemi_1Vtx_dist_smu300_neu280","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu350_neu180 = new TH1F("hData_Hemi_1Vtx_dist_smu350_neu180","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu350_neu200 = new TH1F("hData_Hemi_1Vtx_dist_smu350_neu200","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu350_neu250 = new TH1F("hData_Hemi_1Vtx_dist_smu350_neu250","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu350_neu300 = new TH1F("hData_Hemi_1Vtx_dist_smu350_neu300","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu350_neu330 = new TH1F("hData_Hemi_1Vtx_dist_smu350_neu330","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu400_neu180 = new TH1F("hData_Hemi_1Vtx_dist_smu400_neu180","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu400_neu200 = new TH1F("hData_Hemi_1Vtx_dist_smu400_neu200","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu400_neu250 = new TH1F("hData_Hemi_1Vtx_dist_smu400_neu250","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu400_neu300 = new TH1F("hData_Hemi_1Vtx_dist_smu400_neu300","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu400_neu350 = new TH1F("hData_Hemi_1Vtx_dist_smu400_neu350","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu400_neu380 = new TH1F("hData_Hemi_1Vtx_dist_smu400_neu380","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu450_neu180 = new TH1F("hData_Hemi_1Vtx_dist_smu450_neu180","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu450_neu200 = new TH1F("hData_Hemi_1Vtx_dist_smu450_neu200","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu450_neu250 = new TH1F("hData_Hemi_1Vtx_dist_smu450_neu250","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu450_neu300 = new TH1F("hData_Hemi_1Vtx_dist_smu450_neu300","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu450_neu350 = new TH1F("hData_Hemi_1Vtx_dist_smu450_neu350","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu450_neu400 = new TH1F("hData_Hemi_1Vtx_dist_smu450_neu400","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu450_neu430 = new TH1F("hData_Hemi_1Vtx_dist_smu450_neu430","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu500_neu180 = new TH1F("hData_Hemi_1Vtx_dist_smu500_neu180","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu500_neu200 = new TH1F("hData_Hemi_1Vtx_dist_smu500_neu200","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu500_neu250 = new TH1F("hData_Hemi_1Vtx_dist_smu500_neu250","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu500_neu300 = new TH1F("hData_Hemi_1Vtx_dist_smu500_neu300","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu500_neu350 = new TH1F("hData_Hemi_1Vtx_dist_smu500_neu350","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu500_neu400 = new TH1F("hData_Hemi_1Vtx_dist_smu500_neu400","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu500_neu450 = new TH1F("hData_Hemi_1Vtx_dist_smu500_neu450","",200,0.,200.);
 TH1F* hData_Hemi_1Vtx_dist_smu500_neu480 = new TH1F("hData_Hemi_1Vtx_dist_smu500_neu480","",200,0.,200.);

 TH1F* hSim_Hemi_1Vtx_dist_ping_smu200_neu180 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu200_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu250_neu180 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu250_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu250_neu200 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu250_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu250_neu230 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu250_neu230","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu300_neu180 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu300_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu300_neu200 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu300_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu300_neu250 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu300_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu300_neu280 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu300_neu280","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu350_neu180 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu350_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu350_neu200 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu350_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu350_neu250 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu350_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu350_neu300 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu350_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu350_neu330 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu350_neu330","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu400_neu180 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu400_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu400_neu200 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu400_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu400_neu250 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu400_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu400_neu300 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu400_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu400_neu350 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu400_neu350","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu400_neu380 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu400_neu380","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu450_neu180 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu450_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu450_neu200 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu450_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu450_neu250 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu450_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu450_neu300 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu450_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu450_neu350 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu450_neu350","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu450_neu400 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu450_neu400","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu450_neu430 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu450_neu430","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu500_neu180 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu500_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu500_neu200 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu500_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu500_neu250 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu500_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu500_neu300 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu500_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu500_neu350 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu500_neu350","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu500_neu400 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu500_neu400","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu500_neu450 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu500_neu450","",200,0.,200.);
 TH1F* hSim_Hemi_1Vtx_dist_ping_smu500_neu480 = new TH1F("hSim_Hemi_1Vtx_dist_ping_smu500_neu480","",200,0.,200.);

 TH1F* hData_Hemi_2Vtx_dist_smu200_neu180 = new TH1F("hData_Hemi_2Vtx_dist_smu200_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu250_neu180 = new TH1F("hData_Hemi_2Vtx_dist_smu250_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu250_neu200 = new TH1F("hData_Hemi_2Vtx_dist_smu250_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu250_neu230 = new TH1F("hData_Hemi_2Vtx_dist_smu250_neu230","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu300_neu180 = new TH1F("hData_Hemi_2Vtx_dist_smu300_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu300_neu200 = new TH1F("hData_Hemi_2Vtx_dist_smu300_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu300_neu250 = new TH1F("hData_Hemi_2Vtx_dist_smu300_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu300_neu280 = new TH1F("hData_Hemi_2Vtx_dist_smu300_neu280","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu350_neu180 = new TH1F("hData_Hemi_2Vtx_dist_smu350_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu350_neu200 = new TH1F("hData_Hemi_2Vtx_dist_smu350_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu350_neu250 = new TH1F("hData_Hemi_2Vtx_dist_smu350_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu350_neu300 = new TH1F("hData_Hemi_2Vtx_dist_smu350_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu350_neu330 = new TH1F("hData_Hemi_2Vtx_dist_smu350_neu330","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu400_neu180 = new TH1F("hData_Hemi_2Vtx_dist_smu400_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu400_neu200 = new TH1F("hData_Hemi_2Vtx_dist_smu400_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu400_neu250 = new TH1F("hData_Hemi_2Vtx_dist_smu400_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu400_neu300 = new TH1F("hData_Hemi_2Vtx_dist_smu400_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu400_neu350 = new TH1F("hData_Hemi_2Vtx_dist_smu400_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu400_neu380 = new TH1F("hData_Hemi_2Vtx_dist_smu400_neu380","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu450_neu180 = new TH1F("hData_Hemi_2Vtx_dist_smu450_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu450_neu200 = new TH1F("hData_Hemi_2Vtx_dist_smu450_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu450_neu250 = new TH1F("hData_Hemi_2Vtx_dist_smu450_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu450_neu300 = new TH1F("hData_Hemi_2Vtx_dist_smu450_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu450_neu350 = new TH1F("hData_Hemi_2Vtx_dist_smu450_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu450_neu400 = new TH1F("hData_Hemi_2Vtx_dist_smu450_neu400","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu450_neu430 = new TH1F("hData_Hemi_2Vtx_dist_smu450_neu430","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu500_neu180 = new TH1F("hData_Hemi_2Vtx_dist_smu500_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu500_neu200 = new TH1F("hData_Hemi_2Vtx_dist_smu500_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu500_neu250 = new TH1F("hData_Hemi_2Vtx_dist_smu500_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu500_neu300 = new TH1F("hData_Hemi_2Vtx_dist_smu500_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu500_neu350 = new TH1F("hData_Hemi_2Vtx_dist_smu500_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu500_neu400 = new TH1F("hData_Hemi_2Vtx_dist_smu500_neu400","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu500_neu450 = new TH1F("hData_Hemi_2Vtx_dist_smu500_neu450","",200,0.,200.);
 TH1F* hData_Hemi_2Vtx_dist_smu500_neu480 = new TH1F("hData_Hemi_2Vtx_dist_smu500_neu480","",200,0.,200.);

 TH1F* hSim_Hemi_2Vtx_dist_ping_smu200_neu180 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu200_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu250_neu180 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu250_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu250_neu200 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu250_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu250_neu230 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu250_neu230","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu300_neu180 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu300_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu300_neu200 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu300_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu300_neu250 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu300_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu300_neu280 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu300_neu280","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu350_neu180 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu350_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu350_neu200 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu350_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu350_neu250 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu350_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu350_neu300 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu350_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu350_neu330 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu350_neu330","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu400_neu180 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu400_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu400_neu200 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu400_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu400_neu250 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu400_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu400_neu300 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu400_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu400_neu350 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu400_neu350","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu400_neu380 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu400_neu380","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu450_neu180 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu450_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu450_neu200 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu450_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu450_neu250 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu450_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu450_neu300 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu450_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu450_neu350 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu450_neu350","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu450_neu400 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu450_neu400","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu450_neu430 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu450_neu430","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu500_neu180 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu500_neu180","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu500_neu200 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu500_neu200","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu500_neu250 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu500_neu250","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu500_neu300 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu500_neu300","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu500_neu350 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu500_neu350","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu500_neu400 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu500_neu400","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu500_neu450 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu500_neu450","",200,0.,200.);
 TH1F* hSim_Hemi_2Vtx_dist_ping_smu500_neu480 = new TH1F("hSim_Hemi_2Vtx_dist_ping_smu500_neu480","",200,0.,200.);

 TH1F* hData_Hemi_2VtxIni_dist_smu200_neu180 = new TH1F("hData_Hemi_2VtxIni_dist_smu200_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu250_neu180 = new TH1F("hData_Hemi_2VtxIni_dist_smu250_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu250_neu200 = new TH1F("hData_Hemi_2VtxIni_dist_smu250_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu250_neu230 = new TH1F("hData_Hemi_2VtxIni_dist_smu250_neu230","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu300_neu180 = new TH1F("hData_Hemi_2VtxIni_dist_smu300_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu300_neu200 = new TH1F("hData_Hemi_2VtxIni_dist_smu300_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu300_neu250 = new TH1F("hData_Hemi_2VtxIni_dist_smu300_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu300_neu280 = new TH1F("hData_Hemi_2VtxIni_dist_smu300_neu280","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu350_neu180 = new TH1F("hData_Hemi_2VtxIni_dist_smu350_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu350_neu200 = new TH1F("hData_Hemi_2VtxIni_dist_smu350_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu350_neu250 = new TH1F("hData_Hemi_2VtxIni_dist_smu350_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu350_neu300 = new TH1F("hData_Hemi_2VtxIni_dist_smu350_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu350_neu330 = new TH1F("hData_Hemi_2VtxIni_dist_smu350_neu330","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu400_neu180 = new TH1F("hData_Hemi_2VtxIni_dist_smu400_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu400_neu200 = new TH1F("hData_Hemi_2VtxIni_dist_smu400_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu400_neu250 = new TH1F("hData_Hemi_2VtxIni_dist_smu400_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu400_neu300 = new TH1F("hData_Hemi_2VtxIni_dist_smu400_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu400_neu350 = new TH1F("hData_Hemi_2VtxIni_dist_smu400_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu400_neu380 = new TH1F("hData_Hemi_2VtxIni_dist_smu400_neu380","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu450_neu180 = new TH1F("hData_Hemi_2VtxIni_dist_smu450_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu450_neu200 = new TH1F("hData_Hemi_2VtxIni_dist_smu450_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu450_neu250 = new TH1F("hData_Hemi_2VtxIni_dist_smu450_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu450_neu300 = new TH1F("hData_Hemi_2VtxIni_dist_smu450_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu450_neu350 = new TH1F("hData_Hemi_2VtxIni_dist_smu450_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu450_neu400 = new TH1F("hData_Hemi_2VtxIni_dist_smu450_neu400","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu450_neu430 = new TH1F("hData_Hemi_2VtxIni_dist_smu450_neu430","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu500_neu180 = new TH1F("hData_Hemi_2VtxIni_dist_smu500_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu500_neu200 = new TH1F("hData_Hemi_2VtxIni_dist_smu500_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu500_neu250 = new TH1F("hData_Hemi_2VtxIni_dist_smu500_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu500_neu300 = new TH1F("hData_Hemi_2VtxIni_dist_smu500_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu500_neu350 = new TH1F("hData_Hemi_2VtxIni_dist_smu500_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu500_neu400 = new TH1F("hData_Hemi_2VtxIni_dist_smu500_neu400","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu500_neu450 = new TH1F("hData_Hemi_2VtxIni_dist_smu500_neu450","",200,0.,200.);
 TH1F* hData_Hemi_2VtxIni_dist_smu500_neu480 = new TH1F("hData_Hemi_2VtxIni_dist_smu500_neu480","",200,0.,200.);

 TH1F* hData_Hemi_2VtxClose_dist_smu200_neu180 = new TH1F("hData_Hemi_2VtxClose_dist_smu200_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu250_neu180 = new TH1F("hData_Hemi_2VtxClose_dist_smu250_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu250_neu200 = new TH1F("hData_Hemi_2VtxClose_dist_smu250_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu250_neu230 = new TH1F("hData_Hemi_2VtxClose_dist_smu250_neu230","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu300_neu180 = new TH1F("hData_Hemi_2VtxClose_dist_smu300_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu300_neu200 = new TH1F("hData_Hemi_2VtxClose_dist_smu300_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu300_neu250 = new TH1F("hData_Hemi_2VtxClose_dist_smu300_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu300_neu280 = new TH1F("hData_Hemi_2VtxClose_dist_smu300_neu280","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu350_neu180 = new TH1F("hData_Hemi_2VtxClose_dist_smu350_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu350_neu200 = new TH1F("hData_Hemi_2VtxClose_dist_smu350_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu350_neu250 = new TH1F("hData_Hemi_2VtxClose_dist_smu350_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu350_neu300 = new TH1F("hData_Hemi_2VtxClose_dist_smu350_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu350_neu330 = new TH1F("hData_Hemi_2VtxClose_dist_smu350_neu330","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu400_neu180 = new TH1F("hData_Hemi_2VtxClose_dist_smu400_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu400_neu200 = new TH1F("hData_Hemi_2VtxClose_dist_smu400_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu400_neu250 = new TH1F("hData_Hemi_2VtxClose_dist_smu400_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu400_neu300 = new TH1F("hData_Hemi_2VtxClose_dist_smu400_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu400_neu350 = new TH1F("hData_Hemi_2VtxClose_dist_smu400_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu400_neu380 = new TH1F("hData_Hemi_2VtxClose_dist_smu400_neu380","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu450_neu180 = new TH1F("hData_Hemi_2VtxClose_dist_smu450_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu450_neu200 = new TH1F("hData_Hemi_2VtxClose_dist_smu450_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu450_neu250 = new TH1F("hData_Hemi_2VtxClose_dist_smu450_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu450_neu300 = new TH1F("hData_Hemi_2VtxClose_dist_smu450_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu450_neu350 = new TH1F("hData_Hemi_2VtxClose_dist_smu450_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu450_neu400 = new TH1F("hData_Hemi_2VtxClose_dist_smu450_neu400","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu450_neu430 = new TH1F("hData_Hemi_2VtxClose_dist_smu450_neu430","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu500_neu180 = new TH1F("hData_Hemi_2VtxClose_dist_smu500_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu500_neu200 = new TH1F("hData_Hemi_2VtxClose_dist_smu500_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu500_neu250 = new TH1F("hData_Hemi_2VtxClose_dist_smu500_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu500_neu300 = new TH1F("hData_Hemi_2VtxClose_dist_smu500_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu500_neu350 = new TH1F("hData_Hemi_2VtxClose_dist_smu500_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu500_neu400 = new TH1F("hData_Hemi_2VtxClose_dist_smu500_neu400","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu500_neu450 = new TH1F("hData_Hemi_2VtxClose_dist_smu500_neu450","",200,0.,200.);
 TH1F* hData_Hemi_2VtxClose_dist_smu500_neu480 = new TH1F("hData_Hemi_2VtxClose_dist_smu500_neu480","",200,0.,200.);

 TH1F* hData_Hemi_2VtxMerge_dist_smu200_neu180 = new TH1F("hData_Hemi_2VtxMerge_dist_smu200_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu250_neu180 = new TH1F("hData_Hemi_2VtxMerge_dist_smu250_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu250_neu200 = new TH1F("hData_Hemi_2VtxMerge_dist_smu250_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu250_neu230 = new TH1F("hData_Hemi_2VtxMerge_dist_smu250_neu230","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu300_neu180 = new TH1F("hData_Hemi_2VtxMerge_dist_smu300_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu300_neu200 = new TH1F("hData_Hemi_2VtxMerge_dist_smu300_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu300_neu250 = new TH1F("hData_Hemi_2VtxMerge_dist_smu300_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu300_neu280 = new TH1F("hData_Hemi_2VtxMerge_dist_smu300_neu280","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu350_neu180 = new TH1F("hData_Hemi_2VtxMerge_dist_smu350_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu350_neu200 = new TH1F("hData_Hemi_2VtxMerge_dist_smu350_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu350_neu250 = new TH1F("hData_Hemi_2VtxMerge_dist_smu350_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu350_neu300 = new TH1F("hData_Hemi_2VtxMerge_dist_smu350_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu350_neu330 = new TH1F("hData_Hemi_2VtxMerge_dist_smu350_neu330","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu400_neu180 = new TH1F("hData_Hemi_2VtxMerge_dist_smu400_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu400_neu200 = new TH1F("hData_Hemi_2VtxMerge_dist_smu400_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu400_neu250 = new TH1F("hData_Hemi_2VtxMerge_dist_smu400_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu400_neu300 = new TH1F("hData_Hemi_2VtxMerge_dist_smu400_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu400_neu350 = new TH1F("hData_Hemi_2VtxMerge_dist_smu400_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu400_neu380 = new TH1F("hData_Hemi_2VtxMerge_dist_smu400_neu380","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu450_neu180 = new TH1F("hData_Hemi_2VtxMerge_dist_smu450_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu450_neu200 = new TH1F("hData_Hemi_2VtxMerge_dist_smu450_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu450_neu250 = new TH1F("hData_Hemi_2VtxMerge_dist_smu450_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu450_neu300 = new TH1F("hData_Hemi_2VtxMerge_dist_smu450_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu450_neu350 = new TH1F("hData_Hemi_2VtxMerge_dist_smu450_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu450_neu400 = new TH1F("hData_Hemi_2VtxMerge_dist_smu450_neu400","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu450_neu430 = new TH1F("hData_Hemi_2VtxMerge_dist_smu450_neu430","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu500_neu180 = new TH1F("hData_Hemi_2VtxMerge_dist_smu500_neu180","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu500_neu200 = new TH1F("hData_Hemi_2VtxMerge_dist_smu500_neu200","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu500_neu250 = new TH1F("hData_Hemi_2VtxMerge_dist_smu500_neu250","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu500_neu300 = new TH1F("hData_Hemi_2VtxMerge_dist_smu500_neu300","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu500_neu350 = new TH1F("hData_Hemi_2VtxMerge_dist_smu500_neu350","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu500_neu400 = new TH1F("hData_Hemi_2VtxMerge_dist_smu500_neu400","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu500_neu450 = new TH1F("hData_Hemi_2VtxMerge_dist_smu500_neu450","",200,0.,200.);
 TH1F* hData_Hemi_2VtxMerge_dist_smu500_neu480 = new TH1F("hData_Hemi_2VtxMerge_dist_smu500_neu480","",200,0.,200.);

 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist001 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist001","",100,0.,0.1);
 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist003 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist003","",100,0.,0.1);
 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist010 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist010","",100,0.,0.2);
 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist030 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist030","",100,0.,0.4);
 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist100 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist100","",100,0.,0.5);
 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist200 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist200","",100,0.,1.0);
 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist300 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist300","",100,0.,4.0);
 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist400 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist400","",100,0.,8.0);
 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist500 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist500","",100,0.,10.);
 TH1F* hSim_Hemi_dSV_step12_etaLT15_dist999 = new TH1F("hSim_Hemi_dSV_step12_etaLT15_dist999","",100,0.,20.);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist001 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist001","",100,0.,0.1);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist003 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist003","",100,0.,0.1);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist010 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist010","",100,0.,0.2);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist030 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist030","",100,0.,0.4);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist100 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist100","",100,0.,0.5);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist200 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist200","",100,0.,1.0);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist300 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist300","",100,0.,4.0);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist400 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist400","",100,0.,8.0);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist500 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist500","",100,0.,10.);
 TH1F* hSim_Hemi_dSV_step12_etaGT15_dist999 = new TH1F("hSim_Hemi_dSV_step12_etaGT15_dist999","",100,0.,20.);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist001 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist001","",100,0.,0.4);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist003 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist003","",100,0.,0.4);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist010 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist010","",100,0.,2.0);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist030 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist030","",100,0.,6.0);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist100 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist100","",100,0.,10.);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist200 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist200","",100,0.,20.);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist300 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist300","",100,0.,20.);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist400 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist400","",100,0.,30.);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist500 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist500","",100,0.,30.);
 TH1F* hSim_Hemi_dSV_step34_etaLT15_dist999 = new TH1F("hSim_Hemi_dSV_step34_etaLT15_dist999","",100,0.,60.);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist001 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist001","",100,0.,0.4);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist003 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist003","",100,0.,0.4);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist010 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist010","",100,0.,2.0);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist030 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist030","",100,0.,6.0);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist100 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist100","",100,0.,10.);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist200 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist200","",100,0.,20.);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist300 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist300","",100,0.,20.);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist400 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist400","",100,0.,30.);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist500 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist500","",100,0.,30.);
 TH1F* hSim_Hemi_dSV_step34_etaGT15_dist999 = new TH1F("hSim_Hemi_dSV_step34_etaGT15_dist999","",100,0.,60.);

 TH1F* hSim_LLP_dist = new TH1F("hSim_LLP_dist","",25,0.,50.);
 TH1F* hSim_Hemi_Vtx_dist_ping = new TH1F("hSim_Hemi_Vtx_dist_ping","",25,0.,50.);
 TH1F* hData_Hemi_Vtx_dist = new TH1F("hData_Hemi_Vtx_dist","",25,0.,50.);
 TH1F* hSim_LLP_dist_big = new TH1F("hSim_LLP_dist_big","",50,0.,100.);
 TH1F* hSim_Hemi_Vtx_dist_ping_big = new TH1F("hSim_Hemi_Vtx_dist_ping_big","",50,0.,100.);
 TH1F* hData_Hemi_Vtx_dist_big = new TH1F("hData_Hemi_Vtx_dist_big","",50,0.,100.);

///////////////////////////////////////////////////////////////////

 if (fChain == 0) return;

 Long64_t nentries = fChain->GetEntriesFast();
 Long64_t nentries2 = fChain->GetEntries();
 Long64_t nentries3 = fChain->GetEntriesFast();
 std::cout << "Total Entries : " << nentries << std::endl;
 std::cout << "Total Entries2 : " << nentries2 << std::endl;
 std::cout << "Total Entries3 : " << nentries3 << std::endl;
 Long64_t nentries4 = minitree_njetNOmu->size();
 std::cout << "Total Entries4 : " << nentries4 << std::endl;
  
 Long64_t nbytes = 0, nb = 0;

 int nRecoTracks=0; 
 int nJets=0; 

 int allevents = 0;
 int itest = 0;

////////////////
// Event loop //
////////////////

 for (Long64_t jentry=0; jentry<nentries; jentry++) {
   Long64_t ientry = LoadTree(jentry);
//$$   cout << " jentry ientry " << jentry << " " << ientry << endl;
 if (ientry < 0) break;
   nb = fChain->GetEntry(jentry);   nbytes += nb;
   
   allevents++;
   if ( allevents%1000 == 0 ) std::cout << "events : " << allevents << std::endl;
   
//$$
//   if ( !minitree_Filter ) continue;
//   if ( minitree_njetNOmu->at(0) < 1 ) continue; 
//   if ( abs(minitree_Hemi_eta->at(0)) > 2.4 || abs(minitree_Hemi_eta->at(1)) > 2.4 ) continue;
//   if ( minitree_Hemi_pt->at(0) < 20. || minitree_Hemi_pt->at(1) < 20. ) continue;
//   if ( minitree_Hemi_SecVtx_x->size() < 2 ) continue;
//   if ( minitree_event_MergedVtx_Vtx_dd->at(0) > 0.1 ) continue;
//   if ( minitree_event_MergedVtx_Vtx_step->at(0) != 2 ) continue;
//   if ( minitree_Hemi_Vtx_NChi2->at(0) <= 0 || minitree_Hemi_Vtx_NChi2->at(0) >= 10 ) continue;
//   if ( minitree_Hemi_Vtx_NChi2->at(1) <= 0 || minitree_Hemi_Vtx_NChi2->at(1) >= 10 ) continue;
//   if ( minitree_Hemi_SecVtx_NChi2->at(0) <= 0 || minitree_Hemi_SecVtx_NChi2->at(0) >= 10 ) continue;
//   if ( minitree_Hemi_SecVtx_NChi2->at(1) <= 0 || minitree_Hemi_SecVtx_NChi2->at(1) >= 10 ) continue;
//$$ 

//   itest++;
//   if ( itest > 10 ) break;
//  std::cout << " " << std::endl;
//  std::cout << " Event " << eventNumber << std::endl;

//$$ 
//  std::cout << " GenPV x y z " << minitree_GenPVx << "  " << minitree_GenPVy << "  " << minitree_GenPVz << std::endl;
//  std::cout << "    PV x y z " << minitree_PV_x   << "  " << minitree_PV_y   << "  " << minitree_PV_z << std::endl;
//  std::cout << std::endl;
// 
//  std::cout << " LLP0  eta phi " << minitree_LLP_eta->at(0) << "  " << minitree_LLP_phi->at(0) << "  x y z " << minitree_LLP_x->at(0) << "  " << minitree_LLP_y->at(0) << "  " << minitree_LLP_z->at(0) << std::endl;
//  std::cout << " LLP1  eta phi " << minitree_LLP_eta->at(1) << "  " << minitree_LLP_phi->at(1) << "  x y z " << minitree_LLP_x->at(1) << "  " << minitree_LLP_y->at(1) << "  " << minitree_LLP_z->at(1) << std::endl;
//  std::cout << std::endl;
//  
//  std::cout << " Hemi0 LLP" << minitree_Hemi_LLP->at(0)-1 << "  eta phi " << minitree_Hemi_eta->at(0) << "  " << minitree_Hemi_phi->at(0) << "  x y z " << minitree_Hemi_Vtx_x->at(0) << "  " << minitree_Hemi_Vtx_y->at(0) << "  " << minitree_Hemi_Vtx_z->at(0) << "  ping test " << minitree_Hemi_LLP_ping->at(0) << " " << minitree_event_LLP_ping_test->at(0) << "  dR " << minitree_Hemi_Vtx_dR->at(0) << std::endl;
//  std::cout << " Hemi1 LLP" << minitree_Hemi_LLP->at(1)-1 << "  eta phi " << minitree_Hemi_eta->at(1) << "  " << minitree_Hemi_phi->at(1) << "  x y z " << minitree_Hemi_Vtx_x->at(1) << "  " << minitree_Hemi_Vtx_y->at(1) << "  " << minitree_Hemi_Vtx_z->at(1) << "  ping test " << minitree_Hemi_LLP_ping->at(1) << " " << minitree_event_LLP_ping_test->at(0) << "  dR " << minitree_Hemi_Vtx_dR->at(1) << std::endl;
//  std::cout << std::endl;
//  
//  std::cout << " Hemi" << minitree_Hemi_SecVtx->at(0)-1 << " LLP" << minitree_Hemi_SecLLP->at(0)-1 << "  x y z " << minitree_Hemi_SecVtx_x->at(0) << "  " << minitree_Hemi_SecVtx_y->at(0) << "  " << minitree_Hemi_SecVtx_z->at(0) << "  ping test " << minitree_Hemi_SecLLP_ping->at(0) << " " << minitree_event_LLP_Newping_test->at(0) << "  dR " << minitree_Hemi_SecVtx_dR->at(0) << std::endl;
//  std::cout << " Hemi" << minitree_Hemi_SecVtx->at(1)-1 << " LLP" << minitree_Hemi_SecLLP->at(1)-1 << "  x y z " << minitree_Hemi_SecVtx_x->at(1) << "  " << minitree_Hemi_SecVtx_y->at(1) << "  " << minitree_Hemi_SecVtx_z->at(1) << "  ping test " << minitree_Hemi_SecLLP_ping->at(1) << " " << minitree_event_LLP_Newping_test->at(0) << "  dR " << minitree_Hemi_SecVtx_dR->at(1) << std::endl;
//  std::cout << std::endl;
//$$ 


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////                                                               ////////
////////                          GENERATION                           ////////
////////                                                               ////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   int Msmu = minitree_smu_mass->at(0);
   int Mneu = minitree_neu_mass->at(0);

   if ( abs(Msmu - 200) < 2 )      Msmu = 200;
   else if ( abs(Msmu - 250) < 2 ) Msmu = 250;
   else if ( abs(Msmu - 300) < 2 ) Msmu = 300;
   else if ( abs(Msmu - 350) < 2 ) Msmu = 350;
   else if ( abs(Msmu - 400) < 2 ) Msmu = 400;
   else if ( abs(Msmu - 450) < 2 ) Msmu = 450;
   else if ( abs(Msmu - 500) < 2 ) Msmu = 500;
   else std::cout << " !!! smu mass out of range !!! " << Msmu << std::endl;
   if ( abs(Mneu - 180) < 2 )      Mneu = 180;
   else if ( abs(Mneu - 200) < 2 ) Mneu = 200;
   else if ( abs(Mneu - 230) < 2 ) Mneu = 230;
   else if ( abs(Mneu - 250) < 2 ) Mneu = 250;
   else if ( abs(Mneu - 280) < 2 ) Mneu = 280;
   else if ( abs(Mneu - 300) < 2 ) Mneu = 300;
   else if ( abs(Mneu - 330) < 2 ) Mneu = 330;
   else if ( abs(Mneu - 350) < 2 ) Mneu = 350;
   else if ( abs(Mneu - 380) < 2 ) Mneu = 380;
   else if ( abs(Mneu - 400) < 2 ) Mneu = 400;
   else if ( abs(Mneu - 430) < 2 ) Mneu = 430;
   else if ( abs(Mneu - 450) < 2 ) Mneu = 450;
   else if ( abs(Mneu - 480) < 2 ) Mneu = 480;
   else std::cout << " !!! neu mass out of range !!! " << Mneu << std::endl;

   hGen_Msmu->Fill( minitree_smu_mass->at(0) );
   hGen_Mneu->Fill( minitree_neu_mass->at(0) );

//  std::cout << std::endl; 

//    int ngenpart =  minitree_genParticle_pt->size();
//    for (int i=0; i<ngenpart; i++)    // Loop on GenParticle
//    {
//      float pdgId = minitree_genParticle_pdgId->at(i); 
//      float mother_pdgId = minitree_genParticle_mother_pdgId->at(i); 
//      float ct0 = minitree_genParticle_ct0->at(i);
// 
// // top quark from neutralino
//      if ( abs(pdgId) == 6 && abs(mother_pdgId) == 1000023 ) {
//        hGen_ct0->Fill( ct0 );
//        if ( Msmu == 200 && Mneu == 180 ) hGen_ct0_smu200_neu180->Fill( ct0 );
//        if ( Msmu == 250 && Mneu == 180 ) hGen_ct0_smu250_neu180->Fill( ct0 );
//        if ( Msmu == 250 && Mneu == 200 ) hGen_ct0_smu250_neu200->Fill( ct0 );
//        if ( Msmu == 250 && Mneu == 230 ) hGen_ct0_smu250_neu230->Fill( ct0 );
//        if ( Msmu == 300 && Mneu == 180 ) hGen_ct0_smu300_neu180->Fill( ct0 );
//        if ( Msmu == 300 && Mneu == 200 ) hGen_ct0_smu300_neu200->Fill( ct0 );
//        if ( Msmu == 300 && Mneu == 250 ) hGen_ct0_smu300_neu250->Fill( ct0 );
//        if ( Msmu == 300 && Mneu == 280 ) hGen_ct0_smu300_neu280->Fill( ct0 );
//        if ( Msmu == 350 && Mneu == 180 ) hGen_ct0_smu350_neu180->Fill( ct0 );
//        if ( Msmu == 350 && Mneu == 200 ) hGen_ct0_smu350_neu200->Fill( ct0 );
//        if ( Msmu == 350 && Mneu == 250 ) hGen_ct0_smu350_neu250->Fill( ct0 );
//        if ( Msmu == 350 && Mneu == 300 ) hGen_ct0_smu350_neu300->Fill( ct0 );
//        if ( Msmu == 350 && Mneu == 330 ) hGen_ct0_smu350_neu330->Fill( ct0 );
//        if ( Msmu == 400 && Mneu == 180 ) hGen_ct0_smu400_neu180->Fill( ct0 );
//        if ( Msmu == 400 && Mneu == 200 ) hGen_ct0_smu400_neu200->Fill( ct0 );
//        if ( Msmu == 400 && Mneu == 250 ) hGen_ct0_smu400_neu250->Fill( ct0 );
//        if ( Msmu == 400 && Mneu == 300 ) hGen_ct0_smu400_neu300->Fill( ct0 );
//        if ( Msmu == 400 && Mneu == 350 ) hGen_ct0_smu400_neu350->Fill( ct0 );
//        if ( Msmu == 400 && Mneu == 380 ) hGen_ct0_smu400_neu380->Fill( ct0 );
//        if ( Msmu == 450 && Mneu == 180 ) hGen_ct0_smu450_neu180->Fill( ct0 );
//        if ( Msmu == 450 && Mneu == 200 ) hGen_ct0_smu450_neu200->Fill( ct0 );
//        if ( Msmu == 450 && Mneu == 250 ) hGen_ct0_smu450_neu250->Fill( ct0 );
//        if ( Msmu == 450 && Mneu == 300 ) hGen_ct0_smu450_neu300->Fill( ct0 );
//        if ( Msmu == 450 && Mneu == 350 ) hGen_ct0_smu450_neu350->Fill( ct0 );
//        if ( Msmu == 450 && Mneu == 400 ) hGen_ct0_smu450_neu400->Fill( ct0 );
//        if ( Msmu == 450 && Mneu == 430 ) hGen_ct0_smu450_neu430->Fill( ct0 );
//        if ( Msmu == 500 && Mneu == 180 ) hGen_ct0_smu500_neu180->Fill( ct0 );
//        if ( Msmu == 500 && Mneu == 200 ) hGen_ct0_smu500_neu200->Fill( ct0 );
//        if ( Msmu == 500 && Mneu == 250 ) hGen_ct0_smu500_neu250->Fill( ct0 );
//        if ( Msmu == 500 && Mneu == 300 ) hGen_ct0_smu500_neu300->Fill( ct0 );
//        if ( Msmu == 500 && Mneu == 350 ) hGen_ct0_smu500_neu350->Fill( ct0 );
//        if ( Msmu == 500 && Mneu == 400 ) hGen_ct0_smu500_neu400->Fill( ct0 );
//        if ( Msmu == 500 && Mneu == 450 ) hGen_ct0_smu500_neu450->Fill( ct0 );
//        if ( Msmu == 500 && Mneu == 480 ) hGen_ct0_smu500_neu480->Fill( ct0 );
//      }
// 
//    }    // End Loop on GenParticle
   
   
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////                                                               ////////
////////                         PRESELECTION                          ////////
////////                                                               ////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

   bool Filter = minitree_Filter; // Tight ID and Tight Iso for both muons
//$$
//    Filter = false;
// 
// // eventually reselect the reco muons 
//    TLorentzVector v, v1, v2;
//    float Mmumu = 0.;
//    float mu_mass = 0.1057;
//    int nmu = minitree_muon_pt->size();
// 
//    if ( nmu >= 2 ) {
//    for (int i=0; i<nmu-1; i++) {
//    if ( !minitree_muon_isGlobal->at(i) )      continue;
//    if ( minitree_muon_pt->at(i) < 10. )       continue;
//    if ( abs(minitree_muon_eta->at(i)) > 2.4 ) continue;
//    if ( !minitree_muon_isTight->at(i) )       continue;
//    if ( abs(minitree_muon_dxy->at(i)) > 0.1 ) continue;
//    if ( abs(minitree_muon_dz->at(i)) > 0.2 )  continue;
//    if ( !minitree_muon_MiniIsoTight->at(i) )  continue;
//      v1.SetPtEtaPhiM(minitree_muon_pt->at(i),minitree_muon_eta->at(i),minitree_muon_phi->at(i),mu_mass);
// 
//      for (int j=i+1; j<nmu; j++) {
//      if ( !minitree_muon_isGlobal->at(j) )      continue;
//      if ( minitree_muon_pt->at(j) < 10. )       continue;
//      if ( abs(minitree_muon_eta->at(j)) > 2.4 ) continue;
//      if ( !minitree_muon_isTight->at(j) )       continue;
//      if ( abs(minitree_muon_dxy->at(j)) > 0.1 ) continue;
//      if ( abs(minitree_muon_dz->at(j)) > 0.2 )  continue;
//      if ( !minitree_muon_MiniIsoTight->at(j) )  continue;
//      if ( minitree_muon_charge->at(i) == minitree_muon_charge->at(j) )     continue;
//      if ( minitree_muon_pt->at(i) < 25. && minitree_muon_pt->at(j) < 25. ) continue;
//        v2.SetPtEtaPhiM(minitree_muon_pt->at(j),minitree_muon_eta->at(j),minitree_muon_phi->at(j),mu_mass);
//        v = v1 + v2;
//        if ( v.Mag() > Mmumu ) Mmumu = v.Mag();
//      }
//    }
//    } // >= 2 muons
//    if ( (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v || HLT_IsoMu24_v)
//         && Mmumu > 10. ) Filter = true;
//$$
   if ( minitree_njetNOmu->at(0) < 1 ) Filter = false; 
   
   if ( Msmu == 200 && Mneu == 180 ) hSim_Filter_smu200_neu180->Fill( Filter );
   if ( Msmu == 250 && Mneu == 180 ) hSim_Filter_smu250_neu180->Fill( Filter );
   if ( Msmu == 250 && Mneu == 200 ) hSim_Filter_smu250_neu200->Fill( Filter );
   if ( Msmu == 250 && Mneu == 230 ) hSim_Filter_smu250_neu230->Fill( Filter );
   if ( Msmu == 300 && Mneu == 180 ) hSim_Filter_smu300_neu180->Fill( Filter );
   if ( Msmu == 300 && Mneu == 200 ) hSim_Filter_smu300_neu200->Fill( Filter );
   if ( Msmu == 300 && Mneu == 250 ) hSim_Filter_smu300_neu250->Fill( Filter );
   if ( Msmu == 300 && Mneu == 280 ) hSim_Filter_smu300_neu280->Fill( Filter );
   if ( Msmu == 350 && Mneu == 180 ) hSim_Filter_smu350_neu180->Fill( Filter );
   if ( Msmu == 350 && Mneu == 200 ) hSim_Filter_smu350_neu200->Fill( Filter );
   if ( Msmu == 350 && Mneu == 250 ) hSim_Filter_smu350_neu250->Fill( Filter );
   if ( Msmu == 350 && Mneu == 300 ) hSim_Filter_smu350_neu300->Fill( Filter );
   if ( Msmu == 350 && Mneu == 330 ) hSim_Filter_smu350_neu330->Fill( Filter );
   if ( Msmu == 400 && Mneu == 180 ) hSim_Filter_smu400_neu180->Fill( Filter );
   if ( Msmu == 400 && Mneu == 200 ) hSim_Filter_smu400_neu200->Fill( Filter );
   if ( Msmu == 400 && Mneu == 250 ) hSim_Filter_smu400_neu250->Fill( Filter );
   if ( Msmu == 400 && Mneu == 300 ) hSim_Filter_smu400_neu300->Fill( Filter );
   if ( Msmu == 400 && Mneu == 350 ) hSim_Filter_smu400_neu350->Fill( Filter );
   if ( Msmu == 400 && Mneu == 380 ) hSim_Filter_smu400_neu380->Fill( Filter );
   if ( Msmu == 450 && Mneu == 180 ) hSim_Filter_smu450_neu180->Fill( Filter );
   if ( Msmu == 450 && Mneu == 200 ) hSim_Filter_smu450_neu200->Fill( Filter );
   if ( Msmu == 450 && Mneu == 250 ) hSim_Filter_smu450_neu250->Fill( Filter );
   if ( Msmu == 450 && Mneu == 300 ) hSim_Filter_smu450_neu300->Fill( Filter );
   if ( Msmu == 450 && Mneu == 350 ) hSim_Filter_smu450_neu350->Fill( Filter );
   if ( Msmu == 450 && Mneu == 400 ) hSim_Filter_smu450_neu400->Fill( Filter );
   if ( Msmu == 450 && Mneu == 430 ) hSim_Filter_smu450_neu430->Fill( Filter );
   if ( Msmu == 500 && Mneu == 180 ) hSim_Filter_smu500_neu180->Fill( Filter );
   if ( Msmu == 500 && Mneu == 200 ) hSim_Filter_smu500_neu200->Fill( Filter );
   if ( Msmu == 500 && Mneu == 250 ) hSim_Filter_smu500_neu250->Fill( Filter );
   if ( Msmu == 500 && Mneu == 300 ) hSim_Filter_smu500_neu300->Fill( Filter );
   if ( Msmu == 500 && Mneu == 350 ) hSim_Filter_smu500_neu350->Fill( Filter );
   if ( Msmu == 500 && Mneu == 400 ) hSim_Filter_smu500_neu400->Fill( Filter );
   if ( Msmu == 500 && Mneu == 450 ) hSim_Filter_smu500_neu450->Fill( Filter );
   if ( Msmu == 500 && Mneu == 480 ) hSim_Filter_smu500_neu480->Fill( Filter );

//$$
 if ( !Filter ) continue; 
 if ( minitree_Hemi_pt->size() != 2 ) continue;
//$$

   if ( minitree_nLLP->at(0) > 0 ) hSim_LLP_dist->Fill( minitree_LLP_dist->at(0) );
   if ( minitree_nLLP->at(0) > 1 ) hSim_LLP_dist->Fill( minitree_LLP_dist->at(1) );
   if ( minitree_nLLP->at(0) > 0 ) hSim_LLP_dist_big->Fill( minitree_LLP_dist->at(0) );
   if ( minitree_nLLP->at(0) > 1 ) hSim_LLP_dist_big->Fill( minitree_LLP_dist->at(1) );

//$$
 if ( abs(minitree_Hemi_eta->at(0)) > 2.4 || abs(minitree_Hemi_eta->at(1)) > 2.4 ) continue;
 if ( minitree_Hemi_pt->at(0) < 20. || minitree_Hemi_pt->at(1) < 20. ) continue;
//$$

   bool isHemiVtx1 = false, isHemiVtx2 = false, ping;
   int nVtx = 0, nVtxIni = 0, step;
   float VtxMass = 0., dR, dist, NChi2;

   bool ping0 = minitree_Hemi_LLP_ping->at(0);
   bool ping1 = minitree_Hemi_LLP_ping->at(1);
   int Vtx_step0 = minitree_Hemi_Vtx_step->at(0);
   int Vtx_step1 = minitree_Hemi_Vtx_step->at(1);
   float Vtx_NChi0 = minitree_Hemi_Vtx_NChi2->at(0);
   float Vtx_NChi1 = minitree_Hemi_Vtx_NChi2->at(1);
   float Vtx_Mass0 = minitree_Hemi_Vtx_Mass->at(0);
   float Vtx_Mass1 = minitree_Hemi_Vtx_Mass->at(1);
   float Vtx_dist0 = minitree_Hemi_Vtx_dist->at(0);
   float Vtx_dist1 = minitree_Hemi_Vtx_dist->at(1);
   float dR0 = minitree_Hemi_LLP_dR->at(0);
   float dR1 = minitree_Hemi_LLP_dR->at(1);

   float posx0 = minitree_Hemi_Vtx_x->at(0);
   float posy0 = minitree_Hemi_Vtx_y->at(0);
   float posz0 = minitree_Hemi_Vtx_z->at(0);
   float posx1 = minitree_Hemi_Vtx_x->at(1);
   float posy1 = minitree_Hemi_Vtx_y->at(1);
   float posz1 = minitree_Hemi_Vtx_z->at(1);
   float r0 = TMath::Sqrt( posx0*posx0 + posy0*posy0 );
   float z0 = TMath::Abs( posz0 );
   float r1 = TMath::Sqrt( posx1*posx1 + posy1*posy1 );
   float z1 = TMath::Abs( posz1 );
   
   float LLPx0 = minitree_LLP_x->at(0);
   float LLPy0 = minitree_LLP_y->at(0);
   float LLPz0 = minitree_LLP_z->at(0);
   float LLPx1 = minitree_LLP_x->at(1);
   float LLPy1 = minitree_LLP_y->at(1);
   float LLPz1 = minitree_LLP_z->at(1);
   float ddLLP = TMath::Sqrt( (LLPx0-LLPx1)*(LLPx0-LLPx1) + (LLPy0-LLPy1)*(LLPy0-LLPy1) + (LLPz0-LLPz1)*(LLPz0-LLPz1) );
//$$
//  if ( ddLLP < 0.4 ) continue;
//$$

   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 1 && Vtx_step0 <= 2 ) {
     isHemiVtx1 = true;
   }
   if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 1 && Vtx_step1 <= 2 ) {
     isHemiVtx2 = true;
   }
   if      ( isHemiVtx1 && isHemiVtx2 ) nVtxIni = 2;
   else if ( isHemiVtx1 || isHemiVtx2 ) nVtxIni = 1;

//$$$$
   bool Merging = true; // false if "no merge" or "no close vtx" output
//$$$$
   if ( Merging && minitree_event_MergedVtx_Vtx_step->size() >= 1 
		&& minitree_event_MergedVtx_Vtx_step->at(0) != 0 ) {
     hData_Hemi_SecVtx_size->Fill( minitree_Hemi_SecVtx->size() );
     ping0 = minitree_Hemi_SecLLP_ping->at(0);
     Vtx_step0 = minitree_Hemi_SecVtx_step->at(0);
     Vtx_NChi0 = minitree_Hemi_SecVtx_NChi2->at(0);
     Vtx_Mass0 = minitree_Hemi_SecVtx_Mass->at(0);
     Vtx_dist0 = minitree_Hemi_SecVtx_dist->at(0);
     posx0 = minitree_Hemi_SecVtx_x->at(0);
     posy0 = minitree_Hemi_SecVtx_y->at(0);
     posz0 = minitree_Hemi_SecVtx_z->at(0);
     Vtx_step1 = 0;
     Vtx_NChi1 = -1.;
     Vtx_Mass1 = 0.;
     Vtx_dist1 = 0.;
     ping1 = false;
     if ( minitree_event_MergedVtx_Vtx_step->at(0) > 0 ) {
       ping1 = minitree_Hemi_SecLLP_ping->at(1);
       Vtx_step1 = minitree_Hemi_SecVtx_step->at(1);
       Vtx_NChi1 = minitree_Hemi_SecVtx_NChi2->at(1);
       Vtx_Mass1 = minitree_Hemi_SecVtx_Mass->at(1);
       Vtx_dist1 = minitree_Hemi_SecVtx_dist->at(1);
       posx1 = minitree_Hemi_SecVtx_x->at(1);
       posy1 = minitree_Hemi_SecVtx_y->at(1);
       posz1 = minitree_Hemi_SecVtx_z->at(1);
     }
   }

   isHemiVtx1 = false;
   isHemiVtx2 = false;
   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 1 && Vtx_step0 <= 2 ) {
     isHemiVtx1 = true;
     VtxMass = Vtx_Mass0;
   }
   if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 1 && Vtx_step1 <= 2 ) {
     isHemiVtx2 = true;
     if ( Vtx_Mass1 > VtxMass ) VtxMass = Vtx_Mass1;
   }
   if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;
   else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;
   
//    if ( minitree_Hemi_SecVtx_x->size() == 2 ) {
//      cout << " nVtx " << nVtx << endl;
//      cout << " Vtx0 step chi dist ping " << Vtx_step0 << " " << Vtx_NChi0 << " " << Vtx_dist0 << " " << ping0 << endl;
//      cout << " Vtx1 step chi dist ping " << Vtx_step1 << " " << Vtx_NChi1 << " " << Vtx_dist1 << " " << ping1 << endl;
//      cout << endl;
//    }
   

///////////////////////
// Delta R between hemisphere axis and closest neutralino

   int nHemi = minitree_Hemi_LLP_dR->size();
   for (int i=0; i<nHemi; i++) {   // Loop on hemispheres

     if ( i == 0 ) {
       dR    = dR0;
       dist  = Vtx_dist0;
       NChi2 = Vtx_NChi0;
       step  = Vtx_step0;
       ping  = ping0;
     }
     else if ( i == 1 ) {
       dR    = dR1;
       dist  = Vtx_dist1;
       NChi2 = Vtx_NChi1;
       step  = Vtx_step1;
       ping  = ping1;
     }
     
     if ( Msmu == 200 && Mneu == 180 ) hSim_Hemi_dR_smu200_neu180->Fill( dR );
     if ( Msmu == 250 && Mneu == 180 ) hSim_Hemi_dR_smu250_neu180->Fill( dR );
     if ( Msmu == 250 && Mneu == 200 ) hSim_Hemi_dR_smu250_neu200->Fill( dR );
     if ( Msmu == 250 && Mneu == 230 ) hSim_Hemi_dR_smu250_neu230->Fill( dR );
     if ( Msmu == 300 && Mneu == 180 ) hSim_Hemi_dR_smu300_neu180->Fill( dR );
     if ( Msmu == 300 && Mneu == 200 ) hSim_Hemi_dR_smu300_neu200->Fill( dR );
     if ( Msmu == 300 && Mneu == 250 ) hSim_Hemi_dR_smu300_neu250->Fill( dR );
     if ( Msmu == 300 && Mneu == 280 ) hSim_Hemi_dR_smu300_neu280->Fill( dR );
     if ( Msmu == 350 && Mneu == 180 ) hSim_Hemi_dR_smu350_neu180->Fill( dR );
     if ( Msmu == 350 && Mneu == 200 ) hSim_Hemi_dR_smu350_neu200->Fill( dR );
     if ( Msmu == 350 && Mneu == 250 ) hSim_Hemi_dR_smu350_neu250->Fill( dR );
     if ( Msmu == 350 && Mneu == 300 ) hSim_Hemi_dR_smu350_neu300->Fill( dR );
     if ( Msmu == 350 && Mneu == 330 ) hSim_Hemi_dR_smu350_neu330->Fill( dR );
     if ( Msmu == 400 && Mneu == 180 ) hSim_Hemi_dR_smu400_neu180->Fill( dR );
     if ( Msmu == 400 && Mneu == 200 ) hSim_Hemi_dR_smu400_neu200->Fill( dR );
     if ( Msmu == 400 && Mneu == 250 ) hSim_Hemi_dR_smu400_neu250->Fill( dR );
     if ( Msmu == 400 && Mneu == 300 ) hSim_Hemi_dR_smu400_neu300->Fill( dR );
     if ( Msmu == 400 && Mneu == 350 ) hSim_Hemi_dR_smu400_neu350->Fill( dR );
     if ( Msmu == 400 && Mneu == 380 ) hSim_Hemi_dR_smu400_neu380->Fill( dR );
     if ( Msmu == 450 && Mneu == 180 ) hSim_Hemi_dR_smu450_neu180->Fill( dR );
     if ( Msmu == 450 && Mneu == 200 ) hSim_Hemi_dR_smu450_neu200->Fill( dR );
     if ( Msmu == 450 && Mneu == 250 ) hSim_Hemi_dR_smu450_neu250->Fill( dR );
     if ( Msmu == 450 && Mneu == 300 ) hSim_Hemi_dR_smu450_neu300->Fill( dR );
     if ( Msmu == 450 && Mneu == 350 ) hSim_Hemi_dR_smu450_neu350->Fill( dR );
     if ( Msmu == 450 && Mneu == 400 ) hSim_Hemi_dR_smu450_neu400->Fill( dR );
     if ( Msmu == 450 && Mneu == 430 ) hSim_Hemi_dR_smu450_neu430->Fill( dR );
     if ( Msmu == 500 && Mneu == 180 ) hSim_Hemi_dR_smu500_neu180->Fill( dR );
     if ( Msmu == 500 && Mneu == 200 ) hSim_Hemi_dR_smu500_neu200->Fill( dR );
     if ( Msmu == 500 && Mneu == 250 ) hSim_Hemi_dR_smu500_neu250->Fill( dR );
     if ( Msmu == 500 && Mneu == 300 ) hSim_Hemi_dR_smu500_neu300->Fill( dR );
     if ( Msmu == 500 && Mneu == 350 ) hSim_Hemi_dR_smu500_neu350->Fill( dR );
     if ( Msmu == 500 && Mneu == 400 ) hSim_Hemi_dR_smu500_neu400->Fill( dR );
     if ( Msmu == 500 && Mneu == 450 ) hSim_Hemi_dR_smu500_neu450->Fill( dR );
     if ( Msmu == 500 && Mneu == 480 ) hSim_Hemi_dR_smu500_neu480->Fill( dR );

   if ( NChi2 > 0. && NChi2 < 10. && (step == 1 || step == 2) ) {
     if ( Msmu == 200 && Mneu == 180 ) hData_Hemi_Vtx_dist_smu200_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 180 ) hData_Hemi_Vtx_dist_smu250_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 200 ) hData_Hemi_Vtx_dist_smu250_neu200->Fill( dist );
     if ( Msmu == 250 && Mneu == 230 ) hData_Hemi_Vtx_dist_smu250_neu230->Fill( dist );
     if ( Msmu == 300 && Mneu == 180 ) hData_Hemi_Vtx_dist_smu300_neu180->Fill( dist );
     if ( Msmu == 300 && Mneu == 200 ) hData_Hemi_Vtx_dist_smu300_neu200->Fill( dist );
     if ( Msmu == 300 && Mneu == 250 ) hData_Hemi_Vtx_dist_smu300_neu250->Fill( dist );
     if ( Msmu == 300 && Mneu == 280 ) hData_Hemi_Vtx_dist_smu300_neu280->Fill( dist );
     if ( Msmu == 350 && Mneu == 180 ) hData_Hemi_Vtx_dist_smu350_neu180->Fill( dist );
     if ( Msmu == 350 && Mneu == 200 ) hData_Hemi_Vtx_dist_smu350_neu200->Fill( dist );
     if ( Msmu == 350 && Mneu == 250 ) hData_Hemi_Vtx_dist_smu350_neu250->Fill( dist );
     if ( Msmu == 350 && Mneu == 300 ) hData_Hemi_Vtx_dist_smu350_neu300->Fill( dist );
     if ( Msmu == 350 && Mneu == 330 ) hData_Hemi_Vtx_dist_smu350_neu330->Fill( dist );
     if ( Msmu == 400 && Mneu == 180 ) hData_Hemi_Vtx_dist_smu400_neu180->Fill( dist );
     if ( Msmu == 400 && Mneu == 200 ) hData_Hemi_Vtx_dist_smu400_neu200->Fill( dist );
     if ( Msmu == 400 && Mneu == 250 ) hData_Hemi_Vtx_dist_smu400_neu250->Fill( dist );
     if ( Msmu == 400 && Mneu == 300 ) hData_Hemi_Vtx_dist_smu400_neu300->Fill( dist );
     if ( Msmu == 400 && Mneu == 350 ) hData_Hemi_Vtx_dist_smu400_neu350->Fill( dist );
     if ( Msmu == 400 && Mneu == 380 ) hData_Hemi_Vtx_dist_smu400_neu380->Fill( dist );
     if ( Msmu == 450 && Mneu == 180 ) hData_Hemi_Vtx_dist_smu450_neu180->Fill( dist );
     if ( Msmu == 450 && Mneu == 200 ) hData_Hemi_Vtx_dist_smu450_neu200->Fill( dist );
     if ( Msmu == 450 && Mneu == 250 ) hData_Hemi_Vtx_dist_smu450_neu250->Fill( dist );
     if ( Msmu == 450 && Mneu == 300 ) hData_Hemi_Vtx_dist_smu450_neu300->Fill( dist );
     if ( Msmu == 450 && Mneu == 350 ) hData_Hemi_Vtx_dist_smu450_neu350->Fill( dist );
     if ( Msmu == 450 && Mneu == 400 ) hData_Hemi_Vtx_dist_smu450_neu400->Fill( dist );
     if ( Msmu == 450 && Mneu == 430 ) hData_Hemi_Vtx_dist_smu450_neu430->Fill( dist );
     if ( Msmu == 500 && Mneu == 180 ) hData_Hemi_Vtx_dist_smu500_neu180->Fill( dist );
     if ( Msmu == 500 && Mneu == 200 ) hData_Hemi_Vtx_dist_smu500_neu200->Fill( dist );
     if ( Msmu == 500 && Mneu == 250 ) hData_Hemi_Vtx_dist_smu500_neu250->Fill( dist );
     if ( Msmu == 500 && Mneu == 300 ) hData_Hemi_Vtx_dist_smu500_neu300->Fill( dist );
     if ( Msmu == 500 && Mneu == 350 ) hData_Hemi_Vtx_dist_smu500_neu350->Fill( dist );
     if ( Msmu == 500 && Mneu == 400 ) hData_Hemi_Vtx_dist_smu500_neu400->Fill( dist );
     if ( Msmu == 500 && Mneu == 450 ) hData_Hemi_Vtx_dist_smu500_neu450->Fill( dist );
     if ( Msmu == 500 && Mneu == 480 ) hData_Hemi_Vtx_dist_smu500_neu480->Fill( dist );
     hData_Hemi_Vtx_dist->Fill( dist );
     hData_Hemi_Vtx_dist_big->Fill( dist );
   }
    
   if ( NChi2 > 0. && NChi2 < 10. && (step == 1 || step == 2) && ping ) {
     if ( Msmu == 200 && Mneu == 180 ) hSim_Hemi_Vtx_dist_ping_smu200_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 180 ) hSim_Hemi_Vtx_dist_ping_smu250_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 200 ) hSim_Hemi_Vtx_dist_ping_smu250_neu200->Fill( dist );
     if ( Msmu == 250 && Mneu == 230 ) hSim_Hemi_Vtx_dist_ping_smu250_neu230->Fill( dist );
     if ( Msmu == 300 && Mneu == 180 ) hSim_Hemi_Vtx_dist_ping_smu300_neu180->Fill( dist );
     if ( Msmu == 300 && Mneu == 200 ) hSim_Hemi_Vtx_dist_ping_smu300_neu200->Fill( dist );
     if ( Msmu == 300 && Mneu == 250 ) hSim_Hemi_Vtx_dist_ping_smu300_neu250->Fill( dist );
     if ( Msmu == 300 && Mneu == 280 ) hSim_Hemi_Vtx_dist_ping_smu300_neu280->Fill( dist );
     if ( Msmu == 350 && Mneu == 180 ) hSim_Hemi_Vtx_dist_ping_smu350_neu180->Fill( dist );
     if ( Msmu == 350 && Mneu == 200 ) hSim_Hemi_Vtx_dist_ping_smu350_neu200->Fill( dist );
     if ( Msmu == 350 && Mneu == 250 ) hSim_Hemi_Vtx_dist_ping_smu350_neu250->Fill( dist );
     if ( Msmu == 350 && Mneu == 300 ) hSim_Hemi_Vtx_dist_ping_smu350_neu300->Fill( dist );
     if ( Msmu == 350 && Mneu == 330 ) hSim_Hemi_Vtx_dist_ping_smu350_neu330->Fill( dist );
     if ( Msmu == 400 && Mneu == 180 ) hSim_Hemi_Vtx_dist_ping_smu400_neu180->Fill( dist );
     if ( Msmu == 400 && Mneu == 200 ) hSim_Hemi_Vtx_dist_ping_smu400_neu200->Fill( dist );
     if ( Msmu == 400 && Mneu == 250 ) hSim_Hemi_Vtx_dist_ping_smu400_neu250->Fill( dist );
     if ( Msmu == 400 && Mneu == 300 ) hSim_Hemi_Vtx_dist_ping_smu400_neu300->Fill( dist );
     if ( Msmu == 400 && Mneu == 350 ) hSim_Hemi_Vtx_dist_ping_smu400_neu350->Fill( dist );
     if ( Msmu == 400 && Mneu == 380 ) hSim_Hemi_Vtx_dist_ping_smu400_neu380->Fill( dist );
     if ( Msmu == 450 && Mneu == 180 ) hSim_Hemi_Vtx_dist_ping_smu450_neu180->Fill( dist );
     if ( Msmu == 450 && Mneu == 200 ) hSim_Hemi_Vtx_dist_ping_smu450_neu200->Fill( dist );
     if ( Msmu == 450 && Mneu == 250 ) hSim_Hemi_Vtx_dist_ping_smu450_neu250->Fill( dist );
     if ( Msmu == 450 && Mneu == 300 ) hSim_Hemi_Vtx_dist_ping_smu450_neu300->Fill( dist );
     if ( Msmu == 450 && Mneu == 350 ) hSim_Hemi_Vtx_dist_ping_smu450_neu350->Fill( dist );
     if ( Msmu == 450 && Mneu == 400 ) hSim_Hemi_Vtx_dist_ping_smu450_neu400->Fill( dist );
     if ( Msmu == 450 && Mneu == 430 ) hSim_Hemi_Vtx_dist_ping_smu450_neu430->Fill( dist );
     if ( Msmu == 500 && Mneu == 180 ) hSim_Hemi_Vtx_dist_ping_smu500_neu180->Fill( dist );
     if ( Msmu == 500 && Mneu == 200 ) hSim_Hemi_Vtx_dist_ping_smu500_neu200->Fill( dist );
     if ( Msmu == 500 && Mneu == 250 ) hSim_Hemi_Vtx_dist_ping_smu500_neu250->Fill( dist );
     if ( Msmu == 500 && Mneu == 300 ) hSim_Hemi_Vtx_dist_ping_smu500_neu300->Fill( dist );
     if ( Msmu == 500 && Mneu == 350 ) hSim_Hemi_Vtx_dist_ping_smu500_neu350->Fill( dist );
     if ( Msmu == 500 && Mneu == 400 ) hSim_Hemi_Vtx_dist_ping_smu500_neu400->Fill( dist );
     if ( Msmu == 500 && Mneu == 450 ) hSim_Hemi_Vtx_dist_ping_smu500_neu450->Fill( dist );
     if ( Msmu == 500 && Mneu == 480 ) hSim_Hemi_Vtx_dist_ping_smu500_neu480->Fill( dist );
     hSim_Hemi_Vtx_dist_ping->Fill( dist );
     hSim_Hemi_Vtx_dist_ping_big->Fill( dist );
   }
     
   if ( nVtx == 1 ) {
   if ( NChi2 > 0. && NChi2 < 10. && (step == 1 || step == 2) ) {
     if ( Msmu == 200 && Mneu == 180 ) hData_Hemi_1Vtx_dist_smu200_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 180 ) hData_Hemi_1Vtx_dist_smu250_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 200 ) hData_Hemi_1Vtx_dist_smu250_neu200->Fill( dist );
     if ( Msmu == 250 && Mneu == 230 ) hData_Hemi_1Vtx_dist_smu250_neu230->Fill( dist );
     if ( Msmu == 300 && Mneu == 180 ) hData_Hemi_1Vtx_dist_smu300_neu180->Fill( dist );
     if ( Msmu == 300 && Mneu == 200 ) hData_Hemi_1Vtx_dist_smu300_neu200->Fill( dist );
     if ( Msmu == 300 && Mneu == 250 ) hData_Hemi_1Vtx_dist_smu300_neu250->Fill( dist );
     if ( Msmu == 300 && Mneu == 280 ) hData_Hemi_1Vtx_dist_smu300_neu280->Fill( dist );
     if ( Msmu == 350 && Mneu == 180 ) hData_Hemi_1Vtx_dist_smu350_neu180->Fill( dist );
     if ( Msmu == 350 && Mneu == 200 ) hData_Hemi_1Vtx_dist_smu350_neu200->Fill( dist );
     if ( Msmu == 350 && Mneu == 250 ) hData_Hemi_1Vtx_dist_smu350_neu250->Fill( dist );
     if ( Msmu == 350 && Mneu == 300 ) hData_Hemi_1Vtx_dist_smu350_neu300->Fill( dist );
     if ( Msmu == 350 && Mneu == 330 ) hData_Hemi_1Vtx_dist_smu350_neu330->Fill( dist );
     if ( Msmu == 400 && Mneu == 180 ) hData_Hemi_1Vtx_dist_smu400_neu180->Fill( dist );
     if ( Msmu == 400 && Mneu == 200 ) hData_Hemi_1Vtx_dist_smu400_neu200->Fill( dist );
     if ( Msmu == 400 && Mneu == 250 ) hData_Hemi_1Vtx_dist_smu400_neu250->Fill( dist );
     if ( Msmu == 400 && Mneu == 300 ) hData_Hemi_1Vtx_dist_smu400_neu300->Fill( dist );
     if ( Msmu == 400 && Mneu == 350 ) hData_Hemi_1Vtx_dist_smu400_neu350->Fill( dist );
     if ( Msmu == 400 && Mneu == 380 ) hData_Hemi_1Vtx_dist_smu400_neu380->Fill( dist );
     if ( Msmu == 450 && Mneu == 180 ) hData_Hemi_1Vtx_dist_smu450_neu180->Fill( dist );
     if ( Msmu == 450 && Mneu == 200 ) hData_Hemi_1Vtx_dist_smu450_neu200->Fill( dist );
     if ( Msmu == 450 && Mneu == 250 ) hData_Hemi_1Vtx_dist_smu450_neu250->Fill( dist );
     if ( Msmu == 450 && Mneu == 300 ) hData_Hemi_1Vtx_dist_smu450_neu300->Fill( dist );
     if ( Msmu == 450 && Mneu == 350 ) hData_Hemi_1Vtx_dist_smu450_neu350->Fill( dist );
     if ( Msmu == 450 && Mneu == 400 ) hData_Hemi_1Vtx_dist_smu450_neu400->Fill( dist );
     if ( Msmu == 450 && Mneu == 430 ) hData_Hemi_1Vtx_dist_smu450_neu430->Fill( dist );
     if ( Msmu == 500 && Mneu == 180 ) hData_Hemi_1Vtx_dist_smu500_neu180->Fill( dist );
     if ( Msmu == 500 && Mneu == 200 ) hData_Hemi_1Vtx_dist_smu500_neu200->Fill( dist );
     if ( Msmu == 500 && Mneu == 250 ) hData_Hemi_1Vtx_dist_smu500_neu250->Fill( dist );
     if ( Msmu == 500 && Mneu == 300 ) hData_Hemi_1Vtx_dist_smu500_neu300->Fill( dist );
     if ( Msmu == 500 && Mneu == 350 ) hData_Hemi_1Vtx_dist_smu500_neu350->Fill( dist );
     if ( Msmu == 500 && Mneu == 400 ) hData_Hemi_1Vtx_dist_smu500_neu400->Fill( dist );
     if ( Msmu == 500 && Mneu == 450 ) hData_Hemi_1Vtx_dist_smu500_neu450->Fill( dist );
     if ( Msmu == 500 && Mneu == 480 ) hData_Hemi_1Vtx_dist_smu500_neu480->Fill( dist );
   }
    
   if ( NChi2 > 0. && NChi2 < 10. && (step == 1 || step == 2) && ping ) {
     if ( Msmu == 200 && Mneu == 180 ) hSim_Hemi_1Vtx_dist_ping_smu200_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 180 ) hSim_Hemi_1Vtx_dist_ping_smu250_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 200 ) hSim_Hemi_1Vtx_dist_ping_smu250_neu200->Fill( dist );
     if ( Msmu == 250 && Mneu == 230 ) hSim_Hemi_1Vtx_dist_ping_smu250_neu230->Fill( dist );
     if ( Msmu == 300 && Mneu == 180 ) hSim_Hemi_1Vtx_dist_ping_smu300_neu180->Fill( dist );
     if ( Msmu == 300 && Mneu == 200 ) hSim_Hemi_1Vtx_dist_ping_smu300_neu200->Fill( dist );
     if ( Msmu == 300 && Mneu == 250 ) hSim_Hemi_1Vtx_dist_ping_smu300_neu250->Fill( dist );
     if ( Msmu == 300 && Mneu == 280 ) hSim_Hemi_1Vtx_dist_ping_smu300_neu280->Fill( dist );
     if ( Msmu == 350 && Mneu == 180 ) hSim_Hemi_1Vtx_dist_ping_smu350_neu180->Fill( dist );
     if ( Msmu == 350 && Mneu == 200 ) hSim_Hemi_1Vtx_dist_ping_smu350_neu200->Fill( dist );
     if ( Msmu == 350 && Mneu == 250 ) hSim_Hemi_1Vtx_dist_ping_smu350_neu250->Fill( dist );
     if ( Msmu == 350 && Mneu == 300 ) hSim_Hemi_1Vtx_dist_ping_smu350_neu300->Fill( dist );
     if ( Msmu == 350 && Mneu == 330 ) hSim_Hemi_1Vtx_dist_ping_smu350_neu330->Fill( dist );
     if ( Msmu == 400 && Mneu == 180 ) hSim_Hemi_1Vtx_dist_ping_smu400_neu180->Fill( dist );
     if ( Msmu == 400 && Mneu == 200 ) hSim_Hemi_1Vtx_dist_ping_smu400_neu200->Fill( dist );
     if ( Msmu == 400 && Mneu == 250 ) hSim_Hemi_1Vtx_dist_ping_smu400_neu250->Fill( dist );
     if ( Msmu == 400 && Mneu == 300 ) hSim_Hemi_1Vtx_dist_ping_smu400_neu300->Fill( dist );
     if ( Msmu == 400 && Mneu == 350 ) hSim_Hemi_1Vtx_dist_ping_smu400_neu350->Fill( dist );
     if ( Msmu == 400 && Mneu == 380 ) hSim_Hemi_1Vtx_dist_ping_smu400_neu380->Fill( dist );
     if ( Msmu == 450 && Mneu == 180 ) hSim_Hemi_1Vtx_dist_ping_smu450_neu180->Fill( dist );
     if ( Msmu == 450 && Mneu == 200 ) hSim_Hemi_1Vtx_dist_ping_smu450_neu200->Fill( dist );
     if ( Msmu == 450 && Mneu == 250 ) hSim_Hemi_1Vtx_dist_ping_smu450_neu250->Fill( dist );
     if ( Msmu == 450 && Mneu == 300 ) hSim_Hemi_1Vtx_dist_ping_smu450_neu300->Fill( dist );
     if ( Msmu == 450 && Mneu == 350 ) hSim_Hemi_1Vtx_dist_ping_smu450_neu350->Fill( dist );
     if ( Msmu == 450 && Mneu == 400 ) hSim_Hemi_1Vtx_dist_ping_smu450_neu400->Fill( dist );
     if ( Msmu == 450 && Mneu == 430 ) hSim_Hemi_1Vtx_dist_ping_smu450_neu430->Fill( dist );
     if ( Msmu == 500 && Mneu == 180 ) hSim_Hemi_1Vtx_dist_ping_smu500_neu180->Fill( dist );
     if ( Msmu == 500 && Mneu == 200 ) hSim_Hemi_1Vtx_dist_ping_smu500_neu200->Fill( dist );
     if ( Msmu == 500 && Mneu == 250 ) hSim_Hemi_1Vtx_dist_ping_smu500_neu250->Fill( dist );
     if ( Msmu == 500 && Mneu == 300 ) hSim_Hemi_1Vtx_dist_ping_smu500_neu300->Fill( dist );
     if ( Msmu == 500 && Mneu == 350 ) hSim_Hemi_1Vtx_dist_ping_smu500_neu350->Fill( dist );
     if ( Msmu == 500 && Mneu == 400 ) hSim_Hemi_1Vtx_dist_ping_smu500_neu400->Fill( dist );
     if ( Msmu == 500 && Mneu == 450 ) hSim_Hemi_1Vtx_dist_ping_smu500_neu450->Fill( dist );
     if ( Msmu == 500 && Mneu == 480 ) hSim_Hemi_1Vtx_dist_ping_smu500_neu480->Fill( dist );
   }
   }
     
   if ( nVtx == 2 ) {
   if ( NChi2 > 0. && NChi2 < 10. && (step == 1 || step == 2) ) {
     if ( Msmu == 200 && Mneu == 180 ) hData_Hemi_2Vtx_dist_smu200_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 180 ) hData_Hemi_2Vtx_dist_smu250_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 200 ) hData_Hemi_2Vtx_dist_smu250_neu200->Fill( dist );
     if ( Msmu == 250 && Mneu == 230 ) hData_Hemi_2Vtx_dist_smu250_neu230->Fill( dist );
     if ( Msmu == 300 && Mneu == 180 ) hData_Hemi_2Vtx_dist_smu300_neu180->Fill( dist );
     if ( Msmu == 300 && Mneu == 200 ) hData_Hemi_2Vtx_dist_smu300_neu200->Fill( dist );
     if ( Msmu == 300 && Mneu == 250 ) hData_Hemi_2Vtx_dist_smu300_neu250->Fill( dist );
     if ( Msmu == 300 && Mneu == 280 ) hData_Hemi_2Vtx_dist_smu300_neu280->Fill( dist );
     if ( Msmu == 350 && Mneu == 180 ) hData_Hemi_2Vtx_dist_smu350_neu180->Fill( dist );
     if ( Msmu == 350 && Mneu == 200 ) hData_Hemi_2Vtx_dist_smu350_neu200->Fill( dist );
     if ( Msmu == 350 && Mneu == 250 ) hData_Hemi_2Vtx_dist_smu350_neu250->Fill( dist );
     if ( Msmu == 350 && Mneu == 300 ) hData_Hemi_2Vtx_dist_smu350_neu300->Fill( dist );
     if ( Msmu == 350 && Mneu == 330 ) hData_Hemi_2Vtx_dist_smu350_neu330->Fill( dist );
     if ( Msmu == 400 && Mneu == 180 ) hData_Hemi_2Vtx_dist_smu400_neu180->Fill( dist );
     if ( Msmu == 400 && Mneu == 200 ) hData_Hemi_2Vtx_dist_smu400_neu200->Fill( dist );
     if ( Msmu == 400 && Mneu == 250 ) hData_Hemi_2Vtx_dist_smu400_neu250->Fill( dist );
     if ( Msmu == 400 && Mneu == 300 ) hData_Hemi_2Vtx_dist_smu400_neu300->Fill( dist );
     if ( Msmu == 400 && Mneu == 350 ) hData_Hemi_2Vtx_dist_smu400_neu350->Fill( dist );
     if ( Msmu == 400 && Mneu == 380 ) hData_Hemi_2Vtx_dist_smu400_neu380->Fill( dist );
     if ( Msmu == 450 && Mneu == 180 ) hData_Hemi_2Vtx_dist_smu450_neu180->Fill( dist );
     if ( Msmu == 450 && Mneu == 200 ) hData_Hemi_2Vtx_dist_smu450_neu200->Fill( dist );
     if ( Msmu == 450 && Mneu == 250 ) hData_Hemi_2Vtx_dist_smu450_neu250->Fill( dist );
     if ( Msmu == 450 && Mneu == 300 ) hData_Hemi_2Vtx_dist_smu450_neu300->Fill( dist );
     if ( Msmu == 450 && Mneu == 350 ) hData_Hemi_2Vtx_dist_smu450_neu350->Fill( dist );
     if ( Msmu == 450 && Mneu == 400 ) hData_Hemi_2Vtx_dist_smu450_neu400->Fill( dist );
     if ( Msmu == 450 && Mneu == 430 ) hData_Hemi_2Vtx_dist_smu450_neu430->Fill( dist );
     if ( Msmu == 500 && Mneu == 180 ) hData_Hemi_2Vtx_dist_smu500_neu180->Fill( dist );
     if ( Msmu == 500 && Mneu == 200 ) hData_Hemi_2Vtx_dist_smu500_neu200->Fill( dist );
     if ( Msmu == 500 && Mneu == 250 ) hData_Hemi_2Vtx_dist_smu500_neu250->Fill( dist );
     if ( Msmu == 500 && Mneu == 300 ) hData_Hemi_2Vtx_dist_smu500_neu300->Fill( dist );
     if ( Msmu == 500 && Mneu == 350 ) hData_Hemi_2Vtx_dist_smu500_neu350->Fill( dist );
     if ( Msmu == 500 && Mneu == 400 ) hData_Hemi_2Vtx_dist_smu500_neu400->Fill( dist );
     if ( Msmu == 500 && Mneu == 450 ) hData_Hemi_2Vtx_dist_smu500_neu450->Fill( dist );
     if ( Msmu == 500 && Mneu == 480 ) hData_Hemi_2Vtx_dist_smu500_neu480->Fill( dist );
   }
    
   if ( NChi2 > 0. && NChi2 < 10. && (step == 1 || step == 2) && ping ) {
     if ( Msmu == 200 && Mneu == 180 ) hSim_Hemi_2Vtx_dist_ping_smu200_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 180 ) hSim_Hemi_2Vtx_dist_ping_smu250_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 200 ) hSim_Hemi_2Vtx_dist_ping_smu250_neu200->Fill( dist );
     if ( Msmu == 250 && Mneu == 230 ) hSim_Hemi_2Vtx_dist_ping_smu250_neu230->Fill( dist );
     if ( Msmu == 300 && Mneu == 180 ) hSim_Hemi_2Vtx_dist_ping_smu300_neu180->Fill( dist );
     if ( Msmu == 300 && Mneu == 200 ) hSim_Hemi_2Vtx_dist_ping_smu300_neu200->Fill( dist );
     if ( Msmu == 300 && Mneu == 250 ) hSim_Hemi_2Vtx_dist_ping_smu300_neu250->Fill( dist );
     if ( Msmu == 300 && Mneu == 280 ) hSim_Hemi_2Vtx_dist_ping_smu300_neu280->Fill( dist );
     if ( Msmu == 350 && Mneu == 180 ) hSim_Hemi_2Vtx_dist_ping_smu350_neu180->Fill( dist );
     if ( Msmu == 350 && Mneu == 200 ) hSim_Hemi_2Vtx_dist_ping_smu350_neu200->Fill( dist );
     if ( Msmu == 350 && Mneu == 250 ) hSim_Hemi_2Vtx_dist_ping_smu350_neu250->Fill( dist );
     if ( Msmu == 350 && Mneu == 300 ) hSim_Hemi_2Vtx_dist_ping_smu350_neu300->Fill( dist );
     if ( Msmu == 350 && Mneu == 330 ) hSim_Hemi_2Vtx_dist_ping_smu350_neu330->Fill( dist );
     if ( Msmu == 400 && Mneu == 180 ) hSim_Hemi_2Vtx_dist_ping_smu400_neu180->Fill( dist );
     if ( Msmu == 400 && Mneu == 200 ) hSim_Hemi_2Vtx_dist_ping_smu400_neu200->Fill( dist );
     if ( Msmu == 400 && Mneu == 250 ) hSim_Hemi_2Vtx_dist_ping_smu400_neu250->Fill( dist );
     if ( Msmu == 400 && Mneu == 300 ) hSim_Hemi_2Vtx_dist_ping_smu400_neu300->Fill( dist );
     if ( Msmu == 400 && Mneu == 350 ) hSim_Hemi_2Vtx_dist_ping_smu400_neu350->Fill( dist );
     if ( Msmu == 400 && Mneu == 380 ) hSim_Hemi_2Vtx_dist_ping_smu400_neu380->Fill( dist );
     if ( Msmu == 450 && Mneu == 180 ) hSim_Hemi_2Vtx_dist_ping_smu450_neu180->Fill( dist );
     if ( Msmu == 450 && Mneu == 200 ) hSim_Hemi_2Vtx_dist_ping_smu450_neu200->Fill( dist );
     if ( Msmu == 450 && Mneu == 250 ) hSim_Hemi_2Vtx_dist_ping_smu450_neu250->Fill( dist );
     if ( Msmu == 450 && Mneu == 300 ) hSim_Hemi_2Vtx_dist_ping_smu450_neu300->Fill( dist );
     if ( Msmu == 450 && Mneu == 350 ) hSim_Hemi_2Vtx_dist_ping_smu450_neu350->Fill( dist );
     if ( Msmu == 450 && Mneu == 400 ) hSim_Hemi_2Vtx_dist_ping_smu450_neu400->Fill( dist );
     if ( Msmu == 450 && Mneu == 430 ) hSim_Hemi_2Vtx_dist_ping_smu450_neu430->Fill( dist );
     if ( Msmu == 500 && Mneu == 180 ) hSim_Hemi_2Vtx_dist_ping_smu500_neu180->Fill( dist );
     if ( Msmu == 500 && Mneu == 200 ) hSim_Hemi_2Vtx_dist_ping_smu500_neu200->Fill( dist );
     if ( Msmu == 500 && Mneu == 250 ) hSim_Hemi_2Vtx_dist_ping_smu500_neu250->Fill( dist );
     if ( Msmu == 500 && Mneu == 300 ) hSim_Hemi_2Vtx_dist_ping_smu500_neu300->Fill( dist );
     if ( Msmu == 500 && Mneu == 350 ) hSim_Hemi_2Vtx_dist_ping_smu500_neu350->Fill( dist );
     if ( Msmu == 500 && Mneu == 400 ) hSim_Hemi_2Vtx_dist_ping_smu500_neu400->Fill( dist );
     if ( Msmu == 500 && Mneu == 450 ) hSim_Hemi_2Vtx_dist_ping_smu500_neu450->Fill( dist );
     if ( Msmu == 500 && Mneu == 480 ) hSim_Hemi_2Vtx_dist_ping_smu500_neu480->Fill( dist );
   }
   }
     
   if ( nVtxIni == 2 &&
        NChi2 > 0. && NChi2 < 10. && (step == 1 || step == 2) ) {
     if ( Msmu == 200 && Mneu == 180 ) hData_Hemi_2VtxIni_dist_smu200_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 180 ) hData_Hemi_2VtxIni_dist_smu250_neu180->Fill( dist );
     if ( Msmu == 250 && Mneu == 200 ) hData_Hemi_2VtxIni_dist_smu250_neu200->Fill( dist );
     if ( Msmu == 250 && Mneu == 230 ) hData_Hemi_2VtxIni_dist_smu250_neu230->Fill( dist );
     if ( Msmu == 300 && Mneu == 180 ) hData_Hemi_2VtxIni_dist_smu300_neu180->Fill( dist );
     if ( Msmu == 300 && Mneu == 200 ) hData_Hemi_2VtxIni_dist_smu300_neu200->Fill( dist );
     if ( Msmu == 300 && Mneu == 250 ) hData_Hemi_2VtxIni_dist_smu300_neu250->Fill( dist );
     if ( Msmu == 300 && Mneu == 280 ) hData_Hemi_2VtxIni_dist_smu300_neu280->Fill( dist );
     if ( Msmu == 350 && Mneu == 180 ) hData_Hemi_2VtxIni_dist_smu350_neu180->Fill( dist );
     if ( Msmu == 350 && Mneu == 200 ) hData_Hemi_2VtxIni_dist_smu350_neu200->Fill( dist );
     if ( Msmu == 350 && Mneu == 250 ) hData_Hemi_2VtxIni_dist_smu350_neu250->Fill( dist );
     if ( Msmu == 350 && Mneu == 300 ) hData_Hemi_2VtxIni_dist_smu350_neu300->Fill( dist );
     if ( Msmu == 350 && Mneu == 330 ) hData_Hemi_2VtxIni_dist_smu350_neu330->Fill( dist );
     if ( Msmu == 400 && Mneu == 180 ) hData_Hemi_2VtxIni_dist_smu400_neu180->Fill( dist );
     if ( Msmu == 400 && Mneu == 200 ) hData_Hemi_2VtxIni_dist_smu400_neu200->Fill( dist );
     if ( Msmu == 400 && Mneu == 250 ) hData_Hemi_2VtxIni_dist_smu400_neu250->Fill( dist );
     if ( Msmu == 400 && Mneu == 300 ) hData_Hemi_2VtxIni_dist_smu400_neu300->Fill( dist );
     if ( Msmu == 400 && Mneu == 350 ) hData_Hemi_2VtxIni_dist_smu400_neu350->Fill( dist );
     if ( Msmu == 400 && Mneu == 380 ) hData_Hemi_2VtxIni_dist_smu400_neu380->Fill( dist );
     if ( Msmu == 450 && Mneu == 180 ) hData_Hemi_2VtxIni_dist_smu450_neu180->Fill( dist );
     if ( Msmu == 450 && Mneu == 200 ) hData_Hemi_2VtxIni_dist_smu450_neu200->Fill( dist );
     if ( Msmu == 450 && Mneu == 250 ) hData_Hemi_2VtxIni_dist_smu450_neu250->Fill( dist );
     if ( Msmu == 450 && Mneu == 300 ) hData_Hemi_2VtxIni_dist_smu450_neu300->Fill( dist );
     if ( Msmu == 450 && Mneu == 350 ) hData_Hemi_2VtxIni_dist_smu450_neu350->Fill( dist );
     if ( Msmu == 450 && Mneu == 400 ) hData_Hemi_2VtxIni_dist_smu450_neu400->Fill( dist );
     if ( Msmu == 450 && Mneu == 430 ) hData_Hemi_2VtxIni_dist_smu450_neu430->Fill( dist );
     if ( Msmu == 500 && Mneu == 180 ) hData_Hemi_2VtxIni_dist_smu500_neu180->Fill( dist );
     if ( Msmu == 500 && Mneu == 200 ) hData_Hemi_2VtxIni_dist_smu500_neu200->Fill( dist );
     if ( Msmu == 500 && Mneu == 250 ) hData_Hemi_2VtxIni_dist_smu500_neu250->Fill( dist );
     if ( Msmu == 500 && Mneu == 300 ) hData_Hemi_2VtxIni_dist_smu500_neu300->Fill( dist );
     if ( Msmu == 500 && Mneu == 350 ) hData_Hemi_2VtxIni_dist_smu500_neu350->Fill( dist );
     if ( Msmu == 500 && Mneu == 400 ) hData_Hemi_2VtxIni_dist_smu500_neu400->Fill( dist );
     if ( Msmu == 500 && Mneu == 450 ) hData_Hemi_2VtxIni_dist_smu500_neu450->Fill( dist );
     if ( Msmu == 500 && Mneu == 480 ) hData_Hemi_2VtxIni_dist_smu500_neu480->Fill( dist );
     if ( minitree_Hemi_SecVtx->size() >= 1 ) {
       if ( Msmu == 200 && Mneu == 180 ) hData_Hemi_2VtxClose_dist_smu200_neu180->Fill( dist );
       if ( Msmu == 250 && Mneu == 180 ) hData_Hemi_2VtxClose_dist_smu250_neu180->Fill( dist );
       if ( Msmu == 250 && Mneu == 200 ) hData_Hemi_2VtxClose_dist_smu250_neu200->Fill( dist );
       if ( Msmu == 250 && Mneu == 230 ) hData_Hemi_2VtxClose_dist_smu250_neu230->Fill( dist );
       if ( Msmu == 300 && Mneu == 180 ) hData_Hemi_2VtxClose_dist_smu300_neu180->Fill( dist );
       if ( Msmu == 300 && Mneu == 200 ) hData_Hemi_2VtxClose_dist_smu300_neu200->Fill( dist );
       if ( Msmu == 300 && Mneu == 250 ) hData_Hemi_2VtxClose_dist_smu300_neu250->Fill( dist );
       if ( Msmu == 300 && Mneu == 280 ) hData_Hemi_2VtxClose_dist_smu300_neu280->Fill( dist );
       if ( Msmu == 350 && Mneu == 180 ) hData_Hemi_2VtxClose_dist_smu350_neu180->Fill( dist );
       if ( Msmu == 350 && Mneu == 200 ) hData_Hemi_2VtxClose_dist_smu350_neu200->Fill( dist );
       if ( Msmu == 350 && Mneu == 250 ) hData_Hemi_2VtxClose_dist_smu350_neu250->Fill( dist );
       if ( Msmu == 350 && Mneu == 300 ) hData_Hemi_2VtxClose_dist_smu350_neu300->Fill( dist );
       if ( Msmu == 350 && Mneu == 330 ) hData_Hemi_2VtxClose_dist_smu350_neu330->Fill( dist );
       if ( Msmu == 400 && Mneu == 180 ) hData_Hemi_2VtxClose_dist_smu400_neu180->Fill( dist );
       if ( Msmu == 400 && Mneu == 200 ) hData_Hemi_2VtxClose_dist_smu400_neu200->Fill( dist );
       if ( Msmu == 400 && Mneu == 250 ) hData_Hemi_2VtxClose_dist_smu400_neu250->Fill( dist );
       if ( Msmu == 400 && Mneu == 300 ) hData_Hemi_2VtxClose_dist_smu400_neu300->Fill( dist );
       if ( Msmu == 400 && Mneu == 350 ) hData_Hemi_2VtxClose_dist_smu400_neu350->Fill( dist );
       if ( Msmu == 400 && Mneu == 380 ) hData_Hemi_2VtxClose_dist_smu400_neu380->Fill( dist );
       if ( Msmu == 450 && Mneu == 180 ) hData_Hemi_2VtxClose_dist_smu450_neu180->Fill( dist );
       if ( Msmu == 450 && Mneu == 200 ) hData_Hemi_2VtxClose_dist_smu450_neu200->Fill( dist );
       if ( Msmu == 450 && Mneu == 250 ) hData_Hemi_2VtxClose_dist_smu450_neu250->Fill( dist );
       if ( Msmu == 450 && Mneu == 300 ) hData_Hemi_2VtxClose_dist_smu450_neu300->Fill( dist );
       if ( Msmu == 450 && Mneu == 350 ) hData_Hemi_2VtxClose_dist_smu450_neu350->Fill( dist );
       if ( Msmu == 450 && Mneu == 400 ) hData_Hemi_2VtxClose_dist_smu450_neu400->Fill( dist );
       if ( Msmu == 450 && Mneu == 430 ) hData_Hemi_2VtxClose_dist_smu450_neu430->Fill( dist );
       if ( Msmu == 500 && Mneu == 180 ) hData_Hemi_2VtxClose_dist_smu500_neu180->Fill( dist );
       if ( Msmu == 500 && Mneu == 200 ) hData_Hemi_2VtxClose_dist_smu500_neu200->Fill( dist );
       if ( Msmu == 500 && Mneu == 250 ) hData_Hemi_2VtxClose_dist_smu500_neu250->Fill( dist );
       if ( Msmu == 500 && Mneu == 300 ) hData_Hemi_2VtxClose_dist_smu500_neu300->Fill( dist );
       if ( Msmu == 500 && Mneu == 350 ) hData_Hemi_2VtxClose_dist_smu500_neu350->Fill( dist );
       if ( Msmu == 500 && Mneu == 400 ) hData_Hemi_2VtxClose_dist_smu500_neu400->Fill( dist );
       if ( Msmu == 500 && Mneu == 450 ) hData_Hemi_2VtxClose_dist_smu500_neu450->Fill( dist );
       if ( Msmu == 500 && Mneu == 480 ) hData_Hemi_2VtxClose_dist_smu500_neu480->Fill( dist );
     }
     if ( minitree_Hemi_SecVtx->size() == 2 && nVtx == 2 ) {
       if ( Msmu == 200 && Mneu == 180 ) hData_Hemi_2VtxMerge_dist_smu200_neu180->Fill( dist );
       if ( Msmu == 250 && Mneu == 180 ) hData_Hemi_2VtxMerge_dist_smu250_neu180->Fill( dist );
       if ( Msmu == 250 && Mneu == 200 ) hData_Hemi_2VtxMerge_dist_smu250_neu200->Fill( dist );
       if ( Msmu == 250 && Mneu == 230 ) hData_Hemi_2VtxMerge_dist_smu250_neu230->Fill( dist );
       if ( Msmu == 300 && Mneu == 180 ) hData_Hemi_2VtxMerge_dist_smu300_neu180->Fill( dist );
       if ( Msmu == 300 && Mneu == 200 ) hData_Hemi_2VtxMerge_dist_smu300_neu200->Fill( dist );
       if ( Msmu == 300 && Mneu == 250 ) hData_Hemi_2VtxMerge_dist_smu300_neu250->Fill( dist );
       if ( Msmu == 300 && Mneu == 280 ) hData_Hemi_2VtxMerge_dist_smu300_neu280->Fill( dist );
       if ( Msmu == 350 && Mneu == 180 ) hData_Hemi_2VtxMerge_dist_smu350_neu180->Fill( dist );
       if ( Msmu == 350 && Mneu == 200 ) hData_Hemi_2VtxMerge_dist_smu350_neu200->Fill( dist );
       if ( Msmu == 350 && Mneu == 250 ) hData_Hemi_2VtxMerge_dist_smu350_neu250->Fill( dist );
       if ( Msmu == 350 && Mneu == 300 ) hData_Hemi_2VtxMerge_dist_smu350_neu300->Fill( dist );
       if ( Msmu == 350 && Mneu == 330 ) hData_Hemi_2VtxMerge_dist_smu350_neu330->Fill( dist );
       if ( Msmu == 400 && Mneu == 180 ) hData_Hemi_2VtxMerge_dist_smu400_neu180->Fill( dist );
       if ( Msmu == 400 && Mneu == 200 ) hData_Hemi_2VtxMerge_dist_smu400_neu200->Fill( dist );
       if ( Msmu == 400 && Mneu == 250 ) hData_Hemi_2VtxMerge_dist_smu400_neu250->Fill( dist );
       if ( Msmu == 400 && Mneu == 300 ) hData_Hemi_2VtxMerge_dist_smu400_neu300->Fill( dist );
       if ( Msmu == 400 && Mneu == 350 ) hData_Hemi_2VtxMerge_dist_smu400_neu350->Fill( dist );
       if ( Msmu == 400 && Mneu == 380 ) hData_Hemi_2VtxMerge_dist_smu400_neu380->Fill( dist );
       if ( Msmu == 450 && Mneu == 180 ) hData_Hemi_2VtxMerge_dist_smu450_neu180->Fill( dist );
       if ( Msmu == 450 && Mneu == 200 ) hData_Hemi_2VtxMerge_dist_smu450_neu200->Fill( dist );
       if ( Msmu == 450 && Mneu == 250 ) hData_Hemi_2VtxMerge_dist_smu450_neu250->Fill( dist );
       if ( Msmu == 450 && Mneu == 300 ) hData_Hemi_2VtxMerge_dist_smu450_neu300->Fill( dist );
       if ( Msmu == 450 && Mneu == 350 ) hData_Hemi_2VtxMerge_dist_smu450_neu350->Fill( dist );
       if ( Msmu == 450 && Mneu == 400 ) hData_Hemi_2VtxMerge_dist_smu450_neu400->Fill( dist );
       if ( Msmu == 450 && Mneu == 430 ) hData_Hemi_2VtxMerge_dist_smu450_neu430->Fill( dist );
       if ( Msmu == 500 && Mneu == 180 ) hData_Hemi_2VtxMerge_dist_smu500_neu180->Fill( dist );
       if ( Msmu == 500 && Mneu == 200 ) hData_Hemi_2VtxMerge_dist_smu500_neu200->Fill( dist );
       if ( Msmu == 500 && Mneu == 250 ) hData_Hemi_2VtxMerge_dist_smu500_neu250->Fill( dist );
       if ( Msmu == 500 && Mneu == 300 ) hData_Hemi_2VtxMerge_dist_smu500_neu300->Fill( dist );
       if ( Msmu == 500 && Mneu == 350 ) hData_Hemi_2VtxMerge_dist_smu500_neu350->Fill( dist );
       if ( Msmu == 500 && Mneu == 400 ) hData_Hemi_2VtxMerge_dist_smu500_neu400->Fill( dist );
       if ( Msmu == 500 && Mneu == 450 ) hData_Hemi_2VtxMerge_dist_smu500_neu450->Fill( dist );
       if ( Msmu == 500 && Mneu == 480 ) hData_Hemi_2VtxMerge_dist_smu500_neu480->Fill( dist );
     }
   }

   if ( minitree_event_nVtx->at(0) >= 1 ) {
   if ( ( i == 0 && ((r0 < 60. && z0 < 108.) || (r0 < 73. && z0 >= 108. && z0 < 188.)) ) || 
        ( i == 1 && ((r1 < 60. && z1 < 108.) || (r1 < 73. && z0 >= 108. && z1 < 188.)) ) ) {
     float dSV0 = TMath::Sqrt( (minitree_Hemi_Vtx_x->at(i)-LLPx0)*(minitree_Hemi_Vtx_x->at(i)-LLPx0)
                             + (minitree_Hemi_Vtx_y->at(i)-LLPy0)*(minitree_Hemi_Vtx_y->at(i)-LLPy0)
                             + (minitree_Hemi_Vtx_z->at(i)-LLPz0)*(minitree_Hemi_Vtx_z->at(i)-LLPz0) );
     float dSV1 = TMath::Sqrt( (minitree_Hemi_Vtx_x->at(i)-LLPx1)*(minitree_Hemi_Vtx_x->at(i)-LLPx1)
                             + (minitree_Hemi_Vtx_y->at(i)-LLPy1)*(minitree_Hemi_Vtx_y->at(i)-LLPy1)
                             + (minitree_Hemi_Vtx_z->at(i)-LLPz1)*(minitree_Hemi_Vtx_z->at(i)-LLPz1) );
     float dSV = dSV0;
     if ( dSV1 < dSV0 ) dSV = dSV1;
     int Vtx_step = minitree_Hemi_Vtx_step->at(i);
     float SVdist = minitree_Hemi_Vtx_dist->at(i);
     if ( Vtx_step >= 1 && Vtx_step <= 2 ) {
       if ( abs(minitree_Hemi_eta->at(i)) < 1.5 ) {                                      // 68%    / 90%
         if (     SVdist < 0.15 ) hSim_Hemi_dSV_step12_etaLT15_dist001->Fill( dSV ); // 0.0045 / 0.04
	 else if ( SVdist < 0.5 ) hSim_Hemi_dSV_step12_etaLT15_dist003->Fill( dSV ); // 0.0045 / 0.06
	 else if ( SVdist < 1.5 ) hSim_Hemi_dSV_step12_etaLT15_dist010->Fill( dSV ); // 0.005  / 0.06
	 else if ( SVdist < 5.0 ) hSim_Hemi_dSV_step12_etaLT15_dist030->Fill( dSV ); // 0.01   / 0.2
	 else if ( SVdist < 15. ) hSim_Hemi_dSV_step12_etaLT15_dist100->Fill( dSV ); // 0.015  / 0.25
	 else if ( SVdist < 25. ) hSim_Hemi_dSV_step12_etaLT15_dist200->Fill( dSV ); // 0.04   / 0.5
	 else if ( SVdist < 35. ) hSim_Hemi_dSV_step12_etaLT15_dist300->Fill( dSV ); // 0.12   / 1.8
	 else if ( SVdist < 45. ) hSim_Hemi_dSV_step12_etaLT15_dist400->Fill( dSV ); // 0.6    / 4.5
	 else if ( SVdist < 55. ) hSim_Hemi_dSV_step12_etaLT15_dist500->Fill( dSV ); // 1.8    / 7.0 
	 else                     hSim_Hemi_dSV_step12_etaLT15_dist999->Fill( dSV ); // 3.5    / 12.0
       }
       else {
         if (     SVdist < 0.15 ) hSim_Hemi_dSV_step12_etaGT15_dist001->Fill( dSV ); // 0.018  / 0.10
	 else if ( SVdist < 0.5 ) hSim_Hemi_dSV_step12_etaGT15_dist003->Fill( dSV ); // 0.011  / 0.10
	 else if ( SVdist < 1.5 ) hSim_Hemi_dSV_step12_etaGT15_dist010->Fill( dSV ); // 0.010  / 0.15 
	 else if ( SVdist < 5.0 ) hSim_Hemi_dSV_step12_etaGT15_dist030->Fill( dSV ); // 0.015  / 0.2
	 else if ( SVdist < 15. ) hSim_Hemi_dSV_step12_etaGT15_dist100->Fill( dSV ); // 0.025  / 0.5
	 else if ( SVdist < 25. ) hSim_Hemi_dSV_step12_etaGT15_dist200->Fill( dSV ); // 0.05   / 1.0
	 else if ( SVdist < 35. ) hSim_Hemi_dSV_step12_etaGT15_dist300->Fill( dSV ); // 0.08   / 0.8
	 else if ( SVdist < 45. ) hSim_Hemi_dSV_step12_etaGT15_dist400->Fill( dSV ); // 0.16   / 1.5
	 else if ( SVdist < 55. ) hSim_Hemi_dSV_step12_etaGT15_dist500->Fill( dSV ); // 0.3    / 4.0
	 else                     hSim_Hemi_dSV_step12_etaGT15_dist999->Fill( dSV ); // 0.8    / 6.0
       }
     }
     else if ( Vtx_step >= 3 && Vtx_step <= 4 ) {
       if ( abs(minitree_Hemi_eta->at(i)) < 1.5 ) {
         if (     SVdist < 0.15 ) hSim_Hemi_dSV_step34_etaLT15_dist001->Fill( dSV ); // low stat
	 else if ( SVdist < 0.5 ) hSim_Hemi_dSV_step34_etaLT15_dist003->Fill( dSV ); // low stat
	 else if ( SVdist < 1.5 ) hSim_Hemi_dSV_step34_etaLT15_dist010->Fill( dSV ); // low stat
	 else if ( SVdist < 5.0 ) hSim_Hemi_dSV_step34_etaLT15_dist030->Fill( dSV ); // 0.8    / 4.5 
	 else if ( SVdist < 15. ) hSim_Hemi_dSV_step34_etaLT15_dist100->Fill( dSV ); // 2.5    / > 10.
	 else if ( SVdist < 25. ) hSim_Hemi_dSV_step34_etaLT15_dist200->Fill( dSV ); // 2.0    / 18.0
	 else if ( SVdist < 35. ) hSim_Hemi_dSV_step34_etaLT15_dist300->Fill( dSV ); // 0.060 
	 else if ( SVdist < 45. ) hSim_Hemi_dSV_step34_etaLT15_dist400->Fill( dSV ); // 0.060 
	 else if ( SVdist < 55. ) hSim_Hemi_dSV_step34_etaLT15_dist500->Fill( dSV ); // 0.060 
	 else                     hSim_Hemi_dSV_step34_etaLT15_dist999->Fill( dSV ); // > 2 
       }
       else {
         if (     SVdist < 0.15 ) hSim_Hemi_dSV_step34_etaGT15_dist001->Fill( dSV ); // low stat
	 else if ( SVdist < 0.5 ) hSim_Hemi_dSV_step34_etaGT15_dist003->Fill( dSV ); // 0.1 low stat
	 else if ( SVdist < 1.5 ) hSim_Hemi_dSV_step34_etaGT15_dist010->Fill( dSV ); // 0.1 low stat 
	 else if ( SVdist < 5.0 ) hSim_Hemi_dSV_step34_etaGT15_dist030->Fill( dSV ); // 0.2
	 else if ( SVdist < 15. ) hSim_Hemi_dSV_step34_etaGT15_dist100->Fill( dSV ); // 0.4
	 else if ( SVdist < 25. ) hSim_Hemi_dSV_step34_etaGT15_dist200->Fill( dSV ); // 0.060 
	 else if ( SVdist < 35. ) hSim_Hemi_dSV_step34_etaGT15_dist300->Fill( dSV ); // 0.060 
	 else if ( SVdist < 45. ) hSim_Hemi_dSV_step34_etaGT15_dist400->Fill( dSV ); // 0.060 
	 else if ( SVdist < 55. ) hSim_Hemi_dSV_step34_etaGT15_dist500->Fill( dSV ); // 0.060 
	 else                     hSim_Hemi_dSV_step34_etaGT15_dist999->Fill( dSV ); // > 2 
       }
     }
   }
   }

   } // end loop on hemispheres

 } // end loop on events 

 std::cout << "number of events  "<< allevents << std::endl; 

//################################################

// Output Postscript

//   TCanvas* c = new TCanvas("c");
  hGen_Msmu -> Draw(); 
//   c->Print("output.ps(");
  
//################################################
  HistogramManager h ;
  
  h.WriteAllHistogramsInFile(filename.Data(),"recreate");
//################################################
}   