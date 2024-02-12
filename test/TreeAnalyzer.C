#define TreeAnalyzer_cxx
#include "TreeAnalyzer.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "./HistogramManager.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime> 
#include "../interface/DeltaFunc.h"
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


void TreeAnalyzer::Loop(TString sample)
{
//   In a ROOT session, you can do:
//      root> .L TreeAnalyzer.C
//      root> TreeAnalyzer t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  TString thesample  = sample;
//   TString filename  =  "BKGEstimation_"+thesample+".root";
  TFile * theoutputfile = new TFile( ("outputroot/BKGEstimation_"+thesample+".root").Data() , "recreate");
  std::ofstream ofs ("outputroot/BKGEstimation_"+thesample+".txt", std::ofstream::out);

//**********************************
// Histograms
//**********************************

 TH1F* hSim_weight = new TH1F("hSim_weight","",100,-10.,10.);

 TH1F* hData_Filter     = new TH1F("hData_Filter","",2,-0.5,1.5);
 TH1F* hData_njetNOmuAll  = new TH1F("hData_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_njetAll      = new TH1F("hData_njetAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_eta     = new TH1F("hData_Hemi_eta","",25,0.,7.5);

 TH1F*     hData_MmumuABCD     = new TH1F("hData_MmumuABCD","",25,0.,500.);
 TH1F*     hData_BDTevtABCD    = new TH1F("hData_BDTevtABCD","",40,-1.,1.);
 TH1F*     hData_BDTevtDY  = new TH1F("hData_BDTevtDY","",40,-1.,1.);
 TH1F*     hData_BDTevtTT  = new TH1F("hData_BDTevtTT","",40,-1.,1.);

 TH1F*     hData_BDTHemi1     = new TH1F("hData_BDTHemi1","",40,-1.,1.);
 TH1F*     hData_BDTHemi1DY   = new TH1F("hData_BDTHemi1DY","",40,-1.,1.);
 TH1F*     hData_BDTHemi1TT   = new TH1F("hData_BDTHemi1TT","",40,-1.,1.);

 TH1F*     hData_BDTHemi2     = new TH1F("hData_BDTHemi2","",40,-1.,1.);
 TH1F*     hData_BDTHemi2DY   = new TH1F("hData_BDTHemi2DY","",40,-1.,1.);
 TH1F*     hData_BDTHemi2TT   = new TH1F("hData_BDTHemi2TT","",40,-1.,1.);

// ABCD
   TH1F* hData_SR_1Vtx_MVtx               = new TH1F("hData_SR_1Vtx_MVtx","",25,0.,100.);
   TH1F* hData_SR_1Vtx_MHemi              = new TH1F("hData_SR_1Vtx_MHemi","",25,0.,500.);
   TH1F* hData_SR_1Vtx_Mmumu              = new TH1F("hData_SR_1Vtx_Mmumu","",50,0.,200.);

   TH1F* hData_CRNoEvt_1Vtx_MVtx          = new TH1F("hData_CRNoEvt_1Vtx_MVtx","",25,0.,100.);
   TH1F* hData_CRNoEvt_1Vtx_MHemi         = new TH1F("hData_CRNoEvt_1Vtx_MHemi","",25,0.,500.);
   TH1F* hData_CRNoEvt_1Vtx_Mmumu         = new TH1F("hData_CRNoEvt_1Vtx_Mmumu","",50,0.,200.); 

   TH1F* hData_CRNoEvtNoVtx_1Vtx_MVtx     = new TH1F("hData_CRNoEvtNoVtx_1Vtx_MVtx","",25,0.,100.);
   TH1F* hData_CRNoEvtNoVtx_1Vtx_MHemi    = new TH1F("hData_CRNoEvtNoVtx_1Vtx_MHemi","",25,0.,500.);
   TH1F* hData_CRNoEvtNoVtx_1Vtx_Mmumu    = new TH1F("hData_CRNoEvtNoVtx_1Vtx_Mmumu","",50,0.,200.); 

   TH1F* hData_CRNoVtx_1Vtx_MVtx          = new TH1F("hData_CRNoVtx_1Vtx_MVtx","",25,0.,100.);
   TH1F* hData_CRNoVtx_1Vtx_MHemi         = new TH1F("hData_CRNoVtx_1Vtx_MHemi","",25,0.,500.);
   TH1F* hData_CRNoVtx_1Vtx_Mmumu         = new TH1F("hData_CRNoVtx_1Vtx_Mmumu","",50,0.,200.); 

   TH1F* hData_SR_2Vtx_MVtx               = new TH1F("hData_SR_2Vtx_MVtx","",25,0.,100.);
   TH1F* hData_SR_2Vtx_MHemi              = new TH1F("hData_SR_2Vtx_MHemi","",25,0.,500.);
   TH1F* hData_SR_2Vtx_Mmumu              = new TH1F("hData_SR_2Vtx_Mmumu","",50,0.,200.);

   TH1F* hData_CRNoEvt_2Vtx_MVtx          = new TH1F("hData_CRNoEvt_2Vtx_MVtx","",25,0.,100.);
   TH1F* hData_CRNoEvt_2Vtx_MHemi         = new TH1F("hData_CRNoEvt_2Vtx_MHemi","",25,0.,500.);
   TH1F* hData_CRNoEvt_2Vtx_Mmumu         = new TH1F("hData_CRNoEvt_2Vtx_Mmumu","",50,0.,200.);

   TH1F* hData_CRNoEvtNoVtx_2Vtx_MVtx     = new TH1F("hData_CRNoEvtNoVtx_2Vtx_MVtx","",25,0.,100.);
   TH1F* hData_CRNoEvtNoVtx_2Vtx_MHemi    = new TH1F("hData_CRNoEvtNoVtx_2Vtx_MHemi","",25,0.,500.);
   TH1F* hData_CRNoEvtNoVtx_2Vtx_Mmumu    = new TH1F("hData_CRNoEvtNoVtx_2Vtx_Mmumu","",50,0.,200.); 

   TH1F* hData_CRNoVtx_2Vtx_MVtx          = new TH1F("hData_CRNoVtx_2Vtx_MVtx","",25,0.,100.);
   TH1F* hData_CRNoVtx_2Vtx_MHemi         = new TH1F("hData_CRNoVtx_2Vtx_MHemi","",25,0.,500.);
   TH1F* hData_CRNoVtx_2Vtx_Mmumu         = new TH1F("hData_CRNoVtx_2Vtx_Mmumu","",50,0.,200.);


   TH1F* hData2_SR_1Vtx_MVtx               = new TH1F("hData2_SR_1Vtx_MVtx","",25,0.,100.);
   TH1F* hData2_SR_1Vtx_MHemi              = new TH1F("hData2_SR_1Vtx_MHemi","",25,0.,500.);
   TH1F* hData2_SR_1Vtx_Mmumu              = new TH1F("hData2_SR_1Vtx_Mmumu","",50,0.,200.);

   TH1F* hData2_CRNoEvt_1Vtx_MVtx          = new TH1F("hData2_CRNoEvt_1Vtx_MVtx","",25,0.,100.);
   TH1F* hData2_CRNoEvt_1Vtx_MHemi         = new TH1F("hData2_CRNoEvt_1Vtx_MHemi","",25,0.,500.);
   TH1F* hData2_CRNoEvt_1Vtx_Mmumu         = new TH1F("hData2_CRNoEvt_1Vtx_Mmumu","",50,0.,200.); 

   TH1F* hData2_CRNoEvtNoVtx_1Vtx_MVtx     = new TH1F("hData2_CRNoEvtNoVtx_1Vtx_MVtx","",25,0.,100.);
   TH1F* hData2_CRNoEvtNoVtx_1Vtx_MHemi    = new TH1F("hData2_CRNoEvtNoVtx_1Vtx_MHemi","",25,0.,500.);
   TH1F* hData2_CRNoEvtNoVtx_1Vtx_Mmumu    = new TH1F("hData2_CRNoEvtNoVtx_1Vtx_Mmumu","",50,0.,200.); 

   TH1F* hData2_CRNoVtx_1Vtx_MVtx          = new TH1F("hData2_CRNoVtx_1Vtx_MVtx","",25,0.,100.);
   TH1F* hData2_CRNoVtx_1Vtx_MHemi         = new TH1F("hData2_CRNoVtx_1Vtx_MHemi","",25,0.,500.);
   TH1F* hData2_CRNoVtx_1Vtx_Mmumu         = new TH1F("hData2_CRNoVtx_1Vtx_Mmumu","",50,0.,200.); 

   TH1F* hData2_SR_2Vtx_MVtx               = new TH1F("hData2_SR_2Vtx_MVtx","",25,0.,100.);
   TH1F* hData2_SR_2Vtx_MHemi              = new TH1F("hData2_SR_2Vtx_MHemi","",25,0.,500.);
   TH1F* hData2_SR_2Vtx_Mmumu              = new TH1F("hData2_SR_2Vtx_Mmumu","",50,0.,200.);

   TH1F* hData2_CRNoEvt_2Vtx_MVtx          = new TH1F("hData2_CRNoEvt_2Vtx_MVtx","",25,0.,100.);
   TH1F* hData2_CRNoEvt_2Vtx_MHemi         = new TH1F("hData2_CRNoEvt_2Vtx_MHemi","",25,0.,500.);
   TH1F* hData2_CRNoEvt_2Vtx_Mmumu         = new TH1F("hData2_CRNoEvt_2Vtx_Mmumu","",50,0.,200.);

   TH1F* hData2_CRNoEvtNoVtx_2Vtx_MVtx     = new TH1F("hData2_CRNoEvtNoVtx_2Vtx_MVtx","",25,0.,100.);
   TH1F* hData2_CRNoEvtNoVtx_2Vtx_MHemi    = new TH1F("hData2_CRNoEvtNoVtx_2Vtx_MHemi","",25,0.,500.);
   TH1F* hData2_CRNoEvtNoVtx_2Vtx_Mmumu    = new TH1F("hData2_CRNoEvtNoVtx_2Vtx_Mmumu","",50,0.,200.); 

   TH1F* hData2_CRNoVtx_2Vtx_MVtx          = new TH1F("hData2_CRNoVtx_2Vtx_MVtx","",25,0.,100.);
   TH1F* hData2_CRNoVtx_2Vtx_MHemi         = new TH1F("hData2_CRNoVtx_2Vtx_MHemi","",25,0.,500.);
   TH1F* hData2_CRNoVtx_2Vtx_Mmumu         = new TH1F("hData2_CRNoVtx_2Vtx_Mmumu","",50,0.,200.);
// END OF ABCD


 TH1F* hData_Mmumu        = new TH1F("hData_Mmumu","",25,0.,500.);
 TH1F* hData_BDTevt       = new TH1F("hData_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_lead_ptmin = new TH1F("hData_Hemi_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_ptmin      = new TH1F("hData_Hemi_ptmin","",25,0.,1000.);
 TH1F* hData_Hemi_nTrks   = new TH1F("hData_Hemi_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_njetNOmu= new TH1F("hData_Hemi_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_mass    = new TH1F("hData_Hemi_mass","",25,0.,500.);
 TH1F* hData_Hemi_massMu  = new TH1F("hData_Hemi_massMu","",25,0.,750.);
 TH1F* hData_Hemi_ptMu    = new TH1F("hData_Hemi_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_dRMu    = new TH1F("hData_Hemi_dRMu","",25,0.,7.5);
 TH1F* hData_CutEvt_Mmumu        = new TH1F("hData_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_CutEvt_lead_ptmin = new TH1F("hData_Hemi_CutEvt_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_CutEvt_ptmin   = new TH1F("hData_Hemi_CutEvt_ptmin","",25,0.,1000.);
 TH1F* hData_Hemi_CutEvt_nTrks   = new TH1F("hData_Hemi_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_CutEvt_njetNOmu= new TH1F("hData_Hemi_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_CutEvt_mass    = new TH1F("hData_Hemi_CutEvt_mass","",25,0.,500.);
 TH1F* hData_Hemi_CutEvt_massMu  = new TH1F("hData_Hemi_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_Hemi_CutEvt_ptMu    = new TH1F("hData_Hemi_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_CutEvt_dRMu    = new TH1F("hData_Hemi_CutEvt_dRMu","",25,0.,7.5);

 TH1F* hData_Hemi_0Vtx_BDTevt  = new TH1F("hData_Hemi_0Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_0Vtx_nTrks   = new TH1F("hData_Hemi_0Vtx_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_0Vtx_mass    = new TH1F("hData_Hemi_0Vtx_mass","",25,0.,500.);
 TH1F* hData_Hemi_0Vtx_Mmumu   = new TH1F("hData_Hemi_0Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_0Vtx_PFMet   = new TH1F("hData_Hemi_0Vtx_PFMet","",25,0.,250.);
 TH1F* hData_Hemi_0Vtx_njetNOmuAll= new TH1F("hData_Hemi_0Vtx_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_0Vtx_lead_ptmin = new TH1F("hData_Hemi_0Vtx_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_0Vtx_ptmin      = new TH1F("hData_Hemi_0Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_Hemi_0Vtx_CutEvt_nTrks       = new TH1F("hData_Hemi_0Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_0Vtx_CutEvt_mass        = new TH1F("hData_Hemi_0Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_Hemi_0Vtx_CutEvt_Mmumu       = new TH1F("hData_Hemi_0Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_0Vtx_CutEvt_njetNOmuAll = new TH1F("hData_Hemi_0Vtx_CutEvt_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_0Vtx_CutEvt_lead_ptmin  = new TH1F("hData_Hemi_0Vtx_CutEvt_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_0Vtx_CutEvt_ptmin       = new TH1F("hData_Hemi_0Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_Hemi_1Vtx_BDTevt  = new TH1F("hData_Hemi_1Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_1Vtx_BDTvtx  = new TH1F("hData_Hemi_1Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_1Vtx_njetNOmu= new TH1F("hData_Hemi_1Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_1Vtx_mass    = new TH1F("hData_Hemi_1Vtx_mass","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_nTrks   = new TH1F("hData_Hemi_1Vtx_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_1Vtx_Mass    = new TH1F("hData_Hemi_1Vtx_Mass","",25,0.,100.);
 TH1F* hData_Hemi_1Vtx_massMu  = new TH1F("hData_Hemi_1Vtx_massMu","",25,0.,750.);
 TH1F* hData_Hemi_1Vtx_ptMu    = new TH1F("hData_Hemi_1Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_1Vtx_dRMu    = new TH1F("hData_Hemi_1Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_Hemi_1Vtx_Mmumu   = new TH1F("hData_Hemi_1Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_PFMet   = new TH1F("hData_Hemi_1Vtx_PFMet","",25,0.,250.);
 TH1F* hData_Hemi_1Vtx_njetNOmuAll= new TH1F("hData_Hemi_1Vtx_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_1Vtx_lead_ptmin = new TH1F("hData_Hemi_1Vtx_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_ptmin      = new TH1F("hData_Hemi_1Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_Hemi_2Vtx_BDTevt  = new TH1F("hData_Hemi_2Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_2Vtx_BDTvtx  = new TH1F("hData_Hemi_2Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_2Vtx_njetNOmu= new TH1F("hData_Hemi_2Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_2Vtx_mass    = new TH1F("hData_Hemi_2Vtx_mass","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_nTrks   = new TH1F("hData_Hemi_2Vtx_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_2Vtx_Mass    = new TH1F("hData_Hemi_2Vtx_Mass","",25,0.,100.);
 TH1F* hData_Hemi_2Vtx_massMu  = new TH1F("hData_Hemi_2Vtx_massMu","",25,0.,750.);
 TH1F* hData_Hemi_2Vtx_ptMu    = new TH1F("hData_Hemi_2Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_2Vtx_dRMu    = new TH1F("hData_Hemi_2Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_Hemi_2Vtx_Mmumu   = new TH1F("hData_Hemi_2Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_PFMet   = new TH1F("hData_Hemi_2Vtx_PFMet","",25,0.,250.);
 TH1F* hData_Hemi_2Vtx_njetNOmuAll= new TH1F("hData_Hemi_2Vtx_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_2Vtx_lead_ptmin = new TH1F("hData_Hemi_2Vtx_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_ptmin	  = new TH1F("hData_Hemi_2Vtx_ptmin","",25,0.,1000.);
 TH1F* hData_Hemi_2VtxAll_BDTvtx  = new TH1F("hData_Hemi_2VtxAll_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_2VtxAll_njetNOmu= new TH1F("hData_Hemi_2VtxAll_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_2VtxAll_mass    = new TH1F("hData_Hemi_2VtxAll_mass","",25,0.,500.);
 TH1F* hData_Hemi_2VtxAll_nTrks   = new TH1F("hData_Hemi_2VtxAll_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_2VtxAll_Mass    = new TH1F("hData_Hemi_2VtxAll_Mass","",25,0.,100.);
 TH1F* hData_Hemi_2VtxAll_massMu  = new TH1F("hData_Hemi_2VtxAll_massMu","",25,0.,750.);
 TH1F* hData_Hemi_2VtxAll_ptMu    = new TH1F("hData_Hemi_2VtxAll_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_2VtxAll_dRMu    = new TH1F("hData_Hemi_2VtxAll_dRMu","",25,0.,7.5);

 TH1F* hData_Hemi_1Vtx_CutEvt_BDTevt  = new TH1F("hData_Hemi_1Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_1Vtx_CutEvt_BDTvtx  = new TH1F("hData_Hemi_1Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_1Vtx_CutEvt_njetNOmu= new TH1F("hData_Hemi_1Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_1Vtx_CutEvt_mass    = new TH1F("hData_Hemi_1Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_CutEvt_nTrks   = new TH1F("hData_Hemi_1Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_1Vtx_CutEvt_Mass    = new TH1F("hData_Hemi_1Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_1Vtx_CutEvt_massMu  = new TH1F("hData_Hemi_1Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_Hemi_1Vtx_CutEvt_ptMu    = new TH1F("hData_Hemi_1Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_1Vtx_CutEvt_dRMu    = new TH1F("hData_Hemi_1Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_Hemi_1Vtx_CutEvt_Mmumu   = new TH1F("hData_Hemi_1Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_CutEvt_njetNOmuAll= new TH1F("hData_Hemi_1Vtx_CutEvt_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_1Vtx_CutEvt_lead_ptmin = new TH1F("hData_Hemi_1Vtx_CutEvt_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_CutEvt_ptmin      = new TH1F("hData_Hemi_1Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_Hemi_2Vtx_CutEvt_BDTevt  = new TH1F("hData_Hemi_2Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_2Vtx_CutEvt_BDTvtx  = new TH1F("hData_Hemi_2Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_2Vtx_CutEvt_njetNOmu= new TH1F("hData_Hemi_2Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_2Vtx_CutEvt_mass    = new TH1F("hData_Hemi_2Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_CutEvt_nTrks   = new TH1F("hData_Hemi_2Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_2Vtx_CutEvt_Mass    = new TH1F("hData_Hemi_2Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_2Vtx_CutEvt_massMu  = new TH1F("hData_Hemi_2Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_Hemi_2Vtx_CutEvt_ptMu    = new TH1F("hData_Hemi_2Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_2Vtx_CutEvt_dRMu    = new TH1F("hData_Hemi_2Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_Hemi_2Vtx_CutEvt_Mmumu   = new TH1F("hData_Hemi_2Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_CutEvt_njetNOmuAll= new TH1F("hData_Hemi_2Vtx_CutEvt_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_2Vtx_CutEvt_lead_ptmin = new TH1F("hData_Hemi_2Vtx_CutEvt_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_CutEvt_ptmin	 = new TH1F("hData_Hemi_2Vtx_CutEvt_ptmin","",25,0.,1000.);
 TH1F* hData_Hemi_2VtxAll_CutEvt_BDTvtx  = new TH1F("hData_Hemi_2VtxAll_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_2VtxAll_CutEvt_njetNOmu= new TH1F("hData_Hemi_2VtxAll_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_2VtxAll_CutEvt_mass    = new TH1F("hData_Hemi_2VtxAll_CutEvt_mass","",25,0.,500.);
 TH1F* hData_Hemi_2VtxAll_CutEvt_nTrks   = new TH1F("hData_Hemi_2VtxAll_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_2VtxAll_CutEvt_Mass    = new TH1F("hData_Hemi_2VtxAll_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_2VtxAll_CutEvt_massMu  = new TH1F("hData_Hemi_2VtxAll_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_Hemi_2VtxAll_CutEvt_ptMu    = new TH1F("hData_Hemi_2VtxAll_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_2VtxAll_CutEvt_dRMu    = new TH1F("hData_Hemi_2VtxAll_CutEvt_dRMu","",25,0.,7.5);

 TH1F* hData_Hemi_1Vtx_CutVtx_BDTevt  = new TH1F("hData_Hemi_1Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_1Vtx_CutVtx_BDTvtx  = new TH1F("hData_Hemi_1Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_1Vtx_CutVtx_njetNOmu= new TH1F("hData_Hemi_1Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_1Vtx_CutVtx_mass    = new TH1F("hData_Hemi_1Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_CutVtx_nTrks   = new TH1F("hData_Hemi_1Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_1Vtx_CutVtx_Mass    = new TH1F("hData_Hemi_1Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_Hemi_1Vtx_CutVtx_massMu  = new TH1F("hData_Hemi_1Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_Hemi_1Vtx_CutVtx_ptMu    = new TH1F("hData_Hemi_1Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_1Vtx_CutVtx_dRMu    = new TH1F("hData_Hemi_1Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_Hemi_1Vtx_CutVtx_Mmumu   = new TH1F("hData_Hemi_1Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_CutVtx_PFMet   = new TH1F("hData_Hemi_1Vtx_CutVtx_PFMet","",25,0.,250.);
 TH1F* hData_Hemi_1Vtx_CutVtx_njetNOmuAll= new TH1F("hData_Hemi_1Vtx_CutVtx_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_1Vtx_CutVtx_lead_ptmin = new TH1F("hData_Hemi_1Vtx_CutVtx_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_CutVtx_ptmin      = new TH1F("hData_Hemi_1Vtx_CutVtx_ptmin","",25,0.,1000.);

 TH1F* hData_Hemi_2Vtx_CutVtx_BDTevt  = new TH1F("hData_Hemi_2Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_2Vtx_CutVtx_BDTvtx  = new TH1F("hData_Hemi_2Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_2Vtx_CutVtx_njetNOmu= new TH1F("hData_Hemi_2Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_2Vtx_CutVtx_mass    = new TH1F("hData_Hemi_2Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_CutVtx_nTrks   = new TH1F("hData_Hemi_2Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_2Vtx_CutVtx_Mass    = new TH1F("hData_Hemi_2Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_Hemi_2Vtx_CutVtx_massMu  = new TH1F("hData_Hemi_2Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_Hemi_2Vtx_CutVtx_ptMu    = new TH1F("hData_Hemi_2Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_2Vtx_CutVtx_dRMu    = new TH1F("hData_Hemi_2Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_Hemi_2Vtx_CutVtx_Mmumu   = new TH1F("hData_Hemi_2Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_CutVtx_PFMet   = new TH1F("hData_Hemi_2Vtx_CutVtx_PFMet","",25,0.,250.);
 TH1F* hData_Hemi_2Vtx_CutVtx_njetNOmuAll= new TH1F("hData_Hemi_2Vtx_CutVtx_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_2Vtx_CutVtx_lead_ptmin = new TH1F("hData_Hemi_2Vtx_CutVtx_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_CutVtx_ptmin	 = new TH1F("hData_Hemi_2Vtx_CutVtx_ptmin","",25,0.,1000.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_BDTvtx  = new TH1F("hData_Hemi_2VtxAll_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_njetNOmu= new TH1F("hData_Hemi_2VtxAll_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_2VtxAll_CutVtx_mass    = new TH1F("hData_Hemi_2VtxAll_CutVtx_mass","",25,0.,500.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_nTrks   = new TH1F("hData_Hemi_2VtxAll_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_Mass    = new TH1F("hData_Hemi_2VtxAll_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_massMu  = new TH1F("hData_Hemi_2VtxAll_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_ptMu    = new TH1F("hData_Hemi_2VtxAll_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_dRMu    = new TH1F("hData_Hemi_2VtxAll_CutVtx_dRMu","",25,0.,7.5);

 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_njetNOmuAll= new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_lead_ptmin = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_1Vtx_CutVtx_CutEvt_ptmin      = new TH1F("hData_Hemi_1Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_njetNOmuAll= new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_lead_ptmin = new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_lead_ptmin","",25,0.,500.);
 TH1F* hData_Hemi_2Vtx_CutVtx_CutEvt_ptmin	= new TH1F("hData_Hemi_2Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_Hemi_2VtxAll_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_CutEvt_njetNOmu= new TH1F("hData_Hemi_2VtxAll_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_Hemi_2VtxAll_CutVtx_CutEvt_mass    = new TH1F("hData_Hemi_2VtxAll_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_CutEvt_nTrks   = new TH1F("hData_Hemi_2VtxAll_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_CutEvt_Mass    = new TH1F("hData_Hemi_2VtxAll_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_CutEvt_massMu  = new TH1F("hData_Hemi_2VtxAll_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_CutEvt_ptMu    = new TH1F("hData_Hemi_2VtxAll_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_Hemi_2VtxAll_CutVtx_CutEvt_dRMu    = new TH1F("hData_Hemi_2VtxAll_CutVtx_CutEvt_dRMu","",25,0.,7.5);

 TH1F* hData_Hemi_2SecVtx_CutEvt_Mass  = new TH1F("hData_Hemi_2SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_2SecVtx_CutVtx_CutEvt_Mass = new TH1F("hData_Hemi_2SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_1SecVtx_CutEvt_Mass = new TH1F("hData_Hemi_1SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_1SecVtx_CutEvt_Merge = new TH1F("hData_Hemi_1SecVtx_CutEvt_Merge","",3,-1,1);
 TH1F* hData_Hemi_1SecVtx_CutVtx_CutEvt_Mass = new TH1F("hData_Hemi_1SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_Hemi_1SecVtx_CutVtx_CutEvt_Merge = new TH1F("hData_Hemi_1SecVtx_CutVtx_CutEvt_Merge","",3,-1,1);

 TH1F* hData_CRfwd_Mmumu        = new TH1F("hData_CRfwd_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_BDTevt       = new TH1F("hData_CRfwd_BDTevt","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_ptmin   = new TH1F("hData_CRfwd_Hemi_ptmin","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_nTrks   = new TH1F("hData_CRfwd_Hemi_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_njetNOmu= new TH1F("hData_CRfwd_Hemi_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_mass    = new TH1F("hData_CRfwd_Hemi_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_massMu  = new TH1F("hData_CRfwd_Hemi_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_ptMu    = new TH1F("hData_CRfwd_Hemi_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_dRMu    = new TH1F("hData_CRfwd_Hemi_dRMu","",25,0.,7.5);
 TH1F* hData_CRfwd_CutEvt_Mmumu        = new TH1F("hData_CRfwd_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_CutEvt_ptmin   = new TH1F("hData_CRfwd_Hemi_CutEvt_ptmin","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_CutEvt_nTrks   = new TH1F("hData_CRfwd_Hemi_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_CutEvt_njetNOmu= new TH1F("hData_CRfwd_Hemi_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_CutEvt_mass    = new TH1F("hData_CRfwd_Hemi_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_CutEvt_massMu  = new TH1F("hData_CRfwd_Hemi_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_CutEvt_ptMu    = new TH1F("hData_CRfwd_Hemi_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_CutEvt_dRMu    = new TH1F("hData_CRfwd_Hemi_CutEvt_dRMu","",25,0.,7.5);

 TH1F* hData_CRfwd_Hemi_1Vtx_BDTevt  = new TH1F("hData_CRfwd_Hemi_1Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_1Vtx_BDTvtx  = new TH1F("hData_CRfwd_Hemi_1Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_1Vtx_njetNOmu= new TH1F("hData_CRfwd_Hemi_1Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_1Vtx_mass    = new TH1F("hData_CRfwd_Hemi_1Vtx_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_1Vtx_nTrks   = new TH1F("hData_CRfwd_Hemi_1Vtx_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_1Vtx_Mass    = new TH1F("hData_CRfwd_Hemi_1Vtx_Mass","",25,0.,100.);
 TH1F* hData_CRfwd_Hemi_1Vtx_massMu  = new TH1F("hData_CRfwd_Hemi_1Vtx_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_1Vtx_ptMu    = new TH1F("hData_CRfwd_Hemi_1Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_1Vtx_dRMu    = new TH1F("hData_CRfwd_Hemi_1Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRfwd_Hemi_1Vtx_Mmumu   = new TH1F("hData_CRfwd_Hemi_1Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_1Vtx_ptmin   = new TH1F("hData_CRfwd_Hemi_1Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRfwd_Hemi_2Vtx_BDTevt  = new TH1F("hData_CRfwd_Hemi_2Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_2Vtx_BDTvtx  = new TH1F("hData_CRfwd_Hemi_2Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_2Vtx_njetNOmu= new TH1F("hData_CRfwd_Hemi_2Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_2Vtx_mass    = new TH1F("hData_CRfwd_Hemi_2Vtx_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_2Vtx_nTrks   = new TH1F("hData_CRfwd_Hemi_2Vtx_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_2Vtx_Mass    = new TH1F("hData_CRfwd_Hemi_2Vtx_Mass","",25,0.,100.);
 TH1F* hData_CRfwd_Hemi_2Vtx_massMu  = new TH1F("hData_CRfwd_Hemi_2Vtx_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_2Vtx_ptMu    = new TH1F("hData_CRfwd_Hemi_2Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_2Vtx_dRMu    = new TH1F("hData_CRfwd_Hemi_2Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRfwd_Hemi_2Vtx_Mmumu   = new TH1F("hData_CRfwd_Hemi_2Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_2Vtx_ptmin	  = new TH1F("hData_CRfwd_Hemi_2Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_BDTevt  = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_BDTvtx  = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_njetNOmu= new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_mass    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_nTrks   = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_Mass    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_massMu  = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_ptMu    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_dRMu    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_Mmumu   = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutEvt_ptmin   = new TH1F("hData_CRfwd_Hemi_1Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_BDTevt  = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_BDTvtx  = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_njetNOmu= new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_mass    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_nTrks   = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_Mass    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_massMu  = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_ptMu    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_dRMu    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_Mmumu   = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutEvt_ptmin    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_BDTevt  = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_BDTvtx  = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_njetNOmu= new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_mass    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_nTrks   = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_Mass    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_massMu  = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_ptMu    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_dRMu    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_Mmumu   = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_ptmin   = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_BDTevt  = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_BDTvtx  = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_njetNOmu= new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_mass    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_nTrks   = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_Mass    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_massMu  = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_ptMu    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_dRMu    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_Mmumu   = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_ptmin   = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_ptmin   = new TH1F("hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_ptmin   = new TH1F("hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F*                       hData_CRfwd_Hemi_2SecVtx_CutEvt_Mass = new TH1F("hData_CRfwd_Hemi_2SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRfwd_Hemi_2SecVtx_CutVtx_CutEvt_Mass = new TH1F("hData_CRfwd_Hemi_2SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRfwd_Hemi_1SecVtx_CutEvt_Mass = new TH1F("hData_CRfwd_Hemi_1SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRfwd_Hemi_1SecVtx_CutEvt_Merge = new TH1F("hData_CRfwd_Hemi_1SecVtx_CutEvt_Merge","",3,-1,1);
 TH1F*                       hData_CRfwd_Hemi_1SecVtx_CutVtx_CutEvt_Mass = new TH1F("hData_CRfwd_Hemi_1SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRfwd_Hemi_1SecVtx_CutVtx_CutEvt_Merge = new TH1F("hData_CRfwd_Hemi_1SecVtx_CutVtx_CutEvt_Merge","",3,-1,1);

 TH1F* hData_CRlowpt_Mmumu        = new TH1F("hData_CRlowpt_Mmumu","",25,0.,500.);
 TH1F* hData_CRlowpt_BDTevt       = new TH1F("hData_CRlowpt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRlowpt_Hemi_ptmin   = new TH1F("hData_CRlowpt_Hemi_ptmin","",25,0.,1000.);
 TH1F* hData_CRlowpt_Hemi_nTrks   = new TH1F("hData_CRlowpt_Hemi_nTrks","",25,0.,50.);
 TH1F* hData_CRlowpt_Hemi_njetNOmu= new TH1F("hData_CRlowpt_Hemi_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRlowpt_Hemi_mass    = new TH1F("hData_CRlowpt_Hemi_mass","",25,0.,500.);
 TH1F* hData_CRlowpt_Hemi_massMu  = new TH1F("hData_CRlowpt_Hemi_massMu","",25,0.,750.);
 TH1F* hData_CRlowpt_Hemi_ptMu    = new TH1F("hData_CRlowpt_Hemi_ptMu","",25,0.,1000.);
 TH1F* hData_CRlowpt_Hemi_dRMu    = new TH1F("hData_CRlowpt_Hemi_dRMu","",25,0.,7.5);

 TH1F* hData_CRlowpt_Hemi_1Vtx_BDTevt  = new TH1F("hData_CRlowpt_Hemi_1Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_BDTvtx  = new TH1F("hData_CRlowpt_Hemi_1Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_njetNOmu= new TH1F("hData_CRlowpt_Hemi_1Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRlowpt_Hemi_1Vtx_mass    = new TH1F("hData_CRlowpt_Hemi_1Vtx_mass","",25,0.,500.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_nTrks   = new TH1F("hData_CRlowpt_Hemi_1Vtx_nTrks","",25,0.,50.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_Mass    = new TH1F("hData_CRlowpt_Hemi_1Vtx_Mass","",25,0.,100.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_massMu  = new TH1F("hData_CRlowpt_Hemi_1Vtx_massMu","",25,0.,750.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_ptMu    = new TH1F("hData_CRlowpt_Hemi_1Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_dRMu    = new TH1F("hData_CRlowpt_Hemi_1Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRlowpt_Hemi_1Vtx_Mmumu   = new TH1F("hData_CRlowpt_Hemi_1Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_ptmin   = new TH1F("hData_CRlowpt_Hemi_1Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRlowpt_Hemi_2Vtx_BDTevt  = new TH1F("hData_CRlowpt_Hemi_2Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_BDTvtx  = new TH1F("hData_CRlowpt_Hemi_2Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_njetNOmu= new TH1F("hData_CRlowpt_Hemi_2Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRlowpt_Hemi_2Vtx_mass    = new TH1F("hData_CRlowpt_Hemi_2Vtx_mass","",25,0.,500.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_nTrks   = new TH1F("hData_CRlowpt_Hemi_2Vtx_nTrks","",25,0.,50.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_Mass    = new TH1F("hData_CRlowpt_Hemi_2Vtx_Mass","",25,0.,100.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_massMu  = new TH1F("hData_CRlowpt_Hemi_2Vtx_massMu","",25,0.,750.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_ptMu    = new TH1F("hData_CRlowpt_Hemi_2Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_dRMu    = new TH1F("hData_CRlowpt_Hemi_2Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRlowpt_Hemi_2Vtx_Mmumu   = new TH1F("hData_CRlowpt_Hemi_2Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_ptmin   = new TH1F("hData_CRlowpt_Hemi_2Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_BDTevt  = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_BDTvtx  = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_njetNOmu= new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_mass    = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_nTrks   = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_Mass    = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_massMu  = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_ptMu    = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_dRMu    = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_Mmumu   = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRlowpt_Hemi_1Vtx_CutVtx_ptmin   = new TH1F("hData_CRlowpt_Hemi_1Vtx_CutVtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_BDTevt  = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_BDTvtx  = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_njetNOmu= new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_mass    = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_nTrks   = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_Mass    = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_massMu  = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_ptMu    = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_dRMu    = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_Mmumu   = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRlowpt_Hemi_2Vtx_CutVtx_ptmin   = new TH1F("hData_CRlowpt_Hemi_2Vtx_CutVtx_ptmin","",25,0.,1000.);
 TH1F*                        hData_CRlowpt_Hemi_2SecVtx_CutEvt_Mass = new TH1F("hData_CRlowpt_Hemi_2SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRlowpt_Hemi_2SecVtx_CutVtx_CutEvt_Mass = new TH1F("hData_CRlowpt_Hemi_2SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRlowpt_Hemi_1SecVtx_CutEvt_Mass = new TH1F("hData_CRlowpt_Hemi_1SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRlowpt_Hemi_1SecVtx_CutEvt_Merge = new TH1F("hData_CRlowpt_Hemi_1SecVtx_CutEvt_Merge","",3,-1,1);
 TH1F*                       hData_CRlowpt_Hemi_1SecVtx_CutVtx_CutEvt_Mass = new TH1F("hData_CRlowpt_Hemi_1SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRlowpt_Hemi_1SecVtx_CutVtx_CutEvt_Merge = new TH1F("hData_CRlowpt_Hemi_1SecVtx_CutVtx_CutEvt_Merge","",3,-1,1);


 TH1F* hData_CRloose_Hemi_0Vtx_BDTevt           = new TH1F("hData_CRloose_Hemi_0Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_0Vtx_nTrks            = new TH1F("hData_CRloose_Hemi_0Vtx_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_0Vtx_mass             = new TH1F("hData_CRloose_Hemi_0Vtx_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_0Vtx_Mmumu            = new TH1F("hData_CRloose_Hemi_0Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_0Vtx_njetNOmuAll      = new TH1F("hData_CRloose_Hemi_0Vtx_njetNOmuAll","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_0Vtx_ptmin            = new TH1F("hData_CRloose_Hemi_0Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRloose_Hemi_0Vtx_CutEvt_nTrks        = new TH1F("hData_CRloose_Hemi_0Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_0Vtx_CutEvt_mass         = new TH1F("hData_CRloose_Hemi_0Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_0Vtx_CutEvt_Mmumu        = new TH1F("hData_CRloose_Hemi_0Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_0Vtx_CutEvt_njetNOmuAll  = new TH1F("hData_CRloose_Hemi_0Vtx_CutEvt_njetNOmuAll","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_0Vtx_CutEvt_ptmin        = new TH1F("hData_CRloose_Hemi_0Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_CRloose_Hemi_1Vtx_BDTevt  = new TH1F("hData_CRloose_Hemi_1Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_1Vtx_BDTvtx  = new TH1F("hData_CRloose_Hemi_1Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_1Vtx_njetNOmu= new TH1F("hData_CRloose_Hemi_1Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_1Vtx_mass    = new TH1F("hData_CRloose_Hemi_1Vtx_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_1Vtx_nTrks   = new TH1F("hData_CRloose_Hemi_1Vtx_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_1Vtx_Mass    = new TH1F("hData_CRloose_Hemi_1Vtx_Mass","",25,0.,100.);
 TH1F* hData_CRloose_Hemi_1Vtx_massMu  = new TH1F("hData_CRloose_Hemi_1Vtx_massMu","",25,0.,750.);
 TH1F* hData_CRloose_Hemi_1Vtx_ptMu    = new TH1F("hData_CRloose_Hemi_1Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRloose_Hemi_1Vtx_dRMu    = new TH1F("hData_CRloose_Hemi_1Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRloose_Hemi_1Vtx_Mmumu   = new TH1F("hData_CRloose_Hemi_1Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_1Vtx_ptmin   = new TH1F("hData_CRloose_Hemi_1Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRloose_Hemi_2Vtx_BDTevt  = new TH1F("hData_CRloose_Hemi_2Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_2Vtx_BDTvtx  = new TH1F("hData_CRloose_Hemi_2Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_2Vtx_njetNOmu= new TH1F("hData_CRloose_Hemi_2Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_2Vtx_mass    = new TH1F("hData_CRloose_Hemi_2Vtx_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_2Vtx_nTrks   = new TH1F("hData_CRloose_Hemi_2Vtx_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_2Vtx_Mass    = new TH1F("hData_CRloose_Hemi_2Vtx_Mass","",25,0.,100.);
 TH1F* hData_CRloose_Hemi_2Vtx_massMu  = new TH1F("hData_CRloose_Hemi_2Vtx_massMu","",25,0.,750.);
 TH1F* hData_CRloose_Hemi_2Vtx_ptMu    = new TH1F("hData_CRloose_Hemi_2Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRloose_Hemi_2Vtx_dRMu    = new TH1F("hData_CRloose_Hemi_2Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRloose_Hemi_2Vtx_Mmumu   = new TH1F("hData_CRloose_Hemi_2Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_2Vtx_ptmin	  = new TH1F("hData_CRloose_Hemi_2Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_BDTevt  = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_BDTvtx  = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_njetNOmu= new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_mass    = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_nTrks   = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_Mass    = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_massMu  = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_ptMu    = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_dRMu    = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_Mmumu   = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutEvt_ptmin   = new TH1F("hData_CRloose_Hemi_1Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_BDTevt  = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_BDTvtx  = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_njetNOmu= new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_mass    = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_nTrks   = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_Mass    = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_massMu  = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_ptMu    = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_dRMu    = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_Mmumu   = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutEvt_ptmin    = new TH1F("hData_CRloose_Hemi_2Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_BDTevt  = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_BDTvtx  = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_njetNOmu= new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_mass    = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_nTrks   = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_Mass    = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_massMu  = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_ptMu    = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_dRMu    = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_Mmumu   = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_ptmin   = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_BDTevt  = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_BDTvtx  = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_njetNOmu= new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_mass    = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_nTrks   = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_Mass    = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_massMu  = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_ptMu    = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_dRMu    = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_Mmumu   = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_ptmin   = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_ptmin","",25,0.,1000.);

 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_ptmin   = new TH1F("hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_ptmin   = new TH1F("hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F*                      hData_CRloose_Hemi_2SecVtx_CutEvt_Mass = new TH1F("hData_CRloose_Hemi_2SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRloose_Hemi_2SecVtx_CutVtx_CutEvt_Mass = new TH1F("hData_CRloose_Hemi_2SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRloose_Hemi_1SecVtx_CutEvt_Mass = new TH1F("hData_CRloose_Hemi_1SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRloose_Hemi_1SecVtx_CutEvt_Merge = new TH1F("hData_CRloose_Hemi_1SecVtx_CutEvt_Merge","",3,-1,1);
 TH1F*                       hData_CRloose_Hemi_1SecVtx_CutVtx_CutEvt_Mass = new TH1F("hData_CRloose_Hemi_1SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*                       hData_CRloose_Hemi_1SecVtx_CutVtx_CutEvt_Merge = new TH1F("hData_CRloose_Hemi_1SecVtx_CutVtx_CutEvt_Merge","",3,-1,1);


 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_BDTevt  = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_BDTvtx  = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_njetNOmu= new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_mass    = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_nTrks   = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_Mass    = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_massMu  = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_ptMu    = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_dRMu    = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_Mmumu   = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutEvt_ptmin   = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_BDTevt  = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_BDTvtx  = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_njetNOmu= new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_mass    = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_nTrks   = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_Mass    = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_massMu  = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_ptMu    = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_dRMu    = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_Mmumu   = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_ptmin   = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_ptmin","",25,0.,1000.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutEvt_dRvtx   = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutEvt_dRvtx","",25,0.,5.);

 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_ptmin   = new TH1F("hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_ptmin   = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);
 TH1F* hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_dRvtx   = new TH1F("hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_dRvtx","",25,0.,5.);

 TH1F*   hData_CRlooselowpt_Hemi_2SecVtx_CutEvt_Mass         = new TH1F("hData_CRlooselowpt_Hemi_2SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*   hData_CRlooselowpt_Hemi_2SecVtx_CutVtx_CutEvt_Mass  = new TH1F("hData_CRlooselowpt_Hemi_2SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*   hData_CRlooselowpt_Hemi_1SecVtx_CutEvt_Mass         = new TH1F("hData_CRlooselowpt_Hemi_1SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*   hData_CRlooselowpt_Hemi_1SecVtx_CutEvt_Merge        = new TH1F("hData_CRlooselowpt_Hemi_1SecVtx_CutEvt_Merge","",3,-1,1);
 TH1F*   hData_CRlooselowpt_Hemi_1SecVtx_CutVtx_CutEvt_Mass  = new TH1F("hData_CRlooselowpt_Hemi_1SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*   hData_CRlooselowpt_Hemi_1SecVtx_CutVtx_CutEvt_Merge = new TH1F("hData_CRlooselowpt_Hemi_1SecVtx_CutVtx_CutEvt_Merge","",3,-1,1);


 TH1F* hData_SameSign_Filter     = new TH1F("hData_SameSign_Filter","",2,-0.5,1.5);
 TH1F* hData_SameSign_njetNOmuAll  = new TH1F("hData_SameSign_njetNOmuAll","",21,-0.5,20.5);
 TH1F* hData_SameSign_njetAll      = new TH1F("hData_SameSign_njetAll","",21,-0.5,20.5);
 TH1F* hData_SameSign_Hemi_eta     = new TH1F("hData_SameSign_Hemi_eta","",25,0.,7.5);

 TH1F* hData_SameSign_Mmumu        = new TH1F("hData_SameSign_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_BDTevt       = new TH1F("hData_SameSign_BDTevt","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_ptmin   = new TH1F("hData_SameSign_Hemi_ptmin","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_nTrks   = new TH1F("hData_SameSign_Hemi_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_njetNOmu= new TH1F("hData_SameSign_Hemi_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_mass    = new TH1F("hData_SameSign_Hemi_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_massMu  = new TH1F("hData_SameSign_Hemi_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_ptMu    = new TH1F("hData_SameSign_Hemi_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_dRMu    = new TH1F("hData_SameSign_Hemi_dRMu","",25,0.,7.5);
 TH1F* hData_SameSign_CutEvt_Mmumu        = new TH1F("hData_SameSign_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_CutEvt_ptmin   = new TH1F("hData_SameSign_Hemi_CutEvt_ptmin","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_CutEvt_nTrks   = new TH1F("hData_SameSign_Hemi_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_CutEvt_njetNOmu= new TH1F("hData_SameSign_Hemi_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_CutEvt_mass    = new TH1F("hData_SameSign_Hemi_CutEvt_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_CutEvt_massMu  = new TH1F("hData_SameSign_Hemi_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_CutEvt_ptMu    = new TH1F("hData_SameSign_Hemi_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_CutEvt_dRMu    = new TH1F("hData_SameSign_Hemi_CutEvt_dRMu","",25,0.,7.5);

 TH1F* hData_SameSign_Hemi_1Vtx_BDTevt  = new TH1F("hData_SameSign_Hemi_1Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_1Vtx_BDTvtx  = new TH1F("hData_SameSign_Hemi_1Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_1Vtx_njetNOmu= new TH1F("hData_SameSign_Hemi_1Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_1Vtx_mass    = new TH1F("hData_SameSign_Hemi_1Vtx_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_1Vtx_nTrks   = new TH1F("hData_SameSign_Hemi_1Vtx_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_1Vtx_Mass    = new TH1F("hData_SameSign_Hemi_1Vtx_Mass","",25,0.,100.);
 TH1F* hData_SameSign_Hemi_1Vtx_massMu  = new TH1F("hData_SameSign_Hemi_1Vtx_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_1Vtx_ptMu    = new TH1F("hData_SameSign_Hemi_1Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_1Vtx_dRMu    = new TH1F("hData_SameSign_Hemi_1Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_SameSign_Hemi_1Vtx_Mmumu   = new TH1F("hData_SameSign_Hemi_1Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_1Vtx_ptmin   = new TH1F("hData_SameSign_Hemi_1Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_SameSign_Hemi_2Vtx_BDTevt  = new TH1F("hData_SameSign_Hemi_2Vtx_BDTevt","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_2Vtx_BDTvtx  = new TH1F("hData_SameSign_Hemi_2Vtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_2Vtx_njetNOmu= new TH1F("hData_SameSign_Hemi_2Vtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_2Vtx_mass    = new TH1F("hData_SameSign_Hemi_2Vtx_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_2Vtx_nTrks   = new TH1F("hData_SameSign_Hemi_2Vtx_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_2Vtx_Mass    = new TH1F("hData_SameSign_Hemi_2Vtx_Mass","",25,0.,100.);
 TH1F* hData_SameSign_Hemi_2Vtx_massMu  = new TH1F("hData_SameSign_Hemi_2Vtx_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_2Vtx_ptMu    = new TH1F("hData_SameSign_Hemi_2Vtx_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_2Vtx_dRMu    = new TH1F("hData_SameSign_Hemi_2Vtx_dRMu","",25,0.,7.5);
 TH1F* hData_SameSign_Hemi_2Vtx_Mmumu   = new TH1F("hData_SameSign_Hemi_2Vtx_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_2Vtx_ptmin	= new TH1F("hData_SameSign_Hemi_2Vtx_ptmin","",25,0.,1000.);

 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_BDTevt  = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_BDTvtx  = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_njetNOmu= new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_mass    = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_nTrks   = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_Mass    = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_massMu  = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_ptMu    = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_dRMu    = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_Mmumu   = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutEvt_ptmin   = new TH1F("hData_SameSign_Hemi_1Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_BDTevt  = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_BDTvtx  = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_njetNOmu= new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_mass    = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_nTrks   = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_Mass    = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_massMu  = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_ptMu    = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_dRMu    = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_Mmumu   = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutEvt_ptmin   = new TH1F("hData_SameSign_Hemi_2Vtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_BDTevt  = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_BDTvtx  = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_njetNOmu= new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_mass    = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_nTrks   = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_Mass    = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_massMu  = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_ptMu    = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_dRMu    = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_Mmumu   = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_ptmin   = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_ptmin","",25,0.,1000.);

 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_BDTevt  = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_BDTevt","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_BDTvtx  = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_BDTvtx","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_njetNOmu= new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_mass    = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_nTrks   = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_Mass    = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_Mass","",25,0.,100.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_massMu  = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_ptMu    = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_dRMu    = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_dRMu","",25,0.,7.5);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_Mmumu   = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_ptmin   = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_ptmin","",25,0.,1000.);

 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_ptmin   = new TH1F("hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_BDTevt  = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_BDTevt","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx  = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx","",40,-1.,1.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu= new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu","",12,-0.5,11.5);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_mass    = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_mass","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_nTrks   = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_nTrks","",25,0.,50.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_Mass    = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_massMu  = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_massMu","",25,0.,750.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_ptMu    = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_ptMu","",25,0.,1000.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_dRMu    = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_dRMu","",25,0.,7.5);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_Mmumu   = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_Mmumu","",25,0.,500.);
 TH1F* hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_ptmin   = new TH1F("hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_ptmin","",25,0.,1000.);

 TH1F*   hData_SameSign_Hemi_2SecVtx_CutEvt_Mass             = new TH1F("hData_SameSign_Hemi_2SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*   hData_SameSign_Hemi_2SecVtx_CutVtx_CutEvt_Mass      = new TH1F("hData_SameSign_Hemi_2SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*   hData_SameSign_Hemi_1SecVtx_CutEvt_Mass             = new TH1F("hData_SameSign_Hemi_1SecVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*   hData_SameSign_Hemi_1SecVtx_CutEvt_Merge            = new TH1F("hData_SameSign_Hemi_1SecVtx_CutEvt_Merge","",3,-1,1);
 TH1F*   hData_SameSign_Hemi_1SecVtx_CutVtx_CutEvt_Mass      = new TH1F("hData_SameSign_Hemi_1SecVtx_CutVtx_CutEvt_Mass","",25,0.,100.);
 TH1F*   hData_SameSign_Hemi_1SecVtx_CutVtx_CutEvt_Merge     = new TH1F("hData_SameSign_Hemi_1SecVtx_CutVtx_CutEvt_Merge","",3,-1,1);
///////////////////////////////////////////////////////////////////

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

  float NormFactor = 1. ;
  float XS = 1;

// XS are given in pb
if (thesample.Contains("DYJetsToLL_M10to50"))                     { XS = 15910.0;   }
if (thesample.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays")) { XS = 10.8707;   }
if (thesample.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     { XS = 10.8908;   }
if (thesample.Contains("TTJets_DiLept"))                          { XS = 53.07;     }
if (thesample.Contains("WWTo2L2Nu_MLL_200To600"))                 { XS = 11.09;     }
if (thesample.Contains("WWTo2L2Nu_MLL_600To1200"))                { XS = 11.09;     }
if (thesample.Contains("WWTo2L2Nu"))                              { XS = 11.09;     }
if (thesample.Contains("WZTo2Q2L_mllmin4p0"))                     { XS = 6.535;     }
if (thesample.Contains("ZZTo2Q2L_mllmin4p0"))                     { XS = 3.676;     }
if (thesample.Contains("TTTo2L2Nu"))                              { XS = 88.3;      }
if (thesample.Contains("DYJetsToLL_M50"))                         { XS = 5379;      }
if (thesample.Contains("ttWJetsToLNu_5f_EWK"))                    { XS = 0.290;     } // not found on XSDB, no file on tier2...approximation
      //Took 0.868 pb (CMS-TOP-21-011)
     // as a starting point and then divided by 3 (lepton universality)
if (thesample.Contains("TTZToLL_5f"))                             { XS = 0.05188;   }//Not found on XSDB => used ana.py macro : 5.188e-02 +- 2.437e-04 pb
if (thesample.Contains("TTToHadronic"))                           { XS = 378.9;     }
if (thesample.Contains("TTWW"))                                   { XS = 0.006992;  }//found on XSDB
if (thesample.Contains("TTToSemiLeptonic") )                      { XS = 365.34;    }

//Signal
if (thesample.Contains("smu200")) { XS = 0.01;   }
if (thesample.Contains("smu250")) { XS = 0.0045; }
if (thesample.Contains("smu300")) { XS = 0.002;  }
if (thesample.Contains("smu350")) { XS = 0.001;  }
if (thesample.Contains("smu400")) { XS = 0.0006; }
if (thesample.Contains("smu450")) { XS = 0.0004; }
if (thesample.Contains("smu500")) { XS = 0.00025;}
if (thesample.Contains("70_test")){ XS = 0.0035; } //this is for 14tev
if (thesample.Contains("50_test")){ XS = 0.0025; } //this is for 14tev 
if (thesample.Contains("30_test")){ XS = 0.002;  } //this is for 14tev 
if (thesample.Contains("10_test")){ XS = 0.0035; } //this is for 14tev
// if (thesample.Contains("smu300_250_240")){ XS = 0.002;  }

  float lumiRun2 = 136000 ; //pb-1
  float Nevent = 0;
//   if (nentries>10000000){nentries = 10000000;}
   // NormFactor =  XS/nentries;

   cout<< "Line : "  << __LINE__ << " " << nentries << endl; 

   bool BlindSR = false;
   int   n_EVT_SR_1Vtx = 0;
   int   n_EVT_CRNoEvt_1Vtx = 0;
   int   n_EVT_CRNoEvtNoVtx_1Vtx = 0;
   int   n_EVT_CRNoVtx_1Vtx = 0;
   int   n_EVT_SR_2Vtx = 0;
   int   n_EVT_CRNoEvt_2Vtx = 0;
   int   n_EVT_CRNoEvtNoVtx_2Vtx = 0;
   int  n_EVT_CRNoVtx_2Vtx = 0;

   int   n_EVT2_SR_1Vtx = 0;
   int   n_EVT2_CRNoEvt_1Vtx = 0;
   int   n_EVT2_CRNoEvtNoVtx_1Vtx = 0;
   int   n_EVT2_CRNoVtx_1Vtx = 0;
   int   n_EVT2_SR_2Vtx = 0;
   int   n_EVT2_CRNoEvt_2Vtx = 0;
   int   n_EVT2_CRNoEvtNoVtx_2Vtx = 0;
   int  n_EVT2_CRNoVtx_2Vtx = 0;

   float EVTSWP = 0;
   float VTXWP = 0.5;

   Long64_t nbytes = 0, nb = 0;
   int TempN = 1000000;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      Nevent++;
      auto tim = std::chrono::system_clock::now();
      std::time_t start = std::chrono::system_clock::to_time_t(tim);
      if (jentry==0.1*nentries) {std::cout<<"10/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.2*nentries) {std::cout<<"20/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.3*nentries) {std::cout<<"30/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.4*nentries) {std::cout<<"40/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.5*nentries) {std::cout<<"50/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.6*nentries) {std::cout<<"60/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.7*nentries) {std::cout<<"70/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.8*nentries) {std::cout<<"80/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.9*nentries) {std::cout<<"90/100 :"<<std::ctime(&start)<<std::endl;}
      if ( jentry%1000 == 0 ) std::cout << "events : " << jentry << std::endl;
      // if ( jentry >= 1000000 ) break;

      float TightWP = 0.85;// for tracks => everytime Tight is mentioned, it is a refrence to this value.
      float LooseWP = 0.; // for tracks => same for loose





      float weight = 1.;
      hSim_weight->Fill( weight );
   
      hData_Filter->Fill( tree_Filter );
      if ( tree_Filter ) hData_njetNOmuAll->Fill( tree_njetNOmu );
      if ( tree_Filter && tree_njetNOmu > 0 ) 
         {
            hData_njetAll->Fill( tree_njet );
            hData_Hemi_eta->Fill( abs(tree_Hemi_eta->at(0)) );
            hData_Hemi_eta->Fill( abs(tree_Hemi_eta->at(1)) );
         }

      hData_SameSign_Filter->Fill( tree_FilterSameSign );
   
      if ( tree_FilterSameSign && !tree_Filter ) hData_SameSign_njetNOmuAll->Fill( tree_njetNOmu );
      if ( tree_FilterSameSign && !tree_Filter && tree_njetNOmu > 0 ) 
         {
            hData_SameSign_njetAll->Fill( tree_njet );
            hData_SameSign_Hemi_eta->Fill( abs(tree_Hemi_eta->at(0)) );
            hData_SameSign_Hemi_eta->Fill( abs(tree_Hemi_eta->at(1)) );
         }

      //$$$$
      if ( !tree_Filter && !tree_FilterSameSign ) continue;
      if ( tree_njetNOmu == 0 ) continue;
      // if ( abs(tree_Hemi_eta->at(0)) > 2.5 || abs(tree_Hemi_eta->at(1)) > 2.5 ) continue;
      //$$$$

      // ABCD //
      hData_MmumuABCD->Fill(tree_Mmumu);
      hData_BDTevtABCD->Fill(tree_Evts_MVAval);
      // hData_BDTevtDY->Fill(tree_Evts_MVAvalDY);
      // hData_BDTevtTT->Fill(tree_Evts_MVAvalTT);

      // hData_BDTHemi1->Fill(tree_Hemi1_MVAval);
      // hData_BDTHemi1DY->Fill(tree_Hemi1_MVAvalDY);
      // hData_BDTHemi1TT->Fill(tree_Hemi1_MVAvalTT);

      // hData_BDTHemi2->Fill(tree_Hemi2_MVAval);
      // hData_BDTHemi2DY->Fill(tree_Hemi2_MVAvalDY);
      // hData_BDTHemi2TT->Fill(tree_Hemi2_MVAvalTT);

      // ABCD - 1  Using Event BDT and the Steps of reconstruction of the BDT 1+2 or 3+4 //
      // ABCD - 2  Using Both Event and Vtx BDT //
      // 1Vtx Region : look at Mmumu mass, mass vertex, hemi mass

      if ((tree_Hemi_Vtx_NChi2->at(0) > 0. && tree_Hemi_Vtx_NChi2->at(0) < 10.)
      || (tree_Hemi_Vtx_NChi2->at(1) > 0. && tree_Hemi_Vtx_NChi2->at(1) < 10.))
            {

               int Vtxidx = 0 ;
               if (tree_Hemi_Vtx_NChi2->at(1) > 0. && tree_Hemi_Vtx_NChi2->at(1) < 10.) Vtxidx = 1;
               //ABCD - 1
               // Signal Region  : EVT BDT > WP && Tight 1+2
               if ( tree_Evts_MVAval > EVTSWP &&  (tree_Hemi_Vtx_step->at(Vtxidx) == 1 || tree_Hemi_Vtx_step->at(Vtxidx) == 2) ) 
                  {
                     hData2_SR_1Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(Vtxidx));
                     hData2_SR_1Vtx_MHemi->Fill(tree_Hemi_mass->at(Vtxidx));
                     hData2_SR_1Vtx_Mmumu->Fill(tree_Mmumu);
                     n_EVT2_SR_1Vtx++;
                  }
                // CR : EVT BDT < WP  && Tight 1+2
               if ( tree_Evts_MVAval < EVTSWP && (tree_Hemi_Vtx_step->at(Vtxidx) == 1 || tree_Hemi_Vtx_step->at(Vtxidx) == 2) )
                  {
                     hData2_CRNoEvt_1Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(Vtxidx));
                     hData2_CRNoEvt_1Vtx_MHemi->Fill(tree_Hemi_mass->at(Vtxidx));
                     hData2_CRNoEvt_1Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT2_CRNoEvt_1Vtx++;
                  }
               // CR : EVT BDT < WP && Loose 3+4
               if ( tree_Evts_MVAval < EVTSWP && (tree_Hemi_Vtx_step->at(Vtxidx) == 3 || tree_Hemi_Vtx_step->at(Vtxidx) == 4)  ) 
                  {
                     hData2_CRNoEvtNoVtx_1Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(Vtxidx));
                     hData2_CRNoEvtNoVtx_1Vtx_MHemi->Fill(tree_Hemi_mass->at(Vtxidx));
                     hData2_CRNoEvtNoVtx_1Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT2_CRNoEvtNoVtx_1Vtx++;
                  }
                // CR : EVT BDT > WP && Loose 3+4
               if ( tree_Evts_MVAval > EVTSWP &&  (tree_Hemi_Vtx_step->at(Vtxidx) == 3 || tree_Hemi_Vtx_step->at(Vtxidx) == 4)  ) 
                  {
                        // tree_Hemi_Vtx_Mass->at(iVtx)
                     hData2_CRNoVtx_1Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(Vtxidx));
                     hData2_CRNoVtx_1Vtx_MHemi->Fill(tree_Hemi_mass->at(Vtxidx));
                     hData2_CRNoVtx_1Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT2_CRNoVtx_1Vtx++;
                  }
              
              //ABCD - 2
              if (tree_Hemi_Vtx_step->at(Vtxidx) < 1 || tree_Hemi_Vtx_step->at(Vtxidx) > 2) continue;
               // Signal Region  : EVT BDT > WP && VTX > WP
               if ( tree_Evts_MVAval > EVTSWP &&  tree_Hemi_Vtx_MVAval_Step1->at(Vtxidx) > VTXWP ) 
                  {
                     hData_SR_1Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(Vtxidx));
                     hData_SR_1Vtx_MHemi->Fill(tree_Hemi_mass->at(Vtxidx));
                     hData_SR_1Vtx_Mmumu->Fill(tree_Mmumu);
                     n_EVT_SR_1Vtx++;
                  }
                // CR : EVT BDT < WP && VTX > WP
               if ( tree_Evts_MVAval < EVTSWP && tree_Hemi_Vtx_MVAval_Step1->at(Vtxidx) > VTXWP )
                  {
                     hData_CRNoEvt_1Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(Vtxidx));
                     hData_CRNoEvt_1Vtx_MHemi->Fill(tree_Hemi_mass->at(Vtxidx));
                     hData_CRNoEvt_1Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT_CRNoEvt_1Vtx++;
                  }
               // CR : EVT BDT < WP && VTX < WP
               if ( tree_Evts_MVAval < EVTSWP && tree_Hemi_Vtx_MVAval_Step1->at(Vtxidx) < VTXWP  ) 
                  {
                     hData_CRNoEvtNoVtx_1Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(Vtxidx));
                     hData_CRNoEvtNoVtx_1Vtx_MHemi->Fill(tree_Hemi_mass->at(Vtxidx));
                     hData_CRNoEvtNoVtx_1Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT_CRNoEvtNoVtx_1Vtx++;
                  }
                // CR : EVT BDT > WP && VTX < WP
               if ( tree_Evts_MVAval > EVTSWP &&  tree_Hemi_Vtx_MVAval_Step1->at(Vtxidx) < VTXWP  ) 
                  {
                        // tree_Hemi_Vtx_Mass->at(iVtx)
                     hData_CRNoVtx_1Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(Vtxidx));
                     hData_CRNoVtx_1Vtx_MHemi->Fill(tree_Hemi_mass->at(Vtxidx));
                     hData_CRNoVtx_1Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT_CRNoVtx_1Vtx++;
                  }
            }



      // 2 Vtx Region : look at Mmumu mass and mass vertex (highest one)
      if ((tree_Hemi_Vtx_NChi2->at(0) > 0. && tree_Hemi_Vtx_NChi2->at(0) < 10.)
      && (tree_Hemi_Vtx_NChi2->at(1) > 0. && tree_Hemi_Vtx_NChi2->at(1) < 10.))
            {


               int MVtxmax = 0;
               int HVtxmax = 0;

               if (tree_Hemi_Vtx_Mass->at(1) > tree_Hemi_Vtx_Mass->at(0)) MVtxmax = 1;
               if (tree_Hemi_mass->at(1) > tree_Hemi_mass->at(0)) HVtxmax = 1;

               // ABCD - 1
               // Signal Region  : EVT BDT > WP && Tight 1+2

               if ( tree_Evts_MVAval > EVTSWP  && 
               (tree_Hemi_Vtx_step->at(0) == 1 || tree_Hemi_Vtx_step->at(0) == 2) && (tree_Hemi_Vtx_step->at(1) == 1 || tree_Hemi_Vtx_step->at(1) == 2) ) 
                  {
                     hData2_SR_2Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(MVtxmax));
                     hData2_SR_2Vtx_MHemi->Fill(tree_Hemi_mass->at(HVtxmax));
                     hData2_SR_2Vtx_Mmumu->Fill(tree_Mmumu);
                     n_EVT2_SR_2Vtx++;
                  }
                // CR : EVT BDT < WP && Tight 1+2
               if ( tree_Evts_MVAval < EVTSWP && 
               (tree_Hemi_Vtx_step->at(0) == 1 || tree_Hemi_Vtx_step->at(0) == 2) && (tree_Hemi_Vtx_step->at(1) == 1 || tree_Hemi_Vtx_step->at(1) == 2) ) 
                  {
                     hData2_CRNoEvt_2Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(MVtxmax));
                     hData2_CRNoEvt_2Vtx_MHemi->Fill(tree_Hemi_mass->at(HVtxmax));
                     hData2_CRNoEvt_2Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT2_CRNoEvt_2Vtx++;
                  }
               // CR : EVT BDT < WP && Loose 3+4
               if ( tree_Evts_MVAval < EVTSWP && 
               (tree_Hemi_Vtx_step->at(0) == 3 || tree_Hemi_Vtx_step->at(0) == 4) && (tree_Hemi_Vtx_step->at(1) ==3 || tree_Hemi_Vtx_step->at(1) == 4) ) 
                  {
                     hData2_CRNoEvtNoVtx_2Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(MVtxmax));
                     hData2_CRNoEvtNoVtx_2Vtx_MHemi->Fill(tree_Hemi_mass->at(HVtxmax));
                     hData2_CRNoEvtNoVtx_2Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT2_CRNoEvtNoVtx_2Vtx++;
                  }
                // CR : EVT BDT > WP && Loose 3+4
               if ( tree_Evts_MVAval > EVTSWP && 
               (tree_Hemi_Vtx_step->at(0) == 3 || tree_Hemi_Vtx_step->at(0) == 4) && (tree_Hemi_Vtx_step->at(1) ==3 || tree_Hemi_Vtx_step->at(1) == 4) ) 
                  {
                     hData2_CRNoVtx_2Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(MVtxmax));
                     hData2_CRNoVtx_2Vtx_MHemi->Fill(tree_Hemi_mass->at(HVtxmax));
                     hData2_CRNoVtx_2Vtx_Mmumu->Fill(tree_Mmumu);    
                     n_EVT2_CRNoVtx_2Vtx++;
                  }

               // ABCD - 2
               // Signal Region  : EVT BDT > WP && VTX > WP
               if ( (tree_Hemi_Vtx_step->at(0) < 1 || tree_Hemi_Vtx_step->at(0) > 2) 
               &&   (tree_Hemi_Vtx_step->at(1) < 1 || tree_Hemi_Vtx_step->at(1) > 2)) continue;

               if ( tree_Evts_MVAval > EVTSWP  && 
               ( tree_Hemi_Vtx_MVAval_Step1->at(0) > VTXWP && tree_Hemi_Vtx_MVAval_Step1->at(1) > VTXWP ) ) 
                  {
                     hData_SR_2Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(MVtxmax));
                     hData_SR_2Vtx_MHemi->Fill(tree_Hemi_mass->at(HVtxmax));
                     hData_SR_2Vtx_Mmumu->Fill(tree_Mmumu);
                     n_EVT_SR_2Vtx++;
                  }
                // CR : EVT BDT < WP && VTX > WP
               if ( tree_Evts_MVAval < EVTSWP && 
               ( tree_Hemi_Vtx_MVAval_Step1->at(0) > VTXWP && tree_Hemi_Vtx_MVAval_Step1->at(1) > VTXWP ) ) 
                  {
                     hData_CRNoEvt_2Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(MVtxmax));
                     hData_CRNoEvt_2Vtx_MHemi->Fill(tree_Hemi_mass->at(HVtxmax));
                     hData_CRNoEvt_2Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT_CRNoEvt_2Vtx++;
                  }
               // CR : EVT BDT < WP && VTX < WP
               if ( tree_Evts_MVAval < EVTSWP && 
               ( tree_Hemi_Vtx_MVAval_Step1->at(0) < VTXWP && tree_Hemi_Vtx_MVAval_Step1->at(1) < VTXWP ) ) 
                  {
                     hData_CRNoEvtNoVtx_2Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(MVtxmax));
                     hData_CRNoEvtNoVtx_2Vtx_MHemi->Fill(tree_Hemi_mass->at(HVtxmax));
                     hData_CRNoEvtNoVtx_2Vtx_Mmumu->Fill(tree_Mmumu); 
                     n_EVT_CRNoEvtNoVtx_2Vtx++;
                  }
                // CR : EVT BDT > WP && VTX < WP
               if ( tree_Evts_MVAval > EVTSWP && 
               ( tree_Hemi_Vtx_MVAval_Step1->at(0) < VTXWP && tree_Hemi_Vtx_MVAval_Step1->at(1) < VTXWP ) ) 
                  {
                     hData_CRNoVtx_2Vtx_MVtx->Fill(tree_Hemi_Vtx_Mass->at(MVtxmax));
                     hData_CRNoVtx_2Vtx_MHemi->Fill(tree_Hemi_mass->at(HVtxmax));
                     hData_CRNoVtx_2Vtx_Mmumu->Fill(tree_Mmumu);    
                     n_EVT_CRNoVtx_2Vtx++;
                  }
            }
      // END OF ABCD //

      // Daniel's Method using hemispheres
      bool isHemiVtx1 = false, isHemiVtx2 = false;
      bool isCutVtx = false, isCutVtx1 = false, isCutVtx2 = false;
      bool isCutEvt = false;
      int nVtx = 0, nTrks = 0;
      int njet = 0, njetNOmu = 0;
      float BDTvtx = -2., BDTvtx1 = -2., BDTvtx2 = -2.;
      float mass = -1., VtxMass = -1.;
      float massMu = -1., ptMu = -1., dRMu = 10.;

      nTrks  = tree_Hemi_nTrks->at(0);
      mass   = tree_Hemi_mass->at(0);
      massMu = tree_HemiMu_mass->at(0);
      ptMu   = tree_HemiMu_pt->at(0);
      dRMu   = tree_HemiMu_dR->at(0);
      njetNOmu = tree_Hemi_njet_nomu->at(0);

      if ( tree_Hemi_nTrks->at(1) < nTrks )    nTrks = tree_Hemi_nTrks->at(1);
      if ( tree_Hemi_mass->at(1)	< mass )     mass = tree_Hemi_mass->at(1);
      if ( tree_HemiMu_mass->at(1) < massMu ) massMu = tree_HemiMu_mass->at(1);
      if ( tree_HemiMu_pt->at(1)	< ptMu )     ptMu = tree_HemiMu_pt->at(1);
      if ( tree_HemiMu_dR->at(1)	> dRMu )     dRMu = tree_HemiMu_dR->at(1);
      if ( tree_Hemi_njet_nomu->at(1) < njetNOmu ) njetNOmu = tree_Hemi_njet_nomu->at(1);

      ////////////////// hemisphere pT
      float  hemi1_pt  = -1.;
      float  hemi2_pt  = -1.;
      ////////////////////////////////

      hemi1_pt = tree_Hemi_pt->at(0);
      hemi2_pt = tree_Hemi_pt->at(1);

      float hemi_ptmin = hemi2_pt;
      if ( hemi2_pt > hemi1_pt ) hemi_ptmin = hemi1_pt;

      //$$$$
      // Signal Region
      if(!BlindSR)
         {
            if ( tree_Filter &&
            abs(tree_Hemi_eta->at(0)) < 2.4 && abs(tree_Hemi_eta->at(1)) < 2.4 ) 
               {

                  hData_Mmumu->Fill( tree_Mmumu );
                  hData_BDTevt->Fill( tree_Evts_MVAval );
   
                  if ( tree_Hemi_Vtx_NChi2->at(0) > 0 && tree_Hemi_Vtx_NChi2->at(0) < 10
                  && tree_Hemi_Vtx_step->at(0) >= 1 && tree_Hemi_Vtx_step->at(0) <= 2 
                  )  {
                        isHemiVtx1 = true;
                        BDTvtx1 = tree_Hemi_Vtx_MVAval_Step1->at(0);
                        BDTvtx  = BDTvtx1;
                        VtxMass = tree_Hemi_Vtx_Mass->at(0);
                     }

                  if ( tree_Hemi_Vtx_NChi2->at(1) > 0 && tree_Hemi_Vtx_NChi2->at(1) < 10
                     && tree_Hemi_Vtx_step->at(1) >= 1 && tree_Hemi_Vtx_step->at(1) <= 2 
                   ) { 
                        isHemiVtx2 = true;
                        BDTvtx2 = tree_Hemi_Vtx_MVAval_Step1->at(1);
                        if ( BDTvtx2 > BDTvtx ) BDTvtx = BDTvtx2;
                        if ( tree_Hemi_Vtx_Mass->at(1)  > VtxMass ) VtxMass = tree_Hemi_Vtx_Mass->at(1);
                     }

                  if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;
                  else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;
   
                  //$$
                  if ( VtxMass > 8. ) isCutVtx = true; 
                  //$$

                  // hData_Hemi_lead_ptmin->Fill( hemi_lead_ptmin );
                  hData_Hemi_ptmin->Fill( hemi_ptmin );
                  hData_Hemi_njetNOmu->Fill(	njetNOmu );
                  hData_Hemi_nTrks->Fill(  nTrks );
                  hData_Hemi_mass->Fill(   mass );
                  hData_Hemi_massMu->Fill( massMu );
                  hData_Hemi_ptMu->Fill(   ptMu );
                  hData_Hemi_dRMu->Fill(   dRMu );

                  //$$
                     if ( hemi_ptmin > 80. ) isCutEvt = true;   
                  //$$

                  if ( isCutEvt ) {
                  hData_CutEvt_Mmumu->Fill( tree_Mmumu );
                  hData_Hemi_CutEvt_njetNOmu->Fill(	njetNOmu );
                  hData_Hemi_CutEvt_nTrks->Fill(  nTrks );
                  hData_Hemi_CutEvt_mass->Fill(   mass );
                  hData_Hemi_CutEvt_massMu->Fill( massMu );
                  hData_Hemi_CutEvt_ptMu->Fill(   ptMu );
                  hData_Hemi_CutEvt_dRMu->Fill(   dRMu );
                  // hData_Hemi_CutEvt_lead_ptmin->Fill( hemi_lead_ptmin );
                  hData_Hemi_CutEvt_ptmin->Fill( hemi_ptmin );
                  }
                  
                  if ( nVtx == 0 ) {
                  hData_Hemi_0Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  hData_Hemi_0Vtx_nTrks->Fill(  nTrks );
                  hData_Hemi_0Vtx_mass->Fill(   mass );
                  hData_Hemi_0Vtx_Mmumu->Fill(  tree_Mmumu );
                  hData_Hemi_0Vtx_PFMet->Fill(  tree_PFMet_et );
                  hData_Hemi_0Vtx_njetNOmuAll->Fill(   tree_njetNOmu );
                  // hData_Hemi_0Vtx_lead_ptmin->Fill( hemi_lead_ptmin );
                  hData_Hemi_0Vtx_ptmin->Fill( hemi_ptmin );
                  }
                  if ( nVtx == 0 && isCutEvt ) {
                     hData_Hemi_0Vtx_CutEvt_nTrks->Fill(  nTrks );
                     hData_Hemi_0Vtx_CutEvt_mass->Fill(   mass );
                     hData_Hemi_0Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                     hData_Hemi_0Vtx_CutEvt_njetNOmuAll->Fill(   tree_njetNOmu );
                     // hData_Hemi_0Vtx_CutEvt_lead_ptmin->Fill( hemi_lead_ptmin );
                     hData_Hemi_0Vtx_CutEvt_ptmin->Fill( hemi_ptmin );
                     }

                  if ( nVtx == 1 ) {
                  hData_Hemi_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  hData_Hemi_1Vtx_BDTvtx->Fill( BDTvtx );
                  hData_Hemi_1Vtx_njetNOmu->Fill(   njetNOmu );
                  hData_Hemi_1Vtx_nTrks->Fill(  nTrks );
                  hData_Hemi_1Vtx_mass->Fill(   mass );
                  hData_Hemi_1Vtx_Mass->Fill(   VtxMass );
                  hData_Hemi_1Vtx_massMu->Fill( massMu );
                  hData_Hemi_1Vtx_ptMu->Fill(   ptMu );
                  hData_Hemi_1Vtx_dRMu->Fill(   dRMu );
                  hData_Hemi_1Vtx_Mmumu->Fill(  tree_Mmumu );
                  hData_Hemi_1Vtx_PFMet->Fill(  tree_PFMet_et );
                  hData_Hemi_1Vtx_njetNOmuAll->Fill( tree_njetNOmu );
                  // hData_Hemi_1Vtx_lead_ptmin->Fill(  hemi_lead_ptmin );
                  hData_Hemi_1Vtx_ptmin->Fill(  hemi_ptmin );
                  }
                  
                  if ( nVtx == 2 ) {
                  hData_Hemi_2Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  hData_Hemi_2Vtx_BDTvtx->Fill( BDTvtx );
                  hData_Hemi_2Vtx_njetNOmu->Fill(   njetNOmu );
                  hData_Hemi_2Vtx_nTrks->Fill(  nTrks );
                  hData_Hemi_2Vtx_mass->Fill(   mass );
                  hData_Hemi_2Vtx_Mass->Fill(   VtxMass );
                  hData_Hemi_2Vtx_massMu->Fill( massMu );
                  hData_Hemi_2Vtx_ptMu->Fill(   ptMu );
                  hData_Hemi_2Vtx_dRMu->Fill(   dRMu );
                  hData_Hemi_2Vtx_Mmumu->Fill(  tree_Mmumu );
                  hData_Hemi_2Vtx_PFMet->Fill(  tree_PFMet_et );
                  // hData_Hemi_2Vtx_lead_ptmin->Fill(  hemi_lead_ptmin );
                  hData_Hemi_2Vtx_ptmin->Fill(  hemi_ptmin );
                  hData_Hemi_2VtxAll_BDTvtx->Fill(   tree_Hemi_Vtx_MVAval_Step1->at(0) );
                  hData_Hemi_2VtxAll_BDTvtx->Fill(   tree_Hemi_Vtx_MVAval_Step1->at(1) );
                  hData_Hemi_2VtxAll_njetNOmu->Fill( tree_Hemi_njet_nomu->at(0) );
                  hData_Hemi_2VtxAll_njetNOmu->Fill( tree_Hemi_njet_nomu->at(1) );
                  hData_Hemi_2VtxAll_nTrks->Fill(    tree_Hemi_nTrks->at(0) );
                  hData_Hemi_2VtxAll_nTrks->Fill(    tree_Hemi_nTrks->at(1) );
                  hData_Hemi_2VtxAll_mass->Fill(     tree_Hemi_mass->at(0) );
                  hData_Hemi_2VtxAll_mass->Fill(     tree_Hemi_mass->at(1) );
                  hData_Hemi_2VtxAll_Mass->Fill(     tree_Hemi_Vtx_Mass->at(0) );
                  hData_Hemi_2VtxAll_Mass->Fill(     tree_Hemi_Vtx_Mass->at(1) );
                  hData_Hemi_2VtxAll_massMu->Fill(   tree_HemiMu_mass->at(0) );
                  hData_Hemi_2VtxAll_massMu->Fill(   tree_HemiMu_mass->at(1) );
                  hData_Hemi_2VtxAll_ptMu->Fill(     tree_HemiMu_pt->at(0) );
                  hData_Hemi_2VtxAll_ptMu->Fill(     tree_HemiMu_pt->at(1) );
                  hData_Hemi_2VtxAll_dRMu->Fill(     tree_HemiMu_dR->at(0) );
                  hData_Hemi_2VtxAll_dRMu->Fill(     tree_HemiMu_dR->at(1) );

                                          //------------------------
            int NewnVtx = -1;
            if (tree_Hemi_Merging->at(0) && tree_Hemi_Merging->at(1)) {NewnVtx = 2;}
            else if (tree_Hemi_Merging->at(0)==0 && tree_Hemi_Merging->at(1) == 0){NewnVtx = 0;}
            else {NewnVtx = 1;}
            float NewVtxMass = -10;
             bool SecTight = false;
            if (NewnVtx == 2)
               {
                  if (tree_Hemi_SecVtx_Mass->at(0)<tree_Hemi_SecVtx_Mass->at(1)) NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                  bool NewisCutVtx = false;
                  if (tree_Hemi_SecVtx_isTight->at(0) == true && tree_Hemi_SecVtx_isTight->at(1) == true) SecTight = true;
                  if (NewVtxMass > 8) NewisCutVtx = true;
                  if (  isCutEvt && SecTight ) 
                     {
                        hData_Hemi_2SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                     }
                  if (NewisCutVtx && isCutEvt)
                     {
                        hData_Hemi_2SecVtx_CutVtx_CutEvt_Mass->Fill(   NewVtxMass );
                     }
               }
            if (NewnVtx == 1)
               {  
                  bool Merge = false ;
                  SecTight = tree_Hemi_SecVtx_isTight->at(1);
                  if (tree_Hemi_Merging->at(0)) Merge = true;
                  
                  if (Merge) {NewVtxMass = tree_Hemi_SecVtx_Mass->at(0);SecTight = tree_Hemi_SecVtx_isTight->at(0);}
                  else NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                  bool NewisCutVtx = false;
                  if (NewVtxMass > 8) NewisCutVtx = true;
                  if (  isCutEvt && SecTight ) 
                     {
                        hData_Hemi_1SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                        hData_Hemi_1SecVtx_CutEvt_Merge->Fill(Merge);
                     }
                  if (NewisCutVtx && isCutEvt && SecTight)
                     {
                        hData_Hemi_1SecVtx_CutVtx_CutEvt_Mass->Fill(  NewVtxMass );
                        hData_Hemi_1SecVtx_CutVtx_CutEvt_Merge->Fill(Merge);
                     }
               }

            //--------------------------------------
                  }

                  if ( nVtx == 1 && isCutEvt ) {
                  hData_Hemi_1Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
                  hData_Hemi_1Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
                  hData_Hemi_1Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
                  hData_Hemi_1Vtx_CutEvt_nTrks->Fill(  nTrks );
                  hData_Hemi_1Vtx_CutEvt_mass->Fill(   mass );
                  hData_Hemi_1Vtx_CutEvt_Mass->Fill(   VtxMass );
                  hData_Hemi_1Vtx_CutEvt_massMu->Fill( massMu );
                  hData_Hemi_1Vtx_CutEvt_ptMu->Fill(   ptMu );
                  hData_Hemi_1Vtx_CutEvt_dRMu->Fill(   dRMu );
                  hData_Hemi_1Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_Hemi_1Vtx_CutEvt_njetNOmuAll->Fill( tree_njetNOmu );
                  // hData_Hemi_1Vtx_CutEvt_lead_ptmin->Fill(  hemi_lead_ptmin );
                  hData_Hemi_1Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
                  }

                  if ( nVtx == 2 && isCutEvt ) {
                  hData_Hemi_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
                  hData_Hemi_2Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
                  hData_Hemi_2Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
                  hData_Hemi_2Vtx_CutEvt_nTrks->Fill(  nTrks );
                  hData_Hemi_2Vtx_CutEvt_mass->Fill(   mass );
                  hData_Hemi_2Vtx_CutEvt_Mass->Fill(   VtxMass );
                  hData_Hemi_2Vtx_CutEvt_massMu->Fill( massMu );
                  hData_Hemi_2Vtx_CutEvt_ptMu->Fill(   ptMu );
                  hData_Hemi_2Vtx_CutEvt_dRMu->Fill(   dRMu );
                  hData_Hemi_2Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_Hemi_2Vtx_CutEvt_njetNOmuAll->Fill( tree_njetNOmu );
                  // hData_Hemi_2Vtx_CutEvt_lead_ptmin->Fill(  hemi_lead_ptmin );
                  hData_Hemi_2Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
                  hData_Hemi_2VtxAll_CutEvt_BDTvtx->Fill(   tree_Hemi_Vtx_MVAval_Step1->at(0) );
                  hData_Hemi_2VtxAll_CutEvt_BDTvtx->Fill(   tree_Hemi_Vtx_MVAval_Step1->at(1) );
                  hData_Hemi_2VtxAll_CutEvt_njetNOmu->Fill( tree_Hemi_njet_nomu->at(0) );
                  hData_Hemi_2VtxAll_CutEvt_njetNOmu->Fill( tree_Hemi_njet_nomu->at(1) );
                  hData_Hemi_2VtxAll_CutEvt_nTrks->Fill(    tree_Hemi_nTrks->at(0) );
                  hData_Hemi_2VtxAll_CutEvt_nTrks->Fill(    tree_Hemi_nTrks->at(1) );
                  hData_Hemi_2VtxAll_CutEvt_mass->Fill(     tree_Hemi_mass->at(0) );
                  hData_Hemi_2VtxAll_CutEvt_mass->Fill(     tree_Hemi_mass->at(1) );
                  hData_Hemi_2VtxAll_CutEvt_Mass->Fill(     tree_Hemi_Vtx_Mass->at(0) );
                  hData_Hemi_2VtxAll_CutEvt_Mass->Fill(     tree_Hemi_Vtx_Mass->at(1) );
                  hData_Hemi_2VtxAll_CutEvt_massMu->Fill(   tree_HemiMu_mass->at(0) );
                  hData_Hemi_2VtxAll_CutEvt_massMu->Fill(   tree_HemiMu_mass->at(1) );
                  hData_Hemi_2VtxAll_CutEvt_ptMu->Fill(     tree_HemiMu_pt->at(0) );
                  hData_Hemi_2VtxAll_CutEvt_ptMu->Fill(     tree_HemiMu_pt->at(1) );
                  hData_Hemi_2VtxAll_CutEvt_dRMu->Fill(     tree_HemiMu_dR->at(0) );
                  hData_Hemi_2VtxAll_CutEvt_dRMu->Fill(     tree_HemiMu_dR->at(1) );
                  }

                  if ( nVtx == 1 && isCutVtx ) {
                  hData_Hemi_1Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
                  hData_Hemi_1Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
                  hData_Hemi_1Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
                  hData_Hemi_1Vtx_CutVtx_nTrks->Fill(  nTrks );
                  hData_Hemi_1Vtx_CutVtx_mass->Fill(   mass );
                  hData_Hemi_1Vtx_CutVtx_Mass->Fill(   VtxMass );
                  hData_Hemi_1Vtx_CutVtx_massMu->Fill( massMu );
                  hData_Hemi_1Vtx_CutVtx_ptMu->Fill(   ptMu );
                  hData_Hemi_1Vtx_CutVtx_dRMu->Fill(   dRMu );
                  hData_Hemi_1Vtx_CutVtx_Mmumu->Fill(  tree_Mmumu );
                  hData_Hemi_1Vtx_CutVtx_PFMet->Fill(  tree_PFMet_et );
                  hData_Hemi_1Vtx_CutVtx_njetNOmuAll->Fill( tree_njetNOmu );
                  // hData_Hemi_1Vtx_CutVtx_lead_ptmin->Fill(  hemi_lead_ptmin );
                  hData_Hemi_1Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
                  }

                  if ( nVtx == 2 && isCutVtx ) {
                  hData_Hemi_2Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
                  hData_Hemi_2Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
                  hData_Hemi_2Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
                  hData_Hemi_2Vtx_CutVtx_nTrks->Fill(  nTrks );
                  hData_Hemi_2Vtx_CutVtx_mass->Fill(   mass );
                  hData_Hemi_2Vtx_CutVtx_Mass->Fill(   VtxMass );
                  hData_Hemi_2Vtx_CutVtx_massMu->Fill( massMu );
                  hData_Hemi_2Vtx_CutVtx_ptMu->Fill(   ptMu );
                  hData_Hemi_2Vtx_CutVtx_dRMu->Fill(   dRMu );
                  hData_Hemi_2Vtx_CutVtx_Mmumu->Fill(  tree_Mmumu );
                  hData_Hemi_2Vtx_CutVtx_PFMet->Fill(  tree_PFMet_et );
                  hData_Hemi_2Vtx_CutVtx_njetNOmuAll->Fill( tree_njetNOmu );
                  // hData_Hemi_2Vtx_CutVtx_lead_ptmin->Fill(  hemi_lead_ptmin );
                  hData_Hemi_2Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
                  hData_Hemi_2VtxAll_CutVtx_BDTvtx->Fill(   tree_Hemi_Vtx_MVAval_Step1->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_BDTvtx->Fill(   tree_Hemi_Vtx_MVAval_Step1->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_njetNOmu->Fill( tree_Hemi_njet_nomu->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_njetNOmu->Fill( tree_Hemi_njet_nomu->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_nTrks->Fill(    tree_Hemi_nTrks->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_nTrks->Fill(    tree_Hemi_nTrks->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_mass->Fill(     tree_Hemi_mass->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_mass->Fill(     tree_Hemi_mass->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_Mass->Fill(     tree_Hemi_Vtx_Mass->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_Mass->Fill(     tree_Hemi_Vtx_Mass->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_massMu->Fill(   tree_HemiMu_mass->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_massMu->Fill(   tree_HemiMu_mass->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_ptMu->Fill(     tree_HemiMu_pt->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_ptMu->Fill(     tree_HemiMu_pt->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_dRMu->Fill(     tree_HemiMu_dR->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_dRMu->Fill(     tree_HemiMu_dR->at(1) );
                  }

                  if ( nVtx == 1 && isCutVtx && isCutEvt ) {
                  hData_Hemi_1Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_mass->Fill(   mass );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_njetNOmuAll->Fill( tree_njetNOmu );
                  // hData_Hemi_1Vtx_CutVtx_CutEvt_lead_ptmin->Fill(  hemi_lead_ptmin );
                  hData_Hemi_1Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
                  }

                  if ( nVtx == 2 && isCutVtx && isCutEvt ) {
                  hData_Hemi_2Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_mass->Fill(   mass );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_njetNOmuAll->Fill( tree_njetNOmu );
                  // hData_Hemi_2Vtx_CutVtx_CutEvt_lead_ptmin->Fill(  hemi_lead_ptmin );
                  hData_Hemi_2Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_BDTvtx->Fill(   tree_Hemi_Vtx_MVAval_Step1->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_BDTvtx->Fill(   tree_Hemi_Vtx_MVAval_Step1->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_njetNOmu->Fill( tree_Hemi_njet_nomu->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_njetNOmu->Fill( tree_Hemi_njet_nomu->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_nTrks->Fill(    tree_Hemi_nTrks->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_nTrks->Fill(    tree_Hemi_nTrks->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_mass->Fill(     tree_Hemi_mass->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_mass->Fill(     tree_Hemi_mass->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_Mass->Fill(     tree_Hemi_Vtx_Mass->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_Mass->Fill(     tree_Hemi_Vtx_Mass->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_massMu->Fill(   tree_HemiMu_mass->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_massMu->Fill(   tree_HemiMu_mass->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_ptMu->Fill(     tree_HemiMu_pt->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_ptMu->Fill(     tree_HemiMu_pt->at(1) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_dRMu->Fill(     tree_HemiMu_dR->at(0) );
                  hData_Hemi_2VtxAll_CutVtx_CutEvt_dRMu->Fill(     tree_HemiMu_dR->at(1) );
                  }
               }
         } // End of Signal Region of Daniel'sMethod


         //$$$$
         // Forward Control Region

      isHemiVtx1 = false; isHemiVtx2 = false;
      isCutVtx = false; isCutVtx1 = false; isCutVtx2 = false;
      isCutEvt = false;
      nVtx = 0; 
      BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      VtxMass = -1.;
      
      if ( tree_Filter &&
         ((abs(tree_Hemi_eta->at(0)) > 2.4 && abs(tree_Hemi_eta->at(0)) < 3.0) ||
         (abs(tree_Hemi_eta->at(1)) > 2.4 && abs(tree_Hemi_eta->at(1)) < 3.0)) ) {

            hData_CRfwd_Mmumu->Fill( tree_Mmumu );
            hData_CRfwd_BDTevt->Fill( tree_Evts_MVAval );
   
         if ( tree_Hemi_Vtx_NChi2->at(0) > 0 && tree_Hemi_Vtx_NChi2->at(0) < 10
         && tree_Hemi_Vtx_step->at(0) >= 1 && tree_Hemi_Vtx_step->at(0) <= 2 
          ) {
            isHemiVtx1 = true;
            BDTvtx1 = tree_Hemi_Vtx_MVAval_Step1->at(0);
            BDTvtx  = BDTvtx1;
            VtxMass = tree_Hemi_Vtx_Mass->at(0);
         }

         if ( tree_Hemi_Vtx_NChi2->at(1) > 0 && tree_Hemi_Vtx_NChi2->at(1) < 10
         && tree_Hemi_Vtx_step->at(1) >= 1 && tree_Hemi_Vtx_step->at(1) <= 2 
         ) { 
            isHemiVtx2 = true;
            BDTvtx2 = tree_Hemi_Vtx_MVAval_Step1->at(1);
            if ( BDTvtx2 > BDTvtx ) BDTvtx = BDTvtx2;
            if ( tree_Hemi_Vtx_Mass->at(1)  > VtxMass )   VtxMass = tree_Hemi_Vtx_Mass->at(1);
         }

         if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;
         else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;
   
         //$$
         if ( VtxMass > 8. ) isCutVtx = true; 
         //$$

         hData_CRfwd_Hemi_ptmin->Fill( hemi_ptmin );
         hData_CRfwd_Hemi_njetNOmu->Fill(	njetNOmu );
         hData_CRfwd_Hemi_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_mass->Fill(   mass );
         hData_CRfwd_Hemi_massMu->Fill( massMu );
         hData_CRfwd_Hemi_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_dRMu->Fill(   dRMu );

         //$$
            if ( hemi_ptmin > 80. ) isCutEvt = true;   
         //$$

         if ( isCutEvt ) {
         hData_CRfwd_CutEvt_Mmumu->Fill( tree_Mmumu );
         hData_CRfwd_Hemi_CutEvt_njetNOmu->Fill(	njetNOmu );
         hData_CRfwd_Hemi_CutEvt_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_CutEvt_mass->Fill(   mass );
         hData_CRfwd_Hemi_CutEvt_massMu->Fill( massMu );
         hData_CRfwd_Hemi_CutEvt_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_CutEvt_dRMu->Fill(   dRMu );
         hData_CRfwd_Hemi_CutEvt_ptmin->Fill( hemi_ptmin );
         }

         if ( nVtx == 1 ) {
         hData_CRfwd_Hemi_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRfwd_Hemi_1Vtx_BDTvtx->Fill( BDTvtx );
         hData_CRfwd_Hemi_1Vtx_njetNOmu->Fill(   njetNOmu );
         hData_CRfwd_Hemi_1Vtx_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_1Vtx_mass->Fill(   mass );
         hData_CRfwd_Hemi_1Vtx_Mass->Fill(   VtxMass );
         hData_CRfwd_Hemi_1Vtx_massMu->Fill( massMu );
         hData_CRfwd_Hemi_1Vtx_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_1Vtx_dRMu->Fill(   dRMu );
         hData_CRfwd_Hemi_1Vtx_Mmumu->Fill(  tree_Mmumu );
         hData_CRfwd_Hemi_1Vtx_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 2 ) {
         hData_CRfwd_Hemi_2Vtx_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRfwd_Hemi_2Vtx_BDTvtx->Fill( BDTvtx );
         hData_CRfwd_Hemi_2Vtx_njetNOmu->Fill(   njetNOmu );
         hData_CRfwd_Hemi_2Vtx_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_2Vtx_mass->Fill(   mass );
         hData_CRfwd_Hemi_2Vtx_Mass->Fill(   VtxMass );
         hData_CRfwd_Hemi_2Vtx_massMu->Fill( massMu );
         hData_CRfwd_Hemi_2Vtx_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_2Vtx_dRMu->Fill(   dRMu );
         hData_CRfwd_Hemi_2Vtx_Mmumu->Fill(  tree_Mmumu );
         hData_CRfwd_Hemi_2Vtx_ptmin->Fill(  hemi_ptmin );


                        //------------------------
            int NewnVtx = -1;
            if (tree_Hemi_Merging->at(0) && tree_Hemi_Merging->at(1)) {NewnVtx = 2;}
            else if (tree_Hemi_Merging->at(0)==0 && tree_Hemi_Merging->at(1) == 0){NewnVtx = 0;}
            else {NewnVtx = 1;}
            float NewVtxMass = -10;
            bool SecTight = false;
 
            if (NewnVtx == 2)
               {
                  if (tree_Hemi_SecVtx_Mass->at(0)<tree_Hemi_SecVtx_Mass->at(1)) NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                  bool NewisCutVtx = false;
                  if (tree_Hemi_SecVtx_isTight->at(0) == true && tree_Hemi_SecVtx_isTight->at(1) == true) SecTight = true;
                  if (NewVtxMass > 8) NewisCutVtx = true;
                  if (  isCutEvt  && SecTight ) 
                     {
                        hData_CRfwd_Hemi_2SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                     }
                  if (NewisCutVtx && isCutEvt  && SecTight)
                     {
                        hData_CRfwd_Hemi_2SecVtx_CutVtx_CutEvt_Mass->Fill(   NewVtxMass );
                     }
               }
            if (NewnVtx == 1)
               {  
                  bool Merge = false ;
                  if (tree_Hemi_Merging->at(0)) Merge = true;
                  SecTight = tree_Hemi_SecVtx_isTight->at(1);
                  if (Merge) {NewVtxMass = tree_Hemi_SecVtx_Mass->at(0);SecTight = tree_Hemi_SecVtx_isTight->at(0);}

                  else NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                  bool NewisCutVtx = false;
                  if (NewVtxMass > 8) NewisCutVtx = true;
                  if (  isCutEvt  && SecTight ) 
                     {
                        hData_CRfwd_Hemi_1SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                        hData_CRfwd_Hemi_1SecVtx_CutEvt_Merge->Fill(Merge);
                     }
                  if (NewisCutVtx && isCutEvt  && SecTight)
                     {
                        hData_CRfwd_Hemi_1SecVtx_CutVtx_CutEvt_Mass->Fill(  NewVtxMass );
                        hData_CRfwd_Hemi_1SecVtx_CutVtx_CutEvt_Merge->Fill(Merge);
                     }
               }

            //--------------------------------------
         }

         if ( nVtx == 1 && isCutEvt ) {
         hData_CRfwd_Hemi_1Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRfwd_Hemi_1Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
         hData_CRfwd_Hemi_1Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
         hData_CRfwd_Hemi_1Vtx_CutEvt_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_1Vtx_CutEvt_mass->Fill(   mass );
         hData_CRfwd_Hemi_1Vtx_CutEvt_Mass->Fill(   VtxMass );
         hData_CRfwd_Hemi_1Vtx_CutEvt_massMu->Fill( massMu );
         hData_CRfwd_Hemi_1Vtx_CutEvt_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_1Vtx_CutEvt_dRMu->Fill(   dRMu );
         hData_CRfwd_Hemi_1Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
         hData_CRfwd_Hemi_1Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 2 && isCutEvt ) {
         hData_CRfwd_Hemi_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRfwd_Hemi_2Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
         hData_CRfwd_Hemi_2Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
         hData_CRfwd_Hemi_2Vtx_CutEvt_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_2Vtx_CutEvt_mass->Fill(   mass );
         hData_CRfwd_Hemi_2Vtx_CutEvt_Mass->Fill(   VtxMass );
         hData_CRfwd_Hemi_2Vtx_CutEvt_massMu->Fill( massMu );
         hData_CRfwd_Hemi_2Vtx_CutEvt_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_2Vtx_CutEvt_dRMu->Fill(   dRMu );
         hData_CRfwd_Hemi_2Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
         hData_CRfwd_Hemi_2Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 1 && isCutVtx ) {
         hData_CRfwd_Hemi_1Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRfwd_Hemi_1Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
         hData_CRfwd_Hemi_1Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_1Vtx_CutVtx_mass->Fill(   mass );
         hData_CRfwd_Hemi_1Vtx_CutVtx_Mass->Fill(   VtxMass );
         hData_CRfwd_Hemi_1Vtx_CutVtx_massMu->Fill( massMu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_dRMu->Fill(   dRMu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_Mmumu->Fill(  tree_Mmumu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 2 && isCutVtx ) {
         hData_CRfwd_Hemi_2Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRfwd_Hemi_2Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
         hData_CRfwd_Hemi_2Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_2Vtx_CutVtx_mass->Fill(   mass );
         hData_CRfwd_Hemi_2Vtx_CutVtx_Mass->Fill(   VtxMass );
         hData_CRfwd_Hemi_2Vtx_CutVtx_massMu->Fill( massMu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_dRMu->Fill(   dRMu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_Mmumu->Fill(  tree_Mmumu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 1 && isCutVtx && isCutEvt ) {
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_mass->Fill(   mass );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
         hData_CRfwd_Hemi_1Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 2 && isCutVtx && isCutEvt ) {
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_mass->Fill(   mass );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
         hData_CRfwd_Hemi_2Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
         }
      } // endif Forward Control Region

      //$$$$
      // Low Pt Control Region

      isHemiVtx1 = false; isHemiVtx2 = false;
      isCutVtx = false; isCutVtx1 = false; isCutVtx2 = false;
      isCutEvt = false;
      nVtx = 0; 
      BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      VtxMass = -1.;
   
      //$$
      if ( tree_Filter &&
            abs(tree_Hemi_eta->at(0)) < 2.4 && abs(tree_Hemi_eta->at(1)) < 2.4
            && hemi_ptmin > 40. && hemi_ptmin < 80. ) {
      //$$

            hData_CRlowpt_Mmumu->Fill( tree_Mmumu );
            hData_CRlowpt_BDTevt->Fill( tree_Evts_MVAval );
   
            if ( tree_Hemi_Vtx_NChi2->at(0) > 0 && tree_Hemi_Vtx_NChi2->at(0) < 10
                  && tree_Hemi_Vtx_step->at(0) >= 1 && tree_Hemi_Vtx_step->at(0) <= 2 
               ) {
                  isHemiVtx1 = true;
                  BDTvtx1 = tree_Hemi_Vtx_MVAval_Step1->at(0);
                  BDTvtx  = BDTvtx1;
                  VtxMass = tree_Hemi_Vtx_Mass->at(0);
               }

            if ( tree_Hemi_Vtx_NChi2->at(1) > 0 && tree_Hemi_Vtx_NChi2->at(1) < 10
               && tree_Hemi_Vtx_step->at(1) >= 1 && tree_Hemi_Vtx_step->at(1) <= 2 
               ) { 
                  isHemiVtx2 = true;
                  BDTvtx2 = tree_Hemi_Vtx_MVAval_Step1->at(1);
                  if ( BDTvtx2 > BDTvtx ) BDTvtx = BDTvtx2;
                  if ( tree_Hemi_Vtx_Mass->at(1)  > VtxMass )   VtxMass = tree_Hemi_Vtx_Mass->at(1);
               }

            if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;
            else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;
   
         //$$
            if ( VtxMass > 8. ) isCutVtx = true; 
         //$$

            hData_CRlowpt_Hemi_njetNOmu->Fill(	njetNOmu );
            hData_CRlowpt_Hemi_nTrks->Fill(  nTrks );
            hData_CRlowpt_Hemi_mass->Fill(   mass );
            hData_CRlowpt_Hemi_massMu->Fill( massMu );
            hData_CRlowpt_Hemi_ptMu->Fill(   ptMu );
            hData_CRlowpt_Hemi_dRMu->Fill(   dRMu );
            hData_CRlowpt_Hemi_ptmin->Fill(  hemi_ptmin );

            if ( nVtx == 1 ) {
            hData_CRlowpt_Hemi_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
            hData_CRlowpt_Hemi_1Vtx_BDTvtx->Fill( BDTvtx );
            hData_CRlowpt_Hemi_1Vtx_njetNOmu->Fill(   njetNOmu );
            hData_CRlowpt_Hemi_1Vtx_nTrks->Fill(  nTrks );
            hData_CRlowpt_Hemi_1Vtx_mass->Fill(   mass );
            hData_CRlowpt_Hemi_1Vtx_Mass->Fill(   VtxMass );
            hData_CRlowpt_Hemi_1Vtx_massMu->Fill( massMu );
            hData_CRlowpt_Hemi_1Vtx_ptMu->Fill(   ptMu );
            hData_CRlowpt_Hemi_1Vtx_dRMu->Fill(   dRMu );
            hData_CRlowpt_Hemi_1Vtx_Mmumu->Fill(  tree_Mmumu );
            hData_CRlowpt_Hemi_1Vtx_ptmin->Fill(  hemi_ptmin );
            }

            if ( nVtx == 2 ) {
            hData_CRlowpt_Hemi_2Vtx_BDTevt->Fill( tree_Evts_MVAval );
            hData_CRlowpt_Hemi_2Vtx_BDTvtx->Fill( BDTvtx );
            hData_CRlowpt_Hemi_2Vtx_njetNOmu->Fill(   njetNOmu );
            hData_CRlowpt_Hemi_2Vtx_nTrks->Fill(  nTrks );
            hData_CRlowpt_Hemi_2Vtx_mass->Fill(   mass );
            hData_CRlowpt_Hemi_2Vtx_Mass->Fill(   VtxMass );
            hData_CRlowpt_Hemi_2Vtx_massMu->Fill( massMu );
            hData_CRlowpt_Hemi_2Vtx_ptMu->Fill(   ptMu );
            hData_CRlowpt_Hemi_2Vtx_dRMu->Fill(   dRMu );
            hData_CRlowpt_Hemi_2Vtx_Mmumu->Fill(  tree_Mmumu );
            hData_CRlowpt_Hemi_2Vtx_ptmin->Fill(  hemi_ptmin );

                        //------------------------
            int NewnVtx = -1;
            if (tree_Hemi_Merging->at(0) && tree_Hemi_Merging->at(1)) {NewnVtx = 2;}
            else if (tree_Hemi_Merging->at(0)==0 && tree_Hemi_Merging->at(1) == 0){NewnVtx = 0;}
            else {NewnVtx = 1;}
            float NewVtxMass = -10;
            bool SecTight = false;
            
            if (NewnVtx == 2)
               {
                  if (tree_Hemi_SecVtx_Mass->at(0)<tree_Hemi_SecVtx_Mass->at(1)) NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                  if (tree_Hemi_SecVtx_isTight->at(0) == true && tree_Hemi_SecVtx_isTight->at(1) == true) SecTight = true;
                  bool NewisCutVtx = false;
                  if (NewVtxMass > 8) NewisCutVtx = true;
                  if (  isCutEvt  &&  SecTight ) 
                     {
                        hData_CRlowpt_Hemi_2SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                     }
                  if (NewisCutVtx && isCutEvt &&  SecTight)
                     {
                        hData_CRlowpt_Hemi_2SecVtx_CutVtx_CutEvt_Mass->Fill(   NewVtxMass );
                     }
               }
            if (NewnVtx == 1)
               {  
                  bool Merge = false ;
                  if (tree_Hemi_Merging->at(0)) Merge = true;
                  SecTight = tree_Hemi_SecVtx_isTight->at(1);
                  if (Merge) {NewVtxMass = tree_Hemi_SecVtx_Mass->at(0);SecTight = tree_Hemi_SecVtx_isTight->at(0);}
                  else NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                  bool NewisCutVtx = false;
                  if (NewVtxMass > 8) NewisCutVtx = true;
                  if (  isCutEvt && SecTight) 
                     {
                        hData_CRlowpt_Hemi_1SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                        hData_CRlowpt_Hemi_1SecVtx_CutEvt_Merge->Fill(Merge);
                     }
                  if (NewisCutVtx && isCutEvt && SecTight)
                     {
                        hData_CRlowpt_Hemi_1SecVtx_CutVtx_CutEvt_Mass->Fill(  NewVtxMass );
                        hData_CRlowpt_Hemi_1SecVtx_CutVtx_CutEvt_Merge->Fill(Merge);
                     }
               }

            //--------------------------------------
            }

            if ( nVtx == 1 && isCutVtx ) {
            hData_CRlowpt_Hemi_1Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_nTrks->Fill(  nTrks );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_mass->Fill(   mass );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_Mass->Fill(   VtxMass );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_massMu->Fill( massMu );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_ptMu->Fill(   ptMu );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_dRMu->Fill(   dRMu );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_Mmumu->Fill(  tree_Mmumu );
            hData_CRlowpt_Hemi_1Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
            }

            if ( nVtx == 2 && isCutVtx ) {
            hData_CRlowpt_Hemi_2Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_nTrks->Fill(  nTrks );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_mass->Fill(   mass );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_Mass->Fill(   VtxMass );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_massMu->Fill( massMu );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_ptMu->Fill(   ptMu );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_dRMu->Fill(   dRMu );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_Mmumu->Fill(  tree_Mmumu );
            hData_CRlowpt_Hemi_2Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
            }

      } // endif Low Pt Control Region
      //$$$$

   //$$$$
   // Loose Vertex Control Region

      isHemiVtx1 = false; isHemiVtx2 = false;
      isCutVtx = false; isCutVtx1 = false; isCutVtx2 = false;
      isCutEvt = false;
      nVtx = 0; 
      BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      VtxMass = -1.;
   
      if ( tree_Filter &&
            abs(tree_Hemi_eta->at(0)) < 2.4 && abs(tree_Hemi_eta->at(1)) < 2.4 ) {

      //$$
      //    bool isHemiVtx1Tight = false, isHemiVtx2Tight = false;
         if ( tree_Hemi_Vtx_NChi2->at(0) > 0 && tree_Hemi_Vtx_NChi2->at(0) < 10
                //  && tree_Hemi_Vtx_step->at(0) >= 1 && tree_Hemi_Vtx_step->at(0) <= 4 
               && tree_Hemi_Vtx_step->at(0) >= 3 && tree_Hemi_Vtx_step->at(0) <= 4 
            ) {
                  isHemiVtx1 = true;
                  BDTvtx1 = tree_Hemi_Vtx_MVAval->at(0);
                  BDTvtx  = BDTvtx1;
                  VtxMass = tree_Hemi_Vtx_Mass->at(0);
                  //      if ( tree_Hemi_Vtx_step->at(0) >= 1 && tree_Hemi_Vtx_step->at(0) <= 2 ) isHemiVtx1Tight = true;
               }

         if ( tree_Hemi_Vtx_NChi2->at(1) > 0 && tree_Hemi_Vtx_NChi2->at(1) < 10
      //         && tree_Hemi_Vtx_step->at(1) >= 1 && tree_Hemi_Vtx_step->at(1) <= 4 
            && tree_Hemi_Vtx_step->at(1) >= 3 && tree_Hemi_Vtx_step->at(1) <= 4 
            ) { 
                  isHemiVtx2 = true;
                  BDTvtx2 = tree_Hemi_Vtx_MVAval->at(1);
                  if ( BDTvtx2 > BDTvtx ) BDTvtx = BDTvtx2;
                  if ( tree_Hemi_Vtx_Mass->at(1)  > VtxMass )   VtxMass = tree_Hemi_Vtx_Mass->at(1);
                  //      if ( tree_Hemi_Vtx_step->at(1) >= 1 && tree_Hemi_Vtx_step->at(1) <= 2 ) isHemiVtx2Tight = true;
              }

         //    if      (  isHemiVtx1 && isHemiVtx2 && !isHemiVtx1Tight && !isHemiVtx2Tight )    nVtx = 2;
         //    else if ( (isHemiVtx1 && !isHemiVtx1Tight) || (isHemiVtx2 && !isHemiVtx2Tight) ) nVtx = 1;
            if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;
            else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;
   
            if ( VtxMass > 8. ) isCutVtx = true; 
         //$$

         //$$
            if ( hemi_ptmin > 80. ) isCutEvt = true;   
         //$$

         if ( nVtx == 0 ) {
         hData_CRloose_Hemi_0Vtx_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRloose_Hemi_0Vtx_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_0Vtx_mass->Fill(   mass );
         hData_CRloose_Hemi_0Vtx_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_0Vtx_njetNOmuAll->Fill(   tree_njetNOmu );
         //   hData_CRloose_Hemi_0Vtx_lead_ptmin->Fill( hemi_lead_ptmin );
         hData_CRloose_Hemi_0Vtx_ptmin->Fill( hemi_ptmin );
         }
         if ( nVtx == 0 && isCutEvt ) {
         hData_CRloose_Hemi_0Vtx_CutEvt_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_0Vtx_CutEvt_mass->Fill(   mass );
         hData_CRloose_Hemi_0Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_0Vtx_CutEvt_njetNOmuAll->Fill(   tree_njetNOmu );
         //   hData_CRloose_Hemi_0Vtx_CutEvt_lead_ptmin->Fill( hemi_lead_ptmin );
         hData_CRloose_Hemi_0Vtx_CutEvt_ptmin->Fill( hemi_ptmin );
         }

         if ( nVtx == 1 ) {
         hData_CRloose_Hemi_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRloose_Hemi_1Vtx_BDTvtx->Fill( BDTvtx );
         hData_CRloose_Hemi_1Vtx_njetNOmu->Fill(   njetNOmu );
         hData_CRloose_Hemi_1Vtx_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_1Vtx_mass->Fill(   mass );
         hData_CRloose_Hemi_1Vtx_Mass->Fill(   VtxMass );
         hData_CRloose_Hemi_1Vtx_massMu->Fill( massMu );
         hData_CRloose_Hemi_1Vtx_ptMu->Fill(   ptMu );
         hData_CRloose_Hemi_1Vtx_dRMu->Fill(   dRMu );
         hData_CRloose_Hemi_1Vtx_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_1Vtx_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 2 ) {
         hData_CRloose_Hemi_2Vtx_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRloose_Hemi_2Vtx_BDTvtx->Fill( BDTvtx );
         hData_CRloose_Hemi_2Vtx_njetNOmu->Fill(   njetNOmu );
         hData_CRloose_Hemi_2Vtx_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_2Vtx_mass->Fill(   mass );
         hData_CRloose_Hemi_2Vtx_Mass->Fill(   VtxMass );
         hData_CRloose_Hemi_2Vtx_massMu->Fill( massMu );
         hData_CRloose_Hemi_2Vtx_ptMu->Fill(   ptMu );
         hData_CRloose_Hemi_2Vtx_dRMu->Fill(   dRMu );
         hData_CRloose_Hemi_2Vtx_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_2Vtx_ptmin->Fill(  hemi_ptmin );
                     //------------------------
            int NewnVtx = -1;
            if (tree_Hemi_Merging->at(0) && tree_Hemi_Merging->at(1)) {NewnVtx = 2;}
            else if (tree_Hemi_Merging->at(0)==0 && tree_Hemi_Merging->at(1) == 0){NewnVtx = 0;}
            else {NewnVtx = 1;}
            float NewVtxMass = -10;
            bool SecTight = false;
            
            if (NewnVtx == 2)
               {
                  if (tree_Hemi_SecVtx_Mass->at(0)<tree_Hemi_SecVtx_Mass->at(1)) NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                  bool NewisCutVtx = false;
                  if (tree_Hemi_SecVtx_isTight->at(0) ==false && tree_Hemi_SecVtx_isTight->at(1) == false) SecTight = false;
                  if (NewVtxMass > 8) NewisCutVtx = true;
                  if (  isCutEvt && !SecTight ) 
                     {
                        hData_CRloose_Hemi_2SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                     }
                  if (NewisCutVtx && isCutEvt && !SecTight)
                     {
                        hData_CRloose_Hemi_2SecVtx_CutVtx_CutEvt_Mass->Fill(   NewVtxMass );
                     }
               }
            if (NewnVtx == 1)
               {  
                  bool Merge = false ;
                  if (tree_Hemi_Merging->at(0)) Merge = true;
                  SecTight = tree_Hemi_SecVtx_isTight->at(1);
                  if (Merge) {NewVtxMass = tree_Hemi_SecVtx_Mass->at(0);SecTight = tree_Hemi_SecVtx_isTight->at(0);}
                  else NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                  bool NewisCutVtx = false;
                  if (NewVtxMass > 8) NewisCutVtx = true;
                  if (  isCutEvt && !SecTight) 
                     {
                        hData_CRloose_Hemi_1SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                        hData_CRloose_Hemi_1SecVtx_CutEvt_Merge->Fill(Merge);
                     }
                  if (NewisCutVtx && isCutEvt && !SecTight)
                     {
                        hData_CRloose_Hemi_1SecVtx_CutVtx_CutEvt_Mass->Fill(  NewVtxMass );
                        hData_CRloose_Hemi_1SecVtx_CutVtx_CutEvt_Merge->Fill(Merge);
                     }
               }

            //--------------------------------------
         }

         if ( nVtx == 1 && isCutEvt ) {
         hData_CRloose_Hemi_1Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRloose_Hemi_1Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
         hData_CRloose_Hemi_1Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
         hData_CRloose_Hemi_1Vtx_CutEvt_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_1Vtx_CutEvt_mass->Fill(   mass );
         hData_CRloose_Hemi_1Vtx_CutEvt_Mass->Fill(   VtxMass );
         hData_CRloose_Hemi_1Vtx_CutEvt_massMu->Fill( massMu );
         hData_CRloose_Hemi_1Vtx_CutEvt_ptMu->Fill(   ptMu );
         hData_CRloose_Hemi_1Vtx_CutEvt_dRMu->Fill(   dRMu );
         hData_CRloose_Hemi_1Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_1Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 2 && isCutEvt ) {
         hData_CRloose_Hemi_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRloose_Hemi_2Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
         hData_CRloose_Hemi_2Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
         hData_CRloose_Hemi_2Vtx_CutEvt_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_2Vtx_CutEvt_mass->Fill(   mass );
         hData_CRloose_Hemi_2Vtx_CutEvt_Mass->Fill(   VtxMass );
         hData_CRloose_Hemi_2Vtx_CutEvt_massMu->Fill( massMu );
         hData_CRloose_Hemi_2Vtx_CutEvt_ptMu->Fill(   ptMu );
         hData_CRloose_Hemi_2Vtx_CutEvt_dRMu->Fill(   dRMu );
         hData_CRloose_Hemi_2Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_2Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 1 && isCutVtx ) {
         hData_CRloose_Hemi_1Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRloose_Hemi_1Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
         hData_CRloose_Hemi_1Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
         hData_CRloose_Hemi_1Vtx_CutVtx_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_1Vtx_CutVtx_mass->Fill(   mass );
         hData_CRloose_Hemi_1Vtx_CutVtx_Mass->Fill(   VtxMass );
         hData_CRloose_Hemi_1Vtx_CutVtx_massMu->Fill( massMu );
         hData_CRloose_Hemi_1Vtx_CutVtx_ptMu->Fill(   ptMu );
         hData_CRloose_Hemi_1Vtx_CutVtx_dRMu->Fill(   dRMu );
         hData_CRloose_Hemi_1Vtx_CutVtx_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_1Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 2 && isCutVtx ) {
         hData_CRloose_Hemi_2Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRloose_Hemi_2Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
         hData_CRloose_Hemi_2Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
         hData_CRloose_Hemi_2Vtx_CutVtx_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_2Vtx_CutVtx_mass->Fill(   mass );
         hData_CRloose_Hemi_2Vtx_CutVtx_Mass->Fill(   VtxMass );
         hData_CRloose_Hemi_2Vtx_CutVtx_massMu->Fill( massMu );
         hData_CRloose_Hemi_2Vtx_CutVtx_ptMu->Fill(   ptMu );
         hData_CRloose_Hemi_2Vtx_CutVtx_dRMu->Fill(   dRMu );
         hData_CRloose_Hemi_2Vtx_CutVtx_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_2Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 1 && isCutVtx && isCutEvt ) {
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_mass->Fill(   mass );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_1Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
         }

         if ( nVtx == 2 && isCutVtx && isCutEvt ) {
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_mass->Fill(   mass );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
         hData_CRloose_Hemi_2Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
         }

         } // endif Loose Vertex Control Region
      //$$$$

//$$$$
// Loose Vertex Low Pt Control Region

   isHemiVtx1 = false; isHemiVtx2 = false;
   isCutVtx = false; isCutVtx1 = false; isCutVtx2 = false;
   isCutEvt = false;
   nVtx = 0; 
   BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
   VtxMass = -1.;
   
 if ( tree_Filter &&
//$$
      abs(tree_Hemi_eta->at(0)) < 2.4 && abs(tree_Hemi_eta->at(1)) < 2.4 ) {
//$$

   if ( tree_Hemi_Vtx_NChi2->at(0) > 0 && tree_Hemi_Vtx_NChi2->at(0) < 10
//$$
        && tree_Hemi_Vtx_step->at(0) >= 3 && tree_Hemi_Vtx_step->at(0) <= 4 
//$$
      ) {
     isHemiVtx1 = true;
     BDTvtx1 = tree_Hemi_Vtx_MVAval->at(0);
     BDTvtx  = BDTvtx1;
     VtxMass = tree_Hemi_Vtx_Mass->at(0);
   }

   if ( tree_Hemi_Vtx_NChi2->at(1) > 0 && tree_Hemi_Vtx_NChi2->at(1) < 10
//$$
        && tree_Hemi_Vtx_step->at(1) >= 3 && tree_Hemi_Vtx_step->at(1) <= 4 
//$$
      ) { 
     isHemiVtx2 = true;
     BDTvtx2 = tree_Hemi_Vtx_MVAval->at(1);
     if ( BDTvtx2 > BDTvtx ) BDTvtx = BDTvtx2;
     if ( tree_Hemi_Vtx_Mass->at(1)  > VtxMass )   VtxMass = tree_Hemi_Vtx_Mass->at(1);
   }

   if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;
   else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;
   
//$$
   if ( VtxMass > 8. ) isCutVtx = true; 
//$$

//$$
   if ( hemi_ptmin > 40. && hemi_ptmin < 80. ) isCutEvt = true;   
//$$

   float dRvtx = -1.;
   if ( nVtx == 2 ) {
     float eta_vtx1 = tree_Hemi_Vtx_eta->at(0);
     float eta_vtx2 = tree_Hemi_Vtx_eta->at(1);
     float phi_vtx1 = TMath::ATan2( tree_Hemi_Vtx_y->at(0), tree_Hemi_Vtx_x->at(0) );
     float phi_vtx2 = TMath::ATan2( tree_Hemi_Vtx_y->at(1), tree_Hemi_Vtx_x->at(1) );
     dRvtx = Deltar( eta_vtx1, phi_vtx1, eta_vtx2, phi_vtx2 );
     if ( dRvtx < 1.5 ) {
       nVtx = 1;
       if ( tree_Hemi_Vtx_Mass->at(0) > tree_Hemi_Vtx_Mass->at(1) ) {
         VtxMass = tree_Hemi_Vtx_Mass->at(0);
         BDTvtx = tree_Hemi_Vtx_MVAval->at(0);
       }
       else {
         VtxMass = tree_Hemi_Vtx_Mass->at(1);
         BDTvtx = tree_Hemi_Vtx_MVAval->at(1);
       } 
     }
     if ( isCutEvt ) hData_CRlooselowpt_Hemi_2Vtx_CutEvt_dRvtx->Fill( dRvtx );
     if ( isCutVtx && isCutEvt ) hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_dRvtx->Fill( dRvtx );

      
      
      //------------------------
      int NewnVtx = -1;
      if (tree_Hemi_Merging->at(0) && tree_Hemi_Merging->at(1)) {NewnVtx = 2;}
      else if (tree_Hemi_Merging->at(0)==0 && tree_Hemi_Merging->at(1) == 0){NewnVtx = 0;}
      else {NewnVtx = 1;}
      float NewVtxMass = -10;
      bool SecTight = false;
      if (NewnVtx == 2)
         {
            if (tree_Hemi_SecVtx_Mass->at(0)<tree_Hemi_SecVtx_Mass->at(1)) NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
            bool NewisCutVtx = false;
            if (tree_Hemi_SecVtx_isTight->at(0) ==false && tree_Hemi_SecVtx_isTight->at(1) == false) SecTight = false;
            if (NewVtxMass > 8) NewisCutVtx = true;
            if (  isCutEvt && !SecTight) 
               {
                  hData_CRlooselowpt_Hemi_2SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
               }
            if (NewisCutVtx && isCutEvt && !SecTight)
               {
                  hData_CRlooselowpt_Hemi_2SecVtx_CutVtx_CutEvt_Mass->Fill(   NewVtxMass );
               }
         }
      if (NewnVtx == 1)
         {  
            bool Merge = false ;
            if (tree_Hemi_Merging->at(0)) Merge = true;
            SecTight = tree_Hemi_SecVtx_isTight->at(1);
            if (Merge) {NewVtxMass = tree_Hemi_SecVtx_Mass->at(0);SecTight = tree_Hemi_SecVtx_isTight->at(0);}
            else NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
            bool NewisCutVtx = false;
            if (NewVtxMass > 8) NewisCutVtx = true;
            if (  isCutEvt && !SecTight ) 
               {
                  hData_CRlooselowpt_Hemi_1SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                  hData_CRlooselowpt_Hemi_1SecVtx_CutEvt_Merge->Fill(Merge);
               }
            if (NewisCutVtx && isCutEvt && !SecTight)
               {
                  hData_CRlooselowpt_Hemi_1SecVtx_CutVtx_CutEvt_Mass->Fill(  NewVtxMass );
                  hData_CRlooselowpt_Hemi_1SecVtx_CutVtx_CutEvt_Merge->Fill(Merge);
               }
         }

      //--------------------------------------

   }

   if ( nVtx == 1 && isCutEvt ) {
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_nTrks->Fill(  nTrks );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_mass->Fill(   mass );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_Mass->Fill(   VtxMass );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_massMu->Fill( massMu );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_ptMu->Fill(   ptMu );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_dRMu->Fill(   dRMu );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
     hData_CRlooselowpt_Hemi_1Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
   }

   if ( nVtx == 2 && isCutEvt ) {
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_nTrks->Fill(  nTrks );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_mass->Fill(   mass );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_Mass->Fill(   VtxMass );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_massMu->Fill( massMu );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_ptMu->Fill(   ptMu );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_dRMu->Fill(   dRMu );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
     hData_CRlooselowpt_Hemi_2Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
   }

   if ( nVtx == 1 && isCutVtx && isCutEvt ) {
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_mass->Fill(   mass );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
     hData_CRlooselowpt_Hemi_1Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
   }

   if ( nVtx == 2 && isCutVtx && isCutEvt ) {
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_mass->Fill(   mass );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
     hData_CRlooselowpt_Hemi_2Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
   }

 } // endif Loose Vertex Low Pt Control Region
//$$$$


      //$$$$
      // Same Sign Control Region

      isHemiVtx1 = false; isHemiVtx2 = false;
      isCutVtx = false; isCutVtx1 = false; isCutVtx2 = false;
      isCutEvt = false;
      nVtx = 0; 
      BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      VtxMass = -1.;
   
      if ( tree_FilterSameSign && !tree_Filter &&
            abs(tree_Hemi_eta->at(0)) < 2.4 && abs(tree_Hemi_eta->at(1)) < 2.4 ) {

            hData_SameSign_Mmumu->Fill( tree_MmumuSameSign );
            hData_SameSign_BDTevt->Fill( tree_Evts_MVAval );
   
            if ( tree_Hemi_Vtx_NChi2->at(0) > 0 && tree_Hemi_Vtx_NChi2->at(0) < 10
                  && tree_Hemi_Vtx_step->at(0) >= 1 && tree_Hemi_Vtx_step->at(0) <= 2 
               ) {
                     isHemiVtx1 = true;
                     BDTvtx1 = tree_Hemi_Vtx_MVAval_Step1->at(0);
                     BDTvtx  = BDTvtx1;
                     VtxMass = tree_Hemi_Vtx_Mass->at(0);
                 }

            if ( tree_Hemi_Vtx_NChi2->at(1) > 0 && tree_Hemi_Vtx_NChi2->at(1) < 10
               && tree_Hemi_Vtx_step->at(1) >= 1 && tree_Hemi_Vtx_step->at(1) <= 2 
               ) { 
                     isHemiVtx2 = true;
                     BDTvtx2 = tree_Hemi_Vtx_MVAval_Step1->at(1);
                     if ( BDTvtx2 > BDTvtx ) BDTvtx = BDTvtx2;
                     if ( tree_Hemi_Vtx_Mass->at(1)  > VtxMass ) VtxMass = tree_Hemi_Vtx_Mass->at(1);
                 }

            if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;
            else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;
            
            //$$
            if ( VtxMass > 8. ) isCutVtx = true; 
            //$$

            hData_SameSign_Hemi_ptmin->Fill( hemi_ptmin );
            hData_SameSign_Hemi_njetNOmu->Fill(	njetNOmu );
            hData_SameSign_Hemi_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_mass->Fill(   mass );
            hData_SameSign_Hemi_massMu->Fill( massMu );
            hData_SameSign_Hemi_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_dRMu->Fill(   dRMu );

            //$$
               if ( hemi_ptmin > 80. ) isCutEvt = true;   
            //$$

            if ( isCutEvt ) {
            hData_SameSign_CutEvt_Mmumu->Fill( tree_MmumuSameSign );
            hData_SameSign_Hemi_CutEvt_njetNOmu->Fill(	njetNOmu );
            hData_SameSign_Hemi_CutEvt_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_CutEvt_mass->Fill(   mass );
            hData_SameSign_Hemi_CutEvt_massMu->Fill( massMu );
            hData_SameSign_Hemi_CutEvt_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_CutEvt_dRMu->Fill(   dRMu );
            hData_SameSign_Hemi_CutEvt_ptmin->Fill( hemi_ptmin );
            }

            if ( nVtx == 1 ) {
            hData_SameSign_Hemi_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
            hData_SameSign_Hemi_1Vtx_BDTvtx->Fill( BDTvtx );
            hData_SameSign_Hemi_1Vtx_njetNOmu->Fill(   njetNOmu );
            hData_SameSign_Hemi_1Vtx_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_1Vtx_mass->Fill(   mass );
            hData_SameSign_Hemi_1Vtx_Mass->Fill(   VtxMass );
            hData_SameSign_Hemi_1Vtx_massMu->Fill( massMu );
            hData_SameSign_Hemi_1Vtx_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_1Vtx_dRMu->Fill(   dRMu );
            hData_SameSign_Hemi_1Vtx_Mmumu->Fill(  tree_MmumuSameSign );
            hData_SameSign_Hemi_1Vtx_ptmin->Fill(  hemi_ptmin );
            }

            if ( nVtx == 2 ) {
            hData_SameSign_Hemi_2Vtx_BDTevt->Fill( tree_Evts_MVAval );
            hData_SameSign_Hemi_2Vtx_BDTvtx->Fill( BDTvtx );
            hData_SameSign_Hemi_2Vtx_njetNOmu->Fill(   njetNOmu );
            hData_SameSign_Hemi_2Vtx_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_2Vtx_mass->Fill(   mass );
            hData_SameSign_Hemi_2Vtx_Mass->Fill(   VtxMass );
            hData_SameSign_Hemi_2Vtx_massMu->Fill( massMu );
            hData_SameSign_Hemi_2Vtx_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_2Vtx_dRMu->Fill(   dRMu );
            hData_SameSign_Hemi_2Vtx_Mmumu->Fill(  tree_MmumuSameSign );
            hData_SameSign_Hemi_2Vtx_ptmin->Fill(  hemi_ptmin );
                           //------------------------
               int NewnVtx = -1;
               if (tree_Hemi_Merging->at(0) && tree_Hemi_Merging->at(1)) {NewnVtx = 2;}
               else if (tree_Hemi_Merging->at(0)==0 && tree_Hemi_Merging->at(1) == 0){NewnVtx = 0;}
               else {NewnVtx = 1;}
               bool SecTight = false;
               
               float NewVtxMass = -10;

               if (NewnVtx == 2)
                  {
                     if (tree_Hemi_SecVtx_Mass->at(0)<tree_Hemi_SecVtx_Mass->at(1)) NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                     bool NewisCutVtx = false;
                     if (tree_Hemi_SecVtx_isTight->at(0) ==true && tree_Hemi_SecVtx_isTight->at(1) == true) SecTight = true;

                     if (NewVtxMass > 8) NewisCutVtx = true;
                     if (  isCutEvt && SecTight) 
                        {
                           hData_SameSign_Hemi_2SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                        }
                     if (NewisCutVtx && isCutEvt && SecTight)
                        {
                           hData_SameSign_Hemi_2SecVtx_CutVtx_CutEvt_Mass->Fill(   NewVtxMass );
                        }
                  }
               if (NewnVtx == 1)
                  {  
                     bool Merge = false ;
                     if (tree_Hemi_Merging->at(0)) Merge = true;
                     SecTight = tree_Hemi_SecVtx_isTight->at(1);
                     if (Merge) {NewVtxMass = tree_Hemi_SecVtx_Mass->at(0);SecTight = tree_Hemi_SecVtx_isTight->at(0);}
                     else NewVtxMass = tree_Hemi_SecVtx_Mass->at(1);
                     bool NewisCutVtx = false;
                     if (NewVtxMass > 8) NewisCutVtx = true;
                     if (  isCutEvt && SecTight) 
                        {
                           hData_SameSign_Hemi_1SecVtx_CutEvt_Mass->Fill(   NewVtxMass );
                           hData_SameSign_Hemi_1SecVtx_CutEvt_Merge->Fill(Merge);
                        }
                     if (NewisCutVtx && isCutEvt && SecTight)
                        {
                           hData_SameSign_Hemi_1SecVtx_CutVtx_CutEvt_Mass->Fill(  NewVtxMass );
                           hData_SameSign_Hemi_1SecVtx_CutVtx_CutEvt_Merge->Fill(Merge);
                        }
                  }

               //--------------------------------------
            }

            if ( nVtx == 1 && isCutEvt ) {
            hData_SameSign_Hemi_1Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
            hData_SameSign_Hemi_1Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
            hData_SameSign_Hemi_1Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
            hData_SameSign_Hemi_1Vtx_CutEvt_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_1Vtx_CutEvt_mass->Fill(   mass );
            hData_SameSign_Hemi_1Vtx_CutEvt_Mass->Fill(   VtxMass );
            hData_SameSign_Hemi_1Vtx_CutEvt_massMu->Fill( massMu );
            hData_SameSign_Hemi_1Vtx_CutEvt_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_1Vtx_CutEvt_dRMu->Fill(   dRMu );
            hData_SameSign_Hemi_1Vtx_CutEvt_Mmumu->Fill(  tree_MmumuSameSign );
            hData_SameSign_Hemi_1Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
            }

            if ( nVtx == 2 && isCutEvt ) {
            hData_SameSign_Hemi_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
            hData_SameSign_Hemi_2Vtx_CutEvt_BDTvtx->Fill( BDTvtx );
            hData_SameSign_Hemi_2Vtx_CutEvt_njetNOmu->Fill(   njetNOmu );
            hData_SameSign_Hemi_2Vtx_CutEvt_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_2Vtx_CutEvt_mass->Fill(   mass );
            hData_SameSign_Hemi_2Vtx_CutEvt_Mass->Fill(   VtxMass );
            hData_SameSign_Hemi_2Vtx_CutEvt_massMu->Fill( massMu );
            hData_SameSign_Hemi_2Vtx_CutEvt_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_2Vtx_CutEvt_dRMu->Fill(   dRMu );
            hData_SameSign_Hemi_2Vtx_CutEvt_Mmumu->Fill(  tree_MmumuSameSign );
            hData_SameSign_Hemi_2Vtx_CutEvt_ptmin->Fill(  hemi_ptmin );
            }

            if ( nVtx == 1 && isCutVtx ) {
            hData_SameSign_Hemi_1Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
            hData_SameSign_Hemi_1Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
            hData_SameSign_Hemi_1Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
            hData_SameSign_Hemi_1Vtx_CutVtx_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_1Vtx_CutVtx_mass->Fill(   mass );
            hData_SameSign_Hemi_1Vtx_CutVtx_Mass->Fill(   VtxMass );
            hData_SameSign_Hemi_1Vtx_CutVtx_massMu->Fill( massMu );
            hData_SameSign_Hemi_1Vtx_CutVtx_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_1Vtx_CutVtx_dRMu->Fill(   dRMu );
            hData_SameSign_Hemi_1Vtx_CutVtx_Mmumu->Fill(  tree_MmumuSameSign );
            hData_SameSign_Hemi_1Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
            }

            if ( nVtx == 2 && isCutVtx ) {
            hData_SameSign_Hemi_2Vtx_CutVtx_BDTevt->Fill( tree_Evts_MVAval );
            hData_SameSign_Hemi_2Vtx_CutVtx_BDTvtx->Fill( BDTvtx );
            hData_SameSign_Hemi_2Vtx_CutVtx_njetNOmu->Fill(   njetNOmu );
            hData_SameSign_Hemi_2Vtx_CutVtx_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_2Vtx_CutVtx_mass->Fill(   mass );
            hData_SameSign_Hemi_2Vtx_CutVtx_Mass->Fill(   VtxMass );
            hData_SameSign_Hemi_2Vtx_CutVtx_massMu->Fill( massMu );
            hData_SameSign_Hemi_2Vtx_CutVtx_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_2Vtx_CutVtx_dRMu->Fill(   dRMu );
            hData_SameSign_Hemi_2Vtx_CutVtx_Mmumu->Fill(  tree_MmumuSameSign );
            hData_SameSign_Hemi_2Vtx_CutVtx_ptmin->Fill(  hemi_ptmin );
            }

            if ( nVtx == 1 && isCutVtx && isCutEvt ) {
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_mass->Fill(   mass );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_MmumuSameSign );
            hData_SameSign_Hemi_1Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
            }

            if ( nVtx == 2 && isCutVtx && isCutEvt ) {
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_BDTvtx->Fill( BDTvtx );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_njetNOmu->Fill(   njetNOmu );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_nTrks->Fill(  nTrks );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_mass->Fill(   mass );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_Mass->Fill(   VtxMass );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_massMu->Fill( massMu );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_ptMu->Fill(   ptMu );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_dRMu->Fill(   dRMu );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_Mmumu->Fill(  tree_MmumuSameSign );
            hData_SameSign_Hemi_2Vtx_CutVtx_CutEvt_ptmin->Fill(  hemi_ptmin );
            }

         } // endif tree_FilterSameSign 
         //$$$$



   }// End of loop on Events


   NormFactor = XS/Nevent;
   // Output Postscript
   ofs<<"//---------------------------------------------------// "<<std::endl;
      ofs<<"Event yields for regions for ABCD 1 : "<<std::endl;
   ofs<<"Nevents: "                    <<Nevent*NormFactor*lumiRun2                     <<std::endl;
   ofs<<" Signal REGION : BDT EVT > "  <<EVTSWP<<" and Tight VTX (step 1+2)"              <<std::endl;
   ofs<<"Two Vertices  "               <<n_EVT2_SR_2Vtx*(NormFactor)*lumiRun2             <<std::endl;
   ofs<<"One Vertex  "                 <<n_EVT2_SR_1Vtx*(NormFactor)*lumiRun2             <<std::endl;
   ofs<<" Control REGION : BDT EVT < " <<EVTSWP<<" and Tight VTX (step 1+2)"              <<std::endl;
   ofs<<"Two vertices "                <<n_EVT2_CRNoEvt_2Vtx*(NormFactor)*lumiRun2        <<std::endl;
   ofs<<"One Vertex "                  <<n_EVT2_CRNoEvt_1Vtx*(NormFactor)*lumiRun2        <<std::endl;
   ofs<<" Control REGION : BDT EVT > " <<EVTSWP<<" and Loose VTX (step 3+4)"              <<std::endl;
   ofs<<"Two Vertices "                <<n_EVT2_CRNoVtx_1Vtx*(NormFactor)*lumiRun2        <<std::endl;
   ofs<<"One Vertex "                  <<n_EVT2_CRNoVtx_2Vtx*(NormFactor)*lumiRun2        <<std::endl;
   ofs<<" Control REGION : BDT EVT < " <<EVTSWP<<" and Loose VTX (step 3+4)"              <<std::endl;
   ofs<<"Two Vertices "                <<n_EVT2_CRNoEvtNoVtx_2Vtx*(NormFactor)*lumiRun2   <<std::endl;
   ofs<<"One Vertex "                  <<n_EVT2_CRNoEvtNoVtx_1Vtx*(NormFactor)*lumiRun2   <<std::endl;
   ofs<<"//---------------------------------------------------// "<<std::endl;
   ofs<<"//---------------------------------------------------// "<<std::endl;
   ofs<<"//---------------------------------------------------// "<<std::endl;
   ofs<<"Event yields for regions for ABCD 2 : "<<std::endl;
   ofs<<"Nevents: "                    <<Nevent*NormFactor*lumiRun2                  <<std::endl;
   ofs<<" Signal REGION : BDT EVT > "  <<EVTSWP<<" and BDT VTX > "<<VTXWP               <<std::endl;
   ofs<<"Two Vertices  "               <<n_EVT_SR_2Vtx*(NormFactor)*lumiRun2           <<std::endl;
   ofs<<"One Vertex  "                 <<n_EVT_SR_1Vtx*(NormFactor)*lumiRun2           <<std::endl;
   ofs<<" Control REGION : BDT EVT < " <<EVTSWP<<" and BDT VTX > "<<VTXWP               <<std::endl;
   ofs<<"Two vertices "                <<n_EVT_CRNoEvt_2Vtx*(NormFactor)*lumiRun2      <<std::endl;
   ofs<<"One Vertex "                  <<n_EVT_CRNoEvt_1Vtx*(NormFactor)*lumiRun2      <<std::endl;
   ofs<<" Control REGION : BDT EVT > " <<EVTSWP<<" and BDT VTX < "<<VTXWP               <<std::endl;
   ofs<<"Two Vertices "                <<n_EVT_CRNoVtx_1Vtx*(NormFactor)*lumiRun2      <<std::endl;
   ofs<<"One Vertex "                  <<n_EVT_CRNoVtx_2Vtx*(NormFactor)*lumiRun2      <<std::endl;
   ofs<<" Control REGION : BDT EVT < " <<EVTSWP<<" and BDT VTX < "<<VTXWP               <<std::endl;
   ofs<<"Two Vertices "                <<n_EVT_CRNoEvtNoVtx_2Vtx*(NormFactor)*lumiRun2 <<std::endl;
   ofs<<"One Vertex "                  <<n_EVT_CRNoEvtNoVtx_1Vtx*(NormFactor)*lumiRun2 <<std::endl;
   ofs<<" //----------------- End of Event yields ---------------\\ "<<std::endl;
   ofs<<" "<<std::endl;

   ofs.close();

   // std::cout << "number of events  "<< allevents << std::endl;
   HistogramManager h ;
  
   h.WriteAllHistogramsInFile(("outputroot/BKGEstimation_"+thesample+".root").Data(),"recreate");
}// En dof TreeAnalyzer::Loop
