#define ABCD_cxx
#include "ABCD.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "HistogramManager.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "TTree.h"
#include <TCanvas.h>


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
 if (!ahisto) { std::cout <<     "!! no histo !!" << std::endl; return ;}
 TDirectory* current = gDirectory ;
 afile->cd();
 ahisto->Write();
 current->cd();
}

void ABCD::Loop(TString sample, TString Production, bool Signal)
{
//   In a ROOT session, you can do:
//      root> .L ABCD.C
//      root> ABCD t
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
//$$
//**********************************
// Histograms
//**********************************
   TH1F* hData_Event_Weight = new TH1F("hData_Event_Weight","",200,0,2);
   TH1F* hData_Filter     = new TH1F("hData_Filter","",2,-0.5,1.5);

   TH1F* hData_Mmumu      = new TH1F("hData_Mmumu","",25,0.,500.);
   TH1F* hData_Evt_MVAval = new TH1F("hData_Evt_MVAval","",100,-1,1);
   //-----------------------------------------------------------//
   // ABCD using Evt and Tight+looseWP 
   //-----------------------------------------------------------//

   //------SR-1vtx----//
      TH1F* hData_EVT12_1Vtx_CutEvt_Mmumu    = new TH1F("hData_EVT12_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVT12_1Vtx_CutEvt_Mass     = new TH1F("hData_EVT12_1Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT12_1Vtx_BDTvtx          = new TH1F("hData_EVT12_1Vtx_BDTvtx","",25,-1,1);
      
      TH1F* hData_EVT12_1Vtx_HMass = new TH1F("hData_EVT12_1Vtx_HMass", "", 20,0,1000);
      TH1F* hData_EVT12_1Vtx_NChi2 = new TH1F("hData_EVT12_1Vtx_NChi2", "", 10, 0, 10);
      TH1F* hData_EVT12_1Vtx_nTrks = new TH1F("hData_EVT12_1Vtx_nTrks", "", 50, 0, 50);
      TH1F* hData_EVT12_1Vtx_z = new TH1F("hData_EVT12_1Vtx_z", "", 50, -50, 50);
      TH1F* hData_EVT12_1Vtx_r = new TH1F("hData_EVT12_1Vtx_r", "", 35, 0, 70);
      TH1F* hData_EVT12_1Vtx_dR = new TH1F("hData_EVT12_1Vtx_dR", "", 50, 0, 10);
      TH1F* hData_EVT12_1Vtx_SumtrackWeight = new TH1F("hData_EVT12_1Vtx_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_EVT12_1Vtx_MeantrackWeight = new TH1F("hData_EVT12_1Vtx_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_EVT12_1Vtx_track_MeanDCA_d = new TH1F("hData_EVT12_1Vtx_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_EVT12_1Vtx_dist = new TH1F("hData_EVT12_1Vtx_dist", "", 20, 0, 100);
   
   //------CR-----//
      TH1F* hData_NoEVT12_1Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVT12_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_Mass   = new TH1F("hData_NoEVT12_1Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT12_1Vtx_BDTvtx        = new TH1F("hData_NoEVT12_1Vtx_BDTvtx","",25,-1,1);

      TH1F* hData_NoEVT12_1Vtx_CutEvt_HMass = new TH1F("hData_NoEVT12_1Vtx_CutEvt_HMass", "", 20,0,1000);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_NChi2 = new TH1F("hData_NoEVT12_1Vtx_CutEvt_NChi2", "", 10, 0, 10);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_nTrks = new TH1F("hData_NoEVT12_1Vtx_CutEvt_nTrks", "", 50, 0, 50);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_z = new TH1F("hData_NoEVT12_1Vtx_CutEvt_z", "", 50, -50, 50);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_r = new TH1F("hData_NoEVT12_1Vtx_CutEvt_r", "", 35, 0, 70);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_dR = new TH1F("hData_NoEVT12_1Vtx_CutEvt_dR", "", 50, 0, 10);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_SumtrackWeight = new TH1F("hData_NoEVT12_1Vtx_CutEvt_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_MeantrackWeight = new TH1F("hData_NoEVT12_1Vtx_CutEvt_MeantrackWeight", "", 10, 0,1);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_track_MeanDCA_d = new TH1F("hData_NoEVT12_1Vtx_CutEvt_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_dist = new TH1F("hData_NoEVT12_1Vtx_CutEvt_dist", "", 20, 0, 100);

   //------CR-----//
      TH1F* hData_EVT34_1Vtx_CutEvt_Mmumu    = new TH1F("hData_EVT34_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVT34_1Vtx_CutEvt_Mass     = new TH1F("hData_EVT34_1Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT34_1Vtx_BDTvtx          = new TH1F("hData_EVT34_1Vtx_BDTvtx","",25,-1,1);

      TH1F* hData_EVT34_1Vtx_CutEvt_HMass = new TH1F("hData_EVT34_1Vtx_CutEvt_HMass", "", 20,0,1000);
      TH1F* hData_EVT34_1Vtx_CutEvt_NChi2 = new TH1F("hData_EVT34_1Vtx_CutEvt_NChi2", "", 10, 0, 10);
      TH1F* hData_EVT34_1Vtx_CutEvt_nTrks = new TH1F("hData_EVT34_1Vtx_CutEvt_nTrks", "", 50, 0, 50);
      TH1F* hData_EVT34_1Vtx_CutEvt_z = new TH1F("hData_EVT34_1Vtx_CutEvt_z", "", 50, -50, 50);
      TH1F* hData_EVT34_1Vtx_CutEvt_r = new TH1F("hData_EVT34_1Vtx_CutEvt_r", "", 35, 0, 70);
      TH1F* hData_EVT34_1Vtx_CutEvt_dR = new TH1F("hData_EVT34_1Vtx_CutEvt_dR", "", 50, 0, 10);
      TH1F* hData_EVT34_1Vtx_CutEvt_SumtrackWeight = new TH1F("hData_EVT34_1Vtx_CutEvt_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_EVT34_1Vtx_CutEvt_MeantrackWeight = new TH1F("hData_EVT34_1Vtx_CutEvt_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_EVT34_1Vtx_CutEvt_track_MeanDCA_d = new TH1F("hData_EVT34_1Vtx_CutEvt_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_EVT34_1Vtx_CutEvt_dist = new TH1F("hData_EVT34_1Vtx_CutEvt_dist", "", 20, 0, 100);


   //------CR-----//
      TH1F* hData_NoEVT34_1Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVT34_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_Mass   = new TH1F("hData_NoEVT34_1Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT34_1Vtx_BDTvtx        = new TH1F("hData_NoEVT34_1Vtx_BDTvtx","",100,-1,1);
   
      TH1F* hData_NoEVT34_1Vtx_CutEvt_HMass = new TH1F("hData_NoEVT34_1Vtx_CutEvt_HMass", "", 20,0,1000);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_NChi2 = new TH1F("hData_NoEVT34_1Vtx_CutEvt_NChi2", "", 10, 0, 10);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_nTrks = new TH1F("hData_NoEVT34_1Vtx_CutEvt_nTrks", "", 50, 0, 50);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_z = new TH1F("hData_NoEVT34_1Vtx_CutEvt_z", "", 50, -50, 50);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_r = new TH1F("hData_NoEVT34_1Vtx_CutEvt_r", "", 35, 0, 70);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_dR = new TH1F("hData_NoEVT34_1Vtx_CutEvt_dR", "", 50, 0, 10);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_SumtrackWeight = new TH1F("hData_NoEVT34_1Vtx_CutEvt_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_MeantrackWeight = new TH1F("hData_NoEVT34_1Vtx_CutEvt_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_track_MeanDCA_d = new TH1F("hData_NoEVT34_1Vtx_CutEvt_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_dist = new TH1F("hData_NoEVT34_1Vtx_CutEvt_dist", "", 20, 0, 100);

   //------SR-2Vtx----//
      TH1F* hData_EVT12_2Vtx_CutEvt_Mmumu    = new TH1F("hData_EVT12_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      
      TH1F* hData_EVT12_2Vtx_CutEvt_Mass     = new TH1F("hData_EVT12_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT12_2Vtx_BDTvtx          = new TH1F("hData_EVT12_2Vtx_BDTvtx","",25,-1,1);
      TH1F* hData_EVT12_2Vtx_HMass = new TH1F("hData_EVT12_2Vtx_HMass", "", 20,0,1000);
      TH1F* hData_EVT12_2Vtx_NChi2 = new TH1F("hData_EVT12_2Vtx_NChi2", "", 10, 0, 10);
      TH1F* hData_EVT12_2Vtx_nTrks = new TH1F("hData_EVT12_2Vtx_nTrks", "", 50, 0, 50);
      TH1F* hData_EVT12_2Vtx_z = new TH1F("hData_EVT12_2Vtx_z", "", 50, -50, 50);
      TH1F* hData_EVT12_2Vtx_r = new TH1F("hData_EVT12_2Vtx_r", "", 35, 0, 70);
      TH1F* hData_EVT12_2Vtx_dR = new TH1F("hData_EVT12_2Vtx_dR", "",50, 0, 10);
      TH1F* hData_EVT12_2Vtx_SumtrackWeight = new TH1F("hData_EVT12_2Vtx_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_EVT12_2Vtx_MeantrackWeight = new TH1F("hData_EVT12_2Vtx_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_EVT12_2Vtx_track_MeanDCA_d = new TH1F("hData_EVT12_2Vtx_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_EVT12_2Vtx_dist = new TH1F("hData_EVT12_2Vtx_dist", "", 20, 0, 100);
      TH1F* hData_EVT12_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_EVT12_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);

            TH1F* hData_EVT12_2VtxAll_CutEvt_Mass     = new TH1F("hData_EVT12_2VtxAll_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT12_2VtxAll_BDTvtx          = new TH1F("hData_EVT12_2VtxAll_BDTvtx","",25,-1,1);
      TH1F* hData_EVT12_2VtxAll_HMass = new TH1F("hData_EVT12_2VtxAll_HMass", "",20,0,1000);
      TH1F* hData_EVT12_2VtxAll_NChi2 = new TH1F("hData_EVT12_2VtxAll_NChi2", "", 10, 0, 10);
      TH1F* hData_EVT12_2VtxAll_nTrks = new TH1F("hData_EVT12_2VtxAll_nTrks", "", 50, 0, 50);
      TH1F* hData_EVT12_2VtxAll_z = new TH1F("hData_EVT12_2VtxAll_z", "", 50, -50, 50);
      TH1F* hData_EVT12_2VtxAll_r = new TH1F("hData_EVT12_2VtxAll_r", "", 35, 0, 70);
      TH1F* hData_EVT12_2VtxAll_dR = new TH1F("hData_EVT12_2VtxAll_dR", "", 50, 0, 10);
      TH1F* hData_EVT12_2VtxAll_SumtrackWeight = new TH1F("hData_EVT12_2VtxAll_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_EVT12_2VtxAll_MeantrackWeight = new TH1F("hData_EVT12_2VtxAll_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_EVT12_2VtxAll_track_MeanDCA_d = new TH1F("hData_EVT12_2VtxAll_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_EVT12_2VtxAll_dist = new TH1F("hData_EVT12_2VtxAll_dist", "",20, 0, 100);
      TH1F* hData_EVT12_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_EVT12_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);

   //------CR-----//
      TH1F* hData_NoEVT12_2Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVT12_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      
      TH1F* hData_NoEVT12_2Vtx_CutEvt_Mass   = new TH1F("hData_NoEVT12_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT12_2Vtx_BDTvtx        = new TH1F("hData_NoEVT12_2Vtx_BDTvtx","",25,-1,1);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_HMass = new TH1F("hData_NoEVT12_2Vtx_CutEvt_HMass", "",20,0,1000);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_NChi2 = new TH1F("hData_NoEVT12_2Vtx_CutEvt_NChi2", "", 10, 0, 10);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_nTrks = new TH1F("hData_NoEVT12_2Vtx_CutEvt_nTrks", "", 50, 0, 50);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_z = new TH1F("hData_NoEVT12_2Vtx_CutEvt_z", "",50, -50, 50);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_r = new TH1F("hData_NoEVT12_2Vtx_CutEvt_r", "", 35, 0, 70);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_dR = new TH1F("hData_NoEVT12_2Vtx_CutEvt_dR", "", 100, 0, 10);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_SumtrackWeight = new TH1F("hData_NoEVT12_2Vtx_CutEvt_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_MeantrackWeight = new TH1F("hData_NoEVT12_2Vtx_CutEvt_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_track_MeanDCA_d = new TH1F("hData_NoEVT12_2Vtx_CutEvt_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_dist = new TH1F("hData_NoEVT12_2Vtx_CutEvt_dist", "", 20, 0, 100);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_NoEVT12_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);
            TH1F* hData_NoEVT12_2VtxAll_CutEvt_Mass   = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT12_2VtxAll_BDTvtx        = new TH1F("hData_NoEVT12_2VtxAll_BDTvtx","",100,-1,1);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_HMass = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_HMass", "", 20,0,1000);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_NChi2 = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_NChi2", "", 10, 0, 10);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_nTrks = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_nTrks", "", 50, 0, 50);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_z = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_z", "", 50, -50, 50);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_r = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_r", "", 35, 0, 70);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_dR = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_dR", "", 50, 0, 10);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_SumtrackWeight = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_MeantrackWeight = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_track_MeanDCA_d = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_dist = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_dist", "", 20, 0, 100);
      TH1F* hData_NoEVT12_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_NoEVT12_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);

   //------CR-----//
      TH1F* hData_EVT34_2Vtx_CutEvt_Mmumu    = new TH1F("hData_EVT34_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      
      TH1F* hData_EVT34_2Vtx_CutEvt_Mass     = new TH1F("hData_EVT34_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT34_2Vtx_BDTvtx          = new TH1F("hData_EVT34_2Vtx_BDTvtx","",25,-1,1);
      TH1F* hData_EVT34_2Vtx_CutEvt_HMass = new TH1F("hData_EVT34_2Vtx_CutEvt_HMass", "", 20,0,1000);
      TH1F* hData_EVT34_2Vtx_CutEvt_NChi2 = new TH1F("hData_EVT34_2Vtx_CutEvt_NChi2", "", 10, 0, 10);
      TH1F* hData_EVT34_2Vtx_CutEvt_nTrks = new TH1F("hData_EVT34_2Vtx_CutEvt_nTrks", "", 50, 0, 50);
      TH1F* hData_EVT34_2Vtx_CutEvt_z = new TH1F("hData_EVT34_2Vtx_CutEvt_z", "", 50, -50, 50);
      TH1F* hData_EVT34_2Vtx_CutEvt_r = new TH1F("hData_EVT34_2Vtx_CutEvt_r", "", 35, 0, 70);
      TH1F* hData_EVT34_2Vtx_CutEvt_dR = new TH1F("hData_EVT34_2Vtx_CutEvt_dR", "", 50, 0, 10);
      TH1F* hData_EVT34_2Vtx_CutEvt_SumtrackWeight = new TH1F("hData_EVT34_2Vtx_CutEvt_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_EVT34_2Vtx_CutEvt_MeantrackWeight = new TH1F("hData_EVT34_2Vtx_CutEvt_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_EVT34_2Vtx_CutEvt_track_MeanDCA_d = new TH1F("hData_EVT34_2Vtx_CutEvt_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_EVT34_2Vtx_CutEvt_dist = new TH1F("hData_EVT34_2Vtx_CutEvt_dist", "", 20, 0, 100);
      TH1F* hData_EVT34_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_EVT34_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);

            TH1F* hData_EVT34_2VtxAll_CutEvt_Mass     = new TH1F("hData_EVT34_2VtxAll_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT34_2VtxAll_BDTvtx          = new TH1F("hData_EVT34_2VtxAll_BDTvtx","",100,-1,1);
      TH1F* hData_EVT34_2VtxAll_CutEvt_HMass = new TH1F("hData_EVT34_2VtxAll_CutEvt_HMass", "", 20,0,1000);
      TH1F* hData_EVT34_2VtxAll_CutEvt_NChi2 = new TH1F("hData_EVT34_2VtxAll_CutEvt_NChi2", "", 10, 0, 10);
      TH1F* hData_EVT34_2VtxAll_CutEvt_nTrks = new TH1F("hData_EVT34_2VtxAll_CutEvt_nTrks", "", 50, 0, 50);
      TH1F* hData_EVT34_2VtxAll_CutEvt_z = new TH1F("hData_EVT34_2VtxAll_CutEvt_z", "", 50, -50, 50);
      TH1F* hData_EVT34_2VtxAll_CutEvt_r = new TH1F("hData_EVT34_2VtxAll_CutEvt_r", "", 35, 0, 70);
      TH1F* hData_EVT34_2VtxAll_CutEvt_dR = new TH1F("hData_EVT34_2VtxAll_CutEvt_dR", "", 50, 0, 10);
      TH1F* hData_EVT34_2VtxAll_CutEvt_SumtrackWeight = new TH1F("hData_EVT34_2VtxAll_CutEvt_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_EVT34_2VtxAll_CutEvt_MeantrackWeight = new TH1F("hData_EVT34_2VtxAll_CutEvt_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_EVT34_2VtxAll_CutEvt_track_MeanDCA_d = new TH1F("hData_EVT34_2VtxAll_CutEvt_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_EVT34_2VtxAll_CutEvt_dist = new TH1F("hData_EVT34_2VtxAll_CutEvt_dist", "", 20, 0, 100);
      TH1F* hData_EVT34_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_EVT34_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);

   //------CR-----//
      TH1F* hData_NoEVT34_2Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVT34_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      
      TH1F* hData_NoEVT34_2Vtx_CutEvt_Mass   = new TH1F("hData_NoEVT34_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT34_2Vtx_BDTvtx        = new TH1F("hData_NoEVT34_2Vtx_BDTvtx","",25,-1,1);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_HMass = new TH1F("hData_NoEVT34_2Vtx_CutEvt_HMass", "", 20,0,1000);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_NChi2 = new TH1F("hData_NoEVT34_2Vtx_CutEvt_NChi2", "", 10, 0, 10);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_nTrks = new TH1F("hData_NoEVT34_2Vtx_CutEvt_nTrks", "", 50, 0, 50);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_z = new TH1F("hData_NoEVT34_2Vtx_CutEvt_z", "", 50, -50, 50);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_r = new TH1F("hData_NoEVT34_2Vtx_CutEvt_r", "", 35, 0, 70);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_dR = new TH1F("hData_NoEVT34_2Vtx_CutEvt_dR", "", 50, 0, 10);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_SumtrackWeight = new TH1F("hData_NoEVT34_2Vtx_CutEvt_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_MeantrackWeight = new TH1F("hData_NoEVT34_2Vtx_CutEvt_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_track_MeanDCA_d = new TH1F("hData_NoEVT34_2Vtx_CutEvt_track_MeanDCA_d", "", 25, 0, 25);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_dist = new TH1F("hData_NoEVT34_2Vtx_CutEvt_dist", "", 20, 0, 100);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_NoEVT34_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);

            TH1F* hData_NoEVT34_2VtxAll_CutEvt_Mass   = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT34_2VtxAll_BDTvtx        = new TH1F("hData_NoEVT34_2VtxAll_BDTvtx","",25,-1,1);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_HMass = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_HMass", "",20,0,1000);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_NChi2 = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_NChi2", "", 10, 0, 10);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_nTrks = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_nTrks", "", 50, 0, 50);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_z = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_z", "", 50, -50, 50);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_r = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_r", "", 35, 0, 70);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_dR = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_dR", "", 50, 0, 10);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_SumtrackWeight = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_MeantrackWeight = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_MeantrackWeight", "", 10, 0, 1);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_track_MeanDCA_d = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_track_MeanDCA_d", "",25, 0, 25);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_dist = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_dist", "", 20, 0, 100);
      TH1F* hData_NoEVT34_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_NoEVT34_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);

   //-----------------------------------------------------------//
   // ABCD using Evt and Vtx BDT
   //-----------------------------------------------------------//
      TH1F* hData_1Vtx_MVAval                = new TH1F("hData_1Vtx_MVAval","",100,-1,1);
      TH1F* hData_2Vtx_MVAval                = new TH1F("hData_2Vtx_MVAval","",100,-1,1);
      TH1F* hData_2VtxAll_MVAval             = new TH1F("hData_2VtxAll_MVAval","",100,-1,1);


   //------SR-1vtx----//
      TH1F* hData_EVTVtx_1Vtx_CutEvt_Mmumu   = new TH1F("hData_EVTVtx_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTVtx_1Vtx_CutEvt_Mass    = new TH1F("hData_EVTVtx_1Vtx_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_NoEVTVtx_1Vtx_CutEvt_Mmumu = new TH1F("hData_NoEVTVtx_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTVtx_1Vtx_CutEvt_Mass  = new TH1F("hData_NoEVTVtx_1Vtx_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_EVTNoVtx_1Vtx_CutEvt_Mmumu = new TH1F("hData_EVTNoVtx_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTNoVtx_1Vtx_CutEvt_Mass  = new TH1F("hData_EVTNoVtx_1Vtx_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_NoEVTNoVtx_1Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVTNoVtx_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTNoVtx_1Vtx_CutEvt_Mass   = new TH1F("hData_NoEVTNoVtx_1Vtx_CutEvt_Mass","",25,0.,100.);

   //------SR-2Vtx----//
      TH1F* hData_EVTVtx_2Vtx_CutEvt_Mmumu      = new TH1F("hData_EVTVtx_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTVtx_2Vtx_CutEvt_Mass       = new TH1F("hData_EVTVtx_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVTVtx_2VtxAll_CutEvt_Mmumu   = new TH1F("hData_EVTVtx_2VtxAll_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTVtx_2VtxAll_CutEvt_Mass    = new TH1F("hData_EVTVtx_2VtxAll_CutEvt_Mass","",25,0.,100.);
         
   //------CR-----//
      TH1F* hData_NoEVTVtx_2Vtx_CutEvt_Mmumu    = new TH1F("hData_NoEVTVtx_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTVtx_2Vtx_CutEvt_Mass     = new TH1F("hData_NoEVTVtx_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVTVtx_2VtxAll_CutEvt_Mmumu = new TH1F("hData_NoEVTVtx_2VtxAll_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTVtx_2VtxAll_CutEvt_Mass  = new TH1F("hData_NoEVTVtx_2VtxAll_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_EVTNoVtx_2Vtx_CutEvt_Mmumu    = new TH1F("hData_EVTNoVtx_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTNoVtx_2Vtx_CutEvt_Mass     = new TH1F("hData_EVTNoVtx_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVTNoVtx_2VtxAll_CutEvt_Mmumu = new TH1F("hData_EVTNoVtx_2VtxAll_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTNoVtx_2VtxAll_CutEvt_Mass  = new TH1F("hData_EVTNoVtx_2VtxAll_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_NoEVTNoVtx_2Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVTNoVtx_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTNoVtx_2Vtx_CutEvt_Mass   = new TH1F("hData_NoEVTNoVtx_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVTNoVtx_2VtxAll_CutEvt_Mmumu  = new TH1F("hData_NoEVTNoVtx_2VtxAll_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTNoVtx_2VtxAll_CutEvt_Mass   = new TH1F("hData_NoEVTNoVtx_2VtxAll_CutEvt_Mass","",25,0.,100.);
         
   //-----------------------------------------------------------//
   // ABCD using Hemisphere pt and Tight+looseWP
   //-----------------------------------------------------------//


//------SR-----//
// add other vtx variables


   TH1F* hData_Hemi_BDTevt = new TH1F("hData_Hemi_BDTevt", "", 25, -1, 1);
   TH1F* hData_Hemi_0Vtx_CutEvt_Mmumu = new TH1F("hData_Hemi_0Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_Hemi_0Vtx_BDTevt = new TH1F("hData_Hemi_0Vtx_BDTevt", "", 25, -1, 1);


   TH1F* hData_Hemi_1Vtx_CutEvt_Mmumu = new TH1F("hData_Hemi_1Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_Hemi_1Vtx_CutEvt_Mass = new TH1F("hData_Hemi_1Vtx_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_Hemi_1Vtx_BDTevt = new TH1F("hData_Hemi_1Vtx_BDTevt", "", 25, -1, 1);
   TH1F* hData_Hemi_1Vtx_BDTvtx = new TH1F("hData_Hemi_1Vtx_BDTvtx", "", 25, -1, 1);
      TH1F* hData_Hemi_1Vtx_HMass = new TH1F("hData_Hemi_1Vtx_HMass", "", 20,0,1000);
   TH1F* hData_Hemi_1Vtx_NChi2 = new TH1F("hData_Hemi_1Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_Hemi_1Vtx_nTrks = new TH1F("hData_Hemi_1Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_Hemi_1Vtx_z = new TH1F("hData_Hemi_1Vtx_z", "", 50, -50, 50);
   TH1F* hData_Hemi_1Vtx_r = new TH1F("hData_Hemi_1Vtx_r", "", 35, 0, 70);
   TH1F* hData_Hemi_1Vtx_dR = new TH1F("hData_Hemi_1Vtx_dR", "", 50, 0, 10);
   TH1F* hData_Hemi_1Vtx_SumtrackWeight = new TH1F("hData_Hemi_1Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_Hemi_1Vtx_MeantrackWeight = new TH1F("hData_Hemi_1Vtx_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_Hemi_1Vtx_track_MeanDCA_d = new TH1F("hData_Hemi_1Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_Hemi_1Vtx_dist = new TH1F("hData_Hemi_1Vtx_dist", "", 20, 0, 100);

   TH1F* hData_Hemi_2Vtx_CutEvt_Mmumu = new TH1F("hData_Hemi_2Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_Hemi_2Vtx_CutEvt_Mass = new TH1F("hData_Hemi_2Vtx_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_Hemi_2Vtx_CutEvt_BDTevt = new TH1F("hData_Hemi_2Vtx_CutEvt_BDTevt", "", 25, -1, 1);
   TH1F* hData_Hemi_2Vtx_CutEvt_MaxBDTvtx = new TH1F("hData_Hemi_2Vtx_CutEvt_MaxBDTvtx", "", 25, -1, 1);
   TH1F* hData_Hemi_2VtxAll_CutEvt_Mass = new TH1F("hData_Hemi_2VtxAll_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_Hemi_2VtxAll_CutEvt_BDTvtx = new TH1F("hData_Hemi_2VtxAll_CutEvt_BDTvtx", "", 25, -1, 1);
      TH1F* hData_Hemi_2Vtx_HMass = new TH1F("hData_Hemi_2Vtx_HMass", "", 20,0,1000);
   TH1F* hData_Hemi_2Vtx_NChi2 = new TH1F("hData_Hemi_2Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_Hemi_2Vtx_nTrks = new TH1F("hData_Hemi_2Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_Hemi_2Vtx_z = new TH1F("hData_Hemi_2Vtx_z", "", 50, -50, 50);
   TH1F* hData_Hemi_2Vtx_r = new TH1F("hData_Hemi_2Vtx_r", "", 35, 0, 70);
   TH1F* hData_Hemi_2Vtx_dR = new TH1F("hData_Hemi_2Vtx_dR", "", 50, 0, 10);
   TH1F* hData_Hemi_2Vtx_SumtrackWeight = new TH1F("hData_Hemi_2Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_Hemi_2Vtx_MeantrackWeight = new TH1F("hData_Hemi_2Vtx_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_Hemi_2Vtx_track_MeanDCA_d = new TH1F("hData_Hemi_2Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_Hemi_2Vtx_dist = new TH1F("hData_Hemi_2Vtx_dist", "", 20, 0, 100);
   TH1F* hData_Hemi_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_Hemi_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);
         TH1F* hData_Hemi_2VtxAll_HMass = new TH1F("hData_Hemi_2VtxAll_HMass", "", 20,0,1000);
   TH1F* hData_Hemi_2VtxAll_NChi2 = new TH1F("hData_Hemi_2VtxAll_NChi2", "", 10, 0, 10);
   TH1F* hData_Hemi_2VtxAll_nTrks = new TH1F("hData_Hemi_2VtxAll_nTrks", "", 50, 0, 50);
   TH1F* hData_Hemi_2VtxAll_z = new TH1F("hData_Hemi_2VtxAll_z", "", 50, -50, 50);
   TH1F* hData_Hemi_2VtxAll_r = new TH1F("hData_Hemi_2VtxAll_r", "", 35, 0, 70);
   TH1F* hData_Hemi_2VtxAll_dR = new TH1F("hData_Hemi_2VtxAll_dR", "", 50, 0, 10);
   TH1F* hData_Hemi_2VtxAll_SumtrackWeight = new TH1F("hData_Hemi_2VtxAll_SumtrackWeight", "",25, 0, 25);
   TH1F* hData_Hemi_2VtxAll_MeantrackWeight = new TH1F("hData_Hemi_2VtxAll_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_Hemi_2VtxAll_track_MeanDCA_d = new TH1F("hData_Hemi_2VtxAll_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_Hemi_2VtxAll_dist = new TH1F("hData_Hemi_2VtxAll_dist", "", 20, 0, 100);
   TH1F* hData_Hemi_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_Hemi_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);

//------low pt-----//
   TH1F* hData_CRlowpt_BDTevt              = new TH1F("hData_CRlowpt_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRlowpt_0Vtx_CutEvt_Mmumu   = new TH1F("hData_CRlowpt_0Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRlowpt_0Vtx_BDTevt         = new TH1F("hData_CRlowpt_0Vtx_BDTevt", "", 25, -1, 1);

   TH1F* hData_CRlowpt_1Vtx_CutEvt_Mass    = new TH1F("hData_CRlowpt_1Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRlowpt_1Vtx_CutEvt_Mmumu   = new TH1F("hData_CRlowpt_1Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRlowpt_1Vtx_BDTevt= new TH1F("hData_CRlowpt_1Vtx_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRlowpt_1Vtx_BDTvtx= new TH1F("hData_CRlowpt_1Vtx_BDTvtx", "", 25, -1, 1);
         TH1F* hData_CRlowpt_1Vtx_HMass = new TH1F("hData_CRlowpt_1Vtx_HMass", "", 20,0,1000);
   TH1F* hData_CRlowpt_1Vtx_NChi2 = new TH1F("hData_CRlowpt_1Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRlowpt_1Vtx_nTrks = new TH1F("hData_CRlowpt_1Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRlowpt_1Vtx_z = new TH1F("hData_CRlowpt_1Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRlowpt_1Vtx_r = new TH1F("hData_CRlowpt_1Vtx_r", "", 35, 0, 70);
   TH1F* hData_CRlowpt_1Vtx_dR = new TH1F("hData_CRlowpt_1Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRlowpt_1Vtx_SumtrackWeight = new TH1F("hData_CRlowpt_1Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRlowpt_1Vtx_MeantrackWeight = new TH1F("hData_CRlowpt_1Vtx_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRlowpt_1Vtx_track_MeanDCA_d = new TH1F("hData_CRlowpt_1Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRlowpt_1Vtx_dist = new TH1F("hData_CRlowpt_1Vtx_dist", "", 20, 0, 100);

   TH1F* hData_CRlowpt_2Vtx_CutEvt_Mass    = new TH1F("hData_CRlowpt_2Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRlowpt_2Vtx_CutEvt_Mmumu   = new TH1F("hData_CRlowpt_2Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRlowpt_2Vtx_CutEvt_BDTevt= new TH1F("hData_CRlowpt_2Vtx_CutEvt_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRlowpt_2Vtx_CutEvt_MaxBDTvtx= new TH1F("hData_CRlowpt_2Vtx_CutEvt_MaxBDTvtx", "",25, -1, 1);
   TH1F* hData_CRlowpt_2VtxAll_CutEvt_Mass = new TH1F("hData_CRlowpt_2VtxAll_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRlowpt_2VtxAll_CutEvt_BDTvtx= new TH1F("hData_CRlowpt_2VtxAll_CutEvt_BDTvtx", "", 25, -1, 1);
      TH1F* hData_CRlowpt_2Vtx_HMass = new TH1F("hData_CRlowpt_2Vtx_HMass", "", 20,0,1000);
   TH1F* hData_CRlowpt_2Vtx_NChi2 = new TH1F("hData_CRlowpt_2Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRlowpt_2Vtx_nTrks = new TH1F("hData_CRlowpt_2Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRlowpt_2Vtx_z = new TH1F("hData_CRlowpt_2Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRlowpt_2Vtx_r = new TH1F("hData_CRlowpt_2Vtx_r", "",35, 0, 70);
   TH1F* hData_CRlowpt_2Vtx_dR = new TH1F("hData_CRlowpt_2Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRlowpt_2Vtx_SumtrackWeight = new TH1F("hData_CRlowpt_2Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRlowpt_2Vtx_MeantrackWeight = new TH1F("hData_CRlowpt_2Vtx_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRlowpt_2Vtx_track_MeanDCA_d = new TH1F("hData_CRlowpt_2Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRlowpt_2Vtx_dist = new TH1F("hData_CRlowpt_2Vtx_dist", "", 20, 0, 100);
   TH1F* hData_CRlowpt_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_CRlowpt_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);
      TH1F* hData_CRlowpt_2VtxAll_HMass = new TH1F("hData_CRlowpt_2VtxAll_HMass", "",20,0,1000);
   TH1F* hData_CRlowpt_2VtxAll_NChi2 = new TH1F("hData_CRlowpt_2VtxAll_NChi2", "", 10, 0, 10);
   TH1F* hData_CRlowpt_2VtxAll_nTrks = new TH1F("hData_CRlowpt_2VtxAll_nTrks", "", 50, 0, 50);
   TH1F* hData_CRlowpt_2VtxAll_z = new TH1F("hData_CRlowpt_2VtxAll_z", "", 50, -50, 50);
   TH1F* hData_CRlowpt_2VtxAll_r = new TH1F("hData_CRlowpt_2VtxAll_r", "",35, 0, 70);
   TH1F* hData_CRlowpt_2VtxAll_dR = new TH1F("hData_CRlowpt_2VtxAll_dR", "", 50, 0, 10);
   TH1F* hData_CRlowpt_2VtxAll_SumtrackWeight = new TH1F("hData_CRlowpt_2VtxAll_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRlowpt_2VtxAll_MeantrackWeight = new TH1F("hData_CRlowpt_2VtxAll_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRlowpt_2VtxAll_track_MeanDCA_d = new TH1F("hData_CRlowpt_2VtxAll_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRlowpt_2VtxAll_dist = new TH1F("hData_CRlowpt_2VtxAll_dist", "", 20, 0, 100);
   TH1F* hData_CRlowpt_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_CRlowpt_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);
//------loose -----//

   TH1F* hData_CRloose_BDTevt = new TH1F("hData_CRloose_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRloose_0Vtx_CutEvt_Mmumu = new TH1F("hData_CRloose_0Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_CRloose_0Vtx_BDTevt = new TH1F("hData_CRloose_0Vtx_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRloose_1Vtx_CutEvt_Mmumu = new TH1F("hData_CRloose_1Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_CRloose_1Vtx_CutEvt_Mass = new TH1F("hData_CRloose_1Vtx_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_CRloose_1Vtx_BDTevt = new TH1F("hData_CRloose_1Vtx_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRloose_1Vtx_BDTvtx = new TH1F("hData_CRloose_1Vtx_BDTvtx", "", 25, -1, 1);
            TH1F* hData_CRloose_1Vtx_HMass = new TH1F("hData_CRloose_1Vtx_HMass", "", 20,0,1000);
   TH1F* hData_CRloose_1Vtx_NChi2 = new TH1F("hData_CRloose_1Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRloose_1Vtx_nTrks = new TH1F("hData_CRloose_1Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRloose_1Vtx_z = new TH1F("hData_CRloose_1Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRloose_1Vtx_r = new TH1F("hData_CRloose_1Vtx_r", "", 35, 0, 70);
   TH1F* hData_CRloose_1Vtx_dR = new TH1F("hData_CRloose_1Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRloose_1Vtx_SumtrackWeight = new TH1F("hData_CRloose_1Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRloose_1Vtx_MeantrackWeight = new TH1F("hData_CRloose_1Vtx_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRloose_1Vtx_track_MeanDCA_d = new TH1F("hData_CRloose_1Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRloose_1Vtx_dist = new TH1F("hData_CRloose_1Vtx_dist", "", 20, 0, 100);

   TH1F* hData_CRloose_2Vtx_CutEvt_Mmumu = new TH1F("hData_CRloose_2Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_CRloose_2Vtx_CutEvt_Mass = new TH1F("hData_CRloose_2Vtx_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_CRloose_2Vtx_CutEvt_BDTevt = new TH1F("hData_CRloose_2Vtx_CutEvt_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRloose_2Vtx_CutEvt_MaxBDTvtx = new TH1F("hData_CRloose_2Vtx_CutEvt_MaxBDTvtx", "", 25, -1, 1);
   TH1F* hData_CRloose_2VtxAll_CutEvt_Mass = new TH1F("hData_CRloose_2VtxAll_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_CRloose_2VtxAll_CutEvt_BDTvtx = new TH1F("hData_CRloose_2VtxAll_CutEvt_BDTvtx", "", 25, -1, 1);
         TH1F* hData_CRloose_2Vtx_HMass = new TH1F("hData_CRloose_2Vtx_HMass", "", 20,0,1000);
   TH1F* hData_CRloose_2Vtx_NChi2 = new TH1F("hData_CRloose_2Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRloose_2Vtx_nTrks = new TH1F("hData_CRloose_2Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRloose_2Vtx_z = new TH1F("hData_CRloose_2Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRloose_2Vtx_r = new TH1F("hData_CRloose_2Vtx_r", "", 35, 0, 70);
   TH1F* hData_CRloose_2Vtx_dR = new TH1F("hData_CRloose_2Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRloose_2Vtx_SumtrackWeight = new TH1F("hData_CRloose_2Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRloose_2Vtx_MeantrackWeight = new TH1F("hData_CRloose_2Vtx_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRloose_2Vtx_track_MeanDCA_d = new TH1F("hData_CRloose_2Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRloose_2Vtx_dist = new TH1F("hData_CRloose_2Vtx_dist", "", 20, 0, 100);
   TH1F* hData_CRloose_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_CRloose_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);
            TH1F* hData_CRloose_2VtxAll_HMass = new TH1F("hData_CRloose_2VtxAll_HMass", "", 20,0,1000);
   TH1F* hData_CRloose_2VtxAll_NChi2 = new TH1F("hData_CRloose_2VtxAll_NChi2", "", 10, 0, 10);
   TH1F* hData_CRloose_2VtxAll_nTrks = new TH1F("hData_CRloose_2VtxAll_nTrks", "", 50, 0, 50);
   TH1F* hData_CRloose_2VtxAll_z = new TH1F("hData_CRloose_2VtxAll_z", "", 50, -50, 50);
   TH1F* hData_CRloose_2VtxAll_r = new TH1F("hData_CRloose_2VtxAll_r", "", 35, 0, 70);
   TH1F* hData_CRloose_2VtxAll_dR = new TH1F("hData_CRloose_2VtxAll_dR", "", 50, 0, 10);
   TH1F* hData_CRloose_2VtxAll_SumtrackWeight = new TH1F("hData_CRloose_2VtxAll_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRloose_2VtxAll_MeantrackWeight = new TH1F("hData_CRloose_2VtxAll_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRloose_2VtxAll_track_MeanDCA_d = new TH1F("hData_CRloose_2VtxAll_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRloose_2VtxAll_dist = new TH1F("hData_CRloose_2VtxAll_dist", "", 20, 0, 100);
   TH1F* hData_CRloose_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_CRloose_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);

//------loose low pt -----//

   TH1F* hData_CRlooselowpt_BDTevt = new TH1F("hData_CRlooselowpt_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRlooselowpt_0Vtx_CutEvt_Mmumu = new TH1F("hData_CRlooselowpt_0Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_CRlooselowpt_1Vtx_BDTevt = new TH1F("hData_CRlooselowpt_1Vtx_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRlooselowpt_1Vtx_BDTvtx = new TH1F("hData_CRlooselowpt_1Vtx_BDTvtx", "", 25, -1, 1);
   TH1F* hData_CRlooselowpt_1Vtx_CutEvt_Mmumu = new TH1F("hData_CRlooselowpt_1Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_CRlooselowpt_1Vtx_CutEvt_Mass = new TH1F("hData_CRlooselowpt_1Vtx_CutEvt_Mass", "", 25, 0., 100.);
      TH1F* hData_CRlooselowpt_1Vtx_HMass = new TH1F("hData_CRlooselowpt_1Vtx_HMass", "", 20,0,1000);
   TH1F* hData_CRlooselowpt_1Vtx_NChi2 = new TH1F("hData_CRlooselowpt_1Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRlooselowpt_1Vtx_nTrks = new TH1F("hData_CRlooselowpt_1Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRlooselowpt_1Vtx_z = new TH1F("hData_CRlooselowpt_1Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRlooselowpt_1Vtx_r = new TH1F("hData_CRlooselowpt_1Vtx_r", "", 35, 0, 70);
   TH1F* hData_CRlooselowpt_1Vtx_dR = new TH1F("hData_CRlooselowpt_1Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRlooselowpt_1Vtx_SumtrackWeight = new TH1F("hData_CRlooselowpt_1Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRlooselowpt_1Vtx_MeantrackWeight = new TH1F("hData_CRlooselowpt_1Vtx_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRlooselowpt_1Vtx_track_MeanDCA_d = new TH1F("hData_CRlooselowpt_1Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRlooselowpt_1Vtx_dist = new TH1F("hData_CRlooselowpt_1Vtx_dist", "", 20, 0, 100);


   TH1F* hData_CRlooselowpt_2Vtx_CutEvt_BDTevt = new TH1F("hData_CRlooselowpt_2Vtx_CutEvt_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRlooselowpt_2Vtx_CutEvt_MaxBDTvtx = new TH1F("hData_CRlooselowpt_2Vtx_CutEvt_MaxBDTvtx", "",25, -1, 1);
   TH1F* hData_CRlooselowpt_2Vtx_CutEvt_Mmumu = new TH1F("hData_CRlooselowpt_2Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_CRlooselowpt_2Vtx_CutEvt_Mass = new TH1F("hData_CRlooselowpt_2Vtx_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_CRlooselowpt_2VtxAll_CutEvt_Mass = new TH1F("hData_CRlooselowpt_2VtxAll_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_CRlooselowpt_2VtxAll_CutEvt_BDTvtx = new TH1F("hData_CRlooselowpt_2VtxAll_CutEvt_BDTvtx", "",25, -1, 1);
      TH1F* hData_CRlooselowpt_2Vtx_HMass = new TH1F("hData_CRlooselowpt_2Vtx_HMass", "",20,0,1000);
   TH1F* hData_CRlooselowpt_2Vtx_NChi2 = new TH1F("hData_CRlooselowpt_2Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRlooselowpt_2Vtx_nTrks = new TH1F("hData_CRlooselowpt_2Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRlooselowpt_2Vtx_z = new TH1F("hData_CRlooselowpt_2Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRlooselowpt_2Vtx_r = new TH1F("hData_CRlooselowpt_2Vtx_r", "", 35, 0, 70);
   TH1F* hData_CRlooselowpt_2Vtx_dR = new TH1F("hData_CRlooselowpt_2Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRlooselowpt_2Vtx_SumtrackWeight = new TH1F("hData_CRlooselowpt_2Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRlooselowpt_2Vtx_MeantrackWeight = new TH1F("hData_CRlooselowpt_2Vtx_MeantrackWeight", "", 10, 0, 1);

   TH1F* hData_CRlooselowpt_2Vtx_track_MeanDCA_d = new TH1F("hData_CRlooselowpt_2Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRlooselowpt_2Vtx_dist = new TH1F("hData_CRlooselowpt_2Vtx_dist", "", 20, 0, 100);
   TH1F* hData_CRlooselowpt_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_CRlooselowpt_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);
         TH1F* hData_CRlooselowpt_2VtxAll_HMass = new TH1F("hData_CRlooselowpt_2VtxAll_HMass", "", 20,0,1000);
   TH1F* hData_CRlooselowpt_2VtxAll_NChi2 = new TH1F("hData_CRlooselowpt_2VtxAll_NChi2", "", 10, 0, 10);
   TH1F* hData_CRlooselowpt_2VtxAll_nTrks = new TH1F("hData_CRlooselowpt_2VtxAll_nTrks", "", 50, 0, 50);
   TH1F* hData_CRlooselowpt_2VtxAll_z = new TH1F("hData_CRlooselowpt_2VtxAll_z", "", 50, -50, 50);
   TH1F* hData_CRlooselowpt_2VtxAll_r = new TH1F("hData_CRlooselowpt_2VtxAll_r", "",35, 0, 70);
   TH1F* hData_CRlooselowpt_2VtxAll_dR = new TH1F("hData_CRlooselowpt_2VtxAll_dR", "", 50, 0, 10);
   TH1F* hData_CRlooselowpt_2VtxAll_SumtrackWeight = new TH1F("hData_CRlooselowpt_2VtxAll_SumtrackWeight", "",25, 0, 25);
   TH1F* hData_CRlooselowpt_2VtxAll_MeantrackWeight = new TH1F("hData_CRlooselowpt_2VtxAll_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRlooselowpt_2VtxAll_track_MeanDCA_d = new TH1F("hData_CRlooselowpt_2VtxAll_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRlooselowpt_2VtxAll_dist = new TH1F("hData_CRlooselowpt_2VtxAll_dist", "", 20, 0, 100);
   TH1F* hData_CRlooselowpt_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_CRlooselowpt_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);

//------ fwd -----//

   TH1F* hData_CRfwd_BDTevt = new TH1F("hData_CRfwd_BDTevt", "",25, -1, 1);
   TH1F* hData_CRfwd_0Vtx_CutEvt_Mmumu = new TH1F("hData_CRfwd_0Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_CRfwd_0Vtx_BDTevt = new TH1F("hData_CRfwd_0Vtx_BDTevt", "",25, -1, 1);
   TH1F* hData_CRfwd_1Vtx_CutEvt_Mmumu = new TH1F("hData_CRfwd_1Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_CRfwd_1Vtx_CutEvt_Mass = new TH1F("hData_CRfwd_1Vtx_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_CRfwd_1Vtx_BDTevt = new TH1F("hData_CRfwd_1Vtx_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRfwd_1Vtx_BDTvtx = new TH1F("hData_CRfwd_1Vtx_BDTvtx", "", 25, -1, 1);
      TH1F* hData_CRfwd_1Vtx_HMass = new TH1F("hData_CRfwd_1Vtx_HMass", "", 20,0,1000);
   TH1F* hData_CRfwd_1Vtx_NChi2 = new TH1F("hData_CRfwd_1Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRfwd_1Vtx_nTrks = new TH1F("hData_CRfwd_1Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRfwd_1Vtx_z = new TH1F("hData_CRfwd_1Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRfwd_1Vtx_r = new TH1F("hData_CRfwd_1Vtx_r", "", 35, 0, 70);
   TH1F* hData_CRfwd_1Vtx_dR = new TH1F("hData_CRfwd_1Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRfwd_1Vtx_SumtrackWeight = new TH1F("hData_CRfwd_1Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRfwd_1Vtx_MeantrackWeight = new TH1F("hData_CRfwd_1Vtx_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRfwd_1Vtx_track_MeanDCA_d = new TH1F("hData_CRfwd_1Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRfwd_1Vtx_dist = new TH1F("hData_CRfwd_1Vtx_dist", "", 20, 0, 100);

   TH1F* hData_CRfwd_2Vtx_CutEvt_Mmumu = new TH1F("hData_CRfwd_2Vtx_CutEvt_Mmumu", "", 25, 0., 500.);
   TH1F* hData_CRfwd_2Vtx_CutEvt_Mass = new TH1F("hData_CRfwd_2Vtx_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_CRfwd_2Vtx_CutEvt_BDTevt = new TH1F("hData_CRfwd_2Vtx_CutEvt_BDTevt", "", 25, -1, 1);
   TH1F* hData_CRfwd_2Vtx_CutEvt_MaxBDTvtx = new TH1F("hData_CRfwd_2Vtx_CutEvt_MaxBDTvtx", "", 25, -1, 1);
   TH1F* hData_CRfwd_2VtxAll_CutEvt_Mass = new TH1F("hData_CRfwd_2VtxAll_CutEvt_Mass", "", 25, 0., 100.);
   TH1F* hData_CRfwd_2VtxAll_CutEvt_BDTvtx = new TH1F("hData_CRfwd_2VtxAll_CutEvt_BDTvtx", "", 25, -1, 1);
         TH1F* hData_CRfwd_2Vtx_HMass = new TH1F("hData_CRfwd_2Vtx_HMass", "", 20,0,1000);
   TH1F* hData_CRfwd_2Vtx_NChi2 = new TH1F("hData_CRfwd_2Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRfwd_2Vtx_nTrks = new TH1F("hData_CRfwd_2Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRfwd_2Vtx_z = new TH1F("hData_CRfwd_2Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRfwd_2Vtx_r = new TH1F("hData_CRfwd_2Vtx_r", "", 35, 0, 70);
   TH1F* hData_CRfwd_2Vtx_dR = new TH1F("hData_CRfwd_2Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRfwd_2Vtx_SumtrackWeight = new TH1F("hData_CRfwd_2Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRfwd_2Vtx_MeantrackWeight = new TH1F("hData_CRfwd_2Vtx_MeantrackWeight", "", 10, 0, 1);

   TH1F* hData_CRfwd_2Vtx_track_MeanDCA_d = new TH1F("hData_CRfwd_2Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRfwd_2Vtx_dist = new TH1F("hData_CRfwd_2Vtx_dist", "", 20, 0, 100);
   TH1F* hData_CRfwd_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_CRfwd_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);
            TH1F* hData_CRfwd_2VtxAll_HMass = new TH1F("hData_CRfwd_2VtxAll_HMass", "", 20,0,1000);
   TH1F* hData_CRfwd_2VtxAll_NChi2 = new TH1F("hData_CRfwd_2VtxAll_NChi2", "", 10, 0, 10);
   TH1F* hData_CRfwd_2VtxAll_nTrks = new TH1F("hData_CRfwd_2VtxAll_nTrks", "", 50, 0, 50);
   TH1F* hData_CRfwd_2VtxAll_z = new TH1F("hData_CRfwd_2VtxAll_z", "", 50, -50, 50);
   TH1F* hData_CRfwd_2VtxAll_r = new TH1F("hData_CRfwd_2VtxAll_r", "", 35, 0, 70);
   TH1F* hData_CRfwd_2VtxAll_dR = new TH1F("hData_CRfwd_2VtxAll_dR", "", 50, 0, 10);
   TH1F* hData_CRfwd_2VtxAll_SumtrackWeight = new TH1F("hData_CRfwd_2VtxAll_SumtrackWeight", "", 25, 0, 25);
      TH1F* hData_CRfwd_2VtxAll_MeantrackWeight = new TH1F("hData_CRfwd_2VtxAll_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRfwd_2VtxAll_track_MeanDCA_d = new TH1F("hData_CRfwd_2VtxAll_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRfwd_2VtxAll_dist = new TH1F("hData_CRfwd_2VtxAll_dist", "", 20, 0, 100);
   TH1F* hData_CRfwd_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_CRfwd_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);

   //----- Samesign ----//  
   TH1F* hData_SameSign_BDTevt               = new TH1F("hData_SameSign_BDTevt","",25, -1, 1);
   TH1F* hData_CRss_0Vtx_CutEvt_Mmumu        = new TH1F("hData_CRss_0Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRss_0Vtx_BDTevt              = new TH1F("hData_CRss_0Vtx_BDTevt","",25, -1, 1);
   TH1F* hData_CRss_1Vtx_CutEvt_Mass         = new TH1F("hData_CRss_1Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRss_1Vtx_CutEvt_Mmumu        = new TH1F("hData_CRss_1Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRss_1Vtx_BDTevt              = new TH1F("hData_CRss_1Vtx_BDTevt","",25, -1, 1);
   TH1F* hData_CRss_1Vtx_BDTvtx              = new TH1F("hData_CRss_1Vtx_BDTvtx","",25, -1, 1);
      TH1F* hData_CRss_1Vtx_HMass = new TH1F("hData_CRss_1Vtx_HMass", "",20,0,1000);
   TH1F* hData_CRss_1Vtx_NChi2 = new TH1F("hData_CRss_1Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRss_1Vtx_nTrks = new TH1F("hData_CRss_1Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRss_1Vtx_z = new TH1F("hData_CRss_1Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRss_1Vtx_r = new TH1F("hData_CRss_1Vtx_r", "", 35, 0, 70);
   TH1F* hData_CRss_1Vtx_dR = new TH1F("hData_CRss_1Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRss_1Vtx_SumtrackWeight = new TH1F("hData_CRss_1Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRss_1Vtx_MeantrackWeight = new TH1F("hData_CRs_1Vtx_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRss_1Vtx_track_MeanDCA_d = new TH1F("hData_CRss_1Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRss_1Vtx_dist = new TH1F("hData_CRss_1Vtx_dist", "", 20, 0, 100);

   TH1F* hData_CRss_2Vtx_CutEvt_Mass         = new TH1F("hData_CRss_2Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRss_2Vtx_CutEvt_Mmumu        = new TH1F("hData_CRss_2Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRss_2Vtx_CutEvt_BDTevt       = new TH1F("hData_CRss_2Vtx_CutEvt_BDTevt","",25, -1, 1);
   TH1F* hData_CRss_2Vtx_CutEvt_MaxBDTvtx    = new TH1F("hData_CRss_2Vtx_CutEvt_MaxBDTvtx","",25, -1, 1);
   TH1F* hData_CRss_2VtxAll_CutEvt_Mass      = new TH1F("hData_CRss_2VtxAll_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRss_2VtxAll_CutEvt_BDTvtx    = new TH1F("hData_CRss_2VtxAll_CutEvt_BDTvtx","",25, -1, 1);
            TH1F* hData_CRss_2Vtx_HMass = new TH1F("hData_CRss_2Vtx_HMass", "",20,0,1000);
   TH1F* hData_CRss_2Vtx_NChi2 = new TH1F("hData_CRss_2Vtx_NChi2", "", 10, 0, 10);
   TH1F* hData_CRss_2Vtx_nTrks = new TH1F("hData_CRss_2Vtx_nTrks", "", 50, 0, 50);
   TH1F* hData_CRss_2Vtx_z = new TH1F("hData_CRss_2Vtx_z", "", 50, -50, 50);
   TH1F* hData_CRss_2Vtx_r = new TH1F("hData_CRss_2Vtx_r", "", 35, 0, 70);
   TH1F* hData_CRss_2Vtx_dR = new TH1F("hData_CRss_2Vtx_dR", "", 50, 0, 10);
   TH1F* hData_CRss_2Vtx_SumtrackWeight = new TH1F("hData_CRss_2Vtx_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRss_2Vtx_MeantrackWeight = new TH1F("hData_CRs_2Vtx_MeantrackWeight", "", 10, 0, 1);

   TH1F* hData_CRss_2Vtx_track_MeanDCA_d = new TH1F("hData_CRss_2Vtx_track_MeanDCA_d", "", 25, 0, 25);
   TH1F* hData_CRss_2Vtx_dist = new TH1F("hData_CRss_2Vtx_dist", "", 20, 0, 100);
   TH1F* hData_CRss_2Vtx_CutEvt_VtxVtxdist = new TH1F("hData_CRss_2Vtx_CutEvt_VtxVtxdist", "", 80, 0, 400);
               TH1F* hData_CRss_2VtxAll_HMass = new TH1F("hData_CRss_2VtxAll_HMass", "", 20,0,1000);
   TH1F* hData_CRss_2VtxAll_NChi2 = new TH1F("hData_CRss_2VtxAll_NChi2", "", 10, 0, 10);
   TH1F* hData_CRss_2VtxAll_nTrks = new TH1F("hData_CRss_2VtxAll_nTrks", "", 50, 0, 50);
   TH1F* hData_CRss_2VtxAll_z = new TH1F("hData_CRss_2VtxAll_z", "", 50, -50, 50);
   TH1F* hData_CRss_2VtxAll_r = new TH1F("hData_CRss_2VtxAll_r", "", 35, 0, 70);
   TH1F* hData_CRss_2VtxAll_dR = new TH1F("hData_CRss_2VtxAll_dR", "", 50, 0, 10);
   TH1F* hData_CRss_2VtxAll_SumtrackWeight = new TH1F("hData_CRss_2VtxAll_SumtrackWeight", "", 25, 0, 25);
   TH1F* hData_CRss_2VtxAll_MeantrackWeight = new TH1F("hData_CRs_2VtxAll_MeantrackWeight", "", 10, 0, 1);
   TH1F* hData_CRss_2VtxAll_track_MeanDCA_d = new TH1F("hData_CRss_2VtxAll_track_MeanDCA_d", "",25, 0, 25);
   TH1F* hData_CRss_2VtxAll_dist = new TH1F("hData_CRss_2VtxAll_dist", "", 20, 0, 100);
   TH1F* hData_CRss_2VtxAll_CutEvt_VtxVtxdist = new TH1F("hData_CRss_2VtxAll_CutEvt_VtxVtxdist", "", 80, 0, 400);

///////////////////////////////////////////////////////////////////


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

  float XS = 1;
   TString thesample  = sample;
   // XS are given in pb
   if (thesample.Contains("DYJetsToLL_M-10to50"))                    { XS = 15910.0;   }
   if (thesample.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays")) { XS = 10.8707;   }
   if (thesample.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     { XS = 10.8908;   }
   if (thesample.Contains("TTJets_DiLept"))                          { XS = 53.07;     }
   if (thesample.Contains("WWTo2L2Nu"))                              { XS = 11.09;     }
   if (thesample.Contains("WZTo2Q2L_mllmin4p0"))                     { XS = 6.535;     }
   if (thesample.Contains("ZZTo2Q2L_mllmin4p0"))                     { XS = 3.676;     }
   if (thesample.Contains("TTTo2L2Nu"))                              { XS = 88.3;      }
   if (thesample.Contains("DYJetsToLL_M-50"))                        { XS = 5379;      }
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

   bool BlindSR = false;
   float Nevent = 0;
   bool signal = Signal;
   double Pref_PU_gen_wt=1;
   //   if (nentries>760000){nentries =760000;}

   cout<< "Line : "  << __LINE__ << " " << nentries << endl; 
   cout<< " XS : "<<XS<<endl;
   Long64_t nbytes = 0, nb = 0;
   int allevents = 0;
   // nentries = 100000;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      allevents++;
      if ( allevents%1000 == 0 ) std::cout << "events : " << allevents << std::endl;

       if ( signal &&  minitree_nLLP->at(0) != 2 ) continue; // protection against rare wrong signal events    // 
  
      //--------------------------------------------------------------//
      if (jentry==0.1*nentries) {std::cout<<"10/100 :"<<std::endl;}
      if (jentry==0.2*nentries) {std::cout<<"20/100 :"<<std::endl;}
      if (jentry==0.3*nentries) {std::cout<<"30/100 :"<<std::endl;}
      if (jentry==0.4*nentries) {std::cout<<"40/100 :"<<std::endl;}
      if (jentry==0.5*nentries) {std::cout<<"50/100 :"<<std::endl;}
      if (jentry==0.6*nentries) {std::cout<<"60/100 :"<<std::endl;}
      if (jentry==0.7*nentries) {std::cout<<"70/100 :"<<std::endl;}
      if (jentry==0.8*nentries) {std::cout<<"80/100 :"<<std::endl;}
      if (jentry==0.9*nentries) {std::cout<<"90/100 :"<<std::endl;}

    
      Pref_PU_gen_wt=minitree_only_gen_wt->at(0)*miniPrefweight->at(0)*miniPUweight->at(0);
      if (signal) Pref_PU_gen_wt = miniPrefweight->at(0)*miniPUweight->at(0);
      hData_Event_Weight->Fill( Pref_PU_gen_wt );
      hData_Filter->Fill( minitree_Filter->at(0),Pref_PU_gen_wt );
      //--------------------------------------------------------------//
      bool isHemiVtx1 = false, isHemiVtx2 = false;
      bool isCutVtx = false, isCutVtx1 = false, isCutVtx2 = false;
      bool isCutEvt = false;
      float BDTvtx = -2., BDTvtx1 = -2., BDTvtx2 = -2.;
      bool ping;
      bool isHemiVtx1Loose = false, isHemiVtx2Loose = false;
      int nVtx = 0, nVtxIni = 0, step;
      int nVtxLoose = 0 , nVtxIniLoose = 0; 
      float VtxMass = 0., dR, dist, NChi2, r, eta;
      float Vtx_HMass = 0;
      float Vtx_nTrks= 0;
      float Vtx_z= 0;
      float Vtx_r= 0;
      float Vtx_dR= 0;
      float Vtx_SumtrackWeight= 0;
      float Vtx_track_MeanDCA_d= 0;
      float Vtx_dist = 0;
      float Vtx_NChi = 0;

      bool ping0 = false;
      bool ping1 = false;
      float dR0 = 0.;
      float dR1 = 0.;
      ////////////////// hemisphere pT
      float  hemi1_pt  = -1.;
      float  hemi2_pt  = -1.;
      ////////////////////////////////

      hData_Mmumu->Fill( minitree_Mmumu->at(0),Pref_PU_gen_wt );
      
      if ( !minitree_Filter->at(0) && !minitree_FilterSameSign->at(0) ) continue;//

      if (minitree_njetNOmu->at(0) < 1) continue;

      bool Filter = minitree_Filter->at(0); // Tight ID and Mini IsoTight for both muons for DImuon Channel
      if ( minitree_njetNOmu->at(0) < 1 ) Filter = false; 
      hemi1_pt = minitree_Hemi_pt->at(0);
      hemi2_pt = minitree_Hemi_pt->at(1);

      float hemi_ptmin = hemi2_pt;
      if ( hemi2_pt > hemi1_pt ) hemi_ptmin = hemi1_pt;


   if (signal)
      {
         ping0 = minitree_Hemi_LLP_ping->at(0);
         ping1 = minitree_Hemi_LLP_ping->at(1);
         dR0 = minitree_Hemi_LLP_dR->at(0);
         dR1 = minitree_Hemi_LLP_dR->at(1);
      }

   int Vtx_step0 = minitree_Hemi_Vtx_step->at(0);
   int Vtx_step1 = minitree_Hemi_Vtx_step->at(1);
   float Vtx_NChi0 = minitree_Hemi_Vtx_NChi2->at(0);
   float Vtx_NChi1 = minitree_Hemi_Vtx_NChi2->at(1);
   float Vtx_Mass0 = minitree_Hemi_Vtx_Mass->at(0);
   float Vtx_Mass1 = minitree_Hemi_Vtx_Mass->at(1);
   float Vtx_dist0 = minitree_Hemi_Vtx_dist->at(0);
   float Vtx_dist1 = minitree_Hemi_Vtx_dist->at(1);

   float Vtx_HMass0 = minitree_Hemi_Vtx_BDT_HMass->at(0);
   float Vtx_HMass1 = minitree_Hemi_Vtx_BDT_HMass->at(1);
   float Vtx_nTrks0 = minitree_Hemi_Vtx_nTrks->at(0);
   float Vtx_nTrks1 = minitree_Hemi_Vtx_nTrks->at(1);
   float Vtx_z0 = minitree_Hemi_Vtx_z->at(0);
   float Vtx_z1 = minitree_Hemi_Vtx_z->at(1);
   float Vtx_r0 = minitree_Hemi_Vtx_r->at(0);
   float Vtx_r1 = minitree_Hemi_Vtx_r->at(1);
   float Vtx_dR0 = minitree_Hemi_Vtx_dR->at(0);
   float Vtx_dR1 = minitree_Hemi_Vtx_dR->at(1);
   float Vtx_SumtrackWeight0 = minitree_Hemi_Vtx_SumtrackWeight->at(0);
   float Vtx_SumtrackWeight1 = minitree_Hemi_Vtx_SumtrackWeight->at(1);

   float Vtx_MeantrackWeight0 = 0;
   float Vtx_MeantrackWeight1 = 0;
   float Vtx_MeantrackWeight = 0;
   if ( minitree_Hemi_Vtx_nTrks->at(0) != 0){Vtx_MeantrackWeight0 = minitree_Hemi_Vtx_SumtrackWeight->at(0)/static_cast< float >(minitree_Hemi_Vtx_nTrks->at(0));}
   if ( minitree_Hemi_Vtx_nTrks->at(1) != 0){Vtx_MeantrackWeight1 = minitree_Hemi_Vtx_SumtrackWeight->at(1)/static_cast< float >(minitree_Hemi_Vtx_nTrks->at(1));}


   float Vtx_track_MeanDCA_d0 = minitree_Hemi_Vtx_track_MeanDCA_d->at(0);
   float Vtx_track_MeanDCA_d1 = minitree_Hemi_Vtx_track_MeanDCA_d->at(1);
   float Vtx_Vtx_dist = 0;


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
   float recX0 = posx0 - minitree_PV_x->at(0);
   float recY0 = posy0 - minitree_PV_y->at(0);
   float recZ0 = posz0 - minitree_PV_z->at(0);
   float recX1 = posx1 - minitree_PV_x->at(0);
   float recY1 = posy1 - minitree_PV_y->at(0);
   float recZ1 = posz1 - minitree_PV_z->at(0);

  float theta_Vtx0 = TMath::ATan2(sqrt(recX0*recX0+recY0*recY0),abs(recZ0)) ;
  float theta_Vtx1 = TMath::ATan2(sqrt(recX1*recX1+recY1*recY1),abs(recZ1)) ;

  float eta_Vtx0 = -TMath::Log(tan(theta_Vtx0/2));
  float eta_Vtx1 = -TMath::Log(tan(theta_Vtx1/2));
  if ( posz0 < 0 ) eta_Vtx0 = -eta_Vtx0;
  if ( posz1 < 0 ) eta_Vtx1 = -eta_Vtx1;

   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 1 && Vtx_step0 <= 2 ) {
     isHemiVtx1 = true;
   }
   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 3 && Vtx_step0 <= 4 ) {
     isHemiVtx1Loose = true;
   }
   if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 1 && Vtx_step1 <= 2 ) {
     isHemiVtx2 = true;
   }
   if ( Vtx_NChi1 > 0 && Vtx_NChi1< 10 && Vtx_step1 >= 3 && Vtx_step1 <= 4 ) {
     isHemiVtx2Loose = true;
   }
   if      ( isHemiVtx1 && isHemiVtx2 ) nVtxIni = 2;
   else if ( isHemiVtx1 || isHemiVtx2 ) nVtxIni = 1;

   if      ( isHemiVtx1Loose && isHemiVtx2Loose ) nVtxIniLoose = 2;
   else if ( isHemiVtx1Loose || isHemiVtx2Loose ) nVtxIniLoose = 1;

   bool Merging = true; // false if "no merge" or "no close vtx" output
   bool Protect = false;
   //$$
   if (signal)
      {
         if (minitree_Hemi_SecLLP_ping->size() >= 1) Protect = true;
      }
   else
      {
         Protect = true;
      } 
    if ( Merging && Protect ) { 
         // protection again
     if ( minitree_Hemi_SecVtx->size() >= 1 ) {
         if (signal)
            {ping0 = minitree_Hemi_SecLLP_ping->at(0);}
         ping0 = false;
         Vtx_step0 = minitree_Hemi_SecVtx_step->at(0);
         Vtx_NChi0 = minitree_Hemi_SecVtx_NChi2->at(0);
         Vtx_Mass0 = minitree_Hemi_SecVtx_Mass->at(0);
         Vtx_dist0 = minitree_Hemi_SecVtx_dist->at(0);
         posx0 = minitree_Hemi_SecVtx_x->at(0);
         posy0 = minitree_Hemi_SecVtx_y->at(0);
         posz0 = minitree_Hemi_SecVtx_z->at(0);
         r0 = minitree_Hemi_SecVtx_r->at(0);
         Vtx_step1 = 0;
         Vtx_NChi1 = -1.;
         Vtx_Mass1 = 0.;
         Vtx_dist1 = 0.;
         ping1 = false;
         r1 = 0;
         float SecrecX0 = posx0 - minitree_PV_x->at(0);
         float SecrecY0 = posy0 - minitree_PV_y->at(0);
         float SecrecZ0 = posz0 - minitree_PV_z->at(0);
         float theta_SecVtx0 = TMath::ATan2(sqrt(SecrecX0*SecrecX0+SecrecY0*SecrecY0),abs(SecrecZ0)) ;

         eta_Vtx0 = -TMath::Log(tan(theta_SecVtx0/2.));
         eta_Vtx1 = 0;
         if ( posz0 < 0 ) eta_Vtx0 = -eta_Vtx0;
         Vtx_nTrks0 = minitree_Hemi_SecVtx_nTrks->at(0);
         Vtx_z0 = minitree_Hemi_SecVtx_z->at(0);
         Vtx_r0 = minitree_Hemi_SecVtx_r->at(0);
         Vtx_dR0 = minitree_Hemi_SecVtx_dR->at(0);
         Vtx_SumtrackWeight0 = minitree_Hemi_SecVtx_SumtrackWeight->at(0);
         if (minitree_Hemi_SecVtx_nTrks->at(0) !=0) Vtx_MeantrackWeight0 = minitree_Hemi_SecVtx_SumtrackWeight->at(0)/static_cast< float >(minitree_Hemi_SecVtx_nTrks->at(0));
         Vtx_track_MeanDCA_d0 = minitree_Hemi_SecVtx_track_MeanDCA_d->at(0);
     }
     if ( minitree_Hemi_SecVtx->size() == 2 ) {
         if (signal)
            {ping1 = minitree_Hemi_SecLLP_ping->at(1);}
         ping1 = false;
         Vtx_step1 = minitree_Hemi_SecVtx_step->at(1);
         Vtx_NChi1 = minitree_Hemi_SecVtx_NChi2->at(1);
         Vtx_Mass1 = minitree_Hemi_SecVtx_Mass->at(1);
         Vtx_dist1 = minitree_Hemi_SecVtx_dist->at(1);
         posx1 = minitree_Hemi_SecVtx_x->at(1);
         posy1 = minitree_Hemi_SecVtx_y->at(1);
         posz1 = minitree_Hemi_SecVtx_z->at(1);
         r1 = minitree_Hemi_SecVtx_r->at(1);
         float SecrecX1 = posx1 - minitree_PV_x->at(0);
         float SecrecY1 = posy1 - minitree_PV_y->at(0);
         float SecrecZ1 = posz1 - minitree_PV_z->at(0);
         float theta_SecVtx1 = TMath::ATan2(sqrt(SecrecX1*SecrecX1+SecrecY1*SecrecY1),abs(SecrecZ1)) ;
         eta_Vtx1 = -TMath::Log(tan(theta_SecVtx1/2.));
         if ( posz1 < 0 ) eta_Vtx1 = -eta_Vtx1;
         Vtx_nTrks1 = minitree_Hemi_SecVtx_nTrks->at(1);
         Vtx_z1 = minitree_Hemi_SecVtx_z->at(1);
         Vtx_r1 = minitree_Hemi_SecVtx_r->at(1);
         Vtx_dR1 = minitree_Hemi_SecVtx_dR->at(1);
         Vtx_SumtrackWeight1 = minitree_Hemi_SecVtx_SumtrackWeight->at(1);
         if (minitree_Hemi_SecVtx_nTrks->at(1) !=0) Vtx_MeantrackWeight1 = minitree_Hemi_SecVtx_SumtrackWeight->at(1)/static_cast< float >(minitree_Hemi_SecVtx_nTrks->at(1));
         Vtx_track_MeanDCA_d1 = minitree_Hemi_SecVtx_track_MeanDCA_d->at(1);

     }
   }// End of Merging information
   isHemiVtx1 = false;
   isHemiVtx2 = false;

   isHemiVtx1Loose = false;
   isHemiVtx2Loose = false;        
   Vtx_Vtx_dist = sqrt((posx1-posx0)*(posx1-posx0)+(posy1-posy0)*(posy1-posy0)+
                        (posz1-posz0)*(posz1-posz0));

      //First Vertex
   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 1 && Vtx_step0 <= 2 ) {
      isHemiVtx1 = true;
      VtxMass = Vtx_Mass0;
      BDTvtx1 = minitree_Hemi_Vtx_MVAval_Tight->at(0);
      BDTvtx  = BDTvtx1;

      Vtx_NChi = Vtx_NChi0;
      Vtx_dist = Vtx_dist0;
      Vtx_HMass = Vtx_HMass0;
      Vtx_nTrks = Vtx_nTrks0;
      Vtx_z = Vtx_z0;
      Vtx_r = Vtx_r0;
      Vtx_dR = Vtx_dR0;
      Vtx_SumtrackWeight = Vtx_SumtrackWeight0;
      Vtx_MeantrackWeight = Vtx_MeantrackWeight0;
      Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d0;

   }
   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 3 && Vtx_step0 <= 4 ) {
      isHemiVtx1Loose = true;
      VtxMass = Vtx_Mass0;
      BDTvtx1 = minitree_Hemi_Vtx_MVAval_Loose->at(0);
      BDTvtx  = BDTvtx1;

      Vtx_NChi = Vtx_NChi0;
      Vtx_dist = Vtx_dist0;
      Vtx_HMass = Vtx_HMass0;
      Vtx_nTrks = Vtx_nTrks0;
      Vtx_z = Vtx_z0;
      Vtx_r = Vtx_r0;
      Vtx_dR = Vtx_dR0;
      Vtx_SumtrackWeight = Vtx_SumtrackWeight0;
      Vtx_MeantrackWeight = Vtx_MeantrackWeight0;
      Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d0;
   }

   //Second Vertex
   if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 1 && Vtx_step1 <= 2 ) {
      isHemiVtx2 = true;
      if ( Vtx_Mass1 > VtxMass ) VtxMass = Vtx_Mass1;
      BDTvtx2 = minitree_Hemi_Vtx_MVAval_Tight->at(1);

      if ( BDTvtx2 > BDTvtx )       BDTvtx      = BDTvtx2;
      if ( Vtx_NChi1 < Vtx_NChi ) Vtx_NChi      = Vtx_NChi1;
      if ( Vtx_dist1 > Vtx_dist ) Vtx_dist      = Vtx_dist1;
      if ( Vtx_HMass1 > Vtx_HMass ) Vtx_HMass   = Vtx_HMass1;
      if ( Vtx_nTrks1 > Vtx_nTrks ) Vtx_nTrks      = Vtx_nTrks1;
      if ( Vtx_r1 > Vtx_r )     Vtx_r           = Vtx_r1;
      if ( Vtx_z1 > Vtx_z )     Vtx_z           = Vtx_z1;
      if ( Vtx_dR1 > Vtx_dR )    Vtx_dR         = Vtx_dR1;
      if ( Vtx_SumtrackWeight1 > Vtx_SumtrackWeight )    Vtx_SumtrackWeight = Vtx_SumtrackWeight1;
      if ( Vtx_MeantrackWeight1 > Vtx_MeantrackWeight )    Vtx_MeantrackWeight = Vtx_MeantrackWeight1;
      if ( Vtx_track_MeanDCA_d1 < Vtx_track_MeanDCA_d ) Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d1;
   }

   if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 3 && Vtx_step1 <= 4 ) {
      isHemiVtx2Loose = true;
      if ( Vtx_Mass1 > VtxMass ) VtxMass = Vtx_Mass1;
      BDTvtx2 = minitree_Hemi_Vtx_MVAval_Loose->at(1);

      if ( BDTvtx2 > BDTvtx )       BDTvtx      = BDTvtx2;
      if ( Vtx_NChi1 < Vtx_NChi ) Vtx_NChi      = Vtx_NChi1;
      if ( Vtx_dist1 > Vtx_dist ) Vtx_dist      = Vtx_dist1;
      if ( Vtx_HMass1 > Vtx_HMass ) Vtx_HMass   = Vtx_HMass1;
      if ( Vtx_nTrks1 > Vtx_nTrks ) Vtx_nTrks      = Vtx_nTrks1;
      if ( Vtx_r1 > Vtx_r )     Vtx_r           = Vtx_r1;
      if ( Vtx_z1 > Vtx_z )     Vtx_z           = Vtx_z1;
      if ( Vtx_dR1 > Vtx_dR )    Vtx_dR         = Vtx_dR1;
      if ( Vtx_SumtrackWeight1 > Vtx_SumtrackWeight )    Vtx_SumtrackWeight = Vtx_SumtrackWeight1;
      if ( Vtx_MeantrackWeight1 > Vtx_MeantrackWeight )    Vtx_MeantrackWeight = Vtx_MeantrackWeight1;
      if ( Vtx_track_MeanDCA_d1 < Vtx_track_MeanDCA_d ) Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d1;
   }


   if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;
   else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;

   if      ( isHemiVtx1Loose && isHemiVtx2Loose ) nVtxLoose = 2;
   else if ( isHemiVtx1Loose || isHemiVtx2Loose ) nVtxLoose = 1;

   int nHemi = minitree_Hemi->size();//_LLP_dR
   for (int i = 0; i < nHemi; i++) 
      {   // Loop on hemispheres
         if ( i == 0 ) {
            dR    = dR0;
            dist  = Vtx_dist0;
            NChi2 = Vtx_NChi0;
            step  = Vtx_step0;
            ping  = ping0;
            r = r0;
            eta = eta_Vtx0;
         }
         else if ( i == 1 ) {
            dR    = dR1;
            dist  = Vtx_dist1;
            NChi2 = Vtx_NChi1;
            step  = Vtx_step1;
            ping  = ping1;
            r = r1;
            eta= eta_Vtx1;
         }
      }// end loop on Hemi

   //-----------------------------------------------------------//
   // ABCD using Evt and Tight+looseWP 
   //-----------------------------------------------------------//
   float EVTSWP = 0.85;
   float VTXWP = 0.85;
   hData_Evt_MVAval->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
   //$$
   if (minitree_Evts_MVAval->at(0) > EVTSWP) isCutEvt = true;

   //$$
   //----------------------//
   if (nVtx == 1 && isCutEvt ) {
      hData_EVT12_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_EVT12_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_BDTvtx->Fill( BDTvtx ,Pref_PU_gen_wt);

      hData_EVT12_1Vtx_HMass->Fill( Vtx_HMass ,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_NChi2->Fill( Vtx_NChi ,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_nTrks->Fill( Vtx_nTrks ,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_z->Fill( Vtx_z ,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_r->Fill( Vtx_r ,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_dR->Fill( Vtx_dR ,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight ,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight ,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d ,Pref_PU_gen_wt);
      hData_EVT12_1Vtx_dist->Fill( Vtx_dist ,Pref_PU_gen_wt);
   }

   if (nVtx == 1 && !isCutEvt) {
      hData_NoEVT12_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_BDTvtx->Fill( BDTvtx ,Pref_PU_gen_wt);

      hData_NoEVT12_1Vtx_CutEvt_HMass->Fill( Vtx_HMass ,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_NChi2->Fill( Vtx_NChi ,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_nTrks->Fill( Vtx_nTrks ,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_z->Fill( Vtx_z ,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_r->Fill( Vtx_r ,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_dR->Fill( Vtx_dR ,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight ,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight ,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d ,Pref_PU_gen_wt);
      hData_NoEVT12_1Vtx_CutEvt_dist->Fill( Vtx_dist ,Pref_PU_gen_wt);
   }

   if (nVtxLoose == 1 && isCutEvt) {
      hData_EVT34_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_EVT34_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      hData_EVT34_1Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );
      
      hData_EVT34_1Vtx_CutEvt_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
      hData_EVT34_1Vtx_CutEvt_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
      hData_EVT34_1Vtx_CutEvt_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
      hData_EVT34_1Vtx_CutEvt_z->Fill( Vtx_z,Pref_PU_gen_wt );
      hData_EVT34_1Vtx_CutEvt_r->Fill( Vtx_r,Pref_PU_gen_wt );
      hData_EVT34_1Vtx_CutEvt_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
      hData_EVT34_1Vtx_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
      hData_EVT34_1Vtx_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
      hData_EVT34_1Vtx_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
      hData_EVT34_1Vtx_CutEvt_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
   }

   if (nVtxLoose == 1 && !isCutEvt) {
      hData_NoEVT34_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_NoEVT34_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      hData_NoEVT34_1Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

      hData_NoEVT34_1Vtx_CutEvt_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
      hData_NoEVT34_1Vtx_CutEvt_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
      hData_NoEVT34_1Vtx_CutEvt_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
      hData_NoEVT34_1Vtx_CutEvt_z->Fill( Vtx_z,Pref_PU_gen_wt );
      hData_NoEVT34_1Vtx_CutEvt_r->Fill( Vtx_r,Pref_PU_gen_wt );
      hData_NoEVT34_1Vtx_CutEvt_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
      hData_NoEVT34_1Vtx_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
      hData_NoEVT34_1Vtx_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
      hData_NoEVT34_1Vtx_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
      hData_NoEVT34_1Vtx_CutEvt_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
   }
   //----------------------//

   if (nVtx == 2 && isCutEvt) {
      hData_EVT12_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);

      hData_EVT12_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      hData_EVT12_2Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
      hData_EVT12_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );

            hData_EVT12_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass0,Pref_PU_gen_wt);
      hData_EVT12_2VtxAll_BDTvtx->Fill( 1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_HMass->Fill( Vtx_HMass0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_NChi2->Fill( Vtx_NChi0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_nTrks->Fill( Vtx_nTrks0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_z->Fill( Vtx_z0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_r->Fill( Vtx_r0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_dR->Fill( Vtx_dR0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_dist->Fill( Vtx_dist0,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );
               hData_EVT12_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass1,Pref_PU_gen_wt);
      hData_EVT12_2VtxAll_BDTvtx->Fill( BDTvtx2,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_HMass->Fill( Vtx_HMass1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_NChi2->Fill( Vtx_NChi1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_nTrks->Fill( Vtx_nTrks1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_z->Fill( Vtx_z1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_r->Fill( Vtx_r1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_dR->Fill( Vtx_dR1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1,Pref_PU_gen_wt );
      hData_EVT12_2VtxAll_dist->Fill( Vtx_dist1,Pref_PU_gen_wt );
   }

   if (nVtx == 2 && !isCutEvt) {
      hData_NoEVT12_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      
      hData_NoEVT12_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      hData_NoEVT12_2Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_z->Fill( Vtx_z,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_r->Fill( Vtx_r,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
      hData_NoEVT12_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );
      

            hData_NoEVT12_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass0,Pref_PU_gen_wt);
      hData_NoEVT12_2VtxAll_BDTvtx->Fill( BDTvtx1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_HMass->Fill( Vtx_HMass0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_NChi2->Fill( Vtx_NChi0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_nTrks->Fill( Vtx_nTrks0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_z->Fill( Vtx_z0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_r->Fill( Vtx_r0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_dR->Fill( Vtx_dR0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_dist->Fill( Vtx_dist0,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );
               hData_NoEVT12_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass1,Pref_PU_gen_wt);
      hData_NoEVT12_2VtxAll_BDTvtx->Fill( BDTvtx2,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_HMass->Fill( Vtx_HMass1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_NChi2->Fill( Vtx_NChi1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_nTrks->Fill( Vtx_nTrks1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_z->Fill( Vtx_z1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_r->Fill( Vtx_r1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_dR->Fill( Vtx_dR1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1,Pref_PU_gen_wt );
      hData_NoEVT12_2VtxAll_CutEvt_dist->Fill( Vtx_dist1,Pref_PU_gen_wt );
   }

   if (nVtxLoose == 2 && isCutEvt) {
      hData_EVT34_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      
      hData_EVT34_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      hData_EVT34_2Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_z->Fill( Vtx_z,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_r->Fill( Vtx_r,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
      hData_EVT34_2Vtx_CutEvt_VtxVtxdist->Fill(  Vtx_Vtx_dist,Pref_PU_gen_wt );

            hData_EVT34_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass0,Pref_PU_gen_wt);
      hData_EVT34_2VtxAll_BDTvtx->Fill( BDTvtx1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_HMass->Fill( Vtx_HMass0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_NChi2->Fill( Vtx_NChi0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_nTrks->Fill( Vtx_nTrks0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_z->Fill( Vtx_z0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_r->Fill( Vtx_r0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_dR->Fill( Vtx_dR0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_dist->Fill( Vtx_dist0,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_VtxVtxdist->Fill(  Vtx_Vtx_dist,Pref_PU_gen_wt );
               hData_EVT34_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass1,Pref_PU_gen_wt);
      hData_EVT34_2VtxAll_BDTvtx->Fill( BDTvtx2,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_HMass->Fill( Vtx_HMass1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_NChi2->Fill( Vtx_NChi1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_nTrks->Fill( Vtx_nTrks1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_z->Fill( Vtx_z1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_r->Fill( Vtx_r1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_dR->Fill( Vtx_dR1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1,Pref_PU_gen_wt );
      hData_EVT34_2VtxAll_CutEvt_dist->Fill( Vtx_dist1,Pref_PU_gen_wt );
   }

   if (nVtxLoose == 2 && !isCutEvt) {
      hData_NoEVT34_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
     
      hData_NoEVT34_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_BDTvtx->Fill( BDTvtx ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_HMass->Fill( Vtx_HMass ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_NChi2->Fill( Vtx_NChi ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_nTrks->Fill( Vtx_nTrks ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_z->Fill( Vtx_z ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_r->Fill( Vtx_r ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_dR->Fill( Vtx_dR ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_dist->Fill( Vtx_dist ,Pref_PU_gen_wt);
      hData_NoEVT34_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist ,Pref_PU_gen_wt);

            hData_NoEVT34_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass0,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_BDTvtx->Fill( BDTvtx1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_HMass->Fill( Vtx_HMass0 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_NChi2->Fill( Vtx_NChi0 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_nTrks->Fill( Vtx_nTrks0 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_z->Fill( Vtx_z0 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_r->Fill( Vtx_r0 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_dR->Fill( Vtx_dR0 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight0 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight0 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_dist->Fill( Vtx_dist0 ,Pref_PU_gen_wt);
         hData_NoEVT34_2VtxAll_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass1,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_BDTvtx->Fill( BDTvtx2 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_HMass->Fill( Vtx_HMass1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_NChi2->Fill( Vtx_NChi1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_nTrks->Fill( Vtx_nTrks1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_z->Fill( Vtx_z1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_r->Fill( Vtx_r1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_dR->Fill( Vtx_dR1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_SumtrackWeight->Fill( Vtx_SumtrackWeight1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_MeantrackWeight->Fill( Vtx_MeantrackWeight1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_dist->Fill( Vtx_dist1 ,Pref_PU_gen_wt);
      hData_NoEVT34_2VtxAll_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist ,Pref_PU_gen_wt);
   }
    //----------------------//
   //$$

   //-----------------------------------------------------------//
   // ABCD using Evt and Vtx BDT
   //-----------------------------------------------------------//
   isCutEvt = false;
   isCutVtx = false;
   if (minitree_Evts_MVAval->at(0) > EVTSWP) isCutEvt = true;
   if (BDTvtx > VTXWP) isCutVtx = true;
   if (BDTvtx1 > VTXWP) isCutVtx1 = true;
   if (BDTvtx2 > VTXWP) isCutVtx2 = true;

   if (nVtx == 1 && isCutEvt)
      {
         hData_1Vtx_MVAval->Fill(BDTvtx,Pref_PU_gen_wt);
      }
   if (nVtx == 2 && isCutEvt)
      {
         hData_2Vtx_MVAval->Fill(BDTvtx,Pref_PU_gen_wt);
         hData_2VtxAll_MVAval->Fill(BDTvtx1,Pref_PU_gen_wt);
         hData_2VtxAll_MVAval->Fill(BDTvtx2,Pref_PU_gen_wt);
      }
   //----------------------//

   if (nVtx == 1 && isCutEvt  && isCutVtx) {
      hData_EVTVtx_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_EVTVtx_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
   }

   if (nVtx == 1 && !isCutEvt  && isCutVtx) {
      hData_NoEVTVtx_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_NoEVTVtx_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
   }

   if (nVtx == 1 && isCutEvt && !isCutVtx) {
      hData_EVTNoVtx_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_EVTNoVtx_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
   }

   if (nVtx == 1 && !isCutEvt && !isCutVtx) {
      hData_NoEVTNoVtx_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_NoEVTNoVtx_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
   }
   //----------------------//

   if (nVtx == 2 && isCutEvt && isCutVtx) {
      hData_EVTVtx_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_EVTVtx_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      if (isCutVtx1 && isCutVtx2)
         {
            hData_EVTVtx_2VtxAll_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
            hData_EVTVtx_2VtxAll_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
         }
   }

   if (nVtx == 2 && !isCutEvt && isCutVtx) {
      hData_NoEVTVtx_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_NoEVTVtx_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      if (isCutVtx1 && isCutVtx2)
         {
            hData_NoEVTVtx_2VtxAll_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
            hData_NoEVTVtx_2VtxAll_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
         }
   }

   if (nVtx == 2 && isCutEvt && !isCutVtx ) {
      hData_EVTNoVtx_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_EVTNoVtx_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      if (!isCutVtx1 && !isCutVtx2)
         {
            hData_EVTNoVtx_2VtxAll_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
            hData_EVTNoVtx_2VtxAll_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
         }
   }

   if (nVtx == 2 && !isCutEvt && !isCutVtx ) {
      hData_NoEVTNoVtx_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
      hData_NoEVTNoVtx_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
      if (!isCutVtx1 && !isCutVtx2)
         {
            hData_NoEVTNoVtx_2VtxAll_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
            hData_NoEVTNoVtx_2VtxAll_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
         }
   }

   //-----------------------------------------------------------//
   // ABCD using Hemipshere pt anf Tight+loose steps of vertexing 
   //-----------------------------------------------------------//
   //--------SR--------//
   
      if (!BlindSR) {
         
         if (Filter && (abs(minitree_Hemi_eta->at(0)) > 2.4 || abs(minitree_Hemi_eta->at(1)) > 2.4)) continue;
         hData_Hemi_BDTevt->Fill( minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);
         //$$
         if (VtxMass > 8.) isCutVtx = true;
         if (hemi_ptmin > 80.) isCutEvt = true;
         //$$
         
         if (nVtx == 0 && isCutEvt) {
            hData_Hemi_0Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
            hData_Hemi_0Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
         }

         
         if (nVtx == 1 && isCutEvt) {
            hData_Hemi_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
            hData_Hemi_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
            hData_Hemi_1Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);
            hData_Hemi_1Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

            hData_Hemi_1Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
            hData_Hemi_1Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
            hData_Hemi_1Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
            hData_Hemi_1Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
            hData_Hemi_1Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
            hData_Hemi_1Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
            hData_Hemi_1Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
            hData_Hemi_1Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
            hData_Hemi_1Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
            hData_Hemi_1Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
         }

         
         if (nVtx == 2 && isCutEvt) {
            hData_Hemi_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
            hData_Hemi_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
            hData_Hemi_2Vtx_CutEvt_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
            hData_Hemi_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

            hData_Hemi_2Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
            hData_Hemi_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );

            hData_Hemi_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass0,Pref_PU_gen_wt);
            hData_Hemi_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass1,Pref_PU_gen_wt);
            hData_Hemi_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1 ,Pref_PU_gen_wt);
            hData_Hemi_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2 ,Pref_PU_gen_wt);

                  hData_Hemi_2VtxAll_HMass->Fill( Vtx_HMass0,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_NChi2->Fill( Vtx_NChi0,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_nTrks->Fill( Vtx_nTrks0,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_z->Fill( Vtx_z0,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_r->Fill( Vtx_r0,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_dR->Fill( Vtx_dR0,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight0,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight0,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_dist->Fill( Vtx_dist0,Pref_PU_gen_wt );
                  hData_Hemi_2VtxAll_HMass->Fill( Vtx_HMass1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_NChi2->Fill( Vtx_NChi1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_nTrks->Fill( Vtx_nTrks1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_z->Fill( Vtx_z1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_r->Fill( Vtx_r1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_dR->Fill( Vtx_dR1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_dist->Fill( Vtx_dist1,Pref_PU_gen_wt );
            hData_Hemi_2VtxAll_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );


         }

      }

      //--------CR : Low Pt --------//

      isCutVtx = false; 
      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;

      if (Filter && hemi_ptmin > 40. && hemi_ptmin < 80.)
         {
             //$$
            if ( abs(minitree_Hemi_eta->at(0)) > 2.4 || abs(minitree_Hemi_eta->at(1)) > 2.4 ) continue;
            hData_CRlowpt_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
            if ( VtxMass > 8. ) isCutVtx = true; 
            // if ( hemi_ptmin > 80. ) isCutEvt = true;   
            //$$
            if (nVtx == 0 )
               {
                  hData_CRlowpt_0Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  hData_CRlowpt_0Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  
               }
            if (nVtx == 1 )
               {
                  hData_CRlowpt_1Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                  hData_CRlowpt_1Vtx_CutEvt_Mass->Fill(   VtxMass,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

                  hData_CRlowpt_1Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
                  hData_CRlowpt_1Vtx_dist->Fill(Vtx_dist,Pref_PU_gen_wt );
               }
            if (nVtx == 2 )
               {
                  hData_CRlowpt_2Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0) );
                  hData_CRlowpt_2Vtx_CutEvt_Mass->Fill(   VtxMass ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_CutEvt_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  hData_CRlowpt_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx ,Pref_PU_gen_wt);

                  hData_CRlowpt_2Vtx_HMass->Fill( Vtx_HMass ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_NChi2->Fill( Vtx_NChi ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_z->Fill( Vtx_z ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_r->Fill( Vtx_r ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_dR->Fill( Vtx_dR ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_dist->Fill( Vtx_dist ,Pref_PU_gen_wt);
                  hData_CRlowpt_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist ,Pref_PU_gen_wt);

                  hData_CRlowpt_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass0,Pref_PU_gen_wt );
                  hData_CRlowpt_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1,Pref_PU_gen_wt );
                  hData_CRlowpt_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2 ,Pref_PU_gen_wt);

                                    hData_CRlowpt_2VtxAll_HMass->Fill( Vtx_HMass0 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_NChi2->Fill( Vtx_NChi0 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_nTrks->Fill( Vtx_nTrks0,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_z->Fill( Vtx_z0 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_r->Fill( Vtx_r0 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_dR->Fill( Vtx_dR0 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight0 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight0 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_dist->Fill( Vtx_dist0 ,Pref_PU_gen_wt);
                                    hData_CRlowpt_2VtxAll_HMass->Fill( Vtx_HMass1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_NChi2->Fill( Vtx_NChi1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_nTrks->Fill( Vtx_nTrks1,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_z->Fill( Vtx_z1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_r->Fill( Vtx_r1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_dR->Fill( Vtx_dR1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_dist->Fill( Vtx_dist1 ,Pref_PU_gen_wt);
                  hData_CRlowpt_2VtxAll_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist ,Pref_PU_gen_wt);
               }

         }
      //--------CR : Loose  --------//

      isCutVtx = false;
      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if ( Filter && abs(minitree_Hemi_eta->at(0)) < 2.4 && abs(minitree_Hemi_eta->at(1)) < 2.4 ) 
         {
            hData_CRloose_BDTevt->Fill( minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);
            //$$
            if ( VtxMass > 8. ) isCutVtx = true; 
            if ( hemi_ptmin > 80. ) isCutEvt = true;
            //$$ 
            if (nVtxLoose == 0 && isCutEvt)
               {
                  hData_CRloose_0Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  hData_CRloose_0Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);
                  
               }
            if (nVtxLoose == 1 && isCutEvt)
               {
                  hData_CRloose_1Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                  hData_CRloose_1Vtx_CutEvt_Mass->Fill(   VtxMass,Pref_PU_gen_wt );

                  hData_CRloose_1Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

                  hData_CRloose_1Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
                  hData_CRloose_1Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );


               }
            if (nVtxLoose == 2 && isCutEvt)
               {
                  hData_CRloose_2Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                  hData_CRloose_2Vtx_CutEvt_Mass->Fill(   VtxMass,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_CutEvt_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

                  hData_CRloose_2Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
                  hData_CRloose_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );

                  hData_CRloose_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass0 ,Pref_PU_gen_wt);
                  hData_CRloose_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass1 ,Pref_PU_gen_wt);
                  hData_CRloose_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2,Pref_PU_gen_wt );

                                    hData_CRloose_2VtxAll_HMass->Fill( Vtx_HMass0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_NChi2->Fill( Vtx_NChi0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_nTrks->Fill( Vtx_nTrks0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_z->Fill( Vtx_z0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_r->Fill( Vtx_r0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_dR->Fill( Vtx_dR0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_dist->Fill( Vtx_dist0,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );

                                    hData_CRloose_2VtxAll_HMass->Fill( Vtx_HMass1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_NChi2->Fill( Vtx_NChi1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_nTrks->Fill( Vtx_nTrks1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_z->Fill( Vtx_z1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_r->Fill( Vtx_r1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_dR->Fill( Vtx_dR1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1,Pref_PU_gen_wt );
                  hData_CRloose_2VtxAll_dist->Fill( Vtx_dist1,Pref_PU_gen_wt );

               }  
         }

      //--------CR : Loose Lowpt  --------//
      isCutVtx = false;
      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if (Filter && hemi_ptmin > 40. && hemi_ptmin < 80. )
         {
            hData_CRlooselowpt_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
            //$$
            if ( VtxMass > 8. ) isCutVtx = true; 
            // if ( hemi_ptmin > 80. ) isCutEvt = true;   
            //$$
            if (nVtxLoose == 0 )
               {
                  hData_CRlooselowpt_0Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_BDTvtx->Fill( BDTvtx ,Pref_PU_gen_wt);
                  
               }
            if (nVtxLoose == 1 )
               {
                  hData_CRlooselowpt_1Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_CutEvt_Mass->Fill(   VtxMass,Pref_PU_gen_wt );
                  hData_CRlooselowpt_2Vtx_CutEvt_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  hData_CRlooselowpt_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

                  hData_CRlooselowpt_1Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_SumtrackWeight->Fill(Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_MeantrackWeight->Fill(Vtx_MeantrackWeight,Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
                  hData_CRlooselowpt_1Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );


               }
            if (nVtxLoose == 2)
               {
                  hData_CRlooselowpt_2Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  hData_CRlooselowpt_2Vtx_CutEvt_Mass->Fill(   VtxMass,Pref_PU_gen_wt );

                  hData_CRlooselowpt_2Vtx_HMass->Fill( Vtx_HMass ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_NChi2->Fill( Vtx_NChi ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_nTrks->Fill( Vtx_nTrks ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_z->Fill( Vtx_z ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_r->Fill( Vtx_r ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_dR->Fill( Vtx_dR ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_dist->Fill( Vtx_dist ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist ,Pref_PU_gen_wt);

                  hData_CRlooselowpt_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1,Pref_PU_gen_wt );
                  hData_CRlooselowpt_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2 ,Pref_PU_gen_wt);

                  hData_CRlooselowpt_2VtxAll_HMass->Fill( Vtx_HMass0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_NChi2->Fill( Vtx_NChi0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_nTrks->Fill( Vtx_nTrks0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_z->Fill( Vtx_z0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_r->Fill( Vtx_r0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_dR->Fill( Vtx_dR0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_dist->Fill( Vtx_dist0 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist ,Pref_PU_gen_wt);

                  hData_CRlooselowpt_2VtxAll_HMass->Fill( Vtx_HMass1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_NChi2->Fill( Vtx_NChi1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_nTrks->Fill( Vtx_nTrks1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_z->Fill( Vtx_z1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_r->Fill( Vtx_r1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_dR->Fill( Vtx_dR1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_SumtrackWeight->Fill( Vtx_SumtrackWeight1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_MeantrackWeight->Fill( Vtx_MeantrackWeight1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1 ,Pref_PU_gen_wt);
                  hData_CRlooselowpt_2VtxAll_dist->Fill( Vtx_dist1 ,Pref_PU_gen_wt);

               }

         }
      //--------CR : fwd  --------//
      if(Filter && ((abs(minitree_Hemi_eta->at(0)) > 2.4 && abs(minitree_Hemi_eta->at(0)) < 3.0) ||
         (abs(minitree_Hemi_eta->at(1)) > 2.4 && abs(minitree_Hemi_eta->at(1)) < 3.0)))
            {
               hData_CRfwd_BDTevt->Fill(minitree_Evts_MVAval->at(0),Pref_PU_gen_wt);
               if ( VtxMass > 8. ) isCutVtx = true; 
               if ( hemi_ptmin > 80. ) isCutEvt = true;   
               //$$
               if (nVtx == 0 && isCutEvt)
                  {
                     hData_CRfwd_0Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                     hData_CRfwd_0Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);

                  }
               if (nVtx == 1 && isCutEvt)
                  {
                     hData_CRfwd_1Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                     hData_CRfwd_1Vtx_CutEvt_Mass->Fill(   VtxMass,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);
                     hData_CRfwd_1Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

                     hData_CRfwd_1Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
                     hData_CRfwd_1Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
                  }
               if (nVtx == 2 && isCutEvt)
                  {
                     hData_CRfwd_2Vtx_CutEvt_Mmumu->Fill(  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_CutEvt_Mass->Fill(   VtxMass,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_CutEvt_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

                     hData_CRfwd_2Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );

                     hData_CRfwd_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass0 ,Pref_PU_gen_wt);
                     hData_CRfwd_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass1 ,Pref_PU_gen_wt);
                     hData_CRfwd_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1,Pref_PU_gen_wt );
                     hData_CRfwd_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2,Pref_PU_gen_wt );

                        hData_CRfwd_2Vtx_HMass->Fill( Vtx_HMass0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_NChi2->Fill( Vtx_NChi0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_nTrks->Fill( Vtx_nTrks0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_z->Fill( Vtx_z0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_r->Fill( Vtx_r0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_dR->Fill( Vtx_dR0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_dist->Fill( Vtx_dist0,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );

                        hData_CRfwd_2Vtx_HMass->Fill( Vtx_HMass1,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_NChi2->Fill( Vtx_NChi1,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_nTrks->Fill( Vtx_nTrks1,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_z->Fill( Vtx_z1,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_r->Fill( Vtx_r1,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_dR->Fill( Vtx_dR1,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight1,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight1,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1,Pref_PU_gen_wt );
                     hData_CRfwd_2Vtx_dist->Fill( Vtx_dist1,Pref_PU_gen_wt );

                  }
            }
         //--------CR : SameSign  --------//
      if (minitree_FilterSameSign->at(0) && !Filter &&
            abs(minitree_Hemi_eta->at(0)) < 2.4 && abs(minitree_Hemi_eta->at(1)) < 2.4)
         {
               hData_SameSign_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
               if ( VtxMass > 8. ) isCutVtx = true; 
               if ( hemi_ptmin > 80. ) isCutEvt = true;   
               //$$
               if (nVtx == 0 && isCutEvt) {
                  hData_CRss_0Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
                  hData_CRss_0Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  
               }

               if (nVtx == 1 && isCutEvt) {
                  hData_CRss_1Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
                  hData_CRss_1Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
                  hData_CRss_1Vtx_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  hData_CRss_1Vtx_BDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

                  hData_CRss_1Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
                  hData_CRss_1Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
                  hData_CRss_1Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
                  hData_CRss_1Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
                  hData_CRss_1Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
                  hData_CRss_1Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
                  hData_CRss_1Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  hData_CRss_1Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
                  hData_CRss_1Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
                  hData_CRss_1Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
               }

               
               if (nVtx == 2 && isCutEvt) {
                  hData_CRss_2Vtx_CutEvt_Mmumu->Fill(minitree_Mmumu->at(0),Pref_PU_gen_wt);
                  hData_CRss_2Vtx_CutEvt_Mass->Fill(VtxMass,Pref_PU_gen_wt);
                  hData_CRss_2Vtx_CutEvt_BDTevt->Fill( minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  hData_CRss_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx,Pref_PU_gen_wt );

                  hData_CRss_2Vtx_HMass->Fill( Vtx_HMass,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_NChi2->Fill( Vtx_NChi,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_nTrks->Fill( Vtx_nTrks,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_z->Fill( Vtx_z,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_r->Fill( Vtx_r,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_dR->Fill( Vtx_dR,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_dist->Fill( Vtx_dist,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );

                  hData_CRss_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass0,Pref_PU_gen_wt);
                  hData_CRss_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass1,Pref_PU_gen_wt);
                  hData_CRss_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1 ,Pref_PU_gen_wt);
                  hData_CRss_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2,Pref_PU_gen_wt );

                     hData_CRss_2Vtx_HMass->Fill( Vtx_HMass0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_NChi2->Fill( Vtx_NChi0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_nTrks->Fill( Vtx_nTrks0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_z->Fill( Vtx_z0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_r->Fill( Vtx_r0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_dR->Fill( Vtx_dR0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_dist->Fill( Vtx_dist0,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_CutEvt_VtxVtxdist->Fill( Vtx_Vtx_dist,Pref_PU_gen_wt );

                     hData_CRss_2Vtx_HMass->Fill( Vtx_HMass1,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_NChi2->Fill( Vtx_NChi1,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_nTrks->Fill( Vtx_nTrks1,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_z->Fill( Vtx_z1,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_r->Fill( Vtx_r1,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_dR->Fill( Vtx_dR1,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_SumtrackWeight->Fill( Vtx_SumtrackWeight1,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_MeantrackWeight->Fill( Vtx_MeantrackWeight1,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_track_MeanDCA_d->Fill( Vtx_track_MeanDCA_d1,Pref_PU_gen_wt );
                  hData_CRss_2Vtx_dist->Fill( Vtx_dist1,Pref_PU_gen_wt );

               }
         }

   } // end LOOP on events

   HistogramManager h ;
   h.WriteAllHistogramsInFile((Production+"/ABCD_"+thesample+".root").Data(),"recreate");

} // end ABCD::Loop
