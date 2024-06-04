#include <iostream>
#include <TROOT.h>
#include "TH1.h"

TCanvas * plot4_ABCD_solve(int method)
{
int stati=0;
bool fit= 1;
bool logy=1;

// number of vertices:
//$$
  int nvtx = 2;
//$$
  float hmin = 1.; // cannot be 0 for logy=1
//$$  float hmax = 1E5;	         // for eta<2.4 pt>80 or CRlowpt
  float hmax = 1E6;	         // for eta<2.4 pt>80 or CRlowpt
  float hmaxBD = 1E6;	         // for eta<2.4 pt>80 or CRlowpt
  if ( nvtx == 1 ) {
    hmax = 1E5;   // for eta<2.4 pt>80  
    hmaxBD = 1E7;
  }
  if ( nvtx == 2 ) {
//     hmax = 1E3;   // for eta<2.4 pt>80  
//     hmaxBD = 1E6;
    hmax = 1E6;   // for eta<2.4 pt>80  
    hmaxBD = 1E9;
  }
TString Prod = "23_04_2024";
 TFile* f1_DY  = new TFile("../"+Prod+"/ABCD_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f2_DY  = new TFile("../"+Prod+"/ABCD_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f1_TT  = new TFile("../"+Prod+"/ABCD_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_TT  = new TFile("../"+Prod+"/ABCD_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_ST  = new TFile("../"+Prod+"/ABCD_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_ST  = new TFile("../"+Prod+"/ABCD_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_TTV = new TFile("../"+Prod+"/ABCD_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8.root");
 TFile* f2_TTV = new TFile("../"+Prod+"/ABCD_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f3_TTV = new TFile("../"+Prod+"/ABCD_TTWW_TuneCP5_13TeV-madgraph-pythia8.root");
 TFile* f1_VV  = new TFile("../"+Prod+"/ABCD_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_VV  = new TFile("../"+Prod+"/ABCD_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 TFile* f3_VV  = new TFile("../"+Prod+"/ABCD_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 TFile* f1_LLP = new TFile("../"+Prod+"/ABCD_RPV_2018_smu250_neu200_ctau100.root");
 TFile* f2_LLP = new TFile("../"+Prod+"/ABCD_RPV_2018_smu300_neu250_ctau100.root");
 TFile* f3_LLP = new TFile("../"+Prod+"/ABCD_RPV_2018_smu400_neu300_ctau100.root");
 TFile* f4_LLP = new TFile("../"+Prod+"/ABCD_RPV_2018_smu500_neu350_ctau100.root");
 
//$$
//  TString htitleA = "hData_2Vtx_hemiPt_VtxQual_A_mVtx"; 
//  TString htitleB = "hData_2Vtx_hemiPt_VtxQual_B_mVtx"; 
//  TString htitleC = "hData_2Vtx_hemiPt_VtxQual_C_mVtx"; 
//  TString htitleD = "hData_2Vtx_hemiPt_VtxQual_D_mVtx"; 
//  TString xtitle = "max. vertex mass [GeV]";
//  if ( nvtx == 1 ) xtitle = "vertex mass [GeV]"; 
//  TString ytitle = "Events / 4 GeV"; 
//  int nbin = 25; 
//  float xmin = 0.;
//  float xmax = 100.;

//  TString htitleA = "hData_2Vtx_hemiPt_VtxQual_A_BDTevt"; 
//  TString htitleB = "hData_2Vtx_hemiPt_VtxQual_B_BDTevt"; 
//  TString htitleC = "hData_2Vtx_hemiPt_VtxQual_C_BDTevt"; 
//  TString htitleD = "hData_2Vtx_hemiPt_VtxQual_D_BDTevt"; 

//  if ( nvtx == 1 ) xtitle = "vertex BDT score"; 
 TString ytitle = "Events"; 
//  int nbin = 26; 
//  float xmin = -1.04;
//  float xmax =  1.04;
//$$
    TString htitleA = "hData_EVT34_1Vtx_BDTvtx";
    TString htitleB = "hData_NoEVT34_1Vtx_BDTvtx";
    TString htitleC = "hData_EVT12_1Vtx_BDTvtx";
    TString htitleD = "hData_NoEVT12_1Vtx_BDTvtx";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "A";
    TString HeaderB = "B";
    TString HeaderC = "C";
    TString HeaderD = " D";
    TString HeaderNVtx = "k Vtx";
    TString xtitle = "var";
  //-----------------------------------------------------------//
  // ABCD using Evt and Tight+looseWP 
  //-----------------------------------------------------------//
if (method == 0)
  {
    TString htitleA = "hData_EVT34_1Vtx_BDTvtx";
    TString htitleB = "hData_NoEVT34_1Vtx_BDTvtx";
    TString htitleC = "hData_EVT12_1Vtx_BDTvtx";
    TString htitleD = "hData_NoEVT12_1Vtx_BDTvtx";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Vtx BDT";

  }

if (method == 1)
  {
    TString htitleA = "hData_EVT34_2Vtx_BDTvtx";
    TString htitleB = "hData_NoEVT34_2Vtx_BDTvtx";
    TString htitleC = "hData_EVT12_2Vtx_BDTvtx";
    TString htitleD = "hData_NoEVT12_2Vtx_BDTvtx";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT<0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx BDT";
  }
//   //-----------------------------------------------------------//
if (method == 2)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_Mass";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_Mass";
    TString htitleC = "hData_EVT12_1Vtx_CutEvt_Mass";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT<0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Vtx Mass";
  }

if (method == 3)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_Mass";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_Mass";
    TString htitleC = "hData_EVT12_2Vtx_CutEvt_Mass";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT<0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx Mass";
  }
//   //-----------------------------------------------------------//
if (method == 4)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_Mmumu";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_Mmumu";
    TString htitleC = "hData_EVT12_1Vtx_CutEvt_Mmumu";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_Mmumu";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  500;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT<0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Mmumu";
  }
  if (method == 5)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_Mmumu";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_Mmumu";
    TString htitleC = "hData_EVT12_2Vtx_CutEvt_Mmumu";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_Mmumu";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  500;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT<0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx Mass";
  }
  
//   //-----------------------------------------------------------//
//   // ABCD using Evt and Vtx BDT
//   //-----------------------------------------------------------//
if (method == 6)
  {
    TString htitleA = "hData_EVTNoVtx_1Vtx_CutEvt_Mass";
    TString htitleB = "hData_NoEVTNoVtx_1Vtx_CutEvt_Mass";
    TString htitleC = "hData_EVTVtx_1Vtx_CutEvt_Mass";
    TString htitleD = "hData_NoEVTVtx_1Vtx_CutEvt_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "EVT + NoVtx";
    TString HeaderB = "NoEVT + NoVtx ";
    TString HeaderC = "EVT + Vtx ";
    TString HeaderD = "NoEVT + Vtx ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Vtx Mass";
  }
  if (method == 7)
  {
    TString htitleA = "hData_EVTNoVtx_2Vtx_CutEvt_Mass";
    TString htitleB = "hData_NoEVTNoVtx_2Vtx_CutEvt_Mass";
    TString htitleC = "hData_EVTVtx_2Vtx_CutEvt_Mass";
    TString htitleD = "hData_NoEVTVtx_2Vtx_CutEvt_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "EVT + NoVtx";
    TString HeaderB = "NoEVT + NoVtx ";
    TString HeaderC = "EVT + Vtx ";
    TString HeaderD = "NoEVT + Vtx ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx Mass";
  }
//    //-----------------------------------------------------------//

if (method == 8)
  {
    TString htitleA = "hData_EVTNoVtx_1Vtx_CutEvt_Mmumu";
    TString htitleB = "hData_NoEVTNoVtx_1Vtx_CutEvt_Mmumu";
    TString htitleC = "hData_EVTVtx_1Vtx_CutEvt_Mmumu";
    TString htitleD = "hData_NoEVTVtx_1Vtx_CutEvt_Mmumu";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  500;
    TString HeaderA = "EVT + NoVtx";
    TString HeaderB = "NoEVT + NoVtx ";
    TString HeaderC = "EVT + Vtx ";
    TString HeaderD = "NoEVT + Vtx ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Mmumu";
  }
  if (method == 9)
  {
    TString htitleA = "hData_EVTNoVtx_2Vtx_CutEvt_Mmumu";
    TString htitleB = "hData_NoEVTNoVtx_2Vtx_CutEvt_Mmumu";
    TString htitleC = "hData_EVTVtx_2Vtx_CutEvt_Mmumu";
    TString htitleD = "hData_NoEVTVtx_2Vtx_CutEvt_Mmumu";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  500;
    TString HeaderA = "EVT + NoVtx";
    TString HeaderB = "NoEVT + NoVtx ";
    TString HeaderC = "EVT + Vtx ";
    TString HeaderD = "NoEVT + Vtx ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Mmumu";
  }
//    //-----------------------------------------------------------//

if (method == 10)
  {
      TString htitleA = "hData_EVTNoVtx_2VtxAll_CutEvt_Mass";
    TString htitleB = "hData_NoEVTNoVtx_2VtxAll_CutEvt_Mass";
    TString htitleC = "hData_EVTVtx_2VtxAll_CutEvt_Mass";
    TString htitleD = "hData_NoEVTVtx_2VtxAll_CutEvt_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "EVT + NoVtx";
    TString HeaderB = "NoEVT + NoVtx ";
    TString HeaderC = "EVT + Vtx ";
    TString HeaderD = "NoEVT + Vtx ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Vtx Mass";
  }
  if (method == 11)
  {
    TString htitleA = "hData_EVTNoVtx_2VtxAll_CutEvt_Mmumu";
    TString htitleB = "hData_NoEVTNoVtx_2VtxAll_CutEvt_Mmumu";
    TString htitleC = "hData_EVTVtx_2VtxAll_CutEvt_Mmumu";
    TString htitleD = "hData_NoEVTVtx_2VtxAll_CutEvt_Mmumu";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  500;
    TString HeaderA = "EVT + NoVtx";
    TString HeaderB = "NoEVT + NoVtx ";
    TString HeaderC = "EVT + Vtx ";
    TString HeaderD = "NoEVT + Vtx ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Mmumu";
  }
//    //-----------------------------------------------------------//

//    //-----------------------------------------------------------//
//    // ABCD using Hemisphere pt and Tight+looseWP
//    //-----------------------------------------------------------//
if (method == 12)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_CutEvt_Mass";
    TString htitleB = "hData_CRlooselowpt_1Vtx_CutEvt_Mass";
    TString htitleC = "hData_Hemi_1Vtx_CutEvt_Mass";
    TString htitleD = "hData_CRloose_1Vtx_CutEvt_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Vtx Mass";
  }
  if (method == 13)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_CutEvt_Mass";
    TString htitleB = "hData_CRlooselowpt_2Vtx_CutEvt_Mass";
    TString htitleC = "hData_Hemi_2Vtx_CutEvt_Mass";
    TString htitleD = "hData_CRloose_2Vtx_CutEvt_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx Mass";
  }
// //-----------------------------------------------------------//

if (method == 14)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_CutEvt_Mmumu";
    TString htitleB = "hData_CRlooselowpt_1Vtx_CutEvt_Mmumu";
    TString htitleC = "hData_Hemi_1Vtx_CutEvt_Mmumu";
    TString htitleD = "hData_CRloose_1Vtx_CutEvt_Mmumu";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  500;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Mmumu";
  }

  if (method == 15)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_CutEvt_Mmumu";
    TString htitleB = "hData_CRlooselowpt_2Vtx_CutEvt_Mmumu";
    TString htitleC = "hData_Hemi_2Vtx_CutEvt_Mmumu";
    TString htitleD = "hData_CRloose_2Vtx_CutEvt_Mmumu";
    TString xtitle = "Mmumu";

    int nbin = 25; 
    float xmin = 0;
    float xmax =  500;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
  }
//   //-----------------------------------------------------------//

if (method == 16)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_BDTvtx";
    TString htitleB = "hData_CRlooselowpt_1Vtx_BDTvtx";
    TString htitleC = "hData_Hemi_1Vtx_BDTvtx";
    TString htitleD = "hData_CRloose_1Vtx_BDTvtx";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Vtx BDT";
  }

if (method == 17)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_CutEvt_MaxBDTvtx";
    TString htitleB = "hData_CRlooselowpt_2Vtx_CutEvt_MaxBDTvtx";
    TString htitleC = "hData_Hemi_2Vtx_CutEvt_MaxBDTvtx";
    TString htitleD = "hData_CRloose_2Vtx_CutEvt_MaxBDTvtx";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx BDT";
  }
//   //-----------------------------------------------------------//

if (method == 18)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_BDTevt";
    TString htitleB = "hData_CRlooselowpt_1Vtx_BDTevt";
    TString htitleC = "hData_Hemi_1Vtx_BDTevt";
    TString htitleD = "hData_CRloose_1Vtx_BDTevt";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "Evt BDT";
  }
  if (method == 19)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_CutEvt_BDTevt";
    TString htitleB = "hData_CRlooselowpt_2Vtx_CutEvt_BDTevt";
    TString htitleC = "hData_Hemi_2Vtx_CutEvt_BDTevt";
    TString htitleD = "hData_CRloose_2Vtx_CutEvt_BDTevt";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Evt BDT";
  }
//     //-----------------------------------------------------------//

if (method == 20)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_CutEvt_Mass";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_CutEvt_Mass";
    TString htitleC = "hData_Hemi_2VtxAll_CutEvt_Mass";
    TString htitleD = "hData_CRloose_2VtxAll_CutEvt_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx"; 
     TString xtitle = "Vtx Mass";
  }
  if (method == 21)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_CutEvt_BDTvtx";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_CutEvt_BDTvtx";
    TString htitleC = "hData_Hemi_2VtxAll_CutEvt_BDTvtx";
    TString htitleD = "hData_CRloose_2VtxAll_CutEvt_BDTvtx";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx BDT";
  }

if (method == 22)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_HMass";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_HMass";
    TString htitleC = "hData_EVT12_1Vtx_HMass";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_HMass";
    int nbin = 20; 
    float xmin = 0;
    float xmax =  1000;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "HMass";
  }
if (method == 23)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_NChi2";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_NChi2";
    TString htitleC = "hData_EVT12_1Vtx_NChi2";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_NChi2";
    int nbin = 10; 
    float xmin = 0;
    float xmax =  10;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "NChi2";
  }
if (method == 24)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_nTrks";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_nTrks";
    TString htitleC = "hData_EVT12_1Vtx_nTrks";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_nTrks";
    int nbin = 50; 
    float xmin = 0;
    float xmax =  50;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "nTrks";
  }
if (method == 25)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_r";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_r";
    TString htitleC = "hData_EVT12_1Vtx_r";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_r";
    int nbin = 35; 
    float xmin = 0;
    float xmax =  70;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "r";
  }
if (method == 26)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_z";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_z";
    TString htitleC = "hData_EVT12_1Vtx_z";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_z";
    int nbin = 50; 
    float xmin = -50;
    float xmax =  50;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "z";
  }
if (method == 27)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_dR";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_dR";
    TString htitleC = "hData_EVT12_1Vtx_dR";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_dR";
    int nbin = 50;
    float xmin = 0;
    float xmax =  10;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "dR";
  }
if (method == 28)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_SumtrackWeight";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_SumtrackWeight";
    TString htitleC = "hData_EVT12_1Vtx_SumtrackWeight";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_SumtrackWeight";
    int nbin = 25;
    float xmin = 0;
    float xmax =  25;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "SumtrackWeight";
  }
if (method == 29)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_track_MeanDCA_d";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_track_MeanDCA_d";
    TString htitleC = "hData_EVT12_1Vtx_track_MeanDCA_d";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_track_MeanDCA_d";
    int nbin = 25;
    float xmin = 0;
    float xmax =  25;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "MeanDCA";
  }

if (method == 30)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_dist";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_dist";
    TString htitleC = "hData_EVT12_1Vtx_dist";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_dist";
    int nbin = 20;
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "dist";
  }


if (method == 31)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_HMass";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_HMass";
    TString htitleC = "hData_EVT12_2Vtx_HMass";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_HMass";
    int nbin = 20; 
    float xmin = 0;
    float xmax =  1000;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "HMass";
  }
if (method == 32)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_NChi2";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_NChi2";
    TString htitleC = "hData_EVT12_2Vtx_NChi2";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_NChi2";
    int nbin = 10; 
    float xmin = 0;
    float xmax =  10;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "NChi2";
  }
if (method == 33)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_nTrks";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_nTrks";
    TString htitleC = "hData_EVT12_2Vtx_nTrks";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_nTrks";
    int nbin = 50; 
    float xmin = 0;
    float xmax =  50;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "nTrks";
  }
if (method == 34)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_r";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_r";
    TString htitleC = "hData_EVT12_2Vtx_r";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_r";
    int nbin = 35; 
    float xmin = 0;
    float xmax =  70;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "r";
  }
if (method == 35)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_z";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_z";
    TString htitleC = "hData_EVT12_2Vtx_z";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_z";
    int nbin = 50; 
    float xmin = -50;
    float xmax =  50;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "z";
  }
if (method == 36)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_dR";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_dR";
    TString htitleC = "hData_EVT12_2Vtx_dR";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_dR";
    int nbin = 50;
    float xmin = 0;
    float xmax =  10;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "dR";
  }
if (method == 37)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_SumtrackWeight";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_SumtrackWeight";
    TString htitleC = "hData_EVT12_2Vtx_SumtrackWeight";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_SumtrackWeight";
    int nbin = 25;
    float xmin = 0;
    float xmax =  25;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "SumtrackWeight";
  }
if (method == 38)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_track_MeanDCA_d";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_track_MeanDCA_d";
    TString htitleC = "hData_EVT12_2Vtx_track_MeanDCA_d";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_track_MeanDCA_d";
    int nbin = 25;
    float xmin = 0;
    float xmax =  25;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeanDCA";
  }

if (method == 39)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_dist";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_dist";
    TString htitleC = "hData_EVT12_2Vtx_dist";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_dist";
    int nbin = 20;
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "dist";
  }
if (method == 40)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_VtxVtxdist";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_VtxVtxdist";
    TString htitleC = "hData_EVT12_2Vtx_CutEvt_VtxVtxdist";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_VtxVtxdist";
    int nbin = 80;
    float xmin = 0;
    float xmax =  400;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx-Vtx dist";
  }

if (method == 41)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_HMass";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_HMass";
    TString htitleC = "hData_EVT12_2VtxAll_HMass";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_HMass";
    int nbin = 20; 
    float xmin = 0;
    float xmax =  1000;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "HMass";
  }
if (method == 42)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_NChi2";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_NChi2";
    TString htitleC = "hData_EVT12_2VtxAll_NChi2";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_NChi2";
    int nbin = 10; 
    float xmin = 0;
    float xmax =  10;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "NChi2";
  }
if (method == 43)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_nTrks";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_nTrks";
    TString htitleC = "hData_EVT12_2VtxAll_nTrks";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_nTrks";
    int nbin = 50; 
    float xmin = 0;
    float xmax =  50;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "nTrks";
  }
if (method == 44)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_r";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_r";
    TString htitleC = "hData_EVT12_2VtxAll_r";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_r";
    int nbin = 35; 
    float xmin = 0;
    float xmax =  70;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "r";
  }
if (method == 45)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_z";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_z";
    TString htitleC = "hData_EVT12_2VtxAll_z";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_z";
    int nbin = 50; 
    float xmin = -50;
    float xmax =  50;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "z";
  }
if (method == 46)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_dR";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_dR";
    TString htitleC = "hData_EVT12_2VtxAll_dR";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_dR";
    int nbin = 50;
    float xmin = 0;
    float xmax =  10;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "dR";
  }
if (method == 47)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_SumtrackWeight";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_SumtrackWeight";
    TString htitleC = "hData_EVT12_2VtxAll_SumtrackWeight";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_SumtrackWeight";
    int nbin = 25;
    float xmin = 0;
    float xmax =  25;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "SumtrackWeight";
  }
if (method == 48)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_track_MeanDCA_d";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_track_MeanDCA_d";
    TString htitleC = "hData_EVT12_2VtxAll_track_MeanDCA_d";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_track_MeanDCA_d";
    int nbin = 25;
    float xmin = 0;
    float xmax =  25;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeanDCA";
  }

if (method == 49)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_dist";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_dist";
    TString htitleC = "hData_EVT12_2VtxAll_dist";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_dist";
    int nbin = 20;
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "dist";
  }
if (method == 50)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_VtxVtxdist";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_VtxVtxdist";
    TString htitleC = "hData_EVT12_2VtxAll_CutEvt_VtxVtxdist";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_VtxVtxdist";
    int nbin = 80;
    float xmin = 0;
    float xmax =  400;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx-Vtx dist";
  }

if (method == 51)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_HMass";
    TString htitleB = "hData_CRlooselowpt_1Vtx_HMass";
    TString htitleC = "hData_Hemi_1Vtx_HMass";
    TString htitleD = "hData_CRloose_1Vtx_HMass";
    int nbin = 20; 
    float xmin = 0;
    float xmax =  1000;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "HMass";
  }
if (method == 52)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_NChi2";
    TString htitleB = "hData_CRlooselowpt_1Vtx_NChi2";
    TString htitleC = "hData_Hemi_1Vtx_NChi2";
    TString htitleD = "hData_CRloose_1Vtx_NChi2";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 10;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "NChi2";
  }
if (method == 53)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_nTrks";
    TString htitleB = "hData_CRlooselowpt_1Vtx_nTrks";
    TString htitleC = "hData_Hemi_1Vtx_nTrks";
    TString htitleD = "hData_CRloose_1Vtx_nTrks";
    int nbin = 50; 
    float xmin = 0;
    float xmax = 50;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "nTrks";
  }
if (method == 54)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_z";
    TString htitleB = "hData_CRlooselowpt_1Vtx_z";
    TString htitleC = "hData_Hemi_1Vtx_z";
    TString htitleD = "hData_CRloose_1Vtx_z";
    int nbin = 50; 
    float xmin = -50;
    float xmax =  50;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "z";
  }
if (method == 55)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_r";
    TString htitleB = "hData_CRlooselowpt_1Vtx_r";
    TString htitleC = "hData_Hemi_1Vtx_r";
    TString htitleD = "hData_CRloose_1Vtx_r";
    int nbin = 35; 
    float xmin = 0;
    float xmax = 70;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "r";
  }
if (method == 56)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_dR";
    TString htitleB = "hData_CRlooselowpt_1Vtx_dR";
    TString htitleC = "hData_Hemi_1Vtx_dR";
    TString htitleD = "hData_CRloose_1Vtx_dR";
    int nbin = 50; 
    float xmin = 0;
    float xmax = 10;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "dR";
  }
if (method == 57)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_SumtrackWeight";
    TString htitleB = "hData_CRlooselowpt_1Vtx_SumtrackWeight";
    TString htitleC = "hData_Hemi_1Vtx_SumtrackWeight";
    TString htitleD = "hData_CRloose_1Vtx_SumtrackWeight";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 25;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "SumtrackWeight";
  }
if (method == 58)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_track_MeanDCA_d";
    TString htitleB = "hData_CRlooselowpt_1Vtx_track_MeanDCA_d";
    TString htitleC = "hData_Hemi_1Vtx_track_MeanDCA_d";
    TString htitleD = "hData_CRloose_1Vtx_track_MeanDCA_d";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 25;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "MeanDCA";
  }
if (method == 59)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_dist";
    TString htitleB = "hData_CRlooselowpt_1Vtx_dist";
    TString htitleC = "hData_Hemi_1Vtx_dist";
    TString htitleD = "hData_CRloose_1Vtx_dist";
    int nbin = 20;
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "dist";
  }


if (method == 60)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_HMass";
    TString htitleB = "hData_CRlooselowpt_2Vtx_HMass";
    TString htitleC = "hData_Hemi_2Vtx_HMass";
    TString htitleD = "hData_CRloose_2Vtx_HMass";
    int nbin = 20; 
    float xmin = 0;
    float xmax =  1000;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "HMass";
  }
if (method == 61)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_NChi2";
    TString htitleB = "hData_CRlooselowpt_2Vtx_NChi2";
    TString htitleC = "hData_Hemi_2Vtx_NChi2";
    TString htitleD = "hData_CRloose_2Vtx_NChi2";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 10;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "NChi2";
  }
if (method == 62)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_nTrks";
    TString htitleB = "hData_CRlooselowpt_2Vtx_nTrks";
    TString htitleC = "hData_Hemi_2Vtx_nTrks";
    TString htitleD = "hData_CRloose_2Vtx_nTrks";
    int nbin = 50; 
    float xmin = 0;
    float xmax = 50;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "nTrks";
  }
if (method == 63)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_z";
    TString htitleB = "hData_CRlooselowpt_2Vtx_z";
    TString htitleC = "hData_Hemi_2Vtx_z";
    TString htitleD = "hData_CRloose_2Vtx_z";
    int nbin = 50; 
    float xmin = -50;
    float xmax =  50;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "z";
  }
if (method == 64)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_r";
    TString htitleB = "hData_CRlooselowpt_2Vtx_r";
    TString htitleC = "hData_Hemi_2Vtx_r";
    TString htitleD = "hData_CRloose_2Vtx_r";
    int nbin = 35; 
    float xmin = 0;
    float xmax = 70;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "r";
  }
if (method == 65)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_dR";
    TString htitleB = "hData_CRlooselowpt_2Vtx_dR";
    TString htitleC = "hData_Hemi_2Vtx_dR";
    TString htitleD = "hData_CRloose_2Vtx_dR";
    int nbin = 50; 
    float xmin = 0;
    float xmax = 10;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "dR";
  }
if (method == 66)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_SumtrackWeight";
    TString htitleB = "hData_CRlooselowpt_2Vtx_SumtrackWeight";
    TString htitleC = "hData_Hemi_2Vtx_SumtrackWeight";
    TString htitleD = "hData_CRloose_2Vtx_SumtrackWeight";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 25;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "SumtrackWeight";
  }
if (method == 67)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_track_MeanDCA_d";
    TString htitleB = "hData_CRlooselowpt_2Vtx_track_MeanDCA_d";
    TString htitleC = "hData_Hemi_2Vtx_track_MeanDCA_d";
    TString htitleD = "hData_CRloose_2Vtx_track_MeanDCA_d";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 25;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeanDCA";
  }
if (method == 68)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_dist";
    TString htitleB = "hData_CRlooselowpt_2Vtx_dist";
    TString htitleC = "hData_Hemi_2Vtx_dist";
    TString htitleD = "hData_CRloose_2Vtx_dist";
    int nbin = 20;
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "dist";
  }
if (method == 69)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_CutEvt_VtxVtxdist";
    TString htitleB = "hData_CRlooselowpt_2Vtx_CutEvt_VtxVtxdist";
    TString htitleC = "hData_Hemi_2Vtx_CutEvt_VtxVtxdist";
    TString htitleD = "hData_CRloose_2Vtx_CutEvt_VtxVtxdist";
    int nbin = 80; 
    float xmin = 0;
    float xmax = 400;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "VtxVtx dist";
  }

if (method == 70)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_HMass";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_HMass";
    TString htitleC = "hData_Hemi_2VtxAll_HMass";
    TString htitleD = "hData_CRloose_2VtxAll_HMass";
    int nbin = 20; 
    float xmin = 0;
    float xmax =  1000;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "HMass";
  }
if (method == 71)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_NChi2";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_NChi2";
    TString htitleC = "hData_Hemi_2VtxAll_NChi2";
    TString htitleD = "hData_CRloose_2VtxAll_NChi2";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 10;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "NChi2";
  }
if (method == 72)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_nTrks";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_nTrks";
    TString htitleC = "hData_Hemi_2VtxAll_nTrks";
    TString htitleD = "hData_CRloose_2VtxAll_nTrks";
    int nbin = 50; 
    float xmin = 0;
    float xmax = 50;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "nTrks";
  }
if (method == 73)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_z";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_z";
    TString htitleC = "hData_Hemi_2VtxAll_z";
    TString htitleD = "hData_CRloose_2VtxAll_z";
    int nbin = 50; 
    float xmin = -50;
    float xmax =  50;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "z";
  }
if (method == 74)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_r";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_r";
    TString htitleC = "hData_Hemi_2VtxAll_r";
    TString htitleD = "hData_CRloose_2VtxAll_r";
    int nbin = 35; 
    float xmin = 0;
    float xmax = 70;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "r";
  }
if (method == 75)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_dR";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_dR";
    TString htitleC = "hData_Hemi_2VtxAll_dR";
    TString htitleD = "hData_CRloose_2VtxAll_dR";
    int nbin = 50; 
    float xmin = 0;
    float xmax = 10;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "dR";
  }
if (method == 76)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_SumtrackWeight";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_SumtrackWeight";
    TString htitleC = "hData_Hemi_2VtxAll_SumtrackWeight";
    TString htitleD = "hData_CRloose_2VtxAll_SumtrackWeight";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 25;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "SumtrackWeight";
  }
if (method == 77)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_track_MeanDCA_d";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_track_MeanDCA_d";
    TString htitleC = "hData_Hemi_2VtxAll_track_MeanDCA_d";
    TString htitleD = "hData_CRloose_2VtxAll_track_MeanDCA_d";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 25;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeanDCA";
  }
if (method == 78)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_dist";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_dist";
    TString htitleC = "hData_Hemi_2VtxAll_dist";
    TString htitleD = "hData_CRloose_2VtxAll_dist";
    int nbin = 20;
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "dist";
  }
if (method == 79)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_CutEvt_VtxVtxdist";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_CutEvt_VtxVtxdist";
    TString htitleC = "hData_Hemi_2VtxAll_CutEvt_VtxVtxdist";
    TString htitleD = "hData_CRloose_2VtxAll_CutEvt_VtxVtxdist";
    int nbin = 80; 
    float xmin = 0;
    float xmax = 400;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "VtxVtx dist";
  }
//---------------//
  if (method == 80)
  {
    TString htitleA = "hData_EVT34_1Vtx_CutEvt_MeantrackWeight";
    TString htitleB = "hData_NoEVT34_1Vtx_CutEvt_MeantrackWeight";
    TString htitleC = "hData_EVT12_1Vtx_MeantrackWeight";
    TString htitleD = "hData_NoEVT12_1Vtx_CutEvt_MeantrackWeight";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 1;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "MeantrackWeight";
  }
  if (method == 81)
  {
    TString htitleA = "hData_EVT34_2Vtx_CutEvt_MeantrackWeight";
    TString htitleB = "hData_NoEVT34_2Vtx_CutEvt_MeantrackWeight";
    TString htitleC = "hData_EVT12_2Vtx_MeantrackWeight";
    TString htitleD = "hData_NoEVT12_2Vtx_CutEvt_MeantrackWeight";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 1;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeantrackWeight";
  }

    if (method == 82)
  {
    TString htitleA = "hData_EVT34_2VtxAll_CutEvt_MeantrackWeight";
    TString htitleB = "hData_NoEVT34_2VtxAll_CutEvt_MeantrackWeight";
    TString htitleC = "hData_EVT12_2VtxAll_MeantrackWeight";
    TString htitleD = "hData_NoEVT12_2VtxAll_CutEvt_MeantrackWeight";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 1;
    TString HeaderA = "Loose + EVTBDT>0";
    TString HeaderB = "Loose + EVTBDT <0 ";
    TString HeaderC = "Tight + EVTBDT>0. ";
    TString HeaderD = "Tight + EVTBDT<0. ";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeantrackWeight";
  }
    if (method == 83)
  {
    TString htitleA = "hData_CRlowpt_1Vtx_MeantrackWeight";
    TString htitleB = "hData_CRlooselowpt_1Vtx_MeantrackWeight";
    TString htitleC = "hData_Hemi_1Vtx_MeantrackWeight";
    TString htitleD = "hData_CRloose_1Vtx_MeantrackWeight";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 1;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "1 Vtx";
    TString xtitle = "MeantrackWeight";
  }

    if (method == 84)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_MeantrackWeight";
    TString htitleB = "hData_CRlooselowpt_2Vtx_MeantrackWeight";
    TString htitleC = "hData_Hemi_2Vtx_MeantrackWeight";
    TString htitleD = "hData_CRloose_2Vtx_MeantrackWeight";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 1;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeantrackWeight";
  }
    if (method == 85)
  {
    TString htitleA = "hData_CRlowpt_2VtxAll_MeantrackWeight";
    TString htitleB = "hData_CRlooselowpt_2VtxAll_MeantrackWeight";
    TString htitleC = "hData_Hemi_2VtxAll_MeantrackWeight";
    TString htitleD = "hData_CRloose_2VtxAll_MeantrackWeight";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 1;
    TString HeaderA = "Tight + Hpt<80 && Hpt>40";
    TString HeaderB = "Loose + Hpt<80 && Hpt>40 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeantrackWeight";
  }




    //-----------------------------------------------------------//

 float xsec_DY[2], xsec_TT[3], xsec_ST[2], xsec_TTV[3], xsec_VV[3]; 
 xsec_DY[0] = 5379.0; xsec_DY[1] = 15910.0;
 xsec_TT[0] = 88.3; xsec_TT[1] = 365.3; xsec_TT[2] = 378.0;
 xsec_ST[0] = 32.51; xsec_ST[1] = 32.45;
 xsec_TTV[0] = 0.29; xsec_TTV[1] = 0.052; xsec_TTV[2] =  0.0070;
 xsec_VV[0] = 11.09; xsec_VV[1] = 6.57; xsec_VV[2] = 3.68;
 
// signal cross section (fb), choose 1 or 10
 int SigXsec = 10;
 if ( nvtx == 2 ) SigXsec = 1;
 float SigXsec1 = 7.55;
 float SigXsec2 = 3.63;
 float SigXsec3 = 1.06;
 float SigXsec4 = 0.38;

 float events = 0.;
 float lumi = 60000.;
 float norm = lumi / 100000.;

  TLegend* leg;
    
// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",200,0,700,700);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

// TPad* pad1 = new TPad("pad1","This is pad1",0.04,0.49,0.48,0.91,21);
// TPad* pad2 = new TPad("pad2","This is pad2",0.52,0.49,0.96,0.91,21);
// TPad* pad3 = new TPad("pad3","This is pad3",0.04,0.04,0.48,0.46,21);
// TPad* pad4 = new TPad("pad4","This is pad4",0.52,0.04,0.96,0.46,21);
// 
TPad* pad1 = new TPad("pad1","This is pad1",0.03,0.48,0.49,0.92,21);
TPad* pad2 = new TPad("pad2","This is pad2",0.51,0.48,0.97,0.92,21);
TPad* pad3 = new TPad("pad3","This is pad3",0.03,0.03,0.49,0.47,21);
TPad* pad4 = new TPad("pad4","This is pad4",0.51,0.03,0.97,0.47,21);

pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.07);
   pad1->SetBottomMargin(0.13);
   pad1->SetRightMargin(0.04);
   pad1->SetLeftMargin(0.16);

pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(logy);
   pad2->SetTopMargin(0.07);
   pad2->SetBottomMargin(0.13);
   pad2->SetRightMargin(0.04);
   pad2->SetLeftMargin(0.16);

pad3->SetFillColor(0);
pad3->SetBorderMode(0);
pad3->SetFrameFillColor(10);
pad3->Draw();
pad3->SetLogy(logy);
   pad3->SetTopMargin(0.07);
   pad3->SetBottomMargin(0.13);
   pad3->SetRightMargin(0.04);
   pad3->SetLeftMargin(0.16);

pad4->SetFillColor(0);
pad4->SetBorderMode(0);
pad4->SetFrameFillColor(10);
pad4->Draw();
pad4->SetLogy(logy);
   pad4->SetTopMargin(0.07);
   pad4->SetBottomMargin(0.13);
   pad4->SetRightMargin(0.04);
   pad4->SetLeftMargin(0.16);

gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleFont(62);
gStyle->SetTitleColor(1);
gStyle->SetTitleTextColor(1);
gStyle->SetTitleFillColor(10);
gStyle->SetTitleFontSize(0.05);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(0.09);
gStyle->SetOptStat(stati);
gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
// gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);

 TH1F* hsolve  = new TH1F("hsolve","",nbin,xmin,xmax);
 hsolve->Sumw2();

 TH1F* htot  = new TH1F("htot","",nbin,xmin,xmax);
 TH1F* e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g1_DY = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g2_DY = (TH1F*)gROOT->FindObject(htitleA);
 TH1F*  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

 TH1F* e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g1_VV = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g2_VV = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g3_VV = (TH1F*)gROOT->FindObject(htitleA);
 TH1F*  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);

 TH1F* e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g1_TVV = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g2_TVV = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g3_TVV = (TH1F*)gROOT->FindObject(htitleA);
 TH1F*  h_TVV = new TH1F("h_TVV","",nbin,xmin,xmax);

 TH1F* e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g1_ST = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g2_ST = (TH1F*)gROOT->FindObject(htitleA);
 TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 TH1F* e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g1_TT = (TH1F*)gROOT->FindObject(htitleA);
 TH1F*  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 TH1F* e1_TTV = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g1_TTV = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* e2_TTV = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g2_TTV = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* e3_TTV = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g3_TTV = (TH1F*)gROOT->FindObject(htitleA);
 TH1F*  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);

 TH1F* e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g1_LLP = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 TH1F* e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g2_LLP = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 TH1F* e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 TH1F* g3_LLP = (TH1F*)gROOT->FindObject(htitleA);
 TH1F* h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);

// *****************************************************************************
  if (method <= 5) {c1->SetTitle("ABCD using Evt BDT and Tight+looseWP ");
      TLatex latex;
        latex.SetTextSize(0.04); // Set the text size
        latex.SetTextAlign(13); // Set text alignment (13 means left adjusted and top aligned)
        // Draw text on the canvas
        latex.DrawLatex(0.15, 0.97, "ABCD using Evt BDT and Tight+looseWP"); // Draw text at position (0.2, 0.8)
    }
  if (method > 5 && method <=11) {c1->SetTitle("ABCD using Evt and Vtx BDT ");
          TLatex latex;
          latex.SetTextSize(0.04); // Set the text size
          latex.SetTextAlign(13); // Set text alignment (13 means left adjusted and top aligned)
          // Draw text on the canvas
          latex.DrawLatex(0.15, 0.97, "ABCD using Evt and Vtx BDT"); // Draw text at position (0.2, 0.8)
        }
  if (method > 11 && method <=21) {c1->SetTitle("ABCD using Hemisphere pt and Tight+looseWP ");
            TLatex latex;
          latex.SetTextSize(0.04); // Set the text size
          latex.SetTextAlign(13); // Set text alignment (13 means left adjusted and top aligned)
          // Draw text on the canvas
          latex.DrawLatex(0.15, 0.97, "ABCD using Hemisphere pt and Tight+looseWP"); // Draw text at position (0.2, 0.8)
  }
  
 pad1->cd();

 htot = new TH1F("htot","",nbin,xmin,xmax);

 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = 60000. / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleA);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = 60000. / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleA);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = 60000. / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleA);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = 60000. / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleA);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = 60000. / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleA);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = 60000. / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleA);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = 60000. / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleA);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = 60000. / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleA);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = 60000. / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleA);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = 60000. / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleA);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = 60000. / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleA);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, xsec_TT[0]*norm,0);

 htot->Add(htot, h_DY, 1, 1);
 htot->Add(htot, h_VV, 1, 1);
 htot->Add(htot, h_TTV, 1, 1);
 htot->Add(htot, h_ST, 1, 1);
 htot->Add(htot, h_TT, 1, 1);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);

 htot->Draw("HE"); 
 htot->SetFillColor(kGreen-6);
 htot->SetLineColor(kGreen-6);
 htot->SetLineStyle(1);
 htot->SetLineWidth(3);
 htot->SetTickLength(0.03, "YZ");
 htot->SetTickLength(0.03,"X");
 htot->SetLabelOffset(0.015,"X");
 htot->SetLabelOffset(0.007,"Y");
 htot->SetLabelSize(0.045, "XYZ");
 htot->SetLabelFont(42, "XYZ"); 
 htot->SetTitleSize(0.045, "XYZ"); 
 htot->SetTitleFont(42, "XYZ");
 htot->SetTitleOffset(1.2,"X"); 
 htot->SetTitleOffset(1.3,"Y");
 htot->GetXaxis()->SetTitle(xtitle);
 htot->GetXaxis()->SetTitleColor(1);
 htot->GetYaxis()->SetTitle(ytitle);
 htot->GetYaxis()->SetTitleColor(1);
 htot->SetNdivisions(509,"XYZ");
 htot->SetMinimum(hmin); 
 htot->SetMaximum(hmax); 

 h_VV->Draw("HEsame"); 
 h_VV->SetFillColor(kOrange-2);
 h_VV->SetLineColor(kOrange-2);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

 h_TTV->Draw("HEsame"); 
 h_TTV->SetFillColor(kAzure+4);
 h_TTV->SetLineColor(kAzure+4);
 h_TTV->SetLineStyle(1);
 h_TTV->SetLineWidth(3);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColor(kAzure+2);
 h_ST->SetLineColor(kAzure+2);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 h_TT->Draw("HEsame"); 
 h_TT->SetFillColor(kAzure+1);
 h_TT->SetLineColor(kAzure+1);
 h_TT->SetLineStyle(1);
 h_TT->SetLineWidth(3);
 h_TT->SetTickLength(0.03, "YZ");
 h_TT->SetTickLength(0.03,"X");

 f1_LLP->cd();
 e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_LLP->Sumw2();
 norm = 0.;
 if ( e1_LLP->Integral(0,3) > 0. ) norm = SigXsec1*60000. / e1_LLP->Integral(0,3);
 g1_LLP = (TH1F*)gROOT->FindObject(htitleA);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*60000. / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleA);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*60000. / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleA);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*60000. / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleA);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 0.001*norm,0);

  
 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(kRed-1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);


 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(kRed-2);
 h2_LLP->SetLineStyle(2);
 h2_LLP->SetLineWidth(2);


 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(kRed-3);
 h3_LLP->SetLineStyle(3);
 h3_LLP->SetLineWidth(2);

//  h4_LLP->Draw("HEsame"); 
//  h4_LLP->SetLineColor(kRed-4);
//  h4_LLP->SetLineStyle(2);
//  h4_LLP->SetLineWidth(4);

//  e1_LLP->Draw("HEsame"); 
//  e1_LLP->SetLineColor(kWhite);
//  e1_LLP->SetLineStyle(1);
//  e1_LLP->SetLineWidth(3);

 hsolve->Add(hsolve, htot, 0., 1.);



  leg = new TLegend(0.17,0.93,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("simulation 2018, 60 fb^{-1} (13 TeV)");
  leg->Draw();

  leg = new TLegend(0.58,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htot, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  leg->AddEntry(h1_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 250/200 GeV","L");
  leg->AddEntry(h2_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 300/250 GeV","L");
  leg->AddEntry(h3_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 400/300 GeV","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.20,0.85,0.45,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  // if ( nvtx == 1 )
   leg->SetHeader(HeaderNVtx);
  // if ( nvtx == 2 ) leg->SetHeader(" 2 tight vertices");
  leg->Draw();
  leg = new TLegend(0.20,0.80,0.45,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  // if ( nvtx == 1 )
   leg->SetHeader(HeaderA);
  // if ( nvtx == 2 ) leg->SetHeader(" 1 hem. p_{T} 40-80 GeV");
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 40-80 GeV");
  leg->Draw();

// *****************************************************************************

 pad2->cd();

 htot = new TH1F("htot","",nbin,xmin,xmax);

 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = 60000. / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleC);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = 60000. / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleC);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = 60000. / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleC);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = 60000. / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleC);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = 60000. / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleC);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = 60000. / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleC);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = 60000. / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleC);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = 60000. / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleC);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = 60000. / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleC);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = 60000. / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleC);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = 60000. / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleC);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, xsec_TT[0]*norm,0);

 htot->Add(htot, h_DY, 1, 1);
 htot->Add(htot, h_VV, 1, 1);
 htot->Add(htot, h_TTV, 1, 1);
 htot->Add(htot, h_ST, 1, 1);
 htot->Add(htot, h_TT, 1, 1);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);

 htot->Draw("HE"); 
 htot->SetFillColor(kGreen-6);
 htot->SetLineColor(kGreen-6);
 htot->SetLineStyle(1);
 htot->SetLineWidth(3);
 htot->SetTickLength(0.03, "YZ");
 htot->SetTickLength(0.03,"X");
 htot->SetLabelOffset(0.015,"X");
 htot->SetLabelOffset(0.007,"Y");
 htot->SetLabelSize(0.045, "XYZ");
 htot->SetLabelFont(42, "XYZ"); 
 htot->SetTitleSize(0.045, "XYZ"); 
 htot->SetTitleFont(42, "XYZ");
 htot->SetTitleOffset(1.2,"X"); 
 htot->SetTitleOffset(1.3,"Y");
 htot->GetXaxis()->SetTitle(xtitle);
 htot->GetXaxis()->SetTitleColor(1);
 htot->GetYaxis()->SetTitle(ytitle);
 htot->GetYaxis()->SetTitleColor(1);
 htot->SetNdivisions(509,"XYZ");
 htot->SetMinimum(hmin); 
 htot->SetMaximum(hmax); 

 h_VV->Draw("HEsame"); 
 h_VV->SetFillColor(kOrange-2);
 h_VV->SetLineColor(kOrange-2);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

 h_TTV->Draw("HEsame"); 
 h_TTV->SetFillColor(kAzure+4);
 h_TTV->SetLineColor(kAzure+4);
 h_TTV->SetLineStyle(1);
 h_TTV->SetLineWidth(3);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColor(kAzure+2);
 h_ST->SetLineColor(kAzure+2);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 h_TT->Draw("HEsame"); 
 h_TT->SetFillColor(kAzure+1);
 h_TT->SetLineColor(kAzure+1);
 h_TT->SetLineStyle(1);
 h_TT->SetLineWidth(3);
 h_TT->SetTickLength(0.03, "YZ");
 h_TT->SetTickLength(0.03,"X");

 f1_LLP->cd();
 e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_LLP->Sumw2();
 norm = 0.;
 if ( e1_LLP->Integral(0,3) > 0. ) norm = SigXsec1*60000. / e1_LLP->Integral(0,3);

 g1_LLP = (TH1F*)gROOT->FindObject(htitleC);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);


 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*60000. / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleC);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*60000. / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleC);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*60000. / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleC);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 0.001*norm,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(kRed-1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);

 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(kRed-2);
 h2_LLP->SetLineStyle(2);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(kRed-3);
 h3_LLP->SetLineStyle(3);
 h3_LLP->SetLineWidth(2);

//  h4_LLP->Draw("HEsame"); 
//  h4_LLP->SetLineColor(kRed-4);
//  h4_LLP->SetLineStyle(2);
//  h4_LLP->SetLineWidth(4);

//  e1_LLP->Draw("HEsame"); 
//  e1_LLP->SetLineColor(kWhite);
//  e1_LLP->SetLineStyle(1);
//  e1_LLP->SetLineWidth(3);

  leg = new TLegend(0.17,0.93,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("simulation 2018, 60 fb^{-1} (13 TeV)");
  leg->Draw();

  leg = new TLegend(0.58,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htot, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  leg->AddEntry(h1_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 250/200 GeV","L");
  leg->AddEntry(h2_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 300/250 GeV","L");
  leg->AddEntry(h3_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 400/300 GeV","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.20,0.85,0.45,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  // if ( nvtx == 1 ) 
  leg->SetHeader(HeaderNVtx);
  // if ( nvtx == 2 ) leg->SetHeader(" 2 tight vertices");
  leg->Draw();
  leg = new TLegend(0.20,0.80,0.45,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderC);
  leg->Draw();

// *****************************************************************************

 pad3->cd();

 htot = new TH1F("htot","",nbin,xmin,xmax);
 hmax = hmaxBD; 

 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = 60000. / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleB);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = 60000. / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleB);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = 60000. / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleB);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = 60000. / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleB);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = 60000. / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleB);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = 60000. / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleB);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = 60000. / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleB);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = 60000. / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleB);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = 60000. / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleB);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = 60000. / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleB);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = 60000. / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleB);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, xsec_TT[0]*norm,0);

 htot->Add(htot, h_DY, 1, 1);
 htot->Add(htot, h_VV, 1, 1);
 htot->Add(htot, h_TTV, 1, 1);
 htot->Add(htot, h_ST, 1, 1);
 htot->Add(htot, h_TT, 1, 1);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);

 htot->Draw("HE"); 
 htot->SetFillColor(kGreen-6);
 htot->SetLineColor(kGreen-6);
 htot->SetLineStyle(1);
 htot->SetLineWidth(3);
 htot->SetTickLength(0.03, "YZ");
 htot->SetTickLength(0.03,"X");
 htot->SetLabelOffset(0.015,"X");
 htot->SetLabelOffset(0.007,"Y");
 htot->SetLabelSize(0.045, "XYZ");
 htot->SetLabelFont(42, "XYZ"); 
 htot->SetTitleSize(0.045, "XYZ"); 
 htot->SetTitleFont(42, "XYZ");
 htot->SetTitleOffset(1.2,"X"); 
 htot->SetTitleOffset(1.3,"Y");
 htot->GetXaxis()->SetTitle(xtitle);
 htot->GetXaxis()->SetTitleColor(1);
 htot->GetYaxis()->SetTitle(ytitle);
 htot->GetYaxis()->SetTitleColor(1);
 htot->SetNdivisions(509,"XYZ");
 htot->SetMinimum(hmin); 
 htot->SetMaximum(hmax); 

 h_VV->Draw("HEsame"); 
 h_VV->SetFillColor(kOrange-2);
 h_VV->SetLineColor(kOrange-2);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

 h_TTV->Draw("HEsame"); 
 h_TTV->SetFillColor(kAzure+4);
 h_TTV->SetLineColor(kAzure+4);
 h_TTV->SetLineStyle(1);
 h_TTV->SetLineWidth(3);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColor(kAzure+2);
 h_ST->SetLineColor(kAzure+2);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 h_TT->Draw("HEsame"); 
 h_TT->SetFillColor(kAzure+1);
 h_TT->SetLineColor(kAzure+1);
 h_TT->SetLineStyle(1);
 h_TT->SetLineWidth(3);
 h_TT->SetTickLength(0.03, "YZ");
 h_TT->SetTickLength(0.03,"X");

 f1_LLP->cd();
 e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_LLP->Sumw2();
 norm = 0.;
 if ( e1_LLP->Integral(0,3) > 0. ) norm = SigXsec1*60000. / e1_LLP->Integral(0,3);
 g1_LLP = (TH1F*)gROOT->FindObject(htitleB);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*60000. / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleB);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*60000. / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleB);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*60000. / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleB);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 0.001*norm,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(kRed-1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);

 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(kRed-2);
 h2_LLP->SetLineStyle(2);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(kRed-3);
 h3_LLP->SetLineStyle(3);
 h3_LLP->SetLineWidth(2);

//  h4_LLP->Draw("HEsame"); 
//  h4_LLP->SetLineColor(kRed-4);
//  h4_LLP->SetLineStyle(2);
//  h4_LLP->SetLineWidth(4);

//  e1_LLP->Draw("HEsame"); 
//  e1_LLP->SetLineColor(kWhite);
//  e1_LLP->SetLineStyle(1);
//  e1_LLP->SetLineWidth(3);

 hsolve->Divide(hsolve, htot, 1., 1.);

  leg = new TLegend(0.17,0.93,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("simulation 2018, 60 fb^{-1} (13 TeV)");
  leg->Draw();

  leg = new TLegend(0.58,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htot, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  leg->AddEntry(h1_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 250/200 GeV","L");
  leg->AddEntry(h2_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 300/250 GeV","L");
  leg->AddEntry(h3_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 400/300 GeV","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.20,0.85,0.45,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  // if ( nvtx == 1 ) 
  leg->SetHeader(HeaderNVtx);
  // if ( nvtx == 2 ) leg->SetHeader(" 2 loose vertices");
  leg->Draw();
  leg = new TLegend(0.20,0.80,0.45,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  // if ( nvtx == 1 ) 
  leg->SetHeader(HeaderB);
  // if ( nvtx == 2 ) leg->SetHeader(" 1 hem. p_{T} 40-80 GeV");
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 40-80 GeV");
  leg->Draw();

// *****************************************************************************

 pad4->cd();
 hmax = hmaxBD; 

 htot = new TH1F("htot","",nbin,xmin,xmax);
 
 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = 60000. / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleD);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = 60000. / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleD);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = 60000. / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleD);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = 60000. / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleD);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = 60000. / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleD);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = 60000. / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleD);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = 60000. / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleD);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = 60000. / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleD);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = 60000. / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleD);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = 60000. / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleD);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = 60000. / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleD);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, xsec_TT[0]*norm,0);

 htot->Add(htot, h_DY, 1, 1);
 htot->Add(htot, h_VV, 1, 1);
 htot->Add(htot, h_TTV, 1, 1);
 htot->Add(htot, h_ST, 1, 1);
 htot->Add(htot, h_TT, 1, 1);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);

 htot->Draw("HE"); 
 htot->SetFillColor(kGreen-6);
 htot->SetLineColor(kGreen-6);
 htot->SetLineStyle(1);
 htot->SetLineWidth(3);
 htot->SetTickLength(0.03, "YZ");
 htot->SetTickLength(0.03,"X");
 htot->SetLabelOffset(0.015,"X");
 htot->SetLabelOffset(0.007,"Y");
 htot->SetLabelSize(0.045, "XYZ");
 htot->SetLabelFont(42, "XYZ"); 
 htot->SetTitleSize(0.045, "XYZ"); 
 htot->SetTitleFont(42, "XYZ");
 htot->SetTitleOffset(1.2,"X"); 
 htot->SetTitleOffset(1.3,"Y");
 htot->GetXaxis()->SetTitle(xtitle);
 htot->GetXaxis()->SetTitleColor(1);
 htot->GetYaxis()->SetTitle(ytitle);
 htot->GetYaxis()->SetTitleColor(1);
 htot->SetNdivisions(509,"XYZ");
 htot->SetMinimum(hmin); 
 htot->SetMaximum(hmax); 

 h_VV->Draw("HEsame"); 
 h_VV->SetFillColor(kOrange-2);
 h_VV->SetLineColor(kOrange-2);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

 h_TTV->Draw("HEsame"); 
 h_TTV->SetFillColor(kAzure+4);
 h_TTV->SetLineColor(kAzure+4);
 h_TTV->SetLineStyle(1);
 h_TTV->SetLineWidth(3);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColor(kAzure+2);
 h_ST->SetLineColor(kAzure+2);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 h_TT->Draw("HEsame"); 
 h_TT->SetFillColor(kAzure+1);
 h_TT->SetLineColor(kAzure+1);
 h_TT->SetLineStyle(1);
 h_TT->SetLineWidth(3);
 h_TT->SetTickLength(0.03, "YZ");
 h_TT->SetTickLength(0.03,"X");

 f1_LLP->cd();
 e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_LLP->Sumw2();
 norm = 0.;
 if ( e1_LLP->Integral(0,3) > 0. ) norm = SigXsec1*60000. / e1_LLP->Integral(0,3);
 g1_LLP = (TH1F*)gROOT->FindObject(htitleD);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*60000. / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleD);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*60000. / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleD);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*60000. / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleD);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 0.001*norm,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(kRed-1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);

 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(kRed-2);
 h2_LLP->SetLineStyle(2);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(kRed-3);
 h3_LLP->SetLineStyle(3);
 h3_LLP->SetLineWidth(2);

//  h4_LLP->Draw("HEsame"); 
//  h4_LLP->SetLineColor(kRed-4);
//  h4_LLP->SetLineStyle(2);
//  h4_LLP->SetLineWidth(4);

//  e1_LLP->Draw("HEsame"); 
//  e1_LLP->SetLineColor(kWhite);
//  e1_LLP->SetLineStyle(1);
//  e1_LLP->SetLineWidth(3);

 hsolve->Multiply(hsolve, htot, 1., 1.);

  leg = new TLegend(0.17,0.93,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("simulation 2018, 60 fb^{-1} (13 TeV)");
  leg->Draw();

  leg = new TLegend(0.58,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htot, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  leg->AddEntry(h1_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 250/200 GeV","L");
  leg->AddEntry(h2_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 300/250 GeV","L");
  leg->AddEntry(h3_LLP," m_{#tilde{#mu}/#tilde{#chi}^{0}} = 400/300 GeV","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.20,0.85,0.45,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  // if ( nvtx == 1 ) 
  leg->SetHeader(HeaderNVtx);
  // if ( nvtx == 2 ) leg->SetHeader(" 2 loose vertices");
  leg->Draw();
  leg = new TLegend(0.20,0.80,0.45,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderD);
  leg->Draw();

// *****************************************************************************

 pad2->cd();
//  hsolve->SetFillColor(kGray+1);
 hsolve->Draw("sameE"); 
 hsolve->SetLineColor(kBlack);
//  hsolve->SetFillStyle(3004);
 hsolve->SetLineStyle(1);
 hsolve->SetLineWidth(2);

// *****************************************************************************

  c1->Update();
  // c1->SaveAs("test"+method+".pdf");
  return c1;
  // delete c1;
  
}
