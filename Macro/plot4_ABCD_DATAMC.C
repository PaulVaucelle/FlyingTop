#include <iostream>
#include <TROOT.h>
#include "TH1.h"

TCanvas * plot(int method, TString Prod, TString Name)
{
int stati=0;
bool fit= 1;
bool logy=1;

//$$
bool DATA = true;
//$$

// number of vertices:
//$$
  int nvtx = 1;
//$$
  float hmin = 0.5; // cannot be 0 for logy=1
//$$  float hmax = 1E5;	         // for eta<2.4 pt>80 or CRlowpt
  float hmax = 1E6;	         // for eta<2.4 pt>80 or CRlowpt
  float hmaxBD = 1E6;	         // for eta<2.4 pt>80 or CRlowpt
  if ( nvtx == 1 ) {
    hmax = 1E5;   // for eta<2.4 pt>80  
    hmaxBD = 1E7;
  }
  if ( nvtx == 2 ) {
    hmax = 1E3;   // for eta<2.4 pt>80  
    f = 1E6;
//     hmax = 1E6;   // for eta<2.4 pt>80  
//     hmaxBD = 1E9;
  }
  // TString Prod = "DATA_EMU_2017";
  //ABCD_EMu2018A 094
  // ABCD_2018A 095
  float ReScaleXS = 1.;
  TString Dmode = "DM_";//DM or SM

// Mumu
 TFile* f1_Data_emu  = new TFile("../DATA_MUMU_2018_20_06_2024/ABCD_"+Dmode+"OS_2p4_DM_Run2018A-UL2018_MiniAODv2-v1.root");
 TFile* f2_Data_emu  = new TFile("../DATA_MUMU_2018_20_06_2024/ABCD_"+Dmode+"OS_2p4_DM_Run2018B-UL2018_MiniAODv2-v1.root");
 TFile* f3_Data_emu  = new TFile("../DATA_MUMU_2018_20_06_2024/ABCD_"+Dmode+"OS_2p4_DM_Run2018C-UL2018_MiniAODv2-v1.root");
 TFile* f4_Data_emu  = new TFile("../DATA_MUMU_2018_20_06_2024/ABCD_"+Dmode+"OS_2p4_DM_Run2018D-UL2018_MiniAODv2-v1.root");

//Emu
//   TFile* f1_Data_emu  = new TFile("../"+Prod+"/ABCD_OS_2p4_DataA.root");
//  TFile* f2_Data_emu  = new TFile("../"+Prod+"/ABCD_OS_2p4_DataB.root");
//  TFile* f3_Data_emu  = new TFile("../"+Prod+"/ABCD_OS_2p4_DataC.root");
//  TFile* f4_Data_emu  = new TFile("../"+Prod+"/ABCD_OS_2p4_DataD.root");

//Mumu
TFile* f1_DY  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f2_DY  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f1_TT  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_TT  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_ST  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_ST  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_TTV = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8.root");
 TFile* f2_TTV = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f3_TTV = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8.root");
 TFile* f1_VV  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_VV  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 TFile* f3_VV  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 



//Emu
//  TFile* f1_DY  = new TFile("../"+Prod+"/ABCD_OS_2p4_DYJetsToLL_M10to50.root");
//  TFile* f2_DY  = new TFile("../"+Prod+"/ABCD_OS_2p4_DYJetsToLL_M_50_v1.root");
//  TFile* f1_TT  = new TFile("../"+Prod+"/ABCD_OS_2p4_TTTo2L2Nu.root");
//  TFile* f2_TT  = new TFile("../"+Prod+"/ABCD_OS_2p4_TTToSemiLeptonic.root");
//  TFile* f1_ST  = new TFile("../"+Prod+"/ABCD_OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays.root");
//  TFile* f2_ST  = new TFile("../"+Prod+"/ABCD_OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays.root");
//  TFile* f1_TTV = new TFile("../"+Prod+"/ABCD_OS_2p4_ttWJetsToLNu_5f_EWK.root");
//  TFile* f2_TTV = new TFile("../"+Prod+"/ABCD_OS_2p4_TTZToLL_5f.root");
//  TFile* f3_TTV = new TFile("../"+Prod+"/ABCD_OS_2p4_TTWW_TuneCP5.root");
//  TFile* f1_VV  = new TFile("../"+Prod+"/ABCD_OS_2p4_WWTo2L2Nu.root");
//  TFile* f2_VV  = new TFile("../"+Prod+"/ABCD_OS_2p4_WZTo2Q2L_mllmin4p0.root");
//  TFile* f3_VV  = new TFile("../"+Prod+"/ABCD_OS_2p4_ZZTo2Q2L_mllmin4p0.root");


//signal
 TFile* f1_LLP = new TFile("../PROD_CSI_10_06_2024/ABCD_"+Dmode+"OS_2p4_RPV_2018_smu200_neu180_ctau100.root");
 TFile* f2_LLP = new TFile("../PROD_CSI_10_06_2024/ABCD_"+Dmode+"OS_2p4_RPV_2018_smu300_neu180_ctau100.root");
 TFile* f3_LLP = new TFile("../PROD_CSI_10_06_2024/ABCD_"+Dmode+"OS_2p4_RPV_2018_smu400_neu300_ctau100.root");
 TFile* f4_LLP = new TFile("../PROD_CSI_10_06_2024/ABCD_"+Dmode+"OS_2p4_RPV_2018_smu500_neu350_ctau100.root");

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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx"; 
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
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
    TString HeaderA = "Tight + 40<Hpt<80";
    TString HeaderB = "Loose + 40<Hpt<80 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeantrackWeight";
  }
  if (method == 86)
  {
    TString htitleA = "hData_CRlowpt_TLVtx_SumtrackWeight";
    TString htitleB = "hData_CRlooselowpt_TLVtx_SumtrackWeight";
    TString htitleC = "hData_Hemi_TLVtx_SumtrackWeight";
    TString htitleD = "hData_CRloose_TLVtx_SumtrackWeight";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 25;
    TString HeaderA = "T/L + 40<Hpt<80";
    TString HeaderB = "T/L + 40<Hpt<80 ";
    TString HeaderC = "T/L + Hpt>80 ";
    TString HeaderD = "T/L + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeantrackWeight";
  }

      if (method == 87)
  {
    TString htitleA = "hData_CRlowpt_TLVtxAll_SumtrackWeight";
    TString htitleB = "hData_CRlooselowpt_TLVtxAll_SumtrackWeight";
    TString htitleC = "hData_Hemi_TLVtxAll_SumtrackWeight";
    TString htitleD = "hData_CRloose_TLVtxAll_SumtrackWeight";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 25;
    TString HeaderA = "T/L + 40<Hpt<80";
    TString HeaderB = "T/L 40<Hpt<80 ";
    TString HeaderC = "T/L Hpt>80 ";
    TString HeaderD = "T/L + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "MeantrackWeight";
  }


    if (method == 88)
  {
    TString htitleA = "hData_CRlowpt_TLVtx_Mass";
    TString htitleB = "hData_CRlooselowpt_TLVtx_Mass";
    TString htitleC = "hData_Hemi_TLVtx_Mass";
    TString htitleD = "hData_CRloose_TLVtx_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 100;
    TString HeaderA = "T/L + 40<Hpt<80";
    TString HeaderB = "T/L + 40<Hpt<80 ";
    TString HeaderC = "T/L + Hpt>80 ";
    TString HeaderD = "T/L + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Mass";
  }

      if (method == 89)
  {
    TString htitleA = "hData_CRlowpt_TLVtxAll_Mass";
    TString htitleB = "hData_CRlooselowpt_TLVtxAll_Mass";
    TString htitleC = "hData_Hemi_TLVtxAll_Mass";
    TString htitleD = "hData_CRloose_TLVtxAll_Mass";
    int nbin = 25; 
    float xmin = 0;
    float xmax = 100;
    TString HeaderA = "T/L + 40<Hpt<80";
    TString HeaderB = "T/L + 40<Hpt<80 ";
    TString HeaderC = "T/L + Hpt>80 ";
    TString HeaderD = "T/L + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Mass";
  }



    //-----------------------------------------------------------//
// xsec in pb
 float xsec_DY[2], xsec_TT[3], xsec_ST[2], xsec_TTV[3], xsec_VV[3]; 

//    xsec_DY[0] = ReScaleXS*18610; xsec_DY[1] = ReScaleXS*6077.0;
//  xsec_TT[0] = ReScaleXS*88.3; xsec_TT[1] = ReScaleXS*366.9; xsec_TT[2] = ReScaleXS*378.0;
//  xsec_ST[0] = ReScaleXS*10.89; xsec_ST[1] = ReScaleXS*10.87;
//  xsec_TTV[0] = ReScaleXS*0.29; xsec_TTV[1] = ReScaleXS*0.052; xsec_TTV[2] =  ReScaleXS*0.0070;
//  xsec_VV[0] = ReScaleXS*11.09; xsec_VV[1] = ReScaleXS*6.57; xsec_VV[2] = ReScaleXS*3.68;
//Mumu

 float RESCALE[12] = {73548220,75845390,125925000,100762000,10692740,11270430,13044000,37746080,268000,8437999,28577000,27448000};//Mumu
//  float RESCALE[12] = {1,101415800,146010000,1,11015960,12704300,29721000,38493120,944000,9994000,28577000,29357940};//EMu

//Mumu

 xsec_DY[0] = 18610; xsec_DY[1] = 6077.22;
 xsec_TT[0] = 88.3/72.; xsec_TT[1] = 366.9/300; xsec_TT[2] = 1.0;
 xsec_ST[0] = 10.89/32.29; xsec_ST[1] = 10.87/32.23;
 xsec_TTV[0] = 0.29/0.016; xsec_TTV[1] = 0.05/1; xsec_TTV[2] =  0.007/1;
 xsec_VV[0] = 12.178/11.07; xsec_VV[1] = 6.57/15.28; xsec_VV[2] = 3.68/8.5;

//   xsec_DY[0] = 1; xsec_DY[1] = 1;
//  xsec_TT[0] = 1.0; xsec_TT[1] = 1.0; xsec_TT[2] = 1.0;
//  xsec_ST[0] = 1.0; xsec_ST[1] = 1.0;
//  xsec_TTV[0] = 1.0; xsec_TTV[1] = 1.0; xsec_TTV[2] =  1.0;
//  xsec_VV[0] = 1.0; xsec_VV[1] = 1.0; xsec_VV[2] = 1.0;


// signal cross section (fb), choose 1 or 10
 int SigXsec = 10;
 if ( nvtx == 2 ) SigXsec = 0;
 float SigXsec1 = 6.25;//250
 float SigXsec2 = 3.06;//300
 float SigXsec3 = 0.94;//400
 float SigXsec4 = 0.34;//500

 float events = 0.;
 float lumi = 59700.;//59700 - 56856
 float norm = 0.;// EMU : 10¨1 MC/data and MUMU ~ 10

  TLegend* leg;
    
// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

TPad* pad1 = new TPad("pad1","This is pad1",0.04,0.65,0.48,0.95,21);
TPad* pad2 = new TPad("pad2","This is pad2",0.52,0.65,0.96,0.95,21);
TPad* pad3 = new TPad("pad3","This is pad3",0.04,0.15,0.48,0.45,21);
TPad* pad4 = new TPad("pad4","This is pad4",0.52,0.15,0.96,0.45,21);
// 

  TPad* pad8 = new TPad("pad8","This is pad8",0.04,0.5,0.48,0.65,21);
  TPad* pad9 = new TPad("pad9","This is pad9",0.52,0.5,0.96,0.65,21);
  TPad* pad10 = new TPad("pad10","This is pad10",0.04,0.05,0.48,0.15,21);
  TPad* pad11 = new TPad("pad11","This is pad11",0.52,0.05,0.96,0.15,21);

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

pad8->SetFillColor(0);
pad8->SetBorderMode(0);
pad8->SetFrameFillColor(10);
pad8->Draw();
pad8->SetLogy(0);
   pad8->SetTopMargin(0.07);
   pad8->SetBottomMargin(0.13);
   pad8->SetRightMargin(0.04);
   pad8->SetLeftMargin(0.16);

pad9->SetFillColor(0);
pad9->SetBorderMode(0);
pad9->SetFrameFillColor(10);
pad9->Draw();
pad9->SetLogy(0);
   pad9->SetTopMargin(0.07);
   pad9->SetBottomMargin(0.13);
   pad9->SetRightMargin(0.04);
   pad9->SetLeftMargin(0.16);

pad10->SetFillColor(0);
pad10->SetBorderMode(0);
pad10->SetFrameFillColor(10);
pad10->Draw();
pad10->SetLogy(0);
   pad10->SetTopMargin(0.07);
   pad10->SetBottomMargin(0.13);
   pad10->SetRightMargin(0.04);
   pad10->SetLeftMargin(0.16);

pad11->SetFillColor(0);
pad11->SetBorderMode(0);
pad11->SetFrameFillColor(10);
pad11->Draw();
pad11->SetLogy(0);
   pad11->SetTopMargin(0.07);
   pad11->SetBottomMargin(0.13);
   pad11->SetRightMargin(0.04);
   pad11->SetLeftMargin(0.16);

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
gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);

 TH1F* hsolve  = new TH1F("hsolve","",nbin,xmin,xmax);
 hsolve->Sumw2();

 TH1F* htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 TH1F* htotData  = new TH1F("htotData","",nbin,xmin,xmax);

TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
 hsolve->Sumw2();
  htotData->Sumw2();
  htotMC->Sumw2();

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

// *****************************************************************************

 pad1->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

 f2_Data_emu->cd();
 g2_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
 htotData->Add(g2_Data_emu, htotData, 1, 1);

 f3_Data_emu->cd();
 g3_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
 htotData->Add(g3_Data_emu, htotData, 1, 1);

 f4_Data_emu->cd();
 g4_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
 htotData->Add(g4_Data_emu, htotData, 1, 1);

//------x);

 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");

 e1_DY->Sumw2();

 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / RESCALE[0]; //e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleA);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 
 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
  std::cout<<"Inte: "<<e2_DY->Integral(0,3)<<std::endl;
  std::cout<<e2_DY->GetBinContent(2)<<std::endl;
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / RESCALE[1]; //e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleA);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);
 std::cout<<"DY2: "<<xsec_DY[1]<<"with norm "<<norm<<std::endl;

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / RESCALE[9]; //e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleA);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / RESCALE[10]; //e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleA);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / RESCALE[11]; //e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleA);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / RESCALE[6]; //e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleA);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi /RESCALE[7]; // e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleA);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / RESCALE[8]; //e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleA);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / RESCALE[4]; //e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleA);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / RESCALE[5]; //e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleA);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / RESCALE[2]; //e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleA);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, xsec_TT[0]*norm,0);
// std::cout<<"TT: "<<xsec_TT[0]*norm<<"with norm "<<norm<<std::endl;
// std::cout<<"tt itnegral "<<h_TT->Integral(0,80)<<std::endl;

 htotMC->Add(htotMC, h_DY, 1, 1);
 htotMC->Add(htotMC, h_VV, 1, 1);
 htotMC->Add(htotMC, h_TTV, 1, 1);
 htotMC->Add(htotMC, h_ST, 1, 1);
 htotMC->Add(htotMC, h_TT, 1, 1);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);

 htotMC->Draw("HE"); 
  htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(kGreen+1, 1);
 htotMC->SetLineColor(kGreen+1);
 htotMC->SetLineStyle(1);
 htotMC->SetLineWidth(1);
 htotMC->SetTickLength(0.03, "YZ");
 htotMC->SetTickLength(0.03,"X");
 htotMC->SetLabelOffset(0.015,"X");
 htotMC->SetLabelOffset(0.007,"Y");
 htotMC->SetLabelSize(0.045, "XYZ");
 htotMC->SetLabelFont(42, "XYZ"); 
 htotMC->SetTitleSize(0.045, "XYZ"); 
 htotMC->SetTitleFont(42, "XYZ");
 htotMC->SetTitleOffset(1.2,"X"); 
 htotMC->SetTitleOffset(1.3,"Y");
 htotMC->GetXaxis()->SetTitle(xtitle);
 htotMC->GetXaxis()->SetTitleColor(1);
 htotMC->GetYaxis()->SetTitle(ytitle);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
 htotMC->SetMinimum(hmin); 
 htotMC->SetMaximum(hmax); 
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);



 h_VV->Draw("HEsame"); 
 h_VV->SetFillColorAlpha(kOrange-2, 1);
 h_VV->SetLineColor(kOrange-2);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

 h_TTV->Draw("HEsame"); 
 h_TTV->SetFillColorAlpha(kAzure+4, 1);
 h_TTV->SetLineColor(kAzure+4);
 h_TTV->SetLineStyle(1);
 h_TTV->SetLineWidth(3);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColorAlpha(kAzure+2, 1);
 h_ST->SetLineColor(kAzure+2);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 h_TT->Draw("HEsame"); 
 h_TT->SetFillColorAlpha(kAzure+1, 1);
 h_TT->SetLineColor(kAzure+1);
 h_TT->SetLineStyle(1);
 h_TT->SetLineWidth(3);
 h_TT->SetTickLength(0.03, "YZ");
 h_TT->SetTickLength(0.03,"X");

htotData->Draw("PE1same");
htotData->SetMarkerStyle(20);
htotData->SetMarkerSize(1);
htotData->SetMarkerColor(kBlack);
htotData->SetLineColor(kBlack);
htotData->SetLineWidth(1);

 f1_LLP->cd();
 e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_LLP->Sumw2();
 norm = 0.;
 if ( e1_LLP->Integral(0,3) > 0. ) norm = SigXsec1*lumi / e1_LLP->Integral(0,3);
 g1_LLP = (TH1F*)gROOT->FindObject(htitleA);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleA);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleA);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
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

 hsolve->Add(hsolve, htotData, 0., 1.);

  leg = new TLegend(0.17,0.93,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("simulation 2018, 60 fb^{-1} (13 TeV)");
  leg->Draw();

  leg = new TLegend(0.64,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " e#mu data","PE1");
    leg->AddEntry(htotMC, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();


  leg = new TLegend(0.33,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.33,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderA);
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 40-80 GeV");
  leg->Draw();

  //Data/MC ---------------------
  pad8->cd();
  
  // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
  hDataMC->Divide(htotData,htotMC,1,1);
  hDataMC->Draw("PE1"); 
  hDataMC->SetFillStyle(1001);
  hDataMC->SetFillColorAlpha(kBlack, 1);
  hDataMC->SetLineColor(kBlack);
  hDataMC->SetLineStyle(1);
  hDataMC->SetLineWidth(1);
  hDataMC->SetTickLength(0.03, "YZ");
  hDataMC->SetTickLength(0.03,"X");
  hDataMC->SetLabelOffset(0.015,"X");
  hDataMC->SetLabelOffset(0.007,"Y");
  hDataMC->SetLabelSize(0.09, "XYZ");
  hDataMC->SetLabelFont(42, "XYZ"); 
  hDataMC->SetTitleSize(0.08, "XYZ"); 
  hDataMC->SetTitleFont(42, "XYZ");
  hDataMC->SetTitleOffset(0.8,"X"); 
  hDataMC->SetTitleOffset(0.8,"Y");
  hDataMC->GetXaxis()->SetTitle(xtitle);
  hDataMC->GetXaxis()->SetTitleColor(1);
  hDataMC->GetYaxis()->SetTitle("Data/MC");
  hDataMC->GetYaxis()->SetTitleColor(1);
  hDataMC->SetNdivisions(509,"XYZ");
  hDataMC->SetMinimum(0); 
  hDataMC->SetMaximum(2); 

// *****************************************************************************

 pad2->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(htitleC);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

 f2_Data_emu->cd();
 g2_Data_emu = (TH1F*)gROOT->FindObject(htitleC);
 htotData->Add(g2_Data_emu, htotData, 1, 1);

 f3_Data_emu->cd();
 g3_Data_emu = (TH1F*)gROOT->FindObject(htitleC);
 htotData->Add(g3_Data_emu, htotData, 1, 1);

 f4_Data_emu->cd();
 g4_Data_emu = (TH1F*)gROOT->FindObject(htitleC);
 htotData->Add(g4_Data_emu, htotData, 1, 1);



 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi /RESCALE[0]; // e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleC);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / RESCALE[1]; //e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleC);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi /RESCALE[9]; // e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleC);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi /RESCALE[10]; // e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleC);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi /RESCALE[11]; // e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleC);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi /RESCALE[6]; // e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleC);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi /RESCALE[7]; // e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleC);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi /RESCALE[8]; // e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleC);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / RESCALE[4]; //e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleC);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / RESCALE[5]; //e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleC);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / RESCALE[2]; //e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleC);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, xsec_TT[0]*norm,0);

 htotMC->Add(htotMC, h_DY, 1, 1);
 htotMC->Add(htotMC, h_VV, 1, 1);
 htotMC->Add(htotMC, h_TTV, 1, 1);
 htotMC->Add(htotMC, h_ST, 1, 1);
 htotMC->Add(htotMC, h_TT, 1, 1);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);


 htotMC->Draw("HE"); 
 htotMC->SetFillColorAlpha(kGreen+1, 1);
 htotMC->SetLineColor(kGreen+1);
 htotMC->SetLineStyle(1);
 htotMC->SetLineWidth(1);
 htotMC->SetTickLength(0.03, "YZ");
 htotMC->SetTickLength(0.03,"X");
 htotMC->SetLabelOffset(0.015,"X");
 htotMC->SetLabelOffset(0.007,"Y");
 htotMC->SetLabelSize(0.045, "XYZ");
 htotMC->SetLabelFont(42, "XYZ"); 
 htotMC->SetTitleSize(0.045, "XYZ"); 
 htotMC->SetTitleFont(42, "XYZ");
 htotMC->SetTitleOffset(1.2,"X"); 
 htotMC->SetTitleOffset(1.3,"Y");
 htotMC->GetXaxis()->SetTitle(xtitle);
 htotMC->GetXaxis()->SetTitleColor(1);
 htotMC->GetYaxis()->SetTitle(ytitle);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
 htotMC->SetMinimum(hmin); 
 htotMC->SetMaximum(hmax); 
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);

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

htotData->Draw("PE1same");
  htotData->SetFillColor(kBlack);
  htotData->SetLineColor(kBlack);
  htotData->SetLineStyle(1);
  htotData->SetLineWidth(1);

 f1_LLP->cd();
 e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_LLP->Sumw2();
 norm = 0.;
 if ( e1_LLP->Integral(0,3) > 0. ) norm = SigXsec1*lumi / e1_LLP->Integral(0,3);
 g1_LLP = (TH1F*)gROOT->FindObject(htitleC);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleC);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleC);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
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

//      hsolve->SetLineColor(kBlack);
// //  hsolve->SetFillStyle(3004);
//  hsolve->SetLineStyle(1);
//  hsolve->SetLineWidth(2);


  leg = new TLegend(0.17,0.93,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("simulation 2018, 60 fb^{-1} (13 TeV)");
  leg->Draw();

  leg = new TLegend(0.64,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " e#mu data","PE1");
    leg->AddEntry(htotMC, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
leg->AddEntry(hsolve," Prediction","FE4");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();


  leg = new TLegend(0.33,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.33,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderC);
  leg->Draw();

  pad9->cd();
  
  // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
  hDataMC->Divide(htotData,htotMC,1,1);
  hDataMC->Draw("PE1"); 
  hDataMC->SetFillStyle(1001);
  hDataMC->SetFillColorAlpha(kBlack, 1);
  hDataMC->SetLineColor(kBlack);
  hDataMC->SetLineStyle(1);
  hDataMC->SetLineWidth(1);
  hDataMC->SetTickLength(0.03, "YZ");
  hDataMC->SetTickLength(0.03,"X");
  hDataMC->SetLabelOffset(0.015,"X");
  hDataMC->SetLabelOffset(0.007,"Y");
  hDataMC->SetLabelSize(0.09, "XYZ");
  hDataMC->SetLabelFont(42, "XYZ"); 
  hDataMC->SetTitleSize(0.08, "XYZ"); 
  hDataMC->SetTitleFont(42, "XYZ");
  hDataMC->SetTitleOffset(0.8,"X"); 
  hDataMC->SetTitleOffset(0.8,"Y");
  hDataMC->GetXaxis()->SetTitle(xtitle);
  hDataMC->GetXaxis()->SetTitleColor(1);
  hDataMC->GetYaxis()->SetTitle(ytitle);
  hDataMC->GetYaxis()->SetTitleColor(1);
  hDataMC->SetNdivisions(509,"XYZ");
  hDataMC->SetMinimum(0); 
  hDataMC->SetMaximum(2);

 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//


  TCanvas *c2 = new TCanvas("c2", "plots",200,0,700,700);
  c2->SetFillColor(10);
  c2->SetFillStyle(4000);
  c2->SetBorderSize(2);
  TPad* pad5 = new TPad("pad5","This is pad5",0.04,0.55,0.75,0.96,21);
  TPad* pad6 = new TPad("pad6","This is pad6",0.04,0.35,0.75,0.55,21);
  TPad* pad7 = new TPad("pad7","This is pad7",0.04,0.1,0.75,0.35,21);



  pad5->SetFillColor(0);
pad5->SetBorderMode(0);
pad5->SetFrameFillColor(10);
pad5->Draw();
pad5->SetLogy(logy);
   pad5->SetTopMargin(0.07);
   pad5->SetBottomMargin(0.13);
   pad5->SetRightMargin(0.04);
   pad5->SetLeftMargin(0.16);

pad6->SetFillColor(0);
pad6->SetBorderMode(0);
pad6->SetFrameFillColor(10);
pad6->Draw();
pad6->SetLogy(0);
   pad6->SetTopMargin(0.07);
   pad6->SetBottomMargin(0.13);
   pad6->SetRightMargin(0.04);
   pad6->SetLeftMargin(0.16);

pad7->SetFillColor(0);
pad7->SetBorderMode(0);
pad7->SetFrameFillColor(10);
pad7->Draw();
pad7->SetLogy(0);
   pad7->SetTopMargin(0.07);
   pad7->SetBottomMargin(0.13);
   pad7->SetRightMargin(0.04);
   pad7->SetLeftMargin(0.16);


c2->cd();
  pad5->cd();
  htotData->Draw("E1"); 
 htotData->SetFillColor(kBlack);
 htotData->SetLineColor(kBlack);
 htotData->SetLineStyle(1);
 htotData->SetLineWidth(1);
 htotData->SetTickLength(0.03, "YZ");
 htotData->SetTickLength(0.03,"X");
 htotData->SetLabelOffset(0.015,"X");
 htotData->SetLabelOffset(0.007,"Y");
 htotData->SetLabelSize(0.045, "XYZ");
 htotData->SetLabelFont(42, "XYZ"); 
 htotData->SetTitleSize(0.045, "XYZ"); 
 htotData->SetTitleFont(42, "XYZ");
 htotData->SetTitleOffset(1.2,"X"); 
 htotData->SetTitleOffset(1.3,"Y");
 htotData->GetXaxis()->SetTitle(xtitle);
 htotData->GetXaxis()->SetTitleColor(1);
 htotData->GetYaxis()->SetTitle(ytitle);
 htotData->GetYaxis()->SetTitleColor(1);
 htotData->SetNdivisions(509,"XYZ");
 htotData->SetMinimum(hmin); 
 htotData->SetMaximum(hmax); 
 htotData->SetMarkerStyle(20);
 htotData->SetMarkerSize(1);

if ( !DATA ) { 
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
}

 f1_LLP->cd();
 e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_LLP->Sumw2();
 norm = 0.;
 if ( e1_LLP->Integral(0,3) > 0. ) norm = SigXsec1*lumi / e1_LLP->Integral(0,3);
 g1_LLP = (TH1F*)gROOT->FindObject(htitleC);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleC);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleC);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
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


  leg = new TLegend(0.17,0.93,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("simulation 2018, 60 fb^{-1} (13 TeV)");
  leg->Draw();

 
  leg = new TLegend(0.64,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " e#mu data","PE1");
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
    leg->AddEntry(htotMC, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  leg->AddEntry(hsolve," Prediction","LPE");
leg->AddEntry(hsolve," Prediction","FE4");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.33,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.33,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderC);
  leg->Draw();

pad6->cd();
//pulls between the predictions and data using the following formula:
//pull = (data - prediction) / sqrt(data+sigma^2_prediction)

// Poissionan errors are added due to low stats

 
TH1F* ratio  = new TH1F("ratio","",nbin,0,25);
ratio->Sumw2();
for (int i = 0; i< nbin; i++)
  {
    double data = htotData->GetBinContent(i);
    double pred = hsolve->GetBinContent(i);
    double sigma = hsolve->GetBinError(i);
    // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
    if (data == 0 || pred == 0) {ratio->AddBinContent(i,0);}    //when data == 0 and there is no prediction => you need this
    else
      {
        // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
        double pull = (data - pred) / sqrt(data + sigma*sigma);
        ratio->AddBinContent(i,pull);
      }
  }
 ratio->Draw("E4"); // PE1 ou E4
 ratio->SetFillColor(kRed);
 ratio->SetFillStyle(3004);
 ratio->SetLineColor(kBlack);
 ratio->SetLineStyle(1);
 ratio->SetLineWidth(1);
 ratio->SetTickLength(0.03, "YZ");
 ratio->SetTickLength(0.03,"X");
 ratio->SetLabelOffset(0.015,"X");
 ratio->SetLabelOffset(0.007,"Y");
 ratio->SetLabelSize(0.045, "XYZ");
 ratio->SetLabelFont(42, "XYZ"); 
 ratio->SetTitleSize(0.045, "XYZ"); 
 ratio->SetTitleFont(42, "XYZ");
 ratio->SetTitleOffset(1.2,"X"); 
 ratio->SetTitleOffset(1.3,"Y");
 ratio->GetXaxis()->SetTitle(xtitle);
 ratio->GetXaxis()->SetTitleColor(1);
 ratio->GetYaxis()->SetTitle("Pulls");
 ratio->GetYaxis()->SetTitleColor(1);
 ratio->SetNdivisions(509,"XYZ");
 ratio->SetMinimum(-5); 
 ratio->SetMaximum(5); 
 ratio->SetMarkerStyle(20);
 ratio->SetMarkerSize(1);

pad7->cd();
// // intégrale à droite
TH1F* Inte  = new TH1F("ratio","",nbin,0,25);
Inte->Sumw2();
for (int i = 0; i< nbin-1; i++)
  {
    double sumPredi = 0;
    double sumData = 0;
    double Interatio = 0;
    for (int j = i+1 ; j < nbin ; j++)
      {
          sumPredi += hsolve->GetBinContent(j);
          sumData += htotData->GetBinContent(j);  
      } 
      if (sumPredi == 0){Inte->AddBinContent(i,0);} 
      else {Inte->AddBinContent(i,sumData/sumPredi);}
 
  }
  Inte->Draw("PE1"); 
  //  hsolve->SetFillStyle(3004);
//  Inte->SetFillColor(kBlack);
 Inte->SetLineColor(kBlack);
 Inte->SetLineStyle(1);
 Inte->SetLineWidth(1);
 Inte->SetTickLength(0.03, "YZ");
 Inte->SetTickLength(0.03,"X");
 Inte->SetLabelOffset(0.015,"X");
 Inte->SetLabelOffset(0.007,"Y");
 Inte->SetLabelSize(0.045, "XYZ");
 Inte->SetLabelFont(42, "XYZ"); 
 Inte->SetTitleSize(0.045, "XYZ"); 
 Inte->SetTitleFont(42, "XYZ");
 Inte->SetTitleOffset(1.2,"X"); 
 Inte->SetTitleOffset(1.3,"Y");
 Inte->GetXaxis()->SetTitle(xtitle);
 Inte->GetXaxis()->SetTitleColor(1);
 Inte->GetYaxis()->SetTitle("Data/Prediction Integral ratio");
 Inte->GetYaxis()->SetTitleColor(1);
 Inte->SetNdivisions(509,"XYZ");
 Inte->SetMinimum(0); 
 Inte->SetMaximum(5); 
 Inte->SetMarkerStyle(20);
 Inte->SetMarkerSize(1);
 
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//


// *****************************************************************************
  c1->cd();
 pad3->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
 hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
 hmax = hmaxBD; 

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(htitleB);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

 f2_Data_emu->cd();
 g2_Data_emu = (TH1F*)gROOT->FindObject(htitleB);
 htotData->Add(g2_Data_emu, htotData, 1, 1);

 f3_Data_emu->cd();
 g3_Data_emu = (TH1F*)gROOT->FindObject(htitleB);
 htotData->Add(g3_Data_emu, htotData, 1, 1);

 f4_Data_emu->cd();
 g4_Data_emu = (TH1F*)gROOT->FindObject(htitleB);
 htotData->Add(g4_Data_emu, htotData, 1, 1);
 
 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / RESCALE[0]; // e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleB);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi /RESCALE[1]; // e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleB);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / RESCALE[9]; //e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleB);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / RESCALE[10]; //e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleB);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi /RESCALE[11]; // e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleB);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / RESCALE[6]; //e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleB);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / RESCALE[7]; //e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleB);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / RESCALE[8]; //e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleB);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi /RESCALE[4]; // e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleB);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / RESCALE[5]; //e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleB);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / RESCALE[2]; //e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleB);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, xsec_TT[0]*norm,0);

 htotMC->Add(htotMC, h_DY, 1, 1);
 htotMC->Add(htotMC, h_VV, 1, 1);
 htotMC->Add(htotMC, h_TTV, 1, 1);
 htotMC->Add(htotMC, h_ST, 1, 1);
 htotMC->Add(htotMC, h_TT, 1, 1);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);


 htotMC->Draw("HE"); 
 htotMC->SetFillColor(kGreen+1);
 htotMC->SetLineColor(kGreen+1);
 htotMC->SetLineStyle(1);
 htotMC->SetLineWidth(1);
 htotMC->SetTickLength(0.03, "YZ");
 htotMC->SetTickLength(0.03,"X");
 htotMC->SetLabelOffset(0.015,"X");
 htotMC->SetLabelOffset(0.007,"Y");
 htotMC->SetLabelSize(0.045, "XYZ");
 htotMC->SetLabelFont(42, "XYZ"); 
 htotMC->SetTitleSize(0.045, "XYZ"); 
 htotMC->SetTitleFont(42, "XYZ");
 htotMC->SetTitleOffset(1.2,"X"); 
 htotMC->SetTitleOffset(1.3,"Y");
 htotMC->GetXaxis()->SetTitle(xtitle);
 htotMC->GetXaxis()->SetTitleColor(1);
 htotMC->GetYaxis()->SetTitle(ytitle);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
 htotMC->SetMinimum(hmin); 
 htotMC->SetMaximum(hmax); 
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);


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


htotData->Draw("PE1same");
htotData->SetMarkerStyle(20);
htotData->SetMarkerSize(1);
htotData->SetMarkerColor(kBlack);
htotData->SetLineColor(kBlack);
htotData->SetLineWidth(1);


 f1_LLP->cd();
 e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_LLP->Sumw2();
 norm = 0.;
 if ( e1_LLP->Integral(0,3) > 0. ) norm = SigXsec1*lumi / e1_LLP->Integral(0,3);
 g1_LLP = (TH1F*)gROOT->FindObject(htitleB);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleB);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleB);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
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

 hsolve->Divide(hsolve, htotData, 1., 1.);

  leg = new TLegend(0.17,0.93,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("simulation 2018, 60 fb^{-1} (13 TeV)");
  leg->Draw();

  leg = new TLegend(0.64,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " e#mu data","PE1");
    leg->AddEntry(htotMC, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.33,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.33,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderB);
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 40-80 GeV");
  leg->Draw();

  pad10->cd();
  
  // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
  hDataMC->Divide(htotData,htotMC,1,1);
  hDataMC->Draw("PE1"); 
  hDataMC->SetFillStyle(1001);
  hDataMC->SetFillColorAlpha(kBlack, 1);
  hDataMC->SetLineColor(kBlack);
  hDataMC->SetLineStyle(1);
  hDataMC->SetLineWidth(1);
  hDataMC->SetTickLength(0.03, "YZ");
  hDataMC->SetTickLength(0.03,"X");
  hDataMC->SetLabelOffset(0.015,"X");
  hDataMC->SetLabelOffset(0.007,"Y");
  hDataMC->SetLabelSize(0.09, "XYZ");
  hDataMC->SetLabelFont(42, "XYZ"); 
  hDataMC->SetTitleSize(0.08, "XYZ"); 
  hDataMC->SetTitleFont(42, "XYZ");
  hDataMC->SetTitleOffset(0.8,"X"); 
  hDataMC->SetTitleOffset(0.8,"Y");
  hDataMC->GetXaxis()->SetTitle(xtitle);
  hDataMC->GetXaxis()->SetTitleColor(1);
  hDataMC->GetYaxis()->SetTitle(ytitle);
  hDataMC->GetYaxis()->SetTitleColor(1);
  hDataMC->SetNdivisions(509,"XYZ");
  hDataMC->SetMinimum(0); 
  hDataMC->SetMaximum(2);
// *****************************************************************************

 pad4->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(htitleD);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

 f2_Data_emu->cd();
 g2_Data_emu = (TH1F*)gROOT->FindObject(htitleD);
 htotData->Add(g2_Data_emu, htotData, 1, 1);

 f3_Data_emu->cd();
 g3_Data_emu = (TH1F*)gROOT->FindObject(htitleD);
 htotData->Add(g3_Data_emu, htotData, 1, 1);

 f4_Data_emu->cd();
 g4_Data_emu = (TH1F*)gROOT->FindObject(htitleD);
 htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / RESCALE[0]; //e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleD);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi /RESCALE[1]; // e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleD);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / RESCALE[9]; //e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleD);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi /RESCALE[10]; // e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleD);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / RESCALE[11]; //e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleD);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / RESCALE[6]; //e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleD);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi /RESCALE[7]; // e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleD);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi /RESCALE[8]; // e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleD);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / RESCALE[4]; //e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleD);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / RESCALE[5]; //e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleD);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / RESCALE[2]; //e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleD);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, xsec_TT[0]*norm,0);

 htotMC->Add(htotMC, h_DY, 1, 1);
 htotMC->Add(htotMC, h_VV, 1, 1);
 htotMC->Add(htotMC, h_TTV, 1, 1);
 htotMC->Add(htotMC, h_ST, 1, 1);
 htotMC->Add(htotMC, h_TT, 1, 1);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);


 htotMC->Draw("HE"); 
 htotMC->SetFillColor(kGreen+1);
 htotMC->SetLineColor(kGreen+1);
 htotMC->SetLineStyle(1);
 htotMC->SetLineWidth(1);
 htotMC->SetTickLength(0.03, "YZ");
 htotMC->SetTickLength(0.03,"X");
 htotMC->SetLabelOffset(0.015,"X");
 htotMC->SetLabelOffset(0.007,"Y");
 htotMC->SetLabelSize(0.045, "XYZ");
 htotMC->SetLabelFont(42, "XYZ"); 
 htotMC->SetTitleSize(0.045, "XYZ"); 
 htotMC->SetTitleFont(42, "XYZ");
 htotMC->SetTitleOffset(1.2,"X"); 
 htotMC->SetTitleOffset(1.3,"Y");
 htotMC->GetXaxis()->SetTitle(xtitle);
 htotMC->GetXaxis()->SetTitleColor(1);
 htotMC->GetYaxis()->SetTitle(ytitle);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
 htotMC->SetMinimum(hmin); 
 htotMC->SetMaximum(hmax); 
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);



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

htotData->Draw("PE1same");
htotData->SetMarkerStyle(20);
htotData->SetMarkerSize(1);
htotData->SetMarkerColor(kBlack);
htotData->SetLineColor(kBlack);
htotData->SetLineWidth(1);


 f1_LLP->cd();
 e1_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_LLP->Sumw2();
 norm = 0.;
 if ( e1_LLP->Integral(0,3) > 0. ) norm = SigXsec1*lumi / e1_LLP->Integral(0,3);
 g1_LLP = (TH1F*)gROOT->FindObject(htitleD);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleD);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleD);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
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

 hsolve->Multiply(hsolve, htotData, 1., 1.);

  leg = new TLegend(0.17,0.93,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("simulation 2018, 60 fb^{-1} (13 TeV)");
  leg->Draw();


  leg = new TLegend(0.64,0.5,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " e#mu data","PE1");
    leg->AddEntry(htotMC, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.33,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.33,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderD);
  leg->Draw();

  pad11->cd();
  
  // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
  hDataMC->Divide(htotData,htotMC,1,1);
  hDataMC->Draw("PE1"); 
  hDataMC->SetFillStyle(1001);
  hDataMC->SetFillColorAlpha(kBlack, 1);
  hDataMC->SetLineColor(kBlack);
  hDataMC->SetLineStyle(1);
  hDataMC->SetLineWidth(1);
  hDataMC->SetTickLength(0.03, "YZ");
  hDataMC->SetTickLength(0.03,"X");
  hDataMC->SetLabelOffset(0.015,"X");
  hDataMC->SetLabelOffset(0.007,"Y");
  hDataMC->SetLabelSize(0.09, "XYZ");
  hDataMC->SetLabelFont(42, "XYZ"); 
  hDataMC->SetTitleSize(0.08, "XYZ"); 
  hDataMC->SetTitleFont(42, "XYZ");
  hDataMC->SetTitleOffset(0.8,"X"); 
  hDataMC->SetTitleOffset(0.8,"Y");
  hDataMC->GetXaxis()->SetTitle(xtitle);
  hDataMC->GetXaxis()->SetTitleColor(1);
  hDataMC->GetYaxis()->SetTitle(ytitle);
  hDataMC->GetYaxis()->SetTitleColor(1);
  hDataMC->SetNdivisions(509,"XYZ");
  hDataMC->SetMinimum(0); 
  hDataMC->SetMaximum(2);
// *****************************************************************************

 pad2->cd();
//  hsolve->SetFillColor(kGray+1);
// const  TH1F*  hsolvetemp(hsolve);
// TGraphErrors TGsolve = new TGraphErrors(hsolvetemp);
//    TGsolve->SetFillColor(6);
//    TGsolve->SetFillStyle(3005);
//    TGsolve->Draw("a4same");
  
//  hsolve->SetLineColor(kBlack);
//  hsolve->SetFillColor(kBlack);
//  hsolve->SetFillStyle(3004);
// //  hsolve->SetLineStyle(1);
// //  hsolve->SetLineWidth(2);
//   hsolve->Draw("E4same");


// *****************************************************************************

  c1->Update();
  c2->cd();
  pad5->cd();

   hsolve->SetLineColor(kBlack);
 hsolve->SetFillColor(kBlack);
 hsolve->SetFillStyle(3004);
//  hsolve->SetLineStyle(1);
//  hsolve->SetLineWidth(2);
  hsolve->Draw("E4same");
  c2->SaveAs("../"+Prod+"/"+Name+"_"+Form("%d",method)+"_v2.pdf");

//   f1_Data_emu->Close();
//   f2_Data_emu->Close();
//   f3_Data_emu->Close();
//   f4_Data_emu->Close();

//  f1_DY->Close();
//   f2_DY->Close();
//   f1_TT->Close();
//   f2_TT->Close();
//   f1_ST->Close();
//   f2_ST->Close();
//   f1_TTV->Close();
//   f2_TTV->Close();
//   f3_TTV->Close();
//   f1_VV->Close();
//   f2_VV->Close();
//   f3_VV->Close();
 
//   f1_LLP->Close();
//   f2_LLP->Close();
//   f3_LLP->Close();
//   f4_LLP->Close();
  return c1;
}