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

//Mumu
//  TFile* f1_Data_emu  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DM_Run2018A-UL2018_MiniAODv2-v1.root");
//  TFile* f2_Data_emu  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DM_Run2018B-UL2018_MiniAODv2-v1.root");
//  TFile* f3_Data_emu  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DM_Run2018C-UL2018_MiniAODv2-v1.root");
//  TFile* f4_Data_emu  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DM_Run2018D-UL2018_MiniAODv2-v1.root");

//Emu
  TFile* f1_Data_emu  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DataA.root");
 TFile* f2_Data_emu  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DataB.root");
 TFile* f3_Data_emu  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DataC.root");
 TFile* f4_Data_emu  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DataD.root");

//Mumu
// TFile* f1_DY  = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
//  TFile* f2_DY  = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
//  TFile* f1_TT  = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f2_TT  = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f1_ST  = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f2_ST  = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f1_TTV = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8.root");
//  TFile* f2_TTV = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8.root");
//  TFile* f3_TTV = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8.root");
//  TFile* f1_VV  = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f2_VV  = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
//  TFile* f3_VV  = new TFile("../DATAMC2018_EMU_10_06_2024/ABCD_"+Dmode+"OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 
//Emu
 TFile* f1_DY  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DYJetsToLL_M10to50.root");
 TFile* f2_DY  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_DYJetsToLL_M_50_v1.root");
 TFile* f1_TT  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_TTTo2L2Nu.root");
 TFile* f2_TT  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_TTToSemiLeptonic.root");
 TFile* f1_ST  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays.root");
 TFile* f2_ST  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays.root");
 TFile* f1_TTV = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_ttWJetsToLNu_5f_EWK.root");
 TFile* f2_TTV = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_TTZToLL_5f.root");
 TFile* f3_TTV = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_TTWW_TuneCP5.root");
 TFile* f1_VV  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_WWTo2L2Nu.root");
 TFile* f2_VV  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_WZTo2Q2L_mllmin4p0.root");
 TFile* f3_VV  = new TFile("../"+Prod+"/ABCD_"+Dmode+"OS_2p4_ZZTo2Q2L_mllmin4p0.root");


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

        TString htitleE = "hData_EVT34_1Vtx_BDTvtx";
    TString htitleF = "hData_NoEVT34_1Vtx_BDTvtx";
    TString htitleG = "hData_EVT12_1Vtx_BDTvtx";
    TString htitleH = "hData_NoEVT12_1Vtx_BDTvtx";
    TString htitleI = "hData_NoEVT12_1Vtx_BDTvtx";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "A";
    TString HeaderB = "B";
    TString HeaderC = "C";
    TString HeaderD = " D";
    TString HeaderE = " E";
    TString HeaderF = " F";
    TString HeaderG = " G";
    TString HeaderH = " H";
    TString HeaderI = " I";
    TString HeaderNVtx = "k Vtx";
    TString xtitle = "var";

//    //-----------------------------------------------------------//
//    // ABCD using Hemisphere pt and Tight+looseWP
//    //-----------------------------------------------------------//

  if (method == 13)
  {
    TString htitleA = "hData_CRlowlowpt_2Vtx_CutEvt_Mass";//ok
    TString htitleB = "hData_CRlowlowpt_TLVtx_CutEvt_Mass";//ok
    TString htitleC = "hData_CRlooselowlowpt_TLVtx_CutEvt_Mass";//ok

    TString htitleD = "hData_CRlowpt_2Vtx_CutEvt_Mass";//ok
    TString htitleE = "hData_CRlowpt_TLVtx_CutEvt_Mass";//ok
    TString htitleF = "hData_CRlooselowpt_2Vtx_CutEvt_Mass";//ok

    TString htitleG = "hData_Hemi_2Vtx_CutEvt_Mass";//ok
    TString htitleH = "hData_Hemi_TLVtx_CutEvt_Mass";//ok
    TString htitleI = "hData_CRloose_2Vtx_CutEvt_Mass";//ok

    int nbin = 25; 
    float xmin = 0;
    float xmax =  100;
    TString HeaderA = "TT + 20<Hpt<80";
    TString HeaderB = "TL + 20<Hpt<80 ";
    TString HeaderC = "LL + Hpt>80 ";

    TString HeaderD = "TT + Hpt1>80 & 20<Hpt2<80";
    TString HeaderE = "TL + Hpt1>80 & 20<Hpt2<80";
    TString HeaderF = "LL + Hpt1>80 & 20<Hpt2<80 ";
    
    TString HeaderG = "TT + Hpt>80 ";
    TString HeaderH = "TL + Hpt>80";
    TString HeaderI = "LL + Hpt>80";

    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx Mass";
  }
// //-----------------------------------------------------------//



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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
  }
//   //-----------------------------------------------------------//


if (method == 17)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_CutEvt_MaxBDTvtx";
    TString htitleB = "hData_CRlooselowpt_2Vtx_CutEvt_MaxBDTvtx";
    TString htitleC = "hData_Hemi_2Vtx_CutEvt_MaxBDTvtx";
    TString htitleD = "hData_CRloose_2Vtx_CutEvt_MaxBDTvtx";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx BDT";
  }
//   //-----------------------------------------------------------//

  if (method == 19)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_CutEvt_BDTevt";
    TString htitleB = "hData_CRlooselowpt_2Vtx_CutEvt_BDTevt";
    TString htitleC = "hData_Hemi_2Vtx_CutEvt_BDTevt";
    TString htitleD = "hData_CRloose_2Vtx_CutEvt_BDTevt";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Vtx BDT";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
    TString HeaderC = "Tight + Hpt>80 ";
    TString HeaderD = "Loose + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "VtxVtx dist";
  }
//---------------//

    if (method == 84)
  {
    TString htitleA = "hData_CRlowpt_2Vtx_MeantrackWeight";
    TString htitleB = "hData_CRlooselowpt_2Vtx_MeantrackWeight";
    TString htitleC = "hData_Hemi_2Vtx_MeantrackWeight";
    TString htitleD = "hData_CRloose_2Vtx_MeantrackWeight";
    int nbin = 10; 
    float xmin = 0;
    float xmax = 1;
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "Tight + 20<Hpt<80";
    TString HeaderB = "Loose + 20<Hpt<80 ";
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
    TString HeaderA = "T/L + 20<Hpt<80";
    TString HeaderB = "T/L + 20<Hpt<80 ";
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
    TString HeaderA = "T/L + 20<Hpt<80";
    TString HeaderB = "T/L 20<Hpt<80 ";
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
    TString HeaderA = "T/L + 20<Hpt<80";
    TString HeaderB = "T/L + 20<Hpt<80 ";
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
    TString HeaderA = "T/L + 20<Hpt<80";
    TString HeaderB = "T/L + 20<Hpt<80 ";
    TString HeaderC = "T/L + Hpt>80 ";
    TString HeaderD = "T/L + Hpt>80";
    TString HeaderNVtx = "2 Vtx";
    TString xtitle = "Mass";
  }



    //-----------------------------------------------------------//
// xsec in pb
 float xsec_DY[2], xsec_TT[3], xsec_ST[2], xsec_TTV[3], xsec_VV[3]; 
 xsec_DY[0] = ReScaleXS*5379.0; xsec_DY[1] = ReScaleXS*15910.0;
 xsec_TT[0] = ReScaleXS*88.3; xsec_TT[1] = ReScaleXS*365.3; xsec_TT[2] = ReScaleXS*378.0;
 xsec_ST[0] = ReScaleXS*32.51; xsec_ST[1] = ReScaleXS*32.45;
 xsec_TTV[0] = ReScaleXS*0.29; xsec_TTV[1] = ReScaleXS*0.052; xsec_TTV[2] =  ReScaleXS*0.0070;
 xsec_VV[0] = ReScaleXS*11.09; xsec_VV[1] = ReScaleXS*6.57; xsec_VV[2] = ReScaleXS*3.68;
 
// signal cross section (fb), choose 1 or 10
 int SigXsec = 10;
 if ( nvtx == 2 ) SigXsec = 0;
 float SigXsec1 = 6.25;//250
 float SigXsec2 = 3.06;//300
 float SigXsec3 = 0.94;//400
 float SigXsec4 = 0.34;//500

 float events = 0.;
 float lumi = 60000.;// pb-1  it's supposed to be 60000 but there is an issue with that... that is not solved and not understood (need to find where that 10^3 comes from)
 float norm = 0.;

  TLegend* leg;
    
// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

TPad* pA = new TPad("pA","This is pad1",0.01,0.66,0.33,0.99,21);
TPad* pB = new TPad("pB","This is pad2",0.01,0.34,0.33,0.65,21);
TPad* pC = new TPad("pC","This is pad3",0.01,0.01,0.33,0.33,21);

TPad* pD = new TPad("pD","This is pad4",0.34,0.66,0.66,0.99,21);
TPad* pE = new TPad("pE","This is pad4",0.34,0.34,0.66,0.65,21);
TPad* pF = new TPad("pF","This is pad1",0.34,0.01,0.66,0.33,21);

TPad* pG = new TPad("pG","This is pad2",0.67,0.66,0.99,0.99,21);
TPad* pH = new TPad("pH","This is pad3",0.67,0.34,0.99,0.65,21);
TPad* pI = new TPad("pI","This is pad4",0.67,0.01,0.99,0.33,21);

// 

  // TPad* pad8 = new TPad("pad8","This is pad8",0.04,0.5,0.48,0.65,21);
  // TPad* pad9 = new TPad("pad9","This is pad9",0.52,0.5,0.96,0.65,21);
  // TPad* pad10 = new TPad("pad10","This is pad10",0.04,0.05,0.48,0.15,21);
  // TPad* pad11 = new TPad("pad11","This is pad11",0.52,0.05,0.96,0.15,21);



pA->SetFillColor(0);
pA->SetBorderMode(0);
pA->SetFrameFillColor(10);
pA->Draw();
pA->SetLogy(logy);
   pA->SetTopMargin(0.07);
   pA->SetBottomMargin(0.13);
   pA->SetRightMargin(0.04);
   pA->SetLeftMargin(0.16);

pB->SetFillColor(0);
pB->SetBorderMode(0);
pB->SetFrameFillColor(10);
pB->Draw();
pB->SetLogy(logy);
   pB->SetTopMargin(0.07);
   pB->SetBottomMargin(0.13);
   pB->SetRightMargin(0.04);
   pB->SetLeftMargin(0.16);

pC->SetFillColor(0);
pC->SetBorderMode(0);
pC->SetFrameFillColor(10);
pC->Draw();
pC->SetLogy(logy);
   pC->SetTopMargin(0.07);
   pC->SetBottomMargin(0.13);
   pC->SetRightMargin(0.04);
   pC->SetLeftMargin(0.16);

pD->SetFillColor(0);
pD->SetBorderMode(0);
pD->SetFrameFillColor(10);
pD->Draw();
pD->SetLogy(logy);
   pD->SetTopMargin(0.07);
   pD->SetBottomMargin(0.13);
   pD->SetRightMargin(0.04);
   pD->SetLeftMargin(0.16);


pE->SetFillColor(0);
pE->SetBorderMode(0);
pE->SetFrameFillColor(10);
pE->Draw();
pE->SetLogy(logy);
   pE->SetTopMargin(0.07);
   pE->SetBottomMargin(0.13);
   pE->SetRightMargin(0.04);
   pE->SetLeftMargin(0.16);

pF->SetFillColor(0);
pF->SetBorderMode(0);
pF->SetFrameFillColor(10);
pF->Draw();
pF->SetLogy(logy);
   pF->SetTopMargin(0.07);
   pF->SetBottomMargin(0.13);
   pF->SetRightMargin(0.04);
   pF->SetLeftMargin(0.16);

pG->SetFillColor(0);
pG->SetBorderMode(0);
pG->SetFrameFillColor(10);
pG->Draw();
pG->SetLogy(logy);
   pG->SetTopMargin(0.07);
   pG->SetBottomMargin(0.13);
   pG->SetRightMargin(0.04);
   pG->SetLeftMargin(0.16);

pH->SetFillColor(0);
pH->SetBorderMode(0);
pH->SetFrameFillColor(10);
pH->Draw();
pH->SetLogy(logy);
   pH->SetTopMargin(0.07);
   pH->SetBottomMargin(0.13);
   pH->SetRightMargin(0.04);
   pH->SetLeftMargin(0.16);


   pI->SetFillColor(0);
pI->SetBorderMode(0);
pI->SetFrameFillColor(10);
pI->Draw();
pI->SetLogy(logy);
   pI->SetTopMargin(0.07);
   pI->SetBottomMargin(0.13);
   pI->SetRightMargin(0.04);
   pI->SetLeftMargin(0.16);

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

 pA->cd();

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
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleA);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 
 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleA);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleA);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleA);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleA);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleA);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleA);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleA);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleA);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleA);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleA);
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
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 20-80 GeV");
  leg->Draw();
// *****************************************************************************

 pB->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

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
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleB);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleB);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleB);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleB);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleB);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleB);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleB);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleB);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleB);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleB);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / e1_TT->Integral(0,3);
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
  leg->SetHeader(HeaderB);
  leg->Draw();


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
 
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//


// *****************************************************************************
  c1->cd();
 pC->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
 hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
 hmax = hmaxBD; 

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
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleC);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleC);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleC);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleC);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleC);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleC);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleC);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleC);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleC);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleC);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / e1_TT->Integral(0,3);
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
  leg->SetHeader(HeaderC);
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 20-80 GeV");
  leg->Draw();

// *****************************************************************************

 pD->cd();
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
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleD);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleD);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleD);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleD);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleD);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleD);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleD);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleD);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleD);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleD);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / e1_TT->Integral(0,3);
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

// *****************************************************************************
// *****************************************************************************

  c1->Update();

// *****************************************************************************

 pE->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

 f2_Data_emu->cd();
 g2_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
 htotData->Add(g2_Data_emu, htotData, 1, 1);

 f3_Data_emu->cd();
 g3_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
 htotData->Add(g3_Data_emu, htotData, 1, 1);

 f4_Data_emu->cd();
 g4_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
 htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleE);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleE);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleE);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleE);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleE);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleE);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleE);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleE);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleE);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleE);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleE);
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
 g1_LLP = (TH1F*)gROOT->FindObject(htitleE);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleE);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleE);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleE);
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
  leg->SetHeader(HeaderE);
  leg->Draw();
//******************************************************************//
//******************************************************************//
pF->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(htitleF);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

 f2_Data_emu->cd();
 g2_Data_emu = (TH1F*)gROOT->FindObject(htitleF);
 htotData->Add(g2_Data_emu, htotData, 1, 1);

 f3_Data_emu->cd();
 g3_Data_emu = (TH1F*)gROOT->FindObject(htitleF);
 htotData->Add(g3_Data_emu, htotData, 1, 1);

 f4_Data_emu->cd();
 g4_Data_emu = (TH1F*)gROOT->FindObject(htitleF);
 htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleF);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleF);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleF);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleF);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleF);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleF);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleF);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleF);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleF);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleF);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleF);
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
 g1_LLP = (TH1F*)gROOT->FindObject(htitleF);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleF);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleF);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleF);
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
  leg->SetHeader(HeaderF);
  leg->Draw();

//*************************************************//
//************************************************//


pG->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(htitleG);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

 f2_Data_emu->cd();
 g2_Data_emu = (TH1F*)gROOT->FindObject(htitleG);
 htotData->Add(g2_Data_emu, htotData, 1, 1);

 f3_Data_emu->cd();
 g3_Data_emu = (TH1F*)gROOT->FindObject(htitleG);
 htotData->Add(g3_Data_emu, htotData, 1, 1);

 f4_Data_emu->cd();
 g4_Data_emu = (TH1F*)gROOT->FindObject(htitleG);
 htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleG);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleG);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleG);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleG);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleG);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleG);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleG);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleG);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleG);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleG);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleG);
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
 g1_LLP = (TH1F*)gROOT->FindObject(htitleG);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleG);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleG);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleG);
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
  leg->SetHeader(HeaderG);
  leg->Draw();


//********************************************µ//
//*********************************************//

pH->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(htitleH);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

 f2_Data_emu->cd();
 g2_Data_emu = (TH1F*)gROOT->FindObject(htitleH);
 htotData->Add(g2_Data_emu, htotData, 1, 1);

 f3_Data_emu->cd();
 g3_Data_emu = (TH1F*)gROOT->FindObject(htitleH);
 htotData->Add(g3_Data_emu, htotData, 1, 1);

 f4_Data_emu->cd();
 g4_Data_emu = (TH1F*)gROOT->FindObject(htitleH);
 htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleH);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleH);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleH);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleH);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleH);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleH);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleH);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleH);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleH);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleH);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleH);
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
 g1_LLP = (TH1F*)gROOT->FindObject(htitleH);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleH);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleH);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleH);
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
  leg->SetHeader(HeaderH);
  leg->Draw();
//************************************************************//
//ù**********************************************************//

pI->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(htitleI);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

 f2_Data_emu->cd();
 g2_Data_emu = (TH1F*)gROOT->FindObject(htitleI);
 htotData->Add(g2_Data_emu, htotData, 1, 1);

 f3_Data_emu->cd();
 g3_Data_emu = (TH1F*)gROOT->FindObject(htitleI);
 htotData->Add(g3_Data_emu, htotData, 1, 1);

 f4_Data_emu->cd();
 g4_Data_emu = (TH1F*)gROOT->FindObject(htitleI);
 htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 e1_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_DY->Sumw2();
 norm = 0.;
 if  ( e1_DY->Integral(0,3) > 0 ) norm = lumi / e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(htitleI);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, xsec_DY[0]*norm,0);

 f2_DY->cd();
 e2_DY = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_DY->Sumw2();
 norm = 0.;
 if  ( e2_DY->Integral(0,3) > 0 ) norm = lumi / e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(htitleI);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, xsec_DY[1]*norm, 1);

 f1_VV->cd();
 e1_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_VV->Sumw2();
 norm = 0.;
 if  ( e1_VV->Integral(0,3) > 0 ) norm = lumi / e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(htitleI);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, xsec_VV[0]*norm,0);

 f2_VV->cd();
 e2_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_VV->Sumw2();
 norm = 0.;
 if  ( e2_VV->Integral(0,3) > 0 ) norm = lumi / e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(htitleI);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, xsec_VV[1]*norm, 1);

 f3_VV->cd();
 e3_VV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_VV->Sumw2();
 norm = 0.;
 if  ( e3_VV->Integral(0,3) > 0 ) norm = lumi / e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(htitleI);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, xsec_VV[2]*norm, 1);

 f1_TTV->cd();
 e1_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TVV->Sumw2();
 norm = 0.;
 if  ( e1_TVV->Integral(0,3) > 0 ) norm = lumi / e1_TVV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(htitleI);
 g1_TVV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, xsec_TTV[0]*norm,0);

 f2_TTV->cd();
 e2_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_TVV->Sumw2();
 norm = 0.;
 if  ( e2_TVV->Integral(0,3) > 0 ) norm = lumi / e2_TVV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(htitleI);
 g2_TVV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, xsec_TTV[1]*norm, 1);

 f3_TTV->cd();
 e3_TVV = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_TVV->Sumw2();
 norm = 0.;
 if  ( e3_TVV->Integral(0,3) > 0 ) norm = lumi / e3_TVV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(htitleI);
 g3_TVV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, xsec_TTV[2]*norm, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 e1_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_ST->Sumw2();
 norm = 0.;
 if  ( e1_ST->Integral(0,3) > 0 ) norm = lumi / e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(htitleI);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, xsec_ST[0]*norm,0);

 f2_ST->cd();
 e2_ST = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_ST->Sumw2();
 norm = 0.;
 if  ( e2_ST->Integral(0,3) > 0 ) norm = lumi / e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(htitleI);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, xsec_ST[1]*norm, 1);

 f1_TT->cd();
 e1_TT = (TH1F*)gROOT->FindObject("hData_Filter");
 e1_TT->Sumw2();
 norm = 0.;
 if  ( e1_TT->Integral(0,3) > 0 ) norm = lumi / e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(htitleI);
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
 g1_LLP = (TH1F*)gROOT->FindObject(htitleI);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0.001*norm,0);

 f2_LLP->cd();
 e2_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e2_LLP->Sumw2();
 norm = 0.;
 if ( e2_LLP->Integral(0,3) > 0. ) norm = SigXsec2*lumi / e2_LLP->Integral(0,3);
 g2_LLP = (TH1F*)gROOT->FindObject(htitleI);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0.001*norm,0);

 f3_LLP->cd();
 e3_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
 e3_LLP->Sumw2();
 norm = 0.;
 if ( e3_LLP->Integral(0,3) > 0. ) norm = SigXsec3*lumi / e3_LLP->Integral(0,3);
 g3_LLP = (TH1F*)gROOT->FindObject(htitleI);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0.001*norm,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = (TH1F*)gROOT->FindObject("hData_Filter");
//            e4_LLP->Sumw2();
// 	   norm = 0.;
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleI);
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
  leg->SetHeader(HeaderI);
  leg->Draw();
//************************************µ//



  c2->cd();
  pad5->cd();

   hsolve->SetLineColor(kBlack);
 hsolve->SetFillColor(kBlack);
 hsolve->SetFillStyle(3004);
//  hsolve->SetLineStyle(1);
//  hsolve->SetLineWidth(2);
  hsolve->Draw("E4same");
  // c2->SaveAs("../"+Prod+"/"+Name+"_"+Form("%d",method)+"_v2.pdf");!
  return c1;
}