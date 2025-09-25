#include <iostream>
#include <TROOT.h>
#include "TH1.h"

TCanvas * plot(int method, TString Prod, TString Name, TString Year, TString Dmode,  TString Plots)
{
int stati=0;
bool fit= 1;
bool logy=1;

//$$
bool DATA = true;
//$$

// number of vertices:
//$$
  int nvtx = 2;
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
    hmaxBD = 1E6;
//     hmax = 1E6;   // for eta<2.4 pt>80  
//     hmaxBD = 1E9;
  }
  // TString Prod = "DATA_EMU_2017";
  //ABCD_EMu2018A 094
  // ABCD_2018A 095
  float ReScaleXS = 1.;//0.0812948*2
  // TString Dmode = "DM";//DM or SM
TString Yearcor = Year;
 if (Year == "2016PRE") Yearcor = "2016preVFP";
 if (Year == "2016POST") Yearcor = "2016";
// Dmu
TFile* f1_Data_emu  = new TFile("../../DATA_EMU_"+Year+"_31_10_2024/histofile_"+Dmode+"_OS_2p4_MuonEG-UL"+Year+"_MiniAODv2_GT36-v1_NOM.root");
    //RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1
    //RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9
    // RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17
    // RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11
    TString Campaign = "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17";
    if (Year == "2016PRE") Campaign = "RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11";
    if (Year == "2016POST") Campaign = "RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17";
    if (Year == "2017") Campaign = "RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9";
    if (Year == "2018") Campaign = "RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1";
//Mumu
 TFile* f1_DY  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_NOM.root");
 TFile* f2_DY  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NOM.root");
 TFile* f1_TT  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_NOM.root");
 TFile* f2_TT  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_NOM.root");
 TFile* f1_ST  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_NOM.root");
 TFile* f2_ST  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_NOM.root");
 TFile* f1_TTV = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8_NOM.root");
 TFile* f2_TTV = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8_NOM.root");
 TFile* f3_TTV = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8_NOM.root");
 TFile* f1_VV  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_NOM.root");
 TFile* f2_VV  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_NOM.root");
 TFile* f3_VV  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/histofile_"+Dmode+"_OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_NOM.root");
 

TString MSUON[3] = {"200","300","400"};
  TString MNEU[3] = {"180","180","300"};
  TString CTAU[3] = {"100","100","100"};

//signal
 TFile* f1_LLP = new TFile("../../Signal_"+Year+"/histofile_DM_OS_2p4_RPV_"+Yearcor+"_smu"+MSUON[0]+"_neu"+MNEU[0]+"_ctau"+CTAU[0]+"_NOM.root");
 TFile* f2_LLP = new TFile("../../Signal_"+Year+"/histofile_DM_OS_2p4_RPV_"+Yearcor+"_smu"+MSUON[1]+"_neu"+MNEU[1]+"_ctau"+CTAU[1]+"_NOM.root");
 TFile* f3_LLP = new TFile("../../Signal_"+Year+"/histofile_DM_OS_2p4_RPV_"+Yearcor+"_smu"+MSUON[2]+"_neu"+MNEU[2]+"_ctau"+CTAU[2]+"_NOM.root");

 TString DATAFILE[1] = {"MuonEG-UL"+Year+"_MiniAODv2_GT36-v1_NOM_"
};
if (Year == "2018") DATAFILE[0] = "MuonEG-UL"+Year+"_MiniAODv2_GT36-v1_NOM_";


TString MCFILE[12] = {"DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_NOM_",
                  "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_NOM_",
                  "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_NOM_",
                  "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_NOM_",
                  "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_NOM_",
                  "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_NOM_",
                  "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8_NOM_",
                  "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8_NOM_",
                  "TTWW_TuneCP5_13TeV-madgraph-pythia8_NOM_",
                  "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_NOM_",
                  "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_NOM_",
                  "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_NOM_",
};



TString LLPFILE[3] = {"RPV_"+Yearcor+"_smu"+MSUON[0]+"_neu"+MNEU[0]+"_ctau"+CTAU[0]+"_NOM_",
                  "RPV_"+Yearcor+"_smu"+MSUON[1]+"_neu"+MNEU[1]+"_ctau"+CTAU[1]+"_NOM_",
                  "RPV_"+Yearcor+"_smu"+MSUON[2]+"_neu"+MNEU[2]+"_ctau"+CTAU[2]+"_NOM_"

};


//  if ( nvtx == 1 ) xtitle = "vertex BDT score"; 
 TString ytitle = "Events"; 
//  int nbin = 26; 
//  float xmin = -1.04;
//  float xmax =  1.04;
//$$
 TString HeaderCMS = "CMS";

 if (Year == "2016") HeaderCMS = "2016                            36.3 fb^{-1} (13 TeV)";
 if (Year == "2017") HeaderCMS = "2017                            41.5 fb^{-1} (13 TeV)";
 if (Year == "2018") HeaderCMS = "2018                            59.8 fb^{-1} (13 TeV)";


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

  
 int Method = method;
if (Method == 0)
  {
    htitleA = "hData_CRlowlowpt_2Vtx_Mass_";
    htitleB = "hData_CRlowlowpt_TLVtx_Mass_";
    htitleC = "hData_CRlooselowlowpt_2Vtx_Mass_";

    htitleD = "hData_CRlowpt_2Vtx_Mass_";
    htitleE = "hData_CRlowpt_TLVtx_Mass_";//ok
    htitleF = "hData_CRlooselowpt_2Vtx_Mass_";//ok

    htitleG = "hData_Hemi_2Vtx_Mass_";//ok
    htitleH = "hData_Hemi_TLVtx_Mass_";//ok
    htitleI = "hData_CRloose_2Vtx_Mass_";//ok
    nbin = 25; 
    xmin = 0;
    xmax =  100;
    HeaderA = "TT + 20<pt<80";
    HeaderB = "TL + 20<pt<80";
    HeaderC = "LL + 20<pt<80";

    HeaderD = "TT + 20<pt_{2}<80 & pt_{1}>80";
    HeaderE = "TL + 20<pt_{2}<80 & pt_{1}>80";
    HeaderF = "LL + 20<pt_{2}<80 & pt_{1}>80";
    
    HeaderG = "TT + pt_{1},pt_{2}>80 ";
    HeaderH = "TL + pt_{1},pt_{2}>80";
    HeaderI = "LL + pt_{1},pt_{2}>80";
    HeaderNVtx = "2 Vtx"; 
     xtitle = "Vtx Mass";
  }
if (Method == 1)
  {
    htitleA = "hData_CRlowlowpt_2VtxAll_Mass_";
    htitleB = "hData_CRlowlowpt_TLVtxAll_Mass_";
    htitleC = "hData_CRlooselowlowpt_2VtxAll_Mass_";

    htitleD = "hData_CRlowpt_2VtxAll_Mass_";
    htitleE = "hData_CRlowpt_TLVtxAll_Mass_";//ok
    htitleF = "hData_CRlooselowpt_2VtxAll_Mass_";//ok

    htitleG = "hData_Hemi_2VtxAll_Mass_";//ok
    htitleH = "hData_Hemi_TLVtxAll_Mass_";//ok
    htitleI = "hData_CRloose_2VtxAll_Mass_";//ok
    nbin = 25; 
    xmin = 0;
    xmax =  100;
    HeaderA = "TT + 20<pt<80";
    HeaderB = "TL + 20<pt<80";
    HeaderC = "LL + 20<pt<80";

    HeaderD = "TT + 20<pt_{2}<80 & pt_{1}>80";
    HeaderE = "TL + 20<pt_{2}<80 & pt_{1}>80";
    HeaderF = "LL + 20<pt_{2}<80 & pt_{1}>80";
    
    HeaderG = "TT + pt_{1},pt_{2}>80 ";
    HeaderH = "TL + pt_{1},pt_{2}>80";
    HeaderI = "LL + pt_{1},pt_{2}>80";
    HeaderNVtx = "2 Vtx All"; 
     xtitle = "Vtx Mass";
  }


if (Method == 2)
  {
    htitleA = "hData_CRlowlowpt_2Vtx_SumtrackWeight_";
    htitleB = "hData_CRlowlowpt_TLVtx_SumtrackWeight_";
    htitleC = "hData_CRlooselowlowpt_2Vtx_SumtrackWeight_";

    htitleD = "hData_CRlowpt_2Vtx_SumtrackWeight_";
    htitleE = "hData_CRlowpt_TLVtx_SumtrackWeight_";//ok
    htitleF = "hData_CRlooselowpt_2Vtx_SumtrackWeight_";//ok

    htitleG = "hData_Hemi_2Vtx_SumtrackWeight_";//ok
    htitleH = "hData_Hemi_TLVtx_SumtrackWeight_";//ok
    htitleI = "hData_CRloose_2Vtx_SumtrackWeight_";//ok
    nbin = 40; 
    xmin = 0;
    xmax = 40;
    HeaderA = "TT + 20<pt<80";
    HeaderB = "TL + 20<pt<80";
    HeaderC = "LL + 20<pt<80";

    HeaderD = "TT + 20<pt_{2}<80 & pt_{1}>80";
    HeaderE = "TL + 20<pt_{2}<80 & pt_{1}>80";
    HeaderF = "LL + 20<pt_{2}<80 & pt_{1}>80";
    
    HeaderG = "TT + pt_{1},pt_{2}>80 ";
    HeaderH = "TL + pt_{1},pt_{2}>80";
    HeaderI = "LL + pt_{1},pt_{2}>80";
    HeaderNVtx = "2 Vtx";
    xtitle = "SumtrackWeight";
  }

if (Method == 3)
  {
    htitleA = "hData_CRlowlowpt_2VtxAll_SumtrackWeight_";
    htitleB = "hData_CRlowlowpt_TLVtxAll_SumtrackWeight_";
    htitleC = "hData_CRlooselowlowpt_2VtxAll_SumtrackWeight_";

    htitleD = "hData_CRlowpt_2VtxAll_SumtrackWeight_";
    htitleE = "hData_CRlowpt_TLVtxAll_SumtrackWeight_";//ok
    htitleF = "hData_CRlooselowpt_2VtxAll_SumtrackWeight_";//ok

    htitleG = "hData_Hemi_2VtxAll_SumtrackWeight_";//ok
    htitleH = "hData_Hemi_TLVtxAll_SumtrackWeight_";//ok
    htitleI = "hData_CRloose_2VtxAll_SumtrackWeight_";//ok
    nbin = 40; 
    xmin = 0;
    xmax = 40;
    HeaderA = "TT + 20<pt<80";
    HeaderB = "TL + 20<pt<80";
    HeaderC = "LL + 20<pt<80";

    HeaderD = "TT + 20<pt_{2}<80 & pt_{1}>80";
    HeaderE = "TL + 20<pt_{2}<80 & pt_{1}>80";
    HeaderF = "LL + 20<pt_{2}<80 & pt_{1}>80";
    
    HeaderG = "TT + pt_{1},pt_{2}>80 ";
    HeaderH = "TL + pt_{1},pt_{2}>80";
    HeaderI = "LL + pt_{1},pt_{2}>80";
    HeaderNVtx = "2 VtxAll";
    xtitle = "SumtrackWeight";
  }



    //-----------------------------------------------------------//
// xsec in pb

  TLegend* leg;
    
// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

  TCanvas *c2 = new TCanvas("c2", "plots",200,0,700,700);
  c2->SetFillColor(10);
  c2->SetFillStyle(4000);
  c2->SetBorderSize(2);

c1->cd();
TPad* pA = new TPad("pA","This is pad1",0.01,0.66,0.33,0.99,21);
TPad* pB = new TPad("pB","This is pad2",0.01,0.34,0.33,0.65,21);
TPad* pC = new TPad("pC","This is pad3",0.01,0.01,0.33,0.33,21);

TPad* pD = new TPad("pD","This is pad4",0.34,0.66,0.66,0.99,21);
TPad* pE = new TPad("pE","This is pE",0.34,0.34,0.66,0.65,21);
TPad* pF = new TPad("pF","This is pF",0.34,0.01,0.66,0.33,21);

TPad* pG = new TPad("pG","This is pG",0.67,0.66,0.99,0.99,21);
TPad* pH = new TPad("pH","This is pH",0.67,0.34,0.99,0.65,21);
TPad* pI = new TPad("pI","This is pI",0.67,0.01,0.99,0.33,21);

// 

  // TPad* pad8 = new TPad("pad8","This is pad8",0.04,0.5,0.48,0.65,21);
  // TPad* pad9 = new TPad("pad9","This is pad9",0.52,0.5,0.96,0.65,21);
  // TPad* pad10 = new TPad("pad10","This is pad10",0.04,0.05,0.48,0.15,21);
  // TPad* pad11 = new TPad("pad11","This is pad11",0.52,0.05,0.96,0.15,21);


  c2->cd();
  TPad* padG = new TPad("padG","This is padG",0.01,0.55,0.33,0.96,21);
  TPad* padGRatio = new TPad("padGRatio","This is padGRatio",0.01,0.35,0.33,0.55,21);
  TPad* padGInte = new TPad("padGInte","This is padGInte",0.01,0.1,0.33,0.35,21);

  TPad* padH1 = new TPad("padH1","This is padH1",0.34,0.66,0.66,0.99,21);
  TPad* padH1Ratio = new TPad("padH1Ratio","This is padH1Ratio",0.34,0.55,0.66,0.66,21);
  TPad* padH2 = new TPad("padH2","This is padH2",0.34,0.22,0.66,0.55,21);
  TPad* padH2Ratio = new TPad("padH2Ratio","This is padH2Ratio",0.34,0.05,0.66,0.21,21);

  TPad* padD1 = new TPad("padD1","This is padD1",0.67,0.66,0.99,0.99,21);
    TPad* padD1Ratio = new TPad("padD1Ratio","This is padD1Ratio",0.67,0.55,0.99,0.66,21);
  TPad* padD2 = new TPad("padD2","This is padD2",0.67,0.22,0.99,0.55,21);
TPad* padD2Ratio = new TPad("padD2Ratio","This is padD2Ratio",0.67,0.05,0.99,0.21,21);
  
  padG->SetFillColor(0);
padG->SetBorderMode(0);
padG->SetFrameFillColor(10);
padG->Draw();
padG->SetLogy(logy);
   padG->SetTopMargin(0.07);
   padG->SetBottomMargin(0.13);
   padG->SetRightMargin(0.04);
   padG->SetLeftMargin(0.16);

padGRatio->SetFillColor(0);
padGRatio->SetBorderMode(0);
padGRatio->SetFrameFillColor(10);
padGRatio->Draw();
padGRatio->SetLogy(0);
   padGRatio->SetTopMargin(0.07);
   padGRatio->SetBottomMargin(0.13);
   padGRatio->SetRightMargin(0.04);
   padGRatio->SetLeftMargin(0.16);

padGInte->SetFillColor(0);
padGInte->SetBorderMode(0);
padGInte->SetFrameFillColor(10);
padGInte->Draw();
padGInte->SetLogy(0);
   padGInte->SetTopMargin(0.07);
   padGInte->SetBottomMargin(0.13);
   padGInte->SetRightMargin(0.04);
   padGInte->SetLeftMargin(0.16);
 


   padH1->SetFillColor(0);
padH1->SetBorderMode(0);
padH1->SetFrameFillColor(10);
padH1->Draw();
padH1->SetLogy(logy);
   padH1->SetTopMargin(0.07);
   padH1->SetBottomMargin(0.13);
   padH1->SetRightMargin(0.04);
   padH1->SetLeftMargin(0.16);

      padH1Ratio->SetFillColor(0);
padH1Ratio->SetBorderMode(0);
padH1Ratio->SetFrameFillColor(10);
padH1Ratio->Draw();
padH1Ratio->SetLogy(0);
   padH1Ratio->SetTopMargin(0.07);
   padH1Ratio->SetBottomMargin(0.13);
   padH1Ratio->SetRightMargin(0.04);
   padH1Ratio->SetLeftMargin(0.16);

padH2->SetFillColor(0);
padH2->SetBorderMode(0);
padH2->SetFrameFillColor(10);
padH2->Draw();
padH2->SetLogy(logy);
   padH2->SetTopMargin(0.07);
   padH2->SetBottomMargin(0.13);
   padH2->SetRightMargin(0.04);
   padH2->SetLeftMargin(0.16);

   padH2Ratio->SetFillColor(0);
padH2Ratio->SetBorderMode(0);
padH2Ratio->SetFrameFillColor(10);
padH2Ratio->Draw();
padH2Ratio->SetLogy(0);
   padH2Ratio->SetTopMargin(0.07);
   padH2Ratio->SetBottomMargin(0.13);
   padH2Ratio->SetRightMargin(0.04);
   padH2Ratio->SetLeftMargin(0.16);

padD1->SetFillColor(0);
padD1->SetBorderMode(0);
padD1->SetFrameFillColor(10);
padD1->Draw();
padD1->SetLogy(logy);
   padD1->SetTopMargin(0.07);
   padD1->SetBottomMargin(0.13);
   padD1->SetRightMargin(0.04);
   padD1->SetLeftMargin(0.16);

   padD1Ratio->SetFillColor(0);
padD1Ratio->SetBorderMode(0);
padD1Ratio->SetFrameFillColor(10);
padD1Ratio->Draw();
padD1Ratio->SetLogy(0);
   padD1Ratio->SetTopMargin(0.07);
   padD1Ratio->SetBottomMargin(0.13);
   padD1Ratio->SetRightMargin(0.04);
   padD1Ratio->SetLeftMargin(0.16);

padD2->SetFillColor(0);
padD2->SetBorderMode(0);
padD2->SetFrameFillColor(10);
padD2->Draw();
padD2->SetLogy(logy);
   padD2->SetTopMargin(0.07);
   padD2->SetBottomMargin(0.13);
   padD2->SetRightMargin(0.04);
   padD2->SetLeftMargin(0.16);

      padD2Ratio->SetFillColor(0);
padD2Ratio->SetBorderMode(0);
padD2Ratio->SetFrameFillColor(10);
padD2Ratio->Draw();
padD2Ratio->SetLogy(0);
   padD2Ratio->SetTopMargin(0.07);
   padD2Ratio->SetBottomMargin(0.13);
   padD2Ratio->SetRightMargin(0.04);
   padD2Ratio->SetLeftMargin(0.16);
   
  c1->cd();


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

 TH1F* hsolveD1  = new TH1F("hsolveD1","",nbin,xmin,xmax);
 TH1F* hsolveD2  = new TH1F("hsolveD2","",nbin,xmin,xmax);
 TH1F* hsolveG   = new TH1F("hsolveG","",nbin,xmin,xmax);
 TH1F* hsolveH1  = new TH1F("hsolveH1","",nbin,xmin,xmax);
 TH1F* hsolveH2  = new TH1F("hsolveH2","",nbin,xmin,xmax);
 hsolveD1->Sumw2();
 hsolveD2->Sumw2();
 hsolveG->Sumw2();
 hsolveH1->Sumw2();
 hsolveH2->Sumw2();

 TH1F* htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 TH1F* htotData  = new TH1F("htotData","",nbin,xmin,xmax);

TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

  htotData->Sumw2();
  htotMC->Sumw2();
  hDataMC->Sumw2();
 
 
 TH1F* g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);//ok
 TH1F* g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);//ok
 TH1F*  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

 
 TH1F* g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);//ok
 TH1F* g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);//ok
 TH1F* g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);//ok
 TH1F*  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 
 TH1F* g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);//ok
 TH1F* g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);//ok
 TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 TH1F* g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);//ok
 TH1F*  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 TH1F* g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);//ok
 TH1F* g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);//ok 
 TH1F* g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);//ok
 TH1F*  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
  
 TH1F* g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleA);//ok
 TH1F* h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 TH1F* g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleA);//ok
 TH1F* h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);  
 TH1F* g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleA);//ok
 TH1F* h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);


// *****************************************************************************

 pA->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 TH1F* g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleA);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  TH1F* g2_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//  TH1F* g3_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//  TH1F* g4_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
//  htotData->Add(g4_Data_emu, htotData, 1, 1);

//------x);

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 
 f2_DY->cd();
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, 1, 1);


 f1_VV->cd();
 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);
 g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);
 g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);
 g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);
// std::cout<<"TT: "<<1<<"with norm "<<norm<<std::endl;
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
   
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleA);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0,0);

 f2_LLP->cd();
  g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleA);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0,0);

 f3_LLP->cd();
  
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleA);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = 
//            e4_LLP->Sumw2();
// 	   
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleA);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

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

  hsolveD1->Add(hsolveD1, htotData, 0., 1.);
 hsolveD2->Add(hsolveD2, htotData, 0., 1.);
 hsolveG->Add(hsolveG, htotData, 0., 1.);

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
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
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");
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

       leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("A");
  leg->Draw();
// *****************************************************************************

 pB->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleB);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  g2_Data_emu = (TH1F*)gROOT->FindObject(htitleB);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//  g3_Data_emu = (TH1F*)gROOT->FindObject(htitleB);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//  g4_Data_emu = (TH1F*)gROOT->FindObject(htitleB);
//  htotData->Add(g4_Data_emu, htotData, 1, 1);



 f1_DY->cd();
 
 
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 
 

 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleB);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
 
 
  //e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleB);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 
 

 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleB);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 
 
 
  // e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleB);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 
 
 
  //e1_TTV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleB);
 g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 
 
  //e2_TTV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleB);
 g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 
 
 
  //e3_TTV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleB);
 g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 
 
  // e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleB);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 //e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleB);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
 
 //e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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
htotData->SetMarkerStyle(20);
  htotData->SetFillColor(kBlack);
  htotData->SetLineColor(kBlack);
  htotData->SetLineStyle(1);
  htotData->SetLineWidth(1);

 f1_LLP->cd();
  
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleB);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0,0);

 f2_LLP->cd();
  
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleB);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0,0);

 f3_LLP->cd();
  
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleB);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = 
//            e4_LLP->Sumw2();
// 	   
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleB);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

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
 hsolveD1->Divide(hsolveD1, htotData, 1., 1.);
 hsolveH2->Add(hsolveD2, htotData, 0., 1.);

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
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
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");

// leg->AddEntry(hsolve," Prediction","FE4");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();


  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderB);
  leg->Draw();

     leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("B");
  leg->Draw();

 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//

// *****************************************************************************
  // c1->cd();
 pC->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
 hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
 hmax = hmaxBD; 

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleC);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  g2_Data_emu = (TH1F*)gROOT->FindObject(htitleC);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//  g3_Data_emu = (TH1F*)gROOT->FindObject(htitleC);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//  g4_Data_emu = (TH1F*)gROOT->FindObject(htitleC);
//  htotData->Add(g4_Data_emu, htotData, 1, 1);
 
 f1_DY->cd();
 
 
 
  // e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 
 
 
  //e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleC);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
 
 
 
  // e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleC);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 
 
 
  // e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleC);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 
 
 
  // e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleC);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 
 
 
  // e1_TTV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleC);
 g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 
 
  // e2_TTV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleC);
 g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 
 
 
  // e3_TTV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleC);
 g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 
 
 
  //e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleC);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 
 
  //e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleC);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
 
 
 
  //e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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

 
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleC);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0,0);

 f2_LLP->cd();
  

 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleC);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0,0);

 f3_LLP->cd();
  
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleC);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = 
//            e4_LLP->Sumw2();
// 	   
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleC);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

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

  hsolveG->Divide(hsolveG, htotData, 1., 1.);
 hsolveD2->Divide(hsolveD2, htotData, 1., 1.);
 hsolveH2->Divide(hsolveH2, htotData, 1., 1.);

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
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
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderC);
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 20-80 GeV");
  leg->Draw();
     leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("C");
  leg->Draw();
// *****************************************************************************

 pD->cd();
 hmax = hmaxBD; 

//Security for blinding
 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleD);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  g2_Data_emu = (TH1F*)gROOT->FindObject(htitleD);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//  g3_Data_emu = (TH1F*)gROOT->FindObject(htitleD);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//  g4_Data_emu = (TH1F*)gROOT->FindObject(htitleD);
//  htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 //e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleD);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 // e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleD);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
 
 
  //e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleD);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
  
 // e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleD);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 
 //e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleD);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 
 //e1_TTV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleD);
 g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 // e2_TTV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleD);
 g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 
 
  // e3_TTV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleD);
 g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 
 
 
  //e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleD);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 
 
  //e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleD);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
 
 //e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleD);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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
  
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleD);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0,0);

 f2_LLP->cd();
  
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleD);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0,0);

 f3_LLP->cd();
  
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleD);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = 
//            e4_LLP->Sumw2();
// 	   
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleD);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

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

//  hsolve->Multiply(hsolve, htotData, 1., 1.);

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
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
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderD);
  leg->Draw();

     leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("D");
  leg->Draw();
  TH1F* htotDataD  = new TH1F("htotDataD","",nbin,xmin,xmax);
htotDataD->Sumw2();
htotDataD->Add(htotData, htotDataD, 1, 0);
// std::cout<<"here 1 "<<std::endl;
  c2->cd();
// c2->ls();
padD1->cd();
  // std::cout<<"here 2 "<<std::endl;
  htotDataD->Draw("E1"); 
 htotDataD->SetFillColor(kBlack);
 htotDataD->SetLineColor(kBlack);
 htotDataD->SetLineStyle(1);
 htotDataD->SetLineWidth(1);
 htotDataD->SetTickLength(0.03, "YZ");
 htotDataD->SetTickLength(0.03,"X");
 htotDataD->SetLabelOffset(0.015,"X");
 htotDataD->SetLabelOffset(0.007,"Y");
 htotDataD->SetLabelSize(0.045, "XYZ");
 htotDataD->SetLabelFont(42, "XYZ"); 
 htotDataD->SetTitleSize(0.045, "XYZ"); 
 htotDataD->SetTitleFont(42, "XYZ");
 htotDataD->SetTitleOffset(1.2,"X"); 
 htotDataD->SetTitleOffset(1.3,"Y");
 htotDataD->GetXaxis()->SetTitle(xtitle);
 htotDataD->GetXaxis()->SetTitleColor(1);
 htotDataD->GetYaxis()->SetTitle(ytitle);
 htotDataD->GetYaxis()->SetTitleColor(1);
 htotDataD->SetNdivisions(509,"XYZ");
 htotDataD->SetMinimum(hmin); 
 htotDataD->SetMaximum(hmax); 
 htotDataD->SetMarkerStyle(20);
 htotDataD->SetMarkerSize(1);





 leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
  leg->Draw();


  leg = new TLegend(0.64,0.7,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotDataD, " e#mu data","PE1");
  leg->AddEntry(hsolveD1, " Prediction","FE4" );

  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderD);
  leg->Draw();

   leg = new TLegend(0.2,0.70,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader("D=A*E/B");
  leg->Draw();



 padD2->cd();
  htotDataD->Draw("E1"); 
 htotDataD->SetFillColor(kBlack);
 htotDataD->SetLineColor(kBlack);
 htotDataD->SetLineStyle(1);
 htotDataD->SetLineWidth(1);
 htotDataD->SetTickLength(0.03, "YZ");
 htotDataD->SetTickLength(0.03,"X");
 htotDataD->SetLabelOffset(0.015,"X");
 htotDataD->SetLabelOffset(0.007,"Y");
 htotDataD->SetLabelSize(0.045, "XYZ");
 htotDataD->SetLabelFont(42, "XYZ"); 
 htotDataD->SetTitleSize(0.045, "XYZ"); 
 htotDataD->SetTitleFont(42, "XYZ");
 htotDataD->SetTitleOffset(1.2,"X"); 
 htotDataD->SetTitleOffset(1.3,"Y");
 htotDataD->GetXaxis()->SetTitle(xtitle);
 htotDataD->GetXaxis()->SetTitleColor(1);
 htotDataD->GetYaxis()->SetTitle(ytitle);
 htotDataD->GetYaxis()->SetTitleColor(1);
 htotDataD->SetNdivisions(509,"XYZ");
 htotDataD->SetMinimum(hmin); 
 htotDataD->SetMaximum(hmax); 
 htotDataD->SetMarkerStyle(20);
 htotDataD->SetMarkerSize(1);

 leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
  leg->Draw();


  leg = new TLegend(0.64,0.7,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotDataD, " e#mu data","PE1");
  leg->AddEntry(hsolveD2, " Prediction","FE4" );

  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderD);
  leg->Draw();

   leg = new TLegend(0.2,0.70,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader("D=A*F/C");
  leg->Draw();

   
 c1->cd();
// *****************************************************************************
// *****************************************************************************

  // c1->Update();

// *****************************************************************************

 pE->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleE);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  g2_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//  g3_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//  g4_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
//  htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 // e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleE);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 
 //e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleE);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
 // e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleE);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 // e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleE);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 
 // e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleE);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 
 // e1_TTV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleE);
 g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 // e2_TTV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleE);
 g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 
 // e3_TTV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleE);
 g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 
 //e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleE);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 //e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleE);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
 
 //e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleE);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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
  
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleE);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0,0);

 f2_LLP->cd();
  
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleE);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0,0);

 f3_LLP->cd();
  
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleE);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = 
//            e4_LLP->Sumw2();
// 	   
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleE);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

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

 hsolveD1->Multiply(hsolveD1, htotData, 1., 1.);
 hsolveH1->Add(hsolveH1, htotData, 0., 1.);

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("2018                                  59.8 fb^{-1} (13 TeV)");
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
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderE);
    leg->Draw();
       leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("E");
  leg->Draw();
//******************************************************************//
//******************************************************************//
pF->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);


 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleF);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  g2_Data_emu = (TH1F*)gROOT->FindObject(htitleF);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//  g3_Data_emu = (TH1F*)gROOT->FindObject(htitleF);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//  g4_Data_emu = (TH1F*)gROOT->FindObject(htitleF);
//  htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 // e1_DY->Integral(0,3);
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleF);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 
 //e2_DY->Integral(0,3);
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleF);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
 
 // e1_VV->Integral(0,3);
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleF);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 
 
 // e2_VV->Integral(0,3);
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleF);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 // e3_VV->Integral(0,3);
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleF);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 
 // e1_TTV->Integral(0,3);
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleF);
 g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 
 // e2_TTV->Integral(0,3);
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleF);
 g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 
 
 
 // e3_TTV->Integral(0,3);
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleF);
 g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 
 
 
 //e1_ST->Integral(0,3);
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleF);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 //e2_ST->Integral(0,3);
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleF);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
 
 
 //e1_TT->Integral(0,3);
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleF);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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
  
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleF);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0,0);

 f2_LLP->cd();
  
 
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleF);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0,0);

 f3_LLP->cd();
  
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleF);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = 
//            e4_LLP->Sumw2();
// 	   
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleF);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

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

 hsolveD2->Multiply(hsolveD2, htotData, 1., 1.);
 hsolveH1->Divide(hsolveH1, htotData, 1., 1.);

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
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
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderF);
  leg->Draw();

     leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("F");
  leg->Draw();
//*************************************************//
//************************************************//


pG->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
//securiryt for blinding
 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleG);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  g2_Data_emu = (TH1F*)gROOT->FindObject(htitleG);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//  g3_Data_emu = (TH1F*)gROOT->FindObject(htitleG);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//  g4_Data_emu = (TH1F*)gROOT->FindObject(htitleG);
//  htotData->Add(g4_Data_emu, htotData, 1, 1);

TH1F * htotDataG = new TH1F("htotDataG","",nbin,xmin,xmax);
htotDataG->Add(htotDataG,htotData,0,1);

 
 f1_DY->cd();
 
 
 
 
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleG);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 
 
 
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleG);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
 
 
 
 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleG);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 
 
 
 
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleG);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 
 
 
 
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleG);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 
 
 
 
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleG);
 g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 
 
 
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleG);
 g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 
 
 
 
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleG);
 g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 
 
 
 
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleG);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 
 
 
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleG);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
 
 
 
 
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleG);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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
  
 
 
 
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleG);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0,0);

 f2_LLP->cd();
  
 
 
 
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleG);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0,0);

 f3_LLP->cd();
  
 
 
 
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleG);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = 
//            e4_LLP->Sumw2();
// 	   
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleG);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

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

//  hsolve->Multiply(hsolve, htotData, 1., 1.);

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
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
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderG);
  leg->Draw();

  leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("G");
  leg->Draw();
c2->cd();
padG->cd();
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

   leg = new TLegend(0.2,0.70,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader("G=A*I/C");
  leg->Draw();
c1->cd();
//*********************************************//

pH->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

//security for blinding
 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleH);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  g2_Data_emu = (TH1F*)gROOT->FindObject(htitleH);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//  g3_Data_emu = (TH1F*)gROOT->FindObject(htitleH);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//  g4_Data_emu = (TH1F*)gROOT->FindObject(htitleH);
//  htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 
 
 
 
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleH);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 
 
 
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleH);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
 
 
 
 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleH);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 
 
 
 
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleH);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 
 
 
 
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleH);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 
 
 
 
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleH);
 g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 
 
 
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleH);
 g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 
 
 
 
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleH);
 g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 
 
 
 
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleH);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 
 
 
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleH);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
 
 
 
 
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleH);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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
  
 
 
 
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleH);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0,0);

 f2_LLP->cd();
  
 
 
 
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleH);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0,0);

 f3_LLP->cd();
  
 
 
 
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleH);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = 
//            e4_LLP->Sumw2();
// 	   
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleH);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

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

//  hsolve->Multiply(hsolve, htotData, 1., 1.);

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
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
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderH);
  leg->Draw();
    leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("H");
  leg->Draw();

  TH1F* htotDataH  = new TH1F("htotDataH","",nbin,xmin,xmax);
htotDataH->Sumw2();
htotDataH->Add(htotData, htotDataH, 1, 0);
  c2->cd();
padH1->cd();
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

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
  leg->Draw();


  leg = new TLegend(0.64,0.7,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " e#mu data","PE1");
  leg->AddEntry(hsolveH1, " Prediction","FE4" );

  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderH);
  leg->Draw();

   leg = new TLegend(0.2,0.70,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader("H=E*I/F");
  leg->Draw();




 padH2->cd();
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
  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
  leg->Draw();


  leg = new TLegend(0.64,0.7,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " e#mu data","PE1");
  leg->AddEntry(hsolveH2, " Prediction","FE4" );

  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderH);
  leg->Draw();

   leg = new TLegend(0.2,0.70,0.5,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader("H=B*I/C");
  leg->Draw();

  
c1->cd();
//************************************************************//
//ù**********************************************************//

pI->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleI);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  g2_Data_emu = (TH1F*)gROOT->FindObject(htitleI);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//  g3_Data_emu = (TH1F*)gROOT->FindObject(htitleI);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//  g4_Data_emu = (TH1F*)gROOT->FindObject(htitleI);
//  htotData->Add(g4_Data_emu, htotData, 1, 1);

 
 f1_DY->cd();
 
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleI);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 
 
 
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleI);
 g2_DY->Sumw2();
 h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
 
 
 
 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleI);
 g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 
 
 
 
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleI);
 g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 
 
 
 
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleI);
 g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 
 
 
 
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleI);
 g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 
 
 
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleI);
 g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 
 
 
 
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleI);
 g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 
 
 
 
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleI);
 g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 
 
 
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleI);
 g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
 
 
 
 
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleI);
 g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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
  
 
 
 
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleI);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 0,0);

 f2_LLP->cd();
  
 
 
 
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleI);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 0,0);

 f3_LLP->cd();
  
 
 
 
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleI);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 0,0);

//            f4_LLP->cd();
//      TH1F* e4_LLP = 
//            e4_LLP->Sumw2();
// 	   
// 	   if ( e4_LLP->Integral(0,3) > 0. ) norm = SigXsec4*lumi / e4_LLP->Integral(0,3);
//      TH1F* g4_LLP = (TH1F*)gROOT->FindObject(htitleI);
//            g4_LLP->Sumw2();
//      TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
//            h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

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

 hsolveG->Multiply(hsolveG, htotData, 1., 1.);
 hsolveH1->Multiply(hsolveH1, htotData, 1., 1.);
 hsolveH2->Multiply(hsolveH2, htotData, 1., 1.);

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
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
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderI);
  leg->Draw();

         leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("I");
  leg->Draw();
//************************************µ//



  // !! -------------------!! //
c2->cd();
padG->cd();
//   htotData->Draw("E1"); 
//  htotData->SetFillColor(kBlack);
//  htotData->SetLineColor(kBlack);
//  htotData->SetLineStyle(1);
//  htotData->SetLineWidth(1);
//  htotData->SetTickLength(0.03, "YZ");
//  htotData->SetTickLength(0.03,"X");
//  htotData->SetLabelOffset(0.015,"X");
//  htotData->SetLabelOffset(0.007,"Y");
//  htotData->SetLabelSize(0.045, "XYZ");
//  htotData->SetLabelFont(42, "XYZ"); 
//  htotData->SetTitleSize(0.045, "XYZ"); 
//  htotData->SetTitleFont(42, "XYZ");
//  htotData->SetTitleOffset(1.2,"X"); 
//  htotData->SetTitleOffset(1.3,"Y");
//  htotData->GetXaxis()->SetTitle(xtitle);
//  htotData->GetXaxis()->SetTitleColor(1);
//  htotData->GetYaxis()->SetTitle(ytitle);
//  htotData->GetYaxis()->SetTitleColor(1);
//  htotData->SetNdivisions(509,"XYZ");
//  htotData->SetMinimum(hmin); 
//  htotData->SetMaximum(hmax); 
//  htotData->SetMarkerStyle(20);
//  htotData->SetMarkerSize(1);


//    hsolveG->SetLineColor(kBlack);
//  hsolveG->SetFillColor(kBlack);
//  hsolveG->SetFillStyle(3004);
// //  hsolve->SetLineStyle(1);
// //  hsolve->SetLineWidth(2);
//   hsolveG->Draw("PE1same");//E4same


 hsolveG->SetLineColor(kAzure+6);
 hsolveG->SetFillColor(kAzure+6);
 hsolveG->SetFillStyle(3001);
 hsolveG->SetLineStyle(1);
 hsolveG->SetLineWidth(1);
 hsolveG->Draw("sameE2"); 

  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

 
  leg = new TLegend(0.64,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " e#mu data","PE1");
  // leg->AddEntry(hsolveE," Prediction","LPE");
leg->AddEntry(hsolveG," Prediction","FE4");
  leg->Draw();

  leg = new TLegend(0.2,0.85,0.5,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  leg = new TLegend(0.2,0.80,0.5,0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderG);
  leg->Draw();

padGRatio->cd();

TH1F* ratio  = new TH1F("ratio","",nbin,xmin,xmax);
ratio->Sumw2();
for (int i = 0; i< nbin; i++)
  {
    double data = htotDataG->GetBinContent(i);
    double pred = hsolveG->GetBinContent(i);
    double sigma = hsolveG->GetBinError(i);
    // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
    if (data == 0 || pred == 0) {ratio->AddBinContent(i,0);}    //when data == 0 and there is no prediction => you need this
    else
      {
        // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
        double pull = (data - pred) / sqrt(data + sigma*sigma);
        // double pull = data/pred;
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
 ratio->SetLabelSize(0.1, "XYZ");
 ratio->SetLabelFont(42, "XYZ"); 
 ratio->SetTitleSize(0.1, "XYZ"); 
 ratio->SetTitleFont(42, "XYZ");
 ratio->SetTitleOffset(1.2,"X"); 
 ratio->SetTitleOffset(0.5,"Y");
 ratio->GetXaxis()->SetTitle(xtitle);
 ratio->GetXaxis()->SetTitleColor(1);
 ratio->GetYaxis()->SetTitle("Pulls");
 ratio->GetYaxis()->SetTitleColor(1);
 ratio->SetNdivisions(509,"XYZ");
 ratio->SetMaximum(3); 
 ratio->SetMinimum(-3); 
 ratio->SetMarkerStyle(20);
 ratio->SetMarkerSize(1);

padGInte->cd();

// // intégrale à droite
TH1F* Inte  = new TH1F("ratio","",nbin,xmin,xmax);
Inte->Sumw2();
for (int i = 0; i< nbin-1; i++)
  {
    double sumPredi = 0;
    double sumData = 0;
    double Interatio = 0;
    for (int j = i+1 ; j < nbin ; j++)
      {
          sumPredi += hsolveG->GetBinContent(j);
          sumData += htotDataG->GetBinContent(j);  
      } 
      if (sumPredi == 0){Inte->AddBinContent(i,0);} 
      else {Inte->AddBinContent(i,sumData/sumPredi);}
 
  }
  Inte->Draw("PE1"); 
  //  hsolveC->SetFillStyle(3004);
//  Inte->SetFillColor(kBlack);
 Inte->SetLineColor(kBlack);
 Inte->SetLineStyle(1);
 Inte->SetLineWidth(1);
 Inte->SetTickLength(0.03, "YZ");
 Inte->SetTickLength(0.03,"X");
 Inte->SetLabelOffset(0.015,"X");
 Inte->SetLabelOffset(0.007,"Y");
 Inte->SetLabelSize(0.1, "XYZ");
 Inte->SetLabelFont(42, "XYZ"); 
 Inte->SetTitleSize(0.08, "XYZ"); 
 Inte->SetTitleFont(42, "XYZ");
 Inte->SetTitleOffset(1.2,"X"); 
 Inte->SetTitleOffset(0.5,"Y");
 Inte->GetXaxis()->SetTitle(xtitle);
 Inte->GetXaxis()->SetTitleColor(1);
 Inte->GetYaxis()->SetTitle("Data/Prediction Integral ratio");
 Inte->GetYaxis()->SetTitleColor(1);
 Inte->SetNdivisions(509,"XYZ");
 Inte->SetMaximum(3); 
 Inte->SetMinimum(0); 
 Inte->SetMarkerStyle(20);
 Inte->SetMarkerSize(1);
 
 c2->Update();


 padD1->cd();
    // hsolveD1->SetLineColor(kBlack);
    // hsolveD1->SetFillColor(kBlack);
    // hsolveD1->SetFillStyle(3004);
    // //  hsolve->SetLineStyle(1);
    // //  hsolve->SetLineWidth(2);
    // hsolveD1->Draw("E4same");

    hsolveD1->SetLineColor(kAzure+6);
    hsolveD1->SetFillColor(kAzure+6);
    hsolveD1->SetFillStyle(3001);
    hsolveD1->SetLineStyle(1);
    hsolveD1->SetLineWidth(1);
    hsolveD1->Draw("sameE2");

    padD1Ratio->cd();

TH1F* ratioD1  = new TH1F("ratioD1","",nbin,xmin,xmax);
ratioD1->Sumw2();
for (int i = 0; i< nbin; i++)
  {
    double data = htotDataD->GetBinContent(i);
    double pred = hsolveD1->GetBinContent(i);
    double sigma = hsolveD1->GetBinError(i);
    // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
    if (data == 0 || pred == 0) {ratioD1->AddBinContent(i,0);}    //when data == 0 and there is no prediction => you need this
    else
      {
        // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
        double pull = (data - pred) / sqrt(data + sigma*sigma);
        // double pull = data/pred;
        ratioD1->AddBinContent(i,pull);
      }
  }
 ratioD1->Draw("E4"); // PE1 ou E4
 ratioD1->SetFillColor(kRed);
 ratioD1->SetFillStyle(3004);
 ratioD1->SetLineColor(kBlack);
 ratioD1->SetLineStyle(1);
 ratioD1->SetLineWidth(1);
 ratioD1->SetTickLength(0.03, "YZ");
 ratioD1->SetTickLength(0.03,"X");
 ratioD1->SetLabelOffset(0.015,"X");
 ratioD1->SetLabelOffset(0.007,"Y");
 ratioD1->SetLabelSize(0.1, "XYZ");
 ratioD1->SetLabelFont(42, "XYZ"); 
 ratioD1->SetTitleSize(0.1, "XYZ"); 
 ratioD1->SetTitleFont(42, "XYZ");
 ratioD1->SetTitleOffset(1.2,"X"); 
 ratioD1->SetTitleOffset(0.5,"Y");
 ratioD1->GetXaxis()->SetTitle(xtitle);
 ratioD1->GetXaxis()->SetTitleColor(1);
 ratioD1->GetYaxis()->SetTitle("Pulls");
 ratioD1->GetYaxis()->SetTitleColor(1);
 ratioD1->SetNdivisions(509,"XYZ");
 ratioD1->SetMaximum(3); 
 ratioD1->SetMinimum(-3); 
 ratioD1->SetMarkerStyle(20);
 ratioD1->SetMarkerSize(1);


 padD2->cd();
    // hsolveD2->SetLineColor(kBlack);
    // hsolveD2->SetFillColor(kBlack);
    // hsolveD2->SetFillStyle(3004);
    // //  hsolve->SetLineStyle(1);
    // //  hsolve->SetLineWidth(2);
    // hsolveD2->Draw("E4same");

     hsolveD2->SetLineColor(kAzure+6);
    hsolveD2->SetFillColor(kAzure+6);
    hsolveD2->SetFillStyle(3001);
    hsolveD2->SetLineStyle(1);
    hsolveD2->SetLineWidth(1);
    hsolveD2->Draw("sameE2");

padD2Ratio->cd();

TH1F* ratioD2  = new TH1F("ratioD2","",nbin,xmin,xmax);
ratioD2->Sumw2();
for (int i = 0; i< nbin; i++)
  {
    double data = htotDataD->GetBinContent(i);
    double pred = hsolveD2->GetBinContent(i);
    double sigma = hsolveD2->GetBinError(i);
    // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
    if (data == 0 || pred == 0) {ratioD2->AddBinContent(i,0);}    //when data == 0 and there is no prediction => you need this
    else
      {
        // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
        double pull = (data - pred) / sqrt(data + sigma*sigma);
        // double pull = data/pred;
        ratioD2->AddBinContent(i,pull);
      }
  }
 ratioD2->Draw("E4"); // PE1 ou E4
 ratioD2->SetFillColor(kRed);
 ratioD2->SetFillStyle(3004);
 ratioD2->SetLineColor(kBlack);
 ratioD2->SetLineStyle(1);
 ratioD2->SetLineWidth(1);
 ratioD2->SetTickLength(0.03, "YZ");
 ratioD2->SetTickLength(0.03,"X");
 ratioD2->SetLabelOffset(0.015,"X");
 ratioD2->SetLabelOffset(0.007,"Y");
 ratioD2->SetLabelSize(0.1, "XYZ");
 ratioD2->SetLabelFont(42, "XYZ"); 
 ratioD2->SetTitleSize(0.1, "XYZ"); 
 ratioD2->SetTitleFont(42, "XYZ");
 ratioD2->SetTitleOffset(1.2,"X"); 
 ratioD2->SetTitleOffset(0.5,"Y");
 ratioD2->GetXaxis()->SetTitle(xtitle);
 ratioD2->GetXaxis()->SetTitleColor(1);
 ratioD2->GetYaxis()->SetTitle("Pulls");
 ratioD2->GetYaxis()->SetTitleColor(1);
 ratioD2->SetNdivisions(509,"XYZ");
 ratioD2->SetMaximum(3); 
 ratioD2->SetMinimum(-3); 
 ratioD2->SetMarkerStyle(20);
 ratioD2->SetMarkerSize(1);
 padH1->cd();
    // hsolveH1->SetLineColor(kBlack);
    // hsolveH1->SetFillColor(kBlack);
    // hsolveH1->SetFillStyle(3004);
    // //  hsolve->SetLineStyle(1);
    // //  hsolve->SetLineWidth(2);
    // hsolveH1->Draw("E4same");

     hsolveH1->SetLineColor(kAzure+6);
    hsolveH1->SetFillColor(kAzure+6);
    hsolveH1->SetFillStyle(3001);
    hsolveH1->SetLineStyle(1);
    hsolveH1->SetLineWidth(1);
    hsolveH1->Draw("sameE2");


padH1Ratio->cd();

TH1F* ratioH1  = new TH1F("ratioH1","",nbin,xmin,xmax);
ratioH1->Sumw2();
for (int i = 0; i< nbin; i++)
  {
    double data = htotDataH->GetBinContent(i);
    double pred = hsolveH1->GetBinContent(i);
    double sigma = hsolveH1->GetBinError(i);
    // std::cout<<"dataH = "<<data<<" predH = "<<pred<<" sigmaH = "<<sigma<<std::endl;
    if (data == 0 || pred == 0) {ratioH1->AddBinContent(i,0);}    //when data == 0 and there is no prediction => you need this
    else
      {
        // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
        double pull = (data - pred) / sqrt(data + sigma*sigma);
        // double pull = data/pred;
        ratioH1->AddBinContent(i,pull);
      }
  }
 ratioH1->Draw("E4"); // PE1 ou E4
 ratioH1->SetFillColor(kRed);
 ratioH1->SetFillStyle(3004);
 ratioH1->SetLineColor(kBlack);
 ratioH1->SetLineStyle(1);
 ratioH1->SetLineWidth(1);
 ratioH1->SetTickLength(0.03, "YZ");
 ratioH1->SetTickLength(0.03,"X");
 ratioH1->SetLabelOffset(0.015,"X");
 ratioH1->SetLabelOffset(0.007,"Y");
 ratioH1->SetLabelSize(0.1, "XYZ");
 ratioH1->SetLabelFont(42, "XYZ"); 
 ratioH1->SetTitleSize(0.1, "XYZ"); 
 ratioH1->SetTitleFont(42, "XYZ");
 ratioH1->SetTitleOffset(1.2,"X"); 
 ratioH1->SetTitleOffset(0.5,"Y");
 ratioH1->GetXaxis()->SetTitle(xtitle);
 ratioH1->GetXaxis()->SetTitleColor(1);
 ratioH1->GetYaxis()->SetTitle("Pulls");
 ratioH1->GetYaxis()->SetTitleColor(1);
 ratioH1->SetNdivisions(509,"XYZ");
 ratioH1->SetMaximum(3); 
 ratioH1->SetMinimum(-3); 
 ratioH1->SetMarkerStyle(20);
 ratioH1->SetMarkerSize(1);
 padH2->cd();
    // hsolveH2->SetLineColor(kBlack);
    // hsolveH2->SetFillColor(kBlack);
    // hsolveH2->SetFillStyle(3004);
    // //  hsolve->SetLineStyle(1);
    // //  hsolve->SetLineWidth(2);
    // hsolveH2->Draw("E4same");

    hsolveH2->SetLineColor(kAzure+6);
    hsolveH2->SetFillColor(kAzure+6);
    hsolveH2->SetFillStyle(3001);
    hsolveH2->SetLineStyle(1);
    hsolveH2->SetLineWidth(1);
    hsolveH2->Draw("sameE2");
padH2Ratio->cd();

TH1F* ratioH2  = new TH1F("ratioH2","",nbin,xmin,xmax);
ratioH2->Sumw2();
for (int i = 0; i< nbin; i++)
  {
    double data = htotDataH->GetBinContent(i);
    double pred = hsolveH2->GetBinContent(i);
    double sigma = hsolveH2->GetBinError(i);
    // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
    if (data == 0 || pred == 0) {ratioH2->AddBinContent(i,0);}    //when data == 0 and there is no prediction => you need this
    else
      {
        // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
        double pull = (data - pred) / sqrt(data + sigma*sigma);
        // double pull = data/pred;
        ratioH2->AddBinContent(i,pull);
      }
  }
 ratioH2->Draw("E4"); // PE1 ou E4
 ratioH2->SetFillColor(kRed);
 ratioH2->SetFillStyle(3004);
 ratioH2->SetLineColor(kBlack);
 ratioH2->SetLineStyle(1);
 ratioH2->SetLineWidth(1);
 ratioH2->SetTickLength(0.03, "YZ");
 ratioH2->SetTickLength(0.03,"X");
 ratioH2->SetLabelOffset(0.015,"X");
 ratioH2->SetLabelOffset(0.007,"Y");
 ratioH2->SetLabelSize(0.1, "XYZ");
 ratioH2->SetLabelFont(42, "XYZ"); 
 ratioH2->SetTitleSize(0.1, "XYZ"); 
 ratioH2->SetTitleFont(42, "XYZ");
 ratioH2->SetTitleOffset(1.2,"X"); 
 ratioH2->SetTitleOffset(0.5,"Y");
 ratioH2->GetXaxis()->SetTitle(xtitle);
 ratioH2->GetXaxis()->SetTitleColor(1);
 ratioH2->GetYaxis()->SetTitle("Pulls");
 ratioH2->GetYaxis()->SetTitleColor(1);
 ratioH2->SetNdivisions(509,"XYZ");
 ratioH2->SetMaximum(3); 
 ratioH2->SetMinimum(-3); 
 ratioH2->SetMarkerStyle(20);
 ratioH2->SetMarkerSize(1);

//********************************************µ//
  c1->cd();

  c2->SaveAs("./"+Name+"_"+Plots+"_v2.pdf");
  return c1;
}