#include <iostream>
#include <TROOT.h>
#include "TH1.h"
// #include "../MCWeights.h"

TCanvas * plot(int method, TString Prod, TString Name, TString Year, TString Dmode,TString Plots)
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
    hmaxBD = 1E6;
//     hmax = 1E6;   // for eta<2.4 pt>80  
//     hmaxBD = 1E9;
  }
  // TString Prod = "DATA_EMU_2017";
  //ABCD_EMu2018A 094
  // ABCD_2018A 095
  float ReScaleXS = 1.;
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


    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "A";
    TString HeaderB = "B";
    TString HeaderC = "C";
    TString HeaderD = " D";
    TString HeaderE = " E";
    TString HeaderF = " F";
    TString HeaderNVtx = "k Vtx";
    TString xtitle = "var";
  //-----------------------------------------------------------//
  // ABCD using Evt and Tight+looseWP 
  //-----------------------------------------------------------//

// //    //-----------------------------------------------------------//

// //    //-----------------------------------------------------------//
// //    // ABCD using Hemisphere pt and Tight+looseWP
// //    //-----------------------------------------------------------//
int Method = method;
  if (Method == 0)
  {
    htitleA = "hData_CRlowlowpt_1Vtx_Mass_";//
    htitleB = "hData_CRlooselowlowpt_1Vtx_Mass_";//

    htitleC = "hData_CRlowpt_1Vtx_Mass_";//
    htitleD = "hData_CRlooselowpt_1Vtx_Mass_";//

    htitleE = "hData_Hemi_1Vtx_Mass_";//*2 mais normal parce que blinding
    htitleF = "hData_CRloose_1Vtx_Mass_";//

    nbin = 25; 
    xmin = 0;
    xmax =  100;
  HeaderA = "T + 20<pt<80";
    HeaderB = "L + 20<pt<80";

    HeaderC = "T + 20<pt_{i}<80 & pt_{j}>80";
    HeaderD = "L + 20<pt_{i}<80 & pt_{j}>80";

    HeaderE = "T + pt>80";
    HeaderF = "L + pt>80";
    HeaderNVtx = "1 Vtx"; 
     xtitle = "Vtx Mass";
  }

if (Method == 1)
  {
    htitleA = "hData_CRlowlowpt_1Vtx_SumtrackWeight_";//
    htitleB = "hData_CRlooselowlowpt_1Vtx_SumtrackWeight_";//

    htitleC = "hData_CRlowpt_1Vtx_SumtrackWeight_";//
    htitleD = "hData_CRlooselowpt_1Vtx_SumtrackWeight_";//

    htitleE = "hData_Hemi_1Vtx_SumtrackWeight_";//*2 mais normal parce que blinding
    htitleF = "hData_CRloose_1Vtx_SumtrackWeight_";//

    nbin = 40; 
    xmin = 0;
    xmax = 40;
  HeaderA = "T + 20<pt<80";
    HeaderB = "L + 20<pt<80";

    HeaderC = "T + 20<pt_{i}<80 & pt_{j}>80";
    HeaderD = "L + 20<pt_{i}<80 & pt_{j}>80";

    HeaderE = "T + pt>80";
    HeaderF = "L + pt>80";

    HeaderNVtx = "1 Vtx";
    xtitle = "SumtrackWeight";
  }



  TLegend* leg;
    
// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

TPad* padA = new TPad("padA","This is padA",0.01,0.65,0.32,0.95,21);
TPad* padAratio = new TPad("padAratio","This is padA",0.01,0.50,0.32,0.65,21);
TPad* padB = new TPad("padB","This is padB",0.01,0.15,0.32,0.45,21);
TPad* padBratio = new TPad("padBratio","This is padB",0.01,0.01,0.32,0.15,21);

TPad* padC = new TPad("padC","This is padC",0.33,0.65,0.65,0.95,21);
TPad* padCratio = new TPad("padCratio","This is padC",0.33,0.50,0.65,0.65,21);
TPad* padD = new TPad("padD","This is padD",0.33,0.15,0.65,0.45,21);
TPad* padDratio = new TPad("padDratio","This is padD",0.33,0.01,0.65,0.15,21);

TPad* padE = new TPad("padE","This is padE",0.66,0.65,0.99,0.95,21);
TPad* padEratio = new TPad("padEratio","This is padE",0.66,0.50,0.99,0.65,21);
TPad* padF = new TPad("padF","This is padF",0.66,0.15,0.99,0.45,21);
TPad* padFratio = new TPad("padFratio","This is padF",0.66,0.01,0.99,0.15,21);
// 
  // //ratio padds
  // TPad* pad8 = new TPad("pad8","This is pad8",0.04,0.5,0.48,0.65,21);
  // TPad* pad9 = new TPad("pad9","This is pad9",0.52,0.5,0.96,0.65,21);
  // TPad* padA0 = new TPad("padA0","This is padA0",0.04,0.05,0.48,0.15,21);
  // TPad* padA1 = new TPad("padA1","This is padA1",0.52,0.05,0.96,0.15,21);

padA->SetFillColor(0);
padA->SetBorderMode(0);
padA->SetFrameFillColor(10);
padA->Draw();
padA->SetLogy(logy);
   padA->SetTopMargin(0.07);
   padA->SetBottomMargin(0.13);
   padA->SetRightMargin(0.04);
   padA->SetLeftMargin(0.16);

padAratio->SetFillColor(0);
padAratio->SetBorderMode(0);
padAratio->SetFrameFillColor(10);
padAratio->Draw();
padAratio->SetLogy(0);
   padAratio->SetTopMargin(0.07);
   padAratio->SetBottomMargin(0.13);
   padAratio->SetRightMargin(0.04);
   padAratio->SetLeftMargin(0.16);

padC->SetFillColor(0);
padC->SetBorderMode(0);
padC->SetFrameFillColor(10);
padC->Draw();
padC->SetLogy(logy);
   padC->SetTopMargin(0.07);
   padC->SetBottomMargin(0.13);
   padC->SetRightMargin(0.04);
   padC->SetLeftMargin(0.16);

   padCratio->SetFillColor(0);
padCratio->SetBorderMode(0);
padCratio->SetFrameFillColor(10);
padCratio->Draw();
padCratio->SetLogy(0);
   padCratio->SetTopMargin(0.07);
   padCratio->SetBottomMargin(0.13);
   padCratio->SetRightMargin(0.04);
   padCratio->SetLeftMargin(0.16);

padB->SetFillColor(0);
padB->SetBorderMode(0);
padB->SetFrameFillColor(10);
padB->Draw();
padB->SetLogy(logy);
   padB->SetTopMargin(0.07);
   padB->SetBottomMargin(0.13);
   padB->SetRightMargin(0.04);
   padB->SetLeftMargin(0.16);

padBratio->SetFillColor(0);
padBratio->SetBorderMode(0);
padBratio->SetFrameFillColor(10);
padBratio->Draw();
padBratio->SetLogy(0);
   padBratio->SetTopMargin(0.07);
   padBratio->SetBottomMargin(0.13);
   padBratio->SetRightMargin(0.04);
   padBratio->SetLeftMargin(0.16);


padD->SetFillColor(0);
padD->SetBorderMode(0);
padD->SetFrameFillColor(10);
padD->Draw();
padD->SetLogy(logy);
   padD->SetTopMargin(0.07);
   padD->SetBottomMargin(0.13);
   padD->SetRightMargin(0.04);
   padD->SetLeftMargin(0.16);

padDratio->SetFillColor(0);
padDratio->SetBorderMode(0);
padDratio->SetFrameFillColor(10);
padDratio->Draw();
padDratio->SetLogy(0);
   padDratio->SetTopMargin(0.07);
   padDratio->SetBottomMargin(0.13);
   padDratio->SetRightMargin(0.04);
   padDratio->SetLeftMargin(0.16);

padE->SetFillColor(0);
padE->SetBorderMode(0);
padE->SetFrameFillColor(10);
padE->Draw();
padE->SetLogy(logy);
   padE->SetTopMargin(0.07);
   padE->SetBottomMargin(0.13);
   padE->SetRightMargin(0.04);
   padE->SetLeftMargin(0.16);

padEratio->SetFillColor(0);
padEratio->SetBorderMode(0);
padEratio->SetFrameFillColor(10);
padEratio->Draw();
padEratio->SetLogy(0);
   padEratio->SetTopMargin(0.07);
   padEratio->SetBottomMargin(0.13);
   padEratio->SetRightMargin(0.04);
   padEratio->SetLeftMargin(0.16);

padF->SetFillColor(0);
padF->SetBorderMode(0);
padF->SetFrameFillColor(10);
padF->Draw();
padF->SetLogy(logy);
   padF->SetTopMargin(0.07);
   padF->SetBottomMargin(0.13);
   padF->SetRightMargin(0.04);
   padF->SetLeftMargin(0.16);

   padFratio->SetFillColor(0);
padFratio->SetBorderMode(0);
padFratio->SetFrameFillColor(10);
padFratio->Draw();
padFratio->SetLogy(0);
   padFratio->SetTopMargin(0.07);
   padFratio->SetBottomMargin(0.13);
   padFratio->SetRightMargin(0.04);
   padFratio->SetLeftMargin(0.16);

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

 TH1F* hsolveC  = new TH1F("hsolveC","",nbin,xmin,xmax);
 TH1F* hsolveE  = new TH1F("hsolveE","",nbin,xmin,xmax);
 hsolveC->Sumw2();
 hsolveE->Sumw2();

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

//////// !! ----------------------------!! //////////////
 padA->cd();
//////// !! ----------------------------!! //////////////

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
  TH1F* g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleA);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//   TH1F*g2_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
//  htotData->Add(g2_Data_emu, htotData, 1, 1);

//  f3_Data_emu->cd();
//   TH1F*g3_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
//  htotData->Add(g3_Data_emu, htotData, 1, 1);

//  f4_Data_emu->cd();
//   TH1F*g4_Data_emu = (TH1F*)gROOT->FindObject(htitleA);
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
 h_VV->Add(g2_VV, h_VV,1, 1);

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



 hsolveC->Add(hsolveC, htotData, 0., 1.);
 hsolveE->Add(hsolveE, htotData, 0., 1.);

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
  //   leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV "+CTAU[1]+"mm","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV "+CTAU[2]+"mm","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
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
  leg->SetHeader(HeaderA);
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 40-80 GeV");
  leg->Draw();

     leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("A");
  leg->Draw();
  //Data/MC ---------------------
  padAratio->cd();
  
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

 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//


// *****************************************************************************
  c1->cd();
  //////// !! ----------------------------!! //////////////
 padB->cd();
 //////// !! ----------------------------!! //////////////

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
 hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
 hmax = hmaxBD; 

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
   hDataMC->GetYaxis()->SetTitle("Data/MC");
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

 hsolveC->Divide(hsolveC, htotData, 1., 1.);
 hsolveE->Divide(hsolveE, htotData, 1., 1.);

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
  leg->SetHeader(HeaderB);
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 40-80 GeV");
  leg->Draw();

     leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("B");
  leg->Draw();
  padBratio->cd();
  
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
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
//////// !! ----------------------------!! //////////////
 padD->cd();
 //////// !! ----------------------------!! //////////////
 hmax = hmaxBD; 

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

 hsolveC->Multiply(hsolveC, htotData, 1., 1.);

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
  padDratio->cd();
  
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
//////// !! ----------------------------!! //////////////
padF->cd();
//////// !! ----------------------------!! //////////////
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

hsolveE->Multiply(hsolveE, htotData, 1., 1.);

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
  leg->SetHeader(HeaderF);
  leg->Draw();

leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("F");
  leg->Draw();

    padFratio->cd();
  
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
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
  //////// !! ----------------------------!! //////////////
 padC->cd();
//////// !! ----------------------------!! //////////////


 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
 hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

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

TH1F * htotDataC = new TH1F("htotData","",nbin,xmin,xmax);
  htotDataC->Add(htotData, htotDataC, 1., 0.);


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

//      hsolveC->SetLineColor(kBlack);
// //  hsolveC->SetFillStyle(3004);
//  hsolveC->SetLineStyle(1);
//  hsolveC->SetLineWidth(2);


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
  leg->Draw();

leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("C");
  leg->Draw();
  padCratio->cd();
  
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
  ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//


  TCanvas *c2 = new TCanvas("c2", "plots",200,0,700,700);
  c2->SetFillColor(10);
  c2->SetFillStyle(4000);
  c2->SetBorderSize(2);
  TPad* pad5 = new TPad("pad5","This is pad5",0.51,0.55,0.99,0.96,21);
  TPad* pad6 = new TPad("pad6","This is pad6",0.51,0.35,0.99,0.55,21);
  TPad* pad7 = new TPad("pad7","This is pad7",0.51,0.1,0.99,0.35,21);

  TPad* pC = new TPad("pC","This is pC",0.01,0.55,0.5,0.96,21);
  TPad* pCRatio = new TPad("pCRatio","This is pCRatio",0.01,0.35,0.5,0.55,21);
TPad* pCInte = new TPad("pCInte","This is pCInte",0.01,0.1,0.5,0.35,21);

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

pC->SetFillColor(0);
pC->SetBorderMode(0);
pC->SetFrameFillColor(10);
pC->Draw();
pC->SetLogy(logy);
   pC->SetTopMargin(0.07);
   pC->SetBottomMargin(0.13);
   pC->SetRightMargin(0.04);
   pC->SetLeftMargin(0.16);

pCRatio->SetFillColor(0);
pCRatio->SetBorderMode(0);
pCRatio->SetFrameFillColor(10);
pCRatio->Draw();
pCRatio->SetLogy(0);
   pCRatio->SetTopMargin(0.07);
   pCRatio->SetBottomMargin(0.13);
   pCRatio->SetRightMargin(0.04);
   pCRatio->SetLeftMargin(0.16);

pCInte->SetFillColor(0);
pCInte->SetBorderMode(0);
pCInte->SetFrameFillColor(10);
pCInte->Draw();
pCInte->SetLogy(0);
   pCInte->SetTopMargin(0.07);
   pCInte->SetBottomMargin(0.13);
   pCInte->SetRightMargin(0.04);
   pCInte->SetLeftMargin(0.16);


  c2->cd();
  pC->cd();
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

//  h_VV->Draw("HEsame"); 
//  h_VV->SetFillColor(kOrange-2);
//  h_VV->SetLineColor(kOrange-2);
//  h_VV->SetLineStyle(1);
//  h_VV->SetLineWidth(3);

//  h_TTV->Draw("HEsame"); 
//  h_TTV->SetFillColor(kAzure+4);
//  h_TTV->SetLineColor(kAzure+4);
//  h_TTV->SetLineStyle(1);
//  h_TTV->SetLineWidth(3);

//  h_ST->Draw("HEsame"); 
//  h_ST->SetFillColor(kAzure+2);
//  h_ST->SetLineColor(kAzure+2);
//  h_ST->SetLineStyle(1);
//  h_ST->SetLineWidth(3);

//  h_TT->Draw("HEsame"); 
//  h_TT->SetFillColor(kAzure+1);
//  h_TT->SetLineColor(kAzure+1);
//  h_TT->SetLineStyle(1);
//  h_TT->SetLineWidth(3);
//  h_TT->SetTickLength(0.03, "YZ");
//  h_TT->SetTickLength(0.03,"X");

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
 
  // leg = new TLegend(0.64,0.50,0.85,0.89);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.035);
  // leg->SetMargin(0.2);
  // leg->AddEntry(htotData, " e#mu data","PE1");
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu}}  = 200  GeV","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu}}  = 300  GeV","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu}}  = 400  GeV","L");
  //   leg->AddEntry(htotMC, " Total MC","F");
  // // leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  // // leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  // // leg->AddEntry(h_ST, " tW","F");
  // // leg->AddEntry(h_TT, " t#bar{t}","F");
  // // leg->AddEntry(hsolveE," Prediction","LPE");
  // leg->AddEntry(hsolveE," Prediction","FE4");
  // //   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  // leg->Draw();
c1->cd();
//////// !! ----------------------------!! //////////////
padE->cd();
//////// !! ----------------------------!! //////////////
 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleE);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//  f2_Data_emu->cd();
//  g2_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
//  htotData->Add(g2_Data_emu, htotData, 0, 0);

//  f3_Data_emu->cd();
//  g3_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
//  htotData->Add(g3_Data_emu, htotData, 0, 0);

//  f4_Data_emu->cd();
//  g4_Data_emu = (TH1F*)gROOT->FindObject(htitleE);
//  htotData->Add(g4_Data_emu, htotData, 0, 0);



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
   padEratio->cd();
  
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
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
    leg->AddEntry(htotMC, " Total MC","F");
  // leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  // leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  // leg->AddEntry(h_ST, " tW","F");
  // leg->AddEntry(h_TT, " t#bar{t}","F");
  // leg->AddEntry(hsolveE," Prediction","LPE");
leg->AddEntry(hsolveE," Prediction","FE4");
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
pad6->cd();
//pulls between the predictions and data using the following formula:
//pull = (data - prediction) / sqrt(data+sigma^2_prediction)

// Poissionan errors are added due to low stats

 
TH1F* ratio  = new TH1F("ratio","",nbin,xmin,xmax);
ratio->Sumw2();
for (int i = 0; i< nbin; i++)
  {
    double data = htotData->GetBinContent(i);
    double pred = hsolveE->GetBinContent(i);
    double sigma = hsolveE->GetBinError(i);
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
 ratio->SetMinimum(-3); 
 ratio->SetMaximum(3); 
 ratio->SetMarkerStyle(20);
 ratio->SetMarkerSize(1);

pad7->cd();
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
          sumPredi += hsolveE->GetBinContent(j);
          sumData += htotData->GetBinContent(j);  
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
 Inte->SetLabelOffset(0.010,"Y");
 Inte->SetLabelSize(0.08, "XYZ");
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
 Inte->SetMinimum(0); 
 Inte->SetMaximum(5); 
 Inte->SetMarkerStyle(20);
 Inte->SetMarkerSize(1);
 
//  padC->cd();
pC->cd();
  htotDataC->Draw("E1same"); 
 htotDataC->SetFillColor(kBlack);
 htotDataC->SetLineColor(kBlack);
 htotDataC->SetLineStyle(1);
 htotDataC->SetLineWidth(1);
 htotDataC->SetTickLength(0.03, "YZ");
 htotDataC->SetTickLength(0.03,"X");
 htotDataC->SetLabelOffset(0.015,"X");
 htotDataC->SetLabelOffset(0.007,"Y");
 htotDataC->SetLabelSize(0.045, "XYZ");
 htotDataC->SetLabelFont(42, "XYZ"); 
 htotDataC->SetTitleSize(0.045, "XYZ"); 
 htotDataC->SetTitleFont(42, "XYZ");
 htotDataC->SetTitleOffset(1.2,"X"); 
 htotDataC->SetTitleOffset(1.3,"Y");
 htotDataC->GetXaxis()->SetTitle(xtitle);
 htotDataC->GetXaxis()->SetTitleColor(1);
 htotDataC->GetYaxis()->SetTitle(ytitle);
 htotDataC->GetYaxis()->SetTitleColor(1);
 htotDataC->SetNdivisions(509,"XYZ");
 htotDataC->SetMinimum(hmin); 
 htotDataC->SetMaximum(hmax); 
 htotDataC->SetMarkerStyle(20);
 htotDataC->SetMarkerSize(1);

   hsolveC->SetLineColor(kGray);
 hsolveC->SetFillColor(kGray);
 hsolveC->SetFillStyle(3001);
 hsolveC->SetLineStyle(1);
 hsolveC->SetLineWidth(1);
 hsolveC->Draw("sameE2"); 


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
  leg->AddEntry(htotDataC, " e#mu data","PE1");
    leg->AddEntry(htotMC, "Total MC","F");
  // leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  // leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  // leg->AddEntry(h_ST, " tW","F");
  // leg->AddEntry(h_TT, " t#bar{t}","F");
  leg->AddEntry(hsolveC," Prediction","FE4");
  //   leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV "+CTAU[0]+"mm","L");
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
  leg->Draw();

     leg = new TLegend(0.50,0.50,0.55,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("C");
  leg->Draw();

pCRatio->cd();

TH1F* ratioC  = new TH1F("ratioC","",nbin,xmin,xmax);
ratioC->Sumw2();
for (int i = 0; i< nbin; i++)
  {
    double data = htotDataC->GetBinContent(i);
    double pred = hsolveC->GetBinContent(i);
    double sigma = hsolveC->GetBinError(i);
    // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
    if (data == 0 || pred == 0) {ratioC->AddBinContent(i,0);}    //when data == 0 and there is no prediction => you need this
    else
      {
        // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
        double pull = (data - pred) / sqrt(data + sigma*sigma);
        ratioC->AddBinContent(i,pull);
      }
  }
 ratioC->Draw("E4"); // PE1 ou E4
 ratioC->SetFillColor(kRed);
 ratioC->SetFillStyle(3004);
 ratioC->SetLineColor(kBlack);
 ratioC->SetLineStyle(1);
 ratioC->SetLineWidth(1);
 ratioC->SetTickLength(0.03, "YZ");
 ratioC->SetTickLength(0.03,"X");
 ratioC->SetLabelOffset(0.015,"X");
 ratioC->SetLabelOffset(0.007,"Y");
 ratioC->SetLabelSize(0.1, "XYZ");
 ratioC->SetLabelFont(42, "XYZ"); 
 ratioC->SetTitleSize(0.1, "XYZ"); 
 ratioC->SetTitleFont(42, "XYZ");
 ratioC->SetTitleOffset(1.2,"X"); 
 ratioC->SetTitleOffset(0.5,"Y");
 ratioC->GetXaxis()->SetTitle(xtitle);
 ratioC->GetXaxis()->SetTitleColor(1);
 ratioC->GetYaxis()->SetTitle("Pulls");
 ratioC->GetYaxis()->SetTitleColor(1);
 ratioC->SetNdivisions(509,"XYZ");
 ratioC->SetMinimum(-3); 
 ratioC->SetMaximum(3); 
 ratioC->SetMarkerStyle(20);
 ratioC->SetMarkerSize(1);


pCInte->cd();
// // intégrale à droite
TH1F* InteC  = new TH1F("ratio","",nbin,xmin,xmax);
InteC->Sumw2();
for (int i = 0; i< nbin-1; i++)
  {
    double sumPredi = 0;
    double sumData = 0;
    double Interatio = 0;
    for (int j = i+1 ; j < nbin ; j++)
      {
          sumPredi += hsolveC->GetBinContent(j);
          sumData += htotDataC->GetBinContent(j);  
      } 
      if (sumPredi == 0){InteC->AddBinContent(i,0);} 
      else {InteC->AddBinContent(i,sumData/sumPredi);}
 
  }
  InteC->Draw("PE1"); 
  //  hsolveC->SetFillStyle(3004);
//  InteC->SetFillColor(kBlack);
 InteC->SetLineColor(kBlack);
 InteC->SetLineStyle(1);
 InteC->SetLineWidth(1);
 InteC->SetTickLength(0.03, "YZ");
 InteC->SetTickLength(0.03,"X");
 InteC->SetLabelOffset(0.015,"X");
 InteC->SetLabelOffset(0.010,"Y");
 InteC->SetLabelSize(0.08, "XYZ");
 InteC->SetLabelFont(42, "XYZ"); 
 InteC->SetTitleSize(0.08, "XYZ"); 
 InteC->SetTitleFont(42, "XYZ");
 InteC->SetTitleOffset(1.2,"X"); 
 InteC->SetTitleOffset(0.5,"Y");
 InteC->GetXaxis()->SetTitle(xtitle);
 InteC->GetXaxis()->SetTitleColor(1);
 InteC->GetYaxis()->SetTitle("Data/Prediction InteCgral ratio");
 InteC->GetYaxis()->SetTitleColor(1);
 InteC->SetNdivisions(509,"XYZ");
 InteC->SetMinimum(0); 
 InteC->SetMaximum(5); 
 InteC->SetMarkerStyle(20);
 InteC->SetMarkerSize(1);
// *****************************************************************************

  c1->Update();
  c2->cd();
  pad5->cd();

   //    hsolveE->SetLineColor(kBlack);
 //  hsolveE->SetFillColor(kBlack);
 //  hsolveE->SetFillStyle(3004);
// //  hsolveC->SetLineStyle(1);
// //  hsolveC->SetLineWidth(2);
  //   hsolveE->Draw("E4same");

   hsolveE->SetLineColor(kGray);
 hsolveE->SetFillColor(kGray);
 hsolveE->SetFillStyle(3001);
 hsolveE->SetLineStyle(1);
 hsolveE->SetLineWidth(1);
 hsolveE->Draw("sameE2"); 

  c2->SaveAs("./"+Name+"_"+Plots+"_v2.pdf");

  return c1;
}