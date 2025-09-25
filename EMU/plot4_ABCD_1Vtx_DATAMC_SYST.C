#include <iostream>
#include <TROOT.h>
#include "TH1.h"
#include "TColor.h"
// #include "../MCWeights.h"

TCanvas * plot(int method, TString Prod, TString Name, TString Year, TString Dmode,  TString Plots, bool Data, bool Mc, bool SIGNAL)
{
int stati=0;
bool fit= 1;
bool logy=1;

bool DATA = Data;
bool MC = Mc;
bool Signal = SIGNAL;

  float hmin = 0.5; // cannot be 0 for logy=1
  float hmax = 1E6;	         // for eta<2.4 pt>80 or CRlowlowpt
  float hmaxBD = 1E6;	         // for eta<2.4 pt>80 or CRlowlowpt
  float ReScaleXS = 1.;

TString Yearcor = Year;
 if (Year == "2016PRE") Yearcor = "2016preVFP";
 if (Year == "2016POST") Yearcor = "2016";
// Dmu
TFile* f1_Data_emu  = new TFile("../../DATA_EMU_"+Year+"_03_02_2025/histofile_"+Dmode+"_OS_2p4_MuonEG-Run2018-UL2018_MiniAODv2_GT36-v1.root");
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
 TFile* f1_DY  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f2_DY  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f1_TT  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_TT  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_ST  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_ST  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f3_ST  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f4_ST  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_TTV = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8.root");
 TFile* f2_TTV = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f3_TTV = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8.root");
 TFile* f1_VV  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_VV  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 TFile* f3_VV  = new TFile("../../MC_EMU_03_02_2025/histofile_"+Dmode+"_OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 

TString MSUON[3] = {"200","300","300"};
  TString MNEU[3] = {"180","180","200"};
  TString CTAU[3] = {"100","100","100"};

//signal
 TFile* f1_LLP = new TFile("../../Signal_"+Year+"_L1/histofile_DM_OS_2p4_RPV_"+Yearcor+"_smu"+MSUON[0]+"_neu"+MNEU[0]+"_ctau"+CTAU[0]+".root");
 TFile* f2_LLP = new TFile("../../Signal_"+Year+"_L1/histofile_DM_OS_2p4_RPV_"+Yearcor+"_smu"+MSUON[1]+"_neu"+MNEU[1]+"_ctau"+CTAU[1]+".root");
 TFile* f3_LLP = new TFile("../../Signal_"+Year+"_L1/histofile_DM_OS_2p4_RPV_"+Yearcor+"_smu"+MSUON[2]+"_neu"+MNEU[2]+"_ctau"+CTAU[2]+".root");

//--SYST--//

std::vector<TFile*> f_SYST_Up;
 std::vector<TFile*> f_SYST_Down;

  bool WLepton = true;
  bool WGen = false;
  bool Signal = true;
  std::vector<TString> SYSTNAME_UP = {"LumiUp","L1Up","TriggerUp","PUUp","TopPtUp","JECUp","JERUp"}; //,,"SFEleUp",

  std::vector<TString> SYSTNAME_DOWN = {"LumiDown","L1Down","TriggerDown","PUDown","TopPtDown","JECDown","JERDown"};//,"SFEleDown",,"RoccorDown",

  if (WGen)
    {
      SYSTNAME_UP.push_back("PDFUp");
      SYSTNAME_UP.push_back("ScaleUp");
      SYSTNAME_DOWN.push_back("PDFDown");
      SYSTNAME_DOWN.push_back("ScaleDown");
    }
  if (WLepton)
    {
      
      // SYSTNAME_UP.push_back("RoccorUp");
      SYSTNAME_UP.push_back("MuonIDUp");
      SYSTNAME_UP.push_back("MuonISOUp");
      SYSTNAME_UP.push_back("EleIDUp");
      SYSTNAME_UP.push_back("EleISOUp");
      // SYSTNAME_DOWN.push_back("RoccorDown");
      SYSTNAME_DOWN.push_back("MuonIDDown");
      SYSTNAME_DOWN.push_back("MuonISODown");
      SYSTNAME_DOWN.push_back("EleIDDown");
      SYSTNAME_DOWN.push_back("EleISODown");
    }

  TFile* f_LumiUp  = new TFile("./SYST/LumiUp_SYST.root");
  TFile* f_L1Up  = new TFile("./SYST/L1Up_SYST.root");
  TFile* f_TriggerUp  = new TFile("./SYST/TriggerUp_SYST.root");
  TFile* f_MuonIDUp = new TFile("./SYST/MuonIDUp_SYST.root");
  TFile* f_MuonISOUp = new TFile("./SYST/MuonISOUp_SYST.root");
  TFile* f_EleIDUp = new TFile("./SYST/EleIDUp_SYST.root");
  TFile* f_EleISOUp = new TFile("./SYST/EleISOUp_SYST.root");
  TFile* f_PUUp  = new TFile("./SYST/PUUp_SYST.root");
  TFile* f_SFEleUp  = new TFile("./SYST/SFEleUp_SYST.root");
  TFile* f_TopPtUp  = new TFile("./SYST/TopPtUp_SYST.root");
  TFile* f_PDFUp  = new TFile("./SYST/PDFUp_SYST.root");
  TFile* f_ScaleUp  = new TFile("./SYST/ScaleUp_SYST.root");
    TFile* f_JECUp = new TFile("./SYST/JECUp_SYST.root");
  TFile* f_JERUp  = new TFile("./SYST/JERUp_SYST.root");
  TFile* f_RoccorUp  = new TFile("./SYST/RoccorUp_SYST.root");

  TFile* f_LumiDown  = new TFile("./SYST/LumiDown_SYST.root");
  TFile* f_L1Down  = new TFile("./SYST/L1Down_SYST.root");
  TFile* f_TriggerDown = new TFile("./SYST/TriggerDown_SYST.root");
  TFile* f_MuonIDDown = new TFile("./SYST/MuonIDDown_SYST.root");
  TFile* f_MuonISODown  = new TFile("./SYST/MuonISODown_SYST.root");
  TFile* f_EleIDDown = new TFile("./SYST/EleIDDown_SYST.root");
  TFile* f_EleISODown  = new TFile("./SYST/EleISODown_SYST.root");
  TFile* f_PUDown  = new TFile("./SYST/PUDown_SYST.root");

  TFile* f_SFEleDown  = new TFile("./SYST/SFEleDown_SYST.root");
  TFile* f_TopPtDown  = new TFile("./SYST/TopPtDown_SYST.root");
  TFile* f_PDFDown  = new TFile("./SYST/PDFDown_SYST.root");
  TFile* f_ScaleDown  = new TFile("./SYST/ScaleDown_SYST.root");
      TFile* f_JECDown = new TFile("./SYST/JECDown_SYST.root");
  TFile* f_JERDown  = new TFile("./SYST/JERDown_SYST.root");
  TFile* f_RoccorDown  = new TFile("./SYST/RoccorDown_SYST.root");


f_SYST_Up.push_back(f_LumiUp);
f_SYST_Up.push_back(f_L1Up);
f_SYST_Up.push_back(f_TriggerUp);
f_SYST_Up.push_back(f_PUUp);
// f_SYST_Up.push_back(f_SFEleUp);
f_SYST_Up.push_back(f_TopPtUp);

f_SYST_Up.push_back(f_JECUp);
f_SYST_Up.push_back(f_JERUp);
// f_SYST_Up.push_back(f_RoccorUp);


if (WGen)
  {
    f_SYST_Up.push_back(f_PDFUp);
    f_SYST_Up.push_back(f_ScaleUp);
  }
if (WLepton)
  {
    // f_SYST_Up.push_back(f_RoccorUp);
    f_SYST_Up.push_back(f_MuonIDUp);
    f_SYST_Up.push_back(f_MuonISOUp);
    f_SYST_Up.push_back(f_EleIDUp);
    f_SYST_Up.push_back(f_EleISOUp);
  }

f_SYST_Down.push_back(f_LumiDown);
f_SYST_Down.push_back(f_L1Down);
f_SYST_Down.push_back(f_TriggerDown);
f_SYST_Down.push_back(f_PUDown);
// f_SYST_Down.push_back(f_SFEleDown);
f_SYST_Down.push_back(f_TopPtDown);

f_SYST_Down.push_back(f_JECDown);
f_SYST_Down.push_back(f_JERDown);
// f_SYST_Down.push_back(f_RoccorDown);
// f_SYST_Down.push_back(f_PDFDown);
// f_SYST_Down.push_back(f_ScaleDown);
if (WGen)
  {
    f_SYST_Down.push_back(f_PDFDown);
    f_SYST_Down.push_back(f_ScaleDown);
  }
if (WLepton)
  {
    // f_SYST_Down.push_back(f_RoccorDown);
    f_SYST_Down.push_back(f_MuonIDDown);
    f_SYST_Down.push_back(f_MuonISODown);
    f_SYST_Down.push_back(f_EleIDDown);
    f_SYST_Down.push_back(f_EleISODown);
  }


//------//






 TString DATAFILE[1] = {"MuonEG-Run2018-UL2018_MiniAODv2_GT36-v1_"
};
if (Year == "2018") DATAFILE[0] = "MuonEG-Run2018-UL2018_MiniAODv2_GT36-v1_";


TString MCFILE[14] = {
  
                  "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_",
                  "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_",
                  "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_",
                  "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_",
                  "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_",
                  "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_",
                  "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8_",
                  "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8_",
                  "TTWW_TuneCP5_13TeV-madgraph-pythia8_",
                  "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_",
                  "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_",
                  "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_",
                  "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_",
                  "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_"
};


TString LLPFILE[3] = {"RPV_"+Yearcor+"_smu"+MSUON[0]+"_neu"+MNEU[0]+"_ctau"+CTAU[0]+"_",
                  "RPV_"+Yearcor+"_smu"+MSUON[1]+"_neu"+MNEU[1]+"_ctau"+CTAU[1]+"_",
                  "RPV_"+Yearcor+"_smu"+MSUON[2]+"_neu"+MNEU[2]+"_ctau"+CTAU[2]+"_"
};

 TString ytitle = "Events"; 
 TString HeaderCMS = "CMS";

if (DATA ) 
  {
      if (Year == "2016") HeaderCMS = "2016                       36.3 fb^{-1} (13 TeV)";
      if (Year == "2017") HeaderCMS = "2017                       41.5 fb^{-1} (13 TeV)";
      if (Year == "2018") HeaderCMS = "2018                       59.8 fb^{-1} (13 TeV)";
  }
if (MC)
  {
      if (Year == "2016") HeaderCMS = "Simulation 2016        36.3 fb^{-1} (13 TeV)";
      if (Year == "2017") HeaderCMS = "Simulation 2017        41.5 fb^{-1} (13 TeV)";
      if (Year == "2018") HeaderCMS = "Simulation 2018        59.8 fb^{-1} (13 TeV)";
  }
if (DATA && MC)
  {
      if (Year == "2016") HeaderCMS = "2016                                                 36.3 fb^{-1} (13 TeV)";
      if (Year == "2017") HeaderCMS = "2017                                                 41.5 fb^{-1} (13 TeV)";
      if (Year == "2018") HeaderCMS = "2018                                                 59.8 fb^{-1} (13 TeV)";
  }
    TString htitleA = "hData_EVT34_1Vtx_BDTvtx";
    TString htitleB = "hData_NoEVT34_1Vtx_BDTvtx";
    TString htitleC = "hData_EVT12_1Vtx_BDTvtx";
    TString htitleD = "hData_NoEVT12_1Vtx_BDTvtx";
    int nbin = 40; 
    float xmin = 0;
    float xmax =  40;
    TString HeaderA = "A";
    TString HeaderB = "B";
    TString HeaderC = "C";
    TString HeaderD = "D";
    TString HeaderNVtx = "k Vtx";
    TString xtitle = "var";
    float rwTT = 1;//0.00923001376;//0.923001376;
  //-----------------------------------------------------------//
  // ABCD using Hemipt and Tight+looseWP 
  //-----------------------------------------------------------//
  int Method = method;
  if (Method == 0)
  {
    htitleA = "hData_CRlowlowpt_1Vtx_Mass_";
    htitleB = "hData_CRlooselowlowpt_1Vtx_Mass_";
    htitleC = "hData_CRTight_1Vtx_Mass_";
    htitleD = "hData_CRLoose_1Vtx_Mass_";

    nbin = 25; 
    xmin = 0;
    xmax =  100;
    HeaderA = "T + 30<pt_{i}<80 & pt_{j}>80";
    HeaderB = "L + 30<pt<80";
    HeaderC = "T + 30<pt_{i}<80 & pt_{j}>80";
    HeaderD = "L + 30<pt<80";

    HeaderNVtx = "1 Vtx"; 
     xtitle = "Vtx Mass [GeV]";
  }

if (Method == 1)
  {
    htitleA = "hData_CRlowlowpt_1Vtx_SumtrackWeight_";
    htitleB = "hData_CRlooselowlowpt_1Vtx_SumtrackWeight_";
    htitleC = "hData_CRTight_1Vtx_SumtrackWeight_";
    htitleD = "hData_CRLoose_1Vtx_SumtrackWeight_";

    nbin = 40; 
    xmin = 0;
    xmax = 40;
    HeaderA = "T + 30<pt_{i}<80 & pt_{j}>80";
    HeaderB = "L + 30<pt<80";
    HeaderC = "T + 30<pt_{i}<80 & pt_{j}>80";
    HeaderD = "L + 30<pt<80";

    HeaderNVtx = "1 Vtx";
    xtitle = "SumtrackWeight";
  }

        // !! LT !! //

if (Method == 2)
  {
    htitleA = "Quality_LT_STW_1Vtx_A_";
    htitleB = "Quality_LT_STW_1Vtx_B_";
    htitleC = "Quality_LT_STW_1Vtx_C_";
    htitleD = "Quality_LT_STW_1Vtx_D_";

    nbin = 40; 
    xmin = 0;
    xmax = 40;
    HeaderA = "T + L_{T}<80";
    HeaderB = "L + L_{T}<80";
    HeaderC = "T + L_{T}>80";
    HeaderD = "L + L_{T}>80";

    HeaderNVtx = "1 Vtx";
    xtitle = "SumtrackWeight";
  }

  TLegend* leg;

// Couleurs pour les histogrammes
    Float_t r1 = 0.246;
Float_t g1 = 0.563;
Float_t b1 = 0.852;
TColor color1 = TColor(301,r1, g1, b1);
// color1.SetRGB(r1, g1, b1);
Int_t ColorBlue = color1.GetNumber();

Float_t r2 = 1.000;
Float_t g2 = 0.661;
Float_t b2 = 0.055;
TColor color2 = TColor(302,r2, g2, b2);
// color2.SetRGB(r2, g2, b2);
Int_t ColorOrange = color2.GetNumber();

Float_t r3 = 0.739;
Float_t g3 = 0.122;
Float_t b3 = 0.004;
TColor color3 = TColor(303,r3, g3, b3);
Int_t ColorRed = color3.GetNumber();

Float_t r4 = 0.578;
Float_t g4 = 0.641;
Float_t b4 = 0.635;
TColor color4 = TColor(304,r4, g4, b4);
Int_t ColorGrey = color4.GetNumber();

Float_t r5 = 0.513;
Float_t g5 = 0.176;
Float_t b5 = 0.713;
TColor color5 = TColor(305,r5, g5, b5);
Int_t ColorDarkPurple = color5.GetNumber();

Float_t r6 = 0.661;
Float_t g6 = 0.418;
Float_t b6 = 0.348;
TColor color6 = TColor(306,r6, g6, b6);
Int_t ColorBrown = color6.GetNumber();

Float_t r7 = 0.905;
Float_t g7 = 0.387;
Float_t b7 = 0.000;
TColor color7 = TColor(307,r7, g7, b7);
Int_t ColorDarkOrange = color7.GetNumber();

Float_t r8 = 0.723;
Float_t g8 = 0.672;
Float_t b8 = 0.438;
TColor color8 = TColor(308,r8, g8, b8);
Int_t ColorNeutral = color8.GetNumber();

Float_t r9 = 0.441;
Float_t g9 = 0.457;
Float_t b9 = 0.504;
TColor color9 = TColor(309,r9, g9, b9);
Int_t ColorDarkGrey = color9.GetNumber();

Float_t r10 = 0.571;
Float_t g10 = 0.852;
Float_t b10 = 0.867;
TColor color10 = TColor(310,r10, g10, b10);
// color10.SetRGB(r10, g10, b10);
Int_t ColorLightBlue = color10.GetNumber();
// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

// !! -------------------
// !! Pads
float x1 = 0.03;
float x2 = 0.49;
float x3 = 0.51;
float x4 = 0.97;

float y1 = 0.62;
float y2 = 0.92;
float y3 = 0.17;
float y4 = 0.47;

// !! raps

float x5 = 0.03;
float x6 = 0.49;
float x7 = 0.51;
float x8 = 0.97;

float y5 = 0.48;
float y6 = 0.61;
float y7 = 0.03;
float y8 = 0.16;
// !! ------------------
if ((DATA && !MC ) || (!DATA && MC))
{

  y1 = 0.51;
  y2 = 0.92;
  y3 = 0.03;
  y4 = 0.44;

  y5 = 0.0;
  y6 = 0.005;
  y7 = 0.0;
  y8 = 0.0005;


}


  TPad* pad1 = new TPad("pad1","This is pad1",x1,y1,x2,y2,21);
TPad* pad2 = new TPad("pad2","This is pad2",x3,y1,x4,y2,21);
TPad* pad3 = new TPad("pad3","This is pad3",x1,y3,x2,y4,21);
TPad* pad4 = new TPad("pad4","This is pad4",x3,y3,x4,y4,21);

  TPad* pad8 = new TPad("pad8","This is pad8",x5,y5,x6,y6,21);
  TPad* pad9 = new TPad("pad9","This is pad9",x7,y5,x8,y6,21);
  TPad* pad10 = new TPad("pad10","This is pad10",x5,y7,x6,y8,21);
  TPad* pad11 = new TPad("pad11","This is pad11",x7,y7,x8,y8,21);

float bottomMargin = 0.05;
if ((DATA && !MC) || ( !DATA && MC )) bottomMargin = 0.15;
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.07);
   pad1->SetBottomMargin(bottomMargin);
   pad1->SetRightMargin(0.04);
   pad1->SetLeftMargin(0.16);

pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(logy);
   pad2->SetTopMargin(0.07);
   pad2->SetBottomMargin(bottomMargin);
   pad2->SetRightMargin(0.04);
   pad2->SetLeftMargin(0.16);

pad3->SetFillColor(0);
pad3->SetBorderMode(0);
pad3->SetFrameFillColor(10);
pad3->Draw();
pad3->SetLogy(logy);
   pad3->SetTopMargin(0.07);
   pad3->SetBottomMargin(bottomMargin);
   pad3->SetRightMargin(0.04);
   pad3->SetLeftMargin(0.16);

pad4->SetFillColor(0);
pad4->SetBorderMode(0);
pad4->SetFrameFillColor(10);
pad4->Draw();
pad4->SetLogy(logy);
   pad4->SetTopMargin(0.07);
   pad4->SetBottomMargin(bottomMargin);
   pad4->SetRightMargin(0.04);
   pad4->SetLeftMargin(0.16);

pad8->SetFillColor(0);
pad8->SetBorderMode(0);
pad8->SetFrameFillColor(10);
pad8->Draw();
pad8->SetLogy(0);
   pad8->SetTopMargin(0.07);
   pad8->SetBottomMargin(0.35);
   pad8->SetRightMargin(0.04);
   pad8->SetLeftMargin(0.16);

pad9->SetFillColor(0);
pad9->SetBorderMode(0);
pad9->SetFrameFillColor(10);
pad9->Draw();
pad9->SetLogy(0);
   pad9->SetTopMargin(0.07);
   pad9->SetBottomMargin(0.35);
   pad9->SetRightMargin(0.04);
   pad9->SetLeftMargin(0.16);

pad10->SetFillColor(0);
pad10->SetBorderMode(0);
pad10->SetFrameFillColor(10);
pad10->Draw();
pad10->SetLogy(0);
   pad10->SetTopMargin(0.07);
   pad10->SetBottomMargin(0.35);
   pad10->SetRightMargin(0.04);
   pad10->SetLeftMargin(0.16);

pad11->SetFillColor(0);
pad11->SetBorderMode(0);
pad11->SetFrameFillColor(10);
pad11->Draw();
pad11->SetLogy(0);
   pad11->SetTopMargin(0.07);
   pad11->SetBottomMargin(0.35);
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
TH1F* g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleA);

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
 TH1F* g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);//ok
  TH1F* g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);//ok
 TH1F* g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);//ok
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
float scaleMC = 0.948;// % of jobs that did not fail for emu data
 pad1->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

 // !! ------------------------------- !!//
 TH1F* hSumQuadratic_Up = new TH1F("hSumQuadratic_Up","",nbin,xmin,xmax);
 TH1F* hSumQuadratic_Down = new TH1F("hSumQuadratic_Down","",nbin,xmin,xmax);

  for (unsigned int i_up = 0 ; i_up < f_SYST_Up.size(); i_up++)
    {

        f_SYST_Up[i_up]->cd();
        TH1F* h_Temp = (TH1F*)gROOT->FindObject(htitleA+"_"+SYSTNAME_UP[i_up]);
        for (int bin = 1; bin <= h_Temp->GetNbinsX(); ++bin) {
            double val = h_Temp->GetBinContent(bin)-1; // !! rescale around 0 ;
            double sumVal = hSumQuadratic_Up->GetBinContent(bin);
            // std::cout<< "valUp = " << val << " sumValUpb4 = " << sumVal << std::endl;
            hSumQuadratic_Up->SetBinContent(bin, sumVal + val * val);
            // std::cout<< "valUp = " << val << " sumValUpafter = " << hSumQuadratic_Up->GetBinContent(bin) << std::endl;
        }

        f_SYST_Down[i_up]->cd();
        h_Temp = (TH1F*)gROOT->FindObject(htitleA+"_"+SYSTNAME_DOWN[i_up]);
        for (int bin = 1; bin <= h_Temp->GetNbinsX(); ++bin) {
            double val = h_Temp->GetBinContent(bin) -1 ; // !! rescale around 0 ;
            double sumVal = hSumQuadratic_Down->GetBinContent(bin);
            hSumQuadratic_Down->SetBinContent(bin, sumVal + val * val);
        }        
    }

    if (hSumQuadratic_Up) {
        for (int bin = 1; bin <= hSumQuadratic_Up->GetNbinsX(); ++bin) {
            hSumQuadratic_Up->SetBinContent(bin, sqrt(hSumQuadratic_Up->GetBinContent(bin)));
            // std::cout<< "sumValUp = " << hSumQuadratic_Up->GetBinContent(bin) << std::endl;
        }
    }
    if (hSumQuadratic_Down) {
        for (int bin = 1; bin <= hSumQuadratic_Down->GetNbinsX(); ++bin) {
            hSumQuadratic_Down->SetBinContent(bin, sqrt(hSumQuadratic_Down->GetBinContent(bin)));
        }
    }


// !! ------------------------------- !!//

if (MC)
  {
    f1_DY->cd();
    g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
    h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
    h_DY->Add(g1_DY, h_DY, 1,0);

    f2_DY->cd();
    g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
    h_DY->Add(g2_DY, h_DY, 1, 1);

    f1_VV->cd();
    g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);
    h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
    h_VV->Add(g1_VV, h_VV, 1,0);

    f2_VV->cd();
    g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);
    h_VV->Add(g2_VV, h_VV,1, 1);

    f3_VV->cd();
    g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);
    h_VV->Add(g3_VV, h_VV, 1, 1);

    f1_TTV->cd();
    g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);
    h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
    h_TTV->Add(g1_TTV, h_TTV, 1,0);

    f2_TTV->cd();
    g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);
    h_TTV->Add(g2_TTV, h_TTV, 1, 1);

    f3_TTV->cd();
    g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);
    h_TTV->Add(g3_TTV, h_TTV, 1, 1);

    h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
    f1_ST->cd();
    g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);

    h_ST->Add(g1_ST, h_ST, 1,0);

    f2_ST->cd();
    g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);
    h_ST->Add(g2_ST, h_ST, 1, 1);

      f3_ST->cd();
    g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);
    h_ST->Add(g3_ST, h_ST, 1, 1);

      f4_ST->cd();
    g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);
    h_ST->Add(g4_ST, h_ST, 1, 1);

    f1_TT->cd();
    g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
    h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
    h_TT->Add(g1_TT, h_TT, rwTT*1,0);

    f2_TT->cd();
    g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);
    h_TT->Add(g2_TT, h_TT, 1,1);

    //!! --------------------------
    h_DY->Scale(1*scaleMC);
    h_VV->Scale(1*scaleMC);
    h_TTV->Scale(1*scaleMC);
    h_ST->Scale(1*scaleMC);
    h_TT->Scale(1*scaleMC);

    // Daniel 
    h_VV->Add(h_DY, h_VV, 1, 1);
    h_ST->Add(h_ST, h_VV, 1, 1);
    h_TTV->Add(h_TTV, h_ST, 1, 1);
    htotMC->Add(h_TTV, htotMC, 1, 0);
    htotMC->Add(htotMC,h_TT, 1, 1);

    htotMC->Draw("HE"); 
      htotMC->SetFillStyle(1001);
    htotMC->SetFillColorAlpha(ColorRed, 1);
    htotMC->SetLineColor(ColorRed);
    htotMC->SetLineColorAlpha(ColorRed, 1);
    htotMC->SetLineStyle(1);
    htotMC->SetLineWidth(1);
    htotMC->SetTickLength(0.03, "YZ");
    htotMC->SetTickLength(0.03,"X");
    htotMC->SetLabelOffset(0.015,"X");
    htotMC->SetLabelOffset(0.007,"Y");
    htotMC->SetLabelSize(0.045, "XYZ");
    htotMC->SetLabelFont(42, "XYZ"); 
    htotMC->SetTitleSize(0.055, "XYZ"); 
    htotMC->SetTitleFont(42, "XYZ");
    htotMC->SetTitleOffset(1.2,"X"); 
    htotMC->SetTitleOffset(1.3,"Y");
    htotMC->GetXaxis()->SetTitle(xtitle);
    htotMC->GetXaxis()->SetTitleColor(1);
    htotMC->GetYaxis()->SetTitle(ytitle);
    htotMC->GetYaxis()->SetTitleColor(1);
    htotMC->SetNdivisions(509,"XYZ");
    //  htotMC->SetMinimum(hmin); 
    //  htotMC->SetMaximum(hmax); 
    //  htotMC->SetMarkerStyle(20);
    //  htotMC->SetMarkerSize(1);
    if (logy)
      {
        htotMC->SetMinimum(1); 
        htotMC->SetMaximum(htotMC->GetMaximum()*1000); 
      }
    else 
      {
        htotMC->SetMinimum(1); 
        htotMC->SetMaximum(htotMC->GetMaximum()*2); 
      }

    //  h_TTV->Draw("HEsame"); 
    //  h_TTV->SetFillColorAlpha(ColorNeutral, 1);
    //  h_TTV->SetLineColor(ColorNeutral);
    //  h_TTV->SetLineStyle(1);
    //  h_TTV->SetLineWidth(3);

    h_ST->Draw("HEsame"); 
    h_ST->SetFillColorAlpha(ColorDarkPurple, 1);
    h_ST->SetLineColor(ColorDarkPurple);
    h_ST->SetLineColorAlpha(ColorDarkPurple, 1);
    h_ST->SetLineStyle(1);
    h_ST->SetLineWidth(3);

    //  h_TT->Draw("HEsame"); 
    //  h_TT->SetFillColorAlpha(ColorRed, 1);
    //  h_TT->SetLineColor(ColorRed);
    //  h_TT->SetLineStyle(1);
    //  h_TT->SetLineWidth(3);
    //  h_TT->SetTickLength(0.03, "YZ");
    //  h_TT->SetTickLength(0.03,"X");

    h_VV->Draw("HEsame"); 
    h_VV->SetFillColorAlpha(ColorOrange, 1);
    h_VV->SetLineColor(ColorOrange);
    h_VV->SetLineColorAlpha(ColorOrange, 1);
    h_VV->SetLineStyle(1);
    h_VV->SetLineWidth(3);

      h_DY->Draw("HEsame"); 
    h_DY->SetFillColorAlpha(ColorBlue,1);
    h_DY->SetLineColor(ColorBlue);
    h_DY->SetLineColorAlpha(ColorBlue,1);
    h_DY->SetLineStyle(1);
    h_DY->SetLineWidth(3);
    //  h_DY->SetTickLength(0.03, "YZ");
    //  h_DY->SetTickLength(0.03,"X");
  }
if (DATA)
  {
    f1_Data_emu->cd();
    g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleA);
    htotData->Add(g1_Data_emu, htotData, 1, 0);

    htotData->Draw("PE1same");
    htotData->SetMarkerStyle(20);
    htotData->SetMarkerSize(1);
    htotData->SetMarkerColor(kBlack);
    htotData->SetLineColor(kBlack);
    htotData->SetLineWidth(1);

    if (DATA && !MC)
      {
        htotData->SetTickLength(0.03, "YZ");
        htotData->SetTickLength(0.03,"X");
        htotData->SetLabelOffset(0.015,"X");
        htotData->SetLabelOffset(0.007,"Y");
        htotData->SetLabelSize(0.045, "XYZ");
        htotData->SetLabelFont(42, "XYZ"); 
        htotData->SetTitleSize(0.055, "XYZ"); 
        htotData->SetTitleFont(42, "XYZ");
        htotData->SetTitleOffset(1.2,"X"); 
        htotData->SetTitleOffset(1.3,"Y");
            htotData->GetXaxis()->SetTitle(xtitle);
        htotData->GetXaxis()->SetTitleColor(1);
        htotData->GetYaxis()->SetTitle(ytitle);
        htotData->GetYaxis()->SetTitleColor(1);
        if (logy)
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*100); 
          }
        else 
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*2); 
          }
      }

  }
if (Signal)
  {
    f1_LLP->cd();
    
    g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleA);
    g1_LLP->Sumw2();
    h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
    h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

    f2_LLP->cd();
    g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleA);
    g2_LLP->Sumw2();
    h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
    h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

    f3_LLP->cd();
    
    g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleA);
    g3_LLP->Sumw2();
    h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
    h3_LLP->Add(g3_LLP, h3_LLP, 1,0);


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
  }
 
 if (DATA )
  {
    hsolve->Add(hsolve, htotData, 0., 1.);
  }
 else if (MC && !DATA)
  {
    hsolve->Add(hsolve, htotMC, 0., 1.);
  }

  leg = new TLegend(0.17,0.94,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  float LEGY1 = 0.50;
  float LEGY2 = 0.89;
  if ( (DATA && !MC) || (!DATA && MC))
    {
      LEGY1 = 0.70;
      LEGY2 = 0.89;
    }
  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
  if (DATA)
    {
      leg->AddEntry(htotData, " e#mu data","PE1");
    }  
  if (MC)
    {
      leg->AddEntry(htotMC, " t#bar{t} + X","F");
      leg->AddEntry(h_ST, " Single t","F");
      leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
      leg->AddEntry(h_DY, " DY","F");
    }
  if (Signal)
    {
      leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
      leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
      leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
    }


  leg->Draw();


  leg = new TLegend(0.2,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();

  leg = new TLegend(0.2,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderA);
  leg->Draw();

  //Data/MC ---------------------
  if (DATA  && MC)
    {
      pad8->cd();
      hDataMC->Divide(htotData,htotMC,1,1);
      hDataMC->Draw("E"); 
      hDataMC->SetLineColor(1);
      hDataMC->SetLineStyle(1);
      hDataMC->SetLineWidth(1);
      hDataMC->SetMarkerColor(kBlack);
      hDataMC->SetMarkerStyle(20);
      hDataMC->SetMarkerSize(0.6);
      hDataMC->SetTickLength(0.10, "X"); hDataMC->SetTickLength(0.05, "YZ");
      hDataMC->SetLabelOffset(0.02,"X");
      hDataMC->SetLabelOffset(0.02,"Y");
      hDataMC->SetLabelSize(0.12, "XY");
      hDataMC->SetLabelFont(42, "XYZ"); 
      hDataMC->SetTitleFont(42, "XYZ");
      hDataMC->SetTitleSize(0.14, "XYZ"); 
      hDataMC->SetTitleOffset(0.9,"X");
      hDataMC->SetTitleOffset(0.5,"Y");
      hDataMC->GetXaxis()->SetTitle(xtitle);
      hDataMC->GetXaxis()->SetTitleColor(1);
      hDataMC->GetXaxis()->SetNdivisions(509);
      hDataMC->GetYaxis()->SetTitle("Data / Sim.");
      hDataMC->GetYaxis()->SetTitleColor(1);
      hDataMC->GetYaxis()->SetNdivisions(509);
      hDataMC->SetNdivisions(509,"XYZ");
      hDataMC->SetMinimum(0.5); 
      hDataMC->SetMaximum(1.5); 

        const TH1F* hRatio_const(hDataMC);

        TGraphAsymmErrors *hRatioUp = new TGraphAsymmErrors(hRatio_const);
        TGraphAsymmErrors *hRatioDown = new TGraphAsymmErrors(hRatio_const);

        for (unsigned int b = 1 ; b < hDataMC->GetNbinsX(); b++)
        {
          float htotMCbin = htotMC->GetBinContent(b);
          float htotMCbinErrorStat = 0;
          if (htotMC->GetBinError(b) != 0)
            {
              htotMCbinErrorStat = htotMC->GetBinError(b);
            }

          float bin =  hDataMC->GetBinContent(b);
          float binStatError =  hDataMC->GetBinError(b);

          float htotDATAbinErrorStat  = 0;
          if (htotData->GetBinError(b) != 0)
            {
              htotDATAbinErrorStat = htotData->GetBinError(b);
            }
            float binSysErrorUp = hSumQuadratic_Up->GetBinContent(b);
          float binSysErrorDown = hSumQuadratic_Down->GetBinContent(b);

          float errorTotalUp = sqrt(binSysErrorUp*binSysErrorUp + binStatError*binStatError);// unccorrelated SYS uncertainties htotMCbinErrorStat*htotMCbinErrorStat
          float errorTotalDown = sqrt(binSysErrorDown*binSysErrorDown + binStatError*binStatError);// unccorrelated SYS uncertainties htotMCbinErrorStat*htotMCbinErrorStat

          float xbin = hDataMC->GetBinCenter(b);
          hRatioUp->SetPoint(b,xbin,1);
          hRatioDown->SetPoint(b,xbin,1);

          hRatioUp->SetPointError(b,0,0,0,errorTotalUp);
          hRatioDown->SetPointError(b,0,0,errorTotalDown,0);
          // hDataMC->SetBinError(b, htotDATAbinErrorStat);


        }

          hRatioUp->Draw("E3same");//HISTSAME
          hRatioUp->SetFillColor(kGray+2);
          hRatioUp->SetFillStyle(3001);
          hRatioUp->SetFillColorAlpha(kGray+2, 0.9);


          hRatioUp->SetLineStyle(1);
          hRatioUp->SetLineWidth(1);
          // hRatioUp->SetLineColor(ColorBlue);
          // hRatioUp->SetLineColorAlpha(ColorBlue,0.8);

          hRatioDown->Draw("E3same");//HISTSAME
          // hRatioDown->SetLineColor(ColorRed);
          // hRatioDown->SetLineColorAlpha(ColorRed,0.8);
          hRatioDown->SetLineStyle(1);
          hRatioDown->SetLineWidth(1);
          hRatioDown->SetFillColorAlpha(kGray+2, 0.9);
          hRatioDown->SetFillColor(kGray+2);
          hRatioDown->SetFillStyle(3001); 

    }


// *****************************************************************************

 pad2->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 
 // !! ------------------------------- !!//
 hSumQuadratic_Up = new TH1F("hSumQuadratic_Up","",nbin,xmin,xmax);
 hSumQuadratic_Down = new TH1F("hSumQuadratic_Down","",nbin,xmin,xmax);

  for (unsigned int i_up = 0 ; i_up < f_SYST_Up.size(); i_up++)
    {

        f_SYST_Up[i_up]->cd();
        TH1F* h_Temp = (TH1F*)gROOT->FindObject(htitleB+"_"+SYSTNAME_UP[i_up]);
        for (int bin = 1; bin <= h_Temp->GetNbinsX(); ++bin) {
            double val = h_Temp->GetBinContent(bin)-1; // !! rescale around 0 ;
            double sumVal = hSumQuadratic_Up->GetBinContent(bin);
            // std::cout<< "valUp = " << val << " sumValUpb4 = " << sumVal << std::endl;
            hSumQuadratic_Up->SetBinContent(bin, sumVal + val * val);
            // std::cout<< "valUp = " << val << " sumValUpafter = " << hSumQuadratic_Up->GetBinContent(bin) << std::endl;
        }

        f_SYST_Down[i_up]->cd();
        h_Temp = (TH1F*)gROOT->FindObject(htitleB+"_"+SYSTNAME_DOWN[i_up]);
        for (int bin = 1; bin <= h_Temp->GetNbinsX(); ++bin) {
            double val = h_Temp->GetBinContent(bin) -1 ; // !! rescale around 0 ;
            double sumVal = hSumQuadratic_Down->GetBinContent(bin);
            hSumQuadratic_Down->SetBinContent(bin, sumVal + val * val);
        }        
    }

    if (hSumQuadratic_Up) {
        for (int bin = 1; bin <= hSumQuadratic_Up->GetNbinsX(); ++bin) {
            hSumQuadratic_Up->SetBinContent(bin, sqrt(hSumQuadratic_Up->GetBinContent(bin)));
            // std::cout<< "sumValUp = " << hSumQuadratic_Up->GetBinContent(bin) << std::endl;
        }
    }
    if (hSumQuadratic_Down) {
        for (int bin = 1; bin <= hSumQuadratic_Down->GetNbinsX(); ++bin) {
            hSumQuadratic_Down->SetBinContent(bin, sqrt(hSumQuadratic_Down->GetBinContent(bin)));
        }
    }
     // !! ------------------------------- !!//

if (MC)
  {
      f1_DY->cd();
      g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);

      h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
      h_DY->Add(g1_DY, h_DY, 1,0);

      f2_DY->cd();
      g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleC);
      h_DY->Add(g2_DY, h_DY, 1, 1);

      f1_VV->cd();
      g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleC);

      h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
      h_VV->Add(g1_VV, h_VV, 1,0);

      f2_VV->cd();
      g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleC);

      h_VV->Add(g2_VV, h_VV, 1, 1);
      f3_VV->cd();
      
      g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleC);
      h_VV->Add(g3_VV, h_VV, 1, 1);

      f1_TTV->cd();
      g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleC);

      h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
      h_TTV->Add(g1_TTV, h_TTV, 1,0);

      f2_TTV->cd();
      g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleC);
      h_TTV->Add(g2_TTV, h_TTV, 1, 1);

      f3_TTV->cd();
      g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleC);
      h_TTV->Add(g3_TTV, h_TTV, 1, 1);

      h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
      f1_ST->cd();
      g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleC);
      h_ST->Add(g1_ST, h_ST, 1,0);

      f2_ST->cd();
      g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleC);
      h_ST->Add(g2_ST, h_ST, 1, 1);

        f3_ST->cd();
      g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleC);
      h_ST->Add(g3_ST, h_ST, 1, 1);

        f4_ST->cd();
      g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleC);
      h_ST->Add(g4_ST, h_ST, 1, 1);


      f1_TT->cd();
      g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);

      h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
      h_TT->Add(g1_TT, h_TT, 1,0);

      f2_TT->cd();
      g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleC);
      h_TT->Add(g2_TT, h_TT, 1,1);

      //!! --------------------------
      h_DY->Scale(1*scaleMC);
      h_VV->Scale(1*scaleMC);
      h_TTV->Scale(1*scaleMC);
      h_ST->Scale(1*scaleMC);
      h_TT->Scale(1*scaleMC);
      // Daniel 
      h_VV->Add(h_DY, h_VV, 1, 1);
      h_ST->Add(h_ST, h_VV, 1, 1);
      h_TTV->Add(h_TTV, h_ST, 1, 1);
      htotMC->Add(h_TTV, htotMC, 1, 0);
      htotMC->Add(htotMC,h_TT, 1, 1);
      // -----

      htotMC->Draw("HE"); 
      htotMC->SetFillStyle(1001);
      htotMC->SetFillColorAlpha(ColorRed, 1);
      htotMC->SetLineColor(ColorRed);
      htotMC->SetLineColorAlpha(ColorRed, 1);
      htotMC->SetLineStyle(1);
      htotMC->SetLineWidth(1);
      htotMC->SetTickLength(0.03, "YZ");
      htotMC->SetTickLength(0.03,"X");
      htotMC->SetLabelOffset(0.015,"X");
      htotMC->SetLabelOffset(0.007,"Y");
      htotMC->SetLabelSize(0.045, "XYZ");
      htotMC->SetLabelFont(42, "XYZ"); 
      htotMC->SetTitleSize(0.055, "XYZ"); 
      htotMC->SetTitleFont(42, "XYZ");
      htotMC->SetTitleOffset(1.2,"X"); 
      htotMC->SetTitleOffset(1.3,"Y");
      htotMC->GetXaxis()->SetTitle(xtitle);
      htotMC->GetXaxis()->SetTitleColor(1);
      htotMC->GetYaxis()->SetTitle(ytitle);
      htotMC->GetYaxis()->SetTitleColor(1);
      htotMC->SetNdivisions(509,"XYZ");

      if (logy)
        {
          htotMC->SetMinimum(1); 
          htotMC->SetMaximum(htotMC->GetMaximum()*1000); 
        }
      else 
        {
          htotMC->SetMinimum(1); 
          htotMC->SetMaximum(htotMC->GetMaximum()*2); 
        }



      h_ST->Draw("HEsame"); 
      h_ST->SetFillColorAlpha(ColorDarkPurple, 1);
      h_ST->SetLineColor(ColorDarkPurple);
      h_ST->SetLineColorAlpha(ColorDarkPurple, 1);
      h_ST->SetLineStyle(1);
      h_ST->SetLineWidth(3);


      h_VV->Draw("HEsame"); 
      h_VV->SetFillColorAlpha(ColorOrange, 1);
      h_VV->SetLineColor(ColorOrange);
      h_VV->SetLineColorAlpha(ColorOrange, 1);
      h_VV->SetLineStyle(1);
      h_VV->SetLineWidth(3);

        h_DY->Draw("HEsame"); 
      h_DY->SetFillColorAlpha(ColorBlue,1);
      h_DY->SetLineColor(ColorBlue);
      h_DY->SetLineColorAlpha(ColorBlue,1);
      h_DY->SetLineStyle(1);
      h_DY->SetLineWidth(3);
      h_DY->SetTickLength(0.03, "YZ");
      h_DY->SetTickLength(0.03,"X");
  }
 if (DATA)
  {
    f1_Data_emu->cd();
    g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleC);
    htotData->Add(g1_Data_emu, htotData, 1, 0);
    htotData->Draw("PE1same");
    htotData->SetMarkerStyle(20);
    htotData->SetMarkerSize(1);
    htotData->SetMarkerColor(kBlack);
    htotData->SetLineColor(kBlack);
      htotData->SetLineWidth(1);
          if (DATA && !MC)
      {
        htotData->SetTickLength(0.03, "YZ");
        htotData->SetTickLength(0.03,"X");
        htotData->SetLabelOffset(0.015,"X");
        htotData->SetLabelOffset(0.007,"Y");
        htotData->SetLabelSize(0.045, "XYZ");
        htotData->SetLabelFont(42, "XYZ"); 
        htotData->SetTitleSize(0.055, "XYZ"); 
        htotData->SetTitleFont(42, "XYZ");
        htotData->SetTitleOffset(1.2,"X"); 
        htotData->SetTitleOffset(1.3,"Y");
            htotData->GetXaxis()->SetTitle(xtitle);
        htotData->GetXaxis()->SetTitleColor(1);
        htotData->GetYaxis()->SetTitle(ytitle);
        htotData->GetYaxis()->SetTitleColor(1);
                if (logy)
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*100); 
          }
        else 
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*2); 
          }
      }
  }


if(Signal)
  {
    f1_LLP->cd();
    
    
    g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleC);
    g1_LLP->Sumw2();
    h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
    h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

    f2_LLP->cd();
    
    
    g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleC);
    g2_LLP->Sumw2();
    h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
    h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

    f3_LLP->cd();
    
    g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleC);
    g3_LLP->Sumw2();
    h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
    h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

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
  }


  leg = new TLegend(0.17,0.94,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
  if(DATA)
    {
      leg->AddEntry(htotData, " e#mu data","PE1");
    }
  if(MC)
    {
      leg->AddEntry(htotMC, " t#bar{t} + X","F");
      leg->AddEntry(h_ST, " Single t","F");
      leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
      leg->AddEntry(h_DY, " DY","F");
    }

  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
  leg->AddEntry(hsolve," Prediction","FE4");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();


  leg = new TLegend(0.2,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();

  leg = new TLegend(0.2,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderC);
  leg->Draw();

if(DATA && MC)
  {
    pad9->cd();
    
    // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
    hDataMC->Divide(htotData,htotMC,1,1);
    hDataMC->Draw("E"); 
    hDataMC->SetLineColor(1);
    hDataMC->SetLineStyle(1);
    hDataMC->SetLineWidth(1);
    hDataMC->SetMarkerColor(kBlack);
    hDataMC->SetMarkerStyle(20);
    hDataMC->SetMarkerSize(0.6);
    hDataMC->SetTickLength(0.10, "X"); hDataMC->SetTickLength(0.05, "YZ");
    hDataMC->SetLabelOffset(0.02,"X");
    hDataMC->SetLabelOffset(0.02,"Y");
    hDataMC->SetLabelSize(0.12, "XY");
    hDataMC->SetLabelFont(42, "XYZ"); 
    hDataMC->SetTitleFont(42, "XYZ");
    hDataMC->SetTitleSize(0.14, "XYZ"); 
    hDataMC->SetTitleOffset(0.9,"X");
    hDataMC->SetTitleOffset(0.5,"Y");
    hDataMC->GetXaxis()->SetTitle(xtitle);
    hDataMC->GetXaxis()->SetTitleColor(1);
    hDataMC->GetXaxis()->SetNdivisions(509);
    hDataMC->GetYaxis()->SetTitle("Data / Sim.");
    hDataMC->GetYaxis()->SetTitleColor(1);
    hDataMC->GetYaxis()->SetNdivisions(509);
    hDataMC->SetNdivisions(509,"XYZ");
    hDataMC->SetMinimum(0.5); 
    hDataMC->SetMaximum(1.5);



    const TH1F* hRatio_constB(hDataMC);

    hRatioUp = new TGraphAsymmErrors(hRatio_constB);
    hRatioDown = new TGraphAsymmErrors(hRatio_constB);

    for (unsigned int b = 1 ; b < hDataMC->GetNbinsX(); b++)
    {
      float htotMCbin = htotMC->GetBinContent(b);
      float htotMCbinErrorStat = htotMC->GetBinError(b);

      float bin =  hDataMC->GetBinContent(b);
      float binStatError =  hDataMC->GetBinError(b);

      float binSysErrorUp = hSumQuadratic_Up->GetBinContent(b);
      float binSysErrorDown = hSumQuadratic_Down->GetBinContent(b);
      // std::cout<<"bin:     "<<bin<<" binStatError: "<<binStatError<<" hSumQuadratic_Up->GetBinContent(b): "<<hSumQuadratic_Up->GetBinContent(b)<<" hSumQuadratic_Down->GetBinContent(b): "<<hSumQuadratic_Down->GetBinContent(b)<<std::endl;
      // std::cout<<"bin "<<bin<<" htotMC->GetBinError(b) : "<<htotMC->GetBinError(b)<<" and hDataMC->GetBinError(b): "<<hDataMC->GetBinError(b)<<std::endl;
      float errorTotalUp = sqrt(binSysErrorUp*binSysErrorUp + binStatError*binStatError);// unccorrelated SYS uncertainties
      float errorTotalDown = sqrt(binSysErrorDown*binSysErrorDown + binStatError*binStatError);// unccorrelated SYS uncertainties
      // std::cout<<"bin: "<<bin<<" binStatError: "<<binStatError<<" binSysErrorUp: "<<binSysErrorUp<<" binSysErrorDown: "<<binSysErrorDown<<" errorTotalUp: "<<errorTotalUp<<" errorTotalDown: "<<errorTotalDown<<std::endl;
      float xbin = hDataMC->GetBinCenter(b);
      hRatioUp->SetPoint(b,xbin,1);
      hRatioDown->SetPoint(b,xbin,1);

      hRatioUp->SetPointError(b,0,0,0,errorTotalUp);
      hRatioDown->SetPointError(b,0,0,errorTotalDown,0);
      // hRatioUp->SetBinContent(b, bin);
      // hRatioDown->SetBinContent(b, bin);

      // hRatioUp->SetBinError(b, errorTotalUp);
      // hRatioDown->SetBinError(b, errorTotalDown);

    }

    hRatioUp->Draw("E3same");//HISTSAME
    hRatioUp->SetFillColor(kGray+2);
    hRatioUp->SetFillStyle(3001);
    hRatioUp->SetFillColorAlpha(kGray+2, 0.9);
      hRatioUp->SetLineStyle(1);
      hRatioUp->SetLineWidth(1);

    hRatioDown->Draw("E3same");//HISTSAME
    hRatioDown->SetLineStyle(1);
    hRatioDown->SetLineWidth(1);
    hRatioDown->SetFillColorAlpha(kGray+2, 0.9);
    hRatioDown->SetFillColor(kGray+2);
    hRatioDown->SetFillStyle(3001); 
  }


 /// !! --------------------------------------------------------------------------//
 /// !! --------------------------------------------------------------------------//
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

 /// !! --------------------------------------------------------------------------//
 /// !! --------------------------------------------------------------------------//
c2->cd();
  pad5->cd();
  if (DATA)
  {
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
  }  

if ( MC) { 
   htotMC->Draw("HESAME"); 
  // htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(ColorRed, 1);
 htotMC->SetLineColor(ColorRed);
 htotMC->SetLineColorAlpha(ColorRed, 1);
 htotMC->SetLineStyle(1);
 htotMC->SetLineWidth(1);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColorAlpha(ColorDarkPurple, 1);
 h_ST->SetLineColor(ColorDarkPurple);
 h_ST->SetLineColorAlpha(ColorDarkPurple, 1);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 //  h_TT->Draw("HEsame"); 
 //  h_TT->SetFillColorAlpha(ColorRed, 1);
 //  h_TT->SetLineColor(ColorRed);
 //  h_TT->SetLineStyle(1);
 //  h_TT->SetLineWidth(3);
 //  h_TT->SetTickLength(0.03, "YZ");
 //  h_TT->SetTickLength(0.03,"X");

 h_VV->Draw("HEsame"); 
 h_VV->SetFillColorAlpha(ColorOrange, 1);
 h_VV->SetLineColor(ColorOrange);
 h_VV->SetLineColorAlpha(ColorOrange, 1);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

  h_DY->Draw("HEsame"); 
 h_DY->SetFillColorAlpha(ColorBlue,1);
 h_DY->SetLineColor(ColorBlue);
 h_DY->SetLineColorAlpha(ColorBlue,1);
 h_DY->SetLineStyle(1);
 h_DY->SetLineWidth(3);
 h_DY->SetTickLength(0.03, "YZ");
 h_DY->SetTickLength(0.03,"X");

}

if (Signal)
  {
    f1_LLP->cd();
    g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleC);
    g1_LLP->Sumw2();
    h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
    h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

    f2_LLP->cd();
    g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleC);
    g2_LLP->Sumw2();
    h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
    h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

    f3_LLP->cd();
    g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleC);
    g3_LLP->Sumw2();
    h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
    h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

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
  }



  leg = new TLegend(0.17,0.94,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

 
  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " e#mu data","PE1");
leg->AddEntry(hsolve," Prediction","FE4");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();

  leg = new TLegend(0.2,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderC);
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
TH1F* Inte  = new TH1F("ratio","",nbin,xmin,xmax);
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
 
 /// !! --------------------------------------------------------------------------//
 /// !! --------------------------------------------------------------------------//

  c1->cd();
 pad3->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
 hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
 hmax = hmaxBD; 
 // !! ------------------------------- !!//
 hSumQuadratic_Up = new TH1F("hSumQuadratic_Up","",nbin,xmin,xmax);
 hSumQuadratic_Down = new TH1F("hSumQuadratic_Down","",nbin,xmin,xmax);

  for (unsigned int i_up = 0 ; i_up < f_SYST_Up.size(); i_up++)
    {

        f_SYST_Up[i_up]->cd();
        TH1F* h_Temp = (TH1F*)gROOT->FindObject(htitleC+"_"+SYSTNAME_UP[i_up]);
        for (int bin = 1; bin <= h_Temp->GetNbinsX(); ++bin) {
            double val = h_Temp->GetBinContent(bin)-1; // !! rescale around 0 ;
            double sumVal = hSumQuadratic_Up->GetBinContent(bin);
            // std::cout<< "valUp = " << val << " sumValUpb4 = " << sumVal << std::endl;
            hSumQuadratic_Up->SetBinContent(bin, sumVal + val * val);
            // std::cout<< "valUp = " << val << " sumValUpafter = " << hSumQuadratic_Up->GetBinContent(bin) << std::endl;
        }

        f_SYST_Down[i_up]->cd();
        h_Temp = (TH1F*)gROOT->FindObject(htitleC+"_"+SYSTNAME_DOWN[i_up]);
        for (int bin = 1; bin <= h_Temp->GetNbinsX(); ++bin) {
            double val = h_Temp->GetBinContent(bin) -1 ; // !! rescale around 0 ;
            double sumVal = hSumQuadratic_Down->GetBinContent(bin);
            hSumQuadratic_Down->SetBinContent(bin, sumVal + val * val);
        }        
    }

    if (hSumQuadratic_Up) {
        for (int bin = 1; bin <= hSumQuadratic_Up->GetNbinsX(); ++bin) {
            hSumQuadratic_Up->SetBinContent(bin, sqrt(hSumQuadratic_Up->GetBinContent(bin)));
            // std::cout<< "sumValUp = " << hSumQuadratic_Up->GetBinContent(bin) << std::endl;
        }
    }
    if (hSumQuadratic_Down) {
        for (int bin = 1; bin <= hSumQuadratic_Down->GetNbinsX(); ++bin) {
            hSumQuadratic_Down->SetBinContent(bin, sqrt(hSumQuadratic_Down->GetBinContent(bin)));
        }
    }
  // !! ---------------------------------------

 if(MC)
  {
      f1_DY->cd();

      g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);
      h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
      h_DY->Add(g1_DY, h_DY, 1,0);

      f2_DY->cd();
      g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleB);
      h_DY->Add(g2_DY, h_DY, 1, 1);

      f1_VV->cd();
      g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleB);
      h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
      h_VV->Add(g1_VV, h_VV, 1,0);

      f2_VV->cd();
      g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleB);
      h_VV->Add(g2_VV, h_VV, 1, 1);

      f3_VV->cd();
      g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleB);
      h_VV->Add(g3_VV, h_VV, 1, 1);

      f1_TTV->cd();
      g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleB);
      h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
      h_TTV->Add(g1_TTV, h_TTV, 1,0);

      f2_TTV->cd();
      g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleB);
      h_TTV->Add(g2_TTV, h_TTV, 1, 1);

      f3_TTV->cd(); 
      g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleB);
      h_TTV->Add(g3_TTV, h_TTV, 1, 1);

      h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
      f1_ST->cd();
      g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleB);
      h_ST->Add(g1_ST, h_ST, 1,0);

      f2_ST->cd();
      g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleB);
      h_ST->Add(g2_ST, h_ST, 1, 1);

        f3_ST->cd();
      g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleB);
      h_ST->Add(g3_ST, h_ST, 1, 1);

        f4_ST->cd();
      g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleB);
      h_ST->Add(g4_ST, h_ST, 1, 1);

      f1_TT->cd();
      g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);

      h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
      h_TT->Add(g1_TT, h_TT, rwTT*1,0);

      f2_TT->cd();
      g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleB);
      h_TT->Add(g2_TT, h_TT, 1,1);

      //!! --------------------------
      h_DY->Scale(1*scaleMC);
      h_VV->Scale(1*scaleMC);
      h_TTV->Scale(1*scaleMC);
      h_ST->Scale(1*scaleMC);
      h_TT->Scale(1*scaleMC);
      // Daniel 
      h_VV->Add(h_DY, h_VV, 1, 1);
      h_ST->Add(h_ST, h_VV, 1, 1);
      h_TTV->Add(h_TTV, h_ST, 1, 1);
      htotMC->Add(h_TTV, htotMC, 1, 0);
      htotMC->Add(htotMC,h_TT, 1, 1);
      // -----

      htotMC->Draw("HE"); 
      htotMC->SetFillStyle(1001);
      htotMC->SetFillColorAlpha(ColorRed, 1);
      htotMC->SetLineColor(ColorRed);
      htotMC->SetLineColorAlpha(ColorRed, 1);
      htotMC->SetLineStyle(1);
      htotMC->SetLineWidth(1);
      htotMC->SetTickLength(0.03, "YZ");
      htotMC->SetTickLength(0.03,"X");
      htotMC->SetLabelOffset(0.015,"X");
      htotMC->SetLabelOffset(0.007,"Y");
      htotMC->SetLabelSize(0.045, "XYZ");
      htotMC->SetLabelFont(42, "XYZ"); 
      htotMC->SetTitleSize(0.055, "XYZ"); 
      htotMC->SetTitleFont(42, "XYZ");
      htotMC->SetTitleOffset(1.2,"X"); 
      htotMC->SetTitleOffset(1.3,"Y");
      htotMC->GetXaxis()->SetTitle(xtitle);
      htotMC->GetXaxis()->SetTitleColor(1);
      htotMC->GetYaxis()->SetTitle(ytitle);
      htotMC->GetYaxis()->SetTitleColor(1);
      htotMC->SetNdivisions(509,"XYZ");
      //  htotMC->SetMinimum(1); 
      //  htotMC->SetMaximum(htotMC->GetMaximum()*2); 
      if (logy)
        {
          htotMC->SetMinimum(1); 
      htotMC->SetMaximum(htotMC->GetMaximum()*1000); 
        }
      else 
        {
          htotMC->SetMinimum(1); 
          htotMC->SetMaximum(htotMC->GetMaximum()*2); 
        }
      //  htotMC->SetMarkerStyle(20);
      //  htotMC->SetMarkerSize(1);

      //  h_TTV->Draw("HEsame"); 
      //  h_TTV->SetFillColorAlpha(ColorNeutral, 1);
      //  h_TTV->SetLineColor(ColorNeutral);
      //  h_TTV->SetLineStyle(1);
      //  h_TTV->SetLineWidth(3);

      h_ST->Draw("HEsame"); 
      h_ST->SetFillColorAlpha(ColorDarkPurple, 1);
      h_ST->SetLineColor(ColorDarkPurple);
        h_ST->SetLineColorAlpha(ColorDarkPurple, 1);
      h_ST->SetLineStyle(1);
      h_ST->SetLineWidth(3);

      //  h_TT->Draw("HEsame"); 
      //  h_TT->SetFillColorAlpha(ColorRed, 1);
      //  h_TT->SetLineColor(ColorRed);
      //  h_TT->SetLineStyle(1);
      //  h_TT->SetLineWidth(3);
      //  h_TT->SetTickLength(0.03, "YZ");
      //  h_TT->SetTickLength(0.03,"X");

      h_VV->Draw("HEsame"); 
      h_VV->SetFillColorAlpha(ColorOrange, 1);
      h_VV->SetLineColor(ColorOrange);
      h_VV->SetLineColorAlpha(ColorOrange, 1);
      h_VV->SetLineStyle(1);
      h_VV->SetLineWidth(3);

        h_DY->Draw("HEsame"); 
      h_DY->SetFillColorAlpha(ColorBlue,1);
      h_DY->SetLineColor(ColorBlue);
      h_DY->SetLineColorAlpha(ColorBlue,1);
      h_DY->SetLineStyle(1);
      h_DY->SetLineWidth(3);
      h_DY->SetTickLength(0.03, "YZ");
      h_DY->SetTickLength(0.03,"X");
  }
 


if (DATA)
  {
    f1_Data_emu->cd();
    g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleB);
    htotData->Add(g1_Data_emu, htotData, 1, 0);
    htotData->Draw("PE1same");
    htotData->SetMarkerStyle(20);
    htotData->SetMarkerSize(1);
    htotData->SetMarkerColor(kBlack);
    htotData->SetLineColor(kBlack);
    htotData->SetLineWidth(1);
        if (DATA && !MC)
      {
        htotData->SetTickLength(0.03, "YZ");
        htotData->SetTickLength(0.03,"X");
        htotData->SetLabelOffset(0.015,"X");
        htotData->SetLabelOffset(0.007,"Y");
        htotData->SetLabelSize(0.045, "XYZ");
        htotData->SetLabelFont(42, "XYZ"); 
        htotData->SetTitleSize(0.055, "XYZ"); 
        htotData->SetTitleFont(42, "XYZ");
        htotData->SetTitleOffset(1.2,"X"); 
        htotData->SetTitleOffset(1.3,"Y");
            htotData->GetXaxis()->SetTitle(xtitle);
        htotData->GetXaxis()->SetTitleColor(1);
        htotData->GetYaxis()->SetTitle(ytitle);
        htotData->GetYaxis()->SetTitleColor(1);
                if (logy)
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*100); 
          }
        else 
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*2); 
          }
      }
  }



if(Signal)
  {
     f1_LLP->cd();
 
      g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleB);
      g1_LLP->Sumw2();
      h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
      h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

      f2_LLP->cd();
      
      g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleB);
      g2_LLP->Sumw2();
      h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
      h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

      f3_LLP->cd();
      
      g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleB);
      g3_LLP->Sumw2();
      h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
      h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

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
  }


//  hsolve->Divide(hsolve, htotData, 1., 1.);
  if (DATA )
  {
    hsolve->Divide(hsolve, htotData, 1., 1.);
  }
 else if (MC && !DATA)
  {
    hsolve->Divide(hsolve, htotMC, 1., 1.);
  }

  leg = new TLegend(0.17,0.94,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
  if (DATA)
    {
      leg->AddEntry(htotData, " e#mu data","PE1");
    }
  if (MC)
    {
      leg->AddEntry(htotMC, " t#bar{t} + X","F");
      leg->AddEntry(h_ST, " Single t","F");
      leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
      leg->AddEntry(h_DY, " DY","F");
    }
if (Signal)
  {
    leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
    leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
    leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
  }


  leg->Draw();

  leg = new TLegend(0.2,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();

  leg = new TLegend(0.2,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderB);
  leg->Draw();

if (DATA && MC)
  {
      pad10->cd();
      
      // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
      hDataMC->Divide(htotData,htotMC,1,1);
      hDataMC->Draw("E"); 
      hDataMC->SetLineColor(1);
      hDataMC->SetLineStyle(1);
      hDataMC->SetLineWidth(1);
      hDataMC->SetMarkerColor(kBlack);
      hDataMC->SetMarkerStyle(20);
      hDataMC->SetMarkerSize(0.6);
      hDataMC->SetTickLength(0.10, "X"); hDataMC->SetTickLength(0.05, "YZ");
      hDataMC->SetLabelOffset(0.02,"X");
      hDataMC->SetLabelOffset(0.02,"Y");
      hDataMC->SetLabelSize(0.12, "XY");
      hDataMC->SetLabelFont(42, "XYZ"); 
      hDataMC->SetTitleFont(42, "XYZ");
      hDataMC->SetTitleSize(0.14, "XYZ"); 
      hDataMC->SetTitleOffset(0.9,"X");
      hDataMC->SetTitleOffset(0.5,"Y");
      hDataMC->GetXaxis()->SetTitle(xtitle);
      hDataMC->GetXaxis()->SetTitleColor(1);
      hDataMC->GetXaxis()->SetNdivisions(509);
      hDataMC->GetYaxis()->SetTitle("Data / Sim.");
      hDataMC->GetYaxis()->SetTitleColor(1);
      hDataMC->GetYaxis()->SetNdivisions(509);
      hDataMC->SetNdivisions(509,"XYZ");
      hDataMC->SetMinimum(0.5); 
      hDataMC->SetMaximum(1.5);


      const TH1F* hRatio_constC(hDataMC);

      hRatioUp = new TGraphAsymmErrors(hRatio_constC);
      hRatioDown = new TGraphAsymmErrors(hRatio_constC);

      for (unsigned int b = 1 ; b < hDataMC->GetNbinsX(); b++)
      {
        float htotMCbin = htotMC->GetBinContent(b);
        float htotMCbinErrorStat = htotMC->GetBinError(b);

        float bin =  hDataMC->GetBinContent(b);
        float binStatError =  hDataMC->GetBinError(b);

        float binSysErrorUp = hSumQuadratic_Up->GetBinContent(b);
        float binSysErrorDown = hSumQuadratic_Down->GetBinContent(b);
        // std::cout<<"bin:     "<<bin<<" binStatError: "<<binStatError<<" hSumQuadratic_Up->GetBinContent(b): "<<hSumQuadratic_Up->GetBinContent(b)<<" hSumQuadratic_Down->GetBinContent(b): "<<hSumQuadratic_Down->GetBinContent(b)<<std::endl;
        // std::cout<<"bin "<<bin<<" htotMC->GetBinError(b) : "<<htotMC->GetBinError(b)<<" and hDataMC->GetBinError(b): "<<hDataMC->GetBinError(b)<<std::endl;
        float errorTotalUp = sqrt(binSysErrorUp*binSysErrorUp + binStatError*binStatError);// unccorrelated SYS uncertainties
        float errorTotalDown = sqrt(binSysErrorDown*binSysErrorDown + binStatError*binStatError);// unccorrelated SYS uncertainties
        // std::cout<<"bin: "<<bin<<" binStatError: "<<binStatError<<" binSysErrorUp: "<<binSysErrorUp<<" binSysErrorDown: "<<binSysErrorDown<<" errorTotalUp: "<<errorTotalUp<<" errorTotalDown: "<<errorTotalDown<<std::endl;
        float xbin = hDataMC->GetBinCenter(b);
        hRatioUp->SetPoint(b,xbin,1);
        hRatioDown->SetPoint(b,xbin,1);

        hRatioUp->SetPointError(b,0,0,0,errorTotalUp);
        hRatioDown->SetPointError(b,0,0,errorTotalDown,0);
        // hRatioUp->SetBinContent(b, bin);
        // hRatioDown->SetBinContent(b, bin);

        // hRatioUp->SetBinError(b, errorTotalUp);
        // hRatioDown->SetBinError(b, errorTotalDown);

      }

      hRatioUp->Draw("E3same");//HISTSAME
      hRatioUp->SetFillColor(kGray+2);
      hRatioUp->SetFillStyle(3001);
      hRatioUp->SetFillColorAlpha(kGray+2, 0.9);


        hRatioUp->SetLineStyle(1);
        hRatioUp->SetLineWidth(1);
        // hRatioUp->SetLineColor(ColorBlue);
        // hRatioUp->SetLineColorAlpha(ColorBlue,0.8);

      hRatioDown->Draw("E3same");//HISTSAME
      // hRatioDown->SetLineColor(ColorRed);
      // hRatioDown->SetLineColorAlpha(ColorRed,0.8);
      hRatioDown->SetLineStyle(1);
      hRatioDown->SetLineWidth(1);
      hRatioDown->SetFillColorAlpha(kGray+2, 0.9);
      hRatioDown->SetFillColor(kGray+2);
      hRatioDown->SetFillStyle(3001); 
  }

// *****************************************************************************

 pad4->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
  // !! ------------------------------- !!//
 hSumQuadratic_Up = new TH1F("hSumQuadratic_Up","",nbin,xmin,xmax);
 hSumQuadratic_Down = new TH1F("hSumQuadratic_Down","",nbin,xmin,xmax);

  for (unsigned int i_up = 0 ; i_up < f_SYST_Up.size(); i_up++)
    {

        f_SYST_Up[i_up]->cd();
        TH1F* h_Temp = (TH1F*)gROOT->FindObject(htitleD+"_"+SYSTNAME_UP[i_up]);
        for (int bin = 1; bin <= h_Temp->GetNbinsX(); ++bin) {
            double val = h_Temp->GetBinContent(bin)-1; // !! rescale around 0 ;
            double sumVal = hSumQuadratic_Up->GetBinContent(bin);
            // std::cout<< "valUp = " << val << " sumValUpb4 = " << sumVal << std::endl;
            hSumQuadratic_Up->SetBinContent(bin, sumVal + val * val);
            // std::cout<< "valUp = " << val << " sumValUpafter = " << hSumQuadratic_Up->GetBinContent(bin) << std::endl;
        }

        f_SYST_Down[i_up]->cd();
        h_Temp = (TH1F*)gROOT->FindObject(htitleD+"_"+SYSTNAME_DOWN[i_up]);
        for (int bin = 1; bin <= h_Temp->GetNbinsX(); ++bin) {
            double val = h_Temp->GetBinContent(bin) -1 ; // !! rescale around 0 ;
            double sumVal = hSumQuadratic_Down->GetBinContent(bin);
            hSumQuadratic_Down->SetBinContent(bin, sumVal + val * val);
        }        
    }

    if (hSumQuadratic_Up) {
        for (int bin = 1; bin <= hSumQuadratic_Up->GetNbinsX(); ++bin) {
            hSumQuadratic_Up->SetBinContent(bin, sqrt(hSumQuadratic_Up->GetBinContent(bin)));
            // std::cout<< "sumValUp = " << hSumQuadratic_Up->GetBinContent(bin) << std::endl;
        }
    }
    if (hSumQuadratic_Down) {
        for (int bin = 1; bin <= hSumQuadratic_Down->GetNbinsX(); ++bin) {
            hSumQuadratic_Down->SetBinContent(bin, sqrt(hSumQuadratic_Down->GetBinContent(bin)));
        }
    }
  // !! ---------------------------------------

// ***************************************************************
if (MC)
  {
    f1_DY->cd();
    g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleD);
    h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
    h_DY->Add(g1_DY, h_DY, 1,0);

    f2_DY->cd();
    g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleD);
    h_DY->Add(g2_DY, h_DY, 1, 1);

    f1_VV->cd();
    g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleD);
    h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
    h_VV->Add(g1_VV, h_VV, 1,0);

    f2_VV->cd();
    g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleD);
    h_VV->Add(g2_VV, h_VV, 1, 1);

    f3_VV->cd();
    g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleD);
    h_VV->Add(g3_VV, h_VV, 1, 1);

    f1_TTV->cd();
    g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleD);
    h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
    h_TTV->Add(g1_TTV, h_TTV, 1,0);

    f2_TTV->cd();
    g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleD);
    h_TTV->Add(g2_TTV, h_TTV, 1, 1);

    f3_TTV->cd();
    g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleD);
    h_TTV->Add(g3_TTV, h_TTV, 1, 1);

    h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
    f1_ST->cd();
    g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleD);
    h_ST->Add(g1_ST, h_ST, 1,0);

    f2_ST->cd();
    g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleD);
    h_ST->Add(g2_ST, h_ST, 1, 1);

      f3_ST->cd();
    g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleD);
    h_ST->Add(g3_ST, h_ST, 1, 1);

      f4_ST->cd();
    g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleD);
    h_ST->Add(g4_ST, h_ST, 1, 1);

    f1_TT->cd();
    g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleD);

    h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
    h_TT->Add(g1_TT, h_TT, rwTT*1,0);

    f2_TT->cd();
    g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleD);
    h_TT->Add(g2_TT, h_TT, 1,1);

    //!! --------------------------
    h_DY->Scale(1*scaleMC);
    h_VV->Scale(1*scaleMC);
    h_TTV->Scale(1*scaleMC);
    h_ST->Scale(1*scaleMC);
    h_TT->Scale(1*scaleMC);
    // Daniel 
    h_VV->Add(h_DY, h_VV, 1, 1);
    h_ST->Add(h_ST, h_VV, 1, 1);
    h_TTV->Add(h_TTV, h_ST, 1, 1);
    htotMC->Add(h_TTV, htotMC, 1, 0);
    htotMC->Add(htotMC,h_TT, 1, 1);
    // -----


    htotMC->Draw("HE"); 
    htotMC->SetFillStyle(1001);
    htotMC->SetFillColorAlpha(ColorRed, 1);
    htotMC->SetLineColor(ColorRed);
    htotMC->SetLineColorAlpha(ColorRed, 1);
    htotMC->SetLineStyle(1);
    htotMC->SetLineWidth(1);
    htotMC->SetTickLength(0.03, "YZ");
    htotMC->SetTickLength(0.03,"X");
    htotMC->SetLabelOffset(0.015,"X");
    htotMC->SetLabelOffset(0.007,"Y");
    htotMC->SetLabelSize(0.045, "XYZ");
    htotMC->SetLabelFont(42, "XYZ"); 
    htotMC->SetTitleSize(0.055, "XYZ"); 
    htotMC->SetTitleFont(42, "XYZ");
    htotMC->SetTitleOffset(1.2,"X"); 
    htotMC->SetTitleOffset(1.3,"Y");
    htotMC->GetXaxis()->SetTitle(xtitle);
    htotMC->GetXaxis()->SetTitleColor(1);
    htotMC->GetYaxis()->SetTitle(ytitle);
    htotMC->GetYaxis()->SetTitleColor(1);
    htotMC->SetNdivisions(509,"XYZ");

    if (logy)
      {
        htotMC->SetMinimum(1); 
    htotMC->SetMaximum(htotMC->GetMaximum()*1000); 
      }
    else 
      {
        htotMC->SetMinimum(1); 
        htotMC->SetMaximum(htotMC->GetMaximum()*2); 
      }


    h_ST->Draw("HEsame"); 
    h_ST->SetFillColorAlpha(ColorDarkPurple, 1);
    h_ST->SetLineColor(ColorDarkPurple);
    h_ST->SetLineColorAlpha(ColorDarkPurple, 1);
    h_ST->SetLineStyle(1);
    h_ST->SetLineWidth(3);

    h_VV->Draw("HEsame"); 
    h_VV->SetFillColorAlpha(ColorOrange, 1);
    h_VV->SetLineColor(ColorOrange);
    h_VV->SetLineColorAlpha(ColorOrange, 1);
    h_VV->SetLineStyle(1);
    h_VV->SetLineWidth(3);

      h_DY->Draw("HEsame"); 
    h_DY->SetFillColorAlpha(ColorBlue,1);
    h_DY->SetLineColor(ColorBlue);
    h_DY->SetLineColorAlpha(ColorBlue,1);
    h_DY->SetLineStyle(1);
    h_DY->SetLineWidth(3);
    h_DY->SetTickLength(0.03, "YZ");
    h_DY->SetTickLength(0.03,"X");
  }

if (DATA)
  {
    f1_Data_emu->cd();
    g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleD);
    htotData->Add(g1_Data_emu, htotData, 1, 0);

    htotData->Draw("PE1same");
    htotData->SetMarkerStyle(20);
    htotData->SetMarkerSize(1);
    htotData->SetMarkerColor(kBlack);
    htotData->SetLineColor(kBlack);
    htotData->SetLineWidth(1);
        if (DATA && !MC)
      {
        htotData->SetTickLength(0.03, "YZ");
        htotData->SetTickLength(0.03,"X");
        htotData->SetLabelOffset(0.015,"X");
        htotData->SetLabelOffset(0.007,"Y");
        htotData->SetLabelSize(0.045, "XYZ");
        htotData->SetLabelFont(42, "XYZ"); 
        htotData->SetTitleSize(0.055, "XYZ"); 
        htotData->SetTitleFont(42, "XYZ");
        htotData->SetTitleOffset(1.2,"X"); 
        htotData->SetTitleOffset(1.3,"Y");
            htotData->GetXaxis()->SetTitle(xtitle);
        htotData->GetXaxis()->SetTitleColor(1);
        htotData->GetYaxis()->SetTitle(ytitle);
        htotData->GetYaxis()->SetTitleColor(1);
                if (logy)
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*100); 
          }
        else 
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*2); 
          }
      }
  }

if (Signal)
  {

    f1_LLP->cd();
      g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleD);
    g1_LLP->Sumw2();
    h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
    h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

    f2_LLP->cd();
    g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleD);
    g2_LLP->Sumw2();
    h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
    h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

    f3_LLP->cd();
    
    g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleD);
    g3_LLP->Sumw2();
    h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
    h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

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

  }


//  hsolve->Multiply(hsolve, htotData, 1., 1.);
   if (DATA )
  {
    hsolve->Multiply(hsolve, htotData, 1., 1.);
  }
 else if (MC && !DATA)
  {
    hsolve->Multiply(hsolve, htotMC, 1., 1.);
  }

  leg = new TLegend(0.17,0.94,0.50,0.98);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  
  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
  if (DATA)
    {
      leg->AddEntry(htotData, " e#mu data","PE1");
    }
  if (MC)
    {
      leg->AddEntry(htotMC, " t#bar{t} + X","F");
      leg->AddEntry(h_ST, " Single t","F");
      leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
      leg->AddEntry(h_DY, " DY","F");
    }
  if (Signal)
    {
      leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
      leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
      leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
    }
  leg->Draw();

  leg = new TLegend(0.2,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();

  leg = new TLegend(0.2,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
 leg->SetTextSize(0.06);
  leg->SetHeader(HeaderD);
  leg->Draw();

if (DATA && MC)
  {
      pad11->cd();
      
      // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
      hDataMC->Divide(htotData,htotMC,1,1);
      hDataMC->Draw("E"); 
      hDataMC->SetLineColor(1);
      hDataMC->SetLineStyle(1);
      hDataMC->SetLineWidth(1);
      hDataMC->SetMarkerColor(kBlack);
      hDataMC->SetMarkerStyle(20);
      hDataMC->SetMarkerSize(0.6);
      hDataMC->SetTickLength(0.10, "X"); hDataMC->SetTickLength(0.05, "YZ");
      hDataMC->SetLabelOffset(0.02,"X");
      hDataMC->SetLabelOffset(0.02,"Y");
      hDataMC->SetLabelSize(0.12, "XY");
      hDataMC->SetLabelFont(42, "XYZ"); 
      hDataMC->SetTitleFont(42, "XYZ");
      hDataMC->SetTitleSize(0.14, "XYZ"); 
      hDataMC->SetTitleOffset(0.9,"X");
      hDataMC->SetTitleOffset(0.5,"Y");
      hDataMC->GetXaxis()->SetTitle(xtitle);
      hDataMC->GetXaxis()->SetTitleColor(1);
      hDataMC->GetXaxis()->SetNdivisions(509);
      hDataMC->GetYaxis()->SetTitle("Data / Sim.");
      hDataMC->GetYaxis()->SetTitleColor(1);
      hDataMC->GetYaxis()->SetNdivisions(509);
      hDataMC->SetNdivisions(509,"XYZ");
      hDataMC->SetMinimum(0.5); 
      hDataMC->SetMaximum(1.5);


      const TH1F* hRatio_constD(hDataMC);

      hRatioUp = new TGraphAsymmErrors(hRatio_constD);
      hRatioDown = new TGraphAsymmErrors(hRatio_constD);

      for (unsigned int b = 1 ; b < hDataMC->GetNbinsX(); b++)
      {
        float htotMCbin = htotMC->GetBinContent(b);
        float htotMCbinErrorStat = htotMC->GetBinError(b);

        float bin =  hDataMC->GetBinContent(b);
        float binStatError =  hDataMC->GetBinError(b);

        float binSysErrorUp = hSumQuadratic_Up->GetBinContent(b);
        float binSysErrorDown = hSumQuadratic_Down->GetBinContent(b);
        // std::cout<<"bin:     "<<bin<<" binStatError: "<<binStatError<<" hSumQuadratic_Up->GetBinContent(b): "<<hSumQuadratic_Up->GetBinContent(b)<<" hSumQuadratic_Down->GetBinContent(b): "<<hSumQuadratic_Down->GetBinContent(b)<<std::endl;
        // std::cout<<"bin "<<bin<<" htotMC->GetBinError(b) : "<<htotMC->GetBinError(b)<<" and hDataMC->GetBinError(b): "<<hDataMC->GetBinError(b)<<std::endl;
        float errorTotalUp = sqrt(binSysErrorUp*binSysErrorUp + binStatError*binStatError);// unccorrelated SYS uncertainties
        float errorTotalDown = sqrt(binSysErrorDown*binSysErrorDown + binStatError*binStatError);// unccorrelated SYS uncertainties
        // std::cout<<"bin: "<<bin<<" binStatError: "<<binStatError<<" binSysErrorUp: "<<binSysErrorUp<<" binSysErrorDown: "<<binSysErrorDown<<" errorTotalUp: "<<errorTotalUp<<" errorTotalDown: "<<errorTotalDown<<std::endl;
        float xbin = hDataMC->GetBinCenter(b);
        hRatioUp->SetPoint(b,xbin,1);
        hRatioDown->SetPoint(b,xbin,1);

        hRatioUp->SetPointError(b,0,0,0,errorTotalUp);
        hRatioDown->SetPointError(b,0,0,errorTotalDown,0);
        // hRatioUp->SetBinContent(b, bin);
        // hRatioDown->SetBinContent(b, bin);

        // hRatioUp->SetBinError(b, errorTotalUp);
        // hRatioDown->SetBinError(b, errorTotalDown);

      }

      hRatioUp->Draw("E3same");//HISTSAME
      hRatioUp->SetFillColor(kGray+2);
      hRatioUp->SetFillStyle(3001);
      hRatioUp->SetFillColorAlpha(kGray+2, 0.9);


        hRatioUp->SetLineStyle(1);
        hRatioUp->SetLineWidth(1);
        // hRatioUp->SetLineColor(ColorBlue);
        // hRatioUp->SetLineColorAlpha(ColorBlue,0.8);

      hRatioDown->Draw("E3same");//HISTSAME
      // hRatioDown->SetLineColor(ColorRed);
      // hRatioDown->SetLineColorAlpha(ColorRed,0.8);
      hRatioDown->SetLineStyle(1);
      hRatioDown->SetLineWidth(1);
      hRatioDown->SetFillColorAlpha(kGray+2, 0.9);
      hRatioDown->SetFillColor(kGray+2);
      hRatioDown->SetFillStyle(3001); 

      // //!!
  }

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

//   hsolve->SetLineColor(kBlack);
//  hsolve->SetFillColor(kBlack);
//  hsolve->SetFillStyle(3004);
// //  hsolve->SetLineStyle(1);
// //  hsolve->SetLineWidth(2);

 if ( DATA && !MC ) {
   hsolve->SetLineColor(kAzure+6);
   hsolve->SetFillColor(kAzure+6);
 }
 else {
   hsolve->SetLineColor(kAzure+6);
   hsolve->SetFillColor(kAzure+6);
 }
 hsolve->SetFillStyle(3001);
 hsolve->SetLineStyle(1);
 hsolve->SetLineWidth(1);
 hsolve->Draw("sameE2"); 

  // hsolve->Draw("E4same");
  hsolve->SaveAs("Tight_1Vtx_"+Plots+".root");
  c1->cd();
  pad2->cd();
  hsolve->Draw("E4same");
  TString namele =  Name+"_"+Plots+"_v2";
  if (DATA && !MC) namele += "_Data";
  if (MC && !DATA) namele += "_MC";
  if (MC && DATA) namele += "_DataMC";

  c2->SaveAs(namele+".pdf");

  return c1;
}