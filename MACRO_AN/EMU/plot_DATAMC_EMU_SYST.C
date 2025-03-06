#include <iostream>
#include <TROOT.h>
#include "TH1.h"
#include "TColor.h"

void plot(int method, TString Prod, TString Name, TString Year, TString Dmode, TString Plots)
{
int stati=0;
bool fit= 1;
bool logy=0;


// number of vertices:
//$$
  int nvtx = 2;

//$$
  float hmin = 0.5; // cannot be 0 for logy=1
//$$  float hmax = 1E5;	         // for eta<2.4 pt>80 or CRlowpt
  float hmax = 1E6;	         // for eta<2.4 pt>80 or CRlowpt
  float hmaxBD = 1E6;	         // for eta<2.4 pt>80 or CRlowpt

  if ( nvtx == 2 ) {
    hmax = 10E9;   // for eta<2.4 pt>80  
    hmaxBD = 1E6;
//     hmax = 1E6;   // for eta<2.4 pt>80  
//     hmaxBD = 1E9;
  }

 TString Yearcor = Year;
 if (Year == "2016PRE") Yearcor = "2016preVFP";
 if (Year == "2016POST") Yearcor = "2016";
//MUMU
TString extension = "_SYST";
TString COMPARE = "";
TFile* f1_Data_emu  = new TFile("../../DATA_EMU_"+Year+"_31_10_2024/DATAMC_"+COMPARE+"MuonEG-UL2018_MiniAODv2_GT36-v1.root");

 TFile* f1_DY  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f2_DY  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f1_TT  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");//DATAMC_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_woTopPt ,DATAMC_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8
 TFile* f2_TT  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_ST  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_ST  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f3_ST  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f4_ST  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_TTV = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8.root");
 TFile* f2_TTV = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f3_TTV = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"TTWW_TuneCP5_13TeV-madgraph-pythia8.root");
 TFile* f1_VV  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_VV  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 TFile* f3_VV  = new TFile("../../MC_EMU_03_02_2025/DATAMC_"+COMPARE+"ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");


 std::vector<TFile*> f_SYST_Up;
 std::vector<TFile*> f_SYST_Down;

  bool WLepton = true;
  bool WGen = false;
  std::vector<TString> SYSTNAME_UP = {"LumiUp","L1Up","TriggerUp","PUUp","SFEleUp","TopPtUp","JECUp","JERUp"}; //,,

  std::vector<TString> SYSTNAME_DOWN = {"LumiDown","L1Down","TriggerDown","PUDown","SFEleDown","TopPtDown","JECDown","JERDown"};//,,"RoccorDown",

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
f_SYST_Up.push_back(f_SFEleUp);
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
f_SYST_Down.push_back(f_SFEleDown);
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


 TString DATAFILE[1] = {"MuonEG-UL2018_MiniAODv2_GT36-v1_"
};
if (Year == "2018") DATAFILE[0] = "MuonEG-UL2018_MiniAODv2_GT36-v1_";


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


 TString ytitle = "Events";
 TString xtitle = "vertex BDT score"; 
 int nbin = 26; 
 float xmin = -1.04;
 float xmax =  1.04;
//$$
 TString HeaderCMS = "CMS";

 if (Year == "2016") HeaderCMS = "2016                                      36.3 fb^{-1} (13 TeV)";
 if (Year == "2017") HeaderCMS = "2017                                      41.5 fb^{-1} (13 TeV)";
 if (Year == "2018") HeaderCMS = "2018                                      59.8 fb^{-1} (13 TeV)";


    TString htitleA = "hData_EVT34_1Vtx_BDTvtx";
    TString htitleB = "hData_EVT34_1Vtx_BDTvtx";
    TString htitleC = "hData_EVT34_1Vtx_BDTvtx";
    TString htitleD = "hData_EVT34_1Vtx_BDTvtx";

    int nbin1 = 25; 
    float xmin1 = -1.0;
    float xmax1 =  1.0;
    int nbin2 = 25; 
    float xmin2 = -1.0;
    float xmax2 =  1.0;
    int nbin3 = 25; 
    float xmin3 = -1.0;
    float xmax3 =  1.0;
    int nbin4 = 25; 
    float xmin4 = -1.0;
    float xmax4 =  1.0;
    TString HeaderA = "A";
    TString HeaderB = "B";
    TString HeaderC = "C";
    TString HeaderD = "D";
    TString HeaderNVtx1 = "k Vtx";
    TString HeaderNVtx2= "k Vtx";
    TString HeaderNVtx3 = "k Vtx";
    TString HeaderNVtx4 = "k Vtx";
    TString xtitle1 = "var";
    TString xtitle2 = "var";
    TString xtitle3 = "var";
    TString xtitle4 = "var";

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
    
 int Method = method;
    if (Method == 0) // !! OK
      {
        
        htitleA = "Tree_filter_";
        nbin1 = 2; 
        xmin1 = 0.;
        xmax1 =  2;

        // HeaderNVtx = "2 Vtx"; 
        xtitle1 = "Filter";


        htitleB = "Tree_Mumu_nosel";
        nbin2 = 150; 
        xmin2 = 0;
        xmax2 = 600;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 VtxAll";
        xtitle2 = "M_{e#mu}";

        htitleC = "Vertices_NoSel";
        nbin3 = 100; 
        xmin3 = 0;
        xmax3 =  100;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 Vtx All"; 
        xtitle3 = "nVtx";

        htitleD = "Vertices_filtercut_";
        nbin4 = 100; 
        xmin4 = 0;
        xmax4 = 100;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 VtxAll";
        xtitle4 = "nVtx_{AfterFilter}";


        HeaderA = "";
        HeaderB = " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" GeV ";
        HeaderC = "";
        HeaderD = "";


      }


    if (Method == 1) // !! OK
      {
        return;
        htitleA = "Dilepton_pt_NoSel";
        nbin1 = 100; 
        xmin1 = 0;
        xmax1 = 300;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 VtxAll";
        xtitle1 = "pt_{ll}";


        htitleB = "Dilepton_eta_NoSel";
        nbin2 = 25; 
        xmin2 = -2.4;
        xmax2 = 2.4;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 VtxAll";
        xtitle2 = "#eta_{ll}";

        htitleC = "Dilepton_phi_NoSel";
        nbin3 = 25; 
        xmin3 = -3.14;
        xmax3 = 3.14;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 VtxAll";
        xtitle3 = "#phi_{ll}";

        htitleD = "Dilepton_mass_NoSel";
        nbin4 = 150; 
        xmin4 = 0;
        xmax4 = 600;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 VtxAll";
        xtitle4 = "M_{ll}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = "";
        HeaderC = "";
        HeaderD = " /"+std::to_string((xmax4-xmin4)/(float)nbin4)+" GeV ";
      }

    if (Method == 2) // !! OK
      {
        htitleA = "leading_muon_pt_reco";
        nbin1 = 100;
        xmin1 = 0;
        xmax1 = 400;
        xtitle1 = "pt_{#mu}";

        htitleB = "leading_lepton_pt_reco";
        nbin2 = 100;
        xmin2 = 0;
        xmax2 = 400;
        xtitle2 = "pt_{e}";

        htitleC = "leading_muon_eta_reco";
        nbin3 = 25;
        xmin3 = -2.4;
        xmax3 = 2.4;
        xtitle3 = "#eta_{#mu}";

        htitleD = "leading_lepton_eta_reco";
        nbin4 = 25;
        xmin4 = -2.4;
        xmax4 = 2.4;
        xtitle4 = "#eta_{e}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" GeV ";
        HeaderC = "";
        HeaderD = "";
      }

    if (Method == 3) // !! OK
      {
        htitleA = "leading_muon_dxy_reco";
        nbin1 = 40;
        xmin1 = -0.2;
        xmax1 = 0.2;
        xtitle1 = "dxy_{#mu}";

        htitleB = "leading_lepton_dxy_reco";
        nbin2 = 200;
        xmin2 = -1;
        xmax2 = 1;
        xtitle2 = "dxy_{e}";

        htitleC = "leading_muon_dz_reco";
        nbin3 = 50;
        xmin3 = -0.5;
        xmax3 = 0.5;
        xtitle3 = "dz_{#mu}";

        htitleD = "leading_lepton_dz_reco";
        nbin4 = 200;
        xmin4 = -1;
        xmax4 = 1;
        xtitle4 = "dz_{e}";

        logy = 1;

        HeaderA = "";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
      }


    if (Method == 4) // !! OK
      {
        htitleA = "hData_jet_pt_";
        nbin1 = 100;
        xmin1 = 0;
        xmax1 = 500;
        xtitle1 = "jet_{pt}";

        htitleB = "hData_jet_eta_";
        nbin2 = 55;
        xmin2 = -5;
        xmax2 = 5;
        xtitle2 = "jet_{#eta}";


        htitleC = "njet_NoSel";
        nbin3 = 20;
        xmin3 = 0;
        xmax3 = 20;
        xtitle3 = "N_{jet}";

        htitleD = "njetNOmu_NoSel";
        nbin4 = 20;
        xmin4 = 0;
        xmax4 = 20;
        xtitle4 = "N_{jet NoLepton}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
      }

    if (Method == 5) // !! OK
      {


        htitleA = "Hemi_pt_";
        nbin1 = 100;
        xmin1 = 0;
        xmax1 = 600;
        xtitle1 = "Hemi_{pt}";

        htitleB = "Hemi_eta_";
        nbin2 = 55;
        xmin2 = -2.5;
        xmax2 = 2.5;
        xtitle2 = "Hemi_{#eta}";


        htitleC = "Hemisphere_leadingpt_";
        nbin3 = 100;
        xmin3 = 0;
        xmax3 = 500;
        xtitle3 = "Hemi_{pt_{1}}";

        htitleD = "Hemisphere_subleadingpt_";
        nbin4 = 100;
        xmin4 = 0;
        xmax4 = 500;
        xtitle4 = "Hemi_{pt_{2}}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = "";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = " /"+std::to_string((xmax4-xmin4)/(float)nbin4)+" GeV ";
      }

// !!----------------------
    if (Method == 6) // !! OK
      {
        htitleA = "leading_jet_pt_";
        nbin1 = 100;
        xmin1 = 0;
        xmax1 = 500;
        xtitle1 = "jet^{1}_{pt}";

        htitleB = "leading_jet_eta_";
        nbin2 = 55;
        xmin2 = -2.5;
        xmax2 = 2.5;
        xtitle2 = "jet^{1}_{#eta}";

        htitleC = "subleading_jet_pt_";
        nbin3 = 100;
        xmin3 = 0;
        xmax3 = 500;
        xtitle3 = "jet^{2}_{pt}";

        htitleD = "subleading_jet_eta_";
        nbin4 = 55;
        xmin4 = -2.5;
        xmax4 = 2.5;
        xtitle4 = "jet^{2}_{#eta}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = "";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = "";
      }


// !!----------------------
    if (Method == 7) // !! OK
      {
        htitleA = "muon_pt_NoSel";
        nbin1 = 100;
        xmin1 = 0;
        xmax1 = 300;
        xtitle1 = "#pt_{#mu}";

        htitleB = "muon_eta_NoSel";
        nbin2 = 25;
        xmin2 = -2.4;
        xmax2 = 2.4;
        xtitle2 = "#eta_{#mu}";

        htitleC = "electron_pt_NoSel";
        nbin3 = 100;
        xmin3 = 0;
        xmax3 = 300;
        xtitle3 = "pt_{e}";

        htitleD = "electron_eta_NoSel";
        nbin4 = 25;
        xmin4 = -2.4;
        xmax4 = 2.4;
        xtitle4 = "#eta_{e}";

        HeaderA = "";
        HeaderB = " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" GeV ";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = "";
      }

// !!----------------------
    if (Method == 8)// !! OK
      {
        htitleA = "LeadingLeptons_dR_";
        nbin1 = 50;
        xmin1 = 0;
        xmax1 = 5;
        xtitle1 = "#Delta R_{l_{1}l_{2}}";

        htitleB = "LeadingLeptons_dPhi_";
        nbin2 = 35;
        xmin2 = 0;
        xmax2 = 3.5;
        xtitle2 = "#Delta #Phi_{l_{1}l_{2}}";

        htitleC = "LeadingJets_dR_";
        nbin3 = 50;
        xmin3 = 0;
        xmax3 = 5;
        xtitle3 = "#Delta R_{j_{1}j_{2}}";

        htitleD = "LeadingJets_dPhi_";
        nbin4 = 35;
        xmin4 = 0;
        xmax4 = 3.5;
        xtitle4 = "#Delta #Phi_{j_{1}j_{2}}";

                HeaderA = std::to_string((xmax1-xmin1)/(float)nbin1);
        HeaderB = std::to_string((xmax2-xmin2)/(float)nbin2);
        HeaderC = std::to_string((xmax3-xmin3)/(float)nbin3);
        HeaderD = std::to_string((xmax4-xmin4)/(float)nbin4);
      }


    // !!----------------------
    if (Method == 9) // !! OK
      {
        htitleA = "LeadingLeptonJet_dRmax_";
        nbin1 = 50;
        xmin1 = 0;
        xmax1 = 5;
        xtitle1 = "#Delta R_{l_{1}j_{1}}";

        htitleB = "LeadingLeptonJet_dRmin_";
        nbin2 = 50;
        xmin2 = 0;
        xmax2 = 5;
        xtitle2 = "#Delta Phi_{l_{1}j_{1}}";

        htitleC = "HemiAxis_Mu_dR_";
        nbin3 = 50;
        xmin3 = 0;
        xmax3 = 5;
        xtitle3 = "#Delta R_{Hemi-Mu}";

        htitleD = "HemiAxis_OpMu_dR_";
        nbin4 = 50;
        xmin4 = 0;
        xmax4 = 5;
        xtitle4 = "#Delta R_{Hemi-OpMu}";

                HeaderA = std::to_string((xmax1-xmin1)/(float)nbin1);
        HeaderB = std::to_string((xmax2-xmin2)/(float)nbin2);
        HeaderC = std::to_string((xmax3-xmin3)/(float)nbin3);
        HeaderD = std::to_string((xmax4-xmin4)/(float)nbin4);
        
      }
    // !!----------------------
    if (Method == 10) // !! OK
      {
        htitleA = "HT_";
        nbin1 = 200;
        xmin1 = 0;
        xmax1 = 800;
        xtitle1 = "H_{T}";

        htitleB = "LT_";
        nbin2 = 100;
        xmin2 = 0;
        xmax2 = 500;
        xtitle2 = "L_{T}";

        htitleC = "nTrks_";
        nbin3 = 100;
        xmin3 = 0;
        xmax3 = 100;
        xtitle3 = "N_{TRK}";

        htitleD = "nLostTracks_";
        nbin4 = 25;
        xmin4 = 0;
        xmax4 = 25;
        xtitle4 = "N_{LostTRK}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" GeV ";
        HeaderC = "";
        HeaderD = "";
      }

    // !!----------------------
    if (Method == 11)  // !! OK
      {
        htitleA = "Hemi_nTrks_";
        nbin1 = 25;
        xmin1 = 0;
        xmax1 = 25;
        xtitle1 = "Hemi_{nTrks}";

        htitleB = "Hemi_Mass_";
        nbin2 = 50;
        xmin2 = 0;
        xmax2 = 300;
        xtitle2 = "Hemi_{Mass}";

        htitleC = "Hemi_nJet_";
        nbin3 = 10;
        xmin3 = 0;
        xmax3 = 10;
        xtitle3 = "Hemi N_{jet}";

        htitleD = "Hemi_nJetNoMu_";
        nbin4 = 10;
        xmin4 = 0;
        xmax4 = 10;
        xtitle4 = "Hemi N_{jetNoMu}";

        HeaderA = "";
        HeaderB = " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" GeV ";
        HeaderC = "";
        HeaderD = "";
      }

    // !!----------------------
    if (Method == 12) // We don't car eanymore
      {
        htitleA = "HemiMu_pt_";
        nbin1 = 100;
        xmin1 = 0;
        xmax1 = 600;
        xtitle1 = "HemiMu_{pt}";

        htitleB = "HemiMu_Mass_";
        nbin2 = 100;
        xmin2 = 0;
        xmax2 = 600;
        xtitle2 = "HemiMu_{mass}";

        htitleC = "Hemi_nTrks_";
        nbin3 = 25;
        xmin3 = 0;
        xmax3 = 25;
        xtitle3 = "Hemi_{nTrks}";

        htitleD = "Hemi_Mass_";
        nbin4 = 50;
        xmin4 = 0;
        xmax4 = 300;
        xtitle4 = "Hemi_{Mass}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" GeV ";
        HeaderC = "";
        HeaderD = " /"+std::to_string((xmax4-xmin4)/(float)nbin4)+" GeV ";
      }
          // !!----------------------
    if (Method == 13) // !! OK 
      {
        htitleA = "K0_mass_";
        nbin1 = 202;
        xmin1 = 0.42;
        xmax1 = 0.58;
        xtitle1 = "K0_{mass}";

        htitleB = "K0_pt_";
        nbin2 = 50;
        xmin2 = 0;
        xmax2 = 50;
        xtitle2 = "K0_{pt}";

        htitleC = "Reco_K0_mass_";
        nbin3 = 202;
        xmin3 = 0.42;
        xmax3 = 0.58;
        xtitle3 = "RecoK0_{mass}";

        htitleD = "Reco_K0_pt_";
        nbin4 = 50;
        xmin4 = 0;
        xmax4 = 50;
        xtitle4 = "RecoK0_{pt}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" GeV ";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = " /"+std::to_string((xmax4-xmin4)/(float)nbin4)+" GeV ";
      }
  // !!----------------------
    if (Method == 14) // !! OK
      {
        htitleA = "L0_mass_";
        nbin1 = 202;
        xmin1 = 1.06;
        xmax1 = 1.18;
        xtitle1 = "L0_{mass}";

        htitleB = "L0_pt_";
        nbin2 = 500;
        xmin2 = 0;
        xmax2 = 50;
        xtitle2 = "L0_{pt}";

        htitleC = "Reco_L0_mass_";
        nbin3 = 202;
        xmin3 = 1.06;
        xmax3 = 1.18;
        xtitle3 = "RecoL0_{mass}";

        htitleD = "Reco_L0_pt_";
        nbin4 = 50;
        xmin4 = 0;
        xmax4 = 50;
        xtitle4 = "RecoL0_{pt}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" GeV ";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = " /"+std::to_string((xmax4-xmin4)/(float)nbin4)+" GeV ";
      }

  // !!----------------------
    if (Method == 15) // !! OK
      {
        htitleA = "SecInt_mass_Selec";
        nbin1 = 20;
        xmin1 = 0;
        xmax1 = 2;
        xtitle1 = "SecInt_{mass}";

        htitleB = "SecInt_drSig_Selec";
        nbin2 = 200;
        xmin2 = 0;
        xmax2 = 2000;
        xtitle2 = "SecInt_{drSig}";

        htitleC = "SecInt_pt_Selec";
        nbin3 = 100;
        xmin3 = 0;
        xmax3 = 100;
        xtitle3 = "SecInt_{pt}";

        htitleD = "SecInt_dzSig_Selec";
        nbin4 = 200;
        xmin4 = 0;
        xmax4 = 2000;
        xtitle4 = "SecInt_dzSig";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = "";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = "";
      }
 
        // !!----------------------
    if (Method == 16) // !! OK
      {
        htitleA = "SecInt_mass_TrackerMatched";
        nbin1 = 20;
        xmin1 = 0;
        xmax1 = 2;
        xtitle1 = "SecInt_{mass}";

        htitleB = "SecInt_drSig_TrackerMatched";
        nbin2 = 200;
        xmin2 = 0;
        xmax2 = 2000;
        xtitle2 = "SecInt IP_{dxy}";

        htitleC = "SecInt_pt_TrackerMatched";
        nbin3 = 100;
        xmin3 = 0;
        xmax3 = 100;
        xtitle3 = "SecInt_{pt}";

        htitleD = "SecInt_dzSig_TrackerMatched";
        nbin4 = 200;
        xmin4 = 0;
        xmax4 = 2000;
        xtitle4 = "SecInt IP_{dz}";       
        
        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = "";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = "";
      }                                                      
                                                            
        // !!----------------------
    if (Method == 17) // !! OK
      {
        htitleA = "Vtx_NChi2_";
        nbin1 = 15;
        xmin1 = 0;
        xmax1 = 15;
        xtitle1 = "Vtx_{#chi^{2}}";

        htitleB = "Vtx_nTrks_";
        nbin2 = 40;
        xmin2 = 0;
        xmax2 = 40;
        xtitle2 = "Vtx_{nTrks}";

        htitleC = "Vtx_Mass_";
        nbin3 = 10;
        xmin3 = 0;
        xmax3 = 100;
        xtitle3 = "Vtx_{mass}";

        htitleD = "Vtx_Dist_";
        nbin4 = 20;
        xmin4 = 0;
        xmax4 = 100;
        xtitle4 = "Vtx_{dist}";

        HeaderA = "";
        HeaderB = "";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = " /"+std::to_string((xmax4-xmin4)/(float)nbin4)+" cm ";
      }      

        // !!----------------------
    if (Method == 18) // !! OK
      {
        htitleA = "SecVtx_NChi2_";
        nbin1 = 15;
        xmin1 = 0;
        xmax1 = 15;
        xtitle1 = "SecVtx_{#chi^{2}}";

        htitleB = "SecVtx_nTrks_";
        nbin2 = 40;
        xmin2 = 0;
        xmax2 = 40;
        xtitle2 = "SecVtx_{nTrks}";

        htitleC = "SecVtx_Mass_";
        nbin3 = 150;
        xmin3 = 0;
        xmax3 = 1500;
        xtitle3 = "SecVtx_{mass}";

        htitleD = "SecVtx_Dist_";
        nbin4 = 100;
        xmin4 = 0;
        xmax4 = 100;
        xtitle4 = "SecVtx_{dist}";
        hmax = 1E3;

        HeaderA = "";
        HeaderB = "";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = " /"+std::to_string((xmax4-xmin4)/(float)nbin4)+" cm ";
      }   

        // !!----------------------
    if (Method == 19)
      {
        htitleA = "Vtx_Step_";
        nbin1 = 4;
        xmin1 = 1;
        xmax1 = 5;
        xtitle1 = "Vtx_{step}";

        htitleB = "Vtx_dR_";
        nbin2 = 50;
        xmin2 = 0;
        xmax2 = 5;
        xtitle2 = "Vtx_{#Delta R}";

        htitleC = "SecVtx_Step_";
        nbin3 = 4;
        xmin3 = 1;
        xmax3 = 5;
        xtitle3 = "SecVtx_{step}";

        htitleD = "SecVtx_dR_";
        nbin4 = 50;
        xmin4 = 0;
        xmax4 = 5;
        xtitle4 = "SecVtx_{#Delta R}";
        HeaderA = "";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
      }  

        // !!----------------------
    if (Method == 20) // !! OK
      {
        htitleA = "FinalVtx_nTrks_";
        nbin1 = 40;
        xmin1 = 0;
        xmax1 = 40;
        xtitle1 = "FinalVtx_{nTrks}";

        htitleB = "FinalVtx_Step_";
        nbin2 = 4;
        xmin2 = 1;
        xmax2 = 5;
        xtitle2 = "FinalVtx_{step}";

        htitleC = "FinalVtx_Mass_";
        nbin3 = 100;
        xmin3 = 0;
        xmax3 = 500;
        xtitle3 = "FinalVtx_{Mass}";

        htitleD = "FinalVtx_HMass_";
        nbin4 = 100;
        xmin4 = 0;
        xmax4 = 500;
        xtitle4 = "FinalVtx_{Hmass}";

        HeaderA = "";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
      }  
//--------------------------------------------

    if (Method == 21) // !! OK
      {
        htitleA = "track_nHitTIB_TRK";
        nbin1 = 15;
        xmin1 = 0;
        xmax1 = 15;
        xtitle1 = "track_{nHitTIB}";

        htitleB = "track_nHitTOB_TRK";
        nbin2 = 15;
        xmin2 = 0;
        xmax2 = 15;
        xtitle2 = "track_{nHitTOB}";

        htitleC = "track_nHitTEC_TRK";
        nbin3 = 15;
        xmin3 = 0;
        xmax3 = 15;
        xtitle3 = "track_{nHitTEC}";

        htitleD = "track_nHitPXB_TRK";
        nbin4 = 15;
        xmin4 = 0;
        xmax4 = 15;
        xtitle4 = "track_{nHitPXB}";

        HeaderA = "";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
      }


    if (Method == 22) // !! OK
      {
        htitleA = "track_nHitPixel_TRK";
        nbin1 = 15;
        xmin1 = 0;
        xmax1 = 15;
        xtitle1 = "track_{nHitPixel}";

        htitleB = "track_nHitPXF_TRK";
        nbin2 = 15;
        xmin2 = 0;
        xmax2 = 15;
        xtitle2 = "track_{nHitPXF}";

        htitleC = "track_isHitPixel_TRK";
        nbin3 = 3000;
        xmin3 = 0;
        xmax3 = 1500;
        xtitle3 = "track_{isHitPixel}";

        htitleD = "track_nLayers_TRK";
        nbin4 = 80;
        xmin4 = 0;
        xmax4 = 20;
        xtitle4 = "track_{nLayers}";
        HeaderA = "";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
    }

    if  (Method == 23) // !! OK
    {
      htitleA = "track_pt_TRK";
      nbin1 = 300;
      xmin1 = 0;
      xmax1 = 300;
      xtitle1 = "track_{pt}";

      htitleB = "track_eta_TRK";
      nbin2 = 80;
      xmin2 = -4;
      xmax2 = 4;
      xtitle2 = "track_{eta}";

      htitleC = "track_NChi2_TRK";
      nbin3 = 5;
      xmin3 = 0;
      xmax3 = 5;
      xtitle3 = "track_{NChi2}";

      htitleD = "track_nhits_TRK";
      nbin4 = 40;
      xmin4 = 0;
      xmax4 = 40;
      xtitle4 = "track_{nhits}";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" GeV ";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
    }


    if  (Method == 24) // !! OK
    {
      htitleA = "track_ntrk10_TRK";
      nbin1 = 100;
      xmin1 = 0;
      xmax1 = 100;
      xtitle1 = "track_{ntrk10}";

      htitleB = "track_ntrk20_TRK";
      nbin2 = 100;
      xmin2 = 0;
      xmax2 = 100;
      xtitle2 = "track_{ntrk20}";

      htitleC = "track_ntrk30_TRK";
      nbin3 = 100;
      xmin3 = 0;
      xmax3 = 100;
      xtitle3 = "track_{ntrk30}";

      htitleD = "track_ntrk40_TRK";
      nbin4 = 100;
      xmin4 = 0;
      xmax4 = 100;
      xtitle4 = "track_{ntrk40}";

        HeaderA = "";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
    }

    if  (Method == 25) // !! OK
    {
      htitleA = "track_Hemi_TRK";//empty => normal
      nbin1 = 6;
      xmin1 = -0.5;
      xmax1 = 5.5;
      xtitle1 = "track_{Hemi}";

      htitleB = "track_lost_TRK";
      nbin2 = 2;
      xmin2 = 0;
      xmax2 = 2;
      xtitle2 = "track_{lost}";

      htitleC = "track_dxy_TRK";
      nbin3 = 100;
      xmin3 = -50;
      xmax3 = 50;
      xtitle3 = "track_{dxy}";

      htitleD = "track_dz_TRK";
      nbin4 = 200;
      xmin4 = -100;
      xmax4 = 100;
      xtitle4 = "track_{dz}";

        HeaderA = "";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
    }


    if  (Method == 26) // !! OK
    {
      htitleA = "track_iJet_TRK";
      nbin1 = 22;
      xmin1 = -2;
      xmax1 = 20;
      xtitle1 = "track_{iJet}";

      htitleB = "track_drSig_TRK";
      nbin2 = 1000;
      xmin2 = 0;
      xmax2 = 1000;
      xtitle2 = "track_{drSig}";

      htitleC = "track_dzSig_TRK";
      nbin3 = 5000;
      xmin3 = 0;
      xmax3 = 5000;
      xtitle3 = "track_{dzSig}";

      htitleD = "track_track_Hemi_dR_TRK";
      nbin4 = 50;
      xmin4 = 0;
      xmax4 = 5;
      xtitle4 = "track_{track_Hemi_dR}";

        HeaderA = "";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
    }

    if  (Method == 27) 
    {
      htitleA = "track_track_Hemi_dRmax_TRK";
      nbin1 = 50;
      xmin1 = 0;
      xmax1 = 5;
      xtitle1 = "track_{track_Hemi_dRmax}";

      htitleB = "track_MVAVal_TRK";
      nbin2 = 101;
      xmin2 = -1;
      xmax2 = 1;
      xtitle2 = "track_{MVAVal}";

      htitleC = "track_Track_firstHit_TRK";
      nbin3 = 4000;
      xmin3 = 0;
      xmax3 = 4000;
      xtitle3 = "track_{Track_firstHit}";

      htitleD = "track_Track_firstHit_TRK";
      nbin4 = 4000;
      xmin4 = 0.;
      xmax4 = 4000;
      xtitle4 = "track_{Track_firstHit}";

        HeaderA = "";
        HeaderB = "";
        HeaderC = "";
        HeaderD = "";
    }

        // !!----------------------
    if (Method == 28) // !! OK
      {
        htitleA = "SecInt_r_Selec";
        nbin1 = 200;
        xmin1 = 0;
        xmax1 = 20;
        xtitle1 = "SecInt r [cm]";

        htitleB = "SecInt_z_Selec";
        nbin2 = 800;
        xmin2 = -200;
        xmax2 = 200;
        xtitle2 = "SecInt z [cm]";

        htitleC = "SecInt_pt_Selec";
        nbin3 = 100;
        xmin3 = 0;
        xmax3 = 100;
        xtitle3 = "SecInt p_{t} [GeV]";

        htitleD = "SecInt_dzSig_Selec";
        nbin4 = 200;
        xmin4 = 0;
        xmax4 = 2000;
        xtitle4 = "SecInt dzSig";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" cm ";
        HeaderB =  " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" cm ";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = "";
      } 
    if (Method == 29) // !! OK
      {
        htitleA = "SecInt_r_TrackerMatched";
        nbin1 = 200;
        xmin1 = 0;
        xmax1 = 20;
        xtitle1 = "SecInt  r [cm]";

        htitleB = "SecInt_z_TrackerMatched";
        nbin2 = 800;
        xmin2 = -200;
        xmax2 = 200;
        xtitle2 = "SecInt z [cm]";

        htitleC = "SecInt_pt_TrackerMatched";
        nbin3 = 100;
        xmin3 = 0;
        xmax3 = 100;
        xtitle3 = "SecInt p_{t} [GeV]";

        htitleD = "SecInt_dzSig_TrackerMatched";
        nbin4 = 200;
        xmin4 = 0;
        xmax4 = 2000;
        xtitle4 = "SecInt dzSig";

        HeaderA = " /"+std::to_string((xmax1-xmin1)/(float)nbin1)+" cm ";
        HeaderB =  " /"+std::to_string((xmax2-xmin2)/(float)nbin2)+" cm ";
        HeaderC = " /"+std::to_string((xmax3-xmin3)/(float)nbin3)+" GeV ";
        HeaderD = "";
      } 

    //-----------------------------------------------------------//
// xsec in rap1

  TLegend* leg;
    
// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",0,0,1600,1500);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);


c1->cd();

TPad* pad1 = new TPad("pad1","This is pad1",0.03,0.62,0.49,0.92,21);
TPad* pad2 = new TPad("pad2","This is pad2",0.51,0.62,0.97,0.92,21);
TPad* pad3 = new TPad("pad3","This is pad3",0.03,0.17,0.49,0.47,21);
TPad* pad4 = new TPad("pad4","This is pad4",0.51,0.17,0.97,0.47,21);

TPad* rap1 = new TPad("rap1","This is rap1",0.03,0.48,0.49,0.61,21);
TPad* rap2 = new TPad("rap2","This is rap2",0.51,0.48,0.97,0.61,21);
TPad* rap3 = new TPad("rap3","This is rap3",0.03,0.03,0.49,0.16,21);
TPad* rap4 = new TPad("rap4","This is rap4",0.51,0.03,0.97,0.16,21);

rap1->SetFillColor(0);
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
pad1->SetTopMargin(0.07);
pad1->SetBottomMargin(0.01);
pad1->SetRightMargin(0.04);
pad1->SetLeftMargin(0.16);

pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(logy);
pad2->SetTopMargin(0.07);
pad2->SetBottomMargin(0.01);
pad2->SetRightMargin(0.04);
pad2->SetLeftMargin(0.16);

pad3->SetFillColor(0);
pad3->SetBorderMode(0);
pad3->SetFrameFillColor(10);
pad3->Draw();
pad3->SetLogy(logy);
pad3->SetTopMargin(0.07);
pad3->SetBottomMargin(0.01);
pad3->SetRightMargin(0.04);
pad3->SetLeftMargin(0.16);

pad4->SetFillColor(0);
pad4->SetBorderMode(0);
pad4->SetFrameFillColor(10);
pad4->Draw();
pad4->SetLogy(logy);
pad4->SetTopMargin(0.07);
pad4->SetBottomMargin(0.01);
pad4->SetRightMargin(0.04);
pad4->SetLeftMargin(0.16);

rap1->SetFillColor(0);
rap1->SetBorderMode(0);
rap1->SetFrameFillColor(10);
rap1->Draw();
rap1->SetLogy(0);
rap1->SetTopMargin(0.00);
rap1->SetBottomMargin(0.30);
rap1->SetRightMargin(0.04);
rap1->SetLeftMargin(0.16);
rap1->SetGrid();

rap2->SetFillColor(0);
rap2->SetBorderMode(0);
rap2->SetFrameFillColor(10);
rap2->Draw();
rap2->SetLogy(0);
rap2->SetTopMargin(0.00);
rap2->SetBottomMargin(0.30);
rap2->SetRightMargin(0.04);
rap2->SetLeftMargin(0.16);
rap2->SetGrid();

rap3->SetFillColor(0);
rap3->SetBorderMode(0);
rap3->SetFrameFillColor(10);
rap3->Draw();
rap3->SetLogy(0);
rap3->SetTopMargin(0.00);
rap3->SetBottomMargin(0.30);
rap3->SetRightMargin(0.04);
rap3->SetLeftMargin(0.16);
rap3->SetGrid();

rap4->SetFillColor(0);
rap4->SetBorderMode(0);
rap4->SetFrameFillColor(10);
rap4->Draw();
rap4->SetLogy(0);
rap4->SetTopMargin(0.00);
rap4->SetBottomMargin(0.30);
rap4->SetRightMargin(0.04);
rap4->SetLeftMargin(0.16);
rap4->SetGrid();

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
gROOT->SetBatch(kTRUE); 



// !! --------------------------- PAD 1 --------------------------------------  !! //
  pad1->cd();
   xtitle = xtitle1;
 nbin   =   nbin1; 
 xmin   =   xmin1;
 xmax   =   xmax1;

 TH1F* g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);//ok
 TH1F* g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);//ok
 TH1F*  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

 
 TH1F* g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);//ok
 TH1F* g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);//ok
 TH1F* g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);//ok
 TH1F*  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 
TH1F* g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);//ok
 TH1F* g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);//ok
 TH1F* g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);//ok
 TH1F* g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);//ok
 TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 TH1F* g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);//ok
 TH1F* g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);//ok
 TH1F*  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 TH1F* g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);//ok
 TH1F* g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);//ok 
 TH1F* g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);//ok
 TH1F*  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
  


// *****************************************************************************

 pad1->cd();

 TH1F* htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 TH1F* htotData  = new TH1F("htotData","",nbin,xmin,xmax);
TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

  htotData->Sumw2();
  htotMC->Sumw2();
  hDataMC->Sumw2();
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
 f1_Data_emu->cd();
 TH1F* g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleA);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//------x);

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 
 f2_DY->cd();
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
//  g2_DY->Sumw2();
//  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g2_DY, h_DY, 1,1);


 f1_VV->cd();
 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);
//  g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);
//  g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);
//  g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);
//  g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);
//  g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);
//  g3_TTV->Sumw2();
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

 float rescale = 0.0954031129 ; //13791474/1.4456*10^{8} = smallnentries/norm * nentries/small nentries
 // = nentries/norm where norm = number of evetn sin Ntuple et nentries = number of events in MiniNtuple

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
//  g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1*rescale,0);

f2_TT->cd();
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);
 h_TT->Add(g2_TT, h_TT, 1,1);
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
 htotMC->SetFillColorAlpha(ColorBlue, 1);
 htotMC->SetLineColor(ColorBlue);
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
 htotMC->GetYaxis()->SetTitle(ytitle + HeaderA );
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
 htotMC->SetMinimum(1); 
 htotMC->SetMaximum(htotMC->GetMaximum()*2); 
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);

 h_VV->Draw("HEsame"); 
 h_VV->SetFillColorAlpha(ColorOrange, 1);
 h_VV->SetLineColor(ColorOrange);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

 h_TTV->Draw("HEsame"); 
 h_TTV->SetFillColorAlpha(ColorNeutral, 1);
 h_TTV->SetLineColor(ColorNeutral);
 h_TTV->SetLineStyle(1);
 h_TTV->SetLineWidth(3);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColorAlpha(ColorDarkPurple, 1);
 h_ST->SetLineColor(ColorDarkPurple);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 h_TT->Draw("HEsame"); 
 h_TT->SetFillColorAlpha(ColorRed, 1);
 h_TT->SetLineColor(ColorRed);
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


  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  leg = new TLegend(0.64,0.60,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " #mu#mu data","PE1");
    leg->AddEntry(h_TT, " t#bar{t}","F");
    leg->AddEntry(htotMC, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");


  leg->Draw();


rap1->cd();

TH1F* hRatio = new TH1F("hRatio","",nbin,xmin,xmax);
hRatio->Sumw2();
hRatio->Divide(htotData, htotMC, 1, 1);
hRatio->Draw("E"); 
hRatio->SetLineColor(1);
hRatio->SetLineStyle(1);
hRatio->SetLineWidth(1);
hRatio->SetMarkerColor(kBlack);
hRatio->SetMarkerStyle(20);
hRatio->SetMarkerSize(0.6);
hRatio->SetTickLength(0.10, "X"); hRatio->SetTickLength(0.05, "YZ");
hRatio->SetLabelOffset(0.02,"X");
hRatio->SetLabelOffset(0.02,"Y");
hRatio->SetLabelSize(0.12, "XY");
hRatio->SetLabelFont(42, "XYZ"); 
hRatio->SetTitleFont(42, "XYZ");
hRatio->SetTitleSize(0.14, "XYZ"); 
hRatio->SetTitleOffset(0.9,"X");
hRatio->SetTitleOffset(0.5,"Y");
hRatio->GetXaxis()->SetTitle(xtitle);
hRatio->GetXaxis()->SetTitleColor(1);
hRatio->GetXaxis()->SetNdivisions(509);
hRatio->GetYaxis()->SetTitle("Data / Sim.");
hRatio->GetYaxis()->SetTitleColor(1);
hRatio->GetYaxis()->SetNdivisions(509);
hRatio->SetNdivisions(509,"XYZ");
hRatio->SetMinimum(0.5); 
hRatio->SetMaximum(1.5); 

const TH1F* hRatio_const(hRatio);

TGraphAsymmErrors *hRatioUp = new TGraphAsymmErrors(hRatio_const);
TGraphAsymmErrors *hRatioDown = new TGraphAsymmErrors(hRatio_const);

for (unsigned int b = 1 ; b < hRatio->GetNbinsX(); b++)
{
  float htotMCbin = htotMC->GetBinContent(b);
  float htotMCbinErrorStat = 0;
  if (htotMC->GetBinError(b) != 0)
    {
      htotMCbinErrorStat = htotMC->GetBinError(b);
    }

  float bin =  hRatio->GetBinContent(b);
  float binStatError =  hRatio->GetBinError(b);

  float htotDATAbinErrorStat  = 0;
  if (htotData->GetBinError(b) != 0)
    {
      htotDATAbinErrorStat = htotData->GetBinError(b);
    }
    float binSysErrorUp = hSumQuadratic_Up->GetBinContent(b);
  float binSysErrorDown = hSumQuadratic_Down->GetBinContent(b);

  float errorTotalUp = sqrt(binSysErrorUp*binSysErrorUp + binStatError*binStatError);// unccorrelated SYS uncertainties htotMCbinErrorStat*htotMCbinErrorStat
  float errorTotalDown = sqrt(binSysErrorDown*binSysErrorDown + binStatError*binStatError);// unccorrelated SYS uncertainties htotMCbinErrorStat*htotMCbinErrorStat

  float xbin = hRatio->GetBinCenter(b);
  hRatioUp->SetPoint(b,xbin,1);
  hRatioDown->SetPoint(b,xbin,1);

  hRatioUp->SetPointError(b,0,0,0,errorTotalUp);
  hRatioDown->SetPointError(b,0,0,errorTotalDown,0);
  // hRatio->SetBinError(b, htotDATAbinErrorStat);


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

// !! --------------------------- PAD 2 --------------------------------------  !! //
  pad2->cd();
   xtitle = xtitle2;
 nbin   =   nbin2; 
 xmin   =   xmin2;
 xmax   =   xmax2;

 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);//ok
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleB);//ok
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleB);//ok
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleB);//ok
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleB);//ok
  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 
g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleB);//ok
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleB);//ok
 g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleB);//ok
 g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleB);//ok
  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);//ok
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleB);//ok
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleB);//ok
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleB);//ok 
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleB);//ok
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
  

// *****************************************************************************

 htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 htotData  = new TH1F("htotData","",nbin,xmin,xmax);
hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

  htotData->Sumw2();
  htotMC->Sumw2();
  hDataMC->Sumw2();
 
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


 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleB);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//------x);

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 
 f2_DY->cd();
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleB);
//  g2_DY->Sumw2();
//  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g2_DY, h_DY, 1,1);


 f1_VV->cd();
 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleB);
//  g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleB);
//  g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleB);
//  g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleB);
//  g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleB);
//  g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleB);
//  g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleB);
//  g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleB);
//  g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

  f3_ST->cd();
 g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleB);
 h_ST->Add(g3_ST, h_ST, 1, 1);

  f4_ST->cd();
 g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleB);
 h_ST->Add(g4_ST, h_ST, 1, 1);

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);
//  g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1*rescale,0);

f2_TT->cd();
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleB);
 h_TT->Add(g2_TT, h_TT, 1,1);
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
 htotMC->SetFillColorAlpha(ColorBlue, 1);
 htotMC->SetLineColor(ColorBlue);
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
 htotMC->GetYaxis()->SetTitle(ytitle+ HeaderB);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
 htotMC->SetMinimum(1); 
 htotMC->SetMaximum(htotMC->GetMaximum()*2); 
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);


 h_VV->Draw("HEsame"); 
 h_VV->SetFillColorAlpha(ColorOrange, 1);
 h_VV->SetLineColor(ColorOrange);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

 h_TTV->Draw("HEsame"); 
 h_TTV->SetFillColorAlpha(ColorNeutral, 1);
 h_TTV->SetLineColor(ColorNeutral);
 h_TTV->SetLineStyle(1);
 h_TTV->SetLineWidth(3);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColorAlpha(ColorDarkPurple, 1);
 h_ST->SetLineColor(ColorDarkPurple);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 h_TT->Draw("HEsame"); 
 h_TT->SetFillColorAlpha(ColorRed, 1);
 h_TT->SetLineColor(ColorRed);
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


  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  leg = new TLegend(0.64,0.60,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " #mu#mu data","PE1");
    leg->AddEntry(h_TT, " t#bar{t}","F");
    leg->AddEntry(htotMC, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");


  leg->Draw();


rap2->cd();

hRatio = new TH1F("hRatio","",nbin,xmin,xmax);
hRatio->Sumw2();
hRatio->Divide(htotData, htotMC, 1, 1);
hRatio->Draw("E"); 
hRatio->SetLineColor(1);
hRatio->SetLineStyle(1);
hRatio->SetLineWidth(1);
hRatio->SetMarkerColor(kBlack);
hRatio->SetMarkerStyle(20);
hRatio->SetMarkerSize(0.6);
hRatio->SetTickLength(0.10, "X"); hRatio->SetTickLength(0.05, "YZ");
hRatio->SetLabelOffset(0.02,"X");
hRatio->SetLabelOffset(0.02,"Y");
hRatio->SetLabelSize(0.12, "XY");
hRatio->SetLabelFont(42, "XYZ"); 
hRatio->SetTitleFont(42, "XYZ");
hRatio->SetTitleSize(0.14, "XYZ"); 
hRatio->SetTitleOffset(0.9,"X");
hRatio->SetTitleOffset(0.5,"Y");
hRatio->GetXaxis()->SetTitle(xtitle);
hRatio->GetXaxis()->SetTitleColor(1);
hRatio->GetXaxis()->SetNdivisions(509);
hRatio->GetYaxis()->SetTitle("Data / Sim.");
hRatio->GetYaxis()->SetTitleColor(1);
hRatio->GetYaxis()->SetNdivisions(509);
hRatio->SetNdivisions(509,"XYZ");
hRatio->SetMinimum(0.5); 
hRatio->SetMaximum(1.5);



const TH1F* hRatio_constB(hRatio);

hRatioUp = new TGraphAsymmErrors(hRatio_constB);
hRatioDown = new TGraphAsymmErrors(hRatio_constB);

for (unsigned int b = 1 ; b < hRatio->GetNbinsX(); b++)
{
  float htotMCbin = htotMC->GetBinContent(b);
  float htotMCbinErrorStat = htotMC->GetBinError(b);

  float bin =  hRatio->GetBinContent(b);
  float binStatError =  hRatio->GetBinError(b);

  float binSysErrorUp = hSumQuadratic_Up->GetBinContent(b);
  float binSysErrorDown = hSumQuadratic_Down->GetBinContent(b);
  // std::cout<<"bin:     "<<bin<<" binStatError: "<<binStatError<<" hSumQuadratic_Up->GetBinContent(b): "<<hSumQuadratic_Up->GetBinContent(b)<<" hSumQuadratic_Down->GetBinContent(b): "<<hSumQuadratic_Down->GetBinContent(b)<<std::endl;
  // std::cout<<"bin "<<bin<<" htotMC->GetBinError(b) : "<<htotMC->GetBinError(b)<<" and hRatio->GetBinError(b): "<<hRatio->GetBinError(b)<<std::endl;
  float errorTotalUp = sqrt(binSysErrorUp*binSysErrorUp + binStatError*binStatError);// unccorrelated SYS uncertainties
  float errorTotalDown = sqrt(binSysErrorDown*binSysErrorDown + binStatError*binStatError);// unccorrelated SYS uncertainties
  // std::cout<<"bin: "<<bin<<" binStatError: "<<binStatError<<" binSysErrorUp: "<<binSysErrorUp<<" binSysErrorDown: "<<binSysErrorDown<<" errorTotalUp: "<<errorTotalUp<<" errorTotalDown: "<<errorTotalDown<<std::endl;
  float xbin = hRatio->GetBinCenter(b);
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
// // !! --------------------------- PAD 3 --------------------------------------  !! //
  pad3->cd();
   xtitle = xtitle3;
 nbin   =   nbin3; 
 xmin   =   xmin3;
 xmax   =   xmax3;

 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);//ok
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleC);//ok
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleC);//ok
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleC);//ok
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleC);//ok
  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 
g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleC);//ok
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleC);//ok
 g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleC);//ok
 g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleC);//ok
  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);//ok
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleC);//ok
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleC);//ok
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleC);//ok 
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleC);//ok
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
  

// *****************************************************************************

 htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 htotData  = new TH1F("htotData","",nbin,xmin,xmax);
hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

  htotData->Sumw2();
  htotMC->Sumw2();
  hDataMC->Sumw2();
 
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
 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleC);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//------x);

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 
 f2_DY->cd();
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleC);
//  g2_DY->Sumw2();
//  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g2_DY, h_DY, 1,1);


 f1_VV->cd();
 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleC);
//  g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleC);
//  g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleC);
//  g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleC);
//  g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleC);
//  g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleC);
//  g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleC);
//  g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleC);
//  g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

  f3_ST->cd();
 g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleC);
 h_ST->Add(g3_ST, h_ST, 1, 1);

  f4_ST->cd();
 g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleC);
 h_ST->Add(g4_ST, h_ST, 1, 1);

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);
//  g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1*rescale,0);

f2_TT->cd();
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleC);
 h_TT->Add(g2_TT, h_TT, 1,1);
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
 htotMC->SetFillColorAlpha(ColorBlue, 1);
 htotMC->SetLineColor(ColorBlue);
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
 htotMC->GetYaxis()->SetTitle(ytitle+ HeaderC);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
 htotMC->SetMinimum(1); 
 htotMC->SetMaximum(htotMC->GetMaximum()*2); 
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);


 h_VV->Draw("HEsame"); 
 h_VV->SetFillColorAlpha(ColorOrange, 1);
 h_VV->SetLineColor(ColorOrange);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

 h_TTV->Draw("HEsame"); 
 h_TTV->SetFillColorAlpha(ColorNeutral, 1);
 h_TTV->SetLineColor(ColorNeutral);
 h_TTV->SetLineStyle(1);
 h_TTV->SetLineWidth(3);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColorAlpha(ColorDarkPurple, 1);
 h_ST->SetLineColor(ColorDarkPurple);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 h_TT->Draw("HEsame"); 
 h_TT->SetFillColorAlpha(ColorRed, 1);
 h_TT->SetLineColor(ColorRed);
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


  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  leg = new TLegend(0.64,0.60,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " #mu#mu data","PE1");
    leg->AddEntry(h_TT, " t#bar{t}","F");
    leg->AddEntry(htotMC, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");


  leg->Draw();


rap3->cd();

hRatio = new TH1F("hRatio","",nbin,xmin,xmax);
hRatio->Sumw2();
hRatio->Divide(htotData, htotMC, 1, 1);
hRatio->Draw("E"); 
hRatio->SetLineColor(1);
hRatio->SetLineStyle(1);
hRatio->SetLineWidth(1);
hRatio->SetMarkerColor(kBlack);
hRatio->SetMarkerStyle(20);
hRatio->SetMarkerSize(0.6);
hRatio->SetTickLength(0.10, "X"); hRatio->SetTickLength(0.05, "YZ");
hRatio->SetLabelOffset(0.02,"X");
hRatio->SetLabelOffset(0.02,"Y");
hRatio->SetLabelSize(0.12, "XY");
hRatio->SetLabelFont(42, "XYZ"); 
hRatio->SetTitleFont(42, "XYZ");
hRatio->SetTitleSize(0.14, "XYZ"); 
hRatio->SetTitleOffset(0.9,"X");
hRatio->SetTitleOffset(0.5,"Y");
hRatio->GetXaxis()->SetTitle(xtitle);
hRatio->GetXaxis()->SetTitleColor(1);
hRatio->GetXaxis()->SetNdivisions(509);
hRatio->GetYaxis()->SetTitle("Data / Sim.");
hRatio->GetYaxis()->SetTitleColor(1);
hRatio->GetYaxis()->SetNdivisions(509);
hRatio->SetNdivisions(509,"XYZ");
hRatio->SetMinimum(0.5); 
hRatio->SetMaximum(1.5);


const TH1F* hRatio_constC(hRatio);

hRatioUp = new TGraphAsymmErrors(hRatio_constC);
hRatioDown = new TGraphAsymmErrors(hRatio_constC);

for (unsigned int b = 1 ; b < hRatio->GetNbinsX(); b++)
{
  float htotMCbin = htotMC->GetBinContent(b);
  float htotMCbinErrorStat = htotMC->GetBinError(b);

  float bin =  hRatio->GetBinContent(b);
  float binStatError =  hRatio->GetBinError(b);

  float binSysErrorUp = hSumQuadratic_Up->GetBinContent(b);
  float binSysErrorDown = hSumQuadratic_Down->GetBinContent(b);
  // std::cout<<"bin:     "<<bin<<" binStatError: "<<binStatError<<" hSumQuadratic_Up->GetBinContent(b): "<<hSumQuadratic_Up->GetBinContent(b)<<" hSumQuadratic_Down->GetBinContent(b): "<<hSumQuadratic_Down->GetBinContent(b)<<std::endl;
  // std::cout<<"bin "<<bin<<" htotMC->GetBinError(b) : "<<htotMC->GetBinError(b)<<" and hRatio->GetBinError(b): "<<hRatio->GetBinError(b)<<std::endl;
  float errorTotalUp = sqrt(binSysErrorUp*binSysErrorUp + binStatError*binStatError);// unccorrelated SYS uncertainties
  float errorTotalDown = sqrt(binSysErrorDown*binSysErrorDown + binStatError*binStatError);// unccorrelated SYS uncertainties
  // std::cout<<"bin: "<<bin<<" binStatError: "<<binStatError<<" binSysErrorUp: "<<binSysErrorUp<<" binSysErrorDown: "<<binSysErrorDown<<" errorTotalUp: "<<errorTotalUp<<" errorTotalDown: "<<errorTotalDown<<std::endl;
  float xbin = hRatio->GetBinCenter(b);
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
// // *****************************************************************************
// // !! --------------------------- PAD4 --------------------------------------  !! //
  pad4->cd();
   xtitle = xtitle4;
 nbin   =   nbin4; 
 xmin   =   xmin4;
 xmax   =   xmax4;

 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleD);//ok
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleD);//ok
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleD);//ok
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleD);//ok
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleD);//ok
  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 
g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleD);//ok
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleD);//ok
 g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleD);//ok
 g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleD);//ok

  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleD);//ok
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleD);//ok
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleD);//ok
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleD);//ok 
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleD);//ok
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
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

// *****************************************************************************

 htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 htotData  = new TH1F("htotData","",nbin,xmin,xmax);
hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

  htotData->Sumw2();
  htotMC->Sumw2();
  hDataMC->Sumw2();
 

 f1_Data_emu->cd();
 g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleD);
 htotData->Add(g1_Data_emu, htotData, 1, 0);

//------x);

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleD);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 
 f2_DY->cd();
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleD);
//  g2_DY->Sumw2();
//  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g2_DY, h_DY, 1,1);


 f1_VV->cd();
 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleD);
//  g1_VV->Sumw2();
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleD);
//  g2_VV->Sumw2();
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleD);
//  g3_VV->Sumw2();
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleD);
//  g1_TTV->Sumw2();
 h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleD);
//  g2_TTV->Sumw2();
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleD);
//  g3_TTV->Sumw2();
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleD);
//  g1_ST->Sumw2();
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleD);
//  g2_ST->Sumw2();
 h_ST->Add(g2_ST, h_ST, 1, 1);

  f3_ST->cd();
 g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleD);
 h_ST->Add(g3_ST, h_ST, 1, 1);

  f4_ST->cd();
 g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleD);
 h_ST->Add(g4_ST, h_ST, 1, 1);

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleD);
//  g1_TT->Sumw2();
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1*rescale,0);

f2_TT->cd();
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleD);
 h_TT->Add(g2_TT, h_TT, 1,1);
// std::cout<<"TT: "<<1<<"with norm "<<norm<<std::endl;
// std::cout<<"tt itnegral "<<h_TT->Integral(0,80)<<std::endl;
 rescale = 1.0;
 htotMC->Add(htotMC, h_DY, 1, rescale);
 htotMC->Add(htotMC, h_VV, 1, rescale);
 htotMC->Add(htotMC, h_TTV, 1, rescale);
 htotMC->Add(htotMC, h_ST, 1, rescale);
 htotMC->Add(htotMC, h_TT, 1, rescale);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);

 htotMC->Draw("HE"); 
  htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(ColorBlue, 1);
 htotMC->SetLineColor(ColorBlue);
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
 htotMC->GetYaxis()->SetTitle(+ HeaderD);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
 htotMC->SetMinimum(1); 
 htotMC->SetMaximum(htotMC->GetMaximum()*2); 
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);


 h_VV->Draw("HEsame"); 
 h_VV->SetFillColorAlpha(ColorOrange, 1);
 h_VV->SetLineColor(ColorOrange);
 h_VV->SetLineStyle(1);
 h_VV->SetLineWidth(3);

 h_TTV->Draw("HEsame"); 
 h_TTV->SetFillColorAlpha(ColorNeutral, 1);
 h_TTV->SetLineColor(ColorNeutral);
 h_TTV->SetLineStyle(1);
 h_TTV->SetLineWidth(3);

 h_ST->Draw("HEsame"); 
 h_ST->SetFillColorAlpha(ColorDarkPurple, 1);
 h_ST->SetLineColor(ColorDarkPurple);
 h_ST->SetLineStyle(1);
 h_ST->SetLineWidth(3);

 h_TT->Draw("HEsame"); 
 h_TT->SetFillColorAlpha(ColorRed, 1);
 h_TT->SetLineColor(ColorRed);
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


  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  leg = new TLegend(0.64,0.60,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotData, " #mu#mu data","PE1");
  leg->AddEntry(h_TT, " t#bar{t}","F");
  leg->AddEntry(htotMC, " DY","F");
  leg->AddEntry(h_VV, " WW, WZ, ZZ","F");
  leg->AddEntry(h_TTV," t#bar{t}W, t#bar{t}Z, t#bar{t}WW","F");
  leg->AddEntry(h_ST, " tW","F");
  

  leg->Draw();


rap4->cd();

hRatio = new TH1F("hRatio","",nbin,xmin,xmax);
hRatio->Sumw2();
hRatio->Divide(htotData, htotMC, 1, 1);
hRatio->Draw("E"); 
hRatio->SetLineColor(1);
hRatio->SetLineStyle(1);
hRatio->SetLineWidth(1);
hRatio->SetMarkerColor(kBlack);
hRatio->SetMarkerStyle(20);
hRatio->SetMarkerSize(0.6);
hRatio->SetTickLength(0.10, "X"); hRatio->SetTickLength(0.05, "YZ");
hRatio->SetLabelOffset(0.02,"X");
hRatio->SetLabelOffset(0.02,"Y");
hRatio->SetLabelSize(0.12, "XY");
hRatio->SetLabelFont(42, "XYZ"); 
hRatio->SetTitleFont(42, "XYZ");
hRatio->SetTitleSize(0.14, "XYZ"); 
hRatio->SetTitleOffset(0.9,"X");
hRatio->SetTitleOffset(0.5,"Y");
hRatio->GetXaxis()->SetTitle(xtitle);
hRatio->GetXaxis()->SetTitleColor(1);
hRatio->GetXaxis()->SetNdivisions(509);
hRatio->GetYaxis()->SetTitle("Data / Sim.");
hRatio->GetYaxis()->SetTitleColor(1);
hRatio->GetYaxis()->SetNdivisions(509);
hRatio->SetNdivisions(509,"XYZ");
hRatio->SetMinimum(0.5); 
hRatio->SetMaximum(1.5);


const TH1F* hRatio_constD(hRatio);

hRatioUp = new TGraphAsymmErrors(hRatio_constD);
hRatioDown = new TGraphAsymmErrors(hRatio_constD);

for (unsigned int b = 1 ; b < hRatio->GetNbinsX(); b++)
{
  float htotMCbin = htotMC->GetBinContent(b);
  float htotMCbinErrorStat = htotMC->GetBinError(b);

  float bin =  hRatio->GetBinContent(b);
  float binStatError =  hRatio->GetBinError(b);

  float binSysErrorUp = hSumQuadratic_Up->GetBinContent(b);
  float binSysErrorDown = hSumQuadratic_Down->GetBinContent(b);
  // std::cout<<"bin:     "<<bin<<" binStatError: "<<binStatError<<" hSumQuadratic_Up->GetBinContent(b): "<<hSumQuadratic_Up->GetBinContent(b)<<" hSumQuadratic_Down->GetBinContent(b): "<<hSumQuadratic_Down->GetBinContent(b)<<std::endl;
  // std::cout<<"bin "<<bin<<" htotMC->GetBinError(b) : "<<htotMC->GetBinError(b)<<" and hRatio->GetBinError(b): "<<hRatio->GetBinError(b)<<std::endl;
  float errorTotalUp = sqrt(binSysErrorUp*binSysErrorUp + binStatError*binStatError);// unccorrelated SYS uncertainties
  float errorTotalDown = sqrt(binSysErrorDown*binSysErrorDown + binStatError*binStatError);// unccorrelated SYS uncertainties
  // std::cout<<"bin: "<<bin<<" binStatError: "<<binStatError<<" binSysErrorUp: "<<binSysErrorUp<<" binSysErrorDown: "<<binSysErrorDown<<" errorTotalUp: "<<errorTotalUp<<" errorTotalDown: "<<errorTotalDown<<std::endl;
  float xbin = hRatio->GetBinCenter(b);
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
// // hRatio->SaveAs("./"+htitleD+".root");
// // !! 



  TString name = htitleA+"_TotalErr";
  if (WLepton)
    {
      name = name+"_WLepton";
    }
  if (WGen)
    {
      name = name+"_WGen";
    }
  c1->SaveAs("./"+name+".pdf");

  f1_Data_emu->Close();
 f1_DY->Close();
 f2_DY->Close();
 f1_TT->Close(); 
 f2_TT->Close();  
 f1_ST->Close();  
 f2_ST->Close();
 f3_ST->Close();
 f4_ST->Close();  
 f1_TTV->Close(); 
 f2_TTV->Close(); 
 f3_TTV->Close(); 
 f1_VV->Close();  
 f2_VV->Close(); 
 f3_VV->Close();  
}