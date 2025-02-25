#include <iostream>
#include <TROOT.h>
#include "TH1.h"
#include "TColor.h"

void plot(int method, TString Year, TString SYST , TString Plots)
{
int stati=0;
bool fit= 1;
bool logy=1;

float hmin = 0.5;
float hmax = 1E6;
float hmaxBD = 1E6;

 TString Yearcor = Year;
 if (Year == "2016PRE") Yearcor = "2016preVFP";
 if (Year == "2016POST") Yearcor = "2016";

TString extension = SYST;
TFile * theoutputfile = new TFile( "./SYST/"+SYST+"_SYST.root" , "UPDATE");

 TFile* f1_DY  = new TFile("../../MC_EMU_03_02_2025/DATAMC_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f2_DY  = new TFile("../../MC_EMU_03_02_2025/DATAMC_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f1_TT  = new TFile("../../MC_EMU_03_02_2025/DATAMC_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_TT  = new TFile("../../MC_EMU_03_02_2025/DATAMC_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_ST  = new TFile("../../MC_EMU_03_02_2025/DATAMC_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_ST  = new TFile("../../MC_EMU_03_02_2025/DATAMC_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_TTV = new TFile("../../MC_EMU_03_02_2025/DATAMC_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8.root");
 TFile* f2_TTV = new TFile("../../MC_EMU_03_02_2025/DATAMC_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f3_TTV = new TFile("../../MC_EMU_03_02_2025/DATAMC_TTWW_TuneCP5_13TeV-madgraph-pythia8.root");
 TFile* f1_VV  = new TFile("../../MC_EMU_03_02_2025/DATAMC_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_VV  = new TFile("../../MC_EMU_03_02_2025/DATAMC_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 TFile* f3_VV  = new TFile("../../MC_EMU_03_02_2025/DATAMC_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 
 TString MCFILE[12] = {
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
                  "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_"
};

// extension is <SYST>Up or <SYST>Down
TString EXTRA = "";
  if (SYST == "JECUp") EXTRA = "_JECUp";
  else if (SYST == "JECDown") EXTRA = "_JECDown";
  else if (SYST == "JERUp" ) EXTRA = "_JERUp";
  else if (SYST == "JERDown" ) EXTRA = "_JERDown";


  TFile* f1_DY_SYST  = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_"+extension+".root");
 TFile* f2_DY_SYST  = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_"+extension+".root");
 TFile* f1_TT_SYST  = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f2_TT_SYST  = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f1_ST_SYST  = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f2_ST_SYST  = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f1_TTV_SYST = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8_"+extension+".root");
 TFile* f2_TTV_SYST = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8_"+extension+".root");
 TFile* f3_TTV_SYST = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_TTWW_TuneCP5_13TeV-madgraph-pythia8_"+extension+".root");
 TFile* f1_VV_SYST  = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f2_VV_SYST  = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_"+extension+".root");
 TFile* f3_VV_SYST  = new TFile("../../MC_EMU_03_02_2025"+EXTRA+"/DATAMC_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_"+extension+".root");


TString MCFILE_SYST[12] = {
                  "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_"+extension+"_",
                  "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_"+extension+"_",
                  "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_"+extension+"_",
                  "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_"+extension+"_",
                  "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+"_",
                  "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+"_",
                  "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8_"+extension+"_",
                  "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8_"+extension+"_",
                  "TTWW_TuneCP5_13TeV-madgraph-pythia8_"+extension+"_",
                  "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_"+extension+"_",
                  "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_"+extension+"_",
                  "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_"+extension+"_"
};


 TString ytitle = "Events";
 TString xtitle = "vertex BDT score"; 
 int nbin = 26; 
 float xmin = -1.04;
 float xmax =  1.04;
//$$
 TString HeaderCMS = "CMS";

 if (Year == "2016") HeaderCMS = "2016                            36.3 fb^{-1} (13 TeV)";
 if (Year == "2017") HeaderCMS = "2017                            41.5 fb^{-1} (13 TeV)";
 if (Year == "2018") HeaderCMS = "2018                            59.8 fb^{-1} (13 TeV)";


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
    if (Method == 0)
      {
        htitleA = "Tree_filter_";
        nbin1 = 2; 
        xmin1 = 0.;
        xmax1 =  2;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 Vtx"; 
        xtitle1 = "Filter";

        htitleB = "Vertices_NoSel";
        nbin2 = 100; 
        xmin2 = 0;
        xmax2 =  100;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 Vtx All"; 
        xtitle2 = "nVtx";

        htitleC = "Vertices_filtercut_";
        nbin3 = 100; 
        xmin3 = 0;
        xmax3 = 100;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 VtxAll";
        xtitle3 = "nVtx";


        htitleD = "Tree_Mumu_nosel";
        nbin4 = 150; 
        xmin4 = 0;
        xmax4 = 600;
        // HeaderA = "TT + 20<pt<80";
        // HeaderNVtx = "2 VtxAll";
        xtitle4 = "M_{e#mu}";
      }


    if (Method == 1)
      {
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
      }

    if (Method == 2)
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
      }

    if (Method == 3)
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
      }


    if (Method == 4)
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

        htitleC = "hData_jet_btag_Deepjet_";
        nbin3 = 50;
        xmin3 = 0;
        xmax3 = 1;
        xtitle3 = "jet_{btag}";

        htitleD = "hData_jet_HadronFlavour_";
        nbin4 = 7;
        xmin4 = -0.5;
        xmax4 = 6.5;
        xtitle4 = "jet_{HadronFlavour}";
      }

    if (Method == 5)
      {
        htitleA = "njet_NoSel";
        nbin1 = 20;
        xmin1 = 0;
        xmax1 = 20;
        xtitle1 = "n_{jet}";

        htitleB = "njetNOmu_NoSel";
        nbin2 = 20;
        xmin2 = 0;
        xmax2 = 20;
        xtitle2 = "n_{jet NoLepton}";

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
      }

    // !!----------------------
    if (Method == 6)
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
      }


  // !!----------------------
    if (Method == 7)
      {
        htitleA = "Nmu_";
        nbin1 = 10;
        xmin1 = 0;
        xmax1 = 10;
        xtitle1 = "n_{#mu}";

        htitleB = "Muon_pt_";
        nbin2 = 100;
        xmin2 = 0;
        xmax2 = 500;
        xtitle2 = "#mu_{pt}";

        htitleC = "Muon_PFIsoLoose_";
        nbin3 = 2;
        xmin3 = 0;
        xmax3 = 2;
        xtitle3 = "PFIsoLoose";

        htitleD = "Muon_MiniIsoTight_";
        nbin4 = 2;
        xmin4 = 0;
        xmax4 = 2;
        xtitle4 = "MiniIsoTight";
      }

  // !!----------------------
    if (Method == 8)
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
      }


    // !!----------------------
    if (Method == 9)
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
        
      }
    // !!----------------------
    if (Method == 10)
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
        xtitle3 = "n_{TRK}";

        htitleD = "nLostTracks_";
        nbin4 = 25;
        xmin4 = 0;
        xmax4 = 25;
        xtitle4 = "n_{LostTRK}";
      }

    // !!----------------------
    if (Method == 11)
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

        htitleC = "Hemi_nJet_";
        nbin3 = 10;
        xmin3 = 0;
        xmax3 = 10;
        xtitle3 = "Hemi n{jet}";

        htitleD = "Hemi_nJetNoMu_";
        nbin4 = 10;
        xmin4 = 0;
        xmax4 = 10;
        xtitle4 = "Hemi n_{jetNoMu}";
      }

    // !!----------------------
    if (Method == 12)
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
        nbin4 = 100;
        xmin4 = 0;
        xmax4 = 500;
        xtitle4 = "Hemi_{Mass}";
      }
          // !!----------------------
    if (Method == 13)
      {
        htitleA = "K0_mass_";
        nbin1 = 202;
        xmin1 = 0.42;
        xmax1 = 0.58;
        xtitle1 = "K0_{mass}";

        htitleB = "K0_pt_";
        nbin2 = 200;
        xmin2 = 0;
        xmax2 = 200;
        xtitle2 = "K0_{pt}";

        htitleC = "Reco_K0_mass_";
        nbin3 = 202;
        xmin3 = 0.42;
        xmax3 = 0.58;
        xtitle3 = "RecoK0_{mass}";

        htitleD = "Reco_K0_pt_";
        nbin4 = 200;
        xmin4 = 0;
        xmax4 = 200;
        xtitle4 = "RecoK0_{pt}";
      }
  // !!----------------------
    if (Method == 14)
      {
        htitleA = "L0_mass_";
        nbin1 = 202;
        xmin1 = 1.06;
        xmax1 = 1.18;
        xtitle1 = "L0_{mass}";

        htitleB = "L0_pt_";
        nbin2 = 200;
        xmin2 = 0;
        xmax2 = 200;
        xtitle2 = "L0_{pt}";

        htitleC = "Reco_L0_mass_";
        nbin3 = 202;
        xmin3 = 1.06;
        xmax3 = 1.18;
        xtitle3 = "RecoL0_{mass}";

        htitleD = "Reco_L0_pt_";
        nbin4 = 200;
        xmin4 = 0;
        xmax4 = 200;
        xtitle4 = "RecoL0_{pt}";
      }

  // !!----------------------
    if (Method == 15)
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
        nbin3 = 2000;
        xmin3 = 0;
        xmax3 = 4000;
        xtitle3 = "SecInt_{pt}";

        htitleD = "SecInt_dzSig_Selec";
        nbin4 = 200;
        xmin4 = 0;
        xmax4 = 2000;
        xtitle4 = "SecInt_dzSig";
      }
 
        // !!----------------------
    if (Method == 16)
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
        nbin3 = 2000;
        xmin3 = 0;
        xmax3 = 4000;
        xtitle3 = "SecInt_{pt}";

        htitleD = "SecInt_dzSig_TrackerMatched";
        nbin4 = 200;
        xmin4 = 0;
        xmax4 = 2000;
        xtitle4 = "SecInt IP_{dz}";
      }                                                      
                                                            
        // !!----------------------
    if (Method == 17)
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
      }      

        // !!----------------------
    if (Method == 18)
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
      }  

        // !!----------------------
    if (Method == 20)
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
      }  
    //--------------------------------------------

    if (Method == 21)
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
      }


    if (Method == 22)
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
    }

    if  (Method == 23)
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
    }


    if  (Method == 24)
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
    }

    if  (Method == 25)
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
    }


    if  (Method == 26)
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
    }

        // !!----------------------
    if (Method == 28)
      {
        htitleA = "SecInt_r_Selec";
        nbin1 = 400;
        xmin1 = 0;
        xmax1 = 200;
        xtitle1 = "SecInt r_{cm}";

        htitleB = "SecInt_z_Selec";
        nbin2 = 800;
        xmin2 = -200;
        xmax2 = 200;
        xtitle2 = "SecInt z_{cm}";

        htitleC = "SecInt_pt_Selec";
        nbin3 = 2000;
        xmin3 = 0;
        xmax3 = 4000;
        xtitle3 = "SecInt p_{t} [GeV]";

        htitleD = "SecInt_dzSig_Selec";
        nbin4 = 200;
        xmin4 = 0;
        xmax4 = 2000;
        xtitle4 = "SecInt dzSig";
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
 TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 TH1F* g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);//ok
 TH1F* g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);//ok
 TH1F*  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 TH1F* g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);//ok
 TH1F* g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);//ok 
 TH1F* g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);//ok
 TH1F*  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
  
  

 TH1F* g1_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleA);//ok
 TH1F* g2_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[1]+htitleA);//ok
 TH1F*  h_DY_SYST = new TH1F("h_DY_SYST","",nbin,xmin,xmax);

 
 TH1F* g1_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[9]+htitleA);//ok
 TH1F* g2_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[10]+htitleA);//ok
 TH1F* g3_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[11]+htitleA);//ok
 TH1F*  h_VV_SYST = new TH1F("h_VV_SYST","",nbin,xmin,xmax);
 
TH1F* g1_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[4]+htitleA);//ok
 TH1F* g2_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[5]+htitleA);//ok
 TH1F*  h_ST_SYST = new TH1F("h_ST_SYST","",nbin,xmin,xmax);

 TH1F* g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleA);//ok
 TH1F* g2_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[3]+htitleA);//ok
 TH1F*  h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);

 TH1F* g1_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[6]+htitleA);//ok
 TH1F* g2_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[7]+htitleA);//ok 
 TH1F* g3_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[8]+htitleA);//ok
 TH1F*  h_TTV_SYST = new TH1F("h_TTV_SYST","",nbin,xmin,xmax);
// *****************************************************************************

 pad1->cd();

 TH1F* htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 TH1F* htotMC_SYST  = new TH1F("htotMC_SYST","",nbin,xmin,xmax);
 TH1F* hRatioSYST_NOM = new TH1F("hRatioSYST_NOM","",nbin,xmin,xmax);

  htotMC_SYST->Sumw2();
  htotMC->Sumw2();
  hRatioSYST_NOM->Sumw2();
 

//------x);

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 
 f2_DY->cd();
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
 h_DY->Add(g2_DY, h_DY, 1,1);


 f1_VV->cd();
 
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);
 h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);
 h_VV->Add(g2_VV, h_VV, 1, 1);

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

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 0.5,0);

f2_TT->cd();
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);
 h_TT->Add(g2_TT, h_TT, 1,1);


 htotMC->Add(htotMC, h_DY, 1, 1);
 htotMC->Add(htotMC, h_VV, 1, 1);
 htotMC->Add(htotMC, h_TTV, 1, 1);
 htotMC->Add(htotMC, h_ST, 1, 1);
 htotMC->Add(htotMC, h_TT, 1, 1);



 htotMC->Draw("PE1"); 
  htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(ColorNeutral, 1);
 htotMC->SetLineColor(ColorNeutral);
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
 htotMC->SetMinimum(htotMC->GetMinimum()/2+1); 
 htotMC->SetMaximum(htotMC->GetMaximum()*2); 

 f1_DY_SYST->cd();
 g1_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleA);
 g1_DY_SYST->Sumw2();
 h_DY_SYST = new TH1F("h_DY_SYST","",nbin,xmin,xmax);
 h_DY_SYST->Add(g1_DY_SYST, h_DY_SYST, 1,0);

 f2_DY_SYST->cd();
 g2_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[1]+htitleA);
 h_DY_SYST->Add(g2_DY_SYST, h_DY_SYST, 1,1);

 f1_VV_SYST->cd();
 g1_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[9]+htitleA);
 h_VV_SYST = new TH1F("h_VV_SYST","",nbin,xmin,xmax);
 h_VV_SYST->Add(g1_VV_SYST, h_VV_SYST, 1,0);

 f2_VV_SYST->cd();
 g2_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[10]+htitleA);
 h_VV_SYST->Add(g2_VV_SYST, h_VV_SYST, 1, 1);

 f3_VV_SYST->cd();
 g3_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[11]+htitleA);
 h_VV_SYST->Add(g3_VV_SYST, h_VV_SYST, 1, 1);

 f1_TTV_SYST->cd();
 g1_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[6]+htitleA);
 h_TTV_SYST = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV_SYST->Add(g1_TTV_SYST, h_TTV_SYST, 1,0);

 f2_TTV_SYST->cd();
 
 g2_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[7]+htitleA);
 h_TTV_SYST->Add(g2_TTV_SYST, h_TTV_SYST, 1, 1);

 f3_TTV_SYST->cd();
 g3_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[8]+htitleA);
 h_TTV_SYST->Add(g3_TTV_SYST, h_TTV_SYST, 1, 1);

 h_ST_SYST = new TH1F("h_ST_SYST","",nbin,xmin,xmax);
 f1_ST_SYST->cd();
 g1_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[4]+htitleA);
 h_ST_SYST->Add(g1_ST_SYST, h_ST_SYST, 1,0);

 f2_ST_SYST->cd();
 g2_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[5]+htitleA);
 h_ST_SYST->Add(g2_ST_SYST, h_ST_SYST, 1, 1);

 f1_TT_SYST->cd();
g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleA);
 h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);
 h_TT_SYST->Add(g1_TT_SYST, h_TT_SYST, 0.5,0);

f2_TT_SYST->cd();
 g2_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[3]+htitleA);
 h_TT_SYST->Add(g2_TT_SYST, h_TT_SYST, 1,1);


 htotMC_SYST->Add(htotMC_SYST, h_DY_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_VV_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_TTV_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_ST_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_TT_SYST, 1, 1);


htotMC_SYST->Draw("PE1same");
htotMC_SYST->SetFillStyle(1001);
 htotMC_SYST->SetFillColorAlpha(ColorRed, 1);
 htotMC_SYST->SetLineColor(ColorRed);

// htotMC_SYST->SetMarkerStyle(20);
// htotMC_SYST->SetMarkerSize(1);
// htotMC_SYST->SetMarkerColor(kBlack);
// htotMC_SYST->SetLineColor(kBlack);
// htotMC_SYST->SetLineWidth(1);


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
    leg->AddEntry(htotMC, " DY","F");
    leg->AddEntry(htotMC_SYST, " DY SYST","F");
  leg->Draw();


rap1->cd();

TH1F* hRatio = new TH1F(htitleA+"_"+SYST,"",nbin,xmin,xmax);
hRatio->Sumw2();
hRatio->Divide(htotMC_SYST, htotMC, 1, 1);
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
hRatio->GetYaxis()->SetTitle("SYST / Sim.");
hRatio->GetYaxis()->SetTitleColor(1);
hRatio->GetYaxis()->SetNdivisions(509);
hRatio->SetNdivisions(509,"XYZ");
hRatio->SetMinimum(0.5); 
hRatio->SetMaximum(1.5); 


        theoutputfile->cd();
        hRatio->Write();

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
  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);//ok
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleB);//ok
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleB);//ok
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleB);//ok 
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleB);//ok
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
  

// // *****************************************************************************
 g1_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleB);//ok
 g2_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[1]+htitleB);//ok
  h_DY_SYST = new TH1F("h_DY_SYST","",nbin,xmin,xmax);

 
 g1_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[9]+htitleB);//ok
 g2_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[10]+htitleB);//ok
 g3_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[11]+htitleB);//ok
  h_VV_SYST = new TH1F("h_VV_SYST","",nbin,xmin,xmax);
 
g1_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[4]+htitleB);//ok
 g2_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[5]+htitleB);//ok
  h_ST_SYST = new TH1F("h_ST_SYST","",nbin,xmin,xmax);

 g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleB);//ok
 g2_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[3]+htitleB);//ok
  h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);

 g1_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[6]+htitleB);//ok
 g2_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[7]+htitleB);//ok 
 g3_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[8]+htitleB);//ok
  h_TTV_SYST = new TH1F("h_TTV_SYST","",nbin,xmin,xmax);
// *****************************************************************************


 htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 htotMC_SYST  = new TH1F("htotMC_SYST","",nbin,xmin,xmax);
 hRatioSYST_NOM = new TH1F("hRatioSYST_NOM","",nbin,xmin,xmax);

  htotMC_SYST->Sumw2();
  htotMC->Sumw2();
  hRatioSYST_NOM->Sumw2();
 

//------x);

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 
 f2_DY->cd();
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleB);
 h_DY->Add(g2_DY, h_DY, 1,1);


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

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 0.5,0);

f2_TT->cd();
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleB);
 h_TT->Add(g2_TT, h_TT, 1,1);


 htotMC->Add(htotMC, h_DY, 1, 1);
 htotMC->Add(htotMC, h_VV, 1, 1);
 htotMC->Add(htotMC, h_TTV, 1, 1);
 htotMC->Add(htotMC, h_ST, 1, 1);
 htotMC->Add(htotMC, h_TT, 1, 1);



 htotMC->Draw("PE1"); 
  htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(ColorNeutral, 1);
 htotMC->SetLineColor(ColorNeutral);
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
 htotMC->SetMinimum(htotMC->GetMinimum()/2+1); 
 htotMC->SetMaximum(htotMC->GetMaximum()*2); 

 f1_DY_SYST->cd();
 g1_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleB);
 g1_DY_SYST->Sumw2();
 h_DY_SYST = new TH1F("h_DY_SYST","",nbin,xmin,xmax);
 h_DY_SYST->Add(g1_DY_SYST, h_DY_SYST, 1,0);

 f2_DY_SYST->cd();
 g2_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[1]+htitleB);
 h_DY_SYST->Add(g2_DY_SYST, h_DY_SYST, 1,1);

 f1_VV_SYST->cd();
 g1_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[9]+htitleB);
 h_VV_SYST = new TH1F("h_VV_SYST","",nbin,xmin,xmax);
 h_VV_SYST->Add(g1_VV_SYST, h_VV_SYST, 1,0);

 f2_VV_SYST->cd();
 g2_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[10]+htitleB);
 h_VV_SYST->Add(g2_VV_SYST, h_VV_SYST, 1, 1);

 f3_VV_SYST->cd();
 g3_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[11]+htitleB);
 h_VV_SYST->Add(g3_VV_SYST, h_VV_SYST, 1, 1);

 f1_TTV_SYST->cd();
 g1_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[6]+htitleB);
 h_TTV_SYST = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV_SYST->Add(g1_TTV_SYST, h_TTV_SYST, 1,0);

 f2_TTV_SYST->cd();
 
 g2_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[7]+htitleB);
 h_TTV_SYST->Add(g2_TTV_SYST, h_TTV_SYST, 1, 1);

 f3_TTV_SYST->cd();
 g3_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[8]+htitleB);
 h_TTV_SYST->Add(g3_TTV_SYST, h_TTV_SYST, 1, 1);

 h_ST_SYST = new TH1F("h_ST_SYST","",nbin,xmin,xmax);
 f1_ST_SYST->cd();
 g1_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[4]+htitleB);
 h_ST_SYST->Add(g1_ST_SYST, h_ST_SYST, 1,0);

 f2_ST_SYST->cd();
 g2_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[5]+htitleB);
 h_ST_SYST->Add(g2_ST_SYST, h_ST_SYST, 1, 1);

 f1_TT_SYST->cd();
g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleB);
 h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);
 h_TT_SYST->Add(g1_TT_SYST, h_TT_SYST, 0.5,0);

f2_TT_SYST->cd();
 g2_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[3]+htitleB);
 h_TT_SYST->Add(g2_TT_SYST, h_TT_SYST, 1,1);


 htotMC_SYST->Add(htotMC_SYST, h_DY_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_VV_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_TTV_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_ST_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_TT_SYST, 1, 1);


htotMC_SYST->Draw("PE1same");
htotMC_SYST->SetFillStyle(1001);
 htotMC_SYST->SetFillColorAlpha(ColorRed, 1);
 htotMC_SYST->SetLineColor(ColorRed);

// htotMC_SYST->SetMarkerStyle(20);
// htotMC_SYST->SetMarkerSize(1);
// htotMC_SYST->SetMarkerColor(kBlack);
// htotMC_SYST->SetLineColor(kBlack);
// htotMC_SYST->SetLineWidth(1);


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
    leg->AddEntry(htotMC, " DY","F");
    leg->AddEntry(htotMC_SYST, " DY SYST","F");
  leg->Draw();


rap2->cd();

 hRatio = new TH1F(htitleB+"_"+SYST,"",nbin,xmin,xmax);
hRatio->Sumw2();
hRatio->Divide(htotMC_SYST, htotMC, 1, 1);
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
hRatio->GetYaxis()->SetTitle("SYST / Sim.");
hRatio->GetYaxis()->SetTitleColor(1);
hRatio->GetYaxis()->SetNdivisions(509);
hRatio->SetNdivisions(509,"XYZ");
hRatio->SetMinimum(0.5); 
hRatio->SetMaximum(1.5); 


        theoutputfile->cd();
        hRatio->Write();
// // *****************************************************************************
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
  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);//ok
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleC);//ok
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleC);//ok
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleC);//ok 
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleC);//ok
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
  

// // *****************************************************************************
 g1_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleC);//ok
 g2_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[1]+htitleC);//ok
  h_DY_SYST = new TH1F("h_DY_SYST","",nbin,xmin,xmax);

 
 g1_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[9]+htitleC);//ok
 g2_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[10]+htitleC);//ok
 g3_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[11]+htitleC);//ok
  h_VV_SYST = new TH1F("h_VV_SYST","",nbin,xmin,xmax);
 
g1_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[4]+htitleC);//ok
 g2_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[5]+htitleC);//ok
  h_ST_SYST = new TH1F("h_ST_SYST","",nbin,xmin,xmax);

 g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleC);//ok
 g2_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[3]+htitleC);//ok
  h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);

 g1_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[6]+htitleC);//ok
 g2_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[7]+htitleC);//ok 
 g3_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[8]+htitleC);//ok
  h_TTV_SYST = new TH1F("h_TTV_SYST","",nbin,xmin,xmax);
// *****************************************************************************


 htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 htotMC_SYST  = new TH1F("htotMC_SYST","",nbin,xmin,xmax);
 hRatioSYST_NOM = new TH1F("hRatioSYST_NOM","",nbin,xmin,xmax);

  htotMC_SYST->Sumw2();
  htotMC->Sumw2();
  hRatioSYST_NOM->Sumw2();
 

//------x);

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 
 f2_DY->cd();
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleC);
 h_DY->Add(g2_DY, h_DY, 1,1);


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

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 0.5,0);

f2_TT->cd();
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleC);
 h_TT->Add(g2_TT, h_TT, 1,1);


 htotMC->Add(htotMC, h_DY, 1, 1);
 htotMC->Add(htotMC, h_VV, 1, 1);
 htotMC->Add(htotMC, h_TTV, 1, 1);
 htotMC->Add(htotMC, h_ST, 1, 1);
 htotMC->Add(htotMC, h_TT, 1, 1);



 htotMC->Draw("PE1"); 
  htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(ColorNeutral, 1);
 htotMC->SetLineColor(ColorNeutral);
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
 htotMC->SetMinimum(htotMC->GetMinimum()/2+1); 
 htotMC->SetMaximum(htotMC->GetMaximum()*2); 

 f1_DY_SYST->cd();
 g1_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleC);
 g1_DY_SYST->Sumw2();
 h_DY_SYST = new TH1F("h_DY_SYST","",nbin,xmin,xmax);
 h_DY_SYST->Add(g1_DY_SYST, h_DY_SYST, 1,0);

 f2_DY_SYST->cd();
 g2_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[1]+htitleC);
 h_DY_SYST->Add(g2_DY_SYST, h_DY_SYST, 1,1);

 f1_VV_SYST->cd();
 g1_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[9]+htitleC);
 h_VV_SYST = new TH1F("h_VV_SYST","",nbin,xmin,xmax);
 h_VV_SYST->Add(g1_VV_SYST, h_VV_SYST, 1,0);

 f2_VV_SYST->cd();
 g2_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[10]+htitleC);
 h_VV_SYST->Add(g2_VV_SYST, h_VV_SYST, 1, 1);

 f3_VV_SYST->cd();
 g3_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[11]+htitleC);
 h_VV_SYST->Add(g3_VV_SYST, h_VV_SYST, 1, 1);

 f1_TTV_SYST->cd();
 g1_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[6]+htitleC);
 h_TTV_SYST = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV_SYST->Add(g1_TTV_SYST, h_TTV_SYST, 1,0);

 f2_TTV_SYST->cd();
 
 g2_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[7]+htitleC);
 h_TTV_SYST->Add(g2_TTV_SYST, h_TTV_SYST, 1, 1);

 f3_TTV_SYST->cd();
 g3_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[8]+htitleC);
 h_TTV_SYST->Add(g3_TTV_SYST, h_TTV_SYST, 1, 1);

 h_ST_SYST = new TH1F("h_ST_SYST","",nbin,xmin,xmax);
 f1_ST_SYST->cd();
 g1_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[4]+htitleC);
 h_ST_SYST->Add(g1_ST_SYST, h_ST_SYST, 1,0);

 f2_ST_SYST->cd();
 g2_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[5]+htitleC);
 h_ST_SYST->Add(g2_ST_SYST, h_ST_SYST, 1, 1);

 f1_TT_SYST->cd();
g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleC);
 h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);
 h_TT_SYST->Add(g1_TT_SYST, h_TT_SYST, 0.5,0);

f2_TT_SYST->cd();
 g2_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[3]+htitleC);
 h_TT_SYST->Add(g2_TT_SYST, h_TT_SYST, 1,1);


 htotMC_SYST->Add(htotMC_SYST, h_DY_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_VV_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_TTV_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_ST_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_TT_SYST, 1, 1);


htotMC_SYST->Draw("PE1same");
htotMC_SYST->SetFillStyle(1001);
 htotMC_SYST->SetFillColorAlpha(ColorRed, 1);
 htotMC_SYST->SetLineColor(ColorRed);

// htotMC_SYST->SetMarkerStyle(20);
// htotMC_SYST->SetMarkerSize(1);
// htotMC_SYST->SetMarkerColor(kBlack);
// htotMC_SYST->SetLineColor(kBlack);
// htotMC_SYST->SetLineWidth(1);


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
    leg->AddEntry(htotMC, " DY","F");
    leg->AddEntry(htotMC_SYST, " DY SYST","F");
  leg->Draw();


rap3->cd();

 hRatio = new TH1F(htitleC+"_"+SYST,"",nbin,xmin,xmax);
hRatio->Sumw2();
hRatio->Divide(htotMC_SYST, htotMC, 1, 1);
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
hRatio->GetYaxis()->SetTitle("SYST / Sim.");
hRatio->GetYaxis()->SetTitleColor(1);
hRatio->GetYaxis()->SetNdivisions(509);
hRatio->SetNdivisions(509,"XYZ");
hRatio->SetMinimum(0.5); 
hRatio->SetMaximum(1.5); 


        theoutputfile->cd();
        hRatio->Write();




// !! ------------- PAD 4 ----------------- !!//



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
  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleD);//ok
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleD);//ok
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleD);//ok
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleD);//ok 
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleD);//ok
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
  

// // *****************************************************************************
// // *****************************************************************************
 g1_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleD);//ok
 g2_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[1]+htitleD);//ok
  h_DY_SYST = new TH1F("h_DY_SYST","",nbin,xmin,xmax);

 
 g1_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[9]+htitleD);//ok
 g2_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[10]+htitleD);//ok
 g3_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[11]+htitleD);//ok
  h_VV_SYST = new TH1F("h_VV_SYST","",nbin,xmin,xmax);
 
g1_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[4]+htitleD);//ok
 g2_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[5]+htitleD);//ok
  h_ST_SYST = new TH1F("h_ST_SYST","",nbin,xmin,xmax);

 g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleD);//ok
 g2_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[3]+htitleD);//ok
  h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);

 g1_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[6]+htitleD);//ok
 g2_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[7]+htitleD);//ok 
 g3_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[8]+htitleD);//ok
  h_TTV_SYST = new TH1F("h_TTV_SYST","",nbin,xmin,xmax);
// *****************************************************************************


 htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 htotMC_SYST  = new TH1F("htotMC_SYST","",nbin,xmin,xmax);
 hRatioSYST_NOM = new TH1F("hRatioSYST_NOM","",nbin,xmin,xmax);

  htotMC_SYST->Sumw2();
  htotMC->Sumw2();
  hRatioSYST_NOM->Sumw2();
 

//------x);

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleD);
 g1_DY->Sumw2();
 h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 
 f2_DY->cd();
 
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleD);
 h_DY->Add(g2_DY, h_DY, 1,1);


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

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleD);
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 0.5,0);

f2_TT->cd();
 g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleD);
 h_TT->Add(g2_TT, h_TT, 1,1);


 htotMC->Add(htotMC, h_DY, 1, 1);
 htotMC->Add(htotMC, h_VV, 1, 1);
 htotMC->Add(htotMC, h_TTV, 1, 1);
 htotMC->Add(htotMC, h_ST, 1, 1);
 htotMC->Add(htotMC, h_TT, 1, 1);



 htotMC->Draw("PE1"); 
  htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(ColorNeutral, 1);
 htotMC->SetLineColor(ColorNeutral);
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
 htotMC->SetMinimum(htotMC->GetMinimum()/2+1); 
 htotMC->SetMaximum(htotMC->GetMaximum()*2); 

 f1_DY_SYST->cd();
 g1_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleD);
 g1_DY_SYST->Sumw2();
 h_DY_SYST = new TH1F("h_DY_SYST","",nbin,xmin,xmax);
 h_DY_SYST->Add(g1_DY_SYST, h_DY_SYST, 1,0);

 f2_DY_SYST->cd();
 g2_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[1]+htitleD);
 h_DY_SYST->Add(g2_DY_SYST, h_DY_SYST, 1,1);

 f1_VV_SYST->cd();
 g1_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[9]+htitleD);
 h_VV_SYST = new TH1F("h_VV_SYST","",nbin,xmin,xmax);
 h_VV_SYST->Add(g1_VV_SYST, h_VV_SYST, 1,0);

 f2_VV_SYST->cd();
 g2_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[10]+htitleD);
 h_VV_SYST->Add(g2_VV_SYST, h_VV_SYST, 1, 1);

 f3_VV_SYST->cd();
 g3_VV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[11]+htitleD);
 h_VV_SYST->Add(g3_VV_SYST, h_VV_SYST, 1, 1);

 f1_TTV_SYST->cd();
 g1_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[6]+htitleD);
 h_TTV_SYST = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV_SYST->Add(g1_TTV_SYST, h_TTV_SYST, 1,0);

 f2_TTV_SYST->cd();
 
 g2_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[7]+htitleD);
 h_TTV_SYST->Add(g2_TTV_SYST, h_TTV_SYST, 1, 1);

 f3_TTV_SYST->cd();
 g3_TTV_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[8]+htitleD);
 h_TTV_SYST->Add(g3_TTV_SYST, h_TTV_SYST, 1, 1);

 h_ST_SYST = new TH1F("h_ST_SYST","",nbin,xmin,xmax);
 f1_ST_SYST->cd();
 g1_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[4]+htitleD);
 h_ST_SYST->Add(g1_ST_SYST, h_ST_SYST, 1,0);

 f2_ST_SYST->cd();
 g2_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[5]+htitleD);
 h_ST_SYST->Add(g2_ST_SYST, h_ST_SYST, 1, 1);

 f1_TT_SYST->cd();
g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleD);
 h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);
 h_TT_SYST->Add(g1_TT_SYST, h_TT_SYST, 0.5,0);

f2_TT_SYST->cd();
 g2_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[3]+htitleD);
 h_TT_SYST->Add(g2_TT_SYST, h_TT_SYST, 1,1);


 htotMC_SYST->Add(htotMC_SYST, h_DY_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_VV_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_TTV_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_ST_SYST, 1, 1);
 htotMC_SYST->Add(htotMC_SYST, h_TT_SYST, 1, 1);


htotMC_SYST->Draw("PE1same");
htotMC_SYST->SetFillStyle(1001);
 htotMC_SYST->SetFillColorAlpha(ColorRed, 1);
 htotMC_SYST->SetLineColor(ColorRed);

// htotMC_SYST->SetMarkerStyle(20);
// htotMC_SYST->SetMarkerSize(1);
// htotMC_SYST->SetMarkerColor(kBlack);
// htotMC_SYST->SetLineColor(kBlack);
// htotMC_SYST->SetLineWidth(1);


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
    leg->AddEntry(htotMC, " DY","F");
    leg->AddEntry(htotMC_SYST, " DY SYST","F");
  leg->Draw();


rap3->cd();

 hRatio = new TH1F(htitleD+"_"+SYST,"",nbin,xmin,xmax);
hRatio->Sumw2();
hRatio->Divide(htotMC_SYST, htotMC, 1, 1);
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
hRatio->GetYaxis()->SetTitle("SYST / Sim.");
hRatio->GetYaxis()->SetTitleColor(1);
hRatio->GetYaxis()->SetNdivisions(509);
hRatio->SetNdivisions(509,"XYZ");
hRatio->SetMinimum(0.5); 
hRatio->SetMaximum(1.5); 


        theoutputfile->cd();
        hRatio->Write();





  TString name = htitleA+"_"+SYST+".pdf";
  c1->SaveAs("./SYST/"+name);


 f1_DY->Close();
 f2_DY->Close();
 f1_TT->Close(); 
 f2_TT->Close();  
 f1_ST->Close();  
 f2_ST->Close();  
 f1_TTV->Close(); 
 f2_TTV->Close(); 
 f3_TTV->Close(); 
 f1_VV->Close();  
 f2_VV->Close(); 
 f3_VV->Close();  

  f1_DY_SYST->Close();
 f2_DY_SYST->Close();
 f1_TT_SYST->Close(); 
 f2_TT_SYST->Close();  
 f1_ST_SYST->Close();  
 f2_ST_SYST->Close();  
 f1_TTV_SYST->Close(); 
 f2_TTV_SYST->Close(); 
 f3_TTV_SYST->Close(); 
 f1_VV_SYST->Close();  
 f2_VV_SYST->Close(); 
 f3_VV_SYST->Close();  

  theoutputfile->Close();
  delete theoutputfile;

  delete c1;
}

