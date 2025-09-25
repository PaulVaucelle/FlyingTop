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

 TFile* f1_DY  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f2_DY  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f1_TT  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_TT  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_ST  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_ST  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f3_ST  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f4_ST  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_TTV = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8.root");
 TFile* f2_TTV = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f3_TTV = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8.root");
 TFile* f1_VV  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_VV  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 TFile* f3_VV  = new TFile("../../MC_MUMU_2018_03_02_2025/histofile_DM_OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 
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
                  "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_",
                  // "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_",
                  // "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_"

};

// extension is <SYST>Up or <SYST>Down
TString EXTRA = "";
  if (SYST == "JECUp") EXTRA = "_JECUp";
  else if (SYST == "JECDown") EXTRA = "_JECDown";
  else if (SYST == "JERUp" ) EXTRA = "_JERUp";
  else if (SYST == "JERDown" ) EXTRA = "_JERDown";
  else if (SYST == "RoccorDown") EXTRA = "_RoccorDown";


  TFile* f1_DY_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_"+extension+".root");
 TFile* f2_DY_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_"+extension+".root");
 TFile* f1_TT_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f2_TT_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f1_ST_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f2_ST_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
//  TFile* f3_ST_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
//  TFile* f4_ST_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f1_TTV_SYST = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8_"+extension+".root");
 TFile* f2_TTV_SYST = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8_"+extension+".root");
 TFile* f3_TTV_SYST = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8_"+extension+".root");
 TFile* f1_VV_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_"+extension+".root");
 TFile* f2_VV_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_"+extension+".root");
 TFile* f3_VV_SYST  = new TFile("../../MC_MUMU_2018_03_02_2025"+EXTRA+"/histofile_DM_OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_"+extension+".root");
 


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
                  "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_"+extension+"_",
                  // "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+"_",
                  // "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_"+extension+"_"
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
        htitleA = "hData_CRlowlowpt_1Vtx_Mass_";
        htitleB = "hData_CRlooselowlowpt_1Vtx_Mass_";
        htitleC = "hData_CRTight_1Vtx_Mass_";
        htitleD = "hData_CRLoose_1Vtx_Mass_";

        nbin = 25; 
        xmin = 0;
        xmax =  100;
        HeaderA = "T + 20<pt_{i}<80 && pt_{j}>80";
        HeaderB = "L + 20<pt<80";
        HeaderC = "T + 20<pt_{i}<80 && pt_{j}>80";
        HeaderD = "L + 20<pt<80";

        HeaderNVtx = "1 Vtx"; 
        xtitle = "Vtx Mass";
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
        HeaderA = "T + 20<pt_{i}<80 && pt_{j}>80";
        HeaderB = "L + 20<pt<80";
        HeaderC = "T + 20<pt_{i}<80 && pt_{j}>80";
        HeaderD = "L + 20<pt<80";


        HeaderNVtx = "1 Vtx";
        xtitle = "SumtrackWeight";
      }

    if (Method == 2)
      {
        htitleA = "Quality_LT_STW_1Vtx_A_";
        htitleB = "Quality_LT_STW_1Vtx_B_";
        htitleC = "Quality_LT_STW_1Vtx_C_";
        htitleD = "Quality_LT_STW_1Vtx_D_";

        nbin = 40; 
        xmin = 0;
        xmax = 40;
        HeaderA = "T + L_{T}c80";
        HeaderB = "L + L_{T}<80";
        HeaderC = "T + L_{T}>80";
        HeaderD = "L + L_{T}>80";


        HeaderNVtx = "1 Vtx";
        xtitle = "SumtrackWeight";
      }

    if (Method == 3)
      {
        htitleA = "hData_CRtightlowlowpt_TLVtx_Mass_";//
        htitleB = "hData_CRlooselooselowlowpt_TLVtx_Mass_";//
        htitleC = "hData_CRtighthighpt_2Vtx_Mass_";//
        htitleD = "hData_CRlooselooselowpt_2Vtx_Mass_";//

        nbin = 25; 
        xmin = 0;
        xmax =  100;
        HeaderA = "TT+TL + 20<pt_{i}<80 && pt_{j}>80";
        HeaderB = "LL + 20<pt<80";
        HeaderC = "TT+TL + 20<pt_{i}<80 && pt_{j}>80";
        HeaderD = "LL + 20<pt<80";

        HeaderNVtx = "2 Vtx"; 
        xtitle = "Vtx Mass";
      }


    if (Method == 4)
      {
        htitleA = "hData_CRtightlowlowpt_TLVtxAll_Mass_";//
        htitleB = "hData_CRlooselooselowlowpt_TLVtxAll_Mass_";//
        htitleC = "hData_CRtighthighpt_2VtxAll_Mass_";//
        htitleD = "hData_CRlooselooselowpt_2VtxAll_Mass_";//

        nbin = 25; 
        xmin = 0;
        xmax = 100;
        HeaderA = "TT+TL + 20<pt_{i}<80 && pt_{j}>80";
        HeaderB = "LL + 20<pt<80";
        HeaderC = "TT+TL + 20<pt_{i}<80 && pt_{j}>80";
        HeaderD = "LL + 20<pt<80";


        HeaderNVtx = "2 Vtx";
        xtitle = "Vtx Mass";
      }

    if (Method == 5)
      {
        htitleA = "hData_CRtightlowlowpt_TLVtx_SumtrackWeight_";//
        htitleB = "hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight_";//
        htitleC = "hData_CRtighthighpt_2Vtx_SumtrackWeight_";//
        htitleD = "hData_CRlooselooselowpt_2Vtx_SumtrackWeight_";//

        nbin = 40; 
        xmin = 0;
        xmax = 40;
        HeaderA = "TT+TL + 20<pt_{i}<80 && pt_{j}>80";
        HeaderB = "LL + 20<pt<80";
        HeaderC = "TT+TL + 20<pt_{i}<80 && pt_{j}>80";
        HeaderD = "LL + 20<pt<80";


        HeaderNVtx = "2 Vtx";
        xtitle = "SumtrackWeight";
      }

    // !!----------------------
    if (Method == 6)
      {
    htitleA = "hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight_";//
    htitleB = "hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight_";//
    htitleC = "hData_CRtighthighpt_2VtxAll_SumtrackWeight_";//
    htitleD = "hData_CRlooselooselowpt_2VtxAll_SumtrackWeight_";//

    nbin = 40; 
    xmin = 0;
    xmax = 40;
    HeaderA = "TT+TL + 20<pt_{i}<80 && pt_{j}>80";
    HeaderB = "LL + 20<pt<80";
    HeaderC = "TT+TL + 20<pt_{i}<80 && pt_{j}>80";
    HeaderD = "LL + 20<pt<80";


    HeaderNVtx = "2 VtxAll";
    xtitle = "SumtrackWeight";
      }


  // !!----------------------
    if (Method == 7)
      {
        htitleA = "Quality_LT_STW_2Vtx_A_";//
        htitleB = "Quality_LT_STW_2Vtx_B_";//
        htitleC = "Quality_LT_STW_2Vtx_C_";//
        htitleD = "Quality_LT_STW_2Vtx_D_";//

        nbin = 40; 
        xmin = 0;
        xmax = 40;
        HeaderA = "TT+TL + L_{T}<80";
        HeaderB = "LL + L_{T}<80";
        HeaderC = "TT+TL + L_{T}>80";
        HeaderD = "LL + L_{T}>80";


        HeaderNVtx = "2 Vtx";
        xtitle = "SumtrackWeight";
      }

  // !!----------------------
    if (Method == 8)
      {
        htitleA = "Quality_LT_STW_2VtxAll_A_";//
        htitleB = "Quality_LT_STW_2VtxAll_B_";//
        htitleC = "Quality_LT_STW_2VtxAll_C_";//
        htitleD = "Quality_LT_STW_2VtxAll_D_";//

        nbin = 40; 
        xmin = 0;
        xmax = 40;
        HeaderA = "TT+TL + L_{T}<80";
        HeaderB = "LL + L_{T}<80";
        HeaderC = "TT+TL + L_{T}>80";
        HeaderD = "LL + L_{T}>80";


        HeaderNVtx = "2 VtxAll";
        xtitle = "SumtrackWeight";
      }


    // !!----------------------

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
//   TH1F* g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);//ok
//  TH1F* g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);//ok
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
//  TH1F* g3_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[12]+htitleA);//ok
//  TH1F* g4_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[13]+htitleA);//ok
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

//   f3_ST->cd();
//  g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);
//  h_ST->Add(g3_ST, h_ST, 1, 1);

//   f4_ST->cd();
//  g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);
//  h_ST->Add(g4_ST, h_ST, 1, 1);


 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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

//    f3_ST_SYST->cd();
//  g3_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[12]+htitleA);
//  h_ST_SYST->Add(g3_ST_SYST, h_ST_SYST, 1, 1);

//   f4_ST_SYST->cd();
//  g4_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[13]+htitleA);
//  h_ST_SYST->Add(g4_ST_SYST, h_ST_SYST, 1, 1);

 f1_TT_SYST->cd();
g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleA);
 h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);
 h_TT_SYST->Add(g1_TT_SYST, h_TT_SYST, 1,0);

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

//   f3_ST->cd();
//  g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleB);
//  h_ST->Add(g3_ST, h_ST, 1, 1);

//   f4_ST->cd();
//  g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleB);
//  h_ST->Add(g4_ST, h_ST, 1, 1);


 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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

//     f3_ST_SYST->cd();
//  g3_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[12]+htitleB);
//  h_ST_SYST->Add(g3_ST_SYST, h_ST_SYST, 1, 1);

//   f4_ST_SYST->cd();
//  g4_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[13]+htitleB);
//  h_ST_SYST->Add(g4_ST_SYST, h_ST_SYST, 1, 1);

 f1_TT_SYST->cd();
g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleB);
 h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);
 h_TT_SYST->Add(g1_TT_SYST, h_TT_SYST, 1,0);

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

//   f3_ST->cd();
//  g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleC);
//  h_ST->Add(g3_ST, h_ST, 1, 1);

//   f4_ST->cd();
//  g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleC);
//  h_ST->Add(g4_ST, h_ST, 1, 1);

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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

//    f3_ST_SYST->cd();
//  g3_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[12]+htitleC);
//  h_ST_SYST->Add(g3_ST_SYST, h_ST_SYST, 1, 1);

//   f4_ST_SYST->cd();
//  g4_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[13]+htitleC);
//  h_ST_SYST->Add(g4_ST_SYST, h_ST_SYST, 1, 1);

 f1_TT_SYST->cd();
g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleC);
 h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);
 h_TT_SYST->Add(g1_TT_SYST, h_TT_SYST, 1,0);

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

//    f3_ST->cd();
//  g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleD);
//  h_ST->Add(g3_ST, h_ST, 1, 1);

//   f4_ST->cd();
//  g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleD);
//  h_ST->Add(g4_ST, h_ST, 1, 1);

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleD);
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

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

//    f3_ST_SYST->cd();
//  g3_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[12]+htitleD);
//  h_ST_SYST->Add(g3_ST_SYST, h_ST_SYST, 1, 1);

//   f4_ST_SYST->cd();
//  g4_ST_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[13]+htitleD);
//  h_ST_SYST->Add(g4_ST_SYST, h_ST_SYST, 1, 1);

 f1_TT_SYST->cd();
g1_TT_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[2]+htitleD);
 h_TT_SYST = new TH1F("h_TT_SYST","",nbin,xmin,xmax);
 h_TT_SYST->Add(g1_TT_SYST, h_TT_SYST, 1,0);

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


rap4->cd();

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
//  f3_ST->Close();  
//  f4_ST->Close(); 
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
//  f3_ST_SYST->Close();  
//  f4_ST_SYST->Close();  
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