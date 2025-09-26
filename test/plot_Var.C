#include <iostream>
#include <TROOT.h>
#include "TH1.h"
#include "TColor.h"
// #include "../MCWeights.h"
#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PlotCMS.h"

TCanvas * plot(int method,  TString Name, TString Year, TString Plots, bool SIGNAL)
{
int stati=0;
bool fit= 1;
bool logy=0;

bool Signal = SIGNAL;
float SumDataEvent = 0;
float SumMCEvent = 0;

TString EXTRA = "HT100_";
TString HTcut = "100";


TString ProdMC_MUMU = "MC_MUMU_2018_03_02_2025";
TString suffixMC_MUMU = "";


TString suffixDATA_MUMU = "_Corr";
TString Prod_MUMU = "DATA_MUMU_2018_03_02_2025";
TString SampleDATA_MUMU = "DoubleMuon_UL2018_MiniAODv2_GT36-v1";
TString SampleDATA_MUMUExtra  = "DoubleMuon_UL2018_MiniAODv2_GT36-v1"+suffixDATA_MUMU;

 if (Year == "2016PRE") 
  {


    ProdMC_MUMU = "MC_MUMU_2016PRE_30_03_2025";

    Prod_MUMU = "DATA_MUMU_2016PRE_30_03_2025";


    SampleDATA_MUMU = "MuonEG_Run2016-HIPM_UL2016_MiniAODv2";
    SampleDATA_MUMUExtra = "MuonEG_Run2016-HIPM_UL2016_MiniAODv2"+suffixDATA_MUMU;
  }
  
 if (Year == "2016POST") 
  {

    ProdMC_MUMU = "MC_MUMU_2016POST_30_03_2025";
    Prod_MUMU = "DATA_MUMU_2016POST_30_03_2025";
    SampleDATA_MUMU = "DoubleMuon_Run2016-UL2016_MiniAODv2";
    SampleDATA_MUMUExtra = "DoubleMuon_Run2016-UL2016_MiniAODv2"+suffixDATA_MUMU;
  }
 if (Year == "2017") 
  {

    ProdMC_MUMU = "MC_MUMU_2017_30_03_2025";
    Prod_MUMU = "DATA_MUMU_2017_30_03_2025";

    SampleDATA_MUMU = "DoubleMuon_Run2017-UL2017_MiniAODv2";
    SampleDATA_MUMUExtra = "DoubleMuon_Run2017-UL2017_MiniAODv2"+suffixDATA_MUMU;
  }
 if (Year == "2018") 
  {

    ProdMC_MUMU = "MC_MUMU_2018_03_02_2025";
    Prod_MUMU = "DATA_MUMU_2018_03_02_2025";
    SampleDATA_MUMU = "DoubleMuon_UL2018_MiniAODv2_GT36-v1";
    SampleDATA_MUMUExtra = "DoubleMuon_UL2018_MiniAODv2_GT36-v1"+suffixDATA_MUMU;
  }

TFile* f1_Data_mumu  = new TFile("../../"+Prod_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_"+SampleDATA_MUMUExtra+".root");// remplacer par la prédi dans SR de MUMU
 
 TFile* f1_DY_MUMU  = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8"+suffixMC_MUMU+".root");
 TFile* f2_DY_MUMU  = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"+suffixMC_MUMU+".root");
 TFile* f1_TT_MUMU  = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"+suffixMC_MUMU+".root");
 TFile* f2_TT_MUMU  = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"+suffixMC_MUMU+".root");
 TFile* f1_ST_MUMU  = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"+suffixMC_MUMU+".root");
 TFile* f2_ST_MUMU  = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"+suffixMC_MUMU+".root");
//  TFile* f3_ST  = new TFile("../../"+ProdMC+"/histofile_"+"_OS_2p4_ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
 //  TFile* f4_ST  = new TFile("../../"+ProdMC+"/histofile_"+"_OS_2p4_ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_TTV_MUMU = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8"+suffixMC_MUMU+".root");
 TFile* f2_TTV_MUMU = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8"+suffixMC_MUMU+".root");
 TFile* f3_TTV_MUMU = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8"+suffixMC_MUMU+".root");
 TFile* f1_VV_MUMU  = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"+suffixMC_MUMU+".root");
 TFile* f2_VV_MUMU  = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8"+suffixMC_MUMU+".root");
 TFile* f3_VV_MUMU  = new TFile("../../"+ProdMC_MUMU+"/histofile_"+EXTRA+"DM_OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8"+suffixMC_MUMU+".root");



TString MSUON[3] = {"200","300","500"};
  TString MNEU[3] = {"180","200","300"};
  TString CTAU[3] = {"100","100","100"};

//signal
 TFile* f1_LLP = new TFile("../../Signal_2018_L1/histofile_"+EXTRA+"DM_OS_2p4_RPV_2018_smu"+MSUON[0]+"_neu"+MNEU[0]+"_ctau"+CTAU[0]+".root");
 TFile* f2_LLP = new TFile("../../Signal_2018_L1/histofile_"+EXTRA+"DM_OS_2p4_RPV_2018_smu"+MSUON[1]+"_neu"+MNEU[1]+"_ctau"+CTAU[1]+".root");
 TFile* f3_LLP = new TFile("../../Signal_2018_L1/histofile_"+EXTRA+"DM_OS_2p4_RPV_2018_smu"+MSUON[2]+"_neu"+MNEU[2]+"_ctau"+CTAU[2]+".root");

 TString DATAFILE[1] = {SampleDATA_MUMU+"_"
};



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

};


TString LLPFILE[3] = {"RPV_2018_smu"+MSUON[0]+"_neu"+MNEU[0]+"_ctau"+CTAU[0]+"_",
                  "RPV_2018_smu"+MSUON[1]+"_neu"+MNEU[1]+"_ctau"+CTAU[1]+"_",
                  "RPV_2018_smu"+MSUON[2]+"_neu"+MNEU[2]+"_ctau"+CTAU[2]+"_"
};

    TString ytitle = "a.u"; 
    TString htitleC = "hData_EVT34_1Vtx_BDTvtx";

    int nbin = 40; 
    float xmin = 0;
    float xmax =  40;
    TString HeaderC = "A";
    TString HeaderCbis = "Abis";
    TString HeaderNVtx = "k Vtx";
    TString xtitle = "var";

  //-----------------------------------------------------------//
  // ABCD using Hemipt and Tight+looseWP 
  //-----------------------------------------------------------//
  int Method = method;


if (Method == 0)
  {
     htitleC = "Hemisphere_leadingpt_";//

    nbin = 30; 
    xmin = 0;
    xmax = 300;
HeaderC = ">= 1 tight";
HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderNVtx = "2 Vertices";
    xtitle = "Hemi Leading pt [GeV]";
  }

    if (Method == 1)
  {


    htitleC = "Hemisphere_subleadingpt_";//

    nbin = 30; 
    xmin = 0;
    xmax = 300;
    HeaderC = ">= 1 tight";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderNVtx = "2 Vertices";
    xtitle = "Hemi SubLeading pt [GeV]";
      }

// xsec in pb


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

  TPad* pad1 = new TPad("pad1","This is pad1",0.01,0.05,0.95,0.99,21);

float bottomMargin = 0.05;
bottomMargin = 0.15;
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.07);
   pad1->SetBottomMargin(bottomMargin);
   pad1->SetRightMargin(0.04);
   pad1->SetLeftMargin(0.16);



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
gROOT->SetBatch(kTRUE);


 TH1F* htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 TH1F* htotData  = new TH1F("htotData","",nbin,xmin,xmax);

  htotData->Sumw2();
  htotMC->Sumw2();


TH1F* g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleC);

 TH1F* g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);//ok
 TH1F* g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleC);//ok
 TH1F*  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);


 TH1F* g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleC);//ok
 TH1F* g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleC);//ok
 TH1F* g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleC);//ok
 TH1F*  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);

TH1F* g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleC);//ok
 TH1F* g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleC);//ok
 TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 TH1F* g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);//ok
TH1F* g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleC);//ok
 TH1F*  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 TH1F* g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleC);//ok
 TH1F* g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleC);//ok
 TH1F* g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleC);//ok
 TH1F*  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);

 TH1F* g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleC);//ok
 TH1F* h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 TH1F* g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleC);//ok
 TH1F* h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 TH1F* g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleC);//ok
 TH1F* h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);

// *****************************************************************************

 pad1->cd();

    f1_DY_MUMU->cd();
    g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);
    h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
    h_DY->Add(g1_DY, h_DY, 1,0);

    f2_DY_MUMU->cd();
    g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleC);
    h_DY->Add(g2_DY, h_DY, 1, 1);

    f1_VV_MUMU->cd();
    g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleC);
    h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
    h_VV->Add(g1_VV, h_VV, 1,0);

    f2_VV_MUMU->cd();
    g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleC);
    h_VV->Add(g2_VV, h_VV, 1, 1);

    f3_VV_MUMU->cd();
    g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleC);
    h_VV->Add(g3_VV, h_VV, 1, 1);

    f1_TTV_MUMU->cd();
    g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleC);
    h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
    h_TTV->Add(g1_TTV, h_TTV, 1,0);

    f2_TTV_MUMU->cd();
    g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleC);
    h_TTV->Add(g2_TTV, h_TTV, 1, 1);

    f3_TTV_MUMU->cd();
    g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleC);
    h_TTV->Add(g3_TTV, h_TTV, 1, 1);

    h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
    f1_ST_MUMU->cd();
    g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleC);
    h_ST->Add(g1_ST, h_ST, 1,0);

    f2_ST_MUMU->cd();
    g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleC);
    h_ST->Add(g2_ST, h_ST, 1, 1);

    f1_TT_MUMU->cd();
    g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);

    h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
    h_TT->Add(g1_TT, h_TT, 1,0);

    f2_TT_MUMU->cd();
    g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleC);
    h_TT->Add(g2_TT, h_TT, 1,1);

h_DY->Scale(1./h_DY->Integral(0,-1));
h_VV->Scale(1./h_VV->Integral(0,-1));
h_TTV->Scale(1./h_TTV->Integral(0,-1));
h_ST->Scale(1./h_ST->Integral(0,-1));
h_TT->Scale(1./h_TT->Integral(0,-1));

 htotData = new TH1F("htotData","",nbin,xmin,xmax);

f1_Data_mumu->cd();
g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleC);
htotData->Add(g1_Data_emu, htotData, 1, 0);

htotData->Draw("PE1same");
htotData->SetMarkerStyle(20);
htotData->SetMarkerSize(0);
htotData->SetMarkerColor(kWhite);
htotData->SetLineColor(kWhite);
htotData->SetLineWidth(1);
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
htotData->Scale(1./htotData->Integral(0,-1));
// htotData->GetXaxis()->SetRangeUser(0,20);
if (logy)
    {
        htotData->SetMinimum(1); 
        htotData->SetMaximum(htotData->GetMaximum()*100); 
    }
else 
    {
        htotData->SetMinimum(0); 
        htotData->SetMaximum(0.6); 
    }



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
    h1_LLP->SetLineStyle(2);
    h1_LLP->SetLineWidth(2);

    h2_LLP->Draw("HEsame"); 
    h2_LLP->SetLineColor(kRed-2);
    h2_LLP->SetLineStyle(3);
    h2_LLP->SetLineWidth(2);

    h3_LLP->Draw("HEsame");
    h3_LLP->SetLineColor(kRed-3);
    h3_LLP->SetLineStyle(4);
    h3_LLP->SetLineWidth(2);

    h1_LLP->Scale(1./h1_LLP->Integral(0,-1));
    h2_LLP->Scale(1./h2_LLP->Integral(0,-1));
    h3_LLP->Scale(1./h3_LLP->Integral(0,-1));

    h_DY->Draw("HEsame"); 
    h_DY->SetLineColor(ColorBlue);
    h_DY->SetLineStyle(1);
    h_DY->SetLineWidth(2);

    h_TT->Draw("HEsame");
    h_TT->SetLineColor(ColorRed);
    h_TT->SetLineStyle(1);
    h_TT->SetLineWidth(2);


  float LEGY1 = 0.50;
  float LEGY2 = 0.89;
  float legsize = 0.06;
  float legsize2 = 0.04;
  LEGY1 = 0.6;
  LEGY2 = 0.89;
  leg = new TLegend(0.4,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
  // leg->AddEntry(htotData, "#mu#mu data SR","PE1");
  leg->AddEntry(h_DY, "DY","L");
  leg->AddEntry(h_TT, "t#bar{t}","L");
  leg->AddEntry(h1_LLP," m^{"+CTAU[0]+" cm}_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[0]+" ("+MNEU[0]+") GeV","L");
  leg->AddEntry(h2_LLP," m^{"+CTAU[1]+" cm}_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[1]+" ("+MNEU[1]+") GeV","L");
  leg->AddEntry(h3_LLP," m^{"+CTAU[2]+" cm}_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSUON[2]+" ("+MNEU[2]+") GeV","L");
  leg->Draw();


PlotCMSv3(pad1,Year,false);

// !! --------------------------
// !! --------------------------
// !! --------------------------

// *****************************************************************************
// *****************************************************************************

  TString namele =  Name+"_"+Plots;
  namele += "_"+Year;
  c1->SaveAs(namele+".pdf");

  return c1;
}