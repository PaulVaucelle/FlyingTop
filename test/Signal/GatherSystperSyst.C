



#include <iostream>
#include <TROOT.h>
#include "TH1.h"
#include "TColor.h"

void plot(TString Prod, int method, TString Year, TString SYST , TString Plots, TString MSMU, TString MNEU, TString CTAU)
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
 
TFile* f1_LLP = new TFile("../../"+Prod+"/histofile_"+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau001.root");

 TString MCFILE[1] = {"RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau001_"};

// extension is <SYST>Up or <SYST>Down
TString EXTRA = "";
  if (SYST == "JECUp") EXTRA = "_JECUp";
  else if (SYST == "JECDown") EXTRA = "_JECDown";
  else if (SYST == "JERUp" ) EXTRA = "_JERUp";
  else if (SYST == "JERDown" ) EXTRA = "_JERDown";
  else if (SYST == "RoccorDown") EXTRA = "_RoccorDown";

  TFile* f1_LLP_SYST = new TFile("../../"+Prod+"/histofile_"+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSMU+"_neu"+MNEU+"_ctau"+CTAU+".root");

TString MCFILE_SYST[1] = {"RPV_"+Yearcor+"_smu"+MSMU+"_neu"+MNEU+"_ctau"+CTAU+"_"+extension+"_"};


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
    //-----------------------------------------------------------//
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

 TH1F* g1_LLP = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);//ok
 TH1F* g1_LLP_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleA);//ok

// *****************************************************************************
 pad1->cd();

 TH1F* htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 TH1F* htotMC_SYST  = new TH1F("htotMC_SYST","",nbin,xmin,xmax);

  htotMC_SYST->Sumw2();
  htotMC->Sumw2();

 f1_LLP->cd();
 g1_LLP = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
 g1_LLP->Sumw2();

htotMC->Add(htotMC, g1_LLP, 0, 1);

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

 f1_LLP_SYST->cd();
 g1_LLP_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleA);
 g1_LLP_SYST->Sumw2();

htotMC_SYST->Add(htotMC_SYST, g1_LLP_SYST, 0, 1);
htotMC_SYST->Draw("PE1same");
htotMC_SYST->SetFillStyle(1001);
htotMC_SYST->SetFillColorAlpha(ColorRed, 1);
htotMC_SYST->SetLineColor(ColorRed);

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
    leg->AddEntry(htotMC, " LLP","F");
    leg->AddEntry(htotMC_SYST, " LLP SYST","F");
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

 g1_LLP = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);//ok
 g1_LLP_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleB);//ok
// *****************************************************************************
 htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 htotMC_SYST  = new TH1F("htotMC_SYST","",nbin,xmin,xmax);

  htotMC_SYST->Sumw2();
  htotMC->Sumw2();

 f1_LLP->cd();
 g1_LLP = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);
 g1_LLP->Sumw2();
htotMC->Add(htotMC, g1_LLP, 0, 1);

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

 f1_LLP_SYST->cd();
 g1_LLP_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleB);
 g1_LLP_SYST->Sumw2();

htotMC_SYST->Add(htotMC_SYST, g1_LLP_SYST, 0, 1);
htotMC_SYST->Draw("PE1same");
htotMC_SYST->SetFillStyle(1001);
 htotMC_SYST->SetFillColorAlpha(ColorRed, 1);
 htotMC_SYST->SetLineColor(ColorRed);

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
    leg->AddEntry(htotMC, " LLP","F");
    leg->AddEntry(htotMC_SYST, " LLP SYST","F");
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

 g1_LLP = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);//ok

 htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 htotMC_SYST  = new TH1F("htotMC_SYST","",nbin,xmin,xmax);

  htotMC_SYST->Sumw2();
  htotMC->Sumw2();

 f1_LLP->cd();
 g1_LLP = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);
 g1_LLP->Sumw2();
htotMC->Add(htotMC, g1_LLP, 0, 1);

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

 f1_LLP_SYST->cd();
 g1_LLP_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleC);
 g1_LLP_SYST->Sumw2();

htotMC_SYST->Add(htotMC_SYST, g1_LLP_SYST, 0, 1);

htotMC_SYST->Draw("PE1same");
htotMC_SYST->SetFillStyle(1001);
 htotMC_SYST->SetFillColorAlpha(ColorRed, 1);
 htotMC_SYST->SetLineColor(ColorRed);

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

 g1_LLP = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleD);//ok
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleD);//ok
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

// // *****************************************************************************
// // *****************************************************************************
 g1_LLP_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleD);//ok
 g2_DY_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[1]+htitleD);//ok
  h_DY_SYST = new TH1F("h_DY_SYST","",nbin,xmin,xmax);

// *****************************************************************************
 htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
 htotMC_SYST  = new TH1F("htotMC_SYST","",nbin,xmin,xmax);

  htotMC_SYST->Sumw2();
  htotMC->Sumw2();

 f1_LLP->cd();
 g1_LLP = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleD);
 g1_LLP->Sumw2();

htotMC->Add(htotMC, g1_LLP, 0, 1);

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

 f1_LLP_SYST->cd();
 g1_LLP_SYST = (TH1F*)gROOT->FindObject(MCFILE_SYST[0]+htitleD);
 g1_LLP_SYST->Sumw2();

htotMC_SYST->Add(htotMC_SYST, g1_LLP_SYST, 0, 1);

htotMC_SYST->Draw("PE1same");
htotMC_SYST->SetFillStyle(1001);
 htotMC_SYST->SetFillColorAlpha(ColorRed, 1);
 htotMC_SYST->SetLineColor(ColorRed);

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

  TString name = htitleA+"_RPV_"+Yearcor+"_smu"+MSMU+"_neu"+MNEU+"_ctau"+CTAU+"_"+SYST+".pdf";
  c1->SaveAs("./SYST/"+name);
 f1_LLP->Close();
 f1_LLP_SYST->Close();
  theoutputfile->Close();
  delete theoutputfile;
  delete c1;
}

