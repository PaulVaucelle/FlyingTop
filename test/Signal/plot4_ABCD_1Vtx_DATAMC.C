#include <iostream>
#include <TROOT.h>
#include "TH1.h"
#include "TColor.h"

TCanvas * plot(int method, TString Prod, TString Name, TString Dmode, TString MSmuon, TString MNeu, TString ctau, TString syst, TString Year)
{
int stati=0;
bool fit= 1;
bool logy=0;

// number of vertices:
//$$
  int nvtx = 2;
//$$
  float hmin = 0.1; // cannot be 0 for logy=1
//$$  float hmax = 1E5;	         // for eta<2.4 pt>80 or CRlowpt
  float hmax = 100;	         // for eta<2.4 pt>80 or CRlowpt
  float hmaxBD = 100;	         // for eta<2.4 pt>80 or CRlowpt
  if ( nvtx == 1 ) {
    hmax = 100;   // for eta<2.4 pt>80  
    hmaxBD = 50;
  }
  if ( nvtx == 2 ) { 
    hmax = 100;   // for eta<2.4 pt>80  
    hmaxBD = 1E6;
//     hmax = 1E6;   // for eta<2.4 pt>80  
//     hmaxBD = 1E9;
  }
  // TString Prod = "DATA_EMU_2017";
  //ABCD_EMu2018A 094
  // ABCD_2018A 095
  float ReScaleXS = 1.;//0.0812948*2
  // TString Dmode = "DM";//DM or SM
  TString EXTRA = "HT100_";
TString HTcut = "100";
 TString Yearcor = Year;
 if (Year == "2016PRE") Yearcor = "2016preVFP";
 if (Year == "2016POST") Yearcor = "2016";

//signal
 TFile* f1_LLP = new TFile("../../"+Prod+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau001.root");
 TFile* f2_LLP = new TFile("../../"+Prod+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau003.root");
 TFile* f3_LLP = new TFile("../../"+Prod+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau010.root");
 TFile* f4_LLP = new TFile("../../"+Prod+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau030.root");
 TFile* f5_LLP = new TFile("../../"+Prod+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau100.root");
 TFile* f6_LLP = new TFile("../../"+Prod+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau300.root");
 TFile* f7_LLP = new TFile("../../"+Prod+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau1000.root");

TString FILE[7] = {"RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau001_",
                  "RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau003_",
                  "RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau010_",
                  "RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau030_",
                  "RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau100_",
                  "RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau300_",
                  "RPV_"+Yearcor+"_smu"+MSmuon+"_neu"+MNeu+"_ctau1000_"};


//  if ( nvtx == 1 ) xtitle = "vertex BDT score"; 
 TString ytitle = "Events"; 
 TString HeaderCMS = "CMS";

 if (Year == "2016") HeaderCMS = "2016                  36.3 fb^{-1} (13 TeV)";
 if (Year == "2017") HeaderCMS = "2017                  41.5 fb^{-1} (13 TeV)";
 if (Year == "2018") HeaderCMS = "2018                  59.8 fb^{-1} (13 TeV)";
    TString htitleA = "hData_EVT34_1Vtx_BDTvtx";
    TString htitleB = "hData_NoEVT34_1Vtx_BDTvtx";
    TString htitleC = "hData_EVT12_1Vtx_BDTvtx";
    TString htitleD = "hData_NoEVT12_1Vtx_BDTvtx";
    int nbin = 25; 
    float xmin = -1.0;
    float xmax =  1.0;
    TString HeaderA = "A";
TString HeaderAbis = "Abis";
    TString HeaderB = "B";
TString HeaderBbis = "Bbis";
    TString HeaderC = "C";
TString HeaderCbis = "Cbis";
    TString HeaderD = "D";
TString HeaderDbis = "Dbis";
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
    HeaderA = "Tight";
    HeaderAbis = "30<pt_{i}<"+HTcut+" & pt_{j}>"+HTcut;
    HeaderB = "Loose";
    HeaderBbis = "30<pt<"+HTcut+"";
    HeaderC = "Tight";
    HeaderCbis = "30<pt_{i}<"+HTcut+" & pt_{j}>"+HTcut+"";
    HeaderD = "Loose";
    HeaderDbis = "30<pt<"+HTcut+"";

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
    HeaderA = "Tight";
    HeaderB = "Loose";
    HeaderC = "Tight";
    HeaderD = "Loose";
    HeaderAbis = "30<pt_{i}<"+HTcut+" & pt_{j}>"+HTcut+"";
    HeaderCbis = "30<pt_{i}<"+HTcut+" & pt_{j}>"+HTcut+"";
    HeaderBbis = "30<pt<"+HTcut+"";
    HeaderDbis = "30<pt<"+HTcut+"";

    HeaderNVtx = "1 Vtx";
    xtitle = "Sum of track weights at Vtx";
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
    HeaderA = "Tight";
    HeaderB = "Loose";
    HeaderC = "Tight";
    HeaderD = "Loose";

    HeaderAbis = "L_{T}<80";
    HeaderCbis = "L_{T}>80";
    HeaderBbis = "L_{T}<80";
    HeaderDbis = "L_{T}>80";

    HeaderNVtx = "1 Vtx";
    xtitle = "Sum of track weights at Vtx";
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


c1->cd();


TPad* pA = new TPad("pA","This is pad1",0.01,0.5,0.5,0.99,21);
TPad* pB = new TPad("pB","This is pad2",0.01,0.01,0.5,0.49,21);
TPad* pC = new TPad("pC","This is pad3",0.51,0.5,0.99,0.99,21);
TPad* pD = new TPad("pD","This is pad4",0.51,0.01,0.99,0.49,21);

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


gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleFont(62);
gStyle->SetTitleColor(1);
gStyle->SetTitleTextColor(1);
gStyle->SetTitleFillColor(10);
gStyle->SetTitleFontSize(0.05);
gROOT->SetBatch(kTRUE);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(0.09);
gStyle->SetOptStat(stati);
gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);


 TH1F* g1_LLP = (TH1F*)gROOT->FindObject(FILE[0]+htitleA);

 TH1F* h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);

 TH1F* g2_LLP = (TH1F*)gROOT->FindObject(FILE[1]+htitleA);
 TH1F* h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);

 TH1F* g3_LLP = (TH1F*)gROOT->FindObject(FILE[2]+htitleA);
 TH1F* h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);

 TH1F* g4_LLP = (TH1F*)gROOT->FindObject(FILE[3]+htitleA);
 TH1F* h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);

 TH1F* g5_LLP = (TH1F*)gROOT->FindObject(FILE[4]+htitleA);
 TH1F* h5_LLP = new TH1F("h5_LLP","",nbin,xmin,xmax);

 TH1F* g6_LLP = (TH1F*)gROOT->FindObject(FILE[5]+htitleA);
 TH1F* h6_LLP = new TH1F("h6_LLP","",nbin,xmin,xmax);

 TH1F* g7_LLP = (TH1F*)gROOT->FindObject(FILE[6]+htitleA);
 TH1F* h7_LLP = new TH1F("h7_LLP","",nbin,xmin,xmax);


 pA->cd();


//------x);

 f1_LLP->cd();

 g1_LLP = (TH1F*)gROOT->FindObject(FILE[0]+htitleA);
 
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax); 
 h1_LLP->Add(g1_LLP, h1_LLP, 1,0);


 f2_LLP->cd();

 g2_LLP = (TH1F*)gROOT->FindObject(FILE[1]+htitleA);
 
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

 f3_LLP->cd();

 g3_LLP = (TH1F*)gROOT->FindObject(FILE[2]+htitleA);
 
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

  f4_LLP->cd();

   g4_LLP = (TH1F*)gROOT->FindObject(FILE[3]+htitleA);
  
   h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
  h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

    f5_LLP->cd();

   g5_LLP = (TH1F*)gROOT->FindObject(FILE[4]+htitleA);
  
   h5_LLP = new TH1F("h5_LLP","",nbin,xmin,xmax);
  h5_LLP->Add(g5_LLP, h5_LLP, 1,0);

    f6_LLP->cd();

   g6_LLP = (TH1F*)gROOT->FindObject(FILE[5]+htitleA);
  
   h6_LLP = new TH1F("h6_LLP","",nbin,xmin,xmax);
  h6_LLP->Add(g6_LLP, h6_LLP, 1,0);

    f7_LLP->cd();

   g7_LLP = (TH1F*)gROOT->FindObject(FILE[6]+htitleA);
  
   h7_LLP = new TH1F("h7_LLP","",nbin,xmin,xmax);
  h7_LLP->Add(g7_LLP, h7_LLP, 1,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(ColorBlue);
 h1_LLP->SetLineColorAlpha(ColorBlue,1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);
 h1_LLP->SetTickLength(0.03, "YZ");
 h1_LLP->SetTickLength(0.03,"X");
 h1_LLP->SetLabelOffset(0.015,"X");
 h1_LLP->SetLabelOffset(0.007,"Y");
 h1_LLP->SetLabelSize(0.045, "XYZ");
 h1_LLP->SetLabelFont(42, "XYZ"); 
 h1_LLP->SetTitleSize(0.045, "XYZ"); 
 h1_LLP->SetTitleFont(42, "XYZ");
 h1_LLP->SetTitleOffset(1.2,"X"); 
 h1_LLP->SetTitleOffset(1.3,"Y");
 h1_LLP->GetXaxis()->SetTitle(xtitle);
 h1_LLP->GetXaxis()->SetTitleColor(1);
 h1_LLP->GetYaxis()->SetTitle(ytitle);
 h1_LLP->GetYaxis()->SetTitleColor(1);
 h1_LLP->SetNdivisions(509,"XYZ");
//  h1_LLP->SetMinimum(hmin); 
//  h1_LLP->SetMaximum(hmax);
 if (logy)
  {
    h1_LLP->SetMinimum(1); 
    h1_LLP->SetMaximum(h1_LLP->GetMaximum()*100); 
  }
else 
  {
    h1_LLP->SetMinimum(1); 
    h1_LLP->SetMaximum(100); 
  }

// Int_t ColorBlue = color1.GetNumber();
// Int_t ColorLightBlue = color10.GetNumber();
// Int_t ColorDarkGrey = color9.GetNumber();
// Int_t ColorGrey = color4.GetNumber();
// Int_t ColorNeutral = color8.GetNumber();
// Int_t ColorDarkOrange = color7.GetNumber();
// Int_t ColorOrange = color2.GetNumber();
// Int_t ColorRed = color3.GetNumber();
// Int_t ColorBrown = color6.GetNumber();
// Int_t ColorDarkPurple = color5.GetNumber();

 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(ColorLightBlue);
 h2_LLP->SetLineColorAlpha(ColorLightBlue,1);
 h2_LLP->SetLineStyle(1);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(ColorDarkGrey);
 h3_LLP->SetLineColorAlpha(ColorDarkGrey,1);
 h3_LLP->SetLineStyle(1);
 h3_LLP->SetLineWidth(2);

 h4_LLP->Draw("HEsame"); 
 h4_LLP->SetLineColor(ColorGrey);
 h4_LLP->SetLineColorAlpha(ColorGrey,1);
 h4_LLP->SetLineStyle(1);
 h4_LLP->SetLineWidth(2);

 h5_LLP->Draw("HEsame"); 
 h5_LLP->SetLineColor(ColorNeutral);
 h5_LLP->SetLineColorAlpha(ColorNeutral,1);
 h5_LLP->SetLineStyle(1);
 h5_LLP->SetLineWidth(2);

  h6_LLP->Draw("HEsame"); 
 h6_LLP->SetLineColor(ColorDarkOrange);
 h6_LLP->SetLineColorAlpha(ColorDarkOrange,1);
 h6_LLP->SetLineStyle(1);
 h6_LLP->SetLineWidth(2);

  h7_LLP->Draw("HEsame"); 
 h7_LLP->SetLineColor(ColorOrange);
 h7_LLP->SetLineColorAlpha(ColorOrange,1);
 h7_LLP->SetLineStyle(1);
 h7_LLP->SetLineWidth(2);


  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  leg = new TLegend(0.51,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
  // leg->AddEntry(htotData, " #mu#mu data","PE1");

  leg->SetHeader(" m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSmuon+" ("+MNeu+") GeV ","L");
  leg->AddEntry(h1_LLP," c#tau = 0.1 cm","L");
  leg->AddEntry(h2_LLP,"c#tau = 0.3 cm","L");
  leg->AddEntry(h3_LLP,"c#tau = 1.0 cm","L");
  leg->AddEntry(h4_LLP,"c#tau = 3.0 cm","L");
  leg->AddEntry(h5_LLP,"c#tau = 10.0 cm","L");
  leg->AddEntry(h6_LLP,"c#tau = 30.0 cm","L");
  leg->AddEntry(h7_LLP,"c#tau = 100.0 cm","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();


  leg = new TLegend(0.2,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  
  leg = new TLegend(0.2,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderA);
 leg->Draw();
   leg = new TLegend(0.2,0.69,0.35,0.74);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderAbis);
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 20-80 GeV");
  leg->Draw();

  leg = new TLegend(0.60,0.30,0.65,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("A");
  leg->Draw();
// *****************************************************************************

 pB->cd();



 f1_LLP->cd();

 g1_LLP = (TH1F*)gROOT->FindObject(FILE[0]+htitleB);
 
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

 f2_LLP->cd();

 g2_LLP = (TH1F*)gROOT->FindObject(FILE[1]+htitleB);
 
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

 f3_LLP->cd();

 g3_LLP = (TH1F*)gROOT->FindObject(FILE[2]+htitleB);
 
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

  f4_LLP->cd();

   g4_LLP = (TH1F*)gROOT->FindObject(FILE[3]+htitleB);
  
   h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
  h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

    f5_LLP->cd();

   g5_LLP = (TH1F*)gROOT->FindObject(FILE[4]+htitleB);
  
   h5_LLP = new TH1F("h5_LLP","",nbin,xmin,xmax);
  h5_LLP->Add(g5_LLP, h5_LLP, 1,0);

    f6_LLP->cd();

   g6_LLP = (TH1F*)gROOT->FindObject(FILE[5]+htitleB);
  
   h6_LLP = new TH1F("h6_LLP","",nbin,xmin,xmax);
  h6_LLP->Add(g6_LLP, h6_LLP, 1,0);

    f7_LLP->cd();

   g7_LLP = (TH1F*)gROOT->FindObject(FILE[6]+htitleB);
  
   h7_LLP = new TH1F("h7_LLP","",nbin,xmin,xmax);
  h7_LLP->Add(g7_LLP, h7_LLP, 1,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(ColorBlue);
 h1_LLP->SetLineColorAlpha(ColorBlue,1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);
 h1_LLP->SetTickLength(0.03, "YZ");
 h1_LLP->SetTickLength(0.03,"X");
 h1_LLP->SetLabelOffset(0.015,"X");
 h1_LLP->SetLabelOffset(0.007,"Y");
 h1_LLP->SetLabelSize(0.045, "XYZ");
 h1_LLP->SetLabelFont(42, "XYZ"); 
 h1_LLP->SetTitleSize(0.045, "XYZ"); 
 h1_LLP->SetTitleFont(42, "XYZ");
 h1_LLP->SetTitleOffset(1.2,"X"); 
 h1_LLP->SetTitleOffset(1.3,"Y");
 h1_LLP->GetXaxis()->SetTitle(xtitle);
 h1_LLP->GetXaxis()->SetTitleColor(1);
 h1_LLP->GetYaxis()->SetTitle(ytitle);
 h1_LLP->GetYaxis()->SetTitleColor(1);
 h1_LLP->SetNdivisions(509,"XYZ");
 if (logy)
  {
    h1_LLP->SetMinimum(1); 
    h1_LLP->SetMaximum(h1_LLP->GetMaximum()*100); 
  }
else 
  {
    h1_LLP->SetMinimum(1); 
    h1_LLP->SetMaximum(100); 
  }




 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(ColorLightBlue);
 h2_LLP->SetLineColorAlpha(ColorLightBlue,1);
 h2_LLP->SetLineStyle(1);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(ColorDarkGrey);
 h3_LLP->SetLineColorAlpha(ColorDarkGrey,1);
 h3_LLP->SetLineStyle(1);
 h3_LLP->SetLineWidth(2);

 h4_LLP->Draw("HEsame"); 
 h4_LLP->SetLineColor(ColorGrey);
 h4_LLP->SetLineColorAlpha(ColorGrey,1);
 h4_LLP->SetLineStyle(1);
 h4_LLP->SetLineWidth(2);

 h5_LLP->Draw("HEsame"); 
 h5_LLP->SetLineColor(ColorNeutral);
 h5_LLP->SetLineColorAlpha(ColorNeutral,1);
 h5_LLP->SetLineStyle(1);
 h5_LLP->SetLineWidth(2);

  h6_LLP->Draw("HEsame"); 
 h6_LLP->SetLineColor(ColorDarkOrange);
 h6_LLP->SetLineColorAlpha(ColorDarkOrange,1);
 h6_LLP->SetLineStyle(1);
 h6_LLP->SetLineWidth(2);

  h7_LLP->Draw("HEsame"); 
 h7_LLP->SetLineColor(ColorOrange);
 h7_LLP->SetLineColorAlpha(ColorOrange,1);
 h7_LLP->SetLineStyle(1);
 h7_LLP->SetLineWidth(2);


  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  leg = new TLegend(0.51,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);

  leg->SetHeader(" m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSmuon+" ("+MNeu+") GeV ","L");
  leg->AddEntry(h1_LLP," c#tau = 0.1 cm","L");
  leg->AddEntry(h2_LLP,"c#tau = 0.3 cm","L");
  leg->AddEntry(h3_LLP,"c#tau = 1.0 cm","L");
  leg->AddEntry(h4_LLP,"c#tau = 3.0 cm","L");
  leg->AddEntry(h5_LLP,"c#tau = 10.0 cm","L");
  leg->AddEntry(h6_LLP,"c#tau = 30.0 cm","L");
  leg->AddEntry(h7_LLP,"c#tau = 100.0 cm","L");
// leg->AddEntry(hsolve," Prediction","FE4");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();


  leg = new TLegend(0.2,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  
  leg = new TLegend(0.2,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderB);
  leg->Draw();

   leg = new TLegend(0.2,0.69,0.35,0.74);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderBbis);
  leg->Draw();

  leg = new TLegend(0.60,0.30,0.65,0.35);
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

 f1_LLP->cd();

 g1_LLP = (TH1F*)gROOT->FindObject(FILE[0]+htitleC);
 
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

 f2_LLP->cd();

 g2_LLP = (TH1F*)gROOT->FindObject(FILE[1]+htitleC);
 
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

 f3_LLP->cd();

 g3_LLP = (TH1F*)gROOT->FindObject(FILE[2]+htitleC);
 
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

  f4_LLP->cd();

   g4_LLP = (TH1F*)gROOT->FindObject(FILE[3]+htitleC);
  
   h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
  h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

    f5_LLP->cd();

   g5_LLP = (TH1F*)gROOT->FindObject(FILE[4]+htitleC);
  
   h5_LLP = new TH1F("h5_LLP","",nbin,xmin,xmax);
  h5_LLP->Add(g5_LLP, h5_LLP, 1,0);

    f6_LLP->cd();

   g6_LLP = (TH1F*)gROOT->FindObject(FILE[5]+htitleC);
  
   h6_LLP = new TH1F("h6_LLP","",nbin,xmin,xmax);
  h6_LLP->Add(g6_LLP, h6_LLP, 1,0);

    f7_LLP->cd();

   g7_LLP = (TH1F*)gROOT->FindObject(FILE[6]+htitleC);
  
   h7_LLP = new TH1F("h7_LLP","",nbin,xmin,xmax);
  h7_LLP->Add(g7_LLP, h7_LLP, 1,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(ColorBlue);
 h1_LLP->SetLineColorAlpha(ColorBlue,1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);
 h1_LLP->SetTickLength(0.03, "YZ");
 h1_LLP->SetTickLength(0.03,"X");
 h1_LLP->SetLabelOffset(0.015,"X");
 h1_LLP->SetLabelOffset(0.007,"Y");
 h1_LLP->SetLabelSize(0.045, "XYZ");
 h1_LLP->SetLabelFont(42, "XYZ"); 
 h1_LLP->SetTitleSize(0.045, "XYZ"); 
 h1_LLP->SetTitleFont(42, "XYZ");
 h1_LLP->SetTitleOffset(1.2,"X"); 
 h1_LLP->SetTitleOffset(1.3,"Y");
 h1_LLP->GetXaxis()->SetTitle(xtitle);
 h1_LLP->GetXaxis()->SetTitleColor(1);
 h1_LLP->GetYaxis()->SetTitle(ytitle);
 h1_LLP->GetYaxis()->SetTitleColor(1);
 h1_LLP->SetNdivisions(509,"XYZ");
 if (logy)
  {
    h1_LLP->SetMinimum(1); 
    h1_LLP->SetMaximum(h1_LLP->GetMaximum()*100); 
  }
else 
  {
    h1_LLP->SetMinimum(1); 
    h1_LLP->SetMaximum(100); 
  }



 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(ColorLightBlue);
 h2_LLP->SetLineColorAlpha(ColorLightBlue,1);
 h2_LLP->SetLineStyle(1);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(ColorDarkGrey);
 h3_LLP->SetLineColorAlpha(ColorDarkGrey,1);
 h3_LLP->SetLineStyle(1);
 h3_LLP->SetLineWidth(2);

 h4_LLP->Draw("HEsame"); 
 h4_LLP->SetLineColor(ColorGrey);
 h4_LLP->SetLineColorAlpha(ColorGrey,1);
 h4_LLP->SetLineStyle(1);
 h4_LLP->SetLineWidth(2);

 h5_LLP->Draw("HEsame"); 
 h5_LLP->SetLineColor(ColorNeutral);
 h5_LLP->SetLineColorAlpha(ColorNeutral,1);
 h5_LLP->SetLineStyle(1);
 h5_LLP->SetLineWidth(2);

  h6_LLP->Draw("HEsame"); 
 h6_LLP->SetLineColor(ColorDarkOrange);
 h6_LLP->SetLineColorAlpha(ColorDarkOrange,1);
 h6_LLP->SetLineStyle(1);
 h6_LLP->SetLineWidth(2);

  h7_LLP->Draw("HEsame"); 
 h7_LLP->SetLineColor(ColorOrange);
 h7_LLP->SetLineColorAlpha(ColorOrange,1);
 h7_LLP->SetLineStyle(1);
 h7_LLP->SetLineWidth(2);


  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderCMS);
  leg->Draw();

  leg = new TLegend(0.51,0.50,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);

  leg->SetHeader(" m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSmuon+" ("+MNeu+") GeV ","L");
  leg->AddEntry(h1_LLP," c#tau = 0.1 cm","L");
  leg->AddEntry(h2_LLP,"c#tau = 0.3 cm","L");
  leg->AddEntry(h3_LLP,"c#tau = 1.0 cm","L");
  leg->AddEntry(h4_LLP,"c#tau = 3.0 cm","L");
  leg->AddEntry(h5_LLP,"c#tau = 10.0 cm","L");
  leg->AddEntry(h6_LLP,"c#tau = 30.0 cm","L");
  leg->AddEntry(h7_LLP,"c#tau = 100.0 cm","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  
  leg = new TLegend(0.2,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderC);
  leg->Draw();

     leg = new TLegend(0.2,0.69,0.35,0.74);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderCbis);
//   if ( nvtx == 2 ) leg->SetHeader(" #geq1 hem. p_{T} 20-80 GeV");
  leg->Draw();
  leg = new TLegend(0.60,0.30,0.65,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("C");
  leg->Draw();
// *****************************************************************************

 pD->cd();
//  hmax = hmaxBD; 

 f1_LLP->cd();

 g1_LLP = (TH1F*)gROOT->FindObject(FILE[0]+htitleD);
 
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

 f2_LLP->cd();

 g2_LLP = (TH1F*)gROOT->FindObject(FILE[1]+htitleD);
 
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

 f3_LLP->cd();

 g3_LLP = (TH1F*)gROOT->FindObject(FILE[2]+htitleD);
 
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

  f4_LLP->cd();

   g4_LLP = (TH1F*)gROOT->FindObject(FILE[3]+htitleD);
  
   h4_LLP = new TH1F("h4_LLP","",nbin,xmin,xmax);
  h4_LLP->Add(g4_LLP, h4_LLP, 1,0);

    f5_LLP->cd();

   g5_LLP = (TH1F*)gROOT->FindObject(FILE[4]+htitleD);
  
   h5_LLP = new TH1F("h5_LLP","",nbin,xmin,xmax);
  h5_LLP->Add(g5_LLP, h5_LLP, 1,0);

    f6_LLP->cd();

   g6_LLP = (TH1F*)gROOT->FindObject(FILE[5]+htitleD);
  
   h6_LLP = new TH1F("h6_LLP","",nbin,xmin,xmax);
  h6_LLP->Add(g6_LLP, h6_LLP, 1,0);

    f7_LLP->cd();

   g7_LLP = (TH1F*)gROOT->FindObject(FILE[6]+htitleD);
  
   h7_LLP = new TH1F("h7_LLP","",nbin,xmin,xmax);
  h7_LLP->Add(g7_LLP, h7_LLP, 1,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(ColorBlue);
 h1_LLP->SetLineColorAlpha(ColorBlue,1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);
 h1_LLP->SetTickLength(0.03, "YZ");
 h1_LLP->SetTickLength(0.03,"X");
 h1_LLP->SetLabelOffset(0.015,"X");
 h1_LLP->SetLabelOffset(0.007,"Y");
 h1_LLP->SetLabelSize(0.045, "XYZ");
 h1_LLP->SetLabelFont(42, "XYZ"); 
 h1_LLP->SetTitleSize(0.045, "XYZ"); 
 h1_LLP->SetTitleFont(42, "XYZ");
 h1_LLP->SetTitleOffset(1.2,"X"); 
 h1_LLP->SetTitleOffset(1.3,"Y");
 h1_LLP->GetXaxis()->SetTitle(xtitle);
 h1_LLP->GetXaxis()->SetTitleColor(1);
 h1_LLP->GetYaxis()->SetTitle(ytitle);
 h1_LLP->GetYaxis()->SetTitleColor(1);
 h1_LLP->SetNdivisions(509,"XYZ");
 if (logy)
  {
    h1_LLP->SetMinimum(1); 
    h1_LLP->SetMaximum(h1_LLP->GetMaximum()*100); 
  }
else 
  {
    h1_LLP->SetMinimum(1); 
    h1_LLP->SetMaximum(100); 
  }
 
 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(ColorLightBlue);
 h2_LLP->SetLineColorAlpha(ColorLightBlue,1);
 h2_LLP->SetLineStyle(1);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(ColorDarkGrey);
 h3_LLP->SetLineColorAlpha(ColorDarkGrey,1);
 h3_LLP->SetLineStyle(1);
 h3_LLP->SetLineWidth(2);

 h4_LLP->Draw("HEsame"); 
 h4_LLP->SetLineColor(ColorGrey);
 h4_LLP->SetLineColorAlpha(ColorGrey,1);
 h4_LLP->SetLineStyle(1);
 h4_LLP->SetLineWidth(2);

 h5_LLP->Draw("HEsame"); 
 h5_LLP->SetLineColor(ColorNeutral);
 h5_LLP->SetLineColorAlpha(ColorNeutral,1);
 h5_LLP->SetLineStyle(1);
 h5_LLP->SetLineWidth(2);

  h6_LLP->Draw("HEsame"); 
 h6_LLP->SetLineColor(ColorDarkOrange);
 h6_LLP->SetLineColorAlpha(ColorDarkOrange,1);
 h6_LLP->SetLineStyle(1);
 h6_LLP->SetLineWidth(2);

  h7_LLP->Draw("HEsame"); 
 h7_LLP->SetLineColor(ColorOrange);
 h7_LLP->SetLineColorAlpha(ColorOrange,1);
 h7_LLP->SetLineStyle(1);
 h7_LLP->SetLineWidth(2);



  leg = new TLegend(0.17,0.94,0.50,0.99);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderCMS);
  leg->Draw();


  leg = new TLegend(0.51,0.5,0.85,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);

  leg->SetHeader(" m_{#tilde{#mu} (#tilde{#chi}^{0})} = "+MSmuon+" ("+MNeu+") GeV ","L");
  leg->AddEntry(h1_LLP," c#tau = 0.1 cm","L");
  leg->AddEntry(h2_LLP,"c#tau = 0.3 cm","L");
  leg->AddEntry(h3_LLP,"c#tau = 1.0 cm","L");
  leg->AddEntry(h4_LLP,"c#tau = 3.0 cm","L");
  leg->AddEntry(h5_LLP,"c#tau = 10.0 cm","L");
  leg->AddEntry(h6_LLP,"c#tau = 30.0 cm","L");
  leg->AddEntry(h7_LLP,"c#tau = 100.0 cm","L");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.2,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  
  leg = new TLegend(0.2,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderD);
  leg->Draw();

   leg = new TLegend(0.2,0.69,0.35,0.74);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetHeader(HeaderDbis);
leg->Draw();
  
     leg = new TLegend(0.60,0.30,0.65,0.35);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.1);
  leg->SetHeader("D");
  leg->Draw();


//************************************µ//
  return c1;
}