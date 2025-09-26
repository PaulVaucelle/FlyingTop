#include <iostream>
#include <TROOT.h>
#include "TVectorD.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TEfficiency.h" 
#include "TMath.h" 
#include "../outputroot/Func.h"

void plot()
{
// *****************************************************************************
// gROOT->LoadMacro("../outputroot/tdrstyle.C");
int stati=0;
bool fit= 0;
bool logy=0;
// TCanvas *c1 = new TCanvas("c1", "plots",200,0,900,500);
TCanvas *c1 = new TCanvas("c1", "plots",200,0,1500,1000);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);


TPad* pad1 = new TPad("pad1","This is pad1",0.01,0.51,0.49,1,21);
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.1);
   pad1->SetBottomMargin(0.15);
   pad1->SetRightMargin(0.06);
   pad1->SetLeftMargin(0.12);


TPad* pad2 = new TPad("pad2","This is pad2",0.01,0.01,0.49,0.49,21);
pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(0);
   pad2->SetTopMargin(0.1);
   pad2->SetBottomMargin(0.20);
   pad2->SetRightMargin(0.05);
   pad2->SetLeftMargin(0.15);


TPad* pad3 = new TPad("pad3","This is pad3",0.53,0.51,0.99,1,21);
pad3->SetFillColor(0);
pad3->SetBorderMode(0);
pad3->SetFrameFillColor(10);
pad3->Draw();
pad3->SetLogy(0);
   pad3->SetTopMargin(0.1);
   pad3->SetBottomMargin(0.15);
   pad3->SetRightMargin(0.05);
   pad3->SetLeftMargin(0.10);

   TPad* pad4 = new TPad("pad4","This is pad4",0.53,0.01,0.99,0.49,21);
pad4->SetFillColor(0);
pad4->SetBorderMode(0);
pad4->SetFrameFillColor(10);
pad4->Draw();
pad4->SetLogy(0);
   pad4->SetTopMargin(0.1);
   pad4->SetBottomMargin(0.20);
   pad4->SetRightMargin(0.05);
   pad4->SetLeftMargin(0.10);


gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleFont(62);
gStyle->SetTitleColor(1);
gStyle->SetTitleTextColor(1);
gStyle->SetTitleFillColor(10);
gStyle->SetTitleFontSize(0.06);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(0.09);
gStyle->SetOptStat(stati);
gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
if (fit) {
gStyle->SetStatW(0.3);
gStyle->SetStatH(0.1);
gStyle->SetOptFit(111);
} else {
gStyle->SetStatW(0.3);
gStyle->SetStatH(0.2);
gStyle->SetOptFit(0);
}

// *****************************************************************************

// setTDRStyle();
  pad1->cd();


// htitle0 = "hData_dR_GenGen";  
// htitle1 = "hData_dR_RecoReco"; 
// htitle2 = "hData_dR_GenReco"; 
TString Prod = "Signal_2018";
const int nMC = 3;
  TString  Datasets[nMC] = {"RPV_2018_smu250_neu180_ctau001","RPV_2018_smu400_neu300_ctau010","RPV_2018_smu500_neu450_ctau100"};


//      TString SignalSet001[34]={"RPV_"+Year+"_smu200_neu180_ctau001","RPV_"+Year+"_smu250_neu180_ctau001","RPV_"+Year+"_smu250_neu200_ctau001",
// "RPV_"+Year+"_smu250_neu230_ctau001","RPV_"+Year+"_smu300_neu180_ctau001","RPV_"+Year+"_smu300_neu200_ctau001","RPV_"+Year+"_smu300_neu280_ctau001","RPV_"+Year+"_smu300_neu250_ctau001",
// "RPV_"+Year+"_smu350_neu180_ctau001","RPV_"+Year+"_smu350_neu200_ctau001","RPV_"+Year+"_smu350_neu250_ctau001","RPV_"+Year+"_smu350_neu300_ctau001","RPV_"+Year+"_smu350_neu330_ctau001",
// "RPV_"+Year+"_smu400_neu180_ctau001","RPV_"+Year+"_smu400_neu200_ctau001","RPV_"+Year+"_smu400_neu250_ctau001","RPV_"+Year+"_smu400_neu300_ctau001","RPV_"+Year+"_smu400_neu350_ctau001",
// "RPV_"+Year+"_smu400_neu380_ctau001","RPV_"+Year+"_smu450_neu180_ctau001","RPV_"+Year+"_smu450_neu200_ctau001","RPV_"+Year+"_smu450_neu250_ctau001","RPV_"+Year+"_smu450_neu300_ctau001",
// "RPV_"+Year+"_smu450_neu350_ctau001","RPV_"+Year+"_smu450_neu400_ctau001","RPV_"+Year+"_smu450_neu430_ctau001","RPV_"+Year+"_smu500_neu180_ctau001","RPV_"+Year+"_smu500_neu200_ctau001",
// "RPV_"+Year+"_smu500_neu250_ctau001","RPV_"+Year+"_smu500_neu300_ctau001","RPV_"+Year+"_smu500_neu350_ctau001","RPV_"+Year+"_smu500_neu400_ctau001","RPV_"+Year+"_smu500_neu450_ctau001",
// "RPV_"+Year+"_smu500_neu480_ctau001"};


//       TString SignalSet001[34]={"RPV_"+Year+"_smu200_neu180","RPV_"+Year+"_smu250_neu180","RPV_"+Year+"_smu250_neu200",
// "RPV_"+Year+"_smu250_neu230","RPV_"+Year+"_smu300_neu180","RPV_"+Year+"_smu300_neu200","RPV_"+Year+"_smu300_neu250","RPV_"+Year+"_smu300_neu280",
// "RPV_"+Year+"_smu350_neu180","RPV_"+Year+"_smu350_neu200","RPV_"+Year+"_smu350_neu250","RPV_"+Year+"_smu350_neu300","RPV_"+Year+"_smu350_neu330",
// "RPV_"+Year+"_smu400_neu180","RPV_"+Year+"_smu400_neu200","RPV_"+Year+"_smu400_neu250","RPV_"+Year+"_smu400_neu300","RPV_"+Year+"_smu400_neu350",
// "RPV_"+Year+"_smu400_neu380","RPV_"+Year+"_smu450_neu180","RPV_"+Year+"_smu450_neu200","RPV_"+Year+"_smu450_neu250","RPV_"+Year+"_smu450_neu300",
// "RPV_"+Year+"_smu450_neu350","RPV_"+Year+"_smu450_neu400","RPV_"+Year+"_smu450_neu430","RPV_"+Year+"_smu500_neu180","RPV_"+Year+"_smu500_neu200",
// "RPV_"+Year+"_smu500_neu250","RPV_"+Year+"_smu500_neu300","RPV_"+Year+"_smu500_neu350","RPV_"+Year+"_smu500_neu400","RPV_"+Year+"_smu500_neu450",
// "RPV_"+Year+"_smu500_neu480"};


  TString  Name[nMC] = {"M_{#tilde{#mu}} = 250 GeV, M_{#tilde{#chi}} = 180 GeV, c#tau = 0.1 cm ","M_{#tilde{#mu}} = 400 GeV, M_{#tilde{#chi}} = 300 GeV, c#tau = 1.0 cm","M_{#tilde{#mu}} = 500 GeV, M_{#tilde{#chi}} = 450 GeV, c#tau = 10.0 cm"};
    TFile *fX[nMC];
    for (int i = 0; i < nMC; i++) 
      {
        TString NAME = "../"+Prod+"/histoGen_"+Datasets[i]+".root";
        fX[i] = new TFile("../"+Prod+"/histoGen_"+Datasets[i]+".root");
        std::cout<<"FileName : "<<NAME<<std::endl;
      }



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

  TString htitle0 = "dEtaNeuNeu_";
  TString xtitle0 = "#Delta #eta_{#chi-#chi}"; 
  TString ytitle0 = "a.u"; 
  TH1F* histograms0[nMC];
  // int Colors[nMC] = {4,1,2,3,33,5,6,7,8,9,30,38,41,43,46,49};
  int Colors0[nMC]={ColorBlue,ColorRed,ColorNeutral};

for (int i = 0; i < nMC; i++) {
    // Colors0[i]= i+1;
    fX[i]->cd();
    TString histogramTitle = Datasets[i]+"_"+htitle0; // Assuming htitle is an array of strings
    std::cout<<"histogramTitle : "<<histogramTitle<<std::endl;
    histograms0[i] = (TH1F*)gROOT->FindObject(histogramTitle);
    // histograms0[i]->GetXaxis()->SetRangeUser(1,5);
}

// /////////
  TLegend* leg0 = new TLegend(0.5,0.6,0.8,0.85);
  leg0->SetHeader("MC #tilde{#mu}^{+}#tilde{#mu}^{-}(j)  with  #tilde{#mu} #rightarrow #mu#tilde{#chi}^{0},  #tilde{#chi}^{0} #rightarrow tds");
  leg0->SetBorderSize(0);
  leg0->SetFillStyle(0);
  leg0->SetFillColor(kWhite);
  leg0->SetTextFont(42);
  leg0->SetTextSize(0.03);

   histograms0[0]->SetLineColor(Colors0[0]);
   histograms0[0]->SetLineWidth(2);
   histograms0[0]->SetLineStyle(1);
   histograms0[0]->SetMarkerStyle(76);
   histograms0[0]->SetMarkerSize(1);
   histograms0[0]->SetMarkerColor(Colors0[0]);
   histograms0[0]->Scale(1./histograms0[0]->Integral());
   histograms0[0]->Draw("HIST");


  histograms0[0]->GetXaxis()->SetTitleSize(0.06);
  histograms0[0]->GetYaxis()->SetTitleSize(0.06);
  histograms0[0]->GetYaxis()->SetTitleOffset(1.5);
  histograms0[0]->SetTickLength(0.03, "XYZ");
  histograms0[0]->SetLabelOffset(0.002,"X");//0.007
  histograms0[0]->SetLabelOffset(0.007,"Y");
  histograms0[0]->SetLabelSize(0.042, "XYZ");
  histograms0[0]->SetLabelFont(42, "XYZ"); 
  histograms0[0]->SetTitleFont(42, "XYZ");
  histograms0[0]->SetTitleSize(0.05, "XYZ"); 
  histograms0[0]->SetTitleOffset(1.,"Y");
  histograms0[0]->SetTitleOffset(0.85,"X");
  histograms0[0]->SetNdivisions(509,"XYZ");

  leg0->AddEntry(histograms0[0],Name[0],"L");

  for (unsigned int k = 1; k < nMC; k++)
  {
    StyleAndDraw(histograms0[k], Colors0[k], 2, 1, Datasets[k]);
    leg0->AddEntry(histograms0[k],Name[k],"L");
  }

  histograms0[0]->SetMaximum(0.15);
  histograms0[0]->SetTitle("");
  histograms0[0]->GetXaxis()->SetTitle(xtitle0);
  // histograms0[0]->GetXaxis()->SetTitleSize(0.06);
  histograms0[0]->GetYaxis()->SetTitle(ytitle0);
  // histograms0[0]->GetYaxis()->SetTitleSize(0.06);
  // histograms0[0]->GetYaxis()->SetTitleOffset(1.5);
 leg0->Draw();


 //------Start of Copy Paste
TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = true;
TString extraText   = "Simulation";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "";//137 fb^{-1}
TString lumi_sqrtS = "2018";
TString lumiText = lumi_13TeV+lumi_sqrtS;
  float H = pad1->GetWh();
  float W = pad1->GetWw();
  float l = pad1->GetLeftMargin();
  float t = pad1->GetTopMargin();
  float r = pad1->GetRightMargin();
  float b = pad1->GetBottomMargin();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);

float posX_=0;
  float posY_=0;
  int iPosX = 3;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  posY_ = 1-t - relPosY*(1-t-b);
  	  if( writeExtraText ) 
	    {
         posX_ =   l +  relPosX*(1-l-r);
         posY_ =   1-t+lumiTextOffset*t;
        int alignY_=3;
         int alignX_=2;
         if( iPosX/10==0 ) alignX_=1;
         if( iPosX==0    ) alignX_=1;
         if( iPosX==0    ) alignY_=1;
         if( iPosX/10==1 ) alignX_=1;
         if( iPosX/10==2 ) alignX_=2;
         if( iPosX/10==3 ) alignX_=3;
         //if( iPosX == 0  ) relPosX = 0.12;
         int align_ = 10*alignX_ + alignY_;

      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(11);
      latex.DrawLatex(posX_+0.08, posY_, extraText);
	    }
 c1->Update();
// *****************************************************************************

  pad2->cd();

  /////////
  TString htitle1 = "RecoRecoAxis_dR_";
  TString xtitle1 = "#Delta R_{H-H}"; 
  TString ytitle1 = "a.u"; 
  TH1F* histograms1[nMC];
  // int Colors[nMC] = {4,1,2,3,33,5,6,7,8,9,30,38,41,43,46,49};
  int Colors1[nMC]={ColorBlue,ColorRed,ColorNeutral};
  
for (int i = 0; i < nMC; i++) {
    
    fX[i]->cd();
    TString histogramTitle = Datasets[i]+"_"+htitle1; // Assuming htitle is an array of strings
    histograms1[i] = (TH1F*)gROOT->FindObject(histogramTitle);
    std::cout<<"histogramTitle : "<<histogramTitle<<std::endl;
    // histograms1[i]->GetXaxis()->SetRangeUser(1,5);
}

/////////
  TLegend* leg1 = new TLegend(0.2,0.6,0.45,0.85);
  leg1->SetHeader("MC #tilde{#mu}^{+}#tilde{#mu}^{-}(j)  with  #tilde{#mu} #rightarrow #mu#tilde{#chi}^{0},  #tilde{#chi}^{0} #rightarrow tds");
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(kWhite);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.03);


   histograms1[0]->SetLineColor(Colors0[0]);
   histograms1[0]->SetLineWidth(2);
   histograms1[0]->SetLineStyle(1);
   histograms1[0]->SetMarkerStyle(76);
   histograms1[0]->SetMarkerSize(1);
   histograms1[0]->SetMarkerColor(Colors0[0]);
   histograms1[0]->Scale(1./histograms1[0]->Integral());
   histograms1[0]->Draw("HIST");


  histograms1[0]->GetXaxis()->SetTitleSize(0.06);
  histograms1[0]->GetYaxis()->SetTitleSize(0.06);
  histograms1[0]->GetYaxis()->SetTitleOffset(1.5);
  histograms1[0]->SetTickLength(0.03, "XYZ");
  histograms1[0]->SetLabelOffset(0.002,"X");//0.007
  histograms1[0]->SetLabelOffset(0.007,"Y");
  histograms1[0]->SetLabelSize(0.042, "XYZ");
  histograms1[0]->SetLabelFont(42, "XYZ"); 
  histograms1[0]->SetTitleFont(42, "XYZ");
  histograms1[0]->SetTitleSize(0.05, "XYZ"); 
  histograms1[0]->SetTitleOffset(1.1,"Y");
  histograms1[0]->SetTitleOffset(0.85,"X");
  histograms1[0]->SetNdivisions(505,"XYZ");


  leg1->AddEntry(histograms1[0],Name[0],"L");

  for (unsigned int k = 1; k < nMC; k++)
  {
    StyleAndDraw(histograms1[k], Colors1[k], 2, 1, Datasets[k]);
    leg1->AddEntry(histograms1[k],Name[k],"L");
  }

  histograms1[0]->SetMaximum(0.25);
  histograms1[0]->SetTitle("");
  histograms1[0]->GetXaxis()->SetTitle(xtitle1);
  // histograms1[0]->GetXaxis()->SetTitleSize(0.06);
  histograms1[0]->GetYaxis()->SetTitle(ytitle1);
  // histograms1[0]->GetYaxis()->SetTitleSize(0.06);
  // histograms1[0]->GetYaxis()->SetTitleOffset(1.5);
 leg1->Draw();


 //------Start of Copy Paste


// text sizes and text offsets with respect to the top frame
// in unit of the top margin size

  H = pad2->GetWh();
  W = pad2->GetWw();
  l = pad2->GetLeftMargin();
  t = pad2->GetTopMargin();
  r = pad2->GetRightMargin();
  b = pad2->GetBottomMargin();

  TLatex latex1;
  latex1.SetNDC();
  latex1.SetTextAngle(0);
  latex1.SetTextColor(kBlack);    

  extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex1.SetTextFont(42);
  latex1.SetTextAlign(31); 
  latex1.SetTextSize(lumiTextSize*t);    
  latex1.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

      latex1.SetTextFont(cmsTextFont);
      latex1.SetTextAlign(11); 
      latex1.SetTextSize(cmsTextSize*t);    
      latex1.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);

  posX_=0;
   posY_=0;
  iPosX = 3;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
   posY_ = 1-t - relPosY*(1-t-b);
  	  if( writeExtraText ) 
	    {
         posX_ =   l +  relPosX*(1-l-r);
         posY_ =   1-t+lumiTextOffset*t;
        int alignY_=3;
         int alignX_=2;
         if( iPosX/10==0 ) alignX_=1;
         if( iPosX==0    ) alignX_=1;
         if( iPosX==0    ) alignY_=1;
         if( iPosX/10==1 ) alignX_=1;
         if( iPosX/10==2 ) alignX_=2;
         if( iPosX/10==3 ) alignX_=3;
         //if( iPosX == 0  ) relPosX = 0.12;
         int align_ = 10*alignX_ + alignY_;
      latex1.SetTextFont(extraTextFont);
      latex1.SetTextSize(extraTextSize*t);
      latex1.SetTextAlign(11);
      latex1.DrawLatex(posX_+0.08, posY_, extraText);
	    }
 c1->Update();
// // *****************************************************************************

  pad3->cd();

   /////////
  TString htitle2 = "dPhiNeuNeu_";
  TString xtitle2 = "#Delta #Phi_{#chi-#chi}"; 
  TString ytitle2 = "a.u"; 
  TH1F* histograms2[nMC];
  // int Colors[nMC] = {4,1,2,3,33,5,6,7,8,9,30,38,41,43,46,49};
  int Colors2[nMC]={ColorBlue,ColorRed,ColorNeutral};

for (int i = 0; i < nMC; i++) {
    // Colors2[i]= i+1;
    fX[i]->cd();
    TString histogramTitle = Datasets[i]+"_"+htitle2; // Assuming htitle is an array of strings
    histograms2[i] = (TH1F*)gROOT->FindObject(histogramTitle);
    // histograms2[i]->GetXaxis()->SetRangeUser(0,5);
}

  TLegend* leg2 = new TLegend(0.15,0.6,0.45,0.85);
  leg2->SetHeader("MC #tilde{#mu}^{+}#tilde{#mu}^{-}(j)  with  #tilde{#mu} #rightarrow #mu#tilde{#chi}^{0},  #tilde{#chi}^{0} #rightarrow tds");
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetFillColor(kWhite);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.03);


   histograms2[0]->SetLineColor(Colors0[0]);
   histograms2[0]->SetLineWidth(2);
   histograms2[0]->SetLineStyle(1);
   histograms2[0]->SetMarkerStyle(76);
   histograms2[0]->SetMarkerSize(1);
   histograms2[0]->SetMarkerColor(Colors0[0]);
   histograms2[0]->Scale(1./histograms2[0]->Integral());
   histograms2[0]->Draw("HIST");


  histograms2[0]->GetXaxis()->SetTitleSize(0.06);
  histograms2[0]->GetYaxis()->SetTitleSize(0.06);
  histograms2[0]->GetYaxis()->SetTitleOffset(1.5);
  histograms2[0]->SetTickLength(0.03, "XYZ");
  histograms2[0]->SetLabelOffset(0.002,"X");//0.007
  histograms2[0]->SetLabelOffset(0.007,"Y");
  histograms2[0]->SetLabelSize(0.042, "XYZ");
  histograms2[0]->SetLabelFont(42, "XYZ"); 
  histograms2[0]->SetTitleFont(42, "XYZ");
  histograms2[0]->SetTitleSize(0.05, "XYZ"); 
  histograms2[0]->SetTitleOffset(1.1,"Y");
  histograms2[0]->SetTitleOffset(0.85,"X");
  histograms2[0]->SetNdivisions(505,"XYZ");

  leg2->AddEntry(histograms2[0],Name[0],"L");

    for (unsigned int k = 1; k < nMC; k++)
    {
      StyleAndDraw(histograms2[k], Colors2[k], 2, 1, Datasets[k]);
      leg2->AddEntry(histograms2[k],Name[k],"L");
    }

  histograms2[0]->SetMaximum(0.35);
  histograms2[0]->SetTitle("");
  histograms2[0]->GetXaxis()->SetTitle(xtitle2);
  histograms2[0]->GetYaxis()->SetTitle(ytitle2);
  // histograms2[0]->GetXaxis()->SetTitleSize(0.06);
  // histograms2[0]->GetYaxis()->SetTitleOffset(1.5);
  // histograms2[0]->GetYaxis()->SetTitleSize(0.06);
 leg2->Draw();
 
 //------Start of Copy Paste

  H = pad3->GetWh();
  W = pad3->GetWw();
  l = pad3->GetLeftMargin();
  t = pad3->GetTopMargin();
  r = pad3->GetRightMargin();
  b = pad3->GetBottomMargin();

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextAngle(0);
  latex2.SetTextColor(kBlack);    

  latex2.SetTextFont(42);
  latex2.SetTextAlign(31); 
  latex2.SetTextSize(lumiTextSize*t);    
  latex2.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

      latex2.SetTextFont(cmsTextFont);
      latex2.SetTextAlign(11); 
      latex2.SetTextSize(cmsTextSize*t);    
      latex2.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);

  posX_=0;
  posY_=0;
  iPosX = 3;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  posY_ = 1-t - relPosY*(1-t-b);
  	  if( writeExtraText ) 
	    {
         posX_ =   l +  relPosX*(1-l-r);
         posY_ =   1-t+lumiTextOffset*t;
        int alignY_=3;
         int alignX_=2;
         if( iPosX/10==0 ) alignX_=1;
         if( iPosX==0    ) alignX_=1;
         if( iPosX==0    ) alignY_=1;
         if( iPosX/10==1 ) alignX_=1;
         if( iPosX/10==2 ) alignX_=2;
         if( iPosX/10==3 ) alignX_=3;
         //if( iPosX == 0  ) relPosX = 0.12;
         int align_ = 10*alignX_ + alignY_;
      latex2.SetTextFont(extraTextFont);
      latex2.SetTextSize(extraTextSize*t);
      latex2.SetTextAlign(11);
      latex2.DrawLatex(posX_+0.08, posY_, extraText);
	    }
 c1->Update();

//  // *****************************************************************************

  pad4->cd();

   /////////
  TString htitle3 = "GenRecoAxis_dRmin_";
  TString xtitle3 = "#Delta R_{#chi-H}"; 
  TString ytitle3 = "a.u"; 
  TH1F* histograms3[nMC];
  // int Colors[nMC] = {4,1,2,3,33,5,6,7,8,9,30,38,41,43,46,49};
  int Colors3[nMC]={ColorBlue,ColorRed,ColorNeutral};

for (int i = 0; i < nMC; i++) {
    // Colors3[i]= i+1;
    fX[i]->cd();
    TString histogramTitle = Datasets[i]+"_"+htitle3; // Assuming htitle is an array of strings
    histograms3[i] = (TH1F*)gROOT->FindObject(histogramTitle);
    histograms3[i]->GetXaxis()->SetRangeUser(0,5);
}


/////////
  TLegend* leg3 = new TLegend(0.5,0.6,0.8,0.85);
  leg3->SetHeader("MC #tilde{#mu}^{+}#tilde{#mu}^{-}(j)  with  #tilde{#mu} #rightarrow #mu#tilde{#chi}^{0},  #tilde{#chi}^{0} #rightarrow tds");
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetFillColor(kWhite);
  leg3->SetTextFont(42);
  leg3->SetTextSize(0.03);

   histograms3[0]->SetLineColor(Colors0[0]);
   histograms3[0]->SetLineWidth(2);
   histograms3[0]->SetLineStyle(1);
   histograms3[0]->SetMarkerStyle(76);
   histograms3[0]->SetMarkerSize(1);
   histograms3[0]->SetMarkerColor(Colors0[0]);
   histograms3[0]->Scale(1./histograms3[0]->Integral());
   histograms3[0]->Draw("HIST");


  histograms3[0]->GetXaxis()->SetTitleSize(0.06);
  histograms3[0]->GetYaxis()->SetTitleSize(0.06);
  histograms3[0]->GetYaxis()->SetTitleOffset(1.5);
  histograms3[0]->SetTickLength(0.03, "XYZ");
  histograms3[0]->SetLabelOffset(0.002,"X");//0.007
  histograms3[0]->SetLabelOffset(0.007,"Y");
  histograms3[0]->SetLabelSize(0.042, "XYZ");
  histograms3[0]->SetLabelFont(42, "XYZ"); 
  histograms3[0]->SetTitleFont(42, "XYZ");
  histograms3[0]->SetTitleSize(0.05, "XYZ"); 
  histograms3[0]->SetTitleOffset(1.1,"Y");
  histograms3[0]->SetTitleOffset(0.85,"X");
  histograms3[0]->SetNdivisions(505,"XYZ");


  leg3->AddEntry(histograms3[0],Name[0],"L");
    for (unsigned int k = 1; k < nMC; k++)
    {
      StyleAndDraw(histograms3[k], Colors3[k], 2, 1, Datasets[k]);
      leg3->AddEntry(histograms3[k],Name[k],"L");
    }

  histograms3[0]->SetMaximum(0.45);
  histograms3[0]->SetTitle("");
  histograms3[0]->GetXaxis()->SetTitle(xtitle3);
  // histograms3[0]->GetXaxis()->SetTitleSize(0.06);
  histograms3[0]->GetYaxis()->SetTitle(ytitle3);
  // histograms3[0]->GetYaxis()->SetTitleOffset(1.5);
  // histograms3[0]->GetYaxis()->SetTitleSize(0.06);
 leg3->Draw();
 
 //------Start of Copy Paste

  H = pad4->GetWh();
  W = pad4->GetWw();
  l = pad4->GetLeftMargin();
  t = pad4->GetTopMargin();
  r = pad4->GetRightMargin();
  b = pad4->GetBottomMargin();

  TLatex latex3;
  latex3.SetNDC();
  latex3.SetTextAngle(0);
  latex3.SetTextColor(kBlack);    



  latex3.SetTextFont(42);
  latex3.SetTextAlign(31); 
  latex3.SetTextSize(lumiTextSize*t);    
  latex3.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

      latex3.SetTextFont(cmsTextFont);
      latex3.SetTextAlign(11); 
      latex3.SetTextSize(cmsTextSize*t);    
      latex3.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);

posX_=0;
  posY_=0;
  iPosX = 3;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  posY_ = 1-t - relPosY*(1-t-b);
  	  if( writeExtraText ) 
	    {
         posX_ =   l +  relPosX*(1-l-r);
         posY_ =   1-t+lumiTextOffset*t;
        int alignY_=3;
         int alignX_=2;
         if( iPosX/10==0 ) alignX_=1;
         if( iPosX==0    ) alignX_=1;
         if( iPosX==0    ) alignY_=1;
         if( iPosX/10==1 ) alignX_=1;
         if( iPosX/10==2 ) alignX_=2;
         if( iPosX/10==3 ) alignX_=3;
         //if( iPosX == 0  ) relPosX = 0.12;
         int align_ = 10*alignX_ + alignY_;
      latex3.SetTextFont(extraTextFont);
      latex3.SetTextSize(extraTextSize*t);
      latex3.SetTextAlign(11);
      latex3.DrawLatex(posX_+0.08, posY_, extraText);
	    }
 c1->Update();
 c1->SaveAs("plot_dR.pdf");

}
