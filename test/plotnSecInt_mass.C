#include <iostream>
#include <TROOT.h>
#include "TVectorD.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TEfficiency.h" 
#include "TMath.h" 
#include "TLatex.h"

void plot(std::vector<TString> Sample , TString YEAR, TString SELECTION,TString PU, TString ctau, TString msmu, TString mneu  )
{
 TFile *f1;
TString Prod = "SecInt";
TString Year = YEAR;
TString Selection = SELECTION; // can also be TrackerMatched + PU25 + PU30 + PU35 + PU40 + PU45 + PU50
TString Pu = PU; // can also be PU30 + PU35 + PU40 + PU45 + PU50
TString Suffix = Selection+Pu;

float xmin = 0;
float xmax = 20;

int stati=0;
bool fit= 0;
bool logy=0;

 TLegend* leg;
 TString htitle0;


// *****************************************************************************

// TCanvas *c1 = new TCanvas("c1", "plots",200,0,900,500);
TCanvas *c1 = new TCanvas("c1", "plots",200,0,1000,1100);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);
c1->SetBatch(kTRUE);
TPad* pad1 = new TPad("pad1","This is pad1",0.01,0.01,0.99,1,21);
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
pad1->SetLogz(true);
   pad1->SetTopMargin(0.1);
   pad1->SetBottomMargin(0.1);
   pad1->SetRightMargin(0.1);
   pad1->SetLeftMargin(0.1);
   


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

  pad1->cd();
 //************************
  // !! ---------------------------------
  const int nMC = Sample.size();
  TFile *fX[nMC];
  TString NameMC[nMC] = {
    "c#tau = 0.1 cm",
    "c#tau = 0.3 cm",
    "c#tau = 1.0 cm",
    "c#tau = 3.0 cm",
    "c#tau = 10.0 cm",
    "c#tau = 30.0 cm",
    "c#tau = 100.0 cm"
  }
  for (int i = 0; i < nMC; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/SECINT_"+Sample[i]+".root";
      fX[i] = new TFile(Path);
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

TString xtitle = "# of SecInt"; 
TString ytitle = "a.u";  
TH1F* histograms[nMC];
int Colors[10] = {ColorBlue,ColorOrange,ColorRed,ColorGrey,ColorDarkPurple,ColorBrown,ColorDarkOrange,ColorNeutral,ColorDarkGrey,ColorLightBlue};

for (int i = 0; i < nMC; i++) {
  fX[i]->cd();
  htitle0 = Sample[i]+"_hData_reco_nSecInt_"+Suffix; 
  histograms[i] = (TH1F*)gROOT->FindObject(htitle0);
}

/////////
  std::cout<<"msu : "<<msmu<<" & "<<mneu<<" "<<" ctau : "<<ctau<<std::endl;

  // histograms[0]->Sumw2(); 
  histograms[0]->Draw("hist"); 
  histograms[0]->SetLineColor(Colors[0]);  
  histograms[0]->SetTitle("");
  histograms[0]->GetXaxis()->SetTitle(xtitle);
  histograms[0]->GetXaxis()->SetTitleSize(0.06);
  histograms[0]->GetYaxis()->SetTitle(ytitle);
  histograms[0]->GetYaxis()->SetTitleSize(0.06);
  histograms[0]->GetYaxis()->SetTitleOffset(1.5);
  histograms[0]->SetLineStyle(1);
  histograms[0]->SetLineWidth(2);//3
  histograms[0]->SetTickLength(0.03, "XYZ");
  histograms[0]->SetLabelOffset(0.002,"X");//0.007
  histograms[0]->SetLabelOffset(0.007,"Y");
  histograms[0]->SetLabelSize(0.032, "XYZ");
  histograms[0]->SetLabelFont(42, "XYZ"); 
  histograms[0]->SetTitleFont(42, "XYZ");
  histograms[0]->SetTitleSize(0.04, "XYZ"); 
  histograms[0]->SetTitleOffset(1.1,"Y");
  histograms[0]->SetTitleOffset(0.85,"X");
  histograms[0]->SetNdivisions(505,"XYZ");
  histograms[0]->Scale(1./histograms[0]->Integral());

  TLegend* leg = new TLegend(0.4,0.57,0.85,0.87);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetHeader("Signal : M_{#tilde{#mu}} = "+msmu+" GeV, M_{#tilde{#chi}} = "+mneu+" GeV ");
  leg->AddEntry(histograms[0],NameMC[0],"L");
    for (unsigned int k = 1; k < nMC; k++)
    {
      histograms[k]->SetLineColor(Colors[k]);
      histograms[k]->SetLineStyle(1);
      histograms[k]->SetLineWidth(2);//3
      // histograms[k]->SetFillColor(Colors[k]);
      histograms[k]->Scale(1./histograms[k]->Integral());
      // hs->Add(histograms[k]);
      histograms[k]->Draw("histsame");
      leg->AddEntry(histograms[k],NameMC[k],"L");
    }
   
  histograms[0]->SetMaximum(1);


  leg->Draw();
  
// *****************************************************************************
 TString cmsText     = "Private Work";//CMS
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = true;
TString extraText   = "";//Privater Work
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.43;//0.6
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.55;//0.75
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "";//137 fb^{-1}
TString lumi_sqrtS = "CMS "+Year+" Simulation";
TString lumiText = lumi_13TeV+lumi_sqrtS;
  float H = c1->GetWh();
  float W = c1->GetWw();
  float l = c1->GetLeftMargin();
  float t = c1->GetTopMargin();
  float r = c1->GetRightMargin();
  float b = c1->GetBottomMargin();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t-0.01,lumiText);

  latex.SetTextFont(cmsTextFont);
  latex.SetTextAlign(11); 
  latex.SetTextSize(cmsTextSize*t);    
  latex.DrawLatex(l,1-t+lumiTextOffset*t-0.01,cmsText);

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
      latex.DrawLatex(posX_+0.14, posY_-0.01, extraText);
	    }

 c1->Update();
 c1->SaveAs( msmu+"_"+mneu+"_fixedmass_fctau_"+Suffix+".pdf");
}
