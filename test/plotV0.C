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

void plot(TString sample, TString YEAR,TString PLAN, TString SELECTION,TString Channel )
{
 TFile *f1;
TString Prod = "SecInt";
TString Year = YEAR;
// TString Sample = "Data_"+Year;
TString Sample = sample; //"DoubleMuon_UL"+Year+"_MiniAODv2_GT36-v1";
TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/SECINT_"+Sample+".root";
f1 = new TFile(Path);
TString plan = PLAN; // can also be rz
TString Selection = SELECTION; // can also be TrackerMatched + PU25 + PU30 + PU35 + PU40 + PU45 + PU50


 TString xtitle = "x [cm]"; 
 TString ytitle = "y [cm]";  


float xmin = -5.;
float xmax = 5.;
float ymin = -5.;
float ymax = 5.;

if (plan == "rz") {
  xmin = 0.;
  xmax = 120.;
  ymin = 0.;
  ymax = 70.;
  xtitle = "z [cm]";
  ytitle = "r [cm]";
}
else if (plan == "xy") {
  xmin = -5;
  xmax = 5 ;
  ymin = -5;
  ymax = 5;
}


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

/////////
 f1->cd();

  htitle0 = Sample+"_hData_"+Selection+"_"+Channel+"_"+plan+"_";


 TH2F* h0 = (TH2F*)gROOT->FindObject(htitle0);
 h0->Sumw2(); 

       h0->Draw("COL"); 
       h0->SetLineColor(kBlack);
       h0->SetLineStyle(1);
       h0->SetLineWidth(1);//3
       h0->SetTickLength(0.03, "XYZ");
       h0->SetLabelOffset(0.002,"X");//0.007
      //  h0->GetXaxis()->SetLabelSize(0.03);
      //  h0->GetXaxis()->SetLabelOffset(0.001);
       h0->SetLabelOffset(0.007,"Y");
       h0->SetLabelSize(0.022, "XYZ");
       h0->SetLabelFont(42, "XYZ"); 
       h0->SetTitleFont(42, "XYZ");
       h0->SetTitleSize(0.03, "XYZ"); 
       h0->SetTitleOffset(1.1,"Y");
       h0->SetTitleOffset(0.85,"X");
       h0->GetYaxis()->SetTitle(ytitle);
       h0->GetYaxis()->SetTitleColor(1);
       h0->SetNdivisions(505,"XYZ");
       h0->GetXaxis()->SetRangeUser(xmin, xmax); 
       h0->GetYaxis()->SetRangeUser(ymin,ymax);
       h0->GetXaxis()->SetTitle(xtitle);
       h0->SetTitle("");
  // leg = new TLegend(0.70,0.75,0.9,0.85);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.025);
  // leg->SetHeader("MC Samples");
  // leg->AddEntry(h0,"Secondary Vertices","LE");

  // leg->Draw();
  
// *****************************************************************************
 TString cmsText     = "Private Work";//CMS
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = true;
TString extraText   = "";//Privater Work
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.53;//0.6
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.65;//0.75
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "";//137 fb^{-1}
TString lumi_sqrtS = "CMS "+Year+" Data";
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
  float posY_ = 1-t - relPosY*(1-t-b);
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

// if (plan == "xy") {
//   TLatex latex;
//   latex.SetNDC();
//   latex.SetTextAngle(0);
//   latex.SetTextColor(kBlack);    

//   float labelSize = 0.022;

//   latex.SetTextFont(42);
//   latex.SetTextAlign(31); 
//   latex.SetTextSize(labelSize);    
//   latex.DrawLatex(0.93,0.09,"25");
// }

 c1->Update();
 c1->SaveAs(htitle0+".pdf");
}
