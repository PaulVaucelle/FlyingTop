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

void plot()
{
 TFile *f1, *f2, *f3, *f4, *f5;
 TGraphAsymmErrors *gr0, *gr1, *gr2, *gr3, *gr4;
 


// Parameters
// TString step = "";
// ""
// "_TightWP"
// "_LooseWP"
TString ctau = "030";
TString Dataset = "RPV_2018_ctau"+ctau;
bool UseEta = false;
bool UseR = false;
bool UseDist = true;
TString dof = "";
//Plotting
int stati=0;
bool fit= 0;
bool logy=0;

TString ytitle = "entries"; 
TString xtitle="";
int nbin = 200; 
float binmin =   0;//0
float binmax = 200;//100

//--------------//
// 65,-6.5,6.5

 TLegend* leg;
 TLegend*  leg2;
 TLegend*  leg3;
 TString htitle0, htitle1, htitle2, htitle3, htitle4, htitle5, htitle6; 

if (UseR)
    {
        htitle0 = "hSim_Hemi_Vtx_r_noSel";
        htitle1 = "hData_Hemi_Vtx_r_GoodRecoVtx"; 
        htitle2 = "hData_Hemi_Vtx_r_Ping";
        htitle3 = "hData_Hemi_Vtx_r_TightPing";
        htitle4 = "hData_Hemi_Vtx_r_LoosePing";
        htitle5 = "hData_Hemi_Vtx_r_TightPing_step1";
        htitle6 = "hData_Hemi_Vtx_r_LoosePing_step3";
        xtitle = "Transverse decay length [cm]";
        nbin = 50;
        binmin = 0;
        binmax = 100;
    }
if (UseEta)
    {
        htitle0 = "hSim_Hemi_Vtx_eta_noSel";
        htitle1 = "hData_Hemi_Vtx_eta_GoodRecoVtx"; 
        htitle2 = "hData_Hemi_Vtx_eta_Ping";
        htitle3 = "hData_Hemi_Vtx_eta_TightPing";
        htitle4 = "hData_Hemi_Vtx_eta_LoosePing";
        htitle5 = "hData_Hemi_Vtx_eta_TightPing_step1";
        htitle6 = "hData_Hemi_Vtx_eta_LoosePing_step3";
        xtitle = "#eta";
        nbin = 26;
        binmin = -6.5;
        binmax = 6.5;
    }
if (UseDist)
    {
        htitle0 = "hSim_Hemi_Vtx_dist_noSel"; 
        htitle1 = "hData_Hemi_Vtx_dist_GoodRecoVtx"; 
        htitle2 = "hData_Hemi_Vtx_dist_Ping";
        htitle3 = "hData_Hemi_Vtx_dist_TightPing";
        htitle4 = "hData_Hemi_Vtx_dist_LoosePing";
        htitle5 = "hData_Hemi_Vtx_dist_TightPing_step1";
        htitle6 = "hData_Hemi_Vtx_dist_LoosePing_step3";
        xtitle  = "decay length [cm]";
        nbin = 50;
        binmin = 0;
        binmax = 100; 
    }


//  htitle0 = "hSim_Hemi_LLP_dist"; 
//  htitle1 = "hSim_Hemi_LLP_dist_chiOK"; 
//  htitle2 = "hSim_Hemi_LLP_dist_ping"; 
// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",200,0,1500,1000);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

TPad* pad1 = new TPad("pad1","This is pad1",0.01,0.01,0.99,1,21);
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.1);
   pad1->SetBottomMargin(0.15);
   pad1->SetRightMargin(0.1);
   pad1->SetLeftMargin(0.1);

// TPad* rap1 = new TPad("rap1","This is rap1",0.04,0.01,0.98,0.45,21);
// rap1->SetFillColor(0);
// rap1->SetBorderMode(0);
// rap1->SetFrameFillColor(10);
// rap1->Draw();
// rap1->SetLogy(0);
//    rap1->SetTopMargin(0.1);
//    rap1->SetBottomMargin(0.1);
//    rap1->SetRightMargin(0.15);
//    rap1->SetLeftMargin(0.15);

TCanvas *c2 = new TCanvas("c2", "plots",200,0,1500,1000);
c2->SetFillColor(10);
c2->SetFillStyle(4000);
c2->SetBorderSize(2);

TPad* pad2 = new TPad("pad2","This is pad2",0.01,0.01,0.99,1,21);
pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(logy);
   pad2->SetTopMargin(0.1);
   pad2->SetBottomMargin(0.16);
   pad2->SetRightMargin(0.1);
   pad2->SetLeftMargin(0.1);

// TPad* rap2 = new TPad("rap2","This is rap2",0.04,0.01,0.98,0.45,21);
// rap2->SetFillColor(0);
// rap2->SetBorderMode(0);
// rap2->SetFrameFillColor(10);
// rap2->Draw();
// rap2->SetLogy(0);
//    rap2->SetTopMargin(0.10);
//    rap2->SetBottomMargin(0.16);
//    rap2->SetRightMargin(0.05);
//    rap2->SetLeftMargin(0.10);


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

Float_t r1 = 0.246;
Float_t g1 = 0.563;
Float_t b1 = 0.852;
TColor color1 = TColor(301,r1, g1, b1);
// color1.SetRGB(r1, g1, b1);
Int_t ColorBlue = color1.GetNumber();
Float_t r3 = 0.739;
Float_t g3 = 0.122;
Float_t b3 = 0.004;
TColor color3 = TColor(303,r3, g3, b3);
Int_t ColorRed = color3.GetNumber();

// *****************************************************************************
  TH1F*    h_purity = new TH1F("h_purity","h_purity",nbin,binmin,binmax);
  TH1F*    h_efficacity = new TH1F("h_efficacity","h_efficacity",nbin,binmin,binmax);

  TH1F*    h_dist = new TH1F("h_dist","h_dist",nbin,binmin,binmax);
  TH1F*    h_chi2 = new TH1F("h_chi2","h_chi2",nbin,binmin,binmax);
  TH1F*    h_ping = new TH1F("h_ping","h_ping",nbin,binmin,binmax);


// if(UseDist){nbin=20;binmin=0;binmax=100;}
// if(UseEta){nbin=260;binmin=-6.5;binmax=6.5;}
// if(UseGen){nbin=10;binmin=0;binmax=100;}
// if(UseR){nbin=20;binmin=0;binmax=100;}


  const int nHisto = 7;
  TString HTITLE[nHisto]={htitle0,htitle1,htitle2,htitle3,htitle4,htitle5,htitle6};       
  TFile *fX = new TFile("../Signal_2018/histoGen_"+Dataset+".root");
        
  TH1F* histograms[nHisto];
  fX->cd();
  for (int j = 0 ; j < nHisto ; j++){
    
    TString histogramTitle = HTITLE[j]; 
    // std::cout<<"histo title : "<<histogramTitle<<std::endl;
    TString NAME = Dataset+"_"+histogramTitle;
    // std::cout<<"histo name : "<<NAME<<std::endl;
    histograms[j] = (TH1F*)gROOT->FindObject(NAME);
    // histograms[j]->Sumw2();
  }
  
  pad1->cd();

  h_dist->Add(h_dist,histograms[0],0,1);
  h_chi2->Add(h_chi2,histograms[1],0,1);
  h_ping->Add(h_ping,histograms[3],0,1);

h_efficacity->Divide(h_ping,h_dist,1,1,"B");
h_purity->Divide(h_ping,h_chi2,1,1,"B");

h_purity->Draw("");
// h_purity->GetXaxis()->SetRangeUser(0,50);
 h_purity->GetXaxis()->SetLabelSize(0.050);
 h_purity->GetYaxis()->SetLabelSize(0.050);
 h_purity->SetLabelFont(42, "XYZ"); 
 h_purity->SetTitleFont(42, "XYZ");
 h_purity->GetXaxis()->SetTitleColor(1);
 h_purity->GetYaxis()->SetTitleColor(1);
 h_purity->GetXaxis()->SetTitleSize(0.06);
 h_purity->GetYaxis()->SetTitleSize(0.06);
 h_purity->GetXaxis()->SetTitle(xtitle);
 h_purity->GetYaxis()->SetTitle("");
 h_purity->SetTitleOffset(0.9,"X");
 h_purity->SetTitleOffset(0.6,"Y");
 h_purity->SetLabelOffset(0.01,"Y");
 h_purity->SetNdivisions(509, "XYZ"); 
 h_purity->SetLineColor(kRed);
 h_purity->SetMarkerColor(kRed);
 h_purity->SetLineStyle(1);
 h_purity->SetLineWidth(3);
 h_purity->SetMinimum(0.);
 h_purity->SetMaximum(1.);
 h_purity->SetTitle("");

 h_efficacity->Draw("same"); 
 h_efficacity->SetLineColor(ColorBlue);
 h_efficacity->SetMarkerColor(ColorBlue);
 h_efficacity->SetLineStyle(1);
 h_efficacity->SetLineWidth(3);


  leg = new TLegend(0.6,0.6,0.8,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetHeader("MC"+dof+" #tilde{#chi}^{ 0} #rightarrow tqq");
  leg->AddEntry(h_efficacity,"Vtx reco efficiency","LE");
  leg->AddEntry(h_purity,"#tilde{#chi}^{ 0} Vtx purity","LE");
  leg->Draw();

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

//  rap1->cd();


  // c1->Update();


//---------------------------//
  c2->cd();
  pad2->cd();

//   TH1F*    h_efficacity4 = new TH1F("h_efficacity4","h_efficacity4",nbin,binmin,binmax);
//   TH1F*    h_efficacity5 = new TH1F("h_efficacity5","h_efficacity5",nbin,binmin,binmax);

//   TH1F*    h_dist3 = new TH1F("h_dist3","h_dist3",nbin,binmin,binmax);
//   TH1F*    h_ping4 = new TH1F("h_ping4","h_ping4",nbin,binmin,binmax);
//   TH1F*    h_ping5 = new TH1F("h_ping5","h_ping5",nbin,binmin,binmax);

//   h_dist3->Add(h_dist3,histograms[2],0,1);
//   h_ping4->Add(h_ping4,histograms[5],0,1);
//   h_ping5->Add(h_ping5,histograms[6],0,1);

//   h_efficacity4->Divide(h_ping4,h_dist3,1,1,"B");
//   h_efficacity5->Divide(h_ping5,h_dist3,1,1,"B");
//   // h_purity->Divide(h_ping,h_chi2,1,1,"B");
// h_efficacity4->Draw("");
//  h_efficacity4->GetXaxis()->SetLabelSize(0.050);
//  h_efficacity4->GetYaxis()->SetLabelSize(0.050);
//  h_efficacity4->SetLabelFont(42, "XYZ"); 
//  h_efficacity4->SetTitleFont(42, "XYZ");
//  h_efficacity4->GetXaxis()->SetTitleColor(1);
//  h_efficacity4->GetYaxis()->SetTitleColor(1);
//  h_efficacity4->GetXaxis()->SetTitleSize(0.06);
//  h_efficacity4->GetYaxis()->SetTitleSize(0.06);
//  h_efficacity4->GetXaxis()->SetTitle(xtitle);
//  h_efficacity4->GetYaxis()->SetTitle("");
//  h_efficacity4->SetTitleOffset(0.9,"X");
//  h_efficacity4->SetTitleOffset(0.6,"Y");
//  h_efficacity4->SetLabelOffset(0.01,"Y");
//  h_efficacity4->SetNdivisions(509, "XYZ"); 
//  h_efficacity4->SetLineColor(ColorBlue);
//  h_efficacity4->SetLineStyle(1);
//  h_efficacity4->SetLineWidth(3);
//  h_efficacity4->SetMinimum(0.);
//  h_efficacity4->SetMaximum(1.);
//  h_efficacity4->SetTitle("");

//  h_efficacity5->Draw("same"); 
//  h_efficacity5->SetLineColor(kRed);
//  h_efficacity5->SetLineStyle(1);
//  h_efficacity5->SetLineWidth(3);

//   leg3 = new TLegend(0.15,0.2,0.25,0.45);
//   leg3->SetBorderSize(0);
//   leg3->SetFillColor(kWhite);
//   leg3->SetTextFont(42);
//   leg3->SetTextSize(0.035);
//   leg3->SetHeader("MC"+dof+" #tilde{#chi}^{0} #rightarrow tqq efficiency");
//   leg3->AddEntry(h_efficacity4," Tight Step 1","L");
//   leg3->AddEntry(h_efficacity5," Loose Step 3","L");
//   // leg2->AddEntry(h_purity," purity","L");
//   leg3->Draw();

  TH1F*    h_efficacity2 = new TH1F("h_efficacity2","h_efficacity2",nbin,binmin,binmax);
  TH1F*    h_efficacity3 = new TH1F("h_efficacity3","h_efficacity3",nbin,binmin,binmax);

  TH1F*    h_dist2 = new TH1F("h_dist2","h_dist2",nbin,binmin,binmax);
  TH1F*    h_ping2 = new TH1F("h_ping2","h_ping2",nbin,binmin,binmax);
  TH1F*    h_ping3 = new TH1F("h_ping3","h_ping3",nbin,binmin,binmax);

  h_dist2->Add(h_dist2,histograms[2],0,1);
  h_ping2->Add(h_ping2,histograms[3],0,1);
  h_ping3->Add(h_ping3,histograms[4],0,1);

  h_efficacity2->Divide(h_ping2,h_dist2,1,1,"B");
  h_efficacity3->Divide(h_ping3,h_dist2,1,1,"B");
  // h_purity->Divide(h_ping,h_chi2,1,1,"B");
h_efficacity2->Draw("");
 h_efficacity2->GetXaxis()->SetLabelSize(0.050);
 h_efficacity2->GetYaxis()->SetLabelSize(0.050);
 h_efficacity2->SetLabelFont(42, "XYZ"); 
 h_efficacity2->SetTitleFont(42, "XYZ");
 h_efficacity2->GetXaxis()->SetTitleColor(1);
 h_efficacity2->GetYaxis()->SetTitleColor(1);
 h_efficacity2->GetXaxis()->SetTitleSize(0.06);
 h_efficacity2->GetYaxis()->SetTitleSize(0.06);
 h_efficacity2->GetXaxis()->SetTitle(xtitle);
 h_efficacity2->GetYaxis()->SetTitle("");
 h_efficacity2->SetTitleOffset(0.9,"X");
 h_efficacity2->SetTitleOffset(0.6,"Y");
 h_efficacity2->SetLabelOffset(0.01,"Y");
 h_efficacity2->SetNdivisions(509, "XYZ"); 
 h_efficacity2->SetLineColor(ColorBlue);
 h_efficacity2->SetMarkerColor(ColorBlue);
 h_efficacity2->SetLineStyle(1);
 h_efficacity2->SetLineWidth(3);
 h_efficacity2->SetMinimum(0.);
 h_efficacity2->SetMaximum(1.);
 h_efficacity2->SetTitle("");

 h_efficacity3->Draw("same"); 
 h_efficacity3->SetLineColor(kRed);
 h_efficacity3->SetMarkerColor(kRed);
 h_efficacity3->SetLineStyle(1);
 h_efficacity3->SetLineWidth(3);

  leg2 = new TLegend(0.15,0.4,0.25,0.6);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(kWhite);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.035);
  leg2->SetHeader("Step contribution");
  leg2->AddEntry(h_efficacity2," Tight","L");
  leg2->AddEntry(h_efficacity3," Loose","L");
  // leg2->AddEntry(h_purity," purity","L");
  leg2->Draw();


   H = c2->GetWh();
   W = c2->GetWw();
   l = c2->GetLeftMargin();
   t = c2->GetTopMargin();
   r = c2->GetRightMargin();
   b = c2->GetBottomMargin();

  TLatex latex2;
  latex2.SetNDC();
  latex2.SetTextAngle(0);
  latex2.SetTextColor(kBlack);    

   extraTextSize = extraOverCmsTextSize*cmsTextSize;

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
      latex2.SetTextFont(extraTextFont);
      latex2.SetTextSize(extraTextSize*t);
      latex2.SetTextAlign(11);
      latex2.DrawLatex(posX_+0.08, posY_, extraText);
	    }
c1->SaveAs("./plotVtxRecoEffi_"+ctau+".pdf");
c2->SaveAs("./plotVtxRecoEffi_"+ctau+"_StepContribution.pdf");
}
