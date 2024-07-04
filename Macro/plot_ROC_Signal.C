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
 TFile *f1;
 TGraphAsymmErrors *gr0, *gr1, *gr2, *gr3, *gr4;
 
//  f1 = new TFile("../output/h_UDD_bgctau50_smu275_snu225_Neu.root");
 f1 = new TFile("../output/h_UDD_bgctau50_smu275_snu225_TT_63k.root");

// h_UDD_bgctau10_smu250_snu200.root");
// h_UDD_bgctau30_smu300_snu250.root");
// h_UDD_bgctau50_smu275_snu225_TTNI.root");
// h_UDD_bgctau50_smu275_snu225_Neu.root");
// h_UDD_bgctau70_smu250_snu200.root");
int stati=0;
bool fit= 0;
bool logy=0;

 TLegend* leg;
 TString htitle0, htitle1, htitle2, htitle3;
 TString dof = "50";

 TString xtitle = "signal efficiency"; 
 TString ytitle = "Background rejection"; 
 int nbin = 201; 
 float binmin =   0.;
 float binmax = 6.;

// *****************************************************************************

// TCanvas *c1 = new TCanvas("c1", "plots",200,0,900,500);
TCanvas *c1 = new TCanvas("c1", "plots",200,0,1000,1000);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

TPad* pad1 = new TPad("pad1","This is pad1",0.04,0.4,0.5,1,21);
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.1);
   pad1->SetBottomMargin(0.16);
   pad1->SetRightMargin(0.05);
   pad1->SetLeftMargin(0.15);

TPad* pad2 = new TPad("pad2","This is pad2",0.51,0.4,0.95,1,21);
pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(logy);
   pad2->SetTopMargin(0.1);
   pad2->SetBottomMargin(0.16);
   pad2->SetRightMargin(0.05);
   pad2->SetLeftMargin(0.2);

TPad* pad3 = new TPad("pad3","This is pad3",0.51,0.01,0.95,0.4,21);
pad3->SetFillColor(0);
pad3->SetBorderMode(0);
pad3->SetFrameFillColor(10);
pad3->Draw();
pad3->SetLogy(logy);
   pad3->SetTopMargin(0.1);
   pad3->SetBottomMargin(0.16);
   pad3->SetRightMargin(0.05);
   pad3->SetLeftMargin(0.2);

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
htitle0 = "hBkfEff_SigEff";
htitle1 = "hSignalEff";
htitle2 = "hBkgEff";
htitle3 = "h_sig";
////////
 f1->cd();
 TH2F* h0 = (TH2F*)gROOT->FindObject(htitle0);
 h0->Sumw2();
 TH1F* h1 = (TH1F*)gROOT->FindObject(htitle1);
 h1->Sumw2();
 TH1F* h2 = (TH1F*)gROOT->FindObject(htitle2);
 h2->Sumw2();
 TH1F* h3 = (TH1F*)gROOT->FindObject(htitle3);
 h3->Sumw2();
//  TH2F* h4 = new TH2F("h4","",201,-1,1);
//       h4->Add(h0,h4,1/h0->Integral(1,201+1),0);
       h0->Draw(""); 
       h0->SetLineColor(kBlack);
       h0->SetLineStyle(1);
       h0->SetLineWidth(3);
       h0->SetTickLength(0.03, "XYZ");
       h0->SetLabelOffset(0.007,"X");
       h0->SetLabelOffset(0.007,"Y");
       h0->SetLabelSize(0.045, "XYZ");
       h0->SetLabelFont(42, "XYZ"); 
       h0->SetTitleFont(42, "XYZ");
       h0->SetTitleSize(0.055, "XYZ"); 
       h0->SetTitleOffset(1.,"Y");
       h0->GetYaxis()->SetTitle(ytitle);
       h0->GetYaxis()->SetTitleColor(1);
       h0->SetNdivisions(509,"XYZ");
       h0->SetMinimum(0.); 
       h0->GetXaxis()->SetTitle(xtitle);
//        h0->SetMaximum(0.06); 

  leg = new TLegend(0.2,0.25,0.35,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("ROC Curve");
  leg->AddEntry(h0,"MC"+dof+"cm","LE");
leg->Draw();

pad2->cd();

       h1->Draw("HIST"); 
       h1->SetLineColor(kBlue);
       h1->SetLineStyle(1);
       h1->SetLineWidth(3);
       h1->SetTickLength(0.03, "XYZ");
       h1->SetLabelOffset(0.007,"X");
       h1->SetLabelOffset(0.007,"Y");
       h1->SetLabelSize(0.045, "XYZ");
       h1->SetLabelFont(42, "XYZ"); 
       h1->SetTitleFont(42, "XYZ");
       h1->SetTitleSize(0.055, "XYZ"); 
       h1->SetTitleOffset(1.,"Y");
       h1->GetYaxis()->SetTitle("efficiency");
       h1->GetYaxis()->SetTitleColor(1);
       h1->SetNdivisions(509,"XYZ");
       h1->SetMinimum(0.); 
       h1->GetXaxis()->SetTitle("MVA score");

       h2->Draw("HISTsame"); 
       h2->SetLineColor(kRed);
       h2->SetLineStyle(1);
       h2->SetLineWidth(3);

         leg = new TLegend(0.6,0.75,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  // leg->SetHeader("ROC Curve");
  leg->AddEntry(h1,"Signal Eff","L");
  leg->AddEntry(h2,"Bkg eff","L");
leg->Draw();


pad3->cd();
       
       h3->Draw("HIST"); 
       h3->SetLineColor(kBlack);
       h3->SetLineStyle(1);
       h3->SetLineWidth(3);
       h3->SetTickLength(0.03, "XYZ");
       h3->SetLabelOffset(0.007,"X");
       h3->SetLabelOffset(0.007,"Y");
       h3->SetLabelSize(0.045, "XYZ");
       h3->SetLabelFont(42, "XYZ"); 
       h3->SetTitleFont(42, "XYZ");
       h3->SetTitleSize(0.055, "XYZ"); 
       h3->SetTitleOffset(1.,"Y");
       h3->GetYaxis()->SetTitle("#frac{S}{sqrt(S+B)}");
       h3->GetYaxis()->SetTitleColor(1);
       h3->SetNdivisions(509,"XYZ");
       h3->SetMinimum(0.); 
       h3->GetXaxis()->SetTitle("MVA score");


  leg = new TLegend(0.25,0.25,0.35,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  // leg->SetHeader("Significance");
  leg->AddEntry(h3,"Significance");
leg->Draw();
  // leg = new TLegend(0.40,0.64,0.80,0.89);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.05);
  // leg->SetHeader("MC"+dof+" #tilde{#chi}^{0} #rightarrow tqq");
  // leg->AddEntry(h0," all","LE");
  // leg->AddEntry(h1," reco Vtx (0 < #chi^{2}/ndf < 10)","LE");
  // leg->AddEntry(h2," reco&matched Vtx (#delta(L)/L  < 0.1)","LE");
  // leg->Draw();



  // leg->AddEntry(h2,"#Delta R GenReco","LE");
  leg->Draw();
  
// *****************************************************************************

 c1->Update();
}
