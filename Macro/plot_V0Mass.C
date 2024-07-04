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
#include "TLatex.h"

void plot()
{
 TFile *f1, *f2, *f3, *f4;
 TGraphAsymmErrors *gr0, *gr1, *gr2, *gr3, *gr4;
 TString Prod = "";
 f1 = new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/histofile_DataA.root");

    gROOT->LoadMacro("./tdrstyle.C");
gROOT->LoadMacro("./CMS_lumi.C");
int stati=0;
bool fit= 0;
bool logy=0;

 TLegend* leg;
 TString htitle0, htitle1, htitle2, htitle3;
 TString dof = "10";

//  TString htitle = "hSim_Hemi_LLP_dist_ping"; 
 TString xtitle = "decay length [cm]"; 
 TString ytitle = "# of entries"; 
 int nbin = 201; 
 float binmin =   0.;
 float binmax = 6.;

// *****************************************************************************

// TCanvas *c1 = new TCanvas("c1", "plots",200,0,900,500);
TCanvas *c1 = new TCanvas("c1", "plots",200,0,1500,1000);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

TPad* pad1 = new TPad("pad1","This is pad1",0.04,0.5,0.49,1,21);
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.15);
   pad1->SetBottomMargin(0.16);
   pad1->SetRightMargin(0.05);
   pad1->SetLeftMargin(0.10);

TPad* pad2 = new TPad("pad2","This is pad2",0.51,0.5,0.99,1,21);
pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(0);
   pad2->SetTopMargin(0.15);
   pad2->SetBottomMargin(0.16);
   pad2->SetRightMargin(0.05);
   pad2->SetLeftMargin(0.10);

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
setTDRStyle();
htitle2 = "hData_reco_K0_mass_noSel__DataA";  
htitle3 = "hData_reco_L0_mass_noSel__DataA";

/////////
 f1->cd();

 TH1F* h2 = (TH1F*)gROOT->FindObject(htitle2);
 h2->Sumw2(); 

       h2->Draw("LE"); 
       h2->SetLineColor(kBlack);
       h2->SetLineStyle(1);
       h2->SetLineWidth(3);
       h2->SetTickLength(0.03, "XYZ");
       h2->SetLabelOffset(0.007,"X");
       h2->SetLabelOffset(0.007,"Y");
       h2->SetLabelSize(0.03, "XYZ");
       h2->SetLabelFont(42, "XYZ"); 
       h2->SetTitleFont(42, "XYZ");
       h2->SetTitleSize(0.055, "XYZ"); 
       h2->SetTitleOffset(1.,"Y");
       h2->GetYaxis()->SetTitle("Entries");
       h2->GetYaxis()->SetTitleColor(1);
       h2->GetXaxis()->SetTitleSize(0.06);
       h2->SetNdivisions(509,"XYZ");
       h2->SetMinimum(0.); 
       h2->GetXaxis()->SetTitle(" #pi^{+}#pi^{-} Mass [GeV]");
//        h0->SetMaximum(0.06); 


  // leg = new TLegend(0.75,0.80,0.85,0.85);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.04);
  // // leg->AddEntry(h2,"Data","LE");
  // leg->Draw();

        TLatex *t = new TLatex(0.42,15900,"Private Work");
        t->SetTextFont(61);
        t->SetTextAlign(11);
        float fac = c1->GetTopMargin();
        t->SetTextSize(fac*0.5);
        t->Draw("same");

        TLatex *t3 = new TLatex(0.46,15900,"CMS 2018 Data");
        t3->SetTextFont(52);
        t3->SetTextAlign(11);
        t3->SetTextSize(fac*0.3);
        t3->Draw();

 // *****************************************************************************
 pad2->cd();
setTDRStyle();
 TH1F* h3 = (TH1F*)gROOT->FindObject(htitle3);
 h3->Sumw2(); 

       h3->Draw("LE"); 
       h3->SetLineColor(kBlack);
       h3->SetLineStyle(1);
       h3->SetLineWidth(3);
       h3->SetTickLength(0.03, "XYZ");
       h3->SetLabelOffset(0.007,"X");
       h3->SetLabelOffset(0.007,"Y");
       h3->SetLabelSize(0.03, "XYZ");
       h3->SetLabelFont(42, "XYZ"); 
       h3->SetTitleFont(42, "XYZ");
       h3->SetTitleSize(0.055, "XYZ"); 
       h3->SetTitleOffset(1.,"Y");
       h3->GetXaxis()->SetTitle(" p#pi^{-} / #bar{p}#pi^{+} Mass [GeV]");
       h3->GetYaxis()->SetTitle("Entries");

       h3->SetNdivisions(509,"XYZ");
       h3->SetMinimum(0.); 



  //        leg = new TLegend(0.75,0.80,0.85,0.85);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.04);
  // // leg->SetHeader("MC Samples");
  // // leg->AddEntry(h3,"Data","LE");
  // leg->Draw();
 

        TLatex *t4 = new TLatex(1.06,8000,"Private Work");
        t4->SetTextFont(61);
        t4->SetTextAlign(11);

        t4->SetTextSize(fac*0.5);
        t4->Draw();


        TLatex *t5 = new TLatex(1.09,7950,"CMS 2018 Data");
        t5->SetTextFont(52);
        t5->SetTextAlign(11);
        t5->SetTextSize(fac*0.3);
        t5->Draw();

  // *****************************************************************************
 c1->Update();
 c1->SaveAs("V0Mass.pdf");
}
