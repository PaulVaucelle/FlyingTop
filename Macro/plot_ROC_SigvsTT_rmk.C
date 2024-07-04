#include <iostream>
#include <TROOT.h>
#include "TVectorD.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TEfficiency.h" 
#include "TMath.h" 

void plot()
{
 TFile *f1;
 TGraphAsymmErrors *gr0, *gr1, *gr2, *gr3, *gr4;
 TString Prod = "PROD_CSI_10_06_2024";
 TString Sample = "TrackAna_RPV_2018_ctau100.root";
 TString Sample2 = "TrackAna_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root";
 f1 = new TFile("../"+Prod+"/"+Sample);//signal
   f2 = new TFile("../"+Prod+"/"+Sample2);//signal

 TString SaveFile = "plot_ROC_SigvsTT_rmk";
 
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

 TString xtitle = "Signal efficiency"; 
 TString ytitle = "Background rejection"; 
 int nbin = 201; 
 float binmin =   0.;
 float binmax = 6.;

// *****************************************************************************

// TCanvas *c1 = new TCanvas("c1", "plots",200,0,900,500);
TCanvas *c1 = new TCanvas("c1", "plots",200,0,1400,1000);
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

   TPad* pad4 = new TPad("pad4","This is pad4",0.04,0.01,0.5,0.4,21);
pad4->SetFillColor(0);
pad4->SetBorderMode(0);
pad4->SetFrameFillColor(10);
pad4->Draw();
pad4->SetLogy(logy);
   pad4->SetTopMargin(0.1);
   pad4->SetBottomMargin(0.16);
   pad4->SetRightMargin(0.05);
   pad4->SetLeftMargin(0.2);

// TPad* pad4 = new TPad("pad4","This is pad4",0.04,0.01,0.5,0.4,21);
// pad4->SetFillColor(0);
// pad4->SetBorderMode(0);
// pad4->SetFrameFillColor(10);
// pad4->Draw();
// pad4->SetLogy(logy);
//    pad4->SetTopMargin(0.1);
//    pad4->SetBottomMargin(0.16);
//    pad4->SetRightMargin(0.05);
//    pad4->SetLeftMargin(0.2);


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
pad1->SetGrid();
htitle0 = "hBkfEff_SigEff";
htitle1 = "hSignalEff";//hSignalEff
htitle2 = "hBkgEff";//hBkgEff
htitle3 = "h_sig";
////////
 f1->cd();
//  TH2F* h0 = (TH2F*)gROOT->FindObject(htitle0);
//  h0->Sumw2();
//  TH1F* h3 = (TH1F*)gROOT->FindObject(htitle3);
//  h3->Sumw2();
 TH1F* h1 = (TH1F*)gROOT->FindObject(htitle1);
 h1->Sumw2();
//  f2->cd();
 TH1F* h2 = (TH1F*)gROOT->FindObject(htitle2);//htitle2
 h2->Sumw2();
   h1->Scale(1./31.);
   h2->Scale(1./31.);

 

//---------------------pad1-------------------------//
TH2F*    hBkfEff_SigEff = new TH2F("hBkfEff_SigEff","hBkfEff_SigEff",201,0,1,201,0,1);
TH2F*    h_sig = new TH2F("h_sig","h_sig",201,-1,1,201,0,1);
int size = 200;
float sig = 0;
float max = 0.;
float idxmax = 0.;
float bdtcut2 = 0.;
float bdtcut3 = 0;
bool letBdf =  true;
bool letS = true;
for (int k =1; k < size+1 ; k++)
    {
        float  Sbin = h1->GetBinContent(k);
        float  Bbin = h2->GetBinContent(k);
        hBkfEff_SigEff->Fill(Sbin,1-Bbin);
        if (Sbin+Bbin>0)
            {
               sig = Sbin/sqrt(Sbin+Bbin);
            }
        float index = -1+k*(2./201.);
        if (Bbin<0.001 && letBdf){bdtcut2 = index;letBdf = false;}
        if ((1-Sbin)>0.001 && letS){bdtcut3 = index;letS = false;}
        if (sig>max){max=sig;idxmax=index;}
        h_sig->Fill(index,sig);
    }
    
   //   std::cout<<"//-----------------------------------\\ "<<std::endl;
   //   std::cout<<"// Significance max : "<<max<<std::endl;
   //   std::cout<<"// BDT cut : "<<idxmax<<std::endl;

      // std::cout<<"// BDT cut  for Seff of 0.001 "<<bdtcut3<<std::endl;
     std::cout<<"\\-----------------------------------// "<<std::endl;
       hBkfEff_SigEff->Draw("colz"); 
       hBkfEff_SigEff->SetLineColor(kBlack);
       hBkfEff_SigEff->SetLineStyle(1);
       hBkfEff_SigEff->SetLineWidth(3);
       hBkfEff_SigEff->SetTickLength(0.03, "XYZ");
       hBkfEff_SigEff->SetLabelOffset(0.007,"X");
       hBkfEff_SigEff->SetLabelOffset(0.007,"Y");
       hBkfEff_SigEff->SetLabelSize(0.045, "XYZ");
       hBkfEff_SigEff->SetLabelFont(42, "XYZ"); 
       hBkfEff_SigEff->SetTitleFont(42, "XYZ");
       hBkfEff_SigEff->SetTitleSize(0.055, "XYZ"); 
       hBkfEff_SigEff->SetTitleOffset(1.,"Y");
       hBkfEff_SigEff->GetYaxis()->SetTitle(ytitle);
       hBkfEff_SigEff->GetYaxis()->SetTitleColor(1);
       hBkfEff_SigEff->SetNdivisions(509,"XYZ");
       hBkfEff_SigEff->SetMinimum(0.); 
       hBkfEff_SigEff->GetXaxis()->SetTitle(xtitle);


//   leg = new TLegend(0.2,0.25,0.35,0.45);
//   leg->SetBorderSize(0);
//   leg->SetFillColor(kWhite);
//   leg->SetTextFont(42);
//   leg->SetTextSize(0.05);
//   leg->SetHeader("ROC Curve");
//   leg->AddEntry(hBkfEff_SigEff,"c#tau = 10.0 cm","LE");
// leg->Draw();


//----------------------pad2--------------------------//
pad2->cd();
pad2->SetGrid();
//  f2->cd();
 TH1F* h30 = (TH1F*)gROOT->FindObject(htitle2);
   h30->Sumw2();
f2->cd();
TH1F* h20 = (TH1F*)gROOT->FindObject(htitle2);
   // h30->Scale(1./31.);
float temp_Sbin0 = 0 ;
float temp_Bbin0 = 0 ;

float temp_Sbin85 = 0 ;
float temp_Bbin85 = 0 ;

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
for (int k =0; k < size ; k++)
    {
        float index = -1+k*(2./200.);
        if (index==0){temp_Sbin0=h1->GetBinContent(k);temp_Bbin0=h30->GetBinContent(k);}
        if (index>0.938 && index <0.942){temp_Sbin85=h1->GetBinContent(k);temp_Bbin85=h30->GetBinContent(k);}

    }

      //  h2->Draw("HISTsame"); 
      //  h2->SetLineColor(kRed);
      //  h2->SetLineStyle(1);
      //  h2->SetLineWidth(3);

       h30->Draw("HISTsame"); 
       h30->SetLineColor(kRed+2);
       h30->SetLineStyle(1);
       h30->SetLineWidth(3);

      h20->Draw("HISTsame");
      h20->SetLineColor(kRed);
      h30->SetLineStyle(1);
      h30->SetLineWidth(3);

         leg = new TLegend(0.68,0.7,0.87,0.88);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.02);
  // leg->SetHeader("ROC Curve");
  leg->AddEntry(h1,"Signal Tracks","L");
   leg->AddEntry(h30,"Bkg Tracks in Signal Events","L");
leg->AddEntry(h20,"Bkg Tracks in TTTo2L2Nu Events","L");

leg->Draw();

     std::cout<<"//-----------------------------------\\ "<<std::endl;
     std::cout<<"// BDT score 0 SignalEff and BkgEff : "<< temp_Sbin0<<" \\ "<<temp_Bbin0<<std::endl;
     std::cout<<"// BDT score 0.94 SignalEff and BkgEff : "<< temp_Sbin85<<" \\ "<<temp_Bbin85<<std::endl;
     std::cout<<"//-----------------------------------\\ "<<std::endl;
//----------------------pad3--------------------------//
pad3->cd();
       pad3->SetGrid();
       h_sig->Draw("COLZ"); 
       h_sig->SetLineColor(kBlue);
       h_sig->SetLineStyle(1);
       h_sig->SetLineWidth(3);
       h_sig->SetTickLength(0.03, "XYZ");
       h_sig->SetLabelOffset(0.007,"X");
       h_sig->SetLabelOffset(0.007,"Y");
       h_sig->SetLabelSize(0.045, "XYZ");
       h_sig->SetLabelFont(42, "XYZ"); 
       h_sig->SetTitleFont(42, "XYZ");
       h_sig->SetTitleSize(0.055, "XYZ"); 
       h_sig->SetTitleOffset(1.,"Y");
       h_sig->GetYaxis()->SetTitle("#frac{S}{sqrt(S+B)}");
       h_sig->GetYaxis()->SetTitleColor(1);
       h_sig->SetNdivisions(509,"XYZ");
       h_sig->SetMinimum(0.); 
       h_sig->GetXaxis()->SetTitle("MVA score");


  leg = new TLegend(0.25,0.25,0.35,0.45);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  // leg->SetHeader("Significance");
  leg->AddEntry(h_sig,"Significance");
leg->Draw();


//----------------------pad4--------------------------//

pad4->cd();
pad4->SetGrid();

// htitle1 = "hSignalEff";
// htitle2 = "hBkgEff";

////////
//  f3->cd();
// //  TH2F* h0 = (TH2F*)gROOT->FindObject(htitle0);
// //  h0->Sumw2();
//  TH1F* h30 = (TH1F*)gROOT->FindObject(htitle1);
//  h30->Sumw2();
//  f4->cd();
//  TH1F* h31 = (TH1F*)gROOT->FindObject(htitle2);
//  h31->Sumw2();
TH2F*    h_Sdiff = new TH2F("h_Sdiff","h_Sdiff",201,-1,1,201,0,1);
TH2F*    h_Bdiff = new TH2F("h_Bdiff","h_Bdiff",201,-1,1,201,0,1);
for (int k =1; k < size+1 ; k++)
    {
        float  Sbin = h1->GetBinContent(k);
        float  Bbin = h30->GetBinContent(k);
      //   float  SbinVeto = h30->GetBinContent(k);
      //   float  BbinVeto = h31->GetBinContent(k);
        float  index = -1+k*(2./201.);
        float Sdiff = Sbin-Bbin;
      //   float Bdiff = Bbin-BbinVeto;
        h_Sdiff->Fill(index,Sdiff);
      //   h_Bdiff->Fill(index,Bdiff);
    }

       h_Sdiff->Draw("COLZ"); 
       h_Sdiff->SetLineColor(kBlue);
       h_Sdiff->SetLineStyle(1);
       h_Sdiff->SetLineWidth(3);
       h_Sdiff->SetTickLength(0.03, "XYZ");
       h_Sdiff->SetLabelOffset(0.007,"X");
       h_Sdiff->SetLabelOffset(0.007,"Y");
       h_Sdiff->SetLabelSize(0.045, "XYZ");
       h_Sdiff->SetLabelFont(42, "XYZ"); 
       h_Sdiff->SetTitleFont(42, "XYZ");
       h_Sdiff->SetTitleSize(0.055, "XYZ"); 
       h_Sdiff->SetTitleOffset(1.,"Y");
       h_Sdiff->GetYaxis()->SetTitle("efficiency");
       h_Sdiff->GetYaxis()->SetTitleColor(1);
       h_Sdiff->SetNdivisions(509,"XYZ");
       h_Sdiff->SetMinimum(0.); 
       h_Sdiff->GetXaxis()->SetTitle("MVA score");

      //  h_Bdiff->Draw("HISTsame"); 
      //  h_Bdiff->SetLineColor(kRed);
      //  h_Bdiff->SetLineStyle(1);
      //  h_Bdiff->SetLineWidth(3);

         leg = new TLegend(0.6,0.75,0.75,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetHeader("Diif w/ vs w/o  veto");
  leg->AddEntry(h_Sdiff,"Signal diff","L");
  leg->AddEntry(h_Bdiff,"Bkg diff","L");
leg->Draw();
// *****************************************************************************

 c1->Update();
 c1->SaveAs("./"+SaveFile+".png");
}
