#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PlotCMS.h"

void plot(int method,TString Year, TString File)
{
    
int stati=0;
bool fit= 1;
bool logy=0;
 TString Yearcor = Year;
 if (Year == "2016PRE") Yearcor = "2016preVFP";
 if (Year == "2016POST") Yearcor = "2016";
  float hmin = 0.;
  float hmax = 3;	   
TString Dmode = "DM";//DM or EM

// TFile* f1_Data = new TFile("../Signal_"+Year+"/g1_Dataofile_DM_OS_2p4_RPV_"+Yearcor+"_NOM.root");
float rwTT = 1.0;
float scaleMC = 1.0;

 TFile* f1_DY  = new TFile("../../Signal_2018_L1/histofile_HT100_"+Dmode+"_OS_2p4_"+File+".root");

TString MCFILE[1] = {
                 File+"_"
};

float SF = 1.0;
if (File == "RPV_2018_smu200_neu180_ctau100") SF = 534,912;
if (File == "RPV_2018_smu300_neu200_ctau100") SF = 108,654;
if (File == "RPV_2018_smu400_neu200_ctau100") SF = 31,41;
if (File == "RPV_2018_smu500_neu200_ctau100") SF = 11,343;


 TString ytitle = "Tight/Loose Vtx ratio"; 
 TString HeaderCMS = "CMS";

 if (Year == "2016") HeaderCMS = "2016                                        36.3 fb^{-1} (13 TeV)";
 if (Year == "2017") HeaderCMS = "2017                                        41.5 fb^{-1} (13 TeV)";
 if (Year == "2018") HeaderCMS = "2018                                        59.8 fb^{-1} (13 TeV)";

    TString htitleA = "hData_VtxQualityTight_Hemi1pt_2Vtx";
    TString htitleB = "hData_VtxQualityLoose_Hemi1pt_2Vtx";
    int nbin = 8; 
    float xmin = 0;
    float xmax =  8;
    TString HeaderA = "A";
    TString HeaderNVtx = "k Vtx";
    TString SaveFile = "QualityRatio1D";
    TString xtitle = "Hemi_{pt} [GeV]";

int Method = method;
  if (Method == 0)
  {
    htitleA = "StepEffi_";
    HeaderA = "";
    HeaderNVtx = "Cutflow"; 
    SaveFile = "StepEffi_";
  }

TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

TPad* pad1 = new TPad("pad1","This is pad1",0.04,0.01,0.96,0.99,21);
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.1);
   pad1->SetBottomMargin(0.15);
   pad1->SetRightMargin(0.05);
   pad1->SetLeftMargin(0.15);



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

//  g1_Data->Sumw2();
TH1F* htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);

 TH1F* g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);//ok
 TH1F*  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

  c1->cd();
pad1->cd();

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

htotMC->Add(h_DY, htotMC, 1,0);
  htotMC->SetFillStyle(1001);
//  htotMC->SetFillColorAlpha(kGreen+1, 1);
 htotMC->SetLineColor(kBlack);
 htotMC->Draw("PE1");
 htotMC->SetMarkerStyle(20);
 htotMC->SetMarkerSize(1.5);
 htotMC->SetMarkerColor(kBlack);
 htotMC->SetLineColor(kBlack);
 htotMC->SetLineWidth(1);
 htotMC->SetTickLength(0.03, "YZ");
 htotMC->SetTickLength(0.03,"X");
 htotMC->SetLabelOffset(0.01,"X");
 htotMC->SetLabelOffset(0.007,"Y");
 htotMC->SetLabelSize(0.035, "XYZ");
 htotMC->SetLabelFont(42, "XYZ"); 
 htotMC->SetTitleSize(0.045, "XYZ"); 
 htotMC->SetTitleFont(42, "XYZ");
 htotMC->SetTitleOffset(1.2,"X"); 
 htotMC->SetTitleOffset(1.5,"Y");
 htotMC->GetXaxis()->SetTitle(xtitle);
 htotMC->GetXaxis()->SetTitleColor(1);
 htotMC->GetXaxis()->SetRangeUser(0,260);
 htotMC->GetYaxis()->SetTitle(ytitle);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
 htotMC->SetMinimum(0); 
htotMC->Scale(SF/htotMC->Integral());

//  htotMC->SetMaximum(htotMC->GetMaximum()*2.5); 
 htotMC->SetTitle(""); 


//  htotMC->SetMinimum(0); 
// //  htotMC->SetMaximum(hmax);

// !! ----------------------------------------------------!!//
// Fit par une droite pour avoir la pente => corrélation

TString extraLeg = " #mu#mu";

  TLegend* leg = new TLegend(0.7,0.75,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotMC,"MC "+extraLeg,"PE1");
  leg->Draw();

 leg = new TLegend(0.7,0.70,0.8,0.75);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();

   leg = new TLegend(0.7,0.65,0.8,0.7);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->SetHeader(HeaderA);
  leg->Draw();

  PlotCMSv2(pad1,Year,false);

      SaveFile = SaveFile + "_"+File;

  c1->SaveAs(SaveFile+".pdf");
  c1->SaveAs(SaveFile+".root");
  delete c1;
}
