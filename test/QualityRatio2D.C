#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PlotCMS.h"

void plot(int method,TString Year, TString CHANNEL)
{
    
int stati=0;
bool fit= 1;
bool logy=0;
 TString Yearcor = Year;
 if (Year == "2016PRE") Yearcor = "2016preVFP";
 if (Year == "2016POST") Yearcor = "2016";
  float hmin = 0.;
  float hmax = 3;	   
TString Channel = CHANNEL;//MUMU ou EMU
TString Dmode = "DM";//DM or EM
if (Channel == "EMU") Dmode = "EM";
// TFile* f1_Data = new TFile("../Signal_"+Year+"/g1_Dataofile_DM_OS_2p4_RPV_"+Yearcor+"_NOM.root");
TString Prod = "DATA_EMU_2018_03_02_2025";
TString ProdSignal = "Signal_2018_L1";
TString suffixDATA = "_BDT100";
TString SampleDATA = "MuonEG-Run2018-UL2018_MiniAODv2_GT36-v1";
TString SampleDATAExtra  = "MuonEG-Run2018-UL2018_MiniAODv2_GT36-v1"+suffixDATA;

 if (Year == "2016PRE") 
  {
    Yearcor = "2016preVFP";
    Prod = "DATA_EMU_2016PRE_30_03_2025";
    SampleDATA = "MuonEG_Run2016-HIPM_UL2016_MiniAODv2";
    SampleDATAExtra = "MuonEG_Run2016-HIPM_UL2016_MiniAODv2"+suffixDATA;
  }
  
 if (Year == "2016POST") 
  {
    Yearcor = "2016";
    Prod = "DATA_EMU_2016POST_30_03_2025";
    SampleDATA = "MuonEG_Run2016-UL2016_MiniAODv2";
    SampleDATAExtra = "MuonEG_Run2016-UL2016_MiniAODv2"+suffixDATA;

  }
 if (Year == "2017") 
  {
    Yearcor = "2017";
    Prod = "DATA_EMU_2017_30_03_2025";
    SampleDATA = "MuonEG_Run2017-UL2017_MiniAODv2";
    SampleDATAExtra = "MuonEG_Run2017-UL2017_MiniAODv2"+suffixDATA;

  }
 if (Year == "2018") 
  {
    Yearcor = "2018";
    Prod = "DATA_EMU_2018_03_02_2025";
  }
  

TFile* f1_Data = new TFile("../../"+Prod+"/histofile_HT100_EM_OS_2p4_"+SampleDATAExtra+".root");


TString FILE[1] = { SampleDATA+"_"};

 TString ytitle = "Hemi p_{t} [GeV]"; 
 TString HeaderCMS = "CMS";

 if (Year == "2016") HeaderCMS = "2016                                        36.3 fb^{-1} (13 TeV)";
 if (Year == "2017") HeaderCMS = "2017                                        41.5 fb^{-1} (13 TeV)";
 if (Year == "2018") HeaderCMS = "2018                                        59.8 fb^{-1} (13 TeV)";

    TString htitleA = "hData_VtxQualityTight_Hemi1pt_2Vtx";
    TString htitleB = "hData_VtxQualityLoose_Hemi1pt_2Vtx";
    int nbin = 50; 
    float xmin = 0;
    float xmax =  1000;
    TString HeaderA = "A";
    TString HeaderNVtx = "k Vtx";
    TString SaveFile = "QualityRatio1D";
    TString xtitle = "Vtx BDT Score";

int Method = method;
  if (Method == 0)
  {
    htitleA = "hData_Hemipt_VtxBDT_2Vtx";
    HeaderA = "Hemi_{1}";
    HeaderNVtx = "2 Vtx"; 
    SaveFile = "Hemipt_VtxBDT_2Vtx";
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
   pad1->SetRightMargin(0.105);
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

  f1_Data->cd();
 TH2F* g1_Data = (TH2F*)gROOT->FindObject(FILE[0]+htitleA);//ok


//  g1_Data->Sumw2();

  c1->cd();
pad1->cd();
  // g1_Data->Divide(g1_Data,g2_Data,1,1);

  g1_Data->SetFillStyle(1001);
//  g1_Data->SetFillColorAlpha(kGreen+1, 1);
 g1_Data->SetLineColor(kBlack);
 g1_Data->Draw("COLZ");
 g1_Data->SetMarkerStyle(20);
 g1_Data->SetMarkerSize(1.5);
 g1_Data->SetMarkerColor(kBlack);
 g1_Data->SetLineColor(kBlack);
 g1_Data->SetLineWidth(1);
 g1_Data->SetTickLength(0.03, "YZ");
 g1_Data->SetTickLength(0.03,"X");
 g1_Data->SetLabelOffset(0.01,"X");
 g1_Data->SetLabelOffset(0.007,"Y");
 g1_Data->SetLabelSize(0.035, "XYZ");
 g1_Data->SetLabelFont(42, "XYZ"); 
 g1_Data->SetTitleSize(0.045, "XYZ"); 
 g1_Data->SetTitleFont(42, "XYZ");
 g1_Data->SetTitleOffset(1.2,"X"); 
 g1_Data->SetTitleOffset(1.5,"Y");
 g1_Data->GetXaxis()->SetTitle(xtitle);
 g1_Data->GetXaxis()->SetTitleColor(1);
//  g1_Data->GetXaxis()->SetRangeUser(0,260);
 g1_Data->GetYaxis()->SetTitle(ytitle);
 g1_Data->GetYaxis()->SetTitleColor(1);
 g1_Data->SetNdivisions(509,"XYZ");
//  g1_Data->SetMinimum(hmin); 
//  g1_Data->SetMaximum(g1_Data->GetMaximum()*2.5); 
 g1_Data->SetTitle(""); 

 if(Method == 1)
  {

  }
//  g1_Data->SetMinimum(0); 
//  g1_Data->SetMaximum(hmax);

  PlotCMSv4(pad1,Year,true);

   
  if (Channel == "MUMU")
    {
      SaveFile = SaveFile + "_MUMU_DATA";
    }
  else
    {
      SaveFile = SaveFile + "_EMU_DATA";
    }
  c1->SaveAs(SaveFile+"_"+Year+".pdf");
  pad1->SaveAs(SaveFile+"_"+Year+".root");
  delete c1;
}