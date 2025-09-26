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
float rwTT = 1.0;
float scaleMC = 1.0;
TString ProdMC = "MC_EMU_2018_03_02_2025";
TString suffixDATA = "_BDT100";
TString suffixMC = "";

 if (Year == "2016PRE") 
  {
    Yearcor = "2016preVFP";
    scaleMC = 0.5372;
    ProdMC = "MC_EMU_2016PRE_30_03_2025";
    if (Channel == "MUMU") ProdMC = "MC_MUMU_2016PRE_30_03_2025";

  }
 if (Year == "2016POST") 
  {
    Yearcor = "2016";
    scaleMC = 0.4628*0.964;
    ProdMC = "MC_EMU_2016POST_30_03_2025";
    if (Channel == "MUMU") {ProdMC = "MC_MUMU_2016POST_30_03_2025";scaleMC = 0.4628*1;}
  }
 if (Year == "2017") 
  {
    Yearcor = "2017";
    scaleMC = 1.;
    ProdMC = "MC_EMU_2017_30_03_2025";
    if (Channel == "MUMU") ProdMC = "MC_MUMU_2017_30_03_2025";

  }
 if (Year == "2018") 
  {
    Yearcor = "2018";
    scaleMC = 0.948;
    ProdMC = "MC_EMU_2018_03_02_2025";
    if (Channel == "MUMU") ProdMC = "MC_MUMU_2018_03_02_2025";
  }
  

 TFile* f1_DY  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8"+suffixMC+".root");
 TFile* f2_DY  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"+suffixMC+".root");
 TFile* f1_TT  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"+suffixMC+".root");
 TFile* f2_TT  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"+suffixMC+".root");
 TFile* f1_ST  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"+suffixMC+".root");
 TFile* f2_ST  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"+suffixMC+".root");
 TFile* f3_ST  = new TFile("../../"+ProdMC+"/histofile_HT100_EM_OS_2p4_ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8"+suffixMC+".root");// only emu
 TFile* f4_ST  = new TFile("../../"+ProdMC+"/histofile_HT100_EM_OS_2p4_ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8"+suffixMC+".root");//only emu
 TFile* f1_TTV = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8"+suffixMC+".root");
 TFile* f2_TTV = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8"+suffixMC+".root");
 TFile* f3_TTV = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8"+suffixMC+".root");
 TFile* f1_VV  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"+suffixMC+".root");
 TFile* f2_VV  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8"+suffixMC+".root");
 TFile* f3_VV  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8"+suffixMC+".root");
 


TString MCFILE[14] = {
  
                  "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_",
                  "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_",
                  "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_",
                  "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_",
                  "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_",
                  "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_",
                  "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8_",
                  "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8_",
                  "TTWW_TuneCP5_13TeV-madgraph-pythia8_",
                  "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_",
                  "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_",
                  "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_",
                  "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_",
                  "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_"
};

 TString ytitle = "Tight/Loose Vtx ratio"; 
 TString HeaderCMS = "CMS";

 if (Year == "2016") HeaderCMS = "2016                                        36.3 fb^{-1} (13 TeV)";
 if (Year == "2017") HeaderCMS = "2017                                        41.5 fb^{-1} (13 TeV)";
 if (Year == "2018") HeaderCMS = "2018                                        59.8 fb^{-1} (13 TeV)";

    TString htitleA = "hData_VtxQualityTight_Hemi1pt_2Vtx";
    TString htitleB = "hData_VtxQualityLoose_Hemi1pt_2Vtx";
TString htitleC = "hData_VtxQualityTight_Hemi1pt_2Vtx";
    TString htitleD = "hData_VtxQualityLoose_Hemi1pt_2Vtx";

    int nbin = 50; 
    float xmin = 0;
    float xmax =  1000;
    TString HeaderA = "A";
    TString HeaderNVtx = "k Vtx";
    TString SaveFile = "QualityRatio1D";
    TString xtitle = "Hemi_{pt} [GeV]";

int Method = method;
  if (Method == 0)
    {
      htitleA = "hData_VtxQualityTight_Hemipt_2Vtx";
      htitleB = "hData_VtxQualityLoose_Hemipt_2Vtx";
      HeaderA = "";
      HeaderNVtx = "2 Vtx"; 
      SaveFile = "VtxQuality_Hemipt_2Vtx";
      xtitle = " pt_{Hemi} [GeV]";
    }



TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

TPad* pad1 = new TPad("pad1","This is pad1",0.04,0.3,0.96,0.99,21);
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.1);
   pad1->SetBottomMargin(0.15);
   pad1->SetRightMargin(0.05);
   pad1->SetLeftMargin(0.15);

TPad* rap1 = new TPad("rap1","This is rap1",0.04,0.01,0.96,0.3,21);
rap1->SetFillColor(0);
rap1->SetBorderMode(0);
rap1->SetFrameFillColor(10);
rap1->Draw();
rap1->SetLogy(0);
   rap1->SetTopMargin(0.1);
   rap1->SetBottomMargin(0.15);
   rap1->SetRightMargin(0.05);
   rap1->SetLeftMargin(0.15);


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
TH1F* htotMC_B  = new TH1F("htotMC_B","",nbin,xmin,xmax);

 TH1F* g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);//ok
 TH1F* g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);//ok
 TH1F*  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

  TH1F* g1_DY_B = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);//ok
  TH1F* g2_DY_B = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleB);//ok
  TH1F*  h_DY_B = new TH1F("h_DY_B","",nbin,xmin,xmax);

 TH1F* g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);//ok
 TH1F* g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);//ok
 TH1F* g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);//ok
 TH1F*  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);

  TH1F* g1_VV_B = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleB);//ok
  TH1F* g2_VV_B = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleB);//ok
  TH1F* g3_VV_B = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleB);//ok
  TH1F*  h_VV_B = new TH1F("h_VV_B","",nbin,xmin,xmax);

TH1F* g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);//ok
 TH1F* g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);//ok
   TH1F* g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);//ok
 TH1F* g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);//ok
 TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

  TH1F* g1_ST_B = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleB);//ok
  TH1F* g2_ST_B = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleB);//ok
  TH1F* g3_ST_B = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleB);//ok
  TH1F* g4_ST_B = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleB);//ok
  TH1F*  h_ST_B = new TH1F("h_ST_B","",nbin,xmin,xmax);

 TH1F* g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);//ok
TH1F* g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);//ok
 TH1F*  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

  TH1F* g1_TT_B = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);//ok
  TH1F* g2_TT_B = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleB);//ok
  TH1F*  h_TT_B = new TH1F("h_TT_B","",nbin,xmin,xmax);

 TH1F* g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);//ok
 TH1F* g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);//ok
 TH1F* g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);//ok
 TH1F*  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);

  TH1F* g1_TTV_B = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleB);//ok
  TH1F* g2_TTV_B = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleB);//ok
  TH1F* g3_TTV_B = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleB);//ok
  TH1F*  h_TTV_B = new TH1F("h_TTV_B","",nbin,xmin,xmax);

  c1->cd();
pad1->cd();

 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 g1_DY_B = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);
  h_DY_B = new TH1F("h_DY_B","",nbin,xmin,xmax);
 h_DY_B->Add(g1_DY_B, h_DY_B, 1,0);


 f2_DY->cd();
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
  h_DY->Add(g2_DY, h_DY, 1, 1);

 g2_DY_B = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleB);
  h_DY_B->Add(g2_DY_B, h_DY_B, 1, 1);

 f1_VV->cd();
  g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);
  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

  g1_VV_B = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleB);
  h_VV_B = new TH1F("h_VV_B","",nbin,xmin,xmax);
 h_VV_B->Add(g1_VV_B, h_VV_B, 1,0);


 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);
  h_VV->Add(g2_VV, h_VV,1, 1);

 g2_VV_B = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleB);
  h_VV_B->Add(g2_VV_B, h_VV_B,1, 1);


 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);
  h_VV->Add(g3_VV, h_VV, 1, 1);

 g3_VV_B = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleB);
  h_VV_B->Add(g3_VV_B, h_VV_B, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 g1_TTV_B = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleB);
  h_TTV_B = new TH1F("h_TTV_B","",nbin,xmin,xmax);
 h_TTV_B->Add(g1_TTV_B, h_TTV_B, 1,0);


 f2_TTV->cd();
  g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);
  h_TTV->Add(g2_TTV, h_TTV, 1, 1);

  g2_TTV_B = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleB);
  h_TTV_B->Add(g2_TTV_B, h_TTV_B, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);
  h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 g3_TTV_B = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleB);
  h_TTV_B->Add(g3_TTV_B, h_TTV_B, 1, 1);

 f1_ST->cd();
  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);
 h_ST->Add(g1_ST, h_ST, 1,0);

  h_ST_B = new TH1F("h_ST_B","",nbin,xmin,xmax);
 g1_ST_B = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleB);
 h_ST_B->Add(g1_ST_B, h_ST_B, 1,0);

 f2_ST->cd();
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);
  h_ST->Add(g2_ST, h_ST, 1, 1);

 g2_ST_B = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleB);
  h_ST_B->Add(g2_ST_B, h_ST_B, 1, 1);

  if (Channel=="EMU")
    {
      f3_ST->cd();
      g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);
      h_ST->Add(g3_ST, h_ST, 1, 1);

      g3_ST_B = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleB);
      h_ST_B->Add(g3_ST_B, h_ST_B, 1, 1);

      f4_ST->cd();
      g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);
      h_ST->Add(g4_ST, h_ST, 1, 1);
      g4_ST_B = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleB);
      h_ST_B->Add(g4_ST_B, h_ST_B, 1, 1);
    }

 f1_TT->cd();
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, rwTT*1,0);

  g1_TT_B = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);
  h_TT_B = new TH1F("h_TT_B","",nbin,xmin,xmax);
  h_TT_B->Add(g1_TT_B, h_TT_B, rwTT*1,0);

    f2_TT->cd();
    g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);
    h_TT->Add(g2_TT, h_TT, 1,1);

    g2_TT_B = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleB);
    h_TT_B->Add(g2_TT_B, h_TT_B, 1,1);
//!! --------------------------
    h_DY->Scale(1*scaleMC);
    h_VV->Scale(1*scaleMC);
    h_TTV->Scale(1*scaleMC);
    h_ST->Scale(1*scaleMC);
    h_TT->Scale(1*scaleMC);

    h_DY_B->Scale(1*scaleMC);
    h_VV_B->Scale(1*scaleMC);
    h_TTV_B->Scale(1*scaleMC);
    h_ST_B->Scale(1*scaleMC);
    h_TT_B->Scale(1*scaleMC);

// Daniel
 h_ST->Add(h_ST, h_VV, 1, 1);
 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TT->Add(h_TT, h_TTV, 1, 1);
 htotMC->Add(h_DY, htotMC, 1, 0);
htotMC->Add(h_TT, htotMC, 1, 1);

 h_ST_B->Add(h_ST_B, h_VV_B, 1, 1);
 h_TTV_B->Add(h_TTV_B, h_ST_B, 1, 1);
 h_TT_B->Add(h_TT_B, h_TTV_B, 1, 1);
 htotMC_B->Add(h_DY_B, htotMC_B, 1, 0);
htotMC_B->Add(h_TT_B, htotMC_B, 1, 1);


  htotMC->Divide(htotMC,htotMC_B,1,1);

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
 htotMC->SetMinimum(hmin); 
 htotMC->SetMaximum(htotMC->GetMaximum()*2.5); 
 htotMC->SetTitle(""); 

 
 htotMC->SetMinimum(0); 
//  htotMC->SetMaximum(hmax);


TF1 *CONST = new TF1("CONST", "[0]",100 , htotMC->GetXaxis()->GetXmax());//htotMC->GetXaxis()->GetXmin()
CONST->SetParameter(0,0.04);
htotMC->Fit(CONST, "R");


double c = CONST->GetParameter(0);
double cerr = CONST->GetParError(0);



// !! --------------------------- !!//

TString extraLeg = " #mu#mu";
if (Channel == "EMU") extraLeg = " e#mu";
  TLegend* leg = new TLegend(0.7,0.65,0.8,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(htotMC,"MC "+extraLeg,"PE1");
  // leg->AddEntry(Droite,"Linear Fit","L");
  // leg->AddEntry(DroiteUp,"FitUp","L");
  // leg->AddEntry(DroiteDown,"FitDown","L");
  leg->AddEntry(CONST,"Const Fit","L");
leg->AddEntry(CONSTUP,"Up Fit","L");
  leg->AddEntry(CONSTDOWN,"Down Fit","L");
  leg->Draw();


  if (Method >= 8)
    {
      leg = new TLegend(0.5,0.55,0.85,0.65);
      leg->SetBorderSize(0);
      leg->SetFillColor(kWhite);
      leg->SetTextFont(42);
      leg->SetTextSize(0.035);
      leg->SetMargin(0.2);
      // double a = Droite->GetParameter(0);
      // double b = Droite->GetParameter(1);
      // double aerr = Droite->GetParError(0);
      // double berr = Droite->GetParError(1);
      // leg->SetHeader("Slope = "+TString::Format("%.2e #pm %.2e",a,aerr));
      double c = CONST->GetParameter(0);
      double cerr = CONST->GetParError(0);
      leg->AddEntry(CONST,TString::Format("Const = %.2e #pm %.2e",c,cerr),"L");
      leg->Draw();
    }

  PlotCMSv2(pad1,Year,false);

    rap1->cd();
  TH1F* h_ratio = (TH1F*)htotMC->Clone("h_ratio");
  TH1F* h_ratioUp = (TH1F*)htotMC->Clone("h_ratioUp");

  h_ratio->Reset();
  h_ratioUp->Reset();

  // h_ratio->SetTitle("Ratio: data / fit extrapolation");

  int nBins = htotMC->GetNbinsX();
  for (int i = 1; i <= nBins; ++i) {
    double x = htotMC->GetBinCenter(i);
    if (x >= 30 && x <= 100) {

      
      // double y_data = htotMC->GetBinContent(i);
      // double y_fit = Droite->Eval(x);
      // double y_fit_up = DroiteUp->Eval(x);
      // double y_fit_down = DroiteDown->Eval(x);
      // if (y_fit != 0) {
      //   h_ratio->SetBinContent(i, y_data / y_fit);
      //   h_ratio->SetBinError(i, htotMC->GetBinError(i) / y_fit); // propagation simple
      // }
      // if (y_fit_up != 0) {
      //   h_ratioUp->SetBinContent(i, y_data / y_fit_up);
      //   h_ratioUp->SetBinError(i, htotMC->GetBinError(i) / y_fit_up); // propagation simple
      // }
      // if (y_fit_down != 0) {
      //   h_ratioDown->SetBinContent(i, y_data / y_fit_down);
      //   h_ratioDown->SetBinError(i, htotMC->GetBinError(i) / y_fit_down); // propagation simple
      // }


            double y_data = htotMC->GetBinContent(i);
      double y_fit = CONST->Eval(x);
      double y_fit_up = CONSTUP->Eval(x);
      double y_fit_down = CONSTDOWN->Eval(x);
      if (y_fit != 0) {
        h_ratio->SetBinContent(i, y_data / y_fit);
        h_ratio->SetBinError(i, htotMC->GetBinError(i) / y_fit); // propagation simple
      }
      if (y_fit_up != 0) {
        h_ratioUp->SetBinContent(i, y_data / y_fit_up);
        h_ratioUp->SetBinError(i, htotMC->GetBinError(i) / y_fit_up); // propagation simple
      }
      if (y_fit_down != 0) {
        h_ratioDown->SetBinContent(i, y_data / y_fit_down);
        h_ratioDown->SetBinError(i, htotMC->GetBinError(i) / y_fit_down); // propagation simple
      }
    }
   else if (x >= 100)
    {
      double y_data = htotMC->GetBinContent(i);
      double y_fit = CONST->Eval(x);
      double y_fit_up = CONSTUP->Eval(x);
      double y_fit_down = CONSTDOWN->Eval(x);
      if (y_fit != 0) {
        h_ratio->SetBinContent(i, y_data / y_fit);
        h_ratio->SetBinError(i, htotMC->GetBinError(i) / y_fit); // propagation simple
      }
      if (y_fit_up != 0) {
        h_ratioUp->SetBinContent(i, y_data / y_fit_up);
        h_ratioUp->SetBinError(i, htotMC->GetBinError(i) / y_fit_up); // propagation simple
      }
      if (y_fit_down != 0) {
        h_ratioDown->SetBinContent(i, y_data / y_fit_down);
        h_ratioDown->SetBinError(i, htotMC->GetBinError(i) / y_fit_down); // propagation simple
      }
    }
  }
  h_ratio->SetFillStyle(1001);
//  h_ratio->SetFillColorAlpha(kGreen+1, 1);
 h_ratio->SetLineColor(kRed);
 h_ratio->Draw("PE1");
 h_ratio->SetMarkerStyle(20);
 h_ratio->SetMarkerSize(1.5);
 h_ratio->SetMarkerColor(kRed);
 h_ratio->SetLineColor(kRed);
 h_ratio->SetLineWidth(1);
 h_ratio->SetTickLength(0.03, "YZ");
 h_ratio->SetTickLength(0.03,"X");
 h_ratio->SetLabelOffset(0.01,"X");
 h_ratio->SetLabelOffset(0.007,"Y");
 h_ratio->SetLabelSize(0.06, "XYZ");
 h_ratio->SetLabelFont(42, "XYZ"); 
 h_ratio->SetTitleSize(0.085, "XYZ"); 
 h_ratio->SetTitleFont(42, "XYZ");
 h_ratio->SetTitleOffset(1.2,"X"); 
 h_ratio->SetTitleOffset(0.8,"Y");
 h_ratio->GetXaxis()->SetTitle(xtitle);
 h_ratio->GetXaxis()->SetTitleColor(1);
 h_ratio->GetXaxis()->SetRangeUser(0,260);
//  h_ratio->GetYaxis()->SetTitle(ytitle);
 h_ratio->GetYaxis()->SetTitleColor(1);
 h_ratio->SetNdivisions(509,"XYZ");
 h_ratio->SetMinimum(0.); 
 h_ratio->SetMaximum(2.); 
 h_ratio->SetTitle(""); 

  h_ratio->GetYaxis()->SetTitle("MC / Fit");

  h_ratioUp->SetFillStyle(1001);
  h_ratioUp->SetFillColorAlpha(kBlue+3, 1);
  h_ratioUp->SetLineColor(kBlue+3);
  h_ratioUp->SetLineWidth(1);
  h_ratioUp->Draw("PE1 same");



  leg = new TLegend(0.7,0.65,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.055);
  leg->SetMargin(0.2);

  leg->AddEntry(h_ratio,"Fit","L");

  leg->Draw();

  c1->SaveAs("Compare_Corr.pdf");
//   rap1->SaveAs(SaveFile+"_"+Year+".root");
  delete c1;
}