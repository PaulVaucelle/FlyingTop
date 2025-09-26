#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PlotCMS.h"

void plot(int method,TString Year, TString CHANNEL, TString Sample)
{
    
int stati=0;
bool fit= 1;
bool logy=0;
 TString Yearcor = Year;
 if (Year == "2016PRE") Yearcor = "2016preVFP";
 if (Year == "2016POST") Yearcor = "2016";
  float hmin = 0.;
  float hmax = 3;	   
    double norm1 = 1.;
    double norm2 = 1.;
TString Channel = CHANNEL;//MUMU ou EMU
TString Dmode = "DM";//DM or EM
if (Channel == "EMU") Dmode = "EM";
// TFile* f1_Data = new TFile("../Signal_"+Year+"/g1_Dataofile_DM_OS_2p4_RPV_"+Yearcor+"_NOM.root");
float rwTT = 1.0;
float scaleMC = 1.0;

 TFile* f1_DY  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f2_DY  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f1_TT  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_TT  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_ST  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_ST  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f3_ST  = new TFile("../../MC_EMU_"+Year+"_03_02_2025/histofile_HT100_EM_OS_2p4_ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f4_ST  = new TFile("../../MC_EMU_"+Year+"_03_02_2025/histofile_HT100_EM_OS_2p4_ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f1_TTV = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8.root");
 TFile* f2_TTV = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8.root");
 TFile* f3_TTV = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8.root");
 TFile* f1_VV  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
 TFile* f2_VV  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 TFile* f3_VV  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
TString hNorma = "hEvents_with_gen_wt";
float SF1 = 1.0;
float SF2 = 1.0;
if (Sample == "TT") {SF1 = 59700*88.5;SF2 = 59700*366.5;}
if (Sample == "DY") {SF1 = 59700*22635;SF2 = 59700*6225.4;}


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
htotMC->Sumw2();
 TH1F* g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);//ok
 TH1F* g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);//ok
 TH1F*  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

TFile* f1_DY_Norm  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
TFile* f2_DY_Norm  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
TFile* f1_TT_Norm  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
TFile* f2_TT_Norm  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");

 TH1F* g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);//ok
 TH1F* g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);//ok
 TH1F* g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);//ok
 TH1F*  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);

TH1F* g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);//ok
 TH1F* g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);//ok
   TH1F* g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);//ok
 TH1F* g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);//ok
 TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 TH1F* g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);//ok
  std::cout<<" path : "<<MCFILE[2]+htitleA<<std::endl;
  std::cout<<" g1_TT : "<<g1_TT<<std::endl;
TH1F* g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);//ok
 TH1F*  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 TH1F* g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);//ok

 TH1F* g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);//ok
 TH1F* g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);//ok
 TH1F*  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);

  c1->cd();
pad1->cd();

 f1_DY->cd();
 if (Sample == "DY")
  {

    f1_DY_Norm->cd("");
    TDirectory* dir1 = f1_DY_Norm->GetDirectory("FlyingTop");
    if (dir1) {
      dir1->cd();  
      TH1D*  e1_DY = (TH1D*)gROOT->FindObject(hNorma);
      e1_DY->Sumw2();
      if  ( e1_DY->GetEntries() > 0 ) norm1 =  e1_DY->GetEntries();

    }
    f2_DY_Norm->cd("");
    TDirectory* dir2 = f2_DY_Norm->GetDirectory("FlyingTop");
    if (dir2) {
      dir2->cd();  
      TH1D*  e1_DY = (TH1D*)gROOT->FindObject(hNorma);
      e1_DY->Sumw2();
      if  ( e1_DY->GetEntries() > 0 ) norm2 =  e1_DY->GetEntries();

    }

    f1_DY->cd();
    g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
      h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
    h_DY->Add(g1_DY, h_DY, SF1/norm1,0);

    f2_DY->cd();
    g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
      h_DY->Add(g2_DY, h_DY, SF2/norm2, 1);
  }
else
  {
       f1_DY->cd();
      g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
      h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
      h_DY->Add(g1_DY, h_DY, 1,0);

      f2_DY->cd();
      g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
      h_DY->Add(g2_DY, h_DY, 1, 1);
  }

 f1_VV->cd();
  g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);
  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);
  h_VV->Add(g2_VV, h_VV,1, 1);

 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);
  h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
  g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);
  h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);
  h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 f1_ST->cd();
  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);
  h_ST->Add(g2_ST, h_ST, 1, 1);

  if (Channel=="EMU")
    {
      f3_ST->cd();
      g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);
      h_ST->Add(g3_ST, h_ST, 1, 1);

      f4_ST->cd();
      g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);
      h_ST->Add(g4_ST, h_ST, 1, 1);
    }

if (Sample == "TT")
  {
        f1_TT_Norm->cd("");
    TDirectory* dir1 = f1_TT_Norm->GetDirectory("FlyingTop");
    if (dir1) {
      dir1->cd();  
      TH1D*  e1_DY = (TH1D*)gROOT->FindObject(hNorma);
      e1_DY->Sumw2();
      if  ( e1_DY->GetEntries() > 0 ) norm1 =  e1_DY->GetEntries();

    }
    f2_TT_Norm->cd("");
    TDirectory* dir2 = f2_TT_Norm->GetDirectory("FlyingTop");
    if (dir2) {
      dir2->cd();  
      TH1D*  e1_DY = (TH1D*)gROOT->FindObject(hNorma);
      e1_DY->Sumw2();
      if  ( e1_DY->GetEntries() > 0 ) norm2 =  e1_DY->GetEntries();

    }


     f1_TT->cd();
    g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
      h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
    h_TT->Add(g1_TT, h_TT, rwTT*SF2/norm1,0);

    f2_TT->cd();
    g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);
    h_TT->Add(g2_TT, h_TT, rwTT*SF2/norm2, 1);
  }
else
  {
    f1_TT->cd();
    g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
      h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
    h_TT->Add(g1_TT, h_TT, rwTT*1,0);

        f2_TT->cd();
        g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);
        h_TT->Add(g2_TT, h_TT, 1,1);
  }

//!! --------------------------
    h_DY->Scale(1*scaleMC);
    h_VV->Scale(1*scaleMC);
    h_TTV->Scale(1*scaleMC);
    h_ST->Scale(1*scaleMC);
    h_TT->Scale(1*scaleMC);
  if (Sample == "DY")
    {
      htotMC->Add(h_DY, htotMC, 1, 0);
    }
  else if (Sample == "TT")
    {
      htotMC->Add(h_TT, htotMC, 1, 0);
    }
  else
    {
      std::cout<<" Sample : "<<Sample<<std::endl;
      htotMC->Add(h_DY, htotMC, 1, 0);
      htotMC->Add(h_VV, htotMC, 1, 1);
      htotMC->Add(h_TTV, htotMC, 1, 1);
      htotMC->Add(h_ST, htotMC, 1, 1);
      htotMC->Add(h_TT, htotMC, 1, 1);
    }

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
//  htotMC->SetMinimum(hmin); 
//  htotMC->SetMaximum(htotMC->GetMaximum()*2.5); 
// htotMC->Scale(SF/htotMC->Integral());
 htotMC->SetTitle(""); 


//  htotMC->SetMinimum(0); 
// //  htotMC->SetMaximum(hmax);

// !! ----------------------------------------------------!!//
// Fit par une droite pour avoir la pente => corrélation

TString extraLeg = " #mu#mu";
if (Channel == "EMU") extraLeg = " e#mu";
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
  if (Channel == "MUMU")
    {
      SaveFile = SaveFile + "_"+Sample+"_MUMU_MC";
    }
  else
    {
      SaveFile = SaveFile + "_"+Sample+"_EMU_MC";
    }
  c1->SaveAs(SaveFile+".pdf");
  c1->SaveAs(SaveFile+".root");
  delete c1;
}

// #include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PlotCMS.h"

// void plot(int method,TString Year, TString CHANNEL)
// {
    
// int stati=0;
// bool fit= 1;
// bool logy=0;
//  TString Yearcor = Year;
//  if (Year == "2016PRE") Yearcor = "2016preVFP";
//  if (Year == "2016POST") Yearcor = "2016";
//   float hmin = 0.;
//   float hmax = 3;	   
// TString Channel = CHANNEL;//MUMU ou EMU
// TString Dmode = "DM";//DM or EM
// if (Channel == "EMU") Dmode = "EM";
// // TFile* f1_Data = new TFile("../Signal_"+Year+"/g1_Dataofile_DM_OS_2p4_RPV_"+Yearcor+"_NOM.root");
// float rwTT = 1.0;
// float scaleMC = 1.0;

//  TFile* f1_DY  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
//  TFile* f2_DY  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
//  TFile* f1_TT  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f2_TT  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f1_ST  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f2_ST  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f3_ST  = new TFile("../../MC_EMU_"+Year+"_03_02_2025/histofile_HT100_EM_OS_2p4_ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f4_ST  = new TFile("../../MC_EMU_"+Year+"_03_02_2025/histofile_HT100_EM_OS_2p4_ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f1_TTV = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8.root");
//  TFile* f2_TTV = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8.root");
//  TFile* f3_TTV = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_TTWW_TuneCP5_13TeV-madgraph-pythia8.root");
//  TFile* f1_VV  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
//  TFile* f2_VV  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
//  TFile* f3_VV  = new TFile("../../MC_"+Channel+"_"+Year+"_03_02_2025/histofile_HT100_"+Dmode+"_OS_2p4_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.root");
 


// TString MCFILE[14] = {
  
//                   "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_",
//                   "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_",
//                   "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_",
//                   "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_",
//                   "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_",
//                   "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_",
//                   "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8_",
//                   "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8_",
//                   "TTWW_TuneCP5_13TeV-madgraph-pythia8_",
//                   "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_",
//                   "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_",
//                   "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_",
//                   "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_",
//                   "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8_"
// };

//  TString ytitle = "Tight/Loose Vtx ratio"; 
//  TString HeaderCMS = "CMS";

//  if (Year == "2016") HeaderCMS = "2016                                        36.3 fb^{-1} (13 TeV)";
//  if (Year == "2017") HeaderCMS = "2017                                        41.5 fb^{-1} (13 TeV)";
//  if (Year == "2018") HeaderCMS = "2018                                        59.8 fb^{-1} (13 TeV)";

//     TString htitleA = "hData_VtxQualityTight_Hemi1pt_2Vtx";
//     TString htitleB = "hData_VtxQualityLoose_Hemi1pt_2Vtx";
//     int nbin = 8; 
//     float xmin = 0;
//     float xmax =  8;
//     TString HeaderA = "A";
//     TString HeaderNVtx = "k Vtx";
//     TString SaveFile = "QualityRatio1D";
//     TString xtitle = "Hemi_{pt} [GeV]";

// int Method = method;
//   if (Method == 0)
//   {
//     htitleA = "StepEffi_";
//     HeaderA = "";
//     HeaderNVtx = "Cutflow"; 
//     SaveFile = "StepEffi_";
//   }

// TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
// c1->SetFillColor(10);
// c1->SetFillStyle(4000);
// c1->SetBorderSize(2);

// TPad* pad1 = new TPad("pad1","This is pad1",0.04,0.01,0.96,0.99,21);
// pad1->SetFillColor(0);
// pad1->SetBorderMode(0);
// pad1->SetFrameFillColor(10);
// pad1->Draw();
// pad1->SetLogy(logy);
//    pad1->SetTopMargin(0.1);
//    pad1->SetBottomMargin(0.15);
//    pad1->SetRightMargin(0.05);
//    pad1->SetLeftMargin(0.15);



//    gStyle->SetOptDate(0);
// gStyle->SetStatColor(0);
// gStyle->SetTitleFont(62);
// gStyle->SetTitleColor(1);
// gStyle->SetTitleTextColor(1);
// gStyle->SetTitleFillColor(10);
// gStyle->SetTitleFontSize(0.05);
// gStyle->SetTitleW(0.4);
// gStyle->SetTitleH(0.09);
// gStyle->SetOptStat(stati);
// gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
// gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);
// gROOT->SetBatch(kTRUE);

// //  g1_Data->Sumw2();
// TH1F* htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);

//  TH1F* g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);//ok
//  TH1F* g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);//ok
//  TH1F*  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);

//  TH1F* g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);//ok
//  TH1F* g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);//ok
//  TH1F* g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);//ok
//  TH1F*  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);

// TH1F* g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);//ok
//  TH1F* g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);//ok
//    TH1F* g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);//ok
//  TH1F* g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);//ok
//  TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

//  TH1F* g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);//ok
// TH1F* g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);//ok
//  TH1F*  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

//  TH1F* g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);//ok
//  TH1F* g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);//ok
//  TH1F* g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);//ok
//  TH1F*  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);

//   c1->cd();
// pad1->cd();

//  f1_DY->cd();
//  g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
//   h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
//  h_DY->Add(g1_DY, h_DY, 1,0);

//  f2_DY->cd();
//  g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
//   h_DY->Add(g2_DY, h_DY, 1, 1);

//  f1_VV->cd();
//   g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);
//   h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
//  h_VV->Add(g1_VV, h_VV, 1,0);

//  f2_VV->cd();
//  g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);
//   h_VV->Add(g2_VV, h_VV,1, 1);

//  f3_VV->cd();
//  g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);
//   h_VV->Add(g3_VV, h_VV, 1, 1);

//  f1_TTV->cd();
//  g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);
//   h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
//  h_TTV->Add(g1_TTV, h_TTV, 1,0);

//  f2_TTV->cd();
//   g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);
//   h_TTV->Add(g2_TTV, h_TTV, 1, 1);

//  f3_TTV->cd();
//  g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);
//   h_TTV->Add(g3_TTV, h_TTV, 1, 1);

//  f1_ST->cd();
//   h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
//  g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);
//  h_ST->Add(g1_ST, h_ST, 1,0);

//  f2_ST->cd();
//  g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);
//   h_ST->Add(g2_ST, h_ST, 1, 1);

//   if (Channel=="EMU")
//     {
//       f3_ST->cd();
//       g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);
//       h_ST->Add(g3_ST, h_ST, 1, 1);

//       f4_ST->cd();
//       g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);
//       h_ST->Add(g4_ST, h_ST, 1, 1);
//     }

//  f1_TT->cd();
//  g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
//   h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
//  h_TT->Add(g1_TT, h_TT, rwTT*1,0);

//     f2_TT->cd();
//     g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);
//     h_TT->Add(g2_TT, h_TT, 1,1);
// //!! --------------------------
//     h_DY->Scale(1*scaleMC);
//     h_VV->Scale(1*scaleMC);
//     h_TTV->Scale(1*scaleMC);
//     h_ST->Scale(1*scaleMC);
//     h_TT->Scale(1*scaleMC);



// // Daniel
// //  h_ST->Add(h_ST, h_VV, 1, 1);
// //  h_TTV->Add(h_TTV, h_ST, 1, 1);
// //  h_TT->Add(h_TT, h_TTV, 1, 1);
// //  htotMC->Add(h_DY, htotMC, 1, 0);
// // htotMC->Add(h_TT, htotMC, 1, 1);

// //   htotMC->SetFillStyle(1001);
// // //  htotMC->SetFillColorAlpha(kGreen+1, 1);
// //  htotMC->SetLineColor(kBlack);
// //  htotMC->Draw("P");
// //  htotMC->SetMarkerStyle(20);
// //  htotMC->SetMarkerSize(1.5);
// //  htotMC->SetMarkerColor(kBlack);
// //  htotMC->SetLineColor(kBlack);
// //  htotMC->SetLineWidth(1);
// //  htotMC->SetTickLength(0.03, "YZ");
// //  htotMC->SetTickLength(0.03,"X");
// //  htotMC->SetLabelOffset(0.01,"X");
// //  htotMC->SetLabelOffset(0.007,"Y");
// //  htotMC->SetLabelSize(0.035, "XYZ");
// //  htotMC->SetLabelFont(42, "XYZ"); 
// //  htotMC->SetTitleSize(0.045, "XYZ"); 
// //  htotMC->SetTitleFont(42, "XYZ");
// //  htotMC->SetTitleOffset(1.2,"X"); 
// //  htotMC->SetTitleOffset(1.5,"Y");
// //  htotMC->GetXaxis()->SetTitle(xtitle);
// //  htotMC->GetXaxis()->SetTitleColor(1);
// //  htotMC->GetXaxis()->SetRangeUser(0,260);
// //  htotMC->GetYaxis()->SetTitle(ytitle);
// //  htotMC->GetYaxis()->SetTitleColor(1);
// //  htotMC->SetNdivisions(509,"XYZ");
// //  htotMC->SetMinimum(hmin); 
// //  htotMC->SetMaximum(htotMC->GetMaximum()*2.5); 
// //  htotMC->SetTitle(""); 

//  h_TT->SetLineColor(kBlack);
//  h_TT->Draw("P");
//  h_TT->SetMarkerStyle(20);
//  h_TT->SetMarkerSize(1.5);
//  h_TT->SetMarkerColor(kBlack);
//  h_TT->SetLineColor(kBlack);
//  h_TT->SetLineWidth(1);
//  h_TT->SetTickLength(0.03, "YZ");
//  h_TT->SetTickLength(0.03,"X");
//  h_TT->SetLabelOffset(0.01,"X");
//  h_TT->SetLabelOffset(0.007,"Y");
//  h_TT->SetLabelSize(0.035, "XYZ");
//  h_TT->SetLabelFont(42, "XYZ"); 
//  h_TT->SetTitleSize(0.045, "XYZ"); 
//  h_TT->SetTitleFont(42, "XYZ");
//  h_TT->SetTitleOffset(1.2,"X"); 
//  h_TT->SetTitleOffset(1.5,"Y");
//  h_TT->GetXaxis()->SetTitle(xtitle);
//  h_TT->GetXaxis()->SetTitleColor(1);
//  h_TT->GetXaxis()->SetRangeUser(0,260);
//  h_TT->GetYaxis()->SetTitle(ytitle);
//  h_TT->GetYaxis()->SetTitleColor(1);
//  h_TT->SetNdivisions(509,"XYZ");
//  h_TT->SetMinimum(hmin); 
//  h_TT->SetMaximum(h_TT->GetMaximum()*2.5); 
//  h_TT->SetTitle(""); 


 
// //  htotMC->SetMinimum(0); 
// // //  htotMC->SetMaximum(hmax);

// // !! ----------------------------------------------------!!//
// // Fit par une droite pour avoir la pente => corrélation

// TString extraLeg = " #mu#mu";
// if (Channel == "EMU") extraLeg = " e#mu";
//   TLegend* leg = new TLegend(0.7,0.75,0.8,0.85);
//   leg->SetBorderSize(0);
//   leg->SetFillColor(kWhite);
//   leg->SetTextFont(42);
//   leg->SetTextSize(0.035);
//   leg->SetMargin(0.2);
//   leg->AddEntry(htotMC,"MC "+extraLeg,"PE1");
//   leg->Draw();

//  leg = new TLegend(0.7,0.70,0.8,0.75);
//   leg->SetBorderSize(0);
//   leg->SetFillColor(kWhite);
//   leg->SetTextFont(42);
//   leg->SetTextSize(0.035);
//   leg->SetMargin(0.2);
//   leg->SetHeader(HeaderNVtx);
//   leg->Draw();

//    leg = new TLegend(0.7,0.65,0.8,0.7);
//   leg->SetBorderSize(0);
//   leg->SetFillColor(kWhite);
//   leg->SetTextFont(42);
//   leg->SetTextSize(0.035);
//   leg->SetMargin(0.2);
//   leg->SetHeader(HeaderA);
//   leg->Draw();

//   PlotCMSv2(pad1,Year,false);
//   if (Channel == "MUMU")
//     {
//       SaveFile = SaveFile + "_MUMU_MC";
//     }
//   else
//     {
//       SaveFile = SaveFile + "_EMU_MC";
//     }
//   c1->SaveAs(SaveFile+".pdf");
//   c1->SaveAs(SaveFile+".root");
//   delete c1;
// }