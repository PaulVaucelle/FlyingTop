#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PlotCMS.h"

void plot(int method,TString Year, TString CHANNEL)
{
    
int stati=0;
bool fit= 1;
bool logy=0;
TString Yearcor = "2022";
TString CorBin = "Bin-";
if (Year == "2023A" || Year == "2023B" )Yearcor = "2023";
if (Year == "2024") {Yearcor = "2024";CorBin = "Bin-";}
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
TString DYName1 = "DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8";
TString DYName2 = "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8";


 if (Year == "2022A") 
  {
    Yearcor = "2022A";
    if (Channel == "EMU")
      {
        ProdMC = "MC_EMU_2022_23_04_2025";
        scaleMC = 0.9142*0.3966*1 ;// !! rescale to good lumi * era fraction lumi * jobs succeeded
      }
    else if (Channel == "MUMU")
      {
        ProdMC = "MC_MUMU_2022_23_04_2025";
        scaleMC = 0.9142*0.3966*1;
      }

  }
 if (Year == "2022B") 
  {
    Yearcor = "2022B";
    if (Channel == "EMU")
      {
        ProdMC = "MC_EMU_2022_EFG_23_04_2025";
        scaleMC = 0.9142*0.6004*1 ;
      }
    else if (Channel == "MUMU")
      {
        ProdMC = "MC_MUMU_2022_EFG_23_04_2025";
        scaleMC = 0.9142*0.6004*1;
      }
  }
 if (Year == "2023A") 
  {
    Yearcor = "2023A";
    if (Channel == "EMU")
      {
        ProdMC = "MC_EMU_2023_C_23_04_2025";
        scaleMC =  0.9256*0.653*0.81875;
      }
    else if (Channel == "MUMU")
      {
        ProdMC = "MC_MUMU_2023_C_23_04_2025";
        scaleMC = 0.9256*0.653*0.894;
      }

  }
 if (Year == "2023B") 
  {
    Yearcor = "2023B";
    if (Channel == "EMU")
      {
        ProdMC = "MC_EMU_2023_D_23_04_2025";
        scaleMC =  0.9256*0.347*1  ;
      }
    else if (Channel == "MUMU")
      {
        ProdMC = "MC_MUMU_2023_D_23_04_2025";
        scaleMC = 0.9256*0.347*1;
      }
  }
 if (Year == "2024") 
  {
    Yearcor = "2024";
    if (Channel == "EMU")
      {
        ProdMC = "MC_EMU_2024_23_04_2025";
        scaleMC = 0.99225;
        TString DYName1 = "DYto2L-4Jets_MLL-10to50_TuneCP5_13p6TeV_madgraphMLM-pythia8";
        TString DYName2 = "DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8";
      }
    else if (Channel == "MUMU")
      {
        ProdMC = "MC_MUMU_2024_23_04_2025";
        scaleMC = 0.9424;
        TString DYName1 = "DYto2L-4Jets_MLL-10to50_TuneCP5_13p6TeV_madgraphMLM-pythia8";
        TString DYName2 = "DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8";
      }
  }
  
// !! ---------

 TFile* f1_DY  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"+suffixMC+".root");
 TFile* f2_DY  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"+suffixMC+".root");
if (Year == "2024")
  {
    if (Channel == "EMU") 
      {
        f1_DY  = new TFile("../../MC_EMU_2023_D_23_04_2025_DY/histofile_HT100_"+Dmode+"_OS_2p4_DYto2Tau_"+CorBin+"MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8"+suffixMC+".root");
        f2_DY  = new TFile("../../MC_EMU_2023_D_23_04_2025_DY/histofile_HT100_"+Dmode+"_OS_2p4_DYto2Tau_"+CorBin+"MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8"+suffixMC+".root");
      }
    else if (Channel == "MUMU")
      {
        f1_DY  = new TFile("../../MC_MUMU_2023_D_23_04_2025_DY/histofile_HT100_"+Dmode+"_OS_2p4_DYto2L-4Jets_MLL-10to50_TuneCP5_13p6TeV_madgraphMLM-pythia8.root");
        f2_DY  = new TFile("../../MC_MUMU_2023_D_23_04_2025_DY/histofile_HT100_"+Dmode+"_OS_2p4_DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8.root");

      }
  }
 TFile* f1_TT  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8"+suffixMC+".root");
 TFile* f2_TT  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8"+suffixMC+".root");
 TFile* f1_ST  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8"+suffixMC+".root");
 TFile* f2_ST  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8"+suffixMC+".root");
 TFile* f1_TTV = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8"+suffixMC+".root");
 TFile* f2_TTV = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8"+suffixMC+".root");
 TFile* f3_TTV = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"+suffixMC+".root");
 TFile* f1_VV  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8"+suffixMC+".root");
 TFile* f2_VV  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8"+suffixMC+".root");
 TFile* f3_VV  = new TFile("../../"+ProdMC+"/histofile_HT100_"+Dmode+"_OS_2p4_ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8"+suffixMC+".root");

// !! 


// !! -------------------------------------------------

// !! -------------------------------------------

TString MCFILE[12] = {
  
                    DYName1+"_",
                    DYName2+"_",
                  "TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_",
                  "TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_",

                  "TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_",
                  "TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_",

                
                  "TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8_",
                  "TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8_",
                  "TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_",

                  "WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_",
                  "WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8_",
                  "ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8_",

};

 TString ytitle = "Tight/Loose Vtx ratio"; 

    TString htitleA = "hData_VtxQualityTight_Hemi1pt_2Vtx";
    TString htitleB = "hData_VtxQualityLoose_Hemi1pt_2Vtx";
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
    htitleA = "hData_VtxQualityTight_Hemi1pt_2Vtx";
    htitleB = "hData_VtxQualityLoose_Hemi1pt_2Vtx";
    HeaderA = "Hemi_{1}";
    HeaderNVtx = "2 Vtx"; 
    SaveFile = "Hemi1pt_2Vtx";
  }
    if (Method == 1)
  {
    htitleA = "hData_VtxQualityTight_Hemi2pt_2Vtx";
    htitleB = "hData_VtxQualityLoose_Hemi2pt_2Vtx";
    HeaderA = "Hemi_{2}";
    HeaderNVtx = "2 Vtx"; 
    SaveFile = "Hemi2pt_2Vtx";
  }
  if (Method == 2)
  {
    htitleA = "hData_VtxQualityTight_Hemi1pt_1Vtx";
    htitleB = "hData_VtxQualityLoose_Hemi1pt_1Vtx";
    HeaderA = "Hemi_{1}";
    HeaderNVtx = "1 Vtx"; 
    SaveFile = "Hemi1pt_1Vtx";
  }
    if (Method == 3)
  {
    htitleA = "hData_VtxQualityTight_Hemi2pt_1Vtx";
    htitleB = "hData_VtxQualityLoose_Hemi2pt_1Vtx";
    HeaderA = "Hemi_{2}";
    HeaderNVtx = "1 Vtx"; 
    SaveFile = "Hemi2pt_1Vtx";
  }
    if (Method == 4)
    {
    htitleA = "LT_2Vtx_PromptTT";
    htitleB = "LT_2Vtx_PromptTL";
    HeaderA = "Sum of lepton p_{T}";
    HeaderNVtx = "2 Vtx"; 
    SaveFile = "PromptLT_2Vtx_TTTL";
    xtitle = "L_{T} [GeV]";
  }
    if (Method == 5)
  {
    htitleA = "LT_2Vtx_PromptTL";
    htitleB = "LT_2Vtx_PromptLL";
    HeaderA = "Sum of lepton p_{T}";
    HeaderNVtx = "2 Vtx"; 
    SaveFile = "PromptLT_2Vtx_TLLL";
    xtitle = "L_{T} [GeV]";
}
  if (Method == 6)
    {
      htitleA = "hData_VtxQualityTight_Hemileadingpt_2Vtx";
      htitleB = "hData_VtxQualityLoose_Hemileadingpt_2Vtx";
      HeaderA = "Hemi_{1}";
      HeaderNVtx = "2 Vtx"; 
      SaveFile = "Tight_Hemileadingpt_2Vtx";
      xtitle = " Leading Hemi p_{T} [GeV]";
    }
  if (Method == 7)
    {
      htitleA = "hData_VtxQualityTight_Hemisubleadingpt_2Vtx";
      htitleB = "hData_VtxQualityLoose_Hemisubleadingpt_2Vtx";
      HeaderA = "Hemi_{2}";
      HeaderNVtx = "2 Vtx"; 
      SaveFile = "Tight_Hemisubleadingpt_2Vtx";
      xtitle = " SubLeading Hemi p_{T} [GeV]";
    }
  if (Method == 8)
    {
      htitleA = "hData_VtxQualityTight_HemiAveragept_2Vtx";
      htitleB = "hData_VtxQualityLoose_HemiAveragept_2Vtx";
      HeaderA = "Average Hemi";
      HeaderNVtx = "2 Vtx"; 
      SaveFile = "VtxQuality_HemiAveragept_2Vtx";
      xtitle = " p_{T} [GeV]";
return ;
    }
  if (Method == 9)
    {
      htitleA = "hData_VtxQualityTight_Hemipt_2Vtx";
      htitleB = "hData_VtxQualityLoose_Hemipt_2Vtx";
      HeaderA = "";
      HeaderNVtx = "2 Vtx"; 
      SaveFile = "VtxQuality_Hemipt_2Vtx";
      xtitle = " pt_{Hemi} [GeV]";
    }
if (Method == 10)
  {
      htitleA = "hData_VtxQualityTight_VtxBDT_2Vtx";
      htitleB = "hData_VtxQualityLoose_VtxBDT_2Vtx";
      HeaderA = "";
      HeaderNVtx = "2 Vtx"; 
      SaveFile = "VtxQuality_VtxBT_2Vtx";
      xtitle = " Vtx BDT Score"; 
      return ;
  }

if (Method == 11)
  {
      htitleA = "hData_VtxQualityTight_Hemipt_Focus_2Vtx";
      htitleB = "hData_VtxQualityLoose_Hemipt_Focus_2Vtx";
      HeaderA = "";
      HeaderNVtx = "2 Vtx"; 
      SaveFile = "VtxQuality_Hemipt_Focus_2Vtx";
      xtitle = " Vtx BDT Score"; 

      nbin = 2; 
      xmin = 30;
      xmax =  50;
  }
  if (Method == 12)
    {
      htitleA = "hData_VtxQualityTight_Hemileadingpt_Focus_2Vtx";
      htitleB = "hData_VtxQualityLoose_Hemileadingpt_Focus_2Vtx";
      HeaderA = "Hemi_{1}";
      HeaderNVtx = "2 Vtx"; 
      SaveFile = "Tight_Hemileadingpt_Focus_2Vtx";
      xtitle = " p_{T} [GeV]";
            nbin = 2; 
      xmin = 30;
      xmax =  50;
    }
  if (Method == 13)
    {
      htitleA = "hData_VtxQualityTight_Hemisubleadingpt_Focus_2Vtx";
      htitleB = "hData_VtxQualityLoose_Hemisubleadingpt_Focus_2Vtx";
      HeaderA = "Hemi_{2}";
      HeaderNVtx = "2 Vtx"; 
      SaveFile = "Tight_Hemisubleadingpt_Focus_2Vtx";
      xtitle = " p_{T} [GeV]";
            nbin = 2; 
      xmin = 30;
      xmax =  50;
    }

TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);
c1->SetTicks(1, 1);

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

pad1->SetTicks(0,2);


TPad* rap1 = new TPad("rap1","This is rap1",0.04,0.02,0.96,0.32,21);
rap1->SetFillColor(0);
rap1->SetBorderMode(0);
rap1->SetFrameFillColor(10);
rap1->Draw();
rap1->SetLogy(0);
   rap1->SetTopMargin(0.1);
   rap1->SetBottomMargin(0.25);
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
// gStyle->SetPadTickX(1); 
// gStyle->SetPadTickY(1);
// gStyle->SetPadTicksy(1)
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
//    TH1F* g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);//ok
//  TH1F* g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);//ok
 TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

  TH1F* g1_ST_B = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleB);//ok
  TH1F* g2_ST_B = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleB);//ok
  // TH1F* g3_ST_B = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleB);//ok
  // TH1F* g4_ST_B = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleB);//ok
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
//  std::cout<< "g1_DY = " << MCFILE[0]+htitleA << std::endl;
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 g1_DY_B = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);
  // std::cout<< "g1_DY_B = " << MCFILE[0]+htitleB << std::endl;
  h_DY_B = new TH1F("h_DY_B","",nbin,xmin,xmax);
 h_DY_B->Add(g1_DY_B, h_DY_B, 1,0);


 f2_DY->cd();
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
//  std::cout<< "g2_DY = " << MCFILE[1]+htitleA << std::endl;
  h_DY->Add(g2_DY, h_DY, 1, 1);

 g2_DY_B = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleB);
  // std::cout<< "g2_DY_B = " << MCFILE[1]+htitleB << std::endl;
  h_DY_B->Add(g2_DY_B, h_DY_B, 1, 1);

 f1_VV->cd();
  g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);
  // std::cout<< "g1_VV = " << MCFILE[9]+htitleA << std::endl;
  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

  g1_VV_B = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleB);
  // std::cout<< "g1_VV_B = " << MCFILE[9]+htitleB << std::endl;
  h_VV_B = new TH1F("h_VV_B","",nbin,xmin,xmax);
 h_VV_B->Add(g1_VV_B, h_VV_B, 1,0);


 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);
//  std::cout<< "g2_VV = " << MCFILE[10]+htitleA << std::endl;
  h_VV->Add(g2_VV, h_VV,1, 1);

 g2_VV_B = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleB);
//  std::cout<< "g2_VV_B = " << MCFILE[10]+htitleB << std::endl;
  h_VV_B->Add(g2_VV_B, h_VV_B,1, 1);


 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);
//  std::cout<< "g3_VV = " << MCFILE[11]+htitleA << std::endl; 
  h_VV->Add(g3_VV, h_VV, 1, 1);

 g3_VV_B = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleB);
//  std::cout<< "g3_VV_B = " << MCFILE[11]+htitleB << std::endl;
  h_VV_B->Add(g3_VV_B, h_VV_B, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);
//  std::cout<< "g1_TTV = " << MCFILE[6]+htitleA << std::endl;
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 g1_TTV_B = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleB);
//  std::cout<< "g1_TTV_B = " << MCFILE[6]+htitleB << std::endl;
  h_TTV_B = new TH1F("h_TTV_B","",nbin,xmin,xmax);
 h_TTV_B->Add(g1_TTV_B, h_TTV_B, 1,0);


 f2_TTV->cd();
  g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);
  // // std::cout<< "g2_TTV = " << MCFILE[7]+htitleA << std::endl;
  h_TTV->Add(g2_TTV, h_TTV, 1, 1);

  g2_TTV_B = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleB);
  // std::cout<< "g2_TTV_B = " << MCFILE[7]+htitleB << std::endl;
  h_TTV_B->Add(g2_TTV_B, h_TTV_B, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);
//  std::cout<< "g3_TTV = " << MCFILE[8]+htitleA << std::endl;
  h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 g3_TTV_B = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleB);
//  std::cout<< "g3_TTV_B = " << MCFILE[8]+htitleB << std::endl;
  h_TTV_B->Add(g3_TTV_B, h_TTV_B, 1, 1);

 f1_ST->cd();
  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);
  //  std::cout<< "g1_ST = " << MCFILE[4]+htitleA << std::endl;
 h_ST->Add(g1_ST, h_ST, 1,0);

  h_ST_B = new TH1F("h_ST_B","",nbin,xmin,xmax);
 g1_ST_B = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleB);
//  std::cout<< "g1_ST_B = " << MCFILE[4]+htitleB << std::endl;
 h_ST_B->Add(g1_ST_B, h_ST_B, 1,0);

 f2_ST->cd();
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);
//  std::cout<< "g2_ST = " << MCFILE[5]+htitleA << std::endl;
  h_ST->Add(g2_ST, h_ST, 1, 1);

 g2_ST_B = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleB);
//  std::cout<< "g2_ST_B = " << MCFILE[5]+htitleB << std::endl;
  h_ST_B->Add(g2_ST_B, h_ST_B, 1, 1);

  // if (Channel=="EMU")
  //   {
  //     f3_ST->cd();
  //     g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);
  //     h_ST->Add(g3_ST, h_ST, 1, 1);

  //     g3_ST_B = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleB);
  //     h_ST_B->Add(g3_ST_B, h_ST_B, 1, 1);

  //     f4_ST->cd();
  //     g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);
  //     h_ST->Add(g4_ST, h_ST, 1, 1);
  //     g4_ST_B = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleB);
  //     h_ST_B->Add(g4_ST_B, h_ST_B, 1, 1);
  //   }

 f1_TT->cd();
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
//  std::cout<< "g1_TT = " << MCFILE[2]+htitleA << std::endl;
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, rwTT*1,0);

  g1_TT_B = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);
  // std::cout<< "g1_TT_B = " << MCFILE[2]+htitleB << std::endl;
  h_TT_B = new TH1F("h_TT_B","",nbin,xmin,xmax);
  h_TT_B->Add(g1_TT_B, h_TT_B, rwTT*1,0);

    f2_TT->cd();
    g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);
    // std::cout<< "g2_TT = " << MCFILE[3]+htitleA << std::endl;
    h_TT->Add(g2_TT, h_TT, 1,1);

    g2_TT_B = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleB);
    // std::cout<< "g2_TT_B = " << MCFILE[3]+htitleB << std::endl;
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
 std::cout<<" here 1"<<std::endl;
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

// !! ----------------------------------------------------!!//
// Fit par une droite pour avoir la pente => corrélation

// TF1 *Droite = new TF1("Droite", "[0]*x+[1]",100 , htotMC->GetXaxis()->GetXmax());//htotMC->GetXaxis()->GetXmin()
// Droite->SetParameters(0., 0.04);
// htotMC->Fit(Droite, "R");
// Droite->SetLineColor(kBlack);
// htotMC->Draw("same");

// double a = Droite->GetParameter(0);
// double b = Droite->GetParameter(1);
// double aerr = Droite->GetParError(0);
// double berr = Droite->GetParError(1);
// TF1 *DroiteUp = new TF1("DroiteUp", "[0]*x+[1]",100 , htotMC->GetXaxis()->GetXmax());//htotMC->GetXaxis()->GetXmin()
// DroiteUp->SetParameter(0, a-aerr); // has the higher [1]
// DroiteUp->SetParameter(1, b+berr); 
// DroiteUp->SetLineColor(kBlue);
// DroiteUp->Draw("same");

// TF1 *DroiteDown = new TF1("DroiteDown", "[0]*x+[1]",100 , htotMC->GetXaxis()->GetXmax());//htotMC->GetXaxis()->GetXmin()
// DroiteDown->SetParameter(0, a+aerr); //has the lower [1]
// DroiteDown->SetParameter(1, b-berr); //
// DroiteDown->SetLineColor(kRed);
// DroiteDown->Draw("same");

TF1 *CONST = new TF1("CONST", "[0]",100 , htotMC->GetXaxis()->GetXmax());//htotMC->GetXaxis()->GetXmin()
CONST->SetParameter(0,0.04);
htotMC->Fit(CONST, "R");
// CONST->SetLineColor(kGreen);
// CONST->SetLineWidth(1);
// CONST->Draw("same Y+");

double c = CONST->GetParameter(0);
double cerr = CONST->GetParError(0);

TF1 *CONSTUP = new TF1("CONSTUP", "[0]",100 , htotMC->GetXaxis()->GetXmax());//htotMC->GetXaxis()->GetXmin()
CONSTUP->SetParameter(0,c+cerr);
// htotMC->Fit(CONSTUP, "R");
CONSTUP->SetLineColor(kBlue+3);
CONSTUP->SetLineWidth(1);
CONSTUP->Draw("same");

TF1 *CONSTDOWN = new TF1("CONSTDOWN", "[0]",100 , htotMC->GetXaxis()->GetXmax());//htotMC->GetXaxis()->GetXmin()
CONSTDOWN->SetParameter(0,c-cerr);
// htotMC->Fit(CONSTUP, "R");
CONSTDOWN->SetLineColor(kBlue-9);
CONSTDOWN->SetLineWidth(1);
CONSTDOWN->Draw("same");

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

//  leg = new TLegend(0.7,0.70,0.8,0.75);
//   leg->SetBorderSize(0);
//   leg->SetFillColor(kWhite);
//   leg->SetTextFont(42);
//   leg->SetTextSize(0.035);
//   leg->SetMargin(0.2);
//   leg->SetHeader(HeaderNVtx);
//   leg->Draw();

  //  leg = new TLegend(0.7,0.65,0.8,0.7);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.035);
  // leg->SetMargin(0.2);
  // leg->SetHeader(HeaderA);
  // leg->Draw();

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
// gStyle->SetPadGridX(true);
    rap1->SetGrid();
  TH1F* h_ratio = (TH1F*)htotMC->Clone("h_ratio");
  TH1F* h_ratioUp = (TH1F*)htotMC->Clone("h_ratioUp");
  TH1F* h_ratioDown = (TH1F*)htotMC->Clone("h_ratioDown");
  h_ratio->Reset();
  h_ratioUp->Reset();
  h_ratioDown->Reset();
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

    h_ratioDown->SetFillStyle(1001);
  h_ratioDown->SetFillColorAlpha(kBlue-9, 1);
  h_ratioDown->SetLineColor(kBlue-9);
  h_ratioDown->SetLineWidth(1);
  h_ratioDown->Draw("PE1 same");

  leg = new TLegend(0.7,0.65,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.055);
  leg->SetMargin(0.2);

  leg->AddEntry(h_ratio,"Fit","L");
  leg->AddEntry(h_ratioUp,"FitUp","L");
  leg->AddEntry(h_ratioDown,"FitDown","L");
  leg->Draw();

  if (Channel == "MUMU")
    {
      SaveFile = SaveFile + "_MUMU_MC";
    }
  else
    {
      SaveFile = SaveFile + "_EMU_MC";
    }
  c1->SaveAs(SaveFile+"_"+Year+".pdf");
  rap1->SaveAs(SaveFile+"_"+Year+".root");
  delete c1;
}