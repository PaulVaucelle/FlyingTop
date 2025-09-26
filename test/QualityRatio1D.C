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
TString suffixDATA = "_BDT100_Corr";//
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


  if (Channel == "MUMU") 
  {
    suffixDATA = "_Corr";
        if (Year == "2016PRE") 
      {
         suffixDATA = "";
        Yearcor = "2016preVFP";
        Prod = "DATA_MUMU_2016PRE_30_03_2025";
        SampleDATA = "DoubleMuon_Run2016-HIPM_UL2016_MiniAODv2";
        SampleDATAExtra = SampleDATA+suffixDATA;
      }
      
    if (Year == "2016POST") 
      {
         suffixDATA = "";
        Yearcor = "2016";
        Prod = "DATA_MUMU_2016POST_30_03_2025";
        SampleDATA = "DoubleMuon_Run2016-UL2016_MiniAODv2";
        SampleDATAExtra = SampleDATA+suffixDATA;

      }
    if (Year == "2017") 
      {
         suffixDATA = "";
        Yearcor = "2017";
        Prod = "DATA_MUMU_2017_30_03_2025";
        SampleDATA = "DoubleMuon_Run2017-UL2017_MiniAODv2";
        SampleDATAExtra = SampleDATA +suffixDATA;

      }
    if (Year == "2018") 
      {
        Prod = "DATA_MUMU_2018_03_02_2025";
        SampleDATA = "DoubleMuon_UL2018_MiniAODv2_GT36-v1";
        SampleDATAExtra = SampleDATA+suffixDATA;
       
      }
  }




TFile* f1_Data = new TFile("../../"+Prod+"/histofile_HT100_"+Dmode+"_OS_2p4_"+SampleDATAExtra+".root");

std::cout<<"Opening file: ../../"+Prod+"/histofile_HT100_"+Dmode+"_OS_2p4_"+SampleDATAExtra+".root"<<std::endl;
TString FILE[1] = { SampleDATA+"_"};

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
      xtitle = " Hemi p_{T} [GeV]"; 

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
      xtitle = " Hemi Leading p_{T} [GeV]";
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
      xtitle = " Hemi SubLeading p_{T} [GeV]";
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
gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);
gROOT->SetBatch(kTRUE);

  f1_Data->cd();
 TH1F* g1_Data = (TH1F*)gROOT->FindObject(FILE[0]+htitleA);//ok
 TH1F* g2_Data = (TH1F*)gROOT->FindObject(FILE[0]+htitleB);//ok

//  g1_Data->Sumw2();

  c1->cd();
pad1->cd();
  g1_Data->Divide(g1_Data,g2_Data,1,1);

  g1_Data->SetFillStyle(1001);
//  g1_Data->SetFillColorAlpha(kGreen+1, 1);
 g1_Data->SetLineColor(kBlack);
 g1_Data->Draw("PE1");
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
 g1_Data->GetXaxis()->SetRangeUser(0,260);
 g1_Data->GetYaxis()->SetTitle(ytitle);
 g1_Data->GetYaxis()->SetTitleColor(1);
 g1_Data->SetNdivisions(509,"XYZ");
 g1_Data->SetMinimum(hmin); 
 g1_Data->SetMaximum(g1_Data->GetMaximum()*2.5); 
 g1_Data->SetTitle(""); 

 
 g1_Data->SetMinimum(0); 
//  g1_Data->SetMaximum(hmax);

//----------------------------------------------------//
// Fit par une droite pour avoir la pente => corrélation

// // TF1 *Droite = new TF1("Droite", "[0]*x+[1]",0, g1_Data->GetXaxis()->GetXmax());
// // Droite->SetParameters(0., 0.04);
// // g1_Data->Fit(Droite, "R");
// // Droite->SetLineColor(kBlack);
// // g1_Data->Draw("same");

// // double a = Droite->GetParameter(0);
// // double b = Droite->GetParameter(1);
// // double aerr = Droite->GetParError(0);
// // double berr = Droite->GetParError(1);
// // TF1 *DroiteUp = new TF1("DroiteUp", "[0]*x+[1]",0 , g1_Data->GetXaxis()->GetXmax());//g1_Data->GetXaxis()->GetXmin()
// // DroiteUp->SetParameter(0, a-aerr); // has the higher [1]
// // DroiteUp->SetParameter(1, b+berr); 
// // DroiteUp->SetLineColor(kBlue);
// // DroiteUp->Draw("same");

// // TF1 *DroiteDown = new TF1("DroiteDown", "[0]*x+[1]",0 , g1_Data->GetXaxis()->GetXmax());//g1_Data->GetXaxis()->GetXmin()
// // DroiteDown->SetParameter(0, a+aerr); //has the lower [1]
// // DroiteDown->SetParameter(1, b-berr); //
// // DroiteDown->SetLineColor(kRed);
// // DroiteDown->Draw("same");

TF1 *CONST = new TF1("CONST", "[0]",100 , g1_Data->GetXaxis()->GetXmax());//g1_Data->GetXaxis()->GetXmin()
CONST->SetParameter(0,0.04);
// g1_Data->Fit(CONST, "R");
// CONST->SetLineColor(kGray);
// CONST->SetLineWidth(1);
// CONST->Draw("same");

double c = CONST->GetParameter(0);
double cerr = CONST->GetParError(0);

TF1 *CONSTUP = new TF1("CONSTUP", "[0]",100 , g1_Data->GetXaxis()->GetXmax());//g1_Data->GetXaxis()->GetXmin()
CONSTUP->SetParameter(0,c+cerr);
// g1_Data->Fit(CONSTUP, "R");
CONSTUP->SetLineColor(kBlue+3);
CONSTUP->SetLineWidth(1);
// CONSTUP->Draw("same");

TF1 *CONSTDOWN = new TF1("CONSTDOWN", "[0]",100 , g1_Data->GetXaxis()->GetXmax());//g1_Data->GetXaxis()->GetXmin()
CONSTDOWN->SetParameter(0,c-cerr);
// g1_Data->Fit(CONSTUP, "R");
CONSTDOWN->SetLineColor(kBlue-9);
CONSTDOWN->SetLineWidth(1);
// CONSTDOWN->Draw("same");


// !! --------------------------- !!//

TString extraLeg = " #mu#mu";
if (Channel == "EMU") extraLeg = " e#mu";
  TLegend* leg = new TLegend(0.7,0.65,0.8,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(g1_Data,"data "+extraLeg,"PE1");
  // // leg->AddEntry(Droite,"Linear Fit","L");
  // // leg->AddEntry(DroiteUp,"FitUp","L");
  // // leg->AddEntry(DroiteDown,"FitDown","L");
  // leg->AddEntry(CONST,"Const Fit","L");
  // leg->AddEntry(CONSTUP,"Up Fit","L");
  // leg->AddEntry(CONSTDOWN,"Down Fit","L");
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

  // if (Method >= 8)
  //   {
  //     leg = new TLegend(0.5,0.55,0.85,0.65);
  //     leg->SetBorderSize(0);
  //     leg->SetFillColor(kWhite);
  //     leg->SetTextFont(42);
  //     leg->SetTextSize(0.035);
  //     leg->SetMargin(0.2);
  //     // double a = Droite->GetParameter(0);
  //     // double b = Droite->GetParameter(1);
  //     // double aerr = Droite->GetParError(0);
  //     // double berr = Droite->GetParError(1);
  //     // leg->SetHeader("Slope = "+TString::Format("%.2e #pm %.2e",a,aerr));
  //     double c = CONST->GetParameter(0);
  //     double cerr = CONST->GetParError(0);
  //     leg->AddEntry(CONST,TString::Format("Const = %.2e #pm %.2e",c,cerr),"L");
  //     leg->Draw();
  //   }

  PlotCMSv2(pad1,Year,true);

//     rap1->cd();
//   TH1F* h_ratio = (TH1F*)g1_Data->Clone("h_ratio");
//   TH1F* h_ratioUp = (TH1F*)g1_Data->Clone("h_ratioUp");
//   TH1F* h_ratioDown = (TH1F*)g1_Data->Clone("h_ratioDown");
//   h_ratio->Reset();
//   h_ratioUp->Reset();
//   h_ratioDown->Reset();
//   // h_ratio->SetTitle("Ratio: data / fit extrapolation");

//   int nBins = g1_Data->GetNbinsX();
//   for (int i = 1; i <= nBins; ++i) {
//     double x = g1_Data->GetBinCenter(i);
//     if (x >= 30 && x <= 100) {

      
//       // double y_data = g1_Data->GetBinContent(i);
//       // double y_fit = Droite->Eval(x);
//       // double y_fit_up = DroiteUp->Eval(x);
//       // double y_fit_down = DroiteDown->Eval(x);
//       // if (y_fit != 0) {
//       //   h_ratio->SetBinContent(i, y_data / y_fit);
//       //   h_ratio->SetBinError(i, g1_Data->GetBinError(i) / y_fit); // propagation simple
//       // }
//       // if (y_fit_up != 0) {
//       //   h_ratioUp->SetBinContent(i, y_data / y_fit_up);
//       //   h_ratioUp->SetBinError(i, g1_Data->GetBinError(i) / y_fit_up); // propagation simple
//       // }
//       // if (y_fit_down != 0) {
//       //   h_ratioDown->SetBinContent(i, y_data / y_fit_down);
//       //   h_ratioDown->SetBinError(i, g1_Data->GetBinError(i) / y_fit_down); // propagation simple
//       // }


//             double y_data = g1_Data->GetBinContent(i);
//       double y_fit = CONST->Eval(x);
//       double y_fit_up = CONSTUP->Eval(x);
//       double y_fit_down = CONSTDOWN->Eval(x);
//       if (y_fit != 0) {
//         h_ratio->SetBinContent(i, y_data / y_fit);
//         h_ratio->SetBinError(i, g1_Data->GetBinError(i) / y_fit); // propagation simple
//       }
//       if (y_fit_up != 0) {
//         h_ratioUp->SetBinContent(i, y_data / y_fit_up);
//         h_ratioUp->SetBinError(i, g1_Data->GetBinError(i) / y_fit_up); // propagation simple
//       }
//       if (y_fit_down != 0) {
//         h_ratioDown->SetBinContent(i, y_data / y_fit_down);
//         h_ratioDown->SetBinError(i, g1_Data->GetBinError(i) / y_fit_down); // propagation simple
//       }
//     }
//    else if (x >= 100)
//     {
//       double y_data = g1_Data->GetBinContent(i);
//       double y_fit = CONST->Eval(x);
//       double y_fit_up = CONSTUP->Eval(x);
//       double y_fit_down = CONSTDOWN->Eval(x);
//       if (y_fit != 0) {
//         h_ratio->SetBinContent(i, y_data / y_fit);
//         h_ratio->SetBinError(i, g1_Data->GetBinError(i) / y_fit); // propagation simple
//       }
//       if (y_fit_up != 0) {
//         h_ratioUp->SetBinContent(i, y_data / y_fit_up);
//         h_ratioUp->SetBinError(i, g1_Data->GetBinError(i) / y_fit_up); // propagation simple
//       }
//       if (y_fit_down != 0) {
//         h_ratioDown->SetBinContent(i, y_data / y_fit_down);
//         h_ratioDown->SetBinError(i, g1_Data->GetBinError(i) / y_fit_down); // propagation simple
//       }
//     }
//   }
//   h_ratio->SetFillStyle(1001);
// //  h_ratio->SetFillColorAlpha(kGreen+1, 1);
//  h_ratio->SetLineColor(kRed);
//  h_ratio->Draw("PE1");
//  h_ratio->SetMarkerStyle(20);
//  h_ratio->SetMarkerSize(1.5);
//  h_ratio->SetMarkerColor(kRed);
//  h_ratio->SetLineColor(kRed);
//  h_ratio->SetLineWidth(1);
//  h_ratio->SetTickLength(0.03, "YZ");
//  h_ratio->SetTickLength(0.03,"X");
//  h_ratio->SetLabelOffset(0.01,"X");
//  h_ratio->SetLabelOffset(0.007,"Y");
//  h_ratio->SetLabelSize(0.06, "XYZ");
//  h_ratio->SetLabelFont(42, "XYZ"); 
//  h_ratio->SetTitleSize(0.085, "XYZ"); 
//  h_ratio->SetTitleFont(42, "XYZ");
//  h_ratio->SetTitleOffset(1.2,"X"); 
//  h_ratio->SetTitleOffset(0.8,"Y");
//  h_ratio->GetXaxis()->SetTitle(xtitle);
//  h_ratio->GetXaxis()->SetTitleColor(1);
//  h_ratio->GetXaxis()->SetRangeUser(0,260);
// //  h_ratio->GetYaxis()->SetTitle(ytitle);
//  h_ratio->GetYaxis()->SetTitleColor(1);
//  h_ratio->SetNdivisions(509,"XYZ");
//  h_ratio->SetMinimum(0.); 
//  h_ratio->SetMaximum(2.); 
//  h_ratio->SetTitle(""); 

//   h_ratio->GetYaxis()->SetTitle("DATA / Fit");

//   h_ratioUp->SetFillStyle(1001);
//   h_ratioUp->SetFillColorAlpha(kBlue+3, 1);
//   h_ratioUp->SetLineColor(kBlue+3);
//   h_ratioUp->SetLineWidth(1);
//   h_ratioUp->Draw("PE1 same");

//     h_ratioDown->SetFillStyle(1001);
//   h_ratioDown->SetFillColorAlpha(kBlue-9, 1);
//   h_ratioDown->SetLineColor(kBlue-9);
//   h_ratioDown->SetLineWidth(1);
//   h_ratioDown->Draw("PE1 same");

//   leg = new TLegend(0.7,0.65,0.8,0.85);
//   leg->SetBorderSize(0);
//   leg->SetFillColor(kWhite);
//   leg->SetTextFont(42);
//   leg->SetTextSize(0.055);
//   leg->SetMargin(0.2);

//   leg->AddEntry(h_ratio,"Fit","L");
//   leg->AddEntry(h_ratioUp,"FitUp","L");
//   leg->AddEntry(h_ratioDown,"FitDown","L");
//   leg->Draw();

  if (Channel == "MUMU")
    {
      SaveFile = SaveFile + "_MUMU_DATA";
    }
  else
    {
      SaveFile = SaveFile + "_EMU_DATA";
    }
  c1->SaveAs(SaveFile+"_"+Year+".pdf");
  rap1->SaveAs(SaveFile+"_"+Year+".root");
  delete c1;
}