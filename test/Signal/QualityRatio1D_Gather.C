#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PlotCMS.h"

void plot() {
    int stati=0;
    bool fit= 1;
    bool logy=0;

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

    // Noms des fichiers et des histogrammes
    const char* filenames[5] = {
        "VtxQuality_Hemipt_2Vtx_EMU_DATA_2024.root",
        "VtxQuality_Hemipt_2Vtx_EMU_DATA_2023A.root",
        "VtxQuality_Hemipt_2Vtx_EMU_DATA_2023B.root",
        "VtxQuality_Hemipt_2Vtx_EMU_DATA_2022A.root",
        "VtxQuality_Hemipt_2Vtx_EMU_DATA_2022B.root"
    };

    const char* histonames[5] = {
        "MuonEG_Run2024_hData_VtxQualityTight_Hemipt_2Vtx",
        "MuonEG_Run2023C-22Sep2023_hData_VtxQualityTight_Hemipt_2Vtx",
        "MuonEG_Run2023D-22Sep2023_hData_VtxQualityTight_Hemipt_2Vtx",
        "MuonEG_Run2022-CDE-22Sep2023_hData_VtxQualityTight_Hemipt_2Vtx",
        "MuonEG_Run2022-FG-22Sep2023_hData_VtxQualityTight_Hemipt_2Vtx"
    };

    TH1F* hSum = nullptr;

    for (int i = 0; i < 5; ++i) {
        TFile* file = TFile::Open(filenames[i]);
        if (!file || file->IsZombie()) {
            std::cerr << "Erreur à l'ouverture du fichier : " << filenames[i] << std::endl;
            continue;
        }

        TPad* pad1 = (TPad*)file->Get("pad1");
        if (!pad1) {
            std::cerr << "TPad 'pad1' introuvable dans " << filenames[i] << std::endl;
            continue;
        }

        pad1->cd(); // S'assurer qu'on est dans le bon contexte
        TH1F* h = (TH1F*)pad1->FindObject(histonames[i]);
        if (!h) {
            std::cerr << "TH1F '" << histonames[i] << "' introuvable dans " << filenames[i] << std::endl;
            continue;
        }

        if (!hSum) {
            hSum = (TH1F*)h->Clone("hSum");
            hSum->SetDirectory(0); // Déconnecter du fichier
        } else {
            hSum->Add(h);
        }

        file->Close();
    }

    if (!hSum) {
        std::cerr << "Aucun histogramme chargé !" << std::endl;
        return;
    }

    // Création du canvas avec 2 pads
    TCanvas* c = new TCanvas("c", "Somme des histos VtxQuality", 0,0,1300,1200);
    c->SetFillColor(10);
    c->SetFillStyle(4000);
    c->SetBorderSize(2);

    TPad* padTop = new TPad("padTop", "padTop", 0.0, 0.3, 1.0, 1.0);
    padTop->SetFillColor(0);
    padTop->SetBorderMode(0);
    padTop->SetFrameFillColor(10);
    padTop->Draw();
    padTop->SetLogy(0);
    padTop->SetTopMargin(0.1);
    padTop->SetBottomMargin(0.15);
    padTop->SetRightMargin(0.05);
    padTop->SetLeftMargin(0.15);

    TPad* padBot = new TPad("padBot", "padBot", 0.0, 0.0, 1.0, 0.3);
    padBot->SetFillColor(0);
    padBot->SetBorderMode(0);
    padBot->SetFrameFillColor(10);
    padBot->Draw();
    padBot->SetLogy(logy);
    padBot->SetTopMargin(0.1);
    padBot->SetBottomMargin(0.15);
    padBot->SetRightMargin(0.05);
    padBot->SetLeftMargin(0.15);

    padTop->Draw();
    padBot->Draw();

    padTop->cd();
    hSum->Draw("PE1"); // Avec erreurs
    hSum->SetMinimum(0.055); 
    hSum->SetMaximum(hSum->GetMaximum()*2.5); 
    TF1 *CONST = new TF1("CONST", "[0]",0 , hSum->GetXaxis()->GetXmax());//hSum->GetXaxis()->GetXmin()
    CONST->SetParameter(0,0.04);
    hSum->Fit(CONST, "R");
    CONST->SetLineColor(kGreen);
    CONST->SetLineWidth(1);
    CONST->Draw("same");
    double par = CONST->GetParameter(0);
    double parerr = CONST->GetParError(0);



    TF1 *CONSTUP = new TF1("CONSTUP", "[0]",0 , hSum->GetXaxis()->GetXmax());//hSum->GetXaxis()->GetXmin()
    CONSTUP->SetParameter(0,par+parerr);
    // hSum->Fit(CONSTUP, "R");
    CONSTUP->SetLineColor(kGreen+5);
    CONSTUP->SetLineWidth(1);
    CONSTUP->Draw("same");

    TF1 *CONSTDOWN = new TF1("CONSTDOWN", "[0]",0 , hSum->GetXaxis()->GetXmax());//hSum->GetXaxis()->GetXmin()
    CONSTDOWN->SetParameter(0,par-parerr);
    // hSum->Fit(CONSTUP, "R");
    CONSTDOWN->SetLineColor(kGreen-5);
    CONSTDOWN->SetLineWidth(1);
    CONSTDOWN->Draw("same");


TLegend* leg = new TLegend(0.7,0.65,0.8,0.89);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  leg->SetMargin(0.2);
  leg->AddEntry(hSum,"data e#mu","PE1");
  leg->AddEntry(CONST,"Const Fit","L");
  leg->AddEntry(CONSTUP,"Const FitUp","L");
  leg->AddEntry(CONSTDOWN,"Const FitDown","L");
  leg->Draw();

leg = new TLegend(0.5,0.55,0.85,0.65);
leg->SetBorderSize(0);
leg->SetFillColor(kWhite);
leg->SetTextFont(42);
leg->SetTextSize(0.035);
leg->SetMargin(0.2);

leg->AddEntry(CONST,TString::Format("Const = %.2e #pm %.2e",par,parerr),"L");
leg->Draw();
    

    PlotCMSv2(padTop,"2018",false);

    // À toi de jouer dans padBot !
    padBot->cd();
    // Tu peux dessiner ici ton ratio, résidu ou tout ce que tu veux

  TH1F* h_ratio = (TH1F*)hSum->Clone("h_ratio");
  TH1F* h_ratioUp = (TH1F*)hSum->Clone("h_ratioUp");
  TH1F* h_ratioDown = (TH1F*)hSum->Clone("h_ratioDown");
  h_ratio->Reset();
  h_ratioUp->Reset();
  h_ratioDown->Reset();
  // h_ratio->SetTitle("Ratio: data / fit extrapolation");

  int nBins = hSum->GetNbinsX();
  for (int i = 1; i <= nBins; ++i) {
    double x = hSum->GetBinCenter(i);
    if (x >= 30 && x <= 100) {

      
      // double y_data = hSum->GetBinContent(i);
      // double y_fit = Droite->Eval(x);
      // double y_fit_up = DroiteUp->Eval(x);
      // double y_fit_down = DroiteDown->Eval(x);
      // if (y_fit != 0) {
      //   h_ratio->SetBinContent(i, y_data / y_fit);
      //   h_ratio->SetBinError(i, hSum->GetBinError(i) / y_fit); // propagation simple
      // }
      // if (y_fit_up != 0) {
      //   h_ratioUp->SetBinContent(i, y_data / y_fit_up);
      //   h_ratioUp->SetBinError(i, hSum->GetBinError(i) / y_fit_up); // propagation simple
      // }
      // if (y_fit_down != 0) {
      //   h_ratioDown->SetBinContent(i, y_data / y_fit_down);
      //   h_ratioDown->SetBinError(i, hSum->GetBinError(i) / y_fit_down); // propagation simple
      // }


            double y_data = hSum->GetBinContent(i);
      double y_fit = CONST->Eval(x);
      double y_fit_up = CONSTUP->Eval(x);
      double y_fit_down = CONSTDOWN->Eval(x);
      if (y_fit != 0) {
        h_ratio->SetBinContent(i, y_data / y_fit);
        h_ratio->SetBinError(i, hSum->GetBinError(i) / y_fit); // propagation simple
      }
      if (y_fit_up != 0) {
        h_ratioUp->SetBinContent(i, y_data / y_fit_up);
        h_ratioUp->SetBinError(i, hSum->GetBinError(i) / y_fit_up); // propagation simple
      }
      if (y_fit_down != 0) {
        h_ratioDown->SetBinContent(i, y_data / y_fit_down);
        h_ratioDown->SetBinError(i, hSum->GetBinError(i) / y_fit_down); // propagation simple
      }
    }
   else if (x >= 100)
    {
      double y_data = hSum->GetBinContent(i);
      double y_fit = CONST->Eval(x);
      double y_fit_up = CONSTUP->Eval(x);
      double y_fit_down = CONSTDOWN->Eval(x);
      if (y_fit != 0) {
        h_ratio->SetBinContent(i, y_data / y_fit);
        h_ratio->SetBinError(i, hSum->GetBinError(i) / y_fit); // propagation simple
      }
      if (y_fit_up != 0) {
        h_ratioUp->SetBinContent(i, y_data / y_fit_up);
        h_ratioUp->SetBinError(i, hSum->GetBinError(i) / y_fit_up); // propagation simple
      }
      if (y_fit_down != 0) {
        h_ratioDown->SetBinContent(i, y_data / y_fit_down);
        h_ratioDown->SetBinError(i, hSum->GetBinError(i) / y_fit_down); // propagation simple
      }
    }
  }
  h_ratio->SetFillStyle(1001);
//  h_ratio->SetFillColorAlpha(kGreen+1, 1);
 h_ratio->SetLineColor(kBlack);
 h_ratio->Draw("PE1");
 h_ratio->SetMarkerStyle(20);
 h_ratio->SetMarkerSize(1.5);
 h_ratio->SetMarkerColor(kBlack);
 h_ratio->SetLineColor(kBlack);
 h_ratio->SetLineWidth(1);
 h_ratio->SetTickLength(0.03, "YZ");
 h_ratio->SetTickLength(0.03,"X");
 h_ratio->SetLabelOffset(0.01,"X");
 h_ratio->SetLabelOffset(0.007,"Y");
 h_ratio->SetLabelSize(0.055, "XYZ");
 h_ratio->SetLabelFont(42, "XYZ"); 
 h_ratio->SetTitleSize(0.065, "XYZ"); 
 h_ratio->SetTitleFont(42, "XYZ");
 h_ratio->SetTitleOffset(1.2,"X"); 
 h_ratio->SetTitleOffset(1.,"Y");
//  h_ratio->GetXaxis()->SetTitle(xtitle);
 h_ratio->GetXaxis()->SetTitleColor(1);
 h_ratio->GetXaxis()->SetRangeUser(0,260);
//  h_ratio->GetYaxis()->SetTitle(ytitle);
 h_ratio->GetYaxis()->SetTitleColor(1);
 h_ratio->SetNdivisions(509,"XYZ");
 h_ratio->SetMinimum(0.75); 
 h_ratio->SetMaximum(1.8); 
 h_ratio->SetTitle(""); 

  h_ratio->GetYaxis()->SetTitle("MC / Fit");

  h_ratioUp->SetFillStyle(1001);
  h_ratioUp->SetFillColorAlpha(kBlue, 1);
  h_ratioUp->SetLineColor(kBlue);
  h_ratioUp->SetLineWidth(1);
  h_ratioUp->Draw("PE1 same");

    h_ratioDown->SetFillStyle(1001);
  h_ratioDown->SetFillColorAlpha(kRed, 1);
  h_ratioDown->SetLineColor(kRed);
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

  // !! -- !! //

    c->SaveAs("QualityRatio1D_Gather.pdf");
     c->SaveAs("QualityRatio1D_Gather.root");
}
