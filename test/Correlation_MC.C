#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PlotCMS.h"

void plot() {
    int stati=0;
  // Ouvre les fichiers ROOT
  TFile* fileEMU = TFile::Open("VtxQuality_Hemipt_2Vtx_EMU_MC.root");
  TFile* fileMUMU = TFile::Open("VtxQuality_Hemipt_2Vtx_MUMU_MC.root");
    bool logy = 0;
  if (!fileEMU || fileEMU->IsZombie()) {
    std::cerr << "Erreur: Impossible d'ouvrir le fichier EMU" << std::endl;
    return;
  }

  if (!fileMUMU || fileMUMU->IsZombie()) {
    std::cerr << "Erreur: Impossible d'ouvrir le fichier MUMU" << std::endl;
    return;
  }

  // Récupère les canvas c1 dans chaque fichier
  TCanvas* c1EMU = (TCanvas*)fileEMU->Get("c1");
  TCanvas* c1MUMU = (TCanvas*)fileMUMU->Get("c1");

  if (!c1EMU || !c1MUMU) {
    std::cerr << "Erreur: Canvas 'c1' introuvable dans un des fichiers" << std::endl;
    return;
  }

  // Récupère les histos "htotMC" à partir des canvas
  TH1F* h_EMU = (TH1F*)c1EMU->FindObject("htotMC");
  TH1F* h_MUMU = (TH1F*)c1MUMU->FindObject("htotMC");

  if (!h_EMU || !h_MUMU) {
    std::cerr << "Erreur: Histogramme 'htotMC' introuvable dans un des canvas" << std::endl;
    return;
  }

  // Crée un nouveau canvas pour superposer les deux histos
  TCanvas *c_comp = new TCanvas("c_comp", "plots",0,0,1300,1200);
    c_comp->SetFillColor(10);
    c_comp->SetFillStyle(4000);
    c_comp->SetBorderSize(2);
    c_comp->cd();
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

    TPad* pad2 = new TPad("pad2","This is pad2",0.04,0.05,0.96,0.3,21);
    pad2->SetFillColor(0);
    pad2->SetBorderMode(0);
    pad2->SetFrameFillColor(10);
    pad2->Draw();
    pad2->SetLogy(logy);
    pad2->SetTopMargin(0.1);
    pad2->SetBottomMargin(0.15);
    pad2->SetRightMargin(0.05);
    pad2->SetLeftMargin(0.15);
    pad2->SetGrid();
    pad1->cd();

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


  h_EMU->SetMarkerStyle(4);
  h_MUMU->SetMarkerStyle(20);
h_EMU->SetMarkerSize(1.5);
h_MUMU->SetMarkerSize(1.5);
h_EMU->SetMarkerColor(kBlack);
h_MUMU->SetMarkerColor(kBlack);

  h_EMU->Draw("LP");
  h_MUMU->Draw("LP SAME");

  // Légende pour faire joli
  TLegend* leg = new TLegend(0.7, 0.75, 0.95, 0.90);
  leg->AddEntry(h_EMU, "e#mu MC", "lp");
  leg->AddEntry(h_MUMU, "#mu#mu MC", "lp");
  leg->Draw();
  PlotCMSv2(pad1,"2018",false);
  pad2->cd();
    TH1F* h_ratio = (TH1F*)h_MUMU->Clone("h_ratio");
    h_ratio->Divide(h_EMU);
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetMarkerColor(kBlack);
    h_ratio->SetLineWidth(2);
    h_ratio->SetMarkerStyle(4);
    h_ratio->SetMinimum(0.5);
    h_ratio->SetMaximum(1.5);
    h_ratio->GetXaxis()->SetTitle("x-axis title");
    h_ratio->GetYaxis()->SetTitle("#mu#mu / e#mu");
    h_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_ratio->GetYaxis()->SetTitleSize(0.1);
    h_ratio->GetXaxis()->SetTitleSize(0.1);
    h_ratio->GetXaxis()->SetLabelSize(0.1);
    h_ratio->GetYaxis()->SetLabelSize(0.1);
    h_ratio->GetXaxis()->SetLabelFont(42);
    h_ratio->GetYaxis()->SetLabelFont(42);
    h_ratio->GetXaxis()->SetTitleFont(42);
    h_ratio->GetYaxis()->SetTitleFont(42);
    // h_ratio->GetXaxis()->SetTickLength(0.03, "YZ");
    // h_ratio->GetYaxis()->SetTickLength(0.03,"X");
    h_ratio->GetYaxis()->SetNdivisions(509,"XYZ");
    h_ratio->Draw("LP");
  c_comp->Update();
  c_comp->SaveAs("Correlation_MC_Ratio.pdf");
  c_comp->SaveAs("Correlation_MC_Ratio.root");
}
