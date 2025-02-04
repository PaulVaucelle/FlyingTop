#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <cmath>
#include "TPad.h"
// Fonction pour calculer les erreurs asymétriques
void calculate_errors(const std::vector<double> &y, const std::vector<double> &PDF_up_percent, const std::vector<double> &PDF_down_percent,
                      const std::vector<double> &SCALE_up_percent, const std::vector<double> &SCALE_down_percent,
                      std::vector<double> &err_up, std::vector<double> &err_down) {
    for (size_t i = 0; i < y.size(); ++i) {
        double PDF_up = y[i] * PDF_up_percent[i] / 100.0;
        double PDF_down = y[i] * PDF_down_percent[i] / 100.0;
        double SCALE_up = y[i] * SCALE_up_percent[i] / 100.0;
        double SCALE_down = y[i] * SCALE_down_percent[i] / 100.0;

        err_up[i] = std::sqrt(PDF_up * PDF_up + SCALE_up * SCALE_up);
        err_down[i] = std::sqrt(PDF_down * PDF_down + SCALE_down * SCALE_down);
    }
}

// Fonction pour calculer les ratios avec erreurs
void calculate_ratio(const std::vector<double> &num, const std::vector<double> &den,
                     const std::vector<double> &num_err_up, const std::vector<double> &den_err_up,
                     const std::vector<double> &num_err_down, const std::vector<double> &den_err_down,
                     std::vector<double> &ratio, std::vector<double> &ratio_err_up, std::vector<double> &ratio_err_down) {
    for (size_t i = 0; i < num.size(); ++i) {
        ratio[i] = num[i] / den[i];
        ratio_err_up[i] = ratio[i] * std::sqrt(std::pow(num_err_up[i] / num[i], 2) + std::pow(den_err_up[i] / den[i], 2));
        ratio_err_down[i] = ratio[i] * std::sqrt(std::pow(num_err_down[i] / num[i], 2) + std::pow(den_err_down[i] / den[i], 2));
    }
}


void plot_XSRR_KFactor_AN() {

        int stati=0;
    bool fit= 0;
    bool logy=0;

    // Points sur l'axe x (identiques pour tous les graphes)
    std::vector<double> x = {200, 250, 300, 350, 400, 450, 500};

    // 13 TeV Values
    // Left-Handed
    // 6 vecteurs y (résultats fictifs pour illustration)
    std::vector<double> y1 = {12.3,5.4,2.66,1.43,0.82,0.49,0.30}; // LO + 1 Jet this Analysis
    std::vector<double> y2 = {7.66,3.25,1.57,0.83,0.46,0.27,0.17}; // LO this analysis
    std::vector<double> y3 = {8.15,3.47,1.68,0.89,0.50,0.30,0.18}; // NLO+NLL B.Fucks
    std::vector<double> y4 = {7.03,3.03,1.48,0.78,0.44,0.26,0.16}; // LO B.Fucks


    // Erreurs PDF (symétriques) pour chaque vecteur y en %
    std::vector<double> y1_PDF_up_percent = {2.1,2.3,2.6,2.8,3.2,3.3,3.7};// LO + 1 Jet this Analysis
    std::vector<double> y1_PDF_down_percent = {2.1,2.3,2.6,2.8,3.2,3.3,3.7};// LO + 1 Jet this Analysis
    std::vector<double> y2_PDF_up_percent = {2.1,2.3,2.5,2.8,2.9,3.2,3.5};// LO this analysis
    std::vector<double> y2_PDF_down_percent = {2.1,2.3,2.5,2.8,2.9,3.2,3.5};// LO this analysis
    std::vector<double> y3_PDF_up_percent = {3.9,4.2,4.5,4.9,5.2,5.7,6.2};// NLO+NLL B.Fucks
    std::vector<double> y3_PDF_down_percent = {4.9,5.3,5.7,6.1,6.7,7.0,7.4};// NLO+NLL B.Fucks
    std::vector<double> y4_PDF_up_percent = {1.3,2.8,3.9,4.8,5.7,6.4,7.0};// LO B.Fucks : no pdf variation provided
    std::vector<double> y4_PDF_down_percent = {1.6,2.9,3.8,4.5,5.2,5.7,6.2};// LO B.Fucks: no pdf variation provided


    // Erreurs SCALE (asymétriques) pour chaque vecteur y en %
    std::vector<double> y1_SCALE_up_percent = {5.7,5.9,6.0,6.2,6.2,6.3,6.7};// LO + 1 Jet this Analysis
    std::vector<double> y1_SCALE_down_percent = {4.4,4.6,4.6,4.8,5.0,5.5,6.0};// LO + 1 Jet this Analysis
    std::vector<double> y2_SCALE_up_percent = {1.7,2.9,3.9,4.7,5.4,6.0,6.6};// LO this analysis
    std::vector<double> y2_SCALE_down_percent = {1.9,2.9,3.7,4.3,4.9,5.4,5.9};// LO this analysis    std::vector<double> y3_SCALE_up_percent = {};// NLO+NLL B.Fucks
    std::vector<double> y3_SCALE_up_percent = {0.1,0.3,0.1,0.0,0.0,0.0,0.0};// NLO+NLL B.Fucks
    std::vector<double> y3_SCALE_down_percent = {3.6,3.9,4.2,4.6,5.0,5.4,5.9};// NLO+NLL B.Fucks
    std::vector<double> y4_SCALE_up_percent = {1.3,2.8,4.0,4.9,5.7,6.4,7.1};// LO B.Fucks
    std::vector<double> y4_SCALE_down_percent = {1.6,2.9,3.8,4.6,5.2,5.8,6.3};// LO B.Fucks

    // Calcul des erreurs asymétriques pour tous les vecteurs y
    std::vector<double> y_err_up1(x.size()), y_err_down1(x.size());
    std::vector<double> y_err_up2(x.size()), y_err_down2(x.size());
    std::vector<double> y_err_up3(x.size()), y_err_down3(x.size());
    std::vector<double> y_err_up4(x.size()), y_err_down4(x.size());


    calculate_errors(y1, y1_PDF_up_percent,y1_PDF_down_percent, y1_SCALE_up_percent, y1_SCALE_down_percent, y_err_up1, y_err_down1);
    calculate_errors(y2, y2_PDF_up_percent,y2_PDF_down_percent, y2_SCALE_up_percent, y2_SCALE_down_percent, y_err_up2, y_err_down2);
    calculate_errors(y3, y3_PDF_up_percent,y3_PDF_down_percent, y3_SCALE_up_percent, y3_SCALE_down_percent, y_err_up3, y_err_down3);
    calculate_errors(y4, y4_PDF_up_percent,y4_PDF_down_percent, y4_SCALE_up_percent, y4_SCALE_down_percent, y_err_up4, y_err_down4);

    // Ratios 1/2, 3/4, 5/6
    std::vector<double> ratio12(x.size()), ratio_err_up12(x.size()), ratio_err_down12(x.size());
    std::vector<double> ratio34(x.size()), ratio_err_up34(x.size()), ratio_err_down34(x.size());

    calculate_ratio(y1, y2, y_err_up1, y_err_up2, y_err_down1, y_err_down2, ratio12, ratio_err_up12, ratio_err_down12);
    calculate_ratio(y3, y4, y_err_up3, y_err_up4, y_err_down3, y_err_down4, ratio34, ratio_err_up34, ratio_err_down34);

    // Création des graphes de ratios
    TGraphAsymmErrors *gr_ratio12 = new TGraphAsymmErrors(x.size(), &x[0], &ratio12[0], nullptr, nullptr, &ratio_err_down12[0], &ratio_err_up12[0]);
    TGraphAsymmErrors *gr_ratio34 = new TGraphAsymmErrors(x.size(), &x[0], &ratio34[0], nullptr, nullptr, &ratio_err_down34[0], &ratio_err_up34[0]);

        Float_t r1 = 0.246;
Float_t g1 = 0.563;
Float_t b1 = 0.852;
TColor color1 = TColor(301,r1, g1, b1);
// color1.SetRGB(r1, g1, b1);
Int_t ColorBlue = color1.GetNumber();

Float_t r2 = 1.000;
Float_t g2 = 0.661;
Float_t b2 = 0.055;
TColor color2 = TColor(302,r2, g2, b2);
// color2.SetRGB(r2, g2, b2);
Int_t ColorOrange = color2.GetNumber();

Float_t r3 = 0.739;
Float_t g3 = 0.122;
Float_t b3 = 0.004;
TColor color3 = TColor(303,r3, g3, b3);
Int_t ColorRed = color3.GetNumber();

Float_t r4 = 0.578;
Float_t g4 = 0.641;
Float_t b4 = 0.635;
TColor color4 = TColor(304,r4, g4, b4);
Int_t ColorGrey = color4.GetNumber();

Float_t r5 = 0.513;
Float_t g5 = 0.176;
Float_t b5 = 0.713;
TColor color5 = TColor(305,r5, g5, b5);
Int_t ColorDarkPurple = color5.GetNumber();

Float_t r6 = 0.661;
Float_t g6 = 0.418;
Float_t b6 = 0.348;
TColor color6 = TColor(306,r6, g6, b6);
Int_t ColorBrown = color6.GetNumber();

Float_t r7 = 0.905;
Float_t g7 = 0.387;
Float_t b7 = 0.000;
TColor color7 = TColor(307,r7, g7, b7);
Int_t ColorDarkOrange = color7.GetNumber();

Float_t r8 = 0.723;
Float_t g8 = 0.672;
Float_t b8 = 0.438;
TColor color8 = TColor(308,r8, g8, b8);
Int_t ColorNeutral = color8.GetNumber();

Float_t r9 = 0.441;
Float_t g9 = 0.457;
Float_t b9 = 0.504;
TColor color9 = TColor(309,r9, g9, b9);
Int_t ColorDarkGrey = color9.GetNumber();

Float_t r10 = 0.571;
Float_t g10 = 0.852;
Float_t b10 = 0.867;
TColor color10 = TColor(310,r10, g10, b10);
// color10.SetRGB(r10, g10, b10);
Int_t ColorLightBlue = color10.GetNumber();

    gr_ratio12->SetMarkerStyle(20);
    gr_ratio12->SetMarkerColor(ColorBlue);
    gr_ratio12->SetLineColor(ColorBlue);

    gr_ratio34->SetMarkerStyle(21);
    gr_ratio34->SetMarkerColor(ColorRed);
    gr_ratio34->SetLineColor(ColorRed);


    gStyle->SetOptDate(0);
    gStyle->SetStatColor(0);
    gStyle->SetTitleFont(62);
    gStyle->SetTitleColor(1);
    gStyle->SetTitleTextColor(1);
    gStyle->SetTitleFillColor(10);
    gStyle->SetTitleFontSize(0.06);
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
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);




    // Canvas et légende
    TCanvas *c = new TCanvas("c", "Ratios", 800, 600);
    gr_ratio12->Draw("AP");
    gr_ratio34->Draw("P SAME");

    gr_ratio12->GetXaxis()->SetTitle("M_{#tilde{#mu}} [GeV]");
    gr_ratio12->GetYaxis()->SetTitle("Ratio");
    gr_ratio12->SetTitle("");
    gr_ratio12->GetYaxis()->SetRangeUser(0.5, 2.8);

    gr_ratio12->GetXaxis()->SetTitleSize(0.06);
    gr_ratio12->GetYaxis()->SetTitleSize(0.06);
    gr_ratio12->GetYaxis()->SetTitleOffset(1);
    // Suite de la macro à partir de TLegend
    TLegend *legend = new TLegend(0.55, 0.7, 0.9, 0.85);
        legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetFillColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(gr_ratio12, "LO+1J/LO", "p");
    legend->AddEntry(gr_ratio34, "NLO+NLL/LO B.Fucks", "p");
    legend->Draw();

    // Affichage et sauvegarde du canvas
    c->SaveAs("plot_XSRR_KFactor_AN.pdf");
    // c->SaveAs("plot_XSLL_KFactor_AN.png");

    // Nettoyage de la mémoire
    delete gr_ratio12;
    delete gr_ratio34;
    delete legend;
    delete c;
}
