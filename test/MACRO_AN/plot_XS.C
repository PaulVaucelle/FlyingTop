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


void plot() {
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
     std::ofstream ofs ("./plot_XS_extra.txt", std::ofstream::out);

    int stati=0;
    bool fit= 0;
    bool logy=0;

    // Points sur l'axe x (identiques pour tous les graphes)
    std::vector<double> x = {200, 250, 300, 350, 400, 450, 500};

    // 13 TeV Values
    // 6 vecteurs y (résultats fictifs pour illustration)

    std::vector<double> y2 = {8.96,3.78,1.82,0.95,0.53,0.31,0.19}; // LO this analysis

    // Erreurs PDF (symétriques) pour chaque vecteur y en %
    std::vector<double> y2_PDF_up_percent = {1.8,2.1,2.4,2.7,2.9,3.2,3.5};// LO this analysis
    std::vector<double> y2_PDF_down_percent = {1.8,2.1,2.4,2.7,2.9,3.2,3.5};// LO this analysis

    // Erreurs SCALE (asymétriques) pour chaque vecteur y en %

    std::vector<double> y2_SCALE_up_percent = {1.7,2.9,3.9,4.7,5.4,6.0,6.6};// LO this analysis
    std::vector<double> y2_SCALE_down_percent = {1.9,2.9,3.7,4.3,4.9,5.4,5.9};// LO this analysis   

    // Calcul des erreurs asymétriques pour tous les vecteurs y

    std::vector<double> y_err_up2(x.size()), y_err_down2(x.size());
    calculate_errors(y2, y2_PDF_up_percent,y2_PDF_down_percent, y2_SCALE_up_percent, y2_SCALE_down_percent, y_err_up2, y_err_down2);
    
    // Création des graphes de ratios
    TGraphAsymmErrors *gr_ratio12 = new TGraphAsymmErrors(x.size(), &x[0], &y2[0], nullptr, nullptr, &y_err_down2[0], &y_err_up2[0]);


    gr_ratio12->SetMarkerStyle(20);
    gr_ratio12->SetMarkerColor(ColorBlue);
    gr_ratio12->SetLineColor(ColorBlue);

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
    c->cd();
    gr_ratio12->Draw("AP");
    gr_ratio12->GetXaxis()->SetTitle("m_{#tilde{#mu}} [GeV]");
    gr_ratio12->GetYaxis()->SetTitle("#sigma_{#tilde{#mu_{L}}#tilde{#mu}}^{LO + 1J}");
    gr_ratio12->SetTitle("");
    gr_ratio12->GetYaxis()->SetRangeUser(0., 15.);

    gr_ratio12->GetXaxis()->SetTitleSize(0.06);
    gr_ratio12->GetYaxis()->SetTitleSize(0.06);
    gr_ratio12->GetYaxis()->SetTitleOffset(1);
    gr_ratio12->GetXaxis()->SetTitleFont(42);
    gr_ratio12->GetYaxis()->SetTitleFont(42);
    // Suite de la macro à partir de TLegend

    TF1* fitFunc = new TF1("fitFunc", "[0]*exp(-x^{[2]}/[1])", 100, 500);
    fitFunc->SetParameters(1000,200,1);// STW
    gr_ratio12->Fit(fitFunc, "R");
    fitFunc->SetLineColor(kBlack);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("same");

    double a = fitFunc->GetParameter(0);
    double b = fitFunc->GetParameter(1);
    double d = fitFunc->GetParameter(2);
    double aerr = fitFunc->GetParError(0);
    double berr = fitFunc->GetParError(1);
    double derr = fitFunc->GetParError(2);

    TF1 *fitFuncUp = new TF1("fitFuncUp", "[0]*exp(-x^{[2]}/[1])",100 , 500);//htotMC->GetXaxis()->GetXmin()
    fitFuncUp->SetParameter(0, a); // 
    fitFuncUp->SetParameter(1, b+berr); 
     fitFuncUp->SetParameter(2, d+derr); 
    fitFuncUp->SetLineColor(kBlue);
    fitFuncUp->Draw("same");

    TF1 *fitFuncDown = new TF1("fitFuncDown", "[0]*exp(-x^{[2]}/[1])",100 , 500);//htotMC->GetXaxis()->GetXmin()
    fitFuncDown->SetParameter(0, a); //
    fitFuncDown->SetParameter(1, b-berr); //
    fitFuncDown->SetParameter(2, d-derr); //
    fitFuncDown->SetLineColor(kRed);
    fitFuncDown->Draw("same");

    TLegend *legend = new TLegend(0.45, 0.7, 0.9, 0.85);
        legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetFillColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    legend->AddEntry(gr_ratio12, "LO", "p");
    legend->AddEntry(fitFunc, "Fit", "l");
    legend->AddEntry(fitFuncUp, "FitUp", "l");
    legend->AddEntry(fitFuncDown, "FitDown", "l");
    legend->Draw();
    c->SaveAs("plot_XS_fit.pdf");
    // TGraphAsymmErrors* gr = new TGraphAsymmErrors(extraMsmu.size(), &extraMsmu[0], &y2[0], nullptr, nullptr, &y_err_down2[0], &y_err_up2[0]);;
    //  gr->Reset();
    TCanvas *c2 = new TCanvas("c2", "Ratios", 800, 600);
    c2->cd();
    gPad->SetLogy(1);
    // !! - -
    std::vector<double> extraMsmu = { 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
    std::vector<double> extraY2(extraMsmu.size());
    std::vector<double> extraY2Up(extraMsmu.size());
    std::vector<double> extraY2Down(extraMsmu.size());
    for (unsigned int i = 0; i < extraMsmu.size(); i++) {
        double x = extraMsmu[i];
        extraY2[i] = fitFunc->Eval(x);
        ofs<<x<<" GeV /  : XS [fb] : "<<extraY2[i]<<std::endl;
        extraY2Up[i] = fitFuncUp->Eval(x);
        extraY2Down[i] = fitFuncDown->Eval(x);
        ofs<<x<<" GeV /  : XS Up [fb] : "<<extraY2Up[i]<<std::endl;
        ofs<<x<<" GeV /  : XS Down [fb] : "<<extraY2Down[i]<<std::endl;
    }
    ofs.close();
    // !! -- 
    std::vector<double> FullMsmu = {200,250,300,350,400,450,500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};
    std::vector<double> FullY2(FullMsmu.size());
    std::vector<double> FullY2Up(FullMsmu.size());
    std::vector<double> FullY2Down(FullMsmu.size());

    for (unsigned int i = 0; i < x.size() + extraMsmu.size(); i++) {
        if (i < x.size()) {
            FullY2[i] = y2[i];
            FullY2Up[i] = y2[i];
            FullY2Down[i] = y2[i];
        } else {
            FullY2[i] = extraY2[i - x.size()];
            FullY2Up[i] = extraY2Up[i - x.size()];
            FullY2Down[i] = extraY2Down[i - x.size()];
        }
    }
    
    TGraph * GR = new TGraph(FullMsmu.size(), &FullMsmu[0], &FullY2[0]);
    GR->SetMarkerStyle(20);
    GR->SetMarkerColor(kBlack);
    GR->SetLineColor(kBlack);
    GR->SetLineWidth(1);
    GR->SetMarkerSize(1.);
    GR->SetTitle("");
    GR->GetXaxis()->SetTitle("m_{#tilde{#mu}} [GeV]");
    GR->GetYaxis()->SetTitle("#sigma_{#tilde{#mu_{L}}#tilde{#mu}}^{LO + 1J}");
    GR->GetYaxis()->SetRangeUser(0.001, 20.);
    GR->GetXaxis()->SetTitleSize(0.06);
    GR->GetYaxis()->SetTitleSize(0.06);
    GR->GetYaxis()->SetTitleOffset(1);
    GR->GetXaxis()->SetTitleFont(42);
    GR->GetYaxis()->SetTitleFont(42);
    GR->GetXaxis()->SetLabelFont(42);
    GR->GetYaxis()->SetLabelFont(42);
    GR->GetXaxis()->SetLabelSize(0.06);
    GR->GetYaxis()->SetLabelSize(0.06);
    GR->GetXaxis()->SetTickLength(0.03);
    GR->GetYaxis()->SetTickLength(0.03);
    GR->GetXaxis()->SetRangeUser(0, 1100);
    GR->Draw("APsame");

    TGraph * GRUp = new TGraph(FullMsmu.size(), &FullMsmu[0], &FullY2Up[0]);
    GRUp->SetMarkerStyle(20);
    GRUp->SetMarkerColor(ColorRed);
    GRUp->SetLineColor(ColorRed);
    GRUp->SetLineWidth(1);
    GRUp->SetMarkerSize(1.);
    GRUp->Draw("Psame");

    TGraph * GRDown = new TGraph(FullMsmu.size(), &FullMsmu[0], &FullY2Down[0]);
    GRDown->SetMarkerStyle(20);
    GRDown->SetMarkerColor(ColorBlue);
    GRDown->SetLineColor(ColorBlue);
    GRDown->SetLineWidth(1);
    GRDown->SetMarkerSize(1.);
    GRDown->Draw("Psame");

    legend = new TLegend(0.45, 0.7, 0.9, 0.85);
        legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetFillColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    legend->AddEntry(GR, "Extrapolation", "l");
    legend->AddEntry(GRUp, "ExtrapolationUp", "l");
    legend->AddEntry(GRDown, "ExtrapolationDown", "l");
    legend->Draw();
    // Affichage et sauvegarde du canvas
    c2->SaveAs("plot_XS_extra.pdf");
    // c->SaveAs("plot_XSLL_KFactor_AN.png");

    // Nettoyage de la mémoire
    delete gr_ratio12;
    delete legend;
    delete c;
}
