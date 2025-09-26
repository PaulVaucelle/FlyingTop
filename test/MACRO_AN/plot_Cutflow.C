#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <string>

void plotCutflow(const std::vector<TString>& fileNames,const std::vector<TString>& Names) {

    int stati=0;
    bool fit= 0;
    bool logy=0;

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

    // Couleurs pour les histogrammes
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

    std::vector<int> colors = {
        ColorBlue,
        ColorOrange,
        ColorRed,
        ColorGrey,
        ColorDarkPurple,
        ColorBrown,
        ColorDarkOrange,
        ColorNeutral,
        ColorDarkGrey,
        ColorLightBlue
    };

    // Create a canvas to plot the histograms
    TCanvas* canvas = new TCanvas("canvas", "Step Efficiency", 800, 600);
    canvas->SetFillColor(10);
    canvas->SetFillStyle(4000);
    canvas->SetBorderSize(2);

    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    // Create a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetFillColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);

    // Prepare to hold histograms
    std::vector<TH1*> histograms;
    std::vector<TH1*> histogramsNorm;
    std::vector<TString> binLabels = {"Step 1", "Step 2", "Step 3", "Step 4", "Step 5","step 6","step7","step8","step9"};

    // Loop over the input file names
    for (size_t i = 0; i < fileNames.size(); ++i) {
        // Open the ROOT file
        TFile* file = TFile::Open("../Signal_2018/histofile_DM_OS_2p4_"+fileNames[i]+".root", "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Cannot open file " << fileNames[i] << std::endl;
            continue;
        }

        // Get the histogram
        TString HISTO = fileNames[i]+"_hData_StepEff_";
        TString HISTONorm = fileNames[i]+"_hData_StepEff_NonNormalized";
        TH1* histo = dynamic_cast<TH1*>(file->Get(HISTO));
        if (!histo) {
            std::cerr << "Error: Histogram "<<HISTO <<" not found in " << fileNames[i] << std::endl;
            file->Close();
            continue;
        }

        TH1* histo2 = dynamic_cast<TH1*>(file->Get(HISTO2));
        if (!histo2) {
            std::cerr << "Error: Histogram "<<HISTO2 <<" not found in " << fileNames[i] << std::endl;
            file->Close();
            continue;
        }

        // Clone the histogram to avoid issues when closing the file
        histograms.push_back((TH1*)histo->Clone());
        histogramsNorm.push_back((TH1*)histo2->Clone());

        histograms.back()->SetDirectory(0); // Detach from the file
        histogramsNorm.back()->SetDirectory(0); // Detach from the file

        // Close the file
        file->Close();

        // Style the histogram
        histograms.back()->SetLineColor(colors[i]); // Use different colors for each histogram
        histograms.back()->SetLineWidth(2);


        // Add to legend
        legend->AddEntry(histograms.back(), Names[i], "l");

        for (int bin = 1; bin <= histograms.back()->GetNbinsX(); ++bin) {
            double content = histograms2.back()->GetBinContent(bin);
            double customError = sqrt(content); // Exemple : définir une erreur comme sqrt(content)
            histograms.back()->SetBinError(bin, customError);
        }

        // Draw the histogram
        if (i == 0) {
            histograms.back()->Draw("E HIST"); // Draw the first histogram : Ajouter option E pour les erreurs
            histograms.back()->SetTitle("");
        } else {
            histograms.back()->Draw("E HIST SAME"); // Overlay the rest : Ajouter option E pour les erreurs
        }
    }

    for (int bin = 1; bin <= histograms[0]->GetNbinsX() && bin <= binLabels.size(); ++bin) {
        histograms[0]->GetXaxis()->SetBinLabel(bin, binLabels[bin - 1]);
    }
    histograms[0]->GetXaxis()->SetRangeUser(0,9);
    histograms[0]->GetYaxis()->SetRangeUser(0,1.2);
    // Draw the legend
    legend->Draw();


//------Start of Copy Paste
TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = true;
TString extraText   = "Simulation";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "";//137 fb^{-1}
TString lumi_sqrtS = "2018";
TString lumiText = lumi_13TeV+lumi_sqrtS;
  float H = canvas->GetWh();
  float W = canvas->GetWw();
  float l = canvas->GetLeftMargin();
  float t = canvas->GetTopMargin();
  float r = canvas->GetRightMargin();
  float b = canvas->GetBottomMargin();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);

float posX_=0;
  float posY_=0;
  int iPosX = 3;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
   posY_ = 1-t - relPosY*(1-t-b);
  	  if( writeExtraText ) 
	    {
         posX_ =   l +  relPosX*(1-l-r);
         posY_ =   1-t+lumiTextOffset*t;
        int alignY_=3;
         int alignX_=2;
         if( iPosX/10==0 ) alignX_=1;
         if( iPosX==0    ) alignX_=1;
         if( iPosX==0    ) alignY_=1;
         if( iPosX/10==1 ) alignX_=1;
         if( iPosX/10==2 ) alignX_=2;
         if( iPosX/10==3 ) alignX_=3;
         //if( iPosX == 0  ) relPosX = 0.12;
         int align_ = 10*alignX_ + alignY_;
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(11);
      latex.DrawLatex(posX_+0.08, posY_, extraText);
	    }
    canvas->Update();
    // Sauvegarder le canvas (optionnel)
    // canvas->SaveAs("output.pdf");
    // canvas->SaveAs(histName+extraTXT+".pdf");


    // Save the canvas to a file
    canvas->SaveAs("plot_Cutflow.pdf");

    // Clean up
    for (auto hist : histograms) {
        delete hist;
    }
    delete legend;
    delete canvas;
}
