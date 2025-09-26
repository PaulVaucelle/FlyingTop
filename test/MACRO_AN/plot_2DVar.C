#include <iostream>
#include <TROOT.h>
#include "TVectorD.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TEfficiency.h" 
#include "TMath.h" 
#include "string.h"

void plot2D( int VAR , TString fileName, TString NAME, int SaveInt) {
    // Liste des fichiers ROOT

    int stati=0;
    bool fit= 0;
    bool logy=0;
    bool logz=0;

    TString extraTXT = "_";
    TString Legheader = "";

    int var = VAR;
    TString s = std::to_string(SaveInt);
    TString ytitle = "entries";
    TString xtitle="p^{j1}_{T} [GeV]";
    TString ztitle="a.u";
    // Nom de l'histogramme à récupérer dans chaque fichier
    // Variable à modifier si nécessaire
    // const char* histName = "myHistogram";
    TString histName = "hBkgEff_SigEff";

    if (var == 0)
        {
            ztitle = "a.u";
            ytitle = "Bkg Efficiency";
            xtitle="Signal Efficiency";
            histName = "hBkgEff_SigEff";
        }
    else if (var == 1)
        {
            ytitle = "STW";
            xtitle="nTrks";
            ztitle="a.u";
            histName = "hData_Hemi_Vtx_SumtrackWeight_vs_Ntrks";
            logz=1;
        }

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
    gStyle->SetPalette(kViridis);
    // Couleurs pour les histogrammes


    // Vecteur pour stocker les histogrammes
    std::vector<TH2*> histograms;
    // Boucle sur les fichiers pour récupérer les histogrammes

        TFile* file = TFile::Open(fileName);
        if (!file || file->IsZombie()) {
            std::cerr << "Erreur : Impossible d'ouvrir le fichier " << fileName << std::endl;
            return;
        }

        TH2* hist = dynamic_cast<TH2*>(file->Get(histName));
        if (!hist) {
            std::cerr << "Erreur : Histogramme '" << histName << "' non trouvé dans " << fileName << std::endl;
            file->Close();
            return;
        }

        hist->SetDirectory(0); // Détacher l'histogramme du fichier
        histograms.push_back(hist);

        // Appliquer une couleur à l'histogramme
        // int color = colors[i % colors.size()];
        // hist->SetLineColor(color);

        file->Close();
    

    // Vérifier si au moins un histogramme a été récupéré
    if (histograms.empty()) {
        std::cerr << "Erreur : Aucun histogramme valide n'a été récupéré." << std::endl;
        return;
    }

    // Créer un canvas pour dessiner les histogrammes
    TCanvas* canvas = new TCanvas("canvas", "Var MC Signal vs Bkg", 800, 600);
    canvas->SetFillColor(10);
    canvas->SetFillStyle(4000);
    canvas->SetBorderSize(2);

    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.15);
    gPad->SetBottomMargin(0.15);
    // Trouver les limites pour l'axe Y
    double yMax = 0.0;

    // Dessiner les histogrammes sur le même canvas
    for (size_t i = 0; i < histograms.size(); ++i) { //size is actually 1 :D
        
        if (i == 0) {
            
            // histograms[i]->SetLineColor(colors[i]);
            histograms[i]->SetLineWidth(2);
            histograms[i]->SetLineStyle(1);
            histograms[i]->SetMarkerStyle(76);
            histograms[i]->SetMarkerSize(1);
            // histograms[i]->SetMarkerColor(colors[i]);
            histograms[i]->Scale(1./histograms[i]->Integral());
            histograms[i]->GetXaxis()->SetTitleSize(0.06);
            histograms[i]->GetYaxis()->SetTitleSize(0.06);
            histograms[i]->GetYaxis()->SetTitleOffset(1.5);
            histograms[i]->SetTickLength(0.03, "XYZ");
            histograms[i]->SetLabelOffset(0.002,"X");//0.007
            histograms[i]->SetLabelOffset(0.007,"Y");
            histograms[i]->SetLabelSize(0.042, "XYZ");
            histograms[i]->SetLabelFont(42, "XYZ"); 
            histograms[i]->SetTitleFont(42, "XYZ");
            histograms[i]->SetTitleSize(0.05, "XYZ"); 
            histograms[i]->SetTitleOffset(1.3,"Y");
            histograms[i]->SetTitleOffset(1.1,"X");
            histograms[i]->SetNdivisions(509,"XYZ");
            histograms[i]->SetTitle("");
            histograms[i]->GetXaxis()->SetTitle(xtitle);
            histograms[i]->GetYaxis()->SetTitle(ytitle);
            histograms[i]->GetZaxis()->SetTitle(ztitle);
            histograms[i]->Draw("COLZ"); // Le premier histogramme
            yMax = histograms[i]->GetMaximum();

        } else {

            // histograms[i]->SetLineColor(colors[i]);
            histograms[i]->SetLineWidth(2);
            histograms[i]->SetLineStyle(1);
            histograms[i]->SetMarkerStyle(76);
            histograms[i]->SetMarkerSize(1);
            // histograms[i]->SetMarkerColor(colors[i]);
            histograms[i]->Scale(1./histograms[i]->Integral());
            histograms[i]->Draw("COLZ SAME"); // Les suivants
            double histMax = histograms[i]->GetMaximum();
            if (histMax > yMax) yMax = histMax;
        }
    }
    std::cout<<"ymax : _"<<yMax<<std::endl;
    gPad->SetLogz(logz);
    histograms[0]->SetMaximum(1.2 * yMax); // Ajouter un peu de marge
    // Ajouter une légende
    TLegend* legend = new TLegend(0.1, 0.8, 0.5, 0.9);
    legend->SetHeader(Legheader);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetFillColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry(histograms[0], NAME, "");
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
    canvas->SaveAs(histName+extraTXT+s+".pdf");

}

