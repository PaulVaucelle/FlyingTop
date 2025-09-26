#include <iostream>
#include <TROOT.h>
#include "TVectorD.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TEfficiency.h" 
#include "TMath.h" 


void plot(int choice , int VAR , TString CT) {
    // Liste des fichiers ROOT
    std::vector<TString> fileNames;
    std::vector<TString> Names;
    int stati=0;
    bool fit= 0;
    bool logy=0;

    TString extraTXT = "";
    TString Legheader = "";
    int CHOICE = choice;
    int var = VAR;
    TString CTAU =  CT;
    int Channel = 2;// 1 = EMu and 2 = MuMu
    double yMin = 0.0;
    float fac = 1.2;
    if (CHOICE == 0)
        {
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau001.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau003.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau010.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau030.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau100.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau300.root");
            Names.push_back("c#tau = 0.1 cm");
            Names.push_back("c#tau = 0.3 cm");
            Names.push_back("c#tau = 1.0 cm");
            Names.push_back("c#tau = 10.0 cm");
            Names.push_back("c#tau = 30.0 cm");
            Names.push_back("c#tau = 100.0 cm");
            extraTXT = "_ctau_";
        }
    else if (CHOICE == 1)
        {
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu200_neu180.root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu250_neu180.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu300_neu180.root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu350_neu180.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu400_neu180.root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu450_neu180.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu500_neu180.root");
            Names.push_back("M_{#tilde{#mu}} = 200 GeV, M_{#tilde{#chi}} = 180 GeV");
            // Names.push_back("M_{#tilde{#mu}} = 250 GeV, M_{#tilde{#chi}} = 180 GeV");
            Names.push_back("M_{#tilde{#mu}} = 300 GeV, M_{#tilde{#chi}} = 180 GeV");
            // Names.push_back( "M_{#tilde{#mu}} = 350 GeV, M_{#tilde{#chi}} = 180 GeV");
            Names.push_back("M_{#tilde{#mu}} = 400 GeV, M_{#tilde{#chi}} = 180 GeV");
            // Names.push_back("M_{#tilde{#mu}} = 450 GeV, M_{#tilde{#chi}} = 180 GeV");
            Names.push_back( "M_{#tilde{#mu}} = 500 GeV, M_{#tilde{#chi}} = 180 GeV");  
            Legheader= "M_{#tilde{#chi}} = 180 GeV";
            extraTXT = "_MNEU180_ctau_";
        }
    else if (CHOICE == 2)
        {
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu200_neu180_ctau"+CTAU+".root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu250_neu180_ctau"+CTAU+".root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu300_neu180_ctau"+CTAU+".root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu350_neu180_ctau"+CTAU+".root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu400_neu180_ctau"+CTAU+".root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu450_neu180_ctau"+CTAU+".root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu500_neu180_ctau"+CTAU+".root");
            Names.push_back("M_{#tilde{#mu}} = 200 GeV, M_{#tilde{#chi}} = 180 GeV");
            // Names.push_back("M_{#tilde{#mu}} = 250 GeV, M_{#tilde{#chi}} = 180 GeV");
            Names.push_back("M_{#tilde{#mu}} = 300 GeV, M_{#tilde{#chi}} = 180 GeV");
            // Names.push_back("M_{#tilde{#mu}} = 350 GeV, M_{#tilde{#chi}} = 180 GeV");
            Names.push_back("M_{#tilde{#mu}} = 400 GeV, M_{#tilde{#chi}} = 180 GeV");
            // Names.push_back("M_{#tilde{#mu}} = 450 GeV, M_{#tilde{#chi}} = 180 GeV");
            Names.push_back("M_{#tilde{#mu}} = 500 GeV, M_{#tilde{#chi}} = 180 GeV");
            Legheader="c#tau = 0.1 cm";
            extraTXT = "_MNEU180_ctau"+CTAU+"_";
        }
    else if  (CHOICE == 3) 
        {
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu200_neu180.root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu250_neu230.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu300_neu280.root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu350_neu330.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu400_neu380.root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu450_neu430.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu500_neu480.root");
            Names.push_back("M_{#tilde{#mu}} = 200 GeV, M_{#tilde{#chi}} = 180 GeV");
            // Names.push_back("M_{#tilde{#mu}} = 250 GeV, M_{#tilde{#chi}} = 230 GeV");
            Names.push_back("M_{#tilde{#mu}} = 300 GeV, M_{#tilde{#chi}} = 280 GeV");
            // Names.push_back( "M_{#tilde{#mu}} = 350 GeV, M_{#tilde{#chi}} = 330 GeV");
            Names.push_back("M_{#tilde{#mu}} = 400 GeV, M_{#tilde{#chi}} = 380 GeV");
            // Names.push_back("M_{#tilde{#mu}} = 450 GeV, M_{#tilde{#chi}} = 430 GeV");
            Names.push_back( "M_{#tilde{#mu}} = 500 GeV, M_{#tilde{#chi}} = 480 GeV"); 
            Legheader="#Delta M_{#tilde{#mu}-#tilde{#chi}} = 20 GeV";
            extraTXT = "_DM20_ctau_";
        }
    else if  (CHOICE == 4) 
        {
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu200_neu180_ctau"+CTAU+".root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu250_neu230_ctau"+CTAU+".root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu300_neu280_ctau"+CTAU+".root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu350_neu330_ctau"+CTAU+".root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu400_neu380_ctau"+CTAU+".root");
            // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu450_neu430_ctau"+CTAU+".root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu500_neu480_ctau"+CTAU+".root");
            Names.push_back("M_{#tilde{#mu}} = 200 GeV, M_{#tilde{#chi}} = 180 GeV");
            // Names.push_back("M_{#tilde{#mu}} = 250 GeV, M_{#tilde{#chi}} = 230 GeV");
            Names.push_back("M_{#tilde{#mu}} = 300 GeV, M_{#tilde{#chi}} = 280 GeV");
            // Names.push_back( "M_{#tilde{#mu}} = 350 GeV, M_{#tilde{#chi}} = 330 GeV");
            Names.push_back("M_{#tilde{#mu}} = 400 GeV, M_{#tilde{#chi}} = 380 GeV");
            // Names.push_back("M_{#tilde{#mu}} = 450 GeV, M_{#tilde{#chi}} = 430 GeV");
            Names.push_back( "M_{#tilde{#mu}} = 500 GeV, M_{#tilde{#chi}} = 480 GeV"); 
            Legheader= "#Delta M_{#tilde{#mu}-#tilde{#chi}} = 20 GeV, c#tau = 0.1 cm";
            extraTXT = "_DM20_ctau"+CTAU+"_";
        }
    else if  (CHOICE == 5) 
        {
            if (Channel == 1)
                {
                    fileNames.push_back("../Signal_2018/TrackAna_EMU_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
                    fileNames.push_back("../Signal_2018/TrackAna_EMU_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");

                }
            else if(Channel == 2)
                {
                    fileNames.push_back("../Signal_2018/TrackAna_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
                    fileNames.push_back("../Signal_2018/TrackAna_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");

                }
            else // Mumu by default
                {
                    fileNames.push_back("../Signal_2018/TrackAna_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
                    fileNames.push_back("../Signal_2018/TrackAna_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");

                }
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu300_neu280.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu400_neu300.root");
            fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu500_neu300.root");
            Names.push_back("DYM50");
            Names.push_back("t#bar{t}");
            Names.push_back("M_{#tilde{#mu}} = 300 GeV, M_{#tilde{#chi}} = 280 GeV");
            Names.push_back("M_{#tilde{#mu}} = 400 GeV, M_{#tilde{#chi}} = 300 GeV");
            Names.push_back("M_{#tilde{#mu}} = 500 GeV, M_{#tilde{#chi}} = 300 GeV");
            extraTXT = "_SvsBKG_";
        }

    TString ytitle = "entries";
    TString xtitle="p^{j1}_{T} [GeV]";
    // Nom de l'histogramme à récupérer dans chaque fichier
    // Variable à modifier si nécessaire
    // const char* histName = "myHistogram";
    TString histName = "hData_jet_leadingpt";

    if (var == 0)
        {
            ytitle = "a.u";
            xtitle="p^{j1}_{T} [GeV]";
            histName = "hData_jet_leadingpt";
        }
    else if (var == 1)
        {
            ytitle = "a.u";
            xtitle="p^{j2}_{T} [GeV]";
            histName = "hData_jet_leadingpt2";
        }
    else if (var == 2)
        {
            ytitle = "a.u";
            xtitle="p^{#mu 1}_{T} [GeV]";
            histName = "hData_reco_lepton_leadingpt";
        }
    else if (var == 3)
        {
            ytitle = "a.u";
            xtitle="p^{#mu 2}_{T} [GeV]";
            histName = "hData_reco_lepton_leadingpt2";
        }
    else if (var == 4)
        {
            ytitle = "a.u";
            xtitle="#Delta #phi_{#mu#mu}";
            histName = "hData_lepton_lepton_dPhi";
        }
    else if (var == 5)
        {
            ytitle = "a.u";
            xtitle="#Delta #eta_{#mu#mu}";
            histName = "hData_lepton_lepton_dEta";
        }
    else if (var == 6)
        {
            ytitle = "a.u";
            xtitle="#Delta R_{#mu#mu}";
            histName = "hData_lepton_lepton_dR";
        }
    else if (var == 7)
        {
            return;
        }
    else if (var == 8)
        {
            ytitle = "a.u";
            xtitle="#Delta R_{jet_{1-2}}";
            histName = "hData_jet_jet_dR";
        }
    else if (var == 9)  
        {
            ytitle = "a.u";
            xtitle="trk p_{T} [GeV]";
            histName = "tree_track_pt";
        }
    else if (var == 10)  
        {
            ytitle = "a.u";
            xtitle="trk #eta ";
            histName = "tree_track_eta";
            fac= 1.8;
        }
    else if (var == 11)  
        {
            ytitle = "a.u";
            xtitle="trk nHits ";
            histName = "tree_track_nHit";
        }
    else if (var == 12)  
        {
            ytitle = "a.u";
            xtitle="nTrk10 ";
            histName = "tree_track_ntrk10";
        }
    else if (var == 13)  
        {
            ytitle = "a.u";
            xtitle="nTrk20 ";
            histName = "tree_track_ntrk20";
        }
    else if (var == 14)  
        {
            ytitle = "a.u";
            xtitle="nTrk30 ";
            histName = "tree_track_ntrk30";
        }
    else if (var == 15)  
        {
            ytitle = "a.u";
            xtitle="nTrk40 ";
            histName = "tree_track_ntrk40";
        }
    else if (var == 16)  
        {
            ytitle = "a.u";
            xtitle="Trk Sig_{dz} ";
            histName = "tree_track_dzSig";
        }
    else if (var == 17)  
        {
            ytitle = "a.u";
            xtitle="Trk Sig_{dr} ";
            histName = "tree_track_drSig";
        }
    else if (var == 18)  
        {
            ytitle = "a.u";
            xtitle="Trk InJet ";
            histName = "tree_track_isInJet";
        }
    else if (var == 19)  
        {
            ytitle = "a.u";
            xtitle="Hemi-track #Delta R_{min}";
            histName = "tree_track_Hemi_dR";
        }
    else if (var == 20)  
        {
            ytitle = "a.u";
            xtitle="Hemi-track #Delta R_{max}";
            histName = "tree_track_Hemi_dRmax";
        }
    else if (var == 21)  
        {
            ytitle = "a.u";
            xtitle="Trk isLost";
            histName = "tree_track_isLost";
        }
    else if (var == 22)  
        {
            ytitle = "a.u";
            xtitle="TRK BDT Score";
            histName = "tree_track_MVAval";
            logy = 1;
            yMin = 0.0001;
        }
    else if (var == 23)  
        {
            ytitle = "Signal Efficiency";
            xtitle="TRK BDT Score ";
            histName = "hSignalEff";
        }
    else if (var == 24)  
        {
            ytitle = "BKG Efficiency";
            xtitle="TRK BDT Score";
            histName = "hBkgEff";
        }
    else if (var == 25)  
        {
            ytitle = "Significance";
            xtitle="TRK BDT Score";
            histName = "h_sig";
        }
    else if (var == 26)  
        {
            ytitle = "a.u";
            xtitle="Trk Weight";
            histName = "hData_Hemi_Vtx_trackWeight";
            logy = 1;
            yMin = 0.0001;
        }  
    else if (var == 27)  
        {
            ytitle = "a.u";
            xtitle="Sum Trk Weight";
            histName = "hData_Hemi_Vtx_SumtrackWeight";
            logy = 1;
            yMin = 0.0001;
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

    // Vecteur pour stocker les histogrammes
    std::vector<TH1*> histograms;

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

    // Boucle sur les fichiers pour récupérer les histogrammes
    for (size_t i = 0; i < fileNames.size(); ++i) {
        TFile* file = TFile::Open(fileNames[i]);
        if (!file || file->IsZombie()) {
            std::cerr << "Erreur : Impossible d'ouvrir le fichier " << fileNames[i] << std::endl;
            continue;
        }

        TH1* hist = dynamic_cast<TH1*>(file->Get(histName));
        if (!hist) {
            std::cerr << "Erreur : Histogramme '" << histName << "' non trouvé dans " << fileNames[i] << std::endl;
            file->Close();
            continue;
        }

        hist->SetDirectory(0); // Détacher l'histogramme du fichier
        histograms.push_back(hist);

        // Appliquer une couleur à l'histogramme

        hist->SetLineColor(colors[i]);

        file->Close();
    }

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
    gPad->SetBottomMargin(0.15);

    
    // Trouver les limites pour l'axe Y
    double yMax = 0.0;
 
    // Dessiner les histogrammes sur le même canvas
    for (size_t i = 0; i < histograms.size(); ++i) {
        
        if (i == 0) {
            
            histograms[i]->SetLineColor(colors[i]);
            histograms[i]->SetLineWidth(2);
            histograms[i]->SetLineStyle(1);
            histograms[i]->SetMarkerStyle(76);
            histograms[i]->SetMarkerSize(1);
            histograms[i]->SetMarkerColor(colors[i]);
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
            histograms[i]->Draw("HIST"); // Le premier histogramme
            yMax = histograms[i]->GetMaximum();
            

        } else {

            histograms[i]->SetLineColor(colors[i]);
            histograms[i]->SetLineWidth(2);
            histograms[i]->SetLineStyle(1);
            histograms[i]->SetMarkerStyle(76);
            histograms[i]->SetMarkerSize(1);
            histograms[i]->SetMarkerColor(colors[i]);
            histograms[i]->Scale(1./histograms[i]->Integral());
            histograms[i]->Draw("HIST SAME"); // Les suivants
            double histMax = histograms[i]->GetMaximum();
            if (histMax > yMax) yMax = histMax;
        }
    }
    histograms[0]->SetMinimum(yMin);
    
    gPad->SetLogy(logy);
    std::cout<<"ymax : _"<<yMax<<std::endl;
    histograms[0]->SetMaximum(fac * yMax); // Ajouter un peu de marge
    // Ajouter une légende

    TLegend* legend;
    if (var == 4 || var == 6 || var == 8 || var == 11  || var == 20)
        {
            legend = new TLegend(0.2, 0.6, 0.5, 0.85);
        } 
    else if ( var == 21 )
        {
            legend = new TLegend(0.55, 0.6, 0.8, 0.85);
        }
    else if ( var == 18)
        {
            legend = new TLegend(0.17, 0.6, 0.45, 0.85);
        }  
    else if ( (var == 23 || var == 25) && choice == 0)
        {
            legend = new TLegend(0.35, 0.3, 0.65, 0.55);
        } 
    else 
        {
            legend = new TLegend(0.5, 0.6, 0.8, 0.85);
        }

    legend->SetHeader(Legheader);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetFillColor(kWhite);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    for (size_t i = 0; i < histograms.size(); ++i) {
        legend->AddEntry(histograms[i], Names[i], "l");
    }
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
    canvas->SaveAs(histName+extraTXT+".pdf");

}
