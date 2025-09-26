#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <vector>

void plot() {
    std::vector<TString> filenames = {
        "TrackAna_RPV_2018_ctau001.root", "TrackAna_RPV_2018_ctau003.root", "TrackAna_RPV_2018_ctau010.root", "TrackAna_RPV_2018_ctau030.root", 
        "TrackAna_RPV_2018_ctau100.root", "TrackAna_RPV_2018_ctau300.root", "TrackAna_RPV_2018_ctau1000.root"
    };
    std::vector<TString> NAMES = {
        "0.1", "0.3", "1.0", "3.0", 
        "10.0", "30.0", "100.0"
    };


    TString histoname = "tree_LLPratio"; // Nom du TH1F dans chaque fichier

    TCanvas *c1 = new TCanvas("c1", "Histogram Means", 800, 600);
    TLegend *leg = new TLegend(0.1, 0.6, 0.3, 0.9);
    leg->SetTextSize(0.04);
    leg->SetHeader("c#tau_{#chi} [cm]");
    bool first = true;
    
    for (size_t i = 0; i < filenames.size(); i++) {
        TFile *file = TFile::Open(filenames[i], "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filenames[i] << std::endl;
            continue;
        }
        
        TH1F *hist = (TH1F*)file->Get(histoname);
        if (!hist) {
            std::cerr << "Histogram " << histoname << " not found in " << filenames[i] << std::endl;
            file->Close();
            continue;
        }
        
        double mean = hist->GetMean();
        double rms = hist->GetMeanError();
        std::cout << "Mean of " << filenames[i] << ": " << mean << std::endl;
        
        TH1F *meanHist = new TH1F("meanHist","Fraction of Vetoed LLP tracks", 7, 0, 7);
        meanHist->SetBinContent(i+1, mean);
        meanHist->SetBinError(i+1, rms);
        meanHist->SetMarkerStyle(20);
        meanHist->SetMarkerColor(i + 1);
        meanHist->SetStats(0);

        if (first) {
            meanHist->Draw("PE1");
            meanHist->GetYaxis()->SetRangeUser(0, 0.25);
            first = false;
        } else {
            meanHist->Draw("PE1 SAME");
        }

        leg->AddEntry(meanHist,  NAMES[i], "p");
        c1->Update();

        // file->Close();
    }

    leg->Draw();

    c1->SaveAs("PlotMeanVetoTrack.pdf");
}
