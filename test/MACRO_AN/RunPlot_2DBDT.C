#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot_2DBDT() {
    std::vector<TString> fileNames;
    std::vector<TString> Names;


    TString SMUON = "500";
    TString NEU = "180";
    TString CTAU = "300";
    TString legCTAU = "30";

    //---------------------------------------------

    // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu"+SMUON+"_neu"+NEU+"_ctau001.root");
    // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu"+SMUON+"_neu"+NEU+"_ctau003.root");
    // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu"+SMUON+"_neu"+NEU+"_ctau010.root");
    // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu"+SMUON+"_neu"+NEU+"_ctau030.root");
    // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu"+SMUON+"_neu"+NEU+"_ctau100.root");
    // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu"+SMUON+"_neu"+NEU+"_ctau300.root");
    // fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu"+SMUON+"_neu"+NEU+"_ctau1000.root");

    // Names.push_back("M_{#tilde{#mu}} = "+SMUON+" GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = 0.1 cm");
    // Names.push_back("M_{#tilde{#mu}} = "+SMUON+" GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = 0.3 cm");
    // Names.push_back("M_{#tilde{#mu}} = "+SMUON+" GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = 1.0 cm");
    // Names.push_back("M_{#tilde{#mu}} = "+SMUON+" GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = 3.0 cm");
    // Names.push_back("M_{#tilde{#mu}} = "+SMUON+" GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = 10.0 cm");
    // Names.push_back("M_{#tilde{#mu}} = "+SMUON+" GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = 30.0 cm");
    // Names.push_back("M_{#tilde{#mu}} = "+SMUON+" GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = 100.0 cm");

    //---------------------------------------------

    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu200_neu"+NEU+"_ctau"+CTAU+".root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu250_neu"+NEU+"_ctau"+CTAU+".root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu300_neu"+NEU+"_ctau"+CTAU+".root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu350_neu"+NEU+"_ctau"+CTAU+".root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu400_neu"+NEU+"_ctau"+CTAU+".root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu450_neu"+NEU+"_ctau"+CTAU+".root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu500_neu"+NEU+"_ctau"+CTAU+".root");

    Names.push_back("M_{#tilde{#mu}} = 200 GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = "+legCTAU+" cm");
    Names.push_back("M_{#tilde{#mu}} = 250 GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = "+legCTAU+" cm");
    Names.push_back("M_{#tilde{#mu}} = 300 GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = "+legCTAU+" cm");
    Names.push_back("M_{#tilde{#mu}} = 350 GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = "+legCTAU+" cm");
    Names.push_back("M_{#tilde{#mu}} = 400 GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = "+legCTAU+" cm");
    Names.push_back("M_{#tilde{#mu}} = 450 GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = "+legCTAU+" cm");
    Names.push_back("M_{#tilde{#mu}} = 500 GeV, M_{#tilde{#chi}} = "+NEU+" GeV, c#tau = "+legCTAU+" cm");

    plot2D(fileNames, Names);



}