#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot_2DVar() {
    std::vector<TString> fileNames;
    std::vector<TString> Names;
    fileNames.push_back("../Signal_2018/TrackAna_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
    fileNames.push_back("../Signal_2018/TrackAna_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu300_neu280.root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu400_neu300.root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_smu500_neu300.root");
    Names.push_back("DYM50");
    Names.push_back("t#bar{t}");
    Names.push_back("M_{#tilde{#mu}} = 300 GeV, M_{#tilde{#chi}} = 280 GeV");
    Names.push_back("M_{#tilde{#mu}} = 400 GeV, M_{#tilde{#chi}} = 300 GeV");
    Names.push_back("M_{#tilde{#mu}} = 500 GeV, M_{#tilde{#chi}} = 300 GeV");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau001.root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau003.root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau010.root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau030.root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau100.root");
    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau300.root");

    Names.push_back("c#tau = 0.1 cm");
    Names.push_back("c#tau = 0.3 cm");
    Names.push_back("c#tau = 1.0 cm");
    Names.push_back("c#tau = 3.0 cm");
    Names.push_back("c#tau = 10.0 cm");
    Names.push_back("c#tau = 30.0 cm");

    fileNames.push_back("../Signal_2018/TrackAna_RPV_2018_ctau1000.root");
    Names.push_back("c#tau = 100.0 cm");

    for (unsigned int i = 0 ; i < fileNames.size();  i++)
        {
            for (unsigned int j = 0; j < 2; ++j) {//Var : 2

                plot2D(j, fileNames[i], Names[i], i);
            }
        }

}