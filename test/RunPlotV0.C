#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlotV0() {
    // Load the macro
    gROOT->LoadMacro("plotV0.C");
    // Loop to call the plot function with different integer arguments
    TString YEAR[3] = {"2016","2017","2018"};
    TString PLAN[4] = {"xy","rz"};
    TString SELECTION[2] = {"reco","CMSSW"};
    TString Channel[2] = {"V0","Yc"};
    // TString Sample[2]={"DoubleMuon_UL2018_MiniAODv2_GT36-v1","TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"};
    TString Sample[1]={"MCBKG"};

    for (unsigned int u = 0 ; u < 1 ; u++) { // Loop on samples
        for (unsigned int j = 0; j < 2; ++j) { // loop on the plan (xy, rz)
            plot(Sample[u],2018, PLAN[j], SELECTION[0], Channel[0]);
            plot(Sample[u],2018, PLAN[j], SELECTION[1], Channel[0]);
            plot(Sample[u],2018, PLAN[j], SELECTION[1], Channel[1]);
        }
    }


//    addHisto2D("hData_reco_V0_xy","", sample.Data(), 500,-5,5,500,-5,5);
//    addHisto2D("hData_reco_V0_rz","", sample.Data(), 1200,0.,120.,700,0.,70.);
//    addHisto2D("hData_CMSSW_V0_xy","", sample.Data(),500,-5,5,500,-5,5);
//    addHisto2D("hData_CMSSW_V0_rz","", sample.Data(), 1200,0.,120.,700,0.,70.);
//    addHisto2D("hData_CMSSW_Yc_xy","", sample.Data(), 500,-5,5,500,-5,5);
//    addHisto2D("hData_CMSSW_Yc_rz","", sample.Data(), 1200,0.,120.,700,0.,70.);


}

