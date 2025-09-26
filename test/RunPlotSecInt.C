#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlotSecInt() {
    // Load the macro
    gROOT->LoadMacro("plotSecInt.C");
    // Loop to call the plot function with different integer arguments
    TString YEAR[4] = {"2016POST","2016PRE","2017","2018"};
    TString PLAN[4] = {"xy","rz","xy_Inner","rz_Inner"};
    TString SELECTION[2] = {"Selec","TrackerMatched"};
    TString PU[7] = {"","_PU25","_PU30","_PU35","_PU40","_PU45","_PU50"};
    TString Sample[4]={"DoubleMuon_UL2016POST_MiniAODv2","DoubleMuon_UL2016PRE_MiniAODv2","DoubleMuon_UL2017_MiniAODv2","DoubleMuon_UL2018_MiniAODv2_GT36-v1"};
    // TString Sample[1]={"MCBKG"};//TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8, DoubleMuon_UL2018_MiniAODv2_GT36-v1
    // for (unsigned int u = 0 ; u<4; u++) // Loop on samples
    //     {
            for (unsigned int i = 3; i < 4; ++i) { // loop on the years (relevant only for the data)
                for (unsigned int j = 2; j < 4; ++j) { // loop on the plan (xy, rz, xy_Inner, rz_Inner)
                    for (unsigned int k = 0; k < 2; ++k) { // loop on the selection (Selec, TrackerMatched) 
                        for (unsigned int l = 0; l < 1; ++l) { // loop on the PU (empty, _PU25, _PU30, _PU35, _PU40, _PU45, _PU50)
                            plot(Sample[i],YEAR[i], PLAN[j], SELECTION[k], PU[l]);
                        }
                    }
                }
            }
        // }


}

