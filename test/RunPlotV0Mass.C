#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlotV0Mass() {
    // Load the macro
    gROOT->LoadMacro("plotV0Mass.C");
    // Loop to call the plot function with different integer arguments
    TString YEAR[3] = {"2016","2017","2018"};
    TString V0[2] = {"K0","L0"};


    for (unsigned int i = 2; i < 3; ++i) {
        for (unsigned int j = 0; j < 2; ++j) {

            plot(YEAR[i], V0[j]);

        }
    }

}

