#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlotVtxEffi() {
    // Load the macro
    // gROOT->LoadMacro("plot_Var.C");

    TString CTAU [7] = {"001","003","010","030","100","300","1000"};

    for (unsigned int n = 0; n < 7; ++n) {// 7 
                plot(CTAU[n]);
    }

}