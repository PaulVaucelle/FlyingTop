#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot_Var() {
    // Load the macro
    // gROOT->LoadMacro("plot_Var.C");
    // Loop to call the plot function with different integer arguments
    // TString YEAR[3] = {"2016","2017","2018"};
    // TString V0[2] = {"K0","L0"};
    TString CTAU [7] = {"001","003","010","030","100","300","1000"};

    for (unsigned int n = 0; n < 1; ++n) {// 7 
        for (unsigned int i = 0; i < 5; ++i) {//CHoice
            for (unsigned int j = 0; j < 2; ++j) {//Var : 26

                plot(i,j, CTAU[n]);

            }
        }
    }

}