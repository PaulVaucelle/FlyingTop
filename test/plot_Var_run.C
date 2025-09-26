#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void plot_Var_run() {
    // Load the macro
    // gROOT->LoadMacro("plot4_ABCD_2Vtx_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    // !! Parameters to be changed !! //
    TString Name = "Plot_";
    TString YEAR = "2018";
    bool SIGNAL = true;
    TString Plots[2] = {
    "HemiLeadingPt",
    "HemiSubLeadingPt",
    };
   

    // !! -------------------------!! //

    for (int i =0; i < 2; ++i) { 
        // Plots[i]+= "_Corr";

        c1 = plot(i, Name, YEAR, Plots[i], SIGNAL);
        // TString name = Name+Plots[i]+DMODE;
        // if (SIGNAL) {name += "_SIGNAL";}
        // name += type;
        // c1->SaveAs("./"+name);
    }
}