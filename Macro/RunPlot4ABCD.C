#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCD() {
    // Load the macro
    gROOT->LoadMacro("plot4_ABCD_solve.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    TString Name = "ABCD_MC_MuMu_";
    TString type = ".pdf";
    for (int i = 0; i < 60; ++i) { //85
        c1 = plot4_ABCD_solve(i); // Call the plot function from the macro with argument i
        TString name = Name+Form("%d",i)+type;
        c1->SaveAs(name);
    }
}