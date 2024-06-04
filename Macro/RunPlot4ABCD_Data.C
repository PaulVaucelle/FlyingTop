#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCD_Data() {
    // Load the macro
    gROOT->LoadMacro("plot4_ABCD.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    TString Name = "ABCD_DATA_EMU18_095_";
    TString type = ".pdf";
    for (int i = 12; i < 13; ++i) { //85
        c1 = plot(i); // Call the plot function from the macro with argument i
        TString name = Name+Form("%d",i)+type;
        c1->SaveAs(name);
    }
}