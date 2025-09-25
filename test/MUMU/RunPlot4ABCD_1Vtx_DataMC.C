#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCD_1Vtx_DataMC() {
    // Load the macro
    // gROOT->LoadMacro("plot4_ABCD_1Vtx_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    // !! Parameters to be changed !! //
    TString Name = "ABCD_";
    TString Prod = "DATA_MUMU_2018_19_08_2024";
    TString YEAR = "2018";
    TString DMODE = "DM";
    bool DATA = true ;
    bool MC = false;
    bool SIGNAL = true;
    TString Plots[3] = {"VtxMass_1Vtx_","STW_1Vtx_","LT_1Vtx"};
    TString type = ".pdf";
    bool blind = true;

    // !! -------------------------!! //

    int mixing = 0;
    for (int i =0; i < 3; ++i) { //87
        c1 = plot(i, Prod, Name, YEAR,DMODE , blind, Plots[i],DATA,MC,SIGNAL); // Call the plot function from the macro with argument i
        TString name = Name+Plots[i]+DMODE;
        if (DATA && !MC) {name += "_DATA";}
        if (MC && !DATA) {name += "_MC";}
        if (MC && DATA) {name += "_DATAMC";}
        if (SIGNAL) {name += "_SIGNAL";}
        name += type;
        c1->SaveAs("./"+name);
    }
}