#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCD_2Vtx_DataMC_SYST() {
    // Load the macro
    // gROOT->LoadMacro("plot4_ABCD_2Vtx_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    // !! Parameters to be changed !! //
    TString Name = "ABCD_EMU_";
    TString Prod = "DATA_EMU_2018_19_08_2024";
    
    TString YEAR = "2018";
    TString DMODE = "EM";
    bool DATA = false ;
    bool MC = true;
    bool SIGNAL = false;
    TString Plots[6] = {"VtxMass_2Vtx","VtxMass_2VtxAll","STW_2Vtx","STW_2VtxAll","LT_2Vtx","LT_2VtxAll"};
    TString type = ".pdf";

    // !! -------------------------!! //

    int mixing = 0;
    for (int i =0; i < 6; ++i) { //87 
        c1 = plot(i, Prod, Name, YEAR,DMODE , Plots[i], DATA, MC,SIGNAL); // Call the plot function from the macro with argument i
        TString name = Name+Plots[i]+DMODE;
        if (DATA && !MC) {name += "_DATA_SYST";}
        if (MC && !DATA) {name += "_MC_SYST";}
        if (MC && DATA) {name += "_DATAMC_SYST";}
        name += type;
        c1->SaveAs("./"+name);
    }
}