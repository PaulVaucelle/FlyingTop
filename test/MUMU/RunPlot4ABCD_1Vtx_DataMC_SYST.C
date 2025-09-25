#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCD_1Vtx_DataMC_SYST() {
    // Load the macro
    // gROOT->LoadMacro("plot4_ABCD_1Vtx_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    // !! Parameters to be changed !! //
    TString Name = "ABCD_";
    TString Prod = "DATA_MUMU_2018_19_08_2024";//pour Emu DATAMC2018_EMU_10_06_2024 //Mumu PROD_CSI_10_06_2024 //DATA_MUMU_2018_20_06_2024
    TString YEAR = "2018";
    TString DMODE = "DM";
    bool DATA = true ;
    bool MC = false;
    bool SIGNAL = false;
    TString Plots[3] = {"VtxMass_1Vtx_","STW_1Vtx_","LT_1Vtx"};
    TString type = ".pdf";
    bool blind = true;

    // !! -------------------------!! //

    int mixing = 0;
    for (int i =0; i < 3; ++i) { //87
        c1 = plot(i, Prod, Name, YEAR,DMODE , blind, Plots[i],DATA,MC,SIGNAL); // Call the plot function from the macro with argument i
        TString name = Name+Plots[i]+DMODE+type;
        c1->SaveAs("./"+name);
    }
}