#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCDEF_1Vtx_DataMC() {
    // Load the macro
    // gROOT->LoadMacro("plot4_ABCDEF_1Vtx_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    // !! Parameters to be changed !! //
    TString Name = "DATAMC_ABCDEF_EMU_";
    TString Prod = "DATA_EMU_2018_19_08_2024";
    //RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1
    //RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9
    // RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17
    // RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11
    TString YEAR = "2018";
    TString DMODE = "EM";
    TString Plots[2] = {"VtxMass_1Vtx_","STW_1Vtx_"};
    TString type = ".pdf";

    // !! -------------------------!! //

    int mixing = 0;
    for (int i =0; i < 2; ++i) { 
        c1 = plot(i, Prod, Name, YEAR,DMODE ,  Plots[i]); // Call the plot function from the macro with argument i
        TString name = Name+Plots[i]+DMODE+type;
        c1->SaveAs("./"+name);
    }
}