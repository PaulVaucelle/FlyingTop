#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCDEFGHI_2Vtx_DataMC() {
    // Load the macro
    gROOT->LoadMacro("plot4_ABCDEFGHI_2Vtx_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    TString Name = "DATAMC_ABCDEFGHI_EMU_";
    TString YEAR = "2018";
    TString Prod = "DATA_EMU_2018_19_08_2024";
    //RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1
    //RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9
    // RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17
    // RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11
    TString Dmode = "EM";
    TString type = ".pdf";
        TString Plots[4] = {"VtxMass_2Vtx","VtxMass_2VtxAll","STW_2Vtx","STW_2VtxAll"};

    for (int i =0; i < 4  ; ++i) { 
        c1 = plot(i, Prod, Name, YEAR,Dmode, Plots[i] ); // Call the plot function from the macro with argument i
        TString name = Name+Plots[i]+Dmode+type;
        c1->SaveAs("./"+name);
    }
}

