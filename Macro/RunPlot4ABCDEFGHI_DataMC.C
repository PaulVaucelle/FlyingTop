#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCDEFGHI_DataMC() {
    // Load the macro
    gROOT->LoadMacro("plot4_ABCDEFGHI_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    TString Name = "DATAMC_ABCDEFGHI_";
    TString Prod = "DATAMC2018_EMU_10_06_2024";//pour Emu DATAMC2018_EMU_10_06_2024 //Mumu PROD_CSI_10_06_2024
    TString type = ".pdf";
    for (int i =13; i < 14  ; ++i) { //89
        c1 = plot(i, Prod, Name); // Call the plot function from the macro with argument i
        TString name = Name+Form("%d",i)+type;
        c1->SaveAs("../"+Prod+"/"+name);
    }
}