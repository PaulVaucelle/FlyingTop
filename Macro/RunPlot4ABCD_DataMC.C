#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCD_DataMC() {
    // Load the macro
    gROOT->LoadMacro("plot4_ABCD_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    TString Name = "DATAMC_ABCD_";
    TString Prod = "PROD_CSI_10_06_2024";//pour Emu DATAMC2018_EMU_10_06_2024 //Mumu PROD_CSI_10_06_2024 //DATA_MUMU_2018_20_06_2024
    TString type = ".pdf";
    for (int i =56; i < 67; ++i) { //89
        c1 = plot(i, Prod, Name); // Call the plot function from the macro with argument i
        TString name = Name+Form("%d",i)+type;
        c1->SaveAs("../"+Prod+"/"+name);
    }
}