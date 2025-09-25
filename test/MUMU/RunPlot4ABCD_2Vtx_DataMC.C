#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot4ABCD_2Vtx_DataMC() {
    // Load the macro
    // gROOT->LoadMacro("plot4_ABCD_2Vtx_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    // !! Parameters to be changed !! //
    TString Name = "ABCD_";
    TString Prod = "DATA_MUMU_2018_19_08_2024";//pour Emu DATAMC2018_EMU_10_06_2024 //Mumu PROD_CSI_10_06_2024 //DATA_MUMU_2018_20_06_2024
    TString YEAR = "2022B";
    TString DMODE = "DM";
    bool DATA = false ;
    bool MC = true;
    bool SIGNAL = false;
    TString Plots[17] = {
    "STW",
    "HTL_VtxBDT",
    "HTL_VtxBDT_6Bins",
    "HTL_VtxBDT_7Bins",
    "HTL_VtxBDT_8Bins",
    "HTL_STW_8Bins",
    "HTL_STW_7Bins",
    "HTL_STW_6Bins",
    "HTL_VtxBDT_Sum",
    "HTL_VtxBDT_Ave",
    "HTL_STW_Sum",//10
    "HTL_STW_Ave",
    "STW_Closure",
    "EventBDT_Closure",
    "AveVtxBDT_Closure", //14
    "HTL_EventBDT"//15
    ,"HTL_EventBDT_Focus"
    };

    TString type = ".pdf";
    bool blind = true;
    // !! -------------------------!! //

    int mixing = 0;
    for (int i = 13; i < 17; ++i) { //87 
    std::cout<<"Plotting: " << Plots[i] << "with i = "<<i<<std::endl; 
        Plots[i]+= "_Corr"; 
        c1 = plot(i, Prod, Name, YEAR,DMODE, blind, Plots[i],DATA,MC,SIGNAL ); // Call the plot function from the macro with argument i
        TString name = Name+Plots[i]+DMODE;
                if (DATA && !MC) {name += "_DATA";}
        if (MC && !DATA) {name += "_MC";}
        if (MC && DATA) {name += "_DATAMC";}
        if (SIGNAL) {name += "_SIGNAL";}
        name+= "_"+YEAR;
        name += type;
        c1->SaveAs("./"+name);
    }
}