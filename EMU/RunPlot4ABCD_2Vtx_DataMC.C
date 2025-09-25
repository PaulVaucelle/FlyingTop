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
    TString Name = "ABCD_EMU_";
    TString Prod = "DATA_EMU_2018_19_08_2024";
    
    TString YEAR = "2024";
    TString DMODE = "EM";
    bool DATA = true ;
    bool MC = false;
    bool SIGNAL = false;
    TString Plots[16] = {
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
    };

    TString type = ".pdf";  

    // !! -------------------------!! //

    int mixing = 0;
    for (int i =0; i < 16; ++i) { //87 
    std::cout<<"Plotting: " << Plots[i] << "with i = "<<i<<std::endl;
        Plots[i]+= "_Corr";
        c1 = plot(i, Prod, Name, YEAR,DMODE , Plots[i], DATA, MC,SIGNAL); // Call the plot function from the macro with argument i
        TString name = Name+Plots[i]+DMODE;
        if (DATA && !MC) {name += "_DATA";}
        if (MC && !DATA) {name += "_MC";}
        if (MC && DATA) {name += "_DATAMC";}
        name+= "_"+YEAR;
        name += type;
        c1->SaveAs("./"+name);
    }
}