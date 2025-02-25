#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void GatherSystperSyst() {

    // !! Parameters to be changed !! //
    TString Name = "DATAMC_EMU_";
    TString YEAR = "2018";
    // !! -------------------------!! //
    std::vector<TString> SYSTNAME = {"LumiUp","LumiDown","L1Up","L1Down","TriggerUp","TriggerDown","LepIDUp",
    "LepIDDown","LepISOUp","LepISODown","PUUp","PUDown","JECUp","JECDown","JERUp","JERDown"};

    for (int i = 0; i < 29; ++i) //29
        { 
            std::cout << "i: " << i << std::endl;
            for (unsigned int j = 0 ; j < SYSTNAME.size()  ; j++) 
                {
                    plot(i, YEAR, SYSTNAME[j],"test"); // Call the plot function from the macro with argument i  
                }
            
        }
}