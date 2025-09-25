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
    std::vector<TString> SYSTNAME = {"LumiUp","LumiDown","L1Up","L1Down","TriggerUp","TriggerDown",
    "MuonIDUp","MuonIDDown","MuonISOUp","MuonISODown","EleIDUp","EleIDDown","EleISOUp","EleISODown",
    "PUUp","PUDown","TopPtUp","TopPtDown",
    "ScaleUp","ScaleDown","JECUp","JECDown","JERUp","JERDown","RoccorUp","RoccorDown","PDFUp","PDFDown"};
    // 

    for (int i = 0; i < 8; ++i) //29
        { 
            std::cout << "i: " << i << std::endl;
            for (unsigned int j = 0 ; j < SYSTNAME.size()  ; j++) //SYSTNAME.size()
                {
                    plot(i, YEAR, SYSTNAME[j],"test"); // Call the plot function from the macro with argument i  
                }
            
        }
}