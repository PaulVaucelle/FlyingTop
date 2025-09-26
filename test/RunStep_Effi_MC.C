#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunStep_Effi_MC() {
    // !! Parameters to be changed !! //
    TString YEAR = "2018";
    TString CHANNEL = "EMU"; // "MUMU" or "EMU"
    TString Sample = "DY"; // TT or DY or nothing
    for (int i = 0 ; i < 1; ++i) { 
        plot(i, YEAR,CHANNEL,Sample); 
    }
}