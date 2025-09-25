#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunQualityRatio1D_MC() {
    // !! Parameters to be changed !! //
    TString YEAR = "2022B";
    TString CHANNEL = "MUMU"; // "MUMU" or "EMU"

    for (int i = 6 ; i < 12; ++i) { 
        plot(i, YEAR,CHANNEL); 
    }
}