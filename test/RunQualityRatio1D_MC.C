#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunQualityRatio1D_MC() {
    // !! Parameters to be changed !! //
    TString YEAR = "2016POST";
    TString CHANNEL = "MUMU"; // "MUMU" or "EMU"

    for (int i = 6 ; i < 14; ++i) { 
        plot(i, YEAR,CHANNEL); 
    }
}