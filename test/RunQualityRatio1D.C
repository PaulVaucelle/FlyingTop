#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunQualityRatio1D() {
    // !! Parameters to be changed !! //
    TString YEAR = "2018";
    TString CHANNEL = "MUMU"; // "MUMU" or "EMU"
    // !! -------------------------!! ///RPV_2018_smu500_neu480.root


    for (int i = 6 ; i < 14; ++i) { 
        plot(i, YEAR,CHANNEL); 
    }



}