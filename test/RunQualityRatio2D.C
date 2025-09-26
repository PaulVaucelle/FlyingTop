#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunQualityRatio2D() {
    // !! Parameters to be changed !! //
    TString YEAR = "2018";
    TString CHANNEL = "EMU"; // "MUMU" or "EMU"
    // !! -------------------------!! ///RPV_2018_smu500_neu480.root


    for (int i = 0 ; i < 1; ++i) { 
        plot(i, YEAR,CHANNEL); 
    }



}