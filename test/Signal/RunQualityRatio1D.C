#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunQualityRatio1D() {
    // !! Parameters to be changed !! //
    TString YEAR = "2022A"; // 2022A, 2022B, 2023A, 2023B, 2024
    TString CHANNEL = "MUMU"; // "MUMU" or "EMU"
    // !! -------------------------!! ///RPV_2018_smu500_neu480.root

    for (int i = 6 ; i < 12; ++i) { 
        plot(i, YEAR,CHANNEL); 
    }



}