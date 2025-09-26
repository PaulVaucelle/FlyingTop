#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void ToysGenRun() {
    // !! Parameters to be changed !! //
    TString YEAR[4] = {"2018", "2017", "2016POST", "2016PRE"}; // Years to run the toys generation
    TString CHANNEL = "MUMU"; // "MUMU" or "EMU"
    TString VAR[3]= {"HTL_VtxBDT_Ave_Corr", "HTL_EventBDT_Corr", "HTL_STW_6Bins_Corr"};
    // !!Run

    for (int i = 1; i < 2; ++i) { 
        for (int j = 0; j < 1; ++j) { 
            ToysGen(VAR[i], CHANNEL, YEAR[j]); 
        }
    }


}

  // HTL_VtxBDT_Ave_Corr 
  // HTL_EventBDT_Corr
  // HTL_STW_6Bins_Corr 