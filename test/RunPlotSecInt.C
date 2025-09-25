#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>
#include "plotSecInt.C"

void RunPlotSecInt() {
    // Load the macro
    // gROOT->LoadMacro("plotSecInt.C");
    // Loop to call the plot function with different integer arguments
    TString YEAR[3] = {"2022","2023","2024"};
    TString PLAN[4] = {"xy","rz","xy_Inner","rz_Inner"};
    TString SELECTION[2] = {"Selec","TrackerMatched"};
    TString PU[7] = {"","_PU25","_PU30","_PU35","_PU40","_PU45","_PU50"};
    TString Sample[3]={"DoubleMuon_2022","DoubleMuon_2023","DoubleMuon_2024"};


    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 4; ++j) {
            for (unsigned int k = 0; k < 2; ++k) {
                for (unsigned int l = 0; l < 7; ++l) {
                    plot(Sample[i],YEAR[i], PLAN[j], SELECTION[k], PU[l]);
                }
            }
        }
    }
// }


}

int main() {
    RunPlotSecInt();
    return 0;
}

