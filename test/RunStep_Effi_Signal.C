#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunStep_Effi_MC() {
    // !! Parameters to be changed !! //
    TString YEAR = "2018";

    std::vector<TString> FILENAMES;
    FILENAMES.push_back("RPV_2018_smu200_neu180_ctau100");
    FILENAMES.push_back("RPV_2018_smu300_neu200_ctau100");
    FILENAMES.push_back("RPV_2018_smu400_neu200_ctau100");
    FILENAMES.push_back("RPV_2018_smu500_neu200_ctau100");

    for (unsigned int i = 0 ; i < FILENAMES.size(); ++i) {
        plot(0, YEAR,FILENAMES[i]); 
    }

}