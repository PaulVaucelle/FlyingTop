#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>


void GatherSystperSyst() {

    // !! Parameters to be changed !! //
    TString YEAR = "2018";
    TString Prod = "Signal_2018_L1";
    // !! -------------------------!! //
    std::vector<TString> SYSTNAME = {"LumiUp","LumiDown","L1Up","L1Down","TriggerUp","TriggerDown",
    "MuonIDUp","MuonIDDown","MuonISOUp","MuonISODown","PUUp","PUDown","TopPtUp","TopPtDown",
    "JECUp","JECDown","JERUp","JERDown","XSUp","XSDown"};//"RoccorUp","RoccorDown",
    // 

    TString MSMU = "500";
    TString MNEU = "300";
    TString CTAU = "001";
    // Définition des listes
    std::vector<int> list1 = {180, 200, 230, 250, 280, 300, 330, 350, 380, 400, 430, 450, 480}; // Masses des neutralinos
    std::vector<int> list2 = {200, 250, 300, 350, 400, 450, 500}; // Masses des smuons
    std::vector<TString> list3 = {"001", "003", "010", "030", "100", "300", "1000"};    // ctaus

    // for (int i = 0; i < 9; ++i) //9 : nombre de variables
    //     { 
    //         std::cout << "i: " << i << std::endl;
    //         for (unsigned int j = 0 ; j < SYSTNAME.size()  ; j++) //SYSTNAME.size()
    //             {
    //                 plot(i, YEAR, SYSTNAME[j],"test", MSMU, MNEU, CTAU); // Call the plot function from the macro with argument i  
    //             }
            
    //     }

    // Boucles imbriquées
    for (TString ctau : list3) 
        {
            for (int msmu : list2) 
                {
                    for (int mneu : list1) 
                        {
                            if (mneu < msmu && (mneu ==180 || msmu-mneu == 20 || (msmu-mneu)% 50 ==0))
                                {

                                    MSMU = std::to_string(msmu) ;
                                    MNEU = std::to_string(mneu) ;
                                    CTAU = ctau;

                                    for (int i = 0; i < 9; ++i) //9 : nombre de variables
                                        { 
                                            std::cout << "i: " << i << std::endl;
                                            for (unsigned int j = 0 ; j < SYSTNAME.size()  ; j++) //SYSTNAME.size()
                                                {
                                                    plot(i, YEAR, SYSTNAME[j],"test", MSMU, MNEU, CTAU); // Call the plot function from the macro with argument i  
                                                }
                                        }
                                }
                        }
                }
        }

}