#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TString.h>
#include <stdio.h>


void RunPlot4ABCD_2Vtx_DataMC() {
    // Load the macro
    // gROOT->LoadMacro("plot4_ABCD_2Vtx_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    TString Year = "2018";
    TString Name = "Signal_"+Year+"_ABCD_";
    TString Prod = "Signal_"+Year+"_L1";//pour Emu DATAMC2018_EMU_10_06_2024 //Mumu PROD_CSI_10_06_2024 // PROD_ANNIVERSAIRE_2024
    TString Dmode = "DM";
    int MSmuon = 200;
    int MNeu = 180;
    TString ctau = "100";
    TString syst = "NOM";
    TString Plots[5] = {"diffBin_BestBDT"  ,"diffBin_SumSTW", "diffBin_AveSTW", "diffBin_STW","diffBin_EvtBDT"};
    TString type = ".pdf";
    TString NEUMASSES[13] = {"180","200","230","250","280","300","330","350","380","400","430","450","480"};
    TString SMUMASSES[7] = {"200","250","300","350","400","450","500"};
    int INTNEUMASSES[13] = {180,200,230,250,280,300,330,350,380,400,430,450,480};
    int INTSMUMASSES[7] = {200,250,300,350,400,450,500};
    for (unsigned int k = 0 ; k < 7 ; k++)//7
        {
            MSmuon = INTSMUMASSES[k];
            for (unsigned int j = 0 ; j < 13 ; j++)//13
                {
                    MNeu = INTNEUMASSES[j];
                    if (MNeu >= MSmuon) continue;
                    if (( ((MSmuon - MNeu) == 70)  || ((MSmuon - MNeu) == 120) || ((MSmuon - MNeu) == 170) || ((MSmuon - MNeu) == 220)  || ((MSmuon - MNeu) == 270)) && MNeu != 180) continue;
                    // if (MNeu == "230" && MSmuon == "250") continue;
                    // if (MNeu == "200" && MSmuon == "300") continue;//300_200_300
                    for (int i = 4; i < 5  ; ++i) // 4 to 28
                        { //89
                            c1 = plot(i, Prod, Name, Dmode,SMUMASSES[k], NEUMASSES[j], ctau,syst,Year); // Call the plot function from the macro with argument i
                            TString name = Name+SMUMASSES[k]+"_"+NEUMASSES[j]+"_"+Plots[i]+type;
                            c1->SaveAs("./"+name);
                        }
                }

        }

}