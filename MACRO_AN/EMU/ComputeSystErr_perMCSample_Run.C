#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunComputeSystErr_perMCSample() {
    // Load the macro
    // gROOT->LoadMacro("plot4_ABCD_1Vtx_DATAMC.C");
    TCanvas *c1 ;
    // Loop to call the plot function with different integer arguments
    // !! Parameters to be changed !! //
    TString Name = "DATAMC_EMU_";
    TString Prod = "MC_EMU_03_02_2025";
    //RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1
    //RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9
    // RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17
    // RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11
    TString YEAR = "2018";
    TString DMODE = "EM";
    
    TString type = ".pdf";
    TString Sample [12]= {
        "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8",
        "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
        "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
        "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8",
        "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
        "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
        "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8",
        "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8",
        "TTWW_TuneCP5_13TeV-madgraph-pythia8",
        "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
        "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8",
        "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8"        
        };
    

    // !! -------------------------!! //

    int mixing = 0;

    for (unsigned int j = 2 ; j < 3 ; j++)
    {
        for (int i = 0; i < 29; ++i) { // 29
            std::cout << "i: " << i << std::endl;
                plot(i, Prod, Name, YEAR,DMODE , "test", Sample[j]); // Call the plot function from the macro with argument i
            }
    }

}