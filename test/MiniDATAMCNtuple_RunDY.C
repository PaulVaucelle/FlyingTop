 #include "TChain.h"
#include "MiniDATAMCNtuple.h"
// C++ includes
#include <iostream>
#include <fstream>
#include "TROOT.h"

int main(int argc, char **argv){
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniDATAMCNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniDATAMCNtuple")==0) return 0;
 
 TChain c("FlyingTop/ttree");

// TString Prod = "MC_EMU_2022_EFG_23_04_2025";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
// MC_EMU_2022_23_04_2025
 // MC_EMU_2022_EFG_23_04_2025
 // MC_EMU_2023_C_23_04_2025
 //- MC_EMU_2023_D_23_04_2025

// TString Prod[4] = {"MC_MUMU_2022_23_04_2025","MC_MUMU_2022_EFG_23_04_2025","MC_MUMU_2023_C_23_04_2025","MC_MUMU_2023_D_23_04_2025"};
// TString Prod[1] = {"MC_MUMU_2024_23_06_2025_DY"};

bool Signal = false;
////////////////////////////////////////////////////////////////////////////////
     // //--------------------Background Mumu MiniDATAMCNtuples ------------//
  // TString BKGSet[2]={"DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8","DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"
  // };

  // TString BKGSet[1]={"DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8"
  // };

// TString BKGSet[10]={
//   "DYto2Mu_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8",
//   "DYto2Mu_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8",
//   "DYto2Mu_Bin-MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8",
//   "DYto2Mu_Bin-MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8",
//   "DYto2Mu_Bin-MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8",
//   "DYto2Tau_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8",
//   "DYto2Tau_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8",
//   "DYto2Tau_Bin-MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8",
//   "DYto2Tau_Bin-MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8",
//   "DYto2Tau_Bin-MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8"

//   // "DYto2Mu_MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8",
//   // "DYto2Mu_MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8",
//   // "DYto2Mu_MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8",
//   // "DYto2Mu_MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8",
//   // "DYto2Mu_MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8",
//   // "DYto2Tau_MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8",
//   // "DYto2Tau_MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8",
//   // "DYto2Tau_MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8",
//   // "DYto2Tau_MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8",
//   // "DYto2Tau_MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8"
//   };

TString Prod[5] = {
  "MC_EMU_2023_D_23_04_2025_DY","MC_EMU_2023_D_23_04_2025_DY",
  "MC_EMU_2024_23_06_2025_DY",
  "MC_MUMU_2023_D_23_04_2025_DY","MC_MUMU_2023_D_23_04_2025_DY"
};

  TString BKGSet[5]={
    "DYto2L-4Jets_MLL-10to50_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2Tau-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2L-4Jets_MLL-10to50_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    "DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8"
    // "DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8","DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8"
  // "DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8"
  //     "DYto2Mu_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8",
  // "DYto2Mu_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8",
  // "DYto2Mu_Bin-MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8",
  // "DYto2Mu_Bin-MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8",
  // "DYto2Mu_Bin-MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8",
  // "DYto2Tau_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8",
  // "DYto2Tau_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8",
  // "DYto2Tau_Bin-MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8",
  // "DYto2Tau_Bin-MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8",
  // "DYto2Tau_Bin-MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8"
  };




    for (int i = 1 ; i< 5 ; i++) 
        {
          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod[i]+"/"+BKGSet[i]+".root";
          c.Reset();
          c.Add(Path);
          MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
          t->Loop(BKGSet[i],Prod[i],Signal);
          delete t;
        }


    std::cout<<"----------------- 2016Pre Ended -------------------"<<std::endl;
////////////////////////////////////////////////////////////////////////////////
 return 0;
}
