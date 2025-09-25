#include "TChain.h"
#include "MiniNtuple.h"
// C++ includes
#include <iostream>
#include <fstream>
#include "TROOT.h"

int main(int argc, char **argv)
{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniNtuple")==0) return 0;
 
 TChain c("FlyingTop/ttree");

TString Prod[5] = {
  "MC_EMU_2023_D_23_04_2025_DY","MC_EMU_2023_D_23_04_2025_DY",
  "MC_EMU_2024_23_06_2025_DY",
  "MC_MUMU_2023_D_23_04_2025_DY","MC_MUMU_2023_D_23_04_2025_DY"
};

// MC_EMU_2022_23_04_2025
 // MC_EMU_2022_EFG_23_04_2025
 // MC_EMU_2023_C_23_04_2025
 //- MC_EMU_2023_D_23_04_2025
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////
                 
     // //--------------------Background Mumu MiniNtuples ------------//
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


      for (int i = 0 ; i< 5 ; i++) 
         {
                  // /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024
            // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";

            TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod[i]+"/"+BKGSet[i]+".root";
            c.Reset();
            c.Add(Path);
            MiniNtuple* t = new MiniNtuple(&c);
            t->Loop(BKGSet[i],Prod[i],Signal);
         }
  
////////////////////////////////////////////////////////////////////////////////
   return 0;
}
