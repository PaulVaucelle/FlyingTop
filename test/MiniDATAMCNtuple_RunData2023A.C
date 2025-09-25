#include "TChain.h"
#include "MiniDATAMCNtuple.h"
// C++ includes
#include <iostream>
#include <fstream>
#include "TROOT.h"

int main(int argc, char **argv)
{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniDATAMCNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniDATAMCNtuple")==0) return 0;
 
 TChain c("FlyingTop/ttree");

TString Prod = "DATA_EMU_2023_C_23_04_2025";
bool Signal = false;

// //--------------------DATA 2018 MiniDATAMCNtuples ------------//
  TString BKGProd[8]=
  { 
    "DATA_MUMU_2023_C_28_06_2025",
    "DATA_MUMU_2023_C_28_06_2025",
    "DATA_MUMU_2023_C_28_06_2025",
    "DATA_MUMU_2023_C_28_06_2025_v2",
    "DATA_MUMU_2023_C_28_06_2025_v2",
    "DATA_MUMU_2023_C_28_06_2025_v3",
    "DATA_MUMU_2023_C_28_06_2025_v3",
    "DATA_MUMU_2023_C_28_06_2025_v3",

    };// DM


  TString BKGSet[8]=
  {
    "Muon0_Run2023C-22Sep2023_v2",
    "Muon0_Run2023C-22Sep2023_v3",
    "Muon0_Run2023C-22Sep2023_v4",
    "Muon1_Run2023C-22Sep2023_v3",
    "Muon1_Run2023C-22Sep2023_v4",
    "Muon0_Run2023C-22Sep2023_v1",
    "Muon1_Run2023C-22Sep2023_v1",
    "Muon1_Run2023C-22Sep2023_v2"
   
    };// DM
  
for (int i = 0 ; i< 8 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+BKGProd[i]+"/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
       t->Loop(BKGSet[i],BKGProd[i],Signal);
      delete t;
    }
  return 0;
}
