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

// TString Prod = "DATA_MUMU_2022_25_09_2024";
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////

// //--------------------DATA 2018 MiniNtuples ------------//
TString Prod[8]= {"DATA_EMU_2022_FG_23_04_2025","DATA_EMU_2022_FG_23_04_2025",  
"DATA_MUMU_2022_CDE_23_04_2025","DATA_MUMU_2022_CDE_23_04_2025", "DATA_MUMU_2022_CDE_23_04_2025",

"DATA_MUMU_2022_FG_23_04_2025","DATA_MUMU_2022_FG_23_04_2025",

"DATA_MUMU_2023_C_23_04_2025"  
};

  TString BKGSet[8]={ "MuonEG_Run2022F-22Sep2023-v1","MuonEG_Run2022G-22Sep2023-v1",
  "Muon_Run2022C-22Sep2023","Muon_Run2022D-22Sep2023",   "Muon_Run2022E-22Sep2023",
  "Muon_Run2022F-22Sep2023","Muon_Run2022G-22Sep2023",

  "Muon_Run2023C-22Sep2023"
  };

 for (int i = 0 ; i< 8 ; i++) 
    {
      
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod[i]+"/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(BKGSet[i],Prod[i],Signal);
    }
 return 0;
}
