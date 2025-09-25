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

TString Prod = "DATA_EMU_2024_03_06_2025";
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////

// //--------------------DATA 2018 MiniNtuples ------------//
  // TString BKGSet[1]={ "MuonEG_Run2024"
  // };
  TString BKGSet[8]={ 
    "MuonEG_Run2024C-2024CDEReprocessing-v1",
    "MuonEG_Run2024D-2024CDEReprocessing-v1",
    "MuonEG_Run2024E-2024CDEReprocessing-v1"
    // "MuonEG_Run2024F-PromptReco-v1",
    // "MuonEG_Run2024G-PromptReco-v1",
    // "MuonEG_Run2024H-PromptReco-v1",
    // "MuonEG_Run2024I-PromptReco-v1",
    // "MuonEG_Run2024I-PromptReco-v2"

  };
 for (int i = 0 ; i< 3 ; i++) 
    {
      
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(BKGSet[i],Prod,Signal);
    }
 return 0;
}
