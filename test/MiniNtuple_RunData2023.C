
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

TString Prod = "DATA_MUMU_2023_25_09_2024";
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////

// //--------------------DATA 2018 MiniNtuples ------------//
  TString BKGSet[1]={ "DoubleMuon_2023"
  };
  // "SM_Run2018A-UL2018_MiniAODv2-v3","SM_Run2018B-UL2018_MiniAODv2-v2","SM_Run2018C-UL2018_MiniAODv2-v2","SM_Run2018D-UL2018_MiniAODv2-v3"
// "emu_2018A",
// "emu_2018B",
 for (int i = 0 ; i< 1 ; i++) 
    {
      
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(BKGSet[i],Prod,Signal);
    }
  return 0;
}
