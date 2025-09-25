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

TString Prod = "DATA_MUMU_2017_30_03_2025";
bool Signal = false;
// //--------------------DATA 2018 MiniDATAMCNtuples ------------//
  TString BKGProd[2]={ "DATA_EMU_2023_D_23_04_2025","DATA_MUMU_2023_D_23_04_2025" };// DM
  TString BKGSet[2]={"MuonEG_Run2023D-22Sep2023","Muon_Run2023D-22Sep2023" };// DM
  
for (int i = 0 ; i< 2 ; i++) 
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
