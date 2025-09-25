
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

TString Prod[2] = {"MC_MUMU_2023_C_23_04_2025","MC_MUMU_2023_D_23_04_2025"};


// MC_EMU_2022_23_04_2025
 // MC_EMU_2022_EFG_23_04_2025
 // MC_EMU_2023_C_23_04_2025
 //- MC_EMU_2023_D_23_04_2025
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////


// //--------------------Background Mumu MiniNtuples ------------//
   TString BKGSet[8]={ 
  "TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
  "TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
  "WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
  "WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8",
  "ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8",
  "TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8",
  "TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8",
  "TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
  };
  
 for (int j = 1 ; j < 2 ; j++)
  {     
    for (int i = 0 ; i< 8 ; i++) 
        {
          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod[j]+"/"+BKGSet[i]+".root";
          c.Reset();
          c.Add(Path);
          MiniNtuple* t = new MiniNtuple(&c);
          t->Loop(BKGSet[i],Prod[j],Signal);
        }
        
      }
    return 0;
}
