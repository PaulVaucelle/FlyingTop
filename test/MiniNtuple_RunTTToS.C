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

TString Prod[4] = {"MC_EMU_2024_23_06_2025","MC_MUMU_2024_23_06_2025"};



bool Signal = false;
////////////////////////////////////////////////////////////////////////////////

// //--------------------Background Emu MiniNtuples ------------//
  TString BKGSet[1]={"TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8"
  };

for (int j = 0 ; j < 2 ; j++)
  { 
    for (int i = 0 ; i< 1 ; i++) 
        {
          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod[j]+"/"+BKGSet[i]+".root";
          c.Reset();
          c.Add(Path);
          MiniNtuple* t = new MiniNtuple(&c);
          t->Loop(BKGSet[i],Prod[j],Signal);
        }
  }
////////////////////////////////////////////////////////////////////////////////
return 0;
}
