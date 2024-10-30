
#include "TChain.h"
#include "MiniSecIntNtuple.h"
// C++ includes
#include <iostream>
#include <fstream>
#include "TROOT.h"

int main(int argc, char **argv)
{
 gROOT->Reset() ; 

 // Compile user's analysis class //
  //  gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniSecIntNtuple.C+g") ;
  std::cout<<" Compiling MiniSecIntNtuple.C"<<std::endl;
 if (gROOT->GetClass("MiniSecIntNtuple")==0) return 0;
 
 TChain c("FlyingTop/ttree");

TString Prod = "SecInt";
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////
// /opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY
// /Data_2016pre_emu/240603_bdt94/Data*.root
// /Data_2016post_emu/240603_bdt94/Data*.root
// /Data_2017_emu/240603_bdt94/Data*.root
// //--------------------DATA 2018 MiniSecIntNtuples ------------//
  TString BKGSet[1]={"DoubleMuon_2023"};


// "emu_2018A",
// "emu_2018B",
 for (int i = 0 ; i< 3 ; i++) 
    {
      // TString Path = "/opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY/Data_2023_emu/240819/"+BKGSet[i]+".root";
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/DATA_MUMU_2023_25_09_2024/"+BKGSet[i]+".root";
      std::cout<<" Running on "<<Path<<std::endl;
      c.Reset();
      c.Add(Path);
       MiniSecIntNtuple* t = new MiniSecIntNtuple(&c);
       t->Loop(BKGSet[i],Prod,Signal);
    }
return 0;
}
