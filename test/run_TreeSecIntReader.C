#include "TChain.h"
#include "TreeSecIntReader.h"
// C++ includes
#include <iostream>
#include <fstream>
#include "TROOT.h"

int main(int argc, char **argv)
{
 gROOT->ProcessLine(".L TreeSecIntReader.C+");
 std::cout<<" Compiling TreeSecInt.C"<<std::endl;
 if (gROOT->GetClass("TreeSecIntReader")==0) return 0;
TChain c("ttree");


      // TString SignalSetq[3]={"RPV_2018_smu200_neu180_ctau001","RPV_2018_smu400_neu300_ctau010","RPV_2018_smu500_neu450_ctau100"};
    TString BKGSet[3]={"DoubleMuon_2022","DoubleMuon_2023","DoubleMuon_2024"};
    TString Prod = "SecInt";
 //issue  with RPV_2018_smu500_neu400_ctau001 500 450 001
 for (int i = 2 ; i< 3 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+BKGSet[i]+".root";
     std::cout<<" Running on "<<Path<<std::endl;
      c.Reset();
      c.Add(Path);
      TreeSecIntReader* t = new TreeSecIntReader(&c,Prod, BKGSet[i]);
      t->Loop(Prod, BKGSet[i]);
      delete t;
      
    }
 return 0;
}
