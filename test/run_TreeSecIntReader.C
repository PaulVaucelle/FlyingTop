// #include "TChain.h"
// #include "TreeSecIntReader.h"
// // C++ includes
// #include <iostream>
// #include <fstream>
// #include "TROOT.h"

// int run_TreeSecIntReader()
{
 gROOT->Reset() ; 
 gROOT->ProcessLine(".L TreeSecIntReader.C+g");

 if (gROOT->GetClass("TreeSecIntReader")==0) return ;
  std::cout<<" Compiling TreeSecInt.C"<<std::endl;
TChain c("ttree");


      // TString SignalSetq[3]={"RPV_2018_smu200_neu180_ctau001","RPV_2018_smu400_neu300_ctau010","RPV_2018_smu500_neu450_ctau100"};
    TString BKGSet[4]={"DoubleMuon_UL2018_MiniAODv2_GT36-v1","DoubleMuon_UL2017_MiniAODv2","DoubleMuon_UL2016PRE_MiniAODv2","DoubleMuon_UL2016POST_MiniAODv2"};
        // TString BKGSet[1]={"DoubleMuon_UL2018_MiniAODv2_GT36-v1"};
    // TString BKGSet[3]={"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8","DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8","TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"};
        // TString BKGSet[1]={"MCBKG"};

    TString Prod = "SecInt";
 //issue  with RPV_2018_smu500_neu400_ctau001 500 450 001
 for (int i = 0 ; i< 1 ; i++) 
    {
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+BKGSet[i]+".root";
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeSecIntReader* t = new TreeSecIntReader(&c,Prod, BKGSet[i]);
      t->Loop(Prod, BKGSet[i]);
      delete t;
      
    }
//  return 0 ;
}

