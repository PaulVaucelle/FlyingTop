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
// TString BKGSet[4]={"DoubleMuon_UL2018_MiniAODv2_GT36-v1","DoubleMuon_UL2017_MiniAODv2","DoubleMuon_UL2016PRE_MiniAODv2","DoubleMuon_UL2016POST_MiniAODv2"};
TString Prod = "Signal_2018";
TString Year = "2018";

//  for (int i = 0 ; i< 1 ; i++) 
//     {
//       // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+BKGSet[i]+".root";
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//       TreeSecIntReader* t = new TreeSecIntReader(&c,Prod, BKGSet[i]);
//       t->Loop(Prod, BKGSet[i]);
//       delete t;
      
//     }

      TString SignalSetctau[34]={"RPV_"+Year+"_smu200_neu180","RPV_"+Year+"_smu250_neu180","RPV_"+Year+"_smu250_neu200",
"RPV_"+Year+"_smu250_neu230","RPV_"+Year+"_smu300_neu180","RPV_"+Year+"_smu300_neu200","RPV_"+Year+"_smu300_neu280","RPV_"+Year+"_smu300_neu250",
"RPV_"+Year+"_smu350_neu180","RPV_"+Year+"_smu350_neu200","RPV_"+Year+"_smu350_neu250","RPV_"+Year+"_smu350_neu300","RPV_"+Year+"_smu350_neu330",
"RPV_"+Year+"_smu400_neu180","RPV_"+Year+"_smu400_neu200","RPV_"+Year+"_smu400_neu250","RPV_"+Year+"_smu400_neu300","RPV_"+Year+"_smu400_neu350",
"RPV_"+Year+"_smu400_neu380","RPV_"+Year+"_smu450_neu180","RPV_"+Year+"_smu450_neu200","RPV_"+Year+"_smu450_neu250","RPV_"+Year+"_smu450_neu300",
"RPV_"+Year+"_smu450_neu350","RPV_"+Year+"_smu450_neu400","RPV_"+Year+"_smu450_neu430","RPV_"+Year+"_smu500_neu180","RPV_"+Year+"_smu500_neu200",
"RPV_"+Year+"_smu500_neu250","RPV_"+Year+"_smu500_neu300","RPV_"+Year+"_smu500_neu350","RPV_"+Year+"_smu500_neu400","RPV_"+Year+"_smu500_neu450",
"RPV_"+Year+"_smu500_neu480"};     
 
 for (int i = 0 ; i< 34 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+SignalSetctau[i]+".root";
      c.Reset();
      c.Add(Path);
       TreeSecIntReader* t = new TreeSecIntReader(&c,Prod,SignalSetctau[i]);
       t->Loop(Prod,SignalSetctau[i]);
    }

//       TString SignalSet001[34]={"RPV_"+Year+"_smu200_neu180_ctau001","RPV_"+Year+"_smu250_neu180_ctau001","RPV_"+Year+"_smu250_neu200_ctau001",
// "RPV_"+Year+"_smu250_neu230_ctau001","RPV_"+Year+"_smu300_neu180_ctau001","RPV_"+Year+"_smu300_neu200_ctau001","RPV_"+Year+"_smu300_neu280_ctau001","RPV_"+Year+"_smu300_neu250_ctau001",
// "RPV_"+Year+"_smu350_neu180_ctau001","RPV_"+Year+"_smu350_neu200_ctau001","RPV_"+Year+"_smu350_neu250_ctau001","RPV_"+Year+"_smu350_neu300_ctau001","RPV_"+Year+"_smu350_neu330_ctau001",
// "RPV_"+Year+"_smu400_neu180_ctau001","RPV_"+Year+"_smu400_neu200_ctau001","RPV_"+Year+"_smu400_neu250_ctau001","RPV_"+Year+"_smu400_neu300_ctau001","RPV_"+Year+"_smu400_neu350_ctau001",
// "RPV_"+Year+"_smu400_neu380_ctau001","RPV_"+Year+"_smu450_neu180_ctau001","RPV_"+Year+"_smu450_neu200_ctau001","RPV_"+Year+"_smu450_neu250_ctau001","RPV_"+Year+"_smu450_neu300_ctau001",
// "RPV_"+Year+"_smu450_neu350_ctau001","RPV_"+Year+"_smu450_neu400_ctau001","RPV_"+Year+"_smu450_neu430_ctau001","RPV_"+Year+"_smu500_neu180_ctau001","RPV_"+Year+"_smu500_neu200_ctau001",
// "RPV_"+Year+"_smu500_neu250_ctau001","RPV_"+Year+"_smu500_neu300_ctau001","RPV_"+Year+"_smu500_neu350_ctau001","RPV_"+Year+"_smu500_neu400_ctau001","RPV_"+Year+"_smu500_neu450_ctau001",
// "RPV_"+Year+"_smu500_neu480_ctau001"};     
 
//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+SignalSet001[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TreeSecIntReader* t = new TreeSecIntReader(&c,Prod,SignalSet001[i]);
//        t->Loop(Prod,SignalSet001[i]);
//     }

// // // // ////////////////////////////////////////////////////////////////////////////////

//       TString SignalSet003[34]={"RPV_"+Year+"_smu200_neu180_ctau003","RPV_"+Year+"_smu250_neu180_ctau003","RPV_"+Year+"_smu250_neu200_ctau003",
// "RPV_"+Year+"_smu250_neu230_ctau003","RPV_"+Year+"_smu300_neu180_ctau003","RPV_"+Year+"_smu300_neu200_ctau003","RPV_"+Year+"_smu300_neu280_ctau003","RPV_"+Year+"_smu300_neu250_ctau003",
// "RPV_"+Year+"_smu350_neu180_ctau003","RPV_"+Year+"_smu350_neu200_ctau003","RPV_"+Year+"_smu350_neu250_ctau003","RPV_"+Year+"_smu350_neu300_ctau003","RPV_"+Year+"_smu350_neu330_ctau003",
// "RPV_"+Year+"_smu400_neu180_ctau003","RPV_"+Year+"_smu400_neu200_ctau003","RPV_"+Year+"_smu400_neu250_ctau003","RPV_"+Year+"_smu400_neu300_ctau003","RPV_"+Year+"_smu400_neu350_ctau003",
// "RPV_"+Year+"_smu400_neu380_ctau003","RPV_"+Year+"_smu450_neu180_ctau003","RPV_"+Year+"_smu450_neu200_ctau003","RPV_"+Year+"_smu450_neu250_ctau003","RPV_"+Year+"_smu450_neu300_ctau003",
// "RPV_"+Year+"_smu450_neu350_ctau003","RPV_"+Year+"_smu450_neu400_ctau003","RPV_"+Year+"_smu450_neu430_ctau003","RPV_"+Year+"_smu500_neu180_ctau003","RPV_"+Year+"_smu500_neu200_ctau003",
// "RPV_"+Year+"_smu500_neu250_ctau003","RPV_"+Year+"_smu500_neu300_ctau003","RPV_"+Year+"_smu500_neu350_ctau003","RPV_"+Year+"_smu500_neu400_ctau003","RPV_"+Year+"_smu500_neu450_ctau003",
// "RPV_"+Year+"_smu500_neu480_ctau003"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+SignalSet003[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TreeSecIntReader* t = new TreeSecIntReader(&c,Prod,SignalSet003[i]);
//        t->Loop(Prod,SignalSet003[i]);
//     }

//       TString SignalSet010[34]={"RPV_"+Year+"_smu200_neu180_ctau010","RPV_"+Year+"_smu250_neu180_ctau010","RPV_"+Year+"_smu250_neu200_ctau010",
// "RPV_"+Year+"_smu250_neu230_ctau010","RPV_"+Year+"_smu300_neu180_ctau010","RPV_"+Year+"_smu300_neu200_ctau010","RPV_"+Year+"_smu300_neu280_ctau010","RPV_"+Year+"_smu300_neu250_ctau010",
// "RPV_"+Year+"_smu350_neu180_ctau010","RPV_"+Year+"_smu350_neu200_ctau010","RPV_"+Year+"_smu350_neu250_ctau010","RPV_"+Year+"_smu350_neu300_ctau010","RPV_"+Year+"_smu350_neu330_ctau010",
// "RPV_"+Year+"_smu400_neu180_ctau010","RPV_"+Year+"_smu400_neu200_ctau010","RPV_"+Year+"_smu400_neu250_ctau010","RPV_"+Year+"_smu400_neu300_ctau010","RPV_"+Year+"_smu400_neu350_ctau010",
// "RPV_"+Year+"_smu400_neu380_ctau010","RPV_"+Year+"_smu450_neu180_ctau010","RPV_"+Year+"_smu450_neu200_ctau010","RPV_"+Year+"_smu450_neu250_ctau010","RPV_"+Year+"_smu450_neu300_ctau010",
// "RPV_"+Year+"_smu450_neu350_ctau010","RPV_"+Year+"_smu450_neu400_ctau010","RPV_"+Year+"_smu450_neu430_ctau010","RPV_"+Year+"_smu500_neu180_ctau010","RPV_"+Year+"_smu500_neu200_ctau010",
// "RPV_"+Year+"_smu500_neu250_ctau010","RPV_"+Year+"_smu500_neu300_ctau010","RPV_"+Year+"_smu500_neu350_ctau010","RPV_"+Year+"_smu500_neu400_ctau010","RPV_"+Year+"_smu500_neu450_ctau010",
// "RPV_"+Year+"_smu500_neu480_ctau010"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+SignalSet010[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TreeSecIntReader* t = new TreeSecIntReader(&c,Prod,SignalSet010[i]);
//        t->Loop(Prod,SignalSet010[i]);
//     }

//       TString SignalSet030[34]={"RPV_"+Year+"_smu200_neu180_ctau030","RPV_"+Year+"_smu250_neu180_ctau030","RPV_"+Year+"_smu250_neu200_ctau030",
// "RPV_"+Year+"_smu250_neu230_ctau030","RPV_"+Year+"_smu300_neu180_ctau030","RPV_"+Year+"_smu300_neu200_ctau030","RPV_"+Year+"_smu300_neu280_ctau030","RPV_"+Year+"_smu300_neu250_ctau030",
// "RPV_"+Year+"_smu350_neu180_ctau030","RPV_"+Year+"_smu350_neu200_ctau030","RPV_"+Year+"_smu350_neu250_ctau030","RPV_"+Year+"_smu350_neu300_ctau030","RPV_"+Year+"_smu350_neu330_ctau030",
// "RPV_"+Year+"_smu400_neu180_ctau030","RPV_"+Year+"_smu400_neu200_ctau030","RPV_"+Year+"_smu400_neu250_ctau030","RPV_"+Year+"_smu400_neu300_ctau030","RPV_"+Year+"_smu400_neu350_ctau030",
// "RPV_"+Year+"_smu400_neu380_ctau030","RPV_"+Year+"_smu450_neu180_ctau030","RPV_"+Year+"_smu450_neu200_ctau030","RPV_"+Year+"_smu450_neu250_ctau030","RPV_"+Year+"_smu450_neu300_ctau030",
// "RPV_"+Year+"_smu450_neu350_ctau030","RPV_"+Year+"_smu450_neu400_ctau030","RPV_"+Year+"_smu450_neu430_ctau030","RPV_"+Year+"_smu500_neu180_ctau030","RPV_"+Year+"_smu500_neu200_ctau030",
// "RPV_"+Year+"_smu500_neu250_ctau030","RPV_"+Year+"_smu500_neu300_ctau030","RPV_"+Year+"_smu500_neu350_ctau030","RPV_"+Year+"_smu500_neu400_ctau030","RPV_"+Year+"_smu500_neu450_ctau030",
// "RPV_"+Year+"_smu500_neu480_ctau030"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+SignalSet030[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TreeSecIntReader* t = new TreeSecIntReader(&c,Prod,SignalSet030[i]);
//        t->Loop(Prod,SignalSet030[i]);
//     }

//     TString SignalSet100[34]={"RPV_"+Year+"_smu200_neu180_ctau100","RPV_"+Year+"_smu250_neu180_ctau100","RPV_"+Year+"_smu250_neu200_ctau100",
// "RPV_"+Year+"_smu250_neu230_ctau100","RPV_"+Year+"_smu300_neu180_ctau100","RPV_"+Year+"_smu300_neu200_ctau100","RPV_"+Year+"_smu300_neu250_ctau100","RPV_"+Year+"_smu300_neu280_ctau100",
// "RPV_"+Year+"_smu350_neu180_ctau100","RPV_"+Year+"_smu350_neu200_ctau100","RPV_"+Year+"_smu350_neu250_ctau100","RPV_"+Year+"_smu350_neu300_ctau100","RPV_"+Year+"_smu350_neu330_ctau100",
// "RPV_"+Year+"_smu400_neu180_ctau100","RPV_"+Year+"_smu400_neu200_ctau100","RPV_"+Year+"_smu400_neu250_ctau100","RPV_"+Year+"_smu400_neu300_ctau100","RPV_"+Year+"_smu400_neu350_ctau100",
// "RPV_"+Year+"_smu400_neu380_ctau100","RPV_"+Year+"_smu450_neu180_ctau100","RPV_"+Year+"_smu450_neu200_ctau100","RPV_"+Year+"_smu450_neu250_ctau100","RPV_"+Year+"_smu450_neu300_ctau100",
// "RPV_"+Year+"_smu450_neu350_ctau100","RPV_"+Year+"_smu450_neu400_ctau100","RPV_"+Year+"_smu450_neu430_ctau100","RPV_"+Year+"_smu500_neu180_ctau100","RPV_"+Year+"_smu500_neu200_ctau100",
// "RPV_"+Year+"_smu500_neu250_ctau100","RPV_"+Year+"_smu500_neu300_ctau100","RPV_"+Year+"_smu500_neu350_ctau100","RPV_"+Year+"_smu500_neu400_ctau100","RPV_"+Year+"_smu500_neu450_ctau100",
// "RPV_"+Year+"_smu500_neu480_ctau100"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+SignalSet100[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TreeSecIntReader* t = new TreeSecIntReader(&c,Prod,SignalSet100[i]);
//        t->Loop(Prod,SignalSet100[i]);
//     }

//       TString SignalSet300[34]={"RPV_"+Year+"_smu200_neu180_ctau300","RPV_"+Year+"_smu250_neu180_ctau300","RPV_"+Year+"_smu250_neu200_ctau300",
// "RPV_"+Year+"_smu250_neu230_ctau300","RPV_"+Year+"_smu300_neu180_ctau300","RPV_"+Year+"_smu300_neu200_ctau300","RPV_"+Year+"_smu300_neu280_ctau300","RPV_"+Year+"_smu300_neu250_ctau300",
// "RPV_"+Year+"_smu350_neu180_ctau300","RPV_"+Year+"_smu350_neu200_ctau300","RPV_"+Year+"_smu350_neu250_ctau300","RPV_"+Year+"_smu350_neu300_ctau300","RPV_"+Year+"_smu350_neu330_ctau300",
// "RPV_"+Year+"_smu400_neu180_ctau300","RPV_"+Year+"_smu400_neu200_ctau300","RPV_"+Year+"_smu400_neu250_ctau300","RPV_"+Year+"_smu400_neu300_ctau300","RPV_"+Year+"_smu400_neu350_ctau300",
// "RPV_"+Year+"_smu400_neu380_ctau300","RPV_"+Year+"_smu450_neu180_ctau300","RPV_"+Year+"_smu450_neu200_ctau300","RPV_"+Year+"_smu450_neu250_ctau300","RPV_"+Year+"_smu450_neu300_ctau300",
// "RPV_"+Year+"_smu450_neu350_ctau300","RPV_"+Year+"_smu450_neu400_ctau300","RPV_"+Year+"_smu450_neu430_ctau300","RPV_"+Year+"_smu500_neu180_ctau300","RPV_"+Year+"_smu500_neu200_ctau300",
// "RPV_"+Year+"_smu500_neu250_ctau300","RPV_"+Year+"_smu500_neu300_ctau300","RPV_"+Year+"_smu500_neu350_ctau300","RPV_"+Year+"_smu500_neu400_ctau300","RPV_"+Year+"_smu500_neu450_ctau300",
// "RPV_"+Year+"_smu500_neu480_ctau300"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+SignalSet300[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TreeSecIntReader* t = new TreeSecIntReader(&c,Prod,SignalSet300[i]);
//        t->Loop(Prod,SignalSet300[i]);
//     }

//           TString SignalSet1000[34]={"RPV_"+Year+"_smu200_neu180_ctau1000","RPV_"+Year+"_smu250_neu180_ctau1000","RPV_"+Year+"_smu250_neu200_ctau1000",
// "RPV_"+Year+"_smu250_neu230_ctau1000","RPV_"+Year+"_smu300_neu180_ctau1000","RPV_"+Year+"_smu300_neu200_ctau1000","RPV_"+Year+"_smu300_neu280_ctau1000","RPV_"+Year+"_smu300_neu250_ctau1000",
// "RPV_"+Year+"_smu350_neu180_ctau1000","RPV_"+Year+"_smu350_neu200_ctau1000","RPV_"+Year+"_smu350_neu250_ctau1000","RPV_"+Year+"_smu350_neu300_ctau1000","RPV_"+Year+"_smu350_neu330_ctau1000",
// "RPV_"+Year+"_smu400_neu180_ctau1000","RPV_"+Year+"_smu400_neu200_ctau1000","RPV_"+Year+"_smu400_neu250_ctau1000","RPV_"+Year+"_smu400_neu300_ctau1000","RPV_"+Year+"_smu400_neu350_ctau1000",
// "RPV_"+Year+"_smu400_neu380_ctau1000","RPV_"+Year+"_smu450_neu180_ctau1000","RPV_"+Year+"_smu450_neu200_ctau1000","RPV_"+Year+"_smu450_neu250_ctau1000","RPV_"+Year+"_smu450_neu300_ctau1000",
// "RPV_"+Year+"_smu450_neu350_ctau1000","RPV_"+Year+"_smu450_neu400_ctau1000","RPV_"+Year+"_smu450_neu430_ctau1000","RPV_"+Year+"_smu500_neu180_ctau1000","RPV_"+Year+"_smu500_neu200_ctau1000",
// "RPV_"+Year+"_smu500_neu250_ctau1000","RPV_"+Year+"_smu500_neu300_ctau1000","RPV_"+Year+"_smu500_neu350_ctau1000","RPV_"+Year+"_smu500_neu400_ctau1000","RPV_"+Year+"_smu500_neu450_ctau1000",
// "RPV_"+Year+"_smu500_neu480_ctau1000"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+SignalSet1000[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TreeSecIntReader* t = new TreeSecIntReader(&c,Prod,SignalSet1000[i]);
//        t->Loop(Prod,SignalSet1000[i]);
//     }

}

