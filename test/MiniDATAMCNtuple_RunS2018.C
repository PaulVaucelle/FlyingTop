{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniDATAMCNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniDATAMCNtuple")==0) return;
 
 TChain c("FlyingTop/ttree");

TString Year = "2018";
TString Prod = "Signal_2018";//Signal_2018

bool Signal = true;
////////////////////////////////////////////////////////////////////////////////
   


//       TString SignalSet001[34]={"RPV_"+Year+"_smu200_neu180_ctau001","RPV_"+Year+"_smu250_neu180_ctau001","RPV_"+Year+"_smu250_neu200_ctau001",
// "RPV_"+Year+"_smu250_neu230_ctau001","RPV_"+Year+"_smu300_neu180_ctau001","RPV_"+Year+"_smu300_neu200_ctau001","RPV_"+Year+"_smu300_neu280_ctau001",
// "RPV_"+Year+"_smu350_neu180_ctau001","RPV_"+Year+"_smu350_neu200_ctau001","RPV_"+Year+"_smu350_neu250_ctau001","RPV_"+Year+"_smu350_neu300_ctau001","RPV_"+Year+"_smu350_neu330_ctau001",
// "RPV_"+Year+"_smu400_neu180_ctau001","RPV_"+Year+"_smu400_neu200_ctau001","RPV_"+Year+"_smu400_neu250_ctau001","RPV_"+Year+"_smu400_neu300_ctau001","RPV_"+Year+"_smu400_neu350_ctau001",
// "RPV_"+Year+"_smu400_neu380_ctau001","RPV_"+Year+"_smu450_neu180_ctau001","RPV_"+Year+"_smu450_neu200_ctau001","RPV_"+Year+"_smu450_neu250_ctau001","RPV_"+Year+"_smu450_neu300_ctau001",
// "RPV_"+Year+"_smu450_neu350_ctau001","RPV_"+Year+"_smu450_neu400_ctau001","RPV_"+Year+"_smu450_neu430_ctau001","RPV_"+Year+"_smu500_neu180_ctau001","RPV_"+Year+"_smu500_neu200_ctau001",
// "RPV_"+Year+"_smu500_neu250_ctau001","RPV_"+Year+"_smu500_neu300_ctau001","RPV_"+Year+"_smu500_neu350_ctau001","RPV_"+Year+"_smu500_neu400_ctau001","RPV_"+Year+"_smu500_neu450_ctau001",
// "RPV_"+Year+"_smu500_neu480_ctau001"};     
 
//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet001[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
//        t->Loop(SignalSet001[i],Prod,Signal);
//     }

// // // // ////////////////////////////////////////////////////////////////////////////////

//       TString SignalSet003[34]={"RPV_"+Year+"_smu200_neu180_ctau003","RPV_"+Year+"_smu250_neu180_ctau003","RPV_"+Year+"_smu250_neu200_ctau003",
// "RPV_"+Year+"_smu250_neu230_ctau003","RPV_"+Year+"_smu300_neu180_ctau003","RPV_"+Year+"_smu300_neu200_ctau003","RPV_"+Year+"_smu300_neu280_ctau003",
// "RPV_"+Year+"_smu350_neu180_ctau003","RPV_"+Year+"_smu350_neu200_ctau003","RPV_"+Year+"_smu350_neu250_ctau003"/*,"RPV_"+Year+"_smu350_neu300_ctau003",*/,"RPV_"+Year+"_smu350_neu330_ctau003",
// "RPV_"+Year+"_smu400_neu180_ctau003","RPV_"+Year+"_smu400_neu200_ctau003","RPV_"+Year+"_smu400_neu250_ctau003","RPV_"+Year+"_smu400_neu300_ctau003","RPV_"+Year+"_smu400_neu350_ctau003",
// "RPV_"+Year+"_smu400_neu380_ctau003","RPV_"+Year+"_smu450_neu180_ctau003","RPV_"+Year+"_smu450_neu200_ctau003","RPV_"+Year+"_smu450_neu250_ctau003","RPV_"+Year+"_smu450_neu300_ctau003",
// "RPV_"+Year+"_smu450_neu350_ctau003","RPV_"+Year+"_smu450_neu400_ctau003","RPV_"+Year+"_smu450_neu430_ctau003","RPV_"+Year+"_smu500_neu180_ctau003","RPV_"+Year+"_smu500_neu200_ctau003",
// "RPV_"+Year+"_smu500_neu250_ctau003","RPV_"+Year+"_smu500_neu300_ctau003","RPV_"+Year+"_smu500_neu350_ctau003","RPV_"+Year+"_smu500_neu400_ctau003","RPV_"+Year+"_smu500_neu450_ctau003",
// "RPV_"+Year+"_smu500_neu480_ctau003"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet003[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
//        t->Loop(SignalSet003[i],Prod,Signal);
//     }

//       TString SignalSet010[34]={"RPV_"+Year+"_smu200_neu180_ctau010","RPV_"+Year+"_smu250_neu180_ctau010","RPV_"+Year+"_smu250_neu200_ctau010",
// "RPV_"+Year+"_smu250_neu230_ctau010","RPV_"+Year+"_smu300_neu180_ctau010","RPV_"+Year+"_smu300_neu200_ctau010","RPV_"+Year+"_smu300_neu280_ctau010",
// "RPV_"+Year+"_smu350_neu180_ctau010","RPV_"+Year+"_smu350_neu200_ctau010","RPV_"+Year+"_smu350_neu250_ctau010","RPV_"+Year+"_smu350_neu300_ctau010","RPV_"+Year+"_smu350_neu330_ctau010",
// "RPV_"+Year+"_smu400_neu180_ctau010","RPV_"+Year+"_smu400_neu200_ctau010","RPV_"+Year+"_smu400_neu250_ctau010","RPV_"+Year+"_smu400_neu300_ctau010","RPV_"+Year+"_smu400_neu350_ctau010",
// "RPV_"+Year+"_smu400_neu380_ctau010","RPV_"+Year+"_smu450_neu180_ctau010","RPV_"+Year+"_smu450_neu200_ctau010","RPV_"+Year+"_smu450_neu250_ctau010","RPV_"+Year+"_smu450_neu300_ctau010",
// "RPV_"+Year+"_smu450_neu350_ctau010","RPV_"+Year+"_smu450_neu400_ctau010","RPV_"+Year+"_smu450_neu430_ctau010","RPV_"+Year+"_smu500_neu180_ctau010","RPV_"+Year+"_smu500_neu200_ctau010",
// "RPV_"+Year+"_smu500_neu250_ctau010","RPV_"+Year+"_smu500_neu300_ctau010","RPV_"+Year+"_smu500_neu350_ctau010","RPV_"+Year+"_smu500_neu400_ctau010","RPV_"+Year+"_smu500_neu450_ctau010",
// "RPV_"+Year+"_smu500_neu480_ctau010"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet010[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
//        t->Loop(SignalSet010[i],Prod,Signal);
//     }

//       TString SignalSet030[34]={"RPV_"+Year+"_smu200_neu180_ctau030","RPV_"+Year+"_smu250_neu180_ctau030","RPV_"+Year+"_smu250_neu200_ctau030",
// "RPV_"+Year+"_smu250_neu230_ctau030","RPV_"+Year+"_smu300_neu180_ctau030","RPV_"+Year+"_smu300_neu200_ctau030","RPV_"+Year+"_smu300_neu280_ctau030",
// "RPV_"+Year+"_smu350_neu180_ctau030","RPV_"+Year+"_smu350_neu200_ctau030","RPV_"+Year+"_smu350_neu250_ctau030","RPV_"+Year+"_smu350_neu300_ctau030","RPV_"+Year+"_smu350_neu330_ctau030",
// "RPV_"+Year+"_smu400_neu180_ctau030","RPV_"+Year+"_smu400_neu200_ctau030","RPV_"+Year+"_smu400_neu250_ctau030","RPV_"+Year+"_smu400_neu300_ctau030","RPV_"+Year+"_smu400_neu350_ctau030",
// "RPV_"+Year+"_smu400_neu380_ctau030","RPV_"+Year+"_smu450_neu180_ctau030","RPV_"+Year+"_smu450_neu200_ctau030","RPV_"+Year+"_smu450_neu250_ctau030","RPV_"+Year+"_smu450_neu300_ctau030",
// "RPV_"+Year+"_smu450_neu350_ctau030","RPV_"+Year+"_smu450_neu400_ctau030","RPV_"+Year+"_smu450_neu430_ctau030","RPV_"+Year+"_smu500_neu180_ctau030","RPV_"+Year+"_smu500_neu200_ctau030",
// "RPV_"+Year+"_smu500_neu250_ctau030","RPV_"+Year+"_smu500_neu300_ctau030","RPV_"+Year+"_smu500_neu350_ctau030","RPV_"+Year+"_smu500_neu400_ctau030","RPV_"+Year+"_smu500_neu450_ctau030",
// "RPV_"+Year+"_smu500_neu480_ctau030"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet030[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
//        t->Loop(SignalSet030[i],Prod,Signal);
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
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet100[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
//        t->Loop(SignalSet100[i],Prod,Signal);
//     }

//       TString SignalSet300[34]={"RPV_"+Year+"_smu200_neu180_ctau300","RPV_"+Year+"_smu250_neu180_ctau300","RPV_"+Year+"_smu250_neu200_ctau300",
// "RPV_"+Year+"_smu250_neu230_ctau300","RPV_"+Year+"_smu300_neu180_ctau300","RPV_"+Year+"_smu300_neu200_ctau300","RPV_"+Year+"_smu300_neu280_ctau300",
// "RPV_"+Year+"_smu350_neu180_ctau300","RPV_"+Year+"_smu350_neu200_ctau300","RPV_"+Year+"_smu350_neu250_ctau300","RPV_"+Year+"_smu350_neu300_ctau300","RPV_"+Year+"_smu350_neu330_ctau300",
// "RPV_"+Year+"_smu400_neu180_ctau300","RPV_"+Year+"_smu400_neu200_ctau300","RPV_"+Year+"_smu400_neu250_ctau300","RPV_"+Year+"_smu400_neu300_ctau300","RPV_"+Year+"_smu400_neu350_ctau300",
// "RPV_"+Year+"_smu400_neu380_ctau300","RPV_"+Year+"_smu450_neu180_ctau300","RPV_"+Year+"_smu450_neu200_ctau300","RPV_"+Year+"_smu450_neu250_ctau300","RPV_"+Year+"_smu450_neu300_ctau300",
// "RPV_"+Year+"_smu450_neu350_ctau300","RPV_"+Year+"_smu450_neu400_ctau300","RPV_"+Year+"_smu450_neu430_ctau300","RPV_"+Year+"_smu500_neu180_ctau300","RPV_"+Year+"_smu500_neu200_ctau300",
// "RPV_"+Year+"_smu500_neu250_ctau300","RPV_"+Year+"_smu500_neu300_ctau300","RPV_"+Year+"_smu500_neu350_ctau300","RPV_"+Year+"_smu500_neu400_ctau300","RPV_"+Year+"_smu500_neu450_ctau300",
// "RPV_"+Year+"_smu500_neu480_ctau300"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet300[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
//        t->Loop(SignalSet300[i],Prod,Signal);
//     }

//           TString SignalSet1000[34]={"RPV_"+Year+"_smu200_neu180_ctau1000","RPV_"+Year+"_smu250_neu180_ctau1000","RPV_"+Year+"_smu250_neu200_ctau1000",
// "RPV_"+Year+"_smu250_neu230_ctau1000","RPV_"+Year+"_smu300_neu180_ctau1000","RPV_"+Year+"_smu300_neu200_ctau1000","RPV_"+Year+"_smu300_neu280_ctau1000",
// "RPV_"+Year+"_smu350_neu180_ctau1000","RPV_"+Year+"_smu350_neu200_ctau1000","RPV_"+Year+"_smu350_neu250_ctau1000","RPV_"+Year+"_smu350_neu300_ctau1000","RPV_"+Year+"_smu350_neu330_ctau1000",
// "RPV_"+Year+"_smu400_neu180_ctau1000","RPV_"+Year+"_smu400_neu200_ctau1000","RPV_"+Year+"_smu400_neu250_ctau1000","RPV_"+Year+"_smu400_neu300_ctau1000","RPV_"+Year+"_smu400_neu350_ctau1000",
// "RPV_"+Year+"_smu400_neu380_ctau1000","RPV_"+Year+"_smu450_neu180_ctau1000","RPV_"+Year+"_smu450_neu200_ctau1000","RPV_"+Year+"_smu450_neu250_ctau1000","RPV_"+Year+"_smu450_neu300_ctau1000",
// "RPV_"+Year+"_smu450_neu350_ctau1000","RPV_"+Year+"_smu450_neu400_ctau1000","RPV_"+Year+"_smu450_neu430_ctau1000","RPV_"+Year+"_smu500_neu180_ctau1000","RPV_"+Year+"_smu500_neu200_ctau1000",
// "RPV_"+Year+"_smu500_neu250_ctau1000","RPV_"+Year+"_smu500_neu300_ctau1000","RPV_"+Year+"_smu500_neu350_ctau1000","RPV_"+Year+"_smu500_neu400_ctau1000","RPV_"+Year+"_smu500_neu450_ctau1000",
// "RPV_"+Year+"_smu500_neu480_ctau1000"};

//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet1000[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
//        t->Loop(SignalSet1000[i],Prod,Signal);
//     }

// // //////////////////////////////////////////////////////////////////////////////
TString SignalSet[4]={"RPV_2018_smu300_neu180_ctau100","RPV_2018_smu400_neu300_ctau100","RPV_2018_smu200_neu180_ctau100","RPV_2018_smu500_neu350_ctau100"};


 for (int i = 0 ; i< 5 ; i++) 
    {
      std::cout<<"SignalSet[i] "<<SignalSet[i]<<std::endl;
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
       t->Loop(SignalSet[i],Prod,Signal);
    }

}
