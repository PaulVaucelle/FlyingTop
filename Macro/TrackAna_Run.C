{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L ../HistogramManager.C+g") ;
   gROOT->ProcessLine(".L TrackAna.C+g") ;
  
 if (gROOT->GetClass("TrackAna")==0) return;//LLTopTree
 
 TChain c("FlyingTop/ttree");
TString Prod = "PROD_CSI_10_06_2024";
// TString Sample = "RPV_2018_smu300_neu200_ctau100.root";//RPV_2018_smu300_neu200_ctau100.root //TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root
// TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+Sample;
// c.Add(Path);

//  TrackAna* t = new TrackAna(&c);
//  t->Loop(Prod, Sample);


TString BkgSet = "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8";
TString path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+BkgSet+".root";
       
      c.Add(path);
       TrackAna* t = new TrackAna(&c);
       t->Loop(Prod, BkgSet);

//       TString SignalSet001[32]={"RPV_2018_smu200_neu180_ctau001","RPV_2018_smu250_neu180_ctau001","RPV_2018_smu250_neu200_ctau001",
// "RPV_2018_smu250_neu230_ctau001","RPV_2018_smu300_neu180_ctau001","RPV_2018_smu300_neu200_ctau001","RPV_2018_smu300_neu280_ctau001",
// "RPV_2018_smu350_neu180_ctau001","RPV_2018_smu350_neu200_ctau001","RPV_2018_smu350_neu250_ctau001","RPV_2018_smu350_neu300_ctau001","RPV_2018_smu350_neu330_ctau001",
// "RPV_2018_smu400_neu180_ctau001","RPV_2018_smu400_neu200_ctau001","RPV_2018_smu400_neu250_ctau001","RPV_2018_smu400_neu300_ctau001","RPV_2018_smu400_neu350_ctau001",
// "RPV_2018_smu400_neu380_ctau001","RPV_2018_smu450_neu180_ctau001","RPV_2018_smu450_neu200_ctau001","RPV_2018_smu450_neu250_ctau001","RPV_2018_smu450_neu300_ctau001",
// "RPV_2018_smu450_neu350_ctau001","RPV_2018_smu450_neu400_ctau001","RPV_2018_smu450_neu430_ctau001","RPV_2018_smu500_neu180_ctau001","RPV_2018_smu500_neu200_ctau001",
// "RPV_2018_smu500_neu250_ctau001","RPV_2018_smu500_neu300_ctau001","RPV_2018_smu500_neu350_ctau001",/*"RPV_2018_smu500_neu400_ctau001","RPV_2018_smu500_neu450_ctau001",*/
// "RPV_2018_smu500_neu480_ctau001"};     
 
//  for (int i = 0 ; i< 31 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet001[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TrackAna* t = new TrackAna(&c);
//        t->Loop(Prod, SignalSet001[i]);
//     }

// ////////////////////////////////////////////////////////////////////////////////

//       TString SignalSet003[33]={"RPV_2018_smu200_neu180_ctau003","RPV_2018_smu250_neu180_ctau003","RPV_2018_smu250_neu200_ctau003",
// "RPV_2018_smu250_neu230_ctau003","RPV_2018_smu300_neu180_ctau003","RPV_2018_smu300_neu200_ctau003","RPV_2018_smu300_neu280_ctau003",
// "RPV_2018_smu350_neu180_ctau003","RPV_2018_smu350_neu200_ctau003","RPV_2018_smu350_neu250_ctau003","RPV_2018_smu350_neu300_ctau003","RPV_2018_smu350_neu330_ctau003",
// "RPV_2018_smu400_neu180_ctau003","RPV_2018_smu400_neu200_ctau003","RPV_2018_smu400_neu250_ctau003","RPV_2018_smu400_neu300_ctau003","RPV_2018_smu400_neu350_ctau003",
// "RPV_2018_smu400_neu380_ctau003","RPV_2018_smu450_neu180_ctau003","RPV_2018_smu450_neu200_ctau003","RPV_2018_smu450_neu250_ctau003","RPV_2018_smu450_neu300_ctau003",
// "RPV_2018_smu450_neu350_ctau003","RPV_2018_smu450_neu400_ctau003","RPV_2018_smu450_neu430_ctau003","RPV_2018_smu500_neu180_ctau003","RPV_2018_smu500_neu200_ctau003",
// "RPV_2018_smu500_neu250_ctau003","RPV_2018_smu500_neu300_ctau003","RPV_2018_smu500_neu350_ctau003",/*"RPV_2018_smu500_neu400_ctau003",*/"RPV_2018_smu500_neu450_ctau003",
// "RPV_2018_smu500_neu480_ctau003"};

//  for (int i = 0 ; i< 32 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet003[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TrackAna* t = new TrackAna(&c);
//        t->Loop(Prod, SignalSet003[i]);
//     }

//       TString SignalSet010[32]={"RPV_2018_smu200_neu180_ctau010","RPV_2018_smu250_neu180_ctau010","RPV_2018_smu250_neu200_ctau010",
// "RPV_2018_smu250_neu230_ctau010","RPV_2018_smu300_neu180_ctau010","RPV_2018_smu300_neu200_ctau010","RPV_2018_smu300_neu280_ctau010",
// "RPV_2018_smu350_neu180_ctau010","RPV_2018_smu350_neu200_ctau010","RPV_2018_smu350_neu250_ctau010","RPV_2018_smu350_neu300_ctau010","RPV_2018_smu350_neu330_ctau010",
// "RPV_2018_smu400_neu180_ctau010","RPV_2018_smu400_neu200_ctau010","RPV_2018_smu400_neu250_ctau010","RPV_2018_smu400_neu300_ctau010","RPV_2018_smu400_neu350_ctau010",
// "RPV_2018_smu400_neu380_ctau010","RPV_2018_smu450_neu180_ctau010","RPV_2018_smu450_neu200_ctau010","RPV_2018_smu450_neu250_ctau010","RPV_2018_smu450_neu300_ctau010",
// "RPV_2018_smu450_neu350_ctau010","RPV_2018_smu450_neu400_ctau010","RPV_2018_smu450_neu430_ctau010","RPV_2018_smu500_neu180_ctau010","RPV_2018_smu500_neu200_ctau010",
// "RPV_2018_smu500_neu250_ctau010","RPV_2018_smu500_neu300_ctau010",/*"RPV_2018_smu500_neu350_ctau010","RPV_2018_smu500_neu400_ctau010",*/"RPV_2018_smu500_neu450_ctau010",
// "RPV_2018_smu500_neu480_ctau010"};

//  for (int i = 0 ; i< 31 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet010[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TrackAna* t = new TrackAna(&c);
//        t->Loop(Prod, SignalSet010[i]);
//     }

//       TString SignalSet030[34]={"RPV_2018_smu200_neu180_ctau030","RPV_2018_smu250_neu180_ctau030","RPV_2018_smu250_neu200_ctau030",
// "RPV_2018_smu250_neu230_ctau030","RPV_2018_smu300_neu180_ctau030","RPV_2018_smu300_neu200_ctau030","RPV_2018_smu300_neu280_ctau030",
// "RPV_2018_smu350_neu180_ctau030","RPV_2018_smu350_neu200_ctau030","RPV_2018_smu350_neu250_ctau030","RPV_2018_smu350_neu300_ctau030","RPV_2018_smu350_neu330_ctau030",
// "RPV_2018_smu400_neu180_ctau030","RPV_2018_smu400_neu200_ctau030","RPV_2018_smu400_neu250_ctau030","RPV_2018_smu400_neu300_ctau030","RPV_2018_smu400_neu350_ctau030",
// "RPV_2018_smu400_neu380_ctau030","RPV_2018_smu450_neu180_ctau030","RPV_2018_smu450_neu200_ctau030","RPV_2018_smu450_neu250_ctau030","RPV_2018_smu450_neu300_ctau030",
// "RPV_2018_smu450_neu350_ctau030","RPV_2018_smu450_neu400_ctau030","RPV_2018_smu450_neu430_ctau030","RPV_2018_smu500_neu180_ctau030","RPV_2018_smu500_neu200_ctau030",
// "RPV_2018_smu500_neu250_ctau030","RPV_2018_smu500_neu300_ctau030","RPV_2018_smu500_neu350_ctau030","RPV_2018_smu500_neu400_ctau030","RPV_2018_smu500_neu450_ctau030",
// "RPV_2018_smu500_neu480_ctau030"};

//  for (int i = 0 ; i< 33 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet030[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TrackAna* t = new TrackAna(&c);
//        t->Loop(Prod, SignalSet030[i]);
//     }

//     TString SignalSet100[32]={"RPV_2018_smu200_neu180_ctau100","RPV_2018_smu250_neu180_ctau100","RPV_2018_smu250_neu200_ctau100",
// "RPV_2018_smu250_neu230_ctau100","RPV_2018_smu300_neu180_ctau100",/*"RPV_2018_smu300_neu200_ctau100","RPV_2018_smu300_neu250_ctau100_v3",*/"RPV_2018_smu300_neu280_ctau100",
// "RPV_2018_smu350_neu180_ctau100","RPV_2018_smu350_neu200_ctau100","RPV_2018_smu350_neu250_ctau100","RPV_2018_smu350_neu300_ctau100","RPV_2018_smu350_neu330_ctau100",
// "RPV_2018_smu400_neu180_ctau100","RPV_2018_smu400_neu200_ctau100","RPV_2018_smu400_neu250_ctau100","RPV_2018_smu400_neu300_ctau100","RPV_2018_smu400_neu350_ctau100",
// "RPV_2018_smu400_neu380_ctau100","RPV_2018_smu450_neu180_ctau100","RPV_2018_smu450_neu200_ctau100","RPV_2018_smu450_neu250_ctau100","RPV_2018_smu450_neu300_ctau100",
// "RPV_2018_smu450_neu350_ctau100","RPV_2018_smu450_neu400_ctau100","RPV_2018_smu450_neu430_ctau100","RPV_2018_smu500_neu180_ctau100","RPV_2018_smu500_neu200_ctau100",
// "RPV_2018_smu500_neu250_ctau100","RPV_2018_smu500_neu300_ctau100","RPV_2018_smu500_neu350_ctau100","RPV_2018_smu500_neu400_ctau100","RPV_2018_smu500_neu450_ctau100",
// "RPV_2018_smu500_neu480_ctau100"};

//  for (int i = 0 ; i< 31 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet100[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TrackAna* t = new TrackAna(&c);
//        t->Loop(Prod, SignalSet100[i]);
//     }

//       TString SignalSet300[34]={"RPV_2018_smu200_neu180_ctau300","RPV_2018_smu250_neu180_ctau300","RPV_2018_smu250_neu200_ctau300",
// "RPV_2018_smu250_neu230_ctau300","RPV_2018_smu300_neu180_ctau300","RPV_2018_smu300_neu200_ctau300","RPV_2018_smu300_neu280_ctau300",
// "RPV_2018_smu350_neu180_ctau300","RPV_2018_smu350_neu200_ctau300","RPV_2018_smu350_neu250_ctau300","RPV_2018_smu350_neu300_ctau300","RPV_2018_smu350_neu330_ctau300",
// "RPV_2018_smu400_neu180_ctau300","RPV_2018_smu400_neu200_ctau300","RPV_2018_smu400_neu250_ctau300","RPV_2018_smu400_neu300_ctau300","RPV_2018_smu400_neu350_ctau300",
// "RPV_2018_smu400_neu380_ctau300","RPV_2018_smu450_neu180_ctau300","RPV_2018_smu450_neu200_ctau300","RPV_2018_smu450_neu250_ctau300","RPV_2018_smu450_neu300_ctau300",
// "RPV_2018_smu450_neu350_ctau300","RPV_2018_smu450_neu400_ctau300","RPV_2018_smu450_neu430_ctau300","RPV_2018_smu500_neu180_ctau300","RPV_2018_smu500_neu200_ctau300",
// "RPV_2018_smu500_neu250_ctau300","RPV_2018_smu500_neu300_ctau300","RPV_2018_smu500_neu350_ctau300","RPV_2018_smu500_neu400_ctau300","RPV_2018_smu500_neu450_ctau300",
// "RPV_2018_smu500_neu480_ctau300"};

//  for (int i = 0 ; i< 33 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet300[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TrackAna* t = new TrackAna(&c);
//        t->Loop(Prod, SignalSet300[i]);
//     }

//           TString SignalSet1000[33]={"RPV_2018_smu200_neu180_ctau1000","RPV_2018_smu250_neu180_ctau1000","RPV_2018_smu250_neu200_ctau1000",
// /*"RPV_2018_smu250_neu230_ctau1000",*/"RPV_2018_smu300_neu180_ctau1000","RPV_2018_smu300_neu200_ctau1000","RPV_2018_smu300_neu280_ctau1000",
// "RPV_2018_smu350_neu180_ctau1000","RPV_2018_smu350_neu200_ctau1000","RPV_2018_smu350_neu250_ctau1000","RPV_2018_smu350_neu300_ctau1000","RPV_2018_smu350_neu330_ctau1000",
// "RPV_2018_smu400_neu180_ctau1000","RPV_2018_smu400_neu200_ctau1000","RPV_2018_smu400_neu250_ctau1000","RPV_2018_smu400_neu300_ctau1000","RPV_2018_smu400_neu350_ctau1000",
// "RPV_2018_smu400_neu380_ctau1000","RPV_2018_smu450_neu180_ctau1000","RPV_2018_smu450_neu200_ctau1000","RPV_2018_smu450_neu250_ctau1000","RPV_2018_smu450_neu300_ctau1000",
// "RPV_2018_smu450_neu350_ctau1000","RPV_2018_smu450_neu400_ctau1000","RPV_2018_smu450_neu430_ctau1000","RPV_2018_smu500_neu180_ctau1000","RPV_2018_smu500_neu200_ctau1000",
// "RPV_2018_smu500_neu250_ctau1000","RPV_2018_smu500_neu300_ctau1000","RPV_2018_smu500_neu350_ctau1000","RPV_2018_smu500_neu400_ctau1000","RPV_2018_smu500_neu450_ctau1000",
// "RPV_2018_smu500_neu480_ctau1000"};

//  for (int i = 0 ; i< 32 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet1000[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        TrackAna* t = new TrackAna(&c);
//        t->Loop(Prod, SignalSet1000[i]);
//     }

}
