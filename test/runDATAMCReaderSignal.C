{
//  gROOT->ProcessLine(".L DATAMCReader.C+");
  std::vector<TString > systlist;
  systlist.push_back("");
  int YEAR = 2018;
  TString Prod = "Signal_2018";// Signal : PROD_CSI_10_06_2024 // BKG : DATAMC2018_EMU_10_06_2024
  bool isPostAPV = false;
  bool DoubleMuon = false;
  int Channel = 0;
  bool isMC = true;
  TChain c("ttree");//FlyingTop/
             

// // //--------------------Test Signal Ntuples ------------//
TString SignalSet[4]={"RPV_2018_smu300_neu180_ctau100","RPV_2018_smu400_neu300_ctau100","RPV_2018_smu200_neu180_ctau100","RPV_2018_smu500_neu350_ctau100"};

 for (int i = 0 ; i< 2 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniDATAMC_"+SignalSet[i]+".root";
            // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
      DATAMCReader* t = new DATAMCReader(&c,SignalSet[i],systlist);
      float mean = 1;
      t->Loop(isMC,Prod,SignalSet[i],true,YEAR, isPostAPV, mean,Channel ,DoubleMuon, systlist);
      
    }

//       TString SignalSet001[32]={"RPV_2018_smu200_neu180_ctau001","RPV_2018_smu250_neu180_ctau001","RPV_2018_smu250_neu200_ctau001",
// "RPV_2018_smu250_neu230_ctau001","RPV_2018_smu300_neu180_ctau001","RPV_2018_smu300_neu200_ctau001","RPV_2018_smu300_neu280_ctau001",
// "RPV_2018_smu350_neu180_ctau001","RPV_2018_smu350_neu200_ctau001","RPV_2018_smu350_neu250_ctau001","RPV_2018_smu350_neu300_ctau001","RPV_2018_smu350_neu330_ctau001",
// "RPV_2018_smu400_neu180_ctau001","RPV_2018_smu400_neu200_ctau001","RPV_2018_smu400_neu250_ctau001","RPV_2018_smu400_neu300_ctau001","RPV_2018_smu400_neu350_ctau001",
// "RPV_2018_smu400_neu380_ctau001","RPV_2018_smu450_neu180_ctau001","RPV_2018_smu450_neu200_ctau001","RPV_2018_smu450_neu250_ctau001","RPV_2018_smu450_neu300_ctau001",
// "RPV_2018_smu450_neu350_ctau001","RPV_2018_smu450_neu400_ctau001","RPV_2018_smu450_neu430_ctau001","RPV_2018_smu500_neu180_ctau001","RPV_2018_smu500_neu200_ctau001",
// "RPV_2018_smu500_neu250_ctau001","RPV_2018_smu500_neu300_ctau001","RPV_2018_smu500_neu350_ctau001",/*"RPV_2018_smu500_neu400_ctau001","RPV_2018_smu500_neu450_ctau001",*/
// "RPV_2018_smu500_neu480_ctau001"};
//  //issue  with RPV_2018_smu500_neu400_ctau001 500 450 001
//  for (int i = 0 ; i< 31 ; i++) 
//     {
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/"+SignalSet001[i]+".root";
      //       // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      // c.Reset();
      // c.Add(Path);
      // DATAMCReader* t = new DATAMCReader(&c,SignalSet[i],systlist);
      // float mean = t->MeanGenWeight();
      // t->Loop(isMC,ProdSignalSet[i],true,YEAR, isPostAPV, mean,systlist);
//     }



//       TString SignalSet003[33]={"RPV_2018_smu200_neu180_ctau003","RPV_2018_smu250_neu180_ctau003","RPV_2018_smu250_neu200_ctau003",
// "RPV_2018_smu250_neu230_ctau003","RPV_2018_smu300_neu180_ctau003","RPV_2018_smu300_neu200_ctau003","RPV_2018_smu300_neu280_ctau003",
// "RPV_2018_smu350_neu180_ctau003","RPV_2018_smu350_neu200_ctau003","RPV_2018_smu350_neu250_ctau003","RPV_2018_smu350_neu300_ctau003","RPV_2018_smu350_neu330_ctau003",
// "RPV_2018_smu400_neu180_ctau003","RPV_2018_smu400_neu200_ctau003","RPV_2018_smu400_neu250_ctau003","RPV_2018_smu400_neu300_ctau003","RPV_2018_smu400_neu350_ctau003",
// "RPV_2018_smu400_neu380_ctau003","RPV_2018_smu450_neu180_ctau003","RPV_2018_smu450_neu200_ctau003","RPV_2018_smu450_neu250_ctau003","RPV_2018_smu450_neu300_ctau003",
// "RPV_2018_smu450_neu350_ctau003","RPV_2018_smu450_neu400_ctau003","RPV_2018_smu450_neu430_ctau003","RPV_2018_smu500_neu180_ctau003","RPV_2018_smu500_neu200_ctau003",
// "RPV_2018_smu500_neu250_ctau003","RPV_2018_smu500_neu300_ctau003","RPV_2018_smu500_neu350_ctau003",/*"RPV_2018_smu500_neu400_ctau003",*/"RPV_2018_smu500_neu450_ctau003",
// "RPV_2018_smu500_neu480_ctau003"};
//  //issue RPV_2018_smu500_neu400_ctau003
//  for (int i = 0 ; i< 32 ; i++) 
//     {
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/"+SignalSet003[i]+".root";
      //       // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      // c.Reset();
      // c.Add(Path);
      // DATAMCReader* t = new DATAMCReader(&c,SignalSet[i],systlist);
      // float mean = t->MeanGenWeight();
      // t->Loop(isMC,ProdSignalSet[i],true,YEAR, isPostAPV, mean,systlist);
//     }

//       TString SignalSet010[32]={"RPV_2018_smu200_neu180_ctau010","RPV_2018_smu250_neu180_ctau010","RPV_2018_smu250_neu200_ctau010",
// "RPV_2018_smu250_neu230_ctau010","RPV_2018_smu300_neu180_ctau010","RPV_2018_smu300_neu200_ctau010","RPV_2018_smu300_neu280_ctau010",
// "RPV_2018_smu350_neu180_ctau010","RPV_2018_smu350_neu200_ctau010","RPV_2018_smu350_neu250_ctau010","RPV_2018_smu350_neu300_ctau010","RPV_2018_smu350_neu330_ctau010",
// "RPV_2018_smu400_neu180_ctau010","RPV_2018_smu400_neu200_ctau010","RPV_2018_smu400_neu250_ctau010","RPV_2018_smu400_neu300_ctau010","RPV_2018_smu400_neu350_ctau010",
// "RPV_2018_smu400_neu380_ctau010","RPV_2018_smu450_neu180_ctau010","RPV_2018_smu450_neu200_ctau010","RPV_2018_smu450_neu250_ctau010","RPV_2018_smu450_neu300_ctau010",
// "RPV_2018_smu450_neu350_ctau010","RPV_2018_smu450_neu400_ctau010","RPV_2018_smu450_neu430_ctau010","RPV_2018_smu500_neu180_ctau010","RPV_2018_smu500_neu200_ctau010",
// "RPV_2018_smu500_neu250_ctau010","RPV_2018_smu500_neu300_ctau010",/*"RPV_2018_smu500_neu350_ctau010","RPV_2018_smu500_neu400_ctau010",*/"RPV_2018_smu500_neu450_ctau010",
// "RPV_2018_smu500_neu480_ctau010"};
//  //issue with 500 350 10 // 500 400 10
//  for (int i = 0 ; i< 31 ; i++) 
//     {
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/"+SignalSet010[i]+".root";
      //       // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      // c.Reset();
      // c.Add(Path);
      // DATAMCReader* t = new DATAMCReader(&c,SignalSet[i],systlist);
      // float mean = t->MeanGenWeight();
      // t->Loop(isMC,ProdSignalSet[i],true,YEAR, isPostAPV, mean,systlist);
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
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/"+SignalSet030[i]+".root";
      //       // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      // c.Reset();
      // c.Add(Path);
      // DATAMCReader* t = new DATAMCReader(&c,SignalSet[i],systlist);
      // float mean = t->MeanGenWeight();
      // t->Loop(isMC,ProdSignalSet[i],true,YEAR, isPostAPV, mean,systlist);
//     }

//   TString SignalSet100[32]={"RPV_2018_smu200_neu180_ctau100","RPV_2018_smu250_neu180_ctau100","RPV_2018_smu250_neu200_ctau100",
// "RPV_2018_smu250_neu230_ctau100","RPV_2018_smu300_neu180_ctau100",/*"RPV_2018_smu300_neu200_ctau100","RPV_2018_smu300_neu250_ctau100_v3",*/"RPV_2018_smu300_neu280_ctau100",
// "RPV_2018_smu350_neu180_ctau100","RPV_2018_smu350_neu200_ctau100","RPV_2018_smu350_neu250_ctau100","RPV_2018_smu350_neu300_ctau100","RPV_2018_smu350_neu330_ctau100",
// "RPV_2018_smu400_neu180_ctau100","RPV_2018_smu400_neu200_ctau100","RPV_2018_smu400_neu250_ctau100","RPV_2018_smu400_neu300_ctau100","RPV_2018_smu400_neu350_ctau100",
// "RPV_2018_smu400_neu380_ctau100","RPV_2018_smu450_neu180_ctau100","RPV_2018_smu450_neu200_ctau100","RPV_2018_smu450_neu250_ctau100","RPV_2018_smu450_neu300_ctau100",
// "RPV_2018_smu450_neu350_ctau100","RPV_2018_smu450_neu400_ctau100","RPV_2018_smu450_neu430_ctau100","RPV_2018_smu500_neu180_ctau100","RPV_2018_smu500_neu200_ctau100",
// "RPV_2018_smu500_neu250_ctau100","RPV_2018_smu500_neu300_ctau100","RPV_2018_smu500_neu350_ctau100","RPV_2018_smu500_neu400_ctau100","RPV_2018_smu500_neu450_ctau100",
// "RPV_2018_smu500_neu480_ctau100"};
 
//  for (int i = 0 ; i< 31 ; i++) 
//     {
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/"+SignalSet100[i]+".root";
      //       // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      // c.Reset();
      // c.Add(Path);
      // DATAMCReader* t = new DATAMCReader(&c,SignalSet[i],systlist);
      // float mean = t->MeanGenWeight();
      // t->Loop(isMC,ProdSignalSet[i],true,YEAR, isPostAPV, mean,systlist);
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
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/"+SignalSet300[i]+".root";
      //       // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      // c.Reset();
      // c.Add(Path);
      // DATAMCReader* t = new DATAMCReader(&c,SignalSet[i],systlist);
      // float mean = t->MeanGenWeight();
      // t->Loop(isMC,ProdSignalSet[i],true,YEAR, isPostAPV, mean,systlist);
//     }


//       TString SignalSet1000[33]={"RPV_2018_smu200_neu180_ctau1000","RPV_2018_smu250_neu180_ctau1000","RPV_2018_smu250_neu200_ctau1000",
// /*"RPV_2018_smu250_neu230_ctau1000",*/"RPV_2018_smu300_neu180_ctau1000","RPV_2018_smu300_neu200_ctau1000","RPV_2018_smu300_neu280_ctau1000",
// "RPV_2018_smu350_neu180_ctau1000","RPV_2018_smu350_neu200_ctau1000","RPV_2018_smu350_neu250_ctau1000","RPV_2018_smu350_neu300_ctau1000","RPV_2018_smu350_neu330_ctau1000",
// "RPV_2018_smu400_neu180_ctau1000","RPV_2018_smu400_neu200_ctau1000","RPV_2018_smu400_neu250_ctau1000","RPV_2018_smu400_neu300_ctau1000","RPV_2018_smu400_neu350_ctau1000",
// "RPV_2018_smu400_neu380_ctau1000","RPV_2018_smu450_neu180_ctau1000","RPV_2018_smu450_neu200_ctau1000","RPV_2018_smu450_neu250_ctau1000","RPV_2018_smu450_neu300_ctau1000",
// "RPV_2018_smu450_neu350_ctau1000","RPV_2018_smu450_neu400_ctau1000","RPV_2018_smu450_neu430_ctau1000","RPV_2018_smu500_neu180_ctau1000","RPV_2018_smu500_neu200_ctau1000",
// "RPV_2018_smu500_neu250_ctau1000","RPV_2018_smu500_neu300_ctau1000","RPV_2018_smu500_neu350_ctau1000","RPV_2018_smu500_neu400_ctau1000","RPV_2018_smu500_neu450_ctau1000",
// "RPV_2018_smu500_neu480_ctau1000"};
 
//  for (int i = 0 ; i< 32 ; i++) 
//     {
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/"+SignalSet1000[i]+".root";
      //       // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      // c.Reset();
      // c.Add(Path);
      // DATAMCReader* t = new DATAMCReader(&c,SignalSet[i],systlist);
      // float mean = t->MeanGenWeight();
      // t->Loop(isMC,ProdSignalSet[i],true,YEAR, isPostAPV, mean,systlist);
//     }

}
