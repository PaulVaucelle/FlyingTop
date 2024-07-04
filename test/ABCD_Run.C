{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L ABCD.C+g") ;
  
 if (gROOT->GetClass("ABCD")==0) return;
 
 TChain c("ttree");

TString Prod = "DATAMC2018_EMU_10_06_2024";//PROD_CSI_10_06_2024 for data // DATAMC2018_EMU_10_06_2024 for MC // DATA_MUMU_2018_20_06_2024
// TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+Sample+".root";
// bool Signal = true;
bool SameSign = false;
bool Forward = false; 
bool DoubleMuon = true;
////////////////////////////////////////////////////////////////////////////////
// //----DATA----//
// TString DataSet[4]={"DataA", "DataB","DataC","DataD"};//EMU  Data PROD_CSI_10_06_2024
// TString DataSet[4]={"SM_Run2018A-UL2018_MiniAODv2-v3","SM_Run2018B-UL2018_MiniAODv2-v2","SM_Run2018C-UL2018_MiniAODv2-v2","SM_Run2018D-UL2018_MiniAODv2-v3"};
// // // "DM_Run2018A-UL2018_MiniAODv2-v1","DM_Run2018B-UL2018_MiniAODv2-v1","DM_Run2018C-UL2018_MiniAODv2-v1","DM_Run2018D-UL2018_MiniAODv2-v1"
// //                     // "SM_Run2018A-UL2018_MiniAODv2-v3","SM_Run2018B-UL2018_MiniAODv2-v2","SM_Run2018C-UL2018_MiniAODv2-v2","SM_Run2018D-UL2018_MiniAODv2-v3"};
// // // TString DataSet[4]={"M_Run2018A-UL2018", "M_Run2018B-UL2018","M_Run2018C-UL2018","M_Run2018D-UL2018"};//Mumu Dm+SM 2018 20_06_2024 DATA_MUMU_2018_20_06_2024
//  for (int i = 0 ; i< 4 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+DataSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//       ABCD* t = new ABCD(&c);
//       t->Loop(DataSet[i],Prod,false, SameSign, Forward, DoubleMuon);
//     }


// //--------------------Background MiniNtuples ------------//
//MC mumu
  // TString BKGSet[13]={
  //                     "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8","TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8","TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"
  //                     ,"ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8","ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
  //                     "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
  //                     "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8","WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8","ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8",
  //                     "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8","TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8","TTWW_TuneCP5_13TeV-madgraph-pythia8"
  // };

//mc emu
  TString BKGSet[13]={
                      "TTTo2L2Nu","TTToHadronic","TTToSemiLeptonic"
                      ,"ST_tW_top_5f_NoFullyHadronicDecays","ST_tW_antitop_5f_NoFullyHadronicDecays",
                      "DYJetsToLL_M10to50","DYJetsToLL_M_50_v1",
                      "WWTo2L2Nu","WZTo2Q2L_mllmin4p0","ZZTo2Q2L_mllmin4p0",
                      "ttWJetsToLNu_5f_EWK","TTZToLL_5f","TTWW_TuneCP5"
  };

                            
 for (int i = 0 ; i< 13 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
      ABCD* t = new ABCD(&c);
      t->Loop(BKGSet[i],Prod,false, SameSign, Forward, DoubleMuon);
    }

// ////////////////////////////////////////////////////////////////////////////////


// // ////////////////////////////////////////////////////////////////////////////////
// // //--------------------Test Signal MiniNtuples ------------//
// TString SignalSet[4]={"RPV_2018_smu300_neu180_ctau100","RPV_2018_smu400_neu300_ctau100","RPV_2018_smu200_neu180_ctau100","RPV_2018_smu500_neu350_ctau100"};


//  for (int i = 0 ; i< 4 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        ABCD* t = new ABCD(&c);
//        t->Loop(SignalSet[i],Prod,true, SameSign, Forward, DoubleMuon);
//     }

// //--------------------------------------------------//


// // ////////////////////////////////////////////////////////////////////////////////
// // //--------------------Test DATA-----------//
// TString SignalSet[4]={"RPV_2018_smu300_neu180_ctau100","RPV_2018_smu400_neu300_ctau100","RPV_2018_smu250_neu200_ctau100","RPV_2018_smu500_neu350_ctau100"};


//  for (int i = 0 ; i< 4 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        ABCD* t = new ABCD(&c);
//        t->Loop(SignalSet[i],Prod,true, SameSign, Forward);
//     }

//--------------------------------------------------//


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
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        ABCD* t = new ABCD(&c);
//        t->Loop(SignalSet[i],Prod,Signal, SameSign, Forward);
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
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        ABCD* t = new ABCD(&c);
//        t->Loop(SignalSet[i],Prod,Signal, SameSign, Forward);
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
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        ABCD* t = new ABCD(&c);
//        t->Loop(SignalSet[i],Prod,Signal, SameSign, Forward);
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
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        ABCD* t = new ABCD(&c);
//        t->Loop(SignalSet[i],Prod,Signal, SameSign, Forward);
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
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        ABCD* t = new ABCD(&c);
//        t->Loop(SignalSet[i],Prod,Signal, SameSign, Forward);
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
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        ABCD* t = new ABCD(&c);
//        t->Loop(SignalSet[i],Prod,Signal, SameSign, Forward);
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
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        ABCD* t = new ABCD(&c);
//        t->Loop(SignalSet[i],Prod,Signal, SameSign, Forward);
//     }


////////////////////////////////////////////////////////////////////////////////

}
