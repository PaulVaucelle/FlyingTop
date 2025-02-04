{//  gROOT->ProcessLine(".L DATAMCReader.C+");
  std::vector<TString > systlist;
  systlist.push_back("");
  TString Prod2018 = "MC_EMU_2018_31_10_2024_v2";//
  TString Prod2017 = "MC_EMU_2017_31_10_2024";//
  TString Prod2016POST = "MC_EMU_2016POST_31_10_2024";// 
  TString Prod2016PRE = "MC_EMU_2016PRE_31_10_2024";// 
    int YEAR = 2018;
  bool isPostAPV = false;
    bool isMC = true;
bool DoubleMuon = true;
  int channel = 2;
  TChain c("ttree");//FlyingTop/
   // /mc mumu
  //  RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11
  //  RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17
  // RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9
  // RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1
  TString BKGSet[5]={
                      /*"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8","TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8","TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8,"*/
                      "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8","ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
                      /*"DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",*/
                      "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8","WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8","ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8",
  
  };
  std::vector<double> W;
  W={32.4473,32.5092,11.1262,15.3793,8.51008};
 for (int i = 0 ; i< 5 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/MiniDATAMC_"+BKGSet[i]+".root";
            // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
      DATAMCReader* t = new DATAMCReader(&c,BKGSet[i],systlist);
      // W.push_back(t->MeanGenWeight(BKGSet[i], Prod2018));//
      double mean = W[i];//W[i];//;
      t->Loop(isMC,Prod2018,BKGSet[i],false,YEAR, isPostAPV, mean,channel,DoubleMuon,systlist);
    }
std::cout<<"2018 done"<<std::endl;
// W.clear();
//                 // --                        2017                      --//
//  for (int i = 0 ; i< 8 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2017+"/"+BKGSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//       DATAMCReader* t = new DATAMCReader(&c,BKGSet[i],systlist);
//       W.push_back(t->MeanGenWeight(BKGSet[i], Prod2017));//
//       double mean = W[i];//W[i];//;
//       t->Loop(isMC,Prod2017,BKGSet[i],false,2017, isPostAPV, mean,channel,DoubleMuon,systlist);
//     }

// W.clear();
//                   // --                        2016Post                     --//
//      for (int i = 0 ; i< 8 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2016POST+"/"+BKGSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//       DATAMCReader* t = new DATAMCReader(&c,BKGSet[i],systlist);
//       W.push_back(t->MeanGenWeight(BKGSet[i], Prod2016POST));//
//       double mean = W[i];//W[i];//;
//       t->Loop(isMC,Prod2016POST,BKGSet[i],false,2016, true, mean,channel,DoubleMuon,systlist);
//     }

// W.clear();
//                 // --                        2016PRE                      --//
//      for (int i = 0 ; i< 8 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2016PRE+"/"+BKGSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//       DATAMCReader* t = new DATAMCReader(&c,BKGSet[i],systlist);
//       W.push_back(t->MeanGenWeight(BKGSet[i], Prod2016PRE));//
//       double mean = W[i];//W[i];//;
//       t->Loop(isMC,Prod2016PRE,BKGSet[i],false,2016, isPostAPV, mean,channel,DoubleMuon,systlist);
//     }
}
