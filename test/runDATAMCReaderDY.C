{
//  gROOT->ProcessLine(".L DATAMCReader.C+");
  std::vector<TString > systlist;
  systlist.push_back("");
  TString Prod2018 = "MC_EMU_2018_31_10_2024_v2";//
  TString Prod2017 = "MC_EMU_2017_31_10_2024";//
  TString Prod2016POST = "MC_EMU_2016POST_31_10_2024";// 
  TString Prod2016PRE = "MC_EMU_2016PRE_31_10_2024";// 
  int YEAR = 2018;
  bool isPostAPV = false;
  bool isMC = true;
  bool DoubleMuon = false;
  int channel = 0;
  TChain c("ttree");//FlyingTop/
   // /mc mumu
     //  RunIISummer20UL16MiniAODAPVv2-106X_mcRun2_asymptotic_preVFP_v11
  //  RunIISummer20UL16MiniAODv2-106X_mcRun2_asymptotic_v17
  // RunIISummer20UL17MiniAODv2-106X_mc2017_realistic_v9
// RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1
  TString BKGSet[2]={
                      /*"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8","TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8","TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"
                      ,"ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8","ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
                      */"DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
                      /*"WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8","WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8","ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8",
                      "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8","TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8","TTWW_TuneCP5_13TeV-madgraph-pythia8"
  */
  };

                // --                        2018                      --//
  std::vector<double> W;
 for (int i = 0 ; i< 2 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/MiniDATAMC_"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
      DATAMCReader* t = new DATAMCReader(&c,BKGSet[i],systlist);
      W.push_back(t->MeanGenWeight(BKGSet[i], Prod2018));//
      double mean = 1.;//W[i];//;
      t->Loop(isMC,Prod2018,BKGSet[i],false,2018, isPostAPV, mean,channel,DoubleMuon,systlist);
       
    }
 std::cout<<"End of 2018"<<std::endl;
// W.clear();
//                 // --                        2017                      --//
//  for (int i = 0 ; i< 2 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2017+"/"+BKGSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//       DATAMCReader* t = new DATAMCReader(&c,BKGSet[i],systlist);
//       W.push_back(t->MeanGenWeight(BKGSet[i], Prod2017));//
//       double mean = W[i];//W[i];//;
//       t->Loop(isMC,Prod2017,BKGSet[i],false,2017, isPostAPV, mean,channel,DoubleMuon,systlist);
        // std::cout<<"End of 2017"<<std::endl;
//     }

// W.clear();
//                   // --                        2016Post                     --//
//      for (int i = 0 ; i< 2 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2016POST+"/"+BKGSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//       DATAMCReader* t = new DATAMCReader(&c,BKGSet[i],systlist);
//       W.push_back(t->MeanGenWeight(BKGSet[i], Prod2016POST));//
//       double mean = W[i];//W[i];//;
//       t->Loop(isMC,Prod2016POST,BKGSet[i],false,2016, true, mean,channel,DoubleMuon,systlist);
//       
//     }
// std::cout<<"End of 2016POST"<<std::endl;
// W.clear();
//                 // --                        2016PRE                      --//
//      for (int i = 0 ; i< 2 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2016PRE+"/"+BKGSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//       DATAMCReader* t = new DATAMCReader(&c,BKGSet[i],systlist);
//       W.push_back(t->MeanGenWeight(BKGSet[i], Prod2016PRE));//
//       double mean = W[i];//W[i];//;
//       t->Loop(isMC,Prod2016PRE,BKGSet[i],false,2016, isPostAPV, mean,channel,DoubleMuon,systlist);

//     }

//       std::cout<<"End of 2016PRE"<<std::endl;

}
