{
//  gROOT->ProcessLine(".L DATAMCReader.C+");
  std::vector<TString > systlist;
  systlist.push_back("");
  // systlist.push_back("LumiUp");
  // systlist.push_back("LumiDown");
  // // systlist.push_back("XSUp");
  // // systlist.push_back("XSDown");
  // systlist.push_back("L1Up");
  // systlist.push_back("L1Down");
  // systlist.push_back("TriggerDown");
  // systlist.push_back("TriggerUp");
  // systlist.push_back("MuonIDUp");
  // systlist.push_back("MuonIDDown");
  // systlist.push_back("MuonISOUp");
  // systlist.push_back("MuonISODown");
  // systlist.push_back("EleIDUp");
  // systlist.push_back("EleIDDown");
  // systlist.push_back("EleISOUp");
  // systlist.push_back("EleISODown");
  // systlist.push_back("PUUp");
  // systlist.push_back("PUDown");
  // systlist.push_back("JECUp");
  // systlist.push_back("JECDown");
  // systlist.push_back("JERUp");
  // systlist.push_back("JERDown");
  // systlist.push_back("SFEleUp");
  // systlist.push_back("SFEleDown");
  // systlist.push_back("TopPtUp");
  // systlist.push_back("TopPtDown");
  // systlist.push_back("RoccorUp");
  // systlist.push_back("RoccorDown");
  // systlist.push_back("PDFUp");
  // systlist.push_back("PDFDown");
  // systlist.push_back("ScaleUp");
  // systlist.push_back("ScaleDown");


  TString Prod2018 = "MC_EMU_03_02_2025";//
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
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {

           Prod2018 = "MC_EMU_03_02_2025";
            if (systlist.size() == 0) return;
            else if (systlist[j] == "JECUp")
              {
                Prod2018 = "MC_EMU_03_02_2025_JECUp";
              }
            else if (systlist[j] == "JECDown")
              {
                Prod2018 = "MC_EMU_03_02_2025_JECDown";
              }
            else if (systlist[j] == "JERUp")
              {
                Prod2018 = "MC_EMU_03_02_2025_JERUp";
              }
            else if (systlist[j] == "JERDown")
              {
                Prod2018 = "MC_EMU_03_02_2025_JERDown";
              }
            else if (systlist[j] == "RoccorDown")
              {
                Prod2018 = "MC_EMU_03_02_2025_RoccorDown";
              }
            else 
              {
                Prod2018 = "MC_EMU_03_02_2025";
              }

          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/MiniDATAMC_"+BKGSet[i]+".root";
          c.Reset();
          c.Add(Path);
          DATAMCReader* t = new DATAMCReader(&c,BKGSet[i],systlist[j]);
          // W.push_back(t->MeanGenWeight(BKGSet[i], Prod2018));//
          double mean = 1.;//W[i];//;
          t->Loop(isMC,Prod2018,BKGSet[i],false,2018, isPostAPV, mean,channel,DoubleMuon,systlist[j]);
          delete t;
   
        }

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
