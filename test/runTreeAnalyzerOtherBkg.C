{
//   gROOT->ProcessLine(".L HistogramManager.C");
//  gROOT->ProcessLine(".L TreeAnalyzer.C");
TString Production = "Test1SignalSample10cm";
  TString BKGSet[8]={  "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8","ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
  "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8", "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8", "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8",
  "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8", "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8", "TTWW_TuneCP5_13TeV-madgraph-pythia8"
  };
 
              
 for (int i = 0 ; i< 8 ; i++) 
    {
      TTree* treeB =0;
      TreeAnalyzer * treeB_ = new TreeAnalyzer(treeB, BKGSet[i],Production);
      treeB_->Loop(BKGSet[i],Production);
      delete treeB_;
      delete treeB;
    }

}
