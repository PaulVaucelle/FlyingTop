{
//   gROOT->ProcessLine(".L HistogramManager.C");
//  gROOT->ProcessLine(".L TreeAnalyzer.C");
TString Production = "BDT090";
  TString BKGSet[2]={/*"TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8","TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"
                      ,"ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8","ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",*/
                      "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
  };
 
              
 for (int i = 0 ; i< 2 ; i++) 
    {
      TTree* treeB =0;
      TreeAnalyzer * treeB_ = new TreeAnalyzer(treeB, BKGSet[i],Production);
      treeB_->Loop(BKGSet[i],Production);
      delete treeB_;
      delete treeB;
    }

}
