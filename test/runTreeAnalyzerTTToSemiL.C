{
//   gROOT->ProcessLine(".L HistogramManager.C");
//  gROOT->ProcessLine(".L TreeAnalyzer.C");
TString Production = "Test1SignalSample10cm";
  TString BKGSet[1]={"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"
  };
 
              
 for (int i = 0 ; i< 1 ; i++) 
    {
      TTree* treeB =0;
      TreeAnalyzer * treeB_ = new TreeAnalyzer(treeB, BKGSet[i],Production);
      treeB_->Loop(BKGSet[i],Production);
      delete treeB_;
      delete treeB;
    }
}
