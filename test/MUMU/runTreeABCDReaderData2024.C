{
//   gROOT->ProcessLine(".L ../HistogramManager.C+");
//  gROOT->ProcessLine(".L TreeABCDReader.C+");
  std::vector<TString > systlist;
  systlist.push_back("");

  TString ProdEMU = "DATA_MUMU_2024_23_04_2025"; // "DATA_MUMU_2018_19_08_2024"; // "DATA_EMU_2018_19_08_2024";
  int YEAR = 2024;
  bool isPostAPV = false;
  bool Signal = false;
  bool SameSign = false;
  bool Forward = false; 
  bool DoubleMuon = true;
  bool CorrectCorrelation = true;
  bool isMC = false;
  int mixing = 0; //-1 : Full left, 0 LR+RL , 1 Full RIght
  int Channel = 2; // 0 : EMu, 1 : SM, 2 : DM

  TChain c("ttree");
//  // Data 

  TString BKGSetEMU[1]={ "Muon_Run2024"
  };

   for (int i = 0 ; i< 1 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+ProdEMU+"/Mini"+BKGSetEMU[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeABCDReader* t = new TreeABCDReader(&c,ProdEMU, BKGSetEMU[i],systlist[0]);
      float mean = 1.;//t->MeanGenWeight()
      t->Loop(isMC, ProdEMU, BKGSetEMU[i],false, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[0],CorrectCorrelation);
      delete t;
    }
          
    
}
