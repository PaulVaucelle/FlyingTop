{
//   gROOT->ProcessLine(".L ../HistogramManager.C+");
//  gROOT->ProcessLine(".L TreeABCDReader.C+");
  std::vector<TString > systlist;
  systlist.push_back("");

  int YEAR = 2022;
  TString Prod = "MC_MUMU_2022_23_04_2025"; 
  bool isPostAPV = false;
  bool Signal = false;
  bool SameSign = false;
  bool Forward = false; 
  bool CorrectCorrelation = false;
  bool DoubleMuon = true;
  bool isMC = true;
  int mixing = 0; //-1 : Full left, 0 LR+RL , 1 Full RIght
  int Channel = 2; // 0 : EMu, 1 : SM, 2 : DM

  TChain c("ttree");


  float WGT[12]={
    46928.5,26780,0.869,3.806,3.799,0.07107,0.1659,11.79,7.576,6.831,81.1,337
  };

  TString BKGSet[12]={

    "DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
    "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",

    "TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
    "TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
    "TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",

    "TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8",
    "TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8",

    "WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
    "WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8",
    "ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8",

    "TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
    "TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8"

  };

   for (int i = 0 ; i< 12 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
              Prod = "MC_MUMU_2022_23_04_2025";
              if (systlist.size() == 0) return;
              else if (systlist[j] == "JECUp")
                {
                  Prod = "MC_MUMU_2022_23_04_2025_JECUp";
                }
              else if (systlist[j] == "JECDown")
                {
                  Prod = "MC_MUMU_2022_23_04_2025_JECDown";
                }
              else if (systlist[j] == "JERUp")
                {
                  Prod = "MC_MUMU_2022_23_04_2025_JERUp";
                }
              else if (systlist[j] == "JERDown")
                {
                  Prod = "MC_MUMU_2022_23_04_2025_JERDown";
                }
              else
                {
                  Prod = "MC_MUMU_2022_23_04_2025";
                }
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeABCDReader* t = new TreeABCDReader(&c,Prod, BKGSet[i],systlist[j]);
      // W.push_back(t->MeanGenWeight(BKGSet[i], Prod));//
      double mean = WGT[i];//W[i];//;
      t->Loop(isMC, Prod, BKGSet[i],false, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
            delete t;
        }
    }
    

}
