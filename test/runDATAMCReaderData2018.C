{
//  gROOT->ProcessLine(".L DATAMCReader.C+");
  std::vector<TString > systlist;
  systlist.push_back("");
  int YEAR = 2018;
  TString Prod = "DATA_EMU_2018_31_10_2024";// Signal : PROD_CSI_10_06_2024 // BKG : DATAMC2018_EMU_10_06_2024
  bool isPostAPV = false;
  bool isMC = false;
  bool DoubleMuon = false;
  int channel = 0;
  // Emu = 0 ; SingleMuon = 1; DiMuon = 2
  TChain c("ttree");
 //----DATA----//

 TString DataSet[1]={"MuonEG-UL2018_MiniAODv2_GT36-v1"};
          

 for (int i = 0; i< 1 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniDATAMC_"+DataSet[i]+".root";
            // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+DataSet[i]+".root";
      c.Reset();
      c.Add(Path);
      DATAMCReader* t = new DATAMCReader(&c,DataSet[i],systlist);
      float mean = 1;
      t->Loop(isMC,Prod,DataSet[i],false,YEAR, isPostAPV, mean,channel,DoubleMuon,systlist);
          }
}
