{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniNtuple")==0) return;
 
 TChain c("FlyingTop/ttree");

TString Prod = "PROD_CSI_10_06_2024";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////

// //--------------------Background Emu MiniNtuples ------------//
  //  TString BKGSet[8]={ 
  //    "ST_tW_top_5f_NoFullyHadronicDecays","ST_tW_antitop_5f_NoFullyHadronicDecays",
  // "WWTo2L2Nu", "WZTo2Q2L_mllmin4p0", "ZZTo2Q2L_mllmin4p0",
  // "ttWJetsToLNu_5f_EWK", "TTZToLL_5f", "TTWW_TuneCP5"
  // "TTToHadronic"

  // };

// //--------------------Background Mumu MiniNtuples ------------//
   TString BKGSet[8]={ 
     "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8","ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",
  "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8", "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8", "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8",
  "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8", "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8", "TTWW_TuneCP5_13TeV-madgraph-pythia8"
  "TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8"

  };
  
     
 for (int i = 0 ; i< 1 ; i++) 
    {
      // /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024
      // /opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(BKGSet[i],Prod,Signal);
    }

}
