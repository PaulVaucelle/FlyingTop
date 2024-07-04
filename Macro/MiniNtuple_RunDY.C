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
  // TString BKGSet[3]={"DYJetsToLL_M10to50","DYJetsToLL_M_50_v1",,"DYJetsToLL_M_50_v2"
  // };

     // //--------------------Background Mumu MiniNtuples ------------//
  TString BKGSet[2]={"DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"};

 for (int i = 0 ; i< 2 ; i++) 
    {
            // /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024
      // /opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/

      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(BKGSet[i],Prod,Signal);
    }

////////////////////////////////////////////////////////////////////////////////

}
