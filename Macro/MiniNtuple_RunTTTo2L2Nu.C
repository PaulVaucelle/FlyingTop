{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniNtuple")==0) return;
 
 TChain c("FlyingTop/ttree");



TString Prod = "DATAMC2018_EMU_10_06_2024";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////

// //--------------------Background MiniNtuples ------------//
  // TString BKGSet[1]={"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"
  // };

  TString BKGSet[1]={"TTTo2L2Nu"
  };
     
 for (int i = 0 ; i< 1 ; i++) 
    {
            // /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024
      // /opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/

      TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(BKGSet[i],Prod,Signal);
    }

////////////////////////////////////////////////////////////////////////////////

}
