{
 gROOT->Reset() ; 

 // Compile user's analysis class //
  //  gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniDATAMCNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniDATAMCNtuple")==0) return;
 
 TChain c("FlyingTop/ttree");
// TString Prod2018 = "MC_EMU_2018_31_10_2024_v2";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
// TString Prod2017 = "MC_EMU_2017_31_10_2024";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
// TString Prod2016POST = "MC_EMU_2016POST_31_10_2024";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
// TString Prod2016PRE = "MC_EMU_2016PRE_31_10_2024";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024

TString Prod2018 = "MC_MUMU_2018_19_08_2024";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
TString Prod2017 = "MC_MUMU_2017_19_08_2024";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
TString Prod2016POST = "MC_MUMU_2016POST_19_08_2024";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
TString Prod2016PRE = "MC_MUMU_2016PRE_19_08_2024";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024

bool Signal = false;
////////////////////////////////////////////////////////////////////////////////

// //--------------------Background MiniDATAMCNtuples ------------//
  // TString BKGSet[1]={"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"
  // };

    TString BKGSet[1]={"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_eleIPcut_EMu"
  };


     
 for (int i = 0 ; i< 1 ; i++) 
    {
            // /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024
      // /opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/

      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/"+BKGSet[i]+".root";
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TT_eleIPcut_EMu/"+BKGSet[i]+".root";

      // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
       t->Loop(BKGSet[i],Prod2018,Signal);
    }



//  for (int i = 0 ; i< 1 ; i++) 
//     {
//             // /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024
//       // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";

//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2017+"/"+BKGSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
//        t->Loop(BKGSet[i],Prod2017,Signal);
//     }

    //  for (int i = 0 ; i< 1 ; i++) 
    // {
    //         // /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024
    //   // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";

    //   TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2016POST+"/"+BKGSet[i]+".root";
    //   c.Reset();
    //   c.Add(Path);
    //    MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
    //    t->Loop(BKGSet[i],Prod2016POST,Signal);
    // }

  //  for (int i = 0 ; i< 1 ; i++) 
  //   {
  //           // /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PROD_CSI_10_06_2024
  //     // TString Path = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+BKGSet[i]+".root";

  //     TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2016PRE+"/"+BKGSet[i]+".root";
  //     c.Reset();
  //     c.Add(Path);
  //      MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
  //      t->Loop(BKGSet[i],Prod2016PRE,Signal);
  //   }

////////////////////////////////////////////////////////////////////////////////
}

