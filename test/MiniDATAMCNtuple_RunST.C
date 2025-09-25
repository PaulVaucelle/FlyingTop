{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniDATAMCNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniDATAMCNtuple")==0) return;
 
 TChain c("FlyingTop/ttree");

TString Prod2018 = "MC_EMU_2018_30_03_2025";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
TString Prod2017 = "MC_EMU_2017_30_03_2025";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
TString Prod2016POST = "MC_EMU_2016POST_30_03_2025";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024
TString Prod2016PRE = "MC_EMU_2016PRE_30_03_2025";//PROD_CSI_10_06_2024 // DATAMC2018_EMU_10_06_2024

// TString Prod2018 = "MC_EMU_03_02_2025_RoccorDown";//JECUp // JECDown // JERUp // JERDown

bool Signal = false;

     // //--------------------Background Mumu MiniDATAMCNtuples ------------//
  TString BKGSet[2]={"ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8","ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8"
  };

//  for (int i = 0 ; i< 2 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/"+BKGSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
//        t->Loop(BKGSet[i],Prod2018,Signal);
//        delete t;
//     }

 for (int i = 0 ; i< 2 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2017+"/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
       t->Loop(BKGSet[i],Prod2017,Signal);
          delete t;
    }

     for (int i = 0 ; i< 2 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2016POST+"/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
       t->Loop(BKGSet[i],Prod2016POST,Signal);
              delete t;
    }

   for (int i = 0 ; i< 2 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2016PRE+"/"+BKGSet[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniDATAMCNtuple* t = new MiniDATAMCNtuple(&c);
       t->Loop(BKGSet[i],Prod2016PRE,Signal);
             delete t;
    }
////////////////////////////////////////////////////////////////////////////////

}
