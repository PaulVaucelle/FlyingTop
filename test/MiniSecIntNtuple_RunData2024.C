{
 gROOT->Reset() ; 

 // Compile user's analysis class //
  //  gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniSecIntNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniSecIntNtuple")==0) return;
 
 TChain c("FlyingTop/ttree");

TString Prod = "SecInt";
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////
// /opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY
// /Data_2016pre_emu/240603_bdt94/Data*.root
// /Data_2016post_emu/240603_bdt94/Data*.root
// /Data_2017_emu/240603_bdt94/Data*.root
// //--------------------DATA 2018 MiniSecIntNtuples ------------//

  TString BKGSet[5]={"Data_2024B",  "Data_2024C",  "Data_2024D" , "Data_2024E",  "Data_2024F" };

// "emu_2018A",
// "emu_2018B",
 for (int i = 2 ; i<5 ; i++) 
    {
      TString Path = "/opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY/Data_2024_emu/240819/"+BKGSet[i]+".root";
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+BKGSet[i]+".root";

      c.Reset();
      c.Add(Path);
       MiniSecIntNtuple* t = new MiniSecIntNtuple(&c);
       t->Loop(BKGSet[i],Prod,Signal);
    }

}
