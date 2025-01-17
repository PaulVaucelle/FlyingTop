{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniNtuple")==0) return;
 
 TChain c("FlyingTop/ttree");

TString Year = "2018";
TString Prod = "Signal_2018_JECDown";//Signal_2018

bool Signal = true;
////////////////////////////////////////////////////////////////////////////////
   

c

// // //////////////////////////////////////////////////////////////////////////////
// TString SignalSet[5]={"RPV_"+Year+"_smu300_neu250_ctau001","RPV_"+Year+"_smu300_neu250_ctau003","RPV_"+Year+"_smu300_neu250_ctau300","RPV_"+Year+"_smu300_neu250_ctau1000","RPV_"+Year+"_smu350_neu300_ctau003"};


//  for (int i = 0 ; i< 2 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet[i]+".root";
//       c.Reset();
//       c.Add(Path);
//        MiniNtuple* t = new MiniNtuple(&c);
//        t->Loop(SignalSet[i],Prod,Signal);
//     }

}
