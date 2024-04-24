{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L LLTopEffi.C+g") ;
  
 if (gROOT->GetClass("LLTopTree")==0) return;
 
 TChain c("FlyingTop/ttree");

TString Sample = "RPV_2018_smu200to500_ctau001to100";
TString Prod = "24_03_2024";
////////////////////////////////////////////////////////////////////////////////

// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/240119_fromPaul/RPV_2018_smu200to500_ctau300.root");
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/231212_fromPaul/RPV/RPV_2018_smu200to500_ctau300.root");
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/231206/RPV_2018_smu200to500_ctau300.root");
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/231212/RPV_2018_smu200to500_ctau300.root");
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/2018_231127/RPV_2018_smu200to500_ctau300.root");
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/2018_231112/RPV_2018_smu200to500_ctau300_noSecInt.root");
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/2018_231120/RPV_2018_smu200to500_ctau300.root");

// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/231206/RPV_2018_smu200to500_ctau300.root");
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/240201_fromPaul/RPV_2018_smu200to500_ctau300.root");
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/240207_skipNegIpTk/RPV_2018_smu200to500_ctau300.root");
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC/240211_swapNegIpTk/RPV_2018_smu200to500_ctau300.root");

// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2017/240208_fromPaul/RPV_2017_smu200to500_ctau300.root"); // BDTtk trained on all signal samples and all SM backgrounds: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/2017_ProdTest/
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240209_fromPaul/RPV_2018_smu200to500_ctau300.root"); // BDTtk trained only on smu350_neu300_ctau300 and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Test1SignalSample10cm/
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240212_fromPaul/RPV_2018_smu200to500_ctau300.root"); // BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240216_fromPaul/RPV_2018_smu200to500_ctau300.root"); // BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/
// c.Add("/opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY/CMSSW_10_6_30/MC_2018/"+Prod+"/"+Sample+".root"); //BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/


c.Add("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+Sample+".root"); //BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/

//  MiniIso
// /opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240311
// RPV_2018_smu200to500_ctau001.root
// RPV_2018_smu200to500_ctau003.root
// RPV_2018_smu200to500_ctau010.root
// RPV_2018_smu200to500_ctau030.root
// RPV_2018_smu200to500_ctau100.root
// RPV_2018_smu200to500_ctau300.root
// RPV_2018_smu200to500_ctau100.root

// "RPV_2018_smu200_neu180_ctau300",
//  "RPV_2018_smu250_neu180_ctau300",
//  "RPV_2018_smu250_neu200_ctau300",
//  "RPV_2018_smu250_neu230_ctau300",
//  "RPV_2018_smu300_neu180_ctau300",
//  "RPV_2018_smu300_neu200_ctau300",
// //  "RPV_2018_smu300_neu250_ctau300",
//  "RPV_2018_smu300_neu280_ctau300",
//  "RPV_2018_smu350_neu180_ctau300",
//  "RPV_2018_smu350_neu200_ctau300",
//  "RPV_2018_smu350_neu250_ctau300",
//  "RPV_2018_smu350_neu300_ctau300",
//  "RPV_2018_smu350_neu330_ctau300",
//  "RPV_2018_smu400_neu180_ctau300",
//  "RPV_2018_smu400_neu200_ctau300",
//  "RPV_2018_smu400_neu250_ctau300",
//  "RPV_2018_smu400_neu300_ctau300",
//  "RPV_2018_smu400_neu350_ctau300",
//  "RPV_2018_smu400_neu380_ctau300",
//  "RPV_2018_smu450_neu180_ctau300",
//  "RPV_2018_smu450_neu200_ctau300",
//  "RPV_2018_smu450_neu250_ctau300",
//  "RPV_2018_smu450_neu300_ctau300",
//  "RPV_2018_smu450_neu350_ctau300",
//  "RPV_2018_smu450_neu400_ctau300",
//  "RPV_2018_smu450_neu430_ctau300",
//  "RPV_2018_smu500_neu180_ctau300",
//  "RPV_2018_smu500_neu200_ctau300",
//  "RPV_2018_smu500_neu250_ctau300",
//  "RPV_2018_smu500_neu300_ctau300",
//  "RPV_2018_smu500_neu350_ctau300",
//  "RPV_2018_smu500_neu400_ctau300",
//  "RPV_2018_smu500_neu450_ctau300",
//  "RPV_2018_smu500_neu480_ctau300"


//-  TkIsoTight
// /opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240310/
// RPV_2018_smu200to500_ctau001.root
// RPV_2018_smu200to500_ctau010.root
// RPV_2018_smu200to500_ctau100.root
// RPV_2018_smu200to500_ctau300.root

// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240309/RPV_2018_smu300_neu200_ctau300.root");

////////////////////////////////////////////////////////////////////////////////

 LLTopTree* t = new LLTopTree(&c);

// t->Loop(0,180.,20.,2999.,0.0,10.,0.,1,"output/h_RPV_2018_smu200to500_ctau300_240201_2vtxloose.root","");
// t->Loop(0,180.,20.,2999.,0.0,10.,0.,1,"output/h_RPV_2018_smu200to500_ctau300_240211_swapNegIpTk_2vtx.root","");

t->Loop(0,180.,20.,2999.,0.0,10.,0.,1,"h_"+Sample+"_"+Prod+".root","");
// t->Loop(0,180.,20.,2999.,0.0,10.,0.,1,"output/h_RPV_2018_smu200to500_ctau300_240309_noCloseVtx.root","");

////////////////////////////////////////////////////////////////////////////////

}
