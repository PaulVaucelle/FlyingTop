{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L ABCD.C+g") ;
  
 if (gROOT->GetClass("ABCD")==0) return;
 
 TChain c("FlyingTop/ttree");

TString Prod = "24_03_2024";
TString Sample = "RPV_2018_smu200to500_ctau1000";
TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+Sample+".root";
bool Signal = false;
////////////////////////////////////////////////////////////////////////////////

// RPV_2018_smu200to500_ctau100.root
// RPV_2018_smu200to500_ctau300.root
// RPV_2018_smu200to500_ctau1000.root
// RPV_2018_smu200to500_ctau010.root
// RPV_2018_smu200to500_ctau030.root
// RPV_2018_smu200to500_ctau001.root
// RPV_2018_smu200to500_ctau003.root

// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2017/240208_fromPaul/RPV_2017_smu200to500_ctau300.root"); // BDTtk trained on all signal samples and all SM backgrounds: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/2017_ProdTest/
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240209_fromPaul/RPV_2018_smu200to500_ctau300.root"); // BDTtk trained only on smu350_neu300_ctau300 and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Test1SignalSample10cm/
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240212_fromPaul/RPV_2018_smu200to500_ctau300.root"); // BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/
// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240216_fromPaul/RPV_2018_smu200to500_ctau300.root"); // BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/
// c.Add("/opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240311/RPV_2018_smu250_neu200_ctau100.root"); //BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/

// RPV_2018_smu250_neu200_ctau100
// RPV_2018_smu300_neu180_ctau100
// RPV_2018_smu400_neu300_ctau100
// RPV_2018_smu500_neu350_ctau100

c.Add(Path); //BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/
// c.Add("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/BDT090_Mini/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root"); //BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/
// c.Add("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/BDT090_Mini/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root"); //BDTtk trained on all ctau300 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/

// DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8
// DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8
// TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8


//  MiniIso
// /opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240311
// RPV_2018_smu200to500_ctau001.root
// RPV_2018_smu200to500_ctau010.root
// RPV_2018_smu200to500_ctau100.root
// RPV_2018_smu200to500_ctau300.root

//-  TkIsoTight
// /opt/sbg/cms/ui2_data1/blochd/NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240310/
// RPV_2018_smu200to500_ctau001.root
// RPV_2018_smu200to500_ctau010.root
// RPV_2018_smu200to500_ctau100.root
// RPV_2018_smu200to500_ctau300.root

// c.Add("../NTUPLES_FLY/CMSSW_10_6_30/MC_2018/240309/RPV_2018_smu300_neu200_ctau300.root");

////////////////////////////////////////////////////////////////////////////////

 ABCD* t = new ABCD(&c);

// t->Loop(0,180.,20.,2999.,0.0,10.,0.,1,"output/h_RPV_2018_smu200to500_ctau300_240201_2vtxloose.root","");
// t->Loop(0,180.,20.,2999.,0.0,10.,0.,1,"output/h_RPV_2018_smu200to500_ctau300_240211_swapNegIpTk_2vtx.root","");

t->Loop(0,180.,20.,2999.,0.0,10.,0.,1,"./"+Prod+"/h_"+Sample+"_REF.root","",Sample,Prod,Signal);
// t->Loop(0,180.,20.,2999.,0.0,10.,0.,1,"output/h_RPV_2018_smu200to500_ctau300_240309_noCloseVtx.root","");

////////////////////////////////////////////////////////////////////////////////

}
