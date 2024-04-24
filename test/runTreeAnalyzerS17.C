{
//   gROOT->ProcessLine(".L HistogramManager.C");
//  gROOT->ProcessLine(".L TreeAnalyzer.C");

TString Production = "2017_ProdTest";
//   TString BKGSet[10]={/*"TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8","TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"
//                       ,"ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8","ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8",*/
//                       "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8","WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8",
//                       "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8","ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8","ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8","TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8","TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8","TTWW_TuneCP5_13TeV-madgraph-pythia8","TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"
//   };
 
              
//  for (int i = 0 ; i< 10 ; i++) 
//     {
//       TTree* treeB =0;
//       TreeAnalyzer * treeB_ = new TreeAnalyzer(treeB, BKGSet[i],Production, );
//       treeB_->Loop(BKGSet[i],Production,false);
//       delete treeB_;
//       delete treeB;
//     }


//   TString BKGSet[14]={"TTJets_DiLept","TTTo2L2Nu","ST_tW_top_5f_NoFullyHadronicDecays","ST_tW_antitop_5f_NoFullyHadronicDecays",
//                       "DYJetsToLL_M10to50","DYJetsToLL_M50","WWTo2L2Nu",
//                       "WZTo2Q2L_mllmin4p0","ZZTo2Q2L_mllmin4p0","ttWJetsToLNu_5f_EWK","TTZToLL_5f","TTToHadronic","TTWW","TTToSemiLeptonic"
//   };
 
//  for (int i = 0 ; i< 14 ; i++) 
//     {
//       TTree* treeB =0;
//       TreeAnalyzer * treeB_ = new TreeAnalyzer(treeB, BKGSet[i],Production, );
//       treeB_->Loop(BKGSet[i],Production,false);
//       delete treeB_;
//       delete treeB;
//     }

//---45 min for all the signal samples

      TString SignalSet001[32]={"RPV_2017_smu200_neu180_ctau001","RPV_2017_smu250_neu180_ctau001","RPV_2017_smu250_neu200_ctau001",
"RPV_2017_smu250_neu230_ctau001","RPV_2017_smu300_neu180_ctau001","RPV_2017_smu300_neu200_ctau001","RPV_2017_smu300_neu280_ctau001",
"RPV_2017_smu350_neu180_ctau001","RPV_2017_smu350_neu200_ctau001","RPV_2017_smu350_neu250_ctau001","RPV_2017_smu350_neu300_ctau001","RPV_2017_smu350_neu330_ctau001",
"RPV_2017_smu400_neu180_ctau001","RPV_2017_smu400_neu200_ctau001","RPV_2017_smu400_neu250_ctau001","RPV_2017_smu400_neu300_ctau001","RPV_2017_smu400_neu350_ctau001",
"RPV_2017_smu400_neu380_ctau001","RPV_2017_smu450_neu180_ctau001","RPV_2017_smu450_neu200_ctau001","RPV_2017_smu450_neu250_ctau001","RPV_2017_smu450_neu300_ctau001",
"RPV_2017_smu450_neu350_ctau001","RPV_2017_smu450_neu400_ctau001","RPV_2017_smu450_neu430_ctau001","RPV_2017_smu500_neu180_ctau001","RPV_2017_smu500_neu200_ctau001",
"RPV_2017_smu500_neu250_ctau001","RPV_2017_smu500_neu300_ctau001","RPV_2017_smu500_neu350_ctau001",/*"RPV_2017_smu500_neu400_ctau001","RPV_2017_smu500_neu450_ctau001",*/
"RPV_2017_smu500_neu480_ctau001"};
 //issue  with RPV_2017_smu500_neu400_ctau001 500 450 001
 for (int i = 0 ; i< 31 ; i++) 
    {
      TTree* tree =0;
      TreeAnalyzer * tree_ = new TreeAnalyzer(tree,SignalSet001[i],Production);
      tree_->Loop(SignalSet001[i],Production);
      delete tree_;
      delete tree;
    }



      TString SignalSet003[33]={"RPV_2017_smu200_neu180_ctau003","RPV_2017_smu250_neu180_ctau003","RPV_2017_smu250_neu200_ctau003",
"RPV_2017_smu250_neu230_ctau003","RPV_2017_smu300_neu180_ctau003","RPV_2017_smu300_neu200_ctau003","RPV_2017_smu300_neu280_ctau003",
"RPV_2017_smu350_neu180_ctau003","RPV_2017_smu350_neu200_ctau003","RPV_2017_smu350_neu250_ctau003","RPV_2017_smu350_neu300_ctau003","RPV_2017_smu350_neu330_ctau003",
"RPV_2017_smu400_neu180_ctau003","RPV_2017_smu400_neu200_ctau003","RPV_2017_smu400_neu250_ctau003","RPV_2017_smu400_neu300_ctau003","RPV_2017_smu400_neu350_ctau003",
"RPV_2017_smu400_neu380_ctau003","RPV_2017_smu450_neu180_ctau003","RPV_2017_smu450_neu200_ctau003","RPV_2017_smu450_neu250_ctau003","RPV_2017_smu450_neu300_ctau003",
"RPV_2017_smu450_neu350_ctau003","RPV_2017_smu450_neu400_ctau003","RPV_2017_smu450_neu430_ctau003","RPV_2017_smu500_neu180_ctau003","RPV_2017_smu500_neu200_ctau003",
"RPV_2017_smu500_neu250_ctau003","RPV_2017_smu500_neu300_ctau003","RPV_2017_smu500_neu350_ctau003",/*"RPV_2017_smu500_neu400_ctau003",*/"RPV_2017_smu500_neu450_ctau003",
"RPV_2017_smu500_neu480_ctau003"};
 //issue RPV_2017_smu500_neu400_ctau003
 for (int i = 0 ; i< 32 ; i++) 
    {
      TTree* tree =0;
      TreeAnalyzer * tree_ = new TreeAnalyzer(tree, SignalSet003[i],Production );
      tree_->Loop(SignalSet003[i],Production );
      delete tree_;
      delete tree;
    }

      TString SignalSet010[32]={"RPV_2017_smu200_neu180_ctau010","RPV_2017_smu250_neu180_ctau010","RPV_2017_smu250_neu200_ctau010",
"RPV_2017_smu250_neu230_ctau010","RPV_2017_smu300_neu180_ctau010","RPV_2017_smu300_neu200_ctau010","RPV_2017_smu300_neu280_ctau010",
"RPV_2017_smu350_neu180_ctau010","RPV_2017_smu350_neu200_ctau010","RPV_2017_smu350_neu250_ctau010","RPV_2017_smu350_neu300_ctau010","RPV_2017_smu350_neu330_ctau010",
"RPV_2017_smu400_neu180_ctau010","RPV_2017_smu400_neu200_ctau010","RPV_2017_smu400_neu250_ctau010","RPV_2017_smu400_neu300_ctau010","RPV_2017_smu400_neu350_ctau010",
"RPV_2017_smu400_neu380_ctau010","RPV_2017_smu450_neu180_ctau010","RPV_2017_smu450_neu200_ctau010","RPV_2017_smu450_neu250_ctau010","RPV_2017_smu450_neu300_ctau010",
"RPV_2017_smu450_neu350_ctau010","RPV_2017_smu450_neu400_ctau010","RPV_2017_smu450_neu430_ctau010","RPV_2017_smu500_neu180_ctau010","RPV_2017_smu500_neu200_ctau010",
"RPV_2017_smu500_neu250_ctau010","RPV_2017_smu500_neu300_ctau010",/*"RPV_2017_smu500_neu350_ctau010","RPV_2017_smu500_neu400_ctau010",*/"RPV_2017_smu500_neu450_ctau010",
"RPV_2017_smu500_neu480_ctau010"};
 //issue with 500 350 10 // 500 400 10
 for (int i = 0 ; i< 31 ; i++) 
    {
      TTree* tree =0;
      TreeAnalyzer * tree_ = new TreeAnalyzer(tree, SignalSet010[i],Production );
      tree_->Loop(SignalSet010[i],Production );
      delete tree_;
      delete tree;
    }


      TString SignalSet030[34]={"RPV_2017_smu200_neu180_ctau030","RPV_2017_smu250_neu180_ctau030","RPV_2017_smu250_neu200_ctau030",
"RPV_2017_smu250_neu230_ctau030","RPV_2017_smu300_neu180_ctau030","RPV_2017_smu300_neu200_ctau030","RPV_2017_smu300_neu280_ctau030",
"RPV_2017_smu350_neu180_ctau030","RPV_2017_smu350_neu200_ctau030","RPV_2017_smu350_neu250_ctau030","RPV_2017_smu350_neu300_ctau030","RPV_2017_smu350_neu330_ctau030",
"RPV_2017_smu400_neu180_ctau030","RPV_2017_smu400_neu200_ctau030","RPV_2017_smu400_neu250_ctau030","RPV_2017_smu400_neu300_ctau030","RPV_2017_smu400_neu350_ctau030",
"RPV_2017_smu400_neu380_ctau030","RPV_2017_smu450_neu180_ctau030","RPV_2017_smu450_neu200_ctau030","RPV_2017_smu450_neu250_ctau030","RPV_2017_smu450_neu300_ctau030",
"RPV_2017_smu450_neu350_ctau030","RPV_2017_smu450_neu400_ctau030","RPV_2017_smu450_neu430_ctau030","RPV_2017_smu500_neu180_ctau030","RPV_2017_smu500_neu200_ctau030",
"RPV_2017_smu500_neu250_ctau030","RPV_2017_smu500_neu300_ctau030","RPV_2017_smu500_neu350_ctau030","RPV_2017_smu500_neu400_ctau030","RPV_2017_smu500_neu450_ctau030",
"RPV_2017_smu500_neu480_ctau030"};
 
 for (int i = 0 ; i< 33 ; i++) 
    {
      TTree* tree =0;
      TreeAnalyzer * tree_ = new TreeAnalyzer(tree, SignalSet030[i],Production );
      tree_->Loop(SignalSet030[i],Production );
      delete tree_;
      delete tree;
    }

  TString SignalSet100[32]={"RPV_2017_smu200_neu180_ctau100","RPV_2017_smu250_neu180_ctau100","RPV_2017_smu250_neu200_ctau100",
"RPV_2017_smu250_neu230_ctau100","RPV_2017_smu300_neu180_ctau100",/*"RPV_2017_smu300_neu200_ctau100","RPV_2017_smu300_neu250_ctau100_v3",*/"RPV_2017_smu300_neu280_ctau100",
"RPV_2017_smu350_neu180_ctau100","RPV_2017_smu350_neu200_ctau100","RPV_2017_smu350_neu250_ctau100","RPV_2017_smu350_neu300_ctau100","RPV_2017_smu350_neu330_ctau100",
"RPV_2017_smu400_neu180_ctau100","RPV_2017_smu400_neu200_ctau100","RPV_2017_smu400_neu250_ctau100","RPV_2017_smu400_neu300_ctau100","RPV_2017_smu400_neu350_ctau100",
"RPV_2017_smu400_neu380_ctau100","RPV_2017_smu450_neu180_ctau100","RPV_2017_smu450_neu200_ctau100","RPV_2017_smu450_neu250_ctau100","RPV_2017_smu450_neu300_ctau100",
"RPV_2017_smu450_neu350_ctau100","RPV_2017_smu450_neu400_ctau100","RPV_2017_smu450_neu430_ctau100","RPV_2017_smu500_neu180_ctau100","RPV_2017_smu500_neu200_ctau100",
"RPV_2017_smu500_neu250_ctau100","RPV_2017_smu500_neu300_ctau100","RPV_2017_smu500_neu350_ctau100","RPV_2017_smu500_neu400_ctau100","RPV_2017_smu500_neu450_ctau100",
"RPV_2017_smu500_neu480_ctau100"};
 
 for (int i = 0 ; i< 31 ; i++) 
    {
      TTree* tree =0;
      TreeAnalyzer * tree_ = new TreeAnalyzer(tree, SignalSet100[i],Production);
      tree_->Loop(SignalSet100[i],Production );
      delete tree_;
      delete tree;
    }

      TString SignalSet300[34]={"RPV_2017_smu200_neu180_ctau300","RPV_2017_smu250_neu180_ctau300","RPV_2017_smu250_neu200_ctau300",
"RPV_2017_smu250_neu230_ctau300","RPV_2017_smu300_neu180_ctau300","RPV_2017_smu300_neu200_ctau300","RPV_2017_smu300_neu280_ctau300",
"RPV_2017_smu350_neu180_ctau300","RPV_2017_smu350_neu200_ctau300","RPV_2017_smu350_neu250_ctau300","RPV_2017_smu350_neu300_ctau300","RPV_2017_smu350_neu330_ctau300",
"RPV_2017_smu400_neu180_ctau300","RPV_2017_smu400_neu200_ctau300","RPV_2017_smu400_neu250_ctau300","RPV_2017_smu400_neu300_ctau300","RPV_2017_smu400_neu350_ctau300",
"RPV_2017_smu400_neu380_ctau300","RPV_2017_smu450_neu180_ctau300","RPV_2017_smu450_neu200_ctau300","RPV_2017_smu450_neu250_ctau300","RPV_2017_smu450_neu300_ctau300",
"RPV_2017_smu450_neu350_ctau300","RPV_2017_smu450_neu400_ctau300","RPV_2017_smu450_neu430_ctau300","RPV_2017_smu500_neu180_ctau300","RPV_2017_smu500_neu200_ctau300",
"RPV_2017_smu500_neu250_ctau300","RPV_2017_smu500_neu300_ctau300","RPV_2017_smu500_neu350_ctau300","RPV_2017_smu500_neu400_ctau300","RPV_2017_smu500_neu450_ctau300",
"RPV_2017_smu500_neu480_ctau300"};
 
 for (int i = 0 ; i< 33 ; i++) 
    {
      TTree* tree =0;
      TreeAnalyzer * tree_ = new TreeAnalyzer(tree, SignalSet300[i],Production );
      tree_->Loop(SignalSet300[i],Production );
      delete tree_;
      delete tree;
    }


      TString SignalSet1000[33]={"RPV_2017_smu200_neu180_ctau1000","RPV_2017_smu250_neu180_ctau1000","RPV_2017_smu250_neu200_ctau1000",
/*"RPV_2017_smu250_neu230_ctau1000",*/"RPV_2017_smu300_neu180_ctau1000","RPV_2017_smu300_neu200_ctau1000","RPV_2017_smu300_neu280_ctau1000",
"RPV_2017_smu350_neu180_ctau1000","RPV_2017_smu350_neu200_ctau1000","RPV_2017_smu350_neu250_ctau1000","RPV_2017_smu350_neu300_ctau1000","RPV_2017_smu350_neu330_ctau1000",
"RPV_2017_smu400_neu180_ctau1000","RPV_2017_smu400_neu200_ctau1000","RPV_2017_smu400_neu250_ctau1000","RPV_2017_smu400_neu300_ctau1000","RPV_2017_smu400_neu350_ctau1000",
"RPV_2017_smu400_neu380_ctau1000","RPV_2017_smu450_neu180_ctau1000","RPV_2017_smu450_neu200_ctau1000","RPV_2017_smu450_neu250_ctau1000","RPV_2017_smu450_neu300_ctau1000",
"RPV_2017_smu450_neu350_ctau1000","RPV_2017_smu450_neu400_ctau1000","RPV_2017_smu450_neu430_ctau1000","RPV_2017_smu500_neu180_ctau1000","RPV_2017_smu500_neu200_ctau1000",
"RPV_2017_smu500_neu250_ctau1000","RPV_2017_smu500_neu300_ctau1000","RPV_2017_smu500_neu350_ctau1000","RPV_2017_smu500_neu400_ctau1000","RPV_2017_smu500_neu450_ctau1000",
"RPV_2017_smu500_neu480_ctau1000"};
 
 for (int i = 0 ; i< 32 ; i++) 
    {
      TTree* tree =0;
      TreeAnalyzer * tree_ = new TreeAnalyzer(tree, SignalSet1000[i],Production );
      tree_->Loop(SignalSet1000[i],Production );
      delete tree_;
      delete tree;
    }


//----------------------------//

  //   TTree* tree300_250_240=0;
  // TreeAnalyzer * tree_300_250_240 = new TreeAnalyzer(tree300_250_240, "300_250_24", );
  // tree_300_250_240->Loop("smu300_250_24" );
  // delete tree_300_250_240;


  // TTree* tree70=0;
  // TreeAnalyzer * tree_70_test = new TreeAnalyzer(tree70, "70_test", );
  // tree_70_test->Loop("70_test" );
  // delete tree_70_test;

  // TTree* tree50=0;
  // TreeAnalyzer * tree_50_test = new TreeAnalyzer(tree50, "50_test", );
  // tree_50_test->Loop("50_test" );
  // delete tree_50_test;

  // TTree* tree30=0;
  // TreeAnalyzer * tree_30_test = new TreeAnalyzer(tree30, "30_test", );
  // tree_30_test->Loop("30_test" );
  // delete tree_30_test;

  // TTree* tree10=0;
  // TreeAnalyzer * tree_10_test = new TreeAnalyzer(tree10, "10_test", );
  // tree_10_test->Loop("10_test" );
  // delete tree_10_test;



  //   TTree* treesmu200_neu180_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu200_neu180_lamE_2_ctau20 = new TreeAnalyzer(treesmu200_neu180_lamE_2_ctau20, "smu200_neu180_lamE_2_ctau20", );
  // tree_smu200_neu180_lamE_2_ctau20->Loop("smu200_neu180_lamE_2_ctau20" );
  // delete tree_smu200_neu180_lamE_2_ctau20;

  //   TTree* treesmu250_neu200_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu250_neu200_lamE_2_ctau20 = new TreeAnalyzer(treesmu250_neu200_lamE_2_ctau20, "smu250_neu200_lamE_2_ctau20", );
  // tree_smu250_neu200amE_2_ctau20->Loop("smu250_neu200_lamE_2_ctau20" );
  // delete tree_smu250_neu200_lamE_2_ctau20;

  //   TTree* treesmu300_neu180_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu300_neu180_lamE_2_ctau20 = new TreeAnalyzer(treesmu300_neu180_lamE_2_ctau20, "smu300_neu180_lamE_2_ctau20", );
  // tree_smu300_neu180_lamE_2_ctau20->Loop("smu300_neu180_lamE_2_ctau20" );
  // delete tree_smu300_neu180_lamE_2_ctau20;

  //   TTree* treesmu300_neu200_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu300_neu200_lamE_2_ctau20 = new TreeAnalyzer(treesmu300_neu200_lamE_2_ctau20, "smu300_neu200_lamE_2_ctau20", );
  // tree_smu300_neu200_lamE_2_ctau20->Loop("smu300_neu200_lamE_2_ctau20" );
  // delete tree_smu300_neu200_lamE_2_ctau20;

  //   TTree* treesmu300_neu250_lamE_2_ctau01=0;
  // TreeAnalyzer * tree_smu300_neu250_lamE_2_ctau01 = new TreeAnalyzer(treesmu300_neu250_lamE_2_ctau01, "smu300_neu250_lamE_2_ctau01", );
  // tree_smu300_neu250_lamE_2_ctau01->Loop("smu300_neu250_lamE_2_ctau01" );
  // delete tree_smu300_neu250_lamE_2_ctau01;

  //   TTree* treesmu300_neu250_lamE_2_ctau03=0;
  // TreeAnalyzer * tree_smu300_neu250_lamE_2_ctau03 = new TreeAnalyzer(treesmu300_neu250_lamE_2_ctau03, "smu300_neu250_lamE_2_ctau03", );
  // tree_smu300_neu250_lamE_2_ctau03->Loop("smu300_neu250_lamE_2_ctau03" );
  // delete tree_smu300_neu250_lamE_2_ctau03;



  //   TTree* treeNtuple_smu300_neu250_lamE_2_ctau10=0;
  // TreeAnalyzer * tree_Ntuple_smu300_neu250_lamE_2_ctau10 = new TreeAnalyzer(treeNtuple_smu300_neu250_lamE_2_ctau10, "smu300_neu250_lamE_2_ctau10", );
  // tree_Ntuple_smu300_neu250_lamE_2_ctau10->Loop("smu300_neu250_lamE_2_ctau10" );
  // delete tree_Ntuple_smu300_neu250_lamE_2_ctau10;

  // TTree* treesmu300_neu250_lamE_2_ctau100=0;
  // TreeAnalyzer * tree_smu300_neu250_lamE_2_ctau100 = new TreeAnalyzer(treesmu300_neu250_lamE_2_ctau100, "smu300_neu250_lamE_2_ctau100", );
  // tree_smu300_neu250_lamE_2_ctau100->Loop("smu300_neu250_lamE_2_ctau100" );
  // delete tree_smu300_neu250_lamE_2_ctau100;

  //   TTree* treesmu300_neu250_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu300_neu250_lamE_2_ctau20 = new TreeAnalyzer(treesmu300_neu250_lamE_2_ctau20, "smu300_neu250_lamE_2_ctau20", );
  // tree_smu300_neu250_lamE_2_ctau20->Loop("smu300_neu250_lamE_2_ctau20" );
  // delete tree_smu300_neu250_lamE_2_ctau20;

  //   TTree* treesmu300_neu250_lamE_2_ctau30=0;
  // TreeAnalyzer * tree_smu300_neu250_lamE_2_ctau30 = new TreeAnalyzer(treesmu300_neu250_lamE_2_ctau30, "smu300_neu250_lamE_2_ctau30", );
  // tree_smu300_neu250_lamE_2_ctau30->Loop("smu300_neu250_lamE_2_ctau30" );
  // delete tree_smu300_neu250_lamE_2_ctau30;

  //   TTree* treesmu300_neu250_lamE_2_ctaup1=0;
  // TreeAnalyzer * tree_smu300_neu250_lamE_2_ctaup1 = new TreeAnalyzer(treesmu300_neu250_lamE_2_ctaup1, "smu300_neu250_lamE_2_ctaup1", );
  // tree_smu300_neu250_lamE_2_ctaup1->Loop("smu300_neu250_lamE_2_ctaup1" );
  // delete tree_smu300_neu250_lamE_2_ctaup1;

  //   TTree* treesmu300_neu250_lamE_2_ctaup3=0;
  // TreeAnalyzer * tree_smu300_neu250_lamE_2_ctaup3 = new TreeAnalyzer(treesmu300_neu250_lamE_2_ctaup3, "smu300_neu250_lamE_2_ctaup3", );
  // tree_smu300_neu250_lamE_2_ctaup3->Loop("smu300_neu250_lamE_2_ctaup3" );
  // delete tree_smu300_neu250_lamE_2_ctaup3;

  //   TTree* treesmu300_neu280_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu300_neu280_lamE_2_ctau20 = new TreeAnalyzer(treesmu300_neu280_lamE_2_ctau20, "smu300_neu280_lamE_2_ctau20", );
  // tree_smu300_neu280_lamE_2_ctau20->Loop("smu300_neu280_lamE_2_ctau20" );
  // delete tree_smu300_neu280_lamE_2_ctau20;


  //   TTree* treesmu400_neu180_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu400_neu180_lamE_2_ctau20 = new TreeAnalyzer(treesmu400_neu180_lamE_2_ctau20, "smu400_neu180_lamE_2_ctau20", );
  // tree_smu400_neu180_lamE_2_ctau20->Loop("smu400_neu180_lamE_2_ctau20" );
  // delete tree_smu400_neu180_lamE_2_ctau20;

  //   TTree* treesmu400_neu200_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu400_neu200_lamE_2_ctau20 = new TreeAnalyzer(treesmu400_neu200_lamE_2_ctau20, "smu400_neu200_lamE_2_ctau20", );
  // tree_smu400_neu200_lamE_2_ctau20->Loop("smu400_neu200_lamE_2_ctau20" );
  // delete tree_smu400_neu200_lamE_2_ctau20;


  //   TTree* treesmu400_neu250_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu400_neu250_lamE_2_ctau20 = new TreeAnalyzer(treesmu400_neu250_lamE_2_ctau20, "smu400_neu250_lamE_2_ctau20", );
  // tree_smu400_neu250_lamE_2_ctau20->Loop("smu400_neu250_lamE_2_ctau20" );
  // delete tree_smu400_neu250_lamE_2_ctau20;

  //   TTree* treesmu400_neu300_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu400_neu300_lamE_2_ctau20 = new TreeAnalyzer(treesmu400_neu300_lamE_2_ctau20, "smu400_neu300_lamE_2_ctau20", );
  // tree_smu400_neu300_lamE_2_ctau20->Loop("smu400_neu300_lamE_2_ctau20" );
  // delete tree_smu400_neu300_lamE_2_ctau20;

  //   TTree* treesmu400_neu350_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu400_neu350_lamE_2_ctau20 = new TreeAnalyzer(treesmu400_neu350_lamE_2_ctau20, "smu400_neu350_lamE_2_ctau20", );
  // tree_smu400_neu350_lamE_2_ctau20->Loop("smu400_neu350_lamE_2_ctau20" );
  // delete tree_smu400_neu350_lamE_2_ctau20;


  //   TTree* treesmu400_neu380_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu400_neu380_lamE_2_ctau20 = new TreeAnalyzer(treesmu400_neu380_lamE_2_ctau20, "smu400_neu380_lamE_2_ctau20", );
  // tree_smu400_neu380_lamE_2_ctau20->Loop("smu400_neu380_lamE_2_ctau20" );
  // delete tree_smu400_neu380_lamE_2_ctau20;

  //   TTree* treesmu500_neu180_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu500_neu180_lamE_2_ctau20 = new TreeAnalyzer(treesmu500_neu180_lamE_2_ctau20, "smu500_neu180_lamE_2_ctau20", );
  // tree_smu500_neu180_lamE_2_ctau20->Loop("smu500_neu180_lamE_2_ctau20" );
  // delete tree_smu500_neu180_lamE_2_ctau20;

  //   TTree* treesmu500_neu200_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu500_neu200_lamE_2_ctau20 = new TreeAnalyzer(treesmu500_neu200_lamE_2_ctau20, "smu500_neu200_lamE_2_ctau20", );
  // tree_smu500_neu200_lamE_2_ctau20->Loop("smu500_neu200_lamE_2_ctau20" );
  // delete tree_smu500_neu200_lamE_2_ctau20;

  //   TTree* treesmu500_neu250_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu500_neu250_lamE_2_ctau20 = new TreeAnalyzer(treesmu500_neu250_lamE_2_ctau20, "smu500_neu250_lamE_2_ctau20", );
  // tree_smu500_neu250_lamE_2_ctau20->Loop("smu500_neu250_lamE_2_ctau20" );
  // delete tree_smu500_neu250_lamE_2_ctau20;

  //   TTree* treesmu500_neu300_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu500_neu300_lamE_2_ctau20 = new TreeAnalyzer(treesmu500_neu300_lamE_2_ctau20, "smu500_neu300_lamE_2_ctau20", );
  // tree_smu500_neu300_lamE_2_ctau20->Loop("smu500_neu300_lamE_2_ctau20" );
  // delete tree_smu500_neu300_lamE_2_ctau20;

  //   TTree* treesmu500_neu350_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu500_neu350_lamE_2_ctau20 = new TreeAnalyzer(treesmu500_neu350_lamE_2_ctau20, "smu500_neu350_lamE_2_ctau20", );
  // tree_smu500_neu350_lamE_2_ctau20->Loop("smu500_neu350_lamE_2_ctau20" );
  // delete tree_smu500_neu350_lamE_2_ctau20;

  //   TTree* treesmu500_neu400_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu500_neu400_lamE_2_ctau20 = new TreeAnalyzer(treesmu500_neu400_lamE_2_ctau20, "smu500_neu400_lamE_2_ctau20", );
  // tree_smu500_neu400_lamE_2_ctau20->Loop("smu500_neu400_lamE_2_ctau20" );
  // delete tree_smu500_neu400_lamE_2_ctau20;

  //   TTree* treesmu500_neu450_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu500_neu450_lamE_2_ctau20 = new TreeAnalyzer(treesmu500_neu450_lamE_2_ctau20, "smu500_neu450_lamE_2_ctau20", );
  // tree_smu500_neu450_lamE_2_ctau20->Loop("smu500_neu450_lamE_2_ctau20" );
  // delete tree_smu500_neu450_lamE_2_ctau20;

  //   TTree* treesmu500_neu480_lamE_2_ctau20=0;
  // TreeAnalyzer * tree_smu500_neu480_lamE_2_ctau20 = new TreeAnalyzer(treesmu500_neu480_lamE_2_ctau20, "smu500_neu480_lamE_2_ctau20", );
  // tree_smu500_neu480_lamE_2_ctau20->Loop("smu500_neu480_lamE_2_ctau20" );
  // delete tree_smu500_neu480_lamE_2_ctau20;
  
  //   TTree* treeBKG=0;
  // TreeAnalyzer * tree_BKG = new TreeAnalyzer(treeBKG, "BKG", );
  // tree_BKG->Loop("BKG", false);
  // delete tree_BKG;

  // TTree* treeWZTo2Q2L_mllmin4p0=0;
  // TreeAnalyzer * tree_WZTo2Q2L_mllmin4p0 = new TreeAnalyzer(treeWZTo2Q2L_mllmin4p0, "WZTo2Q2L_mllmin4p0", );
  // tree_WZTo2Q2L_mllmin4p0->Loop("WZTo2Q2L_mllmin4p0", false);
  // delete tree_WZTo2Q2L_mllmin4p0;

  // TTree* treeWWTo2L2Nu=0;
  // TreeAnalyzer * tree_WWTo2L2Nu = new TreeAnalyzer(treeWWTo2L2Nu, "WWTo2L2Nu", );
  // tree_WWTo2L2Nu->Loop("WWTo2L2Nu", false);
  // delete tree_WWTo2L2Nu;

  // TTree* treeWWTo2L2Nu_MLL_600To1200=0;
  // TreeAnalyzer * tree_WWTo2L2Nu_MLL_600To1200 = new TreeAnalyzer(treeWWTo2L2Nu_MLL_600To1200, "WWTo2L2Nu_MLL_600To1200", );
  // tree_WWTo2L2Nu_MLL_600To1200->Loop("WWTo2L2Nu_MLL_600To1200", false);
  // delete tree_WWTo2L2Nu_MLL_600To1200;

  // TTree* treeWWTo2L2Nu_MLL_200To600=0;
  // TreeAnalyzer * tree_WWTo2L2Nu_MLL_200To600 = new TreeAnalyzer(treeWWTo2L2Nu_MLL_200To600, "WWTo2L2Nu_MLL_200To600", );
  // tree_WWTo2L2Nu_MLL_200To600->Loop("WWTo2L2Nu_MLL_200To600", false);
  // delete tree_WWTo2L2Nu_MLL_200To600;

  // TTree* treeTTJets_DiLept=0;
  // TreeAnalyzer * tree_TTJets_DiLept = new TreeAnalyzer(treeTTJets_DiLept, "TTJets_DiLept", );
  // tree_TTJets_DiLept->Loop("TTJets_DiLept", false);
  // delete tree_TTJets_DiLept;

  // TTree* treeZZTo2Q2L_mllmin4p0=0;
  // TreeAnalyzer * tree_ZZTo2Q2L_mllmin4p0 = new TreeAnalyzer(treeZZTo2Q2L_mllmin4p0, "ZZTo2Q2L_mllmin4p0", );
  // tree_ZZTo2Q2L_mllmin4p0->Loop("ZZTo2Q2L_mllmin4p0", false);
  // delete tree_ZZTo2Q2L_mllmin4p0;

  // TTree* treeST_tW_top_5f_NoFullyHadronicDecays=0;
  // TreeAnalyzer * tree_ST_tW_top_5f_NoFullyHadronicDecays = new TreeAnalyzer(treeST_tW_top_5f_NoFullyHadronicDecays, "ST_tW_top_5f_NoFullyHadronicDecays", );
  // tree_ST_tW_top_5f_NoFullyHadronicDecays->Loop("ST_tW_top_5f_NoFullyHadronicDecays", false);
  // delete tree_ST_tW_top_5f_NoFullyHadronicDecays;


  //   TTree* treeDYJetsToLL_M10to50=0;
  // TreeAnalyzer * tree_DYJetsToLL_M10to50 = new TreeAnalyzer(treeDYJetsToLL_M10to50, "DYJetsToLL_M10to50", );
  // tree_DYJetsToLL_M10to50->Loop("DYJetsToLL_M10to50", false);
  // delete tree_DYJetsToLL_M10to50;

  // TTree* treeDYJetsToLL_M50=0;
  // TreeAnalyzer * tree_DYJetsToLL_M50 = new TreeAnalyzer(treeDYJetsToLL_M50, "DYJetsToLL_M50", );
  // tree_DYJetsToLL_M50->Loop("DYJetsToLL_M50", false);
  // delete tree_DYJetsToLL_M50;

  // TTree* treeST_tW_antitop_5f_NoFullyHadronicDecays=0;
  // TreeAnalyzer * tree_ST_tW_antitop_5f_NoFullyHadronicDecays = new TreeAnalyzer(treeST_tW_antitop_5f_NoFullyHadronicDecays, "ST_tW_antitop_5f_NoFullyHadronicDecays", );
  // tree_ST_tW_antitop_5f_NoFullyHadronicDecays->Loop("ST_tW_antitop_5f_NoFullyHadronicDecays", false);
  // delete tree_ST_tW_antitop_5f_NoFullyHadronicDecays;

  // TTree* treeTTTo2L2Nu=0;
  // TreeAnalyzer * tree_TTTo2L2Nu = new TreeAnalyzer(treeTTTo2L2Nu, "TTTo2L2Nu", );
  // tree_TTTo2L2Nu->Loop("TTTo2L2Nu", false);
  // delete tree_TTTo2L2Nu;

}
