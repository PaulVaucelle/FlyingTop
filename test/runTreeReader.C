{
//  gROOT->ProcessLine(".L TreeReader.C+");
  std::vector<TString > systlist;
  systlist.push_back("");

//   TString BKGSet[16]={"TTJets_DiLept","TTTo2L2Nu","ST_tW_top_5f_NoFullyHadronicDecays","ST_tW_antitop_5f_NoFullyHadronicDecays",
//                       "DYJetsToLL_M10to50","DYJetsToLL_M50","WWTo2L2Nu","WWTo2L2Nu_MLL_200To600","WWTo2L2Nu_MLL_600To1200",
//                       "WZTo2Q2L_mllmin4p0","ZZTo2Q2L_mllmin4p0","ttWJetsToLNu_5f_EWK","TTZToLL_5f","TTToHadronic","TTWW","TTToSemiLeptonic"
//   };
 
//  for (int i = 0 ; i< 15 ; i++) 
//     {
//       TTree* treeB =0;
//       TreeReader * treeB_ = new TreeReader(treeB, BKGSet[i], systlist);
//       treeB_->Loop(BKGSet[i],false);
//       delete treeB_;
//       delete treeB;
//     }


//       TString SignalSet001[32]={"RPV_2018_smu200_neu180_ctau001","RPV_2018_smu250_neu180_ctau001","RPV_2018_smu250_neu200_ctau001",
// "RPV_2018_smu250_neu230_ctau001","RPV_2018_smu300_neu180_ctau001","RPV_2018_smu300_neu200_ctau001","RPV_2018_smu300_neu280_ctau001",
// "RPV_2018_smu350_neu180_ctau001","RPV_2018_smu350_neu200_ctau001","RPV_2018_smu350_neu250_ctau001","RPV_2018_smu350_neu300_ctau001","RPV_2018_smu350_neu330_ctau001",
// "RPV_2018_smu400_neu180_ctau001","RPV_2018_smu400_neu200_ctau001","RPV_2018_smu400_neu250_ctau001","RPV_2018_smu400_neu300_ctau001","RPV_2018_smu400_neu350_ctau001",
// "RPV_2018_smu400_neu380_ctau001","RPV_2018_smu450_neu180_ctau001","RPV_2018_smu450_neu200_ctau001","RPV_2018_smu450_neu250_ctau001","RPV_2018_smu450_neu300_ctau001",
// "RPV_2018_smu450_neu350_ctau001","RPV_2018_smu450_neu400_ctau001","RPV_2018_smu450_neu430_ctau001","RPV_2018_smu500_neu180_ctau001","RPV_2018_smu500_neu200_ctau001",
// "RPV_2018_smu500_neu250_ctau001","RPV_2018_smu500_neu300_ctau001","RPV_2018_smu500_neu350_ctau001",/*"RPV_2018_smu500_neu400_ctau001","RPV_2018_smu500_neu450_ctau001",*/
// "RPV_2018_smu500_neu480_ctau001"};
//  //issue  with RPV_2018_smu500_neu400_ctau001 500 450 001
//  for (int i = 0 ; i< 31 ; i++) 
//     {
//       TTree* tree =0;
//       TreeReader * tree_ = new TreeReader(tree, SignalSet001[i], systlist);
//       tree_->Loop(SignalSet001[i],true);
//       delete tree_;
//       delete tree;
//     }



//       TString SignalSet003[33]={"RPV_2018_smu200_neu180_ctau003","RPV_2018_smu250_neu180_ctau003","RPV_2018_smu250_neu200_ctau003",
// "RPV_2018_smu250_neu230_ctau003","RPV_2018_smu300_neu180_ctau003","RPV_2018_smu300_neu200_ctau003","RPV_2018_smu300_neu280_ctau003",
// "RPV_2018_smu350_neu180_ctau003","RPV_2018_smu350_neu200_ctau003","RPV_2018_smu350_neu250_ctau003","RPV_2018_smu350_neu300_ctau003","RPV_2018_smu350_neu330_ctau003",
// "RPV_2018_smu400_neu180_ctau003","RPV_2018_smu400_neu200_ctau003","RPV_2018_smu400_neu250_ctau003","RPV_2018_smu400_neu300_ctau003","RPV_2018_smu400_neu350_ctau003",
// "RPV_2018_smu400_neu380_ctau003","RPV_2018_smu450_neu180_ctau003","RPV_2018_smu450_neu200_ctau003","RPV_2018_smu450_neu250_ctau003","RPV_2018_smu450_neu300_ctau003",
// "RPV_2018_smu450_neu350_ctau003","RPV_2018_smu450_neu400_ctau003","RPV_2018_smu450_neu430_ctau003","RPV_2018_smu500_neu180_ctau003","RPV_2018_smu500_neu200_ctau003",
// "RPV_2018_smu500_neu250_ctau003","RPV_2018_smu500_neu300_ctau003","RPV_2018_smu500_neu350_ctau003",/*"RPV_2018_smu500_neu400_ctau003",*/"RPV_2018_smu500_neu450_ctau003",
// "RPV_2018_smu500_neu480_ctau003"};
//  //issue RPV_2018_smu500_neu400_ctau003
//  for (int i = 0 ; i< 32 ; i++) 
//     {
//       TTree* tree =0;
//       TreeReader * tree_ = new TreeReader(tree, SignalSet003[i], systlist);
//       tree_->Loop(SignalSet003[i],true);
//       delete tree_;
//       delete tree;
//     }

//       TString SignalSet010[32]={"RPV_2018_smu200_neu180_ctau010","RPV_2018_smu250_neu180_ctau010","RPV_2018_smu250_neu200_ctau010",
// "RPV_2018_smu250_neu230_ctau010","RPV_2018_smu300_neu180_ctau010","RPV_2018_smu300_neu200_ctau010","RPV_2018_smu300_neu280_ctau010",
// "RPV_2018_smu350_neu180_ctau010","RPV_2018_smu350_neu200_ctau010","RPV_2018_smu350_neu250_ctau010","RPV_2018_smu350_neu300_ctau010","RPV_2018_smu350_neu330_ctau010",
// "RPV_2018_smu400_neu180_ctau010","RPV_2018_smu400_neu200_ctau010","RPV_2018_smu400_neu250_ctau010","RPV_2018_smu400_neu300_ctau010","RPV_2018_smu400_neu350_ctau010",
// "RPV_2018_smu400_neu380_ctau010","RPV_2018_smu450_neu180_ctau010","RPV_2018_smu450_neu200_ctau010","RPV_2018_smu450_neu250_ctau010","RPV_2018_smu450_neu300_ctau010",
// "RPV_2018_smu450_neu350_ctau010","RPV_2018_smu450_neu400_ctau010","RPV_2018_smu450_neu430_ctau010","RPV_2018_smu500_neu180_ctau010","RPV_2018_smu500_neu200_ctau010",
// "RPV_2018_smu500_neu250_ctau010","RPV_2018_smu500_neu300_ctau010",/*"RPV_2018_smu500_neu350_ctau010","RPV_2018_smu500_neu400_ctau010",*/"RPV_2018_smu500_neu450_ctau010",
// "RPV_2018_smu500_neu480_ctau010"};
//  //issue with 500 350 10 // 500 400 10
//  for (int i = 0 ; i< 31 ; i++) 
//     {
//       TTree* tree =0;
//       TreeReader * tree_ = new TreeReader(tree, SignalSet010[i], systlist);
//       tree_->Loop(SignalSet010[i],true);
//       delete tree_;
//       delete tree;
//     }


//       TString SignalSet030[34]={"RPV_2018_smu200_neu180_ctau030","RPV_2018_smu250_neu180_ctau030","RPV_2018_smu250_neu200_ctau030",
// "RPV_2018_smu250_neu230_ctau030","RPV_2018_smu300_neu180_ctau030","RPV_2018_smu300_neu200_ctau030","RPV_2018_smu300_neu280_ctau030",
// "RPV_2018_smu350_neu180_ctau030","RPV_2018_smu350_neu200_ctau030","RPV_2018_smu350_neu250_ctau030","RPV_2018_smu350_neu300_ctau030","RPV_2018_smu350_neu330_ctau030",
// "RPV_2018_smu400_neu180_ctau030","RPV_2018_smu400_neu200_ctau030","RPV_2018_smu400_neu250_ctau030","RPV_2018_smu400_neu300_ctau030","RPV_2018_smu400_neu350_ctau030",
// "RPV_2018_smu400_neu380_ctau030","RPV_2018_smu450_neu180_ctau030","RPV_2018_smu450_neu200_ctau030","RPV_2018_smu450_neu250_ctau030","RPV_2018_smu450_neu300_ctau030",
// "RPV_2018_smu450_neu350_ctau030","RPV_2018_smu450_neu400_ctau030","RPV_2018_smu450_neu430_ctau030","RPV_2018_smu500_neu180_ctau030","RPV_2018_smu500_neu200_ctau030",
// "RPV_2018_smu500_neu250_ctau030","RPV_2018_smu500_neu300_ctau030","RPV_2018_smu500_neu350_ctau030","RPV_2018_smu500_neu400_ctau030","RPV_2018_smu500_neu450_ctau030",
// "RPV_2018_smu500_neu480_ctau030"};
 
//  for (int i = 0 ; i< 33 ; i++) 
//     {
//       TTree* tree =0;
//       TreeReader * tree_ = new TreeReader(tree, SignalSet030[i], systlist);
//       tree_->Loop(SignalSet030[i],true);
//       delete tree_;
//       delete tree;
//     }

//   TString SignalSet100[32]={"RPV_2018_smu200_neu180_ctau100","RPV_2018_smu250_neu180_ctau100","RPV_2018_smu250_neu200_ctau100",
// "RPV_2018_smu250_neu230_ctau100","RPV_2018_smu300_neu180_ctau100",/*"RPV_2018_smu300_neu200_ctau100","RPV_2018_smu300_neu250_ctau100_v3",*/"RPV_2018_smu300_neu280_ctau100",
// "RPV_2018_smu350_neu180_ctau100","RPV_2018_smu350_neu200_ctau100","RPV_2018_smu350_neu250_ctau100","RPV_2018_smu350_neu300_ctau100","RPV_2018_smu350_neu330_ctau100",
// "RPV_2018_smu400_neu180_ctau100","RPV_2018_smu400_neu200_ctau100","RPV_2018_smu400_neu250_ctau100","RPV_2018_smu400_neu300_ctau100","RPV_2018_smu400_neu350_ctau100",
// "RPV_2018_smu400_neu380_ctau100","RPV_2018_smu450_neu180_ctau100","RPV_2018_smu450_neu200_ctau100","RPV_2018_smu450_neu250_ctau100","RPV_2018_smu450_neu300_ctau100",
// "RPV_2018_smu450_neu350_ctau100","RPV_2018_smu450_neu400_ctau100","RPV_2018_smu450_neu430_ctau100","RPV_2018_smu500_neu180_ctau100","RPV_2018_smu500_neu200_ctau100",
// "RPV_2018_smu500_neu250_ctau100","RPV_2018_smu500_neu300_ctau100","RPV_2018_smu500_neu350_ctau100","RPV_2018_smu500_neu400_ctau100","RPV_2018_smu500_neu450_ctau100",
// "RPV_2018_smu500_neu480_ctau100"};
 
//  for (int i = 0 ; i< 31 ; i++) 
//     {
//       TTree* tree =0;
//       TreeReader * tree_ = new TreeReader(tree, SignalSet100[i], systlist);
//       tree_->Loop(SignalSet100[i],true);
//       delete tree_;
//       delete tree;
//     }

//       TString SignalSet300[34]={"RPV_2018_smu200_neu180_ctau300","RPV_2018_smu250_neu180_ctau300","RPV_2018_smu250_neu200_ctau300",
// "RPV_2018_smu250_neu230_ctau300","RPV_2018_smu300_neu180_ctau300","RPV_2018_smu300_neu200_ctau300","RPV_2018_smu300_neu280_ctau300",
// "RPV_2018_smu350_neu180_ctau300","RPV_2018_smu350_neu200_ctau300","RPV_2018_smu350_neu250_ctau300","RPV_2018_smu350_neu300_ctau300","RPV_2018_smu350_neu330_ctau300",
// "RPV_2018_smu400_neu180_ctau300","RPV_2018_smu400_neu200_ctau300","RPV_2018_smu400_neu250_ctau300","RPV_2018_smu400_neu300_ctau300","RPV_2018_smu400_neu350_ctau300",
// "RPV_2018_smu400_neu380_ctau300","RPV_2018_smu450_neu180_ctau300","RPV_2018_smu450_neu200_ctau300","RPV_2018_smu450_neu250_ctau300","RPV_2018_smu450_neu300_ctau300",
// "RPV_2018_smu450_neu350_ctau300","RPV_2018_smu450_neu400_ctau300","RPV_2018_smu450_neu430_ctau300","RPV_2018_smu500_neu180_ctau300","RPV_2018_smu500_neu200_ctau300",
// "RPV_2018_smu500_neu250_ctau300","RPV_2018_smu500_neu300_ctau300","RPV_2018_smu500_neu350_ctau300","RPV_2018_smu500_neu400_ctau300","RPV_2018_smu500_neu450_ctau300",
// "RPV_2018_smu500_neu480_ctau300"};
 
//  for (int i = 0 ; i< 33 ; i++) 
//     {
//       TTree* tree =0;
//       TreeReader * tree_ = new TreeReader(tree, SignalSet300[i], systlist);
//       tree_->Loop(SignalSet300[i],true);
//       delete tree_;
//       delete tree;
//     }


//       TString SignalSet1000[33]={"RPV_2018_smu200_neu180_ctau1000","RPV_2018_smu250_neu180_ctau1000","RPV_2018_smu250_neu200_ctau1000",
// /*"RPV_2018_smu250_neu230_ctau1000",*/"RPV_2018_smu300_neu180_ctau1000","RPV_2018_smu300_neu200_ctau1000","RPV_2018_smu300_neu280_ctau1000",
// "RPV_2018_smu350_neu180_ctau1000","RPV_2018_smu350_neu200_ctau1000","RPV_2018_smu350_neu250_ctau1000","RPV_2018_smu350_neu300_ctau1000","RPV_2018_smu350_neu330_ctau1000",
// "RPV_2018_smu400_neu180_ctau1000","RPV_2018_smu400_neu200_ctau1000","RPV_2018_smu400_neu250_ctau1000","RPV_2018_smu400_neu300_ctau1000","RPV_2018_smu400_neu350_ctau1000",
// "RPV_2018_smu400_neu380_ctau1000","RPV_2018_smu450_neu180_ctau1000","RPV_2018_smu450_neu200_ctau1000","RPV_2018_smu450_neu250_ctau1000","RPV_2018_smu450_neu300_ctau1000",
// "RPV_2018_smu450_neu350_ctau1000","RPV_2018_smu450_neu400_ctau1000","RPV_2018_smu450_neu430_ctau1000","RPV_2018_smu500_neu180_ctau1000","RPV_2018_smu500_neu200_ctau1000",
// "RPV_2018_smu500_neu250_ctau1000","RPV_2018_smu500_neu300_ctau1000","RPV_2018_smu500_neu350_ctau1000","RPV_2018_smu500_neu400_ctau1000","RPV_2018_smu500_neu450_ctau1000",
// "RPV_2018_smu500_neu480_ctau1000"};
 
//  for (int i = 0 ; i< 32 ; i++) 
//     {
//       TTree* tree =0;
//       TreeReader * tree_ = new TreeReader(tree, SignalSet1000[i], systlist);
//       tree_->Loop(SignalSet1000[i],true);
//       delete tree_;
//       delete tree;
//     }

  //   TTree* tree300_250_240=0;
  // TreeReader * tree_300_250_240 = new TreeReader(tree300_250_240, "300_250_24", systlist);
  // tree_300_250_240->Loop("smu300_250_24",true);
  // delete tree_300_250_240;


  // TTree* tree70=0;
  // TreeReader * tree_70_test = new TreeReader(tree70, "70_test", systlist);
  // tree_70_test->Loop("70_test",true);
  // delete tree_70_test;

  TTree* tree50=0;
  TreeReader * tree_50_test = new TreeReader(tree50, "50_test", systlist);
  tree_50_test->Loop("50_test",true);
  delete tree_50_test;

  // TTree* tree30=0;
  // TreeReader * tree_30_test = new TreeReader(tree30, "30_test", systlist);
  // tree_30_test->Loop("30_test",true);
  // delete tree_30_test;

  // TTree* tree10=0;
  // TreeReader * tree_10_test = new TreeReader(tree10, "10_test", systlist);
  // tree_10_test->Loop("10_test",true);
  // delete tree_10_test;



  //   TTree* treesmu200_neu180_lamE_2_ctau20=0;
  // TreeReader * tree_smu200_neu180_lamE_2_ctau20 = new TreeReader(treesmu200_neu180_lamE_2_ctau20, "smu200_neu180_lamE_2_ctau20", systlist);
  // tree_smu200_neu180_lamE_2_ctau20->Loop("smu200_neu180_lamE_2_ctau20",true);
  // delete tree_smu200_neu180_lamE_2_ctau20;

  //   TTree* treesmu250_neu200_lamE_2_ctau20=0;
  // TreeReader * tree_smu250_neu200_lamE_2_ctau20 = new TreeReader(treesmu250_neu200_lamE_2_ctau20, "smu250_neu200_lamE_2_ctau20", systlist);
  // tree_smu250_neu200amE_2_ctau20->Loop("smu250_neu200_lamE_2_ctau20",true);
  // delete tree_smu250_neu200_lamE_2_ctau20;

  //   TTree* treesmu300_neu180_lamE_2_ctau20=0;
  // TreeReader * tree_smu300_neu180_lamE_2_ctau20 = new TreeReader(treesmu300_neu180_lamE_2_ctau20, "smu300_neu180_lamE_2_ctau20", systlist);
  // tree_smu300_neu180_lamE_2_ctau20->Loop("smu300_neu180_lamE_2_ctau20",true);
  // delete tree_smu300_neu180_lamE_2_ctau20;

  //   TTree* treesmu300_neu200_lamE_2_ctau20=0;
  // TreeReader * tree_smu300_neu200_lamE_2_ctau20 = new TreeReader(treesmu300_neu200_lamE_2_ctau20, "smu300_neu200_lamE_2_ctau20", systlist);
  // tree_smu300_neu200_lamE_2_ctau20->Loop("smu300_neu200_lamE_2_ctau20",true);
  // delete tree_smu300_neu200_lamE_2_ctau20;

  //   TTree* treesmu300_neu250_lamE_2_ctau01=0;
  // TreeReader * tree_smu300_neu250_lamE_2_ctau01 = new TreeReader(treesmu300_neu250_lamE_2_ctau01, "smu300_neu250_lamE_2_ctau01", systlist);
  // tree_smu300_neu250_lamE_2_ctau01->Loop("smu300_neu250_lamE_2_ctau01",true);
  // delete tree_smu300_neu250_lamE_2_ctau01;

  //   TTree* treesmu300_neu250_lamE_2_ctau03=0;
  // TreeReader * tree_smu300_neu250_lamE_2_ctau03 = new TreeReader(treesmu300_neu250_lamE_2_ctau03, "smu300_neu250_lamE_2_ctau03", systlist);
  // tree_smu300_neu250_lamE_2_ctau03->Loop("smu300_neu250_lamE_2_ctau03",true);
  // delete tree_smu300_neu250_lamE_2_ctau03;



  //   TTree* treeNtuple_smu300_neu250_lamE_2_ctau10=0;
  // TreeReader * tree_Ntuple_smu300_neu250_lamE_2_ctau10 = new TreeReader(treeNtuple_smu300_neu250_lamE_2_ctau10, "smu300_neu250_lamE_2_ctau10", systlist);
  // tree_Ntuple_smu300_neu250_lamE_2_ctau10->Loop("smu300_neu250_lamE_2_ctau10",true);
  // delete tree_Ntuple_smu300_neu250_lamE_2_ctau10;

  // TTree* treesmu300_neu250_lamE_2_ctau100=0;
  // TreeReader * tree_smu300_neu250_lamE_2_ctau100 = new TreeReader(treesmu300_neu250_lamE_2_ctau100, "smu300_neu250_lamE_2_ctau100", systlist);
  // tree_smu300_neu250_lamE_2_ctau100->Loop("smu300_neu250_lamE_2_ctau100",true);
  // delete tree_smu300_neu250_lamE_2_ctau100;

  //   TTree* treesmu300_neu250_lamE_2_ctau20=0;
  // TreeReader * tree_smu300_neu250_lamE_2_ctau20 = new TreeReader(treesmu300_neu250_lamE_2_ctau20, "smu300_neu250_lamE_2_ctau20", systlist);
  // tree_smu300_neu250_lamE_2_ctau20->Loop("smu300_neu250_lamE_2_ctau20",true);
  // delete tree_smu300_neu250_lamE_2_ctau20;

  //   TTree* treesmu300_neu250_lamE_2_ctau30=0;
  // TreeReader * tree_smu300_neu250_lamE_2_ctau30 = new TreeReader(treesmu300_neu250_lamE_2_ctau30, "smu300_neu250_lamE_2_ctau30", systlist);
  // tree_smu300_neu250_lamE_2_ctau30->Loop("smu300_neu250_lamE_2_ctau30",true);
  // delete tree_smu300_neu250_lamE_2_ctau30;

  //   TTree* treesmu300_neu250_lamE_2_ctaup1=0;
  // TreeReader * tree_smu300_neu250_lamE_2_ctaup1 = new TreeReader(treesmu300_neu250_lamE_2_ctaup1, "smu300_neu250_lamE_2_ctaup1", systlist);
  // tree_smu300_neu250_lamE_2_ctaup1->Loop("smu300_neu250_lamE_2_ctaup1",true);
  // delete tree_smu300_neu250_lamE_2_ctaup1;

  //   TTree* treesmu300_neu250_lamE_2_ctaup3=0;
  // TreeReader * tree_smu300_neu250_lamE_2_ctaup3 = new TreeReader(treesmu300_neu250_lamE_2_ctaup3, "smu300_neu250_lamE_2_ctaup3", systlist);
  // tree_smu300_neu250_lamE_2_ctaup3->Loop("smu300_neu250_lamE_2_ctaup3",true);
  // delete tree_smu300_neu250_lamE_2_ctaup3;

  //   TTree* treesmu300_neu280_lamE_2_ctau20=0;
  // TreeReader * tree_smu300_neu280_lamE_2_ctau20 = new TreeReader(treesmu300_neu280_lamE_2_ctau20, "smu300_neu280_lamE_2_ctau20", systlist);
  // tree_smu300_neu280_lamE_2_ctau20->Loop("smu300_neu280_lamE_2_ctau20",true);
  // delete tree_smu300_neu280_lamE_2_ctau20;


  //   TTree* treesmu400_neu180_lamE_2_ctau20=0;
  // TreeReader * tree_smu400_neu180_lamE_2_ctau20 = new TreeReader(treesmu400_neu180_lamE_2_ctau20, "smu400_neu180_lamE_2_ctau20", systlist);
  // tree_smu400_neu180_lamE_2_ctau20->Loop("smu400_neu180_lamE_2_ctau20",true);
  // delete tree_smu400_neu180_lamE_2_ctau20;

  //   TTree* treesmu400_neu200_lamE_2_ctau20=0;
  // TreeReader * tree_smu400_neu200_lamE_2_ctau20 = new TreeReader(treesmu400_neu200_lamE_2_ctau20, "smu400_neu200_lamE_2_ctau20", systlist);
  // tree_smu400_neu200_lamE_2_ctau20->Loop("smu400_neu200_lamE_2_ctau20",true);
  // delete tree_smu400_neu200_lamE_2_ctau20;


  //   TTree* treesmu400_neu250_lamE_2_ctau20=0;
  // TreeReader * tree_smu400_neu250_lamE_2_ctau20 = new TreeReader(treesmu400_neu250_lamE_2_ctau20, "smu400_neu250_lamE_2_ctau20", systlist);
  // tree_smu400_neu250_lamE_2_ctau20->Loop("smu400_neu250_lamE_2_ctau20",true);
  // delete tree_smu400_neu250_lamE_2_ctau20;

  //   TTree* treesmu400_neu300_lamE_2_ctau20=0;
  // TreeReader * tree_smu400_neu300_lamE_2_ctau20 = new TreeReader(treesmu400_neu300_lamE_2_ctau20, "smu400_neu300_lamE_2_ctau20", systlist);
  // tree_smu400_neu300_lamE_2_ctau20->Loop("smu400_neu300_lamE_2_ctau20",true);
  // delete tree_smu400_neu300_lamE_2_ctau20;

  //   TTree* treesmu400_neu350_lamE_2_ctau20=0;
  // TreeReader * tree_smu400_neu350_lamE_2_ctau20 = new TreeReader(treesmu400_neu350_lamE_2_ctau20, "smu400_neu350_lamE_2_ctau20", systlist);
  // tree_smu400_neu350_lamE_2_ctau20->Loop("smu400_neu350_lamE_2_ctau20",true);
  // delete tree_smu400_neu350_lamE_2_ctau20;


  //   TTree* treesmu400_neu380_lamE_2_ctau20=0;
  // TreeReader * tree_smu400_neu380_lamE_2_ctau20 = new TreeReader(treesmu400_neu380_lamE_2_ctau20, "smu400_neu380_lamE_2_ctau20", systlist);
  // tree_smu400_neu380_lamE_2_ctau20->Loop("smu400_neu380_lamE_2_ctau20",true);
  // delete tree_smu400_neu380_lamE_2_ctau20;

  //   TTree* treesmu500_neu180_lamE_2_ctau20=0;
  // TreeReader * tree_smu500_neu180_lamE_2_ctau20 = new TreeReader(treesmu500_neu180_lamE_2_ctau20, "smu500_neu180_lamE_2_ctau20", systlist);
  // tree_smu500_neu180_lamE_2_ctau20->Loop("smu500_neu180_lamE_2_ctau20",true);
  // delete tree_smu500_neu180_lamE_2_ctau20;

  //   TTree* treesmu500_neu200_lamE_2_ctau20=0;
  // TreeReader * tree_smu500_neu200_lamE_2_ctau20 = new TreeReader(treesmu500_neu200_lamE_2_ctau20, "smu500_neu200_lamE_2_ctau20", systlist);
  // tree_smu500_neu200_lamE_2_ctau20->Loop("smu500_neu200_lamE_2_ctau20",true);
  // delete tree_smu500_neu200_lamE_2_ctau20;

  //   TTree* treesmu500_neu250_lamE_2_ctau20=0;
  // TreeReader * tree_smu500_neu250_lamE_2_ctau20 = new TreeReader(treesmu500_neu250_lamE_2_ctau20, "smu500_neu250_lamE_2_ctau20", systlist);
  // tree_smu500_neu250_lamE_2_ctau20->Loop("smu500_neu250_lamE_2_ctau20",true);
  // delete tree_smu500_neu250_lamE_2_ctau20;

  //   TTree* treesmu500_neu300_lamE_2_ctau20=0;
  // TreeReader * tree_smu500_neu300_lamE_2_ctau20 = new TreeReader(treesmu500_neu300_lamE_2_ctau20, "smu500_neu300_lamE_2_ctau20", systlist);
  // tree_smu500_neu300_lamE_2_ctau20->Loop("smu500_neu300_lamE_2_ctau20",true);
  // delete tree_smu500_neu300_lamE_2_ctau20;

  //   TTree* treesmu500_neu350_lamE_2_ctau20=0;
  // TreeReader * tree_smu500_neu350_lamE_2_ctau20 = new TreeReader(treesmu500_neu350_lamE_2_ctau20, "smu500_neu350_lamE_2_ctau20", systlist);
  // tree_smu500_neu350_lamE_2_ctau20->Loop("smu500_neu350_lamE_2_ctau20",true);
  // delete tree_smu500_neu350_lamE_2_ctau20;

  //   TTree* treesmu500_neu400_lamE_2_ctau20=0;
  // TreeReader * tree_smu500_neu400_lamE_2_ctau20 = new TreeReader(treesmu500_neu400_lamE_2_ctau20, "smu500_neu400_lamE_2_ctau20", systlist);
  // tree_smu500_neu400_lamE_2_ctau20->Loop("smu500_neu400_lamE_2_ctau20",true);
  // delete tree_smu500_neu400_lamE_2_ctau20;

  //   TTree* treesmu500_neu450_lamE_2_ctau20=0;
  // TreeReader * tree_smu500_neu450_lamE_2_ctau20 = new TreeReader(treesmu500_neu450_lamE_2_ctau20, "smu500_neu450_lamE_2_ctau20", systlist);
  // tree_smu500_neu450_lamE_2_ctau20->Loop("smu500_neu450_lamE_2_ctau20",true);
  // delete tree_smu500_neu450_lamE_2_ctau20;

  //   TTree* treesmu500_neu480_lamE_2_ctau20=0;
  // TreeReader * tree_smu500_neu480_lamE_2_ctau20 = new TreeReader(treesmu500_neu480_lamE_2_ctau20, "smu500_neu480_lamE_2_ctau20", systlist);
  // tree_smu500_neu480_lamE_2_ctau20->Loop("smu500_neu480_lamE_2_ctau20",true);
  // delete tree_smu500_neu480_lamE_2_ctau20;
  
  //   TTree* treeBKG=0;
  // TreeReader * tree_BKG = new TreeReader(treeBKG, "BKG", systlist);
  // tree_BKG->Loop("BKG", false);
  // delete tree_BKG;

  // TTree* treeWZTo2Q2L_mllmin4p0=0;
  // TreeReader * tree_WZTo2Q2L_mllmin4p0 = new TreeReader(treeWZTo2Q2L_mllmin4p0, "WZTo2Q2L_mllmin4p0", systlist);
  // tree_WZTo2Q2L_mllmin4p0->Loop("WZTo2Q2L_mllmin4p0", false);
  // delete tree_WZTo2Q2L_mllmin4p0;

  // TTree* treeWWTo2L2Nu=0;
  // TreeReader * tree_WWTo2L2Nu = new TreeReader(treeWWTo2L2Nu, "WWTo2L2Nu", systlist);
  // tree_WWTo2L2Nu->Loop("WWTo2L2Nu", false);
  // delete tree_WWTo2L2Nu;

  // TTree* treeWWTo2L2Nu_MLL_600To1200=0;
  // TreeReader * tree_WWTo2L2Nu_MLL_600To1200 = new TreeReader(treeWWTo2L2Nu_MLL_600To1200, "WWTo2L2Nu_MLL_600To1200", systlist);
  // tree_WWTo2L2Nu_MLL_600To1200->Loop("WWTo2L2Nu_MLL_600To1200", false);
  // delete tree_WWTo2L2Nu_MLL_600To1200;

  // TTree* treeWWTo2L2Nu_MLL_200To600=0;
  // TreeReader * tree_WWTo2L2Nu_MLL_200To600 = new TreeReader(treeWWTo2L2Nu_MLL_200To600, "WWTo2L2Nu_MLL_200To600", systlist);
  // tree_WWTo2L2Nu_MLL_200To600->Loop("WWTo2L2Nu_MLL_200To600", false);
  // delete tree_WWTo2L2Nu_MLL_200To600;

  // TTree* treeTTJets_DiLept=0;
  // TreeReader * tree_TTJets_DiLept = new TreeReader(treeTTJets_DiLept, "TTJets_DiLept", systlist);
  // tree_TTJets_DiLept->Loop("TTJets_DiLept", false);
  // delete tree_TTJets_DiLept;

  // TTree* treeZZTo2Q2L_mllmin4p0=0;
  // TreeReader * tree_ZZTo2Q2L_mllmin4p0 = new TreeReader(treeZZTo2Q2L_mllmin4p0, "ZZTo2Q2L_mllmin4p0", systlist);
  // tree_ZZTo2Q2L_mllmin4p0->Loop("ZZTo2Q2L_mllmin4p0", false);
  // delete tree_ZZTo2Q2L_mllmin4p0;

  // TTree* treeST_tW_top_5f_NoFullyHadronicDecays=0;
  // TreeReader * tree_ST_tW_top_5f_NoFullyHadronicDecays = new TreeReader(treeST_tW_top_5f_NoFullyHadronicDecays, "ST_tW_top_5f_NoFullyHadronicDecays", systlist);
  // tree_ST_tW_top_5f_NoFullyHadronicDecays->Loop("ST_tW_top_5f_NoFullyHadronicDecays", false);
  // delete tree_ST_tW_top_5f_NoFullyHadronicDecays;


  //   TTree* treeDYJetsToLL_M10to50=0;
  // TreeReader * tree_DYJetsToLL_M10to50 = new TreeReader(treeDYJetsToLL_M10to50, "DYJetsToLL_M10to50", systlist);
  // tree_DYJetsToLL_M10to50->Loop("DYJetsToLL_M10to50", false);
  // delete tree_DYJetsToLL_M10to50;

  // TTree* treeDYJetsToLL_M50=0;
  // TreeReader * tree_DYJetsToLL_M50 = new TreeReader(treeDYJetsToLL_M50, "DYJetsToLL_M50", systlist);
  // tree_DYJetsToLL_M50->Loop("DYJetsToLL_M50", false);
  // delete tree_DYJetsToLL_M50;

  // TTree* treeST_tW_antitop_5f_NoFullyHadronicDecays=0;
  // TreeReader * tree_ST_tW_antitop_5f_NoFullyHadronicDecays = new TreeReader(treeST_tW_antitop_5f_NoFullyHadronicDecays, "ST_tW_antitop_5f_NoFullyHadronicDecays", systlist);
  // tree_ST_tW_antitop_5f_NoFullyHadronicDecays->Loop("ST_tW_antitop_5f_NoFullyHadronicDecays", false);
  // delete tree_ST_tW_antitop_5f_NoFullyHadronicDecays;

  // TTree* treeTTTo2L2Nu=0;
  // TreeReader * tree_TTTo2L2Nu = new TreeReader(treeTTTo2L2Nu, "TTTo2L2Nu", systlist);
  // tree_TTTo2L2Nu->Loop("TTTo2L2Nu", false);
  // delete tree_TTTo2L2Nu;

}
