{
  gROOT->ProcessLine(".L TreeReader.C+");
  
  
  std::vector<TString > systlist;
  systlist.push_back("");
  
  //Signal

  TTree* tree70=0;
  TreeReader * tree_70_test = new TreeReader(tree70, "70_test", systlist);
  tree_70_test->Loop("70_test");
  delete tree_70_test;

  TTree* tree50=0;
  TreeReader * tree_50_test = new TreeReader(tree50, "50_test", systlist);
  tree_50_test->Loop("50_test");
  delete tree_50_test;

  TTree* tree30=0;
  TreeReader * tree_30_test = new TreeReader(tree30, "30_test", systlist);
  tree_30_test->Loop("30_test");
  delete tree_30_test;

  TTree* tree10=0;
  TreeReader * tree_10_test = new TreeReader(tree10, "10_test", systlist);
  tree_10_test->Loop("10_test");
  delete tree_10_test;


//----------//


    TTree* treesmu200_neu180_lamE_2_ctau20=0;
  TreeReader * tree_smu200_neu180_lamE_2_ctau20 = new TreeReader(treesmu200_neu180_lamE_2_ctau20, "smu200_neu180_lamE_2_ctau20", systlist);
  tree_smu200_neu180_lamE_2_ctau20->Loop("smu200_neu180_lamE_2_ctau20");
  delete tree_smu200_neu180_lamE_2_ctau20;

    TTree* treesmu250_neu200_lamE_2_ctau20=0;
  TreeReader * tree_smu250_neu200_lamE_2_ctau20 = new TreeReader(treesmu250_neu200_lamE_2_ctau20, "smu250_neu200_lamE_2_ctau20", systlist);
  tree_smu250_neu200_lamE_2_ctau20->Loop("smu250_neu200_lamE_2_ctau20");
  delete tree_smu250_neu200_lamE_2_ctau20;

    TTree* treesmu300_neu180_lamE_2_ctau20=0;
  TreeReader * tree_smu300_neu180_lamE_2_ctau20 = new TreeReader(treesmu300_neu180_lamE_2_ctau20, "smu300_neu180_lamE_2_ctau20", systlist);
  tree_smu300_neu180_lamE_2_ctau20->Loop("smu300_neu180_lamE_2_ctau20");
  delete tree_smu300_neu180_lamE_2_ctau20;

    TTree* treesmu300_neu200_lamE_2_ctau20=0;
  TreeReader * tree_smu300_neu200_lamE_2_ctau20 = new TreeReader(treesmu300_neu200_lamE_2_ctau20, "smu300_neu200_lamE_2_ctau20", systlist);
  tree_smu300_neu200_lamE_2_ctau20->Loop("smu300_neu200_lamE_2_ctau20");
  delete tree_smu300_neu200_lamE_2_ctau20;

    TTree* treesmu300_neu250_lamE_2_ctau01=0;
  TreeReader * tree_smu300_neu250_lamE_2_ctau01 = new TreeReader(treesmu300_neu250_lamE_2_ctau01, "smu300_neu250_lamE_2_ctau01", systlist);
  tree_smu300_neu250_lamE_2_ctau01->Loop("smu300_neu250_lamE_2_ctau01");
  delete tree_smu300_neu250_lamE_2_ctau01;

    TTree* treesmu300_neu250_lamE_2_ctau03=0;
  TreeReader * tree_smu300_neu250_lamE_2_ctau03 = new TreeReader(treesmu300_neu250_lamE_2_ctau03, "smu300_neu250_lamE_2_ctau03", systlist);
  tree_smu300_neu250_lamE_2_ctau03->Loop("smu300_neu250_lamE_2_ctau03");
  delete tree_smu300_neu250_lamE_2_ctau03;



    TTree* treeNtuple_smu300_neu250_lamE_2_ctau10=0;
  TreeReader * tree_Ntuple_smu300_neu250_lamE_2_ctau10 = new TreeReader(treeNtuple_smu300_neu250_lamE_2_ctau10, "Ntuple_smu300_neu250_lamE_2_ctau10", systlist);
  tree_Ntuple_smu300_neu250_lamE_2_ctau10->Loop("Ntuple_smu300_neu250_lamE_2_ctau10");
  delete tree_Ntuple_smu300_neu250_lamE_2_ctau10;

  TTree* treesmu300_neu250_lamE_2_ctau100=0;
  TreeReader * tree_smu300_neu250_lamE_2_ctau100 = new TreeReader(treesmu300_neu250_lamE_2_ctau100, "smu300_neu250_lamE_2_ctau100", systlist);
  tree_smu300_neu250_lamE_2_ctau100->Loop("smu300_neu250_lamE_2_ctau100");
  delete tree_smu300_neu250_lamE_2_ctau100;

    TTree* treesmu300_neu250_lamE_2_ctau20=0;
  TreeReader * tree_smu300_neu250_lamE_2_ctau20 = new TreeReader(treesmu300_neu250_lamE_2_ctau20, "smu300_neu250_lamE_2_ctau20", systlist);
  tree_smu300_neu250_lamE_2_ctau20->Loop("smu300_neu250_lamE_2_ctau20");
  delete tree_smu300_neu250_lamE_2_ctau20;

    TTree* treesmu300_neu250_lamE_2_ctau30=0;
  TreeReader * tree_smu300_neu250_lamE_2_ctau30 = new TreeReader(treesmu300_neu250_lamE_2_ctau30, "smu300_neu250_lamE_2_ctau30", systlist);
  tree_smu300_neu250_lamE_2_ctau30->Loop("smu300_neu250_lamE_2_ctau30");
  delete tree_smu300_neu250_lamE_2_ctau30;

    TTree* treesmu300_neu250_lamE_2_ctaup1=0;
  TreeReader * tree_smu300_neu250_lamE_2_ctaup1 = new TreeReader(treesmu300_neu250_lamE_2_ctaup1, "smu300_neu250_lamE_2_ctaup1", systlist);
  tree_smu300_neu250_lamE_2_ctaup1->Loop("smu300_neu250_lamE_2_ctaup1");
  delete tree_smu300_neu250_lamE_2_ctaup1;

    TTree* treesmu300_neu250_lamE_2_ctaup3=0;
  TreeReader * tree_smu300_neu250_lamE_2_ctaup3 = new TreeReader(treesmu300_neu250_lamE_2_ctaup3, "smu300_neu250_lamE_2_ctaup3", systlist);
  tree_smu300_neu250_lamE_2_ctaup3->Loop("smu300_neu250_lamE_2_ctaup3");
  delete tree_smu300_neu250_lamE_2_ctaup3;

    TTree* treesmu300_neu280_lamE_2_ctau20=0;
  TreeReader * tree_smu300_neu280_lamE_2_ctau20 = new TreeReader(treesmu300_neu280_lamE_2_ctau20, "smu300_neu280_lamE_2_ctau20", systlist);
  tree_smu300_neu280_lamE_2_ctau20->Loop("smu300_neu280_lamE_2_ctau20");
  delete tree_smu300_neu280_lamE_2_ctau20;


    TTree* treesmu400_neu180_lamE_2_ctau20=0;
  TreeReader * tree_smu400_neu180_lamE_2_ctau20 = new TreeReader(treesmu400_neu180_lamE_2_ctau20, "smu400_neu180_lamE_2_ctau20", systlist);
  tree_smu400_neu180_lamE_2_ctau20->Loop("smu400_neu180_lamE_2_ctau20");
  delete tree_smu400_neu180_lamE_2_ctau20;

    TTree* treesmu400_neu200_lamE_2_ctau20=0;
  TreeReader * tree_smu400_neu200_lamE_2_ctau20 = new TreeReader(treesmu400_neu200_lamE_2_ctau20, "smu400_neu200_lamE_2_ctau20", systlist);
  tree_smu400_neu200_lamE_2_ctau20->Loop("smu400_neu200_lamE_2_ctau20");
  delete tree_smu400_neu200_lamE_2_ctau20;



    TTree* treesmu400_neu250_lamE_2_ctau20=0;
  TreeReader * tree_smu400_neu250_lamE_2_ctau20 = new TreeReader(treesmu400_neu250_lamE_2_ctau20, "smu400_neu250_lamE_2_ctau20", systlist);
  tree_smu400_neu250_lamE_2_ctau20->Loop("smu400_neu250_lamE_2_ctau20");
  delete tree_smu400_neu250_lamE_2_ctau20;

    TTree* treesmu400_neu300_lamE_2_ctau20=0;
  TreeReader * tree_smu400_neu300_lamE_2_ctau20 = new TreeReader(treesmu400_neu300_lamE_2_ctau20, "smu400_neu300_lamE_2_ctau20", systlist);
  tree_smu400_neu300_lamE_2_ctau20->Loop("smu400_neu300_lamE_2_ctau20");
  delete tree_smu400_neu300_lamE_2_ctau20;

    TTree* treesmu400_neu350_lamE_2_ctau20=0;
  TreeReader * tree_smu400_neu350_lamE_2_ctau20 = new TreeReader(treesmu400_neu350_lamE_2_ctau20, "smu400_neu350_lamE_2_ctau20", systlist);
  tree_smu400_neu350_lamE_2_ctau20->Loop("smu400_neu350_lamE_2_ctau20");
  delete tree_smu400_neu350_lamE_2_ctau20;



    TTree* treesmu400_neu380_lamE_2_ctau20=0;
  TreeReader * tree_smu400_neu380_lamE_2_ctau20 = new TreeReader(treesmu400_neu380_lamE_2_ctau20, "smu400_neu380_lamE_2_ctau20", systlist);
  tree_smu400_neu380_lamE_2_ctau20->Loop("smu400_neu380_lamE_2_ctau20");
  delete tree_smu400_neu380_lamE_2_ctau20;

    TTree* treesmu500_neu180_lamE_2_ctau20=0;
  TreeReader * tree_smu500_neu180_lamE_2_ctau20 = new TreeReader(smu500_neu180_lamE_2_ctau20, "smu500_neu180_lamE_2_ctau20", systlist);
  tree_smu500_neu180_lamE_2_ctau20->Loop("smu500_neu180_lamE_2_ctau20");
  delete tree_smu500_neu180_lamE_2_ctau20;

    TTree* treesmu500_neu200_lamE_2_ctau20=0;
  TreeReader * tree_smu500_neu200_lamE_2_ctau20 = new TreeReader(treesmu500_neu200_lamE_2_ctau20, "smu500_neu200_lamE_2_ctau20", systlist);
  tree_smu500_neu200_lamE_2_ctau20->Loop("smu500_neu200_lamE_2_ctau20");
  delete tree_smu500_neu200_lamE_2_ctau20;

    TTree* treesmu500_neu250_lamE_2_ctau20=0;
  TreeReader * tree_smu500_neu250_lamE_2_ctau20 = new TreeReader(treesmu500_neu250_lamE_2_ctau20, "smu500_neu250_lamE_2_ctau20", systlist);
  tree_smu500_neu250_lamE_2_ctau20->Loop("smu500_neu250_lamE_2_ctau20");
  delete tree_smu500_neu250_lamE_2_ctau20;

    TTree* treesmu500_neu300_lamE_2_ctau20=0;
  TreeReader * tree_smu500_neu300_lamE_2_ctau20 = new TreeReader(treesmu500_neu300_lamE_2_ctau20, "smu500_neu300_lamE_2_ctau20", systlist);
  tree_smu500_neu300_lamE_2_ctau20->Loop("smu500_neu300_lamE_2_ctau20");
  delete tree_smu500_neu300_lamE_2_ctau20;

    TTree* treesmu500_neu350_lamE_2_ctau20=0;
  TreeReader * tree_smu500_neu350_lamE_2_ctau20 = new TreeReader(treesmu500_neu350_lamE_2_ctau20, "smu500_neu350_lamE_2_ctau20", systlist);
  tree_smu500_neu350_lamE_2_ctau20->Loop("smu500_neu350_lamE_2_ctau20");
  delete tree_smu500_neu350_lamE_2_ctau20;

    TTree* treesmu500_neu400_lamE_2_ctau20=0;
  TreeReader * tree_smu500_neu400_lamE_2_ctau20 = new TreeReader(treesmu500_neu400_lamE_2_ctau20, "smu500_neu400_lamE_2_ctau20", systlist);
  tree_smu500_neu400_lamE_2_ctau20->Loop("smu500_neu400_lamE_2_ctau20");
  delete tree_smu500_neu400_lamE_2_ctau20;

    TTree* treesmu500_neu450_lamE_2_ctau20=0;
  TreeReader * tree_smu500_neu450_lamE_2_ctau20 = new TreeReader(treesmu500_neu450_lamE_2_ctau20, "smu500_neu450_lamE_2_ctau20", systlist);
  tree_smu500_neu450_lamE_2_ctau20->Loop("smu500_neu450_lamE_2_ctau20");
  delete tree_smu500_neu450_lamE_2_ctau20;

    TTree* treesmu500_neu480_lamE_2_ctau20=0;
  TreeReader * tree_smu500_neu480_lamE_2_ctau20 = new TreeReader(treesmu500_neu480_lamE_2_ctau20, "smu500_neu480_lamE_2_ctau20", systlist);
  tree_smu500_neu480_lamE_2_ctau20->Loop("smu500_neu480_lamE_2_ctau20");
  delete tree_smu500_neu480_lamE_2_ctau20;
  
  //Background 

  TTree* treeWZTo2Q2L_mllmin4p0=0;
  TreeReader * tree_WZTo2Q2L_mllmin4p0 = new TreeReader(treeWZTo2Q2L_mllmin4p0, "WZTo2Q2L_mllmin4p0", systlist);
  tree_WZTo2Q2L_mllmin4p0->Loop("WZTo2Q2L_mllmin4p0");
  delete tree_WZTo2Q2L_mllmin4p0;

  TTree* treeWWTo2L2Nu=0;
  TreeReader * tree_WWTo2L2Nu = new TreeReader(treeWWTo2L2Nu, "WWTo2L2Nu", systlist);
  tree_WWTo2L2Nu->Loop("WWTo2L2Nu");
  delete tree_WWTo2L2Nu;

  TTree* treeWWTo2L2Nu_MLL_600To1200=0;
  TreeReader * tree_WWTo2L2Nu_MLL_600To1200 = new TreeReader(treeWWTo2L2Nu_MLL_600To1200, "WWTo2L2Nu_MLL_600To1200", systlist);
  tree_WWTo2L2Nu_MLL_600To1200->Loop("WWTo2L2Nu_MLL_600To1200");
  delete tree_WWTo2L2Nu_MLL_600To1200;

  TTree* treeWWTo2L2Nu_MLL_200To600=0;
  TreeReader * tree_WWTo2L2Nu_MLL_200To600 = new TreeReader(treeWWTo2L2Nu_MLL_200To600, "WWTo2L2Nu_MLL_200To600", systlist);
  tree_WWTo2L2Nu_MLL_200To600->Loop("WWTo2L2Nu_MLL_200To600");
  delete tree_WWTo2L2Nu_MLL_200To600;

  TTree* treeTTJets_DiLept=0;
  TreeReader * tree_TTJets_DiLept = new TreeReader(treeTTJets_DiLept, "TTJets_DiLept", systlist);
  tree_TTJets_DiLept->Loop("TTJets_DiLept");
  delete tree_TTJets_DiLept;

  TTree* treeZZTo2Q2L_mllmin4p0=0;
  TreeReader * tree_ZZTo2Q2L_mllmin4p0 = new TreeReader(treeZZTo2Q2L_mllmin4p0, "ZZTo2Q2L_mllmin4p0", systlist);
  tree_ZZTo2Q2L_mllmin4p0->Loop("ZZTo2Q2L_mllmin4p0");
  delete tree_ZZTo2Q2L_mllmin4p0;

  TTree* treeST_tW_top_5f_NoFullyHadronicDecays=0;
  TreeReader * tree_ST_tW_top_5f_NoFullyHadronicDecays = new TreeReader(treeST_tW_top_5f_NoFullyHadronicDecays, "ST_tW_top_5f_NoFullyHadronicDecays", systlist);
  tree_ST_tW_top_5f_NoFullyHadronicDecays->Loop("ST_tW_top_5f_NoFullyHadronicDecays");
  delete tree_ST_tW_top_5f_NoFullyHadronicDecays;
  
}
