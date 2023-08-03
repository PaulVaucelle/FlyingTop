{
  gROOT->ProcessLine(".L TreeReader.C+");
  
  
  std::vector<TString > systlist;
  systlist.push_back("");
  
  // TTree* tree70=0;
  // TreeReader * tree_70_test = new TreeReader(tree70, "70_test", systlist);
  // tree_70_test->Loop("70_test");
  // delete tree_70_test;

  TTree* tree50=0;
  TreeReader * tree_50_test = new TreeReader(tree50, "50_test", systlist);
  tree_50_test->Loop("50_test");
  delete tree_50_test;

  // TTree* tree30=0;
  // TreeReader * tree_30_test = new TreeReader(tree30, "30_test", systlist);
  // tree_30_test->Loop("30_test");
  // delete tree_30_test;

  // TTree* tree10=0;
  // TreeReader * tree_10_test = new TreeReader(tree10, "10_test", systlist);
  // tree_10_test->Loop("10_test");
  // delete tree_10_test;
  
}
