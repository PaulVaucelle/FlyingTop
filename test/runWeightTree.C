{
    gROOT->ProcessLine(".L WeightTree.C");
    // WeightTree("10_test");
    // WeightTree("30_test");
    // WeightTree("50_test");
    // WeightTree("70_test");
    WeightTree("smu500_neu200_lamE_2_ctau20");

    // WeightTree("ST_tW_antitop_5f_NoFullyHadronicDecays");
    // WeightTree("ST_tW_top_5f_NoFullyHadronicDecays");
    // WeightTree("DYJetsToLL_M10to50");
    // WeightTree("DYJetsToLL_M50");
    // WeightTree("WWTo2L2Nu_MLL_200To600");
    // WeightTree("WWTo2L2Nu_MLL_600To1200");
    // WeightTree("WWTo2L2Nu");
    // WeightTree("TTJets_DiLept");
    // WeightTree("TTTo2L2Nu");
    // WeightTree("TTWW");
    // WeightTree("ttWJetsToLNu_5f_EWK");
    // WeightTree("TTZToLL_5f");
    // WeightTree("WZTo2Q2L_mllmin4p0");
    // WeightTree("ZZTo2Q2L_mllmin4p0");

    // gROOT->ProcessLine("hadd -f Ntuple_Signal_test.root Ntuple_10_test_weighted.root Ntuple_30_test_weighted.root Ntuple_50_test_weighted.root Ntuple_70_test_weighted.root");

}
