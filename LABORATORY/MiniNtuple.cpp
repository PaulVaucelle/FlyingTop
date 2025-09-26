
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> <output_file.root>" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];

    // ROOT::EnableVerboseLogging(true);
    ROOT::EnableImplicitMT(2); // Multi-threading
    ROOT::RDataFrame df("FlyingTop/ttree", input_file);
    TFile output(output_file.c_str(), "RECREATE");
// g++ -o MiniDataMCNtuple MiniNtuple.cpp $(root-config --cflags --libs) -lROOTDataFrame
// ./MiniDataMCNtuple <<inputfile.root>> <name of outputfile.root>

// g++ -o MiniNtuple MiniNtuple.cpp $(root-config --cflags --libs) -lROOTDataFrame
// ./MiniNtuple /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/MC_EMU_03_02_2025/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root ./MiniDATAMC_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root
//-------tree_Filter-------//
 
// auto filter_df = df.Filter("tree_Filter == true","Filtre simple");
// auto hData_njetNOmu_Filter = filter_df.Histo1D({"hData_njetNOmu_Filter", "Histogramme original;nJetNomu;Counts", 15, 0, 15}, "tree_njetNOmu");
// hData_njetNOmu_Filter->Write();

//------tree_FIlterSamSign-----//
// auto filterSameSign_df = df.Filter("tree_FilterSameSign == true","Filtre simple");

//-------Filtre Global--------//

    // filtered_df.Snapshot("ttree", output_file);
    // auto filtered_df_small = df.Range(1000000);
    auto filtered_df = df.Filter("(tree_Filter == true  && tree_njetNOmu > 0) ", "Filtre global combine");

    filtered_df.Snapshot("ttree", output_file,
     {
        // "eventNumber",
         "tree_LHE_Weights", 
        // "tree_MCEvt_weight", 
        "tree_only_gen_wt", 
    //  "tree_event_weight", 
     "tree_genTop_Weight", "tree_gen_top_pt", 
     "PUweight", "PUweight_Up", "PUweight_Down", 
      "Prefweight_Down","Prefweight","Prefweight_Up",

    //  "tree_Filter", "tree_FilterSameSign", "tree_trigger_doublelepton", "tree_trigger_singlelepton", 
    //     "tree_Good_PV", "tree_nPV", 
    //   "tree_HT", "tree_TRACK_SIZE", "tree_nTracks", "tree_nLostTracks", 
    //   "tree_all_nmu", "tree_nmu", "tree_LT", 
    "tree_Mmumu", 
    // "tree_MmumuSameSign",
    //   "tree_muon_pt",  "tree_muon_eta", "tree_muon_dxy", "tree_muon_dz", 
        "tree_lepton_leadingpt", "tree_lepton_leadingpt2", "tree_lepton_leadingeta", 
      "tree_lepton_leadingeta2",
    //    "tree_lepton_leadingphi", "tree_lepton_leadingphi2",
    //    "tree_lepton_lepton_dR", 
      "tree_PFMet_et"
    //   "tree_lepton_lepton_dPhi", "tree_lepton_lepton_dEta", "tree_lepton_leadingdxy", "tree_lepton_leadingdxy2", 
    //   "tree_lepton_leadingdz", "tree_lepton_leadingdz2", "tree_all_nel", "tree_electron_nEle", "tree_electron_isPrompt", 
    //   "tree_electron_pt", "tree_electron_eta", 
    //    "tree_njet", "tree_njetNOmu", "tree_jet_pt", "tree_jet_eta", "tree_jet_leadingpt", "tree_jet_leadingpt2", "tree_jet_leadingeta", 
    //    "tree_jet_leadingeta2", "tree_jet_jet_dR", "tree_jet_jet_dPhi", "tree_jet_jet_dEta", "tree_muon_jet_dRmin",
    //     "tree_muon_jet_dRmax", "tree_Hemi", 
    //     "tree_Hemi_njet", "tree_Hemi_njet_nomu", "tree_Hemi_pt", "tree_Hemi_eta", "tree_Hemi_phi", 
    //     "tree_Hemi_nTrks",  "tree_Hemi_mass", "tree_HemiMu_mass", 
    //     "tree_HemiMu_pt", "tree_HemiMu_dR", "tree_HemiMuOp_mass", "tree_HemiMuOp_pt", "tree_HemiMuOp_dR", 
    //     "tree_Hemi_dR12", "tree_Hemi_Vtx_step", "tree_Hemi_Vtx_isTight", "tree_Hemi_Vtx_NChi2", "tree_Hemi_Vtx_nTrks",
    //     "tree_Hemi_Vtx_r", "tree_Hemi_Vtx_dR", "tree_Hemi_Vtx_SumtrackWeight",
    //     "tree_Hemi_Vtx_track_MeanDCA_d", "tree_Hemi_Vtx_Mass", "tree_Hemi_Vtx_dist",  "tree_Hemi_SecVtx", 
    //     "tree_Hemi_SecVtx_step", "tree_Hemi_SecVtx_x", "tree_Hemi_SecVtx_y", "tree_Hemi_SecVtx_z", 
    //     "tree_Hemi_SecVtx_r", "tree_Hemi_SecVtx_dR", "tree_Hemi_SecVtx_nTrks", "tree_Hemi_SecVtx_NChi2", 
    //     "tree_Hemi_SecVtx_dist", "tree_Hemi_SecVtx_track_MeanDCA_d", "tree_Hemi_SecVtx_SumtrackWeight", 
    //     "tree_Hemi_SecVtx_Mass", "tree_Hemi_Vtx_BDT_nTrks", "tree_Hemi_Vtx_BDT_NChi2", 
    //     "tree_Hemi_Vtx_BDT_step", "tree_Hemi_Vtx_BDT_STW", "tree_Hemi_Vtx_BDT_Mass", "tree_Hemi_Vtx_BDT_HMass",
    //     "tree_Hemi_Vtx_BDT_ntrk10", "tree_Hemi_Vtx_BDT_ntrk20", "tree_Hemi_Vtx_BDT_MeanDCA", "tree_Hemi_Vtx_MVAval_Loose", "tree_Hemi_Vtx_MVAval_Tight",
           }
           
           );

std::cout << "Analyse terminee. TTree sauvegarde dans " << output_file << std::endl;

    return 0;
}
