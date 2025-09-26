
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

    // ROOT::EnableImplicitMT(8); // Multi-threading
    // ROOT::RDF::PerfStats::Enable();


    ROOT::RDataFrame df("FlyingTop/ttree", input_file);
    TFile output(output_file.c_str(), "RECREATE");

//-------tree_Filter-------//
 
// auto filter_df = df.Filter("tree_Filter == true","Filtre simple");

// auto hData_njetNOmu_Filter = filter_df.Histo1D({"hData_njetNOmu_Filter", "Histogramme original;nJetNomu;Counts", 15, 0, 15}, "tree_njetNOmu");
// hData_njetNOmu_Filter->Write();

//------tree_FIlterSamSign-----//
// auto filterSameSign_df = df.Filter("tree_FilterSameSign == true","Filtre simple");

//-------Filtre Global--------//
    auto filtered_df = df.Filter("(tree_Filter == true  && tree_njetNOmu > 0) ", "Filtre global combine");
    filtered_df.Snapshot("ttree", output_file);


    //-- Renaming of the branches to match the other script that generaters minintuples
    // auto renamed_df = filtered_df

    // renamed_df.Snapshot("ttree", output_file, {"minieventNumber", "minilumiBlock", "minitree_LHE_Weights", "minitree_MCEvt_weight", "minitree_only_gen_wt", "minitree_event_weight", "minitree_genTop_Weight", "minitree_gen_top_pt", "minitree_gen_top_rw_pt", "miniPUweight", "miniPUweight_Up", "miniPUweight_Down", "miniPrefweight", "miniPU_events", "minitree_Filter", "minitree_FilterSameSign", "minitree_trigger_doublelepton", "minitree_trigger_singlelepton", "minitree_GenPVx", "minitree_GenPVy", "minitree_GenPVz", "minitree_smu_mass", "minitree_neu_mass", "minitree_neu_ctau", "minitree_Good_PV", "minitree_nPV", "minitree_PV_x", "minitree_PV_y", "minitree_PV_z", "minitree_PV_ez", "minitree_PV_NChi2", "minitree_PV_ndf", "minitree_PFMet_et", "minitree_PFMet_phi", "minitree_HT", "minitree_TRACK_SIZE", "minitree_nTracks", "minitree_nLostTracks", "minitree_muon_GenRecoTriggerMatched", "minitree_all_nmu", "minitree_nmu", "minitree_LT", "minitree_Mmumu", "minitree_MmumuSameSign", "minitree_muon_isPrompt", "minitree_muon_pt", "minitree_muon_SF", "minitree_muon_eta", "minitree_muon_phi", "minitree_muon_dxy", "minitree_muon_dz", "minitree_muon_charge", "minitree_muon_correction", "minitree_muon_gen", "minitree_muon_dxyError", "minitree_muon_dzError", "minitree_muon_isLoose", "minitree_muon_isMedium", "minitree_muon_isTight", "minitree_muon_isGlobal", "minitree_muon_PFIsoVeryLoose", "minitree_muon_PFIsoLoose", "minitree_muon_PFIsoMedium", "minitree_muon_PFIsoTight", "minitree_muon_TkIsoLoose", "minitree_muon_TkIsoTight", "minitree_muon_MiniIsoLoose", "minitree_muon_MiniIsoMedium", "minitree_muon_MiniIsoTight", "minitree_lepton_leadingpt", "minitree_lepton_leadingpt2", "minitree_lepton_leadingeta", "minitree_lepton_leadingeta2", "minitree_lepton_leadingphi", "minitree_lepton_leadingphi2", "minitree_lepton_lepton_dR", "minitree_lepton_lepton_dPhi", "minitree_lepton_lepton_dEta", "minitree_lepton_leadingdxy", "minitree_lepton_leadingdxy2", "minitree_lepton_leadingdz", "minitree_lepton_leadingdz2", "minitree_all_nel", "minitree_electron_nEle", "minitree_electron_isPrompt", "minitree_electron_pt", "minitree_electron_eta", "minitree_electron_phi", "minitree_electron_charge", "minitree_electron_dxy", "minitree_electron_dz", "minitree_electron_gen", "minitree_electron_energy", "minitree_electron_et", "minitree_electron_ecal_trk_postcorr", "minitree_electron_isoR4", "minitree_electron_IsLoose", "minitree_electron_IsMedium", "minitree_electron_IsTight", "minitree_njet", "minitree_njetNOmu", "minitree_jet_pt", "minitree_jet_eta", "minitree_jet_phi", "minitree_jet_HadronFlavour", "minitree_jet_btag_DeepJet", "minitree_jet_E", "minitree_jet_leadingpt", "minitree_jet_leadingpt2", "minitree_jet_leadingeta", "minitree_jet_leadingeta2", "minitree_jet_jet_dR", "minitree_jet_jet_dPhi", "minitree_jet_jet_dEta", "minitree_muon_jet_dRmin", "minitree_muon_jet_dRmax", "minitree_Evts_MVAval", "minitree_Evts_MVAvalDY", "minitree_Evts_MVAvalTT", "minitree_Hemi", "minitree_Hemi_njet", "minitree_Hemi_njet_nomu", "minitree_Hemi_pt", "minitree_Hemi_eta", "minitree_Hemi_phi", "minitree_Hemi_nTrks", "minitree_Hemi_nTrks_sig", "minitree_Hemi_nTrks_bad", "minitree_Hemi_mass", "minitree_HemiMu_mass", "minitree_HemiMu_pt", "minitree_HemiMu_dR", "minitree_HemiMuOp_mass", "minitree_HemiMuOp_pt", "minitree_HemiMuOp_dR", "minitree_Hemi_dR12", "minitree_Hemi_Vtx_step", "minitree_Hemi_Vtx_isTight", "minitree_Hemi_Vtx_NChi2", "minitree_Hemi_Vtx_nTrks", "minitree_Hemi_Vtx_nTrks_sig", "minitree_Hemi_Vtx_nTrks_bad", "minitree_Hemi_Vtx_x", "minitree_Hemi_Vtx_y", "minitree_Hemi_Vtx_z", "minitree_Hemi_Vtx_r", "minitree_Hemi_Vtx_dR", "minitree_Hemi_Vtx_SumtrackWeight", "minitree_Hemi_Vtx_track_MeanDCA_d", "minitree_Hemi_Vtx_Mass", "minitree_Hemi_Vtx_dist", "minitree_event_nVtx", "minitree_event_Vtx_Vtx_dr", "minitree_event_Vtx_Vtx_dz", "minitree_event_Vtx_Vtx_dd", "minitree_event_Vtx_Vtx_reldd", "minitree_event_Vtx_Vtx_dR", "minitree_event_Vtx_Vtx_step", "minitree_Hemi_SecLLP", "minitree_Hemi_LLP_SecVtx_dz", "minitree_Hemi_LLP_SecVtx_dr", "minitree_Hemi_SecLLP_ping", "minitree_event_SecLLP_ping", "minitree_Hemi_SecVtx", "minitree_Hemi_SecVtx_step", "minitree_Hemi_SecVtx_x", "minitree_Hemi_SecVtx_y", "minitree_Hemi_SecVtx_z", "minitree_Hemi_SecVtx_r", "minitree_Hemi_SecVtx_dR", "minitree_Hemi_SecVtx_nTrks", "minitree_Hemi_SecVtx_NChi2", "minitree_Hemi_SecVtx_dist", "minitree_Hemi_SecVtx_track_MeanDCA_d", "minitree_Hemi_SecVtx_SumtrackWeight", "minitree_Hemi_SecVtx_Mass", "minitree_event_MergedVtx_Vtx_dr", "minitree_event_MergedVtx_Vtx_dz", "minitree_event_MergedVtx_Vtx_dd", "minitree_event_MergedVtx_Vtx_reldd", "minitree_event_MergedVtx_Vtx_dR", "minitree_event_MergedVtx_Vtx_step", "minitree_Hemi_Vtx_BDT_nTrks", "minitree_Hemi_Vtx_BDT_NChi2", "minitree_Hemi_Vtx_BDT_step", "minitree_Hemi_Vtx_BDT_STW", "minitree_Hemi_Vtx_BDT_Mass", "minitree_Hemi_Vtx_BDT_HMass", "minitree_Hemi_Vtx_BDT_ntrk10", "minitree_Hemi_Vtx_BDT_ntrk20", "minitree_Hemi_Vtx_BDT_MeanDCA", "minitree_Hemi_Vtx_MVAval_Loose", "minitree_Hemi_Vtx_MVAval_Tight"});

    std::cout << "Analyse terminee. TTree sauvegarde dans " << output_file << std::endl;

            // auto df_small = df.Range(10000);
            // df_small.Snapshot("ttree", output_file);
    return 0;
}
