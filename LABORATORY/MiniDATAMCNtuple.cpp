
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

    ROOT::EnableVerboseLogging(true);
    ROOT::EnableImplicitMT(16); // Multi-threading
    ROOT::RDataFrame df("FlyingTop/ttree", input_file);
    TFile output(output_file.c_str(), "RECREATE");
// g++ -o MiniDataMCNtuple MiniNtuple.cpp $(root-config --cflags --libs) -lROOTDataFrame
// ./MiniDataMCNtuple <<inputfile.root>> <name of outputfile.root>

// g++ -o MiniDataMCNtuple MiniNtuple.cpp $(root-config --cflags --libs) -lROOTDataFrame
// ./MiniDataMCNtuple /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/MC_EMU_03_02_2025/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root ./MiniDATAMC_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root

    // filtered_df.Snapshot("ttree", output_file);
    // auto filtered_df_small = df.Range(1000000);
    auto filtered_df = df.Filter("(tree_Filter == true  && tree_njetNOmu > 0) ", "Filtre global combine");


    // filtered_df.Range(1000000);
//-- Renaming of the branches to match the other script that generaters minintuples
// auto renamed_df = filtered_df
//     .Define("minieventNumber", "eventNumber")
//     .Define("minilumiBlock", "lumiBlock")

//     .Define("miniPrefweight_Up","Prefweight_Up")
//     .Define("miniPrefweight_Down","Prefweight_Down")
//     .Define("minitree_LHE_Weights", "tree_LHE_Weights")
//     .Define("minitree_MCEvt_weight", "tree_MCEvt_weight")
//     .Define("minitree_only_gen_wt", "tree_only_gen_wt")
//     .Define("minitree_event_weight", "tree_event_weight")
//     .Define("minitree_genTop_Weight", "tree_genTop_Weight")
//     .Define("minitree_gen_top_pt", "tree_gen_top_pt")
//     .Define("minitree_gen_top_rw_pt", "tree_gen_top_rw_pt")
//     .Define("miniPUweight", "PUweight")
//     .Define("miniPUweight_Up", "PUweight_Up")
//     .Define("miniPUweight_Down", "PUweight_Down")
//     .Define("miniPrefweight", "Prefweight")
//     .Define("miniPU_events", "PU_events")
//     .Define("minitree_Filter", "tree_Filter")
//     .Define("minitree_FilterSameSign", "tree_FilterSameSign")
//     .Define("minitree_trigger_doublelepton", "tree_trigger_doublelepton")
//     .Define("minitree_trigger_singlelepton", "tree_trigger_singlelepton")
//     .Define("minitree_Good_PV", "tree_Good_PV")
//     .Define("minitree_nPV", "tree_nPV")
//     .Define("minitree_PV_x", "tree_PV_x")
//     .Define("minitree_PV_y", "tree_PV_y")
//     .Define("minitree_PV_z", "tree_PV_z")
//     .Define("minitree_PV_ez", "tree_PV_ez")
//     .Define("minitree_PV_NChi2", "tree_PV_NChi2")
//     .Define("minitree_PV_ndf", "tree_PV_ndf")
//     .Define("minitree_PFMet_et", "tree_PFMet_et")
//     .Define("minitree_PFMet_phi", "tree_PFMet_phi")
//     .Define("minitree_HT", "tree_HT")
//     .Define("minitree_TRACK_SIZE", "tree_TRACK_SIZE")
//     .Define("minitree_nTracks", "tree_nTracks")
//     .Define("minitree_nLostTracks", "tree_nLostTracks")
//     .Define("minitree_muon_GenRecoTriggerMatched", "tree_muon_GenRecoTriggerMatched")
//     .Define("minitree_all_nmu", "tree_all_nmu")
//     .Define("minitree_nmu", "tree_nmu")
//     .Define("minitree_LT", "tree_LT")
//     .Define("minitree_Mmumu", "tree_Mmumu")
//     .Define("minitree_MmumuSameSign", "tree_MmumuSameSign")
//     .Define("minitree_muon_isPrompt", "tree_muon_isPrompt")
//     .Define("minitree_muon_pt", "tree_muon_pt")
//     .Define("minitree_muon_SF", "tree_muon_SF")
//     .Define("minitree_muon_eta", "tree_muon_eta")
//     .Define("minitree_muon_phi", "tree_muon_phi")
//     .Define("minitree_muon_dxy", "tree_muon_dxy")
//     .Define("minitree_muon_dz", "tree_muon_dz")
//     .Define("minitree_muon_charge", "tree_muon_charge")
//     .Define("minitree_muon_correction", "tree_muon_correction")
//     .Define("minitree_muon_gen", "tree_muon_gen")
//     .Define("minitree_muon_dxyError", "tree_muon_dxyError")
//     .Define("minitree_muon_dzError", "tree_muon_dzError")
//     .Define("minitree_muon_isLoose", "tree_muon_isLoose")
//     .Define("minitree_muon_isMedium", "tree_muon_isMedium")
//     .Define("minitree_muon_isTight", "tree_muon_isTight")
//     .Define("minitree_muon_isGlobal", "tree_muon_isGlobal")
//     .Define("minitree_muon_PFIsoVeryLoose", "tree_muon_PFIsoVeryLoose")
//     .Define("minitree_muon_PFIsoLoose", "tree_muon_PFIsoLoose")
//     .Define("minitree_muon_PFIsoMedium", "tree_muon_PFIsoMedium")
//     .Define("minitree_muon_PFIsoTight", "tree_muon_PFIsoTight")
//     .Define("minitree_muon_TkIsoLoose", "tree_muon_TkIsoLoose")
//     .Define("minitree_muon_TkIsoTight", "tree_muon_TkIsoTight")
//     .Define("minitree_muon_MiniIsoLoose", "tree_muon_MiniIsoLoose")
//     .Define("minitree_muon_MiniIsoMedium", "tree_muon_MiniIsoMedium")
//     .Define("minitree_muon_MiniIsoTight", "tree_muon_MiniIsoTight")
//     .Define("minitree_lepton_leadingpt", "tree_lepton_leadingpt")
//     .Define("minitree_lepton_leadingpt2", "tree_lepton_leadingpt2")
//     .Define("minitree_lepton_leadingeta", "tree_lepton_leadingeta")
//     .Define("minitree_lepton_leadingeta2", "tree_lepton_leadingeta2")
//     .Define("minitree_lepton_leadingphi", "tree_lepton_leadingphi")
//     .Define("minitree_lepton_leadingphi2", "tree_lepton_leadingphi2")
//     .Define("minitree_lepton_lepton_dR", "tree_lepton_lepton_dR")
//     .Define("minitree_lepton_lepton_dPhi", "tree_lepton_lepton_dPhi")
//     .Define("minitree_lepton_lepton_dEta", "tree_lepton_lepton_dEta")
//     .Define("minitree_lepton_leadingdxy", "tree_lepton_leadingdxy")
//     .Define("minitree_lepton_leadingdxy2", "tree_lepton_leadingdxy2")
//     .Define("minitree_lepton_leadingdz", "tree_lepton_leadingdz")
//     .Define("minitree_lepton_leadingdz2", "tree_lepton_leadingdz2")
//     .Define("minitree_all_nel", "tree_all_nel")
//     .Define("minitree_electron_nEle", "tree_electron_nEle")
//     .Define("minitree_electron_isPrompt", "tree_electron_isPrompt")
//     .Define("minitree_electron_pt", "tree_electron_pt")
//     .Define("minitree_electron_eta", "tree_electron_eta")
//     .Define("minitree_electron_phi", "tree_electron_phi")
//     .Define("minitree_electron_charge", "tree_electron_charge")
//     .Define("minitree_electron_dxy", "tree_electron_dxy")
//     .Define("minitree_electron_dz", "tree_electron_dz")
//     .Define("minitree_electron_gen", "tree_electron_gen")
//     .Define("minitree_electron_energy", "tree_electron_energy")
//     .Define("minitree_electron_et", "tree_electron_et")
//     .Define("minitree_electron_ecal_trk_postcorr", "tree_electron_ecal_trk_postcorr")
//     .Define("minitree_electron_isoR4", "tree_electron_isoR4")
//     .Define("minitree_electron_IsLoose", "tree_electron_IsLoose")
//     .Define("minitree_electron_IsMedium", "tree_electron_IsMedium")
//     .Define("minitree_electron_IsTight", "tree_electron_IsTight")
//     .Define("minitree_njet", "tree_njet")
//     .Define("minitree_njetNOmu", "tree_njetNOmu")
//     .Define("minitree_jet_pt", "tree_jet_pt")
//     .Define("minitree_jet_eta", "tree_jet_eta")
//     .Define("minitree_jet_phi", "tree_jet_phi")
//     .Define("minitree_jet_HadronFlavour", "tree_jet_HadronFlavour")
//     .Define("minitree_jet_btag_DeepJet", "tree_jet_btag_DeepJet")
//     .Define("minitree_jet_E", "tree_jet_E")
//     .Define("minitree_jet_leadingpt", "tree_jet_leadingpt")
//     .Define("minitree_jet_leadingpt2", "tree_jet_leadingpt2")
//     .Define("minitree_jet_leadingeta", "tree_jet_leadingeta")
//     .Define("minitree_jet_leadingeta2", "tree_jet_leadingeta2")
//     .Define("minitree_jet_jet_dR", "tree_jet_jet_dR")
//     .Define("minitree_jet_jet_dPhi", "tree_jet_jet_dPhi")
//     .Define("minitree_jet_jet_dEta", "tree_jet_jet_dEta")
//     .Define("minitree_muon_jet_dRmin", "tree_muon_jet_dRmin")
//     .Define("minitree_muon_jet_dRmax", "tree_muon_jet_dRmax")
//     .Define("minitree_Evts_MVAval", "tree_Evts_MVAval")
//     .Define("minitree_Evts_MVAvalDY", "tree_Evts_MVAvalDY")
//     .Define("minitree_Evts_MVAvalTT", "tree_Evts_MVAvalTT")

//     .Define("minitree_Hemi", "tree_Hemi")
//     .Define("minitree_Hemi_njet", "tree_Hemi_njet")
//     .Define("minitree_Hemi_njet_nomu", "tree_Hemi_njet_nomu")
//     .Define("minitree_Hemi_pt", "tree_Hemi_pt")
//     .Define("minitree_Hemi_eta", "tree_Hemi_eta")
//     .Define("minitree_Hemi_phi", "tree_Hemi_phi")
//     .Define("minitree_Hemi_nTrks", "tree_Hemi_nTrks")
//     .Define("minitree_Hemi_nTrks_sig", "tree_Hemi_nTrks_sig")
//     .Define("minitree_Hemi_nTrks_bad", "tree_Hemi_nTrks_bad")
//     .Define("minitree_Hemi_mass", "tree_Hemi_mass")
//     .Define("minitree_HemiMu_mass", "tree_HemiMu_mass")
//     .Define("minitree_HemiMu_pt", "tree_HemiMu_pt")
//     .Define("minitree_HemiMu_dR", "tree_HemiMu_dR")
//     .Define("minitree_HemiMuOp_mass", "tree_HemiMuOp_mass")
//     .Define("minitree_HemiMuOp_pt", "tree_HemiMuOp_pt")
//     .Define("minitree_HemiMuOp_dR", "tree_HemiMuOp_dR")
//     .Define("minitree_Hemi_dR12", "tree_Hemi_dR12")
//     .Define("minitree_Hemi_Vtx_step", "tree_Hemi_Vtx_step")
//     .Define("minitree_Hemi_Vtx_isTight", "tree_Hemi_Vtx_isTight")
//     .Define("minitree_Hemi_Vtx_NChi2", "tree_Hemi_Vtx_NChi2")
//     .Define("minitree_Hemi_Vtx_nTrks", "tree_Hemi_Vtx_nTrks")
//     .Define("minitree_Hemi_Vtx_nTrks_sig", "tree_Hemi_Vtx_nTrks_sig")
//     .Define("minitree_Hemi_Vtx_nTrks_bad", "tree_Hemi_Vtx_nTrks_bad")
//     .Define("minitree_Hemi_Vtx_x", "tree_Hemi_Vtx_x")
//     .Define("minitree_Hemi_Vtx_y", "tree_Hemi_Vtx_y")
//     .Define("minitree_Hemi_Vtx_z", "tree_Hemi_Vtx_z")
//     .Define("minitree_Hemi_Vtx_r", "tree_Hemi_Vtx_r")
//     .Define("minitree_Hemi_Vtx_dR", "tree_Hemi_Vtx_dR")
//     .Define("minitree_Hemi_Vtx_SumtrackWeight", "tree_Hemi_Vtx_SumtrackWeight")
//     .Define("minitree_Hemi_Vtx_track_MeanDCA_d", "tree_Hemi_Vtx_track_MeanDCA_d")
//     .Define("minitree_Hemi_Vtx_Mass", "tree_Hemi_Vtx_Mass")
//     .Define("minitree_Hemi_Vtx_dist", "tree_Hemi_Vtx_dist")
//     .Define("minitree_event_nVtx", "tree_event_nVtx")
//     .Define("minitree_event_Vtx_Vtx_dr", "tree_event_Vtx_Vtx_dr")
//     .Define("minitree_event_Vtx_Vtx_dz", "tree_event_Vtx_Vtx_dz")
//     .Define("minitree_event_Vtx_Vtx_dd", "tree_event_Vtx_Vtx_dd")
//     .Define("minitree_event_Vtx_Vtx_reldd", "tree_event_Vtx_Vtx_reldd")
//     .Define("minitree_event_Vtx_Vtx_dR", "tree_event_Vtx_Vtx_dR")
//     .Define("minitree_event_Vtx_Vtx_step", "tree_event_Vtx_Vtx_step")
//     .Define("minitree_Hemi_SecLLP", "tree_Hemi_SecLLP")
//     .Define("minitree_Hemi_LLP_SecVtx_dz", "tree_Hemi_LLP_SecVtx_dz")
//     .Define("minitree_Hemi_LLP_SecVtx_dr", "tree_Hemi_LLP_SecVtx_dr")
//     .Define("minitree_Hemi_SecLLP_ping", "tree_Hemi_SecLLP_ping")
//     .Define("minitree_event_SecLLP_ping", "tree_event_SecLLP_ping")
//     .Define("minitree_Hemi_SecVtx", "tree_Hemi_SecVtx")
//     .Define("minitree_Hemi_SecVtx_step", "tree_Hemi_SecVtx_step")
//     .Define("minitree_Hemi_SecVtx_x", "tree_Hemi_SecVtx_x")
//     .Define("minitree_Hemi_SecVtx_y", "tree_Hemi_SecVtx_y")
//     .Define("minitree_Hemi_SecVtx_z", "tree_Hemi_SecVtx_z")
//     .Define("minitree_Hemi_SecVtx_r", "tree_Hemi_SecVtx_r")
//     .Define("minitree_Hemi_SecVtx_dR", "tree_Hemi_SecVtx_dR")
//     .Define("minitree_Hemi_SecVtx_nTrks", "tree_Hemi_SecVtx_nTrks")
//     .Define("minitree_Hemi_SecVtx_NChi2", "tree_Hemi_SecVtx_NChi2")
//     .Define("minitree_Hemi_SecVtx_dist", "tree_Hemi_SecVtx_dist")
//     .Define("minitree_Hemi_SecVtx_track_MeanDCA_d", "tree_Hemi_SecVtx_track_MeanDCA_d")
//     .Define("minitree_Hemi_SecVtx_SumtrackWeight", "tree_Hemi_SecVtx_SumtrackWeight")
//     .Define("minitree_Hemi_SecVtx_Mass", "tree_Hemi_SecVtx_Mass")
//     .Define("minitree_event_MergedVtx_Vtx_dr", "tree_event_MergedVtx_Vtx_dr")
//     .Define("minitree_event_MergedVtx_Vtx_dz", "tree_event_MergedVtx_Vtx_dz")
//     .Define("minitree_event_MergedVtx_Vtx_dd", "tree_event_MergedVtx_Vtx_dd")
//     .Define("minitree_event_MergedVtx_Vtx_reldd", "tree_event_MergedVtx_Vtx_reldd")
//     .Define("minitree_event_MergedVtx_Vtx_dR", "tree_event_MergedVtx_Vtx_dR")
//     .Define("minitree_event_MergedVtx_Vtx_step", "tree_event_MergedVtx_Vtx_step")
//     .Define("minitree_Hemi_Vtx_BDT_nTrks", "tree_Hemi_Vtx_BDT_nTrks")
//     .Define("minitree_Hemi_Vtx_BDT_NChi2", "tree_Hemi_Vtx_BDT_NChi2")
//     .Define("minitree_Hemi_Vtx_BDT_step", "tree_Hemi_Vtx_BDT_step")
//     .Define("minitree_Hemi_Vtx_BDT_STW", "tree_Hemi_Vtx_BDT_STW")
//     .Define("minitree_Hemi_Vtx_BDT_Mass", "tree_Hemi_Vtx_BDT_Mass")
//     .Define("minitree_Hemi_Vtx_BDT_HMass", "tree_Hemi_Vtx_BDT_HMass")
//     .Define("minitree_Hemi_Vtx_BDT_ntrk10", "tree_Hemi_Vtx_BDT_ntrk10")
//     .Define("minitree_Hemi_Vtx_BDT_ntrk20", "tree_Hemi_Vtx_BDT_ntrk20")
//     .Define("minitree_Hemi_Vtx_BDT_MeanDCA", "tree_Hemi_Vtx_BDT_MeanDCA")
//     .Define("minitree_Hemi_Vtx_MVAval_Loose", "tree_Hemi_Vtx_MVAval_Loose")
//     .Define("minitree_Hemi_Vtx_MVAval_Tight", "tree_Hemi_Vtx_MVAval_Tight")

//    .Define("minitree_track_ipc","tree_track_ipc")
//    .Define("minitree_track_lost","tree_track_lost")
//    .Define("minitree_track_px","tree_track_px")
//    .Define("minitree_track_py","tree_track_py")
//    .Define("minitree_track_pz","tree_track_pz")
//    .Define("minitree_track_pt","tree_track_pt")
//    .Define("minitree_track_eta","tree_track_eta")
//    .Define("minitree_track_phi","tree_track_phi")
//    .Define("minitree_track_charge","tree_track_charge")
//    .Define("minitree_track_NChi2","tree_track_NChi2")
//    .Define("minitree_track_isHighPurity","tree_track_isHighPurity")
//    .Define("minitree_track_dxy","tree_track_dxy")
//    .Define("minitree_track_dxyError","tree_track_dxyError")
//    .Define("minitree_track_drSig","tree_track_drSig")
//    .Define("minitree_track_dz","tree_track_dz")
//    .Define("minitree_track_dzError","tree_track_dzError")
//    .Define("minitree_track_dzSig","tree_track_dzSig")
//    .Define("minitree_track_nHit","tree_track_nHit")
//    .Define("minitree_track_nHitPixel","tree_track_nHitPixel")
//    .Define("minitree_track_nHitTIB","tree_track_nHitTIB")
//    .Define("minitree_track_nHitTID","tree_track_nHitTID")
//    .Define("minitree_track_nHitTOB","tree_track_nHitTOB")
//    .Define("minitree_track_nHitTEC","tree_track_nHitTEC")
//    .Define("minitree_track_nHitPXB","tree_track_nHitPXB")
//    .Define("minitree_track_nHitPXF","tree_track_nHitPXF")
//    .Define("minitree_track_isHitPixel","tree_track_isHitPixel")
//    .Define("minitree_track_nLayers","tree_track_nLayers")
//    .Define("minitree_track_nLayersPixel","tree_track_nLayersPixel")
//    .Define("minitree_track_x","tree_track_x")
//    .Define("minitree_track_y","tree_track_y")
//    .Define("minitree_track_z","tree_track_z")
//    .Define("minitree_track_firstHit","tree_track_firstHit")
//    .Define("minitree_track_region","tree_track_region")
//    .Define("minitree_track_firstHit_x","tree_track_firstHit_x")
//    .Define("minitree_track_firstHit_y","tree_track_firstHit_y")
//    .Define("minitree_track_firstHit_z","tree_track_firstHit_z")
//    .Define("minitree_track_iJet","tree_track_iJet")
//    .Define("minitree_track_ntrk10","tree_track_ntrk10")
//    .Define("minitree_track_ntrk20","tree_track_ntrk20")
//    .Define("minitree_track_ntrk30","tree_track_ntrk30")
//    .Define("minitree_track_ntrk40","tree_track_ntrk40")
//    .Define("minitree_track_MVAval","tree_track_MVAval")

//    .Define("minitree_track_Hemi_dR","tree_track_Hemi_dR")
//    .Define("minitree_track_Hemi_dRmax","tree_track_Hemi_dRmax")

//    .Define("minitree_K0_mass","tree_K0_mass")
//    .Define("minitree_K0_pt","tree_K0_mass")

//    .Define("minitree_L0_mass","tree_L0_mass")
//    .Define("minitree_L0_pt","tree_L0_pt")

//    .Define("minitree_V0_reco_mass","tree_V0_reco_mass")
//    .Define("minitree_V0_reco_pt","tree_V0_reco_pt")
//    .Define("minitree_V0_reco_source","tree_V0_reco_source")

//    .Define("minitree_SecInt_mass","tree_SecInt_mass")
//    .Define("minitree_SecInt_pt","tree_SecInt_pt")
//    .Define("minitree_SecInt_drSig","tree_SecInt_drSig")
//    .Define("minitree_SecInt_dzSig","tree_SecInt_dzSig")
//    .Define("minitree_SecInt_layer","tree_SecInt_layer")
//    .Define("minitree_SecInt_selec","tree_SecInt_selec")
//    .Define("minitree_SecInt_r","tree_SecInt_r")
//    .Define("minitree_SecInt_z","tree_SecInt_z");



    filtered_df.Snapshot("ttree", output_file
    
    ,
     {
        "eventNumber", "tree_LHE_Weights", "tree_MCEvt_weight", "tree_only_gen_wt", 
     "tree_event_weight", "tree_genTop_Weight", "tree_gen_top_pt", 
     "PUweight", "PUweight_Up", "PUweight_Down", 
      "Prefweight_Down","Prefweight","Prefweight_Up",

     "tree_Filter", "tree_FilterSameSign", "tree_trigger_doublelepton", "tree_trigger_singlelepton", 
        "tree_Good_PV", "tree_nPV", 
      "tree_HT", "tree_TRACK_SIZE", "tree_nTracks", "tree_nLostTracks", 
      "tree_all_nmu", "tree_nmu", "tree_LT", "tree_Mmumu", "tree_MmumuSameSign",
      "tree_muon_pt",  "tree_muon_eta", "tree_muon_dxy", "tree_muon_dz", 
        "tree_lepton_leadingpt", "tree_lepton_leadingpt2", "tree_lepton_leadingeta", 
      "tree_lepton_leadingeta2", "tree_lepton_leadingphi", "tree_lepton_leadingphi2", "tree_lepton_lepton_dR", 
      "tree_lepton_lepton_dPhi", "tree_lepton_lepton_dEta", "tree_lepton_leadingdxy", "tree_lepton_leadingdxy2", 
      "tree_lepton_leadingdz", "tree_lepton_leadingdz2", "tree_all_nel", "tree_electron_nEle", "tree_electron_isPrompt", 
      "tree_electron_pt", "tree_electron_eta", 
       "tree_njet", "tree_njetNOmu", "tree_jet_pt", "tree_jet_eta", "tree_jet_leadingpt", "tree_jet_leadingpt2", "tree_jet_leadingeta", 
       "tree_jet_leadingeta2", "tree_jet_jet_dR", "tree_jet_jet_dPhi", "tree_jet_jet_dEta", "tree_muon_jet_dRmin",
        "tree_muon_jet_dRmax", "tree_Hemi", 
        "tree_Hemi_njet", "tree_Hemi_njet_nomu", "tree_Hemi_pt", "tree_Hemi_eta", "tree_Hemi_phi", 
        "tree_Hemi_nTrks",  "tree_Hemi_mass", "tree_HemiMu_mass"
        
        , 

        "tree_HemiMu_pt", "tree_HemiMu_dR", "tree_HemiMuOp_mass", "tree_HemiMuOp_pt", "tree_HemiMuOp_dR", 
        "tree_Hemi_dR12", "tree_Hemi_Vtx_step", "tree_Hemi_Vtx_isTight", "tree_Hemi_Vtx_NChi2", "tree_Hemi_Vtx_nTrks",
        "tree_Hemi_Vtx_r", "tree_Hemi_Vtx_dR", "tree_Hemi_Vtx_SumtrackWeight",
        "tree_Hemi_Vtx_track_MeanDCA_d", "tree_Hemi_Vtx_Mass", "tree_Hemi_Vtx_dist",  "tree_Hemi_SecVtx", 
        "tree_Hemi_SecVtx_step", "tree_Hemi_SecVtx_x", "tree_Hemi_SecVtx_y", "tree_Hemi_SecVtx_z", 
        "tree_Hemi_SecVtx_r", "tree_Hemi_SecVtx_dR", "tree_Hemi_SecVtx_nTrks", "tree_Hemi_SecVtx_NChi2", 
        "tree_Hemi_SecVtx_dist", "tree_Hemi_SecVtx_track_MeanDCA_d", "tree_Hemi_SecVtx_SumtrackWeight", 
        "tree_Hemi_SecVtx_Mass", "tree_Hemi_Vtx_BDT_nTrks", "tree_Hemi_Vtx_BDT_NChi2", 
        "tree_Hemi_Vtx_BDT_step", "tree_Hemi_Vtx_BDT_STW", "tree_Hemi_Vtx_BDT_Mass", "tree_Hemi_Vtx_BDT_HMass",
        "tree_Hemi_Vtx_BDT_ntrk10", "tree_Hemi_Vtx_BDT_ntrk20", "tree_Hemi_Vtx_BDT_MeanDCA", "tree_Hemi_Vtx_MVAval_Loose", "tree_Hemi_Vtx_MVAval_Tight"

        
        "tree_track_ipc", "tree_track_lost", "tree_track_px", "tree_track_py", "tree_track_pz",  "tree_track_pt",
        "tree_track_eta", "tree_track_phi", "tree_track_charge", "tree_track_NChi2","tree_track_isHighPurity", "tree_track_dxy",
        "tree_track_dxyError",  "tree_track_drSig",  "tree_track_dz",  "tree_track_dzError",  "tree_track_dzSig", "tree_track_nHit",
        "tree_track_nHitPixel",  "tree_track_nHitTIB",  "tree_track_nHitTID",  "tree_track_nHitTOB",  "tree_track_nHitTEC",  "tree_track_nHitPXB",
        "tree_track_nHitPXF", "tree_track_isHitPixel",  "tree_track_nLayers",  "tree_track_nLayersPixel",  "tree_track_x",  "tree_track_y",
        "tree_track_z",  "tree_track_firstHit",  "tree_track_region",  "tree_track_firstHit_x",  "tree_track_firstHit_y",  "tree_track_firstHit_z",
        "tree_track_iJet",  "tree_track_ntrk10",  "tree_track_ntrk20", "tree_track_ntrk30", "tree_track_ntrk40", "tree_track_MVAval", "tree_track_Hemi_dR",
        "tree_track_Hemi_dRmax", "tree_K0_mass", "tree_K0_pt","tree_L0_mass", "tree_L0_pt", "tree_V0_reco_mass","tree_V0_reco_pt", "tree_V0_reco_source",
        "tree_SecInt_mass", "tree_SecInt_pt","tree_SecInt_drSig", "tree_SecInt_dzSig", "tree_SecInt_layer", "tree_SecInt_selec","tree_SecInt_r",  "tree_SecInt_z"

           
           }
           
           );

std::cout << "Analyse terminee. TTree sauvegarde dans " << output_file << std::endl;

    return 0;
}
