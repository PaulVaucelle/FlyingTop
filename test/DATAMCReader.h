//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 20 09:45:42 2024 by ROOT version 5.34/36
// from TTree ttree/summary information
// found on file: MiniDATAMC_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root
//////////////////////////////////////////////////////////

#ifndef DATAMCReader_h
#define DATAMCReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
// Header file for the classes stored in the TTree if any.
#include <vector>


// Fixed size dimensions of array or collections stored in the TTree if any.

class DATAMCReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

      // Declaration of leaf types
   vector<double>  *minitree_only_gen_wt;
   vector<double>  *minitree_genTop_Weight;
   vector<double>  *minitree_gen_top_pt;
   vector<double>  *minitree_gen_top_rw_pt;
   vector<double>  *miniPUweight;
   vector<double>  *miniPUweight_Up;
   vector<double>  *miniPUweight_Down;
   vector<double>  *miniPrefweight;
   vector<double>  *miniPrefweight_Up;
   vector<double>  *miniPrefweight_Down;
   vector<bool>    *minitree_Filter;
   vector<bool>    *minitree_FilterSameSign;
   vector<int>     *minitree_nPV;
   vector<bool>    *minitree_trigger_doublelepton;
   vector<bool>    *minitree_trigger_singlelepton;
   vector<bool>    *minitree_Good_PV;
   vector<int>     *minitree_smu_mass;
   vector<int>     *minitree_neu_mass;
   vector<float>   *minitree_neu_ctau;
   vector<float>   *minitree_HT;
   vector<int>     *minitree_TRACK_SIZE;
   vector<int>     *minitree_nTracks;
   vector<int>     *minitree_nLostTracks;
   vector<float>   *minitree_all_nmu;
   vector<float>   *minitree_nmu;
   vector<float>   *minitree_LT;
   vector<float>   *minitree_Mmumu;
   vector<float>   *minitree_MmumuSameSign;
   vector<bool>    *minitree_muon_isPrompt;
   vector<float>   *minitree_muon_pt;
   vector<float>   *minitree_muon_SF;
   vector<float>   *minitree_muon_eta;
   vector<float>   *minitree_muon_phi;
   vector<float>   *minitree_muon_dxy;
   vector<float>   *minitree_muon_dz;
   vector<int>     *minitree_muon_charge;
   vector<float>   *minitree_muon_correction;
   vector<float>   *minitree_muon_dxyError;
   vector<float>   *minitree_muon_dzError;
   vector<bool>    *minitree_muon_isLoose;
   vector<bool>    *minitree_muon_isMedium;
   vector<bool>    *minitree_muon_isTight;
   vector<bool>    *minitree_muon_isGlobal;
   vector<bool>    *minitree_muon_PFIsoVeryLoose;
   vector<bool>    *minitree_muon_PFIsoLoose;
   vector<bool>    *minitree_muon_PFIsoMedium;
   vector<bool>    *minitree_muon_PFIsoTight;
   vector<bool>    *minitree_muon_TkIsoLoose;
   vector<bool>    *minitree_muon_TkIsoTight;
   vector<bool>    *minitree_muon_MiniIsoLoose;
   vector<bool>    *minitree_muon_MiniIsoMedium;
   vector<bool>    *minitree_muon_MiniIsoTight;
   vector<float>   *minitree_lepton_leadingpt;
   vector<float>   *minitree_lepton_leadingpt2;
   vector<float>   *minitree_lepton_leadingeta;
   vector<float>   *minitree_lepton_leadingeta2;
   vector<float>   *minitree_lepton_leadingphi;
   vector<float>   *minitree_lepton_leadingphi2;
   vector<float>   *minitree_lepton_lepton_dR;
   vector<float>   *minitree_lepton_lepton_dPhi;
   vector<float>   *minitree_lepton_lepton_dEta;
   vector<float>   *minitree_lepton_leadingdxy;
   vector<float>   *minitree_lepton_leadingdxy2;
   vector<float>   *minitree_lepton_leadingdz;
   vector<float>   *minitree_lepton_leadingdz2;
   vector<float>   *minitree_reco_lepton_leadingpt;
   vector<float>   *minitree_reco_lepton_leadingpt2;
   vector<float>   *minitree_reco_lepton_leadingeta;
   vector<float>   *minitree_reco_lepton_leadingeta2;
   vector<float>   *minitree_reco_lepton_leadingphi;
   vector<float>   *minitree_reco_lepton_leadingphi2;
   vector<float>   *minitree_lepton_b4trigger_leadingpt;
   vector<float>   *minitree_lepton_b4trigger_leadingpt2;
   vector<int>     *minitree_all_nel;
   vector<int>     *minitree_electron_nEle;
   vector<bool>    *minitree_electron_isPrompt;
   vector<float>   *minitree_electron_pt;
   vector<float>   *minitree_electron_eta;
   vector<float>   *minitree_electron_phi;
   vector<int>     *minitree_electron_charge;
   vector<float>   *minitree_electron_dxy;
   vector<float>   *minitree_electron_dz;
   vector<int>     *minitree_electron_gen;
   vector<float>   *minitree_electron_energy;
   vector<float>   *minitree_electron_et;
   vector<float>   *minitree_electron_ecal_trk_postcorr;
   vector<float>   *minitree_electron_isoR4;
   vector<bool>    *minitree_electron_IsLoose;
   vector<bool>    *minitree_electron_IsMedium;
   vector<bool>    *minitree_electron_IsTight;
   vector<int>     *minitree_njet;
   vector<int>     *minitree_njetNOmu;
   vector<float>   *minitree_jet_pt;
   vector<float>   *minitree_jet_eta;
   vector<float>   *minitree_jet_phi;
   vector<float>   *minitree_jet_HadronFlavour;
   vector<float>   *minitree_jet_btag_DeepJet;
   vector<float>   *minitree_jet_E;
   vector<float>   *minitree_jet_leadingpt;
   vector<float>   *minitree_jet_leadingpt2;
   vector<float>   *minitree_jet_leadingeta;
   vector<float>   *minitree_jet_leadingeta2;
   vector<float>   *minitree_jet_jet_dR;
   vector<float>   *minitree_jet_jet_dPhi;
   vector<float>   *minitree_jet_jet_dEta;
   vector<float>   *minitree_muon_jet_dRmin;
   vector<float>   *minitree_muon_jet_dRmax;
   vector<float>   *minitree_elemu_jet_dRmin;
   vector<float>   *minitree_elemu_jet_dRmax;
   vector<float>   *minitree_ele_jet_dRmin;
   vector<float>   *minitree_ele_jet_dRmax;
   vector<int>     *minitree_Hemi;
   vector<int>     *minitree_Hemi_njet;
   vector<int>     *minitree_Hemi_njet_nomu;
   vector<float>   *minitree_Hemi_pt;
   vector<float>   *minitree_Hemi_eta;
   vector<float>   *minitree_Hemi_phi;
   vector<int>     *minitree_Hemi_nTrks;
   vector<int>     *minitree_Hemi_nTrks_sig;
   vector<int>     *minitree_Hemi_nTrks_bad;
   vector<float>   *minitree_Hemi_mass;
   vector<float>   *minitree_HemiMu_mass;
   vector<float>   *minitree_HemiMu_pt;
   vector<float>   *minitree_HemiMu_dR;
   vector<float>   *minitree_HemiMuOp_mass;
   vector<float>   *minitree_HemiMuOp_pt;
   vector<float>   *minitree_HemiMuOp_dR;
   vector<float>   *minitree_Hemi_dR12;
   vector<float>   *minitree_ll_pt;
   vector<float>   *minitree_ll_eta;
   vector<float>   *minitree_ll_phi;
   vector<float>   *minitree_ll_px;
   vector<float>   *minitree_ll_py;
   vector<float>   *minitree_ll_pz;
   vector<float>   *minitree_ll_energy;
   vector<float>   *minitree_ll_mass;
   vector<unsigned int> *minitree_track_ipc;
   vector<bool>    *minitree_track_lost;
   vector<float>   *minitree_track_px;
   vector<float>   *minitree_track_py;
   vector<float>   *minitree_track_pz;
   vector<float>   *minitree_track_pt;
   vector<float>   *minitree_track_eta;
   vector<float>   *minitree_track_phi;
   vector<int>     *minitree_track_charge;
   vector<float>   *minitree_track_NChi2;
   vector<bool>    *minitree_track_isHighPurity;
   vector<float>   *minitree_track_dxy;
   vector<float>   *minitree_track_dxyError;
   vector<float>   *minitree_track_drSig;
   vector<float>   *minitree_track_dz;
   vector<float>   *minitree_track_dzError;
   vector<float>   *minitree_track_dzSig;
   vector<int>     *minitree_track_nHit;
   vector<int>     *minitree_track_nHitPixel;
   vector<int>     *minitree_track_nHitTIB;
   vector<int>     *minitree_track_nHitTID;
   vector<int>     *minitree_track_nHitTOB;
   vector<int>     *minitree_track_nHitTEC;
   vector<int>     *minitree_track_nHitPXB;
   vector<int>     *minitree_track_nHitPXF;
   vector<int>     *minitree_track_isHitPixel;
   vector<int>     *minitree_track_nLayers;
   vector<int>     *minitree_track_nLayersPixel;
   vector<float>   *minitree_track_x;
   vector<float>   *minitree_track_y;
   vector<float>   *minitree_track_z;
   vector<int>     *minitree_track_firstHit;
   vector<float>   *minitree_track_region;
   vector<float>   *minitree_track_firstHit_x;
   vector<float>   *minitree_track_firstHit_y;
   vector<float>   *minitree_track_firstHit_z;
   vector<int>     *minitree_track_iJet;
   vector<float>   *minitree_track_ntrk10;
   vector<float>   *minitree_track_ntrk20;
   vector<float>   *minitree_track_ntrk30;
   vector<float>   *minitree_track_ntrk40;
   vector<double>  *minitree_track_MVAval;
   vector<double>  *minitree_track_Hemi_dR;
   vector<double>  *minitree_track_Hemi_dRmax;
   vector<float>   *minitree_K0_mass;
   vector<float>   *minitree_K0_pt;
   vector<float>   *minitree_L0_mass;
   vector<float>   *minitree_L0_pt;
   vector<float>   *minitree_V0_reco_mass;
   vector<float>   *minitree_V0_reco_pt;
   vector<int>     *minitree_V0_reco_source;
   vector<float>   *minitree_SecInt_mass;
   vector<float>   *minitree_SecInt_pt;
   vector<float>   *minitree_SecInt_r;
   vector<float>   *minitree_SecInt_z;
   vector<float>   *minitree_SecInt_drSig;
   vector<float>   *minitree_SecInt_dzSig;
   vector<int>     *minitree_SecInt_layer;
   vector<bool>    *minitree_SecInt_selec;
   vector<int>     *minitree_Hemi_Vtx_step;
   vector<float>   *minitree_Hemi_Vtx_NChi2;
   vector<int>     *minitree_Hemi_Vtx_nTrks;
   vector<float>   *minitree_Hemi_Vtx_dR;
   vector<float>   *minitree_Hemi_Vtx_xError;
   vector<float>   *minitree_Hemi_Vtx_yError;
   vector<float>   *minitree_Hemi_Vtx_zError;
   vector<float>   *minitree_Hemi_Vtx_trackWeight;
   vector<float>   *minitree_Hemi_Vtx_SumtrackWeight;
   vector<float>   *minitree_Hemi_Vtx_track_MeanDCA_d;
   vector<float>   *minitree_Hemi_Vtx_Mass;
   vector<float>   *minitree_Hemi_Vtx_dist;
   vector<int>     *minitree_Hemi_Vtx_ntrk10;
   vector<int>     *minitree_Hemi_Vtx_ntrk20;
   vector<int>     *minitree_Hemi_SecVtx_step;
   vector<float>   *minitree_Hemi_SecVtx_NChi2;
   vector<float>   *minitree_Hemi_SecVtx_dR;
   vector<float>   *minitree_Hemi_SecVtx_nTrks;
   vector<float>   *minitree_Hemi_SecVtx_dist;
   vector<float>   *minitree_Hemi_SecVtx_track_MeanDCA_d;
   vector<float>   *minitree_Hemi_SecVtx_SumtrackWeight;
   vector<float>   *minitree_Hemi_SecVtx_trackWeight;
   vector<float>   *minitree_Hemi_SecVtx_Mass;
   vector<float>   *minitree_Hemi_Vtx_BDT_nTrks;
   vector<float>   *minitree_Hemi_Vtx_BDT_NChi2;
   vector<float>   *minitree_Hemi_Vtx_BDT_step;
   vector<float>   *minitree_Hemi_Vtx_BDT_Mass;
   vector<float>   *minitree_Hemi_Vtx_BDT_HMass;

   // List of branches
   TBranch        *b_minitree_only_gen_wt;   //!
   TBranch        *b_minitree_genTop_Weight;   //!
   TBranch        *b_minitree_gen_top_pt;   //!
   TBranch        *b_minitree_gen_top_rw_pt;   //!
   TBranch        *b_miniPUweight;   //!
   TBranch        *b_miniPUweight_Up;   //!
   TBranch        *b_miniPUweight_Down;   //!
   TBranch        *b_miniPrefweight;   //!
   TBranch        *b_miniPrefweight_Up;   //!
   TBranch        *b_miniPrefweight_Down;   //!
   TBranch        *b_minitree_Filter;   //!
   TBranch        *b_minitree_FilterSameSign;   //!
   TBranch        *b_minitree_nPV;   //!
   TBranch        *b_minitree_trigger_doublelepton;   //!
   TBranch        *b_minitree_trigger_singlelepton;   //!
   TBranch        *b_minitree_Good_PV;   //!
   TBranch        *b_minitree_smu_mass;   //!
   TBranch        *b_minitree_neu_mass;   //!
   TBranch        *b_minitree_neu_ctau;   //!
   TBranch        *b_minitree_HT;   //!
   TBranch        *b_minitree_TRACK_SIZE;   //!
   TBranch        *b_minitree_nTracks;   //!
   TBranch        *b_minitree_nLostTracks;   //!
   TBranch        *b_minitree_all_nmu;   //!
   TBranch        *b_minitree_nmu;   //!
   TBranch        *b_minitree_LT;   //!
   TBranch        *b_minitree_Mmumu;   //!
   TBranch        *b_minitree_MmumuSameSign;   //!
   TBranch        *b_minitree_muon_isPrompt;   //!
   TBranch        *b_minitree_muon_pt;   //!
   TBranch        *b_minitree_muon_SF;   //!
   TBranch        *b_minitree_muon_eta;   //!
   TBranch        *b_minitree_muon_phi;   //!
   TBranch        *b_minitree_muon_dxy;   //!
   TBranch        *b_minitree_muon_dz;   //!
   TBranch        *b_minitree_muon_charge;   //!
   TBranch        *b_minitree_muon_correction;   //!
   TBranch        *b_minitree_muon_dxyError;   //!
   TBranch        *b_minitree_muon_dzError;   //!
   TBranch        *b_minitree_muon_isLoose;   //!
   TBranch        *b_minitree_muon_isMedium;   //!
   TBranch        *b_minitree_muon_isTight;   //!
   TBranch        *b_minitree_muon_isGlobal;   //!
   TBranch        *b_minitree_muon_PFIsoVeryLoose;   //!
   TBranch        *b_minitree_muon_PFIsoLoose;   //!
   TBranch        *b_minitree_muon_PFIsoMedium;   //!
   TBranch        *b_minitree_muon_PFIsoTight;   //!
   TBranch        *b_minitree_muon_TkIsoLoose;   //!
   TBranch        *b_minitree_muon_TkIsoTight;   //!
   TBranch        *b_minitree_muon_MiniIsoLoose;   //!
   TBranch        *b_minitree_muon_MiniIsoMedium;   //!
   TBranch        *b_minitree_muon_MiniIsoTight;   //!
   TBranch        *b_minitree_lepton_leadingpt;   //!
   TBranch        *b_minitree_lepton_leadingpt2;   //!
   TBranch        *b_minitree_lepton_leadingeta;   //!
   TBranch        *b_minitree_lepton_leadingeta2;   //!
   TBranch        *b_minitree_lepton_leadingphi;   //!
   TBranch        *b_minitree_lepton_leadingphi2;   //!
   TBranch        *b_minitree_lepton_lepton_dR;   //!
   TBranch        *b_minitree_lepton_lepton_dPhi;   //!
   TBranch        *b_minitree_lepton_lepton_dEta;   //!
   TBranch        *b_minitree_lepton_leadingdxy;   //!
   TBranch        *b_minitree_lepton_leadingdxy2;   //!
   TBranch        *b_minitree_lepton_leadingdz;   //!
   TBranch        *b_minitree_lepton_leadingdz2;   //!
   TBranch        *b_minitree_reco_lepton_leadingpt;   //!
   TBranch        *b_minitree_reco_lepton_leadingpt2;   //!
   TBranch        *b_minitree_reco_lepton_leadingeta;   //!
   TBranch        *b_minitree_reco_lepton_leadingeta2;   //!
   TBranch        *b_minitree_reco_lepton_leadingphi;   //!
   TBranch        *b_minitree_reco_lepton_leadingphi2;   //!
   TBranch        *b_minitree_lepton_b4trigger_leadingpt;   //!
   TBranch        *b_minitree_lepton_b4trigger_leadingpt2;   //!
   TBranch        *b_minitree_all_nel;   //!
   TBranch        *b_minitree_electron_nEle;   //!
   TBranch        *b_minitree_electron_isPrompt;   //!
   TBranch        *b_minitree_electron_pt;   //!
   TBranch        *b_minitree_electron_eta;   //!
   TBranch        *b_minitree_electron_phi;   //!
   TBranch        *b_minitree_electron_charge;   //!
   TBranch        *b_minitree_electron_dxy;   //!
   TBranch        *b_minitree_electron_dz;   //!
   TBranch        *b_minitree_electron_gen;   //!
   TBranch        *b_minitree_electron_energy;   //!
   TBranch        *b_minitree_electron_et;   //!
   TBranch        *b_minitree_electron_ecal_trk_postcorr;   //!
   TBranch        *b_minitree_electron_isoR4;   //!
   TBranch        *b_minitree_electron_IsLoose;   //!
   TBranch        *b_minitree_electron_IsMedium;   //!
   TBranch        *b_minitree_electron_IsTight;   //!
   TBranch        *b_minitree_njet;   //!
   TBranch        *b_minitree_njetNOmu;   //!
   TBranch        *b_minitree_jet_pt;   //!
   TBranch        *b_minitree_jet_eta;   //!
   TBranch        *b_minitree_jet_phi;   //!
   TBranch        *b_minitree_jet_HadronFlavour;   //!
   TBranch        *b_minitree_jet_btag_DeepJet;   //!
   TBranch        *b_minitree_jet_E;   //!
   TBranch        *b_minitree_jet_leadingpt;   //!
   TBranch        *b_minitree_jet_leadingpt2;   //!
   TBranch        *b_minitree_jet_leadingeta;   //!
   TBranch        *b_minitree_jet_leadingeta2;   //!
   TBranch        *b_minitree_jet_jet_dR;   //!
   TBranch        *b_minitree_jet_jet_dPhi;   //!
   TBranch        *b_minitree_jet_jet_dEta;   //!
   TBranch        *b_minitree_muon_jet_dRmin;   //!
   TBranch        *b_minitree_muon_jet_dRmax;   //!
   TBranch        *b_minitree_elemu_jet_dRmin;   //!
   TBranch        *b_minitree_elemu_jet_dRmax;   //!
   TBranch        *b_minitree_ele_jet_dRmin;   //!
   TBranch        *b_minitree_ele_jet_dRmax;   //!
   TBranch        *b_minitree_Hemi;   //!
   TBranch        *b_minitree_Hemi_njet;   //!
   TBranch        *b_minitree_Hemi_njet_nomu;   //!
   TBranch        *b_minitree_Hemi_pt;   //!
   TBranch        *b_minitree_Hemi_eta;   //!
   TBranch        *b_minitree_Hemi_phi;   //!
   TBranch        *b_minitree_Hemi_nTrks;   //!
   TBranch        *b_minitree_Hemi_nTrks_sig;   //!
   TBranch        *b_minitree_Hemi_nTrks_bad;   //!
   TBranch        *b_minitree_Hemi_mass;   //!
   TBranch        *b_minitree_HemiMu_mass;   //!
   TBranch        *b_minitree_HemiMu_pt;   //!
   TBranch        *b_minitree_HemiMu_dR;   //!
   TBranch        *b_minitree_HemiMuOp_mass;   //!
   TBranch        *b_minitree_HemiMuOp_pt;   //!
   TBranch        *b_minitree_HemiMuOp_dR;   //!
   TBranch        *b_minitree_Hemi_dR12;   //!
   TBranch        *b_minitree_ll_pt;   //!
   TBranch        *b_minitree_ll_eta;   //!
   TBranch        *b_minitree_ll_phi;   //!
   TBranch        *b_minitree_ll_px;   //!
   TBranch        *b_minitree_ll_py;   //!
   TBranch        *b_minitree_ll_pz;   //!
   TBranch        *b_minitree_ll_energy;   //!
   TBranch        *b_minitree_ll_mass;   //!
   TBranch        *b_minitree_track_ipc;   //!
   TBranch        *b_minitree_track_lost;   //!
   TBranch        *b_minitree_track_px;   //!
   TBranch        *b_minitree_track_py;   //!
   TBranch        *b_minitree_track_pz;   //!
   TBranch        *b_minitree_track_pt;   //!
   TBranch        *b_minitree_track_eta;   //!
   TBranch        *b_minitree_track_phi;   //!
   TBranch        *b_minitree_track_charge;   //!
   TBranch        *b_minitree_track_NChi2;   //!
   TBranch        *b_minitree_track_isHighPurity;   //!
   TBranch        *b_minitree_track_dxy;   //!
   TBranch        *b_minitree_track_dxyError;   //!
   TBranch        *b_minitree_track_drSig;   //!
   TBranch        *b_minitree_track_dz;   //!
   TBranch        *b_minitree_track_dzError;   //!
   TBranch        *b_minitree_track_dzSig;   //!
   TBranch        *b_minitree_track_nHit;   //!
   TBranch        *b_minitree_track_nHitPixel;   //!
   TBranch        *b_minitree_track_nHitTIB;   //!
   TBranch        *b_minitree_track_nHitTID;   //!
   TBranch        *b_minitree_track_nHitTOB;   //!
   TBranch        *b_minitree_track_nHitTEC;   //!
   TBranch        *b_minitree_track_nHitPXB;   //!
   TBranch        *b_minitree_track_nHitPXF;   //!
   TBranch        *b_minitree_track_isHitPixel;   //!
   TBranch        *b_minitree_track_nLayers;   //!
   TBranch        *b_minitree_track_nLayersPixel;   //!
   TBranch        *b_minitree_track_x;   //!
   TBranch        *b_minitree_track_y;   //!
   TBranch        *b_minitree_track_z;   //!
   TBranch        *b_minitree_track_firstHit;   //!
   TBranch        *b_minitree_track_region;   //!
   TBranch        *b_minitree_track_firstHit_x;   //!
   TBranch        *b_minitree_track_firstHit_y;   //!
   TBranch        *b_minitree_track_firstHit_z;   //!
   TBranch        *b_minitree_track_iJet;   //!
   TBranch        *b_minitree_track_ntrk10;   //!
   TBranch        *b_minitree_track_ntrk20;   //!
   TBranch        *b_minitree_track_ntrk30;   //!
   TBranch        *b_minitree_track_ntrk40;   //!
   TBranch        *b_minitree_track_MVAval;   //!
   TBranch        *b_minitree_track_Hemi_dR;   //!
   TBranch        *b_minitree_track_Hemi_dRmax;   //!
   TBranch        *b_minitree_K0_mass;   //!
   TBranch        *b_minitree_K0_pt;   //!
   TBranch        *b_minitree_L0_mass;   //!
   TBranch        *b_minitree_L0_pt;   //!
   TBranch        *b_minitree_V0_reco_mass;   //!
   TBranch        *b_minitree_V0_reco_pt;   //!
   TBranch        *b_minitree_V0_reco_source;   //!
   TBranch        *b_minitree_SecInt_mass;   //!
   TBranch        *b_minitree_SecInt_pt;   //!
   TBranch        *b_minitree_SecInt_r;
   TBranch        *b_minitree_SecInt_z;
   TBranch        *b_minitree_SecInt_drSig;   //!
   TBranch        *b_minitree_SecInt_dzSig;   //!
   TBranch        *b_minitree_SecInt_layer;   //!
   TBranch        *b_minitree_SecInt_selec;   //!
   TBranch        *b_minitree_Hemi_Vtx_step;   //!
   TBranch        *b_minitree_Hemi_Vtx_NChi2;   //!
   TBranch        *b_minitree_Hemi_Vtx_nTrks;   //!
   TBranch        *b_minitree_Hemi_Vtx_dR;   //!
   TBranch        *b_minitree_Hemi_Vtx_xError;   //!
   TBranch        *b_minitree_Hemi_Vtx_yError;   //!
   TBranch        *b_minitree_Hemi_Vtx_zError;   //!
   TBranch        *b_minitree_Hemi_Vtx_trackWeight;   //!
   TBranch        *b_minitree_Hemi_Vtx_SumtrackWeight;   //!
   TBranch        *b_minitree_Hemi_Vtx_track_MeanDCA_d;   //!
   TBranch        *b_minitree_Hemi_Vtx_Mass;   //!
   TBranch        *b_minitree_Hemi_Vtx_dist;   //!
   TBranch        *b_minitree_Hemi_Vtx_ntrk10;   //!
   TBranch        *b_minitree_Hemi_Vtx_ntrk20;   //!
   TBranch        *b_minitree_Hemi_SecVtx_step;   //!
   TBranch        *b_minitree_Hemi_SecVtx_NChi2;   //!
   TBranch        *b_minitree_Hemi_SecVtx_dR;   //!
   TBranch        *b_minitree_Hemi_SecVtx_nTrks;   //!
   TBranch        *b_minitree_Hemi_SecVtx_dist;   //!
   TBranch        *b_minitree_Hemi_SecVtx_track_MeanDCA_d;   //!
   TBranch        *b_minitree_Hemi_SecVtx_SumtrackWeight;   //!
   TBranch        *b_minitree_Hemi_SecVtx_trackWeight;   //!
   TBranch        *b_minitree_Hemi_SecVtx_Mass;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_nTrks;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_NChi2;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_step;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_Mass;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_HMass;   //!
   DATAMCReader(TTree *tree=0, TString sample="", std::vector<TString> thesystlist = std::vector<TString>());
   virtual ~DATAMCReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(bool isMC, TString Prod, TString sample, bool Signal, int Year, bool IsPostAPV, float MeanGenW , int Channel,bool DoubleMuon,  std::vector<TString> thesystlist);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

         //add functions to create and fill histograms
 
   void initializeHisto(TString sample, bool isfirstset);
   void addHisto( TString var, TString selstep, TString sample, int nbins, float min, float max);
   void fillHisto(TString var, TString selstep, TString sample, float val, float weight);
   void addHisto2D(TString var, TString selstep, TString sample, int nxbins, float xmin, float xmax, int nybins, float ymin, float ymax);
   void fillHisto2D( TString var, TString selstep, TString sample, float xval,float yval, float weight);
   virtual float    MeanGenWeight(TString thesample, TString Prod);

   std::vector<TH1F*> histo_list_;
   std::vector<TH2F*> histo_list_2D_;
   std::map<std::string, int> histo_map_;
   std::map<std::string, int> histo_map_2D_;


   int numb_histo;
   int numb_histo_2D_;
   void deleteHisto();


   std::vector<TString> systlist;
};

#endif

#ifdef DATAMCReader_cxx
DATAMCReader::DATAMCReader(TTree *tree, TString sample, std::vector<TString> thesystlist) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(("/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+sample+".root").Data());
      if (!f || !f->IsOpen()) {
         f = new TFile(("/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+sample+".root").Data());
      }
      TDirectory * dir = (TDirectory*)f->Get("/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Macro_new/Ntuple_03_06_24/2018/"+sample+".root:/FlyingTop");
      dir->GetObject("ttree",tree);
      systlist = thesystlist; 
   }
systlist = thesystlist;
   Init(tree);
}

DATAMCReader::~DATAMCReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DATAMCReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DATAMCReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DATAMCReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

 // Set object pointer
   minitree_only_gen_wt = 0;
   minitree_genTop_Weight = 0;
   minitree_gen_top_pt = 0;
   minitree_gen_top_rw_pt = 0;
   miniPUweight = 0;
   miniPUweight_Up = 0;
   miniPUweight_Down = 0;
   miniPrefweight = 0;
   miniPrefweight_Up = 0 ;
   miniPrefweight_Down = 0;
   minitree_Filter = 0;
   minitree_FilterSameSign = 0;
   minitree_nPV = 0;
   minitree_trigger_doublelepton = 0;
   minitree_trigger_singlelepton = 0;
   minitree_Good_PV = 0;
   minitree_smu_mass = 0;
   minitree_neu_mass = 0;
   minitree_neu_ctau = 0;
   minitree_HT = 0;
   minitree_TRACK_SIZE = 0;
   minitree_nTracks = 0;
   minitree_nLostTracks = 0;
   minitree_all_nmu = 0;
   minitree_nmu = 0;
   minitree_LT = 0;
   minitree_Mmumu = 0;
   minitree_MmumuSameSign = 0;
   minitree_muon_isPrompt = 0;
   minitree_muon_pt = 0;
   minitree_muon_SF = 0;
   minitree_muon_eta = 0;
   minitree_muon_phi = 0;
   minitree_muon_dxy = 0;
   minitree_muon_dz = 0;
   minitree_muon_charge = 0;
   minitree_muon_correction = 0;
   minitree_muon_dxyError = 0;
   minitree_muon_dzError = 0;
   minitree_muon_isLoose = 0;
   minitree_muon_isMedium = 0;
   minitree_muon_isTight = 0;
   minitree_muon_isGlobal = 0;
   minitree_muon_PFIsoVeryLoose = 0;
   minitree_muon_PFIsoLoose = 0;
   minitree_muon_PFIsoMedium = 0;
   minitree_muon_PFIsoTight = 0;
   minitree_muon_TkIsoLoose = 0;
   minitree_muon_TkIsoTight = 0;
   minitree_muon_MiniIsoLoose = 0;
   minitree_muon_MiniIsoMedium = 0;
   minitree_muon_MiniIsoTight = 0;
   minitree_lepton_leadingpt = 0;
   minitree_lepton_leadingpt2 = 0;
   minitree_lepton_leadingeta = 0;
   minitree_lepton_leadingeta2 = 0;
   minitree_lepton_leadingphi = 0;
   minitree_lepton_leadingphi2 = 0;
   minitree_lepton_lepton_dR = 0;
   minitree_lepton_lepton_dPhi = 0;
   minitree_lepton_lepton_dEta = 0;
   minitree_lepton_leadingdxy = 0;
   minitree_lepton_leadingdxy2 = 0;
   minitree_lepton_leadingdz = 0;
   minitree_lepton_leadingdz2 = 0;
   minitree_reco_lepton_leadingpt = 0;
   minitree_reco_lepton_leadingpt2 = 0;
   minitree_reco_lepton_leadingeta = 0;
   minitree_reco_lepton_leadingeta2 = 0;
   minitree_reco_lepton_leadingphi = 0;
   minitree_reco_lepton_leadingphi2 = 0;
   minitree_lepton_b4trigger_leadingpt = 0;
   minitree_lepton_b4trigger_leadingpt2 = 0;
   minitree_all_nel = 0;
   minitree_electron_nEle = 0;
   minitree_electron_isPrompt = 0;
   minitree_electron_pt = 0;
   minitree_electron_eta = 0;
   minitree_electron_phi = 0;
   minitree_electron_charge = 0;
   minitree_electron_dxy = 0;
   minitree_electron_dz = 0;
   minitree_electron_gen = 0;
   minitree_electron_energy = 0;
   minitree_electron_et = 0;
   minitree_electron_ecal_trk_postcorr = 0;
   minitree_electron_isoR4 = 0;
   minitree_electron_IsLoose = 0;
   minitree_electron_IsMedium = 0;
   minitree_electron_IsTight = 0;
   minitree_njet = 0;
   minitree_njetNOmu = 0;
   minitree_jet_pt = 0;
   minitree_jet_eta = 0;
   minitree_jet_phi = 0;
   minitree_jet_HadronFlavour = 0;
   minitree_jet_btag_DeepJet = 0;
   minitree_jet_E = 0;
   minitree_jet_leadingpt = 0;
   minitree_jet_leadingpt2 = 0;
   minitree_jet_leadingeta = 0;
   minitree_jet_leadingeta2 = 0;
   minitree_jet_jet_dR = 0;
   minitree_jet_jet_dPhi = 0;
   minitree_jet_jet_dEta = 0;
   minitree_muon_jet_dRmin = 0;
   minitree_muon_jet_dRmax = 0;
   minitree_elemu_jet_dRmin = 0;
   minitree_elemu_jet_dRmax = 0;
   minitree_ele_jet_dRmin = 0;
   minitree_ele_jet_dRmax = 0;
   minitree_Hemi = 0;
   minitree_Hemi_njet = 0;
   minitree_Hemi_njet_nomu = 0;
   minitree_Hemi_pt = 0;
   minitree_Hemi_eta = 0;
   minitree_Hemi_phi = 0;
   minitree_Hemi_nTrks = 0;
   minitree_Hemi_nTrks_sig = 0;
   minitree_Hemi_nTrks_bad = 0;
   minitree_Hemi_mass = 0;
   minitree_HemiMu_mass = 0;
   minitree_HemiMu_pt = 0;
   minitree_HemiMu_dR = 0;
   minitree_HemiMuOp_mass = 0;
   minitree_HemiMuOp_pt = 0;
   minitree_HemiMuOp_dR = 0;
   minitree_Hemi_dR12 = 0;
   minitree_ll_pt = 0;
   minitree_ll_eta = 0;
   minitree_ll_phi = 0;
   minitree_ll_px = 0;
   minitree_ll_py = 0;
   minitree_ll_pz = 0;
   minitree_ll_energy = 0;
   minitree_ll_mass = 0;
   minitree_track_ipc = 0;
   minitree_track_lost = 0;
   minitree_track_px = 0;
   minitree_track_py = 0;
   minitree_track_pz = 0;
   minitree_track_pt = 0;
   minitree_track_eta = 0;
   minitree_track_phi = 0;
   minitree_track_charge = 0;
   minitree_track_NChi2 = 0;
   minitree_track_isHighPurity = 0;
   minitree_track_dxy = 0;
   minitree_track_dxyError = 0;
   minitree_track_drSig = 0;
   minitree_track_dz = 0;
   minitree_track_dzError = 0;
   minitree_track_dzSig = 0;
   minitree_track_nHit = 0;
   minitree_track_nHitPixel = 0;
   minitree_track_nHitTIB = 0;
   minitree_track_nHitTID = 0;
   minitree_track_nHitTOB = 0;
   minitree_track_nHitTEC = 0;
   minitree_track_nHitPXB = 0;
   minitree_track_nHitPXF = 0;
   minitree_track_isHitPixel = 0;
   minitree_track_nLayers = 0;
   minitree_track_nLayersPixel = 0;
   minitree_track_x = 0;
   minitree_track_y = 0;
   minitree_track_z = 0;
   minitree_track_firstHit = 0;
   minitree_track_region = 0;
   minitree_track_firstHit_x = 0;
   minitree_track_firstHit_y = 0;
   minitree_track_firstHit_z = 0;
   minitree_track_iJet = 0;
   minitree_track_ntrk10 = 0;
   minitree_track_ntrk20 = 0;
   minitree_track_ntrk30 = 0;
   minitree_track_ntrk40 = 0;
   minitree_track_MVAval = 0;
   minitree_track_Hemi_dR = 0;
   minitree_track_Hemi_dRmax = 0;
   minitree_K0_mass = 0;
   minitree_K0_pt = 0;
   minitree_L0_mass = 0;
   minitree_L0_pt = 0;
   minitree_V0_reco_mass = 0;
   minitree_V0_reco_pt = 0;
   minitree_V0_reco_source = 0;
   minitree_SecInt_mass = 0;
   minitree_SecInt_pt = 0;
   minitree_SecInt_r = 0;
   minitree_SecInt_z = 0;
   minitree_SecInt_drSig = 0;
   minitree_SecInt_dzSig = 0;
   minitree_SecInt_layer = 0;
   minitree_SecInt_selec = 0;
   minitree_Hemi_Vtx_step = 0;
   minitree_Hemi_Vtx_NChi2 = 0;
   minitree_Hemi_Vtx_nTrks = 0;
   minitree_Hemi_Vtx_dR = 0;
   minitree_Hemi_Vtx_xError = 0;
   minitree_Hemi_Vtx_yError = 0;
   minitree_Hemi_Vtx_zError = 0;
   minitree_Hemi_Vtx_trackWeight = 0;
   minitree_Hemi_Vtx_SumtrackWeight = 0;
   minitree_Hemi_Vtx_track_MeanDCA_d = 0;
   minitree_Hemi_Vtx_Mass = 0;
   minitree_Hemi_Vtx_dist = 0;
   minitree_Hemi_Vtx_ntrk10 = 0;
   minitree_Hemi_Vtx_ntrk20 = 0;
   minitree_Hemi_SecVtx_step = 0;
   minitree_Hemi_SecVtx_NChi2 = 0;
   minitree_Hemi_SecVtx_dR = 0;
   minitree_Hemi_SecVtx_nTrks = 0;
   minitree_Hemi_SecVtx_dist = 0;
   minitree_Hemi_SecVtx_track_MeanDCA_d = 0;
   minitree_Hemi_SecVtx_SumtrackWeight = 0;
   minitree_Hemi_SecVtx_trackWeight = 0;
   minitree_Hemi_SecVtx_Mass = 0;
   minitree_Hemi_Vtx_BDT_nTrks = 0;
   minitree_Hemi_Vtx_BDT_NChi2 = 0;
   minitree_Hemi_Vtx_BDT_step = 0;
   minitree_Hemi_Vtx_BDT_Mass = 0;
   minitree_Hemi_Vtx_BDT_HMass = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("minitree_only_gen_wt", &minitree_only_gen_wt, &b_minitree_only_gen_wt);
   fChain->SetBranchAddress("minitree_genTop_Weight", &minitree_genTop_Weight, &b_minitree_genTop_Weight);
   fChain->SetBranchAddress("minitree_gen_top_pt", &minitree_gen_top_pt, &b_minitree_gen_top_pt);
   fChain->SetBranchAddress("minitree_gen_top_rw_pt", &minitree_gen_top_rw_pt, &b_minitree_gen_top_rw_pt);
   fChain->SetBranchAddress("miniPUweight", &miniPUweight, &b_miniPUweight);
   fChain->SetBranchAddress("miniPUweight_Up", &miniPUweight_Up, &b_miniPUweight_Up);
   fChain->SetBranchAddress("miniPUweight_Down", &miniPUweight_Down, &b_miniPUweight_Down);
   fChain->SetBranchAddress("miniPrefweight", &miniPrefweight, &b_miniPrefweight);
   fChain->SetBranchAddress("miniPrefweight_Up",&miniPrefweight_Up,&b_miniPrefweight_Up);
   fChain->SetBranchAddress("miniPrefweight_Down",&miniPrefweight_Down,&b_miniPrefweight_Down);
   fChain->SetBranchAddress("minitree_Filter", &minitree_Filter, &b_minitree_Filter);
   fChain->SetBranchAddress("minitree_FilterSameSign", &minitree_FilterSameSign, &b_minitree_FilterSameSign);
   fChain->SetBranchAddress("minitree_nPV", &minitree_nPV, &b_minitree_nPV);
   fChain->SetBranchAddress("minitree_trigger_doublelepton", &minitree_trigger_doublelepton, &b_minitree_trigger_doublelepton);
   fChain->SetBranchAddress("minitree_trigger_singlelepton", &minitree_trigger_singlelepton, &b_minitree_trigger_singlelepton);
   fChain->SetBranchAddress("minitree_Good_PV", &minitree_Good_PV, &b_minitree_Good_PV);
   fChain->SetBranchAddress("minitree_smu_mass", &minitree_smu_mass, &b_minitree_smu_mass);
   fChain->SetBranchAddress("minitree_neu_mass", &minitree_neu_mass, &b_minitree_neu_mass);
   fChain->SetBranchAddress("minitree_neu_ctau", &minitree_neu_ctau, &b_minitree_neu_ctau);
   fChain->SetBranchAddress("minitree_HT", &minitree_HT, &b_minitree_HT);
   fChain->SetBranchAddress("minitree_TRACK_SIZE", &minitree_TRACK_SIZE, &b_minitree_TRACK_SIZE);
   fChain->SetBranchAddress("minitree_nTracks", &minitree_nTracks, &b_minitree_nTracks);
   fChain->SetBranchAddress("minitree_nLostTracks", &minitree_nLostTracks, &b_minitree_nLostTracks);
   fChain->SetBranchAddress("minitree_all_nmu", &minitree_all_nmu, &b_minitree_all_nmu);
   fChain->SetBranchAddress("minitree_nmu", &minitree_nmu, &b_minitree_nmu);
   fChain->SetBranchAddress("minitree_LT", &minitree_LT, &b_minitree_LT);
   fChain->SetBranchAddress("minitree_Mmumu", &minitree_Mmumu, &b_minitree_Mmumu);
   fChain->SetBranchAddress("minitree_MmumuSameSign", &minitree_MmumuSameSign, &b_minitree_MmumuSameSign);
   fChain->SetBranchAddress("minitree_muon_isPrompt", &minitree_muon_isPrompt, &b_minitree_muon_isPrompt);
   fChain->SetBranchAddress("minitree_muon_pt", &minitree_muon_pt, &b_minitree_muon_pt);
   fChain->SetBranchAddress("minitree_muon_SF", &minitree_muon_SF, &b_minitree_muon_SF);
   fChain->SetBranchAddress("minitree_muon_eta", &minitree_muon_eta, &b_minitree_muon_eta);
   fChain->SetBranchAddress("minitree_muon_phi", &minitree_muon_phi, &b_minitree_muon_phi);
   fChain->SetBranchAddress("minitree_muon_dxy", &minitree_muon_dxy, &b_minitree_muon_dxy);
   fChain->SetBranchAddress("minitree_muon_dz", &minitree_muon_dz, &b_minitree_muon_dz);
   fChain->SetBranchAddress("minitree_muon_charge", &minitree_muon_charge, &b_minitree_muon_charge);
   fChain->SetBranchAddress("minitree_muon_correction", &minitree_muon_correction, &b_minitree_muon_correction);
   fChain->SetBranchAddress("minitree_muon_dxyError", &minitree_muon_dxyError, &b_minitree_muon_dxyError);
   fChain->SetBranchAddress("minitree_muon_dzError", &minitree_muon_dzError, &b_minitree_muon_dzError);
   fChain->SetBranchAddress("minitree_muon_isLoose", &minitree_muon_isLoose, &b_minitree_muon_isLoose);
   fChain->SetBranchAddress("minitree_muon_isMedium", &minitree_muon_isMedium, &b_minitree_muon_isMedium);
   fChain->SetBranchAddress("minitree_muon_isTight", &minitree_muon_isTight, &b_minitree_muon_isTight);
   fChain->SetBranchAddress("minitree_muon_isGlobal", &minitree_muon_isGlobal, &b_minitree_muon_isGlobal);
   fChain->SetBranchAddress("minitree_muon_PFIsoVeryLoose", &minitree_muon_PFIsoVeryLoose, &b_minitree_muon_PFIsoVeryLoose);
   fChain->SetBranchAddress("minitree_muon_PFIsoLoose", &minitree_muon_PFIsoLoose, &b_minitree_muon_PFIsoLoose);
   fChain->SetBranchAddress("minitree_muon_PFIsoMedium", &minitree_muon_PFIsoMedium, &b_minitree_muon_PFIsoMedium);
   fChain->SetBranchAddress("minitree_muon_PFIsoTight", &minitree_muon_PFIsoTight, &b_minitree_muon_PFIsoTight);
   fChain->SetBranchAddress("minitree_muon_TkIsoLoose", &minitree_muon_TkIsoLoose, &b_minitree_muon_TkIsoLoose);
   fChain->SetBranchAddress("minitree_muon_TkIsoTight", &minitree_muon_TkIsoTight, &b_minitree_muon_TkIsoTight);
   fChain->SetBranchAddress("minitree_muon_MiniIsoLoose", &minitree_muon_MiniIsoLoose, &b_minitree_muon_MiniIsoLoose);
   fChain->SetBranchAddress("minitree_muon_MiniIsoMedium", &minitree_muon_MiniIsoMedium, &b_minitree_muon_MiniIsoMedium);
   fChain->SetBranchAddress("minitree_muon_MiniIsoTight", &minitree_muon_MiniIsoTight, &b_minitree_muon_MiniIsoTight);
   fChain->SetBranchAddress("minitree_lepton_leadingpt", &minitree_lepton_leadingpt, &b_minitree_lepton_leadingpt);
   fChain->SetBranchAddress("minitree_lepton_leadingpt2", &minitree_lepton_leadingpt2, &b_minitree_lepton_leadingpt2);
   fChain->SetBranchAddress("minitree_lepton_leadingeta", &minitree_lepton_leadingeta, &b_minitree_lepton_leadingeta);
   fChain->SetBranchAddress("minitree_lepton_leadingeta2", &minitree_lepton_leadingeta2, &b_minitree_lepton_leadingeta2);
   fChain->SetBranchAddress("minitree_lepton_leadingphi", &minitree_lepton_leadingphi, &b_minitree_lepton_leadingphi);
   fChain->SetBranchAddress("minitree_lepton_leadingphi2", &minitree_lepton_leadingphi2, &b_minitree_lepton_leadingphi2);
   fChain->SetBranchAddress("minitree_lepton_lepton_dR", &minitree_lepton_lepton_dR, &b_minitree_lepton_lepton_dR);
   fChain->SetBranchAddress("minitree_lepton_lepton_dPhi", &minitree_lepton_lepton_dPhi, &b_minitree_lepton_lepton_dPhi);
   fChain->SetBranchAddress("minitree_lepton_lepton_dEta", &minitree_lepton_lepton_dEta, &b_minitree_lepton_lepton_dEta);
   fChain->SetBranchAddress("minitree_lepton_leadingdxy", &minitree_lepton_leadingdxy, &b_minitree_lepton_leadingdxy);
   fChain->SetBranchAddress("minitree_lepton_leadingdxy2", &minitree_lepton_leadingdxy2, &b_minitree_lepton_leadingdxy2);
   fChain->SetBranchAddress("minitree_lepton_leadingdz", &minitree_lepton_leadingdz, &b_minitree_lepton_leadingdz);
   fChain->SetBranchAddress("minitree_lepton_leadingdz2", &minitree_lepton_leadingdz2, &b_minitree_lepton_leadingdz2);
   fChain->SetBranchAddress("minitree_reco_lepton_leadingpt", &minitree_reco_lepton_leadingpt, &b_minitree_reco_lepton_leadingpt);
   fChain->SetBranchAddress("minitree_reco_lepton_leadingpt2", &minitree_reco_lepton_leadingpt2, &b_minitree_reco_lepton_leadingpt2);
   fChain->SetBranchAddress("minitree_reco_lepton_leadingeta", &minitree_reco_lepton_leadingeta, &b_minitree_reco_lepton_leadingeta);
   fChain->SetBranchAddress("minitree_reco_lepton_leadingeta2", &minitree_reco_lepton_leadingeta2, &b_minitree_reco_lepton_leadingeta2);
   fChain->SetBranchAddress("minitree_reco_lepton_leadingphi", &minitree_reco_lepton_leadingphi, &b_minitree_reco_lepton_leadingphi);
   fChain->SetBranchAddress("minitree_reco_lepton_leadingphi2", &minitree_reco_lepton_leadingphi2, &b_minitree_reco_lepton_leadingphi2);
   fChain->SetBranchAddress("minitree_lepton_b4trigger_leadingpt", &minitree_lepton_b4trigger_leadingpt, &b_minitree_lepton_b4trigger_leadingpt);
   fChain->SetBranchAddress("minitree_lepton_b4trigger_leadingpt2", &minitree_lepton_b4trigger_leadingpt2, &b_minitree_lepton_b4trigger_leadingpt2);
   fChain->SetBranchAddress("minitree_all_nel", &minitree_all_nel, &b_minitree_all_nel);
   fChain->SetBranchAddress("minitree_electron_nEle", &minitree_electron_nEle, &b_minitree_electron_nEle);
   fChain->SetBranchAddress("minitree_electron_isPrompt", &minitree_electron_isPrompt, &b_minitree_electron_isPrompt);
   fChain->SetBranchAddress("minitree_electron_pt", &minitree_electron_pt, &b_minitree_electron_pt);
   fChain->SetBranchAddress("minitree_electron_eta", &minitree_electron_eta, &b_minitree_electron_eta);
   fChain->SetBranchAddress("minitree_electron_phi", &minitree_electron_phi, &b_minitree_electron_phi);
   fChain->SetBranchAddress("minitree_electron_charge", &minitree_electron_charge, &b_minitree_electron_charge);
   fChain->SetBranchAddress("minitree_electron_dxy", &minitree_electron_dxy, &b_minitree_electron_dxy);
   fChain->SetBranchAddress("minitree_electron_dz", &minitree_electron_dz, &b_minitree_electron_dz);
   fChain->SetBranchAddress("minitree_electron_gen", &minitree_electron_gen, &b_minitree_electron_gen);
   fChain->SetBranchAddress("minitree_electron_energy", &minitree_electron_energy, &b_minitree_electron_energy);
   fChain->SetBranchAddress("minitree_electron_et", &minitree_electron_et, &b_minitree_electron_et);
   fChain->SetBranchAddress("minitree_electron_ecal_trk_postcorr", &minitree_electron_ecal_trk_postcorr, &b_minitree_electron_ecal_trk_postcorr);
   fChain->SetBranchAddress("minitree_electron_isoR4", &minitree_electron_isoR4, &b_minitree_electron_isoR4);
   fChain->SetBranchAddress("minitree_electron_IsLoose", &minitree_electron_IsLoose, &b_minitree_electron_IsLoose);
   fChain->SetBranchAddress("minitree_electron_IsMedium", &minitree_electron_IsMedium, &b_minitree_electron_IsMedium);
   fChain->SetBranchAddress("minitree_electron_IsTight", &minitree_electron_IsTight, &b_minitree_electron_IsTight);
   fChain->SetBranchAddress("minitree_njet", &minitree_njet, &b_minitree_njet);
   fChain->SetBranchAddress("minitree_njetNOmu", &minitree_njetNOmu, &b_minitree_njetNOmu);
   fChain->SetBranchAddress("minitree_jet_pt", &minitree_jet_pt, &b_minitree_jet_pt);
   fChain->SetBranchAddress("minitree_jet_eta", &minitree_jet_eta, &b_minitree_jet_eta);
   fChain->SetBranchAddress("minitree_jet_phi", &minitree_jet_phi, &b_minitree_jet_phi);
   fChain->SetBranchAddress("minitree_jet_HadronFlavour", &minitree_jet_HadronFlavour, &b_minitree_jet_HadronFlavour);
   fChain->SetBranchAddress("minitree_jet_btag_DeepJet", &minitree_jet_btag_DeepJet, &b_minitree_jet_btag_DeepJet);
   fChain->SetBranchAddress("minitree_jet_E", &minitree_jet_E, &b_minitree_jet_E);
   fChain->SetBranchAddress("minitree_jet_leadingpt", &minitree_jet_leadingpt, &b_minitree_jet_leadingpt);
   fChain->SetBranchAddress("minitree_jet_leadingpt2", &minitree_jet_leadingpt2, &b_minitree_jet_leadingpt2);
   fChain->SetBranchAddress("minitree_jet_leadingeta", &minitree_jet_leadingeta, &b_minitree_jet_leadingeta);
   fChain->SetBranchAddress("minitree_jet_leadingeta2", &minitree_jet_leadingeta2, &b_minitree_jet_leadingeta2);
   fChain->SetBranchAddress("minitree_jet_jet_dR", &minitree_jet_jet_dR, &b_minitree_jet_jet_dR);
   fChain->SetBranchAddress("minitree_jet_jet_dPhi", &minitree_jet_jet_dPhi, &b_minitree_jet_jet_dPhi);
   fChain->SetBranchAddress("minitree_jet_jet_dEta", &minitree_jet_jet_dEta, &b_minitree_jet_jet_dEta);
   fChain->SetBranchAddress("minitree_muon_jet_dRmin", &minitree_muon_jet_dRmin, &b_minitree_muon_jet_dRmin);
   fChain->SetBranchAddress("minitree_muon_jet_dRmax", &minitree_muon_jet_dRmax, &b_minitree_muon_jet_dRmax);
   fChain->SetBranchAddress("minitree_elemu_jet_dRmin", &minitree_elemu_jet_dRmin, &b_minitree_elemu_jet_dRmin);
   fChain->SetBranchAddress("minitree_elemu_jet_dRmax", &minitree_elemu_jet_dRmax, &b_minitree_elemu_jet_dRmax);
   fChain->SetBranchAddress("minitree_ele_jet_dRmin", &minitree_ele_jet_dRmin, &b_minitree_ele_jet_dRmin);
   fChain->SetBranchAddress("minitree_ele_jet_dRmax", &minitree_ele_jet_dRmax, &b_minitree_ele_jet_dRmax);
   fChain->SetBranchAddress("minitree_Hemi", &minitree_Hemi, &b_minitree_Hemi);
   fChain->SetBranchAddress("minitree_Hemi_njet", &minitree_Hemi_njet, &b_minitree_Hemi_njet);
   fChain->SetBranchAddress("minitree_Hemi_njet_nomu", &minitree_Hemi_njet_nomu, &b_minitree_Hemi_njet_nomu);
   fChain->SetBranchAddress("minitree_Hemi_pt", &minitree_Hemi_pt, &b_minitree_Hemi_pt);
   fChain->SetBranchAddress("minitree_Hemi_eta", &minitree_Hemi_eta, &b_minitree_Hemi_eta);
   fChain->SetBranchAddress("minitree_Hemi_phi", &minitree_Hemi_phi, &b_minitree_Hemi_phi);
   fChain->SetBranchAddress("minitree_Hemi_nTrks", &minitree_Hemi_nTrks, &b_minitree_Hemi_nTrks);
   fChain->SetBranchAddress("minitree_Hemi_nTrks_sig", &minitree_Hemi_nTrks_sig, &b_minitree_Hemi_nTrks_sig);
   fChain->SetBranchAddress("minitree_Hemi_nTrks_bad", &minitree_Hemi_nTrks_bad, &b_minitree_Hemi_nTrks_bad);
   fChain->SetBranchAddress("minitree_Hemi_mass", &minitree_Hemi_mass, &b_minitree_Hemi_mass);
   fChain->SetBranchAddress("minitree_HemiMu_mass", &minitree_HemiMu_mass, &b_minitree_HemiMu_mass);
   fChain->SetBranchAddress("minitree_HemiMu_pt", &minitree_HemiMu_pt, &b_minitree_HemiMu_pt);
   fChain->SetBranchAddress("minitree_HemiMu_dR", &minitree_HemiMu_dR, &b_minitree_HemiMu_dR);
   fChain->SetBranchAddress("minitree_HemiMuOp_mass", &minitree_HemiMuOp_mass, &b_minitree_HemiMuOp_mass);
   fChain->SetBranchAddress("minitree_HemiMuOp_pt", &minitree_HemiMuOp_pt, &b_minitree_HemiMuOp_pt);
   fChain->SetBranchAddress("minitree_HemiMuOp_dR", &minitree_HemiMuOp_dR, &b_minitree_HemiMuOp_dR);
   fChain->SetBranchAddress("minitree_Hemi_dR12", &minitree_Hemi_dR12, &b_minitree_Hemi_dR12);
   fChain->SetBranchAddress("minitree_ll_pt", &minitree_ll_pt, &b_minitree_ll_pt);
   fChain->SetBranchAddress("minitree_ll_eta", &minitree_ll_eta, &b_minitree_ll_eta);
   fChain->SetBranchAddress("minitree_ll_phi", &minitree_ll_phi, &b_minitree_ll_phi);
   fChain->SetBranchAddress("minitree_ll_px", &minitree_ll_px, &b_minitree_ll_px);
   fChain->SetBranchAddress("minitree_ll_py", &minitree_ll_py, &b_minitree_ll_py);
   fChain->SetBranchAddress("minitree_ll_pz", &minitree_ll_pz, &b_minitree_ll_pz);
   fChain->SetBranchAddress("minitree_ll_energy", &minitree_ll_energy, &b_minitree_ll_energy);
   fChain->SetBranchAddress("minitree_ll_mass", &minitree_ll_mass, &b_minitree_ll_mass);
   fChain->SetBranchAddress("minitree_track_ipc", &minitree_track_ipc, &b_minitree_track_ipc);
   fChain->SetBranchAddress("minitree_track_lost", &minitree_track_lost, &b_minitree_track_lost);
   fChain->SetBranchAddress("minitree_track_px", &minitree_track_px, &b_minitree_track_px);
   fChain->SetBranchAddress("minitree_track_py", &minitree_track_py, &b_minitree_track_py);
   fChain->SetBranchAddress("minitree_track_pz", &minitree_track_pz, &b_minitree_track_pz);
   fChain->SetBranchAddress("minitree_track_pt", &minitree_track_pt, &b_minitree_track_pt);
   fChain->SetBranchAddress("minitree_track_eta", &minitree_track_eta, &b_minitree_track_eta);
   fChain->SetBranchAddress("minitree_track_phi", &minitree_track_phi, &b_minitree_track_phi);
   fChain->SetBranchAddress("minitree_track_charge", &minitree_track_charge, &b_minitree_track_charge);
   fChain->SetBranchAddress("minitree_track_NChi2", &minitree_track_NChi2, &b_minitree_track_NChi2);
   fChain->SetBranchAddress("minitree_track_isHighPurity", &minitree_track_isHighPurity, &b_minitree_track_isHighPurity);
   fChain->SetBranchAddress("minitree_track_dxy", &minitree_track_dxy, &b_minitree_track_dxy);
   fChain->SetBranchAddress("minitree_track_dxyError", &minitree_track_dxyError, &b_minitree_track_dxyError);
   fChain->SetBranchAddress("minitree_track_drSig", &minitree_track_drSig, &b_minitree_track_drSig);
   fChain->SetBranchAddress("minitree_track_dz", &minitree_track_dz, &b_minitree_track_dz);
   fChain->SetBranchAddress("minitree_track_dzError", &minitree_track_dzError, &b_minitree_track_dzError);
   fChain->SetBranchAddress("minitree_track_dzSig", &minitree_track_dzSig, &b_minitree_track_dzSig);
   fChain->SetBranchAddress("minitree_track_nHit", &minitree_track_nHit, &b_minitree_track_nHit);
   fChain->SetBranchAddress("minitree_track_nHitPixel", &minitree_track_nHitPixel, &b_minitree_track_nHitPixel);
   fChain->SetBranchAddress("minitree_track_nHitTIB", &minitree_track_nHitTIB, &b_minitree_track_nHitTIB);
   fChain->SetBranchAddress("minitree_track_nHitTID", &minitree_track_nHitTID, &b_minitree_track_nHitTID);
   fChain->SetBranchAddress("minitree_track_nHitTOB", &minitree_track_nHitTOB, &b_minitree_track_nHitTOB);
   fChain->SetBranchAddress("minitree_track_nHitTEC", &minitree_track_nHitTEC, &b_minitree_track_nHitTEC);
   fChain->SetBranchAddress("minitree_track_nHitPXB", &minitree_track_nHitPXB, &b_minitree_track_nHitPXB);
   fChain->SetBranchAddress("minitree_track_nHitPXF", &minitree_track_nHitPXF, &b_minitree_track_nHitPXF);
   fChain->SetBranchAddress("minitree_track_isHitPixel", &minitree_track_isHitPixel, &b_minitree_track_isHitPixel);
   fChain->SetBranchAddress("minitree_track_nLayers", &minitree_track_nLayers, &b_minitree_track_nLayers);
   fChain->SetBranchAddress("minitree_track_nLayersPixel", &minitree_track_nLayersPixel, &b_minitree_track_nLayersPixel);
   fChain->SetBranchAddress("minitree_track_x", &minitree_track_x, &b_minitree_track_x);
   fChain->SetBranchAddress("minitree_track_y", &minitree_track_y, &b_minitree_track_y);
   fChain->SetBranchAddress("minitree_track_z", &minitree_track_z, &b_minitree_track_z);
   fChain->SetBranchAddress("minitree_track_firstHit", &minitree_track_firstHit, &b_minitree_track_firstHit);
   fChain->SetBranchAddress("minitree_track_region", &minitree_track_region, &b_minitree_track_region);
   fChain->SetBranchAddress("minitree_track_firstHit_x", &minitree_track_firstHit_x, &b_minitree_track_firstHit_x);
   fChain->SetBranchAddress("minitree_track_firstHit_y", &minitree_track_firstHit_y, &b_minitree_track_firstHit_y);
   fChain->SetBranchAddress("minitree_track_firstHit_z", &minitree_track_firstHit_z, &b_minitree_track_firstHit_z);
   fChain->SetBranchAddress("minitree_track_iJet", &minitree_track_iJet, &b_minitree_track_iJet);
   fChain->SetBranchAddress("minitree_track_ntrk10", &minitree_track_ntrk10, &b_minitree_track_ntrk10);
   fChain->SetBranchAddress("minitree_track_ntrk20", &minitree_track_ntrk20, &b_minitree_track_ntrk20);
   fChain->SetBranchAddress("minitree_track_ntrk30", &minitree_track_ntrk30, &b_minitree_track_ntrk30);
   fChain->SetBranchAddress("minitree_track_ntrk40", &minitree_track_ntrk40, &b_minitree_track_ntrk40);
   fChain->SetBranchAddress("minitree_track_MVAval", &minitree_track_MVAval, &b_minitree_track_MVAval);
   fChain->SetBranchAddress("minitree_track_Hemi_dR", &minitree_track_Hemi_dR, &b_minitree_track_Hemi_dR);
   fChain->SetBranchAddress("minitree_track_Hemi_dRmax", &minitree_track_Hemi_dRmax, &b_minitree_track_Hemi_dRmax);
   fChain->SetBranchAddress("minitree_K0_mass", &minitree_K0_mass, &b_minitree_K0_mass);
   fChain->SetBranchAddress("minitree_K0_pt", &minitree_K0_pt, &b_minitree_K0_pt);
   fChain->SetBranchAddress("minitree_L0_mass", &minitree_L0_mass, &b_minitree_L0_mass);
   fChain->SetBranchAddress("minitree_L0_pt", &minitree_L0_pt, &b_minitree_L0_pt);
   fChain->SetBranchAddress("minitree_V0_reco_mass", &minitree_V0_reco_mass, &b_minitree_V0_reco_mass);
   fChain->SetBranchAddress("minitree_V0_reco_pt", &minitree_V0_reco_pt, &b_minitree_V0_reco_pt);
   fChain->SetBranchAddress("minitree_V0_reco_source", &minitree_V0_reco_source, &b_minitree_V0_reco_source);
   fChain->SetBranchAddress("minitree_SecInt_mass", &minitree_SecInt_mass, &b_minitree_SecInt_mass);
   fChain->SetBranchAddress("minitree_SecInt_pt", &minitree_SecInt_pt, &b_minitree_SecInt_pt);
   fChain->SetBranchAddress("minitree_SecInt_r", &minitree_SecInt_r, &b_minitree_SecInt_r);
   fChain->SetBranchAddress("minitree_SecInt_z", &minitree_SecInt_z, &b_minitree_SecInt_z);
   fChain->SetBranchAddress("minitree_SecInt_drSig", &minitree_SecInt_drSig, &b_minitree_SecInt_drSig);
   fChain->SetBranchAddress("minitree_SecInt_dzSig", &minitree_SecInt_dzSig, &b_minitree_SecInt_dzSig);
   fChain->SetBranchAddress("minitree_SecInt_layer", &minitree_SecInt_layer, &b_minitree_SecInt_layer);
   fChain->SetBranchAddress("minitree_SecInt_selec", &minitree_SecInt_selec, &b_minitree_SecInt_selec);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_step", &minitree_Hemi_Vtx_step, &b_minitree_Hemi_Vtx_step);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_NChi2", &minitree_Hemi_Vtx_NChi2, &b_minitree_Hemi_Vtx_NChi2);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_nTrks", &minitree_Hemi_Vtx_nTrks, &b_minitree_Hemi_Vtx_nTrks);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_dR", &minitree_Hemi_Vtx_dR, &b_minitree_Hemi_Vtx_dR);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_xError", &minitree_Hemi_Vtx_xError, &b_minitree_Hemi_Vtx_xError);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_yError", &minitree_Hemi_Vtx_yError, &b_minitree_Hemi_Vtx_yError);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_zError", &minitree_Hemi_Vtx_zError, &b_minitree_Hemi_Vtx_zError);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_trackWeight", &minitree_Hemi_Vtx_trackWeight, &b_minitree_Hemi_Vtx_trackWeight);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_SumtrackWeight", &minitree_Hemi_Vtx_SumtrackWeight, &b_minitree_Hemi_Vtx_SumtrackWeight);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_track_MeanDCA_d", &minitree_Hemi_Vtx_track_MeanDCA_d, &b_minitree_Hemi_Vtx_track_MeanDCA_d);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_Mass", &minitree_Hemi_Vtx_Mass, &b_minitree_Hemi_Vtx_Mass);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_dist", &minitree_Hemi_Vtx_dist, &b_minitree_Hemi_Vtx_dist);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_ntrk10", &minitree_Hemi_Vtx_ntrk10, &b_minitree_Hemi_Vtx_ntrk10);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_ntrk20", &minitree_Hemi_Vtx_ntrk20, &b_minitree_Hemi_Vtx_ntrk20);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_step", &minitree_Hemi_SecVtx_step, &b_minitree_Hemi_SecVtx_step);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_NChi2", &minitree_Hemi_SecVtx_NChi2, &b_minitree_Hemi_SecVtx_NChi2);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_dR", &minitree_Hemi_SecVtx_dR, &b_minitree_Hemi_SecVtx_dR);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_nTrks", &minitree_Hemi_SecVtx_nTrks, &b_minitree_Hemi_SecVtx_nTrks);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_dist", &minitree_Hemi_SecVtx_dist, &b_minitree_Hemi_SecVtx_dist);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_track_MeanDCA_d", &minitree_Hemi_SecVtx_track_MeanDCA_d, &b_minitree_Hemi_SecVtx_track_MeanDCA_d);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_SumtrackWeight", &minitree_Hemi_SecVtx_SumtrackWeight, &b_minitree_Hemi_SecVtx_SumtrackWeight);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_trackWeight", &minitree_Hemi_SecVtx_trackWeight, &b_minitree_Hemi_SecVtx_trackWeight);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_Mass", &minitree_Hemi_SecVtx_Mass, &b_minitree_Hemi_SecVtx_Mass);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_nTrks", &minitree_Hemi_Vtx_BDT_nTrks, &b_minitree_Hemi_Vtx_BDT_nTrks);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_NChi2", &minitree_Hemi_Vtx_BDT_NChi2, &b_minitree_Hemi_Vtx_BDT_NChi2);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_step", &minitree_Hemi_Vtx_BDT_step, &b_minitree_Hemi_Vtx_BDT_step);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_Mass", &minitree_Hemi_Vtx_BDT_Mass, &b_minitree_Hemi_Vtx_BDT_Mass);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_HMass", &minitree_Hemi_Vtx_BDT_HMass, &b_minitree_Hemi_Vtx_BDT_HMass);
   Notify();
}

Bool_t DATAMCReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DATAMCReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DATAMCReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DATAMCReader_cxx
