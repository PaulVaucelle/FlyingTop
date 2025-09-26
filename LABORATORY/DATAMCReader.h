//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar  5 14:12:15 2025 by ROOT version 6.14/09
// from TTree ttree/ttree
// found on file: miniDYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root
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

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventNumber;
   Int_t           lumiBlock;
   vector<float>   *tree_LHE_Weights;
   Float_t         tree_MCEvt_weight;
   Double_t        tree_only_gen_wt;
   Double_t        tree_event_weight;
   Double_t        tree_genTop_Weight;
   vector<float>   *tree_gen_top_pt;
   vector<float>   *tree_gen_top_rw_pt;
   Double_t        PUweight;
   Double_t        PUweight_Up;
   Double_t        PUweight_Down;
   Double_t        Prefweight;
   Double_t        Prefweight_Up;
   Double_t        Prefweight_Down;
   Int_t           PU_events;
   Bool_t          tree_Filter;
   Bool_t          tree_FilterSameSign;
   Bool_t          tree_trigger_doublelepton;
   Bool_t          tree_trigger_singlelepton;
   Bool_t          tree_Good_PV;
   Int_t           tree_nPV;
   Float_t         tree_PV_x;
   Float_t         tree_PV_y;
   Float_t         tree_PV_z;
   Float_t         tree_PV_ez;
   Float_t         tree_PV_NChi2;
   Int_t           tree_PV_ndf;
   Float_t         tree_PFMet_et;
   Float_t         tree_PFMet_phi;
   Float_t         tree_HT;
   Int_t           tree_TRACK_SIZE;
   Int_t           tree_nTracks;
   Int_t           tree_nLostTracks;
   Int_t           tree_muon_GenRecoTriggerMatched;
   Int_t           tree_all_nmu;
   Int_t           tree_nmu;
   Float_t         tree_LT;
   Float_t         tree_Mmumu;
   Float_t         tree_MmumuSameSign;
   vector<bool>    *tree_muon_isPrompt;
   vector<float>   *tree_muon_pt;
   vector<float>   *tree_muon_SF;
   vector<float>   *tree_muon_eta;
   vector<float>   *tree_muon_phi;
   vector<float>   *tree_muon_dxy;
   vector<float>   *tree_muon_dz;
   vector<int>     *tree_muon_charge;
   vector<float>   *tree_muon_correction;
   vector<int>     *tree_muon_gen;
   vector<float>   *tree_muon_dxyError;
   vector<float>   *tree_muon_dzError;
   vector<bool>    *tree_muon_isLoose;
   vector<bool>    *tree_muon_isMedium;
   vector<bool>    *tree_muon_isTight;
   vector<bool>    *tree_muon_isGlobal;
   vector<bool>    *tree_muon_PFIsoVeryLoose;
   vector<bool>    *tree_muon_PFIsoLoose;
   vector<bool>    *tree_muon_PFIsoMedium;
   vector<bool>    *tree_muon_PFIsoTight;
   vector<bool>    *tree_muon_TkIsoLoose;
   vector<bool>    *tree_muon_TkIsoTight;
   vector<float>   *tree_lepton_leadingpt;
   vector<float>   *tree_lepton_leadingpt2;
   vector<float>   *tree_lepton_leadingeta;
   vector<float>   *tree_lepton_leadingeta2;
   vector<float>   *tree_lepton_leadingphi;
   vector<float>   *tree_lepton_leadingphi2;
   vector<float>   *tree_lepton_lepton_dR;
   vector<float>   *tree_lepton_lepton_dPhi;
   vector<float>   *tree_lepton_lepton_dEta;
   vector<float>   *tree_lepton_leadingdxy;
   vector<float>   *tree_lepton_leadingdxy2;
   vector<float>   *tree_lepton_leadingdz;
   vector<float>   *tree_lepton_leadingdz2;
   Int_t           tree_all_nel;
   Int_t           tree_electron_nEle;
   vector<bool>    *tree_electron_isPrompt;
   vector<float>   *tree_electron_pt;
   vector<float>   *tree_electron_eta;
   vector<float>   *tree_electron_phi;
   vector<int>     *tree_electron_charge;
   vector<float>   *tree_electron_dxy;
   vector<float>   *tree_electron_dz;
   vector<int>     *tree_electron_gen;
   vector<float>   *tree_electron_energy;
   vector<float>   *tree_electron_et;
   vector<float>   *tree_electron_ecal_trk_postcorr;
   vector<float>   *tree_electron_isoR4;
   vector<bool>    *tree_electron_IsLoose;
   vector<bool>    *tree_electron_IsMedium;
   vector<bool>    *tree_electron_IsTight;
   Int_t           tree_njet;
   Int_t           tree_njetNOmu;
   vector<float>   *tree_jet_pt;
   vector<float>   *tree_jet_eta;
   vector<float>   *tree_jet_phi;
   vector<float>   *tree_jet_HadronFlavour;
   vector<float>   *tree_jet_btag_DeepJet;
   vector<float>   *tree_jet_E;
   vector<float>   *tree_jet_leadingpt;
   vector<float>   *tree_jet_leadingpt2;
   vector<float>   *tree_jet_leadingeta;
   vector<float>   *tree_jet_leadingeta2;
   vector<float>   *tree_jet_jet_dR;
   vector<float>   *tree_jet_jet_dPhi;
   vector<float>   *tree_jet_jet_dEta;
   vector<float>   *tree_muon_jet_dRmin;
   vector<float>   *tree_muon_jet_dRmax;
   Float_t         tree_Evts_MVAval;
   Float_t         tree_Evts_MVAvalDY;
   Float_t         tree_Evts_MVAvalTT;
   vector<int>     *tree_Hemi;
   vector<int>     *tree_Hemi_njet;
   vector<int>     *tree_Hemi_njet_nomu;
   vector<float>   *tree_Hemi_pt;
   vector<float>   *tree_Hemi_eta;
   vector<float>   *tree_Hemi_phi;
   vector<int>     *tree_Hemi_nTrks;
   vector<int>     *tree_Hemi_nTrks_sig;
   vector<int>     *tree_Hemi_nTrks_bad;
   vector<float>   *tree_Hemi_mass;
   vector<float>   *tree_HemiMu_mass;
   vector<float>   *tree_HemiMu_pt;
   vector<float>   *tree_HemiMu_dR;
   vector<float>   *tree_HemiMuOp_mass;
   vector<float>   *tree_HemiMuOp_pt;
   vector<float>   *tree_HemiMuOp_dR;
   vector<float>   *tree_Hemi_dR12;
   vector<int>     *tree_Hemi_Vtx_step;
   vector<bool>    *tree_Hemi_Vtx_isTight;
   vector<float>   *tree_Hemi_Vtx_NChi2;
   vector<int>     *tree_Hemi_Vtx_nTrks;
   vector<int>     *tree_Hemi_Vtx_nTrks_sig;
   vector<int>     *tree_Hemi_Vtx_nTrks_bad;
   vector<float>   *tree_Hemi_Vtx_x;
   vector<float>   *tree_Hemi_Vtx_y;
   vector<float>   *tree_Hemi_Vtx_z;
   vector<float>   *tree_Hemi_Vtx_r;
   vector<float>   *tree_Hemi_Vtx_dR;
   vector<float>   *tree_Hemi_Vtx_SumtrackWeight;
   vector<float>   *tree_Hemi_Vtx_track_MeanDCA_d;
   vector<float>   *tree_Hemi_Vtx_Mass;
   vector<float>   *tree_Hemi_Vtx_dist;
   vector<int>     *tree_event_nVtx;
   vector<float>   *tree_event_Vtx_Vtx_dr;
   vector<float>   *tree_event_Vtx_Vtx_dz;
   vector<float>   *tree_event_Vtx_Vtx_dd;
   vector<float>   *tree_event_Vtx_Vtx_reldd;
   vector<float>   *tree_event_Vtx_Vtx_dR;
   vector<int>     *tree_event_Vtx_Vtx_step;
   vector<float>   *tree_Hemi_SecLLP;
   vector<float>   *tree_Hemi_LLP_SecVtx_dz;
   vector<float>   *tree_Hemi_LLP_SecVtx_dr;
   vector<bool>    *tree_Hemi_SecLLP_ping;
   vector<int>     *tree_event_SecLLP_ping;
   vector<int>     *tree_Hemi_SecVtx;
   vector<int>     *tree_Hemi_SecVtx_step;
   vector<float>   *tree_Hemi_SecVtx_x;
   vector<float>   *tree_Hemi_SecVtx_y;
   vector<float>   *tree_Hemi_SecVtx_z;
   vector<float>   *tree_Hemi_SecVtx_r;
   vector<float>   *tree_Hemi_SecVtx_dR;
   vector<float>   *tree_Hemi_SecVtx_nTrks;
   vector<float>   *tree_Hemi_SecVtx_NChi2;
   vector<float>   *tree_Hemi_SecVtx_dist;
   vector<float>   *tree_Hemi_SecVtx_track_MeanDCA_d;
   vector<float>   *tree_Hemi_SecVtx_SumtrackWeight;
   vector<float>   *tree_Hemi_SecVtx_Mass;
   vector<float>   *tree_event_MergedVtx_Vtx_dr;
   vector<float>   *tree_event_MergedVtx_Vtx_dz;
   vector<float>   *tree_event_MergedVtx_Vtx_dd;
   vector<float>   *tree_event_MergedVtx_Vtx_reldd;
   vector<float>   *tree_event_MergedVtx_Vtx_dR;
   vector<int>     *tree_event_MergedVtx_Vtx_step;
   vector<float>   *tree_Hemi_Vtx_BDT_nTrks;
   vector<float>   *tree_Hemi_Vtx_BDT_NChi2;
   vector<float>   *tree_Hemi_Vtx_BDT_step;
   vector<float>   *tree_Hemi_Vtx_BDT_STW;
   vector<float>   *tree_Hemi_Vtx_BDT_Mass;
   vector<float>   *tree_Hemi_Vtx_BDT_HMass;
   vector<float>   *tree_Hemi_Vtx_BDT_ntrk10;
   vector<float>   *tree_Hemi_Vtx_BDT_ntrk20;
   vector<float>   *tree_Hemi_Vtx_BDT_MeanDCA;
   vector<float>   *tree_Hemi_Vtx_MVAval_Loose;
   vector<float>   *tree_Hemi_Vtx_MVAval_Tight;
   vector<unsigned int> *tree_track_ipc;
   vector<bool>    *tree_track_lost;
   vector<float>   *tree_track_px;
   vector<float>   *tree_track_py;
   vector<float>   *tree_track_pz;
   vector<float>   *tree_track_pt;
   vector<float>   *tree_track_eta;
   vector<float>   *tree_track_phi;
   vector<int>     *tree_track_charge;
   vector<float>   *tree_track_NChi2;
   vector<bool>    *tree_track_isHighPurity;
   vector<float>   *tree_track_dxy;
   vector<float>   *tree_track_dxyError;
   vector<float>   *tree_track_drSig;
   vector<float>   *tree_track_dz;
   vector<float>   *tree_track_dzError;
   vector<float>   *tree_track_dzSig;
   vector<int>     *tree_track_nHit;
   vector<int>     *tree_track_nHitPixel;
   vector<int>     *tree_track_nHitTIB;
   vector<int>     *tree_track_nHitTID;
   vector<int>     *tree_track_nHitTOB;
   vector<int>     *tree_track_nHitTEC;
   vector<int>     *tree_track_nHitPXB;
   vector<int>     *tree_track_nHitPXF;
   vector<int>     *tree_track_isHitPixel;
   vector<int>     *tree_track_nLayers;
   vector<int>     *tree_track_nLayersPixel;
   vector<float>   *tree_track_x;
   vector<float>   *tree_track_y;
   vector<float>   *tree_track_z;
   vector<int>     *tree_track_firstHit;
   vector<float>   *tree_track_region;
   vector<float>   *tree_track_firstHit_x;
   vector<float>   *tree_track_firstHit_y;
   vector<float>   *tree_track_firstHit_z;
   vector<int>     *tree_track_iJet;
   vector<float>   *tree_track_ntrk10;
   vector<float>   *tree_track_ntrk20;
   vector<float>   *tree_track_ntrk30;
   vector<float>   *tree_track_ntrk40;
   vector<double>  *tree_track_MVAval;
   vector<float>   *tree_track_Hemi_dR;
   vector<float>   *tree_track_Hemi_dRmax;
   vector<float>   *tree_K0_mass;
   vector<float>   *tree_K0_pt;
   vector<float>   *tree_L0_mass;
   vector<float>   *tree_L0_pt;
   vector<float>   *tree_V0_reco_mass;
   vector<float>   *tree_V0_reco_pt;
   vector<int>     *tree_V0_reco_source;
   vector<float>   *tree_SecInt_mass;
   vector<float>   *tree_SecInt_pt;
   vector<float>   *tree_SecInt_drSig;
   vector<float>   *tree_SecInt_dzSig;
   vector<int>     *tree_SecInt_layer;
   vector<bool>    *tree_SecInt_selec;
   vector<float>   *tree_SecInt_r;
   vector<float>   *tree_SecInt_z;

   // List of branches
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_tree_LHE_Weights;   //!
   TBranch        *b_tree_MCEvt_weight;   //!
   TBranch        *b_tree_only_gen_wt;   //!
   TBranch        *b_tree_event_weight;   //!
   TBranch        *b_tree_genTop_Weight;   //!
   TBranch        *b_tree_gen_top_pt;   //!
   TBranch        *b_tree_gen_top_rw_pt;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_PUweight_Up;   //!
   TBranch        *b_PUweight_Down;   //!
   TBranch        *b_Prefweight;   //!
   TBranch        *b_Prefweight_Up;   //!
   TBranch        *b_Prefweight_Down;   //!
   TBranch        *b_PU_events;   //!
   TBranch        *b_tree_Filter;   //!
   TBranch        *b_tree_FilterSameSign;   //!
   TBranch        *b_tree_trigger_doublelepton;   //!
   TBranch        *b_tree_trigger_singlelepton;   //!
   TBranch        *b_tree_Good_PV;   //!
   TBranch        *b_tree_nPV;   //!
   TBranch        *b_tree_PV_x;   //!
   TBranch        *b_tree_PV_y;   //!
   TBranch        *b_tree_PV_z;   //!
   TBranch        *b_tree_PV_ez;   //!
   TBranch        *b_tree_PV_NChi2;   //!
   TBranch        *b_tree_PV_ndf;   //!
   TBranch        *b_tree_PFMet_et;   //!
   TBranch        *b_tree_PFMet_phi;   //!
   TBranch        *b_tree_HT;   //!
   TBranch        *b_tree_TRACK_SIZE;   //!
   TBranch        *b_tree_nTracks;   //!
   TBranch        *b_tree_nLostTracks;   //!
   TBranch        *b_tree_muon_GenRecoTriggerMatched;   //!
   TBranch        *b_tree_all_nmu;   //!
   TBranch        *b_tree_nmu;   //!
   TBranch        *b_tree_LT;   //!
   TBranch        *b_tree_Mmumu;   //!
   TBranch        *b_tree_MmumuSameSign;   //!
   TBranch        *b_tree_muon_isPrompt;   //!
   TBranch        *b_tree_muon_pt;   //!
   TBranch        *b_tree_muon_SF;   //!
   TBranch        *b_tree_muon_eta;   //!
   TBranch        *b_tree_muon_phi;   //!
   TBranch        *b_tree_muon_dxy;   //!
   TBranch        *b_tree_muon_dz;   //!
   TBranch        *b_tree_muon_charge;   //!
   TBranch        *b_tree_muon_correction;   //!
   TBranch        *b_tree_muon_gen;   //!
   TBranch        *b_tree_muon_dxyError;   //!
   TBranch        *b_tree_muon_dzError;   //!
   TBranch        *b_tree_muon_isLoose;   //!
   TBranch        *b_tree_muon_isMedium;   //!
   TBranch        *b_tree_muon_isTight;   //!
   TBranch        *b_tree_muon_isGlobal;   //!
   TBranch        *b_tree_muon_PFIsoVeryLoose;   //!
   TBranch        *b_tree_muon_PFIsoLoose;   //!
   TBranch        *b_tree_muon_PFIsoMedium;   //!
   TBranch        *b_tree_muon_PFIsoTight;   //!
   TBranch        *b_tree_muon_TkIsoLoose;   //!
   TBranch        *b_tree_muon_TkIsoTight;   //!
   TBranch        *b_tree_lepton_leadingpt;   //!
   TBranch        *b_tree_lepton_leadingpt2;   //!
   TBranch        *b_tree_lepton_leadingeta;   //!
   TBranch        *b_tree_lepton_leadingeta2;   //!
   TBranch        *b_tree_lepton_leadingphi;   //!
   TBranch        *b_tree_lepton_leadingphi2;   //!
   TBranch        *b_tree_lepton_lepton_dR;   //!
   TBranch        *b_tree_lepton_lepton_dPhi;   //!
   TBranch        *b_tree_lepton_lepton_dEta;   //!
   TBranch        *b_tree_lepton_leadingdxy;   //!
   TBranch        *b_tree_lepton_leadingdxy2;   //!
   TBranch        *b_tree_lepton_leadingdz;   //!
   TBranch        *b_tree_lepton_leadingdz2;   //!
   TBranch        *b_tree_all_nel;   //!
   TBranch        *b_tree_electron_nEle;   //!
   TBranch        *b_tree_electron_isPrompt;   //!
   TBranch        *b_tree_electron_pt;   //!
   TBranch        *b_tree_electron_eta;   //!
   TBranch        *b_tree_electron_phi;   //!
   TBranch        *b_tree_electron_charge;   //!
   TBranch        *b_tree_electron_dxy;   //!
   TBranch        *b_tree_electron_dz;   //!
   TBranch        *b_tree_electron_gen;   //!
   TBranch        *b_tree_electron_energy;   //!
   TBranch        *b_tree_electron_et;   //!
   TBranch        *b_tree_electron_ecal_trk_postcorr;   //!
   TBranch        *b_tree_electron_isoR4;   //!
   TBranch        *b_tree_electron_IsLoose;   //!
   TBranch        *b_tree_electron_IsMedium;   //!
   TBranch        *b_tree_electron_IsTight;   //!
   TBranch        *b_tree_njet;   //!
   TBranch        *b_tree_njetNOmu;   //!
   TBranch        *b_tree_jet_pt;   //!
   TBranch        *b_tree_jet_eta;   //!
   TBranch        *b_tree_jet_phi;   //!
   TBranch        *b_tree_jet_HadronFlavour;   //!
   TBranch        *b_tree_jet_btag_DeepJet;   //!
   TBranch        *b_tree_jet_E;   //!
   TBranch        *b_tree_jet_leadingpt;   //!
   TBranch        *b_tree_jet_leadingpt2;   //!
   TBranch        *b_tree_jet_leadingeta;   //!
   TBranch        *b_tree_jet_leadingeta2;   //!
   TBranch        *b_tree_jet_jet_dR;   //!
   TBranch        *b_tree_jet_jet_dPhi;   //!
   TBranch        *b_tree_jet_jet_dEta;   //!
   TBranch        *b_tree_muon_jet_dRmin;   //!
   TBranch        *b_tree_muon_jet_dRmax;   //!
   TBranch        *b_tree_Evts_MVAval;   //!
   TBranch        *b_tree_Evts_MVAvalDY;   //!
   TBranch        *b_tree_Evts_MVAvalTT;   //!
   TBranch        *b_tree_Hemi;   //!
   TBranch        *b_tree_Hemi_njet;   //!
   TBranch        *b_tree_Hemi_njet_nomu;   //!
   TBranch        *b_tree_Hemi_pt;   //!
   TBranch        *b_tree_Hemi_eta;   //!
   TBranch        *b_tree_Hemi_phi;   //!
   TBranch        *b_tree_Hemi_nTrks;   //!
   TBranch        *b_tree_Hemi_nTrks_sig;   //!
   TBranch        *b_tree_Hemi_nTrks_bad;   //!
   TBranch        *b_tree_Hemi_mass;   //!
   TBranch        *b_tree_HemiMu_mass;   //!
   TBranch        *b_tree_HemiMu_pt;   //!
   TBranch        *b_tree_HemiMu_dR;   //!
   TBranch        *b_tree_HemiMuOp_mass;   //!
   TBranch        *b_tree_HemiMuOp_pt;   //!
   TBranch        *b_tree_HemiMuOp_dR;   //!
   TBranch        *b_tree_Hemi_dR12;   //!
   TBranch        *b_tree_Hemi_Vtx_step;   //!
   TBranch        *b_tree_Hemi_Vtx_isTight;   //!
   TBranch        *b_tree_Hemi_Vtx_NChi2;   //!
   TBranch        *b_tree_Hemi_Vtx_nTrks;   //!
   TBranch        *b_tree_Hemi_Vtx_nTrks_sig;   //!
   TBranch        *b_tree_Hemi_Vtx_nTrks_bad;   //!
   TBranch        *b_tree_Hemi_Vtx_x;   //!
   TBranch        *b_tree_Hemi_Vtx_y;   //!
   TBranch        *b_tree_Hemi_Vtx_z;   //!
   TBranch        *b_tree_Hemi_Vtx_r;   //!
   TBranch        *b_tree_Hemi_Vtx_dR;   //!
   TBranch        *b_tree_Hemi_Vtx_SumtrackWeight;   //!
   TBranch        *b_tree_Hemi_Vtx_track_MeanDCA_d;   //!
   TBranch        *b_tree_Hemi_Vtx_Mass;   //!
   TBranch        *b_tree_Hemi_Vtx_dist;   //!
   TBranch        *b_tree_event_nVtx;   //!
   TBranch        *b_tree_event_Vtx_Vtx_dr;   //!
   TBranch        *b_tree_event_Vtx_Vtx_dz;   //!
   TBranch        *b_tree_event_Vtx_Vtx_dd;   //!
   TBranch        *b_tree_event_Vtx_Vtx_reldd;   //!
   TBranch        *b_tree_event_Vtx_Vtx_dR;   //!
   TBranch        *b_tree_event_Vtx_Vtx_step;   //!
   TBranch        *b_tree_Hemi_SecLLP;   //!
   TBranch        *b_tree_Hemi_LLP_SecVtx_dz;   //!
   TBranch        *b_tree_Hemi_LLP_SecVtx_dr;   //!
   TBranch        *b_tree_Hemi_SecLLP_ping;   //!
   TBranch        *b_tree_event_SecLLP_ping;   //!
   TBranch        *b_tree_Hemi_SecVtx;   //!
   TBranch        *b_tree_Hemi_SecVtx_step;   //!
   TBranch        *b_tree_Hemi_SecVtx_x;   //!
   TBranch        *b_tree_Hemi_SecVtx_y;   //!
   TBranch        *b_tree_Hemi_SecVtx_z;   //!
   TBranch        *b_tree_Hemi_SecVtx_r;   //!
   TBranch        *b_tree_Hemi_SecVtx_dR;   //!
   TBranch        *b_tree_Hemi_SecVtx_nTrks;   //!
   TBranch        *b_tree_Hemi_SecVtx_NChi2;   //!
   TBranch        *b_tree_Hemi_SecVtx_dist;   //!
   TBranch        *b_tree_Hemi_SecVtx_track_MeanDCA_d;   //!
   TBranch        *b_tree_Hemi_SecVtx_SumtrackWeight;   //!
   TBranch        *b_tree_Hemi_SecVtx_Mass;   //!
   TBranch        *b_tree_event_MergedVtx_Vtx_dr;   //!
   TBranch        *b_tree_event_MergedVtx_Vtx_dz;   //!
   TBranch        *b_tree_event_MergedVtx_Vtx_dd;   //!
   TBranch        *b_tree_event_MergedVtx_Vtx_reldd;   //!
   TBranch        *b_tree_event_MergedVtx_Vtx_dR;   //!
   TBranch        *b_tree_event_MergedVtx_Vtx_step;   //!
   TBranch        *b_tree_Hemi_Vtx_BDT_nTrks;   //!
   TBranch        *b_tree_Hemi_Vtx_BDT_NChi2;   //!
   TBranch        *b_tree_Hemi_Vtx_BDT_step;   //!
   TBranch        *b_tree_Hemi_Vtx_BDT_STW;   //!
   TBranch        *b_tree_Hemi_Vtx_BDT_Mass;   //!
   TBranch        *b_tree_Hemi_Vtx_BDT_HMass;   //!
   TBranch        *b_tree_Hemi_Vtx_BDT_ntrk10;   //!
   TBranch        *b_tree_Hemi_Vtx_BDT_ntrk20;   //!
   TBranch        *b_tree_Hemi_Vtx_BDT_MeanDCA;   //!
   TBranch        *b_tree_Hemi_Vtx_MVAval_Loose;   //!
   TBranch        *b_tree_Hemi_Vtx_MVAval_Tight;   //!
   TBranch        *b_tree_track_ipc;   //!
   TBranch        *b_tree_track_lost;   //!
   TBranch        *b_tree_track_px;   //!
   TBranch        *b_tree_track_py;   //!
   TBranch        *b_tree_track_pz;   //!
   TBranch        *b_tree_track_pt;   //!
   TBranch        *b_tree_track_eta;   //!
   TBranch        *b_tree_track_phi;   //!
   TBranch        *b_tree_track_charge;   //!
   TBranch        *b_tree_track_NChi2;   //!
   TBranch        *b_tree_track_isHighPurity;   //!
   TBranch        *b_tree_track_dxy;   //!
   TBranch        *b_tree_track_dxyError;   //!
   TBranch        *b_tree_track_drSig;   //!
   TBranch        *b_tree_track_dz;   //!
   TBranch        *b_tree_track_dzError;   //!
   TBranch        *b_tree_track_dzSig;   //!
   TBranch        *b_tree_track_nHit;   //!
   TBranch        *b_tree_track_nHitPixel;   //!
   TBranch        *b_tree_track_nHitTIB;   //!
   TBranch        *b_tree_track_nHitTID;   //!
   TBranch        *b_tree_track_nHitTOB;   //!
   TBranch        *b_tree_track_nHitTEC;   //!
   TBranch        *b_tree_track_nHitPXB;   //!
   TBranch        *b_tree_track_nHitPXF;   //!
   TBranch        *b_tree_track_isHitPixel;   //!
   TBranch        *b_tree_track_nLayers;   //!
   TBranch        *b_tree_track_nLayersPixel;   //!
   TBranch        *b_tree_track_x;   //!
   TBranch        *b_tree_track_y;   //!
   TBranch        *b_tree_track_z;   //!
   TBranch        *b_tree_track_firstHit;   //!
   TBranch        *b_tree_track_region;   //!
   TBranch        *b_tree_track_firstHit_x;   //!
   TBranch        *b_tree_track_firstHit_y;   //!
   TBranch        *b_tree_track_firstHit_z;   //!
   TBranch        *b_tree_track_iJet;   //!
   TBranch        *b_tree_track_ntrk10;   //!
   TBranch        *b_tree_track_ntrk20;   //!
   TBranch        *b_tree_track_ntrk30;   //!
   TBranch        *b_tree_track_ntrk40;   //!
   TBranch        *b_tree_track_MVAval;   //!
   TBranch        *b_tree_track_Hemi_dR;   //!
   TBranch        *b_tree_track_Hemi_dRmax;   //!
   TBranch        *b_tree_K0_mass;   //!
   TBranch        *b_tree_K0_pt;   //!
   TBranch        *b_tree_L0_mass;   //!
   TBranch        *b_tree_L0_pt;   //!
   TBranch        *b_tree_V0_reco_mass;   //!
   TBranch        *b_tree_V0_reco_pt;   //!
   TBranch        *b_tree_V0_reco_source;   //!
   TBranch        *b_tree_SecInt_mass;   //!
   TBranch        *b_tree_SecInt_pt;   //!
   TBranch        *b_tree_SecInt_drSig;   //!
   TBranch        *b_tree_SecInt_dzSig;   //!
   TBranch        *b_tree_SecInt_layer;   //!
   TBranch        *b_tree_SecInt_selec;   //!
   TBranch        *b_tree_SecInt_r;   //!
   TBranch        *b_tree_SecInt_z;   //!

   DATAMCReader(TTree *tree=0, TString sample="", TString thesystlist = "");
   virtual ~DATAMCReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(bool isMC, TString Prod, TString sample, bool Signal, int Year, bool IsPostAPV, float MeanGenW , int Channel,bool DoubleMuon,  TString thesystlist);
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


   TString systlist;
};

#endif

#ifdef DATAMCReader_cxx
DATAMCReader::DATAMCReader(TTree *tree, TString sample, TString thesystlist) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("miniDYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("miniDYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
      }
      f->GetObject("ttree",tree);
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
   tree_LHE_Weights = 0;
   tree_gen_top_pt = 0;
   tree_gen_top_rw_pt = 0;
   tree_muon_isPrompt = 0;
   tree_muon_pt = 0;
   tree_muon_SF = 0;
   tree_muon_eta = 0;
   tree_muon_phi = 0;
   tree_muon_dxy = 0;
   tree_muon_dz = 0;
   tree_muon_charge = 0;
   tree_muon_correction = 0;
   tree_muon_gen = 0;
   tree_muon_dxyError = 0;
   tree_muon_dzError = 0;
   tree_muon_isLoose = 0;
   tree_muon_isMedium = 0;
   tree_muon_isTight = 0;
   tree_muon_isGlobal = 0;
   tree_muon_PFIsoVeryLoose = 0;
   tree_muon_PFIsoLoose = 0;
   tree_muon_PFIsoMedium = 0;
   tree_muon_PFIsoTight = 0;
   tree_muon_TkIsoLoose = 0;
   tree_muon_TkIsoTight = 0;
   tree_lepton_leadingpt = 0;
   tree_lepton_leadingpt2 = 0;
   tree_lepton_leadingeta = 0;
   tree_lepton_leadingeta2 = 0;
   tree_lepton_leadingphi = 0;
   tree_lepton_leadingphi2 = 0;
   tree_lepton_lepton_dR = 0;
   tree_lepton_lepton_dPhi = 0;
   tree_lepton_lepton_dEta = 0;
   tree_lepton_leadingdxy = 0;
   tree_lepton_leadingdxy2 = 0;
   tree_lepton_leadingdz = 0;
   tree_lepton_leadingdz2 = 0;
   tree_electron_isPrompt = 0;
   tree_electron_pt = 0;
   tree_electron_eta = 0;
   tree_electron_phi = 0;
   tree_electron_charge = 0;
   tree_electron_dxy = 0;
   tree_electron_dz = 0;
   tree_electron_gen = 0;
   tree_electron_energy = 0;
   tree_electron_et = 0;
   tree_electron_ecal_trk_postcorr = 0;
   tree_electron_isoR4 = 0;
   tree_electron_IsLoose = 0;
   tree_electron_IsMedium = 0;
   tree_electron_IsTight = 0;
   tree_jet_pt = 0;
   tree_jet_eta = 0;
   tree_jet_phi = 0;
   tree_jet_HadronFlavour = 0;
   tree_jet_btag_DeepJet = 0;
   tree_jet_E = 0;
   tree_jet_leadingpt = 0;
   tree_jet_leadingpt2 = 0;
   tree_jet_leadingeta = 0;
   tree_jet_leadingeta2 = 0;
   tree_jet_jet_dR = 0;
   tree_jet_jet_dPhi = 0;
   tree_jet_jet_dEta = 0;
   tree_muon_jet_dRmin = 0;
   tree_muon_jet_dRmax = 0;
   tree_Hemi = 0;
   tree_Hemi_njet = 0;
   tree_Hemi_njet_nomu = 0;
   tree_Hemi_pt = 0;
   tree_Hemi_eta = 0;
   tree_Hemi_phi = 0;
   tree_Hemi_nTrks = 0;
   tree_Hemi_nTrks_sig = 0;
   tree_Hemi_nTrks_bad = 0;
   tree_Hemi_mass = 0;
   tree_HemiMu_mass = 0;
   tree_HemiMu_pt = 0;
   tree_HemiMu_dR = 0;
   tree_HemiMuOp_mass = 0;
   tree_HemiMuOp_pt = 0;
   tree_HemiMuOp_dR = 0;
   tree_Hemi_dR12 = 0;
   tree_Hemi_Vtx_step = 0;
   tree_Hemi_Vtx_isTight = 0;
   tree_Hemi_Vtx_NChi2 = 0;
   tree_Hemi_Vtx_nTrks = 0;
   tree_Hemi_Vtx_nTrks_sig = 0;
   tree_Hemi_Vtx_nTrks_bad = 0;
   tree_Hemi_Vtx_x = 0;
   tree_Hemi_Vtx_y = 0;
   tree_Hemi_Vtx_z = 0;
   tree_Hemi_Vtx_r = 0;
   tree_Hemi_Vtx_dR = 0;
   tree_Hemi_Vtx_SumtrackWeight = 0;
   tree_Hemi_Vtx_track_MeanDCA_d = 0;
   tree_Hemi_Vtx_Mass = 0;
   tree_Hemi_Vtx_dist = 0;
   tree_event_nVtx = 0;
   tree_event_Vtx_Vtx_dr = 0;
   tree_event_Vtx_Vtx_dz = 0;
   tree_event_Vtx_Vtx_dd = 0;
   tree_event_Vtx_Vtx_reldd = 0;
   tree_event_Vtx_Vtx_dR = 0;
   tree_event_Vtx_Vtx_step = 0;
   tree_Hemi_SecLLP = 0;
   tree_Hemi_LLP_SecVtx_dz = 0;
   tree_Hemi_LLP_SecVtx_dr = 0;
   tree_Hemi_SecLLP_ping = 0;
   tree_event_SecLLP_ping = 0;
   tree_Hemi_SecVtx = 0;
   tree_Hemi_SecVtx_step = 0;
   tree_Hemi_SecVtx_x = 0;
   tree_Hemi_SecVtx_y = 0;
   tree_Hemi_SecVtx_z = 0;
   tree_Hemi_SecVtx_r = 0;
   tree_Hemi_SecVtx_dR = 0;
   tree_Hemi_SecVtx_nTrks = 0;
   tree_Hemi_SecVtx_NChi2 = 0;
   tree_Hemi_SecVtx_dist = 0;
   tree_Hemi_SecVtx_track_MeanDCA_d = 0;
   tree_Hemi_SecVtx_SumtrackWeight = 0;
   tree_Hemi_SecVtx_Mass = 0;
   tree_event_MergedVtx_Vtx_dr = 0;
   tree_event_MergedVtx_Vtx_dz = 0;
   tree_event_MergedVtx_Vtx_dd = 0;
   tree_event_MergedVtx_Vtx_reldd = 0;
   tree_event_MergedVtx_Vtx_dR = 0;
   tree_event_MergedVtx_Vtx_step = 0;
   tree_Hemi_Vtx_BDT_nTrks = 0;
   tree_Hemi_Vtx_BDT_NChi2 = 0;
   tree_Hemi_Vtx_BDT_step = 0;
   tree_Hemi_Vtx_BDT_STW = 0;
   tree_Hemi_Vtx_BDT_Mass = 0;
   tree_Hemi_Vtx_BDT_HMass = 0;
   tree_Hemi_Vtx_BDT_ntrk10 = 0;
   tree_Hemi_Vtx_BDT_ntrk20 = 0;
   tree_Hemi_Vtx_BDT_MeanDCA = 0;
   tree_Hemi_Vtx_MVAval_Loose = 0;
   tree_Hemi_Vtx_MVAval_Tight = 0;
   tree_track_ipc = 0;
   tree_track_lost = 0;
   tree_track_px = 0;
   tree_track_py = 0;
   tree_track_pz = 0;
   tree_track_pt = 0;
   tree_track_eta = 0;
   tree_track_phi = 0;
   tree_track_charge = 0;
   tree_track_NChi2 = 0;
   tree_track_isHighPurity = 0;
   tree_track_dxy = 0;
   tree_track_dxyError = 0;
   tree_track_drSig = 0;
   tree_track_dz = 0;
   tree_track_dzError = 0;
   tree_track_dzSig = 0;
   tree_track_nHit = 0;
   tree_track_nHitPixel = 0;
   tree_track_nHitTIB = 0;
   tree_track_nHitTID = 0;
   tree_track_nHitTOB = 0;
   tree_track_nHitTEC = 0;
   tree_track_nHitPXB = 0;
   tree_track_nHitPXF = 0;
   tree_track_isHitPixel = 0;
   tree_track_nLayers = 0;
   tree_track_nLayersPixel = 0;
   tree_track_x = 0;
   tree_track_y = 0;
   tree_track_z = 0;
   tree_track_firstHit = 0;
   tree_track_region = 0;
   tree_track_firstHit_x = 0;
   tree_track_firstHit_y = 0;
   tree_track_firstHit_z = 0;
   tree_track_iJet = 0;
   tree_track_ntrk10 = 0;
   tree_track_ntrk20 = 0;
   tree_track_ntrk30 = 0;
   tree_track_ntrk40 = 0;
   tree_track_MVAval = 0;
   tree_track_Hemi_dR = 0;
   tree_track_Hemi_dRmax = 0;
   tree_K0_mass = 0;
   tree_K0_pt = 0;
   tree_L0_mass = 0;
   tree_L0_pt = 0;
   tree_V0_reco_mass = 0;
   tree_V0_reco_pt = 0;
   tree_V0_reco_source = 0;
   tree_SecInt_mass = 0;
   tree_SecInt_pt = 0;
   tree_SecInt_drSig = 0;
   tree_SecInt_dzSig = 0;
   tree_SecInt_layer = 0;
   tree_SecInt_selec = 0;
   tree_SecInt_r = 0;
   tree_SecInt_z = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("tree_LHE_Weights", &tree_LHE_Weights, &b_tree_LHE_Weights);
   fChain->SetBranchAddress("tree_MCEvt_weight", &tree_MCEvt_weight, &b_tree_MCEvt_weight);
   fChain->SetBranchAddress("tree_only_gen_wt", &tree_only_gen_wt, &b_tree_only_gen_wt);
   fChain->SetBranchAddress("tree_event_weight", &tree_event_weight, &b_tree_event_weight);
   fChain->SetBranchAddress("tree_genTop_Weight", &tree_genTop_Weight, &b_tree_genTop_Weight);
   fChain->SetBranchAddress("tree_gen_top_pt", &tree_gen_top_pt, &b_tree_gen_top_pt);
   fChain->SetBranchAddress("tree_gen_top_rw_pt", &tree_gen_top_rw_pt, &b_tree_gen_top_rw_pt);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("PUweight_Up", &PUweight_Up, &b_PUweight_Up);
   fChain->SetBranchAddress("PUweight_Down", &PUweight_Down, &b_PUweight_Down);
   fChain->SetBranchAddress("Prefweight", &Prefweight, &b_Prefweight);
   fChain->SetBranchAddress("Prefweight_Up", &Prefweight_Up, &b_Prefweight_Up);
   fChain->SetBranchAddress("Prefweight_Down", &Prefweight_Down, &b_Prefweight_Down);
   fChain->SetBranchAddress("PU_events", &PU_events, &b_PU_events);
   fChain->SetBranchAddress("tree_Filter", &tree_Filter, &b_tree_Filter);
   fChain->SetBranchAddress("tree_FilterSameSign", &tree_FilterSameSign, &b_tree_FilterSameSign);
   fChain->SetBranchAddress("tree_trigger_doublelepton", &tree_trigger_doublelepton, &b_tree_trigger_doublelepton);
   fChain->SetBranchAddress("tree_trigger_singlelepton", &tree_trigger_singlelepton, &b_tree_trigger_singlelepton);
   fChain->SetBranchAddress("tree_Good_PV", &tree_Good_PV, &b_tree_Good_PV);
   fChain->SetBranchAddress("tree_nPV", &tree_nPV, &b_tree_nPV);
   fChain->SetBranchAddress("tree_PV_x", &tree_PV_x, &b_tree_PV_x);
   fChain->SetBranchAddress("tree_PV_y", &tree_PV_y, &b_tree_PV_y);
   fChain->SetBranchAddress("tree_PV_z", &tree_PV_z, &b_tree_PV_z);
   fChain->SetBranchAddress("tree_PV_ez", &tree_PV_ez, &b_tree_PV_ez);
   fChain->SetBranchAddress("tree_PV_NChi2", &tree_PV_NChi2, &b_tree_PV_NChi2);
   fChain->SetBranchAddress("tree_PV_ndf", &tree_PV_ndf, &b_tree_PV_ndf);
   fChain->SetBranchAddress("tree_PFMet_et", &tree_PFMet_et, &b_tree_PFMet_et);
   fChain->SetBranchAddress("tree_PFMet_phi", &tree_PFMet_phi, &b_tree_PFMet_phi);
   fChain->SetBranchAddress("tree_HT", &tree_HT, &b_tree_HT);
   fChain->SetBranchAddress("tree_TRACK_SIZE", &tree_TRACK_SIZE, &b_tree_TRACK_SIZE);
   fChain->SetBranchAddress("tree_nTracks", &tree_nTracks, &b_tree_nTracks);
   fChain->SetBranchAddress("tree_nLostTracks", &tree_nLostTracks, &b_tree_nLostTracks);
   fChain->SetBranchAddress("tree_muon_GenRecoTriggerMatched", &tree_muon_GenRecoTriggerMatched, &b_tree_muon_GenRecoTriggerMatched);
   fChain->SetBranchAddress("tree_all_nmu", &tree_all_nmu, &b_tree_all_nmu);
   fChain->SetBranchAddress("tree_nmu", &tree_nmu, &b_tree_nmu);
   fChain->SetBranchAddress("tree_LT", &tree_LT, &b_tree_LT);
   fChain->SetBranchAddress("tree_Mmumu", &tree_Mmumu, &b_tree_Mmumu);
   fChain->SetBranchAddress("tree_MmumuSameSign", &tree_MmumuSameSign, &b_tree_MmumuSameSign);
   fChain->SetBranchAddress("tree_muon_isPrompt", &tree_muon_isPrompt, &b_tree_muon_isPrompt);
   fChain->SetBranchAddress("tree_muon_pt", &tree_muon_pt, &b_tree_muon_pt);
   fChain->SetBranchAddress("tree_muon_SF", &tree_muon_SF, &b_tree_muon_SF);
   fChain->SetBranchAddress("tree_muon_eta", &tree_muon_eta, &b_tree_muon_eta);
   fChain->SetBranchAddress("tree_muon_phi", &tree_muon_phi, &b_tree_muon_phi);
   fChain->SetBranchAddress("tree_muon_dxy", &tree_muon_dxy, &b_tree_muon_dxy);
   fChain->SetBranchAddress("tree_muon_dz", &tree_muon_dz, &b_tree_muon_dz);
   fChain->SetBranchAddress("tree_muon_charge", &tree_muon_charge, &b_tree_muon_charge);
   fChain->SetBranchAddress("tree_muon_correction", &tree_muon_correction, &b_tree_muon_correction);
   fChain->SetBranchAddress("tree_muon_gen", &tree_muon_gen, &b_tree_muon_gen);
   fChain->SetBranchAddress("tree_muon_dxyError", &tree_muon_dxyError, &b_tree_muon_dxyError);
   fChain->SetBranchAddress("tree_muon_dzError", &tree_muon_dzError, &b_tree_muon_dzError);
   fChain->SetBranchAddress("tree_muon_isLoose", &tree_muon_isLoose, &b_tree_muon_isLoose);
   fChain->SetBranchAddress("tree_muon_isMedium", &tree_muon_isMedium, &b_tree_muon_isMedium);
   fChain->SetBranchAddress("tree_muon_isTight", &tree_muon_isTight, &b_tree_muon_isTight);
   fChain->SetBranchAddress("tree_muon_isGlobal", &tree_muon_isGlobal, &b_tree_muon_isGlobal);
   fChain->SetBranchAddress("tree_muon_PFIsoVeryLoose", &tree_muon_PFIsoVeryLoose, &b_tree_muon_PFIsoVeryLoose);
   fChain->SetBranchAddress("tree_muon_PFIsoLoose", &tree_muon_PFIsoLoose, &b_tree_muon_PFIsoLoose);
   fChain->SetBranchAddress("tree_muon_PFIsoMedium", &tree_muon_PFIsoMedium, &b_tree_muon_PFIsoMedium);
   fChain->SetBranchAddress("tree_muon_PFIsoTight", &tree_muon_PFIsoTight, &b_tree_muon_PFIsoTight);
   fChain->SetBranchAddress("tree_muon_TkIsoLoose", &tree_muon_TkIsoLoose, &b_tree_muon_TkIsoLoose);
   fChain->SetBranchAddress("tree_muon_TkIsoTight", &tree_muon_TkIsoTight, &b_tree_muon_TkIsoTight);
   fChain->SetBranchAddress("tree_lepton_leadingpt", &tree_lepton_leadingpt, &b_tree_lepton_leadingpt);
   fChain->SetBranchAddress("tree_lepton_leadingpt2", &tree_lepton_leadingpt2, &b_tree_lepton_leadingpt2);
   fChain->SetBranchAddress("tree_lepton_leadingeta", &tree_lepton_leadingeta, &b_tree_lepton_leadingeta);
   fChain->SetBranchAddress("tree_lepton_leadingeta2", &tree_lepton_leadingeta2, &b_tree_lepton_leadingeta2);
   fChain->SetBranchAddress("tree_lepton_leadingphi", &tree_lepton_leadingphi, &b_tree_lepton_leadingphi);
   fChain->SetBranchAddress("tree_lepton_leadingphi2", &tree_lepton_leadingphi2, &b_tree_lepton_leadingphi2);
   fChain->SetBranchAddress("tree_lepton_lepton_dR", &tree_lepton_lepton_dR, &b_tree_lepton_lepton_dR);
   fChain->SetBranchAddress("tree_lepton_lepton_dPhi", &tree_lepton_lepton_dPhi, &b_tree_lepton_lepton_dPhi);
   fChain->SetBranchAddress("tree_lepton_lepton_dEta", &tree_lepton_lepton_dEta, &b_tree_lepton_lepton_dEta);
   fChain->SetBranchAddress("tree_lepton_leadingdxy", &tree_lepton_leadingdxy, &b_tree_lepton_leadingdxy);
   fChain->SetBranchAddress("tree_lepton_leadingdxy2", &tree_lepton_leadingdxy2, &b_tree_lepton_leadingdxy2);
   fChain->SetBranchAddress("tree_lepton_leadingdz", &tree_lepton_leadingdz, &b_tree_lepton_leadingdz);
   fChain->SetBranchAddress("tree_lepton_leadingdz2", &tree_lepton_leadingdz2, &b_tree_lepton_leadingdz2);
   fChain->SetBranchAddress("tree_all_nel", &tree_all_nel, &b_tree_all_nel);
   fChain->SetBranchAddress("tree_electron_nEle", &tree_electron_nEle, &b_tree_electron_nEle);
   fChain->SetBranchAddress("tree_electron_isPrompt", &tree_electron_isPrompt, &b_tree_electron_isPrompt);
   fChain->SetBranchAddress("tree_electron_pt", &tree_electron_pt, &b_tree_electron_pt);
   fChain->SetBranchAddress("tree_electron_eta", &tree_electron_eta, &b_tree_electron_eta);
   fChain->SetBranchAddress("tree_electron_phi", &tree_electron_phi, &b_tree_electron_phi);
   fChain->SetBranchAddress("tree_electron_charge", &tree_electron_charge, &b_tree_electron_charge);
   fChain->SetBranchAddress("tree_electron_dxy", &tree_electron_dxy, &b_tree_electron_dxy);
   fChain->SetBranchAddress("tree_electron_dz", &tree_electron_dz, &b_tree_electron_dz);
   fChain->SetBranchAddress("tree_electron_gen", &tree_electron_gen, &b_tree_electron_gen);
   fChain->SetBranchAddress("tree_electron_energy", &tree_electron_energy, &b_tree_electron_energy);
   fChain->SetBranchAddress("tree_electron_et", &tree_electron_et, &b_tree_electron_et);
   fChain->SetBranchAddress("tree_electron_ecal_trk_postcorr", &tree_electron_ecal_trk_postcorr, &b_tree_electron_ecal_trk_postcorr);
   fChain->SetBranchAddress("tree_electron_isoR4", &tree_electron_isoR4, &b_tree_electron_isoR4);
   fChain->SetBranchAddress("tree_electron_IsLoose", &tree_electron_IsLoose, &b_tree_electron_IsLoose);
   fChain->SetBranchAddress("tree_electron_IsMedium", &tree_electron_IsMedium, &b_tree_electron_IsMedium);
   fChain->SetBranchAddress("tree_electron_IsTight", &tree_electron_IsTight, &b_tree_electron_IsTight);
   fChain->SetBranchAddress("tree_njet", &tree_njet, &b_tree_njet);
   fChain->SetBranchAddress("tree_njetNOmu", &tree_njetNOmu, &b_tree_njetNOmu);
   fChain->SetBranchAddress("tree_jet_pt", &tree_jet_pt, &b_tree_jet_pt);
   fChain->SetBranchAddress("tree_jet_eta", &tree_jet_eta, &b_tree_jet_eta);
   fChain->SetBranchAddress("tree_jet_phi", &tree_jet_phi, &b_tree_jet_phi);
   fChain->SetBranchAddress("tree_jet_HadronFlavour", &tree_jet_HadronFlavour, &b_tree_jet_HadronFlavour);
   fChain->SetBranchAddress("tree_jet_btag_DeepJet", &tree_jet_btag_DeepJet, &b_tree_jet_btag_DeepJet);
   fChain->SetBranchAddress("tree_jet_E", &tree_jet_E, &b_tree_jet_E);
   fChain->SetBranchAddress("tree_jet_leadingpt", &tree_jet_leadingpt, &b_tree_jet_leadingpt);
   fChain->SetBranchAddress("tree_jet_leadingpt2", &tree_jet_leadingpt2, &b_tree_jet_leadingpt2);
   fChain->SetBranchAddress("tree_jet_leadingeta", &tree_jet_leadingeta, &b_tree_jet_leadingeta);
   fChain->SetBranchAddress("tree_jet_leadingeta2", &tree_jet_leadingeta2, &b_tree_jet_leadingeta2);
   fChain->SetBranchAddress("tree_jet_jet_dR", &tree_jet_jet_dR, &b_tree_jet_jet_dR);
   fChain->SetBranchAddress("tree_jet_jet_dPhi", &tree_jet_jet_dPhi, &b_tree_jet_jet_dPhi);
   fChain->SetBranchAddress("tree_jet_jet_dEta", &tree_jet_jet_dEta, &b_tree_jet_jet_dEta);
   fChain->SetBranchAddress("tree_muon_jet_dRmin", &tree_muon_jet_dRmin, &b_tree_muon_jet_dRmin);
   fChain->SetBranchAddress("tree_muon_jet_dRmax", &tree_muon_jet_dRmax, &b_tree_muon_jet_dRmax);
   fChain->SetBranchAddress("tree_Evts_MVAval", &tree_Evts_MVAval, &b_tree_Evts_MVAval);
   fChain->SetBranchAddress("tree_Evts_MVAvalDY", &tree_Evts_MVAvalDY, &b_tree_Evts_MVAvalDY);
   fChain->SetBranchAddress("tree_Evts_MVAvalTT", &tree_Evts_MVAvalTT, &b_tree_Evts_MVAvalTT);
   fChain->SetBranchAddress("tree_Hemi", &tree_Hemi, &b_tree_Hemi);
   fChain->SetBranchAddress("tree_Hemi_njet", &tree_Hemi_njet, &b_tree_Hemi_njet);
   fChain->SetBranchAddress("tree_Hemi_njet_nomu", &tree_Hemi_njet_nomu, &b_tree_Hemi_njet_nomu);
   fChain->SetBranchAddress("tree_Hemi_pt", &tree_Hemi_pt, &b_tree_Hemi_pt);
   fChain->SetBranchAddress("tree_Hemi_eta", &tree_Hemi_eta, &b_tree_Hemi_eta);
   fChain->SetBranchAddress("tree_Hemi_phi", &tree_Hemi_phi, &b_tree_Hemi_phi);
   fChain->SetBranchAddress("tree_Hemi_nTrks", &tree_Hemi_nTrks, &b_tree_Hemi_nTrks);
   fChain->SetBranchAddress("tree_Hemi_nTrks_sig", &tree_Hemi_nTrks_sig, &b_tree_Hemi_nTrks_sig);
   fChain->SetBranchAddress("tree_Hemi_nTrks_bad", &tree_Hemi_nTrks_bad, &b_tree_Hemi_nTrks_bad);
   fChain->SetBranchAddress("tree_Hemi_mass", &tree_Hemi_mass, &b_tree_Hemi_mass);
   fChain->SetBranchAddress("tree_HemiMu_mass", &tree_HemiMu_mass, &b_tree_HemiMu_mass);
   fChain->SetBranchAddress("tree_HemiMu_pt", &tree_HemiMu_pt, &b_tree_HemiMu_pt);
   fChain->SetBranchAddress("tree_HemiMu_dR", &tree_HemiMu_dR, &b_tree_HemiMu_dR);
   fChain->SetBranchAddress("tree_HemiMuOp_mass", &tree_HemiMuOp_mass, &b_tree_HemiMuOp_mass);
   fChain->SetBranchAddress("tree_HemiMuOp_pt", &tree_HemiMuOp_pt, &b_tree_HemiMuOp_pt);
   fChain->SetBranchAddress("tree_HemiMuOp_dR", &tree_HemiMuOp_dR, &b_tree_HemiMuOp_dR);
   fChain->SetBranchAddress("tree_Hemi_dR12", &tree_Hemi_dR12, &b_tree_Hemi_dR12);
   fChain->SetBranchAddress("tree_Hemi_Vtx_step", &tree_Hemi_Vtx_step, &b_tree_Hemi_Vtx_step);
   fChain->SetBranchAddress("tree_Hemi_Vtx_isTight", &tree_Hemi_Vtx_isTight, &b_tree_Hemi_Vtx_isTight);
   fChain->SetBranchAddress("tree_Hemi_Vtx_NChi2", &tree_Hemi_Vtx_NChi2, &b_tree_Hemi_Vtx_NChi2);
   fChain->SetBranchAddress("tree_Hemi_Vtx_nTrks", &tree_Hemi_Vtx_nTrks, &b_tree_Hemi_Vtx_nTrks);
   fChain->SetBranchAddress("tree_Hemi_Vtx_nTrks_sig", &tree_Hemi_Vtx_nTrks_sig, &b_tree_Hemi_Vtx_nTrks_sig);
   fChain->SetBranchAddress("tree_Hemi_Vtx_nTrks_bad", &tree_Hemi_Vtx_nTrks_bad, &b_tree_Hemi_Vtx_nTrks_bad);
   fChain->SetBranchAddress("tree_Hemi_Vtx_x", &tree_Hemi_Vtx_x, &b_tree_Hemi_Vtx_x);
   fChain->SetBranchAddress("tree_Hemi_Vtx_y", &tree_Hemi_Vtx_y, &b_tree_Hemi_Vtx_y);
   fChain->SetBranchAddress("tree_Hemi_Vtx_z", &tree_Hemi_Vtx_z, &b_tree_Hemi_Vtx_z);
   fChain->SetBranchAddress("tree_Hemi_Vtx_r", &tree_Hemi_Vtx_r, &b_tree_Hemi_Vtx_r);
   fChain->SetBranchAddress("tree_Hemi_Vtx_dR", &tree_Hemi_Vtx_dR, &b_tree_Hemi_Vtx_dR);
   fChain->SetBranchAddress("tree_Hemi_Vtx_SumtrackWeight", &tree_Hemi_Vtx_SumtrackWeight, &b_tree_Hemi_Vtx_SumtrackWeight);
   fChain->SetBranchAddress("tree_Hemi_Vtx_track_MeanDCA_d", &tree_Hemi_Vtx_track_MeanDCA_d, &b_tree_Hemi_Vtx_track_MeanDCA_d);
   fChain->SetBranchAddress("tree_Hemi_Vtx_Mass", &tree_Hemi_Vtx_Mass, &b_tree_Hemi_Vtx_Mass);
   fChain->SetBranchAddress("tree_Hemi_Vtx_dist", &tree_Hemi_Vtx_dist, &b_tree_Hemi_Vtx_dist);
   fChain->SetBranchAddress("tree_event_nVtx", &tree_event_nVtx, &b_tree_event_nVtx);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_dr", &tree_event_Vtx_Vtx_dr, &b_tree_event_Vtx_Vtx_dr);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_dz", &tree_event_Vtx_Vtx_dz, &b_tree_event_Vtx_Vtx_dz);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_dd", &tree_event_Vtx_Vtx_dd, &b_tree_event_Vtx_Vtx_dd);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_reldd", &tree_event_Vtx_Vtx_reldd, &b_tree_event_Vtx_Vtx_reldd);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_dR", &tree_event_Vtx_Vtx_dR, &b_tree_event_Vtx_Vtx_dR);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_step", &tree_event_Vtx_Vtx_step, &b_tree_event_Vtx_Vtx_step);
   fChain->SetBranchAddress("tree_Hemi_SecLLP", &tree_Hemi_SecLLP, &b_tree_Hemi_SecLLP);
   fChain->SetBranchAddress("tree_Hemi_LLP_SecVtx_dz", &tree_Hemi_LLP_SecVtx_dz, &b_tree_Hemi_LLP_SecVtx_dz);
   fChain->SetBranchAddress("tree_Hemi_LLP_SecVtx_dr", &tree_Hemi_LLP_SecVtx_dr, &b_tree_Hemi_LLP_SecVtx_dr);
   fChain->SetBranchAddress("tree_Hemi_SecLLP_ping", &tree_Hemi_SecLLP_ping, &b_tree_Hemi_SecLLP_ping);
   fChain->SetBranchAddress("tree_event_SecLLP_ping", &tree_event_SecLLP_ping, &b_tree_event_SecLLP_ping);
   fChain->SetBranchAddress("tree_Hemi_SecVtx", &tree_Hemi_SecVtx, &b_tree_Hemi_SecVtx);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_step", &tree_Hemi_SecVtx_step, &b_tree_Hemi_SecVtx_step);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_x", &tree_Hemi_SecVtx_x, &b_tree_Hemi_SecVtx_x);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_y", &tree_Hemi_SecVtx_y, &b_tree_Hemi_SecVtx_y);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_z", &tree_Hemi_SecVtx_z, &b_tree_Hemi_SecVtx_z);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_r", &tree_Hemi_SecVtx_r, &b_tree_Hemi_SecVtx_r);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_dR", &tree_Hemi_SecVtx_dR, &b_tree_Hemi_SecVtx_dR);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_nTrks", &tree_Hemi_SecVtx_nTrks, &b_tree_Hemi_SecVtx_nTrks);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_NChi2", &tree_Hemi_SecVtx_NChi2, &b_tree_Hemi_SecVtx_NChi2);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_dist", &tree_Hemi_SecVtx_dist, &b_tree_Hemi_SecVtx_dist);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_track_MeanDCA_d", &tree_Hemi_SecVtx_track_MeanDCA_d, &b_tree_Hemi_SecVtx_track_MeanDCA_d);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_SumtrackWeight", &tree_Hemi_SecVtx_SumtrackWeight, &b_tree_Hemi_SecVtx_SumtrackWeight);
   fChain->SetBranchAddress("tree_Hemi_SecVtx_Mass", &tree_Hemi_SecVtx_Mass, &b_tree_Hemi_SecVtx_Mass);
   fChain->SetBranchAddress("tree_event_MergedVtx_Vtx_dr", &tree_event_MergedVtx_Vtx_dr, &b_tree_event_MergedVtx_Vtx_dr);
   fChain->SetBranchAddress("tree_event_MergedVtx_Vtx_dz", &tree_event_MergedVtx_Vtx_dz, &b_tree_event_MergedVtx_Vtx_dz);
   fChain->SetBranchAddress("tree_event_MergedVtx_Vtx_dd", &tree_event_MergedVtx_Vtx_dd, &b_tree_event_MergedVtx_Vtx_dd);
   fChain->SetBranchAddress("tree_event_MergedVtx_Vtx_reldd", &tree_event_MergedVtx_Vtx_reldd, &b_tree_event_MergedVtx_Vtx_reldd);
   fChain->SetBranchAddress("tree_event_MergedVtx_Vtx_dR", &tree_event_MergedVtx_Vtx_dR, &b_tree_event_MergedVtx_Vtx_dR);
   fChain->SetBranchAddress("tree_event_MergedVtx_Vtx_step", &tree_event_MergedVtx_Vtx_step, &b_tree_event_MergedVtx_Vtx_step);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BDT_nTrks", &tree_Hemi_Vtx_BDT_nTrks, &b_tree_Hemi_Vtx_BDT_nTrks);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BDT_NChi2", &tree_Hemi_Vtx_BDT_NChi2, &b_tree_Hemi_Vtx_BDT_NChi2);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BDT_step", &tree_Hemi_Vtx_BDT_step, &b_tree_Hemi_Vtx_BDT_step);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BDT_STW", &tree_Hemi_Vtx_BDT_STW, &b_tree_Hemi_Vtx_BDT_STW);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BDT_Mass", &tree_Hemi_Vtx_BDT_Mass, &b_tree_Hemi_Vtx_BDT_Mass);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BDT_HMass", &tree_Hemi_Vtx_BDT_HMass, &b_tree_Hemi_Vtx_BDT_HMass);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BDT_ntrk10", &tree_Hemi_Vtx_BDT_ntrk10, &b_tree_Hemi_Vtx_BDT_ntrk10);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BDT_ntrk20", &tree_Hemi_Vtx_BDT_ntrk20, &b_tree_Hemi_Vtx_BDT_ntrk20);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BDT_MeanDCA", &tree_Hemi_Vtx_BDT_MeanDCA, &b_tree_Hemi_Vtx_BDT_MeanDCA);
   fChain->SetBranchAddress("tree_Hemi_Vtx_MVAval_Loose", &tree_Hemi_Vtx_MVAval_Loose, &b_tree_Hemi_Vtx_MVAval_Loose);
   fChain->SetBranchAddress("tree_Hemi_Vtx_MVAval_Tight", &tree_Hemi_Vtx_MVAval_Tight, &b_tree_Hemi_Vtx_MVAval_Tight);
   fChain->SetBranchAddress("tree_track_ipc", &tree_track_ipc, &b_tree_track_ipc);
   fChain->SetBranchAddress("tree_track_lost", &tree_track_lost, &b_tree_track_lost);
   fChain->SetBranchAddress("tree_track_px", &tree_track_px, &b_tree_track_px);
   fChain->SetBranchAddress("tree_track_py", &tree_track_py, &b_tree_track_py);
   fChain->SetBranchAddress("tree_track_pz", &tree_track_pz, &b_tree_track_pz);
   fChain->SetBranchAddress("tree_track_pt", &tree_track_pt, &b_tree_track_pt);
   fChain->SetBranchAddress("tree_track_eta", &tree_track_eta, &b_tree_track_eta);
   fChain->SetBranchAddress("tree_track_phi", &tree_track_phi, &b_tree_track_phi);
   fChain->SetBranchAddress("tree_track_charge", &tree_track_charge, &b_tree_track_charge);
   fChain->SetBranchAddress("tree_track_NChi2", &tree_track_NChi2, &b_tree_track_NChi2);
   fChain->SetBranchAddress("tree_track_isHighPurity", &tree_track_isHighPurity, &b_tree_track_isHighPurity);
   fChain->SetBranchAddress("tree_track_dxy", &tree_track_dxy, &b_tree_track_dxy);
   fChain->SetBranchAddress("tree_track_dxyError", &tree_track_dxyError, &b_tree_track_dxyError);
   fChain->SetBranchAddress("tree_track_drSig", &tree_track_drSig, &b_tree_track_drSig);
   fChain->SetBranchAddress("tree_track_dz", &tree_track_dz, &b_tree_track_dz);
   fChain->SetBranchAddress("tree_track_dzError", &tree_track_dzError, &b_tree_track_dzError);
   fChain->SetBranchAddress("tree_track_dzSig", &tree_track_dzSig, &b_tree_track_dzSig);
   fChain->SetBranchAddress("tree_track_nHit", &tree_track_nHit, &b_tree_track_nHit);
   fChain->SetBranchAddress("tree_track_nHitPixel", &tree_track_nHitPixel, &b_tree_track_nHitPixel);
   fChain->SetBranchAddress("tree_track_nHitTIB", &tree_track_nHitTIB, &b_tree_track_nHitTIB);
   fChain->SetBranchAddress("tree_track_nHitTID", &tree_track_nHitTID, &b_tree_track_nHitTID);
   fChain->SetBranchAddress("tree_track_nHitTOB", &tree_track_nHitTOB, &b_tree_track_nHitTOB);
   fChain->SetBranchAddress("tree_track_nHitTEC", &tree_track_nHitTEC, &b_tree_track_nHitTEC);
   fChain->SetBranchAddress("tree_track_nHitPXB", &tree_track_nHitPXB, &b_tree_track_nHitPXB);
   fChain->SetBranchAddress("tree_track_nHitPXF", &tree_track_nHitPXF, &b_tree_track_nHitPXF);
   fChain->SetBranchAddress("tree_track_isHitPixel", &tree_track_isHitPixel, &b_tree_track_isHitPixel);
   fChain->SetBranchAddress("tree_track_nLayers", &tree_track_nLayers, &b_tree_track_nLayers);
   fChain->SetBranchAddress("tree_track_nLayersPixel", &tree_track_nLayersPixel, &b_tree_track_nLayersPixel);
   fChain->SetBranchAddress("tree_track_x", &tree_track_x, &b_tree_track_x);
   fChain->SetBranchAddress("tree_track_y", &tree_track_y, &b_tree_track_y);
   fChain->SetBranchAddress("tree_track_z", &tree_track_z, &b_tree_track_z);
   fChain->SetBranchAddress("tree_track_firstHit", &tree_track_firstHit, &b_tree_track_firstHit);
   fChain->SetBranchAddress("tree_track_region", &tree_track_region, &b_tree_track_region);
   fChain->SetBranchAddress("tree_track_firstHit_x", &tree_track_firstHit_x, &b_tree_track_firstHit_x);
   fChain->SetBranchAddress("tree_track_firstHit_y", &tree_track_firstHit_y, &b_tree_track_firstHit_y);
   fChain->SetBranchAddress("tree_track_firstHit_z", &tree_track_firstHit_z, &b_tree_track_firstHit_z);
   fChain->SetBranchAddress("tree_track_iJet", &tree_track_iJet, &b_tree_track_iJet);
   fChain->SetBranchAddress("tree_track_ntrk10", &tree_track_ntrk10, &b_tree_track_ntrk10);
   fChain->SetBranchAddress("tree_track_ntrk20", &tree_track_ntrk20, &b_tree_track_ntrk20);
   fChain->SetBranchAddress("tree_track_ntrk30", &tree_track_ntrk30, &b_tree_track_ntrk30);
   fChain->SetBranchAddress("tree_track_ntrk40", &tree_track_ntrk40, &b_tree_track_ntrk40);
   fChain->SetBranchAddress("tree_track_MVAval", &tree_track_MVAval, &b_tree_track_MVAval);
   fChain->SetBranchAddress("tree_track_Hemi_dR", &tree_track_Hemi_dR, &b_tree_track_Hemi_dR);
   fChain->SetBranchAddress("tree_track_Hemi_dRmax", &tree_track_Hemi_dRmax, &b_tree_track_Hemi_dRmax);
   fChain->SetBranchAddress("tree_K0_mass", &tree_K0_mass, &b_tree_K0_mass);
   fChain->SetBranchAddress("tree_K0_pt", &tree_K0_pt, &b_tree_K0_pt);
   fChain->SetBranchAddress("tree_L0_mass", &tree_L0_mass, &b_tree_L0_mass);
   fChain->SetBranchAddress("tree_L0_pt", &tree_L0_pt, &b_tree_L0_pt);
   fChain->SetBranchAddress("tree_V0_reco_mass", &tree_V0_reco_mass, &b_tree_V0_reco_mass);
   fChain->SetBranchAddress("tree_V0_reco_pt", &tree_V0_reco_pt, &b_tree_V0_reco_pt);
   fChain->SetBranchAddress("tree_V0_reco_source", &tree_V0_reco_source, &b_tree_V0_reco_source);
   fChain->SetBranchAddress("tree_SecInt_mass", &tree_SecInt_mass, &b_tree_SecInt_mass);
   fChain->SetBranchAddress("tree_SecInt_pt", &tree_SecInt_pt, &b_tree_SecInt_pt);
   fChain->SetBranchAddress("tree_SecInt_drSig", &tree_SecInt_drSig, &b_tree_SecInt_drSig);
   fChain->SetBranchAddress("tree_SecInt_dzSig", &tree_SecInt_dzSig, &b_tree_SecInt_dzSig);
   fChain->SetBranchAddress("tree_SecInt_layer", &tree_SecInt_layer, &b_tree_SecInt_layer);
   fChain->SetBranchAddress("tree_SecInt_selec", &tree_SecInt_selec, &b_tree_SecInt_selec);
   fChain->SetBranchAddress("tree_SecInt_r", &tree_SecInt_r, &b_tree_SecInt_r);
   fChain->SetBranchAddress("tree_SecInt_z", &tree_SecInt_z, &b_tree_SecInt_z);
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
