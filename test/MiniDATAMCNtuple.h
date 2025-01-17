//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue May 14 11:25:59 2024 by ROOT version 6.14/09
// from TTree ttree/ttree
// found on file: Ntuple_cms.root
//////////////////////////////////////////////////////////

#ifndef MiniDATAMCNtuple_h
#define MiniDATAMCNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include <TMath.h>





class MiniDATAMCNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiBlock;
   Double_t        tree_only_gen_wt;
   Double_t        tree_event_weight;
   Double_t        tree_genTop_Weight;
   Double_t        PUweight;
   Double_t        PUweight_Up;
   Double_t        PUweight_Down;
   Double_t        Prefweight;
   Int_t           PU_events;
   Bool_t          tree_only_tigger_filter;
   Bool_t          tree_trigger_doublelepton;
   Bool_t          tree_trigger_singlelepton;
   Bool_t          tree_Filter;
   Bool_t          tree_FilterSameSign;
   Bool_t          tree_Good_PV;
   Int_t           tree_muon_GenRecoTriggerMatched;
   Float_t         tree_Evts_MVAval;
   Float_t         tree_Evts_MVAvalDY;
   Float_t         tree_Evts_MVAvalTT;
   Float_t         tree_bs_PosX;
   Float_t         tree_bs_PosY;
   Float_t         tree_bs_PosZ;
   Int_t           tree_nPV;
   Float_t         tree_PV_x;
   Float_t         tree_PV_y;
   Float_t         tree_PV_z;
   Float_t         tree_PV_ez;
   Float_t         tree_PV_NChi2;
   Int_t           tree_PV_ndf;
   Float_t         tree_PV_rho;
   Int_t           tree_all_nmu;
   Int_t           tree_nmu;
   Float_t         tree_LT;
   Float_t         tree_Mmumu;
   Float_t         tree_MmumuSameSign;
   std::vector<bool>    *tree_muon_isPrompt;
   std::vector<float>   *tree_muon_pt;
   std::vector<float>   *tree_muon_SF;
   std::vector<float>   *tree_muon_eta;
   std::vector<float>   *tree_muon_phi;
   std::vector<float>   *tree_muon_x;
   std::vector<float>   *tree_muon_y;
   std::vector<float>   *tree_muon_z;
   std::vector<float>   *tree_muon_dxy;
   std::vector<float>   *tree_muon_dxyError;
   std::vector<float>   *tree_muon_dz;
   std::vector<float>   *tree_muon_dzError;
   std::vector<int>     *tree_muon_charge;
   std::vector<bool>    *tree_muon_isLoose;
   std::vector<bool>    *tree_muon_isMedium;
   std::vector<bool>    *tree_muon_isTight;
   std::vector<bool>    *tree_muon_isGlobal;
   std::vector<float>   *tree_muon_isoR3;
   std::vector<bool>    *tree_muon_trigger_dimu;
   std::vector<bool>    *tree_muon_trigger_isomu;
   std::vector<bool>    *tree_muon_PFIsoVeryLoose;
   std::vector<bool>    *tree_muon_PFIsoLoose;
   std::vector<bool>    *tree_muon_PFIsoMedium;
   std::vector<bool>    *tree_muon_PFIsoTight;
   std::vector<bool>    *tree_muon_MiniIsoLoose;
   std::vector<bool>    *tree_muon_MiniIsoMedium;
   std::vector<bool>    *tree_muon_MiniIsoTight;
   std::vector<bool>    *tree_muon_TkIsoLoose;
   std::vector<bool>    *tree_muon_TkIsoTight;
   std::vector<float>   *tree_muon_trkLayers;
   std::vector<float>   *tree_muon_miniIso;
   std::vector<float>   *tree_muon_correction;
   std::vector<int>     *tree_muon_gen;
   std::vector<float>   *tree_reco_muon_leadingpt;
   std::vector<float>   *tree_reco_electron_leadingpt2;
   std::vector<float>   *tree_reco_muon_leadingeta;
   std::vector<float>   *tree_reco_electron_leadingeta2;
   std::vector<float>   *tree_reco_muon_leadingphi;
   std::vector<float>   *tree_reco_electron_leadingphi2;
   std::vector<float>   *tree_trig_muon_leadingpt;
   std::vector<float>   *tree_trig_electron_leadingpt2;
   std::vector<float>   *tree_trig_muon_leadingeta;
   std::vector<float>   *tree_trig_electron_leadingeta2;
   std::vector<float>   *tree_trig_muon_leadingphi;
   std::vector<float>   *tree_trig_electron_leadingphi2;

   std::vector<float>   *tree_lepton_b4trigger_leadingpt;
   std::vector<float>   *tree_lepton_b4trigger_leadingpt2;


   std::vector<float>   *tree_reco_lepton_leadingpt;
   std::vector<float>   *tree_reco_lepton_leadingpt2;
   std::vector<float>   *tree_reco_lepton_leadingeta;
   std::vector<float>   *tree_reco_lepton_leadingeta2;
   std::vector<float>   *tree_reco_lepton_leadingphi;
   std::vector<float>   *tree_reco_lepton_leadingphi2;
   std::vector<float>   *tree_trig_lepton_leadingpt;
   std::vector<float>   *tree_trig_lepton_leadingpt2;
   std::vector<float>   *tree_trig_lepton_leadingeta;
   std::vector<float>   *tree_trig_lepton_leadingeta2;
   std::vector<float>   *tree_trig_lepton_leadingphi;
   std::vector<float>   *tree_trig_lepton_leadingphi2;
   std::vector<float>   *tree_lepton_leadingpt;
   std::vector<float>   *tree_lepton_leadingpt2;
   std::vector<float>   *tree_lepton_leadingeta;
   std::vector<float>   *tree_lepton_leadingeta2;
   std::vector<float>   *tree_lepton_leadingphi;
   std::vector<float>   *tree_lepton_leadingphi2;
   std::vector<float>   *tree_lepton_leadingdxy;
   std::vector<float>   *tree_lepton_leadingdxy2;
   std::vector<float>   *tree_lepton_leadingdz;
   std::vector<float>   *tree_lepton_leadingdz2;
   std::vector<float>   *tree_lepton_lepton_dR;
   std::vector<float>   *tree_lepton_lepton_dPhi;
   std::vector<float>   *tree_lepton_lepton_dEta;
   std::vector<float>   *tree_ll_pt;
   std::vector<float>   *tree_ll_eta;
   std::vector<float>   *tree_ll_phi;
   std::vector<float>   *tree_ll_px;
   std::vector<float>   *tree_ll_py;
   std::vector<float>   *tree_ll_pz;
   std::vector<float>   *tree_ll_energy;
   std::vector<float>   *tree_ll_mass;
   Int_t           tree_all_nel;
   Int_t           tree_electron_nEle;
   std::vector<bool>    *tree_electron_isPrompt;
   std::vector<float>   *tree_electron_pt;
   std::vector<float>   *tree_electron_eta;
   std::vector<float>   *tree_electron_phi;
   std::vector<float>   *tree_electron_x;
   std::vector<float>   *tree_electron_y;
   std::vector<float>   *tree_electron_z;
   std::vector<float>   *tree_electron_energy;
   std::vector<float>   *tree_electron_et;
   std::vector<float>   *tree_electron_ecal_trk_postcorr;
   std::vector<int>     *tree_electron_charge;
   std::vector<float>   *tree_electron_isoR4;
   std::vector<bool>    *tree_electron_IsLoose;
   std::vector<bool>    *tree_electron_IsMedium;
   std::vector<bool>    *tree_electron_IsTight;
   std::vector<float>   *tree_electron_dxy;
   std::vector<float>   *tree_electron_dz;
   std::vector<int>     *tree_electron_gen;
   Float_t         tree_PFMet_et;
   Float_t         tree_PFMet_phi;
   Float_t         tree_PFMet_sig;
   Float_t         tree_PFMet_pt;
   Int_t           tree_njet;
   Int_t           tree_njetNOmu;
   std::vector<float>   *tree_jet_pt;
   std::vector<float>   *tree_jet_eta;
   std::vector<float>   *tree_jet_phi;
   std::vector<float>   *tree_jet_px;
   std::vector<float>   *tree_jet_py;
   std::vector<float>   *tree_jet_pz;
   std::vector<float>   *tree_jet_E;
   std::vector<bool>    *tree_jet_tightid_LepVeto;
   std::vector<bool>    *tree_jet_tightid;
   std::vector<bool>    *tree_jet_TightJetIDLepVeto;
   std::vector<bool>    *tree_jet_TightJetID;
   std::vector<float>   *tree_jet_HadronFlavour;
   std::vector<float>   *tree_jet_pileupID;
   std::vector<float>   *tree_jet_btag_DeepCSV;
   std::vector<float>   *tree_jet_btag_DeepJet;
   std::vector<float>   *tree_jet_leadingpt;
   std::vector<float>   *tree_jet_leadingpt2;
   std::vector<float>   *tree_jet_leadingeta;
   std::vector<float>   *tree_jet_leadingeta2;
   std::vector<float>   *tree_jet_leadingMuon_dR;
   std::vector<float>   *tree_jet_leadingMuon2_dR;
   std::vector<float>   *tree_jet_jet_dR;
   std::vector<float>   *tree_jet_jet_dPhi;
   std::vector<float>   *tree_jet_jet_dEta;
   std::vector<float>   *tree_muon_jet_dRmin;
   std::vector<float>   *tree_muon_jet_dRmax;
   std::vector<float>   *tree_elemu_jet_dRmin;
   std::vector<float>   *tree_elemu_jet_dRmax;
   std::vector<float>   *tree_ele_jet_dRmin;
   std::vector<float>   *tree_ele_jet_dRmax;
   Float_t         tree_HT;
   Int_t           tree_nK0;
   std::vector<float>   *tree_K0_x;
   std::vector<float>   *tree_K0_y;
   std::vector<float>   *tree_K0_z;
   std::vector<float>   *tree_K0_r;
   std::vector<float>   *tree_K0_NChi2;
   std::vector<float>   *tree_K0_ndf;
   std::vector<float>   *tree_K0_mass;
   std::vector<float>   *tree_K0_pt;
   std::vector<float>   *tree_K0_eta;
   std::vector<float>   *tree_K0_phi;
   std::vector<unsigned int> *tree_K0_nDaughters;
   Int_t           tree_nLambda;
   std::vector<float>   *tree_L0_x;
   std::vector<float>   *tree_L0_y;
   std::vector<float>   *tree_L0_z;
   std::vector<float>   *tree_L0_r;
   std::vector<float>   *tree_L0_NChi2;
   std::vector<float>   *tree_L0_ndf;
   std::vector<unsigned int> *tree_L0_nDaughters;
   std::vector<float>   *tree_L0_mass;
   std::vector<float>   *tree_L0_pt;
   std::vector<float>   *tree_L0_eta;
   std::vector<float>   *tree_L0_phi;
   Int_t           tree_nV0_reco;
   std::vector<float>   *tree_V0_reco_x;
   std::vector<float>   *tree_V0_reco_y;
   std::vector<float>   *tree_V0_reco_z;
   std::vector<float>   *tree_V0_reco_r;
   std::vector<float>   *tree_V0_reco_drSig;
   std::vector<float>   *tree_V0_reco_dzSig;
   std::vector<float>   *tree_V0_reco_angleXY;
   std::vector<float>   *tree_V0_reco_angleZ;
   std::vector<float>   *tree_V0_reco_NChi2;
   std::vector<float>   *tree_V0_reco_ndf;
   std::vector<float>   *tree_V0_reco_mass;
   std::vector<float>   *tree_V0_reco_pt;
   std::vector<float>   *tree_V0_reco_eta;
   std::vector<float>   *tree_V0_reco_phi;
   std::vector<int>     *tree_V0_reco_source;
   std::vector<bool>    *tree_V0_reco_badTkHit;
   std::vector<float>   *tree_V0_reco_dca;
   Int_t           tree_nSecInt;
   std::vector<float>   *tree_SecInt_x;
   std::vector<float>   *tree_SecInt_y;
   std::vector<float>   *tree_SecInt_z;
   std::vector<float>   *tree_SecInt_r;
   std::vector<float>   *tree_SecInt_d;
   std::vector<float>   *tree_SecInt_drSig;
   std::vector<float>   *tree_SecInt_dzSig;
   std::vector<float>   *tree_SecInt_angleXY;
   std::vector<float>   *tree_SecInt_angleZ;
   std::vector<float>   *tree_SecInt_NChi2;
   std::vector<float>   *tree_SecInt_ndf;
   std::vector<float>   *tree_SecInt_mass;
   std::vector<float>   *tree_SecInt_pt;
   std::vector<float>   *tree_SecInt_eta;
   std::vector<float>   *tree_SecInt_phi;
   std::vector<int>     *tree_SecInt_charge;
   std::vector<bool>    *tree_SecInt_badTkHit;
   std::vector<float>   *tree_SecInt_dca;
   std::vector<bool>    *tree_SecInt_selec;
   std::vector<int>     *tree_SecInt_layer;
   std::vector<int>     *tree_SecInt_LLP;
   std::vector<float>   *tree_SecInt_LLP_dr;
   std::vector<float>   *tree_SecInt_LLP_dz;
   std::vector<float>   *tree_SecInt_LLP_dd;
   std::vector<unsigned int> *tree_SecInt_tk1;
   std::vector<unsigned int> *tree_SecInt_tk2;
   Int_t           tree_nYConv;
   std::vector<float>   *tree_Yc_x;
   std::vector<float>   *tree_Yc_y;
   std::vector<float>   *tree_Yc_z;
   std::vector<float>   *tree_Yc_r;
   std::vector<float>   *tree_Yc_dr0;
   std::vector<float>   *tree_Yc_dr1;
   std::vector<float>   *tree_Yc_dz0;
   std::vector<float>   *tree_Yc_dz1;
   std::vector<float>   *tree_Yc_costheta;
   std::vector<int>     *tree_Yc_layer;
   std::vector<float>   *tree_Yc_NChi2;
   std::vector<float>   *tree_Yc_ndf;
   std::vector<unsigned int> *tree_Yc_nDaughters;
   std::vector<float>   *tree_Yc_pt;
   std::vector<float>   *tree_Yc_eta;
   std::vector<float>   *tree_Yc_phi;
   std::vector<float>   *tree_Yc_mass;
   Int_t           tree_Yc_ntracks;
   std::vector<int>     *tree_Yc_tracks_index;
   std::vector<int>     *tree_Yc_tracks_charge;
   std::vector<float>   *tree_Yc_tracks_pt;
   std::vector<float>   *tree_Yc_tracks_eta;
   std::vector<float>   *tree_Yc_tracks_phi;
   std::vector<float>   *tree_Yc_tracks_phi0;
   Int_t           tree_TRACK_SIZE;
   Int_t           tree_nTracks;
   Int_t           tree_nLostTracks;
   std::vector<unsigned int> *tree_track_ipc;
   std::vector<bool>    *tree_track_lost;
   std::vector<float>   *tree_track_px;
   std::vector<float>   *tree_track_py;
   std::vector<float>   *tree_track_pz;
   std::vector<float>   *tree_track_pt;
   std::vector<float>   *tree_track_eta;
   std::vector<float>   *tree_track_phi;
   std::vector<int>     *tree_track_charge;
   std::vector<float>   *tree_track_NChi2;
   std::vector<bool>    *tree_track_isHighPurity;
   std::vector<float>   *tree_track_dxy;
   std::vector<float>   *tree_track_dxyError;
   std::vector<float>   *tree_track_drSig;
   std::vector<float>   *tree_track_dz;
   std::vector<float>   *tree_track_dzError;
   std::vector<float>   *tree_track_dzSig;
   std::vector<int>     *tree_track_nHit;
   std::vector<int>     *tree_track_nHitPixel;
   std::vector<int>     *tree_track_nHitTIB;
   std::vector<int>     *tree_track_nHitTID;
   std::vector<int>     *tree_track_nHitTOB;
   std::vector<int>     *tree_track_nHitTEC;
   std::vector<int>     *tree_track_nHitPXB;
   std::vector<int>     *tree_track_nHitPXF;
   std::vector<int>     *tree_track_isHitPixel;
   std::vector<int>     *tree_track_nLayers;
   std::vector<int>     *tree_track_nLayersPixel;
   std::vector<float>   *tree_track_x;
   std::vector<float>   *tree_track_y;
   std::vector<float>   *tree_track_z;
   std::vector<int>     *tree_track_firstHit;
   std::vector<float>   *tree_track_region;
   std::vector<float>   *tree_track_firstHit_x;
   std::vector<float>   *tree_track_firstHit_y;
   std::vector<float>   *tree_track_firstHit_z;
   std::vector<int>     *tree_track_iJet;
   std::vector<float>   *tree_track_ntrk10;
   std::vector<float>   *tree_track_ntrk20;
   std::vector<float>   *tree_track_ntrk30;
   std::vector<float>   *tree_track_ntrk40;
   std::vector<double>  *tree_track_MVAval;
   std::vector<float>   *tree_track_btag;
   std::vector<float>   *tree_track_energy;
   std::vector<int>     *tree_track_Hemi;
   std::vector<float>   *tree_track_Hemi_dR;
   std::vector<float>   *tree_track_Hemi_dRmax;
   std::vector<float>   *tree_track_Hemi_mva_NChi2;
   std::vector<bool>    *tree_track_Hemi_ping;
   std::vector<float>   *tree_track_Hemi_dFirstVtx;
   std::vector<int>     *tree_track_Hemi_LLP;
   std::vector<int>     *tree_track_Hemi_isjet;
   std::vector<int>     *tree_track_sim_LLP;
   std::vector<bool>    *tree_track_sim_isFromB;
   std::vector<bool>    *tree_track_sim_isFromC;
   std::vector<float>   *tree_track_sim_pt;
   std::vector<float>   *tree_track_sim_eta;
   std::vector<float>   *tree_track_sim_phi;
   std::vector<int>     *tree_track_sim_charge;
   std::vector<int>     *tree_track_sim_pdgId;
   std::vector<float>   *tree_track_sim_mass;
   std::vector<float>   *tree_track_sim_x;
   std::vector<float>   *tree_track_sim_y;
   std::vector<float>   *tree_track_sim_z;
   std::vector<float>   *tree_track_sim_dFirstGen;
   std::vector<float>   *tree_track_sim_LLP_r;
   std::vector<float>   *tree_track_sim_LLP_z;
   std::vector<bool>    *tree_V0_track_isFromV0;
   std::vector<bool>    *tree_V0_track_isFromSI;
   std::vector<bool>    *tree_V0_track_lost;
   std::vector<float>   *tree_V0_track_pt;
   std::vector<float>   *tree_V0_track_eta;
   std::vector<float>   *tree_V0_track_phi;
   std::vector<int>     *tree_V0_track_charge;
   std::vector<float>   *tree_V0_track_NChi2;
   std::vector<float>   *tree_V0_track_dxy;
   std::vector<float>   *tree_V0_track_drSig;
   std::vector<float>   *tree_V0_track_dz;
   std::vector<float>   *tree_V0_track_dzSig;
   std::vector<int>     *tree_V0_track_nHit;
   std::vector<int>     *tree_V0_track_nHitPixel;
   std::vector<int>     *tree_V0_track_firstHit;
   std::vector<float>   *tree_V0_track_firstHit_x;
   std::vector<float>   *tree_V0_track_firstHit_y;
   std::vector<float>   *tree_V0_track_firstHit_z;
   std::vector<int>     *tree_V0_track_iJet;
   std::vector<float>   *tree_V0_track_ntrk10;
   std::vector<float>   *tree_V0_track_ntrk20;
   std::vector<float>   *tree_V0_track_ntrk30;
   std::vector<float>   *tree_V0_track_ntrk40;
   std::vector<int>     *tree_V0_track_Hemi;
   std::vector<float>   *tree_V0_track_Hemi_dR;
   std::vector<float>   *tree_V0_track_Hemi_dRmax;
   Float_t         tree_GenPVx;
   Float_t         tree_GenPVy;
   Float_t         tree_GenPVz;
   Int_t           tree_smu_mass;
   Int_t           tree_neu_mass;
   Float_t         tree_neu_ctau;
   std::vector<float>   *tree_genParticle_pt;
   std::vector<float>   *tree_genParticle_eta;
   std::vector<float>   *tree_genParticle_phi;
   std::vector<float>   *tree_genParticle_charge;
   std::vector<int>     *tree_genParticle_pdgId;
   std::vector<float>   *tree_genParticle_mass;
   std::vector<float>   *tree_genParticle_x;
   std::vector<float>   *tree_genParticle_y;
   std::vector<float>   *tree_genParticle_z;
   std::vector<float>   *tree_genParticle_px;
   std::vector<float>   *tree_genParticle_py;
   std::vector<float>   *tree_genParticle_pz;
   std::vector<float>   *tree_genParticle_energy;
   std::vector<bool>    *tree_genParticle_isPromptFinalState;
   std::vector<int>     *tree_genParticle_statusCode;
   std::vector<int>     *tree_genParticle_mother_pdgId;
   std::vector<int>     *tree_genParticle_LLP;
   std::vector<float>   *tree_genParticle_ct;
   std::vector<float>   *tree_genParticle_ct0;
   Int_t           tree_ngenPackPart;
   std::vector<float>   *tree_genPackPart_pt;
   std::vector<float>   *tree_genPackPart_eta;
   std::vector<float>   *tree_genPackPart_phi;
   std::vector<float>   *tree_genPackPart_charge;
   std::vector<int>     *tree_genPackPart_pdgId;
   std::vector<float>   *tree_genPackPart_mass;
   std::vector<float>   *tree_genPackPart_x;
   std::vector<float>   *tree_genPackPart_y;
   std::vector<float>   *tree_genPackPart_z;
   std::vector<int>     *tree_genPackPart_mother_pdgId;
   std::vector<bool>    *tree_genPackPart_isFromB;
   std::vector<bool>    *tree_genPackPart_isFromC;
   Int_t           tree_ngenFromLLP;
   std::vector<int>     *tree_genFromLLP_LLP;
   std::vector<float>   *tree_genFromLLP_pt;
   std::vector<float>   *tree_genFromLLP_eta;
   std::vector<float>   *tree_genFromLLP_phi;
   std::vector<float>   *tree_genFromLLP_charge;
   std::vector<int>     *tree_genFromLLP_pdgId;
   std::vector<float>   *tree_genFromLLP_mass;
   std::vector<float>   *tree_genFromLLP_x;
   std::vector<float>   *tree_genFromLLP_y;
   std::vector<float>   *tree_genFromLLP_z;
   std::vector<int>     *tree_genFromLLP_mother_pdgId;
   std::vector<bool>    *tree_genFromLLP_isFromB;
   std::vector<bool>    *tree_genFromLLP_isFromC;
   std::vector<float>   *tree_gen_top_pt;
   std::vector<float>   *tree_gen_top_rw_pt;
   std::vector<float>   *tree_genAxis_dRneuneu;
   std::vector<float>   *tree_genAxis_dPhineuneu;
   std::vector<float>   *tree_genAxis_dEtaneuneu;
   std::vector<float>   *tree_GenAxes_Mass;
   std::vector<float>   *tree_GenAxes_CombinedHemiLeptonMass;
   std::vector<float>   *tree_GenAxis_Neu_dRmin;
   std::vector<float>   *tree_GenAxis_Neu_dRmax;
   std::vector<float>   *tree_GenAxis_RecoAxis_dRmin;
   std::vector<float>   *tree_GenAxis_RecoAxis_dRmax;
   Int_t           tree_nFromC;
   std::vector<float>   *tree_genFromC_pt;
   std::vector<float>   *tree_genFromC_eta;
   std::vector<float>   *tree_genFromC_phi;
   std::vector<float>   *tree_genFromC_charge;
   std::vector<int>     *tree_genFromC_pdgId;
   std::vector<float>   *tree_genFromC_x;
   std::vector<float>   *tree_genFromC_y;
   std::vector<float>   *tree_genFromC_z;
   std::vector<int>     *tree_genFromC_mother_pdgId;
   Int_t           tree_nFromB;
   std::vector<float>   *tree_genFromB_pt;
   std::vector<float>   *tree_genFromB_eta;
   std::vector<float>   *tree_genFromB_phi;
   std::vector<float>   *tree_genFromB_charge;
   std::vector<int>     *tree_genFromB_pdgId;
   std::vector<float>   *tree_genFromB_x;
   std::vector<float>   *tree_genFromB_y;
   std::vector<float>   *tree_genFromB_z;
   std::vector<int>     *tree_genFromB_mother_pdgId;
   std::vector<float>   *tree_genFromB_dd;
   std::vector<float>   *tree_genFromB_dr;
   std::vector<float>   *tree_genFromB_dz;
   std::vector<float>   *tree_genJet_pt;
   std::vector<float>   *tree_genJet_eta;
   std::vector<float>   *tree_genJet_phi;
   std::vector<float>   *tree_genJet_mass;
   std::vector<float>   *tree_genJet_energy;
   Int_t           tree_nLLP;
   std::vector<int>     *tree_LLP;
   std::vector<float>   *tree_LLP_pt;
   std::vector<float>   *tree_LLP_eta;
   std::vector<float>   *tree_LLP_phi;
   std::vector<float>   *tree_LLP_x;
   std::vector<float>   *tree_LLP_y;
   std::vector<float>   *tree_LLP_z;
   std::vector<float>   *tree_LLP_r;
   std::vector<float>   *tree_LLP_dist;
   std::vector<int>     *tree_LLP_nTrks;
   std::vector<int>     *tree_LLP_Vtx_nTrks;
   std::vector<float>   *tree_LLP_Vtx_NChi2;
   std::vector<float>   *tree_LLP_Vtx_dx;
   std::vector<float>   *tree_LLP_Vtx_dy;
   std::vector<float>   *tree_LLP_Vtx_dz;
   std::vector<float>   *tree_LLP_Vtx_dist;
   std::vector<float>   *tree_LLP_Vtx_dd;
   std::vector<float>   *tree_LLP12_dR;
   std::vector<float>   *tree_LLP12_deta;
   std::vector<float>   *tree_LLP12_dphi;
   std::vector<float>   *tree_LLP_Mass;
   std::vector<int>     *tree_Hemi;
   std::vector<int>     *tree_Hemi_njet;
   std::vector<int>     *tree_Hemi_njet_nomu;
   std::vector<float>   *tree_Hemi_pt;
   std::vector<float>   *tree_Hemi_eta;
   std::vector<float>   *tree_Hemi_phi;
   std::vector<int>     *tree_Hemi_nTrks;
   std::vector<int>     *tree_Hemi_nTrks_sig;
   std::vector<int>     *tree_Hemi_nTrks_bad;
   std::vector<float>   *tree_Hemi_mass;
   std::vector<float>   *tree_HemiMu_mass;
   std::vector<float>   *tree_HemiMu_pt;
   std::vector<float>   *tree_HemiMu_dR;
   std::vector<float>   *tree_HemiMuOp_mass;
   std::vector<float>   *tree_HemiMuOp_pt;
   std::vector<float>   *tree_HemiMuOp_dR;
   std::vector<int>     *tree_Hemi_LooseBTag_axes;
   std::vector<int>     *tree_Hemi_MediumBTag_axes;
   std::vector<int>     *tree_Hemi_TightBTag_axes;
   std::vector<float>   *tree_Hemi_dR12;
   std::vector<int>     *tree_Hemi_LLP;
   std::vector<float>   *tree_Hemi_LLP_pt;
   std::vector<float>   *tree_Hemi_LLP_eta;
   std::vector<float>   *tree_Hemi_LLP_phi;
   std::vector<float>   *tree_Hemi_LLP_dist;
   std::vector<float>   *tree_Hemi_LLP_x;
   std::vector<float>   *tree_Hemi_LLP_y;
   std::vector<float>   *tree_Hemi_LLP_z;
   std::vector<float>   *tree_Hemi_LLP_dR;
   std::vector<int>     *tree_Hemi_LLP_mother;
   std::vector<float>   *tree_Hemi_LLP_Vtx_dx;
   std::vector<float>   *tree_Hemi_LLP_Vtx_dy;
   std::vector<float>   *tree_Hemi_LLP_Vtx_dz;
   std::vector<float>   *tree_Hemi_LLP_Vtx_dr;
   std::vector<float>   *tree_Hemi_LLP_muOK_dR;
   std::vector<float>   *tree_Hemi_LLP_muOK_pt;
   std::vector<float>   *tree_Hemi_LLP_muOK_mass;
   std::vector<float>   *tree_Hemi_LLP_muNO_dR;
   std::vector<float>   *tree_Hemi_LLP_muNO_pt;
   std::vector<float>   *tree_Hemi_LLP_muNO_mass;
   std::vector<float>   *tree_Hemi_LLP_dR12;
   std::vector<bool>    *tree_Hemi_LLP_ping;
   std::vector<int>     *tree_event_LLP_ping;
   std::vector<int>     *tree_Hemi_Vtx_step;
   std::vector<bool>    *tree_Hemi_Vtx_isTight;
   std::vector<float>   *tree_Hemi_Vtx_NChi2;
   std::vector<int>     *tree_Hemi_Vtx_nTrks;
   std::vector<int>     *tree_Hemi_Vtx_nTrks_sig;
   std::vector<int>     *tree_Hemi_Vtx_nTrks_bad;
   std::vector<float>   *tree_Hemi_Vtx_x;
   std::vector<float>   *tree_Hemi_Vtx_y;
   std::vector<float>   *tree_Hemi_Vtx_z;
   std::vector<float>   *tree_Hemi_Vtx_r;
   std::vector<float>   *tree_Hemi_Vtx_dR;
   std::vector<float>   *tree_Hemi_Vtx_xError;
   std::vector<float>   *tree_Hemi_Vtx_yError;
   std::vector<float>   *tree_Hemi_Vtx_zError;
   std::vector<float>   *tree_Hemi_Vtx_BTag;
   std::vector<float>   *tree_Hemi_Vtx_trackWeight;
   std::vector<float>   *tree_Hemi_Vtx_SumtrackWeight;
   std::vector<float>   *tree_Hemi_Vtx_track_MeanDCA_d;
   std::vector<float>   *tree_Hemi_Vtx_Mass;
   std::vector<float>   *tree_Hemi_Vtx_dist;
   std::vector<int>     *tree_Hemi_Vtx_ntrk10;
   std::vector<int>     *tree_Hemi_Vtx_ntrk20;
   std::vector<int>     *tree_event_nVtx;
   std::vector<float>   *tree_event_Vtx_Vtx_dr;
   std::vector<float>   *tree_event_Vtx_Vtx_dz;
   std::vector<float>   *tree_event_Vtx_Vtx_dd;
   std::vector<float>   *tree_event_Vtx_Vtx_reldd;
   std::vector<float>   *tree_event_Vtx_Vtx_dR;
   std::vector<int>     *tree_event_Vtx_Vtx_step;
   std::vector<float>   *tree_Hemi_SecLLP;
   std::vector<float>   *tree_Hemi_LLP_SecVtx_dx;
   std::vector<float>   *tree_Hemi_LLP_SecVtx_dy;
   std::vector<float>   *tree_Hemi_LLP_SecVtx_dz;
   std::vector<float>   *tree_Hemi_LLP_SecVtx_dr;
   std::vector<bool>    *tree_Hemi_SecLLP_ping;
   std::vector<int>     *tree_event_SecLLP_ping;
   std::vector<int>     *tree_Hemi_SecVtx;
   std::vector<int>     *tree_Hemi_SecVtx_step;
   std::vector<float>   *tree_Hemi_SecVtx_x;
   std::vector<float>   *tree_Hemi_SecVtx_y;
   std::vector<float>   *tree_Hemi_SecVtx_z;
   std::vector<float>   *tree_Hemi_SecVtx_r;
   std::vector<float>   *tree_Hemi_SecVtx_dR;
   std::vector<float>   *tree_Hemi_SecVtx_nTrks;
   std::vector<float>   *tree_Hemi_SecVtx_NChi2;
   std::vector<float>   *tree_Hemi_SecVtx_dist;
   std::vector<float>   *tree_Hemi_SecVtx_track_MeanDCA_d;
   std::vector<float>   *tree_Hemi_SecVtx_SumtrackWeight;
   std::vector<float>   *tree_Hemi_SecVtx_trackWeight;
   std::vector<float>   *tree_Hemi_SecVtx_Mass;
   std::vector<float>   *tree_event_MergedVtx_Vtx_dr;
   std::vector<float>   *tree_event_MergedVtx_Vtx_dz;
   std::vector<float>   *tree_event_MergedVtx_Vtx_dd;
   std::vector<float>   *tree_event_MergedVtx_Vtx_reldd;
   std::vector<float>   *tree_event_MergedVtx_Vtx_dR;
   std::vector<int>     *tree_event_MergedVtx_Vtx_step;
   std::vector<float>   *tree_Hemi_Vtx_BDT_nTrks;
   std::vector<float>   *tree_Hemi_Vtx_BDT_NChi2;
   std::vector<float>   *tree_Hemi_Vtx_BDT_step;
   std::vector<float>   *tree_Hemi_Vtx_BDT_STW;
   std::vector<float>   *tree_Hemi_Vtx_BDT_Mass;
   std::vector<float>   *tree_Hemi_Vtx_BDT_HMass;
   std::vector<float>   *tree_Hemi_Vtx_BDT_ntrk10;
   std::vector<float>   *tree_Hemi_Vtx_BDT_ntrk20;
   std::vector<float>   *tree_Hemi_Vtx_BDT_MeanDCA;
   std::vector<float>   *tree_Hemi_Vtx_MVAval_Loose;
   std::vector<float>   *tree_Hemi_Vtx_MVAval_Tight;

   Bool_t          HLT_IsoMu24_v;
   Bool_t          HLT_IsoMu27_v;
   Bool_t          HLT_IsoTkMu24_v;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v;
   Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v;
   Bool_t          HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Bool_t          HLT_Ele27_WPTight_Gsf_v;
   Bool_t          HLT_Ele32_WPTight_Gsf_v;
   Bool_t          HLT_Ele35_WPTight_Gsf_v;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_v;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_PFHT60_v;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;
   Bool_t          HLT_PFMET250_HBHECleaned_v;
   Bool_t          HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;
   Bool_t          HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;    
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v;    
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v;  


   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_tree_only_gen_wt;   //!
   TBranch        *b_tree_event_weight;   //!
   TBranch        *b_tree_genTop_Weight;   //!
   TBranch        *b_PUweight;   //!
   TBranch        *b_PUweight_Up;   //!
   TBranch        *b_PUweight_Down;   //!
   TBranch        *b_Prefweight;   //!
   TBranch        *b_PU_events;   //!
   TBranch        *b_tree_only_tigger_filter;   //!
TBranch        *b_tree_trigger_doublelepton;   //!
   TBranch        *b_tree_trigger_singlelepton;   //!
   TBranch        *b_tree_Filter;   //!
   TBranch        *b_tree_FilterSameSign;   //!
   TBranch        *b_tree_Good_PV;   //!
   TBranch        *b_tree_muon_GenRecoTriggerMatched;   //!
   TBranch        *b_tree_Evts_MVAval;   //!
   TBranch        *b_tree_Evts_MVAvalDY;   //!
   TBranch        *b_tree_Evts_MVAvalTT;   //!
   TBranch        *b_tree_bs_PosX;   //!
   TBranch        *b_tree_bs_PosY;   //!
   TBranch        *b_tree_bs_PosZ;   //!
   TBranch        *b_tree_nPV;   //!
   TBranch        *b_tree_PV_x;   //!
   TBranch        *b_tree_PV_y;   //!
   TBranch        *b_tree_PV_z;   //!
   TBranch        *b_tree_PV_ez;   //!
   TBranch        *b_tree_PV_NChi2;   //!
   TBranch        *b_tree_PV_ndf;   //!
   TBranch        *b_tree_PV_rho;   //!
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
   TBranch        *b_tree_muon_x;   //!
   TBranch        *b_tree_muon_y;   //!
   TBranch        *b_tree_muon_z;   //!
   TBranch        *b_tree_muon_dxy;   //!
   TBranch        *b_tree_muon_dxyError;   //!
   TBranch        *b_tree_muon_dz;   //!
   TBranch        *b_tree_muon_dzError;   //!
   TBranch        *b_tree_muon_charge;   //!
   TBranch        *b_tree_muon_isLoose;   //!
   TBranch        *b_tree_muon_isMedium;   //!
   TBranch        *b_tree_muon_isTight;   //!
   TBranch        *b_tree_muon_isGlobal;   //!
   TBranch        *b_tree_muon_isoR3;   //!
   TBranch        *b_tree_muon_trigger_dimu;   //!
   TBranch        *b_tree_muon_trigger_isomu;   //!
   TBranch        *b_tree_muon_PFIsoVeryLoose;   //!
   TBranch        *b_tree_muon_PFIsoLoose;   //!
   TBranch        *b_tree_muon_PFIsoMedium;   //!
   TBranch        *b_tree_muon_PFIsoTight;   //!
   TBranch        *b_tree_muon_MiniIsoLoose;   //!
   TBranch        *b_tree_muon_MiniIsoMedium;   //!
   TBranch        *b_tree_muon_MiniIsoTight;   //!
   TBranch        *b_tree_muon_TkIsoLoose;   //!
   TBranch        *b_tree_muon_TkIsoTight;   //!
   TBranch        *b_tree_muon_trkLayers;   //!
   TBranch        *b_tree_muon_miniIso;   //!
   TBranch        *b_tree_muon_correction;   //!
   TBranch        *b_tree_muon_gen;   //!
   TBranch        *b_tree_reco_muon_leadingpt;   //!
   TBranch        *b_tree_reco_electron_leadingpt2;   //!
   TBranch        *b_tree_reco_muon_leadingeta;   //!
   TBranch        *b_tree_reco_electron_leadingeta2;   //!
   TBranch        *b_tree_reco_muon_leadingphi;   //!
   TBranch        *b_tree_reco_electron_leadingphi2;   //!
   TBranch        *b_tree_trig_muon_leadingpt;   //!
   TBranch        *b_tree_trig_electron_leadingpt2;   //!
   TBranch        *b_tree_trig_muon_leadingeta;   //!
   TBranch        *b_tree_trig_electron_leadingeta2;   //!
   TBranch        *b_tree_trig_muon_leadingphi;   //!
   TBranch        *b_tree_trig_electron_leadingphi2;   //!

   TBranch        *b_tree_lepton_b4trigger_leadingpt;
   TBranch        *b_tree_lepton_b4trigger_leadingpt2;


   TBranch        *b_tree_reco_lepton_leadingpt;   //!
   TBranch        *b_tree_reco_lepton_leadingpt2;   //!
   TBranch        *b_tree_reco_lepton_leadingeta;   //!
   TBranch        *b_tree_reco_lepton_leadingeta2;   //!
   TBranch        *b_tree_reco_lepton_leadingphi;   //!
   TBranch        *b_tree_reco_lepton_leadingphi2;   //!
   TBranch        *b_tree_trig_lepton_leadingpt;   //!
   TBranch        *b_tree_trig_lepton_leadingpt2;   //!
   TBranch        *b_tree_trig_lepton_leadingeta;   //!
   TBranch        *b_tree_trig_lepton_leadingeta2;   //!
   TBranch        *b_tree_trig_lepton_leadingphi;   //!
   TBranch        *b_tree_trig_lepton_leadingphi2;   //!
   TBranch        *b_tree_lepton_leadingpt;   //!
   TBranch        *b_tree_lepton_leadingpt2;   //!
   TBranch        *b_tree_lepton_leadingeta;   //!
   TBranch        *b_tree_lepton_leadingeta2;   //!
   TBranch        *b_tree_lepton_leadingphi;   //!
   TBranch        *b_tree_lepton_leadingphi2;   //!
TBranch        *b_tree_lepton_leadingdxy;   //!
   TBranch        *b_tree_lepton_leadingdxy2;   //!
   TBranch        *b_tree_lepton_leadingdz;   //!
   TBranch        *b_tree_lepton_leadingdz2;   //!
   TBranch        *b_tree_lepton_lepton_dR;   //!
   TBranch        *b_tree_lepton_lepton_dPhi;   //!
   TBranch        *b_tree_lepton_lepton_dEta;   //!
   TBranch        *b_tree_ll_pt;   //!
   TBranch        *b_tree_ll_eta;   //!
   TBranch        *b_tree_ll_phi;   //!
   TBranch        *b_tree_ll_px;   //!
   TBranch        *b_tree_ll_py;   //!
   TBranch        *b_tree_ll_pz;   //!
   TBranch        *b_tree_ll_energy;   //!
   TBranch        *b_tree_ll_mass;   //!
   TBranch        *b_tree_all_nel;   //!
   TBranch        *b_tree_electron_nEle;   //!
   TBranch        *b_tree_electron_isPrompt;   //!
   TBranch        *b_tree_electron_pt;   //!
   TBranch        *b_tree_electron_eta;   //!
   TBranch        *b_tree_electron_phi;   //!
   TBranch        *b_tree_electron_x;   //!
   TBranch        *b_tree_electron_y;   //!
   TBranch        *b_tree_electron_z;   //!
   TBranch        *b_tree_electron_energy;   //!
   TBranch        *b_tree_electron_et;   //!
   TBranch        *b_tree_electron_ecal_trk_postcorr;   //!
   TBranch        *b_tree_electron_charge;   //!
   TBranch        *b_tree_electron_isoR4;   //!
   TBranch        *b_tree_electron_IsLoose;   //!
   TBranch        *b_tree_electron_IsMedium;   //!
   TBranch        *b_tree_electron_IsTight;   //!
   TBranch        *b_tree_electron_dxy;   //!
   TBranch        *b_tree_electron_dz;   //!
   TBranch        *b_tree_electron_gen;   //!
   TBranch        *b_tree_PFMet_et;   //!
   TBranch        *b_tree_PFMet_phi;   //!
   TBranch        *b_tree_PFMet_sig;   //!
   TBranch        *b_tree_PFMet_pt;   //!
   TBranch        *b_tree_njet;   //!
   TBranch        *b_tree_njetNOmu;   //!
   TBranch        *b_tree_jet_pt;   //!
   TBranch        *b_tree_jet_eta;   //!
   TBranch        *b_tree_jet_phi;   //!
   TBranch        *b_tree_jet_px;   //!
   TBranch        *b_tree_jet_py;   //!
   TBranch        *b_tree_jet_pz;   //!
   TBranch        *b_tree_jet_E;   //!
   TBranch        *b_tree_jet_tightid_LepVeto;   //!
   TBranch        *b_tree_jet_tightid;   //!
   TBranch        *b_tree_jet_TightJetIDLepVeto;   //!
   TBranch        *b_tree_jet_TightJetID;   //!
   TBranch        *b_tree_jet_HadronFlavour;   //!
   TBranch        *b_tree_jet_pileupID;   //!
   TBranch        *b_tree_jet_btag_DeepCSV;   //!
   TBranch        *b_tree_jet_btag_DeepJet;   //!
   TBranch        *b_tree_jet_leadingpt;   //!
   TBranch        *b_tree_jet_leadingpt2;   //!
   TBranch        *b_tree_jet_leadingeta;   //!
   TBranch        *b_tree_jet_leadingeta2;   //!
   TBranch        *b_tree_jet_leadingMuon_dR;   //!
   TBranch        *b_tree_jet_leadingMuon2_dR;   //!
   TBranch        *b_tree_jet_jet_dR;   //!
   TBranch        *b_tree_jet_jet_dPhi;   //!
   TBranch        *b_tree_jet_jet_dEta;   //!
   TBranch        *b_tree_muon_jet_dRmin;   //!
   TBranch        *b_tree_muon_jet_dRmax;   //!
   TBranch        *b_tree_elemu_jet_dRmin;   //!
   TBranch        *b_tree_elemu_jet_dRmax;   //!
   TBranch        *b_tree_ele_jet_dRmin;   //!
   TBranch        *b_tree_ele_jet_dRmax;   //!
   TBranch        *b_tree_HT;   //!
   TBranch        *b_tree_nK0;   //!
   TBranch        *b_tree_K0_x;   //!
   TBranch        *b_tree_K0_y;   //!
   TBranch        *b_tree_K0_z;   //!
   TBranch        *b_tree_K0_r;   //!
   TBranch        *b_tree_K0_NChi2;   //!
   TBranch        *b_tree_K0_ndf;   //!
   TBranch        *b_tree_K0_mass;   //!
   TBranch        *b_tree_K0_pt;   //!
   TBranch        *b_tree_K0_eta;   //!
   TBranch        *b_tree_K0_phi;   //!
   TBranch        *b_tree_K0_nDaughters;   //!
   TBranch        *b_tree_nLambda;   //!
   TBranch        *b_tree_L0_x;   //!
   TBranch        *b_tree_L0_y;   //!
   TBranch        *b_tree_L0_z;   //!
   TBranch        *b_tree_L0_r;   //!
   TBranch        *b_tree_L0_NChi2;   //!
   TBranch        *b_tree_L0_ndf;   //!
   TBranch        *b_tree_L0_nDaughters;   //!
   TBranch        *b_tree_L0_mass;   //!
   TBranch        *b_tree_L0_pt;   //!
   TBranch        *b_tree_L0_eta;   //!
   TBranch        *b_tree_L0_phi;   //!
   TBranch        *b_tree_nV0_reco;   //!
   TBranch        *b_tree_V0_reco_x;   //!
   TBranch        *b_tree_V0_reco_y;   //!
   TBranch        *b_tree_V0_reco_z;   //!
   TBranch        *b_tree_V0_reco_r;   //!
   TBranch        *b_tree_V0_reco_drSig;   //!
   TBranch        *b_tree_V0_reco_dzSig;   //!
   TBranch        *b_tree_V0_reco_angleXY;   //!
   TBranch        *b_tree_V0_reco_angleZ;   //!
   TBranch        *b_tree_V0_reco_NChi2;   //!
   TBranch        *b_tree_V0_reco_ndf;   //!
   TBranch        *b_tree_V0_reco_mass;   //!
   TBranch        *b_tree_V0_reco_pt;   //!
   TBranch        *b_tree_V0_reco_eta;   //!
   TBranch        *b_tree_V0_reco_phi;   //!
   TBranch        *b_tree_V0_reco_source;   //!
   TBranch        *b_tree_V0_reco_badTkHit;   //!
   TBranch        *b_tree_V0_reco_dca;   //!
   TBranch        *b_tree_nSecInt;   //!
   TBranch        *b_tree_SecInt_x;   //!
   TBranch        *b_tree_SecInt_y;   //!
   TBranch        *b_tree_SecInt_z;   //!
   TBranch        *b_tree_SecInt_r;   //!
   TBranch        *b_tree_SecInt_d;   //!
   TBranch        *b_tree_SecInt_drSig;   //!
   TBranch        *b_tree_SecInt_dzSig;   //!
   TBranch        *b_tree_SecInt_angleXY;   //!
   TBranch        *b_tree_SecInt_angleZ;   //!
   TBranch        *b_tree_SecInt_NChi2;   //!
   TBranch        *b_tree_SecInt_ndf;   //!
   TBranch        *b_tree_SecInt_mass;   //!
   TBranch        *b_tree_SecInt_pt;   //!
   TBranch        *b_tree_SecInt_eta;   //!
   TBranch        *b_tree_SecInt_phi;   //!
   TBranch        *b_tree_SecInt_charge;   //!
   TBranch        *b_tree_SecInt_badTkHit;   //!
   TBranch        *b_tree_SecInt_dca;   //!
   TBranch        *b_tree_SecInt_selec;   //!
   TBranch        *b_tree_SecInt_layer;   //!
   TBranch        *b_tree_SecInt_LLP;   //!
   TBranch        *b_tree_SecInt_LLP_dr;   //!
   TBranch        *b_tree_SecInt_LLP_dz;   //!
   TBranch        *b_tree_SecInt_LLP_dd;   //!
   TBranch        *b_tree_SecInt_tk1;   //!
   TBranch        *b_tree_SecInt_tk2;   //!
   TBranch        *b_tree_nYConv;   //!
   TBranch        *b_tree_Yc_x;   //!
   TBranch        *b_tree_Yc_y;   //!
   TBranch        *b_tree_Yc_z;   //!
   TBranch        *b_tree_Yc_r;   //!
   TBranch        *b_tree_Yc_dr0;   //!
   TBranch        *b_tree_Yc_dr1;   //!
   TBranch        *b_tree_Yc_dz0;   //!
   TBranch        *b_tree_Yc_dz1;   //!
   TBranch        *b_tree_Yc_costheta;   //!
   TBranch        *b_tree_Yc_layer;   //!
   TBranch        *b_tree_Yc_NChi2;   //!
   TBranch        *b_tree_Yc_ndf;   //!
   TBranch        *b_tree_Yc_nDaughters;   //!
   TBranch        *b_tree_Yc_pt;   //!
   TBranch        *b_tree_Yc_eta;   //!
   TBranch        *b_tree_Yc_phi;   //!
   TBranch        *b_tree_Yc_mass;   //!
   TBranch        *b_tree_Yc_ntracks;   //!
   TBranch        *b_tree_Yc_tracks_index;   //!
   TBranch        *b_tree_Yc_tracks_charge;   //!
   TBranch        *b_tree_Yc_tracks_pt;   //!
   TBranch        *b_tree_Yc_tracks_eta;   //!
   TBranch        *b_tree_Yc_tracks_phi;   //!
   TBranch        *b_tree_Yc_tracks_phi0;   //!
   TBranch        *b_tree_TRACK_SIZE;   //!
   TBranch        *b_tree_nTracks;   //!
   TBranch        *b_tree_nLostTracks;   //!
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
   TBranch        *b_tree_track_btag;   //!
   TBranch        *b_tree_track_energy;   //!
   TBranch        *b_tree_track_Hemi;   //!
   TBranch        *b_tree_track_Hemi_dR;   //!
   TBranch        *b_tree_track_Hemi_dRmax;   //!
   TBranch        *b_tree_track_Hemi_mva_NChi2;   //!
   TBranch        *b_tree_track_Hemi_ping;   //!
   TBranch        *b_tree_track_Hemi_dFirstVtx;   //!
   TBranch        *b_tree_track_Hemi_LLP;   //!
   TBranch        *b_tree_track_Hemi_isjet;   //!
   TBranch        *b_tree_track_sim_LLP;   //!
   TBranch        *b_tree_track_sim_isFromB;   //!
   TBranch        *b_tree_track_sim_isFromC;   //!
   TBranch        *b_tree_track_sim_pt;   //!
   TBranch        *b_tree_track_sim_eta;   //!
   TBranch        *b_tree_track_sim_phi;   //!
   TBranch        *b_tree_track_sim_charge;   //!
   TBranch        *b_tree_track_sim_pdgId;   //!
   TBranch        *b_tree_track_sim_mass;   //!
   TBranch        *b_tree_track_sim_x;   //!
   TBranch        *b_tree_track_sim_y;   //!
   TBranch        *b_tree_track_sim_z;   //!
   TBranch        *b_tree_track_sim_dFirstGen;   //!
   TBranch        *b_tree_track_sim_LLP_r;   //!
   TBranch        *b_tree_track_sim_LLP_z;   //!
   TBranch        *b_tree_V0_track_isFromV0;   //!
   TBranch        *b_tree_V0_track_isFromSI;   //!
   TBranch        *b_tree_V0_track_lost;   //!
   TBranch        *b_tree_V0_track_pt;   //!
   TBranch        *b_tree_V0_track_eta;   //!
   TBranch        *b_tree_V0_track_phi;   //!
   TBranch        *b_tree_V0_track_charge;   //!
   TBranch        *b_tree_V0_track_NChi2;   //!
   TBranch        *b_tree_V0_track_dxy;   //!
   TBranch        *b_tree_V0_track_drSig;   //!
   TBranch        *b_tree_V0_track_dz;   //!
   TBranch        *b_tree_V0_track_dzSig;   //!
   TBranch        *b_tree_V0_track_nHit;   //!
   TBranch        *b_tree_V0_track_nHitPixel;   //!
   TBranch        *b_tree_V0_track_firstHit;   //!
   TBranch        *b_tree_V0_track_firstHit_x;   //!
   TBranch        *b_tree_V0_track_firstHit_y;   //!
   TBranch        *b_tree_V0_track_firstHit_z;   //!
   TBranch        *b_tree_V0_track_iJet;   //!
   TBranch        *b_tree_V0_track_ntrk10;   //!
   TBranch        *b_tree_V0_track_ntrk20;   //!
   TBranch        *b_tree_V0_track_ntrk30;   //!
   TBranch        *b_tree_V0_track_ntrk40;   //!
   TBranch        *b_tree_V0_track_Hemi;   //!
   TBranch        *b_tree_V0_track_Hemi_dR;   //!
   TBranch        *b_tree_V0_track_Hemi_dRmax;   //!
   TBranch        *b_tree_GenPVx;   //!
   TBranch        *b_tree_GenPVy;   //!
   TBranch        *b_tree_GenPVz;   //!
   TBranch        *b_tree_smu_mass;   //!
   TBranch        *b_tree_neu_mass;   //!
   TBranch        *b_tree_neu_ctau;   //!
   TBranch        *b_tree_genParticle_pt;   //!
   TBranch        *b_tree_genParticle_eta;   //!
   TBranch        *b_tree_genParticle_phi;   //!
   TBranch        *b_tree_genParticle_charge;   //!
   TBranch        *b_tree_genParticle_pdgId;   //!
   TBranch        *b_tree_genParticle_mass;   //!
   TBranch        *b_tree_genParticle_x;   //!
   TBranch        *b_tree_genParticle_y;   //!
   TBranch        *b_tree_genParticle_z;   //!
   TBranch        *b_tree_genParticle_px;   //!
   TBranch        *b_tree_genParticle_py;   //!
   TBranch        *b_tree_genParticle_pz;   //!
   TBranch        *b_tree_genParticle_energy;   //!
   TBranch        *b_tree_genParticle_isPromptFinalState;   //!
   TBranch        *b_tree_genParticle_statusCode;   //!
   TBranch        *b_tree_genParticle_mother_pdgId;   //!
   TBranch        *b_tree_genParticle_LLP;   //!
   TBranch        *b_tree_genParticle_ct;   //!
   TBranch        *b_tree_genParticle_ct0;   //!
   TBranch        *b_tree_ngenPackPart;   //!
   TBranch        *b_tree_genPackPart_pt;   //!
   TBranch        *b_tree_genPackPart_eta;   //!
   TBranch        *b_tree_genPackPart_phi;   //!
   TBranch        *b_tree_genPackPart_charge;   //!
   TBranch        *b_tree_genPackPart_pdgId;   //!
   TBranch        *b_tree_genPackPart_mass;   //!
   TBranch        *b_tree_genPackPart_x;   //!
   TBranch        *b_tree_genPackPart_y;   //!
   TBranch        *b_tree_genPackPart_z;   //!
   TBranch        *b_tree_genPackPart_mother_pdgId;   //!
   TBranch        *b_tree_genPackPart_isFromB;   //!
   TBranch        *b_tree_genPackPart_isFromC;   //!
   TBranch        *b_tree_ngenFromLLP;   //!
   TBranch        *b_tree_genFromLLP_LLP;   //!
   TBranch        *b_tree_genFromLLP_pt;   //!
   TBranch        *b_tree_genFromLLP_eta;   //!
   TBranch        *b_tree_genFromLLP_phi;   //!
   TBranch        *b_tree_genFromLLP_charge;   //!
   TBranch        *b_tree_genFromLLP_pdgId;   //!
   TBranch        *b_tree_genFromLLP_mass;   //!
   TBranch        *b_tree_genFromLLP_x;   //!
   TBranch        *b_tree_genFromLLP_y;   //!
   TBranch        *b_tree_genFromLLP_z;   //!
   TBranch        *b_tree_genFromLLP_mother_pdgId;   //!
   TBranch        *b_tree_genFromLLP_isFromB;   //!
   TBranch        *b_tree_genFromLLP_isFromC;   //!
TBranch        *b_tree_gen_top_pt;   //!
   TBranch        *b_tree_gen_top_rw_pt;   //!
   TBranch        *b_tree_genAxis_dRneuneu;   //!
   TBranch        *b_tree_genAxis_dPhineuneu;   //!
   TBranch        *b_tree_genAxis_dEtaneuneu;   //!
   TBranch        *b_tree_GenAxes_Mass;   //!
   TBranch        *b_tree_GenAxes_CombinedHemiLeptonMass;   //!
   TBranch        *b_tree_GenAxis_Neu_dRmin;   //!
   TBranch        *b_tree_GenAxis_Neu_dRmax;   //!
   TBranch        *b_tree_GenAxis_RecoAxis_dRmin;   //!
   TBranch        *b_tree_GenAxis_RecoAxis_dRmax;   //!
   TBranch        *b_tree_nFromC;   //!
   TBranch        *b_tree_genFromC_pt;   //!
   TBranch        *b_tree_genFromC_eta;   //!
   TBranch        *b_tree_genFromC_phi;   //!
   TBranch        *b_tree_genFromC_charge;   //!
   TBranch        *b_tree_genFromC_pdgId;   //!
   TBranch        *b_tree_genFromC_x;   //!
   TBranch        *b_tree_genFromC_y;   //!
   TBranch        *b_tree_genFromC_z;   //!
   TBranch        *b_tree_genFromC_mother_pdgId;   //!
   TBranch        *b_tree_nFromB;   //!
   TBranch        *b_tree_genFromB_pt;   //!
   TBranch        *b_tree_genFromB_eta;   //!
   TBranch        *b_tree_genFromB_phi;   //!
   TBranch        *b_tree_genFromB_charge;   //!
   TBranch        *b_tree_genFromB_pdgId;   //!
   TBranch        *b_tree_genFromB_x;   //!
   TBranch        *b_tree_genFromB_y;   //!
   TBranch        *b_tree_genFromB_z;   //!
   TBranch        *b_tree_genFromB_mother_pdgId;   //!
   TBranch        *b_tree_genFromB_dd;   //!
   TBranch        *b_tree_genFromB_dr;   //!
   TBranch        *b_tree_genFromB_dz;   //!
   TBranch        *b_tree_genJet_pt;   //!
   TBranch        *b_tree_genJet_eta;   //!
   TBranch        *b_tree_genJet_phi;   //!
   TBranch        *b_tree_genJet_mass;   //!
   TBranch        *b_tree_genJet_energy;   //!
   TBranch        *b_tree_nLLP;   //!
   TBranch        *b_tree_LLP;   //!
   TBranch        *b_tree_LLP_pt;   //!
   TBranch        *b_tree_LLP_eta;   //!
   TBranch        *b_tree_LLP_phi;   //!
   TBranch        *b_tree_LLP_x;   //!
   TBranch        *b_tree_LLP_y;   //!
   TBranch        *b_tree_LLP_z;   //!
   TBranch        *b_tree_LLP_r;   //!
   TBranch        *b_tree_LLP_dist;   //!
   TBranch        *b_tree_LLP_nTrks;   //!
   TBranch        *b_tree_LLP_Vtx_nTrks;   //!
   TBranch        *b_tree_LLP_Vtx_NChi2;   //!
   TBranch        *b_tree_LLP_Vtx_dx;   //!
   TBranch        *b_tree_LLP_Vtx_dy;   //!
   TBranch        *b_tree_LLP_Vtx_dz;   //!
   TBranch        *b_tree_LLP_Vtx_dist;   //!
   TBranch        *b_tree_LLP_Vtx_dd;   //!
   TBranch        *b_tree_LLP12_dR;   //!
   TBranch        *b_tree_LLP12_deta;   //!
   TBranch        *b_tree_LLP12_dphi;   //!
   TBranch        *b_tree_LLP_Mass;   //!
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
   TBranch        *b_tree_Hemi_LooseBTag_axes;   //!
   TBranch        *b_tree_Hemi_MediumBTag_axes;   //!
   TBranch        *b_tree_Hemi_TightBTag_axes;   //!
   TBranch        *b_tree_Hemi_dR12;   //!
   TBranch        *b_tree_Hemi_LLP;   //!
   TBranch        *b_tree_Hemi_LLP_pt;   //!
   TBranch        *b_tree_Hemi_LLP_eta;   //!
   TBranch        *b_tree_Hemi_LLP_phi;   //!
   TBranch        *b_tree_Hemi_LLP_dist;   //!
   TBranch        *b_tree_Hemi_LLP_x;   //!
   TBranch        *b_tree_Hemi_LLP_y;   //!
   TBranch        *b_tree_Hemi_LLP_z;   //!
   TBranch        *b_tree_Hemi_LLP_dR;   //!
   TBranch        *b_tree_Hemi_LLP_mother;   //!
   TBranch        *b_tree_Hemi_LLP_Vtx_dx;   //!
   TBranch        *b_tree_Hemi_LLP_Vtx_dy;   //!
   TBranch        *b_tree_Hemi_LLP_Vtx_dz;   //!
   TBranch        *b_tree_Hemi_LLP_Vtx_dr;   //!
   TBranch        *b_tree_Hemi_LLP_muOK_dR;   //!
   TBranch        *b_tree_Hemi_LLP_muOK_pt;   //!
   TBranch        *b_tree_Hemi_LLP_muOK_mass;   //!
   TBranch        *b_tree_Hemi_LLP_muNO_dR;   //!
   TBranch        *b_tree_Hemi_LLP_muNO_pt;   //!
   TBranch        *b_tree_Hemi_LLP_muNO_mass;   //!
   TBranch        *b_tree_Hemi_LLP_dR12;   //!
   TBranch        *b_tree_Hemi_LLP_ping;   //!
   TBranch        *b_tree_event_LLP_ping;   //!
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
   TBranch        *b_tree_Hemi_Vtx_xError;   //!
   TBranch        *b_tree_Hemi_Vtx_yError;   //!
   TBranch        *b_tree_Hemi_Vtx_zError;   //!
   TBranch        *b_tree_Hemi_Vtx_BTag;   //!
   TBranch        *b_tree_Hemi_Vtx_trackWeight;   //!
   TBranch        *b_tree_Hemi_Vtx_SumtrackWeight;   //!
   TBranch        *b_tree_Hemi_Vtx_track_MeanDCA_d;   //!
   TBranch        *b_tree_Hemi_Vtx_Mass;   //!
   TBranch        *b_tree_Hemi_Vtx_dist;   //!
   TBranch        *b_tree_Hemi_Vtx_ntrk10;   //!
   TBranch        *b_tree_Hemi_Vtx_ntrk20;   //!
   TBranch        *b_tree_event_nVtx;   //!
   TBranch        *b_tree_event_Vtx_Vtx_dr;   //!
   TBranch        *b_tree_event_Vtx_Vtx_dz;   //!
   TBranch        *b_tree_event_Vtx_Vtx_dd;   //!
   TBranch        *b_tree_event_Vtx_Vtx_reldd;   //!
   TBranch        *b_tree_event_Vtx_Vtx_dR;   //!
   TBranch        *b_tree_event_Vtx_Vtx_step;   //!
   TBranch        *b_tree_Hemi_SecLLP;   //!
   TBranch        *b_tree_Hemi_LLP_SecVtx_dx;   //!
   TBranch        *b_tree_Hemi_LLP_SecVtx_dy;   //!
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
   TBranch        *b_tree_Hemi_SecVtx_trackWeight;   //!
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
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v;   //!
TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;   //!
TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_Ele27_WPTight_Gsf_v;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf_v;   //!
   TBranch        *b_HLT_Ele35_WPTight_Gsf_v;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_v;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;   //!
   TBranch        *b_HLT_PFMET250_HBHECleaned_v;   //!
   TBranch        *b_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;    
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v;    
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v;  
   TBranch        *b_HLT_IsoMu24_v;   //!
   TBranch        *b_HLT_IsoMu27_v;   //!
   TBranch        *b_HLT_IsoTkMu24_v;   //!

   MiniDATAMCNtuple(TTree *tree=0);
   virtual ~MiniDATAMCNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString sample ="", TString Production="",bool Signal = true);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MiniDATAMCNtuple_cxx
MiniDATAMCNtuple::MiniDATAMCNtuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Ntuple_cms.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Ntuple_cms.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Ntuple_cms.root:/FlyingTop");
      dir->GetObject("ttree",tree);

   }
   Init(tree);
}

MiniDATAMCNtuple::~MiniDATAMCNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MiniDATAMCNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MiniDATAMCNtuple::LoadTree(Long64_t entry)
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

void MiniDATAMCNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   tree_muon_isPrompt = 0;
   tree_muon_pt = 0;
   tree_muon_SF = 0;
   tree_muon_eta = 0;
   tree_muon_phi = 0;
   tree_muon_x = 0;
   tree_muon_y = 0;
   tree_muon_z = 0;
   tree_muon_dxy = 0;
   tree_muon_dxyError = 0;
   tree_muon_dz = 0;
   tree_muon_dzError = 0;
   tree_muon_charge = 0;
   tree_muon_isLoose = 0;
   tree_muon_isMedium = 0;
   tree_muon_isTight = 0;
   tree_muon_isGlobal = 0;
   tree_muon_isoR3 = 0;
   tree_muon_trigger_dimu = 0;
   tree_muon_trigger_isomu = 0;
   tree_muon_PFIsoVeryLoose = 0;
   tree_muon_PFIsoLoose = 0;
   tree_muon_PFIsoMedium = 0;
   tree_muon_PFIsoTight = 0;
   tree_muon_MiniIsoLoose = 0;
   tree_muon_MiniIsoMedium = 0;
   tree_muon_MiniIsoTight = 0;
   tree_muon_TkIsoLoose = 0;
   tree_muon_TkIsoTight = 0;
   tree_muon_trkLayers = 0;
   tree_muon_miniIso = 0;
   tree_muon_correction = 0;
   tree_muon_gen = 0;
   tree_reco_muon_leadingpt = 0;
   tree_reco_electron_leadingpt2 = 0;
   tree_reco_muon_leadingeta = 0;
   tree_reco_electron_leadingeta2 = 0;
   tree_reco_muon_leadingphi = 0;
   tree_reco_electron_leadingphi2 = 0;
   tree_trig_muon_leadingpt = 0;
   tree_trig_electron_leadingpt2 = 0;
   tree_trig_muon_leadingeta = 0;
   tree_trig_electron_leadingeta2 = 0;
   tree_trig_muon_leadingphi = 0;
   tree_trig_electron_leadingphi2 = 0;
   tree_lepton_b4trigger_leadingpt = 0;
   tree_lepton_b4trigger_leadingpt2 = 0;
   tree_reco_lepton_leadingpt = 0;
   tree_reco_lepton_leadingpt2 = 0;
   tree_reco_lepton_leadingeta = 0;
   tree_reco_lepton_leadingeta2 = 0;
   tree_reco_lepton_leadingphi = 0;
   tree_reco_lepton_leadingphi2 = 0;
   tree_trig_lepton_leadingpt = 0;
   tree_trig_lepton_leadingpt2 = 0;
   tree_trig_lepton_leadingeta = 0;
   tree_trig_lepton_leadingeta2 = 0;
   tree_trig_lepton_leadingphi = 0;
   tree_trig_lepton_leadingphi2 = 0;
   tree_lepton_leadingpt = 0;
   tree_lepton_leadingpt2 = 0;
   tree_lepton_leadingeta = 0;
   tree_lepton_leadingeta2 = 0;
   tree_lepton_leadingphi = 0;
   tree_lepton_leadingphi2 = 0;
tree_lepton_leadingdxy = 0;
   tree_lepton_leadingdxy2 = 0;
   tree_lepton_leadingdz = 0;
   tree_lepton_leadingdz2 = 0;
   tree_lepton_lepton_dR = 0;
   tree_lepton_lepton_dPhi = 0;
   tree_lepton_lepton_dEta = 0;
   tree_ll_pt = 0;
   tree_ll_eta = 0;
   tree_ll_phi = 0;
   tree_ll_px = 0;
   tree_ll_py = 0;
   tree_ll_pz = 0;
   tree_ll_energy = 0;
   tree_ll_mass = 0;
   tree_electron_isPrompt = 0;
   tree_electron_pt = 0;
   tree_electron_eta = 0;
   tree_electron_phi = 0;
   tree_electron_x = 0;
   tree_electron_y = 0;
   tree_electron_z = 0;
   tree_electron_energy = 0;
   tree_electron_et = 0;
   tree_electron_ecal_trk_postcorr = 0;
   tree_electron_charge = 0;
   tree_electron_isoR4 = 0;
   tree_electron_IsLoose = 0;
   tree_electron_IsMedium = 0;
   tree_electron_IsTight = 0;
   tree_electron_dxy = 0;
   tree_electron_dz = 0;
   tree_electron_gen = 0;
   tree_jet_pt = 0;
   tree_jet_eta = 0;
   tree_jet_phi = 0;
   tree_jet_px = 0;
   tree_jet_py = 0;
   tree_jet_pz = 0;
   tree_jet_E = 0;
   tree_jet_tightid_LepVeto = 0;
   tree_jet_tightid = 0;
   tree_jet_TightJetIDLepVeto = 0;
   tree_jet_TightJetID = 0;
   tree_jet_HadronFlavour = 0;
   tree_jet_pileupID = 0;
   tree_jet_btag_DeepCSV = 0;
   tree_jet_btag_DeepJet = 0;
   tree_jet_leadingpt = 0;
   tree_jet_leadingpt2 = 0;
   tree_jet_leadingeta = 0;
   tree_jet_leadingeta2 = 0;
   tree_jet_leadingMuon_dR = 0;
   tree_jet_leadingMuon2_dR = 0;
   tree_jet_jet_dR = 0;
   tree_jet_jet_dPhi = 0;
   tree_jet_jet_dEta = 0;
   tree_muon_jet_dRmin = 0;
   tree_muon_jet_dRmax = 0;
   tree_elemu_jet_dRmin = 0;
   tree_elemu_jet_dRmax = 0;
   tree_ele_jet_dRmin = 0;
   tree_ele_jet_dRmax = 0;
   tree_K0_x = 0;
   tree_K0_y = 0;
   tree_K0_z = 0;
   tree_K0_r = 0;
   tree_K0_NChi2 = 0;
   tree_K0_ndf = 0;
   tree_K0_mass = 0;
   tree_K0_pt = 0;
   tree_K0_eta = 0;
   tree_K0_phi = 0;
   tree_K0_nDaughters = 0;
   tree_L0_x = 0;
   tree_L0_y = 0;
   tree_L0_z = 0;
   tree_L0_r = 0;
   tree_L0_NChi2 = 0;
   tree_L0_ndf = 0;
   tree_L0_nDaughters = 0;
   tree_L0_mass = 0;
   tree_L0_pt = 0;
   tree_L0_eta = 0;
   tree_L0_phi = 0;
   tree_V0_reco_x = 0;
   tree_V0_reco_y = 0;
   tree_V0_reco_z = 0;
   tree_V0_reco_r = 0;
   tree_V0_reco_drSig = 0;
   tree_V0_reco_dzSig = 0;
   tree_V0_reco_angleXY = 0;
   tree_V0_reco_angleZ = 0;
   tree_V0_reco_NChi2 = 0;
   tree_V0_reco_ndf = 0;
   tree_V0_reco_mass = 0;
   tree_V0_reco_pt = 0;
   tree_V0_reco_eta = 0;
   tree_V0_reco_phi = 0;
   tree_V0_reco_source = 0;
   tree_V0_reco_badTkHit = 0;
   tree_V0_reco_dca = 0;
   tree_SecInt_x = 0;
   tree_SecInt_y = 0;
   tree_SecInt_z = 0;
   tree_SecInt_r = 0;
   tree_SecInt_d = 0;
   tree_SecInt_drSig = 0;
   tree_SecInt_dzSig = 0;
   tree_SecInt_angleXY = 0;
   tree_SecInt_angleZ = 0;
   tree_SecInt_NChi2 = 0;
   tree_SecInt_ndf = 0;
   tree_SecInt_mass = 0;
   tree_SecInt_pt = 0;
   tree_SecInt_eta = 0;
   tree_SecInt_phi = 0;
   tree_SecInt_charge = 0;
   tree_SecInt_badTkHit = 0;
   tree_SecInt_dca = 0;
   tree_SecInt_selec = 0;
   tree_SecInt_layer = 0;
   tree_SecInt_LLP = 0;
   tree_SecInt_LLP_dr = 0;
   tree_SecInt_LLP_dz = 0;
   tree_SecInt_LLP_dd = 0;
   tree_SecInt_tk1 = 0;
   tree_SecInt_tk2 = 0;
   tree_Yc_x = 0;
   tree_Yc_y = 0;
   tree_Yc_z = 0;
   tree_Yc_r = 0;
   tree_Yc_dr0 = 0;
   tree_Yc_dr1 = 0;
   tree_Yc_dz0 = 0;
   tree_Yc_dz1 = 0;
   tree_Yc_costheta = 0;
   tree_Yc_layer = 0;
   tree_Yc_NChi2 = 0;
   tree_Yc_ndf = 0;
   tree_Yc_nDaughters = 0;
   tree_Yc_pt = 0;
   tree_Yc_eta = 0;
   tree_Yc_phi = 0;
   tree_Yc_mass = 0;
   tree_Yc_tracks_index = 0;
   tree_Yc_tracks_charge = 0;
   tree_Yc_tracks_pt = 0;
   tree_Yc_tracks_eta = 0;
   tree_Yc_tracks_phi = 0;
   tree_Yc_tracks_phi0 = 0;
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
   tree_track_btag = 0;
   tree_track_energy = 0;
   tree_track_Hemi = 0;
   tree_track_Hemi_dR = 0;
   tree_track_Hemi_dRmax = 0;
   tree_track_Hemi_mva_NChi2 = 0;
   tree_track_Hemi_ping = 0;
   tree_track_Hemi_dFirstVtx = 0;
   tree_track_Hemi_LLP = 0;
   tree_track_Hemi_isjet = 0;
   tree_track_sim_LLP = 0;
   tree_track_sim_isFromB = 0;
   tree_track_sim_isFromC = 0;
   tree_track_sim_pt = 0;
   tree_track_sim_eta = 0;
   tree_track_sim_phi = 0;
   tree_track_sim_charge = 0;
   tree_track_sim_pdgId = 0;
   tree_track_sim_mass = 0;
   tree_track_sim_x = 0;
   tree_track_sim_y = 0;
   tree_track_sim_z = 0;
   tree_track_sim_dFirstGen = 0;
   tree_track_sim_LLP_r = 0;
   tree_track_sim_LLP_z = 0;
   tree_V0_track_isFromV0 = 0;
   tree_V0_track_isFromSI = 0;
   tree_V0_track_lost = 0;
   tree_V0_track_pt = 0;
   tree_V0_track_eta = 0;
   tree_V0_track_phi = 0;
   tree_V0_track_charge = 0;
   tree_V0_track_NChi2 = 0;
   tree_V0_track_dxy = 0;
   tree_V0_track_drSig = 0;
   tree_V0_track_dz = 0;
   tree_V0_track_dzSig = 0;
   tree_V0_track_nHit = 0;
   tree_V0_track_nHitPixel = 0;
   tree_V0_track_firstHit = 0;
   tree_V0_track_firstHit_x = 0;
   tree_V0_track_firstHit_y = 0;
   tree_V0_track_firstHit_z = 0;
   tree_V0_track_iJet = 0;
   tree_V0_track_ntrk10 = 0;
   tree_V0_track_ntrk20 = 0;
   tree_V0_track_ntrk30 = 0;
   tree_V0_track_ntrk40 = 0;
   tree_V0_track_Hemi = 0;
   tree_V0_track_Hemi_dR = 0;
   tree_V0_track_Hemi_dRmax = 0;
   tree_genParticle_pt = 0;
   tree_genParticle_eta = 0;
   tree_genParticle_phi = 0;
   tree_genParticle_charge = 0;
   tree_genParticle_pdgId = 0;
   tree_genParticle_mass = 0;
   tree_genParticle_x = 0;
   tree_genParticle_y = 0;
   tree_genParticle_z = 0;
   tree_genParticle_px = 0;
   tree_genParticle_py = 0;
   tree_genParticle_pz = 0;
   tree_genParticle_energy = 0;
   tree_genParticle_isPromptFinalState = 0;
   tree_genParticle_statusCode = 0;
   tree_genParticle_mother_pdgId = 0;
   tree_genParticle_LLP = 0;
   tree_genParticle_ct = 0;
   tree_genParticle_ct0 = 0;
   tree_genPackPart_pt = 0;
   tree_genPackPart_eta = 0;
   tree_genPackPart_phi = 0;
   tree_genPackPart_charge = 0;
   tree_genPackPart_pdgId = 0;
   tree_genPackPart_mass = 0;
   tree_genPackPart_x = 0;
   tree_genPackPart_y = 0;
   tree_genPackPart_z = 0;
   tree_genPackPart_mother_pdgId = 0;
   tree_genPackPart_isFromB = 0;
   tree_genPackPart_isFromC = 0;
   tree_genFromLLP_LLP = 0;
   tree_genFromLLP_pt = 0;
   tree_genFromLLP_eta = 0;
   tree_genFromLLP_phi = 0;
   tree_genFromLLP_charge = 0;
   tree_genFromLLP_pdgId = 0;
   tree_genFromLLP_mass = 0;
   tree_genFromLLP_x = 0;
   tree_genFromLLP_y = 0;
   tree_genFromLLP_z = 0;
   tree_genFromLLP_mother_pdgId = 0;
   tree_genFromLLP_isFromB = 0;
   tree_genFromLLP_isFromC = 0;
tree_gen_top_pt = 0;
   tree_gen_top_rw_pt = 0;
   tree_genAxis_dRneuneu = 0;
   tree_genAxis_dPhineuneu = 0;
   tree_genAxis_dEtaneuneu = 0;
   tree_GenAxes_Mass = 0;
   tree_GenAxes_CombinedHemiLeptonMass = 0;
   tree_GenAxis_Neu_dRmin = 0;
   tree_GenAxis_Neu_dRmax = 0;
   tree_GenAxis_RecoAxis_dRmin = 0;
   tree_GenAxis_RecoAxis_dRmax = 0;
   tree_genFromC_pt = 0;
   tree_genFromC_eta = 0;
   tree_genFromC_phi = 0;
   tree_genFromC_charge = 0;
   tree_genFromC_pdgId = 0;
   tree_genFromC_x = 0;
   tree_genFromC_y = 0;
   tree_genFromC_z = 0;
   tree_genFromC_mother_pdgId = 0;
   tree_genFromB_pt = 0;
   tree_genFromB_eta = 0;
   tree_genFromB_phi = 0;
   tree_genFromB_charge = 0;
   tree_genFromB_pdgId = 0;
   tree_genFromB_x = 0;
   tree_genFromB_y = 0;
   tree_genFromB_z = 0;
   tree_genFromB_mother_pdgId = 0;
   tree_genFromB_dd = 0;
   tree_genFromB_dr = 0;
   tree_genFromB_dz = 0;
   tree_genJet_pt = 0;
   tree_genJet_eta = 0;
   tree_genJet_phi = 0;
   tree_genJet_mass = 0;
   tree_genJet_energy = 0;
   tree_LLP = 0;
   tree_LLP_pt = 0;
   tree_LLP_eta = 0;
   tree_LLP_phi = 0;
   tree_LLP_x = 0;
   tree_LLP_y = 0;
   tree_LLP_z = 0;
   tree_LLP_r = 0;
   tree_LLP_dist = 0;
   tree_LLP_nTrks = 0;
   tree_LLP_Vtx_nTrks = 0;
   tree_LLP_Vtx_NChi2 = 0;
   tree_LLP_Vtx_dx = 0;
   tree_LLP_Vtx_dy = 0;
   tree_LLP_Vtx_dz = 0;
   tree_LLP_Vtx_dist = 0;
   tree_LLP_Vtx_dd = 0;
   tree_LLP12_dR = 0;
   tree_LLP12_deta = 0;
   tree_LLP12_dphi = 0;
   tree_LLP_Mass = 0;
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
   tree_Hemi_LooseBTag_axes = 0;
   tree_Hemi_MediumBTag_axes = 0;
   tree_Hemi_TightBTag_axes = 0;
   tree_Hemi_dR12 = 0;
   tree_Hemi_LLP = 0;
   tree_Hemi_LLP_pt = 0;
   tree_Hemi_LLP_eta = 0;
   tree_Hemi_LLP_phi = 0;
   tree_Hemi_LLP_dist = 0;
   tree_Hemi_LLP_x = 0;
   tree_Hemi_LLP_y = 0;
   tree_Hemi_LLP_z = 0;
   tree_Hemi_LLP_dR = 0;
   tree_Hemi_LLP_mother = 0;
   tree_Hemi_LLP_Vtx_dx = 0;
   tree_Hemi_LLP_Vtx_dy = 0;
   tree_Hemi_LLP_Vtx_dz = 0;
   tree_Hemi_LLP_Vtx_dr = 0;
   tree_Hemi_LLP_muOK_dR = 0;
   tree_Hemi_LLP_muOK_pt = 0;
   tree_Hemi_LLP_muOK_mass = 0;
   tree_Hemi_LLP_muNO_dR = 0;
   tree_Hemi_LLP_muNO_pt = 0;
   tree_Hemi_LLP_muNO_mass = 0;
   tree_Hemi_LLP_dR12 = 0;
   tree_Hemi_LLP_ping = 0;
   tree_event_LLP_ping = 0;
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
   tree_Hemi_Vtx_xError = 0;
   tree_Hemi_Vtx_yError = 0;
   tree_Hemi_Vtx_zError = 0;
   tree_Hemi_Vtx_BTag = 0;
   tree_Hemi_Vtx_trackWeight = 0;
   tree_Hemi_Vtx_SumtrackWeight = 0;
   tree_Hemi_Vtx_track_MeanDCA_d = 0;
   tree_Hemi_Vtx_Mass = 0;
   tree_Hemi_Vtx_dist = 0;
   tree_Hemi_Vtx_ntrk10 = 0;
   tree_Hemi_Vtx_ntrk20 = 0;
   tree_event_nVtx = 0;
   tree_event_Vtx_Vtx_dr = 0;
   tree_event_Vtx_Vtx_dz = 0;
   tree_event_Vtx_Vtx_dd = 0;
   tree_event_Vtx_Vtx_reldd = 0;
   tree_event_Vtx_Vtx_dR = 0;
   tree_event_Vtx_Vtx_step = 0;
   tree_Hemi_SecLLP = 0;
   tree_Hemi_LLP_SecVtx_dx = 0;
   tree_Hemi_LLP_SecVtx_dy = 0;
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
   tree_Hemi_SecVtx_trackWeight = 0;
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
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("tree_only_gen_wt", &tree_only_gen_wt, &b_tree_only_gen_wt);
   fChain->SetBranchAddress("tree_event_weight", &tree_event_weight, &b_tree_event_weight);
   fChain->SetBranchAddress("tree_genTop_Weight", &tree_genTop_Weight, &b_tree_genTop_Weight);
   fChain->SetBranchAddress("PUweight", &PUweight, &b_PUweight);
   fChain->SetBranchAddress("PUweight_Up", &PUweight_Up, &b_PUweight_Up);
   fChain->SetBranchAddress("PUweight_Down", &PUweight_Down, &b_PUweight_Down);
   fChain->SetBranchAddress("Prefweight", &Prefweight, &b_Prefweight);
   fChain->SetBranchAddress("PU_events", &PU_events, &b_PU_events);
   fChain->SetBranchAddress("tree_only_tigger_filter", &tree_only_tigger_filter, &b_tree_only_tigger_filter);
fChain->SetBranchAddress("tree_trigger_doublelepton", &tree_trigger_doublelepton, &b_tree_trigger_doublelepton);
   fChain->SetBranchAddress("tree_trigger_singlelepton", &tree_trigger_singlelepton, &b_tree_trigger_singlelepton);
   fChain->SetBranchAddress("tree_Filter", &tree_Filter, &b_tree_Filter);
   fChain->SetBranchAddress("tree_FilterSameSign", &tree_FilterSameSign, &b_tree_FilterSameSign);
   fChain->SetBranchAddress("tree_Good_PV", &tree_Good_PV, &b_tree_Good_PV);
   fChain->SetBranchAddress("tree_muon_GenRecoTriggerMatched", &tree_muon_GenRecoTriggerMatched, &b_tree_muon_GenRecoTriggerMatched);
   fChain->SetBranchAddress("tree_Evts_MVAval", &tree_Evts_MVAval, &b_tree_Evts_MVAval);
   fChain->SetBranchAddress("tree_Evts_MVAvalDY", &tree_Evts_MVAvalDY, &b_tree_Evts_MVAvalDY);
   fChain->SetBranchAddress("tree_Evts_MVAvalTT", &tree_Evts_MVAvalTT, &b_tree_Evts_MVAvalTT);
   fChain->SetBranchAddress("tree_bs_PosX", &tree_bs_PosX, &b_tree_bs_PosX);
   fChain->SetBranchAddress("tree_bs_PosY", &tree_bs_PosY, &b_tree_bs_PosY);
   fChain->SetBranchAddress("tree_bs_PosZ", &tree_bs_PosZ, &b_tree_bs_PosZ);
   fChain->SetBranchAddress("tree_nPV", &tree_nPV, &b_tree_nPV);
   fChain->SetBranchAddress("tree_PV_x", &tree_PV_x, &b_tree_PV_x);
   fChain->SetBranchAddress("tree_PV_y", &tree_PV_y, &b_tree_PV_y);
   fChain->SetBranchAddress("tree_PV_z", &tree_PV_z, &b_tree_PV_z);
   fChain->SetBranchAddress("tree_PV_ez", &tree_PV_ez, &b_tree_PV_ez);
   fChain->SetBranchAddress("tree_PV_NChi2", &tree_PV_NChi2, &b_tree_PV_NChi2);
   fChain->SetBranchAddress("tree_PV_ndf", &tree_PV_ndf, &b_tree_PV_ndf);
   fChain->SetBranchAddress("tree_PV_rho", &tree_PV_rho, &b_tree_PV_rho);
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
   fChain->SetBranchAddress("tree_muon_x", &tree_muon_x, &b_tree_muon_x);
   fChain->SetBranchAddress("tree_muon_y", &tree_muon_y, &b_tree_muon_y);
   fChain->SetBranchAddress("tree_muon_z", &tree_muon_z, &b_tree_muon_z);
   fChain->SetBranchAddress("tree_muon_dxy", &tree_muon_dxy, &b_tree_muon_dxy);
   fChain->SetBranchAddress("tree_muon_dxyError", &tree_muon_dxyError, &b_tree_muon_dxyError);
   fChain->SetBranchAddress("tree_muon_dz", &tree_muon_dz, &b_tree_muon_dz);
   fChain->SetBranchAddress("tree_muon_dzError", &tree_muon_dzError, &b_tree_muon_dzError);
   fChain->SetBranchAddress("tree_muon_charge", &tree_muon_charge, &b_tree_muon_charge);
   fChain->SetBranchAddress("tree_muon_isLoose", &tree_muon_isLoose, &b_tree_muon_isLoose);
   fChain->SetBranchAddress("tree_muon_isMedium", &tree_muon_isMedium, &b_tree_muon_isMedium);
   fChain->SetBranchAddress("tree_muon_isTight", &tree_muon_isTight, &b_tree_muon_isTight);
   fChain->SetBranchAddress("tree_muon_isGlobal", &tree_muon_isGlobal, &b_tree_muon_isGlobal);
   fChain->SetBranchAddress("tree_muon_isoR3", &tree_muon_isoR3, &b_tree_muon_isoR3);
   fChain->SetBranchAddress("tree_muon_trigger_dimu", &tree_muon_trigger_dimu, &b_tree_muon_trigger_dimu);
   fChain->SetBranchAddress("tree_muon_trigger_isomu", &tree_muon_trigger_isomu, &b_tree_muon_trigger_isomu);
   fChain->SetBranchAddress("tree_muon_PFIsoVeryLoose", &tree_muon_PFIsoVeryLoose, &b_tree_muon_PFIsoVeryLoose);
   fChain->SetBranchAddress("tree_muon_PFIsoLoose", &tree_muon_PFIsoLoose, &b_tree_muon_PFIsoLoose);
   fChain->SetBranchAddress("tree_muon_PFIsoMedium", &tree_muon_PFIsoMedium, &b_tree_muon_PFIsoMedium);
   fChain->SetBranchAddress("tree_muon_PFIsoTight", &tree_muon_PFIsoTight, &b_tree_muon_PFIsoTight);
   fChain->SetBranchAddress("tree_muon_MiniIsoLoose", &tree_muon_MiniIsoLoose, &b_tree_muon_MiniIsoLoose);
   fChain->SetBranchAddress("tree_muon_MiniIsoMedium", &tree_muon_MiniIsoMedium, &b_tree_muon_MiniIsoMedium);
   fChain->SetBranchAddress("tree_muon_MiniIsoTight", &tree_muon_MiniIsoTight, &b_tree_muon_MiniIsoTight);
   fChain->SetBranchAddress("tree_muon_TkIsoLoose", &tree_muon_TkIsoLoose, &b_tree_muon_TkIsoLoose);
   fChain->SetBranchAddress("tree_muon_TkIsoTight", &tree_muon_TkIsoTight, &b_tree_muon_TkIsoTight);
   fChain->SetBranchAddress("tree_muon_trkLayers", &tree_muon_trkLayers, &b_tree_muon_trkLayers);
   fChain->SetBranchAddress("tree_muon_miniIso", &tree_muon_miniIso, &b_tree_muon_miniIso);
   fChain->SetBranchAddress("tree_muon_correction", &tree_muon_correction, &b_tree_muon_correction);
   fChain->SetBranchAddress("tree_muon_gen", &tree_muon_gen, &b_tree_muon_gen);
   fChain->SetBranchAddress("tree_reco_muon_leadingpt", &tree_reco_muon_leadingpt, &b_tree_reco_muon_leadingpt);
   fChain->SetBranchAddress("tree_reco_electron_leadingpt2", &tree_reco_electron_leadingpt2, &b_tree_reco_electron_leadingpt2);
   fChain->SetBranchAddress("tree_reco_muon_leadingeta", &tree_reco_muon_leadingeta, &b_tree_reco_muon_leadingeta);
   fChain->SetBranchAddress("tree_reco_electron_leadingeta2", &tree_reco_electron_leadingeta2, &b_tree_reco_electron_leadingeta2);
   fChain->SetBranchAddress("tree_reco_muon_leadingphi", &tree_reco_muon_leadingphi, &b_tree_reco_muon_leadingphi);
   fChain->SetBranchAddress("tree_reco_electron_leadingphi2", &tree_reco_electron_leadingphi2, &b_tree_reco_electron_leadingphi2);
   fChain->SetBranchAddress("tree_trig_muon_leadingpt", &tree_trig_muon_leadingpt, &b_tree_trig_muon_leadingpt);
   fChain->SetBranchAddress("tree_trig_electron_leadingpt2", &tree_trig_electron_leadingpt2, &b_tree_trig_electron_leadingpt2);
   fChain->SetBranchAddress("tree_trig_muon_leadingeta", &tree_trig_muon_leadingeta, &b_tree_trig_muon_leadingeta);
   fChain->SetBranchAddress("tree_trig_electron_leadingeta2", &tree_trig_electron_leadingeta2, &b_tree_trig_electron_leadingeta2);
   fChain->SetBranchAddress("tree_trig_muon_leadingphi", &tree_trig_muon_leadingphi, &b_tree_trig_muon_leadingphi);
   fChain->SetBranchAddress("tree_trig_electron_leadingphi2", &tree_trig_electron_leadingphi2, &b_tree_trig_electron_leadingphi2);

      fChain->SetBranchAddress("tree_lepton_b4trigger_leadingpt",&tree_lepton_b4trigger_leadingpt,&b_tree_lepton_b4trigger_leadingpt);
   fChain->SetBranchAddress("tree_lepton_b4trigger_leadingpt2",&tree_lepton_b4trigger_leadingpt2,&b_tree_lepton_b4trigger_leadingpt2);


   fChain->SetBranchAddress("tree_reco_lepton_leadingpt", &tree_reco_lepton_leadingpt, &b_tree_reco_lepton_leadingpt);
   fChain->SetBranchAddress("tree_reco_lepton_leadingpt2", &tree_reco_lepton_leadingpt2, &b_tree_reco_lepton_leadingpt2);
   fChain->SetBranchAddress("tree_reco_lepton_leadingeta", &tree_reco_lepton_leadingeta, &b_tree_reco_lepton_leadingeta);
   fChain->SetBranchAddress("tree_reco_lepton_leadingeta2", &tree_reco_lepton_leadingeta2, &b_tree_reco_lepton_leadingeta2);
   fChain->SetBranchAddress("tree_reco_lepton_leadingphi", &tree_reco_lepton_leadingphi, &b_tree_reco_lepton_leadingphi);
   fChain->SetBranchAddress("tree_reco_lepton_leadingphi2", &tree_reco_lepton_leadingphi2, &b_tree_reco_lepton_leadingphi2);
   fChain->SetBranchAddress("tree_trig_lepton_leadingpt", &tree_trig_lepton_leadingpt, &b_tree_trig_lepton_leadingpt);
   fChain->SetBranchAddress("tree_trig_lepton_leadingpt2", &tree_trig_lepton_leadingpt2, &b_tree_trig_lepton_leadingpt2);
   fChain->SetBranchAddress("tree_trig_lepton_leadingeta", &tree_trig_lepton_leadingeta, &b_tree_trig_lepton_leadingeta);
   fChain->SetBranchAddress("tree_trig_lepton_leadingeta2", &tree_trig_lepton_leadingeta2, &b_tree_trig_lepton_leadingeta2);
   fChain->SetBranchAddress("tree_trig_lepton_leadingphi", &tree_trig_lepton_leadingphi, &b_tree_trig_lepton_leadingphi);
   fChain->SetBranchAddress("tree_trig_lepton_leadingphi2", &tree_trig_lepton_leadingphi2, &b_tree_trig_lepton_leadingphi2);
   fChain->SetBranchAddress("tree_lepton_leadingpt", &tree_lepton_leadingpt, &b_tree_lepton_leadingpt);
   fChain->SetBranchAddress("tree_lepton_leadingpt2", &tree_lepton_leadingpt2, &b_tree_lepton_leadingpt2);
   fChain->SetBranchAddress("tree_lepton_leadingeta", &tree_lepton_leadingeta, &b_tree_lepton_leadingeta);
   fChain->SetBranchAddress("tree_lepton_leadingeta2", &tree_lepton_leadingeta2, &b_tree_lepton_leadingeta2);
   fChain->SetBranchAddress("tree_lepton_leadingphi", &tree_lepton_leadingphi, &b_tree_lepton_leadingphi);
   fChain->SetBranchAddress("tree_lepton_leadingphi2", &tree_lepton_leadingphi2, &b_tree_lepton_leadingphi2);
fChain->SetBranchAddress("tree_lepton_leadingdxy", &tree_lepton_leadingdxy, &b_tree_lepton_leadingdxy);
   fChain->SetBranchAddress("tree_lepton_leadingdxy2", &tree_lepton_leadingdxy2, &b_tree_lepton_leadingdxy2);
   fChain->SetBranchAddress("tree_lepton_leadingdz", &tree_lepton_leadingdz, &b_tree_lepton_leadingdz);
   fChain->SetBranchAddress("tree_lepton_leadingdz2", &tree_lepton_leadingdz2, &b_tree_lepton_leadingdz2);
   fChain->SetBranchAddress("tree_lepton_lepton_dR", &tree_lepton_lepton_dR, &b_tree_lepton_lepton_dR);
   fChain->SetBranchAddress("tree_lepton_lepton_dPhi", &tree_lepton_lepton_dPhi, &b_tree_lepton_lepton_dPhi);
   fChain->SetBranchAddress("tree_lepton_lepton_dEta", &tree_lepton_lepton_dEta, &b_tree_lepton_lepton_dEta);
   fChain->SetBranchAddress("tree_ll_pt", &tree_ll_pt, &b_tree_ll_pt);
   fChain->SetBranchAddress("tree_ll_eta", &tree_ll_eta, &b_tree_ll_eta);
   fChain->SetBranchAddress("tree_ll_phi", &tree_ll_phi, &b_tree_ll_phi);
   fChain->SetBranchAddress("tree_ll_px", &tree_ll_px, &b_tree_ll_px);
   fChain->SetBranchAddress("tree_ll_py", &tree_ll_py, &b_tree_ll_py);
   fChain->SetBranchAddress("tree_ll_pz", &tree_ll_pz, &b_tree_ll_pz);
   fChain->SetBranchAddress("tree_ll_energy", &tree_ll_energy, &b_tree_ll_energy);
   fChain->SetBranchAddress("tree_ll_mass", &tree_ll_mass, &b_tree_ll_mass);
   fChain->SetBranchAddress("tree_all_nel", &tree_all_nel, &b_tree_all_nel);
   fChain->SetBranchAddress("tree_electron_nEle", &tree_electron_nEle, &b_tree_electron_nEle);
   fChain->SetBranchAddress("tree_electron_isPrompt", &tree_electron_isPrompt, &b_tree_electron_isPrompt);
   fChain->SetBranchAddress("tree_electron_pt", &tree_electron_pt, &b_tree_electron_pt);
   fChain->SetBranchAddress("tree_electron_eta", &tree_electron_eta, &b_tree_electron_eta);
   fChain->SetBranchAddress("tree_electron_phi", &tree_electron_phi, &b_tree_electron_phi);
   fChain->SetBranchAddress("tree_electron_x", &tree_electron_x, &b_tree_electron_x);
   fChain->SetBranchAddress("tree_electron_y", &tree_electron_y, &b_tree_electron_y);
   fChain->SetBranchAddress("tree_electron_z", &tree_electron_z, &b_tree_electron_z);
   fChain->SetBranchAddress("tree_electron_energy", &tree_electron_energy, &b_tree_electron_energy);
   fChain->SetBranchAddress("tree_electron_et", &tree_electron_et, &b_tree_electron_et);
   fChain->SetBranchAddress("tree_electron_ecal_trk_postcorr", &tree_electron_ecal_trk_postcorr, &b_tree_electron_ecal_trk_postcorr);
   fChain->SetBranchAddress("tree_electron_charge", &tree_electron_charge, &b_tree_electron_charge);
   fChain->SetBranchAddress("tree_electron_isoR4", &tree_electron_isoR4, &b_tree_electron_isoR4);
   fChain->SetBranchAddress("tree_electron_IsLoose", &tree_electron_IsLoose, &b_tree_electron_IsLoose);
   fChain->SetBranchAddress("tree_electron_IsMedium", &tree_electron_IsMedium, &b_tree_electron_IsMedium);
   fChain->SetBranchAddress("tree_electron_IsTight", &tree_electron_IsTight, &b_tree_electron_IsTight);
   fChain->SetBranchAddress("tree_electron_dxy", &tree_electron_dxy, &b_tree_electron_dxy);
   fChain->SetBranchAddress("tree_electron_dz", &tree_electron_dz, &b_tree_electron_dz);
   fChain->SetBranchAddress("tree_electron_gen", &tree_electron_gen, &b_tree_electron_gen);
   fChain->SetBranchAddress("tree_PFMet_et", &tree_PFMet_et, &b_tree_PFMet_et);
   fChain->SetBranchAddress("tree_PFMet_phi", &tree_PFMet_phi, &b_tree_PFMet_phi);
   fChain->SetBranchAddress("tree_PFMet_sig", &tree_PFMet_sig, &b_tree_PFMet_sig);
   fChain->SetBranchAddress("tree_PFMet_pt", &tree_PFMet_pt, &b_tree_PFMet_pt);
   fChain->SetBranchAddress("tree_njet", &tree_njet, &b_tree_njet);
   fChain->SetBranchAddress("tree_njetNOmu", &tree_njetNOmu, &b_tree_njetNOmu);
   fChain->SetBranchAddress("tree_jet_pt", &tree_jet_pt, &b_tree_jet_pt);
   fChain->SetBranchAddress("tree_jet_eta", &tree_jet_eta, &b_tree_jet_eta);
   fChain->SetBranchAddress("tree_jet_phi", &tree_jet_phi, &b_tree_jet_phi);
   fChain->SetBranchAddress("tree_jet_px", &tree_jet_px, &b_tree_jet_px);
   fChain->SetBranchAddress("tree_jet_py", &tree_jet_py, &b_tree_jet_py);
   fChain->SetBranchAddress("tree_jet_pz", &tree_jet_pz, &b_tree_jet_pz);
   fChain->SetBranchAddress("tree_jet_E", &tree_jet_E, &b_tree_jet_E);
   fChain->SetBranchAddress("tree_jet_tightid_LepVeto", &tree_jet_tightid_LepVeto, &b_tree_jet_tightid_LepVeto);
   fChain->SetBranchAddress("tree_jet_tightid", &tree_jet_tightid, &b_tree_jet_tightid);
   fChain->SetBranchAddress("tree_jet_TightJetIDLepVeto", &tree_jet_TightJetIDLepVeto, &b_tree_jet_TightJetIDLepVeto);
   fChain->SetBranchAddress("tree_jet_TightJetID", &tree_jet_TightJetID, &b_tree_jet_TightJetID);
   fChain->SetBranchAddress("tree_jet_HadronFlavour", &tree_jet_HadronFlavour, &b_tree_jet_HadronFlavour);
   fChain->SetBranchAddress("tree_jet_pileupID", &tree_jet_pileupID, &b_tree_jet_pileupID);
   fChain->SetBranchAddress("tree_jet_btag_DeepCSV", &tree_jet_btag_DeepCSV, &b_tree_jet_btag_DeepCSV);
   fChain->SetBranchAddress("tree_jet_btag_DeepJet", &tree_jet_btag_DeepJet, &b_tree_jet_btag_DeepJet);
   fChain->SetBranchAddress("tree_jet_leadingpt", &tree_jet_leadingpt, &b_tree_jet_leadingpt);
   fChain->SetBranchAddress("tree_jet_leadingpt2", &tree_jet_leadingpt2, &b_tree_jet_leadingpt2);
   fChain->SetBranchAddress("tree_jet_leadingeta", &tree_jet_leadingeta, &b_tree_jet_leadingeta);
   fChain->SetBranchAddress("tree_jet_leadingeta2", &tree_jet_leadingeta2, &b_tree_jet_leadingeta2);
   fChain->SetBranchAddress("tree_jet_leadingMuon_dR", &tree_jet_leadingMuon_dR, &b_tree_jet_leadingMuon_dR);
   fChain->SetBranchAddress("tree_jet_leadingMuon2_dR", &tree_jet_leadingMuon2_dR, &b_tree_jet_leadingMuon2_dR);
   fChain->SetBranchAddress("tree_jet_jet_dR", &tree_jet_jet_dR, &b_tree_jet_jet_dR);
   fChain->SetBranchAddress("tree_jet_jet_dPhi", &tree_jet_jet_dPhi, &b_tree_jet_jet_dPhi);
   fChain->SetBranchAddress("tree_jet_jet_dEta", &tree_jet_jet_dEta, &b_tree_jet_jet_dEta);
   fChain->SetBranchAddress("tree_muon_jet_dRmin", &tree_muon_jet_dRmin, &b_tree_muon_jet_dRmin);
   fChain->SetBranchAddress("tree_muon_jet_dRmax", &tree_muon_jet_dRmax, &b_tree_muon_jet_dRmax);
   fChain->SetBranchAddress("tree_elemu_jet_dRmin", &tree_elemu_jet_dRmin, &b_tree_elemu_jet_dRmin);
   fChain->SetBranchAddress("tree_elemu_jet_dRmax", &tree_elemu_jet_dRmax, &b_tree_elemu_jet_dRmax);
   fChain->SetBranchAddress("tree_ele_jet_dRmin", &tree_ele_jet_dRmin, &b_tree_ele_jet_dRmin);
   fChain->SetBranchAddress("tree_ele_jet_dRmax", &tree_ele_jet_dRmax, &b_tree_ele_jet_dRmax);
   fChain->SetBranchAddress("tree_HT", &tree_HT, &b_tree_HT);
   fChain->SetBranchAddress("tree_nK0", &tree_nK0, &b_tree_nK0);
   fChain->SetBranchAddress("tree_K0_x", &tree_K0_x, &b_tree_K0_x);
   fChain->SetBranchAddress("tree_K0_y", &tree_K0_y, &b_tree_K0_y);
   fChain->SetBranchAddress("tree_K0_z", &tree_K0_z, &b_tree_K0_z);
   fChain->SetBranchAddress("tree_K0_r", &tree_K0_r, &b_tree_K0_r);
   fChain->SetBranchAddress("tree_K0_NChi2", &tree_K0_NChi2, &b_tree_K0_NChi2);
   fChain->SetBranchAddress("tree_K0_ndf", &tree_K0_ndf, &b_tree_K0_ndf);
   fChain->SetBranchAddress("tree_K0_mass", &tree_K0_mass, &b_tree_K0_mass);
   fChain->SetBranchAddress("tree_K0_pt", &tree_K0_pt, &b_tree_K0_pt);
   fChain->SetBranchAddress("tree_K0_eta", &tree_K0_eta, &b_tree_K0_eta);
   fChain->SetBranchAddress("tree_K0_phi", &tree_K0_phi, &b_tree_K0_phi);
   fChain->SetBranchAddress("tree_K0_nDaughters", &tree_K0_nDaughters, &b_tree_K0_nDaughters);
   fChain->SetBranchAddress("tree_nLambda", &tree_nLambda, &b_tree_nLambda);
   fChain->SetBranchAddress("tree_L0_x", &tree_L0_x, &b_tree_L0_x);
   fChain->SetBranchAddress("tree_L0_y", &tree_L0_y, &b_tree_L0_y);
   fChain->SetBranchAddress("tree_L0_z", &tree_L0_z, &b_tree_L0_z);
   fChain->SetBranchAddress("tree_L0_r", &tree_L0_r, &b_tree_L0_r);
   fChain->SetBranchAddress("tree_L0_NChi2", &tree_L0_NChi2, &b_tree_L0_NChi2);
   fChain->SetBranchAddress("tree_L0_ndf", &tree_L0_ndf, &b_tree_L0_ndf);
   fChain->SetBranchAddress("tree_L0_nDaughters", &tree_L0_nDaughters, &b_tree_L0_nDaughters);
   fChain->SetBranchAddress("tree_L0_mass", &tree_L0_mass, &b_tree_L0_mass);
   fChain->SetBranchAddress("tree_L0_pt", &tree_L0_pt, &b_tree_L0_pt);
   fChain->SetBranchAddress("tree_L0_eta", &tree_L0_eta, &b_tree_L0_eta);
   fChain->SetBranchAddress("tree_L0_phi", &tree_L0_phi, &b_tree_L0_phi);
   fChain->SetBranchAddress("tree_nV0_reco", &tree_nV0_reco, &b_tree_nV0_reco);
   fChain->SetBranchAddress("tree_V0_reco_x", &tree_V0_reco_x, &b_tree_V0_reco_x);
   fChain->SetBranchAddress("tree_V0_reco_y", &tree_V0_reco_y, &b_tree_V0_reco_y);
   fChain->SetBranchAddress("tree_V0_reco_z", &tree_V0_reco_z, &b_tree_V0_reco_z);
   fChain->SetBranchAddress("tree_V0_reco_r", &tree_V0_reco_r, &b_tree_V0_reco_r);
   fChain->SetBranchAddress("tree_V0_reco_drSig", &tree_V0_reco_drSig, &b_tree_V0_reco_drSig);
   fChain->SetBranchAddress("tree_V0_reco_dzSig", &tree_V0_reco_dzSig, &b_tree_V0_reco_dzSig);
   fChain->SetBranchAddress("tree_V0_reco_angleXY", &tree_V0_reco_angleXY, &b_tree_V0_reco_angleXY);
   fChain->SetBranchAddress("tree_V0_reco_angleZ", &tree_V0_reco_angleZ, &b_tree_V0_reco_angleZ);
   fChain->SetBranchAddress("tree_V0_reco_NChi2", &tree_V0_reco_NChi2, &b_tree_V0_reco_NChi2);
   fChain->SetBranchAddress("tree_V0_reco_ndf", &tree_V0_reco_ndf, &b_tree_V0_reco_ndf);
   fChain->SetBranchAddress("tree_V0_reco_mass", &tree_V0_reco_mass, &b_tree_V0_reco_mass);
   fChain->SetBranchAddress("tree_V0_reco_pt", &tree_V0_reco_pt, &b_tree_V0_reco_pt);
   fChain->SetBranchAddress("tree_V0_reco_eta", &tree_V0_reco_eta, &b_tree_V0_reco_eta);
   fChain->SetBranchAddress("tree_V0_reco_phi", &tree_V0_reco_phi, &b_tree_V0_reco_phi);
   fChain->SetBranchAddress("tree_V0_reco_source", &tree_V0_reco_source, &b_tree_V0_reco_source);
   fChain->SetBranchAddress("tree_V0_reco_badTkHit", &tree_V0_reco_badTkHit, &b_tree_V0_reco_badTkHit);
   fChain->SetBranchAddress("tree_V0_reco_dca", &tree_V0_reco_dca, &b_tree_V0_reco_dca);
   fChain->SetBranchAddress("tree_nSecInt", &tree_nSecInt, &b_tree_nSecInt);
   fChain->SetBranchAddress("tree_SecInt_x", &tree_SecInt_x, &b_tree_SecInt_x);
   fChain->SetBranchAddress("tree_SecInt_y", &tree_SecInt_y, &b_tree_SecInt_y);
   fChain->SetBranchAddress("tree_SecInt_z", &tree_SecInt_z, &b_tree_SecInt_z);
   fChain->SetBranchAddress("tree_SecInt_r", &tree_SecInt_r, &b_tree_SecInt_r);
   fChain->SetBranchAddress("tree_SecInt_d", &tree_SecInt_d, &b_tree_SecInt_d);
   fChain->SetBranchAddress("tree_SecInt_drSig", &tree_SecInt_drSig, &b_tree_SecInt_drSig);
   fChain->SetBranchAddress("tree_SecInt_dzSig", &tree_SecInt_dzSig, &b_tree_SecInt_dzSig);
   fChain->SetBranchAddress("tree_SecInt_angleXY", &tree_SecInt_angleXY, &b_tree_SecInt_angleXY);
   fChain->SetBranchAddress("tree_SecInt_angleZ", &tree_SecInt_angleZ, &b_tree_SecInt_angleZ);
   fChain->SetBranchAddress("tree_SecInt_NChi2", &tree_SecInt_NChi2, &b_tree_SecInt_NChi2);
   fChain->SetBranchAddress("tree_SecInt_ndf", &tree_SecInt_ndf, &b_tree_SecInt_ndf);
   fChain->SetBranchAddress("tree_SecInt_mass", &tree_SecInt_mass, &b_tree_SecInt_mass);
   fChain->SetBranchAddress("tree_SecInt_pt", &tree_SecInt_pt, &b_tree_SecInt_pt);
   fChain->SetBranchAddress("tree_SecInt_eta", &tree_SecInt_eta, &b_tree_SecInt_eta);
   fChain->SetBranchAddress("tree_SecInt_phi", &tree_SecInt_phi, &b_tree_SecInt_phi);
   fChain->SetBranchAddress("tree_SecInt_charge", &tree_SecInt_charge, &b_tree_SecInt_charge);
   fChain->SetBranchAddress("tree_SecInt_badTkHit", &tree_SecInt_badTkHit, &b_tree_SecInt_badTkHit);
   fChain->SetBranchAddress("tree_SecInt_dca", &tree_SecInt_dca, &b_tree_SecInt_dca);
   fChain->SetBranchAddress("tree_SecInt_selec", &tree_SecInt_selec, &b_tree_SecInt_selec);
   fChain->SetBranchAddress("tree_SecInt_layer", &tree_SecInt_layer, &b_tree_SecInt_layer);
   fChain->SetBranchAddress("tree_SecInt_LLP", &tree_SecInt_LLP, &b_tree_SecInt_LLP);
   fChain->SetBranchAddress("tree_SecInt_LLP_dr", &tree_SecInt_LLP_dr, &b_tree_SecInt_LLP_dr);
   fChain->SetBranchAddress("tree_SecInt_LLP_dz", &tree_SecInt_LLP_dz, &b_tree_SecInt_LLP_dz);
   fChain->SetBranchAddress("tree_SecInt_LLP_dd", &tree_SecInt_LLP_dd, &b_tree_SecInt_LLP_dd);
   fChain->SetBranchAddress("tree_SecInt_tk1", &tree_SecInt_tk1, &b_tree_SecInt_tk1);
   fChain->SetBranchAddress("tree_SecInt_tk2", &tree_SecInt_tk2, &b_tree_SecInt_tk2);
   fChain->SetBranchAddress("tree_nYConv", &tree_nYConv, &b_tree_nYConv);
   fChain->SetBranchAddress("tree_Yc_x", &tree_Yc_x, &b_tree_Yc_x);
   fChain->SetBranchAddress("tree_Yc_y", &tree_Yc_y, &b_tree_Yc_y);
   fChain->SetBranchAddress("tree_Yc_z", &tree_Yc_z, &b_tree_Yc_z);
   fChain->SetBranchAddress("tree_Yc_r", &tree_Yc_r, &b_tree_Yc_r);
   fChain->SetBranchAddress("tree_Yc_dr0", &tree_Yc_dr0, &b_tree_Yc_dr0);
   fChain->SetBranchAddress("tree_Yc_dr1", &tree_Yc_dr1, &b_tree_Yc_dr1);
   fChain->SetBranchAddress("tree_Yc_dz0", &tree_Yc_dz0, &b_tree_Yc_dz0);
   fChain->SetBranchAddress("tree_Yc_dz1", &tree_Yc_dz1, &b_tree_Yc_dz1);
   fChain->SetBranchAddress("tree_Yc_costheta", &tree_Yc_costheta, &b_tree_Yc_costheta);
   fChain->SetBranchAddress("tree_Yc_layer", &tree_Yc_layer, &b_tree_Yc_layer);
   fChain->SetBranchAddress("tree_Yc_NChi2", &tree_Yc_NChi2, &b_tree_Yc_NChi2);
   fChain->SetBranchAddress("tree_Yc_ndf", &tree_Yc_ndf, &b_tree_Yc_ndf);
   fChain->SetBranchAddress("tree_Yc_nDaughters", &tree_Yc_nDaughters, &b_tree_Yc_nDaughters);
   fChain->SetBranchAddress("tree_Yc_pt", &tree_Yc_pt, &b_tree_Yc_pt);
   fChain->SetBranchAddress("tree_Yc_eta", &tree_Yc_eta, &b_tree_Yc_eta);
   fChain->SetBranchAddress("tree_Yc_phi", &tree_Yc_phi, &b_tree_Yc_phi);
   fChain->SetBranchAddress("tree_Yc_mass", &tree_Yc_mass, &b_tree_Yc_mass);
   fChain->SetBranchAddress("tree_Yc_ntracks", &tree_Yc_ntracks, &b_tree_Yc_ntracks);
   fChain->SetBranchAddress("tree_Yc_tracks_index", &tree_Yc_tracks_index, &b_tree_Yc_tracks_index);
   fChain->SetBranchAddress("tree_Yc_tracks_charge", &tree_Yc_tracks_charge, &b_tree_Yc_tracks_charge);
   fChain->SetBranchAddress("tree_Yc_tracks_pt", &tree_Yc_tracks_pt, &b_tree_Yc_tracks_pt);
   fChain->SetBranchAddress("tree_Yc_tracks_eta", &tree_Yc_tracks_eta, &b_tree_Yc_tracks_eta);
   fChain->SetBranchAddress("tree_Yc_tracks_phi", &tree_Yc_tracks_phi, &b_tree_Yc_tracks_phi);
   fChain->SetBranchAddress("tree_Yc_tracks_phi0", &tree_Yc_tracks_phi0, &b_tree_Yc_tracks_phi0);
   fChain->SetBranchAddress("tree_TRACK_SIZE", &tree_TRACK_SIZE, &b_tree_TRACK_SIZE);
   fChain->SetBranchAddress("tree_nTracks", &tree_nTracks, &b_tree_nTracks);
   fChain->SetBranchAddress("tree_nLostTracks", &tree_nLostTracks, &b_tree_nLostTracks);
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
   fChain->SetBranchAddress("tree_track_btag", &tree_track_btag, &b_tree_track_btag);
   fChain->SetBranchAddress("tree_track_energy", &tree_track_energy, &b_tree_track_energy);
   fChain->SetBranchAddress("tree_track_Hemi", &tree_track_Hemi, &b_tree_track_Hemi);
   fChain->SetBranchAddress("tree_track_Hemi_dR", &tree_track_Hemi_dR, &b_tree_track_Hemi_dR);
   fChain->SetBranchAddress("tree_track_Hemi_dRmax", &tree_track_Hemi_dRmax, &b_tree_track_Hemi_dRmax);
   fChain->SetBranchAddress("tree_track_Hemi_mva_NChi2", &tree_track_Hemi_mva_NChi2, &b_tree_track_Hemi_mva_NChi2);
   fChain->SetBranchAddress("tree_track_Hemi_ping", &tree_track_Hemi_ping, &b_tree_track_Hemi_ping);
   fChain->SetBranchAddress("tree_track_Hemi_dFirstVtx", &tree_track_Hemi_dFirstVtx, &b_tree_track_Hemi_dFirstVtx);
   fChain->SetBranchAddress("tree_track_Hemi_LLP", &tree_track_Hemi_LLP, &b_tree_track_Hemi_LLP);
   fChain->SetBranchAddress("tree_track_Hemi_isjet", &tree_track_Hemi_isjet, &b_tree_track_Hemi_isjet);
   fChain->SetBranchAddress("tree_track_sim_LLP", &tree_track_sim_LLP, &b_tree_track_sim_LLP);
   fChain->SetBranchAddress("tree_track_sim_isFromB", &tree_track_sim_isFromB, &b_tree_track_sim_isFromB);
   fChain->SetBranchAddress("tree_track_sim_isFromC", &tree_track_sim_isFromC, &b_tree_track_sim_isFromC);
   fChain->SetBranchAddress("tree_track_sim_pt", &tree_track_sim_pt, &b_tree_track_sim_pt);
   fChain->SetBranchAddress("tree_track_sim_eta", &tree_track_sim_eta, &b_tree_track_sim_eta);
   fChain->SetBranchAddress("tree_track_sim_phi", &tree_track_sim_phi, &b_tree_track_sim_phi);
   fChain->SetBranchAddress("tree_track_sim_charge", &tree_track_sim_charge, &b_tree_track_sim_charge);
   fChain->SetBranchAddress("tree_track_sim_pdgId", &tree_track_sim_pdgId, &b_tree_track_sim_pdgId);
   fChain->SetBranchAddress("tree_track_sim_mass", &tree_track_sim_mass, &b_tree_track_sim_mass);
   fChain->SetBranchAddress("tree_track_sim_x", &tree_track_sim_x, &b_tree_track_sim_x);
   fChain->SetBranchAddress("tree_track_sim_y", &tree_track_sim_y, &b_tree_track_sim_y);
   fChain->SetBranchAddress("tree_track_sim_z", &tree_track_sim_z, &b_tree_track_sim_z);
   fChain->SetBranchAddress("tree_track_sim_dFirstGen", &tree_track_sim_dFirstGen, &b_tree_track_sim_dFirstGen);
   fChain->SetBranchAddress("tree_track_sim_LLP_r", &tree_track_sim_LLP_r, &b_tree_track_sim_LLP_r);
   fChain->SetBranchAddress("tree_track_sim_LLP_z", &tree_track_sim_LLP_z, &b_tree_track_sim_LLP_z);
   fChain->SetBranchAddress("tree_V0_track_isFromV0", &tree_V0_track_isFromV0, &b_tree_V0_track_isFromV0);
   fChain->SetBranchAddress("tree_V0_track_isFromSI", &tree_V0_track_isFromSI, &b_tree_V0_track_isFromSI);
   fChain->SetBranchAddress("tree_V0_track_lost", &tree_V0_track_lost, &b_tree_V0_track_lost);
   fChain->SetBranchAddress("tree_V0_track_pt", &tree_V0_track_pt, &b_tree_V0_track_pt);
   fChain->SetBranchAddress("tree_V0_track_eta", &tree_V0_track_eta, &b_tree_V0_track_eta);
   fChain->SetBranchAddress("tree_V0_track_phi", &tree_V0_track_phi, &b_tree_V0_track_phi);
   fChain->SetBranchAddress("tree_V0_track_charge", &tree_V0_track_charge, &b_tree_V0_track_charge);
   fChain->SetBranchAddress("tree_V0_track_NChi2", &tree_V0_track_NChi2, &b_tree_V0_track_NChi2);
   fChain->SetBranchAddress("tree_V0_track_dxy", &tree_V0_track_dxy, &b_tree_V0_track_dxy);
   fChain->SetBranchAddress("tree_V0_track_drSig", &tree_V0_track_drSig, &b_tree_V0_track_drSig);
   fChain->SetBranchAddress("tree_V0_track_dz", &tree_V0_track_dz, &b_tree_V0_track_dz);
   fChain->SetBranchAddress("tree_V0_track_dzSig", &tree_V0_track_dzSig, &b_tree_V0_track_dzSig);
   fChain->SetBranchAddress("tree_V0_track_nHit", &tree_V0_track_nHit, &b_tree_V0_track_nHit);
   fChain->SetBranchAddress("tree_V0_track_nHitPixel", &tree_V0_track_nHitPixel, &b_tree_V0_track_nHitPixel);
   fChain->SetBranchAddress("tree_V0_track_firstHit", &tree_V0_track_firstHit, &b_tree_V0_track_firstHit);
   fChain->SetBranchAddress("tree_V0_track_firstHit_x", &tree_V0_track_firstHit_x, &b_tree_V0_track_firstHit_x);
   fChain->SetBranchAddress("tree_V0_track_firstHit_y", &tree_V0_track_firstHit_y, &b_tree_V0_track_firstHit_y);
   fChain->SetBranchAddress("tree_V0_track_firstHit_z", &tree_V0_track_firstHit_z, &b_tree_V0_track_firstHit_z);
   fChain->SetBranchAddress("tree_V0_track_iJet", &tree_V0_track_iJet, &b_tree_V0_track_iJet);
   fChain->SetBranchAddress("tree_V0_track_ntrk10", &tree_V0_track_ntrk10, &b_tree_V0_track_ntrk10);
   fChain->SetBranchAddress("tree_V0_track_ntrk20", &tree_V0_track_ntrk20, &b_tree_V0_track_ntrk20);
   fChain->SetBranchAddress("tree_V0_track_ntrk30", &tree_V0_track_ntrk30, &b_tree_V0_track_ntrk30);
   fChain->SetBranchAddress("tree_V0_track_ntrk40", &tree_V0_track_ntrk40, &b_tree_V0_track_ntrk40);
   fChain->SetBranchAddress("tree_V0_track_Hemi", &tree_V0_track_Hemi, &b_tree_V0_track_Hemi);
   fChain->SetBranchAddress("tree_V0_track_Hemi_dR", &tree_V0_track_Hemi_dR, &b_tree_V0_track_Hemi_dR);
   fChain->SetBranchAddress("tree_V0_track_Hemi_dRmax", &tree_V0_track_Hemi_dRmax, &b_tree_V0_track_Hemi_dRmax);
   fChain->SetBranchAddress("tree_GenPVx", &tree_GenPVx, &b_tree_GenPVx);
   fChain->SetBranchAddress("tree_GenPVy", &tree_GenPVy, &b_tree_GenPVy);
   fChain->SetBranchAddress("tree_GenPVz", &tree_GenPVz, &b_tree_GenPVz);
   fChain->SetBranchAddress("tree_smu_mass", &tree_smu_mass, &b_tree_smu_mass);
   fChain->SetBranchAddress("tree_neu_mass", &tree_neu_mass, &b_tree_neu_mass);
   fChain->SetBranchAddress("tree_neu_ctau", &tree_neu_ctau, &b_tree_neu_ctau);
   fChain->SetBranchAddress("tree_genParticle_pt", &tree_genParticle_pt, &b_tree_genParticle_pt);
   fChain->SetBranchAddress("tree_genParticle_eta", &tree_genParticle_eta, &b_tree_genParticle_eta);
   fChain->SetBranchAddress("tree_genParticle_phi", &tree_genParticle_phi, &b_tree_genParticle_phi);
   fChain->SetBranchAddress("tree_genParticle_charge", &tree_genParticle_charge, &b_tree_genParticle_charge);
   fChain->SetBranchAddress("tree_genParticle_pdgId", &tree_genParticle_pdgId, &b_tree_genParticle_pdgId);
   fChain->SetBranchAddress("tree_genParticle_mass", &tree_genParticle_mass, &b_tree_genParticle_mass);
   fChain->SetBranchAddress("tree_genParticle_x", &tree_genParticle_x, &b_tree_genParticle_x);
   fChain->SetBranchAddress("tree_genParticle_y", &tree_genParticle_y, &b_tree_genParticle_y);
   fChain->SetBranchAddress("tree_genParticle_z", &tree_genParticle_z, &b_tree_genParticle_z);
   fChain->SetBranchAddress("tree_genParticle_px", &tree_genParticle_px, &b_tree_genParticle_px);
   fChain->SetBranchAddress("tree_genParticle_py", &tree_genParticle_py, &b_tree_genParticle_py);
   fChain->SetBranchAddress("tree_genParticle_pz", &tree_genParticle_pz, &b_tree_genParticle_pz);
   fChain->SetBranchAddress("tree_genParticle_energy", &tree_genParticle_energy, &b_tree_genParticle_energy);
   fChain->SetBranchAddress("tree_genParticle_isPromptFinalState", &tree_genParticle_isPromptFinalState, &b_tree_genParticle_isPromptFinalState);
   fChain->SetBranchAddress("tree_genParticle_statusCode", &tree_genParticle_statusCode, &b_tree_genParticle_statusCode);
   fChain->SetBranchAddress("tree_genParticle_mother_pdgId", &tree_genParticle_mother_pdgId, &b_tree_genParticle_mother_pdgId);
   fChain->SetBranchAddress("tree_genParticle_LLP", &tree_genParticle_LLP, &b_tree_genParticle_LLP);
   fChain->SetBranchAddress("tree_genParticle_ct", &tree_genParticle_ct, &b_tree_genParticle_ct);
   fChain->SetBranchAddress("tree_genParticle_ct0", &tree_genParticle_ct0, &b_tree_genParticle_ct0);
   fChain->SetBranchAddress("tree_ngenPackPart", &tree_ngenPackPart, &b_tree_ngenPackPart);
   fChain->SetBranchAddress("tree_genPackPart_pt", &tree_genPackPart_pt, &b_tree_genPackPart_pt);
   fChain->SetBranchAddress("tree_genPackPart_eta", &tree_genPackPart_eta, &b_tree_genPackPart_eta);
   fChain->SetBranchAddress("tree_genPackPart_phi", &tree_genPackPart_phi, &b_tree_genPackPart_phi);
   fChain->SetBranchAddress("tree_genPackPart_charge", &tree_genPackPart_charge, &b_tree_genPackPart_charge);
   fChain->SetBranchAddress("tree_genPackPart_pdgId", &tree_genPackPart_pdgId, &b_tree_genPackPart_pdgId);
   fChain->SetBranchAddress("tree_genPackPart_mass", &tree_genPackPart_mass, &b_tree_genPackPart_mass);
   fChain->SetBranchAddress("tree_genPackPart_x", &tree_genPackPart_x, &b_tree_genPackPart_x);
   fChain->SetBranchAddress("tree_genPackPart_y", &tree_genPackPart_y, &b_tree_genPackPart_y);
   fChain->SetBranchAddress("tree_genPackPart_z", &tree_genPackPart_z, &b_tree_genPackPart_z);
   fChain->SetBranchAddress("tree_genPackPart_mother_pdgId", &tree_genPackPart_mother_pdgId, &b_tree_genPackPart_mother_pdgId);
   fChain->SetBranchAddress("tree_genPackPart_isFromB", &tree_genPackPart_isFromB, &b_tree_genPackPart_isFromB);
   fChain->SetBranchAddress("tree_genPackPart_isFromC", &tree_genPackPart_isFromC, &b_tree_genPackPart_isFromC);
   fChain->SetBranchAddress("tree_ngenFromLLP", &tree_ngenFromLLP, &b_tree_ngenFromLLP);
   fChain->SetBranchAddress("tree_genFromLLP_LLP", &tree_genFromLLP_LLP, &b_tree_genFromLLP_LLP);
   fChain->SetBranchAddress("tree_genFromLLP_pt", &tree_genFromLLP_pt, &b_tree_genFromLLP_pt);
   fChain->SetBranchAddress("tree_genFromLLP_eta", &tree_genFromLLP_eta, &b_tree_genFromLLP_eta);
   fChain->SetBranchAddress("tree_genFromLLP_phi", &tree_genFromLLP_phi, &b_tree_genFromLLP_phi);
   fChain->SetBranchAddress("tree_genFromLLP_charge", &tree_genFromLLP_charge, &b_tree_genFromLLP_charge);
   fChain->SetBranchAddress("tree_genFromLLP_pdgId", &tree_genFromLLP_pdgId, &b_tree_genFromLLP_pdgId);
   fChain->SetBranchAddress("tree_genFromLLP_mass", &tree_genFromLLP_mass, &b_tree_genFromLLP_mass);
   fChain->SetBranchAddress("tree_genFromLLP_x", &tree_genFromLLP_x, &b_tree_genFromLLP_x);
   fChain->SetBranchAddress("tree_genFromLLP_y", &tree_genFromLLP_y, &b_tree_genFromLLP_y);
   fChain->SetBranchAddress("tree_genFromLLP_z", &tree_genFromLLP_z, &b_tree_genFromLLP_z);
   fChain->SetBranchAddress("tree_genFromLLP_mother_pdgId", &tree_genFromLLP_mother_pdgId, &b_tree_genFromLLP_mother_pdgId);
   fChain->SetBranchAddress("tree_genFromLLP_isFromB", &tree_genFromLLP_isFromB, &b_tree_genFromLLP_isFromB);
   fChain->SetBranchAddress("tree_genFromLLP_isFromC", &tree_genFromLLP_isFromC, &b_tree_genFromLLP_isFromC);
fChain->SetBranchAddress("tree_gen_top_pt", &tree_gen_top_pt, &b_tree_gen_top_pt);
   fChain->SetBranchAddress("tree_gen_top_rw_pt", &tree_gen_top_rw_pt, &b_tree_gen_top_rw_pt);
   fChain->SetBranchAddress("tree_genAxis_dRneuneu", &tree_genAxis_dRneuneu, &b_tree_genAxis_dRneuneu);
   fChain->SetBranchAddress("tree_genAxis_dPhineuneu", &tree_genAxis_dPhineuneu, &b_tree_genAxis_dPhineuneu);
   fChain->SetBranchAddress("tree_genAxis_dEtaneuneu", &tree_genAxis_dEtaneuneu, &b_tree_genAxis_dEtaneuneu);
   fChain->SetBranchAddress("tree_GenAxes_Mass", &tree_GenAxes_Mass, &b_tree_GenAxes_Mass);
   fChain->SetBranchAddress("tree_GenAxes_CombinedHemiLeptonMass", &tree_GenAxes_CombinedHemiLeptonMass, &b_tree_GenAxes_CombinedHemiLeptonMass);
   fChain->SetBranchAddress("tree_GenAxis_Neu_dRmin", &tree_GenAxis_Neu_dRmin, &b_tree_GenAxis_Neu_dRmin);
   fChain->SetBranchAddress("tree_GenAxis_Neu_dRmax", &tree_GenAxis_Neu_dRmax, &b_tree_GenAxis_Neu_dRmax);
   fChain->SetBranchAddress("tree_GenAxis_RecoAxis_dRmin", &tree_GenAxis_RecoAxis_dRmin, &b_tree_GenAxis_RecoAxis_dRmin);
   fChain->SetBranchAddress("tree_GenAxis_RecoAxis_dRmax", &tree_GenAxis_RecoAxis_dRmax, &b_tree_GenAxis_RecoAxis_dRmax);
   fChain->SetBranchAddress("tree_nFromC", &tree_nFromC, &b_tree_nFromC);
   fChain->SetBranchAddress("tree_genFromC_pt", &tree_genFromC_pt, &b_tree_genFromC_pt);
   fChain->SetBranchAddress("tree_genFromC_eta", &tree_genFromC_eta, &b_tree_genFromC_eta);
   fChain->SetBranchAddress("tree_genFromC_phi", &tree_genFromC_phi, &b_tree_genFromC_phi);
   fChain->SetBranchAddress("tree_genFromC_charge", &tree_genFromC_charge, &b_tree_genFromC_charge);
   fChain->SetBranchAddress("tree_genFromC_pdgId", &tree_genFromC_pdgId, &b_tree_genFromC_pdgId);
   fChain->SetBranchAddress("tree_genFromC_x", &tree_genFromC_x, &b_tree_genFromC_x);
   fChain->SetBranchAddress("tree_genFromC_y", &tree_genFromC_y, &b_tree_genFromC_y);
   fChain->SetBranchAddress("tree_genFromC_z", &tree_genFromC_z, &b_tree_genFromC_z);
   fChain->SetBranchAddress("tree_genFromC_mother_pdgId", &tree_genFromC_mother_pdgId, &b_tree_genFromC_mother_pdgId);
   fChain->SetBranchAddress("tree_nFromB", &tree_nFromB, &b_tree_nFromB);
   fChain->SetBranchAddress("tree_genFromB_pt", &tree_genFromB_pt, &b_tree_genFromB_pt);
   fChain->SetBranchAddress("tree_genFromB_eta", &tree_genFromB_eta, &b_tree_genFromB_eta);
   fChain->SetBranchAddress("tree_genFromB_phi", &tree_genFromB_phi, &b_tree_genFromB_phi);
   fChain->SetBranchAddress("tree_genFromB_charge", &tree_genFromB_charge, &b_tree_genFromB_charge);
   fChain->SetBranchAddress("tree_genFromB_pdgId", &tree_genFromB_pdgId, &b_tree_genFromB_pdgId);
   fChain->SetBranchAddress("tree_genFromB_x", &tree_genFromB_x, &b_tree_genFromB_x);
   fChain->SetBranchAddress("tree_genFromB_y", &tree_genFromB_y, &b_tree_genFromB_y);
   fChain->SetBranchAddress("tree_genFromB_z", &tree_genFromB_z, &b_tree_genFromB_z);
   fChain->SetBranchAddress("tree_genFromB_mother_pdgId", &tree_genFromB_mother_pdgId, &b_tree_genFromB_mother_pdgId);
   fChain->SetBranchAddress("tree_genFromB_dd", &tree_genFromB_dd, &b_tree_genFromB_dd);
   fChain->SetBranchAddress("tree_genFromB_dr", &tree_genFromB_dr, &b_tree_genFromB_dr);
   fChain->SetBranchAddress("tree_genFromB_dz", &tree_genFromB_dz, &b_tree_genFromB_dz);
   fChain->SetBranchAddress("tree_genJet_pt", &tree_genJet_pt, &b_tree_genJet_pt);
   fChain->SetBranchAddress("tree_genJet_eta", &tree_genJet_eta, &b_tree_genJet_eta);
   fChain->SetBranchAddress("tree_genJet_phi", &tree_genJet_phi, &b_tree_genJet_phi);
   fChain->SetBranchAddress("tree_genJet_mass", &tree_genJet_mass, &b_tree_genJet_mass);
   fChain->SetBranchAddress("tree_genJet_energy", &tree_genJet_energy, &b_tree_genJet_energy);
   fChain->SetBranchAddress("tree_nLLP", &tree_nLLP, &b_tree_nLLP);
   fChain->SetBranchAddress("tree_LLP", &tree_LLP, &b_tree_LLP);
   fChain->SetBranchAddress("tree_LLP_pt", &tree_LLP_pt, &b_tree_LLP_pt);
   fChain->SetBranchAddress("tree_LLP_eta", &tree_LLP_eta, &b_tree_LLP_eta);
   fChain->SetBranchAddress("tree_LLP_phi", &tree_LLP_phi, &b_tree_LLP_phi);
   fChain->SetBranchAddress("tree_LLP_x", &tree_LLP_x, &b_tree_LLP_x);
   fChain->SetBranchAddress("tree_LLP_y", &tree_LLP_y, &b_tree_LLP_y);
   fChain->SetBranchAddress("tree_LLP_z", &tree_LLP_z, &b_tree_LLP_z);
   fChain->SetBranchAddress("tree_LLP_r", &tree_LLP_r, &b_tree_LLP_r);
   fChain->SetBranchAddress("tree_LLP_dist", &tree_LLP_dist, &b_tree_LLP_dist);
   fChain->SetBranchAddress("tree_LLP_nTrks", &tree_LLP_nTrks, &b_tree_LLP_nTrks);
   fChain->SetBranchAddress("tree_LLP_Vtx_nTrks", &tree_LLP_Vtx_nTrks, &b_tree_LLP_Vtx_nTrks);
   fChain->SetBranchAddress("tree_LLP_Vtx_NChi2", &tree_LLP_Vtx_NChi2, &b_tree_LLP_Vtx_NChi2);
   fChain->SetBranchAddress("tree_LLP_Vtx_dx", &tree_LLP_Vtx_dx, &b_tree_LLP_Vtx_dx);
   fChain->SetBranchAddress("tree_LLP_Vtx_dy", &tree_LLP_Vtx_dy, &b_tree_LLP_Vtx_dy);
   fChain->SetBranchAddress("tree_LLP_Vtx_dz", &tree_LLP_Vtx_dz, &b_tree_LLP_Vtx_dz);
   fChain->SetBranchAddress("tree_LLP_Vtx_dist", &tree_LLP_Vtx_dist, &b_tree_LLP_Vtx_dist);
   fChain->SetBranchAddress("tree_LLP_Vtx_dd", &tree_LLP_Vtx_dd, &b_tree_LLP_Vtx_dd);
   fChain->SetBranchAddress("tree_LLP12_dR", &tree_LLP12_dR, &b_tree_LLP12_dR);
   fChain->SetBranchAddress("tree_LLP12_deta", &tree_LLP12_deta, &b_tree_LLP12_deta);
   fChain->SetBranchAddress("tree_LLP12_dphi", &tree_LLP12_dphi, &b_tree_LLP12_dphi);
   fChain->SetBranchAddress("tree_LLP_Mass", &tree_LLP_Mass, &b_tree_LLP_Mass);
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
   fChain->SetBranchAddress("tree_Hemi_LooseBTag_axes", &tree_Hemi_LooseBTag_axes, &b_tree_Hemi_LooseBTag_axes);
   fChain->SetBranchAddress("tree_Hemi_MediumBTag_axes", &tree_Hemi_MediumBTag_axes, &b_tree_Hemi_MediumBTag_axes);
   fChain->SetBranchAddress("tree_Hemi_TightBTag_axes", &tree_Hemi_TightBTag_axes, &b_tree_Hemi_TightBTag_axes);
   fChain->SetBranchAddress("tree_Hemi_dR12", &tree_Hemi_dR12, &b_tree_Hemi_dR12);
   fChain->SetBranchAddress("tree_Hemi_LLP", &tree_Hemi_LLP, &b_tree_Hemi_LLP);
   fChain->SetBranchAddress("tree_Hemi_LLP_pt", &tree_Hemi_LLP_pt, &b_tree_Hemi_LLP_pt);
   fChain->SetBranchAddress("tree_Hemi_LLP_eta", &tree_Hemi_LLP_eta, &b_tree_Hemi_LLP_eta);
   fChain->SetBranchAddress("tree_Hemi_LLP_phi", &tree_Hemi_LLP_phi, &b_tree_Hemi_LLP_phi);
   fChain->SetBranchAddress("tree_Hemi_LLP_dist", &tree_Hemi_LLP_dist, &b_tree_Hemi_LLP_dist);
   fChain->SetBranchAddress("tree_Hemi_LLP_x", &tree_Hemi_LLP_x, &b_tree_Hemi_LLP_x);
   fChain->SetBranchAddress("tree_Hemi_LLP_y", &tree_Hemi_LLP_y, &b_tree_Hemi_LLP_y);
   fChain->SetBranchAddress("tree_Hemi_LLP_z", &tree_Hemi_LLP_z, &b_tree_Hemi_LLP_z);
   fChain->SetBranchAddress("tree_Hemi_LLP_dR", &tree_Hemi_LLP_dR, &b_tree_Hemi_LLP_dR);
   fChain->SetBranchAddress("tree_Hemi_LLP_mother", &tree_Hemi_LLP_mother, &b_tree_Hemi_LLP_mother);
   fChain->SetBranchAddress("tree_Hemi_LLP_Vtx_dx", &tree_Hemi_LLP_Vtx_dx, &b_tree_Hemi_LLP_Vtx_dx);
   fChain->SetBranchAddress("tree_Hemi_LLP_Vtx_dy", &tree_Hemi_LLP_Vtx_dy, &b_tree_Hemi_LLP_Vtx_dy);
   fChain->SetBranchAddress("tree_Hemi_LLP_Vtx_dz", &tree_Hemi_LLP_Vtx_dz, &b_tree_Hemi_LLP_Vtx_dz);
   fChain->SetBranchAddress("tree_Hemi_LLP_Vtx_dr", &tree_Hemi_LLP_Vtx_dr, &b_tree_Hemi_LLP_Vtx_dr);
   fChain->SetBranchAddress("tree_Hemi_LLP_muOK_dR", &tree_Hemi_LLP_muOK_dR, &b_tree_Hemi_LLP_muOK_dR);
   fChain->SetBranchAddress("tree_Hemi_LLP_muOK_pt", &tree_Hemi_LLP_muOK_pt, &b_tree_Hemi_LLP_muOK_pt);
   fChain->SetBranchAddress("tree_Hemi_LLP_muOK_mass", &tree_Hemi_LLP_muOK_mass, &b_tree_Hemi_LLP_muOK_mass);
   fChain->SetBranchAddress("tree_Hemi_LLP_muNO_dR", &tree_Hemi_LLP_muNO_dR, &b_tree_Hemi_LLP_muNO_dR);
   fChain->SetBranchAddress("tree_Hemi_LLP_muNO_pt", &tree_Hemi_LLP_muNO_pt, &b_tree_Hemi_LLP_muNO_pt);
   fChain->SetBranchAddress("tree_Hemi_LLP_muNO_mass", &tree_Hemi_LLP_muNO_mass, &b_tree_Hemi_LLP_muNO_mass);
   fChain->SetBranchAddress("tree_Hemi_LLP_dR12", &tree_Hemi_LLP_dR12, &b_tree_Hemi_LLP_dR12);
   fChain->SetBranchAddress("tree_Hemi_LLP_ping", &tree_Hemi_LLP_ping, &b_tree_Hemi_LLP_ping);
   fChain->SetBranchAddress("tree_event_LLP_ping", &tree_event_LLP_ping, &b_tree_event_LLP_ping);
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
   fChain->SetBranchAddress("tree_Hemi_Vtx_xError", &tree_Hemi_Vtx_xError, &b_tree_Hemi_Vtx_xError);
   fChain->SetBranchAddress("tree_Hemi_Vtx_yError", &tree_Hemi_Vtx_yError, &b_tree_Hemi_Vtx_yError);
   fChain->SetBranchAddress("tree_Hemi_Vtx_zError", &tree_Hemi_Vtx_zError, &b_tree_Hemi_Vtx_zError);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BTag", &tree_Hemi_Vtx_BTag, &b_tree_Hemi_Vtx_BTag);
   fChain->SetBranchAddress("tree_Hemi_Vtx_trackWeight", &tree_Hemi_Vtx_trackWeight, &b_tree_Hemi_Vtx_trackWeight);
   fChain->SetBranchAddress("tree_Hemi_Vtx_SumtrackWeight", &tree_Hemi_Vtx_SumtrackWeight, &b_tree_Hemi_Vtx_SumtrackWeight);
   fChain->SetBranchAddress("tree_Hemi_Vtx_track_MeanDCA_d", &tree_Hemi_Vtx_track_MeanDCA_d, &b_tree_Hemi_Vtx_track_MeanDCA_d);
   fChain->SetBranchAddress("tree_Hemi_Vtx_Mass", &tree_Hemi_Vtx_Mass, &b_tree_Hemi_Vtx_Mass);
   fChain->SetBranchAddress("tree_Hemi_Vtx_dist", &tree_Hemi_Vtx_dist, &b_tree_Hemi_Vtx_dist);
   fChain->SetBranchAddress("tree_Hemi_Vtx_ntrk10", &tree_Hemi_Vtx_ntrk10, &b_tree_Hemi_Vtx_ntrk10);
   fChain->SetBranchAddress("tree_Hemi_Vtx_ntrk20", &tree_Hemi_Vtx_ntrk20, &b_tree_Hemi_Vtx_ntrk20);
   fChain->SetBranchAddress("tree_event_nVtx", &tree_event_nVtx, &b_tree_event_nVtx);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_dr", &tree_event_Vtx_Vtx_dr, &b_tree_event_Vtx_Vtx_dr);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_dz", &tree_event_Vtx_Vtx_dz, &b_tree_event_Vtx_Vtx_dz);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_dd", &tree_event_Vtx_Vtx_dd, &b_tree_event_Vtx_Vtx_dd);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_reldd", &tree_event_Vtx_Vtx_reldd, &b_tree_event_Vtx_Vtx_reldd);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_dR", &tree_event_Vtx_Vtx_dR, &b_tree_event_Vtx_Vtx_dR);
   fChain->SetBranchAddress("tree_event_Vtx_Vtx_step", &tree_event_Vtx_Vtx_step, &b_tree_event_Vtx_Vtx_step);
   fChain->SetBranchAddress("tree_Hemi_SecLLP", &tree_Hemi_SecLLP, &b_tree_Hemi_SecLLP);
   fChain->SetBranchAddress("tree_Hemi_LLP_SecVtx_dx", &tree_Hemi_LLP_SecVtx_dx, &b_tree_Hemi_LLP_SecVtx_dx);
   fChain->SetBranchAddress("tree_Hemi_LLP_SecVtx_dy", &tree_Hemi_LLP_SecVtx_dy, &b_tree_Hemi_LLP_SecVtx_dy);
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
   fChain->SetBranchAddress("tree_Hemi_SecVtx_trackWeight", &tree_Hemi_SecVtx_trackWeight, &b_tree_Hemi_SecVtx_trackWeight);
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
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v);
fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v", &HLT_Ele27_WPTight_Gsf_v, &b_HLT_Ele27_WPTight_Gsf_v);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_v", &HLT_Ele32_WPTight_Gsf_v, &b_HLT_Ele32_WPTight_Gsf_v);
   fChain->SetBranchAddress("HLT_Ele35_WPTight_Gsf_v", &HLT_Ele35_WPTight_Gsf_v, &b_HLT_Ele35_WPTight_Gsf_v);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v", &HLT_PFMET120_PFMHT120_IDTight_v, &b_HLT_PFMET120_PFMHT120_IDTight_v);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_v, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
   fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned_v", &HLT_PFMET250_HBHECleaned_v, &b_HLT_PFMET250_HBHECleaned_v);
   fChain->SetBranchAddress("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v", &HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v, &b_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v,&b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);    
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v,&b_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v);    
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v,&b_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v);  
   fChain->SetBranchAddress("HLT_IsoMu24_v", &HLT_IsoMu24_v, &b_HLT_IsoMu24_v);
   fChain->SetBranchAddress("HLT_IsoMu27_v", &HLT_IsoMu27_v, &b_HLT_IsoMu27_v);
   fChain->SetBranchAddress("HLT_IsoTkMu24_v", &HLT_IsoTkMu24_v, &b_HLT_IsoTkMu24_v);
   Notify();
}

Bool_t MiniDATAMCNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MiniDATAMCNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MiniDATAMCNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

//$$
double DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double DeltaPhi = TMath::Abs(phi2 - phi1);
  if (DeltaPhi > 3.141593 ) DeltaPhi = 2.*3.141593 - DeltaPhi;
  return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}
double DeltaPhi(double phi1, double phi2) {
  double DeltaPhi = phi1 - phi2;
  if (abs(DeltaPhi) > 3.141593 ) {
    DeltaPhi = 2.*3.141593 - abs(DeltaPhi);
    DeltaPhi = -DeltaPhi * (phi1 - phi2) / abs(phi1 - phi2);
  }
  return DeltaPhi;
}
#endif // #ifdef MiniDATAMCNtuple_cxx
