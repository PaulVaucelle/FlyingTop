//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 31 09:09:59 2023 by ROOT version 6.14/09
// from TTree ttree/ttree
// found on file: Ntuple_50_test.root
//////////////////////////////////////////////////////////

#ifndef TreeReader_h
#define TreeReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "ROOT/RVec.hxx"

class TreeReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runNumber;
   Int_t           eventNumber;
   Int_t           lumiBlock;
   Int_t           tree_nPV;
   Float_t         tree_PV_x;
   Float_t         tree_PV_y;
   Float_t         tree_PV_z;
   Float_t         tree_PV_ez;
   Float_t         tree_PV_NChi2;
   Int_t           tree_PV_ndf;
   vector<float>   *tree_Evts_MVAval;
   vector<int>     *tree_allPV_i;
   vector<float>   *tree_allPV_x;
   vector<float>   *tree_allPV_y;
   vector<float>   *tree_allPV_z;
   vector<float>   *tree_allPV_ex;
   vector<float>   *tree_allPV_ey;
   vector<float>   *tree_allPV_ez;
   vector<float>   *tree_allPV_NChi2;
   vector<int>     *tree_allPV_ndf;
   Float_t         tree_bs_PosX;
   Float_t         tree_bs_PosY;
   Float_t         tree_bs_PosZ;
   Int_t           tree_NbrOfZCand;
   Bool_t          tree_Filter;
   Int_t           tree_nK0;
   vector<float>   *tree_K0_x;
   vector<float>   *tree_K0_y;
   vector<float>   *tree_K0_z;
   vector<float>   *tree_K0_r;
   vector<float>   *tree_K0_NChi2;
   vector<float>   *tree_K0_ndf;
   vector<float>   *tree_K0_mass;
   vector<float>   *tree_K0_pt;
   vector<float>   *tree_K0_eta;
   vector<float>   *tree_K0_phi;
   vector<unsigned int> *tree_K0_nDaughters;
   Int_t           tree_nLambda;
   vector<float>   *tree_L0_x;
   vector<float>   *tree_L0_y;
   vector<float>   *tree_L0_z;
   vector<float>   *tree_L0_r;
   vector<float>   *tree_L0_NChi2;
   vector<float>   *tree_L0_ndf;
   vector<unsigned int> *tree_L0_nDaughters;
   vector<float>   *tree_L0_mass;
   vector<float>   *tree_L0_pt;
   vector<float>   *tree_L0_eta;
   vector<float>   *tree_L0_phi;
   Int_t           tree_nV0_reco;
   vector<float>   *tree_V0_reco_x;
   vector<float>   *tree_V0_reco_y;
   vector<float>   *tree_V0_reco_z;
   vector<float>   *tree_V0_reco_r;
   vector<float>   *tree_V0_reco_drSig;
   vector<float>   *tree_V0_reco_dzSig;
   vector<float>   *tree_V0_reco_angleXY;
   vector<float>   *tree_V0_reco_angleZ;
   vector<float>   *tree_V0_reco_NChi2;
   vector<float>   *tree_V0_reco_ndf;
   vector<float>   *tree_V0_reco_mass;
   vector<float>   *tree_V0_reco_pt;
   vector<float>   *tree_V0_reco_eta;
   vector<float>   *tree_V0_reco_phi;
   vector<int>     *tree_V0_reco_source;
   vector<bool>    *tree_V0_reco_badTkHit;
   vector<float>   *tree_V0_reco_dca;
   Int_t           tree_nSecInt;
   vector<float>   *tree_SecInt_x;
   vector<float>   *tree_SecInt_y;
   vector<float>   *tree_SecInt_z;
   vector<float>   *tree_SecInt_r;
   vector<float>   *tree_SecInt_d;
   vector<float>   *tree_SecInt_drSig;
   vector<float>   *tree_SecInt_dzSig;
   vector<float>   *tree_SecInt_angleXY;
   vector<float>   *tree_SecInt_angleZ;
   vector<float>   *tree_SecInt_NChi2;
   vector<float>   *tree_SecInt_ndf;
   vector<float>   *tree_SecInt_mass;
   vector<float>   *tree_SecInt_pt;
   vector<float>   *tree_SecInt_eta;
   vector<float>   *tree_SecInt_phi;
   vector<int>     *tree_SecInt_charge;
   vector<bool>    *tree_SecInt_badTkHit;
   vector<float>   *tree_SecInt_dca;
   vector<bool>    *tree_SecInt_selec;
   vector<int>     *tree_SecInt_layer;
   vector<int>     *tree_SecInt_LLP;
   vector<float>   *tree_SecInt_LLP_dr;
   vector<float>   *tree_SecInt_LLP_dz;
   vector<float>   *tree_SecInt_LLP_dd;
   Int_t           tree_nYConv;
   vector<float>   *tree_Yc_x;
   vector<float>   *tree_Yc_y;
   vector<float>   *tree_Yc_z;
   vector<float>   *tree_Yc_r;
   vector<float>   *tree_Yc_dr0;
   vector<float>   *tree_Yc_dr1;
   vector<float>   *tree_Yc_dz0;
   vector<float>   *tree_Yc_dz1;
   vector<float>   *tree_Yc_costheta;
   vector<int>     *tree_Yc_layer;
   vector<float>   *tree_Yc_NChi2;
   vector<float>   *tree_Yc_ndf;
   vector<unsigned int> *tree_Yc_nDaughters;
   vector<float>   *tree_Yc_pt;
   vector<float>   *tree_Yc_eta;
   vector<float>   *tree_Yc_phi;
   vector<float>   *tree_Yc_mass;
   Int_t           tree_Yc_ntracks;
   vector<int>     *tree_Yc_tracks_index;
   vector<int>     *tree_Yc_tracks_charge;
   vector<float>   *tree_Yc_tracks_pt;
   vector<float>   *tree_Yc_tracks_eta;
   vector<float>   *tree_Yc_tracks_phi;
   vector<float>   *tree_Yc_tracks_phi0;
   Float_t         tree_PFMet_et;
   Float_t         tree_PFMet_phi;
   Float_t         tree_PFMet_sig;
   Int_t           tree_njet;
   vector<float>   *tree_jet_E;
   vector<float>   *tree_jet_pt;
   vector<float>   *tree_jet_eta;
   vector<float>   *tree_jet_phi;
   vector<float>   *tree_jet_btag_DeepCSV;
   vector<float>   *tree_jet_btag_DeepJet;
   vector<float>   *tree_jet_leadingpt;
   vector<float>   *tree_jet_leadingpt2;
   vector<float>   *tree_jet_leadingMuon_dR;
   vector<float>   *tree_jet_leadingMuon2_dR;
   vector<float>   *tree_jet_jet_dR;
   vector<float>   *tree_jet_jet_dPhi;
   vector<float>   *tree_jet_jet_dEta;
   vector<float>   *tree_muon_jet_dRmin;
   vector<float>   *tree_muon_jet_dRmax;
   Float_t         tree_HT;
   Int_t           tree_electron_nEle;
   vector<float>   *tree_electron_pt;
   vector<float>   *tree_electron_eta;
   vector<float>   *tree_electron_phi;
   vector<float>   *tree_electron_x;
   vector<float>   *tree_electron_y;
   vector<float>   *tree_electron_z;
   vector<float>   *tree_electron_px;
   vector<float>   *tree_electron_py;
   vector<float>   *tree_electron_pz;
   vector<float>   *tree_electron_energy;
   vector<int>     *tree_electron_charge;
   vector<float>   *tree_electron_isoR4;
   vector<bool>    *tree_electron_IsLoose;
   vector<bool>    *tree_electron_IsMedium;
   vector<bool>    *tree_electron_IsTight;
   vector<bool>    *tree_electron_trigger_Ele;
   vector<bool>    *tree_electron_trigger_diEle;
   vector<float>   *tree_electron_dxy;
   vector<float>   *tree_electron_dz;
   vector<float>   *tree_ST;
   Float_t         tree_Mmumu;
   vector<float>   *tree_muon_pt;
   vector<float>   *tree_muon_eta;
   vector<float>   *tree_muon_phi;
   vector<float>   *tree_muon_x;
   vector<float>   *tree_muon_y;
   vector<float>   *tree_muon_z;
   vector<float>   *tree_muon_px;
   vector<float>   *tree_muon_py;
   vector<float>   *tree_muon_pz;
   vector<float>   *tree_muon_energy;
   vector<float>   *tree_muon_dxy;
   vector<float>   *tree_muon_dxyError;
   vector<float>   *tree_muon_dz;
   vector<float>   *tree_muon_dzError;
   vector<int>     *tree_muon_charge;
   vector<bool>    *tree_muon_isLoose;
   vector<bool>    *tree_muon_isTight;
   vector<bool>    *tree_muon_isGlobal;
   vector<float>   *tree_muon_isoR3;
   vector<bool>    *tree_muon_trigger_dimu;
   vector<bool>    *tree_muon_trigger_isomu;
   vector<bool>    *tree_muon_PFIsoLoose;
   vector<bool>    *tree_muon_PFIsoMedium;
   vector<bool>    *tree_muon_PFIsoTight;
   vector<float>   *tree_muon_nmu;
   vector<float>   *tree_muon_leadingpt;
   vector<float>   *tree_muon_leadingpt2;
   vector<float>   *tree_muon_muon_dR;
   vector<float>   *tree_muon_muon_dPhi;
   vector<float>   *tree_muon_muon_dEta;
   Int_t           tree_TRACK_SIZE;
   Int_t           tree_nTracks;
   Int_t           tree_nLostTracks;
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
   vector<float>   *tree_track_dzTOpu;
   vector<float>   *tree_track_dzSigTOpu;
   vector<unsigned int> *tree_track_algo;
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
   vector<float>   *tree_track_btag;
   vector<float>   *tree_track_energy;
   vector<int>     *tree_track_Hemi;
   vector<float>   *tree_track_Hemi_dR;
   vector<float>   *tree_track_Hemi_dRmax;
   vector<float>   *tree_track_Hemi_mva_NChi2;
   vector<bool>    *tree_track_Hemi_ping;
   vector<float>   *tree_track_Hemi_dFirstVtx;
   vector<int>     *tree_track_Hemi_LLP;
   vector<int>     *tree_track_sim_LLP;
   vector<bool>    *tree_track_sim_isFromB;
   vector<bool>    *tree_track_sim_isFromC;
   vector<float>   *tree_track_sim_pt;
   vector<float>   *tree_track_sim_eta;
   vector<float>   *tree_track_sim_phi;
   vector<int>     *tree_track_sim_charge;
   vector<int>     *tree_track_sim_pdgId;
   vector<float>   *tree_track_sim_mass;
   vector<float>   *tree_track_sim_x;
   vector<float>   *tree_track_sim_y;
   vector<float>   *tree_track_sim_z;
   vector<float>   *tree_track_sim_dFirstGen;
   vector<float>   *tree_track_sim_LLP_r;
   vector<float>   *tree_track_sim_LLP_z;
   Float_t         tree_GenPVx;
   Float_t         tree_GenPVy;
   Float_t         tree_GenPVz;
   Int_t           tree_smu_mass;
   Int_t           tree_neu_mass;
   vector<float>   *tree_genParticle_pt;
   vector<float>   *tree_genParticle_eta;
   vector<float>   *tree_genParticle_phi;
   vector<float>   *tree_genParticle_charge;
   vector<int>     *tree_genParticle_pdgId;
   vector<float>   *tree_genParticle_mass;
   vector<float>   *tree_genParticle_x;
   vector<float>   *tree_genParticle_y;
   vector<float>   *tree_genParticle_z;
   vector<float>   *tree_genParticle_px;
   vector<float>   *tree_genParticle_py;
   vector<float>   *tree_genParticle_pz;
   vector<float>   *tree_genParticle_energy;
   vector<bool>    *tree_genParticle_isPromptFinalState;
   vector<int>     *tree_genParticle_statusCode;
   vector<int>     *tree_genParticle_mother_pdgId;
   vector<int>     *tree_genParticle_LLP;
   vector<float>   *tree_genParticle_ct0;
   Int_t           tree_ngenPackPart;
   vector<float>   *tree_genPackPart_pt;
   vector<float>   *tree_genPackPart_eta;
   vector<float>   *tree_genPackPart_phi;
   vector<float>   *tree_genPackPart_charge;
   vector<int>     *tree_genPackPart_pdgId;
   vector<float>   *tree_genPackPart_mass;
   vector<float>   *tree_genPackPart_x;
   vector<float>   *tree_genPackPart_y;
   vector<float>   *tree_genPackPart_z;
   vector<int>     *tree_genPackPart_mother_pdgId;
   vector<bool>    *tree_genPackPart_isFromB;
   vector<bool>    *tree_genPackPart_isFromC;
   Int_t           tree_ngenFromLLP;
   vector<int>     *tree_genFromLLP_LLP;
   vector<float>   *tree_genFromLLP_pt;
   vector<float>   *tree_genFromLLP_eta;
   vector<float>   *tree_genFromLLP_phi;
   vector<float>   *tree_genFromLLP_charge;
   vector<int>     *tree_genFromLLP_pdgId;
   vector<float>   *tree_genFromLLP_mass;
   vector<float>   *tree_genFromLLP_x;
   vector<float>   *tree_genFromLLP_y;
   vector<float>   *tree_genFromLLP_z;
   vector<int>     *tree_genFromLLP_mother_pdgId;
   vector<bool>    *tree_genFromLLP_isFromB;
   vector<bool>    *tree_genFromLLP_isFromC;
   vector<float>   *tree_genAxis_dRneuneu;
   vector<float>   *tree_genAxis_dPhineuneu;
   vector<float>   *tree_genAxis_dEtaneuneu;
   vector<float>   *tree_GenAxes_Mass;
   vector<float>   *tree_GenAxis_Neu_dRmin;
   vector<float>   *tree_GenAxis_Neu_dRmax;
   vector<float>   *tree_GenAxis_RecoAxis_dRmin;
   vector<float>   *tree_GenAxis_RecoAxis_dRmax;
   Int_t           tree_nFromC;
   vector<float>   *tree_genFromC_pt;
   vector<float>   *tree_genFromC_eta;
   vector<float>   *tree_genFromC_phi;
   vector<float>   *tree_genFromC_charge;
   vector<int>     *tree_genFromC_pdgId;
   vector<float>   *tree_genFromC_x;
   vector<float>   *tree_genFromC_y;
   vector<float>   *tree_genFromC_z;
   vector<int>     *tree_genFromC_mother_pdgId;
   Int_t           tree_nFromB;
   vector<float>   *tree_genFromB_pt;
   vector<float>   *tree_genFromB_eta;
   vector<float>   *tree_genFromB_phi;
   vector<float>   *tree_genFromB_charge;
   vector<int>     *tree_genFromB_pdgId;
   vector<float>   *tree_genFromB_x;
   vector<float>   *tree_genFromB_y;
   vector<float>   *tree_genFromB_z;
   vector<int>     *tree_genFromB_mother_pdgId;
   vector<float>   *tree_genFromB_dd;
   vector<float>   *tree_genFromB_dr;
   vector<float>   *tree_genFromB_dz;
   vector<float>   *tree_genJet_pt;
   vector<float>   *tree_genJet_eta;
   vector<float>   *tree_genJet_phi;
   vector<float>   *tree_genJet_mass;
   vector<float>   *tree_genJet_energy;
   Int_t           tree_nLLP;
   vector<int>     *tree_LLP;
   vector<float>   *tree_LLP_pt;
   vector<float>   *tree_LLP_eta;
   vector<float>   *tree_LLP_phi;
   vector<float>   *tree_LLP_x;
   vector<float>   *tree_LLP_y;
   vector<float>   *tree_LLP_z;
   vector<float>   *tree_LLP_r;
   vector<float>   *tree_LLP_dist;
   vector<int>     *tree_LLP_nTrks;
   vector<int>     *tree_LLP_Vtx_nTrks;
   vector<float>   *tree_LLP_Vtx_NChi2;
   vector<float>   *tree_LLP_Vtx_dx;
   vector<float>   *tree_LLP_Vtx_dy;
   vector<float>   *tree_LLP_Vtx_dz;
   vector<float>   *tree_LLP_Vtx_dist;
   vector<float>   *tree_LLP_Vtx_dd;
   vector<float>   *tree_LLP12_dR;
   vector<float>   *tree_LLP12_deta;
   vector<float>   *tree_LLP12_dphi;
   vector<float>   *tree_LLP_Mass;
   vector<int>     *tree_Hemi;
   vector<float>   *tree_Hemi_Mass;
   vector<int>     *tree_Hemi_njet;
   vector<float>   *tree_Hemi_eta;
   vector<float>   *tree_Hemi_phi;
   vector<float>   *tree_Hemi_dR;
   vector<int>     *tree_Hemi_nTrks;
   vector<int>     *tree_Hemi_nTrks_sig;
   vector<int>     *tree_Hemi_nTrks_bad;
   vector<int>     *tree_Hemi_LLP;
   vector<float>   *tree_Hemi_LLP_pt;
   vector<float>   *tree_Hemi_LLP_eta;
   vector<float>   *tree_Hemi_LLP_phi;
   vector<float>   *tree_Hemi_LLP_dist;
   vector<float>   *tree_Hemi_LLP_x;
   vector<float>   *tree_Hemi_LLP_y;
   vector<float>   *tree_Hemi_LLP_z;
   vector<int>     *tree_Hemi_Vtx_step;
   vector<float>   *tree_Hemi_Vtx_NChi2;
   vector<int>     *tree_Hemi_Vtx_nTrks;
   vector<int>     *tree_Hemi_Vtx_nTrks_sig;
   vector<int>     *tree_Hemi_Vtx_nTrks_bad;
   vector<float>   *tree_Hemi_Vtx_x;
   vector<float>   *tree_Hemi_Vtx_y;
   vector<float>   *tree_Hemi_Vtx_r;
   vector<float>   *tree_Hemi_Vtx_z;
   vector<float>   *tree_Hemi_Vtx_xError;
   vector<float>   *tree_Hemi_Vtx_yError;
   vector<float>   *tree_Hemi_Vtx_zError;
   vector<float>   *tree_Hemi_Vtx_eta;
   vector<float>   *tree_Hemi_Vtx_Vtx_dr;
   vector<float>   *tree_Hemi_Vtx_Vtx_dz;
   vector<float>   *tree_Hemi_Vtx_Vtx_dd;
   vector<float>   *tree_Hemi_Vtx_BTag;
   vector<int>     *tree_Hemi_Vtx_nVtx;
   vector<float>   *tree_Hemi_Vtx_trackWeight;
   vector<float>   *tree_Hemi_Vtx_MeantrackWeight;
   vector<float>   *tree_Hemi_Vtx_track_DCA_x;
   vector<float>   *tree_Hemi_Vtx_track_DCA_y;
   vector<float>   *tree_Hemi_Vtx_track_DCA_z;
   vector<float>   *tree_Hemi_Vtx_track_DCA_r;
   vector<float>   *tree_Hemi_Vtx_track_DCA_d;
   vector<float>   *tree_Hemi_Vtx_track_MeanDCA_d;
   vector<float>   *tree_Hemi_Vtx_Mass;
   vector<float>   *tree_Hemi_Vtx_MVAval;
   vector<float>   *tree_Hemi_Vtx_MVAval_Step1;
   vector<float>   *tree_Hemi_Vtx_TVtx_dx;
   vector<float>   *tree_Hemi_Vtx_TVtx_dy;
   vector<float>   *tree_Hemi_Vtx_TVtx_dz;
   vector<float>   *tree_Hemi_Vtx_TVtx_NChi2;
   vector<float>   *tree_Hemi_Vtx_dist;
   vector<float>   *tree_Hemi_Vtx_dx;
   vector<float>   *tree_Hemi_Vtx_dy;
   vector<float>   *tree_Hemi_Vtx_dz;
   vector<float>   *tree_Hemi_Vtx_dr;
   vector<float>   *tree_Hemi_Vtx_dd;
   vector<float>   *tree_Hemi_dR12;
   vector<float>   *tree_Hemi_LLP_dR12;
   vector<float>   *tree_Hemi_Vtx_ddbad;
   vector<int>     *tree_Hemi_Vtx_ntrk10;
   vector<int>     *tree_Hemi_Vtx_ntrk20;
   vector<int>     *tree_track_Hemi_isjet;
   vector<float>   *tree_Hemi_Vtx_ddToBkg;
   vector<bool>    *tree_Hemi_LLP_ping;
   vector<int>     *tree_event_LLP_ping;
   vector<int>     *tree_Hemi_LooseBTag_axes;
   vector<int>     *tree_Hemi_MediumBTag_axes;
   vector<int>     *tree_Hemi_TightBTag_axes;
   Bool_t          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Bool_t          HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;
   Bool_t          HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   Bool_t          HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Bool_t          HLT_Ele27_WPTight_Gsf_v;
   Bool_t          HLT_Ele32_WPTight_Gsf_v;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   Bool_t          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_v;
   Bool_t          HLT_PFMET120_PFMHT120_IDTight_PFHT60_v;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;
   Bool_t          HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;
   Bool_t          HLT_PFMET250_HBHECleaned_v;
   Bool_t          HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;
   Bool_t          HLT_IsoMu24_v;
   Bool_t          HLT_IsoMu27_v;
   Bool_t          HLT_PFHT180_v;
   Bool_t          HLT_PFHT250_v;
   Bool_t          HLT_PFHT370_v;
   Bool_t          HLT_PFHT430_v;
   Bool_t          HLT_PFHT510_v;
   Bool_t          HLT_PFHT590_v;
   Bool_t          HLT_PFHT680_v;
   Bool_t          HLT_PFHT500_PFMET100_PFMHT100_IDTight_v;
   Bool_t          HLT_PFHT700_PFMET85_PFMHT85_IDTight_v;
   Bool_t          HLT_PFHT800_PFMET75_PFMHT75_IDTight_v;

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_lumiBlock;   //!
   TBranch        *b_tree_nPV;   //!
   TBranch        *b_tree_PV_x;   //!
   TBranch        *b_tree_PV_y;   //!
   TBranch        *b_tree_PV_z;   //!
   TBranch        *b_tree_PV_ez;   //!
   TBranch        *b_tree_PV_NChi2;   //!
   TBranch        *b_tree_PV_ndf;   //!
   TBranch        *b_tree_Evts_MVAval;
   TBranch        *b_tree_allPV_i;   //!
   TBranch        *b_tree_allPV_x;   //!
   TBranch        *b_tree_allPV_y;   //!
   TBranch        *b_tree_allPV_z;   //!
   TBranch        *b_tree_allPV_ex;   //!
   TBranch        *b_tree_allPV_ey;   //!
   TBranch        *b_tree_allPV_ez;   //!
   TBranch        *b_tree_allPV_NChi2;   //!
   TBranch        *b_tree_allPV_ndf;   //!
   TBranch        *b_tree_bs_PosX;   //!
   TBranch        *b_tree_bs_PosY;   //!
   TBranch        *b_tree_bs_PosZ;   //!
   TBranch        *b_tree_NbrOfZCand;   //!
   TBranch        *b_tree_Filter;   //!
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
   TBranch        *b_tree_PFMet_et;   //!
   TBranch        *b_tree_PFMet_phi;   //!
   TBranch        *b_tree_PFMet_sig;   //!
   TBranch        *b_tree_njet;   //!
   TBranch        *b_tree_jet_E;   //!
   TBranch        *b_tree_jet_pt;   //!
   TBranch        *b_tree_jet_eta;   //!
   TBranch        *b_tree_jet_phi;   //!
   TBranch        *b_tree_jet_btag_DeepCSV;   //!
   TBranch        *b_tree_jet_btag_DeepJet;   //!
   TBranch        *b_tree_jet_leadingpt;   //!
   TBranch        *b_tree_jet_leadingpt2;   //!
   TBranch        *b_tree_jet_leadingMuon_dR;   //!
   TBranch        *b_tree_jet_leadingMuon2_dR;   //!
   TBranch        *b_tree_jet_jet_dR;   //!
   TBranch        *b_tree_jet_jet_dPhi;   //!
   TBranch        *b_tree_jet_jet_dEta;   //!
   TBranch        *b_tree_muon_jet_dRmin;   //!
   TBranch        *b_tree_muon_jet_dRmax;   //!
   TBranch        *b_tree_HT;   //!
   TBranch        *b_tree_electron_nEle;   //!
   TBranch        *b_tree_electron_pt;   //!
   TBranch        *b_tree_electron_eta;   //!
   TBranch        *b_tree_electron_phi;   //!
   TBranch        *b_tree_electron_x;   //!
   TBranch        *b_tree_electron_y;   //!
   TBranch        *b_tree_electron_z;   //!
   TBranch        *b_tree_electron_px;   //!
   TBranch        *b_tree_electron_py;   //!
   TBranch        *b_tree_electron_pz;   //!
   TBranch        *b_tree_electron_energy;   //!
   TBranch        *b_tree_electron_charge;   //!
   TBranch        *b_tree_electron_isoR4;   //!
   TBranch        *b_tree_electron_IsLoose;   //!
   TBranch        *b_tree_electron_IsMedium;   //!
   TBranch        *b_tree_electron_IsTight;   //!
   TBranch        *b_tree_electron_trigger_Ele;   //!
   TBranch        *b_tree_electron_trigger_diEle;   //!
   TBranch        *b_tree_electron_dxy;   //!
   TBranch        *b_tree_electron_dz;   //!
   TBranch        *b_tree_ST;   //!
   TBranch        *b_tree_Mmumu;   //!
   TBranch        *b_tree_muon_pt;   //!
   TBranch        *b_tree_muon_eta;   //!
   TBranch        *b_tree_muon_phi;   //!
   TBranch        *b_tree_muon_x;   //!
   TBranch        *b_tree_muon_y;   //!
   TBranch        *b_tree_muon_z;   //!
   TBranch        *b_tree_muon_px;   //!
   TBranch        *b_tree_muon_py;   //!
   TBranch        *b_tree_muon_pz;   //!
   TBranch        *b_tree_muon_energy;   //!
   TBranch        *b_tree_muon_dxy;   //!
   TBranch        *b_tree_muon_dxyError;   //!
   TBranch        *b_tree_muon_dz;   //!
   TBranch        *b_tree_muon_dzError;   //!
   TBranch        *b_tree_muon_charge;   //!
   TBranch        *b_tree_muon_isLoose;   //!
   TBranch        *b_tree_muon_isTight;   //!
   TBranch        *b_tree_muon_isGlobal;   //!
   TBranch        *b_tree_muon_isoR3;   //!
   TBranch        *b_tree_muon_trigger_dimu;   //!
   TBranch        *b_tree_muon_trigger_isomu;   //!
   TBranch        *b_tree_muon_PFIsoLoose;   //!
   TBranch        *b_tree_muon_PFIsoMedium;   //!
   TBranch        *b_tree_muon_PFIsoTight;   //!
   TBranch        *b_tree_muon_nmu;   //!
   TBranch        *b_tree_muon_leadingpt;   //!
   TBranch        *b_tree_muon_leadingpt2;   //!
   TBranch        *b_tree_muon_muon_dR;   //!
   TBranch        *b_tree_muon_muon_dPhi;   //!
   TBranch        *b_tree_muon_muon_dEta;   //!
   TBranch        *b_tree_TRACK_SIZE;
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
   TBranch        *b_tree_track_dzTOpu;   //!
   TBranch        *b_tree_track_dzSigTOpu;   //!
   TBranch        *b_tree_track_algo;   //!
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
   TBranch        *b_tree_GenPVx;   //!
   TBranch        *b_tree_GenPVy;   //!
   TBranch        *b_tree_GenPVz;   //!
   TBranch        *b_tree_smu_mass;
   TBranch        *b_tree_neu_mass;
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
   TBranch        *b_tree_genParticle_ct0;
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
   TBranch        *b_tree_genAxis_dRneuneu;   //!
   TBranch        *b_tree_genAxis_dPhineuneu;   //!
   TBranch        *b_tree_genAxis_dEtaneuneu;   //!
   TBranch        *b_tree_GenAxes_Mass;   //!
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
   TBranch        *b_tree_Hemi_Mass;   //!
   TBranch        *b_tree_Hemi_njet;   //!
   TBranch        *b_tree_Hemi_eta;   //!
   TBranch        *b_tree_Hemi_phi;   //!
   TBranch        *b_tree_Hemi_dR;   //!
   TBranch        *b_tree_Hemi_nTrks;   //!
   TBranch        *b_tree_Hemi_nTrks_sig;   //!
   TBranch        *b_tree_Hemi_nTrks_bad;   //!
   TBranch        *b_tree_Hemi_LLP;   //!
   TBranch        *b_tree_Hemi_LLP_pt;   //!
   TBranch        *b_tree_Hemi_LLP_eta;   //!
   TBranch        *b_tree_Hemi_LLP_phi;   //!
   TBranch        *b_tree_Hemi_LLP_dist;   //!
   TBranch        *b_tree_Hemi_LLP_x;   //!
   TBranch        *b_tree_Hemi_LLP_y;   //!
   TBranch        *b_tree_Hemi_LLP_z;   //!
   TBranch        *b_tree_Hemi_Vtx_step;   //!
   TBranch        *b_tree_Hemi_Vtx_NChi2;   //!
   TBranch        *b_tree_Hemi_Vtx_nTrks;   //!
   TBranch        *b_tree_Hemi_Vtx_nTrks_sig;   //!
   TBranch        *b_tree_Hemi_Vtx_nTrks_bad;   //!
   TBranch        *b_tree_Hemi_Vtx_x;   //!
   TBranch        *b_tree_Hemi_Vtx_y;   //!
   TBranch        *b_tree_Hemi_Vtx_r;   //!
   TBranch        *b_tree_Hemi_Vtx_z;   //!
   TBranch        *b_tree_Hemi_Vtx_xError;   //!
   TBranch        *b_tree_Hemi_Vtx_yError;   //!
   TBranch        *b_tree_Hemi_Vtx_zError;   //!
   TBranch        *b_tree_Hemi_Vtx_eta;   //!
   TBranch        *b_tree_Hemi_Vtx_Vtx_dr;   //!
   TBranch        *b_tree_Hemi_Vtx_Vtx_dz;   //!
   TBranch        *b_tree_Hemi_Vtx_Vtx_dd;   //!
   TBranch        *b_tree_Hemi_Vtx_BTag;   //!
   TBranch        *b_tree_Hemi_Vtx_nVtx;   //!
   TBranch        *b_tree_Hemi_Vtx_trackWeight;   //!
   TBranch        *b_tree_Hemi_Vtx_MeantrackWeight;   //!
   TBranch        *b_tree_Hemi_Vtx_track_DCA_x;   //!
   TBranch        *b_tree_Hemi_Vtx_track_DCA_y;   //!
   TBranch        *b_tree_Hemi_Vtx_track_DCA_z;   //!
   TBranch        *b_tree_Hemi_Vtx_track_DCA_r;   //!
   TBranch        *b_tree_Hemi_Vtx_track_DCA_d;   //!
   TBranch        *b_tree_Hemi_Vtx_track_MeanDCA_d;   //!
   TBranch        *b_tree_Hemi_Vtx_Mass;   //!
   TBranch        *b_tree_Hemi_Vtx_MVAval;   //!
   TBranch        *b_tree_Hemi_Vtx_MVAval_Step1;   //!
   TBranch        *b_tree_Hemi_Vtx_TVtx_dx;   //!
   TBranch        *b_tree_Hemi_Vtx_TVtx_dy;   //!
   TBranch        *b_tree_Hemi_Vtx_TVtx_dz;   //!
   TBranch        *b_tree_Hemi_Vtx_TVtx_NChi2;   //!
   TBranch        *b_tree_Hemi_Vtx_dist;   //!
   TBranch        *b_tree_Hemi_Vtx_dx;   //!
   TBranch        *b_tree_Hemi_Vtx_dy;   //!
   TBranch        *b_tree_Hemi_Vtx_dz;   //!
   TBranch        *b_tree_Hemi_Vtx_dr;   //!
   TBranch        *b_tree_Hemi_Vtx_dd;   //!
   TBranch        *b_tree_Hemi_dR12;   //!
   TBranch        *b_tree_Hemi_LLP_dR12;   //!
   TBranch        *b_tree_Hemi_Vtx_ddbad;   //!
   TBranch        *b_tree_Hemi_Vtx_ntrk10;   //!
   TBranch        *b_tree_Hemi_Vtx_ntrk20;   //!
   TBranch        *b_tree_track_Hemi_isjet;   //!
   TBranch        *b_tree_Hemi_Vtx_ddToBkg;   //!
   TBranch        *b_tree_Hemi_LLP_ping;   //!
   TBranch        *b_tree_event_LLP_ping;   //!
   TBranch        *b_tree_Hemi_LooseBTag_axes;   //!
   TBranch        *b_tree_Hemi_MediumBTag_axes;   //!
   TBranch        *b_tree_Hemi_TightBTag_axes;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_Ele27_WPTight_Gsf_v;   //!
   TBranch        *b_HLT_Ele32_WPTight_Gsf_v;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_v;   //!
   TBranch        *b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;   //!
   TBranch        *b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;   //!
   TBranch        *b_HLT_PFMET250_HBHECleaned_v;   //!
   TBranch        *b_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;   //!
   TBranch        *b_HLT_IsoMu24_v;   //!
   TBranch        *b_HLT_IsoMu27_v;   //!
   TBranch        *b_HLT_PFHT180_v;   //!
   TBranch        *b_HLT_PFHT250_v;   //!
   TBranch        *b_HLT_PFHT370_v;   //!
   TBranch        *b_HLT_PFHT430_v;   //!
   TBranch        *b_HLT_PFHT510_v;   //!
   TBranch        *b_HLT_PFHT590_v;   //!
   TBranch        *b_HLT_PFHT680_v;   //!
   TBranch        *b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v;   //!
   TBranch        *b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v;   //!
   TBranch        *b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v;   //!

   // TreeReader(TTree *tree=0);
      TreeReader(TTree *tree=0, TString sample="", std::vector<TString> thesystlist = std::vector<TString>());
   virtual ~TreeReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   // virtual void     Loop();
   virtual void     Loop(TString sample, bool Signal);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

      //add functions to create and fill histograms
 
   void initializeHisto(TString sample, bool isfirstset);
   void addHisto( TString var, TString selstep, TString sample, int nbins, float min, float max);
   void fillHisto(TString var, TString selstep, TString sample, float val, float weight);


   std::vector<TH1F*> histo_list_;
   std::map<std::string, int> histo_map_;


   int numb_histo;
   void deleteHisto();


   std::vector<TString> systlist;

};

#endif

#ifdef TreeReader_cxx
TreeReader::TreeReader(TTree *tree, TString sample, std::vector<TString> thesystlist) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(   ("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/Ntuple_"+sample+".root").Data()   );
      if (!f || !f->IsOpen()) {
         f = new TFile(  ("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/Ntuple_"+sample+".root").Data() );
      }
      //f->GetObject("events",tree);
      TDirectory * dir = (TDirectory*)f->Get("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/Ntuple_"+sample+".root:/FlyingTop");
      dir->GetObject("ttree",tree);
      systlist = thesystlist;
   }
   Init(tree);
}

// TreeReader::TreeReader(TTree *tree) : fChain(0) 
// {
// // if parameter tree is not specified (or zero), connect the file
// // used to generate this class and read the Tree.
//    if (tree == 0) {
//       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Ntuple_50_test.root");
//       if (!f || !f->IsOpen()) {
//          f = new TFile("Ntuple_50_test.root");
//       }
//       TDirectory * dir = (TDirectory*)f->Get("Ntuple_50_test.root:/FlyingTop");
//       dir->GetObject("ttree",tree);

//    }
//    Init(tree);
// }

TreeReader::~TreeReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeReader::LoadTree(Long64_t entry)
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

void TreeReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   tree_allPV_i = 0;
   tree_allPV_x = 0;
   tree_allPV_y = 0;
   tree_allPV_z = 0;
   tree_allPV_ex = 0;
   tree_allPV_ey = 0;
   tree_allPV_ez = 0;
   tree_allPV_NChi2 = 0;
   tree_allPV_ndf = 0;
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
   tree_jet_E = 0;
   tree_jet_pt = 0;
   tree_jet_eta = 0;
   tree_jet_phi = 0;
   tree_jet_btag_DeepCSV = 0;
   tree_jet_btag_DeepJet = 0;
   tree_jet_leadingpt = 0;
   tree_jet_leadingpt2 = 0;
   tree_jet_leadingMuon_dR = 0;
   tree_jet_leadingMuon2_dR = 0;
   tree_jet_jet_dR = 0;
   tree_jet_jet_dPhi = 0;
   tree_jet_jet_dEta = 0;
   tree_muon_jet_dRmin = 0;
   tree_muon_jet_dRmax = 0;
   tree_electron_pt = 0;
   tree_electron_eta = 0;
   tree_electron_phi = 0;
   tree_electron_x = 0;
   tree_electron_y = 0;
   tree_electron_z = 0;
   tree_electron_px = 0;
   tree_electron_py = 0;
   tree_electron_pz = 0;
   tree_electron_energy = 0;
   tree_electron_charge = 0;
   tree_electron_isoR4 = 0;
   tree_electron_IsLoose = 0;
   tree_electron_IsMedium = 0;
   tree_electron_IsTight = 0;
   tree_electron_trigger_Ele = 0;
   tree_electron_trigger_diEle = 0;
   tree_electron_dxy = 0;
   tree_electron_dz = 0;
   tree_ST = 0;
   tree_muon_pt = 0;
   tree_muon_eta = 0;
   tree_muon_phi = 0;
   tree_muon_x = 0;
   tree_muon_y = 0;
   tree_muon_z = 0;
   tree_muon_px = 0;
   tree_muon_py = 0;
   tree_muon_pz = 0;
   tree_muon_energy = 0;
   tree_muon_dxy = 0;
   tree_muon_dxyError = 0;
   tree_muon_dz = 0;
   tree_muon_dzError = 0;
   tree_muon_charge = 0;
   tree_muon_isLoose = 0;
   tree_muon_isTight = 0;
   tree_muon_isGlobal = 0;
   tree_muon_isoR3 = 0;
   tree_muon_trigger_dimu = 0;
   tree_muon_trigger_isomu = 0;
   tree_muon_PFIsoLoose = 0;
   tree_muon_PFIsoMedium = 0;
   tree_muon_PFIsoTight = 0;
   tree_muon_nmu = 0;
   tree_muon_leadingpt = 0;
   tree_muon_leadingpt2 = 0;
   tree_muon_muon_dR = 0;
   tree_muon_muon_dPhi = 0;
   tree_muon_muon_dEta = 0;
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
   tree_track_dzTOpu = 0;
   tree_track_dzSigTOpu = 0;
   tree_track_algo = 0;
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
   tree_genAxis_dRneuneu = 0;
   tree_genAxis_dPhineuneu = 0;
   tree_genAxis_dEtaneuneu = 0;
   tree_GenAxes_Mass = 0;
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
   tree_Hemi_Mass = 0;
   tree_Hemi_njet = 0;
   tree_Hemi_eta = 0;
   tree_Hemi_phi = 0;
   tree_Hemi_dR = 0;
   tree_Hemi_nTrks = 0;
   tree_Hemi_nTrks_sig = 0;
   tree_Hemi_nTrks_bad = 0;
   tree_Hemi_LLP = 0;
   tree_Hemi_LLP_pt = 0;
   tree_Hemi_LLP_eta = 0;
   tree_Hemi_LLP_phi = 0;
   tree_Hemi_LLP_dist = 0;
   tree_Hemi_LLP_x = 0;
   tree_Hemi_LLP_y = 0;
   tree_Hemi_LLP_z = 0;
   tree_Hemi_Vtx_step = 0;
   tree_Hemi_Vtx_NChi2 = 0;
   tree_Hemi_Vtx_nTrks = 0;
   tree_Hemi_Vtx_nTrks_sig = 0;
   tree_Hemi_Vtx_nTrks_bad = 0;
   tree_Hemi_Vtx_x = 0;
   tree_Hemi_Vtx_y = 0;
   tree_Hemi_Vtx_r = 0;
   tree_Hemi_Vtx_z = 0;
   tree_Hemi_Vtx_xError = 0;
   tree_Hemi_Vtx_yError = 0;
   tree_Hemi_Vtx_zError = 0;
   tree_Hemi_Vtx_eta = 0;
   tree_Hemi_Vtx_Vtx_dr = 0;
   tree_Hemi_Vtx_Vtx_dz = 0;
   tree_Hemi_Vtx_Vtx_dd = 0;
   tree_Hemi_Vtx_BTag = 0;
   tree_Hemi_Vtx_nVtx = 0;
   tree_Hemi_Vtx_trackWeight = 0;
   tree_Hemi_Vtx_MeantrackWeight = 0;
   tree_Hemi_Vtx_track_DCA_x = 0;
   tree_Hemi_Vtx_track_DCA_y = 0;
   tree_Hemi_Vtx_track_DCA_z = 0;
   tree_Hemi_Vtx_track_DCA_r = 0;
   tree_Hemi_Vtx_track_DCA_d = 0;
   tree_Hemi_Vtx_track_MeanDCA_d = 0;
   tree_Hemi_Vtx_Mass = 0;
   tree_Hemi_Vtx_MVAval = 0;
   tree_Hemi_Vtx_MVAval_Step1 = 0;
   tree_Hemi_Vtx_TVtx_dx = 0;
   tree_Hemi_Vtx_TVtx_dy = 0;
   tree_Hemi_Vtx_TVtx_dz = 0;
   tree_Hemi_Vtx_TVtx_NChi2 = 0;
   tree_Hemi_Vtx_dist = 0;
   tree_Hemi_Vtx_dx = 0;
   tree_Hemi_Vtx_dy = 0;
   tree_Hemi_Vtx_dz = 0;
   tree_Hemi_Vtx_dr = 0;
   tree_Hemi_Vtx_dd = 0;
   tree_Hemi_dR12 = 0;
   tree_Hemi_LLP_dR12 = 0;
   tree_Hemi_Vtx_ddbad = 0;
   tree_Hemi_Vtx_ntrk10 = 0;
   tree_Hemi_Vtx_ntrk20 = 0;
   tree_track_Hemi_isjet = 0;
   tree_Hemi_Vtx_ddToBkg = 0;
   tree_Hemi_LLP_ping = 0;
   tree_event_LLP_ping = 0;
   tree_Hemi_LooseBTag_axes = 0;
   tree_Hemi_MediumBTag_axes = 0;
   tree_Hemi_TightBTag_axes = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("lumiBlock", &lumiBlock, &b_lumiBlock);
   fChain->SetBranchAddress("tree_nPV", &tree_nPV, &b_tree_nPV);
   fChain->SetBranchAddress("tree_PV_x", &tree_PV_x, &b_tree_PV_x);
   fChain->SetBranchAddress("tree_PV_y", &tree_PV_y, &b_tree_PV_y);
   fChain->SetBranchAddress("tree_PV_z", &tree_PV_z, &b_tree_PV_z);
   fChain->SetBranchAddress("tree_PV_ez", &tree_PV_ez, &b_tree_PV_ez);
   fChain->SetBranchAddress("tree_PV_NChi2", &tree_PV_NChi2, &b_tree_PV_NChi2);
   fChain->SetBranchAddress("tree_PV_ndf", &tree_PV_ndf, &b_tree_PV_ndf);
   fChain->SetBranchAddress("tree_Evts_MVAval",&tree_Evts_MVAval, &b_tree_Evts_MVAval);
   fChain->SetBranchAddress("tree_allPV_i", &tree_allPV_i, &b_tree_allPV_i);
   fChain->SetBranchAddress("tree_allPV_x", &tree_allPV_x, &b_tree_allPV_x);
   fChain->SetBranchAddress("tree_allPV_y", &tree_allPV_y, &b_tree_allPV_y);
   fChain->SetBranchAddress("tree_allPV_z", &tree_allPV_z, &b_tree_allPV_z);
   fChain->SetBranchAddress("tree_allPV_ex", &tree_allPV_ex, &b_tree_allPV_ex);
   fChain->SetBranchAddress("tree_allPV_ey", &tree_allPV_ey, &b_tree_allPV_ey);
   fChain->SetBranchAddress("tree_allPV_ez", &tree_allPV_ez, &b_tree_allPV_ez);
   fChain->SetBranchAddress("tree_allPV_NChi2", &tree_allPV_NChi2, &b_tree_allPV_NChi2);
   fChain->SetBranchAddress("tree_allPV_ndf", &tree_allPV_ndf, &b_tree_allPV_ndf);
   fChain->SetBranchAddress("tree_bs_PosX", &tree_bs_PosX, &b_tree_bs_PosX);
   fChain->SetBranchAddress("tree_bs_PosY", &tree_bs_PosY, &b_tree_bs_PosY);
   fChain->SetBranchAddress("tree_bs_PosZ", &tree_bs_PosZ, &b_tree_bs_PosZ);
   fChain->SetBranchAddress("tree_NbrOfZCand", &tree_NbrOfZCand, &b_tree_NbrOfZCand);
   fChain->SetBranchAddress("tree_Filter", &tree_Filter, &b_tree_Filter);
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
   fChain->SetBranchAddress("tree_PFMet_et", &tree_PFMet_et, &b_tree_PFMet_et);
   fChain->SetBranchAddress("tree_PFMet_phi", &tree_PFMet_phi, &b_tree_PFMet_phi);
   fChain->SetBranchAddress("tree_PFMet_sig", &tree_PFMet_sig, &b_tree_PFMet_sig);
   fChain->SetBranchAddress("tree_njet", &tree_njet, &b_tree_njet);
   fChain->SetBranchAddress("tree_jet_E", &tree_jet_E, &b_tree_jet_E);
   fChain->SetBranchAddress("tree_jet_pt", &tree_jet_pt, &b_tree_jet_pt);
   fChain->SetBranchAddress("tree_jet_eta", &tree_jet_eta, &b_tree_jet_eta);
   fChain->SetBranchAddress("tree_jet_phi", &tree_jet_phi, &b_tree_jet_phi);
   fChain->SetBranchAddress("tree_jet_btag_DeepCSV", &tree_jet_btag_DeepCSV, &b_tree_jet_btag_DeepCSV);
   fChain->SetBranchAddress("tree_jet_btag_DeepJet", &tree_jet_btag_DeepJet, &b_tree_jet_btag_DeepJet);
   fChain->SetBranchAddress("tree_jet_leadingpt", &tree_jet_leadingpt, &b_tree_jet_leadingpt);
   fChain->SetBranchAddress("tree_jet_leadingpt2", &tree_jet_leadingpt2, &b_tree_jet_leadingpt2);
   fChain->SetBranchAddress("tree_jet_leadingMuon_dR", &tree_jet_leadingMuon_dR, &b_tree_jet_leadingMuon_dR);
   fChain->SetBranchAddress("tree_jet_leadingMuon2_dR", &tree_jet_leadingMuon2_dR, &b_tree_jet_leadingMuon2_dR);
   fChain->SetBranchAddress("tree_jet_jet_dR", &tree_jet_jet_dR, &b_tree_jet_jet_dR);
   fChain->SetBranchAddress("tree_jet_jet_dPhi", &tree_jet_jet_dPhi, &b_tree_jet_jet_dPhi);
   fChain->SetBranchAddress("tree_jet_jet_dEta", &tree_jet_jet_dEta, &b_tree_jet_jet_dEta);
   fChain->SetBranchAddress("tree_muon_jet_dRmin", &tree_muon_jet_dRmin, &b_tree_muon_jet_dRmin);
   fChain->SetBranchAddress("tree_muon_jet_dRmax", &tree_muon_jet_dRmax, &b_tree_muon_jet_dRmax);
   fChain->SetBranchAddress("tree_HT", &tree_HT, &b_tree_HT);
   fChain->SetBranchAddress("tree_electron_nEle", &tree_electron_nEle, &b_tree_electron_nEle);
   fChain->SetBranchAddress("tree_electron_pt", &tree_electron_pt, &b_tree_electron_pt);
   fChain->SetBranchAddress("tree_electron_eta", &tree_electron_eta, &b_tree_electron_eta);
   fChain->SetBranchAddress("tree_electron_phi", &tree_electron_phi, &b_tree_electron_phi);
   fChain->SetBranchAddress("tree_electron_x", &tree_electron_x, &b_tree_electron_x);
   fChain->SetBranchAddress("tree_electron_y", &tree_electron_y, &b_tree_electron_y);
   fChain->SetBranchAddress("tree_electron_z", &tree_electron_z, &b_tree_electron_z);
   fChain->SetBranchAddress("tree_electron_px", &tree_electron_px, &b_tree_electron_px);
   fChain->SetBranchAddress("tree_electron_py", &tree_electron_py, &b_tree_electron_py);
   fChain->SetBranchAddress("tree_electron_pz", &tree_electron_pz, &b_tree_electron_pz);
   fChain->SetBranchAddress("tree_electron_energy", &tree_electron_energy, &b_tree_electron_energy);
   fChain->SetBranchAddress("tree_electron_charge", &tree_electron_charge, &b_tree_electron_charge);
   fChain->SetBranchAddress("tree_electron_isoR4", &tree_electron_isoR4, &b_tree_electron_isoR4);
   fChain->SetBranchAddress("tree_electron_IsLoose", &tree_electron_IsLoose, &b_tree_electron_IsLoose);
   fChain->SetBranchAddress("tree_electron_IsMedium", &tree_electron_IsMedium, &b_tree_electron_IsMedium);
   fChain->SetBranchAddress("tree_electron_IsTight", &tree_electron_IsTight, &b_tree_electron_IsTight);
   fChain->SetBranchAddress("tree_electron_trigger_Ele", &tree_electron_trigger_Ele, &b_tree_electron_trigger_Ele);
   fChain->SetBranchAddress("tree_electron_trigger_diEle", &tree_electron_trigger_diEle, &b_tree_electron_trigger_diEle);
   fChain->SetBranchAddress("tree_electron_dxy", &tree_electron_dxy, &b_tree_electron_dxy);
   fChain->SetBranchAddress("tree_electron_dz", &tree_electron_dz, &b_tree_electron_dz);
   fChain->SetBranchAddress("tree_ST", &tree_ST, &b_tree_ST);
   fChain->SetBranchAddress("tree_Mmumu", &tree_Mmumu, &b_tree_Mmumu);
   fChain->SetBranchAddress("tree_muon_pt", &tree_muon_pt, &b_tree_muon_pt);
   fChain->SetBranchAddress("tree_muon_eta", &tree_muon_eta, &b_tree_muon_eta);
   fChain->SetBranchAddress("tree_muon_phi", &tree_muon_phi, &b_tree_muon_phi);
   fChain->SetBranchAddress("tree_muon_x", &tree_muon_x, &b_tree_muon_x);
   fChain->SetBranchAddress("tree_muon_y", &tree_muon_y, &b_tree_muon_y);
   fChain->SetBranchAddress("tree_muon_z", &tree_muon_z, &b_tree_muon_z);
   fChain->SetBranchAddress("tree_muon_px", &tree_muon_px, &b_tree_muon_px);
   fChain->SetBranchAddress("tree_muon_py", &tree_muon_py, &b_tree_muon_py);
   fChain->SetBranchAddress("tree_muon_pz", &tree_muon_pz, &b_tree_muon_pz);
   fChain->SetBranchAddress("tree_muon_energy", &tree_muon_energy, &b_tree_muon_energy);
   fChain->SetBranchAddress("tree_muon_dxy", &tree_muon_dxy, &b_tree_muon_dxy);
   fChain->SetBranchAddress("tree_muon_dxyError", &tree_muon_dxyError, &b_tree_muon_dxyError);
   fChain->SetBranchAddress("tree_muon_dz", &tree_muon_dz, &b_tree_muon_dz);
   fChain->SetBranchAddress("tree_muon_dzError", &tree_muon_dzError, &b_tree_muon_dzError);
   fChain->SetBranchAddress("tree_muon_charge", &tree_muon_charge, &b_tree_muon_charge);
   fChain->SetBranchAddress("tree_muon_isLoose", &tree_muon_isLoose, &b_tree_muon_isLoose);
   fChain->SetBranchAddress("tree_muon_isTight", &tree_muon_isTight, &b_tree_muon_isTight);
   fChain->SetBranchAddress("tree_muon_isGlobal", &tree_muon_isGlobal, &b_tree_muon_isGlobal);
   fChain->SetBranchAddress("tree_muon_isoR3", &tree_muon_isoR3, &b_tree_muon_isoR3);
   fChain->SetBranchAddress("tree_muon_trigger_dimu", &tree_muon_trigger_dimu, &b_tree_muon_trigger_dimu);
   fChain->SetBranchAddress("tree_muon_trigger_isomu", &tree_muon_trigger_isomu, &b_tree_muon_trigger_isomu);
   fChain->SetBranchAddress("tree_muon_PFIsoLoose", &tree_muon_PFIsoLoose, &b_tree_muon_PFIsoLoose);
   fChain->SetBranchAddress("tree_muon_PFIsoMedium", &tree_muon_PFIsoMedium, &b_tree_muon_PFIsoMedium);
   fChain->SetBranchAddress("tree_muon_PFIsoTight", &tree_muon_PFIsoTight, &b_tree_muon_PFIsoTight);
   fChain->SetBranchAddress("tree_muon_nmu", &tree_muon_nmu, &b_tree_muon_nmu);
   fChain->SetBranchAddress("tree_muon_leadingpt", &tree_muon_leadingpt, &b_tree_muon_leadingpt);
   fChain->SetBranchAddress("tree_muon_leadingpt2", &tree_muon_leadingpt2, &b_tree_muon_leadingpt2);
   fChain->SetBranchAddress("tree_muon_muon_dR", &tree_muon_muon_dR, &b_tree_muon_muon_dR);
   fChain->SetBranchAddress("tree_muon_muon_dPhi", &tree_muon_muon_dPhi, &b_tree_muon_muon_dPhi);
   fChain->SetBranchAddress("tree_muon_muon_dEta", &tree_muon_muon_dEta, &b_tree_muon_muon_dEta);
   fChain->SetBranchAddress("tree_TRACK_SIZE",&tree_TRACK_SIZE,&b_tree_TRACK_SIZE);
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
   fChain->SetBranchAddress("tree_track_dzTOpu", &tree_track_dzTOpu, &b_tree_track_dzTOpu);
   fChain->SetBranchAddress("tree_track_dzSigTOpu", &tree_track_dzSigTOpu, &b_tree_track_dzSigTOpu);
   fChain->SetBranchAddress("tree_track_algo", &tree_track_algo, &b_tree_track_algo);
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
   fChain->SetBranchAddress("tree_GenPVx", &tree_GenPVx, &b_tree_GenPVx);
   fChain->SetBranchAddress("tree_GenPVy", &tree_GenPVy, &b_tree_GenPVy);
   fChain->SetBranchAddress("tree_GenPVz", &tree_GenPVz, &b_tree_GenPVz);
   fChain->SetBranchAddress("tree_smu_mass",&tree_smu_mass,&b_tree_smu_mass);
   fChain->SetBranchAddress("tree_neu_mass",&tree_neu_mass,&b_tree_neu_mass);
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
   fChain->SetBranchAddress("tree_genAxis_dRneuneu", &tree_genAxis_dRneuneu, &b_tree_genAxis_dRneuneu);
   fChain->SetBranchAddress("tree_genAxis_dPhineuneu", &tree_genAxis_dPhineuneu, &b_tree_genAxis_dPhineuneu);
   fChain->SetBranchAddress("tree_genAxis_dEtaneuneu", &tree_genAxis_dEtaneuneu, &b_tree_genAxis_dEtaneuneu);
   fChain->SetBranchAddress("tree_GenAxes_Mass", &tree_GenAxes_Mass, &b_tree_GenAxes_Mass);
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
   fChain->SetBranchAddress("tree_Hemi_Mass", &tree_Hemi_Mass, &b_tree_Hemi_Mass);
   fChain->SetBranchAddress("tree_Hemi_njet", &tree_Hemi_njet, &b_tree_Hemi_njet);
   fChain->SetBranchAddress("tree_Hemi_eta", &tree_Hemi_eta, &b_tree_Hemi_eta);
   fChain->SetBranchAddress("tree_Hemi_phi", &tree_Hemi_phi, &b_tree_Hemi_phi);
   fChain->SetBranchAddress("tree_Hemi_dR", &tree_Hemi_dR, &b_tree_Hemi_dR);
   fChain->SetBranchAddress("tree_Hemi_nTrks", &tree_Hemi_nTrks, &b_tree_Hemi_nTrks);
   fChain->SetBranchAddress("tree_Hemi_nTrks_sig", &tree_Hemi_nTrks_sig, &b_tree_Hemi_nTrks_sig);
   fChain->SetBranchAddress("tree_Hemi_nTrks_bad", &tree_Hemi_nTrks_bad, &b_tree_Hemi_nTrks_bad);
   fChain->SetBranchAddress("tree_Hemi_LLP", &tree_Hemi_LLP, &b_tree_Hemi_LLP);
   fChain->SetBranchAddress("tree_Hemi_LLP_pt", &tree_Hemi_LLP_pt, &b_tree_Hemi_LLP_pt);
   fChain->SetBranchAddress("tree_Hemi_LLP_eta", &tree_Hemi_LLP_eta, &b_tree_Hemi_LLP_eta);
   fChain->SetBranchAddress("tree_Hemi_LLP_phi", &tree_Hemi_LLP_phi, &b_tree_Hemi_LLP_phi);
   fChain->SetBranchAddress("tree_Hemi_LLP_dist", &tree_Hemi_LLP_dist, &b_tree_Hemi_LLP_dist);
   fChain->SetBranchAddress("tree_Hemi_LLP_x", &tree_Hemi_LLP_x, &b_tree_Hemi_LLP_x);
   fChain->SetBranchAddress("tree_Hemi_LLP_y", &tree_Hemi_LLP_y, &b_tree_Hemi_LLP_y);
   fChain->SetBranchAddress("tree_Hemi_LLP_z", &tree_Hemi_LLP_z, &b_tree_Hemi_LLP_z);
   fChain->SetBranchAddress("tree_Hemi_Vtx_step", &tree_Hemi_Vtx_step, &b_tree_Hemi_Vtx_step);
   fChain->SetBranchAddress("tree_Hemi_Vtx_NChi2", &tree_Hemi_Vtx_NChi2, &b_tree_Hemi_Vtx_NChi2);
   fChain->SetBranchAddress("tree_Hemi_Vtx_nTrks", &tree_Hemi_Vtx_nTrks, &b_tree_Hemi_Vtx_nTrks);
   fChain->SetBranchAddress("tree_Hemi_Vtx_nTrks_sig", &tree_Hemi_Vtx_nTrks_sig, &b_tree_Hemi_Vtx_nTrks_sig);
   fChain->SetBranchAddress("tree_Hemi_Vtx_nTrks_bad", &tree_Hemi_Vtx_nTrks_bad, &b_tree_Hemi_Vtx_nTrks_bad);
   fChain->SetBranchAddress("tree_Hemi_Vtx_x", &tree_Hemi_Vtx_x, &b_tree_Hemi_Vtx_x);
   fChain->SetBranchAddress("tree_Hemi_Vtx_y", &tree_Hemi_Vtx_y, &b_tree_Hemi_Vtx_y);
   fChain->SetBranchAddress("tree_Hemi_Vtx_r", &tree_Hemi_Vtx_r, &b_tree_Hemi_Vtx_r);
   fChain->SetBranchAddress("tree_Hemi_Vtx_z", &tree_Hemi_Vtx_z, &b_tree_Hemi_Vtx_z);
   fChain->SetBranchAddress("tree_Hemi_Vtx_xError", &tree_Hemi_Vtx_xError, &b_tree_Hemi_Vtx_xError);
   fChain->SetBranchAddress("tree_Hemi_Vtx_yError", &tree_Hemi_Vtx_yError, &b_tree_Hemi_Vtx_yError);
   fChain->SetBranchAddress("tree_Hemi_Vtx_zError", &tree_Hemi_Vtx_zError, &b_tree_Hemi_Vtx_zError);
   fChain->SetBranchAddress("tree_Hemi_Vtx_eta", &tree_Hemi_Vtx_eta, &b_tree_Hemi_Vtx_eta);
   fChain->SetBranchAddress("tree_Hemi_Vtx_Vtx_dr", &tree_Hemi_Vtx_Vtx_dr, &b_tree_Hemi_Vtx_Vtx_dr);
   fChain->SetBranchAddress("tree_Hemi_Vtx_Vtx_dz", &tree_Hemi_Vtx_Vtx_dz, &b_tree_Hemi_Vtx_Vtx_dz);
   fChain->SetBranchAddress("tree_Hemi_Vtx_Vtx_dd", &tree_Hemi_Vtx_Vtx_dd, &b_tree_Hemi_Vtx_Vtx_dd);
   fChain->SetBranchAddress("tree_Hemi_Vtx_BTag", &tree_Hemi_Vtx_BTag, &b_tree_Hemi_Vtx_BTag);
   fChain->SetBranchAddress("tree_Hemi_Vtx_nVtx", &tree_Hemi_Vtx_nVtx, &b_tree_Hemi_Vtx_nVtx);
   fChain->SetBranchAddress("tree_Hemi_Vtx_trackWeight", &tree_Hemi_Vtx_trackWeight, &b_tree_Hemi_Vtx_trackWeight);
   fChain->SetBranchAddress("tree_Hemi_Vtx_MeantrackWeight", &tree_Hemi_Vtx_MeantrackWeight, &b_tree_Hemi_Vtx_MeantrackWeight);
   fChain->SetBranchAddress("tree_Hemi_Vtx_track_DCA_x", &tree_Hemi_Vtx_track_DCA_x, &b_tree_Hemi_Vtx_track_DCA_x);
   fChain->SetBranchAddress("tree_Hemi_Vtx_track_DCA_y", &tree_Hemi_Vtx_track_DCA_y, &b_tree_Hemi_Vtx_track_DCA_y);
   fChain->SetBranchAddress("tree_Hemi_Vtx_track_DCA_z", &tree_Hemi_Vtx_track_DCA_z, &b_tree_Hemi_Vtx_track_DCA_z);
   fChain->SetBranchAddress("tree_Hemi_Vtx_track_DCA_r", &tree_Hemi_Vtx_track_DCA_r, &b_tree_Hemi_Vtx_track_DCA_r);
   fChain->SetBranchAddress("tree_Hemi_Vtx_track_DCA_d", &tree_Hemi_Vtx_track_DCA_d, &b_tree_Hemi_Vtx_track_DCA_d);
   fChain->SetBranchAddress("tree_Hemi_Vtx_track_MeanDCA_d", &tree_Hemi_Vtx_track_MeanDCA_d, &b_tree_Hemi_Vtx_track_MeanDCA_d);
   fChain->SetBranchAddress("tree_Hemi_Vtx_Mass", &tree_Hemi_Vtx_Mass, &b_tree_Hemi_Vtx_Mass);
   fChain->SetBranchAddress("tree_Hemi_Vtx_MVAval", &tree_Hemi_Vtx_MVAval, &b_tree_Hemi_Vtx_MVAval);
   fChain->SetBranchAddress("tree_Hemi_Vtx_MVAval_Step1", &tree_Hemi_Vtx_MVAval_Step1, &b_tree_Hemi_Vtx_MVAval_Step1);
   fChain->SetBranchAddress("tree_Hemi_Vtx_TVtx_dx", &tree_Hemi_Vtx_TVtx_dx, &b_tree_Hemi_Vtx_TVtx_dx);
   fChain->SetBranchAddress("tree_Hemi_Vtx_TVtx_dy", &tree_Hemi_Vtx_TVtx_dy, &b_tree_Hemi_Vtx_TVtx_dy);
   fChain->SetBranchAddress("tree_Hemi_Vtx_TVtx_dz", &tree_Hemi_Vtx_TVtx_dz, &b_tree_Hemi_Vtx_TVtx_dz);
   fChain->SetBranchAddress("tree_Hemi_Vtx_TVtx_NChi2", &tree_Hemi_Vtx_TVtx_NChi2, &b_tree_Hemi_Vtx_TVtx_NChi2);
   fChain->SetBranchAddress("tree_Hemi_Vtx_dist", &tree_Hemi_Vtx_dist, &b_tree_Hemi_Vtx_dist);
   fChain->SetBranchAddress("tree_Hemi_Vtx_dx", &tree_Hemi_Vtx_dx, &b_tree_Hemi_Vtx_dx);
   fChain->SetBranchAddress("tree_Hemi_Vtx_dy", &tree_Hemi_Vtx_dy, &b_tree_Hemi_Vtx_dy);
   fChain->SetBranchAddress("tree_Hemi_Vtx_dz", &tree_Hemi_Vtx_dz, &b_tree_Hemi_Vtx_dz);
   fChain->SetBranchAddress("tree_Hemi_Vtx_dr", &tree_Hemi_Vtx_dr, &b_tree_Hemi_Vtx_dr);
   fChain->SetBranchAddress("tree_Hemi_Vtx_dd", &tree_Hemi_Vtx_dd, &b_tree_Hemi_Vtx_dd);
   fChain->SetBranchAddress("tree_Hemi_dR12", &tree_Hemi_dR12, &b_tree_Hemi_dR12);
   fChain->SetBranchAddress("tree_Hemi_LLP_dR12", &tree_Hemi_LLP_dR12, &b_tree_Hemi_LLP_dR12);
   fChain->SetBranchAddress("tree_Hemi_Vtx_ddbad", &tree_Hemi_Vtx_ddbad, &b_tree_Hemi_Vtx_ddbad);
   fChain->SetBranchAddress("tree_Hemi_Vtx_ntrk10", &tree_Hemi_Vtx_ntrk10, &b_tree_Hemi_Vtx_ntrk10);
   fChain->SetBranchAddress("tree_Hemi_Vtx_ntrk20", &tree_Hemi_Vtx_ntrk20, &b_tree_Hemi_Vtx_ntrk20);
   fChain->SetBranchAddress("tree_track_Hemi_isjet", &tree_track_Hemi_isjet, &b_tree_track_Hemi_isjet);
   fChain->SetBranchAddress("tree_Hemi_Vtx_ddToBkg", &tree_Hemi_Vtx_ddToBkg, &b_tree_Hemi_Vtx_ddToBkg);
   fChain->SetBranchAddress("tree_Hemi_LLP_ping", &tree_Hemi_LLP_ping, &b_tree_Hemi_LLP_ping);
   fChain->SetBranchAddress("tree_event_LLP_ping", &tree_event_LLP_ping, &b_tree_event_LLP_ping);
   fChain->SetBranchAddress("tree_Hemi_LooseBTag_axes", &tree_Hemi_LooseBTag_axes, &b_tree_Hemi_LooseBTag_axes);
   fChain->SetBranchAddress("tree_Hemi_MediumBTag_axes", &tree_Hemi_MediumBTag_axes, &b_tree_Hemi_MediumBTag_axes);
   fChain->SetBranchAddress("tree_Hemi_TightBTag_axes", &tree_Hemi_TightBTag_axes, &b_tree_Hemi_TightBTag_axes);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf_v", &HLT_Ele27_WPTight_Gsf_v, &b_HLT_Ele27_WPTight_Gsf_v);
   fChain->SetBranchAddress("HLT_Ele32_WPTight_Gsf_v", &HLT_Ele32_WPTight_Gsf_v, &b_HLT_Ele32_WPTight_Gsf_v);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_v", &HLT_PFMET120_PFMHT120_IDTight_v, &b_HLT_PFMET120_PFMHT120_IDTight_v);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v", &HLT_PFMET120_PFMHT120_IDTight_PFHT60_v, &b_HLT_PFMET120_PFMHT120_IDTight_PFHT60_v);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
   fChain->SetBranchAddress("HLT_PFMET250_HBHECleaned_v", &HLT_PFMET250_HBHECleaned_v, &b_HLT_PFMET250_HBHECleaned_v);
   fChain->SetBranchAddress("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v", &HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v, &b_HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v);
   fChain->SetBranchAddress("HLT_IsoMu24_v", &HLT_IsoMu24_v, &b_HLT_IsoMu24_v);
   fChain->SetBranchAddress("HLT_IsoMu27_v", &HLT_IsoMu27_v, &b_HLT_IsoMu27_v);
   fChain->SetBranchAddress("HLT_PFHT180_v", &HLT_PFHT180_v, &b_HLT_PFHT180_v);
   fChain->SetBranchAddress("HLT_PFHT250_v", &HLT_PFHT250_v, &b_HLT_PFHT250_v);
   fChain->SetBranchAddress("HLT_PFHT370_v", &HLT_PFHT370_v, &b_HLT_PFHT370_v);
   fChain->SetBranchAddress("HLT_PFHT430_v", &HLT_PFHT430_v, &b_HLT_PFHT430_v);
   fChain->SetBranchAddress("HLT_PFHT510_v", &HLT_PFHT510_v, &b_HLT_PFHT510_v);
   fChain->SetBranchAddress("HLT_PFHT590_v", &HLT_PFHT590_v, &b_HLT_PFHT590_v);
   fChain->SetBranchAddress("HLT_PFHT680_v", &HLT_PFHT680_v, &b_HLT_PFHT680_v);
   fChain->SetBranchAddress("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v", &HLT_PFHT500_PFMET100_PFMHT100_IDTight_v, &b_HLT_PFHT500_PFMET100_PFMHT100_IDTight_v);
   fChain->SetBranchAddress("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v", &HLT_PFHT700_PFMET85_PFMHT85_IDTight_v, &b_HLT_PFHT700_PFMET85_PFMHT85_IDTight_v);
   fChain->SetBranchAddress("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v", &HLT_PFHT800_PFMET75_PFMHT75_IDTight_v, &b_HLT_PFHT800_PFMET75_PFMHT75_IDTight_v);
   Notify();
}

Bool_t TreeReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TreeReader_cxx
