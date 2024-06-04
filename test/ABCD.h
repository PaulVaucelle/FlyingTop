//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 26 08:15:26 2024 by ROOT version 6.14/09
// from TTree ttree/summary information
// found on file: MiniDYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root
//////////////////////////////////////////////////////////

#ifndef ABCD_h
#define ABCD_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

//$$
#include <TMath.h>
using namespace std;
//$$

class ABCD {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<int>     *minirunNumber;
   vector<int>     *minieventNumber;
   vector<int>     *minilumiBlock;
   vector<float>   *minitree_LHE_Weights;
   vector<float>   *minitree_MCEvt_weight;
   vector<double>  *minitree_only_gen_wt;
   vector<double>  *minitree_event_weight;
   vector<double>  *minitree_genTop_Weight;
   vector<double>  *miniPUweight;
   vector<double>  *miniPrefweight;
   vector<int>     *miniPU_events;
   vector<int>     *miniAllPU_events_weight;
   vector<bool>    *minitree_Filter;
   vector<bool>    *minitree_FilterSameSign;
   vector<float>   *minitree_GenPVx;
   vector<float>   *minitree_GenPVy;
   vector<float>   *minitree_GenPVz;
   vector<int>     *minitree_smu_mass;
   vector<int>     *minitree_neu_mass;
   vector<int>     *minitree_neu_ctau;
   vector<bool>    *minitree_Good_PV;
   vector<int>     *minitree_nPV;
   vector<float>   *minitree_PV_x;
   vector<float>   *minitree_PV_y;
   vector<float>   *minitree_PV_z;
   vector<float>   *minitree_PV_ez;
   vector<float>   *minitree_PV_NChi2;
   vector<int>     *minitree_PV_ndf;
   vector<float>   *minitree_PFMet_et;
   vector<float>   *minitree_PFMet_phi;
   vector<float>   *minitree_HT;
   vector<int>     *minitree_TRACK_SIZE;
   vector<int>     *minitree_nTracks;
   vector<int>     *minitree_nLostTracks;
   vector<int>     *minitree_muon_GenRecoTriggerMatched;
   vector<int>     *minitree_all_nmu;
   vector<int>     *minitree_nmu;
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
   vector<int>     *minitree_muon_gen;
   vector<float>   *minitree_lepton_leadingpt;
   vector<float>   *minitree_lepton_leadingpt2;
   vector<float>   *minitree_lepton_lepton_dR;
   vector<float>   *minitree_lepton_lepton_dPhi;
   vector<float>   *minitree_lepton_lepton_dEta;
   vector<float>   *minitree_ll_pt;
   vector<float>   *minitree_ll_eta;
   vector<float>   *minitree_ll_phi;
   vector<float>   *minitree_ll_px;
   vector<float>   *minitree_ll_py;
   vector<float>   *minitree_ll_pz;
   vector<float>   *minitree_ll_energy;
   vector<float>   *minitree_ll_mass;
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
   vector<float>   *minitree_Evts_MVAval;
   vector<float>   *minitree_Evts_MVAvalDY;
   vector<float>   *minitree_Evts_MVAvalTT;
   vector<int>     *minitree_nLLP;
   vector<int>     *minitree_LLP;
   vector<float>   *minitree_LLP_pt;
   vector<float>   *minitree_LLP_eta;
   vector<float>   *minitree_LLP_phi;
   vector<float>   *minitree_LLP_x;
   vector<float>   *minitree_LLP_y;
   vector<float>   *minitree_LLP_z;
   vector<float>   *minitree_LLP_r;
   vector<float>   *minitree_LLP_dist;
   vector<int>     *minitree_LLP_nTrks;
   vector<float>   *minitree_LLP12_dR;
   vector<float>   *minitree_LLP12_deta;
   vector<float>   *minitree_LLP12_dphi;
   vector<float>   *minitree_LLP_Mass;
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
   vector<int>     *minitree_Hemi_LLP;
   vector<float>   *minitree_Hemi_LLP_pt;
   vector<float>   *minitree_Hemi_LLP_eta;
   vector<float>   *minitree_Hemi_LLP_phi;
   vector<float>   *minitree_Hemi_LLP_dist;
   vector<float>   *minitree_Hemi_LLP_x;
   vector<float>   *minitree_Hemi_LLP_y;
   vector<float>   *minitree_Hemi_LLP_z;
   vector<float>   *minitree_Hemi_LLP_dR;
   vector<int>     *minitree_Hemi_LLP_mother;
   vector<float>   *minitree_Hemi_LLP_Vtx_dx;
   vector<float>   *minitree_Hemi_LLP_Vtx_dy;
   vector<float>   *minitree_Hemi_LLP_Vtx_dz;
   vector<float>   *minitree_Hemi_LLP_Vtx_dr;
   vector<float>   *minitree_Hemi_LLP_muOK_dR;
   vector<float>   *minitree_Hemi_LLP_muOK_pt;
   vector<float>   *minitree_Hemi_LLP_muOK_mass;
   vector<float>   *minitree_Hemi_LLP_muNO_dR;
   vector<float>   *minitree_Hemi_LLP_muNO_pt;
   vector<float>   *minitree_Hemi_LLP_muNO_mass;
   vector<float>   *minitree_Hemi_LLP_dR12;
   vector<bool>    *minitree_Hemi_LLP_ping;
   vector<int>     *minitree_event_LLP_ping;
   vector<int>     *minitree_Hemi_Vtx_step;
   vector<bool>    *minitree_Hemi_Vtx_isTight;
   vector<float>   *minitree_Hemi_Vtx_NChi2;
   vector<int>     *minitree_Hemi_Vtx_nTrks;
   vector<int>     *minitree_Hemi_Vtx_nTrks_sig;
   vector<int>     *minitree_Hemi_Vtx_nTrks_bad;
   vector<float>   *minitree_Hemi_Vtx_x;
   vector<float>   *minitree_Hemi_Vtx_y;
   vector<float>   *minitree_Hemi_Vtx_z;
   vector<float>   *minitree_Hemi_Vtx_r;
   vector<float>   *minitree_Hemi_Vtx_dR;
   vector<float>   *minitree_Hemi_Vtx_SumtrackWeight;
   vector<float>   *minitree_Hemi_Vtx_track_MeanDCA_d;
   vector<float>   *minitree_Hemi_Vtx_Mass;
   vector<float>   *minitree_Hemi_Vtx_dist;
   vector<int>     *minitree_event_nVtx;
   vector<float>   *minitree_event_Vtx_Vtx_dr;
   vector<float>   *minitree_event_Vtx_Vtx_dz;
   vector<float>   *minitree_event_Vtx_Vtx_dd;
   vector<float>   *minitree_event_Vtx_Vtx_reldd;
   vector<float>   *minitree_event_Vtx_Vtx_dR;
   vector<int>     *minitree_event_Vtx_Vtx_step;
   vector<float>   *minitree_Hemi_SecLLP;
   vector<float>   *minitree_Hemi_LLP_SecVtx_dz;
   vector<float>   *minitree_Hemi_LLP_SecVtx_dr;
   vector<bool>    *minitree_Hemi_SecLLP_ping;
   vector<int>     *minitree_event_SecLLP_ping;
   vector<int>     *minitree_Hemi_SecVtx;
   vector<int>     *minitree_Hemi_SecVtx_step;
   vector<float>   *minitree_Hemi_SecVtx_x;
   vector<float>   *minitree_Hemi_SecVtx_y;
   vector<float>   *minitree_Hemi_SecVtx_z;
   vector<float>   *minitree_Hemi_SecVtx_r;
   vector<float>   *minitree_Hemi_SecVtx_dR;
   vector<float>   *minitree_Hemi_SecVtx_nTrks;
   vector<float>   *minitree_Hemi_SecVtx_NChi2;
   vector<float>   *minitree_Hemi_SecVtx_dist;
   vector<float>   *minitree_Hemi_SecVtx_track_MeanDCA_d;
   vector<float>   *minitree_Hemi_SecVtx_SumtrackWeight;
   vector<float>   *minitree_Hemi_SecVtx_Mass;
   vector<float>   *minitree_event_MergedVtx_Vtx_dr;
   vector<float>   *minitree_event_MergedVtx_Vtx_dz;
   vector<float>   *minitree_event_MergedVtx_Vtx_dd;
   vector<float>   *minitree_event_MergedVtx_Vtx_reldd;
   vector<float>   *minitree_event_MergedVtx_Vtx_dR;
   vector<int>     *minitree_event_MergedVtx_Vtx_step;
   vector<float>   *minitree_Hemi_Vtx_BDT_nTrks;
   vector<float>   *minitree_Hemi_Vtx_BDT_NChi2;
   vector<float>   *minitree_Hemi_Vtx_BDT_step;
   vector<float>   *minitree_Hemi_Vtx_BDT_STW;
   vector<float>   *minitree_Hemi_Vtx_BDT_Mass;
   vector<float>   *minitree_Hemi_Vtx_BDT_HMass;
   vector<float>   *minitree_Hemi_Vtx_BDT_ntrk10;
   vector<float>   *minitree_Hemi_Vtx_BDT_ntrk20;
   vector<float>   *minitree_Hemi_Vtx_BDT_MeanDCA;
   vector<float>   *minitree_Hemi_Vtx_MVAval_Loose;
   vector<float>   *minitree_Hemi_Vtx_MVAval_Tight;
   vector<bool>    *miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;
   vector<bool>    *miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v;
   vector<bool>    *miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v;
   vector<bool>    *miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   vector<bool>    *miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;
   vector<bool>    *miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   vector<bool>    *miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   vector<bool>    *miniHLT_Ele27_WPTight_Gsf_v;
   vector<bool>    *miniHLT_Ele32_WPTight_Gsf_v;
   vector<bool>    *miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   vector<bool>    *miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   vector<bool>    *miniHLT_PFMET120_PFMHT120_IDTight_v;
   vector<bool>    *miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v;
   vector<bool>    *miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;
   vector<bool>    *miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;
   vector<bool>    *miniHLT_PFMET250_HBHECleaned_v;
   vector<bool>    *miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;
   vector<bool>    *miniHLT_IsoMu24_v;
   vector<bool>    *miniHLT_IsoMu27_v;

   // List of branches
   TBranch        *b_minirunNumber;   //!
   TBranch        *b_minieventNumber;   //!
   TBranch        *b_minilumiBlock;   //!
   TBranch        *b_minitree_LHE_Weights;   //!
   TBranch        *b_minitree_MCEvt_weight;   //!
   TBranch        *b_minitree_only_gen_wt;   //!
   TBranch        *b_minitree_event_weight;   //!
   TBranch        *b_minitree_genTop_Weight;   //!
   TBranch        *b_miniPUweight;   //!
   TBranch        *b_miniPrefweight;   //!
   TBranch        *b_miniPU_events;   //!
   TBranch        *b_miniAllPU_events_weight;   //!
   TBranch        *b_minitree_Filter;   //!
   TBranch        *b_minitree_FilterSameSign;   //!
   TBranch        *b_minitree_GenPVx;   //!
   TBranch        *b_minitree_GenPVy;   //!
   TBranch        *b_minitree_GenPVz;   //!
   TBranch        *b_minitree_smu_mass;   //!
   TBranch        *b_minitree_neu_mass;   //!
   TBranch        *b_minitree_neu_ctau;   //!
   TBranch        *b_minitree_Good_PV;   //!
   TBranch        *b_minitree_nPV;   //!
   TBranch        *b_minitree_PV_x;   //!
   TBranch        *b_minitree_PV_y;   //!
   TBranch        *b_minitree_PV_z;   //!
   TBranch        *b_minitree_PV_ez;   //!
   TBranch        *b_minitree_PV_NChi2;   //!
   TBranch        *b_minitree_PV_ndf;   //!
   TBranch        *b_minitree_PFMet_et;   //!
   TBranch        *b_minitree_PFMet_phi;   //!
   TBranch        *b_minitree_HT;   //!
   TBranch        *b_minitree_TRACK_SIZE;   //!
   TBranch        *b_minitree_nTracks;   //!
   TBranch        *b_minitree_nLostTracks;   //!
   TBranch        *b_minitree_muon_GenRecoTriggerMatched;   //!
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
   TBranch        *b_minitree_muon_gen;   //!
   TBranch        *b_minitree_lepton_leadingpt;   //!
   TBranch        *b_minitree_lepton_leadingpt2;   //!
   TBranch        *b_minitree_lepton_lepton_dR;   //!
   TBranch        *b_minitree_lepton_lepton_dPhi;   //!
   TBranch        *b_minitree_lepton_lepton_dEta;   //!
   TBranch        *b_minitree_ll_pt;   //!
   TBranch        *b_minitree_ll_eta;   //!
   TBranch        *b_minitree_ll_phi;   //!
   TBranch        *b_minitree_ll_px;   //!
   TBranch        *b_minitree_ll_py;   //!
   TBranch        *b_minitree_ll_pz;   //!
   TBranch        *b_minitree_ll_energy;   //!
   TBranch        *b_minitree_ll_mass;   //!
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
   TBranch        *b_minitree_Evts_MVAval;   //!
   TBranch        *b_minitree_Evts_MVAvalDY;   //!
   TBranch        *b_minitree_Evts_MVAvalTT;   //!
   TBranch        *b_minitree_nLLP;   //!
   TBranch        *b_minitree_LLP;   //!
   TBranch        *b_minitree_LLP_pt;   //!
   TBranch        *b_minitree_LLP_eta;   //!
   TBranch        *b_minitree_LLP_phi;   //!
   TBranch        *b_minitree_LLP_x;   //!
   TBranch        *b_minitree_LLP_y;   //!
   TBranch        *b_minitree_LLP_z;   //!
   TBranch        *b_minitree_LLP_r;   //!
   TBranch        *b_minitree_LLP_dist;   //!
   TBranch        *b_minitree_LLP_nTrks;   //!
   TBranch        *b_minitree_LLP12_dR;   //!
   TBranch        *b_minitree_LLP12_deta;   //!
   TBranch        *b_minitree_LLP12_dphi;   //!
   TBranch        *b_minitree_LLP_Mass;   //!
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
   TBranch        *b_minitree_Hemi_LLP;   //!
   TBranch        *b_minitree_Hemi_LLP_pt;   //!
   TBranch        *b_minitree_Hemi_LLP_eta;   //!
   TBranch        *b_minitree_Hemi_LLP_phi;   //!
   TBranch        *b_minitree_Hemi_LLP_dist;   //!
   TBranch        *b_minitree_Hemi_LLP_x;   //!
   TBranch        *b_minitree_Hemi_LLP_y;   //!
   TBranch        *b_minitree_Hemi_LLP_z;   //!
   TBranch        *b_minitree_Hemi_LLP_dR;   //!
   TBranch        *b_minitree_Hemi_LLP_mother;   //!
   TBranch        *b_minitree_Hemi_LLP_Vtx_dx;   //!
   TBranch        *b_minitree_Hemi_LLP_Vtx_dy;   //!
   TBranch        *b_minitree_Hemi_LLP_Vtx_dz;   //!
   TBranch        *b_minitree_Hemi_LLP_Vtx_dr;   //!
   TBranch        *b_minitree_Hemi_LLP_muOK_dR;   //!
   TBranch        *b_minitree_Hemi_LLP_muOK_pt;   //!
   TBranch        *b_minitree_Hemi_LLP_muOK_mass;   //!
   TBranch        *b_minitree_Hemi_LLP_muNO_dR;   //!
   TBranch        *b_minitree_Hemi_LLP_muNO_pt;   //!
   TBranch        *b_minitree_Hemi_LLP_muNO_mass;   //!
   TBranch        *b_minitree_Hemi_LLP_dR12;   //!
   TBranch        *b_minitree_Hemi_LLP_ping;   //!
   TBranch        *b_minitree_event_LLP_ping;   //!
   TBranch        *b_minitree_Hemi_Vtx_step;   //!
   TBranch        *b_minitree_Hemi_Vtx_isTight;   //!
   TBranch        *b_minitree_Hemi_Vtx_NChi2;   //!
   TBranch        *b_minitree_Hemi_Vtx_nTrks;   //!
   TBranch        *b_minitree_Hemi_Vtx_nTrks_sig;   //!
   TBranch        *b_minitree_Hemi_Vtx_nTrks_bad;   //!
   TBranch        *b_minitree_Hemi_Vtx_x;   //!
   TBranch        *b_minitree_Hemi_Vtx_y;   //!
   TBranch        *b_minitree_Hemi_Vtx_z;   //!
   TBranch        *b_minitree_Hemi_Vtx_r;   //!
   TBranch        *b_minitree_Hemi_Vtx_dR;   //!
   TBranch        *b_minitree_Hemi_Vtx_SumtrackWeight;   //!
   TBranch        *b_minitree_Hemi_Vtx_track_MeanDCA_d;   //!
   TBranch        *b_minitree_Hemi_Vtx_Mass;   //!
   TBranch        *b_minitree_Hemi_Vtx_dist;   //!
   TBranch        *b_minitree_event_nVtx;   //!
   TBranch        *b_minitree_event_Vtx_Vtx_dr;   //!
   TBranch        *b_minitree_event_Vtx_Vtx_dz;   //!
   TBranch        *b_minitree_event_Vtx_Vtx_dd;   //!
   TBranch        *b_minitree_event_Vtx_Vtx_reldd;   //!
   TBranch        *b_minitree_event_Vtx_Vtx_dR;   //!
   TBranch        *b_minitree_event_Vtx_Vtx_step;   //!
   TBranch        *b_minitree_Hemi_SecLLP;   //!
   TBranch        *b_minitree_Hemi_LLP_SecVtx_dz;   //!
   TBranch        *b_minitree_Hemi_LLP_SecVtx_dr;   //!
   TBranch        *b_minitree_Hemi_SecLLP_ping;   //!
   TBranch        *b_minitree_event_SecLLP_ping;   //!
   TBranch        *b_minitree_Hemi_SecVtx;   //!
   TBranch        *b_minitree_Hemi_SecVtx_step;   //!
   TBranch        *b_minitree_Hemi_SecVtx_x;   //!
   TBranch        *b_minitree_Hemi_SecVtx_y;   //!
   TBranch        *b_minitree_Hemi_SecVtx_z;   //!
   TBranch        *b_minitree_Hemi_SecVtx_r;   //!
   TBranch        *b_minitree_Hemi_SecVtx_dR;   //!
   TBranch        *b_minitree_Hemi_SecVtx_nTrks;   //!
   TBranch        *b_minitree_Hemi_SecVtx_NChi2;   //!
   TBranch        *b_minitree_Hemi_SecVtx_dist;   //!
   TBranch        *b_minitree_Hemi_SecVtx_track_MeanDCA_d;   //!
   TBranch        *b_minitree_Hemi_SecVtx_SumtrackWeight;   //!
   TBranch        *b_minitree_Hemi_SecVtx_Mass;   //!
   TBranch        *b_minitree_event_MergedVtx_Vtx_dr;   //!
   TBranch        *b_minitree_event_MergedVtx_Vtx_dz;   //!
   TBranch        *b_minitree_event_MergedVtx_Vtx_dd;   //!
   TBranch        *b_minitree_event_MergedVtx_Vtx_reldd;   //!
   TBranch        *b_minitree_event_MergedVtx_Vtx_dR;   //!
   TBranch        *b_minitree_event_MergedVtx_Vtx_step;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_nTrks;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_NChi2;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_step;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_STW;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_Mass;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_HMass;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_ntrk10;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_ntrk20;   //!
   TBranch        *b_minitree_Hemi_Vtx_BDT_MeanDCA;   //!
   TBranch        *b_minitree_Hemi_Vtx_MVAval_Loose;   //!
   TBranch        *b_minitree_Hemi_Vtx_MVAval_Tight;   //!
   TBranch        *b_miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;   //!
   TBranch        *b_miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v;   //!
   TBranch        *b_miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v;   //!
   TBranch        *b_miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_miniHLT_Ele27_WPTight_Gsf_v;   //!
   TBranch        *b_miniHLT_Ele32_WPTight_Gsf_v;   //!
   TBranch        *b_miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;   //!
   TBranch        *b_miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;   //!
   TBranch        *b_miniHLT_PFMET120_PFMHT120_IDTight_v;   //!
   TBranch        *b_miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v;   //!
   TBranch        *b_miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;   //!
   TBranch        *b_miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;   //!
   TBranch        *b_miniHLT_PFMET250_HBHECleaned_v;   //!
   TBranch        *b_miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;   //!
   TBranch        *b_miniHLT_IsoMu24_v;   //!
   TBranch        *b_miniHLT_IsoMu27_v;   //!

   ABCD(TTree *tree=0);
   virtual ~ABCD();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString sample ="", TString Production="",bool Signal = true);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ABCD_cxx
ABCD::ABCD(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MiniDYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MiniDYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root");
      }
      f->GetObject("ttree",tree);

   }
   Init(tree);
}

ABCD::~ABCD()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ABCD::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ABCD::LoadTree(Long64_t entry)
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

void ABCD::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   minirunNumber = 0;
   minieventNumber = 0;
   minilumiBlock = 0;
   minitree_LHE_Weights = 0;
   minitree_MCEvt_weight = 0;
   minitree_only_gen_wt = 0;
   minitree_event_weight = 0;
   minitree_genTop_Weight = 0;
   miniPUweight = 0;
   miniPrefweight = 0;
   miniPU_events = 0;
   miniAllPU_events_weight = 0;
   minitree_Filter = 0;
   minitree_FilterSameSign = 0;
   minitree_GenPVx = 0;
   minitree_GenPVy = 0;
   minitree_GenPVz = 0;
   minitree_smu_mass = 0;
   minitree_neu_mass = 0;
   minitree_neu_ctau = 0;
   minitree_Good_PV = 0;
   minitree_nPV = 0;
   minitree_PV_x = 0;
   minitree_PV_y = 0;
   minitree_PV_z = 0;
   minitree_PV_ez = 0;
   minitree_PV_NChi2 = 0;
   minitree_PV_ndf = 0;
   minitree_PFMet_et = 0;
   minitree_PFMet_phi = 0;
   minitree_HT = 0;
   minitree_TRACK_SIZE = 0;
   minitree_nTracks = 0;
   minitree_nLostTracks = 0;
   minitree_muon_GenRecoTriggerMatched = 0;
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
   minitree_muon_gen = 0;
   minitree_lepton_leadingpt = 0;
   minitree_lepton_leadingpt2 = 0;
   minitree_lepton_lepton_dR = 0;
   minitree_lepton_lepton_dPhi = 0;
   minitree_lepton_lepton_dEta = 0;
   minitree_ll_pt = 0;
   minitree_ll_eta = 0;
   minitree_ll_phi = 0;
   minitree_ll_px = 0;
   minitree_ll_py = 0;
   minitree_ll_pz = 0;
   minitree_ll_energy = 0;
   minitree_ll_mass = 0;
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
   minitree_Evts_MVAval = 0;
   minitree_Evts_MVAvalDY = 0;
   minitree_Evts_MVAvalTT = 0;
   minitree_nLLP = 0;
   minitree_LLP = 0;
   minitree_LLP_pt = 0;
   minitree_LLP_eta = 0;
   minitree_LLP_phi = 0;
   minitree_LLP_x = 0;
   minitree_LLP_y = 0;
   minitree_LLP_z = 0;
   minitree_LLP_r = 0;
   minitree_LLP_dist = 0;
   minitree_LLP_nTrks = 0;
   minitree_LLP12_dR = 0;
   minitree_LLP12_deta = 0;
   minitree_LLP12_dphi = 0;
   minitree_LLP_Mass = 0;
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
   minitree_Hemi_LLP = 0;
   minitree_Hemi_LLP_pt = 0;
   minitree_Hemi_LLP_eta = 0;
   minitree_Hemi_LLP_phi = 0;
   minitree_Hemi_LLP_dist = 0;
   minitree_Hemi_LLP_x = 0;
   minitree_Hemi_LLP_y = 0;
   minitree_Hemi_LLP_z = 0;
   minitree_Hemi_LLP_dR = 0;
   minitree_Hemi_LLP_mother = 0;
   minitree_Hemi_LLP_Vtx_dx = 0;
   minitree_Hemi_LLP_Vtx_dy = 0;
   minitree_Hemi_LLP_Vtx_dz = 0;
   minitree_Hemi_LLP_Vtx_dr = 0;
   minitree_Hemi_LLP_muOK_dR = 0;
   minitree_Hemi_LLP_muOK_pt = 0;
   minitree_Hemi_LLP_muOK_mass = 0;
   minitree_Hemi_LLP_muNO_dR = 0;
   minitree_Hemi_LLP_muNO_pt = 0;
   minitree_Hemi_LLP_muNO_mass = 0;
   minitree_Hemi_LLP_dR12 = 0;
   minitree_Hemi_LLP_ping = 0;
   minitree_event_LLP_ping = 0;
   minitree_Hemi_Vtx_step = 0;
   minitree_Hemi_Vtx_isTight = 0;
   minitree_Hemi_Vtx_NChi2 = 0;
   minitree_Hemi_Vtx_nTrks = 0;
   minitree_Hemi_Vtx_nTrks_sig = 0;
   minitree_Hemi_Vtx_nTrks_bad = 0;
   minitree_Hemi_Vtx_x = 0;
   minitree_Hemi_Vtx_y = 0;
   minitree_Hemi_Vtx_z = 0;
   minitree_Hemi_Vtx_r = 0;
   minitree_Hemi_Vtx_dR = 0;
   minitree_Hemi_Vtx_SumtrackWeight = 0;
   minitree_Hemi_Vtx_track_MeanDCA_d = 0;
   minitree_Hemi_Vtx_Mass = 0;
   minitree_Hemi_Vtx_dist = 0;
   minitree_event_nVtx = 0;
   minitree_event_Vtx_Vtx_dr = 0;
   minitree_event_Vtx_Vtx_dz = 0;
   minitree_event_Vtx_Vtx_dd = 0;
   minitree_event_Vtx_Vtx_reldd = 0;
   minitree_event_Vtx_Vtx_dR = 0;
   minitree_event_Vtx_Vtx_step = 0;
   minitree_Hemi_SecLLP = 0;
   minitree_Hemi_LLP_SecVtx_dz = 0;
   minitree_Hemi_LLP_SecVtx_dr = 0;
   minitree_Hemi_SecLLP_ping = 0;
   minitree_event_SecLLP_ping = 0;
   minitree_Hemi_SecVtx = 0;
   minitree_Hemi_SecVtx_step = 0;
   minitree_Hemi_SecVtx_x = 0;
   minitree_Hemi_SecVtx_y = 0;
   minitree_Hemi_SecVtx_z = 0;
   minitree_Hemi_SecVtx_r = 0;
   minitree_Hemi_SecVtx_dR = 0;
   minitree_Hemi_SecVtx_nTrks = 0;
   minitree_Hemi_SecVtx_NChi2 = 0;
   minitree_Hemi_SecVtx_dist = 0;
   minitree_Hemi_SecVtx_track_MeanDCA_d = 0;
   minitree_Hemi_SecVtx_SumtrackWeight = 0;
   minitree_Hemi_SecVtx_Mass = 0;
   minitree_event_MergedVtx_Vtx_dr = 0;
   minitree_event_MergedVtx_Vtx_dz = 0;
   minitree_event_MergedVtx_Vtx_dd = 0;
   minitree_event_MergedVtx_Vtx_reldd = 0;
   minitree_event_MergedVtx_Vtx_dR = 0;
   minitree_event_MergedVtx_Vtx_step = 0;
   minitree_Hemi_Vtx_BDT_nTrks = 0;
   minitree_Hemi_Vtx_BDT_NChi2 = 0;
   minitree_Hemi_Vtx_BDT_step = 0;
   minitree_Hemi_Vtx_BDT_STW = 0;
   minitree_Hemi_Vtx_BDT_Mass = 0;
   minitree_Hemi_Vtx_BDT_HMass = 0;
   minitree_Hemi_Vtx_BDT_ntrk10 = 0;
   minitree_Hemi_Vtx_BDT_ntrk20 = 0;
   minitree_Hemi_Vtx_BDT_MeanDCA = 0;
   minitree_Hemi_Vtx_MVAval_Loose = 0;
   minitree_Hemi_Vtx_MVAval_Tight = 0;
   miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v = 0;
   miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v = 0;
   miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v = 0;
   miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = 0;
   miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = 0;
   miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = 0;
   miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = 0;
   miniHLT_Ele27_WPTight_Gsf_v = 0;
   miniHLT_Ele32_WPTight_Gsf_v = 0;
   miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = 0;
   miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v = 0;
   miniHLT_PFMET120_PFMHT120_IDTight_v = 0;
   miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v = 0;
   miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v = 0;
   miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v = 0;
   miniHLT_PFMET250_HBHECleaned_v = 0;
   miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v = 0;
   miniHLT_IsoMu24_v = 0;
   miniHLT_IsoMu27_v = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("minirunNumber", &minirunNumber, &b_minirunNumber);
   fChain->SetBranchAddress("minieventNumber", &minieventNumber, &b_minieventNumber);
   fChain->SetBranchAddress("minilumiBlock", &minilumiBlock, &b_minilumiBlock);
   fChain->SetBranchAddress("minitree_LHE_Weights", &minitree_LHE_Weights, &b_minitree_LHE_Weights);
   fChain->SetBranchAddress("minitree_MCEvt_weight", &minitree_MCEvt_weight, &b_minitree_MCEvt_weight);
   fChain->SetBranchAddress("minitree_only_gen_wt", &minitree_only_gen_wt, &b_minitree_only_gen_wt);
   fChain->SetBranchAddress("minitree_event_weight", &minitree_event_weight, &b_minitree_event_weight);
   fChain->SetBranchAddress("minitree_genTop_Weight", &minitree_genTop_Weight, &b_minitree_genTop_Weight);
   fChain->SetBranchAddress("miniPUweight", &miniPUweight, &b_miniPUweight);
   fChain->SetBranchAddress("miniPrefweight", &miniPrefweight, &b_miniPrefweight);
   fChain->SetBranchAddress("miniPU_events", &miniPU_events, &b_miniPU_events);
   fChain->SetBranchAddress("miniAllPU_events_weight", &miniAllPU_events_weight, &b_miniAllPU_events_weight);
   fChain->SetBranchAddress("minitree_Filter", &minitree_Filter, &b_minitree_Filter);
   fChain->SetBranchAddress("minitree_FilterSameSign", &minitree_FilterSameSign, &b_minitree_FilterSameSign);
   fChain->SetBranchAddress("minitree_GenPVx", &minitree_GenPVx, &b_minitree_GenPVx);
   fChain->SetBranchAddress("minitree_GenPVy", &minitree_GenPVy, &b_minitree_GenPVy);
   fChain->SetBranchAddress("minitree_GenPVz", &minitree_GenPVz, &b_minitree_GenPVz);
   fChain->SetBranchAddress("minitree_smu_mass", &minitree_smu_mass, &b_minitree_smu_mass);
   fChain->SetBranchAddress("minitree_neu_mass", &minitree_neu_mass, &b_minitree_neu_mass);
   fChain->SetBranchAddress("minitree_neu_ctau", &minitree_neu_ctau, &b_minitree_neu_ctau);
   fChain->SetBranchAddress("minitree_Good_PV", &minitree_Good_PV, &b_minitree_Good_PV);
   fChain->SetBranchAddress("minitree_nPV", &minitree_nPV, &b_minitree_nPV);
   fChain->SetBranchAddress("minitree_PV_x", &minitree_PV_x, &b_minitree_PV_x);
   fChain->SetBranchAddress("minitree_PV_y", &minitree_PV_y, &b_minitree_PV_y);
   fChain->SetBranchAddress("minitree_PV_z", &minitree_PV_z, &b_minitree_PV_z);
   fChain->SetBranchAddress("minitree_PV_ez", &minitree_PV_ez, &b_minitree_PV_ez);
   fChain->SetBranchAddress("minitree_PV_NChi2", &minitree_PV_NChi2, &b_minitree_PV_NChi2);
   fChain->SetBranchAddress("minitree_PV_ndf", &minitree_PV_ndf, &b_minitree_PV_ndf);
   fChain->SetBranchAddress("minitree_PFMet_et", &minitree_PFMet_et, &b_minitree_PFMet_et);
   fChain->SetBranchAddress("minitree_PFMet_phi", &minitree_PFMet_phi, &b_minitree_PFMet_phi);
   fChain->SetBranchAddress("minitree_HT", &minitree_HT, &b_minitree_HT);
   fChain->SetBranchAddress("minitree_TRACK_SIZE", &minitree_TRACK_SIZE, &b_minitree_TRACK_SIZE);
   fChain->SetBranchAddress("minitree_nTracks", &minitree_nTracks, &b_minitree_nTracks);
   fChain->SetBranchAddress("minitree_nLostTracks", &minitree_nLostTracks, &b_minitree_nLostTracks);
   fChain->SetBranchAddress("minitree_muon_GenRecoTriggerMatched", &minitree_muon_GenRecoTriggerMatched, &b_minitree_muon_GenRecoTriggerMatched);
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
   fChain->SetBranchAddress("minitree_muon_gen", &minitree_muon_gen, &b_minitree_muon_gen);
   fChain->SetBranchAddress("minitree_lepton_leadingpt", &minitree_lepton_leadingpt, &b_minitree_lepton_leadingpt);
   fChain->SetBranchAddress("minitree_lepton_leadingpt2", &minitree_lepton_leadingpt2, &b_minitree_lepton_leadingpt2);
   fChain->SetBranchAddress("minitree_lepton_lepton_dR", &minitree_lepton_lepton_dR, &b_minitree_lepton_lepton_dR);
   fChain->SetBranchAddress("minitree_lepton_lepton_dPhi", &minitree_lepton_lepton_dPhi, &b_minitree_lepton_lepton_dPhi);
   fChain->SetBranchAddress("minitree_lepton_lepton_dEta", &minitree_lepton_lepton_dEta, &b_minitree_lepton_lepton_dEta);
   fChain->SetBranchAddress("minitree_ll_pt", &minitree_ll_pt, &b_minitree_ll_pt);
   fChain->SetBranchAddress("minitree_ll_eta", &minitree_ll_eta, &b_minitree_ll_eta);
   fChain->SetBranchAddress("minitree_ll_phi", &minitree_ll_phi, &b_minitree_ll_phi);
   fChain->SetBranchAddress("minitree_ll_px", &minitree_ll_px, &b_minitree_ll_px);
   fChain->SetBranchAddress("minitree_ll_py", &minitree_ll_py, &b_minitree_ll_py);
   fChain->SetBranchAddress("minitree_ll_pz", &minitree_ll_pz, &b_minitree_ll_pz);
   fChain->SetBranchAddress("minitree_ll_energy", &minitree_ll_energy, &b_minitree_ll_energy);
   fChain->SetBranchAddress("minitree_ll_mass", &minitree_ll_mass, &b_minitree_ll_mass);
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
   fChain->SetBranchAddress("minitree_Evts_MVAval", &minitree_Evts_MVAval, &b_minitree_Evts_MVAval);
   fChain->SetBranchAddress("minitree_Evts_MVAvalDY", &minitree_Evts_MVAvalDY, &b_minitree_Evts_MVAvalDY);
   fChain->SetBranchAddress("minitree_Evts_MVAvalTT", &minitree_Evts_MVAvalTT, &b_minitree_Evts_MVAvalTT);
   fChain->SetBranchAddress("minitree_nLLP", &minitree_nLLP, &b_minitree_nLLP);
   fChain->SetBranchAddress("minitree_LLP", &minitree_LLP, &b_minitree_LLP);
   fChain->SetBranchAddress("minitree_LLP_pt", &minitree_LLP_pt, &b_minitree_LLP_pt);
   fChain->SetBranchAddress("minitree_LLP_eta", &minitree_LLP_eta, &b_minitree_LLP_eta);
   fChain->SetBranchAddress("minitree_LLP_phi", &minitree_LLP_phi, &b_minitree_LLP_phi);
   fChain->SetBranchAddress("minitree_LLP_x", &minitree_LLP_x, &b_minitree_LLP_x);
   fChain->SetBranchAddress("minitree_LLP_y", &minitree_LLP_y, &b_minitree_LLP_y);
   fChain->SetBranchAddress("minitree_LLP_z", &minitree_LLP_z, &b_minitree_LLP_z);
   fChain->SetBranchAddress("minitree_LLP_r", &minitree_LLP_r, &b_minitree_LLP_r);
   fChain->SetBranchAddress("minitree_LLP_dist", &minitree_LLP_dist, &b_minitree_LLP_dist);
   fChain->SetBranchAddress("minitree_LLP_nTrks", &minitree_LLP_nTrks, &b_minitree_LLP_nTrks);
   fChain->SetBranchAddress("minitree_LLP12_dR", &minitree_LLP12_dR, &b_minitree_LLP12_dR);
   fChain->SetBranchAddress("minitree_LLP12_deta", &minitree_LLP12_deta, &b_minitree_LLP12_deta);
   fChain->SetBranchAddress("minitree_LLP12_dphi", &minitree_LLP12_dphi, &b_minitree_LLP12_dphi);
   fChain->SetBranchAddress("minitree_LLP_Mass", &minitree_LLP_Mass, &b_minitree_LLP_Mass);
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
   fChain->SetBranchAddress("minitree_Hemi_LLP", &minitree_Hemi_LLP, &b_minitree_Hemi_LLP);
   fChain->SetBranchAddress("minitree_Hemi_LLP_pt", &minitree_Hemi_LLP_pt, &b_minitree_Hemi_LLP_pt);
   fChain->SetBranchAddress("minitree_Hemi_LLP_eta", &minitree_Hemi_LLP_eta, &b_minitree_Hemi_LLP_eta);
   fChain->SetBranchAddress("minitree_Hemi_LLP_phi", &minitree_Hemi_LLP_phi, &b_minitree_Hemi_LLP_phi);
   fChain->SetBranchAddress("minitree_Hemi_LLP_dist", &minitree_Hemi_LLP_dist, &b_minitree_Hemi_LLP_dist);
   fChain->SetBranchAddress("minitree_Hemi_LLP_x", &minitree_Hemi_LLP_x, &b_minitree_Hemi_LLP_x);
   fChain->SetBranchAddress("minitree_Hemi_LLP_y", &minitree_Hemi_LLP_y, &b_minitree_Hemi_LLP_y);
   fChain->SetBranchAddress("minitree_Hemi_LLP_z", &minitree_Hemi_LLP_z, &b_minitree_Hemi_LLP_z);
   fChain->SetBranchAddress("minitree_Hemi_LLP_dR", &minitree_Hemi_LLP_dR, &b_minitree_Hemi_LLP_dR);
   fChain->SetBranchAddress("minitree_Hemi_LLP_mother", &minitree_Hemi_LLP_mother, &b_minitree_Hemi_LLP_mother);
   fChain->SetBranchAddress("minitree_Hemi_LLP_Vtx_dx", &minitree_Hemi_LLP_Vtx_dx, &b_minitree_Hemi_LLP_Vtx_dx);
   fChain->SetBranchAddress("minitree_Hemi_LLP_Vtx_dy", &minitree_Hemi_LLP_Vtx_dy, &b_minitree_Hemi_LLP_Vtx_dy);
   fChain->SetBranchAddress("minitree_Hemi_LLP_Vtx_dz", &minitree_Hemi_LLP_Vtx_dz, &b_minitree_Hemi_LLP_Vtx_dz);
   fChain->SetBranchAddress("minitree_Hemi_LLP_Vtx_dr", &minitree_Hemi_LLP_Vtx_dr, &b_minitree_Hemi_LLP_Vtx_dr);
   fChain->SetBranchAddress("minitree_Hemi_LLP_muOK_dR", &minitree_Hemi_LLP_muOK_dR, &b_minitree_Hemi_LLP_muOK_dR);
   fChain->SetBranchAddress("minitree_Hemi_LLP_muOK_pt", &minitree_Hemi_LLP_muOK_pt, &b_minitree_Hemi_LLP_muOK_pt);
   fChain->SetBranchAddress("minitree_Hemi_LLP_muOK_mass", &minitree_Hemi_LLP_muOK_mass, &b_minitree_Hemi_LLP_muOK_mass);
   fChain->SetBranchAddress("minitree_Hemi_LLP_muNO_dR", &minitree_Hemi_LLP_muNO_dR, &b_minitree_Hemi_LLP_muNO_dR);
   fChain->SetBranchAddress("minitree_Hemi_LLP_muNO_pt", &minitree_Hemi_LLP_muNO_pt, &b_minitree_Hemi_LLP_muNO_pt);
   fChain->SetBranchAddress("minitree_Hemi_LLP_muNO_mass", &minitree_Hemi_LLP_muNO_mass, &b_minitree_Hemi_LLP_muNO_mass);
   fChain->SetBranchAddress("minitree_Hemi_LLP_dR12", &minitree_Hemi_LLP_dR12, &b_minitree_Hemi_LLP_dR12);
   fChain->SetBranchAddress("minitree_Hemi_LLP_ping", &minitree_Hemi_LLP_ping, &b_minitree_Hemi_LLP_ping);
   fChain->SetBranchAddress("minitree_event_LLP_ping", &minitree_event_LLP_ping, &b_minitree_event_LLP_ping);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_step", &minitree_Hemi_Vtx_step, &b_minitree_Hemi_Vtx_step);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_isTight", &minitree_Hemi_Vtx_isTight, &b_minitree_Hemi_Vtx_isTight);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_NChi2", &minitree_Hemi_Vtx_NChi2, &b_minitree_Hemi_Vtx_NChi2);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_nTrks", &minitree_Hemi_Vtx_nTrks, &b_minitree_Hemi_Vtx_nTrks);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_nTrks_sig", &minitree_Hemi_Vtx_nTrks_sig, &b_minitree_Hemi_Vtx_nTrks_sig);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_nTrks_bad", &minitree_Hemi_Vtx_nTrks_bad, &b_minitree_Hemi_Vtx_nTrks_bad);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_x", &minitree_Hemi_Vtx_x, &b_minitree_Hemi_Vtx_x);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_y", &minitree_Hemi_Vtx_y, &b_minitree_Hemi_Vtx_y);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_z", &minitree_Hemi_Vtx_z, &b_minitree_Hemi_Vtx_z);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_r", &minitree_Hemi_Vtx_r, &b_minitree_Hemi_Vtx_r);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_dR", &minitree_Hemi_Vtx_dR, &b_minitree_Hemi_Vtx_dR);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_SumtrackWeight", &minitree_Hemi_Vtx_SumtrackWeight, &b_minitree_Hemi_Vtx_SumtrackWeight);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_track_MeanDCA_d", &minitree_Hemi_Vtx_track_MeanDCA_d, &b_minitree_Hemi_Vtx_track_MeanDCA_d);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_Mass", &minitree_Hemi_Vtx_Mass, &b_minitree_Hemi_Vtx_Mass);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_dist", &minitree_Hemi_Vtx_dist, &b_minitree_Hemi_Vtx_dist);
   fChain->SetBranchAddress("minitree_event_nVtx", &minitree_event_nVtx, &b_minitree_event_nVtx);
   fChain->SetBranchAddress("minitree_event_Vtx_Vtx_dr", &minitree_event_Vtx_Vtx_dr, &b_minitree_event_Vtx_Vtx_dr);
   fChain->SetBranchAddress("minitree_event_Vtx_Vtx_dz", &minitree_event_Vtx_Vtx_dz, &b_minitree_event_Vtx_Vtx_dz);
   fChain->SetBranchAddress("minitree_event_Vtx_Vtx_dd", &minitree_event_Vtx_Vtx_dd, &b_minitree_event_Vtx_Vtx_dd);
   fChain->SetBranchAddress("minitree_event_Vtx_Vtx_reldd", &minitree_event_Vtx_Vtx_reldd, &b_minitree_event_Vtx_Vtx_reldd);
   fChain->SetBranchAddress("minitree_event_Vtx_Vtx_dR", &minitree_event_Vtx_Vtx_dR, &b_minitree_event_Vtx_Vtx_dR);
   fChain->SetBranchAddress("minitree_event_Vtx_Vtx_step", &minitree_event_Vtx_Vtx_step, &b_minitree_event_Vtx_Vtx_step);
   fChain->SetBranchAddress("minitree_Hemi_SecLLP", &minitree_Hemi_SecLLP, &b_minitree_Hemi_SecLLP);
   fChain->SetBranchAddress("minitree_Hemi_LLP_SecVtx_dz", &minitree_Hemi_LLP_SecVtx_dz, &b_minitree_Hemi_LLP_SecVtx_dz);
   fChain->SetBranchAddress("minitree_Hemi_LLP_SecVtx_dr", &minitree_Hemi_LLP_SecVtx_dr, &b_minitree_Hemi_LLP_SecVtx_dr);
   fChain->SetBranchAddress("minitree_Hemi_SecLLP_ping", &minitree_Hemi_SecLLP_ping, &b_minitree_Hemi_SecLLP_ping);
   fChain->SetBranchAddress("minitree_event_SecLLP_ping", &minitree_event_SecLLP_ping, &b_minitree_event_SecLLP_ping);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx", &minitree_Hemi_SecVtx, &b_minitree_Hemi_SecVtx);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_step", &minitree_Hemi_SecVtx_step, &b_minitree_Hemi_SecVtx_step);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_x", &minitree_Hemi_SecVtx_x, &b_minitree_Hemi_SecVtx_x);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_y", &minitree_Hemi_SecVtx_y, &b_minitree_Hemi_SecVtx_y);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_z", &minitree_Hemi_SecVtx_z, &b_minitree_Hemi_SecVtx_z);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_r", &minitree_Hemi_SecVtx_r, &b_minitree_Hemi_SecVtx_r);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_dR", &minitree_Hemi_SecVtx_dR, &b_minitree_Hemi_SecVtx_dR);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_nTrks", &minitree_Hemi_SecVtx_nTrks, &b_minitree_Hemi_SecVtx_nTrks);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_NChi2", &minitree_Hemi_SecVtx_NChi2, &b_minitree_Hemi_SecVtx_NChi2);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_dist", &minitree_Hemi_SecVtx_dist, &b_minitree_Hemi_SecVtx_dist);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_track_MeanDCA_d", &minitree_Hemi_SecVtx_track_MeanDCA_d, &b_minitree_Hemi_SecVtx_track_MeanDCA_d);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_SumtrackWeight", &minitree_Hemi_SecVtx_SumtrackWeight, &b_minitree_Hemi_SecVtx_SumtrackWeight);
   fChain->SetBranchAddress("minitree_Hemi_SecVtx_Mass", &minitree_Hemi_SecVtx_Mass, &b_minitree_Hemi_SecVtx_Mass);
   fChain->SetBranchAddress("minitree_event_MergedVtx_Vtx_dr", &minitree_event_MergedVtx_Vtx_dr, &b_minitree_event_MergedVtx_Vtx_dr);
   fChain->SetBranchAddress("minitree_event_MergedVtx_Vtx_dz", &minitree_event_MergedVtx_Vtx_dz, &b_minitree_event_MergedVtx_Vtx_dz);
   fChain->SetBranchAddress("minitree_event_MergedVtx_Vtx_dd", &minitree_event_MergedVtx_Vtx_dd, &b_minitree_event_MergedVtx_Vtx_dd);
   fChain->SetBranchAddress("minitree_event_MergedVtx_Vtx_reldd", &minitree_event_MergedVtx_Vtx_reldd, &b_minitree_event_MergedVtx_Vtx_reldd);
   fChain->SetBranchAddress("minitree_event_MergedVtx_Vtx_dR", &minitree_event_MergedVtx_Vtx_dR, &b_minitree_event_MergedVtx_Vtx_dR);
   fChain->SetBranchAddress("minitree_event_MergedVtx_Vtx_step", &minitree_event_MergedVtx_Vtx_step, &b_minitree_event_MergedVtx_Vtx_step);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_nTrks", &minitree_Hemi_Vtx_BDT_nTrks, &b_minitree_Hemi_Vtx_BDT_nTrks);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_NChi2", &minitree_Hemi_Vtx_BDT_NChi2, &b_minitree_Hemi_Vtx_BDT_NChi2);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_step", &minitree_Hemi_Vtx_BDT_step, &b_minitree_Hemi_Vtx_BDT_step);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_STW", &minitree_Hemi_Vtx_BDT_STW, &b_minitree_Hemi_Vtx_BDT_STW);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_Mass", &minitree_Hemi_Vtx_BDT_Mass, &b_minitree_Hemi_Vtx_BDT_Mass);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_HMass", &minitree_Hemi_Vtx_BDT_HMass, &b_minitree_Hemi_Vtx_BDT_HMass);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_ntrk10", &minitree_Hemi_Vtx_BDT_ntrk10, &b_minitree_Hemi_Vtx_BDT_ntrk10);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_ntrk20", &minitree_Hemi_Vtx_BDT_ntrk20, &b_minitree_Hemi_Vtx_BDT_ntrk20);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_BDT_MeanDCA", &minitree_Hemi_Vtx_BDT_MeanDCA, &b_minitree_Hemi_Vtx_BDT_MeanDCA);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_MVAval_Loose", &minitree_Hemi_Vtx_MVAval_Loose, &b_minitree_Hemi_Vtx_MVAval_Loose);
   fChain->SetBranchAddress("minitree_Hemi_Vtx_MVAval_Tight", &minitree_Hemi_Vtx_MVAval_Tight, &b_minitree_Hemi_Vtx_MVAval_Tight);
   fChain->SetBranchAddress("miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v", &miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v, &b_miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
   fChain->SetBranchAddress("miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v", &miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v, &b_miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v);
   fChain->SetBranchAddress("miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v", &miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v, &b_miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v);
   fChain->SetBranchAddress("miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", &miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v", &miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v, &b_miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v", &miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v", &miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("miniHLT_Ele27_WPTight_Gsf_v", &miniHLT_Ele27_WPTight_Gsf_v, &b_miniHLT_Ele27_WPTight_Gsf_v);
   fChain->SetBranchAddress("miniHLT_Ele32_WPTight_Gsf_v", &miniHLT_Ele32_WPTight_Gsf_v, &b_miniHLT_Ele32_WPTight_Gsf_v);
   fChain->SetBranchAddress("miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v, &b_miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
   fChain->SetBranchAddress("miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v", &miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v, &b_miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   fChain->SetBranchAddress("miniHLT_PFMET120_PFMHT120_IDTight_v", &miniHLT_PFMET120_PFMHT120_IDTight_v, &b_miniHLT_PFMET120_PFMHT120_IDTight_v);
   fChain->SetBranchAddress("miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v", &miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v, &b_miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v);
   fChain->SetBranchAddress("miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v", &miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v, &b_miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v);
   fChain->SetBranchAddress("miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v", &miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v, &b_miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
   fChain->SetBranchAddress("miniHLT_PFMET250_HBHECleaned_v", &miniHLT_PFMET250_HBHECleaned_v, &b_miniHLT_PFMET250_HBHECleaned_v);
   fChain->SetBranchAddress("miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v", &miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v, &b_miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v);
   fChain->SetBranchAddress("miniHLT_IsoMu24_v", &miniHLT_IsoMu24_v, &b_miniHLT_IsoMu24_v);
   fChain->SetBranchAddress("miniHLT_IsoMu27_v", &miniHLT_IsoMu27_v, &b_miniHLT_IsoMu27_v);
   Notify();
}

Bool_t ABCD::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ABCD::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ABCD::Cut(Long64_t entry)
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
#endif // #ifdef ABCD_cxx
