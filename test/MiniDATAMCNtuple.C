#define MiniDATAMCNtuple_cxx
#include "MiniDATAMCNtuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "HistogramManager.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "TTree.h"
using namespace std;



void MiniDATAMCNtuple::Loop(TString sample , TString Production,bool Signal )
{
//   In a ROOT session, you can do:
//      root> .L MiniDATAMCNtuple.C
//      root> MiniDATAMCNtuple t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   TFile * myFile = new TFile( (Production+"/MiniDATAMC_"+sample+".root").Data(), "recreate");
   TTree *smalltree = new TTree("ttree", "summary information");
  

   std::vector<double> minitree_only_gen_wt;
   std::vector<double> minitree_genTop_Weight;

   std::vector<double> minitree_gen_top_pt;
   std::vector<double> minitree_gen_top_rw_pt;

   std::vector<double> miniPUweight;
   std::vector<double> miniPUweight_Up;
   std::vector<double> miniPUweight_Down;
   std::vector<double> miniPrefweight;

   std::vector<bool>   minitree_Filter;
   std::vector<bool>   minitree_FilterSameSign;
   std::vector<int>    minitree_nPV;
   std::vector<bool>   minitree_trigger_doublelepton;
   std::vector<bool>   minitree_trigger_singlelepton;
   std::vector<bool>   minitree_Good_PV;
   std::vector<int>    minitree_smu_mass; // used for signal only
   std::vector<int>    minitree_neu_mass; // used for signal only
   std::vector<float>  minitree_neu_ctau; // used for signal only

   std::vector<float>  minitree_Mmumu;
   std::vector<float>  minitree_MmumuSameSign;

   std::vector<float>  minitree_ll_pt;
   std::vector<float>  minitree_ll_eta;
   std::vector<float>  minitree_ll_phi;
   std::vector<float>  minitree_ll_px;
   std::vector<float>  minitree_ll_py;
   std::vector<float>  minitree_ll_pz;
   std::vector<float>  minitree_ll_energy;
   std::vector<float>  minitree_ll_mass;
   std::vector<float>  minitree_all_nmu;
   std::vector<float>  minitree_nmu;
   std::vector<bool>   minitree_muon_isPrompt;
   std::vector<float>  minitree_muon_pt;
   std::vector<float>  minitree_muon_SF;
   std::vector<float>  minitree_muon_eta;
   std::vector<float>  minitree_muon_phi;
   std::vector<float>  minitree_muon_dxy;
   std::vector<float>  minitree_muon_dz;
   std::vector<int>    minitree_muon_charge;
   std::vector<float>  minitree_muon_correction;
   std::vector<float>  minitree_muon_dxyError;
   std::vector<float>  minitree_muon_dzError;
   std::vector<bool>   minitree_muon_isLoose;
   std::vector<bool>   minitree_muon_isMedium;
   std::vector<bool>   minitree_muon_isTight;
   std::vector<bool>   minitree_muon_isGlobal;
   std::vector<bool>   minitree_muon_PFIsoVeryLoose;
   std::vector<bool>   minitree_muon_PFIsoLoose;
   std::vector<bool>   minitree_muon_PFIsoMedium;
   std::vector<bool>   minitree_muon_PFIsoTight;
   std::vector<bool>   minitree_muon_TkIsoLoose;
   std::vector<bool>   minitree_muon_TkIsoTight;
   std::vector<bool>   minitree_muon_MiniIsoLoose;
   std::vector<bool>   minitree_muon_MiniIsoMedium;
   std::vector<bool>   minitree_muon_MiniIsoTight;

   std::vector<float>  minitree_lepton_leadingpt;
   std::vector<float>  minitree_lepton_leadingpt2;
   std::vector<float>  minitree_lepton_leadingeta;
   std::vector<float>  minitree_lepton_leadingeta2;
   std::vector<float>  minitree_lepton_leadingphi;
   std::vector<float>  minitree_lepton_leadingphi2;
   std::vector<float>  minitree_lepton_lepton_dR;
   std::vector<float>  minitree_lepton_lepton_dPhi;
   std::vector<float>  minitree_lepton_lepton_dEta;
   
   std::vector<float>  minitree_lepton_leadingdxy;
   std::vector<float>  minitree_lepton_leadingdxy2;
   std::vector<float>  minitree_lepton_leadingdz;
   std::vector<float>  minitree_lepton_leadingdz2;

   std::vector<float>   minitree_reco_lepton_leadingpt;
   std::vector<float>   minitree_reco_lepton_leadingpt2;
   std::vector<float>   minitree_reco_lepton_leadingeta;
   std::vector<float>   minitree_reco_lepton_leadingeta2;
   std::vector<float>   minitree_reco_lepton_leadingphi;
   std::vector<float>   minitree_reco_lepton_leadingphi2;

   std::vector<float>  minitree_lepton_b4trigger_leadingpt;
   std::vector<float>  minitree_lepton_b4trigger_leadingpt2;

   std::vector<int>    minitree_all_nel; 
   std::vector<int>    minitree_electron_nEle; 
   std::vector<bool>   minitree_electron_isPrompt;
   std::vector<float>  minitree_electron_pt;
   std::vector<float>  minitree_electron_eta;
   std::vector<float>  minitree_electron_phi;
   std::vector<int>    minitree_electron_charge;
   std::vector<float>  minitree_electron_dxy;
   std::vector<float>  minitree_electron_dz;
   std::vector<int>    minitree_electron_gen;
   std::vector<float>  minitree_electron_energy;
   std::vector<float>  minitree_electron_et;
   std::vector<float>  minitree_electron_ecal_trk_postcorr;
   std::vector<float>  minitree_electron_isoR4;
   std::vector<bool>   minitree_electron_IsLoose;
   std::vector<bool>   minitree_electron_IsMedium;
   std::vector<bool>   minitree_electron_IsTight;

   std::vector<int>    minitree_njet;
   std::vector<int>    minitree_njetNOmu; // in BDT_EVT
   std::vector<float>  minitree_jet_pt;
   std::vector<float>  minitree_jet_eta;
   std::vector<float>  minitree_jet_phi;
   std::vector<float>  minitree_jet_HadronFlavour;
   std::vector<float>  minitree_jet_btag_DeepJet;
   std::vector<float>  minitree_jet_E;

   std::vector<float>  minitree_HT;
   std::vector<float>  minitree_LT;

   std::vector<float>  minitree_jet_leadingpt; // in BDT_EVT
   std::vector<float>  minitree_jet_leadingpt2; // in BDT_EVT
   std::vector<float>  minitree_jet_leadingeta; // in BDT_EVT
   std::vector<float>  minitree_jet_leadingeta2; // in BDT_EVT
   std::vector<float>  minitree_jet_jet_dR; // in BDT_EVT
   std::vector<float>  minitree_jet_jet_dPhi; // in BDT_EVT
   std::vector<float>  minitree_jet_jet_dEta; // in BDT_EVT


   std::vector<float>  minitree_muon_jet_dRmin;
   std::vector<float>  minitree_muon_jet_dRmax;
   std::vector<float>  minitree_elemu_jet_dRmin;
   std::vector<float>  minitree_elemu_jet_dRmax;
   std::vector<float>  minitree_ele_jet_dRmin;
   std::vector<float>  minitree_ele_jet_dRmax;

    std::vector< int >   minitree_Hemi;
    std::vector< int >   minitree_Hemi_njet;
    std::vector< int >   minitree_Hemi_njet_nomu;
    std::vector< float > minitree_Hemi_pt;
    std::vector< float > minitree_Hemi_eta;
    std::vector< float > minitree_Hemi_phi;
    std::vector< int >   minitree_Hemi_nTrks;
    std::vector< int >   minitree_Hemi_nTrks_sig;
    std::vector< int >   minitree_Hemi_nTrks_bad;
    std::vector< float > minitree_Hemi_mass;
    std::vector< float > minitree_HemiMu_mass;
    std::vector< float > minitree_HemiMu_pt;
    std::vector< float > minitree_HemiMu_dR;
    std::vector< float > minitree_HemiMuOp_mass;
    std::vector< float > minitree_HemiMuOp_pt;
    std::vector< float > minitree_HemiMuOp_dR;
    std::vector< int >   minitree_Hemi_LooseBTag_axes;
    std::vector< int >   minitree_Hemi_MediumBTag_axes;
    std::vector< int >   minitree_Hemi_TightBTag_axes;
    std::vector< float > minitree_Hemi_dR12;

   std::vector<int>  minitree_TRACK_SIZE;
   std::vector<int>  minitree_nTracks;
   std::vector<int>  minitree_nLostTracks;

   std::vector<unsigned int> minitree_track_ipc;
   std::vector<bool>    minitree_track_lost;
   std::vector<float>   minitree_track_px;
   std::vector<float>   minitree_track_py;
   std::vector<float>   minitree_track_pz;
   std::vector<float>   minitree_track_pt;
   std::vector<float>   minitree_track_eta;
   std::vector<float>   minitree_track_phi;
   std::vector<int>     minitree_track_charge;
   std::vector<float>   minitree_track_NChi2;
   std::vector<bool>    minitree_track_isHighPurity;
   std::vector<float>   minitree_track_dxy;
   std::vector<float>   minitree_track_dxyError;
   std::vector<float>   minitree_track_drSig;
   std::vector<float>   minitree_track_dz;
   std::vector<float>   minitree_track_dzError;
   std::vector<float>   minitree_track_dzSig;
   std::vector<int>     minitree_track_nHit;
   std::vector<int>     minitree_track_nHitPixel;
   std::vector<int>     minitree_track_nHitTIB;
   std::vector<int>     minitree_track_nHitTID;
   std::vector<int>     minitree_track_nHitTOB;
   std::vector<int>     minitree_track_nHitTEC;
   std::vector<int>     minitree_track_nHitPXB;
   std::vector<int>     minitree_track_nHitPXF;
   std::vector<int>     minitree_track_isHitPixel;
   std::vector<int>     minitree_track_nLayers;
   std::vector<int>     minitree_track_nLayersPixel;
   std::vector<float>   minitree_track_x;
   std::vector<float>   minitree_track_y;
   std::vector<float>   minitree_track_z;
   std::vector<int>     minitree_track_firstHit;
   std::vector<float>   minitree_track_region;
   std::vector<float>   minitree_track_firstHit_x;
   std::vector<float>   minitree_track_firstHit_y;
   std::vector<float>   minitree_track_firstHit_z;
   std::vector<int>     minitree_track_iJet;
   std::vector<float>   minitree_track_ntrk10;
   std::vector<float>   minitree_track_ntrk20;
   std::vector<float>   minitree_track_ntrk30;
   std::vector<float>   minitree_track_ntrk40;
   std::vector<double>  minitree_track_MVAval;
   std::vector<double>  minitree_track_Hemi_dR;
   std::vector<double>  minitree_track_Hemi_dRmax;
     
   std::vector<float>   minitree_K0_mass;
   std::vector<float>   minitree_K0_pt;
   std::vector<float>   minitree_L0_mass;
   std::vector<float>   minitree_L0_pt;
   std::vector<float>   minitree_V0_reco_mass;
   std::vector<float>   minitree_V0_reco_pt;
   std::vector<int>     minitree_V0_reco_source;

   std::vector<float>   minitree_SecInt_drSig;
   std::vector<float>   minitree_SecInt_dzSig;
   std::vector<float>   minitree_SecInt_mass;
   std::vector<float>   minitree_SecInt_pt;
   std::vector<bool>    minitree_SecInt_selec;
   std::vector<int>     minitree_SecInt_layer;
   std::vector<float>   minitree_SecInt_r;
   std::vector<float>   minitree_SecInt_z;


   std::vector<int>     minitree_Hemi_Vtx_step;
   std::vector<float>   minitree_Hemi_Vtx_NChi2;
   std::vector<int>     minitree_Hemi_Vtx_nTrks;
   std::vector<float>   minitree_Hemi_Vtx_dR;
   std::vector<float>   minitree_Hemi_Vtx_xError;
   std::vector<float>   minitree_Hemi_Vtx_yError;
   std::vector<float>   minitree_Hemi_Vtx_zError;
   std::vector<float>   minitree_Hemi_Vtx_trackWeight;
   std::vector<float>   minitree_Hemi_Vtx_SumtrackWeight;
   std::vector<float>   minitree_Hemi_Vtx_track_MeanDCA_d;
   std::vector<float>   minitree_Hemi_Vtx_Mass;
   std::vector<float>   minitree_Hemi_Vtx_dist;
   std::vector<int>     minitree_Hemi_Vtx_ntrk10;
   std::vector<int>     minitree_Hemi_Vtx_ntrk20;

   std::vector<int>     minitree_Hemi_SecVtx_step;
   std::vector<float>   minitree_Hemi_SecVtx_dR;
   std::vector<float>   minitree_Hemi_SecVtx_nTrks;
   std::vector<float>   minitree_Hemi_SecVtx_NChi2;
   std::vector<float>   minitree_Hemi_SecVtx_dist;
   std::vector<float>   minitree_Hemi_SecVtx_track_MeanDCA_d;
   std::vector<float>   minitree_Hemi_SecVtx_SumtrackWeight;
   std::vector<float>   minitree_Hemi_SecVtx_trackWeight;
   std::vector<float>   minitree_Hemi_SecVtx_Mass;

   std::vector<float>   minitree_Hemi_Vtx_BDT_nTrks;
   std::vector<float>   minitree_Hemi_Vtx_BDT_NChi2;
   std::vector<float>   minitree_Hemi_Vtx_BDT_step;
   std::vector<float>   minitree_Hemi_Vtx_BDT_Mass;
   std::vector<float>   minitree_Hemi_Vtx_BDT_HMass;


   smalltree->Branch("minitree_only_gen_wt",&minitree_only_gen_wt);

   smalltree->Branch("minitree_genTop_Weight",&minitree_genTop_Weight);
   smalltree->Branch("minitree_gen_top_pt",&minitree_gen_top_pt);
   smalltree->Branch("minitree_gen_top_rw_pt",&minitree_gen_top_rw_pt);

   smalltree->Branch("miniPUweight",&miniPUweight);
   smalltree->Branch("miniPUweight_Up",&miniPUweight_Up);
   smalltree->Branch("miniPUweight_Down",&miniPUweight_Down);
   smalltree->Branch("miniPrefweight",&miniPrefweight);

   smalltree->Branch("minitree_Filter",        &minitree_Filter);
   smalltree->Branch("minitree_FilterSameSign",&minitree_FilterSameSign);

   smalltree->Branch("minitree_nPV",&minitree_nPV);
   smalltree->Branch("minitree_trigger_doublelepton",&minitree_trigger_doublelepton);
   smalltree->Branch("minitree_trigger_singlelepton",&minitree_trigger_singlelepton);

   smalltree->Branch("minitree_Good_PV",&minitree_Good_PV);

   smalltree->Branch("minitree_smu_mass",&minitree_smu_mass); // used for signal only
   smalltree->Branch("minitree_neu_mass",&minitree_neu_mass); // used for signal only
   smalltree->Branch("minitree_neu_ctau",&minitree_neu_ctau); // used for signal only


   smalltree->Branch("minitree_HT",&minitree_HT); // in BDT_EVT

   smalltree->Branch("minitree_TRACK_SIZE",&minitree_TRACK_SIZE); // in BDT_EVT
   smalltree->Branch("minitree_nTracks",&minitree_nTracks);
   smalltree->Branch("minitree_nLostTracks",&minitree_nLostTracks);

   // smalltree->Branch("minitree_muon_GenRecoTriggerMatched",&minitree_muon_GenRecoTriggerMatched);
   smalltree->Branch("minitree_all_nmu",&minitree_all_nmu); 
   smalltree->Branch("minitree_nmu",&minitree_nmu);    
   smalltree->Branch("minitree_LT",&minitree_LT);
   smalltree->Branch("minitree_Mmumu",&minitree_Mmumu);
   smalltree->Branch("minitree_MmumuSameSign",&minitree_MmumuSameSign);

   smalltree->Branch("minitree_muon_isPrompt",&minitree_muon_isPrompt);
   smalltree->Branch("minitree_muon_pt",&minitree_muon_pt);
   smalltree->Branch("minitree_muon_SF",&minitree_muon_SF);
   smalltree->Branch("minitree_muon_eta",&minitree_muon_eta);
   smalltree->Branch("minitree_muon_phi",&minitree_muon_phi);
   smalltree->Branch("minitree_muon_dxy",&minitree_muon_dxy);
   smalltree->Branch("minitree_muon_dz",&minitree_muon_dz);
   smalltree->Branch("minitree_muon_charge",&minitree_muon_charge);
   smalltree->Branch("minitree_muon_correction",&minitree_muon_correction);
   // smalltree->Branch("minitree_muon_gen",&minitree_muon_gen);
   smalltree->Branch("minitree_muon_dxyError",&minitree_muon_dxyError);
   smalltree->Branch("minitree_muon_dzError",&minitree_muon_dzError);
   smalltree->Branch("minitree_muon_isLoose",&minitree_muon_isLoose);
   smalltree->Branch("minitree_muon_isMedium",&minitree_muon_isMedium);
   smalltree->Branch("minitree_muon_isTight",&minitree_muon_isTight);
   smalltree->Branch("minitree_muon_isGlobal",&minitree_muon_isGlobal);
   smalltree->Branch("minitree_muon_PFIsoVeryLoose",&minitree_muon_PFIsoVeryLoose);
   smalltree->Branch("minitree_muon_PFIsoLoose",&minitree_muon_PFIsoLoose);
   smalltree->Branch("minitree_muon_PFIsoMedium",&minitree_muon_PFIsoMedium);
   smalltree->Branch("minitree_muon_PFIsoTight",&minitree_muon_PFIsoTight);
   smalltree->Branch("minitree_muon_TkIsoLoose",&minitree_muon_TkIsoLoose);
   smalltree->Branch("minitree_muon_TkIsoTight",&minitree_muon_TkIsoTight);
   smalltree->Branch("minitree_muon_MiniIsoLoose",&minitree_muon_MiniIsoLoose);
   smalltree->Branch("minitree_muon_MiniIsoMedium",&minitree_muon_MiniIsoMedium);
   smalltree->Branch("minitree_muon_MiniIsoTight",&minitree_muon_MiniIsoTight);

   smalltree->Branch("minitree_lepton_leadingpt",&minitree_lepton_leadingpt);
   smalltree->Branch("minitree_lepton_leadingpt2",&minitree_lepton_leadingpt2);
   smalltree->Branch("minitree_lepton_leadingeta",&minitree_lepton_leadingeta);
   smalltree->Branch("minitree_lepton_leadingeta2",&minitree_lepton_leadingeta2);
   smalltree->Branch("minitree_lepton_leadingphi",&minitree_lepton_leadingphi);
   smalltree->Branch("minitree_lepton_leadingphi2",&minitree_lepton_leadingphi2);

   smalltree->Branch("minitree_lepton_lepton_dR",&minitree_lepton_lepton_dR);
   smalltree->Branch("minitree_lepton_lepton_dPhi",&minitree_lepton_lepton_dPhi);
   smalltree->Branch("minitree_lepton_lepton_dEta",&minitree_lepton_lepton_dEta);
   smalltree->Branch("minitree_lepton_leadingdxy",&minitree_lepton_leadingdxy);
   smalltree->Branch("minitree_lepton_leadingdxy2",&minitree_lepton_leadingdxy2);
   smalltree->Branch("minitree_lepton_leadingdz",&minitree_lepton_leadingdz);
   smalltree->Branch("minitree_lepton_leadingdz2",&minitree_lepton_leadingdz2);


   // smalltree->Branch("minitree_reco_lepton_leadingpt",&minitree_reco_lepton_leadingpt);
   // smalltree->Branch("minitree_reco_lepton_leadingpt2",&minitree_reco_lepton_leadingpt2);
   // smalltree->Branch("minitree_reco_lepton_leadingeta",&minitree_reco_lepton_leadingeta);
   // smalltree->Branch("minitree_reco_lepton_leadingeta2",&minitree_reco_lepton_leadingeta2);
   // smalltree->Branch("minitree_reco_lepton_leadingphi",&minitree_reco_lepton_leadingphi);
   // smalltree->Branch("minitree_reco_lepton_leadingphi2",&minitree_reco_lepton_leadingphi2);

   // smalltree->Branch("minitree_lepton_b4trigger_leadingpt",&minitree_lepton_b4trigger_leadingpt);
   // smalltree->Branch("minitree_lepton_b4trigger_leadingpt2",&minitree_lepton_b4trigger_leadingpt2);


   smalltree->Branch("minitree_all_nel",&minitree_all_nel); 
   smalltree->Branch("minitree_electron_nEle",&minitree_electron_nEle); 
   smalltree->Branch("minitree_electron_isPrompt",&minitree_electron_isPrompt);
   smalltree->Branch("minitree_electron_pt",&minitree_electron_pt);
   smalltree->Branch("minitree_electron_eta",&minitree_electron_eta);
   smalltree->Branch("minitree_electron_phi",&minitree_electron_phi);
   smalltree->Branch("minitree_electron_charge",&minitree_electron_charge);
   smalltree->Branch("minitree_electron_dxy",&minitree_electron_dxy);
   smalltree->Branch("minitree_electron_dz",&minitree_electron_dz);
   smalltree->Branch("minitree_electron_gen",&minitree_electron_gen);
   smalltree->Branch("minitree_electron_energy",&minitree_electron_energy);
   smalltree->Branch("minitree_electron_et",&minitree_electron_et);
   smalltree->Branch("minitree_electron_ecal_trk_postcorr",&minitree_electron_ecal_trk_postcorr);
   smalltree->Branch("minitree_electron_isoR4",&minitree_electron_isoR4);
   smalltree->Branch("minitree_electron_IsLoose",&minitree_electron_IsLoose);
   smalltree->Branch("minitree_electron_IsMedium",&minitree_electron_IsMedium);
   smalltree->Branch("minitree_electron_IsTight",&minitree_electron_IsTight);

   smalltree->Branch("minitree_njet",&minitree_njet);
   smalltree->Branch("minitree_njetNOmu",&minitree_njetNOmu); // in BDT_EVT
   smalltree->Branch("minitree_jet_pt",&minitree_jet_pt);
   smalltree->Branch("minitree_jet_eta",&minitree_jet_eta);
   smalltree->Branch("minitree_jet_phi",&minitree_jet_phi);
   smalltree->Branch("minitree_jet_HadronFlavour",&minitree_jet_HadronFlavour);
   smalltree->Branch("minitree_jet_btag_DeepJet",&minitree_jet_btag_DeepJet);
   smalltree->Branch("minitree_jet_E",&minitree_jet_E);

   smalltree->Branch("minitree_jet_leadingpt",&minitree_jet_leadingpt); // in BDT_EVT
   smalltree->Branch("minitree_jet_leadingpt2",&minitree_jet_leadingpt2); // in BDT_EVT
   smalltree->Branch("minitree_jet_leadingeta",&minitree_jet_leadingeta); // in BDT_EVT
   smalltree->Branch("minitree_jet_leadingeta2",&minitree_jet_leadingeta2); // in BDT_EVT
   smalltree->Branch("minitree_jet_jet_dR",&minitree_jet_jet_dR); // in BDT_EVT
   smalltree->Branch("minitree_jet_jet_dPhi",&minitree_jet_jet_dPhi); // in BDT_EVT
   smalltree->Branch("minitree_jet_jet_dEta",&minitree_jet_jet_dEta); // in BDT_EVT
   smalltree->Branch("minitree_muon_jet_dRmin",&minitree_muon_jet_dRmin);
   smalltree->Branch("minitree_muon_jet_dRmax",&minitree_muon_jet_dRmax);
   smalltree->Branch("minitree_elemu_jet_dRmin",&minitree_elemu_jet_dRmin);
   smalltree->Branch("minitree_elemu_jet_dRmax",&minitree_elemu_jet_dRmax);
   smalltree->Branch("minitree_ele_jet_dRmin",&minitree_ele_jet_dRmin); // empty, usefull ???
   smalltree->Branch("minitree_ele_jet_dRmax",&minitree_ele_jet_dRmax); // empty, usefull ???

   // smalltree->Branch("minitree_Evts_MVAval",   &minitree_Evts_MVAval);
   // smalltree->Branch("minitree_Evts_MVAvalDY",   &minitree_Evts_MVAvalDY);
   // smalltree->Branch("minitree_Evts_MVAvalTT",   &minitree_Evts_MVAvalTT);
   smalltree->Branch("minitree_Hemi",&minitree_Hemi);
   smalltree->Branch("minitree_Hemi_njet",&minitree_Hemi_njet);
   smalltree->Branch("minitree_Hemi_njet_nomu",&minitree_Hemi_njet_nomu);
   smalltree->Branch("minitree_Hemi_pt",&minitree_Hemi_pt);
   smalltree->Branch("minitree_Hemi_eta",&minitree_Hemi_eta);
   smalltree->Branch("minitree_Hemi_phi",&minitree_Hemi_phi);
   smalltree->Branch("minitree_Hemi_nTrks",&minitree_Hemi_nTrks);
   smalltree->Branch("minitree_Hemi_nTrks_sig",&minitree_Hemi_nTrks_sig);
   smalltree->Branch("minitree_Hemi_nTrks_bad",&minitree_Hemi_nTrks_bad);
   smalltree->Branch("minitree_Hemi_mass",&minitree_Hemi_mass);
   smalltree->Branch("minitree_HemiMu_mass",&minitree_HemiMu_mass);
   smalltree->Branch("minitree_HemiMu_pt",&minitree_HemiMu_pt);
   smalltree->Branch("minitree_HemiMu_dR",&minitree_HemiMu_dR);
   smalltree->Branch("minitree_HemiMuOp_mass",&minitree_HemiMuOp_mass);
   smalltree->Branch("minitree_HemiMuOp_pt",&minitree_HemiMuOp_pt);
   smalltree->Branch("minitree_HemiMuOp_dR",&minitree_HemiMuOp_dR);
   smalltree->Branch("minitree_Hemi_dR12",&minitree_Hemi_dR12);

   smalltree->Branch("minitree_ll_pt",&minitree_ll_pt);
   smalltree->Branch("minitree_ll_eta",&minitree_ll_eta);
   smalltree->Branch("minitree_ll_phi",&minitree_ll_phi);
   smalltree->Branch("minitree_ll_px",&minitree_ll_px);
   smalltree->Branch("minitree_ll_py",&minitree_ll_py);
   smalltree->Branch("minitree_ll_pz",&minitree_ll_pz);
   smalltree->Branch("minitree_ll_energy",&minitree_ll_energy);
   smalltree->Branch("minitree_ll_mass",&minitree_ll_mass);


   smalltree->Branch("minitree_track_ipc",&minitree_track_ipc);
   smalltree->Branch("minitree_track_lost",&minitree_track_lost);
   smalltree->Branch("minitree_track_px",&minitree_track_px);
   smalltree->Branch("minitree_track_py",&minitree_track_py);
   smalltree->Branch("minitree_track_pz",&minitree_track_pz);
   smalltree->Branch("minitree_track_pt",&minitree_track_pt);
   smalltree->Branch("minitree_track_eta",&minitree_track_eta);
   smalltree->Branch("minitree_track_phi",&minitree_track_phi);
   smalltree->Branch("minitree_track_charge",&minitree_track_charge);
   smalltree->Branch("minitree_track_NChi2",&minitree_track_NChi2);
   smalltree->Branch("minitree_track_isHighPurity",&minitree_track_isHighPurity);
   smalltree->Branch("minitree_track_dxy",&minitree_track_dxy);
   smalltree->Branch("minitree_track_dxyError",&minitree_track_dxyError);
   smalltree->Branch("minitree_track_drSig",&minitree_track_drSig);
   smalltree->Branch("minitree_track_dz",&minitree_track_dz);
   smalltree->Branch("minitree_track_dzError",&minitree_track_dzError);
   smalltree->Branch("minitree_track_dzSig",&minitree_track_dzSig);
   smalltree->Branch("minitree_track_nHit",&minitree_track_nHit);
   smalltree->Branch("minitree_track_nHitPixel",&minitree_track_nHitPixel);
   smalltree->Branch("minitree_track_nHitTIB",&minitree_track_nHitTIB);
   smalltree->Branch("minitree_track_nHitTID",&minitree_track_nHitTID);
   smalltree->Branch("minitree_track_nHitTOB",&minitree_track_nHitTOB);
   smalltree->Branch("minitree_track_nHitTEC",&minitree_track_nHitTEC);
   smalltree->Branch("minitree_track_nHitPXB",&minitree_track_nHitPXB);
   smalltree->Branch("minitree_track_nHitPXF",&minitree_track_nHitPXF);
   smalltree->Branch("minitree_track_isHitPixel",&minitree_track_isHitPixel);
   smalltree->Branch("minitree_track_nLayers",&minitree_track_nLayers);
   smalltree->Branch("minitree_track_nLayersPixel",&minitree_track_nLayersPixel);
   smalltree->Branch("minitree_track_x",&minitree_track_x);
   smalltree->Branch("minitree_track_y",&minitree_track_y);
   smalltree->Branch("minitree_track_z",&minitree_track_z);
   smalltree->Branch("minitree_track_firstHit",&minitree_track_firstHit);
   smalltree->Branch("minitree_track_region",&minitree_track_region);
   smalltree->Branch("minitree_track_firstHit_x",&minitree_track_firstHit_x);
   smalltree->Branch("minitree_track_firstHit_y",&minitree_track_firstHit_y);
   smalltree->Branch("minitree_track_firstHit_z",&minitree_track_firstHit_z);
   smalltree->Branch("minitree_track_iJet",&minitree_track_iJet);
   smalltree->Branch("minitree_track_ntrk10",&minitree_track_ntrk10);
   smalltree->Branch("minitree_track_ntrk20",&minitree_track_ntrk20);
   smalltree->Branch("minitree_track_ntrk30",&minitree_track_ntrk30);
   smalltree->Branch("minitree_track_ntrk40",&minitree_track_ntrk40);
   smalltree->Branch("minitree_track_MVAval",&minitree_track_MVAval);

// !! 
   smalltree->Branch("minitree_track_Hemi_dR",&minitree_track_Hemi_dR);
   smalltree->Branch("minitree_track_Hemi_dRmax",&minitree_track_Hemi_dRmax);


   smalltree->Branch("minitree_K0_mass",&minitree_K0_mass);
   smalltree->Branch("minitree_K0_pt",&minitree_K0_mass);

   smalltree->Branch("minitree_L0_mass",&minitree_L0_mass);
   smalltree->Branch("minitree_L0_pt",&minitree_L0_pt);

   smalltree->Branch("minitree_V0_reco_mass",&minitree_V0_reco_mass);
   smalltree->Branch("minitree_V0_reco_pt",&minitree_V0_reco_pt);
   smalltree->Branch("minitree_V0_reco_source",&minitree_V0_reco_source);

   smalltree->Branch("minitree_SecInt_mass",&minitree_SecInt_mass);
   smalltree->Branch("minitree_SecInt_pt",&minitree_SecInt_pt);
   smalltree->Branch("minitree_SecInt_drSig",&minitree_SecInt_drSig);
   smalltree->Branch("minitree_SecInt_dzSig",&minitree_SecInt_dzSig);
   smalltree->Branch("minitree_SecInt_layer",&minitree_SecInt_layer);
   smalltree->Branch("minitree_SecInt_selec",&minitree_SecInt_selec);
   smalltree->Branch("minitree_SecInt_r",&minitree_SecInt_r);
   smalltree->Branch("minitree_SecInt_z",&minitree_SecInt_z);

   smalltree->Branch("minitree_Hemi_Vtx_step",&minitree_Hemi_Vtx_step);
   smalltree->Branch("minitree_Hemi_Vtx_NChi2",&minitree_Hemi_Vtx_NChi2);
   smalltree->Branch("minitree_Hemi_Vtx_nTrks",&minitree_Hemi_Vtx_nTrks);
   smalltree->Branch("minitree_Hemi_Vtx_dR",&minitree_Hemi_Vtx_dR);
   smalltree->Branch("minitree_Hemi_Vtx_xError",&minitree_Hemi_Vtx_xError);
   smalltree->Branch("minitree_Hemi_Vtx_yError",&minitree_Hemi_Vtx_yError);
   smalltree->Branch("minitree_Hemi_Vtx_zError",&minitree_Hemi_Vtx_zError);
   smalltree->Branch("minitree_Hemi_Vtx_trackWeight",&minitree_Hemi_Vtx_trackWeight);
   smalltree->Branch("minitree_Hemi_Vtx_SumtrackWeight",&minitree_Hemi_Vtx_SumtrackWeight);
   smalltree->Branch("minitree_Hemi_Vtx_track_MeanDCA_d",&minitree_Hemi_Vtx_track_MeanDCA_d);
   smalltree->Branch("minitree_Hemi_Vtx_Mass",&minitree_Hemi_Vtx_Mass);
   smalltree->Branch("minitree_Hemi_Vtx_dist",&minitree_Hemi_Vtx_dist);
   smalltree->Branch("minitree_Hemi_Vtx_ntrk10",&minitree_Hemi_Vtx_ntrk10);
   smalltree->Branch("minitree_Hemi_Vtx_ntrk20",&minitree_Hemi_Vtx_ntrk20);

   smalltree->Branch("minitree_Hemi_SecVtx_step",&minitree_Hemi_SecVtx_step);
   smalltree->Branch("minitree_Hemi_SecVtx_NChi2",&minitree_Hemi_SecVtx_NChi2);
   smalltree->Branch("minitree_Hemi_SecVtx_dR",&minitree_Hemi_SecVtx_dR);
   smalltree->Branch("minitree_Hemi_SecVtx_nTrks",&minitree_Hemi_SecVtx_nTrks);
   smalltree->Branch("minitree_Hemi_SecVtx_dist",&minitree_Hemi_SecVtx_dist);
   smalltree->Branch("minitree_Hemi_SecVtx_track_MeanDCA_d",&minitree_Hemi_SecVtx_track_MeanDCA_d);
   smalltree->Branch("minitree_Hemi_SecVtx_SumtrackWeight",&minitree_Hemi_SecVtx_SumtrackWeight);
   smalltree->Branch("minitree_Hemi_SecVtx_trackWeight",&minitree_Hemi_SecVtx_trackWeight);
   smalltree->Branch("minitree_Hemi_SecVtx_Mass",&minitree_Hemi_SecVtx_Mass);

   smalltree->Branch("minitree_Hemi_Vtx_BDT_nTrks",&minitree_Hemi_Vtx_BDT_nTrks);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_NChi2",&minitree_Hemi_Vtx_BDT_NChi2);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_step",&minitree_Hemi_Vtx_BDT_step);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_Mass",&minitree_Hemi_Vtx_BDT_Mass);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_HMass",&minitree_Hemi_Vtx_BDT_HMass);

   std::cout<<"//-------------------------//"<<std::endl;
   std::cout<<"// "<<Production<<" : "<<sample<<"  //"<<std::endl;
   std::cout<<"//-------------------------//"<<std::endl;

   bool debug = false;
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   Long64_t nentries2 = fChain->GetEntriesFast();
   std::cout << "Total Entries : " << nentries << std::endl;
   std::cout << "Total Entries2 : " << nentries2 << std::endl;
   Long64_t nbytes = 0, nb = 0;

   int allevents = 0;
   int itest = 0;

   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      allevents++;
      // if ( allevents%10000 == 0 ) std::cout << "events : " << allevents << std::endl;

      // if (Cut(ientry) < 0) continue;
      if ( jentry%10000 == 0 ) std::cout << "events : " << jentry << std::endl;
      // if (jentry < 3300000) continue;
      // std::cout<<"jentry : "<<jentry<<std::endl;
      // if (jentry==0.1*nentries) {std::cout<<"10/100 :"<<std::endl;}
      // if (jentry==0.2*nentries) {std::cout<<"20/100 :"<<std::endl;}
      // if (jentry==0.3*nentries) {std::cout<<"30/100 :"<<std::endl;}
      // if (jentry==0.4*nentries) {std::cout<<"40/100 :"<<std::endl;}
      // if (jentry==0.5*nentries) {std::cout<<"50/100 :"<<std::endl;}
      // if (jentry==0.6*nentries) {std::cout<<"60/100 :"<<std::endl;}
      // if (jentry==0.7*nentries) {std::cout<<"70/100 :"<<std::endl;}
      // if (jentry==0.8*nentries) {std::cout<<"80/100 :"<<std::endl;}
      // if (jentry==0.9*nentries) {std::cout<<"90/100 :"<<std::endl;}

//$$ 
      itest++;
//    if ( itest > 1000000 ) break;
//$$

      //----------------------//
      //        Event         //
      //----------------------//



   // if (tree_Mmumu >10)
   //    {
   //       minitree_lepton_b4trigger_leadingpt.push_back(tree_lepton_b4trigger_leadingpt->at(0));
   //       minitree_lepton_b4trigger_leadingpt2.push_back(tree_lepton_b4trigger_leadingpt2->at(0));
   //    }

//$$
   
    // -------------- --------------------------------------------------------------------//
    // -----------------------------------------------------------------------------------//
    if ( !((tree_Filter || tree_FilterSameSign) && tree_njetNOmu>0) && !Signal ) continue;
    //------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------//
//$$
      minitree_nPV.push_back(tree_nPV);
      minitree_trigger_doublelepton.push_back(tree_trigger_doublelepton);
      minitree_trigger_singlelepton.push_back(tree_trigger_singlelepton);

      minitree_Good_PV.push_back(tree_Good_PV);

      for(unsigned int ill = 0; ill <tree_ll_pt->size(); ill ++)//test                                                                                                        
          {
            minitree_ll_pt.push_back(tree_ll_pt->at(ill));
            minitree_ll_eta.push_back(tree_ll_eta->at(ill));
            minitree_ll_phi.push_back(tree_ll_phi->at(ill));
            minitree_ll_mass.push_back(tree_ll_mass->at(ill));
         }
         // std::cout<< "tracks   "<<std::endl;
      for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++)
         {
            //--- Track pre-selected ---------//  
            

            minitree_track_ipc.push_back(tree_track_ipc->at(iTrk));
            minitree_track_lost.push_back(tree_track_lost->at(iTrk));
            minitree_track_px.push_back(tree_track_px->at(iTrk));
            minitree_track_py.push_back(tree_track_py->at(iTrk));
            minitree_track_pz.push_back(tree_track_pz->at(iTrk));
            minitree_track_pt.push_back(tree_track_pt->at(iTrk));
            minitree_track_eta.push_back(tree_track_eta->at(iTrk));
            minitree_track_phi.push_back(tree_track_phi->at(iTrk));
            minitree_track_charge.push_back(tree_track_charge->at(iTrk));
            minitree_track_NChi2.push_back(tree_track_NChi2->at(iTrk));
            minitree_track_isHighPurity.push_back(tree_track_isHighPurity->at(iTrk));
            minitree_track_dxy.push_back(tree_track_dxy->at(iTrk));
            minitree_track_dxyError.push_back(tree_track_dxyError->at(iTrk));
            minitree_track_drSig.push_back(tree_track_drSig->at(iTrk));
            minitree_track_dz.push_back(tree_track_dz->at(iTrk));
            minitree_track_dzError.push_back(tree_track_dzError->at(iTrk));
            minitree_track_dzSig.push_back(tree_track_dzSig->at(iTrk));
            minitree_track_nHit.push_back(tree_track_nHit->at(iTrk));
            minitree_track_nHitPixel.push_back(tree_track_nHitPixel->at(iTrk));
            minitree_track_nHitTIB.push_back(tree_track_nHitTIB->at(iTrk));
            minitree_track_nHitTID.push_back(tree_track_nHitTID->at(iTrk));
            minitree_track_nHitTOB.push_back(tree_track_nHitTOB->at(iTrk));
            minitree_track_nHitTEC.push_back(tree_track_nHitTEC->at(iTrk));
            minitree_track_nHitPXB.push_back(tree_track_nHitPXB->at(iTrk));
            minitree_track_nHitPXF.push_back(tree_track_nHitPXF->at(iTrk));
            minitree_track_isHitPixel.push_back(tree_track_isHitPixel->at(iTrk));
            minitree_track_nLayers.push_back(tree_track_nLayers->at(iTrk));
            minitree_track_nLayersPixel.push_back(tree_track_nLayersPixel->at(iTrk));
            minitree_track_x.push_back(tree_track_x->at(iTrk));
            minitree_track_y.push_back(tree_track_y->at(iTrk));
            minitree_track_z.push_back(tree_track_z->at(iTrk));
            minitree_track_firstHit.push_back(tree_track_firstHit->at(iTrk));
            minitree_track_region.push_back(tree_track_region->at(iTrk));
            minitree_track_firstHit_x.push_back(tree_track_firstHit_x->at(iTrk));
            minitree_track_firstHit_y.push_back(tree_track_firstHit_y->at(iTrk));
            minitree_track_firstHit_z.push_back(tree_track_firstHit_z->at(iTrk));
            minitree_track_iJet.push_back(tree_track_iJet->at(iTrk));
            minitree_track_ntrk10.push_back(tree_track_ntrk10->at(iTrk));
            minitree_track_ntrk20.push_back(tree_track_ntrk20->at(iTrk));
            minitree_track_ntrk30.push_back(tree_track_ntrk30->at(iTrk));
            minitree_track_ntrk40.push_back(tree_track_ntrk40->at(iTrk));
            minitree_track_MVAval.push_back(tree_track_MVAval->at(iTrk));


            minitree_track_Hemi_dR.push_back(tree_track_Hemi_dR->at(iTrk));
            minitree_track_Hemi_dRmax.push_back(tree_track_Hemi_dRmax->at(iTrk));

         }//end loop on tracks


      minitree_Filter.push_back(tree_Filter);
      minitree_FilterSameSign.push_back(tree_FilterSameSign);
      // if ( tree_Filter && tree_njetNOmu > 0 ) hData_nPV->Fill( tree_nPV );
         // std::cout<< "muon  "<<std::endl;
      if (debug)  {std::cout<<"muon  "<<std::endl;}
      for (unsigned int i = 0; i < tree_muon_pt->size(); i++)
         {
            minitree_muon_isPrompt.push_back(tree_muon_isPrompt->at(i));
            minitree_muon_pt.push_back(tree_muon_pt->at(i));
            minitree_muon_SF.push_back(tree_muon_SF->at(i));
            minitree_muon_eta.push_back(tree_muon_eta->at(i));
            minitree_muon_phi.push_back(tree_muon_phi->at(i));
            minitree_muon_dxy.push_back(tree_muon_dxy->at(i));
            minitree_muon_dz.push_back(tree_muon_dz->at(i));
            minitree_muon_charge.push_back(tree_muon_charge->at(i));
            minitree_muon_correction.push_back(tree_muon_correction->at(i));
            // minitree_muon_gen.push_back(tree_muon_gen->at(i));
            minitree_muon_dxyError.push_back( tree_muon_dxyError->at(i));
            minitree_muon_dzError.push_back( tree_muon_dzError->at(i));
            minitree_muon_isLoose.push_back( tree_muon_isLoose->at(i));
            minitree_muon_isMedium.push_back( tree_muon_isMedium->at(i));
            minitree_muon_isTight.push_back(  tree_muon_isTight->at(i));
            minitree_muon_isGlobal.push_back( tree_muon_isGlobal->at(i));
            minitree_muon_PFIsoVeryLoose.push_back(tree_muon_PFIsoVeryLoose->at(i));
            minitree_muon_PFIsoLoose.push_back(tree_muon_PFIsoLoose->at(i));
            minitree_muon_PFIsoMedium.push_back(tree_muon_PFIsoMedium->at(i));
            minitree_muon_PFIsoTight.push_back( tree_muon_PFIsoTight->at(i));
            minitree_muon_TkIsoLoose.push_back( tree_muon_TkIsoLoose->at(i));
            minitree_muon_TkIsoTight.push_back( tree_muon_TkIsoTight->at(i));
            minitree_muon_MiniIsoLoose.push_back( tree_muon_MiniIsoLoose->at(i));  //Wrongly stored
            minitree_muon_MiniIsoMedium.push_back(tree_muon_MiniIsoMedium->at(i)); //Wrongly stored
            minitree_muon_MiniIsoTight.push_back( tree_muon_MiniIsoTight->at(i));  //Wrongly stored
         }
      // std::cout<< "electron   "<<std::endl;
      if (debug)  {std::cout<<"electron "<<std::endl;}
      for (unsigned int i = 0; i < tree_electron_pt->size(); i++)
         {
            minitree_electron_isPrompt.push_back(tree_electron_isPrompt->at(i));
            minitree_electron_pt.push_back(tree_electron_pt->at(i));
            minitree_electron_eta.push_back(tree_electron_eta->at(i));
            minitree_electron_phi.push_back(tree_electron_phi->at(i));
            minitree_electron_charge.push_back(tree_electron_charge->at(i));
            minitree_electron_dxy.push_back(tree_electron_dxy->at(i));
            minitree_electron_dz.push_back(tree_electron_dz->at(i));
            minitree_electron_gen.push_back(tree_electron_gen->at(i));
            minitree_electron_energy.push_back(tree_electron_energy->at(i));
            minitree_electron_et.push_back(tree_electron_et->at(i));
            minitree_electron_ecal_trk_postcorr.push_back(tree_electron_ecal_trk_postcorr->at(i));
            minitree_electron_isoR4.push_back(tree_electron_isoR4->at(i));
            minitree_electron_IsLoose.push_back(tree_electron_IsLoose->at(i));
            minitree_electron_IsMedium.push_back(tree_electron_IsMedium->at(i));
            minitree_electron_IsTight.push_back(tree_electron_IsTight->at(i));
         }

      // hData_njetNOmu->Fill( tree_njetNOmu );

      if (debug)  {std::cout<<"Weights  "<<std::endl;}
      // minirunNumber.push_back(runNumber);
      // minieventNumber.push_back(eventNumber);
      // minilumiBlock.push_back(lumiBlock);
      // std::cout<< "weight   "<<std::endl;
      minitree_only_gen_wt.push_back(tree_only_gen_wt);
      minitree_genTop_Weight.push_back(tree_genTop_Weight);
      if (tree_gen_top_pt->size() > 0)
         {
            minitree_gen_top_pt.push_back(tree_gen_top_pt->at(0));
            minitree_gen_top_rw_pt.push_back(tree_gen_top_rw_pt->at(0));
         }

      miniPUweight.push_back(PUweight);
      double PileUp_Up = 1;
      double PileUp_Down = 1;

      PileUp_Up = PUweight_Up;
      PileUp_Down = PUweight_Down;

      miniPUweight_Up.push_back(PileUp_Up);
      miniPUweight_Down.push_back(PileUp_Down);
      miniPrefweight.push_back(Prefweight);
      // miniPU_events.push_back(PU_events);
      // miniAllPU_events_weight.push_back(AllPU_events_weight);

      if (Signal)
         {
            if (debug)  {std::cout<<"Signal parameters "<<std::endl;}
            minitree_smu_mass.push_back(tree_smu_mass); // used for signal only
            minitree_neu_mass.push_back(tree_neu_mass); // used for signal only
            minitree_neu_ctau.push_back(tree_neu_ctau); // used for signal only
         }

      // minitree_muon_GenRecoTriggerMatched.push_back(tree_muon_GenRecoTriggerMatched);

      if (debug)  {std::cout<<"PV (reco and Gen)  "<<std::endl;}


      minitree_njet.push_back(tree_njet);
      minitree_njetNOmu.push_back(tree_njetNOmu); // in BDT_EVT

      minitree_HT.push_back(tree_HT); // in BDT_EVT

      minitree_TRACK_SIZE.push_back(tree_TRACK_SIZE); // in BDT_EVT
      minitree_nTracks.push_back(tree_nTracks);
      minitree_nLostTracks.push_back(tree_nLostTracks);

      minitree_all_nmu.push_back(tree_all_nmu); 
      minitree_nmu.push_back(tree_nmu);    
      minitree_all_nel.push_back(tree_all_nel); 
      minitree_electron_nEle.push_back(tree_electron_nEle); 
      minitree_LT.push_back(tree_LT);
      minitree_Mmumu.push_back(tree_Mmumu);
      minitree_MmumuSameSign.push_back(tree_MmumuSameSign);


      // minitree_Evts_MVAval.push_back(tree_Evts_MVAval);
      // minitree_Evts_MVAvalDY.push_back(tree_Evts_MVAvalDY);
      // minitree_Evts_MVAvalTT.push_back(tree_Evts_MVAvalTT);


      // std::cout<< "jet   "<<std::endl;
      if (debug)  {std::cout<<"jet "<<std::endl;}
      for (unsigned int i = 0; i < tree_jet_pt->size(); i++)
         {
            minitree_jet_pt.push_back(tree_jet_pt->at(i));
            minitree_jet_eta.push_back(tree_jet_eta->at(i));
            minitree_jet_phi.push_back(tree_jet_phi->at(i));
            minitree_jet_HadronFlavour.push_back(tree_jet_HadronFlavour->at(i));
            minitree_jet_btag_DeepJet.push_back(tree_jet_btag_DeepJet->at(i));
            minitree_jet_E.push_back(tree_jet_E->at(i));
         }
      if ((tree_Filter || tree_FilterSameSign) && tree_njetNOmu>0)
         {
            minitree_lepton_leadingpt.push_back(tree_lepton_leadingpt->at(0));
            minitree_lepton_leadingpt2.push_back(tree_lepton_leadingpt2->at(0));
            minitree_lepton_leadingeta.push_back(tree_lepton_leadingeta->at(0));
            minitree_lepton_leadingeta2.push_back(tree_lepton_leadingeta2->at(0));
            minitree_lepton_leadingphi.push_back(tree_lepton_leadingphi->at(0));
            minitree_lepton_leadingphi2.push_back(tree_lepton_leadingphi2->at(0));

            minitree_lepton_lepton_dR.push_back(tree_lepton_lepton_dR->at(0));
            minitree_lepton_lepton_dPhi.push_back(tree_lepton_lepton_dPhi->at(0));
            minitree_lepton_lepton_dEta.push_back(tree_lepton_lepton_dEta->at(0));
            
            minitree_lepton_leadingdxy.push_back(tree_lepton_leadingdxy->at(0));
            minitree_lepton_leadingdxy2.push_back(tree_lepton_leadingdxy2->at(0));
            minitree_lepton_leadingdz.push_back(tree_lepton_leadingdz->at(0));
            minitree_lepton_leadingdz2.push_back(tree_lepton_leadingdz2->at(0));
            // std::cout<< "good pv  "<<std::endl;
            // if (tree_Good_PV)
            //    {
            //       minitree_reco_lepton_leadingpt.push_back(tree_reco_lepton_leadingpt->at(0) );
            //       minitree_reco_lepton_leadingpt2.push_back(tree_reco_lepton_leadingpt2->at(0) );
            //       minitree_reco_lepton_leadingeta.push_back(tree_reco_lepton_leadingeta->at(0) );
            //       minitree_reco_lepton_leadingeta2.push_back(tree_reco_lepton_leadingeta2->at(0) );
            //       minitree_reco_lepton_leadingphi.push_back(tree_reco_lepton_leadingphi->at(0) );
            //       minitree_reco_lepton_leadingphi2.push_back(tree_reco_lepton_leadingphi2->at(0) );

            //    }
         }
      if (debug)  {std::cout<<"Lepton-jet "<<std::endl;}
       if ((tree_Filter || tree_FilterSameSign) && tree_njetNOmu>0)
         {
            minitree_jet_leadingpt.push_back(tree_jet_leadingpt->at(0) ); // in BDT_EVT
            minitree_jet_leadingpt2.push_back(tree_jet_leadingpt2->at(0) ); // in BDT_EVT
            minitree_jet_leadingeta.push_back(tree_jet_leadingeta->at(0) ); // in BDT_EVT
            minitree_jet_leadingeta2.push_back(tree_jet_leadingeta2->at(0) ); // in BDT_EVT
            minitree_jet_jet_dR.push_back(tree_jet_jet_dR->at(0) ); // in BDT_EVT
            minitree_jet_jet_dPhi.push_back( tree_jet_jet_dPhi->at(0) ); // in BDT_EVT
            minitree_jet_jet_dEta.push_back( tree_jet_jet_dEta->at(0) ); // in BDT_EVT
            minitree_muon_jet_dRmin.push_back(tree_muon_jet_dRmin->at(0) );
            minitree_muon_jet_dRmax.push_back(tree_muon_jet_dRmax->at(0) );
         }

      // std::cout<< "Hemi   "<<std::endl;
      for (unsigned int i = 0; i < tree_Hemi->size(); i++)
         {
            minitree_Hemi.push_back(tree_Hemi->at(i));
            minitree_Hemi_njet.push_back(tree_Hemi_njet->at(i));
            minitree_Hemi_njet_nomu.push_back(tree_Hemi_njet_nomu->at(i));
            minitree_Hemi_pt.push_back(tree_Hemi_pt->at(i));
            minitree_Hemi_eta.push_back(tree_Hemi_eta->at(i));
            minitree_Hemi_phi.push_back(tree_Hemi_phi->at(i));
            minitree_Hemi_nTrks.push_back(tree_Hemi_nTrks->at(i));
            minitree_Hemi_nTrks_sig.push_back(tree_Hemi_nTrks_sig->at(i));
            minitree_Hemi_nTrks_bad.push_back(tree_Hemi_nTrks_bad->at(i));
            minitree_Hemi_mass.push_back(tree_Hemi_mass->at(i));
            minitree_HemiMu_mass.push_back(tree_HemiMu_mass->at(i));
            minitree_HemiMu_pt.push_back(tree_HemiMu_pt->at(i));
            minitree_HemiMu_dR.push_back(tree_HemiMu_dR->at(i));
            minitree_HemiMuOp_mass.push_back(tree_HemiMuOp_mass->at(i));
            minitree_HemiMuOp_pt.push_back(tree_HemiMuOp_pt->at(i));
            minitree_HemiMuOp_dR.push_back(tree_HemiMuOp_dR->at(i));
            minitree_Hemi_dR12.push_back(tree_Hemi_dR12->at(i));

         }
      // std::cout<< "V0 "<<std::endl;
      for (unsigned int i = 0 ; i < tree_K0_mass->size(); i++)
         {
               minitree_K0_mass.push_back(tree_K0_mass->at(i));
               minitree_K0_pt.push_back(tree_K0_pt->at(i));
         }
      for (unsigned int i = 0 ; i < tree_L0_mass->size(); i++)
         {
               minitree_L0_mass.push_back(tree_L0_mass->at(i));
               minitree_L0_pt.push_back(tree_L0_pt->at(i));
         }

      for (unsigned int i = 0 ; i < tree_V0_reco_mass->size(); i++)
         {
               minitree_V0_reco_mass.push_back(tree_V0_reco_mass->at(i));
               minitree_V0_reco_pt.push_back(tree_V0_reco_pt->at(i));
               minitree_V0_reco_source.push_back(tree_V0_reco_source->at(i));
         }

      for (unsigned int i = 0 ; i < tree_SecInt_mass->size() ; i++)
         {
               minitree_SecInt_mass.push_back(tree_SecInt_mass->at(i));
               minitree_SecInt_pt.push_back(tree_SecInt_pt->at(i));
               minitree_SecInt_drSig.push_back(tree_SecInt_drSig->at(i));
               minitree_SecInt_dzSig.push_back(tree_SecInt_dzSig->at(i));
               minitree_SecInt_layer.push_back(tree_SecInt_layer->at(i));
               minitree_SecInt_selec.push_back(tree_SecInt_selec->at(i));
               minitree_SecInt_r.push_back(tree_SecInt_r->at(i));
               minitree_SecInt_z.push_back(tree_SecInt_z->at(i));
         }
      // std::cout<< "Vtx   "<<std::endl;

      for (unsigned int i = 0; i < tree_Hemi_Vtx_step->size(); i++)
         {

               minitree_Hemi_Vtx_step.push_back(tree_Hemi_Vtx_step->at(i));
               minitree_Hemi_Vtx_NChi2.push_back(tree_Hemi_Vtx_NChi2->at(i));
               minitree_Hemi_Vtx_nTrks.push_back(tree_Hemi_Vtx_nTrks->at(i));
               minitree_Hemi_Vtx_dR.push_back(tree_Hemi_Vtx_dR->at(i));
               minitree_Hemi_Vtx_xError.push_back(tree_Hemi_Vtx_xError->at(i));
               minitree_Hemi_Vtx_yError.push_back(tree_Hemi_Vtx_yError->at(i));
               minitree_Hemi_Vtx_zError.push_back(tree_Hemi_Vtx_zError->at(i));
               // minitree_Hemi_Vtx_trackWeight.push_back(tree_Hemi_Vtx_trackWeight->at(i));
               minitree_Hemi_Vtx_SumtrackWeight.push_back(tree_Hemi_Vtx_SumtrackWeight->at(i));
               minitree_Hemi_Vtx_track_MeanDCA_d.push_back(tree_Hemi_Vtx_track_MeanDCA_d->at(i));
               minitree_Hemi_Vtx_Mass.push_back(tree_Hemi_Vtx_Mass->at(i));
               minitree_Hemi_Vtx_dist.push_back(tree_Hemi_Vtx_dist->at(i));
               minitree_Hemi_Vtx_ntrk10.push_back(tree_Hemi_Vtx_ntrk10->at(i));
               minitree_Hemi_Vtx_ntrk20.push_back(tree_Hemi_Vtx_ntrk20->at(i));
         }
      // std::cout<< "SecVtx   "<<std::endl;
      for (unsigned int i = 0; i < tree_Hemi_SecVtx_step->size(); i++)
         {
               minitree_Hemi_SecVtx_step.push_back(tree_Hemi_SecVtx_step->at(i));
               minitree_Hemi_SecVtx_NChi2.push_back(tree_Hemi_SecVtx_NChi2->at(i));
               minitree_Hemi_SecVtx_dR.push_back(tree_Hemi_SecVtx_dR->at(i));
               minitree_Hemi_SecVtx_nTrks.push_back(tree_Hemi_SecVtx_nTrks->at(i));
               minitree_Hemi_SecVtx_dist.push_back(tree_Hemi_SecVtx_dist->at(i));
               minitree_Hemi_SecVtx_track_MeanDCA_d.push_back(tree_Hemi_SecVtx_track_MeanDCA_d->at(i));
               minitree_Hemi_SecVtx_SumtrackWeight.push_back(tree_Hemi_SecVtx_SumtrackWeight->at(i));
               // minitree_Hemi_SecVtx_trackWeight.push_back(tree_Hemi_SecVtx_trackWeight->at(i));
               minitree_Hemi_SecVtx_Mass.push_back(tree_Hemi_SecVtx_Mass->at(i));
         }

      for (unsigned int i = 0; i < tree_Hemi_Vtx_BDT_nTrks->size(); i++)
         {
               minitree_Hemi_Vtx_BDT_nTrks.push_back(tree_Hemi_Vtx_BDT_nTrks->at(i));
               minitree_Hemi_Vtx_BDT_NChi2.push_back(tree_Hemi_Vtx_BDT_NChi2->at(i));
               minitree_Hemi_Vtx_BDT_step.push_back(tree_Hemi_Vtx_BDT_step->at(i));
               minitree_Hemi_Vtx_BDT_Mass.push_back(tree_Hemi_Vtx_BDT_Mass->at(i));
               minitree_Hemi_Vtx_BDT_HMass.push_back(tree_Hemi_Vtx_BDT_HMass->at(i));
         }


      //////////// FILLING //////////

      smalltree->Fill();

      //////////// CLEARING //////////

      minitree_only_gen_wt.clear();

      minitree_genTop_Weight.clear();
      minitree_gen_top_pt.clear();
      minitree_gen_top_rw_pt.clear();
      miniPUweight.clear();
      miniPUweight_Up.clear();
      miniPUweight_Down.clear();
      miniPrefweight.clear();
      minitree_Filter.clear();
      minitree_FilterSameSign.clear();
      minitree_nPV.clear();
      minitree_trigger_doublelepton.clear();
      minitree_trigger_singlelepton.clear();

      minitree_Good_PV.clear();
      minitree_smu_mass.clear();
      minitree_neu_mass.clear();
      minitree_neu_ctau.clear();

      minitree_HT.clear();

      minitree_TRACK_SIZE.clear();
      minitree_nTracks.clear();
      minitree_nLostTracks.clear();

      // minitree_muon_GenRecoTriggerMatched.clear();
      minitree_all_nmu.clear();
      minitree_nmu.clear();    
      minitree_LT.clear();
      minitree_Mmumu.clear();
      minitree_MmumuSameSign.clear();

      minitree_muon_isPrompt.clear();
      minitree_muon_pt.clear();
      minitree_muon_SF.clear();
      minitree_muon_eta.clear();
      minitree_muon_phi.clear();
      minitree_muon_dxy.clear();
      minitree_muon_dz.clear();
      minitree_muon_charge.clear();
      minitree_muon_correction.clear();
      // minitree_muon_gen.clear();
      minitree_muon_dxyError.clear();
      minitree_muon_dzError.clear();
      minitree_muon_isLoose.clear();
      minitree_muon_isMedium.clear();
      minitree_muon_isTight.clear();
      minitree_muon_isGlobal.clear();
      minitree_muon_PFIsoVeryLoose.clear();
      minitree_muon_PFIsoLoose.clear();
      minitree_muon_PFIsoMedium.clear();
      minitree_muon_PFIsoTight.clear();
      minitree_muon_TkIsoLoose.clear();
      minitree_muon_TkIsoTight.clear();
      minitree_muon_MiniIsoLoose.clear();
      minitree_muon_MiniIsoMedium.clear();
      minitree_muon_MiniIsoTight.clear();

      minitree_lepton_leadingpt.clear();
      minitree_lepton_leadingpt2.clear();
      minitree_lepton_leadingeta.clear();
      minitree_lepton_leadingeta2.clear();
      minitree_lepton_leadingphi.clear();
      minitree_lepton_leadingphi2.clear();

      minitree_lepton_lepton_dR.clear();
      minitree_lepton_lepton_dPhi.clear();
      minitree_lepton_lepton_dEta.clear();

      minitree_lepton_leadingdxy.clear();
      minitree_lepton_leadingdxy2.clear();
      minitree_lepton_leadingdz.clear();
      minitree_lepton_leadingdz2.clear();

      // minitree_reco_lepton_leadingpt.clear();
      // minitree_reco_lepton_leadingpt2.clear();
      // minitree_reco_lepton_leadingeta.clear();
      // minitree_reco_lepton_leadingeta2.clear();
      // minitree_reco_lepton_leadingphi.clear();
      // minitree_reco_lepton_leadingphi2.clear();

      // minitree_lepton_b4trigger_leadingpt.clear();
      // minitree_lepton_b4trigger_leadingpt2.clear();

      minitree_all_nel.clear(); 
      minitree_electron_nEle.clear(); 
      minitree_electron_isPrompt.clear();
      minitree_electron_pt.clear();
      minitree_electron_eta.clear();
      minitree_electron_phi.clear();
      minitree_electron_charge.clear();
      minitree_electron_dxy.clear();
      minitree_electron_dz.clear();
      minitree_electron_gen.clear();

      minitree_electron_energy.clear();
      minitree_electron_et.clear();
      minitree_electron_ecal_trk_postcorr.clear();
      minitree_electron_isoR4.clear();
      minitree_electron_IsLoose.clear();
      minitree_electron_IsMedium.clear();
      minitree_electron_IsTight.clear();
   
      minitree_njet.clear();
      minitree_njetNOmu.clear();
      minitree_jet_pt.clear();
      minitree_jet_eta.clear();
      minitree_jet_phi.clear();
      minitree_jet_E.clear();
      minitree_jet_HadronFlavour.clear();
      minitree_jet_btag_DeepJet.clear();

      minitree_jet_leadingpt.clear();
      minitree_jet_leadingpt2.clear();
      minitree_jet_leadingeta.clear();
      minitree_jet_leadingeta2.clear();
      minitree_jet_jet_dR.clear();
      minitree_jet_jet_dPhi.clear();
      minitree_jet_jet_dEta.clear();
      minitree_muon_jet_dRmin.clear();
      minitree_muon_jet_dRmax.clear();
      minitree_elemu_jet_dRmin.clear();
      minitree_elemu_jet_dRmax.clear();
      minitree_ele_jet_dRmin.clear();
      minitree_ele_jet_dRmax.clear();

      minitree_Hemi.clear();
      minitree_Hemi_njet.clear();
      minitree_Hemi_njet_nomu.clear();
      minitree_Hemi_pt.clear();
      minitree_Hemi_eta.clear();
      minitree_Hemi_phi.clear();
      minitree_Hemi_nTrks.clear();
      minitree_Hemi_nTrks_sig.clear();
      minitree_Hemi_nTrks_bad.clear();
      minitree_Hemi_mass.clear();
      minitree_HemiMu_mass.clear();
      minitree_HemiMu_pt.clear();
      minitree_HemiMu_dR.clear();
      minitree_HemiMuOp_mass.clear();
      minitree_HemiMuOp_pt.clear();
      minitree_HemiMuOp_dR.clear();
      minitree_Hemi_dR12.clear();

      // minitree_Evts_MVAval.clear();
      // minitree_Evts_MVAvalDY.clear();
      // minitree_Evts_MVAvalTT.clear();
    
      minitree_ll_pt.clear();
      minitree_ll_eta.clear();
      minitree_ll_phi.clear();
      minitree_ll_px.clear();
      minitree_ll_py.clear();
      minitree_ll_pz.clear();
      minitree_ll_energy.clear();
      minitree_ll_mass.clear();


      minitree_track_ipc.clear();
      minitree_track_lost.clear();
      minitree_track_px.clear();
      minitree_track_py.clear();
      minitree_track_pz.clear();
      minitree_track_pt.clear();
      minitree_track_eta.clear();
      minitree_track_phi.clear();
      minitree_track_charge.clear();
      minitree_track_NChi2.clear();
      minitree_track_isHighPurity.clear();
      minitree_track_dxy.clear();
      minitree_track_dxyError.clear();
      minitree_track_drSig.clear();
      minitree_track_dz.clear();
      minitree_track_dzError.clear();
      minitree_track_dzSig.clear();
      minitree_track_nHit.clear();
      minitree_track_nHitPixel.clear();
      minitree_track_nHitTIB.clear();
      minitree_track_nHitTID.clear();
      minitree_track_nHitTOB.clear();
      minitree_track_nHitTEC.clear();
      minitree_track_nHitPXB.clear();
      minitree_track_nHitPXF.clear();
      minitree_track_isHitPixel.clear();
      minitree_track_nLayers.clear();
      minitree_track_nLayersPixel.clear();
      minitree_track_x.clear();
      minitree_track_y.clear();
      minitree_track_z.clear();
      minitree_track_firstHit.clear();
      minitree_track_region.clear();
      minitree_track_firstHit_x.clear();
      minitree_track_firstHit_y.clear();
      minitree_track_firstHit_z.clear();
      minitree_track_iJet.clear();
      minitree_track_ntrk10.clear();
      minitree_track_ntrk20.clear();
      minitree_track_ntrk30.clear();
      minitree_track_ntrk40.clear();
      minitree_track_MVAval.clear();
      minitree_track_Hemi_dR.clear();
      minitree_track_Hemi_dRmax.clear();

      minitree_K0_mass.clear();
      minitree_K0_pt.clear();
      minitree_L0_mass.clear();
      minitree_L0_pt.clear();
      minitree_V0_reco_mass.clear();
      minitree_V0_reco_pt.clear();
      minitree_V0_reco_source.clear();
      minitree_SecInt_mass.clear();
      minitree_SecInt_pt.clear();
      minitree_SecInt_drSig.clear();
      minitree_SecInt_dzSig.clear();
      minitree_SecInt_layer.clear();
      minitree_SecInt_selec.clear();
      minitree_SecInt_r.clear();
      minitree_SecInt_z.clear();
      minitree_Hemi_Vtx_step.clear();
      minitree_Hemi_Vtx_NChi2.clear();
      minitree_Hemi_Vtx_nTrks.clear();
      minitree_Hemi_Vtx_dR.clear();
      minitree_Hemi_Vtx_xError.clear();
      minitree_Hemi_Vtx_yError.clear();
      minitree_Hemi_Vtx_zError.clear();
      minitree_Hemi_Vtx_trackWeight.clear();
      minitree_Hemi_Vtx_SumtrackWeight.clear();
      minitree_Hemi_Vtx_track_MeanDCA_d.clear();
      minitree_Hemi_Vtx_Mass.clear();
      minitree_Hemi_Vtx_dist.clear();
      minitree_Hemi_Vtx_ntrk10.clear();
      minitree_Hemi_Vtx_ntrk20.clear();
      minitree_Hemi_SecVtx_step.clear();
      minitree_Hemi_SecVtx_NChi2.clear();
      minitree_Hemi_SecVtx_dR.clear();
      minitree_Hemi_SecVtx_nTrks.clear();
      minitree_Hemi_SecVtx_dist.clear();
      minitree_Hemi_SecVtx_track_MeanDCA_d.clear();
      minitree_Hemi_SecVtx_SumtrackWeight.clear();
      minitree_Hemi_SecVtx_trackWeight.clear();
      minitree_Hemi_SecVtx_Mass.clear();
      minitree_Hemi_Vtx_BDT_nTrks.clear();
      minitree_Hemi_Vtx_BDT_NChi2.clear();
      minitree_Hemi_Vtx_BDT_step.clear();
      minitree_Hemi_Vtx_BDT_Mass.clear();
      minitree_Hemi_Vtx_BDT_HMass.clear();
     
   } // end loop on events 

   std::cout << "number of events  "<< allevents << std::endl;


   myFile->Write();
   delete myFile;

   // HistogramManager h ;
   // h.WriteAllHistogramsInFile((Production+"/Mini"+sample+".root").Data(),"recreate");
}

