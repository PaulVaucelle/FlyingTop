#define MiniNtuple_cxx
#include "MiniNtuple.h"
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



void MiniNtuple::Loop(TString sample , TString Production,bool Signal )
{
//   In a ROOT session, you can do:
//      root> .L MiniNtuple.C
//      root> MiniNtuple t
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
   TFile * myFile = new TFile( (Production+"/Mini"+sample+".root").Data(), "recreate");
   TTree *smalltree = new TTree("ttree", "summary information");
  
   TH1F* hData_Good_PV = new TH1F("hData_Good_PV","",2,-0.5,1.5);
   TH1F* hData_Filter = new TH1F("hData_Filter","",2,-0.5,1.5);
   TH1F* hData_FilterSameSign = new TH1F("hData_FilterSameSign","",2,-0.5,1.5);
   TH1F* hData_njetNOmu_Filter = new TH1F("hData_njetNOmu_Filter","",31,-0.5,30.5);
   TH1F* hData_njetNOmu_FilterSameSign = new TH1F("hData_njetNOmu_FilterSameSign","",31,-0.5,30.5);
   TH1F* hData_njetNOmu = new TH1F("hData_njetNOmu","",31,-0.5,30.5);

   TH1F* hData_nPV          = new TH1F("hData_nPV","",81,-0.5,80.5);
   TH1F* hTk_MVA            = new TH1F("hTk_MVA","",200,-1.,1.);
   TH1F* hData_Hemi         = new TH1F("hData_Hemi","",2,0.5,2.5);
   TH1F* hData_Hemi_stepGE1 = new TH1F("hData_Hemi_stepGE1","",2,0.5,2.5);
   TH1F* hData_Hemi_step12  = new TH1F("hData_Hemi_step12","",2,0.5,2.5);

   std::vector<int>    minirunNumber;
   std::vector<int>    minieventNumber;
   std::vector<int>    minilumiBlock;
   std::vector<float>  minitree_LHE_Weights;
   std::vector<float>  minitree_MCEvt_weight;
   std::vector<double> minitree_only_gen_wt;
   std::vector<double> minitree_event_weight; // Gen event wt
   std::vector<double> minitree_genTop_Weight;
   std::vector<double> minitree_gen_top_pt;
   std::vector<double> minitree_gen_top_rw_pt;

   std::vector<double> miniPUweight;
   std::vector<double> miniPUweight_Up;
   std::vector<double> miniPUweight_Down;
   std::vector<double> miniPrefweight;
   std::vector<int>    miniPU_events;
   std::vector<int>    miniAllPU_events_weight;
   std::vector<bool>   minitree_Filter;
   std::vector<bool>   minitree_FilterSameSign;
   std::vector<bool>   minitree_trigger_doublelepton;
   std::vector<bool>   minitree_trigger_singlelepton;
    
   std::vector<float>  minitree_GenPVx;
   std::vector<float>  minitree_GenPVy;
   std::vector<float>  minitree_GenPVz;
   std::vector<int>    minitree_smu_mass; // used for signal only
   std::vector<int>    minitree_neu_mass; // used for signal only
   std::vector<int>    minitree_neu_ctau; // used for signal only
   
   std::vector<bool>   minitree_Good_PV;
   std::vector<int>    minitree_nPV;
   std::vector<float>  minitree_PV_x;
   std::vector<float>  minitree_PV_y;
   std::vector<float>  minitree_PV_z;
   std::vector<float>  minitree_PV_ez;
   std::vector<float>  minitree_PV_NChi2;
   std::vector<int>    minitree_PV_ndf;

   std::vector<float>  minitree_PFMet_et; // in BDT_EVT
   std::vector<float>  minitree_PFMet_phi;

   std::vector<float>  minitree_HT; // in BDT_EVT

   std::vector<int>    minitree_TRACK_SIZE; // in BDT_EVT
   std::vector<int>    minitree_nTracks;
   std::vector<int>    minitree_nLostTracks;

   std::vector<int>    minitree_muon_GenRecoTriggerMatched;
   std::vector<int>    minitree_all_nmu; // count all muons
   std::vector<int>    minitree_nmu;	 // count prompt muons
   std::vector<float>  minitree_LT;
   std::vector<float>  minitree_Mmumu;
   std::vector<float>  minitree_MmumuSameSign;

   std::vector<bool>   minitree_muon_isPrompt;
   std::vector<float>  minitree_muon_pt;
   std::vector<float>  minitree_muon_SF;
   std::vector<float>  minitree_muon_eta;
   std::vector<float>  minitree_muon_phi;
   std::vector<float>  minitree_muon_dxy;
   std::vector<float>  minitree_muon_dz;
   std::vector<int>    minitree_muon_charge;
   std::vector<float>  minitree_muon_correction;
   std::vector<int>    minitree_muon_gen;
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
   std::vector<float>  minitree_ele_jet_dRmin; // empty, usefull ???
   std::vector<float>  minitree_ele_jet_dRmax; // empty, usefull ???

   std::vector<float>  minitree_Evts_MVAval;
   std::vector<float>  minitree_Evts_MVAvalDY;
   std::vector<float>  minitree_Evts_MVAvalTT;

    std::vector<int>   minitree_nLLP;
    std::vector<int>   minitree_LLP;
    std::vector<float> minitree_LLP_pt;
    std::vector<float> minitree_LLP_eta;
    std::vector<float> minitree_LLP_phi;
    std::vector<float> minitree_LLP_x;
    std::vector<float> minitree_LLP_y;
    std::vector<float> minitree_LLP_z;
    std::vector<float> minitree_LLP_r;
    std::vector<float> minitree_LLP_dist;
    std::vector<int>   minitree_LLP_nTrks;
    std::vector<float> minitree_LLP12_dR;
    std::vector<float> minitree_LLP12_deta;
    std::vector<float> minitree_LLP12_dphi;
    std::vector<float> minitree_LLP_Mass;

    std::vector<int>   minitree_Hemi;
    std::vector<int>   minitree_Hemi_njet;
    std::vector<int>   minitree_Hemi_njet_nomu;
    std::vector<float> minitree_Hemi_pt;
    std::vector<float> minitree_Hemi_eta;
    std::vector<float> minitree_Hemi_phi;
    std::vector<int>   minitree_Hemi_nTrks;
    std::vector<int>   minitree_Hemi_nTrks_sig;
    std::vector<int>   minitree_Hemi_nTrks_bad;
    std::vector<float> minitree_Hemi_mass;
    std::vector<float> minitree_HemiMu_mass;
    std::vector<float> minitree_HemiMu_pt;
    std::vector<float> minitree_HemiMu_dR;
    std::vector<float> minitree_HemiMuOp_mass;
    std::vector<float> minitree_HemiMuOp_pt;
    std::vector<float> minitree_HemiMuOp_dR;
    std::vector<float> minitree_Hemi_dR12;

    std::vector<int>   minitree_Hemi_LLP;
    std::vector<float> minitree_Hemi_LLP_pt;
    std::vector<float> minitree_Hemi_LLP_eta;
    std::vector<float> minitree_Hemi_LLP_phi;
    std::vector<float> minitree_Hemi_LLP_dist;
    std::vector<float> minitree_Hemi_LLP_x;
    std::vector<float> minitree_Hemi_LLP_y;
    std::vector<float> minitree_Hemi_LLP_z;
    std::vector<float> minitree_Hemi_LLP_dR;
    std::vector<int>   minitree_Hemi_LLP_mother;
    std::vector<float> minitree_Hemi_LLP_Vtx_dx;
    std::vector<float> minitree_Hemi_LLP_Vtx_dy;
    std::vector<float> minitree_Hemi_LLP_Vtx_dz;
    std::vector<float> minitree_Hemi_LLP_Vtx_dr;
    std::vector<float> minitree_Hemi_LLP_muOK_dR;
    std::vector<float> minitree_Hemi_LLP_muOK_pt;
    std::vector<float> minitree_Hemi_LLP_muOK_mass;
    std::vector<float> minitree_Hemi_LLP_muNO_dR;
    std::vector<float> minitree_Hemi_LLP_muNO_pt;
    std::vector<float> minitree_Hemi_LLP_muNO_mass;
    std::vector<float> minitree_Hemi_LLP_dR12;
    std::vector<bool>  minitree_Hemi_LLP_ping; 
    std::vector<int>   minitree_event_LLP_ping;
    
    std::vector<int>   minitree_Hemi_Vtx_step;
    std::vector<bool>  minitree_Hemi_Vtx_isTight;
    std::vector<float> minitree_Hemi_Vtx_NChi2;
    std::vector<int>   minitree_Hemi_Vtx_nTrks;
    std::vector<int>   minitree_Hemi_Vtx_nTrks_sig;
    std::vector<int>   minitree_Hemi_Vtx_nTrks_bad;
    std::vector<float> minitree_Hemi_Vtx_x;
    std::vector<float> minitree_Hemi_Vtx_y;
    std::vector<float> minitree_Hemi_Vtx_z;
    std::vector<float> minitree_Hemi_Vtx_r;
    std::vector<float> minitree_Hemi_Vtx_dR;
    std::vector<float> minitree_Hemi_Vtx_SumtrackWeight;//Vertx selection variable for the BDT
    std::vector<float> minitree_Hemi_Vtx_Mass;
    std::vector<float> minitree_Hemi_Vtx_track_MeanDCA_d;//Veertex selection BDT
    std::vector<float> minitree_Hemi_Vtx_dist;
    std::vector<int>   minitree_event_nVtx;
    std::vector<float> minitree_event_Vtx_Vtx_dr;
    std::vector<float> minitree_event_Vtx_Vtx_dz;
    std::vector<float> minitree_event_Vtx_Vtx_dd;
    std::vector<float> minitree_event_Vtx_Vtx_reldd;
    std::vector<float> minitree_event_Vtx_Vtx_dR;
    std::vector<int>   minitree_event_Vtx_Vtx_step;

    std::vector<float> minitree_Hemi_SecLLP;
    std::vector<float> minitree_Hemi_LLP_SecVtx_dz;
    std::vector<float> minitree_Hemi_LLP_SecVtx_dr;
    std::vector<bool>  minitree_Hemi_SecLLP_ping;
    std::vector<int>   minitree_event_SecLLP_ping;

    std::vector<int>   minitree_Hemi_SecVtx;      // Hemi (1 or 2) if merging
    std::vector<int>   minitree_Hemi_SecVtx_step; // vertex step for this Hemi if merging
    std::vector<float> minitree_Hemi_SecVtx_x;
    std::vector<float> minitree_Hemi_SecVtx_y;
    std::vector<float> minitree_Hemi_SecVtx_z;
    std::vector<float> minitree_Hemi_SecVtx_r;
    std::vector<float> minitree_Hemi_SecVtx_dR;
    std::vector<float> minitree_Hemi_SecVtx_nTrks;
    std::vector<float> minitree_Hemi_SecVtx_NChi2;
    std::vector<float> minitree_Hemi_SecVtx_dist;
    std::vector<float> minitree_Hemi_SecVtx_track_MeanDCA_d;
    std::vector<float> minitree_Hemi_SecVtx_SumtrackWeight;
    std::vector<float> minitree_Hemi_SecVtx_Mass;
    std::vector<float> minitree_event_MergedVtx_Vtx_dr;
    std::vector<float> minitree_event_MergedVtx_Vtx_dz;
    std::vector<float> minitree_event_MergedVtx_Vtx_dd;
    std::vector<float> minitree_event_MergedVtx_Vtx_reldd;
    std::vector<float> minitree_event_MergedVtx_Vtx_dR;
    std::vector<int>   minitree_event_MergedVtx_Vtx_step;
    
    std::vector<float> minitree_Hemi_Vtx_BDT_nTrks;
    std::vector<float> minitree_Hemi_Vtx_BDT_NChi2;
    std::vector<float> minitree_Hemi_Vtx_BDT_step;
    std::vector<float> minitree_Hemi_Vtx_BDT_STW;
    std::vector<float> minitree_Hemi_Vtx_BDT_Mass;
    std::vector<float> minitree_Hemi_Vtx_BDT_HMass;
    std::vector<float> minitree_Hemi_Vtx_BDT_ntrk10;
    std::vector<float> minitree_Hemi_Vtx_BDT_ntrk20;
    std::vector<float> minitree_Hemi_Vtx_BDT_MeanDCA;
    std::vector<float> minitree_Hemi_Vtx_MVAval_Loose;
    std::vector<float> minitree_Hemi_Vtx_MVAval_Tight;//TIght WP
    
   std::vector<bool>  miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;
   std::vector<bool>  miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v;
   std::vector<bool>  miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v;
   std::vector<bool>  miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   std::vector<bool>  miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;
   std::vector<bool>  miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   std::vector<bool>  miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   std::vector<bool>  miniHLT_Ele27_WPTight_Gsf_v;
   std::vector<bool>  miniHLT_Ele32_WPTight_Gsf_v;
   std::vector<bool>  miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   std::vector<bool>  miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   std::vector<bool>  miniHLT_PFMET120_PFMHT120_IDTight_v;
   std::vector<bool>  miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v;
   std::vector<bool>  miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;
   std::vector<bool>  miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;
   std::vector<bool>  miniHLT_PFMET250_HBHECleaned_v;
   std::vector<bool>  miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;
   std::vector<bool>  miniHLT_IsoMu24_v;
   std::vector<bool>  miniHLT_IsoMu27_v;
   std::vector<bool>  miniHLT_IsoTkMu24_v; 
   std::vector<bool>  miniHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v;  
   std::vector<bool> miniHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;    
   std::vector<bool> miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v;    
   std::vector<bool> miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v;  
   std::vector<bool> miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
 

     
   smalltree->Branch("minirunNumber",&minirunNumber);
   smalltree->Branch("minieventNumber",&minieventNumber);
   smalltree->Branch("minilumiBlock",&minilumiBlock);
   smalltree->Branch("minitree_LHE_Weights",&minitree_LHE_Weights);
   smalltree->Branch("minitree_MCEvt_weight",&minitree_MCEvt_weight);
   smalltree->Branch("minitree_only_gen_wt",&minitree_only_gen_wt);
   smalltree->Branch("minitree_event_weight",&minitree_event_weight);
   smalltree->Branch("minitree_genTop_Weight",&minitree_genTop_Weight);
   smalltree->Branch("minitree_gen_top_pt",&minitree_gen_top_pt);
   smalltree->Branch("minitree_gen_top_rw_pt",&minitree_gen_top_rw_pt);

   smalltree->Branch("miniPUweight",&miniPUweight);
   smalltree->Branch("miniPUweight_Up",&miniPUweight_Up);
   smalltree->Branch("miniPUweight_Down",&miniPUweight_Down);
   smalltree->Branch("miniPrefweight",&miniPrefweight);
   smalltree->Branch("miniPU_events",&miniPU_events);
   smalltree->Branch("miniAllPU_events_weight",&miniAllPU_events_weight);
   smalltree->Branch("minitree_Filter",        &minitree_Filter);
   smalltree->Branch("minitree_FilterSameSign",&minitree_FilterSameSign);
   smalltree->Branch("minitree_trigger_doublelepton",&minitree_trigger_doublelepton);
   smalltree->Branch("minitree_trigger_singlelepton",&minitree_trigger_singlelepton);

   smalltree->Branch("minitree_GenPVx",&minitree_GenPVx);
   smalltree->Branch("minitree_GenPVy",&minitree_GenPVy);
   smalltree->Branch("minitree_GenPVz",&minitree_GenPVz);
   smalltree->Branch("minitree_smu_mass",&minitree_smu_mass); // used for signal only
   smalltree->Branch("minitree_neu_mass",&minitree_neu_mass); // used for signal only
   smalltree->Branch("minitree_neu_ctau",&minitree_neu_ctau); // used for signal only

   smalltree->Branch("minitree_Good_PV",&minitree_Good_PV);
   smalltree->Branch("minitree_nPV",&minitree_nPV);
   smalltree->Branch("minitree_PV_x",&minitree_PV_x);
   smalltree->Branch("minitree_PV_y",&minitree_PV_y);
   smalltree->Branch("minitree_PV_z",&minitree_PV_z);
   smalltree->Branch("minitree_PV_ez",&minitree_PV_ez);
   smalltree->Branch("minitree_PV_NChi2",&minitree_PV_NChi2);
   smalltree->Branch("minitree_PV_ndf",&minitree_PV_ndf);

   smalltree->Branch("minitree_PFMet_et",&minitree_PFMet_et); // in BDT_EVT
   smalltree->Branch("minitree_PFMet_phi",&minitree_PFMet_phi);

   smalltree->Branch("minitree_HT",&minitree_HT); // in BDT_EVT

   smalltree->Branch("minitree_TRACK_SIZE",&minitree_TRACK_SIZE); // in BDT_EVT
   smalltree->Branch("minitree_nTracks",&minitree_nTracks);
   smalltree->Branch("minitree_nLostTracks",&minitree_nLostTracks);

   smalltree->Branch("minitree_muon_GenRecoTriggerMatched",&minitree_muon_GenRecoTriggerMatched);
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
   smalltree->Branch("minitree_muon_gen",&minitree_muon_gen);
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

   smalltree->Branch("minitree_Evts_MVAval",   &minitree_Evts_MVAval);
   smalltree->Branch("minitree_Evts_MVAvalDY",   &minitree_Evts_MVAvalDY);
   smalltree->Branch("minitree_Evts_MVAvalTT",   &minitree_Evts_MVAvalTT);

   smalltree->Branch("minitree_nLLP",          &minitree_nLLP);
   smalltree->Branch("minitree_LLP",           &minitree_LLP);
   smalltree->Branch("minitree_LLP_pt" ,       &minitree_LLP_pt);
   smalltree->Branch("minitree_LLP_eta",       &minitree_LLP_eta);
   smalltree->Branch("minitree_LLP_phi",       &minitree_LLP_phi);
   smalltree->Branch("minitree_LLP_x",         &minitree_LLP_x);
   smalltree->Branch("minitree_LLP_y",         &minitree_LLP_y);
   smalltree->Branch("minitree_LLP_z",         &minitree_LLP_z);
   smalltree->Branch("minitree_LLP_r",         &minitree_LLP_r);
   smalltree->Branch("minitree_LLP_dist",      &minitree_LLP_dist);
   smalltree->Branch("minitree_LLP_nTrks",     &minitree_LLP_nTrks);
   smalltree->Branch("minitree_LLP12_dR",      &minitree_LLP12_dR);
   smalltree->Branch("minitree_LLP12_deta",    &minitree_LLP12_deta);
   smalltree->Branch("minitree_LLP12_dphi",    &minitree_LLP12_dphi);
   smalltree->Branch("minitree_LLP_Mass",      &minitree_LLP_Mass);

   smalltree->Branch("minitree_Hemi",       &minitree_Hemi);
   smalltree->Branch("minitree_Hemi_njet",  &minitree_Hemi_njet);
   smalltree->Branch("minitree_Hemi_njet_nomu",  &minitree_Hemi_njet_nomu);
   smalltree->Branch("minitree_Hemi_pt",    &minitree_Hemi_pt);
   smalltree->Branch("minitree_Hemi_eta",   &minitree_Hemi_eta);
   smalltree->Branch("minitree_Hemi_phi",   &minitree_Hemi_phi);
   smalltree->Branch("minitree_Hemi_nTrks", &minitree_Hemi_nTrks);
   smalltree->Branch("minitree_Hemi_nTrks_sig", &minitree_Hemi_nTrks_sig);
   smalltree->Branch("minitree_Hemi_nTrks_bad", &minitree_Hemi_nTrks_bad);
   smalltree->Branch("minitree_Hemi_mass",     &minitree_Hemi_mass);
   smalltree->Branch("minitree_HemiMu_mass",   &minitree_HemiMu_mass);
   smalltree->Branch("minitree_HemiMu_pt",     &minitree_HemiMu_pt);
   smalltree->Branch("minitree_HemiMu_dR",     &minitree_HemiMu_dR);
   smalltree->Branch("minitree_HemiMuOp_mass", &minitree_HemiMuOp_mass);
   smalltree->Branch("minitree_HemiMuOp_pt",   &minitree_HemiMuOp_pt);
   smalltree->Branch("minitree_HemiMuOp_dR",   &minitree_HemiMuOp_dR);

   smalltree->Branch("minitree_Hemi_dR12",      &minitree_Hemi_dR12);

   smalltree->Branch("minitree_Hemi_LLP",       &minitree_Hemi_LLP);
   smalltree->Branch("minitree_Hemi_LLP_pt",    &minitree_Hemi_LLP_pt);
   smalltree->Branch("minitree_Hemi_LLP_eta",   &minitree_Hemi_LLP_eta);
   smalltree->Branch("minitree_Hemi_LLP_phi",   &minitree_Hemi_LLP_phi);
   smalltree->Branch("minitree_Hemi_LLP_dist",  &minitree_Hemi_LLP_dist);
   smalltree->Branch("minitree_Hemi_LLP_x",     &minitree_Hemi_LLP_x);
   smalltree->Branch("minitree_Hemi_LLP_y",     &minitree_Hemi_LLP_y);
   smalltree->Branch("minitree_Hemi_LLP_z",     &minitree_Hemi_LLP_z);
   smalltree->Branch("minitree_Hemi_LLP_dR",    &minitree_Hemi_LLP_dR);
   smalltree->Branch("minitree_Hemi_LLP_mother",&minitree_Hemi_LLP_mother);
   smalltree->Branch("minitree_Hemi_LLP_Vtx_dx",     &minitree_Hemi_LLP_Vtx_dx);
   smalltree->Branch("minitree_Hemi_LLP_Vtx_dy",     &minitree_Hemi_LLP_Vtx_dy);
   smalltree->Branch("minitree_Hemi_LLP_Vtx_dz",     &minitree_Hemi_LLP_Vtx_dz);
   smalltree->Branch("minitree_Hemi_LLP_Vtx_dr",     &minitree_Hemi_LLP_Vtx_dr);
   smalltree->Branch("minitree_Hemi_LLP_muOK_dR",    &minitree_Hemi_LLP_muOK_dR);
   smalltree->Branch("minitree_Hemi_LLP_muOK_pt",    &minitree_Hemi_LLP_muOK_pt);
   smalltree->Branch("minitree_Hemi_LLP_muOK_mass",  &minitree_Hemi_LLP_muOK_mass);
   smalltree->Branch("minitree_Hemi_LLP_muNO_dR",    &minitree_Hemi_LLP_muNO_dR);
   smalltree->Branch("minitree_Hemi_LLP_muNO_pt",    &minitree_Hemi_LLP_muNO_pt);
   smalltree->Branch("minitree_Hemi_LLP_muNO_mass",  &minitree_Hemi_LLP_muNO_mass);
   smalltree->Branch("minitree_Hemi_LLP_dR12",  &minitree_Hemi_LLP_dR12);
   smalltree->Branch("minitree_Hemi_LLP_ping",  &minitree_Hemi_LLP_ping);
   smalltree->Branch("minitree_event_LLP_ping", &minitree_event_LLP_ping);
    
   smalltree->Branch("minitree_Hemi_Vtx_step",  &minitree_Hemi_Vtx_step);
   smalltree->Branch("minitree_Hemi_Vtx_isTight",&minitree_Hemi_Vtx_isTight);
   smalltree->Branch("minitree_Hemi_Vtx_NChi2", &minitree_Hemi_Vtx_NChi2);
   smalltree->Branch("minitree_Hemi_Vtx_nTrks", &minitree_Hemi_Vtx_nTrks);
   smalltree->Branch("minitree_Hemi_Vtx_nTrks_sig", &minitree_Hemi_Vtx_nTrks_sig);
   smalltree->Branch("minitree_Hemi_Vtx_nTrks_bad", &minitree_Hemi_Vtx_nTrks_bad);
   smalltree->Branch("minitree_Hemi_Vtx_x",     &minitree_Hemi_Vtx_x);
   smalltree->Branch("minitree_Hemi_Vtx_y",     &minitree_Hemi_Vtx_y);
   smalltree->Branch("minitree_Hemi_Vtx_z",     &minitree_Hemi_Vtx_z);
   smalltree->Branch("minitree_Hemi_Vtx_r",     &minitree_Hemi_Vtx_r);
   smalltree->Branch("minitree_Hemi_Vtx_dR",    &minitree_Hemi_Vtx_dR);
   smalltree->Branch("minitree_Hemi_Vtx_SumtrackWeight",&minitree_Hemi_Vtx_SumtrackWeight);
   smalltree->Branch("minitree_Hemi_Vtx_track_MeanDCA_d",&minitree_Hemi_Vtx_track_MeanDCA_d);
   smalltree->Branch("minitree_Hemi_Vtx_Mass", &minitree_Hemi_Vtx_Mass);
   smalltree->Branch("minitree_Hemi_Vtx_dist",  &minitree_Hemi_Vtx_dist);
   smalltree->Branch("minitree_event_nVtx",      &minitree_event_nVtx);
   smalltree->Branch("minitree_event_Vtx_Vtx_dr",&minitree_event_Vtx_Vtx_dr);
   smalltree->Branch("minitree_event_Vtx_Vtx_dz",&minitree_event_Vtx_Vtx_dz);
   smalltree->Branch("minitree_event_Vtx_Vtx_dd",&minitree_event_Vtx_Vtx_dd);
   smalltree->Branch("minitree_event_Vtx_Vtx_reldd",&minitree_event_Vtx_Vtx_reldd);
   smalltree->Branch("minitree_event_Vtx_Vtx_dR",&minitree_event_Vtx_Vtx_dR);
   smalltree->Branch("minitree_event_Vtx_Vtx_step",&minitree_event_Vtx_Vtx_step);

   smalltree->Branch("minitree_Hemi_SecLLP",&minitree_Hemi_SecLLP);
   smalltree->Branch("minitree_Hemi_LLP_SecVtx_dz",&minitree_Hemi_LLP_SecVtx_dz);
   smalltree->Branch("minitree_Hemi_LLP_SecVtx_dr",&minitree_Hemi_LLP_SecVtx_dr);
   smalltree->Branch("minitree_Hemi_SecLLP_ping",&minitree_Hemi_SecLLP_ping);
   smalltree->Branch("minitree_event_SecLLP_ping",&minitree_event_SecLLP_ping);

   smalltree->Branch("minitree_Hemi_SecVtx",     &minitree_Hemi_SecVtx);
   smalltree->Branch("minitree_Hemi_SecVtx_step",&minitree_Hemi_SecVtx_step);
   smalltree->Branch("minitree_Hemi_SecVtx_x",&minitree_Hemi_SecVtx_x);
   smalltree->Branch("minitree_Hemi_SecVtx_y",&minitree_Hemi_SecVtx_y);
   smalltree->Branch("minitree_Hemi_SecVtx_z",&minitree_Hemi_SecVtx_z);
   smalltree->Branch("minitree_Hemi_SecVtx_r",&minitree_Hemi_SecVtx_r);
   smalltree->Branch("minitree_Hemi_SecVtx_dR",&minitree_Hemi_SecVtx_dR);
   smalltree->Branch("minitree_Hemi_SecVtx_nTrks",&minitree_Hemi_SecVtx_nTrks);
   smalltree->Branch("minitree_Hemi_SecVtx_NChi2", &minitree_Hemi_SecVtx_NChi2);
   smalltree->Branch("minitree_Hemi_SecVtx_dist",&minitree_Hemi_SecVtx_dist);
   smalltree->Branch("minitree_Hemi_SecVtx_track_MeanDCA_d",&minitree_Hemi_SecVtx_track_MeanDCA_d);
   smalltree->Branch("minitree_Hemi_SecVtx_SumtrackWeight",&minitree_Hemi_SecVtx_SumtrackWeight);
   smalltree->Branch("minitree_Hemi_SecVtx_Mass",&minitree_Hemi_SecVtx_Mass);
   smalltree->Branch("minitree_event_MergedVtx_Vtx_dr",&minitree_event_MergedVtx_Vtx_dr);
   smalltree->Branch("minitree_event_MergedVtx_Vtx_dz",&minitree_event_MergedVtx_Vtx_dz);
   smalltree->Branch("minitree_event_MergedVtx_Vtx_dd",&minitree_event_MergedVtx_Vtx_dd);
   smalltree->Branch("minitree_event_MergedVtx_Vtx_reldd",&minitree_event_MergedVtx_Vtx_reldd);
   smalltree->Branch("minitree_event_MergedVtx_Vtx_dR",&minitree_event_MergedVtx_Vtx_dR);
   smalltree->Branch("minitree_event_MergedVtx_Vtx_step",&minitree_event_MergedVtx_Vtx_step);

   smalltree->Branch("minitree_Hemi_Vtx_BDT_nTrks",&minitree_Hemi_Vtx_BDT_nTrks);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_NChi2",&minitree_Hemi_Vtx_BDT_NChi2);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_step",&minitree_Hemi_Vtx_BDT_step);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_STW",&minitree_Hemi_Vtx_BDT_STW);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_Mass",&minitree_Hemi_Vtx_BDT_Mass);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_HMass",&minitree_Hemi_Vtx_BDT_HMass);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_ntrk10",&minitree_Hemi_Vtx_BDT_ntrk10);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_ntrk20",&minitree_Hemi_Vtx_BDT_ntrk20);
   smalltree->Branch("minitree_Hemi_Vtx_BDT_MeanDCA",&minitree_Hemi_Vtx_BDT_MeanDCA);
   smalltree->Branch("minitree_Hemi_Vtx_MVAval_Loose", &minitree_Hemi_Vtx_MVAval_Loose);
   smalltree->Branch("minitree_Hemi_Vtx_MVAval_Tight",&minitree_Hemi_Vtx_MVAval_Tight);

   smalltree->Branch("miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",&miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
   smalltree->Branch("miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",&miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v);
   smalltree->Branch("miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",&miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v);
   smalltree->Branch("miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   smalltree->Branch("miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",&miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
   smalltree->Branch("miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",&miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   smalltree->Branch("miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   smalltree->Branch("miniHLT_Ele27_WPTight_Gsf_v",&miniHLT_Ele27_WPTight_Gsf_v);
   smalltree->Branch("miniHLT_Ele32_WPTight_Gsf_v",&miniHLT_Ele32_WPTight_Gsf_v);
   smalltree->Branch("miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
   smalltree->Branch("miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",&miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   smalltree->Branch("miniHLT_PFMET120_PFMHT120_IDTight_v",&miniHLT_PFMET120_PFMHT120_IDTight_v);
   smalltree->Branch("miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v",&miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v);
   smalltree->Branch("miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",&miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v);
   smalltree->Branch("miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",&miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
   smalltree->Branch("miniHLT_PFMET250_HBHECleaned_v",&miniHLT_PFMET250_HBHECleaned_v);
   smalltree->Branch("miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",&miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v);
   smalltree->Branch("miniHLT_IsoMu24_v",&miniHLT_IsoMu24_v);
   smalltree->Branch("miniHLT_IsoMu27_v",&miniHLT_IsoMu27_v);
   smalltree->Branch("miniHLT_IsoTkMu24_v",&miniHLT_IsoTkMu24_v); 
   smalltree->Branch("miniHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",&miniHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v); 
   smalltree->Branch("miniHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",&miniHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);    
   smalltree->Branch("miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",&miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v);      
   smalltree->Branch("miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v",&miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v);   
   smalltree->Branch("miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);  


   std::cout<<"//-------------------------//"<<std::endl;
   std::cout<<"// "<<Production<<" : "<<sample<<"  //"<<std::endl;
   std::cout<<"//-------------------------//"<<std::endl;

   bool debug = false;
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "Total Entries : " << nentries << std::endl;

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

      hData_Good_PV->Fill( tree_Good_PV );
      hData_Filter->Fill( tree_Filter );
      hData_FilterSameSign->Fill( tree_FilterSameSign );
      if ( tree_Filter ) hData_njetNOmu_Filter->Fill( tree_njetNOmu );
      if ( tree_FilterSameSign ) hData_njetNOmu_FilterSameSign->Fill( tree_njetNOmu );

//$$

    // -------------- --------------------------------------------------------------------//
    // -----------------------------------------------------------------------------------//
    if ( !((tree_Filter || tree_FilterSameSign) && tree_njetNOmu>0) && !Signal ) continue;
    //------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------//
//$$

      minitree_trigger_doublelepton.push_back(tree_trigger_doublelepton);
      minitree_trigger_singlelepton.push_back(tree_trigger_singlelepton);

      minitree_Good_PV.push_back(tree_Good_PV);


      minitree_Filter.push_back(tree_Filter);
      minitree_FilterSameSign.push_back(tree_FilterSameSign);
      if ( tree_Filter && tree_njetNOmu > 0 ) hData_nPV->Fill( tree_nPV );

      miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v.push_back(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
      miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v.push_back(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v);
      miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v.push_back(HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v);
      miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.push_back(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
      miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.push_back(HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
      miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.push_back(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
      miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.push_back(HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
      miniHLT_Ele27_WPTight_Gsf_v.push_back(HLT_Ele27_WPTight_Gsf_v);
      miniHLT_Ele32_WPTight_Gsf_v.push_back(HLT_Ele32_WPTight_Gsf_v);
      miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.push_back(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
      miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v.push_back(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);
      miniHLT_PFMET120_PFMHT120_IDTight_v.push_back(HLT_PFMET120_PFMHT120_IDTight_v);
      miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v.push_back(HLT_PFMET120_PFMHT120_IDTight_PFHT60_v);
      miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v.push_back(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v);
      miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v.push_back(HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
      miniHLT_PFMET250_HBHECleaned_v.push_back(HLT_PFMET250_HBHECleaned_v);
      miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v.push_back(HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v);
      miniHLT_IsoMu24_v.push_back(HLT_IsoMu24_v);
      miniHLT_IsoMu27_v.push_back(HLT_IsoMu27_v);
      miniHLT_IsoTkMu24_v.push_back(HLT_IsoTkMu24_v); 
      miniHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.push_back(HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v);   
      miniHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.push_back(HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);    
      miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.push_back(HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v);    
      miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v.push_back(HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v);  
      miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.push_back(HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);


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
            minitree_muon_gen.push_back(tree_muon_gen->at(i));
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
            // minitree_electron_ecal_trk_postcorr.push_back(tree_electron_ecal_trk_postcorr->at(i));
            minitree_electron_isoR4.push_back(tree_electron_isoR4->at(i));
            minitree_electron_IsLoose.push_back(tree_electron_IsLoose->at(i));
            minitree_electron_IsMedium.push_back(tree_electron_IsMedium->at(i));
            minitree_electron_IsTight.push_back(tree_electron_IsTight->at(i));
         }

      hData_njetNOmu->Fill( tree_njetNOmu );

      if (debug)  {std::cout<<"Weights  "<<std::endl;}
      minirunNumber.push_back(runNumber);
      minieventNumber.push_back(eventNumber);
      minilumiBlock.push_back(lumiBlock);
      // for (unsigned int i = 0; i < tree_LHE_Weights->size(); i++)
      //    {
      //       minitree_LHE_Weights.push_back(tree_LHE_Weights->at(i));
      //    }
      
      // minitree_MCEvt_weight.push_back(tree_MCEvt_weight);

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
      miniPU_events.push_back(PU_events);
      // miniAllPU_events_weight.push_back(AllPU_events_weight);

      if (Signal)
         {
            if (debug)  {std::cout<<"Signal parameters "<<std::endl;}
            minitree_smu_mass.push_back(tree_smu_mass); // used for signal only
            minitree_neu_mass.push_back(tree_neu_mass); // used for signal only
            minitree_neu_ctau.push_back(tree_neu_ctau); // used for signal only
         }

      minitree_muon_GenRecoTriggerMatched.push_back(tree_muon_GenRecoTriggerMatched);

      if (debug)  {std::cout<<"PV (reco and Gen)  "<<std::endl;}
      minitree_GenPVx.push_back(tree_GenPVx);
      minitree_GenPVy.push_back(tree_GenPVy);
      minitree_GenPVz.push_back(tree_GenPVz);

      minitree_nPV.push_back(tree_nPV);
      minitree_PV_x.push_back(tree_PV_x);
      minitree_PV_y.push_back(tree_PV_y);
      minitree_PV_z.push_back(tree_PV_z);
      minitree_PV_ez.push_back(tree_PV_ez);
      minitree_PV_NChi2.push_back(tree_PV_NChi2);
      minitree_PV_ndf.push_back(tree_PV_ndf);

      if (debug)  {std::cout<<"Event Level info  "<<std::endl;}
      minitree_PFMet_et.push_back(tree_PFMet_et); // in BDT_EVT
      minitree_PFMet_phi.push_back(tree_PFMet_phi);

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


      minitree_Evts_MVAval.push_back(tree_Evts_MVAval);
      minitree_Evts_MVAvalDY.push_back(tree_Evts_MVAvalDY);
      minitree_Evts_MVAvalTT.push_back(tree_Evts_MVAvalTT);
      minitree_nLLP.push_back(tree_nLLP);
      minitree_event_weight.push_back(tree_event_weight);

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

      // sALUT DADA, JE SUIS TON ESPRIT. NOUBLIE PAS DE BIEN PRENDRE SOIN DE TON PAUL, CET ETRE PRECIEUX
      // depends on the decay channel :(
      // minitree_elemu_jet_dRmin.push_back(tree_elemu_jet_dRmin->at(0) );
      // minitree_elemu_jet_dRmax.push_back( tree_elemu_jet_dRmax->at(0));
      // minitree_ele_jet_dRmin.push_back( tree_ele_jet_dRmin->at(0)); // empty, usefull ???
      // minitree_ele_jet_dRmax.push_back( tree_ele_jet_dRmax->at(0)); // empty, usefull ???

      if (debug)  {std::cout<<"Tracks "<<std::endl;}
      if ( tree_Filter && tree_njetNOmu>0 ) {
        for (unsigned int i = 0; i < tree_track_pt->size(); i++)
         {
            hTk_MVA->Fill( tree_track_MVAval->at(i) );
         }
      }

      // std::cout<<tree_Filter.capacity()<<std::endl;
      if (debug)  {std::cout<<"Event"<<std::endl;}
      for (unsigned int i = 0; i < tree_event_nVtx->size(); i++)
         {
            minitree_event_nVtx.push_back( tree_event_nVtx->at(i));
            if (Signal)
               {
                  minitree_event_LLP_ping.push_back( tree_event_LLP_ping->at(i));
                  minitree_LLP12_dR.push_back(tree_LLP12_dR->at(i));                        
                  minitree_LLP12_deta.push_back(tree_LLP12_deta->at(i));                      
                  minitree_LLP12_dphi.push_back(tree_LLP12_dphi->at(i)); 
               }

         }

      for (unsigned int i = 0; i < tree_event_MergedVtx_Vtx_step->size(); i++)
         {
            minitree_event_MergedVtx_Vtx_step.push_back( tree_event_MergedVtx_Vtx_step->at(i));
         }
         

      for (unsigned int i = 0; i < tree_event_Vtx_Vtx_dr->size(); i++)
         {
            minitree_event_Vtx_Vtx_dr.push_back( tree_event_Vtx_Vtx_dr->at(i));
            minitree_event_Vtx_Vtx_dz.push_back( tree_event_Vtx_Vtx_dz->at(i));
            minitree_event_Vtx_Vtx_dd.push_back( tree_event_Vtx_Vtx_dd->at(i));
            minitree_event_Vtx_Vtx_reldd.push_back( tree_event_Vtx_Vtx_reldd->at(i));
            minitree_event_Vtx_Vtx_dR.push_back( tree_event_Vtx_Vtx_dR->at(i));
            minitree_event_Vtx_Vtx_step.push_back( tree_event_Vtx_Vtx_step->at(i));
         }
      for (unsigned int i = 0; i < tree_event_SecLLP_ping->size(); i++)
         {
            minitree_event_SecLLP_ping.push_back( tree_event_SecLLP_ping->at(i));
         }

      for (unsigned int i = 0; i < tree_event_MergedVtx_Vtx_dr->size(); i++)
         {
            minitree_event_MergedVtx_Vtx_dr.push_back( tree_event_MergedVtx_Vtx_dr->at(i));
            minitree_event_MergedVtx_Vtx_dz.push_back( tree_event_MergedVtx_Vtx_dz->at(i));
            minitree_event_MergedVtx_Vtx_dd.push_back( tree_event_MergedVtx_Vtx_dd->at(i));
            minitree_event_MergedVtx_Vtx_reldd.push_back( tree_event_MergedVtx_Vtx_reldd->at(i));
            minitree_event_MergedVtx_Vtx_dR.push_back( tree_event_MergedVtx_Vtx_dR->at(i));
         }

      //----------------------//
      //         LLP          //
      //----------------------//
      if (debug)  {std::cout<<"LLP"<<std::endl;}
      for (unsigned int i = 0; i < tree_Hemi->size(); i++)
         {
            if (Signal)
               {
                  minitree_LLP.push_back(tree_LLP->at(i));                          
                  minitree_LLP_pt.push_back(tree_LLP_pt->at(i));                      
                  minitree_LLP_eta.push_back(tree_LLP_eta->at(i));                  
                  minitree_LLP_phi.push_back(tree_LLP_phi->at(i));                      
                  minitree_LLP_x.push_back(tree_LLP_x->at(i));                        
                  minitree_LLP_y.push_back(tree_LLP_y->at(i));                           
                  minitree_LLP_z.push_back(tree_LLP_z->at(i));                           
                  minitree_LLP_r.push_back(tree_LLP_r->at(i));                         
                  minitree_LLP_dist.push_back(tree_LLP_dist->at(i));                       
                  minitree_LLP_nTrks.push_back(tree_LLP_nTrks->at(i));                                         
                  minitree_LLP_Mass.push_back(tree_LLP_Mass->at(i));   


                     //----------------------//
                     //         Hemi         //
                     //----------------------//
                  minitree_Hemi_LLP.push_back(tree_Hemi_LLP->at(i));
                  minitree_Hemi_LLP_pt.push_back(tree_Hemi_LLP_pt->at(i));
                  minitree_Hemi_LLP_eta.push_back(tree_Hemi_LLP_eta->at(i));
                  minitree_Hemi_LLP_phi.push_back(tree_Hemi_LLP_phi->at(i));
                  minitree_Hemi_LLP_dist.push_back(tree_Hemi_LLP_dist->at(i));
                  minitree_Hemi_LLP_x.push_back(tree_Hemi_LLP_x->at(i));
                  minitree_Hemi_LLP_y.push_back(tree_Hemi_LLP_y->at(i));
                  minitree_Hemi_LLP_z.push_back(tree_Hemi_LLP_z->at(i));
                  minitree_Hemi_LLP_dR.push_back(tree_Hemi_LLP_dR->at(i));
                  minitree_Hemi_LLP_mother.push_back(tree_Hemi_LLP_mother->at(i));
                  minitree_Hemi_LLP_Vtx_dx.push_back(tree_Hemi_LLP_Vtx_dx->at(i));
                  minitree_Hemi_LLP_Vtx_dy.push_back(tree_Hemi_LLP_Vtx_dy->at(i));
                  minitree_Hemi_LLP_Vtx_dz.push_back(tree_Hemi_LLP_Vtx_dz->at(i));
                  minitree_Hemi_LLP_Vtx_dr.push_back(tree_Hemi_LLP_Vtx_dr->at(i));
                  minitree_Hemi_LLP_muOK_dR.push_back(tree_Hemi_LLP_dR->at(i));
                  minitree_Hemi_LLP_muOK_pt.push_back(tree_Hemi_LLP_pt->at(i));
                  minitree_Hemi_LLP_muOK_mass.push_back(tree_Hemi_LLP_muOK_mass->at(i));
                  minitree_Hemi_LLP_muNO_dR.push_back(tree_Hemi_LLP_muNO_dR->at(i));
                  minitree_Hemi_LLP_muNO_pt.push_back(tree_Hemi_LLP_muNO_pt->at(i));
                  minitree_Hemi_LLP_muNO_mass.push_back(tree_Hemi_LLP_muNO_mass->at(i));
                  minitree_Hemi_LLP_dR12.push_back(tree_Hemi_LLP_dR12->at(i));
                  minitree_Hemi_LLP_ping.push_back(tree_Hemi_LLP_ping->at(i));
               }


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
               //----------------------//
               //         Vtx          //
               //----------------------//
            minitree_Hemi_Vtx_step.push_back(tree_Hemi_Vtx_step->at(i));
            minitree_Hemi_Vtx_isTight.push_back(tree_Hemi_Vtx_isTight->at(i));
            minitree_Hemi_Vtx_NChi2.push_back(tree_Hemi_Vtx_NChi2->at(i));
            minitree_Hemi_Vtx_nTrks.push_back(tree_Hemi_Vtx_nTrks->at(i));
            minitree_Hemi_Vtx_nTrks_sig.push_back(tree_Hemi_Vtx_nTrks_sig->at(i));
            minitree_Hemi_Vtx_nTrks_bad.push_back(tree_Hemi_Vtx_nTrks_bad->at(i));
            minitree_Hemi_Vtx_x.push_back(tree_Hemi_Vtx_x->at(i));
            minitree_Hemi_Vtx_y.push_back(tree_Hemi_Vtx_y->at(i));
            minitree_Hemi_Vtx_z.push_back(tree_Hemi_Vtx_z->at(i));
            minitree_Hemi_Vtx_r.push_back(tree_Hemi_Vtx_r->at(i));
            minitree_Hemi_Vtx_dR.push_back(tree_Hemi_Vtx_dR->at(i));
            minitree_Hemi_Vtx_SumtrackWeight.push_back(tree_Hemi_Vtx_SumtrackWeight->at(i));
            minitree_Hemi_Vtx_track_MeanDCA_d.push_back(tree_Hemi_Vtx_track_MeanDCA_d->at(i));
            minitree_Hemi_Vtx_Mass.push_back(tree_Hemi_Vtx_Mass->at(i));
            minitree_Hemi_Vtx_dist.push_back(tree_Hemi_Vtx_dist->at(i));

            minitree_Hemi_Vtx_BDT_nTrks.push_back(tree_Hemi_Vtx_BDT_nTrks->at(i));
            minitree_Hemi_Vtx_BDT_NChi2.push_back(tree_Hemi_Vtx_BDT_NChi2->at(i));
            minitree_Hemi_Vtx_BDT_step.push_back(tree_Hemi_Vtx_BDT_step->at(i));
            minitree_Hemi_Vtx_BDT_STW.push_back(tree_Hemi_Vtx_BDT_STW->at(i));
            minitree_Hemi_Vtx_BDT_Mass.push_back(tree_Hemi_Vtx_BDT_Mass->at(i));
            minitree_Hemi_Vtx_BDT_HMass.push_back(tree_Hemi_Vtx_BDT_HMass->at(i));
            minitree_Hemi_Vtx_BDT_ntrk10.push_back(tree_Hemi_Vtx_BDT_ntrk10->at(i));
            minitree_Hemi_Vtx_BDT_ntrk20.push_back(tree_Hemi_Vtx_BDT_ntrk20->at(i));
            minitree_Hemi_Vtx_BDT_MeanDCA.push_back(tree_Hemi_Vtx_BDT_MeanDCA->at(i));
            minitree_Hemi_Vtx_MVAval_Loose.push_back(tree_Hemi_Vtx_MVAval_Loose->at(i));
            minitree_Hemi_Vtx_MVAval_Tight.push_back(tree_Hemi_Vtx_MVAval_Tight->at(i));

            if ( tree_Filter && tree_njetNOmu > 0 ) {
              hData_Hemi->Fill( tree_Hemi->at(i) );
              if ( tree_Hemi_Vtx_step->at(i) >= 1 ) hData_Hemi_stepGE1->Fill( tree_Hemi->at(i) );
              if ( tree_Hemi_Vtx_step->at(i) >= 1 && tree_Hemi_Vtx_step->at(i) <= 2 ) hData_Hemi_step12->Fill( tree_Hemi->at(i) );
            }
         }

      if (debug)  {std::cout<<"SecVtxLLP"<<std::endl;}
      for (unsigned int i = 0; i < tree_Hemi_SecLLP->size(); i++)
      {
         minitree_Hemi_SecLLP.push_back(tree_Hemi_SecLLP->at(i));
         minitree_Hemi_LLP_SecVtx_dz.push_back(tree_Hemi_LLP_SecVtx_dz->at(i));
         minitree_Hemi_LLP_SecVtx_dr.push_back(tree_Hemi_LLP_SecVtx_dr->at(i));
         minitree_Hemi_SecLLP_ping.push_back(tree_Hemi_SecLLP_ping->at(i));
      }

      if (debug)  {std::cout<<"SecVtxLLP"<<std::endl;}
      for (unsigned int i = 0; i < tree_Hemi_SecVtx->size(); i++) 
      {
         minitree_Hemi_SecVtx.push_back(tree_Hemi_SecVtx->at(i));
         minitree_Hemi_SecVtx_step.push_back(tree_Hemi_SecVtx_step->at(i));
         minitree_Hemi_SecVtx_x.push_back(tree_Hemi_SecVtx_x->at(i));
         minitree_Hemi_SecVtx_y.push_back(tree_Hemi_SecVtx_y->at(i));
         minitree_Hemi_SecVtx_z.push_back(tree_Hemi_SecVtx_z->at(i));
         minitree_Hemi_SecVtx_r.push_back(tree_Hemi_SecVtx_r->at(i));
         minitree_Hemi_SecVtx_dR.push_back(tree_Hemi_SecVtx_dR->at(i));
         minitree_Hemi_SecVtx_nTrks.push_back(tree_Hemi_SecVtx_nTrks->at(i));
         minitree_Hemi_SecVtx_NChi2.push_back(tree_Hemi_SecVtx_NChi2->at(i));
         minitree_Hemi_SecVtx_dist.push_back(tree_Hemi_SecVtx_dist->at(i));
         minitree_Hemi_SecVtx_track_MeanDCA_d.push_back(tree_Hemi_SecVtx_track_MeanDCA_d->at(i));
         minitree_Hemi_SecVtx_SumtrackWeight.push_back(tree_Hemi_SecVtx_SumtrackWeight->at(i));
         minitree_Hemi_SecVtx_Mass.push_back(tree_Hemi_SecVtx_Mass->at(i));
      }

      //////////// FILLING //////////

      smalltree->Fill();

      //////////// CLEARING //////////

      minirunNumber.clear();
      minieventNumber.clear();
      minilumiBlock.clear();
      minitree_LHE_Weights.clear();
      minitree_MCEvt_weight.clear();
      minitree_only_gen_wt.clear();
      minitree_event_weight.clear();
      minitree_genTop_Weight.clear();
      minitree_gen_top_pt.clear();
      minitree_gen_top_rw_pt.clear();
      miniPUweight.clear();
      miniPUweight_Up.clear();
      miniPUweight_Down.clear();
      miniPrefweight.clear();
      miniPU_events.clear();
      miniAllPU_events_weight.clear();
      minitree_Filter.clear();
      minitree_FilterSameSign.clear();
      minitree_trigger_doublelepton.clear();
      minitree_trigger_singlelepton.clear();

      minitree_GenPVx.clear();
      minitree_GenPVy.clear();
      minitree_GenPVz.clear();
      minitree_smu_mass.clear();
      minitree_neu_mass.clear();
      minitree_neu_ctau.clear();

      minitree_Good_PV.clear();
      minitree_nPV.clear();
      minitree_PV_x.clear();
      minitree_PV_y.clear();
      minitree_PV_z.clear();
      minitree_PV_ez.clear();
      minitree_PV_NChi2.clear();
      minitree_PV_ndf.clear();

      minitree_PFMet_et.clear();
      minitree_PFMet_phi.clear();

      minitree_HT.clear();

      minitree_TRACK_SIZE.clear();
      minitree_nTracks.clear();
      minitree_nLostTracks.clear();

      minitree_muon_GenRecoTriggerMatched.clear();
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
      minitree_muon_gen.clear();
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

      minitree_Evts_MVAval.clear();
      minitree_Evts_MVAvalDY.clear();
      minitree_Evts_MVAvalTT.clear();
    
      minitree_nLLP.clear();
      minitree_LLP.clear();
      minitree_LLP_pt.clear();
      minitree_LLP_eta.clear();
      minitree_LLP_phi.clear();
      minitree_LLP_x.clear();
      minitree_LLP_y.clear();
      minitree_LLP_z.clear();
      minitree_LLP_r.clear();
      minitree_LLP_dist.clear();
      minitree_LLP_nTrks.clear();
      minitree_LLP12_dR.clear();
      minitree_LLP12_deta.clear();
      minitree_LLP12_dphi.clear();
      minitree_LLP_Mass.clear();

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

      minitree_Hemi_LLP.clear();
      minitree_Hemi_LLP_pt.clear();
      minitree_Hemi_LLP_eta.clear();
      minitree_Hemi_LLP_phi.clear();
      minitree_Hemi_LLP_dist.clear();
      minitree_Hemi_LLP_x.clear();
      minitree_Hemi_LLP_y.clear();
      minitree_Hemi_LLP_z.clear();
      minitree_Hemi_LLP_dR.clear();
      minitree_Hemi_LLP_mother.clear();
      minitree_Hemi_LLP_Vtx_dx.clear();
      minitree_Hemi_LLP_Vtx_dy.clear();
      minitree_Hemi_LLP_Vtx_dz.clear();
      minitree_Hemi_LLP_Vtx_dr.clear();
      minitree_Hemi_LLP_muOK_dR.clear();
      minitree_Hemi_LLP_muOK_pt.clear();
      minitree_Hemi_LLP_muOK_mass.clear();
      minitree_Hemi_LLP_muNO_dR.clear();
      minitree_Hemi_LLP_muNO_pt.clear();
      minitree_Hemi_LLP_muNO_mass.clear();
      minitree_Hemi_LLP_dR12.clear();
      minitree_Hemi_LLP_ping.clear();
      minitree_event_LLP_ping.clear();
    
      minitree_Hemi_Vtx_step.clear();
      minitree_Hemi_Vtx_isTight.clear();
      minitree_Hemi_Vtx_NChi2.clear();
      minitree_Hemi_Vtx_nTrks.clear();
      minitree_Hemi_Vtx_nTrks_sig.clear();
      minitree_Hemi_Vtx_nTrks_bad.clear();
      minitree_Hemi_Vtx_x.clear();
      minitree_Hemi_Vtx_y.clear();
      minitree_Hemi_Vtx_z.clear();
      minitree_Hemi_Vtx_r.clear();
      minitree_Hemi_Vtx_dR.clear();
      minitree_Hemi_Vtx_SumtrackWeight.clear();
      minitree_Hemi_Vtx_track_MeanDCA_d.clear();
      minitree_Hemi_Vtx_Mass.clear();
      minitree_Hemi_Vtx_dist.clear();

      minitree_event_nVtx.clear();
      minitree_event_Vtx_Vtx_dr.clear();
      minitree_event_Vtx_Vtx_dz.clear();
      minitree_event_Vtx_Vtx_dd.clear();
      minitree_event_Vtx_Vtx_reldd.clear();
      minitree_event_Vtx_Vtx_dR.clear();
      minitree_event_Vtx_Vtx_step.clear();

      minitree_Hemi_SecLLP.clear();
      minitree_Hemi_LLP_SecVtx_dz.clear();
      minitree_Hemi_LLP_SecVtx_dr.clear();
      minitree_Hemi_SecLLP_ping.clear();
      minitree_event_SecLLP_ping.clear();

      minitree_Hemi_SecVtx.clear();
      minitree_Hemi_SecVtx_step.clear();
      minitree_Hemi_SecVtx_x.clear();
      minitree_Hemi_SecVtx_y.clear();
      minitree_Hemi_SecVtx_z.clear();
      minitree_Hemi_SecVtx_r.clear();
      minitree_Hemi_SecVtx_dR.clear();
      minitree_Hemi_SecVtx_nTrks.clear();
      minitree_Hemi_SecVtx_NChi2.clear();
      minitree_Hemi_SecVtx_dist.clear();
      minitree_Hemi_SecVtx_track_MeanDCA_d.clear();
      minitree_Hemi_SecVtx_SumtrackWeight.clear();
      minitree_Hemi_SecVtx_Mass.clear();

      minitree_event_MergedVtx_Vtx_dr.clear();
      minitree_event_MergedVtx_Vtx_dz.clear();
      minitree_event_MergedVtx_Vtx_dd.clear();
      minitree_event_MergedVtx_Vtx_reldd.clear();
      minitree_event_MergedVtx_Vtx_dR.clear();
      minitree_event_MergedVtx_Vtx_step.clear();
        
      minitree_Hemi_Vtx_BDT_nTrks.clear();
      minitree_Hemi_Vtx_BDT_NChi2.clear();
      minitree_Hemi_Vtx_BDT_step.clear();
      minitree_Hemi_Vtx_BDT_STW.clear();
      minitree_Hemi_Vtx_BDT_Mass.clear();
      minitree_Hemi_Vtx_BDT_HMass.clear();
      minitree_Hemi_Vtx_BDT_ntrk10.clear();
      minitree_Hemi_Vtx_BDT_ntrk20.clear();
      minitree_Hemi_Vtx_BDT_MeanDCA.clear();
      minitree_Hemi_Vtx_MVAval_Loose.clear();
      minitree_Hemi_Vtx_MVAval_Tight.clear();

      miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v.clear();
      miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v.clear();
      miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v.clear();
      miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.clear();
      miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v.clear();
      miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.clear();
      miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v.clear();
      miniHLT_Ele27_WPTight_Gsf_v.clear();
      miniHLT_Ele32_WPTight_Gsf_v.clear();
      miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.clear();
      miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v.clear();
      miniHLT_PFMET120_PFMHT120_IDTight_v.clear();
      miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v.clear();
      miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v.clear();
      miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v.clear();
      miniHLT_PFMET250_HBHECleaned_v.clear();
      miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v.clear();
      miniHLT_IsoMu24_v.clear();
      miniHLT_IsoMu27_v.clear();

      miniHLT_IsoTkMu24_v.clear(); 
      miniHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v.clear();   
 
      miniHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v.clear();    
      miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v.clear();    
      miniHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v.clear();  
      miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v.clear();


     
   } // end loop on events 

   std::cout << "number of events  "<< allevents << std::endl;
   hData_Good_PV -> Draw(); 
   hData_Filter -> Draw(); 
   hData_FilterSameSign -> Draw(); 
   hData_njetNOmu_Filter -> Draw(); 
   hData_njetNOmu_FilterSameSign -> Draw(); 
   hData_njetNOmu -> Draw(); 

   hData_nPV -> Draw();
   hTk_MVA -> Draw();
   hData_Hemi -> Draw();
   hData_Hemi_stepGE1 -> Draw();
   hData_Hemi_step12 -> Draw();

   myFile->Write();
   delete myFile;

   // HistogramManager h ;
   // h.WriteAllHistogramsInFile((Production+"/Mini"+sample+".root").Data(),"recreate");
}

