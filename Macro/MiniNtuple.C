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
   TTree *tree = new TTree("T", "summary information");
  
   std::vector<int>              minirunNumber;
   std::vector<int>              minieventNumber;
   std::vector<int>              minilumiBlock;
   std::vector<float>            minitree_LHE_Weights;
   std::vector< Float_t >        minitree_MCEvt_weight;
   std::vector< Double_t >       minitree_only_gen_wt;
   std::vector< Double_t >       minitree_genTop_Weight;
   std::vector< Double_t >       miniPUweight;
   std::vector< Double_t >       miniPrefweight;
   std::vector<int>              miniPU_events;
   std::vector<int>              miniAllPU_events_weight;
    
   vector<bool>                  minitree_Good_PV;
   std::vector<int>              minitree_muon_GenRecoTriggerMatched;

   std::vector< float >         minitree_GenPVx;
   std::vector< float >         minitree_GenPVy;
   std::vector< float >         minitree_GenPVz;
   std::vector< int >           minitree_smu_mass; // used for signal only
   std::vector< int >           minitree_neu_mass; // used for signal only
   std::vector< int >           minitree_neu_ctau; // used for signal only
   
   vector<int>           minitree_nPV;
   vector<float>         minitree_PV_x;
   vector<float>         minitree_PV_y;
   vector<float>         minitree_PV_z;
   vector<float>         minitree_PV_ez;
   vector<float>         minitree_PV_NChi2;
   vector<int>           minitree_PV_ndf;

   vector<float>         minitree_PFMet_et; // in BDT_EVT
   vector<float>         minitree_PFMet_phi;

   vector<int>           minitree_njet;
   vector<int>           minitree_njetNOmu; // in BDT_EVT

   vector<float>         minitree_HT; // in BDT_EVT

   vector<int>           minitree_TRACK_SIZE; // in BDT_EVT
   vector<int>           minitree_nTracks;
   vector<int>           minitree_nLostTracks;

   vector<bool>    minitree_muon_isPrompt;
   vector<float>   minitree_muon_pt;
   vector<float>   minitree_muon_SF;
   vector<float>   minitree_muon_eta;
   vector<float>   minitree_muon_phi;
   vector<float>   minitree_muon_dxy;
   vector<float>   minitree_muon_dz;
   vector<int>     minitree_muon_charge;
   vector<float>   minitree_muon_correction;
   vector<int>     minitree_muon_gen;

   vector<float>   minitree_lepton_leadingpt;
   vector<float>   minitree_lepton_leadingpt2;
   vector<float>   minitree_lepton_lepton_dR;
   vector<float>   minitree_lepton_lepton_dPhi;
   vector<float>   minitree_lepton_lepton_dEta;

   vector<float>   minitree_ll_pt;
   vector<float>   minitree_ll_eta;
   vector<float>   minitree_ll_phi;
   vector<float>   minitree_ll_px;
   vector<float>   minitree_ll_py;
   vector<float>   minitree_ll_pz;
   vector<float>   minitree_ll_energy;
   vector<float>   minitree_ll_mass;

   std::vector < Int_t >           minitree_all_nel; 
   std::vector < Int_t >            minitree_electron_nEle; 
   vector<bool>    minitree_electron_isPrompt;
   vector<float>   minitree_electron_pt;
   vector<float>   minitree_electron_eta;
   vector<float>   minitree_electron_phi;
   vector<int>     minitree_electron_charge;
   vector<float>   minitree_electron_dxy;
   vector<float>   minitree_electron_dz;
   vector<int>     minitree_electron_gen;

   std::vector<float> minitree_jet_pt;
   std::vector<float> minitree_jet_eta;
   std::vector<float> minitree_jet_phi;
   std::vector<float> minitree_jet_HadronFlavour;
   std::vector<float> minitree_jet_btag_DeepJet;
   std::vector<float> minitree_jet_E;

   std::vector<float>   minitree_jet_leadingpt; // in BDT_EVT
   std::vector<float>   minitree_jet_leadingpt2; // in BDT_EVT
   std::vector<float>   minitree_jet_leadingeta; // in BDT_EVT
   std::vector<float>   minitree_jet_leadingeta2; // in BDT_EVT
   std::vector<float>   minitree_jet_jet_dR; // in BDT_EVT
   std::vector<float>   minitree_jet_jet_dPhi; // in BDT_EVT
   std::vector<float>   minitree_jet_jet_dEta; // in BDT_EVT
   std::vector<float>   minitree_muon_jet_dRmin;
   std::vector<float>   minitree_muon_jet_dRmax;
   std::vector<float>   minitree_elemu_jet_dRmin;
   std::vector<float>   minitree_elemu_jet_dRmax;
   std::vector<float>   minitree_ele_jet_dRmin; // empty, usefull ???
   std::vector<float>   minitree_ele_jet_dRmax; // empty, usefull ???

    std::vector<Bool_t> minitree_Filter;
    std::vector<Bool_t> minitree_FilterSameSign;

    std::vector<Float_t> minitree_Evts_MVAval;
    std::vector<Float_t> minitree_Evts_MVAvalDY;
    std::vector<Float_t> minitree_Evts_MVAvalTT;
      //---------------  Gen event wt------------------------//

    std::vector<double> minitree_event_weight;

    //--------------------------------
    // muons infos -------
    //--------------------------------
    std::vector<int> minitree_all_nmu; // count all muons
    std::vector<int> minitree_nmu;     // count prompt muons
    std::vector< float > minitree_LT;
    std::vector< float > minitree_Mmumu;
    std::vector< float > minitree_MmumuSameSign;

    std::vector<int>     minitree_nLLP;
    std::vector< int >   minitree_LLP;
    std::vector< float > minitree_LLP_pt;
    std::vector< float > minitree_LLP_eta;
    std::vector< float > minitree_LLP_phi;
    std::vector< float > minitree_LLP_x;
    std::vector< float > minitree_LLP_y;
    std::vector< float > minitree_LLP_z;
    std::vector< float > minitree_LLP_r;
    std::vector< float > minitree_LLP_dist;
    std::vector< int >   minitree_LLP_nTrks;
    std::vector< float > minitree_LLP12_dR;
    std::vector< float > minitree_LLP12_deta;
    std::vector< float > minitree_LLP12_dphi;
    std::vector< float > minitree_LLP_Mass;

    //-----------------------
    // - Vertices information
    //-----------------------

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
    std::vector< float > minitree_Hemi_dR12;

    std::vector< int >   minitree_Hemi_LLP;
    std::vector< float > minitree_Hemi_LLP_pt;
    std::vector< float > minitree_Hemi_LLP_eta;
    std::vector< float > minitree_Hemi_LLP_phi;
    std::vector< float > minitree_Hemi_LLP_dist;
    std::vector< float > minitree_Hemi_LLP_x;
    std::vector< float > minitree_Hemi_LLP_y;
    std::vector< float > minitree_Hemi_LLP_z;
    std::vector< float > minitree_Hemi_LLP_dR;
    std::vector< int >   minitree_Hemi_LLP_mother;
    std::vector< float > minitree_Hemi_LLP_Vtx_dx;
    std::vector< float > minitree_Hemi_LLP_Vtx_dy;
    std::vector< float > minitree_Hemi_LLP_Vtx_dz;
    std::vector< float > minitree_Hemi_LLP_Vtx_dr;
    std::vector< float > minitree_Hemi_LLP_muOK_dR;
    std::vector< float > minitree_Hemi_LLP_muOK_pt;
    std::vector< float > minitree_Hemi_LLP_muOK_mass;
    std::vector< float > minitree_Hemi_LLP_muNO_dR;
    std::vector< float > minitree_Hemi_LLP_muNO_pt;
    std::vector< float > minitree_Hemi_LLP_muNO_mass;
    std::vector< float > minitree_Hemi_LLP_dR12;
    std::vector< bool >  minitree_Hemi_LLP_ping; 
    std::vector< int >   minitree_event_LLP_ping;
    
    std::vector< int >   minitree_Hemi_Vtx_step;
    std::vector< bool >  minitree_Hemi_Vtx_isTight;
    std::vector< float > minitree_Hemi_Vtx_NChi2;
    std::vector< int >   minitree_Hemi_Vtx_nTrks;
    std::vector< int >   minitree_Hemi_Vtx_nTrks_sig;
    std::vector< int >   minitree_Hemi_Vtx_nTrks_bad;
    std::vector< float > minitree_Hemi_Vtx_x;
    std::vector< float > minitree_Hemi_Vtx_y;
    std::vector< float > minitree_Hemi_Vtx_z;
    std::vector< float > minitree_Hemi_Vtx_r;
    std::vector< float > minitree_Hemi_Vtx_dR;
    std::vector< float > minitree_Hemi_Vtx_SumtrackWeight;//Vertx selection variable for the BDT
    std::vector< float > minitree_Hemi_Vtx_Mass;
    std::vector< float > minitree_Hemi_Vtx_track_MeanDCA_d;//Veertex selection BDT
    std::vector< float > minitree_Hemi_Vtx_dist;
    std::vector< int >   minitree_event_nVtx;
    std::vector< float > minitree_event_Vtx_Vtx_dr;
    std::vector< float > minitree_event_Vtx_Vtx_dz;
    std::vector< float > minitree_event_Vtx_Vtx_dd;
    std::vector< float > minitree_event_Vtx_Vtx_reldd;
    std::vector< float > minitree_event_Vtx_Vtx_dR;
    std::vector< int >   minitree_event_Vtx_Vtx_step;

    std::vector< float > minitree_Hemi_SecLLP;
    std::vector< float > minitree_Hemi_LLP_SecVtx_dz;
    std::vector< float > minitree_Hemi_LLP_SecVtx_dr;
    std::vector< bool >  minitree_Hemi_SecLLP_ping;
    std::vector< int >   minitree_event_SecLLP_ping;

    std::vector< int >   minitree_Hemi_SecVtx;      // Hemi (1 or 2) if merging
    std::vector< int >   minitree_Hemi_SecVtx_step; // vertex step for this Hemi if merging
    std::vector< float > minitree_Hemi_SecVtx_x;
    std::vector< float > minitree_Hemi_SecVtx_y;
    std::vector< float > minitree_Hemi_SecVtx_z;
    std::vector< float > minitree_Hemi_SecVtx_r;
    std::vector< float > minitree_Hemi_SecVtx_dR;
    std::vector< float > minitree_Hemi_SecVtx_nTrks;
    std::vector< float > minitree_Hemi_SecVtx_NChi2;
    std::vector< float > minitree_Hemi_SecVtx_dist;
    std::vector< float > minitree_Hemi_SecVtx_track_MeanDCA_d;
    std::vector< float > minitree_Hemi_SecVtx_SumtrackWeight;
    std::vector< float > minitree_Hemi_SecVtx_Mass;
    
    std::vector< float > minitree_event_MergedVtx_Vtx_dr;
    std::vector< float > minitree_event_MergedVtx_Vtx_dz;
    std::vector< float > minitree_event_MergedVtx_Vtx_dd;
    std::vector< float > minitree_event_MergedVtx_Vtx_reldd;
    std::vector< float > minitree_event_MergedVtx_Vtx_dR;
    std::vector< int >   minitree_event_MergedVtx_Vtx_step;
    
    std::vector< float > minitree_Hemi_Vtx_BDT_nTrks;
    std::vector< float > minitree_Hemi_Vtx_BDT_NChi2;
    std::vector< float > minitree_Hemi_Vtx_BDT_step;
    std::vector< float > minitree_Hemi_Vtx_BDT_STW;
    std::vector< float > minitree_Hemi_Vtx_BDT_Mass;
    std::vector< float > minitree_Hemi_Vtx_BDT_HMass;
    std::vector< float > minitree_Hemi_Vtx_BDT_ntrk10;
    std::vector< float > minitree_Hemi_Vtx_BDT_ntrk20;
    std::vector< float > minitree_Hemi_Vtx_BDT_MeanDCA;
    std::vector< float > minitree_Hemi_Vtx_MVAval_Loose;
    std::vector< float > minitree_Hemi_Vtx_MVAval_Tight;//TIght WP
    
   std::vector< Bool_t >  miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;
   std::vector< Bool_t >  miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v;
   std::vector< Bool_t >  miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v;
   std::vector< Bool_t >  miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   std::vector< Bool_t >  miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;
   std::vector< Bool_t >  miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   std::vector< Bool_t >  miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;
   std::vector< Bool_t >  miniHLT_Ele27_WPTight_Gsf_v;
   std::vector< Bool_t >  miniHLT_Ele32_WPTight_Gsf_v;
   std::vector< Bool_t >  miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
   std::vector< Bool_t >  miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;
   std::vector< Bool_t >  miniHLT_PFMET120_PFMHT120_IDTight_v;
   std::vector< Bool_t >  miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v;
   std::vector< Bool_t >  miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;
   std::vector< Bool_t >  miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;
   std::vector< Bool_t >  miniHLT_PFMET250_HBHECleaned_v;
   std::vector< Bool_t >  miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;
   std::vector< Bool_t >  miniHLT_IsoMu24_v;
   std::vector< Bool_t >  miniHLT_IsoMu27_v;
   
  
   tree->Branch("minirunNumber",&minirunNumber);
   tree->Branch("minieventNumber",&minieventNumber);
   tree->Branch("minilumiBlock",&minilumiBlock);
   tree->Branch("minitree_LHE_Weights",&minitree_LHE_Weights);
   tree->Branch("minitree_MCEvt_weight",&minitree_MCEvt_weight);
   tree->Branch("minitree_only_gen_wt",&minitree_only_gen_wt);
   tree->Branch("minitree_genTop_Weight",&minitree_genTop_Weight);
   tree->Branch("miniPUweight",&miniPUweight);
   tree->Branch("miniPrefweight",&miniPrefweight);
   tree->Branch("miniPU_events",&miniPU_events);
   tree->Branch("miniAllPU_events_weight",&miniAllPU_events_weight);

   tree->Branch("minitree_Good_PV",&minitree_Good_PV);
   tree->Branch("minitree_muon_GenRecoTriggerMatched",&minitree_muon_GenRecoTriggerMatched);

   tree->Branch("minitree_GenPVx",&minitree_GenPVx);
   tree->Branch("minitree_GenPVy",&minitree_GenPVy);
   tree->Branch("minitree_GenPVz",&minitree_GenPVz);
   tree->Branch("minitree_smu_mass",&minitree_smu_mass); // used for signal only
   tree->Branch("minitree_neu_mass",&minitree_neu_mass); // used for signal only
   tree->Branch("minitree_neu_ctau",&minitree_neu_ctau); // used for signal only

   tree->Branch("minitree_nPV",&minitree_nPV);
   tree->Branch("minitree_PV_x",&minitree_PV_x);
   tree->Branch("minitree_PV_y",&minitree_PV_y);
   tree->Branch("minitree_PV_z",&minitree_PV_z);
   tree->Branch("minitree_PV_ez",&minitree_PV_ez);
   tree->Branch("minitree_PV_NChi2",&minitree_PV_NChi2);
   tree->Branch("minitree_PV_ndf",&minitree_PV_ndf);

   tree->Branch("minitree_PFMet_et",&minitree_PFMet_et); // in BDT_EVT
   tree->Branch("minitree_PFMet_phi",&minitree_PFMet_phi);

   tree->Branch("minitree_njet",&minitree_njet);
   tree->Branch("minitree_njetNOmu",&minitree_njetNOmu); // in BDT_EVT

   tree->Branch("minitree_HT",&minitree_HT); // in BDT_EVT

   tree->Branch("minitree_TRACK_SIZE",&minitree_TRACK_SIZE); // in BDT_EVT
   tree->Branch("minitree_nTracks",&minitree_nTracks);
   tree->Branch("minitree_nLostTracks",&minitree_nLostTracks);

   tree->Branch("minitree_all_nmu",&minitree_all_nmu); 
   tree->Branch("minitree_nmu",&minitree_nmu);    
   tree->Branch("minitree_LT",&minitree_LT);
   tree->Branch("minitree_Mmumu",&minitree_Mmumu);
   tree->Branch("minitree_MmumuSameSign",&minitree_MmumuSameSign);

   tree->Branch("minitree_muon_isPrompt",&minitree_muon_isPrompt);
   tree->Branch("minitree_muon_pt",&minitree_muon_pt);
   tree->Branch("minitree_muon_SF",&minitree_muon_SF);
   tree->Branch("minitree_muon_eta",&minitree_muon_eta);
   tree->Branch("minitree_muon_phi",&minitree_muon_phi);
   tree->Branch("minitree_muon_dxy",&minitree_muon_dxy);
   tree->Branch("minitree_muon_dz",&minitree_muon_dz);
   tree->Branch("minitree_muon_charge",&minitree_muon_charge);
   tree->Branch("minitree_muon_correction",&minitree_muon_correction);
   tree->Branch("minitree_muon_gen",&minitree_muon_gen);

   tree->Branch("minitree_lepton_leadingpt",&minitree_lepton_leadingpt);
   tree->Branch("minitree_lepton_leadingpt2",&minitree_lepton_leadingpt2);
   tree->Branch("minitree_lepton_lepton_dR",&minitree_lepton_lepton_dR);
   tree->Branch("minitree_lepton_lepton_dPhi",&minitree_lepton_lepton_dPhi);
   tree->Branch("minitree_lepton_lepton_dEta",&minitree_lepton_lepton_dEta);


   tree->Branch("minitree_ll_pt",&minitree_ll_pt);
   tree->Branch("minitree_ll_eta",&minitree_ll_eta);
   tree->Branch("minitree_ll_phi",&minitree_ll_phi);
   tree->Branch("minitree_ll_px",&minitree_ll_px);
   tree->Branch("minitree_ll_py",&minitree_ll_py);
   tree->Branch("minitree_ll_pz",&minitree_ll_pz);
   tree->Branch("minitree_ll_energy",&minitree_ll_energy);
   tree->Branch("minitree_ll_mass",&minitree_ll_mass);

   tree->Branch("minitree_all_nel",&minitree_all_nel); 
   tree->Branch("minitree_electron_nEle",&minitree_electron_nEle); 
   tree->Branch("minitree_electron_isPrompt",&minitree_electron_isPrompt);
   tree->Branch("minitree_electron_pt",&minitree_electron_pt);
   tree->Branch("minitree_electron_eta",&minitree_electron_eta);
   tree->Branch("minitree_electron_phi",&minitree_electron_phi);
   tree->Branch("minitree_electron_charge",&minitree_electron_charge);
   tree->Branch("minitree_electron_dxy",&minitree_electron_dxy);
   tree->Branch("minitree_electron_dz",&minitree_electron_dz);
   tree->Branch("minitree_electron_gen",&minitree_electron_gen);

   tree->Branch("minitree_jet_pt",&minitree_jet_pt);
   tree->Branch("minitree_jet_eta",&minitree_jet_eta);
   tree->Branch("minitree_jet_phi",&minitree_jet_phi);
   tree->Branch("minitree_jet_HadronFlavour",&minitree_jet_HadronFlavour);
   tree->Branch("minitree_jet_btag_DeepJet",&minitree_jet_btag_DeepJet);
   tree->Branch("minitree_jet_E",&minitree_jet_E);

   tree->Branch("minitree_jet_leadingpt",&minitree_jet_leadingpt); // in BDT_EVT
   tree->Branch("minitree_jet_leadingpt2",&minitree_jet_leadingpt2); // in BDT_EVT
   tree->Branch("minitree_jet_leadingeta",&minitree_jet_leadingeta); // in BDT_EVT
   tree->Branch("minitree_jet_leadingeta2",&minitree_jet_leadingeta2); // in BDT_EVT
   tree->Branch("minitree_jet_jet_dR",&minitree_jet_jet_dR); // in BDT_EVT
   tree->Branch("minitree_jet_jet_dPhi",&minitree_jet_jet_dPhi); // in BDT_EVT
   tree->Branch("minitree_jet_jet_dEta",&minitree_jet_jet_dEta); // in BDT_EVT
   tree->Branch("minitree_muon_jet_dRmin",&minitree_muon_jet_dRmin);
   tree->Branch("minitree_muon_jet_dRmax",&minitree_muon_jet_dRmax);
   tree->Branch("minitree_elemu_jet_dRmin",&minitree_elemu_jet_dRmin);
   tree->Branch("minitree_elemu_jet_dRmax",&minitree_elemu_jet_dRmax);
   tree->Branch("minitree_ele_jet_dRmin",&minitree_ele_jet_dRmin); // empty, usefull ???
   tree->Branch("minitree_ele_jet_dRmax",&minitree_ele_jet_dRmax); // empty, usefull ???

    tree->Branch("minitree_event_weight",&minitree_event_weight);
    tree->Branch("minitree_Filter",        &minitree_Filter);
    tree->Branch("minitree_FilterSameSign",&minitree_FilterSameSign);
    tree->Branch("minitree_Evts_MVAval",   &minitree_Evts_MVAval);
    tree->Branch("minitree_Evts_MVAvalDY",   &minitree_Evts_MVAvalDY);
    tree->Branch("minitree_Evts_MVAvalTT",   &minitree_Evts_MVAvalTT);

    tree->Branch("minitree_nLLP",          &minitree_nLLP);
    tree->Branch("minitree_LLP",           &minitree_LLP);
    tree->Branch("minitree_LLP_pt" ,       &minitree_LLP_pt);
    tree->Branch("minitree_LLP_eta",       &minitree_LLP_eta);
    tree->Branch("minitree_LLP_phi",       &minitree_LLP_phi);
    tree->Branch("minitree_LLP_x",         &minitree_LLP_x);
    tree->Branch("minitree_LLP_y",         &minitree_LLP_y);
    tree->Branch("minitree_LLP_z",         &minitree_LLP_z);
    tree->Branch("minitree_LLP_r",         &minitree_LLP_r);
    tree->Branch("minitree_LLP_dist",      &minitree_LLP_dist);
    tree->Branch("minitree_LLP_nTrks",     &minitree_LLP_nTrks);
    tree->Branch("minitree_LLP12_dR",      &minitree_LLP12_dR);
    tree->Branch("minitree_LLP12_deta",    &minitree_LLP12_deta);
    tree->Branch("minitree_LLP12_dphi",    &minitree_LLP12_dphi);
    tree->Branch("minitree_LLP_Mass",      &minitree_LLP_Mass);

    tree->Branch("minitree_Hemi",       &minitree_Hemi);
    tree->Branch("minitree_Hemi_njet",  &minitree_Hemi_njet);
    tree->Branch("minitree_Hemi_njet_nomu",  &minitree_Hemi_njet_nomu);
    tree->Branch("minitree_Hemi_pt",    &minitree_Hemi_pt);
    tree->Branch("minitree_Hemi_eta",   &minitree_Hemi_eta);
    tree->Branch("minitree_Hemi_phi",   &minitree_Hemi_phi);
    tree->Branch("minitree_Hemi_nTrks", &minitree_Hemi_nTrks);
    tree->Branch("minitree_Hemi_nTrks_sig", &minitree_Hemi_nTrks_sig);
    tree->Branch("minitree_Hemi_nTrks_bad", &minitree_Hemi_nTrks_bad);
    tree->Branch("minitree_Hemi_mass",     &minitree_Hemi_mass);
    tree->Branch("minitree_HemiMu_mass",   &minitree_HemiMu_mass);
    tree->Branch("minitree_HemiMu_pt",     &minitree_HemiMu_pt);
    tree->Branch("minitree_HemiMu_dR",     &minitree_HemiMu_dR);
    tree->Branch("minitree_HemiMuOp_mass", &minitree_HemiMuOp_mass);
    tree->Branch("minitree_HemiMuOp_pt",   &minitree_HemiMuOp_pt);
    tree->Branch("minitree_HemiMuOp_dR",   &minitree_HemiMuOp_dR);

    tree->Branch("minitree_Hemi_dR12",      &minitree_Hemi_dR12);

    tree->Branch("minitree_Hemi_LLP",       &minitree_Hemi_LLP);
    tree->Branch("minitree_Hemi_LLP_pt",    &minitree_Hemi_LLP_pt);
    tree->Branch("minitree_Hemi_LLP_eta",   &minitree_Hemi_LLP_eta);
    tree->Branch("minitree_Hemi_LLP_phi",   &minitree_Hemi_LLP_phi);
    tree->Branch("minitree_Hemi_LLP_dist",  &minitree_Hemi_LLP_dist);
    tree->Branch("minitree_Hemi_LLP_x",     &minitree_Hemi_LLP_x);
    tree->Branch("minitree_Hemi_LLP_y",     &minitree_Hemi_LLP_y);
    tree->Branch("minitree_Hemi_LLP_z",     &minitree_Hemi_LLP_z);
    tree->Branch("minitree_Hemi_LLP_dR",    &minitree_Hemi_LLP_dR);
    tree->Branch("minitree_Hemi_LLP_mother",&minitree_Hemi_LLP_mother);
    tree->Branch("minitree_Hemi_LLP_Vtx_dx",     &minitree_Hemi_LLP_Vtx_dx);
    tree->Branch("minitree_Hemi_LLP_Vtx_dy",     &minitree_Hemi_LLP_Vtx_dy);
    tree->Branch("minitree_Hemi_LLP_Vtx_dz",     &minitree_Hemi_LLP_Vtx_dz);
    tree->Branch("minitree_Hemi_LLP_Vtx_dr",     &minitree_Hemi_LLP_Vtx_dr);
    tree->Branch("minitree_Hemi_LLP_muOK_dR",    &minitree_Hemi_LLP_muOK_dR);
    tree->Branch("minitree_Hemi_LLP_muOK_pt",    &minitree_Hemi_LLP_muOK_pt);
    tree->Branch("minitree_Hemi_LLP_muOK_mass",  &minitree_Hemi_LLP_muOK_mass);
    tree->Branch("minitree_Hemi_LLP_muNO_dR",    &minitree_Hemi_LLP_muNO_dR);
    tree->Branch("minitree_Hemi_LLP_muNO_pt",    &minitree_Hemi_LLP_muNO_pt);
    tree->Branch("minitree_Hemi_LLP_muNO_mass",  &minitree_Hemi_LLP_muNO_mass);
    tree->Branch("minitree_Hemi_LLP_dR12",  &minitree_Hemi_LLP_dR12);
    tree->Branch("minitree_Hemi_LLP_ping",  &minitree_Hemi_LLP_ping);
    tree->Branch("minitree_event_LLP_ping", &minitree_event_LLP_ping);
    
    tree->Branch("minitree_Hemi_Vtx_step",  &minitree_Hemi_Vtx_step);
    tree->Branch("minitree_Hemi_Vtx_isTight",&minitree_Hemi_Vtx_isTight);
    tree->Branch("minitree_Hemi_Vtx_NChi2", &minitree_Hemi_Vtx_NChi2);
    tree->Branch("minitree_Hemi_Vtx_nTrks", &minitree_Hemi_Vtx_nTrks);
    tree->Branch("minitree_Hemi_Vtx_nTrks_sig", &minitree_Hemi_Vtx_nTrks_sig);
    tree->Branch("minitree_Hemi_Vtx_nTrks_bad", &minitree_Hemi_Vtx_nTrks_bad);
    tree->Branch("minitree_Hemi_Vtx_x",     &minitree_Hemi_Vtx_x);
    tree->Branch("minitree_Hemi_Vtx_y",     &minitree_Hemi_Vtx_y);
    tree->Branch("minitree_Hemi_Vtx_z",     &minitree_Hemi_Vtx_z);
    tree->Branch("minitree_Hemi_Vtx_r",     &minitree_Hemi_Vtx_r);
    tree->Branch("minitree_Hemi_Vtx_dR",    &minitree_Hemi_Vtx_dR);
    tree->Branch("minitree_Hemi_Vtx_SumtrackWeight",&minitree_Hemi_Vtx_SumtrackWeight);
    tree->Branch("minitree_Hemi_Vtx_track_MeanDCA_d",&minitree_Hemi_Vtx_track_MeanDCA_d);
    tree->Branch("minitree_Hemi_Vtx_Mass", &minitree_Hemi_Vtx_Mass);
    tree->Branch("minitree_Hemi_Vtx_dist",  &minitree_Hemi_Vtx_dist);
    tree->Branch("minitree_event_nVtx",      &minitree_event_nVtx);
    tree->Branch("minitree_event_Vtx_Vtx_dr",&minitree_event_Vtx_Vtx_dr);
    tree->Branch("minitree_event_Vtx_Vtx_dz",&minitree_event_Vtx_Vtx_dz);
    tree->Branch("minitree_event_Vtx_Vtx_dd",&minitree_event_Vtx_Vtx_dd);
    tree->Branch("minitree_event_Vtx_Vtx_reldd",&minitree_event_Vtx_Vtx_reldd);
    tree->Branch("minitree_event_Vtx_Vtx_dR",&minitree_event_Vtx_Vtx_dR);
    tree->Branch("minitree_event_Vtx_Vtx_step",&minitree_event_Vtx_Vtx_step);

    tree->Branch("minitree_Hemi_SecLLP",&minitree_Hemi_SecLLP);
    tree->Branch("minitree_Hemi_LLP_SecVtx_dz",&minitree_Hemi_LLP_SecVtx_dz);
    tree->Branch("minitree_Hemi_LLP_SecVtx_dr",&minitree_Hemi_LLP_SecVtx_dr);
    tree->Branch("minitree_Hemi_SecLLP_ping",&minitree_Hemi_SecLLP_ping);
    tree->Branch("minitree_event_SecLLP_ping",&minitree_event_SecLLP_ping);

    tree->Branch("minitree_Hemi_SecVtx",     &minitree_Hemi_SecVtx);
    tree->Branch("minitree_Hemi_SecVtx_step",&minitree_Hemi_SecVtx_step);
    tree->Branch("minitree_Hemi_SecVtx_x",&minitree_Hemi_SecVtx_x);
    tree->Branch("minitree_Hemi_SecVtx_y",&minitree_Hemi_SecVtx_y);
    tree->Branch("minitree_Hemi_SecVtx_z",&minitree_Hemi_SecVtx_z);
    tree->Branch("minitree_Hemi_SecVtx_r",&minitree_Hemi_SecVtx_r);
    tree->Branch("minitree_Hemi_SecVtx_dR",&minitree_Hemi_SecVtx_dR);
    tree->Branch("minitree_Hemi_SecVtx_nTrks",&minitree_Hemi_SecVtx_nTrks);
    tree->Branch("minitree_Hemi_SecVtx_NChi2", &minitree_Hemi_SecVtx_NChi2);
    tree->Branch("minitree_Hemi_SecVtx_dist",&minitree_Hemi_SecVtx_dist);
    tree->Branch("minitree_Hemi_SecVtx_track_MeanDCA_d",&minitree_Hemi_SecVtx_track_MeanDCA_d);
    tree->Branch("minitree_Hemi_SecVtx_SumtrackWeight",&minitree_Hemi_SecVtx_SumtrackWeight);
    tree->Branch("minitree_Hemi_SecVtx_Mass",&minitree_Hemi_SecVtx_Mass);
    tree->Branch("minitree_event_MergedVtx_Vtx_dr",&minitree_event_MergedVtx_Vtx_dr);
    tree->Branch("minitree_event_MergedVtx_Vtx_dz",&minitree_event_MergedVtx_Vtx_dz);
    tree->Branch("minitree_event_MergedVtx_Vtx_dd",&minitree_event_MergedVtx_Vtx_dd);
    tree->Branch("minitree_event_MergedVtx_Vtx_reldd",&minitree_event_MergedVtx_Vtx_reldd);
    tree->Branch("minitree_event_MergedVtx_Vtx_dR",&minitree_event_MergedVtx_Vtx_dR);
    tree->Branch("minitree_event_MergedVtx_Vtx_step",&minitree_event_MergedVtx_Vtx_step);


    tree->Branch("minitree_Hemi_Vtx_BDT_nTrks",&minitree_Hemi_Vtx_BDT_nTrks);
    tree->Branch("minitree_Hemi_Vtx_BDT_NChi2",&minitree_Hemi_Vtx_BDT_NChi2);
    tree->Branch("minitree_Hemi_Vtx_BDT_step",&minitree_Hemi_Vtx_BDT_step);
    tree->Branch("minitree_Hemi_Vtx_BDT_STW",&minitree_Hemi_Vtx_BDT_STW);
    tree->Branch("minitree_Hemi_Vtx_BDT_Mass",&minitree_Hemi_Vtx_BDT_Mass);
    tree->Branch("minitree_Hemi_Vtx_BDT_HMass",&minitree_Hemi_Vtx_BDT_HMass);
    tree->Branch("minitree_Hemi_Vtx_BDT_ntrk10",&minitree_Hemi_Vtx_BDT_ntrk10);
    tree->Branch("minitree_Hemi_Vtx_BDT_ntrk20",&minitree_Hemi_Vtx_BDT_ntrk20);
    tree->Branch("minitree_Hemi_Vtx_BDT_MeanDCA",&minitree_Hemi_Vtx_BDT_MeanDCA);
    tree->Branch("minitree_Hemi_Vtx_MVAval_Loose", &minitree_Hemi_Vtx_MVAval_Loose);
    tree->Branch("minitree_Hemi_Vtx_MVAval_Tight",&minitree_Hemi_Vtx_MVAval_Tight);

   tree->Branch("miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",&miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
   tree->Branch("miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",&miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v);
   tree->Branch("miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",&miniHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v);
   tree->Branch("miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   tree->Branch("miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",&miniHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
   tree->Branch("miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",&miniHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   tree->Branch("miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&miniHLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
   tree->Branch("miniHLT_Ele27_WPTight_Gsf_v",&miniHLT_Ele27_WPTight_Gsf_v);
   tree->Branch("miniHLT_Ele32_WPTight_Gsf_v",&miniHLT_Ele32_WPTight_Gsf_v);
   tree->Branch("miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
   tree->Branch("miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",&miniHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);
   tree->Branch("miniHLT_PFMET120_PFMHT120_IDTight_v",&miniHLT_PFMET120_PFMHT120_IDTight_v);
   tree->Branch("miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v",&miniHLT_PFMET120_PFMHT120_IDTight_PFHT60_v);
   tree->Branch("miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",&miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v);
   tree->Branch("miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",&miniHLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
   tree->Branch("miniHLT_PFMET250_HBHECleaned_v",&miniHLT_PFMET250_HBHECleaned_v);
   tree->Branch("miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",&miniHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v);
   tree->Branch("miniHLT_IsoMu24_v",&miniHLT_IsoMu24_v);
   tree->Branch("miniHLT_IsoMu27_v",&miniHLT_IsoMu27_v);

   minitree_Filter.reserve(100000000000);// reserve memory for the vector :O(2^37) due to ttbar large sampel size
   minitree_FilterSameSign.reserve(100000000000);// reserve memory for the vector :O(2^37)
   // std::cout<<minitree_Filter.max_size()<<std::endl;
   // std::cout<<minitree_Filter.capacity()<<std::endl;
   std::cout<<"//-------------------------//"<<std::endl;
   std::cout<<"// "<<Production<<" : "<<sample<<"  //"<<std::endl;
   std::cout<<"//-------------------------//"<<std::endl;

   bool debug = false;
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   nentries = 100000;
   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if ( jentry%1000 == 0 ) std::cout << "events : " << jentry << std::endl;
      if (jentry==0.1*nentries) {std::cout<<"10/100 :"<<std::endl;}
      if (jentry==0.2*nentries) {std::cout<<"20/100 :"<<std::endl;}
      if (jentry==0.3*nentries) {std::cout<<"30/100 :"<<std::endl;}
      if (jentry==0.4*nentries) {std::cout<<"40/100 :"<<std::endl;}
      if (jentry==0.5*nentries) {std::cout<<"50/100 :"<<std::endl;}
      if (jentry==0.6*nentries) {std::cout<<"60/100 :"<<std::endl;}
      if (jentry==0.7*nentries) {std::cout<<"70/100 :"<<std::endl;}
      if (jentry==0.8*nentries) {std::cout<<"80/100 :"<<std::endl;}
      if (jentry==0.9*nentries) {std::cout<<"90/100 :"<<std::endl;}

      //----------------------//
      //        Event         //
      //----------------------//

      minitree_Filter.push_back(tree_Filter);
      minitree_FilterSameSign.push_back(tree_FilterSameSign);

      //------------------------------------------------------------------------------------//
      //------------------------------------------------------------------------------------//
      if ( !((tree_Filter || tree_FilterSameSign) && tree_njetNOmu>0) && !Signal)  continue;
      //------------------------------------------------------------------------------------//
      //------------------------------------------------------------------------------------//

      if (debug)  {std::cout<<"Weights  "<<std::endl;}
      minirunNumber.push_back(runNumber);
      minieventNumber.push_back(eventNumber);
      minilumiBlock.push_back(lumiBlock);
      for (unsigned int i = 0; i < tree_LHE_Weights->size(); i++)
         {
            minitree_LHE_Weights.push_back(tree_LHE_Weights->at(i));
         }
      
      minitree_MCEvt_weight.push_back(tree_MCEvt_weight);
      minitree_only_gen_wt.push_back(tree_only_gen_wt);
      minitree_genTop_Weight.push_back(tree_genTop_Weight);
      miniPUweight.push_back(PUweight);
      miniPrefweight.push_back(Prefweight);
      miniPU_events.push_back(PU_events);
      miniAllPU_events_weight.push_back(AllPU_events_weight);

      if (Signal)
         {
            if (debug)  {std::cout<<"Signal parameters "<<std::endl;}
            minitree_smu_mass.push_back(tree_smu_mass); // used for signal only
            minitree_neu_mass.push_back(tree_neu_mass); // used for signal only
            minitree_neu_ctau.push_back(tree_neu_ctau); // used for signal only
         }

      minitree_Good_PV.push_back(tree_Good_PV);
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

      if (debug)  {std::cout<<"dilepton invariant quantities  "<<std::endl;}
            if ((tree_Filter || tree_FilterSameSign) && tree_njetNOmu>0)
         {
            minitree_ll_pt.push_back(tree_ll_pt->at(0));
            minitree_ll_eta.push_back(tree_ll_eta->at(0));
            minitree_ll_phi.push_back(tree_ll_phi->at(0));
            minitree_ll_px.push_back(tree_ll_px->at(0));
            minitree_ll_py.push_back(tree_ll_py->at(0));
            minitree_ll_pz.push_back(tree_ll_pz->at(0));
            minitree_ll_energy.push_back(tree_ll_energy->at(0));
            minitree_ll_mass.push_back(tree_ll_mass->at(0));
         }
      minitree_Evts_MVAval.push_back(tree_Evts_MVAval);
      minitree_Evts_MVAvalDY.push_back(tree_Evts_MVAvalDY);
      minitree_Evts_MVAvalTT.push_back(tree_Evts_MVAvalTT);
      minitree_nLLP.push_back(tree_nLLP);
      minitree_event_weight.push_back(tree_event_weight);

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
         }
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
            minitree_lepton_lepton_dR.push_back(tree_lepton_lepton_dR->at(0));
            minitree_lepton_lepton_dPhi.push_back(tree_lepton_lepton_dPhi->at(0));
            minitree_lepton_lepton_dEta.push_back(tree_lepton_lepton_dEta->at(0));
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


      // depends on the decay channel :(
      // minitree_elemu_jet_dRmin.push_back(tree_elemu_jet_dRmin->at(0) );
      // minitree_elemu_jet_dRmax.push_back( tree_elemu_jet_dRmax->at(0));
      // minitree_ele_jet_dRmin.push_back( tree_ele_jet_dRmin->at(0)); // empty, usefull ???
      // minitree_ele_jet_dRmax.push_back( tree_ele_jet_dRmax->at(0)); // empty, usefull ???

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
               // tree_event_Vtx_Vtx_reldd.push_back( tree_event_Vtx_Vtx_reldd->at(i));
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
               // tree_event_MergedVtx_Vtx_reldd.push_back( tree_event_MergedVtx_Vtx_reldd->at(i));
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


   }
   tree->Fill();
   myFile->Write();
   delete myFile;
   // HistogramManager h ;
   // h.WriteAllHistogramsInFile((Production+"/Mini"+sample+".root").Data(),"recreate");

}
