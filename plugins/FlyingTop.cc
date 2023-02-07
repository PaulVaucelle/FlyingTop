// system include files
#include <memory>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <map>
#include <utility>
#include <TNtuple.h>
#include <bitset>

// user include files
#include "TH2F.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TVector3.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "boost/functional/hash.hpp"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
              //CMSSW13 interface//
// #include "../interface/PackedCandidate.h"
//------------------------------------
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"//!!!!!!!!!!!!!!!!!!!!
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//!!!!!!
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/VecArray.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/Math/interface/liblogintpack.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
//!!!!

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FlyingTop/FlyingTop/interface/Proto.h"
#include "FlyingTop/FlyingTop/interface/DeltaFunc.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
//
//---------------------------------Paul-----------------------------//
              //-----------Transient Track/Vtx--------//
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "FWCore/Framework/interface/ESTransientHandle.h"
              //-------------Propagators------------//
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GeomPropagators/interface/StraightLinePlaneCrossing.h"
              //-------------Surfaces---------------//
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
              //----------------?-----------------//
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
              //----------------BField--------------//
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"
              //----------------New interface----------------------//
#include "../interface/PropaHitPattern.h"
              //----------------Trigger---------------------------//
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include <TEfficiency.h>
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//------------------------------End of Paul------------------------//


//
// class declaration
//

// skeleton from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#4_7_MiniAOD_Analysis_Documentati
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class FlyingTopAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
// class FlyingTopAnalyzer : public edm::EDAnalyzer {
  public:
    explicit FlyingTopAnalyzer(const edm::ParameterSet&);
    ~FlyingTopAnalyzer() {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------

    void clearVariables();

    std::string weightFile_;
    std::string weightFileLost_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<pat::JetCollection> jetToken_;
    edm::EDGetTokenT<edm::View<reco::GenJet> > genJetToken_;
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pcToken_;
    edm::EDGetTokenT<pat::PackedCandidateCollection>  lostpcToken_; //LOST

    
    std::string parametersDefinerName_;
    
  ///////////////
  // Ntuple info

    TTree *smalltree;
    
    edm::Service<TFileService> fs;
    
//     std::string ttrhbuilder_;
    
    edm::ESHandle<MagneticField> bField;
    
    edm::ParameterSet kvfPSet;
    //trig
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
    //trig
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> K0Token_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> LambdaToken_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> SVToken_;
    edm::EDGetTokenT<reco::ConversionCollection> PhotonToken_;
    // edm::EDGetTokenT<pat::PackedTriggerPrescales> PrescaleToken_;
    int runNumber, eventNumber, lumiBlock;
    int  tree_NbrOfZCand;
    bool tree_Filter;
    int  tree_nTracks, tree_nLostTracks; 
    int  nBC = 0, tree_nFromC = 0, tree_nFromB = 0; 
    int nEvent;
    
    float LLP1_pt, LLP1_eta, LLP1_phi, LLP2_pt, LLP2_eta, LLP2_phi;
    float LLP1_x, LLP1_y, LLP1_z, LLP2_x, LLP2_y, LLP2_z;
    float LLP1_dist, LLP2_dist;
    int   LLP1_nTrks = 0, LLP2_nTrks = 0;

//     //-----------------------
//     // trigger variable
// //     std::vector<string > tree_trigger_names;
// //     std::vector<bool >   tree_trigger_bits;
//     std::vector<int>    tree_trigger_size;
//     std::vector<int>    tree_passesTrigger;
//     std::vector<string> tree_passesTriggerName;
//     std::vector<string> tree_Trigger_Muon;//+ dilepton channel emu
//     std::vector<string> tree_Trigger_Ele;
//     std::vector<string> tree_Trigger_DoubleMu;
//     std::vector<string> tree_Trigger_DoubleEle;
//     std::vector<string> tree_Trigger_Dimuon0;
//     std::vector<string> tree_Trigger_PFMET;
//     std::vector<string> tree_Trigger_HT;
//     std::vector<string> tree_Trigger_AK4;
//     std::vector<string> tree_Trigger_PFJet;
//     std::vector<string> tree_Trigger_DoublePFJets;
//     std::vector<string> tree_Trigger_DiPFJet;
//     std::vector<string> tree_Trigger_QuadPFJet;
//     std::vector<string> tree_Trigger_BTagMu;
    
//$$
//The BDT variables are declared here to reduce computation time
    float pt, eta, NChi, nhits, ntrk10, drSig, isinjet;
    TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );

    bool NewCovMat = true;//Allow for Covariance Matrix correction due to the MiniAOD dataformat apporixmation
    bool IterAVF = true; // Activate IAVF step of the vertexing 
    bool ActivateTrigger = true;
    bool TrackMatchingToV0 = true; //Matching between tracks from V0Candidates collection and pfcandidate
    int index[1000];
    double MVAval[1000];
//$$

    //--------------------------------
    // primary vertex infos -------
    //--------------------------------
  
    int   tree_nPV;
    std::vector<float> tree_PV_x;
    std::vector<float> tree_PV_y;
    std::vector<float> tree_PV_z;
    std::vector<float> tree_PV_ez;
    std::vector<float> tree_PV_NChi2;
    std::vector<float> tree_PV_ndf;

    int tree_nSV;
    std::vector<float> tree_SV_x;
    std::vector<float> tree_SV_y;
    std::vector<float> tree_SV_z;
    std::vector<float> tree_SV_r;
    std::vector<float> tree_SV_NChi2;
    std::vector<float> tree_SV_ndf;
    std::vector<int>   tree_SV_nDaughters;
    std::vector<float> tree_SV_mass;
    std::vector<int>   tree_SV_pdgid;
    std::vector<float> tree_SV_pt;
    std::vector<float> tree_SV_eta;
    std::vector<float> tree_SV_phi;
    std::vector<bool>  tree_SV_IsConvertedPhoton;
//$$$$    std::vector<int>   tree_SV_daughter_pdgid;

    int tree_nK0;
    std::vector<float>     tree_K0_x;
    std::vector<float>     tree_K0_y;
    std::vector<float>     tree_K0_z;
    std::vector<float>     tree_K0_r;
    std::vector<float>     tree_K0_NChi2;
    std::vector<float>     tree_K0_ndf;
    std::vector<float>     tree_K0_mass;
    std::vector<float>     tree_K0_pt;
    std::vector<float>     tree_K0_eta;
    std::vector<float>     tree_K0_phi;
    std::vector<unsigned int> tree_K0_nDaughters;
    std::vector<float>     tree_K0_daughter_pt;
    std::vector<float>     tree_K0_daughter_eta;
    std::vector<float>     tree_K0_daughter_phi;

    int tree_nLambda;
    std::vector<float>     tree_L0_x;
    std::vector<float>     tree_L0_y;
    std::vector<float>     tree_L0_z;
    std::vector<float>     tree_L0_r;
    std::vector<float>     tree_L0_NChi2;
    std::vector<float>     tree_L0_ndf;
    std::vector<unsigned int> tree_L0_nDaughters;
    std::vector<float>     tree_L0_mass;
    std::vector<float>     tree_L0_pt;
    std::vector<float>     tree_L0_eta;
    std::vector<float>     tree_L0_phi;
    std::vector<float>     tree_L0_daughter_pt;
    std::vector<float>     tree_L0_daughter_eta;
    std::vector<float>     tree_L0_daughter_phi;

    int tree_nYConv;
    std::vector<float>     tree_Yc_x; 
    std::vector<float>     tree_Yc_y;
    std::vector<float>     tree_Yc_z;
    std::vector<float>     tree_Yc_r;
    std::vector<float>     tree_Yc_Vtx_NChi2;
    std::vector<float>     tree_Yc_Vtx_ndf;
    std::vector<float>     tree_Yc_mass;
    std::vector<float>     tree_Yc_tracks_pt;
    std::vector<float>     tree_Yc_tracks_eta;
    std::vector<float>     tree_Yc_tracks_phi;
    std::vector<math::XYZVectorF>     tree_Yc_tracks_sum3p;
    std::vector<float>     tree_Yc_tracks_InPosx;
    std::vector<float>     tree_Yc_tracks_InPosy;
    std::vector<float>     tree_Yc_tracks_InPosz;
    std::vector<float>     tree_Yc_tracks_InPx;
    std::vector<float>     tree_Yc_tracks_InPy;
    std::vector<float>     tree_Yc_tracks_InPz;
    std::vector<float>     tree_Yc_tracks_OutPx;
    std::vector<float>     tree_Yc_tracks_OutPy;
    std::vector<float>     tree_Yc_tracks_OutPz;

    std::vector<float>     tree_track_bkg_pt;
    std::vector<float>     tree_track_bkg_eta;
    std::vector<float>     tree_track_bkg_phi;
    std::vector<float>     tree_track_bkg_charge;
    std::vector<int>       tree_track_bkg_source;
    std::vector<float>     tree_track_dpt;
    std::vector<float>     tree_track_deta;
    std::vector<float>     tree_track_dphi;
    std::vector<bool>      tree_track_bkg;
    std::vector<float>     tree_Hemi_Vtx_bkg_x;
    std::vector<float>     tree_Hemi_Vtx_bkg_y;
    std::vector<float>     tree_Hemi_Vtx_bkg_z;

    std::vector<float> tree_vtx_PosX;
    std::vector<float> tree_vtx_PosY;
    std::vector<float> tree_vtx_PosZ;
    std::vector<float> tree_vtx_NChi2;
    std::vector<float> tree_vtx_PosXError;
    std::vector<float> tree_vtx_PosYError;
    std::vector<float> tree_vtx_PosZError;
    
    //--------------------------------
    // met infos -------
    //--------------------------------
    float tree_PFMet_et;
    float tree_PFMet_phi;
    float tree_PFMet_sig;
    
    //--------------------------------
    // jet infos -------
    //--------------------------------
    
    int tree_njet;
    std::vector<float> tree_jet_E;
    std::vector<float> tree_jet_pt;
    std::vector<float> tree_jet_eta;
    std::vector<float> tree_jet_phi;
    float tree_HT;
    
    //--------------------------------
    // electrons infos -------
    //--------------------------------
    std::vector<float> tree_electron_pt;
    std::vector<float> tree_electron_eta;
    std::vector<float> tree_electron_phi;
    std::vector<float> tree_electron_x;
    std::vector<float> tree_electron_y;
    std::vector<float> tree_electron_z;
    std::vector<float> tree_electron_energy;
    std::vector< int > tree_electron_charge;
    std::vector<float> tree_electron_isoR4;

    //--------------------------------
    // muons infos -------
    //--------------------------------
    float tree_Mmumu;
    std::vector<float> tree_muon_pt;
    std::vector<float> tree_muon_eta;
    std::vector<float> tree_muon_phi;
    std::vector<float> tree_muon_x;
    std::vector<float> tree_muon_y;
    std::vector<float> tree_muon_z;
    std::vector<float> tree_muon_energy;
    std::vector<float> tree_muon_dxy;
    std::vector<float> tree_muon_dxyError;
    std::vector<float> tree_muon_dz;
    std::vector<float> tree_muon_dzError;
    std::vector< int > tree_muon_charge;
    std::vector<bool>  tree_muon_isLoose;
    std::vector<bool>  tree_muon_isTight;
    std::vector<bool>  tree_muon_isGlobal;
    std::vector<float> tree_muon_isoR3;

    //-----------------------
    // per track
    //-----------------------
//     std::vector<bool>     tree_passesTrkPtr;
    std::vector<unsigned int> tree_track_ipc;
    std::vector<bool>     tree_track_lost;
    std::vector<float>    tree_track_pt;
    std::vector<float>    tree_track_eta;
    std::vector<float>    tree_track_phi;
    std::vector<int>      tree_track_charge;
    std::vector<float>    tree_track_NChi2;
    std::vector<bool>     tree_track_isHighPurity;
    std::vector<float>    tree_track_dxy; // with respect to PV
    std::vector<float>    tree_track_dxyError;
    std::vector<float>    tree_track_drSig;
    std::vector<float>    tree_track_dz;  // with respect to PV
    std::vector<float>    tree_track_dzError;
//     http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_1_3/doc/html/d8/df2/classreco_1_1TrackBase.html#aca7611bd1a33d535cefc72b6e497ece8
    std::vector<unsigned int>    tree_track_algo;
    std::vector<int>      tree_track_nHit;
    std::vector<int>      tree_track_nHitPixel;
    std::vector<int>      tree_track_nHitTIB;
    std::vector<int>      tree_track_nHitTID;
    std::vector<int>      tree_track_nHitTOB;
    std::vector<int>      tree_track_nHitTEC;
    std::vector<int>      tree_track_nHitPXB;
    std::vector<int>      tree_track_nHitPXF;
    std::vector<int>      tree_track_isHitPixel;
    std::vector<int>      tree_track_nLayers;
    std::vector<int>      tree_track_nLayersPixel;
//     std::vector<int>     tree_track_nLostHit;
    
    std::vector< float >  tree_track_x;
    std::vector< float >  tree_track_y;
    std::vector< float >  tree_track_z;
    std::vector< int >    tree_track_firstHit;
    std::vector< float >  tree_track_firstHit_x;
    std::vector< float >  tree_track_firstHit_y;
    std::vector< float >  tree_track_firstHit_z;
    std::vector< int >    tree_track_iJet;
//     std::vector<int>    tree_track_recoVertex_idx;
    std::vector<float>    tree_track_region;
    std::vector<float>    tree_track_ntrk10;
    std::vector<float>    tree_track_ntrk20;
    std::vector<float>    tree_track_ntrk30;
    std::vector<float>    tree_track_ntrk40;
    std::vector< double > tree_track_MVAval;
    
    std::vector< int >    tree_track_Hemi;
    std::vector< double > tree_track_Hemi_dR;
    std::vector< double > tree_track_Hemi_mva_NChi2;
    std::vector< bool >   tree_track_Hemi_ping;
    std::vector< float >  tree_track_Hemi_dFirstVtx;
    std::vector< int >    tree_track_Hemi_LLP;
    
    std::vector< int >    tree_track_sim_LLP;
    std::vector< bool >   tree_track_sim_isFromB;
    std::vector< bool >   tree_track_sim_isFromC;
    std::vector< float >  tree_track_sim_pt;
    std::vector< float >  tree_track_sim_eta  ;
    std::vector< float >  tree_track_sim_phi  ;
    std::vector< int >    tree_track_sim_charge;
    std::vector< int >    tree_track_sim_pdgId;
    std::vector< float >  tree_track_sim_mass  ;
    std::vector< float >  tree_track_sim_x;
    std::vector< float >  tree_track_sim_y;
    std::vector< float >  tree_track_sim_z;
    std::vector< float >  tree_track_sim_dFirstGen;
    std::vector< float >  tree_track_sim_LLP_r;
    std::vector< float >  tree_track_sim_LLP_z;
   
    //--------------------------------
    // gen infos -------
    //--------------------------------
    float tree_GenPVx;
    float tree_GenPVy;
    float tree_GenPVz;
    
    std::vector< float > tree_genParticle_pt;
    std::vector< float > tree_genParticle_eta;
    std::vector< float > tree_genParticle_phi;
    std::vector< float > tree_genParticle_charge;
    std::vector< int >   tree_genParticle_pdgId;
    std::vector< float > tree_genParticle_mass;
    std::vector< float > tree_genParticle_x;
    std::vector< float > tree_genParticle_y;
    std::vector< float > tree_genParticle_z;
    std::vector< int >   tree_genParticle_statusCode;
    std::vector< int >   tree_genParticle_mother_pdgId;
    std::vector< int >   tree_genParticle_LLP;

    int tree_ngenPackPart;
    std::vector< float > tree_genPackPart_pt;
    std::vector< float > tree_genPackPart_eta;
    std::vector< float > tree_genPackPart_phi;
    std::vector< float > tree_genPackPart_charge;
    std::vector< int >   tree_genPackPart_pdgId;
    std::vector< float > tree_genPackPart_mass;
    std::vector< float > tree_genPackPart_x;
    std::vector< float > tree_genPackPart_y;
    std::vector< float > tree_genPackPart_z;
    std::vector< int >   tree_genPackPart_mother_pdgId;
    std::vector< bool >  tree_genPackPart_isFromB;
    std::vector< bool >  tree_genPackPart_isFromC;

    int tree_ngenFromLLP;
    std::vector< int >   tree_genFromLLP_LLP;
    std::vector< float > tree_genFromLLP_pt;
    std::vector< float > tree_genFromLLP_eta;
    std::vector< float > tree_genFromLLP_phi;
    std::vector< float > tree_genFromLLP_charge;
    std::vector< int >   tree_genFromLLP_pdgId;
    std::vector< float > tree_genFromLLP_mass;
    std::vector< float > tree_genFromLLP_x;
    std::vector< float > tree_genFromLLP_y;
    std::vector< float > tree_genFromLLP_z;
    std::vector< int >   tree_genFromLLP_mother_pdgId;
    std::vector< bool >  tree_genFromLLP_isFromB;
    std::vector< bool >  tree_genFromLLP_isFromC;

    std::vector< float > tree_genAxis_dRneuneu;

    std::vector< float > tree_genFromC_pt;
    std::vector< float > tree_genFromC_eta;
    std::vector< float > tree_genFromC_phi;
    std::vector< float > tree_genFromC_charge;
    std::vector< int >   tree_genFromC_pdgId;
    std::vector< float > tree_genFromC_x;
    std::vector< float > tree_genFromC_y;
    std::vector< float > tree_genFromC_z;
    std::vector< int >   tree_genFromC_mother_pdgId;

    std::vector< float > tree_genFromB_pt;
    std::vector< float > tree_genFromB_eta;
    std::vector< float > tree_genFromB_phi;
    std::vector< float > tree_genFromB_charge;
    std::vector< int >   tree_genFromB_pdgId;
    std::vector< float > tree_genFromB_x;
    std::vector< float > tree_genFromB_y;
    std::vector< float > tree_genFromB_z;
    std::vector< int >   tree_genFromB_mother_pdgId;
   
    //--------------------------------
    // gen jet infos -------
    //--------------------------------
    std::vector<float> tree_genJet_pt;
    std::vector<float> tree_genJet_eta;
    std::vector<float> tree_genJet_phi;
    std::vector<float> tree_genJet_mass;
    std::vector<float> tree_genJet_energy;
    
    //--------------------------------
    // gen event info -------
    //--------------------------------
    
    //--------------------------------
    // lhe event infos -------
    //--------------------------------
    
    //--------------------------------
    // PF infos -------
    //--------------------------------
    
    //-----------------------
    // generated LLPs 
    //-----------------------
    int   tree_nLLP = -1;
    std::vector< int >   tree_LLP;
    std::vector< float > tree_LLP_pt;
    std::vector< float > tree_LLP_eta;
    std::vector< float > tree_LLP_phi;
    std::vector< float > tree_LLP_x;
    std::vector< float > tree_LLP_y;
    std::vector< float > tree_LLP_z;
    std::vector< float > tree_LLP_dist;
    std::vector< int >   tree_LLP_nTrks;
    std::vector< int >   tree_LLP_Vtx_nTrks;
    std::vector< float > tree_LLP_Vtx_NChi2;
    std::vector< float > tree_LLP_Vtx_dx;
    std::vector< float > tree_LLP_Vtx_dy;
    std::vector< float > tree_LLP_Vtx_dz;
    std::vector< float > tree_LLP_Vtx_dist;
    std::vector< float > tree_LLP_Vtx_dd;
    
    //-----------------------
    //Analysis with the two hemispheres
    //-----------------------
    std::vector< int >   tree_Hemi_Vtx_Layer;
    std::vector< int >   tree_Hemi_Vtx_evt2vtx;
    std::vector< int >   tree_Hemi;
    std::vector< int >   tree_Hemi_njet;
    std::vector< float > tree_Hemi_eta;
    std::vector< float > tree_Hemi_phi;
    std::vector< float > tree_Hemi_dR;
    std::vector< int >   tree_Hemi_nTrks;
    std::vector< int >   tree_Hemi_nTrks_sig;
    std::vector< int >   tree_Hemi_nTrks_bad;
    std::vector< int >   tree_Hemi_LLP;
    std::vector< float > tree_Hemi_LLP_pt;
    std::vector< float > tree_Hemi_LLP_eta;
    std::vector< float > tree_Hemi_LLP_phi;
    std::vector< float > tree_Hemi_LLP_dist;
    std::vector< float > tree_Hemi_LLP_x;
    std::vector< float > tree_Hemi_LLP_y;
    std::vector< float > tree_Hemi_LLP_z;
    std::vector< int >   tree_Hemi_Vtx_step;
    std::vector< float > tree_Hemi_Vtx_NChi2;
    std::vector< int >   tree_Hemi_Vtx_nTrks;
    std::vector< int >   tree_Hemi_Vtx_nTrks_sig;
    std::vector< int >   tree_Hemi_Vtx_nTrks_bad;
    std::vector< float > tree_Hemi_Vtx_x;
    std::vector< float > tree_Hemi_Vtx_y;
    std::vector< float > tree_Hemi_Vtx_r;
    std::vector< float > tree_Hemi_Vtx_z;
    std::vector< float > tree_Hemi_Vtx_dist;
    std::vector< float > tree_Hemi_Vtx_dx;
    std::vector< float > tree_Hemi_Vtx_dy;
    std::vector< float > tree_Hemi_Vtx_dz;
    std::vector< float > tree_Hemi_Vtx_dr;
    std::vector< float > tree_Hemi_Vtx_dd;
    std::vector< float > tree_Hemi_dR12;
    std::vector< float > tree_Hemi_LLP_dR12;
    std::vector< float > tree_Hemi_Vtx_ddbad;
    std::vector< int >   tree_Hemi_Vtx_ntrk10;
    std::vector< int >   tree_Hemi_Vtx_ntrk20;
    std::vector< float > tree_Hemi_Vtx_ddToBkg;
//     std::vector< float > tree_Hemi_Vtx_trackWeight;
    std::vector< bool >  tree_Hemi_LLP_ping;
    std::vector< int >   tree_event_LLP_ping;
    std::vector< bool >    tree_Hemi_Vtx_K0;
    std::vector< bool >    tree_Hemi_Vtx_L0;
    std::vector< bool >    tree_Hemi_Vtx_V0;
    std::vector< bool >    tree_Hemi_Vtx_Yc;

    //All preselected triggers
// ----------------Trigger Muon + dilepton-------------
    bool HLT_Mu27_Ele37_CaloIdL_MW_v;
    bool HLT_Mu37_Ele27_CaloIdL_MW_v;
    bool HLT_Mu37_TkMu27_v;
    bool HLT_Mu3_PFJet40_v;
    bool HLT_Mu7p5_L2Mu2_Jpsi_v;
    bool HLT_Mu7p5_L2Mu2_Upsilon_v;
    bool HLT_Mu7p5_Track2_Jpsi_v;
    bool HLT_Mu7p5_Track3p5_Jpsi_v;
    bool HLT_Mu7p5_Track7_Jpsi_v;
    bool HLT_Mu7p5_Track2_Upsilon_v;
    bool HLT_Mu7p5_Track3p5_Upsilon_v;
    bool HLT_Mu7p5_Track7_Upsilon_v;
    bool HLT_Mu3_L1SingleMu5orSingleMu7_v;
    bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;  // USED in 2016-2018
    bool HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_v;
    bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v;
    bool HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v;
    bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v;
    bool HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v;
    bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v;
    bool HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v;
    bool HLT_Mu25_TkMu0_Onia_v;
    bool HLT_Mu30_TkMu0_Psi_v;
    bool HLT_Mu30_TkMu0_Upsilon_v;
    bool HLT_Mu20_TkMu0_Phi_v;
    bool HLT_Mu25_TkMu0_Phi_v;
    bool HLT_Mu12_v;
    bool HLT_Mu15_v;
    bool HLT_Mu20_v;
    bool HLT_Mu27_v;
    bool HLT_Mu50_v;
    bool HLT_Mu55_v;
    bool HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_v;
    bool HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_v;
    bool HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_v;
    bool HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_v;
    bool HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v;
    bool HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v;
    bool HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v;
    bool HLT_Mu8_TrkIsoVVL_v;
    bool HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v;
    bool HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v;
    bool HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v;
    bool HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v;
    bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;  // USED in 2016-2018
    bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_v;
    bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_v;
    bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_v;
    bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_v;
    bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;     // USED in 2016-2018
    bool HLT_Mu17_TrkIsoVVL_v;
    bool HLT_Mu19_TrkIsoVVL_v;
    bool HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
    bool HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;    // USED in 2016-2018
    bool HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;
    bool HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v; // USED in 2016-2018
    bool HLT_Mu12_DoublePhoton20_v;
    bool HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v;
    bool HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_v;
    bool HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_v;
    bool HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_v;
    bool HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v;
    bool HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v;
    bool HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v;
    bool HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v;
    bool HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v;
    bool HLT_Mu15_IsoVVVL_PFHT450_v;
    bool HLT_Mu50_IsoVVVL_PFHT450_v;
    bool HLT_Mu15_IsoVVVL_PFHT600_v;
    bool HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_v;
    bool HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_v;
    bool HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_v;
    bool HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_v;
    bool HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_v;
    bool HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_v;
    bool HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_v;
    bool HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_v;
    bool HLT_Mu8_v;
    bool HLT_Mu17_v;
    bool HLT_Mu19_v;
    bool HLT_Mu17_Photon30_IsoCaloId_v;
    bool HLT_Mu18_Mu9_SameSign_v;
    bool HLT_Mu18_Mu9_SameSign_DZ_v;
    bool HLT_Mu18_Mu9_v;
    bool HLT_Mu18_Mu9_DZ_v;
    bool HLT_Mu20_Mu10_SameSign_v;
    bool HLT_Mu20_Mu10_SameSign_DZ_v;
    bool HLT_Mu20_Mu10_v;
    bool HLT_Mu20_Mu10_DZ_v;
    bool HLT_Mu23_Mu12_SameSign_v;
    bool HLT_Mu23_Mu12_SameSign_DZ_v;
    bool HLT_Mu23_Mu12_v;
    bool HLT_Mu23_Mu12_DZ_v;
    bool HLT_Mu12_IP6_part0_v;
    bool HLT_Mu12_IP6_part1_v;
    bool HLT_Mu12_IP6_part2_v;
    bool HLT_Mu12_IP6_part3_v;
    bool HLT_Mu12_IP6_part4_v;
    bool HLT_Mu9_IP5_part0_v;
    bool HLT_Mu9_IP5_part1_v;
    bool HLT_Mu9_IP5_part2_v;
    bool HLT_Mu9_IP5_part3_v;
    bool HLT_Mu9_IP5_part4_v;
    bool HLT_Mu7_IP4_part0_v;
    bool HLT_Mu7_IP4_part1_v;
    bool HLT_Mu7_IP4_part2_v;
    bool HLT_Mu7_IP4_part3_v;
    bool HLT_Mu7_IP4_part4_v;
    bool HLT_Mu9_IP4_part0_v;
    bool HLT_Mu9_IP4_part1_v;
    bool HLT_Mu9_IP4_part2_v;
    bool HLT_Mu9_IP4_part3_v;
    bool HLT_Mu9_IP4_part4_v;
    bool HLT_Mu8_IP5_part0_v;
    bool HLT_Mu8_IP5_part1_v;
    bool HLT_Mu8_IP5_part2_v;
    bool HLT_Mu8_IP5_part3_v;
    bool HLT_Mu8_IP5_part4_v;
    bool HLT_Mu8_IP6_part0_v;
    bool HLT_Mu8_IP6_part1_v;
    bool HLT_Mu8_IP6_part2_v;
    bool HLT_Mu8_IP6_part3_v;
    bool HLT_Mu8_IP6_part4_v;
    bool HLT_Mu9_IP6_part0_v;
    bool HLT_Mu9_IP6_part1_v;
    bool HLT_Mu9_IP6_part2_v;
    bool HLT_Mu9_IP6_part3_v;
    bool HLT_Mu9_IP6_part4_v;
    bool HLT_Mu8_IP3_part0_v;
    bool HLT_Mu8_IP3_part1_v;
    bool HLT_Mu8_IP3_part2_v;
    bool HLT_Mu8_IP3_part3_v;
    bool HLT_Mu8_IP3_part4_v;

// ----------------Trigger Electron-------------
    bool HLT_Ele27_Ele37_CaloIdL_MW_v;
    bool HLT_Ele20_WPTight_Gsf_v;
    bool HLT_Ele15_WPLoose_Gsf_v;
    bool HLT_Ele17_WPLoose_Gsf_v;
    bool HLT_Ele20_WPLoose_Gsf_v;
    bool HLT_Ele20_eta2p1_WPLoose_Gsf_v;
    bool HLT_Ele27_WPTight_Gsf_v;   // USED in 2016
    bool HLT_Ele28_WPTight_Gsf_v;
    bool HLT_Ele30_WPTight_Gsf_v;
    bool HLT_Ele32_WPTight_Gsf_v;   // USED in 2017-2018
    bool HLT_Ele35_WPTight_Gsf_v;
    bool HLT_Ele35_WPTight_Gsf_L1EGMT_v;
    bool HLT_Ele38_WPTight_Gsf_v;
    bool HLT_Ele40_WPTight_Gsf_v;
    bool HLT_Ele32_WPTight_Gsf_L1DoubleEG_v;
    bool HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v;
    bool HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_v;
    bool HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_v;
    bool HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v;
    bool HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v;
    bool HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v;
    bool HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v;
    bool HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v; // USED in 2016-2018
    bool HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;    // USED in 2016-2018
    bool HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v;
    bool HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v;
    bool HLT_Ele28_HighEta_SC20_Mass55_v;
    bool HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v;
    bool HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v;
    bool HLT_Ele15_IsoVVVL_PFHT450_v;
    bool HLT_Ele50_IsoVVVL_PFHT450_v;
    bool HLT_Ele15_IsoVVVL_PFHT600_v;
    bool HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v;
    bool HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v;
    bool HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v;
    bool HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v;
    bool HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v;
    bool HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v;
    bool HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v;
    bool HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v;
    bool HLT_Ele115_CaloIdVT_GsfTrkIdT_v;
    bool HLT_Ele135_CaloIdVT_GsfTrkIdT_v;
    bool HLT_Ele145_CaloIdVT_GsfTrkIdT_v;
    bool HLT_Ele200_CaloIdVT_GsfTrkIdT_v;
    bool HLT_Ele250_CaloIdVT_GsfTrkIdT_v;
    bool HLT_Ele300_CaloIdVT_GsfTrkIdT_v;
    bool HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v;

// ----------------Trigger DoubleMu-------------
    bool HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v;
    bool HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v;
    bool HLT_DoubleMu4_3_Bs_v;
    bool HLT_DoubleMu4_3_Jpsi_v;
    bool HLT_DoubleMu4_JpsiTrk_Displaced_v;
    bool HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v;
    bool HLT_DoubleMu3_Trk_Tau3mu_v;
    bool HLT_DoubleMu3_TkMu_DsTau3Mu_v;
    bool HLT_DoubleMu4_PsiPrimeTrk_Displaced_v;
    bool HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v;
    bool HLT_DoubleMu3_DZ_PFMET50_PFMHT60_v;
    bool HLT_DoubleMu3_DZ_PFMET70_PFMHT70_v;
    bool HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v;
    bool HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v;
    bool HLT_DoubleMu4_Jpsi_Displaced_v;
    bool HLT_DoubleMu4_Jpsi_NoVertexing_v;
    bool HLT_DoubleMu4_JpsiTrkTrk_Displaced_v;
    bool HLT_DoubleMu43NoFiltersNoVtx_v;
    bool HLT_DoubleMu48NoFiltersNoVtx_v;
    bool HLT_DoubleMu33NoFiltersNoVtxDisplaced_v;
    bool HLT_DoubleMu40NoFiltersNoVtxDisplaced_v;
    bool HLT_DoubleMu20_7_Mass0to30_L1_DM4_v;
    bool HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v;
    bool HLT_DoubleMu20_7_Mass0to30_Photon23_v;
    bool HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_v;
    bool HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v;
    bool HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v;

// ----------------Trigger DoubleEle-------------
    bool HLT_DoubleEle25_CaloIdL_MW_v;
    bool HLT_DoubleEle27_CaloIdL_MW_v;
    bool HLT_DoubleEle33_CaloIdL_MW_v;
    bool HLT_DoubleEle24_eta2p1_WPTight_Gsf_v;
    bool HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v;
    bool HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v;

// ----------------Trigger Dimuon0-------------
    bool HLT_Dimuon0_Jpsi_L1_NoOS_v;
    bool HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v;
    bool HLT_Dimuon0_Jpsi_v;
    bool HLT_Dimuon0_Jpsi_NoVertexing_v;
    bool HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v;
    bool HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v;
    bool HLT_Dimuon0_Jpsi3p5_Muon2_v;
    bool HLT_Dimuon0_Upsilon_L1_4p5_v;
    bool HLT_Dimuon0_Upsilon_L1_5_v;
    bool HLT_Dimuon0_Upsilon_L1_4p5NoOS_v;
    bool HLT_Dimuon0_Upsilon_L1_4p5er2p0_v;
    bool HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v;
    bool HLT_Dimuon0_Upsilon_NoVertexing_v;
    bool HLT_Dimuon0_Upsilon_L1_5M_v;
    bool HLT_Dimuon0_LowMass_L1_0er1p5R_v;
    bool HLT_Dimuon0_LowMass_L1_0er1p5_v;
    bool HLT_Dimuon0_LowMass_v;
    bool HLT_Dimuon0_LowMass_L1_4_v;
    bool HLT_Dimuon0_LowMass_L1_4R_v;
    bool HLT_Dimuon0_LowMass_L1_TM530_v;
    bool HLT_Dimuon0_Upsilon_Muon_L1_TM0_v;
    bool HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v;

// ----------------Trigger PFMET-------------
    bool HLT_PFMET110_PFMHT110_IDTight_v;
    bool HLT_PFMET120_PFMHT120_IDTight_v;   // USED
    bool HLT_PFMET130_PFMHT130_IDTight_v;
    bool HLT_PFMET140_PFMHT140_IDTight_v;
    bool HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_v;
    bool HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_v;
    bool HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_v;
    bool HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_v;
    bool HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_v;
    bool HLT_PFMET120_PFMHT120_IDTight_PFHT60_v;   // USED
    bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;   // USED
    bool HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v;
    bool HLT_PFMETTypeOne110_PFMHT110_IDTight_v;
    bool HLT_PFMETTypeOne120_PFMHT120_IDTight_v;
    bool HLT_PFMETTypeOne130_PFMHT130_IDTight_v;
    bool HLT_PFMETTypeOne140_PFMHT140_IDTight_v;
    bool HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v;
    bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;   // USED
    bool HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v;
    bool HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v;
    bool HLT_PFMET200_NotCleaned_v;
    bool HLT_PFMET200_HBHECleaned_v;
    bool HLT_PFMET250_HBHECleaned_v;   // USED
    bool HLT_PFMET300_HBHECleaned_v;
    bool HLT_PFMET200_HBHE_BeamHaloCleaned_v;
    bool HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;   // USED
    bool HLT_PFMET100_PFMHT100_IDTight_PFHT60_v;
    bool HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v;
    bool HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v;

// ----------------Trigger HT-------------
    bool HLT_HT450_Beamspot_v;
    bool HLT_HT300_Beamspot_v;
    bool HLT_HT425_v;
    bool HLT_HT430_DisplacedDijet40_DisplacedTrack_v;
    bool HLT_HT500_DisplacedDijet40_DisplacedTrack_v;
    bool HLT_HT430_DisplacedDijet60_DisplacedTrack_v;
    bool HLT_HT400_DisplacedDijet40_DisplacedTrack_v;
    bool HLT_HT650_DisplacedDijet60_Inclusive_v;
    bool HLT_HT550_DisplacedDijet60_Inclusive_v;

// ----------------Trigger AK4-------------
    bool HLT_AK4CaloJet30_v;
    bool HLT_AK4CaloJet40_v;
    bool HLT_AK4CaloJet50_v;
    bool HLT_AK4CaloJet80_v;
    bool HLT_AK4CaloJet100_v;
    bool HLT_AK4CaloJet120_v;
    bool HLT_AK4PFJet30_v;
    bool HLT_AK4PFJet50_v;
    bool HLT_AK4PFJet80_v;
    bool HLT_AK4PFJet100_v;
    bool HLT_AK4PFJet120_v;

// ----------------Trigger PFJet-------------
    bool HLT_PFJet15_v;
    bool HLT_PFJet25_v;
    bool HLT_PFJet40_v;
    bool HLT_PFJet60_v;
    bool HLT_PFJet80_v;
    bool HLT_PFJet140_v;
    bool HLT_PFJet200_v;
    bool HLT_PFJet260_v;
    bool HLT_PFJet320_v;
    bool HLT_PFJet400_v;
    bool HLT_PFJet450_v;
    bool HLT_PFJet500_v;
    bool HLT_PFJet550_v;
    bool HLT_PFJetFwd15_v;
    bool HLT_PFJetFwd25_v;
    bool HLT_PFJetFwd40_v;
    bool HLT_PFJetFwd60_v;
    bool HLT_PFJetFwd80_v;
    bool HLT_PFJetFwd140_v;
    bool HLT_PFJetFwd200_v;
    bool HLT_PFJetFwd260_v;
    bool HLT_PFJetFwd320_v;
    bool HLT_PFJetFwd400_v;
    bool HLT_PFJetFwd450_v;
    bool HLT_PFJetFwd500_v;

//-------Trigger DiPFJetAve-------//
    bool HLT_DiPFJetAve40_v;
    bool HLT_DiPFJetAve60_v;
    bool HLT_DiPFJetAve80_v;
    bool HLT_DiPFJetAve140_v;
    bool HLT_DiPFJetAve200_v;
    bool HLT_DiPFJetAve260_v;
    bool HLT_DiPFJetAve320_v;
    bool HLT_DiPFJetAve400_v;
    bool HLT_DiPFJetAve500_v;
    bool HLT_DiPFJetAve60_HFJEC_v;
    bool HLT_DiPFJetAve80_HFJEC_v;
    bool HLT_DiPFJetAve100_HFJEC_v;
    bool HLT_DiPFJetAve160_HFJEC_v;
    bool HLT_DiPFJetAve220_HFJEC_v;
    bool HLT_DiPFJetAve300_HFJEC_v;

//-------Trigger DoublePFJets-------//
    bool HLT_DoublePFJets40_CaloBTagDeepCSV_p71_v;
    bool HLT_DoublePFJets100_CaloBTagDeepCSV_p71_v;
    bool HLT_DoublePFJets200_CaloBTagDeepCSV_p71_v;
    bool HLT_DoublePFJets350_CaloBTagDeepCSV_p71_v;
    bool HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v;
    bool HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v;

//-------Trigger BTagMu-------//
    bool HLT_BTagMu_AK4DiJet20_Mu5_v;
    bool HLT_BTagMu_AK4DiJet40_Mu5_v;
    bool HLT_BTagMu_AK4DiJet70_Mu5_v;
    bool HLT_BTagMu_AK4DiJet110_Mu5_v;
    bool HLT_BTagMu_AK4DiJet170_Mu5_v;
    bool HLT_BTagMu_AK4Jet300_Mu5_v;
    bool HLT_BTagMu_AK8DiJet170_Mu5_v;
    bool HLT_BTagMu_AK8Jet170_DoubleMu5_v;
    bool HLT_BTagMu_AK8Jet300_Mu5_v;
    bool HLT_BTagMu_AK4DiJet20_Mu5_noalgo_v;
    bool HLT_BTagMu_AK4DiJet40_Mu5_noalgo_v;
    bool HLT_BTagMu_AK4DiJet70_Mu5_noalgo_v;
    bool HLT_BTagMu_AK4DiJet110_Mu5_noalgo_v;
    bool HLT_BTagMu_AK4DiJet170_Mu5_noalgo_v;
    bool HLT_BTagMu_AK4Jet300_Mu5_noalgo_v;
    bool HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v;
    bool HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_v;
    bool HLT_BTagMu_AK8Jet300_Mu5_noalgo_v;

//-------Trigger QuadPFJet-------//
    bool HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v;
    bool HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v;
    bool HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v;
    bool HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_v;
    bool HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_v;
    bool HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_v;
    bool HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_v;
    bool HLT_QuadPFJet98_83_71_15_v;
    bool HLT_QuadPFJet103_88_75_15_v;
    bool HLT_QuadPFJet105_88_76_15_v;
    bool HLT_QuadPFJet111_90_80_15_v;
    bool HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v;

//-------Trigger IsoMu-------//
    bool HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v;
    bool HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v;
    bool HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v;
    bool HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v;
    bool HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v;
    bool HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v;
    bool HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v;
    bool HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v;
    bool HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v;
    bool HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v;
    bool HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v;
    bool HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v;
    bool HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v;
    bool HLT_IsoMu20_v;
    bool HLT_IsoMu24_v;     // USED in 2016 and 2018
    bool HLT_IsoMu24_eta2p1_v;
    bool HLT_IsoMu27_v;     // USED in 2017
    bool HLT_IsoMu30_v;
    bool HLT_IsoMu24_TwoProngs35_v;
    bool HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v;
    bool HLT_IsoMu27_MET90_v;

//------------Trigger PFHT--------------//
    bool HLT_PFHT180_v;
    bool HLT_PFHT250_v;
    bool HLT_PFHT370_v;
    bool HLT_PFHT430_v;
    bool HLT_PFHT510_v;
    bool HLT_PFHT590_v;
    bool HLT_PFHT680_v;
    bool HLT_PFHT780_v;
    bool HLT_PFHT890_v;
    bool HLT_PFHT1050_v;
    bool HLT_PFHT500_PFMET100_PFMHT100_IDTight_v;   // USED
    bool HLT_PFHT500_PFMET110_PFMHT110_IDTight_v;
    bool HLT_PFHT700_PFMET85_PFMHT85_IDTight_v;   // USED
    bool HLT_PFHT700_PFMET95_PFMHT95_IDTight_v;
    bool HLT_PFHT800_PFMET75_PFMHT75_IDTight_v;   // USED
    bool HLT_PFHT800_PFMET85_PFMHT85_IDTight_v;
    bool HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v;
    bool HLT_PFHT330PT30_QuadPFJet_75_60_45_40_v;
    bool HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_v;
    bool HLT_PFHT400_SixPFJet32_v;
    bool HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_v;
    bool HLT_PFHT450_SixPFJet36_v;
    bool HLT_PFHT350_v;
    bool HLT_PFHT350MinPFJet15_v;

    // Trigger plots if needed
    TH2F* test  = new TH2F("test","test",200,0,1000,2,0,1);
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
FlyingTopAnalyzer::FlyingTopAnalyzer(const edm::ParameterSet& iConfig):

    weightFile_( iConfig.getUntrackedParameter<std::string>("weightFileMVA") ),   

    prunedGenToken_(consumes<edm::View<reco::GenParticle> >(     iConfig.getParameter<edm::InputTag>("genpruned"))),
    packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genpacked"))),
    vertexToken_(   consumes<reco::VertexCollection>(            iConfig.getParameter<edm::InputTag>("vertices"))),
    metToken_(      consumes<pat::METCollection>(                iConfig.getParameter<edm::InputTag>("mets"))),
    jetToken_(      consumes<pat::JetCollection>(                iConfig.getParameter<edm::InputTag>("jets"))),
    genJetToken_(   consumes<edm::View<reco::GenJet>>(           iConfig.getParameter<edm::InputTag>("genjets"))),
    electronToken_( consumes<pat::ElectronCollection>(           iConfig.getParameter<edm::InputTag>("electrons"))),
    muonToken_(     consumes<pat::MuonCollection>(               iConfig.getParameter<edm::InputTag>("muons"))),
    pcToken_(       consumes<pat::PackedCandidateCollection>(    iConfig.getParameter<edm::InputTag>("pfCands"))),
    lostpcToken_(   consumes<pat::PackedCandidateCollection>(    iConfig.getParameter<edm::InputTag>("lostpfCands"))), //LOST
    triggerResultsToken_(consumes<edm::TriggerResults>(edm::InputTag(std::string("TriggerResults"),std::string(""),std::string("HLT"))) )//trig
    ,K0Token_(      consumes<reco::VertexCompositePtrCandidateCollection>(            iConfig.getParameter<edm::InputTag>("Kshorts"))),  
    LambdaToken_(   consumes<reco::VertexCompositePtrCandidateCollection>(            iConfig.getParameter<edm::InputTag>("Lambda")))
    , SVToken_(   consumes<reco::VertexCompositePtrCandidateCollection>(            iConfig.getParameter<edm::InputTag>("SV")))
    ,PhotonToken_(  consumes<reco::ConversionCollection>(edm::InputTag(std::string("reducedEgamma"),std::string("reducedConversions")))) //,std::string("PAT")TTbar: "PAT" _____ Neu: "RECO" 
    // , PrescaleToken_( consumes<pat::PackedTriggerPrescales>(edm::InputTag(std::string("patTrigger"),std::string("")))  )
{
   //now do what ever initialization is needed
    nEvent = 0;
    usesResource("TFileService");
    
    smalltree = fs->make<TTree>("ttree", "ttree");
    
    // event info
    smalltree->Branch("runNumber",  &runNumber,  "runNumber/I");
    smalltree->Branch("eventNumber",&eventNumber,"eventNumber/I");
    smalltree->Branch("lumiBlock"  ,&lumiBlock,  "lumiBlock/I");
    
    // primary vertex info
    smalltree->Branch("tree_nPV", &tree_nPV);
    smalltree->Branch("tree_PV_x",     &tree_PV_x);
    smalltree->Branch("tree_PV_y",     &tree_PV_y);
    smalltree->Branch("tree_PV_z",     &tree_PV_z);    
    smalltree->Branch("tree_PV_ez",    &tree_PV_ez);    
    smalltree->Branch("tree_PV_NChi2", &tree_PV_NChi2);    
    smalltree->Branch("tree_PV_ndf",   &tree_PV_ndf);

    smalltree->Branch("tree_nSV",&tree_nSV);
    smalltree->Branch("tree_SV_x",&tree_SV_x);
    smalltree->Branch("tree_SV_y",&tree_SV_y);
    smalltree->Branch("tree_SV_z",&tree_SV_z);
    smalltree->Branch("tree_SV_r",&tree_SV_r);
    smalltree->Branch("tree_SV_NChi2",&tree_SV_NChi2);
    smalltree->Branch("tree_SV_ndf",&tree_SV_ndf);
    smalltree->Branch("tree_SV_nDaughters",&tree_SV_nDaughters);
    smalltree->Branch("tree_SV_mass",&tree_SV_mass);
    smalltree->Branch("tree_SV_pdgid",&tree_SV_pdgid);
    smalltree->Branch("tree_SV_pt",&tree_SV_pt);
    smalltree->Branch("tree_SV_eta",&tree_SV_eta);
    smalltree->Branch("tree_SV_phi",&tree_SV_phi);
    smalltree->Branch("tree_SV_IsConvertedPhoton",&tree_SV_IsConvertedPhoton);
//$$$$    smalltree->Branch("tree_SV_daughter_pdgid",&tree_SV_daughter_pdgid);

    smalltree->Branch("tree_nK0",           &tree_nK0);
    smalltree->Branch("tree_K0_x",          &tree_K0_x); // l'index 0 donne le PV!
    smalltree->Branch("tree_K0_y",          &tree_K0_y);
    smalltree->Branch("tree_K0_z",          &tree_K0_z);
    smalltree->Branch("tree_K0_r",          &tree_K0_r);
    smalltree->Branch("tree_K0_NChi2",      &tree_K0_NChi2);
    smalltree->Branch("tree_K0_ndf",        &tree_K0_ndf);
    smalltree->Branch("tree_K0_mass",       &tree_K0_mass);
    smalltree->Branch("tree_K0_pt",         &tree_K0_pt);
    smalltree->Branch("tree_K0_eta",        &tree_K0_eta);
    smalltree->Branch("tree_K0_phi",        &tree_K0_phi);
    smalltree->Branch("tree_K0_nDaughters", &tree_K0_nDaughters);
    smalltree->Branch("tree_K0_daughter_pt",&tree_K0_daughter_pt);
    smalltree->Branch("tree_K0_daughter_eta",&tree_K0_daughter_eta);
    smalltree->Branch("tree_K0_daughter_phi",&tree_K0_daughter_phi);

    smalltree->Branch("tree_nLambda",       &tree_nLambda);
    smalltree->Branch("tree_L0_x",          &tree_L0_x); // l'index 0 donne le PV!
    smalltree->Branch("tree_L0_y",          &tree_L0_y);
    smalltree->Branch("tree_L0_z",          &tree_L0_z);
    smalltree->Branch("tree_L0_r",          &tree_L0_r);
    smalltree->Branch("tree_L0_NChi2",      &tree_L0_NChi2);
    smalltree->Branch("tree_L0_ndf",        &tree_L0_ndf);
    smalltree->Branch("tree_L0_nDaughters", &tree_L0_nDaughters);
    smalltree->Branch("tree_L0_mass",       &tree_L0_mass);
    smalltree->Branch("tree_L0_pt",         &tree_L0_pt);
    smalltree->Branch("tree_L0_eta",        &tree_L0_eta);
    smalltree->Branch("tree_L0_phi",        &tree_L0_phi);
    smalltree->Branch("tree_L0_daughter_pt",&tree_L0_daughter_pt);
    smalltree->Branch("tree_L0_daughter_eta",&tree_L0_daughter_eta);
    smalltree->Branch("tree_L0_daughter_phi",&tree_L0_daughter_phi);

    smalltree->Branch("tree_nYConv",        &tree_nYConv);
    smalltree->Branch("tree_Yc_x",          &tree_Yc_x); 
    smalltree->Branch("tree_Yc_y",          &tree_Yc_y);
    smalltree->Branch("tree_Yc_z",          &tree_Yc_z);
    smalltree->Branch("tree_Yc_r",          &tree_Yc_r);
    smalltree->Branch("tree_Yc_Vtx_NChi2",  &tree_Yc_Vtx_NChi2);
    smalltree->Branch("tree_Yc_Vtx_ndf",    &tree_Yc_Vtx_ndf);
    smalltree->Branch("tree_Yc_mass",       &tree_Yc_mass);
    smalltree->Branch("tree_Yc_tracks_pt",  &tree_Yc_tracks_pt);
    smalltree->Branch("tree_Yc_tracks_eta", &tree_Yc_tracks_eta);
    smalltree->Branch("tree_Yc_tracks_phi", &tree_Yc_tracks_phi);
    smalltree->Branch("tree_Yc_tracks_sum3p",&tree_Yc_tracks_sum3p);
    smalltree->Branch("tree_Yc_tracks_InPosx",&tree_Yc_tracks_InPosx);
    smalltree->Branch("tree_Yc_tracks_InPosy",&tree_Yc_tracks_InPosy);
    smalltree->Branch("tree_Yc_tracks_InPosz",&tree_Yc_tracks_InPosz);
    smalltree->Branch("tree_Yc_tracks_InPx",&tree_Yc_tracks_InPx);
    smalltree->Branch("tree_Yc_tracks_InPy",&tree_Yc_tracks_InPy);
    smalltree->Branch("tree_Yc_tracks_InPz",&tree_Yc_tracks_InPz);
    smalltree->Branch("tree_Yc_tracks_OutPx",&tree_Yc_tracks_OutPx);
    smalltree->Branch("tree_Yc_tracks_OutPy",&tree_Yc_tracks_OutPy);
    smalltree->Branch("tree_Yc_tracks_OutPz",&tree_Yc_tracks_OutPz);

    smalltree->Branch("tree_track_bkg_pt",   &tree_track_bkg_pt);
    smalltree->Branch("tree_track_bkg_eta",  &tree_track_bkg_eta);
    smalltree->Branch("tree_track_bkg_phi",  &tree_track_bkg_phi);
    smalltree->Branch("tree_track_bkg_charge",&tree_track_bkg_charge);
    smalltree->Branch("tree_track_bkg_source",&tree_track_bkg_source);
    smalltree->Branch("tree_track_dpt",&tree_track_dpt);
    smalltree->Branch("tree_track_deta",&tree_track_deta);
    smalltree->Branch("tree_track_dphi",&tree_track_dphi);
    smalltree->Branch("tree_track_bkg",  &tree_track_bkg);
    smalltree->Branch("tree_Hemi_Vtx_bkg_x",&tree_Hemi_Vtx_bkg_x);
    smalltree->Branch("tree_Hemi_Vtx_bkg_y",&tree_Hemi_Vtx_bkg_y);
    smalltree->Branch("tree_Hemi_Vtx_bkg_z",&tree_Hemi_Vtx_bkg_z);
    
        
//     smalltree->Branch("tree_vtx_PosX", &tree_vtx_PosX);
//     smalltree->Branch("tree_vtx_PosY", &tree_vtx_PosY);
//     smalltree->Branch("tree_vtx_PosZ", &tree_vtx_PosZ);
//     smalltree->Branch("tree_vtx_NChi2", &tree_vtx_NChi2);
//     smalltree->Branch("tree_vtx_PosXError", &tree_vtx_PosXError);
//     smalltree->Branch("tree_vtx_PosYError", &tree_vtx_PosYError);
//     smalltree->Branch("tree_vtx_PosZError", &tree_vtx_PosZError);
    
//     // trigger info
// //     smalltree->Branch("tree_trigger_names", &tree_trigger_names);
// //     smalltree->Branch("tree_trigger_bits",  &tree_trigger_bits);
//     smalltree->Branch("tree_trigger_size", &tree_trigger_size);
//     smalltree->Branch("tree_passesTrigger", &tree_passesTrigger);
//     smalltree->Branch("tree_passesTriggerName", &tree_passesTriggerName);
//     smalltree->Branch("tree_Trigger_Muon",&tree_Trigger_Muon);//+ dilepton channel emu
//     smalltree->Branch("tree_Trigger_Ele",&tree_Trigger_Ele);
//     smalltree->Branch("tree_Trigger_DoubleMu",&tree_Trigger_DoubleMu);
//     smalltree->Branch("tree_Trigger_DoubleEle",&tree_Trigger_DoubleEle);
//     smalltree->Branch("tree_Trigger_Dimuon0",&tree_Trigger_Dimuon0);
//     smalltree->Branch("tree_Trigger_PFMET",&tree_Trigger_PFMET);
//     smalltree->Branch("tree_Trigger_HT",&tree_Trigger_HT);
//     smalltree->Branch("tree_Trigger_AK4",&tree_Trigger_AK4);
//     smalltree->Branch("tree_Trigger_PFJet",&tree_Trigger_PFJet);
//     smalltree->Branch("tree_Trigger_DoublePFJets",&tree_Trigger_DoublePFJets);
//     smalltree->Branch("tree_Trigger_DiPFJet",&tree_Trigger_DiPFJet);
//     smalltree->Branch("tree_Trigger_QuadPFJet",&tree_Trigger_QuadPFJet);
//     smalltree->Branch("tree_Trigger_BTagMu",&tree_Trigger_BTagMu);

    smalltree->Branch("tree_NbrOfZCand",  &tree_NbrOfZCand,  "tree_NbrOfZCand/I");
    smalltree->Branch("tree_Filter", &tree_Filter);
    
    // met info
    smalltree->Branch("tree_PFMet_et" ,  &tree_PFMet_et);
    smalltree->Branch("tree_PFMet_phi" , &tree_PFMet_phi);
    smalltree->Branch("tree_PFMet_sig" , &tree_PFMet_sig);
    
    // jet info
    smalltree->Branch("tree_njet"  ,        &tree_njet);
    smalltree->Branch("tree_jet_E"  ,       &tree_jet_E);
    smalltree->Branch("tree_jet_pt"  ,      &tree_jet_pt);
    smalltree->Branch("tree_jet_eta" ,      &tree_jet_eta);
    smalltree->Branch("tree_jet_phi" ,      &tree_jet_phi);
    smalltree->Branch("tree_HT"  ,          &tree_HT);
    
    // electrons info
    smalltree->Branch("tree_electron_pt"  ,   &tree_electron_pt);
    smalltree->Branch("tree_electron_eta" ,   &tree_electron_eta);
    smalltree->Branch("tree_electron_phi" ,   &tree_electron_phi);
    smalltree->Branch("tree_electron_x"  ,    &tree_electron_x);
    smalltree->Branch("tree_electron_y" ,     &tree_electron_y);
    smalltree->Branch("tree_electron_z" ,     &tree_electron_z);
    smalltree->Branch("tree_electron_energy", &tree_electron_energy);
    smalltree->Branch("tree_electron_charge", &tree_electron_charge);
    smalltree->Branch("tree_electron_isoR4",&tree_electron_isoR4);

    // muons info
    smalltree->Branch("tree_Mmumu"  ,         &tree_Mmumu);
    smalltree->Branch("tree_muon_pt"  ,       &tree_muon_pt);
    smalltree->Branch("tree_muon_eta" ,       &tree_muon_eta);
    smalltree->Branch("tree_muon_phi" ,       &tree_muon_phi);
    smalltree->Branch("tree_muon_x"  ,        &tree_muon_x);
    smalltree->Branch("tree_muon_y" ,         &tree_muon_y);
    smalltree->Branch("tree_muon_z" ,         &tree_muon_z);
    smalltree->Branch("tree_muon_energy",     &tree_muon_energy);
    smalltree->Branch("tree_muon_dxy",        &tree_muon_dxy);
    smalltree->Branch("tree_muon_dxyError",   &tree_muon_dxyError);
    smalltree->Branch("tree_muon_dz",         &tree_muon_dz);
    smalltree->Branch("tree_muon_dzError",    &tree_muon_dzError);
    smalltree->Branch("tree_muon_charge",     &tree_muon_charge);
    smalltree->Branch("tree_muon_isLoose",    &tree_muon_isLoose);
    smalltree->Branch("tree_muon_isTight",    &tree_muon_isTight);
    smalltree->Branch("tree_muon_isGlobal",   &tree_muon_isGlobal);
    smalltree->Branch("tree_muon_isoR3",&tree_muon_isoR3);
    // track
    smalltree->Branch("tree_nTracks",                &tree_nTracks, "tree_nTracks/I"); 
    smalltree->Branch("tree_nLostTracks",            &tree_nLostTracks, "tree_nLostTracks/I"); 
//     smalltree->Branch("tree_passesTrkPtr",           &tree_passesTrkPtr);
    smalltree->Branch("tree_track_ipc",              &tree_track_ipc);
    smalltree->Branch("tree_track_lost",             &tree_track_lost);
    smalltree->Branch("tree_track_pt",               &tree_track_pt);
    smalltree->Branch("tree_track_eta",              &tree_track_eta );
    smalltree->Branch("tree_track_phi",              &tree_track_phi );
    smalltree->Branch("tree_track_charge",           &tree_track_charge );
    smalltree->Branch("tree_track_NChi2",            &tree_track_NChi2);
    smalltree->Branch("tree_track_isHighPurity",     &tree_track_isHighPurity);
    smalltree->Branch("tree_track_dxy",              &tree_track_dxy );
    smalltree->Branch("tree_track_dxyError",         &tree_track_dxyError);
    smalltree->Branch("tree_track_drSig",            &tree_track_drSig);
    smalltree->Branch("tree_track_dz",               &tree_track_dz);
    smalltree->Branch("tree_track_dzError",          &tree_track_dzError  );
    smalltree->Branch("tree_track_algo",             &tree_track_algo);
    smalltree->Branch("tree_track_nHit",         &tree_track_nHit);
    smalltree->Branch("tree_track_nHitPixel",    &tree_track_nHitPixel);
    smalltree->Branch("tree_track_nHitTIB",      &tree_track_nHitTIB);
    smalltree->Branch("tree_track_nHitTID",      &tree_track_nHitTID);
    smalltree->Branch("tree_track_nHitTOB",      &tree_track_nHitTOB);
    smalltree->Branch("tree_track_nHitTEC",      &tree_track_nHitTEC);
    smalltree->Branch("tree_track_nHitPXB",      &tree_track_nHitPXB);
    smalltree->Branch("tree_track_nHitPXF",      &tree_track_nHitPXF);
    smalltree->Branch("tree_track_isHitPixel",   &tree_track_isHitPixel);
    smalltree->Branch("tree_track_nLayers",      &tree_track_nLayers);
    smalltree->Branch("tree_track_nLayersPixel", &tree_track_nLayersPixel);
    
    smalltree->Branch("tree_track_x",            &tree_track_x );
    smalltree->Branch("tree_track_y",            &tree_track_y );
    smalltree->Branch("tree_track_z",            &tree_track_z );
    smalltree->Branch("tree_track_firstHit",     &tree_track_firstHit);
    smalltree->Branch("tree_track_region",       &tree_track_region);
    smalltree->Branch("tree_track_firstHit_x",   &tree_track_firstHit_x);
    smalltree->Branch("tree_track_firstHit_y",   &tree_track_firstHit_y);
    smalltree->Branch("tree_track_firstHit_z",   &tree_track_firstHit_z);
    smalltree->Branch("tree_track_iJet",         &tree_track_iJet);
    smalltree->Branch("tree_track_ntrk10",       &tree_track_ntrk10);
    smalltree->Branch("tree_track_ntrk20",       &tree_track_ntrk20);
    smalltree->Branch("tree_track_ntrk30",       &tree_track_ntrk30);
    smalltree->Branch("tree_track_ntrk40",       &tree_track_ntrk40);
    smalltree->Branch("tree_track_MVAval",       &tree_track_MVAval);
    
    smalltree->Branch("tree_track_Hemi",           &tree_track_Hemi);
    smalltree->Branch("tree_track_Hemi_dR",        &tree_track_Hemi_dR);
    smalltree->Branch("tree_track_Hemi_mva_NChi2", &tree_track_Hemi_mva_NChi2);
    smalltree->Branch("tree_track_Hemi_ping",      &tree_track_Hemi_ping);
    smalltree->Branch("tree_track_Hemi_dFirstVtx", &tree_track_Hemi_dFirstVtx);
    smalltree->Branch("tree_track_Hemi_LLP",       &tree_track_Hemi_LLP);
        
    // info about the simulated track from LLP matched to the reco track
    smalltree->Branch("tree_track_sim_LLP",        &tree_track_sim_LLP );
    smalltree->Branch("tree_track_sim_isFromB",    &tree_track_sim_isFromB );
    smalltree->Branch("tree_track_sim_isFromC",    &tree_track_sim_isFromC );
    smalltree->Branch("tree_track_sim_pt",         &tree_track_sim_pt );
    smalltree->Branch("tree_track_sim_eta",        &tree_track_sim_eta );
    smalltree->Branch("tree_track_sim_phi",        &tree_track_sim_phi );
    smalltree->Branch("tree_track_sim_charge",     &tree_track_sim_charge );
    smalltree->Branch("tree_track_sim_pdgId",      &tree_track_sim_pdgId );
    smalltree->Branch("tree_track_sim_mass",       &tree_track_sim_mass );
    smalltree->Branch("tree_track_sim_x",          &tree_track_sim_x );
    smalltree->Branch("tree_track_sim_y",          &tree_track_sim_y );
    smalltree->Branch("tree_track_sim_z",          &tree_track_sim_z );    
    smalltree->Branch("tree_track_sim_dFirstGen",  &tree_track_sim_dFirstGen );
    smalltree->Branch("tree_track_sim_LLP_r",      &tree_track_sim_LLP_r );
    smalltree->Branch("tree_track_sim_LLP_z",      &tree_track_sim_LLP_z );

    // gen info
    smalltree->Branch("tree_GenPVx" ,  &tree_GenPVx);
    smalltree->Branch("tree_GenPVy" ,  &tree_GenPVy);
    smalltree->Branch("tree_GenPVz" ,  &tree_GenPVz);
    
    smalltree->Branch("tree_genParticle_pt"  ,          &tree_genParticle_pt);
    smalltree->Branch("tree_genParticle_eta" ,          &tree_genParticle_eta);
    smalltree->Branch("tree_genParticle_phi" ,          &tree_genParticle_phi);
    smalltree->Branch("tree_genParticle_charge" ,       &tree_genParticle_charge);
    smalltree->Branch("tree_genParticle_pdgId" ,        &tree_genParticle_pdgId);
    smalltree->Branch("tree_genParticle_mass" ,         &tree_genParticle_mass);
    smalltree->Branch("tree_genParticle_x"  ,	        &tree_genParticle_x);
    smalltree->Branch("tree_genParticle_y" ,	        &tree_genParticle_y);
    smalltree->Branch("tree_genParticle_z" ,	        &tree_genParticle_z);
    smalltree->Branch("tree_genParticle_statusCode",    &tree_genParticle_statusCode);
    smalltree->Branch("tree_genParticle_mother_pdgId" , &tree_genParticle_mother_pdgId);
    smalltree->Branch("tree_genParticle_LLP" ,          &tree_genParticle_LLP);

    smalltree->Branch("tree_ngenPackPart"  ,            &tree_ngenPackPart);
    smalltree->Branch("tree_genPackPart_pt"  ,          &tree_genPackPart_pt);
    smalltree->Branch("tree_genPackPart_eta" ,          &tree_genPackPart_eta);
    smalltree->Branch("tree_genPackPart_phi" ,          &tree_genPackPart_phi);
    smalltree->Branch("tree_genPackPart_charge" ,       &tree_genPackPart_charge);
    smalltree->Branch("tree_genPackPart_pdgId" ,        &tree_genPackPart_pdgId);
    smalltree->Branch("tree_genPackPart_mass" ,         &tree_genPackPart_mass);
    smalltree->Branch("tree_genPackPart_x"  ,	        &tree_genPackPart_x);
    smalltree->Branch("tree_genPackPart_y" ,	        &tree_genPackPart_y);
    smalltree->Branch("tree_genPackPart_z" ,	        &tree_genPackPart_z);
    smalltree->Branch("tree_genPackPart_mother_pdgId" , &tree_genPackPart_mother_pdgId);
    smalltree->Branch("tree_genPackPart_isFromB" ,      &tree_genPackPart_isFromB);
    smalltree->Branch("tree_genPackPart_isFromC" ,      &tree_genPackPart_isFromC);

    smalltree->Branch("tree_ngenFromLLP"  ,            &tree_ngenFromLLP);
    smalltree->Branch("tree_genFromLLP_LLP"  ,         &tree_genFromLLP_LLP);
    smalltree->Branch("tree_genFromLLP_pt"  ,          &tree_genFromLLP_pt);
    smalltree->Branch("tree_genFromLLP_eta" ,          &tree_genFromLLP_eta);
    smalltree->Branch("tree_genFromLLP_phi" ,          &tree_genFromLLP_phi);
    smalltree->Branch("tree_genFromLLP_charge" ,       &tree_genFromLLP_charge);
    smalltree->Branch("tree_genFromLLP_pdgId" ,        &tree_genFromLLP_pdgId);
    smalltree->Branch("tree_genFromLLP_mass" ,         &tree_genFromLLP_mass);
    smalltree->Branch("tree_genFromLLP_x"  ,	       &tree_genFromLLP_x);
    smalltree->Branch("tree_genFromLLP_y" ,	       &tree_genFromLLP_y);
    smalltree->Branch("tree_genFromLLP_z" ,	       &tree_genFromLLP_z);
    smalltree->Branch("tree_genFromLLP_mother_pdgId" , &tree_genFromLLP_mother_pdgId);
    smalltree->Branch("tree_genFromLLP_isFromB" ,      &tree_genFromLLP_isFromB);
    smalltree->Branch("tree_genFromLLP_isFromC" ,      &tree_genFromLLP_isFromC);

    smalltree->Branch("tree_genAxis_dRneuneu",       &tree_genAxis_dRneuneu);

    smalltree->Branch("tree_nFromC",                 &tree_nFromC,  "tree_nFromC/I");
    smalltree->Branch("tree_genFromC_pt"  ,          &tree_genFromC_pt);
    smalltree->Branch("tree_genFromC_eta" ,          &tree_genFromC_eta);
    smalltree->Branch("tree_genFromC_phi" ,          &tree_genFromC_phi);
    smalltree->Branch("tree_genFromC_charge" ,       &tree_genFromC_charge);
    smalltree->Branch("tree_genFromC_pdgId" ,        &tree_genFromC_pdgId);
    smalltree->Branch("tree_genFromC_x"  ,	     &tree_genFromC_x);
    smalltree->Branch("tree_genFromC_y" ,	     &tree_genFromC_y);
    smalltree->Branch("tree_genFromC_z" ,	     &tree_genFromC_z);
    smalltree->Branch("tree_genFromC_mother_pdgId" , &tree_genFromC_mother_pdgId);

    smalltree->Branch("tree_nFromB",                 &tree_nFromB,  "tree_nFromB/I");
    smalltree->Branch("tree_genFromB_pt"  ,	     &tree_genFromB_pt);
    smalltree->Branch("tree_genFromB_eta" ,	     &tree_genFromB_eta);
    smalltree->Branch("tree_genFromB_phi" ,	     &tree_genFromB_phi);
    smalltree->Branch("tree_genFromB_charge" ,	     &tree_genFromB_charge);
    smalltree->Branch("tree_genFromB_pdgId" ,	     &tree_genFromB_pdgId);
    smalltree->Branch("tree_genFromB_x"  ,	     &tree_genFromB_x);
    smalltree->Branch("tree_genFromB_y" ,	     &tree_genFromB_y);
    smalltree->Branch("tree_genFromB_z" ,	     &tree_genFromB_z);
    smalltree->Branch("tree_genFromB_mother_pdgId" , &tree_genFromB_mother_pdgId);
    
    // genJet info
    smalltree->Branch("tree_genJet_pt"  ,   &tree_genJet_pt);
    smalltree->Branch("tree_genJet_eta" ,   &tree_genJet_eta);
    smalltree->Branch("tree_genJet_phi" ,   &tree_genJet_phi);
    smalltree->Branch("tree_genJet_mass",   &tree_genJet_mass);
    smalltree->Branch("tree_genJet_energy", &tree_genJet_energy);
    
    smalltree->Branch("tree_nLLP",          &tree_nLLP);
    smalltree->Branch("tree_LLP",           &tree_LLP);
    smalltree->Branch("tree_LLP_pt" ,       &tree_LLP_pt);
    smalltree->Branch("tree_LLP_eta",       &tree_LLP_eta);
    smalltree->Branch("tree_LLP_phi",       &tree_LLP_phi);
    smalltree->Branch("tree_LLP_x",         &tree_LLP_x);
    smalltree->Branch("tree_LLP_y",         &tree_LLP_y);
    smalltree->Branch("tree_LLP_z",         &tree_LLP_z);
    smalltree->Branch("tree_LLP_dist",      &tree_LLP_dist);
    smalltree->Branch("tree_LLP_nTrks",     &tree_LLP_nTrks);
    smalltree->Branch("tree_LLP_Vtx_nTrks", &tree_LLP_Vtx_nTrks);
    smalltree->Branch("tree_LLP_Vtx_NChi2", &tree_LLP_Vtx_NChi2);
    smalltree->Branch("tree_LLP_Vtx_dx",    &tree_LLP_Vtx_dx);
    smalltree->Branch("tree_LLP_Vtx_dy",    &tree_LLP_Vtx_dy);
    smalltree->Branch("tree_LLP_Vtx_dz",    &tree_LLP_Vtx_dz);
    smalltree->Branch("tree_LLP_Vtx_dist",  &tree_LLP_Vtx_dist);
    smalltree->Branch("tree_LLP_Vtx_dd",    &tree_LLP_Vtx_dd);

    smalltree->Branch("tree_Hemi_Vtx_Layer",&tree_Hemi_Vtx_Layer);
    smalltree->Branch("tree_Hemi_Vtx_evt2vtx",&tree_Hemi_Vtx_evt2vtx);
    smalltree->Branch("tree_Hemi",       &tree_Hemi);
    smalltree->Branch("tree_Hemi_njet",  &tree_Hemi_njet);
    smalltree->Branch("tree_Hemi_eta",   &tree_Hemi_eta);
    smalltree->Branch("tree_Hemi_phi",   &tree_Hemi_phi);
    smalltree->Branch("tree_Hemi_dR",    &tree_Hemi_dR);
    smalltree->Branch("tree_Hemi_nTrks", &tree_Hemi_nTrks);
    smalltree->Branch("tree_Hemi_nTrks_sig", &tree_Hemi_nTrks_sig);
    smalltree->Branch("tree_Hemi_nTrks_bad", &tree_Hemi_nTrks_bad);
    smalltree->Branch("tree_Hemi_LLP",       &tree_Hemi_LLP);
    smalltree->Branch("tree_Hemi_LLP_pt",    &tree_Hemi_LLP_pt);
    smalltree->Branch("tree_Hemi_LLP_eta",   &tree_Hemi_LLP_eta);
    smalltree->Branch("tree_Hemi_LLP_phi",   &tree_Hemi_LLP_phi);
    smalltree->Branch("tree_Hemi_LLP_dist",  &tree_Hemi_LLP_dist);
    smalltree->Branch("tree_Hemi_LLP_x",     &tree_Hemi_LLP_x);
    smalltree->Branch("tree_Hemi_LLP_y",     &tree_Hemi_LLP_y);
    smalltree->Branch("tree_Hemi_LLP_z",     &tree_Hemi_LLP_z);
    smalltree->Branch("tree_Hemi_Vtx_step",  &tree_Hemi_Vtx_step);
    smalltree->Branch("tree_Hemi_Vtx_NChi2", &tree_Hemi_Vtx_NChi2);
    smalltree->Branch("tree_Hemi_Vtx_nTrks", &tree_Hemi_Vtx_nTrks);
    smalltree->Branch("tree_Hemi_Vtx_nTrks_sig", &tree_Hemi_Vtx_nTrks_sig);
    smalltree->Branch("tree_Hemi_Vtx_nTrks_bad", &tree_Hemi_Vtx_nTrks_bad);
    smalltree->Branch("tree_Hemi_Vtx_x",     &tree_Hemi_Vtx_x);
    smalltree->Branch("tree_Hemi_Vtx_y",     &tree_Hemi_Vtx_y);
    smalltree->Branch("tree_Hemi_Vtx_r",     &tree_Hemi_Vtx_r);
    smalltree->Branch("tree_Hemi_Vtx_z",     &tree_Hemi_Vtx_z);
    smalltree->Branch("tree_Hemi_Vtx_dist",  &tree_Hemi_Vtx_dist);
    smalltree->Branch("tree_Hemi_Vtx_dx",    &tree_Hemi_Vtx_dx);
    smalltree->Branch("tree_Hemi_Vtx_dy",    &tree_Hemi_Vtx_dy);
    smalltree->Branch("tree_Hemi_Vtx_dz",    &tree_Hemi_Vtx_dz);
    smalltree->Branch("tree_Hemi_Vtx_dr",    &tree_Hemi_Vtx_dr);
    smalltree->Branch("tree_Hemi_Vtx_dd",    &tree_Hemi_Vtx_dd);
    smalltree->Branch("tree_Hemi_dR12",      &tree_Hemi_dR12);
    smalltree->Branch("tree_Hemi_LLP_dR12",  &tree_Hemi_LLP_dR12);
    smalltree->Branch("tree_Hemi_Vtx_ddbad", &tree_Hemi_Vtx_ddbad);
    smalltree->Branch("tree_Hemi_Vtx_ntrk10",&tree_Hemi_Vtx_ntrk10);
    smalltree->Branch("tree_Hemi_Vtx_ntrk20",&tree_Hemi_Vtx_ntrk20);
    smalltree->Branch("tree_Hemi_Vtx_ddToBkg",&tree_Hemi_Vtx_ddToBkg);
//     smalltree->Branch("tree_Hemi_Vtx_trackWeight", &tree_Hemi_Vtx_trackWeight);
    smalltree->Branch("tree_Hemi_LLP_ping",  &tree_Hemi_LLP_ping);
    smalltree->Branch("tree_event_LLP_ping", &tree_event_LLP_ping);
    smalltree->Branch("tree_Hemi_Vtx_K0",&tree_Hemi_Vtx_K0);
    smalltree->Branch("tree_Hemi_Vtx_L0",&tree_Hemi_Vtx_L0);
    smalltree->Branch("tree_Hemi_Vtx_V0",&tree_Hemi_Vtx_V0);
    smalltree->Branch("tree_Hemi_Vtx_Yc",&tree_Hemi_Vtx_Yc);
// // ----------------Trigger Muon + dilepton-------------
// smalltree->Branch("HLT_Mu27_Ele37_CaloIdL_MW_v",&HLT_Mu27_Ele37_CaloIdL_MW_v);
// smalltree->Branch("HLT_Mu37_Ele27_CaloIdL_MW_v",&HLT_Mu37_Ele27_CaloIdL_MW_v);
// smalltree->Branch("HLT_Mu37_TkMu27_v",&HLT_Mu37_TkMu27_v);
// smalltree->Branch("HLT_Mu3_PFJet40_v",&HLT_Mu3_PFJet40_v);
// smalltree->Branch("HLT_Mu7p5_L2Mu2_Jpsi_v",&HLT_Mu7p5_L2Mu2_Jpsi_v);
// smalltree->Branch("HLT_Mu7p5_L2Mu2_Upsilon_v",&HLT_Mu7p5_L2Mu2_Upsilon_v);
// smalltree->Branch("HLT_Mu7p5_Track2_Jpsi_v",&HLT_Mu7p5_Track2_Jpsi_v);
// smalltree->Branch("HLT_Mu7p5_Track3p5_Jpsi_v",&HLT_Mu7p5_Track3p5_Jpsi_v);
// smalltree->Branch("HLT_Mu7p5_Track7_Jpsi_v",&HLT_Mu7p5_Track7_Jpsi_v);
// smalltree->Branch("HLT_Mu7p5_Track2_Upsilon_v",&HLT_Mu7p5_Track2_Upsilon_v);
// smalltree->Branch("HLT_Mu7p5_Track3p5_Upsilon_v",&HLT_Mu7p5_Track3p5_Upsilon_v);
// smalltree->Branch("HLT_Mu7p5_Track7_Upsilon_v",&HLT_Mu7p5_Track7_Upsilon_v);
// smalltree->Branch("HLT_Mu3_L1SingleMu5orSingleMu7_v",&HLT_Mu3_L1SingleMu5orSingleMu7_v);
smalltree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
// smalltree->Branch("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_v",&HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_v);
// smalltree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v);
// smalltree->Branch("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v",&HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v);
// smalltree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v);
// smalltree->Branch("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v",&HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v);
// smalltree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v);
// smalltree->Branch("HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v",&HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v);
// smalltree->Branch("HLT_Mu25_TkMu0_Onia_v",&HLT_Mu25_TkMu0_Onia_v);
// smalltree->Branch("HLT_Mu30_TkMu0_Psi_v",&HLT_Mu30_TkMu0_Psi_v);
// smalltree->Branch("HLT_Mu30_TkMu0_Upsilon_v",&HLT_Mu30_TkMu0_Upsilon_v);
// smalltree->Branch("HLT_Mu20_TkMu0_Phi_v",&HLT_Mu20_TkMu0_Phi_v);
// smalltree->Branch("HLT_Mu25_TkMu0_Phi_v",&HLT_Mu25_TkMu0_Phi_v);
// smalltree->Branch("HLT_Mu12_v",&HLT_Mu12_v);
// smalltree->Branch("HLT_Mu15_v",&HLT_Mu15_v);
// smalltree->Branch("HLT_Mu20_v",&HLT_Mu20_v);
// smalltree->Branch("HLT_Mu27_v",&HLT_Mu27_v);
// smalltree->Branch("HLT_Mu50_v",&HLT_Mu50_v);
// smalltree->Branch("HLT_Mu55_v",&HLT_Mu55_v);
// smalltree->Branch("HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_v",&HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_v",&HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_v",&HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_v",&HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v",&HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v",&HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v",&HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_Mu8_TrkIsoVVL_v",&HLT_Mu8_TrkIsoVVL_v);
// smalltree->Branch("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v",&HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v);
// smalltree->Branch("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v",&HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v);
// smalltree->Branch("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v",&HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v);
// smalltree->Branch("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v",&HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v);
smalltree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
// smalltree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_v);
// smalltree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_v);
// smalltree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_v);
// smalltree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_v);
smalltree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
// smalltree->Branch("HLT_Mu17_TrkIsoVVL_v",&HLT_Mu17_TrkIsoVVL_v);
// smalltree->Branch("HLT_Mu19_TrkIsoVVL_v",&HLT_Mu19_TrkIsoVVL_v);
// smalltree->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
smalltree->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
// smalltree->Branch("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
smalltree->Branch("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
// smalltree->Branch("HLT_Mu12_DoublePhoton20_v",&HLT_Mu12_DoublePhoton20_v);
// smalltree->Branch("HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v",&HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v);
// smalltree->Branch("HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_v",&HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_v);
// smalltree->Branch("HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_v",&HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_v);
// smalltree->Branch("HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_v",&HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_v);
// smalltree->Branch("HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v",&HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v);
// smalltree->Branch("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v",&HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v);
// smalltree->Branch("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v",&HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v);
// smalltree->Branch("HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v",&HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v);
// smalltree->Branch("HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v",&HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v);
// smalltree->Branch("HLT_Mu15_IsoVVVL_PFHT450_v",&HLT_Mu15_IsoVVVL_PFHT450_v);
// smalltree->Branch("HLT_Mu50_IsoVVVL_PFHT450_v",&HLT_Mu50_IsoVVVL_PFHT450_v);
// smalltree->Branch("HLT_Mu15_IsoVVVL_PFHT600_v",&HLT_Mu15_IsoVVVL_PFHT600_v);
// smalltree->Branch("HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_v",&HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_v);
// smalltree->Branch("HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_v",&HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_v);
// smalltree->Branch("HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_v",&HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_v);
// smalltree->Branch("HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_v",&HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_v);
// smalltree->Branch("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_v",&HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_v);
// smalltree->Branch("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_v",&HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_v);
// smalltree->Branch("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_v",&HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_v);
// smalltree->Branch("HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_v",&HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_v);
// smalltree->Branch("HLT_Mu8_v",&HLT_Mu8_v);
// smalltree->Branch("HLT_Mu17_v",&HLT_Mu17_v);
// smalltree->Branch("HLT_Mu19_v",&HLT_Mu19_v);
// smalltree->Branch("HLT_Mu17_Photon30_IsoCaloId_v",&HLT_Mu17_Photon30_IsoCaloId_v);
// smalltree->Branch("HLT_Mu18_Mu9_SameSign_v",&HLT_Mu18_Mu9_SameSign_v);
// smalltree->Branch("HLT_Mu18_Mu9_SameSign_DZ_v",&HLT_Mu18_Mu9_SameSign_DZ_v);
// smalltree->Branch("HLT_Mu18_Mu9_v",&HLT_Mu18_Mu9_v);
// smalltree->Branch("HLT_Mu18_Mu9_DZ_v",&HLT_Mu18_Mu9_DZ_v);
// smalltree->Branch("HLT_Mu20_Mu10_SameSign_v",&HLT_Mu20_Mu10_SameSign_v);
// smalltree->Branch("HLT_Mu20_Mu10_SameSign_DZ_v",&HLT_Mu20_Mu10_SameSign_DZ_v);
// smalltree->Branch("HLT_Mu20_Mu10_v",&HLT_Mu20_Mu10_v);
// smalltree->Branch("HLT_Mu20_Mu10_DZ_v",&HLT_Mu20_Mu10_DZ_v);
// smalltree->Branch("HLT_Mu23_Mu12_SameSign_v",&HLT_Mu23_Mu12_SameSign_v);
// smalltree->Branch("HLT_Mu23_Mu12_SameSign_DZ_v",&HLT_Mu23_Mu12_SameSign_DZ_v);
// smalltree->Branch("HLT_Mu23_Mu12_v",&HLT_Mu23_Mu12_v);
// smalltree->Branch("HLT_Mu23_Mu12_DZ_v",&HLT_Mu23_Mu12_DZ_v);
// smalltree->Branch("HLT_Mu12_IP6_part0_v",&HLT_Mu12_IP6_part0_v);
// smalltree->Branch("HLT_Mu12_IP6_part1_v",&HLT_Mu12_IP6_part1_v);
// smalltree->Branch("HLT_Mu12_IP6_part2_v",&HLT_Mu12_IP6_part2_v);
// smalltree->Branch("HLT_Mu12_IP6_part3_v",&HLT_Mu12_IP6_part3_v);
// smalltree->Branch("HLT_Mu12_IP6_part4_v",&HLT_Mu12_IP6_part4_v);
// smalltree->Branch("HLT_Mu9_IP5_part0_v",&HLT_Mu9_IP5_part0_v);
// smalltree->Branch("HLT_Mu9_IP5_part1_v",&HLT_Mu9_IP5_part1_v);
// smalltree->Branch("HLT_Mu9_IP5_part2_v",&HLT_Mu9_IP5_part2_v);
// smalltree->Branch("HLT_Mu9_IP5_part3_v",&HLT_Mu9_IP5_part3_v);
// smalltree->Branch("HLT_Mu9_IP5_part4_v",&HLT_Mu9_IP5_part4_v);
// smalltree->Branch("HLT_Mu7_IP4_part0_v",&HLT_Mu7_IP4_part0_v);
// smalltree->Branch("HLT_Mu7_IP4_part1_v",&HLT_Mu7_IP4_part1_v);
// smalltree->Branch("HLT_Mu7_IP4_part2_v",&HLT_Mu7_IP4_part2_v);
// smalltree->Branch("HLT_Mu7_IP4_part3_v",&HLT_Mu7_IP4_part3_v);
// smalltree->Branch("HLT_Mu7_IP4_part4_v",&HLT_Mu7_IP4_part4_v);
// smalltree->Branch("HLT_Mu9_IP4_part0_v",&HLT_Mu9_IP4_part0_v);
// smalltree->Branch("HLT_Mu9_IP4_part1_v",&HLT_Mu9_IP4_part1_v);
// smalltree->Branch("HLT_Mu9_IP4_part2_v",&HLT_Mu9_IP4_part2_v);
// smalltree->Branch("HLT_Mu9_IP4_part3_v",&HLT_Mu9_IP4_part3_v);
// smalltree->Branch("HLT_Mu9_IP4_part4_v",&HLT_Mu9_IP4_part4_v);
// smalltree->Branch("HLT_Mu8_IP5_part0_v",&HLT_Mu8_IP5_part0_v);
// smalltree->Branch("HLT_Mu8_IP5_part1_v",&HLT_Mu8_IP5_part1_v);
// smalltree->Branch("HLT_Mu8_IP5_part2_v",&HLT_Mu8_IP5_part2_v);
// smalltree->Branch("HLT_Mu8_IP5_part3_v",&HLT_Mu8_IP5_part3_v);
// smalltree->Branch("HLT_Mu8_IP5_part4_v",&HLT_Mu8_IP5_part4_v);
// smalltree->Branch("HLT_Mu8_IP6_part0_v",&HLT_Mu8_IP6_part0_v);
// smalltree->Branch("HLT_Mu8_IP6_part1_v",&HLT_Mu8_IP6_part1_v);
// smalltree->Branch("HLT_Mu8_IP6_part2_v",&HLT_Mu8_IP6_part2_v);
// smalltree->Branch("HLT_Mu8_IP6_part3_v",&HLT_Mu8_IP6_part3_v);
// smalltree->Branch("HLT_Mu8_IP6_part4_v",&HLT_Mu8_IP6_part4_v);
// smalltree->Branch("HLT_Mu9_IP6_part0_v",&HLT_Mu9_IP6_part0_v);
// smalltree->Branch("HLT_Mu9_IP6_part1_v",&HLT_Mu9_IP6_part1_v);
// smalltree->Branch("HLT_Mu9_IP6_part2_v",&HLT_Mu9_IP6_part2_v);
// smalltree->Branch("HLT_Mu9_IP6_part3_v",&HLT_Mu9_IP6_part3_v);
// smalltree->Branch("HLT_Mu9_IP6_part4_v",&HLT_Mu9_IP6_part4_v);
// smalltree->Branch("HLT_Mu8_IP3_part0_v",&HLT_Mu8_IP3_part0_v);
// smalltree->Branch("HLT_Mu8_IP3_part1_v",&HLT_Mu8_IP3_part1_v);
// smalltree->Branch("HLT_Mu8_IP3_part2_v",&HLT_Mu8_IP3_part2_v);
// smalltree->Branch("HLT_Mu8_IP3_part3_v",&HLT_Mu8_IP3_part3_v);
// smalltree->Branch("HLT_Mu8_IP3_part4_v",&HLT_Mu8_IP3_part4_v);
// // ----------------Trigger Electron-------------
// smalltree->Branch("HLT_Ele27_Ele37_CaloIdL_MW_v",&HLT_Ele27_Ele37_CaloIdL_MW_v);
// smalltree->Branch("HLT_Ele20_WPTight_Gsf_v",&HLT_Ele20_WPTight_Gsf_v);
// smalltree->Branch("HLT_Ele15_WPLoose_Gsf_v",&HLT_Ele15_WPLoose_Gsf_v);
// smalltree->Branch("HLT_Ele17_WPLoose_Gsf_v",&HLT_Ele17_WPLoose_Gsf_v);
// smalltree->Branch("HLT_Ele20_WPLoose_Gsf_v",&HLT_Ele20_WPLoose_Gsf_v);
// smalltree->Branch("HLT_Ele20_eta2p1_WPLoose_Gsf_v",&HLT_Ele20_eta2p1_WPLoose_Gsf_v);
smalltree->Branch("HLT_Ele27_WPTight_Gsf_v",&HLT_Ele27_WPTight_Gsf_v);
// smalltree->Branch("HLT_Ele28_WPTight_Gsf_v",&HLT_Ele28_WPTight_Gsf_v);
// smalltree->Branch("HLT_Ele30_WPTight_Gsf_v",&HLT_Ele30_WPTight_Gsf_v);
smalltree->Branch("HLT_Ele32_WPTight_Gsf_v",&HLT_Ele32_WPTight_Gsf_v);
// smalltree->Branch("HLT_Ele35_WPTight_Gsf_v",&HLT_Ele35_WPTight_Gsf_v);
// smalltree->Branch("HLT_Ele35_WPTight_Gsf_L1EGMT_v",&HLT_Ele35_WPTight_Gsf_L1EGMT_v);
// smalltree->Branch("HLT_Ele38_WPTight_Gsf_v",&HLT_Ele38_WPTight_Gsf_v);
// smalltree->Branch("HLT_Ele40_WPTight_Gsf_v",&HLT_Ele40_WPTight_Gsf_v);
// smalltree->Branch("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v",&HLT_Ele32_WPTight_Gsf_L1DoubleEG_v);
// smalltree->Branch("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v",&HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v);
// smalltree->Branch("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_v",&HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_v);
// smalltree->Branch("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_v",&HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_v);
// smalltree->Branch("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v",&HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v);
// smalltree->Branch("HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v",&HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v);
// smalltree->Branch("HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v",&HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v);
// smalltree->Branch("HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v",&HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v);
smalltree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
smalltree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);
// smalltree->Branch("HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v",&HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v);
// smalltree->Branch("HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v",&HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v);
// smalltree->Branch("HLT_Ele28_HighEta_SC20_Mass55_v",&HLT_Ele28_HighEta_SC20_Mass55_v);
// smalltree->Branch("HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v",&HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v);
// smalltree->Branch("HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v",&HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v);
// smalltree->Branch("HLT_Ele15_IsoVVVL_PFHT450_v",&HLT_Ele15_IsoVVVL_PFHT450_v);
// smalltree->Branch("HLT_Ele50_IsoVVVL_PFHT450_v",&HLT_Ele50_IsoVVVL_PFHT450_v);
// smalltree->Branch("HLT_Ele15_IsoVVVL_PFHT600_v",&HLT_Ele15_IsoVVVL_PFHT600_v);
// smalltree->Branch("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v",&HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v);
// smalltree->Branch("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v",&HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v);
// smalltree->Branch("HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v",&HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v);
// smalltree->Branch("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v",&HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v);
// smalltree->Branch("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v",&HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v);
// smalltree->Branch("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v",&HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v);
// smalltree->Branch("HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v",&HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v);
// smalltree->Branch("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v",&HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v);
// smalltree->Branch("HLT_Ele115_CaloIdVT_GsfTrkIdT_v",&HLT_Ele115_CaloIdVT_GsfTrkIdT_v);
// smalltree->Branch("HLT_Ele135_CaloIdVT_GsfTrkIdT_v",&HLT_Ele135_CaloIdVT_GsfTrkIdT_v);
// smalltree->Branch("HLT_Ele145_CaloIdVT_GsfTrkIdT_v",&HLT_Ele145_CaloIdVT_GsfTrkIdT_v);
// smalltree->Branch("HLT_Ele200_CaloIdVT_GsfTrkIdT_v",&HLT_Ele200_CaloIdVT_GsfTrkIdT_v);
// smalltree->Branch("HLT_Ele250_CaloIdVT_GsfTrkIdT_v",&HLT_Ele250_CaloIdVT_GsfTrkIdT_v);
// smalltree->Branch("HLT_Ele300_CaloIdVT_GsfTrkIdT_v",&HLT_Ele300_CaloIdVT_GsfTrkIdT_v);
// smalltree->Branch("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v",&HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v);
// // ----------------Trigger DoubleMu-------------
// smalltree->Branch("HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v",&HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v);
// smalltree->Branch("HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v",&HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v);
// smalltree->Branch("HLT_DoubleMu4_3_Bs_v",&HLT_DoubleMu4_3_Bs_v);
// smalltree->Branch("HLT_DoubleMu4_3_Jpsi_v",&HLT_DoubleMu4_3_Jpsi_v);
// smalltree->Branch("HLT_DoubleMu4_JpsiTrk_Displaced_v",&HLT_DoubleMu4_JpsiTrk_Displaced_v);
// smalltree->Branch("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v",&HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v);
// smalltree->Branch("HLT_DoubleMu3_Trk_Tau3mu_v",&HLT_DoubleMu3_Trk_Tau3mu_v);
// smalltree->Branch("HLT_DoubleMu3_TkMu_DsTau3Mu_v",&HLT_DoubleMu3_TkMu_DsTau3Mu_v);
// smalltree->Branch("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v",&HLT_DoubleMu4_PsiPrimeTrk_Displaced_v);
// smalltree->Branch("HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v",&HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v);
// smalltree->Branch("HLT_DoubleMu3_DZ_PFMET50_PFMHT60_v",&HLT_DoubleMu3_DZ_PFMET50_PFMHT60_v);
// smalltree->Branch("HLT_DoubleMu3_DZ_PFMET70_PFMHT70_v",&HLT_DoubleMu3_DZ_PFMET70_PFMHT70_v);
// smalltree->Branch("HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v",&HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v);
// smalltree->Branch("HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v",&HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v);
// smalltree->Branch("HLT_DoubleMu4_Jpsi_Displaced_v",&HLT_DoubleMu4_Jpsi_Displaced_v);
// smalltree->Branch("HLT_DoubleMu4_Jpsi_NoVertexing_v",&HLT_DoubleMu4_Jpsi_NoVertexing_v);
// smalltree->Branch("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v",&HLT_DoubleMu4_JpsiTrkTrk_Displaced_v);
// smalltree->Branch("HLT_DoubleMu43NoFiltersNoVtx_v",&HLT_DoubleMu43NoFiltersNoVtx_v);
// smalltree->Branch("HLT_DoubleMu48NoFiltersNoVtx_v",&HLT_DoubleMu48NoFiltersNoVtx_v);
// smalltree->Branch("HLT_DoubleMu33NoFiltersNoVtxDisplaced_v",&HLT_DoubleMu33NoFiltersNoVtxDisplaced_v);
// smalltree->Branch("HLT_DoubleMu40NoFiltersNoVtxDisplaced_v",&HLT_DoubleMu40NoFiltersNoVtxDisplaced_v);
// smalltree->Branch("HLT_DoubleMu20_7_Mass0to30_L1_DM4_v",&HLT_DoubleMu20_7_Mass0to30_L1_DM4_v);
// smalltree->Branch("HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v",&HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v);
// smalltree->Branch("HLT_DoubleMu20_7_Mass0to30_Photon23_v",&HLT_DoubleMu20_7_Mass0to30_Photon23_v);
// smalltree->Branch("HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_v",&HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_v);
// smalltree->Branch("HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v",&HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v);
// smalltree->Branch("HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v",&HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v);
// // ----------------Trigger DoubleEle-------------
// smalltree->Branch("HLT_DoubleEle25_CaloIdL_MW_v",&HLT_DoubleEle25_CaloIdL_MW_v);
// smalltree->Branch("HLT_DoubleEle27_CaloIdL_MW_v",&HLT_DoubleEle27_CaloIdL_MW_v);
// smalltree->Branch("HLT_DoubleEle33_CaloIdL_MW_v",&HLT_DoubleEle33_CaloIdL_MW_v);
// smalltree->Branch("HLT_DoubleEle24_eta2p1_WPTight_Gsf_v",&HLT_DoubleEle24_eta2p1_WPTight_Gsf_v);
// smalltree->Branch("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v",&HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v);
// smalltree->Branch("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v",&HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v);
// // ----------------Trigger Dimuon0-------------
// smalltree->Branch("HLT_Dimuon0_Jpsi_L1_NoOS_v",&HLT_Dimuon0_Jpsi_L1_NoOS_v);
// smalltree->Branch("HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v",&HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v);
// smalltree->Branch("HLT_Dimuon0_Jpsi_v",&HLT_Dimuon0_Jpsi_v);
// smalltree->Branch("HLT_Dimuon0_Jpsi_NoVertexing_v",&HLT_Dimuon0_Jpsi_NoVertexing_v);
// smalltree->Branch("HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v",&HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v);
// smalltree->Branch("HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v",&HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v);
// smalltree->Branch("HLT_Dimuon0_Jpsi3p5_Muon2_v",&HLT_Dimuon0_Jpsi3p5_Muon2_v);
// smalltree->Branch("HLT_Dimuon0_Upsilon_L1_4p5_v",&HLT_Dimuon0_Upsilon_L1_4p5_v);
// smalltree->Branch("HLT_Dimuon0_Upsilon_L1_5_v",&HLT_Dimuon0_Upsilon_L1_5_v);
// smalltree->Branch("HLT_Dimuon0_Upsilon_L1_4p5NoOS_v",&HLT_Dimuon0_Upsilon_L1_4p5NoOS_v);
// smalltree->Branch("HLT_Dimuon0_Upsilon_L1_4p5er2p0_v",&HLT_Dimuon0_Upsilon_L1_4p5er2p0_v);
// smalltree->Branch("HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v",&HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v);
// smalltree->Branch("HLT_Dimuon0_Upsilon_NoVertexing_v",&HLT_Dimuon0_Upsilon_NoVertexing_v);
// smalltree->Branch("HLT_Dimuon0_Upsilon_L1_5M_v",&HLT_Dimuon0_Upsilon_L1_5M_v);
// smalltree->Branch("HLT_Dimuon0_LowMass_L1_0er1p5R_v",&HLT_Dimuon0_LowMass_L1_0er1p5R_v);
// smalltree->Branch("HLT_Dimuon0_LowMass_L1_0er1p5_v",&HLT_Dimuon0_LowMass_L1_0er1p5_v);
// smalltree->Branch("HLT_Dimuon0_LowMass_v",&HLT_Dimuon0_LowMass_v);
// smalltree->Branch("HLT_Dimuon0_LowMass_L1_4_v",&HLT_Dimuon0_LowMass_L1_4_v);
// smalltree->Branch("HLT_Dimuon0_LowMass_L1_4R_v",&HLT_Dimuon0_LowMass_L1_4R_v);
// smalltree->Branch("HLT_Dimuon0_LowMass_L1_TM530_v",&HLT_Dimuon0_LowMass_L1_TM530_v);
// smalltree->Branch("HLT_Dimuon0_Upsilon_Muon_L1_TM0_v",&HLT_Dimuon0_Upsilon_Muon_L1_TM0_v);
// smalltree->Branch("HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v",&HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v);
// // ----------------Trigger PFMET-------------
// smalltree->Branch("HLT_PFMET110_PFMHT110_IDTight_v",&HLT_PFMET110_PFMHT110_IDTight_v);
smalltree->Branch("HLT_PFMET120_PFMHT120_IDTight_v",&HLT_PFMET120_PFMHT120_IDTight_v);
// smalltree->Branch("HLT_PFMET130_PFMHT130_IDTight_v",&HLT_PFMET130_PFMHT130_IDTight_v);
// smalltree->Branch("HLT_PFMET140_PFMHT140_IDTight_v",&HLT_PFMET140_PFMHT140_IDTight_v);
// smalltree->Branch("HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_v",&HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_v);
// smalltree->Branch("HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_v",&HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_v);
// smalltree->Branch("HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_v",&HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_v);
// smalltree->Branch("HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_v",&HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_v);
// smalltree->Branch("HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_v",&HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_v);
smalltree->Branch("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",&HLT_PFMET120_PFMHT120_IDTight_PFHT60_v);
smalltree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v);
// smalltree->Branch("HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v",&HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v);
// smalltree->Branch("HLT_PFMETTypeOne110_PFMHT110_IDTight_v",&HLT_PFMETTypeOne110_PFMHT110_IDTight_v);
// smalltree->Branch("HLT_PFMETTypeOne120_PFMHT120_IDTight_v",&HLT_PFMETTypeOne120_PFMHT120_IDTight_v);
// smalltree->Branch("HLT_PFMETTypeOne130_PFMHT130_IDTight_v",&HLT_PFMETTypeOne130_PFMHT130_IDTight_v);
// smalltree->Branch("HLT_PFMETTypeOne140_PFMHT140_IDTight_v",&HLT_PFMETTypeOne140_PFMHT140_IDTight_v);
// smalltree->Branch("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v",&HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v);
smalltree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
// smalltree->Branch("HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v",&HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v);
// smalltree->Branch("HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v",&HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v);
// smalltree->Branch("HLT_PFMET200_NotCleaned_v",&HLT_PFMET200_NotCleaned_v);
// smalltree->Branch("HLT_PFMET200_HBHECleaned_v",&HLT_PFMET200_HBHECleaned_v);
smalltree->Branch("HLT_PFMET250_HBHECleaned_v",&HLT_PFMET250_HBHECleaned_v);
// smalltree->Branch("HLT_PFMET300_HBHECleaned_v",&HLT_PFMET300_HBHECleaned_v);
// smalltree->Branch("HLT_PFMET200_HBHE_BeamHaloCleaned_v",&HLT_PFMET200_HBHE_BeamHaloCleaned_v);
smalltree->Branch("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",&HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v);
// smalltree->Branch("HLT_PFMET100_PFMHT100_IDTight_PFHT60_v",&HLT_PFMET100_PFMHT100_IDTight_PFHT60_v);
// smalltree->Branch("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v",&HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v);
// smalltree->Branch("HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v",&HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v);
// // ----------------Trigger HT-------------
// smalltree->Branch("HLT_HT450_Beamspot_v",&HLT_HT450_Beamspot_v);
// smalltree->Branch("HLT_HT300_Beamspot_v",&HLT_HT300_Beamspot_v);
// smalltree->Branch("HLT_HT425_v",&HLT_HT425_v);
// smalltree->Branch("HLT_HT430_DisplacedDijet40_DisplacedTrack_v",&HLT_HT430_DisplacedDijet40_DisplacedTrack_v);
// smalltree->Branch("HLT_HT500_DisplacedDijet40_DisplacedTrack_v",&HLT_HT500_DisplacedDijet40_DisplacedTrack_v);
// smalltree->Branch("HLT_HT430_DisplacedDijet60_DisplacedTrack_v",&HLT_HT430_DisplacedDijet60_DisplacedTrack_v);
// smalltree->Branch("HLT_HT400_DisplacedDijet40_DisplacedTrack_v",&HLT_HT400_DisplacedDijet40_DisplacedTrack_v);
// smalltree->Branch("HLT_HT650_DisplacedDijet60_Inclusive_v",&HLT_HT650_DisplacedDijet60_Inclusive_v);
// smalltree->Branch("HLT_HT550_DisplacedDijet60_Inclusive_v",&HLT_HT550_DisplacedDijet60_Inclusive_v);
// // ----------------Trigger AK4-------------
// smalltree->Branch("HLT_AK4CaloJet30_v",&HLT_AK4CaloJet30_v);
// smalltree->Branch("HLT_AK4CaloJet40_v",&HLT_AK4CaloJet40_v);
// smalltree->Branch("HLT_AK4CaloJet50_v",&HLT_AK4CaloJet50_v);
// smalltree->Branch("HLT_AK4CaloJet80_v",&HLT_AK4CaloJet80_v);
// smalltree->Branch("HLT_AK4CaloJet100_v",&HLT_AK4CaloJet100_v);
// smalltree->Branch("HLT_AK4CaloJet120_v",&HLT_AK4CaloJet120_v);
// smalltree->Branch("HLT_AK4PFJet30_v",&HLT_AK4PFJet30_v);
// smalltree->Branch("HLT_AK4PFJet50_v",&HLT_AK4PFJet50_v);
// smalltree->Branch("HLT_AK4PFJet80_v",&HLT_AK4PFJet80_v);
// smalltree->Branch("HLT_AK4PFJet100_v",&HLT_AK4PFJet100_v);
// smalltree->Branch("HLT_AK4PFJet120_v",&HLT_AK4PFJet120_v);
// // ----------------Trigger PFJet-------------
// smalltree->Branch("HLT_PFJet15_v",&HLT_PFJet15_v);
// smalltree->Branch("HLT_PFJet25_v",&HLT_PFJet25_v);
// smalltree->Branch("HLT_PFJet40_v",&HLT_PFJet40_v);
// smalltree->Branch("HLT_PFJet60_v",&HLT_PFJet60_v);
// smalltree->Branch("HLT_PFJet80_v",&HLT_PFJet80_v);
// smalltree->Branch("HLT_PFJet140_v",&HLT_PFJet140_v);
// smalltree->Branch("HLT_PFJet200_v",&HLT_PFJet200_v);
// smalltree->Branch("HLT_PFJet260_v",&HLT_PFJet260_v);
// smalltree->Branch("HLT_PFJet320_v",&HLT_PFJet320_v);
// smalltree->Branch("HLT_PFJet400_v",&HLT_PFJet400_v);
// smalltree->Branch("HLT_PFJet450_v",&HLT_PFJet450_v);
// smalltree->Branch("HLT_PFJet500_v",&HLT_PFJet500_v);
// smalltree->Branch("HLT_PFJet550_v",&HLT_PFJet550_v);
// smalltree->Branch("HLT_PFJetFwd15_v",&HLT_PFJetFwd15_v);
// smalltree->Branch("HLT_PFJetFwd25_v",&HLT_PFJetFwd25_v);
// smalltree->Branch("HLT_PFJetFwd40_v",&HLT_PFJetFwd40_v);
// smalltree->Branch("HLT_PFJetFwd60_v",&HLT_PFJetFwd60_v);
// smalltree->Branch("HLT_PFJetFwd80_v",&HLT_PFJetFwd80_v);
// smalltree->Branch("HLT_PFJetFwd140_v",&HLT_PFJetFwd140_v);
// smalltree->Branch("HLT_PFJetFwd200_v",&HLT_PFJetFwd200_v);
// smalltree->Branch("HLT_PFJetFwd260_v",&HLT_PFJetFwd260_v);
// smalltree->Branch("HLT_PFJetFwd320_v",&HLT_PFJetFwd320_v);
// smalltree->Branch("HLT_PFJetFwd400_v",&HLT_PFJetFwd400_v);
// smalltree->Branch("HLT_PFJetFwd450_v",&HLT_PFJetFwd450_v);
// smalltree->Branch("HLT_PFJetFwd500_v",&HLT_PFJetFwd500_v);
// // ----------------Trigger DiPFJet-------------
// smalltree->Branch("HLT_DiPFJetAve40_v",&HLT_DiPFJetAve40_v);
// smalltree->Branch("HLT_DiPFJetAve60_v",&HLT_DiPFJetAve60_v);
// smalltree->Branch("HLT_DiPFJetAve80_v",&HLT_DiPFJetAve80_v);
// smalltree->Branch("HLT_DiPFJetAve140_v",&HLT_DiPFJetAve140_v);
// smalltree->Branch("HLT_DiPFJetAve200_v",&HLT_DiPFJetAve200_v);
// smalltree->Branch("HLT_DiPFJetAve260_v",&HLT_DiPFJetAve260_v);
// smalltree->Branch("HLT_DiPFJetAve320_v",&HLT_DiPFJetAve320_v);
// smalltree->Branch("HLT_DiPFJetAve400_v",&HLT_DiPFJetAve400_v);
// smalltree->Branch("HLT_DiPFJetAve500_v",&HLT_DiPFJetAve500_v);
// smalltree->Branch("HLT_DiPFJetAve60_HFJEC_v",&HLT_DiPFJetAve60_HFJEC_v);
// smalltree->Branch("HLT_DiPFJetAve80_HFJEC_v",&HLT_DiPFJetAve80_HFJEC_v);
// smalltree->Branch("HLT_DiPFJetAve100_HFJEC_v",&HLT_DiPFJetAve100_HFJEC_v);
// smalltree->Branch("HLT_DiPFJetAve160_HFJEC_v",&HLT_DiPFJetAve160_HFJEC_v);
// smalltree->Branch("HLT_DiPFJetAve220_HFJEC_v",&HLT_DiPFJetAve220_HFJEC_v);
// smalltree->Branch("HLT_DiPFJetAve300_HFJEC_v",&HLT_DiPFJetAve300_HFJEC_v);
// // ----------------Trigger DoublePFJet-------------
// smalltree->Branch("HLT_DoublePFJets40_CaloBTagDeepCSV_p71_v",&HLT_DoublePFJets40_CaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_DoublePFJets100_CaloBTagDeepCSV_p71_v",&HLT_DoublePFJets100_CaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_DoublePFJets200_CaloBTagDeepCSV_p71_v",&HLT_DoublePFJets200_CaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_DoublePFJets350_CaloBTagDeepCSV_p71_v",&HLT_DoublePFJets350_CaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v",&HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v);
// smalltree->Branch("HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v",&HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v);
// // ----------------Trigger BTagMu-------------
// smalltree->Branch("HLT_BTagMu_AK4DiJet20_Mu5_v",&HLT_BTagMu_AK4DiJet20_Mu5_v);
// smalltree->Branch("HLT_BTagMu_AK4DiJet40_Mu5_v",&HLT_BTagMu_AK4DiJet40_Mu5_v);
// smalltree->Branch("HLT_BTagMu_AK4DiJet70_Mu5_v",&HLT_BTagMu_AK4DiJet70_Mu5_v);
// smalltree->Branch("HLT_BTagMu_AK4DiJet110_Mu5_v",&HLT_BTagMu_AK4DiJet110_Mu5_v);
// smalltree->Branch("HLT_BTagMu_AK4DiJet170_Mu5_v",&HLT_BTagMu_AK4DiJet170_Mu5_v);
// smalltree->Branch("HLT_BTagMu_AK4Jet300_Mu5_v",&HLT_BTagMu_AK4Jet300_Mu5_v);
// smalltree->Branch("HLT_BTagMu_AK8DiJet170_Mu5_v",&HLT_BTagMu_AK8DiJet170_Mu5_v);
// smalltree->Branch("HLT_BTagMu_AK8Jet170_DoubleMu5_v",&HLT_BTagMu_AK8Jet170_DoubleMu5_v);
// smalltree->Branch("HLT_BTagMu_AK8Jet300_Mu5_v",&HLT_BTagMu_AK8Jet300_Mu5_v);
// smalltree->Branch("HLT_BTagMu_AK4DiJet20_Mu5_noalgo_v",&HLT_BTagMu_AK4DiJet20_Mu5_noalgo_v);
// smalltree->Branch("HLT_BTagMu_AK4DiJet40_Mu5_noalgo_v",&HLT_BTagMu_AK4DiJet40_Mu5_noalgo_v);
// smalltree->Branch("HLT_BTagMu_AK4DiJet70_Mu5_noalgo_v",&HLT_BTagMu_AK4DiJet70_Mu5_noalgo_v);
// smalltree->Branch("HLT_BTagMu_AK4DiJet110_Mu5_noalgo_v",&HLT_BTagMu_AK4DiJet110_Mu5_noalgo_v);
// smalltree->Branch("HLT_BTagMu_AK4DiJet170_Mu5_noalgo_v",&HLT_BTagMu_AK4DiJet170_Mu5_noalgo_v);
// smalltree->Branch("HLT_BTagMu_AK4Jet300_Mu5_noalgo_v",&HLT_BTagMu_AK4Jet300_Mu5_noalgo_v);
// smalltree->Branch("HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v",&HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v);
// smalltree->Branch("HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_v",&HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_v);
// smalltree->Branch("HLT_BTagMu_AK8Jet300_Mu5_noalgo_v",&HLT_BTagMu_AK8Jet300_Mu5_noalgo_v);
// // ----------------Trigger QuadPFJet-------------
// smalltree->Branch("HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v",&HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v);
// smalltree->Branch("HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v",&HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v);
// smalltree->Branch("HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v",&HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v);
// smalltree->Branch("HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_v",&HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_v);
// smalltree->Branch("HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_v",&HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_v);
// smalltree->Branch("HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_v",&HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_v);
// smalltree->Branch("HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_v",&HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_v);
// smalltree->Branch("HLT_QuadPFJet98_83_71_15_v",&HLT_QuadPFJet98_83_71_15_v);
// smalltree->Branch("HLT_QuadPFJet103_88_75_15_v",&HLT_QuadPFJet103_88_75_15_v);
// smalltree->Branch("HLT_QuadPFJet105_88_76_15_v",&HLT_QuadPFJet105_88_76_15_v);
// smalltree->Branch("HLT_QuadPFJet111_90_80_15_v",&HLT_QuadPFJet111_90_80_15_v);
// smalltree->Branch("HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v",&HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v);
// // ----------------Trigger IsoMu-------------
// smalltree->Branch("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v",&HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v);
// smalltree->Branch("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v",&HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v);
// smalltree->Branch("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v",&HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v);
// smalltree->Branch("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v",&HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v);
// smalltree->Branch("HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v",&HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v);
// smalltree->Branch("HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v",&HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v);
// smalltree->Branch("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v",&HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v);
// smalltree->Branch("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v",&HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v);
// smalltree->Branch("HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v",&HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v);
// smalltree->Branch("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v",&HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v);
// smalltree->Branch("HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v",&HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v);
// smalltree->Branch("HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v",&HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v);
// smalltree->Branch("HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v",&HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v);
// smalltree->Branch("HLT_IsoMu20_v",&HLT_IsoMu20_v);
smalltree->Branch("HLT_IsoMu24_v",&HLT_IsoMu24_v);
// smalltree->Branch("HLT_IsoMu24_eta2p1_v",&HLT_IsoMu24_eta2p1_v);
smalltree->Branch("HLT_IsoMu27_v",&HLT_IsoMu27_v);
// smalltree->Branch("HLT_IsoMu30_v",&HLT_IsoMu30_v);
// smalltree->Branch("HLT_IsoMu24_TwoProngs35_v",&HLT_IsoMu24_TwoProngs35_v);
// smalltree->Branch("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v",&HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v);
// smalltree->Branch("HLT_IsoMu27_MET90_v",&HLT_IsoMu27_MET90_v);
// // ----------------Trigger PFHT-------------
// smalltree->Branch("HLT_PFHT180_v",&HLT_PFHT180_v);
// smalltree->Branch("HLT_PFHT250_v",&HLT_PFHT250_v);
// smalltree->Branch("HLT_PFHT370_v",&HLT_PFHT370_v);
// smalltree->Branch("HLT_PFHT430_v",&HLT_PFHT430_v);
// smalltree->Branch("HLT_PFHT510_v",&HLT_PFHT510_v);
// smalltree->Branch("HLT_PFHT590_v",&HLT_PFHT590_v);
// smalltree->Branch("HLT_PFHT680_v",&HLT_PFHT680_v);
// smalltree->Branch("HLT_PFHT780_v",&HLT_PFHT780_v);
// smalltree->Branch("HLT_PFHT890_v",&HLT_PFHT890_v);
// smalltree->Branch("HLT_PFHT1050_v",&HLT_PFHT1050_v);
smalltree->Branch("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v",&HLT_PFHT500_PFMET100_PFMHT100_IDTight_v);
// smalltree->Branch("HLT_PFHT500_PFMET110_PFMHT110_IDTight_v",&HLT_PFHT500_PFMET110_PFMHT110_IDTight_v);
smalltree->Branch("HLT_PFHT700_PFMET85_PFMHT85_IDTight_v",&HLT_PFHT700_PFMET85_PFMHT85_IDTight_v);
// smalltree->Branch("HLT_PFHT700_PFMET95_PFMHT95_IDTight_v",&HLT_PFHT700_PFMET95_PFMHT95_IDTight_v);
smalltree->Branch("HLT_PFHT800_PFMET75_PFMHT75_IDTight_v",&HLT_PFHT800_PFMET75_PFMHT75_IDTight_v);
// smalltree->Branch("HLT_PFHT800_PFMET85_PFMHT85_IDTight_v",&HLT_PFHT800_PFMET85_PFMHT85_IDTight_v);
// smalltree->Branch("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v",&HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v);
// smalltree->Branch("HLT_PFHT330PT30_QuadPFJet_75_60_45_40_v",&HLT_PFHT330PT30_QuadPFJet_75_60_45_40_v);
// smalltree->Branch("HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_v",&HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_v);
// smalltree->Branch("HLT_PFHT400_SixPFJet32_v",&HLT_PFHT400_SixPFJet32_v);
// smalltree->Branch("HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_v",&HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_v);
// smalltree->Branch("HLT_PFHT450_SixPFJet36_v",&HLT_PFHT450_SixPFJet36_v);
// smalltree->Branch("HLT_PFHT350_v",&HLT_PFHT350_v);
// smalltree->Branch("HLT_PFHT350MinPFJet15_v",&HLT_PFHT350MinPFJet15_v);

//$$
    //add the variables from my BDT (Paul)
    reader->AddVariable( "mva_track_pt", &pt );
    reader->AddVariable( "mva_track_eta", &eta );
    reader->AddVariable( "mva_track_nchi2", &NChi );
    reader->AddVariable( "mva_track_nhits", &nhits );
    reader->AddVariable( "mva_ntrk10", &ntrk10);
    reader->AddVariable( "mva_drSig", &drSig); /*!*/
    reader->AddVariable( "mva_track_isinjet", &isinjet); /*!*/
    reader->BookMVA( "BDTG", weightFile_ ); // root 6.14/09, care compatiblity of versions for tmva
//$$
}


// FlyingTopAnalyzer::~FlyingTopAnalyzer()
// {
//    // do anything here that needs to be done at destruction time
//    // (e.g. close files, deallocate resources etc.)
// }


bool FlyingTopAnalyzer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
//particle is already the ancestor
  if ( ancestor == particle ) return true;

//otherwise loop on mothers, if any and return true if the ancestor is found
  for (size_t i=0; i < particle->numberOfMothers(); i++)
  {
    if ( isAncestor(ancestor,particle->mother(i)) ) return true;
  }
//if we did not return yet, then particle and ancestor are not relatives
  return false;
}


//
// member functions
//

// ------------ method called for each event  ------------
void FlyingTopAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  clearVariables();

  runNumber   = iEvent.id().run();
  eventNumber = iEvent.id().event();
  lumiBlock   = iEvent.luminosityBlock();

  using namespace edm;
  using namespace reco;
  using namespace pat;

  bool runOnData_ = false;

  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > pruned;
  if ( !runOnData_ ) iEvent.getByToken(prunedGenToken_, pruned);

  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  Handle<edm::View<pat::PackedGenParticle> > packed;
  if ( !runOnData_ ) iEvent.getByToken(packedGenToken_, packed);

  edm::Handle<edm::View<reco::GenJet>> genJets;
  if ( !runOnData_ ) iEvent.getByToken(genJetToken_, genJets);

  edm::Handle<reco::VertexCollection> primaryVertex;
  iEvent.getByToken(vertexToken_, primaryVertex);
  edm::Handle<pat::METCollection> PFMETs;
  iEvent.getByToken(metToken_, PFMETs);
  edm::Handle<pat::JetCollection> jets;
  iEvent.getByToken(jetToken_, jets);
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);
  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  edm::Handle<pat::PackedCandidateCollection> pcs;
  iEvent.getByToken(pcToken_, pcs);
  const pat::PackedCandidateCollection* pc = pcs.product();    

  //LOST
  edm::Handle<pat::PackedCandidateCollection> lostpcs;
  iEvent.getByToken(lostpcToken_, lostpcs);
  const pat::PackedCandidateCollection* lostpc = lostpcs.product(); 

  //V0 particles Collections
  //Kshort
    edm::Handle<reco::VertexCompositePtrCandidateCollection> KshortVertex;
  iEvent.getByToken(K0Token_, KshortVertex);
  //Lambda
    edm::Handle<reco::VertexCompositePtrCandidateCollection> LambdaVertex;
  iEvent.getByToken(LambdaToken_, LambdaVertex);
  //SV
    edm::Handle<reco::VertexCompositePtrCandidateCollection> SVertex;
  iEvent.getByToken(SVToken_, SVertex);

  //Photon Conversion
  edm::Handle<reco::ConversionCollection> PhotonConversion;
  iEvent.getByToken(PhotonToken_,PhotonConversion);

  //Prescales
  // edm::Handle<pat::PackedTriggerPrescales> prescale;
  // iEvent.getByToken(PrescaleToken_,prescale);
 
//--------------------------------------//
//               Trigger                //
//--------------------------------------//
 //HLT trigger we want to keep:
  //__________________________________//
  //HLT_Mu* (includes dilepton channel HLT_MuXX_EleXX)
  //HLT_Ele*
  //HLT_DoubleEle*
  //HLT_DoubleMu*
  //HLT_Dimuon0
  //HLT_PFMET*
  //HLT_HTXX
  //HLT_AK4
  //HLT_PFJet

  //HLT_DoublePFJets
  //HLT_DiPFJet
  //HLT_QuadPFJet
  //HLT_BTagMu
  //__________________________________//

  std::vector<std::string> TriggerCheck;
  std::vector<std::string> VTrigger;
  std::vector<bool> PassTrigger;

  // ######################
  const edm::Handle<edm::TriggerResults> triggerH = iEvent.getHandle(triggerResultsToken_);
  const auto triggerNames = iEvent.triggerNames(*triggerH);
  std::string TName;
  // if (eventNumber==1)
  // {
  //   for (unsigned int i = 0; i < triggerH->size(); i++) 
  //     {
  //       if (TName.substr(0,6) !="HLT_Mu" && TName.substr(0,7) !="HLT_Ele" && TName.substr(0,12) !="HLT_DoubleMu" && TName.substr(0,13) !="HLT_DoubleEle" &&TName.substr(0,11) !="HLT_Dimuon0" &&  TName.substr(0,9) !="HLT_PFMET" &&TName.substr(0,6) !="HLT_HT" && TName.substr(0,7) !="HLT_AK4" && TName.substr(0,9) !="HLT_PFJet")
  //         {
  //           std::cout<<triggerNames.triggerName(i)<<std::endl;
  //         }
  //     }
  // }
  for (unsigned int i = 0; i < triggerH->size(); i++) 
    {
      TName=triggerNames.triggerName(i);
      TriggerCheck.push_back(TName);
//        if (triggerH->accept(i))
//         {
//           tree_passesTrigger.push_back(i);
//           tree_passesTriggerName.push_back(TName);
//                  // std::cout<<" trigger passed : "<<triggerNames.triggerName(i)<<std::endl;
//           // if(TName.substr(0,6) =="HLT_Mu"){tree_Trigger_Muon.push_back(TName);}//+ dilepton channel emu
//           // if(TName.substr(0,7) =="HLT_Ele"){tree_Trigger_Ele.push_back(TName);}
//           // if(TName.substr(0,12) =="HLT_DoubleMu"){tree_Trigger_DoubleMu.push_back(TName);}
//           // if(TName.substr(0,13) =="HLT_DoubleEle"){tree_Trigger_DoubleEle.push_back(TName);}
//           // if(TName.substr(0,11) =="HLT_Dimuon0"){tree_Trigger_Dimuon0.push_back(TName);}
//           // if(TName.substr(0,9) =="HLT_PFMET"){tree_Trigger_PFMET.push_back(TName);}
//           // if(TName.substr(0,6) =="HLT_HT"){tree_Trigger_HT.push_back(TName);}
//           // if(TName.substr(0,7) =="HLT_AK4"){tree_Trigger_AK4.push_back(TName);}
//           // if(TName.substr(0,9) =="HLT_PFJet"){tree_Trigger_PFJet.push_back(TName);}
//           // if(TName.substr(0,16) =="HLT_DoublePFJets"){tree_Trigger_DoublePFJets.push_back(TName);}
//           // if(TName.substr(0,11) =="HLT_DiPFJet"){tree_Trigger_DiPFJet.push_back(TName);}
//           // if(TName.substr(0,13) =="HLT_QuadPFJet"){tree_Trigger_QuadPFJet.push_back(TName);}
//           // if(TName.substr(0,10) =="HLT_BTagMu"){tree_Trigger_BTagMu.push_back(TName);}
//         }

//************************************************************************************************************************//



//   ||||||||||              ||||||||||     |||||||                    |||||||||||||||||||||||||||||||||||||||        
//   ||||||||||              ||||||||||     |||||||                    |||||||||||||||||||||||||||||||||||||||
//   ||||||||||              ||||||||||     |||||||                    |||||||||||||||||||||||||||||||||||||||
//   ||||||||||              ||||||||||     |||||||                    |||||||||||||||||||||||||||||||||||||||
//   ||||||||||              ||||||||||     |||||||                                |||||||||||||||
//   ||||||||||              ||||||||||     |||||||                                |||||||||||||||
//   ||||||||||              ||||||||||     |||||||                                |||||||||||||||
//   ||||||||||||||||||||||||||||||||||     |||||||                                |||||||||||||||
//   ||||||||||||||||||||||||||||||||||     |||||||                                |||||||||||||||
//   ||||||||||||||||||||||||||||||||||     |||||||                                |||||||||||||||
//   ||||||||||              ||||||||||     |||||||                                |||||||||||||||
//   ||||||||||              ||||||||||     |||||||                                |||||||||||||||
//   ||||||||||              ||||||||||     |||||||                                |||||||||||||||
//   ||||||||||              ||||||||||     |||||||                                |||||||||||||||
//   ||||||||||              ||||||||||     |||||||                                |||||||||||||||
//   ||||||||||              ||||||||||     ||||||||||||||||||||||||||||           |||||||||||||||
//   ||||||||||              ||||||||||     ||||||||||||||||||||||||||||           |||||||||||||||
//   ||||||||||              ||||||||||     ||||||||||||||||||||||||||||           |||||||||||||||

//************************************************************************************************************************//

//************************************************************************************************************************//

// // ----------------Trigger Muon + dilepton-------------
//    if (strstr(TName.c_str(),"HLT_Mu27_Ele37_CaloIdL_MW_v") &&  triggerH->accept(i)){HLT_Mu27_Ele37_CaloIdL_MW_v = true;} else if (strstr(TName.c_str(),"HLT_Mu27_Ele37_CaloIdL_MW_v") && !triggerH->accept(i)){HLT_Mu27_Ele37_CaloIdL_MW_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu37_Ele27_CaloIdL_MW_v") &&  triggerH->accept(i)){HLT_Mu37_Ele27_CaloIdL_MW_v = true;} else if (strstr(TName.c_str(),"HLT_Mu37_Ele27_CaloIdL_MW_v") && !triggerH->accept(i)){HLT_Mu37_Ele27_CaloIdL_MW_v = false;};
// /* 55%*/if (strstr(TName.c_str(),"HLT_Mu37_TkMu27_v") &&  triggerH->accept(i)){HLT_Mu37_TkMu27_v = true;} else if (strstr(TName.c_str(),"HLT_Mu37_TkMu27_v") && !triggerH->accept(i)){HLT_Mu37_TkMu27_v = false;};
//   if (strstr(TName.c_str(),"HLT_Mu3_PFJet40_v") &&  triggerH->accept(i)){HLT_Mu3_PFJet40_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3_PFJet40_v") && !triggerH->accept(i)){HLT_Mu3_PFJet40_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu7p5_L2Mu2_Jpsi_v") &&  triggerH->accept(i)){HLT_Mu7p5_L2Mu2_Jpsi_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7p5_L2Mu2_Jpsi_v") && !triggerH->accept(i)){HLT_Mu7p5_L2Mu2_Jpsi_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu7p5_L2Mu2_Upsilon_v") &&  triggerH->accept(i)){HLT_Mu7p5_L2Mu2_Upsilon_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7p5_L2Mu2_Upsilon_v") && !triggerH->accept(i)){HLT_Mu7p5_L2Mu2_Upsilon_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu7p5_Track2_Jpsi_v") &&  triggerH->accept(i)){HLT_Mu7p5_Track2_Jpsi_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7p5_Track2_Jpsi_v") && !triggerH->accept(i)){HLT_Mu7p5_Track2_Jpsi_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu7p5_Track3p5_Jpsi_v") &&  triggerH->accept(i)){HLT_Mu7p5_Track3p5_Jpsi_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7p5_Track3p5_Jpsi_v") && !triggerH->accept(i)){HLT_Mu7p5_Track3p5_Jpsi_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu7p5_Track7_Jpsi_v") &&  triggerH->accept(i)){HLT_Mu7p5_Track7_Jpsi_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7p5_Track7_Jpsi_v") && !triggerH->accept(i)){HLT_Mu7p5_Track7_Jpsi_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu7p5_Track2_Upsilon_v") &&  triggerH->accept(i)){HLT_Mu7p5_Track2_Upsilon_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7p5_Track2_Upsilon_v") && !triggerH->accept(i)){HLT_Mu7p5_Track2_Upsilon_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu7p5_Track3p5_Upsilon_v") &&  triggerH->accept(i)){HLT_Mu7p5_Track3p5_Upsilon_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7p5_Track3p5_Upsilon_v") && !triggerH->accept(i)){HLT_Mu7p5_Track3p5_Upsilon_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu7p5_Track7_Upsilon_v") &&  triggerH->accept(i)){HLT_Mu7p5_Track7_Upsilon_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7p5_Track7_Upsilon_v") && !triggerH->accept(i)){HLT_Mu7p5_Track7_Upsilon_v = false;};
// /* 98.5%*/if (strstr(TName.c_str(),"HLT_Mu3_L1SingleMu5orSingleMu7_v") &&  triggerH->accept(i)){HLT_Mu3_L1SingleMu5orSingleMu7_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3_L1SingleMu5orSingleMu7_v") && !triggerH->accept(i)){HLT_Mu3_L1SingleMu5orSingleMu7_v = false;};
/* 78%*/ if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") && !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v = false;};//std::cout<<"triggerName / prescale : "<<TName<<" / "<<prescale->getPrescaleForInder(i)<<std::endl;
// /* 78%*/if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_v") &&  triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_v") && !triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_v = false;};
// /* 75%*/if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") && !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v = false;};
// /* 75%*/if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v") && !triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_v = false;};
// /* 75%*/if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v") && !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v = false;};
// /* 75%*/if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v") &&  triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v = true;} else if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v") && !triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass8_v = false;};
// /* 76%*/if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") && !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v = false;};
// /* 75%*/if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v") &&  triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v = true;} else if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v") && !triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_Mu9_TrkIsoVVL_DZ_Mass3p8_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu25_TkMu0_Onia_v") &&  triggerH->accept(i)){HLT_Mu25_TkMu0_Onia_v = true;} else if (strstr(TName.c_str(),"HLT_Mu25_TkMu0_Onia_v") && !triggerH->accept(i)){HLT_Mu25_TkMu0_Onia_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu30_TkMu0_Psi_v") &&  triggerH->accept(i)){HLT_Mu30_TkMu0_Psi_v = true;} else if (strstr(TName.c_str(),"HLT_Mu30_TkMu0_Psi_v") && !triggerH->accept(i)){HLT_Mu30_TkMu0_Psi_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu30_TkMu0_Upsilon_v") &&  triggerH->accept(i)){HLT_Mu30_TkMu0_Upsilon_v = true;} else if (strstr(TName.c_str(),"HLT_Mu30_TkMu0_Upsilon_v") && !triggerH->accept(i)){HLT_Mu30_TkMu0_Upsilon_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu20_TkMu0_Phi_v") &&  triggerH->accept(i)){HLT_Mu20_TkMu0_Phi_v = true;} else if (strstr(TName.c_str(),"HLT_Mu20_TkMu0_Phi_v") && !triggerH->accept(i)){HLT_Mu20_TkMu0_Phi_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu25_TkMu0_Phi_v") &&  triggerH->accept(i)){HLT_Mu25_TkMu0_Phi_v = true;} else if (strstr(TName.c_str(),"HLT_Mu25_TkMu0_Phi_v") && !triggerH->accept(i)){HLT_Mu25_TkMu0_Phi_v = false;};
// /* 97%*/if (strstr(TName.c_str(),"HLT_Mu12_v") &&  triggerH->accept(i)){HLT_Mu12_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_v") && !triggerH->accept(i)){HLT_Mu12_v = false;};
// /* 97%*/if (strstr(TName.c_str(),"HLT_Mu15_v") &&  triggerH->accept(i)){HLT_Mu15_v = true;} else if (strstr(TName.c_str(),"HLT_Mu15_v") && !triggerH->accept(i)){HLT_Mu15_v = false;};
// /* 95%*/if (strstr(TName.c_str(),"HLT_Mu20_v") &&  triggerH->accept(i)){HLT_Mu20_v = true;} else if (strstr(TName.c_str(),"HLT_Mu20_v") && !triggerH->accept(i)){HLT_Mu20_v = false;};
// /* 91%*/if (strstr(TName.c_str(),"HLT_Mu27_v") &&  triggerH->accept(i)){HLT_Mu27_v = true;} else if (strstr(TName.c_str(),"HLT_Mu27_v") && !triggerH->accept(i)){HLT_Mu27_v = false;};
// /* 74%*/if (strstr(TName.c_str(),"HLT_Mu50_v") &&  triggerH->accept(i)){HLT_Mu50_v = true;} else if (strstr(TName.c_str(),"HLT_Mu50_v") && !triggerH->accept(i)){HLT_Mu50_v = false;};
// /* 67*/if (strstr(TName.c_str(),"HLT_Mu55_v") &&  triggerH->accept(i)){HLT_Mu55_v = true;} else if (strstr(TName.c_str(),"HLT_Mu55_v") && !triggerH->accept(i)){HLT_Mu55_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_v") &&  triggerH->accept(i)){HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_v") && !triggerH->accept(i)){HLT_Mu12_DoublePFJets40_CaloBTagDeepCSV_p71_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_v") &&  triggerH->accept(i)){HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_v") && !triggerH->accept(i)){HLT_Mu12_DoublePFJets100_CaloBTagDeepCSV_p71_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_v") &&  triggerH->accept(i)){HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_v") && !triggerH->accept(i)){HLT_Mu12_DoublePFJets200_CaloBTagDeepCSV_p71_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_v") &&  triggerH->accept(i)){HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_v") && !triggerH->accept(i)){HLT_Mu12_DoublePFJets350_CaloBTagDeepCSV_p71_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") &&  triggerH->accept(i)){HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") && !triggerH->accept(i)){HLT_Mu12_DoublePFJets40MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") &&  triggerH->accept(i)){HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") && !triggerH->accept(i)){HLT_Mu12_DoublePFJets54MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") &&  triggerH->accept(i)){HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") && !triggerH->accept(i)){HLT_Mu12_DoublePFJets62MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = false;};
// /* 98%*/if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v") &&  triggerH->accept(i)){HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v") && !triggerH->accept(i)){HLT_Mu8_DiEle12_CaloIdL_TrackIdL_DZ_v = false;};
// /* 1Null Efficacity*/if (strstr(TName.c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v") &&  triggerH->accept(i)){HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v") && !triggerH->accept(i)){HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v") &&  triggerH->accept(i)){HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v") && !triggerH->accept(i)){HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v") &&  triggerH->accept(i)){HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v") && !triggerH->accept(i)){HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_v = false;};
/* 1%*/if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_PFDiJet30_PFBtagDeepCSV_1p5_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_CaloDiJet30_CaloBtagDeepCSV_1p5_v = false;};
/* 1%*/if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = false;};
// /* 97%*/if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_v") && !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_v = false;};
// /* 96*/if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_v") &&  triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu19_TrkIsoVVL_v") && !triggerH->accept(i)){HLT_Mu19_TrkIsoVVL_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") && !triggerH->accept(i)){HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = false;};
/* 1%*/if (strstr(TName.c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") &&  triggerH->accept(i)){HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") && !triggerH->accept(i)){HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") &&  triggerH->accept(i)){HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") && !triggerH->accept(i)){HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = false;};
/* 1%*/if (strstr(TName.c_str(),"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") && !triggerH->accept(i)){HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = false;};
// /* 10%*/if (strstr(TName.c_str(),"HLT_Mu12_DoublePhoton20_v") &&  triggerH->accept(i)){HLT_Mu12_DoublePhoton20_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_DoublePhoton20_v") && !triggerH->accept(i)){HLT_Mu12_DoublePhoton20_v = false;};
// /* 11%*/if (strstr(TName.c_str(),"HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v") &&  triggerH->accept(i)){HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v") && !triggerH->accept(i)){HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL_v = false;};
// /* 11%*/if (strstr(TName.c_str(),"HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_v") &&  triggerH->accept(i)){HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_v") && !triggerH->accept(i)){HLT_Mu48NoFiltersNoVtx_Photon48_CaloIdL_v = false;};
// /* 2%*/if (strstr(TName.c_str(),"HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_v") &&  triggerH->accept(i)){HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_v") && !triggerH->accept(i)){HLT_Mu38NoFiltersNoVtxDisplaced_Photon38_CaloIdL_v = false;};
// /* 2%*/if (strstr(TName.c_str(),"HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_v") &&  triggerH->accept(i)){HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_v") && !triggerH->accept(i)){HLT_Mu43NoFiltersNoVtxDisplaced_Photon43_CaloIdL_v = false;};
// /* 2%*/if (strstr(TName.c_str(),"HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v") &&  triggerH->accept(i)){HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v = true;} else if (strstr(TName.c_str(),"HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v") && !triggerH->accept(i)){HLT_Mu4_TrkIsoVVL_DiPFJet90_40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v = false;};
// /* 2%*/if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60_v = false;};
// /* 2%*/if (strstr(TName.c_str(),"HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v") &&  triggerH->accept(i)){HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v = true;} else if (strstr(TName.c_str(),"HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v") && !triggerH->accept(i)){HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60_v = false;};
// /* 15%*/if (strstr(TName.c_str(),"HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v") &&  triggerH->accept(i)){HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v = true;} else if (strstr(TName.c_str(),"HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v") && !triggerH->accept(i)){HLT_Mu15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v = false;};
// /* 30%*/if (strstr(TName.c_str(),"HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v") &&  triggerH->accept(i)){HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v = true;} else if (strstr(TName.c_str(),"HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v") && !triggerH->accept(i)){HLT_Mu15_IsoVVVL_PFHT450_PFMET50_v = false;};
// /* 52%*/if (strstr(TName.c_str(),"HLT_Mu15_IsoVVVL_PFHT450_v") &&  triggerH->accept(i)){HLT_Mu15_IsoVVVL_PFHT450_v = true;} else if (strstr(TName.c_str(),"HLT_Mu15_IsoVVVL_PFHT450_v") && !triggerH->accept(i)){HLT_Mu15_IsoVVVL_PFHT450_v = false;};
// /* 45%*/if (strstr(TName.c_str(),"HLT_Mu50_IsoVVVL_PFHT450_v") &&  triggerH->accept(i)){HLT_Mu50_IsoVVVL_PFHT450_v = true;} else if (strstr(TName.c_str(),"HLT_Mu50_IsoVVVL_PFHT450_v") && !triggerH->accept(i)){HLT_Mu50_IsoVVVL_PFHT450_v = false;};
// /* 33%*/if (strstr(TName.c_str(),"HLT_Mu15_IsoVVVL_PFHT600_v") &&  triggerH->accept(i)){HLT_Mu15_IsoVVVL_PFHT600_v = true;} else if (strstr(TName.c_str(),"HLT_Mu15_IsoVVVL_PFHT600_v") && !triggerH->accept(i)){HLT_Mu15_IsoVVVL_PFHT600_v = false;};
// /* 19%*/if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_v") &&  triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_v") && !triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMET70_PFMHT70_IDTight_v = false;};
// /* 15%*/if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_v") &&  triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_v") && !triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMET80_PFMHT80_IDTight_v = false;};
// /* 71%*/if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_v") &&  triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_v") && !triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMET90_PFMHT90_IDTight_v = false;};
// /* 10%*/if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_v") &&  triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_v") && !triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMET100_PFMHT100_IDTight_v = false;};
// /* 31%*/if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_v") &&  triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_v") && !triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu70_PFMHTNoMu70_IDTight_v = false;};
// /* 27%*/if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_v") &&  triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_v") && !triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu80_PFMHTNoMu80_IDTight_v = false;};
// /* 72%*/if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_v") &&  triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_v") && !triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu90_PFMHTNoMu90_IDTight_v = false;};
// /* 19%*/if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_v") &&  triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_v") && !triggerH->accept(i)){HLT_Mu3er1p5_PFJet100er2p5_PFMETNoMu100_PFMHTNoMu100_IDTight_v = false;};
// /* 98%*/if (strstr(TName.c_str(),"HLT_Mu8_v") &&  triggerH->accept(i)){HLT_Mu8_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_v") && !triggerH->accept(i)){HLT_Mu8_v = false;};
// /* 97%*/if (strstr(TName.c_str(),"HLT_Mu17_v") &&  triggerH->accept(i)){HLT_Mu17_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_v") && !triggerH->accept(i)){HLT_Mu17_v = false;};
// /* 96%*/if (strstr(TName.c_str(),"HLT_Mu19_v") &&  triggerH->accept(i)){HLT_Mu19_v = true;} else if (strstr(TName.c_str(),"HLT_Mu19_v") && !triggerH->accept(i)){HLT_Mu19_v = false;};
// /* 15%*/if (strstr(TName.c_str(),"HLT_Mu17_Photon30_IsoCaloId_v") &&  triggerH->accept(i)){HLT_Mu17_Photon30_IsoCaloId_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_Photon30_IsoCaloId_v") && !triggerH->accept(i)){HLT_Mu17_Photon30_IsoCaloId_v = false;};
// /* 10%*/if (strstr(TName.c_str(),"HLT_Mu18_Mu9_SameSign_v") &&  triggerH->accept(i)){HLT_Mu18_Mu9_SameSign_v = true;} else if (strstr(TName.c_str(),"HLT_Mu18_Mu9_SameSign_v") && !triggerH->accept(i)){HLT_Mu18_Mu9_SameSign_v = false;};
// /* 9%*/if (strstr(TName.c_str(),"HLT_Mu18_Mu9_SameSign_DZ_v") &&  triggerH->accept(i)){HLT_Mu18_Mu9_SameSign_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu18_Mu9_SameSign_DZ_v") && !triggerH->accept(i)){HLT_Mu18_Mu9_SameSign_DZ_v = false;};
// /* 78%*/if (strstr(TName.c_str(),"HLT_Mu18_Mu9_v") &&  triggerH->accept(i)){HLT_Mu18_Mu9_v = true;} else if (strstr(TName.c_str(),"HLT_Mu18_Mu9_v") && !triggerH->accept(i)){HLT_Mu18_Mu9_v = false;};
// /* 76%*/if (strstr(TName.c_str(),"HLT_Mu18_Mu9_DZ_v") &&  triggerH->accept(i)){HLT_Mu18_Mu9_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu18_Mu9_DZ_v") && !triggerH->accept(i)){HLT_Mu18_Mu9_DZ_v = false;};
// /* 10%*/if (strstr(TName.c_str(),"HLT_Mu20_Mu10_SameSign_v") &&  triggerH->accept(i)){HLT_Mu20_Mu10_SameSign_v = true;} else if (strstr(TName.c_str(),"HLT_Mu20_Mu10_SameSign_v") && !triggerH->accept(i)){HLT_Mu20_Mu10_SameSign_v = false;};
// /* 9%*/if (strstr(TName.c_str(),"HLT_Mu20_Mu10_SameSign_DZ_v") &&  triggerH->accept(i)){HLT_Mu20_Mu10_SameSign_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu20_Mu10_SameSign_DZ_v") && !triggerH->accept(i)){HLT_Mu20_Mu10_SameSign_DZ_v = false;};
// /* 77%*/if (strstr(TName.c_str(),"HLT_Mu20_Mu10_v") &&  triggerH->accept(i)){HLT_Mu20_Mu10_v = true;} else if (strstr(TName.c_str(),"HLT_Mu20_Mu10_v") && !triggerH->accept(i)){HLT_Mu20_Mu10_v = false;};
// /* 75%*/if (strstr(TName.c_str(),"HLT_Mu20_Mu10_DZ_v") &&  triggerH->accept(i)){HLT_Mu20_Mu10_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu20_Mu10_DZ_v") && !triggerH->accept(i)){HLT_Mu20_Mu10_DZ_v = false;};
// /* 10%*/if (strstr(TName.c_str(),"HLT_Mu23_Mu12_SameSign_v") &&  triggerH->accept(i)){HLT_Mu23_Mu12_SameSign_v = true;} else if (strstr(TName.c_str(),"HLT_Mu23_Mu12_SameSign_v") && !triggerH->accept(i)){HLT_Mu23_Mu12_SameSign_v = false;};
// /* 9%*/if (strstr(TName.c_str(),"HLT_Mu23_Mu12_SameSign_DZ_v") &&  triggerH->accept(i)){HLT_Mu23_Mu12_SameSign_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu23_Mu12_SameSign_DZ_v") && !triggerH->accept(i)){HLT_Mu23_Mu12_SameSign_DZ_v = false;};
// /* 77%*/if (strstr(TName.c_str(),"HLT_Mu23_Mu12_v") &&  triggerH->accept(i)){HLT_Mu23_Mu12_v = true;} else if (strstr(TName.c_str(),"HLT_Mu23_Mu12_v") && !triggerH->accept(i)){HLT_Mu23_Mu12_v = false;};
// /* 75%*/if (strstr(TName.c_str(),"HLT_Mu23_Mu12_DZ_v") &&  triggerH->accept(i)){HLT_Mu23_Mu12_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu23_Mu12_DZ_v") && !triggerH->accept(i)){HLT_Mu23_Mu12_DZ_v = false;};
// /* 9%*/if (strstr(TName.c_str(),"HLT_Mu12_IP6_part0_v") &&  triggerH->accept(i)){HLT_Mu12_IP6_part0_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_IP6_part0_v") && !triggerH->accept(i)){HLT_Mu12_IP6_part0_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu12_IP6_part1_v") &&  triggerH->accept(i)){HLT_Mu12_IP6_part1_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_IP6_part1_v") && !triggerH->accept(i)){HLT_Mu12_IP6_part1_v = false;};
// /* 7.4%*/if (strstr(TName.c_str(),"HLT_Mu12_IP6_part2_v") &&  triggerH->accept(i)){HLT_Mu12_IP6_part2_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_IP6_part2_v") && !triggerH->accept(i)){HLT_Mu12_IP6_part2_v = false;};
// /* 7.2%*/if (strstr(TName.c_str(),"HLT_Mu12_IP6_part3_v") &&  triggerH->accept(i)){HLT_Mu12_IP6_part3_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_IP6_part3_v") && !triggerH->accept(i)){HLT_Mu12_IP6_part3_v = false;};
// /* 7.5*/if (strstr(TName.c_str(),"HLT_Mu12_IP6_part4_v") &&  triggerH->accept(i)){HLT_Mu12_IP6_part4_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_IP6_part4_v") && !triggerH->accept(i)){HLT_Mu12_IP6_part4_v = false;};
// /* 8.6*/if (strstr(TName.c_str(),"HLT_Mu9_IP5_part0_v") &&  triggerH->accept(i)){HLT_Mu9_IP5_part0_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP5_part0_v") && !triggerH->accept(i)){HLT_Mu9_IP5_part0_v = false;};
// /* 8.6*/if (strstr(TName.c_str(),"HLT_Mu9_IP5_part1_v") &&  triggerH->accept(i)){HLT_Mu9_IP5_part1_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP5_part1_v") && !triggerH->accept(i)){HLT_Mu9_IP5_part1_v = false;};
// /* 8.6*/if (strstr(TName.c_str(),"HLT_Mu9_IP5_part2_v") &&  triggerH->accept(i)){HLT_Mu9_IP5_part2_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP5_part2_v") && !triggerH->accept(i)){HLT_Mu9_IP5_part2_v = false;};
// /* 8.6*/if (strstr(TName.c_str(),"HLT_Mu9_IP5_part3_v") &&  triggerH->accept(i)){HLT_Mu9_IP5_part3_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP5_part3_v") && !triggerH->accept(i)){HLT_Mu9_IP5_part3_v = false;};
// /* 8.6*/if (strstr(TName.c_str(),"HLT_Mu9_IP5_part4_v") &&  triggerH->accept(i)){HLT_Mu9_IP5_part4_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP5_part4_v") && !triggerH->accept(i)){HLT_Mu9_IP5_part4_v = false;};
// /* 10*/if (strstr(TName.c_str(),"HLT_Mu7_IP4_part0_v") &&  triggerH->accept(i)){HLT_Mu7_IP4_part0_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7_IP4_part0_v") && !triggerH->accept(i)){HLT_Mu7_IP4_part0_v = false;};
// /* 10*/if (strstr(TName.c_str(),"HLT_Mu7_IP4_part1_v") &&  triggerH->accept(i)){HLT_Mu7_IP4_part1_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7_IP4_part1_v") && !triggerH->accept(i)){HLT_Mu7_IP4_part1_v = false;};
// /* 10*/if (strstr(TName.c_str(),"HLT_Mu7_IP4_part2_v") &&  triggerH->accept(i)){HLT_Mu7_IP4_part2_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7_IP4_part2_v") && !triggerH->accept(i)){HLT_Mu7_IP4_part2_v = false;};
// /* 10*/if (strstr(TName.c_str(),"HLT_Mu7_IP4_part3_v") &&  triggerH->accept(i)){HLT_Mu7_IP4_part3_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7_IP4_part3_v") && !triggerH->accept(i)){HLT_Mu7_IP4_part3_v = false;};
// /* 10*/if (strstr(TName.c_str(),"HLT_Mu7_IP4_part4_v") &&  triggerH->accept(i)){HLT_Mu7_IP4_part4_v = true;} else if (strstr(TName.c_str(),"HLT_Mu7_IP4_part4_v") && !triggerH->accept(i)){HLT_Mu7_IP4_part4_v = false;};
// /* 9%*/if (strstr(TName.c_str(),"HLT_Mu9_IP4_part0_v") &&  triggerH->accept(i)){HLT_Mu9_IP4_part0_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP4_part0_v") && !triggerH->accept(i)){HLT_Mu9_IP4_part0_v = false;};
// /* 9%*/if (strstr(TName.c_str(),"HLT_Mu9_IP4_part1_v") &&  triggerH->accept(i)){HLT_Mu9_IP4_part1_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP4_part1_v") && !triggerH->accept(i)){HLT_Mu9_IP4_part1_v = false;};
// /* 9%*/if (strstr(TName.c_str(),"HLT_Mu9_IP4_part2_v") &&  triggerH->accept(i)){HLT_Mu9_IP4_part2_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP4_part2_v") && !triggerH->accept(i)){HLT_Mu9_IP4_part2_v = false;};
// /* 9%*/if (strstr(TName.c_str(),"HLT_Mu9_IP4_part3_v") &&  triggerH->accept(i)){HLT_Mu9_IP4_part3_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP4_part3_v") && !triggerH->accept(i)){HLT_Mu9_IP4_part3_v = false;};
// /* 9%*/if (strstr(TName.c_str(),"HLT_Mu9_IP4_part4_v") &&  triggerH->accept(i)){HLT_Mu9_IP4_part4_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP4_part4_v") && !triggerH->accept(i)){HLT_Mu9_IP4_part4_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP5_part0_v") &&  triggerH->accept(i)){HLT_Mu8_IP5_part0_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP5_part0_v") && !triggerH->accept(i)){HLT_Mu8_IP5_part0_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP5_part1_v") &&  triggerH->accept(i)){HLT_Mu8_IP5_part1_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP5_part1_v") && !triggerH->accept(i)){HLT_Mu8_IP5_part1_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP5_part2_v") &&  triggerH->accept(i)){HLT_Mu8_IP5_part2_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP5_part2_v") && !triggerH->accept(i)){HLT_Mu8_IP5_part2_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP5_part3_v") &&  triggerH->accept(i)){HLT_Mu8_IP5_part3_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP5_part3_v") && !triggerH->accept(i)){HLT_Mu8_IP5_part3_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP5_part4_v") &&  triggerH->accept(i)){HLT_Mu8_IP5_part4_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP5_part4_v") && !triggerH->accept(i)){HLT_Mu8_IP5_part4_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP6_part0_v") &&  triggerH->accept(i)){HLT_Mu8_IP6_part0_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP6_part0_v") && !triggerH->accept(i)){HLT_Mu8_IP6_part0_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP6_part1_v") &&  triggerH->accept(i)){HLT_Mu8_IP6_part1_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP6_part1_v") && !triggerH->accept(i)){HLT_Mu8_IP6_part1_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP6_part2_v") &&  triggerH->accept(i)){HLT_Mu8_IP6_part2_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP6_part2_v") && !triggerH->accept(i)){HLT_Mu8_IP6_part2_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP6_part3_v") &&  triggerH->accept(i)){HLT_Mu8_IP6_part3_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP6_part3_v") && !triggerH->accept(i)){HLT_Mu8_IP6_part3_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu8_IP6_part4_v") &&  triggerH->accept(i)){HLT_Mu8_IP6_part4_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP6_part4_v") && !triggerH->accept(i)){HLT_Mu8_IP6_part4_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu9_IP6_part0_v") &&  triggerH->accept(i)){HLT_Mu9_IP6_part0_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP6_part0_v") && !triggerH->accept(i)){HLT_Mu9_IP6_part0_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu9_IP6_part1_v") &&  triggerH->accept(i)){HLT_Mu9_IP6_part1_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP6_part1_v") && !triggerH->accept(i)){HLT_Mu9_IP6_part1_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu9_IP6_part2_v") &&  triggerH->accept(i)){HLT_Mu9_IP6_part2_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP6_part2_v") && !triggerH->accept(i)){HLT_Mu9_IP6_part2_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu9_IP6_part3_v") &&  triggerH->accept(i)){HLT_Mu9_IP6_part3_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP6_part3_v") && !triggerH->accept(i)){HLT_Mu9_IP6_part3_v = false;};
// /* 8%*/if (strstr(TName.c_str(),"HLT_Mu9_IP6_part4_v") &&  triggerH->accept(i)){HLT_Mu9_IP6_part4_v = true;} else if (strstr(TName.c_str(),"HLT_Mu9_IP6_part4_v") && !triggerH->accept(i)){HLT_Mu9_IP6_part4_v = false;};
// /* 13%*/if (strstr(TName.c_str(),"HLT_Mu8_IP3_part0_v") &&  triggerH->accept(i)){HLT_Mu8_IP3_part0_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP3_part0_v") && !triggerH->accept(i)){HLT_Mu8_IP3_part0_v = false;};
// /* 13%*/if (strstr(TName.c_str(),"HLT_Mu8_IP3_part1_v") &&  triggerH->accept(i)){HLT_Mu8_IP3_part1_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP3_part1_v") && !triggerH->accept(i)){HLT_Mu8_IP3_part1_v = false;};
// /* 13%*/if (strstr(TName.c_str(),"HLT_Mu8_IP3_part2_v") &&  triggerH->accept(i)){HLT_Mu8_IP3_part2_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP3_part2_v") && !triggerH->accept(i)){HLT_Mu8_IP3_part2_v = false;};
// /* 13%*/if (strstr(TName.c_str(),"HLT_Mu8_IP3_part3_v") &&  triggerH->accept(i)){HLT_Mu8_IP3_part3_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP3_part3_v") && !triggerH->accept(i)){HLT_Mu8_IP3_part3_v = false;};
// /* 13%*/if (strstr(TName.c_str(),"HLT_Mu8_IP3_part4_v") &&  triggerH->accept(i)){HLT_Mu8_IP3_part4_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_IP3_part4_v") && !triggerH->accept(i)){HLT_Mu8_IP3_part4_v = false;};
// // ----------------Trigger Electron-------------
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele27_Ele37_CaloIdL_MW_v") &&  triggerH->accept(i)){HLT_Ele27_Ele37_CaloIdL_MW_v = true;} else if (strstr(TName.c_str(),"HLT_Ele27_Ele37_CaloIdL_MW_v") && !triggerH->accept(i)){HLT_Ele27_Ele37_CaloIdL_MW_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele20_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele20_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele20_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele20_WPTight_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele15_WPLoose_Gsf_v") &&  triggerH->accept(i)){HLT_Ele15_WPLoose_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele15_WPLoose_Gsf_v") && !triggerH->accept(i)){HLT_Ele15_WPLoose_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele17_WPLoose_Gsf_v") &&  triggerH->accept(i)){HLT_Ele17_WPLoose_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele17_WPLoose_Gsf_v") && !triggerH->accept(i)){HLT_Ele17_WPLoose_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele20_WPLoose_Gsf_v") &&  triggerH->accept(i)){HLT_Ele20_WPLoose_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele20_WPLoose_Gsf_v") && !triggerH->accept(i)){HLT_Ele20_WPLoose_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele20_eta2p1_WPLoose_Gsf_v") &&  triggerH->accept(i)){HLT_Ele20_eta2p1_WPLoose_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele20_eta2p1_WPLoose_Gsf_v") && !triggerH->accept(i)){HLT_Ele20_eta2p1_WPLoose_Gsf_v = false;};
/* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele27_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele27_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele27_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele27_WPTight_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele28_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele28_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele28_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele28_WPTight_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele30_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele30_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele30_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele30_WPTight_Gsf_v = false;};
/* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele32_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele32_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele32_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele32_WPTight_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele35_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele35_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele35_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele35_WPTight_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele35_WPTight_Gsf_L1EGMT_v") &&  triggerH->accept(i)){HLT_Ele35_WPTight_Gsf_L1EGMT_v = true;} else if (strstr(TName.c_str(),"HLT_Ele35_WPTight_Gsf_L1EGMT_v") && !triggerH->accept(i)){HLT_Ele35_WPTight_Gsf_L1EGMT_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele38_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele38_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele38_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele38_WPTight_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele40_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele40_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele40_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele40_WPTight_Gsf_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele32_WPTight_Gsf_L1DoubleEG_v") &&  triggerH->accept(i)){HLT_Ele32_WPTight_Gsf_L1DoubleEG_v = true;} else if (strstr(TName.c_str(),"HLT_Ele32_WPTight_Gsf_L1DoubleEG_v") && !triggerH->accept(i)){HLT_Ele32_WPTight_Gsf_L1DoubleEG_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v") &&  triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v") && !triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_v") &&  triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_v") && !triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_CrossL1_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_v") &&  triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_v") && !triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_CrossL1_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v") &&  triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v") && !triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v") &&  triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v") && !triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_MediumChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v") &&  triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v") && !triggerH->accept(i)){HLT_Ele24_eta2p1_WPTight_Gsf_TightChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v") &&  triggerH->accept(i)){HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v = true;} else if (strstr(TName.c_str(),"HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v") && !triggerH->accept(i)){HLT_Ele15_Ele8_CaloIdL_TrackIdL_IsoVL_v = false;};
/* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") &&  triggerH->accept(i)){HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") && !triggerH->accept(i)){HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = false;};
/* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v") &&  triggerH->accept(i)){HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v = true;} else if (strstr(TName.c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v") && !triggerH->accept(i)){HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v") &&  triggerH->accept(i)){HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v = true;} else if (strstr(TName.c_str(),"HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v") && !triggerH->accept(i)){HLT_Ele30_eta2p1_WPTight_Gsf_CentralPFJet35_EleCleaned_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v") &&  triggerH->accept(i)){HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v = true;} else if (strstr(TName.c_str(),"HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v") && !triggerH->accept(i)){HLT_Ele28_eta2p1_WPTight_Gsf_HT150_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele28_HighEta_SC20_Mass55_v") &&  triggerH->accept(i)){HLT_Ele28_HighEta_SC20_Mass55_v = true;} else if (strstr(TName.c_str(),"HLT_Ele28_HighEta_SC20_Mass55_v") && !triggerH->accept(i)){HLT_Ele28_HighEta_SC20_Mass55_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v") &&  triggerH->accept(i)){HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v = true;} else if (strstr(TName.c_str(),"HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v") && !triggerH->accept(i)){HLT_Ele15_IsoVVVL_PFHT450_CaloBTagDeepCSV_4p5_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v") &&  triggerH->accept(i)){HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v = true;} else if (strstr(TName.c_str(),"HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v") && !triggerH->accept(i)){HLT_Ele15_IsoVVVL_PFHT450_PFMET50_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele15_IsoVVVL_PFHT450_v") &&  triggerH->accept(i)){HLT_Ele15_IsoVVVL_PFHT450_v = true;} else if (strstr(TName.c_str(),"HLT_Ele15_IsoVVVL_PFHT450_v") && !triggerH->accept(i)){HLT_Ele15_IsoVVVL_PFHT450_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele50_IsoVVVL_PFHT450_v") &&  triggerH->accept(i)){HLT_Ele50_IsoVVVL_PFHT450_v = true;} else if (strstr(TName.c_str(),"HLT_Ele50_IsoVVVL_PFHT450_v") && !triggerH->accept(i)){HLT_Ele50_IsoVVVL_PFHT450_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele15_IsoVVVL_PFHT600_v") &&  triggerH->accept(i)){HLT_Ele15_IsoVVVL_PFHT600_v = true;} else if (strstr(TName.c_str(),"HLT_Ele15_IsoVVVL_PFHT600_v") && !triggerH->accept(i)){HLT_Ele15_IsoVVVL_PFHT600_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") &&  triggerH->accept(i)){HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v = true;} else if (strstr(TName.c_str(),"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && !triggerH->accept(i)){HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") &&  triggerH->accept(i)){HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v = true;} else if (strstr(TName.c_str(),"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && !triggerH->accept(i)){HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v") &&  triggerH->accept(i)){HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v = true;} else if (strstr(TName.c_str(),"HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && !triggerH->accept(i)){HLT_Ele15_CaloIdL_TrackIdL_IsoVL_PFJet30_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") &&  triggerH->accept(i)){HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v = true;} else if (strstr(TName.c_str(),"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v") && !triggerH->accept(i)){HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v") &&  triggerH->accept(i)){HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v = true;} else if (strstr(TName.c_str(),"HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v") && !triggerH->accept(i)){HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v") &&  triggerH->accept(i)){HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v = true;} else if (strstr(TName.c_str(),"HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v") && !triggerH->accept(i)){HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v") &&  triggerH->accept(i)){HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v = true;} else if (strstr(TName.c_str(),"HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v") && !triggerH->accept(i)){HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v") &&  triggerH->accept(i)){HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v = true;} else if (strstr(TName.c_str(),"HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v") && !triggerH->accept(i)){HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele115_CaloIdVT_GsfTrkIdT_v") &&  triggerH->accept(i)){HLT_Ele115_CaloIdVT_GsfTrkIdT_v = true;} else if (strstr(TName.c_str(),"HLT_Ele115_CaloIdVT_GsfTrkIdT_v") && !triggerH->accept(i)){HLT_Ele115_CaloIdVT_GsfTrkIdT_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele135_CaloIdVT_GsfTrkIdT_v") &&  triggerH->accept(i)){HLT_Ele135_CaloIdVT_GsfTrkIdT_v = true;} else if (strstr(TName.c_str(),"HLT_Ele135_CaloIdVT_GsfTrkIdT_v") && !triggerH->accept(i)){HLT_Ele135_CaloIdVT_GsfTrkIdT_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele145_CaloIdVT_GsfTrkIdT_v") &&  triggerH->accept(i)){HLT_Ele145_CaloIdVT_GsfTrkIdT_v = true;} else if (strstr(TName.c_str(),"HLT_Ele145_CaloIdVT_GsfTrkIdT_v") && !triggerH->accept(i)){HLT_Ele145_CaloIdVT_GsfTrkIdT_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele200_CaloIdVT_GsfTrkIdT_v") &&  triggerH->accept(i)){HLT_Ele200_CaloIdVT_GsfTrkIdT_v = true;} else if (strstr(TName.c_str(),"HLT_Ele200_CaloIdVT_GsfTrkIdT_v") && !triggerH->accept(i)){HLT_Ele200_CaloIdVT_GsfTrkIdT_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele250_CaloIdVT_GsfTrkIdT_v") &&  triggerH->accept(i)){HLT_Ele250_CaloIdVT_GsfTrkIdT_v = true;} else if (strstr(TName.c_str(),"HLT_Ele250_CaloIdVT_GsfTrkIdT_v") && !triggerH->accept(i)){HLT_Ele250_CaloIdVT_GsfTrkIdT_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele300_CaloIdVT_GsfTrkIdT_v") &&  triggerH->accept(i)){HLT_Ele300_CaloIdVT_GsfTrkIdT_v = true;} else if (strstr(TName.c_str(),"HLT_Ele300_CaloIdVT_GsfTrkIdT_v") && !triggerH->accept(i)){HLT_Ele300_CaloIdVT_GsfTrkIdT_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v") &&  triggerH->accept(i)){HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v = true;} else if (strstr(TName.c_str(),"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v") && !triggerH->accept(i)){HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v = false;};
// // ----------------Trigger DoubleMu-------------
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v") &&  triggerH->accept(i)){HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v") && !triggerH->accept(i)){HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v") &&  triggerH->accept(i)){HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v") && !triggerH->accept(i)){HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu4_3_Bs_v") &&  triggerH->accept(i)){HLT_DoubleMu4_3_Bs_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu4_3_Bs_v") && !triggerH->accept(i)){HLT_DoubleMu4_3_Bs_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu4_3_Jpsi_v") &&  triggerH->accept(i)){HLT_DoubleMu4_3_Jpsi_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu4_3_Jpsi_v") && !triggerH->accept(i)){HLT_DoubleMu4_3_Jpsi_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu4_JpsiTrk_Displaced_v") &&  triggerH->accept(i)){HLT_DoubleMu4_JpsiTrk_Displaced_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu4_JpsiTrk_Displaced_v") && !triggerH->accept(i)){HLT_DoubleMu4_JpsiTrk_Displaced_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v") &&  triggerH->accept(i)){HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v") && !triggerH->accept(i)){HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu3_Trk_Tau3mu_v") &&  triggerH->accept(i)){HLT_DoubleMu3_Trk_Tau3mu_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu3_Trk_Tau3mu_v") && !triggerH->accept(i)){HLT_DoubleMu3_Trk_Tau3mu_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu3_TkMu_DsTau3Mu_v") &&  triggerH->accept(i)){HLT_DoubleMu3_TkMu_DsTau3Mu_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu3_TkMu_DsTau3Mu_v") && !triggerH->accept(i)){HLT_DoubleMu3_TkMu_DsTau3Mu_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu4_PsiPrimeTrk_Displaced_v") &&  triggerH->accept(i)){HLT_DoubleMu4_PsiPrimeTrk_Displaced_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu4_PsiPrimeTrk_Displaced_v") && !triggerH->accept(i)){HLT_DoubleMu4_PsiPrimeTrk_Displaced_v = false;};
// /* 59%*/if (strstr(TName.c_str(),"HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v") &&  triggerH->accept(i)){HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v") && !triggerH->accept(i)){HLT_DoubleMu4_Mass3p8_DZ_PFHT350_v = false;};
// /* 7.1%*/if (strstr(TName.c_str(),"HLT_DoubleMu3_DZ_PFMET50_PFMHT60_v") &&  triggerH->accept(i)){HLT_DoubleMu3_DZ_PFMET50_PFMHT60_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu3_DZ_PFMET50_PFMHT60_v") && !triggerH->accept(i)){HLT_DoubleMu3_DZ_PFMET50_PFMHT60_v = false;};
// /* 4.5%*/if (strstr(TName.c_str(),"HLT_DoubleMu3_DZ_PFMET70_PFMHT70_v") &&  triggerH->accept(i)){HLT_DoubleMu3_DZ_PFMET70_PFMHT70_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu3_DZ_PFMET70_PFMHT70_v") && !triggerH->accept(i)){HLT_DoubleMu3_DZ_PFMET70_PFMHT70_v = false;};
// /* 3%*/if (strstr(TName.c_str(),"HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v") &&  triggerH->accept(i)){HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v") && !triggerH->accept(i)){HLT_DoubleMu3_DZ_PFMET90_PFMHT90_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v") &&  triggerH->accept(i)){HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v") && !triggerH->accept(i)){HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu4_Jpsi_Displaced_v") &&  triggerH->accept(i)){HLT_DoubleMu4_Jpsi_Displaced_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu4_Jpsi_Displaced_v") && !triggerH->accept(i)){HLT_DoubleMu4_Jpsi_Displaced_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu4_Jpsi_NoVertexing_v") &&  triggerH->accept(i)){HLT_DoubleMu4_Jpsi_NoVertexing_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu4_Jpsi_NoVertexing_v") && !triggerH->accept(i)){HLT_DoubleMu4_Jpsi_NoVertexing_v = false;};
// /* Null Efficacity*/if (strstr(TName.c_str(),"HLT_DoubleMu4_JpsiTrkTrk_Displaced_v") &&  triggerH->accept(i)){HLT_DoubleMu4_JpsiTrkTrk_Displaced_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu4_JpsiTrkTrk_Displaced_v") && !triggerH->accept(i)){HLT_DoubleMu4_JpsiTrkTrk_Displaced_v = false;};
// /* 36%*/if (strstr(TName.c_str(),"HLT_DoubleMu43NoFiltersNoVtx_v") &&  triggerH->accept(i)){HLT_DoubleMu43NoFiltersNoVtx_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu43NoFiltersNoVtx_v") && !triggerH->accept(i)){HLT_DoubleMu43NoFiltersNoVtx_v = false;};
// /* 31%*/if (strstr(TName.c_str(),"HLT_DoubleMu48NoFiltersNoVtx_v") &&  triggerH->accept(i)){HLT_DoubleMu48NoFiltersNoVtx_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu48NoFiltersNoVtx_v") && !triggerH->accept(i)){HLT_DoubleMu48NoFiltersNoVtx_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_DoubleMu33NoFiltersNoVtxDisplaced_v") &&  triggerH->accept(i)){HLT_DoubleMu33NoFiltersNoVtxDisplaced_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu33NoFiltersNoVtxDisplaced_v") && !triggerH->accept(i)){HLT_DoubleMu33NoFiltersNoVtxDisplaced_v = false;};
// /* 1%*/if (strstr(TName.c_str(),"HLT_DoubleMu40NoFiltersNoVtxDisplaced_v") &&  triggerH->accept(i)){HLT_DoubleMu40NoFiltersNoVtxDisplaced_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu40NoFiltersNoVtxDisplaced_v") && !triggerH->accept(i)){HLT_DoubleMu40NoFiltersNoVtxDisplaced_v = false;};
// /* 4.4%*/if (strstr(TName.c_str(),"HLT_DoubleMu20_7_Mass0to30_L1_DM4_v") &&  triggerH->accept(i)){HLT_DoubleMu20_7_Mass0to30_L1_DM4_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu20_7_Mass0to30_L1_DM4_v") && !triggerH->accept(i)){HLT_DoubleMu20_7_Mass0to30_L1_DM4_v = false;};
// /* 4.5%*/if (strstr(TName.c_str(),"HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v") &&  triggerH->accept(i)){HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v") && !triggerH->accept(i)){HLT_DoubleMu20_7_Mass0to30_L1_DM4EG_v = false;};
// /* 1.3%*/if (strstr(TName.c_str(),"HLT_DoubleMu20_7_Mass0to30_Photon23_v") &&  triggerH->accept(i)){HLT_DoubleMu20_7_Mass0to30_Photon23_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu20_7_Mass0to30_Photon23_v") && !triggerH->accept(i)){HLT_DoubleMu20_7_Mass0to30_Photon23_v = false;};
// /* Null Efficacity%*/if (strstr(TName.c_str(),"HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_v") &&  triggerH->accept(i)){HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_v") && !triggerH->accept(i)){HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05_v = false;};
// /* Null Efficacity%*/if (strstr(TName.c_str(),"HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v") &&  triggerH->accept(i)){HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v") && !triggerH->accept(i)){HLT_DoubleMu2_Jpsi_DoubleTkMu0_Phi_v = false;};
// /* 10%*/if (strstr(TName.c_str(),"HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v") &&  triggerH->accept(i)){HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v") && !triggerH->accept(i)){HLT_DoubleMu3_DCA_PFMET50_PFMHT60_v = false;};
// // ----------------Trigger DoubleEle-------------
// if (strstr(TName.c_str(),"HLT_DoubleEle25_CaloIdL_MW_v") &&  triggerH->accept(i)){HLT_DoubleEle25_CaloIdL_MW_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleEle25_CaloIdL_MW_v") && !triggerH->accept(i)){HLT_DoubleEle25_CaloIdL_MW_v = false;};
// if (strstr(TName.c_str(),"HLT_DoubleEle27_CaloIdL_MW_v") &&  triggerH->accept(i)){HLT_DoubleEle27_CaloIdL_MW_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleEle27_CaloIdL_MW_v") && !triggerH->accept(i)){HLT_DoubleEle27_CaloIdL_MW_v = false;};
// if (strstr(TName.c_str(),"HLT_DoubleEle33_CaloIdL_MW_v") &&  triggerH->accept(i)){HLT_DoubleEle33_CaloIdL_MW_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleEle33_CaloIdL_MW_v") && !triggerH->accept(i)){HLT_DoubleEle33_CaloIdL_MW_v = false;};
// if (strstr(TName.c_str(),"HLT_DoubleEle24_eta2p1_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_DoubleEle24_eta2p1_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleEle24_eta2p1_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_DoubleEle24_eta2p1_WPTight_Gsf_v = false;};
// if (strstr(TName.c_str(),"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v") &&  triggerH->accept(i)){HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v") && !triggerH->accept(i)){HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_DZ_PFHT350_v = false;};
// if (strstr(TName.c_str(),"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v") &&  triggerH->accept(i)){HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v = true;} else if (strstr(TName.c_str(),"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v") && !triggerH->accept(i)){HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350_v = false;};
// // ----------------Trigger Dimuon0-------------
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_L1_NoOS_v") &&  triggerH->accept(i)){HLT_Dimuon0_Jpsi_L1_NoOS_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_L1_NoOS_v") && !triggerH->accept(i)){HLT_Dimuon0_Jpsi_L1_NoOS_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v") &&  triggerH->accept(i)){HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v") && !triggerH->accept(i)){HLT_Dimuon0_Jpsi_NoVertexing_NoOS_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_v") &&  triggerH->accept(i)){HLT_Dimuon0_Jpsi_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_v") && !triggerH->accept(i)){HLT_Dimuon0_Jpsi_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_NoVertexing_v") &&  triggerH->accept(i)){HLT_Dimuon0_Jpsi_NoVertexing_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_NoVertexing_v") && !triggerH->accept(i)){HLT_Dimuon0_Jpsi_NoVertexing_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v") &&  triggerH->accept(i)){HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v") && !triggerH->accept(i)){HLT_Dimuon0_Jpsi_L1_4R_0er1p5R_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v") &&  triggerH->accept(i)){HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v") && !triggerH->accept(i)){HLT_Dimuon0_Jpsi_NoVertexing_L1_4R_0er1p5R_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi3p5_Muon2_v") &&  triggerH->accept(i)){HLT_Dimuon0_Jpsi3p5_Muon2_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Jpsi3p5_Muon2_v") && !triggerH->accept(i)){HLT_Dimuon0_Jpsi3p5_Muon2_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_4p5_v") &&  triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_4p5_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_4p5_v") && !triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_4p5_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_5_v") &&  triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_5_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_5_v") && !triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_5_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_4p5NoOS_v") &&  triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_4p5NoOS_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_4p5NoOS_v") && !triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_4p5NoOS_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_4p5er2p0_v") &&  triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_4p5er2p0_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_4p5er2p0_v") && !triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_4p5er2p0_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v") &&  triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v") && !triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_4p5er2p0M_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_NoVertexing_v") &&  triggerH->accept(i)){HLT_Dimuon0_Upsilon_NoVertexing_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_NoVertexing_v") && !triggerH->accept(i)){HLT_Dimuon0_Upsilon_NoVertexing_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_5M_v") &&  triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_5M_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_L1_5M_v") && !triggerH->accept(i)){HLT_Dimuon0_Upsilon_L1_5M_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_0er1p5R_v") &&  triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_0er1p5R_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_0er1p5R_v") && !triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_0er1p5R_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_0er1p5_v") &&  triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_0er1p5_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_0er1p5_v") && !triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_0er1p5_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_v") &&  triggerH->accept(i)){HLT_Dimuon0_LowMass_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_v") && !triggerH->accept(i)){HLT_Dimuon0_LowMass_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_4_v") &&  triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_4_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_4_v") && !triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_4_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_4R_v") &&  triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_4R_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_4R_v") && !triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_4R_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_TM530_v") &&  triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_TM530_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_LowMass_L1_TM530_v") && !triggerH->accept(i)){HLT_Dimuon0_LowMass_L1_TM530_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_Muon_L1_TM0_v") &&  triggerH->accept(i)){HLT_Dimuon0_Upsilon_Muon_L1_TM0_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_Muon_L1_TM0_v") && !triggerH->accept(i)){HLT_Dimuon0_Upsilon_Muon_L1_TM0_v = false;};
// /*Null efficacity*/if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v") &&  triggerH->accept(i)){HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v = true;} else if (strstr(TName.c_str(),"HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v") && !triggerH->accept(i)){HLT_Dimuon0_Upsilon_Muon_NoL1Mass_v = false;};
// // ----------------Trigger PFMET-------------
// /* 8.8%*/if (strstr(TName.c_str(),"HLT_PFMET110_PFMHT110_IDTight_v") &&  triggerH->accept(i)){HLT_PFMET110_PFMHT110_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET110_PFMHT110_IDTight_v") && !triggerH->accept(i)){HLT_PFMET110_PFMHT110_IDTight_v = false;};
/* 7.6%*/if (strstr(TName.c_str(),"HLT_PFMET120_PFMHT120_IDTight_v") &&  triggerH->accept(i)){HLT_PFMET120_PFMHT120_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET120_PFMHT120_IDTight_v") && !triggerH->accept(i)){HLT_PFMET120_PFMHT120_IDTight_v = false;};
// /* 6.5%*/if (strstr(TName.c_str(),"HLT_PFMET130_PFMHT130_IDTight_v") &&  triggerH->accept(i)){HLT_PFMET130_PFMHT130_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET130_PFMHT130_IDTight_v") && !triggerH->accept(i)){HLT_PFMET130_PFMHT130_IDTight_v = false;};
// /* 5.5%*/if (strstr(TName.c_str(),"HLT_PFMET140_PFMHT140_IDTight_v") &&  triggerH->accept(i)){HLT_PFMET140_PFMHT140_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET140_PFMHT140_IDTight_v") && !triggerH->accept(i)){HLT_PFMET140_PFMHT140_IDTight_v = false;};
// /* 2.6%*/if (strstr(TName.c_str(),"HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_v") &&  triggerH->accept(i)){HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_v") && !triggerH->accept(i)){HLT_PFMET100_PFMHT100_IDTight_CaloBTagDeepCSV_3p1_v = false;};
// /*2.3%*/if (strstr(TName.c_str(),"HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_v") &&  triggerH->accept(i)){HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_v") && !triggerH->accept(i)){HLT_PFMET110_PFMHT110_IDTight_CaloBTagDeepCSV_3p1_v = false;};
// /* 2.1%*/if (strstr(TName.c_str(),"HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_v") &&  triggerH->accept(i)){HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_v") && !triggerH->accept(i)){HLT_PFMET120_PFMHT120_IDTight_CaloBTagDeepCSV_3p1_v = false;};
// /* 1.9%*/if (strstr(TName.c_str(),"HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_v") &&  triggerH->accept(i)){HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_v") && !triggerH->accept(i)){HLT_PFMET130_PFMHT130_IDTight_CaloBTagDeepCSV_3p1_v = false;};
// /* 1.5%*/if (strstr(TName.c_str(),"HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_v") &&  triggerH->accept(i)){HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_v") && !triggerH->accept(i)){HLT_PFMET140_PFMHT140_IDTight_CaloBTagDeepCSV_3p1_v = false;};
/* 7.8%*/if (strstr(TName.c_str(),"HLT_PFMET120_PFMHT120_IDTight_PFHT60_v") &&  triggerH->accept(i)){HLT_PFMET120_PFMHT120_IDTight_PFHT60_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET120_PFMHT120_IDTight_PFHT60_v") && !triggerH->accept(i)){HLT_PFMET120_PFMHT120_IDTight_PFHT60_v = false;};
/* 15%*/if (strstr(TName.c_str(),"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v") &&  triggerH->accept(i)){HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v") && !triggerH->accept(i)){HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v = false;};
// /*8.8%*/if (strstr(TName.c_str(),"HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v") &&  triggerH->accept(i)){HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v") && !triggerH->accept(i)){HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v = false;};
// /* 10%*/if (strstr(TName.c_str(),"HLT_PFMETTypeOne110_PFMHT110_IDTight_v") &&  triggerH->accept(i)){HLT_PFMETTypeOne110_PFMHT110_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETTypeOne110_PFMHT110_IDTight_v") && !triggerH->accept(i)){HLT_PFMETTypeOne110_PFMHT110_IDTight_v = false;};
// /* 8.5%*/if (strstr(TName.c_str(),"HLT_PFMETTypeOne120_PFMHT120_IDTight_v") &&  triggerH->accept(i)){HLT_PFMETTypeOne120_PFMHT120_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETTypeOne120_PFMHT120_IDTight_v") && !triggerH->accept(i)){HLT_PFMETTypeOne120_PFMHT120_IDTight_v = false;};
// /* 7.5%*/if (strstr(TName.c_str(),"HLT_PFMETTypeOne130_PFMHT130_IDTight_v") &&  triggerH->accept(i)){HLT_PFMETTypeOne130_PFMHT130_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETTypeOne130_PFMHT130_IDTight_v") && !triggerH->accept(i)){HLT_PFMETTypeOne130_PFMHT130_IDTight_v = false;};
// /* 6%*/if (strstr(TName.c_str(),"HLT_PFMETTypeOne140_PFMHT140_IDTight_v") &&  triggerH->accept(i)){HLT_PFMETTypeOne140_PFMHT140_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETTypeOne140_PFMHT140_IDTight_v") && !triggerH->accept(i)){HLT_PFMETTypeOne140_PFMHT140_IDTight_v = false;};
// /* 17%*/if (strstr(TName.c_str(),"HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") &&  triggerH->accept(i)){HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v") && !triggerH->accept(i)){HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v = false;};
/* 14.6%*/if (strstr(TName.c_str(),"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") &&  triggerH->accept(i)){HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") && !triggerH->accept(i)){HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v = false;};
// /* 12%*/if (strstr(TName.c_str(),"HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v") &&  triggerH->accept(i)){HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v") && !triggerH->accept(i)){HLT_PFMETNoMu130_PFMHTNoMu130_IDTight_v = false;};
// /* 11%*/if (strstr(TName.c_str(),"HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v") &&  triggerH->accept(i)){HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v") && !triggerH->accept(i)){HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v = false;};
// /* 3.6%*/if (strstr(TName.c_str(),"HLT_PFMET200_NotCleaned_v") &&  triggerH->accept(i)){HLT_PFMET200_NotCleaned_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET200_NotCleaned_v") && !triggerH->accept(i)){HLT_PFMET200_NotCleaned_v = false;};
// /* 3.6%*/if (strstr(TName.c_str(),"HLT_PFMET200_HBHECleaned_v") &&  triggerH->accept(i)){HLT_PFMET200_HBHECleaned_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET200_HBHECleaned_v") && !triggerH->accept(i)){HLT_PFMET200_HBHECleaned_v = false;};
/* 1%*/if (strstr(TName.c_str(),"HLT_PFMET250_HBHECleaned_v") &&  triggerH->accept(i)){HLT_PFMET250_HBHECleaned_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET250_HBHECleaned_v") && !triggerH->accept(i)){HLT_PFMET250_HBHECleaned_v = false;};
// /* <1%*/if (strstr(TName.c_str(),"HLT_PFMET300_HBHECleaned_v") &&  triggerH->accept(i)){HLT_PFMET300_HBHECleaned_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET300_HBHECleaned_v") && !triggerH->accept(i)){HLT_PFMET300_HBHECleaned_v = false;};
// /* 3.6%*/if (strstr(TName.c_str(),"HLT_PFMET200_HBHE_BeamHaloCleaned_v") &&  triggerH->accept(i)){HLT_PFMET200_HBHE_BeamHaloCleaned_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET200_HBHE_BeamHaloCleaned_v") && !triggerH->accept(i)){HLT_PFMET200_HBHE_BeamHaloCleaned_v = false;};
/* 4%*/if (strstr(TName.c_str(),"HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v") &&  triggerH->accept(i)){HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v") && !triggerH->accept(i)){HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v = false;};
// /* 10%*/if (strstr(TName.c_str(),"HLT_PFMET100_PFMHT100_IDTight_PFHT60_v") &&  triggerH->accept(i)){HLT_PFMET100_PFMHT100_IDTight_PFHT60_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET100_PFMHT100_IDTight_PFHT60_v") && !triggerH->accept(i)){HLT_PFMET100_PFMHT100_IDTight_PFHT60_v = false;};
// /* 20%*/if (strstr(TName.c_str(),"HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v") &&  triggerH->accept(i)){HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v") && !triggerH->accept(i)){HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_PFHT60_v = false;};
// /* 12%*/if (strstr(TName.c_str(),"HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v") &&  triggerH->accept(i)){HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v") && !triggerH->accept(i)){HLT_PFMETTypeOne100_PFMHT100_IDTight_PFHT60_v = false;};
// // // ----------------Trigger HT-------------
// /* 47%*/if (strstr(TName.c_str(),"HLT_HT450_Beamspot_v") &&  triggerH->accept(i)){HLT_HT450_Beamspot_v = true;} else if (strstr(TName.c_str(),"HLT_HT450_Beamspot_v") && !triggerH->accept(i)){HLT_HT450_Beamspot_v = false;};
// /* 72%*/if (strstr(TName.c_str(),"HLT_HT300_Beamspot_v") &&  triggerH->accept(i)){HLT_HT300_Beamspot_v = true;} else if (strstr(TName.c_str(),"HLT_HT300_Beamspot_v") && !triggerH->accept(i)){HLT_HT300_Beamspot_v = false;};
// /* 50%*/if (strstr(TName.c_str(),"HLT_HT425_v") &&  triggerH->accept(i)){HLT_HT425_v = true;} else if (strstr(TName.c_str(),"HLT_HT425_v") && !triggerH->accept(i)){HLT_HT425_v = false;};
// /* <25%*/if (strstr(TName.c_str(),"HLT_HT430_DisplacedDijet40_DisplacedTrack_v") &&  triggerH->accept(i)){HLT_HT430_DisplacedDijet40_DisplacedTrack_v = true;} else if (strstr(TName.c_str(),"HLT_HT430_DisplacedDijet40_DisplacedTrack_v") && !triggerH->accept(i)){HLT_HT430_DisplacedDijet40_DisplacedTrack_v = false;};
// /* <25%*/if (strstr(TName.c_str(),"HLT_HT500_DisplacedDijet40_DisplacedTrack_v") &&  triggerH->accept(i)){HLT_HT500_DisplacedDijet40_DisplacedTrack_v = true;} else if (strstr(TName.c_str(),"HLT_HT500_DisplacedDijet40_DisplacedTrack_v") && !triggerH->accept(i)){HLT_HT500_DisplacedDijet40_DisplacedTrack_v = false;};
// /* <25%*/if (strstr(TName.c_str(),"HLT_HT430_DisplacedDijet60_DisplacedTrack_v") &&  triggerH->accept(i)){HLT_HT430_DisplacedDijet60_DisplacedTrack_v = true;} else if (strstr(TName.c_str(),"HLT_HT430_DisplacedDijet60_DisplacedTrack_v") && !triggerH->accept(i)){HLT_HT430_DisplacedDijet60_DisplacedTrack_v = false;};
// /* <25%*/if (strstr(TName.c_str(),"HLT_HT400_DisplacedDijet40_DisplacedTrack_v") &&  triggerH->accept(i)){HLT_HT400_DisplacedDijet40_DisplacedTrack_v = true;} else if (strstr(TName.c_str(),"HLT_HT400_DisplacedDijet40_DisplacedTrack_v") && !triggerH->accept(i)){HLT_HT400_DisplacedDijet40_DisplacedTrack_v = false;};
// /* <25%*/if (strstr(TName.c_str(),"HLT_HT650_DisplacedDijet60_Inclusive_v") &&  triggerH->accept(i)){HLT_HT650_DisplacedDijet60_Inclusive_v = true;} else if (strstr(TName.c_str(),"HLT_HT650_DisplacedDijet60_Inclusive_v") && !triggerH->accept(i)){HLT_HT650_DisplacedDijet60_Inclusive_v = false;};
// /* <25%*/if (strstr(TName.c_str(),"HLT_HT550_DisplacedDijet60_Inclusive_v") &&  triggerH->accept(i)){HLT_HT550_DisplacedDijet60_Inclusive_v = true;} else if (strstr(TName.c_str(),"HLT_HT550_DisplacedDijet60_Inclusive_v") && !triggerH->accept(i)){HLT_HT550_DisplacedDijet60_Inclusive_v = false;};
// // // ----------------Trigger AK4-------------
// /* 99.9%*/if (strstr(TName.c_str(),"HLT_AK4CaloJet30_v") &&  triggerH->accept(i)){HLT_AK4CaloJet30_v = true;} else if (strstr(TName.c_str(),"HLT_AK4CaloJet30_v") && !triggerH->accept(i)){HLT_AK4CaloJet30_v = false;};
// /* 99.9%*/if (strstr(TName.c_str(),"HLT_AK4CaloJet40_v") &&  triggerH->accept(i)){HLT_AK4CaloJet40_v = true;} else if (strstr(TName.c_str(),"HLT_AK4CaloJet40_v") && !triggerH->accept(i)){HLT_AK4CaloJet40_v = false;};
// /* 99.5%*/if (strstr(TName.c_str(),"HLT_AK4CaloJet50_v") &&  triggerH->accept(i)){HLT_AK4CaloJet50_v = true;} else if (strstr(TName.c_str(),"HLT_AK4CaloJet50_v") && !triggerH->accept(i)){HLT_AK4CaloJet50_v = false;};
// /* 93%*/if (strstr(TName.c_str(),"HLT_AK4CaloJet80_v") &&  triggerH->accept(i)){HLT_AK4CaloJet80_v = true;} else if (strstr(TName.c_str(),"HLT_AK4CaloJet80_v") && !triggerH->accept(i)){HLT_AK4CaloJet80_v = false;};
// /* 84%*/if (strstr(TName.c_str(),"HLT_AK4CaloJet100_v") &&  triggerH->accept(i)){HLT_AK4CaloJet100_v = true;} else if (strstr(TName.c_str(),"HLT_AK4CaloJet100_v") && !triggerH->accept(i)){HLT_AK4CaloJet100_v = false;};
// /* 73%*/if (strstr(TName.c_str(),"HLT_AK4CaloJet120_v") &&  triggerH->accept(i)){HLT_AK4CaloJet120_v = true;} else if (strstr(TName.c_str(),"HLT_AK4CaloJet120_v") && !triggerH->accept(i)){HLT_AK4CaloJet120_v = false;};
// /* 99.96%*/if (strstr(TName.c_str(),"HLT_AK4PFJet30_v") &&  triggerH->accept(i)){HLT_AK4PFJet30_v = true;} else if (strstr(TName.c_str(),"HLT_AK4PFJet30_v") && !triggerH->accept(i)){HLT_AK4PFJet30_v = false;};
// /* 98.2%*/if (strstr(TName.c_str(),"HLT_AK4PFJet50_v") &&  triggerH->accept(i)){HLT_AK4PFJet50_v = true;} else if (strstr(TName.c_str(),"HLT_AK4PFJet50_v") && !triggerH->accept(i)){HLT_AK4PFJet50_v = false;};
// /* 87%*/if (strstr(TName.c_str(),"HLT_AK4PFJet80_v") &&  triggerH->accept(i)){HLT_AK4PFJet80_v = true;} else if (strstr(TName.c_str(),"HLT_AK4PFJet80_v") && !triggerH->accept(i)){HLT_AK4PFJet80_v = false;};
// /* 75%*/if (strstr(TName.c_str(),"HLT_AK4PFJet100_v") &&  triggerH->accept(i)){HLT_AK4PFJet100_v = true;} else if (strstr(TName.c_str(),"HLT_AK4PFJet100_v") && !triggerH->accept(i)){HLT_AK4PFJet100_v = false;};
// /* 65%*/if (strstr(TName.c_str(),"HLT_AK4PFJet120_v") &&  triggerH->accept(i)){HLT_AK4PFJet120_v = true;} else if (strstr(TName.c_str(),"HLT_AK4PFJet120_v") && !triggerH->accept(i)){HLT_AK4PFJet120_v = false;};
// // ----------------Trigger PFJet-------------
// /* 100%*/if (strstr(TName.c_str(),"HLT_PFJet15_v") &&  triggerH->accept(i)){HLT_PFJet15_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet15_v") && !triggerH->accept(i)){HLT_PFJet15_v = false;};
// /* 99.99%*/if (strstr(TName.c_str(),"HLT_PFJet25_v") &&  triggerH->accept(i)){HLT_PFJet25_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet25_v") && !triggerH->accept(i)){HLT_PFJet25_v = false;};
// /* 99.65%*/if (strstr(TName.c_str(),"HLT_PFJet40_v") &&  triggerH->accept(i)){HLT_PFJet40_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet40_v") && !triggerH->accept(i)){HLT_PFJet40_v = false;};
// /* 96%*/if (strstr(TName.c_str(),"HLT_PFJet60_v") &&  triggerH->accept(i)){HLT_PFJet60_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet60_v") && !triggerH->accept(i)){HLT_PFJet60_v = false;};
// /* 86%*/if (strstr(TName.c_str(),"HLT_PFJet80_v") &&  triggerH->accept(i)){HLT_PFJet80_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet80_v") && !triggerH->accept(i)){HLT_PFJet80_v = false;};
// /* 53%*/ if (strstr(TName.c_str(),"HLT_PFJet140_v") &&  triggerH->accept(i)){HLT_PFJet140_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet140_v") && !triggerH->accept(i)){HLT_PFJet140_v = false;};
// /* 20-50%*/ if (strstr(TName.c_str(),"HLT_PFJet200_v") &&  triggerH->accept(i)){HLT_PFJet200_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet200_v") && !triggerH->accept(i)){HLT_PFJet200_v = false;};
// /* 20-50%*/ if (strstr(TName.c_str(),"HLT_PFJet260_v") &&  triggerH->accept(i)){HLT_PFJet260_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet260_v") && !triggerH->accept(i)){HLT_PFJet260_v = false;};
// /* 20-50%*/ if (strstr(TName.c_str(),"HLT_PFJet320_v") &&  triggerH->accept(i)){HLT_PFJet320_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet320_v") && !triggerH->accept(i)){HLT_PFJet320_v = false;};
// /* 20-50%*/ if (strstr(TName.c_str(),"HLT_PFJet400_v") &&  triggerH->accept(i)){HLT_PFJet400_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet400_v") && !triggerH->accept(i)){HLT_PFJet400_v = false;};
// /* 20-50%*/ if (strstr(TName.c_str(),"HLT_PFJet450_v") &&  triggerH->accept(i)){HLT_PFJet450_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet450_v") && !triggerH->accept(i)){HLT_PFJet450_v = false;};
// /* 20-50%*/ if (strstr(TName.c_str(),"HLT_PFJet500_v") &&  triggerH->accept(i)){HLT_PFJet500_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet500_v") && !triggerH->accept(i)){HLT_PFJet500_v = false;};
// /* 20-50%*/ if (strstr(TName.c_str(),"HLT_PFJet550_v") &&  triggerH->accept(i)){HLT_PFJet550_v = true;} else if (strstr(TName.c_str(),"HLT_PFJet550_v") && !triggerH->accept(i)){HLT_PFJet550_v = false;};
// /* 69%*/ if (strstr(TName.c_str(),"HLT_PFJetFwd15_v") &&  triggerH->accept(i)){HLT_PFJetFwd15_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd15_v") && !triggerH->accept(i)){HLT_PFJetFwd15_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd25_v") &&  triggerH->accept(i)){HLT_PFJetFwd25_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd25_v") && !triggerH->accept(i)){HLT_PFJetFwd25_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd40_v") &&  triggerH->accept(i)){HLT_PFJetFwd40_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd40_v") && !triggerH->accept(i)){HLT_PFJetFwd40_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd60_v") &&  triggerH->accept(i)){HLT_PFJetFwd60_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd60_v") && !triggerH->accept(i)){HLT_PFJetFwd60_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd80_v") &&  triggerH->accept(i)){HLT_PFJetFwd80_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd80_v") && !triggerH->accept(i)){HLT_PFJetFwd80_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd140_v") &&  triggerH->accept(i)){HLT_PFJetFwd140_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd140_v") && !triggerH->accept(i)){HLT_PFJetFwd140_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd200_v") &&  triggerH->accept(i)){HLT_PFJetFwd200_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd200_v") && !triggerH->accept(i)){HLT_PFJetFwd200_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd260_v") &&  triggerH->accept(i)){HLT_PFJetFwd260_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd260_v") && !triggerH->accept(i)){HLT_PFJetFwd260_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd320_v") &&  triggerH->accept(i)){HLT_PFJetFwd320_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd320_v") && !triggerH->accept(i)){HLT_PFJetFwd320_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd400_v") &&  triggerH->accept(i)){HLT_PFJetFwd400_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd400_v") && !triggerH->accept(i)){HLT_PFJetFwd400_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd450_v") &&  triggerH->accept(i)){HLT_PFJetFwd450_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd450_v") && !triggerH->accept(i)){HLT_PFJetFwd450_v = false;};
// /*Null efficacity*/ if (strstr(TName.c_str(),"HLT_PFJetFwd500_v") && triggerH->accept(i)){HLT_PFJetFwd500_v = true;} else if (strstr(TName.c_str(),"HLT_PFJetFwd500_v") && !triggerH->accept(i)){HLT_PFJetFwd500_v = false;}
// // ----------------Trigger DiPFJet-------------
// /*99.79%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve40_v") && triggerH->accept(i)){HLT_DiPFJetAve40_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve40_v") &&!triggerH->accept(i)){HLT_DiPFJetAve40_v = false;};
// /*96.5%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve60_v") && triggerH->accept(i)){HLT_DiPFJetAve60_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve60_v") &&!triggerH->accept(i)){HLT_DiPFJetAve60_v = false;};
// /*86%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve80_v") && triggerH->accept(i)){HLT_DiPFJetAve80_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve80_v") &&!triggerH->accept(i)){HLT_DiPFJetAve80_v = false;};
// /*46%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve140_v") && triggerH->accept(i)){HLT_DiPFJetAve140_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve140_v") &&!triggerH->accept(i)){HLT_DiPFJetAve140_v = false;};
// /*27%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve200_v") && triggerH->accept(i)){HLT_DiPFJetAve200_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve200_v") &&!triggerH->accept(i)){HLT_DiPFJetAve200_v = false;};
// /*17%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve260_v") && triggerH->accept(i)){HLT_DiPFJetAve260_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve260_v") &&!triggerH->accept(i)){HLT_DiPFJetAve260_v = false;};
// /*12%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve320_v") && triggerH->accept(i)){HLT_DiPFJetAve320_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve320_v") &&!triggerH->accept(i)){HLT_DiPFJetAve320_v = false;};
// /*7%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve400_v") && triggerH->accept(i)){HLT_DiPFJetAve400_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve400_v") &&!triggerH->accept(i)){HLT_DiPFJetAve400_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve500_v") && triggerH->accept(i)){HLT_DiPFJetAve500_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve500_v") &&!triggerH->accept(i)){HLT_DiPFJetAve500_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve60_HFJEC_v") && triggerH->accept(i)){HLT_DiPFJetAve60_HFJEC_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve60_HFJEC_v") &&!triggerH->accept(i)){HLT_DiPFJetAve60_HFJEC_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve80_HFJEC_v") && triggerH->accept(i)){HLT_DiPFJetAve80_HFJEC_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve80_HFJEC_v") &&!triggerH->accept(i)){HLT_DiPFJetAve80_HFJEC_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve100_HFJEC_v") && triggerH->accept(i)){HLT_DiPFJetAve100_HFJEC_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve100_HFJEC_v") &&!triggerH->accept(i)){HLT_DiPFJetAve100_HFJEC_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve160_HFJEC_v") && triggerH->accept(i)){HLT_DiPFJetAve160_HFJEC_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve160_HFJEC_v") &&!triggerH->accept(i)){HLT_DiPFJetAve160_HFJEC_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve220_HFJEC_v") && triggerH->accept(i)){HLT_DiPFJetAve220_HFJEC_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve220_HFJEC_v") &&!triggerH->accept(i)){HLT_DiPFJetAve220_HFJEC_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DiPFJetAve300_HFJEC_v") && triggerH->accept(i)){HLT_DiPFJetAve300_HFJEC_v = true;} else if (strstr(TName.c_str(),"HLT_DiPFJetAve300_HFJEC_v") &&!triggerH->accept(i)){HLT_DiPFJetAve300_HFJEC_v = false;};
// // ----------------Trigger DoublePFJets-------------
// /*8%*/if (strstr(TName.c_str(),"HLT_DoublePFJets40_CaloBTagDeepCSV_p71_v") && triggerH->accept(i)){HLT_DoublePFJets40_CaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_DoublePFJets40_CaloBTagDeepCSV_p71_v") &&!triggerH->accept(i)){HLT_DoublePFJets40_CaloBTagDeepCSV_p71_v = false;};
// /*5%*/if (strstr(TName.c_str(),"HLT_DoublePFJets100_CaloBTagDeepCSV_p71_v") && triggerH->accept(i)){HLT_DoublePFJets100_CaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_DoublePFJets100_CaloBTagDeepCSV_p71_v") &&!triggerH->accept(i)){HLT_DoublePFJets100_CaloBTagDeepCSV_p71_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DoublePFJets200_CaloBTagDeepCSV_p71_v") && triggerH->accept(i)){HLT_DoublePFJets200_CaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_DoublePFJets200_CaloBTagDeepCSV_p71_v") &&!triggerH->accept(i)){HLT_DoublePFJets200_CaloBTagDeepCSV_p71_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DoublePFJets350_CaloBTagDeepCSV_p71_v") && triggerH->accept(i)){HLT_DoublePFJets350_CaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_DoublePFJets350_CaloBTagDeepCSV_p71_v") &&!triggerH->accept(i)){HLT_DoublePFJets350_CaloBTagDeepCSV_p71_v = false;};
// if (strstr(TName.c_str(),"HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") && triggerH->accept(i)){HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") &&!triggerH->accept(i)){HLT_DoublePFJets116MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") && triggerH->accept(i)){HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = true;} else if (strstr(TName.c_str(),"HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v") &&!triggerH->accept(i)){HLT_DoublePFJets128MaxDeta1p6_DoubleCaloBTagDeepCSV_p71_v = false;};
// // ----------------Trigger BTagMu-------------
// /*15%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet20_Mu5_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet20_Mu5_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet20_Mu5_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet20_Mu5_v = false;};
// /*15%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet40_Mu5_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet40_Mu5_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet40_Mu5_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet40_Mu5_v = false;};
// /*10%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet70_Mu5_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet70_Mu5_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet70_Mu5_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet70_Mu5_v = false;};
// /*6%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet110_Mu5_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet110_Mu5_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet110_Mu5_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet110_Mu5_v = false;};
// /*4%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet170_Mu5_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet170_Mu5_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet170_Mu5_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet170_Mu5_v = false;};
// /*2%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4Jet300_Mu5_v") && triggerH->accept(i)){HLT_BTagMu_AK4Jet300_Mu5_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4Jet300_Mu5_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4Jet300_Mu5_v = false;};
// /*6%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK8DiJet170_Mu5_v") && triggerH->accept(i)){HLT_BTagMu_AK8DiJet170_Mu5_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK8DiJet170_Mu5_v") &&!triggerH->accept(i)){HLT_BTagMu_AK8DiJet170_Mu5_v = false;};
// /*4%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK8Jet170_DoubleMu5_v") && triggerH->accept(i)){HLT_BTagMu_AK8Jet170_DoubleMu5_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK8Jet170_DoubleMu5_v") &&!triggerH->accept(i)){HLT_BTagMu_AK8Jet170_DoubleMu5_v = false;};
// /*4%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK8Jet300_Mu5_v") && triggerH->accept(i)){HLT_BTagMu_AK8Jet300_Mu5_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK8Jet300_Mu5_v") &&!triggerH->accept(i)){HLT_BTagMu_AK8Jet300_Mu5_v = false;};
// /*40%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet20_Mu5_noalgo_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet20_Mu5_noalgo_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet20_Mu5_noalgo_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet20_Mu5_noalgo_v = false;};
// /*36%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet40_Mu5_noalgo_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet40_Mu5_noalgo_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet40_Mu5_noalgo_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet40_Mu5_noalgo_v = false;};
// /*28%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet70_Mu5_noalgo_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet70_Mu5_noalgo_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet70_Mu5_noalgo_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet70_Mu5_noalgo_v = false;};
// /*20%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet110_Mu5_noalgo_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet110_Mu5_noalgo_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet110_Mu5_noalgo_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet110_Mu5_noalgo_v = false;};
// /*10%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet170_Mu5_noalgo_v") && triggerH->accept(i)){HLT_BTagMu_AK4DiJet170_Mu5_noalgo_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4DiJet170_Mu5_noalgo_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4DiJet170_Mu5_noalgo_v = false;};
// /*10%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK4Jet300_Mu5_noalgo_v") && triggerH->accept(i)){HLT_BTagMu_AK4Jet300_Mu5_noalgo_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK4Jet300_Mu5_noalgo_v") &&!triggerH->accept(i)){HLT_BTagMu_AK4Jet300_Mu5_noalgo_v = false;};
// /*10%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v") && triggerH->accept(i)){HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v") &&!triggerH->accept(i)){HLT_BTagMu_AK8DiJet170_Mu5_noalgo_v = false;};
// /*10%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_v") && triggerH->accept(i)){HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_v") &&!triggerH->accept(i)){HLT_BTagMu_AK8Jet170_DoubleMu5_noalgo_v = false;};
// /*10%*/if (strstr(TName.c_str(),"HLT_BTagMu_AK8Jet300_Mu5_noalgo_v") && triggerH->accept(i)){HLT_BTagMu_AK8Jet300_Mu5_noalgo_v = true;} else if (strstr(TName.c_str(),"HLT_BTagMu_AK8Jet300_Mu5_noalgo_v") &&!triggerH->accept(i)){HLT_BTagMu_AK8Jet300_Mu5_noalgo_v = false;};
// // ----------------Trigger QuadPFJet-------------
// /*1%*/if (strstr(TName.c_str(),"HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v") && triggerH->accept(i)){HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v") &&!triggerH->accept(i)){HLT_QuadPFJet98_83_71_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v") && triggerH->accept(i)){HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v") &&!triggerH->accept(i)){HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v") && triggerH->accept(i)){HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v") &&!triggerH->accept(i)){HLT_QuadPFJet111_90_80_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_v") && triggerH->accept(i)){HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_v") &&!triggerH->accept(i)){HLT_QuadPFJet98_83_71_15_PFBTagDeepCSV_1p3_VBF2_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_v") && triggerH->accept(i)){HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_v") &&!triggerH->accept(i)){HLT_QuadPFJet103_88_75_15_PFBTagDeepCSV_1p3_VBF2_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_v") && triggerH->accept(i)){HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_v") &&!triggerH->accept(i)){HLT_QuadPFJet105_88_76_15_PFBTagDeepCSV_1p3_VBF2_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_v") && triggerH->accept(i)){HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_v") &&!triggerH->accept(i)){HLT_QuadPFJet111_90_80_15_PFBTagDeepCSV_1p3_VBF2_v = false;};
// /*20%*/if (strstr(TName.c_str(),"HLT_QuadPFJet98_83_71_15_v") && triggerH->accept(i)){HLT_QuadPFJet98_83_71_15_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet98_83_71_15_v") &&!triggerH->accept(i)){HLT_QuadPFJet98_83_71_15_v = false;};
// /*20%*/if (strstr(TName.c_str(),"HLT_QuadPFJet103_88_75_15_v") && triggerH->accept(i)){HLT_QuadPFJet103_88_75_15_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet103_88_75_15_v") &&!triggerH->accept(i)){HLT_QuadPFJet103_88_75_15_v = false;};
// /*20%*/if (strstr(TName.c_str(),"HLT_QuadPFJet105_88_76_15_v") && triggerH->accept(i)){HLT_QuadPFJet105_88_76_15_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet105_88_76_15_v") &&!triggerH->accept(i)){HLT_QuadPFJet105_88_76_15_v = false;};
// /*20%*/if (strstr(TName.c_str(),"HLT_QuadPFJet111_90_80_15_v") && triggerH->accept(i)){HLT_QuadPFJet111_90_80_15_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet111_90_80_15_v") &&!triggerH->accept(i)){HLT_QuadPFJet111_90_80_15_v = false;};
// /*1%*/if (strstr(TName.c_str(),"HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v") && triggerH->accept(i)){HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v = true;} else if (strstr(TName.c_str(),"HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v") &&!triggerH->accept(i)){HLT_QuadPFJet105_88_76_15_DoublePFBTagDeepCSV_1p3_7p7_VBF1_v = false;};
// // ----------------Trigger IsoMu-------------
// if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu20_eta2p1_MediumChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu20_eta2p1_TightChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu24_eta2p1_TightChargedIsoPFTauHPS35_Trk1_TightID_eta2p1_Reg_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v") && triggerH->accept(i)){HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v") &&!triggerH->accept(i)){HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v") && triggerH->accept(i)){HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v") &&!triggerH->accept(i)){HLT_IsoMu27_LooseChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v") && triggerH->accept(i)){HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v") &&!triggerH->accept(i)){HLT_IsoMu27_MediumChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v") && triggerH->accept(i)){HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v") &&!triggerH->accept(i)){HLT_IsoMu27_TightChargedIsoPFTauHPS20_Trk1_eta2p1_SingleL1_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu20_v") && triggerH->accept(i)){HLT_IsoMu20_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu20_v") &&!triggerH->accept(i)){HLT_IsoMu20_v = false;};
if (strstr(TName.c_str(),"HLT_IsoMu24_v") && triggerH->accept(i)){HLT_IsoMu24_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu24_v") &&!triggerH->accept(i)){HLT_IsoMu24_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_v") && triggerH->accept(i)){HLT_IsoMu24_eta2p1_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_v") &&!triggerH->accept(i)){HLT_IsoMu24_eta2p1_v = false;};
if (strstr(TName.c_str(),"HLT_IsoMu27_v") && triggerH->accept(i)){HLT_IsoMu27_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu27_v") &&!triggerH->accept(i)){HLT_IsoMu27_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu30_v") && triggerH->accept(i)){HLT_IsoMu30_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu30_v") &&!triggerH->accept(i)){HLT_IsoMu30_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu24_TwoProngs35_v") && triggerH->accept(i)){HLT_IsoMu24_TwoProngs35_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu24_TwoProngs35_v") &&!triggerH->accept(i)){HLT_IsoMu24_TwoProngs35_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v") && triggerH->accept(i)){HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v") &&!triggerH->accept(i)){HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_v = false;};
// if (strstr(TName.c_str(),"HLT_IsoMu27_MET90_v") && triggerH->accept(i)){HLT_IsoMu27_MET90_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu27_MET90_v") &&!triggerH->accept(i)){HLT_IsoMu27_MET90_v = false;};
// // ----------------Trigger PFHT-------------
// if (strstr(TName.c_str(),"HLT_PFHT180_v") && triggerH->accept(i)){HLT_PFHT180_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT180_v") &&!triggerH->accept(i)){HLT_PFHT180_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT250_v") && triggerH->accept(i)){HLT_PFHT250_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT250_v") &&!triggerH->accept(i)){HLT_PFHT250_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT370_v") && triggerH->accept(i)){HLT_PFHT370_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT370_v") &&!triggerH->accept(i)){HLT_PFHT370_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT430_v") && triggerH->accept(i)){HLT_PFHT430_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT430_v") &&!triggerH->accept(i)){HLT_PFHT430_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT510_v") && triggerH->accept(i)){HLT_PFHT510_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT510_v") &&!triggerH->accept(i)){HLT_PFHT510_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT590_v") && triggerH->accept(i)){HLT_PFHT590_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT590_v") &&!triggerH->accept(i)){HLT_PFHT590_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT680_v") && triggerH->accept(i)){HLT_PFHT680_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT680_v") &&!triggerH->accept(i)){HLT_PFHT680_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT780_v") && triggerH->accept(i)){HLT_PFHT780_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT780_v") &&!triggerH->accept(i)){HLT_PFHT780_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT890_v") && triggerH->accept(i)){HLT_PFHT890_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT890_v") &&!triggerH->accept(i)){HLT_PFHT890_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT1050_v") && triggerH->accept(i)){HLT_PFHT1050_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT1050_v") &&!triggerH->accept(i)){HLT_PFHT1050_v = false;};
if (strstr(TName.c_str(),"HLT_PFHT500_PFMET100_PFMHT100_IDTight_v") && triggerH->accept(i)){HLT_PFHT500_PFMET100_PFMHT100_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT500_PFMET100_PFMHT100_IDTight_v") &&!triggerH->accept(i)){HLT_PFHT500_PFMET100_PFMHT100_IDTight_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT500_PFMET110_PFMHT110_IDTight_v") && triggerH->accept(i)){HLT_PFHT500_PFMET110_PFMHT110_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT500_PFMET110_PFMHT110_IDTight_v") &&!triggerH->accept(i)){HLT_PFHT500_PFMET110_PFMHT110_IDTight_v = false;};
if (strstr(TName.c_str(),"HLT_PFHT700_PFMET85_PFMHT85_IDTight_v") && triggerH->accept(i)){HLT_PFHT700_PFMET85_PFMHT85_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT700_PFMET85_PFMHT85_IDTight_v") &&!triggerH->accept(i)){HLT_PFHT700_PFMET85_PFMHT85_IDTight_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT700_PFMET95_PFMHT95_IDTight_v") && triggerH->accept(i)){HLT_PFHT700_PFMET95_PFMHT95_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT700_PFMET95_PFMHT95_IDTight_v") &&!triggerH->accept(i)){HLT_PFHT700_PFMET95_PFMHT95_IDTight_v = false;};
if (strstr(TName.c_str(),"HLT_PFHT800_PFMET75_PFMHT75_IDTight_v") && triggerH->accept(i)){HLT_PFHT800_PFMET75_PFMHT75_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT800_PFMET75_PFMHT75_IDTight_v") &&!triggerH->accept(i)){HLT_PFHT800_PFMET75_PFMHT75_IDTight_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT800_PFMET85_PFMHT85_IDTight_v") && triggerH->accept(i)){HLT_PFHT800_PFMET85_PFMHT85_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT800_PFMET85_PFMHT85_IDTight_v") &&!triggerH->accept(i)){HLT_PFHT800_PFMET85_PFMHT85_IDTight_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v") && triggerH->accept(i)){HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v") &&!triggerH->accept(i)){HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepCSV_4p5_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT330PT30_QuadPFJet_75_60_45_40_v") && triggerH->accept(i)){HLT_PFHT330PT30_QuadPFJet_75_60_45_40_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT330PT30_QuadPFJet_75_60_45_40_v") &&!triggerH->accept(i)){HLT_PFHT330PT30_QuadPFJet_75_60_45_40_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_v") && triggerH->accept(i)){HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_v") &&!triggerH->accept(i)){HLT_PFHT400_SixPFJet32_DoublePFBTagDeepCSV_2p94_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT400_SixPFJet32_v") && triggerH->accept(i)){HLT_PFHT400_SixPFJet32_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT400_SixPFJet32_v") &&!triggerH->accept(i)){HLT_PFHT400_SixPFJet32_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_v") && triggerH->accept(i)){HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_v") &&!triggerH->accept(i)){HLT_PFHT450_SixPFJet36_PFBTagDeepCSV_1p59_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT450_SixPFJet36_v") && triggerH->accept(i)){HLT_PFHT450_SixPFJet36_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT450_SixPFJet36_v") &&!triggerH->accept(i)){HLT_PFHT450_SixPFJet36_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT350_v") && triggerH->accept(i)){HLT_PFHT350_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT350_v") &&!triggerH->accept(i)){HLT_PFHT350_v = false;};
// if (strstr(TName.c_str(),"HLT_PFHT350MinPFJet15_v") && triggerH->accept(i)){HLT_PFHT350MinPFJet15_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT350MinPFJet15_v") &&!triggerH->accept(i)){HLT_PFHT350MinPFJet15_v = false;};
//************************************************************************************************************************//

}

//   tree_trigger_size.push_back(TriggerCheck.size()); // has to be the same for each evnet!!

//---Test Patate----//

//------Test One plot for One trigger (For one file only)
// TFile* pFile = new TFile("triggerTest.root","recreate");
// TRandom3 rand3;
// TEfficiency* EffvsObs2;
// TDirectory* drt2 = new TDirectory();
// EffvsObs2 = new TEfficiency("Eff2","Efficiency2;patate2;#epsilon",20,0,10);
// EffvsObs2->SetDirectory(drt2);
// bool bPassed;
// double PATATE;
// for (int i = 0 ; i <10000 ; i++)
//   {
//     PATATE = rand3.Uniform(10);//variable
//     bPassed = rand3.Rndm() < TMath::Gaus(PATATE,5,4);//triggerH->accept(0)//trigger boolean
//     EffvsObs2->TEfficiency::Fill(bPassed,PATATE);
//   }

//trig

//------Test One plot for each trigger
// TEfficiency** EffvsObs3;
// TDirectory* drt3 = new TDirectory();

// for (unsigned int j=0;j< triggerH->size();j++)
//   {
//     EffvsObs3[j] = new TEfficiency("Eff3","Efficiency3;MET;#epsilon",20,0,10);//crash
//     EffvsObs3[j]->SetDirectory(drt3);
//     for (int i = 0 ; i <1000 ; i++)
//       {
//         PATATE = rand3.Uniform(10);//variable
//         bPassed = rand3.Rndm() < TMath::Gaus(PATATE,5,4);//triggerH->accept(0)//trigger boolean
//         EffvsObs3[j]->TEfficiency::Fill(bPassed,PATATE);//triggerH->accept(j)
//       }
//       EffvsObs3[j]->Write();
//   }

//         for (unsigned int i=0;i<tree_passesTriggerName.size();i++)
//             {
//               if (strstr(tree_passesTriggerName[i].c_str(),"HLT_Mu15_v"))
//                 {
//                   test->Fill(tk_pt,triggerH->accept(i));
//                 }
//             }


  //////////////////////////////////
  //////////////////////////////////
  //////    Primary Vertex   ///////
  //////////////////////////////////
  //////////////////////////////////
  
//   for (unsigned int i = 0; i< primaryVertex->size() ; ++i) {
//     tree_vtx_PosX.push_back((*primaryVertex)[i].x());
//     tree_vtx_PosY.push_back((*primaryVertex)[i].y());
//     tree_vtx_PosZ.push_back((*primaryVertex)[i].z());
//     tree_vtx_NChi2.push_back((*primaryVertex)[i].normalizedChi2());
//     tree_vtx_PosXError.push_back((*primaryVertex)[i].xError());
//     tree_vtx_PosYError.push_back((*primaryVertex)[i].yError());
//     tree_vtx_PosZError.push_back((*primaryVertex)[i].zError());
//   }
  tree_nPV = primaryVertex->size();
  if ( !primaryVertex->empty() ) {
    tree_PV_x.push_back(     (*primaryVertex)[0].x()); // l'index 0 donne le PV!
    tree_PV_y.push_back(     (*primaryVertex)[0].y());
    tree_PV_z.push_back(     (*primaryVertex)[0].z());
    tree_PV_ez.push_back(    (*primaryVertex)[0].zError());
    tree_PV_NChi2.push_back( (*primaryVertex)[0].normalizedChi2());
    tree_PV_ndf.push_back(   (*primaryVertex)[0].ndof());
  }
  const reco::Vertex &PV = primaryVertex->front();


  //////////////////////////////////
  //////////////////////////////////
  //////          SV         /////// to replace Ycandidate
  //////////////////////////////////
  //////////////////////////////////

  tree_nSV = SVertex->size();
  if ( !SVertex->empty() ) {
    for (int j=0; j<tree_nSV ; j++)
    {
      float SVx = (*SVertex)[j].position().x();
      float SVy = (*SVertex)[j].position().y();
      tree_SV_x.push_back(     SVx);
      tree_SV_y.push_back(     SVy);
      tree_SV_z.push_back(     (*SVertex)[j].position().z());
      tree_SV_r.push_back(     TMath::Sqrt(SVx*SVx+SVy*SVy));
      tree_SV_NChi2.push_back( (*SVertex)[j].vertexNormalizedChi2());
      tree_SV_ndf.push_back(   (*SVertex)[j].vertexNdof());
      tree_SV_nDaughters.push_back((*SVertex)[j].numberOfDaughters());
      tree_SV_mass.push_back(  (*SVertex)[j].mass());
      tree_SV_pdgid.push_back( (*SVertex)[j].pdgId());//doesnot work
      tree_SV_pt.push_back(    (*SVertex)[j].pt());
      tree_SV_eta.push_back(   (*SVertex)[j].eta());
      tree_SV_phi.push_back(   (*SVertex)[j].phi());
      tree_SV_IsConvertedPhoton.push_back((*SVertex)[j].isConvertedPhoton()); // does not work
//$$$$
//       for (unsigned int i = 0 ; i < (*SVertex)[j].numberOfDaughters() ; i++)
//       {
//         tree_SV_daughter_pdgid.push_back((*SVertex)[j].daughter(i)->pdgId()); // does not work
//       }
//$$$$
      //Potential loop on daughter Candidates : pdgid, status, vx, vy ,vz, phi ;eta, p, mt,charge
      // to match with the tracks of packPFCandidate
    }
  }


  //////////////////////////////////
  //////////////////////////////////
  //////    V0 Candidates     //////
  //////////////////////////////////
  //////////////////////////////////
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideV0Producer

  tree_nK0= KshortVertex->size();
  if ( !KshortVertex->empty() ) {
    for (int j = 0; j<tree_nK0 ; j++ )
    {
      float K0x = (*KshortVertex)[j].position().x();
      float K0y = (*KshortVertex)[j].position().y();
      float K0z = (*KshortVertex)[j].position().z(); 
      tree_K0_x.push_back(	K0x);
      tree_K0_y.push_back(	K0y);
      tree_K0_z.push_back(	K0z);
      tree_Hemi_Vtx_bkg_x.push_back(K0x);
      tree_Hemi_Vtx_bkg_y.push_back(K0y);
      tree_Hemi_Vtx_bkg_z.push_back(K0z);
      tree_K0_r.push_back(	TMath::Sqrt(K0x*K0x+K0y*K0y));
      tree_K0_NChi2.push_back(  (*KshortVertex)[j].vertexNormalizedChi2());
      tree_K0_ndf.push_back(	(*KshortVertex)[j].vertexNdof());
      tree_K0_mass.push_back(	(*KshortVertex)[j].mass());
      tree_K0_pt.push_back(	(*KshortVertex)[j].pt());
      tree_K0_eta.push_back(	(*KshortVertex)[j].eta());
      tree_K0_phi.push_back(	(*KshortVertex)[j].phi());
//$$$$
      if ( (*KshortVertex)[j].vertexNormalizedChi2() < 5. &&
           (*KshortVertex)[j].mass() > 0.48 && (*KshortVertex)[j].mass() < 0.52 ) {
        tree_K0_nDaughters.push_back((*KshortVertex)[j].numberOfDaughters());
//$$$$      
        for (unsigned int i =0 ; i< (*KshortVertex)[j].numberOfDaughters() ; i++)
        {
          tree_K0_daughter_pt.push_back((*KshortVertex)[j].daughter(i)->pt());
          tree_K0_daughter_eta.push_back((*KshortVertex)[j].daughter(i)->eta());
          tree_K0_daughter_phi.push_back((*KshortVertex)[j].daughter(i)->phi());
          tree_track_bkg_pt.push_back((*KshortVertex)[j].daughter(i)->pt());
          tree_track_bkg_eta.push_back((*KshortVertex)[j].daughter(i)->eta());
          tree_track_bkg_phi.push_back((*KshortVertex)[j].daughter(i)->phi());
          tree_track_bkg_charge.push_back((*KshortVertex)[j].daughter(i)->charge());
          tree_track_bkg_source.push_back(1);
        }
//$$$$      
      }
      else  tree_K0_nDaughters.push_back(0);     
//$$$$      
      //Potential loop on daughter Candidates : pdgid, status, vx, vy ,vz, phi ;eta, p, mt,charge
      // to match with the ttrakcs of packPFCandidate
    }
  }

  tree_nLambda = LambdaVertex->size();
  if ( !LambdaVertex->empty() ) {
    for (int j = 0; j<tree_nLambda ; j++ )
    {
      float L0x = (*LambdaVertex)[j].position().x();
      float L0y = (*LambdaVertex)[j].position().y();
      float L0z = (*LambdaVertex)[j].position().z();
      tree_L0_x.push_back(         L0x); 
      tree_L0_y.push_back(         L0y);
      tree_L0_z.push_back(         L0z);
      tree_Hemi_Vtx_bkg_x.push_back(L0x);
      tree_Hemi_Vtx_bkg_y.push_back(L0y);
      tree_Hemi_Vtx_bkg_z.push_back(L0z);
      tree_L0_r.push_back(         TMath::Sqrt(L0x*L0x+L0y*L0y));
      tree_L0_NChi2.push_back(     (*LambdaVertex)[j].vertexNormalizedChi2());
      tree_L0_ndf.push_back(       (*LambdaVertex)[j].vertexNdof());
      tree_L0_mass.push_back(      (*LambdaVertex)[j].mass());
      tree_L0_pt.push_back(        (*LambdaVertex)[j].pt());
      tree_L0_eta.push_back(       (*LambdaVertex)[j].eta());
      tree_L0_phi.push_back(       (*LambdaVertex)[j].phi());
//$$$$
      if ( (*LambdaVertex)[j].vertexNormalizedChi2() < 5. &&
           (*LambdaVertex)[j].mass() > 1.111 && (*LambdaVertex)[j].mass() < 1.121 ) {
        tree_L0_nDaughters.push_back((*LambdaVertex)[j].numberOfDaughters());
//$$$$      
        for (unsigned int i =0 ; i< (*LambdaVertex)[j].numberOfDaughters() ; i++)
        {
          tree_L0_daughter_pt.push_back((*LambdaVertex)[j].daughter(i)->pt());
          tree_L0_daughter_eta.push_back((*LambdaVertex)[j].daughter(i)->eta());
          tree_L0_daughter_phi.push_back((*LambdaVertex)[j].daughter(i)->phi());
          tree_track_bkg_pt.push_back((*LambdaVertex)[j].daughter(i)->pt());
          tree_track_bkg_eta.push_back((*LambdaVertex)[j].daughter(i)->eta());
          tree_track_bkg_phi.push_back((*LambdaVertex)[j].daughter(i)->phi());
          tree_track_bkg_charge.push_back((*LambdaVertex)[j].daughter(i)->charge());
          tree_track_bkg_source.push_back(2);
        }
//$$$$      
      }
      else  tree_L0_nDaughters.push_back(0);     
//$$$$      
      //Potential loop on daughter Candidates : pdgid, status, vx, vy ,vz, phi ;eta, p, mt,charge
      // to match with the ttrakcs of packPFCandidate
    }
  }


  //////////////////////////////////
  //////////////////////////////////
  //////    YConversion      ///////
  //////////////////////////////////
  //////////////////////////////////
// https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/DataFormats/EgammaCandidates/interface/Conversion.h

  tree_nYConv = PhotonConversion->size();
  std::vector<reco::Track> YcTracks;
  if ( !PhotonConversion->empty() ) {
    for (int j = 0; j<tree_nYConv ; j++ )
    {
      if ((*PhotonConversion)[j].isConverted()) // tracksize>0 // Number of tracks= 0,1,2 unsigned int nTracks() const {return  tracks().size(); }
      {
        float Yx = (*PhotonConversion)[j].conversionVertex().position().x();
        float Yy = (*PhotonConversion)[j].conversionVertex().position().y();
        float Yz = (*PhotonConversion)[j].conversionVertex().position().z();
        tree_Yc_x.push_back(        Yx); 
        tree_Yc_y.push_back(        Yy);
        tree_Yc_z.push_back(        Yz);
        tree_Hemi_Vtx_bkg_x.push_back(Yx);
        tree_Hemi_Vtx_bkg_y.push_back(Yy);
        tree_Hemi_Vtx_bkg_z.push_back(Yz);
        tree_Yc_r.push_back(        TMath::Sqrt(Yx*Yx+Yy*Yy));
        tree_Yc_Vtx_NChi2.push_back((*PhotonConversion)[j].conversionVertex().normalizedChi2());
        tree_Yc_Vtx_ndf.push_back(   (*PhotonConversion)[j].conversionVertex().ndof());        
        tree_Yc_tracks_sum3p.push_back((*PhotonConversion)[j].refittedPairMomentum());

        std::vector<math::XYZPointF>  inPos = (*PhotonConversion)[j].tracksInnerPosition();
        std::vector<math::XYZVectorF> inP   = (*PhotonConversion)[j].tracksPin();
        std::vector<math::XYZVectorF> outP  = (*PhotonConversion)[j].tracksPout();
        std::vector<reco::Track> YcVertexTracks = (*PhotonConversion)[j].conversionVertex().refittedTracks();
        for (unsigned int k = 0; k<(*PhotonConversion)[j].nTracks(); k ++)
        {
          YcTracks.push_back(YcVertexTracks[k]);
          tree_Yc_tracks_InPosx.push_back(inPos[k].x());
          tree_Yc_tracks_InPosy.push_back(inPos[k].y());
          tree_Yc_tracks_InPosz.push_back(inPos[k].z());
          tree_Yc_tracks_InPx.push_back(  inP[k].x());
          tree_Yc_tracks_InPy.push_back(  inP[k].y());
          tree_Yc_tracks_InPz.push_back(  inP[k].z());
          tree_Yc_tracks_OutPx.push_back( outP[k].x());
          tree_Yc_tracks_OutPy.push_back( outP[k].y());
          tree_Yc_tracks_OutPz.push_back( outP[k].z());
        }

        if ((*PhotonConversion)[j].nTracks()==2 && (*PhotonConversion)[j].conversionVertex().normalizedChi2()<5)
          {
            if ((*PhotonConversion)[j].pairInvariantMass()<1 )
              {
                  for (unsigned int i=0 ; i < YcTracks.size() ; i++)
                    {
                      tree_Yc_tracks_pt.push_back( YcTracks[i].pt());
                      tree_Yc_tracks_eta.push_back(YcTracks[i].eta());
                      tree_Yc_tracks_phi.push_back(YcTracks[i].phi());
                      tree_track_bkg_pt.push_back( YcTracks[i].pt());
                      tree_track_bkg_eta.push_back(YcTracks[i].eta());
                      tree_track_bkg_phi.push_back(YcTracks[i].phi());
                      tree_track_bkg_charge.push_back(YcTracks[i].charge());
                      tree_track_bkg_source.push_back(3);
                      //we have access to the same amount of information for the covnersion as the tracks in pfcandidate
                      //aka TrackBase.h
                    }
                tree_Yc_mass.push_back((*PhotonConversion)[j].pairInvariantMass());
              }
          }
        else tree_Yc_mass.push_back(-1);
      }
    }
  }


  //////////////////////////////////
  //////////////////////////////////
  /////////   Simulation   /////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_GenPVx = -1.;
  tree_GenPVy = -1.;
  tree_GenPVz = -20.;

  int nLLP = 0;
  int nllp = 0;
  nBC = 0; 
  tree_nFromC = 0; 
  tree_nFromB = 0;
  bool showlog=false;
  // Gen Information  for event axis //
  float  Gen_neu1_eta=-10, Gen_neu1_phi=-10;
  float  Gen_neu2_eta=-10, Gen_neu2_phi=-10;
  int  neu[2];
  int  nneu = 0;
  TLorentzVector vneu[2];
  
  float dRneuneu = 0.;
  
  for (int k=0; k<2; k++) {
    neu[k] = -1;
  }

  if ( !runOnData_ ) {
  
    // cout << endl; cout << endl; cout << endl;
    int genParticle_idx=0;
    for (size_t i=0; i<pruned->size(); i++)
    {
      const GenParticle & genIt = (*pruned)[i];
      const Candidate * mom   = genIt.mother();
      unsigned int nDaughters = genIt.numberOfDaughters();
      genParticle_idx++;

      int ID = abs(genIt.pdgId());
      float Gen_pt  = genIt.pt();
      float Gen_eta = genIt.eta();
      float Gen_phi = genIt.phi();
      float Gen_m   = genIt.mass();

      // primary incoming partons
      if ( genIt.status() == 21  ) {
	tree_GenPVx = genIt.vx();
	tree_GenPVy = genIt.vy();
	tree_GenPVz = genIt.vz();
// 	if ( i > 3 ) {
// 	  cout << " !!! incoming parton to be checked !!! " << endl;
//           cout << i << " status " << genIt.status() << " pt eta phi id " 
// 	       << Gen_pt << " " << Gen_eta << " " << Gen_phi << " " << genIt.pdgId() 
//                << " x y z " << genIt.vx() << " " << genIt.vy() << " " << genIt.vz() << " " << endl; 
// 	}
      }
      
      // neutralino from smuon
      if ( ID == 1000023 && abs(mom->pdgId()) == 1000013 ) {
	nLLP++;
        // cout << " neutralino" << nLLP << " pt eta phi " << Gen_pt << " " << Gen_eta << " " << Gen_phi << endl;
	if ( nLLP == 1 ) {
	  LLP1_pt  = Gen_pt;
	  LLP1_eta = Gen_eta;
	  LLP1_phi = Gen_phi;
	}
	if ( nLLP == 2 ) {
	  LLP2_pt  = Gen_pt;
	  LLP2_eta = Gen_eta;
	  LLP2_phi = Gen_phi;
	}
	if ( neu[0] < 0 ) {
	  neu[0] = genParticle_idx;
	  vneu[0].SetPtEtaPhiM( Gen_pt, Gen_eta, Gen_phi, Gen_m );
	  Gen_neu1_eta = Gen_eta;
	  Gen_neu1_phi = Gen_phi;
	}
	else if ( neu[1] < 0 ) {
	  neu[1] = genParticle_idx;
	  vneu[1].SetPtEtaPhiM( Gen_pt, Gen_eta, Gen_phi, Gen_m );
	  Gen_neu2_eta = Gen_eta;
	  Gen_neu2_phi = Gen_phi;
	}
	nneu++;
      }
      
      if ( nneu == 2 ) {
	dRneuneu = Deltar( Gen_neu1_eta, Gen_neu1_phi, Gen_neu2_eta, Gen_neu2_phi );
        tree_genAxis_dRneuneu.push_back(dRneuneu);
      }
      
      // quarks from neutralino
      if ( ID >= 1 && ID <= 6 && abs(mom->pdgId()) == 1000023 ) {
	if ( nllp >= 2 ) {
	  float dV1 = (genIt.vx() - LLP1_x)*(genIt.vx() - LLP1_x)
	            + (genIt.vy() - LLP1_y)*(genIt.vy() - LLP1_y)
	            + (genIt.vz() - LLP1_z)*(genIt.vz() - LLP1_z); // dV1 is equal to dV from nllp==1
	  float dV2 = (genIt.vx() - LLP2_x)*(genIt.vx() - LLP2_x)
	            + (genIt.vy() - LLP2_y)*(genIt.vy() - LLP2_y)
	            + (genIt.vz() - LLP2_z)*(genIt.vz() - LLP2_z);
	  if ( dV1 > 0.01 && dV2 > 0.01 ) nllp++; // should be == 2, so just to check : dV2 is always equal to 0 here
	}
	if ( nllp == 1 ) {
	  float dV = (genIt.vx() - LLP1_x)*(genIt.vx() - LLP1_x)
	           + (genIt.vy() - LLP1_y)*(genIt.vy() - LLP1_y)
	           + (genIt.vz() - LLP1_z)*(genIt.vz() - LLP1_z);
	  if ( dV > 0.01 ) {
	    nllp = 2;
	    LLP2_x = genIt.vx();
	    LLP2_y = genIt.vy();
	    LLP2_z = genIt.vz();
	    LLP2_dist = TMath::Sqrt( (LLP2_x - tree_GenPVx)*(LLP2_x - tree_GenPVx) 
				   + (LLP2_y - tree_GenPVy)*(LLP2_y - tree_GenPVy) 
				   + (LLP2_z - tree_GenPVz)*(LLP2_z - tree_GenPVz) ); 
	  }
	}
	if ( nllp == 0 ) {
	  nllp = 1;
	  LLP1_x = genIt.vx();
	  LLP1_y = genIt.vy();
	  LLP1_z = genIt.vz();
	  LLP1_dist = TMath::Sqrt( (LLP1_x - tree_GenPVx)*(LLP1_x - tree_GenPVx) 
				 + (LLP1_y - tree_GenPVy)*(LLP1_y - tree_GenPVy) 
				 + (LLP1_z - tree_GenPVz)*(LLP1_z - tree_GenPVz) ); 
	}
        // cout << " quark " << genIt.pdgId() << " from " << mom->pdgId() 
        //      << " pt eta phi " << Gen_pt << " " << Gen_eta << " " << Gen_phi 
        //      << " x y z " << genIt.vx() << " " << genIt.vy() << " " << genIt.vz() 
        //      << endl;
      }
      
      // Final c Hadron and get all its final charged particles
      bool isFinalD = false;
      if ( (ID/100)%10 == 4 || (ID/1000)%10 == 4 ) {
        isFinalD = true;
        for (unsigned int d1=0; d1<nDaughters; d1++) {
          const Candidate* gen1 = genIt.daughter(d1);
          int ID1 = abs(gen1->pdgId());
          if ( (ID1/100)%10 == 4 || (ID1/1000)%10 == 4 ) isFinalD = false;
        }
      }
      if ( isFinalD && abs(genIt.eta()) < 4. ) {
        const Candidate * Ancestor = &genIt;
        for (size_t j=0; j<packed->size(); j++) 
        {
        if ( (*packed)[j].pt() < 0.9 || fabs((*packed)[j].eta()) > 3.0 || (*packed)[j].charge() == 0 ) continue;
          //get the pointer to the first survived ancestor of a given packed GenParticle in the prunedCollection
  	  const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
  	if ( !(motherInPrunedCollection != nullptr && isAncestor( Ancestor , motherInPrunedCollection)) ) continue;
          tree_nFromC++;
          tree_genFromC_pt.push_back(	 (*packed)[j].pt());
          tree_genFromC_eta.push_back(   (*packed)[j].eta());
          tree_genFromC_phi.push_back(   (*packed)[j].phi());
          tree_genFromC_charge.push_back((*packed)[j].charge());
          tree_genFromC_pdgId.push_back( (*packed)[j].pdgId());
          tree_genFromC_mother_pdgId.push_back( genIt.pdgId());
	  if ( nDaughters > 0 ) {
            const Candidate* gen2 = genIt.daughter(0);
            tree_genFromC_x.push_back(gen2->vx());
            tree_genFromC_y.push_back(gen2->vy());
            tree_genFromC_z.push_back(gen2->vz());
	  }
	  else { // never happens a priori
            tree_genFromC_x.push_back(-10);
            tree_genFromC_y.push_back(-10);
            tree_genFromC_z.push_back(-10);
	  }
        }
      } // final c hadron

      // Final b Hadron and get all its final charged particles
      bool isFinalB = false;
      if ( (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) {
        isFinalB = true;
        for (unsigned int d1=0; d1<nDaughters; d1++) {
          const Candidate* gen1 = genIt.daughter(d1);
          int ID1 = abs(gen1->pdgId());
          if ( (ID1/100)%10 == 5 || (ID1/1000)%10 == 5 ) isFinalB = false;
        }
      }
      if ( isFinalB && abs(genIt.eta()) < 4. ) {
        const Candidate * Ancestor = &genIt;
        for (size_t j=0; j<packed->size(); j++) 
        {
        if ( (*packed)[j].pt() < 0.9 || fabs((*packed)[j].eta()) > 3.0 || (*packed)[j].charge() == 0 ) continue;
          //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
  	  const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
  	if ( !(motherInPrunedCollection != nullptr && isAncestor( Ancestor , motherInPrunedCollection)) ) continue;
          tree_nFromB++;
          tree_genFromB_pt.push_back(	 (*packed)[j].pt());
          tree_genFromB_eta.push_back(   (*packed)[j].eta());
          tree_genFromB_phi.push_back(   (*packed)[j].phi());
          tree_genFromB_charge.push_back((*packed)[j].charge());
          tree_genFromB_pdgId.push_back( (*packed)[j].pdgId());
          tree_genFromB_mother_pdgId.push_back( genIt.pdgId());
	  if ( nDaughters > 0 ) {
            const Candidate* gen2 = genIt.daughter(0);
            tree_genFromB_x.push_back(gen2->vx());
            tree_genFromB_y.push_back(gen2->vy());
            tree_genFromB_z.push_back(gen2->vz());
	  }
	  else { // never happens a priori
            tree_genFromB_x.push_back(-10);
            tree_genFromB_y.push_back(-10);
            tree_genFromB_z.push_back(-10);
	  }
        }
      } // final b hadron
    
      float dV0 = (genIt.vx() - tree_GenPVx)*(genIt.vx() - tree_GenPVx)
        	+ (genIt.vy() - tree_GenPVy)*(genIt.vy() - tree_GenPVy)
        	+ (genIt.vz() - tree_GenPVz)*(genIt.vz() - tree_GenPVz);
      float dV1 = (genIt.vx() - LLP1_x)*(genIt.vx() - LLP1_x)
        	+ (genIt.vy() - LLP1_y)*(genIt.vy() - LLP1_y)
        	+ (genIt.vz() - LLP1_z)*(genIt.vz() - LLP1_z);
      float dV2 = (genIt.vx() - LLP2_x)*(genIt.vx() - LLP2_x)
        	+ (genIt.vy() - LLP2_y)*(genIt.vy() - LLP2_y)
        	+ (genIt.vz() - LLP2_z)*(genIt.vz() - LLP2_z);
      int fromLLP = -1;
      if      ( dV1 < dV2 && dV1 < 0.01 ) fromLLP = 1;
      else if ( dV2 < dV1 && dV2 < 0.01 ) fromLLP = 2;
      else if ( dV0 < 0.01 )		  fromLLP = 0;

    if ( genIt.pt() < 0.9 || fabs(genIt.eta()) > 4.0 ) continue;
      
      tree_genParticle_pt.push_back(        genIt.pt());
      tree_genParticle_eta.push_back(       genIt.eta());
      tree_genParticle_phi.push_back(       genIt.phi());
      tree_genParticle_charge.push_back(    genIt.charge());
      tree_genParticle_pdgId.push_back(     genIt.pdgId());
      tree_genParticle_mass.push_back(      genIt.mass());
      tree_genParticle_x.push_back(	    genIt.vx());
      tree_genParticle_y.push_back(	    genIt.vy());
      tree_genParticle_z.push_back(	    genIt.vz());
      tree_genParticle_statusCode.push_back(genIt.status());
      tree_genParticle_mother_pdgId.push_back( mom ? mom->pdgId() :  -10 );
      tree_genParticle_LLP.push_back(fromLLP);

    } // end loop on pruned genparticles

    tree_nLLP = nllp;
    // cout << endl;

    // second pass to recover the final particles from LLP decay
    int nLLPbis = 0;
    tree_ngenFromLLP = 0;

    for (size_t i=0; i<pruned->size(); i++) // loop on pruned genparticles
    {
      const GenParticle & genIt = (*pruned)[i];
      const Candidate * mom = genIt.mother();
      int pdgid     = genIt.pdgId();

    // neutralino
    if ( !(pdgid == 1000023 && abs(mom->pdgId()) == 1000013) ) continue;
      nLLPbis++;
      // if ( nLLPbis == 2 ) cout << endl;
      const Candidate * Neutralino = &genIt;
      for (size_t j=0; j<packed->size(); j++) // loop on packed genparticles
      {
      if ( (*packed)[j].pt() < 0.9 || fabs((*packed)[j].eta()) > 3.0 || (*packed)[j].charge() == 0 ) continue;
        //get the pointer to the first survived ancestor of a given packed GenParticle in the prunedCollection
        const Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
      if ( !(motherInPrunedCollection != nullptr && isAncestor( Neutralino , motherInPrunedCollection)) ) continue;
        tree_ngenFromLLP++;
        tree_genFromLLP_LLP.push_back(       nLLPbis);
        float pack_pt  = (*packed)[j].pt();
        float pack_eta = (*packed)[j].eta();
        float pack_phi = (*packed)[j].phi();
        float pack_pdgId = (*packed)[j].pdgId();

        tree_genFromLLP_pt.push_back(	     pack_pt);
        tree_genFromLLP_eta.push_back(       pack_eta);
        tree_genFromLLP_phi.push_back(       pack_phi);
        tree_genFromLLP_charge.push_back(    (*packed)[j].charge());
        tree_genFromLLP_pdgId.push_back(     pack_pdgId);
        tree_genFromLLP_mass.push_back(      (*packed)[j].mass());
        const Candidate * momj =	     (*packed)[j].mother(0);

        int mom_pdgid = -9999;
	float vx = -10., vy = -10., vz = -10.;
        if ( momj ) {
	  mom_pdgid = momj->pdgId();
	  vx = momj->vx();
	  vy = momj->vy();
	  vz = momj->vz();
          if ( momj->numberOfDaughters() > 0 ) { // always the case a priori
	    vx = momj->daughter(0)->vx();
	    vy = momj->daughter(0)->vy();
	    vz = momj->daughter(0)->vz();
	  }
	}
        tree_genFromLLP_mother_pdgId.push_back(mom_pdgid);
        tree_genFromLLP_x.push_back( vx );
        tree_genFromLLP_y.push_back( vy );
        tree_genFromLLP_z.push_back( vz );

        // match to final b hadron
	bool matchB = false;
	for (int k = 0; k < tree_nFromB; k++)
	{
        if ( pack_pdgId != tree_genFromB_pdgId[k] ) continue;
	  float dpt  = abs( pack_pt / tree_genFromB_pt[k] - 1. );
	  float deta = abs( pack_eta - tree_genFromB_eta[k] );
  	  float dphi = abs( Deltaphi( pack_phi, tree_genFromB_phi[k] ) );
          if ( abs(deta) < 0.01 && abs(dphi) < 0.01 && abs(dpt) < 0.01 ) {
	    matchB = true;
	    break;
	  }
	}
        tree_genFromLLP_isFromB.push_back(matchB);

        // match to final c hadron
	bool matchC = false;
	for (int k = 0; k < tree_nFromC; k++)
	{
        if ( pack_pdgId != tree_genFromC_pdgId[k] ) continue;
	  float dpt  = abs( pack_pt / tree_genFromC_pt[k] - 1. );
	  float deta = abs( pack_eta - tree_genFromC_eta[k] );
  	  float dphi = abs( Deltaphi( pack_phi, tree_genFromC_phi[k] ) );
          if ( abs(deta) < 0.01 && abs(dphi) < 0.01 && abs(dpt) < 0.01 ) {
	    matchC = true;
	    break;
	  }
	}
        tree_genFromLLP_isFromC.push_back(matchC);
        // cout << " gentk " << pack_pdgId << " from " << mom_pdgid 
        //      << " LLP " << nLLPbis << " BC " << matchB << matchC
        //      << " pt eta phi " << pack_pt << " " << pack_eta << " " << pack_phi 
        //      << " x y z " << vx << " " << vy << " " << vz 
        //      << endl;

      } // end loop on packed genparticles
      if ( nLLPbis == 2 ) break;

    } // end loop on pruned genparticles

    // packed genparticles (final particles)
    tree_ngenPackPart = 0;
    for (size_t i=0; i<packed->size(); i++) 
    {
    if ( (*packed)[i].pt() < 0.9 || fabs((*packed)[i].eta()) > 3.0 || (*packed)[i].charge() == 0 ) continue;
      const Candidate * mom = (*packed)[i].mother(0);
      tree_genPackPart_pt.push_back(        (*packed)[i].pt());
      tree_genPackPart_eta.push_back(       (*packed)[i].eta());
      tree_genPackPart_phi.push_back(       (*packed)[i].phi());
      tree_genPackPart_charge.push_back(    (*packed)[i].charge());
      tree_genPackPart_pdgId.push_back(     (*packed)[i].pdgId());
      tree_genPackPart_mass.push_back(      (*packed)[i].mass());
      tree_genPackPart_mother_pdgId.push_back( mom ? mom->pdgId() :  -10 );

      float vx = -10., vy = -10., vz = -10.;
      if ( mom != nullptr ) {
        if ( mom->vx() == 0 && mom->vy() == 0 && mom->vz() == 0 ) {
	  vx = tree_GenPVx;
	  vy = tree_GenPVy;
	  vz = tree_GenPVz;
	}
	else {
	  vx = mom->vx();
	  vy = mom->vy();
	  vz = mom->vz();
	}
      }
      tree_genPackPart_x.push_back( vx );
      tree_genPackPart_y.push_back( vy );
      tree_genPackPart_z.push_back( vz );

      // match to final b hadron
      bool matchB = false;
      for (int k = 0; k < tree_nFromB; k++)
      {
      if ( (*packed)[i].pdgId() != tree_genFromB_pdgId[k] ) continue;
        float dpt  = abs( (*packed)[i].pt() / tree_genFromB_pt[k] - 1. );
        float deta = abs( (*packed)[i].eta() - tree_genFromB_eta[k] );
        float dphi = abs( Deltaphi( (*packed)[i].phi(), tree_genFromB_phi[k] ) );
        if ( abs(deta) < 0.01 && abs(dphi) < 0.01 && abs(dpt) < 0.01 ) {
          matchB = true;
          break;
        }
      }
      tree_genPackPart_isFromB.push_back(matchB);

      // match to final c hadron
      bool matchC = false;
      for (int k = 0; k < tree_nFromC; k++)
      {
      if ( (*packed)[i].pdgId() != tree_genFromC_pdgId[k] ) continue;
        float dpt  = abs( (*packed)[i].pt() / tree_genFromC_pt[k] - 1. );
        float deta = abs( (*packed)[i].eta() - tree_genFromC_eta[k] );
        float dphi = abs( Deltaphi( (*packed)[i].phi(), tree_genFromC_phi[k] ) );
        if ( abs(deta) < 0.01 && abs(dphi) < 0.01 && abs(dpt) < 0.01 ) {
          matchC = true;
          break;
        }
      }
      tree_genPackPart_isFromC.push_back(matchC);
      tree_ngenPackPart++;
    }

    // gen jets
    for (auto const & genJet : *genJets)
    {
    if ( genJet.pt() < 20. ) continue;
      tree_genJet_pt.push_back(genJet.pt());
      tree_genJet_eta.push_back(genJet.eta());
      tree_genJet_phi.push_back(genJet.phi());
      tree_genJet_mass.push_back(genJet.mass());
      tree_genJet_energy.push_back(genJet.energy());
    }
    
  } // endif simulation


  //////////////////////////////////
  //////////////////////////////////
  ///////////   MET   //////////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_PFMet_et  = -10.;
  tree_PFMet_phi = -10.;
  tree_PFMet_sig = -10.;
  if ( PFMETs->size() > 0 ) {
    const pat::MET &themet = PFMETs->front();
    tree_PFMet_et  = themet.et();
    tree_PFMet_phi = themet.phi();
    tree_PFMet_sig = themet.significance();
  }

  //////////////////////////////////
  //////////////////////////////////
  ///////////	Jets   /////////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_njet = 0;
  float HT_val = 0;
  float jet_pt_min = 20.;
  for (const pat::Jet &jet : *jets) {
  if ( jet.pt() < jet_pt_min ) continue;
    tree_jet_E.push_back(jet.energy());
    tree_jet_pt.push_back(jet.pt());
    tree_jet_eta.push_back(jet.eta());
    tree_jet_phi.push_back(jet.phi());
    tree_njet++;
    if ( abs(jet.eta()) < 2.4 ) HT_val += jet.pt(); // used in HT filter !
  }
    tree_HT = HT_val;

  //////////////////////////////////
  //////////////////////////////////
  ////////   Electrons   ///////////
  //////////////////////////////////
  //////////////////////////////////
    //One could consider the emu channel, or even ee channel
  // float Eiso=0;  
  for (const pat::Electron &el: *electrons)
  {
  if ( el.pt() < 5. ) continue;
    tree_electron_pt.push_back(     el.pt());
    tree_electron_eta.push_back(    el.eta());
    tree_electron_phi.push_back(    el.phi());
    tree_electron_x.push_back(      el.vx());
    tree_electron_y.push_back(      el.vy());
    tree_electron_z.push_back(      el.vz());
    tree_electron_energy.push_back( el.energy());
    tree_electron_charge.push_back(el.charge());
    tree_electron_isoR4.push_back(el.trackIso());//returns the value of the summed track pt in a cone of deltaR<0.4
    //     float ecalIso() const { return dr04EcalRecHitSumEt(); }
    // /// Overload of pat::Lepton::hcalIso(); returns the value of the summed Et of all caloTowers in the hcal in a cone of deltaR<0.4
    // float hcalIso() const { return dr04HcalTowerSumEt(); }
    // /// Overload of pat::Lepton::caloIso(); returns the sum of ecalIso() and hcalIso
    // float caloIso() const { return ecalIso() + hcalIso(); }
    // std::cout<<"Electron isolation : "<<Eiso<<std::endl;
  }

  //////////////////////////////////
  //////////////////////////////////
  ///////////   Muons   ////////////
  //////////////////////////////////
  //////////////////////////////////
  
  int nmu = 0;
   // float isoR3=0;
  for (const pat::Muon &mu : *muons)
  {
  if ( mu.pt() < 3. ) continue;
  if ( abs(mu.eta()) > 2.4 ) continue;  // muon acceptance
    tree_muon_pt.push_back(       mu.pt());
    tree_muon_eta.push_back(      mu.eta());
    tree_muon_phi.push_back(      mu.phi());
    tree_muon_x.push_back(        mu.vx());
    tree_muon_y.push_back(        mu.vy());
    tree_muon_z.push_back(        mu.vz());
    tree_muon_energy.push_back(   mu.energy());
    tree_muon_dxy.push_back(	  mu.muonBestTrack()->dxy(PV.position()));
    tree_muon_dxyError.push_back( mu.muonBestTrack()->dxyError());
    tree_muon_dz.push_back(       mu.muonBestTrack()->dz(PV.position()));
    tree_muon_dzError.push_back(  mu.muonBestTrack()->dzError());
    tree_muon_charge.push_back(   mu.charge());
    tree_muon_isLoose.push_back(  mu.isLooseMuon());
    tree_muon_isTight.push_back(  mu.isTightMuon(PV));
    tree_muon_isGlobal.push_back( mu.isGlobalMuon());
    tree_muon_isoR3.push_back(mu.trackIso());
        // isoR3 = mu.trackIso();//the summed track pt in a cone of deltaR<0.3
    //     /// Overload of pat::Lepton::trackIso(); returns the value of
    // /// the summed Et of all recHits in the ecal in a cone of
    // /// deltaR<0.3
    // float ecalIso() const { return isolationR03().emEt; }
    // /// Overload of pat::Lepton::trackIso(); returns the value of
    // /// the summed Et of all caloTowers in the hcal in a cone of
    // /// deltaR<0.4
    // float hcalIso() const { return isolationR03().hadEt; }
    // /// Overload of pat::Lepton::trackIso(); returns the sum of
    // /// ecalIso() and hcalIso
    // float caloIso() const { return ecalIso() + hcalIso(); }
    // std::cout<<"Muon isolation : "<<isoR3<<std::endl;
    nmu++;
  }
    
  int imu1 = -1, imu2 = -1;
  float mupt1, mueta1, muphi1, mupt2, mueta2, muphi2;
  float mu_mass = 0.1057;
  TLorentzVector v1, v2, v;
  tree_Mmumu = 0.;
  
  if ( nmu >= 2 ) {
    for (int mu = 0; mu < nmu-1; mu++)
    { 
    if ( !tree_muon_isGlobal[mu] ) continue;
      mupt1  = tree_muon_pt[mu];
    if ( mupt1 < 10. ) continue; // Zmu filter
    if ( abs(tree_muon_dxy[mu]) > 0.1 || abs(tree_muon_dz[mu]) > 0.2 ) continue; // muons closed to PV
      mueta1 = tree_muon_eta[mu];
      muphi1 = tree_muon_phi[mu];
      v1.SetPtEtaPhiM(mupt1,mueta1,muphi1,mu_mass);
      for ( int mu2=mu+1; mu2<nmu; mu2++) 
      {	    
      if ( !tree_muon_isGlobal[mu2] ) continue;
      if ( tree_muon_charge[mu] == tree_muon_charge[mu2] ) continue;
      if ( abs(tree_muon_dxy[mu2]) > 0.1 || abs(tree_muon_dz[mu2]) > 0.2 ) continue;
        mupt2  = tree_muon_pt[mu2];
      if ( mupt2 < 10. ) continue;
      if ( mupt1 < 25. && mupt2 < 25. ) continue; // Zmu Filter
        mueta2 = tree_muon_eta[mu2];
        muphi2 = tree_muon_phi[mu2];
        v2.SetPtEtaPhiM(mupt2,mueta2,muphi2,mu_mass);
        v = v1 + v2;
        if ( v.Mag() > tree_Mmumu )
        { // Mag pour masse invariante (magnitude)
          tree_Mmumu = v.Mag();
          imu1 = mu;
          imu2 = mu2;
        }
      }
    } // end loop on muons
  }

  if ( imu1 >= 0 && imu2 >= 0 && tree_muon_pt[imu2] > tree_muon_pt[imu1] ) {
    int imu0 = imu2;
    imu2 = imu1; // muons reco with imu1 having the highest pt
    imu1 = imu0;
  }

  //////////////////////////////////
  //////////////////////////////////
  ////////   FILTER CHECK  /////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_NbrOfZCand = 0;
  tree_Filter = false;
  tree_nTracks = 0;
  tree_nLostTracks = 0;

  if ( tree_Mmumu > 60. )                  tree_NbrOfZCand = 1;
//   if ( tree_Mmumu > 60. && HT_val > 180. ) tree_Filter = true;
//$$
  if ( (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v || HLT_IsoMu24_v) 
       && tree_Mmumu > 10. ) tree_Filter = true;
//$$
 
  edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder); // Asking for reco collection of PV..
  vector<reco::TransientTrack> BestTracks;
  std::vector<std::pair<uint16_t,float> > Players;
  int count =0;
  std::map<size_t , int > trackToAK4SlimmedJetMap;

 if ( tree_Filter && ActivateTrigger) {


  //////////////////////////////////
  //////////////////////////////////
  //////////   Tracks   ////////////
  //////////////////////////////////
  //////////////////////////////////
  
    float pt_Cut = 1.;//1. for MC Signal
    float NChi2_Cut = 5.; // 5. for MC Signal
    float drSig_Cut = 5.;// 5. for MC Signal

    // track variables declaration //
    float tk_nHit   ;
    float tk_charge ;
    float tk_pt ;
    float tk_eta ;
    float tk_phi ;
    float tk_NChi2 ;
    float tk_drSig = -1;
    float tk_dxy ;
    float tk_dxyError ;
    float tk_dz ;
    float tk_dzError;
    float tk_vx ;
    float tk_vy ;
    float tk_vz ;
    float tk_px ;
    float tk_py ;
    float tk_pz ;
    HitPattern tk_HitPattern;

    vector <pat::PackedCandidateRef> MINIgeneralTracks;
    for (unsigned int ipc = 0; ipc < pc->size(); ipc++) 
    {
      MINIgeneralTracks.push_back(pat::PackedCandidateRef(pcs, ipc));
    }
    for (unsigned int k=0; k<lostpc->size();k++)
    {
      MINIgeneralTracks.push_back(pat::PackedCandidateRef(lostpcs, k));
    }

    // loop on all packedPFCandidates + lostTracks

    // from /PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc
    // and
    // from /DQM/TrackingMonitor/src/PackedCandidateTrackValidator.cc

    for (unsigned int ipc = 0; ipc < pc->size()+lostpc->size(); ipc++) {
      pat::PackedCandidateRef pcref = MINIgeneralTracks[ipc];
      const reco::Track *trackPcPtr = pcref->bestTrack(); // const
//       if ( !trackPcPtr) tree_passesTrkPtr.push_back(0);
    if ( !trackPcPtr) continue;
//       tree_passesTrkPtr.push_back(1);
      reco::Track tk;
      
        tk = *trackPcPtr;
        tk_nHit   = tk.hitPattern().numberOfValidHits();
        tk_charge = tk.charge();
        tk_pt  = tk.pt();
        tk_eta = tk.eta();
        tk_phi = tk.phi();
        tk_NChi2 = tk.normalizedChi2();
        tk_dxy = tk.dxy(PV.position());
        tk_dxyError = tk.dxyError();
        tk_dz = tk.dz(PV.position());
        tk_dzError = tk.dzError();
        if ( tk_dxyError > 0 ) 
          tk_drSig = abs(tk_dxy) / tk_dxyError; // from Paul
        tk_vx = tk.vx();
        tk_vy = tk.vy();
        tk_vz = tk.vz();
        tk_px = tk.px();
        tk_py = tk.py();
        tk_pz = tk.pz();
        tk_HitPattern = tk.hitPattern();

    if ( tk_nHit == 0 ) continue;
    if ( tk_charge == 0 ) continue;
//$$
    if ( !(tk_pt > pt_Cut && tk_NChi2 < NChi2_Cut && tk_drSig > drSig_Cut) ) continue; // preselection
//$$

// Tracks from V0Candidates and  photon conversion
  bool matchTObkg = false;
  for (unsigned int i = 0; i< tree_track_bkg_pt.size() ; i++ ) 
    {
      float dpt = tree_track_bkg_pt[i]-tk_pt;
      float deta = tree_track_bkg_eta[i]-tk_eta;
      float dphi = tree_track_bkg_phi[i]-tk_phi;
      float bkg_charge = tree_track_bkg_charge[i];
      tree_track_dpt.push_back(dpt);
      tree_track_deta.push_back(deta);
      tree_track_dphi.push_back(dphi);
      if ( abs(deta) < 0.01 && abs(dphi) < 0.01 && abs(dpt) < 0.01 && bkg_charge == tk_charge && TrackMatchingToV0)   { matchTObkg = true;}
    }
    tree_track_bkg.push_back(matchTObkg);
    if (matchTObkg)continue;

      tree_nTracks++;
      tree_track_ipc.push_back(ipc);
      if ( ipc < pc->size() ) {
        tree_track_lost.push_back(0);
      }
      else {
        tree_track_lost.push_back(1);
        tree_nLostTracks++;
      }  
      tree_track_pt.push_back           (tk_pt);
      tree_track_eta.push_back          (tk_eta);
      tree_track_phi.push_back          (tk_phi);
      tree_track_charge.push_back       (tk_charge);
      tree_track_NChi2.push_back        (tk_NChi2);
      tree_track_dxy.push_back          (tk_dxy);
      tree_track_dxyError.push_back     (tk_dxyError);
      tree_track_drSig.push_back        (tk_drSig); 
      tree_track_dz.push_back           (tk_dz);
      tree_track_dzError.push_back      (tk_dzError);
//       tree_track_algo.push_back         (tk_temp.algo());
      tree_track_isHighPurity.push_back (static_cast<int>(tk.quality(reco::TrackBase::highPurity)));
      tree_track_nHit.push_back         (tk_nHit);
      tree_track_nHitPixel.push_back    (tk_HitPattern.numberOfValidPixelHits());
      tree_track_nHitTIB.push_back      (tk_HitPattern.numberOfValidStripTIBHits());
      tree_track_nHitTID.push_back      (tk_HitPattern.numberOfValidStripTIDHits());
      tree_track_nHitTOB.push_back      (tk_HitPattern.numberOfValidStripTOBHits());
      tree_track_nHitTEC.push_back      (tk_HitPattern.numberOfValidStripTECHits());
      tree_track_nHitPXB.push_back      (tk_HitPattern.numberOfValidPixelBarrelHits());
      tree_track_nHitPXF.push_back      (tk_HitPattern.numberOfValidPixelEndcapHits());

      int hitPixelLayer = 0;
      if ( tk_HitPattern.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) ) hitPixelLayer += 1;
      if ( tk_HitPattern.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2) ) hitPixelLayer += 10;
      if ( tk_HitPattern.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) ) hitPixelLayer += 100;
      if ( tk_HitPattern.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 4) ) hitPixelLayer += 1000;
      if ( tk_HitPattern.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1) ) hitPixelLayer += 2;
      if ( tk_HitPattern.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2) ) hitPixelLayer += 20;
      if ( tk_HitPattern.hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 3) ) hitPixelLayer += 200;
      tree_track_isHitPixel.push_back( hitPixelLayer );

      tree_track_nLayers.push_back      (tk_HitPattern.trackerLayersWithMeasurement());
      tree_track_nLayersPixel.push_back (tk_HitPattern.pixelLayersWithMeasurement());
      tree_track_x.push_back            (tk_vx);
      tree_track_y.push_back            (tk_vy);
      tree_track_z.push_back            (tk_vz);

                  //----------------MINIAOD_Firsthit----------//
                  //-----------------IMPORTANT----------------//
                  // TSOS is said to be better for the -------//
                  // propagators (see Propagator.h)...--------//
                  // ../interface/PropaHitPattern.h           //
                  //------------------------------------------//
      //-----hitpattern -> Database ---/
      const HitPattern hp = tk_HitPattern;
      uint16_t firsthit = hp.getHitPattern(HitPattern::HitCategory::TRACK_HITS,0);
//$$      Approximation for the lostTrack since the hitpattern information is not available (only 1160)      
      if ( ipc >= pc->size() ) {
        if ( abs(tk_eta) < 1. ) firsthit = 1184; // PIXBL4 in barrel
        else                    firsthit = 1296; // PIXFD2 in forward
      }      
//$$      
      tree_track_firstHit.push_back(firsthit);

      //---Creating State to propagate from  TT---//
      BestTracks.push_back(theTransientTrackBuilder->build(tk));
      const MagneticField* B = BestTracks[count].field(); // 3.8T
      reco::TransientTrack TT (tk,BestTracks[count].field());
      // const FreeTrajectoryState Freetraj = TT.initialFreeState(); // Propagator in the barrel can also use FTS (WARNING: the so-called reference point (where the propagation starts might be different from the first vtx, a check should be done))
      GlobalPoint vert (tk_vx,tk_vy,tk_vz); // Point where the propagation will start (Reference Point)
      const TrajectoryStateOnSurface Surtraj = TT.stateOnSurface(vert); // TSOS of this point
      AnalyticalPropagator* Prop = new AnalyticalPropagator(B); // Propagator that will be used for barrel, crashes in the disks when using Plane
      Basic3DVector<float> P3D2(tk_vx,tk_vy,tk_vz);  // global frame
      Basic3DVector<float> B3DV (tk_px,tk_py,tk_pz); // global frame 
      float Eta = tk_eta;
      float Phi = tk_phi;
      float vz  = tk_vz;
      //------Propagation with new interface --> See ../interface/PropaHitPattern.h-----//
      PropaHitPattern* PHP = new PropaHitPattern();
      std::pair<int,GloballyPositioned<float>::PositionType> FHPosition = PHP->Main(firsthit,Prop,Surtraj,Eta,Phi,vz,P3D2,B3DV);
   //returns 0 if in Barrel, 1 if disks with the position of the firsthit. Different porpagators are used between 
      // barrel (StraightLinePlaneCrossing/geometry if propagator does not work)
      // and disks (HelixPlaneCrossing) 

      float xFirst = FHPosition.second.x();
      float yFirst = FHPosition.second.y();
      float zFirst = FHPosition.second.z();
      tree_track_firstHit_x.push_back(xFirst);
      tree_track_firstHit_y.push_back(yFirst);
      tree_track_firstHit_z.push_back(zFirst);
      tree_track_region.push_back(FHPosition.first);
      count++;
      //-----------------------END OF MINIAOD firsthit-----------------------//

      // track association to jet
      int iJet = 0;
      bool matchTOjet = false;
      for (const pat::Jet &jet : *jets) {
      if ( jet.pt() < jet_pt_min ) continue;
        float dR = Deltar( jet.eta(), jet.phi(), tk_eta, tk_phi );
        if ( dR < 0.4 ) {
          matchTOjet = true;
          break;
        }
        else iJet++;
      }
      if ( matchTOjet ) tree_track_iJet.push_back (iJet);
      else              tree_track_iJet.push_back (-1);

      // match to gen particle from LLP decay
      int      kmatch = -1;
      float    dFirstGenMin = 1000000.;
      int      track_sim_LLP = -1;
      bool     track_sim_isFromB = 0;
      bool     track_sim_isFromC = 0;
      float    track_sim_pt = 0;
      float    track_sim_eta = 0;
      float    track_sim_phi = 0;
      int      track_sim_charge = 0;
      int      track_sim_pdgId = 0;
      float    track_sim_mass = 0;
      float    track_sim_x = 0;
      float    track_sim_y = 0;
      float    track_sim_z = 0;
      float    track_sim_LLP_r = 0.;
      float    track_sim_LLP_z = 0.;

      float xPV = tree_GenPVx;
      float yPV = tree_GenPVy;
      for (int k = 0; k < tree_ngenFromLLP; k++) // loop on final gen part from LLP
      {
      if ( tk_charge != tree_genFromLLP_charge[k] ) continue;

        float qGen   = tree_genFromLLP_charge[k];
        float ptGen  = tree_genFromLLP_pt[k];
        float etaGen = tree_genFromLLP_eta[k];
        float phiGen = tree_genFromLLP_phi[k]; // given at production point
        float xGen   = tree_genFromLLP_x[k];
        float yGen   = tree_genFromLLP_y[k];
        float zGen   = tree_genFromLLP_z[k];
        float qR = qGen * ptGen * 100 / 0.3 / 3.8;

        float dpt  = (tk_pt - ptGen) / tk_pt;
        float deta = tk_eta - etaGen;
	
	// compute phi0 at dca(PV) for the gen particle (instead of production point)
        float sin0 = qR * sin( phiGen ) + (xGen - xPV);
        float cos0 = qR * cos( phiGen ) - (yGen - yPV);
        float phi0 = TMath::ATan2( sin0, cos0 ); // but note that it can be wrong by +_pi ! 
        float dphi = tk_phi - phi0;
        if      ( dphi < -3.14159 / 2. ) dphi += 3.14159;
        else if ( dphi >  3.14159 / 2. ) dphi -= 3.14159;

        // resolution depend on the number of hits... (here select 97% of signal tracks)
        bool matchTOgen = false;
	if ( tk_nHit <= 10 ) {
          if ( abs(dpt) < 0.70 && abs(deta) < 0.30 && abs(dphi) < 0.08 ) matchTOgen = true; 
        }
        else if ( tk_nHit <= 13 ) {
          if ( abs(dpt) < 0.20 && abs(deta) < 0.12 && abs(dphi) < 0.05 ) matchTOgen = true; 
        }
        else if ( tk_nHit <= 17 ) {
          if ( abs(dpt) < 0.08 && abs(deta) < 0.04 && abs(dphi) < 0.03 ) matchTOgen = true; 
        }
        else {
          if ( abs(dpt) < 0.07 && abs(deta) < 0.02 && abs(dphi) < 0.02 ) matchTOgen = true; 
        }

	if ( matchTOgen ) {
	  float dFirstGen = (xFirst-xGen)*(xFirst-xGen) + (yFirst-yGen)*(yFirst-yGen) + (zFirst-zGen)*(zFirst-zGen);
	  if ( dFirstGen < dFirstGenMin ) {
	    kmatch = k;
	    dFirstGenMin = dFirstGen;
	  }
	}
      } // end loop on final gen part from LLP

      if ( kmatch >= 0 ) {
        track_sim_LLP =     tree_genFromLLP_LLP[kmatch];
        track_sim_isFromB = tree_genFromLLP_isFromB[kmatch];
        track_sim_isFromC = tree_genFromLLP_isFromC[kmatch];
        track_sim_pt  =     tree_genFromLLP_pt[kmatch];
        track_sim_eta =     tree_genFromLLP_eta[kmatch];
        track_sim_phi =     tree_genFromLLP_phi[kmatch];
        track_sim_charge =  tree_genFromLLP_charge[kmatch];
        track_sim_pdgId =   tree_genFromLLP_pdgId[kmatch];
        track_sim_mass =    tree_genFromLLP_mass[kmatch];
        track_sim_x =	    tree_genFromLLP_x[kmatch];
        track_sim_y =	    tree_genFromLLP_y[kmatch];
        track_sim_z =	    tree_genFromLLP_z[kmatch];

        float dSign = 1.;
        if ( xFirst*track_sim_x+yFirst*track_sim_y+zFirst*track_sim_z < 0. ) dSign = -1.;
	dFirstGenMin = TMath::Sqrt(dFirstGenMin)*dSign;
        if ( track_sim_LLP == 1 ) {
          track_sim_LLP_r = TMath::Sqrt( LLP1_x*LLP1_x + LLP1_y*LLP1_y );
          track_sim_LLP_z = abs(LLP1_z);
        }
        if ( track_sim_LLP == 2 ) {
          track_sim_LLP_r = TMath::Sqrt( LLP2_x*LLP2_x + LLP2_y*LLP2_y );
          track_sim_LLP_z = abs(LLP2_z);
        }
      }

      else {
        int omatch = -1;
        for (int k = 0; k < tree_ngenPackPart; k++) // loop on the other final gen particles 
        {
        if ( tk_charge != tree_genPackPart_charge[k] ) continue;
          float ptGen  = tree_genPackPart_pt[k];
          float etaGen = tree_genPackPart_eta[k];
          float phiGen = tree_genPackPart_phi[k]; // given at production point
          float dpt  = (tk_pt - ptGen) / tk_pt;
          float deta = tk_eta - etaGen;
  	  float dphi = Deltaphi( tk_phi, phiGen );
          // resolution depend on the number of hits... (here select 97% of signal tracks)
          bool matchTOgen = false;
	  if ( tk_nHit <= 10 ) {
            if ( abs(dpt) < 0.70 && abs(deta) < 0.30 && abs(dphi) < 0.08 ) matchTOgen = true; 
          }
          else if ( tk_nHit <= 13 ) {
            if ( abs(dpt) < 0.20 && abs(deta) < 0.12 && abs(dphi) < 0.05 ) matchTOgen = true; 
          }
          else if ( tk_nHit <= 17 ) {
            if ( abs(dpt) < 0.08 && abs(deta) < 0.04 && abs(dphi) < 0.03 ) matchTOgen = true; 
          }
          else {
            if ( abs(dpt) < 0.07 && abs(deta) < 0.02 && abs(dphi) < 0.02 ) matchTOgen = true; 
          }
	  if ( matchTOgen ) omatch = k;
        } // end loop on the other final gen particles

        if ( omatch >= 0 ) {
          track_sim_pt  =     tree_genPackPart_pt[omatch];
          track_sim_eta =     tree_genPackPart_eta[omatch];
          track_sim_phi =     tree_genPackPart_phi[omatch];
          track_sim_charge =  tree_genPackPart_charge[omatch];
          track_sim_pdgId =   tree_genPackPart_pdgId[omatch];
          track_sim_mass =    tree_genPackPart_mass[omatch];
          track_sim_x =	      tree_genPackPart_x[omatch];
          track_sim_y =	      tree_genPackPart_y[omatch];
          track_sim_z =	      tree_genPackPart_z[omatch];
          track_sim_isFromB = tree_genPackPart_isFromB[omatch];
          track_sim_isFromC = tree_genPackPart_isFromC[omatch];
        }
      }

      tree_track_sim_LLP.push_back(	  track_sim_LLP );
      tree_track_sim_isFromB.push_back(   track_sim_isFromB );
      tree_track_sim_isFromC.push_back(   track_sim_isFromC );
      tree_track_sim_pt.push_back(	  track_sim_pt );
      tree_track_sim_eta.push_back(	  track_sim_eta );
      tree_track_sim_phi.push_back(	  track_sim_phi );
      tree_track_sim_charge.push_back(    track_sim_charge );
      tree_track_sim_pdgId.push_back(	  track_sim_pdgId );
      tree_track_sim_mass.push_back(	  track_sim_mass );
      tree_track_sim_x.push_back(	  track_sim_x );
      tree_track_sim_y.push_back(	  track_sim_y );
      tree_track_sim_z.push_back(	  track_sim_z );
      tree_track_sim_dFirstGen.push_back( dFirstGenMin );
      tree_track_sim_LLP_r.push_back(     track_sim_LLP_r );
      tree_track_sim_LLP_z.push_back(     track_sim_LLP_z );

    } // end loop on all track candidates


    /////////////////////////////////////////////////////////
    //-------------------------------------------------------
    // Jets for event axes                                 
    //-------------------------------------------------------
    /////////////////////////////////////////////////////////

    int njet = 0, njet1 = 0, njet2 = 0;
    bool isjet[99], isjet1[99], isjet2[99];
    TLorentzVector vaxis1, vaxis2, vjet[99];
    float PtMin = 20;   // (GeV) minimum jet pt is optimum
    float EtaMax = 10.; // no cut on eta is optimum
    int jetidx = 0; // : May be in the loop/ not sure it changes anything
    for (const pat::Jet &jet : *jets)    // Loop on jet
    {
      float jet_pt  = jet.pt();
      float jet_eta = jet.eta();
      float jet_phi = jet.phi();
      isjet[jetidx]  = false;
      isjet1[jetidx] = false; // first neutralino jets
      isjet2[jetidx] = false; // second neutralino jets
      v.SetPtEtaPhiM( jet_pt, jet_eta, jet_phi, 0. ); //set the axis
      
    if ( jet_pt < PtMin ) continue;
    if ( abs(jet_eta) > EtaMax ) continue;
      
      // look if prompt muon inside
      float deltaR1 = 1000., deltaR2 = 1000.;
      if ( imu1 >= 0 ) deltaR1 = Deltar( jet_eta, jet_phi, tree_muon_eta[imu1], tree_muon_phi[imu1] );
      if ( imu2 >= 0 ) deltaR2 = Deltar( jet_eta, jet_phi, tree_muon_eta[imu2], tree_muon_phi[imu2] );
      if ( deltaR1 < 0.4 || deltaR2 < 0.4 )
      {
        if ( deltaR1 < 0.4 )
        { //if muon is inside, we remove the muons infomation from the jet
          v1.SetPtEtaPhiM( tree_muon_pt[imu1],
          		  tree_muon_eta[imu1],
          		  tree_muon_phi[imu1],
          		  0 );
          v -= v1; //v TLorentzFactor being just above, defined by jet data
        }
        if ( deltaR2 < 0.4 )
        {
          v2.SetPtEtaPhiM( tree_muon_pt[imu2],
          		  tree_muon_eta[imu2],
          		  tree_muon_phi[imu2],
          		  0 );
          v -= v2;
        }
        jet_pt  = v.Pt(); //Update jet data by removing the muons information (muons that could be in the jet)
        jet_eta = v.Eta(); //+ we do not want muons data to build the two axis since they come from the PV
        jet_phi = v.Phi();
      }
      
      njet++;
      isjet[jetidx] = true;
      vjet[jetidx] = v; // Only jet data (with  possible muons being removed)
      if ( njet1 == 0 && jet_pt > PtMin && abs(jet_eta) < EtaMax )
      {
        njet1 = 1;
        isjet1[jetidx] = true;
        vaxis1 = v;
      }
      jetidx++;
    } // End Loop on jets

    /////////////////////////////////////////////////////////
    //-------------------------------------------------------
    // Event Axes
    //-------------------------------------------------------
    /////////////////////////////////////////////////////////

    float dR, dR1 = 10., dR2 = 10.;
    float dRcut_hemis  = 1.5; // subjective choice
    float dRcut_tracks = 10.; // no cut is better (could bias low track pT and high LLP ct) 
     
    for (int i=0; i<jetidx; i++) // Loop on jet
    {
    if ( !isjet[i] ) continue;
      // float jet_pt  = vjet[i].Pt();
      float jet_eta = vjet[i].Eta();
      float jet_phi = vjet[i].Phi();
      if ( njet1 > 0 ) dR1 = Deltar( jet_eta, jet_phi, vaxis1.Eta(), vaxis1.Phi() );
      if ( njet2 > 0 ) dR2 = Deltar( jet_eta, jet_phi, vaxis2.Eta(), vaxis2.Phi() );
      // axis 1
      if ( njet1 > 0 && !isjet2[i]  && dR1 < dRcut_hemis) {
        njet1++;
        vaxis1 += vjet[i];
        isjet1[i] = true;
      }
      // axis 2
      if ( njet2 == 0 && !isjet1[i] ) {
        njet2 = 1;
        vaxis2 = vjet[i];
        isjet2[i] = true;
      }
      else if ( njet2 > 0 && !isjet1[i] && !isjet2[i] && dR2 < dRcut_hemis ) {//
        njet2++;
        vaxis2 += vjet[i];
        isjet2[i] = true;
      }
    }       // end Loop on jet
    
//$$
//     // force the axes to the true LLP
//     vaxis1 = vneu[0];
//     vaxis2 = vneu[1];
//$$

    ///////////////////////////////
    // compare with neutralino axis
    ///////////////////////////////
    
    int iLLPrec1 = 1, iLLPrec2 = 2;
    float axis1_eta = vaxis1.Eta();
    float axis1_phi = vaxis1.Phi();
    if ( neu[0] >= 0 ) dR1 = Deltar( axis1_eta, axis1_phi, Gen_neu1_eta, Gen_neu1_phi ); //dR between reco axis of jets and gen neutralino
    if ( neu[1] >= 0 ) dR2 = Deltar( axis1_eta, axis1_phi, Gen_neu2_eta, Gen_neu2_phi );
    dR = dR1;
    if ( dR2 < dR1 )
    { // make sure that the reco axis defined matches well with the axis of the gen neutralino, if not it is swapped
      iLLPrec1 = 2;
      iLLPrec2 = 1;
      dR = dR2;
    }
    float axis1_dR = dR;
    float axis2_eta = vaxis2.Eta();
    float axis2_phi = vaxis2.Phi();
    if ( njet2 == 0 )
    {  // compute an axis 2 even without jet, by taking the opposite in phi to axis 1
      axis2_eta = axis1_eta;
      axis2_phi = axis1_phi - 3.14159;
      if ( axis1_phi < 0 ) axis2_phi = axis1_phi + 3.14159;
    }
    if ( iLLPrec2 == 1 ) dR = Deltar( axis2_eta, axis2_phi, Gen_neu1_eta, Gen_neu1_phi );
    else                 dR = Deltar( axis2_eta, axis2_phi, Gen_neu2_eta, Gen_neu2_phi );
    float axis2_dR = dR;

    float dR_axis12 = Deltar(axis1_eta,axis1_phi,axis2_eta,axis2_phi);

        // cout << " njet1 " << njet1 << " and njet2" << njet2 << endl;
        // cout << " axis1_eta " << axis1_eta << " and axis2_eta" << axis2_eta << endl;
        // cout << " axis1_phi " << axis1_phi << " and axis2_phi" << axis2_phi << endl;
        // cout << " axis1_dR " << axis1_dR << " and axis2_dR" << axis2_dR << endl;
        // cout << " dR_axis12 " << dR_axis12 << endl;


    ///////////////////////////////////////////////////////
    //-----------------------------------------------------
    // selection of displaced tracks
    //-----------------------------------------------------
    ///////////////////////////////////////////////////////

    int jet; /*!*/
    float ntrk20, ntrk30, ntrk40; /*!*/
    float firsthit_X, firsthit_Y, firsthit_Z, phi;
    double bdtval = -100.;

    int nTrks_axis1 = 0;
    int nTrks_axis1_sig=0, nTrks_axis1_bad=0;
    int nTrks_axis2 = 0;
    int nTrks_axis2_sig=0, nTrks_axis2_bad=0;
    int nTrks_axis1_sig_mva=0, nTrks_axis1_bad_mva=0;
    int nTrks_axis2_sig_mva=0, nTrks_axis2_bad_mva=0;
    
    LLP1_nTrks = 0;
    LLP2_nTrks = 0;

//$$
//     double bdtcut = -0.0401; // for TMVAbgctau50withnhits.xml BDToldrecoavecalgo
//     double bdtcut = -0.0815; // for TMVAClassification_BDTG50sansalgo.weights.xml BDToldreco
//     double bdtcut =  0.0327; // for TMVAClassification_BDTG50cm_NewSignal.weights.xml BDTreco
//     double bdtcut = -0.1456; // for TMVAClassification_BDTG50cm_HighPurity.weights.xml BDTrecohp
//     double bdtcut = -0.1083; // for TMVAClassification_BDTG_FromBC.weights.xml from BDTminipf
//     double bdtcut = -0.0067; // for TMVAClassification_BDTG50cm_sansntrk10_avecHP.weights.xml BDTrecohpsansntrk10
//     double bdtcut = -10.; // no BDT cut
//$$
    double bdtcut = 0.5; 
    double bdtcut_step2 = -0.1456; // for TMVAClassification_BDTG50cm_HighPurity.weights.xml BDTrecohp
//$$

    //---------------------------//
    // if (tree_Filter){

    for (int counter_track = 0; counter_track < tree_nTracks; counter_track++) 
    {
      firsthit_X = tree_track_firstHit_x[counter_track];
      firsthit_Y = tree_track_firstHit_y[counter_track];
      firsthit_Z = tree_track_firstHit_z[counter_track];
      pt   = tree_track_pt[counter_track];
      eta  = tree_track_eta[counter_track];
      phi	 = tree_track_phi[counter_track];
      NChi	   = tree_track_NChi2[counter_track];
      nhits	   = tree_track_nHit[counter_track];
      // algo	 = tree_track_algo[counter_track];
      drSig	 = tree_track_drSig[counter_track];
      ntrk10 = 0;
      float ntrk10_lost = 0, ntrk20_lost = 0, ntrk30_lost = 0, ntrk40_lost = 0;
      isinjet = 0.;
      ntrk20 = 0;
      ntrk30 = 0;
      ntrk40 = 0;
      bdtval = -10.;
      dR = -1.;
      int tracks_axis = 0; // flag to check which axis is the closest from the track

      jet = tree_track_iJet[counter_track];
      if ( jet >= 0 ) isinjet = 1.; /*!*/
      int isFromLLP = tree_track_sim_LLP[counter_track];

      // check the dR between the tracks and the second axis (without any selection on the tracks)
      float dR1  = Deltar( eta, phi, axis1_eta, axis1_phi ); // axis1_phi and axis1_eta for the first axis
      float dR2  = Deltar( eta, phi, axis2_eta, axis2_phi );
      tracks_axis = 1;
      dR = dR1;
      if ( dR2 < dR1 ) { // a restriction could be added on the value of dR to assign the value Tracks_axis  (avoid some background???)
        tracks_axis = 2;
        dR = dR2;
      }

      //Computation of the distances needed for the BDT
      for (int counter_othertrack = 0; counter_othertrack < tree_nTracks; counter_othertrack++) 
      {
      if ( counter_othertrack == counter_track ) continue;
        float x2 = tree_track_firstHit_x[counter_othertrack];
        float y2 = tree_track_firstHit_y[counter_othertrack];
        float z2 = tree_track_firstHit_z[counter_othertrack];
        float dist = TMath::Sqrt( (firsthit_X-x2)*(firsthit_X-x2) + (firsthit_Y-y2)*(firsthit_Y-y2) + (firsthit_Z-z2)*(firsthit_Z-z2) ); // pour chaque reconstruite, on regarde les autres tracks
        if ( dist < 10. ) ntrk10++; 
        if ( dist < 20. ) ntrk20++;
        if ( dist < 30. ) ntrk30++;
        if ( dist < 40. ) ntrk40++;
	if ( tree_track_lost[counter_track] && tree_track_lost[counter_othertrack] ) {
          if ( dist < 10. ) ntrk10_lost++; 
          if ( dist < 20. ) ntrk20_lost++;
          if ( dist < 30. ) ntrk30_lost++;
          if ( dist < 40. ) ntrk40_lost++;
	}
      }  // end Loop on other Tracks
      if ( tree_track_lost[counter_track] ) {
        ntrk10 = ntrk10_lost;
        ntrk20 = ntrk20_lost;
        ntrk30 = ntrk30_lost;
        ntrk40 = ntrk40_lost;
      }

      if ( dR < dRcut_tracks ) 
      {
        if ( isFromLLP == 1 ) LLP1_nTrks++;
        if ( isFromLLP == 2 ) LLP2_nTrks++;

        bdtval = reader->EvaluateMVA( "BDTG" ); //default value = -10 (no -10 observed and -999 comes from EvaluateMVA)

        if ( tracks_axis == 1 ) {
          nTrks_axis1++;
          if ( isFromLLP == iLLPrec1 ) nTrks_axis1_sig++;
          else if ( isFromLLP >= 1 )   nTrks_axis1_bad++;
          if ( bdtval > bdtcut ) {
            if ( isFromLLP == iLLPrec1 ) nTrks_axis1_sig_mva++;
            else if ( isFromLLP >= 1 )   nTrks_axis1_bad_mva++;
          }
        }
      
        if ( tracks_axis == 2 ) {
          nTrks_axis2++;
          if ( isFromLLP == iLLPrec2 ) nTrks_axis2_sig++;
          else if ( isFromLLP >= 1 )   nTrks_axis2_bad++;
          if ( bdtval > bdtcut ) {
            if ( isFromLLP == iLLPrec2 ) nTrks_axis2_sig_mva++;
            else if ( isFromLLP >= 1 )   nTrks_axis2_bad_mva++;
          }
        }
      }
     
      tree_track_ntrk10.push_back(ntrk10);
      tree_track_ntrk20.push_back(ntrk20);
      tree_track_ntrk30.push_back(ntrk30);
      tree_track_ntrk40.push_back(ntrk40);
      tree_track_MVAval.push_back(bdtval);
      tree_track_Hemi.push_back(tracks_axis);
      tree_track_Hemi_dR.push_back(dR);
      if      ( tracks_axis == 1 ) tree_track_Hemi_LLP.push_back(iLLPrec1);
      else if ( tracks_axis == 2 ) tree_track_Hemi_LLP.push_back(iLLPrec2);
      else		           tree_track_Hemi_LLP.push_back(0);
      
    } //End loop on all the tracks
    

    /////////////////////////////////////////
    // Sort tracks by decreasing BDT value //
    /////////////////////////////////////////

    for (int counter_track = 0; counter_track < tree_nTracks; counter_track++)
      MVAval[counter_track] = tree_track_MVAval[counter_track];
    if ( tree_nTracks < 1000 ) {
      TMath::Sort(tree_nTracks, MVAval, index);
    }
    else cout << " !!!!!! ERROR tree_nTracks " << tree_nTracks << " ERROR !!!!!! " << endl;


    ///////////////////////////////
    // Fill the transient tracks //
    ///////////////////////////////

    vector<reco::TransientTrack> displacedTracks_llp1_mva, displacedTracks_llp2_mva; // Control Tracks
    vector<reco::TransientTrack> displacedTracks_Hemi1_mva, displacedTracks_Hemi2_mva; // Tracks selected wrt the hemisphere
    vector<reco::TransientTrack> displacedTracks_step2_Hemi1, displacedTracks_step2_Hemi2;

    for (int k = 0; k < tree_nTracks; k++)
    {
      int counter_track = index[k];
      unsigned int ipc = tree_track_ipc[counter_track];
      pat::PackedCandidateRef pcref = MINIgeneralTracks[ipc];
      const reco::Track *trackPcPtr = pcref->bestTrack();
      reco::Track tk;

      //PseudoDefTrack forces the covariance matrix to be positive definite. Can be negative in MINIAod with bestTrack/pseudoTrack methods
      //-----------------------------------------------------
      const reco::Track& tk_temp = *trackPcPtr;//const ..;reco::Track
      reco::TrackBase::CovarianceMatrix m_ = tk_temp.covariance(); //math::Error<5>::type
      double det=0;
      //Covariance Matrix Correction//
      bool notPosDef = (!(m_).Sub<AlgebraicSymMatrix22>(0, 0).Det(det) || det < 0) ||
                   (!(m_).Sub<AlgebraicSymMatrix33>(0, 0).Det(det) || det < 0) ||
                   (!(m_).Sub<AlgebraicSymMatrix44>(0, 0).Det(det) || det < 0) || (!(m_).Det(det) || det < 0);
      if ( notPosDef && NewCovMat ) {//flag to allow for covariance matrix corrections
        reco::TrackBase::CovarianceMatrix m(m_);
        // if not positive-definite, alter values to allow for pos-def
        TMatrixDSym eigenCov(5);
        for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
            if (std::isnan((m)(i, j)) || std::isinf((m)(i, j)))
              eigenCov(i, j) = 1e-6;
            else
              eigenCov(i, j) = (m)(i, j);
          }
        }
        TVectorD eigenValues(5);
        eigenCov.EigenVectors(eigenValues);
        double minEigenValue = eigenValues.Min();
        double delta = 1e-6;
        if (minEigenValue < 0) {
          for (int i = 0; i < 5; i++)
            m(i, i) += delta - minEigenValue;
        }
        // make a track object with pos def covariance matrix  
        tk = reco::Track(tk_temp.normalizedChi2() * tk_temp.ndof(),tk_temp.ndof(),tk_temp.TrackBase::referencePoint(),tk_temp.momentum(),tk_temp.charge(),m,reco::TrackBase::undefAlgorithm,reco::TrackBase::loose);
        // for(int i = 0 ; i<4 ;i++ )
        //   {
        //     for(int j = 0 ; j<4 ;j++ )
        //	 {
        //	     std::cout<<" Apres Covcor m["<<i<<"]["<<j<<"]="<<m[i][j]<<std::endl;
        //	 }
        //   }
      //-------------------end covariance matrix correction-----------//
      }
      else tk = *trackPcPtr;

      int isFromLLP   = tree_track_sim_LLP[counter_track];
      int tracks_axis = tree_track_Hemi[counter_track];
      double bdtval   = tree_track_MVAval[counter_track];
      if ( bdtval > bdtcut ) {
        if ( isFromLLP == 1 )
          displacedTracks_llp1_mva.push_back(theTransientTrackBuilder->build(trackPcPtr));
        if ( isFromLLP == 2 )
          displacedTracks_llp2_mva.push_back(theTransientTrackBuilder->build(trackPcPtr));
        if ( tracks_axis == 1 )
          displacedTracks_Hemi1_mva.push_back(theTransientTrackBuilder->build(&tk));
        if ( tracks_axis == 2 )
          displacedTracks_Hemi2_mva.push_back(theTransientTrackBuilder->build(&tk));
      }
      if ( bdtval > bdtcut_step2 ) {
        if ( tracks_axis == 1 )
          displacedTracks_step2_Hemi1.push_back(theTransientTrackBuilder->build(&tk));
        if ( tracks_axis == 2 )
          displacedTracks_step2_Hemi2.push_back(theTransientTrackBuilder->build(&tk));
      }
    }

  
    
        // cout << " displaced tracks LLP1 " << LLP1_nTrks << " and with mva" << displacedTracks_llp1_mva.size() << endl;
        // cout << " displaced tracks LLP2 " << LLP2_nTrks << " and with mva" << displacedTracks_llp2_mva.size() << endl;
        // cout << " displaced tracks Hemi1 " << nTrks_axis1 << " and with mva" << displacedTracks_Hemi1_mva.size() << endl;
        // cout << " displaced tracks Hemi2 " << nTrks_axis2 << " and with mva" << displacedTracks_Hemi2_mva.size() << endl;

    ///////////////////////////////////////////////////////
    //-----------------------------------------------------
    // Vertex fitting 
    //-----------------------------------------------------
    ///////////////////////////////////////////////////////

    int   Vtx_ntk = 0, Vtx_step = 0;
    float Vtx_x = 0., Vtx_y = 0., Vtx_z= 0., Vtx_chi = -10.;
    float recX, recY, recZ, dSV, recD;

//$$
    // parameters for the Adaptive Vertex Fitter (AVF)
    double maxshift        = 0.0001;
    unsigned int maxstep   = 30;
    double maxlpshift      = 0.1;
    double weightThreshold = 0.001;
    double sigmacut        = 3.;
    double Tini            = 256.;
    double ratio           = 0.25;
//$$

//------------------------------- FIRST LLP WITH MVA ----------------------------------//

    static AdaptiveVertexFitter 
    theFitter_vertex_llp1_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_vertex_llp1_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
    
    Vtx_ntk = displacedTracks_llp1_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -10.;
    
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_llp1_mva = theFitter_vertex_llp1_mva.vertex(displacedTracks_llp1_mva); // fitted vertex
      // std::cout<< "displacedVertex_llp1_mva is built" << std::endl;
      if ( displacedVertex_llp1_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      {
        Vtx_ntk = displacedTracks_llp1_mva.size();
        Vtx_x = displacedVertex_llp1_mva.position().x();
        Vtx_y = displacedVertex_llp1_mva.position().y();
        Vtx_z = displacedVertex_llp1_mva.position().z();
        Vtx_chi = displacedVertex_llp1_mva.normalisedChiSquared();
      }
    }

    tree_LLP.push_back(1);
    tree_LLP_pt.push_back(   LLP1_pt);
    tree_LLP_eta.push_back(  LLP1_eta);
    tree_LLP_phi.push_back(  LLP1_phi);
    tree_LLP_x.push_back(    LLP1_x);
    tree_LLP_y.push_back(    LLP1_y);
    tree_LLP_z.push_back(    LLP1_z);
    tree_LLP_dist.push_back( LLP1_dist);
    tree_LLP_nTrks.push_back(LLP1_nTrks);
    tree_LLP_Vtx_nTrks.push_back(Vtx_ntk);
    tree_LLP_Vtx_NChi2.push_back(Vtx_chi);
    tree_LLP_Vtx_dx.push_back(Vtx_x - LLP1_x);
    tree_LLP_Vtx_dy.push_back(Vtx_y - LLP1_y);
    tree_LLP_Vtx_dz.push_back(Vtx_z - LLP1_z);

    dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
    recX = Vtx_x - tree_PV_x[0];
    recY = Vtx_y - tree_PV_y[0];
    recZ = Vtx_z - tree_PV_z[0];
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_LLP_Vtx_dist.push_back( recD );
    tree_LLP_Vtx_dd.push_back( TMath::Sqrt(dSV)/LLP1_dist );
    
//&&&&&
// //     bool dump = false;
// //     if ( Vtx_chi < 0 && LLP1_nTrks > 1 ) dump = true;
//     cout << endl;
//     cout << endl;
//     cout << " &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& " << endl;
//     cout << " run event " << runNumber << " " << eventNumber << endl;
//     cout << " LLP " << 1 << " pt eta phi " << LLP1_pt << " " << LLP1_eta << " " << LLP1_phi << " x y z " << LLP1_x << " " << LLP1_y << " " << LLP1_z << " nTrks " << LLP1_nTrks << endl;
//     cout << "	  Vtx Chi2 " << Vtx_chi << " dx dy dz " << Vtx_x - LLP1_x << " " << Vtx_y - LLP1_y  << " " << Vtx_z - LLP1_z << " nTrks " << Vtx_ntk << endl;
//&&&&&


    //-------------------------- SECOND LLP WITH MVA -------------------------------------//

    static AdaptiveVertexFitter 
    theFitter_vertex_llp2_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_vertex_llp2_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
    
    Vtx_ntk = displacedTracks_llp2_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -10.;
    
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_llp2_mva = theFitter_vertex_llp2_mva.vertex(displacedTracks_llp2_mva); // fitted vertex
      
      if ( displacedVertex_llp2_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      {
        Vtx_x = displacedVertex_llp2_mva.position().x();
        Vtx_y = displacedVertex_llp2_mva.position().y();
        Vtx_z = displacedVertex_llp2_mva.position().z();
        Vtx_chi = displacedVertex_llp2_mva.normalisedChiSquared();
      }
    }

    tree_LLP.push_back(2);
    tree_LLP_pt.push_back(   LLP2_pt);
    tree_LLP_eta.push_back(  LLP2_eta);
    tree_LLP_phi.push_back(  LLP2_phi);
    tree_LLP_x.push_back(    LLP2_x);
    tree_LLP_y.push_back(    LLP2_y);
    tree_LLP_z.push_back(    LLP2_z);
    tree_LLP_dist.push_back( LLP2_dist);
    tree_LLP_nTrks.push_back(LLP2_nTrks);
    tree_LLP_Vtx_nTrks.push_back(Vtx_ntk);
    tree_LLP_Vtx_NChi2.push_back(Vtx_chi);
    tree_LLP_Vtx_dx.push_back(Vtx_x - LLP2_x);
    tree_LLP_Vtx_dy.push_back(Vtx_y - LLP2_y);
    tree_LLP_Vtx_dz.push_back(Vtx_z - LLP2_z);

    dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
    recX = Vtx_x - tree_PV_x[0];
    recY = Vtx_y - tree_PV_y[0];
    recZ = Vtx_z - tree_PV_z[0];
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_LLP_Vtx_dist.push_back( recD );
    tree_LLP_Vtx_dd.push_back( TMath::Sqrt(dSV)/LLP2_dist );
    
//&&&&&
// //     if ( Vtx_chi < 0 && LLP2_nTrks > 1 ) dump = true;
//     cout << endl;
//     cout << " LLP " << 2 << " pt eta phi " << LLP2_pt << " " << LLP2_eta << " " << LLP2_phi << " x y z " << LLP2_x << " " << LLP2_y << " " << LLP2_z << " nTrks " << LLP2_nTrks << endl;
//     cout << "	  Vtx Chi2 " << Vtx_chi << " dx dy dz " << Vtx_x - LLP2_x << " " << Vtx_y - LLP2_y  << " " << Vtx_z - LLP2_z << " nTrks " << Vtx_ntk << endl;
//     cout << " &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& " << endl;
// //     cout << endl;
//&&&&&

     
    //--------------------------- FIRST HEMISPHERE WITH MVA -------------------------------------//

    static AdaptiveVertexFitter 
    theFitter_Vertex_Hemi1_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_Vertex_Hemi1_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
    
    Vtx_ntk = displacedTracks_Hemi1_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -10.;
    Vtx_step = 0;

    TransientVertex displacedVertex_Hemi1_mva;
	
    if ( Vtx_ntk > 1 )
    {
      displacedVertex_Hemi1_mva = theFitter_Vertex_Hemi1_mva.vertex(displacedTracks_Hemi1_mva); // fitted vertex
      if ( displacedVertex_Hemi1_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      { 
        Vtx_x = displacedVertex_Hemi1_mva.position().x();
        Vtx_y = displacedVertex_Hemi1_mva.position().y();
        Vtx_z = displacedVertex_Hemi1_mva.position().z();
        Vtx_chi = displacedVertex_Hemi1_mva.normalisedChiSquared();
	Vtx_step = 1;
//         for (int p = 0; p < Vtx_ntk; p++)
//         {
//           tree_Hemi_Vtx_trackWeight.push_back(displacedVertex_Hemi1_mva.trackWeight(displacedTracks_Hemi1_mva[p]));
//           // if (showlog) std::cout<<"1st vtx_chi / weight / chi2 / ndof / NCHi2: "<<Vtx_chi<<" / " <<displacedVertex_Hemi1_mva.trackWeight(displacedTracks_Hemi1_mva[p])<<" / "<<displacedTracks_Hemi1_mva[p].ndof()<<" / "<<displacedTracks_Hemi1_mva[p].normalizedChi2()<<std::endl;  	  
//         }
      }
    }

    // step 2
    TransientVertex displacedVertex_step2_Hemi1;
    static AdaptiveVertexFitter theFitter_Vertex_step2_Hemi1(
    	       GeometricAnnealing ( sigmacut, Tini, ratio ), 
    	       DefaultLinearizationPointFinder(),
    	       KalmanVertexUpdator<5>(), 
    	       KalmanVertexTrackCompatibilityEstimator<5>(), 
    	       KalmanVertexSmoother() );
    theFitter_Vertex_step2_Hemi1.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
    bool badVtx = false;
    if ( Vtx_chi < 0. || Vtx_chi > 10. ) badVtx = true;
    if ( badVtx && displacedTracks_step2_Hemi1.size() > 1 ) {
      Vtx_ntk = displacedTracks_step2_Hemi1.size();
      Vtx_chi = -10.;
      displacedVertex_step2_Hemi1 = theFitter_Vertex_step2_Hemi1.vertex(displacedTracks_step2_Hemi1);
      if ( displacedVertex_step2_Hemi1.isValid() )
      { 
        Vtx_x	= displacedVertex_step2_Hemi1.position().x();
        Vtx_y	= displacedVertex_step2_Hemi1.position().y();
        Vtx_z	= displacedVertex_step2_Hemi1.position().z();
        Vtx_chi = displacedVertex_step2_Hemi1.normalisedChiSquared();
        Vtx_step = 2;
      }
    }

    //----------------------------------------IAVF-----------------------------------------------//
    //                           Iterative Adaptive Vertex Fitter                                //
    //Input : Collections od displaced Tracks ordered by decreasing values of BDT => to have the //
    //        the best tracks at the beginning of the collection                                 //
    //                                                                                           //
    //Process : Within the collection of tracks, look for the first good seed made of 2 tracks,  //
    //          (good seed : vertex valid with a chi2 within a certain range (LowerLimit and     //
    //          UpperLimit)). From this seed, tracks are added one by one, and a vertex is built //
    //          at each step. The vertex is either valid or not. If valid in a certain range of  //
    //          Chi2, we keep going until there are no more tracks.                              //
    //                                                                                           //
    //Degrees of freedom: -LowerLimit (Chi2 restriction)                                         //
    //                    -UpperLimit (Chi2 restriction)                                         //
    //                    -MeanBDtValue (Restriction on the BDT value of the tracks that formed  //
    //                     the vertex)                                                           //
    //                                                                                           //
    //Why? : Address the drop in efficiency due to MiniAOD datatier, turns out it improves the   //
    //       basic AVF implementation (even in RECO/AOD)                                         //
    //                                                                                           //
    //PS: Maximum efficiency is reached for MiniAOD when using the covariance matrix correction  //
    //-------------------------------------------------------------------------------------------//

    // step 3
    int ntracks    = -2;
    float tempchi2 = -10.;
    float tempx = -100.;
    float tempy = -100.;
    float tempz = -100.;
    badVtx = false;
    if ( Vtx_chi < 0. || Vtx_chi > 10. ) badVtx = true;
    if ( badVtx && IterAVF && Vtx_ntk > 1 ) {
      bool success = false;
      std::vector<TransientTrack> vTT;
      tempchi2 = -10.;
      for ( int p = 1; p < Vtx_ntk; p++ )
      {
        for  ( int k = 0; k < p; k++ ) // take pairs of tracks of highest BDT value
        {
          vTT.push_back(displacedTracks_step2_Hemi1[p]);
          vTT.push_back(displacedTracks_step2_Hemi1[k]);
          ntracks = 2;
          TransientVertex TV = theFitter_Vertex_step2_Hemi1.vertex(vTT); // We take the first "good-looking" seed to start
          if ( TV.isValid() && TV.normalisedChiSquared()>0 && TV.normalisedChiSquared()<10 ) {
	    tempchi2 = TV.normalisedChiSquared();
	    tempx = TV.position().x();
	    tempy = TV.position().y();
	    tempz = TV.position().z();
	  }
          if ( TV.isValid() && TV.normalisedChiSquared()>0 && TV.normalisedChiSquared()<10 )
          {
	    success = true; 
            if (showlog) std::cout<<"1st LLP seed is created for k and p : "<<k<<" / "<<p<<" with chi2: "<<TV.normalisedChiSquared()<<std::endl;
            for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
            {
            if (m == k || m == p) continue;
              ntracks++;
              vTT.push_back(displacedTracks_step2_Hemi1[m]);
              TransientVertex updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
              if (showlog) std::cout<<"m = "<<m<<"/ chi2 of the vetex constructed with this seed : "<<updatedTV.normalisedChiSquared()<<" originally : "<<Vtx_chi<<std::endl;
              if ( updatedTV.isValid() ) tempchi2 = updatedTV.normalisedChiSquared();
//               DOF = updatedTV.degreesOfFreedom();
//               temptChi2 = updatedTV.totalChiSquared();
//               if ( !updatedTV.isValid() || tempchi2 < 0. || tempchi2 > 10. ) 
              if ( !updatedTV.isValid() ) 
	      {  
	    	vTT.pop_back();
	    	ntracks--;
	    	updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
	    	tempchi2 = updatedTV.normalisedChiSquared();
	    	tempx=updatedTV.position().x();
	    	tempy=updatedTV.position().y();
        	tempz=updatedTV.position().z(); 
// 	    	DOF = updatedTV.degreesOfFreedom();
// 	    	temptChi2 = updatedTV.totalChiSquared();
	    	continue;
	      } 
              else
	      {
	    	tempx=updatedTV.position().x();
	    	tempy=updatedTV.position().y();
	    	tempz=updatedTV.position().z();
	      }
            } // end loop on the other tracks
            Vtx_ntk = ntracks;
            Vtx_chi = tempchi2;
            Vtx_x = tempx;
            Vtx_y = tempy;
            Vtx_z = tempz;
            Vtx_step = 3;
            // if (showlog) std::cout << "OriginalTracks : " << Vtx_ntk << " vs Actual tracks : "<<ntracks<<"and Final chi2 of : "<<Vtx_chi<<std::endl; 
          }
          else
          { // If not valid : we build a new seed
            ntracks = 0;
            vTT.clear();
          }
          if ( success ) break;
        }
        if ( success ) break;
      }
      //--------------------ENDOF IAVF--------------------------//
    }

//$$$$   float Vtx_r = sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y);
//$$$$   PropaHitPattern* NI = new PropaHitPattern();
   int VtxLayerNI = -1;
//$$$$  VtxLayerNI = NI->VertexBelongsToBarrelLayer(Vtx_r); // Only for Barrel atm
//$$$$   if ( VtxLayerNI == 0 ) VtxLayerNI=NI->VertexBelongsToDiskLayer(Vtx_z);
   tree_Hemi_Vtx_Layer.push_back(VtxLayerNI);
   
     bool L0_Vtx = false;
  if ( !KshortVertex->empty() ) {
  for (int j =0 ; j<tree_nLambda ; j++)
    {
      float L0x = (*LambdaVertex)[j].position().x();
      float L0y = (*LambdaVertex)[j].position().y();
      float L0z = (*LambdaVertex)[j].position().z();     
      float dd_L0 = sqrt((Vtx_x-L0x)*(Vtx_x-L0x)+(Vtx_y-L0y)*(Vtx_y-L0y)+(Vtx_z-L0z)*(Vtx_z-L0z));
      if (dd_L0<0.5 && Vtx_ntk <= 3) // cm
        {
          L0_Vtx = true;
          break;
        }
    }
  }
  tree_Hemi_Vtx_L0.push_back(L0_Vtx);
  bool K0_Vtx = false;
  if ( !KshortVertex->empty() && L0_Vtx==false ) {
  for (int j =0 ; j<tree_nK0 ; j++)
    {
      float K0x = (*KshortVertex)[j].position().x();
      float K0y = (*KshortVertex)[j].position().y();
      float K0z = (*KshortVertex)[j].position().z();     
      float dd_K0 = sqrt((Vtx_x-K0x)*(Vtx_x-K0x)+(Vtx_y-K0y)*(Vtx_y-K0y)+(Vtx_z-K0z)*(Vtx_z-K0z));
      if (dd_K0<0.5 && Vtx_ntk <= 3) // cm
        {
          K0_Vtx = true;
          break;
        }
    }
  }
  tree_Hemi_Vtx_K0.push_back(K0_Vtx);
  if (K0_Vtx || L0_Vtx ) tree_Hemi_Vtx_V0.push_back(true);
  else  tree_Hemi_Vtx_V0.push_back(false);

  bool Yc_Vtx = false;
  if ( !PhotonConversion->empty() && !K0_Vtx && !L0_Vtx ) {
    for (int j = 0; j<tree_nYConv ; j++ )
    {
      if((*PhotonConversion)[j].isConverted())//tracksize>0       /// Number of tracks= 0,1,2      unsigned int nTracks() const {return  tracks().size(); }
      {
        float Yx = (*PhotonConversion)[j].conversionVertex().position().x();
        float Yy = (*PhotonConversion)[j].conversionVertex().position().y();
        float Yz = (*PhotonConversion)[j].conversionVertex().position().z();
        float dd_Yc = sqrt((Vtx_x-Yx)*(Vtx_x-Yx)+(Vtx_y-Yy)*(Vtx_y-Yy)+(Vtx_z-Yz)*(Vtx_z-Yz));
        if (dd_Yc<0.5 && Vtx_ntk <= 3) // cm
        {
          Yc_Vtx = true;
          break;
        }
      }
    }
  }
  tree_Hemi_Vtx_Yc.push_back(Yc_Vtx);

    float Vtx_chi1 = Vtx_chi;
    int nVtx = 0;
    if (abs(Vtx_chi)<10 && Vtx_ntk>1){nVtx++;}
    tree_Hemi.push_back(1);
    tree_Hemi_njet.push_back(njet1);
    tree_Hemi_eta.push_back(axis1_eta);
    tree_Hemi_phi.push_back(axis1_phi);
    tree_Hemi_dR.push_back(axis1_dR);
    tree_Hemi_nTrks.push_back(nTrks_axis1);
    tree_Hemi_nTrks_sig.push_back(nTrks_axis1_sig);
    tree_Hemi_nTrks_bad.push_back(nTrks_axis1_bad);
    tree_Hemi_Vtx_step.push_back(Vtx_step);
    tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
    tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
    tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis1_sig_mva);
    tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis1_bad_mva);
    tree_Hemi_Vtx_x.push_back(Vtx_x);
    tree_Hemi_Vtx_y.push_back(Vtx_y);
    tree_Hemi_Vtx_r.push_back(sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y));
    tree_Hemi_Vtx_z.push_back(Vtx_z);
    recX = Vtx_x - tree_PV_x[0];
    recY = Vtx_y - tree_PV_y[0];
    recZ = Vtx_z - tree_PV_z[0];
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_Hemi_Vtx_dist.push_back( recD );

    float ddok, ddbad;
    float ping1 = 0;

    if ( iLLPrec1 == 1 ) {
      tree_Hemi_LLP_pt.push_back( LLP1_pt);
      tree_Hemi_LLP_eta.push_back(LLP1_eta);
      tree_Hemi_LLP_phi.push_back(LLP1_phi);
      tree_Hemi_LLP_x.push_back(LLP1_x);
      tree_Hemi_LLP_y.push_back(LLP1_y);
      tree_Hemi_LLP_z.push_back(LLP1_z);
      tree_Hemi_LLP_dist.push_back(LLP1_dist);
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP1_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP1_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP1_z);
      tree_Hemi_Vtx_dr.push_back(TMath::Sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y) - TMath::Sqrt(LLP1_x*LLP1_x+LLP1_y*LLP1_y));

      dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
      ddok = TMath::Sqrt(dSV)/LLP1_dist;
      tree_Hemi_Vtx_dd.push_back( ddok );
      dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
      ddbad = TMath::Sqrt(dSV)/LLP2_dist;
      tree_Hemi_Vtx_ddbad.push_back( ddbad );
      if      ( Vtx_chi >= 0. && Vtx_chi < 10. && ddok  < 0.1 ) ping1 = 1;
      else if ( Vtx_chi >= 0. && Vtx_chi < 10. && ddbad < 0.1 ) ping1 = 2;
    }
    else {
      tree_Hemi_LLP_pt.push_back( LLP2_pt);
      tree_Hemi_LLP_eta.push_back(LLP2_eta);
      tree_Hemi_LLP_phi.push_back(LLP2_phi);
      tree_Hemi_LLP_x.push_back(LLP2_x);
      tree_Hemi_LLP_y.push_back(LLP2_y);
      tree_Hemi_LLP_z.push_back(LLP2_z);
      tree_Hemi_LLP_dist.push_back(LLP2_dist);
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP2_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP2_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP2_z);
      tree_Hemi_Vtx_dr.push_back(TMath::Sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y) - TMath::Sqrt(LLP2_x*LLP2_x+LLP2_y*LLP2_y));

      dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
      ddok = TMath::Sqrt(dSV)/LLP2_dist;
      tree_Hemi_Vtx_dd.push_back( ddok );
      dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
      ddbad = TMath::Sqrt(dSV)/LLP1_dist;
      tree_Hemi_Vtx_ddbad.push_back( ddbad );
      if      ( Vtx_chi >= 0. && Vtx_chi < 10. && ddok  < 0.1 ) ping1 = 2;
      else if ( Vtx_chi >= 0. && Vtx_chi < 10. && ddbad < 0.1 ) ping1 = 1;
    }
    tree_Hemi_LLP.push_back(iLLPrec1);
         

    //--------------------------------------------------------------------------------------------//
    //--------------------------- SECOND HEMISPHERE WITH MVA -------------------------------------//
    //--------------------------------------------------------------------------------------------//

    static AdaptiveVertexFitter 
    theFitter_Vertex_Hemi2_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_Vertex_Hemi2_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );

    Vtx_ntk = displacedTracks_Hemi2_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -10.;
    Vtx_step = 0;
    TransientVertex displacedVertex_Hemi2_mva;

    if ( Vtx_ntk > 1 )
    {
      displacedVertex_Hemi2_mva = theFitter_Vertex_Hemi2_mva.vertex(displacedTracks_Hemi2_mva); // fitted vertex
      
      if ( displacedVertex_Hemi2_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      {
        Vtx_x = displacedVertex_Hemi2_mva.position().x();
        Vtx_y = displacedVertex_Hemi2_mva.position().y();
        Vtx_z = displacedVertex_Hemi2_mva.position().z();
        Vtx_chi = displacedVertex_Hemi2_mva.normalisedChiSquared();
	Vtx_step = 1;
//         for (int p =0; p<Vtx_ntk; p++)
//         {
//           tree_Hemi_Vtx_trackWeight.push_back(displacedVertex_Hemi2_mva.trackWeight(displacedTracks_Hemi2_mva[p]));
//           //  if(showlog) std::cout<<"2nd vtx_chi  / weight / chi2 / ndof / NChi2 : "<<Vtx_chi<<" / " <<displacedVertex_Hemi2_mva.trackWeight(displacedTracks_Hemi2_mva[p])<<" / "<<displacedTracks_Hemi2_mva[p].chi2()<<" / "<<displacedTracks_Hemi2_mva[p].ndof()<<" / "<<displacedTracks_Hemi2_mva[p].normalizedChi2()<<std::endl;
//         }
      }
    }    

    // step 2
    TransientVertex displacedVertex_step2_Hemi2;
    static AdaptiveVertexFitter theFitter_Vertex_step2_Hemi2(
    	       GeometricAnnealing ( sigmacut, Tini, ratio ), 
    	       DefaultLinearizationPointFinder(),
    	       KalmanVertexUpdator<5>(), 
    	       KalmanVertexTrackCompatibilityEstimator<5>(), 
    	       KalmanVertexSmoother() );
    theFitter_Vertex_step2_Hemi2.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
    badVtx = false;
    if ( Vtx_chi < 0. || Vtx_chi > 10. ) badVtx = true;
    if ( badVtx && displacedTracks_step2_Hemi2.size() > 1 ) {
      Vtx_ntk = displacedTracks_step2_Hemi2.size();
      Vtx_chi = -10.;
      displacedVertex_step2_Hemi2 = theFitter_Vertex_step2_Hemi2.vertex(displacedTracks_step2_Hemi2);
      if ( displacedVertex_step2_Hemi2.isValid() )
      { 
        Vtx_x	= displacedVertex_step2_Hemi2.position().x();
        Vtx_y	= displacedVertex_step2_Hemi2.position().y();
        Vtx_z	= displacedVertex_step2_Hemi2.position().z();
        Vtx_chi = displacedVertex_step2_Hemi2.normalisedChiSquared();
        Vtx_step = 2;
      }
    }

    // step 3
    ntracks   = -2;
    tempchi2  = -10.;
    tempx = -100.;
    tempy = -100.;
    tempz = -100.;
    badVtx = false;
    if ( Vtx_chi < 0. || Vtx_chi > 10. ) badVtx = true;
    if ( badVtx && IterAVF && Vtx_ntk > 1 ) {
      bool success = false;
      std::vector<TransientTrack> vTT;
      tempchi2 = -10.;
      for ( int p = 1; p < Vtx_ntk; p++ )
      {
        for  ( int k = 0; k < p; k++ ) // take pairs of tracks of highest BDT value
        {
          vTT.push_back(displacedTracks_step2_Hemi2[p]);
          vTT.push_back(displacedTracks_step2_Hemi2[k]);
          ntracks = 2;
          TransientVertex TV = theFitter_Vertex_step2_Hemi2.vertex(vTT); // We take the first "good-looking" seed to start
          if ( TV.isValid() && TV.normalisedChiSquared()>0 && TV.normalisedChiSquared()<10 ) {
	    tempchi2 = TV.normalisedChiSquared();
	    tempx = TV.position().x();
	    tempy = TV.position().y();
	    tempz = TV.position().z();
	  }
          if ( TV.isValid() && TV.normalisedChiSquared()>0 && TV.normalisedChiSquared()<10 )
          { 
            success = true;
            if (showlog) std::cout<<"1st LLP seed is created for k and p : "<<k<<" / "<<p<<" with chi2: "<<TV.normalisedChiSquared()<<std::endl;
            for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
            {
            if (m == k || m == p) continue;
              ntracks++;
              vTT.push_back(displacedTracks_step2_Hemi2[m]);
              TransientVertex updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
              if (showlog) std::cout<<"m = "<<m<<"/ chi2 of the vetex constructed with this seed : "<<updatedTV.normalisedChiSquared()<<" originally : "<<Vtx_chi<<std::endl;
              if ( updatedTV.isValid() ) tempchi2 = updatedTV.normalisedChiSquared();
//               DOF = updatedTV.degreesOfFreedom();
//               temptChi2 = updatedTV.totalChiSquared();
//               if ( !updatedTV.isValid() || tempchi2 < 0. || tempchi2 > 10. ) 
              if ( !updatedTV.isValid() ) 
	      {  
	    	vTT.pop_back();
	    	ntracks--;
	    	updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
	    	tempchi2 = updatedTV.normalisedChiSquared();
	    	tempx=updatedTV.position().x();
	    	tempy=updatedTV.position().y();
        	tempz=updatedTV.position().z(); 
// 	    	DOF = updatedTV.degreesOfFreedom();
// 	    	temptChi2 = updatedTV.totalChiSquared();
	    	continue;
	      } 
              else
	      {
	    	tempx=updatedTV.position().x();
	    	tempy=updatedTV.position().y();
	    	tempz=updatedTV.position().z();
	      }
            } // end loop on the other tracks
            Vtx_ntk = ntracks;
            Vtx_chi = tempchi2;
            Vtx_x = tempx;
            Vtx_y = tempy;
            Vtx_z = tempz;
            Vtx_step = 3;
            // if (showlog) std::cout << "OriginalTracks : " << Vtx_ntk << " vs Actual tracks : "<<ntracks<<"and Final chi2 of : "<<Vtx_chi<<std::endl; 
          }
          else
          { // If not valid : we build a new seed
            ntracks = 0;
            vTT.clear();
          }
          if ( success ) break;
        }
        if ( success ) break;
      }
      //--------------------ENDOF IAVF-------------------------//
    }
//$$$$   Vtx_r = sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y);
   VtxLayerNI = -1;
//$$$$   VtxLayerNI = NI->VertexBelongsToBarrelLayer(Vtx_r);//Only for Barrel atm
//$$$$   if (VtxLayerNI == 0) VtxLayerNI=NI->VertexBelongsToDiskLayer(Vtx_z);
   tree_Hemi_Vtx_Layer.push_back(VtxLayerNI);

   L0_Vtx = false;
  if ( !KshortVertex->empty() ) {
  for (int j =0 ; j<tree_nLambda ; j++)
    {
      float L0x = (*LambdaVertex)[j].position().x();
      float L0y = (*LambdaVertex)[j].position().y();
      float L0z = (*LambdaVertex)[j].position().z();     
      float dd_L0 = sqrt((Vtx_x-L0x)*(Vtx_x-L0x)+(Vtx_y-L0y)*(Vtx_y-L0y)+(Vtx_z-L0z)*(Vtx_z-L0z));
      if ( dd_L0<0.5 && Vtx_ntk<=3 ) // cm
        {
          L0_Vtx = true;
          break;
        }
    }
  }
  tree_Hemi_Vtx_L0.push_back(L0_Vtx);

  K0_Vtx = false;
  if ( !KshortVertex->empty() && L0_Vtx==false ) {
  for (int j =0 ; j<tree_nK0 ; j++)
    {
      float K0x = (*KshortVertex)[j].position().x();
      float K0y = (*KshortVertex)[j].position().y();
      float K0z = (*KshortVertex)[j].position().z();     
      float dd_K0 = sqrt((Vtx_x-K0x)*(Vtx_x-K0x)+(Vtx_y-K0y)*(Vtx_y-K0y)+(Vtx_z-K0z)*(Vtx_z-K0z));
      if (dd_K0<0.5 && Vtx_ntk<=3 ) //cm
        {
          K0_Vtx=true;
          break;
        }
    } 
  }
  tree_Hemi_Vtx_K0.push_back(K0_Vtx);
  if (K0_Vtx || L0_Vtx) tree_Hemi_Vtx_V0.push_back(true);
  else tree_Hemi_Vtx_V0.push_back(false);

  Yc_Vtx = false;
  if ( !PhotonConversion->empty() && !K0_Vtx && !L0_Vtx) {
    for (int j = 0; j<tree_nYConv ; j++ )
    {
      if((*PhotonConversion)[j].isConverted())//tracksize>0       /// Number of tracks= 0,1,2      unsigned int nTracks() const {return  tracks().size(); }
      {
        float Yx = (*PhotonConversion)[j].conversionVertex().position().x();
        float Yy = (*PhotonConversion)[j].conversionVertex().position().y();
        float Yz = (*PhotonConversion)[j].conversionVertex().position().z();
        float dd_Yc = sqrt((Vtx_x-Yx)*(Vtx_x-Yx)+(Vtx_y-Yy)*(Vtx_y-Yy)+(Vtx_z-Yz)*(Vtx_z-Yz));
        if (dd_Yc<0.5 && Vtx_ntk<=3 ) //cm 
        {
          Yc_Vtx = true;
          break;
        }
      }
    }
  }
  tree_Hemi_Vtx_Yc.push_back(Yc_Vtx);

    //---------------End of Secondary Interactions----------------// 

    float Vtx_chi2 = Vtx_chi;
    if (abs(Vtx_chi)<10 && Vtx_ntk>1){nVtx++;}
    if (nVtx==2){tree_Hemi_Vtx_evt2vtx.push_back(nVtx);}
    tree_Hemi.push_back(2);
    tree_Hemi_njet.push_back(njet2);
    tree_Hemi_eta.push_back(axis2_eta);
    tree_Hemi_phi.push_back(axis2_phi);
    tree_Hemi_dR.push_back(axis2_dR);
    tree_Hemi_nTrks.push_back(nTrks_axis2);
    tree_Hemi_nTrks_sig.push_back(nTrks_axis2_sig);
    tree_Hemi_nTrks_bad.push_back(nTrks_axis2_bad);
    tree_Hemi_Vtx_step.push_back(Vtx_step);
    tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
    tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
    tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis2_sig_mva);
    tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis2_bad_mva);
    tree_Hemi_Vtx_x.push_back(Vtx_x);
    tree_Hemi_Vtx_y.push_back(Vtx_y);
    tree_Hemi_Vtx_r.push_back(sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y));
    tree_Hemi_Vtx_z.push_back(Vtx_z);
    recX = Vtx_x - tree_PV_x[0];
    recY = Vtx_y - tree_PV_y[0];
    recZ = Vtx_z - tree_PV_z[0];
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_Hemi_Vtx_dist.push_back( recD );

    float ping2 = 0;
    if ( iLLPrec2 == 1 ) {
      tree_Hemi_LLP_pt.push_back( LLP1_pt);
      tree_Hemi_LLP_eta.push_back(LLP1_eta);
      tree_Hemi_LLP_phi.push_back(LLP1_phi);
      tree_Hemi_LLP_x.push_back(LLP1_x);
      tree_Hemi_LLP_y.push_back(LLP1_y);
      tree_Hemi_LLP_z.push_back(LLP1_z);
      tree_Hemi_LLP_dist.push_back(LLP1_dist);
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP1_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP1_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP1_z);
      tree_Hemi_Vtx_dr.push_back(TMath::Sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y) - TMath::Sqrt(LLP1_x*LLP1_x+LLP1_y*LLP1_y));

      dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
      ddok = TMath::Sqrt(dSV)/LLP1_dist;
      tree_Hemi_Vtx_dd.push_back( ddok );
      dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
      ddbad = TMath::Sqrt(dSV)/LLP2_dist;
      tree_Hemi_Vtx_ddbad.push_back( ddbad );
      if      ( Vtx_chi >= 0. && Vtx_chi < 10. && ddok  < 0.1 ) ping2 = 1;
      else if ( Vtx_chi >= 0. && Vtx_chi < 10. && ddbad < 0.1 ) ping2 = 2;
    }
    else {
      tree_Hemi_LLP_pt.push_back( LLP2_pt);
      tree_Hemi_LLP_eta.push_back(LLP2_eta);
      tree_Hemi_LLP_phi.push_back(LLP2_phi);
      tree_Hemi_LLP_x.push_back(LLP2_x);
      tree_Hemi_LLP_y.push_back(LLP2_y);
      tree_Hemi_LLP_z.push_back(LLP2_z);
      tree_Hemi_LLP_dist.push_back(LLP2_dist);
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP2_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP2_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP2_z);
      tree_Hemi_Vtx_dr.push_back(TMath::Sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y) - TMath::Sqrt(LLP2_x*LLP2_x+LLP2_y*LLP2_y));

      dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
      ddok = TMath::Sqrt(dSV)/LLP2_dist;
      tree_Hemi_Vtx_dd.push_back( ddok );
      dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
      ddbad = TMath::Sqrt(dSV)/LLP1_dist;
      tree_Hemi_Vtx_ddbad.push_back( ddbad );
      if      ( Vtx_chi >= 0. && Vtx_chi < 10. && ddok  < 0.1 ) ping2 = 2;
      else if ( Vtx_chi >= 0. && Vtx_chi < 10. && ddbad < 0.1 ) ping2 = 1;
    }
    tree_Hemi_LLP.push_back(iLLPrec2);
    
    tree_Hemi_dR12.push_back(dR_axis12);
    tree_Hemi_dR12.push_back(dR_axis12);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);

    bool ping_Hemi1 = false, ping_Hemi2 = false;
    if ( ping1 == iLLPrec1 ) ping_Hemi1 = true;
    if ( ping2 == iLLPrec2 ) ping_Hemi2 = true;
    if ( ping1 == iLLPrec2 && ping2 == iLLPrec1 ) {
      ping_Hemi1 = true;
      ping_Hemi2 = true;
    }
    if ( ping2 == 0 && ping1 == iLLPrec2 ) ping_Hemi1 = true;
    if ( ping1 == 0 && ping2 == iLLPrec1 ) ping_Hemi2 = true;
    tree_Hemi_LLP_ping.push_back( ping_Hemi1 );
    tree_Hemi_LLP_ping.push_back( ping_Hemi2 );
    int ping_event = 0;
    if      ( ping_Hemi1 && ping_Hemi2 ) ping_event = 2;
    else if ( ping_Hemi1 || ping_Hemi2 ) ping_event = 1;
    tree_event_LLP_ping.push_back( ping_event );
      
    // some informations for tracks in their hemisphere
    int ntrk10_vtx_hemi1 = 0., ntrk10_vtx_hemi2 = 0.;
    int ntrk20_vtx_hemi1 = 0., ntrk20_vtx_hemi2 = 0.;
    for (int counter_track = 0; counter_track < tree_nTracks; counter_track++) 
    {
      int hemi      = tree_track_Hemi[counter_track];
      double MVAval = tree_track_MVAval[counter_track];
      bool ping = false;
      Vtx_chi = -10.;
      float dist = -100.;
      if ( MVAval > bdtcut ) {
        if      ( hemi == 1 ) Vtx_chi = Vtx_chi1;
        else if ( hemi == 2 ) Vtx_chi = Vtx_chi2;
        if ( hemi == 1 && ping_Hemi1 ) ping = true;
        if ( hemi == 2 && ping_Hemi2 ) ping = true;
        if ( abs(Vtx_chi) < 10. ) {
          float x1 = tree_track_firstHit_x[counter_track] - tree_PV_x[0];
          float y1 = tree_track_firstHit_y[counter_track] - tree_PV_y[0];
          float z1 = tree_track_firstHit_z[counter_track] - tree_PV_z[0];
          float vtx_x = tree_Hemi_Vtx_x[hemi-1] - tree_PV_x[0];
          float vtx_y = tree_Hemi_Vtx_y[hemi-1] - tree_PV_y[0];
          float vtx_z = tree_Hemi_Vtx_z[hemi-1] - tree_PV_z[0];
          dist = TMath::Sqrt( (x1-vtx_x)*(x1-vtx_x) + (y1-vtx_y)*(y1-vtx_y) + (z1-vtx_z)*(z1-vtx_z) );
	  if ( x1*vtx_x + y1*vtx_y + z1*vtx_z < 0. ) dist = -dist;
	  if ( dist > 0. && dist < 10. && hemi == 1 ) ntrk10_vtx_hemi1++;
	  if ( dist > 0. && dist < 20. && hemi == 1 ) ntrk20_vtx_hemi1++;
	  if ( dist > 0. && dist < 10. && hemi == 2 ) ntrk10_vtx_hemi2++;
	  if ( dist > 0. && dist < 20. && hemi == 2 ) ntrk20_vtx_hemi2++;
	}
      }
      tree_track_Hemi_mva_NChi2.push_back(Vtx_chi);
      tree_track_Hemi_ping.push_back(ping);
      tree_track_Hemi_dFirstVtx.push_back( dist );
    } // End loop on tracks
    tree_Hemi_Vtx_ntrk10.push_back(ntrk10_vtx_hemi1);
    tree_Hemi_Vtx_ntrk10.push_back(ntrk10_vtx_hemi2);
    tree_Hemi_Vtx_ntrk20.push_back(ntrk20_vtx_hemi1);
    tree_Hemi_Vtx_ntrk20.push_back(ntrk20_vtx_hemi2);

  //Compute distance distributino bewteen vertices (including V0 candidates and Yconversion)


        for (unsigned int i= 0; i<tree_Hemi_Vtx_x.size();i++)
          {
            float vtx_x = tree_Hemi_Vtx_x[i];
            float vtx_y = tree_Hemi_Vtx_y[i];
            float vtx_z = tree_Hemi_Vtx_z[i];
            if (abs(tree_Hemi_Vtx_NChi2[i]) < 10. && tree_Hemi_nTrks[i]>1)
              {
                for (unsigned j = 0 ; j < tree_Hemi_Vtx_bkg_x.size() ; j++)
                  {
                    float dVtx_bkg = sqrt((vtx_x-tree_Hemi_Vtx_bkg_x[j])*(vtx_x-tree_Hemi_Vtx_bkg_x[j])+(vtx_y-tree_Hemi_Vtx_bkg_y[j])*(vtx_y-tree_Hemi_Vtx_bkg_y[j])+(vtx_z-tree_Hemi_Vtx_bkg_z[j])*(vtx_z-tree_Hemi_Vtx_bkg_z[j]));
                    tree_Hemi_Vtx_ddToBkg.push_back(dVtx_bkg);
                  }
              }
          }
      
  //////////////////////////////////
  // } // endif Filter
  
    // TFile* pFile2 = new TFile("triggerTest2.root","UPDATE");//RECREATE
  //Vertex Veto//


  } // endif Filter
  smalltree->Fill();
  // pFile->cd();
  // EffvsObs2->Write();
  // test->Write();

}


// ------------ method called once each job just before starting event loop  ------------
void
FlyingTopAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
FlyingTopAnalyzer::endJob()
{
          // std::cout<<"----------------Trigger Muon + dilepton-------------"<<std::endl;
          // for(unsigned int i = 0 ; i < tree_Trigger_Muon.size(); i++)
          //   {
          //     std::cout<<tree_Trigger_Muon[i]<<std::endl;
          //   }
          // std::cout<<"----------------Trigger Electron-------------"<<std::endl;
          // for(unsigned int i = 0 ; i < tree_Trigger_Ele.size(); i++)
          //   {
          //     std::cout<<tree_Trigger_Ele[i]<<std::endl;
          //   }
          // std::cout<<"----------------Trigger DoubleMu-------------"<<std::endl;
          // for(unsigned int i = 0 ; i < tree_Trigger_DoubleMu.size(); i++)
          //   {
          //     std::cout<<tree_Trigger_DoubleMu[i]<<std::endl;
          //   }
          // std::cout<<"----------------Trigger DoubleEle-------------"<<std::endl;
          // for(unsigned int i = 0 ; i < tree_Trigger_DoubleEle.size(); i++)
          //   {
          //     std::cout<<tree_Trigger_DoubleEle[i]<<std::endl;
          //   }
          //   std::cout<<"----------------Trigger Dimuon0-------------"<<std::endl;
          // for(unsigned int i = 0 ; i < tree_Trigger_Dimuon0.size(); i++)
          //   {
          //     std::cout<<tree_Trigger_Dimuon0[i]<<std::endl;
          //   }
          // std::cout<<"----------------Trigger PFMET-------------"<<std::endl;
          // for(unsigned int i = 0 ; i < tree_Trigger_PFMET.size(); i++)
          //   {
          //     std::cout<<tree_Trigger_PFMET[i]<<std::endl;
          //   }
          // std::cout<<"----------------Trigger HT-------------"<<std::endl;
          // for(unsigned int i = 0 ; i < tree_Trigger_HT.size(); i++)
          //   {
          //     std::cout<<tree_Trigger_HT[i]<<std::endl;
          //   }
          // std::cout<<"----------------Trigger AK4-------------"<<std::endl;
          // for(unsigned int i = 0 ; i < tree_Trigger_AK4.size(); i++)
          //   {
          //     std::cout<<tree_Trigger_AK4[i]<<std::endl;
          //   }
          // std::cout<<"----------------Trigger PFJet-------------"<<std::endl;
          // for(unsigned int i = 0 ; i < tree_Trigger_PFJet.size(); i++)
          //   {
          //     std::cout<<tree_Trigger_PFJet[i]<<std::endl;
          //   }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlyingTopAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FlyingTopAnalyzer);

void FlyingTopAnalyzer::clearVariables() {
    
//     tree_vtx_PosX.clear();
//     tree_vtx_PosY.clear();
//     tree_vtx_PosZ.clear();
//     tree_vtx_NChi2.clear();
//     tree_vtx_PosXError.clear();
//     tree_vtx_PosYError.clear();
//     tree_vtx_PosZError.clear();
    
// //     tree_trigger_names.clear();
// //     tree_trigger_bits.clear();
//     tree_trigger_size.clear();
//     tree_passesTrigger.clear();
//     tree_passesTriggerName.clear();
//     tree_Trigger_Muon.clear();//+ dilepton channel emu
//     tree_Trigger_Ele.clear();
//     tree_Trigger_DoubleMu.clear();
//     tree_Trigger_DoubleEle.clear();
//     tree_Trigger_Dimuon0.clear();
//     tree_Trigger_PFMET.clear();
//     tree_Trigger_HT.clear();
//     tree_Trigger_AK4.clear();
//     tree_Trigger_PFJet.clear();
//     tree_Trigger_DoublePFJets.clear();
//     tree_Trigger_DiPFJet.clear();
//     tree_Trigger_QuadPFJet.clear();
//     tree_Trigger_BTagMu.clear();
    
    tree_PV_x.clear(); // l'index 0 donne le PV!
    tree_PV_y.clear();
    tree_PV_z.clear();
    tree_PV_ez.clear();
    tree_PV_NChi2.clear();
    tree_PV_ndf.clear();

      
    tree_SV_x.clear();
    tree_SV_y.clear();
    tree_SV_z.clear();
    tree_SV_r.clear();
    tree_SV_NChi2.clear();
    tree_SV_ndf.clear();
    tree_SV_nDaughters.clear();
    tree_SV_mass.clear();
    tree_SV_pdgid.clear();
    tree_SV_pt.clear();
    tree_SV_eta.clear();
    tree_SV_phi.clear();
    tree_SV_IsConvertedPhoton.clear();
//$$$$    tree_SV_daughter_pdgid.clear();

    tree_K0_x.clear();
    tree_K0_y.clear();
    tree_K0_z.clear();
    tree_K0_r.clear();
    tree_K0_NChi2.clear();
    tree_K0_ndf.clear();
    tree_K0_nDaughters.clear();
    tree_K0_mass.clear();
    tree_K0_pt.clear();
    tree_K0_eta.clear();
    tree_K0_phi.clear();
    tree_K0_daughter_pt.clear();
    tree_K0_daughter_eta.clear();
    tree_K0_daughter_phi.clear();

    tree_L0_x.clear();
    tree_L0_y.clear();
    tree_L0_z.clear();
    tree_L0_r.clear();
    tree_L0_NChi2.clear();
    tree_L0_ndf.clear();
    tree_L0_nDaughters.clear();
    tree_L0_mass.clear();
    tree_L0_pt.clear();
    tree_L0_eta.clear();
    tree_L0_phi.clear();
    tree_L0_daughter_pt.clear();
    tree_L0_daughter_eta.clear();
    tree_L0_daughter_phi.clear();
    

    tree_Yc_x.clear(); 
    tree_Yc_y.clear();
    tree_Yc_z.clear();
    tree_Yc_r.clear();
    tree_Yc_Vtx_NChi2.clear();
    tree_Yc_Vtx_ndf.clear();
    tree_Yc_mass.clear();
    tree_Yc_tracks_pt.clear();
    tree_Yc_tracks_eta.clear();
    tree_Yc_tracks_phi.clear();
    tree_Yc_tracks_sum3p.clear();
    tree_Yc_tracks_InPosx.clear();
    tree_Yc_tracks_InPosy.clear();
    tree_Yc_tracks_InPosz.clear();
    tree_Yc_tracks_InPx.clear();
    tree_Yc_tracks_InPy.clear();
    tree_Yc_tracks_InPz.clear();
    tree_Yc_tracks_OutPx.clear();
    tree_Yc_tracks_OutPy.clear();
    tree_Yc_tracks_OutPz.clear();
    tree_Yc_tracks_pt.clear();
    tree_Yc_tracks_eta.clear();
    tree_Yc_tracks_phi.clear();

    tree_track_bkg_pt.clear();
    tree_track_bkg_eta.clear();
    tree_track_bkg_phi.clear();
    tree_track_bkg_charge.clear();
    tree_track_bkg_source.clear();
    tree_track_dpt.clear();
    tree_track_deta.clear();
    tree_track_dphi.clear();
    tree_track_bkg.clear();
    tree_Hemi_Vtx_bkg_x.clear();
    tree_Hemi_Vtx_bkg_y.clear();
    tree_Hemi_Vtx_bkg_z.clear();

    tree_jet_E.clear();
    tree_jet_pt.clear();
    tree_jet_eta.clear();
    tree_jet_phi.clear();
    
    tree_electron_pt.clear();
    tree_electron_eta.clear();
    tree_electron_phi.clear();
    tree_electron_x.clear();
    tree_electron_y.clear();
    tree_electron_z.clear();
    tree_electron_energy.clear();
    tree_electron_charge.clear();
    tree_electron_isoR4.clear(); 

    tree_muon_pt.clear();
    tree_muon_eta.clear();
    tree_muon_phi.clear();
    tree_muon_x.clear();
    tree_muon_y.clear();
    tree_muon_z.clear();
    tree_muon_energy.clear();
    tree_muon_dxy.clear();
    tree_muon_dxyError.clear();
    tree_muon_dz.clear();
    tree_muon_dzError.clear();
    tree_muon_charge.clear();
    tree_muon_isLoose.clear();
    tree_muon_isTight.clear();
    tree_muon_isGlobal.clear();
    tree_muon_isoR3.clear();
//     tree_passesTrkPtr.clear();

    tree_track_ipc.clear();
    tree_track_lost.clear();
    tree_track_pt.clear();
    tree_track_eta.clear();
    tree_track_phi.clear();
    tree_track_charge.clear();
    tree_track_NChi2.clear();
    tree_track_isHighPurity.clear();
    tree_track_dxy.clear();
    tree_track_dxyError.clear();
    tree_track_drSig.clear();
    tree_track_dz.clear();
    tree_track_dzError.clear();
    tree_track_algo.clear();
    tree_track_nHit.clear();
    tree_track_nHitPixel.clear();
    tree_track_nHitTIB.clear();
    tree_track_nHitTID.clear();
    tree_track_nHitTOB.clear();
    tree_track_nHitTEC.clear();
    tree_track_nHitPXB.clear();
    tree_track_nHitPXF.clear();
    tree_track_isHitPixel.clear();
    tree_track_nLayers.clear();
    tree_track_nLayersPixel.clear();
    tree_track_x.clear();
    tree_track_y.clear();
    tree_track_z.clear();
    tree_track_firstHit.clear();
    tree_track_firstHit_x.clear();
    tree_track_firstHit_y.clear();
    tree_track_firstHit_z.clear();
    tree_track_region.clear();
    tree_track_iJet.clear();
    tree_track_ntrk10.clear();
    tree_track_ntrk20.clear();
    tree_track_ntrk30.clear();
    tree_track_ntrk40.clear();
    tree_track_MVAval.clear();
    
    tree_track_Hemi.clear();
    tree_track_Hemi_dR.clear();
    tree_track_Hemi_mva_NChi2.clear();
    tree_track_Hemi_ping.clear();
    tree_track_Hemi_dFirstVtx.clear();
    tree_track_Hemi_LLP.clear();
    
    tree_track_sim_LLP.clear();
    tree_track_sim_isFromB.clear();
    tree_track_sim_isFromC.clear();
    tree_track_sim_pt.clear();
    tree_track_sim_eta.clear();
    tree_track_sim_phi.clear();
    tree_track_sim_charge.clear();
    tree_track_sim_pdgId.clear();
    tree_track_sim_mass.clear();
    tree_track_sim_x.clear();
    tree_track_sim_y.clear();
    tree_track_sim_z.clear();
    tree_track_sim_dFirstGen.clear();
    tree_track_sim_LLP_r.clear();
    tree_track_sim_LLP_z.clear();

    tree_genParticle_pt.clear();
    tree_genParticle_eta.clear();
    tree_genParticle_phi.clear();
    tree_genParticle_charge.clear();
    tree_genParticle_pdgId.clear();
    tree_genParticle_mass.clear();
    tree_genParticle_x.clear();
    tree_genParticle_y.clear();
    tree_genParticle_z.clear();
    tree_genParticle_statusCode.clear();
    tree_genParticle_mother_pdgId.clear();
    tree_genParticle_LLP.clear();

    tree_genPackPart_pt.clear();
    tree_genPackPart_eta.clear();
    tree_genPackPart_phi.clear();
    tree_genPackPart_charge.clear();
    tree_genPackPart_pdgId.clear();
    tree_genPackPart_mass.clear();

    tree_genPackPart_x.clear();
    tree_genPackPart_y.clear();
    tree_genPackPart_z.clear();
    tree_genPackPart_mother_pdgId.clear();
    tree_genPackPart_isFromB.clear();
    tree_genPackPart_isFromC.clear();

    tree_genFromLLP_LLP.clear();
    tree_genFromLLP_pt.clear();
    tree_genFromLLP_eta.clear();
    tree_genFromLLP_phi.clear();
    tree_genFromLLP_charge.clear();
    tree_genFromLLP_pdgId.clear();
    tree_genFromLLP_mass.clear();
    tree_genFromLLP_x.clear();
    tree_genFromLLP_y.clear();
    tree_genFromLLP_z.clear();
    tree_genFromLLP_mother_pdgId.clear();
    tree_genFromLLP_isFromB.clear();
    tree_genFromLLP_isFromC.clear();

    tree_genAxis_dRneuneu.clear();

    tree_genFromC_pt.clear();
    tree_genFromC_eta.clear();
    tree_genFromC_phi.clear();
    tree_genFromC_charge.clear();
    tree_genFromC_pdgId.clear();
    tree_genFromC_x.clear();
    tree_genFromC_y.clear();
    tree_genFromC_z.clear();
    tree_genFromC_mother_pdgId.clear();

    tree_genFromB_pt.clear();
    tree_genFromB_eta.clear();
    tree_genFromB_phi.clear();
    tree_genFromB_charge.clear();
    tree_genFromB_pdgId.clear();
    tree_genFromB_x.clear();
    tree_genFromB_y.clear();
    tree_genFromB_z.clear();
    tree_genFromB_mother_pdgId.clear();

    tree_genJet_pt.clear();
    tree_genJet_eta.clear();
    tree_genJet_phi.clear();
    tree_genJet_mass.clear();
    tree_genJet_energy.clear();
    
    tree_LLP.clear();
    tree_LLP_pt.clear();
    tree_LLP_eta.clear();
    tree_LLP_phi.clear();
    tree_LLP_x.clear();
    tree_LLP_y.clear();
    tree_LLP_z.clear();
    tree_LLP_dist.clear();
    tree_LLP_nTrks.clear();
    tree_LLP_Vtx_NChi2.clear();
    tree_LLP_Vtx_nTrks.clear();
    tree_LLP_Vtx_dx.clear();
    tree_LLP_Vtx_dy.clear();
    tree_LLP_Vtx_dz.clear();
    tree_LLP_Vtx_dist.clear();
    tree_LLP_Vtx_dd.clear();

    tree_Hemi_Vtx_Layer.clear();
    tree_Hemi_Vtx_evt2vtx.clear();
    tree_Hemi.clear();
    tree_Hemi_njet.clear();
    tree_Hemi_eta.clear();
    tree_Hemi_phi.clear();
    tree_Hemi_dR.clear();
    tree_Hemi_nTrks.clear();
    tree_Hemi_nTrks_sig.clear();
    tree_Hemi_nTrks_bad.clear();
    tree_Hemi_LLP.clear();
    tree_Hemi_LLP_pt.clear();
    tree_Hemi_LLP_eta.clear();
    tree_Hemi_LLP_phi.clear();
    tree_Hemi_LLP_dist.clear();
    tree_Hemi_LLP_x.clear();
    tree_Hemi_LLP_y.clear();
    tree_Hemi_LLP_z.clear();
    tree_Hemi_Vtx_step.clear();
    tree_Hemi_Vtx_NChi2.clear();
    tree_Hemi_Vtx_nTrks.clear();
    tree_Hemi_Vtx_nTrks_sig.clear();
    tree_Hemi_Vtx_nTrks_bad.clear();
    tree_Hemi_Vtx_x.clear();
    tree_Hemi_Vtx_y.clear();
    tree_Hemi_Vtx_r.clear();
    tree_Hemi_Vtx_z.clear();
    tree_Hemi_Vtx_dist.clear();
    tree_Hemi_Vtx_dx.clear();
    tree_Hemi_Vtx_dy.clear();
    tree_Hemi_Vtx_dz.clear();
    tree_Hemi_Vtx_dr.clear();
    tree_Hemi_Vtx_dd.clear();
    tree_Hemi_dR12.clear();
    tree_Hemi_LLP_dR12.clear();
    tree_Hemi_Vtx_ddbad.clear();
    tree_Hemi_Vtx_ntrk10.clear();
    tree_Hemi_Vtx_ntrk20.clear();
    tree_Hemi_Vtx_ddToBkg.clear();
//     tree_Hemi_Vtx_trackWeight.clear();
    tree_Hemi_LLP_ping.clear();
    tree_event_LLP_ping.clear();
    tree_Hemi_Vtx_K0.clear();
    tree_Hemi_Vtx_L0.clear();
    tree_Hemi_Vtx_V0.clear();
    tree_Hemi_Vtx_Yc.clear();
}
