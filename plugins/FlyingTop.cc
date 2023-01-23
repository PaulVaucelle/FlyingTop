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
    int runNumber, eventNumber, lumiBlock;
    int  tree_NbrOfZCand;
    bool tree_passesHTFilter;
    int  tree_nTracks, tree_nLostTracks; 
    int  nBC = 0, tree_nFromC = 0, tree_nFromB = 0; 
    int nEvent;
    
    float LLP1_pt, LLP1_eta, LLP1_phi, LLP2_pt, LLP2_eta, LLP2_phi;
    float LLP1_x, LLP1_y, LLP1_z, LLP2_x, LLP2_y, LLP2_z;
    float LLP1_dist, LLP2_dist;
    int   LLP1_nTrks = 0, LLP2_nTrks = 0;

    //-----------------------
    // trigger variable
//     std::vector<string > tree_trigger_names;
//     std::vector<bool >   tree_trigger_bits;
   //trig
    std::vector<int>     tree_trigger_size;
    std::vector<int>        tree_passesTrigger;
    std::vector<string> tree_passesTriggerName;
    //trig
    
//$$
//The BDT variables are declared here to reduce computation time
    float pt, eta, NChi, nhits, ntrk10, drSig, isinjet;
    TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );

    bool NewCovMat = true;//Allow for Covariance Matrix correction due to the MiniAOD dataformat apporixmation
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
    std::vector< double > tree_track_MVAval;

    std::vector< int >    tree_track_Hemi;
    std::vector< double > tree_track_Hemi_dR;
    std::vector< double > tree_track_Hemi_mva_NChi2;
//$$
    std::vector< bool >   tree_track_Hemi_ping;
    std::vector< float >  tree_track_Hemi_dFirstVtx;
//$$
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

    std::vector< float > tree_genPackPart_pt;
    std::vector< float > tree_genPackPart_eta;
    std::vector< float > tree_genPackPart_phi;
    std::vector< float > tree_genPackPart_charge;
    std::vector< int >   tree_genPackPart_pdgId;
    std::vector< float > tree_genPackPart_mass;
    std::vector< int >   tree_genPackPart_mother_pdgId;

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
    std::vector< float > tree_LLP_Vtx_trackWeight;
    
    //-----------------------
    //Analysis with the two hemispheres
    //-----------------------
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
    std::vector< float > tree_Hemi_Vtx_NChi2;
    std::vector< int >   tree_Hemi_Vtx_nTrks;
    std::vector< int >   tree_Hemi_Vtx_nTrks_sig;
    std::vector< int >   tree_Hemi_Vtx_nTrks_bad;
    std::vector< float > tree_Hemi_Vtx_x;
    std::vector< float > tree_Hemi_Vtx_y;
    std::vector< float > tree_Hemi_Vtx_z;
    std::vector< float > tree_Hemi_Vtx_dist;
    std::vector< float > tree_Hemi_Vtx_dx;
    std::vector< float > tree_Hemi_Vtx_dy;
    std::vector< float > tree_Hemi_Vtx_dz;
    std::vector< float > tree_Hemi_Vtx_dd;
    std::vector< float > tree_Hemi_dR12;
    std::vector< float > tree_Hemi_LLP_dR12;
    std::vector< float > tree_Hemi_Vtx_UpdatedChi2;
    std::vector< int >   tree_Hemi_Vtx_ntrk;
    std::vector< float > tree_Hemi_Vtx_TotalChi2;
    std::vector< float > tree_Hemi_Vtx_DOF;
//$$
    std::vector< float > tree_Hemi_Vtx_ddbad;
    std::vector< int >   tree_Hemi_Vtx_ntrk10;
    std::vector< int >   tree_Hemi_Vtx_ntrk20;
    std::vector< float > tree_Hemi_Vtx_trackWeight;
    std::vector< bool >  tree_Hemi_LLP_ping;
    std::vector< int >   tree_event_LLP_ping;
//$$
    std::vector< float > tree_Hemi_Vtx_MeanBdtVal;
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
    
//     smalltree->Branch("tree_vtx_PosX", &tree_vtx_PosX);
//     smalltree->Branch("tree_vtx_PosY", &tree_vtx_PosY);
//     smalltree->Branch("tree_vtx_PosZ", &tree_vtx_PosZ);
//     smalltree->Branch("tree_vtx_NChi2", &tree_vtx_NChi2);
//     smalltree->Branch("tree_vtx_PosXError", &tree_vtx_PosXError);
//     smalltree->Branch("tree_vtx_PosYError", &tree_vtx_PosYError);
//     smalltree->Branch("tree_vtx_PosZError", &tree_vtx_PosZError);
    
    // trigger info
//     smalltree->Branch("tree_trigger_names", &tree_trigger_names);
//     smalltree->Branch("tree_trigger_bits",  &tree_trigger_bits);
        //trig
    smalltree->Branch("tree_trigger_size", &tree_trigger_size);
    smalltree->Branch("tree_passesTrigger", &tree_passesTrigger);
    smalltree->Branch("tree_passesTriggerName", &tree_passesTriggerName);
    //trig

    smalltree->Branch("tree_NbrOfZCand",  &tree_NbrOfZCand,  "tree_NbrOfZCand/I");
    smalltree->Branch("tree_passesHTFilter", &tree_passesHTFilter);
    
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
    
    // electrons info
    smalltree->Branch("tree_electron_pt"  ,   &tree_electron_pt);
    smalltree->Branch("tree_electron_eta" ,   &tree_electron_eta);
    smalltree->Branch("tree_electron_phi" ,   &tree_electron_phi);
    smalltree->Branch("tree_electron_x"  ,    &tree_electron_x);
    smalltree->Branch("tree_electron_y" ,     &tree_electron_y);
    smalltree->Branch("tree_electron_z" ,     &tree_electron_z);
    smalltree->Branch("tree_electron_energy", &tree_electron_energy);
    smalltree->Branch("tree_electron_charge", &tree_electron_charge);
    
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
    smalltree->Branch("tree_track_MVAval",         &tree_track_MVAval);
    

    smalltree->Branch("tree_track_Hemi",           &tree_track_Hemi);
    smalltree->Branch("tree_track_Hemi_dR",        &tree_track_Hemi_dR);
    smalltree->Branch("tree_track_Hemi_mva_NChi2", &tree_track_Hemi_mva_NChi2);
//$$
    smalltree->Branch("tree_track_Hemi_ping",      &tree_track_Hemi_ping);
    smalltree->Branch("tree_track_Hemi_dFirstVtx", &tree_track_Hemi_dFirstVtx);
//$$
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

    smalltree->Branch("tree_genPackPart_pt"  ,          &tree_genPackPart_pt);
    smalltree->Branch("tree_genPackPart_eta" ,          &tree_genPackPart_eta);
    smalltree->Branch("tree_genPackPart_phi" ,          &tree_genPackPart_phi);
    smalltree->Branch("tree_genPackPart_charge" ,       &tree_genPackPart_charge);
    smalltree->Branch("tree_genPackPart_pdgId" ,        &tree_genPackPart_pdgId);
    smalltree->Branch("tree_genPackPart_mass" ,         &tree_genPackPart_mass);
    smalltree->Branch("tree_genPackPart_mother_pdgId" , &tree_genPackPart_mother_pdgId);

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
    smalltree->Branch("tree_LLP_Vtx_trackWeight", &tree_LLP_Vtx_trackWeight);

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
    smalltree->Branch("tree_Hemi_Vtx_NChi2", &tree_Hemi_Vtx_NChi2);
    smalltree->Branch("tree_Hemi_Vtx_nTrks", &tree_Hemi_Vtx_nTrks);
    smalltree->Branch("tree_Hemi_Vtx_nTrks_sig", &tree_Hemi_Vtx_nTrks_sig);
    smalltree->Branch("tree_Hemi_Vtx_nTrks_bad", &tree_Hemi_Vtx_nTrks_bad);
    smalltree->Branch("tree_Hemi_Vtx_x",     &tree_Hemi_Vtx_x);
    smalltree->Branch("tree_Hemi_Vtx_y",     &tree_Hemi_Vtx_y);
    smalltree->Branch("tree_Hemi_Vtx_z",     &tree_Hemi_Vtx_z);
    smalltree->Branch("tree_Hemi_Vtx_dist",  &tree_Hemi_Vtx_dist);
    smalltree->Branch("tree_Hemi_Vtx_dx",    &tree_Hemi_Vtx_dx);
    smalltree->Branch("tree_Hemi_Vtx_dy",    &tree_Hemi_Vtx_dy);
    smalltree->Branch("tree_Hemi_Vtx_dz",    &tree_Hemi_Vtx_dz);
    smalltree->Branch("tree_Hemi_Vtx_dd",    &tree_Hemi_Vtx_dd);
    smalltree->Branch("tree_Hemi_dR12",      &tree_Hemi_dR12);
    smalltree->Branch("tree_Hemi_LLP_dR12",  &tree_Hemi_LLP_dR12);
    smalltree->Branch("tree_Hemi_Vtx_UpdatedChi2", &tree_Hemi_Vtx_UpdatedChi2);
    smalltree->Branch("tree_Hemi_Vtx_ntrk",        &tree_Hemi_Vtx_ntrk);
    smalltree->Branch("tree_Hemi_Vtx_TotalChi2",   &tree_Hemi_Vtx_TotalChi2);
    smalltree->Branch("tree_Hemi_Vtx_DOF",         &tree_Hemi_Vtx_DOF);
//$$
    smalltree->Branch("tree_Hemi_Vtx_ddbad", &tree_Hemi_Vtx_ddbad);
    smalltree->Branch("tree_Hemi_Vtx_ntrk10",&tree_Hemi_Vtx_ntrk10);
    smalltree->Branch("tree_Hemi_Vtx_ntrk20",&tree_Hemi_Vtx_ntrk20);
    smalltree->Branch("tree_Hemi_Vtx_trackWeight", &tree_Hemi_Vtx_trackWeight);
    smalltree->Branch("tree_Hemi_LLP_ping",  &tree_Hemi_LLP_ping);
    smalltree->Branch("tree_event_LLP_ping", &tree_event_LLP_ping);
    smalltree->Branch("tree_Hemi_Vtx_MeanBdtVal", &tree_Hemi_Vtx_MeanBdtVal);
//$$

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
 
 //trig
std::vector<std::string> VTrigger;
std::vector<bool> PassTrigger;
  // ######################
  const edm::Handle<edm::TriggerResults> triggerH = iEvent.getHandle(triggerResultsToken_);
  const auto triggerNames = iEvent.triggerNames(*triggerH);
  for (unsigned int i = 0; i < triggerH->size(); i++) {
  VTrigger.push_back(triggerNames.triggerName(i));
  //  std::cout << "Trigger name : " << triggerNames.triggerName(i) << " , trigger pass ? " << triggerH->accept(i)<<std::endl;
  if (triggerH->accept(i))
    {
      // std::cout<<" trigger passed : "<<triggerNames.triggerName(i)<<std::endl;
      tree_passesTrigger.push_back(i);
      tree_passesTriggerName.push_back(triggerNames.triggerName(i));
    }
}
tree_trigger_size.push_back(VTrigger.size());//should be the same for each evnet!!



//trig

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

      // smuon
      if ( genIt.pdgId() == 1000013 ) {
	tree_GenPVx = genIt.vx();
	tree_GenPVy = genIt.vy();
	tree_GenPVz = genIt.vz();
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

  //////////////////////////////////
  //////////////////////////////////
  ////////   Electrons   ///////////
  //////////////////////////////////
  //////////////////////////////////
  //One could consider the emu channel, or even ee channel
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
  }

  //////////////////////////////////
  //////////////////////////////////
  ///////////   Muons   ////////////
  //////////////////////////////////
  //////////////////////////////////
  
  int nmu = 0;
  for (const pat::Muon &mu : *muons)
  {
  if ( mu.pt() < 3. ) continue;
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
    nmu++;
  }
    
  int imu1 = -1, imu2 = -1;
  float mupt1, mueta1, muphi1, mupt2, mueta2, muphi2;
  float mu_mass = 0.1057;
  TLorentzVector v1, v2, v;
  tree_Mmumu = 0.;
  
  for ( int mu=0; mu<nmu; mu++)
  { 
  if ( !tree_muon_isGlobal[mu] ) continue;
    mupt1  = tree_muon_pt[mu];
  if ( mupt1 < 10. ) continue; // Zmu filter
//$$  if ( abs(tree_muon_dxy[mu]) > 0.1 || abs(tree_muon_dz[mu]) > 0.2 ) continue; // muons closed to PV
    mueta1 = tree_muon_eta[mu];
    muphi1 = tree_muon_phi[mu];
    v1.SetPtEtaPhiM(mupt1,mueta1,muphi1,mu_mass);
    for ( int mu2=mu+1; mu2<nmu; mu2++) 
    {	    
    if ( !tree_muon_isGlobal[mu2] ) continue;
    if ( tree_muon_charge[mu] == tree_muon_charge[mu2] ) continue;
//$$    if ( abs(tree_muon_dxy[mu2]) > 0.1 || abs(tree_muon_dz[mu2]) > 0.2 ) continue;
      mupt2  = tree_muon_pt[mu2];
    if ( mupt2 < 10. ) continue;
    if ( mupt1 < 28. && mupt2 < 28. ) continue; // Zmu Filter
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
  }

  if ( tree_muon_pt[imu2] > tree_muon_pt[imu1] ) {
    int imu0 = imu2;
    imu2 = imu1; // muons reco with imu1 having the highest pt
    imu1 = imu0;
  }

  //////////////////////////////////
  //////////////////////////////////
  //////// HT FILTER CHECK /////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_NbrOfZCand = 0;
  tree_passesHTFilter = false;
  tree_nTracks = 0;
  tree_nLostTracks = 0;

  if ( tree_Mmumu > 60. )                  tree_NbrOfZCand = 1;
  if ( tree_Mmumu > 60. && HT_val > 180. ) tree_passesHTFilter = true;
  

  edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder); // Asking for reco collection of PV..
  vector<reco::TransientTrack> BestTracks;
  std::vector<std::pair<uint16_t,float> > Players;
  int count =0;
  std::map<size_t , int > trackToAK4SlimmedJetMap;

//$$ // if ( tree_passesHTFilter ) {


  //////////////////////////////////
  //////////////////////////////////
  //////////   Tracks   ////////////
  //////////////////////////////////
  //////////////////////////////////
  
    float pt_Cut = 1.;
    float NChi2_Cut = 5.;
    float drSig_Cut = 5.;

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
      
      // ----------------------------------------------------
      const reco::Track& tk_temp = *trackPcPtr; 
      reco::Track tk;
      reco::TrackBase::CovarianceMatrix m_ = tk_temp.covariance(); //math::Error<5>::type

      // In MiniAOD dataTier, the covariance matrix is approximated => Drop of 10% in the vertex reco efficiency
      // A Correction is given  on the CMSSW_13 release (that is written here)
      double det=0;
      //Covariance Matrix Correction//
      bool notPosDef = (!(m_).Sub<AlgebraicSymMatrix22>(0, 0).Det(det) || det < 0) ||
                       (!(m_).Sub<AlgebraicSymMatrix33>(0, 0).Det(det) || det < 0) ||
                       (!(m_).Sub<AlgebraicSymMatrix44>(0, 0).Det(det) || det < 0) || 
		       (!(m_).Det(det) || det < 0);
      if ( notPosDef && NewCovMat ) {
        reco::TrackBase::CovarianceMatrix m(m_);
        //if not positive-definite, alter values to allow for pos-def
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
          for (int i = 0; i < 5; i++) m(i, i) += delta - minEigenValue;
        }
        // make a track object with pos def covariance matrix  
        //If a correction is applied, we consider the track object with the corrections
        tk = reco::Track(tk_temp.normalizedChi2() * tk_temp.ndof(),tk_temp.ndof(),tk_temp.TrackBase::referencePoint(),tk_temp.momentum(),tk_temp.charge(),m,reco::TrackBase::undefAlgorithm,reco::TrackBase::loose);
	tk_nHit   = tk_temp.hitPattern().numberOfValidHits();
        tk_charge = tk_temp.charge();
        tk_pt  = tk_temp.pt();
        tk_eta = tk_temp.eta();
        tk_phi = tk_temp.phi();
        tk_NChi2 = tk_temp.normalizedChi2();
        if ( tk_temp.dxyError() > 0 ) 
          tk_drSig = abs(tk_temp.dxy(PV.position())) / tk_temp.dxyError(); // from Paul
        tk_dxy = tk_temp.dxy(PV.position());
        tk_dxyError = tk_temp.dxyError();
        tk_dz = tk_temp.dz(PV.position());
        tk_dzError = tk_temp.dzError();
        tk_vx = tk_temp.vx();
        tk_vy = tk_temp.vy();
        tk_vz = tk_temp.vz();
        tk_px = tk_temp.px();
        tk_py = tk_temp.py();
        tk_pz = tk_temp.pz();
        tk_HitPattern = tk_temp.hitPattern();

      }
      else //If no correction is applied, we consider the actual pfcandidate as track object 
      {
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
      }
  //-------------------end covariance matrix correction-----------//

    if ( tk_nHit == 0 ) continue;
    if ( tk_charge == 0 ) continue;
//$$
    if ( !(tk_pt > pt_Cut && tk_NChi2 < NChi2_Cut && tk_drSig > drSig_Cut) ) continue; // preselection
//$$
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
      tree_track_algo.push_back         (tk_temp.algo());
      tree_track_isHighPurity.push_back (static_cast<int>(tk_temp.quality(reco::TrackBase::highPurity)));
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

      float dSign = 1.;
      if ( kmatch >= 0 &&
           xFirst*tree_genFromLLP_x[kmatch]+yFirst*tree_genFromLLP_y[kmatch]+zFirst*tree_genFromLLP_z[kmatch] < 0. ) dSign = -1.;
      tree_track_sim_dFirstGen.push_back( TMath::Sqrt(dFirstGenMin)*dSign );
      float track_sim_LLP_r = 0.;
      float track_sim_LLP_z = 0.;
      if ( track_sim_LLP == 1 ) {
        track_sim_LLP_r = TMath::Sqrt( LLP1_x*LLP1_x + LLP1_y*LLP1_y );
        track_sim_LLP_z = abs(LLP1_z);
      }
      if ( track_sim_LLP == 2 ) {
        track_sim_LLP_r = TMath::Sqrt( LLP2_x*LLP2_x + LLP2_y*LLP2_y );
        track_sim_LLP_z = abs(LLP2_z);
      }
      tree_track_sim_LLP_r.push_back( track_sim_LLP_r );
      tree_track_sim_LLP_z.push_back( track_sim_LLP_z );

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

    vector<reco::TransientTrack> displacedTracks_llp1_mva, displacedTracks_llp2_mva; // Control Tracks
    vector<reco::TransientTrack> displacedTracks_Hemi1_mva, displacedTracks_Hemi2_mva; // Tracks selected wrt the hemisphere
    std::pair<vector<reco::TransientTrack>, std::vector<float>> PairSortedHemi1;
    std::pair<vector<reco::TransientTrack>, std::vector<float>> PairSortedHemi2;

    int jet; /*!*/
    float ntrk20, ntrk30; /*!*/
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
    double bdtcut = -0.1456; // for TMVAClassification_BDTG50cm_HighPurity.weights.xml BDTrecohp
//     double bdtcut = -0.1083; // for TMVAClassification_BDTG_FromBC.weights.xml from BDTminipf
//     double bdtcut = -0.0067; // for TMVAClassification_BDTG50cm_sansntrk10_avecHP.weights.xml BDTrecohpsansntrk10
//     double bdtcut = -10.; // no BDT cut
//$$

    //---------------------------//
    // if (tree_passesHTFilter){

    for (int counter_track = 0; counter_track < tree_nTracks; counter_track++) 
    {
      unsigned int ipc = tree_track_ipc[counter_track];
      pat::PackedCandidateRef pcref = MINIgeneralTracks[ipc];
      const reco::Track *trackPcPtr = pcref->bestTrack();

      reco::Track tk ;
      const reco::Track& tk_temp = *trackPcPtr;
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
      }
      else tk = *trackPcPtr;


        //---------------end of Wip-----------------------//
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
      isinjet = 0.;
      ntrk20 = 0;
      ntrk30 = 0;
      bdtval = -10.;
      dR = -1.;
      int tracks_axis = 0; // flag to check which axis is the closest from the track

      jet = tree_track_iJet[counter_track];
      if ( jet >= 0 ) isinjet = 1.; /*!*/
      int isFromLLP = tree_track_sim_LLP[counter_track];

      //check the dR between the tracks and the second axis (without any selection on the tracks)
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
        if ( dist < 10. )	   {ntrk10++;} // les sctocker les 3 , on teste sur une seule couche quand on regarde vers l'avant
        if ( dist < 20. )	   {ntrk20++;}
        if ( dist < 30. )	   {ntrk30++;}
      }  // end Loop on other Tracks

      if ( dR < dRcut_tracks ) 
      {
        if ( isFromLLP == 1 ) LLP1_nTrks++;
        if ( isFromLLP == 2 ) LLP2_nTrks++;

        bdtval = reader->EvaluateMVA( "BDTG" ); //default value = -10 (no -10 observed and -999 comes from EvaluateMVA)

        if ( tracks_axis == 1 ) {
          nTrks_axis1++;
          if ( isFromLLP == iLLPrec1 ) nTrks_axis1_sig++;
          else if ( isFromLLP >= 1 )   nTrks_axis1_bad++;
        }
      
        if ( tracks_axis == 2 ) {
          nTrks_axis2++;
          if ( isFromLLP == iLLPrec2 ) nTrks_axis2_sig++;
          else if ( isFromLLP >= 1 )   nTrks_axis2_bad++;
        }
      
        if ( bdtval > bdtcut ) {
          ////--------------Control tracks-----------------////
          if ( isFromLLP == 1 )
          {
            displacedTracks_llp1_mva.push_back(theTransientTrackBuilder->build(trackPcPtr));
          }
          if ( isFromLLP == 2 )
          {
            displacedTracks_llp2_mva.push_back(theTransientTrackBuilder->build(trackPcPtr));
          }

          if ( tracks_axis == 1 )
          {
            displacedTracks_Hemi1_mva.push_back(theTransientTrackBuilder->build(&tk));//&tk contains the covariance matrix corrections if needed
            PairSortedHemi1.first.push_back(theTransientTrackBuilder->build(&tk));//&tk contains the covariance matrix corrections if needed
            PairSortedHemi1.second.push_back(bdtval);
            if ( isFromLLP == iLLPrec1 ) nTrks_axis1_sig_mva++;
            else if ( isFromLLP >= 1 )   nTrks_axis1_bad_mva++;
          }
 
          if ( tracks_axis == 2 )
          {
            displacedTracks_Hemi2_mva.push_back(theTransientTrackBuilder->build(&tk));//&tk contains the covariance matrix corrections if needed
            PairSortedHemi2.first.push_back(theTransientTrackBuilder->build(&tk));//&tk contains the covariance matrix corrections if needed
            PairSortedHemi2.second.push_back(bdtval);
            if ( isFromLLP == iLLPrec2 ) nTrks_axis2_sig_mva++;
            else if ( isFromLLP >= 1 )   nTrks_axis2_bad_mva++;
          }
        }
      }
      
      tree_track_ntrk10.push_back(ntrk10);
      tree_track_ntrk20.push_back(ntrk20);
      tree_track_ntrk30.push_back(ntrk30);
      tree_track_MVAval.push_back(bdtval);
      tree_track_Hemi.push_back(tracks_axis);
      tree_track_Hemi_dR.push_back(dR);
      if      ( tracks_axis == 1 ) tree_track_Hemi_LLP.push_back(iLLPrec1);
      else if ( tracks_axis == 2 ) tree_track_Hemi_LLP.push_back(iLLPrec2);
      else		           tree_track_Hemi_LLP.push_back(0);
      
    } //End loop on all the tracks

    ////////////////////////////////////
    // Sort Tracks by their BDT value //
    ////////////////////////////////////
    //--Hemi1--//
    if ( displacedTracks_Hemi1_mva.size() > 1 )
    {
      unsigned int idxmax = 0;
      for (unsigned int p = 0; p<displacedTracks_Hemi1_mva.size(); p++)
      {
        idxmax = p;
        for (unsigned int k = p+1; k <displacedTracks_Hemi1_mva.size(); k++)
        {
          if ( PairSortedHemi1.second[k] > PairSortedHemi1.second[idxmax] )
          {
            std::pair<reco::TransientTrack,float> tempSorted;
            tempSorted.first  = PairSortedHemi1.first[idxmax];
            tempSorted.second = PairSortedHemi1.second[idxmax];
            PairSortedHemi1.first[idxmax]  = PairSortedHemi1.first[k];
            PairSortedHemi1.second[idxmax] = PairSortedHemi1.second[k];
            PairSortedHemi1.first[k]  = tempSorted.first;
            PairSortedHemi1.second[k] = tempSorted.second;
          }
        }
      }
    }
    //--Hemi2--//
    if ( displacedTracks_Hemi2_mva.size() > 1 )
    {
      unsigned int idxmax2 = 0;
      for (unsigned int p = 0; p<displacedTracks_Hemi2_mva.size(); p++)
      {
        idxmax2 = p;
        for (unsigned int k = p+1; k <displacedTracks_Hemi2_mva.size(); k++)
        {
          if ( PairSortedHemi2.second[k] > PairSortedHemi2.second[idxmax2] )
          {
             std::pair<reco::TransientTrack,float> tempSorted;
             tempSorted.first = PairSortedHemi2.first[idxmax2];
             tempSorted.second = PairSortedHemi2.second[idxmax2];
             PairSortedHemi2.first[idxmax2]=PairSortedHemi2.first[k];
             PairSortedHemi2.second[idxmax2]=PairSortedHemi2.second[k];
             PairSortedHemi2.first[k]=tempSorted.first;
             PairSortedHemi2.second[k]=tempSorted.second;
          }
        }
      }
    }

    ///////////////////////////////////////////////////////
    //-----------------------------------------------------
    // Vertex fitting 
    //-----------------------------------------------------
    ///////////////////////////////////////////////////////

    int   Vtx_ntk = 0;
    float Vtx_x = 0., Vtx_y = 0., Vtx_z= 0., Vtx_chi = -10.;
    float recX, recY, recZ, dSV, recD;
    float Vtx_TotalChi=0.;
    float Vtx_DOF=0.;
    float MeanBdtValue = 0.1;
    float UpperLimit = 10.;
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
    int ntracks=-2;
    float tempchi2 = -10.;
    float DOF = -10.;
    float temptChi2 = -10.;
    float tempx = -100.;
    float tempy = -100.;
    float tempz = -100.;
    float mean= -2.;
    bool Vtx_Valid = false;
    if ( Vtx_ntk > 1 )
    {
      // TransientVertex displacedVertex_Hemi1_mva = theFitter_Vertex_Hemi1_mva.vertex(displacedTracks_Hemi1_mva); // fitted vertex
      TransientVertex displacedVertex_Hemi1_mva = theFitter_Vertex_Hemi1_mva.vertex(PairSortedHemi1.first); // fitted vertex

      if ( displacedVertex_Hemi1_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      { 
        Vtx_x = displacedVertex_Hemi1_mva.position().x();
        Vtx_y = displacedVertex_Hemi1_mva.position().y();
        Vtx_z = displacedVertex_Hemi1_mva.position().z();
        Vtx_chi = displacedVertex_Hemi1_mva.normalisedChiSquared();
        Vtx_ntk = displacedTracks_Hemi1_mva.size();
        for (int p = 0; p < Vtx_ntk; p++)
        {
          // tree_LLP_Vtx_trackWeight.push_back(displacedVertex_Hemi1_mva.trackWeight(displacedTracks_Hemi1_mva[p]));
          tree_LLP_Vtx_trackWeight.push_back(displacedVertex_Hemi1_mva.trackWeight(PairSortedHemi1.first[p]));
        }
        mean=0;
        for (int y = 0; y < Vtx_ntk ; y++)
          {
            mean = mean + PairSortedHemi1.second[y];
          }
          mean = mean/Vtx_ntk;
                            //-----------------------IAVF-Hemi1----------------------------//
        std::vector<TransientTrack> vTT;

        if ( Vtx_chi < 0 )
        {
          for ( int p = 0; p < Vtx_ntk; p++ ) // loop on all the TT
          {
            for  ( int k = p+1; k < Vtx_ntk; k++ ) // loop on all the other TT
            {            
              vTT.push_back(PairSortedHemi1.first[p]);
              vTT.push_back(PairSortedHemi1.first[k]);
              TransientVertex TV = theFitter_Vertex_Hemi1_mva.vertex(vTT); // We take the first "good-looking" seed to start but does not discard any track from the collection
              if (TV.isValid() && TV.normalisedChiSquared() > 0  && TV.normalisedChiSquared() < UpperLimit) // && TV.normalisedChiSquared()>0 && abs(TV.normalisedChiSquared())<10
              { 
                mean = (PairSortedHemi1.second[k]+PairSortedHemi1.second[p]);
                for (int m = 0; m < Vtx_ntk; m++)// We then add track by track to the vertex and check the validity of the vertex
                {
                if (m == k || m == p) continue;
                  vTT.push_back(PairSortedHemi1.first[m]);
                  TransientVertex updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                  tempchi2 = updatedTV.normalisedChiSquared();
                  DOF = updatedTV.degreesOfFreedom();
                  temptChi2 = updatedTV.totalChiSquared();
                  if ( !updatedTV.isValid() || updatedTV.normalisedChiSquared() < 0 || updatedTV.normalisedChiSquared() > UpperLimit ) 
		                {  
		   		   		   		    vTT.pop_back();
		    		   		   		  updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
		    		   		   		  tempchi2 = updatedTV.normalisedChiSquared();
		    		   		   		  tempx=updatedTV.position().x();
		    		   		   		  tempy=updatedTV.position().y();
                    		  tempz=updatedTV.position().z(); 
		    		   		   		  DOF = updatedTV.degreesOfFreedom();
		    		   		   		  temptChi2 = updatedTV.totalChiSquared();
		    		   		   		  continue;
		  		   		   } 
                  else
		               {
        		   		    mean = mean + PairSortedHemi1.second[m];        
		                  tempx=updatedTV.position().x();
		                  tempy=updatedTV.position().y();
		                  tempz=updatedTV.position().z();
		               }
                }
                Vtx_chi = tempchi2;
                Vtx_TotalChi = temptChi2;
                Vtx_DOF = DOF;
                Vtx_x = tempx;
                Vtx_y = tempy;
                Vtx_z = tempz;
                ntracks=vTT.size();
                mean = mean/vTT.size();
                Vtx_Valid = true;//Final vertex is valid
                break;  //can potentially give a bad seed	=> usually happens for a low a amount of tracks (few % of the event)   
              }
              else
              { // If not valid : we build a new seed
                ntracks=0;
                vTT.clear();
                Vtx_Valid = false;//FInal Vertex is not valid with the considered seed, so try a new seed
                continue;// can have no output
              }
              // break; //never encountered
            }
            // break;//This has to be changed, can give a bad seed (espcially when there is low number of track)
            if (Vtx_Valid==false || mean <MeanBdtValue) continue ;//MeanBdtValue is arbitrary : ttree->Draw("tree_Hemi_Vtx_dd","abs(tree_Hemi_Vtx_NChi2)<10 && tree_Hemi_Vtx_MeanBdtVal<MeanBdtValue")
            else if (Vtx_Valid==true && mean < MeanBdtValue) continue;
            else if (Vtx_Valid==true && mean > MeanBdtValue) break;// (Vtx_Valid==true && mean > MeanBdtValue)
          }
        }
      }
      else
      {
        std::vector<TransientTrack> vTT;


        for (int p = 0; p < Vtx_ntk; p++) // loop on all the TT
        {
          for  (int k = p+1; k < Vtx_ntk; k++ ) //loop on all the other TT
          {
            vTT.push_back(PairSortedHemi1.first[p]);
            vTT.push_back(PairSortedHemi1.first[k]);
            TransientVertex TV = theFitter_Vertex_Hemi1_mva.vertex(vTT); // We take the first "good-looking" seed to start
            if ( TV.isValid() && TV.normalisedChiSquared() > 0 && TV.normalisedChiSquared() < UpperLimit) // && TV.normalisedChiSquared()>0 && abs(TV.normalisedChiSquared())<10
            { 
              mean = (PairSortedHemi1.second[k]+PairSortedHemi1.second[p]);
              for (int m = 0; m < Vtx_ntk ; m++)// We then add track by track to the vertex and check the validity of the vertex
              {
              if (m == k || m == p) continue;
                vTT.push_back(PairSortedHemi1.first[m]);
                TransientVertex updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                tempchi2 = updatedTV.normalisedChiSquared();
                DOF = updatedTV.degreesOfFreedom();
                temptChi2 = updatedTV.totalChiSquared();
                if (!updatedTV.isValid() || updatedTV.normalisedChiSquared()<0 || updatedTV.normalisedChiSquared()>UpperLimit) 
		              {
		                vTT.pop_back();
		  		  		  	updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
		  		  		    tempchi2=updatedTV.normalisedChiSquared();
		  		  		    tempx=updatedTV.position().x();
		  		  		    tempy=updatedTV.position().y();
                    tempz=updatedTV.position().z();
		  		  		    DOF = updatedTV.degreesOfFreedom();
		  		  		    temptChi2 = updatedTV.totalChiSquared();
		  		  		    continue;
		              }
                else
		              {
                    mean = mean + PairSortedHemi1.second[m];
		                tempx=updatedTV.position().x();
		                tempy=updatedTV.position().y();
		                tempz=updatedTV.position().z();
		              }  
              }
              Vtx_chi = tempchi2;
              Vtx_TotalChi = temptChi2;
              Vtx_DOF=DOF;
              Vtx_x = tempx;
              Vtx_y = tempy;
              Vtx_z = tempz;
              ntracks=vTT.size();
              mean = mean/vTT.size();
              // if (showlog) std::cout << "OriginalTracks : "<<Vtx_ntk<<" vs Actual tracks : "<<ntracks<<"and Final chi2 of : "<<Vtx_chi<<std::endl;
              Vtx_Valid = true;
              break;		 
            }
            else
            {// If not valid : we build a new seed
              ntracks=0;
              Vtx_Valid = false;
              vTT.clear();
              continue;
            }
            // break;//never encountered
          }
          // break;//this need to be changed
            if (Vtx_Valid==false || mean <MeanBdtValue) continue ;//MeanBdtValue is arbitrary : ttree->Draw("tree_Hemi_Vtx_dd","abs(tree_Hemi_Vtx_NChi2)<10 && tree_Hemi_Vtx_MeanBdtVal<MeanBdtValue")
            else if (Vtx_Valid==true && mean < MeanBdtValue) continue;
            else if (Vtx_Valid==true && mean > MeanBdtValue) break;// (Vtx_Valid==true && mean > MeanBdtValue)
        }
      }
                    //--------------------ENDOF IAVF--------------------------//
    }
    std::cout<<"the mean of the tracks used : "<<mean<<std::endl;
    tree_Hemi_Vtx_MeanBdtVal.push_back(mean);
    // !! New information caused by iterative AVF !! //
    tree_Hemi_Vtx_UpdatedChi2.push_back(tempchi2);
    tree_Hemi_Vtx_ntrk.push_back(ntracks);
    tree_Hemi_Vtx_TotalChi2.push_back(Vtx_TotalChi);
    tree_Hemi_Vtx_DOF.push_back(Vtx_DOF);
    //----------------------------------------//   

    float Vtx_chi1 = Vtx_chi;
    tree_Hemi.push_back(1);
    tree_Hemi_njet.push_back(njet1);
    tree_Hemi_eta.push_back(axis1_eta);
    tree_Hemi_phi.push_back(axis1_phi);
    tree_Hemi_dR.push_back(axis1_dR);
    tree_Hemi_nTrks.push_back(nTrks_axis1);
    tree_Hemi_nTrks_sig.push_back(nTrks_axis1_sig);
    tree_Hemi_nTrks_bad.push_back(nTrks_axis1_bad);
    tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
    tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
    tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis1_sig_mva);
    tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis1_bad_mva);
    tree_Hemi_Vtx_x.push_back(Vtx_x);
    tree_Hemi_Vtx_y.push_back(Vtx_y);
    tree_Hemi_Vtx_z.push_back(Vtx_z);
    recX = Vtx_x - tree_PV_x[0];
    recY = Vtx_y - tree_PV_y[0];
    recZ = Vtx_z - tree_PV_z[0];
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_Hemi_Vtx_dist.push_back( recD );
//$$
    float ddok, ddbad;
    float ping1 = 0;
//$$
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
//$$
      dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
      ddok = TMath::Sqrt(dSV)/LLP1_dist;
      tree_Hemi_Vtx_dd.push_back( ddok );
      dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
      ddbad = TMath::Sqrt(dSV)/LLP2_dist;
      tree_Hemi_Vtx_ddbad.push_back( ddbad );
      if      ( abs(Vtx_chi) < 10. && ddok  < 0.1 ) ping1 = 1;
      else if ( abs(Vtx_chi) < 10. && ddbad < 0.1 ) ping1 = 2;
//$$
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
//$$
      dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
      ddok = TMath::Sqrt(dSV)/LLP2_dist;
      tree_Hemi_Vtx_dd.push_back( ddok );
      dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
      ddbad = TMath::Sqrt(dSV)/LLP1_dist;
      tree_Hemi_Vtx_ddbad.push_back( ddbad );
      if      ( abs(Vtx_chi) < 10. && ddok  < 0.1 ) ping1 = 2;
      else if ( abs(Vtx_chi) < 10. && ddbad < 0.1 ) ping1 = 1;
//$$
    }
    tree_Hemi_LLP.push_back(iLLPrec1);
         

    //--------------------------- SECOND HEMISPHERE WITH MVA -------------------------------------//

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
    ntracks=-2;
    tempchi2= -10.;
    mean = -2.;
    Vtx_Valid = false; 
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_Hemi2_mva = theFitter_Vertex_Hemi2_mva.vertex(displacedTracks_Hemi2_mva); // fitted vertex
      
      if ( displacedVertex_Hemi2_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      {
        Vtx_x = displacedVertex_Hemi2_mva.position().x();
        Vtx_y = displacedVertex_Hemi2_mva.position().y();
        Vtx_z = displacedVertex_Hemi2_mva.position().z();
        Vtx_chi = displacedVertex_Hemi2_mva.normalisedChiSquared();
        mean=0;
        for (int y = 0; y < Vtx_ntk ; y++)
          {
            mean = mean + PairSortedHemi2.second[y];
          }
          mean = mean/Vtx_ntk;
        for (int p =0; p<Vtx_ntk; p++)
        {
          // tree_LLP_Vtx_trackWeight.push_back(displacedVertex_Hemi2_mva.trackWeight(displacedTracks_Hemi2_mva[p]));
          tree_LLP_Vtx_trackWeight.push_back(displacedVertex_Hemi2_mva.trackWeight(PairSortedHemi2.first[p]));
          //  if(showlog) std::cout<<"2nd vtx_chi  / weight / chi2 / ndof / NChi2 : "<<Vtx_chi<<" / " <<displacedVertex_Hemi2_mva.trackWeight(displacedTracks_Hemi2_mva[p])<<" / "<<displacedTracks_Hemi2_mva[p].chi2()<<" / "<<displacedTracks_Hemi2_mva[p].ndof()<<" / "<<displacedTracks_Hemi2_mva[p].normalizedChi2()<<std::endl;
        }
                    //-----------------------IAVF-Hemi2---------------------------//
        std::vector<TransientTrack> vTT;
        if (Vtx_chi<0)
        {
          for  (int p = 0; p<Vtx_ntk; p++) //loop on all the TT
          {
            for  (int k = p+1; k < Vtx_ntk; k++) //loop on all the other TT
            {
              // vTT.push_back(displacedTracks_Hemi2_mva[p]);
              // vTT.push_back(displacedTracks_Hemi2_mva[k]);
              vTT.push_back(PairSortedHemi2.first[p]);
              vTT.push_back(PairSortedHemi2.first[k]);
              TransientVertex TV = theFitter_Vertex_Hemi2_mva.vertex(vTT);//We take the first "good-looking" seed to start,
              if (TV.isValid() && TV.normalisedChiSquared()>0 && TV.normalisedChiSquared()<UpperLimit) // && TV.normalisedChiSquared()>0 && abs(TV.normalisedChiSquared())<10
              { 
                mean = (PairSortedHemi2.second[k]+PairSortedHemi2.second[p]);
                if (showlog) std::cout << "2nd LLP seed is created for k and p : "<<k<<" / "<<p<<" with chi2: "<<TV.normalisedChiSquared()<<std::endl;
                for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                {
                if (m == k || m == p) continue;
                  // vTT.push_back(displacedTracks_Hemi2_mva[m]);
                  vTT.push_back(PairSortedHemi2.first[m]);
                  TransientVertex updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                  if (showlog) std::cout << "m = "<<m<<"/ chi2 of the vetex constructed with this seed : "<<updatedTV.normalisedChiSquared()<<" originally : "<<Vtx_chi<<std::endl;
                  tempchi2=updatedTV.normalisedChiSquared();
                  DOF = updatedTV.degreesOfFreedom();
                  temptChi2 = updatedTV.totalChiSquared();
                  // || updatedTV.normalisedChiSquared()<0
                  // || abs(updatedTV.normalisedChiSquared())>10
                  if ( !updatedTV.isValid() || updatedTV.normalisedChiSquared()<0 || updatedTV.normalisedChiSquared()>UpperLimit) 
		  {
		    vTT.pop_back();
		    updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
		    tempchi2=updatedTV.normalisedChiSquared();
		    tempx=updatedTV.position().x();
		    tempy=updatedTV.position().y();
                    tempz=updatedTV.position().z();
		    DOF = updatedTV.degreesOfFreedom();
		    temptChi2 = updatedTV.totalChiSquared();
		    continue;//pas utile
		  } 
                  else
		  {
        mean = mean + PairSortedHemi2.second[m];
		    tempx=updatedTV.position().x();
		    tempy=updatedTV.position().y();
		    tempz=updatedTV.position().z();
		  }
                }       
                Vtx_chi = tempchi2;
                Vtx_TotalChi = temptChi2;
                Vtx_DOF = DOF;
                Vtx_x = tempx;
                Vtx_y = tempy;
                Vtx_z = tempz;
                Vtx_Valid = true;
                ntracks=vTT.size();
                mean = mean/vTT.size();
                if (showlog) std::cout << "OriginalTracks : "<<Vtx_ntk<<" vs Actual tracks : "<<ntracks<<"and Final chi2 of : "<<Vtx_chi<<std::endl;      
                break;	    
              }
              else
              { // If not valid : we build a new seed
                ntracks = 0;
                Vtx_Valid = false;
                vTT.clear();
                continue;
              }
              // break;//nerver encountered
            }
            // break;   
            if (Vtx_Valid==false || mean <MeanBdtValue) continue ;//MeanBdtValue is arbitrary : ttree->Draw("tree_Hemi_Vtx_dd","abs(tree_Hemi_Vtx_NChi2)<10 && tree_Hemi_Vtx_MeanBdtVal<MeanBdtValue")
            else if (Vtx_Valid==true && mean < MeanBdtValue) continue;
            else if (Vtx_Valid==true && mean > MeanBdtValue) break;// (Vtx_Valid==true && mean > MeanBdtValue)
          }
        }
      }
      else
      {
        std::vector<TransientTrack> vTT;
        for  (int p = 0; p < Vtx_ntk; p++) //loop on all the TT
        {
          for  (int k = p+1; k < Vtx_ntk; k++ ) //loop on all the other TT
          {
            vTT.push_back(PairSortedHemi2.first[p]);
            vTT.push_back(PairSortedHemi2.first[k]);
            TransientVertex TV = theFitter_Vertex_Hemi2_mva.vertex(vTT);//We take the first "good-looking" seed to start, 
            if ( TV.isValid() && TV.normalisedChiSquared() > 0 && TV.normalisedChiSquared() < UpperLimit )// && TV.normalisedChiSquared()>0 && abs(TV.normalisedChiSquared())<10
            { 
              mean = PairSortedHemi2.second[k]+PairSortedHemi2.second[p];
              ntracks=2;
              if (showlog) std::cout << " 2nd LLP seed is created for k and p : "<<k<<" / "<<p<<" with chi2: "<<TV.normalisedChiSquared()<<std::endl;
              for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
              {
              if ( m == k || m == p ) continue;
                // vTT.push_back(displacedTracks_Hemi2_mva[m]);
                vTT.push_back(PairSortedHemi2.first[m]);
                TransientVertex updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                if (showlog) std::cout << " m = " << m << "/ chi2 of the vetex constructed with this seed : "<<updatedTV.normalisedChiSquared()<<" originally : "<<Vtx_chi<<std::endl;
                tempchi2 = updatedTV.normalisedChiSquared();
                DOF = updatedTV.degreesOfFreedom();
                temptChi2 = updatedTV.totalChiSquared();
                if ( !updatedTV.isValid() || updatedTV.normalisedChiSquared()<0  || updatedTV.normalisedChiSquared()> UpperLimit) 
		{
		  vTT.pop_back();
		  updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
		  tempchi2=updatedTV.normalisedChiSquared();
		  tempx=updatedTV.position().x();
		  tempy=updatedTV.position().y();
                  tempz=updatedTV.position().z();
		  DOF = updatedTV.degreesOfFreedom();
		  temptChi2 = updatedTV.totalChiSquared();
		  continue;
		}
                else
		{
      mean = mean + PairSortedHemi2.second[m];
		  tempx=updatedTV.position().x();
		  tempy=updatedTV.position().y();
		  tempz=updatedTV.position().z();
		}
              }
              Vtx_chi = tempchi2;
              Vtx_TotalChi = temptChi2;
              Vtx_DOF = DOF;
              Vtx_x = tempx;
              Vtx_y = tempy;
              Vtx_z = tempz;
              Vtx_Valid = true;
              ntracks=vTT.size();
              mean = mean/vTT.size();
              break;		 
            }
            else
            { // If not valid : we build a new seed
              ntracks = 0;
              vTT.clear();
              Vtx_Valid = false;
              continue;
            }
            // break;//never encountered
          }
          // break;//this needs to be changed
            if (Vtx_Valid==false || mean <MeanBdtValue) continue ;//MeanBdtValue is arbitrary : ttree->Draw("tree_Hemi_Vtx_dd","abs(tree_Hemi_Vtx_NChi2)<10 && tree_Hemi_Vtx_MeanBdtVal<MeanBdtValue")
            else if (Vtx_Valid==true && mean < MeanBdtValue) continue;
            else if (Vtx_Valid==true && mean > MeanBdtValue) break;// (Vtx_Valid==true && mean > MeanBdtValue)
        }
      //-----------------------ENFO OF IAVF HEMI2---------------------//
      }
    } 
              std::cout<<"the mean of the tracks used : "<<mean<<std::endl;
    tree_Hemi_Vtx_MeanBdtVal.push_back(mean);  
    // !!New information caused by iterative AVF !! //
    tree_Hemi_Vtx_UpdatedChi2.push_back(tempchi2);
    tree_Hemi_Vtx_ntrk.push_back(ntracks); 
    tree_Hemi_Vtx_TotalChi2.push_back(Vtx_TotalChi);
    tree_Hemi_Vtx_DOF.push_back(Vtx_DOF);
    //----------------------------------------//

    float Vtx_chi2 = Vtx_chi;
    tree_Hemi.push_back(2);
    tree_Hemi_njet.push_back(njet2);
    tree_Hemi_eta.push_back(axis2_eta);
    tree_Hemi_phi.push_back(axis2_phi);
    tree_Hemi_dR.push_back(axis2_dR);
    tree_Hemi_nTrks.push_back(nTrks_axis2);
    tree_Hemi_nTrks_sig.push_back(nTrks_axis2_sig);
    tree_Hemi_nTrks_bad.push_back(nTrks_axis2_bad);
    tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
    tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
    tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis2_sig_mva);
    tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis2_bad_mva);
    tree_Hemi_Vtx_x.push_back(Vtx_x);
    tree_Hemi_Vtx_y.push_back(Vtx_y);
    tree_Hemi_Vtx_z.push_back(Vtx_z);
    recX = Vtx_x - tree_PV_x[0];
    recY = Vtx_y - tree_PV_y[0];
    recZ = Vtx_z - tree_PV_z[0];
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_Hemi_Vtx_dist.push_back( recD );
//$$
    float ping2 = 0;
//$$
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
//$$
      dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
      ddok = TMath::Sqrt(dSV)/LLP1_dist;
      tree_Hemi_Vtx_dd.push_back( ddok );
      dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
      ddbad = TMath::Sqrt(dSV)/LLP2_dist;
      tree_Hemi_Vtx_ddbad.push_back( ddbad );
      if      ( abs(Vtx_chi) < 10. && ddok  < 0.1 ) ping2 = 1;
      else if ( abs(Vtx_chi) < 10. && ddbad < 0.1 ) ping2 = 2;
//$$
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
//$$
      dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
      ddok = TMath::Sqrt(dSV)/LLP2_dist;
      tree_Hemi_Vtx_dd.push_back( ddok );
      dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
      ddbad = TMath::Sqrt(dSV)/LLP1_dist;
      tree_Hemi_Vtx_ddbad.push_back( ddbad );
      if      ( abs(Vtx_chi) < 10. && ddok  < 0.1 ) ping2 = 2;
      else if ( abs(Vtx_chi) < 10. && ddbad < 0.1 ) ping2 = 1;
//$$
    }
    tree_Hemi_LLP.push_back(iLLPrec2);
    
    tree_Hemi_dR12.push_back(dR_axis12);
    tree_Hemi_dR12.push_back(dR_axis12);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);

//$$
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
//$$
      
    // some informations for tracks in their hemisphere
//$$
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
//$$

     
  //////////////////////////////////
  // }//end passes htfilter
  smalltree->Fill();
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
  float countChi2 = 0.;
  float countEff = 0.;
  float effChi2;
  float Vtxeff;
  std::cout<<"size : "<<tree_Hemi.size()<<std::endl;
for(unsigned int i=0;i<20000;i++)
  {
    if (abs(tree_Hemi_Vtx_NChi2[i])<10)
      {
        countChi2=countChi2+1;
        if (tree_Hemi_Vtx_dd[i]<0.1)
          {
            countEff=countEff+1;
          }
      }
  }
  effChi2= countChi2/20000.;
  Vtxeff = countEff/20000.;

  std::cout<<"counts : "<<countChi2<<" / "<<countEff<<std::endl;
  std::cout<<"//----------Efficiencies-----------\\"<<std::endl;
  std::cout<<"||                                 ||"<<std::endl;
  std::cout<<"|| Vertex w/ good Chi2 : "<<effChi2<<"% ||"<<std::endl;
  std::cout<<"|| Vertex RecoEff : "<<Vtxeff<<"%       ||"<<std::endl;
  std::cout<<"||                                  ||"<<std::endl;
  std::cout<<"\\----------------------------------//"<<std::endl;
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
    
//     tree_trigger_names.clear();
//     tree_trigger_bits.clear();
    //trig
    tree_trigger_size.clear();
    tree_passesTrigger.clear();
    tree_passesTriggerName.clear();
    //trig
    tree_PV_x.clear(); // l'index 0 donne le PV!
    tree_PV_y.clear();
    tree_PV_z.clear();
    tree_PV_ez.clear();
    tree_PV_NChi2.clear();
    tree_PV_ndf.clear();

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
    tree_track_MVAval.clear();
    
    tree_track_Hemi.clear();
    tree_track_Hemi_dR.clear();
    tree_track_Hemi_mva_NChi2.clear();
//$$
    tree_track_Hemi_ping.clear();
    tree_track_Hemi_dFirstVtx.clear();
//$$
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
    tree_genPackPart_mother_pdgId.clear();

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
    tree_LLP_Vtx_trackWeight.clear();    

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
    tree_Hemi_Vtx_NChi2.clear();
    tree_Hemi_Vtx_nTrks.clear();
    tree_Hemi_Vtx_nTrks_sig.clear();
    tree_Hemi_Vtx_nTrks_bad.clear();
    tree_Hemi_Vtx_x.clear();
    tree_Hemi_Vtx_y.clear();
    tree_Hemi_Vtx_z.clear();
    tree_Hemi_Vtx_dist.clear();
    tree_Hemi_Vtx_dx.clear();
    tree_Hemi_Vtx_dy.clear();
    tree_Hemi_Vtx_dz.clear();
    tree_Hemi_Vtx_dd.clear();
    tree_Hemi_dR12.clear();
    tree_Hemi_LLP_dR12.clear();
    tree_Hemi_Vtx_UpdatedChi2.clear();
    tree_Hemi_Vtx_ntrk.clear();
    tree_Hemi_Vtx_TotalChi2.clear();
    tree_Hemi_Vtx_DOF.clear();
    tree_Hemi_Vtx_MeanBdtVal.clear();
//$$
    tree_Hemi_Vtx_ddbad.clear();
    tree_Hemi_Vtx_ntrk10.clear();
    tree_Hemi_Vtx_ntrk20.clear();
    tree_Hemi_Vtx_trackWeight.clear();
    tree_Hemi_LLP_ping.clear();
    tree_event_LLP_ping.clear();
//$$
}
