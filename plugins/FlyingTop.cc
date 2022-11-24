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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

//!!!!!!
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/VecArray.h"
#include "FWCore/Utilities/interface/isFinite.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/Math/interface/liblogintpack.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
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
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
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
    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<pat::JetCollection> jetToken_;
    edm::EDGetTokenT<edm::View<reco::GenJet> > genJetToken_;
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
//$$    edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
//$$
    edm::EDGetTokenT<pat::PackedCandidateCollection> pc_;
    edm::EDGetTokenT<edm::Association<reco::PFCandidateCollection>> pc2pf_;
//$$
    
    std::string parametersDefinerName_;

  ///////////////
  // Ntuple info

    TTree *smalltree;
    
    edm::Service<TFileService> fs;
    
//     std::string ttrhbuilder_;
    
    edm::ESHandle<MagneticField> bField;

    edm::ParameterSet kvfPSet;
       
    int runNumber, eventNumber, lumiBlock;
    int  tree_NbrOfZCand;
    bool tree_passesHTFilter;
    int  tree_nTracks; 
    int  nBC = 0, nFromC = 0, nFromB = 0; 
    int nEvent;
    
    //-----------------------
    // trigger variable
//     std::vector<string > tree_trigger_names;
//     std::vector<bool >   tree_trigger_bits;
    
    //--------------------------------
    // primary vertex infos -------
    //--------------------------------
  
    int   tree_nPV;
    float tree_PV_x;
    float tree_PV_y;
    float tree_PV_z;
    float tree_PV_ez;
    float tree_PV_NChi2;
    float tree_PV_ndf;
    
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
    
    int tree_nAK4jet;
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
    std::vector< int >   tree_track_iJet;
    std::vector< float > tree_track_pt;
    std::vector< float > tree_track_outerPt;
    std::vector< float > tree_track_eta;
    std::vector< float > tree_track_phi;
    std::vector<int>     tree_track_charge;
    std::vector<float >  tree_track_NChi2;
    std::vector<bool >   tree_track_isHighQuality;
    std::vector<bool >   tree_track_isLoose;
    std::vector<bool >   tree_track_isTight;
    std::vector< float>  tree_track_dxy; // dxy with respect to beam spot position
    std::vector< float>  tree_track_dxyError;
    std::vector< float>  tree_track_dz;
    std::vector< float>  tree_track_dzError ;
    std::vector<int>     tree_track_numberOfLostHits;
    std::vector<unsigned int>    tree_track_originalAlgo; // definition as comments at the end of the file,
    //http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_1_3/doc/html/d8/df2/classreco_1_1TrackBase.html#aca7611bd1a33d535cefc72b6e497ece8
    std::vector<unsigned int>    tree_track_algo;
    std::vector<unsigned short>  tree_track_stopReason;
    
    std::vector<int>     tree_track_nHit;
    std::vector<int>     tree_track_nHitPixel;
    std::vector<int>     tree_track_nHitTIB;
    std::vector<int>     tree_track_nHitTID;
    std::vector<int>     tree_track_nHitTOB;
    std::vector<int>     tree_track_nHitTEC;
    std::vector<int>     tree_track_nHitPXB;
    std::vector<int>     tree_track_nHitPXF;
    std::vector<int>     tree_track_isHitL1;
    std::vector<int>     tree_track_nLayers;
    std::vector<int>     tree_track_nLayersPixel;
    std::vector<int>     tree_track_stripTECLayersWithMeasurement ;
    std::vector<int>     tree_track_stripTIBLayersWithMeasurement;
    std::vector<int>     tree_track_stripTIDLayersWithMeasurement;
    std::vector<int>     tree_track_stripTOBLayersWithMeasurement;
    std::vector<int>     tree_track_hitpattern;
    std::vector< float >  tree_track_vx;
    std::vector< float >  tree_track_vy;
    std::vector< float >  tree_track_vz;
    std::vector<float>    tree_track_firsthit_X;
    std::vector<float>    tree_track_firsthit_Y;
    std::vector<float>    tree_track_firsthit_Z;
    std::vector<float>    tree_track_region;
    std::vector<float>    tree_track_ntrk10;
    std::vector<float>    tree_track_ntrk20;
    std::vector<float>    tree_track_ntrk30;

    //!!!!
 
    std::vector<float>   tree_track_GeoBarrel_firsthit_x;
    std::vector<float>   tree_track_GeoBarrel_firsthit_y;
    std::vector<float>   tree_track_GeoBarrel_firsthit_z;

    std::vector<float>                         tree_track_PropBarrel_firsthit_x;
    std::vector<float>                         tree_track_PropBarrel_firsthit_y;
    std::vector<float>                         tree_track_PropBarrel_firsthit_z;

    std::vector<float>                         tree_track_GeoDisk_firsthit_x;
    std::vector<float>                         tree_track_GeoDisk_firsthit_y;
    std::vector<float>                         tree_track_GeoDisk_firsthit_z;

    std::vector<float>                         tree_track_PropDisk_SLPC_firsthit_x;
    std::vector<float>                         tree_track_PropDisk_SLPC_firsthit_y;
    std::vector<float>                         tree_track_PropDisk_SLPC_firsthit_z;
    //!!!!

    std::vector<double>    tree_track_MVAval;

    std::vector<int>       tree_track_Hemi;
    std::vector<double>    tree_track_Hemi_dR;
    std::vector<double>    tree_track_Hemi_mva_NChi2;
    std::vector<int>       tree_track_Hemi_LLP;
    
    std::vector<int>    tree_track_recoVertex_idx;
    std::vector<int>    tree_track_recoAK4SlimmedJet_idx;
    std::vector<int>    tree_track_recoAK4PFJet_idx;
    std::vector<int>    tree_track_reco08Jet_idx;
    std::vector<int>    tree_track_recoCaloJet_idx;
    
    std::vector<int>      tree_track_nSimHits;
    std::vector<bool>     tree_track_isSimMatched;
    
    std::vector< int >    tree_track_sim_charge;
    std::vector< float >  tree_track_sim_pt;
    std::vector< float >  tree_track_sim_eta  ;
    std::vector< float >  tree_track_sim_phi  ;
    std::vector<bool>     tree_track_sim_longLived      ;
    // std::vector<int>   tree_track_sim_matchedHit    ;
    std::vector<int>      tree_track_sim_pdgId;
    std::vector<int>      tree_track_sim_numberOfTrackerHits  ;
    std::vector<int>      tree_track_sim_numberOfTrackerLayers;
    std::vector<float>    tree_track_sim_mass  ;
    std::vector<int>      tree_track_sim_status;
    
    std::vector<float>    tree_track_sim_vx;
    std::vector<float>    tree_track_sim_vy;
    std::vector<float>    tree_track_sim_vz;
    std::vector<int>      tree_track_sim_isFromLLP;

    std::vector<int>      tree_track_sim_isFromBC;
    std::vector<int>      tree_track_sim_isFromBC_mother_pdgId;
    std::vector<int>      tree_track_sim_isFromBC_LLP;
    std::vector<float>    tree_track_sim_dV;
    std::vector<float>    tree_track_sim_dist;
    
    //--------------------------------
    // gen infos -------
    //--------------------------------
    float tree_GenPVx;
    float tree_GenPVy;
    float tree_GenPVz;
    
    std::vector< float > tree_genAxis_dRneuneu;
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

    std::vector< float > tree_genPackPart_pt;
    std::vector< float > tree_genPackPart_eta;
    std::vector< float > tree_genPackPart_phi;
    std::vector< float > tree_genPackPart_charge;
    std::vector< int >   tree_genPackPart_pdgId;
    std::vector< float > tree_genPackPart_x;
    std::vector< float > tree_genPackPart_y;
    std::vector< float > tree_genPackPart_z;
    std::vector< float > tree_genPackPart_mass;
    std::vector< int >   tree_genPackPart_statusCode;
    std::vector< int >   tree_genPackPart_mother_pdgId;
    std::vector< int >   tree_genPackPart_nDaughters;

    std::vector< int >   tree_genFromLLP_LLP;
    std::vector< float > tree_genFromLLP_LLPx;
    std::vector< float > tree_genFromLLP_LLPy;
    std::vector< float > tree_genFromLLP_LLPz;
    std::vector< float > tree_genFromLLP_pt;
    std::vector< float > tree_genFromLLP_eta;
    std::vector< float > tree_genFromLLP_phi;
    std::vector< float > tree_genFromLLP_charge;
    std::vector< int >   tree_genFromLLP_pdgId;
    std::vector< float > tree_genFromLLP_x;
    std::vector< float > tree_genFromLLP_y;
    std::vector< float > tree_genFromLLP_z;
    std::vector< float > tree_genFromLLP_mass;
    std::vector< int >   tree_genFromLLP_statusCode;
    std::vector< int >   tree_genFromLLP_mother_pdgId;

    std::vector< float > tree_genFromC_pt;
    std::vector< float > tree_genFromC_eta;
    std::vector< float > tree_genFromC_phi;
    std::vector< float > tree_genFromC_charge;
    std::vector< int >   tree_genFromC_pdgId;
    std::vector< float > tree_genFromC_x;
    std::vector< float > tree_genFromC_y;
    std::vector< float > tree_genFromC_z;
    std::vector< int >   tree_genFromC_mother_pdgId;
    std::vector< int >   tree_genFromC_generation;
    std::vector< int >   tree_genFromC_LLP;

    std::vector< float > tree_genFromB_pt;
    std::vector< float > tree_genFromB_eta;
    std::vector< float > tree_genFromB_phi;
    std::vector< float > tree_genFromB_charge;
    std::vector< int >   tree_genFromB_pdgId;
    std::vector< float > tree_genFromB_x;
    std::vector< float > tree_genFromB_y;
    std::vector< float > tree_genFromB_z;
    std::vector< int >   tree_genFromB_mother_pdgId;
    std::vector< int >   tree_genFromB_generation;
    std::vector< int >   tree_genFromB_LLP;
   
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
    float LLP1_pt, LLP1_eta, LLP1_phi, LLP2_pt, LLP2_eta, LLP2_phi;
    float LLP1_x, LLP1_y, LLP1_z, LLP2_x, LLP2_y, LLP2_z;
    float LLP1_dist, LLP2_dist;
    int   LLP1_nTrks = 0, LLP2_nTrks = 0;
    int   tree_nLLP = -1;

    std::vector< int >   tree_LLP;
    std::vector< float > tree_LLP_pt;
    std::vector< float > tree_LLP_eta;
    std::vector< float > tree_LLP_phi;
    std::vector< float > tree_LLP_x;
    std::vector< float > tree_LLP_y;
    std::vector< float > tree_LLP_z;
    std::vector< int >   tree_LLP_nTrks;
    std::vector< int >   tree_LLP_Vtx_nTrks;
    std::vector< float > tree_LLP_Vtx_NChi2;
    std::vector< float > tree_LLP_Vtx_dx;
    std::vector< float > tree_LLP_Vtx_dy;
    std::vector< float > tree_LLP_Vtx_dz;
    
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
    std::vector< float > tree_Hemi_Vtx_dx;
    std::vector< float > tree_Hemi_Vtx_dy;
    std::vector< float > tree_Hemi_Vtx_dz;
    std::vector< float > tree_Hemi_dR12;
    std::vector< float > tree_Hemi_LLP_dR12;
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
//  
    weightFile_(                                                 iConfig.getUntrackedParameter<std::string>("weightFileMVA") ),
    prunedGenToken_(consumes<edm::View<reco::GenParticle> >(     iConfig.getParameter<edm::InputTag>("genpruned"))),
    packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genpacked"))),
    vertexToken_(   consumes<reco::VertexCollection>(            iConfig.getParameter<edm::InputTag>("vertices"))),
    metToken_(      consumes<pat::METCollection>(                iConfig.getParameter<edm::InputTag>("mets"))),
    jetToken_(      consumes<pat::JetCollection>(                iConfig.getParameter<edm::InputTag>("jets"))),
    genJetToken_(   consumes<edm::View<reco::GenJet>>(           iConfig.getParameter<edm::InputTag>("genjets"))),
    electronToken_( consumes<pat::ElectronCollection>(           iConfig.getParameter<edm::InputTag>("electrons"))),
    muonToken_(     consumes<pat::MuonCollection>(               iConfig.getParameter<edm::InputTag>("muons"))),
//$$    pfToken_(       consumes<pat::PackedCandidateCollection>(    iConfig.getParameter<edm::InputTag>("pfCands")))
//$$
     pc_(           consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
     pc2pf_(        consumes<edm::Association<reco::PFCandidateCollection>>( iConfig.getParameter<edm::InputTag>("pfCands")))
//$$
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
    smalltree->Branch("tree_PV_x", &tree_PV_x);
    smalltree->Branch("tree_PV_y", &tree_PV_y);
    smalltree->Branch("tree_PV_z", &tree_PV_z);    
    smalltree->Branch("tree_PV_ez", &tree_PV_ez);    
    smalltree->Branch("tree_PV_NChi2", &tree_PV_NChi2);    
    smalltree->Branch("tree_PV_ndf", &tree_PV_ndf);    
    
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
    
    smalltree->Branch("tree_NbrOfZCand",  &tree_NbrOfZCand,  "tree_NbrOfZCand/I");
    smalltree->Branch("tree_passesHTFilter", &tree_passesHTFilter);
    
    // met info
    smalltree->Branch("tree_PFMet_et" ,  &tree_PFMet_et);
    smalltree->Branch("tree_PFMet_phi" , &tree_PFMet_phi);
    smalltree->Branch("tree_PFMet_sig" , &tree_PFMet_sig);
    
    // jet info
    smalltree->Branch("tree_nAK4jet"  ,        &tree_nAK4jet);
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
    smalltree->Branch("tree_nTracks", &tree_nTracks, "tree_nTracks/I");
    smalltree->Branch("tree_track_iJet",&tree_track_iJet);
    smalltree->Branch("tree_track_pt",            &tree_track_pt);
    smalltree->Branch("tree_track_outerPt",           &tree_track_outerPt );
    smalltree->Branch("tree_track_eta",                  &tree_track_eta );
    smalltree->Branch("tree_track_phi",                  &tree_track_phi );
    smalltree->Branch("tree_track_charge",            &tree_track_charge );
    smalltree->Branch("tree_track_NChi2",                &tree_track_NChi2);
    smalltree->Branch("tree_track_isHighQuality",     &tree_track_isHighQuality);
    smalltree->Branch("tree_track_isLoose",        &tree_track_isLoose);
    smalltree->Branch("tree_track_isTight",        &tree_track_isTight);
    smalltree->Branch("tree_track_dxy",                  &tree_track_dxy );
    smalltree->Branch("tree_track_dxyError",        &tree_track_dxyError);
    smalltree->Branch("tree_track_dz",            &tree_track_dz);
    smalltree->Branch("tree_track_dzError",        &tree_track_dzError  );
    smalltree->Branch("tree_track_numberOfLostHits",  &tree_track_numberOfLostHits);
    smalltree->Branch("tree_track_originalAlgo",      &tree_track_originalAlgo);
    smalltree->Branch("tree_track_algo",             &tree_track_algo);
    smalltree->Branch("tree_track_stopReason",        &tree_track_stopReason);
    smalltree->Branch("tree_track_isSimMatched",      &tree_track_isSimMatched    );
    
    smalltree->Branch("tree_track_nHit",         &tree_track_nHit);
    smalltree->Branch("tree_track_nHitPixel",    &tree_track_nHitPixel);
    smalltree->Branch("tree_track_nHitTIB",      &tree_track_nHitTIB);
    smalltree->Branch("tree_track_nHitTID",      &tree_track_nHitTID);
    smalltree->Branch("tree_track_nHitTOB",      &tree_track_nHitTOB);
    smalltree->Branch("tree_track_nHitTEC",      &tree_track_nHitTEC);
    smalltree->Branch("tree_track_nHitPXB",      &tree_track_nHitPXB);
    smalltree->Branch("tree_track_nHitPXF",      &tree_track_nHitPXF);
    smalltree->Branch("tree_track_isHitL1",      &tree_track_isHitL1);
    smalltree->Branch("tree_track_nLayers",      &tree_track_nLayers);
    smalltree->Branch("tree_track_nLayersPixel", &tree_track_nLayersPixel);
    smalltree->Branch("tree_track_stripTECLayersWithMeasurement", &tree_track_stripTECLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTIBLayersWithMeasurement", &tree_track_stripTIBLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTIDLayersWithMeasurement", &tree_track_stripTIDLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTOBLayersWithMeasurement", &tree_track_stripTOBLayersWithMeasurement);
    smalltree->Branch("tree_track_hitpattern",&tree_track_hitpattern);
    smalltree->Branch("tree_track_vx",           &tree_track_vx );
    smalltree->Branch("tree_track_vy",           &tree_track_vy );
    smalltree->Branch("tree_track_vz",           &tree_track_vz );
    smalltree->Branch("tree_track_firsthit_X",   &tree_track_firsthit_X);
    smalltree->Branch("tree_track_firsthit_Y",   &tree_track_firsthit_Y);
    smalltree->Branch("tree_track_firsthit_Z",   &tree_track_firsthit_Z);
    smalltree->Branch("tree_track_region",&tree_track_region);
    smalltree->Branch("tree_track_ntrk10",&tree_track_ntrk10);
    smalltree->Branch("tree_track_ntrk20",&tree_track_ntrk20);
    smalltree->Branch("tree_track_ntrk30",&tree_track_ntrk30);

    //!!!!
 
    smalltree->Branch("tree_track_GeoBarrel_firsthit_x",&tree_track_GeoBarrel_firsthit_x);
    smalltree->Branch("tree_track_GeoBarrel_firsthit_y",&tree_track_GeoBarrel_firsthit_y);
    smalltree->Branch("tree_track_GeoBarrel_firsthit_z",&tree_track_GeoBarrel_firsthit_z);

    smalltree->Branch("tree_track_PropBarrel_firsthit_x",&tree_track_PropBarrel_firsthit_x);
    smalltree->Branch("tree_track_PropBarrel_firsthit_y",&tree_track_PropBarrel_firsthit_y);
    smalltree->Branch("tree_track_PropBarrel_firsthit_z",&tree_track_PropBarrel_firsthit_z);

    smalltree->Branch("tree_track_GeoDisk_firsthit_x",&tree_track_GeoDisk_firsthit_x);
    smalltree->Branch("tree_track_GeoDisk_firsthit_y",&tree_track_GeoDisk_firsthit_y);
    smalltree->Branch("tree_track_GeoDisk_firsthit_z",&tree_track_GeoDisk_firsthit_z);
    
    smalltree->Branch("tree_track_PropDisk_SLPC_firsthit_x",&tree_track_PropDisk_SLPC_firsthit_x);
    smalltree->Branch("tree_track_PropDisk_SLPC_firsthit_y",&tree_track_PropDisk_SLPC_firsthit_y);
    smalltree->Branch("tree_track_PropDisk_SLPC_firsthit_z",&tree_track_PropDisk_SLPC_firsthit_z);


    //!!!!
    
    smalltree->Branch("tree_track_MVAval", &tree_track_MVAval);
    smalltree->Branch("tree_track_Hemi", &tree_track_Hemi);
    smalltree->Branch("tree_track_Hemi_dR", &tree_track_Hemi_dR);
    smalltree->Branch("tree_track_Hemi_mva_NChi2", &tree_track_Hemi_mva_NChi2);
    smalltree->Branch("tree_track_Hemi_LLP", &tree_track_Hemi_LLP);
    
    smalltree->Branch("tree_track_recoVertex_idx", &tree_track_recoVertex_idx);
    smalltree->Branch("tree_track_recoAK4SlimmedJet_idx", &tree_track_recoAK4SlimmedJet_idx);
    
    // NOT USED IN MINIAOD info about the simulated track matched to the reco track
    smalltree->Branch("tree_track_sim_charge",                &tree_track_sim_charge );
    smalltree->Branch("tree_track_sim_pt",                    &tree_track_sim_pt );
    smalltree->Branch("tree_track_sim_eta",                   &tree_track_sim_eta  );
    smalltree->Branch("tree_track_sim_phi",                   &tree_track_sim_phi  );
    smalltree->Branch("tree_track_sim_longLived",             &tree_track_sim_longLived );
    smalltree->Branch("tree_track_sim_pdgId",                 &tree_track_sim_pdgId );
    smalltree->Branch("tree_track_sim_numberOfTrackerHits",   &tree_track_sim_numberOfTrackerHits   );
    smalltree->Branch("tree_track_sim_numberOfTrackerLayers", &tree_track_sim_numberOfTrackerLayers );
    smalltree->Branch("tree_track_sim_mass",                  &tree_track_sim_mass   );
    smalltree->Branch("tree_track_sim_status",                &tree_track_sim_status );
    smalltree->Branch("tree_track_sim_vx",                    &tree_track_sim_vx );
    smalltree->Branch("tree_track_sim_vy",                    &tree_track_sim_vy );
    smalltree->Branch("tree_track_sim_vz",                    &tree_track_sim_vz );
    smalltree->Branch("tree_track_sim_isFromLLP",             &tree_track_sim_isFromLLP );
    smalltree->Branch("tree_track_sim_isFromBC",              &tree_track_sim_isFromBC );
    smalltree->Branch("tree_track_sim_isFromBC_mother_pdgId", &tree_track_sim_isFromBC_mother_pdgId );
    smalltree->Branch("tree_track_sim_isFromBC_LLP",          &tree_track_sim_isFromBC_LLP );
    smalltree->Branch("tree_track_sim_dV",                    &tree_track_sim_dV );
    smalltree->Branch("tree_track_sim_dist",                  &tree_track_sim_dist );
    
    // gen info
    smalltree->Branch("tree_GenPVx" ,  &tree_GenPVx);
    smalltree->Branch("tree_GenPVy" ,  &tree_GenPVy);
    smalltree->Branch("tree_GenPVz" ,  &tree_GenPVz);
    
    smalltree->Branch("tree_genAxis_dRneuneu",&tree_genAxis_dRneuneu);
    smalltree->Branch("tree_genParticle_pt"  ,          &tree_genParticle_pt);
    smalltree->Branch("tree_genParticle_eta" ,          &tree_genParticle_eta);
    smalltree->Branch("tree_genParticle_phi" ,          &tree_genParticle_phi);
    smalltree->Branch("tree_genParticle_charge" ,       &tree_genParticle_charge);
    smalltree->Branch("tree_genParticle_pdgId" ,        &tree_genParticle_pdgId);
    smalltree->Branch("tree_genParticle_x"  ,	        &tree_genParticle_x);
    smalltree->Branch("tree_genParticle_y" ,	        &tree_genParticle_y);
    smalltree->Branch("tree_genParticle_z" ,	        &tree_genParticle_z);
    smalltree->Branch("tree_genParticle_mass" ,         &tree_genParticle_mass);
    smalltree->Branch("tree_genParticle_statusCode",    &tree_genParticle_statusCode);
    smalltree->Branch("tree_genParticle_mother_pdgId" , &tree_genParticle_mother_pdgId);

    smalltree->Branch("tree_genPackPart_pt"  ,          &tree_genPackPart_pt);
    smalltree->Branch("tree_genPackPart_eta" ,          &tree_genPackPart_eta);
    smalltree->Branch("tree_genPackPart_phi" ,          &tree_genPackPart_phi);
    smalltree->Branch("tree_genPackPart_charge" ,       &tree_genPackPart_charge);
    smalltree->Branch("tree_genPackPart_pdgId" ,        &tree_genPackPart_pdgId);
    smalltree->Branch("tree_genPackPart_x"  ,	        &tree_genPackPart_x);
    smalltree->Branch("tree_genPackPart_y" ,	        &tree_genPackPart_y);
    smalltree->Branch("tree_genPackPart_z" ,	        &tree_genPackPart_z);
    smalltree->Branch("tree_genPackPart_mass" ,         &tree_genPackPart_mass);
    smalltree->Branch("tree_genPackPart_statusCode",    &tree_genPackPart_statusCode);
    smalltree->Branch("tree_genPackPart_mother_pdgId" , &tree_genPackPart_mother_pdgId);
    smalltree->Branch("tree_genPackPart_nDaughters" ,   &tree_genPackPart_nDaughters);

    smalltree->Branch("tree_genFromLLP_LLP"  ,         &tree_genFromLLP_LLP);
    smalltree->Branch("tree_genFromLLP_LLPx"  ,        &tree_genFromLLP_LLPx);
    smalltree->Branch("tree_genFromLLP_LLPy"  ,        &tree_genFromLLP_LLPy);
    smalltree->Branch("tree_genFromLLP_LLPz"  ,        &tree_genFromLLP_LLPz);
    smalltree->Branch("tree_genFromLLP_pt"  ,          &tree_genFromLLP_pt);
    smalltree->Branch("tree_genFromLLP_eta" ,          &tree_genFromLLP_eta);
    smalltree->Branch("tree_genFromLLP_phi" ,          &tree_genFromLLP_phi);
    smalltree->Branch("tree_genFromLLP_charge" ,       &tree_genFromLLP_charge);
    smalltree->Branch("tree_genFromLLP_pdgId" ,        &tree_genFromLLP_pdgId);
    smalltree->Branch("tree_genFromLLP_x"  ,	       &tree_genFromLLP_x);
    smalltree->Branch("tree_genFromLLP_y" ,	       &tree_genFromLLP_y);
    smalltree->Branch("tree_genFromLLP_z" ,	       &tree_genFromLLP_z);
    smalltree->Branch("tree_genFromLLP_mass" ,         &tree_genFromLLP_mass);
    smalltree->Branch("tree_genFromLLP_statusCode",    &tree_genFromLLP_statusCode);
    smalltree->Branch("tree_genFromLLP_mother_pdgId" , &tree_genFromLLP_mother_pdgId);

    smalltree->Branch("nFromC",  &nFromC,  "nFromC/I");
    smalltree->Branch("tree_genFromC_pt"  ,          &tree_genFromC_pt);
    smalltree->Branch("tree_genFromC_eta" ,          &tree_genFromC_eta);
    smalltree->Branch("tree_genFromC_phi" ,          &tree_genFromC_phi);
    smalltree->Branch("tree_genFromC_charge" ,       &tree_genFromC_charge);
    smalltree->Branch("tree_genFromC_pdgId" ,        &tree_genFromC_pdgId);
    smalltree->Branch("tree_genFromC_x"  ,	     &tree_genFromC_x);
    smalltree->Branch("tree_genFromC_y" ,	     &tree_genFromC_y);
    smalltree->Branch("tree_genFromC_z" ,	     &tree_genFromC_z);
    smalltree->Branch("tree_genFromC_mother_pdgId" , &tree_genFromC_mother_pdgId);
    smalltree->Branch("tree_genFromC_generation" ,   &tree_genFromC_generation);
    smalltree->Branch("tree_genFromC_LLP" ,          &tree_genFromC_LLP);

    smalltree->Branch("nFromB",  &nFromB,  "nFromB/I");
    smalltree->Branch("tree_genFromB_pt"  ,	     &tree_genFromB_pt);
    smalltree->Branch("tree_genFromB_eta" ,	     &tree_genFromB_eta);
    smalltree->Branch("tree_genFromB_phi" ,	     &tree_genFromB_phi);
    smalltree->Branch("tree_genFromB_charge" ,	     &tree_genFromB_charge);
    smalltree->Branch("tree_genFromB_pdgId" ,	     &tree_genFromB_pdgId);
    smalltree->Branch("tree_genFromB_x"  ,	     &tree_genFromB_x);
    smalltree->Branch("tree_genFromB_y" ,	     &tree_genFromB_y);
    smalltree->Branch("tree_genFromB_z" ,	     &tree_genFromB_z);
    smalltree->Branch("tree_genFromB_mother_pdgId" , &tree_genFromB_mother_pdgId);
    smalltree->Branch("tree_genFromB_generation" ,   &tree_genFromB_generation);
    smalltree->Branch("tree_genFromB_LLP" ,          &tree_genFromB_LLP);
    
    // genJet info
    smalltree->Branch("tree_genJet_pt"  ,   &tree_genJet_pt);
    smalltree->Branch("tree_genJet_eta" ,   &tree_genJet_eta);
    smalltree->Branch("tree_genJet_phi" ,   &tree_genJet_phi);
    smalltree->Branch("tree_genJet_mass",   &tree_genJet_mass);
    smalltree->Branch("tree_genJet_energy", &tree_genJet_energy);
    
    smalltree->Branch("tree_nLLP",&tree_nLLP);
    smalltree->Branch("tree_LLP",&tree_LLP);
    smalltree->Branch("tree_LLP_pt" ,&tree_LLP_pt);
    smalltree->Branch("tree_LLP_eta",&tree_LLP_eta);
    smalltree->Branch("tree_LLP_phi",&tree_LLP_phi);
    smalltree->Branch("tree_LLP_x",&tree_LLP_x);
    smalltree->Branch("tree_LLP_y",&tree_LLP_y);
    smalltree->Branch("tree_LLP_z",&tree_LLP_z);
    smalltree->Branch("tree_LLP_nTrks",&tree_LLP_nTrks);
    smalltree->Branch("tree_LLP_Vtx_nTrks",&tree_LLP_Vtx_nTrks);
    smalltree->Branch("tree_LLP_Vtx_NChi2",&tree_LLP_Vtx_NChi2);
    smalltree->Branch("tree_LLP_Vtx_dx",&tree_LLP_Vtx_dx);
    smalltree->Branch("tree_LLP_Vtx_dy",&tree_LLP_Vtx_dy);
    smalltree->Branch("tree_LLP_Vtx_dz",&tree_LLP_Vtx_dz);

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
    smalltree->Branch("tree_Hemi_Vtx_dx",    &tree_Hemi_Vtx_dx);
    smalltree->Branch("tree_Hemi_Vtx_dy",    &tree_Hemi_Vtx_dy);
    smalltree->Branch("tree_Hemi_Vtx_dz",    &tree_Hemi_Vtx_dz);
    smalltree->Branch("tree_Hemi_dR12",      &tree_Hemi_dR12);
    smalltree->Branch("tree_Hemi_LLP_dR12",  &tree_Hemi_LLP_dR12);

    tree_NbrOfZCand= 0;
    
    runNumber = 0;
    eventNumber = 0;
    lumiBlock = 0;
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
//$$   edm::Handle<pat::PackedCandidateCollection> pfs;
//$$   iEvent.getByToken(pfToken_, pfs);
//$$
   edm::Handle<pat::PackedCandidateCollection> pc_h;
   iEvent.getByToken(pc_, pc_h);
   const pat::PackedCandidateCollection* pc = pc_h.product();	
   edm::Handle<edm::Association<reco::PFCandidateCollection>> pc2pf;
   iEvent.getByToken(pc2pf_, pc2pf);
//$$
  iSetup.get<IdealMagneticFieldRecord>().get(bField);
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
    tree_PV_x     = (*primaryVertex)[0].x(); // l'index 0 donne le PV!
    tree_PV_y     = (*primaryVertex)[0].y();
    tree_PV_z     = (*primaryVertex)[0].z();
    tree_PV_ez    = (*primaryVertex)[0].zError();
    tree_PV_NChi2 = (*primaryVertex)[0].normalizedChi2();
    tree_PV_ndf   = (*primaryVertex)[0].ndof();
  }
  const reco::Vertex &PV = primaryVertex->front();


  //////////////////////////////////
  //////////////////////////////////
  ///////   Gen Particles   ////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_nLLP = -1;
  tree_GenPVx = -1.;
  tree_GenPVy = -1.;
  tree_GenPVz = -20.;

  int nLLP = 0;
  int nllp = 0;
  nBC = 0; 
      
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
 
    // GENPRUNED PARTICLES (initial partons and heavy hadrons)
    int genParticle_idx=0;
    for (size_t i=0; i<pruned->size(); i++)
    {
      genParticle_idx++;
      int pdgid     = (*pruned)[i].pdgId();
      float Gen_pt  = (*pruned)[i].pt();
      float Gen_eta = (*pruned)[i].eta();
      float Gen_phi = (*pruned)[i].phi();
      float Gen_m   = (*pruned)[i].mass();
      // smuon
      if ( pdgid == 1000013 ) {
	tree_GenPVx = (*pruned)[i].vx();
	tree_GenPVy = (*pruned)[i].vy();
	tree_GenPVz = (*pruned)[i].vz();
      }
      const Candidate * mom = (*pruned)[i].mother();
      
      // neutralino from smuon
      if ( pdgid == 1000023 && abs(mom->pdgId()) == 1000013 ) {
	nLLP++;
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
      if ( abs(pdgid) >= 1 && abs(pdgid) <= 6 && abs(mom->pdgId()) == 1000023 ) {
	if ( nllp >= 2 ) {
	  float dV1 = ((*pruned)[i].vx() - LLP1_x)*((*pruned)[i].vx() - LLP1_x)
	            + ((*pruned)[i].vy() - LLP1_y)*((*pruned)[i].vy() - LLP1_y)
	            + ((*pruned)[i].vz() - LLP1_z)*((*pruned)[i].vz() - LLP1_z); // dV1 is equal to dV from nllp==1
	  float dV2 = ((*pruned)[i].vx() - LLP2_x)*((*pruned)[i].vx() - LLP2_x)
	            + ((*pruned)[i].vy() - LLP2_y)*((*pruned)[i].vy() - LLP2_y)
	            + ((*pruned)[i].vz() - LLP2_z)*((*pruned)[i].vz() - LLP2_z);
	  if ( dV1 > 0.01 && dV2 > 0.01 ) nllp++; // should be == 2, so just to check : dV2 is always equal to 0 here
	}
	if ( nllp == 1 ) {
	  float dV = ((*pruned)[i].vx() - LLP1_x)*((*pruned)[i].vx() - LLP1_x)
	           + ((*pruned)[i].vy() - LLP1_y)*((*pruned)[i].vy() - LLP1_y)
	           + ((*pruned)[i].vz() - LLP1_z)*((*pruned)[i].vz() - LLP1_z);
	  if ( dV > 0.01 ) {
	    nllp = 2;
	    LLP2_x = (*pruned)[i].vx();
	    LLP2_y = (*pruned)[i].vy();
	    LLP2_z = (*pruned)[i].vz();
	    LLP2_dist = TMath::Sqrt( (LLP2_x - tree_GenPVx)*(LLP2_x - tree_GenPVx) 
				   + (LLP2_y - tree_GenPVy)*(LLP2_y - tree_GenPVy) 
				   + (LLP2_z - tree_GenPVz)*(LLP2_z - tree_GenPVz) ); 
	  }
	}
	if ( nllp == 0 ) {
	  nllp = 1;
	  LLP1_x = (*pruned)[i].vx();
	  LLP1_y = (*pruned)[i].vy();
	  LLP1_z = (*pruned)[i].vz();
	  LLP1_dist = TMath::Sqrt( (LLP1_x - tree_GenPVx)*(LLP1_x - tree_GenPVx) 
				 + (LLP1_y - tree_GenPVy)*(LLP1_y - tree_GenPVy) 
				 + (LLP1_z - tree_GenPVz)*(LLP1_z - tree_GenPVz) ); 
	}
      }
      
    if ((*pruned)[i].pt() < 0.9 || fabs((*pruned)[i].eta()) > 4.0) continue;
      
      tree_genParticle_pt.push_back(        (*pruned)[i].pt());
      tree_genParticle_eta.push_back(       (*pruned)[i].eta());
      tree_genParticle_phi.push_back(       (*pruned)[i].phi());
      tree_genParticle_charge.push_back(    (*pruned)[i].charge());
      tree_genParticle_pdgId.push_back(     (*pruned)[i].pdgId());
      tree_genParticle_x.push_back(	    (*pruned)[i].vx());
      tree_genParticle_y.push_back(	    (*pruned)[i].vy());
      tree_genParticle_z.push_back(	    (*pruned)[i].vz());
      tree_genParticle_mass.push_back(      (*pruned)[i].mass());
      tree_genParticle_statusCode.push_back((*pruned)[i].status());
      tree_genParticle_mother_pdgId.push_back( mom ? mom->pdgId() :  -1 );

    } // end loop on genpruned particles
    
    // second pass to recover the final gen particles from LLP decay
    int nLLPbis = 0;
    for (size_t i=0; i<pruned->size(); i++)     // loop on pruned genparticles
    {
      int pdgid     = (*pruned)[i].pdgId();
      const Candidate * mom = (*pruned)[i].mother();
      if ( pdgid == 1000023 && abs(mom->pdgId()) == 1000013 ) {
	nLLPbis++;
	if ( nLLPbis == 1 || nLLPbis == 2 ) {
          const Candidate * Neutralino = &(*pruned)[i];
          for (size_t j=0; j<packed->size(); j++) {     // loop on packed genparticles
  	    const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
  	    if ( motherInPrunedCollection != nullptr && isAncestor( Neutralino , motherInPrunedCollection)) {
              if (!( (*packed)[j].pt() < 0.9 || fabs((*packed)[j].eta()) > 3.0 || (*packed)[j].charge() == 0 )) {
                tree_genFromLLP_LLP.push_back(nLLPbis);
		if ( nLLPbis == 1 ) {
                  tree_genFromLLP_LLPx.push_back(LLP1_x);
                  tree_genFromLLP_LLPy.push_back(LLP1_y);
                  tree_genFromLLP_LLPz.push_back(LLP1_z);
		}
		else if ( nLLPbis == 2 ) {
                  tree_genFromLLP_LLPx.push_back(LLP2_x);
                  tree_genFromLLP_LLPy.push_back(LLP2_y);
                  tree_genFromLLP_LLPz.push_back(LLP2_z);
		}
                tree_genFromLLP_pt.push_back(	     (*packed)[j].pt());
                tree_genFromLLP_eta.push_back(       (*packed)[j].eta());
                tree_genFromLLP_phi.push_back(       (*packed)[j].phi());
                tree_genFromLLP_charge.push_back(    (*packed)[j].charge());
                tree_genFromLLP_pdgId.push_back(     (*packed)[j].pdgId());
                tree_genFromLLP_mass.push_back(      (*packed)[j].mass());
                tree_genFromLLP_statusCode.push_back((*packed)[j].status());
                const Candidate * momj =             (*packed)[j].mother(0);
                if ( momj ) {
                  tree_genFromLLP_mother_pdgId.push_back(momj->pdgId());
                  tree_genFromLLP_x.push_back(momj->vx());
                  tree_genFromLLP_y.push_back(momj->vy());
                  tree_genFromLLP_z.push_back(momj->vz());
		}
                else {
                  tree_genFromLLP_mother_pdgId.push_back(-9999);
                  tree_genFromLLP_x.push_back(-10);
                  tree_genFromLLP_y.push_back(-10);
                  tree_genFromLLP_z.push_back(-10);
                }
	      }
  	    }
          } // end loop on packed genparticles
	  if ( nLLPbis == 2 ) break;
	}
      }
    } // end loop on pruned genparticles

    tree_nLLP = nllp;

    nFromC = 0; 
    nFromB = 0;
    for (size_t i=0; i<pruned->size(); i++) // loop on genpruned particles
    {
      const GenParticle & genIt = (*pruned)[i];
      int ID = abs(genIt.pdgId());
      unsigned int nDaughters = genIt.numberOfDaughters();

      int fromLLP = -1;
      if ( (ID/100)%10 == 4 || (ID/1000)%10 == 4 || (ID/100)%10 == 5 || (ID/1000)%10 == 5 ) {
	float dV0 = (genIt.vx() - tree_GenPVx)*(genIt.vx() - tree_GenPVx)
	    	  + (genIt.vy() - tree_GenPVy)*(genIt.vy() - tree_GenPVy)
	    	  + (genIt.vz() - tree_GenPVz)*(genIt.vz() - tree_GenPVz);
	float dV1 = (genIt.vx() - LLP1_x)*(genIt.vx() - LLP1_x)
	    	  + (genIt.vy() - LLP1_y)*(genIt.vy() - LLP1_y)
	    	  + (genIt.vz() - LLP1_z)*(genIt.vz() - LLP1_z);
	float dV2 = (genIt.vx() - LLP2_x)*(genIt.vx() - LLP2_x)
	    	  + (genIt.vy() - LLP2_y)*(genIt.vy() - LLP2_y)
	    	  + (genIt.vz() - LLP2_z)*(genIt.vz() - LLP2_z);
        if      ( dV1 < dV2 && dV1 < 0.01 ) fromLLP = 1;
        else if ( dV2 < dV1 && dV2 < 0.01 ) fromLLP = 2;
        else if ( dV0 < 0.01 )              fromLLP = 0;
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
        const Candidate * Ancestor = &(*pruned)[i];
        for (size_t j=0; j<packed->size(); j++) 
        {
          //get the pointer to the first survived ancestor of a given packed GenParticle in the prunedCollection
  	  const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
  	  if ( motherInPrunedCollection != nullptr && isAncestor( Ancestor , motherInPrunedCollection))
	  {
            if (!( (*packed)[j].pt() < 0.9 || fabs((*packed)[j].eta()) > 3.0 || (*packed)[j].charge() == 0 ))
	    {
              nFromC++;
              tree_genFromC_pt.push_back(    (*packed)[j].pt());
              tree_genFromC_eta.push_back(   (*packed)[j].eta());
              tree_genFromC_phi.push_back(   (*packed)[j].phi());
              tree_genFromC_charge.push_back((*packed)[j].charge());
              tree_genFromC_pdgId.push_back( (*packed)[j].pdgId());
              tree_genFromC_mother_pdgId.push_back( genIt.pdgId());
              unsigned int nDaughters2 = genIt.numberOfDaughters();
              tree_genFromC_generation.push_back(nDaughters2);
	      if ( nDaughters2 > 0 ) {
                const Candidate* gen2 = genIt.daughter(0);
                tree_genFromC_x.push_back(gen2->vx());
                tree_genFromC_y.push_back(gen2->vy());
                tree_genFromC_z.push_back(gen2->vz());
	      }
	      else {
                tree_genFromC_x.push_back(-10);
                tree_genFromC_y.push_back(-10);
                tree_genFromC_z.push_back(-10);
	      }
              tree_genFromC_LLP.push_back(fromLLP);
  	    }
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
        const Candidate * Ancestor = &(*pruned)[i];
        for (size_t j=0; j<packed->size(); j++) 
        {
          //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
  	  const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
  	  if ( motherInPrunedCollection != nullptr && isAncestor( Ancestor , motherInPrunedCollection))
	  {
            if (!( (*packed)[j].pt() < 0.9 || fabs((*packed)[j].eta()) > 3.0 || (*packed)[j].charge() == 0 ))
	    { 
              nFromB++;
              tree_genFromB_pt.push_back(    (*packed)[j].pt());
              tree_genFromB_eta.push_back(   (*packed)[j].eta());
              tree_genFromB_phi.push_back(   (*packed)[j].phi());
              tree_genFromB_charge.push_back((*packed)[j].charge());
              tree_genFromB_pdgId.push_back( (*packed)[j].pdgId());
              tree_genFromB_mother_pdgId.push_back( genIt.pdgId());
              unsigned int nDaughters2 = genIt.numberOfDaughters();
              tree_genFromB_generation.push_back(nDaughters2);
	      if ( nDaughters2 > 0 ) {
                const Candidate* gen2 = genIt.daughter(0);
                tree_genFromB_x.push_back(gen2->vx());
                tree_genFromB_y.push_back(gen2->vy());
                tree_genFromB_z.push_back(gen2->vz());
	      }
	      else {
                tree_genFromB_x.push_back(-10);
                tree_genFromB_y.push_back(-10);
                tree_genFromB_z.push_back(-10);
	      }
              tree_genFromB_LLP.push_back(fromLLP);
  	    }
  	  }
        }
      } // final b hadron
    } // end loop on genpruned particles

    // GENPACKED PARTICLES (final particles)
    for (size_t i=0; i<packed->size(); i++) 
    {
    if ( (*packed)[i].pt() < 0.9 || fabs((*packed)[i].eta()) > 3.0 || (*packed)[i].charge() == 0 ) continue;
      const Candidate * mom = (*packed)[i].mother(0);
      tree_genPackPart_pt.push_back(        (*packed)[i].pt());
      tree_genPackPart_eta.push_back(       (*packed)[i].eta());
      tree_genPackPart_phi.push_back(       (*packed)[i].phi());
      tree_genPackPart_charge.push_back(    (*packed)[i].charge());
      tree_genPackPart_pdgId.push_back(     (*packed)[i].pdgId());
      tree_genPackPart_x.push_back(	    (*packed)[i].vx());
      tree_genPackPart_y.push_back(	    (*packed)[i].vy());
      tree_genPackPart_z.push_back(	    (*packed)[i].vz());
      tree_genPackPart_mass.push_back(      (*packed)[i].mass());
      tree_genPackPart_statusCode.push_back((*packed)[i].status());
      tree_genPackPart_mother_pdgId.push_back( mom ? mom->pdgId() :  -10 );
      tree_genPackPart_nDaughters.push_back((*packed)[i].numberOfDaughters());
    }

    // GEN JETS
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
  
  tree_nAK4jet = 0;
  float HT_val = 0;
  for (const pat::Jet &jet : *jets) {
  if ( jet.pt() < 20. ) continue;
    tree_jet_E.push_back(jet.energy());
    tree_jet_pt.push_back(jet.pt());
    tree_jet_eta.push_back(jet.eta());
    tree_jet_phi.push_back(jet.phi());
    tree_nAK4jet++;
    if ( abs(jet.eta()) < 2.4 ) HT_val += jet.pt(); // used in HT filter !
  }
  
  //////////////////////////////////
  //////////////////////////////////
  ////////   Electrons   ///////////
  //////////////////////////////////
  //////////////////////////////////
  
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

  if ( tree_Mmumu > 60. )                  tree_NbrOfZCand = 1;
  if ( tree_Mmumu > 60. && HT_val > 180. ) tree_passesHTFilter = true;


  edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);//Asking for reco colelction of PV..
  vector<reco::TransientTrack> BestTracks;
  std::vector<std::pair<uint16_t,float> > Players;
  int count =0;
  std::map<size_t , int > trackToAK4SlimmedJetMap;

  if ( tree_passesHTFilter ) {
//$$

  //////////////////////////////////
  //////////////////////////////////
  //////////   Tracks   ////////////
  //////////////////////////////////
  //////////////////////////////////
  
    // loop on pf candidates

// from /PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc
// and
// from /DQM/TrackingMonitor/src/PackedCandidateTrackValidator.cc

    for (unsigned int ipc = 0; ipc < pc->size(); ipc++) {
      // const pat::PackedCandidate& pfCand = pc->at(ipc);
      pat::PackedCandidateRef pcref = pat::PackedCandidateRef(pc_h, ipc);

      const reco::Track *trackPcPtr = pcref->bestTrack();
    if( !trackPcPtr ) continue;

      const reco::Track& trackPc = *trackPcPtr;

    //---------------------------Je to track mathcign==> isinjet BDT varaible-----//
      int idxSlimmedJet = 0;
      bool found_match = false;
      for (const pat::Jet &jet : *jets) {
      // for (unsigned int ij=0;ij<ak4slimmedJets->size();ij++) {
      //   const Jet& jet = ak4slimmedJets->at(ij);
      if ( jet.pt() < 20 ) continue;
        TLorentzVector jet4m, track4m;
        jet4m.SetPtEtaPhiM(jet.pt(), jet.eta(), jet.phi(), 0);
        track4m.SetPtEtaPhiM(trackPc.pt(), trackPc.eta(), trackPc.phi(), 0);
        if ( jet4m.DeltaR(track4m) < 0.4 ) {
          found_match = true;
          break;
        }
        else idxSlimmedJet++;
      }
      if (found_match) trackToAK4SlimmedJetMap[ipc] = idxSlimmedJet;
      else	       trackToAK4SlimmedJetMap[ipc] = -1;
      tree_track_iJet.push_back(trackToAK4SlimmedJetMap[ipc]);

      //------------------------End of Jet to track matching---------------//

    if ( trackPc.hitPattern().numberOfValidHits() == 0 ) continue;

    if ( trackPc.charge() == 0 || trackPc.pt() < 1. ) continue;

      tree_nTracks++; 
      tree_track_pt.push_back             (trackPc.pt());
      tree_track_eta.push_back            (trackPc.eta());
      tree_track_phi.push_back            (trackPc.phi());
      tree_track_charge.push_back         (trackPc.charge());
      tree_track_NChi2.push_back          (trackPc.normalizedChi2());
      tree_track_vx.push_back              (trackPc.vx());
      tree_track_vy.push_back              (trackPc.vy());
      tree_track_vz.push_back              (trackPc.vz());
      tree_track_dxy.push_back            (trackPc.dxy((*primaryVertex)[0].position()));
      tree_track_dxyError.push_back       (trackPc.dxyError());
      tree_track_dz.push_back             (trackPc.dz((*primaryVertex)[0].position()));
//       tree_track_isHighQuality               (trackPc.highPurity());//All tracks are high quality in MINIAOD
      tree_track_nHit.push_back           (trackPc.numberOfValidHits());
      tree_track_nHitPixel.push_back      (trackPc.hitPattern().numberOfValidPixelHits());
      tree_track_nHitTIB.push_back        (trackPc.hitPattern().numberOfValidStripTIBHits());
      tree_track_nHitTID.push_back        (trackPc.hitPattern().numberOfValidStripTIDHits());
      tree_track_nHitTOB.push_back        (trackPc.hitPattern().numberOfValidStripTOBHits());
      tree_track_nHitTEC.push_back        (trackPc.hitPattern().numberOfValidStripTECHits());
      tree_track_nHitPXB.push_back        (trackPc.hitPattern().numberOfValidPixelBarrelHits());
      tree_track_nHitPXF.push_back        (trackPc.hitPattern().numberOfValidPixelEndcapHits());
      tree_track_isHitL1.push_back        (trackPc.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1));
      tree_track_nLayers.push_back        (trackPc.hitPattern().trackerLayersWithMeasurement());
      tree_track_nLayersPixel.push_back   (trackPc.hitPattern().pixelLayersWithMeasurement());
      tree_track_algo.push_back(trackPc.algo());
      tree_track_originalAlgo.push_back(trackPc.originalAlgo());
      
      //!!!!
      //----------------MINIAOD_Firsthit-----------//
                  //-----------------IMPORTANT----------------//
                  // BestTracks seems to be destructed somehow//
                  // and therefore cannot be used after-------//
                  // TSOS is said to be better for the -------//
                  // propagators (see Propagator.h)...--------//
                  //------------------------------------------//
      //-----hitpattern -> Database ---/
      const HitPattern hp = trackPc.hitPattern();
      uint16_t firsthit = hp.getHitPattern(HitPattern::HitCategory::TRACK_HITS,0);
      tree_track_hitpattern.push_back(firsthit);

      //---Creating State to propagate from  TT---//
      BestTracks.push_back(theTransientTrackBuilder->build(trackPc));
      const MagneticField* B = BestTracks[count].field();//3.8T
      reco::TransientTrack TT (trackPc,BestTracks[count].field());
      // const FreeTrajectoryState Freetraj = TT.initialFreeState();//Propagator in the barrel can also use FTS (WARNING: the so-called reference point (where the propagation starts might be different from the first vtx, a check should be done))
      GlobalPoint vert (trackPc.vx(),trackPc.vy(),trackPc.vz());//Point where the propagation will start (Reference Point)
      const TrajectoryStateOnSurface Surtraj = TT.stateOnSurface(vert);//TSOS of this point
      AnalyticalPropagator* Prop = new AnalyticalPropagator(B);//Propagator that will be used for barrel, crashes in the disks when using Plane
      Basic3DVector<float> P3D2(trackPc.vx(),trackPc.vy(),trackPc.vz());//global frame
      Basic3DVector<float> B3DV (trackPc.px(),trackPc.py(),trackPc.pz());//global frame 
      float Eta = trackPc.eta();
      float Phi = trackPc.phi();
      float vz  = trackPc.vz();
      // double pz = trackPc.pz();
      //------Propagation with new interface --> See ../interface/PropaHitPattern.h-----//
      PropaHitPattern* PHP = new PropaHitPattern();
      std::pair<int,GloballyPositioned<float>::PositionType> FHPosition = PHP->Main(firsthit,Prop,Surtraj,Eta,Phi,vz,P3D2,B3DV);
      tree_track_firsthit_X.push_back(FHPosition.second.x());
      tree_track_firsthit_Y.push_back(FHPosition.second.y());
      tree_track_firsthit_Z.push_back(FHPosition.second.z());
      tree_track_region.push_back(FHPosition.first);
      count+=1;
      //-----------------------END OF MINIAOD firsthit-----------------------//
      //!!!!


    }

  } // endif HT filter (+Zmu)
//$$
  
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
    int jetidx = 0; //FIXME : May be in the loop/ not sure it changes anything
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
     
    for (int i=0; i<jetidx; i++) // Loop on jet, //FIXME
    {
    if ( !isjet[i] ) continue;
      // float jet_pt  = vjet[i].Pt();
      float jet_eta = vjet[i].Eta();
      float jet_phi = vjet[i].Phi();
      if ( njet1 > 0 ) dR1 = Deltar( jet_eta, jet_phi, vaxis1.Eta(), vaxis1.Phi() );
      if ( njet2 > 0 ) dR2 = Deltar( jet_eta, jet_phi, vaxis2.Eta(), vaxis2.Phi() );
      // axis 1
      if ( njet1 > 0 && !isjet2[i]  && dR1 < dRcut_hemis) {	     //
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
    //selection of displaced tracks
    //-----------------------------------------------------
    ///////////////////////////////////////////////////////

    vector<reco::TransientTrack> displacedTracks_llp1_mva, displacedTracks_llp2_mva; // Control Tracks
    vector<reco::TransientTrack> displacedTracks_Hemi1_mva, displacedTracks_Hemi2_mva; // Tracks selected wrt the hemisphere

    //ajoute par Paul /*!*/
    float drSig, isinjet;
    // float dptSig; /*!*/
    int jet; /*!*/
    float ntrk10, ntrk20, ntrk30; /*!*/
    float firsthit_X, firsthit_Y, firsthit_Z, dxy, dxyError, pt, eta,phi, NChi2, nhits, algo;
    // float  track_dR;
    // float track_dRmax ; /*!*/
    double bdtval = -100.;

    int nTrks_axis1 = 0;
    // int nTrks_axis1_sig=0, nTrks_axis1_bad=0;
    int nTrks_axis2 = 0;
    // int nTrks_axis2_sig=0, nTrks_axis2_bad=0;

    // int nTrks_axis1_sig_mva=0, nTrks_axis1_bad_mva=0;
    // int nTrks_axis2_sig_mva=0, nTrks_axis2_bad_mva=0;
    
    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

    // reader->AddVariable( "mva_track_firstHit_x", &firsthit_X );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_y", &firsthit_Y );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_z", &firsthit_Z );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_dxy", &dxy );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_dxyError", &dxyError );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_dz", &dz );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    // reader->AddVariable( "mva_track_firstHit_dzError", &dzError );//to be exluded if TMVAbgctau50withnhits.xml is chosen
    reader->AddVariable( "mva_track_pt", &pt );
    reader->AddVariable( "mva_track_eta", &eta );
    reader->AddVariable( "mva_track_nchi2", &NChi2 );
    reader->AddVariable( "mva_track_nhits", &nhits );
    // reader->AddVariable( "mva_track_algo", &algo);
    //add the variables from my BDT(Paul)
    reader->AddVariable( "mva_ntrk10", &ntrk10);
    reader->AddVariable( "mva_drSig", &drSig); /*!*/
    reader->AddVariable( "mva_track_isinjet", &isinjet); /*!*/
    //reader->AddVariable("mva_track_dR",&track_dR);
    // reader->AddVariable("mva_track_dRmax",&track_dRmax);

    reader->BookMVA( "BDTG", weightFile_ );//root 6.14/09, care compatiblity of versions for tmva

    //$$
    float pt_Cut = 1.;
    float NChi2_Cut = 5.;
    float drSig_Cut = 5.;
    double bdtcut = -0.0815; 
    //optimal value : -0.0815 for TMVAClassification_BDTG50sansalgo.weights.xml
    // optimal value 0.0057 for TMVAClassification_BDTG50cm.weights.xml
    // optimal value w/o track association to axis: -0.0401: TMVAbgctau50withnhits.xml

//$$    double bdtcut = -10.; // if NO BDT cut !**

    int counter_track = -1;
    //---------------------------//
    for (unsigned int ipc = 0; ipc < pc->size(); ipc++) 
    {
      pat::PackedCandidateRef pcref = pat::PackedCandidateRef(pc_h, ipc);
      const reco::Track *trackPcPtr = pcref->bestTrack();
      if( !trackPcPtr ) continue;//happens (~5%) 
            const reco::Track& trackPc = *trackPcPtr;
      if ( trackPc.hitPattern().numberOfValidHits() == 0 ) continue;

      if ( trackPc.charge() == 0 || trackPc.pt() < 1. ) continue;
    // for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack)  // Loop on all the tracks
    // {
      counter_track++;
      // const auto& itTrack = trackRefs[iTrack];
      firsthit_X = tree_track_firsthit_X[counter_track];
      firsthit_Y = tree_track_firsthit_Y[counter_track];
      firsthit_Z = tree_track_firsthit_Z[counter_track];
      dxy	 = tree_track_dxy[counter_track];
      dxyError   = tree_track_dxyError[counter_track];
      pt	 = tree_track_pt[counter_track];
      // ptError    = tree_track_ptError[counter_track];
      eta	 = tree_track_eta[counter_track];
      phi	 = tree_track_phi[counter_track];
      NChi2	 = tree_track_NChi2[counter_track];
      nhits	 = tree_track_nHit[counter_track];
      algo	 = tree_track_algo[counter_track];
      // float dptSig=-1;
      // if (pt>0) dptSig=ptError/pt;
      
      //Ajoute par Paul /*!*/
      drSig = -1.;
      if ( dxyError > 0 ) drSig = abs(dxy) / dxyError; /*!*/
//       tree_track_drSig.push_back(drSig);
      ntrk10 = 0;
      ntrk20 = 0;
      ntrk30 = 0;
      bdtval = -10.;
      dR = -1.;
      int tracks_axis = 0; // flag to check which axis is the closest from the track
      algo=algo+1;
//$$
      if ( pt > pt_Cut && NChi2 < NChi2_Cut && drSig > drSig_Cut ) // preselection : pt > 1. && NChi2 < 5. && drSig > 5.
//$$
      { // On regarde les track_selec[i] qui sont True donc de potentielles tracks secondaires
        

        //-----------------------------------------------FIXME---------------//
        jet = tree_track_iJet[counter_track]; /*!*/
        isinjet = 0.;
        if ( jet >= 0 ) isinjet = 1.; /*!*/
        // int isFromLLP    = tree_track_sim_LLP[counter_track];//FIXME
        //--------------------------------------------------------------------//


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
        int counter_othertrack = -1; ///*!*/
        // for (size_t iTrack2=0; iTrack2<trackRefs.size(); ++iTrack2)    // Loop on all the other Tracks/*!*/
        // {
              for (unsigned int ipc2 = 0; ipc2 < pc->size(); ipc2++) 
    {
          counter_othertrack++;
          if ( counter_othertrack == counter_track ) continue;
          float pt2  = tree_track_pt[counter_othertrack];
          float dr2 = abs(tree_track_dxy[counter_othertrack]);
          float drSig2 = -1.;
          if ( tree_track_dxyError[counter_othertrack] > 0 ) drSig2 = dr2 / tree_track_dxyError[counter_othertrack];
          NChi2 = tree_track_NChi2[counter_othertrack];
//$$
        if ( !(pt2 > pt_Cut && NChi2 < NChi2_Cut && drSig2 > drSig_Cut ) ) continue; // On regarde les autres track_selec[i] qui sont True donc de potnetielles tracks secondaires
//$$
          float x2 = tree_track_firsthit_X[counter_othertrack];
          float y2 = tree_track_firsthit_Y[counter_othertrack];
          float z2 = tree_track_firsthit_Z[counter_othertrack];
          float dist = TMath::Sqrt( (firsthit_X-x2)*(firsthit_X-x2) + (firsthit_Y-y2)*(firsthit_Y-y2) + (firsthit_Z-z2)*(firsthit_Z-z2) );//pour chaque reconstruite, on regarde les autres tracks,
          if ( dist < 10. )	     {ntrk10++;} // les sctocker les 3 , on teste sur une seule couche quand on regarde vers l'avant
          if ( dist < 20. )	     {ntrk20++;}
          if ( dist < 30. )	     {ntrk30++;}
        }  // end Loop on other Tracks

        if ( dR < dRcut_tracks ) 
	{
	  // if ( isFromLLP == 1 ) LLP1_nTrks++;
	  // if ( isFromLLP == 2 ) LLP2_nTrks++;
	
          bdtval = reader->EvaluateMVA( "BDTG" ); //default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
          //cout << "BDT VAL " << bdtval <<endl;

      //     if ( tracks_axis == 1 ) {
	    // nTrks_axis1++;
      //       if ( isFromLLP == iLLPrec1 ) nTrks_axis1_sig++;
      //       else if ( isFromLLP >= 1 )   nTrks_axis1_bad++;
      //     }
        
      //     if ( tracks_axis == 2 ) {
	    // nTrks_axis2++;
      //       if ( isFromLLP == iLLPrec2 ) nTrks_axis2_sig++;
      //       else if ( isFromLLP >= 1 )   nTrks_axis2_bad++;
      //     }
        
          if ( bdtval > bdtcut )      //optimal cut for 50 cm ,trained on 10k events //Paul, expected to change
          { 
            ////--------------Control tracks-----------------////
            // if ( isFromLLP == 1 )
            // {
            //   displacedTracks_llp1_mva.push_back(theTransientTrackBuilder->build(trackPcPtr));
            // }
            // if ( isFromLLP == 2 )
            // {
            //   displacedTracks_llp2_mva.push_back(theTransientTrackBuilder->build(trackPcPtr));
            // }

            if ( tracks_axis == 1 )
            {
              displacedTracks_Hemi1_mva.push_back(theTransientTrackBuilder->build(trackPcPtr));
              // if ( isFromLLP == iLLPrec1 ) nTrks_axis1_sig_mva++;
              // else if ( isFromLLP >= 1 )   nTrks_axis1_bad_mva++;
            }
 
            if ( tracks_axis == 2 )
            {
              displacedTracks_Hemi2_mva.push_back(theTransientTrackBuilder->build(trackPcPtr));
              // if ( isFromLLP == iLLPrec2 ) nTrks_axis2_sig_mva++;
              // else if ( isFromLLP >= 1 )   nTrks_axis2_bad_mva++;
            }
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

        // cout << " displaced tracks LLP1 " << LLP1_nTrks << " and with mva" << displacedTracks_llp1_mva.size() << endl;
        // cout << " displaced tracks LLP2 " << LLP2_nTrks << " and with mva" << displacedTracks_llp2_mva.size() << endl;
        // cout << " displaced tracks Hemi1 " << nTrks_axis1 << " and with mva" << displacedTracks_Hemi1_mva.size() << endl;
        // cout << " displaced tracks Hemi2 " << nTrks_axis2 << " and with mva" << displacedTracks_Hemi2_mva.size() << endl;

    
    //---------------------------------------------------------------------------------------//

    ///////////////////////////////////////////////////////
    //-----------------------------------------------------
    // Vertex fitting 
    //-----------------------------------------------------
    ///////////////////////////////////////////////////////

    int   Vtx_ntk = 0;
    float Vtx_x = 0., Vtx_y = 0., Vtx_z= 0., Vtx_chi = 0.;
    
// $$
    // parameters for the Adaptive Vertex Fitter (AVF)
    double maxshift        = 0.0001;
    unsigned int maxstep   = 30;
    double maxlpshift      = 0.1;
    double weightThreshold = 0.001;
    double sigmacut        = 5.;
    double Tini            = 256.;
    double ratio           = 0.25;
// $$

//$$    KalmanVertexFitter theFitter_Vertex_Hemi1_mva(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
//$$
    static AdaptiveVertexFitter 
    theFitter_Vertex_Hemi1_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_Vertex_Hemi1_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
//$$
    
    Vtx_ntk = displacedTracks_Hemi1_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -1.;
	
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_Hemi1_mva = theFitter_Vertex_Hemi1_mva.vertex(displacedTracks_Hemi1_mva); // fitted vertex
      if ( displacedVertex_Hemi1_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      { 
        Vtx_x = displacedVertex_Hemi1_mva.position().x();
        Vtx_y = displacedVertex_Hemi1_mva.position().y();
        Vtx_z = displacedVertex_Hemi1_mva.position().z();
        Vtx_chi = displacedVertex_Hemi1_mva.normalisedChiSquared();
        Vtx_ntk = displacedTracks_Hemi1_mva.size();
      }
    }
   
    float Vtx_chi1 = Vtx_chi;
    tree_Hemi.push_back(1);
    tree_Hemi_njet.push_back(njet1);
    tree_Hemi_eta.push_back(axis1_eta);
    tree_Hemi_phi.push_back(axis1_phi);
    tree_Hemi_dR.push_back(axis1_dR);
    tree_Hemi_nTrks.push_back(nTrks_axis1);
    // tree_Hemi_nTrks_sig.push_back(nTrks_axis1_sig);
    // tree_Hemi_nTrks_bad.push_back(nTrks_axis1_bad);
    tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
    tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
    // tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis1_sig_mva);
    // tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis1_bad_mva);
    tree_Hemi_Vtx_x.push_back(Vtx_x);
    tree_Hemi_Vtx_y.push_back(Vtx_y);
    tree_Hemi_Vtx_z.push_back(Vtx_z);
    if ( iLLPrec1 == 1 ) {
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP1_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP1_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP1_z);
    }
    else {
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP2_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP2_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP2_z);
    }

    float dist = 0.;
    if ( iLLPrec1 == 1 ) {
      dist = (LLP1_x - tree_GenPVx)*(LLP1_x - tree_GenPVx)
           + (LLP1_y - tree_GenPVy)*(LLP1_y - tree_GenPVy)
           + (LLP1_z - tree_GenPVz)*(LLP1_z - tree_GenPVz);
      tree_Hemi_LLP_pt.push_back( LLP1_pt);
      tree_Hemi_LLP_eta.push_back(LLP1_eta);
      tree_Hemi_LLP_phi.push_back(LLP1_phi);
      tree_Hemi_LLP_x.push_back(LLP1_x);
      tree_Hemi_LLP_y.push_back(LLP1_y);
      tree_Hemi_LLP_z.push_back(LLP1_z);
    }
    else {
      dist = (LLP2_x - tree_GenPVx)*(LLP2_x - tree_GenPVx)
           + (LLP2_y - tree_GenPVy)*(LLP2_y - tree_GenPVy)
           + (LLP2_z - tree_GenPVz)*(LLP2_z - tree_GenPVz);
      tree_Hemi_LLP_pt.push_back( LLP2_pt);
      tree_Hemi_LLP_eta.push_back(LLP2_eta);
      tree_Hemi_LLP_phi.push_back(LLP2_phi);
      tree_Hemi_LLP_x.push_back(LLP2_x);
      tree_Hemi_LLP_y.push_back(LLP2_y);
      tree_Hemi_LLP_z.push_back(LLP2_z);
    }
    tree_Hemi_LLP.push_back(iLLPrec1);
    tree_Hemi_LLP_dist.push_back(TMath::Sqrt(dist));
    

    //--------------------------- SECOND HEMISPHERE WITH MVA -------------------------------------//
    
//$$    KalmanVertexFitter theFitter_Vertex_Hemi2_mva(kvfPSet); //One or less vertex: either Valid or NotValid /*!*/
//$$
    static AdaptiveVertexFitter 
    theFitter_Vertex_Hemi2_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_Vertex_Hemi2_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
//$$
    
    Vtx_ntk = displacedTracks_Hemi2_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -1.;
    
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_Hemi2_mva = theFitter_Vertex_Hemi2_mva.vertex(displacedTracks_Hemi2_mva); // fitted vertex
      
      if ( displacedVertex_Hemi2_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      {
        Vtx_x = displacedVertex_Hemi2_mva.position().x();
        Vtx_y = displacedVertex_Hemi2_mva.position().y();
        Vtx_z = displacedVertex_Hemi2_mva.position().z();
        Vtx_chi = displacedVertex_Hemi2_mva.normalisedChiSquared();
      }
    }
    
    float Vtx_chi2 = Vtx_chi;
    tree_Hemi.push_back(2);
    tree_Hemi_njet.push_back(njet2);
    tree_Hemi_eta.push_back(axis2_eta);
    tree_Hemi_phi.push_back(axis2_phi);
    tree_Hemi_dR.push_back(axis2_dR);
    tree_Hemi_nTrks.push_back(nTrks_axis2);
    // tree_Hemi_nTrks_sig.push_back(nTrks_axis2_sig);
    // tree_Hemi_nTrks_bad.push_back(nTrks_axis2_bad);
    tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
    tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
    // tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis2_sig_mva);
    // tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis2_bad_mva);
    tree_Hemi_Vtx_x.push_back(Vtx_x);
    tree_Hemi_Vtx_y.push_back(Vtx_y);
    tree_Hemi_Vtx_z.push_back(Vtx_z);
    if ( iLLPrec2 == 1 ) {
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP1_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP1_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP1_z);
    }
    else {
      tree_Hemi_Vtx_dx.push_back(Vtx_x - LLP2_x);
      tree_Hemi_Vtx_dy.push_back(Vtx_y - LLP2_y);
      tree_Hemi_Vtx_dz.push_back(Vtx_z - LLP2_z);
    }

    dist = 0.;
    if ( iLLPrec2 == 1 ) {
      dist = (LLP1_x - tree_GenPVx)*(LLP1_x - tree_GenPVx)
           + (LLP1_y - tree_GenPVy)*(LLP1_y - tree_GenPVy)
           + (LLP1_z - tree_GenPVz)*(LLP1_z - tree_GenPVz);
      tree_Hemi_LLP_pt.push_back( LLP1_pt);
      tree_Hemi_LLP_eta.push_back(LLP1_eta);
      tree_Hemi_LLP_phi.push_back(LLP1_phi);
      tree_Hemi_LLP_x.push_back(LLP1_x);
      tree_Hemi_LLP_y.push_back(LLP1_y);
      tree_Hemi_LLP_z.push_back(LLP1_z);
    }
    else {
      dist = (LLP2_x - tree_GenPVx)*(LLP2_x - tree_GenPVx)
           + (LLP2_y - tree_GenPVy)*(LLP2_y - tree_GenPVy)
           + (LLP2_z - tree_GenPVz)*(LLP2_z - tree_GenPVz);
      tree_Hemi_LLP_pt.push_back( LLP2_pt);
      tree_Hemi_LLP_eta.push_back(LLP2_eta);
      tree_Hemi_LLP_phi.push_back(LLP2_phi);
      tree_Hemi_LLP_x.push_back(LLP2_x);
      tree_Hemi_LLP_y.push_back(LLP2_y);
      tree_Hemi_LLP_z.push_back(LLP2_z);
    }
    tree_Hemi_LLP.push_back(iLLPrec2);
    tree_Hemi_LLP_dist.push_back(TMath::Sqrt(dist));
    

    //     for (unsigned int ipc = 0; ipc < pc->size(); ipc++) 
    // {
    //   pat::PackedCandidateRef pcref = pat::PackedCandidateRef(pc_h, ipc);
    //   const reco::Track *trackPcPtr = pcref->bestTrack();

    counter_track = -1;
    // for (size_t iTrack = 0; iTrack<trackRefs.size(); ++iTrack) {
              for (unsigned int ipc = 0; ipc < pc->size(); ipc++) // Loop on all the tracks
    {       
            pat::PackedCandidateRef pcref = pat::PackedCandidateRef(pcs, ipc);
      const reco::Track *trackPcPtr = pcref->bestTrack();
      if( !trackPcPtr ) continue;//happens (~5%) 
      const reco::Track& trackPc = *trackPcPtr;
      if ( trackPc.hitPattern().numberOfValidHits() == 0 ) continue;

      if ( trackPc.charge() == 0 || trackPc.pt() < 1. ) continue; 
      counter_track++;
      int hemi      = tree_track_Hemi[counter_track];
      double MVAval = tree_track_MVAval[counter_track];
      Vtx_chi = -1000.;
      if      ( hemi == 1 && MVAval > bdtcut ) Vtx_chi = Vtx_chi1;
      else if ( hemi == 2 && MVAval > bdtcut ) Vtx_chi = Vtx_chi2;
      tree_track_Hemi_mva_NChi2.push_back(Vtx_chi);
    } //End loop on all the tracks
    
    tree_Hemi_dR12.push_back(dR_axis12);
    tree_Hemi_dR12.push_back(dR_axis12);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);

  //////////////////////////////////

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
    
    tree_track_charge.clear();
    tree_track_iJet.clear();
    tree_track_pt.clear();
    tree_track_outerPt.clear();
    tree_track_eta.clear();
    tree_track_phi.clear();
    tree_track_charge.clear();
    tree_track_NChi2.clear();
    tree_track_isHighQuality.clear();
    tree_track_isLoose.clear();
    tree_track_isTight.clear();
    tree_track_dxy.clear();
    tree_track_dxyError.clear();
    tree_track_dz.clear();
    tree_track_dzError.clear();
    tree_track_numberOfLostHits.clear();
    tree_track_originalAlgo.clear();
    tree_track_algo.clear();
    tree_track_nHit.clear();
    tree_track_nHitPixel.clear();
    tree_track_nHitTIB.clear();
    tree_track_nHitTID.clear();
    tree_track_nHitTOB.clear();
    tree_track_nHitTEC.clear();
    tree_track_nHitPXB.clear();
    tree_track_nHitPXF.clear();
    tree_track_isHitL1.clear();
    tree_track_nLayers.clear();
    tree_track_nLayersPixel.clear();
    tree_track_stripTECLayersWithMeasurement .clear();
    tree_track_stripTIBLayersWithMeasurement.clear();
    tree_track_stripTIDLayersWithMeasurement.clear();
    tree_track_stripTOBLayersWithMeasurement.clear();
    tree_track_hitpattern.clear();
    tree_track_vx.clear();
    tree_track_vy.clear();
    tree_track_vz.clear();
    tree_track_firsthit_X.clear();
    tree_track_firsthit_Y.clear();
    tree_track_firsthit_Z.clear();
    tree_track_region.clear();
    tree_track_ntrk10.clear();
    tree_track_ntrk20.clear();
    tree_track_ntrk30.clear();
    
    //!!!!

    tree_track_GeoBarrel_firsthit_x.clear();
    tree_track_GeoBarrel_firsthit_y.clear();
    tree_track_GeoBarrel_firsthit_z.clear();

    tree_track_PropBarrel_firsthit_x.clear();
    tree_track_PropBarrel_firsthit_y.clear();
    tree_track_PropBarrel_firsthit_z.clear();

    tree_track_GeoDisk_firsthit_x.clear();
    tree_track_GeoDisk_firsthit_y.clear();
    tree_track_GeoDisk_firsthit_z.clear();

    tree_track_PropDisk_SLPC_firsthit_x.clear();
    tree_track_PropDisk_SLPC_firsthit_y.clear();
    tree_track_PropDisk_SLPC_firsthit_z.clear();
    //!!!!

    tree_track_MVAval.clear();
    tree_track_Hemi.clear();
    tree_track_Hemi_dR.clear();
    tree_track_Hemi_mva_NChi2.clear();
    tree_track_Hemi_LLP.clear();
    
    tree_track_recoVertex_idx.clear();
    tree_track_recoAK4SlimmedJet_idx.clear();
    
    tree_track_nSimHits.clear();
    tree_track_isSimMatched.clear();
    
    tree_track_sim_charge.clear();
    tree_track_sim_pt.clear();
    tree_track_sim_eta.clear();
    tree_track_sim_phi.clear();
    tree_track_sim_longLived.clear();
    //tree_track_sim_matchedHit .clear();
    tree_track_sim_pdgId.clear();
    tree_track_sim_numberOfTrackerHits.clear();
    tree_track_sim_numberOfTrackerLayers.clear();
    tree_track_sim_mass.clear();
    tree_track_sim_status.clear();
    
    tree_track_sim_vx.clear();
    tree_track_sim_vy.clear();
    tree_track_sim_vz.clear();
    tree_track_sim_isFromLLP.clear();
    tree_track_sim_isFromBC.clear();
    tree_track_sim_isFromBC_mother_pdgId.clear();
    tree_track_sim_isFromBC_LLP.clear();
    tree_track_sim_dV.clear();
    tree_track_sim_dist.clear();
    
    tree_genAxis_dRneuneu.clear();
    tree_genParticle_pt.clear();
    tree_genParticle_eta.clear();
    tree_genParticle_phi.clear();
    tree_genParticle_charge.clear();
    tree_genParticle_pdgId.clear();
    tree_genParticle_x.clear();
    tree_genParticle_y.clear();
    tree_genParticle_z.clear();
    tree_genParticle_mass.clear();
    tree_genParticle_statusCode.clear();
    tree_genParticle_mother_pdgId.clear();

    tree_genPackPart_pt.clear();
    tree_genPackPart_eta.clear();
    tree_genPackPart_phi.clear();
    tree_genPackPart_charge.clear();
    tree_genPackPart_pdgId.clear();
    tree_genPackPart_x.clear();
    tree_genPackPart_y.clear();
    tree_genPackPart_z.clear();
    tree_genPackPart_mass.clear();
    tree_genPackPart_statusCode.clear();
    tree_genPackPart_mother_pdgId.clear();
    tree_genPackPart_nDaughters.clear();

    tree_genFromLLP_LLP.clear();
    tree_genFromLLP_LLPx.clear();
    tree_genFromLLP_LLPy.clear();
    tree_genFromLLP_LLPz.clear();
    tree_genFromLLP_pt.clear();
    tree_genFromLLP_eta.clear();
    tree_genFromLLP_phi.clear();
    tree_genFromLLP_charge.clear();
    tree_genFromLLP_pdgId.clear();
    tree_genFromLLP_x.clear();
    tree_genFromLLP_y.clear();
    tree_genFromLLP_z.clear();
    tree_genFromLLP_mass.clear();
    tree_genFromLLP_statusCode.clear();
    tree_genFromLLP_mother_pdgId.clear();

    tree_genFromC_pt.clear();
    tree_genFromC_eta.clear();
    tree_genFromC_phi.clear();
    tree_genFromC_charge.clear();
    tree_genFromC_pdgId.clear();
    tree_genFromC_x.clear();
    tree_genFromC_y.clear();
    tree_genFromC_z.clear();
    tree_genFromC_mother_pdgId.clear();
    tree_genFromC_generation.clear();
    tree_genFromC_LLP.clear();

    tree_genFromB_pt.clear();
    tree_genFromB_eta.clear();
    tree_genFromB_phi.clear();
    tree_genFromB_charge.clear();
    tree_genFromB_pdgId.clear();
    tree_genFromB_x.clear();
    tree_genFromB_y.clear();
    tree_genFromB_z.clear();
    tree_genFromB_mother_pdgId.clear();
    tree_genFromB_generation.clear();
    tree_genFromB_LLP.clear();

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
    tree_LLP_nTrks.clear();
    tree_LLP_Vtx_NChi2.clear();
    tree_LLP_Vtx_nTrks.clear();
    tree_LLP_Vtx_dx.clear();
    tree_LLP_Vtx_dy.clear();
    tree_LLP_Vtx_dz.clear();

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
    tree_Hemi_Vtx_dx.clear();
    tree_Hemi_Vtx_dy.clear();
    tree_Hemi_Vtx_dz.clear();
    tree_Hemi_dR12.clear();
    tree_Hemi_LLP_dR12.clear();
}
