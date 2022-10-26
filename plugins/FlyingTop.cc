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

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FlyingTop/FlyingTop/interface/Proto.h"
#include "FlyingTop/FlyingTop/interface/DeltaFunc.h"


//$$
// #include "CondFormats/GeometryObjects/interface/HcalParameters.h"
// #include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
// #include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
// 
// #include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
// #include "DataFormats/Candidate/interface/Candidate.h"
// #include "DataFormats/Common/interface/MergeableCounter.h"
// #include "DataFormats/Common/interface/TriggerResults.h"
// #include "DataFormats/Common/interface/View.h"
// #include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
// #include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
// #include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
// #include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
// #include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
// #include "DataFormats/HLTReco/interface/TriggerEvent.h"
// #include "DataFormats/HLTReco/interface/TriggerObject.h"
// #include "DataFormats/JetReco/interface/PFJet.h"
// #include "DataFormats/JetReco/interface/PFJetCollection.h"
// #include "DataFormats/L1Trigger/interface/Jet.h"
// #include "DataFormats/L1Trigger/interface/Tau.h"
// #include "DataFormats/Math/interface/deltaR.h"
// #include "DataFormats/Math/interface/LorentzVector.h"
// #include "DataFormats/MuonReco/interface/MuonSelectors.h"
// #include "DataFormats/METReco/interface/CaloMET.h"
// #include "DataFormats/METReco/interface/CaloMETFwd.h"
// #include "DataFormats/METReco/interface/CaloMETCollection.h"
// #include "DataFormats/METReco/interface/MET.h"
// #include "DataFormats/METReco/interface/METFwd.h"
// #include "DataFormats/METReco/interface/PFMET.h"
// #include "DataFormats/METReco/interface/PFMETFwd.h"
// #include "DataFormats/METReco/interface/PFMETCollection.h"
// #include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
// #include "DataFormats/PatCandidates/interface/GenericParticle.h"
// #include "DataFormats/PatCandidates/interface/Tau.h"
// #include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
// #include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
// #include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
// #include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
// #include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
// #include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
// #include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
// #include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
// #include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
// #include "DataFormats/TrackReco/interface/Track.h"
// #include "DataFormats/TrackReco/interface/TrackFwd.h"
// 
// #include "DetectorDescription/Core/interface/DDCompactView.h"
// 
// #include "FWCore/Common/interface/TriggerNames.h"
// #include "FWCore/Framework/interface/eventsetuprecord_registration_macro.h"
// #include "FWCore/Framework/interface/EventSetup.h"
// #include "FWCore/Framework/interface/ESHandle.h"
// #include "FWCore/Framework/interface/ESTransientHandle.h"
// #include "FWCore/Framework/interface/LuminosityBlock.h"
// #include "FWCore/Utilities/interface/InputTag.h"
// #include "FWCore/Utilities/interface/RegexMatch.h"
// 
// #include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"
// 
// #include "Geometry/HcalCommonData/interface/HcalParametersFromDD.h"
// #include "Geometry/Records/interface/HcalParametersRcd.h"
// #include "Geometry/Records/interface/IdealGeometryRecord.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
// #include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
// #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
// 
// #include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
// 
// #include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
// #include "JetMETCorrections/Objects/interface/JetCorrector.h"
// #include "JetMETCorrections/Modules/interface/JetResolution.h"
// 
// #include "RecoLocalTracker/ClusterParameterEstimator/interface/StripClusterParameterEstimator.h"
// #include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
// #include "RecoLocalTracker/SiStripClusterizer/interface/SiStripClusterInfo.h"
// 
// #include "RecoTracker/FinalTrackSelectors/plugins/getBestVertex.h"
// 
// #include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
// #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
// #include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
// #include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
// #include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
// #include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
// #include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
// 
// #include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
// #include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
// #include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
// #include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
// #include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
// #include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
// 
// #include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"
// #include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
// #include "TrackingTools/IPTools/interface/IPTools.h"
// #include "TrackingTools/PatternTools/interface/Trajectory.h"
// #include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
// #include "TrackingTools/Records/interface/TransientTrackRecord.h"
// #include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
// #include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
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
#include "TrackingTools/GeomPropagators/interface/SmartPropagator.h"
              //-------------Surfaces---------------//
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
              //----------------?-----------------//
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
//------------------------------End of Paul------------------------//
//$$

//
// class declaration
//

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
    
    edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
    edm::EDGetTokenT<pat::METCollection> metToken_;
    edm::EDGetTokenT<pat::JetCollection> jetToken_;
    edm::EDGetTokenT<edm::View<reco::GenJet> > genJetToken_;
    edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
    edm::EDGetTokenT<pat::MuonCollection> muonToken_;
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
    
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
    std::vector< float > tree_track_pt;
    std::vector< float > tree_track_outerPt;
    std::vector< float > tree_track_eta;
    std::vector< float > tree_track_phi;
    std::vector<int>     tree_track_charge;
    std::vector<int>     tree_track_nhits;
    std::vector<float >  tree_track_NChi2;
    std::vector<bool >   tree_track_isHighQuality;
    std::vector<bool >   tree_track_isLoose;
    std::vector<bool >   tree_track_isTight;
    std::vector< float>  tree_track_dxy; // dxy with respect to beam spot position
    std::vector< float>  tree_track_dxyError;
    std::vector< float>  tree_track_dz;
    std::vector< float>  tree_track_dzError ;
    std::vector<int>     tree_track_numberOfLostHits;
    std::vector<int>     tree_track_numberOfValidHits;
    std::vector<unsigned int>    tree_track_originalAlgo; // definition as comments at the end of the file,
    //http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_1_3/doc/html/d8/df2/classreco_1_1TrackBase.html#aca7611bd1a33d535cefc72b6e497ece8
    std::vector<unsigned int>    tree_track_algo;
    std::vector<unsigned short>  tree_track_stopReason;
    
    std::vector<int>     tree_track_numberOfValidPixelHits;
    std::vector<int>     tree_track_numberOfValidStripHits;
    std::vector<int>     tree_track_numberOfValidStripTIBHits;
    std::vector<int>     tree_track_numberOfValidStripTIDHits;
    std::vector<int>     tree_track_numberOfValidStripTOBHits;
    std::vector<int>     tree_track_numberOfValidStripTECHits;
    std::vector<int>     tree_track_numberOfValidPixelBarrelHits;
    std::vector<int>     tree_track_numberOfValidPixelEndcapHits;
    std::vector<int>     tree_track_hasValidHitInPixelLayer;
    std::vector<int>     tree_track_trackerLayersWithMeasurement;
    std::vector<int>     tree_track_pixelLayersWithMeasurement;
    std::vector<int>     tree_track_stripTECLayersWithMeasurement ;
    std::vector<int>     tree_track_stripTIBLayersWithMeasurement;
    std::vector<int>     tree_track_stripTIDLayersWithMeasurement;
    std::vector<int>     tree_track_stripTOBLayersWithMeasurement;
    
    std::vector< float >  tree_track_vx;
    std::vector< float >  tree_track_vy;
    std::vector< float >  tree_track_vz;
    std::vector<float>    tree_track_firsthit_X;
    std::vector<float>    tree_track_firsthit_Y;
    std::vector<float>    tree_track_firsthit_Z;
    std::vector<float>    tree_track_ntrk10;
    std::vector<float>    tree_track_ntrk20;
    std::vector<float>    tree_track_ntrk30;
    
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
    
    std::vector< float > tree_genParticle_pt;
    std::vector< float > tree_genParticle_eta;
    std::vector< float > tree_genParticle_phi;
    std::vector< float > tree_genParticle_charge;
    std::vector< int >   tree_genParticle_pdgId;
    std::vector< float > tree_genParticle_x;
    std::vector< float > tree_genParticle_y;
    std::vector< float > tree_genParticle_z;
    std::vector< float > tree_genParticle_mass;
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
    prunedGenToken_(consumes<edm::View<reco::GenParticle> >(     iConfig.getParameter<edm::InputTag>("genpruned"))),
    packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("genpacked"))),
    vertexToken_(   consumes<reco::VertexCollection>(            iConfig.getParameter<edm::InputTag>("vertices"))),
    metToken_(      consumes<pat::METCollection>(                iConfig.getParameter<edm::InputTag>("mets"))),
    jetToken_(      consumes<pat::JetCollection>(                iConfig.getParameter<edm::InputTag>("jets"))),
    genJetToken_(   consumes<edm::View<reco::GenJet>>(           iConfig.getParameter<edm::InputTag>("genjets"))),
    electronToken_( consumes<pat::ElectronCollection>(           iConfig.getParameter<edm::InputTag>("electrons"))),
    muonToken_(     consumes<pat::MuonCollection>(               iConfig.getParameter<edm::InputTag>("muons"))),
    pfToken_(       consumes<pat::PackedCandidateCollection>(    iConfig.getParameter<edm::InputTag>("pfCands")))
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
    smalltree->Branch("tree_track_pt",            &tree_track_pt);
    smalltree->Branch("tree_track_outerPt",           &tree_track_outerPt );
    smalltree->Branch("tree_track_eta",                  &tree_track_eta );
    smalltree->Branch("tree_track_phi",                  &tree_track_phi );
    smalltree->Branch("tree_track_charge",            &tree_track_charge );
    smalltree->Branch("tree_track_nhits",                &tree_track_nhits);
    smalltree->Branch("tree_track_NChi2",                &tree_track_NChi2);
    smalltree->Branch("tree_track_isHighQuality",     &tree_track_isHighQuality);
    smalltree->Branch("tree_track_isLoose",        &tree_track_isLoose);
    smalltree->Branch("tree_track_isTight",        &tree_track_isTight);
    smalltree->Branch("tree_track_dxy",                  &tree_track_dxy );
    smalltree->Branch("tree_track_dxyError",        &tree_track_dxyError);
    smalltree->Branch("tree_track_dz",            &tree_track_dz);
    smalltree->Branch("tree_track_dzError",        &tree_track_dzError  );
    smalltree->Branch("tree_track_numberOfLostHits",  &tree_track_numberOfLostHits);
    smalltree->Branch("tree_track_numberOfValidHits", &tree_track_numberOfValidHits);
    smalltree->Branch("tree_track_originalAlgo",      &tree_track_originalAlgo);
    smalltree->Branch("tree_track_algo",             &tree_track_algo);
    smalltree->Branch("tree_track_stopReason",        &tree_track_stopReason);
    smalltree->Branch("tree_track_isSimMatched",      &tree_track_isSimMatched    );
    
    smalltree->Branch("tree_track_numberOfValidPixelHits",        &tree_track_numberOfValidPixelHits);
    smalltree->Branch("tree_track_numberOfValidStripHits",        &tree_track_numberOfValidStripHits);
    smalltree->Branch("tree_track_numberOfValidStripTIBHits",     &tree_track_numberOfValidStripTIBHits);
    smalltree->Branch("tree_track_numberOfValidStripTIDHits",     &tree_track_numberOfValidStripTIDHits);
    smalltree->Branch("tree_track_numberOfValidStripTOBHits",     &tree_track_numberOfValidStripTOBHits);
    smalltree->Branch("tree_track_numberOfValidStripTECHits",     &tree_track_numberOfValidStripTECHits);
    smalltree->Branch("tree_track_numberOfValidPixelBarrelHits",  &tree_track_numberOfValidPixelBarrelHits);
    smalltree->Branch("tree_track_numberOfValidPixelEndcapHits",  &tree_track_numberOfValidPixelEndcapHits);
    smalltree->Branch("tree_track_hasValidHitInPixelLayer",       &tree_track_hasValidHitInPixelLayer);
    smalltree->Branch("tree_track_trackerLayersWithMeasurement",  &tree_track_trackerLayersWithMeasurement);
    smalltree->Branch("tree_track_pixelLayersWithMeasurement",    &tree_track_pixelLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTECLayersWithMeasurement", &tree_track_stripTECLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTIBLayersWithMeasurement", &tree_track_stripTIBLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTIDLayersWithMeasurement", &tree_track_stripTIDLayersWithMeasurement);
    smalltree->Branch("tree_track_stripTOBLayersWithMeasurement", &tree_track_stripTOBLayersWithMeasurement);
    
    smalltree->Branch("tree_track_vx",           &tree_track_vx );
    smalltree->Branch("tree_track_vy",           &tree_track_vy );
    smalltree->Branch("tree_track_vz",           &tree_track_vz );
    smalltree->Branch("tree_track_firsthit_X",   &tree_track_firsthit_X);
    smalltree->Branch("tree_track_firsthit_Y",   &tree_track_firsthit_Y);
    smalltree->Branch("tree_track_firsthit_Z",   &tree_track_firsthit_Z);
    smalltree->Branch("tree_track_ntrk10",&tree_track_ntrk10);
    smalltree->Branch("tree_track_ntrk20",&tree_track_ntrk20);
    smalltree->Branch("tree_track_ntrk30",&tree_track_ntrk30);
    
    smalltree->Branch("tree_track_MVAval", &tree_track_MVAval);
    smalltree->Branch("tree_track_Hemi", &tree_track_Hemi);
    smalltree->Branch("tree_track_Hemi_dR", &tree_track_Hemi_dR);
    smalltree->Branch("tree_track_Hemi_mva_NChi2", &tree_track_Hemi_mva_NChi2);
    smalltree->Branch("tree_track_Hemi_LLP", &tree_track_Hemi_LLP);
    
    smalltree->Branch("tree_track_recoVertex_idx", &tree_track_recoVertex_idx);
    smalltree->Branch("tree_track_recoAK4SlimmedJet_idx", &tree_track_recoAK4SlimmedJet_idx);
    
    // info about the simulated track matched to the reco track
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

  std::cout<<"event number :"<<eventNumber<<std::endl;
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
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  
  //----------------paul------------//
  iSetup.get<IdealMagneticFieldRecord>().get(bField);

  edm::ESHandle<EcalChannelStatusRcd> chanstat;//useless atm
  iSetup.get<EcalChannelStatusRcd>().get(chanstat);//useless atm

  //-------------end of Paul---------//

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
  // float  Gen_neu1_eta=-10, Gen_neu1_phi=-10;
  // float  Gen_neu2_eta=-10, Gen_neu2_phi=-10;
//   int  neu[2];
  int  nneu = 0;
  TLorentzVector vneu[2];
  
//   float dRneuneu = 0.;
  
//   for (int k=0; k<2; k++) {
//     neu[k] = -1;
//   }
  /////////////////////////////////////////////////////////
  
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
// 	if ( neu[0] < 0 ) {
// 	  neu[0] = genParticle_idx;
// 	  vneu[0].SetPtEtaPhiM( Gen_pt, Gen_eta, Gen_phi, Gen_m );
// 	  Gen_neu1_eta = Gen_eta;
// 	  Gen_neu1_phi = Gen_phi;
// 	}
// 	else if ( neu[1] < 0 ) {
// 	  neu[1] = genParticle_idx;
// 	  vneu[1].SetPtEtaPhiM( Gen_pt, Gen_eta, Gen_phi, Gen_m );
// 	  Gen_neu2_eta = Gen_eta;
// 	  Gen_neu2_phi = Gen_phi;
// 	}
	nneu++;
      }
      
//       if ( nneu == 2 ) {
// 	dRneuneu = Deltar( Gen_neu1_eta, Gen_neu1_phi, Gen_neu2_eta, Gen_neu2_phi );
//       }
      
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
  ////////////   Muons   ///////////
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
    // tree_muon_dxy.push_back(0.5);
    tree_muon_dxyError.push_back( mu.muonBestTrack()->dxyError());
    tree_muon_dz.push_back(       mu.muonBestTrack()->dz(PV.position()));
    // tree_muon_dz.push_back(0.5);
    tree_muon_dzError.push_back(  mu.muonBestTrack()->dzError());
    tree_muon_charge.push_back(   mu.charge());
    tree_muon_isLoose.push_back(  mu.isLooseMuon());
    // tree_muon_isTight.push_back(  mu.isTightMuon(PV));
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
  
  if ( tree_passesHTFilter ) {
//$$
  
  //////////////////////////////////
  //////////////////////////////////
  //////////   Tracks  /////////////
  //////////////////////////////////
  //////////////////////////////////

  edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder);//Asking for reco colelction of PV..
  vector<reco::TransientTrack> BestTracks;

    // loop on pf candidates
    for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
      const pat::PackedCandidate &pf = (*pfs)[i];
    if ( pf.charge() == 0 || pf.pt() < 1. ) continue;
      
      tree_nTracks++; 
      // BestTracks.push_back(theTransientTrackBuilder->build(pf.bestTrack()));
      tree_track_pt.push_back                (pf.pt());
      tree_track_eta.push_back               (pf.eta());
      tree_track_phi.push_back               (pf.phi());
      tree_track_charge.push_back            (pf.charge());
//       tree_track_nhits.push_back             (pf.numberOfValidHits());
      tree_track_NChi2.push_back             (pf.bestTrack()->normalizedChi2());
//       tree_track_numberOfValidHits.push_back (pf.numberOfValidHits());
      tree_track_dxy.push_back               (pf.dxy((*primaryVertex)[0].position()));
      tree_track_dz.push_back                (pf.dz((*primaryVertex)[0].position()));
//       tree_track_isHighQuality               (pf.trackHighPurity());
//       tree_track_numberOfValidPixelHits.push_back       (pf.hitPattern().numberOfValidPixelHits());
//       tree_track_numberOfValidStripHits.push_back       (pf.hitPattern().numberOfValidStripHits());
//       tree_track_numberOfValidStripTIBHits.push_back    (pf.hitPattern().numberOfValidStripTIBHits());
//       tree_track_numberOfValidStripTIDHits.push_back    (pf.hitPattern().numberOfValidStripTIDHits());
//       tree_track_numberOfValidStripTOBHits.push_back    (pf.hitPattern().numberOfValidStripTOBHits());
//       tree_track_numberOfValidStripTECHits.push_back    (pf.hitPattern().numberOfValidStripTECHits());
//       tree_track_numberOfValidPixelBarrelHits.push_back (pf.hitPattern().numberOfValidPixelBarrelHits());
//       tree_track_numberOfValidPixelEndcapHits.push_back (pf.hitPattern().numberOfValidPixelEndcapHits());
//       tree_track_hasValidHitInPixelLayer.push_back      (pf.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1);
//       tree_track_trackerLayersWithMeasurement.push_back (pf.hitPattern().trackerLayersWithMeasurement());
//       tree_track_pixelLayersWithMeasurement.push_back   (pf.hitPattern().pixelLayersWithMeasurement());
      tree_track_numberOfValidPixelHits.push_back       (pf.numberOfPixelHits());
      tree_track_vx.push_back                 (pf.vx());
      tree_track_vy.push_back                 (pf.vy());
      tree_track_vz.push_back                 (pf.vz());

      //----------------------------paul--------------------------------//
      if (pf.firstHit()!=0)
        {
          const FreeTrajectoryState Freetraj = BestTracks[i].initialFreeState();//SegFault
          GlobalPoint vert (pf.vx(),pf.vy(),pf.vz());
          // const TrajectoryStateOnSurface Freetraj = BestTracks[i].stateOnSurface(vert);//SegFault : same as above
          Cylinder Layer1 = Cylinder(3);//radius
          std::cout<<"MINIAOD TT FReeTRaj=> glob X: "<<Freetraj.position().x()<<"|Y: "<<Freetraj.position().y()<<"|Z"<<Freetraj.position().z()<<std::endl;
          // std::cout<<"MINIAOD TT FReeTRaj=> glob X: "<<Freetraj.globalPosition().x()<<"|Y: "<<Freetraj.globalPosition().y()<<"|Z"<<Freetraj.globalPosition().z()<<std::endl;

          // const MagneticField* Bfield = BestTracks[i].field();
          // const Propagator* InsideTk;//to be initialized?
          // const Propagator* OutsideTk;//to be initialized?
          // SmartPropagator SmartProp(*InsideTk,*OutsideTk,Bfield);
          // TrajectoryStateOnSurface tsos = SmartProp.propagate(Freetraj,Layer1);
          // std::cout<<"MINIAOD TT Propagator=> tsos glob X: "<<tsos.globalPosition().x()<<"|Y: "<<tsos.globalPosition().y()<<"|Z"<<tsos.globalPosition().z()<<std::endl;
          // std::cout<<"MINIAOD TT Propagator=> tsos local X: "<<tsos.localPosition().x()<<"|Y: "<<tsos.localPosition().y()<<"|Z"<<tsos.localPosition().z()<<std::endl;
        }      
        //----------------IMPORTANT-----------------//
        //  The use of propagate will be different depending on the "initialfreestate" obtained above.
        // Disk might get considered as plane, so the propagate method will change
        //+ Cylinders and disk will have to be defined for each possible hitpattern (firsthit)
        //----------------------------------------------------------------------------

      //---------------------------End of paul--------------------------//

    }
  }//HTfilter
//$$
  
  //////////////////////////////////
  
  smalltree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void
FlyingTopAnalyzer::beginJob( )
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
    tree_track_pt.clear();
    tree_track_outerPt.clear();
    tree_track_eta.clear();
    tree_track_phi.clear();
    tree_track_charge.clear();
    tree_track_nhits.clear();
    tree_track_NChi2.clear();
    tree_track_isHighQuality.clear();
    tree_track_isLoose.clear();
    tree_track_isTight.clear();
    
    tree_track_dxy.clear();
    tree_track_dxyError.clear();
    tree_track_dz.clear();
    tree_track_dzError.clear();
    tree_track_numberOfLostHits.clear();
    tree_track_numberOfValidHits.clear();
    tree_track_originalAlgo.clear();
    tree_track_algo.clear();
    tree_track_numberOfValidPixelHits.clear();
    tree_track_numberOfValidStripHits.clear();
    tree_track_numberOfValidStripTIBHits.clear();
    tree_track_numberOfValidStripTIDHits.clear();
    tree_track_numberOfValidStripTOBHits.clear();
    tree_track_numberOfValidStripTECHits.clear();
    tree_track_numberOfValidPixelBarrelHits.clear();
    tree_track_numberOfValidPixelEndcapHits.clear();
    tree_track_hasValidHitInPixelLayer.clear();
    tree_track_trackerLayersWithMeasurement.clear();
    tree_track_pixelLayersWithMeasurement.clear();
    tree_track_stripTECLayersWithMeasurement .clear();
    tree_track_stripTIBLayersWithMeasurement.clear();
    tree_track_stripTIDLayersWithMeasurement.clear();
    tree_track_stripTOBLayersWithMeasurement.clear();
    
    tree_track_vx.clear();
    tree_track_vy.clear();
    tree_track_vz.clear();
    tree_track_firsthit_X.clear();
    tree_track_firsthit_Y.clear();
    tree_track_firsthit_Z.clear();
    tree_track_ntrk10.clear();
    tree_track_ntrk20.clear();
    tree_track_ntrk30.clear();
    
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
