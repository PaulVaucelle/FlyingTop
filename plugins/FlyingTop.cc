// system include files

///////////////////////////////////////////////////////////////////
//---------------------------------------------------------------//
//---------------  Displaced Top quark Analysis Code ------------//
//---------------------------------------------------------------//
///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
// This code aims at reconstructing long-lived decays of new     //
// massive particle decaying into a top and a stop quark         //
// https://arxiv.org/pdf/2212.06678.pdf -> see production of     //
// smuon decaying into muon and neutralino. The Neutralino       //
// decaying to a top and a stop with the RPV coupling            //
// lambda'' 312                                                  //
///////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////
//---------------------------------------------------------------//
// The electron and tau channels could also be analyzed but      //
// since their identification/isolation is way different from    //
// the muons (and also harder), the muon channel is chosen even  //
// though the code is also ready for the electron channel.       //
// We will use the prompt muons to trigger.                      //
//  For a signal event, we aim at reconstructing a maximum of 2  //
// vertices but conveners ask for potentially 3 or more vertices //
// so there could be some tests performed to reconstruct other   //
// vertices.                                                     //
//---------------------------------------------------------------//
///////////////////////////////////////////////////////////////////

//Please contact me at : paul.vaucelle@cern.ch (and before the 01/10/25 dd/mm/yyyy as I finish my Ph.d 
// on this date and I will not be able to answer to your questions after this date)


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
#include "RoccoR.h"
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
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>
//#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "PhysicsTools/Utilities/interface/LumiReweightingStandAlone.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"
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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

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
#include "../interface/Vtx.h"
#include "../interface/Filter.h"
              //----------------Trigger---------------------------//
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include <TEfficiency.h>
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//--PU--//
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
//----Others----//
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
//-----------COnversion---------//

#include "RecoEgamma/EgammaPhotonAlgos/interface/ConversionHitChecker.h"
#include "DataFormats/EgammaTrackReco/interface/ConversionTrack.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
 
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "DataFormats/GeometrySurface/interface/BoundCylinder.h"
#include "DataFormats/GeometrySurface/interface/BoundDisk.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleDiskBounds.h"
 
 //-------------CaloCluster---------------------//
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
// #include "../interface/OwnConversion.cc"

//--------------------Muons---------------------//
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//----------------------Generator/LHE infos--------------------------//

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHECommonBlocks.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHERunInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h"

//-------------------------Top pt reweighting ---------------------//
// https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_dat
// #include "AnalysisDataFormats/TopObjects/interface/TtGenEvent.h"
//-----------------------------------------------------------------//
//
// class declaration
//

// skeleton from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#4_7_MiniAOD_Analysis_Documentati
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
// edm::one::EDAnalyzer<> will be mandatory for CMSSW_X with X>=13

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

    RoccoR rc;
  
    bool isMC_;
    int YEAR_ ;
    string RochString;
    std::string weightFile_;
    std::string weightFileEVTS_;
    std::string weightFileEVTSDY_ ;
    std::string weightFileEVTSTT_;
    std::string weightFileHEMI1_;
    std::string weightFileHEMI1DY_;
    std::string weightFileHEMI1TT_;
    std::string weightFileHEMI2_;
    std::string weightFileHEMI2DY_;
    std::string weightFileHEMI2TT_;
    std::string weightFileVtx_;
    std::string weightFileVtxStep1_;
    std::string mcPileupFile_, dataPileupFile_;
    std::string mcPileupPath_, dataPileupPath_;
    std::string MuonEps1File_ ;
    std::string MuonEps1Path_ ;
    std::string MuonEps2File_ ;
    std::string MuonEps2Path_ ;
    std::string MuonEps3File_ ;
    std::string MuonEps3Path_ ;

    const edm::EDGetTokenT<GenEventInfoProduct>           genEventInfoToken_;
    const edm::EDGetTokenT<LHEEventProduct>               LHEEventProductToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> >       prunedGenToken_;
    edm::EDGetTokenT<edm::View<pat::PackedGenParticle> >  packedGenToken_;
    edm::EDGetTokenT<reco::VertexCollection>              vertexToken_;
    edm::EDGetTokenT<pat::METCollection>                  metToken_;
    edm::EDGetTokenT<pat::JetCollection>                  jetToken_;
    edm::EDGetTokenT<edm::View<reco::GenJet> >            genJetToken_;
    edm::EDGetTokenT<pat::ElectronCollection>             electronToken_;
    edm::EDGetTokenT<pat::MuonCollection>                 muonToken_;
    edm::EDGetTokenT<pat::PackedCandidateCollection>      pcToken_;
    edm::EDGetTokenT<pat::PackedCandidateCollection>      lostpcToken_; //LOST

    reweight::LumiReWeighting* lumiWeights_;

  ///////////////
  // Ntuple info

    TTree *smalltree;
    
    edm::Service<TFileService> fs;
    
    edm::ESHandle<MagneticField> bField;
    
    // edm::ParameterSet kvfPSet;
    //trig
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
    //trig
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> K0Token_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> LambdaToken_;
    // edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> SVToken_;

    edm::EDGetTokenT<vector<PileupSummaryInfo>>     puToken_;
    //edm::EDGetTokenT<vector<PileupSummaryInfo>>   thePUTag;
    edm::EDGetTokenT<reco::ConversionCollection>    PhotonToken_;
    const edm::EDGetTokenT<reco::BeamSpot>          beamSpotToken_;

    // edm::EDGetTokenT<pat::PackedTriggerPrescales> PrescaleToken_;
    edm::EDGetTokenT< double >                      prefweight_token;
    edm::EDGetTokenT<double>                        rho_token_;

    int runNumber, eventNumber, lumiBlock;
    double PUweight;
    int PU_events;// AllPU_events_weight;
    double Prefweight;
    bool tree_only_tigger_filter;
    bool tree_Filter;
    bool tree_FilterSameSign;
    bool tree_Good_PV;
    int tree_muon_GenRecoTriggerMatched;
    int  tree_nTracks, tree_nLostTracks, tree_TRACK_SIZE; 
    int  tree_nFromC = 0, tree_nFromB = 0; 
    int  nEvent;
    
    int   LLP1_mother = 0, LLP2_mother = 0;    
    float LLP1_pt, LLP1_eta, LLP1_phi, LLP2_pt, LLP2_eta, LLP2_phi;
    float LLP1_x, LLP1_y, LLP1_z, LLP2_x, LLP2_y, LLP2_z;
    float LLP1_dist, LLP2_dist;
    int   LLP1_nTrks = 0, LLP2_nTrks = 0;    

  float  Evts_muon1_pt;
  float  Evts_muon2_pt;
  float  Evts_muon1_eta;
  float  Evts_muon2_eta;
  float  Evts_muon1_phi;
  float  Evts_muon2_phi;
  float  Evts_muon12_dR;
  float  Evts_muon12_dPhi;
  float  Evts_muon12_dEta;

 
  //The BDT variables are declared here to reduce computation time
  //EVTS level BDT to select signal events
    double Rho = 0;
    float  mva_Evts_MET_et;
    float  mva_Evts_nTrks;
    float  mva_Evts_muon1_pt;
    float  mva_Evts_muon2_pt;
    float  mva_Evts_jet1_pt;
    float  mva_Evts_jet1_eta;
    float  mva_Evts_jet2_pt;
    float  mva_Evts_jet2_eta;
    float  mva_Evts_muon12_dR;
    float  mva_Evts_muon12_dPhi;
    float  mva_Evts_muon12_dEta;
    float  mva_Evts_jet12_dR;
    float  mva_Evts_jet12_dPhi;
    float  mva_Evts_jet12_dEta;
    float  mva_Evts_muon_jet_dRmin0;
    float  mva_Evts_muon_jet_dRmax0;
    float  mva_Evts_muon_jet_dRmin1;
    float  mva_Evts_muon_jet_dRmax1;
    float  mva_Evts_Hemi1_njet_nomu;
    float  mva_Evts_Hemi2_njet_nomu;
    float  mva_Evts_Hemi1_pt;
    float  mva_Evts_Hemi2_pt;
    float  mva_Evts_Hemi1_eta;
    float  mva_Evts_Hemi2_eta;
    float  mva_Evts_Hemi1_phi;
    float  mva_Evts_Hemi2_phi;
    float  mva_Evts_Hemi1_nTrks;
    float  mva_Evts_Hemi2_nTrks;
    float  mva_Evts_Hemi1_Mass;
    float  mva_Evts_Hemi2_Mass;
    float  mva_Evts_nVtx;
    float  mva_HT;
    float  mva_Evts_ST;
    float  mva_Evts_njets;
    float  mva_Evts_nmuon;
    float  mva_Evts_MediumAxes;
    float  mva_Evts_LooseAxes;
    float  mva_Evts_TightAxes;
    float  mva_Evts_Mmumu;
    float  mva_Evts_all_muon;


    TMVA::Reader *readerEvts = new TMVA::Reader("!Color:Silent");
    // TMVA::Reader *readerHemi1 = new TMVA::Reader("!Color:Silent");
    // TMVA::Reader *readerHemi2 = new TMVA::Reader("!Color:Silent");
  //-------------------------------//
    float pt, eta, NChi, nhits, ntrk10, drSig, isinjet, phi;
    float firsthit_X; /*!*/
    float firsthit_Y; /*!*/
    float firsthit_Z; /*!*/
    float dxy; /*!*/
    float dxyError; /*!*/
    float dz; /*!*/
    float dzSig;
    float ntrk20; /*!*/
    float ntrk30; /*!*/
    float ntrk40; /*!*/
    float dR; /*!*/
    float dRmax; /*!*/
//$$$$$$
//     float dzTopu;
//     float dzSigTopu;
//$$$$$$
    float TibHit;
    float TobHit;
    float PixBarHit;
    float TecHit;
    float isLost;
//$$$$$$
//     float ntrk10rel;
//     float ntrk20rel;
//     float ntrk30rel;
//     float ntrk10reltot;
//     float ntrk20reltot;
//     float ntrk30reltot;
//     float ntrk40reltot;
//$$$$$$
    //Track level BDT for displaced track selection
    TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );

    //Vtx level BDT to select signal vertices
    float mva_V_nTrks = 0 ;
    float mva_V_chi = 0;
    float mva_V_step = 0;
    float mva_V_r= 0;
    float mva_V_z = 0;
    float mva_V_MTW = 0;
    float mva_V_Mass =  0;
    float mva_H_Mass = 0;
    float mva_V_dist = 0;
    float mva_V_ntrk10 = 0;
    float mva_V_ntrk20 = 0 ;
    float mva_V_MeanDCA = 0;
    TMVA::Reader *readerVtx = new TMVA::Reader( "!Color:Silent" );
    TMVA::Reader *readerVtxStep1 = new TMVA::Reader( "!Color:Silent" );

    int index[10000];
    double MVAval[10000];

    int index_muon[10000];
    double muon_pt[10000];

    int index_el[10000];
    double el_pt[10000];

    //  ---------------------------------------------------------------- //
    //  ------- Booleans to activate/desactivate part of the code ------ //
    //  ---------------------------------------------------------------- //

    bool showlog            = false;
    
    bool MuonChannel        = true;
    bool ElChannel          = false;
    bool EMuChannel         = false;

    bool AllowDiLeptonSameSign = true;

    bool NewCovMat          = true;//Keep True : Allow for Covariance Matrix correction due to the MiniAOD dataformat apporixmation
          //Vetos to find vertices from different secondary interactions
    bool DetailedMap        = true;// Detailed map of the CMS tracker to apply a veto on the tracks of the vertices that belong to this map
        // Vetos applied on tracks of the vertices belonging to V0Candidates, Photon conversions and Secindary Interactions
    bool ActivateV0Veto     = true;
    bool ActivateYcVeto     = true;// This Veto is not doing anything on RunIISummer20UL18 TTbar and on MC signal
    bool ActivateSecIntVeto = true;
        // Activate steps of the vertexing workflow (better to keep everything true for the development, since we may want to keeep the tight WP => true false true false false)
    bool ActivateStep1      = true; //TIghtwp, standard AVF
    bool IterAVF            = true; // Activate IAVF step of the vertexing 
    bool ActivateStep2      = true; // Tight WP IAVF => Works like this (ActivateStep2 || IterAVF)
    bool ActivateStep3      = true; // LooseWP STep3 => standard AVF
    bool ActivateStep4      = true; // LooseWP STep4  IAVF => (ActivateStep4 || IterAVF)
    bool ActivateMerging    = true;
    bool ActivateONLYAVF    = false; 

    // -- B Tagging related information
    // WorkingPoints : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    // Updated Page by BTV : https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
    bool ActivateBtag       = true ;
    bool ActivateBtagLog    = false;
    float LooseWP  = 0.0490;
    float MediumWP = 0.2783;
    float TightWP  = 0.7100;

    //-Track Selection-------/
    bool IncludeLostTrack= true;// Keep as true // For 2023 and beyong Lost track will have their hitpattern implemented following the request we made
    bool RequestHighPurity = false ; // does not rlly matter since the TRK BDT selects mainly HighPuirty tracks (99%)
    //-END of Track Selection-------/

    //------Vtx Selection ---//
    bool RemoveLostTrackFromVtxSelec = true;// This should be removed for 2023 and beyong since we will have the hitpattern 
    //  ---------------------------------------------------------------- //
 
    //  ---------------------------------------------------------------- //
    //  -------------------- track preselection cuts ------------------- //
    //  ---------------------------------------------------------------- //
    float pt_Cut = 1;    // default 1. 
    float NChi2_Cut = 5; // default 5. 
    float drSig_Cut = 5; // default 5. 
   //  ---------------------------------------------------------------- //

  //    std::vector<float> tree_LHE_Weights;
  // float tree_MCEvt_weight;

    float tree_Evts_MVAval;
    float tree_Evts_MVAvalDY;
    float tree_Evts_MVAvalTT;

    float tree_Hemi1_MVAval ;
    float tree_Hemi1_MVAvalDY ;
    float tree_Hemi1_MVAvalTT ;
    float tree_Hemi2_MVAval ;
    float tree_Hemi2_MVAvalDY ;
    float tree_Hemi2_MVAvalTT ;

    //---------------  Gen event wt------------------------//

    TH1D *hEvents;
    TH1D *hEvents_with_gen_wt;
    double tree_only_gen_wt;
    double tree_event_weight;
    double tree_genTop_Weight;


    //--------------------------------
    // primary vertex infos -------
    //--------------------------------

    float tree_bs_PosX ;
    float tree_bs_PosY ;
    float tree_bs_PosZ ;
  
    int   tree_nPV;
    float tree_PV_x;
    float tree_PV_y;
    float tree_PV_z;
    float tree_PV_ez;
    float tree_PV_NChi2;
    int   tree_PV_ndf;
    float tree_PV_rho;

    //--------------------------------
    // muons infos -------
    //--------------------------------
    int tree_all_nmu; // count all muons
    int tree_nmu;     // count prompt muons
    float tree_LT;
    float tree_Mmumu;
    float tree_MmumuSameSign;

    std::vector<bool>  tree_muon_isPrompt; // prompt candidate 
    std::vector<float> tree_muon_pt;
    std::vector<float> tree_muon_SF;
    std::vector<float> tree_muon_eta;
    std::vector<float> tree_muon_phi;
    std::vector<float> tree_muon_x;
    std::vector<float> tree_muon_y;
    std::vector<float> tree_muon_z;
    std::vector<float> tree_muon_dxy;
    std::vector<float> tree_muon_dxyError;
    std::vector<float> tree_muon_dz;
    std::vector<float> tree_muon_dzError;
    std::vector< int > tree_muon_charge;
    std::vector<bool>  tree_muon_isLoose;
    std::vector<bool>  tree_muon_isMedium;
    std::vector<bool>  tree_muon_isTight;
    std::vector<bool>  tree_muon_isGlobal;
    std::vector<float> tree_muon_isoR3;
    std::vector<bool>  tree_muon_trigger_dimu;  
    std::vector<bool>  tree_muon_trigger_isomu;
    std::vector<bool>  tree_muon_PFIsoVeryLoose;
    std::vector<bool>  tree_muon_PFIsoLoose;
    std::vector<bool>  tree_muon_PFIsoMedium;
    std::vector<bool>  tree_muon_PFIsoTight;
    std::vector<bool>  tree_muon_MiniIsoLoose;
    std::vector<bool>  tree_muon_MiniIsoMedium;
    std::vector<bool>  tree_muon_MiniIsoTight;
    std::vector<bool>  tree_muon_TkIsoLoose;
    std::vector<bool>  tree_muon_TkIsoTight;
    std::vector<float> tree_muon_trkLayers;
    std::vector<float> tree_muon_miniIso;
    std::vector<float> tree_muon_correction;
    std::vector<int>   tree_muon_gen; // generated parent pdgid from reco muon

  std::vector<float> tree_reco_muon_leadingpt;
  std::vector<float> tree_reco_electron_leadingpt2;
  std::vector<float> tree_reco_muon_leadingeta;
  std::vector<float> tree_reco_electron_leadingeta2;
  std::vector<float> tree_reco_muon_leadingphi;
  std::vector<float> tree_reco_electron_leadingphi2;
  
  std::vector<float> tree_trig_muon_leadingpt;
  std::vector<float> tree_trig_electron_leadingpt2;
  std::vector<float> tree_trig_muon_leadingeta;
  std::vector<float> tree_trig_electron_leadingeta2;
  std::vector<float> tree_trig_muon_leadingphi;
  std::vector<float> tree_trig_electron_leadingphi2;

  std::vector<float> tree_reco_lepton_leadingpt;
  std::vector<float> tree_reco_lepton_leadingpt2;
  std::vector<float> tree_reco_lepton_leadingeta;
  std::vector<float> tree_reco_lepton_leadingeta2;
  std::vector<float> tree_reco_lepton_leadingphi;
  std::vector<float> tree_reco_lepton_leadingphi2;
    
  std::vector<float> tree_trig_lepton_leadingpt;
  std::vector<float> tree_trig_lepton_leadingpt2;
  std::vector<float> tree_trig_lepton_leadingeta;
  std::vector<float> tree_trig_lepton_leadingeta2;
  std::vector<float> tree_trig_lepton_leadingphi;
  std::vector<float> tree_trig_lepton_leadingphi2;

    std::vector<float> tree_lepton_leadingpt;
    std::vector<float> tree_lepton_leadingpt2;
    std::vector<float> tree_lepton_leadingeta;
    std::vector<float> tree_lepton_leadingeta2;
    std::vector<float> tree_lepton_leadingphi;
    std::vector<float> tree_lepton_leadingphi2;

    std::vector<float> tree_lepton_lepton_dR;
    std::vector<float> tree_lepton_lepton_dPhi;
    std::vector<float> tree_lepton_lepton_dEta;
    
    std::vector<float> tree_ll_pt;
    std::vector<float> tree_ll_eta;
    std::vector<float> tree_ll_phi;
    std::vector<float> tree_ll_px;
    std::vector<float> tree_ll_py;
    std::vector<float> tree_ll_pz;
    std::vector<float> tree_ll_energy;
    std::vector<float> tree_ll_mass;

    //--------------------------------
    // electrons infos -------
    //--------------------------------

    int                tree_all_nel;         // count all electrons
    int                tree_electron_nEle;   // count prompt electrons
    std::vector<bool>  tree_electron_isPrompt;
    std::vector<float> tree_electron_pt;
    std::vector<float> tree_electron_eta;
    std::vector<float> tree_electron_phi;
    std::vector<float> tree_electron_x;
    std::vector<float> tree_electron_y;
    std::vector<float> tree_electron_z;
    std::vector<float> tree_electron_energy;
    std::vector< int > tree_electron_charge;
    std::vector<float> tree_electron_et;
    std::vector<float> tree_electron_ecal_trk_postcorr;
    std::vector<float> tree_electron_isoR4;
    std::vector<bool>  tree_electron_IsLoose;
    std::vector<bool>  tree_electron_IsMedium;
    std::vector<bool>  tree_electron_IsTight;
    std::vector<float> tree_electron_dxy;
    std::vector<float> tree_electron_dz;
    std::vector<int>   tree_electron_gen; // generated parent pdgid from reco electron

    //--------------------------------
    // met infos -------
    //--------------------------------

    float tree_PFMet_et;
    float tree_PFMet_phi;
    float tree_PFMet_sig;
    float tree_PFMet_pt;
    
    //--------------------------------
    // jet infos -------
    //--------------------------------
    
    int tree_njet; 
    int tree_njetNOmu; // only for counting jets without prompt muon inside !
    std::vector<float> tree_jet_pt;
    std::vector<float> tree_jet_eta;
    std::vector<float> tree_jet_phi;
    std::vector<float> tree_jet_px;
    std::vector<float> tree_jet_py;
    std::vector<float> tree_jet_pz;
    std::vector<float> tree_jet_E;
    std::vector<bool>  tree_jet_tightid_LepVeto;
    std::vector<bool>  tree_jet_tightid;
    std::vector<bool>  tree_jet_TightJetIDLepVeto;
    std::vector<bool>  tree_jet_TightJetID;
    std::vector<float> tree_jet_HadronFlavour;
    std::vector<float> tree_jet_pileupID;
    std::vector<float> tree_jet_btag_DeepCSV;
    std::vector<float> tree_jet_btag_DeepJet;
    std::vector<float> tree_jet_leadingpt;
    std::vector<float> tree_jet_leadingpt2;
    std::vector<float> tree_jet_leadingeta;
    std::vector<float> tree_jet_leadingeta2;
    std::vector<float> tree_jet_leadingMuon_dR;
    std::vector<float> tree_jet_leadingMuon2_dR;
    std::vector<float> tree_jet_jet_dR;
    std::vector<float> tree_jet_jet_dPhi;
    std::vector<float> tree_jet_jet_dEta;
    std::vector<float> tree_muon_jet_dRmin;
    std::vector<float> tree_muon_jet_dRmax;
    std::vector<float> tree_elemu_jet_dRmin;
    std::vector<float> tree_elemu_jet_dRmax;
    std::vector<float> tree_ele_jet_dRmin;
    std::vector<float> tree_ele_jet_dRmax;
    float tree_HT;
    
    //--------------------------------
    // ------ V0 Candidates  ---------
    //--------------------------------

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

    int tree_nLambda;
    std::vector<float>     tree_L0_x;
    std::vector<float>     tree_L0_y;
    std::vector<float>     tree_L0_z;
    std::vector<float>     tree_L0_r;
    std::vector<float>     tree_L0_NChi2;
    std::vector<float>     tree_L0_ndf;
    std::vector<float>     tree_L0_mass;
    std::vector<float>     tree_L0_pt;
    std::vector<float>     tree_L0_eta;
    std::vector<float>     tree_L0_phi;
    std::vector<unsigned int> tree_L0_nDaughters;

    // reconstructed V0
    int tree_nV0_reco;
    std::vector<float>     tree_V0_reco_x;
    std::vector<float>     tree_V0_reco_y;
    std::vector<float>     tree_V0_reco_z;
    std::vector<float>     tree_V0_reco_r;
    std::vector<float>     tree_V0_reco_drSig;
    std::vector<float>     tree_V0_reco_dzSig;
    std::vector<float>     tree_V0_reco_angleXY;
    std::vector<float>     tree_V0_reco_angleZ;
    std::vector<float>     tree_V0_reco_NChi2;
    std::vector<float>     tree_V0_reco_ndf;
    std::vector<float>     tree_V0_reco_mass;
    std::vector<float>     tree_V0_reco_pt;
    std::vector<float>     tree_V0_reco_eta;
    std::vector<float>     tree_V0_reco_phi;
    std::vector<int>       tree_V0_reco_source;
    std::vector<bool>      tree_V0_reco_badTkHit;
    std::vector<float>     tree_V0_reco_dca;

    //--------------------------------
    // ------ Secondary Interactions -
    //--------------------------------

    int tree_nSecInt;
    std::vector<float>     tree_SecInt_x;
    std::vector<float>     tree_SecInt_y;
    std::vector<float>     tree_SecInt_z;
    std::vector<float>     tree_SecInt_r;
    std::vector<float>     tree_SecInt_d;
    std::vector<float>     tree_SecInt_drSig;
    std::vector<float>     tree_SecInt_dzSig;
    std::vector<float>     tree_SecInt_angleXY;
    std::vector<float>     tree_SecInt_angleZ;
    std::vector<float>     tree_SecInt_NChi2;
    std::vector<float>     tree_SecInt_ndf;
    std::vector<float>     tree_SecInt_mass;
    std::vector<float>     tree_SecInt_pt;
    std::vector<float>     tree_SecInt_eta;
    std::vector<float>     tree_SecInt_phi;
    std::vector<int>       tree_SecInt_charge;
    std::vector<bool>      tree_SecInt_badTkHit;
    std::vector<float>     tree_SecInt_dca;
    std::vector<bool>      tree_SecInt_selec;
    std::vector<int>       tree_SecInt_layer;
    std::vector<int>       tree_SecInt_LLP;
    std::vector<float>     tree_SecInt_LLP_dr;
    std::vector<float>     tree_SecInt_LLP_dz;
    std::vector<float>     tree_SecInt_LLP_dd;
    std::vector<unsigned int>     tree_SecInt_tk1;
    std::vector<unsigned int>     tree_SecInt_tk2;

    //---------------------------------------------------------
    // ------ Photon Conversions => Fomm CMSSW collection -----
    //---------------------------------------------------------

    int tree_nYConv;
    std::vector<float>     tree_Yc_x; 
    std::vector<float>     tree_Yc_y;
    std::vector<float>     tree_Yc_z;
    std::vector<float>     tree_Yc_r;
    std::vector<float>     tree_Yc_dr0;
    std::vector<float>     tree_Yc_dr1;
    std::vector<float>     tree_Yc_dz0;
    std::vector<float>     tree_Yc_dz1;
    std::vector<float>     tree_Yc_costheta;
    std::vector<int>       tree_Yc_layer;
    std::vector<float>     tree_Yc_NChi2;
    std::vector<float>     tree_Yc_ndf;
    std::vector<unsigned int> tree_Yc_nDaughters;
    std::vector<float>     tree_Yc_pt;
    std::vector<float>     tree_Yc_eta;
    std::vector<float>     tree_Yc_phi;
    std::vector<float>     tree_Yc_mass;

    int tree_Yc_ntracks;
    std::vector<int>       tree_Yc_tracks_index;
    std::vector<int>       tree_Yc_tracks_charge;
    std::vector<float>     tree_Yc_tracks_pt;
    std::vector<float>     tree_Yc_tracks_eta;
    std::vector<float>     tree_Yc_tracks_phi;
    std::vector<float>     tree_Yc_tracks_phi0;

    //-----------------------
    // track info
    //-----------------------

//     std::vector<bool>     tree_passesTrkPtr;
    std::vector<unsigned int> tree_track_ipc;
    std::vector<bool>     tree_track_lost;
    std::vector<float>    tree_track_px;
    std::vector<float>    tree_track_py;
    std::vector<float>    tree_track_pz;
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
    std::vector<float>    tree_track_dzSig;
//$$$$$$
//     std::vector<float>    tree_track_dzTOpu;  // with respect to clostest PU
//     std::vector<float>    tree_track_dzSigTOpu;
//$$$$$$
//     http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_1_3/doc/html/d8/df2/classreco_1_1TrackBase.html#aca7611bd1a33d535cefc72b6e497ece8

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
    std::vector<float>    tree_track_btag;
    std::vector<float>    tree_track_energy;

    std::vector< int >    tree_track_Hemi;
    std::vector< float >  tree_track_Hemi_dR; // dRmin
    std::vector< float >  tree_track_Hemi_dRmax;
    std::vector< float >  tree_track_Hemi_mva_NChi2;
    std::vector< bool >   tree_track_Hemi_ping;
    std::vector< float >  tree_track_Hemi_dFirstVtx;
    std::vector< int >    tree_track_Hemi_LLP;
//$$$$$$
//     std::vector< float >  tree_track_Hemi_d0;
//     std::vector< float >  tree_track_Hemi_d0Sig;
//     std::vector< float >  tree_track_HemiOp_d0;
//     std::vector< float >  tree_track_HemiOp_d0Sig;
//$$$$$$
    std::vector< int >    tree_track_Hemi_isjet;

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
   
//$$$$
    std::vector<bool>     tree_V0_track_isFromV0;
    std::vector<bool>     tree_V0_track_isFromSI;
    std::vector<bool>     tree_V0_track_lost;
    std::vector<float>    tree_V0_track_pt;
    std::vector<float>    tree_V0_track_eta;
    std::vector<float>    tree_V0_track_phi;
    std::vector<int>      tree_V0_track_charge;
    std::vector<float>    tree_V0_track_NChi2;
    std::vector<float>    tree_V0_track_dxy;
    std::vector<float>    tree_V0_track_drSig;
    std::vector<float>    tree_V0_track_dz;
    std::vector<float>    tree_V0_track_dzSig;
    std::vector<int>      tree_V0_track_nHit;
    std::vector<int>      tree_V0_track_nHitPixel;
    std::vector< int >    tree_V0_track_firstHit;
    std::vector< float >  tree_V0_track_firstHit_x;
    std::vector< float >  tree_V0_track_firstHit_y;
    std::vector< float >  tree_V0_track_firstHit_z;
    std::vector< int >    tree_V0_track_iJet;
    std::vector<float>    tree_V0_track_ntrk10;
    std::vector<float>    tree_V0_track_ntrk20;
    std::vector<float>    tree_V0_track_ntrk30;
    std::vector<float>    tree_V0_track_ntrk40;
    std::vector< int >    tree_V0_track_Hemi;
    std::vector< float >  tree_V0_track_Hemi_dR;
    std::vector< float >  tree_V0_track_Hemi_dRmax;
//$$$$
   
    //--------------------------------
    // gen infos -------
    //--------------------------------

    float tree_GenPVx;
    float tree_GenPVy;
    float tree_GenPVz;
 
    int tree_smu_mass = 0;
    int tree_neu_mass = 0;
    int tree_neu_ctau = 0;
    
    std::vector< float > tree_genParticle_pt;
    std::vector< float > tree_genParticle_eta;
    std::vector< float > tree_genParticle_phi;
    std::vector< float > tree_genParticle_charge;
    std::vector< int >   tree_genParticle_pdgId;
    std::vector< float > tree_genParticle_mass;
    std::vector< float > tree_genParticle_x;
    std::vector< float > tree_genParticle_y;
    std::vector< float > tree_genParticle_z;
    std::vector< float > tree_genParticle_px;
    std::vector< float > tree_genParticle_py;
    std::vector< float > tree_genParticle_pz;
    std::vector< float > tree_genParticle_energy;
    std::vector< bool >  tree_genParticle_isPromptFinalState;
    std::vector< int >   tree_genParticle_statusCode;
    std::vector< int >   tree_genParticle_mother_pdgId;
    std::vector< int >   tree_genParticle_LLP;
    std::vector< float > tree_genParticle_ct;
    std::vector< float > tree_genParticle_ct0;

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
    std::vector< float > tree_genAxis_dPhineuneu;
    std::vector< float > tree_genAxis_dEtaneuneu;
    std::vector< float > tree_GenAxes_Mass;
    std::vector< float > tree_GenAxes_CombinedHemiLeptonMass;
    std::vector< float > tree_GenAxis_Neu_dRmin;
    std::vector< float > tree_GenAxis_Neu_dRmax;
    std::vector< float > tree_GenAxis_RecoAxis_dRmin;//something to look at in the future
    std::vector< float > tree_GenAxis_RecoAxis_dRmax;

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
    std::vector< float > tree_genFromB_dd;
    std::vector< float > tree_genFromB_dr;
    std::vector< float > tree_genFromB_dz;

    //--------------------------------
    // gen jet infos -------
    //--------------------------------

    std::vector<float> tree_genJet_pt;
    std::vector<float> tree_genJet_eta;
    std::vector<float> tree_genJet_phi;
    std::vector<float> tree_genJet_mass;
    std::vector<float> tree_genJet_energy;
    
    //--------------------------------
    // lhe event infos -------
    //--------------------------------
    
    //=> pdf
    
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
    std::vector< float > tree_LLP_r;
    std::vector< float > tree_LLP_dist;
    std::vector< int >   tree_LLP_nTrks;
    std::vector< int >   tree_LLP_Vtx_nTrks;
    std::vector< float > tree_LLP_Vtx_NChi2;
    std::vector< float > tree_LLP_Vtx_dx;
    std::vector< float > tree_LLP_Vtx_dy;
    std::vector< float > tree_LLP_Vtx_dz;
    std::vector< float > tree_LLP_Vtx_dist;
    std::vector< float > tree_LLP_Vtx_dd;
    std::vector< float > tree_LLP12_dR;
    std::vector< float > tree_LLP12_deta;
    std::vector< float > tree_LLP12_dphi;
    std::vector< float > tree_LLP_Mass;

    //-----------------------
    // - Vertices information
    //-----------------------

    std::vector< int >   tree_Hemi;
    std::vector< int >   tree_Hemi_njet;
    std::vector< int >   tree_Hemi_njet_nomu;
    std::vector< float > tree_Hemi_pt;
    std::vector< float > tree_Hemi_eta;
    std::vector< float > tree_Hemi_phi;
    std::vector< int >   tree_Hemi_nTrks;
    std::vector< int >   tree_Hemi_nTrks_sig;
    std::vector< int >   tree_Hemi_nTrks_bad;
    std::vector< float > tree_Hemi_mass;
    std::vector< float > tree_HemiMu_mass;
    std::vector< float > tree_HemiMu_pt;
    std::vector< float > tree_HemiMu_dR;
    std::vector< float > tree_HemiMuOp_mass;
    std::vector< float > tree_HemiMuOp_pt;
    std::vector< float > tree_HemiMuOp_dR;
    std::vector< int >   tree_Hemi_LooseBTag_axes;
    std::vector< int >   tree_Hemi_MediumBTag_axes;
    std::vector< int >   tree_Hemi_TightBTag_axes;
    std::vector< float > tree_Hemi_dR12;

    std::vector< int >   tree_Hemi_LLP;
    std::vector< float > tree_Hemi_LLP_pt;
    std::vector< float > tree_Hemi_LLP_eta;
    std::vector< float > tree_Hemi_LLP_phi;
    std::vector< float > tree_Hemi_LLP_dist;
    std::vector< float > tree_Hemi_LLP_x;
    std::vector< float > tree_Hemi_LLP_y;
    std::vector< float > tree_Hemi_LLP_z;
    std::vector< float > tree_Hemi_LLP_dR;
    std::vector< int >   tree_Hemi_LLP_mother;
    std::vector< float > tree_Hemi_LLP_Vtx_dx;
    std::vector< float > tree_Hemi_LLP_Vtx_dy;
    std::vector< float > tree_Hemi_LLP_Vtx_dz;
    std::vector< float > tree_Hemi_LLP_Vtx_dr;
        std::vector< float > tree_Hemi_LLP_muOK_dR;
    std::vector< float > tree_Hemi_LLP_muOK_pt;
    std::vector< float > tree_Hemi_LLP_muOK_mass;
    std::vector< float > tree_Hemi_LLP_muNO_dR;
    std::vector< float > tree_Hemi_LLP_muNO_pt;
    std::vector< float > tree_Hemi_LLP_muNO_mass;
    std::vector< float > tree_Hemi_LLP_dR12;
    std::vector< bool >  tree_Hemi_LLP_ping; 
    std::vector< int >   tree_event_LLP_ping;
    
    std::vector< int >   tree_Hemi_Vtx_step;
    std::vector< bool >  tree_Hemi_Vtx_isTight;
    std::vector< float > tree_Hemi_Vtx_NChi2;
    std::vector< int >   tree_Hemi_Vtx_nTrks;
    std::vector< int >   tree_Hemi_Vtx_nTrks_sig;
    std::vector< int >   tree_Hemi_Vtx_nTrks_bad;
    std::vector< float > tree_Hemi_Vtx_x;
    std::vector< float > tree_Hemi_Vtx_y;
    std::vector< float > tree_Hemi_Vtx_z;
    std::vector< float > tree_Hemi_Vtx_r;
    std::vector< float > tree_Hemi_Vtx_dR;
    std::vector< float > tree_Hemi_Vtx_xError;
    std::vector< float > tree_Hemi_Vtx_yError;
    std::vector< float > tree_Hemi_Vtx_zError;
    std::vector< float > tree_Hemi_Vtx_BTag;
    std::vector< float > tree_Hemi_Vtx_trackWeight;
    std::vector< float > tree_Hemi_Vtx_SumtrackWeight;//Vertx selection variable for the BDT
    std::vector< float > tree_Hemi_Vtx_Mass;
    std::vector< float > tree_Hemi_Vtx_track_MeanDCA_d;//Veertex selection BDT
    std::vector< float > tree_Hemi_Vtx_dist;
    std::vector< int >   tree_Hemi_Vtx_ntrk10;//Vertex selection variables
    std::vector< int >   tree_Hemi_Vtx_ntrk20;
    std::vector< int >   tree_event_nVtx;
    std::vector< float > tree_event_Vtx_Vtx_dr;
    std::vector< float > tree_event_Vtx_Vtx_dz;
    std::vector< float > tree_event_Vtx_Vtx_dd;
    std::vector< float > tree_event_Vtx_Vtx_reldd;
    std::vector< float > tree_event_Vtx_Vtx_dR;
    std::vector< int >   tree_event_Vtx_Vtx_step;

    std::vector< float > tree_Hemi_SecLLP;
    std::vector< float > tree_Hemi_LLP_SecVtx_dx;
    std::vector< float > tree_Hemi_LLP_SecVtx_dy;
    std::vector< float > tree_Hemi_LLP_SecVtx_dz;
    std::vector< float > tree_Hemi_LLP_SecVtx_dr;
    std::vector< bool >  tree_Hemi_SecLLP_ping;
    std::vector< int >   tree_event_SecLLP_ping;

    // std::vector< float > tree_Hemi_SecVtx_track_DCA_x;
    // std::vector< float > tree_Hemi_SecVtx_track_DCA_y;
    // std::vector< float > tree_Hemi_SecVtx_track_DCA_z;
    // std::vector< float > tree_Hemi_SecVtx_track_DCA_r;
    // std::vector< float > tree_Hemi_SecVtx_track_DCA_d;
    std::vector< int >   tree_Hemi_SecVtx;      // Hemi (1 or 2) if merging
    std::vector< int >   tree_Hemi_SecVtx_step; // vertex step for this Hemi if merging
    std::vector< float > tree_Hemi_SecVtx_x;
    std::vector< float > tree_Hemi_SecVtx_y;
    std::vector< float > tree_Hemi_SecVtx_z;
    std::vector< float > tree_Hemi_SecVtx_r;
    std::vector< float > tree_Hemi_SecVtx_dR;
    std::vector< float > tree_Hemi_SecVtx_nTrks;
    std::vector< float > tree_Hemi_SecVtx_NChi2;
    std::vector< float > tree_Hemi_SecVtx_dist;
    std::vector< float > tree_Hemi_SecVtx_track_MeanDCA_d;
    std::vector< float > tree_Hemi_SecVtx_SumtrackWeight;
    std::vector< float > tree_Hemi_SecVtx_trackWeight;
    std::vector< float > tree_Hemi_SecVtx_Mass;
    
    std::vector< float > tree_event_MergedVtx_Vtx_dr;
    std::vector< float > tree_event_MergedVtx_Vtx_dz;
    std::vector< float > tree_event_MergedVtx_Vtx_dd;
    std::vector< float > tree_event_MergedVtx_Vtx_reldd;
    std::vector< float > tree_event_MergedVtx_Vtx_dR;
    std::vector< int >   tree_event_MergedVtx_Vtx_step;
    
    std::vector< float > tree_Hemi_Vtx_BDT_nTrks;
    std::vector< float > tree_Hemi_Vtx_BDT_NChi2;
    std::vector< float > tree_Hemi_Vtx_BDT_step;
    std::vector< float > tree_Hemi_Vtx_BDT_STW;
    std::vector< float > tree_Hemi_Vtx_BDT_Mass;
    std::vector< float > tree_Hemi_Vtx_BDT_HMass;
    std::vector< float > tree_Hemi_Vtx_BDT_ntrk10;
    std::vector< float > tree_Hemi_Vtx_BDT_ntrk20;
    std::vector< float > tree_Hemi_Vtx_BDT_MeanDCA;
    std::vector< float > tree_Hemi_Vtx_MVAval_Loose;
    std::vector< float > tree_Hemi_Vtx_MVAval_Tight;//TIght WP

    // Used triggers
    //------- Trigger IsoMu -------//
    bool HLT_IsoMu24_v;     // USED in 2016 and 2018
    bool HLT_IsoMu27_v;     // USED in 2017
    bool HLT_IsoTkMu24_v;

    // ---------------- Trigger MuMu + MuEl -------------

  bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v;                   // USED in 2016-2018                                                                                          
  bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v;
  bool HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v;
  bool HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v;
  bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v;
  bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v;        // USED in 2016-2018                                                                                          
  bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v;  // USED in 2016-2018                                                                                           
  bool HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v;     // USED in 2016-2018                                                                                          
  bool HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;    // USED in 2016-2018                                                                                          
  bool HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v;
  bool HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v; // USED in 2016-2018     


    // ---------------- Trigger Electron-------------
  bool HLT_Ele27_WPTight_Gsf_v;                     // USED in 2016                                                                                                         
  bool HLT_Ele32_WPTight_Gsf_v;                     // USED in 2017-2018                                                                                                    
  bool HLT_Ele35_WPTight_Gsf_v;

  bool HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v; // USED in 2016-2018                                                                                                    
  bool HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v;    // USED in 2016-2018                                                                                                    


    // ---------------- Trigger PFMET -------------
    bool HLT_PFMET120_PFMHT120_IDTight_v;   // USED
    bool HLT_PFMET120_PFMHT120_IDTight_PFHT60_v;   // USED
    bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v;   // USED
    bool HLT_PFMETTypeOne120_PFMHT120_IDTight_PFHT60_v;
    bool HLT_PFMETTypeOne120_PFMHT120_IDTight_v;
    bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v;   // USED
    bool HLT_PFMET250_HBHECleaned_v;   // USED
    bool HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v;   // USED


};


//
// constants, enums and typedefs
//
typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

//
// static data member definitions
//

//
// constructors and destructor
//
FlyingTopAnalyzer::FlyingTopAnalyzer(const edm::ParameterSet& iConfig):
    isMC_(iConfig.getParameter<bool>("isMC")),
    YEAR_ (iConfig.getParameter<int>("YEAR")),
    RochString (iConfig.getParameter<std::string>("RochString")),
    weightFile_( iConfig.getUntrackedParameter<std::string>("weightFileMVA") ),
    weightFileEVTS_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_EVTS")),
    weightFileEVTSDY_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_EVTSDY")),
    weightFileEVTSTT_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_EVTSTT")),
    // weightFileHEMI1_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_HEMI1")),
    // weightFileHEMI1DY_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_HEMI1DY")),
    // weightFileHEMI1TT_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_HEMI1TT")),
    // weightFileHEMI2_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_HEMI2")),
    // weightFileHEMI2DY_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_HEMI2DY")),
    // weightFileHEMI2TT_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_HEMI2TT")),
    weightFileVtx_( iConfig.getUntrackedParameter<std::string>("weightFileMVA_VTX") ),
    weightFileVtxStep1_( iConfig.getUntrackedParameter<std::string>("weightFileMVA_VTX_step1") ),
    mcPileupFile_	( iConfig.getParameter<std::string>( "mcpufile" ) ),
    dataPileupFile_	( iConfig.getParameter<std::string>( "datapufile" ) ),
    mcPileupPath_	( iConfig.getParameter<std::string>( "mcpupath" ) ),
    dataPileupPath_	( iConfig.getParameter<std::string>( "datapupath" ) ),
    // MuonEps1File_ (iConfig.getParameter<std::string>( "muoneps1file" ) ),
    // MuonEps1Path_ (iConfig.getParameter<std::string>( "muoneps1path" ) ),
    // MuonEps2File_ (iConfig.getParameter<std::string>( "muoneps2file" ) ),
    // MuonEps2Path_ (iConfig.getParameter<std::string>( "muoneps2path" ) ),
    // MuonEps3File_ (iConfig.getParameter<std::string>( "muoneps3file" ) ),
    // MuonEps3Path_ (iConfig.getParameter<std::string>( "muoneps3path" ) ),
    genEventInfoToken_(    consumes<GenEventInfoProduct>(        iConfig.getParameter<edm::InputTag>("genEventInfoInput"))),
    LHEEventProductToken_( consumes<LHEEventProduct>(            iConfig.getParameter<edm::InputTag>("LHEEventProductInput"))),
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
    LambdaToken_(   consumes<reco::VertexCompositePtrCandidateCollection>(            iConfig.getParameter<edm::InputTag>("Lambda"))),
    // ,puToken_(      consumes<PileupSummaryInfo>(                                iConfig.getParameter<edm::InputTag>("pileup")))
      puToken_(  consumes<vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("puCollection")))
    ,PhotonToken_(  consumes<reco::ConversionCollection>(edm::InputTag(std::string("reducedEgamma"),std::string("reducedConversions")))) 
    ,beamSpotToken_(     consumes<reco::BeamSpot>(               iConfig.getUntrackedParameter<edm::InputTag>("beamSpot")))
    ,  prefweight_token (consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb")))  //working
    ,rho_token_ (consumes<double> (iConfig.getParameter<edm::InputTag>("rhoCollection")))
    // ,Top_token (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genEvt")))

    // , PrescaleToken_( consumes<pat::PackedTriggerPrescales>(edm::InputTag(std::string("patTrigger"),std::string("")))  )
{
   //now do what ever initialization is needed
    nEvent = 0;
    usesResource("TFileService");

    rc.init(edm::FileInPath(RochString).fullPath());
    lumiWeights_ = new reweight::LumiReWeighting( mcPileupFile_, dataPileupFile_, mcPileupPath_, dataPileupPath_ );

    smalltree = fs->make<TTree>("ttree", "ttree");
    
    // event info
    smalltree->Branch("runNumber",        &runNumber,  "runNumber/I");
    smalltree->Branch("eventNumber",      &eventNumber,"eventNumber/I");
    smalltree->Branch("lumiBlock"  ,      &lumiBlock,  "lumiBlock/I");
    //smalltree->Branch("tree_LHE_Weights", &tree_LHE_Weights);
    //smalltree->Branch("tree_MCEvt_weight", &tree_MCEvt_weight, "tree_MCEvt_weight/F");
    
    hEvents = fs->make<TH1D>("hEvents","hEvents",100,0,100);
    hEvents_with_gen_wt = fs->make<TH1D>("hEvents_with_gen_wt","hEvents_with_gen_wt",2,0,2);
    
    smalltree->Branch("tree_only_gen_wt",&tree_only_gen_wt,"tree_only_gen_wt/D");
    smalltree->Branch("tree_event_weight",&tree_event_weight,"tree_event_weight/D");
    smalltree->Branch("tree_genTop_Weight",&tree_genTop_Weight,"tree_genTop_Weight/D");
    
    smalltree->Branch("PUweight",         &PUweight, "PUweight/D");
    smalltree->Branch("Prefweight",       &Prefweight, "Prefweight/D");
    smalltree->Branch("PU_events", &PU_events, "PU_events/I");
    //smalltree->Branch("AllPU_events_weight", &AllPU_events_weight, "AllPU_events_weight/I");
    smalltree->Branch("tree_only_tigger_filter", &tree_only_tigger_filter);
    smalltree->Branch("tree_Filter",        &tree_Filter);
    smalltree->Branch("tree_FilterSameSign",&tree_FilterSameSign);
    smalltree->Branch("tree_Good_PV",       &tree_Good_PV);
    smalltree->Branch("tree_muon_GenRecoTriggerMatched" ,&tree_muon_GenRecoTriggerMatched);
    smalltree->Branch("tree_Evts_MVAval",   &tree_Evts_MVAval);
    smalltree->Branch("tree_Evts_MVAvalDY",   &tree_Evts_MVAvalDY);
    smalltree->Branch("tree_Evts_MVAvalTT",   &tree_Evts_MVAvalTT);

    // smalltree->Branch("tree_Hemi1_MVAval",    &tree_Hemi1_MVAval);
    // smalltree->Branch("tree_Hemi1_MVAvalDY",  &tree_Hemi1_MVAvalDY);
    // smalltree->Branch("tree_Hemi1_MVAvalTT",  &tree_Hemi1_MVAvalTT);
    // smalltree->Branch("tree_Hemi2_MVAval",    &tree_Hemi2_MVAval);
    // smalltree->Branch("tree_Hemi2_MVAvalDY",  &tree_Hemi2_MVAvalDY);
    // smalltree->Branch("tree_Hemi2_MVAvalTT",  &tree_Hemi2_MVAvalTT);
    //Beamspot
    smalltree->Branch("tree_bs_PosX", &tree_bs_PosX) ;
    smalltree->Branch("tree_bs_PosY", &tree_bs_PosY) ;
    smalltree->Branch("tree_bs_PosZ", &tree_bs_PosZ) ;

    // primary vertex info
    smalltree->Branch("tree_nPV",      &tree_nPV);
    smalltree->Branch("tree_PV_x",     &tree_PV_x);
    smalltree->Branch("tree_PV_y",     &tree_PV_y);
    smalltree->Branch("tree_PV_z",     &tree_PV_z);    
    smalltree->Branch("tree_PV_ez",    &tree_PV_ez);    
    smalltree->Branch("tree_PV_NChi2", &tree_PV_NChi2);    
    smalltree->Branch("tree_PV_ndf",   &tree_PV_ndf);
    smalltree->Branch("tree_PV_rho",   &tree_PV_rho);

    // muons info
    smalltree->Branch("tree_all_nmu",         &tree_all_nmu);
    smalltree->Branch("tree_nmu",             &tree_nmu);
    smalltree->Branch("tree_LT",              &tree_LT);
    smalltree->Branch("tree_Mmumu"  ,         &tree_Mmumu);
    smalltree->Branch("tree_MmumuSameSign"  , &tree_MmumuSameSign);
    smalltree->Branch("tree_muon_isPrompt" ,  &tree_muon_isPrompt);
    smalltree->Branch("tree_muon_pt"  ,       &tree_muon_pt);
    smalltree->Branch("tree_muon_SF" ,        &tree_muon_SF);
    smalltree->Branch("tree_muon_eta" ,       &tree_muon_eta);
    smalltree->Branch("tree_muon_phi" ,       &tree_muon_phi);
    smalltree->Branch("tree_muon_x"  ,        &tree_muon_x);
    smalltree->Branch("tree_muon_y" ,         &tree_muon_y);
    smalltree->Branch("tree_muon_z" ,         &tree_muon_z);
    smalltree->Branch("tree_muon_dxy",        &tree_muon_dxy);
    smalltree->Branch("tree_muon_dxyError",   &tree_muon_dxyError);
    smalltree->Branch("tree_muon_dz",         &tree_muon_dz);
    smalltree->Branch("tree_muon_dzError",    &tree_muon_dzError);
    smalltree->Branch("tree_muon_charge",     &tree_muon_charge);
    smalltree->Branch("tree_muon_isLoose",    &tree_muon_isLoose);
    smalltree->Branch("tree_muon_isMedium",   &tree_muon_isMedium);
    smalltree->Branch("tree_muon_isTight",    &tree_muon_isTight);
    smalltree->Branch("tree_muon_isGlobal",   &tree_muon_isGlobal);
    smalltree->Branch("tree_muon_isoR3",      &tree_muon_isoR3);
    smalltree->Branch("tree_muon_trigger_dimu",  &tree_muon_trigger_dimu);  
    smalltree->Branch("tree_muon_trigger_isomu", &tree_muon_trigger_isomu);
    smalltree->Branch("tree_muon_PFIsoVeryLoose",&tree_muon_PFIsoVeryLoose);
    smalltree->Branch("tree_muon_PFIsoLoose",    &tree_muon_PFIsoLoose);
    smalltree->Branch("tree_muon_PFIsoMedium",   &tree_muon_PFIsoMedium);
    smalltree->Branch("tree_muon_PFIsoTight",    &tree_muon_PFIsoTight);
    smalltree->Branch("tree_muon_MiniIsoLoose",  &tree_muon_MiniIsoLoose);
    smalltree->Branch("tree_muon_MiniIsoMedium", &tree_muon_MiniIsoMedium);
    smalltree->Branch("tree_muon_MiniIsoTight",  &tree_muon_MiniIsoTight);
    smalltree->Branch("tree_muon_TkIsoLoose", &tree_muon_TkIsoLoose);
    smalltree->Branch("tree_muon_TkIsoTight", &tree_muon_TkIsoTight);
    smalltree->Branch("tree_muon_trkLayers",  &tree_muon_trkLayers);
    smalltree->Branch("tree_muon_miniIso",    &tree_muon_miniIso);
    smalltree->Branch("tree_muon_correction", &tree_muon_correction);
    smalltree->Branch("tree_muon_gen",        &tree_muon_gen);
    
    smalltree->Branch("tree_reco_muon_leadingpt",&tree_reco_muon_leadingpt);
    smalltree->Branch("tree_reco_electron_leadingpt2",&tree_reco_electron_leadingpt2);
    smalltree->Branch("tree_reco_muon_leadingeta",&tree_reco_muon_leadingeta);
    smalltree->Branch("tree_reco_electron_leadingeta2",&tree_reco_electron_leadingeta2);
    smalltree->Branch("tree_reco_muon_leadingphi",&tree_reco_muon_leadingphi);
    smalltree->Branch("tree_reco_electron_leadingphi2",&tree_reco_electron_leadingphi2);


    smalltree->Branch("tree_trig_muon_leadingpt",&tree_trig_muon_leadingpt);
    smalltree->Branch("tree_trig_electron_leadingpt2",&tree_trig_electron_leadingpt2);
    smalltree->Branch("tree_trig_muon_leadingeta",&tree_trig_muon_leadingeta);
    smalltree->Branch("tree_trig_electron_leadingeta2",&tree_trig_electron_leadingeta2);
    smalltree->Branch("tree_trig_muon_leadingphi",&tree_trig_muon_leadingphi);
    smalltree->Branch("tree_trig_electron_leadingphi2",&tree_trig_electron_leadingphi2);


    smalltree->Branch("tree_reco_lepton_leadingpt",&tree_reco_lepton_leadingpt);
    smalltree->Branch("tree_reco_lepton_leadingpt2",&tree_reco_lepton_leadingpt2);
    smalltree->Branch("tree_reco_lepton_leadingeta",&tree_reco_lepton_leadingeta);
    smalltree->Branch("tree_reco_lepton_leadingeta2",&tree_reco_lepton_leadingeta2);
    smalltree->Branch("tree_reco_lepton_leadingphi",&tree_reco_lepton_leadingphi);
    smalltree->Branch("tree_reco_lepton_leadingphi2",&tree_reco_lepton_leadingphi2);

    smalltree->Branch("tree_trig_lepton_leadingpt",&tree_trig_lepton_leadingpt);
    smalltree->Branch("tree_trig_lepton_leadingpt2",&tree_trig_lepton_leadingpt2);
    smalltree->Branch("tree_trig_lepton_leadingeta",&tree_trig_lepton_leadingeta);
    smalltree->Branch("tree_trig_lepton_leadingeta2",&tree_trig_lepton_leadingeta2);
    smalltree->Branch("tree_trig_lepton_leadingphi",&tree_trig_lepton_leadingphi);
    smalltree->Branch("tree_trig_lepton_leadingphi2",&tree_trig_lepton_leadingphi2);

    smalltree->Branch("tree_lepton_leadingpt",&tree_lepton_leadingpt);
    smalltree->Branch("tree_lepton_leadingpt2",&tree_lepton_leadingpt2);
    smalltree->Branch("tree_lepton_leadingeta",&tree_lepton_leadingeta);
    smalltree->Branch("tree_lepton_leadingeta2",&tree_lepton_leadingeta2);
    smalltree->Branch("tree_lepton_leadingphi",&tree_lepton_leadingphi);
    smalltree->Branch("tree_lepton_leadingphi2",&tree_lepton_leadingphi2);

    smalltree->Branch("tree_lepton_lepton_dR",&tree_lepton_lepton_dR);
    smalltree->Branch("tree_lepton_lepton_dPhi",&tree_lepton_lepton_dPhi);
    smalltree->Branch("tree_lepton_lepton_dEta",&tree_lepton_lepton_dEta);
 
    smalltree->Branch("tree_ll_pt", &tree_ll_pt);
    smalltree->Branch("tree_ll_eta", &tree_ll_eta);
    smalltree->Branch("tree_ll_phi", &tree_ll_phi);
    smalltree->Branch("tree_ll_px", &tree_ll_px);
    smalltree->Branch("tree_ll_py", &tree_ll_py);
    smalltree->Branch("tree_ll_pz", &tree_ll_pz);
    smalltree->Branch("tree_ll_energy", &tree_ll_energy);
    smalltree->Branch("tree_ll_mass", &tree_ll_mass);

    // electrons info
    smalltree->Branch("tree_all_nel",           &tree_all_nel);
    smalltree->Branch("tree_electron_nEle",     &tree_electron_nEle);
    smalltree->Branch("tree_electron_isPrompt", &tree_electron_isPrompt);
    smalltree->Branch("tree_electron_pt"  ,     &tree_electron_pt);
    smalltree->Branch("tree_electron_eta" ,     &tree_electron_eta);
    smalltree->Branch("tree_electron_phi" ,     &tree_electron_phi);
    smalltree->Branch("tree_electron_x"  ,      &tree_electron_x);
    smalltree->Branch("tree_electron_y" ,       &tree_electron_y);
    smalltree->Branch("tree_electron_z" ,       &tree_electron_z);
    smalltree->Branch("tree_electron_energy",   &tree_electron_energy);
    smalltree->Branch("tree_electron_et",       &tree_electron_et);
    smalltree->Branch("tree_electron_ecal_trk_postcorr", &tree_electron_ecal_trk_postcorr);
    smalltree->Branch("tree_electron_charge",   &tree_electron_charge);
    smalltree->Branch("tree_electron_isoR4",    &tree_electron_isoR4);
    smalltree->Branch("tree_electron_IsLoose",  &tree_electron_IsLoose);
    smalltree->Branch("tree_electron_IsMedium", &tree_electron_IsMedium);
    smalltree->Branch("tree_electron_IsTight",  &tree_electron_IsTight);
    smalltree->Branch("tree_electron_dxy",      &tree_electron_dxy);
    smalltree->Branch("tree_electron_dz",       &tree_electron_dz);
    smalltree->Branch("tree_electron_gen",      &tree_electron_gen);

    // met info
    smalltree->Branch("tree_PFMet_et" ,   &tree_PFMet_et);
    smalltree->Branch("tree_PFMet_phi" ,  &tree_PFMet_phi);
    smalltree->Branch("tree_PFMet_sig" ,  &tree_PFMet_sig);
    smalltree->Branch("tree_PFMet_pt",    &tree_PFMet_pt);
    
    // jet info
    smalltree->Branch("tree_njet"  ,                &tree_njet);
    smalltree->Branch("tree_njetNOmu"  ,            &tree_njetNOmu);
    smalltree->Branch("tree_jet_pt"  ,              &tree_jet_pt);
    smalltree->Branch("tree_jet_eta" ,              &tree_jet_eta);
    smalltree->Branch("tree_jet_phi" ,              &tree_jet_phi);
    smalltree->Branch("tree_jet_px"  ,              &tree_jet_px);
    smalltree->Branch("tree_jet_py"  ,              &tree_jet_py);
    smalltree->Branch("tree_jet_pz"  ,              &tree_jet_pz);
    smalltree->Branch("tree_jet_E"  ,               &tree_jet_E);
    smalltree->Branch("tree_jet_tightid_LepVeto",   &tree_jet_tightid_LepVeto);
    smalltree->Branch("tree_jet_tightid",           &tree_jet_tightid);
    smalltree->Branch("tree_jet_TightJetIDLepVeto", &tree_jet_TightJetIDLepVeto);
    smalltree->Branch("tree_jet_TightJetID"        ,&tree_jet_TightJetID);
    smalltree->Branch("tree_jet_HadronFlavour",     &tree_jet_HadronFlavour);
    smalltree->Branch("tree_jet_pileupID",          &tree_jet_pileupID);
    smalltree->Branch("tree_jet_btag_DeepCSV",      &tree_jet_btag_DeepCSV);
    smalltree->Branch("tree_jet_btag_DeepJet",      &tree_jet_btag_DeepJet);
    smalltree->Branch("tree_jet_leadingpt",         &tree_jet_leadingpt);
    smalltree->Branch("tree_jet_leadingpt2",        &tree_jet_leadingpt2);
    smalltree->Branch("tree_jet_leadingeta",        &tree_jet_leadingeta);
    smalltree->Branch("tree_jet_leadingeta2",       &tree_jet_leadingeta2);
    smalltree->Branch("tree_jet_leadingMuon_dR",    &tree_jet_leadingMuon_dR);
    smalltree->Branch("tree_jet_leadingMuon2_dR",   &tree_jet_leadingMuon2_dR);
    smalltree->Branch("tree_jet_jet_dR",            &tree_jet_jet_dR);
    smalltree->Branch("tree_jet_jet_dPhi",          &tree_jet_jet_dPhi);
    smalltree->Branch("tree_jet_jet_dEta",          &tree_jet_jet_dEta);
    smalltree->Branch("tree_muon_jet_dRmin",        &tree_muon_jet_dRmin);
    smalltree->Branch("tree_muon_jet_dRmax",        &tree_muon_jet_dRmax);
    smalltree->Branch("tree_elemu_jet_dRmin",       &tree_elemu_jet_dRmin);
    smalltree->Branch("tree_elemu_jet_dRmax",       &tree_elemu_jet_dRmax);
    smalltree->Branch("tree_ele_jet_dRmin",         &tree_ele_jet_dRmin);
    smalltree->Branch("tree_ele_jet_dRmax",         &tree_ele_jet_dRmax);
    smalltree->Branch("tree_HT"  ,                  &tree_HT);
    
    smalltree->Branch("tree_nK0",           &tree_nK0);
    smalltree->Branch("tree_K0_x",          &tree_K0_x);
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

    smalltree->Branch("tree_nLambda",       &tree_nLambda);
    smalltree->Branch("tree_L0_x",          &tree_L0_x);
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

    // reco V0
    smalltree->Branch("tree_nV0_reco",           &tree_nV0_reco);
    smalltree->Branch("tree_V0_reco_x",          &tree_V0_reco_x); // l'index 0 donne le PV!
    smalltree->Branch("tree_V0_reco_y",          &tree_V0_reco_y);
    smalltree->Branch("tree_V0_reco_z",          &tree_V0_reco_z);
    smalltree->Branch("tree_V0_reco_r",          &tree_V0_reco_r);
    smalltree->Branch("tree_V0_reco_drSig",      &tree_V0_reco_drSig);
    smalltree->Branch("tree_V0_reco_dzSig",      &tree_V0_reco_dzSig);
    smalltree->Branch("tree_V0_reco_angleXY",    &tree_V0_reco_angleXY);
    smalltree->Branch("tree_V0_reco_angleZ",     &tree_V0_reco_angleZ);
    smalltree->Branch("tree_V0_reco_NChi2",      &tree_V0_reco_NChi2);
    smalltree->Branch("tree_V0_reco_ndf",        &tree_V0_reco_ndf);
    smalltree->Branch("tree_V0_reco_mass",       &tree_V0_reco_mass);
    smalltree->Branch("tree_V0_reco_pt",         &tree_V0_reco_pt);
    smalltree->Branch("tree_V0_reco_eta",        &tree_V0_reco_eta);
    smalltree->Branch("tree_V0_reco_phi",        &tree_V0_reco_phi);
    smalltree->Branch("tree_V0_reco_source",     &tree_V0_reco_source);
    smalltree->Branch("tree_V0_reco_badTkHit",   &tree_V0_reco_badTkHit);
    smalltree->Branch("tree_V0_reco_dca",        &tree_V0_reco_dca);

    // reco Secondary Interaction
    smalltree->Branch("tree_nSecInt",           &tree_nSecInt);
    smalltree->Branch("tree_SecInt_x",	        &tree_SecInt_x);
    smalltree->Branch("tree_SecInt_y",	        &tree_SecInt_y);
    smalltree->Branch("tree_SecInt_z",	        &tree_SecInt_z);
    smalltree->Branch("tree_SecInt_r",	        &tree_SecInt_r);
    smalltree->Branch("tree_SecInt_d",	        &tree_SecInt_d);
    smalltree->Branch("tree_SecInt_drSig",      &tree_SecInt_drSig);
    smalltree->Branch("tree_SecInt_dzSig",      &tree_SecInt_dzSig);
    smalltree->Branch("tree_SecInt_angleXY",    &tree_SecInt_angleXY);
    smalltree->Branch("tree_SecInt_angleZ",     &tree_SecInt_angleZ);
    smalltree->Branch("tree_SecInt_NChi2",      &tree_SecInt_NChi2);
    smalltree->Branch("tree_SecInt_ndf",        &tree_SecInt_ndf);
    smalltree->Branch("tree_SecInt_mass",       &tree_SecInt_mass);
    smalltree->Branch("tree_SecInt_pt",         &tree_SecInt_pt);
    smalltree->Branch("tree_SecInt_eta",        &tree_SecInt_eta);
    smalltree->Branch("tree_SecInt_phi",        &tree_SecInt_phi);
    smalltree->Branch("tree_SecInt_charge",     &tree_SecInt_charge);
    smalltree->Branch("tree_SecInt_badTkHit",   &tree_SecInt_badTkHit);
    smalltree->Branch("tree_SecInt_dca",        &tree_SecInt_dca);
    smalltree->Branch("tree_SecInt_selec",      &tree_SecInt_selec);
    smalltree->Branch("tree_SecInt_layer",      &tree_SecInt_layer);
    smalltree->Branch("tree_SecInt_LLP",        &tree_SecInt_LLP);
    smalltree->Branch("tree_SecInt_LLP_dr",     &tree_SecInt_LLP_dr);
    smalltree->Branch("tree_SecInt_LLP_dz",     &tree_SecInt_LLP_dz);
    smalltree->Branch("tree_SecInt_LLP_dd",     &tree_SecInt_LLP_dd);
    smalltree->Branch("tree_SecInt_tk1",        &tree_SecInt_tk1);
    smalltree->Branch("tree_SecInt_tk2",        &tree_SecInt_tk2);

    smalltree->Branch("tree_nYConv",        &tree_nYConv);
    smalltree->Branch("tree_Yc_x",          &tree_Yc_x); 
    smalltree->Branch("tree_Yc_y",          &tree_Yc_y);
    smalltree->Branch("tree_Yc_z",          &tree_Yc_z);
    smalltree->Branch("tree_Yc_r",          &tree_Yc_r);
    smalltree->Branch("tree_Yc_dr0",        &tree_Yc_dr0);
    smalltree->Branch("tree_Yc_dr1",        &tree_Yc_dr1);
    smalltree->Branch("tree_Yc_dz0",        &tree_Yc_dz0);
    smalltree->Branch("tree_Yc_dz1",        &tree_Yc_dz1);
    smalltree->Branch("tree_Yc_costheta",   &tree_Yc_costheta);
    smalltree->Branch("tree_Yc_layer",      &tree_Yc_layer);
    smalltree->Branch("tree_Yc_NChi2",      &tree_Yc_NChi2);
    smalltree->Branch("tree_Yc_ndf",        &tree_Yc_ndf);
    smalltree->Branch("tree_Yc_nDaughters", &tree_Yc_nDaughters);
    smalltree->Branch("tree_Yc_pt",         &tree_Yc_pt);
    smalltree->Branch("tree_Yc_eta",        &tree_Yc_eta);
    smalltree->Branch("tree_Yc_phi",        &tree_Yc_phi);
    smalltree->Branch("tree_Yc_mass",       &tree_Yc_mass);
    smalltree->Branch("tree_Yc_ntracks",       &tree_Yc_ntracks);
    smalltree->Branch("tree_Yc_tracks_index",  &tree_Yc_tracks_index);
    smalltree->Branch("tree_Yc_tracks_charge", &tree_Yc_tracks_charge);
    smalltree->Branch("tree_Yc_tracks_pt",     &tree_Yc_tracks_pt);
    smalltree->Branch("tree_Yc_tracks_eta",    &tree_Yc_tracks_eta);
    smalltree->Branch("tree_Yc_tracks_phi",    &tree_Yc_tracks_phi);
    smalltree->Branch("tree_Yc_tracks_phi0",   &tree_Yc_tracks_phi0);
        
    // track
    smalltree->Branch("tree_TRACK_SIZE",         &tree_TRACK_SIZE, "tree_TRACK_SIZE/I");
    smalltree->Branch("tree_nTracks",            &tree_nTracks, "tree_nTracks/I"); 
    smalltree->Branch("tree_nLostTracks",        &tree_nLostTracks, "tree_nLostTracks/I"); 
//     smalltree->Branch("tree_passesTrkPtr",       &tree_passesTrkPtr);
    smalltree->Branch("tree_track_ipc",          &tree_track_ipc);
    smalltree->Branch("tree_track_lost",         &tree_track_lost);
    smalltree->Branch("tree_track_px",           &tree_track_px);
    smalltree->Branch("tree_track_py",           &tree_track_py);
    smalltree->Branch("tree_track_pz",           &tree_track_pz);
    smalltree->Branch("tree_track_pt",           &tree_track_pt);
    smalltree->Branch("tree_track_eta",          &tree_track_eta );
    smalltree->Branch("tree_track_phi",          &tree_track_phi );
    smalltree->Branch("tree_track_charge",       &tree_track_charge );
    smalltree->Branch("tree_track_NChi2",        &tree_track_NChi2);
    smalltree->Branch("tree_track_isHighPurity", &tree_track_isHighPurity);
    smalltree->Branch("tree_track_dxy",          &tree_track_dxy );
    smalltree->Branch("tree_track_dxyError",     &tree_track_dxyError);
    smalltree->Branch("tree_track_drSig",        &tree_track_drSig);
    smalltree->Branch("tree_track_dz",           &tree_track_dz);
    smalltree->Branch("tree_track_dzError",      &tree_track_dzError  );
    smalltree->Branch("tree_track_dzSig",        &tree_track_dzSig);
//$$$$$$
//     smalltree->Branch("tree_track_dzTOpu",       &tree_track_dzTOpu);
//     smalltree->Branch("tree_track_dzSigTOpu",    &tree_track_dzSigTOpu  );
//$$$$$$
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
    smalltree->Branch("tree_track_btag",         &tree_track_btag);
    smalltree->Branch("tree_track_energy",       &tree_track_energy);

    smalltree->Branch("tree_track_Hemi",           &tree_track_Hemi);
    smalltree->Branch("tree_track_Hemi_dR",        &tree_track_Hemi_dR);
    smalltree->Branch("tree_track_Hemi_dRmax",     &tree_track_Hemi_dRmax);
    smalltree->Branch("tree_track_Hemi_mva_NChi2", &tree_track_Hemi_mva_NChi2);
    smalltree->Branch("tree_track_Hemi_ping",      &tree_track_Hemi_ping);
    smalltree->Branch("tree_track_Hemi_dFirstVtx", &tree_track_Hemi_dFirstVtx);
    smalltree->Branch("tree_track_Hemi_LLP",       &tree_track_Hemi_LLP);
//$$$$$$
//     smalltree->Branch("tree_track_Hemi_d0",        &tree_track_Hemi_d0);
//     smalltree->Branch("tree_track_Hemi_d0Sig",     &tree_track_Hemi_d0Sig);
//     smalltree->Branch("tree_track_HemiOp_d0",      &tree_track_HemiOp_d0);
//     smalltree->Branch("tree_track_HemiOp_d0Sig",   &tree_track_HemiOp_d0Sig);
//$$$$$$
    smalltree->Branch("tree_track_Hemi_isjet",&tree_track_Hemi_isjet);
    
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

//$$$$
    // tracks from V0 (and SI) candidate
    smalltree->Branch("tree_V0_track_isFromV0",     &tree_V0_track_isFromV0);
    smalltree->Branch("tree_V0_track_isFromSI",     &tree_V0_track_isFromSI);
    smalltree->Branch("tree_V0_track_lost",         &tree_V0_track_lost);
    smalltree->Branch("tree_V0_track_pt",           &tree_V0_track_pt);
    smalltree->Branch("tree_V0_track_eta",          &tree_V0_track_eta );
    smalltree->Branch("tree_V0_track_phi",          &tree_V0_track_phi );
    smalltree->Branch("tree_V0_track_charge",       &tree_V0_track_charge );
    smalltree->Branch("tree_V0_track_NChi2",        &tree_V0_track_NChi2);
    smalltree->Branch("tree_V0_track_dxy",          &tree_V0_track_dxy );
    smalltree->Branch("tree_V0_track_drSig",        &tree_V0_track_drSig);
    smalltree->Branch("tree_V0_track_dz",           &tree_V0_track_dz);
    smalltree->Branch("tree_V0_track_dzSig",        &tree_V0_track_dzSig);
    smalltree->Branch("tree_V0_track_nHit",         &tree_V0_track_nHit);
    smalltree->Branch("tree_V0_track_nHitPixel",    &tree_V0_track_nHitPixel);
    smalltree->Branch("tree_V0_track_firstHit",     &tree_V0_track_firstHit);
    smalltree->Branch("tree_V0_track_firstHit_x",   &tree_V0_track_firstHit_x);
    smalltree->Branch("tree_V0_track_firstHit_y",   &tree_V0_track_firstHit_y);
    smalltree->Branch("tree_V0_track_firstHit_z",   &tree_V0_track_firstHit_z);
    smalltree->Branch("tree_V0_track_iJet",         &tree_V0_track_iJet);
    smalltree->Branch("tree_V0_track_ntrk10",       &tree_V0_track_ntrk10);
    smalltree->Branch("tree_V0_track_ntrk20",       &tree_V0_track_ntrk20);
    smalltree->Branch("tree_V0_track_ntrk30",       &tree_V0_track_ntrk30);
    smalltree->Branch("tree_V0_track_ntrk40",       &tree_V0_track_ntrk40);
    smalltree->Branch("tree_V0_track_Hemi",         &tree_V0_track_Hemi);
    smalltree->Branch("tree_V0_track_Hemi_dR",      &tree_V0_track_Hemi_dR);
    smalltree->Branch("tree_V0_track_Hemi_dRmax",   &tree_V0_track_Hemi_dRmax);
//$$$$

    // gen info
    smalltree->Branch("tree_GenPVx" ,  &tree_GenPVx);
    smalltree->Branch("tree_GenPVy" ,  &tree_GenPVy);
    smalltree->Branch("tree_GenPVz" ,  &tree_GenPVz);

    smalltree->Branch("tree_smu_mass" ,  &tree_smu_mass);
    smalltree->Branch("tree_neu_mass" ,  &tree_neu_mass);
    smalltree->Branch("tree_neu_ctau" ,  &tree_neu_ctau);
    
    smalltree->Branch("tree_genParticle_pt"  ,          &tree_genParticle_pt);
    smalltree->Branch("tree_genParticle_eta" ,          &tree_genParticle_eta);
    smalltree->Branch("tree_genParticle_phi" ,          &tree_genParticle_phi);
    smalltree->Branch("tree_genParticle_charge" ,       &tree_genParticle_charge);
    smalltree->Branch("tree_genParticle_pdgId" ,        &tree_genParticle_pdgId);
    smalltree->Branch("tree_genParticle_mass" ,         &tree_genParticle_mass);
    smalltree->Branch("tree_genParticle_x"  ,	        &tree_genParticle_x);
    smalltree->Branch("tree_genParticle_y" ,	        &tree_genParticle_y);
    smalltree->Branch("tree_genParticle_z" ,	        &tree_genParticle_z);
    smalltree->Branch("tree_genParticle_px",            &tree_genParticle_px);
    smalltree->Branch("tree_genParticle_py",            &tree_genParticle_py);
    smalltree->Branch("tree_genParticle_pz",            &tree_genParticle_pz);
    smalltree->Branch("tree_genParticle_energy",        &tree_genParticle_energy);
    smalltree->Branch("tree_genParticle_isPromptFinalState", &tree_genParticle_isPromptFinalState);
    smalltree->Branch("tree_genParticle_statusCode",    &tree_genParticle_statusCode);
    smalltree->Branch("tree_genParticle_mother_pdgId" , &tree_genParticle_mother_pdgId);
    smalltree->Branch("tree_genParticle_LLP" ,          &tree_genParticle_LLP);
    smalltree->Branch("tree_genParticle_ct" ,	        &tree_genParticle_ct);
    smalltree->Branch("tree_genParticle_ct0" ,	        &tree_genParticle_ct0);

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
    smalltree->Branch("tree_genAxis_dPhineuneu",     &tree_genAxis_dPhineuneu);
    smalltree->Branch("tree_genAxis_dEtaneuneu",     &tree_genAxis_dEtaneuneu);
    smalltree->Branch("tree_GenAxes_Mass",           &tree_GenAxes_Mass);
    smalltree->Branch("tree_GenAxes_CombinedHemiLeptonMass",&tree_GenAxes_CombinedHemiLeptonMass );
    smalltree->Branch("tree_GenAxis_Neu_dRmin",      &tree_GenAxis_Neu_dRmin);
    smalltree->Branch("tree_GenAxis_Neu_dRmax",      &tree_GenAxis_Neu_dRmax);
    smalltree->Branch("tree_GenAxis_RecoAxis_dRmin", &tree_GenAxis_RecoAxis_dRmin);
    smalltree->Branch("tree_GenAxis_RecoAxis_dRmax", &tree_GenAxis_RecoAxis_dRmax);

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
    smalltree->Branch("tree_genFromB_dd",&tree_genFromB_dd);
    smalltree->Branch("tree_genFromB_dr",&tree_genFromB_dr);
    smalltree->Branch("tree_genFromB_dz",&tree_genFromB_dz);
    
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
    smalltree->Branch("tree_LLP_r",         &tree_LLP_r);
    smalltree->Branch("tree_LLP_dist",      &tree_LLP_dist);
    smalltree->Branch("tree_LLP_nTrks",     &tree_LLP_nTrks);
    smalltree->Branch("tree_LLP_Vtx_nTrks", &tree_LLP_Vtx_nTrks);
    smalltree->Branch("tree_LLP_Vtx_NChi2", &tree_LLP_Vtx_NChi2);
    smalltree->Branch("tree_LLP_Vtx_dx",    &tree_LLP_Vtx_dx);
    smalltree->Branch("tree_LLP_Vtx_dy",    &tree_LLP_Vtx_dy);
    smalltree->Branch("tree_LLP_Vtx_dz",    &tree_LLP_Vtx_dz);
    smalltree->Branch("tree_LLP_Vtx_dist",  &tree_LLP_Vtx_dist);
    smalltree->Branch("tree_LLP_Vtx_dd",    &tree_LLP_Vtx_dd);
    smalltree->Branch("tree_LLP12_dR",      &tree_LLP12_dR);
    smalltree->Branch("tree_LLP12_deta",    &tree_LLP12_deta);
    smalltree->Branch("tree_LLP12_dphi",    &tree_LLP12_dphi);
    smalltree->Branch("tree_LLP_Mass",      &tree_LLP_Mass);

    smalltree->Branch("tree_Hemi",       &tree_Hemi);
    smalltree->Branch("tree_Hemi_njet",  &tree_Hemi_njet);
    smalltree->Branch("tree_Hemi_njet_nomu",  &tree_Hemi_njet_nomu);
    smalltree->Branch("tree_Hemi_pt",    &tree_Hemi_pt);
    smalltree->Branch("tree_Hemi_eta",   &tree_Hemi_eta);
    smalltree->Branch("tree_Hemi_phi",   &tree_Hemi_phi);
    smalltree->Branch("tree_Hemi_nTrks", &tree_Hemi_nTrks);
    smalltree->Branch("tree_Hemi_nTrks_sig", &tree_Hemi_nTrks_sig);
    smalltree->Branch("tree_Hemi_nTrks_bad", &tree_Hemi_nTrks_bad);
    smalltree->Branch("tree_Hemi_mass",     &tree_Hemi_mass);
    smalltree->Branch("tree_HemiMu_mass",   &tree_HemiMu_mass);
    smalltree->Branch("tree_HemiMu_pt",     &tree_HemiMu_pt);
    smalltree->Branch("tree_HemiMu_dR",     &tree_HemiMu_dR);
    smalltree->Branch("tree_HemiMuOp_mass", &tree_HemiMuOp_mass);
    smalltree->Branch("tree_HemiMuOp_pt",   &tree_HemiMuOp_pt);
    smalltree->Branch("tree_HemiMuOp_dR",   &tree_HemiMuOp_dR);
    smalltree->Branch("tree_Hemi_LooseBTag_axes",&tree_Hemi_LooseBTag_axes);
    smalltree->Branch("tree_Hemi_MediumBTag_axes",&tree_Hemi_MediumBTag_axes);
    smalltree->Branch("tree_Hemi_TightBTag_axes",&tree_Hemi_TightBTag_axes);
    smalltree->Branch("tree_Hemi_dR12",      &tree_Hemi_dR12);

    smalltree->Branch("tree_Hemi_LLP",       &tree_Hemi_LLP);
    smalltree->Branch("tree_Hemi_LLP_pt",    &tree_Hemi_LLP_pt);
    smalltree->Branch("tree_Hemi_LLP_eta",   &tree_Hemi_LLP_eta);
    smalltree->Branch("tree_Hemi_LLP_phi",   &tree_Hemi_LLP_phi);
    smalltree->Branch("tree_Hemi_LLP_dist",  &tree_Hemi_LLP_dist);
    smalltree->Branch("tree_Hemi_LLP_x",     &tree_Hemi_LLP_x);
    smalltree->Branch("tree_Hemi_LLP_y",     &tree_Hemi_LLP_y);
    smalltree->Branch("tree_Hemi_LLP_z",     &tree_Hemi_LLP_z);
    smalltree->Branch("tree_Hemi_LLP_dR",    &tree_Hemi_LLP_dR);
    smalltree->Branch("tree_Hemi_LLP_mother",&tree_Hemi_LLP_mother);
    smalltree->Branch("tree_Hemi_LLP_Vtx_dx",     &tree_Hemi_LLP_Vtx_dx);
    smalltree->Branch("tree_Hemi_LLP_Vtx_dy",     &tree_Hemi_LLP_Vtx_dy);
    smalltree->Branch("tree_Hemi_LLP_Vtx_dz",     &tree_Hemi_LLP_Vtx_dz);
    smalltree->Branch("tree_Hemi_LLP_Vtx_dr",     &tree_Hemi_LLP_Vtx_dr);
    smalltree->Branch("tree_Hemi_LLP_muOK_dR",    &tree_Hemi_LLP_muOK_dR);
    smalltree->Branch("tree_Hemi_LLP_muOK_pt",    &tree_Hemi_LLP_muOK_pt);
    smalltree->Branch("tree_Hemi_LLP_muOK_mass",  &tree_Hemi_LLP_muOK_mass);
    smalltree->Branch("tree_Hemi_LLP_muNO_dR",    &tree_Hemi_LLP_muNO_dR);
    smalltree->Branch("tree_Hemi_LLP_muNO_pt",    &tree_Hemi_LLP_muNO_pt);
    smalltree->Branch("tree_Hemi_LLP_muNO_mass",  &tree_Hemi_LLP_muNO_mass);
    smalltree->Branch("tree_Hemi_LLP_dR12",  &tree_Hemi_LLP_dR12);
    smalltree->Branch("tree_Hemi_LLP_ping",  &tree_Hemi_LLP_ping);
    smalltree->Branch("tree_event_LLP_ping", &tree_event_LLP_ping);
    
    smalltree->Branch("tree_Hemi_Vtx_step",  &tree_Hemi_Vtx_step);
    smalltree->Branch("tree_Hemi_Vtx_isTight",&tree_Hemi_Vtx_isTight);
    smalltree->Branch("tree_Hemi_Vtx_NChi2", &tree_Hemi_Vtx_NChi2);
    smalltree->Branch("tree_Hemi_Vtx_nTrks", &tree_Hemi_Vtx_nTrks);
    smalltree->Branch("tree_Hemi_Vtx_nTrks_sig", &tree_Hemi_Vtx_nTrks_sig);
    smalltree->Branch("tree_Hemi_Vtx_nTrks_bad", &tree_Hemi_Vtx_nTrks_bad);
    smalltree->Branch("tree_Hemi_Vtx_x",     &tree_Hemi_Vtx_x);
    smalltree->Branch("tree_Hemi_Vtx_y",     &tree_Hemi_Vtx_y);
    smalltree->Branch("tree_Hemi_Vtx_z",     &tree_Hemi_Vtx_z);
    smalltree->Branch("tree_Hemi_Vtx_r",     &tree_Hemi_Vtx_r);
    smalltree->Branch("tree_Hemi_Vtx_dR",    &tree_Hemi_Vtx_dR);
    smalltree->Branch("tree_Hemi_Vtx_xError",&tree_Hemi_Vtx_xError);
    smalltree->Branch("tree_Hemi_Vtx_yError",&tree_Hemi_Vtx_yError);
    smalltree->Branch("tree_Hemi_Vtx_zError",&tree_Hemi_Vtx_zError);
    smalltree->Branch("tree_Hemi_Vtx_BTag",  &tree_Hemi_Vtx_BTag);
    smalltree->Branch("tree_Hemi_Vtx_trackWeight",&tree_Hemi_Vtx_trackWeight);
    smalltree->Branch("tree_Hemi_Vtx_SumtrackWeight",&tree_Hemi_Vtx_SumtrackWeight);
    smalltree->Branch("tree_Hemi_Vtx_track_MeanDCA_d",&tree_Hemi_Vtx_track_MeanDCA_d);
    smalltree->Branch("tree_Hemi_Vtx_Mass", &tree_Hemi_Vtx_Mass);
    smalltree->Branch("tree_Hemi_Vtx_dist",  &tree_Hemi_Vtx_dist);
    smalltree->Branch("tree_Hemi_Vtx_ntrk10",&tree_Hemi_Vtx_ntrk10);
    smalltree->Branch("tree_Hemi_Vtx_ntrk20",&tree_Hemi_Vtx_ntrk20);
    smalltree->Branch("tree_event_nVtx",      &tree_event_nVtx);
    smalltree->Branch("tree_event_Vtx_Vtx_dr",&tree_event_Vtx_Vtx_dr);
    smalltree->Branch("tree_event_Vtx_Vtx_dz",&tree_event_Vtx_Vtx_dz);
    smalltree->Branch("tree_event_Vtx_Vtx_dd",&tree_event_Vtx_Vtx_dd);
    smalltree->Branch("tree_event_Vtx_Vtx_reldd",&tree_event_Vtx_Vtx_reldd);
    smalltree->Branch("tree_event_Vtx_Vtx_dR",&tree_event_Vtx_Vtx_dR);
    smalltree->Branch("tree_event_Vtx_Vtx_step",&tree_event_Vtx_Vtx_step);

    smalltree->Branch("tree_Hemi_SecLLP",&tree_Hemi_SecLLP);
    smalltree->Branch("tree_Hemi_LLP_SecVtx_dx",&tree_Hemi_LLP_SecVtx_dx);
    smalltree->Branch("tree_Hemi_LLP_SecVtx_dy",&tree_Hemi_LLP_SecVtx_dy);
    smalltree->Branch("tree_Hemi_LLP_SecVtx_dz",&tree_Hemi_LLP_SecVtx_dz);
    smalltree->Branch("tree_Hemi_LLP_SecVtx_dr",&tree_Hemi_LLP_SecVtx_dr);
    smalltree->Branch("tree_Hemi_SecLLP_ping",&tree_Hemi_SecLLP_ping);
    smalltree->Branch("tree_event_SecLLP_ping",&tree_event_SecLLP_ping);

    // smalltree->Branch("tree_Hemi_SecVtx_track_DCA_x",&tree_Hemi_SecVtx_track_DCA_x);
    // smalltree->Branch("tree_Hemi_SecVtx_track_DCA_y",&tree_Hemi_SecVtx_track_DCA_y);
    // smalltree->Branch("tree_Hemi_SecVtx_track_DCA_z",&tree_Hemi_SecVtx_track_DCA_z);
    // smalltree->Branch("tree_Hemi_SecVtx_track_DCA_r",&tree_Hemi_SecVtx_track_DCA_r);
    // smalltree->Branch("tree_Hemi_SecVtx_track_DCA_d",&tree_Hemi_SecVtx_track_DCA_d);
    smalltree->Branch("tree_Hemi_SecVtx",     &tree_Hemi_SecVtx);
    smalltree->Branch("tree_Hemi_SecVtx_step",&tree_Hemi_SecVtx_step);
    smalltree->Branch("tree_Hemi_SecVtx_x",&tree_Hemi_SecVtx_x);
    smalltree->Branch("tree_Hemi_SecVtx_y",&tree_Hemi_SecVtx_y);
    smalltree->Branch("tree_Hemi_SecVtx_z",&tree_Hemi_SecVtx_z);
    smalltree->Branch("tree_Hemi_SecVtx_r",&tree_Hemi_SecVtx_r);
    smalltree->Branch("tree_Hemi_SecVtx_dR",&tree_Hemi_SecVtx_dR);
    smalltree->Branch("tree_Hemi_SecVtx_nTrks",&tree_Hemi_SecVtx_nTrks);
    smalltree->Branch("tree_Hemi_SecVtx_NChi2", &tree_Hemi_SecVtx_NChi2);
    smalltree->Branch("tree_Hemi_SecVtx_dist",&tree_Hemi_SecVtx_dist);
    smalltree->Branch("tree_Hemi_SecVtx_track_MeanDCA_d",&tree_Hemi_SecVtx_track_MeanDCA_d);
    smalltree->Branch("tree_Hemi_SecVtx_SumtrackWeight",&tree_Hemi_SecVtx_SumtrackWeight);
    smalltree->Branch("tree_Hemi_SecVtx_trackWeight",&tree_Hemi_SecVtx_trackWeight);
    smalltree->Branch("tree_Hemi_SecVtx_Mass",&tree_Hemi_SecVtx_Mass);
    smalltree->Branch("tree_event_MergedVtx_Vtx_dr",&tree_event_MergedVtx_Vtx_dr);
    smalltree->Branch("tree_event_MergedVtx_Vtx_dz",&tree_event_MergedVtx_Vtx_dz);
    smalltree->Branch("tree_event_MergedVtx_Vtx_dd",&tree_event_MergedVtx_Vtx_dd);
    smalltree->Branch("tree_event_MergedVtx_Vtx_reldd",&tree_event_MergedVtx_Vtx_reldd);
    smalltree->Branch("tree_event_MergedVtx_Vtx_dR",&tree_event_MergedVtx_Vtx_dR);
    smalltree->Branch("tree_event_MergedVtx_Vtx_step",&tree_event_MergedVtx_Vtx_step);

        // smalltree->Branch("tree_Hemi_Vtx_track_DCA_x",&tree_Hemi_Vtx_track_DCA_x);
    // smalltree->Branch("tree_Hemi_Vtx_track_DCA_y",&tree_Hemi_Vtx_track_DCA_y);
    // smalltree->Branch("tree_Hemi_Vtx_track_DCA_z",&tree_Hemi_Vtx_track_DCA_z);
    // smalltree->Branch("tree_Hemi_Vtx_track_DCA_r",&tree_Hemi_Vtx_track_DCA_r);
    // smalltree->Branch("tree_Hemi_Vtx_track_DCA_d",&tree_Hemi_Vtx_track_DCA_d);
        smalltree->Branch("tree_Hemi_Vtx_BDT_nTrks",&tree_Hemi_Vtx_BDT_nTrks);
    smalltree->Branch("tree_Hemi_Vtx_BDT_NChi2",&tree_Hemi_Vtx_BDT_NChi2);
    smalltree->Branch("tree_Hemi_Vtx_BDT_step",&tree_Hemi_Vtx_BDT_step);
    smalltree->Branch("tree_Hemi_Vtx_BDT_STW",&tree_Hemi_Vtx_BDT_STW);
    smalltree->Branch("tree_Hemi_Vtx_BDT_Mass",&tree_Hemi_Vtx_BDT_Mass);
    smalltree->Branch("tree_Hemi_Vtx_BDT_HMass",&tree_Hemi_Vtx_BDT_HMass);
    smalltree->Branch("tree_Hemi_Vtx_BDT_ntrk10",&tree_Hemi_Vtx_BDT_ntrk10);
    smalltree->Branch("tree_Hemi_Vtx_BDT_ntrk20",&tree_Hemi_Vtx_BDT_ntrk20);
    smalltree->Branch("tree_Hemi_Vtx_BDT_MeanDCA",&tree_Hemi_Vtx_BDT_MeanDCA);
    smalltree->Branch("tree_Hemi_Vtx_MVAval_Loose", &tree_Hemi_Vtx_MVAval_Loose);
    smalltree->Branch("tree_Hemi_Vtx_MVAval_Tight",&tree_Hemi_Vtx_MVAval_Tight);

    // ----------------Trigger Muon + dilepton-------------
    smalltree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
    smalltree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v);
    smalltree->Branch("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v);
    smalltree->Branch("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v);
    smalltree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v);
    smalltree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v);
    smalltree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
    smalltree->Branch("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);
    smalltree->Branch("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);
    smalltree->Branch("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);
    smalltree->Branch("HLT_Ele27_WPTight_Gsf_v",&HLT_Ele27_WPTight_Gsf_v);
    smalltree->Branch("HLT_Ele32_WPTight_Gsf_v",&HLT_Ele32_WPTight_Gsf_v);
    smalltree->Branch("HLT_Ele35_WPTight_Gsf_v", &HLT_Ele35_WPTight_Gsf_v);
    smalltree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v);
    smalltree->Branch("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
    // ----------------Trigger PFMET-------------                                                                                                                             
    smalltree->Branch("HLT_PFMET120_PFMHT120_IDTight_v",&HLT_PFMET120_PFMHT120_IDTight_v);
    smalltree->Branch("HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",&HLT_PFMET120_PFMHT120_IDTight_PFHT60_v);
    smalltree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v);
    smalltree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",&HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v);
    smalltree->Branch("HLT_PFMET250_HBHECleaned_v",&HLT_PFMET250_HBHECleaned_v);
    smalltree->Branch("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",&HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v);
    smalltree->Branch("HLT_IsoMu24_v",&HLT_IsoMu24_v);
    smalltree->Branch("HLT_IsoMu27_v",&HLT_IsoMu27_v);
    smalltree->Branch("HLT_IsoTkMu24_v", &HLT_IsoTkMu24_v);



    //----------------------------------------
    // - BDT Input Variables -----------------
    //----------------------------------------

    //---------------EVTS--------------------
    readerEvts->AddVariable( "mva_Evts_MET_et",             &mva_Evts_MET_et);
    readerEvts->AddVariable( "mva_Evts_nTrks",              &mva_Evts_nTrks);
    // readerEvts->AddVariable( "mva_Evts_muon1_pt",           &mva_Evts_muon1_pt);
    // readerEvts->AddVariable( "mva_Evts_muon2_pt",           &mva_Evts_muon2_pt);
    readerEvts->AddVariable( "mva_Evts_jet1_pt",            &mva_Evts_jet1_pt);
    readerEvts->AddVariable( "mva_Evts_jet2_pt",            &mva_Evts_jet2_pt);
    readerEvts->AddVariable("mva_Evts_jet1_eta",            &mva_Evts_jet1_eta);
    readerEvts->AddVariable("mva_Evts_jet2_eta",            &mva_Evts_jet2_eta);
    readerEvts->AddVariable( "mva_Evts_jet12_dR",           &mva_Evts_jet12_dR);
    readerEvts->AddVariable( "mva_Evts_jet12_dPhi",         &mva_Evts_jet12_dPhi);
    readerEvts->AddVariable( "mva_Evts_jet12_dEta",         &mva_Evts_jet12_dEta);
    readerEvts->AddVariable("mva_Evts_Hemi1_njet_nomu",     &mva_Evts_Hemi1_njet_nomu);
    readerEvts->AddVariable("mva_Evts_Hemi2_njet_nomu",     &mva_Evts_Hemi2_njet_nomu);
    readerEvts->AddVariable("mva_Evts_Hemi1_pt",            &mva_Evts_Hemi1_pt);
    readerEvts->AddVariable("mva_Evts_Hemi2_pt",            &mva_Evts_Hemi2_pt);
    readerEvts->AddVariable("mva_Evts_Hemi1_eta",           &mva_Evts_Hemi1_eta);
    readerEvts->AddVariable("mva_Evts_Hemi2_eta",           &mva_Evts_Hemi2_eta);
    readerEvts->AddVariable("mva_Evts_Hemi1_phi",           &mva_Evts_Hemi1_phi);
    readerEvts->AddVariable("mva_Evts_Hemi2_phi",           &mva_Evts_Hemi2_phi);
    readerEvts->AddVariable("mva_Evts_Hemi1_nTrks",         &mva_Evts_Hemi1_nTrks);
    readerEvts->AddVariable("mva_Evts_Hemi2_nTrks",         &mva_Evts_Hemi2_nTrks);
    readerEvts->AddVariable("mva_Evts_Hemi1_Mass",          &mva_Evts_Hemi1_Mass);
    readerEvts->AddVariable("mva_Evts_Hemi2_Mass",          &mva_Evts_Hemi2_Mass);
    readerEvts->AddVariable( "mva_HT",                      &mva_HT);
    readerEvts->AddVariable( "mva_Evts_ST",                 &mva_Evts_ST);
    readerEvts->AddVariable( "mva_Evts_njets",              &mva_Evts_njets);
    readerEvts->AddVariable( "mva_Evts_nmuon",              &mva_Evts_nmuon);
    // readerEvts->AddVariable( "mva_Evts_MediumAxes",         &mva_Evts_MediumAxes);
    // readerEvts->AddVariable( "mva_Evts_LooseAxes",          &mva_Evts_LooseAxes);
    // readerEvts->AddVariable( "mva_Evts_TightAxes",          &mva_Evts_TightAxes);
    // readerEvts->AddVariable( "mva_Evts_Mmumu",              &mva_Evts_Mmumu);
    readerEvts->AddVariable("mva_Evts_all_muon",            &mva_Evts_all_muon);
    readerEvts->BookMVA("BDTGALLBKG", weightFileEVTS_ ); // root 6.14/09, care compatiblity of versions for tmva
    readerEvts->BookMVA("BDTGDY",weightFileEVTSDY_);
    readerEvts->BookMVA("BDTGTT",weightFileEVTSTT_);

     //--------------Hemi 1--------------------

    // readerHemi1->AddVariable( "mva_Evts_jet1_pt",            &mva_Evts_jet1_pt);
    // readerHemi1->AddVariable("mva_Evts_jet1_eta",            &mva_Evts_jet1_eta);
    // readerHemi1->AddVariable("mva_Evts_Hemi1_njet_nomu",     &mva_Evts_Hemi1_njet_nomu);
    // readerHemi1->AddVariable("mva_Evts_Hemi1_pt",            &mva_Evts_Hemi1_pt);
    // readerHemi1->AddVariable("mva_Evts_Hemi1_eta",           &mva_Evts_Hemi1_eta);
    // readerHemi1->AddVariable("mva_Evts_Hemi1_phi",           &mva_Evts_Hemi1_phi);
    // readerHemi1->AddVariable("mva_Evts_Hemi1_nTrks",         &mva_Evts_Hemi1_nTrks);
    // readerHemi1->AddVariable( "mva_Evts_Hemi1_Mass",         &mva_Evts_Hemi1_Mass);
    // // readerHemi1->AddVariable( "mva_Evts_nmuon",              &mva_Evts_nmuon);
    // readerHemi1->AddVariable("mva_Evts_all_muon",            &mva_Evts_all_muon);

    // readerHemi1->BookMVA("BDTGALLBKG", weightFileHEMI1_ ); // root 6.14/09, care compatiblity of versions for tmva
    // readerHemi1->BookMVA("BDTGDY",weightFileHEMI1DY_);
    // readerHemi1->BookMVA("BDTGTT",weightFileHEMI1TT_);

      //--------------Hemi 2--------------------

    // readerHemi2->AddVariable("mva_Evts_jet2_pt",            &mva_Evts_jet2_pt);
    // readerHemi2->AddVariable("mva_Evts_jet2_eta",           &mva_Evts_jet2_eta);
    // readerHemi2->AddVariable("mva_Evts_Hemi2_njet_nomu",    &mva_Evts_Hemi2_njet_nomu);
    // readerHemi2->AddVariable("mva_Evts_Hemi2_pt",           &mva_Evts_Hemi2_pt);
    // readerHemi2->AddVariable("mva_Evts_Hemi2_eta",          &mva_Evts_Hemi2_eta);
    // readerHemi2->AddVariable("mva_Evts_Hemi2_phi",          &mva_Evts_Hemi2_phi);
    // readerHemi2->AddVariable("mva_Evts_Hemi2_nTrks",        &mva_Evts_Hemi2_nTrks);
    // readerHemi2->AddVariable("mva_Evts_Hemi2_Mass",         &mva_Evts_Hemi2_Mass);

    // readerHemi2->BookMVA("BDTGALLBKG", weightFileHEMI2_ ); // root 6.14/09, care compatiblity of versions for tmva
    // readerHemi2->BookMVA("BDTGDY",weightFileHEMI2DY_);
    // readerHemi2->BookMVA("BDTGTT",weightFileHEMI2TT_);

    //-------------Tracks---------------------
    //add the variables from my BDT (Paul)
    // reader->AddVariable( "mva_track_firstHit_x", &firsthit_X); /*!*/
    // reader->AddVariable( "mva_track_firstHit_y", &firsthit_Y); /*!*/
    // reader->AddVariable( "mva_track_firstHit_z", &firsthit_Z); /*!*/
//$$$$$$
//     reader->AddVariable( "mva_track_dxy",     &dxy);
//     reader->AddVariable( "mva_track_dz",      &dz);
//     reader->AddVariable( "mva_track_pt",      &pt );
//     reader->AddVariable( "mva_track_eta",     &eta );
//     reader->AddVariable( "mva_track_nchi2",   &NChi );
//     reader->AddVariable( "mva_track_nhits",   &nhits );
//     reader->AddVariable( "mva_ntrk10",        &ntrk10);
//     reader->AddVariable( "mva_ntrk20",        &ntrk20);
//     reader->AddVariable( "mva_ntrk30",        &ntrk30);
//     reader->AddVariable( "mva_ntrk40",        &ntrk40);
//     reader->AddVariable( "mva_dzSig",         &dzSig);
//     reader->AddVariable( "mva_drSig",         &drSig); 
//     reader->AddVariable( "mva_track_isinjet", &isinjet);
//     reader->AddVariable( "mva_track_dR",      &dR);
//     reader->AddVariable( "mva_track_dRmax",   &dRmax);
//     reader->AddVariable(" mva_track_lost", &isLost);
//$$$$$$
    // updated 240510:
//     reader->AddVariable( "mva_track_dxy",     &dxy);
//     reader->AddVariable( "mva_track_dz",      &dz);
    reader->AddVariable( "mva_track_pt",      &pt );
    reader->AddVariable( "mva_track_eta",     &eta );
//     reader->AddVariable( "mva_track_nchi2",   &NChi );
    reader->AddVariable( "mva_track_nhits",   &nhits );
    reader->AddVariable( "mva_ntrk10",        &ntrk10);
    reader->AddVariable( "mva_ntrk20",        &ntrk20);
    reader->AddVariable( "mva_ntrk30",        &ntrk30);
    reader->AddVariable( "mva_ntrk40",        &ntrk40);
    reader->AddVariable( "mva_dzSig",         &dzSig);
    reader->AddVariable( "mva_drSig",         &drSig); 
    reader->AddVariable( "mva_track_isinjet", &isinjet);
    reader->AddVariable( "mva_track_dR",      &dR);
    reader->AddVariable( "mva_track_dRmax",   &dRmax);
    reader->AddVariable(" mva_track_lost",    &isLost);
//$$$$$$

    // reader->AddVariable(" mva_track_dxyError", &dxyError); /*!*/
    // reader->AddVariable(" mva_track_dzTOpu", &dzTopu);//added on 24/03/2023 : if using bdts generated before this date, =>crash
    // reader->AddVariable(" mva_track_dzSigTOpu", &dzSigTopu);//added on 24/03/2023
    // reader->AddVariable(" mva_ValTIBHit", &TibHit);
    // reader->AddVariable(" mva_ValTOBHit", &TobHit);
    // reader->AddVariable(" mva_ValPixBarHit", &PixBarHit);
    // reader->AddVariable(" mva_nValTECHHit", &TecHit);

    reader->BookMVA( "BDTG", weightFile_ ); // root 6.14/09, care compatiblity of versions for tmva

    //-------------Vertex---------------------
    readerVtx->AddVariable( "mva_Vtx_nTrks",   &mva_V_nTrks);
    readerVtx->AddVariable( "mva_Vtx_NChi2",   &mva_V_chi);
    readerVtx->AddVariable( "mva_Vtx_step",    &mva_V_step);
    // readerVtx->AddVariable( "mva_Vtx_r",       &mva_V_r);
    // readerVtx->AddVariable( "mva_Vtx_z",       &mva_V_z);
    readerVtx->AddVariable( "mva_Vtx_MTW",     &mva_V_MTW);
    readerVtx->AddVariable( "mva_Vtx_Mass",    &mva_V_Mass);
    // readerVtx->AddVariable( "mva_Hemi_Mass",   &mva_H_Mass);
    // readerVtx->AddVariable( "mva_Vtx_dist",     &mva_V_dist);
    readerVtx->AddVariable( "mva_Vtx_ntrk10",  &mva_V_ntrk10);
    readerVtx->AddVariable( "mva_Vtx_ntrk20",  &mva_V_ntrk20);
    readerVtx->AddVariable( "mva_Vtx_MeanDCA", &mva_V_MeanDCA);

    readerVtx->BookMVA( "BDTG", weightFileVtx_ ); // root 6.14/09, care compatiblity of versions for tmva
   
    //----------------Vertex Step1----------//
    readerVtxStep1->AddVariable( "mva_Vtx_nTrks",    &mva_V_nTrks);
    readerVtxStep1->AddVariable( "mva_Vtx_NChi2",    &mva_V_chi);
    readerVtxStep1->AddVariable( "mva_Vtx_step",        &mva_V_step);
    // readerVtxStep1->AddVariable( "mva_Vtx_r",        &mva_V_r);
    // readerVtxStep1->AddVariable( "mva_Vtx_z",        &mva_V_z);
    readerVtxStep1->AddVariable( "mva_Vtx_MTW",      &mva_V_MTW);
    readerVtxStep1->AddVariable( "mva_Vtx_Mass",     &mva_V_Mass);
    // readerVtxStep1->AddVariable( "mva_Hemi_Mass",    &mva_H_Mass);
    // readerVtxStep1->AddVariable( "mva_Vtx_dist",     &mva_V_dist);
    readerVtxStep1->AddVariable( "mva_Vtx_ntrk10",   &mva_V_ntrk10);
    readerVtxStep1->AddVariable( "mva_Vtx_ntrk20",   &mva_V_ntrk20);
    readerVtxStep1->AddVariable(" mva_Vtx_MeanDCA",  &mva_V_MeanDCA);

    readerVtxStep1->BookMVA("BDTG",weightFileVtxStep1_);
}


// FlyingTopAnalyzer::~FlyingTopAnalyzer()
// {
//    // do anything here that needs to be done at destruction time
//    // (e.g. close files, deallocate resources etc.)
// }


// - Gen Particule level info about ancestor
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

  //  ---------------------------------------------------------------- //
  //  --------------- beam pipe and detector centers ----------------- //
  //  ------ tuned via a scan in (x,y), see /ui2_data1/blochd/LLTopAna/SecIntAna.C and output/h_SecInt_2017_data.root, ...
  //  ---------------------------------------------------------------- //
  float x_bmp = 0., y_bmp = 0.;
  float x_det = 0., y_det = 0.;
  if ( !isMC_ ) {
    if ( YEAR_ == 2017 ) {
      x_bmp =  0.12; y_bmp = -0.19;
      x_det =  0.12; y_det = -0.12;
    }
    if ( YEAR_ == 2018 ) {
      x_bmp =  0.18; y_bmp = -0.19;
      x_det =  0.09; y_det = -0.12;
    }
  }

  //------------------------------------
/// - Propagators init. ---------------
///------------------------------------

    PropaHitPattern* PHP = new PropaHitPattern(YEAR_);
    PropaHitPattern* NI = new PropaHitPattern(YEAR_);
    PropaHitPattern* posPHP = new PropaHitPattern(YEAR_);
    PropaHitPattern* negPHP = new PropaHitPattern(YEAR_);
  //  ---------------------------------------------------------------- //

  runNumber   = iEvent.id().run();
  eventNumber = iEvent.id().event();
  lumiBlock   = iEvent.luminosityBlock();
  
  tree_only_gen_wt=1.0;
  
  PUweight = 1;
  Prefweight = 1;

  using namespace edm;
  using namespace reco;
  using namespace pat;
  
  //LHE/gen infos

  edm::Handle<GenEventInfoProduct> genEventInfo;
  if(isMC_) iEvent.getByToken(genEventInfoToken_,genEventInfo);
  
  edm::Handle<LHEEventProduct> lheEventProduct;
  if(isMC_) iEvent.getByToken(LHEEventProductToken_,lheEventProduct);
  
  // Pruned particles are the one containing "important" stuff
  Handle<edm::View<reco::GenParticle> > pruned;
  if ( isMC_ ) iEvent.getByToken(prunedGenToken_, pruned);

  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  Handle<edm::View<pat::PackedGenParticle> > packed;
  if ( isMC_) iEvent.getByToken(packedGenToken_, packed);

  edm::Handle<edm::View<reco::GenJet>> genJets;
  if ( isMC_ ) iEvent.getByToken(genJetToken_, genJets);

  if ( isMC_ ) {
    edm::Handle< double > theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight ) ;
    double _prefiringweight =(*theprefweight);
    Prefweight=_prefiringweight;
  }
  //cout<<" _prefiringweight ="<< _prefiringweight<<endl;
  //  tree_prefir_weight= _prefiringweight;

  /*  edm::Handle< double > theprefweightup;
  if ( !runOnData_ )iEvent.getByToken(prefweightup_token, theprefweightup ) ;
  double _prefiringweightup =(*theprefweightup);
  edm::Handle< double > theprefweightdown;
  iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
  if ( !runOnData_ )double _prefiringweightdown =(*theprefweightdown);*/

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
  //PU
  // edm::Handle<PileupSummaryInfo> PU ;
  // iEvent.getByToken(puToken_,PU);
  //Photon Conversion
  edm::Handle<reco::ConversionCollection> PhotonConversion;
  iEvent.getByToken(PhotonToken_,PhotonConversion);

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotToken_, recoBeamSpotHandle);
  
  //Prescales
  // edm::Handle<pat::PackedTriggerPrescales> prescale;
  // iEvent.getByToken(PrescaleToken_,prescale);

  iSetup.get<IdealMagneticFieldRecord>().get(bField);
  const MagneticField* theMagneticField = bField.product();

  edm::ESHandle<TrackerGeometry> trackerGeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get( trackerGeomHandle );

  Handle<double> hRho;
  iEvent.getByToken(rho_token_,hRho);
  // float topweight = 1.0;
  // if ( isMC_ )
  //   {
  //   float topPtLepTrue=genEvt->leptonicDecayTop()->pt();
  //   float topPtHadTrue=genEvt->hadronicDecayTop()->pt();
  //   topweight = sqrt(SF(topPtLepTrue)*SF(topPtHadTrue));
  //   std::cout<<"topweight = "<<topweight<<std::endl;
  //   }

  hEvents->Fill(1);
  if ( isMC_ )
  {
   hEvents_with_gen_wt->Fill(1,genEventInfo->weight());
   tree_only_gen_wt=genEventInfo->weight();
  }
 
  if ( isMC_ ) 
  {
    int TruePUI = -99;
    Handle<std::vector< PileupSummaryInfo > > PupInfo;
    iEvent.getByToken(puToken_,PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
    {
      int BX = PVI->getBunchCrossing();
      if( BX == 0 ) 
      {
        TruePUI = PVI->getTrueNumInteractions();
        continue;
      }
    } // Pileup info loop ends																				    
    PUweight = lumiWeights_->weight(TruePUI);
    PU_events = float(TruePUI);

    //  cout<<" pile up weight ="<<PUweight<<endl;
  }
  tree_event_weight = 1;
  tree_event_weight = tree_only_gen_wt*Prefweight*PUweight;


  //--------------------------------------//
  //               Trigger                //
  //--------------------------------------//
  std::vector<std::string> TriggerCheck;
  const edm::Handle<edm::TriggerResults> triggerH = iEvent.getHandle(triggerResultsToken_);
  const auto triggerNames = iEvent.triggerNames(*triggerH);
  std::string TName;

  for (unsigned int i = 0; i < triggerH->size(); i++) 
  {
    TName=triggerNames.triggerName(i);
    TriggerCheck.push_back(TName);
    if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v") && !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v = false;};//std::cout<<"triggerName / prescale : "<<TName<<" / "<<prescale->getPrescaleForInder(i)<<std::endl;
    
    if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v")&& !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v = false;};//std::cout<<"triggerName / prescale : "<<TName<<" / "<<prescale->getPrescaleForInder(i)<<std::endl;    
    
    if(strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")&& !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v = false;};
    if(strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v")&& !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v = false;};

    if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v") && !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v = false;};
    if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") && !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v = false;};
    if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = false;};
    if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") &&  triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v") && !triggerH->accept(i)){HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v = false;};
    if (strstr(TName.c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") &&  triggerH->accept(i)){HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = true;} else if (strstr(TName.c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v") && !triggerH->accept(i)){HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v = false;};
    //if (strstr(TName.c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") && !triggerH->accept(i)){HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = false;};
    if (strstr(TName.c_str(),"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") &&  triggerH->accept(i)){HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v") && !triggerH->accept(i)){HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v = false;};
    if (strstr(TName.c_str(),"HLT_Ele27_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele27_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele27_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele27_WPTight_Gsf_v = false;};
    if (strstr(TName.c_str(),"HLT_Ele32_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele32_WPTight_Gsf_v = true;} else if (strstr(TName.c_str(),"HLT_Ele32_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele32_WPTight_Gsf_v = false;};
    if (strstr(TName.c_str(),"HLT_Ele35_WPTight_Gsf_v") &&  triggerH->accept(i)){HLT_Ele35_WPTight_Gsf_v = true;} else if (strstr(TName.c_str()," HLT_Ele35_WPTight_Gsf_v") && !triggerH->accept(i)){HLT_Ele35_WPTight_Gsf_v = false;};
    if (strstr(TName.c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") &&  triggerH->accept(i)){HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = true;} else if (strstr(TName.c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v") && !triggerH->accept(i)){HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v = false;};
    if (strstr(TName.c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v") &&  triggerH->accept(i)){HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v = true;} else if (strstr(TName.c_str(),"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v") && !triggerH->accept(i)){HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v = false;};
    if (strstr(TName.c_str(),"HLT_PFMET120_PFMHT120_IDTight_PFHT60_v") &&  triggerH->accept(i)){HLT_PFMET120_PFMHT120_IDTight_PFHT60_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET120_PFMHT120_IDTight_PFHT60_v") && !triggerH->accept(i)){HLT_PFMET120_PFMHT120_IDTight_PFHT60_v = false;};
    if (strstr(TName.c_str(),"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v") &&  triggerH->accept(i)){HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v") && !triggerH->accept(i)){HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v = false;};
    if (strstr(TName.c_str(),"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") &&  triggerH->accept(i)){HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v") && !triggerH->accept(i)){HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v = false;};
    if (strstr(TName.c_str(),"HLT_PFMET250_HBHECleaned_v") &&  triggerH->accept(i)){HLT_PFMET250_HBHECleaned_v = true;} else if (strstr(TName.c_str(),"HLT_PFMET250_HBHECleaned_v") && !triggerH->accept(i)){HLT_PFMET250_HBHECleaned_v = false;};
    if (strstr(TName.c_str(),"HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v") &&  triggerH->accept(i)){HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v = true;} else if (strstr(TName.c_str(),"HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v") && !triggerH->accept(i)){HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v = false;};
    if (strstr(TName.c_str(),"HLT_IsoMu24_v") && triggerH->accept(i)){HLT_IsoMu24_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu24_v") &&!triggerH->accept(i)){HLT_IsoMu24_v = false;};
    if (strstr(TName.c_str(),"HLT_IsoTkMu24_v") && triggerH->accept(i)){HLT_IsoTkMu24_v = true;} else if (strstr(TName.c_str(),"HLT_IsoTkMu24_v") &&!triggerH->accept(i)){HLT_IsoTkMu24_v = false;};
    if (strstr(TName.c_str(),"HLT_IsoMu27_v") && triggerH->accept(i)){HLT_IsoMu27_v = true;} else if (strstr(TName.c_str(),"HLT_IsoMu27_v") &&!triggerH->accept(i)){HLT_IsoMu27_v = false;};
  }


  //////////////////////////////////
  //////////////////////////////////
  //////////    BS     /////////////
  //////////////////////////////////
  //////////////////////////////////
  
  BeamSpot const & bs = *recoBeamSpotHandle;
  tree_bs_PosX = bs.x0();
  tree_bs_PosY = bs.y0();
  tree_bs_PosZ = bs.z0();


  //////////////////////////////////
  //////////////////////////////////
  //////    Primary Vertex   ///////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_nPV = primaryVertex->size();
  if ( !primaryVertex->empty() ) {
    tree_PV_x     = (*primaryVertex)[0].x(); // l'index 0 donne le PV!
    tree_PV_y     = (*primaryVertex)[0].y();
    tree_PV_z     = (*primaryVertex)[0].z();
    tree_PV_ez    = (*primaryVertex)[0].zError();
    tree_PV_NChi2 = (*primaryVertex)[0].normalizedChi2();
    tree_PV_ndf   = (*primaryVertex)[0].ndof();
    tree_PV_rho   = (*primaryVertex)[0].position().Rho();
  }
  const reco::Vertex &PV = primaryVertex->front();

  //T_Rho = *hRho;
  const double rho_val = *(hRho.product());                                                                                                                                                                  
  Rho = rho_val; // note that this has nothing to do with tree_PV_rho !

  bool tree_Good_PV = false;
  if ( tree_PV_ndf > 4 && abs(tree_PV_z) < 24. && abs(tree_PV_rho) < 2. ) 
    tree_Good_PV = true;


  //////////////////////////////////
  //////////////////////////////////
  ///////////   Muons   ////////////
  //////////////////////////////////
  //////////////////////////////////

  int isMatched = 0;
  if ( isMC_ ) 
  {
    for (const pat::Muon &mu : *muons)
    {
      for (size_t i=0; i<pruned->size(); i++)
      {
        const GenParticle & genIt = (*pruned)[i];
                // just to keep the smuon and neutralino masses in case the event is not filtered
        if ( abs(genIt.pdgId()) != 13 ) continue;
        if ( genIt.charge() * mu.charge() < 0. ) continue;
        float Gen_pt  = genIt.pt();
        float Gen_eta = genIt.eta();
        float Gen_phi = genIt.phi();
        if ( abs(Gen_eta) > 2.5 ) continue;
        float dpt  = abs( Gen_pt / mu.pt() - 1. );
        float deta = abs( Gen_eta - mu.eta() );
        float dphi = abs( Deltaphi( Gen_phi, mu.phi() ) );
        if ( abs(deta) < 0.1 && abs(dphi) < 0.1 && abs(dpt) < 0.1 ) {
            if (mu.triggered("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*") )
              {
                isMatched += 10;
              }
            if ( mu.triggered("HLT_IsoMu24_v*"))
              {
                isMatched += 1;
              }
        } // muon matching					
        if ( isMatched > 0 ) break;																	      
      } //  end loop on generated muons 																						      
    }
  } // endif MC 

  // std::cout<<"isMatched "<<isMatched<<std::endl;
  tree_muon_GenRecoTriggerMatched = isMatched;
  int nmu = 0;
  int allnmu = 0;
  tree_nmu = 0;
  tree_all_nmu = 0;
  tree_LT  = 0.;

  for (const pat::Muon &mu : *muons)
  {
    //------------Muon Corrections --------------//
    if ( mu.pt() < 3. ) continue;
    if ( abs(mu.eta()) > 2.4 ) continue;  // muon acceptance
    if ( !mu.isGlobalMuon() ) continue;
    bool isPromptMuon = false;
    // reco::Muon::Selector::TkIsoTight
    if ( abs(mu.muonBestTrack()->dxy(PV.position())) < 0.1 &&
         abs(mu.muonBestTrack()->dz(PV.position()))  < 0.2    ) //mu.isTightMuon(PV) &&
                                                                //  mu.passed(reco::Muon::Selector::MiniIsoTight) &&
      isPromptMuon = true;
    // if ( !mu.isMediumMuon() && !isPromptMuon ) continue;

    float correction = 1.;
    float  smearedPt = 0.;
    int isGen = 0;
    if ( !isMC_ ) 
      correction = rc.kScaleDT( mu.charge(), mu.pt(), mu.eta(), mu.phi(), 0, 0 );
    if ( isMC_ ) 
    {
      for (size_t i=0; i<pruned->size(); i++)
      {
        const GenParticle & genIt = (*pruned)[i];
                // just to keep the smuon and neutralino masses in case the event is not filtered
        if ( abs(genIt.pdgId()) == 1000013 ) tree_smu_mass = genIt.mass();
        if ( abs(genIt.pdgId()) == 1000023 ) tree_neu_mass = genIt.mass();
        if ( abs(genIt.pdgId()) != 13 ) continue;
        if ( genIt.charge() * mu.charge() < 0. ) continue;
        float Gen_pt  = genIt.pt();
        float Gen_eta = genIt.eta();
        float Gen_phi = genIt.phi();
        if ( Gen_pt < 10. ) continue;
        if ( abs(Gen_eta) > 2.5 ) continue;
        float dpt  = abs( Gen_pt / mu.pt() - 1. );
        float deta = abs( Gen_eta - mu.eta() );
        float dphi = abs( Deltaphi( Gen_phi, mu.phi() ) );
        if ( abs(deta) < 0.1 && abs(dphi) < 0.1 && abs(dpt) < 0.1 ) {
          correction = rc.kSpreadMC( mu.charge(), mu.pt(), mu.eta(), mu.phi(), Gen_pt, 0, 0 );
          isGen = 1;
          int motherPdgId = 0;
          const Candidate * mom   = genIt.mother();
          if ( mom ) motherPdgId = mom->pdgId();
          if ( motherPdgId == 23 )      isGen = 23; // Z
          if ( abs(motherPdgId) == 24 ) isGen = 24; // W
          if ( isGen == 13 && motherPdgId != 0 ) isGen = abs(motherPdgId);
        } // muon matching																						      
        if ( isGen != 0 ) break;
      } //  end loop on generated muons 																						      
      if  ( isGen == 0 ) correction = rc.kSmearMC(mu.charge(), mu.pt(), mu.eta(), mu.phi(), mu.innerTrack()->hitPattern().trackerLayersWithMeasurement(), gRandom->Rndm(), 0, 0);
    } // endif MC 
    																					    
    smearedPt = mu.pt() * correction;
    if ( smearedPt < 10. ) continue;
    if (showlog) std::cout<<"smeared pt vs direct pt "<<smearedPt<<" vs "<<mu.pt()<<std::endl;
    																			  
    tree_muon_isPrompt.push_back( isPromptMuon );
    tree_muon_correction.push_back( correction );
    tree_muon_gen.push_back(      isGen );
    tree_muon_pt.push_back(       smearedPt);

    // test muon scale factor Paul
    // MuonEps3File_ (iConfig.getParameter<std::string>( "muoneps3file" ) ),
    // MuonEps3Path_

    // bins
    //   "pt": [15, 20, 25, 30, 40, 50, 60, 120],
    // "abseta": [0, 0.9, 1.2, 2.1, 2.4],

    // //-------------Epsilon 1 trigger/(id and iso) ----------//
    // TFile file1(MuonEps1File_.c_str(), "READ");
    // if (!file1.IsOpen()) {
    //   std::cout << "Failed to open muon eps1 SF file: " << MuonEps1File_;
    //   return;
    // }
    // // Get the histogram
    // TH2F* histogram1 = dynamic_cast<TH2F*>(file1.Get(MuonEps1Path_.c_str()));
    // if (!histogram1) {
    //   std::cout << "Failed to retrieve muon eps3 histogram: " << MuonEps1Path_;
    //   file1.Close();
    //   return;
    // }
    //  //-----------END--Epsilon 1 trigger/(id and iso) ----------//

    // //-------------Epsilon 2 trigger/(id and iso) ----------//
    // TFile file2(MuonEps2File_.c_str(), "READ");
    // if (!file2.IsOpen()) {
    //   std::cout << "Failed to open muon eps2 SF file: " << MuonEps2File_;
    //   return;
    // }
    // // Get the histogram
    // TH2F* histogram2 = dynamic_cast<TH2F*>(file2.Get(MuonEps2Path_.c_str()));
    // if (!histogram2) {
    //   std::cout << "Failed to retrieve muon eps2 histogram2: " << MuonEps2Path_;
    //   file2.Close();
    //   return;
    // }
    //  //-----------END--Epsilon 2 trigger/(id and iso) ----------//


    // //-------------Epsilon 3 trigger/(id and iso) ----------//
    // TFile file3(MuonEps3File_.c_str(), "READ");
    // if (!file3.IsOpen()) {
    //   std::cout << "Failed to open muon eps3 SF file: " << MuonEps3File_;
    //   return;
    // }
    // // Get the histogram
    // TH2F* histogram3 = dynamic_cast<TH2F*>(file3.Get(MuonEps3Path_.c_str()));
    // if (!histogram3) {
    //   std::cout << "Failed to retrieve muon eps3 histogram3: " << MuonEps3Path_;
    //   file3.Close();
    //   return;
    // }
     //-----------END--Epsilon 3 trigger/(id and iso) ----------//

    // // Access data from the histogram and perform analysis
    // int nbinsX = histogram1->GetNbinsX();
    // int nbinsY = histogram1->GetNbinsY();

    float SF = 1;

    // for (int i = 1; i <= nbinsX; ++i) {
    //   for (int j = 1; j <= nbinsY; ++j) {
    //     // double binContent = histogram1->GetBinContent(i, j);
    //     double binlowedgeX = histogram1->GetXaxis()->GetBinLowEdge(i);
    //     double binWidthX = histogram1->GetXaxis()->GetBinWidth(i);
    //     double binlowedgeY = histogram1->GetYaxis()->GetBinLowEdge(j);
    //     double binWidthY = histogram1->GetYaxis()->GetBinWidth(j);
    //     // double binError = histogram1->GetBinError(i, j);
    //     // std::cout<<smearedPt<<" with abs(mu.eta()) : "<<abs(mu.eta())<<std::endl;
    //     // std::cout<<" binlowedgeX and binlowedgeX+binWidthX : "<<binlowedgeX<<" and "<<binlowedgeX+binWidthX<<std::endl;
    //     // std::cout<<" binlowedgeY and binlowedgeY+binWidthY : "<<binlowedgeY<<" and "<<binlowedgeY+binWidthY<<std::endl;
    //     if ( smearedPt >= binlowedgeY && smearedPt <= (binlowedgeY+binWidthY) && abs(mu.eta()) >= binlowedgeX && abs(mu.eta()) <= (binlowedgeX+binWidthX) )
    //       {
    //         SF = histogram1->GetBinContent(i, j)*histogram2->GetBinContent(i, j)*histogram3->GetBinContent(i, j);
    //         // std::cout<<" SF : "<<SF<<std::endl;
    //         // break;
    //       }
    //     // std::cout<<" bin low edge X with binwidth: "<<binlowedgeX<<"  width :"<<binWidthX<<" and Y : "<<binlowedgeY<<" with binWidth : "<<binWidthY<<std::endl;
    //     // Perform your analysis using binContent and binError
    //     // Example: Print bin content and error
    //     // std::cout << "Bin (" << i << ", " << j << "): Content = " << binContent << ", Error = " << binError << std::endl;
    //   }
    // }
    // file1.Close();
    // file2.Close();
    // file3.Close();
    
    tree_muon_SF.push_back(SF);
    tree_muon_eta.push_back(      mu.eta());
    tree_muon_phi.push_back(      mu.phi());
    tree_muon_x.push_back(        mu.vx());
    tree_muon_y.push_back(        mu.vy());
    tree_muon_z.push_back(        mu.vz());
    tree_muon_dxy.push_back(	    mu.muonBestTrack()->dxy(PV.position()));
    tree_muon_dxyError.push_back( mu.muonBestTrack()->dxyError());
    tree_muon_dz.push_back(       mu.muonBestTrack()->dz(PV.position()));
    tree_muon_dzError.push_back(  mu.muonBestTrack()->dzError());
    tree_muon_charge.push_back(   mu.charge());
    tree_muon_isLoose.push_back(  mu.isLooseMuon());
    tree_muon_isMedium.push_back( mu.isMediumMuon());
    tree_muon_isTight.push_back(  mu.isTightMuon(PV));
    tree_muon_isGlobal.push_back( mu.isGlobalMuon());
    tree_muon_isoR3.push_back(    mu.trackIso());
    tree_muon_trigger_dimu.push_back(  mu.triggered("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*"));  
    tree_muon_trigger_isomu.push_back( mu.triggered("HLT_IsoMu24_v*"));
    // reco::TrackRef MuonInnerRef = mu.innerTrack();
    // https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/DataFormats/MuonReco/interface/Muon.h
    tree_muon_PFIsoVeryLoose.push_back(mu.passed(reco::Muon::Selector::PFIsoVeryLoose));
    tree_muon_PFIsoLoose.push_back( mu.passed(reco::Muon::Selector::PFIsoLoose));
    tree_muon_PFIsoMedium.push_back(mu.passed(reco::Muon::Selector::PFIsoMedium));
    tree_muon_PFIsoTight.push_back( mu.passed(reco::Muon::Selector::PFIsoTight));
    tree_muon_TkIsoLoose.push_back( mu.passed(reco::Muon::Selector::TkIsoLoose));
    tree_muon_TkIsoTight.push_back( mu.passed(reco::Muon::Selector::TkIsoTight));
    tree_muon_MiniIsoLoose.push_back( mu.passed(reco::Muon::Selector::MiniIsoLoose));  //Wrongly stored
    tree_muon_MiniIsoMedium.push_back(mu.passed(reco::Muon::Selector::MiniIsoMedium)); //Wrongly stored
    tree_muon_MiniIsoTight.push_back( mu.passed(reco::Muon::Selector::MiniIsoTight));  //Wrongly stored

    double Aeff_Fall17[5] = { 0.0566, 0.0562, 0.0363, 0.0119, 0.0064 };
    double EA;
    auto iso = mu.miniPFIsolation();
    auto chg = iso.chargedHadronIso();
    auto neu = iso.neutralHadronIso();
    auto pho = iso.photonIso();
    if( TMath::Abs(mu.eta()) < 0.8 ) EA = Aeff_Fall17[0];
    else if( TMath::Abs(mu.eta()) < 1.3 ) EA = Aeff_Fall17[1];
    else if( TMath::Abs(mu.eta()) < 2.0 ) EA = Aeff_Fall17[2];
    else if( TMath::Abs(mu.eta()) < 2.2 ) EA = Aeff_Fall17[3];
    else EA = Aeff_Fall17[4];
    float R = 10.0 / std::min( std::max( mu.pt(), 50.0 ), 200.0 );
    EA *= std::pow( R / 0.3, 2 );
    float miniIso = ( chg + TMath::Max( 0.0, neu + pho - (Rho) * EA ) ) / mu.pt();//smaeredPT ???
    tree_muon_miniIso.push_back(miniIso);
    if ( mu.isGlobalMuon() ) {
      //cout<<" muon inner hits="<<mu.innerTrack()->hitPattern().trackerLayersWithMeasurement()<<endl;
      tree_muon_trkLayers.push_back(mu.innerTrack()->hitPattern().trackerLayersWithMeasurement());
    }
    
    allnmu++;
    if ( isPromptMuon ) nmu++;
    tree_LT    += smearedPt;
  } // end loop on muons
  tree_nmu = nmu;
  tree_all_nmu = allnmu;

  //----------------------------- SORTED MUONS ----------------------------------------//
  // Please note that after the loop on muons that have been corrected, they may not be ordered by decreasing value of pt
  // => needs sorting and therefore, accessing the information about muons has to be changed => index_muon[k] instead of k where k is the iterator basically
  int size_muon = 0;
  for (unsigned int counter_muon = 0; counter_muon < tree_muon_pt.size(); counter_muon++) {
    muon_pt[counter_muon] = tree_muon_pt[counter_muon];size_muon++;
  }
  if ( tree_muon_pt.size() < 10000 ) { 
    TMath::Sort(size_muon, muon_pt, index_muon);
  }


  //////////////////////////////////
  //////////////////////////////////
  ////////   Electrons   ///////////
  //////////////////////////////////
  //////////////////////////////////
  // One could consider the emu channel, or even ee channel
 
  int nEl = 0;
  tree_all_nel = 0;
  for (const pat::Electron &el: *electrons)
  {
    if ( el.pt() < 10. ) continue;
    if ( abs(el.eta()) > 2.4 || (abs(el.eta()) > 1.442 && abs(el.eta()) < 1.556)) continue;
    bool isPromptElec = false;
    if (  // for 2018 el.electronID("cutBasedElectronID-Fall17-94X-V2-tight") &&
    (( abs(el.eta()) <= 1.479 && abs(el.gsfTrack()->dxy(PV.position())) < 0.05 && abs(el.gsfTrack()->dz(PV.position())) < 0.10 ) || 
     ( abs(el.eta()) >  1.556 && abs(el.gsfTrack()->dxy(PV.position())) < 0.10 && abs(el.gsfTrack()->dz(PV.position())) < 0.20 )   ) )
      isPromptElec = true; 
    // if ( !el.electronID("cutBasedElectronID-Fall17-94X-V2-medium") && !isPromptElec ) continue;

    float correction = 1;
    TLorentzVector elcor(0.,0.,0.,0);
    TLorentzVector el1(0.,0.,0.,0.);
    el1.SetPxPyPzE(el.px(), el.py(),el.pz(),el.energy());
    tree_electron_ecal_trk_postcorr.push_back(el.userFloat("ecalTrkEnergyPostCorr"));
    if ( el.energy() != 0 )
      correction=el.userFloat("ecalTrkEnergyPostCorr")/el.energy();
    // cout<< " ele pt before CORREC***** = "<<el1.Pt()<<endl;    
    elcor = el1 * correction;
    // cout<< " ele pt AFTER CORRECT = "<<elcor.Pt()<<endl;

    int isGen = 0;
    if ( isMC_ ) 
    {
      for (size_t i=0; i<pruned->size(); i++)
      {
        const GenParticle & genIt = (*pruned)[i];
      if ( abs(genIt.pdgId()) != 11 ) continue;
      if ( genIt.charge() * el.charge() < 0. ) continue;
        float Gen_pt  = genIt.pt();
        float Gen_eta = genIt.eta();
        float Gen_phi = genIt.phi();
      if ( Gen_pt < 10. ) continue;
      if ( abs(Gen_eta) > 2.5 ) continue;
	float dpt  = abs( Gen_pt / elcor.Pt() - 1. );
	float deta = abs( Gen_eta - elcor.Eta() );
  	float dphi = abs( Deltaphi( Gen_phi, elcor.Phi() ) );
        if ( deta < 0.1 && dphi < 0.1 && dpt < 0.1 ) {
          isGen = 1;
          int motherPdgId = 0;
          const Candidate * mom   = genIt.mother();
	  if ( mom ) motherPdgId = mom->pdgId();
	  if ( motherPdgId == 23 )      isGen = 23; // Z
	  if ( abs(motherPdgId) == 24 ) isGen = 24; // W
	  if ( isGen == 1 && motherPdgId != 0 ) isGen = abs(motherPdgId);
        } // electron matching																						      
	if ( isGen != 0 ) break;
      } //  end loop on generated electrons 																						      
    } // endif MC 																					      
    tree_electron_gen.push_back( isGen );

    tree_electron_isPrompt.push_back( isPromptElec );
    tree_electron_IsLoose.push_back(  el.electronID("cutBasedElectronID-Fall17-94X-V2-loose"));  // valid for Run2
    tree_electron_IsMedium.push_back( el.electronID("cutBasedElectronID-Fall17-94X-V2-medium")); //  valid for Run2
    tree_electron_IsTight.push_back(  el.electronID("cutBasedElectronID-Fall17-94X-V2-tight"));  //  valid for Run2																	   
    tree_electron_pt.push_back(     elcor.Pt());
    tree_electron_eta.push_back(    elcor.Eta());//cluster eta
    tree_electron_phi.push_back(    elcor.Phi());
    tree_electron_x.push_back(      el.vx());
    tree_electron_y.push_back(      el.vy());
    tree_electron_z.push_back(      el.vz());
    tree_electron_energy.push_back( elcor.E());
    tree_electron_et.push_back(     el.et()); //TLorentzvector
    tree_electron_charge.push_back( el.charge());
    tree_electron_isoR4.push_back(  el.trackIso());//returns the value of the summed track pt in a cone of deltaR<0.4
    tree_electron_dxy.push_back(    el.gsfTrack()->dxy(PV.position()));
    tree_electron_dz.push_back(     el.gsfTrack()->dz(PV.position()));
    // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_94X_and_later
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCategoryBasedElectronID
    // HLT Ele23Ele12 noDZ v*
    // HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
        //       	barrel 	endcap
    // d0, cm 	0.05 	0.10
    // dz, cm 	0.10 	0.20 

    tree_all_nel++;
    if ( isPromptElec ) nEl++;
  } // end loop on electrons
  tree_electron_nEle = nEl;

  //-----------------------------Sorted ELECTRONS -------------------------------//
  // Please note that after the loop on electrons that have been corrected, they mya not be ordered by decreasing value of pt
  // => needs sortings and therefore, accessing the information about electrons has to be changed => index_el[k] instead of k where k is the iterator basically
  int size_el = 0;
  for (unsigned int counter_el = 0; counter_el < tree_electron_pt.size(); counter_el++) {
    el_pt[counter_el] = tree_electron_pt[counter_el];
    size_el++;
  }
  if ( tree_electron_pt.size() < 10000 ) { 
    TMath::Sort(size_el, el_pt, index_el);
  }
  

  //////////////////////////////////
  //////////////////////////////////
  /////   DiLepton Selection   /////
  //////////////////////////////////
  //////////////////////////////////

  int imu1 = -1, imu2 = -1;
  int imu1_SS = -1, imu2_SS = -1;
  // int Q1 = 0, Q2 = 0;
  //float mupt1, mueta1, muphi1, mupt2, mueta2, muphi2;
  float mu_mass = 0.1057, el_mass = 0.0005;
  float MuonMasses[2] = {mu_mass,mu_mass};//For muon channel
  float MuEMasses[2] = {mu_mass,el_mass};//For EMu channel
  float EMasses[2] = {el_mass,el_mass};//For electron channel
  TLorentzVector v1, v2, v;
  tree_Mmumu = 0.;
  tree_MmumuSameSign = 0.;
  
  //------ Dimuon Channel ------//
  if ( nmu >= 2 && MuonChannel ) 
  {
    // Find the function in ../interface/Filter.h
    std::vector<float> DiLeptonData = DiLeptonMass(AllowDiLeptonSameSign,nmu,
    tree_all_nmu,tree_muon_isTight,tree_muon_MiniIsoTight,tree_muon_isPrompt,
    tree_muon_pt,tree_muon_eta,tree_muon_phi,tree_muon_charge,index_muon,MuonMasses);
    tree_Mmumu = DiLeptonData[0];
    tree_MmumuSameSign = DiLeptonData[1];
    imu1 = DiLeptonData[2];
    imu2 = DiLeptonData[3];
    imu1_SS = DiLeptonData[4];
    imu2_SS = DiLeptonData[5];
  }

  //---------- Dielectron Channel ----------//
  if ( nEl >= 2 && ElChannel ) 
  {
    // Find the function in ../interface/Filter.h
    std::vector<float> DiLeptonData = DiLeptonMass(AllowDiLeptonSameSign,nEl,
    tree_all_nel,tree_electron_IsTight,tree_electron_IsTight,tree_electron_isPrompt,
    tree_electron_pt,tree_electron_eta,tree_electron_phi,tree_electron_charge,index_el,EMasses);//The Iso and ID flags are contained inside the electron ID
    tree_Mmumu = DiLeptonData[0];
    tree_MmumuSameSign = DiLeptonData[1];
    imu1 = DiLeptonData[2];
    imu2 = DiLeptonData[3];
    imu1_SS = DiLeptonData[4];
    imu2_SS = DiLeptonData[5];
  }

  //----------   EleMu Channel  ----------//
  bool LeadingMuon = true;
  if( nmu >=1 && nEl >=1 && EMuChannel )
  {
    // Find the function in ../interface/Filter.h
    std::vector<float> DiLeptonData = EMuMass(AllowDiLeptonSameSign,
    nmu,tree_all_nmu,tree_muon_isTight,tree_muon_MiniIsoLoose,tree_muon_isPrompt,
    tree_muon_pt,tree_muon_eta,tree_muon_phi,tree_muon_charge,index_muon,
     nEl,tree_all_nel,tree_electron_IsTight,tree_electron_IsTight,tree_electron_isPrompt,
    tree_electron_pt,tree_electron_eta,tree_electron_phi,tree_electron_charge,index_el,
    MuEMasses);

    tree_Mmumu = DiLeptonData[0];
    tree_MmumuSameSign = DiLeptonData[1];
    imu1 = DiLeptonData[2];
    imu2 = DiLeptonData[3];
    imu1_SS = DiLeptonData[4];
    imu2_SS = DiLeptonData[5];
  }// end of Emu Channel

  //---------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------/////
  // After that for mumu and ee chanels: call the leading lepton by imu1 index and the sub leading lepton by imu2 index
  // but for the emu channel: imu1 is the muon and imu2 is the electron
  //---------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------/////


  //////////////////////////////////
  //////////////////////////////////
  ////////   FILTER CHECK  /////////
  //////////////////////////////////
  //////////////////////////////////
  tree_only_tigger_filter = false;
  tree_Filter = false;
  tree_FilterSameSign = false;
  tree_nTracks = 0;
  tree_nLostTracks = 0;
  tree_TRACK_SIZE = 0;
  tree_nSecInt = 0;
  tree_nV0_reco = 0;

  // Mu trigger 2018 : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v || HLT_IsoMu24_v
  // Mu trigger 2017 : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v || HLT_IsoMu27_v
  // Mu trigger 2016 : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL(_DZ)_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL(_DZ)_v || HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_(DZ)_v  || HLT_IsoMu24_v || HLT_IsoTkMu24_v 

  // bool HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v;    // USED in 2016-2018                                                                                                             
  // bool HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v; // USED in 2016-2018  


  if (YEAR_ == 2018)
    {
      if( ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v || HLT_IsoMu24_v ) && MuonChannel ) tree_only_tigger_filter = true; // for all data and MC                       
      if( ( HLT_Ele32_WPTight_Gsf_v || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v) && ElChannel ) tree_only_tigger_filter = true;//HLT_DoubleEle25_CaloIdL_MW                
      if( ( HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v ||HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v || HLT_IsoMu24_v || HLT_Ele32_WPTight_Gsf_v) && EMuChannel ) tree_only_tigger_filter = true;
    }
  else if (YEAR_ == 2017)
    {
      if((HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v ||  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v||  HLT_IsoMu27_v ) && MuonChannel ) tree_only_tigger_filter = true;// only for MC      
      //if(( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v ||  HLT_IsoMu27_v ) && MuonChannel ) tree_only_tigger_filter = true;//only for data runs B                              
      // if(( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v ||  HLT_IsoMu27_v ) && MuonChannel ) tree_only_tigger_filter = true;// only for data run C-F                                                                                                                                                  
      if ( (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v || HLT_Ele35_WPTight_Gsf_v) && ElChannel ) tree_only_tigger_filter = true;//HLT_DoubleEle33_CaloIdL_MW              
      if ( (HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v ||HLT_Ele35_WPTight_Gsf_v || HLT_IsoMu27_v) && EMuChannel )  tree_only_tigger_filter = true;//HLT_Ele35_WPTight_Gsf, //HLT_IsoMu27                                                                                                 
    }
  else {
    if(( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v) && MuonChannel ) tree_only_tigger_filter = true;// for MC                                                                                    
                                                                                                                                                                            
    //if((HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v) && MuonChannel ) tree_only_tigger_filter = true;// \only for data B-G          
    //if((  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v) && MuonChannel ) tree_only_tigger_filter = true;// only for data run H  
    if ( (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Ele27_WPTight_Gsf_v)  && ElChannel ) tree_only_tigger_filter = true; //DoubleEle33_CaloIdL_MW, DoubleEle33_CaloIdL_GsfTrkIdVL               
    if ( ( HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v || HLT_Ele27_WPTight_Gsf_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v ) && EMuChannel ) tree_only_tigger_filter = true;  // only for B_G and MC , and HLT_Ele27_WPTight_Gsf, HLT_IsoMu24, HLT_IsoTkMu24 for all data & MC                         
    //if ( (HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Ele27_WPTight_Gsf_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v) && EMuChannel )  tree_only_tigger_filter = true;        // only for H data                                                                                          
  }//else 

if (YEAR_ == 2018)
  {
    if ( ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v || HLT_IsoMu24_v  ) // for all data and MC
        && nmu >= 2 && MuonChannel ) {
      if ( tree_Mmumu > 10. ) tree_Filter = true; 
      if ( AllowDiLeptonSameSign && tree_MmumuSameSign > 10. ) tree_FilterSameSign = true;
    }

    if ( (HLT_Ele32_WPTight_Gsf_v || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v)//HLT_DoubleEle25_CaloIdL_MW  
        && nEl >= 2 && ElChannel ) {
      if ( tree_Mmumu > 10. ) tree_Filter = true; 
      if ( AllowDiLeptonSameSign && tree_MmumuSameSign > 10. ) tree_FilterSameSign = true;
    }
    
    if ( ( HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v || HLT_IsoMu24_v || HLT_Ele32_WPTight_Gsf_v)&& nmu >= 1 && nEl >= 1 && EMuChannel ) {
      if ( tree_Mmumu > 10. ) tree_Filter = true; 
      if ( AllowDiLeptonSameSign && tree_MmumuSameSign > 10. ) tree_FilterSameSign = true;
    }
  }





else if (YEAR_ == 2017)
  {
    if ( ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v ||  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v||  HLT_IsoMu27_v ) 
    	 //if(( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v ||  HLT_IsoMu27_v )//only for data runs B                                                                             
	 // if(( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v ||  HLT_IsoMu27_v )// only for data run C-F  
        && nmu >= 2 && MuonChannel ) {
      if ( tree_Mmumu > 10. ) tree_Filter = true; 
      if ( AllowDiLeptonSameSign && tree_MmumuSameSign > 10. ) tree_FilterSameSign = true;
    }

    if ( (HLT_Ele35_WPTight_Gsf_v || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v)
        && nEl >= 2 && ElChannel ) {
      if ( tree_Mmumu > 10. ) tree_Filter = true; 
      if ( AllowDiLeptonSameSign && tree_MmumuSameSign > 10. ) tree_FilterSameSign = true;
    }
    
    if ( ( HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v ) 
        && nmu >= 1 && nEl >= 1 && EMuChannel ) {
      if ( tree_Mmumu > 10. ) tree_Filter = true; 
      if ( AllowDiLeptonSameSign && tree_MmumuSameSign > 10. ) tree_FilterSameSign = true;
    }
  }






else if (YEAR_ == 2016)
  {
    if ( ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v)//for MC 
             //if((HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v) && nmu >= 2 && MuonChannel ) {//only for data B-G    //if((HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v) && nmu >= 2 && MuonChannel) {//only for data  H                
        && nmu >= 2 && MuonChannel ) {
          if ( tree_Mmumu > 10. ) tree_Filter = true; 
          if ( AllowDiLeptonSameSign && tree_MmumuSameSign > 10. ) tree_FilterSameSign = true;
    }
    
    if ( (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Ele27_WPTight_Gsf_v)&& nEl >= 2 && ElChannel ) { //DoubleEle33_CaloIdL_MW, DoubleEle33_CaloIdL_GsfTrkIdVL     
      if ( tree_Mmumu > 10. ) tree_Filter = true;
      if ( AllowDiLeptonSameSign && tree_MmumuSameSign > 10. ) tree_FilterSameSign = true;
    }
    if ( ( HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v || HLT_Ele27_WPTight_Gsf_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v ) && nmu >= 1 && nEl >= 1 && EMuChannel ) {// only for B_G and MC                             
      // if ( (HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v || HLT_Ele27_WPTight_Gsf_v || HLT_IsoMu24_v || HLT_IsoTkMu24_v) && nmu >= 1 && nEl >= 1 && EMuChannel ) { // only for H data                                                                                                 
      if ( tree_Mmumu > 10. ) tree_Filter = true;
      if ( AllowDiLeptonSameSign && tree_MmumuSameSign > 10. ) tree_FilterSameSign = true;
    }
  }//else       

  if ( !tree_Good_PV ) {
    tree_Filter = false;
    tree_FilterSameSign = false;
  }


  //////////////////////////////////
  //////////////////////////////////
  /////////    FILTER    ///////////
  //////////////////////////////////
  //////////////////////////////////

// //$$
//  if ( tree_Filter || tree_FilterSameSign ) {
// //$$


  ///////////////////////////////////////////////
  /////////    lepton informations    ///////////
  ///////////////////////////////////////////////

  if ( !tree_Filter && tree_FilterSameSign ) {
    imu1 = imu1_SS;
    imu2 = imu2_SS;
  }

  TLorentzVector Vlep1, Vlep2, vll;
  float lep1_pt=0, lep2_pt=0, lep1_eta=0, lep2_eta=0, lep1_phi=0, lep2_phi=0;
  int lep1_Q =0; int lep2_Q =0;
  float lep1_mass=-1, lep2_mass=-1;

  if ( imu1 >= 0 && imu2 >= 0 ) {
  if ( MuonChannel ) {
    lep1_pt  = tree_muon_pt[imu1];
    lep2_pt  = tree_muon_pt[imu2];
    lep1_eta = tree_muon_eta[imu1];
    lep2_eta = tree_muon_eta[imu2];
    lep1_phi = tree_muon_phi[imu1];
    lep2_phi = tree_muon_phi[imu2];
    lep1_Q   = tree_muon_charge[imu1];
    lep2_Q   = tree_muon_charge[imu2];
    lep1_mass = mu_mass;
    lep2_mass = mu_mass;
  }
  if ( ElChannel ) {
    lep1_pt  = tree_electron_pt[imu1];
    lep2_pt  = tree_electron_pt[imu2];
    lep1_eta = tree_electron_eta[imu1];
    lep2_eta = tree_electron_eta[imu2];
    lep1_phi = tree_electron_phi[imu1];
    lep2_phi = tree_electron_phi[imu2];
    lep1_Q   = tree_electron_charge[imu1];
    lep2_Q   = tree_electron_charge[imu2];
    lep1_mass = el_mass;
    lep2_mass = el_mass;
  }
  if ( EMuChannel ) {
          if(tree_Good_PV){
	tree_reco_muon_leadingpt.push_back(tree_muon_pt[imu1] );
	tree_reco_electron_leadingpt2.push_back(tree_electron_pt[imu2]);
	tree_reco_muon_leadingeta.push_back(tree_muon_eta[imu1]);
	tree_reco_electron_leadingeta2.push_back(tree_electron_eta[imu2]);
	tree_reco_muon_leadingphi.push_back(tree_muon_phi[imu1] );
	tree_reco_electron_leadingphi2.push_back(tree_electron_phi[imu2] );
	if(tree_only_tigger_filter){
        tree_trig_muon_leadingpt.push_back(tree_muon_pt[imu1] );
        tree_trig_electron_leadingpt2.push_back(tree_electron_pt[imu2]);
        tree_trig_muon_leadingeta.push_back(tree_muon_eta[imu1]);
        tree_trig_electron_leadingeta2.push_back(tree_electron_eta[imu2]);
        tree_trig_muon_leadingphi.push_back(tree_muon_phi[imu1]);
        tree_trig_electron_leadingphi2.push_back(tree_electron_phi[imu2]);
      }//end of trigger loop                                                                                                                                              
  }//end of Good PV   

    lep1_pt  = tree_muon_pt[imu1];
    lep2_pt  = tree_electron_pt[imu2];
    lep1_eta = tree_muon_eta[imu1];
    lep2_eta = tree_electron_eta[imu2];
    lep1_phi = tree_muon_phi[imu1];
    lep2_phi = tree_electron_phi[imu2];
    lep1_Q   = tree_muon_charge[imu1];
    lep2_Q   = tree_electron_charge[imu2];
    lep1_mass = mu_mass;
    lep2_mass = el_mass;
    if ( tree_electron_pt[imu2] > tree_muon_pt[imu1] ) LeadingMuon = false;

  Evts_muon1_pt    = lep1_pt;
  Evts_muon2_pt    = lep2_pt;
  Evts_muon1_eta   = lep1_eta;
  Evts_muon2_eta   = lep2_eta;
  Evts_muon1_phi   = lep1_phi;
  Evts_muon2_phi   = lep2_phi;
  if ( !LeadingMuon ) {
    Evts_muon1_pt  = lep2_pt;
    Evts_muon2_pt  = lep1_pt;
  }
      if(tree_Good_PV){
	tree_reco_lepton_leadingpt.push_back(Evts_muon1_pt );
	tree_reco_lepton_leadingpt2.push_back(  Evts_muon2_pt );
	tree_reco_lepton_leadingeta.push_back(Evts_muon1_eta );
	tree_reco_lepton_leadingeta2.push_back(Evts_muon2_eta );
	tree_reco_lepton_leadingphi.push_back(Evts_muon1_phi );
	tree_reco_lepton_leadingphi2.push_back(Evts_muon2_phi );
	if(tree_only_tigger_filter){
	  tree_trig_lepton_leadingpt.push_back(Evts_muon1_pt );
	  tree_trig_lepton_leadingpt2.push_back(Evts_muon2_pt);
	  tree_trig_lepton_leadingeta.push_back(Evts_muon1_eta );
	  tree_trig_lepton_leadingeta2.push_back(Evts_muon2_eta);
	  tree_trig_lepton_leadingphi.push_back(Evts_muon1_phi);
	  tree_trig_lepton_leadingphi2.push_back(Evts_muon2_phi);
	} // end trigger check                                                                                                                                              
      }// Good PV 
  }

  Vlep1.SetPtEtaPhiM(lep1_pt,lep1_eta,lep1_phi,lep1_mass);
  Vlep2.SetPtEtaPhiM(lep2_pt,lep2_eta,lep2_phi,lep2_mass);
  vll = Vlep1 + Vlep2;
    if(tree_Good_PV && tree_Mmumu > 10){
      tree_ll_pt.push_back(vll.Pt());
      tree_ll_eta.push_back(vll.Eta());
      tree_ll_phi.push_back(vll.Phi());
      tree_ll_px.push_back(vll.Px());
      tree_ll_py.push_back(vll.Py());
      tree_ll_pz.push_back(vll.Pz());
      tree_ll_energy.push_back(vll.Energy());
      tree_ll_mass.push_back(tree_Mmumu);
    }// Good PV and mass cut 
  }// imu1> =0 && imu2>=0       


  //////////////////////////////////
  //////////////////////////////////
  /////////    FILTER    ///////////
  //////////////////////////////////
  //////////////////////////////////

//$$
 if ( tree_Filter || tree_FilterSameSign ) {
//$$


  ///////////////////////////////////////////////
  /////////    lepton informations    ///////////
  ///////////////////////////////////////////////

  if ( !tree_Filter && tree_FilterSameSign ) {
    imu1 = imu1_SS;
    imu2 = imu2_SS;
  }
  mva_Evts_Mmumu = tree_Mmumu;

  edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder); 
  vector<reco::TransientTrack> BestTracks;
  vector<reco::TransientTrack> posBestTracks;
  vector<reco::TransientTrack> negBestTracks;
  int count =0;
  std::map<size_t , int > trackToAK4SlimmedJetMap;


  //////////////////////////////////
  //////////////////////////////////
  //////////    LHE    /////////////
  //////////////////////////////////
  //////////////////////////////////

//   if ( isMC_ )
//   {
//     std::vector<double> evtWeights = genEventInfo->weights();
//     tree_MCEvt_weight = genEventInfo->weight();

// //     const gen::PdfInfo *PDF = genEventInfo->pdf();
// //     // std::cout<<"scalePDF : "<<PDF->scalePDF<<std::endl;
// //     int id1 = PDF->id.first ;// [-4,-3,-2,-1,1,2,3,4]
// //     int id2 = PDF->id.second;
// //     // std::cout<<"id1 and id2 : "<<id1<<"//"<<id2<<std::endl;
// // 
// //     double x1 = PDF->x.first;
// //     double x2 = PDF->x.second;
// //     //  std::cout<<"x1 and x2 : "<<x1<<"//"<<x2<<std::endl;
// // 
// //     double xPDF1 = PDF->xPDF.first;//==0
// //     double xPDF2 = PDF->xPDF.second;//==0
// //     //  std::cout<<"xPDF1 and xPDF2 : "<<xPDF1<<"//"<<xPDF2<<std::endl;
// // 
// //     unsigned int ProcID = genEventInfo->signalProcessID();//9999
// // 	  // double qscale = genEventInfo->qScale();//sameasPDFscale
// //     double alphaqcd = genEventInfo->alphaQCD();
// //     // std::cout<<"ProcID and qscale and alphaqcd : "<<ProcID<<"//"<<" alphaqcd : "<<alphaqcd<<std::endl;

//     for (unsigned int k = 0 ; k<evtWeights.size() ; k++)
//       tree_LHE_Weights.push_back(evtWeights[k]);
//   }


  //////////////////////////////////
  //////////////////////////////////
  /////////     Pileup     /////////
  //////////////////////////////////
  //////////////////////////////////

   // const int nPU = PU->getPU_NumInteractions();
  // std::vector<float> PUzpos;
  // if (nPU)
  //   {
  //     PUzpos = PU->getPU_zpositions();
  //     for ( int i = 0 ; i < nPU ; i++ )
  //       {dRneuneu
  //         std::cout<<"zpositions of PU : "<<PUzpos[i]<<std::endl;
  //       }
  //   }


  //////////////////////////////////
  //////////////////////////////////
  /////////   Simulation   /////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_GenPVx = -1.;
  tree_GenPVy = -1.;
  tree_GenPVz = -20.;

  TLorentzVector vgen;
  float smumass = 0., neumass = 0.;

  int nLLP = 0;
  int nllp = 0;
  // nBC = 0; !
  tree_nFromC = 0; 
  tree_nFromB = 0;
  
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

  int nGenjet = 0, nGenjet1 = 0, nGenjet2 = 0;
  bool isGenjet[99], isGenjet1[99], isGenjet2[99];//the "real" size is given by the final value of jetidx since non-valid jets are replaced
  TLorentzVector Genvaxis1, Genvaxis2, Genvjet[99];
  TLorentzVector GenVaxis1, GenVaxis2, GenVjet[99];
  TLorentzVector GenV1, GenV2, GenV;
  TLorentzVector Genv1, Genv2, Genv;
  float GenEtaMax = 10; // no cut on eta is optimum
  int Genjetidx = 0; // : May be in the loop/ not sure it changes anything
  float GenPtMin = 20;   // (GeV) minimum jet pt is optimum
  float TopWeight = 1.0;
  int ntop = 0;
  if ( isMC_ ) {
  
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
        // if ( i > 3 ) {
        //   cout << " !!! incoming parton to be checked !!! " << endl;
        // 	cout << i << " status " << genIt.status() << " pt eta phi id " 
        // 	  << Gen_pt << " " << Gen_eta << " " << Gen_phi << " " << genIt.pdgId() 
        // 	     << " x y z " << genIt.vx() << " " << genIt.vy() << " " << genIt.vz() << " " << endl; 
        // }
      }

      if (ID == 6)
        {
          if (genIt.isLastCopy())
            {
              ntop++;
              // tree_gen_top_pt.push_back(Gen_pt);
              if (ntop == 1)
                {
                  TopWeight = 0.103*exp(-0.0118*Gen_pt)-0.000134*Gen_pt+0.973;
                  if (showlog)std::cout<<"TopWeight : "<<TopWeight<<" fro the : "<<ntop<<" top quark"<<std::endl;
                }
              if (ntop == 2)
                {
                  TopWeight = sqrt(TopWeight*(0.103*exp(-0.0118*Gen_pt)-0.000134*Gen_pt+0.973));
                  if (showlog) std::cout<<"TopWeight : "<<TopWeight<<" fro the : "<<ntop<<" top quark"<<std::endl;
                }
              tree_genTop_Weight = TopWeight;
            }
        }
      // neutralino from smuon
      if ( ID == 1000023 && abs(mom->pdgId()) == 1000013 ) {
        nLLP++;
              // cout << " neutralino" << nLLP << " pt eta phi " << Gen_pt << " " << Gen_eta << " " << Gen_phi << endl;
        if ( nLLP == 1 ) {
          LLP1_pt  = Gen_pt;
          LLP1_eta = Gen_eta;
          LLP1_phi = Gen_phi;
          LLP1_mother = mom->pdgId();
          smumass = mom->mass();
          neumass = Gen_m;
        }
        if ( nLLP == 2 ) {
          LLP2_pt  = Gen_pt;
          LLP2_eta = Gen_eta;
          LLP2_phi = Gen_phi;
          LLP2_mother = mom->pdgId();
          smumass += mom->mass();
          neumass += Gen_m;
          tree_smu_mass = (smumass/2. + 0.5);
          tree_neu_mass = (neumass/2. + 0.5);
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
        float dEta = Gen_neu2_eta-Gen_neu1_eta;
        float dPhi = Gen_neu2_phi-Gen_neu1_phi;
        tree_genAxis_dPhineuneu.push_back(dPhi);
        tree_genAxis_dEtaneuneu.push_back(dEta);
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
          if ( dV1 > 0.0001 && dV2 > 0.0001 ) nllp++; // should be == 2, so just to check : dV2 is always equal to 0 here
        }
        if ( nllp == 1 ) {
          float dV = (genIt.vx() - LLP1_x)*(genIt.vx() - LLP1_x)
                   + (genIt.vy() - LLP1_y)*(genIt.vy() - LLP1_y)
                   + (genIt.vz() - LLP1_z)*(genIt.vz() - LLP1_z);
          if ( dV > 0.0001 ) {
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
              tree_genFromB_pt.push_back(	     (*packed)[j].pt());
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
                tree_genFromB_dr.push_back((motherInPrunedCollection->vx()-gen2->vx())*(motherInPrunedCollection->vx()-gen2->vx()) 
                                            +(motherInPrunedCollection->vy()-gen2->vy())*(motherInPrunedCollection->vy()-gen2->vy()) );
                tree_genFromB_dz.push_back((motherInPrunedCollection->vz()-gen2->vz())*(motherInPrunedCollection->vz()-gen2->vz()) );
                tree_genFromB_dd.push_back(sqrt((motherInPrunedCollection->vx()-gen2->vx())*(motherInPrunedCollection->vx()-gen2->vx()) 
                                                    +(motherInPrunedCollection->vy()-gen2->vy())*(motherInPrunedCollection->vy()-gen2->vy()) 
                                                    +(motherInPrunedCollection->vz()-gen2->vz())*(motherInPrunedCollection->vz()-gen2->vz()) 
                                              ));
                // std::cout<<"b xyz : "<<(*packed)[j].vx()<<" - "<<(*packed)[j].vy()<<" - "<<(*packed)[j].vz()<<std::endl;
                // std::cout<<"Fromb xyz : "<<gen2->vx()<<" - "<<gen2->vy()<<" - "<<gen2->vz()<<std::endl;
        }
        else { // never happens a priori
                tree_genFromB_x.push_back(-10);
                tree_genFromB_y.push_back(-10);
                tree_genFromB_z.push_back(-10);
        }
            }
      } // final b hadron

      //------------END of Decay length of  b hadrons --------------//
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

      if ( !(abs(genIt.pdgId()) == 6 && abs(mom->pdgId()) == 1000023) 
         && (genIt.pt() < 0.9 || fabs(genIt.eta()) > 4.0) ) continue;

      // std::cout << mom->pdgId() << " mom pdg id "<<std::endl;
      // <=> ((If a particle is a top and comes from a neutralino) OR (has a valid pt and eta )) we look for the information of the gen particle
      tree_genParticle_pt.push_back(        genIt.pt());
      tree_genParticle_eta.push_back(       genIt.eta());
      tree_genParticle_phi.push_back(       genIt.phi());
      tree_genParticle_charge.push_back(    genIt.charge());
      tree_genParticle_pdgId.push_back(     genIt.pdgId());
      tree_genParticle_mass.push_back(      genIt.mass());
      tree_genParticle_x.push_back(	    genIt.vx());
      tree_genParticle_y.push_back(	    genIt.vy());
      tree_genParticle_z.push_back(	    genIt.vz());
      tree_genParticle_px.push_back(	    genIt.px());
      tree_genParticle_py.push_back(	    genIt.py());
      tree_genParticle_pz.push_back(	    genIt.pz());
      tree_genParticle_energy.push_back(    genIt.energy());
      tree_genParticle_isPromptFinalState.push_back( genIt.isPromptFinalState());
      tree_genParticle_statusCode.push_back(genIt.status());
      tree_genParticle_mother_pdgId.push_back( mom ? mom->pdgId() :  -10 );
      tree_genParticle_LLP.push_back(fromLLP);

      // get generated lifetime of neutralino (stored for each of its top)
      float ct  = -1.;
      float ct0 = -1.;
      if ( abs(genIt.pdgId()) == 6 && abs(mom->pdgId()) == 1000023 ) {
        float dx = genIt.vx()-tree_GenPVx;
        float dy = genIt.vy()-tree_GenPVy;
        float dz = genIt.vz()-tree_GenPVz;
        ct = sqrt( dx*dx + dy*dy + dz*dz ); // cm
        vgen.SetPtEtaPhiM(mom->pt(),mom->eta(),mom->phi(),mom->mass());
        float bg = vgen.P() / mom->mass();
        ct0 = ct / bg;
	      tree_neu_ctau = ct0;
      }
      tree_genParticle_ct.push_back(ct);
      tree_genParticle_ct0.push_back(ct0);
    } // end loop on pruned genparticles

    tree_nLLP = nllp;

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
    } // end of loop on packed particles

    //////////////////////
    // gen jets
    //////////////////////

    // genMuPair
    int nGenLepton = 0;
    for (auto const & genJet : *genJets)
    {
      if ( genJet.pt() < GenPtMin ) continue;
      tree_genJet_pt.push_back(genJet.pt());
      tree_genJet_eta.push_back(genJet.eta());
      tree_genJet_phi.push_back(genJet.phi());
      tree_genJet_mass.push_back(genJet.mass());
      tree_genJet_energy.push_back(genJet.energy());

      float Genjet_pt  = genJet.pt();
      float Genjet_eta = genJet.eta();
      float Genjet_phi = genJet.phi();
      float Genjet_px = genJet.px();
      float Genjet_py = genJet.py();
      float Genjet_pz = genJet.pz();
      float Genjet_e = genJet.energy();

      isGenjet[Genjetidx]  = false;
      isGenjet1[Genjetidx] = false; // first neutralino jets
      isGenjet2[Genjetidx] = false; // second neutralino jets
      Genv.SetPtEtaPhiM( Genjet_pt, Genjet_eta, Genjet_phi, 0. ); //set the axis
      GenV.SetPxPyPzE(Genjet_px,Genjet_py,Genjet_pz,Genjet_e);
      if ( Genjet_pt < GenPtMin ) continue;
      if ( abs(Genjet_eta) > GenEtaMax ) continue;

     //----------------------------------------//
     //-----------GenEvent Axes----------------//
     //----------------------------------------//

      nGenjet++;
      isGenjet[Genjetidx] = true;
      Genvjet[Genjetidx] = Genv; // Only jet data (with  possible muons being removed)
      GenVjet[Genjetidx] = GenV;
      if ( nGenjet1 == 0 && Genjet_pt > GenPtMin && abs(Genjet_eta) < GenEtaMax && Genjetidx==0 )
      {
        nGenjet1 = 1;
        isGenjet1[Genjetidx] = true;
        Genvaxis1 = Genv;
        GenVaxis1 = GenV;
      }
      
      float GendR1 = 10., GendR2 = 10.;
      float GendRcut_hemis  = 1.5; // subjective choice default is 1.5

      if ( !isGenjet[Genjetidx] ) continue;
      float jet_eta = Genvjet[Genjetidx].Eta();
      float jet_phi = Genvjet[Genjetidx].Phi();
      if ( nGenjet1 > 0 ) GendR1 = Deltar( jet_eta, jet_phi, Genvaxis1.Eta(), Genvaxis1.Phi() );
      if ( nGenjet2 > 0 ) GendR2 = Deltar( jet_eta, jet_phi, Genvaxis2.Eta(), Genvaxis2.Phi() );

            // axis 1  <---------------------
      if ( nGenjet1 > 0 && !isGenjet2[Genjetidx]  && GendR1 < GendRcut_hemis) {
        nGenjet1++;
        Genvaxis1 += Genvjet[Genjetidx];
        GenVaxis1 += GenVjet[Genjetidx];
        isGenjet1[Genjetidx] = true;
      }
            // axis 2  ---------->
      if ( nGenjet2 == 0 && !isGenjet1[Genjetidx] ) {
        nGenjet2 = 1;
        Genvaxis2 = Genvjet[Genjetidx];
        GenVaxis2 += GenVjet[Genjetidx];
      }
            // axis 2  --------------------->
      else if ( nGenjet2 > 0 && !isGenjet1[Genjetidx] && !isGenjet2[Genjetidx] && GendR2 < GendRcut_hemis) {//
        nGenjet2++;
        Genvaxis2 += Genvjet[Genjetidx];
        GenVaxis2 += GenVjet[Genjetidx];
        isGenjet2[Genjetidx] = true;
      }
      Genjetidx++;
      // loop over genmuons and remove the contribution of the prompt muosn to the axes building procedure 

      float deltaGenR1 = 1000.;
      float deltaGenR2 = 1000.;

      // int nGenMuon = 0;
      // int nGeneleMuon = 0;//changes meena  
      nGenLepton = 0;
      for (unsigned int k = 0 ; k < tree_genParticle_pdgId.size() ; k++) //loop over genMuons
      {
        if (abs(tree_genParticle_pdgId[k])!=13 && MuonChannel) continue;//pruned collection may be should check also with the packed colelction
        if ( ((abs(tree_genParticle_pdgId[k])!=11) || (abs(tree_genParticle_pdgId[k])!=13)) && (EMuChannel) ) continue;// changes meena // check 
        if ((abs(tree_genParticle_pdgId[k])!=11) && ElChannel) continue;
        //feels like there are sometimes two gen muons that are ony one?? close to having the same pt eta and phi (deltaQuantity  ~ 0.001) => continue => splitting oftracks of muons
        //adressed in reco
        if ( !tree_genParticle_isPromptFinalState[k] ) continue; //is Prompt and Final state (not from hadron, muon, or tau decay))
        if (abs(tree_genParticle_mother_pdgId[k]) != 1000013) continue;
        if ( tree_genParticle_pt[k] < 25) continue;
        // nGenMuon++;
        // nGeneleMuon++;
        nGenLepton++;
        deltaGenR1 = Deltar( Genvaxis1.Eta(), Genvaxis1.Phi(), tree_genParticle_eta[k], tree_genParticle_phi[k] );
        deltaGenR2 = Deltar( Genvaxis2.Eta(), Genvaxis2.Phi(), tree_genParticle_eta[k], tree_genParticle_phi[k] );
        if ( deltaGenR1 < 0.4)
          {
            Genv1.SetPtEtaPhiM( tree_genParticle_pt[k],tree_genParticle_eta[k], tree_genParticle_phi[k], 0 ); 
            Genvaxis1 -= Genv1; // v TLorentzFactor being just above, defined by jet data
            GenV1.SetPxPyPzE(tree_genParticle_px[k],tree_genParticle_py[k],tree_genParticle_pz[k],tree_genParticle_energy[k]);
            GenVaxis1 -= GenV1;
          }
        if ( deltaGenR2 < 0.4 )
          {
            Genv2.SetPtEtaPhiM( tree_genParticle_pt[k],tree_genParticle_eta[k], tree_genParticle_phi[k], 0 );
            GenV2.SetPxPyPzE(tree_genParticle_px[k],tree_genParticle_py[k],tree_genParticle_pz[k],tree_genParticle_energy[k]);
            Genvaxis2 -= Genv2;
            GenVaxis2 -= GenV2;
          }
      } // end loop over gen muons

// std::cout<<"ngen prompt muons = "<<nGenLepton<<std::endl;
    } // end loop over gen jets

    // ///////////////////////////////
    // // Invariant Mass of GenAxes 
    // ///////////////////////////////
    
    float temp_px1 = GenVaxis1.Px();
    float temp_py1 = GenVaxis1.Py();
    float temp_pz1 = GenVaxis1.Pz();
    float temp_e1  = GenVaxis1.E();
    float temp_px2 = GenVaxis2.Px();
    float temp_py2 = GenVaxis2.Py();
    float temp_pz2 = GenVaxis2.Pz();
    float temp_e2  = GenVaxis2.E();

    //... then, add the lepton contribution
    for (unsigned int i = 0 ; i<tree_genParticle_pdgId.size() ; i++)
    {
      if (abs(tree_genParticle_pdgId[i])!=13 && MuonChannel) continue;//pruned collection may be should check also with the packed colelction
      if ( ((abs(tree_genParticle_pdgId[i])!=11) || (abs(tree_genParticle_pdgId[i])!=13)) && (EMuChannel) ) continue;// changes meena // check 
      if ((abs(tree_genParticle_pdgId[i])!=11) && ElChannel) continue;
      if ( tree_genParticle_isPromptFinalState[i] ) continue;//||  !tree_muon_isTight[mu]
      
      // Avoid prompt leptons
      float GendRAxis1muon = Deltar( Genvaxis1.Eta(), Genvaxis1.Phi(), tree_genParticle_eta[i], tree_genParticle_phi[i] );
      float GendRAxis2muon = Deltar( Genvaxis2.Eta(), Genvaxis2.Phi(), tree_genParticle_eta[i], tree_genParticle_phi[i] );
      if ( GendRAxis1muon < 1.5 )
    	{
    	  temp_px1 += tree_genParticle_px[i];
    	  temp_py1 += tree_genParticle_py[i];
    	  temp_pz1 += tree_genParticle_pz[i];
    	  temp_e1  += tree_genParticle_energy[i];
    	}
      if ( GendRAxis2muon < 1.5 )
    	{
    	  temp_px2 += tree_genParticle_px[i];
    	  temp_py2 += tree_genParticle_py[i];
    	  temp_pz2 += tree_genParticle_pz[i];
    	  temp_e2  += tree_genParticle_energy[i];
    	}
    }

    TLorentzVector TLorentzGenAxis1(temp_px1,temp_py1,temp_pz1,temp_e1);
    TLorentzVector TLorentzGenAxis2(temp_px2,temp_py2,temp_pz2,temp_e2);
    tree_GenAxes_Mass.push_back(TLorentzGenAxis1.Mag());
    tree_GenAxes_Mass.push_back(TLorentzGenAxis2.Mag());

    int count = 0 ;
    int idx1 = -1;
    int idx2 = -1;
    for (unsigned int i=0 ; i<tree_genParticle_pdgId.size() ; i++)
    {
      if (abs(tree_genParticle_pdgId[i])!=13 && MuonChannel) continue;//pruned collection may be should check also with the packed colelction
      if ( !tree_genParticle_isPromptFinalState[i] ) continue;//||  !tree_muon_isTight[mu]
      if (abs(tree_genParticle_mother_pdgId[i]) != 1000013) continue;
      if (nGenLepton>=1 && count==0 && tree_genParticle_pt[i] > 25)
        {
          count++;
          idx1 = i;
        }
      else if (nGenLepton>1 && count>0 && tree_genParticle_pt[i] > 10)
        {
          count++;
          idx2=i;
        }
    }
    // std::cout<<"count idx1 and idx2: "<<count<<" // "<<idx1<<" // "<<idx2<<std::endl;
    //add prompt gen muon selection
    float dR1 = Deltar(tree_genParticle_eta[0],tree_genParticle_phi[0],TLorentzGenAxis1.Eta(),TLorentzGenAxis1.Phi());
    float dR2 = Deltar(tree_genParticle_eta[0],tree_genParticle_phi[0],TLorentzGenAxis2.Eta(),TLorentzGenAxis2.Phi());
    // std::cout<<"dR1 and dR2 : "<<dR1<<" : "<<dR2<<std::endl;
    if ( count > 0 ) {
      if ( dR1 < dR2 ) {
        temp_px1+= tree_genParticle_px[idx1];
        temp_py1+= tree_genParticle_py[idx1];
        temp_pz1+= tree_genParticle_pz[idx1];
        temp_e1+= tree_genParticle_energy[idx1];
        if ( count >= 2 ) {
          temp_px2+= tree_genParticle_px[idx2];
          temp_py2+= tree_genParticle_py[idx2];
          temp_pz2+= tree_genParticle_pz[idx2];
          temp_e2+= tree_genParticle_energy[idx2];
        }
      }
      else {
        temp_px1+= tree_genParticle_px[idx2];
        temp_py1+= tree_genParticle_py[idx2];
        temp_pz1+= tree_genParticle_pz[idx2];
        temp_e1+= tree_genParticle_energy[idx2];
        if ( count >= 2 ) {
          temp_px2+= tree_genParticle_px[idx1];
          temp_py2+= tree_genParticle_py[idx1];
          temp_pz2+= tree_genParticle_pz[idx1];
          temp_e2+= tree_genParticle_energy[idx1];
        }
      }
    }
  
    TLorentzVector TLorentzCombinedAxis1(temp_px1,temp_py1,temp_pz1,temp_e1);
    TLorentzVector TLorentzCombinedAxis2(temp_px2,temp_py2,temp_pz2,temp_e2);

    float CombinedMass1 = TLorentzCombinedAxis1.Mag();
    float CombinedMass2 = TLorentzCombinedAxis2.Mag();
    if (isnan(CombinedMass1) || isinf(CombinedMass1)) CombinedMass1 = 0;
    if (isnan(CombinedMass2) || isinf(CombinedMass2)) CombinedMass2 = 0;
    tree_GenAxes_CombinedHemiLeptonMass.push_back(CombinedMass1);
    tree_GenAxes_CombinedHemiLeptonMass.push_back(CombinedMass2);

  } // endif MC


  //////////////////////////////////
  //////////////////////////////////
  ///////////   MET   //////////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_PFMet_et  = -10.;
  tree_PFMet_phi = -10.;
  tree_PFMet_sig = -10.;
  tree_PFMet_pt = -10.;
  if ( PFMETs->size() > 0 ) {
    const pat::MET &themet = PFMETs->front();
    tree_PFMet_et  = themet.et();
    tree_PFMet_phi = themet.phi();
    tree_PFMet_sig = themet.significance();
    float met_x   =  themet.px();
    float met_y   =  themet.py();
    tree_PFMet_pt  = sqrt(met_x*met_x + met_y*met_y);
  }


  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  //////    V0 Candidates  from CMSSW Collection   //////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////

  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideV0Producer
  // We actually do not use these ones => See V0 Producer part
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
      tree_K0_r.push_back(	TMath::Sqrt(K0x*K0x+K0y*K0y));
      tree_K0_NChi2.push_back(  (*KshortVertex)[j].vertexNormalizedChi2());
      tree_K0_ndf.push_back(	(*KshortVertex)[j].vertexNdof());
      tree_K0_mass.push_back(	(*KshortVertex)[j].mass());
      tree_K0_pt.push_back(	(*KshortVertex)[j].pt());
      tree_K0_eta.push_back(	(*KshortVertex)[j].eta());
      tree_K0_phi.push_back(	(*KshortVertex)[j].phi());
      tree_K0_nDaughters.push_back((*KshortVertex)[j].numberOfDaughters());
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
      tree_L0_r.push_back(         TMath::Sqrt(L0x*L0x+L0y*L0y));
      tree_L0_NChi2.push_back(     (*LambdaVertex)[j].vertexNormalizedChi2());
      tree_L0_ndf.push_back(       (*LambdaVertex)[j].vertexNdof());
      tree_L0_mass.push_back(      (*LambdaVertex)[j].mass());
      tree_L0_pt.push_back(        (*LambdaVertex)[j].pt());
      tree_L0_eta.push_back(       (*LambdaVertex)[j].eta());
      tree_L0_phi.push_back(       (*LambdaVertex)[j].phi());
      tree_L0_nDaughters.push_back((*LambdaVertex)[j].numberOfDaughters());
    }
  }


  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  ///      Photon Conversion from CMSSW Collection   ////
  ///////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////
  // https://github.com/cms-sw/cmssw/tree/CMSSW_10_6_X/DataFormats/EgammaCandidates/interface/Conversion.h

  tree_nYConv = PhotonConversion->size();
  tree_Yc_ntracks = 0;
  int iYc = -1;
  if ( !PhotonConversion->empty() ) {
    for (int j = 0; j < tree_nYConv ; j++)
    {
      if ( !((*PhotonConversion)[j].isConverted()) ) continue; // tracksize>0 // Number of tracks= 0,1,2 unsigned int nTracks() const {return  tracks().size(); }
      iYc++;
      float Yx = (*PhotonConversion)[j].conversionVertex().position().x();
      float Yy = (*PhotonConversion)[j].conversionVertex().position().y();
      float Yz = (*PhotonConversion)[j].conversionVertex().position().z();
      float Yr = TMath::Sqrt(Yx*Yx + Yy*Yy);
      tree_Yc_x.push_back(	   Yx); 
      tree_Yc_y.push_back(	   Yy);
      tree_Yc_z.push_back(	   Yz);
      tree_Yc_r.push_back(	   Yr);
      
      int VtxLayerNI = -10;
      if ( DetailedMap ) VtxLayerNI = NI->VertexBelongsToTracker(Yr, Yz);
      else { // deprecated
          VtxLayerNI = NI->VertexBelongsToBarrelLayer(Yr, Yz);
          if ( VtxLayerNI == 0 ) VtxLayerNI = NI->VertexBelongsToDiskLayer(Yr, Yz);
        }
      tree_Yc_layer.push_back(     VtxLayerNI);

      float Ynchi2 = (*PhotonConversion)[j].conversionVertex().normalizedChi2();
      tree_Yc_NChi2.push_back(     Ynchi2);
      tree_Yc_ndf.push_back(	   (*PhotonConversion)[j].conversionVertex().ndof());	   
      tree_Yc_nDaughters.push_back((*PhotonConversion)[j].nTracks());

      std::vector<math::XYZPointF>  inPos = (*PhotonConversion)[j].tracksInnerPosition();
      std::vector<math::XYZVectorF> inP   = (*PhotonConversion)[j].tracksPin();
      std::vector<reco::Track> YcVertexTracks = (*PhotonConversion)[j].conversionVertex().refittedTracks();

      TLorentzVector vY;
      vY.SetXYZM( YcVertexTracks[0].px()+YcVertexTracks[1].px(), 
        	  YcVertexTracks[0].py()+YcVertexTracks[1].py(), 
        	  YcVertexTracks[0].pz()+YcVertexTracks[1].pz(), 0. );
      tree_Yc_pt.push_back(	   vY.Pt()); 
      tree_Yc_eta.push_back(	   vY.Eta()); 
      tree_Yc_phi.push_back(	   vY.Phi()); 
      float Ymass  = (*PhotonConversion)[j].pairInvariantMass();
      tree_Yc_mass.push_back(	   Ymass);	 

      // check if tracks are starting downstream from the conversion vertex
      float Ydr0 = TMath::Sqrt((Yx-inPos[0].x())*(Yx-inPos[0].x()) + (Yy-inPos[0].y())*(Yy-inPos[0].y()));
      float Ydr1 = TMath::Sqrt((Yx-inPos[1].x())*(Yx-inPos[1].x()) + (Yy-inPos[1].y())*(Yy-inPos[1].y()));
      if ( TMath::Sqrt(inPos[0].x()*inPos[0].x() + inPos[0].y()*inPos[0].y()) < Yr ) Ydr0 = -Ydr0;
      if ( TMath::Sqrt(inPos[1].x()*inPos[1].x() + inPos[1].y()*inPos[1].y()) < Yr ) Ydr1 = -Ydr1;
      float Ydz0 = abs( Yz - inPos[0].z() );
      float Ydz1 = abs( Yz - inPos[1].z() );
      if ( abs(inPos[0].z()) < abs(Yz) ) Ydz0 = -Ydz0;
      if ( abs(inPos[1].z()) < abs(Yz) ) Ydz1 = -Ydz1;
      tree_Yc_dr0.push_back(	   Ydr0); 
      tree_Yc_dr1.push_back(	   Ydr1); 
      tree_Yc_dz0.push_back(	   Ydz0); 
      tree_Yc_dz1.push_back(	   Ydz1);
      // angle between the photon direction and its momentum vector
      float costheta = vY.Px()*(Yx-tree_PV_x) + vY.Py()*(Yy-tree_PV_y) + vY.Pz()*(Yz-tree_PV_z);
      costheta = costheta / vY.P() / TMath::Sqrt( (Yx-tree_PV_x)*(Yx-tree_PV_x) + (Yy-tree_PV_y)*(Yy-tree_PV_y) + (Yz-tree_PV_z)*(Yz-tree_PV_z) );
      tree_Yc_costheta.push_back(  costheta);

    if ( !(VtxLayerNI != 0 && Ymass > 0. && Ymass < 1. && Ynchi2 < 10.) ) continue;
 
      for (unsigned int k = 0; k < (*PhotonConversion)[j].nTracks(); k++)
      {
        tree_Yc_ntracks++;
        TLorentzVector vYtk;
        vYtk.SetXYZM( inP[k].x(), inP[k].y(), inP[k].z(), 0. );
        float Ytk_pt  = vYtk.Pt();
        float Ytk_eta = vYtk.Eta();
        float Ytk_phi = vYtk.Phi();
        float Ytk_q   = YcVertexTracks[k].charge();
        float qR = Ytk_q * Ytk_pt * 100 / 0.3 / 3.8;
        float sin0 = qR * sin( Ytk_phi ) + (inPos[k].x() - tree_PV_x);
        float cos0 = qR * cos( Ytk_phi ) - (inPos[k].y() - tree_PV_y);
        float Ytk_phi0 = TMath::ATan2( sin0, cos0 ); // but note that it can be wrong by +_pi ! 
        tree_Yc_tracks_index.push_back(  iYc);
        tree_Yc_tracks_charge.push_back( Ytk_q);
        tree_Yc_tracks_pt.push_back(	 Ytk_pt);
        tree_Yc_tracks_eta.push_back(	 Ytk_eta);
        tree_Yc_tracks_phi.push_back(	 Ytk_phi);
        tree_Yc_tracks_phi0.push_back(   Ytk_phi0);
      }
    }
  }


  //////////////////////////////////
  //////////////////////////////////
  ///////////   Leptons  ///////////
  //////////////////////////////////
  //////////////////////////////////
 
    mva_Evts_muon1_pt    = lep1_pt;
    mva_Evts_muon2_pt    = lep2_pt;
    mva_Evts_muon12_dR   = Deltar(lep1_eta, lep1_phi, lep2_eta, lep2_phi);
    mva_Evts_muon12_dPhi = abs( Deltaphi(lep1_phi, lep2_phi) );
    mva_Evts_muon12_dEta = abs( lep1_eta - lep2_eta );

  tree_lepton_leadingpt.push_back(   mva_Evts_muon1_pt );  
  tree_lepton_leadingpt2.push_back(  mva_Evts_muon2_pt );
  tree_lepton_leadingeta.push_back( lep1_eta );
  tree_lepton_leadingeta2.push_back( lep2_eta );
  tree_lepton_leadingphi.push_back( lep1_phi );
  tree_lepton_leadingphi2.push_back( lep2_phi );


  tree_lepton_lepton_dR.push_back(   mva_Evts_muon12_dR );
  tree_lepton_lepton_dPhi.push_back( mva_Evts_muon12_dPhi );
  tree_lepton_lepton_dEta.push_back( mva_Evts_muon12_dEta );


  //////////////////////////////////
  //////////////////////////////////
  ///////////	Jets   /////////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_njet = 0;      // all selected jets
  tree_njetNOmu = 0;  // selected jets without a prompt muon track inside
  float HT_val = 0;
  float jet_pt_min = 20.;
  int indjet = -1;

  float jet1_pt = 0;
  float jet1_eta = 0;
  float jet1_phi = -10;

  float jet2_pt = 0;
  float jet2_eta = 0; 
  float jet2_phi = -10;

  for (const pat::Jet &jet : *jets) 
  {
    
  if ( jet.pt() < jet_pt_min ) continue;
  if ( !jet.userInt("tightLepVetoId") ) continue;
    indjet++;
    float NHF                 = jet.neutralHadronEnergyFraction();
    float NEMF                = jet.neutralEmEnergyFraction();
    float CHF                 = jet.chargedHadronEnergyFraction();
    float MUF                 = jet.muonEnergyFraction();
    float CEMF                = jet.chargedEmEnergyFraction();
    float NumConst            = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    float NumNeutralParticles = jet.neutralMultiplicity();
    float CHM                 = jet.chargedMultiplicity(); 
    bool TightJetIDLepVeto = false;
    bool TightJetID        = false;
    if ( abs(jet.eta()) <= 2.6 )
    {
      if ( NHF<0.9 && NEMF<0.9 && NumConst>1 && CHF>0 && CHM>0 )
      {
        TightJetID = true;
        if ( MUF<0.8 && CEMF<0.8 ) TightJetIDLepVeto = true; // lepton veto
      }
    }
    else if ( abs(jet.eta()) <= 2.7 )
    {
      if ( NHF<0.9 && NEMF<0.99  && CHM>0 )
      {
        TightJetID = true;
        if ( MUF<0.8 && CEMF<0.8 ) TightJetIDLepVeto = true; // lepton veto
      }
    }
    else if ( abs(jet.eta()) <= 3.0 )
    {
      if ( NEMF<0.99 && NEMF>0.01 && NumNeutralParticles>1 )
      {
        TightJetID = true;
        TightJetIDLepVeto = true;
      }
    }
    else if ( abs(jet.eta()) <= 5.0 )
    {
      if ( NEMF<0.90 && NHF>0.2 && NumNeutralParticles>10 )
      {
        TightJetID = true;
        TightJetIDLepVeto = true;
      }
    }
    tree_jet_TightJetIDLepVeto.push_back(TightJetIDLepVeto);
    tree_jet_TightJetID.push_back(TightJetID);
    tree_jet_tightid.push_back(jet.userInt("tightId"));
    tree_jet_tightid_LepVeto.push_back(jet.userInt("tightLepVetoId"));

    if ( indjet == 0 ) {jet1_pt = jet.pt(); jet1_eta = jet.eta();jet1_phi=jet.phi();}
    if ( indjet == 1 ) {jet2_pt = jet.pt(); jet2_eta = jet.eta();jet2_phi=jet.phi();}
    if (jet1_pt < jet2_pt){float tempjetpt = jet2_pt; jet2_pt = jet1_pt; jet1_pt = tempjetpt; }
    if (indjet==1 && showlog){std::cout<<"delta R between leading and subleadingjet : "<< Deltar(jet1_eta,jet1_phi,jet2_eta,jet2_phi) <<std::endl;}
    tree_jet_pt.push_back(jet.pt());
    tree_jet_eta.push_back(jet.eta());
    tree_jet_phi.push_back(jet.phi());
    tree_jet_px.push_back(jet.px());
    tree_jet_py.push_back(jet.py());
    tree_jet_pz.push_back(jet.pz());
    tree_jet_E.push_back(jet.energy());
    // std::cout<<" jet.n90() : "<<jet.n90()<<" and  n60: "<<jet.n60()<<std::endl;

    tree_jet_HadronFlavour.push_back(jet.hadronFlavour());
    tree_jet_pileupID.push_back(jet.userFloat("pileupJetId:fullDiscriminant"));

    // btag infos :
    // WorkingPoints : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    float DeepCSVb = jet.bDiscriminator("pfDeepCSVJetTags:probb");
    float DeepCSVbb = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
    float DeepFlavourb = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
    float DeepFlavourbb = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
    float DeepFlavourblep = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
    float DeepCSV = DeepCSVb + DeepCSVbb;
    float DeepJet = -10.;
    if ( DeepFlavourb > -5 ) DeepJet = DeepFlavourb + DeepFlavourbb + DeepFlavourblep;
    tree_jet_btag_DeepCSV.push_back( DeepCSV );
    tree_jet_btag_DeepJet.push_back( DeepJet );

    if ( abs(jet.eta()) < 2.4 ) HT_val += jet.pt();
    if (DeepJet > MediumWP && nmu>=1 )
    {
      tree_jet_leadingMuon_dR.push_back(              Deltar( jet.eta(), jet.phi(), lep1_eta, lep1_phi ));
      if ( nmu>1 ) tree_jet_leadingMuon2_dR.push_back(Deltar( jet.eta(), jet.phi(), lep2_eta, lep2_phi ));
      else         tree_jet_leadingMuon2_dR.push_back(0); 
    }
    else           tree_jet_leadingMuon_dR.push_back(0);
    tree_njet++;

    bool PromptLeptonJetMatching = false ; 
        float dRmujet1 = Deltar( jet.eta(), jet.phi(), lep1_eta, lep1_phi );
    float dRmujet2 = Deltar( jet.eta(), jet.phi(), lep2_eta, lep2_phi );
        if ( dRmujet1 < 0.4 || dRmujet2 < 0.4 ) PromptLeptonJetMatching = true;
        if ( !PromptLeptonJetMatching ) tree_njetNOmu++;

  } // end jet loop

  tree_HT = HT_val;

    /////////////////////////////////////////////////////
  ////////////   ONLY EVENTS WITH JETS !   ////////////
  /////////////////////////////////////////////////////

//   if ( tree_njetNOmu == 0 ) // not tree_njet ?
  //   {
//     tree_jet_jet_dR.push_back(0);
    //     tree_jet_jet_dPhi.push_back(0);
    //     tree_jet_jet_dEta.push_back(0);
    //     tree_jet_leadingpt.push_back(0);
    //     tree_jet_leadingpt2.push_back(0);
    //     mva_Evts_jet12_dR = 0;
    //     mva_Evts_jet12_dPhi = 0;
    //     mva_Evts_jet12_dEta  = 0;
    //     mva_Evts_jet1_pt = 0;
    //     mva_Evts_jet2_pt = 0;
  //   }

//$$
  if ( tree_njetNOmu > 0 ) {
//$$

  
  //////////////////////////////////
  //////////////////////////////////
  //     EVTS Selection BDT       //
  //////////////////////////////////
  //////////////////////////////////

  if ( tree_njetNOmu == 1 )
  {
    tree_jet_jet_dR.push_back(0);
    tree_jet_jet_dPhi.push_back(0);
    tree_jet_jet_dEta.push_back(0);
    tree_jet_leadingpt.push_back(jet1_pt);
    tree_jet_leadingpt2.push_back(0);
    tree_jet_leadingeta.push_back(jet1_eta);
    tree_jet_leadingeta2.push_back(-10);
    mva_Evts_jet12_dR = 0;
    mva_Evts_jet12_dPhi = 0;
    mva_Evts_jet12_dEta  = 0;
    mva_Evts_jet1_pt = jet1_pt;
    mva_Evts_jet2_pt = 0;
    mva_Evts_jet1_eta = jet1_eta;
    mva_Evts_jet2_eta = 0;
  }
  
  if ( tree_njetNOmu >= 2 ) 
  {
    tree_jet_jet_dR.push_back(Deltar(jet1_eta,tree_jet_phi[0],jet2_eta,tree_jet_phi[1]));
    tree_jet_jet_dPhi.push_back(abs(Deltaphi(tree_jet_phi[0],tree_jet_phi[1])));
    tree_jet_jet_dEta.push_back(abs(jet1_eta-jet2_eta));
    tree_jet_leadingpt.push_back(jet1_pt);
    tree_jet_leadingpt2.push_back(jet2_pt);
    tree_jet_leadingeta.push_back(jet1_eta);
    tree_jet_leadingeta2.push_back(jet2_eta);
    mva_Evts_jet12_dR = Deltar(jet1_eta,tree_jet_phi[0],jet2_eta,tree_jet_phi[1]);
    mva_Evts_jet12_dPhi = abs(Deltaphi(tree_jet_phi[0],tree_jet_phi[1]));
    mva_Evts_jet12_dEta  = abs(jet1_eta-jet2_eta);
    mva_Evts_jet1_pt = jet1_pt;
    mva_Evts_jet2_pt = jet2_pt;
    mva_Evts_jet1_eta = jet1_eta;
    mva_Evts_jet2_eta = jet2_eta;   
  }

  mva_Evts_MET_et           = tree_PFMet_et;
  // mva_Evts_nVtx             = ;
  mva_HT                    = HT_val;
  mva_Evts_ST               = tree_LT;
  mva_Evts_njets            = tree_njetNOmu;
  mva_Evts_nmuon            = nmu;
  mva_Evts_all_muon         = allnmu;
  // mva_Evts_MediumAxes       = ;
  // mva_Evts_LooseAxes        = ;
  // mva_Evts_TightAxes        = ;
  //the other variables are in loops so they can not be put there
  // but here is the list:

    // mva_Evts_MET_et
    // mva_Evts_nTrks = theFormat.tree_TRACK_SIZE;
    // mva_Evts_muon1_pt= theFormat.tree_lepton_leadingpt->at(0);
    // mva_Evts_muon2_pt= theFormat.tree_lepton_leadingpt2->at(0);
    // mva_Evts_jet1_pt = theFormat.tree_jet_leadingpt->at(0);
    // mva_Evts_jet2_pt = theFormat.tree_jet_leadingpt2->at(0);
    // mva_Evts_muon12_dR = theFormat.tree_lepton_lepton_dR->at(0);
    // mva_Evts_muon12_dPhi = theFormat.tree_lepton_lepton_dPhi->at(0);
    // mva_Evts_muon12_dEta = theFormat.tree_lepton_lepton_dEta->at(0);
    // mva_Evts_jet12_dR = theFormat.tree_jet_jet_dR->at(0);
    // mva_Evts_jet12_dPhi = theFormat.tree_jet_jet_dPhi->at(0);
    // mva_Evts_jet12_dEta  = theFormat.tree_jet_jet_dEta->at(0);
    // mva_Evts_muon_jet_dRmin = theFormat.tree_muon_jet_dRmin->at(0);
    // mva_Evts_muon_jet_dRmax = theFormat.tree_muon_jet_dRmax->at(0);
    // mva_Evts_nVtx = theFormat.tree_Hemi_Vtx_nVtx->at(0);
    // mva_HT = theFormat.tree_HT;
    // mva_Evts_ST = theFormat.tree_LT->at(0);
    // mva_Evts_njets = theFormat.tree_njetNOmu;
    // mva_Evts_nmuon 

  float dRmuon1_jet_min = 0;
  float dRmuon1_jet_max = 0;
  float dRmuon2_jet_min = 0;
  float dRmuon2_jet_max = 0;
  if ( tree_njetNOmu >= 2 ) {
    float dRmuon1_jet0 = Deltar(lep1_eta,lep1_phi,tree_jet_eta[0],tree_jet_phi[0]);
    float dRmuon1_jet1 = Deltar(lep1_eta,lep1_phi,tree_jet_eta[1],tree_jet_phi[1]);
    if (dRmuon1_jet0 < dRmuon1_jet1) {
      dRmuon1_jet_min = dRmuon1_jet0;
      dRmuon1_jet_max = dRmuon1_jet1;
    }
    else {
      dRmuon1_jet_min = dRmuon1_jet1;
      dRmuon1_jet_max = dRmuon1_jet0;
    }
    float dRmuon2_jet0 = Deltar(lep2_eta,lep2_phi,tree_jet_eta[0],tree_jet_phi[0]);
    float dRmuon2_jet1 = Deltar(lep2_eta,lep2_phi,tree_jet_eta[1],tree_jet_phi[1]);
    if (dRmuon2_jet0 < dRmuon2_jet1) {
      dRmuon2_jet_min = dRmuon2_jet0;
      dRmuon2_jet_max = dRmuon2_jet1;
    }
    else {
      dRmuon2_jet_min = dRmuon2_jet1;
      dRmuon2_jet_max = dRmuon2_jet0;
    }
  }

  if (dRmuon1_jet_min < dRmuon2_jet_min) {
      tree_muon_jet_dRmin.push_back(dRmuon1_jet_min);
      tree_muon_jet_dRmin.push_back(dRmuon2_jet_min);
    }
    else {
      tree_muon_jet_dRmin.push_back(dRmuon2_jet_min);
      tree_muon_jet_dRmin.push_back(dRmuon1_jet_min);
    }

    if (dRmuon1_jet_max > dRmuon2_jet_max) {
      tree_muon_jet_dRmax.push_back(dRmuon1_jet_max);
      tree_muon_jet_dRmax.push_back(dRmuon2_jet_max);
    }
    else {
      tree_muon_jet_dRmax.push_back(dRmuon2_jet_max);
      tree_muon_jet_dRmax.push_back(dRmuon1_jet_max);
    } 
    mva_Evts_muon_jet_dRmin0 = tree_muon_jet_dRmin[0];
    mva_Evts_muon_jet_dRmax0 = tree_muon_jet_dRmax[0];
    mva_Evts_muon_jet_dRmin1 = tree_muon_jet_dRmin[1];
    mva_Evts_muon_jet_dRmax1 = tree_muon_jet_dRmax[1];


  /////////////////////////////////////////////////////////
  //-------------------------------------------------------
  // Jets for event axes				 
  //-------------------------------------------------------
  /////////////////////////////////////////////////////////

  int njet1 = 0, njet2 = 0, njet1_nomu = 0, njet2_nomu = 0;
  bool isjet1[99], isjet2[99]; //  isjet[99], the "real" size is given by the final value of jetidx since non-valid jets are replaced
  float btag[99]={0};
  float btag1[99]={0};
  float btag2[99]={0};
  TLorentzVector vaxis1, vaxis2, vjet[99];
  TLorentzVector V1, V2, V;
   
  int index_jetnomu[99];
  double pt_jetnomu[99];

  for (int i=0; i<tree_njet; i++) // Loop on jet
  {
    float jet_pt  = tree_jet_pt[i];
    float jet_eta = tree_jet_eta[i];
    float jet_phi = tree_jet_phi[i];
    V.SetPxPyPzE(tree_jet_px[i], tree_jet_py[i], tree_jet_pz[i], tree_jet_E[i]);
    
    // look if prompt lepton inside
    float deltaR1 = Deltar( jet_eta, jet_phi, lep1_eta, lep1_phi );
    float deltaR2 = Deltar( jet_eta, jet_phi, lep2_eta, lep2_phi );
    if ( deltaR1 < 0.4 ) {
      if ( jet_pt > lep1_pt ) V -= Vlep1;
      else V.SetPtEtaPhiM(lep1_pt/100., lep1_eta, lep1_phi, lep1_mass); // to be safe, this jet is ignored
      jet_pt  = V.Pt();
    }
    if ( deltaR2 < 0.4 ) {
      if ( jet_pt > lep2_pt ) V -= Vlep2;
      else V.SetPtEtaPhiM(lep2_pt/100., lep2_eta, lep2_phi, lep2_mass); // to be safe, this jet is ignored
      jet_pt  = V.Pt();
    }
    vjet[i] = V;
    pt_jetnomu[i] = jet_pt;
    isjet1[i] = false; // first neutralino jets
    isjet2[i] = false; // second neutralino jets

    btag[i]  = tree_jet_btag_DeepJet[i];
    btag1[i] = 0;
    btag2[i] = 0;    
  } // End Loop on jets

  // sort jets (again) by decreasing pT (but after having subtracted the prompt muons)
  TMath::Sort(tree_njet, pt_jetnomu, index_jetnomu);

  // jet seed
  if ( tree_njet > 0 ) {
    njet1 = 1;
    int jetseed = index_jetnomu[0];
    isjet1[jetseed] = true;
    vaxis1 = vjet[jetseed];
    btag1[jetseed] = btag[jetseed];
    float deltaR1 = Deltar( vaxis1.Eta(), vaxis1.Phi(), lep1_eta, lep1_phi );
    float deltaR2 = Deltar( vaxis1.Eta(), vaxis1.Phi(), lep2_eta, lep2_phi );
    if ( deltaR1 >= 0.4 && deltaR2 >= 0.4 ) njet1_nomu = 1;
  }

  /////////////////////////////////////////////////////////
  //-------------------------------------------------------
  // Event Axes
  //-------------------------------------------------------
  /////////////////////////////////////////////////////////

  float dR1 = 10., dR2 = 10.;
  float dRcut_hemis  = 1.5; // subjective choice default is 1.5
  float dRcut_tracks = 10; // no cut is better (could bias low track pT and high LLP ct) 
  int jet_noHemi = 0;
  if ( tree_njet > 1 ) {
    for (int ii=1; ii<tree_njet; ii++) // Loop on jet (but skip the seed)
  {
    int i = index_jetnomu[ii];
 
    // float jet_pt  = vjet[i].Pt();
    float jet_eta = vjet[i].Eta();
    float jet_phi = vjet[i].Phi();
    if ( njet1 > 0 ) dR1 = Deltar( jet_eta, jet_phi, vaxis1.Eta(), vaxis1.Phi() );
    if ( njet2 > 0 ) dR2 = Deltar( jet_eta, jet_phi, vaxis2.Eta(), vaxis2.Phi() );
    // test Paul (in case the conveners ask for three vertices)
    // In case Conveners ask for mroe than 2 vertices for BSM physics, one way we could expect other vertices could be jets not being part of the
    // hemispheres => back to back in eta and phi jets from the jets from the signal. This really should have a minor impact if not null
    // btu we never know  what BSM physics is made of :D .
    if (dR1>dRcut_hemis && dR2>dRcut_hemis )
      {
        jet_noHemi++;
      }

        float deltaR1 = Deltar( jet_eta, jet_phi, lep1_eta, lep1_phi );
    float deltaR2 = Deltar( jet_eta, jet_phi, lep2_eta, lep2_phi );
    // axis 1
    if ( njet1 > 0 && dR1 < dRcut_hemis ) {
      njet1++;
      vaxis1 += vjet[i];
      isjet1[i] = true;
      btag1[i]=btag[i];
      if (btag1[i]>0.7264 && ActivateBtagLog) {std::cout<<"Tight b jet ID"<<std::endl;}
      // 0.2770 : Medium 
      // 0.0494 : loose
      if ( deltaR1 >= 0.4 && deltaR2 >= 0.4 ) njet1_nomu++;
    }
    // axis 2
    if ( njet2 == 0 && !isjet1[i] ) {
      njet2 = 1;
      vaxis2 = vjet[i];
      isjet2[i] = true;
      btag2[i]=btag[i];
      if (btag2[i]>0.7264 && ActivateBtagLog) {std::cout<<"Tight b jet ID"<<std::endl;}
      if ( deltaR1 >= 0.4 && deltaR2 >= 0.4 ) njet2_nomu = 1;
    }
    else if ( njet2 > 0 && !isjet1[i] && !isjet2[i] && dR2 < dRcut_hemis ) {
      njet2++;
      vaxis2 += vjet[i];
      isjet2[i] = true;
      btag2[i]=btag[i];//njet2 insteag of i would also make sense, here there are voids when one of the other conditions above is filled
      if (btag2[i]>0.7264 && ActivateBtagLog) {std::cout<<"Tight b jet ID"<<std::endl;}
      if ( deltaR1 >= 0.4 && deltaR2 >= 0.4 ) njet2_nomu++;
    }
  }	  // end Loop on jet
}
  
  float axis1_eta = vaxis1.Eta();
  float axis1_phi = vaxis1.Phi();
  float axis2_eta = vaxis2.Eta();
  float axis2_phi = vaxis2.Phi();
  if ( njet2 == 0 )
  {  // compute an axis 2 even without jet, by taking the opposite in phi to axis 1
    axis2_eta = axis1_eta;
    axis2_phi = axis1_phi - 3.14159;
    if ( axis1_phi < 0 ) axis2_phi = axis1_phi + 3.14159;
    vaxis2.SetPtEtaPhiM(10., axis2_eta, axis2_phi, 0.);
  }
 
  // std::cout<<"Number of jets not belonging to the two hemispheres : "<<jet_noHemi<<std::endl;
//   // force the axes to the true LLP
//   vaxis1.SetPtEtaPhiM(LLP1_pt, LLP1_eta, LLP1_phi, neumass);
//   vaxis2.SetPtEtaPhiM(LLP2_pt, LLP2_eta, LLP2_phi, neumass);
//   axis1_eta = vaxis1.Eta();
//   axis1_phi = vaxis1.Phi();
//   axis2_eta = vaxis2.Eta();
//   axis2_phi = vaxis2.Phi();

  ///////////////////////////////
  // observed smuon 
  ///////////////////////////////
  
  // associate prompt muons and first hemisphere to get the highest overall pT
  TLorentzVector Vobs1, Vobs2; // retained smuon candidates 
  TLorentzVector VobsOp1, VobsOp2; // other candidate 
  float dRobs1, dRobs2, dRobsOp1, dRobsOp2; 
  V1 = Vlep1 + vaxis1;
  V2 = Vlep2 + vaxis1;
  if ( V1.Pt() > V2.Pt() ) {
    Vobs1 = V1;            
    Vobs2 = Vlep2 + vaxis2;
    VobsOp1 = V2;
    VobsOp2 = Vlep1 + vaxis2;
    dRobs1 = Deltar( lep1_eta, lep1_phi, axis1_eta, axis1_phi );
    dRobs2 = Deltar( lep2_eta, lep2_phi, axis2_eta, axis2_phi );
    dRobsOp1 = Deltar( lep2_eta, lep2_phi, axis1_eta, axis1_phi );
    dRobsOp2 = Deltar( lep1_eta, lep1_phi, axis2_eta, axis2_phi );
  }
  else {
    Vobs1 = V2;
    Vobs2 = Vlep1 + vaxis2;
    VobsOp1 = V1;            
    VobsOp2 = Vlep2 + vaxis2;
    dRobs1 = Deltar( lep2_eta, lep2_phi, axis1_eta, axis1_phi );
    dRobs2 = Deltar( lep1_eta, lep1_phi, axis2_eta, axis2_phi );
    dRobsOp1 = Deltar( lep1_eta, lep1_phi, axis1_eta, axis1_phi );
    dRobsOp2 = Deltar( lep2_eta, lep2_phi, axis2_eta, axis2_phi );
  }

  mva_Evts_Hemi1_njet_nomu = njet1_nomu;
  mva_Evts_Hemi2_njet_nomu = njet2_nomu;
  mva_Evts_Hemi1_pt = vaxis1.Pt(); 
  mva_Evts_Hemi2_pt = vaxis2.Pt();
  mva_Evts_Hemi1_eta = axis1_eta; 
  mva_Evts_Hemi2_eta = axis2_eta;
  mva_Evts_Hemi1_phi = axis1_phi;
  mva_Evts_Hemi2_phi = axis2_phi;

  mva_Evts_Hemi1_Mass = TMath::Max(vaxis1.Mag(),0.);
  mva_Evts_Hemi2_Mass = TMath::Max(vaxis2.Mag(),0.);

  ///////////////////////////////
  // compare with neutralino axis
  ///////////////////////////////

  float axis1_dR, axis2_dR, dR_axis12;
  int iLLPrec1 = 1, iLLPrec2 = 2;
  if ( isMC_ && nLLP == 2 )
  {
    dR1 = Deltar( axis1_eta, axis1_phi, LLP1_eta, LLP1_phi ); //dR between reco axis of jets and gen neutralino
    dR2 = Deltar( axis1_eta, axis1_phi, LLP2_eta, LLP2_phi );
    axis1_dR = dR1;
    if ( dR2 < dR1 ) {
      iLLPrec1 = 2;
      axis1_dR = dR2;
    }
    dR1 = Deltar( axis2_eta, axis2_phi, LLP1_eta, LLP1_phi );
    dR2 = Deltar( axis2_eta, axis2_phi, LLP2_eta, LLP2_phi );
    axis2_dR = dR2;
    if ( dR1 < dR2 ) { 
      iLLPrec2 = 1;
      axis2_dR = dR1;
    }
    if ( axis1_dR < axis2_dR ) {
      if ( iLLPrec1 == 1 ) {
        iLLPrec2 = 2;
	      axis2_dR = dR2; 
      }
      else {
        iLLPrec2 = 1;
	      axis2_dR = dR1; 
      }
    }
    else {
    if ( iLLPrec2 == 1 ) {
        iLLPrec1 = 2;
        axis2_dR = dR1;
        axis1_dR = Deltar( axis1_eta, axis1_phi, LLP2_eta, LLP2_phi );
    }
    else {
        iLLPrec1 = 1;
        axis2_dR = dR2;
        axis1_dR = Deltar( axis1_eta, axis1_phi, LLP1_eta, LLP1_phi );
      }
    }
  }
  dR_axis12 = Deltar(axis1_eta,axis1_phi,axis2_eta,axis2_phi);

  // cout << " njet1 " << njet1 << " and njet2" << njet2 << endl;
  // cout << " axis1_eta " << axis1_eta << " and axis2_eta" << axis2_eta << endl;
  // cout << " axis1_phi " << axis1_phi << " and axis2_phi" << axis2_phi << endl;
  // cout << " axis1_dR " << axis1_dR << " and axis2_dR" << axis2_dR << endl;
  // cout << " dR_axis12 " << dR_axis12 << endl;
  
    ///////////////////////////////////////////////////////
    // Compare between the gen neutralino and the reco axis from gen jets and with the reco axis with reco jets
    ///////////////////////////////////////////////////////
    if ( isMC_ )
    {
      float dRGenAxisNeuMin = -1000;
      float dRGenAxisNeuMax = -1000;
      dRGenAxisNeuMin = Deltar(Gen_neu1_eta,Gen_neu1_phi,Genvaxis1.Eta(),Genvaxis1.Phi() );
      dRGenAxisNeuMax = Deltar(Gen_neu1_eta,Gen_neu1_phi,Genvaxis2.Eta(),Genvaxis2.Phi() );
      if(dRGenAxisNeuMin > dRGenAxisNeuMax )
      {
	      dRGenAxisNeuMin = Deltar(Gen_neu1_eta,Gen_neu1_phi,Genvaxis2.Eta(),Genvaxis2.Phi() );
	      dRGenAxisNeuMax = Deltar(Gen_neu1_eta,Gen_neu1_phi,Genvaxis1.Eta(),Genvaxis1.Phi() );
      }
      tree_GenAxis_Neu_dRmin.push_back(dRGenAxisNeuMin);
      tree_GenAxis_Neu_dRmax.push_back(dRGenAxisNeuMax);

      dRGenAxisNeuMin = Deltar(Gen_neu2_eta,Gen_neu2_phi,Genvaxis2.Eta(),Genvaxis2.Phi());
      dRGenAxisNeuMax = Deltar(Gen_neu2_eta,Gen_neu2_phi,Genvaxis1.Eta(),Genvaxis1.Phi());
      if(dRGenAxisNeuMin > dRGenAxisNeuMax )
      {
	      dRGenAxisNeuMin = Deltar(Gen_neu2_eta,Gen_neu2_phi,Genvaxis1.Eta(),Genvaxis1.Phi());
	      dRGenAxisNeuMax = Deltar(Gen_neu2_eta,Gen_neu2_phi,Genvaxis2.Eta(),Genvaxis2.Phi() );
      }
      tree_GenAxis_Neu_dRmin.push_back(dRGenAxisNeuMin);
      tree_GenAxis_Neu_dRmax.push_back(dRGenAxisNeuMax);

      //-------------------------:/
      float dRGenAxisRecoAxisMin = Deltar(vaxis1.Eta(),vaxis1.Phi(),Genvaxis1.Eta(),Genvaxis1.Phi() );
      float dRGenAxisRecoAxisMax = Deltar(vaxis1.Eta(),vaxis1.Phi(),Genvaxis2.Eta(),Genvaxis2.Phi() );
      if(dRGenAxisRecoAxisMin > dRGenAxisRecoAxisMax )
      {
	      dRGenAxisRecoAxisMin = Deltar(vaxis1.Eta(),vaxis1.Phi(),Genvaxis2.Eta(),Genvaxis2.Phi() );
	      dRGenAxisRecoAxisMax = Deltar(vaxis1.Eta(),vaxis1.Phi(),Genvaxis1.Eta(),Genvaxis1.Phi() );
      }
      tree_GenAxis_RecoAxis_dRmin.push_back(dRGenAxisRecoAxisMin);
      tree_GenAxis_RecoAxis_dRmax.push_back(dRGenAxisRecoAxisMax);

      dRGenAxisRecoAxisMin = Deltar(vaxis2.Eta(),vaxis2.Phi(),Genvaxis2.Eta(),Genvaxis2.Phi());
      dRGenAxisRecoAxisMax = Deltar(vaxis2.Eta(),vaxis2.Phi(),Genvaxis1.Eta(),Genvaxis1.Phi());
      if(dRGenAxisRecoAxisMin > dRGenAxisRecoAxisMax )
      {
	      dRGenAxisRecoAxisMin = Deltar(vaxis2.Eta(),vaxis2.Phi(),Genvaxis1.Eta(),Genvaxis1.Phi());
	      dRGenAxisRecoAxisMax = Deltar(vaxis2.Eta(),vaxis2.Phi(),Genvaxis2.Eta(),Genvaxis2.Phi());
      }
      tree_GenAxis_RecoAxis_dRmin.push_back(dRGenAxisRecoAxisMin);
      tree_GenAxis_RecoAxis_dRmax.push_back(dRGenAxisRecoAxisMax);
    } // endif MC
    
      bool LooseAxesStatus = false ;
      bool MediumAxesStatus = false ;
      bool TightAxesStatus = false ;
   if ( ActivateBtag )
    {
      for (int indice = 0 ; indice < 99 ; indice++)    // Loop on jet
        {
          // std::cout<<"btag1 et btag2 : "<<btag1[indice]<<" / "<<btag2[indice]<<std::endl;
          // if (btag1[indice]<0 ||btag1[indice]>1 ) continue ;
          if (btag1[indice] >LooseWP ) //request at least a b jet in each hemisphere to then continue to potentially veto TTbar
            {
              for (int indice2 = 0 ; indice2 <99 ; indice2++)
                {
                  if (btag2[indice2]>LooseWP)
                    {
                      LooseAxesStatus = true;
                    }
                } 
            }

          if (btag1[indice] >MediumWP) //request at least a b jet in each hemisphere to then continue to potentially veto TTbar
            {
              for (int indice2 = 0 ; indice2 <99 ; indice2++)
                {
                  if ( btag2[indice2]>MediumWP)
                    {
                      MediumAxesStatus = true;
                    }
                }          
            }
          if (btag1[indice] >TightWP) //request at least a b jet in each hemisphere to then continue to potentially veto TTbar
            {
              for (int indice2 = 0 ; indice2 <99 ; indice2++)
                {
                  if ( btag2[indice2]>TightWP)
                    {
                      TightAxesStatus = true;
                    }
                } 
            }
        }
        // tree_Hemi_LooseBTag_axes.push_back(LooseAxesStatus);
        tree_Hemi_LooseBTag_axes.push_back(LooseAxesStatus);
        // tree_Hemi_MediumBTag_axes.push_back(MediumAxesStatus);
        tree_Hemi_MediumBTag_axes.push_back(MediumAxesStatus);
        // tree_Hemi_TightBTag_axes.push_back(TightAxesStatus);
        tree_Hemi_TightBTag_axes.push_back(TightAxesStatus);
    }

    mva_Evts_LooseAxes	    = LooseAxesStatus;
    mva_Evts_MediumAxes	    = MediumAxesStatus;
    mva_Evts_TightAxes	    = TightAxesStatus;


    //////////////////////////////////
    //////////////////////////////////
    //////////   Tracks   ////////////
    //////////////////////////////////
    //////////////////////////////////

    vector <pat::PackedCandidateRef> MINIgeneralTracks;
    unsigned int TRACK_SIZE = pc->size();
    for (unsigned int ipc = 0; ipc < pc->size(); ipc++) 
    {
      MINIgeneralTracks.push_back(pat::PackedCandidateRef(pcs, ipc));
    }
    if ( IncludeLostTrack )
      {
        for (unsigned int k=0; k<lostpc->size();k++)
          {
            MINIgeneralTracks.push_back(pat::PackedCandidateRef(lostpcs, k));  
          }
          TRACK_SIZE = TRACK_SIZE+lostpc->size();
      }   
    tree_TRACK_SIZE = TRACK_SIZE;
    mva_Evts_nTrks = tree_TRACK_SIZE;

    //-------------------------------------//
   
    // track variables declaration //
    float tk_nHit = -1000  ;
    float tk_charge = -1000;
    float tk_pt = -1000;
    float tk_eta= -1000 ;
    float tk_phi= -1000 ;
    float tk_NChi2= -1000 ;
    float tk_drSig = -1000;
    float tk_dxy = -1000;
    float tk_dxyError= -1000 ;
    float tk_dz = -1000;
    float tk_dzError= -1000;
    float tk_dzSig= -1000;
    float tk_vx= -1000 ;
    float tk_vy = -1000;
    float tk_vz = -1000;
    float tk_px = -1000;
    float tk_py = -1000;
    float tk_pz = -1000;
    float tk_e  = -1000;
    HitPattern tk_HitPattern;

    //---------------------------------------------------------------------------------//
    //----------------------- V0 Producer adapted in MiniAOD---------------------------//
    //---------------------------------------------------------------------------------//

    // - Be careful => tree_V0_xxx comes from the collection of CMSSW
    // the following is the V0Producer code adapted in MiniAOD that will be also used to
    // reconstruct secondary interactiosn (photon conversions and nuclear intractions) with slight modification
    // Basically, only cuts and the badtkhit have been changed from the V0Producer code.
    // Keep in mind that we are using the MINIGeneral Tracks (packedPFCandidate and LostTrack)  != cmssw collection with a different pt cut

    const double piMass = 0.13957018;
    const double piMassSquared = piMass*piMass;
    const double protonMass = 0.938272046;
    const double protonMassSquared = protonMass*protonMass;
    const double kShortMass = 0.497614;
    const double lambdaMass = 1.115683;

    bool useBS_ = false;       // false by default; if true use beamspot instead of primary vertex
    bool vertexFitter_ = true; // true  by default: Kalman & False : AVF
    bool useRefTracks_ = true; // true  by default: # use the refitted tracks returned from the KVF for V0Candidate kinematics, # this is automatically set to False if using the AdaptiveVertexFitter

    //V0 to be reconstructed
    bool doKShorts_ = true;
    bool doLambdas_ = true;
    // # Track normalized Chi2 <
    float tkChi2Cut_    = 10.; // 10. by default
    // # Number of valid hits on track >=
    int tkNHitsCut_   = 3;   // 3 by default
    // # Pt of track >
    float tkPtCut_     = 0.35; // 0.35 by default
    // # Track impact parameter significance >
    float tkIPSigXYCut_ = 2.;  // 2. by default
    float tkIPSigZCut_  = -1.;

    //  # -- cuts on the vertex --
    // # Vertex chi2 <
    float vtxChi2Cut_ = 6.63;
    // # XY decay distance significance >
    float vtxDecaySigXYCut_ = 15.;
    // # XYZ decay distance significance >
    float vtxDecaySigXYZCut_ = -1.;
    
    // # -- miscellaneous cuts --
    // # POCA distance between tracks <
    float tkDCACut_ = 1.; // 1. by default
    // # invariant mass of track pair - assuming both tracks are charged pions <
    float mPiPiCut_ = 0.6; // 0.6 by default
    // # check if either track has a hit radially inside the vertex position minus this number times the sigma of the vertex fit
    // # note: Set this to -1 to disable this cut, which MUST be done if you want to run V0Producer on the AOD track collection!
    float innerHitPosCut_ = 4.;
    // # cos(angleXY) between x and p of V0 candidate >
    float cosThetaXYCut_ = 0.998;
    // # cos(angleXYZ) between x and p of V0 candidate >
    float cosThetaXYZCut_ = -2.;
    
    // # -- cuts on the V0 candidate mass --
    // # V0 mass window +- pdg value
 
    float kShortMassCut_ = 0.07;  // 0.07 by default
    float lambdaMassCut_ = 0.05;  // 0.05 by default
    float kShortMassSel_ = 0.022;  // track selection cut 0.476 - 0.520
    float lambdaMassSel_ = 0.0060; // track selection cut /: could go for asymmetric selectino
 
    int temp_nV0_reco = 0;

    //------------------Selection of good tracks for the vertexing---------------------//

    // fill vectors of TransientTracks and TrackRefs after applying preselection cuts

    std::vector<reco::Track> theTrackRefs;
    std::vector<reco::TransientTrack> theTransTracks;

    std::vector<std::pair<bool,float>> idxMGT; // contains the following informations : the size of this vector is the number of selected tracks passing the cuts
    // The bool says "is used at the end to build a V0 candidate", the float keeps in memory the index of the track from MINIGeneralTracks that is used to build the V0 candiates
    std::vector<std::pair<bool,float>> idxSecIntMGT;
    // std::vector<std::pair<bool,float>> idxMuonMGT;

    for (unsigned int ipc = 0; ipc < TRACK_SIZE; ipc++) { // loop on all packedPFCandidates + lostTrackss
      pat::PackedCandidateRef pcref = MINIgeneralTracks[ipc];
      const reco::Track *trackPcPtr = pcref->bestTrack();
      if ( !trackPcPtr ) continue;

      reco::Track tk = *trackPcPtr;
      if ( tk.hitPattern().numberOfValidHits() == 0 ) continue;
      if ( tk.charge() == 0 ) continue;
      double ipsigXY = std::abs(tk.dxy(PV.position()) / tk.dxyError());
      if ( useBS_ ) ipsigXY = std::abs(tk.dxy(bs) / tk.dxyError());
      double ipsigZ = std::abs(tk.dz(PV.position()) / tk.dzError());
      if ( tk.normalizedChi2() < tkChi2Cut_ && tk.hitPattern().numberOfValidHits() >= tkNHitsCut_ &&
           tk.pt() > tkPtCut_ && ipsigXY > tkIPSigXYCut_ && ipsigZ > tkIPSigZCut_ ) 
      {
        reco::Track* tmpRef = &tk;
        idxMGT.push_back(make_pair(false,ipc));
        idxSecIntMGT.push_back(make_pair(false,ipc));
        theTrackRefs.push_back(std::move(*tmpRef));
        reco::TransientTrack tmpTransient(*tmpRef,theMagneticField);
        theTransTracks.push_back(std::move(tmpTransient));
      }
    }
    // good tracks have now been selected for vertexing

//------------------------------ loop over tracks and vertex good charged track pairs
    for (unsigned int trd1 = 0; trd1 < theTrackRefs.size()-1; ++trd1) 
    {
      for (unsigned int trd2 = trd1+1; trd2 < theTrackRefs.size(); ++trd2) 
      {
        reco::Track positiveTrackRef;
        reco::Track negativeTrackRef;
        reco::TransientTrack* posTransTkPtr = nullptr;
        reco::TransientTrack* negTransTkPtr = nullptr;
        int Selec = 0;
        float drsig1 = std::abs(theTrackRefs[trd1].dxy(PV.position())/theTrackRefs[trd1].dxyError());
        float drsig2 = std::abs(theTrackRefs[trd2].dxy(PV.position())/theTrackRefs[trd2].dxyError());
        if ( useBS_ ) {
          drsig1 = std::abs(theTrackRefs[trd1].dxy(bs)/theTrackRefs[trd1].dxyError());
          drsig2 = std::abs(theTrackRefs[trd2].dxy(bs)/theTrackRefs[trd2].dxyError());
	}
 
        if ( !(theTrackRefs[trd1].pt() > pt_Cut && theTrackRefs[trd1].normalizedChi2() < NChi2_Cut && drsig1 > drSig_Cut) 
	  && !(theTrackRefs[trd2].pt() > pt_Cut && theTrackRefs[trd2].normalizedChi2() < NChi2_Cut && drsig2 > drSig_Cut) ) continue;
 
        if ( theTrackRefs[trd1].charge() < 0. && theTrackRefs[trd2].charge() > 0. ) 
        {
          negativeTrackRef = theTrackRefs[trd1];
          positiveTrackRef = theTrackRefs[trd2];
          negTransTkPtr = &theTransTracks[trd1];
          posTransTkPtr = &theTransTracks[trd2];
          Selec = 1;
        } 
        else if ( theTrackRefs[trd1].charge() > 0. && theTrackRefs[trd2].charge() < 0. )
        {
          negativeTrackRef = theTrackRefs[trd2];
          positiveTrackRef = theTrackRefs[trd1];
          negTransTkPtr = &theTransTracks[trd2];
          posTransTkPtr = &theTransTracks[trd1];
          Selec = 2;
        } 
        if ( Selec == 0 ) continue;

        // ---------------measure distance between tracks at their closest approach---------------
        // these two variables are needed to 'pin' the temporary value returned to the stack
        // in order to keep posState and negState from pointing to destructed objects
        auto const& posImpact = posTransTkPtr->impactPointTSCP();
        auto const& negImpact = negTransTkPtr->impactPointTSCP();
        if ( !posImpact.isValid() || !negImpact.isValid() ) continue;

        FreeTrajectoryState const & posState = posImpact.theState();
        FreeTrajectoryState const & negState = negImpact.theState();
        ClosestApproachInRPhi cApp;
        cApp.calculate(posState, negState);
        if ( !cApp.status() ) continue;

        float dca = std::abs(cApp.distance());
        if ( dca > tkDCACut_ ) continue;

        // the POCA should at least be in the sensitive volume
        GlobalPoint cxPt = cApp.crossingPoint();
        if ( sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || std::abs(cxPt.z()) > 300. ) continue;

        // the tracks should at least point in the same quadrant
        TrajectoryStateClosestToPoint const & posTSCP = posTransTkPtr->trajectoryStateClosestToPoint(cxPt);
        TrajectoryStateClosestToPoint const & negTSCP = negTransTkPtr->trajectoryStateClosestToPoint(cxPt);
        if ( !posTSCP.isValid() || !negTSCP.isValid() ) continue;
        if ( posTSCP.momentum().dot(negTSCP.momentum()) < 0. ) continue;
        	  
        // calculate mPiPi
        double totalE = sqrt(posTSCP.momentum().mag2() + piMassSquared) + sqrt(negTSCP.momentum().mag2() + piMassSquared);
        double totalESq = totalE*totalE;
        double totalPSq = (posTSCP.momentum() + negTSCP.momentum()).mag2();
        double mass = TMath::Sqrt(totalESq - totalPSq);
        if ( mass > mPiPiCut_ ) continue; // mPiPi < 0.6 GeV (true for Lambda ?)
 
        // Fill the vector of TransientTracks to send to KVF
        std::vector<reco::TransientTrack> transTracks;
        transTracks.reserve(2);
        transTracks.push_back(*posTransTkPtr);
        transTracks.push_back(*negTransTkPtr);
        
        // create the vertex fitter object and vertex the tracks
        TransientVertex theRecoVertex;
        if ( vertexFitter_ ) // true is recommended, AVF is more likely to be better for high number of tracks vertices
        {
          KalmanVertexFitter theKalmanFitter(useRefTracks_ == 0 ? false : true); 
          // # use the refitted tracks returned from the KVF for V0Candidate kinematics
          // # this is automatically set to False if using the AdaptiveVertexFitter
          theRecoVertex = theKalmanFitter.vertex(transTracks);
        } 
        else if ( !vertexFitter_ ) 
        {
          useRefTracks_ = false;
          AdaptiveVertexFitter theAdaptiveFitter;
          theRecoVertex = theAdaptiveFitter.vertex(transTracks);
        }
        if ( !theRecoVertex.isValid() ) continue;

        reco::Vertex theVtx = theRecoVertex;
        if ( theVtx.normalizedChi2() > vtxChi2Cut_ ) continue; // Vtx NChi2 < 6.63

        GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());
        // 2D decay length significance 
	float distsigXY = 10000.;
        SMatrixSym3D  totalCov = PV.covariance() + theVtx.covariance();
        if ( useBS_ ) totalCov = bs.rotatedCovariance3D() + theVtx.covariance();
        SVector3 distVecXY(vtxPos.x()-PV.x(), vtxPos.y()-PV.y(), 0.);
        double distMagXY = ROOT::Math::Mag(distVecXY);
        double sigmaDistMagXY = sqrt(ROOT::Math::Similarity(totalCov, distVecXY)) / distMagXY;
	if ( sigmaDistMagXY > 0. ) distsigXY = distMagXY / sigmaDistMagXY;
        if ( distsigXY < vtxDecaySigXYCut_ ) continue; // Vtx XY signif > 15.
       	     
        // 3D decay significance
	float distsigXYZ = 10000.;
        SVector3 distVecXYZ(vtxPos.x()-PV.x(), vtxPos.y()-PV.y(), vtxPos.z()-PV.z());
        double distMagXYZ = ROOT::Math::Mag(distVecXYZ);
        double sigmaDistMagXYZ = sqrt(ROOT::Math::Similarity(totalCov, distVecXYZ)) / distMagXYZ;
	if ( sigmaDistMagXYZ > 0. ) distsigXYZ = distMagXYZ / sigmaDistMagXYZ;
        if ( vtxDecaySigXYZCut_ > 0. && distsigXYZ < vtxDecaySigXYZCut_ ) continue; // no cut on Vtx XYZ significance
       	     
        // z significance
	float distsigZ = 10000.;
        SVector3 distVecZ(0, 0, vtxPos.z()-PV.z());
        double distMagZ = TMath::Abs(vtxPos.z()-PV.z());
        double sigmaDistMagZ = sqrt(ROOT::Math::Similarity(theVtx.covariance(), distVecZ)) / distMagZ;
	if ( sigmaDistMagZ > 0. ) distsigZ = distMagZ / sigmaDistMagZ;

        // make sure the vertex radius is within the inner track hit radius
        // The original code is in RECO datatier using .innerOk and innerPosition => asking for the position of the inner hit
        // In miniaod, we use the implemented class PropaHitpattern => introducing a bit more of uncertainty on the position of the firsthit
        // < 1 cm in barrel and ~1-2cm in disks
        bool badTkHit = false;
        if ( innerHitPosCut_ > 0. ) // && positiveTrackRef.innerOk() //innerOk is trackextra => not availeble in MINIaod
        {
          int trdx = trd1;
          if ( Selec == 1 ) trdx = trd2;
          if ( idxMGT[trdx].second < pc->size() ) // not for the lost tracks, but as their default first hit is PIXBL1 it would not hurt...
          {
            const HitPattern hp = positiveTrackRef.hitPattern();//tk_HitPattern;
            uint16_t firsthit = hp.getHitPattern(HitPattern::HitCategory::TRACK_HITS,0);     
             //---Creating State to propagate from  TT---//
            reco::TransientTrack TTpos = theTransientTrackBuilder->build(positiveTrackRef);
            // reco::TransientTrack TT (positiveTrackRef,posBestTracks[countpos].field());
            GlobalPoint vert (positiveTrackRef.vx(),positiveTrackRef.vy(),positiveTrackRef.vz()); // Point where the propagation will start (Reference Point)
            const TrajectoryStateOnSurface Surtraj = TTpos.stateOnSurface(vert); // TSOS of this point
            const MagneticField* B = TTpos.field(); // 3.8T
            AnalyticalPropagator* Prop = new AnalyticalPropagator(B);
            Basic3DVector<float> P3D2(positiveTrackRef.vx(),positiveTrackRef.vy(),positiveTrackRef.vz());  // global frame
            Basic3DVector<float> B3DV (positiveTrackRef.px(),positiveTrackRef.py(),positiveTrackRef.pz()); // global frame 
            float Eta = positiveTrackRef.eta();
            float Phi = positiveTrackRef.phi();
            float vz  = positiveTrackRef.vz();
            
            std::pair<int,GloballyPositioned<float>::PositionType> posFHPosition = posPHP->Main(firsthit,Prop,Surtraj,Eta,Phi,vz,P3D2,B3DV);
            float xF = posFHPosition.second.x() - PV.x();
            float yF = posFHPosition.second.y() - PV.y();
            float rF = TMath::Sqrt( xF*xF + yF*yF );
                  float zF = TMath::Abs( posFHPosition.second.z() - PV.z() );
            // PXB
            if ( rF > 2.7 && rF < 16.5 && zF < 27.0 ) {
              if ( rF <  distMagXY - 0.40 ) badTkHit = true;
            }
            // TIB
            if ( rF > 23.5 && rF < 36.2 && zF < 66.0 ) {
              if ( rF <  distMagXY - 0.84 ) badTkHit = true;
            }
            if ( rF > 40.0 && rF < 52.0 && zF < 66.0 ) {
              if ( rF <  distMagXY - 2.24 ) badTkHit = true;
            }
            // TOB
            if ( rF > 58.4 && rF < 62.9 && zF < 107. ) {
              if ( rF <  distMagXY - 2.8 ) badTkHit = true;
            }
            // PXF
            if ( rF > 4.5 && rF < 16.1 && zF > 29.6 && zF < 52.0 ) {
              if ( zF <  distMagZ - 3.2 ) badTkHit = true;
            }
            // TID
            if ( rF > 22.5 && rF < 50.5 && zF > 74.3 && zF < 109.7 ) {
              if ( zF <  distMagZ - 3.9 ) badTkHit = true;
            }
            // TEC
            if ( rF > 22.0 && rF < 70.0 && zF > 126.4 && zF < 193.7 ) {
              if ( zF <  distMagZ - 6.7 ) badTkHit = true;
            }

          }
        }

        if ( innerHitPosCut_ > 0. && !badTkHit ) // && negativeTrackRef.innerOk()
        {
          int trdx = trd2;
          if ( Selec == 1 ) trdx = trd1;
          if ( idxMGT[trdx].second < pc->size() ) // not for the lost tracks, but as their default first hit is PIXBL1 it would not hurt...
          {
            const HitPattern hp = negativeTrackRef.hitPattern();
            uint16_t firsthit = hp.getHitPattern(HitPattern::HitCategory::TRACK_HITS,0); 
            //---Creating State to propagate from  TT---//
            reco::TransientTrack  TTneg = theTransientTrackBuilder->build(negativeTrackRef);
            GlobalPoint vert (negativeTrackRef.vx(),negativeTrackRef.vy(),negativeTrackRef.vz()); // Point where the propagation will start (Reference Point)
            const TrajectoryStateOnSurface Surtraj = TTneg.stateOnSurface(vert); // TSOS of this point
            const MagneticField* B = TTneg.field(); // 3.8T
            AnalyticalPropagator* Prop = new AnalyticalPropagator(B);
            Basic3DVector<float> P3D2(negativeTrackRef.vx(),negativeTrackRef.vy(),negativeTrackRef.vz());  // global frame
            Basic3DVector<float> B3DV(negativeTrackRef.px(),negativeTrackRef.py(),negativeTrackRef.pz()); // global frame 
            float Eta = negativeTrackRef.eta();
            float Phi = negativeTrackRef.phi();
            float vz  = negativeTrackRef.vz();
            
            std::pair<int,GloballyPositioned<float>::PositionType> negFHPosition = negPHP->Main(firsthit,Prop,Surtraj,Eta,Phi,vz,P3D2,B3DV);
            float xF = negFHPosition.second.x() - PV.x();
            float yF = negFHPosition.second.y() - PV.y();
            float rF = TMath::Sqrt( xF*xF + yF*yF );
                  float zF = TMath::Abs( negFHPosition.second.z() - PV.z() );
            // PXB
            if ( rF > 2.7 && rF < 16.5 && zF < 27.0 ) {
              if ( rF <  distMagXY - 0.40 ) badTkHit = true;
            }
            // TIB
            if ( rF > 23.5 && rF < 36.2 && zF < 66.0 ) {
              if ( rF <  distMagXY - 0.84 ) badTkHit = true;
            }
            if ( rF > 40.0 && rF < 52.0 && zF < 66.0 ) {
              if ( rF <  distMagXY - 2.24 ) badTkHit = true;
            }
            // TOB
            if ( rF > 58.4 && rF < 62.9 && zF < 107. ) {
              if ( rF <  distMagXY - 2.8 ) badTkHit = true;
            }
            // PXF
            if ( rF > 4.5 && rF < 16.1 && zF > 29.6 && zF < 52.0 ) {
              if ( zF <  distMagZ - 3.2 ) badTkHit = true;
            }
            // TID
            if ( rF > 22.5 && rF < 50.5 && zF > 74.3 && zF < 109.7 ) {
              if ( zF <  distMagZ - 3.9 ) badTkHit = true;
            }
            // TEC
            if ( rF > 22.0 && rF < 70.0 && zF > 126.4 && zF < 193.7 ) {
              if ( zF <  distMagZ - 6.7 ) badTkHit = true;
            }

          }
        }
        if ( badTkHit && dca > 0.1 ) continue;
 
        std::unique_ptr<TrajectoryStateClosestToPoint> trajPlus;
        std::unique_ptr<TrajectoryStateClosestToPoint> trajMins;
        std::vector<reco::TransientTrack> theRefTracks;
        if ( theRecoVertex.hasRefittedTracks() ) {
          theRefTracks = theRecoVertex.refittedTracks();
        }
                  
        if ( useRefTracks_ && theRefTracks.size() > 1 ) 
        {
          reco::TransientTrack* thePositiveRefTrack = nullptr;
          reco::TransientTrack* theNegativeRefTrack = nullptr;
          for (std::vector<reco::TransientTrack>::iterator iTrack = theRefTracks.begin(); iTrack != theRefTracks.end(); ++iTrack) 
          {
            if ( iTrack->track().charge() > 0. ) thePositiveRefTrack = &*iTrack;
            else  			   theNegativeRefTrack = &*iTrack;
          }
          if ( thePositiveRefTrack == nullptr || theNegativeRefTrack == nullptr ) continue;
          trajPlus.reset(new TrajectoryStateClosestToPoint(thePositiveRefTrack->trajectoryStateClosestToPoint(vtxPos)));
          trajMins.reset(new TrajectoryStateClosestToPoint(theNegativeRefTrack->trajectoryStateClosestToPoint(vtxPos)));
        } 
        else 
        {
          trajPlus.reset(new TrajectoryStateClosestToPoint(posTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
          trajMins.reset(new TrajectoryStateClosestToPoint(negTransTkPtr->trajectoryStateClosestToPoint(vtxPos)));
        }
        if ( trajPlus.get() == nullptr || trajMins.get() == nullptr || !trajPlus->isValid() || !trajMins->isValid() ) continue;
                   
        GlobalVector positiveP(trajPlus->momentum());
        GlobalVector negativeP(trajMins->momentum());
        GlobalVector totalP(positiveP + negativeP);
     
        // 2D pointing angle
        double dx = theVtx.x()-PV.x();
        double dy = theVtx.y()-PV.y();
        double px = totalP.x();
        double py = totalP.y();
        float  angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
        if ( angleXY < cosThetaXYCut_ ) continue; // cos( dxy , pxy ) > 0.998
     
        // 3D pointing angle
        float angleXYZ = -2.;
        // double dz = theVtx.z()-PV.z();
        dz = theVtx.z()-PV.z();
        double pz = totalP.z();
        angleXYZ = (dx*px+dy*py+dz*pz)/(sqrt(dx*dx+dy*dy+dz*dz)*sqrt(px*px+py*py+pz*pz));
        if ( angleXYZ < cosThetaXYZCut_ ) continue; // not cut by default

        // delta eta pointing
        double dr = sqrt(dx*dx + dy*dy);
        double etaZ = 0.;
        if ( dz != 0. ) etaZ = -TMath::Log(abs(TMath::Tan(TMath::ATan(dr/dz)/2.)));
        if ( dz < 0. )  etaZ = -etaZ;
        double ptot = TMath::Sqrt(px*px + py*py + pz*pz);
        double etaP = 0.;
        if ( abs(pz) < ptot ) etaP = 0.5 * TMath::Log((ptot+pz)/(ptot-pz));
	      float detaPointing = etaP - etaZ;
                   
        // calculate total energy of V0 3 ways: assume it's a Kshort, a Lambda, or a LambdaBar.
        double piPlusE = sqrt(positiveP.mag2() + piMassSquared);
        double piMinusE = sqrt(negativeP.mag2() + piMassSquared);
        double protonE = sqrt(positiveP.mag2() + protonMassSquared);
        double antiProtonE = sqrt(negativeP.mag2() + protonMassSquared);
        double kShortETot = piPlusE + piMinusE;
        double lambdaEtot = protonE + piMinusE;
        double lambdaBarEtot = antiProtonE + piPlusE;
     
        // Create momentum 4-vectors for the 3 candidate types
        const reco::Particle::LorentzVector kShortP4(totalP.x(), totalP.y(), totalP.z(), kShortETot);
        const reco::Particle::LorentzVector lambdaP4(totalP.x(), totalP.y(), totalP.z(), lambdaEtot);
        const reco::Particle::LorentzVector lambdaBarP4(totalP.x(), totalP.y(), totalP.z(), lambdaBarEtot);
     
        reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
        const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
        double vtxChi2(theVtx.chi2());
        double vtxNdof(theVtx.ndof());
      
        // Create the VertexCompositeCandidate object that will be stored in the Event
        reco::VertexCompositeCandidate* theKshort = nullptr;
        reco::VertexCompositeCandidate* theLambda = nullptr;
        reco::VertexCompositeCandidate* theLambdaBar = nullptr;
     
        if ( doKShorts_ ) {
          theKshort = new reco::VertexCompositeCandidate(0, kShortP4, vtx, vtxCov, vtxChi2, vtxNdof);
        }
        if ( doLambdas_ ) {
          if ( positiveP.mag2() > negativeP.mag2() ) {
            theLambda = new reco::VertexCompositeCandidate(0, lambdaP4, vtx, vtxCov, vtxChi2, vtxNdof);
          } 
	        else {
            theLambdaBar = new reco::VertexCompositeCandidate(0, lambdaBarP4, vtx, vtxCov, vtxChi2, vtxNdof);
          }
        }
     
        // Create daughter candidates for the VertexCompositeCandidates
        reco::RecoChargedCandidate thePiPlusCand(
           1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), piPlusE), vtx);
        reco::TrackRef posTrackRef = thePiPlusCand.track();
        thePiPlusCand.setTrack(posTrackRef); // positiveTrackRef
        
        reco::RecoChargedCandidate thePiMinusCand(
          -1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), piMinusE), vtx);
        reco::TrackRef negTrackRef = thePiMinusCand.track();
        thePiMinusCand.setTrack(negTrackRef); // negativeTrackRef
        
        reco::RecoChargedCandidate theProtonCand(
           1, reco::Particle::LorentzVector(positiveP.x(), positiveP.y(), positiveP.z(), protonE), vtx);
        reco::TrackRef pos2TrackRef = theProtonCand.track();
        theProtonCand.setTrack(pos2TrackRef); // positiveTrackRef
     
        reco::RecoChargedCandidate theAntiProtonCand(
          -1, reco::Particle::LorentzVector(negativeP.x(), negativeP.y(), negativeP.z(), antiProtonE), vtx);
        reco::TrackRef neg2TrackRef = theAntiProtonCand.track();
        theAntiProtonCand.setTrack(neg2TrackRef); // negativeTrackRef
     
        AddFourMomenta addp4;
        // Store the daughter Candidates in the VertexCompositeCandidates if they pass mass cuts
        if ( doKShorts_ ) {
          theKshort->addDaughter(thePiPlusCand);
          theKshort->addDaughter(thePiMinusCand);
          theKshort->setPdgId(310);
          // addp4.set(*theKshort);//undefined reference to `AddFourMomenta::set(reco::Candidate&) const even with the include ..
          Candidate::LorentzVector p4(0,0,0,0);
          Candidate::Charge charge = 0;
          p4 += thePiPlusCand.p4();
          p4 += thePiMinusCand.p4();
          theKshort->setP4(p4);
          theKshort->setCharge(charge);
          if ( abs(theKshort->mass() - kShortMass) < kShortMassCut_ ) 
          {
            float K0x = theKshort->vertex().x();
            float K0y = theKshort->vertex().y();
            float K0z = theKshort->vertex().z(); 
            tree_V0_reco_x.push_back(   K0x);
            tree_V0_reco_y.push_back(   K0y);
            tree_V0_reco_z.push_back(   K0z);
            tree_V0_reco_r.push_back(   TMath::Sqrt(K0x*K0x+K0y*K0y));
            tree_V0_reco_drSig.push_back( distsigXY);
            tree_V0_reco_dzSig.push_back( distsigZ);
            tree_V0_reco_angleXY.push_back(  angleXY);
            tree_V0_reco_angleZ.push_back(   detaPointing);
            tree_V0_reco_NChi2.push_back( theKshort->vertexNormalizedChi2());
            tree_V0_reco_ndf.push_back( theKshort->vertexNdof());
            tree_V0_reco_mass.push_back(  theKshort->mass());
            tree_V0_reco_pt.push_back(  theKshort->pt());
            tree_V0_reco_eta.push_back( theKshort->eta());
            tree_V0_reco_phi.push_back( theKshort->phi());
            tree_V0_reco_source.push_back(1);
            tree_V0_reco_badTkHit.push_back(badTkHit);
            tree_V0_reco_dca.push_back( dca);
            temp_nV0_reco++;
            if ( abs(theKshort->mass() - kShortMass) < kShortMassSel_ && ActivateV0Veto) { 
              idxMGT[trd1].first = true;
              idxMGT[trd2].first = true;
            }
          }
        }
        if ( doLambdas_ && theLambda ) {
          theLambda->addDaughter(theProtonCand);
          theLambda->addDaughter(thePiMinusCand);
          theLambda->setPdgId(3122);
          // addp4.set( *theLambda );//undefined reference to `AddFourMomenta::set(reco::Candidate&) const even with the include ..
          Candidate::LorentzVector p4(0,0,0,0);
          Candidate::Charge charge = 0;
          p4 += theProtonCand.p4();
          p4 += thePiMinusCand.p4();
          theLambda->setP4(p4);
          theLambda->setCharge(charge);
          if ( abs(theLambda->mass() - lambdaMass) < lambdaMassCut_ ) 
          {
            float L0x = theLambda->vertex().x();
            float L0y = theLambda->vertex().y();
            float L0z = theLambda->vertex().z(); 
            tree_V0_reco_x.push_back(   L0x);
            tree_V0_reco_y.push_back(   L0y);
            tree_V0_reco_z.push_back(   L0z);
            tree_V0_reco_r.push_back(   TMath::Sqrt(L0x*L0x+L0y*L0y));
            tree_V0_reco_drSig.push_back( distsigXY);
            tree_V0_reco_dzSig.push_back( distsigZ);
            tree_V0_reco_angleXY.push_back(  angleXY);
            tree_V0_reco_angleZ.push_back(   detaPointing);
            tree_V0_reco_NChi2.push_back( theLambda->vertexNormalizedChi2());
            tree_V0_reco_ndf.push_back( theLambda->vertexNdof());
            tree_V0_reco_mass.push_back(  theLambda->mass());
            tree_V0_reco_pt.push_back(  theLambda->pt());
            tree_V0_reco_eta.push_back( theLambda->eta());
            tree_V0_reco_phi.push_back( theLambda->phi());
            tree_V0_reco_source.push_back(2);
            tree_V0_reco_badTkHit.push_back(badTkHit);
            tree_V0_reco_dca.push_back( dca);
            if ( abs(theLambda->mass() - lambdaMass) < lambdaMassSel_ && ActivateV0Veto) { 
              idxMGT[trd1].first = true;
              idxMGT[trd2].first = true;
            }
          }
        } 
	else if ( doLambdas_ && theLambdaBar ) {
          theLambdaBar->addDaughter(theAntiProtonCand);
          theLambdaBar->addDaughter(thePiPlusCand);
          theLambdaBar->setPdgId(-3122);
          // addp4.set(*theLambdaBar);//undefined reference to `AddFourMomenta::set(reco::Candidate&) const even with the include ..
          Candidate::LorentzVector p4(0,0,0,0);
          Candidate::Charge charge = 0;
          p4 += theAntiProtonCand.p4();
          p4 += thePiPlusCand.p4();
          theLambdaBar->setP4(p4);
          theLambdaBar->setCharge(charge);
          if ( abs(theLambdaBar->mass() - lambdaMass) < lambdaMassCut_ ) 
          {
            float L0x = theLambdaBar->vertex().x();
            float L0y = theLambdaBar->vertex().y();
            float L0z = theLambdaBar->vertex().z(); 
            tree_V0_reco_x.push_back(   L0x);
            tree_V0_reco_y.push_back(   L0y);
            tree_V0_reco_z.push_back(   L0z);
            tree_V0_reco_r.push_back(   TMath::Sqrt(L0x*L0x+L0y*L0y));
            tree_V0_reco_drSig.push_back( distsigXY);
            tree_V0_reco_dzSig.push_back( distsigXYZ);
            tree_V0_reco_angleXY.push_back(  angleXY);
            tree_V0_reco_angleZ.push_back(   detaPointing);
            tree_V0_reco_NChi2.push_back( theLambdaBar->vertexNormalizedChi2());
            tree_V0_reco_ndf.push_back( theLambdaBar->vertexNdof());
            tree_V0_reco_mass.push_back(  theLambdaBar->mass());
            tree_V0_reco_pt.push_back(  theLambdaBar->pt());
            tree_V0_reco_eta.push_back( theLambdaBar->eta());
            tree_V0_reco_phi.push_back( theLambdaBar->phi());
            tree_V0_reco_source.push_back(2);
            tree_V0_reco_badTkHit.push_back(badTkHit);
            tree_V0_reco_dca.push_back( dca);
            if ( abs(theLambdaBar->mass() - lambdaMass) < lambdaMassSel_ && ActivateV0Veto) { 
              idxMGT[trd1].first = true;
              idxMGT[trd2].first = true;
            }
          }
        }
        
      }
    }
    tree_nV0_reco = temp_nV0_reco;
    //------------------ END OF V0 reconstruction ---------------------//


    //---------------------------------------------------------------//
    //----------------------- Secondary Interactions ----------------//
    //---------------------------------------------------------------//

    //------------- loop over tracks and vertex good charged track pairs
    for (unsigned int trd1 = 0; trd1 < theTrackRefs.size()-1; ++trd1) 
    {
      if ( idxMGT[trd1].first ) continue;
 
      int iq1 = theTrackRefs[trd1].charge();
      for (unsigned int trd2 = trd1+1; trd2 < theTrackRefs.size(); ++trd2) 
      {
        if ( idxMGT[trd2].first ) continue;
 
        int iq2 = theTrackRefs[trd2].charge();
        float drsig1 = std::abs(theTrackRefs[trd1].dxy(PV.position())/theTrackRefs[trd1].dxyError());
        float drsig2 = std::abs(theTrackRefs[trd2].dxy(PV.position())/theTrackRefs[trd2].dxyError());
        if ( useBS_ ) {
          drsig1 = std::abs(theTrackRefs[trd1].dxy(bs)/theTrackRefs[trd1].dxyError());
          drsig2 = std::abs(theTrackRefs[trd2].dxy(bs)/theTrackRefs[trd2].dxyError());
        }
        if ( !(theTrackRefs[trd1].pt() > pt_Cut && theTrackRefs[trd1].normalizedChi2() < NChi2_Cut && drsig1 > drSig_Cut) 
	  && !(theTrackRefs[trd2].pt() > pt_Cut && theTrackRefs[trd2].normalizedChi2() < NChi2_Cut && drsig2 > drSig_Cut) ) continue;
 
        reco::Track TrackRef1 = theTrackRefs[trd1];
        reco::Track TrackRef2 = theTrackRefs[trd2];
        reco::TransientTrack* TransTkPtr1 = &theTransTracks[trd1];
        reco::TransientTrack* TransTkPtr2 = &theTransTracks[trd2];

        // ---------------measure distance between tracks at their closest approach---------------
        auto const& Impact1 = TransTkPtr1->impactPointTSCP();
        auto const& Impact2 = TransTkPtr2->impactPointTSCP();
        if ( !Impact1.isValid() || !Impact2.isValid() ) continue;

        FreeTrajectoryState const & State1 = Impact1.theState();
        FreeTrajectoryState const & State2 = Impact2.theState();
        ClosestApproachInRPhi cApp;
        cApp.calculate(State1, State2);
        if ( !cApp.status() ) continue;

        float dca = std::abs(cApp.distance());
        if ( dca > tkDCACut_ ) continue;
 
        // the POCA should at least be in the sensitive volume
        GlobalPoint cxPt = cApp.crossingPoint();
        if ( sqrt(cxPt.x()*cxPt.x() + cxPt.y()*cxPt.y()) > 120. || std::abs(cxPt.z()) > 300. ) continue;

        // the tracks should at least point in the same quadrant
        TrajectoryStateClosestToPoint const & TSCP1 = TransTkPtr1->trajectoryStateClosestToPoint(cxPt);
        TrajectoryStateClosestToPoint const & TSCP2 = TransTkPtr2->trajectoryStateClosestToPoint(cxPt);
        if ( !TSCP1.isValid() || !TSCP2.isValid() ) continue;
        if ( TSCP1.momentum().dot(TSCP2.momentum()) < 0. ) continue;
        	  
        // calculate the invariant mass of the pair
        double totalE = sqrt(TSCP1.momentum().mag2()) + sqrt(TSCP2.momentum().mag2());
        double totalESq = totalE*totalE;
        double PtotSq = (TSCP1.momentum() + TSCP2.momentum()).mag2();
        double mass = TMath::Sqrt(totalESq - PtotSq);
        if ( mass > 4. ) continue; 
 
        // Fill the vector of TransientTracks to send to KVF
        std::vector<reco::TransientTrack> transTracks;
        transTracks.reserve(2);
        transTracks.push_back(*TransTkPtr1);
        transTracks.push_back(*TransTkPtr2);
        
        // create the vertex fitter object and vertex the tracks
        TransientVertex theRecoVertex;
        if ( vertexFitter_ ) // true is recommended, AVF is more likely to be better for high number of tracks vertices
        {
          KalmanVertexFitter theKalmanFitter(useRefTracks_ == 0 ? false : true); 
          theRecoVertex = theKalmanFitter.vertex(transTracks);
        } 
        else if ( !vertexFitter_ ) 
        {
          useRefTracks_ = false;
          AdaptiveVertexFitter theAdaptiveFitter;
          theRecoVertex = theAdaptiveFitter.vertex(transTracks);
        }
        if ( !theRecoVertex.isValid() ) continue;

        reco::Vertex theVtx = theRecoVertex;
        if ( theVtx.normalizedChi2() > 20. ) continue;

        GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());
        // 2D decay length significance 
	      float distsigXY = 10000.;
        SMatrixSym3D  totalCov = PV.covariance() + theVtx.covariance();
        if ( useBS_ ) totalCov = bs.rotatedCovariance3D() + theVtx.covariance();
        SVector3 distVecXY(vtxPos.x()-PV.x(), vtxPos.y()-PV.y(), 0.);
        double distMagXY = ROOT::Math::Mag(distVecXY);
        double sigmaDistMagXY = sqrt(ROOT::Math::Similarity(totalCov, distVecXY)) / distMagXY;
        if ( sigmaDistMagXY > 0. ) distsigXY = distMagXY / sigmaDistMagXY;
        if ( distsigXY < 50. ) continue;
   	     
        // z significance
        float distsigZ = 10000.;
        SVector3 distVecZ(0, 0, vtxPos.z()-PV.z());
        double distMagZ = TMath::Abs(vtxPos.z()-PV.z());
        double sigmaDistMagZ = sqrt(ROOT::Math::Similarity(theVtx.covariance(), distVecZ)) / distMagZ;
        if ( sigmaDistMagZ > 0. ) distsigZ = distMagZ / sigmaDistMagZ;
       	     
        // make sure the vertex radius is within the inner track hit radius
        bool badTkHit = false;
        for (int k = 1; k <= 2; k++) 
        {
          unsigned int  trdx = trd1;
          if ( k == 2 ) trdx = trd2;
          if ( idxMGT[trdx].second >= pc->size() ) continue; // not for the lost tracks, but as their default first hit is PIXBL1 it would not hurt...

          reco::Track   TrackRef = TrackRef1;
          if ( k == 2 ) TrackRef = TrackRef2;
          const HitPattern hp = TrackRef.hitPattern();
          uint16_t firsthit = hp.getHitPattern(HitPattern::HitCategory::TRACK_HITS,0);     
           //---Creating State to propagate from  TT---//
          reco::TransientTrack TTrack = theTransientTrackBuilder->build(TrackRef);
          GlobalPoint vert (TrackRef.vx(),TrackRef.vy(),TrackRef.vz()); // Point where the propagation will start (Reference Point)
          const TrajectoryStateOnSurface Surtraj = TTrack.stateOnSurface(vert); // TSOS of this point
          const MagneticField* B = TTrack.field(); // 3.8T
          AnalyticalPropagator* Prop = new AnalyticalPropagator(B);
          Basic3DVector<float> P3D2(TrackRef.vx(),TrackRef.vy(),TrackRef.vz()); // global frame
          Basic3DVector<float> B3DV(TrackRef.px(),TrackRef.py(),TrackRef.pz()); // global frame 
          float Eta = TrackRef.eta();
          float Phi = TrackRef.phi();
          float vz  = TrackRef.vz();
          
          std::pair<int,GloballyPositioned<float>::PositionType> FHPosition = PHP->Main(firsthit,Prop,Surtraj,Eta,Phi,vz,P3D2,B3DV);

          float xF = FHPosition.second.x() - PV.x();
          float yF = FHPosition.second.y() - PV.y();
          float rF = TMath::Sqrt( xF*xF + yF*yF );
          float zF = TMath::Abs( FHPosition.second.z() - PV.z() );
	  // PXB
	  if ( rF > 2.7 && rF < 16.5 && zF < 27.0 ) {
	    if ( rF <  distMagXY - 0.40 ) badTkHit = true;
          }
        	  // TIB
	  if ( rF > 23.5 && rF < 36.2 && zF < 66.0 ) {
	    if ( rF <  distMagXY - 0.84 ) badTkHit = true;
	  }
	  if ( rF > 40.0 && rF < 52.0 && zF < 66.0 ) {
	    if ( rF <  distMagXY - 2.24 ) badTkHit = true;
	  }
	  // TOB
	  if ( rF > 58.4 && rF < 62.9 && zF < 107. ) {
	    if ( rF <  distMagXY - 2.8 ) badTkHit = true;
	  }
	  // PXF
	  if ( rF > 4.5 && rF < 16.1 && zF > 29.6 && zF < 52.0 ) {
	    if ( zF <  distMagZ - 3.2 ) badTkHit = true;
	  }
	  // TID
	  if ( rF > 22.5 && rF < 50.5 && zF > 74.3 && zF < 109.7 ) {
	    if ( zF <  distMagZ - 3.9 ) badTkHit = true;
	  }
	  // TEC
	  if ( rF > 22.0 && rF < 70.0 && zF > 126.4 && zF < 193.7 ) {
	    if ( zF <  distMagZ - 6.7 ) badTkHit = true;
	  }
        }
        if ( badTkHit && dca > 0.1 ) continue;
 
        std::unique_ptr<TrajectoryStateClosestToPoint> traj1;
        std::unique_ptr<TrajectoryStateClosestToPoint> traj2;
        std::vector<reco::TransientTrack> theRefTracks;
        if ( theRecoVertex.hasRefittedTracks() ) theRefTracks = theRecoVertex.refittedTracks();
                  
        if ( useRefTracks_ && theRefTracks.size() > 1 ) 
        {
          reco::TransientTrack* theRefTrack1 = nullptr;
          reco::TransientTrack* theRefTrack2 = nullptr;
	  int idum = 0;
          for (std::vector<reco::TransientTrack>::iterator iTrack = theRefTracks.begin(); iTrack != theRefTracks.end(); ++iTrack) 
          {
	    idum++;
            if ( idum == 1 ) theRefTrack1 = &*iTrack;
            else	     theRefTrack2 = &*iTrack;
          }
          if ( theRefTrack1 == nullptr || theRefTrack2 == nullptr ) continue;
          traj1.reset(new TrajectoryStateClosestToPoint(theRefTrack1->trajectoryStateClosestToPoint(vtxPos)));
          traj2.reset(new TrajectoryStateClosestToPoint(theRefTrack2->trajectoryStateClosestToPoint(vtxPos)));
        } 
        else 
        {
          traj1.reset(new TrajectoryStateClosestToPoint(TransTkPtr1->trajectoryStateClosestToPoint(vtxPos)));
          traj2.reset(new TrajectoryStateClosestToPoint(TransTkPtr2->trajectoryStateClosestToPoint(vtxPos)));
        }
        if ( traj1.get() == nullptr || traj2.get() == nullptr || !traj1->isValid() || !traj2->isValid() ) continue;
                   
        // calculate total energy assuming 0 mass for the tracks
        GlobalVector P1(traj1->momentum());
        GlobalVector P2(traj2->momentum());
        GlobalVector Ptot(P1 + P2);	
        double E1 = sqrt(P1.mag2());
        double E2 = sqrt(P2.mag2());
        double ETot = E1 + E2;
        const reco::Particle::LorentzVector SecIntP4(Ptot.x(), Ptot.y(), Ptot.z(), ETot);
     
        // 2D pointing angle
        double dx = theVtx.x()-PV.x();
        double dy = theVtx.y()-PV.y();
        double px = Ptot.x();
        double py = Ptot.y();
        double pz = Ptot.z();
        double ptot = TMath::Sqrt(px*px + py*py + pz*pz);
        float  angleXY = (dx*px+dy*py)/(sqrt(dx*dx+dy*dy)*sqrt(px*px+py*py));
        if ( angleXY < 0.7 ) continue;
                     
        // delta eta pointing
        double dr = sqrt(dx*dx + dy*dy);
        double etaZ = 0.;
        if ( dz != 0. ) etaZ = -TMath::Log(abs(TMath::Tan(TMath::ATan(dr/dz)/2.)));
        if ( dz < 0. )  etaZ = -etaZ;
        double etaP = 0.;
        if ( abs(pz) < ptot ) etaP = 0.5 * TMath::Log((ptot+pz)/(ptot-pz));
        float detaPointing = etaP - etaZ;

        reco::Particle::Point vtx(theVtx.x(), theVtx.y(), theVtx.z());
        const reco::Vertex::CovarianceMatrix vtxCov(theVtx.covariance());
        double vtxChi2(theVtx.chi2());
        double vtxNdof(theVtx.ndof());
      
        // Create the VertexCompositeCandidate object that will be stored in the Event
        reco::VertexCompositeCandidate* theSecInt = new reco::VertexCompositeCandidate(0, SecIntP4, vtx, vtxCov, vtxChi2, vtxNdof);
     
        // Create daughter candidates for the VertexCompositeCandidates
        reco::RecoChargedCandidate theTkCand1( iq1, reco::Particle::LorentzVector(P1.x(), P1.y(), P1.z(), E1), vtx);
        reco::TrackRef TkRef1 = theTkCand1.track();
        theTkCand1.setTrack(TkRef1);
        
        reco::RecoChargedCandidate theTkCand2( iq2, reco::Particle::LorentzVector(P2.x(), P2.y(), P2.z(), E2), vtx);
        reco::TrackRef TkRef2 = theTkCand2.track();
        theTkCand2.setTrack(TkRef2);
        
        AddFourMomenta addp4;
        // Store the daughter Candidates in the VertexCompositeCandidates if they pass mass cuts
        theSecInt->addDaughter(theTkCand1);
        theSecInt->addDaughter(theTkCand2);
        theSecInt->setPdgId(22);
        Candidate::LorentzVector p4(0,0,0,0);
        Candidate::Charge charge = iq1+iq2;
        p4 += theTkCand1.p4();
        p4 += theTkCand2.p4();
        theSecInt->setP4(p4);
        theSecInt->setCharge(charge);
 
        if ( abs(theSecInt->mass()) > 4. ) continue; 
 
        float SecInt_x = theSecInt->vertex().x();
        float SecInt_y = theSecInt->vertex().y();
        float SecInt_z = theSecInt->vertex().z();
        float SecInt_r = TMath::Sqrt(SecInt_x*SecInt_x + SecInt_y*SecInt_y);
        float SecInt_d = TMath::Sqrt(SecInt_x*SecInt_x + SecInt_y*SecInt_y+SecInt_z*SecInt_z);
        tree_SecInt_x.push_back(	 SecInt_x);
        tree_SecInt_y.push_back(	 SecInt_y);
        tree_SecInt_z.push_back(	 SecInt_z);
        tree_SecInt_r.push_back(	 SecInt_r);
        tree_SecInt_d.push_back(	 SecInt_d);
        tree_SecInt_drSig.push_back(   distsigXY);
        tree_SecInt_dzSig.push_back(   distsigZ);
        tree_SecInt_angleXY.push_back( angleXY);
        tree_SecInt_angleZ.push_back(  detaPointing);
        tree_SecInt_NChi2.push_back(   theSecInt->vertexNormalizedChi2());
        tree_SecInt_ndf.push_back(     theSecInt->vertexNdof());
        tree_SecInt_mass.push_back(    theSecInt->mass());
        tree_SecInt_pt.push_back(      theSecInt->pt());
        tree_SecInt_eta.push_back(     theSecInt->eta());
        tree_SecInt_phi.push_back(     theSecInt->phi());
        tree_SecInt_charge.push_back(  theSecInt->charge());
        tree_SecInt_badTkHit.push_back(badTkHit);
        tree_SecInt_dca.push_back(     dca);
        tree_nSecInt++;

        // get r-z distance to closest generated LLP decay point
        if ( isMC_ )
        {
	  float drLLP, dzLLP, ddLLP1, ddLLP2;
	  ddLLP1 = (SecInt_x - LLP1_x)*(SecInt_x - LLP1_x) + (SecInt_y - LLP1_y)*(SecInt_y - LLP1_y) + (SecInt_z - LLP1_z)*(SecInt_z - LLP1_z);
	  ddLLP2 = (SecInt_x - LLP2_x)*(SecInt_x - LLP2_x) + (SecInt_y - LLP2_y)*(SecInt_y - LLP2_y) + (SecInt_z - LLP2_z)*(SecInt_z - LLP2_z);
          float ddLLP = -100;
	  if ( ddLLP1 < ddLLP2 ) {
            float deltaX = SecInt_x - LLP1_x;
            float deltaY = SecInt_y - LLP1_y;
	    drLLP = TMath::Sqrt( (SecInt_x - LLP1_x)*(SecInt_x - LLP1_x) + (SecInt_y - LLP1_y)*(SecInt_y - LLP1_y) );
            if ((deltaX<0 && deltaY<0) || (deltaX<0 && abs(deltaX)>abs(deltaY)) || (deltaY<0 && abs(deltaY)>abs(deltaX))){drLLP = -drLLP;}
	    dzLLP = SecInt_z - LLP1_z;
            ddLLP = ddLLP1;
	  }
	  else {
            float deltaX = SecInt_x - LLP2_x;
            float deltaY = SecInt_y - LLP2_y;
	    drLLP = TMath::Sqrt( (SecInt_x - LLP2_x)*(SecInt_x - LLP2_x) + (SecInt_y - LLP2_y)*(SecInt_y - LLP2_y) );
            if ((deltaX<0 && deltaY<0) || (deltaX<0 && abs(deltaX)>abs(deltaY)) || (deltaY<0 && abs(deltaY)>abs(deltaX))){drLLP = -drLLP;}
	    dzLLP = SecInt_z - LLP2_z;
            ddLLP = ddLLP2;
	  }
          tree_SecInt_LLP_dr.push_back(  drLLP);
          tree_SecInt_LLP_dz.push_back(  dzLLP);
          tree_SecInt_LLP_dd.push_back(  ddLLP);
	
	  // are the 2 tracks from LLP ?
          float   q1 = theTrackRefs[trd1].charge();
          float  pt1 = theTrackRefs[trd1].pt();
          float eta1 = theTrackRefs[trd1].eta();
          float phi1 = theTrackRefs[trd1].phi();
          int   hit1 = theTrackRefs[trd1].hitPattern().numberOfValidHits();
          float   q2 = theTrackRefs[trd2].charge();
          float  pt2 = theTrackRefs[trd2].pt();
          float eta2 = theTrackRefs[trd2].eta();
          float phi2 = theTrackRefs[trd2].phi();
          int   hit2 = theTrackRefs[trd2].hitPattern().numberOfValidHits();
          int LLPtrd1 = 0, LLPtrd2 = 0;

          for (int k = 0; k < tree_ngenFromLLP; k++) // loop on final gen part from LLP
          {
            float qGen   = tree_genFromLLP_charge[k];
            if ( q1 != qGen && q2 != qGen ) continue;

            int kLLP =     tree_genFromLLP_LLP[k];
            float ptGen  = tree_genFromLLP_pt[k];
            float etaGen = tree_genFromLLP_eta[k];
            float phiGen = tree_genFromLLP_phi[k]; // given at production point
            float xGen   = tree_genFromLLP_x[k];
            float yGen   = tree_genFromLLP_y[k];
            float qR = qGen * ptGen * 100 / 0.3 / 3.8;
	    // compute phi0 at dca(PV) for the gen particle (instead of production point)
            float sin0 = qR * sin( phiGen ) + (xGen - tree_GenPVx);
            float cos0 = qR * cos( phiGen ) - (yGen - tree_GenPVy);
            float phi0 = TMath::ATan2( sin0, cos0 ); // but note that it can be wrong by +_pi ! 
            float dpt= 100., deta = 100., dphi = 100.;
            if ( q1 == qGen && LLPtrd1 == 0 ) {
              dpt  = (pt1 - ptGen) / pt1;
              deta = eta1 - etaGen;
              dphi = phi1 - phi0;
              if	  ( dphi < -3.14159 / 2. ) dphi += 3.14159;
              else if ( dphi >  3.14159 / 2. ) dphi -= 3.14159;
              if ( hit1 <= 10 ) {
                if ( abs(dpt) < 0.70 && abs(deta) < 0.30 && abs(dphi) < 0.08 ) LLPtrd1 = kLLP; 
              }
              else if ( hit1 <= 13 ) {
                if ( abs(dpt) < 0.20 && abs(deta) < 0.12 && abs(dphi) < 0.05 ) LLPtrd1 = kLLP; 
              }
              else if ( hit1 <= 17 ) {
                if ( abs(dpt) < 0.08 && abs(deta) < 0.04 && abs(dphi) < 0.03 ) LLPtrd1 = kLLP; 
              }
              else {
                if ( abs(dpt) < 0.07 && abs(deta) < 0.02 && abs(dphi) < 0.02 ) LLPtrd1 = kLLP; 
              }
	    }
	    if ( q2 == qGen && LLPtrd2 == 0 ) {
              dpt  = (pt2 - ptGen) / pt2;
              deta = eta2 - etaGen;
              dphi = phi2 - phi0;
              if ( dphi < -3.14159 / 2. ) dphi += 3.14159;
              else if ( dphi >  3.14159 / 2. ) dphi -= 3.14159;
	      if ( hit2 <= 10 ) {
                if ( abs(dpt) < 0.70 && abs(deta) < 0.30 && abs(dphi) < 0.08 ) LLPtrd2 = kLLP; 
              }
              else if ( hit2 <= 13 ) {
                if ( abs(dpt) < 0.20 && abs(deta) < 0.12 && abs(dphi) < 0.05 ) LLPtrd2 = kLLP; 
              }
              else if ( hit2 <= 17 ) {
                if ( abs(dpt) < 0.08 && abs(deta) < 0.04 && abs(dphi) < 0.03 ) LLPtrd2 = kLLP; 
              }
              else {
                if ( abs(dpt) < 0.07 && abs(deta) < 0.02 && abs(dphi) < 0.02 ) LLPtrd2 = kLLP; 
              }
	    }
          } // end loop on final gen part from LLP
	  if ( LLPtrd1 != 0 && LLPtrd1 == LLPtrd2 ) {
            tree_SecInt_LLP.push_back(LLPtrd1);
	  }
	  else {
            tree_SecInt_LLP.push_back(0);
	  }
        } // endif MC

	bool SecInt_selec = false;
	if ( distsigXY > 100. && angleXY > 0.9 && theSecInt->mass() < 2. &&
             theSecInt->vertexNormalizedChi2() < 10. ) SecInt_selec = true; //&& dca < 1.
        tree_SecInt_selec.push_back(   SecInt_selec);

        // redefine SecInt_r to account for the displaced detector center in data
        if ( !isMC_ ) SecInt_r = TMath::Sqrt( (SecInt_x-x_det)*(SecInt_x-x_det) + (SecInt_y-y_det)*(SecInt_y-y_det) );

	// tracker active layers
        // PropaHitPattern* NI = new PropaHitPattern();
        // Shift from the (0,0) center of CMS w.r.t the alignement of the tracker.
        // / With this correction => everything is in the (0,0) reference frame

        int VtxLayerNI = -1;
        if (DetailedMap) VtxLayerNI = NI->VertexBelongsToTracker(SecInt_r, SecInt_z);
        else { // deprecated
          VtxLayerNI = NI->VertexBelongsToBarrelLayer(SecInt_r, SecInt_z);
          if ( VtxLayerNI == 0 ) VtxLayerNI = NI->VertexBelongsToDiskLayer(SecInt_r, SecInt_z);
        }
 
	if ( SecInt_selec && ActivateSecIntVeto ) {
	  if ( VtxLayerNI != 0 ) {
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // beam pipe
          float r_bmp = TMath::Sqrt( (SecInt_x-x_bmp)*(SecInt_x-x_bmp) + (SecInt_y-y_bmp)*(SecInt_y-y_bmp) );
	  if ( abs(SecInt_z) < 27. && r_bmp > 2.15 && r_bmp < 2.27 ) {
      	    VtxLayerNI = -1;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // PIXB inner support
          if ( abs(SecInt_z) < 27. && SecInt_r > 2.44 && SecInt_r < 2.55 ) {
	    VtxLayerNI = -2;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // PIXB outer support
          if ( SecInt_r > 21.4 && SecInt_r < 22.1 ) {
	    VtxLayerNI = -3;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // PIXB rails
	  if ( abs(SecInt_z) < 27. && abs(SecInt_x) < 10. && 
	       abs(SecInt_y) > 18.9 && abs(SecInt_y) < 20.9 ) {
	    VtxLayerNI = -4;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // services : r 18 - 19 and z 29 - 200 according to "other" material map
	  if ( abs(SecInt_z) > 29. && SecInt_r > 18.0 && SecInt_r < 20.6 ) {
	    VtxLayerNI = -3;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	}
        tree_SecInt_layer.push_back( VtxLayerNI );
        tree_SecInt_tk1.push_back(   trd1 );
        tree_SecInt_tk2.push_back(   trd2 );
      }
    }
    //---------------------------END OF Sec.Int. reconstruction------------------------------------------//
    

    //---------------------------------------------------------------//
    //------------------------- TRACKS ------------------------------//
    //---------------------------------------------------------------//

    for (unsigned int ipc = 0; ipc < TRACK_SIZE; ipc++) { // loop on all packedPFCandidates + lostTrackss
      pat::PackedCandidateRef pcref = MINIgeneralTracks[ipc];
      float energy = pcref->energy();
      const reco::Track *trackPcPtr = pcref->bestTrack(); // const
      if ( !trackPcPtr ) continue;

      // reject tracks from K0 and Lambda decays and from secondary interactions
      bool Veto_tk = false;
      for (unsigned int j = 0 ; j < idxMGT.size() ; j++) 
      {
        if ( idxSecIntMGT[j].second != idxMGT[j].second ) cout << " ERROR idxSecIntMGT " << endl;
        if ( (idxMGT[j].first || idxSecIntMGT[j].first) && idxMGT[j].second == ipc ) 
        {
          Veto_tk = true;//keep as true
          break;
        }
      }
      if ( Veto_tk ) continue;
 
      reco::Track tk = *trackPcPtr;
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
      if ( tk_dzError > 0 )  tk_dzSig = abs(tk_dz) / tk_dzError; 
      if ( tk_dxyError > 0 ) tk_drSig = abs(tk_dxy) / tk_dxyError; 
      tk_vx = tk.vx();
      tk_vy = tk.vy();
      tk_vz = tk.vz();
      tk_px = tk.px();
      tk_py = tk.py();
      tk_pz = tk.pz();
      tk_e = energy;
      tk_HitPattern = tk.hitPattern();

      if ( tk_nHit == 0 ) continue;
      if ( tk_charge == 0 ) continue;

      // preselection
      if ( RequestHighPurity
             && !( tk_pt > pt_Cut && tk_NChi2 < NChi2_Cut && tk_drSig > drSig_Cut
             && static_cast<int>(tk.quality(reco::TrackBase::highPurity))) ) continue;
      else if ( !( tk_pt > pt_Cut && tk_NChi2 < NChi2_Cut && tk_drSig > drSig_Cut ) ) continue;
 
      // reject prompt muon tracks
      bool PromptMuonVeto = false;
      float dptMu1 = 10., detaMu1 = 10., dphiMu1 = 10.;
      float dptMu2 = 10., detaMu2 = 10., dphiMu2 = 10.;
      int iqMu1 = 0, iqMu2 = 0;
      if ( imu1 >= 0 ) {
        dptMu1 = (tk_pt - lep1_pt) / tk_pt;
        detaMu1 = tk_eta - lep1_eta;
        dphiMu1 = Deltaphi( tk_phi, lep1_phi );
        iqMu1 = lep1_Q;
      }
      if ( imu2 >= 0 ) {
        dptMu2 = (tk_pt - lep2_pt) / tk_pt;
        detaMu2 = tk_eta - lep2_eta;
        dphiMu2 = Deltaphi( tk_phi, lep2_phi );
        iqMu2 =lep2_Q;
      }
      if ( abs(dptMu1) < 0.1 && abs(detaMu1) < 0.1 && abs(dphiMu1) < 0.1 
	   && tk_charge == iqMu1 ) PromptMuonVeto = true;		 
      if ( abs(dptMu2) < 0.1 && abs(detaMu2) < 0.1 && abs(dphiMu2) < 0.1
	   && tk_charge == iqMu2 ) PromptMuonVeto = true;		 
      if ( PromptMuonVeto ) continue;
 
      // reject tracks from conversions
      bool Yc_tk = false;
      for (int k = 0; k < tree_Yc_ntracks; k++) {   // Loop on conversions
        float Yq    = tree_Yc_tracks_charge[k];
        if ( tk_charge*Yq < 0. ) continue;
        float Ypt   = tree_Yc_tracks_pt[k];
        float Yeta  = tree_Yc_tracks_eta[k];
        float Yphi0 = tree_Yc_tracks_phi0[k];
        float dpt = ( tk_pt - Ypt ) / tk_pt;
        float deta = tk_eta - Yeta;
        float dphi0 = tk_phi - Yphi0;
        if      ( dphi0 < -3.14159 / 2. ) dphi0 += 3.14159;
        else if ( dphi0 >  3.14159 / 2. ) dphi0 -= 3.14159;
        if ( abs(dpt) < 0.2 && abs(deta) < 0.1 && abs(dphi0) < 0.1 ) {
          Yc_tk = true;
          break;
        }
      }
      if ( Yc_tk && ActivateYcVeto ) continue;

      tree_nTracks++;
      tree_track_ipc.push_back(ipc);
      if ( ipc < pc->size() ) {
        tree_track_lost.push_back(0);
      }
      else {
        tree_track_lost.push_back(1);
        tree_nLostTracks++;
      }
      tree_track_px.push_back	  (tk_px);
      tree_track_py.push_back	  (tk_py);
      tree_track_pz.push_back	  (tk_pz);
      tree_track_pt.push_back	  (tk_pt);
      tree_track_eta.push_back	  (tk_eta);
      tree_track_phi.push_back	  (tk_phi);
      tree_track_charge.push_back   (tk_charge);
      tree_track_NChi2.push_back    (tk_NChi2);
      tree_track_dxy.push_back	  (tk_dxy);
      tree_track_dxyError.push_back (tk_dxyError);
      tree_track_drSig.push_back    (tk_drSig); 
      tree_track_dz.push_back	  (tk_dz);
      tree_track_dzError.push_back  (tk_dzError);
      tree_track_dzSig.push_back    (tk_dzSig);
      tree_track_energy.push_back   (tk_e);

//$$$$$$
//       // Tracks from pileup
//       float tk_dzmin = 100.;
//       float tk_dzminSig = 100000.;
//       for (unsigned int i = 0; i< primaryVertex->size() ; ++i) {
//         if ( i == 0 ) continue;
//         float dzTOpu = tk_dz + tree_PV_z - (*primaryVertex)[i].z();
//         if ( abs(dzTOpu) < abs(tk_dzmin) ) tk_dzmin = dzTOpu;
//         float dzerror = TMath::Sqrt( tk_dzError*tk_dzError + (*primaryVertex)[i].zError()*(*primaryVertex)[i].zError() );
//         if ( dzerror > 0. ) {
//           if ( abs(dzTOpu / dzerror) < abs(tk_dzminSig) ) tk_dzminSig = dzTOpu / dzerror;
// 	   }
//       }
//       tree_track_dzTOpu.push_back       (tk_dzmin);
//       tree_track_dzSigTOpu.push_back    (tk_dzminSig);
//$$$$$$

      tree_track_isHighPurity.push_back (static_cast<int>(tk.quality(reco::TrackBase::highPurity)));
      tree_track_nHit.push_back	        (tk_nHit);
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
      tree_track_x.push_back	  (tk_vx);
      tree_track_y.push_back	  (tk_vy);
      tree_track_z.push_back	  (tk_vz);

                  //----------------MINIAOD_Firsthit----------//
                  //-----------------IMPORTANT----------------//
                  // TSOS is said to be better for the -------//
                  // propagators (see Propagator.h)...--------//
                  // ../interface/PropaHitPattern.h	    //
                  //------------------------------------------//
      //-----hitpattern -> Database ---/
      const HitPattern hp = tk_HitPattern;
      uint16_t firsthit = hp.getHitPattern(HitPattern::HitCategory::TRACK_HITS,0);
 
      // Approximation for the lostTrack since the hitpattern information is not available (only 1160, tracking POG knows about it but do not seem to care)      
      if ( ipc >= pc->size() ) {
        if ( abs(tk_eta) < 1. && YEAR_ >= 2017) firsthit = 1184; // PIXBL4 in barrel
        else if (abs(tk_eta) < 1. && YEAR_ == 2016) firsthit = 1176;
        else		      firsthit = 1296; // PIXFD2 in forward
      }      
 
      tree_track_firstHit.push_back(firsthit);

      //---Creating State to propagate from  TT---//
      BestTracks.push_back(theTransientTrackBuilder->build(tk));
      const MagneticField* B = BestTracks[count].field(); // 3.8T
      reco::TransientTrack TT (tk,BestTracks[count].field());
      // const FreeTrajectoryState Freetraj = TT.initialFreeState(); // Propagator in the barrel can also use FTS (WARNING: the so-called reference point (where the propagation starts might be different from the first vtx, a check should be done))
      GlobalPoint vert (tk_vx,tk_vy,tk_vz); // Point where the propagation will start (Reference Point)
      const TrajectoryStateOnSurface Surtraj = TT.stateOnSurface(vert); // TSOS of this point
      AnalyticalPropagator* Prop = new AnalyticalPropagator(B);
      Basic3DVector<float> P3D2(tk_vx,tk_vy,tk_vz);  // global frame
      Basic3DVector<float> B3DV (tk_px,tk_py,tk_pz); // global frame 
      float Eta = tk_eta;
      float Phi = tk_phi;
      float vz  = tk_vz;
      //------Propagation with new interface --> See ../interface/PropaHitPattern.h-----//
      
      std::pair<int,GloballyPositioned<float>::PositionType> FHPosition = PHP->Main(firsthit,Prop,Surtraj,Eta,Phi,vz,P3D2,B3DV);
      // returns 0 if in Barrel, 1 if disks with the position of the firsthit. Different propagators are used between 
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
      float btagFromJet = 0;
      bool matchTOjet = false;

      for (const pat::Jet &jet : *jets) {
      if ( jet.pt() < jet_pt_min ) continue;
        float dR_jet = Deltar( jet.eta(), jet.phi(), tk_eta, tk_phi );
        float DeepFlavourb = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
        float DeepFlavourbb = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
        float DeepFlavourblep = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
      if ( DeepFlavourb > -5 ) btagFromJet = DeepFlavourb + DeepFlavourbb + DeepFlavourblep;
        if ( dR_jet < 0.4 ) {
          matchTOjet = true;
          break;
        }
        else iJet++;
      }
      if ( matchTOjet ) {tree_track_iJet.push_back (iJet); tree_track_btag.push_back(btagFromJet);}
      else	        {tree_track_iJet.push_back (-1); tree_track_btag.push_back(-1);}

      // match to gen particle from LLP decay
      if ( isMC_ )
      {
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
        float sin0 = qR * sin( phiGen ) + (xGen - tree_GenPVx);
        float cos0 = qR * cos( phiGen ) - (yGen - tree_GenPVy);
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
            if ( abs(dpt) < 1.00 && abs(deta) < 0.30 && abs(dphi) < 0.09 ) matchTOgen = true; 
          }
          else if ( tk_nHit <= 13 ) {
            if ( abs(dpt) < 0.30 && abs(deta) < 0.13 && abs(dphi) < 0.05 ) matchTOgen = true; 
          }
          else if ( tk_nHit <= 17 ) {
            if ( abs(dpt) < 0.12 && abs(deta) < 0.04 && abs(dphi) < 0.03 ) matchTOgen = true; 
          }
          else {
            if ( abs(dpt) < 0.08 && abs(deta) < 0.02 && abs(dphi) < 0.02 ) matchTOgen = true; 
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
    } // endif MC

    } // end loop on all track candidates


    ///////////////////////////////////////////////////////
    //-----------------------------------------------------
    // selection of displaced tracks
    //-----------------------------------------------------
    ///////////////////////////////////////////////////////

    int jet; /*!*/

    int nTrks_axis1 = 0;
    int nTrks_axis1_sig=0, nTrks_axis1_bad=0;
    int nTrks_axis2 = 0;
    int nTrks_axis2_sig=0, nTrks_axis2_bad=0;
    int nTrks_axis1_sig_mva=0, nTrks_axis1_bad_mva=0;
    int nTrks_axis2_sig_mva=0, nTrks_axis2_bad_mva=0;
    
    LLP1_nTrks = 0;
    LLP2_nTrks = 0;

    // ------------------------------------- //
    // --------- Track Selection BDT WPs --- //
    // ------------------------------------- //

//$$$$$$
    double bdtcut = 0.92; //Tight WP      // ttbar ~ 1E-3 : selection efficiency
//$$$$$$
    double bdtcut_step2 = 0.0; //Loose WP // ttbar ~ 1E-2
 
    //---------------------------//

    for (int counter_track = 0; counter_track < tree_nTracks; counter_track++) 
    {
      firsthit_X = tree_track_firstHit_x[counter_track];
      firsthit_Y = tree_track_firstHit_y[counter_track];
      firsthit_Z = tree_track_firstHit_z[counter_track];
      pt         = tree_track_pt[counter_track];
      eta        = tree_track_eta[counter_track];
      phi        = tree_track_phi[counter_track];
      NChi       = tree_track_NChi2[counter_track];
      nhits      = tree_track_nHit[counter_track];
      dxy        = abs(tree_track_dxy[counter_track]);
      dz         = abs(tree_track_dz[counter_track]);
      drSig      = tree_track_drSig[counter_track];
      dzSig      = tree_track_dzSig[counter_track];

//$$$$$$
//       dzTopu     = tree_track_dzTOpu[counter_track];
//       dzSigTopu  = tree_track_dzSigTOpu[counter_track];
//$$$$$$
      TibHit     = tree_track_nHitTIB[counter_track] ;
      TobHit     = tree_track_nHitTOB[counter_track] ;
      PixBarHit  = tree_track_nHitPXB[counter_track];
      TecHit     = tree_track_nHitTEC[counter_track];

      ntrk10 = 0, ntrk20 = 0, ntrk30 = 0, ntrk40 = 0;
      isLost     = tree_track_lost[counter_track];
      float ntrk10_lost = 0, ntrk20_lost = 0, ntrk30_lost = 0, ntrk40_lost = 0;
      isinjet = 0.;
      double bdtval = -10.;
      dR = -1.;
      int tracks_axis = 0; // flag to check which axis is the closest from the track

      jet = tree_track_iJet[counter_track];
      int isFromLLP = -1;
      if ( jet >= 0 ) isinjet = 1.; /*!*/

      if ( isMC_ ) isFromLLP = tree_track_sim_LLP[counter_track];
              
      // check the dR between the tracks and the second axis (without any selection on the tracks)
      dR1  = Deltar( eta, phi, axis1_eta, axis1_phi ); // axis1_phi and axis1_eta for the first axis
      dR2  = Deltar( eta, phi, axis2_eta, axis2_phi );
      if ( dR1 < dR2 &&  dR1 < dRcut_tracks  ) { // a restriction could be added on the value of dR to assign the value Tracks_axis  (avoid some background???)
        tracks_axis = 1;
        dR = dR1;
        dRmax = dR2;
      }
      if ( dR2 < dR1 &&  dR2 < dRcut_tracks  ) { // a restriction could be added on the value of dR to assign the value Tracks_axis  (avoid some background???)
        tracks_axis = 2;
        dR = dR2;
        dRmax = dR1;
      }

      //Computation of the distances needed for the BDT
      for (int counter_othertrack = 0; counter_othertrack < tree_nTracks; counter_othertrack++) 
      {
        if ( counter_othertrack == counter_track ) continue;

        float x2 = tree_track_firstHit_x[counter_othertrack];
        float y2 = tree_track_firstHit_y[counter_othertrack];
        float z2 = tree_track_firstHit_z[counter_othertrack];
        float dist = TMath::Sqrt( (firsthit_X-x2)*(firsthit_X-x2) + (firsthit_Y-y2)*(firsthit_Y-y2) + (firsthit_Z-z2)*(firsthit_Z-z2) ); // pour chaque reconstruite, on regarde les autres tracks
        if ( dist < 10.  ) ntrk10++; //&& tree_track_lost[counter_othertrack] == false
        if ( dist < 20. ) ntrk20++; //&& tree_track_lost[counter_othertrack] == false
        if ( dist < 30. ) ntrk30++; //&& tree_track_lost[counter_othertrack] == false
        if ( dist < 40.) ntrk40++; // && tree_track_lost[counter_othertrack] == false
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
      
//$$$$$$
//       if ( ntrk40 == 0 ) {
//         ntrk10rel = -0.02;
//         ntrk20rel = -0.02;
//         ntrk30rel = -0.02; 
//       }
//       else {
//         ntrk10rel = ntrk10 / ntrk40;
//         ntrk20rel = ntrk20 / ntrk40;
//         ntrk30rel = ntrk30 / ntrk40; 
//       }
//       if ( tree_nTracks == 0 ) {
//         ntrk10reltot = -0.02;
//         ntrk20reltot = -0.02;
//         ntrk30reltot = -0.02; 
//         ntrk40reltot = -0.02; 
//       }
//       else {
//         ntrk10reltot = ntrk10 / tree_nTracks;
//         ntrk20reltot = ntrk20 / tree_nTracks;
//         ntrk30reltot = ntrk30 / tree_nTracks; 
//         ntrk40reltot = ntrk40 / tree_nTracks; 
//       }
//$$$$$$

      // - Apply BDT -------------------
      bdtval = reader->EvaluateMVA( "BDTG" ); // default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
      // std::cout<<"tk pt : "<<pt<<" tk _eta : "<<eta<<" tk_phi :"<<phi<<" bdt_val :"<<bdtval<<std::endl;

      if ( dR < dRcut_tracks ) 
      {
        if ( isFromLLP == 1 ) LLP1_nTrks++;
        if ( isFromLLP == 2 ) LLP2_nTrks++;

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
      
      tree_track_Hemi_dR.push_back(dR);
      tree_track_Hemi_dRmax.push_back(dRmax);
      tree_track_Hemi.push_back(tracks_axis);
      if      ( tracks_axis == 1 ) tree_track_Hemi_LLP.push_back(iLLPrec1);
      else if ( tracks_axis == 2 ) tree_track_Hemi_LLP.push_back(iLLPrec2);
      else		           tree_track_Hemi_LLP.push_back(0);

      // float x_tk = tree_track_x[counter_track] - tree_PV_x;   
      // float y_tk = tree_track_y[counter_track] - tree_PV_y; 
      // float dr = abs(tree_track_dxy[counter_track]);
      // float drSig = tree_track_drSig[counter_track];
      // float phi_tk = tree_track_phi[counter_track];   

//$$$$$$
//       // linear intercept of track and its hemisphere axis in transverse plane 
//       float xT = (y_tk - x_tk*tan(phi_tk)) / (tan(vaxis1.Phi()) - tan(phi_tk));
//       float yT = xT * tan(vaxis1.Phi());
//       float d1 = dr;
//       float d1Sig = drSig;
//       if ( xT*cos(vaxis1.Phi()) + yT*sin(vaxis1.Phi()) < 0. ) {
//         d1    = -d1;
//         d1Sig = -d1Sig;
//       }
//       xT = (y_tk - x_tk*tan(phi_tk)) / (tan(vaxis2.Phi()) - tan(phi_tk));
//       yT = xT * tan(vaxis2.Phi());
//       float d2 = dr;
//       float d2Sig = drSig;
//       if ( xT*cos(vaxis2.Phi()) + yT*sin(vaxis2.Phi()) < 0. ) {
//         d2    = -d2;
//         d2Sig = -d2Sig;
//       }
//       if (tracks_axis == 1 ) {
//         tree_track_Hemi_d0.push_back(d1);
//         tree_track_Hemi_d0Sig.push_back(d1Sig);
//         tree_track_HemiOp_d0.push_back(d2);
//         tree_track_HemiOp_d0Sig.push_back(d2Sig);
//       }
//       else if (tracks_axis == 2 ) {
//         tree_track_Hemi_d0.push_back(d2);
//         tree_track_Hemi_d0Sig.push_back(d2Sig);
//         tree_track_HemiOp_d0.push_back(d1);
//         tree_track_HemiOp_d0Sig.push_back(d1Sig);
//       }
//       else {
//         tree_track_Hemi_d0.push_back(0);
//         tree_track_Hemi_d0Sig.push_back(0);
//         tree_track_HemiOp_d0.push_back(0);
//         tree_track_HemiOp_d0Sig.push_back(0);
//       }
//$$$$$$
      
    } //End loop on all the tracks
    
    //apply BDT event here or remove nTrks_axis1 and 2 from the BDTs ...

    mva_Evts_Hemi1_nTrks = nTrks_axis1;
    mva_Evts_Hemi2_nTrks = nTrks_axis2;
    //--------------------

    //-------Applying the EVTS BDT selection
    float EVTS_BDTval = -10;
    float EVTS_BDTvalDY = -10;
    float EVTS_BDTvalTT = -10;

    EVTS_BDTval = readerEvts->EvaluateMVA( "BDTGALLBKG" );
    EVTS_BDTvalDY =  readerEvts->EvaluateMVA("BDTGDY");
    EVTS_BDTvalTT =  readerEvts->EvaluateMVA("BDTGTT");

    // To select the events that you want, use  a macro (at the moment)
    tree_Evts_MVAval = EVTS_BDTval;
    tree_Evts_MVAvalDY = EVTS_BDTvalDY;
    tree_Evts_MVAvalTT = EVTS_BDTvalTT;

    // - Applyign the Hemi1 BDT selection
    // float HEMI1_BDTval = -10;
    // float HEMI1_BDTvalDY = -10;
    // float HEMI1_BDTvalTT = -10;

    // HEMI1_BDTval =    readerHemi1->EvaluateMVA("BDTGALLBKG" );//Current wp is -0.58 (optimizing S/sqrt(S+B)) but -0.7 is the point that is the maximum value to keep 90% of the signal)
    // HEMI1_BDTvalDY =  readerHemi1->EvaluateMVA("BDTGDY");
    // HEMI1_BDTvalTT =  readerHemi1->EvaluateMVA("BDTGTT");

    // tree_Hemi1_MVAval = HEMI1_BDTval;
    // tree_Hemi1_MVAvalDY = HEMI1_BDTvalDY;
    // tree_Hemi1_MVAvalTT = HEMI1_BDTvalTT;

    // - Applyign the Hemi2 BDT selection
    // float HEMI2_BDTval = -10;
    // float HEMI2_BDTvalDY = -10;
    // float HEMI2_BDTvalTT = -10;

    // HEMI2_BDTval =    readerHemi2->EvaluateMVA("BDTGALLBKG" );//Current wp is -0.58 (optimizing S/sqrt(S+B)) but -0.7 is the point that is the maximum value to keep 90% of the signal)
    // HEMI2_BDTvalDY =  readerHemi2->EvaluateMVA("BDTGDY");
    // HEMI2_BDTvalTT =  readerHemi2->EvaluateMVA("BDTGTT");

    // tree_Hemi2_MVAval = HEMI2_BDTval;
    // tree_Hemi2_MVAvalDY = HEMI2_BDTvalDY;
    // tree_Hemi2_MVAvalTT = HEMI2_BDTvalTT;


//$$$$
    //---------------------------------------------------------------//
    //----------- STORE TRACKS FROM V0s (AND SEC. INT.) -------------//
    //---------------------------------------------------------------//

    int countV0 = 0;
    for (unsigned int ipc = 0; ipc < TRACK_SIZE; ipc++) {
      pat::PackedCandidateRef pcref = MINIgeneralTracks[ipc];
      const reco::Track *trackPcPtr = pcref->bestTrack();
      if ( !trackPcPtr ) continue;

      // keep only tracks from K0 and Lambda decays
      bool isFromV0 = false, isFromSI = false;
      for (unsigned int j = 0 ; j < idxMGT.size() ; j++) 
      {
        if ( idxMGT[j].first && idxMGT[j].second == ipc ) isFromV0 = true;
        if ( idxSecIntMGT[j].first && idxSecIntMGT[j].second == ipc ) isFromSI = true;
      }
      if ( !isFromV0 && !isFromSI ) continue;
 
      reco::Track tk = *trackPcPtr;
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
      if ( tk_dzError > 0 )  tk_dzSig = abs(tk_dz) / tk_dzError; 
      if ( tk_dxyError > 0 ) tk_drSig = abs(tk_dxy) / tk_dxyError; 
      tk_vx = tk.vx();
      tk_vy = tk.vy();
      tk_vz = tk.vz();
      tk_px = tk.px();
      tk_py = tk.py();
      tk_pz = tk.pz();
      tk_HitPattern = tk.hitPattern();

      if ( tk_nHit == 0 ) continue;
      if ( tk_charge == 0 ) continue;

      // preselection
      if ( RequestHighPurity &&
           !( tk_pt > pt_Cut && tk_NChi2 < NChi2_Cut && tk_drSig > drSig_Cut
        && static_cast<int>(tk.quality(reco::TrackBase::highPurity))) ) continue;
      else if ( !( tk_pt > pt_Cut && tk_NChi2 < NChi2_Cut && tk_drSig > drSig_Cut ) ) continue;
 
      bool isLostTrack;
      if ( ipc < pc->size() ) isLostTrack = false;
      else                    isLostTrack = true;
      tree_V0_track_lost.push_back(     isLostTrack);
      tree_V0_track_isFromV0.push_back (isFromV0);
      tree_V0_track_isFromSI.push_back (isFromSI);

      tree_V0_track_pt.push_back       (tk_pt);
      tree_V0_track_eta.push_back      (tk_eta);
      tree_V0_track_phi.push_back      (tk_phi);
      tree_V0_track_charge.push_back   (tk_charge);
      tree_V0_track_NChi2.push_back    (tk_NChi2);
      tree_V0_track_dxy.push_back      (tk_dxy);
      tree_V0_track_drSig.push_back    (tk_drSig); 
      tree_V0_track_dz.push_back       (tk_dz);
      tree_V0_track_dzSig.push_back    (tk_dzSig);
      tree_V0_track_nHit.push_back     (tk_nHit);
      tree_V0_track_nHitPixel.push_back (tk_HitPattern.numberOfValidPixelHits());

      //----------------MINIAOD_Firsthit----------//
      const HitPattern hp = tk_HitPattern;
      uint16_t firsthit = hp.getHitPattern(HitPattern::HitCategory::TRACK_HITS,0);
 
      // Approximation for the lostTrack since the hitpattern information is not available (only 1160, tracking POG knows about it but do not seem to care)      
      if ( ipc >= pc->size() ) {
        if ( abs(tk_eta) < 1. && YEAR_ >= 2017) firsthit = 1184; // PIXBL4 in barrel
        else if (abs(tk_eta) < 1. && YEAR_ == 2016) firsthit = 1176;
        else		      firsthit = 1296; // PIXFD2 in forward
      }      
      tree_V0_track_firstHit.push_back(firsthit);

      //---Creating State to propagate from  TT---//
      BestTracks.push_back(theTransientTrackBuilder->build(tk));
      const MagneticField* B = BestTracks[countV0].field(); // 3.8T
      reco::TransientTrack TT (tk,BestTracks[countV0].field());
      // const FreeTrajectoryState Freetraj = TT.initialFreeState(); // Propagator in the barrel can also use FTS (WARNING: the so-called reference point (where the propagation starts might be different from the first vtx, a check should be done))
      GlobalPoint vert (tk_vx,tk_vy,tk_vz); // Point where the propagation will start (Reference Point)
      const TrajectoryStateOnSurface Surtraj = TT.stateOnSurface(vert); // TSOS of this point
      AnalyticalPropagator* Prop = new AnalyticalPropagator(B);
      Basic3DVector<float> P3D2(tk_vx,tk_vy,tk_vz);  // global frame
      Basic3DVector<float> B3DV (tk_px,tk_py,tk_pz); // global frame 
      float Eta = tk_eta;
      float Phi = tk_phi;
      float vz  = tk_vz;
      //------Propagation with new interface --> See ../interface/PropaHitPattern.h-----//
      
      std::pair<int,GloballyPositioned<float>::PositionType> FHPosition = PHP->Main(firsthit,Prop,Surtraj,Eta,Phi,vz,P3D2,B3DV);
      // returns 0 if in Barrel, 1 if disks with the position of the firsthit. Different propagators are used between 
      // barrel (StraightLinePlaneCrossing/geometry if propagator does not work)
      // and disks (HelixPlaneCrossing) 

      float xFirst = FHPosition.second.x();
      float yFirst = FHPosition.second.y();
      float zFirst = FHPosition.second.z();
      tree_V0_track_firstHit_x.push_back(xFirst);
      tree_V0_track_firstHit_y.push_back(yFirst);
      tree_V0_track_firstHit_z.push_back(zFirst);
      countV0++;
      //-----------------------END OF MINIAOD firsthit-----------------------//

      // track association to jet
      int iJet = 0;
      bool matchTOjet = false;
      for (const pat::Jet &jet : *jets) {
        if ( jet.pt() < jet_pt_min ) continue;

        float dR_jet = Deltar( jet.eta(), jet.phi(), tk_eta, tk_phi );
        if ( dR_jet < 0.4 ) {
          matchTOjet = true;
          break;
        }
        else iJet++;
      }
      if ( matchTOjet ) tree_V0_track_iJet.push_back (iJet);
      else	        tree_V0_track_iJet.push_back (-1);

      // dR between the track and the hemiqphere axes
      int tracks_axis = 0;
      dR1  = Deltar( tk_eta, tk_phi, axis1_eta, axis1_phi );
      dR2  = Deltar( tk_eta, tk_phi, axis2_eta, axis2_phi );
      if ( dR1 < dR2 &&  dR1 < dRcut_tracks  ) {
        tracks_axis = 1;
        dR = dR1;
        dRmax = dR2;
      }
      if ( dR2 < dR1 &&  dR2 < dRcut_tracks  ) {
        tracks_axis = 2;
        dR = dR2;
        dRmax = dR1;
      }
      tree_V0_track_Hemi_dR.push_back(dR);
      tree_V0_track_Hemi_dRmax.push_back(dRmax);
      tree_V0_track_Hemi.push_back(tracks_axis);

      // neighbour tracks
      ntrk10 = 0, ntrk20 = 0, ntrk30 = 0, ntrk40 = 0;
      float ntrk10_lost = 0, ntrk20_lost = 0, ntrk30_lost = 0, ntrk40_lost = 0;
      for (int counter_othertrack = 0; counter_othertrack < tree_nTracks; counter_othertrack++) 
      {
        if ( tree_track_ipc[counter_othertrack] == ipc ) continue;
        float x2 = tree_track_firstHit_x[counter_othertrack];
        float y2 = tree_track_firstHit_y[counter_othertrack];
        float z2 = tree_track_firstHit_z[counter_othertrack];
        float dist = TMath::Sqrt( (xFirst-x2)*(xFirst-x2) + (yFirst-y2)*(yFirst-y2) + (zFirst-z2)*(zFirst-z2) );
        if ( dist < 10. ) ntrk10++; 
        if ( dist < 20. ) ntrk20++;
        if ( dist < 30. ) ntrk30++;
        if ( dist < 40. ) ntrk40++;
        if ( isLostTrack && tree_track_lost[counter_othertrack] ) {
	  if ( dist < 10. ) ntrk10_lost++; 
  	  if ( dist < 20. ) ntrk20_lost++;
	  if ( dist < 30. ) ntrk30_lost++;
	  if ( dist < 40. ) ntrk40_lost++;
        }
      }  // end Loop on other Tracks
      if ( isLostTrack ) {
        ntrk10 = ntrk10_lost;
        ntrk20 = ntrk20_lost;
        ntrk30 = ntrk30_lost;
        ntrk40 = ntrk40_lost;
      }
      tree_V0_track_ntrk10.push_back(ntrk10);
      tree_V0_track_ntrk20.push_back(ntrk20);
      tree_V0_track_ntrk30.push_back(ntrk30);
      tree_V0_track_ntrk40.push_back(ntrk40);

    } // end loop on V0 (and SI) track candidates
//$$$$


    /////////////////////////////////////////
    // Sort tracks by decreasing BDT value //
    /////////////////////////////////////////

    for (int counter_track = 0; counter_track < tree_nTracks; counter_track++)
      MVAval[counter_track] = tree_track_MVAval[counter_track];
    if ( tree_nTracks < 10000 ) {
      TMath::Sort(tree_nTracks, MVAval, index);
    }
    else cout << " !!!!!! ERROR tree_nTracks " << tree_nTracks << " ERROR !!!!!! " << endl;


    ///////////////////////////////
    // Fill the transient tracks //
    ///////////////////////////////

    vector<reco::TransientTrack> displacedTracks_llp1_mva, displacedTracks_llp2_mva; // Control Tracks
    vector<reco::TransientTrack> displacedTracks_Hemi1_mva, displacedTracks_Hemi2_mva; // Tracks selected wrt the hemisphere
    vector<reco::TransientTrack> displacedTracks_Hemi1_mva_TW, displacedTracks_Hemi2_mva_TW; // Tracks selected wrt the hemisphere
    vector<reco::TransientTrack> displacedTracks_step2_Hemi1, displacedTracks_step2_Hemi2;
    
    vector<TLorentzVector> TrackInfo_llp1_mva;
    vector<TLorentzVector>TrackInfo_llp2_mva;

    vector<std::pair<float, TLorentzVector > > TrackInfo_Hemi1_mva;
    vector<std::pair<float, TLorentzVector > > TrackInfo_Hemi2_mva;
    vector<std::pair<float, TLorentzVector > > TrackInfo_step2_Hemi1;
    vector<std::pair<float, TLorentzVector > > TrackInfo_step2_Hemi2;

    vector<std::pair<bool,TLorentzVector>> Track_FirstHit_Hemi1_mva;
    vector<std::pair<bool,TLorentzVector>> Track_FirstHit_Hemi2_mva;
    vector<std::pair<bool,TLorentzVector>> Track_FirstHit_step2_Hemi1;
    vector<std::pair<bool,TLorentzVector>> Track_FirstHit_step2_Hemi2;

    for (int k = 0; k < tree_nTracks; k++)
    {
      int counter_track = index[k];
      unsigned int ipc = tree_track_ipc[counter_track];
      pat::PackedCandidateRef pcref = MINIgeneralTracks[ipc];
      const reco::Track *trackPcPtr = pcref->bestTrack();
      reco::Track tk;

      //PseudoDefTrack forces the covariance matrix to be positive definite. Can be negative in MINIAod with bestTrack/pseudoTrack methods
      //Correction defined by the BPH group => retrieve our effiency with this
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
        //       {
        //       std::cout<<" Apres Covcor m["<<i<<"]["<<j<<"] <<m[i][j]<<std::endl;
        //       }
        //   }
      //-------------------end covariance matrix correction-----------//
      }
      else tk = *trackPcPtr;

      int isFromLLP   = -1;
      if ( isMC_ ) isFromLLP = tree_track_sim_LLP[counter_track];
              
      int tracks_axis = tree_track_Hemi[counter_track];
      float track_p = sqrt(tree_track_px[counter_track]*tree_track_px[counter_track]+tree_track_py[counter_track]*tree_track_py[counter_track]+tree_track_pz[counter_track]*tree_track_pz[counter_track]);
      float track_e = tree_track_energy[counter_track];
      double bdtval   = tree_track_MVAval[counter_track];
      if ((tree_track_energy[counter_track]*tree_track_energy[counter_track]-track_p*track_p)<0)
        track_e = sqrt(track_p*track_p + 0.01948); // pi+ mass
              TLorentzVector vTrack(tree_track_px[counter_track],tree_track_py[counter_track],tree_track_pz[counter_track],track_e);
	      
	//  float x_IP =   tk.vx();   
	//  float y_IP =   tk.vy(); 
	//  float dPhiHemi1 = Deltaphi(tk.phi(),vaxis1.Phi());
	//  float dPhiHemi2 = Deltaphi(tk.phi(),vaxis2.Phi());
	//  float x_H1 = (y_IP-x_IP*tan(tk.phi()))/(tan(vaxis1.Phi())-tan(tk.phi()));
	//  float y_H1 = x_H1 * tan(vaxis1.Phi());
	//  float x_H2 = (y_IP-x_IP*tan(tk.phi()))/(tan(vaxis2.Phi())-tan(tk.phi()));
	//  float y_H2 = x_H2 * tan(vaxis2.Phi());
	//  float dxy_H1 = x_H1*cos(dPhiHemi1)+y_H1*sin(dPhiHemi1);// -x_H1*sin(dPhiHemi1)+y_H1*cos(dPhiHemi1);//using the CMS convention
	//  float dxy_H2 = x_H2*cos(dPhiHemi2)+y_H2*sin(dPhiHemi2);//-x_H2*sin(dPhiHemi2)+y_H2*cos(dPhiHemi2);//using the CMS convention

  //     //Computation of the IP of the tracks w.r.t the reco Axis
  //     if (tracks_axis == 1 )
  //     	{
	// 	tree_track_Hemi1_x.push_back(x_H1);
	// 	tree_track_Hemi2_x.push_back(x_H2);
	// 	tree_track_Hemi1_y.push_back(y_H1);
	// 	tree_track_Hemi2_y.push_back(y_H2);
	// 	tree_track_Hemi1_dxy.push_back(dxy_H1);
	// 	tree_track_Hemi2_dxy.push_back(dxy_H2);
	// }
	//       else if (tracks_axis == 2)
  //     	{
	// 	tree_track_Hemi1_x.push_back(x_H1);
	// 	tree_track_Hemi2_x.push_back(x_H2);
	// 	tree_track_Hemi1_y.push_back(y_H1);
	// 	tree_track_Hemi2_y.push_back(y_H2);
	// 	tree_track_Hemi1_dxy.push_back(dxy_H1);
	// 	tree_track_Hemi2_dxy.push_back(dxy_H2);
	// }
            //----------------------------------------------------------
      
      if ( bdtval > bdtcut ) { // tight wp
        if ( isFromLLP == 1 )
          {
            displacedTracks_llp1_mva.push_back(theTransientTrackBuilder->build(&tk));
            TLorentzVector SignalvTrack(tree_track_px[counter_track],tree_track_py[counter_track],tree_track_pz[counter_track],track_e);
            TrackInfo_llp1_mva.push_back(SignalvTrack);
          }
        if ( isFromLLP == 2 )
          {
            displacedTracks_llp2_mva.push_back(theTransientTrackBuilder->build(&tk));
            TLorentzVector SignalvTrack(tree_track_px[counter_track],tree_track_py[counter_track],tree_track_pz[counter_track],track_e);
            TrackInfo_llp2_mva.push_back(SignalvTrack);
          }

        if ( tracks_axis == 1 )
          {
            displacedTracks_Hemi1_mva.push_back(theTransientTrackBuilder->build(&tk));
            TrackInfo_Hemi1_mva.push_back(make_pair(tree_track_btag[counter_track],vTrack));
            TLorentzVector TrackFH(tree_track_firstHit_x[counter_track],tree_track_firstHit_y[counter_track],tree_track_firstHit_z[counter_track],0);
            Track_FirstHit_Hemi1_mva.push_back(make_pair(tree_track_lost[counter_track],TrackFH));
                      }
        if ( tracks_axis == 2 )
          {
            displacedTracks_Hemi2_mva.push_back(theTransientTrackBuilder->build(&tk));
            TrackInfo_Hemi2_mva.push_back(make_pair(tree_track_btag[counter_track],vTrack));
            TLorentzVector TrackFH(tree_track_firstHit_x[counter_track],tree_track_firstHit_y[counter_track],tree_track_firstHit_z[counter_track],0);
            Track_FirstHit_Hemi2_mva.push_back(make_pair(tree_track_lost[counter_track],TrackFH));
          }
      }

      if ( bdtval > bdtcut_step2  ) { // loose wp
        if ( tracks_axis == 1 )
          {
            displacedTracks_step2_Hemi1.push_back(theTransientTrackBuilder->build(&tk));
            TrackInfo_step2_Hemi1.push_back(make_pair(tree_track_btag[counter_track],vTrack));
            TLorentzVector TrackFH(tree_track_firstHit_x[counter_track],tree_track_firstHit_y[counter_track],tree_track_firstHit_z[counter_track],0);
            Track_FirstHit_step2_Hemi1.push_back(make_pair(tree_track_lost[counter_track],TrackFH));
          }    
        if ( tracks_axis == 2 )
          {
            displacedTracks_step2_Hemi2.push_back(theTransientTrackBuilder->build(&tk));
            TrackInfo_step2_Hemi2.push_back(make_pair(tree_track_btag[counter_track],vTrack));
            TLorentzVector TrackFH(tree_track_firstHit_x[counter_track],tree_track_firstHit_y[counter_track],tree_track_firstHit_z[counter_track],0);
            Track_FirstHit_step2_Hemi2.push_back(make_pair(tree_track_lost[counter_track],TrackFH));
          } 
      }
    }  // end loop on tracks
  
    // cout << " displaced tracks LLP1 " << LLP1_nTrks << " and with mva" << displacedTracks_llp1_mva.size() << endl;
    // cout << " displaced tracks LLP2 " << LLP2_nTrks << " and with mva" << displacedTracks_llp2_mva.size() << endl;
    // cout << " displaced tracks Hemi1 " << nTrks_axis1 << " and with mva" << displacedTracks_Hemi1_mva.size() << endl;
    // cout << " displaced tracks Hemi2 " << nTrks_axis2 << " and with mva" << displacedTracks_Hemi2_mva.size() << endl;


    ///////////////////////////////////////////////////////
    //-----------------------------------------------------
    // Vertex fitting 
    //-----------------------------------------------------
    ///////////////////////////////////////////////////////

    int Vtx_ntk = 0, Vtx_step = 0;
    float Vtx_x = 0., Vtx_y = 0., Vtx_z= 0., Vtx_chi = -10.;
    float recX, recY, recZ, dSV, recD;
 
    // parameters for the Adaptive Vertex Fitter (AVF)
    double maxshift        = 0.0001;
    unsigned int maxstep   = 30;
    double maxlpshift      = 0.1;
    double weightThreshold = 0.001;
    double sigmacut        = 3.;
    double Tini            = 256.;
    double ratio           = 0.25;
 
  
  if ( isMC_ && tree_nLLP > 0 ) {
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
    TLorentzVector Total4Vector1(0,0,0,0);
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
        float temp_px = 0;
        float temp_py = 0;
        float temp_pz = 0 ;
        float temp_e = 0 ;

        //--------------- B-tagging-----------------------//
            for (unsigned int i =  0 ; i < displacedTracks_llp1_mva.size(); i++)
              {
                temp_px = TrackInfo_llp1_mva[i].Px();
                temp_py = TrackInfo_llp1_mva[i].Py();
                temp_pz = TrackInfo_llp1_mva[i].Pz();
                temp_e  = TrackInfo_llp1_mva[i].E();
                TLorentzVector TLorentzTrack(temp_px,temp_py,temp_pz,temp_e);
                Total4Vector1 +=TLorentzTrack;
              }

    // std::cout<<"invaraint mass of LLP1 with LorentzVector : "<<sqrt(Total4Vector1.Mag2())<<std::endl;
      }
    }
    tree_LLP_Mass.push_back(Total4Vector1.Mag());
    tree_LLP.push_back(1);
    tree_LLP_pt.push_back(   LLP1_pt);
    tree_LLP_eta.push_back(  LLP1_eta);
    tree_LLP_phi.push_back(  LLP1_phi);
    tree_LLP_x.push_back(    LLP1_x);
    tree_LLP_y.push_back(    LLP1_y);
    tree_LLP_z.push_back(    LLP1_z);
    tree_LLP_r.push_back(sqrt(LLP1_x*LLP1_x+LLP1_y*LLP1_y));
    tree_LLP_dist.push_back( LLP1_dist);
    tree_LLP_nTrks.push_back(LLP1_nTrks);
    tree_LLP_Vtx_nTrks.push_back(Vtx_ntk);
    tree_LLP_Vtx_NChi2.push_back(Vtx_chi);
    tree_LLP_Vtx_dx.push_back(Vtx_x - LLP1_x);
    tree_LLP_Vtx_dy.push_back(Vtx_y - LLP1_y);
    tree_LLP_Vtx_dz.push_back(Vtx_z - LLP1_z);

    dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
    recX = Vtx_x - tree_PV_x;
    recY = Vtx_y - tree_PV_y;
    recZ = Vtx_z - tree_PV_z;
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
    Total4Vector1.SetPxPyPzE(0,0,0,0);
    if ( Vtx_ntk > 1 )
    {
      TransientVertex displacedVertex_llp2_mva = theFitter_vertex_llp2_mva.vertex(displacedTracks_llp2_mva); // fitted vertex
      
      if ( displacedVertex_llp2_mva.isValid() ) // NotValid if the max number of steps has been exceded or the fitted position is out of tracker bounds.
      {
        Vtx_x = displacedVertex_llp2_mva.position().x();
        Vtx_y = displacedVertex_llp2_mva.position().y();
        Vtx_z = displacedVertex_llp2_mva.position().z();
        Vtx_chi = displacedVertex_llp2_mva.normalisedChiSquared();

        float temp_px = 0;
        float temp_py = 0;
        float temp_pz = 0 ;
        float temp_e = 0 ;
        
        //--------------- B-tagging-----------------------//
            for (unsigned int i =  0 ; i < displacedTracks_llp2_mva.size(); i++)
              {
                temp_px = TrackInfo_llp2_mva[i].Px();
                temp_py = TrackInfo_llp2_mva[i].Py();
                temp_pz = TrackInfo_llp2_mva[i].Pz();
                temp_e  = TrackInfo_llp2_mva[i].E();
                TLorentzVector TLorentzTrack(temp_px,temp_py,temp_pz,temp_e);
                Total4Vector1 +=TLorentzTrack;
              }

        // std::cout<<"invaraint mass of LLP2 with LorentzVector : "<<sqrt(Total4Vector1.Mag2())<<std::endl;
      }
    }
    
    tree_LLP_Mass.push_back(Total4Vector1.Mag());
    tree_LLP.push_back(2);
    tree_LLP_pt.push_back(   LLP2_pt);
    tree_LLP_eta.push_back(  LLP2_eta);
    tree_LLP_phi.push_back(  LLP2_phi);
    tree_LLP_x.push_back(    LLP2_x);
    tree_LLP_y.push_back(    LLP2_y);
    tree_LLP_z.push_back(    LLP2_z);
    tree_LLP_dist.push_back( LLP2_dist);
    tree_LLP_r.push_back(sqrt(LLP2_x*LLP2_x+LLP2_y*LLP2_y));
    tree_LLP_nTrks.push_back(LLP2_nTrks);
    tree_LLP_Vtx_nTrks.push_back(Vtx_ntk);
    tree_LLP_Vtx_NChi2.push_back(Vtx_chi);
    tree_LLP_Vtx_dx.push_back(Vtx_x - LLP2_x);
    tree_LLP_Vtx_dy.push_back(Vtx_y - LLP2_y);
    tree_LLP_Vtx_dz.push_back(Vtx_z - LLP2_z);

    dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
    recX = Vtx_x - tree_PV_x;
    recY = Vtx_y - tree_PV_y;
    recZ = Vtx_z - tree_PV_z;
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_LLP_Vtx_dist.push_back( recD );
    tree_LLP_Vtx_dd.push_back( TMath::Sqrt(dSV)/LLP2_dist );
    
    float dR_LLP12   = Deltar(LLP1_eta, LLP1_phi, LLP2_eta, LLP2_phi);
    float deta_LLP12 = abs(LLP1_eta - LLP2_eta);
    float dphi_LLP12 = abs(Deltaphi(LLP1_phi, LLP2_phi));
    tree_LLP12_dR.push_back(   dR_LLP12);
    tree_LLP12_deta.push_back( deta_LLP12);
    tree_LLP12_dphi.push_back( dphi_LLP12);

//&&&&&
// //     if ( Vtx_chi < 0 && LLP2_nTrks > 1 ) dump = true;
//     cout << endl;
//     cout << " LLP " << 2 << " pt eta phi " << LLP2_pt << " " << LLP2_eta << " " << LLP2_phi << " x y z " << LLP2_x << " " << LLP2_y << " " << LLP2_z << " nTrks " << LLP2_nTrks << endl;
//     cout << "	  Vtx Chi2 " << Vtx_chi << " dx dy dz " << Vtx_x - LLP2_x << " " << Vtx_y - LLP2_y  << " " << Vtx_z - LLP2_z << " nTrks " << Vtx_ntk << endl;
//     cout << " &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& " << endl;
// //     cout << endl;
//&&&&&
  } // endif MC


      //-----------------------------------Vertexing-----------------------------------------------//
      //                               4 steps Vertexing                                           //
      // These 4 steps can be divided into 2 => (1,2) & (3,4)                                      //
      //        (1,2) : Corresponds to the Tight WP  of the Track level  BDT                       //
      //        (3,4) : Corresponds to the Loose WP of the Track level  BDT                        //
      //     => 1 & 3 are using the classic Adaptive Vertex Fitter (AVF) to build the vertices     //
      //     => 2 & 4 are using the Iterative Adaptive Vertex Fitter (IAVF)to build the vertices   //
      //          => Description of the IAVF at step 2                                             //
      //Process : Takes a collection of displaced track (from the BDT) as an input, then           //
      //          a vertex is built using the AVF or the IAVF if needed. No vertex can be obtained //
      //          out of these 4 steps                                                             //
      //                                                                                           //
      //Reconstruction Criteria : We require the vertex to have a chi2 per D.O.F between 0 and 10  //                                                                    
      //                          and to be matched to a generated vertex when looking at signal MC//
      //                                                                                           //
      //-------------------------------------------------------------------------------------------//
     
    //  Warning :  Sorry for the people reading the vertexing code (Vtx.h), it's not easy to read and to understand. I hope that the comments will be good enough :D
    // If you want to complain about this code, please contact : "Paul Vaucelle" <paul.vaucelle@cern.ch>.
    // If you come from Brittany (West of France), your remarks have higher chances of being accepted, thanks!

    //--------------------------- FIRST HEMISPHERE WITH MVA -------------------------------------//

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 1 step 1-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    Vtx_ntk = displacedTracks_Hemi1_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -10.;
    Vtx_step = 0;
    float SumWeight = 0;
    float DCA_VTX_Meand = 0;
    std::vector<float>  Vtx1_Weights;
    std::vector<unsigned int> Vtx1_index;
    int TightVertex = false;

    GlobalPoint PVPos(PV.x(),PV.y(),PV.z());
    Vtx* VtxHemi1 = new Vtx();
    // void Vertexing(std::vector<reco::TransientTrack> VertexTracks, vector<std::pair<bool,TLorentzVector>> Track_FirstHit, bool ActivateStep = true, bool RequireGoodChi2Seed = false,bool RequireGoodChi2VertexIter = false, float Chi2down = 0., float Chi2up = 10., GlobalPoint *PV = nullptr )
    // Track_FirstHit_Hemi1_mva.first tells you if the track is lost => don't car about the first hit of lost tracks.
    VtxHemi1->IAVFVertexing(displacedTracks_Hemi1_mva,Track_FirstHit_Hemi1_mva,ActivateStep1,false,false,0.,10.,&PVPos);
    if (VtxHemi1->chi2()>0 && VtxHemi1->chi2()<10)
      {
        Vtx_step = 1; // 1 here
        TightVertex = true; 
      }
    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 1 step 2-----------------------------------------//
    //------------------------------------------------------------------------------------------------//
    bool badVtx = false;
    if ( !(VtxHemi1->chi2()>0 && VtxHemi1->chi2()<10)  ) badVtx = true;
    if ( badVtx  && displacedTracks_Hemi1_mva.size() > 1 && (ActivateStep2 || IterAVF))
      {
        VtxHemi1->IAVFVertexing(displacedTracks_Hemi1_mva,Track_FirstHit_Hemi1_mva,ActivateStep2,true,false,0.,10.,&PVPos);
        Vtx_step = 2; 
        TightVertex = true;       
      }
    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 1 step 3-----------------------------------------//
    //------------------------------------------------------------------------------------------------//
    badVtx = false;
    if ( !(VtxHemi1->chi2()>0 && VtxHemi1->chi2()<10) ) badVtx = true;
    if ( badVtx  && displacedTracks_step2_Hemi1.size() > 1 && ActivateStep3)
      {
        VtxHemi1->IAVFVertexing(displacedTracks_step2_Hemi1,Track_FirstHit_step2_Hemi1,ActivateStep3,false,false,0.,10.,&PVPos);
        Vtx_step = 3; 
        TightVertex = false;
      }
    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 1 step 4-----------------------------------------//
    //------------------------------------------------------------------------------------------------//
    badVtx = false;
    if ( !(VtxHemi1->chi2()>0 && VtxHemi1->chi2()<10)  ) badVtx = true;
    if ( badVtx  && displacedTracks_step2_Hemi1.size() > 1 && (ActivateStep4 || IterAVF) )
      {
        VtxHemi1->IAVFVertexing(displacedTracks_step2_Hemi1,Track_FirstHit_step2_Hemi1,ActivateStep4,true,false,0.,10.,&PVPos);
        if ( VtxHemi1->chi2()>0 && VtxHemi1->chi2()<10 ) Vtx_step = 4;
//$$$$$$
        else  Vtx_step = 0;
//$$$$$$ 
        TightVertex = false;
      }

    Vtx_ntk = VtxHemi1->nTrk();
    Vtx_x = VtxHemi1->x();
    Vtx_y = VtxHemi1->y();
    Vtx_z = VtxHemi1->z();
    Vtx_chi = VtxHemi1->chi2();
    SumWeight = VtxHemi1->SumWeight();
    DCA_VTX_Meand = VtxHemi1->MeanDCA();
    Vtx1_Weights =  VtxHemi1->Vtx_TrkWeights();
    Vtx1_index = VtxHemi1->Vtx_TrkIndex();
    std::vector<TransientTrack> Vtx1_Trks = VtxHemi1->GetTTracks();
    GlobalError Vtx1posError = VtxHemi1->Vtx_PosErr();

    if ( VtxHemi1->chi2() == -10 ) Vtx1_Trks.clear();
    TLorentzVector Vtx1Vector(0,0,0,0);
    for (unsigned int i = 0 ; i <Vtx1_Weights.size(); i++)
      {
        tree_Hemi_Vtx_trackWeight.push_back(Vtx1_Weights[i]);
        if (Vtx1_Weights[i]>0.5)
          {
            const Track &tkm = Vtx1_Trks[i].track();            
            float tkpx = tkm.px();
            float tkpy = tkm.py();
            float tkpz = tkm.pz();
            float track_p = sqrt(tkpx*tkpx+tkpy*tkpy+tkpz*tkpz);
            float tke  = sqrt(track_p*track_p + 0.01948); // pi+ mass
            TLorentzVector TLorentzTrack(tkpx,tkpy,tkpz,tke);
            Vtx1Vector += TLorentzTrack;
          }
      }

    float Vtx1Mass = TMath::Max(Vtx1Vector.Mag(),0.);
    tree_Hemi_Vtx_Mass.push_back(Vtx1Mass);

  //------------------------------------------------------------------------------------------------//
  //------------------------------------------------------------------------------------------------//

    float Vtx_chi1 = Vtx_chi;
    tree_Hemi.push_back(1);
    tree_Hemi_njet.push_back(njet1);
    tree_Hemi_njet_nomu.push_back(njet1_nomu);
    tree_Hemi_pt.push_back(vaxis1.Pt());
    tree_Hemi_eta.push_back(axis1_eta);
    tree_Hemi_phi.push_back(axis1_phi);
    tree_Hemi_nTrks.push_back(nTrks_axis1);
    tree_Hemi_nTrks_sig.push_back(nTrks_axis1_sig);
    tree_Hemi_nTrks_bad.push_back(nTrks_axis1_bad);
    tree_Hemi_mass.push_back(  TMath::Max(vaxis1.Mag(),0.) );

    tree_HemiMu_mass.push_back( TMath::Max(Vobs1.Mag(),0.) );
    tree_HemiMu_pt.push_back( Vobs1.Pt() );
    tree_HemiMu_dR.push_back( dRobs1 );
    tree_HemiMuOp_mass.push_back( TMath::Max(VobsOp1.Mag(),0.) );
    tree_HemiMuOp_pt.push_back( VobsOp1.Pt() );
    tree_HemiMuOp_dR.push_back( dRobsOp1 );

    if ( ActivateONLYAVF )
      {
        TightVertex = false;
        VtxHemi1->AVFVertexing(displacedTracks_Hemi1_mva);
        if (VtxHemi1->chi2()>0 && VtxHemi1->chi2()<10)
          {
            Vtx_step = 1; 
            TightVertex = true; 
          }
        if ( !(VtxHemi1->chi2()>0 && VtxHemi1->chi2()<10)  ) badVtx = true;
        if ( badVtx  && displacedTracks_step2_Hemi1.size() > 1 && ActivateStep3)
          {
            VtxHemi1->AVFVertexing(displacedTracks_step2_Hemi1);
            Vtx_step = 3; 
//$$$$$$
	    if ( !(VtxHemi1->chi2()>0 && VtxHemi1->chi2()<10) ) Vtx_step = 0; 
//$$$$$$
            TightVertex = false;       
          }
        Vtx_x = VtxHemi1->x();
        Vtx_y = VtxHemi1->y();
        Vtx_z = VtxHemi1->z();
        Vtx_chi = VtxHemi1->chi2();
      }

    tree_Hemi_Vtx_step.push_back(Vtx_step);
    tree_Hemi_Vtx_isTight.push_back(TightVertex);
    tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
    tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
    tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis1_sig_mva);
    tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis1_bad_mva);
    tree_Hemi_Vtx_x.push_back(Vtx_x);
    tree_Hemi_Vtx_y.push_back(Vtx_y);
    tree_Hemi_Vtx_r.push_back(sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y));
    tree_Hemi_Vtx_z.push_back(Vtx_z);
    tree_Hemi_Vtx_xError.push_back(Vtx1posError.cxx());
    tree_Hemi_Vtx_yError.push_back(Vtx1posError.cyy());
    tree_Hemi_Vtx_zError.push_back(Vtx1posError.czz());

    recX = Vtx_x - tree_PV_x;
    recY = Vtx_y - tree_PV_y;
    recZ = Vtx_z - tree_PV_z;
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_Hemi_Vtx_dist.push_back( recD );
    tree_Hemi_Vtx_SumtrackWeight.push_back(SumWeight);
    tree_Hemi_Vtx_track_MeanDCA_d.push_back(DCA_VTX_Meand);

    int nVertex = 0;
    if ( Vtx_step>0 && Vtx_chi<10 && Vtx_chi>0 ) nVertex++;

  // -------------------------------------------//
    float posx1 = Vtx_x- tree_PV_x;
    float posy1 = Vtx_y- tree_PV_y;
    float posz1 = Vtx_z- tree_PV_z;
    float theta_Vtx1 = TMath::ATan2( sqrt(posx1*posx1+posy1*posy1) , abs(posz1) ) ;
    float eta_Vtx1 = -TMath::Log(tan(theta_Vtx1/2));
    if ( posz1 < 0 ) eta_Vtx1 = -eta_Vtx1;
    float phi1 = TMath::ATan2(posy1,posx1);
    float dRvtx1 = Deltar( eta_Vtx1, phi1, axis1_eta, axis1_phi );
    tree_Hemi_Vtx_dR.push_back(dRvtx1);
    
    int ping1 = 0; // 0 if no match, 1/2 if matched to LLP1/2, 3 if matched to both 
            // selection on distance between vertices for ping and merging
    float dSVcut = 0.1; // absolute cut (cm)
    float ddcut = 0.1;  // relative cut wrt decay length
    
  if ( isMC_ && tree_nLLP > 0 ) {
    tree_Hemi_LLP_dR.push_back(axis1_dR);
    if ( iLLPrec1 == 1 ) {
      tree_Hemi_LLP_pt.push_back( LLP1_pt);
      tree_Hemi_LLP_eta.push_back(LLP1_eta);
      tree_Hemi_LLP_phi.push_back(LLP1_phi);
      tree_Hemi_LLP_x.push_back(LLP1_x);
      tree_Hemi_LLP_y.push_back(LLP1_y);
      tree_Hemi_LLP_z.push_back(LLP1_z);
      tree_Hemi_LLP_dist.push_back(LLP1_dist);
      tree_Hemi_LLP_mother.push_back( LLP1_mother );
      tree_Hemi_LLP_Vtx_dx.push_back(Vtx_x - LLP1_x);
      tree_Hemi_LLP_Vtx_dy.push_back(Vtx_y - LLP1_y);
      tree_Hemi_LLP_Vtx_dz.push_back(Vtx_z - LLP1_z);
      tree_Hemi_LLP_Vtx_dr.push_back(TMath::Sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y) - TMath::Sqrt(LLP1_x*LLP1_x+LLP1_y*LLP1_y));
    }
    else {
      tree_Hemi_LLP_pt.push_back( LLP2_pt);
      tree_Hemi_LLP_eta.push_back(LLP2_eta);
      tree_Hemi_LLP_phi.push_back(LLP2_phi);
      tree_Hemi_LLP_x.push_back(LLP2_x);
      tree_Hemi_LLP_y.push_back(LLP2_y);
      tree_Hemi_LLP_z.push_back(LLP2_z);
      tree_Hemi_LLP_dist.push_back(LLP2_dist);
      tree_Hemi_LLP_mother.push_back( LLP2_mother );
      tree_Hemi_LLP_Vtx_dx.push_back(Vtx_x - LLP2_x);
      tree_Hemi_LLP_Vtx_dy.push_back(Vtx_y - LLP2_y);
      tree_Hemi_LLP_Vtx_dz.push_back(Vtx_z - LLP2_z);
      tree_Hemi_LLP_Vtx_dr.push_back(TMath::Sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y) - TMath::Sqrt(LLP2_x*LLP2_x+LLP2_y*LLP2_y));
    }
    tree_Hemi_LLP.push_back(iLLPrec1);

    dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
    float dSV1 = TMath::Sqrt(dSV)/LLP1_dist;
    if ( Vtx_chi > 0. && Vtx_chi < 10. && 
         (dSV < dSVcut || dSV1 < ddcut) ) ping1 = 1;
    dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
    float dSV2 = TMath::Sqrt(dSV)/LLP2_dist;
    if ( Vtx_chi > 0. && Vtx_chi < 10. && 
         (dSV < dSVcut || dSV2 < ddcut) ) ping1 += 2;
    
    // associate reconstructed muon and neutralino
    TLorentzVector vmuOK, vmuNO;
    float muOK_dR, muOK_pt, muOK_mass;
    float muNO_dR, muNO_pt, muNO_mass;
    int muOK, muNO;
    if ( iLLPrec1 == 1 ) {
      // be cautious: smuon ID 100013 has charge -1 !
      //if ( LLP1_mother * tree_muon_charge[imu1] < 0 ) {
if ( LLP1_mother * lep1_Q < 0 ) {
        muOK = imu1;
        muNO = imu2;
      }
      else {
        muOK = imu2;
        muNO = imu1;
      }
    }
    else {
      // be cautious: smuon ID 100013 has charge -1 !
      //if ( LLP2_mother * tree_muon_charge[imu1] < 0 ) {
if ( LLP2_mother * lep1_Q < 0 ) {
        muOK = imu1;
        muNO = imu2;
      }
      else {
        muOK = imu2;
        muNO = imu1;
      }
    }
    vmuOK.SetPtEtaPhiM(tree_muon_pt[muOK], tree_muon_eta[muOK], tree_muon_phi[muOK], mu_mass);
    vmuNO.SetPtEtaPhiM(tree_muon_pt[muNO], tree_muon_eta[muNO], tree_muon_phi[muNO], mu_mass);
        muOK_dR = Deltar( tree_muon_eta[muOK], tree_muon_phi[muOK], axis1_eta, axis1_phi );
    muNO_dR = Deltar( tree_muon_eta[muNO], tree_muon_phi[muNO], axis1_eta, axis1_phi );
        muOK_pt = (vmuOK+vaxis1).Pt();
    muNO_pt = (vmuNO+vaxis1).Pt();
    muOK_mass = TMath::Max((vmuOK+vaxis1).Mag(),0.);
    muNO_mass = TMath::Max((vmuNO+vaxis1).Mag(),0.);
    tree_Hemi_LLP_muOK_dR.push_back( muOK_dR );
    tree_Hemi_LLP_muNO_dR.push_back( muNO_dR );
    tree_Hemi_LLP_muOK_pt.push_back( muOK_pt );
    tree_Hemi_LLP_muNO_pt.push_back( muNO_pt );
    tree_Hemi_LLP_muOK_mass.push_back( muOK_mass );
    tree_Hemi_LLP_muNO_mass.push_back( muNO_mass );
  } // endif MC  

    // -------------------- End Of Invariant Mass ------------------------//

    float Vtx1_ntk = Vtx_ntk;
    float Vtx1_chi = Vtx_chi;
    float Vtx1_step = Vtx_step;
    // float Vtx1_r = sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y);
    // float Vtx1_z = Vtx_z;
    float Vtx1_STW = SumWeight;
    float Vtx1_Mass = Vtx1Mass ;
    float H1_Mass = Vobs1.Mag();
    // float Vtx1_dist = recD;
    float Vtx1_MeanDCA = DCA_VTX_Meand;
    if (Vtx1_Weights.size() != Vtx1_index.size() && showlog){ std::cout<<"size Vtx1_weights and Vtx1_index and  ntracks and chi and step: "<<Vtx1_Weights.size()<<" and "<<Vtx1_index.size()<<" and "<<Vtx_ntk<<" and "<<Vtx_chi<<" and "<<Vtx_step<<std::endl;}
    

    //--------------------------------------------------------------------------------------------//
    //--------------------------- SECOND HEMISPHERE WITH MVA -------------------------------------//
    //--------------------------------------------------------------------------------------------//

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 2 step 1-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    Vtx_ntk = displacedTracks_Hemi2_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -10.;
    Vtx_step = 0;
    SumWeight = 0;
    DCA_VTX_Meand = 0;
    std::vector<float>  Vtx2_Weights;
    std::vector<unsigned int>    Vtx2_index;
    TightVertex = false;
    // void Vertexing(std::vector<reco::TransientTrack> VertexTracks, vector<std::pair<bool,TLorentzVector>> Track_FirstHit, bool ActivateStep = true, bool RequireGoodChi2Seed = false,bool RequireGoodChi2VertexIter = false, float Chi2down = 0., float Chi2up = 10., GlobalPoint *PV = nullptr )//return type to be chnged
                
    // GlobalPoint PVPos(PV.x(),PV.y(),PV.z());
    Vtx* VtxHemi2 = new Vtx();
    VtxHemi2->IAVFVertexing(displacedTracks_Hemi2_mva,Track_FirstHit_Hemi2_mva,ActivateStep1,false,false,0.,10.,&PVPos);
    if (VtxHemi2->chi2()>0 && VtxHemi2->chi2()<10)
      {
        Vtx_step = 1;
        TightVertex = true; 
      }

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 2 step 2-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    badVtx = false;
    if ( !(VtxHemi2->chi2()>0 && VtxHemi2->chi2()<10)  ) badVtx = true;
    if ( badVtx  && displacedTracks_Hemi2_mva.size() > 1 && (ActivateStep2 || IterAVF))
      {
        VtxHemi2->IAVFVertexing(displacedTracks_Hemi2_mva,Track_FirstHit_Hemi2_mva,ActivateStep2,true,false,0.,10.,&PVPos);
        Vtx_step = 2;
        TightVertex = true;       
      }

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 2 step 3-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    badVtx = false;
    if ( !(VtxHemi2->chi2()>0 && VtxHemi2->chi2()<10) ) badVtx = true;
    if ( badVtx  && displacedTracks_step2_Hemi2.size() > 1 && ActivateStep3)
      {
        VtxHemi2->IAVFVertexing(displacedTracks_step2_Hemi2,Track_FirstHit_step2_Hemi2,ActivateStep3,false,false,0.,10.,&PVPos);
        Vtx_step = 3;
        TightVertex = false;
      }

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 2 step 4-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    badVtx = false;
    if ( !(VtxHemi2->chi2()>0 && VtxHemi2->chi2()<10)  ) badVtx = true;
    if ( badVtx  && displacedTracks_step2_Hemi2.size() > 1 && (ActivateStep4 || IterAVF) )
      {
        VtxHemi2->IAVFVertexing(displacedTracks_step2_Hemi2,Track_FirstHit_step2_Hemi2,ActivateStep4,true,false,0.,10.,&PVPos);
        if ( VtxHemi2->chi2()>0 && VtxHemi2->chi2()<10 ) Vtx_step = 4;
//$$$$$$
        else  Vtx_step = 0;
//$$$$$$ 
        TightVertex = false;
      }

    Vtx_ntk = VtxHemi2->nTrk();
    Vtx_x = VtxHemi2->x();
    Vtx_y = VtxHemi2->y();
    Vtx_z = VtxHemi2->z();
    Vtx_chi = VtxHemi2->chi2();
    SumWeight = VtxHemi2->SumWeight();
    DCA_VTX_Meand = VtxHemi2->MeanDCA();
    Vtx2_Weights =  VtxHemi2->Vtx_TrkWeights();
    Vtx2_index = VtxHemi2->Vtx_TrkIndex();
    std::vector<TransientTrack> Vtx2_Trks = VtxHemi2->GetTTracks();
    GlobalError Vtx2posError = VtxHemi2->Vtx_PosErr();

    if ( VtxHemi2->chi2() == -10 ) Vtx2_Trks.clear();
    TLorentzVector Vtx2Vector(0,0,0,0);
    for (unsigned int i = 0 ; i <Vtx2_Weights.size(); i++)
      {
        tree_Hemi_Vtx_trackWeight.push_back(Vtx2_Weights[i]);
        if (Vtx2_Weights[i]>0.5)
          {
            const Track &tkm = Vtx2_Trks[i].track();            
            float tkpx = tkm.px();
            float tkpy = tkm.py();
            float tkpz = tkm.pz();
            float track_p = sqrt(tkpx*tkpx+tkpy*tkpy+tkpz*tkpz);
            float tke  = sqrt(track_p*track_p + 0.01948); // pi+ mass
            TLorentzVector TLorentzTrack(tkpx,tkpy,tkpz,tke);
            Vtx2Vector += TLorentzTrack;
          }
      }

    float Vtx2Mass = TMath::Max(Vtx2Vector.Mag(),0.);
    tree_Hemi_Vtx_Mass.push_back(Vtx2Mass);

  //------------------------------------------------------------------------------------------------//
  //------------------------------------------------------------------------------------------------//

  float Vtx_chi2 = Vtx_chi;
  tree_Hemi.push_back(2);
  tree_Hemi_njet.push_back(njet2);
  tree_Hemi_njet_nomu.push_back(njet2_nomu);
  tree_Hemi_pt.push_back(vaxis2.Pt());
  tree_Hemi_eta.push_back(axis2_eta);
  tree_Hemi_phi.push_back(axis2_phi);
  tree_Hemi_nTrks.push_back(nTrks_axis2);
  tree_Hemi_nTrks_sig.push_back(nTrks_axis2_sig);
  tree_Hemi_nTrks_bad.push_back(nTrks_axis2_bad);
  tree_Hemi_mass.push_back(  TMath::Max(vaxis2.Mag(),0.) );
  tree_HemiMu_mass.push_back( TMath::Max(Vobs2.Mag(),0.) );
  tree_HemiMu_pt.push_back( Vobs2.Pt() );
  tree_HemiMu_dR.push_back( dRobs2 );
  tree_HemiMuOp_mass.push_back( TMath::Max(VobsOp2.Mag(),0.) );
  tree_HemiMuOp_pt.push_back( VobsOp2.Pt() );
  tree_HemiMuOp_dR.push_back( dRobsOp2 );

  if ( ActivateONLYAVF )
      {
        VtxHemi2->AVFVertexing(displacedTracks_Hemi2_mva);
        TightVertex = false;
        if (VtxHemi2->chi2()>0 && VtxHemi2->chi2()<10)
          {
            Vtx_step = 1; 
            TightVertex = true; 
          }
        if ( !(VtxHemi2->chi2()>0 && VtxHemi2->chi2()<10)  ) badVtx = true;
        if ( badVtx  && displacedTracks_step2_Hemi2.size() > 1 && ActivateStep3)
          {
            VtxHemi2->AVFVertexing(displacedTracks_step2_Hemi2);
            Vtx_step = 3;
//$$$$$$
	    if ( !(VtxHemi2->chi2()>0 && VtxHemi2->chi2()<10) ) Vtx_step = 0; 
//$$$$$$
            TightVertex = false;       
          }
        Vtx_x = VtxHemi2->x();
        Vtx_y = VtxHemi2->y();
        Vtx_z = VtxHemi2->z();
        Vtx_chi = VtxHemi2->chi2();
      }

  tree_Hemi_Vtx_step.push_back(Vtx_step);
  tree_Hemi_Vtx_isTight.push_back(TightVertex);
  tree_Hemi_Vtx_NChi2.push_back(Vtx_chi);
  tree_Hemi_Vtx_nTrks.push_back(Vtx_ntk);
  tree_Hemi_Vtx_nTrks_sig.push_back(nTrks_axis2_sig_mva);
  tree_Hemi_Vtx_nTrks_bad.push_back(nTrks_axis2_bad_mva);
  tree_Hemi_Vtx_x.push_back(Vtx_x);
  tree_Hemi_Vtx_y.push_back(Vtx_y);
  tree_Hemi_Vtx_r.push_back(sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y));
  tree_Hemi_Vtx_z.push_back(Vtx_z);
  tree_Hemi_Vtx_xError.push_back(Vtx2posError.cxx());
  tree_Hemi_Vtx_yError.push_back(Vtx2posError.cyy());
  tree_Hemi_Vtx_zError.push_back(Vtx2posError.czz());
  
  recX = Vtx_x - tree_PV_x;
  recY = Vtx_y - tree_PV_y;
  recZ = Vtx_z - tree_PV_z;
  recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
  tree_Hemi_Vtx_dist.push_back( recD );
  tree_Hemi_Vtx_SumtrackWeight.push_back(SumWeight);
  tree_Hemi_Vtx_track_MeanDCA_d.push_back(DCA_VTX_Meand);

  if ( Vtx_step>0 && Vtx_chi<10 && Vtx_chi>0 ) nVertex++;

  // -------------------------------------------//
  float posx2 = Vtx_x - tree_PV_x;
  float posy2 = Vtx_y - tree_PV_y;
  float posz2 = Vtx_z - tree_PV_z;
  float theta_Vtx2 = TMath::ATan2( sqrt(posx2*posx2+posy2*posy2) , abs(posz2) ) ;
  float eta_Vtx2 = -TMath::Log(tan(theta_Vtx2/2));
  if ( posz2 < 0 ) eta_Vtx2 = -eta_Vtx2;
  float phi2 = TMath::ATan2(posy2,posx2);
  float dRvtx2 = Deltar( eta_Vtx2, phi2, axis2_eta, axis2_phi );
  tree_Hemi_Vtx_dR.push_back(dRvtx2);

  int ping2 = 0; // 0 if no match, 1/2 if matched to LLP1/2, 3 if matched to both
  bool ping_Hemi1 = false, ping_Hemi2 = false;

  if ( isMC_ && tree_nLLP > 0 ) {
    tree_Hemi_LLP_dR.push_back(axis2_dR);
    if ( iLLPrec2 == 1 ) {
      tree_Hemi_LLP_pt.push_back( LLP1_pt);
      tree_Hemi_LLP_eta.push_back(LLP1_eta);
      tree_Hemi_LLP_phi.push_back(LLP1_phi);
      tree_Hemi_LLP_x.push_back(LLP1_x);
      tree_Hemi_LLP_y.push_back(LLP1_y);
      tree_Hemi_LLP_z.push_back(LLP1_z);
      tree_Hemi_LLP_dist.push_back(LLP1_dist);
      tree_Hemi_LLP_mother.push_back( LLP1_mother );
      tree_Hemi_LLP_Vtx_dx.push_back(Vtx_x - LLP1_x);
      tree_Hemi_LLP_Vtx_dy.push_back(Vtx_y - LLP1_y);
      tree_Hemi_LLP_Vtx_dz.push_back(Vtx_z - LLP1_z);
      tree_Hemi_LLP_Vtx_dr.push_back(TMath::Sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y) - TMath::Sqrt(LLP1_x*LLP1_x+LLP1_y*LLP1_y));
    }
    else {
      tree_Hemi_LLP_pt.push_back( LLP2_pt);
      tree_Hemi_LLP_eta.push_back(LLP2_eta);
      tree_Hemi_LLP_phi.push_back(LLP2_phi);
      tree_Hemi_LLP_x.push_back(LLP2_x);
      tree_Hemi_LLP_y.push_back(LLP2_y);
      tree_Hemi_LLP_z.push_back(LLP2_z);
      tree_Hemi_LLP_dist.push_back(LLP2_dist);
      tree_Hemi_LLP_mother.push_back( LLP2_mother );
      tree_Hemi_LLP_Vtx_dx.push_back(Vtx_x - LLP2_x);
      tree_Hemi_LLP_Vtx_dy.push_back(Vtx_y - LLP2_y);
      tree_Hemi_LLP_Vtx_dz.push_back(Vtx_z - LLP2_z);
      tree_Hemi_LLP_Vtx_dr.push_back(TMath::Sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y) - TMath::Sqrt(LLP2_x*LLP2_x+LLP2_y*LLP2_y));
    }
    tree_Hemi_LLP.push_back(iLLPrec2);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);

    // associate reconstructed muon and neutralino
    TLorentzVector vmuOK, vmuNO;
    float muOK_dR, muOK_pt, muOK_mass;
    float muNO_dR, muNO_pt, muNO_mass;
    int muOK, muNO;
    if ( iLLPrec2 == 1 ) {
      // be cautious: smuon ID 100013 has charge -1 !
      if ( LLP1_mother * lep1_Q < 0 ) {
        muOK = imu1;
        muNO = imu2;
      }
      else {
        muOK = imu2;
        muNO = imu1;
      }
    }
    else {
      // be cautious: smuon ID 100013 has charge -1 !
      if ( LLP2_mother * lep1_Q < 0 ) {
        muOK = imu1;
        muNO = imu2;
      }
      else {
        muOK = imu2;
        muNO = imu1;
      }
    }
    vmuOK.SetPtEtaPhiM(tree_muon_pt[muOK], tree_muon_eta[muOK], tree_muon_phi[muOK], mu_mass);
    vmuNO.SetPtEtaPhiM(tree_muon_pt[muNO], tree_muon_eta[muNO], tree_muon_phi[muNO], mu_mass);
        muOK_dR = Deltar( tree_muon_eta[muOK], tree_muon_phi[muOK], axis2_eta, axis2_phi );
    muNO_dR = Deltar( tree_muon_eta[muNO], tree_muon_phi[muNO], axis2_eta, axis2_phi );
        muOK_pt = (vmuOK+vaxis2).Pt();
    muNO_pt = (vmuNO+vaxis2).Pt();
    muOK_mass = TMath::Max((vmuOK+vaxis2).Mag(),0.);
    muNO_mass = TMath::Max((vmuNO+vaxis2).Mag(),0.);
    tree_Hemi_LLP_muOK_dR.push_back( muOK_dR );
    tree_Hemi_LLP_muNO_dR.push_back( muNO_dR );
    tree_Hemi_LLP_muOK_pt.push_back( muOK_pt );
    tree_Hemi_LLP_muNO_pt.push_back( muNO_pt );
    tree_Hemi_LLP_muOK_mass.push_back( muOK_mass );
    tree_Hemi_LLP_muNO_mass.push_back( muNO_mass );

    dSV = (Vtx_x - LLP1_x)*(Vtx_x - LLP1_x) + (Vtx_y - LLP1_y)*(Vtx_y - LLP1_y) + (Vtx_z - LLP1_z)*(Vtx_z - LLP1_z);
    float dSV1 = TMath::Sqrt(dSV)/LLP1_dist;
    if ( Vtx_chi > 0. && Vtx_chi < 10. && 
         (dSV < dSVcut || dSV1 < ddcut) ) ping2 = 1;
    dSV = (Vtx_x - LLP2_x)*(Vtx_x - LLP2_x) + (Vtx_y - LLP2_y)*(Vtx_y - LLP2_y) + (Vtx_z - LLP2_z)*(Vtx_z - LLP2_z);
    float dSV2 = TMath::Sqrt(dSV)/LLP2_dist;
    if ( Vtx_chi > 0. && Vtx_chi < 10. && 
         (dSV < dSVcut || dSV2 < ddcut) ) ping2 += 2;

    if ( (ping1 == 3 || ping2 == 3) && ping1*ping2 > 0 ) {
      ping_Hemi1 = true;
      ping_Hemi2 = true;
    }
    else if ( ping1 > 0 && ping2 == 0 ) {
      ping_Hemi1 = true;
    }
    else if ( ping2 > 0 && ping1 == 0 ) {
      ping_Hemi2 = true;
    }
    else if ( ping1 != ping2 ) {
      ping_Hemi1 = true;
      ping_Hemi2 = true;
    }
    else if ( ping1 == ping2 && ping1*ping2 > 0 ) {
      if      ( ping1 == iLLPrec1 ) {
        ping_Hemi1 = true;
      }
      else if ( ping2 == iLLPrec2 ) {
        ping_Hemi2 = true;
      }
      else if ( (Vtx1_step >= 1 && Vtx1_step <= 2) || Vtx_step > 2 ) {
        ping_Hemi1 = true;
      }
      else {
        ping_Hemi2 = true;
      }
    }

    tree_Hemi_LLP_ping.push_back( ping_Hemi1 );
    tree_Hemi_LLP_ping.push_back( ping_Hemi2 );
    int ping_event = 0;
    if      ( ping_Hemi1 && ping_Hemi2 ) ping_event = 2;
    else if ( ping_Hemi1 || ping_Hemi2 ) ping_event = 1;
    tree_event_LLP_ping.push_back( ping_event );
  } // endif MC
  
  tree_Hemi_dR12.push_back(dR_axis12);
  tree_Hemi_dR12.push_back(dR_axis12);
    
  //------------------------------------------------//

  float Vtx2_ntk = Vtx_ntk;
  float Vtx2_chi = Vtx_chi;
  float Vtx2_step = Vtx_step;
  // float Vtx2_r = sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y);
  // float Vtx2_z = Vtx_z;
  float Vtx2_STW = SumWeight;
  float Vtx2_Mass = Vtx2Mass;
  float H2_Mass = Vobs2.Mag();
  // float Vtx2_dist = recD;
  float Vtx2_MeanDCA = DCA_VTX_Meand;

  
  // -------------------------------------------------//
  // ---------------- CLOSE VERTICES -----------------//
  // -------------------------------------------------//

  float dr_2V = -10.;//quantities that should be lower for the backgrounds/ compared to signal
  float dz_2V = -10.;
  float dd_2V = -10.;
  float dRVtx = -10.;
  
  int MergeHemi = 1, SecHemi = 2;
  int MergeStep = 0, SecStep = -1;
  int pong_Merge = 0, pong_Sec = 0;
  bool pong_Hemi1 = false, pong_Hemi2 = false;

  if ( nVertex == 2 )
    {
        //Trying to deal with the two reconstructed vertices being really close <1mm
      dr_2V = sqrt((posx2-posx1)*(posx2-posx1)+(posy2-posy1)*(posy2-posy1));
      dz_2V = abs(posz2-posz1);
      dd_2V = sqrt(dr_2V*dr_2V + dz_2V*dz_2V);
      float recD1 = TMath::Sqrt(posx1*posx1 + posy1*posy1 + posz1*posz1);
      float recD2 = TMath::Sqrt(posx2*posx2 + posy2*posy2 + posz2*posz2);
      recD = (recD1 + recD2) / 2.;
      if ( Vtx1_step >= 1 && Vtx1_step <= 2 && Vtx2_step >= 3 ) recD = recD1; 
      if ( Vtx2_step >= 1 && Vtx2_step <= 2 && Vtx1_step >= 3 ) recD = recD2;
      
      dRVtx = Deltar(eta_Vtx1,phi1,eta_Vtx2,phi2); 
      tree_event_Vtx_Vtx_dr.push_back(dr_2V);
      tree_event_Vtx_Vtx_dz.push_back(dz_2V);
      tree_event_Vtx_Vtx_dd.push_back(dd_2V);
      tree_event_Vtx_Vtx_reldd.push_back(dd_2V/recD);
      tree_event_Vtx_Vtx_dR.push_back(dRVtx);

      int eventStep = 0;
      if      (  Vtx1_step >= 1 && Vtx1_step <= 2 
             &&  Vtx2_step >= 1 && Vtx2_step <= 2 )  eventStep = 1; 
      else if ( (Vtx1_step >= 1 && Vtx1_step <= 2) 
             || (Vtx2_step >= 1 && Vtx2_step <= 2) ) eventStep = 2;
      else eventStep = 3;

      tree_event_Vtx_Vtx_step.push_back(eventStep);

      if ( ActivateMerging && (dd_2V < dSVcut || dd_2V/recD < ddcut) ) 
      {  // criteria to activate the merging is in cm
         // The criteria that should be applied is that this dd_2V(dz and dr) distance should be below the resolution that we have on the vertices
         // Here are just some low  values for the two resolutions

          float Mergedx = -1000;
          float Mergedy = -1000;
          float Mergedz = -1000;
          float MergeVtx_ntk = 0;
          float MergeVtx_chi = -100;
          float MergedVtxMass = -1000;
          std::vector<float>  MergeVtx_Weights;
          std::vector<unsigned int>  MergeVtx_index;
          SumWeight=0;
                    
          if ( eventStep == 3 ) // both vtx are loose => final vtx is loose
            {
              Mergedx = (posx1+posx2)/2;
              Mergedy = (posy1+posy2)/2;
              Mergedz = (posz1+posz2)/2;
              MergeVtx_ntk = Vtx1_ntk+Vtx2_ntk;
              MergeVtx_chi =  (Vtx1_chi+Vtx2_chi)/2;
              SecStep = -4;
              SumWeight = (Vtx1_STW+Vtx2_STW);
              DCA_VTX_Meand = (Vtx1_MeanDCA+Vtx2_MeanDCA)/2;
              MergedVtxMass = (Vtx1Vector+Vtx2Vector).Mag();
            }
          else if ( eventStep == 2 ) //mix of tight and loose, we keep the tight vtx
            {
              bool isHemi1 = true;
              if ( Vtx1_step > 2 && Vtx2_step > 0 && Vtx2_step <= 2 ) isHemi1 = false;
              if ( isHemi1 )
                {
                  Mergedx = posx1;
                  Mergedy = posy1;
                  Mergedz = posz1;
                  MergeVtx_ntk = Vtx1_ntk;
                  MergeVtx_chi = Vtx1_chi;
                  SecStep = -3;
                  SumWeight = Vtx1_STW ;
                  DCA_VTX_Meand = Vtx1_MeanDCA;
                  MergedVtxMass = Vtx1_Mass;
                }
              else  
                {
                  Mergedx = posx2;
                  Mergedy = posy2;
                  Mergedz = posz2;
                  MergeVtx_ntk = Vtx2_ntk;
                  MergeVtx_chi = Vtx2_chi;
                  SecStep = -2;
                  SumWeight = Vtx2_STW ;
                  DCA_VTX_Meand = Vtx2_MeanDCA;
                  MergedVtxMass = Vtx2_Mass; 
                }
            }
          else if ( eventStep == 1 ) // both vtx are tight=> final vtx is tight
            {
              Mergedx = (posx1+posx2)/2;
              Mergedy = (posy1+posy2)/2;
              Mergedz = (posz1+posz2)/2;
              MergeVtx_ntk = Vtx1_ntk+Vtx2_ntk;
              MergeVtx_chi = (Vtx1_chi+Vtx2_chi)/2;
              SecStep = -1;
              SumWeight = (Vtx1_STW+Vtx2_STW);
              DCA_VTX_Meand = (Vtx1_MeanDCA+Vtx2_MeanDCA)/2;
              MergedVtxMass = (Vtx1Vector+Vtx2Vector).Mag();
            }

                  // Merged vertex is always valid, but just to be cautious...
            if ( MergeVtx_chi > 0. && MergeVtx_chi < 10. ) 
            {
            // std::cout<<"vtx1_chi: "<<Vtx1_chi<<std::endl;
            // std::cout<<"vtx2_chi: "<<Vtx2_chi<<std::endl;
            float thetaMerged = TMath::ATan2( sqrt(Mergedx*Mergedx+Mergedy*Mergedy) , abs(Mergedz) ) ;
            float MergedEta = -TMath::Log(tan(thetaMerged));
            if ( Mergedz < 0 ) MergedEta = -MergedEta;
            float MergedPhi = TMath::ATan2(Mergedy,Mergedx);
            recD = TMath::Sqrt(Mergedx*Mergedx + Mergedy*Mergedy + Mergedz*Mergedz);
            Mergedx += tree_PV_x;
            Mergedy += tree_PV_y;
            Mergedz += tree_PV_z;// becasue we comptued the Merged position takign into acount the PC position !!!
            float Mergedr = sqrt(Mergedx*Mergedx+Mergedy*Mergedy);

            //----------------------------------------------------
            // --  Check the closest Hemisphere to this pseudo-vertex
            //----------------------------------------------------
            float dR_axis1 = Deltar(MergedEta, MergedPhi, axis1_eta, axis1_phi);
            float dR_axis2 = Deltar(MergedEta, MergedPhi, axis2_eta, axis2_phi);
            float dRmerge = dR_axis1;
            bool isAxis1 = true;
            if (dR_axis2 < dR_axis1) isAxis1 = false;
            if ( isAxis1 ) {
              MergeHemi = 1;
              MergeStep = Vtx1_step;
            }
            else {
              MergeHemi = 2;
              MergeStep = Vtx2_step;
              dRmerge = dR_axis2;
            }

            if ( isMC_ && tree_nLLP > 0 ) {
              if ( MergeHemi == 1 ) {
                tree_Hemi_SecLLP.push_back( iLLPrec1 );
	            }
              if ( MergeHemi == 2 ) {
                tree_Hemi_SecLLP.push_back( iLLPrec2 );
	            }
	            if ( ((MergeHemi == 1 && iLLPrec1 == 1) || 
	                  (MergeHemi == 2 && iLLPrec2 == 1)) ) {
                tree_Hemi_LLP_SecVtx_dx.push_back(Mergedx - LLP1_x);
                tree_Hemi_LLP_SecVtx_dy.push_back(Mergedy - LLP1_y);
                tree_Hemi_LLP_SecVtx_dz.push_back(Mergedz - LLP1_z);
                tree_Hemi_LLP_SecVtx_dr.push_back(TMath::Sqrt(Mergedr*Mergedr) - TMath::Sqrt(LLP1_x*LLP1_x+LLP1_y*LLP1_y));
	            }
              if ( ((MergeHemi == 1 && iLLPrec1 == 2) ||
	                  (MergeHemi == 2 && iLLPrec2 == 2)) ) {
                tree_Hemi_LLP_SecVtx_dx.push_back(Mergedx - LLP2_x);
                tree_Hemi_LLP_SecVtx_dy.push_back(Mergedy - LLP2_y);
                tree_Hemi_LLP_SecVtx_dz.push_back(Mergedz - LLP2_z);
                tree_Hemi_LLP_SecVtx_dr.push_back(TMath::Sqrt(Mergedr*Mergedr) - TMath::Sqrt(LLP2_x*LLP2_x+LLP2_y*LLP2_y));
	            }
	          
              dSV = (Mergedx - LLP1_x)*(Mergedx - LLP1_x) + (Mergedy - LLP1_y)*(Mergedy - LLP1_y) + (Mergedz - LLP1_z)*(Mergedz - LLP1_z);
              float dSV1 = TMath::Sqrt(dSV)/LLP1_dist;
              if ( (dSV < dSVcut || dSV1 < ddcut) ) pong_Merge = 1;
              dSV = (Mergedx - LLP2_x)*(Mergedx - LLP2_x) + (Mergedy - LLP2_y)*(Mergedy - LLP2_y) + (Mergedz - LLP2_z)*(Mergedz - LLP2_z);
              float dSV2 = TMath::Sqrt(dSV)/LLP2_dist;
              if ( (dSV < dSVcut || dSV2 < ddcut) ) pong_Merge += 2;
              } // endif MC

                //----------------------------------------------------
                // -- Use the whole set of tracks available from the hemisphere + remove the ones used in the previous vertex 
                //----------------------------------------------------

            unsigned int count_trk = 0;
                // CHeck for step 1 and 2 first ...
            const unsigned int size1 = displacedTracks_Hemi1_mva.size();
            const unsigned int size2 = displacedTracks_Hemi2_mva.size();
            const unsigned int size3 = displacedTracks_step2_Hemi1.size();
            const unsigned int size4 = displacedTracks_step2_Hemi2.size();

              int IndextoRemove1[size1]={0} ;
              int IndextoRemove2[size2]={0};
              int IndextoRemove1_2[size3]={0} ;
              int IndextoRemove2_2[size4]={0} ;
              vector<std::pair<bool,TLorentzVector>> FHNewTTrack;
              bool NewTightVertex = true;
              std::vector<TransientTrack> NewTTracks;
              std::vector<TransientTrack> TempTT;
              
          if ( NewTightVertex )
                {
                  if (Vtx1_step == 1 || Vtx1_step == 2)
                      {
                        //Hemisphere 1

                        //the logic is:
                        // from the whole set of tracks available for the steps 1 aand 2 of Hemi1, we loop over the tracks that are effectively used (Vtx1_index loop)
                        //But one has to becareful since the weight of tracks assoicated to the vtx can be low=> one can use those tracks to form a vtx
                        for(std::vector<TransientTrack>::iterator it = displacedTracks_Hemi1_mva.begin(); it != displacedTracks_Hemi1_mva.end();)// On the whole set of tracks used for steps 1 and 2 for hemi1...
                          {
                            for(unsigned int j = 0; j < Vtx1_index.size(); j++ )//. we look at the ones that are effectively used
                            {
                                // std::cout<<"count_trk : "<<count_trk<<" and Vtx1_index j :"<<Vtx1_index[j]<<std::endl;
                                if (count_trk == Vtx1_index[j])
                                  {
                                    if (Vtx1_Weights[j]<0.05)// keep tracks that did not work for the merged vtx
                                    //We can make the hypothesis that tracks that have a weight > 0.5 would also have aweight above 0.5
                                    // ( not exactly true but one can think that we look for tracks that are not near the merged vtx therefore
                                    // not near the first two vertices)
                                      {
                                        break;
                                      }
                                    IndextoRemove1[count_trk]++;
                                    break;
                                  }
                            }
                            count_trk++;
                            if (count_trk==displacedTracks_Hemi1_mva.size()) break;
                          }
                      }
                    if (Vtx2_step == 1 || Vtx2_step == 2)
                      {
                         // Hemisphere 2 
                        count_trk = 0;
                        for(std::vector<TransientTrack>::iterator it = displacedTracks_Hemi2_mva.begin(); it != displacedTracks_Hemi2_mva.end(); )
                          {
                            for(unsigned int j = 0; j < Vtx2_index.size(); j++ )
                            {
                              // std::cout<<"count_trk : "<<count_trk<<" and Vtx2_index j :"<<Vtx2_index[j]<<std::endl;
                                if (count_trk == Vtx2_index[j])
                                  {
                                    if (Vtx2_Weights[j]<0.05) // keep tracks that di not work for the merged vtx
                                      {
                                        break;
                                      }
                                    IndextoRemove2[count_trk]++ ;
                                    break;
                                  }
                            }
                            count_trk++;
                            if (count_trk==displacedTracks_Hemi2_mva.size()) break;
                          }
                      }

                    if (!isAxis1 && (SecStep == -1 ))// we take the tracks from the hemisphere the further away from the vtx <=> in the other hemisphere
                      {
                        for (unsigned int i = 0 ; i < size1 ; i++)
                          {
                            if(IndextoRemove1[i]==0)
                              {
                                    NewTTracks.push_back(displacedTracks_Hemi1_mva[i]);
                                    FHNewTTrack.push_back(Track_FirstHit_Hemi1_mva[i]);
                              }
                          }
                      }
                  
                    if (isAxis1 && (SecStep == -1))//steps -1  
                      {
                        for (unsigned int i = 0 ; i < size2 ; i++)
                          {
                            if(IndextoRemove2[i]==0)
                              {
                                NewTTracks.push_back(displacedTracks_Hemi2_mva[i]);
                                FHNewTTrack.push_back(Track_FirstHit_Hemi2_mva[i]);
                              }
                          }
                      }
                } // endif newTight vertex

                // GlobalPoint PVPos(PV.x(),PV.y(),PV.z());
                SecStep = 0;
                Vtx* VtxSec = new Vtx();
                // void Vertexing(std::vector<reco::TransientTrack> VertexTracks, vector<std::pair<bool,TLorentzVector>> Track_FirstHit, bool ActivateStep = true, bool RequireGoodChi2Seed,bool RequireGoodChi2VertexIter = false, float Chi2down = 0., float Chi2up = 10.)
                VtxSec->IAVFVertexing(NewTTracks,FHNewTTrack,true,false,false,0.,10.,&PVPos);
                if (VtxSec->chi2() > 0. && VtxSec->chi2() < 10.) 
                  {
                    SecStep = 1;
                  }

                badVtx = false;
                if ( !(VtxSec->chi2() > 0. && VtxSec->chi2() < 10.) ) badVtx = true;
                if ( badVtx  && NewTTracks.size() > 1 && (ActivateStep2 || IterAVF))
                  {
                    VtxSec->IAVFVertexing(NewTTracks,FHNewTTrack,true,true,false,0.,10.,&PVPos);
                    SecStep = 2;
                  }

                // then steps 3 and 4 =< Loose 
                NewTTracks.clear();
                FHNewTTrack.clear();
                if ( (!(VtxSec->chi2() > 0. && VtxSec->chi2() < 10.)) || NewTTracks.size()<2) // || NewTTracks.size()<2
                  {
                    count_trk = 0;
                    if (Vtx1_step == 3 || Vtx1_step == 4)
                      {
                        //Hemisphere 1
                        for(std::vector<TransientTrack>::iterator it = displacedTracks_step2_Hemi1.begin(); it != displacedTracks_step2_Hemi1.end(); )
                          {
                            for(unsigned int j = 0; j < Vtx1_index.size(); j++ )
                            {
                                // std::cout<<"count_trk : "<<count_trk<<" and Vtx1_index j :"<<Vtx1_index[j]<<std::endl;
                                if (count_trk == Vtx1_index[j])
                                  {
                                    if (Vtx1_Weights[j]<0.05)// keep tracks that did not work for the merged vtx
                                    //We can make the hypothesis that tracks that have a weight > 0.5 would also have aweight above 0.5
                                    // ( not exactly true but one can think that we look for tracks that are not near the merged vtx therefore
                                    // not near the first two vertices)
                                      {
                                        break;
                                      }
                                    IndextoRemove1_2[count_trk]++;
                                  }
                            }
                            count_trk++;
                            if (count_trk==displacedTracks_step2_Hemi1.size()) break;
                          }
                      }
                    if (Vtx2_step == 3 || Vtx2_step == 4)
                      {
                          // Hemisphere 2 
                        count_trk = 0;
                        for(std::vector<TransientTrack>::iterator it = displacedTracks_step2_Hemi2.begin(); it != displacedTracks_step2_Hemi2.end();)
                          {
                            for(unsigned int j = 0; j < Vtx2_index.size(); j++ )
                            {
                              // std::cout<<"count_trk : "<<count_trk<<" and Vtx2_index j :"<<Vtx2_index[j]<<std::endl;
                                if (count_trk == Vtx2_index[j])
                                  {
                                    if (Vtx2_Weights[j]<0.05) // keep tracks that di not work for the merged vtx
                                      {
                                        break;
                                      }
                                    IndextoRemove2_2[count_trk]++;
                                  }
                            }
                            count_trk++;
                            if (count_trk==displacedTracks_step2_Hemi2.size()) break;
                          }
                      }

                    if (!isAxis1)// all the SecStep go through this => For SecStep == -1, because we might no have enought tracks from step1 and 2
                      {
                        for (unsigned int i = 0 ; i < size3 ; i++)
                          {
                            if(IndextoRemove1_2[i]==0)
                              {
                                NewTTracks.push_back(displacedTracks_step2_Hemi1[i]);
                                FHNewTTrack.push_back(Track_FirstHit_step2_Hemi1[i]);
                              }
                          }
                      }

                    if (isAxis1 )// all the SecStep go through this => For SecStep == -1, because we might no have enought tracks from step1 and 2
                      {
                        for (unsigned int i = 0 ; i < size4 ; i++)
                          {
                            if(IndextoRemove2_2[i]==0)
                              {
                                NewTTracks.push_back(displacedTracks_step2_Hemi2[i]);
                                FHNewTTrack.push_back(Track_FirstHit_step2_Hemi2[i]);
                              }
                          }
                      }
                  }

                badVtx = false;
                if ( !(VtxSec->chi2() > 0. && VtxSec->chi2() < 10.) ) badVtx = true;
                NewTightVertex = false;
                if ( badVtx  && NewTTracks.size() > 1 && ActivateStep3 && !NewTightVertex)
                  {
                    VtxSec->IAVFVertexing(NewTTracks,FHNewTTrack,true,false,false,0.,10.,&PVPos);
                    SecStep = 3;
                  }
                badVtx = false;
                if ( !(VtxSec->chi2() > 0. && VtxSec->chi2() < 10.) ) badVtx = true;
                if ( badVtx  && NewTTracks.size() > 1 && (ActivateStep4 || IterAVF) && !NewTightVertex)
                  {
                    VtxSec->IAVFVertexing(NewTTracks,FHNewTTrack,true,true,false,0.,10.,&PVPos);
                    //$$
                    if ( VtxSec->chi2() > 0. && VtxSec->chi2() < 10.) SecStep = 4; // missing ?
                    //$$$$$$
                                  else SecStep = 0;
                    //$$$$$$ 
                  }

                  float SecVtx_x = VtxSec->x();
                  float SecVtx_y = VtxSec->y();
                  float SecVtx_z = VtxSec->z();
                  float SecVtx_ntk = VtxSec->nTrk();
                  float SecVtx_chi = VtxSec->chi2();
                  std::vector<TransientTrack> SecVtx_Trks = VtxSec->GetTTracks();
                  std::vector<float>  NewVtx_Weights = VtxSec->Vtx_TrkWeights();
                  std::vector<unsigned int>  NewVtx_index = VtxSec->Vtx_TrkIndex();
                  float SecDCA_VTX_Meand = VtxSec->MeanDCA();
                  float SecSumWeight	 = VtxSec->SumWeight();

                  TLorentzVector NewVtxVector(0,0,0,0);
                  //$$          if ( SecVtx_chi== -10 && SecStep == 4 ) SecVtx_Trks.clear();
                  if ( SecVtx_chi== -10 ) SecVtx_Trks.clear();
                  for (unsigned int i = 0 ; i <SecVtx_Trks.size(); i++)
                    {
                      tree_Hemi_SecVtx_trackWeight.push_back(NewVtx_Weights[i]);
                      if (NewVtx_Weights[i]>0.5)
                        {
                          const Track &tkm = SecVtx_Trks[i].track();            
                          float tkpx = tkm.px();
                          float tkpy = tkm.py();
                          float tkpz = tkm.pz();
                          float track_p = sqrt(tkpx*tkpx+tkpy*tkpy+tkpz*tkpz);
                          float tke  = sqrt(track_p*track_p + 0.01948); // pi+ mass
                          TLorentzVector TLorentzTrack(tkpx,tkpy,tkpz,tke);
                          NewVtxVector += TLorentzTrack;
                        }
                    }
                  
                            float SecrecX = SecVtx_x - tree_PV_x;
                            float SecrecY = SecVtx_y - tree_PV_y;
                            float SecrecZ = SecVtx_z - tree_PV_z;
                            float SecrecD = TMath::Sqrt(SecrecX*SecrecX + SecrecY*SecrecY + SecrecZ*SecrecZ);
                        
                        //-> In this section: we have in fact two newly built vertices with a good chi2
                        //Check that the two newly built vertices are not "clsoe" from each other <=> far from the millimeter scale and that the hopefully,
                        // the delta R between the two vertices is going towards 3.14 :s

                        if ( SecVtx_chi > 0. && SecVtx_chi < 10. )
                          {
                            float New_dr_2V = sqrt((Mergedx-SecVtx_x)*(Mergedx-SecVtx_x)+(Mergedy-SecVtx_y)*(Mergedy-SecVtx_y));
                            float New_dz_2V = abs(Mergedz-SecVtx_z);
                            float New_dd_2V = sqrt(New_dr_2V*New_dr_2V + (Mergedz-SecVtx_z)*(Mergedz-SecVtx_z));
                            tree_event_MergedVtx_Vtx_dr.push_back(New_dr_2V);
                            tree_event_MergedVtx_Vtx_dz.push_back(New_dz_2V);
                            tree_event_MergedVtx_Vtx_dd.push_back(New_dd_2V);
                            tree_event_MergedVtx_Vtx_reldd.push_back(New_dd_2V/SecrecD);
                            float theta_SecVtx = TMath::ATan2(sqrt(SecrecX*SecrecX+SecrecY*SecrecY),abs(SecrecZ)) ;
                            float eta_SecVtx = -TMath::Log(tan(theta_SecVtx/2));
                            if ( SecVtx_z < 0 ) eta_SecVtx = -eta_SecVtx;
                            float SecPhi = TMath::ATan2(SecrecY,SecrecX);
                            float New_dRVtx = Deltar(MergedEta,MergedPhi,eta_SecVtx,SecPhi);
                            tree_event_MergedVtx_Vtx_dR.push_back(New_dRVtx);

                            int eventNewStep = 0;
                            if ( MergeStep >= 1 && MergeStep <= 2 
                                && SecStep >= 1 && SecStep <= 2 ) eventNewStep = 1; 
                            else if ( (MergeStep >= 1 && MergeStep <= 2) 
                                    || (SecStep >= 1 && SecStep <= 2) ) eventNewStep = 2;
                            else eventNewStep = 3;
                            
                            if (New_dd_2V < dSVcut || New_dd_2V/SecrecD < ddcut) // re-merging
                              {
                                  SecStep = 0;
                                  eventNewStep = -eventNewStep;
                                  if ( eventNewStep == -3 ) // both vtx are loose => final vtx is loose
                                    {
                                      Mergedx = (Mergedx+SecVtx_x)/2;
                                      Mergedy = (Mergedy+SecVtx_y)/2;
                                      Mergedz = (Mergedz+SecVtx_z)/2;
                                      MergeVtx_ntk = MergeVtx_ntk+SecVtx_ntk;
                                      MergeVtx_chi = (MergeVtx_chi+SecVtx_chi)/2;
                                      SumWeight = SumWeight+SecSumWeight;
                                      DCA_VTX_Meand = (DCA_VTX_Meand+SecDCA_VTX_Meand)/2;
                                      MergedVtxMass = (Vtx1Vector+Vtx2Vector+NewVtxVector).Mag();
                                                                          }
                                  else if ( eventNewStep == -2 ) // mix of tight and loose, we keep the tight vtx
                                    {
                                      if ( MergeStep > 2 && SecStep > 0 && SecStep <= 2 ) 
                                        {
                                          Mergedx = SecVtx_x;
                                          Mergedy = SecVtx_y;
                                          Mergedz = SecVtx_z;
                                          MergeVtx_ntk = SecVtx_ntk;
                                          MergeVtx_chi = SecVtx_chi;
                                          SumWeight = SecSumWeight ;
                                          DCA_VTX_Meand = SecDCA_VTX_Meand;
                                          MergedVtxMass = NewVtxVector.Mag(); 
                                        }
                                        //$$
                                            // else do nothing (MergeVtx informations still valid)
                                        //$$
                                    }
                                  else if ( eventNewStep == -1 ) // both vtx are tight=> final vtx is tight
                                    {
                                      Mergedx = (Mergedx+SecVtx_x)/2;
                                      Mergedy = (Mergedy+SecVtx_y)/2;
                                      Mergedz = (Mergedz+SecVtx_z)/2;
                                      MergeVtx_ntk = MergeVtx_ntk+SecVtx_ntk;
                                      MergeVtx_chi = (MergeVtx_chi+SecVtx_chi)/2;
                                      SumWeight = SumWeight+SecSumWeight;
                                      DCA_VTX_Meand = (DCA_VTX_Meand+SecDCA_VTX_Meand)/2;
                                      MergedVtxMass = (Vtx1Vector+Vtx2Vector+NewVtxVector).Mag();
                                      }
                              } // endif re-merging

                            tree_event_MergedVtx_Vtx_step.push_back(eventNewStep);
                            } // endif SecVtx is valid

                        //---------Merged Vertex-----//
                        tree_Hemi_SecVtx.push_back(MergeHemi);
                        tree_Hemi_SecVtx_step.push_back(MergeStep);
                        tree_Hemi_SecVtx_x.push_back(Mergedx);
                        tree_Hemi_SecVtx_y.push_back(Mergedy);
                        tree_Hemi_SecVtx_z.push_back(Mergedz);
                        tree_Hemi_SecVtx_r.push_back(Mergedr);
                        tree_Hemi_SecVtx_nTrks.push_back(MergeVtx_ntk);
                        tree_Hemi_SecVtx_NChi2.push_back(MergeVtx_chi);
                        tree_Hemi_SecVtx_dist.push_back( recD );
                        tree_Hemi_SecVtx_track_MeanDCA_d.push_back(DCA_VTX_Meand);
                        tree_Hemi_SecVtx_SumtrackWeight.push_back(SumWeight);         
                        tree_Hemi_SecVtx_Mass.push_back(MergedVtxMass);
                        tree_Hemi_SecVtx_dR.push_back(dRmerge); 

                        //-------New secondary vertex---//
                        SecHemi = 2;
		                    if ( MergeHemi == 2 ) SecHemi = 1;
                        tree_Hemi_SecVtx.push_back(SecHemi);
                        tree_Hemi_SecVtx_step.push_back(SecStep);
                        tree_Hemi_SecVtx_x.push_back(SecVtx_x);
                        tree_Hemi_SecVtx_y.push_back(SecVtx_y);
                        tree_Hemi_SecVtx_z.push_back(SecVtx_z);
                        tree_Hemi_SecVtx_r.push_back(sqrt(SecVtx_x*SecVtx_x+SecVtx_y*SecVtx_y));
                        tree_Hemi_SecVtx_nTrks.push_back(SecVtx_ntk);
                        tree_Hemi_SecVtx_NChi2.push_back(SecVtx_chi);
                                                // std::cout<<"SecVtx_chi: "<<SecVtx_chi<<std::endl;
                        // recX = SecVtx_x - tree_PV_x;
                        // recY = SecVtx_y - tree_PV_y;
                        // recZ = SecVtx_z - tree_PV_z;
                        // recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
                        float theta_SecVtx = TMath::ATan2(sqrt(SecrecX*SecrecX+SecrecY*SecrecY),abs(SecrecZ)) ;
                        float eta_SecVtx = -TMath::Log(tan(theta_SecVtx/2));
                        if ( SecVtx_z < 0 ) eta_SecVtx = -eta_SecVtx;
                        float SecPhi = TMath::ATan2(SecrecY,SecrecX);
	                      float dRsec = Deltar( eta_SecVtx, SecPhi, axis2_eta, axis2_phi );
	    	                if ( !isAxis1 ) dRsec = Deltar( eta_SecVtx, SecPhi, axis1_eta, axis1_phi );
	                      tree_Hemi_SecVtx_dR.push_back(dRsec);
                        tree_Hemi_SecVtx_dist.push_back( SecrecD );
                        tree_Hemi_SecVtx_track_MeanDCA_d.push_back(SecDCA_VTX_Meand);
                        tree_Hemi_SecVtx_SumtrackWeight.push_back(SecSumWeight);
                        float NewVtxMass = TMath::Max(NewVtxVector.Mag(),0.);
                        tree_Hemi_SecVtx_Mass.push_back(NewVtxMass);

                        if ( isMC_ && tree_nLLP > 0 ) {
                          if ( SecHemi == 1 ) tree_Hemi_SecLLP.push_back( iLLPrec1 );
                          if ( SecHemi == 2 ) tree_Hemi_SecLLP.push_back( iLLPrec2 );
                          if ( ((SecHemi == 1 && iLLPrec1 == 1) || 
	                              (SecHemi == 2 && iLLPrec2 == 1)) ) {
                            tree_Hemi_LLP_SecVtx_dx.push_back(SecVtx_x - LLP1_x);
                            tree_Hemi_LLP_SecVtx_dy.push_back(SecVtx_y - LLP1_y);
                            tree_Hemi_LLP_SecVtx_dz.push_back(SecVtx_z - LLP1_z);
                            tree_Hemi_LLP_SecVtx_dr.push_back(sqrt(SecVtx_x*SecVtx_x+SecVtx_y*SecVtx_y) - sqrt(LLP1_x*LLP1_x+LLP1_y*LLP1_y));
	                        }
                          if ( ((SecHemi == 1 && iLLPrec1 == 2) ||
	                              (SecHemi == 2 && iLLPrec2 == 2)) ) {
                            tree_Hemi_LLP_SecVtx_dx.push_back(SecVtx_x - LLP2_x);
                            tree_Hemi_LLP_SecVtx_dy.push_back(SecVtx_y - LLP2_y);
                            tree_Hemi_LLP_SecVtx_dz.push_back(SecVtx_z - LLP2_z);
                            tree_Hemi_LLP_SecVtx_dr.push_back(sqrt(SecVtx_x*SecVtx_x+SecVtx_y*SecVtx_y) - sqrt(LLP2_x*LLP2_x+LLP2_y*LLP2_y));
                            }

                                    dSV = (SecVtx_x - LLP1_x)*(SecVtx_x - LLP1_x) + (SecVtx_y - LLP1_y)*(SecVtx_y - LLP1_y) + (SecVtx_z - LLP1_z)*(SecVtx_z - LLP1_z);
                                    float dSV1 = TMath::Sqrt(dSV)/LLP1_dist;
                                    if ( SecVtx_chi > 0. && SecVtx_chi < 10. && 
                                              (dSV < dSVcut || dSV1 < ddcut) ) pong_Sec = 1;
                                    dSV = (SecVtx_x - LLP2_x)*(SecVtx_x - LLP2_x) + (SecVtx_y - LLP2_y)*(SecVtx_y - LLP2_y) + (SecVtx_z - LLP2_z)*(SecVtx_z - LLP2_z);
                                    float dSV2 = TMath::Sqrt(dSV)/LLP2_dist;
                                    if ( SecVtx_chi > 0. && SecVtx_chi < 10. && 
                                              (dSV < dSVcut || dSV2 < ddcut) ) pong_Sec += 2;
                            
                            if ( (pong_Merge == 3 || pong_Sec == 3) && pong_Merge*pong_Sec > 0 ) {//One of the two new vtx is close to both the gen Vtx 
                              pong_Hemi1 = true;
                              pong_Hemi2 = true;
                            }
                            else if ( pong_Merge > 0 && pong_Sec == 0 ) {//Merged Vtx is close to 
                              pong_Hemi1 = true;
                            }
                            else if ( pong_Sec > 0 && pong_Merge == 0 ) {
                              pong_Hemi2 = true;
                            }
                            else if ( pong_Merge != pong_Sec ) {
                              pong_Hemi1 = true;
                              pong_Hemi2 = true;
                            }
                            else if ( pong_Merge == pong_Sec && pong_Merge*pong_Sec > 0 ) {
                              if      ( pong_Merge == iLLPrec1 ) {
                                pong_Hemi1 = true;
                              }
                              else if ( pong_Sec == iLLPrec2 ) {
                                pong_Hemi2 = true;
                              }
                              else if ( (Vtx1_step >= 1 && Vtx1_step <= 2) || Vtx_step > 2 ) {
                                pong_Hemi1 = true;
                                }
                                else {
                                  pong_Hemi2 = true;
                                }
                            }

                          if ( isAxis1 ) {
                              tree_Hemi_SecLLP_ping.push_back( pong_Hemi1 );
                              tree_Hemi_SecLLP_ping.push_back( pong_Hemi2 );
                          }
                        else {
                          tree_Hemi_SecLLP_ping.push_back( pong_Hemi2 );
			                    tree_Hemi_SecLLP_ping.push_back( pong_Hemi1 );
	                      }
                        int pong_event = 0;
                        if	  ( pong_Hemi1 && pong_Hemi2 ) pong_event = 2;
                        else if ( pong_Hemi1 || pong_Hemi2 ) pong_event = 1;
                        tree_event_SecLLP_ping.push_back( pong_event );
                  } // endif MC
            	     
              //  } // SecVtx is Valid

        //$$
        } // endif merge vertex is valid (just to be cautious)
        //$$

      } // end of critera for merging
      else {
        tree_event_MergedVtx_Vtx_step.push_back(0);
      }

    } // endif nVertex == 2
    else {
      tree_event_MergedVtx_Vtx_step.push_back(0);
    }

    tree_event_nVtx.push_back(nVertex);

    // For the merging flag, if you have 2N vertices to merge together, there is a highchance that you get N1 (= N in theory) merged vertices thanks to the IAVF 
    // (for the signal sample at least). But you should also reconstruct N2 (N2 = N) other vertices to get 2 vertices in each event. However, we may be lacking
    // tracks for the N2 vertex => usually N2 < N1. The way this information is regitered is that the Merging2 value will be set to 0 if 
    //the vertxing does not converge but the merging can be good so Merging1 will be set to 1. (to explain the difference between N1 and N2)
    // -------------------------------------------//
  
    // "mva_Vtx_ntrk10" 
    // "mva_Vtx_ntrk20"
    // "mva_Vtx_MeanDCA
    //done after

    //------------Duplicate for each hemisphere-----------//
    // some informations for tracks in their hemisphere
    int ntrk10_vtx_hemi1 = 0., ntrk10_vtx_hemi2 = 0.;
    int ntrk20_vtx_hemi1 = 0., ntrk20_vtx_hemi2 = 0.;
    int NisjetH = 0;
    bool FromMerging = false;
    for (int counter_track = 0; counter_track < tree_nTracks; counter_track++) 
    {
      int hemi	= tree_track_Hemi[counter_track];
      double MVAval = tree_track_MVAval[counter_track];
      bool ping = false;
      
      Vtx_chi = -10.;
      float dist = -100.;
      if ( MVAval > bdtcut ) {
        if ( hemi == 1 ) {Vtx_chi = Vtx_chi1;}
        else if ( tree_Hemi_SecVtx.size()>= 1 )
          {
            if (tree_Hemi_SecVtx[0] == 1 || tree_Hemi_SecVtx[0] == 2 ) {Vtx_chi = tree_Hemi_SecVtx_NChi2[0];FromMerging = true;}
          } 
        else if ( hemi == 2 ) {Vtx_chi = Vtx_chi2;} 
        else if ( tree_Hemi_SecVtx.size()> 1 )
          {
            if (tree_Hemi_SecVtx[1] == 1 || tree_Hemi_SecVtx[1] == 2){Vtx_chi = tree_Hemi_SecVtx_NChi2[1];FromMerging = true;} 
          }
        if ( hemi == 1 && ping_Hemi1 ) ping = true;
        if ( hemi == 2 && ping_Hemi2 ) ping = true;
        if ( Vtx_chi < 10. && Vtx_chi>0 ) {
          bool isLostTrack = tree_track_lost[counter_track];
          if (isLostTrack && RemoveLostTrackFromVtxSelec ) {continue;}
          float x1 = tree_track_firstHit_x[counter_track] - tree_PV_x;
          float y1 = tree_track_firstHit_y[counter_track] - tree_PV_y;
          float z1 = tree_track_firstHit_z[counter_track] - tree_PV_z;
          float vtx_x =  tree_Hemi_Vtx_x[hemi-1] - tree_PV_x;
          float vtx_y = tree_Hemi_Vtx_y[hemi-1] - tree_PV_y;
          float vtx_z = tree_Hemi_Vtx_z[hemi-1] - tree_PV_z;
          if (FromMerging)
            {
              vtx_x = tree_Hemi_SecVtx_x[hemi-1] - tree_PV_x;
              vtx_y = tree_Hemi_SecVtx_y[hemi-1] - tree_PV_y;
              vtx_z = tree_Hemi_SecVtx_z[hemi-1] - tree_PV_z;
            }
          dist = TMath::Sqrt( (x1-vtx_x)*(x1-vtx_x) + (y1-vtx_y)*(y1-vtx_y) + (z1-vtx_z)*(z1-vtx_z) );
            int ijet = tree_track_iJet[counter_track];
            if ( ijet >= 0 ) NisjetH ++; /*!*/
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
      
    tree_track_Hemi_isjet.push_back(NisjetH);
    tree_track_Hemi_isjet.push_back(NisjetH);
    tree_Hemi_Vtx_ntrk10.push_back(ntrk10_vtx_hemi1);
    tree_Hemi_Vtx_ntrk10.push_back(ntrk10_vtx_hemi2);
    tree_Hemi_Vtx_ntrk20.push_back(ntrk20_vtx_hemi1);
    tree_Hemi_Vtx_ntrk20.push_back(ntrk20_vtx_hemi2);
    // ---------------------------------------------------//


  //////////////////////////////////////////
  ////// Vertex Selection ----------------//
  //////////////////////////////////////////

            //--------------VTX1-------------------//
// float mva_V1_r = Vtx1_r;
// float mva_V1_z = Vtx1_z;
// float mva_V1_dist = Vtx1_dist;
float mva_V1_nTrks = Vtx1_ntk;
float mva_V1_NChi2 = Vtx1_chi;
float mva_V1_step = Vtx1_step;
float mva_V1_STW = Vtx1_STW;
float mva_V1_Mass = Vtx1_Mass;
float mva_V1_HMass = H1_Mass;
float mva_V1_ntrk10 = ntrk10_vtx_hemi1;
float mva_V1_ntrk20 = ntrk20_vtx_hemi1;
float mva_V1_MeanDCA = Vtx1_MeanDCA;

// float mva_V2_r = Vtx2_r;
// float mva_V2_z = Vtx2_z;
// float mva_V2_dist = Vtx2_dist;
float mva_V2_nTrks = Vtx2_ntk;
float mva_V2_NChi2 = Vtx2_chi;
float mva_V2_step = Vtx2_step;
float mva_V2_STW = Vtx2_STW;
float mva_V2_Mass = Vtx2_Mass;
float mva_V2_HMass = H2_Mass;
float mva_V2_ntrk10 = ntrk10_vtx_hemi2;
float mva_V2_ntrk20 = ntrk20_vtx_hemi2;
float mva_V2_MeanDCA = Vtx2_MeanDCA;

if ( tree_Hemi_SecVtx.size() >= 1 ) 
  {
    // mva_V1_r = tree_Hemi_SecVtx_r[0];
    // mva_V1_z = tree_Hemi_SecVtx_z[0];
    // mva_V1_dist = tree_Hemi_SecVtx_dist[0];
    mva_V1_nTrks = tree_Hemi_SecVtx_nTrks[0];
    mva_V1_NChi2 = tree_Hemi_SecVtx_NChi2[0];
    mva_V1_step = tree_Hemi_SecVtx_step[0];
    mva_V1_STW = tree_Hemi_SecVtx_SumtrackWeight[0];
    mva_V1_Mass = tree_Hemi_SecVtx_Mass[0];
    mva_V1_MeanDCA = tree_Hemi_SecVtx_track_MeanDCA_d[0];
  }
if ( tree_Hemi_SecVtx.size() == 2 ) {
    // mva_V2_r = tree_Hemi_SecVtx_r[1];
    // mva_V2_z = tree_Hemi_SecVtx_z[1];
    // mva_V2_dist = tree_Hemi_SecVtx_dist[1];
    mva_V2_nTrks = tree_Hemi_SecVtx_nTrks[1];
    mva_V2_NChi2 = tree_Hemi_SecVtx_NChi2[1];
    mva_V2_step = tree_Hemi_SecVtx_step[1];
    mva_V2_STW = tree_Hemi_SecVtx_SumtrackWeight[1];
    mva_V2_Mass = tree_Hemi_SecVtx_Mass[1];
    mva_V2_MeanDCA = tree_Hemi_SecVtx_track_MeanDCA_d[1];
  }

if ( mva_V1_NChi2 >0 && mva_V1_NChi2 < 10 )
  {
    // mva_V_r    =  mva_V1_r;
    // mva_V_z    =  mva_V1_z;
    // mva_V_dist =  mva_V1_dist;

    mva_V_nTrks   =  mva_V1_nTrks;
    mva_V_chi	    =  mva_V1_NChi2;
    mva_V_step    =  mva_V1_step;
    mva_V_MTW     =  mva_V1_STW;
    mva_V_Mass    =  mva_V1_Mass;
    mva_H_Mass    =  mva_V1_HMass;
    mva_V_ntrk10  =  mva_V1_ntrk10;
    mva_V_ntrk20  =  mva_V1_ntrk20;
    mva_V_MeanDCA =  mva_V1_MeanDCA;
  }
    double Vtx1_bdtVal = -10;
    Vtx1_bdtVal = readerVtx->EvaluateMVA("BDTG");// values at -999 => thishappens if the Hemi_Mass and Hemi_Vtx_Mass are not definite
   // default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
    tree_Hemi_Vtx_MVAval_Loose.push_back(Vtx1_bdtVal);

  if ( (mva_V1_step == 1 || mva_V1_step == 2) &&  mva_V1_NChi2 > 0 && mva_V1_NChi2 <10 )
    {
      // mva_V_r    =   Vtx1_r;
      // mva_V_z    =    Vtx1_z;
      // mva_V_dist = Vtx1_dist;

      mva_V_nTrks   =  mva_V1_nTrks;
      mva_V_chi	    =  mva_V1_NChi2;
      mva_V_step    =  mva_V1_step;
      mva_V_MTW     =  mva_V1_STW;
      mva_V_Mass    =  mva_V1_Mass;
      mva_H_Mass    =  mva_V1_HMass;
      mva_V_ntrk10  =  mva_V1_ntrk10;
      mva_V_ntrk20  =  mva_V1_ntrk20;
      mva_V_MeanDCA =  mva_V1_MeanDCA;
    }
  double Vtx1_bdtVal_Step1 = -10; 
  Vtx1_bdtVal_Step1 = readerVtxStep1->EvaluateMVA("BDTG");
  tree_Hemi_Vtx_BDT_nTrks.push_back(mva_V_nTrks);
  tree_Hemi_Vtx_BDT_NChi2.push_back(mva_V_chi);
  tree_Hemi_Vtx_BDT_step.push_back(mva_V_step);
  tree_Hemi_Vtx_BDT_STW.push_back(mva_V_MTW);
  tree_Hemi_Vtx_BDT_Mass.push_back(mva_V_Mass);
  tree_Hemi_Vtx_BDT_HMass.push_back(mva_H_Mass);
  tree_Hemi_Vtx_BDT_ntrk10.push_back(mva_V_ntrk10);
  tree_Hemi_Vtx_BDT_ntrk20.push_back(mva_V_ntrk20);
  tree_Hemi_Vtx_BDT_MeanDCA.push_back(mva_V_MeanDCA);
  tree_Hemi_Vtx_MVAval_Tight.push_back(Vtx1_bdtVal_Step1);

	//--------------VTX2-------------------//

  if (mva_V2_NChi2>0 && mva_V2_NChi2<10 )
    {
    // mva_V_r    =  mva_V2_r;
    // mva_V_z    =  mva_V2_z;
    // mva_V_dist =  mva_V2_dist;

    mva_V_nTrks   =  mva_V2_nTrks;
    mva_V_chi	    =  mva_V2_NChi2;
    mva_V_step    =  mva_V2_step;
    mva_V_MTW     =  mva_V2_STW;
    mva_V_Mass    =  mva_V2_Mass;
    mva_H_Mass    =  mva_V2_HMass;
    mva_V_ntrk10  =  mva_V2_ntrk10;
    mva_V_ntrk20  =  mva_V2_ntrk20;
    mva_V_MeanDCA =  mva_V2_MeanDCA;
    }
  double Vtx2_bdtVal = -10;
  Vtx2_bdtVal = readerVtx->EvaluateMVA("BDTG");// values at -999
  // default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
  tree_Hemi_Vtx_MVAval_Loose.push_back(Vtx2_bdtVal);

  if ((mva_V2_step == 1 || mva_V2_step == 2) && mva_V2_NChi2>0 && mva_V2_NChi2<10)
  {
    // mva_V_r    =  mva_V2_r;
    // mva_V_z    =  mva_V2_z;
    // mva_V_dist =  mva_V2_dist;

    mva_V_nTrks   =  mva_V2_nTrks;
    mva_V_chi	    =  mva_V2_NChi2;
    mva_V_step    =  mva_V2_step;
    mva_V_MTW     =  mva_V2_STW;
    mva_V_Mass    =  mva_V2_Mass;
    mva_H_Mass    =  mva_V2_HMass;
    mva_V_ntrk10  =  mva_V2_ntrk10;
    mva_V_ntrk20  =  mva_V2_ntrk20;
    mva_V_MeanDCA =  mva_V2_MeanDCA;
  }

  double Vtx2_bdtVal_Step1 = -10; 
  Vtx2_bdtVal_Step1 = readerVtxStep1->EvaluateMVA("BDTG");
  tree_Hemi_Vtx_BDT_nTrks.push_back(mva_V_nTrks);
  tree_Hemi_Vtx_BDT_NChi2.push_back(mva_V_chi);
  tree_Hemi_Vtx_BDT_step.push_back(mva_V_step);
  tree_Hemi_Vtx_BDT_STW.push_back(mva_V_MTW);
  tree_Hemi_Vtx_BDT_Mass.push_back(mva_V_Mass);
  tree_Hemi_Vtx_BDT_HMass.push_back(mva_H_Mass);
  tree_Hemi_Vtx_BDT_ntrk10.push_back(mva_V_ntrk10);
  tree_Hemi_Vtx_BDT_ntrk20.push_back(mva_V_ntrk20);
  tree_Hemi_Vtx_BDT_MeanDCA.push_back(mva_V_MeanDCA);
  tree_Hemi_Vtx_MVAval_Tight.push_back(Vtx2_bdtVal_Step1);

//$$
  } // endif tree_njetNOmu > 0
} // endif tree_Filter
//$$

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

  //tree_LHE_Weights.clear();
    
    tree_muon_isPrompt.clear();
    tree_muon_pt.clear();
    tree_muon_SF.clear();
    tree_muon_eta.clear();
    tree_muon_phi.clear();
    tree_muon_x.clear();
    tree_muon_y.clear();
    tree_muon_z.clear();
    tree_muon_dxy.clear();
    tree_muon_dxyError.clear();
    tree_muon_dz.clear();
    tree_muon_dzError.clear();
    tree_muon_charge.clear();
    tree_muon_isLoose.clear();
    tree_muon_isMedium.clear();
    tree_muon_isTight.clear();
    tree_muon_isGlobal.clear();
    tree_muon_isoR3.clear();
    tree_muon_trigger_dimu.clear();  
    tree_muon_trigger_isomu.clear();
    tree_muon_PFIsoVeryLoose.clear();
    tree_muon_PFIsoLoose.clear();
    tree_muon_PFIsoMedium.clear();
    tree_muon_PFIsoTight.clear();
    tree_muon_TkIsoLoose.clear();
    tree_muon_TkIsoTight.clear();
    tree_muon_MiniIsoLoose.clear();
    tree_muon_MiniIsoMedium.clear();
    tree_muon_MiniIsoTight.clear();
    tree_muon_trkLayers.clear();
    tree_muon_miniIso.clear();
    tree_muon_correction.clear();
    tree_muon_gen.clear();

    tree_reco_muon_leadingpt.clear();
    tree_reco_electron_leadingpt2.clear();
    tree_reco_muon_leadingeta.clear();
    tree_reco_electron_leadingeta2.clear();
    tree_reco_muon_leadingphi.clear();
    tree_reco_electron_leadingphi2.clear();

    tree_trig_muon_leadingpt.clear();
    tree_trig_electron_leadingpt2.clear();
    tree_trig_muon_leadingeta.clear();
    tree_trig_electron_leadingeta2.clear();
    tree_trig_muon_leadingphi.clear();
    tree_trig_electron_leadingphi2.clear();

    tree_reco_lepton_leadingpt.clear();
    tree_reco_lepton_leadingpt2.clear();
    tree_reco_lepton_leadingeta.clear();
    tree_reco_lepton_leadingeta2.clear();
    tree_reco_lepton_leadingphi.clear();
    tree_reco_lepton_leadingphi2.clear();
    tree_trig_lepton_leadingpt.clear();
    tree_trig_lepton_leadingpt2.clear();
    tree_trig_lepton_leadingeta.clear();
    tree_trig_lepton_leadingeta2.clear();
    tree_trig_lepton_leadingphi.clear();
    tree_trig_lepton_leadingphi2.clear();

    tree_lepton_leadingpt.clear();
    tree_lepton_leadingpt2.clear();
    tree_lepton_leadingeta.clear();
    tree_lepton_leadingeta2.clear();
    tree_lepton_leadingphi.clear();
    tree_lepton_leadingphi2.clear();

    tree_lepton_lepton_dR.clear();
    tree_lepton_lepton_dPhi.clear();
    tree_lepton_lepton_dEta.clear();

    tree_ll_pt.clear();
    tree_ll_eta.clear();
    tree_ll_phi.clear();
    tree_ll_px.clear();
    tree_ll_py.clear();
    tree_ll_pz.clear();
    tree_ll_energy.clear();
    tree_ll_mass.clear();

    tree_electron_isPrompt.clear();
    tree_electron_pt.clear();
    tree_electron_eta.clear();
    tree_electron_phi.clear();
    tree_electron_x.clear();
    tree_electron_y.clear();
    tree_electron_z.clear();
    tree_electron_energy.clear();
    tree_electron_et.clear();
    tree_electron_ecal_trk_postcorr.clear();
    tree_electron_charge.clear();
    tree_electron_isoR4.clear();
    tree_electron_IsLoose.clear();
    tree_electron_IsMedium.clear();
    tree_electron_IsTight.clear();
    tree_electron_dxy.clear();
    tree_electron_dz.clear();
    tree_electron_gen.clear();

    tree_jet_pt.clear();
    tree_jet_eta.clear();
    tree_jet_phi.clear();
    tree_jet_px.clear();
    tree_jet_py.clear();
    tree_jet_pz.clear();
    tree_jet_E.clear();
    tree_jet_tightid_LepVeto.clear();
    tree_jet_tightid.clear();
    tree_jet_TightJetIDLepVeto.clear();
    tree_jet_TightJetID.clear();
    tree_jet_HadronFlavour.clear();
    tree_jet_pileupID.clear();
    tree_jet_btag_DeepCSV.clear();
    tree_jet_btag_DeepJet.clear();
    tree_jet_leadingpt.clear();
    tree_jet_leadingpt2.clear();
    tree_jet_leadingeta.clear();
    tree_jet_leadingeta2.clear();
    tree_jet_leadingMuon_dR.clear();
    tree_jet_leadingMuon2_dR.clear();
    tree_jet_jet_dR.clear();
    tree_jet_jet_dPhi.clear();
    tree_jet_jet_dEta.clear();
    tree_muon_jet_dRmin.clear();
    tree_muon_jet_dRmax.clear();
    tree_elemu_jet_dRmin.clear();
    tree_elemu_jet_dRmax.clear();
    tree_ele_jet_dRmin.clear();
    tree_ele_jet_dRmax.clear();
    
    tree_K0_x.clear();
    tree_K0_y.clear();
    tree_K0_z.clear();
    tree_K0_r.clear();
    tree_K0_NChi2.clear();
    tree_K0_ndf.clear();
    tree_K0_mass.clear();
    tree_K0_pt.clear();
    tree_K0_eta.clear();
    tree_K0_phi.clear();
    tree_K0_nDaughters.clear();

    tree_L0_x.clear();
    tree_L0_y.clear();
    tree_L0_z.clear();
    tree_L0_r.clear();
    tree_L0_NChi2.clear();
    tree_L0_ndf.clear();
    tree_L0_mass.clear();
    tree_L0_pt.clear();
    tree_L0_eta.clear();
    tree_L0_phi.clear();
    tree_L0_nDaughters.clear();

    tree_V0_reco_x.clear();
    tree_V0_reco_y.clear();
    tree_V0_reco_z.clear();
    tree_V0_reco_r.clear();
    tree_V0_reco_drSig.clear();
    tree_V0_reco_dzSig.clear();
    tree_V0_reco_angleXY.clear();
    tree_V0_reco_angleZ.clear();
    tree_V0_reco_NChi2.clear();
    tree_V0_reco_ndf.clear();
    tree_V0_reco_mass.clear();
    tree_V0_reco_pt.clear();
    tree_V0_reco_eta.clear();
    tree_V0_reco_phi.clear();
    tree_V0_reco_source.clear();
    tree_V0_reco_badTkHit.clear();
    tree_V0_reco_dca.clear();

    tree_SecInt_x.clear();
    tree_SecInt_y.clear();
    tree_SecInt_z.clear();
    tree_SecInt_r.clear();
    tree_SecInt_d.clear();
    tree_SecInt_drSig.clear();
    tree_SecInt_dzSig.clear();
    tree_SecInt_angleXY.clear();
    tree_SecInt_angleZ.clear();
    tree_SecInt_NChi2.clear();
    tree_SecInt_ndf.clear();
    tree_SecInt_mass.clear();
    tree_SecInt_pt.clear();
    tree_SecInt_eta.clear();
    tree_SecInt_phi.clear();
    tree_SecInt_charge.clear();
    tree_SecInt_badTkHit.clear();
    tree_SecInt_dca.clear();
    tree_SecInt_selec.clear();
    tree_SecInt_layer.clear();
    tree_SecInt_LLP.clear();
    tree_SecInt_LLP_dr.clear();
    tree_SecInt_LLP_dz.clear();
    tree_SecInt_LLP_dd.clear();
    tree_SecInt_tk1.clear();
    tree_SecInt_tk2.clear();

    tree_Yc_x.clear(); 
    tree_Yc_y.clear();
    tree_Yc_z.clear();
    tree_Yc_r.clear();
    tree_Yc_dr0.clear();
    tree_Yc_dr1.clear();
    tree_Yc_dz0.clear();
    tree_Yc_dz1.clear();
    tree_Yc_costheta.clear();
    tree_Yc_layer.clear();
    tree_Yc_NChi2.clear();
    tree_Yc_ndf.clear();
    tree_Yc_nDaughters.clear();
    tree_Yc_pt.clear();
    tree_Yc_eta.clear();
    tree_Yc_phi.clear();
    tree_Yc_mass.clear();
    tree_Yc_tracks_index.clear();
    tree_Yc_tracks_charge.clear();
    tree_Yc_tracks_pt.clear();
    tree_Yc_tracks_eta.clear();
    tree_Yc_tracks_phi.clear();
    tree_Yc_tracks_phi0.clear();

    tree_track_ipc.clear();
    tree_track_lost.clear();
    tree_track_px.clear();
    tree_track_py.clear();
    tree_track_pz.clear(); 
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
    tree_track_dzSig.clear();
//$$$$$$
//     tree_track_dzTOpu.clear();
//     tree_track_dzSigTOpu.clear();
//$$$$$$
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
    tree_track_btag.clear();
    tree_track_energy.clear();

    tree_track_Hemi.clear();
    tree_track_Hemi_dR.clear();
    tree_track_Hemi_dRmax.clear();
    tree_track_Hemi_mva_NChi2.clear();
    tree_track_Hemi_ping.clear();
    tree_track_Hemi_dFirstVtx.clear();
    tree_track_Hemi_LLP.clear();
//$$$$$$
//     tree_track_Hemi_d0.clear();
//     tree_track_Hemi_d0Sig.clear();
//     tree_track_HemiOp_d0.clear();
//     tree_track_HemiOp_d0Sig.clear();
//$$$$$$
    tree_track_Hemi_isjet.clear();
        
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

//$$$$
    tree_V0_track_isFromV0.clear();
    tree_V0_track_isFromSI.clear();
    tree_V0_track_lost.clear();
    tree_V0_track_pt.clear();
    tree_V0_track_eta.clear();
    tree_V0_track_phi.clear();
    tree_V0_track_charge.clear();
    tree_V0_track_NChi2.clear();
    tree_V0_track_dxy.clear();
    tree_V0_track_drSig.clear();
    tree_V0_track_dz.clear();
    tree_V0_track_dzSig.clear();
    tree_V0_track_nHit.clear();
    tree_V0_track_nHitPixel.clear();
    tree_V0_track_firstHit.clear();
    tree_V0_track_firstHit_x.clear();
    tree_V0_track_firstHit_y.clear();
    tree_V0_track_firstHit_z.clear();
    tree_V0_track_iJet.clear();
    tree_V0_track_ntrk10.clear();
    tree_V0_track_ntrk20.clear();
    tree_V0_track_ntrk30.clear();
    tree_V0_track_ntrk40.clear();
    tree_V0_track_Hemi.clear();
    tree_V0_track_Hemi_dR.clear();
    tree_V0_track_Hemi_dRmax.clear();
//$$$$

    tree_genParticle_pt.clear();
    tree_genParticle_eta.clear();
    tree_genParticle_phi.clear();
    tree_genParticle_charge.clear();
    tree_genParticle_pdgId.clear();
    tree_genParticle_mass.clear();
    tree_genParticle_x.clear();
    tree_genParticle_y.clear();
    tree_genParticle_z.clear();
    tree_genParticle_px.clear();
    tree_genParticle_py.clear();
    tree_genParticle_pz.clear();
    tree_genParticle_energy.clear();
    tree_genParticle_isPromptFinalState.clear();
    tree_genParticle_ct.clear();
    tree_genParticle_ct0.clear();
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
    tree_genAxis_dPhineuneu.clear();
    tree_genAxis_dEtaneuneu.clear();
    tree_GenAxes_Mass.clear();
    tree_GenAxes_CombinedHemiLeptonMass.clear();
    tree_GenAxis_Neu_dRmin.clear();
    tree_GenAxis_Neu_dRmax.clear();
    tree_GenAxis_RecoAxis_dRmin.clear();
    tree_GenAxis_RecoAxis_dRmax.clear();

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
    tree_genFromB_dr.clear();
    tree_genFromB_dz.clear();
    tree_genFromB_dd.clear();

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
    tree_LLP_r.clear();
    tree_LLP_dist.clear();
    tree_LLP_nTrks.clear();
    tree_LLP_Vtx_NChi2.clear();
    tree_LLP_Vtx_nTrks.clear();
    tree_LLP_Vtx_dx.clear();
    tree_LLP_Vtx_dy.clear();
    tree_LLP_Vtx_dz.clear();
    tree_LLP_Vtx_dist.clear();
    tree_LLP_Vtx_dd.clear();
    tree_LLP12_dR.clear();
    tree_LLP12_deta.clear();
    tree_LLP12_dphi.clear();
    tree_LLP_Mass.clear();

    tree_Hemi.clear();
    tree_Hemi_njet.clear();
    tree_Hemi_njet_nomu.clear();
    tree_Hemi_pt.clear();
    tree_Hemi_eta.clear();
    tree_Hemi_phi.clear();
    tree_Hemi_nTrks.clear();
    tree_Hemi_nTrks_sig.clear();
    tree_Hemi_nTrks_bad.clear();
    tree_Hemi_mass.clear();
    tree_HemiMu_mass.clear();
    tree_HemiMu_pt.clear();
    tree_HemiMu_dR.clear();
    tree_HemiMuOp_mass.clear();
    tree_HemiMuOp_pt.clear();
    tree_HemiMuOp_dR.clear();
    tree_Hemi_LooseBTag_axes.clear();
    tree_Hemi_MediumBTag_axes.clear();
    tree_Hemi_TightBTag_axes.clear();
    tree_Hemi_dR12.clear();

    tree_Hemi_LLP.clear();
    tree_Hemi_LLP_pt.clear();
    tree_Hemi_LLP_eta.clear();
    tree_Hemi_LLP_phi.clear();
    tree_Hemi_LLP_dist.clear();
    tree_Hemi_LLP_x.clear();
    tree_Hemi_LLP_y.clear();
    tree_Hemi_LLP_z.clear();
    tree_Hemi_LLP_dR.clear();
    tree_Hemi_LLP_mother.clear();
    tree_Hemi_LLP_Vtx_dx.clear();
    tree_Hemi_LLP_Vtx_dy.clear();
    tree_Hemi_LLP_Vtx_dz.clear();
    tree_Hemi_LLP_Vtx_dr.clear();
    tree_Hemi_LLP_muOK_dR.clear();
    tree_Hemi_LLP_muOK_pt.clear();
    tree_Hemi_LLP_muOK_mass.clear();
    tree_Hemi_LLP_muNO_dR.clear();
    tree_Hemi_LLP_muNO_pt.clear();
    tree_Hemi_LLP_muNO_mass.clear();
    tree_Hemi_LLP_dR12.clear();
    tree_Hemi_LLP_ping.clear();
    tree_event_LLP_ping.clear();
    
    tree_Hemi_Vtx_step.clear();
    tree_Hemi_Vtx_isTight.clear();
    tree_Hemi_Vtx_NChi2.clear();
    tree_Hemi_Vtx_nTrks.clear();
    tree_Hemi_Vtx_nTrks_sig.clear();
    tree_Hemi_Vtx_nTrks_bad.clear();
    tree_Hemi_Vtx_x.clear();
    tree_Hemi_Vtx_y.clear();
    tree_Hemi_Vtx_z.clear();
    tree_Hemi_Vtx_r.clear();
    tree_Hemi_Vtx_dR.clear();
    tree_Hemi_Vtx_xError.clear();
    tree_Hemi_Vtx_yError.clear();
    tree_Hemi_Vtx_zError.clear();
    tree_Hemi_Vtx_trackWeight.clear();
    tree_Hemi_Vtx_SumtrackWeight.clear();
    tree_Hemi_Vtx_track_MeanDCA_d.clear();
    tree_Hemi_Vtx_BTag.clear();
    tree_Hemi_Vtx_Mass.clear();
    tree_Hemi_Vtx_dist.clear();
    tree_Hemi_Vtx_ntrk10.clear();
    tree_Hemi_Vtx_ntrk20.clear();
    tree_event_nVtx.clear();
    tree_event_Vtx_Vtx_dr.clear();
    tree_event_Vtx_Vtx_dz.clear();
    tree_event_Vtx_Vtx_dd.clear();
    tree_event_Vtx_Vtx_reldd.clear();
    tree_event_Vtx_Vtx_dR.clear();
    tree_event_Vtx_Vtx_step.clear();

    tree_Hemi_SecLLP.clear();
    tree_Hemi_LLP_SecVtx_dx.clear();
    tree_Hemi_LLP_SecVtx_dy.clear();
    tree_Hemi_LLP_SecVtx_dz.clear();
    tree_Hemi_LLP_SecVtx_dr.clear();
    tree_Hemi_SecLLP_ping.clear();
    tree_event_SecLLP_ping.clear();

    // tree_Hemi_SecVtx_track_DCA_x.clear();
    // tree_Hemi_SecVtx_track_DCA_y.clear();
    // tree_Hemi_SecVtx_track_DCA_z.clear();
    // tree_Hemi_SecVtx_track_DCA_r.clear();
    // tree_Hemi_SecVtx_track_DCA_d.clear();
    tree_Hemi_SecVtx.clear();
    tree_Hemi_SecVtx_step.clear();
    tree_Hemi_SecVtx_x.clear();
    tree_Hemi_SecVtx_y.clear();
    tree_Hemi_SecVtx_z.clear();
    tree_Hemi_SecVtx_r.clear();
    tree_Hemi_SecVtx_dR.clear();
    tree_Hemi_SecVtx_nTrks.clear();
    tree_Hemi_SecVtx_NChi2.clear();
    tree_Hemi_SecVtx_dist.clear();
    tree_Hemi_SecVtx_track_MeanDCA_d.clear();
    tree_Hemi_SecVtx_SumtrackWeight.clear();
    tree_Hemi_SecVtx_trackWeight.clear();
    tree_Hemi_SecVtx_Mass.clear();
    tree_event_MergedVtx_Vtx_dr.clear();
    tree_event_MergedVtx_Vtx_dz.clear();
    tree_event_MergedVtx_Vtx_dd.clear();
    tree_event_MergedVtx_Vtx_reldd.clear();
    tree_event_MergedVtx_Vtx_dR.clear();
    tree_event_MergedVtx_Vtx_step.clear();
        
    tree_Hemi_Vtx_BDT_nTrks.clear();
    tree_Hemi_Vtx_BDT_NChi2.clear();
    tree_Hemi_Vtx_BDT_step.clear();
    tree_Hemi_Vtx_BDT_STW.clear();
    tree_Hemi_Vtx_BDT_Mass.clear();
    tree_Hemi_Vtx_BDT_HMass.clear();
    tree_Hemi_Vtx_BDT_ntrk10.clear();
    tree_Hemi_Vtx_BDT_ntrk20.clear();
    tree_Hemi_Vtx_BDT_MeanDCA.clear();
    tree_Hemi_Vtx_MVAval_Loose.clear();
    tree_Hemi_Vtx_MVAval_Tight.clear();
}
