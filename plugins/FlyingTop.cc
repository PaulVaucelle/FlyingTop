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
#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>

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

    std::string weightFile_;
    std::string weightFileEVTS_;
    std::string weightFileVtx_;
    std::string weightFileVtxStep1_;

    const edm::EDGetTokenT<GenEventInfoProduct>         genEventInfoToken_;
    const edm::EDGetTokenT<LHEEventProduct>             LHEEventProductToken_;
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
    
    // edm::ParameterSet kvfPSet;
    //trig
    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
    //trig
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> K0Token_;
    edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> LambdaToken_;
    // edm::EDGetTokenT<reco::VertexCompositePtrCandidateCollection> SVToken_;
    // edm::EDGetTokenT<PileupSummaryInfo> puToken_ ;
  
    edm::EDGetTokenT<reco::ConversionCollection> PhotonToken_;
    const edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
    edm::EDGetTokenT<reco::CaloClusterCollection>clusterToken_;
    edm::EDGetTokenT<reco::CaloClusterCollection> showerToken_;
    edm::EDGetTokenT<reco::SuperClusterCollection>superclusterToken_;
    // edm::EDGetTokenT<pat::PackedTriggerPrescales> PrescaleToken_;

    int runNumber, eventNumber, lumiBlock;
    int  tree_NbrOfZCand;
    bool tree_Filter;
    
    int  tree_GoodMu1, tree_GoodMu2;

    int  tree_nTracks, tree_nLostTracks, tree_TRACK_SIZE; 
    int  tree_nFromC = 0, tree_nFromB = 0; 
    int nEvent;
    
    float LLP1_pt, LLP1_eta, LLP1_phi, LLP2_pt, LLP2_eta, LLP2_phi;
    float LLP1_x, LLP1_y, LLP1_z, LLP2_x, LLP2_y, LLP2_z;
    float LLP1_dist, LLP2_dist;
    int   LLP1_nTrks = 0, LLP2_nTrks = 0;

//$$
//The BDT variables are declared here to reduce computation time
  //EVTS level BDT to select signal events

     float  mva_Evts_MET_et;
     float  mva_Evts_nTrks;
     float  mva_Evts_muon1_pt;
     float  mva_Evts_muon2_pt;
     float  mva_Evts_jet1_pt;
     float  mva_Evts_jet2_pt;
     float  mva_Evts_muon12_dR;
     float  mva_Evts_muon12_dPhi;
     float  mva_Evts_muon12_dEta;
     float  mva_Evts_jet12_dR;
     float  mva_Evts_jet12_dPhi;
     float  mva_Evts_jet12_dEta;
     float  mva_Evts_muon_jet_dRmin;
     float  mva_Evts_muon_jet_dRmax;
     float  mva_Evts_nVtx;
     float  mva_HT;
     float  mva_Evts_ST;
     float  mva_Evts_njets;
     float  mva_Evts_nmuon;
     float  mva_Evts_MediumAxes;
     float  mva_Evts_LooseAxes;
     float  mva_Evts_TightAxes;
    
    TMVA::Reader *readerEvts = new TMVA::Reader("!Color:Silent");

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
    float dzTopu;
    float dzSigTopu;
    float TibHit;
    float TobHit;
    float PixBarHit;
    float TecHit;
    float isLost;
    //Track level BDT for displaced track selection
    TMVA::Reader *reader = new TMVA::Reader( "!Color:Silent" );

    //Vtx levelt BDT to select signal vertices
    float mva_V_nTrks;
    float mva_V_chi;
    float mva_V_step;
    float mva_V_r;
    float mva_V_z;
    float mva_V_MTW;
    float mva_V_Mass;
    float mva_H_Mass;
    float mva_V_dist;
    float mva_V_MeanDCA;
    float mva_V_ntrk10;
    float mva_V_ntrk20;

    TMVA::Reader *readerVtx = new TMVA::Reader( "!Color:Silent" );
    TMVA::Reader *readerVtxStep1 = new TMVA::Reader( "!Color:Silent" );
    
    int index[10000];
    double MVAval[10000];

    //  ---------------------------------------------------------------- //
    //  ------- Booleans to activate/desactivate part of the code ------ //
    //  ---------------------------------------------------------------- //
    bool showlog=false;//not used atm
    bool isMC = true;//not used yet

    bool MuonChannel        = true; //=> If false, look at the electron channel
    bool NewCovMat          = true;//Keep True : Allow for Covariance Matrix correction due to the MiniAOD dataformat apporixmation
    
          //Vetos to find vertices from different secondary interactions
    bool DetailedMap        = true;// Detailed map of the CMS tracker to apply a veto on the tracks of the vertices that belong to this map
        // Vetos applied on tracks of the vertices belonging to V0Candidates, Photon conversions and Secindary Interactions
    bool ActivateV0Veto     = true;
    bool ActivateYcVeto     = true;// This Veto is not doing anything on RunIISummer20UL18 TTvar and on MC signal
    bool ActivateSecIntVeto = true;
        // Activate steps of the vertexing workflow (better to keep everything true for the development, since we may want to keeep the tight WP => true false true false false)
    bool ActivateStep1      = true;//TIghtwp, standard AVF
    bool IterAVF            = true; // Activate IAVF step of the vertexing 
    bool ActivateStep2      = true;// Tight WP IAVF => Works like this (ActivateStep2 || IterAVF)
    bool ActivateStep3      = true;// LooseWP STep3 => standard AVF
    bool ActivateStep4      = true;// LooseWP STep4  IAVF => (ActivateStep4 || IterAVF)

    // -- B Tagging related information
    // WorkingPoints : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
    // Updated Page by BTV : https://btv-wiki.docs.cern.ch/ScaleFactors/UL2018/
    bool ActivateBtag       = true ;
    bool ActivateBtagLog    = false;
    float LooseWP  = 0.0490;
    float MediumWP = 0.2783;
    float TightWP  = 0.7100;

    //-Track Selection-------/
    bool IncludeLostTrack= true;// Keep as true
    bool RequestHighPurity = false ; // does not rlly matter
    //-END of Track Selection-------/

    //------Vtx Selection ---//
    bool RemoveLostTrackFromVtxSelec = true;
    //  ---------------------------------------------------------------- //
 

    //  ---------------------------------------------------------------- //
    //  -------------------- track preselection cuts ------------------- //
    //  ---------------------------------------------------------------- //
    float pt_Cut = 1;    // default 1. 
    float NChi2_Cut = 5; // default 5. 
    float drSig_Cut = 5; // default 5. 
   //  ---------------------------------------------------------------- //

    float tree_Evt_weight;

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

    std::vector<float> tree_Evts_MVAval;
    std::vector<int>   tree_allPV_i;
    std::vector<float> tree_allPV_x;
    std::vector<float> tree_allPV_y;
    std::vector<float> tree_allPV_z;
    std::vector<float> tree_allPV_ex;
    std::vector<float> tree_allPV_ey;
    std::vector<float> tree_allPV_ez;
    std::vector<float> tree_allPV_NChi2;
    std::vector<int>   tree_allPV_ndf;

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
    std::vector<int>       tree_SecInt_ntrk10;
    std::vector<int>       tree_SecInt_ntrk20;
    std::vector<int>       tree_SecInt_LLP;
    std::vector<float>     tree_SecInt_LLP_dr;
    std::vector<float>     tree_SecInt_LLP_dz;
    std::vector<float>     tree_SecInt_LLP_dd;


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
    std::vector<float> tree_jet_btag_DeepCSV;
    std::vector<float> tree_jet_btag_DeepJet;
    std::vector<float> tree_jet_leadingpt;
    std::vector<float> tree_jet_leadingpt2;
    std::vector<float> tree_jet_leadingMuon_dR;
    std::vector<float> tree_jet_leadingMuon2_dR;
    std::vector<float> tree_jet_jet_dR;
    std::vector<float> tree_jet_jet_dPhi;
    std::vector<float> tree_jet_jet_dEta;
    std::vector<float> tree_muon_jet_dRmin;
    std::vector<float> tree_muon_jet_dRmax;
    float tree_HT;
    
    //--------------------------------
    // electrons infos -------
    //--------------------------------
    int                tree_electron_nEle;
    std::vector<float> tree_electron_pt;
    std::vector<float> tree_electron_eta;
    std::vector<float> tree_electron_phi;
    std::vector<float> tree_electron_x;
    std::vector<float> tree_electron_y;
    std::vector<float> tree_electron_z;
    std::vector<float> tree_electron_px;
    std::vector<float> tree_electron_py;
    std::vector<float> tree_electron_pz;
    std::vector<float> tree_electron_energy;
    std::vector< int > tree_electron_charge;
    std::vector<float> tree_electron_isoR4;
    std::vector<bool>  tree_electron_IsLoose;
    std::vector<bool>  tree_electron_IsMedium;
    std::vector<bool>  tree_electron_IsTight;
    std::vector<bool>  tree_electron_trigger_Ele;
    std::vector<bool>  tree_electron_trigger_diEle;
    std::vector<float> tree_electron_dxy;
    std::vector<float> tree_electron_dz;
    
    //--------------------------------
    // muons infos -------
    //--------------------------------
    float tree_Mmumu;
    std::vector<float> tree_muon_pt;
    std::vector<float> tree_ST;
    std::vector<float> tree_muon_eta;
    std::vector<float> tree_muon_phi;
    std::vector<float> tree_muon_x;
    std::vector<float> tree_muon_y;
    std::vector<float> tree_muon_z;
    std::vector<float> tree_muon_px;
    std::vector<float> tree_muon_py;
    std::vector<float> tree_muon_pz;
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
    std::vector<bool>  tree_muon_trigger_dimu;  
    std::vector<bool>  tree_muon_trigger_isomu;
    std::vector<bool>  tree_muon_PFIsoLoose;
    std::vector<bool>  tree_muon_PFIsoMedium;
    std::vector<bool>  tree_muon_PFIsoTight;
    std::vector<bool>  tree_muon_TkIsoLoose;
    std::vector<bool>  tree_muon_TkIsoTight;
    std::vector<float> tree_muon_nmu;
    std::vector<float> tree_muon_leadingpt;
    std::vector<float> tree_muon_leadingpt2;
    std::vector<float> tree_muon_muon_dR;
    std::vector<float> tree_muon_muon_dPhi;
    std::vector<float> tree_muon_muon_dEta;

    //-----------------------
    // per track
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
    std::vector<float>    tree_track_dzTOpu;  // with respect to clostest PU
    std::vector<float>    tree_track_dzSigTOpu;
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
    std::vector<float>    tree_track_btag;
    std::vector<float>    tree_track_energy;

    std::vector< int >    tree_track_Hemi;
    std::vector< float >  tree_track_Hemi_dR;//dRmin
    std::vector< float >  tree_track_Hemi_dRmax;
    std::vector< float >  tree_track_Hemi_mva_NChi2;
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

    //$$$$
    int tree_smu_mass = 0;
    int tree_neu_mass = 0;
    int tree_neu_ctau = 0;
//$$$$
    
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
    // std::vector< float > tree_genFromb_x;
    // std::vector< float > tree_genFromb_y;
    // std::vector< float > tree_genFromb_z;
    // std::vector< int >  tree_genFromb_pdgid;
    // std::vector< float > tree_genFromb_LLP_dV;
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
    std::vector< float > tree_Hemi_Mass;
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
    std::vector< float > tree_Hemi_Vtx_xError;
    std::vector< float > tree_Hemi_Vtx_yError;
    std::vector< float > tree_Hemi_Vtx_zError;
    std::vector< float > tree_Hemi_Vtx_eta;
    std::vector< float > tree_Hemi_Vtx_Vtx_dr;
    std::vector< float > tree_Hemi_Vtx_Vtx_dz;
    std::vector< float > tree_Hemi_Vtx_Vtx_dd;
    std::vector< float > tree_Hemi_Vtx_BTag;
    std::vector< int >   tree_Hemi_Vtx_nVtx;
    std::vector< float > tree_Hemi_Vtx_trackWeight;
    std::vector< float > tree_Hemi_Vtx_MeantrackWeight;//Vertx selection variable for the BDT
    std::vector< float > tree_Hemi_Vtx_Mass;
    std::vector< float > tree_Hemi_Vtx_MVAval;
    std::vector< float > tree_Hemi_Vtx_MVAval_Step1;//TIght WP

    std::vector< float > tree_Hemi_Vtx_track_DCA_x;
    std::vector< float > tree_Hemi_Vtx_track_DCA_y;
    std::vector< float > tree_Hemi_Vtx_track_DCA_z;
    std::vector< float > tree_Hemi_Vtx_track_DCA_r;
    std::vector< float > tree_Hemi_Vtx_track_DCA_d;
    std::vector< float > tree_Hemi_Vtx_track_MeanDCA_d;//Veertex selection BDT

    std::vector< float > tree_Hemi_Vtx_TVtx_dx;
    std::vector< float > tree_Hemi_Vtx_TVtx_dy;
    std::vector< float > tree_Hemi_Vtx_TVtx_dz;
    std::vector< float > tree_Hemi_Vtx_TVtx_NChi2;
      
    std::vector< float > tree_Hemi_Vtx_dist;
    std::vector< float > tree_Hemi_Vtx_dx;
    std::vector< float > tree_Hemi_Vtx_dy;
    std::vector< float > tree_Hemi_Vtx_dz;
    std::vector< float > tree_Hemi_Vtx_dr;
    std::vector< float > tree_Hemi_Vtx_dd;
    std::vector< float > tree_Hemi_dR12;
    std::vector< float > tree_Hemi_LLP_dR12;
    std::vector< float > tree_Hemi_Vtx_ddbad;
    std::vector< int >   tree_Hemi_Vtx_ntrk10;//Vertex selection variables
    std::vector< int >   tree_Hemi_Vtx_ntrk20;
    std::vector< int >   tree_track_Hemi_isjet;
    std::vector< float > tree_Hemi_Vtx_ddToBkg;
    std::vector< bool >  tree_Hemi_LLP_ping;//set as true 
    std::vector< int >   tree_event_LLP_ping;

    std::vector< int >   tree_Hemi_LooseBTag_axes;
    std::vector< int >   tree_Hemi_MediumBTag_axes;
    std::vector< int >   tree_Hemi_TightBTag_axes;


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

//------------------------------------
// - Propagators init. ---------------
//------------------------------------

    PropaHitPattern* PHP = new PropaHitPattern();
    PropaHitPattern* NI = new PropaHitPattern();
    PropaHitPattern* posPHP = new PropaHitPattern();
    PropaHitPattern* negPHP = new PropaHitPattern();
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

    weightFile_( iConfig.getUntrackedParameter<std::string>("weightFileMVA") ),
    weightFileEVTS_ (iConfig.getUntrackedParameter<std::string>("weightFileMVA_EVTS")), 
    weightFileVtx_( iConfig.getUntrackedParameter<std::string>("weightFileMVA_VTX") ),
    weightFileVtxStep1_( iConfig.getUntrackedParameter<std::string>("weightFileMVA_VTX_step1") ),
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
    LambdaToken_(   consumes<reco::VertexCompositePtrCandidateCollection>(            iConfig.getParameter<edm::InputTag>("Lambda")))
    // ,puToken_(      consumes<PileupSummaryInfo>(                                iConfig.getParameter<edm::InputTag>("pileup")))
    ,PhotonToken_(  consumes<reco::ConversionCollection>(edm::InputTag(std::string("reducedEgamma"),std::string("reducedConversions")))) 
    ,beamSpotToken_(     consumes<reco::BeamSpot>(               iConfig.getUntrackedParameter<edm::InputTag>("beamSpot"))),
    clusterToken_ (consumes<reco::CaloClusterCollection>(edm::InputTag(std::string("reducedEgamma"),std::string("reducedEBEEClusters"),std::string("RECO")))), //enlever les clusters
    showerToken_ (consumes<reco::CaloClusterCollection>(edm::InputTag(std::string("reducedEgamma"),std::string("educedESClusters"),std::string("RECO")))),//enlever les clusters
    superclusterToken_ (consumes<reco::SuperClusterCollection>(edm::InputTag(std::string("reducedEgamma"),std::string("reducedSuperClusters"),std::string("RECO"))))//enlever les clusters
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
    smalltree->Branch("tree_Evt_weight",&tree_Evt_weight, "tree_Evt_weight/I");
    // primary vertex info
    smalltree->Branch("tree_nPV",      &tree_nPV);
    smalltree->Branch("tree_PV_x",     &tree_PV_x);
    smalltree->Branch("tree_PV_y",     &tree_PV_y);
    smalltree->Branch("tree_PV_z",     &tree_PV_z);    
    smalltree->Branch("tree_PV_ez",    &tree_PV_ez);    
    smalltree->Branch("tree_PV_NChi2", &tree_PV_NChi2);    
    smalltree->Branch("tree_PV_ndf",   &tree_PV_ndf);

    smalltree->Branch("tree_Evts_MVAval", &tree_Evts_MVAval);
    smalltree->Branch("tree_allPV_i",     &tree_allPV_i);
    smalltree->Branch("tree_allPV_x",     &tree_allPV_x);
    smalltree->Branch("tree_allPV_y",     &tree_allPV_y);
    smalltree->Branch("tree_allPV_z",     &tree_allPV_z);
    smalltree->Branch("tree_allPV_ex",    &tree_allPV_ex);
    smalltree->Branch("tree_allPV_ey",    &tree_allPV_ey);
    smalltree->Branch("tree_allPV_ez",    &tree_allPV_ez);
    smalltree->Branch("tree_allPV_NChi2", &tree_allPV_NChi2);
    smalltree->Branch("tree_allPV_ndf",   &tree_allPV_ndf);

    //Beamspot
    smalltree->Branch("tree_bs_PosX", &tree_bs_PosX) ;
    smalltree->Branch("tree_bs_PosY", &tree_bs_PosY) ;
    smalltree->Branch("tree_bs_PosZ", &tree_bs_PosZ) ;

    smalltree->Branch("tree_NbrOfZCand",  &tree_NbrOfZCand,  "tree_NbrOfZCand/I");
    smalltree->Branch("tree_Filter", &tree_Filter);
    //$$$$
    smalltree->Branch("tree_GoodMu1", &tree_GoodMu1);
    smalltree->Branch("tree_GoodMu2", &tree_GoodMu2);
//$$$$

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
    smalltree->Branch("tree_SecInt_ntrk10",     &tree_SecInt_ntrk10);
    smalltree->Branch("tree_SecInt_ntrk20",     &tree_SecInt_ntrk20);
    smalltree->Branch("tree_SecInt_LLP",        &tree_SecInt_LLP);
    smalltree->Branch("tree_SecInt_LLP_dr",     &tree_SecInt_LLP_dr);
    smalltree->Branch("tree_SecInt_LLP_dz",     &tree_SecInt_LLP_dz);
    smalltree->Branch("tree_SecInt_LLP_dd",     &tree_SecInt_LLP_dd);

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
    smalltree->Branch("tree_jet_btag_DeepCSV",&tree_jet_btag_DeepCSV);
    smalltree->Branch("tree_jet_btag_DeepJet",&tree_jet_btag_DeepJet);
    smalltree->Branch("tree_jet_leadingpt",&tree_jet_leadingpt);
    smalltree->Branch("tree_jet_leadingpt2",&tree_jet_leadingpt2);
    smalltree->Branch("tree_jet_leadingMuon_dR",&tree_jet_leadingMuon_dR);
    smalltree->Branch("tree_jet_leadingMuon2_dR",&tree_jet_leadingMuon2_dR);
    smalltree->Branch("tree_jet_jet_dR",&tree_jet_jet_dR);
    smalltree->Branch("tree_jet_jet_dPhi",&tree_jet_jet_dPhi);
    smalltree->Branch("tree_jet_jet_dEta",&tree_jet_jet_dEta);
    smalltree->Branch("tree_muon_jet_dRmin",&tree_muon_jet_dRmin);
    smalltree->Branch("tree_muon_jet_dRmax",&tree_muon_jet_dRmax);
    smalltree->Branch("tree_HT"  ,          &tree_HT);
    
    // electrons info
    smalltree->Branch("tree_electron_nEle",   &tree_electron_nEle);
    smalltree->Branch("tree_electron_pt"  ,   &tree_electron_pt);
    smalltree->Branch("tree_electron_eta" ,   &tree_electron_eta);
    smalltree->Branch("tree_electron_phi" ,   &tree_electron_phi);
    smalltree->Branch("tree_electron_x"  ,    &tree_electron_x);
    smalltree->Branch("tree_electron_y" ,     &tree_electron_y);
    smalltree->Branch("tree_electron_z" ,     &tree_electron_z);
    smalltree->Branch("tree_electron_px",     &tree_electron_px);
    smalltree->Branch("tree_electron_py",     &tree_electron_py);
    smalltree->Branch("tree_electron_pz",     &tree_electron_pz);
    smalltree->Branch("tree_electron_energy", &tree_electron_energy);
    smalltree->Branch("tree_electron_charge", &tree_electron_charge);
    smalltree->Branch("tree_electron_isoR4",  &tree_electron_isoR4);
    smalltree->Branch("tree_electron_IsLoose",&tree_electron_IsLoose);
    smalltree->Branch("tree_electron_IsMedium",&tree_electron_IsMedium);
    smalltree->Branch("tree_electron_IsTight",&tree_electron_IsTight);
    smalltree->Branch("tree_electron_trigger_Ele",&tree_electron_trigger_Ele);
    smalltree->Branch("tree_electron_trigger_diEle",&tree_electron_trigger_diEle);
    smalltree->Branch("tree_electron_dxy",    &tree_electron_dxy);
    smalltree->Branch("tree_electron_dz",     &tree_electron_dz);

    smalltree->Branch("tree_ST",&tree_ST);
    // muons info
    smalltree->Branch("tree_Mmumu"  ,         &tree_Mmumu);
    smalltree->Branch("tree_muon_pt"  ,       &tree_muon_pt);
    smalltree->Branch("tree_muon_eta" ,       &tree_muon_eta);
    smalltree->Branch("tree_muon_phi" ,       &tree_muon_phi);
    smalltree->Branch("tree_muon_x"  ,        &tree_muon_x);
    smalltree->Branch("tree_muon_y" ,         &tree_muon_y);
    smalltree->Branch("tree_muon_z" ,         &tree_muon_z);
    smalltree->Branch("tree_muon_px"  ,       &tree_muon_px);
    smalltree->Branch("tree_muon_py" ,        &tree_muon_py);
    smalltree->Branch("tree_muon_pz" ,        &tree_muon_pz);   
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
    smalltree->Branch("tree_muon_trigger_dimu",&tree_muon_trigger_dimu);  
    smalltree->Branch("tree_muon_trigger_isomu",&tree_muon_trigger_isomu);
    smalltree->Branch("tree_muon_PFIsoLoose",&tree_muon_PFIsoLoose);
    smalltree->Branch("tree_muon_PFIsoMedium",&tree_muon_PFIsoMedium);
    smalltree->Branch("tree_muon_PFIsoTight",&tree_muon_PFIsoTight);
    smalltree->Branch("tree_muon_TkIsoLoose", &tree_muon_TkIsoLoose);
    smalltree->Branch("tree_muon_TkIsoTight", &tree_muon_TkIsoTight);
    smalltree->Branch("tree_muon_nmu",&tree_muon_nmu);
    smalltree->Branch("tree_muon_leadingpt",&tree_muon_leadingpt);
    smalltree->Branch("tree_muon_leadingpt2",&tree_muon_leadingpt2);
    smalltree->Branch("tree_muon_muon_dR",&tree_muon_muon_dR);
    smalltree->Branch("tree_muon_muon_dPhi",&tree_muon_muon_dPhi);
    smalltree->Branch("tree_muon_muon_dEta",&tree_muon_muon_dEta);

    // track
    smalltree->Branch("tree_TRACK_SIZE", &tree_TRACK_SIZE, "tree_TRACK_SIZE/I");
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
    smalltree->Branch("tree_track_dzTOpu",       &tree_track_dzTOpu);
    smalltree->Branch("tree_track_dzSigTOpu",    &tree_track_dzSigTOpu  );
    smalltree->Branch("tree_track_algo",         &tree_track_algo);
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

    //$$$$
    smalltree->Branch("tree_smu_mass" ,  &tree_smu_mass);
    smalltree->Branch("tree_neu_mass" ,  &tree_neu_mass);
    smalltree->Branch("tree_neu_ctau" ,  &tree_neu_ctau);
//$$$$
    
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
    // smalltree->Branch("tree_genFromb_x",&tree_genFromb_x);
    // smalltree->Branch("tree_genFromb_y",&tree_genFromb_y);
    // smalltree->Branch("tree_genFromb_z",&tree_genFromb_z);
    // smalltree->Branch("tree_genFromb_pdgid",&tree_genFromb_pdgid);
    // smalltree->Branch("tree_genFromb_LLP_dV",&tree_genFromb_LLP_dV);
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
    smalltree->Branch("tree_Hemi_Mass",  &tree_Hemi_Mass);
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
    smalltree->Branch("tree_Hemi_Vtx_xError",&tree_Hemi_Vtx_xError);
    smalltree->Branch("tree_Hemi_Vtx_yError",&tree_Hemi_Vtx_yError);
    smalltree->Branch("tree_Hemi_Vtx_zError",&tree_Hemi_Vtx_zError);
    smalltree->Branch("tree_Hemi_Vtx_eta",   &tree_Hemi_Vtx_eta);
    smalltree->Branch("tree_Hemi_Vtx_Vtx_dr",&tree_Hemi_Vtx_Vtx_dr);
    smalltree->Branch("tree_Hemi_Vtx_Vtx_dz",&tree_Hemi_Vtx_Vtx_dz);
    smalltree->Branch("tree_Hemi_Vtx_Vtx_dd",&tree_Hemi_Vtx_Vtx_dd);
    smalltree->Branch("tree_Hemi_Vtx_BTag",  &tree_Hemi_Vtx_BTag);
    smalltree->Branch("tree_Hemi_Vtx_nVtx",  &tree_Hemi_Vtx_nVtx);
    smalltree->Branch("tree_Hemi_Vtx_trackWeight",&tree_Hemi_Vtx_trackWeight);
    smalltree->Branch("tree_Hemi_Vtx_MeantrackWeight",&tree_Hemi_Vtx_MeantrackWeight);
    smalltree->Branch("tree_Hemi_Vtx_track_DCA_x",&tree_Hemi_Vtx_track_DCA_x);
    smalltree->Branch("tree_Hemi_Vtx_track_DCA_y",&tree_Hemi_Vtx_track_DCA_y);
    smalltree->Branch("tree_Hemi_Vtx_track_DCA_z",&tree_Hemi_Vtx_track_DCA_z);
    smalltree->Branch("tree_Hemi_Vtx_track_DCA_r",&tree_Hemi_Vtx_track_DCA_r);
    smalltree->Branch("tree_Hemi_Vtx_track_DCA_d",&tree_Hemi_Vtx_track_DCA_d);
    smalltree->Branch("tree_Hemi_Vtx_track_MeanDCA_d",&tree_Hemi_Vtx_track_MeanDCA_d);
    smalltree->Branch("tree_Hemi_Vtx_Mass", &tree_Hemi_Vtx_Mass);
    smalltree->Branch("tree_Hemi_Vtx_MVAval", &tree_Hemi_Vtx_MVAval);
    smalltree->Branch("tree_Hemi_Vtx_MVAval_Step1",&tree_Hemi_Vtx_MVAval_Step1);
    smalltree->Branch("tree_Hemi_Vtx_TVtx_dx",&tree_Hemi_Vtx_TVtx_dx);
    smalltree->Branch("tree_Hemi_Vtx_TVtx_dy",&tree_Hemi_Vtx_TVtx_dy);
    smalltree->Branch("tree_Hemi_Vtx_TVtx_dz",&tree_Hemi_Vtx_TVtx_dz);
    smalltree->Branch("tree_Hemi_Vtx_TVtx_NChi2",&tree_Hemi_Vtx_TVtx_NChi2); 
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
    smalltree->Branch("tree_track_Hemi_isjet",&tree_track_Hemi_isjet);
    smalltree->Branch("tree_Hemi_Vtx_ddToBkg",&tree_Hemi_Vtx_ddToBkg);
    smalltree->Branch("tree_Hemi_LLP_ping",  &tree_Hemi_LLP_ping);
    smalltree->Branch("tree_event_LLP_ping", &tree_event_LLP_ping);
    smalltree->Branch("tree_Hemi_LooseBTag_axes",&tree_Hemi_LooseBTag_axes);
    smalltree->Branch("tree_Hemi_MediumBTag_axes",&tree_Hemi_MediumBTag_axes);
    smalltree->Branch("tree_Hemi_TightBTag_axes",&tree_Hemi_TightBTag_axes);

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
smalltree->Branch("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v);
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
smalltree->Branch("HLT_PFHT180_v",&HLT_PFHT180_v);
smalltree->Branch("HLT_PFHT250_v",&HLT_PFHT250_v);
smalltree->Branch("HLT_PFHT370_v",&HLT_PFHT370_v);
smalltree->Branch("HLT_PFHT430_v",&HLT_PFHT430_v);
smalltree->Branch("HLT_PFHT510_v",&HLT_PFHT510_v);
smalltree->Branch("HLT_PFHT590_v",&HLT_PFHT590_v);
smalltree->Branch("HLT_PFHT680_v",&HLT_PFHT680_v);
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


//----------------------------------------
// - BDT Input Variables -----------------
//----------------------------------------
//---------------EVTS--------------------

readerEvts->AddVariable( "mva_Evts_MET_et",             &mva_Evts_MET_et);
readerEvts->AddVariable( "mva_Evts_nTrks",              &mva_Evts_nTrks);
readerEvts->AddVariable( "mva_Evts_muon1_pt",           &mva_Evts_muon1_pt);
readerEvts->AddVariable( "mva_Evts_muon2_pt",           &mva_Evts_muon2_pt);
readerEvts->AddVariable( "mva_Evts_jet1_pt",            &mva_Evts_jet1_pt);
readerEvts->AddVariable( "mva_Evts_jet2_pt",            &mva_Evts_jet2_pt);
readerEvts->AddVariable( "mva_Evts_muon12_dR",          &mva_Evts_muon12_dR);
readerEvts->AddVariable( "mva_Evts_muon12_dPhi",        &mva_Evts_muon12_dPhi);
readerEvts->AddVariable( "mva_Evts_muon12_dEta",        &mva_Evts_muon12_dEta);
readerEvts->AddVariable( "mva_Evts_jet12_dR",           &mva_Evts_jet12_dR);
readerEvts->AddVariable( "mva_Evts_jet12_dPhi",         &mva_Evts_jet12_dPhi);
readerEvts->AddVariable( "mva_Evts_jet12_dEta",         &mva_Evts_jet12_dEta);
readerEvts->AddVariable( "mva_Evts_muon_jet_dRmin",     &mva_Evts_muon_jet_dRmin);
readerEvts->AddVariable( "mva_Evts_muon_jet_dRmax",     &mva_Evts_muon_jet_dRmax);
// readerEvts->AddVariable( "mva_Evts_nVtx",               &mva_Evts_nVtx);
readerEvts->AddVariable( "mva_HT",                      &mva_HT);
readerEvts->AddVariable( "mva_Evts_ST",                 &mva_Evts_ST);
readerEvts->AddVariable( "mva_Evts_njets",              &mva_Evts_njets);
readerEvts->AddVariable( "mva_Evts_nmuon",              &mva_Evts_nmuon);
readerEvts->AddVariable( "mva_Evts_MediumAxes",         &mva_Evts_MediumAxes);
readerEvts->AddVariable( "mva_Evts_LooseAxes",          &mva_Evts_LooseAxes);
readerEvts->AddVariable( "mva_Evts_TightAxes",          &mva_Evts_TightAxes);
readerEvts->BookMVA( "BDTG", weightFileEVTS_ ); // root 6.14/09, care compatiblity of versions for tmva

//-------------Tracks---------------------
//$$
    //add the variables from my BDT (Paul)
    // reader->AddVariable( "mva_track_firstHit_x", &firsthit_X); /*!*/
    // reader->AddVariable( "mva_track_firstHit_y", &firsthit_Y); /*!*/
    // reader->AddVariable( "mva_track_firstHit_z", &firsthit_Z); /*!*/
    reader->AddVariable( "mva_track_dxy",     &dxy);
    reader->AddVariable( "mva_track_dz",      &dz);
    reader->AddVariable( "mva_track_pt",      &pt );
    reader->AddVariable( "mva_track_eta",     &eta );
    reader->AddVariable( "mva_track_nchi2",   &NChi );
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
    // reader->AddVariable(" mva_track_lost", &isLost);

    // reader->AddVariable(" mva_track_dxyError", &dxyError); /*!*/
    // reader->AddVariable(" mva_track_dzTOpu", &dzTopu);//added on 24/03/2023 : if using bdts generated before this date, =>crash
    // reader->AddVariable(" mva_track_dzSigTOpu", &dzSigTopu);//added on 24/03/2023
    // reader->AddVariable(" mva_ValTIBHit", &TibHit);
    // reader->AddVariable(" mva_ValTOBHit", &TobHit);
    // reader->AddVariable(" mva_ValPixBarHit", &PixBarHit);
    // reader->AddVariable(" mva_nValTECHHit", &TecHit);

    reader->BookMVA( "BDTG", weightFile_ ); // root 6.14/09, care compatiblity of versions for tmva

    //-------------Vertex---------------------
    readerVtx->AddVariable( "mva_Vtx_nTrks",     &mva_V_nTrks);
    readerVtx->AddVariable( "mva_Vtx_NChi2",     &mva_V_chi);
    readerVtx->AddVariable( "mva_Vtx_step",     &mva_V_step);
    // readerVtx->AddVariable( "mva_Vtx_r",     &mva_V_r);
    // readerVtx->AddVariable( "mva_Vtx_z",     &mva_V_z);
    readerVtx->AddVariable( "mva_Vtx_MTW",     &mva_V_MTW);
    readerVtx->AddVariable( "mva_Vtx_Mass",     &mva_V_Mass);
    readerVtx->AddVariable("mva_Hemi_Mass",&mva_H_Mass);
    // readerVtx->AddVariable( "mva_Vtx_dist",     &mva_V_dist);
    readerVtx->AddVariable( "mva_Vtx_ntrk10",     &mva_V_ntrk10);
    readerVtx->AddVariable( "mva_Vtx_ntrk20",     &mva_V_ntrk20);
    readerVtx->AddVariable(" mva_Vtx_MeanDCA",&mva_V_MeanDCA);
    readerVtx->BookMVA( "BDTG", weightFileVtx_ ); // root 6.14/09, care compatiblity of versions for tmva
    //----------------Vertex Step1----------//
    readerVtxStep1->AddVariable( "mva_Vtx_nTrks",     &mva_V_nTrks);
    readerVtxStep1->AddVariable( "mva_Vtx_NChi2",     &mva_V_chi);
    readerVtxStep1->AddVariable( "mva_Vtx_step",     &mva_V_step);
    // readerVtxStep1->AddVariable( "mva_Vtx_r",     &mva_V_r);
    // readerVtxStep1->AddVariable( "mva_Vtx_z",     &mva_V_z);
    readerVtxStep1->AddVariable( "mva_Vtx_MTW",     &mva_V_MTW);
    readerVtxStep1->AddVariable( "mva_Vtx_Mass",     &mva_V_Mass);
    readerVtxStep1->AddVariable("mva_Hemi_Mass",&mva_H_Mass);
    // readerVtxStep1->AddVariable( "mva_Vtx_dist",     &mva_V_dist);
    readerVtxStep1->AddVariable( "mva_Vtx_ntrk10",     &mva_V_ntrk10);
    readerVtxStep1->AddVariable( "mva_Vtx_ntrk20",     &mva_V_ntrk20);
    readerVtxStep1->AddVariable(" mva_Vtx_MeanDCA",&mva_V_MeanDCA);
    readerVtxStep1->BookMVA("BDTG",weightFileVtxStep1_);
//$$
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

  runNumber   = iEvent.id().run();
  eventNumber = iEvent.id().event();
  lumiBlock   = iEvent.luminosityBlock();
  tree_Evt_weight = 1;
// std::cout<<"event : "<<eventNumber<<std::endl;
  using namespace edm;
  using namespace reco;
  using namespace pat;
  
  bool runOnData_ = false;

  //LHE/gen infos

  edm::Handle<GenEventInfoProduct> genEventInfo;
  if(!runOnData_) iEvent.getByToken(genEventInfoToken_,genEventInfo);
  
  edm::Handle<LHEEventProduct> lheEventProduct;
  if(!runOnData_) iEvent.getByToken(LHEEventProductToken_,lheEventProduct);
  
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

  //---calocluster
  edm::Handle<reco::CaloClusterCollection> CaloClusters;
  iEvent.getByToken(clusterToken_,CaloClusters);
 
  edm::Handle<reco::CaloClusterCollection> ESClusters;
  iEvent.getByToken(showerToken_,ESClusters);

 edm::Handle<reco::SuperClusterCollection> SClusters;
 iEvent.getByToken(superclusterToken_,SClusters);

 //----------------------//

 iSetup.get<IdealMagneticFieldRecord>().get(bField);
 const MagneticField* theMagneticField = bField.product();

 edm::ESHandle<TrackerGeometry> trackerGeomHandle;
 iSetup.get<TrackerDigiGeometryRecord>().get( trackerGeomHandle );

// const TrackerGeometry* trackerGeom = trackerGeomHandle.product();

 //--------------------------------------//
 //               Trigger                //
 //--------------------------------------//
  std::vector<std::string> TriggerCheck;
  // std::vector<std::string> VTrigger;
  // std::vector<bool> PassTrigger;

  // ######################
  const edm::Handle<edm::TriggerResults> triggerH = iEvent.getHandle(triggerResultsToken_);
  const auto triggerNames = iEvent.triggerNames(*triggerH);
  std::string TName;
  // if (eventNumber==1)
  // {
  //   for (unsigned int i = 0; i < triggerH->size(); i++) 
  //     {
  //       if (TName.substr(0,6) ! HLT_Mu" && TName.substr(0,7) ! HLT_Ele" && TName.substr(0,12) ! HLT_DoubleMu" && TName.substr(0,13) ! HLT_DoubleEle" &&TName.substr(0,11) ! HLT_Dimuon0" &&  TName.substr(0,9) ! HLT_PFMET" &&TName.substr(0,6) ! HLT_HT" && TName.substr(0,7) ! HLT_AK4" && TName.substr(0,9) ! HLT_PFJet")
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
//           // if(TName.substr(0,6) = HLT_Mu"){tree_Trigger_Muon.push_back(TName);}//+ dilepton channel emu
//           // if(TName.substr(0,7) = HLT_Ele"){tree_Trigger_Ele.push_back(TName);}
//           // if(TName.substr(0,12) = HLT_DoubleMu"){tree_Trigger_DoubleMu.push_back(TName);}
//           // if(TName.substr(0,13) = HLT_DoubleEle"){tree_Trigger_DoubleEle.push_back(TName);}
//           // if(TName.substr(0,11) = HLT_Dimuon0"){tree_Trigger_Dimuon0.push_back(TName);}
//           // if(TName.substr(0,9) = HLT_PFMET"){tree_Trigger_PFMET.push_back(TName);}
//           // if(TName.substr(0,6) = HLT_HT"){tree_Trigger_HT.push_back(TName);}
//           // if(TName.substr(0,7) = HLT_AK4"){tree_Trigger_AK4.push_back(TName);}
//           // if(TName.substr(0,9) = HLT_PFJet"){tree_Trigger_PFJet.push_back(TName);}
//           // if(TName.substr(0,16) = HLT_DoublePFJets"){tree_Trigger_DoublePFJets.push_back(TName);}
//           // if(TName.substr(0,11) = HLT_DiPFJet"){tree_Trigger_DiPFJet.push_back(TName);}
//           // if(TName.substr(0,13) = HLT_QuadPFJet"){tree_Trigger_QuadPFJet.push_back(TName);}
//           // if(TName.substr(0,10) = HLT_BTagMu"){tree_Trigger_BTagMu.push_back(TName);}
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
/* 76%*/if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") &&  triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v = true;} else if (strstr(TName.c_str(),"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v") && !triggerH->accept(i)){HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v = false;};
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
if (strstr(TName.c_str(),"HLT_PFHT180_v") && triggerH->accept(i)){HLT_PFHT180_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT180_v") &&!triggerH->accept(i)){HLT_PFHT180_v = false;};
if (strstr(TName.c_str(),"HLT_PFHT250_v") && triggerH->accept(i)){HLT_PFHT250_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT250_v") &&!triggerH->accept(i)){HLT_PFHT250_v = false;};
if (strstr(TName.c_str(),"HLT_PFHT370_v") && triggerH->accept(i)){HLT_PFHT370_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT370_v") &&!triggerH->accept(i)){HLT_PFHT370_v = false;};
if (strstr(TName.c_str(),"HLT_PFHT430_v") && triggerH->accept(i)){HLT_PFHT430_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT430_v") &&!triggerH->accept(i)){HLT_PFHT430_v = false;};
if (strstr(TName.c_str(),"HLT_PFHT510_v") && triggerH->accept(i)){HLT_PFHT510_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT510_v") &&!triggerH->accept(i)){HLT_PFHT510_v = false;};
if (strstr(TName.c_str(),"HLT_PFHT590_v") && triggerH->accept(i)){HLT_PFHT590_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT590_v") &&!triggerH->accept(i)){HLT_PFHT590_v = false;};
if (strstr(TName.c_str(),"HLT_PFHT680_v") && triggerH->accept(i)){HLT_PFHT680_v = true;} else if (strstr(TName.c_str(),"HLT_PFHT680_v") &&!triggerH->accept(i)){HLT_PFHT680_v = false;};
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

  //////////////////////////////////
  //////////////////////////////////
  //////////    LHE    /////////////
  //////////////////////////////////
  //////////////////////////////////

//----------GEN------------//
  std::vector<double> evtWeights = genEventInfo->weights();// only weight[0] is stored
  double theWeight = genEventInfo->weight();

    const gen::PdfInfo *PDF = genEventInfo->pdf();
    // std::cout<<"scalePDF : "<<PDF->scalePDF<<std::endl;
    int id1 = PDF->id.first ;// [-4,-3,-2,-1,1,2,3,4]
    int id2 = PDF->id.second;
    // std::cout<<"id1 and id2 : "<<id1<<"//"<<id2<<std::endl;

   double x1 = PDF->x.first;
   double x2 = PDF->x.second;
  //  std::cout<<"x1 and x2 : "<<x1<<"//"<<x2<<std::endl;

   double xPDF1 = PDF->xPDF.first;//==0
   double xPDF2 = PDF->xPDF.second;//==0
  //  std::cout<<"xPDF1 and xPDF2 : "<<xPDF1<<"//"<<xPDF2<<std::endl;

  unsigned int ProcID = genEventInfo->signalProcessID();//9999
	// double qscale = genEventInfo->qScale();//sameasPDFscale
  double alphaqcd = genEventInfo->alphaQCD();
// std::cout<<"ProcID and qscale and alphaqcd : "<<ProcID<<"//"<<" alphaqcd : "<<alphaqcd<<std::endl;

// for (unsigned int k = 0 ; k<evtWeights.size() ; k++)
//   {
//     std::cout<<"evtWeights k : "<<evtWeights[k]<<std::endl;
//   }

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
  
  for (unsigned int i = 0; i< primaryVertex->size() ; ++i) {
    tree_allPV_i.push_back(i);
    tree_allPV_x.push_back((*primaryVertex)[i].x());
    tree_allPV_y.push_back((*primaryVertex)[i].y());
    tree_allPV_z.push_back((*primaryVertex)[i].z());
    tree_allPV_ex.push_back((*primaryVertex)[i].xError());
    tree_allPV_ey.push_back((*primaryVertex)[i].yError());
    tree_allPV_ez.push_back((*primaryVertex)[i].zError());
    tree_allPV_NChi2.push_back((*primaryVertex)[i].normalizedChi2());
    tree_allPV_ndf.push_back((*primaryVertex)[i].ndof());
  }

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
if ( DetailedMap ) {VtxLayerNI = NI->VertexBelongsToTracker(Yr, Yz);}
else {
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

//$$
    if ( !(VtxLayerNI != 0 && Ymass > 0. && Ymass < 1. && Ynchi2 < 10.) ) continue;
//$$
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
//$$$$
  float smumass = 0., neumass = 0.;
//$$$$

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
    float GenPtMin = 20;   // (GeV) minimum jet pt is optimum
    float GenEtaMax = 10; // no cut on eta is optimum
    int Genjetidx = 0; // : May be in the loop/ not sure it changes anything

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
    //$$$$
	  smumass = mom->mass();
	  neumass = Gen_m;
//$$$$
	}
	if ( nLLP == 2 ) {
	  LLP2_pt  = Gen_pt;
	  LLP2_eta = Gen_eta;
	  LLP2_phi = Gen_phi;
    //$$$$
	  smumass += mom->mass();
	  neumass += Gen_m;
	  tree_smu_mass = (smumass/2. + 0.5);
	  tree_neu_mass = (neumass/2. + 0.5);
//$$$$
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
            tree_genFromB_dr.push_back((motherInPrunedCollection->vx()-gen2->vx())*(motherInPrunedCollection->vx()-gen2->vx()) 
                                        +(motherInPrunedCollection->vy()-gen2->vy())*(motherInPrunedCollection->vy()-gen2->vy()) );
            tree_genFromB_dz.push_back((motherInPrunedCollection->vz()-gen2->vz())*(motherInPrunedCollection->vz()-gen2->vz()) );
            tree_genFromB_dd.push_back(sqrt(     (motherInPrunedCollection->vx()-gen2->vx())*(motherInPrunedCollection->vx()-gen2->vx()) 
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
      tree_genParticle_energy.push_back(	    genIt.energy());
      tree_genParticle_isPromptFinalState.push_back(   genIt.isPromptFinalState());
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
        float bg = vgen.E() / mom->mass();
        ct0 = ct / bg;
      }
      tree_genParticle_ct.push_back(ct);
      tree_genParticle_ct0.push_back(ct0);

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


    //////////////////////
    // gen jets
    //////////////////////

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

      int nGenMuon = 0;
      for (unsigned int k = 0 ; k < tree_genParticle_pdgId.size() ; k++) //loop over genMuons
        {
          if (abs(tree_genParticle_pdgId[k])!=13) continue;//pruned collection may be should check also with the packed colelction
          //feels like there are sometimes two gen muons that are ony one?? close to having the same pt eta and phi (deltaQuantity  ~ 0.001) => continue
          if ( !tree_genParticle_isPromptFinalState[k] ) continue;
          if ( tree_genParticle_pt[k] < 25) continue;
          nGenMuon++;
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
    }// end loop over gen jets

//     ///////////////////////////////
//     // Invariant Mass of GenAxes 
//     ///////////////////////////////
    
    //The computation of the invariant mass does not look good when using only the tracks
    //So, we try to do it with the jets and the possible letpon coming out of the neutralino
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
        if (abs(tree_genParticle_pdgId[i])!=13) continue;
        if ( tree_genParticle_isPromptFinalState[i] ) continue;//||  !tree_muon_isTight[mu]
        // Avoid prompt leptons
        float GendRAxis1muon = Deltar( Genvaxis1.Eta(), Genvaxis1.Phi(), tree_genParticle_eta[i], tree_genParticle_phi[i] );
        float GendRAxis2muon = Deltar( Genvaxis2.Eta(), Genvaxis2.Phi(), tree_genParticle_eta[i], tree_genParticle_phi[i] );
        if (GendRAxis1muon<1.5  )
          {
            temp_px1 += tree_genParticle_px[i];
            temp_py1 += tree_genParticle_py[i];
            temp_pz1 += tree_genParticle_pz[i];
            temp_e1  += tree_genParticle_energy[i];
          }
        if(GendRAxis2muon<1.5)
          {
            temp_px2 += tree_genParticle_px[i];
            temp_py2 += tree_genParticle_py[i];
            temp_pz2 += tree_genParticle_pz[i];
            temp_e2  += tree_genParticle_energy[i];
          }
      }

    for(unsigned int i = 0 ; i <tree_genParticle_pdgId.size() ; i++)
      {
        if (abs(tree_genParticle_pdgId[i])!=11) continue;
        if(!tree_genParticle_isPromptFinalState[i]) continue; //tree_electron_IsTight
        // Avoid prompt leptons
        float GendRAxis1electron = Deltar( Genvaxis1.Eta(), Genvaxis1.Phi(), tree_genParticle_eta[i], tree_genParticle_phi[i] );
        float GendRAxis2electron = Deltar( Genvaxis2.Eta(), Genvaxis2.Phi(), tree_genParticle_eta[i], tree_genParticle_phi[i] );
        if (GendRAxis1electron<1.5  )
          {
            temp_px1 += tree_genParticle_px[i];
            temp_py1 += tree_genParticle_py[i];
            temp_pz1 += tree_genParticle_pz[i];
            temp_e1  += tree_genParticle_energy[i];
          }
        if(GendRAxis2electron<1.5)
          {
            temp_px2 += tree_genParticle_px[i];
            temp_py2 += tree_genParticle_py[i];
            temp_pz2 += tree_genParticle_pz[i];
            temp_e2  += tree_genParticle_energy[i];
          }
      }
    TLorentzVector TLorentzGenAxis1(temp_px1,temp_py1,temp_pz1,temp_e1);
    TLorentzVector TLorentzGenAxis2(temp_px2,temp_py2,temp_pz2,temp_e2);
    tree_GenAxes_Mass.push_back(sqrt(TLorentzGenAxis1.Mag2()));
    // std::cout<<"Mass GenAxis 1 : "<<sqrt(TLorentzGenAxis1.Mag2())<<std::endl;
    tree_GenAxes_Mass.push_back(sqrt(TLorentzGenAxis2.Mag2()));
    // std::cout<<"Mass GenAxis 2 : "<<sqrt(TLorentzGenAxis2.Mag2()) <<std::endl;

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
  ///////////   Muons   ////////////
  //////////////////////////////////
  //////////////////////////////////
  
  int nmu = 0;
  float LT = 0;
   // float isoR3=0;
    // bool triggered(const char* pathName) const { return triggerObjectMatchByPath(pathName, true, true) != nullptr; }
    int genMuons = 0;
    float muonpt = 0;

  for (const pat::Muon &mu : *muons)
  {
  if ( mu.pt() < 3. ) continue;
  if ( abs(mu.eta()) > 2.5 ) continue;  // muon acceptance
    tree_muon_pt.push_back(       mu.pt());
    // std::cout<<" mu.pt()"<< mu.pt()<<std::endl;
    LT += mu.pt();
    tree_muon_eta.push_back(      mu.eta());
    tree_muon_phi.push_back(      mu.phi());
    tree_muon_x.push_back(        mu.vx());
    tree_muon_y.push_back(        mu.vy());
    tree_muon_z.push_back(        mu.vz());
    tree_muon_px.push_back(       mu.px());
    tree_muon_py.push_back(       mu.py());
    tree_muon_pz.push_back(       mu.pz());
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
    tree_muon_trigger_dimu.push_back(mu.triggered("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*"));  
    tree_muon_trigger_isomu.push_back(mu.triggered("HLT_IsoMu24_v*"));
    // reco::TrackRef MuonInnerRef = mu.innerTrack();

    // https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/DataFormats/MuonReco/interface/Muon.h
    tree_muon_PFIsoLoose.push_back( mu.passed(reco::Muon::Selector::PFIsoLoose));
    tree_muon_PFIsoMedium.push_back(mu.passed(reco::Muon::Selector::PFIsoMedium));
    tree_muon_PFIsoTight.push_back( mu.passed(reco::Muon::Selector::PFIsoTight));
    tree_muon_TkIsoLoose.push_back( mu.passed(reco::Muon::Selector::TkIsoLoose));
    tree_muon_TkIsoTight.push_back( mu.passed(reco::Muon::Selector::TkIsoTight));
        if (muonpt < mu.pt() && nmu >=1 ) { std::cout<<"MAYDAYDYAYDAYDYAYDA1"<<std::endl;}
    muonpt = mu.pt();
    nmu++;

    // Matching to gen muons
    if ( !runOnData_ )
      {
    bool IsMatched = false;
    for (unsigned int k = 0 ; k < tree_genParticle_pdgId.size() ; k++)
      {
        if (abs(tree_genParticle_pdgId[k])!=13) continue;
        if (tree_genParticle_charge[k] != mu.charge()) continue;
        //feels like there are sometimes two gen muons that are ony one?? close to having the same pt eta and phi (deltaQuantity  ~ 0.001) => continue
        float deta = tree_genParticle_eta[k]-mu.eta();
        float dphi = tree_genParticle_phi[k]-mu.phi();
        float dpt  = (tree_genParticle_pt[k]-mu.pt())/tree_genParticle_pt[k];
        if(abs(dpt)<0.1 && abs(dphi)<0.1 && abs(deta)<0.1 )
          {
            IsMatched = true;
          }
          genMuons++;
          if ( IsMatched ) break;
      }
  }
  }

    if ( nmu >= 1 )
      {
        tree_muon_leadingpt.push_back(tree_muon_pt[0]);  
        mva_Evts_muon1_pt = tree_muon_pt[0];
      }
    if ( nmu > 1 )
      {
        tree_muon_leadingpt2.push_back(tree_muon_pt[1]);
        tree_muon_muon_dR.push_back(Deltar( tree_muon_eta[0], tree_muon_phi[0], tree_muon_eta[1], tree_muon_phi[1]));
        tree_muon_muon_dPhi.push_back(abs(Deltaphi(tree_muon_phi[0],tree_muon_phi[1])));
        tree_muon_muon_dEta.push_back(abs( tree_muon_eta[0]-tree_muon_eta[1]));
        mva_Evts_muon2_pt         = tree_muon_pt[1];
        mva_Evts_muon12_dR        = Deltar( tree_muon_eta[0], tree_muon_phi[0], tree_muon_eta[1], tree_muon_phi[1]);
        mva_Evts_muon12_dPhi      = abs(Deltaphi(tree_muon_phi[0],tree_muon_phi[1]));
        mva_Evts_muon12_dEta      = abs( tree_muon_eta[0]-tree_muon_eta[1]);
      }
    if ( nmu == 0){tree_muon_leadingpt.push_back(0);mva_Evts_muon1_pt=0;}
    if ( nmu <= 1)
    {
      tree_muon_leadingpt2.push_back(0);
      tree_muon_muon_dR.push_back(0);
      tree_muon_muon_dPhi.push_back(0);
      tree_muon_muon_dEta.push_back(0);
      mva_Evts_muon2_pt         = 0;
      mva_Evts_muon12_dR        = 0;
      mva_Evts_muon12_dPhi      = 0;
      mva_Evts_muon12_dEta      = 0;
    }
  float ST = 0;
  ST+=LT;
  tree_ST.push_back(ST);
  tree_muon_nmu.push_back(nmu);


  //////////////////////////////////
  //////////////////////////////////
  ///////////	Jets   /////////////
  //////////////////////////////////
  ///////////////////////////////////
  
  tree_njet = 0;
  float HT_val = 0;
  float jet_pt_min = 20.;
  int indjet = -1;
  float jetpt=0;
  float jet1_pt =0;
  float jet2_pt = 0;
  for (const pat::Jet &jet : *jets) 
  {
    indjet++;
    if ( jet.pt() < jet_pt_min ) continue;
    if ( indjet==0) {jet1_pt = jet.pt();}
    if ( indjet==1) {jet2_pt = jet.pt();}
    if (jetpt<jet.pt() && tree_njet>=1){ std::cout<<"MAYDAYDYAYDAYDYAYDA2 "<<tree_njet<<" pt1 and pt+1 : "<<jetpt<<" and "<<jet.pt()<<std::endl;}
    jetpt=jet.pt();
    tree_jet_E.push_back(jet.energy());
    tree_jet_pt.push_back(jet.pt());
    tree_jet_eta.push_back(jet.eta());
    tree_jet_phi.push_back(jet.phi());
  // std::cout<<"jet.pt()"<<jet.pt()<<std::endl;
    // btag infos :
    float DeepCSVb = jet.bDiscriminator("pfDeepCSVJetTags:probb");
    float DeepCSVbb = jet.bDiscriminator("pfDeepCSVJetTags:probbb");
    float DeepFlavourb = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
    float DeepFlavourbb = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
    float DeepFlavourblep = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
    float DeepCSV = DeepCSVb + DeepCSVbb;
    float DeepJet = -10.;
    if ( DeepFlavourb > -5 ) DeepJet = DeepFlavourb + DeepFlavourbb + DeepFlavourblep;
    tree_jet_btag_DeepCSV.push_back(DeepCSV);
    tree_jet_btag_DeepJet.push_back(DeepJet);
    // WorkingPoints : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X

    
    if ( abs(jet.eta()) < 2.4 ) {HT_val += jet.pt();} // used in HT filter !
    if (DeepJet>MediumWP && nmu>=1)
      {
        tree_jet_leadingMuon_dR.push_back(Deltar( jet.eta(), jet.phi(), tree_muon_eta[0], tree_muon_phi[0] ));
        if (nmu>1)
          {
            tree_jet_leadingMuon2_dR.push_back(Deltar( jet.eta(), jet.phi(), tree_muon_eta[1], tree_muon_phi[1] ));
          }
        else
          {
            tree_jet_leadingMuon2_dR.push_back(0);
          }
      }
    else 
    {
      tree_jet_leadingMuon_dR.push_back(0);
    }

    bool MuonJetMatching = false ; 
    //matching with jets
    for (unsigned int mu = 0 ; mu < tree_muon_pt.size() ; mu++)
      {
        float dRmujet = Deltar( jet.eta(), jet.phi(), tree_muon_eta[mu], tree_muon_phi[mu] );
        if (dRmujet <0.4)
          {
            MuonJetMatching = true;
          }
      }
    if (MuonJetMatching){continue;}   
    tree_njet++;
  }

  if (tree_njet==0)
    {
        tree_jet_jet_dR.push_back(0);
        tree_jet_jet_dPhi.push_back(0);
        tree_jet_jet_dEta.push_back(0);
        tree_jet_leadingpt.push_back(0);
        tree_jet_leadingpt2.push_back(0);
        mva_Evts_jet12_dR = 0;
        mva_Evts_jet12_dPhi = 0;
        mva_Evts_jet12_dEta  = 0;
        mva_Evts_jet1_pt = 0;
        mva_Evts_jet2_pt = 0;
    }
  if (tree_njet==1)
    {

        tree_jet_jet_dR.push_back(0);
        tree_jet_jet_dPhi.push_back(0);
        tree_jet_jet_dEta.push_back(0);
        tree_jet_leadingpt.push_back(jet1_pt);
        tree_jet_leadingpt2.push_back(0);
        mva_Evts_jet12_dR = 0;
        mva_Evts_jet12_dPhi = 0;
        mva_Evts_jet12_dEta  = 0;
        mva_Evts_jet1_pt = jet1_pt;
        mva_Evts_jet2_pt = 0;
    }
  
  if(tree_njet>=2)
      {
        tree_jet_jet_dR.push_back(Deltar(tree_jet_eta[0],tree_jet_phi[0],tree_jet_eta[1],tree_jet_phi[1]));
        tree_jet_jet_dPhi.push_back(abs(Deltaphi(tree_jet_phi[0],tree_jet_phi[1])));
        tree_jet_jet_dEta.push_back(abs(tree_jet_eta[0]-tree_jet_eta[1]));
        tree_jet_leadingpt2.push_back(tree_jet_pt[1]);
        mva_Evts_jet12_dR = Deltar(tree_jet_eta[0],tree_jet_phi[0],tree_jet_eta[1],tree_jet_phi[1]);
        mva_Evts_jet12_dPhi = abs(Deltaphi(tree_jet_phi[0],tree_jet_phi[1]));
        mva_Evts_jet12_dEta  = abs(tree_jet_eta[0]-tree_jet_eta[1]);
        mva_Evts_jet1_pt = jet1_pt;
        mva_Evts_jet2_pt = jet2_pt;    
      }

  tree_HT = HT_val;
  // ST+=HT_val;


  //////////////////////////////////
  //////////////////////////////////
  ////////   Electrons   ///////////
  //////////////////////////////////
  //////////////////////////////////
  //One could consider the emu channel, or even ee channel

  // float Eiso=0;  
  int nEl = 0;
  float elpt = 0;
  for (const pat::Electron &el: *electrons)
  {
  if ( el.pt() < 5. ) continue;
  if (abs(el.eta())>2.5 || (abs(el.eta())>1.442 && abs(el.eta())<1.556)) continue;//cluster eta
  //

    if (elpt<el.pt()&& nEl>=1){ std::cout<<"MAYDAYDYAYDAYDYAYDA"<<std::endl;}
    tree_electron_pt.push_back(     el.pt());
    tree_electron_eta.push_back(    el.eta());//cluster eta
    tree_electron_phi.push_back(    el.phi());
    tree_electron_x.push_back(      el.vx());
    tree_electron_y.push_back(      el.vy());
    tree_electron_z.push_back(      el.vz());
    tree_electron_px.push_back(      el.px());
    tree_electron_py.push_back(      el.py());
    tree_electron_pz.push_back(      el.pz());
    tree_electron_energy.push_back( el.energy());
    tree_electron_charge.push_back(el.charge());
    tree_electron_isoR4.push_back(el.trackIso());//returns the value of the summed track pt in a cone of deltaR<0.4
    tree_electron_dxy.push_back(el.gsfTrack()->dxy(PV.position()));
    tree_electron_dz.push_back(el.gsfTrack()->dz(PV.position()));

  // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_94X_and_later
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCategoryBasedElectronID
  // HLT Ele23Ele12 noDZ v*
  // HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
    tree_electron_IsLoose.push_back(el.electronID("cutBasedElectronID-Fall17-94X-V2-loose"));//for 2018
    tree_electron_IsMedium.push_back(el.electronID("cutBasedElectronID-Fall17-94X-V2-medium"));//for 2018
    tree_electron_IsTight.push_back(el.electronID("cutBasedElectronID-Fall17-94X-V2-tight"));//for 2018
    nEl++;

        elpt=el.pt();

    // tree_electron_trigger_Ele.push_back(el.triggered("HLT_Ele27_WPTight_Gsf_v*"));
    // tree_electron_trigger_diEle.push_back(el.triggered("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"));
    if (el.gsfTrack()->dxy(PV.position())>0.05 || el.gsfTrack()->dz(PV.position())>0.1) continue;
    
    if (showlog){std::cout<<"prmpt electron"<<std::endl;}
//       	barrel 	endcap
// d0, cm 	0.05 	0.10
// dz, cm 	0.10 	0.20 
  }
  tree_electron_nEle=nEl;


  //////////////////////////////////
  //////////////////////////////////
  /////   DiMuon Selection   ///////
  //////////////////////////////////
  //////////////////////////////////

  int imu1 = -1, imu2 = -1;
  float mupt1, mueta1, muphi1, mupt2, mueta2, muphi2;
  float mu_mass = 0.1057;
  TLorentzVector v1, v2, v;
  tree_Mmumu = 0.;

  if ( nmu >= 2 ) 
  {
    for (int mu = 0; mu < nmu-1; mu++)
    { 
      if ( !tree_muon_isGlobal[mu] ) continue;
      mupt1  = tree_muon_pt[mu];
      if ( mupt1 < 10. ) continue; // Zmu filter
      // if ( abs(tree_muon_dxy[mu]) > 0.1 || abs(tree_muon_dz[mu]) > 0.2 || (!tree_muon_trigger_dimu[mu] && !tree_muon_trigger_isomu[mu]) || !tree_muon_isLoose[mu] ) continue; // muons closed to PV
      // if ( abs(tree_muon_dxy[mu]) > 0.1 || abs(tree_muon_dz[mu]) > 0.2 ||  !tree_muon_isLoose[mu] ) continue; //|| !tree_muon_PFIsoLoose[mu] muons closed to PV
      if ( abs(tree_muon_dxy[mu]) > 0.1 || abs(tree_muon_dz[mu]) > 0.2 ||  !tree_muon_isLoose[mu] || !tree_muon_TkIsoLoose[mu]) continue;

      mueta1 = tree_muon_eta[mu];
      muphi1 = tree_muon_phi[mu];
      v1.SetPtEtaPhiM(mupt1,mueta1,muphi1,mu_mass);

      for ( int mu2=mu+1; mu2<nmu; mu2++) 
        {	    
          if ( !tree_muon_isGlobal[mu2] ) continue;
          mupt2  = tree_muon_pt[mu2];
          if ( mupt2 < 10. ) continue;
          if ( mupt1 < 25. && mupt2 < 25. ) continue; // Zmu Filter
          if ( tree_muon_charge[mu] == tree_muon_charge[mu2] ) continue;
          // // if ( (tree_muon_trigger_isomu[mu] && (abs(tree_muon_dxy[mu2]) > 0.1 || abs(tree_muon_dz[mu2]) > 0.2 || !tree_muon_isLoose[mu2]))  || (abs(tree_muon_dxy[mu2]) > 0.1 || abs(tree_muon_dz[mu2]) > 0.2 || (!tree_muon_trigger_dimu[mu2] && tree_muon_trigger_dimu[mu] ) || !tree_muon_isLoose[mu2]) ) continue;
          // if ( ( (abs(tree_muon_dxy[mu2]) > 0.1 || abs(tree_muon_dz[mu2]) > 0.2 || !tree_muon_isLoose[mu2] )) ) continue;//|| !tree_muon_PFIsoLoose[mu2]
                if ( ( (abs(tree_muon_dxy[mu2]) > 0.1 || abs(tree_muon_dz[mu2]) > 0.2 || !tree_muon_isLoose[mu2] || !tree_muon_TkIsoLoose[mu2])) ) continue;

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

  //---------------EVTS Selection BDT------------------//
  mva_Evts_MET_et           = tree_PFMet_et;

  // mva_Evts_nVtx             = ;
  mva_HT                    = HT_val;
  mva_Evts_ST               = ST;
  mva_Evts_njets            = tree_njet;
  mva_Evts_nmuon            = nmu;
  // mva_Evts_MediumAxes       = ;
  // mva_Evts_LooseAxes        = ;
  // mva_Evts_TightAxes        = ;
  //the other variables are in loops so they can not be put there
  // but here is the list:

    // mva_Evts_MET_et
    // mva_Evts_nTrks = theFormat.tree_TRACK_SIZE;
    // mva_Evts_muon1_pt= theFormat.tree_muon_leadingpt->at(0);
    // mva_Evts_muon2_pt= theFormat.tree_muon_leadingpt2->at(0);
    // mva_Evts_jet1_pt = theFormat.tree_jet_leadingpt->at(0);
    // mva_Evts_jet2_pt = theFormat.tree_jet_leadingpt2->at(0);
    // mva_Evts_muon12_dR = theFormat.tree_muon_muon_dR->at(0);
    // mva_Evts_muon12_dPhi = theFormat.tree_muon_muon_dPhi->at(0);
    // mva_Evts_muon12_dEta = theFormat.tree_muon_muon_dEta->at(0);
    // mva_Evts_jet12_dR = theFormat.tree_jet_jet_dR->at(0);
    // mva_Evts_jet12_dPhi = theFormat.tree_jet_jet_dPhi->at(0);
    // mva_Evts_jet12_dEta  = theFormat.tree_jet_jet_dEta->at(0);
    // mva_Evts_muon_jet_dRmin = theFormat.tree_muon_jet_dRmin->at(0);
    // mva_Evts_muon_jet_dRmax = theFormat.tree_muon_jet_dRmax->at(0);
    // mva_Evts_nVtx = theFormat.tree_Hemi_Vtx_nVtx->at(0);
    // mva_HT = theFormat.tree_HT;
    // mva_Evts_ST = theFormat.tree_ST->at(0);
    // mva_Evts_njets = theFormat.tree_njet;
    // mva_Evts_nmuon 


  //////////////////////////////////
  //////////////////////////////////
  ////////   FILTER CHECK  /////////
  //////////////////////////////////
  //////////////////////////////////
  
  tree_NbrOfZCand = 0;
  tree_Filter = false;
  tree_nTracks = 0;
  tree_nLostTracks = 0;
  tree_TRACK_SIZE = 0;
  tree_nSecInt = 0;
  tree_nV0_reco = 0;

  if ( tree_Mmumu > 60. )                  tree_NbrOfZCand = 1;
//   if ( tree_Mmumu > 60. && HT_val > 180. ) tree_Filter = true;
//$$

//Mu trigger 2018 : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v || HLT_IsoMu24_v
//Mu trigger 2017 : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v || HLT_IsoMu27_v
// Mu trigger 2016 : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL(_DZ)_v || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL(_DZ)_v || HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_(DZ)_v  || HLT_IsoMu24_v || HLT_IsoTkMu24_v 
  
  if ( ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v || HLT_IsoMu24_v  ) 
       && tree_Mmumu > 10. ) tree_Filter = true;
//$$

//$$$$
    tree_GoodMu1 = -1;
    if ( imu1 >= 0 ) { 
      tree_GoodMu1 = 0;
      if ( tree_muon_isTight[imu1] )    tree_GoodMu1 = 1; 
      if ( tree_muon_TkIsoTight[imu1] ) tree_GoodMu1 += 10;
    }
    tree_GoodMu2 = -1; 
    if ( imu2 >= 0 ) { 
      tree_GoodMu2 = 0;
      if ( tree_muon_isTight[imu2] )    tree_GoodMu2 = 1; 
      if ( tree_muon_TkIsoTight[imu2] ) tree_GoodMu2 += 10;
    }
//$$$$
  
  edm::ESHandle<TransientTrackBuilder> theTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTransientTrackBuilder); // Asking for reco collection of PV..
  vector<reco::TransientTrack> BestTracks;
  vector<reco::TransientTrack> posBestTracks;
  vector<reco::TransientTrack> negBestTracks;
  int count =0;
  std::map<size_t , int > trackToAK4SlimmedJetMap;


  if ( tree_Filter ) 
  {

  float dRmuon1_jet_min = 0;
  float dRmuon1_jet_max = 0;
  float dRmuon2_jet_min = 0;
  float dRmuon2_jet_max = 0;
  if (nmu>=2 && tree_njet>=2)
  	{
		float dRmuon1_jet0=Deltar(tree_muon_eta[imu1],tree_muon_phi[imu1],tree_jet_eta[0],tree_jet_phi[0]);
		float dRmuon1_jet1=Deltar(tree_muon_eta[imu1],tree_muon_phi[imu1],tree_jet_eta[1],tree_jet_phi[1]);
		if (dRmuon1_jet0 < dRmuon1_jet1)
			{
				dRmuon1_jet_min = dRmuon1_jet0;
				dRmuon1_jet_max = dRmuon1_jet1;
			}
		else
			{
				dRmuon1_jet_min = dRmuon1_jet1;
				dRmuon1_jet_max = dRmuon1_jet0;
			}
		float dRmuon2_jet0=Deltar(tree_muon_eta[imu2],tree_muon_phi[imu2],tree_jet_eta[0],tree_jet_phi[0]);
		float dRmuon2_jet1=Deltar(tree_muon_eta[imu2],tree_muon_phi[imu2],tree_jet_eta[1],tree_jet_phi[1]);
		if (dRmuon2_jet0 < dRmuon2_jet1)
			{
				dRmuon2_jet_min = dRmuon2_jet0;
				dRmuon2_jet_max = dRmuon2_jet1;
			}
		else
			{
				dRmuon2_jet_min = dRmuon2_jet1;
				dRmuon2_jet_max = dRmuon2_jet0;
			}
	}

  if (dRmuon1_jet_min<dRmuon2_jet_min)
    {
      tree_muon_jet_dRmin.push_back(dRmuon1_jet_min);
	    tree_muon_jet_dRmin.push_back(dRmuon2_jet_min);
    }
  else
    {
      tree_muon_jet_dRmin.push_back(dRmuon2_jet_min);
	    tree_muon_jet_dRmin.push_back(dRmuon1_jet_min);
    }

  if (dRmuon1_jet_max>dRmuon2_jet_max)
    {
      tree_muon_jet_dRmax.push_back(dRmuon1_jet_max);
	    tree_muon_jet_dRmax.push_back(dRmuon2_jet_max);
    }
  else
    {
      tree_muon_jet_dRmax.push_back(dRmuon2_jet_max);
	    tree_muon_jet_dRmax.push_back(dRmuon1_jet_max);
    }	
   mva_Evts_muon_jet_dRmin = tree_muon_jet_dRmin[0];
   mva_Evts_muon_jet_dRmax = tree_muon_jet_dRmax[0];

       /////////////////////////////////////////////////////////
    //-------------------------------------------------------
    // Jets for event axes                                 
    //-------------------------------------------------------
    /////////////////////////////////////////////////////////

    int njet = 0, njet1 = 0, njet2 = 0;
    bool isjet[99], isjet1[99], isjet2[99];//the "real" size is given by the final value of jetidx since non-valid jets are replaced
    float btag[99]={0};
    float btag1[99]={0};
    float btag2[99]={0};
    TLorentzVector vaxis1, vaxis2, vjet[99];
    TLorentzVector Vaxis1, Vaxis2, Vjet[99];
    TLorentzVector V1, V2, V;
    float PtMin = 20;   // (GeV) minimum jet pt is optimum
    float EtaMax = 10; // no cut on eta is optimum
    int jetidx = 0; // : May be in the loop/ not sure it changes anything
    for (const pat::Jet &jet : *jets)    // Loop on jet
    {
      float jet_pt  = jet.pt();
      float jet_eta = jet.eta();
      float jet_phi = jet.phi();
      float jet_px = jet.px();
      float jet_py = jet.py();
      float jet_pz = jet.pz();
      float jet_e = jet.energy();
      float jetp = sqrt(jet_px*jet_px+jet_py*jet_py+jet_pz*jet_pz);
      // float jetmass = sqrt(jet_e*jet_e-jetp*jetp);
      // std::cout<<"jet mass: "<<jetmass<<std::endl;
      // std::cout<<"tree_jet_pt : "<< jet.pt()<<std::endl;
      float DeepFlavourb = jet.bDiscriminator("pfDeepFlavourJetTags:probb");
      float DeepFlavourbb = jet.bDiscriminator("pfDeepFlavourJetTags:probbb");
      float DeepFlavourblep = jet.bDiscriminator("pfDeepFlavourJetTags:problepb");
      float DeepJet = -10.;
      if ( DeepFlavourb > -5 ) DeepJet = DeepFlavourb + DeepFlavourbb + DeepFlavourblep;

      btag[jetidx]  =  DeepJet;
      btag1[jetidx]  =  0;
      btag2[jetidx]  =  0;    
      isjet[jetidx]  = false;
      isjet1[jetidx] = false; // first neutralino jets
      isjet2[jetidx] = false; // second neutralino jets
      v.SetPtEtaPhiM( jet_pt, jet_eta, jet_phi, 0. ); //set the axis
      V.SetPxPyPzE(jet_px,jet_py,jet_pz,jet_e);
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
          v1.SetPtEtaPhiM( tree_muon_pt[imu1], tree_muon_eta[imu1], tree_muon_phi[imu1], 0 );
          V1.SetPxPyPzE(tree_muon_px[imu1],tree_muon_py[imu1],tree_muon_pz[imu1],tree_muon_energy[imu1]);
          v -= v1; // v TLorentzFactor being just above, defined by jet data
          V-= V1;
        }
        if ( deltaR2 < 0.4 )
        {
          v2.SetPtEtaPhiM( tree_muon_pt[imu2], tree_muon_eta[imu2], tree_muon_phi[imu2], 0 );
          V2.SetPxPyPzE(tree_muon_px[imu2],tree_muon_py[imu2],tree_muon_pz[imu2],tree_muon_energy[imu2]);
          v -= v2;
          V-= V2;
        }
        jet_pt  = v.Pt(); //Update jet data by removing the muons information (muons that could be in the jet)
        jet_eta = v.Eta(); //+ we do not want muons data to build the two axis since they come from the PV except the one that come from a top decay
        jet_phi = v.Phi();
        jet_px = V.Px();
        jet_pz = V.Py();
        jet_py = V.Pz();
        jet_e  = V.E();
      }
      
      njet++;
      isjet[jetidx] = true;
      vjet[jetidx] = v; // Only jet data (with  possible muons being removed)
      Vjet[jetidx] = V;
      if ( njet1 == 0 && jet_pt > PtMin && abs(jet_eta) < EtaMax )
      {
        njet1 = 1;
        isjet1[jetidx] = true;
        vaxis1 = v;
        Vaxis1 = V;
        btag1[jetidx] = btag[jetidx];
      }
      jetidx++;
    } // End Loop on jets

    /////////////////////////////////////////////////////////
    //-------------------------------------------------------
    // Event Axes
    //-------------------------------------------------------
    /////////////////////////////////////////////////////////

    // There should be a dependance on the decay channel of the top. in theory, the consutrction of the axes
    // should be improved by the use of the MET+lepton. The current method is good for the hradronic decay of the top
    float dR_axis = 10., dR1 = 10., dR2 = 10.;
    float dRcut_hemis  = 1.5; // subjective choice default is 1.5
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
        Vaxis1 += Vjet[i];
        isjet1[i] = true;
        btag1[i]=btag[i];
        if (btag1[i]>0.7264 && ActivateBtagLog) {std::cout<<"Tight b jet ID"<<std::endl;}
        // 0.2770 : Medium 
        // 0.0494 : loose
      }
      // axis 2
      if ( njet2 == 0 && !isjet1[i] ) {
        njet2 = 1;
        vaxis2 = vjet[i];
        Vaxis2 += Vjet[i];
        isjet2[i] = true;
        btag2[i]=btag[i];
        if (btag2[i]>0.7264 && ActivateBtagLog) {std::cout<<"Tight b jet ID"<<std::endl;}
      }
      else if ( njet2 > 0 && !isjet1[i] && !isjet2[i] && dR2 < dRcut_hemis) {//
        njet2++;
        vaxis2 += vjet[i];
        Vaxis2 += Vjet[i];
        isjet2[i] = true;
        btag2[i]=btag[i];//njet2 insteag of i would also make sense, here there are voids when one of the other conditions above is filled
        if (btag2[i]>0.7264 && ActivateBtagLog) {std::cout<<"Tight b jet ID"<<std::endl;}
      }
    }       // end Loop on jet
    
//$$
//     // force the axes to the true LLP
//     vaxis1 = vneu[0];
//     vaxis2 = vneu[1];
//$$

    ///////////////////////////////
    // Invariant Mass of Axes 
    ///////////////////////////////
    
    //The computation of the invariant mass does not look good when using only the tracks
    //So, we try to do it with the jets and the possible letpon coming out of the neutralino
        float temp_px1 = Vaxis1.Px();
        float temp_py1 = Vaxis1.Py();
        float temp_pz1 = Vaxis1.Pz();
        float temp_e1  = Vaxis1.E();
        float temp_px2 = Vaxis2.Px();
        float temp_py2 = Vaxis2.Py();
        float temp_pz2 = Vaxis2.Pz();
        float temp_e2  = Vaxis2.E();

    //... then, add the lepton contribution
    // float temp_mass_muon = -1;
    // float temp_mass_e    = -1;
    for (int i = 0 ; i<nmu ; i++)
      {
        // if ( !tree_muon_isGlobal[mu] ) continue;
        if ( abs(tree_muon_dxy[i]) < 0.1 || abs(tree_muon_dz[i]) < 0.2 ) continue;//||  !tree_muon_isTight[mu]
        // Avoid prompt leptons
        float dRAxis1muon = Deltar( vaxis1.Eta(), vaxis1.Phi(), tree_muon_eta[i], tree_muon_phi[i] );
        float dRAxis2muon = Deltar( vaxis2.Eta(), vaxis2.Phi(), tree_muon_eta[i], tree_muon_phi[i] );
        if (dRAxis1muon<1.5  )
          {
            temp_px1 += tree_muon_px[i];
            temp_py1 += tree_muon_py[i];
            temp_pz1 += tree_muon_pz[i];
            temp_e1  += tree_muon_energy[i];
            // temp_mass_muon = sqrt(temp_e1*temp_e1-temp_px1*temp_px1+temp_py1*temp_py1+temp_pz1*temp_pz1);
            // if (temp_mass_muon<0)
            //   {
            //     std::cout<<"mass muon & energy : "<<temp_mass_muon<<" && "<<tree_muon_energy[i]<<std::endl;
            //   }
            
          }
        if(dRAxis2muon<1.5)
          {
            temp_px2 += tree_muon_px[i];
            temp_py2 += tree_muon_py[i];
            temp_pz2 += tree_muon_pz[i];
            temp_e2  += tree_muon_energy[i];
            // temp_mass_muon = sqrt(temp_e2*temp_e2-temp_px2*temp_px2+temp_py2*temp_py2+temp_pz2*temp_pz2);
            // if (temp_mass_muon<0)
            //   {
            //     std::cout<<"mass muon & energy : "<<temp_mass_muon<<" && "<<tree_muon_energy[i]<<std::endl;
            //   }
          }
      }

    for (int i = 0 ; i <nEl ; i++)
      {
        if(abs(tree_electron_dxy[i])< 0.05|| abs(tree_electron_dz[i])<0.1) continue; //tree_electron_IsTight
        // Avoid prompt leptons
        float dRAxis1electron = Deltar( vaxis1.Eta(), vaxis1.Phi(), tree_electron_eta[i], tree_electron_phi[i] );
        float dRAxis2electron = Deltar( vaxis2.Eta(), vaxis2.Phi(), tree_electron_eta[i], tree_electron_phi[i] );
        if (dRAxis1electron<1.5  )
          {
            temp_px1 += tree_electron_px[i];
            temp_py1 += tree_electron_py[i];
            temp_pz1 += tree_electron_pz[i];
            temp_e1  += tree_electron_energy[i];
            // temp_mass_e = sqrt(temp_e1*temp_e1-temp_px1*temp_px1+temp_py1*temp_py1+temp_pz1*temp_pz1);
            // if (temp_mass_e<0)
            //   {
            //     std::cout<<"mass electron & energy : "<<temp_mass_e<<" && "<<tree_electron_energy[i]<<std::endl;
            //   }
            
          }
        if(dRAxis2electron<1.5)
          {
            temp_px2 += tree_electron_px[i];
            temp_py2 += tree_electron_py[i];
            temp_pz2 += tree_electron_pz[i];
            temp_e2  += tree_electron_energy[i];
            // temp_mass_e = sqrt(temp_e2*temp_e2-temp_px2*temp_px2+temp_py2*temp_py2+temp_pz2*temp_pz2);
            // if (temp_mass_e<0)
            //   {
            //     std::cout<<"mass electron & energy : "<<temp_mass_e<<" && "<<tree_electron_energy[i]<<std::endl;
            //   }
            
          }
      }
    TLorentzVector TLorentzAxis1(temp_px1,temp_py1,temp_pz1,temp_e1);
    TLorentzVector TLorentzAxis2(temp_px2,temp_py2,temp_pz2,temp_e2);

    float Mass1 = sqrt(TLorentzAxis1.Mag2());
    float Mass2 = sqrt(TLorentzAxis2.Mag2());
    if (isnan(Mass1) || isinf(Mass1))
      {
        Mass1 = 0;
      }
    if (isnan(Mass2) || isinf(Mass2))
      {
        Mass2 = 0;
      }
    tree_Hemi_Mass.push_back(Mass1);
    // std::cout<<"Mass Axis 1 : "<<sqrt(TLorentzAxis1.Mag2())<<std::endl;
    tree_Hemi_Mass.push_back(Mass2);
    // std::cout<<"Mass Axis 2 : "<<sqrt(TLorentzAxis2.Mag2()) <<std::endl;


    ///////////////////////////////
    // compare with neutralino axis
    ///////////////////////////////

    int iLLPrec1 = 1, iLLPrec2 = 2;
    float axis1_eta = vaxis1.Eta();
    float axis1_phi = vaxis1.Phi();
    if ( !runOnData_ )
      {

    if ( neu[0] >= 0 ) dR1 = Deltar( axis1_eta, axis1_phi, Gen_neu1_eta, Gen_neu1_phi ); //dR between reco axis of jets and gen neutralino
    if ( neu[1] >= 0 ) dR2 = Deltar( axis1_eta, axis1_phi, Gen_neu2_eta, Gen_neu2_phi );
    dR_axis = dR1;
      if ( dR2 < dR1 )
        { // make sure that the reco axis defined matches well with the axis of the gen neutralino, if not it is swapped
          iLLPrec1 = 2;
          iLLPrec2 = 1;
          dR_axis = dR2;
        }
      }
    float axis1_dR = dR_axis;
    float axis2_eta = vaxis2.Eta();
    float axis2_phi = vaxis2.Phi();

    float axis2_dR = dR2;
    float dR_axis12 = 10.;
    if ( !runOnData_ )
      {
        if ( njet2 == 0 )
          {  // compute an axis 2 even without jet, by taking the opposite in phi to axis 1
            axis2_eta = axis1_eta;
            axis2_phi = axis1_phi - 3.14159;
            if ( axis1_phi < 0 ) axis2_phi = axis1_phi + 3.14159;
         }
        if ( iLLPrec2 == 1 ) dR_axis = Deltar( axis2_eta, axis2_phi, Gen_neu1_eta, Gen_neu1_phi );
        else                 dR_axis = Deltar( axis2_eta, axis2_phi, Gen_neu2_eta, Gen_neu2_phi );
        axis2_dR = dR_axis;
        dR_axis12 = Deltar(axis1_eta,axis1_phi,axis2_eta,axis2_phi);
      }
        // cout << " njet1 " << njet1 << " and njet2" << njet2 << endl;
        // cout << " axis1_eta " << axis1_eta << " and axis2_eta" << axis2_eta << endl;
        // cout << " axis1_phi " << axis1_phi << " and axis2_phi" << axis2_phi << endl;
        // cout << " axis1_dR " << axis1_dR << " and axis2_dR" << axis2_dR << endl;
        // cout << " dR_axis12 " << dR_axis12 << endl;
      
    ///////////////////////////////////////////////////////
    // Compare between the gen enutralino and the reco axis from gen jets and with the reco axis with reco jets
    ///////////////////////////////////////////////////////
    if ( !runOnData_ )
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
      }
    
      bool LooseAxesStatus = false ;
      bool MediumAxesStatus = false ;
      bool TightAxesStatus = false ;
   if (ActivateBtag)
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

  mva_Evts_LooseAxes        = LooseAxesStatus;
  mva_Evts_MediumAxes       = MediumAxesStatus;
  mva_Evts_TightAxes        = TightAxesStatus;

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


    //-------Applying the EVTS BDT selection
    float EVTS_BDTval = -10;
    EVTS_BDTval = readerEvts->EvaluateMVA( "BDTG" );//Current wp is -0.58 (optimizing S/sqrt(S+B)) but -0.7 is the point that is the maximum value to keep 90% of the signal)
    // To select the events that you want, use  a macro (at the moment)
    tree_Evts_MVAval.push_back(EVTS_BDTval);
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
    // reconstruct secondary interactiosn (photon conversions an dnuclear intractions) with slight modification
    // Basically, only cuts and the badtkhit have been changed from the V0Producer code.
    // Keep in mind that we are using the MINIGeneral Tracks (packedPFCandidate and LostTrack)  

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
    float tkChi2Cut_	= 10.; // 10. by default
    // # Number of valid hits on track >=
    int   tkNHitsCut_	= 3;   // 3 by default
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
//$$
    float kShortMassCut_ = 0.07;  // 0.07 by default
    float lambdaMassCut_ = 0.05;  // 0.05 by default
    float kShortMassSel_ = 0.022;  // track selection cut 0.476 - 0.520
    float lambdaMassSel_ = 0.0060; // track selection cut /: could go for asymmetric selectino
//$$
    int temp_nV0_reco = 0;

    //------------------Selection of good tracks for the vertexing---------------------//

    // fill vectors of TransientTracks and TrackRefs after applying preselection cuts

    std::vector<reco::Track> theTrackRefs;
    std::vector<reco::TransientTrack> theTransTracks;

    std::vector<std::pair<bool,float>> idxMGT; // contains the following informations : the size of this vector is the number of selected tracks passing the cuts
    // The bool says "is used at the end to build a V0 candidate", the float keeps in memory the index of the track from MINIGeneralTracks that is used to build the V0 candiates
//
    std::vector<std::pair<bool,float>> idxSecIntMGT;
//

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

//-------------------------------- loop over tracks and vertex good charged track pairs
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
//$$
      if ( !(theTrackRefs[trd1].pt() > pt_Cut && theTrackRefs[trd1].normalizedChi2() < NChi2_Cut && drsig1 > drSig_Cut) 
	&& !(theTrackRefs[trd2].pt() > pt_Cut && theTrackRefs[trd2].normalizedChi2() < NChi2_Cut && drsig2 > drSig_Cut) ) continue;
//$$
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
//$$
      if ( badTkHit && dca > 0.1 ) continue;
//$$

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
            else                                 theNegativeRefTrack = &*iTrack;
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
            tree_V0_reco_x.push_back(	  K0x);
            tree_V0_reco_y.push_back(	  K0y);
            tree_V0_reco_z.push_back(	  K0z);
            tree_V0_reco_r.push_back(	  TMath::Sqrt(K0x*K0x+K0y*K0y));
            tree_V0_reco_drSig.push_back( distsigXY);
            tree_V0_reco_dzSig.push_back( distsigZ);
            tree_V0_reco_angleXY.push_back(  angleXY);
            tree_V0_reco_angleZ.push_back(   detaPointing);
            tree_V0_reco_NChi2.push_back( theKshort->vertexNormalizedChi2());
            tree_V0_reco_ndf.push_back(   theKshort->vertexNdof());
            tree_V0_reco_mass.push_back(  theKshort->mass());
            tree_V0_reco_pt.push_back(    theKshort->pt());
            tree_V0_reco_eta.push_back(   theKshort->eta());
            tree_V0_reco_phi.push_back(   theKshort->phi());
            tree_V0_reco_source.push_back(1);
            tree_V0_reco_badTkHit.push_back(badTkHit);
            tree_V0_reco_dca.push_back(   dca);
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
            tree_V0_reco_x.push_back(	  L0x);
            tree_V0_reco_y.push_back(	  L0y);
            tree_V0_reco_z.push_back(	  L0z);
            tree_V0_reco_r.push_back(	  TMath::Sqrt(L0x*L0x+L0y*L0y));
            tree_V0_reco_drSig.push_back( distsigXY);
            tree_V0_reco_dzSig.push_back( distsigZ);
            tree_V0_reco_angleXY.push_back(  angleXY);
            tree_V0_reco_angleZ.push_back(   detaPointing);
            tree_V0_reco_NChi2.push_back( theLambda->vertexNormalizedChi2());
            tree_V0_reco_ndf.push_back(   theLambda->vertexNdof());
            tree_V0_reco_mass.push_back(  theLambda->mass());
            tree_V0_reco_pt.push_back(    theLambda->pt());
            tree_V0_reco_eta.push_back(   theLambda->eta());
            tree_V0_reco_phi.push_back(   theLambda->phi());
            tree_V0_reco_source.push_back(2);
            tree_V0_reco_badTkHit.push_back(badTkHit);
            tree_V0_reco_dca.push_back(   dca);
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
            tree_V0_reco_x.push_back(	  L0x);
            tree_V0_reco_y.push_back(	  L0y);
            tree_V0_reco_z.push_back(	  L0z);
            tree_V0_reco_r.push_back(	  TMath::Sqrt(L0x*L0x+L0y*L0y));
            tree_V0_reco_drSig.push_back( distsigXY);
            tree_V0_reco_dzSig.push_back( distsigXYZ);
            tree_V0_reco_angleXY.push_back(  angleXY);
            tree_V0_reco_angleZ.push_back(   detaPointing);
            tree_V0_reco_NChi2.push_back( theLambdaBar->vertexNormalizedChi2());
            tree_V0_reco_ndf.push_back(   theLambdaBar->vertexNdof());
            tree_V0_reco_mass.push_back(  theLambdaBar->mass());
            tree_V0_reco_pt.push_back(    theLambdaBar->pt());
            tree_V0_reco_eta.push_back(   theLambdaBar->eta());
            tree_V0_reco_phi.push_back(   theLambdaBar->phi());
            tree_V0_reco_source.push_back(2);
            tree_V0_reco_badTkHit.push_back(badTkHit);
            tree_V0_reco_dca.push_back(   dca);
            if ( abs(theLambdaBar->mass() - lambdaMass) < lambdaMassSel_ && ActivateV0Veto) { 
              idxMGT[trd1].first = true;
              idxMGT[trd2].first = true;
            }
          }
        }
        
      }
    }
    tree_nV0_reco = temp_nV0_reco;
//-------------------- END OF V0 reconstruction ---------------------//


    //---------------------------------------------------------------//
    //----------------------- Secondary Interactions ----------------//
    //---------------------------------------------------------------//
int SecInt_ntrk10 = 0;
int SecInt_ntrk20 = 0;
//----------------------- loop over tracks and vertex good charged track pairs
    for (unsigned int trd1 = 0; trd1 < theTrackRefs.size()-1; ++trd1) 
    {
//$$
    if ( idxMGT[trd1].first ) continue;
//$$

      int iq1 = theTrackRefs[trd1].charge();
      for (unsigned int trd2 = trd1+1; trd2 < theTrackRefs.size(); ++trd2) 
      {
//$$
      if ( idxMGT[trd2].first ) continue;
//$$
        int iq2 = theTrackRefs[trd2].charge();
        float drsig1 = std::abs(theTrackRefs[trd1].dxy(PV.position())/theTrackRefs[trd1].dxyError());
        float drsig2 = std::abs(theTrackRefs[trd2].dxy(PV.position())/theTrackRefs[trd2].dxyError());
        if ( useBS_ ) {
          drsig1 = std::abs(theTrackRefs[trd1].dxy(bs)/theTrackRefs[trd1].dxyError());
          drsig2 = std::abs(theTrackRefs[trd2].dxy(bs)/theTrackRefs[trd2].dxyError());
	}
//$$
      if ( !(theTrackRefs[trd1].pt() > pt_Cut && theTrackRefs[trd1].normalizedChi2() < NChi2_Cut && drsig1 > drSig_Cut) 
	&& !(theTrackRefs[trd2].pt() > pt_Cut && theTrackRefs[trd2].normalizedChi2() < NChi2_Cut && drsig2 > drSig_Cut) ) continue;
//$$

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
//$$
      if ( dca > tkDCACut_ ) continue;
//$$

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
//$$
      if ( mass > 4. ) continue; 
//$$
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
//$$
      if ( theVtx.normalizedChi2() > 20. ) continue;
//$$
        GlobalPoint vtxPos(theVtx.x(), theVtx.y(), theVtx.z());

        // 2D decay length significance 
	float distsigXY = 10000.;
        SMatrixSym3D  totalCov = PV.covariance() + theVtx.covariance();
        if ( useBS_ ) totalCov = bs.rotatedCovariance3D() + theVtx.covariance();
        SVector3 distVecXY(vtxPos.x()-PV.x(), vtxPos.y()-PV.y(), 0.);
        double distMagXY = ROOT::Math::Mag(distVecXY);
        double sigmaDistMagXY = sqrt(ROOT::Math::Similarity(totalCov, distVecXY)) / distMagXY;
	if ( sigmaDistMagXY > 0. ) distsigXY = distMagXY / sigmaDistMagXY;
//$$
      if ( distsigXY < 50. ) continue;
//$$  	     
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

          //Check track first hit density here => There are too many SecInt reconstructed in Signal MC
          // for (unsigned int trd = 0; trd < theTrackRefs.size(); ++trd)
          //   {
          //     const HitPattern hp2 = theTrackRefs[trd].hitPattern();
          //     uint16_t firsthit2 = hp2.getHitPattern(HitPattern::HitCategory::TRACK_HITS,0);     
          //     //---Creating State to propagate from  TT---//
          //     reco::TransientTrack TTrack2 = theTransientTrackBuilder->build(theTrackRefs[trd]);
          //     GlobalPoint vert2 (theTrackRefs[trd].vx(),theTrackRefs[trd].vy(),theTrackRefs[trd].vz()); // Point where the propagation will start (Reference Point)
          //     const TrajectoryStateOnSurface Surtraj2 = TTrack.stateOnSurface(vert2); // TSOS of this point
          //     const MagneticField* B2 = TTrack2.field(); // 3.8T
          //     AnalyticalPropagator* Prop = new AnalyticalPropagator(B2);
          //     Basic3DVector<float> P3D22(theTrackRefs[trd].vx(),theTrackRefs[trd].vy(),theTrackRefs[trd].vz()); // global frame
          //     Basic3DVector<float> B3DV2(theTrackRefs[trd].px(),theTrackRefs[trd].py(),theTrackRefs[trd].pz()); // global frame 
          //     float Eta2 = theTrackRefs[trd].eta();
          //     float Phi2 = theTrackRefs[trd].phi();
          //     float vz2  = theTrackRefs[trd].vz();
          
          //     std::pair<int,GloballyPositioned<float>::PositionType> FHPosition2 = PHP->Main(firsthit2,Prop,Surtraj2,Eta2,Phi2,vz2,P3D22,B3DV2);
          //     if (sqrt((FHPosition2.second.x()-FHPosition.second.x())*(FHPosition2.second.x()-FHPosition.second.x())
          //               +(FHPosition2.second.y()-FHPosition.second.y())*(FHPosition2.second.y()-FHPosition.second.y())
          //               +(FHPosition2.second.z()-FHPosition.second.z())*(FHPosition2.second.z()-FHPosition.second.z())  )<10)
          //                 {
          //                   SecInt_ntrk10++;
          //                 }
          //     else if (sqrt((FHPosition2.second.x()-FHPosition.second.x())*(FHPosition2.second.x()-FHPosition.second.x())
          //               +(FHPosition2.second.y()-FHPosition.second.y())*(FHPosition2.second.y()-FHPosition.second.y())
          //               +(FHPosition2.second.z()-FHPosition.second.z())*(FHPosition2.second.z()-FHPosition.second.z())  )<20)
          //                 {
          //                   SecInt_ntrk20++;
          //                 }
          //   }
          // tree_SecInt_ntrk10.push_back(SecInt_ntrk10);
          // tree_SecInt_ntrk20.push_back(SecInt_ntrk20);
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
//$$
      if ( badTkHit && dca > 0.1 ) continue;
//$$
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
            else             theRefTrack2 = &*iTrack;
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
//$$
      if ( angleXY < 0.7 ) continue;
//$$                    
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
//$$
      if ( abs(theSecInt->mass()) > 4. ) continue; 
//$$
        float SecInt_x = theSecInt->vertex().x();
        float SecInt_y = theSecInt->vertex().y();
        float SecInt_z = theSecInt->vertex().z();
        float SecInt_r = TMath::Sqrt(SecInt_x*SecInt_x + SecInt_y*SecInt_y);
        float SecInt_d = TMath::Sqrt(SecInt_x*SecInt_x + SecInt_y*SecInt_y+SecInt_z*SecInt_z);
	tree_SecInt_x.push_back(       SecInt_x);
	tree_SecInt_y.push_back(       SecInt_y);
	tree_SecInt_z.push_back(       SecInt_z);
	tree_SecInt_r.push_back(       SecInt_r);
  tree_SecInt_d.push_back(       SecInt_d);
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
        if (!runOnData_)
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
            if      ( dphi < -3.14159 / 2. ) dphi += 3.14159;
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
            if      ( dphi < -3.14159 / 2. ) dphi += 3.14159;
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
}
	bool SecInt_selec = false;
//$$
	if ( distsigXY > 100. && angleXY > 0.9 && theSecInt->mass() < 2. &&
             theSecInt->vertexNormalizedChi2() < 10. ) SecInt_selec = true; //&& dca < 1.
//$$
        tree_SecInt_selec.push_back(   SecInt_selec);

	// tracker active layers
        // PropaHitPattern* NI = new PropaHitPattern();
        int VtxLayerNI = -1;
        if (DetailedMap) {VtxLayerNI = NI->VertexBelongsToTracker(SecInt_r, SecInt_z);}
        else {
          VtxLayerNI = NI->VertexBelongsToBarrelLayer(SecInt_r, SecInt_z);
          if ( VtxLayerNI == 0 ) VtxLayerNI = NI->VertexBelongsToDiskLayer(SecInt_r, SecInt_z);
        }
//$$
	if ( SecInt_selec && ActivateSecIntVeto ) {
	  if ( VtxLayerNI != 0 ) {
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // beam pipe
	  if ( abs(SecInt_z) < 27. && SecInt_r > 2.16 && SecInt_r < 2.26 ) {
	    VtxLayerNI = -1;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // PIXB inner support
	  if ( abs(SecInt_z) < 27. && SecInt_r > 2.49 && SecInt_r < 2.54 ) {
	    VtxLayerNI = -2;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // PIXB outer support
	  if ( SecInt_r > 21.55 && SecInt_r < 21.85 ) {
	    VtxLayerNI = -3;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // PIXB rails
	  if ( abs(SecInt_z) < 27. && abs(SecInt_x) < 10. && 
	       abs(SecInt_y) > 18.9 && abs(SecInt_y) < 19.4 ) {
	    VtxLayerNI = -4;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	  // services : r 18 - 19 and z 29 - 200 according to "other" material map
	  if ( abs(SecInt_z) > 29. && SecInt_r > 18.0 && SecInt_r < 19.0 ) {
	    VtxLayerNI = -3;
            idxSecIntMGT[trd1].first = true;
            idxSecIntMGT[trd2].first = true;
	  }
	}
//$$
        tree_SecInt_layer.push_back(   VtxLayerNI);
      }
    }
//-----------------------------END OF Sec.Int. reconstruction------------------------------------------//
    
//$$$$
//     // std::vector<std::pair<bool,float>> idxMuonMGT;
//     for (unsigned int ipc = 0; ipc < TRACK_SIZE; ipc++) { // loop on all packedPFCandidates + lostTrackss
//       pat::PackedCandidateRef pcref = MINIgeneralTracks[ipc];
//       const reco::Track *trackPcPtr = pcref->bestTrack();
//     if ( !trackPcPtr ) continue;
// 
//       reco::Track tk = *trackPcPtr;
//     if ( tk.hitPattern().numberOfValidHits() == 0 ) continue;
//     if ( tk.charge() == 0 ) continue;
//     // double ipsigXY = std::abs(tk.dxy(PV.position()) / tk.dxyError());
//     // if ( useBS_ ) ipsigXY = std::abs(tk.dxy(bs) / tk.dxyError());
//     // double ipsigZ = std::abs(tk.dz(PV.position()) / tk.dzError());
// 
//   if ( imu1 >= 0 && imu2 >= 0 && tree_muon_pt[imu2] > tree_muon_pt[imu1] ) {
//  
//       for (unsigned int k = 0 ; k<tree_muon_pt.size(); k++)
//         {
//           float dptMu = (tk.pt()-tree_muon_pt[k])/tk.pt();
//           float detaMu = tk.eta()-tree_muon_eta[k];
//           float dphiMu = Deltaphi(tk.phi(),tree_muon_phi[k] );
//           if (abs(dptMu)<0.1 && abs(detaMu)<0.1 && abs(dphiMu)<0.1 && (tree_muon_trigger_isomu[k] || tree_muon_trigger_dimu[k]) && abs(tree_muon_dxy[k]) < 0.1 && abs(tree_muon_dz[k]) < 0.2 &&  tree_muon_isLoose[k] && tree_muon_TkIsoLoose[k])
//             {
//               idxMuonMGT.push_back(make_pair(true,ipc));              
//             }
//         }
//     }
//$$$$


    //---------------------------------------------------------------//
    //------------------------- TRACKS ------------------------------//
    //---------------------------------------------------------------//

    for (unsigned int ipc = 0; ipc < TRACK_SIZE; ipc++) { // loop on all packedPFCandidates + lostTrackss
      pat::PackedCandidateRef pcref = MINIgeneralTracks[ipc];
      float energy = pcref->energy();
      const reco::Track *trackPcPtr = pcref->bestTrack(); // const
//       if ( !trackPcPtr) tree_passesTrkPtr.push_back(0);
    if ( !trackPcPtr ) continue;
//       tree_passesTrkPtr.push_back(1);

      // reject tracks from K0 and Lambda decays and from secondary interactions
      bool Veto_tk = false;
      for (unsigned int j = 0 ; j < idxMGT.size() ; j++) {
          if ( idxSecIntMGT[j].second != idxMGT[j].second ) cout << " ERROR idxSecIntMGT " << endl;
        if ( (idxMGT[j].first || idxSecIntMGT[j].first) && idxMGT[j].second == ipc ) {
	  Veto_tk = true;//keep as true
	  break;
        }
      }
//$$
    if ( Veto_tk ) continue;
//$$

      reco::Track tk = *trackPcPtr;
      tk_nHit	= tk.hitPattern().numberOfValidHits();
      tk_charge = tk.charge();
      tk_pt  = tk.pt();
      tk_eta = tk.eta();
      tk_phi = tk.phi();
      tk_NChi2 = tk.normalizedChi2();
      tk_dxy = tk.dxy(PV.position());
      tk_dxyError = tk.dxyError();
      tk_dz = tk.dz(PV.position());
      tk_dzError = tk.dzError();
      if ( tk_dzError > 0 ) 
        tk_dzSig = abs(tk_dz) / tk_dzError; 
      if ( tk_dxyError > 0 ) 
        tk_drSig = abs(tk_dxy) / tk_dxyError; 
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

//$$
if (RequestHighPurity) // preselection
  {
    if ( !(tk_pt > pt_Cut && tk_NChi2 < NChi2_Cut && tk_drSig > drSig_Cut
     && static_cast<int>(tk.quality(reco::TrackBase::highPurity))) ) continue;
  }
else  
  {
    if ( !(tk_pt > pt_Cut && tk_NChi2 < NChi2_Cut && tk_drSig > drSig_Cut ) ) continue;
  }
//$$

//$$$$
      // reject prompt muon tracks
      bool PromptMuonVeto = false;
      float dptMu1 = 10., detaMu1 = 10., dphiMu1 = 10.;
      float dptMu2 = 10., detaMu2 = 10., dphiMu2 = 10.;
      int iqMu1 = 0, iqMu2 = 0;
      if ( imu1 >= 0 ) {
        dptMu1 = (tk_pt - tree_muon_pt[imu1]) / tk_pt;
        detaMu1 = tk_eta - tree_muon_eta[imu1];
        dphiMu1 = Deltaphi( tk_phi, tree_muon_phi[imu1] );
	iqMu1 = tree_muon_charge[imu1];
      }
      if ( imu2 >= 0 ) {
        dptMu2 = (tk_pt - tree_muon_pt[imu2]) / tk_pt;
        detaMu2 = tk_eta - tree_muon_eta[imu2];
        dphiMu2 = Deltaphi( tk_phi, tree_muon_phi[imu2] );
	iqMu2 = tree_muon_charge[imu2];
      }
      if ( abs(dptMu1)<0.1 && abs(detaMu1)<0.1 && abs(dphiMu1)<0.1 
	   && tk_charge == iqMu1 ) PromptMuonVeto = true;		 
      if ( abs(dptMu2)<0.1 && abs(detaMu2)<0.1 && abs(dphiMu2)<0.1
	   && tk_charge == iqMu2 ) PromptMuonVeto = true;		 
//$$
    if ( PromptMuonVeto ) continue;
//$$
//$$$$

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
//$$
    if ( Yc_tk && ActivateYcVeto ) continue;
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
      tree_track_px.push_back           (tk_px);
      tree_track_py.push_back           (tk_py);
      tree_track_pz.push_back           (tk_pz);
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
      tree_track_dzSig.push_back        (tk_dzSig);
      tree_track_energy.push_back       (tk_e);

// Tracks from pileup
      float tk_dzmin = 100.;
      float tk_dzminSig = 100000.;
      for (unsigned int i = 0; i< primaryVertex->size() ; ++i) {
      if ( i == 0 ) continue;
        float dzTOpu = tk_dz + tree_PV_z - (*primaryVertex)[i].z();
        if ( abs(dzTOpu) < abs(tk_dzmin) ) tk_dzmin = dzTOpu;
        float dzerror = TMath::Sqrt( tk_dzError*tk_dzError + (*primaryVertex)[i].zError()*(*primaryVertex)[i].zError() );
        if ( dzerror > 0. ) {
          if ( abs(dzTOpu / dzerror) < abs(tk_dzminSig) ) tk_dzminSig = dzTOpu / dzerror;
	}
      }
      tree_track_dzTOpu.push_back       (tk_dzmin);
      tree_track_dzSigTOpu.push_back    (tk_dzminSig);

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
//$$
// //     Approximation for the lostTrack since the hitpattern information is not available (only 1160, tracking POG knows about it but do not seem to care)      
//       if ( ipc >= pc->size() ) {
//         if ( abs(tk_eta) < 1. ) firsthit = 1184; // PIXBL4 in barrel
//         else                    firsthit = 1296; // PIXFD2 in forward
//       }      
// //$$
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
      if ( matchTOjet ) {tree_track_iJet.push_back (iJet);tree_track_btag.push_back(btagFromJet);}
      else              {tree_track_iJet.push_back (-1);tree_track_btag.push_back(-1);}

      // match to gen particle from LLP decay

      if (!runOnData_)
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
        }
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

//$$
//     double _pt = -0.33; (0.22)   // for TMVAClassification_BDTG50cm_TT_WTrigger.weights.xml BDTMiniaod
//     double bdtcut = -0.35; (0.7241)   // for TMVAClassification_BDTG50cm_TT.weights.xml BDTMiniaod
//     double bdtcut = -0.0401; // for TMVAbgctau50withnhits.xml BDToldrecoavecalgo
//     double bdtcut = -0.0815; // for TMVAClassification_BDTG50sansalgo.weights.xml BDToldreco
//     double bdtcut =  0.0327; // for TMVAClassification_BDTG50cm_NewSignal.weights.xml BDTreco
//     double bdtcut = -0.1456; // for TMVAClassification_BDTG50cm_HighPurity.weights.xml BDTrecohp
//     double bdtcut = -0.1083; // for TMVAClassification_BDTG_FromBC.weights.xml from BDTminipf
//     double bdtcut = -0.0067; // for TMVAClassification_BDTG50cm_sansntrk10_avecHP.weights.xml BDTrecohpsansntrk10
//     double bdtcut = -10.; // no BDT cut
//     double bdtcut = 0.1624; // for TMVAClassification_BDTG50cm_wVeto.weights.xml
// TMVAClassification_BDTG50cm_V0Veto.weights.xml"),  # BDTMiniAOD //-0.0090
// TMVAClassification_BDTG50cm_V0_YcVeto.weights.xml,  # BDTMiniAOD //0.0372
// TMVAClassification_BDTG50cm_NoVeto.weights.xml,  # BDTMiniAOD ///0.1270
//$$
    double bdtcut = 0.85;        // ttbar ~ 1E-3
    double bdtcut_step2 = 0.0;   // ttbar ~ 1E-2
//$$

    //---------------------------//

    for (int counter_track = 0; counter_track < tree_nTracks; counter_track++) 
    {
      firsthit_X = tree_track_firstHit_x[counter_track];
      firsthit_Y = tree_track_firstHit_y[counter_track];
      firsthit_Z = tree_track_firstHit_z[counter_track];
      pt         = tree_track_pt[counter_track];
      eta        = tree_track_eta[counter_track];
      phi	 = tree_track_phi[counter_track];
      NChi	 = tree_track_NChi2[counter_track];
      nhits	 = tree_track_nHit[counter_track];
      dxy	 = abs(tree_track_dxy[counter_track]);
      dz	 = abs(tree_track_dz[counter_track]);
      drSig	 = tree_track_drSig[counter_track];
      dzSig      = tree_track_dzSig[counter_track];

      dzTopu     = tree_track_dzTOpu[counter_track];
      dzSigTopu  = tree_track_dzSigTOpu[counter_track];
      TibHit     = tree_track_nHitTIB[counter_track] ;
      TobHit     = tree_track_nHitTOB[counter_track] ;
      PixBarHit  = tree_track_nHitPXB[counter_track];
      TecHit     = tree_track_nHitTEC[counter_track];
      // algo	    = tree_track_algo[counter_track];

      ntrk10 = 0, ntrk20 = 0, ntrk30 = 0, ntrk40 = 0;
      isLost     = tree_track_lost[counter_track];
      float ntrk10_lost = 0, ntrk20_lost = 0, ntrk30_lost = 0, ntrk40_lost = 0;
      isinjet = 0.;
      double bdtval = -10.;
      dR = -1.;
      int tracks_axis = 0; // flag to check which axis is the closest from the track

      jet = tree_track_iJet[counter_track];
      if ( jet >= 0 ) isinjet = 1.; /*!*/

      int isFromLLP = tree_track_sim_LLP[counter_track];

      // check the dR between the tracks and the second axis (without any selection on the tracks)
      dR1  = Deltar( eta, phi, axis1_eta, axis1_phi ); // axis1_phi and axis1_eta for the first axis
      dR2  = Deltar( eta, phi, axis2_eta, axis2_phi );
      tracks_axis = 1;
      dR = dR1;
      dRmax = dR2;
      if ( dR2 < dR1 ) { // a restriction could be added on the value of dR to assign the value Tracks_axis  (avoid some background???)
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
      else		tree_track_Hemi_LLP.push_back(0);
      
    } //End loop on all the tracks
    

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
        //	 {
        //	     std::cout<<" Apres Covcor m["<<i<<"]["<<j<<"] <<m[i][j]<<std::endl;
        //	 }
        //   }
      //-------------------end covariance matrix correction-----------//
      }
      else tk = *trackPcPtr;

      int isFromLLP   = tree_track_sim_LLP[counter_track];
      int tracks_axis = tree_track_Hemi[counter_track];
      float track_p = sqrt(tree_track_px[counter_track]*tree_track_px[counter_track]+tree_track_py[counter_track]*tree_track_py[counter_track]+tree_track_pz[counter_track]*tree_track_pz[counter_track]);
      float track_e = tree_track_energy[counter_track];
      double bdtval   = tree_track_MVAval[counter_track];
      if ((tree_track_energy[counter_track]*tree_track_energy[counter_track]-track_p*track_p)<0)
        {
          track_e = sqrt(track_p*track_p+0.01946);//pi0 mass
        }
      TLorentzVector vTrack(tree_track_px[counter_track],tree_track_py[counter_track],tree_track_pz[counter_track],track_e);



      if ( bdtval > bdtcut ) {//tight wp
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

      if ( bdtval > bdtcut_step2 ) {//loose wp
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
 if (!runOnData_)
  {

  
  if ( tree_nLLP > 0 ) {
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
    tree_LLP_Mass.push_back(sqrt(Total4Vector1.Mag2()));
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
    
    tree_LLP_Mass.push_back(sqrt(Total4Vector1.Mag2()));
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
  } // tree_LLP>0
  } // !runonData

      //-----------------------------------Vertexing-----------------------------------------------//
      //                               4 steps Vertexing                                           //
      // These 4 steps can be divided into 2 => (1,2) & (3,4)                                      //
      //        (1,2) : Corresponds to the Tight WP  of the Track level  BDT                       //
      //        (3,4) : Corresponds to the Loose WP of the Track level  BDT 
      //     => 1 & 3 are using the classic Adaptive Vertex Fitter (AVF) to build the vertices     //
      //     => 2 & 4 are using the Iterative Adaptive Vertex Fitter (IAVF)to build the vertices   //
      //          => Description of the IAVF at step 2                                             //
      //Process : Takes a collection of displaced track (from the BDT) as an input, then           //
      //          a vertex is built using the AVF or the IAVF if needed. No vertex can be obtained //
      //          out of these 4 steps                                                             //
      //                                                                                           //
      //Reconstruction Criteria : We require the vertex to have a chi2 per D.O.F between 0 and 10  //                                                                    //
      //                          and to be matched to a generated vertex when looking at signal MC//
      //                                                                                           //
      //-------------------------------------------------------------------------------------------//
     
    //  Warning :  Sorry for the people reading the vertexing code, it's not easy to read and understand. I hope that the comments will be enough :D

    //--------------------------- FIRST HEMISPHERE WITH MVA -------------------------------------//

    static AdaptiveVertexFitter 
    theFitter_Vertex_Hemi1_mva(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
    theFitter_Vertex_Hemi1_mva.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 1 step 1-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    
    Vtx_ntk = displacedTracks_Hemi1_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -10.;
    Vtx_step = 0;
    GlobalError posError;
    float MeanWeight =0;
    TransientVertex displacedVertex_Hemi1_mva;
    std::vector<float> Vtx1_Weights;
    std::vector<unsigned int> Vtx1_index;
    float DCA_VTX_Meand = 0;
    int badtkhit_index = -1;
    float tempMeanWeight=0;
    int ntracks = 0;
    float tempchi2 = -10.;
    float tempx = -100.;
    float tempy = -100.;
    float tempz = -100.;      

    if ( Vtx_ntk > 1 && ActivateStep1)
    {
      DCA_VTX_Meand = 0;
      badtkhit_index = -1;
      bool success = false;
      MeanWeight=0;
      tempMeanWeight=0;
      std::vector<TransientTrack> vTT;
          for (int k = 0 ; k <Vtx_ntk-1;k++)
            {
              for (int p = k+1 ; p < Vtx_ntk ; p++)
                {
                  vTT.push_back(displacedTracks_Hemi1_mva[p]);
                  vTT.push_back(displacedTracks_Hemi1_mva[k]);  
                  ntracks = 2;
                  TransientVertex TV = theFitter_Vertex_Hemi1_mva.vertex(vTT); // We take the first "good-looking" seed to start
                  if ( TV.isValid())
                    {
                      for (int m = 0; m < ntracks; m++) // we check thaat both tracks have their first hit "after" the vertex
                        {
                          if(Track_FirstHit_Hemi1_mva[m].first == true) continue; //first hit of lost track is biaised
                          float PosFH = sqrt((Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z()));
                          float PosVtx1 = sqrt((TV.position().x()-PV.x())*(TV.position().x()-PV.x())+(TV.position().y()-PV.y())*(TV.position().y()-PV.y())+(TV.position().z()-PV.z())*(TV.position().z()-PV.z()));
                          if (PosFH>PosVtx1) 
                            {
                              success = true; 
                              tempchi2 = TV.normalisedChiSquared();
	    	                      tempx=TV.position().x();
	    	                      tempy=TV.position().y();
	    	                      tempz=TV.position().z();
                              posError = TV.positionError();
                              continue;//continue not useful
                            }
                          else{badtkhit_index = m;}// we keep in memory the index of the track that does not have a godd first hit
                        }

                      if (Vtx_ntk == 2 && !success){break;}// removing 1 track gives no other option
                      else if (success)
                        {
                          Vtx1_index.clear();
                          Vtx1_index.push_back(k);
                          Vtx1_index.push_back(p);
                          for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                            {
                              if (m == k || m == p) continue;
                              ntracks++;
                              vTT.push_back(displacedTracks_Hemi1_mva[m]);
                              TransientVertex updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                              if ( !updatedTV.isValid() ) 
	                              {  
	    	                          vTT.pop_back();
                                  ntracks--;
	    	                          updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  tempMeanWeight=0;
                                  Vtx1_Weights.clear();
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	        continue;
	                              }
                              if ( updatedTV.isValid() ) 
	                              {
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  float TPosFH = sqrt((Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z()));
                                  float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                                  if (TPosFH>TPosVtx1 || Track_FirstHit_Hemi1_mva[m].first == true) {success = true;}//continue not useful
                                  else  
                                    {
                                      vTT.pop_back();
                                      ntracks--;
	    	                              updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                                      tempchi2 = updatedTV.normalisedChiSquared();
	    	                              tempx=updatedTV.position().x();
	    	                              tempy=updatedTV.position().y();
	    	                              tempz=updatedTV.position().z();
                                      posError = updatedTV.positionError();

                                      tempMeanWeight=0;
                                      Vtx1_Weights.clear();
                                      
                                      for(int i = 0; i<ntracks;i++)
                                        {
                                          tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                          Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                          if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                        }
                        	            continue;
                                    }
                                  tempMeanWeight=0;
                                  Vtx1_Weights.clear();
                                  Vtx1_index.push_back(m);
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
	                              }
                            }
                        }
                        //----------------------------//
                      else if (Vtx_ntk > 2 && !success)
                        {
                          Vtx1_index.clear();
                          if (badtkhit_index == k) {vTT.erase(vTT.begin());ntracks--;Vtx1_index.push_back(p);}
                          else if (badtkhit_index == p) {vTT.erase(vTT.end());ntracks--;Vtx1_index.push_back(k);}
                          else {Vtx1_index.push_back(k);Vtx1_index.push_back(p);}
                          for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                            {
                              if (m == k || m == p || m == badtkhit_index) continue;// we take care to not take into account the track with a wrong first hit
                              ntracks++;//++
                              tempMeanWeight=0;
                             
                              vTT.push_back(displacedTracks_Hemi1_mva[m]);
                             
                              TransientVertex updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                              if ( !updatedTV.isValid() ) 
	                              {  
	    	                          vTT.pop_back();
                                  ntracks--;
                                  if (vTT.size()<2) continue;
	    	                          updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx1_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  if (ntracks<2){success = false;}
                                  if (ntracks>=2){success = true;}
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	        continue;
	                              }
                              if ( updatedTV.isValid() ) 
	                              {
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  float TPosFH = sqrt((Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z()));
                                  float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                                  if (TPosFH>TPosVtx1 || Track_FirstHit_Hemi1_mva[m].first == true ) {success = true;}//continue not useful
                                  else  
                                    {
                                      vTT.pop_back();
                                      ntracks--;
                                      if (vTT.size()<2) continue;
	    	                              updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                                      tempchi2 = updatedTV.normalisedChiSquared();
	    	                              tempx=updatedTV.position().x();
	    	                              tempy=updatedTV.position().y();
	    	                              tempz=updatedTV.position().z();
                                      posError = updatedTV.positionError();
                                      Vtx1_Weights.clear();
                                      tempMeanWeight=0;
                                      
                                      if (ntracks<2){success = false;}
                                      if (ntracks>=2){success = true;}
                                      for(int i = 0; i<ntracks;i++)
                                        {
                                          tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                          Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                          if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                        }
                        	            continue;
                                    }

                                   Vtx1_Weights.clear();
                                   tempMeanWeight=0;
                                  Vtx1_index.push_back(m);
                                  success = true;
                                for(int i = 0; i<ntracks;i++)
                                  {
                                    tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                    Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                    if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                  }
	                              }
                            }
                        }

                        // We should have a Vertex after these conditions
                        Vtx_ntk = ntracks;
                        Vtx_chi = tempchi2;
                        Vtx_x = tempx;
                        Vtx_y = tempy;
                        Vtx_z = tempz;
                        Vtx_step = 1;
                        GlobalPoint RECOvtxPos(Vtx_x, Vtx_y, Vtx_z);
                        float DCA_VTX_Meand = 0;
                        for (int k = 0; k< Vtx_ntk; k++)
                          {

                            TrajectoryStateClosestToPoint DCA_Vtx = vTT[k].trajectoryStateClosestToPoint(RECOvtxPos);
                            if ( DCA_Vtx.isValid() ) // Be careful, all tracks are considered when looking at the DCA,
                            //but one could look of the wieghts of the track at the same time
                              {// The positions are given in the Global frame
                                float pca_Vtx_x = DCA_Vtx.position().x();
                                float pca_Vtx_y = DCA_Vtx.position().y();
                                float pca_Vtx_z = DCA_Vtx.position().z();
                                float refPoint_x = DCA_Vtx.referencePoint().x();
                                float refPoint_y = DCA_Vtx.referencePoint().y();
                                float refPoint_z = DCA_Vtx.referencePoint().z();
                                float DCA_Vtx_x = refPoint_x-pca_Vtx_x;
                                float DCA_Vtx_y = refPoint_y-pca_Vtx_y;
                                float DCA_Vtx_z = refPoint_z-pca_Vtx_z;
                                float DCA_VTX_r = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y);
                                float DCA_VTX_d = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y+DCA_Vtx_z*DCA_Vtx_z);
                                DCA_VTX_Meand+=DCA_VTX_d;
                                tree_Hemi_Vtx_track_DCA_x.push_back(DCA_Vtx_x);
                                tree_Hemi_Vtx_track_DCA_y.push_back(DCA_Vtx_y);
                                tree_Hemi_Vtx_track_DCA_z.push_back(DCA_Vtx_z);
                                tree_Hemi_Vtx_track_DCA_r.push_back(DCA_VTX_r);
                                tree_Hemi_Vtx_track_DCA_d.push_back(DCA_VTX_d);
                      
                              }
                          }
                        DCA_VTX_Meand = DCA_VTX_Meand/(float)(Vtx_ntk);
                        if (MeanWeight==0)//<=>only two tracks in the valid vertex
                          {
                            Vtx1_Weights.clear();
                            for(int i = 0; i<ntracks;i++)
                              {
                                MeanWeight+=TV.trackWeight(vTT[i]);
                                Vtx1_Weights.push_back(TV.trackWeight(vTT[i]));
                              }
                          }
                    }
                    else
                      {
                        ntracks=0;
                        vTT.clear();
                      }
                    if ( success ) break;
                }
              if(success){break;}
            }
                if (showlog){std::cout<<"success Hemi1 step 1 : "<<success<<std::endl;}

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
      //          Chi2, we keep going until there are no more tracks left                          //
      //                                                                                           //
      //Degrees of freedom: -LowerLimit (Chi2 restriction)                                         //
      //                    -UpperLimit (Chi2 restriction)                                         //
      //                    -BDtValue (Restriction on the BDT value of the tracks that formed  //
      //                     the vertex)                                                           //
      //                                                                                           //
      //Why? : Address the drop in efficiency due to MiniAOD datatier, turns out it improves the   //
      //       basic AVF implementation (even in RECO/AOD)                                         //
      //                                                                                           //
      //PS: Maximum efficiency is reached for MiniAOD when using the covariance matrix correction  //
      //-------------------------------------------------------------------------------------------//


    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 1 step 2-----------------------------------------//
    //------------------------------------------------------------------------------------------------// 
    ntracks    = -2;
    tempchi2 = -10.;
    tempx = -100.;
    tempy = -100.;
    tempz = -100.;
    bool badVtx = false;
    
    tempMeanWeight=0;
    if ( (Vtx_chi < 0. || Vtx_chi > 10.) ) badVtx = true;
    if ( badVtx  && displacedTracks_Hemi1_mva.size() > 1 && (IterAVF || ActivateStep2) ) 
    {
      MeanWeight=0;
      DCA_VTX_Meand = 0;
      bool success = false;
      std::vector<TransientTrack> vTT;
      tempchi2 = -10.;
      tempMeanWeight=0;
      Vtx_ntk = displacedTracks_Hemi1_mva.size();
      for ( int p = 1; p < Vtx_ntk; p++ )
      {
        for  ( int k = 0; k < p; k++ ) // take pairs of tracks of highest BDT value
        {
          vTT.push_back(displacedTracks_Hemi1_mva[p]);
          vTT.push_back(displacedTracks_Hemi1_mva[k]);
          ntracks = 2;
          TransientVertex TV = theFitter_Vertex_Hemi1_mva.vertex(vTT); // We take the first "good-looking" seed to start
          

          if ( TV.isValid() && TV.normalisedChiSquared()>0 && TV.normalisedChiSquared()<10 ) 
            {
	            tempchi2 = TV.normalisedChiSquared();
	            tempx = TV.position().x();
	            tempy = TV.position().y();
	            tempz = TV.position().z();
              posError = TV.positionError();
	            success = true;

              for (int m = 0; m < ntracks; m++) // We then add track by track to the vertex and check the validity of the vertex
                {
                  if(Track_FirstHit_Hemi1_mva[m].first == true) continue; //first hit of lost track is biaised
                  float PosFH = sqrt((Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z()));
                  float PosVtx1 = sqrt((TV.position().x()-PV.x())*(TV.position().x()-PV.x())+(TV.position().y()-PV.y())*(TV.position().y()-PV.y())+(TV.position().z()-PV.z())*(TV.position().z()-PV.z()));
                  if (PosFH>PosVtx1) 
                    {
                      success = true; 
                      tempchi2 = TV.normalisedChiSquared();
	    	              tempx=TV.position().x();
	    	              tempy=TV.position().y();
	    	              tempz=TV.position().z();
                      posError = TV.positionError();
                      continue;
                    }//continue not useful
                  else{badtkhit_index = m;}
                }
                if (Vtx_ntk == 2 && !success){break;}// removing 1 track gives no other option
              //fin construction seed 
              else if (success)
                {
                  Vtx1_index.clear();
                  Vtx1_index.push_back(k);
                  Vtx1_index.push_back(p);
                  for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                    {
                      if (m == k || m == p) continue;
                      ntracks++;
                      tempMeanWeight=0;
                      vTT.push_back(displacedTracks_Hemi1_mva[m]);
                      TransientVertex updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                      if ( !updatedTV.isValid() ) 
	                      {  
	    	                  vTT.pop_back();
	    	                  ntracks--;
	    	                  updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
	    	                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                  tempx=updatedTV.position().x();
	    	                  tempy=updatedTV.position().y();
        	                tempz=updatedTV.position().z();
                          posError = updatedTV.positionError();
                          tempMeanWeight=0;
                          Vtx1_Weights.clear();
                                    
                          for(int i = 0; i<ntracks;i++)
                            {
                              tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                              Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                              if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                            }
	    	                  continue;
	                      } 
                      if ( updatedTV.isValid() ) 
	                      {
                          tempchi2 = updatedTV.normalisedChiSquared();
	    	                  tempx=updatedTV.position().x();
	    	                  tempy=updatedTV.position().y();
	    	                  tempz=updatedTV.position().z();
                          posError = updatedTV.positionError();
                          float TPosFH = sqrt((Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z()));
                          float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                          if (TPosFH>TPosVtx1 || Track_FirstHit_Hemi1_mva[m].first == true) {success = true;}//continue not useful
                          else  
                            {
                              vTT.pop_back();
                              ntracks--;
	    	                      updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                              tempchi2 = updatedTV.normalisedChiSquared();
	    	                      tempx=updatedTV.position().x();
	    	                      tempy=updatedTV.position().y();
	    	                      tempz=updatedTV.position().z();
                              posError = updatedTV.positionError();
                              Vtx1_Weights.clear();
                              
                              tempMeanWeight=0;
                              for(int i = 0; i<ntracks;i++)
                                {
                                  tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                  Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                  if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                }
                        	    continue;
                            }
                          Vtx1_Weights.clear();
                          Vtx1_index.push_back(m);
                          tempMeanWeight=0;
                          for(int i = 0; i<ntracks;i++)
                            {
                              tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                              Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                              if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                            }
	                      }
                    } // end loop on the other tracks
                } // end of success
              else if (Vtx_ntk > 2 && !success)
                {
                  Vtx1_index.clear();
                  if (badtkhit_index == k) {vTT.erase(vTT.begin());ntracks--;Vtx1_index.push_back(p);}
                  else if (badtkhit_index == p) {vTT.erase(vTT.end());ntracks--;Vtx1_index.push_back(k);}
                  else {Vtx1_index.push_back(k);Vtx1_index.push_back(p);}
                  for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                    {
                      if (m == k || m == p || m == badtkhit_index) continue;
                      ntracks++;
                      tempMeanWeight=0;
                      
                      vTT.push_back(displacedTracks_Hemi1_mva[m]);
                      if (vTT.size()<2) continue;//should not hapen
                      TransientVertex updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                      if ( !updatedTV.isValid() ) 
	                      {  
	    	                  vTT.pop_back();
                          ntracks--;
                           if (vTT.size()<2) continue;
	    	                  updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                          tempchi2 = updatedTV.normalisedChiSquared();
	    	                  tempx=updatedTV.position().x();
	    	                  tempy=updatedTV.position().y();
	    	                  tempz=updatedTV.position().z();
                          posError = updatedTV.positionError();
                          Vtx1_Weights.clear();
                          tempMeanWeight=0;
                          
                          if (ntracks<2){success = false;}
                          if (ntracks>=2){success = true;}
                          for(int i = 0; i<ntracks;i++)
                            {
                              tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                              Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                              if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                            }
                          continue;
	                      }
                      if ( updatedTV.isValid() ) 
	                      {
                          tempchi2 = updatedTV.normalisedChiSquared();
	    	                  tempx=updatedTV.position().x();
	    	                  tempy=updatedTV.position().y();
	    	                  tempz=updatedTV.position().z();
                          posError = updatedTV.positionError();
                          float TPosFH = sqrt((Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi1_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi1_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi1_mva[m].second.Z()-PV.z()));
                          float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                          if (TPosFH>TPosVtx1 || Track_FirstHit_Hemi1_mva[m].first == true ) {success = true;}//continue not useful
                          else  
                            {
                              vTT.pop_back();
                              ntracks--;
                              if (vTT.size()<2) continue;
	    	                      updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                              tempchi2 = updatedTV.normalisedChiSquared();
	    	                      tempx=updatedTV.position().x();
	    	                      tempy=updatedTV.position().y();
	    	                      tempz=updatedTV.position().z();
                              posError = updatedTV.positionError();
                              Vtx1_Weights.clear();
                              tempMeanWeight=0;
                              
                              if (ntracks<2){success = false;}
                              if (ntracks>=2){success = true;}
                              for(int i = 0; i<ntracks;i++)
                                {
                                  tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                  Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                  if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                }
                        	    continue;
                            }
                          Vtx1_Weights.clear();
                          Vtx1_index.push_back(m);
                          tempMeanWeight = 0;
                          success = true;
                          for(int i = 0; i<ntracks;i++)
                            {
                              tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                              Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                              if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                            }
	                      }
                    }
                }
              Vtx_ntk = ntracks;
              Vtx_chi = tempchi2;
              Vtx_x = tempx;
              Vtx_y = tempy;
              Vtx_z = tempz;
              Vtx_step = 2;
              GlobalPoint RECOvtxPos(Vtx_x, Vtx_y, Vtx_z);
              float DCA_VTX_Meand = 0;
              for (int k = 0; k< Vtx_ntk; k++)
                {
                  // MeanWeight+=displacedVertex_Hemi1_mva.trackWeight(displacedTracks_Hemi1_mva[p]);
                  // Vtx1_Weights.push_back(displacedVertex_Hemi1_mva.trackWeight(displacedTracks_Hemi1_mva[p]));
                  // tree_Hemi_Vtx_trackWeight.push_back(displacedVertex_Hemi1_mva.trackWeight(displacedTracks_Hemi1_mva[p]));
                  TrajectoryStateClosestToPoint DCA_Vtx = vTT[k].trajectoryStateClosestToPoint(RECOvtxPos);
                  if ( DCA_Vtx.isValid() ) // Be careful, all tracks are considered when looking at the DCA,
                  //but one should look of the wieghts of the track at the same time
                    {// The positions are given in the Global frame
                      float pca_Vtx_x = DCA_Vtx.position().x();
                      float pca_Vtx_y = DCA_Vtx.position().y();
                      float pca_Vtx_z = DCA_Vtx.position().z();
                      // to check 
                      float refPoint_x = DCA_Vtx.referencePoint().x();
                      float refPoint_y = DCA_Vtx.referencePoint().y();
                      float refPoint_z = DCA_Vtx.referencePoint().z();
                      float DCA_Vtx_x = refPoint_x-pca_Vtx_x;
                      float DCA_Vtx_y = refPoint_y-pca_Vtx_y;
                      float DCA_Vtx_z = refPoint_z-pca_Vtx_z;
                      float DCA_VTX_r = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y);
                      float DCA_VTX_d = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y+DCA_Vtx_z*DCA_Vtx_z);
                      DCA_VTX_Meand+=DCA_VTX_d;
                      tree_Hemi_Vtx_track_DCA_x.push_back(DCA_Vtx_x);
                      tree_Hemi_Vtx_track_DCA_y.push_back(DCA_Vtx_y);
                      tree_Hemi_Vtx_track_DCA_z.push_back(DCA_Vtx_z);
                      tree_Hemi_Vtx_track_DCA_r.push_back(DCA_VTX_r);
                      tree_Hemi_Vtx_track_DCA_d.push_back(DCA_VTX_d);
                      
                    }
                }
              DCA_VTX_Meand = DCA_VTX_Meand/(float)(Vtx_ntk);
              // tree_Hemi_Vtx_track_MeanDCA_d.push_back(DCA_VTX_Meand);
              if (MeanWeight==0)//<=>only two tracks in the valid vertex
                  {
                    Vtx1_Weights.clear();
                    for(int i = 0; i<ntracks;i++)
                      {
                        MeanWeight+=TV.trackWeight(vTT[i]);
                        Vtx1_Weights.push_back(TV.trackWeight(vTT[i]));
                      }
                  }
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
            if (showlog){std::cout<<"success Hemi1 step 2 : "<<success<<std::endl;}

    }

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 1 step 3-----------------------------------------//
    //------------------------------------------------------------------------------------------------// 

    TransientVertex displacedVertex_step2_Hemi1;
    static AdaptiveVertexFitter theFitter_Vertex_step2_Hemi1(
    	       GeometricAnnealing ( sigmacut, Tini, ratio ), 
    	       DefaultLinearizationPointFinder(),
    	       KalmanVertexUpdator<5>(), 
    	       KalmanVertexTrackCompatibilityEstimator<5>(), 
    	       KalmanVertexSmoother() );
    theFitter_Vertex_step2_Hemi1.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
    badVtx = false;

    if ( (Vtx_chi < 0. || Vtx_chi > 10.) && ActivateStep3 ) badVtx = true;
    if ( badVtx && displacedTracks_step2_Hemi1.size() > 1)
    {
      DCA_VTX_Meand = 0;
      badtkhit_index = -1;
      bool success = false;
      MeanWeight=0;
      tempMeanWeight=0;
      std::vector<TransientTrack> vTT;
      Vtx_ntk = displacedTracks_step2_Hemi1.size();
      Vtx_chi = -10.;
      // displacedVertex_step2_Hemi1 = theFitter_Vertex_step2_Hemi1.vertex(displacedTracks_step2_Hemi1);
          for (int k = 0 ; k <Vtx_ntk-1;k++)
            {
              for (int p = k+1 ; p < Vtx_ntk ; p++)
                {
                  vTT.push_back(displacedTracks_step2_Hemi1[p]);
                  vTT.push_back(displacedTracks_step2_Hemi1[k]);
                  ntracks = 2;
                  TransientVertex TV = theFitter_Vertex_step2_Hemi1.vertex(vTT); // We take the first "good-looking" seed to start
                  if ( TV.isValid())
                    {
                      for (int m = 0; m < ntracks; m++) // We then add track by track to the vertex and check the validity of the vertex
                        {
                          if(Track_FirstHit_step2_Hemi1[m].first == true) continue; //first hit of lost track is biaised
                          float PosFH = sqrt((Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z()));
                          float PosVtx1 = sqrt((TV.position().x()-PV.x())*(TV.position().x()-PV.x())+(TV.position().y()-PV.y())*(TV.position().y()-PV.y())+(TV.position().z()-PV.z())*(TV.position().z()-PV.z()));
                          if (PosFH>PosVtx1) 
                            {
                              success = true; 
                              tempchi2 = TV.normalisedChiSquared();
	    	                      tempx=TV.position().x();
	    	                      tempy=TV.position().y();
	    	                      tempz=TV.position().z();
                              posError = TV.positionError();
                              continue;
                            }//continue not useful
                          else{badtkhit_index = m;}
                        }

                      if (Vtx_ntk == 2 && !success){break;}// removing 1 track gives no other option
                      else if (success)
                        {
                          Vtx1_index.clear();
                          Vtx1_index.push_back(k);
                          Vtx1_index.push_back(p);
                          for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                            {
                              if (m == k || m == p) continue;
                              ntracks++;
                              vTT.push_back(displacedTracks_step2_Hemi1[m]);
                              TransientVertex updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                              if ( !updatedTV.isValid() ) 
	                              {  
	    	                          vTT.pop_back();
                                  ntracks--;
	    	                          updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx1_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	        continue;
	                              }
                              if ( updatedTV.isValid() ) 
	                              {
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  float TPosFH = sqrt((Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z()));
                                  float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                                  if (TPosFH>TPosVtx1 || Track_FirstHit_step2_Hemi1[m].first == true) {success = true;}//continue not useful
                                  else  
                                    {
                                      vTT.pop_back();
                                      ntracks--;
	    	                              updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                                      tempchi2 = updatedTV.normalisedChiSquared();
	    	                              tempx=updatedTV.position().x();
	    	                              tempy=updatedTV.position().y();
	    	                              tempz=updatedTV.position().z();
                                      posError = updatedTV.positionError();
                                      Vtx1_Weights.clear();
                                      tempMeanWeight=0;
                                      
                                      for(int i = 0; i<ntracks;i++)
                                        {
                                          tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                          Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                          if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                        }
                        	            continue;
                                    }
                                  tempMeanWeight=0;
                                  Vtx1_Weights.clear();
                                  Vtx1_index.push_back(m);
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
	                              }
                            }
                        }
                      else if (Vtx_ntk > 2 && !success)
                        {
                          Vtx1_index.clear();
                          if (badtkhit_index == k) {vTT.erase(vTT.begin());ntracks--;Vtx1_index.push_back(p);}
                          else if (badtkhit_index == p) {vTT.erase(vTT.end());ntracks--;Vtx1_index.push_back(k);}
                          else {Vtx1_index.push_back(k);Vtx1_index.push_back(p);}
                          for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                            {
                              if (m == k || m == p || m == badtkhit_index) continue;
                              ntracks++;
                              tempMeanWeight=0;
                              
                              vTT.push_back(displacedTracks_step2_Hemi1[m]);
                              TransientVertex updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                              if ( !updatedTV.isValid() ) 
	                              {  
	    	                          vTT.pop_back();
                                  ntracks--;
                                   if (vTT.size()<2) continue;
	    	                          updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx1_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  if (ntracks<2){success = false;}
                                  if (ntracks>=2){success = true;}
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	        continue;
	                              }
                              if ( updatedTV.isValid() ) 
	                              {
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  float TPosFH = sqrt((Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z()));
                                  float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                                  if (TPosFH>TPosVtx1 || Track_FirstHit_step2_Hemi1[m].first == true ) {success = true;}//continue not useful
                                  else  
                                    {
                                      vTT.pop_back();
                                      ntracks--;
                                       if (vTT.size()<2) continue;
	    	                              updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                                      tempchi2 = updatedTV.normalisedChiSquared();
	    	                              tempx=updatedTV.position().x();
	    	                              tempy=updatedTV.position().y();
	    	                              tempz=updatedTV.position().z();
                                      posError = updatedTV.positionError();
                                      Vtx1_Weights.clear();
                                      tempMeanWeight=0;
                                      
                                      if (ntracks<2){success = false;}
                                      if (ntracks>=2){success = true;}
                                      for(int i = 0; i<ntracks;i++)
                                        {
                                          tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                          Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                          if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                        }
                        	            continue;
                                    }
                                  Vtx1_Weights.clear();
                                  Vtx1_index.push_back(m);
                                  success = true;
                                for(int i = 0; i<ntracks;i++)
                                  {
                                    tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                    Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                    if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                  }
	                              }
                            }
                        }

                        // We should have a Vertex after these conditions
                        Vtx_ntk = ntracks;
                        Vtx_chi = tempchi2;
                        Vtx_x = tempx;
                        Vtx_y = tempy;
                        Vtx_z = tempz;
                        Vtx_step = 3;
                        GlobalPoint RECOvtxPos(Vtx_x, Vtx_y, Vtx_z);
                        float DCA_VTX_Meand = 0;
                        for (int k = 0; k< Vtx_ntk; k++)
                          {
                            // MeanWeight+=displacedVertex_step2_Hemi1.trackWeight(displacedTracks_step2_Hemi1[p]);
                            // Vtx1_Weights.push_back(displacedVertex_step2_Hemi1.trackWeight(displacedTracks_step2_Hemi1[p]));
                            // tree_Hemi_Vtx_trackWeight.push_back(displacedVertex_step2_Hemi1.trackWeight(displacedTracks_step2_Hemi1[p]));
                            TrajectoryStateClosestToPoint DCA_Vtx = vTT[k].trajectoryStateClosestToPoint(RECOvtxPos);
                            if ( DCA_Vtx.isValid() ) // Be careful, all tracks are considered when looking at the DCA,
                            //but one could look of the wieghts of the track at the same time
                              {// The positions are given in the Global frame
                                float pca_Vtx_x = DCA_Vtx.position().x();
                                float pca_Vtx_y = DCA_Vtx.position().y();
                                float pca_Vtx_z = DCA_Vtx.position().z();
                                float refPoint_x = DCA_Vtx.referencePoint().x();
                                float refPoint_y = DCA_Vtx.referencePoint().y();
                                float refPoint_z = DCA_Vtx.referencePoint().z();
                                float DCA_Vtx_x = refPoint_x-pca_Vtx_x;
                                float DCA_Vtx_y = refPoint_y-pca_Vtx_y;
                                float DCA_Vtx_z = refPoint_z-pca_Vtx_z;
                                float DCA_VTX_r = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y);
                                float DCA_VTX_d = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y+DCA_Vtx_z*DCA_Vtx_z);
                                DCA_VTX_Meand+=DCA_VTX_d;
                                tree_Hemi_Vtx_track_DCA_x.push_back(DCA_Vtx_x);
                                tree_Hemi_Vtx_track_DCA_y.push_back(DCA_Vtx_y);
                                tree_Hemi_Vtx_track_DCA_z.push_back(DCA_Vtx_z);
                                tree_Hemi_Vtx_track_DCA_r.push_back(DCA_VTX_r);
                                tree_Hemi_Vtx_track_DCA_d.push_back(DCA_VTX_d);
                      
                              }
                          }
                        DCA_VTX_Meand = DCA_VTX_Meand/(float)(Vtx_ntk);
                        if (MeanWeight==0)//<=>only two tracks in the valid vertex
                          {
                            Vtx1_Weights.clear();
                            for(int i = 0; i<ntracks;i++)
                              {
                                MeanWeight+=TV.trackWeight(vTT[i]);
                                Vtx1_Weights.push_back(TV.trackWeight(vTT[i]));
                              }
                          }
                    }
                    else
                      {
                        ntracks=0;
                        vTT.clear();
                      }
                      if ( success ) break;
                }
              if(success){break;}
            }
                if (showlog){std::cout<<"success Hemi1 step 3 : "<<success<<std::endl;}

    }

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 1 step 4-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    ntracks    = -2;
    tempchi2 = -10.;
    tempx = -100.;
    tempy = -100.;
    tempz = -100.;
    badVtx = false;
    
    if ( (Vtx_chi < 0. || Vtx_chi > 10.)  ) badVtx = true;
    if ( badVtx  && displacedTracks_step2_Hemi1.size() > 1 && (ActivateStep4 || IterAVF) ) //IterAVF est redondant avec AcitvateStep3
    {
      MeanWeight=0;
      bool success = false;
      std::vector<TransientTrack> vTT;
      tempchi2 = -10.;
      tempMeanWeight=0;
      DCA_VTX_Meand = 0;
      Vtx_ntk = displacedTracks_step2_Hemi1.size();
      for ( int p = 1; p < Vtx_ntk; p++ )
      {
        for  ( int k = 0; k < p; k++ ) // take pairs of tracks of highest BDT value
        {
          vTT.push_back(displacedTracks_step2_Hemi1[p]);
          vTT.push_back(displacedTracks_step2_Hemi1[k]);
          ntracks = 2;
          TransientVertex TV = theFitter_Vertex_step2_Hemi1.vertex(vTT); // We take the first "good-looking" seed to start
          
          if ( TV.isValid() && TV.normalisedChiSquared()>0 && TV.normalisedChiSquared()<10 ) 
            {
	            tempchi2 = TV.normalisedChiSquared();
	            tempx = TV.position().x();
	            tempy = TV.position().y();
	            tempz = TV.position().z();
              posError = TV.positionError();
	            success = true; 

              for (int m = 0; m < ntracks; m++) // We then add track by track to the vertex and check the validity of the vertex
                {
                  if(Track_FirstHit_step2_Hemi1[m].first == true) continue; //first hit of lost track is biaised
                  float PosFH = sqrt((Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z()));
                  float PosVtx1 = sqrt((TV.position().x()-PV.x())*(TV.position().x()-PV.x())+(TV.position().y()-PV.y())*(TV.position().y()-PV.y())+(TV.position().z()-PV.z())*(TV.position().z()-PV.z()));
                  if (PosFH>PosVtx1) 
                    {
                      success = true; 
                      tempchi2 = TV.normalisedChiSquared();
	    	              tempx=TV.position().x();
	    	              tempy=TV.position().y();
	    	              tempz=TV.position().z();
                      posError = TV.positionError();
                      continue;
                    }//continue not useful
                  else{badtkhit_index = m;}
                }
            if (Vtx_ntk == 2 && !success){break;}// removing 1 track gives no other option
              //fin construction seed 
            else if (success)
              {
                Vtx1_index.clear();
                Vtx1_index.push_back(k);
                Vtx1_index.push_back(p);
                for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                  {
                    if (m == k || m == p) continue;
                    ntracks++;
                    tempMeanWeight=0;
                    
                    vTT.push_back(displacedTracks_step2_Hemi1[m]);
                    TransientVertex updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                    if ( !updatedTV.isValid() ) 
	                    {  
	    	                vTT.pop_back();
	    	                ntracks--;
	    	                updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
	    	                tempchi2 = updatedTV.normalisedChiSquared();
	    	                tempx=updatedTV.position().x();
	    	                tempy=updatedTV.position().y();
        	              tempz=updatedTV.position().z();
                        posError = updatedTV.positionError();
                        Vtx1_Weights.clear();
                        tempMeanWeight=0;
                        
                        for(int i = 0; i<ntracks;i++)
                          {
                            tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                            Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                            if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                          }
	    	                continue;
	                    } 
                    if ( updatedTV.isValid() ) 
	                    {
                        tempchi2 = updatedTV.normalisedChiSquared();
	    	                tempx=updatedTV.position().x();
	    	                tempy=updatedTV.position().y();
	    	                tempz=updatedTV.position().z();
                        posError = updatedTV.positionError();
                        float TPosFH = sqrt(( Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())*( Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())+( Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())*( Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z())*( Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z()));
                        float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                        if (TPosFH>TPosVtx1 ||  Track_FirstHit_step2_Hemi1[m].first == true) {success = true;}//continue not useful
                          else  
                            {
                              vTT.pop_back();
                              ntracks--;
	    	                      updatedTV = theFitter_Vertex_Hemi1_mva.vertex(vTT);
                                   tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx1_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	    continue;
                            }
                        Vtx1_Weights.clear();
                        Vtx1_index.push_back(m);
                        
                        for(int i = 0; i<ntracks;i++)
                          {
                            tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                            Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                            if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                          }
	                    }
                  } // end loop on the other tracks
              }

            else if (Vtx_ntk > 2 && !success)
                {
                  Vtx1_index.clear();
                  if (badtkhit_index == k) {vTT.erase(vTT.begin());ntracks--;Vtx1_index.push_back(p);}
                  else if (badtkhit_index == p) {vTT.erase(vTT.end());ntracks--;Vtx1_index.push_back(k);}
                  else {Vtx1_index.push_back(k);Vtx1_index.push_back(p);}
                  for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                    {
                      if (m == k || m == p || m == badtkhit_index) continue;
                      ntracks++;
                      tempMeanWeight=0;
                      
                      vTT.push_back(displacedTracks_step2_Hemi1[m]);
                      
                      TransientVertex updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                      if ( !updatedTV.isValid() ) 
	                      {  
	    	                  vTT.pop_back();
                          ntracks--;
                          if (vTT.size()<2) continue;
                          updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                          tempchi2 = updatedTV.normalisedChiSquared();
	    	                  tempx=updatedTV.position().x();
	    	                  tempy=updatedTV.position().y();
	    	                  tempz=updatedTV.position().z();
                          posError = updatedTV.positionError();
                          Vtx1_Weights.clear();
                          tempMeanWeight=0;
                          
                          if (ntracks<2){success = false;}
                          if (ntracks>=2){success = true;}
                          for(int i = 0; i<ntracks;i++)
                            {
                              tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                              Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                              if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                            }
                          continue;
	                      }
                      if ( updatedTV.isValid() ) 
	                      {
                          tempchi2 = updatedTV.normalisedChiSquared();
	    	                  tempx=updatedTV.position().x();
	    	                  tempy=updatedTV.position().y();
	    	                  tempz=updatedTV.position().z();
                          posError = updatedTV.positionError();
                          float TPosFH = sqrt((Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi1[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi1[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi1[m].second.Z()-PV.z()));
                          float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                          if (TPosFH>TPosVtx1 || Track_FirstHit_step2_Hemi1[m].first == true ) {success = true;}//c
                          else  
                            {
                              vTT.pop_back();
                              ntracks--;
                              if (vTT.size()<2) continue;
	    	                      updatedTV = theFitter_Vertex_step2_Hemi1.vertex(vTT);
                              tempchi2 = updatedTV.normalisedChiSquared();
	    	                      tempx=updatedTV.position().x();
	    	                      tempy=updatedTV.position().y();
	    	                      tempz=updatedTV.position().z();
                              posError = updatedTV.positionError();
                              Vtx1_Weights.clear();
                              tempMeanWeight=0;
                              
                              if (ntracks<2){success = false;}
                              if (ntracks>=2){success = true;}
                              for(int i = 0; i<ntracks;i++)
                                {
                                  tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                  Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                  if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                }
                        	    continue;
                            }
                          Vtx1_Weights.clear();
                          Vtx1_index.push_back(m);
                          success = true;
                          for(int i = 0; i<ntracks;i++)
                            {
                              tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                              Vtx1_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                              if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                            }
	                      }
                    }
                }
              Vtx_ntk = ntracks;
              Vtx_chi = tempchi2;
              Vtx_x = tempx;
              Vtx_y = tempy;
              Vtx_z = tempz;
              Vtx_step = 4;
               GlobalPoint RECOvtxPos(Vtx_x, Vtx_y, Vtx_z);
              //  float DCA_VTX_Meand = 0;
              // float totalpx = 0 ;
              // float totalpy = 0 ;
              // float totalpz = 0 ;
              // float totale = 0 ;
              for (int k = 0; k< Vtx_ntk; k++)
                {
                  // MeanWeight+=displacedVertex_Hemi1_mva.trackWeight(displacedTracks_step2_Hemi1[p]);
                  // Vtx1_Weights.push_back(displacedVertex_Hemi1_mva.trackWeight(displacedTracks_step2_Hemi1[p]));
                  // tree_Hemi_Vtx_trackWeight.push_back(displacedVertex_Hemi1_mva.trackWeight(displacedTracks_step2_Hemi1[p]));
                  TrajectoryStateClosestToPoint DCA_Vtx = vTT[k].trajectoryStateClosestToPoint(RECOvtxPos);
                  if ( DCA_Vtx.isValid() ) // Be careful, all tracks are considered when looking at the DCA,
                  //but one should look of the wieghts of the track at the same time
                    {// The positions are given in the Global frame
                      float pca_Vtx_x = DCA_Vtx.position().x();
                      float pca_Vtx_y = DCA_Vtx.position().y();
                      float pca_Vtx_z = DCA_Vtx.position().z();
                      // to check 
                      float refPoint_x = DCA_Vtx.referencePoint().x();
                      float refPoint_y = DCA_Vtx.referencePoint().y();
                      float refPoint_z = DCA_Vtx.referencePoint().z();
                      float DCA_Vtx_x = refPoint_x-pca_Vtx_x;
                      float DCA_Vtx_y = refPoint_y-pca_Vtx_y;
                      float DCA_Vtx_z = refPoint_z-pca_Vtx_z;
                      float DCA_VTX_r = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y);
                      float DCA_VTX_d = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y+DCA_Vtx_z*DCA_Vtx_z);
                      DCA_VTX_Meand+=DCA_VTX_d;
                      tree_Hemi_Vtx_track_DCA_x.push_back(DCA_Vtx_x);
                      tree_Hemi_Vtx_track_DCA_y.push_back(DCA_Vtx_y);
                      tree_Hemi_Vtx_track_DCA_z.push_back(DCA_Vtx_z);
                      tree_Hemi_Vtx_track_DCA_r.push_back(DCA_VTX_r);
                      tree_Hemi_Vtx_track_DCA_d.push_back(DCA_VTX_d);
                    }
                }
              DCA_VTX_Meand = DCA_VTX_Meand/(float)(Vtx_ntk);
              if (MeanWeight==0)//<=>only two tracks in the valid vertex
                  {
                    Vtx1_Weights.clear();
                    for(int i = 0; i<ntracks;i++)
                      {
                        MeanWeight+=TV.trackWeight(vTT[i]);
                        Vtx1_Weights.push_back(TV.trackWeight(vTT[i]));
                      }
                  }
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
      if (showlog){std::cout<<"success Hemi1 step 4 : "<<success<<std::endl;}
    }

    float Vtx_chi1 = Vtx_chi;
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
    tree_Hemi_Vtx_xError.push_back(posError.cxx());
    tree_Hemi_Vtx_yError.push_back(posError.cyy());
    tree_Hemi_Vtx_zError.push_back(posError.czz());
    float theta_Vtx = tan(sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y)/abs(Vtx_z)) ;
    float eta_Vtx = -TMath::Log(tan(theta_Vtx/2));
    if (Vtx_z<0){eta_Vtx = -eta_Vtx;}
    tree_Hemi_Vtx_eta.push_back(eta_Vtx);
    recX = Vtx_x - tree_PV_x;
    recY = Vtx_y - tree_PV_y;
    recZ = Vtx_z - tree_PV_z;
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_Hemi_Vtx_dist.push_back( recD );
    tree_Hemi_Vtx_MeantrackWeight.push_back(MeanWeight);
    tree_Hemi_Vtx_track_MeanDCA_d.push_back(DCA_VTX_Meand);
    int nVertex = 0;
    float posx1 = Vtx_x;
    float posy1 = Vtx_y;
    float posz1 = Vtx_z;
    if ( Vtx_step>0 && Vtx_chi<10 && Vtx_chi>0 ){nVertex++;}

    float ddok, ddbad;
    float ping1 = 0;

  if (!runOnData_)
    {

    
  if ( tree_nLLP > 0 ) {
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
  }  
    }     
    // Vertex Analysis Step
    // float LooseWP  = 0.0494;
    // float MediumWP = 0.2770;
    // float TightWP  = 0.7264;
    int BtagGood_Hemi1 = 0; // Number of b-tagged hets over Medium WP in the Hemi1
    float temp_px = 0;
    float temp_py = 0;
    float temp_pz = 0 ;
    float temp_e = 0 ;

// ROOT::Math::PxPyPzEVector
    TLorentzVector Total4Vector1(0,0,0,0);
      //--------------- B-tagging-----------------------//
      //Vtx1_weights Vtx1_index TrackInfo_Hemi1_mva
    if(Vtx_step==1 || Vtx_step==2)
      {
        // if(Vtx_chi>0 && Vtx_chi<10)
        //   {
            for (unsigned int i = 0 ; i <TrackInfo_Hemi1_mva.size(); i++)
              {
                for (unsigned int j = 0 ; j < Vtx1_index.size(); j++)
                  {
                    if (i == Vtx1_index[j])
                      {
                        if (Vtx1_Weights[j]>0.5)
                          {
                            temp_px = TrackInfo_Hemi1_mva[i].second.Px();
                            temp_py = TrackInfo_Hemi1_mva[i].second.Py();
                            temp_pz = TrackInfo_Hemi1_mva[i].second.Pz();
                            temp_e  = TrackInfo_Hemi1_mva[i].second.E();
                            TLorentzVector TLorentzTrack(temp_px,temp_py,temp_pz,temp_e);
                            Total4Vector1 +=TLorentzTrack;
                          }
                      }
                  }
              }
          // }
      }

    else if (Vtx_step == 3 || Vtx_step==4)
      {

            for (unsigned int i = 0 ; i <TrackInfo_Hemi1_mva.size(); i++)
              {
                for (unsigned int j = 0 ; j < Vtx1_index.size(); j++)
                  {
                    if (i == Vtx1_index[j])
                      {
                        if (Vtx1_Weights[j]>0.5)
                          {
                            temp_px = TrackInfo_step2_Hemi1[i].second.Px();
                            temp_py = TrackInfo_step2_Hemi1[i].second.Py();
                            temp_pz = TrackInfo_step2_Hemi1[i].second.Pz();
                            temp_e  = TrackInfo_step2_Hemi1[i].second.E();
                            TLorentzVector TLorentzTrack(temp_px,temp_py,temp_pz,temp_e);
                            Total4Vector1 +=TLorentzTrack;
                          }
                      }
                  }
              }
      }


    tree_Hemi_Vtx_BTag.push_back(BtagGood_Hemi1);
      // -------------------- End  of  B-Tagging --------------------//
    // ROOT::Math::PxPyPzEVector for later version of root (later than 6.14 at least)
    tree_Hemi_Vtx_Mass.push_back(sqrt(Total4Vector1.Mag2()));

    // -------------------- End Of Invariant Mass ------------------------//

    float Vtx1_ntk = Vtx_ntk;
    float Vtx1_chi = Vtx_chi;
    float Vtx1_step = Vtx_step;
    float Vtx1_r = sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y);
    float Vtx1_z = Vtx_z;
    float Vtx1_MTW = MeanWeight;
    float Vtx1_Mass = sqrt(Total4Vector1.Mag2());
    float H1_Mass = sqrt(TLorentzAxis1.Mag2());
    float Vtx1_dist = recD;

    float Vtx1_MeanDCA = DCA_VTX_Meand;
if (Vtx1_Weights.size() != Vtx1_index.size()){std::cout<<"size Vtx1_weights and Vtx1_index and  ntracks and chi and step: "<<Vtx1_Weights.size()<<" and "<<Vtx1_index.size()<<" and "<<Vtx_ntk<<" and "<<Vtx_chi<<" and "<<Vtx_step<<std::endl;
}
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

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 2 step 1-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    Vtx_ntk = displacedTracks_Hemi2_mva.size();
    Vtx_x = -100.;
    Vtx_y = -100.;
    Vtx_z = -100.;
    Vtx_chi = -10.;
    Vtx_step = 0;
    MeanWeight = 0;
    DCA_VTX_Meand = 0;
    TransientVertex displacedVertex_Hemi2_mva;
    std::vector<float>  Vtx2_Weights;
    std::vector<unsigned int>    Vtx2_index;
    ntracks = 0;
    tempchi2 = -10.;
    tempx = -100.;
    tempy = -100.;
    tempz = -100.;    
//-------------test Paul----------------//
    if ( displacedTracks_Hemi2_mva.size() > 1 && ActivateStep1)
    {
      DCA_VTX_Meand = 0;
      badtkhit_index = -1;
      bool success = false;
      MeanWeight=0;
      tempMeanWeight=0;
      std::vector<TransientTrack> vTT;
          for (int k = 0 ; k <Vtx_ntk-1;k++)
            {
              for (int p = k+1 ; p < Vtx_ntk ; p++)
                {
                  vTT.push_back(displacedTracks_Hemi2_mva[p]);
                  vTT.push_back(displacedTracks_Hemi2_mva[k]);
                  ntracks = 2;
                  
                  TransientVertex TV = theFitter_Vertex_Hemi2_mva.vertex(vTT); // We take the first "good-looking" seed to start
                  if ( TV.isValid())
                    {
                      for (int m = 0; m < ntracks; m++) // We then add track by track to the vertex and check the validity of the vertex
                        {
                          if(Track_FirstHit_Hemi2_mva[m].first == true) continue; //first hit of lost track is biaised
                          float PosFH = sqrt((Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z()));
                          float PosVtx1 = sqrt((TV.position().x()-PV.x())*(TV.position().x()-PV.x())+(TV.position().y()-PV.y())*(TV.position().y()-PV.y())+(TV.position().z()-PV.z())*(TV.position().z()-PV.z()));
                          if (PosFH>PosVtx1) 
                            {
                              success = true; 
                              tempchi2 = TV.normalisedChiSquared();
	    	                      tempx=TV.position().x();
	    	                      tempy=TV.position().y();
	    	                      tempz=TV.position().z();
                              posError = TV.positionError();
                              continue;
                            }//continue not useful
                          else{badtkhit_index = m;}
                        }
                      // std::cout<<"success "<<success<<std::endl; //bug when success is false
                      if (Vtx_ntk == 2 && !success){break;}// removing 1 track gives no other option
                      else if (success)
                        {
                          Vtx2_index.clear();
                          Vtx2_index.push_back(k);
                          Vtx2_index.push_back(p);
                          for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                            {
                              if (m == k || m == p) continue;
                              ntracks++;
                              vTT.push_back(displacedTracks_Hemi2_mva[m]);
                              TransientVertex updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                              if ( !updatedTV.isValid() ) 
	                              {  
	    	                          vTT.pop_back();
                                  ntracks--;
	    	                          updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx2_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	        continue;
	                              }
                              if ( updatedTV.isValid() ) 
	                              {
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  float TPosFH = sqrt((Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z()));
                                  float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                                  if (TPosFH>TPosVtx1 || Track_FirstHit_Hemi2_mva[m].first == true  ) {success = true;}//continue not useful
                                  else  
                                    {
                                      vTT.pop_back();
                                      ntracks--;
                                      updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                                      tempchi2 = updatedTV.normalisedChiSquared();
	    	                              tempx=updatedTV.position().x();
	    	                              tempy=updatedTV.position().y();
	    	                              tempz=updatedTV.position().z();
                                      posError = updatedTV.positionError(); 
                                      Vtx2_Weights.clear();
                                      tempMeanWeight=0;
                                      
                                      for(int i = 0; i<ntracks;i++)
                                        {
                                          tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                          Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                          if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                        }
                        	            continue;
                                    }
                                  tempMeanWeight=0;
                                  Vtx2_Weights.clear();
                                  Vtx2_index.push_back(m);
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
	                              }
                            }
                        }
                      else if (Vtx_ntk > 2 && !success)
                        {
                          Vtx2_index.clear();
                          if (badtkhit_index == k) {vTT.erase(vTT.begin());ntracks--;Vtx2_index.push_back(p);}
                          else if (badtkhit_index == p) {vTT.erase(vTT.end());ntracks--;Vtx2_index.push_back(k);}
                          else {Vtx2_index.push_back(k);Vtx2_index.push_back(p);}
                          for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                            {
                              if (m == k || m == p || m == badtkhit_index) continue;
                              ntracks++;
                              tempMeanWeight=0;
                              
                              vTT.push_back(displacedTracks_Hemi2_mva[m]);
                              TransientVertex updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                              if ( !updatedTV.isValid() ) 
	                              {  
	    	                          vTT.pop_back();
                                  ntracks--;
                                  if (vTT.size()<2) continue;
	    	                          updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx2_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  if (ntracks<2){success = false;}
                                  if (ntracks>=2){success = true;}
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	        continue;
	                              }
                              if ( updatedTV.isValid() ) 
	                              {
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  float TPosFH = sqrt((Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z()));
                                  float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                                  if (TPosFH>TPosVtx1 || Track_FirstHit_Hemi2_mva[m].first == true ) {success = true;}//continue not useful
                                  else  
                                    {
                                      vTT.pop_back();
                                      ntracks--;
                                      if (vTT.size()<2) continue;
	    	                              updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                                      tempchi2 = updatedTV.normalisedChiSquared();
	    	                              tempx=updatedTV.position().x();
	    	                              tempy=updatedTV.position().y();
	    	                              tempz=updatedTV.position().z();
                                      posError = updatedTV.positionError();
                                      Vtx2_Weights.clear();
                                      tempMeanWeight=0;
                                      
                                      if (ntracks<2){success = false;}
                                      if (ntracks>=2){success = true;}
                                      for(int i = 0; i<ntracks;i++)
                                        {
                                          tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                          Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                          if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                        }
                        	            continue;
                                    }
                                    Vtx2_Weights.clear();
                                  Vtx2_index.push_back(m);
                                  tempMeanWeight=0;
                                  success = true;
                                for(int i = 0; i<ntracks;i++)
                                  {
                                    tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                    Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                    if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                  }
	                              }
                            }
                        }

                        // We should have a Vertex after these conditions
                        Vtx_ntk = ntracks;
                        Vtx_chi = tempchi2;
                        Vtx_x = tempx;
                        Vtx_y = tempy;
                        Vtx_z = tempz;
                        Vtx_step = 1;
                        GlobalPoint RECOvtxPos(Vtx_x, Vtx_y, Vtx_z);
                        float DCA_VTX_Meand = 0;
                        for (int k = 0; k< Vtx_ntk; k++)
                          {
                            // MeanWeight+=displacedVertex_Hemi2_mva.trackWeight(displacedTracks_Hemi2_mva[p]);
                            // Vtx2_Weights.push_back(displacedVertex_Hemi2_mva.trackWeight(displacedTracks_Hemi2_mva[p]));
                            // tree_Hemi_Vtx_trackWeight.push_back(displacedVertex_Hemi2_mva.trackWeight(displacedTracks_Hemi2_mva[p]));
                            TrajectoryStateClosestToPoint DCA_Vtx = vTT[k].trajectoryStateClosestToPoint(RECOvtxPos);
                            if ( DCA_Vtx.isValid() ) // Be careful, all tracks are considered when looking at the DCA,
                            //but one could look of the wieghts of the track at the same time
                              {// The positions are given in the Global frame
                                float pca_Vtx_x = DCA_Vtx.position().x();
                                float pca_Vtx_y = DCA_Vtx.position().y();
                                float pca_Vtx_z = DCA_Vtx.position().z();
                                float refPoint_x = DCA_Vtx.referencePoint().x();
                                float refPoint_y = DCA_Vtx.referencePoint().y();
                                float refPoint_z = DCA_Vtx.referencePoint().z();
                                float DCA_Vtx_x = refPoint_x-pca_Vtx_x;
                                float DCA_Vtx_y = refPoint_y-pca_Vtx_y;
                                float DCA_Vtx_z = refPoint_z-pca_Vtx_z;
                                float DCA_VTX_r = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y);
                                float DCA_VTX_d = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y+DCA_Vtx_z*DCA_Vtx_z);
                                DCA_VTX_Meand+=DCA_VTX_d;
                                tree_Hemi_Vtx_track_DCA_x.push_back(DCA_Vtx_x);
                                tree_Hemi_Vtx_track_DCA_y.push_back(DCA_Vtx_y);
                                tree_Hemi_Vtx_track_DCA_z.push_back(DCA_Vtx_z);
                                tree_Hemi_Vtx_track_DCA_r.push_back(DCA_VTX_r);
                                tree_Hemi_Vtx_track_DCA_d.push_back(DCA_VTX_d);
                      
                              }
                          }
                        DCA_VTX_Meand = DCA_VTX_Meand/(float)(Vtx_ntk);
                        if (MeanWeight==0)//<=>only two tracks in the valid vertex
                          {
                            Vtx2_Weights.clear();
                            for(int i = 0; i<ntracks;i++)
                              {
                                MeanWeight+=TV.trackWeight(vTT[i]);
                                Vtx2_Weights.push_back(TV.trackWeight(vTT[i]));
                              }
                          }
                    }
                    else
                      {
                        ntracks=0;
                        vTT.clear();
                      }
                      if ( success ) break;
                }
              if(success){break;}
            }
            if (showlog){std::cout<<"success Hemi2 step 1 : "<<success<<std::endl;}
    }

    

    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 2 step 2-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    ntracks   = -2;
    tempchi2  = -10.;
    tempx = -100.;
    tempy = -100.;
    tempz = -100.;
    badVtx = false;
    //---------------------//
    if ( (Vtx_chi < 0. || Vtx_chi > 10.)  ) badVtx = true;
    if ( badVtx  && displacedTracks_Hemi2_mva.size() > 1 && (ActivateStep2 || IterAVF)) 
      {
        MeanWeight=0;
        tempMeanWeight=0;
        DCA_VTX_Meand = 0;
        bool success = false;
        std::vector<TransientTrack> vTT;
        tempchi2 = -10.;
        Vtx_ntk = displacedTracks_Hemi2_mva.size();
        for ( int p = 1; p < Vtx_ntk; p++ )
          {
            for  ( int k = 0; k < p; k++ ) // take pairs of tracks of highest BDT value
              {
                vTT.push_back(displacedTracks_Hemi2_mva[p]);
                vTT.push_back(displacedTracks_Hemi2_mva[k]);
                ntracks = 2;
                
                TransientVertex TV = theFitter_Vertex_Hemi2_mva.vertex(vTT); // We take the first "good-looking" seed to start
                if ( TV.isValid() && TV.normalisedChiSquared()>0 && TV.normalisedChiSquared()<10 ) 
                  {
	                  tempchi2 = TV.normalisedChiSquared();
	                  tempx = TV.position().x();
	                  tempy = TV.position().y();
	                  tempz = TV.position().z();
                    posError =TV.positionError();
                    success = true;
                    for (int m = 0; m < ntracks; m++) // We then add track by track to the vertex and check the validity of the vertex
                      {
                        if(Track_FirstHit_Hemi2_mva[m].first == true) continue; //first hit of lost track is biaised
                        float PosFH = sqrt((Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z()));
                        float PosVtx1 = sqrt((TV.position().x()-PV.x())*(TV.position().x()-PV.x())+(TV.position().y()-PV.y())*(TV.position().y()-PV.y())+(TV.position().z()-PV.z())*(TV.position().z()-PV.z()));
                        if (PosFH>PosVtx1) 
                          {
                            success = true; 
                            tempchi2 = TV.normalisedChiSquared();
	    	                    tempx=TV.position().x();
	    	                    tempy=TV.position().y();
	    	                    tempz=TV.position().z();
                            posError = TV.positionError();
                            continue;
                          }//continue not useful
                        else{badtkhit_index = m;}
                      }
                if (Vtx_ntk == 2 && !success){break;}// removing 1 track gives no other option
                else if (success)
                  {
                    Vtx2_index.clear();
                    Vtx2_index.push_back(k);
                    Vtx2_index.push_back(p);
                    for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                      {
                        if (m == k || m == p) continue;
                        ntracks++;
                        tempMeanWeight=0;
                        
                        vTT.push_back(displacedTracks_Hemi2_mva[m]);
                        TransientVertex updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                        if ( !updatedTV.isValid() ) 
	                        {  
	    	                    vTT.pop_back();
	    	                    ntracks--;
	    	                    updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
	    	                    tempchi2 = updatedTV.normalisedChiSquared();
	    	                    tempx=updatedTV.position().x();
	    	                    tempy=updatedTV.position().y();
        	                  tempz=updatedTV.position().z();
                            posError = updatedTV.positionError();
                            Vtx2_Weights.clear();
                            tempMeanWeight=0;
                            
                            for(int i = 0; i<ntracks;i++)
                              {
                                tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                              }
	    	                    continue;
	                        } 
                        if ( updatedTV.isValid() ) 
	                        {
                            tempchi2 = updatedTV.normalisedChiSquared();
	    	                    tempx=updatedTV.position().x();
	    	                    tempy=updatedTV.position().y();
	    	                    tempz=updatedTV.position().z();
                            posError = updatedTV.positionError();
                            float TPosFH = sqrt((Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z()));
                          float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                          if (TPosFH>TPosVtx1 || Track_FirstHit_Hemi2_mva[m].first == true) {success = true;}//continue not useful
                          else  
                            {
                              vTT.pop_back();
                              ntracks--;
	    	                      updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                                                                tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx2_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  if (ntracks<2){success = false;}
                                  if (ntracks>=2){success = true;}
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	    continue;
                            }
                            Vtx2_Weights.clear();
                            Vtx2_index.push_back(m);
                            for(int i = 0; i<ntracks;i++)
                              {
                                tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                              }
	                        }
                      } // end loop on the other tracks
                  }
                else if (Vtx_ntk > 2 && !success)
                  {
                    Vtx2_index.clear();
                  if (badtkhit_index == k) {vTT.erase(vTT.begin());ntracks--;Vtx2_index.push_back(p);}
                  else if (badtkhit_index == p) {vTT.erase(vTT.end());ntracks--;Vtx2_index.push_back(k);}
                  else {Vtx2_index.push_back(k);Vtx2_index.push_back(p);}
                  for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                    {
                      if (m == k || m == p || m == badtkhit_index) continue;
                      ntracks++;
                      tempMeanWeight=0;
                      
                      vTT.push_back(displacedTracks_Hemi2_mva[m]);
                      
                      TransientVertex updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                      if ( !updatedTV.isValid() ) 
	                      {  
	    	                  vTT.pop_back();
                          ntracks--;
                          if (vTT.size()<2) continue;
	    	                  updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                                                            tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx2_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  if (ntracks<2){success = false;}
                                  if (ntracks>=2){success = true;}
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                          continue;
	                      }
                      if ( updatedTV.isValid() ) 
	                      {
                          tempchi2 = updatedTV.normalisedChiSquared();
	    	                  tempx=updatedTV.position().x();
	    	                  tempy=updatedTV.position().y();
	    	                  tempz=updatedTV.position().z();
                          posError = updatedTV.positionError();
                          float TPosFH = sqrt((Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())*(Track_FirstHit_Hemi2_mva[m].second.X()-PV.x())+(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())*(Track_FirstHit_Hemi2_mva[m].second.Y()-PV.y())+(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z())*(Track_FirstHit_Hemi2_mva[m].second.Z()-PV.z()));
                          float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                          if (TPosFH>TPosVtx1 || Track_FirstHit_Hemi2_mva[m].first == true) {success = true;}//continue not useful
                          else  
                            {
                              vTT.pop_back();
                              ntracks--;
                              if (vTT.size()<2) continue;
	    	                      updatedTV = theFitter_Vertex_Hemi2_mva.vertex(vTT);
                                                                tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx2_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  if (ntracks<2){success = false;}
                                  if (ntracks>=2){success = true;}
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	    continue;
                            }
                            Vtx2_Weights.clear();
                          Vtx2_index.push_back(m);
                          success = true;
                          for(int i = 0; i<ntracks;i++)
                            {
                              tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                              Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                              if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                            }
	                      }
                    }
                }


                    Vtx_ntk = ntracks;
                    Vtx_chi = tempchi2;
                    Vtx_x = tempx;
                    Vtx_y = tempy;
                    Vtx_z = tempz;
                    Vtx_step = 2;
                    GlobalPoint RECOvtxPos(Vtx_x, Vtx_y, Vtx_z);
                    float DCA_VTX_Meand = 0;
                    for (int k = 0; k< Vtx_ntk; k++)
                      {
                        TrajectoryStateClosestToPoint DCA_Vtx = vTT[k].trajectoryStateClosestToPoint(RECOvtxPos);
                        if ( DCA_Vtx.isValid() ) // Be careful, all tracks are considered when looking at the DCA,
                        //but one should look of the wieghts of the track at the same time
                          {// The positions are given in the Global frame
                            float pca_Vtx_x = DCA_Vtx.position().x();
                            float pca_Vtx_y = DCA_Vtx.position().y();
                            float pca_Vtx_z = DCA_Vtx.position().z();
                            float refPoint_x = DCA_Vtx.referencePoint().x();
                            float refPoint_y = DCA_Vtx.referencePoint().y();
                            float refPoint_z = DCA_Vtx.referencePoint().z();
                            float DCA_Vtx_x = refPoint_x-pca_Vtx_x;
                            float DCA_Vtx_y = refPoint_y-pca_Vtx_y;
                            float DCA_Vtx_z = refPoint_z-pca_Vtx_z;
                            float DCA_VTX_r = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y);
                            float DCA_VTX_d = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y+DCA_Vtx_z*DCA_Vtx_z);
                            DCA_VTX_Meand+=DCA_VTX_d;
                            tree_Hemi_Vtx_track_DCA_x.push_back(DCA_Vtx_x);
                            tree_Hemi_Vtx_track_DCA_y.push_back(DCA_Vtx_y);
                            tree_Hemi_Vtx_track_DCA_z.push_back(DCA_Vtx_z);
                            tree_Hemi_Vtx_track_DCA_r.push_back(DCA_VTX_r);
                            tree_Hemi_Vtx_track_DCA_d.push_back(DCA_VTX_d);
                          }
                      }
                    if (MeanWeight==0)//<=>only two tracks in the valid vertex
                      {
                        Vtx2_Weights.clear();
                        for(int i = 0; i<ntracks;i++)
                          {
                            MeanWeight+=TV.trackWeight(vTT[i]);
                            Vtx2_Weights.push_back(TV.trackWeight(vTT[i]));
                          }
                      }
                    DCA_VTX_Meand = DCA_VTX_Meand/(float)(Vtx_ntk);
                    // tree_Hemi_Vtx_track_MeanDCA_d.push_back(DCA_VTX_Meand);
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
if (showlog){std::cout<<"success Hemi2 step 2 : "<<success<<std::endl;}
      }
    


    TransientVertex displacedVertex_step2_Hemi2;
    static AdaptiveVertexFitter theFitter_Vertex_step2_Hemi2(
    	       GeometricAnnealing ( sigmacut, Tini, ratio ), 
    	       DefaultLinearizationPointFinder(),
    	       KalmanVertexUpdator<5>(), 
    	       KalmanVertexTrackCompatibilityEstimator<5>(), 
    	       KalmanVertexSmoother() );
    theFitter_Vertex_step2_Hemi2.setParameters ( maxshift, maxlpshift, maxstep, weightThreshold );
    badVtx = false;


    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 2 step 3-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

  if ( (Vtx_chi < 0. || Vtx_chi > 10.) && ActivateStep2 ) badVtx = true;
  if ( badVtx && displacedTracks_step2_Hemi2.size() > 1)
    {
      DCA_VTX_Meand = 0;
      badtkhit_index = -1;
      bool success = false;
      MeanWeight=0;
      tempMeanWeight=0;
      std::vector<TransientTrack> vTT;
      Vtx_ntk = displacedTracks_step2_Hemi2.size();
      Vtx_chi = -10.;
      for (int k = 0 ; k <Vtx_ntk-1;k++)
        {
          for (int p = k+1 ; p < Vtx_ntk ; p++)
            {
              vTT.push_back(displacedTracks_step2_Hemi2[p]);
              vTT.push_back(displacedTracks_step2_Hemi2[k]);
              
              ntracks = 2;
              TransientVertex TV = theFitter_Vertex_step2_Hemi2.vertex(vTT); // We take the first "good-looking" seed to start
              if ( TV.isValid())
                {
                  for (int m = 0; m < ntracks; m++) // We then add track by track to the vertex and check the validity of the vertex
                    {
                      if(Track_FirstHit_step2_Hemi2[m].first == true) continue; //first hit of lost track is biaised
                      float PosFH = sqrt((Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z()));
                      float PosVtx1 = sqrt((TV.position().x()-PV.x())*(TV.position().x()-PV.x())+(TV.position().y()-PV.y())*(TV.position().y()-PV.y())+(TV.position().z()-PV.z())*(TV.position().z()-PV.z()));
                      if (PosFH>PosVtx1) 
                        {
                          success = true; 
                          tempchi2 = TV.normalisedChiSquared();
	    	                  tempx=TV.position().x();
	    	                  tempy=TV.position().y();
	    	                  tempz=TV.position().z();
                          posError = TV.positionError();
                          continue;
                        }//continue not useful
                      else{badtkhit_index = m;}
                    }

                    if (Vtx_ntk == 2 && !success){break;}// removing 1 track gives no other option
                    else if (success)
                      {
                        Vtx2_index.clear();
                          Vtx2_index.push_back(k);
                          Vtx2_index.push_back(p);
                        for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                          {
                            if (m == k || m == p) continue;
                            ntracks++;
                            vTT.push_back(displacedTracks_step2_Hemi2[m]);
                            TransientVertex updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                            if ( !updatedTV.isValid() ) 
	                            {  
	    	                        vTT.pop_back();
                                ntracks--;
	    	                        updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                                tempchi2 = updatedTV.normalisedChiSquared();
	    	                        tempx=updatedTV.position().x();
	    	                        tempy=updatedTV.position().y();
	    	                        tempz=updatedTV.position().z();
                                posError = updatedTV.positionError();
                                Vtx2_Weights.clear();
                                tempMeanWeight=0;
                                
                                for(int i = 0; i<ntracks;i++)
                                  {
                                    tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                    Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                    if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                  }
                        	      continue;
	                            }
                            if ( updatedTV.isValid() ) 
	                            {
                                tempchi2 = updatedTV.normalisedChiSquared();
	    	                        tempx=updatedTV.position().x();
	    	                        tempy=updatedTV.position().y();
	    	                        tempz=updatedTV.position().z();
                                posError = updatedTV.positionError();
                                float TPosFH = sqrt((Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z()));
                                float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                                if (TPosFH>TPosVtx1 || Track_FirstHit_step2_Hemi2[m].first == true) {success = true;}//continue not useful
                                else  
                                  {
                                    vTT.pop_back();
                                    ntracks--;
	    	                            updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                                    tempchi2 = updatedTV.normalisedChiSquared();
	    	                            tempx=updatedTV.position().x();
	    	                            tempy=updatedTV.position().y();
	    	                            tempz=updatedTV.position().z();
                                    posError = updatedTV.positionError();
                                    Vtx2_Weights.clear();
                                    tempMeanWeight=0;
                                    
                                    if (ntracks<2){success = false;}
                                    if (ntracks>=2){success = true;}
                                    for(int i = 0; i<ntracks;i++)
                                      {
                                        tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                        Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                        if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                      }
                        	          continue;
                                  }
                                tempMeanWeight=0;
                                Vtx2_Weights.clear();
                                Vtx2_index.push_back(m);
                                for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
	                            }
                          }
                      }
                    else if (Vtx_ntk > 2 && !success)
                      {
                        Vtx2_index.clear();
                        if (badtkhit_index == k) {vTT.erase(vTT.begin());ntracks--;Vtx2_index.push_back(p);}
                        else if (badtkhit_index == p) {vTT.erase(vTT.end());ntracks--;Vtx2_index.push_back(k);}
                        else {Vtx2_index.push_back(k);Vtx2_index.push_back(p);}
                        for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                          {
                            if (m == k || m == p || m == badtkhit_index) continue;
                            ntracks++;
                            tempMeanWeight=0;
                            
                            vTT.push_back(displacedTracks_step2_Hemi2[m]);
                            
                            TransientVertex updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                            if ( !updatedTV.isValid() ) 
	                            {  
	    	                        vTT.pop_back();
                                ntracks--;
                                if (vTT.size()<2) continue;
	    	                        updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                                tempchi2 = updatedTV.normalisedChiSquared();
	    	                        tempx=updatedTV.position().x();
	    	                        tempy=updatedTV.position().y();
	    	                        tempz=updatedTV.position().z();
                                posError = updatedTV.positionError();
                                Vtx2_Weights.clear();
                                tempMeanWeight=0;
                                
                                if (ntracks<2){success = false;}
                                if (ntracks>=2){success = true;}
                                for(int i = 0; i<ntracks;i++)
                                  {
                                    tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                    Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                    if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                  }
                        	      continue;
	                            }
                            if ( updatedTV.isValid() ) 
	                            {
                                tempchi2 = updatedTV.normalisedChiSquared();
	    	                        tempx=updatedTV.position().x();
	    	                        tempy=updatedTV.position().y();
	    	                        tempz=updatedTV.position().z();
                                posError = updatedTV.positionError();
                                float TPosFH = sqrt((Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z()));
                                float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                                if (TPosFH>TPosVtx1 || Track_FirstHit_step2_Hemi2[m].first == true) {success = true;}//continue not useful
                                else  
                                  {
                                    vTT.pop_back();
                                    ntracks--;
                                    if (vTT.size()<2) continue;
	    	                            updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                                    tempchi2 = updatedTV.normalisedChiSquared();
	    	                            tempx=updatedTV.position().x();
	    	                            tempy=updatedTV.position().y();
	    	                            tempz=updatedTV.position().z();
                                    posError = updatedTV.positionError();
                                    Vtx2_Weights.clear();
                                    tempMeanWeight=0;
                                    
                                    if (ntracks<2){success = false;}
                                    if (ntracks>=2){success = true;}
                                    for(int i = 0; i<ntracks;i++)
                                      {
                                        tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                        Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                        if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                      }
                        	          continue;
                                  }
                                  Vtx2_Weights.clear();
                                  Vtx2_index.push_back(m);
                                  success = true;
                                for(int i = 0; i<ntracks;i++)
                                  {
                                    tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                    Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                    if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                  }
	                            }
                          }
                      }

                        // We should have a Vertex after these conditions
                        Vtx_ntk = ntracks;
                        Vtx_chi = tempchi2;
                        Vtx_x = tempx;
                        Vtx_y = tempy;
                        Vtx_z = tempz;
                        Vtx_step = 3;
                        GlobalPoint RECOvtxPos(Vtx_x, Vtx_y, Vtx_z);
                        float DCA_VTX_Meand = 0;
                        for (int k = 0; k< Vtx_ntk; k++)
                          {
                            // MeanWeight+=displacedVertex_step2_Hemi2.trackWeight(displacedTracks_step2_Hemi2[p]);
                            // Vtx1_Weights.push_back(displacedVertex_step2_Hemi2.trackWeight(displacedTracks_step2_Hemi2[p]));
                            // tree_Hemi_Vtx_trackWeight.push_back(displacedVertex_step2_Hemi2.trackWeight(displacedTracks_step2_Hemi2[p]));
                            TrajectoryStateClosestToPoint DCA_Vtx = vTT[k].trajectoryStateClosestToPoint(RECOvtxPos);
                            if ( DCA_Vtx.isValid() ) // Be careful, all tracks are considered when looking at the DCA,
                            //but one could look of the wieghts of the track at the same time
                              {// The positions are given in the Global frame
                                float pca_Vtx_x = DCA_Vtx.position().x();
                                float pca_Vtx_y = DCA_Vtx.position().y();
                                float pca_Vtx_z = DCA_Vtx.position().z();
                                float refPoint_x = DCA_Vtx.referencePoint().x();
                                float refPoint_y = DCA_Vtx.referencePoint().y();
                                float refPoint_z = DCA_Vtx.referencePoint().z();
                                float DCA_Vtx_x = refPoint_x-pca_Vtx_x;
                                float DCA_Vtx_y = refPoint_y-pca_Vtx_y;
                                float DCA_Vtx_z = refPoint_z-pca_Vtx_z;
                                float DCA_VTX_r = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y);
                                float DCA_VTX_d = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y+DCA_Vtx_z*DCA_Vtx_z);
                                DCA_VTX_Meand+=DCA_VTX_d;
                                tree_Hemi_Vtx_track_DCA_x.push_back(DCA_Vtx_x);
                                tree_Hemi_Vtx_track_DCA_y.push_back(DCA_Vtx_y);
                                tree_Hemi_Vtx_track_DCA_z.push_back(DCA_Vtx_z);
                                tree_Hemi_Vtx_track_DCA_r.push_back(DCA_VTX_r);
                                tree_Hemi_Vtx_track_DCA_d.push_back(DCA_VTX_d);
                      
                              }
                          }
                        DCA_VTX_Meand = DCA_VTX_Meand/(float)(Vtx_ntk);
                        if (MeanWeight==0)//<=>only two tracks in the valid vertex
                          {
                            Vtx2_Weights.clear();
                            for(int i = 0; i<ntracks;i++)
                              {
                                MeanWeight+=TV.trackWeight(vTT[i]);
                                Vtx2_Weights.push_back(TV.trackWeight(vTT[i]));
                              }
                          }
                }
              else
                {
                  ntracks=0;
                  vTT.clear();
                }
              if ( success ) break;
            }
          if(success){break;}
        }
        if (showlog){std::cout<<"success Hemi2 step 3 : "<<success<<std::endl;}
    }



    //------------------------------------------------------------------------------------------------//
    //------------------------------------------Hemi 2 step 4-----------------------------------------//
    //------------------------------------------------------------------------------------------------//

    ntracks   = -2;
    tempchi2  = -10.;
    tempx = -100.;
    tempy = -100.;
    tempz = -100.;
    badVtx = false;
    if ( (Vtx_chi < 0. || Vtx_chi > 10.)  ) badVtx = true;
    if ( badVtx && IterAVF && displacedTracks_step2_Hemi2.size() > 1 && (ActivateStep4 || IterAVF) ) 
      {
        MeanWeight=0;
        bool success = false;
        std::vector<TransientTrack> vTT;
        tempchi2 = -10.;
        DCA_VTX_Meand = 0;
        tempMeanWeight=0;
        Vtx_ntk = displacedTracks_step2_Hemi2.size();
        for ( int p = 1; p < Vtx_ntk; p++ )
          {
            for  ( int k = 0; k < p; k++ ) // take pairs of tracks of highest BDT value
              {
                vTT.push_back(displacedTracks_step2_Hemi2[p]);
                vTT.push_back(displacedTracks_step2_Hemi2[k]);
                ntracks = 2;
                
                TransientVertex TV = theFitter_Vertex_step2_Hemi2.vertex(vTT); // We take the first "good-looking" seed to start
                if ( TV.isValid() && TV.normalisedChiSquared()>0 && TV.normalisedChiSquared()<10 ) 
                  {
	                  tempchi2 = TV.normalisedChiSquared();
	                  tempx = TV.position().x();
	                  tempy = TV.position().y();
	                  tempz = TV.position().z();
                    posError =TV.positionError();
                    success = true;

                    for (int m = 0; m < ntracks; m++) // We then add track by track to the vertex and check the validity of the vertex
                      {
                        if(Track_FirstHit_step2_Hemi2[m].first == true) continue; //first hit of lost track is biaised
                        float PosFH = sqrt((Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z()));
                        float PosVtx1 = sqrt((TV.position().x()-PV.x())*(TV.position().x()-PV.x())+(TV.position().y()-PV.y())*(TV.position().y()-PV.y())+(TV.position().z()-PV.z())*(TV.position().z()-PV.z()));
                        if (PosFH>PosVtx1) 
                          {
                            success = true; 
                            tempchi2 = TV.normalisedChiSquared();
	    	                    tempx=TV.position().x();
	    	                    tempy=TV.position().y();
	    	                    tempz=TV.position().z();
                            posError = TV.positionError();
                            continue;
                          }//continue not useful
                        else{badtkhit_index = m;}
                      }

                if (Vtx_ntk == 2 && !success){break;}// removing 1 track gives no other option
                else if (success)
                  {
                    Vtx2_index.clear();
                    Vtx2_index.push_back(k);
                    Vtx2_index.push_back(p);
                    for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                      {
                        if (m == k || m == p) continue;
                        ntracks++;
                        tempMeanWeight=0;
                        
                        vTT.push_back(displacedTracks_step2_Hemi2[m]);
                        TransientVertex updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                        if ( !updatedTV.isValid() ) 
	                        {  
	    	                    vTT.pop_back();
	    	                    ntracks--;
	    	                    updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
	    	                    tempchi2 = updatedTV.normalisedChiSquared();
	    	                    tempx=updatedTV.position().x();
	    	                    tempy=updatedTV.position().y();
        	                  tempz=updatedTV.position().z();
                            posError = updatedTV.positionError();
                            Vtx2_Weights.clear();
                            tempMeanWeight=0;
                            
                           for(int i = 0; i<ntracks;i++)
                              {
                                tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                              }
	    	                    continue;
	                        } 
                        if ( updatedTV.isValid() ) 
	                        {
                            tempchi2 = updatedTV.normalisedChiSquared();
	    	                    tempx=updatedTV.position().x();
	    	                    tempy=updatedTV.position().y();
	    	                    tempz=updatedTV.position().z();
                            posError = updatedTV.positionError();
                            float TPosFH = sqrt((Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z()));
                          float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                          if (TPosFH>TPosVtx1 || Track_FirstHit_step2_Hemi2[m].first == true) {success = true;}//continue not useful
                          else  
                            {
                              vTT.pop_back();
                              ntracks--;
	    	                      updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                                                                tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx2_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  if (ntracks<2){success = false;}
                                  if (ntracks>=2){success = true;}
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	    continue;
                            }
                          Vtx2_Weights.clear();
                          Vtx2_index.push_back(m);
                            for(int i = 0; i<ntracks;i++)
                              {
                                tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                              }
	                        }
                      } // end loop on the other tracks
                  }


                else if (Vtx_ntk > 2 && !success)
                  {
                    Vtx2_index.clear();
                    if (badtkhit_index == k) {vTT.erase(vTT.begin());ntracks--;Vtx2_index.push_back(p);}
                    else if (badtkhit_index == p) {vTT.erase(vTT.end());ntracks--;Vtx2_index.push_back(k);}
                    else {Vtx2_index.push_back(k);Vtx2_index.push_back(p);}
                  for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                    {
                      if (m == k || m == p || m == badtkhit_index) continue;
                      ntracks++;
                      tempMeanWeight=0;
                      
                      vTT.push_back(displacedTracks_step2_Hemi2[m]);
                      TransientVertex updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                      if ( !updatedTV.isValid() ) 
	                      {  
	    	                  vTT.pop_back();
                          ntracks--;
                          if (vTT.size()<2) continue;
	    	                  updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                                                            tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx2_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  if (ntracks<2){success = false;}
                                  if (ntracks>=2){success = true;}
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                          continue;
	                      }
                      if ( updatedTV.isValid() ) 
	                      {
                          tempchi2 = updatedTV.normalisedChiSquared();
	    	                  tempx=updatedTV.position().x();
	    	                  tempy=updatedTV.position().y();
	    	                  tempz=updatedTV.position().z();
                          posError = updatedTV.positionError();
                          
                          float TPosFH = sqrt((Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())*(Track_FirstHit_step2_Hemi2[m].second.X()-PV.x())+(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())*(Track_FirstHit_step2_Hemi2[m].second.Y()-PV.y())+(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z())*(Track_FirstHit_step2_Hemi2[m].second.Z()-PV.z()));
                          float TPosVtx1 = sqrt((tempx-PV.x())*(tempx-PV.x())+(tempy-PV.y())*(tempy-PV.y())+(tempz-PV.z())*(tempz-PV.z()));
                          if (TPosFH>TPosVtx1 || Track_FirstHit_step2_Hemi2[m].first == true) {success = true;}//continue not useful
                          else  
                            {
                              vTT.pop_back();
                              ntracks--;
                              if (vTT.size()<2) continue;
	    	                      updatedTV = theFitter_Vertex_step2_Hemi2.vertex(vTT);
                                  tempchi2 = updatedTV.normalisedChiSquared();
	    	                          tempx=updatedTV.position().x();
	    	                          tempy=updatedTV.position().y();
	    	                          tempz=updatedTV.position().z();
                                  posError = updatedTV.positionError();
                                  Vtx2_Weights.clear();
                                  tempMeanWeight=0;
                                  
                                  if (ntracks<2){success = false;}
                                  if (ntracks>=2){success = true;}
                                  for(int i = 0; i<ntracks;i++)
                                    {
                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                      Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                      if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                    }
                        	    continue;
                            }
                            Vtx2_Weights.clear();
                          Vtx2_index.push_back(m);
                          success = true;
                          for(int i = 0; i<ntracks;i++)
                            {
                              tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                              Vtx2_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                              if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                            }
	                      }
                    }
                  }

                    Vtx_ntk = ntracks;
                    Vtx_chi = tempchi2;
                    Vtx_x = tempx;
                    Vtx_y = tempy;
                    Vtx_z = tempz;
                    Vtx_step = 4;
                    GlobalPoint RECOvtxPos(Vtx_x, Vtx_y, Vtx_z);
                    DCA_VTX_Meand = 0;
                    for (int k = 0; k< Vtx_ntk; k++)
                      {
                        TrajectoryStateClosestToPoint DCA_Vtx = vTT[k].trajectoryStateClosestToPoint(RECOvtxPos);
                        if ( DCA_Vtx.isValid() ) // Be careful, all tracks are considered when looking at the DCA,
                        //but one should look of the wieghts of the track at the same time
                          {// The positions are given in the Global frame
                            float pca_Vtx_x = DCA_Vtx.position().x();
                            float pca_Vtx_y = DCA_Vtx.position().y();
                            float pca_Vtx_z = DCA_Vtx.position().z();
                            float refPoint_x = DCA_Vtx.referencePoint().x();
                            float refPoint_y = DCA_Vtx.referencePoint().y();
                            float refPoint_z = DCA_Vtx.referencePoint().z();
                            float DCA_Vtx_x = refPoint_x-pca_Vtx_x;
                            float DCA_Vtx_y = refPoint_y-pca_Vtx_y;
                            float DCA_Vtx_z = refPoint_z-pca_Vtx_z;
                            float DCA_VTX_r = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y);
                            float DCA_VTX_d = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y+DCA_Vtx_z*DCA_Vtx_z);
                            DCA_VTX_Meand+=DCA_VTX_d;
                            tree_Hemi_Vtx_track_DCA_x.push_back(DCA_Vtx_x);
                            tree_Hemi_Vtx_track_DCA_y.push_back(DCA_Vtx_y);
                            tree_Hemi_Vtx_track_DCA_z.push_back(DCA_Vtx_z);
                            tree_Hemi_Vtx_track_DCA_r.push_back(DCA_VTX_r);
                            tree_Hemi_Vtx_track_DCA_d.push_back(DCA_VTX_d);
                          }
                      }
                    if (MeanWeight==0)//<=>only two tracks in the valid vertex
                      {
                        Vtx2_Weights.clear();
                        for(int i = 0; i<ntracks;i++)
                          {
                            MeanWeight+=TV.trackWeight(vTT[i]);
                            Vtx2_Weights.push_back(TV.trackWeight(vTT[i]));
                          }
                      }
                    DCA_VTX_Meand = DCA_VTX_Meand/(float)(Vtx_ntk);
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
        if (showlog){std::cout<<"success Hemi2 step 4 : "<<success<<std::endl;}
      }


    float Vtx_chi2 = Vtx_chi;
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
    tree_Hemi_Vtx_xError.push_back(posError.cxx());
    tree_Hemi_Vtx_yError.push_back(posError.cyy());
    tree_Hemi_Vtx_zError.push_back(posError.czz());
    theta_Vtx = tan(sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y)/abs(Vtx_z)) ;
    eta_Vtx = -TMath::Log(tan(theta_Vtx/2));
    if (Vtx_z<0){eta_Vtx = -eta_Vtx;}
    tree_Hemi_Vtx_eta.push_back(eta_Vtx);
    recX = Vtx_x - tree_PV_x;
    recY = Vtx_y - tree_PV_y;
    recZ = Vtx_z - tree_PV_z;
    recD = TMath::Sqrt(recX*recX + recY*recY + recZ*recZ);
    tree_Hemi_Vtx_dist.push_back( recD );
    tree_Hemi_Vtx_MeantrackWeight.push_back(MeanWeight);
    tree_Hemi_Vtx_track_MeanDCA_d.push_back(DCA_VTX_Meand);
    float posx2 = Vtx_x;
    float posy2 = Vtx_y;
    float posz2 = Vtx_z;
    if(Vtx_step>0 && Vtx_chi<10 && Vtx_chi>0){nVertex++;}

    float dr_2V = -10.;//quantities that should be lower for the backgrounds/ compared to signal
    float dz_2V = -10.;
    float dd_2V = -10.;
    if(nVertex==2)
      {
        dr_2V = sqrt((posx2-posx1)*(posx2-posx1)+(posy2-posy1)*(posy2-posy1));
        dz_2V = sqrt((posz2-posz1)*(posz2-posz1));
        dd_2V = sqrt((posx2-posx1)*(posx2-posx1)+(posy2-posy1)*(posy2-posy1)+(posz2-posz1)*(posz2-posz1));
      }
    tree_Hemi_Vtx_nVtx.push_back(nVertex);
    // tree_Hemi_Vtx_nVtx.push_back(nVertex);
    tree_Hemi_Vtx_Vtx_dr.push_back(dr_2V);
    tree_Hemi_Vtx_Vtx_dr.push_back(dr_2V);
    tree_Hemi_Vtx_Vtx_dz.push_back(dz_2V);
    tree_Hemi_Vtx_Vtx_dz.push_back(dz_2V);
    tree_Hemi_Vtx_Vtx_dd.push_back(dd_2V);
    tree_Hemi_Vtx_Vtx_dd.push_back(dd_2V);

    // Vertex Analysis Step

        // 0.7264 : Tight
        // 0.2770 : Medium 
        // 0.0494 : loose
    // pass BtagGood_Hemi1 
    int BtagGood_Hemi2 = 0 ;
    temp_px = 0;
    temp_py = 0;
    temp_pz = 0 ;
    temp_e = 0 ;
TLorentzVector Total4Vector2(0,0,0,0);
if (Vtx2_index.size() != Vtx2_Weights.size()) {std::cout<<"size Vtx2_weights and Vtx2_index and ntracks and  chi and step : "<<Vtx2_Weights.size()<<" and "<<Vtx2_index.size()<<" and "<<Vtx_ntk<<" and "<<Vtx_chi<<" and "<<Vtx_step<<std::endl;}



    if(Vtx_step==1 || Vtx_step==2)
      {
        // if(Vtx_chi>0 && Vtx_chi<10)
        //   {
            for (unsigned int i = 0 ; i <TrackInfo_Hemi2_mva.size(); i++)
              {
                for (unsigned int j = 0 ; j < Vtx2_index.size(); j++)
                  {
                    if (i == Vtx2_index[j])
                      {
                        if (Vtx2_Weights[j]>0.5)
                          {
                            temp_px = TrackInfo_Hemi2_mva[i].second.Px();
                            temp_py = TrackInfo_Hemi2_mva[i].second.Py();
                            temp_pz = TrackInfo_Hemi2_mva[i].second.Pz();
                            temp_e  = TrackInfo_Hemi2_mva[i].second.E();
                            TLorentzVector TLorentzTrack(temp_px,temp_py,temp_pz,temp_e);
                            Total4Vector2 +=TLorentzTrack;
                          }
                      }
                  }
              }
          // }
      }

    else if (Vtx_step == 3 || Vtx_step==4)
      {

            for (unsigned int i = 0 ; i <TrackInfo_Hemi2_mva.size(); i++)
              {
                for (unsigned int j = 0 ; j < Vtx2_index.size(); j++)
                  {
                    if (i == Vtx2_index[j])
                      {
                        if (Vtx2_Weights[j]>0.5)
                          {
                            temp_px = TrackInfo_step2_Hemi2[i].second.Px();
                            temp_py = TrackInfo_step2_Hemi2[i].second.Py();
                            temp_pz = TrackInfo_step2_Hemi2[i].second.Pz();
                            temp_e  = TrackInfo_step2_Hemi2[i].second.E();
                            TLorentzVector TLorentzTrack(temp_px,temp_py,temp_pz,temp_e);
                            Total4Vector2 +=TLorentzTrack;
                          }
                      }
                  }
              }
      }
    tree_Hemi_Vtx_BTag.push_back(BtagGood_Hemi2 );
    tree_Hemi_Vtx_Mass.push_back(sqrt(Total4Vector2.Mag2()));


//------------------------------------------------//
    float Vtx2_ntk = Vtx_ntk;
    float Vtx2_chi = Vtx_chi;
    float Vtx2_step = Vtx_step;
    float Vtx2_r = sqrt(Vtx_x*Vtx_x+Vtx_y*Vtx_y);
    float Vtx2_z = Vtx_z;
    float Vtx2_MTW = MeanWeight;
    float Vtx2_Mass = sqrt(Total4Vector2.Mag2());
    float H2_Mass = sqrt(TLorentzAxis2.Mag2());
    float Vtx2_dist = recD;
    float Vtx2_MeanDCA = DCA_VTX_Meand;


    // -------------------------------------------//
  
    float ping2 = 0;
    bool ping_Hemi1 = false, ping_Hemi2 = false;
      if (!runOnData_)
        {

  if ( tree_nLLP > 0 ) {
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
    tree_Hemi_LLP_dR12.push_back(dRneuneu);
    tree_Hemi_LLP_dR12.push_back(dRneuneu);

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
  }
  else  {
        tree_Hemi_LLP_ping.push_back( false );
        tree_Hemi_LLP_ping.push_back( false );
        }
    
    tree_Hemi_dR12.push_back(dR_axis12);
    tree_Hemi_dR12.push_back(dR_axis12);
        }

//------------Duplicate for each hemisphere-----------//
    // some informations for tracks in their hemisphere
    int ntrk10_vtx_hemi1 = 0., ntrk10_vtx_hemi2 = 0.;
    int ntrk20_vtx_hemi1 = 0., ntrk20_vtx_hemi2 = 0.;
    int NisjetH = 0;
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
        if ( Vtx_chi < 10. && Vtx_chi>0 ) {
          bool isLostTrack = tree_track_lost[counter_track];
          if (isLostTrack && RemoveLostTrackFromVtxSelec ) {continue;}
          float x1 = tree_track_firstHit_x[counter_track] - tree_PV_x;
          float y1 = tree_track_firstHit_y[counter_track] - tree_PV_y;
          float z1 = tree_track_firstHit_z[counter_track] - tree_PV_z;
          float vtx_x = tree_Hemi_Vtx_x[hemi-1] - tree_PV_x;
          float vtx_y = tree_Hemi_Vtx_y[hemi-1] - tree_PV_y;
          float vtx_z = tree_Hemi_Vtx_z[hemi-1] - tree_PV_z;
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
float Vtx1_ntrk10 =0;
float Vtx1_ntrk20 =0;
        //--------------VTX1-------------------//
        if (Vtx1_chi>0 && Vtx1_chi<10)
          {
            Vtx1_ntrk10 = ntrk10_vtx_hemi1 ;
            Vtx1_ntrk20 = ntrk20_vtx_hemi1;
            mva_V_nTrks  =  Vtx1_ntk;
            mva_V_chi    =  Vtx1_chi;
            mva_V_step   =  Vtx1_step;
            // mva_V_r  =   Vtx1_r;
            // mva_V_z  =    Vtx1_z;
            mva_V_MTW    =  Vtx1_MTW;
            mva_V_Mass   =  Vtx1_Mass;
            mva_H_Mass   =  H1_Mass;

            mva_V_ntrk10 = Vtx1_ntrk10;
            mva_V_ntrk20 = Vtx1_ntrk20;
            mva_V_MeanDCA = Vtx1_MeanDCA;
            // mva_V_dist  = Vtx1_dist;
          }
    double Vtx1_bdtVal = -10;
    Vtx1_bdtVal = readerVtx->EvaluateMVA("BDTG");// values at -999 => thishappens if the Hemi_Mass and Hemi_Vtx_Mass are not definite
   // default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
    tree_Hemi_Vtx_MVAval.push_back(Vtx1_bdtVal);

      if ((Vtx1_step==1 || Vtx1_step==2) && Vtx1_chi>0 && Vtx1_chi<10)
    {
      Vtx1_ntrk10 = ntrk10_vtx_hemi2 ;
      Vtx1_ntrk20 = ntrk20_vtx_hemi2;

      mva_V_nTrks  =  Vtx1_ntk;
      mva_V_chi    =  Vtx1_chi;
      mva_V_step   =  Vtx1_step;
      // mva_V_r  =   Vtx1_r;
      // mva_V_z  =    Vtx1_z;
      mva_V_MTW    =  Vtx1_MTW;
      mva_V_Mass   =  Vtx1_Mass;
      mva_H_Mass   =  H1_Mass;

      mva_V_ntrk10 = Vtx1_ntrk10;
      mva_V_ntrk20 = Vtx1_ntrk20;
      mva_V_MeanDCA = Vtx1_MeanDCA;

    }
      double Vtx1_bdtVal_Step1 = -10; 
      Vtx1_bdtVal_Step1 = readerVtxStep1->EvaluateMVA("BDTG");
      tree_Hemi_Vtx_MVAval_Step1.push_back(Vtx1_bdtVal_Step1);
float Vtx2_ntrk10=0;
float Vtx2_ntrk20=0;
        //--------------VTX2-------------------//
        if (Vtx2_chi>0 && Vtx2_chi<10)
          {
            Vtx2_ntrk10 = ntrk10_vtx_hemi2 ;
            Vtx2_ntrk20 = ntrk20_vtx_hemi2;

            mva_V_nTrks  =  Vtx2_ntk;
            mva_V_chi    =  Vtx2_chi;
            mva_V_step   =  Vtx2_step;
            // mva_V_r  =   Vtx2_r;
            // mva_V_z  =    Vtx2_z;
            mva_V_MTW    =  Vtx2_MTW;
            mva_V_Mass   =  Vtx2_Mass;
            mva_H_Mass   =  H2_Mass;
            mva_V_ntrk10 = Vtx2_ntrk10;
            mva_V_ntrk20 = Vtx2_ntrk20;
            mva_V_MeanDCA = Vtx2_MeanDCA;
            // mva_V_dist  = Vtx1_dist;
          }
    double Vtx2_bdtVal = -10;
    Vtx2_bdtVal = readerVtx->EvaluateMVA("BDTG");// values at -999
   // default value = -10 (no -10 observed and -999 comes from EvaluateMVA)
    tree_Hemi_Vtx_MVAval.push_back(Vtx2_bdtVal);

  if ((Vtx2_step==1 || Vtx2_step==2) && Vtx2_chi>0 && Vtx2_chi<10)
    {
      Vtx2_ntrk10 = ntrk10_vtx_hemi2 ;
      Vtx2_ntrk20 = ntrk20_vtx_hemi2;

      mva_V_nTrks  =  Vtx2_ntk;
      mva_V_chi    =  Vtx2_chi;
      mva_V_step   =  Vtx2_step;
      // mva_V_r  =   Vtx2_r;
      // mva_V_z  =    Vtx2_z;
      mva_V_MTW    =  Vtx2_MTW;
      mva_V_Mass   =  Vtx2_Mass;
      mva_H_Mass   =  H2_Mass;

      mva_V_ntrk10 = Vtx2_ntrk10;
      mva_V_ntrk20 = Vtx2_ntrk20;
      mva_V_MeanDCA = Vtx2_MeanDCA;

    }

      double Vtx2_bdtVal_Step1 = -10; 
      Vtx2_bdtVal_Step1 = readerVtxStep1->EvaluateMVA("BDTG");
      tree_Hemi_Vtx_MVAval_Step1.push_back(Vtx2_bdtVal_Step1);
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
    
    tree_Evts_MVAval.clear();

    tree_allPV_i.clear();
    tree_allPV_x.clear();
    tree_allPV_y.clear();
    tree_allPV_z.clear();
    tree_allPV_ex.clear();
    tree_allPV_ey.clear();
    tree_allPV_ez.clear();
    tree_allPV_NChi2.clear();
    tree_allPV_ndf.clear();
    
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
    tree_SecInt_ntrk10.clear();
    tree_SecInt_ntrk20.clear();
    tree_SecInt_LLP.clear();
    tree_SecInt_LLP_dr.clear();
    tree_SecInt_LLP_dz.clear();
    tree_SecInt_LLP_dd.clear();

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

    tree_jet_E.clear();
    tree_jet_pt.clear();
    tree_jet_eta.clear();
    tree_jet_phi.clear();
    tree_jet_btag_DeepCSV.clear();
    tree_jet_btag_DeepJet.clear();
    tree_jet_leadingpt.clear();
    tree_jet_leadingpt2.clear();
    tree_jet_leadingMuon_dR.clear();
    tree_jet_leadingMuon2_dR.clear();
    tree_jet_jet_dR.clear();
    tree_jet_jet_dPhi.clear();
    tree_jet_jet_dEta.clear();
    tree_muon_jet_dRmin.clear();
    tree_muon_jet_dRmax.clear();
    
    tree_electron_pt.clear();
    tree_electron_eta.clear();
    tree_electron_phi.clear();
    tree_electron_x.clear();
    tree_electron_y.clear();
    tree_electron_z.clear();
    tree_electron_px.clear();
    tree_electron_py.clear();
    tree_electron_pz.clear();
    tree_electron_energy.clear();
    tree_electron_charge.clear();
    tree_electron_isoR4.clear();
    tree_electron_IsLoose.clear();
    tree_electron_IsMedium.clear();
    tree_electron_IsTight.clear();
    tree_electron_trigger_Ele.clear();
    tree_electron_trigger_diEle.clear();
    tree_electron_dxy.clear();
    tree_electron_dz.clear();

    tree_ST.clear();
    tree_muon_pt.clear();
    tree_muon_eta.clear();
    tree_muon_phi.clear();
    tree_muon_x.clear();
    tree_muon_y.clear();
    tree_muon_z.clear();
    tree_muon_px.clear();
    tree_muon_py.clear();
    tree_muon_pz.clear();
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
    tree_muon_trigger_dimu.clear();  
    tree_muon_trigger_isomu.clear();
    tree_muon_PFIsoLoose.clear();
    tree_muon_PFIsoMedium.clear();
    tree_muon_PFIsoTight.clear();
    tree_muon_TkIsoLoose.clear();
    tree_muon_TkIsoTight.clear();
    tree_muon_nmu.clear();
    tree_muon_leadingpt.clear();
    tree_muon_leadingpt2.clear();
    tree_muon_muon_dR.clear();
    tree_muon_muon_dPhi.clear();
    tree_muon_muon_dEta.clear();
//     tree_passesTrkPtr.clear();

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
    tree_track_dzTOpu.clear();
    tree_track_dzSigTOpu.clear();
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
    tree_track_btag.clear();
    tree_track_energy.clear();

    tree_track_Hemi.clear();
    tree_track_Hemi_dR.clear();
    tree_track_Hemi_dRmax.clear();
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

    // tree_genFromb_x.clear();
    // tree_genFromb_y.clear();
    // tree_genFromb_z.clear();
    // tree_genFromb_pdgid.clear();
    // tree_genFromb_LLP_dV.clear();
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
    tree_Hemi_Mass.clear();
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
    tree_Hemi_Vtx_xError.clear();
    tree_Hemi_Vtx_yError.clear();
    tree_Hemi_Vtx_zError.clear();
    tree_Hemi_Vtx_eta.clear();
    tree_Hemi_Vtx_Vtx_dr.clear();
    tree_Hemi_Vtx_Vtx_dz.clear();
    tree_Hemi_Vtx_Vtx_dd.clear();
    tree_Hemi_Vtx_nVtx.clear();
    tree_Hemi_Vtx_trackWeight.clear();
    tree_Hemi_Vtx_MeantrackWeight.clear();
    tree_Hemi_Vtx_track_DCA_x.clear();
    tree_Hemi_Vtx_track_DCA_y.clear();
    tree_Hemi_Vtx_track_DCA_z.clear();
    tree_Hemi_Vtx_track_DCA_r.clear();
    tree_Hemi_Vtx_track_DCA_d.clear();
    tree_Hemi_Vtx_track_MeanDCA_d.clear();
    tree_Hemi_Vtx_BTag.clear();
    tree_Hemi_Vtx_Mass.clear();
    tree_Hemi_Vtx_MVAval.clear();
    tree_Hemi_Vtx_MVAval_Step1.clear();
    tree_Hemi_Vtx_TVtx_dx.clear();
    tree_Hemi_Vtx_TVtx_dy.clear();
    tree_Hemi_Vtx_TVtx_dz.clear();
    tree_Hemi_Vtx_TVtx_NChi2.clear();

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
    tree_track_Hemi_isjet.clear();
    tree_Hemi_Vtx_ddToBkg.clear();
    tree_Hemi_LLP_ping.clear();
    tree_event_LLP_ping.clear();

    tree_Hemi_LooseBTag_axes.clear();
    tree_Hemi_MediumBTag_axes.clear();
    tree_Hemi_TightBTag_axes.clear();

}
