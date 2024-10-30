import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
process = cms.Process("FlyingTop",Run3_2023)

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#Prompt_Reco


process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
##----------------------paul--------------------------##
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")


## JeC JER for systematics ###########################
process.load("JetMETCorrections.Configuration.JetCorrectors_cff")
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
from CondCore.CondDB.CondDB_cfi import *


# https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_dat 
###--------------------------
IsMC = False
#$$
year = 2023
#$$

from Configuration.AlCa.GlobalTag import GlobalTag
# Global Tags:

isPost = True
ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2018UL.txt"
EGERA = '2022-Prompt'
GT = '140X_dataRun3_v3' 
TIGHTJETIDERA = 'RUN3CHSruns2022FGruns2023CD'  #NOTE: RUN3CHSrunsBCDEprompt, RUN3CHSruns2022FGruns2023CD, RUN2ULCHS
L1PREFERA = '20172018'
DATAPUFILE = 'PU_Run2023_data.root'
DATAPUFILEUP = 'pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-72400ub-100bins.root'
DATAPUFILEDOWN = 'pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-66000ub-100bins.root'
MCPUFILE   = 'PU_Run2023_BPix_MC.root'
if isPost:
    MCPUFILE   = 'PU_Run2023_MC.root'   

if IsMC:                                                                                                                                                                                     
    process.GlobalTag = GlobalTag(process.GlobalTag, GT, '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, GT, '') 


#$$
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
#$$

process.GoodVertexFilter = cms.EDFilter("VertexSelector",
                                        src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                                        filter = cms.bool(True)   # otherwise it won't filter the events, just produce an empty vertex collection.                                        
                                       )

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
  # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/TTTo2L2Nu_2024.root'
#$$
  '/store/data/Run2023C/MuonEG/MINIAOD/PromptReco-v1/000/367/112/00000/7d9bff86-5cd0-4b11-b42e-8f7a7949a028.root'
#$$
)
)


##########################################################
#
# Setup AK4 jets to be used for analysis
#
##########################################################
#
# Setup JEC factors
#
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
process.jetCorrFactors = patJetCorrFactors.clone(src='slimmedJetsPuppi',
    levels = cms.vstring('L1FastJet',
        'L2Relative',
        'L3Absolute',
        'L2L3Residual'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
)
#
# This module will take the JEC factors and update them on slimmedJets
#
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import *
process.updatedJets = updatedPatJets.clone(
    addBTagInfo=False,
    jetSource='slimmedJetsPuppi',
    jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactors") ),
)
#
#  Module to calculate JetID for "Tight" working point
#
process.tightJetId = cms.EDProducer("PatJetIDValueMapProducer",
    filterParams=cms.PSet(
        version = cms.string(TIGHTJETIDERA),
        quality = cms.string('TIGHT'),
    ),
    src = cms.InputTag("updatedJets")
)

process.tightLepVetoJetId = cms.EDProducer("PatJetIDValueMapProducer",
  filterParams=cms.PSet(
    version = cms.string(TIGHTJETIDERA),
    quality = cms.string('TIGHTLEPVETO'),
  ),
  src = cms.InputTag("updatedJets")
)
#
# Module to calculate Pileup Jet ID
# Might not be useful for RUn 3 since we have PUPPI jets
# !!  from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_140X_23
# process.load("RecoJets.JetProducers.PileupJetID_cfi")
# process.pileupJetIdUpdated = process.pileupJetId.clone(
#     jets=cms.InputTag("updatedJets"),# JEC corrected jets
#     inputIsCorrected=True,
#     applyJec=False,
#     vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
#     algos = cms.VPSet(_chsalgos_140X_23),
# )
#
# Embed the Jet ID and Pileup Jet ID variables in the jets.
#
process.updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
    src = cms.InputTag("updatedJets"),
    # userFloats = cms.PSet(
    #     puIdDisc = cms.InputTag('pileupJetIdUpdated:fullDiscriminant'),
    # ),
    userInts = cms.PSet(
        # puId = cms.InputTag('pileupJetIdUpdated:fullId'),
        tightId = cms.InputTag("tightJetId"),
        tightLepVetoId = cms.InputTag("tightLepVetoJetId"),
    ),
)

################################################""
# process.jer = cms.ESSource("PoolDBESSource",
#         CondDBSetup,
#         toGet = cms.VPSet(
#             # Resolution
#             cms.PSet(
#                 record = cms.string('JetResolutionRcd'),
#                 # JR_dataRun2_25nsV1b_V3b_V7b_106X_DATA_SF_AK4PF
#                 # JR_dataRun2_25nsV1b_V3b_V7b_106X_DATA_PtResolution_AK4PF
#                 # JR_Summer15_25nsV6_MC_PtResolution_AK4PF
#                 tag    = cms.string('JR_Summer15_25nsV6_MC_PtResolution_AK4PFchs'),
#                 label  = cms.untracked.string('AK4PFchs_pt')
#                 ),

#             # Scale factors
#             cms.PSet(
#                 record = cms.string('JetResolutionScaleFactorRcd'),
#                 tag    = cms.string('JR_Summer15_25nsV6_MC_SF_AK4PFchs'),
#                 label  = cms.untracked.string('AK4PFchs')
#                 ),
#             ),
#         connect = cms.string('sqlite:Summer16_25nsV1b_DATA.db')
#         )
################################################""

# from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
# process.prefiringweight = l1PrefiringWeightProducer.clone(
#     #TheJets = cms.InputTag("updatedPatJetsUpdatedJEC"), #this should be the slimmedJets collection with up to date JECs !                                                               
#     TheJets= cms.InputTag('slimmedJets'),
#     DataEraECAL = cms.string("None"),
#     DataEraMuon = cms.string(L1PREFERA),
#     UseJetEMPt = cms.bool(False),
#     PrefiringRateSystematicUnctyECAL = cms.double(0.2),
#     PrefiringRateSystematicUnctyMuon = cms.double(0.2)
# )


process.options = cms.untracked.PSet( )
process.FlyingTop = cms.EDAnalyzer("FlyingTopAnalyzer",
       #    RochString = cms.string("./FlyingTop/FlyingTop/test/"), 
       #    Roccor = cms.FileInPath("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"),
           DATASET = cms.untracked.vstring(process.source.fileNames),
           isMC =cms.bool(IsMC), 
           YEAR = cms.int32(year),
           ERA2016 = cms.bool(isPost),
           RochString = cms.string(ROCCORPATH),
           weightFileMVA = cms.untracked.string( "BDT_TRK_240510_ctau100vsEMUdata_NOchi2NOdxyNOdz.xml"),#track selection => previous :BDT_TRK_CTAU10cm_vs_TT.xml //BDT_TRK_ALLSIGvsDYTT_30_01_2024.xml/ BDT_TRK_ALLSIGvsALLBKG.xml // TMVAClassification_BDTG_TRKSEL_.weights.xml
           weightFileMVA_EVTS = cms.untracked.string("BDT_EVT_ALLSIGvsALLBKG.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
           weightFileMVA_EVTSDY = cms.untracked.string("BDT_EVT_ALLSIGvsDYM50.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
           weightFileMVA_EVTSTT = cms.untracked.string("BDT_EVT_ALLSIGvsTTTo2L2Nu.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
           weightFileMVA_VTX = cms.untracked.string("BDT_VTX_ALLSTEPS.xml"),#vtx selection : TMVAClassification_BDTG_VTXSEL_.weights.xml
           weightFileMVA_VTX_step1 = cms.untracked.string("BDT_VTX_STEP12.xml"),#vtx selection :  TMVAClassification_BDTG_VTXSel_TIGHTWP.weights.xml
           mcpufile = cms.string(MCPUFILE),
           mcpupath = cms.string("pileup"),
           datapufile = cms.string(DATAPUFILE),
           datapileupfileup = cms.string(DATAPUFILEUP),
           datapileupfiledown = cms.string(DATAPUFILEDOWN),
           datapupath = cms.string("pileup"),
       #    muoneps1file   = cms.string("NUM_TightID_DEN_TrackerMuons_abseta_pt.root"),
       #    muoneps1path   = cms.string("NUM_TightID_DEN_TrackerMuons_abseta_pt"),
       #    muoneps2file   = cms.string("NUM_TkIsoTight_DEN_TightID_TrackerMuons_abseta_pt.root"),
       #    muoneps2path   = cms.string("NUM_TkIsoTight_DEN_TightID_TrackerMuons_abseta_pt"),
       #    muoneps3file   = cms.string("NUM_IsoMu24_or_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_DEN_TkIsoTight_and_TightID_abseta_pt.root"),
       #    muoneps3path   = cms.string("NUM_IsoMu24_or_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_DEN_TkIsoTight_and_TightID_abseta_pt"),
           genEventInfoInput	    = cms.InputTag("generator"),
           LHEEventProductInput     = cms.InputTag("externalLHEProducer"),#source or externalLHEProducer
           genpruned   = cms.InputTag('prunedGenParticles'),
           genpacked   = cms.InputTag('packedGenParticles'),
           vertices    = cms.InputTag('offlineSlimmedPrimaryVertices'),
           mets        = cms.InputTag("slimmedMETs"),
           jets        = cms.InputTag("updatedJetsWithUserData"),
           genjets     = cms.InputTag("slimmedGenJets"),
           electrons   = cms.InputTag("slimmedElectrons"),
           muons       = cms.InputTag("slimmedMuons"),
           pfCands     = cms.InputTag("packedPFCandidates"),
           lostpfCands = cms.InputTag("lostTracks"),
           Kshorts     = cms.InputTag("slimmedKshortVertices"),#recoVertexCompositePtrCandidates_slimmedKshortVertices__PAT
           Lambda      = cms.InputTag("slimmedLambdaVertices"),#recoVertexCompositePtrCandidates_slimmedLambdaVertices__PAT
           beamSpot    = cms.untracked.InputTag('offlineBeamSpot'),
           puCollection = cms.InputTag("slimmedAddPileupInfo"),
           rhoCollection = cms.InputTag("fixedGridRhoFastjetAll"),
       )

#process.tsk = cms.Task()
#for mod in process.producers_().itervalues():
#    process.tsk.add(mod)
#    for mod in process.filters_().itervalues():
#        process.tsk.add(mod)
#-------------------------------------
#process.p = cms.Path(process.prefiringweight* process.egammaPostRecoSeq* process.updatedPatJetsTransientCorrectedNewDFTraining* process.FlyingTop,process.tsk)
process.p = cms.Path(
    process.GoodVertexFilter*
    process.jetCorrFactors*
    process.updatedJets*
    process.tightJetId*
    process.tightLepVetoJetId *
    # !!  process.pileupJetIdUpdated* # Using Puppi jets so not needed
    process.updatedJetsWithUserData*
    # !! process.prefiringweight* # I have to figure out why it is not working for run 3
    # !! process.egammaPostRecoSeq* #I have to figure out why it is not working for run 3
    process.FlyingTop
)
# //jet energy corrections
from PhysicsTools.PatAlgos.tools.helpers  import getPatAlgosToolsTask
process.patAlgosToolsTask = getPatAlgosToolsTask(process)
process.pathRunPatAlgos = cms.Path(process.patAlgosToolsTask)
########## output of ntuple

process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple.root") )

#$$
process.options.numberOfThreads=cms.untracked.uint32(4)
#$$

