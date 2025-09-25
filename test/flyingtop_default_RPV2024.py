import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2024_cff import Run3_2024
process = cms.Process("FlyingTop",Run3_2024)

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
#$$
IsMC = True
# Year = 2024
# isPost = True
#$$

from Configuration.AlCa.GlobalTag import GlobalTag

ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2018UL.txt"
EGERA = '2022-Prompt'
L1PREFERA = '20172018'

GT = '140X_mcRun3_2024_realistic_v26'
TIGHTJETIDERA = 'RUN3CHSruns2022FGruns2023CD'
MCPUFILE   = 'Run2024_MC.root'
DATAPUFILE = 'Run2024_DATA.root'
DATAPUFILEUP   = 'pileupHistogram-Cert_Collisions2024_378891_386951_GoldenJson-13p6TeV-72400ub-100bins.root'
DATAPUFILEDOWN = 'pileupHistogram-Cert_Collisions2024_378891_386951_GoldenJson-13p6TeV-66000ub-100bins.root'
JECUNCDATA = 'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt' 
JECUNCMC =   'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt'
JERMC =      'Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt' 
JERSFMC =    'Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFPuppi.txt'
JERDATA =    'Summer23BPixPrompt23_RunD_JRV1_DATA_PtResolution_AK4PFchs.txt' 
JERSFDATA =  'Summer23BPixPrompt23_RunD_JRV1_DATA_SF_AK4PFPuppi.txt' 
TRKBDT = "BDT_TRK_250527_2024_ctau100vsEMUdata.xml"

process.GlobalTag = GlobalTag(process.GlobalTag, GT, '')

#$$
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#$$

process.GoodVertexFilter = cms.EDFilter("VertexSelector",
                                        src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                                        filter = cms.bool(True)   # otherwise it won't filter the events, just produce an empty vertex collection.                                        
                                       )

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#$$
#  '/store/mc/Run3Summer22MiniAODv4/Zto2SmuTo2Mu2ChiTo2Top2Stop_Msmu-450_Mchi-430_ct-030_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2530000/ee97de7d-7c7b-4f01-b48f-2d928ee08fbf.root'
#  'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2022A/RPV_2022A_smu400_neu300_ctau100/MINIAODSIM_1.root',
#  '/store/mc/Run3Summer22MiniAODv4/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5_ext1-v2/50000/1561f73e-71d7-46c8-872f-40e92b7227d1.root'
 '/store/mc/RunIII2024Summer24MiniAOD/Zto2SmuTo2Mu2ChiTo2T2Stop_Par-ct-100-MChi-250-MSmu-400_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/140X_mcRun3_2024_realistic_v26-v3/2540000/b064b486-5ffa-4c59-986e-5e7d6eb831db.root'
#$$
)
)

##########################################################
#
# Setup AK4 jets to be used for analysis
#
##########################################################
#
# Setup JEC factors, see https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/jetsAK4_Puppi_cff.py
#
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import *
process.jetCorrFactors = patJetCorrFactors.clone(src='slimmedJets',
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
    jetSource='slimmedJets',
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
#$$    src = cms.InputTag("updatedJets")
    src = cms.InputTag("slimmedJets") # for 2024
)
process.tightLepVetoJetId = cms.EDProducer("PatJetIDValueMapProducer",
  filterParams=cms.PSet(
    version = cms.string(TIGHTJETIDERA),
    quality = cms.string('TIGHTLEPVETO'),
  ),
#$$    src = cms.InputTag("updatedJets")
    src = cms.InputTag("slimmedJets") # for 2024
)
#
# Module to calculate Pileup Jet ID
# _chsalgos_106X_UL16
# _chsalgos_106X_UL17
# _chsalgos_106X_UL18
# from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL18
# process.load("RecoJets.JetProducers.PileupJetID_cfi")
# process.pileupJetIdUpdated = process.pileupJetId.clone(
#     jets=cms.InputTag("updatedJets"),# JEC corrected jets
#     inputIsCorrected=True,
#     applyJec=False,
#     vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
#     algos = cms.VPSet(_chsalgos_106X_UL18),
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
           DATASET = cms.untracked.vstring(process.source.fileNames),
           isMC =cms.bool(IsMC), 
           YEAR = cms.int32(Year),
           ERA2016 = cms.bool(isPost),
           RochString = cms.string(ROCCORPATH),
           weightFileMVA = cms.untracked.string(TRKBDT),
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
           genEventInfoInput	    = cms.InputTag("generator"),
           LHEEventProductInput     = cms.InputTag("externalLHEProducer"),#source or externalLHEProducer
           genpruned   = cms.InputTag('prunedGenParticles'),
           genpacked   = cms.InputTag('packedGenParticles'),
           vertices    = cms.InputTag('offlineSlimmedPrimaryVertices'),
           mets        = cms.InputTag("slimmedMETs"),
           jets        = cms.InputTag("slimmedJets"), # updatedJetsWithUserData"
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
           jetjecuncdata = cms.string(JECUNCDATA),
           jetjecuncmc = cms.string(JECUNCMC),
           jetjerdata = cms.string(JERDATA),
           jetjersfdata = cms.string(JERSFDATA),
           jetjermc = cms.string(JERMC),
           jetjersfmc = cms.string(JERSFMC)
       )

#-------------------------------------
#process.p = cms.Path(process.prefiringweight* process.egammaPostRecoSeq* process.updatedPatJetsTransientCorrectedNewDFTraining* process.FlyingTop,process.tsk)
process.p = cms.Path(
    process.GoodVertexFilter*
    # process.jetCorrFactors*
    # process.updatedJets*
    # process.tightJetId*
    # process.tightLepVetoJetId *
    # !!  process.pileupJetIdUpdated* # Using Puppi jets so not needed
    # process.updatedJetsWithUserData*
    # !! process.prefiringweight* # I have to figure out why it is not working for run 3
    # !! process.egammaPostRecoSeq* #I have to figure out why it is not working for run 3
    process.FlyingTop
)
# //jet energy corrections
from PhysicsTools.PatAlgos.tools.helpers  import getPatAlgosToolsTask
process.patAlgosToolsTask = getPatAlgosToolsTask(process)
process.pathRunPatAlgos = cms.Path(process.patAlgosToolsTask)

########## output of ntuple
#$$
if isPost:
    process.TFileService = cms.Service("TFileService", fileName = cms.string("RPV_YearB_inputSample.root") )
else:
    process.TFileService = cms.Service("TFileService", fileName = cms.string("RPV_Year_inputSample.root") )
#$$

#$$
process.options.numberOfThreads=cms.untracked.uint32(2)
#$$

