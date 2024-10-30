import FWCore.ParameterSet.Config as cms

#$$
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
process = cms.Process("FlyingTop",Run2_2018)
IsMC = True
year = 2018
#$$

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
##----------------------paul--------------------------##
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
## JeC JER for systematics ###########################
process.load("JetMETCorrections.Configuration.JetCorrectors_cff")
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
from CondCore.CondDB.CondDB_cfi import *
######################################################""
#$$CondCore.CondDB.CondDB_cfi
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag

# Global Tags:

#     Data: 106X_dataRun2_v37
#     MC 2016APV: 106X_mcRun2_asymptotic_preVFP_v11
#     MC 2016: 106X_mcRun2_asymptotic_v17
#     MC 2017: 106X_mc2017_realistic_v10
#     MC 2018: 106X_upgrade2018_realistic_v16_L1v1 
# /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/MC/2016
# mcpufile = cms.string("Pileup_MC2018UL_bin100.root"),
#  datapufile = cms.string("MyDataPileupHistogram_bin100.root"),


isPostAPV = False
ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2018UL.txt"
EGERA = '2018-UL'
GT = '106X_upgrade2018_realistic_v16_L1v1'
TIGHTJETIDERA = 'RUN2ULCHS'  #NOTE: Use "RUN2UL16CHS" for UL2016 eras
L1PREFERA = '20172018'
DATAPUFILE = '/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PileupHistogram-goldenJSON-13tev-2018-69200ub-100bins.root'
DATAPUFILEUP = '/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PileupHistogram-goldenJSON-13tev-2018-72400ub-100bins.root'
DATAPUFILEDOWN = '/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PileupHistogram-goldenJSON-13tev-2018-66000ub-100bins.root'
MCPUFILE   = '/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Pileup_MC2018UL_bin100.root'



if year == 2018 :
    ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2018UL.txt"
    GT = '106X_upgrade2018_realistic_v16_L1v1'
    EGERA = '2018-UL'
    TIGHTJETIDERA = 'RUN2ULCHS'
    L1PREFERA = '20172018'

if year == 2017 or year == 2016:
    ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2017UL.txt"
    GT = '106X_mc2017_realistic_v10'
    EGERA = '2017-UL'
    TIGHTJETIDERA = 'RUN2ULCHS' 
    L1PREFERA = '20172018'
    DATAPUFILE = 'MyDataPileupHistogram2017.root'
    MCPUFILE   = 'Pileup_MC2017UL_bin100.root'

# if year == 2016 :
#     TIGHTJETIDERA = 'RUN2UL16CHS'  #NOTE: Use "RUN2UL16CHS" for UL2016 eras
#     DATAPUFILE = 'MyDataPileupHistogram2016.root'
#     MCPUFILE   = 'Pileup_MC2016UL_bin100.root'
#     if isPostAPV :
#         ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2016aUL.txt"
#         GT = '106X_mcRun2_asymptotic_preVFP_v11'
#         EGERA = '2016-UL'
#         L1PREFERA = '2016'

        
#     else :
#         ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2016bUL.txt"
#         GT = '106X_mcRun2_asymptotic_v17'
#         EGERA = '2016-UL'
#         L1PREFERA = '2016'


if IsMC:                                                                                                                                                                                     
    process.GlobalTag = GlobalTag(process.GlobalTag, GT, '')##106X_upgrade2018_realistic_v16_L1v1 default one signal sample                               
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018Rereco_v1', '') ## for 2018 MuonEG    
    
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1','')                                                                                              

#$$                                                                                                                                                                                       
## 102X_dataRun2_Sep2018Rereco_v1 => /MuonEG/Run2018B-17Sep2018-v1/MINIAOD                                                                                                                
# FlyingTopAnalyzer                                                                                          
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
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_1.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_2.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_3.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_4.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_5.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_6.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_7.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_9.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_8.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_9.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/inputFile/MINIAODSIM_v16_L1v1_10.root'
#$$
    )
)


    # 2016ULpreVFP : era='2016preVFP-UL'
    # 2016ULpostVFP : era='2016postVFP-UL'
    # 2017 UL, era='2017-UL'
    # 2018 UL, era='2018-UL' 
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era=EGERA
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
        version = cms.string(TIGHTJETIDERA), #NOTE: Use "RUN2UL16CHS" for UL2016 eras
        quality = cms.string('TIGHT'),
    ),
    src = cms.InputTag("updatedJets")
)

process.tightLepVetoJetId = cms.EDProducer("PatJetIDValueMapProducer",
  filterParams=cms.PSet(
    version = cms.string(TIGHTJETIDERA), #NOTE: Use "RUN2UL16CHS" for UL2016 eras^M
    quality = cms.string('TIGHTLEPVETO'),
  ),
  src = cms.InputTag("updatedJets")
)
#
# Module to calculate Pileup Jet ID
# _chsalgos_106X_UL16
# _chsalgos_106X_UL17
# _chsalgos_106X_UL18
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL18
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
    jets=cms.InputTag("updatedJets"),# JEC corrected jets
    inputIsCorrected=True,
    applyJec=False,
    vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
    algos = cms.VPSet(_chsalgos_106X_UL18),
)
#
# Embed the Jet ID and Pileup Jet ID variables in the jets.
#
process.updatedJetsWithUserData = cms.EDProducer("PATJetUserDataEmbedder",
    src = cms.InputTag("updatedJets"),
    userFloats = cms.PSet(
        puIdDisc = cms.InputTag('pileupJetIdUpdated:fullDiscriminant'),
    ),
    userInts = cms.PSet(
        puId = cms.InputTag('pileupJetIdUpdated:fullId'),
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

from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
    #TheJets = cms.InputTag("updatedPatJetsUpdatedJEC"), #this should be the slimmedJets collection with up to date JECs !                                                               
    TheJets= cms.InputTag('slimmedJets'),
    DataEraECAL = cms.string("None"),
    DataEraMuon = cms.string(L1PREFERA),
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUnctyECAL = cms.double(0.2),
    PrefiringRateSystematicUnctyMuon = cms.double(0.2)
)




process.options = cms.untracked.PSet( )
process.FlyingTop = cms.EDAnalyzer("FlyingTopAnalyzer",
                                #    RochString = cms.string("./FlyingTop/FlyingTop/test/"), 
                                #    Roccor = cms.FileInPath("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"),
                                    DATASET = cms.untracked.vstring(process.source.fileNames),
                                    isMC =cms.bool(IsMC), 
                                    YEAR = cms.int32(year),
                                    ERA2016 = cms.bool(isPostAPV),
                                    RochString = cms.string(ROCCORPATH),
                                    weightFileMVA = cms.untracked.string( "BDT_TRK_240510_ctau100vsEMUdata_NOchi2NOdxyNOdz.xml"),#track selection => previous :BDT_TRK_CTAU10cm_vs_TT.xml //BDT_TRK_ALLSIGvsDYTT_30_01_2024.xml/ BDT_TRK_ALLSIGvsALLBKG.xml // TMVAClassification_BDTG_TRKSEL_.weights.xml
                                    weightFileMVA_EVTS = cms.untracked.string("BDT_EVT_ALLSIGvsALLBKG.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
                                    weightFileMVA_EVTSDY = cms.untracked.string("BDT_EVT_ALLSIGvsDYM50.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
                                    weightFileMVA_EVTSTT = cms.untracked.string("BDT_EVT_ALLSIGvsTTTo2L2Nu.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
                                    # weightFileMVA_HEMI1 = cms.untracked.string("BDT_HEMI1_ALLSIGvsALLBKG.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
                                    # weightFileMVA_HEMI1DY = cms.untracked.string("BDT_HEMI1_ALLSIGvsDYM50.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
                                    # weightFileMVA_HEMI1TT = cms.untracked.string("BDT_HEMI1_ALLSIGvsTTTo2L2Nu.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
                                    # weightFileMVA_HEMI2 = cms.untracked.string("BDT_HEMI2_ALLSIGvsALLBKG.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
                                    # weightFileMVA_HEMI2DY = cms.untracked.string("BDT_HEMI2_ALLSIGvsDYM50.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
                                    # weightFileMVA_HEMI2TT = cms.untracked.string("BDT_HEMI2_ALLSIGvsTTTo2L2Nu.xml"),#evts selection => previous :  BDT_TRK_ALLSignal.xml
                                    weightFileMVA_VTX = cms.untracked.string("BDT_VTX_ALLSTEPS.xml"),#vtx selection : TMVAClassification_BDTG_VTXSEL_.weights.xml
                                    weightFileMVA_VTX_step1 = cms.untracked.string("BDT_VTX_STEP12.xml"),#vtx selection :  TMVAClassification_BDTG_VTXSel_TIGHTWP.weights.xml
                                    mcpufile = cms.string("Pileup_MC2018UL_bin100.root"),
                                    mcpupath = cms.string("pileup"),
                                    datapufile = cms.string("MyDataPileupHistogram2018.root"),
                                    datapileupfileup = cms.string(DATAPUFILEUP),
                                    datapileupfiledown = cms.string(DATAPUFILEDOWN),
                                    datapupath = cms.string("pileup"),
                                #    muoneps1file   = cms.string("NUM_TightID_DEN_TrackerMuons_abseta_pt.root"),
                                #    muoneps1path   = cms.string("NUM_TightID_DEN_TrackerMuons_abseta_pt"),
                                #    muoneps2file   = cms.string("NUM_TkIsoTight_DEN_TightID_TrackerMuons_abseta_pt.root"),
                                #    muoneps2path   = cms.string("NUM_TkIsoTight_DEN_TightID_TrackerMuons_abseta_pt"),
                                #    muoneps3file   = cms.string("NUM_IsoMu24_or_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_DEN_TkIsoTight_and_TightID_abseta_pt.root"),
                                #    muoneps3path   = cms.string("NUM_IsoMu24_or_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_DEN_TkIsoTight_and_TightID_abseta_pt"),
                                    genEventInfoInput        = cms.InputTag("generator"),
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
                                    rhoCollection = cms.InputTag("fixedGridRhoFastjetAll")
                                )

#process.tsk = cms.Task()
#for mod in process.producers_().itervalues():
#    process.tsk.add(mod)
#    for mod in process.filters_().itervalues():
#        process.tsk.add(mod)
#-------------------------------------
# process.p = cms.Path(process.prefiringweight* process.egammaPostRecoSeq* process.updatedPatJetsTransientCorrectedNewDFTraining* process.FlyingTop,process.tsk)
# process.p = cms.Path(process.prefiringweight* process.egammaPostRecoSeq* process.FlyingTop)

#$$
process.p = cms.Path(
    process.GoodVertexFilter*
    process.jetCorrFactors*
    process.updatedJets*
    process.tightJetId*
    process.tightLepVetoJetId *
    process.pileupJetIdUpdated*
    process.updatedJetsWithUserData*
    process.prefiringweight*
    process.egammaPostRecoSeq*
    process.FlyingTop
)
# //jet energy corrections
from PhysicsTools.PatAlgos.tools.helpers  import getPatAlgosToolsTask
process.patAlgosToolsTask = getPatAlgosToolsTask(process)
process.pathRunPatAlgos = cms.Path(process.patAlgosToolsTask)
########## output of ntuple
#$$
process.TFileService = cms.Service("TFileService", fileName = cms.string("inputFile.root") )
#$$

process.options.numberOfThreads=cms.untracked.uint32(4)

# #Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
# from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
# process = customiseLogErrorHarvesterUsingOutputCommands(process)

# # Add early deletion of temporary data products to reduce peak memory need
# from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
# process = customiseEarlyDelete(process)
# # End adding early deletion

# Ntuple_TT_highStat
#$$
