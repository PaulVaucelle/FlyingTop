import FWCore.ParameterSet.Config as cms

#$$
# from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
# process = cms.Process("FlyingTop",Run2_2018)

from Configuration.Eras.Era_Run2_2017_cff import Run2_2016
process = cms.Process('FlyingTop',Run2_2016)
IsMC = True
year = 2016
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
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag

# Global Tags:

#     Data: 106X_dataRun2_v37
#     MC 2016APV: 106X_mcRun2_asymptotic_preVFP_v11
#     MC 2016: 106X_mcRun2_asymptotic_v17
#     MC 2017: 106X_mc2017_realistic_v10
#     MC 2018: 106X_upgrade2018_realistic_v16_L1v1 
isPostAPV = False
if year == 2016 :
    TIGHTJETIDERA = 'RUN2UL16CHS'  #NOTE: Use "RUN2UL16CHS" for UL2016 eras
    MCPUFILE   = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/MC/2016/Pileup_MC2016UL_bin100.root"
    DATAPUFILE = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/Data/2016/MyDataPileupHistogram2016.root"
    DATAPUFILEUP = '/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PileupHistogram-goldenJSON-13tev-2016-72400ub-100bins.root'
    DATAPUFILEDOWN = '/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PileupHistogram-goldenJSON-13tev-2016-66000ub-100bins.root'

    if isPostAPV :
        ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2016bUL.txt"
        GT = '106X_dataRun2_v37'
        EGERA = '2016postVFP-UL'
        L1PREFERAECAL = 'UL2016postVFP'
        L1PREFERAMUON= '2016postVFP'
        #chsalgos_106X_UL = _chsalgos_106X_UL16
    else :
        ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2016aUL.txt"
        GT = '106X_dataRun2_v37'
        EGERA = '2016preVFP-UL'
        L1PREFERAECAL = 'UL2016preVFP'
        L1PREFERAMUON= '2016preVFP'
if IsMC:                                                                                                                                                                                     
    process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_preVFP_v11', '') ## 106X_upgrade2018_realistic_v16_L1v1 default one signal sample                               
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018Rereco_v1', '') ## for 2018 MuonEG    
    
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
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_1.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_2.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_3.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_4.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_5.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_6.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_7.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_9.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_8.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_9.root',
'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2016preVFP/inputFile/MINIAODSIM_10.root'
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
##############################There is an issue when running with 2016preVFP-Ul where : 
# FileInPath unable to find file EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2016_UltraLegacy_preVFP_RunFineEtaR9Gain_v3_scales.dat anywhere in the search path

#
# Setup JEC factors for AK4 jets
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
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import *
process.updatedJets = updatedPatJets.clone(
    addBTagInfo=False,
    jetSource='slimmedJets',
    jetCorrFactorsSource=cms.VInputTag(cms.InputTag("jetCorrFactors") ),
)
#
#  Module to calculate JetID for "Tight" working point
process.tightJetId = cms.EDProducer("PatJetIDValueMapProducer",
    filterParams=cms.PSet(
        version = cms.string('RUN2UL16CHS'), #NOTE: Use "RUN2UL16CHS" for UL2016 eras
        quality = cms.string('TIGHT'),
    ),
    src = cms.InputTag("updatedJets")
)

process.tightLepVetoJetId = cms.EDProducer("PatJetIDValueMapProducer",
  filterParams=cms.PSet(
    version = cms.string('RUN2UL16CHS'), #NOTE: Use "RUN2UL16CHS" for UL2016 eras^M
    quality = cms.string('TIGHTLEPVETO'),
  ),
  src = cms.InputTag("updatedJets")
)
#
# Module to calculate Pileup Jet ID
# _chsalgos_106X_UL16
# _chsalgos_106X_UL17
# _chsalgos_106X_UL18
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL16
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
    jets=cms.InputTag("updatedJets"),# JEC corrected jets
    inputIsCorrected=True,
    applyJec=False,
    vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
    algos = cms.VPSet(_chsalgos_106X_UL16),
)
#
# Embed the Jet ID and Pileup Jet ID variables in the jets.
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

# from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
# updateJetCollection(
#     process,
#     jetSource = cms.InputTag('slimmedJets'),
#     jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#     btagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb',
#                           'pfDeepFlavourJetTags:probb','pfDeepFlavourJetTags:probbb','pfDeepFlavourJetTags:problepb'], ## to add discriminators
#     btagPrefix = 'TEST'
# )

from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
    #TheJets = cms.InputTag("updatedPatJetsUpdatedJEC"), #this should be the slimmedJets collection with up to date JECs !                                                               
    TheJets= cms.InputTag('slimmedJets'),
    DataEraECAL = cms.string(L1PREFERAECAL),
    DataEraMuon = cms.string(L1PREFERAMUON),
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
                                   weightFileMVA = cms.untracked.string( "BDT_TRK_240603_2016pre_ctau100vsEMUdata.xml"),#track selection => previous :BDT_TRK_ALLSIGvsDYTT_30_01_2024.xml/ BDT_TRK_ALLSIGvsALLBKG.xml // TMVAClassification_BDTG_TRKSEL_.weights.xml
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
                                   #mcpufile = cms.string("/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/MC/2017/Pileup_MC2017UL_bin100.root"),
                                   #datapufile = cms.string("/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/Data/2017/MyDataPileupHistogram2017.root"),
                                   #mcpufile = cms.string(MCPUFILE),
                                   #datapufile = cms.string(DATAPUFILE),
                                   mcpufile = cms.string("Pileup_MC2016UL_bin100.root"),
                                   mcpupath = cms.string("pileup"),
                                   datapufile = cms.string("MyDataPileupHistogram2016.root"),
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

