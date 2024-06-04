import FWCore.ParameterSet.Config as cms

IsMC=False
year = 2016
isPostAPV = True

if year == 2018 :
    from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
    process = cms.Process("FlyingTop",Run2_2018)
if  year == 2017:
    from Configuration.Eras.Era_Run2_2017_cff import Run2_2017
    process = cms.Process("FlyingTop",Run2_2017)
if  year == 2016:
    if isPostAPV :
        from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
        process = cms.Process("FlyingTop",Run2_2016)
    else:
        from Configuration.Eras.Era_Run2_2016_HIPM_cff import Run2_2016_HIPM
        process = cms.Process("FlyingTop",Run2_2016_HIPM)
        

#Run2_2016     this is pre HIPM, pre APV
#Run2_2016_HIPM    this is post APV if no HIPM
#Run2_2017
#Run3

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
##----------------------paul--------------------------##
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

#$$
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

# https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_dat 
###--------------------------

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



#ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2018UL.txt"
#EGERA = '2018-UL'
#GT = '106X_upgrade2018_realistic_v16_L1v1'
#TIGHTJETIDERA = 'RUN2ULCHS'  #NOTE: Use "RUN2UL16CHS" for UL2016 eras
#L1PREFERA = '20172018'
#DATAPUFILE = 'MyDataPileupHistogram2018.root'
#MCPUFILE   = 'Pileup_MC2018UL_bin100.root'

if year == 2018 :
    ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2018UL.txt"
    GT = '106X_dataRun2_v37'
    EGERA = '2018-UL'
    TIGHTJETIDERA = 'RUN2ULCHS'
    L1PREFERAECAL = 'None'
    L1PREFERAMUON= '20172018'
    #chsalgos_106X_UL = _chsalgos_106X_UL18
    MCPUFILE   = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/MC/2016/Pileup_MC2018UL_bin100.root"
    DATAPUFILE = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/Data/2018/MyDataPileupHistogram2018_bin100.root"
if year == 2017 :
    ROCCORPATH = "FlyingTop/FlyingTop/data/RoccoR2017UL.txt"
    GT = '106X_dataRun2_v37'
    EGERA = '2017-UL'
    TIGHTJETIDERA = 'RUN2ULCHS' 
    L1PREFERAECAL = 'UL2017BtoF'
    L1PREFERAMUON= '20172018'
    MCPUFILE   = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/MC/2017/Pileup_MC2017UL_bin100.root"
    DATAPUFILE = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/Data/2017/MyDataPileupHistogram2017.root"
if year == 2016 :
    TIGHTJETIDERA = 'RUN2UL16CHS'  #NOTE: Use "RUN2UL16CHS" for UL2016 eras
    MCPUFILE   = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/MC/2016/Pileup_MC2016UL_bin100.root"
    DATAPUFILE = "/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/Data/2016/MyDataPileupHistogram2016.root"
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
	#chsalgos_106X_UL = _chsalgos_106X_UL16

if IsMC:                                                                                                                                                                                     
    process.GlobalTag = GlobalTag(process.GlobalTag, GT, '')##106X_upgrade2018_realistic_v16_L1v1 default one signal sample                               
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, GT, '') ## for 2018 MuonEG    
    

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                            #'/store/data/Run2016B/MuonEG/MINIAOD/ver1_HIPM_UL2016_MiniAODv2-v2/120000/229FBC5A-7299-8642-881A-D7DE49FE5AF3.root' #pre VFP
                                '/store/data/Run2016H/MuonEG/MINIAOD/UL2016_MiniAODv2-v2/140000/11B2462D-A84A-7944-9880-5DBEE38EF19B.root'       # post VFP
                                #'/store/data/Run2017B/MuonEG/MINIAOD/UL2017_MiniAODv2-v1/270000/17AB9CA6-CD61-1149-B492-A487F2B570B9.root'        #2017                          
                                #'/store/data/Run2018A/MuonEG/MINIAOD/UL2018_MiniAODv2-v1/100000/05BD6361-1727-7248-B956-62956572AE23.root'      # 2018
                            )
                        )

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
if year == 2018 :
    from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL18
    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    process.pileupJetIdUpdated = process.pileupJetId.clone(
        jets=cms.InputTag("updatedJets"),# JEC corrected jets                                                                   
        inputIsCorrected=True,
        applyJec=False,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos = cms.VPSet(_chsalgos_106X_UL18),
    )

if year == 2017 :
    from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL17
    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    process.pileupJetIdUpdated = process.pileupJetId.clone(
        jets=cms.InputTag("updatedJets"),# JEC corrected jets                                                                    
        inputIsCorrected=True,
        applyJec=False,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos = cms.VPSet(_chsalgos_106X_UL17),
    )

if year == 2016 :
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


'''
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
'''


process.options = cms.untracked.PSet( )
process.FlyingTop = cms.EDAnalyzer("FlyingTopAnalyzer",
                                   #    RochString = cms.string("./FlyingTop/FlyingTop/test/"), 
                                   #    Roccor = cms.FileInPath("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"),
                                   isMC =cms.bool(IsMC), 
                                   YEAR = cms.int32(year),
                                   RochString = cms.string(ROCCORPATH),
                                   weightFileMVA = cms.untracked.string( "BDT_TRK_240510_ctau100vsEMUdata_NOchi2NOdxyNOdz.xml"),#track selection => previous :BDT_TRK_ALLSIGvsDYTT_30_01_2024.xml/ BDT_TRK_ALLSIGvsALLBKG.xml // TMVAClassification_BDTG_TRKSEL_.weights.xml
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
                                   datapufile = cms.string("MyDataPileupHistogram2016.root"),
                                   #mcpufile = cms.string("Pileup_MC2017UL_bin100.root"),
                                   #datapufile = cms.string("MyDataPileupHistogram2017.root"),
                                   mcpupath = cms.string("pileup"),
                                   datapupath = cms.string("pileup"),
                                   #DATAPUFILE = '/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/Data/2017/MyDataPileupHistogram2017.root'
                                   #MCPUFILE   = '/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/PU/MC/2017/Pileup_MC2017UL_bin100.root'
                                   
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
    #process.GoodVertexFilter*
    process.jetCorrFactors*
    process.updatedJets*
    process.tightJetId*
    process.tightLepVetoJetId *
    process.pileupJetIdUpdated*
    process.updatedJetsWithUserData*
    #process.prefiringweight*
    process.egammaPostRecoSeq*
    process.FlyingTop
)

########## output of ntuple

process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple_Mini_data.root") )

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
