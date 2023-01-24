import FWCore.ParameterSet.Config as cms

#$$
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
#$$

process = cms.Process("FlyingTop",Run2_2018)

# import of standard configurations
#$$
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
##----------------------paul--------------------------##
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#$$

# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#$$
#       'file:MINIAODSIM_v16_L1v1.root'
#       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_1.root'

# 70cm -------------------------------------------------------------------------------
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_1.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_2.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_3.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_4.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_5.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_6.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_7.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_8.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_9.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_10.root'

# 50cm---------------------------------------------------
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_1.root',
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_2.root',
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_3.root',
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_4.root',
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_5.root',
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_6.root',
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_7.root',
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_8.root',
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_9.root',
       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_10.root'

#30cm-------------------------
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_1.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_2.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_3.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_4.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_5.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_6.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_7.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_8.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_9.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu300_snu250/MINIAODSIM_v16_L1v1_10.root'

#10cm-----------
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_1.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_2.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_3.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_4.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_5.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_6.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_7.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_8.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_9.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_10.root'
#$$
    )
)

process.options = cms.untracked.PSet(
)

#-------------------------------------
## Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
#$$
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')
#$$

# FlyingTopAnalyzer
process.FlyingTop = cms.EDAnalyzer("FlyingTopAnalyzer",
#$$
#    weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50sansalgo.weights.xml"),  # BDToldreco
    weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm_HighPurity.weights.xml"),  # BDTrecohp
#    weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG_FromBC.weights.xml"),  # BDTminipf
#    weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm_sansntrk10_avecHP.weights.xml"),  # BDTrecohpsansntrk10
#    weightFileMVALost = cms.untracked.string( "TMVAClassification_BDTG50cm_sansntrk10_avecHP.weights.xml"),  # BDTrecohpsansntrk10
#$$
    genpruned = cms.InputTag('prunedGenParticles'),
    genpacked = cms.InputTag('packedGenParticles'),
    vertices  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    mets      = cms.InputTag("slimmedMETs"),
    jets      = cms.InputTag("slimmedJets"),
    genjets   = cms.InputTag("slimmedGenJets"),
    electrons = cms.InputTag("slimmedElectrons"),
    muons     = cms.InputTag("slimmedMuons"),
    pfCands   = cms.InputTag("packedPFCandidates"),
    lostpfCands   = cms.InputTag("lostTracks"),
)

#-------------------------------------
process.p = cms.Path(process.FlyingTop)

########## output of ntuple
#$$
process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple_50_TriggerTest.root") )
#$$

#$$ 
process.options.numberOfThreads=cms.untracked.uint32(4)
#$$ 

#$$
# #Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
# from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
# process = customiseLogErrorHarvesterUsingOutputCommands(process)

# # Add early deletion of temporary data products to reduce peak memory need
# from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
# process = customiseEarlyDelete(process)
# # End adding early deletion
#$$
