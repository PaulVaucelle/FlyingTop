import FWCore.ParameterSet.Config as cms

##----------------------paul--------------------------##
from Configuration.AlCa.GlobalTag import GlobalTag


##----------------end of paul------------------------##
process = cms.Process("FlyingTop")

process.load("FWCore.MessageService.MessageLogger_cfi")
##----------------------paul--------------------------##
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
# process.load('Configuration.Geometry.GeometryIdeal_cff') # it may happen that this is needed instead of the one above
##https://twiki.cern.ch/twiki/bin/view/Sandbox/MyRootMakerFrom72XTo74X#DDVectorGetter_vectors_are_empty
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
##----------------end of paul------------------------##


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#        '/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
# 'file:MINIAODSIM_v16_L1v1.root'
        'file:MINIAODSIM_v16_L1v1_1.root',
        # 'file:MINIAODSIM_v16_L1v1_2.root',
        # 'file:MINIAODSIM_v16_L1v1_3.root',
        # 'file:MINIAODSIM_v16_L1v1_4.root',
        # 'file:MINIAODSIM_v16_L1v1_5.root',
        # 'file:MINIAODSIM_v16_L1v1_6.root',
        # 'file:MINIAODSIM_v16_L1v1_7.root',
        # 'file:MINIAODSIM_v16_L1v1_8.root',
        # 'file:MINIAODSIM_v16_L1v1_9.root',
        # 'file:MINIAODSIM_v16_L1v1_10.root'
        # '/UDD_bgctau50_smu275_snu225/blochd-2018_step2HLT-b403a189a2d057e62e59ed092120c7f4/USER'
    )
)
# FlyingTopAnalyzer
process.FlyingTop = cms.EDAnalyzer("FlyingTopAnalyzer",
    # weightFileMVA = cms.untracked.string( "TMVAbgctau50withnhits.xml" ),
    # weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm.weights.xml"),  
    weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50sansalgo.weights.xml"),  

    genpruned = cms.InputTag('prunedGenParticles'),
    genpacked = cms.InputTag('packedGenParticles'),
    vertices  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    mets      = cms.InputTag("slimmedMETs"),
    jets      = cms.InputTag("slimmedJets"),
    genjets   = cms.InputTag("slimmedGenJets"),
    electrons = cms.InputTag("slimmedElectrons"),
    muons     = cms.InputTag("slimmedMuons"),
    pfCands   = cms.InputTag("packedPFCandidates"),
)

#-------------------------------------
## Global tag
#$$
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')
#$$
########## output of ntuple
#$$
process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple_Daniel.root") )
#$$

process.p = cms.Path(process.FlyingTop)

# #$$ 
# process.options.numberOfThreads=cms.untracked.uint32(8)
# #$$ 
