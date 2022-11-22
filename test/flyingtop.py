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
        'file:MINIAODSIM_v16_L1v1.root'
    )
)
# FlyingTopAnalyzer
process.FlyingTop = cms.EDAnalyzer("FlyingTopAnalyzer",
    # weightFileMVA = cms.untracked.string( "TMVAbgctau50withnhits.xml" ),
    weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm.weights.xml"),  
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

########## output of ntuple
#$$
process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple.root") )
#$$

process.p = cms.Path(process.FlyingTop)
