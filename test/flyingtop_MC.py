import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_2023_cff import Run3_2023
process = cms.Process("FlyingTop",Run3_2023)

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVRun3Analysis#Prompt_Reco


process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
##----------------------paul--------------------------##
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
# process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")


# process.load("Geometry.GEMGeometryBuilder.gemGeometryDB_cfi")
# process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

# process.es_prefer_GEMGeometry = cms.ESPrefer("GEMGeometryESModule", "gemGeometry")

# process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")

# process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")

# process.load("Alignment.MuonAlignment.muonGeometryDBConverter_cfi")



# process.load('Configuration.StandardSequences.Services_cff')
# process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
# process.load("FWCore.MessageService.MessageLogger_cfi")
# process.load("Configuration.StandardSequences.MagneticField_cff")
# process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
# process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
# process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# process.load('Configuration.StandardSequences.Reconstruction_cff')
# process.load('RecoTracker.Configuration.RecoTracker_cff')
# process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')


## JeC JER for systematics ###########################
process.load("JetMETCorrections.Configuration.JetCorrectors_cff")
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
## !! process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff") Not available for Run3
from CondCore.CondDB.CondDB_cfi import *
######################################################""
#$$CondCore.CondDB.CondDB_cfi
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")

# https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_dat 
###--------------------------
IsMC=True
year = 2018
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
EGERA = '2022-Prompt'
GT = '140X_mcRun3_2024_realistic_v14' #140X_mcRun3_2024_realistic_v14  130X_mcRun3_2023_realistic_v14-v2 106X_upgrade2018_realistic_v16_L1v1
TIGHTJETIDERA = 'RUN2ULCHS'  #NOTE: Use "RUN2UL16CHS" for UL2016 eras
L1PREFERA = '20172018'
DATAPUFILE = '/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/MyDataPileupHistogram2018.root'
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
    process.GlobalTag = GlobalTag(process.GlobalTag, '140X_mcRun3_2024_realistic_v14', '')##106X_upgrade2018_realistic_v16_L1v1 default one signal sample                               
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '140X_dataRun3_v3', '') ## for 2018 MuonEG    
    
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1','')                                                                                              

#$$                                                                                                                                                                                       
## 102X_dataRun2_Sep2018Rereco_v1 => /MuonEG/Run2018B-17Sep2018-v1/MINIAOD                                                                                                                
# FlyingTopAnalyzer                                                                                          

# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.GoodVertexFilter = cms.EDFilter("VertexSelector",
                                        src = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                        cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
                                        filter = cms.bool(True)   # otherwise it won't filter the events, just produce an empty vertex collection.                                        
                                        
                                    )

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                #  '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/03B96C60-30E8-A449-BF9B-7D4BF8E11222.root'

                                # '/store/mc/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2430000/03B96C60-30E8-A449-BF9B-7D4BF8E11222.root'
                                #'/store/data/Run2018D/SingleMuon/MINIAOD/22Jan2019-v2/110000/03549C1A-345E-8244-9D47-19752B13FAC8.root'
                                #'/store/data/Run2018B/MuonEG/MINIAOD/17Sep2018-v1/00000/1044E92E-2236-7547-A380-E4E28961E076.root'
                                #'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/564C53A2-0646-A24D-B580-AEBF43B22A7B.root'
                                #$$
                                #       'file:MINIAODSIM_v16_L1v1.root'
                                # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_1.root'
# 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/TTTo2L2Nu_2024.root'
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/TTTo2L2Nu_2023.root'
# 'file:/store/mc/Run3Summer23MiniAODv4/TTto2L2Nu_TuneCP5CR1_13p6TeV_powheg-pythia8/MINIAODSIM/130X_mcRun3_2023_realistic_v14-v2/2560000/08d7e38d-8b18-4f25-822b-92c92d8f6a73.root'

#        'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_1.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_2.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_3.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_4.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_5.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_6.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_7.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_9.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_8.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_9.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2018/RPV_2018_smu300_neu200_ctau010/MINIAODSIM_v16_L1v1_10.root'


# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_1.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_2.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_3.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_4.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_5.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_6.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_7.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_8.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_9.root',
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC_2017/RPV_2017_smu200_neu180_ctau010/MINIAODSIM_10.root',

)
)


# !! from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
# from EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
# setupEgammaPostRecoSeq(process,
#                        runEnergyCorrections=True,
#                        runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
#                        era=EGERA
#                         ,eleIDModules=[ 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_iso_V1_cff',
#                         'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_RunIIIWinter22_noIso_V1_cff',
#                         'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Winter22_122X_V1_cff']
# )    


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
                                    # # # datapileupfileup = cms.string(DATAPUFILEUP),
                                    # # # datapileupfiledown = cms.string(DATAPUFILEDOWN),
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

process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple_Test.root") )

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
