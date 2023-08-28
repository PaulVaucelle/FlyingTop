
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

#$$
# from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
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
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
#$$
# setupEgammaPostRecoSeq(process,era='2018-Prompt')
# process.p=cms.Path(process.egammaPostRecoSeq)
#$$

# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000000) )

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#$$
#       'file:MINIAODSIM_v16_L1v1.root'
#       'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_1.root'

# # 70cm -------------------------------------------------------------------------------
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_1.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_2.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_3.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_4.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_5.root',
    #   'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_6.root',
    #   'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_7.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_8.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_9.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau70_smu250_snu200/MINIAODSIM_v16_L1v1_10.root'

# # # 50cm---------------------------------------------------
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_1.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_2.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_3.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_4.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_5.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_6.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_7.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_8.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_9.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau50_smu275_snu225/MINIAODSIM_v16_L1v1_10.root'

# # #30cm-------------------------
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_1.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_2.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_3.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_4.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_5.root',
    #   'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_6.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_7.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_8.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_9.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_10.root'

# # 10cm-----------
    #      'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_1.root',
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_2.root',
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_3.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_4.root',
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_5.root',
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_6.root',
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_7.root',
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_8.root',
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_9.root',
    #     'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau10_smu250_snu200/MINIAODSIM_v16_L1v1_10.root'
# $$
# 'file:/store/mc/RunIISummer20UL18MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/04A0B676-D63A-6D41-B47F-F4CF8CBE7DB8.root'



##TT samples
####
# # #Everything takes abotu 3 hours in interactve mode
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/564C53A2-0646-A24D-B580-AEBF43B22A7B.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/6C51F452-445D-7547-B585-BDA1278294AE.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/0023597E-6F49-214D-A518-CC9330A404BD.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/00B564A1-86ED-2C45-8159-5FE8CA552C71.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/041BD979-6754-BC42-8A91-8698F5744F27.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/0938183D-C346-BE42-AD41-EAAB4DD5F61C.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/11E8D8AA-5485-2B43-937C-2137A3A95624.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/1207B318-4CEB-AF41-B68C-2BD3598339F0.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/13B67860-F6B0-EF4F-AB59-2646C0AAA906.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/1C224BAE-6F90-DB41-A495-1A976B13F267.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/210F0D88-DC80-484E-BC6D-6929507FDEE2.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/25720137-0046-D346-8451-1F3EC899F9FE.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/30E5F4EF-74DF-4A40-9CE8-17E329305027.root',
    # # # #328k events (30min)
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/3440172F-14FF-9640-845A-8A7505EF06F7.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/37134442-447E-DE48-9B51-6EDE08191165.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/3813D538-73AF-CE48-BD3C-DCEFBFE88F9E.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/3AF7AAD5-3376-3C45-B6F8-4DBB8CC5056D.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/43183F74-7A2D-8941-8317-6FEAA34987BE.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/45B2708A-78DC-A340-A7DE-AE6174133BDD.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/4EBF1ED8-105C-B549-99FA-5C757EB79470.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/500F6E0A-B2D7-CC45-80BA-480F8BE85552.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/564C53A2-0646-A24D-B580-AEBF43B22A7B.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/57913C40-2F64-8E45-9BBD-FBF0622A5771.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/582F5B79-F9CF-D543-84A6-F0B5953D270D.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/58F7CA73-DEE8-F846-9469-D62EF845C11C.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/5DD6C708-8AAA-0C4E-BD44-8F9F762869CF.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/6800A177-83B7-4F40-8397-C8F90D9AE37C.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/6A74FF77-02A9-9147-B126-5FC44AFC05B6.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/6C51F452-445D-7547-B585-BDA1278294AE.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/6E8C7F63-2E85-AD4A-8079-441882EF0CB2.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/71482C68-8212-BB43-8A72-CA4DB92B1716.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/71D9D8E8-2579-0F40-8F2B-91585B2BE98B.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/72470E8D-E4D3-6342-8974-B8617CC3AE74.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/7983AC65-6B64-AE4B-B5E8-5D2D7A5D93F4.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/7B6D4654-9443-A944-98EF-779F1C5A6D1B.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/7D6AFFF8-DFB0-DA4D-BB48-5FD6FFD2218A.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/7EEAEA3B-E8B3-F548-9749-63BA57A3A0A6.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/810DFBF5-5DFC-F141-98D9-54EA8E0561C7.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/829F1307-B3CB-4348-8134-B261BAEB603E.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/82AA756A-C136-5E42-A95B-B840C517ECD9.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/83FF1415-978D-5B45-9976-2F6388E3CE44.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/8A52BDD4-78C8-9B4E-AE94-D3FBA76C5600.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/8ACBAF44-7B32-D346-BD40-7861614DF836.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/8CFEDB6B-E310-294B-A023-B3F6E1370084.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/8E7B2251-5E93-4F4C-9B18-1C51D3902672.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/900EBC0E-2C1D-E14A-AA72-2B9C1807ADD8.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/9602191D-73F4-E242-8447-F599912FCE0A.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/9741B695-0711-344C-8260-B0BCD71FAE38.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/A07147D5-0CB7-6F4C-9BC5-607B41BD1AEC.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/A3864815-A3E2-7F4A-9BCE-A059C762DFA3.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/ADA73145-983F-7B4D-B121-A4C58925EBA3.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/AE736E74-4B48-1F48-B9BB-9718AE6953D8.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/B21F05AE-CE6E-0841-BDE2-CC0E5269043F.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/C014E57E-71AF-BF46-8F7E-B2170531D329.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/C4286A37-51D5-AD4B-969F-EE50E7A28244.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/C4476F07-0641-F74B-8D6E-90CA536441CC.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/CE7010E8-D38D-794E-9461-455B79D44D66.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/D387C571-45F1-1844-BFF0-853020F9F91F.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/D5D2F8C4-0ABC-734B-BA5A-1828C8A32B4C.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/D86E503C-53E3-4649-9125-BDC48EC8BBC3.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/D8F91B35-84C3-6D4E-AC2E-26A56C88A862.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/DD79E24E-DF1D-9342-A9E9-018FA0441F96.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/DDCF1232-9E92-444E-8115-998CB9052A8A.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/E37DB1D7-18C2-0F4E-8DDD-7184CAF818D7.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/E43C3B2D-3CE0-8741-AF69-3CA53F701A9D.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/E5DCEF5C-F2D0-074B-A20A-96756C2FD5D6.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/E70064F0-E134-5D4C-9CAC-F2CCBF124854.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/EBD623A1-84F8-734E-B912-835D23555E5E.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/EC59C293-5805-D046-97F1-367668D45BB5.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/F4EABEBA-0405-254F-BBF9-5B8B6B727973.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/F8DA96F9-7958-AB44-87AF-914D435896CA.root',
    # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/FB2A28D0-3315-EE42-A78D-DC38DA4B4B2E.root'
    

    # # //TTTo2L2Nu
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/04A0B676-D63A-6D41-B47F-F4CF8CBE7DB8.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/04ED7864-0830-774F-A499-D26AC997401C.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/093341BD-F5AE-404D-BF2C-157BAF2B9AC5.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/0B7957F4-6494-7B4C-BD85-E682DA1EFA4F.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/0C6623EF-B101-694A-8904-D7578B1093C8.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/0CA220FA-8ED3-824A-8445-0599C8754362.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/15D2D06A-6B01-894D-9A2B-A9CBBBA7FA2B.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/18D0B062-0B22-1C42-8437-73179AAC2A1B.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/21B996BC-8CB3-1A40-ACE3-64DDECADDDE8.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/25840049-8F8B-B449-8A69-57D711B71239.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/27A87F71-7D93-D547-AD16-31084E549053.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/2B6DCD9D-A589-0F48-91D5-06C355FFBDA1.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/2CF51C43-067F-3A4C-A47E-508DB497B5C8.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/315847FE-0611-E144-B82E-EC6E26C66490.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/388FA3C8-8312-5F44-A675-77C8E63522EA.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/3BA34C6E-E559-DE49-AFC2-C52294644E27.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/3BBDC332-6E2D-EB4F-896E-C49B3485DB78.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/3D051F63-4F6D-D147-8367-3054535A1766.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/3EE11089-6746-7049-A5B3-6F66A5162607.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/479256D0-E6B6-454B-81B1-DAA37E532B87.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/4C0D4F9E-176A-704C-B2F4-A6850456738C.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/4DF9E033-6718-5045-A9BB-F65C1C3BEE16.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/51E0D3B6-BF85-0F4F-9F44-B0DAA9805912.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/54F97890-7FF0-DE41-98D5-E41A56892631.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/57DBB937-4C72-6B45-9CD4-1B2C5671EBDE.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/604BEB4A-B3B4-D846-9FFA-1DE1335888C5.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/60BED2E3-37BA-7843-A629-3E4C6BB76A99.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/62128EF3-C14B-2845-A1A5-76C03E2FE10F.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/65F17326-5D17-3B41-9992-95D46B49A583.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/708D8761-95E4-244B-A678-7FAA965B15E1.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/710DC2B8-40C1-7445-B026-2BB2CC9F5515.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/73BE428C-D34E-DD4F-B147-85299F1F914B.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/743EBE98-9E3F-154C-AEEA-2A9951E027C4.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/791F8D83-E992-F645-AC75-087EB7292E20.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/83DEFE89-905A-5448-BE2B-C7DD049A5659.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/8498CA05-C30C-3B4A-A800-C77B3E19503E.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/88A15057-5447-E64F-AC43-DE9723DFB0B7.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/91CD1A2F-85BF-9144-91A2-32FB5C466CBD.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/95965DF4-2D86-DB4C-9818-096E61BCFDFF.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/978450D8-3B00-D64C-BAFB-3CDBCF951DD7.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/9D704B6A-B565-5A43-B5DD-B73A5A014582.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/A4508DAD-6AD8-2C4F-BC7F-D31BAEF54BE6.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/A6DA2625-8B3F-464E-8E64-E9F24A6217EF.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/AC687BC7-76F7-B048-9367-9C29E9FF57F6.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/B42E2585-848A-A64A-A458-C70B31114D98.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/B4EEF49C-E3C5-8144-B774-6D69AF5A9341.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/B70F044C-491F-1C4C-A039-625493131004.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/B7C1E24E-6CCC-3C45-AB21-3563D770B863.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/B86867C3-7DC0-1045-80E5-9F060A4B0547.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/B9535C20-A60D-5D47-864F-6A9D40B34C72.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/BC2272E7-79E8-E942-B321-18B9858C35CC.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/BCB59DB8-3B7B-824B-AC5C-C6EF6A16B44D.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/BDEC421C-D547-7648-9786-B1628BA027EE.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/BE3415D7-2F4D-6042-A53B-58AD5FE731E7.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/BF58DC22-5591-6F45-A8D7-9B061C069EC3.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C18EF78E-6C7D-F547-BCA3-6CAF0BF65355.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C1946A27-F0C6-C64F-9BE8-9128A95ECDCC.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C3730241-E766-0440-88D1-EC5A195DD0F7.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C398B198-8A24-414E-9EC9-59A160096573.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C48A0239-487B-BC4D-AD01-C8F3A6E7D6D9.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C5824C5E-68ED-7247-9129-A63F1C458669.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C76FF9AD-F4FF-FA4F-AFBB-E73A4B417ECD.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C78CD641-3E2E-3B43-B707-BDD2140FC320.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C7E5F05D-2DD8-3B40-8EFF-2B2115122CB4.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/C7F868D2-7A62-4C43-8FC3-E0D270D4035C.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/CF741A25-D0DF-3642-9BA0-B9C3668092D5.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/D71DDA43-8E09-0945-80AA-CE3D5297FD70.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/D731C7FD-7B6E-AB41-A753-A7A599400922.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/D83203E1-714F-FB45-8067-4C3BA07CDF61.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/D97C4A43-5CBE-5342-BE3C-AEA8A6BE10E8.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/DB5D0369-42BC-DE4B-B4FA-04DE84AAADB0.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/DBFB0B2F-86F0-CC4F-B2C4-C09DD6A12D81.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/E108A05B-B5DD-3845-A6A6-9944CD63B126.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/E27C9AEE-3D9C-BB41-9BBF-F53E146A3E25.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/E687E8AB-8966-1443-B034-8657BE7D1FFF.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/E6926812-F42F-AD4B-8989-B11831C73DBA.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/E7DE5B21-4D83-2644-85C8-9C1637674143.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/E97BCD93-C443-2842-97AF-1CA7EB491B84.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/E9BD7694-9EE0-944E-ADEE-115C434B6093.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/EB1BCA1E-1872-2A4D-8163-40EBDB22D532.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/EB8A4EDB-6F5C-7F4F-A837-F311125E0B88.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/EC712038-1B86-114F-BE1E-63CBEB16DBFF.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/F0EE0AD1-7BDB-154D-82CD-4FAEA1C39442.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/F42CCCB1-DD90-0A44-AD86-FEA804C212B3.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/F44E0965-65FA-1149-A833-BB14D305C286.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/F48F3794-2C9E-0A44-BF3E-EAC1369F91DE.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/F4D0BED7-5979-4441-9201-DF8D12C9396C.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/FA822740-577A-8C44-9352-8ADFFFDBFCB3.root',
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/00000/FE8E6D38-B6D2-1540-8DD4-C36B5984A34C.root'

# DYJEtsToLL
    # 'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYLLtoJets_M50/001C8DDF-599C-5E45-BF2C-76F887C9ADE9.root'

   'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/001C8DDF-599C-5E45-BF2C-76F887C9ADE9.root',
    'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/01522709-8919-C542-91B2-2262F3995F48.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/015753DA-CD2E-F546-9A7B-9DD451DEA159.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/01C7FCCE-F23B-2242-BF72-73AA8BAF4C47.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/02991CFE-6B1C-6148-8837-3D68EB56C8C9.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/02FD7B88-3EDC-C64D-8E18-F9A8F9E7E7DF.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/03311968-AA82-2F46-8949-F39CBA6E33CC.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/035BDD75-657D-3E4D-BB95-06F715874873.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/03A89FD2-3EA0-5F47-9877-CA74F8A01950.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/059DAA5D-FF55-C64B-BA61-2F4D084E3A1E.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/05EC91D0-CAED-2A4B-B25A-4189C5C951F6.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/069E946D-802D-8A4B-BA0B-E73E3CDC7BCD.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/06D81BB9-375D-C746-AF50-67BD16501B82.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/06FEF028-669E-3642-81EE-27B06A97BC9C.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/077533C3-85C2-1C48-A029-34E755781A26.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0789E1C9-5C74-A44D-AFB9-F46EC27DEB10.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/08B421E5-B5B0-3940-92D7-E7421CE78AA8.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/093E3572-4578-D342-B2A6-D3FD4E049C0E.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0997B6BB-EB0D-CD4E-8F7E-7F145E7A9AEE.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0A3663AA-EB69-6344-A787-E622693E26A6.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0B706D8E-F68B-DD45-A986-E4A3689A4531.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0BC4DB6D-A6D2-464A-AA77-BC1F7E4E52F2.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0CEE97C2-FC2E-3848-A9AA-6EE57219F97C.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0D720AB7-12F3-144B-B81A-5B5D603C2BC7.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0DB0D698-946B-7149-AEFA-200F94986EBB.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0DC91181-1A95-D24C-B58A-C3FE6FD64F94.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0E9A1DE3-8B33-8948-8ED7-9134D88E6004.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0EA4DA82-4D8B-3B49-9AB3-07128CAE0296.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0F5DDB72-13E4-0945-A34E-959D60A70868.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/0FC94000-7AF7-4A47-A3EE-721BDDC32DF9.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/11103D7B-D884-7744-BF07-536EA418F6F0.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/112B54FA-9CF7-0C42-A228-4BFE04BDAA8D.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/11827714-A478-2D4A-A968-D5C340DC2335.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/11C2016B-D8C8-1641-95B4-42930ADC4335.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1218A43F-2F71-B649-9963-19166FC68509.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/12A5308C-C8FA-DC42-9813-B45BBDF4A56F.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/13B2279B-8ADC-E945-990A-581F3408C35E.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1578F2E4-D79D-7B4A-839F-D037933AB67D.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/16C21BAB-0411-ED48-996E-37BB2195B011.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/16E6ECD2-560E-6344-93FE-4A727A195BA1.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/16F638B5-48A5-5D48-892A-7552C7457620.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/17803CB6-AA01-FB4A-951D-A090C273654C.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/17AAEF93-B0AE-9440-BC3B-C7511157A5B1.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1819A02B-9630-1E40-BD85-232FAF014F3E.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1873AC42-133F-BA44-BC24-B4315F1B38C1.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/18916CB9-2107-9647-835C-CB2B6D45E5C3.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/193B4ACF-4FDC-5746-B6E1-A5FB9BEA6CC5.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/19841D9D-54C4-9640-B42D-AB7A6A84FA0C.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/19A4BD1F-EBAA-CC48-A0F1-0E628BF1B028.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/19A51DDC-742F-EF47-91D9-D1010AB742A0.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/19A91FDE-44E5-F248-A07B-B4CFD30BF149.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/19DCEA8C-7F8F-F34F-A7A5-517D27D65DAF.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1A00818F-24B4-5342-BDB7-F2DEB6686091.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1A4E8B78-2535-5A4C-809D-8E3F5F0010AD.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1AE43FCE-69A2-DB45-8277-3B6A476D0F0B.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1BF69BD0-39B3-4642-AD14-277EF956D3EB.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1C431755-105A-464D-B099-416FF564D5D5.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1CBF78A6-A4FC-3F46-A7E5-FB7A407D4D9F.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1D88B92C-DCCC-C442-8323-0407EBEEEE77.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1E0DDC2C-F8BB-374D-A00B-77627AFB92F0.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/1E8FF42F-58B6-E941-AC74-8D8131B4BFDE.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/20147561-35B0-5E45-AD55-91BFEFE491D1.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/206C1852-FD74-2344-8DFF-E5792CB2EAAB.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/20F2105A-F471-4147-91CE-C0EC1F893686.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/22242799-25ED-8146-87E8-471CF72D4334.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/230D9833-C988-3F4D-8539-08443DB5D74B.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/266573BB-CF15-3B4C-88E4-710E8FAD70E7.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/2792BBE6-FB96-1E4C-BBA7-0B8EB6ED06BA.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/2821C780-9B3A-324F-9B18-F2B9C49E73BE.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/284B9F60-CE95-9142-9ABC-49DE67762CBF.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/28BBD339-BEAA-B148-A0A0-FC85E6C5451D.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/29D45569-4826-2344-9193-F4071A720A9F.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/2AC9D3D8-B47C-184E-8FC2-71CF5D59CB36.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/2BC4FCA4-3F34-2141-8D3E-B5553E7582F4.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/2D6B952F-122C-B440-8398-643F148FD554.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/2DBAAA93-1984-CF43-9F1E-110E9AB37FD0.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/2F1EC487-5FE5-C049-A504-9B091D2F6E9F.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3106808A-A2CE-864B-BA4A-B5260FD012B8.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/31C04865-65FD-BC41-9675-7E21CCAEC6B6.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/35557DDE-C2E0-A54A-A092-1C8984910B31.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/35E7C1BB-84DF-CC4E-B9E0-5931B0C13DE7.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/369F8A2E-11FB-D848-8FD7-F2C737640DEA.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/36C7C2F8-F785-3242-8427-7ABFBCF7684E.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/371B67C7-3399-FD41-801E-64F5C1C21EB3.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/386BF912-8016-E644-8D46-4D72AB950C96.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3A282FF1-2262-E142-A3DE-A3128577160B.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3BFCCAFE-37F7-F841-8066-10E1B298F2AC.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3C362AD9-969F-2146-A8D6-A4E8E26CB4B4.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3C5559A3-C5BF-E246-8635-A569B89E474C.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3C5F5250-8AE5-9E4D-B7A5-9C158212EC82.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3C76F160-EB37-6A46-8253-203CDDAD34CF.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3CEF106A-8B25-5A42-AAE9-860D1ECBDB34.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3DC47A03-48B3-F945-BF60-72C09D37B38A.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3DDA7F11-6FAE-AE4E-AE40-3D05F2A1620D.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3DF4BA2F-444A-2F46-8DF8-E8074D251C14.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3DFD0B21-8C09-E646-9423-F944E550199D.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3E00A6A7-2AD9-6D4E-9075-61E93E04B88D.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/3E32083D-F30A-984B-BAEB-725C945B25B2.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/402D73C7-FE83-C14D-A3BD-C55A37D0E45D.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4039836C-844E-E24E-B130-550A901F0349.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/40667CCA-1E86-D44A-9067-6E449A87530A.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/41C2071D-EE5A-5849-80BB-A2614A3866C1.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/42BA2532-0C5D-A744-987A-7B3A6C04F071.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/431192D1-936F-CD49-977A-97A598759694.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/454F59F4-6AB3-C74F-964B-CD8081FB4586.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/470467C4-F735-784E-AD4F-5601C8C13169.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/47861CE9-2547-9646-A872-47315C8C8871.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/47C7108E-DAE2-3B45-A7E0-8297B0B4DE68.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/48A79C40-A538-E743-A1EE-E383F67E7256.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/48E86767-BDFD-1E47-B7FE-415960CCBF6B.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4A55B8AC-8A4E-6C49-82EE-EB5D2AC8BAF1.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4AA50662-AE98-3243-ACCB-BF8DCA3C82CB.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4AB9B611-E012-A044-BE10-970850D30AC1.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4B536C2B-F075-9C49-B032-C577417ABE79.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4B7084B1-C7E2-4B48-B473-EC9C046A249B.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4DC81A62-66BC-6248-B3C0-20938EFFDB3D.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4E63A500-7C8C-B54B-977F-644C5577E214.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4F132011-38DE-764C-B52A-0983AA61CC3B.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/4FD27CE6-0888-064F-A26E-673BC9B42E79.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/508C54AF-9980-A246-888A-3DD3D72289FF.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/50F68336-FA50-374F-9418-F2183B3A2A48.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/516635C2-DA22-FD43-9737-D7D6DF50B1E0.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/5169D9A1-30F3-C647-B23B-A54B5EF73968.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/517F989E-AC8F-E342-9452-8D1346C447A0.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/51A8D7EA-4A3C-2142-A701-8FC6497515A9.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/51BDE0E6-5214-A548-A228-8122F8C35BF3.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/534C146E-E0FD-8147-85B9-787294A22681.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/53ED64BD-405D-8E43-BC48-3155096AAAE6.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/551AE7FD-5F1B-654C-9C6B-83B0B32BD78B.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/55353174-C59F-ED48-9AF6-81A4805916C2.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/55A76F90-A2F7-D24A-BE73-F56CA50AB7FD.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/56D22CD0-B03E-DB4F-A246-3C330112C2AF.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/579DD16E-2DCC-DB4A-AD84-5541AA9B9653.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/57E6619D-8997-744F-AB62-B5FF98247557.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/581F25D2-5A6E-8E46-9D6A-EC109FB42051.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/5845CE12-830B-B545-9D8D-254DB201CA5D.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/594EF710-6AEE-D346-8013-A5687819947A.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/5A5CEC5D-0D44-C645-91BC-D6C05E29F5AF.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/5C76A0F3-DA3E-B743-B615-63C9065A0428.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/5CCCF4E1-81E7-8445-9173-9AA6FA3A8F52.root',
'file:/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/MC/DYJetsTo_LL_M50/5EE47880-81FA-CD42-A9E7-FC556E8CC9E1.root'

    
    )
)

# //btagging
# from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
# patAlgosToolsTask = getPatAlgosToolsTask(process)
# process.outpath = cms.EndPath(process, patAlgosToolsTask)
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
#$$$$    btagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probc','pfDeepCSVJetTags:probudsg','pfDeepCSVJetTags:probbb','pfDeepFlavourJetTags:probb'], ## to add discriminators
    btagDiscriminators = ['pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb',
                          'pfDeepFlavourJetTags:probb','pfDeepFlavourJetTags:probbb','pfDeepFlavourJetTags:problepb'], ## to add discriminators
    btagPrefix = 'TEST'
)

process.options = cms.untracked.PSet(

)

#-------------------------------------
## Global tag
from Configuration.AlCa.GlobalTag import GlobalTag
#$$
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')##106X_upgrade2018_realistic_v16_L1v1 default one 
#$$
## 102X_dataRun2_Sep2018Rereco_v1 => /MuonEG/Run2018B-17Sep2018-v1/MINIAOD
# FlyingTopAnalyzer
process.FlyingTop = cms.EDAnalyzer("FlyingTopAnalyzer",
#$$
    weightFileMVA = cms.untracked.string( "BDTG_SIG50vsBKGtt_filter_veto_15var.xml"),#track selection

# BDTG_SIG50vsBKGtt_filter_veto_15var.xml => Daniel => version sans Yc Veto bug

    weightFileMVA_EVTS = cms.untracked.string("BDTG_EvtsSel_SIG50vsAll.xml"),#evts selection

    weightFileMVA_VTX = cms.untracked.string("BDTG_VtxSel_SIG50vsAll.xml"),#vtx selection
    weightFileMVA_VTX_step1 = cms.untracked.string("BDTG_VtxSel_SIG50vsAll_Tight.xml"),#vtx selection
    genpruned   = cms.InputTag('prunedGenParticles'),
    genpacked   = cms.InputTag('packedGenParticles'),
    vertices    = cms.InputTag('offlineSlimmedPrimaryVertices'),
    mets        = cms.InputTag("slimmedMETs"),
    jets        = cms.InputTag("slimmedJets"),
    genjets     = cms.InputTag("slimmedGenJets"),
    electrons   = cms.InputTag("slimmedElectrons"),
    muons       = cms.InputTag("slimmedMuons"),
    pfCands     = cms.InputTag("packedPFCandidates"),
    lostpfCands = cms.InputTag("lostTracks"),
    Kshorts     = cms.InputTag("slimmedKshortVertices"),#recoVertexCompositePtrCandidates_slimmedKshortVertices__PAT
    Lambda      = cms.InputTag("slimmedLambdaVertices"),#recoVertexCompositePtrCandidates_slimmedLambdaVertices__PAT
    beamSpot    = cms.untracked.InputTag('offlineBeamSpot')
    #    kcvFitter_ =  KinematicConstrainedVertexFitter(),
    # pileup = cms.InputTag("slimmedAddPileupInfo")
    # SV = cms.InputTag("slimmedSecondaryVertices")
)

#-------------------------------------
process.p = cms.Path(process.FlyingTop)
# getattr(process,'slimmedJets').addTagInfos = cms.bool(True)
########## output of ntuple

process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple_.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("UDD_bgctau50_smu275_snu225.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8_10k_NOdrSigCut.root") )

# process.options.numberOfThreads=cms.untracked.uint32(4)

# #Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
# from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
# process = customiseLogErrorHarvesterUsingOutputCommands(process)

# # Add early deletion of temporary data products to reduce peak memory need
# from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
# process = customiseEarlyDelete(process)
# # End adding early deletion

# Ntuple_TT_highStat
#$$
