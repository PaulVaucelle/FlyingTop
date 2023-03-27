
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

#30cm-------------------------
        # 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_1.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_2.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_3.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_4.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_5.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_6.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_7.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_8.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_9.root',
    #    'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/CMSSW_10_6_20/UDD_bgctau30_smu300_snu250/MINIAODSIM_v16_L1v1_10.root'

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
# $$

##TT samples
####

'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/564C53A2-0646-A24D-B580-AEBF43B22A7B.root',
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
# 'file:/opt/sbg/cms/ui2_data1/blochd/MINIAODSIM/MC/RunIISummer20UL18MiniAODv2/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/50000/FB2A28D0-3315-EE42-A78D-DC38DA4B4B2E.root*
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

process.options = cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound')

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
    # weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm_HighPurity.weights.xml"),  # BDTrecohp
    # weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm_TT_WTrigger.weights.xml"),  # BDTMiniAOD
    # weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm_wVeto.weights.xml"),  # BDTMiniAOD //0.1624
    weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm_V0Veto.weights.xml"),  # BDTMiniAOD //-0.0090
    # weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm_V0_YcVeto.weights.xml"),  # BDTMiniAOD //0.0372
    # weightFileMVA = cms.untracked.string( "TMVAClassification_BDTG50cm_NoVeto.weights.xml"),  # BDTMiniAOD //0.1270

#$$
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
    beamSpot    = cms.untracked.InputTag('offlineBeamSpot'),
    KVFParameters = cms.PSet( 
       maxDelta_ = cms.double(0.01),
       maxReducedChiSq_ = cms.double(225.),
       minChiSqImprovement_  = cms.double(50.),
       maxNbrOfIterations_ = cms.int32(40)
    #    kcvFitter_ =  KinematicConstrainedVertexFitter(),
      )

    # pileup = cms.InputTag("slimmedAddPileupInfo")
    # SV = cms.InputTag("slimmedSecondaryVertices")
)

#-------------------------------------
process.p = cms.Path(process.FlyingTop)
# getattr(process,'slimmedJets').addTagInfos = cms.bool(True)
########## output of ntuple
#$$
# process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple.root") )
process.TFileService = cms.Service("TFileService", fileName = cms.string("Ntuple_test.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("UDD_bgctau50_smu275_snu225.root") )
# process.TFileService = cms.Service("TFileService", fileName = cms.string("TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8_10k_NOdrSigCut.root") )
#$$
# Ntuple_50cm_TTNI_Treader_input.root
# Ntuple_50cm_Neu_Treader_input.root
#$$ 
process.options.numberOfThreads=cms.untracked.uint32(4)

# #Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
# from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
# process = customiseLogErrorHarvesterUsingOutputCommands(process)

# # Add early deletion of temporary data products to reduce peak memory need
# from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
# process = customiseEarlyDelete(process)
# # End adding early deletion
#$$
