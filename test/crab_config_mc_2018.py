from WMCore.Configuration import Configuration
# More details here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial

config = Configuration()

config.section_("General")
config.General.requestName = '' # output logs directory
#config.General.workArea = 'crab_reco_step2_benchmark_bgctau10cm'
#config.General.transferLogs = True

## Specific option of the job type
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'flyingtop.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.maxMemoryMB = 8000
config.JobType.numCores = 4
# config.JobType.maxJobRuntimeMin = 2630
config.JobType.maxJobRuntimeMin = 720
#$$
# config.JobType.inputFiles = ['TMVAClassification_BDTG50cm_HighPurity.weights.xml']
config.JobType.inputFiles = ['BDTG_SIG50vsBKGtt_filter_NOveto_19varLOST.xml']
#$$

## Specific data options
config.section_("Data")
#$$
# config.Data.inputDataset = '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM'
config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'
#$$
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
config.Data.totalUnits =  20 
# config.Data.totalUnits = -1

#$$
config.Data.outLFNDirBase = '/store/user/blochd/CMSSW_10_6_20_FLY/MC/2018'
#$$

config.Data.publication = False
#$$
# config.Data.outputDatasetTag = '2018_step3_221228'
#$$

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

# config.Site.whitelist = ['T2_FR_IPHC']
# config.Site.blacklist = ['T3_UK_SGrid_Oxford']
