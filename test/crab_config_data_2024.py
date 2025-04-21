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
config.JobType.psetName = 'flyingtop_data_2024.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.maxMemoryMB = 5000
config.JobType.numCores = 4
# config.JobType.maxJobRuntimeMin = 2630
# config.JobType.maxJobRuntimeMin = 720
config.JobType.inputFiles = ['BDT_TRK_240510_ctau100vsEMUdata_NOchi2NOdxyNOdz.xml'
                            ,'BDT_EVT_ALLSIGvsALLBKG.xml'
                            ,'BDT_EVT_ALLSIGvsDYM50.xml'
                            ,'BDT_EVT_ALLSIGvsTTTo2L2Nu.xml'
							,'BDT_VTX_ALLSTEPS.xml'
							,'BDT_VTX_STEP12.xml'
							,'PU_Run2023_MC.root'
							,'PU_Run2023_data.root'
							,'RoccoR2018UL.txt'
							
							
							]

#$$
## Specific data options
config.section_("Data")
#$$

##############------------------------------Monte-Carlos----------------###############
# config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'

##############------------------------------Data------------------------###############
# config.Data.inputDataset = '/MuonEG/Run2024B-PromptReco-v1/MINIAOD'
# config.Data.inputDataset = '/MuonEG/Run2024C-PromptReco-v1/MINIAOD'
# config.Data.inputDataset = '/MuonEG/Run2024D-PromptReco-v1/MINIAOD'
# config.Data.inputDataset = '/MuonEG/Run2024E-PromptReco-v1/MINIAOD'
# config.Data.inputDataset = '/MuonEG/Run2024F-PromptReco-v1/MINIAOD'
# config.Data.inputDataset = '/MuonEG/Run2024G-PromptReco-v1/MINIAOD'
# config.Data.inputDataset = '/MuonEG/Run2024H-PromptReco-v1/MINIAOD'
# config.Data.inputDataset = '/MuonEG/Run2024I-PromptReco-v1/MINIAOD'


# config.Data.inputDataset ='/Muon0/Run2024A-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon0/Run2024B-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon0/Run2024C-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon0/Run2024D-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon0/Run2024E-PromptReco-v2/MINIAOD'
# config.Data.inputDataset ='/Muon0/Run2024F-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon0/Run2024G-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon0/Run2024H-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon0/Run2024I-PromptReco-v1/MINIAOD'

# config.Data.inputDataset ='/Muon1/Run2024A-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon1/Run2024B-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon1/Run2024C-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon1/Run2024D-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon1/Run2024E-PromptReco-v2/MINIAOD'
# config.Data.inputDataset ='/Muon1/Run2024F-PromptReco-v1/MINIAOD'
# config.Data.inputDataset ='/Muon1/Run2024G-PromptReco-v1/MINIAOD'
config.Data.inputDataset ='/Muon1/Run2024H-PromptReco-v1/MINIAOD'
config.Data.inputDataset ='/Muon1/Run2024I-PromptReco-v1/MINIAOD'

config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 200

#$$ config.Data.totalUnits = 100 
config.Data.totalUnits = -1 

config.Data.lumiMask = 'Cert_Collisions2024_378981_384380_Golden.json'

#$$
config.Data.outLFNDirBase = '/store/user/pvaucell/DATA_MUMU_2024_25_09_2024'
#$$

config.Data.publication = False
#$$
# config.Data.outputDatasetTag = '2018_step3_221228'
#$$

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

# config.Site.whitelist = ['T2_FR_IPHC','T2_US_Florida', 'T2_IT_Legnaro', 'T1_FR_CCIN2P3', 'T2_IT_Bari', 'T2_DE_DESY']
# config.Site.blacklist = ['T1_US_FNAL']
