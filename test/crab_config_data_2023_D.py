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
#$$
config.JobType.psetName = 'flyingtop_data_2023_D.py'
#$$
config.JobType.allowUndistributedCMSSW = True

config.JobType.maxMemoryMB = 5000
config.JobType.numCores = 4
# config.JobType.maxJobRuntimeMin = 720
# config.JobType.inputFiles = ['BDT_TRK_240510_ctau100vsEMUdata_NOchi2NOdxyNOdz.xml'
#$$
config.JobType.inputFiles = [
							'BDT_TRK_241119_2023B_ctau100vsEMUdata.xml'
                            ,'BDT_EVT_ALLSIGvsALLBKG.xml'
                            ,'BDT_EVT_ALLSIGvsDYM50.xml'
                            ,'BDT_EVT_ALLSIGvsTTTo2L2Nu.xml'
							,'BDT_VTX_ALLSTEPS.xml'
							,'BDT_VTX_STEP12.xml'
							,'PU_Run2023_BPix_MC.root'
							,'pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-69200ub-100bins.root'
							,'pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-72400ub-100bins.root'
							,'pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-66000ub-100bins.root'
							,'RoccoR2018UL.txt'
							,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt' #./JECDatabase/textFiles/Summer23BPixPrompt23_RunD_V1_DATA/Summer23BPixPrompt23_RunD_V1_DATA_Uncertainty_AK4PFchs.txt
							,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt'#/JECDatabase/textFiles/Summer23BPixPrompt23_V1_MC/Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFchs.txt
							,'Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_MC/Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFchs.txt
							,'Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFPuppi.txt'#./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_MC/Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFchs.txt
							,'Summer23BPixPrompt23_RunD_JRV1_DATA_PtResolution_AK4PFchs.txt' #./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_DATA/
							,'Summer23BPixPrompt23_RunD_JRV1_DATA_SF_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_DATA/Summer23BPixPrompt23_RunD_JRV1_DATA_SF_AK4PFchs.txt
				]

## Specific data options
config.section_("Data")

##############------------------------------Monte-Carlos----------------###############
# config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM'

##############------------------------------Data------------------------###############
# config.Data.inputDataset = '/MuonEG/Run2023D-22Sep2023_v1-v1/MINIAOD'
config.Data.inputDataset = '/MuonEG/Run2023D-22Sep2023_v2-v1/MINIAOD'


# config.Data.inputDataset =           '/Muon0/Run2023D-22Sep2023_v1-v1/MINIAOD'
# config.Data.inputDataset =           '/Muon0/Run2023D-22Sep2023_v2-v1/MINIAOD'
            


# config.Data.inputDataset =          '/Muon1/Run2023D-22Sep2023_v1-v1/MINIAOD'
# config.Data.inputDataset =          '/Muon1/Run2023D-22Sep2023_v2-v1/MINIAOD'

config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 200


#$$ config.Data.totalUnits = 100 
config.Data.totalUnits = -1 

config.Data.lumiMask = 'Cert_Collisions2023_366442_370790_Golden.json'

#$$
config.Data.outLFNDirBase = '/store/user/pvaucell/DATA_EMU_2023D_v2_04_02_2025'
#$$

config.Data.publication = False
#$$
# config.Data.outputDatasetTag = '2018_step3_221228'
#$$

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

# config.Site.whitelist = ['T2_FR_IPHC']
config.Site.blacklist = ['T1_US_FNAL']
