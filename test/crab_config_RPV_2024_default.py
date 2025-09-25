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
config.JobType.psetName = 'flyingtop_temp.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 2
config.JobType.maxJobRuntimeMin = 60
#$$
config.JobType.inputFiles = ['BDT_TRK_250527_2024_ctau100vsEMUdata.xml'
                            ,'BDT_EVT_ALLSIGvsALLBKG.xml'
                            ,'BDT_EVT_ALLSIGvsDYM50.xml'
                            ,'BDT_EVT_ALLSIGvsTTTo2L2Nu.xml'
 			    ,'BDT_VTX_ALLSTEPS.xml'
			    ,'BDT_VTX_STEP12.xml'
			    ,'Run2024_MC.root'
			    ,'Run2024_DATA.root'
			    ,'pileupHistogram-Cert_Collisions2024_378891_386951_GoldenJson-13p6TeV-66000ub-100bins.root'
			    ,'pileupHistogram-Cert_Collisions2024_378891_386951_GoldenJson-13p6TeV-72400ub-100bins.root'
    		            ,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt'
    		            ,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt'
    		            ,'Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt'
   		            ,'Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFPuppi.txt'
   		            ,'Summer23BPixPrompt23_RunD_JRV1_DATA_PtResolution_AK4PFchs.txt'
   		            ,'Summer23BPixPrompt23_RunD_JRV1_DATA_SF_AK4PFPuppi.txt'
			    ,'RoccoR2018UL.txt']
#$$

## Specific data options
config.section_("Data")
#$$

##############------------------------------Monte-Carlos------------------------###############
config.Data.inputDataset = '/Zto2SmuTo2Mu2ChiTo2T2Stop_inputSample_TuneCP5_13p6TeV_madgraphMLM-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v3/MINIAODSIM'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.totalUnits = -1

#$$
config.Data.outLFNDirBase = '/store/user/pvaucell/RPV_2024/250527/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26'
#$$

config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

# config.Site.whitelist = ['T2_FR_IPHC']
config.Site.blacklist = ['T1_RU_JINR']
