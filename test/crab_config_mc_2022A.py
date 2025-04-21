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
config.JobType.psetName = 'flyingtop_default_MC2022A.py'
config.JobType.allowUndistributedCMSSW = True

config.JobType.maxMemoryMB = 5000
config.JobType.numCores = 4
# config.JobType.maxJobRuntimeMin = 2630
config.JobType.maxJobRuntimeMin = 720
#$$
# config.JobType.inputFiles = ['TMVAClassification_BDTG50cm_HighPurity.weights.xml']
config.JobType.inputFiles = ['BDT_TRK_241028_2022A_ctau100vsEMUdata.xml'

                            ,'BDT_EVT_ALLSIGvsALLBKG.xml'
                            ,'BDT_EVT_ALLSIGvsDYM50.xml'
                            ,'BDT_EVT_ALLSIGvsTTTo2L2Nu.xml'
 			    ,'BDT_VTX_ALLSTEPS.xml'
			    ,'BDT_VTX_STEP12.xml'
			    ,'PU_Run2022_MC.root'
				,'PU_Run2022EE_MC.root'
				,'pileupHistogram-Cert_Collisions2022_355100_362760_GoldenJson-13p6TeV-69200ub-100bins.root'
				,'pileupHistogram-Cert_Collisions2022_355100_362760_GoldenJson-13p6TeV-72400ub-100bins.root'
				,'pileupHistogram-Cert_Collisions2022_355100_362760_GoldenJson-13p6TeV-66000ub-100bins.root'
			    ,'RoccoR2018UL.txt',
                'Summer22_22Sep2023_RunCD_V2_DATA_Uncertainty_AK4PFchs.txt'
                'Summer22_22Sep2023_V2_MC_Uncertainty_AK4PFchs.txt'

                'Summer22_22Sep2023_JRV1_MC_PtResolution_AK4PFchs.txt'
                'Summer22_22Sep2023_JRV1_MC_SF_AK4PFchs.txt'
                'Summer22_22Sep2023_JRV1_DATA_PtResolution_AK4PFchs.txt'
                'Summer22_22Sep2023_JRV1_DATA_SF_AK4PFchs.txt'
            ]
#$$

## Specific data options
config.section_("Data")
#$$

##############------------------------------Data------------------------###############
config.Data.inputDataset = 
                        '/TTTo2L2Nu_TuneCP5_13p6TeV-powheg-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                         '/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13p6TeV-powheg-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                         '/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13p6TeV-powheg-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                         '/DYJetsToLL_M-10to50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                         '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                        #  '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                         '/WWTo2L2Nu_TuneCP5_13p6TeV-powheg-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                         '/WZTo2Q2L_mllmin4p0_TuneCP5_13p6TeV-amcatnloFXFX-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                         '/ZZTo2Q2L_mllmin4p0_TuneCP5_13p6TeV-amcatnloFXFX-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                         '/ttWJetsToLNu_5f_EWK_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
			             '/TTZToLL_5f_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
			            #  '/TTToHadronic_TuneCP5CR1_13p6TeV-powheg-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
			             '/TTWW_TuneCP5_13p6TeV-madgraph-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM',
                        '/TTToSemiLeptonic_TuneCP5_13p6TeV-powheg-pythia8/Run3Summer22MiniAODv4-130X_mcRun3_2022_realistic_v5-v2/MINIAODSIM'



#$$                          
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 200
# config.Data.totalUnits =  20 
config.Data.totalUnits = -1

#$$
config.Data.outLFNDirBase = '/store/user/pvaucell/'
#$$

config.Data.publication = False
#$$
# config.Data.outputDatasetTag = '2018_step3_221228'
#$$

config.section_("Site")
config.Site.storageSite = 'T2_FR_IPHC'

# config.Site.whitelist = ['T2_FR_IPHC']
# config.Site.blacklist = ['T3_UK_SGrid_Oxford']
