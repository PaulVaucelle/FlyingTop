#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
"""
import os
from optparse import OptionParser


from datetime import date


import CRABClient
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
# from httplib import HTTPException


def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = '',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = '',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def main():

    options = getOptions()

    # The submit command needs special treatment.
    if options.crabCmd == 'submit':
        today = date.today()
        d1 = today.strftime("%d_%m_%Y")
        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()

        config.General.requestName = None
        config.General.workArea = 'MC_EMU_2024_23_06_2025_DY'  #d1

        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'flyingtop_mc_2024.py'
        
        config.JobType.inputFiles = [
            
        'BDT_TRK_250428_2024_ctau100vsEMUdata.xml',
        'BDT_EVT_ALLSIGvsALLBKG.xml',
        'BDT_EVT_ALLSIGvsDYM50.xml',
        'BDT_EVT_ALLSIGvsTTTo2L2Nu.xml',
        'BDT_VTX_ALLSTEPS.xml',
        'BDT_VTX_STEP12.xml',
        'Run2024_MC.root',
        'Run2024_DATA.root',
        'RoccoR2018UL.txt',
        'pileupHistogram-Cert_Collisions2024_378891_386951_GoldenJson-13p6TeV-72400ub-100bins.root',
        'pileupHistogram-Cert_Collisions2024_378891_386951_GoldenJson-13p6TeV-66000ub-100bins.root'
        ,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt' #./JECDatabase/textFiles/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_Uncertainty_AK4PFchs.txt
        ,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt' #./JECDatabase/textFiles/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_Uncertainty_AK4PFchs.txt
        ,'Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23Prompt23_RunCv1234_JRV1_MC/
        ,'Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23Prompt23_RunCv1234_JRV1_MC/Summer23Prompt23_RunCv1234_JRV1_MC_SF_AK4PFchs.txt
        ,'Summer23BPixPrompt23_RunD_JRV1_DATA_PtResolution_AK4PFchs.txt' #./JRDatabase/textFiles/Summer23Prompt23_RunCv1234_JRV1_DATA/
        ,'Summer23BPixPrompt23_RunD_JRV1_DATA_SF_AK4PFPuppi.txt' #RunCv1234_JRV1_DATA_SF_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23Prompt23_RunCv1234_JRV1_DATA/Summer23Prompt23_RunCv1234_JRV1_DATA_SF_AK4PFchs.txt
    
        ]

        config.JobType.maxMemoryMB = 5000
        config.JobType.numCores = 4

        config.Data.inputDataset = None
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 50
        # config.Data.totalUnits = 20 # problem with data/ MC if we want tto quick check, we can uncomment this
        config.Data.outLFNDirBase = '/store/user/pvaucell/MC_EMU_2024_23_06_2025_DY'
        config.Data.outputDatasetTag = 'MC'
        config.Data.inputDBS = 'global' #for signal : phys03 // for bkg :  global 
        config.Data.ignoreLocality = False
        config.Site.storageSite = 'T2_FR_IPHC' # Choose your site. 

        config.General.transferOutputs = True

        

        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDatasets = [
            
            # '/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # '/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # #MUMU
            # '/DYto2Mu_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',        
            # '/DYto2Mu_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',       
            # '/DYto2Mu_Bin-MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # '/DYto2Mu_Bin-MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # '/DYto2Mu_Bin-MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # #EMU
            # '/DYto2Tau_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',       
            # '/DYto2Tau_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # '/DYto2Tau_Bin-MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # '/DYto2Tau_Bin-MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # '/DYto2Tau_Bin-MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',

            # '/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # '/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM',
            # '/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v2/MINIAODSIM'


            "/DYto2Tau-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v3/MINIAODSIM"
            # "/DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/RunIII2024Summer24MiniAOD-140X_mcRun3_2024_realistic_v26-v3/MINIAODSIM"

            		    # '/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/MINIAODSIM', # !! ok
						# '/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v4/MINIAODSIM', # !! ok
						#  '/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/MINIAODSIM', # !! ok'
                        #  '/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM',# !! ok 
                        #  '/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM', # !! ok 
                        #  '/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v4/MINIAODSIM',# !! ok
                        #  '/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3/MINIAODSIM',# !! ok 
                        #  '/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3/MINIAODSIM',# !! ok
						#  '/TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2/MINIAODSIM', # !! ok
						#  '/TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2/MINIAODSIM', # !! ok
						# '/TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2/MINIAODSIM', # !! ok
                        # '/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/MINIAODSIM' # !! ok : TTToSemiLeptonic

			            ]

# /DYto2Mu-2Jets_Bin-MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v6/MINIAODSIM
# /DYto2Mu_Bin-MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /DYto2Mu_Bin-MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /DYto2Mu_Bin-MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /DYto2Mu_Bin-MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /DYto2Mu_Bin-MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v3/MINIAODSIM
# /TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM
# /TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM

        for inDS in inputDatasets:
            # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            config.Data.publication = False
            config.General.requestName = inDS.split('/')[1]
            config.Data.inputDataset = inDS
            config.Data.outputDatasetTag = '%s_%s' % (d1, config.General.requestName)
            # Submit.
            try:
                print("Submitting for input dataset %s" % (inDS))
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print("Submission for input dataset %s failed: %s" % (inDS, hte.headers))
            except ClientException as cle:
                print("Submission for input dataset %s failed: %s" % (inDS, cle))

    # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print("-"*len(msg))
            print(msg)
            print("-"*len(msg))
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print("Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers))
            except ClientException as cle:
                print("Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle))


if __name__ == '__main__':
    main()
                                                                                                                      
