#!/usr/bin/env python3
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
        config.General.workArea = 'MC_MUMU_2023_C_23_04_2025'  #d1

        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'flyingtop_mc_2023_C.py'
        config.JobType.inputFiles = [
							'BDT_TRK_241122_2023A_ctau100vsEMUdata.xml'
                            ,'BDT_EVT_ALLSIGvsALLBKG.xml'
                            ,'BDT_EVT_ALLSIGvsDYM50.xml'
                            ,'BDT_EVT_ALLSIGvsTTTo2L2Nu.xml'
							,'BDT_VTX_ALLSTEPS.xml'
							,'BDT_VTX_STEP12.xml'
							,'PU_Run2023_MC.root'
							,'PU_Run2023_BPix_MC.root'
							,'pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-69200ub-100bins.root'
							,'pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-72400ub-100bins.root'
							,'pileupHistogram-Cert_Collisions2023_366442_370790_GoldenJson-13p6TeV-66000ub-100bins.root'
							,'RoccoR2018UL.txt'
							,'Summer23Prompt23_V1_MC_Uncertainty_AK4PFchs.txt' #./JECDatabase/textFiles/Summer23Prompt23_RunCv4_V1_DATA/Summer23Prompt23_RunCv4_V1_DATA_Uncertainty_AK4PFchs.txt
							,'Summer23Prompt23_V1_MC_Uncertainty_AK4PFPuppi.txt' #./JECDatabase/textFiles/Summer23Prompt23_V1_MC/Summer23Prompt23_V1_MC_Uncertainty_AK4PFchs.txt
							,'Summer23Prompt23_RunCv1234_JRV1_MC_PtResolution_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23Prompt23_RunCv1234_JRV1_MC/
							,'Summer23Prompt23_RunCv1234_JRV1_MC_SF_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23Prompt23_RunCv1234_JRV1_MC/Summer23Prompt23_RunCv1234_JRV1_MC_SF_AK4PFchs.txt
							,'Summer23Prompt23_RunCv1234_JRV1_DATA_PtResolution_AK4PFchs.txt' #./JRDatabase/textFiles/Summer23Prompt23_RunCv1234_JRV1_DATA/
							,'Summer23Prompt23_RunCv1234_JRV1_DATA_SF_AK4PFPuppi.txt' #RunCv1234_JRV1_DATA_SF_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23Prompt23_RunCv1234_JRV1_DATA/Summer23Prompt23_RunCv1234_JRV1_DATA_SF_AK4PFchs.txt
            ]
        config.JobType.maxMemoryMB = 5000
        config.JobType.numCores = 4

        config.Data.inputDataset = None
        config.Data.splitting = 'FileBased'
        config.Data.unitsPerJob = 100
        # config.Data.totalUnits = 20 # problem with data/ MC if we want tto quick check, we can uncomment this
        config.Data.outLFNDirBase = '/store/user/pvaucell/MC_MUMU_2023_C_23_04_20254'
        config.Data.outputDatasetTag = 'MC'
        config.Data.inputDBS = 'global' #for signal : phys03 // for bkg :  global 
        config.Data.ignoreLocality = False
        config.Site.storageSite = 'T2_FR_IPHC' # Choose your site. 

        config.General.transferOutputs = True

        

        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDatasets = [
            		    '/TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/MINIAODSIM', # !! ok
						'/TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v4/MINIAODSIM', # !! ok
						 '/TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/MINIAODSIM', # !! ok'
                         '/DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM',# !! ok 
                         '/DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v2/MINIAODSIM', # !! ok 
                         '/WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v4/MINIAODSIM',# !! ok
                         '/WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3/MINIAODSIM',# !! ok 
                         '/ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v3/MINIAODSIM',# !! ok
						 '/TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2/MINIAODSIM', # !! ok
						 '/TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2/MINIAODSIM', # !! ok
						'/TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v15-v2/MINIAODSIM', # !! ok
                        '/TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8/Run3Summer23MiniAODv4-130X_mcRun3_2023_realistic_v14-v2/MINIAODSIM' # !! ok : TTToSemiLeptonic

			            ]


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
                                                                                                                      
