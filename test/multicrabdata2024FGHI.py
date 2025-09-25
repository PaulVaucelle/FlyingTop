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
        config.General.workArea = "DATA_MUMU_2024_04_06_2025_FGHI"

        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'flyingtop_data_2024_FGHI.py'
        config.JobType.inputFiles = [
                'BDT_TRK_250527_2024_ctau100vsEMUdata.xml'
                ,'BDT_EVT_ALLSIGvsALLBKG.xml'
                ,'BDT_EVT_ALLSIGvsDYM50.xml'
                ,'BDT_EVT_ALLSIGvsTTTo2L2Nu.xml'
                ,'BDT_VTX_ALLSTEPS.xml'
                ,'BDT_VTX_STEP12.xml'
                ,'Run2024_MC.root'
                ,'Run2024_DATA.root'
                ,'pileupHistogram-Cert_Collisions2024_378891_386951_GoldenJson-13p6TeV-72400ub-100bins.root'
                ,'pileupHistogram-Cert_Collisions2024_378891_386951_GoldenJson-13p6TeV-66000ub-100bins.root'
                ,'RoccoR2018UL.txt'

                ,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt' #./JECDatabase/textFiles/Summer23BPixPrompt23_RunD_V1_DATA/Summer23BPixPrompt23_RunD_V1_DATA_Uncertainty_AK4PFchs.txt
                ,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt'#/JECDatabase/textFiles/Summer23BPixPrompt23_V1_MC/Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFchs.txt
                ,'Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_MC/Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFchs.txt
                ,'Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFPuppi.txt'#./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_MC/Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFchs.txt
                ,'Summer23BPixPrompt23_RunD_JRV1_DATA_PtResolution_AK4PFchs.txt' #./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_DATA/
                ,'Summer23BPixPrompt23_RunD_JRV1_DATA_SF_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_DATA/Summer23BPixPrompt23_RunD_JRV1_DATA_SF_AK4PFchs.txt
            ]
        config.JobType.maxMemoryMB = 5000
        config.JobType.numCores = 4

        config.Data.inputDataset = None
        config.Data.splitting = 'LumiBased'
        config.Data.unitsPerJob = 100
        #config.Data.unitsPerJob = 50
        # config.Data.totalUnits = 800

        config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions24/Cert_Collisions2024_378981_386951_Golden.json'
  
        config.Data.outLFNDirBase = '/store/user/pvaucell/DATA_MUMU_2024_04_06_2025_FGHI'
        config.Data.outputDatasetTag = 'data'
        config.Data.inputDBS = 'global' #for signal : phys03 // for bkg :  global 
        config.Data.ignoreLocality = False
        config.Site.storageSite = 'T2_FR_IPHC' # Choose your site. 
        config.Site.whitelist = ['T2_FR_IPHC']

        config.General.transferOutputs = True

        

        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDatasets = [
           
          
                # '/Muon0/Run2024C-2024CDEReprocessing-v1/MINIAOD',
                # '/Muon0/Run2024D-2024CDEReprocessing-v1/MINIAOD',
                # '/Muon0/Run2024E-2024CDEReprocessing-v1/MINIAOD',

                # '/Muon0/Run2024F-PromptReco-v1/MINIAOD',
                # '/Muon0/Run2024G-PromptReco-v1/MINIAOD',
                # '/Muon0/Run2024H-PromptReco-v1/MINIAOD',
                # '/Muon0/Run2024I-PromptReco-v1/MINIAOD',

                # '/Muon1/Run2024C-2024CDEReprocessing-v1/MINIAOD',
                # '/Muon1/Run2024D-2024CDEReprocessing-v1/MINIAOD',
                # '/Muon1/Run2024E-2024CDEReprocessing-v1/MINIAOD'

                # '/Muon1/Run2024F-PromptReco-v1/MINIAOD',
                '/Muon1/Run2024G-PromptReco-v1/MINIAOD'
                # '/Muon1/Run2024H-PromptReco-v1/MINIAOD',
                # '/Muon1/Run2024I-PromptReco-v2/MINIAOD'


                # '/MuonEG/Run2024C-2024CDEReprocessing-v1/MINIAOD',
                # '/MuonEG/Run2024D-2024CDEReprocessing-v1/MINIAOD',
                # '/MuonEG/Run2024E-2024CDEReprocessing-v1/MINIAOD',
                # '/MuonEG/Run2024F-PromptReco-v1/MINIAOD',
                # '/MuonEG/Run2024G-PromptReco-v1/MINIAOD',
                # '/MuonEG/Run2024H-PromptReco-v1/MINIAOD',
                # '/MuonEG/Run2024I-PromptReco-v1/MINIAOD',
                # '/MuonEG/Run2024I-PromptReco-v2/MINIAOD'
        ]

        for inDS in inputDatasets:
            # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            config.Data.publication = False
            # MuonDATA = inDS.split('/')[1:3]
            config.General.requestName = '_'.join(inDS.split('/')[1:3])
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
                                                                                                                      
