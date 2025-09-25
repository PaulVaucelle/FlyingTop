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
        config.General.workArea = "DATA_MUMU_2023_D_23_04_2025_v2"

        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'flyingtop_data_2023_D.py'
        config.JobType.inputFiles =  [
                'BDT_TRK_241119_2023B_ctau100vsEMUdata.xml'
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
                ,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt' #./JECDatabase/textFiles/Summer23BPixPrompt23_RunD_V1_DATA/Summer23BPixPrompt23_RunD_V1_DATA_Uncertainty_AK4PFchs.txt
                ,'Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFPuppi.txt'#/JECDatabase/textFiles/Summer23BPixPrompt23_V1_MC/Summer23BPixPrompt23_V1_MC_Uncertainty_AK4PFchs.txt
                ,'Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_MC/Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFchs.txt
                ,'Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFPuppi.txt'#./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_MC/Summer23BPixPrompt23_RunD_JRV1_MC_SF_AK4PFchs.txt
                ,'Summer23BPixPrompt23_RunD_JRV1_DATA_PtResolution_AK4PFchs.txt' #./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_DATA/
                ,'Summer23BPixPrompt23_RunD_JRV1_DATA_SF_AK4PFPuppi.txt' #./JRDatabase/textFiles/Summer23BPixPrompt23_RunD_JRV1_DATA/Summer23BPixPrompt23_RunD_JRV1_DATA_SF_AK4PFchs.txt
        ]

        config.JobType.maxMemoryMB = 4000
        config.JobType.numCores = 4

        config.Data.inputDataset = None
        config.Data.splitting = 'LumiBased'
        config.Data.unitsPerJob = 100
        #config.Data.unitsPerJob = 50
        # config.Data.totalUnits = 800

        config.Data.lumiMask = 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Golden.json'

        #config.Data.lumiMask = '/opt/sbg/cms/ui2_data1/mmeena/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Crab_23_11_23/crab_MuonEG/results/notFinishedLumis.txt'
        config.Data.outLFNDirBase = '/store/user/pvaucell/DATA_MUMU_D_2023_23_04_2025_v2'
        config.Data.outputDatasetTag = 'data'
        config.Data.inputDBS = 'global' #for signal : phys03 // for bkg :  global 
        config.Data.ignoreLocality = False
        config.Site.storageSite = 'T2_FR_IPHC' # Choose your site. 
        config.General.transferOutputs = True
        
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDatasets = [
            # '/Muon0/Run2023D-22Sep2023_v1-v1/MINIAOD',
            # '/Muon0/Run2023D-22Sep2023_v2-v1/MINIAOD',
            # '/Muon1/Run2023D-22Sep2023_v1-v1/MINIAOD'
            # '/Muon1/Run2023D-22Sep2023_v2-v1/MINIAOD'

            # '/MuonEG/Run2023D-22Sep2023_v1-v1/MINIAOD',
            # '/MuonEG/Run2023D-22Sep2023_v2-v1/MINIAOD'

        ]

        for inDS in inputDatasets:
            # inDS is of the form /A/B/C. Since B is unique for each inDS, use this in the CRAB request name.
            config.Data.publication = False
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
                                                                                                                      
