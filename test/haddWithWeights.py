import ROOT, sys, os, time, re, numpy, os.path
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog outputDir")
from tqdm import tqdm

outputDir = sys.argv[1]
#date = name of working  directory

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

WorkingDir = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/"

inputDatasets = [
'/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',
'/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM',
'/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM',
'/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM',
'/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM',
'/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',
'/WWTo2L2Nu_MLL_200To600_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1/MINIAODSIM',
'/WWTo2L2Nu_MLL_600To1200_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v3/MINIAODSIM',
'/WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',
'/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',
]

BackgroundSamples = [
WorkingDir+outputDir+"crab_"+inputDatasets[0].split('/')[1]+"/results/",
WorkingDir+outputDir+"crab_"+inputDatasets[1].split('/')[1]+"/results/",
WorkingDir+outputDir+"crab_"+inputDatasets[2].split('/')[1]+"/results/",
WorkingDir+outputDir+"crab_"+inputDatasets[3].split('/')[1]+"/results/",
WorkingDir+outputDir+"crab_"+inputDatasets[4].split('/')[1]+"/results/",
WorkingDir+outputDir+"crab_"+inputDatasets[5].split('/')[1]+"/results/",
WorkingDir+outputDir+"crab_"+inputDatasets[6].split('/')[1]+"/results/",
WorkingDir+outputDir+"crab_"+inputDatasets[7].split('/')[1]+"/results/",
WorkingDir+outputDir+"crab_"+inputDatasets[8].split('/')[1]+"/results/",
WorkingDir+outputDir+"crab_"+inputDatasets[9].split('/')[1]+"/results/"
]

intLumi = 10.0 # 137.0

count = 0
fileInArray = []

for inDS in BackgroundSamples:
    # for path in os.listdir(inDS):
    #     # check if current path is a file
    #     if os.path.isfile(os.path.join(inDS, path)):
    #         if path.endswith(".root"):
    #             print('Filename : ',path)
    dir1 = inDS.split('/')[12]
    dir2  = dir1.split('_Tune')[0]
    dir3 = dir2.split('_',1)[1]#retrieve the name of the background
    os.chdir(inDS)
    os.system("hadd -f Ntuple_"+dir3+".root Ntuple_*.root")
    os.system("mv Ntuple_"+dir3+".root ../../..")
    # fileInArray.append(ROOT.TFile.Open("Ntuple_"+dir3+".root","UPDATE"))# update the root file

os.chdir("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/")
os.system("mv Ntuple_DYJetsToLL_M-10to50.root Ntuple_DYJetsToLL_M10to50.root")
os.system("mv Ntuple_DYJetsToLL_M-50.root Ntuple_DYJetsToLL_M50.root")



