import ROOT, sys, os, time, re, numpy, os.path
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog outputDir")
from tqdm import tqdm

# outputDir = sys.argv[1]
#date = name of working  directory

WorkingDir = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/17_04_2024"

Samples =[
    'DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8',
    'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8',
    'TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8',
    'TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8',
    'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8'
]

inputDatasets = [
"/240418_123823/0000",
"/240418_123839/0000",
"/240418_123737/0000",
"/240418_123955/0000",
"/240418_124041/0000"
]


BackgroundSamples = [
 WorkingDir+inputDatasets[0],
 WorkingDir+inputDatasets[1],
 WorkingDir+inputDatasets[2],
 WorkingDir+inputDatasets[3],
 WorkingDir+inputDatasets[4]
]

count = 0
for inDS in BackgroundSamples:
    # dir1 = inDS.split('/')[12]
    # dir2  = dir1.split('_Tune')[0]
    # # print('Filename : ',dir1)
    # print('Filename : ',dir2)

    # dir3 = dir2.split('_',1)[1]#retrieve the name of the background
    # print('Filename : ',dir3)
    for i in range(50):
        os.chdir(inDS)
        os.system("hadd -f "+Samples[count]+"_"+str(i)+".root NtupleS_Testavec_"+str(i)+"?.root") #use dir2 : fro background samples ,dir3 can be useful in some cases :os.system("hadd -f Ntuple_"+dir2+".root Ntuple_*.root")
        os.system("mv "+Samples[count]+"_"+str(i)+".root ../..") #os.system("mv Ntuple_"+dir2+".root ../../../../..")
    #os.system("mv "+inDS+".root "+Samples[count]) #os.system("mv Ntuple_"+dir2+".root ../../../../..")

    # fileInArray.append(ROOT.TFile.Open("Ntuple_"+dir3+".root","UPDATE"))# update the root file
    count= count +1
# os.chdir("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/")




