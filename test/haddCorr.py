import ROOT, sys, os, time, re, numpy, os.path
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog outputDir")
from tqdm import tqdm

# outputDir = sys.argv[1]
#date = name of working  directory

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

WorkingDir = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"

inputDatasets = [
# 'MC/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/19_10_2023_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/231019_213703/0000/',
# 'MC/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/19_10_2023_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/231019_213649/0000/',
# 'MC/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/19_10_2023_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/231019_213636/0000/',
# 'MC/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/19_10_2023_TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/231019_213601/0000/',
# 'MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/19_10_2023_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/231019_213623/0000/',
# 'MC/TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8/19_10_2023_TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8/231019_213900/0000/',
# 'MC/ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8/19_10_2023_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8/231019_213833/0000/',
# 'MC/TTWW_TuneCP5_13TeV-madgraph-pythia8/19_10_2023_TTWW_TuneCP5_13TeV-madgraph-pythia8/231019_213913/0000/', 
# 'MC/TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8/19_10_2023_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8/231019_213846/0000/',
# 'MC/WWTo2L2Nu_MLL_200To600_TuneCP5_13TeV-powheg-pythia8/19_10_2023_WWTo2L2Nu_MLL_200To600_TuneCP5_13TeV-powheg-pythia8/231019_213741/0000/',
# 'MC/WWTo2L2Nu_MLL_600To1200_TuneCP5_13TeV-powheg-pythia8/19_10_2023_WWTo2L2Nu_MLL_600To1200_TuneCP5_13TeV-powheg-pythia8/231019_213754/0000/',
# 'MC/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/19_10_2023_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/231019_213729/0000/',
# 'MC/WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/19_10_2023_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/231019_213807/0000/',
# 'MC/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/19_10_2023_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/231019_213820/0000/',

'MCHighStat/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/0000/',
'MCHighStat/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/0000/',
'MCHighStat/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/0000/',
'MCHighStat/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/0000/',
'MCHighStat/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/0000/',
'MCHighStat/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/0000/',
'MCHighStat/TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8/0000/',
'MCHighStat/ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8/0000/',
'MCHighStat/TTWW_TuneCP5_13TeV-madgraph-pythia8/0000/', 
'MCHighStat/TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8/0000/',
'MCHighStat/WWTo2L2Nu_MLL_200To600_TuneCP5_13TeV-powheg-pythia8/0000/',
'MCHighStat/WWTo2L2Nu_MLL_600To1200_TuneCP5_13TeV-powheg-pythia8/0000/',
'MCHighStat/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/0000/',
'MCHighStat/WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/0000/',
'MCHighStat/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/0000/'


# 'MC/RPV_2018_smu200_neu180_ctau001/19_10_2023_RPV_2018_smu200_neu180_ctau001/231019_210831/0000/',
# 'MC/RPV_2018_smu200_neu180_ctau100/19_10_2023_RPV_2018_smu200_neu180_ctau100/231019_211005/0000/',
# 'MC/RPV_2018_smu250_neu180_ctau100/19_10_2023_RPV_2018_smu250_neu180_ctau100/231019_211018/0000/',
# 'MC/RPV_2018_smu250_neu180_ctau001/19_10_2023_RPV_2018_smu250_neu180_ctau001/231019_210845/0000/',
# 'MC/RPV_2018_smu250_neu200_ctau001/19_10_2023_RPV_2018_smu250_neu200_ctau001/231019_210858/0000/',
# 'MC/RPV_2018_smu250_neu200_ctau100/19_10_2023_RPV_2018_smu250_neu200_ctau100/231019_211031/0000/',
# 'MC/RPV_2018_smu250_neu230_ctau001/19_10_2023_RPV_2018_smu250_neu230_ctau001/231019_210912/0000/',
# 'MC/RPV_2018_smu250_neu230_ctau100/19_10_2023_RPV_2018_smu250_neu230_ctau100/231019_211044/0000/',
# 'MC/RPV_2018_smu300_neu180_ctau100/19_10_2023_RPV_2018_smu300_neu180_ctau100/231019_211058/0000/',
# 'MC/RPV_2018_smu300_neu200_ctau100/19_10_2023_RPV_2018_smu300_neu200_ctau100/231019_211111/0000/',
# 'MC/RPV_2018_smu300_neu250_ctau100_v3/19_10_2023_RPV_2018_smu300_neu250_ctau100_v3/231019_211138/0000/',
# 'MC/RPV_2018_smu300_neu280_ctau100/19_10_2023_RPV_2018_smu300_neu280_ctau100/231019_211205/0000/',
# 'MC/RPV_2018_smu350_neu180_ctau100/19_10_2023_RPV_2018_smu350_neu180_ctau100/231019_211218/0000/',
# 'MC/RPV_2018_smu350_neu200_ctau100/19_10_2023_RPV_2018_smu350_neu200_ctau100/231019_211231/0000/',
# 'MC/RPV_2018_smu350_neu250_ctau100/19_10_2023_RPV_2018_smu350_neu250_ctau100/231019_211244/0000/',
# 'MC/RPV_2018_smu350_neu300_ctau100/19_10_2023_RPV_2018_smu350_neu300_ctau100/231019_211257/0000/',
# 'MC/RPV_2018_smu350_neu330_ctau100/19_10_2023_RPV_2018_smu350_neu330_ctau100/231019_211310/0000/',
# 'MC/RPV_2018_smu400_neu180_ctau100/19_10_2023_RPV_2018_smu400_neu180_ctau100/231019_211323/0000/',
# 'MC/RPV_2018_smu400_neu200_ctau100/19_10_2023_RPV_2018_smu400_neu200_ctau100/231019_211336/0000/',
# 'MC/RPV_2018_smu400_neu250_ctau100/19_10_2023_RPV_2018_smu400_neu250_ctau100/231019_211349/0000/',
# 'MC/RPV_2018_smu400_neu300_ctau100/19_10_2023_RPV_2018_smu400_neu300_ctau100/231019_211402/0000/',
# 'MC/RPV_2018_smu400_neu350_ctau100/19_10_2023_RPV_2018_smu400_neu350_ctau100/231019_211415/0000/',
# 'MC/RPV_2018_smu400_neu380_ctau100/19_10_2023_RPV_2018_smu400_neu380_ctau100/231019_211428/0000/',
# 'MC/RPV_2018_smu450_neu180_ctau100/19_10_2023_RPV_2018_smu450_neu180_ctau100/231019_211442/0000/',
# 'MC/RPV_2018_smu450_neu200_ctau100/19_10_2023_RPV_2018_smu450_neu200_ctau100/231019_211455/0000/',
# 'MC/RPV_2018_smu450_neu250_ctau100/19_10_2023_RPV_2018_smu450_neu250_ctau100/231019_211508/0000/',
# 'MC/RPV_2018_smu450_neu300_ctau100/19_10_2023_RPV_2018_smu450_neu300_ctau100/231019_211521/0000/',
# 'MC/RPV_2018_smu450_neu350_ctau100/19_10_2023_RPV_2018_smu450_neu350_ctau100/231019_211534/0000/',
# 'MC/RPV_2018_smu450_neu400_ctau100/19_10_2023_RPV_2018_smu450_neu400_ctau100/231019_211548/0000/',
# 'MC/RPV_2018_smu450_neu430_ctau100/19_10_2023_RPV_2018_smu450_neu430_ctau100/231019_211601/0000/',
# 'MC/RPV_2018_smu500_neu180_ctau100/19_10_2023_RPV_2018_smu500_neu180_ctau100/231019_211614/0000/',
# 'MC/RPV_2018_smu500_neu200_ctau100/19_10_2023_RPV_2018_smu500_neu200_ctau100/231019_211627/0000/',
# 'MC/RPV_2018_smu500_neu250_ctau100/19_10_2023_RPV_2018_smu500_neu250_ctau100/231019_211640/0000/',
# 'MC/RPV_2018_smu500_neu300_ctau100/19_10_2023_RPV_2018_smu500_neu300_ctau100/231019_211653/0000/',
# 'MC/RPV_2018_smu500_neu350_ctau100/19_10_2023_RPV_2018_smu500_neu350_ctau100/231019_211707/0000/',
# 'MC/RPV_2018_smu500_neu400_ctau100/19_10_2023_RPV_2018_smu500_neu400_ctau100/231019_211720/0000/',
# 'MC/RPV_2018_smu500_neu450_ctau100/19_10_2023_RPV_2018_smu500_neu450_ctau100/231019_211733/0000/',
# 'MC/RPV_2018_smu500_neu480_ctau100/19_10_2023_RPV_2018_smu500_neu480_ctau100/231019_211747/0000/'
]


BackgroundSamples = [
WorkingDir+inputDatasets[0],
WorkingDir+inputDatasets[1],
WorkingDir+inputDatasets[2],
WorkingDir+inputDatasets[3],
WorkingDir+inputDatasets[4],
WorkingDir+inputDatasets[5],
WorkingDir+inputDatasets[6],
WorkingDir+inputDatasets[7],
WorkingDir+inputDatasets[8],
WorkingDir+inputDatasets[9],
WorkingDir+inputDatasets[10],
WorkingDir+inputDatasets[11],
WorkingDir+inputDatasets[12],
WorkingDir+inputDatasets[13],
WorkingDir+inputDatasets[14],
# WorkingDir+inputDatasets[15],
# WorkingDir+inputDatasets[16],
# WorkingDir+inputDatasets[17],
# WorkingDir+inputDatasets[18],
# WorkingDir+inputDatasets[19],
# WorkingDir+inputDatasets[20],
# WorkingDir+inputDatasets[21],
# WorkingDir+inputDatasets[22],
# WorkingDir+inputDatasets[23],
# WorkingDir+inputDatasets[24],
# WorkingDir+inputDatasets[25],
# WorkingDir+inputDatasets[26],
# WorkingDir+inputDatasets[27],
# WorkingDir+inputDatasets[28],
# WorkingDir+inputDatasets[29],
# WorkingDir+inputDatasets[30],
# WorkingDir+inputDatasets[31],
# WorkingDir+inputDatasets[32],
# WorkingDir+inputDatasets[33],
# WorkingDir+inputDatasets[34],
# WorkingDir+inputDatasets[35],
# WorkingDir+inputDatasets[36],
# WorkingDir+inputDatasets[37],
# WorkingDir+inputDatasets[38],
# WorkingDir+inputDatasets[39],
# WorkingDir+inputDatasets[40],
# WorkingDir+inputDatasets[41],
# WorkingDir+inputDatasets[42],
# WorkingDir+inputDatasets[43],
# WorkingDir+inputDatasets[44],
# WorkingDir+inputDatasets[45],
# WorkingDir+inputDatasets[46],
# WorkingDir+inputDatasets[47],
# WorkingDir+inputDatasets[48],
# WorkingDir+inputDatasets[49],
# WorkingDir+inputDatasets[50],
# WorkingDir+inputDatasets[51],

]


for inDS in BackgroundSamples:
    dir1 = inDS.split('/')[12]
    dir2  = dir1.split('_Tune')[0]
    # print('Filename : ',dir1)
    print('Filename : ',dir2)

    # dir3 = dir2.split('_',1)[1]#retrieve the name of the background
    # print('Filename : ',dir3)
    os.chdir(inDS)
    os.system("hadd -f "+dir2+".root Ntuple_*.root") #use dir2 : fro background samples ,dir3 can be useful in some cases :os.system("hadd -f Ntuple_"+dir2+".root Ntuple_*.root")
    os.system("mv "+dir2+".root ../../..") #os.system("mv Ntuple_"+dir2+".root ../../../../..")
    # fileInArray.append(ROOT.TFile.Open("Ntuple_"+dir3+".root","UPDATE"))# update the root file

# os.chdir("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/")




