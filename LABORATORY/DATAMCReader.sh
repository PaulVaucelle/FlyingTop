#!/bin/bash

# Compilation
g++ -o DATAMCReader DATAMCReader.cpp $(root-config --cflags --libs) -lROOTDataFrame
if [[ $? -ne 0 ]]; then
    echo "Erreur lors de la compilation."
    # exit 1
fi

# Liste des fichiers d'entrée
input_files=(
            #  "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8" 
            #  "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8" 
            #  "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8" 
             "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8" 
            #  "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"
            #  "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8" 
            #  "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8" 
            #  "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8" 
            #  "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8" 
            #  "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8"
            #  "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8" 
            #  "TTWW_TuneCP5_13TeV-madgraph-pythia8"
            #     "ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8"
            #     "ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8"
            )

# The scripts works for all Ntuples exept the TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8 and TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8, because of the size i guess
# so this script is not really useful for the moment

    # ./DATAMCReader "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/MC_EMU_03_02_2025/MiniDATAMC_$infile.root" "./EMU/DATAMC_$infile.root"  MC_EMU_2018_03_02_2025  2018 0 0 0

# Exécution de MiniDataMCNtuple pour chaque fichier
for infile in "${input_files[@]}"; do
    echo "Processing: $infile -> ./EMU/$infile"
    ./DATAMCReader $infile $infile MC_EMU_2018_03_02_2025 2018 0 0 0
done
