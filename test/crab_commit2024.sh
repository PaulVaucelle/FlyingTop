#!/bin/bash


Sample=("/Muon0/Run2024A-PromptReco-v1/MINIAOD" "/Muon0/Run2024B-PromptReco-v1/MINIAOD" "/Muon0/Run2024C-PromptReco-v1/MINIAOD" "/Muon0/Run2024D-PromptReco-v1/MINIAOD" 
"/Muon0/Run2024E-PromptReco-v2/MINIAOD" "/Muon0/Run2024F-PromptReco-v1/MINIAOD" "/Muon0/Run2024G-PromptReco-v1/MINIAOD" "/Muon0/Run2024H-PromptReco-v1/MINIAOD"
"/Muon1/Run2024A-PromptReco-v1/MINIAOD" "/Muon1/Run2024B-PromptReco-v1/MINIAOD" "/Muon1/Run2024C-PromptReco-v1/MINIAOD" "/Muon1/Run2024D-PromptReco-v1/MINIAOD"
"/Muon1/Run2024E-PromptReco-v2/MINIAOD" "/Muon1/Run2024F-PromptReco-v1/MINIAOD" "/Muon1/Run2024G-PromptReco-v1/MINIAOD" "/Muon1/Run2024H-PromptReco-v1/MINIAOD"
  )

command=submit 
submit_file=crab_config_data_2024.py
# Loop over each directory
for dir in "${Sample[@]}"; do
  echo "Sumitting  for: $dir"
  python ChangeDataset2024.py "${Sample[@]}" '/store/user/pvaucell/DATA_MUMU_2024_25_09_2024'
  crab $command -c ./"$submit_file"
  echo ""
done
