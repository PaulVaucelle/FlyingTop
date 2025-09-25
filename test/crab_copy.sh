#!/bin/bash

# List of work directories

work_directories=(


# "MC_EMU_2022_23_04_2025"
# "MC_EMU_2022_EFG_23_04_2025"
# "DATA_EMU_2024_17_05_2025"

# "DATA_MUMU_2024_04_06_2025"
# "DATA_MUMU_2024_16_06_2025"
MC/RPV_2022/250428/
 )

# work_directories=("DATA_MUMU_2022_25_09_2024/Muon")

# Loop over each directory
for dir in "${work_directories[@]}"; do
  echo "Copying: $dir"
  mkdir -p $dir
  cd $dir
  gfal-copy -r --timeout 0 -p davs://sbgdcache.in2p3.fr/cms/phedex/store/user/pvaucell/$dir ./
  echo "Ended copy of $dir"
  cd ..
done



