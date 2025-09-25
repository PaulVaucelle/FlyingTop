#!/bin/bash

# List of work directories
# work_directories=("crab_20240924_143601" "crab_20240924_192448" "crab_20240924_192609"
#                   "crab_20240924_192705" "crab_20240924_194401" "crab_20240924_194432" "crab_20240924_194800" 
#                   "crab_20240924_194831" "crab_20240924_195014" "crab_20240924_195122" "crab_20240924_195145" 
#                   "crab_20240924_195223" "crab_20240924_195330" "crab_20240924_195407" "crab_20240924_195431" 
#                   "crab_20240924_195554" "crab_20240924_195657" "crab_20240924_195722" "crab_20240924_195828" 
#                   "crab_20240924_195852" "crab_20240924_200003" "crab_20240924_203034" "crab_20240924_203059" 
#                   "crab_20240924_203133" "crab_20240924_203204" "crab_20240924_203227"  )  # Add your directory names here
# work_directories=$(ls -d crab_*/)
work_directories=(
  "crab_20250204_133008" 
  "crab_20250205_090809" 

   )

# work_directories=("crab_20240924_143601" "crab_20240924_192448" "crab_20240924_192609" "crab_20240924_192705")
command=${1}
# Loop over each directory
for dir in "${work_directories[@]}"; do
  echo "Checking status for: $dir"
  crab $command -d ./"$dir"
  echo ""
done