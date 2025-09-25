import ROOT
import os

list_arguments = [

## 3 hours - 3h30min per benchmarck of ctau
                      
                        "RPV_2022A_smu500_neu200_ctau001",
                      
                        "RPV_2022A_smu500_neu200_ctau003",
                       
                        "RPV_2022A_smu500_neu200_ctau010",
                       
                        "RPV_2022A_smu500_neu200_ctau030",
                        
                        "RPV_2022A_smu500_neu200_ctau100"
                       
                        
]

for i in range(len(list_arguments)):
    command = "rm flyingtop_temp.py"
    os.system(command) 
    sedFile = "sed -e 's#inputFile#" + list_arguments[i]+"#' flyingtop_default_2022.py > flyingtop_temp.py"
    os.system(sedFile) 
    command = "cmsRun flyingtop_temp.py"
    os.system(command) 

