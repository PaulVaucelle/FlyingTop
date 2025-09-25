import os

list_arguments = [
  "Msmu-400_Mchi-250_ct-100"
 ,"Msmu-450_Mchi-200_ct-100"
]

year = "2023"
ispost = "True"
#$$ chsalgos = "_chsalgos_106X_UL18"

for i in range(len(list_arguments)):
    command = "rm flyingtop_temp0.py"
    os.system(command) 
    command = "rm flyingtop_temp1.py"
    os.system(command) 
#$$     command = "rm flyingtop_temp2.py"
#$$     os.system(command) 
    command = "rm flyingtop_temp.py"
    os.system(command) 
    command = "rm crab_config_temp.py"
    os.system(command) 
    sedFile0 = "sed -e 's#Year#" + year+"#' flyingtop_RPV_default.py > flyingtop_temp0.py"
    os.system(sedFile0) 
    sedFile1 = "sed -e 's#isPost#" + ispost+"#' flyingtop_temp0.py > flyingtop_temp1.py"
    os.system(sedFile1) 
#$$    sedFile2 = "sed -e 's#CHSALGOS#" + chsalgos+"#' flyingtop_temp1.py > flyingtop_temp2.py"
#$$    os.system(sedFile2) 
#$$    sedFileZ = "sed -e 's#inputSample#" + list_arguments[i]+"#' flyingtop_temp2.py > flyingtop_temp.py"
    sedFileZ = "sed -e 's#inputSample#" + list_arguments[i]+"#' flyingtop_temp1.py > flyingtop_temp.py"
    os.system(sedFileZ) 
    sedFile = "sed -e 's#inputSample#" + list_arguments[i]+"#' crab_config_RPV_2023B_default_v3.py > crab_config_temp.py"
    os.system(sedFile) 
    command = "crab-pre submit -c crab_config_temp.py"
    os.system(command) 

