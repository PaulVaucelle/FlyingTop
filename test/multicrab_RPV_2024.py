import os

list_arguments = [
     "Par-ct-001-MChi-180-MSmu-200",
"Par-ct-001-MChi-180-MSmu-250", "Par-ct-001-MChi-200-MSmu-250"
,"Par-ct-001-MChi-180-MSmu-300", "Par-ct-001-MChi-200-MSmu-300", "Par-ct-001-MChi-250-MSmu-300", "Par-ct-001-MChi-280-MSmu-300" 
,"Par-ct-001-MChi-180-MSmu-350", "Par-ct-001-MChi-200-MSmu-350", "Par-ct-001-MChi-250-MSmu-350", "Par-ct-001-MChi-300-MSmu-350", "Par-ct-001-MChi-330-MSmu-350" 
,"Par-ct-001-MChi-180-MSmu-400", "Par-ct-001-MChi-200-MSmu-400", "Par-ct-001-MChi-250-MSmu-400", "Par-ct-001-MChi-300-MSmu-400", "Par-ct-001-MChi-350-MSmu-400", "Par-ct-001-MChi-380-MSmu-400" 
,"Par-ct-001-MChi-180-MSmu-450", "Par-ct-001-MChi-200-MSmu-450", "Par-ct-001-MChi-250-MSmu-450", "Par-ct-001-MChi-300-MSmu-450", "Par-ct-001-MChi-350-MSmu-450", "Par-ct-001-MChi-400-MSmu-450", "Par-ct-001-MChi-430-MSmu-450" 
,"Par-ct-001-MChi-180-MSmu-500", "Par-ct-001-MChi-200-MSmu-500", "Par-ct-001-MChi-250-MSmu-500", "Par-ct-001-MChi-300-MSmu-500", "Par-ct-001-MChi-350-MSmu-500", "Par-ct-001-MChi-400-MSmu-500", "Par-ct-001-MChi-450-MSmu-500", 
"Par-ct-001-MChi-480-MSmu-500",

 "Par-ct-003-MChi-180-MSmu-200",
"Par-ct-003-MChi-180-MSmu-250", "Par-ct-003-MChi-200-MSmu-250"
,"Par-ct-003-MChi-180-MSmu-300", "Par-ct-003-MChi-200-MSmu-300", "Par-ct-003-MChi-250-MSmu-300", "Par-ct-003-MChi-280-MSmu-300" 
,"Par-ct-003-MChi-180-MSmu-350", "Par-ct-003-MChi-200-MSmu-350", "Par-ct-003-MChi-250-MSmu-350", "Par-ct-003-MChi-300-MSmu-350", "Par-ct-003-MChi-330-MSmu-350" 
,"Par-ct-003-MChi-180-MSmu-400", "Par-ct-003-MChi-200-MSmu-400", "Par-ct-003-MChi-250-MSmu-400", "Par-ct-003-MChi-300-MSmu-400", "Par-ct-003-MChi-350-MSmu-400", "Par-ct-003-MChi-380-MSmu-400" 
,"Par-ct-003-MChi-180-MSmu-450", "Par-ct-003-MChi-200-MSmu-450", "Par-ct-003-MChi-250-MSmu-450", "Par-ct-003-MChi-300-MSmu-450", "Par-ct-003-MChi-350-MSmu-450", "Par-ct-003-MChi-400-MSmu-450", "Par-ct-003-MChi-430-MSmu-450" 
,"Par-ct-003-MChi-180-MSmu-500", "Par-ct-003-MChi-200-MSmu-500", "Par-ct-003-MChi-250-MSmu-500", "Par-ct-003-MChi-300-MSmu-500", "Par-ct-003-MChi-350-MSmu-500", "Par-ct-003-MChi-400-MSmu-500", "Par-ct-003-MChi-450-MSmu-500", 
"Par-ct-003-MChi-480-MSmu-500" ,

 "Par-ct-010-MChi-180-MSmu-200",
"Par-ct-010-MChi-180-MSmu-250", "Par-ct-010-MChi-200-MSmu-250"
,"Par-ct-010-MChi-180-MSmu-300", "Par-ct-010-MChi-200-MSmu-300", "Par-ct-010-MChi-250-MSmu-300", "Par-ct-010-MChi-280-MSmu-300" 
,"Par-ct-010-MChi-180-MSmu-350", "Par-ct-010-MChi-200-MSmu-350", "Par-ct-010-MChi-250-MSmu-350", "Par-ct-010-MChi-300-MSmu-350", "Par-ct-010-MChi-330-MSmu-350" 
,"Par-ct-010-MChi-180-MSmu-400", "Par-ct-010-MChi-200-MSmu-400", "Par-ct-010-MChi-250-MSmu-400", "Par-ct-010-MChi-300-MSmu-400", "Par-ct-010-MChi-350-MSmu-400", "Par-ct-010-MChi-380-MSmu-400" 
,"Par-ct-010-MChi-180-MSmu-450", "Par-ct-010-MChi-200-MSmu-450", "Par-ct-010-MChi-250-MSmu-450", "Par-ct-010-MChi-300-MSmu-450", "Par-ct-010-MChi-350-MSmu-450", "Par-ct-010-MChi-400-MSmu-450", "Par-ct-010-MChi-430-MSmu-450" 
,"Par-ct-010-MChi-180-MSmu-500", "Par-ct-010-MChi-200-MSmu-500", "Par-ct-010-MChi-250-MSmu-500", "Par-ct-010-MChi-300-MSmu-500", "Par-ct-010-MChi-350-MSmu-500", "Par-ct-010-MChi-400-MSmu-500", "Par-ct-010-MChi-450-MSmu-500", 
"Par-ct-010-MChi-480-MSmu-500" ,

 "Par-ct-030-MChi-180-MSmu-200",
"Par-ct-030-MChi-180-MSmu-250", "Par-ct-030-MChi-200-MSmu-250"
,"Par-ct-030-MChi-180-MSmu-300", "Par-ct-030-MChi-200-MSmu-300", "Par-ct-030-MChi-250-MSmu-300", "Par-ct-030-MChi-280-MSmu-300" 
,"Par-ct-030-MChi-180-MSmu-350", "Par-ct-030-MChi-200-MSmu-350", "Par-ct-030-MChi-250-MSmu-350", "Par-ct-030-MChi-300-MSmu-350", "Par-ct-030-MChi-330-MSmu-350" 
,"Par-ct-030-MChi-180-MSmu-400", "Par-ct-030-MChi-200-MSmu-400", "Par-ct-030-MChi-250-MSmu-400", "Par-ct-030-MChi-300-MSmu-400", "Par-ct-030-MChi-350-MSmu-400", "Par-ct-030-MChi-380-MSmu-400" 
,"Par-ct-030-MChi-180-MSmu-450", "Par-ct-030-MChi-200-MSmu-450", "Par-ct-030-MChi-250-MSmu-450", "Par-ct-030-MChi-300-MSmu-450", "Par-ct-030-MChi-350-MSmu-450", "Par-ct-030-MChi-400-MSmu-450", "Par-ct-030-MChi-430-MSmu-450" 
,"Par-ct-030-MChi-180-MSmu-500", "Par-ct-030-MChi-200-MSmu-500", "Par-ct-030-MChi-250-MSmu-500", "Par-ct-030-MChi-300-MSmu-500", "Par-ct-030-MChi-350-MSmu-500", "Par-ct-030-MChi-400-MSmu-500", "Par-ct-030-MChi-450-MSmu-500", 
"Par-ct-030-MChi-480-MSmu-500" ,


 "Par-ct-100-MChi-180-MSmu-200",
"Par-ct-100-MChi-180-MSmu-250", "Par-ct-100-MChi-200-MSmu-250"
,"Par-ct-100-MChi-180-MSmu-300", "Par-ct-100-MChi-200-MSmu-300", "Par-ct-100-MChi-250-MSmu-300", "Par-ct-100-MChi-280-MSmu-300" 
,"Par-ct-100-MChi-180-MSmu-350", "Par-ct-100-MChi-200-MSmu-350", "Par-ct-100-MChi-250-MSmu-350", "Par-ct-100-MChi-300-MSmu-350", "Par-ct-100-MChi-330-MSmu-350" 
,"Par-ct-100-MChi-180-MSmu-400", "Par-ct-100-MChi-200-MSmu-400", "Par-ct-100-MChi-250-MSmu-400", "Par-ct-100-MChi-300-MSmu-400", "Par-ct-100-MChi-350-MSmu-400", "Par-ct-100-MChi-380-MSmu-400" 
,"Par-ct-100-MChi-180-MSmu-450", "Par-ct-100-MChi-200-MSmu-450", "Par-ct-100-MChi-250-MSmu-450", "Par-ct-100-MChi-300-MSmu-450", "Par-ct-100-MChi-350-MSmu-450", "Par-ct-100-MChi-400-MSmu-450", "Par-ct-100-MChi-430-MSmu-450" 
,"Par-ct-100-MChi-180-MSmu-500", "Par-ct-100-MChi-200-MSmu-500", "Par-ct-100-MChi-250-MSmu-500", "Par-ct-100-MChi-300-MSmu-500", "Par-ct-100-MChi-350-MSmu-500", "Par-ct-100-MChi-400-MSmu-500", "Par-ct-100-MChi-450-MSmu-500",
 "Par-ct-100-MChi-480-MSmu-500" ,

 "Par-ct-300-MChi-180-MSmu-200",
"Par-ct-300-MChi-180-MSmu-250", "Par-ct-300-MChi-200-MSmu-250"
,"Par-ct-300-MChi-180-MSmu-300", "Par-ct-300-MChi-200-MSmu-300", "Par-ct-300-MChi-250-MSmu-300", "Par-ct-300-MChi-280-MSmu-300" 
,"Par-ct-300-MChi-180-MSmu-350", "Par-ct-300-MChi-200-MSmu-350", "Par-ct-300-MChi-250-MSmu-350", "Par-ct-300-MChi-300-MSmu-350", "Par-ct-300-MChi-330-MSmu-350" 
,"Par-ct-300-MChi-180-MSmu-400", "Par-ct-300-MChi-200-MSmu-400", "Par-ct-300-MChi-250-MSmu-400", "Par-ct-300-MChi-300-MSmu-400", "Par-ct-300-MChi-350-MSmu-400", "Par-ct-300-MChi-380-MSmu-400" 
,"Par-ct-300-MChi-180-MSmu-450", "Par-ct-300-MChi-200-MSmu-450", "Par-ct-300-MChi-250-MSmu-450", "Par-ct-300-MChi-300-MSmu-450", "Par-ct-300-MChi-350-MSmu-450", "Par-ct-300-MChi-400-MSmu-450", "Par-ct-300-MChi-430-MSmu-450" 
,"Par-ct-300-MChi-180-MSmu-500", "Par-ct-300-MChi-200-MSmu-500", "Par-ct-300-MChi-250-MSmu-500", "Par-ct-300-MChi-300-MSmu-500", "Par-ct-300-MChi-350-MSmu-500", "Par-ct-300-MChi-400-MSmu-500", "Par-ct-300-MChi-450-MSmu-500", 
"Par-ct-300-MChi-480-MSmu-500" ,


 "Par-ct-1000-MChi-180-MSmu-200",
"Par-ct-1000-MChi-180-MSmu-250", "Par-ct-1000-MChi-200-MSmu-250"
,"Par-ct-1000-MChi-180-MSmu-300", "Par-ct-1000-MChi-200-MSmu-300", "Par-ct-1000-MChi-250-MSmu-300", "Par-ct-1000-MChi-280-MSmu-300" 
,"Par-ct-1000-MChi-180-MSmu-350", "Par-ct-1000-MChi-200-MSmu-350", "Par-ct-1000-MChi-250-MSmu-350", "Par-ct-1000-MChi-300-MSmu-350", "Par-ct-1000-MChi-330-MSmu-350" 
,"Par-ct-1000-MChi-180-MSmu-400", "Par-ct-1000-MChi-200-MSmu-400", "Par-ct-1000-MChi-250-MSmu-400", "Par-ct-1000-MChi-300-MSmu-400", "Par-ct-1000-MChi-350-MSmu-400", "Par-ct-1000-MChi-380-MSmu-400" 
,"Par-ct-1000-MChi-180-MSmu-450", "Par-ct-1000-MChi-200-MSmu-450", "Par-ct-1000-MChi-250-MSmu-450", "Par-ct-1000-MChi-300-MSmu-450", "Par-ct-1000-MChi-350-MSmu-450", "Par-ct-1000-MChi-400-MSmu-450", "Par-ct-1000-MChi-430-MSmu-450" 
,"Par-ct-1000-MChi-180-MSmu-500", "Par-ct-1000-MChi-200-MSmu-500", "Par-ct-1000-MChi-250-MSmu-500", "Par-ct-1000-MChi-300-MSmu-500", "Par-ct-1000-MChi-350-MSmu-500", "Par-ct-1000-MChi-400-MSmu-500", "Par-ct-1000-MChi-450-MSmu-500", 
"Par-ct-1000-MChi-480-MSmu-500" 


]
year = "2024"
ispost = "False"
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
    sedFile0 = "sed -e 's#Year#" + year+"#' flyingtop_default_RPV2024.py > flyingtop_temp0.py"
    os.system(sedFile0) 
    sedFile1 = "sed -e 's#isPost#" + ispost+"#' flyingtop_temp0.py > flyingtop_temp1.py"
    os.system(sedFile1) 
#$$    sedFile2 = "sed -e 's#CHSALGOS#" + chsalgos+"#' flyingtop_temp1.py > flyingtop_temp2.py"
#$$    os.system(sedFile2) 
#$$    sedFileZ = "sed -e 's#inputSample#" + list_arguments[i]+"#' flyingtop_temp2.py > flyingtop_temp.py"
    sedFileZ = "sed -e 's#inputSample#" + list_arguments[i]+"#' flyingtop_temp1.py > flyingtop_temp.py"
    os.system(sedFileZ) 
    sedFile = "sed -e 's#inputSample#" + list_arguments[i]+"#' crab_config_RPV_2024_default.py > crab_config_temp.py"
    os.system(sedFile) 
    command = "crab-pre submit -c crab_config_temp.py"
    os.system(command) 
