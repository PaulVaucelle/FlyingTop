import os

list_arguments = [
  "Msmu-200_Mchi-180_ct-001"
 "Msmu-250_Mchi-180_ct-001","Msmu-250_Mchi-200_ct-001"
 ,"Msmu-300_Mchi-180_ct-001","Msmu-300_Mchi-200_ct-001","Msmu-300_Mchi-250_ct-001","Msmu-300_Mchi-280_ct-001"
 ,"Msmu-350_Mchi-180_ct-001","Msmu-350_Mchi-200_ct-001","Msmu-350_Mchi-250_ct-001","Msmu-350_Mchi-300_ct-001","Msmu-350_Mchi-330_ct-001"
 ,"Msmu-400_Mchi-180_ct-001","Msmu-400_Mchi-200_ct-001","Msmu-400_Mchi-250_ct-001","Msmu-400_Mchi-300_ct-001","Msmu-400_Mchi-350_ct-001","Msmu-400_Mchi-380_ct-001"
 ,"Msmu-450_Mchi-180_ct-001","Msmu-450_Mchi-200_ct-001","Msmu-450_Mchi-250_ct-001","Msmu-450_Mchi-300_ct-001","Msmu-450_Mchi-350_ct-001","Msmu-450_Mchi-400_ct-001","Msmu-450_Mchi-430_ct-001"
 ,"Msmu-500_Mchi-180_ct-001","Msmu-500_Mchi-200_ct-001","Msmu-500_Mchi-250_ct-001","Msmu-500_Mchi-300_ct-001","Msmu-500_Mchi-350_ct-001","Msmu-500_Mchi-400_ct-001","Msmu-500_Mchi-450_ct-001","Msmu-500_Mchi-480_ct-001"
 ,"Msmu-200_Mchi-180_ct-003"
 ,"Msmu-250_Mchi-180_ct-003","Msmu-250_Mchi-200_ct-003"
 ,"Msmu-300_Mchi-180_ct-003","Msmu-300_Mchi-200_ct-003","Msmu-300_Mchi-250_ct-003","Msmu-300_Mchi-280_ct-003"
 ,"Msmu-350_Mchi-180_ct-003","Msmu-350_Mchi-200_ct-003","Msmu-350_Mchi-250_ct-003","Msmu-350_Mchi-300_ct-003","Msmu-350_Mchi-330_ct-003"
 ,"Msmu-400_Mchi-180_ct-003","Msmu-400_Mchi-200_ct-003","Msmu-400_Mchi-250_ct-003","Msmu-400_Mchi-300_ct-003",
 "Msmu-400_Mchi-350_ct-003","Msmu-400_Mchi-380_ct-003"
 ,"Msmu-450_Mchi-180_ct-003","Msmu-450_Mchi-200_ct-003","Msmu-450_Mchi-250_ct-003","Msmu-450_Mchi-300_ct-003","Msmu-450_Mchi-350_ct-003","Msmu-450_Mchi-400_ct-003","Msmu-450_Mchi-430_ct-003"
 ,"Msmu-500_Mchi-180_ct-003","Msmu-500_Mchi-200_ct-003","Msmu-500_Mchi-250_ct-003","Msmu-500_Mchi-300_ct-003","Msmu-500_Mchi-350_ct-003","Msmu-500_Mchi-400_ct-003","Msmu-500_Mchi-450_ct-003","Msmu-500_Mchi-480_ct-003"
 ,"Msmu-200_Mchi-180_ct-010"
 ,"Msmu-250_Mchi-180_ct-010","Msmu-250_Mchi-200_ct-010"
 ,"Msmu-300_Mchi-180_ct-010","Msmu-300_Mchi-200_ct-010","Msmu-300_Mchi-250_ct-010","Msmu-300_Mchi-280_ct-010"
 ,"Msmu-350_Mchi-180_ct-010","Msmu-350_Mchi-200_ct-010","Msmu-350_Mchi-250_ct-010","Msmu-350_Mchi-300_ct-010","Msmu-350_Mchi-330_ct-010"
 ,"Msmu-400_Mchi-180_ct-010","Msmu-400_Mchi-200_ct-010","Msmu-400_Mchi-250_ct-010","Msmu-400_Mchi-300_ct-010","Msmu-400_Mchi-350_ct-010","Msmu-400_Mchi-380_ct-010"
 ,"Msmu-450_Mchi-180_ct-010","Msmu-450_Mchi-200_ct-010","Msmu-450_Mchi-250_ct-010","Msmu-450_Mchi-300_ct-010","Msmu-450_Mchi-350_ct-010","Msmu-450_Mchi-400_ct-010","Msmu-450_Mchi-430_ct-010"
 ,"Msmu-500_Mchi-180_ct-010","Msmu-500_Mchi-200_ct-010","Msmu-500_Mchi-250_ct-010","Msmu-500_Mchi-300_ct-010","Msmu-500_Mchi-350_ct-010","Msmu-500_Mchi-400_ct-010","Msmu-500_Mchi-450_ct-010","Msmu-500_Mchi-480_ct-010"
 ,"Msmu-200_Mchi-180_ct-030"
 ,"Msmu-250_Mchi-180_ct-030","Msmu-250_Mchi-200_ct-030"
 ,"Msmu-300_Mchi-180_ct-030","Msmu-300_Mchi-200_ct-030","Msmu-300_Mchi-250_ct-030","Msmu-300_Mchi-280_ct-030"
 ,"Msmu-350_Mchi-180_ct-030","Msmu-350_Mchi-200_ct-030","Msmu-350_Mchi-250_ct-030","Msmu-350_Mchi-300_ct-030","Msmu-350_Mchi-330_ct-030"
 ,"Msmu-400_Mchi-180_ct-030","Msmu-400_Mchi-200_ct-030","Msmu-400_Mchi-250_ct-030","Msmu-400_Mchi-300_ct-030","Msmu-400_Mchi-350_ct-030","Msmu-400_Mchi-380_ct-030"
 ,"Msmu-450_Mchi-180_ct-030","Msmu-450_Mchi-200_ct-030","Msmu-450_Mchi-250_ct-030","Msmu-450_Mchi-300_ct-030","Msmu-450_Mchi-350_ct-030","Msmu-450_Mchi-400_ct-030","Msmu-450_Mchi-430_ct-030"
 ,"Msmu-500_Mchi-180_ct-030","Msmu-500_Mchi-200_ct-030","Msmu-500_Mchi-250_ct-030","Msmu-500_Mchi-300_ct-030","Msmu-500_Mchi-350_ct-030","Msmu-500_Mchi-400_ct-030","Msmu-500_Mchi-450_ct-030","Msmu-500_Mchi-480_ct-030"
 ,"Msmu-200_Mchi-180_ct-100"
 ,"Msmu-250_Mchi-180_ct-100","Msmu-250_Mchi-200_ct-100"
 ,"Msmu-300_Mchi-180_ct-100","Msmu-300_Mchi-200_ct-100","Msmu-300_Mchi-250_ct-100","Msmu-300_Mchi-280_ct-100"
 ,"Msmu-350_Mchi-180_ct-100","Msmu-350_Mchi-200_ct-100","Msmu-350_Mchi-250_ct-100","Msmu-350_Mchi-300_ct-100","Msmu-350_Mchi-330_ct-100"
 ,"Msmu-400_Mchi-180_ct-100","Msmu-400_Mchi-200_ct-100","Msmu-400_Mchi-250_ct-100","Msmu-400_Mchi-300_ct-100","Msmu-400_Mchi-350_ct-100","Msmu-400_Mchi-380_ct-100"
 ,"Msmu-450_Mchi-180_ct-100","Msmu-450_Mchi-200_ct-100","Msmu-450_Mchi-250_ct-100","Msmu-450_Mchi-300_ct-100","Msmu-450_Mchi-350_ct-100","Msmu-450_Mchi-400_ct-100","Msmu-450_Mchi-430_ct-100"
 ,"Msmu-500_Mchi-180_ct-100","Msmu-500_Mchi-200_ct-100","Msmu-500_Mchi-250_ct-100","Msmu-500_Mchi-300_ct-100","Msmu-500_Mchi-350_ct-100","Msmu-500_Mchi-400_ct-100","Msmu-500_Mchi-450_ct-100","Msmu-500_Mchi-480_ct-100"
 ,"Msmu-200_Mchi-180_ct-300"
 ,"Msmu-250_Mchi-180_ct-300","Msmu-250_Mchi-200_ct-300"
 ,"Msmu-300_Mchi-180_ct-300","Msmu-300_Mchi-200_ct-300","Msmu-300_Mchi-250_ct-300","Msmu-300_Mchi-280_ct-300"
 ,"Msmu-350_Mchi-180_ct-300","Msmu-350_Mchi-200_ct-300","Msmu-350_Mchi-250_ct-300","Msmu-350_Mchi-300_ct-300","Msmu-350_Mchi-330_ct-300"
 ,"Msmu-400_Mchi-180_ct-300","Msmu-400_Mchi-200_ct-300","Msmu-400_Mchi-250_ct-300","Msmu-400_Mchi-300_ct-300","Msmu-400_Mchi-350_ct-300","Msmu-400_Mchi-380_ct-300"
 ,"Msmu-450_Mchi-180_ct-300","Msmu-450_Mchi-200_ct-300","Msmu-450_Mchi-250_ct-300","Msmu-450_Mchi-300_ct-300","Msmu-450_Mchi-350_ct-300","Msmu-450_Mchi-400_ct-300","Msmu-450_Mchi-430_ct-300"
 ,"Msmu-500_Mchi-180_ct-300","Msmu-500_Mchi-200_ct-300","Msmu-500_Mchi-250_ct-300","Msmu-500_Mchi-300_ct-300","Msmu-500_Mchi-350_ct-300","Msmu-500_Mchi-400_ct-300","Msmu-500_Mchi-450_ct-300","Msmu-500_Mchi-480_ct-300"
 ,"Msmu-200_Mchi-180_ct-1000"
 ,"Msmu-250_Mchi-180_ct-1000","Msmu-250_Mchi-200_ct-1000"
 ,"Msmu-300_Mchi-180_ct-1000","Msmu-300_Mchi-200_ct-1000","Msmu-300_Mchi-250_ct-1000","Msmu-300_Mchi-280_ct-1000"
 ,"Msmu-350_Mchi-180_ct-1000","Msmu-350_Mchi-200_ct-1000","Msmu-350_Mchi-250_ct-1000","Msmu-350_Mchi-300_ct-1000","Msmu-350_Mchi-330_ct-1000"
 ,"Msmu-400_Mchi-180_ct-1000","Msmu-400_Mchi-200_ct-1000","Msmu-400_Mchi-250_ct-1000","Msmu-400_Mchi-300_ct-1000","Msmu-400_Mchi-350_ct-1000","Msmu-400_Mchi-380_ct-1000"
 ,"Msmu-450_Mchi-180_ct-1000","Msmu-450_Mchi-200_ct-1000","Msmu-450_Mchi-250_ct-1000","Msmu-450_Mchi-300_ct-1000","Msmu-450_Mchi-350_ct-1000","Msmu-450_Mchi-400_ct-1000","Msmu-450_Mchi-430_ct-1000"
 ,"Msmu-500_Mchi-180_ct-1000","Msmu-500_Mchi-200_ct-1000","Msmu-500_Mchi-250_ct-1000","Msmu-500_Mchi-300_ct-1000","Msmu-500_Mchi-350_ct-1000","Msmu-500_Mchi-400_ct-1000","Msmu-500_Mchi-450_ct-1000","Msmu-500_Mchi-480_ct-1000"
]

year = "2022"
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
    sedFile = "sed -e 's#inputSample#" + list_arguments[i]+"#' crab_config_RPV_2022B_default.py > crab_config_temp.py"
    os.system(sedFile) 
    command = "crab-pre submit -c crab_config_temp.py"
    os.system(command) 

