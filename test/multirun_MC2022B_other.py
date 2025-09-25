import os

list_arguments = [
  "RPV_2022B_smu350_neu180_ctau001","RPV_2022B_smu350_neu200_ctau001","RPV_2022B_smu350_neu250_ctau001","RPV_2022B_smu350_neu300_ctau001","RPV_2022B_smu350_neu330_ctau001"
 ,"RPV_2022B_smu450_neu180_ctau001","RPV_2022B_smu450_neu200_ctau001","RPV_2022B_smu450_neu250_ctau001","RPV_2022B_smu450_neu300_ctau001","RPV_2022B_smu450_neu350_ctau001","RPV_2022B_smu450_neu400_ctau001","RPV_2022B_smu450_neu430_ctau001"
 ,"RPV_2022B_smu350_neu180_ctau010","RPV_2022B_smu350_neu200_ctau010","RPV_2022B_smu350_neu250_ctau010","RPV_2022B_smu350_neu300_ctau010","RPV_2022B_smu350_neu330_ctau010"
 ,"RPV_2022B_smu450_neu180_ctau010","RPV_2022B_smu450_neu200_ctau010","RPV_2022B_smu450_neu250_ctau010","RPV_2022B_smu450_neu300_ctau010","RPV_2022B_smu450_neu350_ctau010","RPV_2022B_smu450_neu400_ctau010","RPV_2022B_smu450_neu430_ctau010"
 ,"RPV_2022B_smu350_neu180_ctau100","RPV_2022B_smu350_neu200_ctau100","RPV_2022B_smu350_neu250_ctau100","RPV_2022B_smu350_neu300_ctau100","RPV_2022B_smu350_neu330_ctau100"
 ,"RPV_2022B_smu450_neu180_ctau100","RPV_2022B_smu450_neu200_ctau100","RPV_2022B_smu450_neu250_ctau100","RPV_2022B_smu450_neu300_ctau100","RPV_2022B_smu450_neu350_ctau100","RPV_2022B_smu450_neu400_ctau100","RPV_2022B_smu450_neu430_ctau100"
 ,"RPV_2022B_smu350_neu180_ctau300","RPV_2022B_smu350_neu200_ctau300","RPV_2022B_smu350_neu250_ctau300","RPV_2022B_smu350_neu300_ctau300","RPV_2022B_smu350_neu330_ctau300"
 ,"RPV_2022B_smu450_neu180_ctau300","RPV_2022B_smu450_neu200_ctau300","RPV_2022B_smu450_neu250_ctau300","RPV_2022B_smu450_neu300_ctau300","RPV_2022B_smu450_neu350_ctau300","RPV_2022B_smu450_neu400_ctau300","RPV_2022B_smu450_neu430_ctau300"
 ,"RPV_2022B_smu200_neu180_ctau003"
 ,"RPV_2022B_smu250_neu180_ctau003","RPV_2022B_smu250_neu200_ctau003","RPV_2022B_smu250_neu230_ctau003"
 ,"RPV_2022B_smu300_neu180_ctau003","RPV_2022B_smu300_neu200_ctau003","RPV_2022B_smu300_neu250_ctau003","RPV_2022B_smu300_neu280_ctau003"
 ,"RPV_2022B_smu350_neu180_ctau003","RPV_2022B_smu350_neu200_ctau003","RPV_2022B_smu350_neu250_ctau003","RPV_2022B_smu350_neu300_ctau003","RPV_2022B_smu350_neu330_ctau003"
 ,"RPV_2022B_smu400_neu180_ctau003","RPV_2022B_smu400_neu200_ctau003","RPV_2022B_smu400_neu250_ctau003","RPV_2022B_smu400_neu300_ctau003","RPV_2022B_smu400_neu350_ctau003","RPV_2022B_smu400_neu380_ctau003"
 ,"RPV_2022B_smu450_neu180_ctau003","RPV_2022B_smu450_neu200_ctau003","RPV_2022B_smu450_neu250_ctau003","RPV_2022B_smu450_neu300_ctau003","RPV_2022B_smu450_neu350_ctau003","RPV_2022B_smu450_neu400_ctau003","RPV_2022B_smu450_neu430_ctau003"
 ,"RPV_2022B_smu500_neu180_ctau003","RPV_2022B_smu500_neu200_ctau003","RPV_2022B_smu500_neu250_ctau003","RPV_2022B_smu500_neu300_ctau003","RPV_2022B_smu500_neu350_ctau003","RPV_2022B_smu500_neu400_ctau003","RPV_2022B_smu500_neu450_ctau003","RPV_2022B_smu500_neu480_ctau003"
 ,"RPV_2022B_smu200_neu180_ctau030"
 ,"RPV_2022B_smu250_neu180_ctau030","RPV_2022B_smu250_neu200_ctau030","RPV_2022B_smu250_neu230_ctau030"
 ,"RPV_2022B_smu300_neu180_ctau030","RPV_2022B_smu300_neu200_ctau030","RPV_2022B_smu300_neu250_ctau030","RPV_2022B_smu300_neu280_ctau030"
 ,"RPV_2022B_smu350_neu180_ctau030","RPV_2022B_smu350_neu200_ctau030","RPV_2022B_smu350_neu250_ctau030","RPV_2022B_smu350_neu300_ctau030","RPV_2022B_smu350_neu330_ctau030"
 ,"RPV_2022B_smu400_neu180_ctau030","RPV_2022B_smu400_neu200_ctau030","RPV_2022B_smu400_neu250_ctau030","RPV_2022B_smu400_neu300_ctau030","RPV_2022B_smu400_neu350_ctau030","RPV_2022B_smu400_neu380_ctau030"
 ,"RPV_2022B_smu450_neu180_ctau030","RPV_2022B_smu450_neu200_ctau030","RPV_2022B_smu450_neu250_ctau030","RPV_2022B_smu450_neu300_ctau030","RPV_2022B_smu450_neu350_ctau030","RPV_2022B_smu450_neu400_ctau030","RPV_2022B_smu450_neu430_ctau030"
 ,"RPV_2022B_smu500_neu180_ctau030","RPV_2022B_smu500_neu200_ctau030","RPV_2022B_smu500_neu250_ctau030","RPV_2022B_smu500_neu300_ctau030","RPV_2022B_smu500_neu350_ctau030","RPV_2022B_smu500_neu400_ctau030","RPV_2022B_smu500_neu450_ctau030","RPV_2022B_smu500_neu480_ctau030"
 ,"RPV_2022B_smu200_neu180_ctau1000"
 ,"RPV_2022B_smu250_neu180_ctau1000","RPV_2022B_smu250_neu200_ctau1000","RPV_2022B_smu250_neu230_ctau1000"
 ,"RPV_2022B_smu300_neu180_ctau1000","RPV_2022B_smu300_neu200_ctau1000","RPV_2022B_smu300_neu250_ctau1000","RPV_2022B_smu300_neu280_ctau1000"
 ,"RPV_2022B_smu350_neu180_ctau1000","RPV_2022B_smu350_neu200_ctau1000","RPV_2022B_smu350_neu250_ctau1000","RPV_2022B_smu350_neu300_ctau1000","RPV_2022B_smu350_neu330_ctau1000"
 ,"RPV_2022B_smu400_neu180_ctau1000","RPV_2022B_smu400_neu200_ctau1000","RPV_2022B_smu400_neu250_ctau1000","RPV_2022B_smu400_neu300_ctau1000","RPV_2022B_smu400_neu350_ctau1000","RPV_2022B_smu400_neu380_ctau1000"
 ,"RPV_2022B_smu450_neu180_ctau1000","RPV_2022B_smu450_neu200_ctau1000","RPV_2022B_smu450_neu250_ctau1000","RPV_2022B_smu450_neu300_ctau1000","RPV_2022B_smu450_neu350_ctau1000","RPV_2022B_smu450_neu400_ctau1000","RPV_2022B_smu450_neu430_ctau1000"
 ,"RPV_2022B_smu500_neu180_ctau1000","RPV_2022B_smu500_neu200_ctau1000","RPV_2022B_smu500_neu250_ctau1000","RPV_2022B_smu500_neu300_ctau1000","RPV_2022B_smu500_neu350_ctau1000","RPV_2022B_smu500_neu400_ctau1000","RPV_2022B_smu500_neu450_ctau1000","RPV_2022B_smu500_neu480_ctau1000"
]

for i in range(len(list_arguments)):
    command = "rm flyingtop_temp.py"
    os.system(command) 
    sedFile = "sed -e 's#inputFile#" + list_arguments[i]+"#' flyingtop_default_MC2022B.py > flyingtop_temp.py"
    os.system(sedFile) 
    command = "cmsRun flyingtop_temp.py"
    os.system(command) 

