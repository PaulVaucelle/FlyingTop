#include "TChain.h"
#include "MiniNtuple.h"
// C++ includes
#include <iostream>
#include <fstream>
#include "TROOT.h"

int main(int argc, char **argv)
{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniNtuple.C+g") ;
  
 if (gROOT->GetClass("MiniNtuple")==0) return 0;
 
 TChain c("FlyingTop/ttree");

TString Prod = "RPV_2022B_JERUp";
bool Signal = true;
TString GlobalPath = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"; 
////////////////////////////////////////////////////////////////////////////////

       TString SignalSet001[34]={
         "RPV_2022B_Msmu-200_Mchi-180_ct-001","RPV_2022B_Msmu-250_Mchi-180_ct-001","RPV_2022B_Msmu-250_Mchi-200_ct-001",
"RPV_2022B_Msmu-250_Mchi-230_ct-001","RPV_2022B_Msmu-300_Mchi-180_ct-001","RPV_2022B_Msmu-300_Mchi-200_ct-001","RPV_2022B_Msmu-300_Mchi-250_ct-001","RPV_2022B_Msmu-300_Mchi-280_ct-001",
"RPV_2022B_Msmu-350_Mchi-180_ct-001","RPV_2022B_Msmu-350_Mchi-200_ct-001","RPV_2022B_Msmu-350_Mchi-250_ct-001","RPV_2022B_Msmu-350_Mchi-300_ct-001","RPV_2022B_Msmu-350_Mchi-330_ct-001",
"RPV_2022B_Msmu-400_Mchi-180_ct-001","RPV_2022B_Msmu-400_Mchi-200_ct-001","RPV_2022B_Msmu-400_Mchi-250_ct-001","RPV_2022B_Msmu-400_Mchi-300_ct-001","RPV_2022B_Msmu-400_Mchi-350_ct-001",
"RPV_2022B_Msmu-400_Mchi-380_ct-001","RPV_2022B_Msmu-450_Mchi-180_ct-001","RPV_2022B_Msmu-450_Mchi-200_ct-001","RPV_2022B_Msmu-450_Mchi-250_ct-001","RPV_2022B_Msmu-450_Mchi-300_ct-001",
"RPV_2022B_Msmu-450_Mchi-350_ct-001","RPV_2022B_Msmu-450_Mchi-400_ct-001","RPV_2022B_Msmu-450_Mchi-430_ct-001","RPV_2022B_Msmu-500_Mchi-180_ct-001",
"RPV_2022B_Msmu-500_Mchi-200_ct-001"
,
"RPV_2022B_Msmu-500_Mchi-250_ct-001","RPV_2022B_Msmu-500_Mchi-300_ct-001","RPV_2022B_Msmu-500_Mchi-350_ct-001","RPV_2022B_Msmu-500_Mchi-400_ct-001","RPV_2022B_Msmu-500_Mchi-450_ct-001",
"RPV_2022B_Msmu-500_Mchi-480_ct-001"
}; 
 for (int i = 0; i< 34 ; i++) 
    {
      TString Path = GlobalPath+"/"+Prod+"/"+SignalSet001[i]+".root";
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet001[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(SignalSet001[i],Prod,Signal);
    }


      TString SignalSet003[34]={
         "RPV_2022B_Msmu-200_Mchi-180_ct-003","RPV_2022B_Msmu-250_Mchi-180_ct-003","RPV_2022B_Msmu-250_Mchi-200_ct-003",
"RPV_2022B_Msmu-250_Mchi-230_ct-003","RPV_2022B_Msmu-300_Mchi-180_ct-003","RPV_2022B_Msmu-300_Mchi-200_ct-003","RPV_2022B_Msmu-300_Mchi-250_ct-003","RPV_2022B_Msmu-300_Mchi-280_ct-003",
"RPV_2022B_Msmu-350_Mchi-180_ct-003","RPV_2022B_Msmu-350_Mchi-200_ct-003","RPV_2022B_Msmu-350_Mchi-250_ct-003","RPV_2022B_Msmu-350_Mchi-300_ct-003","RPV_2022B_Msmu-350_Mchi-330_ct-003",
"RPV_2022B_Msmu-400_Mchi-180_ct-003","RPV_2022B_Msmu-400_Mchi-200_ct-003","RPV_2022B_Msmu-400_Mchi-250_ct-003","RPV_2022B_Msmu-400_Mchi-300_ct-003","RPV_2022B_Msmu-400_Mchi-350_ct-003",
"RPV_2022B_Msmu-400_Mchi-380_ct-003","RPV_2022B_Msmu-450_Mchi-180_ct-003","RPV_2022B_Msmu-450_Mchi-200_ct-003","RPV_2022B_Msmu-450_Mchi-250_ct-003","RPV_2022B_Msmu-450_Mchi-300_ct-003",
"RPV_2022B_Msmu-450_Mchi-350_ct-003","RPV_2022B_Msmu-450_Mchi-400_ct-003","RPV_2022B_Msmu-450_Mchi-430_ct-003","RPV_2022B_Msmu-500_Mchi-180_ct-003",
"RPV_2022B_Msmu-500_Mchi-200_ct-003"
,
"RPV_2022B_Msmu-500_Mchi-250_ct-003","RPV_2022B_Msmu-500_Mchi-300_ct-003","RPV_2022B_Msmu-500_Mchi-350_ct-003","RPV_2022B_Msmu-500_Mchi-400_ct-003","RPV_2022B_Msmu-500_Mchi-450_ct-003",
"RPV_2022B_Msmu-500_Mchi-480_ct-003"
};
 for (int i = 3; i< 4 ; i++) 
    {
      
      TString Path = GlobalPath+"/"+Prod+"/"+SignalSet003[i]+".root";
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet003[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(SignalSet003[i],Prod,Signal);
    }

      TString SignalSet010[34]={
         "RPV_2022B_Msmu-200_Mchi-180_ct-010","RPV_2022B_Msmu-250_Mchi-180_ct-010","RPV_2022B_Msmu-250_Mchi-200_ct-010",
"RPV_2022B_Msmu-250_Mchi-230_ct-010","RPV_2022B_Msmu-300_Mchi-180_ct-010","RPV_2022B_Msmu-300_Mchi-200_ct-010","RPV_2022B_Msmu-300_Mchi-250_ct-010","RPV_2022B_Msmu-300_Mchi-280_ct-010",
"RPV_2022B_Msmu-350_Mchi-180_ct-010","RPV_2022B_Msmu-350_Mchi-200_ct-010","RPV_2022B_Msmu-350_Mchi-250_ct-010","RPV_2022B_Msmu-350_Mchi-300_ct-010","RPV_2022B_Msmu-350_Mchi-330_ct-010",
"RPV_2022B_Msmu-400_Mchi-180_ct-010","RPV_2022B_Msmu-400_Mchi-200_ct-010","RPV_2022B_Msmu-400_Mchi-250_ct-010","RPV_2022B_Msmu-400_Mchi-300_ct-010","RPV_2022B_Msmu-400_Mchi-350_ct-010",
"RPV_2022B_Msmu-400_Mchi-380_ct-010","RPV_2022B_Msmu-450_Mchi-180_ct-010","RPV_2022B_Msmu-450_Mchi-200_ct-010","RPV_2022B_Msmu-450_Mchi-250_ct-010","RPV_2022B_Msmu-450_Mchi-300_ct-010",
"RPV_2022B_Msmu-450_Mchi-350_ct-010","RPV_2022B_Msmu-450_Mchi-400_ct-010","RPV_2022B_Msmu-450_Mchi-430_ct-010","RPV_2022B_Msmu-500_Mchi-180_ct-010",
"RPV_2022B_Msmu-500_Mchi-200_ct-010"
,
"RPV_2022B_Msmu-500_Mchi-250_ct-010","RPV_2022B_Msmu-500_Mchi-300_ct-010","RPV_2022B_Msmu-500_Mchi-350_ct-010","RPV_2022B_Msmu-500_Mchi-400_ct-010","RPV_2022B_Msmu-500_Mchi-450_ct-010",
"RPV_2022B_Msmu-500_Mchi-480_ct-010"
};
 for (int i = 0; i< 34 ; i++) 
    {
      TString Path = GlobalPath+"/"+Prod+"/"+SignalSet010[i]+".root";
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet010[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(SignalSet010[i],Prod,Signal);
    }
      TString SignalSet030[34]={
         "RPV_2022B_Msmu-200_Mchi-180_ct-030","RPV_2022B_Msmu-250_Mchi-180_ct-030","RPV_2022B_Msmu-250_Mchi-200_ct-030",
"RPV_2022B_Msmu-250_Mchi-230_ct-030","RPV_2022B_Msmu-300_Mchi-180_ct-030","RPV_2022B_Msmu-300_Mchi-200_ct-030","RPV_2022B_Msmu-300_Mchi-250_ct-030","RPV_2022B_Msmu-300_Mchi-280_ct-030",
"RPV_2022B_Msmu-350_Mchi-180_ct-030","RPV_2022B_Msmu-350_Mchi-200_ct-030","RPV_2022B_Msmu-350_Mchi-250_ct-030","RPV_2022B_Msmu-350_Mchi-300_ct-030","RPV_2022B_Msmu-350_Mchi-330_ct-030",
"RPV_2022B_Msmu-400_Mchi-180_ct-030","RPV_2022B_Msmu-400_Mchi-200_ct-030","RPV_2022B_Msmu-400_Mchi-250_ct-030","RPV_2022B_Msmu-400_Mchi-300_ct-030","RPV_2022B_Msmu-400_Mchi-350_ct-030",
"RPV_2022B_Msmu-400_Mchi-380_ct-030","RPV_2022B_Msmu-450_Mchi-180_ct-030","RPV_2022B_Msmu-450_Mchi-200_ct-030","RPV_2022B_Msmu-450_Mchi-250_ct-030","RPV_2022B_Msmu-450_Mchi-300_ct-030",
"RPV_2022B_Msmu-450_Mchi-350_ct-030","RPV_2022B_Msmu-450_Mchi-400_ct-030","RPV_2022B_Msmu-450_Mchi-430_ct-030","RPV_2022B_Msmu-500_Mchi-180_ct-030",
"RPV_2022B_Msmu-500_Mchi-200_ct-030"
,
"RPV_2022B_Msmu-500_Mchi-250_ct-030","RPV_2022B_Msmu-500_Mchi-300_ct-030","RPV_2022B_Msmu-500_Mchi-350_ct-030","RPV_2022B_Msmu-500_Mchi-400_ct-030","RPV_2022B_Msmu-500_Mchi-450_ct-030",
"RPV_2022B_Msmu-500_Mchi-480_ct-030"
};
 for (int i = 0 ; i< 34 ; i++) 
    {
      TString Path = GlobalPath+"/"+Prod+"/"+SignalSet030[i]+".root";
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet030[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(SignalSet030[i],Prod,Signal);
    }

    TString SignalSet100[34]={
      "RPV_2022B_Msmu-200_Mchi-180_ct-100","RPV_2022B_Msmu-250_Mchi-180_ct-100","RPV_2022B_Msmu-250_Mchi-200_ct-100",
"RPV_2022B_Msmu-250_Mchi-230_ct-100","RPV_2022B_Msmu-300_Mchi-180_ct-100","RPV_2022B_Msmu-300_Mchi-200_ct-100","RPV_2022B_Msmu-300_Mchi-250_ct-100","RPV_2022B_Msmu-300_Mchi-280_ct-100",
"RPV_2022B_Msmu-350_Mchi-180_ct-100","RPV_2022B_Msmu-350_Mchi-200_ct-100","RPV_2022B_Msmu-350_Mchi-250_ct-100","RPV_2022B_Msmu-350_Mchi-300_ct-100","RPV_2022B_Msmu-350_Mchi-330_ct-100",
"RPV_2022B_Msmu-400_Mchi-180_ct-100","RPV_2022B_Msmu-400_Mchi-200_ct-100","RPV_2022B_Msmu-400_Mchi-250_ct-100","RPV_2022B_Msmu-400_Mchi-300_ct-100","RPV_2022B_Msmu-400_Mchi-350_ct-100",
"RPV_2022B_Msmu-400_Mchi-380_ct-100","RPV_2022B_Msmu-450_Mchi-180_ct-100","RPV_2022B_Msmu-450_Mchi-200_ct-100","RPV_2022B_Msmu-450_Mchi-250_ct-100","RPV_2022B_Msmu-450_Mchi-300_ct-100",
"RPV_2022B_Msmu-450_Mchi-350_ct-100","RPV_2022B_Msmu-450_Mchi-400_ct-100","RPV_2022B_Msmu-450_Mchi-430_ct-100","RPV_2022B_Msmu-500_Mchi-180_ct-100",
"RPV_2022B_Msmu-500_Mchi-200_ct-100"
,
"RPV_2022B_Msmu-500_Mchi-250_ct-100","RPV_2022B_Msmu-500_Mchi-300_ct-100","RPV_2022B_Msmu-500_Mchi-350_ct-100","RPV_2022B_Msmu-500_Mchi-400_ct-100","RPV_2022B_Msmu-500_Mchi-450_ct-100",
"RPV_2022B_Msmu-500_Mchi-480_ct-100"
};

 for (int i = 0; i< 34 ; i++) 
    {
      TString Path = GlobalPath+"/"+Prod+"/"+SignalSet100[i]+".root";
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet100[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(SignalSet100[i],Prod,Signal);
    }


      TString SignalSet300[34]={
         "RPV_2022B_Msmu-200_Mchi-180_ct-300","RPV_2022B_Msmu-250_Mchi-180_ct-300","RPV_2022B_Msmu-250_Mchi-200_ct-300",
"RPV_2022B_Msmu-250_Mchi-230_ct-300","RPV_2022B_Msmu-300_Mchi-180_ct-300","RPV_2022B_Msmu-300_Mchi-200_ct-300","RPV_2022B_Msmu-300_Mchi-250_ct-300","RPV_2022B_Msmu-300_Mchi-280_ct-300",
"RPV_2022B_Msmu-350_Mchi-180_ct-300","RPV_2022B_Msmu-350_Mchi-200_ct-300","RPV_2022B_Msmu-350_Mchi-250_ct-300","RPV_2022B_Msmu-350_Mchi-300_ct-300","RPV_2022B_Msmu-350_Mchi-330_ct-300",
"RPV_2022B_Msmu-400_Mchi-180_ct-300","RPV_2022B_Msmu-400_Mchi-200_ct-300","RPV_2022B_Msmu-400_Mchi-250_ct-300","RPV_2022B_Msmu-400_Mchi-300_ct-300","RPV_2022B_Msmu-400_Mchi-350_ct-300",
"RPV_2022B_Msmu-400_Mchi-380_ct-300","RPV_2022B_Msmu-450_Mchi-180_ct-300","RPV_2022B_Msmu-450_Mchi-200_ct-300","RPV_2022B_Msmu-450_Mchi-250_ct-300","RPV_2022B_Msmu-450_Mchi-300_ct-300",
"RPV_2022B_Msmu-450_Mchi-350_ct-300","RPV_2022B_Msmu-450_Mchi-400_ct-300","RPV_2022B_Msmu-450_Mchi-430_ct-300","RPV_2022B_Msmu-500_Mchi-180_ct-300",
"RPV_2022B_Msmu-500_Mchi-200_ct-300"
,
"RPV_2022B_Msmu-500_Mchi-250_ct-300","RPV_2022B_Msmu-500_Mchi-300_ct-300","RPV_2022B_Msmu-500_Mchi-350_ct-300","RPV_2022B_Msmu-500_Mchi-400_ct-300","RPV_2022B_Msmu-500_Mchi-450_ct-300",
"RPV_2022B_Msmu-500_Mchi-480_ct-300"
};
 for (int i = 0 ; i< 34 ; i++) 
    {
      TString Path = GlobalPath+"/"+Prod+"/"+SignalSet300[i]+".root";
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet300[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(SignalSet300[i],Prod,Signal);
    }

          TString SignalSet1000[34]={
            "RPV_2022B_Msmu-200_Mchi-180_ct-1000","RPV_2022B_Msmu-250_Mchi-180_ct-1000","RPV_2022B_Msmu-250_Mchi-200_ct-1000",
"RPV_2022B_Msmu-250_Mchi-230_ct-1000","RPV_2022B_Msmu-300_Mchi-180_ct-1000","RPV_2022B_Msmu-300_Mchi-200_ct-1000","RPV_2022B_Msmu-300_Mchi-250_ct-1000","RPV_2022B_Msmu-300_Mchi-280_ct-1000",
"RPV_2022B_Msmu-350_Mchi-180_ct-1000","RPV_2022B_Msmu-350_Mchi-200_ct-1000","RPV_2022B_Msmu-350_Mchi-250_ct-1000","RPV_2022B_Msmu-350_Mchi-300_ct-1000","RPV_2022B_Msmu-350_Mchi-330_ct-1000",
"RPV_2022B_Msmu-400_Mchi-180_ct-1000","RPV_2022B_Msmu-400_Mchi-200_ct-1000","RPV_2022B_Msmu-400_Mchi-250_ct-1000","RPV_2022B_Msmu-400_Mchi-300_ct-1000","RPV_2022B_Msmu-400_Mchi-350_ct-1000",
"RPV_2022B_Msmu-400_Mchi-380_ct-1000","RPV_2022B_Msmu-450_Mchi-180_ct-1000","RPV_2022B_Msmu-450_Mchi-200_ct-1000","RPV_2022B_Msmu-450_Mchi-250_ct-1000","RPV_2022B_Msmu-450_Mchi-300_ct-1000",
"RPV_2022B_Msmu-450_Mchi-350_ct-1000","RPV_2022B_Msmu-450_Mchi-400_ct-1000","RPV_2022B_Msmu-450_Mchi-430_ct-1000","RPV_2022B_Msmu-500_Mchi-180_ct-1000",
"RPV_2022B_Msmu-500_Mchi-200_ct-1000"
,
"RPV_2022B_Msmu-500_Mchi-250_ct-1000","RPV_2022B_Msmu-500_Mchi-300_ct-1000","RPV_2022B_Msmu-500_Mchi-350_ct-1000","RPV_2022B_Msmu-500_Mchi-400_ct-1000","RPV_2022B_Msmu-500_Mchi-450_ct-1000",
"RPV_2022B_Msmu-500_Mchi-480_ct-1000"
};
 for (int i = 0 ; i< 34 ; i++) 
    {
      TString Path = GlobalPath+"/"+Prod+"/"+SignalSet1000[i]+".root";
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet1000[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(SignalSet1000[i],Prod,Signal);
    }
return 0;
}
