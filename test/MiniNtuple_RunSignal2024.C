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

TString Prod = "RPV_2024_JECUp";
bool Signal = true;
TString GlobalPath = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"; 
////////////////////////////////////////////////////////////////////////////////
 
      TString SignalSet001[34]={
         "RPV_2024_Par-ct-001-MChi-180-MSmu-200","RPV_2024_Par-ct-001-MChi-180-MSmu-250","RPV_2024_Par-ct-001-MChi-200-MSmu-250","RPV_2024_Par-ct-001-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-300","RPV_2024_Par-ct-001-MChi-200-MSmu-300","RPV_2024_Par-ct-001-MChi-250-MSmu-300","RPV_2024_Par-ct-001-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-350","RPV_2024_Par-ct-001-MChi-200-MSmu-350","RPV_2024_Par-ct-001-MChi-250-MSmu-350","RPV_2024_Par-ct-001-MChi-300-MSmu-350","RPV_2024_Par-ct-001-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-400","RPV_2024_Par-ct-001-MChi-200-MSmu-400","RPV_2024_Par-ct-001-MChi-250-MSmu-400","RPV_2024_Par-ct-001-MChi-300-MSmu-400","RPV_2024_Par-ct-001-MChi-350-MSmu-400","RPV_2024_Par-ct-001-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-450","RPV_2024_Par-ct-001-MChi-200-MSmu-450","RPV_2024_Par-ct-001-MChi-250-MSmu-450","RPV_2024_Par-ct-001-MChi-300-MSmu-450","RPV_2024_Par-ct-001-MChi-350-MSmu-450","RPV_2024_Par-ct-001-MChi-400-MSmu-450","RPV_2024_Par-ct-001-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-500","RPV_2024_Par-ct-001-MChi-200-MSmu-500","RPV_2024_Par-ct-001-MChi-250-MSmu-500","RPV_2024_Par-ct-001-MChi-300-MSmu-500","RPV_2024_Par-ct-001-MChi-350-MSmu-500","RPV_2024_Par-ct-001-MChi-400-MSmu-500","RPV_2024_Par-ct-001-MChi-450-MSmu-500","RPV_2024_Par-ct-001-MChi-480-MSmu-500"
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
         "RPV_2024_Par-ct-003-MChi-180-MSmu-200","RPV_2024_Par-ct-003-MChi-180-MSmu-250","RPV_2024_Par-ct-003-MChi-200-MSmu-250","RPV_2024_Par-ct-003-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-300","RPV_2024_Par-ct-003-MChi-200-MSmu-300","RPV_2024_Par-ct-003-MChi-250-MSmu-300","RPV_2024_Par-ct-003-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-350","RPV_2024_Par-ct-003-MChi-200-MSmu-350","RPV_2024_Par-ct-003-MChi-250-MSmu-350","RPV_2024_Par-ct-003-MChi-300-MSmu-350","RPV_2024_Par-ct-003-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-400","RPV_2024_Par-ct-003-MChi-200-MSmu-400","RPV_2024_Par-ct-003-MChi-250-MSmu-400","RPV_2024_Par-ct-003-MChi-300-MSmu-400","RPV_2024_Par-ct-003-MChi-350-MSmu-400","RPV_2024_Par-ct-003-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-450","RPV_2024_Par-ct-003-MChi-200-MSmu-450","RPV_2024_Par-ct-003-MChi-250-MSmu-450","RPV_2024_Par-ct-003-MChi-300-MSmu-450","RPV_2024_Par-ct-003-MChi-350-MSmu-450","RPV_2024_Par-ct-003-MChi-400-MSmu-450","RPV_2024_Par-ct-003-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-500","RPV_2024_Par-ct-003-MChi-200-MSmu-500","RPV_2024_Par-ct-003-MChi-250-MSmu-500","RPV_2024_Par-ct-003-MChi-300-MSmu-500","RPV_2024_Par-ct-003-MChi-350-MSmu-500","RPV_2024_Par-ct-003-MChi-400-MSmu-500","RPV_2024_Par-ct-003-MChi-450-MSmu-500","RPV_2024_Par-ct-003-MChi-480-MSmu-500"

};
 for (int i = 0; i< 34 ; i++) 
    {
      
      TString Path = GlobalPath+"/"+Prod+"/"+SignalSet003[i]+".root";
      // TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+SignalSet003[i]+".root";
      c.Reset();
      c.Add(Path);
       MiniNtuple* t = new MiniNtuple(&c);
       t->Loop(SignalSet003[i],Prod,Signal);
    }

      TString SignalSet010[34]={
         "RPV_2024_Par-ct-010-MChi-180-MSmu-200","RPV_2024_Par-ct-010-MChi-180-MSmu-250","RPV_2024_Par-ct-010-MChi-200-MSmu-250","RPV_2024_Par-ct-010-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-300","RPV_2024_Par-ct-010-MChi-200-MSmu-300","RPV_2024_Par-ct-010-MChi-250-MSmu-300","RPV_2024_Par-ct-010-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-350","RPV_2024_Par-ct-010-MChi-200-MSmu-350","RPV_2024_Par-ct-010-MChi-250-MSmu-350","RPV_2024_Par-ct-010-MChi-300-MSmu-350","RPV_2024_Par-ct-010-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-400","RPV_2024_Par-ct-010-MChi-200-MSmu-400","RPV_2024_Par-ct-010-MChi-250-MSmu-400","RPV_2024_Par-ct-010-MChi-300-MSmu-400","RPV_2024_Par-ct-010-MChi-350-MSmu-400","RPV_2024_Par-ct-010-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-450","RPV_2024_Par-ct-010-MChi-200-MSmu-450","RPV_2024_Par-ct-010-MChi-250-MSmu-450","RPV_2024_Par-ct-010-MChi-300-MSmu-450","RPV_2024_Par-ct-010-MChi-350-MSmu-450","RPV_2024_Par-ct-010-MChi-400-MSmu-450","RPV_2024_Par-ct-010-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-500","RPV_2024_Par-ct-010-MChi-200-MSmu-500","RPV_2024_Par-ct-010-MChi-250-MSmu-500","RPV_2024_Par-ct-010-MChi-300-MSmu-500","RPV_2024_Par-ct-010-MChi-350-MSmu-500","RPV_2024_Par-ct-010-MChi-400-MSmu-500","RPV_2024_Par-ct-010-MChi-450-MSmu-500","RPV_2024_Par-ct-010-MChi-480-MSmu-500"

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
         "RPV_2024_Par-ct-030-MChi-180-MSmu-200","RPV_2024_Par-ct-030-MChi-180-MSmu-250","RPV_2024_Par-ct-030-MChi-200-MSmu-250","RPV_2024_Par-ct-030-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-030-MChi-180-MSmu-300","RPV_2024_Par-ct-030-MChi-200-MSmu-300","RPV_2024_Par-ct-030-MChi-250-MSmu-300","RPV_2024_Par-ct-030-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-030-MChi-180-MSmu-350","RPV_2024_Par-ct-030-MChi-200-MSmu-350","RPV_2024_Par-ct-030-MChi-250-MSmu-350","RPV_2024_Par-ct-030-MChi-300-MSmu-350","RPV_2024_Par-ct-030-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-030-MChi-180-MSmu-400","RPV_2024_Par-ct-030-MChi-200-MSmu-400","RPV_2024_Par-ct-030-MChi-250-MSmu-400","RPV_2024_Par-ct-030-MChi-300-MSmu-400","RPV_2024_Par-ct-030-MChi-350-MSmu-400","RPV_2024_Par-ct-030-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-030-MChi-180-MSmu-450","RPV_2024_Par-ct-030-MChi-200-MSmu-450","RPV_2024_Par-ct-030-MChi-250-MSmu-450","RPV_2024_Par-ct-030-MChi-300-MSmu-450","RPV_2024_Par-ct-030-MChi-350-MSmu-450","RPV_2024_Par-ct-030-MChi-400-MSmu-450","RPV_2024_Par-ct-030-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-030-MChi-180-MSmu-500","RPV_2024_Par-ct-030-MChi-200-MSmu-500","RPV_2024_Par-ct-030-MChi-250-MSmu-500","RPV_2024_Par-ct-030-MChi-300-MSmu-500","RPV_2024_Par-ct-030-MChi-350-MSmu-500","RPV_2024_Par-ct-030-MChi-400-MSmu-500","RPV_2024_Par-ct-030-MChi-450-MSmu-500","RPV_2024_Par-ct-030-MChi-480-MSmu-500"

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
         "RPV_2024_Par-ct-100-MChi-180-MSmu-200","RPV_2024_Par-ct-100-MChi-180-MSmu-250","RPV_2024_Par-ct-100-MChi-200-MSmu-250","RPV_2024_Par-ct-100-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-300","RPV_2024_Par-ct-100-MChi-200-MSmu-300","RPV_2024_Par-ct-100-MChi-250-MSmu-300","RPV_2024_Par-ct-100-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-350","RPV_2024_Par-ct-100-MChi-200-MSmu-350","RPV_2024_Par-ct-100-MChi-250-MSmu-350","RPV_2024_Par-ct-100-MChi-300-MSmu-350","RPV_2024_Par-ct-100-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-400","RPV_2024_Par-ct-100-MChi-200-MSmu-400","RPV_2024_Par-ct-100-MChi-250-MSmu-400","RPV_2024_Par-ct-100-MChi-300-MSmu-400","RPV_2024_Par-ct-100-MChi-350-MSmu-400","RPV_2024_Par-ct-100-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-450","RPV_2024_Par-ct-100-MChi-200-MSmu-450","RPV_2024_Par-ct-100-MChi-250-MSmu-450","RPV_2024_Par-ct-100-MChi-300-MSmu-450","RPV_2024_Par-ct-100-MChi-350-MSmu-450","RPV_2024_Par-ct-100-MChi-400-MSmu-450","RPV_2024_Par-ct-100-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-500","RPV_2024_Par-ct-100-MChi-200-MSmu-500","RPV_2024_Par-ct-100-MChi-250-MSmu-500","RPV_2024_Par-ct-100-MChi-300-MSmu-500","RPV_2024_Par-ct-100-MChi-350-MSmu-500","RPV_2024_Par-ct-100-MChi-400-MSmu-500","RPV_2024_Par-ct-100-MChi-450-MSmu-500","RPV_2024_Par-ct-100-MChi-480-MSmu-500"

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
         "RPV_2024_Par-ct-300-MChi-180-MSmu-200","RPV_2024_Par-ct-300-MChi-180-MSmu-250","RPV_2024_Par-ct-300-MChi-200-MSmu-250","RPV_2024_Par-ct-300-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-300-MChi-180-MSmu-300","RPV_2024_Par-ct-300-MChi-200-MSmu-300","RPV_2024_Par-ct-300-MChi-250-MSmu-300","RPV_2024_Par-ct-300-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-300-MChi-180-MSmu-350","RPV_2024_Par-ct-300-MChi-200-MSmu-350","RPV_2024_Par-ct-300-MChi-250-MSmu-350","RPV_2024_Par-ct-300-MChi-300-MSmu-350","RPV_2024_Par-ct-300-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-300-MChi-180-MSmu-400","RPV_2024_Par-ct-300-MChi-200-MSmu-400","RPV_2024_Par-ct-300-MChi-250-MSmu-400","RPV_2024_Par-ct-300-MChi-300-MSmu-400","RPV_2024_Par-ct-300-MChi-350-MSmu-400","RPV_2024_Par-ct-300-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-300-MChi-180-MSmu-450","RPV_2024_Par-ct-300-MChi-200-MSmu-450","RPV_2024_Par-ct-300-MChi-250-MSmu-450","RPV_2024_Par-ct-300-MChi-300-MSmu-450","RPV_2024_Par-ct-300-MChi-350-MSmu-450","RPV_2024_Par-ct-300-MChi-400-MSmu-450","RPV_2024_Par-ct-300-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-300-MChi-180-MSmu-500","RPV_2024_Par-ct-300-MChi-200-MSmu-500","RPV_2024_Par-ct-300-MChi-250-MSmu-500","RPV_2024_Par-ct-300-MChi-300-MSmu-500","RPV_2024_Par-ct-300-MChi-350-MSmu-500","RPV_2024_Par-ct-300-MChi-400-MSmu-500","RPV_2024_Par-ct-300-MChi-450-MSmu-500","RPV_2024_Par-ct-300-MChi-480-MSmu-500"

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
         "RPV_2024_Par-ct-1000-MChi-180-MSmu-200","RPV_2024_Par-ct-1000-MChi-180-MSmu-250","RPV_2024_Par-ct-1000-MChi-200-MSmu-250","RPV_2024_Par-ct-1000-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-1000-MChi-180-MSmu-300","RPV_2024_Par-ct-1000-MChi-200-MSmu-300","RPV_2024_Par-ct-1000-MChi-250-MSmu-300","RPV_2024_Par-ct-1000-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-1000-MChi-180-MSmu-350","RPV_2024_Par-ct-1000-MChi-200-MSmu-350","RPV_2024_Par-ct-1000-MChi-250-MSmu-350","RPV_2024_Par-ct-1000-MChi-300-MSmu-350","RPV_2024_Par-ct-1000-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-1000-MChi-180-MSmu-400","RPV_2024_Par-ct-1000-MChi-200-MSmu-400","RPV_2024_Par-ct-1000-MChi-250-MSmu-400","RPV_2024_Par-ct-1000-MChi-300-MSmu-400","RPV_2024_Par-ct-1000-MChi-350-MSmu-400","RPV_2024_Par-ct-1000-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-1000-MChi-180-MSmu-450","RPV_2024_Par-ct-1000-MChi-200-MSmu-450","RPV_2024_Par-ct-1000-MChi-250-MSmu-450","RPV_2024_Par-ct-1000-MChi-300-MSmu-450","RPV_2024_Par-ct-1000-MChi-350-MSmu-450","RPV_2024_Par-ct-1000-MChi-400-MSmu-450","RPV_2024_Par-ct-1000-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-1000-MChi-180-MSmu-500","RPV_2024_Par-ct-1000-MChi-200-MSmu-500","RPV_2024_Par-ct-1000-MChi-250-MSmu-500","RPV_2024_Par-ct-1000-MChi-300-MSmu-500","RPV_2024_Par-ct-1000-MChi-350-MSmu-500","RPV_2024_Par-ct-1000-MChi-400-MSmu-500","RPV_2024_Par-ct-1000-MChi-450-MSmu-500","RPV_2024_Par-ct-1000-MChi-480-MSmu-500"

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
