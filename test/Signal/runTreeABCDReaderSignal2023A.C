{
//   gROOT->ProcessLine(".L ../HistogramManager.C+");
//  gROOT->ProcessLine(".L TreeABCDReader.C+");
  std::vector<TString > systlist;
  //   systlist.push_back("");
  // systlist.push_back("Lumi2023AUp");
  // systlist.push_back("Lumi2023ADown");
  // systlist.push_back("L12023AUp");
  // systlist.push_back("L12023ADown");
  // systlist.push_back("Trigger2023ADown");
  // systlist.push_back("Trigger2023AUp");
  // systlist.push_back("MuonID2023AUp");
  // systlist.push_back("MuonID2023ADown");
  // systlist.push_back("MuonISO2023AUp");
  // systlist.push_back("MuonISO2023ADown");
  // systlist.push_back("PU2023AUp");
  // systlist.push_back("PU2023ADown");
  systlist.push_back("JEC2023AUp");
  systlist.push_back("JEC2023ADown");
  systlist.push_back("JER2023AUp");
  systlist.push_back("JER2023ADown");
  // systlist.push_back("Vtx2023AUp");
  // systlist.push_back("Vtx2023ADown");
  // // // systlist.push_back("RoccorUp");
  // // // systlist.push_back("RoccorDown");
  //   systlist.push_back("XSUp");
  // systlist.push_back("XSDown");
  // systlist.push_back("TopPtUp");
  // systlist.push_back("TopPtDown");


  // systlist.push_back("PDFUp");
  // systlist.push_back("PDFDown");
  // systlist.push_back("ScaleUp");
  // systlist.push_back("ScaleDown");


  int YEAR = 2023;
  TString Year = "2023A";
  TString Prod = "RPV_2023A"; //DATAMC2018_EMU_10_06_2024 //  EMU: SYST_EMU_CTAU100/Prod100_EMU_JER_up !! M: SYST_CTAU100/Prod100_JER_up // BKG MC and Data PROD_ANNIVERSAIRE_2024
  bool isPostAPV = false;
  bool Signal = true;
  bool SameSign = false;
  bool Forward = false; 
  bool DoubleMuon = true;
bool CorrectCorrelation = false;
  bool isMC = true;
  int mixing = 0; //-1 : Full left, 0 LR+RL , 1 Full RIght
  int Channel = 2; // 0 : EMu, 1 : SM, 2 : DM
  
  TChain c("ttree");
//  // Signal
   

      TString SignalSet001[34]={         "RPV_2023A_Msmu-200_Mchi-180_ct-001","RPV_2023A_Msmu-250_Mchi-180_ct-001","RPV_2023A_Msmu-250_Mchi-200_ct-001",
"RPV_2023A_Msmu-250_Mchi-230_ct-001","RPV_2023A_Msmu-300_Mchi-180_ct-001","RPV_2023A_Msmu-300_Mchi-200_ct-001","RPV_2023A_Msmu-300_Mchi-250_ct-001","RPV_2023A_Msmu-300_Mchi-280_ct-001",
"RPV_2023A_Msmu-350_Mchi-180_ct-001","RPV_2023A_Msmu-350_Mchi-200_ct-001","RPV_2023A_Msmu-350_Mchi-250_ct-001","RPV_2023A_Msmu-350_Mchi-300_ct-001","RPV_2023A_Msmu-350_Mchi-330_ct-001",
"RPV_2023A_Msmu-400_Mchi-180_ct-001","RPV_2023A_Msmu-400_Mchi-200_ct-001","RPV_2023A_Msmu-400_Mchi-250_ct-001","RPV_2023A_Msmu-400_Mchi-300_ct-001","RPV_2023A_Msmu-400_Mchi-350_ct-001",
"RPV_2023A_Msmu-400_Mchi-380_ct-001","RPV_2023A_Msmu-450_Mchi-180_ct-001","RPV_2023A_Msmu-450_Mchi-200_ct-001","RPV_2023A_Msmu-450_Mchi-250_ct-001","RPV_2023A_Msmu-450_Mchi-300_ct-001",
"RPV_2023A_Msmu-450_Mchi-350_ct-001","RPV_2023A_Msmu-450_Mchi-400_ct-001","RPV_2023A_Msmu-450_Mchi-430_ct-001","RPV_2023A_Msmu-500_Mchi-180_ct-001",
"RPV_2023A_Msmu-500_Mchi-200_ct-001"
,
"RPV_2023A_Msmu-500_Mchi-250_ct-001","RPV_2023A_Msmu-500_Mchi-300_ct-001","RPV_2023A_Msmu-500_Mchi-350_ct-001","RPV_2023A_Msmu-500_Mchi-400_ct-001","RPV_2023A_Msmu-500_Mchi-450_ct-001",
"RPV_2023A_Msmu-500_Mchi-480_ct-001"
}; 
//  //issue  with RPV_"+Year+"_smu500_neu400_ctau001 500 450 001
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
            Prod = "RPV_2023A";
            if (systlist.size() == 0) return;
            else if (systlist[j] == "JEC"+Year+"Up")
              {
                Prod = "RPV_2023A_JECUp";
              }
            else if (systlist[j] == "JEC"+Year+"Down")
              {
                Prod = "RPV_2023A_JECDown";
              }
            else if (systlist[j] == "JER"+Year+"Up")
              {
                Prod = "RPV_2023A_JERUp";
              }
            else if (systlist[j] == "JER"+Year+"Down")
              {
                Prod = "RPV_2023A_JERDown";
              }
            else
              {
                Prod = "RPV_2023A";
              }
            TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet001[i]+".root";
            c.Reset();
            c.Add(Path);
            TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet001[i],systlist[j]);
            float mean = 1.;//t->MeanGenWeight()
            t->Loop(isMC, Prod, SignalSet001[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
            delete t;
        }
    }



      TString SignalSet003[34]={         "RPV_2023A_Msmu-200_Mchi-180_ct-003","RPV_2023A_Msmu-250_Mchi-180_ct-003","RPV_2023A_Msmu-250_Mchi-200_ct-003",
"RPV_2023A_Msmu-250_Mchi-230_ct-003","RPV_2023A_Msmu-300_Mchi-180_ct-003","RPV_2023A_Msmu-300_Mchi-200_ct-003","RPV_2023A_Msmu-300_Mchi-250_ct-003","RPV_2023A_Msmu-300_Mchi-280_ct-003",
"RPV_2023A_Msmu-350_Mchi-180_ct-003","RPV_2023A_Msmu-350_Mchi-200_ct-003","RPV_2023A_Msmu-350_Mchi-250_ct-003","RPV_2023A_Msmu-350_Mchi-300_ct-003","RPV_2023A_Msmu-350_Mchi-330_ct-003",
"RPV_2023A_Msmu-400_Mchi-180_ct-003","RPV_2023A_Msmu-400_Mchi-200_ct-003","RPV_2023A_Msmu-400_Mchi-250_ct-003","RPV_2023A_Msmu-400_Mchi-300_ct-003","RPV_2023A_Msmu-400_Mchi-350_ct-003",
"RPV_2023A_Msmu-400_Mchi-380_ct-003","RPV_2023A_Msmu-450_Mchi-180_ct-003","RPV_2023A_Msmu-450_Mchi-200_ct-003","RPV_2023A_Msmu-450_Mchi-250_ct-003","RPV_2023A_Msmu-450_Mchi-300_ct-003",
"RPV_2023A_Msmu-450_Mchi-350_ct-003","RPV_2023A_Msmu-450_Mchi-400_ct-003","RPV_2023A_Msmu-450_Mchi-430_ct-003","RPV_2023A_Msmu-500_Mchi-180_ct-003",
"RPV_2023A_Msmu-500_Mchi-200_ct-003"
,
"RPV_2023A_Msmu-500_Mchi-250_ct-003","RPV_2023A_Msmu-500_Mchi-300_ct-003","RPV_2023A_Msmu-500_Mchi-350_ct-003","RPV_2023A_Msmu-500_Mchi-400_ct-003","RPV_2023A_Msmu-500_Mchi-450_ct-003",
"RPV_2023A_Msmu-500_Mchi-480_ct-003"
};

// // //  issue RPV_"+Year+"_smu500_neu400_ctau003
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2023A";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2023A";
                  }
              TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet003[i]+".root";
                c.Reset();
                c.Add(Path);
                TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet003[i],systlist[j]);
                float mean = 1.;//t->MeanGenWeight()
                t->Loop(isMC, Prod, SignalSet003[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
                delete t;
        }
    }

      TString SignalSet010[34]={         "RPV_2023A_Msmu-200_Mchi-180_ct-010","RPV_2023A_Msmu-250_Mchi-180_ct-010","RPV_2023A_Msmu-250_Mchi-200_ct-010",
"RPV_2023A_Msmu-250_Mchi-230_ct-010","RPV_2023A_Msmu-300_Mchi-180_ct-010","RPV_2023A_Msmu-300_Mchi-200_ct-010","RPV_2023A_Msmu-300_Mchi-250_ct-010","RPV_2023A_Msmu-300_Mchi-280_ct-010",
"RPV_2023A_Msmu-350_Mchi-180_ct-010","RPV_2023A_Msmu-350_Mchi-200_ct-010","RPV_2023A_Msmu-350_Mchi-250_ct-010","RPV_2023A_Msmu-350_Mchi-300_ct-010","RPV_2023A_Msmu-350_Mchi-330_ct-010",
"RPV_2023A_Msmu-400_Mchi-180_ct-010","RPV_2023A_Msmu-400_Mchi-200_ct-010","RPV_2023A_Msmu-400_Mchi-250_ct-010","RPV_2023A_Msmu-400_Mchi-300_ct-010","RPV_2023A_Msmu-400_Mchi-350_ct-010",
"RPV_2023A_Msmu-400_Mchi-380_ct-010","RPV_2023A_Msmu-450_Mchi-180_ct-010","RPV_2023A_Msmu-450_Mchi-200_ct-010","RPV_2023A_Msmu-450_Mchi-250_ct-010","RPV_2023A_Msmu-450_Mchi-300_ct-010",
"RPV_2023A_Msmu-450_Mchi-350_ct-010","RPV_2023A_Msmu-450_Mchi-400_ct-010","RPV_2023A_Msmu-450_Mchi-430_ct-010","RPV_2023A_Msmu-500_Mchi-180_ct-010",
"RPV_2023A_Msmu-500_Mchi-200_ct-010"
,
"RPV_2023A_Msmu-500_Mchi-250_ct-010","RPV_2023A_Msmu-500_Mchi-300_ct-010","RPV_2023A_Msmu-500_Mchi-350_ct-010","RPV_2023A_Msmu-500_Mchi-400_ct-010","RPV_2023A_Msmu-500_Mchi-450_ct-010",
"RPV_2023A_Msmu-500_Mchi-480_ct-010"};
 //issue with 500 350 10 // 500 400 10
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2023A";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2023A";
                  }
          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet010[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet010[i],systlist[j]);
      float mean = 1.;//t->MeanGenWeight()
      t->Loop(isMC, Prod, SignalSet010[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
      delete t;
        }
    }


      TString SignalSet030[34]={       "RPV_2023A_Msmu-200_Mchi-180_ct-030","RPV_2023A_Msmu-250_Mchi-180_ct-030","RPV_2023A_Msmu-250_Mchi-200_ct-030",
"RPV_2023A_Msmu-250_Mchi-230_ct-030","RPV_2023A_Msmu-300_Mchi-180_ct-030","RPV_2023A_Msmu-300_Mchi-200_ct-030","RPV_2023A_Msmu-300_Mchi-250_ct-030","RPV_2023A_Msmu-300_Mchi-280_ct-030",
"RPV_2023A_Msmu-350_Mchi-180_ct-030","RPV_2023A_Msmu-350_Mchi-200_ct-030","RPV_2023A_Msmu-350_Mchi-250_ct-030","RPV_2023A_Msmu-350_Mchi-300_ct-030","RPV_2023A_Msmu-350_Mchi-330_ct-030",
"RPV_2023A_Msmu-400_Mchi-180_ct-030","RPV_2023A_Msmu-400_Mchi-200_ct-030","RPV_2023A_Msmu-400_Mchi-250_ct-030","RPV_2023A_Msmu-400_Mchi-300_ct-030","RPV_2023A_Msmu-400_Mchi-350_ct-030",
"RPV_2023A_Msmu-400_Mchi-380_ct-030","RPV_2023A_Msmu-450_Mchi-180_ct-030","RPV_2023A_Msmu-450_Mchi-200_ct-030","RPV_2023A_Msmu-450_Mchi-250_ct-030","RPV_2023A_Msmu-450_Mchi-300_ct-030",
"RPV_2023A_Msmu-450_Mchi-350_ct-030","RPV_2023A_Msmu-450_Mchi-400_ct-030","RPV_2023A_Msmu-450_Mchi-430_ct-030","RPV_2023A_Msmu-500_Mchi-180_ct-030",
"RPV_2023A_Msmu-500_Mchi-200_ct-030"
,
"RPV_2023A_Msmu-500_Mchi-250_ct-030","RPV_2023A_Msmu-500_Mchi-300_ct-030","RPV_2023A_Msmu-500_Mchi-350_ct-030","RPV_2023A_Msmu-500_Mchi-400_ct-030","RPV_2023A_Msmu-500_Mchi-450_ct-030",
"RPV_2023A_Msmu-500_Mchi-480_ct-030"};
 
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2023A";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2023A";
                  }
          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet030[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet030[i],systlist[j]);
      float mean = 1.;//t->MeanGenWeight()
      t->Loop(isMC, Prod, SignalSet030[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
          delete t;
        }
    }

  TString SignalSet100[34]=      {"RPV_2023A_Msmu-200_Mchi-180_ct-100","RPV_2023A_Msmu-250_Mchi-180_ct-100","RPV_2023A_Msmu-250_Mchi-200_ct-100",
"RPV_2023A_Msmu-250_Mchi-230_ct-100","RPV_2023A_Msmu-300_Mchi-180_ct-100","RPV_2023A_Msmu-300_Mchi-200_ct-100","RPV_2023A_Msmu-300_Mchi-250_ct-100","RPV_2023A_Msmu-300_Mchi-280_ct-100",
"RPV_2023A_Msmu-350_Mchi-180_ct-100","RPV_2023A_Msmu-350_Mchi-200_ct-100","RPV_2023A_Msmu-350_Mchi-250_ct-100","RPV_2023A_Msmu-350_Mchi-300_ct-100","RPV_2023A_Msmu-350_Mchi-330_ct-100",
"RPV_2023A_Msmu-400_Mchi-180_ct-100","RPV_2023A_Msmu-400_Mchi-200_ct-100","RPV_2023A_Msmu-400_Mchi-250_ct-100","RPV_2023A_Msmu-400_Mchi-300_ct-100","RPV_2023A_Msmu-400_Mchi-350_ct-100",
"RPV_2023A_Msmu-400_Mchi-380_ct-100","RPV_2023A_Msmu-450_Mchi-180_ct-100","RPV_2023A_Msmu-450_Mchi-200_ct-100","RPV_2023A_Msmu-450_Mchi-250_ct-100","RPV_2023A_Msmu-450_Mchi-300_ct-100",
"RPV_2023A_Msmu-450_Mchi-350_ct-100","RPV_2023A_Msmu-450_Mchi-400_ct-100","RPV_2023A_Msmu-450_Mchi-430_ct-100","RPV_2023A_Msmu-500_Mchi-180_ct-100",
"RPV_2023A_Msmu-500_Mchi-200_ct-100"
,
"RPV_2023A_Msmu-500_Mchi-250_ct-100","RPV_2023A_Msmu-500_Mchi-300_ct-100","RPV_2023A_Msmu-500_Mchi-350_ct-100","RPV_2023A_Msmu-500_Mchi-400_ct-100","RPV_2023A_Msmu-500_Mchi-450_ct-100",
"RPV_2023A_Msmu-500_Mchi-480_ct-100"};
 
 for (int i = 0 ; i< 34 ; i++) 
    {
            for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2023A";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2023A";
                  }
          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet100[i]+".root";
          c.Reset();
          c.Add(Path);
          TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet100[i],systlist[j]);
          float mean = 1.;//t->MeanGenWeight()
          t->Loop(isMC, Prod, SignalSet100[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
              delete t;
        }
    }

      TString SignalSet300[34]={         "RPV_2023A_Msmu-200_Mchi-180_ct-300","RPV_2023A_Msmu-250_Mchi-180_ct-300","RPV_2023A_Msmu-250_Mchi-200_ct-300",
"RPV_2023A_Msmu-250_Mchi-230_ct-300","RPV_2023A_Msmu-300_Mchi-180_ct-300","RPV_2023A_Msmu-300_Mchi-200_ct-300","RPV_2023A_Msmu-300_Mchi-250_ct-300","RPV_2023A_Msmu-300_Mchi-280_ct-300",
"RPV_2023A_Msmu-350_Mchi-180_ct-300","RPV_2023A_Msmu-350_Mchi-200_ct-300","RPV_2023A_Msmu-350_Mchi-250_ct-300","RPV_2023A_Msmu-350_Mchi-300_ct-300","RPV_2023A_Msmu-350_Mchi-330_ct-300",
"RPV_2023A_Msmu-400_Mchi-180_ct-300","RPV_2023A_Msmu-400_Mchi-200_ct-300","RPV_2023A_Msmu-400_Mchi-250_ct-300","RPV_2023A_Msmu-400_Mchi-300_ct-300","RPV_2023A_Msmu-400_Mchi-350_ct-300",
"RPV_2023A_Msmu-400_Mchi-380_ct-300","RPV_2023A_Msmu-450_Mchi-180_ct-300","RPV_2023A_Msmu-450_Mchi-200_ct-300","RPV_2023A_Msmu-450_Mchi-250_ct-300","RPV_2023A_Msmu-450_Mchi-300_ct-300",
"RPV_2023A_Msmu-450_Mchi-350_ct-300","RPV_2023A_Msmu-450_Mchi-400_ct-300","RPV_2023A_Msmu-450_Mchi-430_ct-300","RPV_2023A_Msmu-500_Mchi-180_ct-300",
"RPV_2023A_Msmu-500_Mchi-200_ct-300"
,
"RPV_2023A_Msmu-500_Mchi-250_ct-300","RPV_2023A_Msmu-500_Mchi-300_ct-300","RPV_2023A_Msmu-500_Mchi-350_ct-300","RPV_2023A_Msmu-500_Mchi-400_ct-300","RPV_2023A_Msmu-500_Mchi-450_ct-300",
"RPV_2023A_Msmu-500_Mchi-480_ct-300"};
 
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2023A";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2023A";
                  }
            TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet300[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet300[i],systlist[j]);
      float mean = 1.;//t->MeanGenWeight()
      t->Loop(isMC, Prod, SignalSet300[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
delete t;
        }
    }



          TString SignalSet1000[34]={            "RPV_2023A_Msmu-200_Mchi-180_ct-1000","RPV_2023A_Msmu-250_Mchi-180_ct-1000","RPV_2023A_Msmu-250_Mchi-200_ct-1000",
"RPV_2023A_Msmu-250_Mchi-230_ct-1000","RPV_2023A_Msmu-300_Mchi-180_ct-1000","RPV_2023A_Msmu-300_Mchi-200_ct-1000","RPV_2023A_Msmu-300_Mchi-250_ct-1000","RPV_2023A_Msmu-300_Mchi-280_ct-1000",
"RPV_2023A_Msmu-350_Mchi-180_ct-1000","RPV_2023A_Msmu-350_Mchi-200_ct-1000","RPV_2023A_Msmu-350_Mchi-250_ct-1000","RPV_2023A_Msmu-350_Mchi-300_ct-1000","RPV_2023A_Msmu-350_Mchi-330_ct-1000",
"RPV_2023A_Msmu-400_Mchi-180_ct-1000","RPV_2023A_Msmu-400_Mchi-200_ct-1000","RPV_2023A_Msmu-400_Mchi-250_ct-1000","RPV_2023A_Msmu-400_Mchi-300_ct-1000","RPV_2023A_Msmu-400_Mchi-350_ct-1000",
"RPV_2023A_Msmu-400_Mchi-380_ct-1000","RPV_2023A_Msmu-450_Mchi-180_ct-1000","RPV_2023A_Msmu-450_Mchi-200_ct-1000","RPV_2023A_Msmu-450_Mchi-250_ct-1000","RPV_2023A_Msmu-450_Mchi-300_ct-1000",
"RPV_2023A_Msmu-450_Mchi-350_ct-1000","RPV_2023A_Msmu-450_Mchi-400_ct-1000","RPV_2023A_Msmu-450_Mchi-430_ct-1000","RPV_2023A_Msmu-500_Mchi-180_ct-1000",
"RPV_2023A_Msmu-500_Mchi-200_ct-1000"
,
"RPV_2023A_Msmu-500_Mchi-250_ct-1000","RPV_2023A_Msmu-500_Mchi-300_ct-1000","RPV_2023A_Msmu-500_Mchi-350_ct-1000","RPV_2023A_Msmu-500_Mchi-400_ct-1000","RPV_2023A_Msmu-500_Mchi-450_ct-1000",
"RPV_2023A_Msmu-500_Mchi-480_ct-1000"
};
 
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2023A";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2023A_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2023A_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2023A";
                  }
            TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet1000[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet1000[i],systlist[j]);
      float mean = 1.;//t->MeanGenWeight()
      t->Loop(isMC, Prod, SignalSet1000[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
            delete t;
    }
}



//       TString SignalSetM[34]={"RPV_"+Year+"_smu200_neu180","RPV_"+Year+"_smu250_neu180","RPV_"+Year+"_smu250_neu200",
// "RPV_"+Year+"_smu250_neu230","RPV_"+Year+"_smu300_neu180","RPV_"+Year+"_smu300_neu200","RPV_"+Year+"_smu300_neu250","RPV_"+Year+"_smu300_neu280",
// "RPV_"+Year+"_smu350_neu180","RPV_"+Year+"_smu350_neu200","RPV_"+Year+"_smu350_neu250","RPV_"+Year+"_smu350_neu300","RPV_"+Year+"_smu350_neu330",
// "RPV_"+Year+"_smu400_neu180","RPV_"+Year+"_smu400_neu200","RPV_"+Year+"_smu400_neu250","RPV_"+Year+"_smu400_neu300","RPV_"+Year+"_smu400_neu350",
// "RPV_"+Year+"_smu400_neu380","RPV_"+Year+"_smu450_neu180","RPV_"+Year+"_smu450_neu200","RPV_"+Year+"_smu450_neu250","RPV_"+Year+"_smu450_neu300",
// "RPV_"+Year+"_smu450_neu350","RPV_"+Year+"_smu450_neu400","RPV_"+Year+"_smu450_neu430","RPV_"+Year+"_smu500_neu180","RPV_"+Year+"_smu500_neu200",
// "RPV_"+Year+"_smu500_neu250","RPV_"+Year+"_smu500_neu300","RPV_"+Year+"_smu500_neu350","RPV_"+Year+"_smu500_neu400","RPV_"+Year+"_smu500_neu450",
// "RPV_"+Year+"_smu500_neu480"};
 
//  for (int i = 0 ; i< 34 ; i++) 
//     {
//             for (unsigned int j = 0 ; j < systlist.size(); j++)
//         {
//                 Prod = "RPV_2023A";
//                 if (systlist.size() == 0) return;
//                 else if (systlist[j] == "JECUp")
//                   {
//                     Prod = "RPV_2023A_JECUp";
//                   }
//                 else if (systlist[j] == "JECDown")
//                   {
//                     Prod = "RPV_2023A_JECDown";
//                   }
//                 else if (systlist[j] == "JERUp")
//                   {
//                     Prod = "RPV_2023A_JERUp";
//                   }
//                 else if (systlist[j] == "JERDown")
//                   {
//                     Prod = "RPV_2023A_JERDown";
//                   }
//                 else
//                   {
//                     Prod = "RPV_2023A";
//                   }
//             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSetM[i]+".root";
//             c.Reset();
//             c.Add(Path);
//             TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSetM[i],systlist[j]);
//             float mean = 1.;//t->MeanGenWeight()
//             t->Loop(isMC, Prod, SignalSetM[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j]);
//             delete t;
//         }
//     }




//       TString SignalSetctau[7]={"RPV_ctau001","RPV_ctau003","RPV_ctau010",
// "RPV_ctau030","RPV_ctau100","RPV_ctau300","RPV_ctau1000"
// };
 
//  for (int i = 0 ; i< 7 ; i++) 
//     {
//             for (unsigned int j = 0 ; j < systlist.size(); j++)
//         {
//                 Prod = "RPV_2023A";
//                 if (systlist.size() == 0) return;
//                 else if (systlist[j] == "JECUp")
//                   {
//                     Prod = "RPV_2023A_JECUp";
//                   }
//                 else if (systlist[j] == "JECDown")
//                   {
//                     Prod = "RPV_2023A_JECDown";
//                   }
//                 else if (systlist[j] == "JERUp")
//                   {
//                     Prod = "RPV_2023A_JERUp";
//                   }
//                 else if (systlist[j] == "JERDown")
//                   {
//                     Prod = "RPV_2023A_JERDown";
//                   }
//                 else
//                   {
//                     Prod = "RPV_2023A";
//                   }
//             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSetctau[i]+".root";
//             c.Reset();
//             c.Add(Path);
//             TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSetctau[i],systlist[j]);
//             float mean = 1.;//t->MeanGenWeight()
//             t->Loop(isMC, Prod, SignalSetctau[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j]);
//             delete t;
//         }
//     }


}
