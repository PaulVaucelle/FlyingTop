{
//   gROOT->ProcessLine(".L ../HistogramManager.C+");
//  gROOT->ProcessLine(".L TreeABCDReader.C+");
  std::vector<TString > systlist;
  // systlist.push_back("");
  // systlist.push_back("Lumi2017Up");
  // systlist.push_back("Lumi2017Down");
  // systlist.push_back("L12017Up");
  // systlist.push_back("L12017Down");
  // systlist.push_back("Trigger2017Down");
  // systlist.push_back("Trigger2017Up");
  // systlist.push_back("MuonID2017Up");
  // systlist.push_back("MuonID2017Down");
  // systlist.push_back("MuonISO2017Up");
  // systlist.push_back("MuonISO2017Down");
  // systlist.push_back("PU2017Up");
  // systlist.push_back("PU2017Down");
  systlist.push_back("JEC2017Up");
  systlist.push_back("JEC2017Down");
  systlist.push_back("JER2017Up");
  systlist.push_back("JER2017Down");
  // systlist.push_back("Vtx2017Up");
  // systlist.push_back("Vtx2017Down");
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


  int YEAR = 2017;
  TString Year = "2017";
  TString Prod = "Signal_2017"; //DATAMC2018_EMU_10_06_2024 //  EMU: SYST_EMU_CTAU100/Prod100_EMU_JER_up !! M: SYST_CTAU100/Prod100_JER_up // BKG MC and Data PROD_ANNIVERSAIRE_2024
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
   

      TString SignalSet001[34]={"RPV_"+Year+"_smu200_neu180_ctau001","RPV_"+Year+"_smu250_neu180_ctau001","RPV_"+Year+"_smu250_neu200_ctau001",
"RPV_"+Year+"_smu250_neu230_ctau001","RPV_"+Year+"_smu300_neu180_ctau001","RPV_"+Year+"_smu300_neu200_ctau001","RPV_"+Year+"_smu300_neu280_ctau001","RPV_"+Year+"_smu300_neu250_ctau001",
"RPV_"+Year+"_smu350_neu180_ctau001","RPV_"+Year+"_smu350_neu200_ctau001","RPV_"+Year+"_smu350_neu250_ctau001","RPV_"+Year+"_smu350_neu300_ctau001","RPV_"+Year+"_smu350_neu330_ctau001",
"RPV_"+Year+"_smu400_neu180_ctau001","RPV_"+Year+"_smu400_neu200_ctau001","RPV_"+Year+"_smu400_neu250_ctau001","RPV_"+Year+"_smu400_neu300_ctau001","RPV_"+Year+"_smu400_neu350_ctau001",
"RPV_"+Year+"_smu400_neu380_ctau001","RPV_"+Year+"_smu450_neu180_ctau001","RPV_"+Year+"_smu450_neu200_ctau001","RPV_"+Year+"_smu450_neu250_ctau001","RPV_"+Year+"_smu450_neu300_ctau001",
"RPV_"+Year+"_smu450_neu350_ctau001","RPV_"+Year+"_smu450_neu400_ctau001","RPV_"+Year+"_smu450_neu430_ctau001","RPV_"+Year+"_smu500_neu180_ctau001","RPV_"+Year+"_smu500_neu200_ctau001",
"RPV_"+Year+"_smu500_neu250_ctau001","RPV_"+Year+"_smu500_neu300_ctau001","RPV_"+Year+"_smu500_neu350_ctau001","RPV_"+Year+"_smu500_neu400_ctau001","RPV_"+Year+"_smu500_neu450_ctau001",
"RPV_"+Year+"_smu500_neu480_ctau001"
}; 
//  //issue  with RPV_"+Year+"_smu500_neu400_ctau001 500 450 001
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
            Prod = "Signal_2017";
            if (systlist.size() == 0) return;
            else if (systlist[j] == "JEC2017Up")
              {
                Prod = "Signal_2017_JECUp";
              }
            else if (systlist[j] == "JEC2017Down")
              {
                Prod = "Signal_2017_JECDown";
              }
            else if (systlist[j] == "JER2017Up")
              {
                Prod = "Signal_2017_JERUp";
              }
            else if (systlist[j] == "JER2017Down")
              {
                Prod = "Signal_2017_JERDown";
              }
            else
              {
                Prod = "Signal_2017";
              }
            TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet001[i]+".root";
            c.Reset();
            c.Add(Path);
            TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet001[i],systlist[j]);
            float mean = 1.;//t->MeanGenWeight()
            t->Loop(isMC, Prod, SignalSet001[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
            delete t;
        }
    }



      TString SignalSet003[34]={
        "RPV_"+Year+"_smu200_neu180_ctau003","RPV_"+Year+"_smu250_neu180_ctau003","RPV_"+Year+"_smu250_neu200_ctau003",
"RPV_"+Year+"_smu250_neu230_ctau003","RPV_"+Year+"_smu300_neu180_ctau003","RPV_"+Year+"_smu300_neu200_ctau003","RPV_"+Year+"_smu300_neu280_ctau003","RPV_"+Year+"_smu300_neu250_ctau003",
"RPV_"+Year+"_smu350_neu180_ctau003","RPV_"+Year+"_smu350_neu200_ctau003","RPV_"+Year+"_smu350_neu250_ctau003","RPV_"+Year+"_smu350_neu300_ctau003","RPV_"+Year+"_smu350_neu330_ctau003",
"RPV_"+Year+"_smu400_neu180_ctau003","RPV_"+Year+"_smu400_neu200_ctau003","RPV_"+Year+"_smu400_neu250_ctau003","RPV_"+Year+"_smu400_neu300_ctau003","RPV_"+Year+"_smu400_neu350_ctau003",
"RPV_"+Year+"_smu400_neu380_ctau003","RPV_"+Year+"_smu450_neu180_ctau003","RPV_"+Year+"_smu450_neu200_ctau003","RPV_"+Year+"_smu450_neu250_ctau003","RPV_"+Year+"_smu450_neu300_ctau003",
"RPV_"+Year+"_smu450_neu350_ctau003","RPV_"+Year+"_smu450_neu400_ctau003","RPV_"+Year+"_smu450_neu430_ctau003","RPV_"+Year+"_smu500_neu180_ctau003","RPV_"+Year+"_smu500_neu200_ctau003",
"RPV_"+Year+"_smu500_neu250_ctau003","RPV_"+Year+"_smu500_neu300_ctau003","RPV_"+Year+"_smu500_neu350_ctau003","RPV_"+Year+"_smu500_neu400_ctau003","RPV_"+Year+"_smu500_neu450_ctau003",
"RPV_"+Year+"_smu500_neu480_ctau003"
};

// // //  issue RPV_"+Year+"_smu500_neu400_ctau003
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "Signal_2017";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC2017Up")
                  {
                    Prod = "Signal_2017_JECUp";
                  }
                else if (systlist[j] == "JEC2017Down")
                  {
                    Prod = "Signal_2017_JECDown";
                  }
                else if (systlist[j] == "JER2017Up")
                  {
                    Prod = "Signal_2017_JERUp";
                  }
                else if (systlist[j] == "JER2017Down")
                  {
                    Prod = "Signal_2017_JERDown";
                  }
                else
                  {
                    Prod = "Signal_2017";
                  }
              TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet003[i]+".root";
                c.Reset();
                c.Add(Path);
                TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet003[i],systlist[j]);
                float mean = 1.;//t->MeanGenWeight()
                t->Loop(isMC, Prod, SignalSet003[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
                delete t;
        }
    }

      TString SignalSet010[34]={"RPV_"+Year+"_smu200_neu180_ctau010","RPV_"+Year+"_smu250_neu180_ctau010","RPV_"+Year+"_smu250_neu200_ctau010",
"RPV_"+Year+"_smu250_neu230_ctau010","RPV_"+Year+"_smu300_neu180_ctau010","RPV_"+Year+"_smu300_neu200_ctau010","RPV_"+Year+"_smu300_neu250_ctau010","RPV_"+Year+"_smu300_neu280_ctau010",
"RPV_"+Year+"_smu350_neu180_ctau010","RPV_"+Year+"_smu350_neu200_ctau010","RPV_"+Year+"_smu350_neu250_ctau010","RPV_"+Year+"_smu350_neu300_ctau010","RPV_"+Year+"_smu350_neu330_ctau010",
"RPV_"+Year+"_smu400_neu180_ctau010","RPV_"+Year+"_smu400_neu200_ctau010","RPV_"+Year+"_smu400_neu250_ctau010","RPV_"+Year+"_smu400_neu300_ctau010","RPV_"+Year+"_smu400_neu350_ctau010",
"RPV_"+Year+"_smu400_neu380_ctau010","RPV_"+Year+"_smu450_neu180_ctau010","RPV_"+Year+"_smu450_neu200_ctau010","RPV_"+Year+"_smu450_neu250_ctau010","RPV_"+Year+"_smu450_neu300_ctau010",
"RPV_"+Year+"_smu450_neu350_ctau010","RPV_"+Year+"_smu450_neu400_ctau010","RPV_"+Year+"_smu450_neu430_ctau010","RPV_"+Year+"_smu500_neu180_ctau010","RPV_"+Year+"_smu500_neu200_ctau010",
"RPV_"+Year+"_smu500_neu250_ctau010","RPV_"+Year+"_smu500_neu300_ctau010","RPV_"+Year+"_smu500_neu350_ctau010","RPV_"+Year+"_smu500_neu400_ctau010","RPV_"+Year+"_smu500_neu450_ctau010",
"RPV_"+Year+"_smu500_neu480_ctau010"};
 //issue with 500 350 10 // 500 400 10
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "Signal_2017";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC2017Up")
                  {
                    Prod = "Signal_2017_JECUp";
                  }
                else if (systlist[j] == "JEC2017Down")
                  {
                    Prod = "Signal_2017_JECDown";
                  }
                else if (systlist[j] == "JER2017Up")
                  {
                    Prod = "Signal_2017_JERUp";
                  }
                else if (systlist[j] == "JER2017Down")
                  {
                    Prod = "Signal_2017_JERDown";
                  }
                else
                  {
                    Prod = "Signal_2017";
                  }
          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet010[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet010[i],systlist[j]);
      float mean = 1.;//t->MeanGenWeight()
      t->Loop(isMC, Prod, SignalSet010[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
      delete t;
        }
    }


      TString SignalSet030[34]={"RPV_"+Year+"_smu200_neu180_ctau030","RPV_"+Year+"_smu250_neu180_ctau030","RPV_"+Year+"_smu250_neu200_ctau030",
"RPV_"+Year+"_smu250_neu230_ctau030","RPV_"+Year+"_smu300_neu180_ctau030","RPV_"+Year+"_smu300_neu200_ctau030","RPV_"+Year+"_smu300_neu250_ctau030","RPV_"+Year+"_smu300_neu280_ctau030",
"RPV_"+Year+"_smu350_neu180_ctau030","RPV_"+Year+"_smu350_neu200_ctau030","RPV_"+Year+"_smu350_neu250_ctau030","RPV_"+Year+"_smu350_neu300_ctau030","RPV_"+Year+"_smu350_neu330_ctau030",
"RPV_"+Year+"_smu400_neu180_ctau030","RPV_"+Year+"_smu400_neu200_ctau030","RPV_"+Year+"_smu400_neu250_ctau030","RPV_"+Year+"_smu400_neu300_ctau030","RPV_"+Year+"_smu400_neu350_ctau030",
"RPV_"+Year+"_smu400_neu380_ctau030","RPV_"+Year+"_smu450_neu180_ctau030","RPV_"+Year+"_smu450_neu200_ctau030","RPV_"+Year+"_smu450_neu250_ctau030","RPV_"+Year+"_smu450_neu300_ctau030",
"RPV_"+Year+"_smu450_neu350_ctau030","RPV_"+Year+"_smu450_neu400_ctau030","RPV_"+Year+"_smu450_neu430_ctau030","RPV_"+Year+"_smu500_neu180_ctau030","RPV_"+Year+"_smu500_neu200_ctau030",
"RPV_"+Year+"_smu500_neu250_ctau030","RPV_"+Year+"_smu500_neu300_ctau030","RPV_"+Year+"_smu500_neu350_ctau030","RPV_"+Year+"_smu500_neu400_ctau030","RPV_"+Year+"_smu500_neu450_ctau030",
"RPV_"+Year+"_smu500_neu480_ctau030"};
 
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "Signal_2017";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC2017Up")
                  {
                    Prod = "Signal_2017_JECUp";
                  }
                else if (systlist[j] == "JEC2017Down")
                  {
                    Prod = "Signal_2017_JECDown";
                  }
                else if (systlist[j] == "JER2017Up")
                  {
                    Prod = "Signal_2017_JERUp";
                  }
                else if (systlist[j] == "JER2017Down")
                  {
                    Prod = "Signal_2017_JERDown";
                  }
                else
                  {
                    Prod = "Signal_2017";
                  }
          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet030[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet030[i],systlist[j]);
      float mean = 1.;//t->MeanGenWeight()
      t->Loop(isMC, Prod, SignalSet030[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
          delete t;
        }
    }

  TString SignalSet100[34]={"RPV_"+Year+"_smu200_neu180_ctau100","RPV_"+Year+"_smu250_neu180_ctau100","RPV_"+Year+"_smu250_neu200_ctau100",
"RPV_"+Year+"_smu250_neu230_ctau100","RPV_"+Year+"_smu300_neu180_ctau100","RPV_"+Year+"_smu300_neu200_ctau100","RPV_"+Year+"_smu300_neu250_ctau100","RPV_"+Year+"_smu300_neu280_ctau100",
"RPV_"+Year+"_smu350_neu180_ctau100","RPV_"+Year+"_smu350_neu200_ctau100","RPV_"+Year+"_smu350_neu250_ctau100","RPV_"+Year+"_smu350_neu300_ctau100","RPV_"+Year+"_smu350_neu330_ctau100",
"RPV_"+Year+"_smu400_neu180_ctau100","RPV_"+Year+"_smu400_neu200_ctau100","RPV_"+Year+"_smu400_neu250_ctau100","RPV_"+Year+"_smu400_neu300_ctau100","RPV_"+Year+"_smu400_neu350_ctau100",
"RPV_"+Year+"_smu400_neu380_ctau100","RPV_"+Year+"_smu450_neu180_ctau100","RPV_"+Year+"_smu450_neu200_ctau100","RPV_"+Year+"_smu450_neu250_ctau100","RPV_"+Year+"_smu450_neu300_ctau100",
"RPV_"+Year+"_smu450_neu350_ctau100","RPV_"+Year+"_smu450_neu400_ctau100","RPV_"+Year+"_smu450_neu430_ctau100","RPV_"+Year+"_smu500_neu180_ctau100","RPV_"+Year+"_smu500_neu200_ctau100",
"RPV_"+Year+"_smu500_neu250_ctau100","RPV_"+Year+"_smu500_neu300_ctau100","RPV_"+Year+"_smu500_neu350_ctau100","RPV_"+Year+"_smu500_neu400_ctau100","RPV_"+Year+"_smu500_neu450_ctau100",
"RPV_"+Year+"_smu500_neu480_ctau100"};
 
 for (int i = 0 ; i< 34 ; i++) 
    {
            for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "Signal_2017";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC2017Up")
                  {
                    Prod = "Signal_2017_JECUp";
                  }
                else if (systlist[j] == "JEC2017Down")
                  {
                    Prod = "Signal_2017_JECDown";
                  }
                else if (systlist[j] == "JER2017Up")
                  {
                    Prod = "Signal_2017_JERUp";
                  }
                else if (systlist[j] == "JER2017Down")
                  {
                    Prod = "Signal_2017_JERDown";
                  }
                else
                  {
                    Prod = "Signal_2017";
                  }
          TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet100[i]+".root";
          c.Reset();
          c.Add(Path);
          TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet100[i],systlist[j]);
          float mean = 1.;//t->MeanGenWeight()
          t->Loop(isMC, Prod, SignalSet100[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
              delete t;
        }
    }

      TString SignalSet300[34]={"RPV_"+Year+"_smu200_neu180_ctau300","RPV_"+Year+"_smu250_neu180_ctau300","RPV_"+Year+"_smu250_neu200_ctau300",
"RPV_"+Year+"_smu250_neu230_ctau300","RPV_"+Year+"_smu300_neu180_ctau300","RPV_"+Year+"_smu300_neu200_ctau300","RPV_"+Year+"_smu300_neu250_ctau300","RPV_"+Year+"_smu300_neu280_ctau300",
"RPV_"+Year+"_smu350_neu180_ctau300","RPV_"+Year+"_smu350_neu200_ctau300","RPV_"+Year+"_smu350_neu250_ctau300","RPV_"+Year+"_smu350_neu300_ctau300","RPV_"+Year+"_smu350_neu330_ctau300",
"RPV_"+Year+"_smu400_neu180_ctau300","RPV_"+Year+"_smu400_neu200_ctau300","RPV_"+Year+"_smu400_neu250_ctau300","RPV_"+Year+"_smu400_neu300_ctau300","RPV_"+Year+"_smu400_neu350_ctau300",
"RPV_"+Year+"_smu400_neu380_ctau300","RPV_"+Year+"_smu450_neu180_ctau300","RPV_"+Year+"_smu450_neu200_ctau300","RPV_"+Year+"_smu450_neu250_ctau300","RPV_"+Year+"_smu450_neu300_ctau300",
"RPV_"+Year+"_smu450_neu350_ctau300","RPV_"+Year+"_smu450_neu400_ctau300","RPV_"+Year+"_smu450_neu430_ctau300","RPV_"+Year+"_smu500_neu180_ctau300","RPV_"+Year+"_smu500_neu200_ctau300",
"RPV_"+Year+"_smu500_neu250_ctau300","RPV_"+Year+"_smu500_neu300_ctau300","RPV_"+Year+"_smu500_neu350_ctau300","RPV_"+Year+"_smu500_neu400_ctau300","RPV_"+Year+"_smu500_neu450_ctau300",
"RPV_"+Year+"_smu500_neu480_ctau300"};
 
 for (int i = 0; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "Signal_2017";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC2017Up")
                  {
                    Prod = "Signal_2017_JECUp";
                  }
                else if (systlist[j] == "JEC2017Down")
                  {
                    Prod = "Signal_2017_JECDown";
                  }
                else if (systlist[j] == "JER2017Up")
                  {
                    Prod = "Signal_2017_JERUp";
                  }
                else if (systlist[j] == "JER2017Down")
                  {
                    Prod = "Signal_2017_JERDown";
                  }
                else
                  {
                    Prod = "Signal_2017";
                  }
            TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet300[i]+".root";
      c.Reset();
      c.Add(Path);
      TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSet300[i],systlist[j]);
      float mean = 1.;//t->MeanGenWeight()
      t->Loop(isMC, Prod, SignalSet300[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
delete t;
        }
    }



          TString SignalSet1000[34]={"RPV_"+Year+"_smu200_neu180_ctau1000","RPV_"+Year+"_smu250_neu180_ctau1000","RPV_"+Year+"_smu250_neu200_ctau1000",
"RPV_"+Year+"_smu250_neu230_ctau1000","RPV_"+Year+"_smu300_neu180_ctau1000","RPV_"+Year+"_smu300_neu200_ctau1000","RPV_"+Year+"_smu300_neu280_ctau1000","RPV_"+Year+"_smu300_neu250_ctau1000",
"RPV_"+Year+"_smu350_neu180_ctau1000","RPV_"+Year+"_smu350_neu200_ctau1000","RPV_"+Year+"_smu350_neu250_ctau1000","RPV_"+Year+"_smu350_neu300_ctau1000","RPV_"+Year+"_smu350_neu330_ctau1000",
"RPV_"+Year+"_smu400_neu180_ctau1000","RPV_"+Year+"_smu400_neu200_ctau1000","RPV_"+Year+"_smu400_neu250_ctau1000","RPV_"+Year+"_smu400_neu300_ctau1000","RPV_"+Year+"_smu400_neu350_ctau1000",
"RPV_"+Year+"_smu400_neu380_ctau1000","RPV_"+Year+"_smu450_neu180_ctau1000","RPV_"+Year+"_smu450_neu200_ctau1000","RPV_"+Year+"_smu450_neu250_ctau1000","RPV_"+Year+"_smu450_neu300_ctau1000",
"RPV_"+Year+"_smu450_neu350_ctau1000","RPV_"+Year+"_smu450_neu400_ctau1000","RPV_"+Year+"_smu450_neu430_ctau1000","RPV_"+Year+"_smu500_neu180_ctau1000","RPV_"+Year+"_smu500_neu200_ctau1000",
"RPV_"+Year+"_smu500_neu250_ctau1000","RPV_"+Year+"_smu500_neu300_ctau1000","RPV_"+Year+"_smu500_neu350_ctau1000","RPV_"+Year+"_smu500_neu400_ctau1000","RPV_"+Year+"_smu500_neu450_ctau1000",
"RPV_"+Year+"_smu500_neu480_ctau1000"
};
 
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "Signal_2017";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC2017Up")
                  {
                    Prod = "Signal_2017_JECUp";
                  }
                else if (systlist[j] == "JEC2017Down")
                  {
                    Prod = "Signal_2017_JECDown";
                  }
                else if (systlist[j] == "JER2017Up")
                  {
                    Prod = "Signal_2017_JERUp";
                  }
                else if (systlist[j] == "JER2017Down")
                  {
                    Prod = "Signal_2017_JERDown";
                  }
                else
                  {
                    Prod = "Signal_2017";
                  }
            TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSet1000[i]+".root";
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
//                 Prod = "Signal_2017";
//                 if (systlist.size() == 0) return;
//                 else if (systlist[j] == "JECUp")
//                   {
//                     Prod = "Signal_2017_JECUp";
//                   }
//                 else if (systlist[j] == "JECDown")
//                   {
//                     Prod = "Signal_2017_JECDown";
//                   }
//                 else if (systlist[j] == "JERUp")
//                   {
//                     Prod = "Signal_2017_JERUp";
//                   }
//                 else if (systlist[j] == "JERDown")
//                   {
//                     Prod = "Signal_2017_JERDown";
//                   }
//                 else
//                   {
//                     Prod = "Signal_2017";
//                   }
//             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSetM[i]+".root";
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
//                 Prod = "Signal_2017";
//                 if (systlist.size() == 0) return;
//                 else if (systlist[j] == "JECUp")
//                   {
//                     Prod = "Signal_2017_JECUp";
//                   }
//                 else if (systlist[j] == "JECDown")
//                   {
//                     Prod = "Signal_2017_JECDown";
//                   }
//                 else if (systlist[j] == "JERUp")
//                   {
//                     Prod = "Signal_2017_JERUp";
//                   }
//                 else if (systlist[j] == "JERDown")
//                   {
//                     Prod = "Signal_2017_JERDown";
//                   }
//                 else
//                   {
//                     Prod = "Signal_2017";
//                   }
//             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSetctau[i]+".root";
//             c.Reset();
//             c.Add(Path);
//             TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSetctau[i],systlist[j]);
//             float mean = 1.;//t->MeanGenWeight()
//             t->Loop(isMC, Prod, SignalSetctau[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j]);
//             delete t;
//         }
//     }


}
