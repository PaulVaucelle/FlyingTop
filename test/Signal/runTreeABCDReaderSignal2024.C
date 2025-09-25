{
//   gROOT->ProcessLine(".L ../HistogramManager.C+");
//  gROOT->ProcessLine(".L TreeABCDReader.C+");
  std::vector<TString > systlist;
  //   systlist.push_back("");
  // systlist.push_back("Lumi2024Up");
  // systlist.push_back("Lumi2024Down");
  // systlist.push_back("L12024Up");
  // systlist.push_back("L12024Down");
  // systlist.push_back("Trigger2024Down");
  // systlist.push_back("Trigger2024Up");
  // systlist.push_back("MuonID2024Up");
  // systlist.push_back("MuonID2024Down");
  // systlist.push_back("MuonISO2024Up");
  // systlist.push_back("MuonISO2024Down");
  // systlist.push_back("PU2024Up");
  // systlist.push_back("PU2024Down");
  systlist.push_back("JEC2024Up");
  systlist.push_back("JEC2024Down");
  systlist.push_back("JER2024Up");
  systlist.push_back("JER2024Down");
  // systlist.push_back("Vtx2024Up");
  // systlist.push_back("Vtx2024Down");
  // // // systlist.push_back("RoccorUp");
  // // // systlist.push_back("RoccorDown");
  //   systlist.push_back("XSUp");
  // systlist.push_back("XSDown");
  // systlist.push_back("TopPtUp");
  // systlist.push_back("TopPtDown");




// lumi_2016
// lumi_2017
// lumi_2018
// lumi_13p6TeV_2022
// lumi_13p6TeV_2023
// lumi_13p6TeV_2024
// CMS_pileup_(2016preVFP|2016postVFP|2017|2018|2022|2022EE|2023|2023BPix)
// CMS_eff_e_id_(2016|2017|2018|2022|2022EE|2023|2023BPix)
// CMS_eff_e_reco_(2016|2017|2018|2022|2022EE|2023|2023BPix)
// CMS_scale_m_(2016|2017|2018|2022|2022EE|2023|2023BPix)

// CMS_eff_m_id_stat_(2016preVFP|2016postVFP|2017|2018|2022|2022EE|2023|2023BPix)
// CMS_eff_m_iso_stat_(2016preVFP|2016postVFP|2017|2018|2022|2022EE|2023|2023BPix)
// CMS_eff_m_trigger_(2016|2017|2018|2022|2022EE|2023|2023BPix)
// CMS_scale_j_(2016|2016preVFP|2016postVFP|2017|2018|2022|2022EE|2023|2023BPix)
//  CMS_res_j_(2016|2016preVFP|2016postVFP|2017|2018|2022|2022EE|2023|2023BPix)
//  top_pt_reweighting
// cross_section_
// CMS_muon_prefiring_2016
// CMS_muon_prefiring_2017
// CMS_muon_prefiring_2018

  // systlist.push_back("PDFUp");
  // systlist.push_back("PDFDown");
  // systlist.push_back("ScaleUp");
  // systlist.push_back("ScaleDown");

  int YEAR = 2024;
  TString Year = "2024";
  TString Prod = "RPV_2024"; //DATAMC2018_EMU_10_06_2024 //  EMU: SYST_EMU_CTAU100/Prod100_EMU_JER_up !! M: SYST_CTAU100/Prod100_JER_up // BKG MC and Data PROD_ANNIVERSAIRE_2024

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


      TString SignalSet001[34]={
         "RPV_2024_Par-ct-001-MChi-180-MSmu-200","RPV_2024_Par-ct-001-MChi-180-MSmu-250","RPV_2024_Par-ct-001-MChi-200-MSmu-250","RPV_2024_Par-ct-001-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-300","RPV_2024_Par-ct-001-MChi-200-MSmu-300","RPV_2024_Par-ct-001-MChi-250-MSmu-300","RPV_2024_Par-ct-001-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-350","RPV_2024_Par-ct-001-MChi-200-MSmu-350","RPV_2024_Par-ct-001-MChi-250-MSmu-350","RPV_2024_Par-ct-001-MChi-300-MSmu-350","RPV_2024_Par-ct-001-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-400","RPV_2024_Par-ct-001-MChi-200-MSmu-400","RPV_2024_Par-ct-001-MChi-250-MSmu-400","RPV_2024_Par-ct-001-MChi-300-MSmu-400","RPV_2024_Par-ct-001-MChi-350-MSmu-400","RPV_2024_Par-ct-001-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-450","RPV_2024_Par-ct-001-MChi-200-MSmu-450","RPV_2024_Par-ct-001-MChi-250-MSmu-450","RPV_2024_Par-ct-001-MChi-300-MSmu-450","RPV_2024_Par-ct-001-MChi-350-MSmu-450","RPV_2024_Par-ct-001-MChi-400-MSmu-450","RPV_2024_Par-ct-001-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-001-MChi-180-MSmu-500","RPV_2024_Par-ct-001-MChi-200-MSmu-500","RPV_2024_Par-ct-001-MChi-250-MSmu-500","RPV_2024_Par-ct-001-MChi-300-MSmu-500","RPV_2024_Par-ct-001-MChi-350-MSmu-500","RPV_2024_Par-ct-001-MChi-400-MSmu-500","RPV_2024_Par-ct-001-MChi-450-MSmu-500","RPV_2024_Par-ct-001-MChi-480-MSmu-500"
}; 
//  issue  with RPV_"+Year+"_smu500_neu400_ctau001 500 450 001
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j <systlist.size(); j++)//systlist.size()
        {
            Prod = "RPV_2024";
            if (systlist.size() == 0) return;
            else if (systlist[j] == "JEC"+Year+"Up")
              {
                Prod = "RPV_2024_JECUp";
              }
            else if (systlist[j] == "JEC"+Year+"Down")
              {
                Prod = "RPV_2024_JECDown";
              }
            else if (systlist[j] == "JER"+Year+"Up")
              {
                Prod = "RPV_2024_JERUp";
              }
            else if (systlist[j] == "JER"+Year+"Down")
              {
                Prod = "RPV_2024_JERDown";
              }
            else
              {
                Prod = "RPV_2024";
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



      TString SignalSet003[34]={
         "RPV_2024_Par-ct-003-MChi-180-MSmu-200","RPV_2024_Par-ct-003-MChi-180-MSmu-250","RPV_2024_Par-ct-003-MChi-200-MSmu-250","RPV_2024_Par-ct-003-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-300","RPV_2024_Par-ct-003-MChi-200-MSmu-300","RPV_2024_Par-ct-003-MChi-250-MSmu-300","RPV_2024_Par-ct-003-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-350","RPV_2024_Par-ct-003-MChi-200-MSmu-350","RPV_2024_Par-ct-003-MChi-250-MSmu-350","RPV_2024_Par-ct-003-MChi-300-MSmu-350","RPV_2024_Par-ct-003-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-400","RPV_2024_Par-ct-003-MChi-200-MSmu-400","RPV_2024_Par-ct-003-MChi-250-MSmu-400","RPV_2024_Par-ct-003-MChi-300-MSmu-400","RPV_2024_Par-ct-003-MChi-350-MSmu-400","RPV_2024_Par-ct-003-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-450","RPV_2024_Par-ct-003-MChi-200-MSmu-450","RPV_2024_Par-ct-003-MChi-250-MSmu-450","RPV_2024_Par-ct-003-MChi-300-MSmu-450","RPV_2024_Par-ct-003-MChi-350-MSmu-450","RPV_2024_Par-ct-003-MChi-400-MSmu-450","RPV_2024_Par-ct-003-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-003-MChi-180-MSmu-500","RPV_2024_Par-ct-003-MChi-200-MSmu-500","RPV_2024_Par-ct-003-MChi-250-MSmu-500","RPV_2024_Par-ct-003-MChi-300-MSmu-500","RPV_2024_Par-ct-003-MChi-350-MSmu-500","RPV_2024_Par-ct-003-MChi-400-MSmu-500","RPV_2024_Par-ct-003-MChi-450-MSmu-500","RPV_2024_Par-ct-003-MChi-480-MSmu-500"

};
//  issue RPV_"+Year+"_smu500_neu400_ctau003
 for (int i = 0 ; i< 34 ; i++) 
    {
      for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2024";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2024_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2024_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2024_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2024_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2024";
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


      TString SignalSet010[34]={
         "RPV_2024_Par-ct-010-MChi-180-MSmu-200","RPV_2024_Par-ct-010-MChi-180-MSmu-250","RPV_2024_Par-ct-010-MChi-200-MSmu-250","RPV_2024_Par-ct-010-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-300","RPV_2024_Par-ct-010-MChi-200-MSmu-300","RPV_2024_Par-ct-010-MChi-250-MSmu-300","RPV_2024_Par-ct-010-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-350","RPV_2024_Par-ct-010-MChi-200-MSmu-350","RPV_2024_Par-ct-010-MChi-250-MSmu-350","RPV_2024_Par-ct-010-MChi-300-MSmu-350","RPV_2024_Par-ct-010-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-400","RPV_2024_Par-ct-010-MChi-200-MSmu-400","RPV_2024_Par-ct-010-MChi-250-MSmu-400","RPV_2024_Par-ct-010-MChi-300-MSmu-400","RPV_2024_Par-ct-010-MChi-350-MSmu-400","RPV_2024_Par-ct-010-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-450","RPV_2024_Par-ct-010-MChi-200-MSmu-450","RPV_2024_Par-ct-010-MChi-250-MSmu-450","RPV_2024_Par-ct-010-MChi-300-MSmu-450","RPV_2024_Par-ct-010-MChi-350-MSmu-450","RPV_2024_Par-ct-010-MChi-400-MSmu-450","RPV_2024_Par-ct-010-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-010-MChi-180-MSmu-500","RPV_2024_Par-ct-010-MChi-200-MSmu-500","RPV_2024_Par-ct-010-MChi-250-MSmu-500","RPV_2024_Par-ct-010-MChi-300-MSmu-500","RPV_2024_Par-ct-010-MChi-350-MSmu-500","RPV_2024_Par-ct-010-MChi-400-MSmu-500","RPV_2024_Par-ct-010-MChi-450-MSmu-500","RPV_2024_Par-ct-010-MChi-480-MSmu-500"

};
//  issue with 500 350 10 // 500 400 10
 for (int i = 0 ; i< 34 ; i++) 
    {
            for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2024";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2024_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2024_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2024_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2024_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2024";
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
            for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2024";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2024_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2024_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2024_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2024_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2024";
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

    TString SignalSet100[34]={
         "RPV_2024_Par-ct-100-MChi-180-MSmu-200","RPV_2024_Par-ct-100-MChi-180-MSmu-250","RPV_2024_Par-ct-100-MChi-200-MSmu-250","RPV_2024_Par-ct-100-MChi-230-MSmu-250",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-300","RPV_2024_Par-ct-100-MChi-200-MSmu-300","RPV_2024_Par-ct-100-MChi-250-MSmu-300","RPV_2024_Par-ct-100-MChi-280-MSmu-300",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-350","RPV_2024_Par-ct-100-MChi-200-MSmu-350","RPV_2024_Par-ct-100-MChi-250-MSmu-350","RPV_2024_Par-ct-100-MChi-300-MSmu-350","RPV_2024_Par-ct-100-MChi-330-MSmu-350",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-400","RPV_2024_Par-ct-100-MChi-200-MSmu-400","RPV_2024_Par-ct-100-MChi-250-MSmu-400","RPV_2024_Par-ct-100-MChi-300-MSmu-400","RPV_2024_Par-ct-100-MChi-350-MSmu-400","RPV_2024_Par-ct-100-MChi-380-MSmu-400",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-450","RPV_2024_Par-ct-100-MChi-200-MSmu-450","RPV_2024_Par-ct-100-MChi-250-MSmu-450","RPV_2024_Par-ct-100-MChi-300-MSmu-450","RPV_2024_Par-ct-100-MChi-350-MSmu-450","RPV_2024_Par-ct-100-MChi-400-MSmu-450","RPV_2024_Par-ct-100-MChi-430-MSmu-450",
         "RPV_2024_Par-ct-100-MChi-180-MSmu-500","RPV_2024_Par-ct-100-MChi-200-MSmu-500","RPV_2024_Par-ct-100-MChi-250-MSmu-500","RPV_2024_Par-ct-100-MChi-300-MSmu-500","RPV_2024_Par-ct-100-MChi-350-MSmu-500","RPV_2024_Par-ct-100-MChi-400-MSmu-500","RPV_2024_Par-ct-100-MChi-450-MSmu-500","RPV_2024_Par-ct-100-MChi-480-MSmu-500"

};
 
 for (int i = 0 ; i< 34 ; i++) 
    {
            for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2024";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2024_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2024_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2024_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2024_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2024";
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
            for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2024";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2024_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2024_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2024_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2024_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2024";
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
            for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod = "RPV_2024";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC"+Year+"Up")
                  {
                    Prod = "RPV_2024_JECUp";
                  }
                else if (systlist[j] == "JEC"+Year+"Down")
                  {
                    Prod = "RPV_2024_JECDown";
                  }
                else if (systlist[j] == "JER"+Year+"Up")
                  {
                    Prod = "RPV_2024_JERUp";
                  }
                else if (systlist[j] == "JER"+Year+"Down")
                  {
                    Prod = "RPV_2024_JERDown";
                  }
                else
                  {
                    Prod = "RPV_2024";
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



// //       TString SignalSetM[34]={"RPV_"+Year+"_smu200_neu180","RPV_"+Year+"_smu250_neu180","RPV_"+Year+"_smu250_neu200",
// // "RPV_"+Year+"_smu250_neu230","RPV_"+Year+"_smu300_neu180","RPV_"+Year+"_smu300_neu200","RPV_"+Year+"_smu300_neu250","RPV_"+Year+"_smu300_neu280",
// // "RPV_"+Year+"_smu350_neu180","RPV_"+Year+"_smu350_neu200","RPV_"+Year+"_smu350_neu250","RPV_"+Year+"_smu350_neu300","RPV_"+Year+"_smu350_neu330",
// // "RPV_"+Year+"_smu400_neu180","RPV_"+Year+"_smu400_neu200","RPV_"+Year+"_smu400_neu250","RPV_"+Year+"_smu400_neu300","RPV_"+Year+"_smu400_neu350",
// // "RPV_"+Year+"_smu400_neu380","RPV_"+Year+"_smu450_neu180","RPV_"+Year+"_smu450_neu200","RPV_"+Year+"_smu450_neu250","RPV_"+Year+"_smu450_neu300",
// // "RPV_"+Year+"_smu450_neu350","RPV_"+Year+"_smu450_neu400","RPV_"+Year+"_smu450_neu430","RPV_"+Year+"_smu500_neu180","RPV_"+Year+"_smu500_neu200",
// // "RPV_"+Year+"_smu500_neu250","RPV_"+Year+"_smu500_neu300","RPV_"+Year+"_smu500_neu350","RPV_"+Year+"_smu500_neu400","RPV_"+Year+"_smu500_neu450",
// // "RPV_"+Year+"_smu500_neu480"};
 
// //  for (int i = 0 ; i< 34 ; i++) 
// //     {
// //             for (unsigned int j = 0 ; j < systlist.size(); j++)
// //         {
// //                 Prod = "RPV_2024";
// //                 if (systlist.size() == 0) return;
// //                 else if (systlist[j] == "JECUp")
// //                   {
// //                     Prod = "RPV_2024_JECUp";
// //                   }
// //                 else if (systlist[j] == "JECDown")
// //                   {
// //                     Prod = "RPV_2024_JECDown";
// //                   }
// //                 else if (systlist[j] == "JERUp")
// //                   {
// //                     Prod = "RPV_2024_JERUp";
// //                   }
// //                 else if (systlist[j] == "JERDown")
// //                   {
// //                     Prod = "RPV_2024_JERDown";
// //                   }
// //                 else
// //                   {
// //                     Prod = "RPV_2024";
// //                   }
// //             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSetM[i]+".root";
// //             c.Reset();
// //             c.Add(Path);
// //             TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSetM[i],systlist[j]);
// //             float mean = 1.;//t->MeanGenWeight()
// //             t->Loop(isMC, Prod, SignalSetM[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
// //             delete t;
// //         }
// //     }




// //       TString SignalSetctau[7]={"RPV_ctau001","RPV_ctau003","RPV_ctau010",
// // "RPV_ctau030","RPV_ctau100","RPV_ctau300","RPV_ctau1000"
// // };
 
// //  for (int i = 0 ; i< 7 ; i++) 
// //     {
// //             for (unsigned int j = 0 ; j < systlist.size(); j++)
// //         {
// //                 Prod = "RPV_2024";
// //                 if (systlist.size() == 0) return;
// //                 else if (systlist[j] == "JECUp")
// //                   {
// //                     Prod = "RPV_2024_JECUp";
// //                   }
// //                 else if (systlist[j] == "JECDown")
// //                   {
// //                     Prod = "RPV_2024_JECDown";
// //                   }
// //                 else if (systlist[j] == "JERUp")
// //                   {
// //                     Prod = "RPV_2024_JERUp";
// //                   }
// //                 else if (systlist[j] == "JERDown")
// //                   {
// //                     Prod = "RPV_2024_JERDown";
// //                   }
// //                 else
// //                   {
// //                     Prod = "RPV_2024";
// //                   }
// //             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+SignalSetctau[i]+".root";
// //             c.Reset();
// //             c.Add(Path);
// //             TreeABCDReader* t = new TreeABCDReader(&c,Prod, SignalSetctau[i],systlist[j]);
// //             float mean = 1.;//t->MeanGenWeight()
// //             t->Loop(isMC, Prod, SignalSetctau[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j]);
// //             delete t;
// //         }
// //     }


}
