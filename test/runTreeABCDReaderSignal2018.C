{
//   gROOT->ProcessLine(".L ../HistogramManager.C+");
//  gROOT->ProcessLine(".L TreeABCDReader.C+");
  std::vector<TString > systlist;
    // systlist.push_back("");
  systlist.push_back("Lumi2018Up");
  systlist.push_back("Lumi2018Down");
  systlist.push_back("L12018Up");
  systlist.push_back("L12018Down");
  systlist.push_back("Trigger2018Down");
  systlist.push_back("Trigger2018Up");
  systlist.push_back("MuonID2018Up");
  systlist.push_back("MuonID2018Down");
  systlist.push_back("MuonISO2018Up");
  systlist.push_back("MuonISO2018Down");
  systlist.push_back("PU2018Up");
  systlist.push_back("PU2018Down");
  // systlist.push_back("JEC2018Up");
  // systlist.push_back("JEC2018Down");
  // systlist.push_back("JER2018Up");
  // systlist.push_back("JER2018Down");
  systlist.push_back("Vtx2018Up");
  systlist.push_back("Vtx2018Down");
  // // systlist.push_back("RoccorUp");
  // // systlist.push_back("RoccorDown");
    systlist.push_back("XSUp");
  systlist.push_back("XSDown");
  systlist.push_back("TopPtUp");
  systlist.push_back("TopPtDown");




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
// CMS_l1_muon_prefiring_2016
// CMS_l1_muon_prefiring_2017
// CMS_l1_muon_prefiring_2018

  // systlist.push_back("PDFUp");
  // systlist.push_back("PDFDown");
  // systlist.push_back("ScaleUp");
  // systlist.push_back("ScaleDown");

  int YEAR = 2018;
  TString Year = "2018";
  TString Prod2018 = "Signal_2018_L1";
  bool isPostAPV = false;
  bool Signal = true;
  bool SameSign = false;
  bool Forward = false; 
  bool DoubleMuon = true;
  bool CorrectCorrelation = true;
  bool isMC = true;
  int mixing = 0; //-1 : Full left, 0 LR+RL , 1 Full RIght
  int Channel = 2; // 0 : EMu, 1 : SM, 2 : DM

  TChain c("ttree");


// TString SignalSetCTAU[7]={"RPV_"+Year+"_smu200","RPV_"+Year+"_smu250","RPV_"+Year+"_smu300",
// "RPV_"+Year+"_smu350","RPV_"+Year+"_smu400","RPV_"+Year+"_smu450","RPV_"+Year+"_smu500"

// };

//  for (int i = 0 ; i< 7 ; i++) 
//     {
//       TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSetCTAU[i]+".root";
//       c.Reset();
//       c.Add(Path);
//       TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSetCTAU[i],systlist);
//       float mean = 1.;//t->MeanGenWeight()
//       t->Loop(isMC, Prod2018, SignalSetCTAU[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV);
//     }


//       TString SignalSet001[34]={"RPV_"+Year+"_smu200_neu180_ctau001","RPV_"+Year+"_smu250_neu180_ctau001","RPV_"+Year+"_smu250_neu200_ctau001",
// "RPV_"+Year+"_smu250_neu230_ctau001","RPV_"+Year+"_smu300_neu180_ctau001","RPV_"+Year+"_smu300_neu200_ctau001","RPV_"+Year+"_smu300_neu250_ctau001","RPV_"+Year+"_smu300_neu280_ctau001",
// "RPV_"+Year+"_smu350_neu180_ctau001","RPV_"+Year+"_smu350_neu200_ctau001","RPV_"+Year+"_smu350_neu250_ctau001","RPV_"+Year+"_smu350_neu300_ctau001","RPV_"+Year+"_smu350_neu330_ctau001",
// "RPV_"+Year+"_smu400_neu180_ctau001","RPV_"+Year+"_smu400_neu200_ctau001","RPV_"+Year+"_smu400_neu250_ctau001","RPV_"+Year+"_smu400_neu300_ctau001","RPV_"+Year+"_smu400_neu350_ctau001",
// "RPV_"+Year+"_smu400_neu380_ctau001","RPV_"+Year+"_smu450_neu180_ctau001","RPV_"+Year+"_smu450_neu200_ctau001","RPV_"+Year+"_smu450_neu250_ctau001","RPV_"+Year+"_smu450_neu300_ctau001",
// "RPV_"+Year+"_smu450_neu350_ctau001","RPV_"+Year+"_smu450_neu400_ctau001","RPV_"+Year+"_smu450_neu430_ctau001","RPV_"+Year+"_smu500_neu180_ctau001","RPV_"+Year+"_smu500_neu200_ctau001",
// "RPV_"+Year+"_smu500_neu250_ctau001","RPV_"+Year+"_smu500_neu300_ctau001","RPV_"+Year+"_smu500_neu350_ctau001","RPV_"+Year+"_smu500_neu400_ctau001","RPV_"+Year+"_smu500_neu450_ctau001",
// "RPV_"+Year+"_smu500_neu480_ctau001"};
// //  issue  with RPV_"+Year+"_smu500_neu400_ctau001 500 450 001
//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       for (unsigned int j = 0 ; j < systlist.size(); j++)
//         {
//             Prod2018 = "Signal_2018_L1";
//             if (systlist.size() == 0) return;
//             else if (systlist[j] == "JEC2018Up")
//               {
//                 Prod2018 = "Signal_2018_JECUp";
//               }
//             else if (systlist[j] == "JEC2018Down")
//               {
//                 Prod2018 = "Signal_2018_JECDown";
//               }
//             else if (systlist[j] == "JER2018Up")
//               {
//                 Prod2018 = "Signal_2018_JERUp";
//               }
//             else if (systlist[j] == "JER2018Down")
//               {
//                 Prod2018 = "Signal_2018_JERDown";
//               }
//             else
//               {
//                 Prod2018 = "Signal_2018_L1";
//               }
//             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSet001[i]+".root";
//             c.Reset();
//             c.Add(Path);
//             TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSet001[i],systlist[j]);
//             float mean = 1.;//t->MeanGenWeight()
//             t->Loop(isMC, Prod2018, SignalSet001[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
//             delete t;
//         }
//     }



//       TString SignalSet003[34]={"RPV_"+Year+"_smu200_neu180_ctau003","RPV_"+Year+"_smu250_neu180_ctau003","RPV_"+Year+"_smu250_neu200_ctau003",
// "RPV_"+Year+"_smu250_neu230_ctau003","RPV_"+Year+"_smu300_neu180_ctau003","RPV_"+Year+"_smu300_neu200_ctau003","RPV_"+Year+"_smu300_neu250_ctau003","RPV_"+Year+"_smu300_neu280_ctau003",
// "RPV_"+Year+"_smu350_neu180_ctau003","RPV_"+Year+"_smu350_neu200_ctau003","RPV_"+Year+"_smu350_neu250_ctau003","RPV_"+Year+"_smu350_neu300_ctau003","RPV_"+Year+"_smu350_neu330_ctau003",
// "RPV_"+Year+"_smu400_neu180_ctau003","RPV_"+Year+"_smu400_neu200_ctau003","RPV_"+Year+"_smu400_neu250_ctau003","RPV_"+Year+"_smu400_neu300_ctau003","RPV_"+Year+"_smu400_neu350_ctau003",
// "RPV_"+Year+"_smu400_neu380_ctau003","RPV_"+Year+"_smu450_neu180_ctau003","RPV_"+Year+"_smu450_neu200_ctau003","RPV_"+Year+"_smu450_neu250_ctau003","RPV_"+Year+"_smu450_neu300_ctau003",
// "RPV_"+Year+"_smu450_neu350_ctau003","RPV_"+Year+"_smu450_neu400_ctau003","RPV_"+Year+"_smu450_neu430_ctau003","RPV_"+Year+"_smu500_neu180_ctau003","RPV_"+Year+"_smu500_neu200_ctau003",
// "RPV_"+Year+"_smu500_neu250_ctau003","RPV_"+Year+"_smu500_neu300_ctau003","RPV_"+Year+"_smu500_neu350_ctau003","RPV_"+Year+"_smu500_neu400_ctau003","RPV_"+Year+"_smu500_neu450_ctau003",
// "RPV_"+Year+"_smu500_neu480_ctau003"};
// //  issue RPV_"+Year+"_smu500_neu400_ctau003
//  for (int i = 0 ; i< 34 ; i++) 
//     {
//       for (unsigned int j = 0 ; j < systlist.size(); j++)
//         {
//                 Prod2018 = "Signal_2018_L1";
//                 if (systlist.size() == 0) return;
//                 else if (systlist[j] == "JEC2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JECUp";
//                   }
//                 else if (systlist[j] == "JEC2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JECDown";
//                   }
//                 else if (systlist[j] == "JER2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JERUp";
//                   }
//                 else if (systlist[j] == "JER2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JERDown";
//                   }
//                 else
//                   {
//                     Prod2018 = "Signal_2018_L1";
//                   }
//           TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSet003[i]+".root";
//           c.Reset();
//           c.Add(Path);
//           TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSet003[i],systlist[j]);
//           float mean = 1.;//t->MeanGenWeight()
//           t->Loop(isMC, Prod2018, SignalSet003[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
//           delete t;
//         }
//     }

//       TString SignalSet010[34]={"RPV_"+Year+"_smu200_neu180_ctau010","RPV_"+Year+"_smu250_neu180_ctau010","RPV_"+Year+"_smu250_neu200_ctau010",
// "RPV_"+Year+"_smu250_neu230_ctau010","RPV_"+Year+"_smu300_neu180_ctau010","RPV_"+Year+"_smu300_neu200_ctau010","RPV_"+Year+"_smu300_neu250_ctau010","RPV_"+Year+"_smu300_neu280_ctau010",
// "RPV_"+Year+"_smu350_neu180_ctau010","RPV_"+Year+"_smu350_neu200_ctau010","RPV_"+Year+"_smu350_neu250_ctau010","RPV_"+Year+"_smu350_neu300_ctau010","RPV_"+Year+"_smu350_neu330_ctau010",
// "RPV_"+Year+"_smu400_neu180_ctau010","RPV_"+Year+"_smu400_neu200_ctau010","RPV_"+Year+"_smu400_neu250_ctau010","RPV_"+Year+"_smu400_neu300_ctau010","RPV_"+Year+"_smu400_neu350_ctau010",
// "RPV_"+Year+"_smu400_neu380_ctau010","RPV_"+Year+"_smu450_neu180_ctau010","RPV_"+Year+"_smu450_neu200_ctau010","RPV_"+Year+"_smu450_neu250_ctau010","RPV_"+Year+"_smu450_neu300_ctau010",
// "RPV_"+Year+"_smu450_neu350_ctau010","RPV_"+Year+"_smu450_neu400_ctau010","RPV_"+Year+"_smu450_neu430_ctau010","RPV_"+Year+"_smu500_neu180_ctau010","RPV_"+Year+"_smu500_neu200_ctau010",
// "RPV_"+Year+"_smu500_neu250_ctau010","RPV_"+Year+"_smu500_neu300_ctau010","RPV_"+Year+"_smu500_neu350_ctau010","RPV_"+Year+"_smu500_neu400_ctau010","RPV_"+Year+"_smu500_neu450_ctau010",
// "RPV_"+Year+"_smu500_neu480_ctau010"};
// //  issue with 500 350 10 // 500 400 10
//  for (int i = 0 ; i< 34 ; i++) 
//     {
//             for (unsigned int j = 0 ; j < systlist.size(); j++)
//         {
//                 Prod2018 = "Signal_2018_L1";
//                 if (systlist.size() == 0) return;
//                 else if (systlist[j] == "JEC2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JECUp";
//                   }
//                 else if (systlist[j] == "JEC2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JECDown";
//                   }
//                 else if (systlist[j] == "JER2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JERUp";
//                   }
//                 else if (systlist[j] == "JER2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JERDown";
//                   }
//                 else
//                   {
//                     Prod2018 = "Signal_2018_L1";
//                   }
//           TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSet010[i]+".root";
//           c.Reset();
//           c.Add(Path);
//           TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSet010[i],systlist[j]);
//           float mean = 1.;//t->MeanGenWeight()
//           t->Loop(isMC, Prod2018, SignalSet010[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
//           delete t;
//         }
//     }


//       TString SignalSet030[34]={"RPV_"+Year+"_smu200_neu180_ctau030","RPV_"+Year+"_smu250_neu180_ctau030","RPV_"+Year+"_smu250_neu200_ctau030",
// "RPV_"+Year+"_smu250_neu230_ctau030","RPV_"+Year+"_smu300_neu180_ctau030","RPV_"+Year+"_smu300_neu200_ctau030","RPV_"+Year+"_smu300_neu250_ctau030","RPV_"+Year+"_smu300_neu280_ctau030",
// "RPV_"+Year+"_smu350_neu180_ctau030","RPV_"+Year+"_smu350_neu200_ctau030","RPV_"+Year+"_smu350_neu250_ctau030","RPV_"+Year+"_smu350_neu300_ctau030","RPV_"+Year+"_smu350_neu330_ctau030",
// "RPV_"+Year+"_smu400_neu180_ctau030","RPV_"+Year+"_smu400_neu200_ctau030","RPV_"+Year+"_smu400_neu250_ctau030","RPV_"+Year+"_smu400_neu300_ctau030","RPV_"+Year+"_smu400_neu350_ctau030",
// "RPV_"+Year+"_smu400_neu380_ctau030","RPV_"+Year+"_smu450_neu180_ctau030","RPV_"+Year+"_smu450_neu200_ctau030","RPV_"+Year+"_smu450_neu250_ctau030","RPV_"+Year+"_smu450_neu300_ctau030",
// "RPV_"+Year+"_smu450_neu350_ctau030","RPV_"+Year+"_smu450_neu400_ctau030","RPV_"+Year+"_smu450_neu430_ctau030","RPV_"+Year+"_smu500_neu180_ctau030","RPV_"+Year+"_smu500_neu200_ctau030",
// "RPV_"+Year+"_smu500_neu250_ctau030","RPV_"+Year+"_smu500_neu300_ctau030","RPV_"+Year+"_smu500_neu350_ctau030","RPV_"+Year+"_smu500_neu400_ctau030","RPV_"+Year+"_smu500_neu450_ctau030",
// "RPV_"+Year+"_smu500_neu480_ctau030"};
 
//  for (int i = 0 ; i< 34 ; i++) 
//     {
//             for (unsigned int j = 0 ; j < systlist.size(); j++)
//         {
//                 Prod2018 = "Signal_2018_L1";
//                 if (systlist.size() == 0) return;
//                 else if (systlist[j] == "JEC2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JECUp";
//                   }
//                 else if (systlist[j] == "JEC2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JECDown";
//                   }
//                 else if (systlist[j] == "JER2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JERUp";
//                   }
//                 else if (systlist[j] == "JER2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JERDown";
//                   }
//                 else
//                   {
//                     Prod2018 = "Signal_2018_L1";
//                   }
//           TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSet030[i]+".root";
//           c.Reset();
//           c.Add(Path);
//           TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSet030[i],systlist[j]);
//           float mean = 1.;//t->MeanGenWeight()
//           t->Loop(isMC, Prod2018, SignalSet030[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
//           delete t;
//         }
//     }

//   TString SignalSet100[34]={"RPV_"+Year+"_smu200_neu180_ctau100","RPV_"+Year+"_smu250_neu180_ctau100","RPV_"+Year+"_smu250_neu200_ctau100",
// "RPV_"+Year+"_smu250_neu230_ctau100","RPV_"+Year+"_smu300_neu180_ctau100","RPV_"+Year+"_smu300_neu200_ctau100","RPV_"+Year+"_smu300_neu250_ctau100","RPV_"+Year+"_smu300_neu280_ctau100",
// "RPV_"+Year+"_smu350_neu180_ctau100","RPV_"+Year+"_smu350_neu200_ctau100","RPV_"+Year+"_smu350_neu250_ctau100","RPV_"+Year+"_smu350_neu300_ctau100","RPV_"+Year+"_smu350_neu330_ctau100",
// "RPV_"+Year+"_smu400_neu180_ctau100","RPV_"+Year+"_smu400_neu200_ctau100","RPV_"+Year+"_smu400_neu250_ctau100","RPV_"+Year+"_smu400_neu300_ctau100","RPV_"+Year+"_smu400_neu350_ctau100",
// "RPV_"+Year+"_smu400_neu380_ctau100","RPV_"+Year+"_smu450_neu180_ctau100","RPV_"+Year+"_smu450_neu200_ctau100","RPV_"+Year+"_smu450_neu250_ctau100","RPV_"+Year+"_smu450_neu300_ctau100",
// "RPV_"+Year+"_smu450_neu350_ctau100","RPV_"+Year+"_smu450_neu400_ctau100","RPV_"+Year+"_smu450_neu430_ctau100","RPV_"+Year+"_smu500_neu180_ctau100","RPV_"+Year+"_smu500_neu200_ctau100",
// "RPV_"+Year+"_smu500_neu250_ctau100","RPV_"+Year+"_smu500_neu300_ctau100","RPV_"+Year+"_smu500_neu350_ctau100","RPV_"+Year+"_smu500_neu400_ctau100","RPV_"+Year+"_smu500_neu450_ctau100",
// "RPV_"+Year+"_smu500_neu480_ctau100"};
 
//  for (int i = 0 ; i< 34 ; i++) 
//     {
//             for (unsigned int j = 0 ; j < systlist.size(); j++)
//         {
//                 Prod2018 = "Signal_2018_L1";
//                 if (systlist.size() == 0) return;
//                 else if (systlist[j] == "JEC2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JECUp";
//                   }
//                 else if (systlist[j] == "JEC2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JECDown";
//                   }
//                 else if (systlist[j] == "JER2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JERUp";
//                   }
//                 else if (systlist[j] == "JER2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JERDown";
//                   }
//                 else
//                   {
//                     Prod2018 = "Signal_2018_L1";
//                   }
//           TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSet100[i]+".root";
//           c.Reset();
//           c.Add(Path);
//           TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSet100[i],systlist[j]);
//           float mean = 1.;//t->MeanGenWeight()
//           t->Loop(isMC, Prod2018, SignalSet100[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
//           delete t;
//         }
//     }

//       TString SignalSet300[34]={"RPV_"+Year+"_smu200_neu180_ctau300","RPV_"+Year+"_smu250_neu180_ctau300","RPV_"+Year+"_smu250_neu200_ctau300",
// "RPV_"+Year+"_smu250_neu230_ctau300","RPV_"+Year+"_smu300_neu180_ctau300","RPV_"+Year+"_smu300_neu200_ctau300","RPV_"+Year+"_smu300_neu250_ctau300","RPV_"+Year+"_smu300_neu280_ctau300",
// "RPV_"+Year+"_smu350_neu180_ctau300","RPV_"+Year+"_smu350_neu200_ctau300","RPV_"+Year+"_smu350_neu250_ctau300","RPV_"+Year+"_smu350_neu300_ctau300","RPV_"+Year+"_smu350_neu330_ctau300",
// "RPV_"+Year+"_smu400_neu180_ctau300","RPV_"+Year+"_smu400_neu200_ctau300","RPV_"+Year+"_smu400_neu250_ctau300","RPV_"+Year+"_smu400_neu300_ctau300","RPV_"+Year+"_smu400_neu350_ctau300",
// "RPV_"+Year+"_smu400_neu380_ctau300","RPV_"+Year+"_smu450_neu180_ctau300","RPV_"+Year+"_smu450_neu200_ctau300","RPV_"+Year+"_smu450_neu250_ctau300","RPV_"+Year+"_smu450_neu300_ctau300",
// "RPV_"+Year+"_smu450_neu350_ctau300","RPV_"+Year+"_smu450_neu400_ctau300","RPV_"+Year+"_smu450_neu430_ctau300","RPV_"+Year+"_smu500_neu180_ctau300","RPV_"+Year+"_smu500_neu200_ctau300",
// "RPV_"+Year+"_smu500_neu250_ctau300","RPV_"+Year+"_smu500_neu300_ctau300","RPV_"+Year+"_smu500_neu350_ctau300","RPV_"+Year+"_smu500_neu400_ctau300","RPV_"+Year+"_smu500_neu450_ctau300",
// "RPV_"+Year+"_smu500_neu480_ctau300"};
 
//  for (int i = 0 ; i< 34 ; i++) 
//     {
//             for (unsigned int j = 0 ; j < systlist.size(); j++)
//         {
//                 Prod2018 = "Signal_2018_L1";
//                 if (systlist.size() == 0) return;
//                 else if (systlist[j] == "JEC2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JECUp";
//                   }
//                 else if (systlist[j] == "JEC2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JECDown";
//                   }
//                 else if (systlist[j] == "JER2018Up")
//                   {
//                     Prod2018 = "Signal_2018_JERUp";
//                   }
//                 else if (systlist[j] == "JER2018Down")
//                   {
//                     Prod2018 = "Signal_2018_JERDown";
//                   }
//                 else
//                   {
//                     Prod2018 = "Signal_2018_L1";
//                   }
//             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSet300[i]+".root";
//             c.Reset();
//             c.Add(Path);
//             TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSet300[i],systlist[j]);
//             float mean = 1.;//t->MeanGenWeight()
//             t->Loop(isMC, Prod2018, SignalSet300[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
//             delete t;
//         }
//     }


      TString SignalSet1000[34]={"RPV_"+Year+"_smu200_neu180_ctau1000","RPV_"+Year+"_smu250_neu180_ctau1000","RPV_"+Year+"_smu250_neu200_ctau1000",
"RPV_"+Year+"_smu250_neu230_ctau1000","RPV_"+Year+"_smu300_neu180_ctau1000","RPV_"+Year+"_smu300_neu200_ctau1000","RPV_"+Year+"_smu300_neu250_ctau1000","RPV_"+Year+"_smu300_neu280_ctau1000",
"RPV_"+Year+"_smu350_neu180_ctau1000","RPV_"+Year+"_smu350_neu200_ctau1000","RPV_"+Year+"_smu350_neu250_ctau1000","RPV_"+Year+"_smu350_neu300_ctau1000","RPV_"+Year+"_smu350_neu330_ctau1000",
"RPV_"+Year+"_smu400_neu180_ctau1000","RPV_"+Year+"_smu400_neu200_ctau1000","RPV_"+Year+"_smu400_neu250_ctau1000","RPV_"+Year+"_smu400_neu300_ctau1000","RPV_"+Year+"_smu400_neu350_ctau1000",
"RPV_"+Year+"_smu400_neu380_ctau1000","RPV_"+Year+"_smu450_neu180_ctau1000","RPV_"+Year+"_smu450_neu200_ctau1000","RPV_"+Year+"_smu450_neu250_ctau1000","RPV_"+Year+"_smu450_neu300_ctau1000",
"RPV_"+Year+"_smu450_neu350_ctau1000","RPV_"+Year+"_smu450_neu400_ctau1000","RPV_"+Year+"_smu450_neu430_ctau1000","RPV_"+Year+"_smu500_neu180_ctau1000","RPV_"+Year+"_smu500_neu200_ctau1000",
"RPV_"+Year+"_smu500_neu250_ctau1000","RPV_"+Year+"_smu500_neu300_ctau1000","RPV_"+Year+"_smu500_neu350_ctau1000","RPV_"+Year+"_smu500_neu400_ctau1000","RPV_"+Year+"_smu500_neu450_ctau1000",
"RPV_"+Year+"_smu500_neu480_ctau1000"};
 
 for (int i = 19 ; i< 20 ; i++) 
    {
            for (unsigned int j = 0 ; j < systlist.size(); j++)
        {
                Prod2018 = "Signal_2018_L1";
                if (systlist.size() == 0) return;
                else if (systlist[j] == "JEC2018Up")
                  {
                    Prod2018 = "Signal_2018_JECUp";
                  }
                else if (systlist[j] == "JEC2018Down")
                  {
                    Prod2018 = "Signal_2018_JECDown";
                  }
                else if (systlist[j] == "JER2018Up")
                  {
                    Prod2018 = "Signal_2018_JERUp";
                  }
                else if (systlist[j] == "JER2018Down")
                  {
                    Prod2018 = "Signal_2018_JERDown";
                  }
                else
                  {
                    Prod2018 = "Signal_2018_L1";
                  }
            TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSet1000[i]+".root";
            c.Reset();
            c.Add(Path);
            TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSet1000[i],systlist[j]);
            float mean = 1.;//t->MeanGenWeight()
            t->Loop(isMC, Prod2018, SignalSet1000[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
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
// //                 Prod2018 = "Signal_2018_L1";
// //                 if (systlist.size() == 0) return;
// //                 else if (systlist[j] == "JECUp")
// //                   {
// //                     Prod2018 = "Signal_2018_JECUp";
// //                   }
// //                 else if (systlist[j] == "JECDown")
// //                   {
// //                     Prod2018 = "Signal_2018_JECDown";
// //                   }
// //                 else if (systlist[j] == "JERUp")
// //                   {
// //                     Prod2018 = "Signal_2018_JERUp";
// //                   }
// //                 else if (systlist[j] == "JERDown")
// //                   {
// //                     Prod2018 = "Signal_2018_JERDown";
// //                   }
// //                 else
// //                   {
// //                     Prod2018 = "Signal_2018_L1";
// //                   }
// //             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSetM[i]+".root";
// //             c.Reset();
// //             c.Add(Path);
// //             TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSetM[i],systlist[j]);
// //             float mean = 1.;//t->MeanGenWeight()
// //             t->Loop(isMC, Prod2018, SignalSetM[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j],CorrectCorrelation);
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
// //                 Prod2018 = "Signal_2018_L1";
// //                 if (systlist.size() == 0) return;
// //                 else if (systlist[j] == "JECUp")
// //                   {
// //                     Prod2018 = "Signal_2018_JECUp";
// //                   }
// //                 else if (systlist[j] == "JECDown")
// //                   {
// //                     Prod2018 = "Signal_2018_JECDown";
// //                   }
// //                 else if (systlist[j] == "JERUp")
// //                   {
// //                     Prod2018 = "Signal_2018_JERUp";
// //                   }
// //                 else if (systlist[j] == "JERDown")
// //                   {
// //                     Prod2018 = "Signal_2018_JERDown";
// //                   }
// //                 else
// //                   {
// //                     Prod2018 = "Signal_2018_L1";
// //                   }
// //             TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod2018+"/Mini"+SignalSetctau[i]+".root";
// //             c.Reset();
// //             c.Add(Path);
// //             TreeABCDReader* t = new TreeABCDReader(&c,Prod2018, SignalSetctau[i],systlist[j]);
// //             float mean = 1.;//t->MeanGenWeight()
// //             t->Loop(isMC, Prod2018, SignalSetctau[i],Signal, SameSign, Forward, YEAR, mean, DoubleMuon, mixing, Channel, isPostAPV,systlist[j]);
// //             delete t;
// //         }
// //     }


}
