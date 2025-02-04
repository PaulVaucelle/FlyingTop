#define TreeABCDReader_cxx
#include "TreeABCDReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime> 
#include "TTree.h"
#include "../HistogramManager.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1F.h"
#include <TDirectory.h>
#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/MCWeights.h"

double TreeABCDReader::MeanGenWeight(TString thesample, TString Prod)
   {
      double mean = 0;
      double ratio = 1.;
      // std::cout<<"number of entries : "<<fChain->GetEntries()<<std::endl;
      std::ofstream ofs ("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Efficacity_"+thesample+".txt", std::ofstream::out);
      if (fChain->GetEntries()==0) {return 1.;}

      TString hNorma = "hEvents";
      TFile* f1_DY= new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+thesample+".root");
      f1_DY->cd("");
      double norm = 0.;
      TDirectory* dir = f1_DY->GetDirectory("FlyingTop");
      if (dir) {
         dir->cd();  // Change to this directory if needed
         // Now you can access histograms, trees, etc. from the directory
            TH1D*  e1_DY = (TH1D*)gROOT->FindObject(hNorma);
      e1_DY->Sumw2();
         if  ( e1_DY->GetEntries() > 0 ) norm =  e1_DY->GetEntries();
      std::cout<<"norm : "<<norm<<std::endl;
      }
      int nentry = norm;
      for (int i = 0; i < nentry; i++)
      {
         fChain->GetEntry(i);
         mean += abs(minitree_only_gen_wt->at(0));
      }
      ratio = mean/norm;
      ofs<<thesample<<" with mean weight : "<<ratio;
      ofs.close();
      std::cout<<thesample<<" with other mean weight : "<<mean/fChain->GetEntries()<<std::endl;
      return ratio;
   }

void TreeABCDReader::Loop(bool isMC, TString Prod, TString sample, bool Signal, bool SS, bool FWD, int Year, float mean, bool DoubleMuon, int mixing, int Channel,  bool isPostAPV)
{

  TString thesample  = sample;

   TString ADD_Text = "OS_2p4";
   if (SS)
      {
         ADD_Text = "SS_2p4";
      }
   else if (FWD)
      {
         ADD_Text = "OS_3p0";
      }
   else if (SS && FWD)
      {
         ADD_Text = "SS_3p0";
      }
   TString CHANNEL = "DM";
   if (Channel == 0)
      {
         CHANNEL = "EM";
      }
   else if (Channel == 1)
      {
         CHANNEL = "SM";
      }
double norm = 0.;
   if (!Signal)
      {
         TString hNorma = "hEvents";
         TFile* f1_DY= new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+sample+".root");
         f1_DY->cd("");
         
         TDirectory* dir = f1_DY->GetDirectory("FlyingTop");
         if (dir) {
            dir->cd();  // Change to this directory if needed
            // Now you can access histograms, trees, etc. from the directory
               TH1D*  e1_DY = (TH1D*)gROOT->FindObject(hNorma);
         e1_DY->Sumw2();
            if  ( e1_DY->GetEntries() > 0 ) norm =  e1_DY->GetEntries();
         std::cout<<"norm : "<<norm<<std::endl;
         }
      }

   bool firstinit = false;
   TString samplename = thesample;
   for(unsigned int i=0; i< systlist.size(); i++){
     
     if( systlist[i]== "") samplename = thesample;
     else                  samplename = thesample+"_"+systlist[i];
     
     
      firstinit = true;

   TFile * theoutputfile = new TFile( ("../"+Prod+"/histofile_"+CHANNEL+"_"+ADD_Text+"_"+samplename+".root").Data() , "recreate");
   theoutputfile->cd();
   initializeHisto(samplename, firstinit);
   if (fChain == 0) return;





   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries2 = fChain->GetEntries();
   std::cout<<"number of entries : "<<nentries<<std::endl;
   std::cout<<"number of entries2 : "<<nentries2<<std::endl;
   // TFile* FiltreTest_File =  new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+thesample+".root");
   // FiltreTest_File->cd();
   // TH1F *Filtretest = (TH1F*)(FiltreTest_File->Get("hData_Filter"));//ok
   // int GoodNentries = Filtretest->GetEntries();
   // // std::cout<<"test nentries : "<<GoodNentries<<std::endl;
   
// Normalisation factor (XS)

  float XS = 1.;
  float XS_up = 1.;
   float XS_down = 1.;
if (thesample.Contains("DYJetsToLL_M-10to50"))                     { XS = 22635;   }//15910
if (thesample.Contains("DYJetsToLL_M-50"))                         { XS = 6225.4;      }//5379.0;   }
if (thesample.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays")) { XS = 21.6;   }//ok
if (thesample.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     { XS = 21.6;   }//ok
if (thesample.Contains("TTJets_DiLept"))                          { XS = 53.07;     }
if (thesample.Contains("WWTo2L2Nu_MLL_200To600"))                 { XS = 11.09;     }//ok
if (thesample.Contains("WWTo2L2Nu_MLL_600To1200"))                { XS = 11.09;     }//ok
if (thesample.Contains("WWTo2L2Nu"))                              { XS = 11.09;     }//ok
if (thesample.Contains("WZTo2Q2L_mllmin4p0"))                     { XS = 6.535;     }//ok
if (thesample.Contains("ZZTo2Q2L_mllmin4p0"))                     { XS = 3.676;     }//ok
if (thesample.Contains("TTTo2L2Nu"))                              { XS = 88.5;      }//ok
if (thesample.Contains("ttWJetsToLNu_5f_EWK"))                    { XS = 0.290;     }//ok // not found on XSDB, no file on tier2...approximation
      //Took 0.868 pb (CMS-TOP-21-011)
     // as a starting point and then divided by 3 (lepton universality)
if (thesample.Contains("TTZToLL_5f"))                             { XS = 0.05188;   }//ok//Not found on XSDB => used ana.py macro : 5.188e-02 +- 2.437e-04 pb
if (thesample.Contains("TTToHadronic"))                           { XS = 377.9;     }//ok
if (thesample.Contains("TTWW"))                                   { XS = 0.006992;  }//found on XSDB
if (thesample.Contains("TTToSemiLeptonic")  )                      { XS = 365.34;    }//ok


// float ReScaleXS = 1.0;
//  float xsec_DY[2], xsec_TT[3], xsec_ST[2], xsec_TTV[3], xsec_VV[3]; 
//  xsec_DY[0] = ReScaleXS*18610; xsec_DY[1] = ReScaleXS*6077.0;
//  xsec_TT[0] = ReScaleXS*88.3; xsec_TT[1] = ReScaleXS*365.3; xsec_TT[2] = ReScaleXS*378.0;
//  xsec_ST[0] = ReScaleXS*32.51; xsec_ST[1] = ReScaleXS*32.45;
//  xsec_TTV[0] = ReScaleXS*0.29; xsec_TTV[1] = ReScaleXS*0.052; xsec_TTV[2] =  ReScaleXS*0.0070;
//  xsec_VV[0] = ReScaleXS*11.09; xsec_VV[1] = ReScaleXS*6.57; xsec_VV[2] = ReScaleXS*3.68;

// Signla XS are given in fb
PDFWeight PDFW(mixing);

if (thesample.Contains("smu200"))
    {        
      XS = 0.001*PDFW.GetXS(200);// the parmaeter is the mass of the smuon
      XS_up = 0.001*PDFW.GetXSUp(200);
      XS_down = 0.001*PDFW.GetXSDown(200);
    }
    if (thesample.Contains("smu250"))
    {        
      XS = 0.001*PDFW.GetXS(250);// the parmaeter is the mass of the smuon
      XS_up = 0.001*PDFW.GetXSUp(250);
      XS_down = 0.001*PDFW.GetXSDown(250);
    }
    if (thesample.Contains("smu300"))
    {        
      XS = 0.001*PDFW.GetXS(300);// the parmaeter is the mass of the smuon
      XS_up = 0.001*PDFW.GetXSUp(300);
      XS_down = 0.001*PDFW.GetXSDown(300);
    }

    if (thesample.Contains("smu350"))
    {        
      XS = 0.001*PDFW.GetXS(350);// the parmaeter is the mass of the smuon
      XS_up = 0.001*PDFW.GetXSUp(350);
      XS_down = 0.001*PDFW.GetXSDown(350);
    }

    if (thesample.Contains("smu400"))
    {        
      XS = 0.001*PDFW.GetXS(400);// the parmaeter is the mass of the smuon
      XS_up = 0.001*PDFW.GetXSUp(400);
      XS_down = 0.001*PDFW.GetXSDown(400);
    }

    if (thesample.Contains("smu450"))
    {        
      XS = 0.001*PDFW.GetXS(450);// the parmaeter is the mass of the smuon
      XS_up = 0.001*PDFW.GetXSUp(450);
      XS_down = 0.001*PDFW.GetXSDown(450);
    }

    if (thesample.Contains("smu500"))
    {        
      XS = 0.001*PDFW.GetXS(500);// the parmaeter is the mass of the smuon
      XS_up = 0.001*PDFW.GetXSUp(500);
      XS_down = 0.001*PDFW.GetXSDown(500);
    }


   float NormFactor1 = 1;

   double NormFactor=1;
   double top_pt_wt=1;
   double Pref_PU_gen_wt=1;

   double NormFactor_mu=1;
   double NormFactor_mu2 =1;
   double NormFactor_ele=1;

   float NormFactorLumiUp = 1.;
   float NormFactorLumi = 1.;
   float NormFactorLumiDown = 1.;   

   float NormFactorXSUp = 1.;
   float NormFactorXS = 1.;
   float NormFactorXSDown = 1.;

   float Lumi = 59700.;
   float LumiUp = 66000.0;
   float LumiDown = 55000.0;
   LumiWeights LUMI(Year);

   // the year is taken into account at the creation of the instance of LumiWeights
   // Uncertainties are taken year by year (and not for the full run 2, see LUMIPOG)
   Lumi = LUMI.GetLumi();
   LumiUp = LUMI.GetLumiUp();
   LumiDown = LUMI.GetLumiDown();

  if (!isMC)//<=> Data
    {
      NormFactorLumi = 1;
      NormFactorLumiUp = 1.;  
      NormFactorLumiDown = 1.;    

      NormFactorXS = 1;  
      NormFactorXSUp = 1.;  
      NormFactorXSDown = 1.;                                                                                                                                                               
    }
  else{
      NormFactorLumi =  XS*Lumi;  
      NormFactorLumiUp = XS*LumiUp;
      NormFactorLumiDown = XS*LumiDown;

      NormFactorXS =  XS*Lumi;  
      NormFactorXSUp = XS_up*Lumi;
      NormFactorXSDown = XS_down*Lumi;
   }

  // !! -- Parameters -- !! //
   // // std::cout<<"Mixing : "<<mixing<<std::endl;
   // // std::cout<<"Msmuon : "<<MSmuon<<std::endl;
   bool BlindSR = false;
   float Nevent = 0;
   bool signal = Signal;
   int allevents = 0;
   double lowHpt = 20.;
   double LT = 0.;
   // nentries = 1;

   double nFilterEvt = 0;
   double nFilterJet = 0;
   double nFilternHemi = 0;
   double nFilterHpt = 0;

   double nRecoVertex = 0;
   double nReco1Vertex = 0;
   double nReco1TightVertex = 0;
   double nReco2Vertex = 0;
   double nReco2TightVertex = 0;
   double nRecoTightVertex = 0;
   double nRecoLooseVertex = 0;

   bool showoutput = false;
  // !! -- END OF Parameters -- !! //
  

     
      //---------------------------------------------------------------------------//
      // !! ----------------------------- Scale Factors---------------------------------//
      //---------------------------------------------------------------------------//

     // Where L1 is always a muon and L2 is either a muon for the dimuon channel or a
   // an electron for the Emu channel
   TFile*            fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
   TFile*            fL1_ID_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root");
   TFile*           fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root");
   TFile*           fL2_ISO_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root");
   // $$$$$ Changer les noms des fichiers
   TFile*           fL2_Reco_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root");
   TFile*            fL2_Reco_SF2= new TFile("../../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2018.root");
   TFile*            fL2_ID_SF= new TFile("../../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_EGM2D_U18.root");
   TFile*            fL1L2_TRG_SF= new TFile("../../Scale_factors/Electron/Top_trigger_group/TriggerSF_2018_ULv2.root");
   TFile*            fL1_TRG_SF = new TFile("../../Scale_factors/Muon/DoubleMuon_2018/NUM_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_DEN_MiniIsoTight_and_TightID_abseta_pt.root");
   
   
   TFile*            fL1L2_TRG_SFerr = new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Scale_factors/Trig_2018/trigSF_2Derr.root");



   TFile*            fEle_SF = new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/DATAMCREADER/EMU/hRatioCor.root");
   TH1F*             hEle_SF = (TH1F*)(fEle_SF->Get("hRatio"));

   TFile*            fEle_2DSF = new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/DATAMCREADER/EMU/2DSFele.root");
   TH2F*             hEle_2DSF = (TH2F*)(fEle_2DSF->Get("h_ratio"));




   TH2F *fh2DL1_ID_SF1;//
   TH2F *fh2DL2_ID_SF1;
   TH2F *fh2DL1_ID_SF1err;//
   TH2F *fh2DL2_ID_SF1err;


   TH2F *fh2DL1_ISO_SF1;
   TH2F *fh2DL2_ISO_SF1;
   TH2F *fh2DL1_ISO_SF1err;
   TH2F *fh2DL2_ISO_SF1err;

   TH2F *fh2DL1_Reco_SF1;
   TH2F *fh2DL2_Reco_SF1;

   
   TH2F *fhL1L2_TRG_SF; // x- axiss pt and y axis pt  
   TH2F *fhL1L2_TRG_SFerr; // x- axiss pt and y axis pt 


   int MuMu = Channel; // Emu = 0 ; SingleMuon = 1; DiMuon = 2
   
   if (MuMu == 0 && isMC)
      {
         if (Year == 2018)
            {
               //cout<< " 2018 year="<<endl;
               fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root");
               fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root");
               
               fL2_Reco_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root");
               fL2_ISO_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2018.root");
               fL2_ID_SF= new TFile("../../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_EGM2D_U18.root");

               // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2018_ULv2.root");
               fL1L2_TRG_SF= new TFile("../Scale_factors/Trig_2018/trigSF_2D.root");
               fL1L2_TRG_SFerr = new TFile("../Scale_factors/Trig_2018/trigSF_2Derr.root");
            }
         if(Year == 2017)
            {
               //cout<<" year 2017"<<endl;
               fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2017_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root");
               fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root");
               
               fL2_Reco_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2017.root");
               fL2_ISO_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2017.root");
               fL2_ID_SF= new TFile("../../Scale_factors/Electron/egammaEffi.txt_EGM2D_Tight_UL17.root");

               // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2018_ULv2.root");
               fL1L2_TRG_SF= new TFile("../Scale_factors/Trig_2018/trigSF_2D.root");
               fL1L2_TRG_SFerr = new TFile("../Scale_factors/Trig_2017/trigSF_2Derr.root");
            }
         if(Year == 2016)
            {
               if (isPostAPV){
                  //cout<<" post 2016="<<endl;
                  fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root");
                  fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root");

                  fL2_Reco_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root");
                  fL2_ISO_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root");
                  fL2_ID_SF= new TFile("../../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_postVFP_EGM2D.root");


                  // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2016postVFP_ULv2.root");
                  fL1L2_TRG_SF= new TFile("../Scale_factors/Trig_2016/trigSF_2D.root");
                  fL1L2_TRG_SFerr = new TFile("../Scale_factors/Trig_2016/trigSF_2Derr.root");

               }
               else{
                  // cout<<" pre  2016="<<endl;
                  fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root");
                  fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root");
                  
                  fL2_Reco_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root");
                  fL2_ISO_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root");
                  fL2_ID_SF= new TFile("../../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_preVFP_EGM2D.root");


                  // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2016preVFP_ULv2.root");
                  fL1L2_TRG_SF= new TFile("../Scale_factors/Trig_2016_pre/trigSF_2D.root");
                  fL1L2_TRG_SFerr = new TFile("../Scale_factors/Trig_2016_pre/trigSF_2Derr.root");
               }
               
            }// else 2016

         fh2DL1_Reco_SF1 = (TH2F*)(fL1_Reco_SF->Get("NUM_TrackerMuons_DEN_genTracks"));
         fh2DL1_ID_SF1 = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt"));//
         fh2DL1_ISO_SF1 = (TH2F*)(fL1_ISO_SF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"));
         // fh2DL2_ISO_SF1 = (TH2F*)(fL2_ISO_SF->Get("NUM_probe_MiniIsoTight_DEN_probe_isTight_abseta_pt"));
         //$$$$ 
         fh2DL2_Reco_SF1 = (TH2F*)(fL2_Reco_SF->Get("EGamma_SF2D"));
         fh2DL2_ISO_SF1 = (TH2F*)(fL2_ISO_SF->Get("EGamma_SF2D"));
         fh2DL2_ID_SF1 = (TH2F*)(fL2_ID_SF->Get("EGamma_SF2D"));
         

         fh2DL1_ID_SF1err  = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_syst"));//
         fh2DL2_ID_SF1err  = (TH2F*)(fL2_ID_SF->Get("statMC"));
         fh2DL1_ISO_SF1err  = (TH2F*)(fL1_ISO_SF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_syst"));
         fh2DL2_ISO_SF1err  = (TH2F*)(fL2_ISO_SF->Get("statMC"));



         // fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("h2D_SF_emu_lepABpt_FullError")); // x- axiss pt and y axis pt 
         fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("c1")); // x- axiss pt and y axis pt 
         fhL1L2_TRG_SFerr = (TH2F*)(fL1L2_TRG_SFerr->Get("c1")); // x- axiss pt and y axis pt 
      }
   else if (MuMu > 0 && isMC) // SingleMuon or DiMuon
      {
         if (Year == 2018)
            {
               //cout<< " 2018 year="<<endl;
               fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../../Scale_factors/Muon/Run2018_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Run2018_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
               
               fL2_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL2_ID_SF= new TFile("../../Scale_factors/Muon/Run2018_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL2_ISO_SF= new TFile("../../Scale_factors/Muon/Run2018_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

               fL1L2_TRG_SF = new TFile("../../Scale_factors/Muon/Run2018_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
               fL1L2_TRG_SFerr = new TFile("../../Scale_factors/Muon/Run2018_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");

               // histos associated to files
               fh2DL1_Reco_SF1 = (TH2F*)(fL1_Reco_SF->Get("NUM_TrackerMuons_DEN_genTracks"));
               fh2DL2_Reco_SF1 = (TH2F*)(fL2_Reco_SF->Get("NUM_TrackerMuons_DEN_genTracks"));

               fh2DL1_ID_SF1 = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt"));//
               fh2DL2_ID_SF1 = (TH2F*)(fL2_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt"));
               fh2DL1_ID_SF1err  = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_combined_syst"));//
               fh2DL2_ID_SF1err  = (TH2F*)(fL2_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_combined_syst"));


               fh2DL1_ISO_SF1 = (TH2F*)(fL1_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt"));
               fh2DL2_ISO_SF1 = (TH2F*)(fL2_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt"));
               fh2DL1_ISO_SF1err  = (TH2F*)(fL1_ISO_SF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_combined_syst"));
               fh2DL2_ISO_SF1err  = (TH2F*)(fL2_ISO_SF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_combined_syst"));
         
               fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt")); // x- axiss pt and y axis pt 
               fhL1L2_TRG_SFerr = (TH2F*)(fL1L2_TRG_SFerr->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt_combined_syst"));

            }
         if(Year == 2017)
            {
               //cout<<" year 2017"<<endl;
               fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
               
               fL2_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL2_ID_SF= new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL2_ISO_SF= new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

               fL1L2_TRG_SF = new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
               fL1L2_TRG_SFerr = new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
               // histos associated to files

               fh2DL1_Reco_SF1 = (TH2F*)(fL1_Reco_SF->Get("NUM_TrackerMuons_DEN_genTracks"));
               fh2DL2_Reco_SF1 = (TH2F*)(fL2_Reco_SF->Get("NUM_TrackerMuons_DEN_genTracks"));

               fh2DL1_ID_SF1 = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt"));//
               fh2DL2_ID_SF1 = (TH2F*)(fL2_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt"));
               fh2DL1_ID_SF1err  = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_combined_syst"));//
               fh2DL2_ID_SF1err  = (TH2F*)(fL2_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_combined_syst"));

               fh2DL1_ISO_SF1 = (TH2F*)(fL1_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt"));
               fh2DL2_ISO_SF1 = (TH2F*)(fL2_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt"));
               fh2DL1_ISO_SF1err  = (TH2F*)(fL1_ISO_SF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_combined_syst"));
               fh2DL2_ISO_SF1err  = (TH2F*)(fL2_ISO_SF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_combined_syst"));

               fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt")); // x- axiss pt and y axis pt 
               //$$$$ 

               fhL1L2_TRG_SFerr = (TH2F*)(fL1L2_TRG_SFerr->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt_combined_syst"));



            }
         if(Year == 2016)
            {
               if (isPostAPV){
                  //cout<<" post 2016="<<endl;
                  fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
                  
                  fL2_Reco_SF =  new TFile("../../Scale_factors/Muon/Run2016_UL/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
                  fL2_ID_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL2_ISO_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

                  fL1L2_TRG_SF = new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
                  fL1L2_TRG_SFerr = new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");



               // histos associated to files

               }
               else{
                  // cout<<" pre  2016="<<endl;
                  fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
                  
                  fL2_Reco_SF =  new TFile("../../Scale_factors/Muon/Run2016_UL/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
                  fL2_ID_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL2_ISO_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

                  fL1L2_TRG_SF = new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
                  fL1L2_TRG_SFerr = new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
                  // histos associated to files

               }
               
               // histos associated to files
               fh2DL1_Reco_SF1 = (TH2F*)(fL1_Reco_SF->Get("NUM_TrackerMuons_DEN_genTracks"));
               fh2DL2_Reco_SF1 = (TH2F*)(fL2_Reco_SF->Get("NUM_TrackerMuons_DEN_genTracks"));

               fh2DL1_ID_SF1 = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt"));//
               fh2DL2_ID_SF1 = (TH2F*)(fL2_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt"));
               fh2DL1_ID_SF1err  = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_combined_syst"));//
               fh2DL2_ID_SF1err  = (TH2F*)(fL2_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_combined_syst"));

               fh2DL1_ISO_SF1 = (TH2F*)(fL1_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt"));
               fh2DL2_ISO_SF1 = (TH2F*)(fL2_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt"));
               fh2DL1_ISO_SF1err  = (TH2F*)(fL1_ISO_SF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_combined_syst"));
               fh2DL2_ISO_SF1err  = (TH2F*)(fL2_ISO_SF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt_combined_syst"));


               fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt")); // x- axiss pt and y axis pt 
               fhL1L2_TRG_SFerr = (TH2F*)(fL1L2_TRG_SFerr->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt_combined_syst"));
               //$$$$ 
            }// else 2016
      }

      //---------------------------------------------------------------------------//
      // !! ----------------------------- End of Scale Factors---------------------------------//
      //---------------------------------------------------------------------------//

std::vector<TString> NJET = {"0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"};

  Long64_t nbytes = 0, nb = 0;
  cout<< "Line : "  << __LINE__ << " " << nentries2 << endl; 
  cout<< " XS : "<<XS<<endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {//nentries2
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      allevents++;
      if ( allevents%10000 == 0 )  std::cout << "events : " << allevents << std::endl;
      // cout<< " count "<<jentry<<endl;
      // if (jentry > 100000) break;
            
      if ( signal  ) 
         {
          if (  minitree_nLLP->at(0) != 2)  continue; // protection against rare wrong signal events    //
         } 

      if (minitree_Mmumu->at(0) < 10) continue;
 
      if ( !minitree_Filter->at(0) && !minitree_FilterSameSign->at(0) ) continue;//

      nFilterEvt++;
      if (minitree_njetNOmu->at(0) < 1) continue;
      nFilterJet++;
      if ( minitree_Hemi_pt->size() != 2 ) continue;
      nFilternHemi++;
      if ( minitree_Hemi_pt->at(0) < lowHpt || minitree_Hemi_pt->at(1) < lowHpt ) continue;
      nFilterHpt++;
      fillHisto("hData_Filter","",samplename, minitree_Filter->at(0),1 );
            
      // --------- Lepton SF ------//
      double Mu_SF=1, Ele_SF=1, Mu_SF2 = 1;
      double EMu_trig_SF1=1;
      double triggerSF = 1;
      double triggerSFerr=0;
      double  Mu_ID_SF1=1, Mu_ISO_SF1=1, Mu_Reco_SF1=1,Mu_ID_SF2=1, Mu_ISO_SF2=1, Mu_Reco_SF2=1, Ele_Reco_SF1=1, Ele_ID_SF1=1, Ele_Reco_SF2=1;//$$$$
      double  Mu_ID_SF1err=0, Mu_ISO_SF1err=0, Mu_ID_SF2err=0, Mu_ISO_SF2err=0,  Ele_ID_SF1err=0, Ele_Reco_SF1err = 0;//$$$$
      float mu1_pt=0;
      float mu2_pt = 0;
	   float ele1_pt=0;
	   float ele1_eta=0;
      float SFele_pt_eta = 1.;
      float SFele_eta = 1.;

      //  MuMu = Channel; // Emu = 0 ; SingleMuon = 1; DiMuon = 2

            for (unsigned int iMuon = 0; iMuon < minitree_lepton_leadingpt->size(); iMuon++)
               {
                  if (minitree_lepton_leadingpt->at(0) == 0   )
                     {
                        Mu_Reco_SF1 = 1;
                        Mu_ID_SF1 = 1;	    
                        Mu_ISO_SF1 = 1 ;
                        continue;
                     }
                  if (MuMu == 0 && isMC)
                     {// leading letpn is muon by default

                              mu1_pt = minitree_lepton_leadingpt->at(iMuon);
                              Mu_Reco_SF1= fh2DL1_Reco_SF1->GetBinContent(fh2DL1_Reco_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_Reco_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));
                              Mu_ID_SF1 = 1;//fh2DL1_ID_SF1->GetBinContent(fh2DL2_ID_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ID_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));	    
                              Mu_ID_SF1err = 0.;
                              Mu_ISO_SF1 = fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));//fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));
                              Mu_ISO_SF1err = fh2DL1_ISO_SF1err->GetBinContent(fh2DL1_ISO_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));//fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));
                              if(minitree_lepton_leadingpt->at(iMuon) <= 15. ) {   Mu_Reco_SF1=1.; Mu_ID_SF1 = 1.;  Mu_ISO_SF1=1.; Mu_ID_SF1err = 0.;  Mu_ISO_SF1err=0.; } //muon ID, reco,ISO, Y axis  pt range 15-200,  x axis eta 0-2.4
                              if(minitree_lepton_leadingpt->at(iMuon) > 40  ) { Mu_Reco_SF1 = 1. ;} 
                              if(minitree_lepton_leadingpt->at(iMuon) >= 120. ) { Mu_ID_SF1 = 1.; Mu_ISO_SF1=1.;Mu_ID_SF1err = 0.; Mu_ISO_SF1err=0.;}

                     }
                  else if (MuMu > 0 && isMC) // $$$$$
                     {
                        
                        mu1_pt = minitree_lepton_leadingpt->at(iMuon);
                        mu2_pt = minitree_lepton_leadingpt2->at(iMuon);

                        Mu_Reco_SF1= fh2DL1_Reco_SF1->GetBinContent(fh2DL1_Reco_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_Reco_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));
                        Mu_ID_SF1 = fh2DL1_ID_SF1->GetBinContent(fh2DL1_ID_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ID_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));	    
                        Mu_ISO_SF1 = fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));
                        Mu_ID_SF1err = fh2DL1_ID_SF1err->GetBinContent(fh2DL1_ID_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ID_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));	    
                        Mu_ISO_SF1err = fh2DL1_ISO_SF1err->GetBinContent(fh2DL1_ISO_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt->at(iMuon)));

                        
                        Mu_Reco_SF2= fh2DL1_Reco_SF1->GetBinContent(fh2DL1_Reco_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_Reco_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iMuon)));
                        Mu_ID_SF2 = fh2DL1_ID_SF1->GetBinContent(fh2DL1_ID_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_ID_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iMuon)));	    
                        Mu_ISO_SF2 = fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iMuon)));
                        Mu_ID_SF2err = fh2DL1_ID_SF1err->GetBinContent(fh2DL1_ID_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_ID_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iMuon)));	    
                        Mu_ISO_SF2err = fh2DL1_ISO_SF1err->GetBinContent(fh2DL1_ISO_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_ISO_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iMuon)));


                        if(minitree_lepton_leadingpt->at(iMuon) <= 15. ) {   Mu_Reco_SF1=1.; Mu_ID_SF1 = 1.;  Mu_ISO_SF1=1.;  Mu_ID_SF1err = 0.;  Mu_ISO_SF1err=0.;  } //muon ID, reco,ISO, Y axis  pt range 15-200,  x axis eta 0-2.4
                        if(minitree_lepton_leadingpt->at(iMuon) > 40  ) { Mu_Reco_SF1 = 1. ;} 
                        if(minitree_lepton_leadingpt->at(iMuon) >= 120. ) { Mu_ID_SF1 = 1.; Mu_ISO_SF1=1.;Mu_ID_SF1err = 0.; Mu_ISO_SF1err=0.;}

                        
                        if(minitree_lepton_leadingpt2->at(iMuon) <= 15. ) {   Mu_Reco_SF2=1.; Mu_ID_SF2 = 1.;  Mu_ISO_SF2=1.; Mu_ID_SF2err = 0.;  Mu_ISO_SF2err=0.; } //muon ID, reco,ISO, Y axis  pt range 15-200,  x axis eta 0-2.4
                        if(minitree_lepton_leadingpt2->at(iMuon) > 40  ) { Mu_Reco_SF2 = 1. ;} 
                        if(minitree_lepton_leadingpt2->at(iMuon) >= 120. ) { Mu_ID_SF2 = 1.; Mu_ISO_SF2=1.;Mu_ID_SF2err = 0.; Mu_ISO_SF2err=0.;}
                        
                     }

               }
            if (MuMu == 0 && isMC)
               {
                     for(unsigned int iel = 0; iel <minitree_lepton_leadingpt2->size(); iel ++)//test    after trigger and lepton selection cut                                    
                        {
                          if (minitree_lepton_leadingpt2->at(0) == 0  )
                           {
                              Ele_Reco_SF1 = 1;
                              Ele_ID_SF1 = 1;	
                              Ele_ID_SF1err = 0;   
                              Ele_Reco_SF1err = 0; 
                              continue;
                           }

                           ele1_pt = minitree_lepton_leadingpt2->at(iel);
                           ele1_eta = fabs(minitree_lepton_leadingeta2->at(iel));
                           if ( ele1_pt < 20)
                              {  
                                 Ele_Reco_SF1 = fh2DL2_ISO_SF1->GetBinContent(fh2DL2_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_ISO_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                                 Ele_Reco_SF1err = fh2DL2_ISO_SF1err->GetBinContent(fh2DL2_ISO_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_ISO_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                              }
                           else
                              {
                                 Ele_Reco_SF1 = fh2DL2_Reco_SF1->GetBinContent(fh2DL2_Reco_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_Reco_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                                 Ele_Reco_SF1err = fh2DL2_Reco_SF1->GetBinContent(fh2DL2_Reco_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_Reco_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                              }
                           Ele_ID_SF1 = fh2DL2_ID_SF1->GetBinContent(fh2DL2_ID_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_ID_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                           Ele_ID_SF1err = fh2DL2_ID_SF1err->GetBinContent(fh2DL2_ID_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_ID_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));

                           if(minitree_lepton_leadingpt2->at(iel) <= 10. ) {  Ele_ID_SF1=1.;Ele_ID_SF1err=0.;Ele_Reco_SF1err = 0; }// ID, 10-500, eta -2.5 to 2.5, x axis = eta, y axis=pt   
                           if(minitree_lepton_leadingpt2->at(iel) <= 10. ) {  Ele_Reco_SF1 = 1.;  }// Reco , 10-500, eta -2.5 to 2.5, x axis = eta, y axis=pt
                           if(minitree_lepton_leadingpt2->at(iel) >= 500. ) { Ele_Reco_SF1 = 1.; Ele_ID_SF1err=0.; Ele_Reco_SF1err = 0;}

                        }
               }                        


            //cout<<" mu pt inside ="<<mu1_pt<<" ele pt outside ="<<ele1_pt<<endl;
         if (isMC && MuMu == 0) {   
                                  triggerSF  = fhL1L2_TRG_SF->GetBinContent(fhL1L2_TRG_SF->GetXaxis()->FindBin(ele1_pt),fhL1L2_TRG_SF->GetYaxis()->FindBin(mu1_pt));
                                 triggerSFerr = fhL1L2_TRG_SFerr->GetBinContent(fhL1L2_TRG_SFerr->GetXaxis()->FindBin(ele1_pt),fhL1L2_TRG_SFerr->GetYaxis()->FindBin(mu1_pt));
                                 if ((ele1_pt < 15. && mu1_pt < 15.) ||  (ele1_pt > 200 && mu1_pt > 200)){triggerSF=1.;triggerSFerr=0.;} // x and y axis 15 to 500
                                 }
         else if (isMC && MuMu >0 ) 
                  {
                     triggerSF = fhL1L2_TRG_SF->GetBinContent(fhL1L2_TRG_SF->GetXaxis()->FindBin(mu1_pt),fhL1L2_TRG_SF->GetYaxis()->FindBin(mu2_pt)); 
                     triggerSFerr = fhL1L2_TRG_SFerr->GetBinContent(fhL1L2_TRG_SFerr->GetXaxis()->FindBin(mu1_pt),fhL1L2_TRG_SFerr->GetYaxis()->FindBin(mu2_pt));
                     if ((mu1_pt < 15. && mu2_pt < 15.) ||  (mu1_pt > 500 && mu2_pt > 500)){triggerSF=1.;triggerSFerr=0.;}
                  }// !! 

            // std::cout<<" here 1 "<<std::endl;
      //---------------------------------------------------------------------------//
      if (MuMu >= 1 )
         {
            LT= mu1_pt + mu2_pt;
         }
      else if (MuMu == 0 )
         {
            LT= mu1_pt + ele1_pt;
         }
      //---------------------------------------------------------------------------//
      // !! -----------------------------Event Weight ---------------------------------//
      //---------------------------------------------------------------------------//  
      //---------------------------------------------------------------------------//


      // if (showoutput) // std::cout<<"GenWeight : "<<minitree_only_gen_wt->at(0)<<std::endl;
      float GenWeight = minitree_only_gen_wt->at(0)/mean;
      if (Signal || !isMC)
         {
            GenWeight = 1.;
         }
      fillHisto("GenWeight","", samplename, GenWeight,1);
      float top_pt_wt = 1.;
      float Prefweight = miniPrefweight->at(0);
      float PrefweightUp = 1.;
      float PrefweightDown = 1.;
      float PUweight = miniPUweight->at(0);
      float PUweightUp = 1.;
      float PUweightDown = 1.;
      float TriggerSyst = triggerSF;
      Ele_SF =Ele_Reco_SF1*Ele_ID_SF1;
      Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * Mu_ISO_SF1;
      Mu_SF2 = Mu_Reco_SF2* Mu_ID_SF2 * Mu_ISO_SF2;
      if (Signal)
         {
            PUweightUp = miniPUweight_Up->at(0);
            PUweightDown = miniPUweight_Down->at(0);
         }
       
      double NormFactor=1;
      double NormFactor_wo_top_pt=1;
      double NormFactor_full=1;
      double NormFactor_mu=1;
      double NormFactor_ele=1;
      NormFactorLumi =  XS*Lumi;  
      NormFactorLumiUp = XS*LumiUp;
      NormFactorLumiDown = XS*LumiDown;

      NormFactorXS =  XS*Lumi;  
      NormFactorXSUp = XS_up*Lumi;
      NormFactorXSDown = XS_down*Lumi;

      // rajouter l1 prefiring
      float NormFactorSYST = NormFactorLumi/norm;
         PUweight = miniPUweight->at(0);
         PUweightUp = miniPUweight_Up->at(0);
         PUweightDown = miniPUweight_Down->at(0);
         Prefweight = miniPrefweight->at(0);
         PrefweightUp = miniPrefweight_Up->at(0);
         PrefweightDown = miniPrefweight_Down->at(0);
         NormFactorSYST =  NormFactorLumi/norm;
         // It is supposed to run for several systematics sequentially, however the code oes not work when lloping on  different systematics ...
         // so we have to put the 0 instead of the index of the syst in the vector systlist
         if (systlist[0].Contains("LumiUp"))
            {
               NormFactorSYST =  NormFactorLumiUp/norm;
            }
         else if (systlist[0].Contains("LumiDown"))
            {
               NormFactorSYST =  NormFactorLumiDown/norm;
            }
         else if (systlist[0].Contains("XSUp"))
            {
               NormFactorSYST =  NormFactorXSUp/norm;
            }
         else if (systlist[0].Contains("XSDown"))
            {
               NormFactorSYST =  NormFactorXSDown/norm;
            }
         else if (systlist[0].Contains("L1Up"))
            {
               Prefweight = PrefweightUp;
            }
         else if (systlist[0].Contains("L1Down"))
            {
               Prefweight = PrefweightDown;
            }
         else if (systlist[0].Contains("TriggerDown"))
            {
               TriggerSyst = triggerSF-triggerSFerr;
               
            }
         else if (systlist[0].Contains("TriggerUp"))
            {
               TriggerSyst = triggerSF+triggerSFerr;
            }
         else if (systlist[0].Contains("LepIDUp"))
            {
               if (MuMu == 0)
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1+abs(Mu_ID_SF1-Mu_ID_SF1err)) * Mu_ISO_SF1;
                     Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1+Ele_ID_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1+Mu_ID_SF1err) * Mu_ISO_SF1;
                     Mu_SF2 = Mu_Reco_SF2* (Mu_ID_SF2+Mu_ID_SF2err) * Mu_ISO_SF2;
                  }
            }    
         else if (systlist[0].Contains("LepIDDown"))
            {
               if (MuMu == 0)
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1- abs(Mu_ID_SF1-Mu_ID_SF1err)) * Mu_ISO_SF1;
                     Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1-Ele_ID_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1-Mu_ID_SF1err) * Mu_ISO_SF1;
                     Mu_SF2 = Mu_Reco_SF2* (Mu_ID_SF2-Mu_ID_SF2err) * Mu_ISO_SF2;
                  }
            }  
         else if (systlist[0].Contains("LepISOUp"))
            {
               if (MuMu == 0)
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1+abs(Mu_ISO_SF1-Mu_ISO_SF1err));
                     Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1+Ele_Reco_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1+Mu_ISO_SF1err);
                     Mu_SF2 = Mu_Reco_SF2* Mu_ID_SF2 * (Mu_ISO_SF2+Mu_ISO_SF2err);
                  }
            } 
         else if (systlist[0].Contains("LepISODown"))
            {
               if (MuMu == 0)
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1-abs(Mu_ISO_SF1-Mu_ISO_SF1err));
                     Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1-Ele_Reco_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1-Mu_ISO_SF1err);
                     Mu_SF2 = Mu_Reco_SF2* Mu_ID_SF2 * (Mu_ISO_SF2-Mu_ISO_SF2err);
                  }
            } 
         else  if (systlist[0].Contains("PUUp"))
            {
               PUweight = PUweightUp;
            }
         else if (systlist[0].Contains("PUDown"))
            {
               PUweight = PUweightDown;
            }
         else
            {
               NormFactorSYST =  NormFactorLumi/norm;
            }


      if(thesample.Contains("TTTo2L2Nu") || thesample.Contains("TTToHadronic") ||  thesample.Contains("TTToSemiLeptonic")){
         top_pt_wt=minitree_genTop_Weight->at(0);
         //top_pt_wt=1;
      }
      else{
	      top_pt_wt=1;
      }
      if (!isMC)
         {
            NormFactor = 1;//found on XSDB                                                                                                                                   
            NormFactor_mu = 1;
            NormFactor_mu2 = 1;
            NormFactor_ele = 1;
            Pref_PU_gen_wt=1;
            top_pt_wt=1;

            Mu_SF =1;
            Ele_SF =1;
            EMu_trig_SF1=1;
            triggerSFerr = 0;
            Ele_ID_SF1err = 0;
            Mu_ID_SF1err = 0;
            Mu_ISO_SF1err = 0;
            Mu_ID_SF2err = 0;
            Mu_ISO_SF2err = 0;
         }
      else{
         
            Pref_PU_gen_wt= NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt;
            NormFactor_mu =  SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Mu_SF*TriggerSyst;
            NormFactor_mu2 =  SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Mu_SF2*TriggerSyst;
            NormFactor_ele=  SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Ele_SF*TriggerSyst;
            NormFactor =  SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Mu_SF*Ele_SF*TriggerSyst;
         }
   if (Signal)
      {
            Pref_PU_gen_wt= NormFactorSYST*Prefweight*PUweight*top_pt_wt;
            NormFactor_mu =  SFele_eta*NormFactorSYST*Prefweight*PUweight*top_pt_wt*Mu_SF*TriggerSyst;
            NormFactor_mu2 =  SFele_eta*NormFactorSYST*Prefweight*PUweight*top_pt_wt*Mu_SF2*TriggerSyst;
            NormFactor_ele=  SFele_eta*NormFactorSYST*Prefweight*PUweight*top_pt_wt*Ele_SF*TriggerSyst;
            NormFactor =  SFele_eta*NormFactorSYST*Prefweight*PUweight*top_pt_wt*Mu_SF*Ele_SF*TriggerSyst;
      }
      fillHisto("hData_Event_Weight","",samplename, Pref_PU_gen_wt,1 );

      // std::cout<<" here 2 "<<std::endl;

      //--------------------------------------------------------------//
      bool isHemiVtx1 = false, isHemiVtx2 = false;
      bool isCutVtx = false, isCutVtx1 = false, isCutVtx2 = false;
      bool isCutEvt = false;
      float BDTvtx = -2., BDTvtx1 = -2., BDTvtx2 = -2.;
      bool ping;
      bool isHemiVtx1Loose = false, isHemiVtx2Loose = false;
      int nVtx = 0, nVtxIni = 0, step;
      int nVtxLoose = 0 , nVtxIniLoose = 0; 
      float VtxMass = 0., dR, dist, NChi2, r, eta;
      float Vtx_HMass = 0;
      float Vtx_nTrks= 0;
      float Vtx_z= 0;
      float Vtx_r= 0;
      float Vtx_dR= 0;
      float Vtx_SumtrackWeight= 0;
      float Vtx_track_MeanDCA_d= 0;
      float Vtx_dist = 0;
      float Vtx_NChi = 0;

      bool ping0 = false;
      bool ping1 = false;
      float dR0 = 0.;
      float dR1 = 0.;
      ////////////////// hemisphere pT
      float  hemi1_pt  = -1.;
      float  hemi2_pt  = -1.;
      ////////////////////////////////
      
      fillHisto("hData_Mmumu","",samplename, minitree_Mmumu->at(0),Pref_PU_gen_wt );

      int FilterSample = -1;
      // if (showoutput) // std::cout<<"minitree_Mmumu->at(0) : "<<minitree_Mmumu->at(0)<<std::endl;
      if (MuMu == 2 || MuMu == 1)
         {
            // // std::cout<<"minitree_trigger_doublelepton->at(0) : "<<minitree_trigger_doublelepton->at(0)<<std::endl;
            // // std::cout<<"minitree_trigger_singlelepton->at(0) : "<<minitree_trigger_singlelepton->at(0)<<std::endl;
            if (DoubleMuon && minitree_trigger_doublelepton->at(0)  ) //&& !minitree_trigger_singlelepton->at(0)
               {
                  FilterSample = 2;
                        // // std::cout<<"FS 2  "<<std::endl;
               }
            else if  (!DoubleMuon && minitree_trigger_singlelepton->at(0) &&  !minitree_trigger_doublelepton->at(0) ) 
               {
                  FilterSample = 1;
                        // // std::cout<<"FS 1 "<<std::endl;
               }
         }
      if (MuMu == 0)
         {
               FilterSample = 0;
         }
      // // std::cout<<"FilterSample : "<<FilterSample<<std::endl;
      if (DoubleMuon && FilterSample != 2) continue;
      else if (!DoubleMuon && FilterSample != 1 && MuMu ==1) continue;
      else if (!DoubleMuon && FilterSample != 0 && MuMu ==0) continue;

         // std::cout<<" here 3 "<<std::endl;

      //--------------------------------------------------------//
      //--------------------------------------------------------//
      bool isSS = SS;// decide if you want SS category or not
      bool isFWD = FWD;// decide if you want FWD category or not
      //**//
      bool Filter = false;
      if (isSS)
         {
           Filter = minitree_FilterSameSign->at(0);
         }
      else  
         {
            Filter = minitree_Filter->at(0);
         }

      if (isFWD)
         {
            if ( (abs(minitree_Hemi_eta->at(0)) < 2.4 || abs(minitree_Hemi_eta->at(0)) > 3.0)  && ( abs(minitree_Hemi_eta->at(1)) < 2.4 || abs(minitree_Hemi_eta->at(1)) > 3.0) ) continue; 
         }
      else
         {
            if ( abs(minitree_Hemi_eta->at(0)) > 2.4 && abs(minitree_Hemi_eta->at(1)) > 2.4 ) continue;
         }

      //--------------------------------------------------------//
      if ( minitree_njetNOmu->at(0) < 1 ) Filter = false; 
      hemi1_pt = minitree_Hemi_pt->at(0);
      hemi2_pt = minitree_Hemi_pt->at(1);



      //--------------------------------------------------------//
      //--------------------------------------------------------//

      // Forward region
      //       if(Filter && ((abs(minitree_Hemi_eta->at(0)) > 2.4 && abs(minitree_Hemi_eta->at(0)) < 3.0) ||
      //    (abs(minitree_Hemi_eta->at(1)) > 2.4 && abs(minitree_Hemi_eta->at(1)) < 3.0)))
      // std::cout<<" here 4 "<<std::endl;

      float hemi_ptmin = hemi2_pt;
      float hemi_ptmax = hemi1_pt;
      if ( hemi2_pt > hemi1_pt ) {hemi_ptmin = hemi1_pt;hemi_ptmax = hemi2_pt;}
      fillHisto("Hemisphere_leadingpt","", samplename, hemi_ptmax,Pref_PU_gen_wt);
      fillHisto("Hemisphere_subleadingpt","", samplename, hemi_ptmin,Pref_PU_gen_wt);

   if (signal)
      {
         ping0 = minitree_Hemi_LLP_ping->at(0);
         ping1 = minitree_Hemi_LLP_ping->at(1);
         dR0 = minitree_Hemi_LLP_dR->at(0);
         dR1 = minitree_Hemi_LLP_dR->at(1);
      }

   int Vtx_step0 = minitree_Hemi_Vtx_step->at(0);
   int Vtx_step1 = minitree_Hemi_Vtx_step->at(1);
   float Vtx_NChi0 = minitree_Hemi_Vtx_NChi2->at(0);
   float Vtx_NChi1 = minitree_Hemi_Vtx_NChi2->at(1);
   float Vtx_Mass0 = minitree_Hemi_Vtx_Mass->at(0);
   float Vtx_Mass1 = minitree_Hemi_Vtx_Mass->at(1);
   float Vtx_dist0 = minitree_Hemi_Vtx_dist->at(0);
   float Vtx_dist1 = minitree_Hemi_Vtx_dist->at(1);

   float Vtx_HMass0 = minitree_Hemi_Vtx_BDT_HMass->at(0);
   float Vtx_HMass1 = minitree_Hemi_Vtx_BDT_HMass->at(1);
   float Vtx_nTrks0 = minitree_Hemi_Vtx_nTrks->at(0);
   float Vtx_nTrks1 = minitree_Hemi_Vtx_nTrks->at(1);
   float Vtx_z0 = minitree_Hemi_Vtx_z->at(0);
   float Vtx_z1 = minitree_Hemi_Vtx_z->at(1);
   float Vtx_r0 = minitree_Hemi_Vtx_r->at(0);
   float Vtx_r1 = minitree_Hemi_Vtx_r->at(1);
   float Vtx_dR0 = minitree_Hemi_Vtx_dR->at(0);
   float Vtx_dR1 = minitree_Hemi_Vtx_dR->at(1);
   float Vtx_SumtrackWeight0 = minitree_Hemi_Vtx_SumtrackWeight->at(0);
   float Vtx_SumtrackWeight1 = minitree_Hemi_Vtx_SumtrackWeight->at(1);

   float Vtx_MeantrackWeight0 = 0;
   float Vtx_MeantrackWeight1 = 0;
   float Vtx_MeantrackWeight = 0;
   if ( minitree_Hemi_Vtx_nTrks->at(0) != 0){Vtx_MeantrackWeight0 = minitree_Hemi_Vtx_SumtrackWeight->at(0)/static_cast< float >(minitree_Hemi_Vtx_nTrks->at(0));}
   if ( minitree_Hemi_Vtx_nTrks->at(1) != 0){Vtx_MeantrackWeight1 = minitree_Hemi_Vtx_SumtrackWeight->at(1)/static_cast< float >(minitree_Hemi_Vtx_nTrks->at(1));}


   float Vtx_track_MeanDCA_d0 = minitree_Hemi_Vtx_track_MeanDCA_d->at(0);
   float Vtx_track_MeanDCA_d1 = minitree_Hemi_Vtx_track_MeanDCA_d->at(1);
   float Vtx_Vtx_dist = 0;


   float posx0 = minitree_Hemi_Vtx_x->at(0);
   float posy0 = minitree_Hemi_Vtx_y->at(0);
   float posz0 = minitree_Hemi_Vtx_z->at(0);
   float posx1 = minitree_Hemi_Vtx_x->at(1);
   float posy1 = minitree_Hemi_Vtx_y->at(1);
   float posz1 = minitree_Hemi_Vtx_z->at(1);
   float r0 = TMath::Sqrt( posx0*posx0 + posy0*posy0 );
   float z0 = TMath::Abs( posz0 );
   float r1 = TMath::Sqrt( posx1*posx1 + posy1*posy1 );
   float z1 = TMath::Abs( posz1 );
   float recX0 = posx0 - minitree_PV_x->at(0);
   float recY0 = posy0 - minitree_PV_y->at(0);
   float recZ0 = posz0 - minitree_PV_z->at(0);
   float recX1 = posx1 - minitree_PV_x->at(0);
   float recY1 = posy1 - minitree_PV_y->at(0);
   float recZ1 = posz1 - minitree_PV_z->at(0);

  float theta_Vtx0 = TMath::ATan2(sqrt(recX0*recX0+recY0*recY0),abs(recZ0)) ;
  float theta_Vtx1 = TMath::ATan2(sqrt(recX1*recX1+recY1*recY1),abs(recZ1)) ;

  float eta_Vtx0 = -TMath::Log(tan(theta_Vtx0/2));
  float eta_Vtx1 = -TMath::Log(tan(theta_Vtx1/2));
  if ( posz0 < 0 ) eta_Vtx0 = -eta_Vtx0;
  if ( posz1 < 0 ) eta_Vtx1 = -eta_Vtx1;
   // // std::cout<<"Vtxvar : "<<std::endl;
   bool Merging = true; // false if "no merge" or "no close vtx" output
   bool Protect = false;
   //$$
   if (signal)
      {
         if (minitree_Hemi_SecLLP_ping->size() >= 1) Protect = true;
      }
   else
      {
         Protect = true;
      } 
    if ( Merging && Protect ) { 
            // protection again
        if ( minitree_Hemi_SecVtx->size() >= 1 ) {
          ping0 = false;
            if (signal)
                {ping0 = minitree_Hemi_SecLLP_ping->at(0);}
            
            Vtx_step0 = minitree_Hemi_SecVtx_step->at(0);
            Vtx_NChi0 = minitree_Hemi_SecVtx_NChi2->at(0);
            Vtx_Mass0 = minitree_Hemi_SecVtx_Mass->at(0);
            Vtx_dist0 = minitree_Hemi_SecVtx_dist->at(0);
            posx0 = minitree_Hemi_SecVtx_x->at(0);
            posy0 = minitree_Hemi_SecVtx_y->at(0);
            posz0 = minitree_Hemi_SecVtx_z->at(0);
            r0 = minitree_Hemi_SecVtx_r->at(0);
            Vtx_step1 = 0;
            Vtx_NChi1 = -1.;
            Vtx_Mass1 = 0.;
            Vtx_dist1 = 0.;
            ping1 = false;
            r1 = 0;
            float SecrecX0 = posx0 - minitree_PV_x->at(0);
            float SecrecY0 = posy0 - minitree_PV_y->at(0);
            float SecrecZ0 = posz0 - minitree_PV_z->at(0);
            float theta_SecVtx0 = TMath::ATan2(sqrt(SecrecX0*SecrecX0+SecrecY0*SecrecY0),abs(SecrecZ0)) ;

            eta_Vtx0 = -TMath::Log(tan(theta_SecVtx0/2.));
            eta_Vtx1 = 0;
            if ( posz0 < 0 ) eta_Vtx0 = -eta_Vtx0;
            Vtx_nTrks0 = minitree_Hemi_SecVtx_nTrks->at(0);
            Vtx_z0 = minitree_Hemi_SecVtx_z->at(0);
            Vtx_r0 = minitree_Hemi_SecVtx_r->at(0);
            Vtx_dR0 = minitree_Hemi_SecVtx_dR->at(0);
            Vtx_SumtrackWeight0 = minitree_Hemi_SecVtx_SumtrackWeight->at(0);
            if (minitree_Hemi_SecVtx_nTrks->at(0) !=0) Vtx_MeantrackWeight0 = minitree_Hemi_SecVtx_SumtrackWeight->at(0)/static_cast< float >(minitree_Hemi_SecVtx_nTrks->at(0));
            Vtx_track_MeanDCA_d0 = minitree_Hemi_SecVtx_track_MeanDCA_d->at(0);
        }
        if ( minitree_Hemi_SecVtx->size() == 2 ) {
            ping1 = false;
            if (signal)
                {ping1 = minitree_Hemi_SecLLP_ping->at(1);}
            
            Vtx_step1 = minitree_Hemi_SecVtx_step->at(1);
            Vtx_NChi1 = minitree_Hemi_SecVtx_NChi2->at(1);
            Vtx_Mass1 = minitree_Hemi_SecVtx_Mass->at(1);
            Vtx_dist1 = minitree_Hemi_SecVtx_dist->at(1);
            posx1 = minitree_Hemi_SecVtx_x->at(1);
            posy1 = minitree_Hemi_SecVtx_y->at(1);
            posz1 = minitree_Hemi_SecVtx_z->at(1);
            r1 = minitree_Hemi_SecVtx_r->at(1);
            float SecrecX1 = posx1 - minitree_PV_x->at(0);
            float SecrecY1 = posy1 - minitree_PV_y->at(0);
            float SecrecZ1 = posz1 - minitree_PV_z->at(0);
            float theta_SecVtx1 = TMath::ATan2(sqrt(SecrecX1*SecrecX1+SecrecY1*SecrecY1),abs(SecrecZ1)) ;
            eta_Vtx1 = -TMath::Log(tan(theta_SecVtx1/2.));
            if ( posz1 < 0 ) eta_Vtx1 = -eta_Vtx1;
            Vtx_nTrks1 = minitree_Hemi_SecVtx_nTrks->at(1);
            Vtx_z1 = minitree_Hemi_SecVtx_z->at(1);
            Vtx_r1 = minitree_Hemi_SecVtx_r->at(1);
            Vtx_dR1 = minitree_Hemi_SecVtx_dR->at(1);
            Vtx_SumtrackWeight1 = minitree_Hemi_SecVtx_SumtrackWeight->at(1);
            if (minitree_Hemi_SecVtx_nTrks->at(1) !=0) Vtx_MeantrackWeight1 = minitree_Hemi_SecVtx_SumtrackWeight->at(1)/static_cast< float >(minitree_Hemi_SecVtx_nTrks->at(1));
            Vtx_track_MeanDCA_d1 = minitree_Hemi_SecVtx_track_MeanDCA_d->at(1);

        }
      }// End of Merging information
      // // std::cout<<"End of Merging : "<<std::endl;
      isHemiVtx1 = false;
      isHemiVtx2 = false;
      isHemiVtx1Loose = false;
      isHemiVtx2Loose = false;  

      Vtx_Vtx_dist = sqrt((posx1-posx0)*(posx1-posx0)+(posy1-posy0)*(posy1-posy0)+
                            (posz1-posz0)*(posz1-posz0));

          //First Vertex
      if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 1 && Vtx_step0 <= 2 ) {
          isHemiVtx1 = true;
          VtxMass = Vtx_Mass0;
          BDTvtx1 = minitree_Hemi_Vtx_MVAval_Tight->at(0);
          BDTvtx  = BDTvtx1;

          Vtx_NChi = Vtx_NChi0;
          Vtx_dist = Vtx_dist0;
          Vtx_HMass = Vtx_HMass0;
          Vtx_nTrks = Vtx_nTrks0;
          Vtx_z = Vtx_z0;
          Vtx_r = Vtx_r0;
          Vtx_dR = Vtx_dR0;
          Vtx_SumtrackWeight = Vtx_SumtrackWeight0;
          Vtx_MeantrackWeight = Vtx_MeantrackWeight0;
          Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d0;


      }
      if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 3 && Vtx_step0 <= 4 ) {
          isHemiVtx1Loose = true;
          VtxMass = Vtx_Mass0;
          BDTvtx1 = minitree_Hemi_Vtx_MVAval_Loose->at(0);
          BDTvtx  = BDTvtx1;

          Vtx_NChi = Vtx_NChi0;
          Vtx_dist = Vtx_dist0;
          Vtx_HMass = Vtx_HMass0;
          Vtx_nTrks = Vtx_nTrks0;
          Vtx_z = Vtx_z0;
          Vtx_r = Vtx_r0;
          Vtx_dR = Vtx_dR0;
          Vtx_SumtrackWeight = Vtx_SumtrackWeight0;
          Vtx_MeantrackWeight = Vtx_MeantrackWeight0;
          Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d0;
      }

      //Second Vertex
      if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 1 && Vtx_step1 <= 2 ) {
          isHemiVtx2 = true;
          if ( Vtx_Mass1 > VtxMass ) VtxMass = Vtx_Mass1;
          BDTvtx2 = minitree_Hemi_Vtx_MVAval_Tight->at(1);

          if ( BDTvtx2 > BDTvtx )       BDTvtx      = BDTvtx2;
          if ( Vtx_NChi1 < Vtx_NChi ) Vtx_NChi      = Vtx_NChi1;
          if ( Vtx_dist1 > Vtx_dist ) Vtx_dist      = Vtx_dist1;
          if ( Vtx_HMass1 > Vtx_HMass ) Vtx_HMass   = Vtx_HMass1;
          if ( Vtx_nTrks1 > Vtx_nTrks ) Vtx_nTrks      = Vtx_nTrks1;
          if ( Vtx_r1 > Vtx_r )     Vtx_r           = Vtx_r1;
          if ( Vtx_z1 > Vtx_z )     Vtx_z           = Vtx_z1;
          if ( Vtx_dR1 > Vtx_dR )    Vtx_dR         = Vtx_dR1;
          if ( Vtx_SumtrackWeight1 > Vtx_SumtrackWeight )    Vtx_SumtrackWeight = Vtx_SumtrackWeight1;
          if ( Vtx_MeantrackWeight1 > Vtx_MeantrackWeight )    Vtx_MeantrackWeight = Vtx_MeantrackWeight1;
          if ( Vtx_track_MeanDCA_d1 < Vtx_track_MeanDCA_d ) Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d1;
      }

      if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 3 && Vtx_step1 <= 4 ) {
          isHemiVtx2Loose = true;
          if ( Vtx_Mass1 > VtxMass ) VtxMass = Vtx_Mass1;
          BDTvtx2 = minitree_Hemi_Vtx_MVAval_Loose->at(1);

          if ( BDTvtx2 > BDTvtx )       BDTvtx      = BDTvtx2;
          if ( Vtx_NChi1 < Vtx_NChi ) Vtx_NChi      = Vtx_NChi1;
          if ( Vtx_dist1 > Vtx_dist ) Vtx_dist      = Vtx_dist1;
          if ( Vtx_HMass1 > Vtx_HMass ) Vtx_HMass   = Vtx_HMass1;
          if ( Vtx_nTrks1 > Vtx_nTrks ) Vtx_nTrks      = Vtx_nTrks1;
          if ( Vtx_r1 > Vtx_r )     Vtx_r           = Vtx_r1;
          if ( Vtx_z1 > Vtx_z )     Vtx_z           = Vtx_z1;
          if ( Vtx_dR1 > Vtx_dR )    Vtx_dR         = Vtx_dR1;
          if ( Vtx_SumtrackWeight1 > Vtx_SumtrackWeight )    Vtx_SumtrackWeight = Vtx_SumtrackWeight1;
          if ( Vtx_MeantrackWeight1 > Vtx_MeantrackWeight )    Vtx_MeantrackWeight = Vtx_MeantrackWeight1;
          if ( Vtx_track_MeanDCA_d1 < Vtx_track_MeanDCA_d ) Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d1;
      }


      if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;// Tight TIght
      else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;// Tight and ???Loose or nothing

      if      ( isHemiVtx1Loose && isHemiVtx2Loose ) nVtxLoose = 2;// Loose and Loose
      else if ( isHemiVtx1Loose || isHemiVtx2Loose ) nVtxLoose = 1;// Loose and Tight or nothing

   nRecoTightVertex= nVtx;
   nRecoLooseVertex= nVtxLoose;

   float RATIOVTX = 0.;
   // We can say that a Tight vertex counts for 2, to increase the value that we cna observe in the ratio :D

   

   int nVtx1 = 0;
   int nVtx2 = 0;
   int nVtxLoose1 = 0;
   int nVtxLoose2 = 0;

   if (isHemiVtx1){nVtx1=1;}
   if (isHemiVtx2){nVtx2=1;}

   if (isHemiVtx1Loose){nVtxLoose1=1;}
   if (isHemiVtx2Loose){nVtxLoose2=1;}

   if (nVtx == 0 && nVtxLoose == 0)
      {
         RATIOVTX = 0;
         
      }
   if (nVtx == 2 && nVtxLoose == 0)
      {
         RATIOVTX = 4;
       
        fillHisto("hData_VtxQualityTight_Hemi1pt","2Vtx",samplename,hemi1_pt,1);
        fillHisto("hData_VtxQualityTight_Hemi2pt","2Vtx",samplename,hemi2_pt,1);

      }
   if (nVtx == 1 && nVtxLoose == 0)
      {
         RATIOVTX = 2;

         if (isHemiVtx1)
            {
               fillHisto("hData_VtxQualityTight_Hemi1pt","1Vtx",samplename,hemi1_pt,1);
            }
         else if (isHemiVtx2)
            {
               fillHisto("hData_VtxQualityTight_Hemi2pt","1Vtx",samplename,hemi2_pt,1);
            }
      }
   if (nVtx == 0 && nVtxLoose == 1)
      {
         RATIOVTX = 0.5;

         if (isHemiVtx1Loose)
            {
               fillHisto("hData_VtxQualityLoose_Hemi1pt","1Vtx",samplename,hemi1_pt,1);

            }
         else if (isHemiVtx2Loose)
            {
               fillHisto("hData_VtxQualityLoose_Hemi2pt","1Vtx",samplename,hemi2_pt,1);
            }
      }
   if (nVtx == 0 && nVtxLoose == 2)
      {
         RATIOVTX = 1;
         fillHisto("hData_VtxQualityLoose_Hemi1pt","2Vtx",samplename,hemi1_pt,1);
         fillHisto("hData_VtxQualityLoose_Hemi2pt","2Vtx",samplename,hemi2_pt,1);
      }
   if (nVtx == 1 && nVtxLoose == 1)
      {
         RATIOVTX = 3;


         if (isHemiVtx1)
            {
              fillHisto("hData_VtxQualityTight_Hemi1pt","2Vtx",samplename,hemi1_pt,1);
               // else if ( hemi1_pt > lowHpt && hemi1_pt< 80 ){fillHisto("hData_VtxQualityTight_LowHemi1pt","2Vtx",samplename,hemi1_pt,1);}
            }
         else if (isHemiVtx2)
            {
               fillHisto("hData_VtxQualityTight_Hemi2pt","2Vtx",samplename,hemi2_pt,1);
               // else if ( hemi2_pt > lowHpt && hemi2_pt< 80 ){fillHisto("hData_VtxQualityTight_LowHemi2pt","2Vtx",samplename,hemi2_pt,1);}
            }
         if (isHemiVtx1Loose)
            {
               fillHisto("hData_VtxQualityLoose_Hemi1pt","2Vtx",samplename,hemi1_pt,1);
               // else if ( hemi1_pt > lowHpt && hemi1_pt< 80 ){fillHisto("hData_VtxQualityLoose_LowHemi1pt","2Vtx",samplename,hemi1_pt,1);}
            }
         else if (isHemiVtx2Loose)
            {
               fillHisto("hData_VtxQualityLoose_Hemi2pt","2Vtx",samplename,hemi2_pt,1);
               // else if ( hemi2_pt > lowHpt && hemi2_pt< 80 ){fillHisto("hData_VtxQualityLoose_LowHemi2pt","2Vtx",samplename,hemi2_pt,1);}
            }

      }



   //-------------- //    
   float HemiAveragePt = (hemi_ptmax+hemi_ptmin)/2.;
   fillHisto("hData_VtxQualityTight_Hemi1pt","NoSel",samplename,hemi1_pt,nVtx1);
   fillHisto("hData_VtxQualityTight_Hemi2pt","NoSel",samplename,hemi2_pt,nVtx2);

   fillHisto("hData_VtxQualityLoose_Hemi1pt","NoSel",samplename,hemi1_pt,nVtxLoose1);
   fillHisto("hData_VtxQualityLoose_Hemi2pt","NoSel",samplename,hemi2_pt,nVtxLoose2);

   //-------------- //
   fillHisto("hData_VtxQualityTight_Hemileadingpt","NoSel",samplename,hemi_ptmax,nRecoTightVertex);
   fillHisto("hData_VtxQualityTight_Hemisubleadingpt","NoSel",samplename,hemi_ptmin,nRecoTightVertex);
  
   fillHisto("hData_VtxQualityLoose_Hemileadingpt","NoSel",samplename,hemi_ptmax,nRecoLooseVertex);
   fillHisto("hData_VtxQualityLoose_Hemisubleadingpt","NoSel",samplename,hemi_ptmin,nRecoLooseVertex);
   
   fillHisto("hData_VtxQualityTight_HemiAveragept","NoSel",samplename,HemiAveragePt,nRecoTightVertex);
   fillHisto("hData_VtxQualityLoose_HemiAveragept","NoSel",samplename,HemiAveragePt,nRecoLooseVertex);

   fillHisto("hData_VtxQualityRatio_Hemileadingpt","NoSel",samplename,hemi_ptmax,RATIOVTX);
   fillHisto("hData_VtxQualityRatio_Hemisubleadingpt","NoSel",samplename,hemi_ptmin,RATIOVTX);
   fillHisto("hData_VtxQualityRatio_HemiAveragept","NoSel",samplename,HemiAveragePt,RATIOVTX);


   fillHisto("hData_VtxQualityTight_Mmumu","NoSel",samplename,minitree_Mmumu->at(0),nRecoTightVertex);
   fillHisto("hData_VtxQualityLoose_Mmumu","NoSel",samplename,minitree_Mmumu->at(0),nRecoLooseVertex);

   fillHisto("hData_VtxQualityTight_Njet1","NoSel",samplename,minitree_Hemi_njet_nomu->at(0),nVtx1);
   fillHisto("hData_VtxQualityLoose_Njet1","NoSel",samplename,minitree_Hemi_njet_nomu->at(0),nVtxLoose1);

   fillHisto("hData_VtxQualityTight_Njet2","NoSel",samplename,minitree_Hemi_njet_nomu->at(1),nVtx2);
   fillHisto("hData_VtxQualityLoose_Njet2","NoSel",samplename,minitree_Hemi_njet_nomu->at(1),nVtxLoose2);

   fillHisto2D("hData_njetvsHemi1pt","NoSel",samplename,minitree_Hemi_njet_nomu->at(0),hemi1_pt,1);
   fillHisto2D("hData_njetvsHemi2pt","NoSel",samplename,minitree_Hemi_njet_nomu->at(1),hemi2_pt,1);



      if ((nVtx == 1 && nVtxLoose == 0)|| (nVtx == 0 && nVtxLoose == 1) )
      {
          fillHisto2D("hData_Hemi_1Vtx_STW_Ntrks","NoSel",samplename,Vtx_SumtrackWeight,Vtx_nTrks,1);
          fillHisto("hData_Hemi_1Vtx_NChi2","",samplename,Vtx_NChi,1);
          fillHisto("LT_1Vtx","",samplename,minitree_LT->at(0),1);
          fillHisto("LT_1Vtx","Prompt",samplename,LT,1);
          nReco1Vertex++;
          if ((nVtx == 1 && nVtxLoose == 0))
            {
               fillHisto("LT_1Vtx","Tight",samplename,minitree_LT->at(0),1);
               fillHisto("LT_1Vtx","PromptTight",samplename,LT,1);
               nReco1TightVertex++;
            }
          if ((nVtx == 0 && nVtxLoose == 1))
            {
               fillHisto("LT_1Vtx","Loose",samplename,minitree_LT->at(0),1); 
               fillHisto("LT_1Vtx","PromptLoose",samplename,LT,1);
            }
      }

          if ((nVtx == 2 && nVtxLoose == 0)|| (nVtx == 0 && nVtxLoose == 2) || (nVtx == 1 && nVtxLoose == 1))
      {
          fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","NoSel",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,1);
          fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","NoSel",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,1);
          fillHisto("hData_Hemi_2VtxAll_NChi2","",samplename,Vtx_NChi0,1);
          fillHisto("hData_Hemi_2VtxAll_NChi2","",samplename,Vtx_NChi1,1);
          fillHisto("LT_2Vtx","",samplename,minitree_LT->at(0),1);
          fillHisto("LT_2Vtx","Prompt",samplename,LT,1);
          nReco2Vertex++;
          if ((nVtx == 2 && nVtxLoose == 0))
            {
                fillHisto("LT_2Vtx","TT",samplename,minitree_LT->at(0),1);
                fillHisto("LT_2Vtx","PromptTT",samplename,LT,1);
                nReco2TightVertex++;
            }
         if ((nVtx == 0 && nVtxLoose == 2))
            {
               fillHisto("LT_2Vtx","LL",samplename,minitree_LT->at(0),1);
               fillHisto("LT_2Vtx","PromptLL",samplename,LT,1);
            }
         if ((nVtx == 1 && nVtxLoose == 1))
            {
               fillHisto("LT_2Vtx","TL",samplename,minitree_LT->at(0),1);
               fillHisto("LT_2Vtx","PromptTL",samplename,LT,1);
            }
      }

      //-----------------------------------------------------------//
      // ABCD using Hemipshere pt anf Tight+loose steps of vertexing 
      //-----------------------------------------------------------//

      // !! --------SR-------- !! //

    if (minitree_Filter->at(0))// CHANGE
      {
         if ((nVtx == 1 && nVtxLoose == 0) || (nVtx == 0 && nVtxLoose == 1)) 
            {
               int Quality = 0 ;// 0 is tight and 1 is loose
               if (nVtxLoose == 1)  Quality = 1;
               fillHisto2D("hData_ABCD_1Vtx_TL_Hemipt","",samplename, hemi_ptmin, Quality, Vtx_SumtrackWeight);
            }
         if (nVtx == 2  || nVtxLoose == 2) 
            {
               int Quality = 0;
               if (nVtxLoose == 2)  Quality = 1;
               fillHisto2D("hData_ABCD_TLVtxAll_Hemipt","",samplename, hemi_ptmin,Quality, Vtx_SumtrackWeight);
               fillHisto2D("hData_ABCD_TLVtxAll_Hemipt","",samplename, hemi1_pt, Quality, Vtx_SumtrackWeight0);
               fillHisto2D("hData_ABCD_TLVtxAll_Hemipt","",samplename, hemi2_pt, Quality, Vtx_SumtrackWeight1);
            }

      }    


                        // !! -------------- Beginning OF LT ------------  !! //

float LTcut = 80.; //GeV
// can try both LT (prompt muons ) and minitree_LT->at(0) (all muons) 
isCutEvt = false;
if (LT > LTcut){isCutEvt = true;}

                           // ---        1 Vertex        ---//
         if (nVtx == 1 && nVtxLoose == 0 ) 
            {
               if (isCutEvt) // Signal Region C
                  {
                     if (!BlindSR) 
                        {
                           fillHisto("Quality_LT_STW_1Vtx_C","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                        }
                     else
                        {
                           fillHisto("Quality_LT_STW_1Vtx_C","",samplename,0,Pref_PU_gen_wt);
                        }
                     
                  }
               else // Contol region A
                  {
                     fillHisto("Quality_LT_STW_1Vtx_A","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                  }
            }
         if (nVtx == 0 && nVtxLoose == 1 ) 
            {
               if (isCutEvt) // Control Region D
                  {
                     fillHisto("Quality_LT_STW_1Vtx_D","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                  }
               else // Contol region B
                  {
                     fillHisto("Quality_LT_STW_1Vtx_B","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                  } 
            }


                        // ---        2 Vertex        ---//
         if (nVtx == 1 && nVtxLoose == 1 ) 
            {
               if (isCutEvt) // Signal Region C
                  {
                     if (!BlindSR) 
                        {
                           fillHisto("Quality_LT_STW_2Vtx_C","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                           fillHisto("Quality_LT_STW_2VtxAll_C","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                           fillHisto("Quality_LT_STW_2VtxAll_C","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                        }
                     else
                        {
                           fillHisto("Quality_LT_STW_2Vtx_C","",samplename,0,Pref_PU_gen_wt);
                           fillHisto("Quality_LT_STW_2VtxAll_C","",samplename,0,Pref_PU_gen_wt);
                           fillHisto("Quality_LT_STW_2VtxAll_C","",samplename,0,Pref_PU_gen_wt);
                        }

                  }
               else // Contol region A
                  {
                     fillHisto("Quality_LT_STW_2Vtx_A","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                     fillHisto("Quality_LT_STW_2VtxAll_A","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                     fillHisto("Quality_LT_STW_2VtxAll_A","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  } 
            }
         if (nVtx == 2 && nVtxLoose == 0 ) 
            {
               if (isCutEvt) // Signal Region C
                  {
                     if (!BlindSR) 
                        {
                           fillHisto("Quality_LT_STW_2Vtx_C","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                           fillHisto("Quality_LT_STW_2VtxAll_C","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                           fillHisto("Quality_LT_STW_2VtxAll_C","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                        }
                     else
                        {
                           fillHisto("Quality_LT_STW_2Vtx_C","",samplename,0,Pref_PU_gen_wt);
                           fillHisto("Quality_LT_STW_2VtxAll_C","",samplename,0,Pref_PU_gen_wt);
                           fillHisto("Quality_LT_STW_2VtxAll_C","",samplename,0,Pref_PU_gen_wt);
                        }

                  }
               else // Contol region A
                  {
                     fillHisto("Quality_LT_STW_2Vtx_A","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                     fillHisto("Quality_LT_STW_2VtxAll_A","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                     fillHisto("Quality_LT_STW_2VtxAll_A","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  } 
            }
         if (nVtx == 0 && nVtxLoose == 2 ) 
            {
               if (isCutEvt) // Control Region D
                  {
                     fillHisto("Quality_LT_STW_2Vtx_D","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                     fillHisto("Quality_LT_STW_2VtxAll_D","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                     fillHisto("Quality_LT_STW_2VtxAll_D","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  }
               else // Contol region B
                  {
                     fillHisto("Quality_LT_STW_2Vtx_B","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                     fillHisto("Quality_LT_STW_2VtxAll_B","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                     fillHisto("Quality_LT_STW_2VtxAll_B","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  } 
            }

                        // !! -------------- END OF LT ------------  !! //


   isCutEvt = false;
    if (!BlindSR) {
         
         // if (showoutput) // std::cout<<"SR "<<std::endl;
         fillHisto("hData_Hemi_BDTevt","",samplename, minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);
         //$$
   
         if (hemi_ptmin > 80.) isCutEvt = true;
            //    //$$
            //          hemi1_pt = minitree_Hemi_pt->at(0);
            // hemi2_pt = minitree_Hemi_pt->at(1);

         if (nVtx == 0 && nVtxLoose == 0 && isCutEvt) {
            fillHisto("hData_Hemi_0Vtx_Mmumu","",samplename,minitree_Mmumu->at(0),Pref_PU_gen_wt);
            fillHisto("hData_Hemi_0Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
         }

         
         if (nVtx == 1 && nVtxLoose==0 && isCutEvt ) {
            fillHisto("hData_Hemi_1Vtx_Mmumu","",samplename,minitree_Mmumu->at(0),Pref_PU_gen_wt);
            fillHisto("hData_Hemi_1Vtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_1Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_1Vtx_BDTvtx","",samplename, BDTvtx,Pref_PU_gen_wt );
            fillHisto("hData_Hemi_1Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight,Pref_PU_gen_wt );
            fillHisto("hData_Hemi_1Vtx_dist","",samplename, Vtx_dist,Pref_PU_gen_wt );
            fillHisto2D("hData_Hemi_1Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);

         }


      if (nVtx == 1 && nVtxLoose == 1 && isCutEvt) { // same as Loose since TL = LT
         fillHisto("hData_Hemi_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
         fillHisto("hData_Hemi_TLVtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);

         fillHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
         fillHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
         fillHisto("hData_Hemi_TLVtxAll_Mass","",samplename,Vtx_Mass0,Pref_PU_gen_wt);
         fillHisto("hData_Hemi_TLVtxAll_Mass","",samplename,Vtx_Mass1,Pref_PU_gen_wt);
         fillHisto2D("hData_Hemi_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
         fillHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
         fillHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

      }
         
         if (nVtx == 2 && isCutEvt ) {
            fillHisto("hData_Hemi_2Vtx_Mmumu","",samplename,minitree_Mmumu->at(0),Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2Vtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
            fillHisto("hData_Hemi_2Vtx_MaxBDTvtx","",samplename, BDTvtx,Pref_PU_gen_wt );
            fillHisto("hData_Hemi_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight,Pref_PU_gen_wt );
            fillHisto2D("hData_Hemi_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);

            fillHisto("hData_Hemi_2VtxAll_Mass","",samplename,Vtx_Mass0,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2VtxAll_Mass","",samplename,Vtx_Mass1,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2VtxAll_BDTvtx","",samplename, BDTvtx1 ,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0,Pref_PU_gen_wt );
            fillHisto("hData_Hemi_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1,Pref_PU_gen_wt );
            fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
            fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);


         }

      }
   else {
         
         // if (showoutput) // std::cout<<"SR "<<std::endl;
         fillHisto("hData_Hemi_BDTevt","",samplename, 0 ,Pref_PU_gen_wt);
         //$$
   
         if (hemi_ptmin > 80.) isCutEvt = true;
            //    //$$
            //          hemi1_pt = minitree_Hemi_pt->at(0);
            // hemi2_pt = minitree_Hemi_pt->at(1);

         if (nVtx == 0 && nVtxLoose == 0 && isCutEvt) {
            fillHisto("hData_Hemi_0Vtx_Mmumu","",samplename,0,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_0Vtx_BDTevt","",samplename, 0,Pref_PU_gen_wt );
         }

         
         if (nVtx == 1 && nVtxLoose==0 && isCutEvt ) {
            fillHisto("hData_Hemi_1Vtx_Mmumu","",samplename,0,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_1Vtx_Mass","",samplename,0,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_1Vtx_BDTevt","",samplename, 0 ,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_1Vtx_BDTvtx","",samplename, 0,Pref_PU_gen_wt );
            fillHisto("hData_Hemi_1Vtx_SumtrackWeight","",samplename, 0,Pref_PU_gen_wt );
            fillHisto("hData_Hemi_1Vtx_dist","",samplename, 0,Pref_PU_gen_wt );
            fillHisto2D("hData_Hemi_1Vtx_STW_Ntrks","",samplename,0,0,Pref_PU_gen_wt);

         }


      if (nVtx == 1 && nVtxLoose == 1 && isCutEvt ) { // same as Loose since TL = LT
         fillHisto("hData_Hemi_TLVtx_SumtrackWeight","",samplename,0,Pref_PU_gen_wt);
         fillHisto("hData_Hemi_TLVtx_Mass","",samplename,0,Pref_PU_gen_wt);

         fillHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",samplename,0,Pref_PU_gen_wt);
         fillHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",samplename,0,Pref_PU_gen_wt);
         fillHisto("hData_Hemi_TLVtxAll_Mass","",samplename,0,Pref_PU_gen_wt);
         fillHisto("hData_Hemi_TLVtxAll_Mass","",samplename,0,Pref_PU_gen_wt);
         fillHisto2D("hData_Hemi_TLVtx_STW_Ntrks","",samplename,0,0,Pref_PU_gen_wt);
         fillHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",samplename,0,0,Pref_PU_gen_wt);
         fillHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",samplename,0,0,Pref_PU_gen_wt);

      }
         
         if (nVtx == 2 && isCutEvt ) {
            fillHisto("hData_Hemi_2Vtx_Mmumu","",samplename,minitree_Mmumu->at(0),Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2Vtx_Mass","",samplename,0,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2Vtx_BDTevt","",samplename, 0,Pref_PU_gen_wt );
            fillHisto("hData_Hemi_2Vtx_MaxBDTvtx","",samplename, 0,Pref_PU_gen_wt );
            fillHisto("hData_Hemi_2Vtx_SumtrackWeight","",samplename, 0,Pref_PU_gen_wt );
            fillHisto2D("hData_Hemi_2Vtx_STW_Ntrks","",samplename,0,Vtx_nTrks,Pref_PU_gen_wt);

            fillHisto("hData_Hemi_2VtxAll_Mass","",samplename,0,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2VtxAll_Mass","",samplename,0,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2VtxAll_BDTvtx","",samplename, 0 ,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,Pref_PU_gen_wt);
            fillHisto("hData_Hemi_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0,Pref_PU_gen_wt );
            fillHisto("hData_Hemi_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1,Pref_PU_gen_wt );
            fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
            fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);


         }

      }

      //--------CR : Low Pt --------//

      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      // hemi1_pt = minitree_Hemi_pt->at(0);
      // hemi2_pt = minitree_Hemi_pt->at(1);
      if (Filter &&  ( (hemi1_pt >= lowHpt && hemi1_pt < 80. && hemi2_pt >= 80.) ||
             (hemi2_pt >= lowHpt && hemi2_pt < 80. && hemi1_pt >= 80.) ))//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {
             //$$
            fillHisto("hData_CRlowpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );

            // if ( hemi_ptmin > 80. ) isCutEvt = true;   
            //$$
            if (nVtx == 0   && nVtxLoose == 0 )
               {
                  fillHisto("hData_CRlowpt_0Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlowpt_0Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  
               }
            if (nVtx == 1  && nVtxLoose==0 )
               {
                  fillHisto("hData_CRlowpt_1Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowpt_1Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );
                  fillHisto("hData_CRlowpt_1Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlowpt_1Vtx_BDTvtx","",samplename, BDTvtx,Pref_PU_gen_wt );
                  fillHisto("hData_CRlowpt_1Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  fillHisto2D("hData_CRlowpt_1Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
               }


            if (nVtx == 1 && nVtxLoose == 1  ) { // same as Loose low pt since TL = LT
                 fillHisto("hData_CRlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                 fillHisto("hData_CRlowpt_TLVtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);

                 fillHisto("hData_CRlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                 fillHisto("hData_CRlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                 fillHisto("hData_CRlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,Pref_PU_gen_wt);
                 fillHisto("hData_CRlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,Pref_PU_gen_wt);
                 fillHisto2D("hData_CRlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                 fillHisto2D("hData_CRlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                 fillHisto2D("hData_CRlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

            }

            if (nVtx == 2 )
               {
                  fillHisto("hData_CRlowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowpt_2Vtx_Mass","",samplename,   VtxMass ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowpt_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlowpt_2Vtx_MaxBDTvtx","",samplename, BDTvtx ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlowpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  
                  fillHisto("hData_CRlowpt_2VtxAll_Mass","",samplename, Vtx_Mass0,Pref_PU_gen_wt );
                  fillHisto("hData_CRlowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,Pref_PU_gen_wt );
                  fillHisto("hData_CRlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,Pref_PU_gen_wt);         
                  fillHisto2D("hData_CRlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);


               }

         }

         //--------CR : Both vertices have Low Pt --------//
 
      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      // hemi1_pt = minitree_Hemi_pt->at(0);
      // hemi2_pt = minitree_Hemi_pt->at(1);
      if (Filter &&   hemi1_pt >= lowHpt && hemi1_pt < 80. && lowHpt < hemi2_pt < 80. )//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {
             //$$

            fillHisto("hData_CRlowlowpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );

            // if ( hemi_ptmin > 80. ) isCutEvt = true;   
            //$$
            if (nVtx == 0   && nVtxLoose == 0 )
               {
                  fillHisto("hData_CRlowlowpt_0Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlowlowpt_0Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  
               }
            if (nVtx == 1  && nVtxLoose==0 )
               {
                  fillHisto("hData_CRlowlowpt_1Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_1Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );
                  fillHisto("hData_CRlowlowpt_1Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlowlowpt_1Vtx_BDTvtx","",samplename, BDTvtx,Pref_PU_gen_wt );

                  fillHisto("hData_CRlowlowpt_1Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  fillHisto2D("hData_CRlowlowpt_1Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
               }


            if (nVtx == 1 && nVtxLoose == 1  ) { // same as LooseLowLowpt since TL = LT
                  fillHisto("hData_CRlowlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_TLVtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlowlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);

                  fillHisto("hData_CRlowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

            }

            if (nVtx == 2 )
               {
                  fillHisto("hData_CRlowlowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlowlowpt_2Vtx_Mass","",samplename,   VtxMass ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlowlowpt_2Vtx_MaxBDTvtx","",samplename, BDTvtx ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_2VtxAll_Mass","",samplename, Vtx_Mass0,Pref_PU_gen_wt );
                  fillHisto("hData_CRlowlowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,Pref_PU_gen_wt );
                  fillHisto("hData_CRlowlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlowlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,Pref_PU_gen_wt);

                  fillHisto2D("hData_CRlowlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlowlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);


               }

         }
      //--------CR : Loose  --------//

      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if ( Filter ) 
         {
            fillHisto("hData_CRloose_BDTevt","",samplename, minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);
            //$$
            if ( hemi_ptmin > 80. ) isCutEvt = true;
            //$$ 
            if (nVtxLoose == 0  && nVtx == 0 && isCutEvt)
               {
                  fillHisto("hData_CRloose_0Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRloose_0Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0) ,Pref_PU_gen_wt);
                  
               }
            if (nVtxLoose == 1 && nVtx==0 && isCutEvt)
               {
                  fillHisto("hData_CRloose_1Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_1Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );
                  fillHisto("hData_CRloose_1Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRloose_1Vtx_BDTvtx","",samplename, BDTvtx,Pref_PU_gen_wt );
                  fillHisto("hData_CRloose_1Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  fillHisto2D("hData_CRloose_1Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);

               }

            if (nVtx == 1 && nVtxLoose == 1 && isCutEvt ) { 
                  fillHisto("hData_CRloose_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_TLVtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRloose_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_TLVtxAll_Mass","",samplename,Vtx_Mass0,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_TLVtxAll_Mass","",samplename,Vtx_Mass1,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRloose_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRloose_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);
            }

            if (nVtxLoose == 2 && isCutEvt)
               {
                  fillHisto("hData_CRloose_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0) ,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_2Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );
                  fillHisto("hData_CRloose_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRloose_2Vtx_MaxBDTvtx","",samplename, BDTvtx,Pref_PU_gen_wt );
                  fillHisto("hData_CRloose_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  fillHisto2D("hData_CRloose_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRloose_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRloose_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_2VtxAll_Mass","",samplename, Vtx_Mass0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_2VtxAll_Mass","",samplename, Vtx_Mass1 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRloose_2VtxAll_BDTvtx","",samplename, BDTvtx1,Pref_PU_gen_wt );
                  fillHisto("hData_CRloose_2VtxAll_BDTvtx","",samplename, BDTvtx2,Pref_PU_gen_wt );           
                  fillHisto("hData_CRloose_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0,Pref_PU_gen_wt );
                  fillHisto("hData_CRloose_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1,Pref_PU_gen_wt );

               }  
         }

      //--------CR : Loose Lowpt  --------//

      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if (Filter && ( (hemi1_pt >= lowHpt && hemi1_pt < 80. && hemi2_pt >= 80.) ||
             (hemi2_pt >= lowHpt && hemi2_pt < 80. && hemi1_pt >= 80.) )  )//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {
            fillHisto("hData_CRlooselowpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );

            if (nVtxLoose == 0  && nVtx == 0)
               {
                  fillHisto("hData_CRlooselowpt_0Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowpt_1Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowpt_1Vtx_BDTvtx","",samplename, BDTvtx ,Pref_PU_gen_wt);
               }
            if (nVtxLoose == 1 && nVtx==0 )
               {
                  fillHisto("hData_CRlooselowpt_1Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowpt_1Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowpt_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowpt_2Vtx_MaxBDTvtx","",samplename, BDTvtx,Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowpt_1Vtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  fillHisto2D("hData_CRlooselowpt_1Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
               }

            if (nVtx == 1 && nVtxLoose == 1  ) {
                  fillHisto("hData_CRlooselowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowpt_TLVtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);

                  fillHisto("hData_CRlooselowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

            }
            if (nVtxLoose == 2)
               {
                  fillHisto("hData_CRlooselowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowpt_2Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );

                  fillHisto("hData_CRlooselowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,Pref_PU_gen_wt);

                  fillHisto2D("hData_CRlooselowpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

                  fillHisto("hData_CRlooselowpt_2VtxAll_Mass","",samplename, Vtx_Mass0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,Pref_PU_gen_wt);

                  fillHisto("hData_CRlooselowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,Pref_PU_gen_wt);

               }

         }

               //--------CR : LooseLoose  LowLowpt  --------//

      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if (Filter &&  hemi1_pt >= lowHpt && hemi1_pt < 80. &&  lowHpt <= hemi2_pt  &&  hemi2_pt <= 80.   )//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {
            fillHisto("hData_CRlooselowlowpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
            if (nVtxLoose == 0  && nVtx == 0)
               {
                  fillHisto("hData_CRlooselowlowpt_0Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowlowpt_1Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowlowpt_1Vtx_BDTvtx","",samplename, BDTvtx ,Pref_PU_gen_wt);
                  
               }
            if (nVtxLoose == 1 && nVtx==0 )
               {
                  fillHisto("hData_CRlooselowlowpt_1Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowlowpt_1Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowlowpt_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowlowpt_2Vtx_MaxBDTvtx","",samplename, BDTvtx,Pref_PU_gen_wt );

                  fillHisto("hData_CRlooselowlowpt_1Vtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt );
                  fillHisto2D("hData_CRlooselowlowpt_1Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);

               }

            if (nVtx == 1 && nVtxLoose == 1  ) {
                  fillHisto("hData_CRlooselowlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowlowpt_TLVtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);


                  fillHisto2D("hData_CRlooselowlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

                  fillHisto("hData_CRlooselowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,Pref_PU_gen_wt);
            }
            if (nVtxLoose == 2)
               {
                  fillHisto("hData_CRlooselowlowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowlowpt_2Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );

                  fillHisto("hData_CRlooselowlowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,Pref_PU_gen_wt);

                  fillHisto2D("hData_CRlooselowlowpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselowlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselowlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

                  fillHisto("hData_CRlooselowlowpt_2VtxAll_Mass","",samplename, Vtx_Mass0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowlowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselowlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,Pref_PU_gen_wt);

                  fillHisto("hData_CRlooselowlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselowlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,Pref_PU_gen_wt);


               }

         }


         // !! rename la première catégorie !! //
         // Histograms for 2Vertices but TT and TL are combined together as well as HighLow and LowLow
         // Gives more stats in one box for the signal region
         // ends up in a new ABCD method and not ABCDEFGHI


      // this is new B, same as old C with LooseLooseLowLowpt
      if (Filter &&  hemi1_pt >= lowHpt && hemi1_pt < 80. &&  lowHpt <= hemi2_pt  &&  hemi2_pt <= 80.   )//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {
            fillHisto("hData_CRlooselooselowlowpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
            if ( nVtx == 0 && nVtxLoose == 2 ) {
                  fillHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);


                  fillHisto2D("hData_CRlooselooselowlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselooselowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselooselowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,Pref_PU_gen_wt);
            }
         }


//--------CR : TightLoose + ITghtTight for low lowpt --------//
// Old A and B combined to form new A
      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if (Filter &&  hemi1_pt >= lowHpt && hemi1_pt < 80. &&  lowHpt <= hemi2_pt  &&  hemi2_pt <= 80. )//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {
            fillHisto("hData_CRtightlowlowpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
            if ( (nVtx == 1 && nVtxLoose == 1 ) || (nVtx == 2 && nVtxLoose == 0) ){
                  fillHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt);
                  fillHisto("hData_CRtightlowlowpt_TLVtx_Mass","",samplename,VtxMass,Pref_PU_gen_wt);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,Pref_PU_gen_wt);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,Pref_PU_gen_wt);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,Pref_PU_gen_wt);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRtightlowlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRtightlowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRtightlowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);
            }
         }


//--------CR : LooseLoose for highhighpt or highlowpt --------//
// Old F and I combined to form new D


      if (Filter &&  ((hemi1_pt >= 80 && hemi2_pt >= 80) || ( (hemi1_pt >= lowHpt && hemi1_pt < 80. && hemi2_pt >= 80.) ||
             (hemi2_pt >= lowHpt && hemi2_pt < 80. && hemi1_pt >= 80.) )  ))//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {
            fillHisto("hData_CRlooselooselowpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );
            if (nVtxLoose == 1 && nVtx == 0)
               {
                  fillHisto("hData_CRLoose_1Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );
                  fillHisto("hData_CRLoose_1Vtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt );
               }
            if (nVtxLoose == 0 && nVtx == 1)
               {
                  fillHisto("hData_CRTight_1Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );
                  fillHisto("hData_CRTight_1Vtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,Pref_PU_gen_wt );
               }

            if (nVtxLoose == 2 && nVtx == 0)
               {
                  fillHisto("hData_CRlooselooselowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );

                  fillHisto("hData_CRlooselooselowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,Pref_PU_gen_wt);

                  fillHisto2D("hData_CRlooselooselowpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselooselowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRlooselooselowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_Mass","",samplename, Vtx_Mass0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,Pref_PU_gen_wt );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,Pref_PU_gen_wt);

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,Pref_PU_gen_wt);
               }

         }

      //--------SR : TIghtLoose+ TIghtTight for highhighpt or highlowpt --------//
// Old F and I combined to form new C


      if (Filter &&  ((hemi1_pt >= 80 && hemi2_pt >= 80) || ( (hemi1_pt >= lowHpt && hemi1_pt < 80. && hemi2_pt >= 80.) ||
             (hemi2_pt >= lowHpt && hemi2_pt < 80. && hemi1_pt >= 80.) )  ))//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {
            fillHisto("hData_CRtighthighpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),Pref_PU_gen_wt );

            if (nVtx == 2 || (nVtx == 1 && nVtxLoose == 1 ) ){
   
                  fillHisto("hData_CRtighthighpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),Pref_PU_gen_wt );
                  fillHisto("hData_CRtighthighpt_2Vtx_Mass","",samplename,   VtxMass,Pref_PU_gen_wt );

                  fillHisto("hData_CRtighthighpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,Pref_PU_gen_wt);

                  fillHisto2D("hData_CRtighthighpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRtighthighpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,Pref_PU_gen_wt);
                  fillHisto2D("hData_CRtighthighpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,Pref_PU_gen_wt);

                  fillHisto("hData_CRtighthighpt_2VtxAll_Mass","",samplename, Vtx_Mass0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRtighthighpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRtighthighpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,Pref_PU_gen_wt );
                  fillHisto("hData_CRtighthighpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,Pref_PU_gen_wt);

                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,Pref_PU_gen_wt);
                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,Pref_PU_gen_wt);
               }

         }

         // std::cout<<" here 4 "<<std::endl;
//----------------END of njet caeroisation ---------------------------//


   }// end of event loop

   //add end accolade
      NormFactor = NormFactor/allevents;
   fillHisto("hData_StepEff","",samplename,0.,allevents*NormFactor);
   fillHisto("hData_StepEff","",samplename,1.,nFilterEvt*NormFactor);//
   fillHisto("hData_StepEff","",samplename,2.,nFilterJet*NormFactor);
   fillHisto("hData_StepEff","",samplename,3.,nFilternHemi*NormFactor);
   fillHisto("hData_StepEff","",samplename,4.,nFilterHpt*NormFactor);
   fillHisto("hData_StepEff","",samplename,5.,nReco1Vertex*NormFactor);
   fillHisto("hData_StepEff","",samplename,6.,nReco2Vertex*NormFactor/2.);
   fillHisto("hData_StepEff","",samplename,7.,nReco1TightVertex*NormFactor);
   fillHisto("hData_StepEff","",samplename,8.,nReco2TightVertex*NormFactor/2.);

   fillHisto("hData_StepEff","NonNormalized",samplename,0.,allevents);
   fillHisto("hData_StepEff","NonNormalized",samplename,1.,nFilterEvt);//
   fillHisto("hData_StepEff","NonNormalized",samplename,2.,nFilterJet);
   fillHisto("hData_StepEff","NonNormalized",samplename,3.,nFilternHemi);
   fillHisto("hData_StepEff","NonNormalized",samplename,4.,nFilterHpt);
   fillHisto("hData_StepEff","NonNormalized",samplename,5.,nReco1Vertex);
   fillHisto("hData_StepEff","NonNormalized",samplename,6.,nReco2Vertex);
   fillHisto("hData_StepEff","NonNormalized",samplename,7.,nReco1TightVertex);
   fillHisto("hData_StepEff","NonNormalized",samplename,8.,nReco2TightVertex);

       theoutputfile->Write();
       // std::cout<<"File has been written: "<<theoutputfile->GetName()<<std::endl;

   // std::cout<<"1"<<std::endl;
   fL1_Reco_SF->Close();
      // // std::cout<<"2"<<std::endl;
// std::cout<<"2"<<std::endl;
   fL1_ID_SF->Close();
      // // std::cout<<"3"<<std::endl;
// std::cout<<"3"<<std::endl;
   fL1_ISO_SF->Close();
      // // std::cout<<"4"<<std::endl;
// std::cout<<"4"<<std::endl;
   fL2_ISO_SF->Close();
      // // std::cout<<"5"<<std::endl;
// std::cout<<"5"<<std::endl;
   fL2_Reco_SF->Close();
      // // std::cout<<"6"<<std::endl;
// std::cout<<"6"<<std::endl;
   // fL2_Reco_SF2->Close();
      // // std::cout<<"7"<<std::endl;

   fL2_ID_SF->Close();
      // // std::cout<<"8"<<std::endl;
// std::cout<<"7"<<std::endl;
   fL1L2_TRG_SF->Close();
      // // std::cout<<"9"<<std::endl;
//  std::cout<<"8"<<std::endl;
   //deleteHisto();
   theoutputfile->Close();
      // // std::cout<<"10"<<std::endl;

   // delete theoutputfile;


}// end of systematics
}// end of loop method
void TreeABCDReader::initializeHisto(TString sample, bool isfirstset){


  cout << "#####################################" << endl;
  cout << "#####################################" << endl;
  cout << " initialize histograms of sample :  " <<sample<< endl;
  cout << "#####################################" << endl;
  cout << "#####################################" << endl;
  

   if(isfirstset){
      numb_histo = 0;
      numb_histo_2D_ = 0;
      TH1F * first_emptyHisto = new TH1F("first_emptyHisto", "first_emptyHisto", 100, 0, 1000);
      TH2F * first_emptyHisto_2D = new TH2F("first_emptyHisto_2D", "first_emptyHisto_2D", 100, 0, 1000,100,0,1000);
      histo_list_.push_back(first_emptyHisto);
      histo_list_2D_.push_back(first_emptyHisto_2D);
      numb_histo++;
      numb_histo_2D_++;
   }

//-----------AddHisto ------//
   addHisto("GenWeight","", sample.Data(), 300,0,3);
   addHisto("hData_Event_Weight","",sample.Data(), 1000,0,1000 );
   addHisto("hData_Filter","",sample.Data(),2,-0.5,1.5 );
   addHisto("hData_Mmumu","",sample.Data(), 150, 0 ,1500);

   addHisto2D("hData_TriggerEff2D","pass",sample.Data(), 20, 0 ,100,20, 0 ,100);
   addHisto2D("hData_TriggerEff2D","fail",sample.Data(), 20, 0 ,100,20, 0 ,100);
   addHisto("hData_TriggerEff","pass",sample.Data(), 50, 0 ,100);
   addHisto("hData_TriggerEff","fail",sample.Data(), 50, 0 ,100);

   addHisto2D("hData_2TriggerEff2D","pass",sample.Data(), 20, 0 ,100,20, 0 ,100);
   addHisto2D("hData_2TriggerEff2D","fail",sample.Data(), 20, 0 ,100,20, 0 ,100);
   addHisto("hData_2TriggerEff","pass",sample.Data(), 50, 0 ,100);
   addHisto("hData_2TriggerEff","fail",sample.Data(), 50, 0 ,100);
   //------- LT --------//

   addHisto("Quality_LT_STW_1Vtx_C","",sample.Data(),40, 0, 40);
   addHisto("Quality_LT_STW_1Vtx_A","",sample.Data(),40, 0, 40);
   addHisto("Quality_LT_STW_1Vtx_D","",sample.Data(),40, 0, 40);
   addHisto("Quality_LT_STW_1Vtx_B","",sample.Data(),40, 0, 40);

   addHisto("Quality_LT_STW_2Vtx_C","",sample.Data(),40, 0, 40);
   addHisto("Quality_LT_STW_2VtxAll_C","",sample.Data(),40, 0, 40);


   addHisto("Quality_LT_STW_2Vtx_A","",sample.Data(),40, 0, 40);
   addHisto("Quality_LT_STW_2VtxAll_A","",sample.Data(),40, 0, 40);


   addHisto("Quality_LT_STW_2Vtx_D","",sample.Data(),40, 0, 40);
   addHisto("Quality_LT_STW_2VtxAll_D","",sample.Data(),40, 0, 40);


   addHisto("Quality_LT_STW_2Vtx_B","",sample.Data(),40, 0, 40);
   addHisto("Quality_LT_STW_2VtxAll_B","",sample.Data(),40, 0, 40);


   //---------Tight/Loose------//

addHisto("hData_VtxQualityTight_Hemi1pt","2Vtx",sample.Data(),50,0,1000);
addHisto("hData_VtxQualityTight_Hemi2pt","2Vtx",sample.Data(),50,0,1000);
addHisto("hData_VtxQualityTight_Hemi1pt","1Vtx",sample.Data(),50,0,1000);
addHisto("hData_VtxQualityTight_Hemi2pt","1Vtx",sample.Data(),50,0,1000);


addHisto("hData_VtxQualityLoose_Hemi1pt","1Vtx",sample.Data(),50,0,1000);
addHisto("hData_VtxQualityLoose_Hemi2pt","1Vtx",sample.Data(),50,0,1000);
addHisto("hData_VtxQualityLoose_Hemi1pt","2Vtx",sample.Data(),50,0,1000);
addHisto("hData_VtxQualityLoose_Hemi2pt","2Vtx",sample.Data(),50,0,1000);

   addHisto2D("hData_njetvsHemi1pt","NoSel",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi2pt","NoSel",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi1pt","1jet",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi2pt","1jet",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi1pt","2jet",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi2pt","2jet",sample.Data(),20, 0 ,20,50,0,1000);
      addHisto2D("hData_njetvsHemi1pt","3jet",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi2pt","3jet",sample.Data(),20, 0 ,20,50,0,1000);
      addHisto2D("hData_njetvsHemi1pt","4jet",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi2pt","4jet",sample.Data(),20, 0 ,20,50,0,1000);
      addHisto2D("hData_njetvsHemi1pt","5jet",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi2pt","5jet",sample.Data(),20, 0 ,20,50,0,1000);
      addHisto2D("hData_njetvsHemi1pt","6jet",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi2pt","6jet",sample.Data(),20, 0 ,20,50,0,1000);
         addHisto2D("hData_njetvsHemi1pt","7jet",sample.Data(),20, 0 ,20,50,0,1000);
   addHisto2D("hData_njetvsHemi2pt","7jet",sample.Data(),20, 0 ,20,50,0,1000);
   //------------------//
   addHisto("hData_VtxQualityRatio_Hemileadingpt","NoSel",sample.Data(), 50,0,1000);
   addHisto("hData_VtxQualityRatio_Hemisubleadingpt","NoSel",sample.Data(), 50,0,1000);
   addHisto("hData_VtxQualityRatio_HemiAveragept","NoSel",sample.Data(), 50,0,1000);

   addHisto("hData_VtxQualityTight_Hemi1pt","NoSel",sample.Data(),50,0,1000);
   addHisto("hData_VtxQualityTight_Hemi2pt","NoSel",sample.Data(),50,0,1000);
   addHisto("hData_VtxQualityLoose_Hemi1pt","NoSel",sample.Data(),50,0,1000);
   addHisto("hData_VtxQualityLoose_Hemi2pt","NoSel",sample.Data(),50,0,1000);

   addHisto("hData_VtxQualityTight_Hemileadingpt","NoSel",sample.Data(), 50,0,1000);
   addHisto("hData_VtxQualityTight_Hemisubleadingpt","NoSel",sample.Data(), 50,0,1000);

   addHisto("hData_VtxQualityLoose_Hemileadingpt","NoSel",sample.Data(), 50,0,1000);
   addHisto("hData_VtxQualityLoose_Hemisubleadingpt","NoSel",sample.Data(), 50,0,1000);

   addHisto("hData_VtxQualityTight_HemiAveragept","NoSel",sample.Data(), 50,0,1000);
   addHisto("hData_VtxQualityLoose_HemiAveragept","NoSel",sample.Data(), 50,0,1000);


   addHisto("hData_VtxQualityTight_Mmumu","NoSel",sample.Data(),100, 0 ,200);
   addHisto("hData_VtxQualityLoose_Mmumu","NoSel",sample.Data(),100, 0 ,200);

   addHisto("hData_VtxQualityTight_Njet1","NoSel",sample.Data(),20, 0 ,20);
   addHisto("hData_VtxQualityLoose_Njet1","NoSel",sample.Data(),20, 0 ,20);

   addHisto("hData_VtxQualityTight_Njet2","NoSel",sample.Data(),20, 0 ,20);
   addHisto("hData_VtxQualityLoose_Njet2","NoSel",sample.Data(),20, 0 ,20);

   
 //-----------------//
   addHisto("LT_1Vtx","",sample.Data(),50,0,1000);
   addHisto("LT_1Vtx","Tight",sample.Data(),50,0,1000);
   addHisto("LT_1Vtx","Loose",sample.Data(),50,0,1000); 
   addHisto("LT_2Vtx","",sample.Data(),50,0,1000);
   addHisto("LT_2Vtx","TT",sample.Data(),50,0,1000);
   addHisto("LT_2Vtx","LL",sample.Data(),50,0,1000);
   addHisto("LT_2Vtx","TL",sample.Data(),50,0,1000);
   addHisto("LT_1Vtx","Prompt",sample.Data(),50,0,1000);
   addHisto("LT_1Vtx","PromptTight",sample.Data(),50,0,1000);
   addHisto("LT_1Vtx","PromptLoose",sample.Data(),50,0,1000); 
   addHisto("LT_2Vtx","Prompt",sample.Data(),50,0,1000);
   addHisto("LT_2Vtx","PromptTT",sample.Data(),50,0,1000);
   addHisto("LT_2Vtx","PromptLL",sample.Data(),50,0,1000);
   addHisto("LT_2Vtx","PromptTL",sample.Data(),50,0,1000);
   //-----------------//
   addHisto2D("hData_Hemi_1Vtx_STW_Ntrks","NoSel",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_Hemi_1Vtx_NChi2","",sample.Data(),64,-1,15);
   addHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","NoSel",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_Hemi_2VtxAll_NChi2","",sample.Data(),64,-1,15);

   addHisto2D("hData_ABCD_1Vtx_TL_Hemipt","",sample.Data(), 300, 0, 500, 2, 0, 2);
   addHisto2D("hData_ABCD_TLVtxAll_Hemipt","",sample.Data(), 300, 0, 500, 2, 0, 2);

   addHisto("Hemisphere_leadingpt","", sample.Data(), 50,0,1000);
   addHisto("Hemisphere_subleadingpt","", sample.Data(), 50,0,1000);

   addHisto("hData_Hemi_BDTevt","",sample.Data(), 25,-1,1);
   addHisto("hData_Hemi_0Vtx_Mmumu","",sample.Data(),25,0.,500.);
   addHisto("hData_Hemi_0Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_Hemi_1Vtx_Mmumu","",sample.Data(),25,0.,500.);
   addHisto("hData_Hemi_1Vtx_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_Hemi_1Vtx_BDTevt","",sample.Data(), 25,-1,1);
   addHisto("hData_Hemi_1Vtx_BDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_Hemi_1Vtx_SumtrackWeight","",sample.Data(), 40, 0, 40 );
   addHisto("hData_Hemi_1Vtx_dist","",sample.Data(), 25,-1,1 );
   addHisto2D("hData_Hemi_1Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_Hemi_TLVtx_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_Hemi_TLVtx_Mass","",sample.Data(),25,0.,100.);

   addHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",sample.Data(),40, 0, 40);

   addHisto("hData_Hemi_TLVtxAll_Mass","",sample.Data(),25,0.,100.);

   addHisto2D("hData_Hemi_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_Hemi_2Vtx_Mmumu","",sample.Data(),25,0.,500.);
   addHisto("hData_Hemi_2Vtx_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_Hemi_2Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_Hemi_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_Hemi_2Vtx_SumtrackWeight","",sample.Data(), 40, 0, 40 );
   addHisto2D("hData_Hemi_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);

   addHisto("hData_Hemi_2VtxAll_Mass","",sample.Data(),25,0.,100.);

   addHisto("hData_Hemi_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1);

   addHisto("hData_Hemi_2VtxAll_SumtrackWeight","",sample.Data(),40, 0, 40 );

   addHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);


   addHisto("hData_CRlowpt_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowpt_0Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlowpt_0Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowpt_1Vtx_Mmumu","",sample.Data(),  25,0.,500.);
   addHisto("hData_CRlowpt_1Vtx_Mass","",sample.Data(),   25,0.,100. );
   addHisto("hData_CRlowpt_1Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowpt_1Vtx_BDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowpt_1Vtx_SumtrackWeight","",sample.Data(),40, 0, 40 );
   addHisto2D("hData_CRlowpt_1Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlowpt_TLVtx_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_CRlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);

   addHisto("hData_CRlowpt_TLVtxAll_SumtrackWeight","",sample.Data(),40, 0, 40);

   addHisto("hData_CRlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);

   addHisto2D("hData_CRlowpt_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlowpt_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlowpt_2Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlowpt_2Vtx_Mass","",sample.Data(),   25,0.,100.);
   addHisto("hData_CRlowpt_2Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowpt_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1);
   addHisto("hData_CRlowpt_2Vtx_SumtrackWeight","",sample.Data(), 40, 0, 40);
   addHisto2D("hData_CRlowpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);

   addHisto("hData_CRlowpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);

   addHisto("hData_CRlowpt_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1 );

   addHisto("hData_CRlowpt_2VtxAll_SumtrackWeight","",sample.Data(), 40, 0, 40);
   
   addHisto2D("hData_CRlowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);

   addHisto("hData_CRlowlowpt_BDTevt","",sample.Data(), 25,-1,1);
   addHisto("hData_CRlowlowpt_0Vtx_Mmumu","",sample.Data(), 25,0.,500. );
   addHisto("hData_CRlowlowpt_0Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowlowpt_1Vtx_Mmumu","",sample.Data(),  25,0.,500.);
   addHisto("hData_CRlowlowpt_1Vtx_Mass","",sample.Data(),   25,0.,100.);
   addHisto("hData_CRlowlowpt_1Vtx_BDTevt","",sample.Data(), 25,-1,1);
   addHisto("hData_CRlowlowpt_1Vtx_BDTvtx","",sample.Data(),25,-1,1);

   addHisto("hData_CRlowlowpt_1Vtx_SumtrackWeight","",sample.Data(), 40, 0, 40 );
   addHisto2D("hData_CRlowlowpt_1Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlowlowpt_TLVtx_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_CRlowlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRlowlowpt_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);

   addHisto("hData_CRlowlowpt_TLVtxAll_SumtrackWeight","",sample.Data(),40, 0, 40);

   addHisto("hData_CRlowlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);

   addHisto2D("hData_CRlowlowpt_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlowlowpt_2Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlowlowpt_2Vtx_Mass","",sample.Data(),   25,0.,100.);
   addHisto("hData_CRlowlowpt_2Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowlowpt_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1);
   addHisto("hData_CRlowlowpt_2Vtx_SumtrackWeight","",sample.Data(), 40, 0, 40);

   addHisto("hData_CRlowlowpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);

   addHisto("hData_CRlowlowpt_2VtxAll_BDTvtx","",sample.Data(),25,-1,1);

   addHisto("hData_CRlowlowpt_2VtxAll_SumtrackWeight","",sample.Data(),40, 0, 40);

   addHisto2D("hData_CRlowlowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRloose_BDTevt","",sample.Data(), 25,-1,1);
   addHisto("hData_CRloose_0Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRloose_0Vtx_BDTevt","",sample.Data(), 25,-1,1);
   addHisto("hData_CRloose_1Vtx_Mmumu","",sample.Data(),  25,0.,500.);
   addHisto("hData_CRloose_1Vtx_Mass","",sample.Data(),   25,0.,100. );
   addHisto("hData_CRloose_1Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRloose_1Vtx_BDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRloose_1Vtx_SumtrackWeight","",sample.Data(),40, 0, 40 );
   addHisto2D("hData_CRloose_1Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRloose_TLVtx_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_CRloose_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRloose_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);

   addHisto("hData_CRloose_TLVtxAll_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_CRloose_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRloose_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRloose_2Vtx_Mmumu","",sample.Data(),  25,0.,500.);
   addHisto("hData_CRloose_2Vtx_Mass","",sample.Data(),   25,0.,100. );
   addHisto("hData_CRloose_2Vtx_BDTevt","",sample.Data(), 25,-1,1);
   addHisto("hData_CRloose_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRloose_2Vtx_SumtrackWeight","",sample.Data(), 40, 0, 40 );
   addHisto2D("hData_CRloose_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRloose_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRloose_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRloose_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1 );           
   addHisto("hData_CRloose_2VtxAll_SumtrackWeight","",sample.Data(), 40, 0, 40 );

   addHisto("hData_CRlooselowpt_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlooselowpt_0Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlooselowpt_1Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlooselowpt_1Vtx_BDTvtx","",sample.Data(), 25,-1,1);
   addHisto("hData_CRlooselowpt_1Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlooselowpt_1Vtx_Mass","",sample.Data(),   25,0.,100. );
   addHisto("hData_CRlooselowpt_2Vtx_BDTevt","",sample.Data(), 25,-1,1);
   addHisto("hData_CRlooselowpt_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlooselowpt_1Vtx_SumtrackWeight","",sample.Data(),40, 0, 40 );
   addHisto2D("hData_CRlooselowpt_1Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlooselowpt_TLVtx_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_CRlooselowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);


   addHisto("hData_CRlooselowpt_TLVtxAll_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_CRlooselowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRlooselowpt_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselowpt_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlooselowpt_2Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlooselowpt_2Vtx_Mass","",sample.Data(),   25,0.,100. );

   addHisto("hData_CRlooselowpt_2Vtx_SumtrackWeight","",sample.Data(), 40, 0, 40);

   addHisto2D("hData_CRlooselowpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);

   addHisto("hData_CRlooselowpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRlooselowpt_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1 );


   addHisto("hData_CRlooselowpt_2VtxAll_SumtrackWeight","",sample.Data(), 40, 0, 40);

   addHisto("hData_CRlooselowlowpt_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlooselowlowpt_0Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlooselowlowpt_1Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlooselowlowpt_1Vtx_BDTvtx","",sample.Data(), 25,-1,1);
   addHisto("hData_CRlooselowlowpt_1Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlooselowlowpt_1Vtx_Mass","",sample.Data(),   25,0.,100.);
   addHisto("hData_CRlooselowlowpt_2Vtx_BDTevt","",sample.Data(),25,-1,1 );
   addHisto("hData_CRlooselowlowpt_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1 );

   addHisto("hData_CRlooselowlowpt_1Vtx_SumtrackWeight","",sample.Data(),40, 0, 40 );
   addHisto2D("hData_CRlooselowlowpt_1Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlooselowlowpt_TLVtx_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_CRlooselowlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);


   addHisto2D("hData_CRlooselowlowpt_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselowlowpt_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);


   addHisto("hData_CRlooselowlowpt_TLVtxAll_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_CRlooselowlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_CRlooselowlowpt_2Vtx_Mmumu","",sample.Data(), 25,0.,500. );
   addHisto("hData_CRlooselowlowpt_2Vtx_Mass","",sample.Data(),   25,0.,100. );

   addHisto("hData_CRlooselowlowpt_2Vtx_SumtrackWeight","",sample.Data(),40, 0, 40);

   addHisto2D("hData_CRlooselowlowpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselowlowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);


   addHisto("hData_CRlooselowlowpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRlooselowlowpt_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1 );


   addHisto("hData_CRlooselowlowpt_2VtxAll_SumtrackWeight","",sample.Data(), 40, 0, 40);


   addHisto("hData_StepEff","",sample.Data(),9, 0, 9);
   addHisto("hData_StepEff","NonNormalized",sample.Data(),9, 0, 9);
   
    // --------------- zzzzzzzzzzzzzzzzzzzzzzzzzzzz ------------------------------------//

//A
   addHisto("hData_CRtightlowlowpt_BDTevt","",sample.Data(),  25,-1,1  );
   addHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","",sample.Data(), 40, 0, 40);
   addHisto("hData_CRtightlowlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","",sample.Data(), 40, 0, 40);
   addHisto("hData_CRtightlowlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRtightlowlowpt_TLVtx_STW_Ntrks","",sample.Data(), 25,0,25,25,0,25);
   addHisto2D("hData_CRtightlowlowpt_TLVtxAll_STW_Ntrks","",sample.Data() ,25,0,25,25,0,25);

   //B
   addHisto("hData_CRlooselooselowlowpt_BDTevt","",sample.Data(),  25,-1,1  );
   addHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","",sample.Data(), 40, 0, 40);
   addHisto("hData_CRlooselooselowlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRlooselooselowlowpt_TLVtx_STW_Ntrks","",sample.Data() ,25,0,25,25,0,25);
   addHisto2D("hData_CRlooselooselowlowpt_TLVtxAll_STW_Ntrks","",sample.Data() ,25,0,25,25,0,25);
   addHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","",sample.Data(), 40, 0, 40);
   addHisto("hData_CRlooselooselowlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);


//C
   addHisto("hData_CRtighthighpt_BDTevt","",sample.Data(),25,-1,1  );
   addHisto("hData_CRtighthighpt_2Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRtighthighpt_2Vtx_Mass","",sample.Data(),   25,0.,100.);
   addHisto("hData_CRtighthighpt_2Vtx_SumtrackWeight","",sample.Data(), 40, 0, 40);
   addHisto2D("hData_CRtighthighpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRtighthighpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRtighthighpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRtighthighpt_2VtxAll_BDTvtx","",sample.Data(),  25,-1,1 );
   addHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","",sample.Data(), 40, 0, 40);
   
//D
   addHisto("hData_CRlooselooselowpt_BDTevt","",sample.Data(), 25,-1,1  );
   addHisto("hData_CRlooselooselowpt_2Vtx_Mmumu","",sample.Data(),   25,0.,500.);
   addHisto("hData_CRlooselooselowpt_2Vtx_Mass","",sample.Data(),  25,0.,100. );
   addHisto("hData_CRlooselooselowpt_2Vtx_SumtrackWeight","",sample.Data(), 40, 0, 40);
   addHisto2D("hData_CRlooselooselowpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselooselowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlooselooselowpt_2VtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_CRlooselooselowpt_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1  );
   addHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","",sample.Data(), 40, 0, 40);

//----- if using ABCD for 1 Vtx with combined regions ------//

   addHisto("hData_CRLoose_1Vtx_Mass","",sample.Data(),   25,0.,100. );
   addHisto("hData_CRLoose_1Vtx_SumtrackWeight","",sample.Data(),40, 0, 40);
   addHisto("hData_CRTight_1Vtx_Mass","",sample.Data(),   25,0.,100. );
   addHisto("hData_CRTight_1Vtx_SumtrackWeight","",sample.Data(),40, 0, 40);

}

//-------------------------------------------------------------
//instantiate and add
//first parameter is the variable name,
// second parameter is the selection step (like "afterleptsel")
//third parameter is the sample name (like "Z)
//others are TH1F binning
//creates one histograms per channel
//-------------------------------------------------------------
void TreeABCDReader::addHisto(TString var, TString selstep, TString sample, int nbins, float min, float max){
 
   TString name =  sample+"_"+var+"_"+selstep;
  TH1F * thehisto = new TH1F(name,name,nbins,min,max);
  thehisto->Sumw2();
  thehisto->SetOption("HIST");
   // // std::cout<<"adding histo with name : "<<name<<std::endl;
  histo_list_.push_back(thehisto);
  histo_map_[name.Data()] = numb_histo;
  numb_histo++;
}

void TreeABCDReader::addHisto2D(TString var, TString selstep, TString sample, int nxbins, float xmin, float xmax, int nybins, float ymin, float ymax){
 
  TString name =  sample+"_"+var+"_"+selstep;
  TH2F * thehisto = new TH2F(name,name,nxbins,xmin,xmax,nybins,ymin,ymax);
  // thehisto->Sumw2();
  thehisto->SetOption("COL");
   // // std::cout<<"adding histo with name : "<<name<<std::endl;
  histo_list_2D_.push_back(thehisto);
  histo_map_2D_[name.Data()] = numb_histo_2D_;
  numb_histo_2D_++;
}

//-------------------------------------------------------------
//fill histograms
//first parameter is the channel,
//second parameter is the variable name,
//third parameter is the selection step (like "afterleptsel")
//forths parameter is the sample name (like "Z)
//others are value and weight
//-------------------------------------------------------------
void TreeABCDReader::fillHisto( TString var, TString selstep,TString sample, float val, float weight){
  TString name = sample+"_"+var+"_"+selstep;

   // // std::cout<<"filling histo with name : "<<name<<std::endl;
  if(histo_map_[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histograms " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }else  histo_list_[histo_map_[name.Data()]]->Fill(val, weight);
  
}


void TreeABCDReader::fillHisto2D( TString var, TString selstep, TString sample, float xval,float yval, float weight){
  TString name = sample+"_"+var+"_"+selstep;

// // std::cout<<"filling histo with name : "<<name<<std::endl;
  if(histo_map_2D_[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histograms " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }else  histo_list_2D_[histo_map_2D_[name.Data()]]->Fill(xval,yval, weight);
  
}


void TreeABCDReader::deleteHisto(){
   cout << __LINE__ << endl;

   /*for(unsigned int i=0; i<histo_list_mmm.size(); i++){
     
     delete  histo_list_mmm[i];
     delete  histo_list_mme[i];
     delete  histo_list_eem[i];
     delete  histo_list_eee[i];
     
     
   }*/
   cout << __LINE__ << endl;
  //delete TheTree;
   cout << __LINE__ << endl;
  
  
}