#define DATAMCReader_cxx
#include "DATAMCReader.h"
#include <TH2.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>
// #include <chrono>
#include <ctime> 
#include "TTree.h"
#include "../../HistogramManager.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1F.h"
#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/MCWeights.h"

float DATAMCReader::MeanGenWeight(TString thesample, TString Prod)
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
   
void DATAMCReader::Loop(bool isMC, TString Prod, TString sample, bool Signal, int  Year, bool isPostAPV, float MeanGenW, int Channel ,bool DoubleMuon, std::vector<TString> thesystlist)
{

  TString thesample  = sample;
   TFile * theoutputfile = new TFile( ("../../"+Prod+"/DATAMC_"+thesample+".root").Data() , "recreate");
   // std::ofstream ofs (Prod+"/Efficacity_"+thesample+".txt", std::ofstream::out);
   // std::cout<<"syslist size "<<systlist.size()<<std::endl;
   for(unsigned int i=0; i< systlist.size(); i++){
     TString samplename = "";
     if( systlist[i]== "") samplename = thesample;
     else                  samplename = thesample+"__"+systlist[i];
     
     bool firstinit = false;
     if(i==0) firstinit = true;
     //cout << " iiii " << i << endl;
     //cout << samplename << endl;
     initializeHisto(samplename, firstinit);
   }
   // const long double sysTime = time(0);


   // Where L1 is always a muon and L2 is either a muon for the dimuon channel or a
   // an electron for the Emu channel
     
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
   // fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("h2D_SF_emu_lepABpt_FullError")); // x- axiss pt and y axis pt  
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

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   std::cout<<"Number of entries in the tree = "<<nentries<<std::endl;

   float XS = 1;
   float XS_up = 1.;
   float XS_down = 1.;

   float Nevent = nentries;

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

   int MSmuon = 0;
   if (thesample.Contains("DYJetsToLL_M-10to50"))                     { XS = 22635;   }//15910 :https://doi.org/10.48550/arXiv.2402.08486
   if (thesample.Contains("DYJetsToLL_M-50"))                         { XS = 6225.4;      }//5379.0;   }https://doi.org/10.48550/arXiv.2402.08486
   if (thesample.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays")) { XS = 21.6;   }
   if (thesample.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     { XS = 21.6;   }
   if (thesample.Contains("TTJets_DiLept"))                          { XS = 53.07;     }
   if (thesample.Contains("WWTo2L2Nu"))                              { XS = 11.09;     }
   if (thesample.Contains("WZTo2Q2L_mllmin4p0"))                     { XS = 6.535;     }
   if (thesample.Contains("ZZTo2Q2L_mllmin4p0"))                     { XS = 3.676;     }
   if (thesample.Contains("TTTo2L2Nu"))                              { XS = 88.5*2;      }
   if (thesample.Contains("ttWJetsToLNu_5f_EWK"))                    { XS = 0.290;     } // not found on XSDB, no file on tier2...approximation
         //Took 0.868 pb (CMS-TOP-21-011)
      // as a starting point and then divided by 3 (lepton universality)
   if (thesample.Contains("TTZToLL_5f"))                             { XS = 0.05188;   }//Not found on XSDB => used ana.py macro : 5.188e-02 +- 2.437e-04 pb
   if (thesample.Contains("TTToHadronic"))                           { XS = 378.9;     }
   if (thesample.Contains("TTWW"))                                   { XS = 0.006992;  }//found on XSDB
   if (thesample.Contains("TTToSemiLeptonic") )                      { XS = 365.34;    }

   //Signal
   if (thesample.Contains("smu200")) { XS = 0.01; MSmuon = 200;   }
   if (thesample.Contains("smu250")) { XS = 0.0045; MSmuon = 250; }
   if (thesample.Contains("smu300")) { XS = 0.002; MSmuon = 300; }
   if (thesample.Contains("smu350")) { XS = 0.001; MSmuon = 350; }
   if (thesample.Contains("smu400")) { XS = 0.0006; MSmuon = 400;}
   if (thesample.Contains("smu450")) { XS = 0.0004; MSmuon = 450;}
   if (thesample.Contains("smu500")) { XS = 0.00025;MSmuon = 500;}

// Signla XS are given in fb
int mixing = 0;
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
   // std::cout<<"here 2"<<std::endl;



   Long64_t nbytes = 0, nb = 0;
   float allevents = 0;
   // nentries = 100;
   float tempNentries = nentries;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      // if (jentry > 100000) break;
      allevents++;
      // if ( allevents%1000 == 0 ) std::cout << "events : " << allevents << std::endl;
      if ( allevents/tempNentries > 0.0999 && allevents/tempNentries < 0.1001 ) std::cout << "10pt events done" << std::endl;
      if (  allevents/tempNentries > 0.1999 && allevents/tempNentries < 0.2001 ) std::cout << "20pt of events done" << std::endl;
      if ( allevents/tempNentries > 0.2999 && allevents/tempNentries < 0.3001) std::cout << "30pt f events done" << std::endl;
      if (  allevents/tempNentries > 0.3999 && allevents/tempNentries < 0.4001) std::cout << "40pt  events done" << std::endl;
      if (  allevents/tempNentries > 0.4999 && allevents/tempNentries < 0.5001 ) std::cout << "50pt  events done" << std::endl;
      if (  allevents/tempNentries > 0.5999 && allevents/tempNentries < 0.6001) std::cout << "60pt  events done" << std::endl;
      if (  allevents/tempNentries > 0.6999 && allevents/tempNentries < 0.7001) std::cout << "70pt  events done" << std::endl;
      if (  allevents/tempNentries > 0.7999 && allevents/tempNentries < 0.8001) std::cout << "80pt  events done" << std::endl;
      if (  allevents/tempNentries > 0.8999 && allevents/tempNentries < 0.9001) std::cout << "90pt  events done" << std::endl;
      if (  allevents/tempNentries > 0.9499 && allevents/tempNentries < 0.9501) std::cout << "95pt  events done" << std::endl;
      // std::cout<<"allevents : "<<allevents<<std::endl;
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


      if ( minitree_Hemi_pt->size() != 2 ) continue;
      if ( abs(minitree_Hemi_eta->at(0)) > 2.4 || abs(minitree_Hemi_eta->at(1)) > 2.4 ) continue;
      if ( minitree_Hemi_pt->at(0) < 20. || minitree_Hemi_pt->at(1) < 20. ) continue;

      // if (abs(minitree_lepton_leadingeta2->at(0))>1.5) continue;
      // !!
      if (abs(minitree_lepton_leadingeta2->at(0))>1.5  && MuMu==0) 
         {
            SFele_eta = hEle_SF->GetBinContent(hEle_SF->GetXaxis()->FindBin(minitree_lepton_leadingeta2->at(0)));
            // std::cout<<"SFele_eta : "<<SFele_eta<<std::endl;
            // SFele_pt_eta = hEle_2DSF->GetBinContent(hEle_2DSF->GetXaxis()->FindBin(minitree_lepton_leadingpt2->at(0)),hEle_2DSF->GetYaxis()->FindBin(minitree_lepton_leadingeta2->at(0)));
         }
      // !! 
      if (minitree_Filter->at(0))
         {

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
                       
         }// End of minitree_filter
      // --------- End of Lepton SF ------//
   
      // --------- Event SF ------//
   float GenWeight = minitree_only_gen_wt->at(0)/MeanGenW;
   float PUweight = 1.;
   float PUweightUp = 1.;
   float PUweightDown = 1.;
   
   float Prefweight = 1.;
   float PrefweightUp = 1.;
   float PrefweightDown = 1.;

   float TriggerSyst = triggerSF;
   Ele_SF =Ele_Reco_SF1*Ele_ID_SF1;
   Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * Mu_ISO_SF1;
   Mu_SF2 = Mu_Reco_SF2* Mu_ID_SF2 * Mu_ISO_SF2;


   float NormFactorSYST = NormFactorLumi/norm;
   if (Signal) GenWeight = 1;
   if (!isMC)//<=> Data
      {
         NormFactorSYST = 1;                                                                                                                                                               
      }
   else{
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
            // std::cout<<"NormFactorSYST="<<NormFactorSYST<<" GenWeight="<<GenWeight<<" miniPrefweight="<<Prefweight<<" miniPUweight="<<PUweight<<" top_pt_wt="<<top_pt_wt<<" Mu_SF="<<Mu_SF<<" Ele_SF="<<Ele_SF<<" EMu_trig_SF1="<<EMu_trig_SF1<<std::endl;
            // std::cout<<" XS = "<<XS<<" Lumi = "<<Lumi<<" nentries = "<<nentries<<std::endl;
            // std::cout<<"MeanGenW : "<<MeanGenW<<"and genweight : "<<minitree_only_gen_wt->at(0)<<std::endl;
            //  std::cout<<"SFele_eta = "<<SFele_eta<<" and the eta of leadingletpno 2 : "<< minitree_lepton_leadingeta2->at(0)<<std::endl;
            // std::cout<<"NormFactor e= "<<NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Ele_SF*EMu_trig_SF1<<" NormFactor_e cor = "<<SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Ele_SF*EMu_trig_SF1<<std::endl;
            // std::cout<<"NormFactor = "<<NormFactor<<" NormFactor_mu = "<<NormFactor_mu<<" NormFactor_ele = "<<NormFactor_ele<<" Pref_PU_gen_wt = "<<Pref_PU_gen_wt<<std::endl;
      // !! ---------------------------------------------------------------------------- !! //
      // --------- End of Event SF ------//

      fillHisto("Tree_filter", "", thesample, minitree_Filter->at(0), NormFactor);// Filter                                                                                                                          
      fillHisto("Vertices", "NoSel", thesample,minitree_nPV->at(0), Pref_PU_gen_wt);// no sel                                                                                                                              

      if (minitree_Filter->at(0) || minitree_FilterSameSign->at(0)){                                                                                                                                                                     
            
            fillHisto("Vertices_filtercut", "", thesample,minitree_nPV->at(0), NormFactor);// no sel               
            fillHisto("njet","NoSel",  thesample,minitree_njet->at(0), NormFactor);
            fillHisto("njetNOmu","NoSel",  thesample,minitree_njetNOmu->at(0), NormFactor);
         }

      fillHisto("Tree_Mumu","nosel", thesample,  minitree_Mmumu->at(0), NormFactor);/// offline only                                                                                                              

      if (MuMu == 0)
         {
            for(unsigned int iMuon = 0; iMuon <minitree_lepton_leadingpt->size(); iMuon++)//test    after trigger and lepton selection cut                                       
               {

                  fillHisto("leading_muon_pt", "reco",  thesample, minitree_lepton_leadingpt->at(iMuon), NormFactor_mu);
                  fillHisto("leading_muon_eta", "reco",  thesample, minitree_lepton_leadingeta->at(iMuon), NormFactor_mu);
                  fillHisto("leading_muon_phi", "reco",  thesample, minitree_lepton_leadingphi->at(iMuon), NormFactor_mu);
                  fillHisto("leading_muon_dxy", "reco", thesample,minitree_lepton_leadingdxy->at(iMuon), NormFactor_mu );
                  fillHisto("leading_muon_dz", "reco", thesample,minitree_lepton_leadingdz->at(iMuon), NormFactor_mu);

                  fillHisto("leading_lepton_pt", "reco", thesample, minitree_lepton_leadingpt2->at(iMuon), NormFactor_ele);
                  fillHisto("leading_lepton_eta", "reco", thesample, minitree_lepton_leadingeta2->at(iMuon), NormFactor_ele);
                  fillHisto("leading_lepton_phi", "reco", thesample, minitree_lepton_leadingphi2->at(iMuon), NormFactor_ele);
                  fillHisto("leading_lepton_dxy", "reco", thesample,minitree_lepton_leadingdxy2->at(iMuon), NormFactor_ele);
                  fillHisto("leading_lepton_dz", "reco", thesample,minitree_lepton_leadingdz2->at(iMuon), NormFactor_ele);
                  fillHisto2D("leading_lepton_pt_eta", "reco", thesample, minitree_lepton_leadingpt2->at(iMuon), minitree_lepton_leadingeta2->at(iMuon), NormFactor_ele);
               }
               
         }
      else if (MuMu > 0)
         {
            for(unsigned int iMuon = 0; iMuon <minitree_lepton_leadingpt->size(); iMuon ++)//test    after trigger and lepton selection cut                                       
               {
                  //if(minitree_muon_dxy->at(iMuon)< 0.1  && minitree_muon_dz->at(iMuon) < 0.2 ){ 
                  fillHisto("leading_muon_pt", "reco",  thesample, minitree_lepton_leadingpt->at(iMuon), NormFactor_mu);
                  fillHisto("leading_muon_eta", "reco",  thesample, minitree_lepton_leadingeta->at(iMuon), NormFactor_mu);
                  fillHisto("leading_muon_phi", "reco",  thesample, minitree_lepton_leadingphi->at(iMuon), NormFactor_mu);
                  fillHisto("leading_muon_dxy", "reco", thesample,minitree_lepton_leadingdxy->at(iMuon), NormFactor_mu );
                  fillHisto("leading_muon_dz", "reco", thesample,minitree_lepton_leadingdz->at(iMuon), NormFactor_mu);
                  
                  //}//muon impact
               }
               
               for(unsigned int iel = 0; iel <minitree_lepton_leadingpt2->size(); iel ++)//test    after trigger and lepton selection cut                                      
               {
                  //if (  ( fabs(minitree_reco_electron_leadingeta2->at(iel)) <= 1.479 && abs(minitree_electron_dxy->at(iel)) < 0.05 && abs(minitree_electron_dz->at(iel) ) < 0.10 ) || ( fabs(minitree_reco_electron_leadingeta2->at(iel)) > 1.556 && abs(minitree_electron_dxy->at(iel)) < 0.10 && abs(minitree_electron_dz->at(iel)) < 0.20 )){

                  fillHisto("leading_lepton_pt", "reco", thesample, minitree_lepton_leadingpt2->at(iel), NormFactor_mu2);
                  fillHisto("leading_lepton_eta", "reco", thesample, minitree_lepton_leadingeta2->at(iel),NormFactor_mu2);
                  fillHisto("leading_lepton_phi", "reco", thesample, minitree_lepton_leadingphi2->at(iel), NormFactor_mu2);
                  fillHisto("leading_lepton_dxy", "reco", thesample,minitree_lepton_leadingdxy2->at(iel), NormFactor_mu2);
                  fillHisto("leading_lepton_dz", "reco", thesample,minitree_lepton_leadingdz2->at(iel), NormFactor_mu2);
                  fillHisto2D("leading_lepton_pt_eta", "reco", thesample, minitree_lepton_leadingpt2->at(iel), minitree_lepton_leadingeta2->at(iel), NormFactor_mu2);
                  //}//ele impact para
               }
         }

      for(unsigned int ill = 0; ill <minitree_ll_pt->size(); ill ++)//test                                                                                                        
          {
            fillHisto("Dilepton_pt", "NoSel", thesample, minitree_ll_pt->at(ill), NormFactor);
            fillHisto("Dilepton_eta", "NoSel", thesample, minitree_ll_eta->at(ill), NormFactor);
            fillHisto("Dilepton_phi", "NoSel", thesample, minitree_ll_phi->at(ill), NormFactor);
            fillHisto("Dilepton_mass", "NoSel", thesample, minitree_ll_mass->at(ill), NormFactor);
         }

      //----------------------------------------------------------------//
      int FilterSample = -1;
      if (MuMu == 2 || MuMu == 1)
         {
            if (DoubleMuon && minitree_trigger_doublelepton->at(0)  ) //&& !minitree_trigger_singlelepton
               {
                  FilterSample = 2;
               }
            else if  (!DoubleMuon && minitree_trigger_singlelepton->at(0) &&  !minitree_trigger_doublelepton->at(0) ) 
               {
                  FilterSample = 1;
               }
         }
      if (MuMu == 0)
         {
               FilterSample = 0;
         }
      // std::cout<<"FilterSample : "<<FilterSample<<std::endl;
      if (DoubleMuon && FilterSample != 2) continue;
      else if (!DoubleMuon && FilterSample != 1 && MuMu ==1) continue;
      else if (!DoubleMuon && FilterSample != 0 && MuMu ==0) continue;

      //----------------------------------------------------------------//

      if (  minitree_Filter->at(0) && minitree_njetNOmu->at(0)>= 1) // !!!!! pour EMu: tester minitree_lepton_leadingeta2->at(iel)<1.5 (barrel region)
         {
            //   std::cout<<" here 0  "<<std::endl;
            fillHisto("LeadingLeptons_dR","",               thesample,  minitree_lepton_lepton_dR->at(0),NormFactor);//filter *1
            fillHisto("LeadingLeptons_dPhi","",             thesample,  minitree_lepton_lepton_dPhi->at(0),NormFactor);//filter *1
            fillHisto("LeadingJets_dR","",                  thesample,  minitree_jet_jet_dR->at(0),NormFactor);//jetnoMu>1 *1
            fillHisto("LeadingJets_dPhi","",                thesample,  minitree_jet_jet_dPhi->at(0),NormFactor);//jetnoMu>1 * 1
            // std::cout<<" here 0A  "<<std::endl;
            fillHisto("LeadingLeptonJet_dRmax","",           thesample,  minitree_muon_jet_dRmax->at(0),NormFactor);//jetnoMu>1 *2
            fillHisto("LeadingLeptonJet_dRmin","",           thesample,  minitree_muon_jet_dRmin->at(0),NormFactor);//jetnoMu>1 *2
            fillHisto("HemiAxis_Mu_dR","",                   thesample,  minitree_HemiMu_dR->at(0),NormFactor);//jetnoMu>1 *2
            fillHisto("HemiAxis_OpMu_dR","",                 thesample,  minitree_HemiMuOp_dR->at(0),NormFactor);//jetnoMu>1 *2

            fillHisto("Nmu","",                 thesample, minitree_nmu->at(0),NormFactor);

            for (unsigned int i = 0 ; i < minitree_muon_pt->size(); i++)
               {
                  fillHisto("Muon_pt","",             thesample, minitree_muon_pt->at(0),NormFactor);
                  fillHisto("Muon_PFIsoLoose","",     thesample, minitree_muon_PFIsoLoose->at(0),NormFactor);
                  fillHisto("Muon_MiniIsoTight","",   thesample, minitree_muon_MiniIsoTight->at(0),NormFactor);
               }

            fillHisto("HT","",            thesample,  minitree_HT->at(0),NormFactor);
            fillHisto("LT","",            thesample,  minitree_LT->at(0),NormFactor);
            fillHisto("nTrks","",         thesample,  minitree_nTracks->at(0),NormFactor);
            fillHisto("nLostTracks","",   thesample,  minitree_nLostTracks->at(0),NormFactor);

            // --- CMSSW Version of V0 Candidates
            for (unsigned int i = 0 ; i < minitree_K0_mass->size(); i++ )
               {
                  fillHisto("K0_mass","", thesample, minitree_K0_mass->at(i),NormFactor);
                  fillHisto("K0_pt","", thesample, minitree_K0_pt->at(i),NormFactor);
               }
            for (unsigned int i = 0 ; i < minitree_L0_mass->size(); i++ )
               {
                  fillHisto("L0_mass","", thesample,minitree_L0_mass->at(i),NormFactor);
                  fillHisto("L0_pt","", thesample, minitree_L0_pt->at(i),NormFactor);
               }

            // -- Our own collection of V0 Candidates
            for (unsigned int i = 0 ; i <minitree_V0_reco_source->size(); i++ )
               {
                     if (minitree_V0_reco_source->at(i) == 1) 
                        {
                           fillHisto("Reco_K0_mass","",thesample,minitree_V0_reco_mass->at(i),NormFactor);
                           fillHisto("Reco_K0_pt","", thesample, minitree_V0_reco_pt->at(i),NormFactor);
                        }
                     else if (minitree_V0_reco_source->at(i) == 2) 
                        {
                           fillHisto("Reco_L0_mass","", thesample,minitree_V0_reco_mass->at(i),NormFactor);
                           fillHisto("Reco_L0_pt","",thesample, minitree_V0_reco_pt->at(i),NormFactor);
                        }

               }

            fillHisto("njet","NoSel",  thesample,minitree_njet->at(0),NormFactor);
            fillHisto("njetNOmu","NoSel", thesample,minitree_njetNOmu->at(0),NormFactor);

            for (unsigned int k = 0 ; k < minitree_jet_pt->size(); k++)
               {
                  fillHisto("hData_jet_pt","",            thesample, minitree_jet_pt->at(k),NormFactor);
                  fillHisto("hData_jet_eta","",           thesample, minitree_jet_eta->at(k),NormFactor);
                  fillHisto("hData_jet_btag_Deepjet","",  thesample, minitree_jet_btag_DeepJet->at(k),NormFactor);
                  fillHisto("hData_jet_HadronFlavour","", thesample, minitree_jet_HadronFlavour->at(k),NormFactor);
               }
            fillHisto("leading_jet_pt","",            thesample,  minitree_jet_leadingpt->at(0),NormFactor);
            fillHisto("leading_jet_eta","",           thesample,  minitree_jet_leadingeta->at(0),NormFactor);
            fillHisto("subleading_jet_pt","",         thesample,  minitree_jet_leadingpt2->at(0),NormFactor);
            fillHisto("subleading_jet_eta","",        thesample,  minitree_jet_leadingeta2->at(0),NormFactor);

            for(unsigned int iTrk = 0; iTrk <minitree_track_MVAval->size(); iTrk ++)
               {
                  //--- Track pre-selected ---------//                                   
                  //get total number of tracks bfore any selection with minitree_TRACK_SIZE                                                                                                      
                  // fillHisto("track_nTracks","TRK", thesample,   minitree_nTracks->at(0),NormFactor); 
                  // fillHisto("track_nLostTracks","TRK", thesample, minitree_nLostTracks->at(0),NormFactor);
                  fillHisto("track_ipc","TRK", thesample,  minitree_track_ipc->at(iTrk),NormFactor);
                  
                  //fillHisto("track_charge","TRK", thesample, minitree_track_charge->at(iTrk),NormFactor);
                  //fillHisto("track_isHighPurity","TRK", thesample, minitree_track_isHighPurity->at(iTrk),NormFactor);

                  fillHisto("track_dxyError","TRK", thesample, minitree_track_dxyError->at(iTrk),NormFactor);

                  fillHisto("track_nHitPixel","TRK", thesample, minitree_track_nHitPixel->at(iTrk),NormFactor);
                  fillHisto("track_nHitTIB","TRK", thesample, minitree_track_nHitTIB->at(iTrk),NormFactor);
                  fillHisto("track_nHitTOB","TRK", thesample, minitree_track_nHitTOB->at(iTrk),NormFactor);
                  fillHisto("track_nHitTEC","TRK", thesample, minitree_track_nHitTEC->at(iTrk),NormFactor);
                  fillHisto("track_nHitPXB","TRK", thesample, minitree_track_nHitPXB->at(iTrk),NormFactor);
                  fillHisto("track_nHitPXF","TRK", thesample, minitree_track_nHitPXF->at(iTrk),NormFactor);
                  fillHisto("track_isHitPixel","TRK", thesample, minitree_track_isHitPixel->at(iTrk),NormFactor);
                  
                  fillHisto("track_nLayers","TRK", thesample, minitree_track_nLayers->at(iTrk),NormFactor);
                  fillHisto("track_nLayersPixel","TRK", thesample, minitree_track_nLayersPixel->at(iTrk),NormFactor);
                  fillHisto("track_region","TRK", thesample, minitree_track_region->at(iTrk),NormFactor);
                  // fillHisto("track_btag","TRK", thesample, minitree_track_btag->at(iTrk),NormFactor);
                  //fillHisto("track_energy","TRK", thesample, minitree_track_energy->at(iTrk),NormFactor);
                  // fillHisto("track_Hemi","TRK", thesample, minitree_track_Hemi->at(iTrk),NormFactor);
                     
                  // fillHisto("track_TRACK_SIZE","TRK", thesample,minitree_TRACK_SIZE->at(0),NormFactor);
                  

                  fillHisto("track_dxy","TRK", thesample,minitree_track_dxy->at(iTrk),NormFactor);
                  fillHisto("track_dz","TRK", thesample,minitree_track_dz->at(iTrk),NormFactor);
                  fillHisto("track_pt","TRK", thesample,minitree_track_pt->at(iTrk),NormFactor);
                  fillHisto("track_eta","TRK", thesample,minitree_track_eta->at(iTrk),NormFactor);
                  fillHisto("track_NChi2","TRK", thesample,minitree_track_NChi2->at(iTrk),NormFactor);
                  fillHisto("track_nhits","TRK", thesample,minitree_track_nHit->at(iTrk),NormFactor);
                  fillHisto("track_ntrk10","TRK", thesample,minitree_track_ntrk10->at(iTrk),NormFactor);
                  fillHisto("track_ntrk20","TRK", thesample,minitree_track_ntrk20->at(iTrk),NormFactor);
                  fillHisto("track_ntrk30","TRK", thesample,minitree_track_ntrk30->at(iTrk),NormFactor);
                  fillHisto("track_ntrk40","TRK", thesample,minitree_track_ntrk40->at(iTrk),NormFactor);
                  fillHisto("track_drSig","TRK", thesample,minitree_track_drSig->at(iTrk),NormFactor);
                  fillHisto("track_dzSig","TRK", thesample,minitree_track_dzSig->at(iTrk),NormFactor);
                  fillHisto("track_iJet","TRK", thesample,minitree_track_iJet->at(iTrk),NormFactor);// tracks matched with jets dR< 0.4
                  fillHisto("track_track_Hemi_dR","TRK", thesample,minitree_track_Hemi_dR->at(iTrk),NormFactor);// dR b/w tracks and jets
                  fillHisto("track_track_Hemi_dRmax","TRK", thesample,minitree_track_Hemi_dRmax->at(iTrk),NormFactor);
                  fillHisto("track_lost","TRK", thesample,minitree_track_lost->at(iTrk),NormFactor);
                  fillHisto("track_MVAVal","TRK", thesample,minitree_track_MVAval->at(iTrk),NormFactor);

                  fillHisto("track_Track_firstHit","TRK", thesample,minitree_track_firstHit->at(iTrk),NormFactor);
                  fillHisto("track_Track_firstHit_x","TRK", thesample,minitree_track_firstHit_x->at(iTrk),NormFactor);
                  fillHisto("track_Track_firstHit_y","TRK", thesample,minitree_track_firstHit_y->at(iTrk),NormFactor);
                  fillHisto("track_Track_firstHit_z","TRK", thesample,minitree_track_firstHit_z->at(iTrk),NormFactor);
                  fillHisto2D("track_Track_firstHit_X_Y",  "TRK", thesample,minitree_track_firstHit_x->at(iTrk),minitree_track_firstHit_y->at(iTrk),NormFactor);
                  
                  float firstHit_x =minitree_track_firstHit_x->at(iTrk);
                  float firstHit_y =minitree_track_firstHit_y->at(iTrk);
                  float firstHit_r=sqrt(firstHit_x*firstHit_x+firstHit_y*firstHit_y);
                  
                  fillHisto("track_Track_firstHit_r","TRK", thesample,firstHit_r,NormFactor);
                                    
                  float mva_ntrk10                  =  minitree_track_ntrk10->at(iTrk);
                  float mva_ntrk20                  =  minitree_track_ntrk20->at(iTrk);
                  float mva_ntrk30                  =  minitree_track_ntrk30->at(iTrk);
                  float mva_ntrk40                  =  minitree_track_ntrk40->at(iTrk);
                  
                  float mva_ntrk10rel;
                  float mva_ntrk20rel;
                  float mva_ntrk30rel;
                  float mva_ntrk10reltot;
                  float mva_ntrk20reltot;
                  float mva_ntrk30reltot;
                  float mva_ntrk40reltot;
                  
                  if ( mva_ntrk40 == 0 ) {
                  mva_ntrk10rel = -0.02;
                  mva_ntrk20rel = -0.02;
                  mva_ntrk30rel = -0.02; 
                  }
                  else {
                  mva_ntrk10rel = mva_ntrk10 / mva_ntrk40;
                  mva_ntrk20rel = mva_ntrk20 / mva_ntrk40;
                  mva_ntrk30rel = mva_ntrk30 / mva_ntrk40;
                  }
                  
               }//end loop on tracks

         float Hemisphere1_pt = 0;
         float Hemisphere2_pt = 0;

         Hemisphere1_pt = minitree_Hemi_pt->at(0);
         Hemisphere2_pt = minitree_Hemi_pt->at(1);
         float hemi_ptmin = Hemisphere2_pt;
         float hemi_ptmax = Hemisphere1_pt;
         if ( Hemisphere2_pt > Hemisphere1_pt ) {hemi_ptmin = Hemisphere1_pt;hemi_ptmax = Hemisphere2_pt;}
         fillHisto("Hemisphere_leadingpt","", thesample, hemi_ptmax,NormFactor);
         fillHisto("Hemisphere_subleadingpt","", thesample, hemi_ptmin,NormFactor);
         for (unsigned int i = 0 ; i < minitree_Hemi_pt->size(); i++)
            {
               fillHisto("Hemi_pt","",       thesample,  minitree_Hemi_pt->at(i),NormFactor);
               fillHisto("Hemi_eta","",      thesample,  minitree_Hemi_eta->at(i),NormFactor);
               fillHisto("Hemi_nJet","",     thesample,  minitree_Hemi_njet->at(i),NormFactor);
               fillHisto("Hemi_nJetNoMu","",thesample,   minitree_Hemi_njet_nomu->at(i),NormFactor);

               fillHisto("HemiMu_pt","",    thesample,  minitree_HemiMu_pt->at(i),NormFactor);
               fillHisto("HemiMu_Mass","",  thesample,   minitree_HemiMu_mass->at(i),NormFactor);
               fillHisto("Hemi_nTrks","",   thesample,  minitree_Hemi_nTrks->at(i),NormFactor);
               fillHisto("Hemi_Mass","",    thesample,  minitree_Hemi_mass->at(i),NormFactor);
            }

         for (unsigned int i = 0 ; i < minitree_SecInt_pt->size(); i++)
            {
               if (minitree_SecInt_selec->at(i))
                  {
                     fillHisto("SecInt_mass","Selec", thesample, minitree_SecInt_mass->at(i), NormFactor);
                     fillHisto("SecInt_drSig","Selec", thesample, minitree_SecInt_drSig->at(i), NormFactor);
                     fillHisto("SecInt_pt","Selec", thesample, minitree_SecInt_pt->at(i), NormFactor);
                     fillHisto("SecInt_dzSig","Selec", thesample, minitree_SecInt_dzSig->at(i), NormFactor);
                     fillHisto("SecInt_selec","Selec", thesample, minitree_SecInt_selec->at(i), NormFactor);
                     fillHisto("SecInt_layer","Selec", thesample, minitree_SecInt_layer->at(i), NormFactor);
                     fillHisto("SecInt_r","Selec", thesample, minitree_SecInt_r->at(i), NormFactor);
                     fillHisto("SecInt_z","Selec", thesample, minitree_SecInt_z->at(i), NormFactor);
                     if ( minitree_SecInt_layer->at(i) )
                        {
                           fillHisto("SecInt_mass","TrackerMatched", thesample, minitree_SecInt_mass->at(i), NormFactor);
                           fillHisto("SecInt_drSig","TrackerMatched", thesample, minitree_SecInt_drSig->at(i), NormFactor);
                           fillHisto("SecInt_pt","TrackerMatched", thesample, minitree_SecInt_pt->at(i), NormFactor);
                           fillHisto("SecInt_dzSig","TrackerMatched", thesample, minitree_SecInt_dzSig->at(i), NormFactor);
                           fillHisto("SecInt_selec","TrackerMatched", thesample, minitree_SecInt_selec->at(i), NormFactor);
                           fillHisto("SecInt_layer","TrackerMatched", thesample, minitree_SecInt_layer->at(i), NormFactor);
                        }
                  }
            }
         // std::cout<<" here 5  "<<std::endl;
         for (unsigned int i = 0 ; i < minitree_Hemi_Vtx_Mass->size(); i++)
            {
               fillHisto("Vtx_NChi2","",            thesample,  minitree_Hemi_Vtx_NChi2->at(i),NormFactor);
               fillHisto("Vtx_nTrks","",           thesample,  minitree_Hemi_Vtx_nTrks->at(i),NormFactor);
               fillHisto("Vtx_Mass","",         thesample,  minitree_Hemi_Vtx_Mass->at(i),NormFactor);
               fillHisto("Vtx_Dist","",        thesample,   minitree_Hemi_Vtx_dist->at(i),NormFactor);

               fillHisto("Vtx_Step","",            thesample,  minitree_Hemi_Vtx_step->at(i),NormFactor);
               fillHisto("Vtx_dR","",           thesample,  minitree_Hemi_Vtx_dR->at(i),NormFactor);
            }
         // std::cout<<" here 6  "<<std::endl;
         for (unsigned int i = 0 ; i < minitree_Hemi_SecVtx_Mass->size(); i++)
            {
               fillHisto("SecVtx_NChi2","",            thesample,  minitree_Hemi_SecVtx_NChi2->at(i),NormFactor);
               fillHisto("SecVtx_nTrks","",           thesample,  minitree_Hemi_SecVtx_nTrks->at(i),NormFactor);
               fillHisto("SecVtx_Mass","",         thesample,  minitree_Hemi_SecVtx_Mass->at(i),NormFactor);
               fillHisto("SecVtx_Dist","",        thesample,   minitree_Hemi_SecVtx_dist->at(i),NormFactor);

               fillHisto("SecVtx_Step","",         thesample,  minitree_Hemi_SecVtx_step->at(i),NormFactor);
               fillHisto("SecVtx_dR","",        thesample,minitree_Hemi_SecVtx_dR->at(i),NormFactor);
            }

         // std::cout<<" here 7  "<<std::endl;
         for (unsigned int i = 0 ; i < minitree_Hemi_Vtx_BDT_nTrks->size() ; i++)
            {
               fillHisto("FinalVtx_nTrks","",           thesample,  minitree_Hemi_Vtx_BDT_nTrks->at(i),NormFactor);
               fillHisto("FinalVtx_Step","",          thesample, minitree_Hemi_Vtx_BDT_step->at(i),NormFactor);
               fillHisto("FinalVtx_Mass","",         thesample,  minitree_Hemi_Vtx_BDT_HMass->at(i),NormFactor);
               fillHisto("FinalVtx_HMass","",        thesample,  minitree_Hemi_Vtx_BDT_HMass->at(i),NormFactor);
            }
         // std::cout<<" here 8  "<<std::endl;

      }// end filter

   }// end loop on events
   // ofs.close();
   theoutputfile->Write();
   //deleteHisto();
   theoutputfile->Close();
   delete theoutputfile;
} // end loop method

void DATAMCReader::initializeHisto(TString sample, bool isfirstset){


  std::cout << "#####################################" << std::endl;
  std::cout << "#####################################" << std::endl;
  std::cout << " initialize histograms of sample :  " <<sample<< std::endl;
  std::cout << "#####################################" << std::endl;
  std::cout << "#####################################" << std::endl;
  

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

   // ---------------------------- Filter (Online+Offline Selection)--------------
// Dimuon : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
// OR
// Single Muon : HLT_IsoMu24_v
// Offline selection of muons : Two Loose ID and PFIsoLoose muons required
// |dxy | < 0.1 cm & |dz | < 0.2 cm (prompt)
// pt1 > 25GeV & pt2 > 10GeV
// Mμμ > 10 GeV (remove low-resonances)

//----------------//

   addHisto("Tree_filter", "", sample.Data(), 2,0,2);
   addHisto("Vertices", "NoSel", sample.Data(), 100,0,100);
   addHisto("Vertices_filtercut","", sample.Data(), 100,0,100); 
   addHisto("Tree_Mumu","nosel", sample.Data(),  150, 0 ,600);
   
   addHisto("Dilepton_pt","NoSel", sample.Data(),100, 0 ,500 );
   addHisto("Dilepton_eta","NoSel",sample.Data(),25,-2.4,2.4);
   addHisto("Dilepton_phi","NoSel",sample.Data(),25,-3.14,3.14);
   addHisto("Dilepton_mass","NoSel",sample.Data(),150, 0 ,600);

   addHisto("leading_muon_pt", "reco", sample.Data(),100, 0 ,400 );
   addHisto("leading_lepton_pt", "reco", sample.Data(),100, 0 ,400);
   addHisto("leading_muon_eta", "reco", sample.Data(),25,-2.4,2.4);
   addHisto("leading_lepton_eta", "reco", sample.Data(),25,-2.4,2.4);

   addHisto("leading_muon_phi", "reco", sample.Data(),25,-3.14,3.14);
   addHisto("leading_lepton_phi", "reco", sample.Data(),25,-3.14,3.14);

   addHisto2D("leading_lepton_pt_eta", "reco", sample.Data(),40,0,200,25,-2.4,2.4);

   // !! --------------------------- being tested 
   addHisto("leading_muon_dxy", "reco", sample.Data(),40,-0.2,0.2 );
   addHisto("leading_lepton_dxy", "reco", sample.Data(),200,-1,1);
   addHisto("leading_muon_dz", "reco", sample.Data(),50,-0.5,0.5);
   addHisto("leading_lepton_dz", "reco", sample.Data(),200,-1,1);


   addHisto("njet","NoSel",  sample.Data(),20,0,20);
   addHisto("njetNOmu","NoSel",  sample.Data(),20,0,20);
   addHisto("Hemisphere_leadingpt","", sample.Data(),100 ,0,500);
   addHisto("Hemisphere_subleadingpt","", sample.Data(), 100 ,0,500);
   // !!  ------------------------------end 

   // ------------------------ Jets   -------------------------------------------
   addHisto("hData_jet_pt","",            sample.Data(),  100, 0 ,500);
   addHisto("hData_jet_eta","",           sample.Data(),  55,-5,5);
   addHisto("hData_jet_btag_Deepjet","",  sample.Data(),  50,0,1);
   addHisto("hData_jet_HadronFlavour","", sample.Data(),  7,-0.5,6.5);

   // !! ------------------- bgin 

   addHisto("leading_jet_pt","",            sample.Data(),  100, 0 ,500);
   addHisto("leading_jet_eta","",           sample.Data(),  55,-2.5,2.5);
   addHisto("subleading_jet_pt","",         sample.Data(),  100, 0 ,500);
   addHisto("subleading_jet_eta","",        sample.Data(),  55,-2.5,2.5);

   addHisto("Nmu","",            sample.Data(),  10,0,10);
   addHisto("Muon_pt","",           sample.Data(),  100,0,500);
   addHisto("Muon_PFIsoLoose","",         sample.Data(),  2,0,2);
   addHisto("Muon_MiniIsoTight","",        sample.Data(),  2,0,2);

   addHisto("LeadingLeptons_dR","",            sample.Data(),  50,0,5);
   addHisto("LeadingLeptons_dPhi","",            sample.Data(),  35,0,3.5);
   addHisto("LeadingJets_dR","",            sample.Data(),  50,0,5);
   addHisto("LeadingJets_dPhi","",            sample.Data(),  35,0,3.5);

   addHisto("LeadingLeptonJet_dRmax","",            sample.Data(),  50,0,5);
   addHisto("LeadingLeptonJet_dRmin","",            sample.Data(),  50,0,5);
   addHisto("HemiAxis_Mu_dR","",                    sample.Data(),  50,0,5);
   addHisto("HemiAxis_OpMu_dR","",                  sample.Data(),  50,0,5);

   addHisto("HT","",            sample.Data(),  200, 0 ,800);
   addHisto("LT","",           sample.Data(),  100, 0 ,500);
   addHisto("nTrks","",         sample.Data(),  100,0,100);
   addHisto("nLostTracks","",        sample.Data(),  25,0,25);

   addHisto("Hemi_pt","",            sample.Data(),  100, 0 ,600);
   addHisto("Hemi_eta","",           sample.Data(),  55,-2.5,2.5);
   addHisto("Hemi_nJet","",         sample.Data(),  10,0,10);
   addHisto("Hemi_nJetNoMu","",        sample.Data(),  10,0,10);

   addHisto("HemiMu_pt","",            sample.Data(),  100, 0 ,600);
   addHisto("HemiMu_Mass","",           sample.Data(), 100, 0 ,600);
   addHisto("Hemi_nTrks","",         sample.Data(),  25,0,25);
   addHisto("Hemi_Mass","",        sample.Data(),  100, 0 ,500);

   // !! ------------------- end

   addHisto("track_TRACK_SIZE","EVT",           sample.Data(),5001,-0.5,5000.5);
   addHisto("track_nTracks", "TRK",           sample.Data(),150,0,150 );
   addHisto("track_nLostTracks", "TRK",           sample.Data(),70,0,70 );
   addHisto("track_ipc", "TRK",           sample.Data(), 50,0,5000);
   addHisto("track_dxyError", "TRK",           sample.Data(), 41,-20,20);
   addHisto("track_nHitPixel", "TRK",         sample.Data(), 15,0,15);
   addHisto("track_nHitTIB", "TRK",           sample.Data(), 15,0,15);
   addHisto("track_nHitTOB", "TRK",           sample.Data(), 15,0,15);
   addHisto("track_nHitTEC", "TRK",           sample.Data(), 15,0,15);
   addHisto("track_nHitPXB", "TRK",           sample.Data(), 15,0,15);
   addHisto("track_nHitPXF", "TRK",           sample.Data(), 15,0,15);
   addHisto("track_isHitPixel", "TRK",           sample.Data(), 3000,0,1500);

   addHisto("track_nLayers", "TRK",           sample.Data(), 80,0,20);
   addHisto("track_nLayersPixel", "TRK",           sample.Data(), 80,0,20);
   addHisto("track_region", "TRK",           sample.Data(),6,-0.5,5.5);
   addHisto("track_btag", "TRK",           sample.Data(), 22,-0.5,10.5);
   addHisto("track_Hemi", "TRK",           sample.Data(), 6,-0.5,5.5);
   addHisto("track_TRACK_SIZE","TRK",           sample.Data(),5001,-0.5,5000.5);



   addHisto("track_lost","TRK",                 sample.Data(),2,0,2);
   addHisto("track_dxy","TRK",                  sample.Data(),100,-50,50);
   addHisto("track_dz","TRK",                   sample.Data(),200,-100,100);
   addHisto("track_pt","TRK",                   sample.Data(),300, 0 ,300);
   addHisto("track_eta","TRK",                  sample.Data(),80,-4,4);
   addHisto("track_NChi2","TRK",                sample.Data(),5,0,5);
   addHisto("track_nhits","TRK",                sample.Data(),40,0,40);
   addHisto("track_iJet","TRK",                 sample.Data(),22,-2,20);
   addHisto("track_drSig","TRK",                sample.Data(),1000,0,1000);
   addHisto("track_dzSig","TRK",                sample.Data(),5000,0,5000);
   addHisto("track_ntrk10","TRK",               sample.Data(),100,0,100);
   addHisto("track_ntrk20","TRK",               sample.Data(),100,0,100);
   addHisto("track_ntrk30","TRK",               sample.Data(),100,0,100);
   addHisto("track_ntrk40","TRK",               sample.Data(),100,0,100);
   addHisto("track_track_Hemi_dR","TRK",              sample.Data(),50,0,5);
   addHisto("track_track_Hemi_dRmax","TRK",           sample.Data(),50,0,5);
   addHisto("track_MVAVal","TRK",                sample.Data(),101,-1,1);
   
   addHisto("track_Track_firstHit","TRK",sample.Data(),4000,0,4000);
   addHisto("track_Track_firstHit_x","TRK",sample.Data(),201,-100,100);// multiple of 401 
   addHisto("track_Track_firstHit_y","TRK",sample.Data(),201,-100,100);
   addHisto("track_Track_firstHit_r","TRK",sample.Data(),600,0,200);
   addHisto("track_Track_firstHit_z","TRK",sample.Data(),801,-200,200);
   addHisto2D("track_Track_firstHit_X_Y","TRK",sample.Data(),201,-100,100,201,-100,100);


// !! bgni
   addHisto("K0_mass","", sample.Data(),202,0.42,0.58);
   addHisto("K0_pt","", sample.Data(), 200, 0 ,200);

   addHisto("L0_mass","", sample.Data(),202,1.06,1.18);
   addHisto("L0_pt","", sample.Data(), 200, 0 ,200);

   addHisto("Reco_K0_mass","", sample.Data(),202,0.42,0.58);
   addHisto("Reco_K0_pt","", sample.Data(), 200, 0 ,200);
   addHisto("Reco_L0_mass","", sample.Data(),202,1.06,1.18);
   addHisto("Reco_L0_pt","", sample.Data(), 200, 0 ,200);


   addHisto("SecInt_mass","Selec", sample.Data(),20,0,2);
   addHisto("SecInt_drSig","Selec", sample.Data(),200,0,2000);                                                                                                                  //  addHisto("SecInt_dzSig","Selec", sample.Data(),100,0,500);                                                                                                          
   addHisto("SecInt_pt","Selec", sample.Data(),2000, 0 ,4000);
   addHisto("SecInt_dzSig","Selec", sample.Data(),200,0,2000);                                                              
   addHisto("SecInt_selec","Selec", sample.Data(),5,0,5);
   addHisto("SecInt_layer","Selec", sample.Data(),100,0,50);
addHisto("SecInt_r","Selec", sample.Data(), 400,0,200);
   addHisto("SecInt_z","Selec", sample.Data(), 800,-200,200);


   addHisto("SecInt_mass","TrackerMatched", sample.Data(),20,0,2);
   addHisto("SecInt_drSig","TrackerMatched", sample.Data(),200,0,2000);                                                                                                                  //  addHisto("SecInt_dzSig","Selec", sample.Data(),100,0,500);                                                                                                          
   addHisto("SecInt_pt","TrackerMatched", sample.Data(),2000, 0 ,4000);
   addHisto("SecInt_dzSig","TrackerMatched", sample.Data(),200,0,2000);                                                              
   addHisto("SecInt_selec","TrackerMatched", sample.Data(),5,0,5);
   addHisto("SecInt_layer","TrackerMatched", sample.Data(),100,0,50);
// !! end 

// !! bgin

   addHisto("Vtx_NChi2","",            sample.Data(),  15,0,15);
   addHisto("Vtx_nTrks","",           sample.Data(), 40,0,40);
   addHisto("Vtx_Mass","",         sample.Data(),  10,0,100);
   addHisto("Vtx_Dist","",        sample.Data(),  20,0,100);

   addHisto("SecVtx_NChi2","",            sample.Data(),  15,0,15);
   addHisto("SecVtx_nTrks","",           sample.Data(), 40,0,40);
   addHisto("SecVtx_Mass","",         sample.Data(),  150,0,1500);
   addHisto("SecVtx_Dist","",        sample.Data(),  100,0,100);

   addHisto("Vtx_Step","",            sample.Data(),  4,1,5);
   addHisto("Vtx_dR","",           sample.Data(), 50,0,5);
   addHisto("SecVtx_Step","",         sample.Data(),  4,1,5);
   addHisto("SecVtx_dR","",        sample.Data(), 50,0,5);


   addHisto("FinalVtx_nTrks","",            sample.Data(),  40,0,40);
   addHisto("FinalVtx_Step","",           sample.Data(), 4,1,5);
   addHisto("FinalVtx_Mass","",         sample.Data(),  100,0,500);
   addHisto("FinalVtx_HMass","",        sample.Data(),  100,0,500);

// !! end 

    // ----------------------------------------------//
}

//-------------------------------------------------------------
//instantiate and add
//first parameter is the variable name,
// second parameter is the selection step (like "afterleptsel")
//third parameter is the sample name (like "Z)
//others are TH1F binning
//creates one histograms per channel
//-------------------------------------------------------------
void DATAMCReader::addHisto(TString var, TString selstep, TString sample, int nbins, float min, float max){
 
  TString name =  sample+"_"+var+"_"+selstep;
  TH1F * thehisto = new TH1F(name,name,nbins,min,max);
//   std::cout << "   adding histo " << name << std::endl;
  thehisto->Sumw2();
  thehisto->SetOption("HIST");

  histo_list_.push_back(thehisto);
  histo_map_[name.Data()] = numb_histo;
  numb_histo++;
}

void DATAMCReader::addHisto2D(TString var, TString selstep, TString sample, int nxbins, float xmin, float xmax, int nybins, float ymin, float ymax){
 
  TString name =  sample+"_"+var+"_"+selstep;
  TH2F * thehisto = new TH2F(name,name,nxbins,xmin,xmax,nybins,ymin,ymax);
  // thehisto->Sumw2();
  thehisto->SetOption("COL");

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
void DATAMCReader::fillHisto( TString var, TString selstep, TString sample, float val, float weight){
  TString name = sample+"_"+var+"_"+selstep;


  if(histo_map_[name.Data()] == 0) {
    std::cout << "   WARNING trying to fill a non existing histograms " << std::endl;
    std::cout << "   please check the naming conventions " << std::endl;
    std::cout << "   histo name "  << name << std::endl;
  }else  histo_list_[histo_map_[name.Data()]]->Fill(val, weight);
  
}


void DATAMCReader::fillHisto2D( TString var, TString selstep, TString sample, float xval,float yval, float weight){
  TString name = sample+"_"+var+"_"+selstep;


  if(histo_map_2D_[name.Data()] == 0) {
    std::cout << "   WARNING trying to fill a non existing histograms " << std::endl;
    std::cout << "   please check the naming conventions " << std::endl;
    std::cout << "   histo name "  << name << std::endl;
  }else  histo_list_2D_[histo_map_2D_[name.Data()]]->Fill(xval,yval, weight);
  
}


void DATAMCReader::deleteHisto(){
   std::cout << __LINE__ << std::endl;

   /*for(unsigned int i=0; i<histo_list_mmm.size(); i++){
     
     delete  histo_list_mmm[i];
     delete  histo_list_mme[i];
     delete  histo_list_eem[i];
     delete  histo_list_eee[i];
     
     
   }*/
   std::cout << __LINE__ << std::endl;
  //delete TheTree;
   std::cout << __LINE__ << std::endl;
  
  
}
