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
#include "../../HistogramManager.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1F.h"
#include <TDirectory.h>
#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/MCWeights.h"
#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/DATAMCREADER/EMU/PDFSCALEVar.h"

// !! new 
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

// !! new 

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

void TreeABCDReader::Loop(bool isMC, TString Prod, TString sample, bool Signal, bool SS, bool FWD, int Year, float mean, bool DoubleMuon, int mixing, int Channel,  bool isPostAPV, TString thesystlist)
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


   bool firstinit = false;
   TString samplename = thesample;
     
     if( systlist== "") samplename = thesample;
     else                  samplename = thesample+"_"+systlist;
     
     
      firstinit = true;

   TFile * theoutputfile = new TFile( ("../../"+Prod+"/histofile_HT100_"+CHANNEL+"_"+ADD_Text+"_"+samplename+"_BDT100.root").Data() , "recreate");
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
   
  // !! -- Parameters -- !! //
   // // std::cout<<"Mixing : "<<mixing<<std::endl;
   // // std::cout<<"Msmuon : "<<MSmuon<<std::endl;
   bool BlindSR = false;
   float Nevent = 0;
Nevent = nentries;
   bool signal = Signal;
   int allevents = 0;
   double lowHpt = 30.;
   double LT = 0.;
   
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
  
  // !! new Vertex BDT 
   float mva_V_nTrks = 0 ;
    float mva_V_chi = 0;
    float mva_V_step = 0;
    float mva_V_r= 0;
    float mva_V_z = 0;
    float mva_V_MTW = 0;
    float mva_V_Mass =  0;
    float mva_H_Mass = 0;
    float mva_V_dist = 0;
    float mva_V_ntrk10 = 0;
    float mva_V_ntrk20 = 0 ;
    float mva_V_MeanDCA = 0;


      float mva_V1_nTrks = 0 ;
      float mva_V1_chi = 0;
      float mva_V1_step = 0;
      float mva_V1_r= 0;
      float mva_V1_z = 0;
      float mva_V1_MTW = 0;
      float mva_V1_Mass =  0;
      float mva_V1_dist = 0;
      float mva_V1_ntrk10 = 0;
      float mva_V1_ntrk20 = 0 ;
      float mva_V1_MeanDCA = 0;

      float mva_V2_nTrks = 0 ;
      float mva_V2_chi = 0;
      float mva_V2_step = 0;
      float mva_V2_r= 0;
      float mva_V2_z = 0;
      float mva_V2_MTW = 0;
      float mva_V2_Mass =  0;
      float mva_V2_dist = 0;
      float mva_V2_ntrk10 = 0;
      float mva_V2_ntrk20 = 0 ;
      float mva_V2_MeanDCA = 0;
    
    
   

    TMVA::Reader *readerVtx = new TMVA::Reader( "!Color:Silent" );

   // readerVtx->AddVariable( "mva_Vtx_nTrks",   &mva_V_nTrks);
   readerVtx->AddVariable( "mva_Vtx_NChi2",   &mva_V_chi);
   // readerVtx->AddVariable( "mva_Vtx_step",    &mva_V_step);
   readerVtx->AddVariable( "mva_Vtx_r",       &mva_V_r);
   readerVtx->AddVariable( "mva_Vtx_z",       &mva_V_z);
   readerVtx->AddVariable( "mva_Vtx_MTW",     &mva_V_MTW);
   readerVtx->AddVariable( "mva_Vtx_Mass",    &mva_V_Mass);
   // readerVtx->AddVariable( "mva_Hemi_Mass",   &mva_H_Mass);
   // readerVtx->AddVariable( "mva_Vtx_dist",     &mva_V_dist);
   // readerVtx->AddVariable( "mva_Vtx_ntrk10",  &mva_V_ntrk10);
   // readerVtx->AddVariable( "mva_Vtx_ntrk20",  &mva_V_ntrk20);
   readerVtx->AddVariable( "mva_Vtx_MeanDCA", &mva_V_MeanDCA);

   readerVtx->BookMVA( "BDTG", "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/BDT_2VTX_18_06_2025_TrainingvF_ctau100.xml"); // root 6.14/09, care compatiblity of versions for tmva

   // -----------------------
    TMVA::Reader *readerVtx_EVT = new TMVA::Reader( "!Color:Silent" );
   // readerVtx_EVT->AddVariable( "mva_Vtx1_step",    &mva_V1_step);
   readerVtx_EVT->AddVariable( "mva_Vtx1_NChi2",   &mva_V1_chi);
   readerVtx_EVT->AddVariable( "mva_Vtx1_r",       &mva_V1_r);
   readerVtx_EVT->AddVariable( "mva_Vtx1_z",       &mva_V1_z);
   readerVtx_EVT->AddVariable( "mva_Vtx1_MTW",     &mva_V1_MTW);
   readerVtx_EVT->AddVariable( "mva_Vtx1_Mass",    &mva_V1_Mass);
   readerVtx_EVT->AddVariable( "mva_Vtx1_MeanDCA", &mva_V1_MeanDCA);
   // readerVtx_EVT->AddVariable( "mva_Vtx2_step",    &mva_V2_step);
   readerVtx_EVT->AddVariable( "mva_Vtx2_NChi2",   &mva_V2_chi);
   readerVtx_EVT->AddVariable( "mva_Vtx2_r",       &mva_V2_r);
   readerVtx_EVT->AddVariable( "mva_Vtx2_z",       &mva_V2_z);
   readerVtx_EVT->AddVariable( "mva_Vtx2_MTW",     &mva_V2_MTW);
   readerVtx_EVT->AddVariable( "mva_Vtx2_Mass",    &mva_V2_Mass);
   readerVtx_EVT->AddVariable( "mva_Vtx2_MeanDCA", &mva_V2_MeanDCA);

   readerVtx_EVT->BookMVA( "BDTG", "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/BDT_EVT2VTX_18_06_2025_TrainingvF_ctau100_NoStep.xml"); // root 6.14/09, care compatiblity of versions for tmva

    
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


fL1_Reco_SF->Close();
fL1_ID_SF->Close();
fL1_ISO_SF->Close();
fL2_ISO_SF->Close();
fL2_Reco_SF->Close();
fL2_ID_SF->Close();
fL1L2_TRG_SFerr->Close();
fL1L2_TRG_SF->Close();

   TH2F *fh2DL1_ID_SF1;//
   TH2F *fh2DL2_ID_SF1;
   TH2F *fh2DL1_ID_SF1err;//
   TH2F *fh2DL2_ID_SF1err;


   TH2F *fh2DL1_ISO_SF1;
   TH2F *fh2DL2_ISO_SF1;
   TH2F *fh2DL1_ISO_SF1err;
   TH2F *fh2DL2_ISO_SF1err;
TH2F *fh2DL2_Reco_SF1err;
   TH2F *fh2DL1_Reco_SF1;
   TH2F *fh2DL2_Reco_SF1;

   
   TH2F *fhL1L2_TRG_SF; // x- axiss pt and y axis pt  
   TH2F *fhL1L2_TRG_SFerr; // x- axiss pt and y axis pt 


   int MuMu = Channel; // Emu = 0 ; SingleMuon = 1; DiMuon = 2
   

   if (MuMu == 0 && isMC)
      {
         if (Year == 2018)
            {
               // cout<< " 2018 year="<<endl;
               fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../../Scale_factors/Muon/Run2018_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Run2018_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

               fL2_Reco_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root");
               fL2_ISO_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2018.root");
               fL2_ID_SF= new TFile("../../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_EGM2D_U18.root");

               fL1L2_TRG_SF= new TFile("../../Scale_factors/Trig_2018/trigSF_2D.root");
               fL1L2_TRG_SFerr = new TFile("../../Scale_factors/Trig_2018/trigSF_2Derr.root");

            }
         if(Year == 2017)
            {
               //cout<<" year 2017"<<endl;
               fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Run2017_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
               
               fL2_Reco_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2017.root");
               fL2_ISO_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2017.root");
               fL2_ID_SF= new TFile("../../Scale_factors/Electron/egammaEffi.txt_EGM2D_Tight_UL17.root");

               // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2018_ULv2.root");
               fL1L2_TRG_SF= new TFile("../../Scale_factors/Trig_2018/trigSF_2D.root");
               fL1L2_TRG_SFerr = new TFile("../../Scale_factors/Trig_2017/trigSF_2Derr.root");
            }
         if(Year == 2016)
            {
               if (isPostAPV){
                  //cout<<" post 2016="<<endl;
                  // fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root");
                  // fL1_ID_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root");
                  // fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root");

                  fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

                  fL2_Reco_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root");
                  fL2_ISO_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root");
                  fL2_ID_SF= new TFile("../../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_postVFP_EGM2D.root");


                  // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2016postVFP_ULv2.root");
                  fL1L2_TRG_SF= new TFile("../../Scale_factors/Trig_2016/trigSF_2D.root");
                  fL1L2_TRG_SFerr = new TFile("../../Scale_factors/Trig_2016/trigSF_2Derr.root");

               }
               else{
                  // cout<<" pre  2016="<<endl;
                  // fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root");
                  // fL1_ID_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root");
                  // fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root");

                  fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL1_ISO_SF= new TFile("../../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
                  
                  fL2_Reco_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root");
                  fL2_ISO_SF= new TFile("../../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root");
                  fL2_ID_SF= new TFile("../../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_preVFP_EGM2D.root");


                  // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2016preVFP_ULv2.root");
                  fL1L2_TRG_SF= new TFile("../../Scale_factors/Trig_2016/trigSF_2D.root");
                  fL1L2_TRG_SFerr = new TFile("../../Scale_factors/Trig_2016/trigSF_2Derr.root");
               }
               
            }// else 2016

         fh2DL1_Reco_SF1 = (TH2F*)(fL1_Reco_SF->Get("NUM_TrackerMuons_DEN_genTracks"));
         fh2DL1_ID_SF1 = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt"));//
         
         fh2DL1_ISO_SF1 = (TH2F*)(fL1_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt"));
         //$$$$ 
         fh2DL2_Reco_SF1 = (TH2F*)(fL2_Reco_SF->Get("EGamma_SF2D"));
         fh2DL2_ISO_SF1 = (TH2F*)(fL2_ISO_SF->Get("EGamma_SF2D"));
         fh2DL2_ID_SF1 = (TH2F*)(fL2_ID_SF->Get("EGamma_SF2D"));
         
         fh2DL1_ID_SF1err  = (TH2F*)(fL1_ID_SF->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_combined_syst"));//
         fh2DL1_ISO_SF1err  = (TH2F*)(fL1_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt_combined_syst"));
         
         fh2DL2_ID_SF1err  = (TH2F*)(fL2_ID_SF->Get("statMC"));
         fh2DL2_ISO_SF1err  = (TH2F*)(fL2_ISO_SF->Get("statMC"));
         fh2DL2_Reco_SF1err = (TH2F*)(fL2_Reco_SF->Get("statMC"));


         // fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("h2D_SF_emu_lepABpt_FullError")); // x- axiss pt and y axis pt 
         TCanvas* c1 =  (TCanvas*)fL1L2_TRG_SF->Get("c1");
         TCanvas* c2 = (TCanvas*)fL1L2_TRG_SFerr->Get("c1");


         if (c1) {
            // Récupérer le premier TPad (ou cherche un spécifique si tu connais son nom)
            TPad *pad = (TPad*)c1->GetListOfPrimitives()->FindObject("pad1");

            if (pad) {
               // Récupérer le TH2F à l'intérieur du TPad
               fhL1L2_TRG_SF = (TH2F*)pad->GetListOfPrimitives()->FindObject("hSF");
               //  pad->ls();
               if (fhL1L2_TRG_SF) {
                     std::cout << "TH2F trouvé: " << fhL1L2_TRG_SF->GetName() << std::endl;
                     // h2->Draw("COLZ");  // Exemple d'utilisation
               } else {
                     std::cout << "Aucun TH2F trouvé dans le TPad." << std::endl;
               }
            } else {
               std::cout << "TPad introuvable dans le canvas." << std::endl;
            }
         } else {
            std::cout << "Canvas introuvable." << std::endl;
         }

         if (c2) {
            // Récupérer le premier TPad (ou cherche un spécifique si tu connais son nom)
            TPad *pad = (TPad*)c2->GetListOfPrimitives()->FindObject("pad1");
            // pad->ls();
            if (pad) {
               // Récupérer le TH2F à l'intérieur du TPad
               fhL1L2_TRG_SFerr = (TH2F*)pad->GetListOfPrimitives()->FindObject("hSFerr");

               if (fhL1L2_TRG_SFerr) {
                     std::cout << "TH2F 2 trouvé: " << fhL1L2_TRG_SFerr->GetName() << std::endl;
                     // h2->Draw("COLZ");  // Exemple d'utilisation
               } else {
                     std::cout << "Aucun TH2F trouvé dans le TPad2." << std::endl;
               }
            } else {
               std::cout << "TPad2 introuvable dans le canvas2." << std::endl;
            }
         } else {
            std::cout << "Canvas 2 introuvable." << std::endl;
         }

         // delete c1; delete c2;
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
               fh2DL1_ISO_SF1err  = (TH2F*)(fL1_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt_combined_syst"));
               fh2DL2_ISO_SF1err  = (TH2F*)(fL2_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt_combined_syst"));
         
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
               fh2DL1_ISO_SF1err  = (TH2F*)(fL1_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt_combined_syst"));
               fh2DL2_ISO_SF1err  = (TH2F*)(fL2_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt_combined_syst"));

               fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt")); // x- axiss pt and y axis pt 
               //$$$$ 

               fhL1L2_TRG_SFerr = (TH2F*)(fL1L2_TRG_SFerr->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt_combined_syst"));



            }
         if(Year == 2016)
            {
               if (isPostAPV){
                  //cout<<" post 2016="<<endl;
                  fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root");
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
                  fL1_Reco_SF =  new TFile("../../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root");
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
               fh2DL1_ISO_SF1err  = (TH2F*)(fL1_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt_combined_syst"));
               fh2DL2_ISO_SF1err  = (TH2F*)(fL2_ISO_SF->Get("NUM_MiniIsoTight_DEN_TightID_abseta_pt_combined_syst"));


               fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt")); // x- axiss pt and y axis pt 
               fhL1L2_TRG_SFerr = (TH2F*)(fL1L2_TRG_SFerr->Get("NUM_Trigger_DEN_MiniIsoTight_abseta_pt_combined_syst"));
               //$$$$ 
            }// else 2016
      }
      //---------------------------------------------------------------------------//
      // !! ----------------------------- End of Scale Factors---------------------------------//
      //---------------------------------------------------------------------------//
   double norm = 1.;
   if (!Signal && isMC)
      {
         TString hNorma = "hEvents_with_gen_wt";
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
         f1_DY->Close();
      }

  Long64_t nbytes = 0, nb = 0;
  cout<< "Line : "  << __LINE__ << " " << nentries2 << endl; 
  // Normalisation factor (XS)

   float XS = 1;
   float XS_up = 1.;
   float XS_down = 1.;

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
   float FracEffEvent = 1.;
   if (thesample.Contains("DYJetsToLL_M-10to50"))                     { XS = 22635; FracEffEvent = 1.;   }//ok : 15910 :https://doi.org/10.48550/arXiv.2402.08486
   if (thesample.Contains("DYJetsToLL_M-50"))                         { XS = 6225.4;  FracEffEvent = 1.;     }//5379.0;  ok  }https://doi.org/10.48550/arXiv.2402.08486
   if (thesample.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays")) { XS = 21.6;  FracEffEvent = 1.;  }//21.6 ok
   if (thesample.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     { XS = 21.6;  FracEffEvent = 1.;  }//21.6 ok
   if (thesample.Contains("TTJets_DiLept"))                          { XS = 53.07;     }//don't care
   if (thesample.Contains("TTJets"))                          { XS = 831.76;  FracEffEvent = 0.362;   }//don't care
   if (thesample.Contains("TTTo2L2Nu"))                              { XS = 88.5;  FracEffEvent = 0.993;    }// TOP 20 -006 : 88.9
   if (thesample.Contains("TTToSemiLeptonic") )                      { XS = 366.3;  FracEffEvent = 0.993;  }// TOP 20 -006 :366.6
      if (thesample.Contains("TTToHadronic"))                           { XS = 378.9;     }// TOP 20 -006 :377.6
   if (thesample.Contains("WWTo2L2Nu"))                              { XS = 11.09;  FracEffEvent = 0.997;    }
   if (thesample.Contains("WZTo2Q2L_mllmin4p0"))                     { XS = 6.535;  FracEffEvent = 0.583;   }
   if (thesample.Contains("ZZTo2Q2L_mllmin4p0"))                     { XS = 3.676;  FracEffEvent = 0.612;   }
   if (thesample.Contains("ttWJetsToLNu_5f_EWK"))                    { XS = 0.290;  FracEffEvent = 0.551;   } // not found on XSDB, no file on tier2...approximation
      //Took 0.868 pb (CMS-TOP-21-011)
   // as a starting point and then divided by 3 (lepton universality)
   if (thesample.Contains("TTZToLL_5f"))                             { XS = 0.253 ;FracEffEvent = 1.;   }// !!  0.253 Not found on XSDB => used ana.py macro : 0.05188  +- 2.437e-04 pb
   if (thesample.Contains("TTWW"))                                   { XS = 0.006992; FracEffEvent = 1.; }//found on XSDB
   if (thesample.Contains("ST_t-channel_antitop_5f_InclusiveDecays")) {XS = 80.95; FracEffEvent = 0.995;}
   if (thesample.Contains("ST_t-channel_top_5f_InclusiveDecays")) {XS = 136.02;FracEffEvent = 0.995;}
   //Signal
   if (thesample.Contains("smu200")) { XS = 0.01; MSmuon = 200;   }
   if (thesample.Contains("smu250")) { XS = 0.0045; MSmuon = 250; }
   if (thesample.Contains("smu300")) { XS = 0.002; MSmuon = 300; }
   if (thesample.Contains("smu350")) { XS = 0.001; MSmuon = 350; }
   if (thesample.Contains("smu400")) { XS = 0.0006; MSmuon = 400;}
   if (thesample.Contains("smu450")) { XS = 0.0004; MSmuon = 450;}
   if (thesample.Contains("smu500")) { XS = 0.00025;MSmuon = 500;}

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
// if (nentries > 10000000){nentries = 10000000; norm = 10000000.;}
   for (Long64_t jentry=0; jentry<nentries;jentry++) {//nentries2
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      allevents++;
      if ( allevents%1000000 == 0 )  std::cout << "events : " << allevents << std::endl;
      // cout<< " count "<<jentry<<endl;
      // if (jentry > 100000) break;
            
      if ( signal  ) 
         {
          if (  minitree_nLLP->at(0) != 2)  continue; // protection against rare wrong signal events    //
         } 
       fillHisto("StepEffi","",samplename, 0,1 );
      if (minitree_Mmumu->at(0) < 20) continue;
            float L1 = 1; 
      float L2 = 1;
      if ( !minitree_Filter->at(0) ) continue;//
      nFilterEvt++;
fillHisto("StepEffi","",samplename, 1,1 );
      if (minitree_njetNOmu->at(0) < 1) continue;
      nFilterJet++;
fillHisto("StepEffi","",samplename, 2,1 );
      if ( minitree_Hemi_pt->size() != 2 ) continue;
      nFilternHemi++;
fillHisto("StepEffi","",samplename, 3,1 );
      if ( minitree_Hemi_pt->at(0) < lowHpt || minitree_Hemi_pt->at(1) < lowHpt ) continue;
      nFilterHpt++;
fillHisto("StepEffi","",samplename, 4,1 );

      if (minitree_lepton_leadingpt->at(0) > minitree_lepton_leadingpt2->at(0)) {L1 = minitree_lepton_leadingpt->at(0); L2 = minitree_lepton_leadingpt2->at(0);}
      else {L1 = minitree_lepton_leadingpt2->at(0); L2 = minitree_lepton_leadingpt->at(0);}
      if ( L1 < 25 || L2 < 14) continue;

      // if ( minitree_Hemi_pt->size() != 2 ) continue;

      // if ( abs(minitree_Hemi_eta->at(0)) > 2.4 || abs(minitree_Hemi_eta->at(1)) > 2.4 ) continue;
      fillHisto("hData_Filter","",samplename, minitree_Filter->at(0),1 );
            
      // --------- Lepton SF ------//
      double Mu_SF=1, Ele_SF=1, Mu_SF2 = 1;
      double EMu_trig_SF1=1;
      double triggerSF = 1;
      double triggerSFerr=0;
      double  Mu_ID_SF1=1, Mu_ISO_SF1=1, Mu_Reco_SF1=1,Mu_ID_SF2=1, Mu_ISO_SF2=1, Mu_Reco_SF2=1, Ele_Reco_SF1=1, Ele_ID_SF1=1, Ele_Reco_SF2=1;//$$$$
      double  Mu_ID_SF1err=0, Mu_ISO_SF1err=0, Mu_ID_SF2err=0, Mu_ISO_SF2err=0,  Ele_ID_SF1err=0, Ele_Reco_SF1err = 0;//$$$$
      float mu1_pt=0;
      float mu1_eta = 0;
      float mu2_pt = 0;
	   float ele1_pt=0;
	   float ele1_eta=0;
      float SFele_pt_eta = 1.;
      float SFele_eta = 1.;
      float SFele_eta_err = 0.;
      float ele1_trigger_pt = 0;

      //  MuMu = Channel; // Emu = 0 ; SingleMuon = 1; DiMuon = 2
      // !!
      // std::cout<<" here 0 "<<std::endl;
      // std::cout << "minitree_lepton_leadingeta2.size() : " <<minitree_lepton_leadingeta2->size()<< std::endl;
      if (abs(minitree_lepton_leadingeta2->at(0))>1.5  && MuMu==0) 
         {
            SFele_eta = hEle_SF->GetBinContent(hEle_SF->GetXaxis()->FindBin(minitree_lepton_leadingeta2->at(0)));
            SFele_eta_err = hEle_SF->GetBinError(hEle_SF->GetXaxis()->FindBin(minitree_lepton_leadingeta2->at(0)));
         }
      // !! 
      if (!isMC)
         {
            mu1_pt = minitree_lepton_leadingpt->at(0);
            if (MuMu == 0)
               {
                  ele1_pt = minitree_lepton_leadingpt2->at(0);
               }
            else
               {
                  mu2_pt = minitree_lepton_leadingpt2->at(0);
               }
         }
            for (unsigned int iMuon = 0; iMuon < minitree_lepton_leadingpt->size(); iMuon++)
               {
                  if (minitree_lepton_leadingpt->at(0) == 0   )
                     {
                        Mu_Reco_SF1 = 1;
                        Mu_ID_SF1 = 1;	    
                        Mu_ISO_SF1 = 1 ;
                        //  std::cout<<"here2A : "<<std::endl;
                        continue;
                     }
                  if (MuMu == 0 && isMC)
                     {// leading letpn is muon by default
                        // leading letpn is muon by default , 
                           // !! ne pas mettre les valeurs par défaut à 1 quand on est en dehors des ranges de la SF => 15 1 et 119 etc
                              //  std::cout<<"here2B : "<<std::endl;
                              mu1_pt = minitree_lepton_leadingpt->at(iMuon);
                              float mu1_pt_reco = mu1_pt;
                              if(minitree_lepton_leadingpt->at(iMuon) <= 15. ) {  mu1_pt = 15.1 ;} 
                              if(minitree_lepton_leadingpt->at(iMuon) >= 120. ) { mu1_pt = 119. ;}
                              if (minitree_lepton_leadingpt->at(iMuon) >= 40.) {mu1_pt_reco = 39.9;}
                              Mu_Reco_SF1= fh2DL1_Reco_SF1->GetBinContent(fh2DL1_Reco_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_Reco_SF1->GetYaxis()->FindBin(mu1_pt_reco));
                              Mu_ID_SF1 = fh2DL1_ID_SF1->GetBinContent(fh2DL1_ID_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ID_SF1->GetYaxis()->FindBin(mu1_pt));	    
                              Mu_ID_SF1err = fh2DL1_ID_SF1err->GetBinContent(fh2DL1_ID_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ID_SF1err->GetYaxis()->FindBin(mu1_pt));
                              Mu_ISO_SF1 = fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(mu1_pt));//fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(mu1_pt));
                              Mu_ISO_SF1err = fh2DL1_ISO_SF1err->GetBinContent(fh2DL1_ISO_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1err->GetYaxis()->FindBin(mu1_pt));//fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(mu1_pt));

                     }
                  else if (MuMu > 0 && isMC) // $$$$$
                     {

                        mu1_pt = minitree_lepton_leadingpt->at(iMuon);
                        mu1_eta = minitree_lepton_leadingeta->at(iMuon);
                        mu2_pt = minitree_lepton_leadingpt2->at(iMuon);
                        float mu1_pt_reco = mu1_pt;
                        float mu2_pt_reco = mu2_pt;
                        if (minitree_lepton_leadingpt->at(iMuon) >= 40.) {mu1_pt_reco = 39.9;}
                        if (minitree_lepton_leadingpt2->at(iMuon) >= 40.) {mu2_pt_reco = 39.9;}

                        if(minitree_lepton_leadingpt->at(iMuon) <= 15. ) {  mu1_pt = 15.1 ;} 
                        if(minitree_lepton_leadingpt->at(iMuon) >= 120. ) { mu1_pt = 119. ;}
                        if(minitree_lepton_leadingpt2->at(iMuon) <= 15. ) {  mu2_pt = 15.1 ;} 
                        if(minitree_lepton_leadingpt2->at(iMuon) >= 120. ) { mu2_pt = 119. ;}
                        Mu_Reco_SF1= fh2DL1_Reco_SF1->GetBinContent(fh2DL1_Reco_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_Reco_SF1->GetYaxis()->FindBin(mu1_pt_reco));
                        Mu_ID_SF1 = fh2DL1_ID_SF1->GetBinContent(fh2DL1_ID_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ID_SF1->GetYaxis()->FindBin(mu1_pt));	    
                        Mu_ISO_SF1 = fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(mu1_pt));
                        Mu_ID_SF1err = fh2DL1_ID_SF1err->GetBinContent(fh2DL1_ID_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ID_SF1err->GetYaxis()->FindBin(mu1_pt));	    
                        Mu_ISO_SF1err = fh2DL1_ISO_SF1err->GetBinContent(fh2DL1_ISO_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta->at(iMuon))),fh2DL1_ISO_SF1err->GetYaxis()->FindBin(mu1_pt));
                        Mu_Reco_SF2= fh2DL1_Reco_SF1->GetBinContent(fh2DL1_Reco_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_Reco_SF1->GetYaxis()->FindBin(mu2_pt_reco));
                        Mu_ID_SF2 = fh2DL1_ID_SF1->GetBinContent(fh2DL1_ID_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_ID_SF1->GetYaxis()->FindBin(mu2_pt));	    
                        Mu_ISO_SF2 = fh2DL1_ISO_SF1->GetBinContent(fh2DL1_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_ISO_SF1->GetYaxis()->FindBin(mu2_pt));
                        Mu_ID_SF2err = fh2DL1_ID_SF1err->GetBinContent(fh2DL1_ID_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_ID_SF1err->GetYaxis()->FindBin(mu2_pt));	    
                        Mu_ISO_SF2err = fh2DL1_ISO_SF1err->GetBinContent(fh2DL1_ISO_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iMuon))),fh2DL1_ISO_SF1err->GetYaxis()->FindBin(mu2_pt));
                     }

               }
            if (MuMu == 0 && isMC)
               {
                     for(unsigned int iel = 0; iel <minitree_lepton_leadingpt2->size(); iel ++)//test    after trigger and lepton selection cut                                    
                        {
                           //  std::cout<<"here4 : "<<std::endl;
                          if (minitree_lepton_leadingpt2->at(0) == 0  )
                           {
                              Ele_Reco_SF1 = 1;
                              Ele_ID_SF1 = 1;	
                              Ele_ID_SF1err = 0;   
                              Ele_Reco_SF1err = 0; 
                              //  std::cout<<"here5A : "<<std::endl;
                              continue;
                           }

                           ele1_pt = minitree_lepton_leadingpt2->at(iel);
                            ele1_trigger_pt = minitree_lepton_leadingpt2->at(iel);
                           ele1_eta = fabs(minitree_lepton_leadingeta2->at(iel));

                           if(minitree_lepton_leadingpt2->at(iel) <= 10. ) {  ele1_pt = 11; }// ID, 10-500, eta -2.5 to 2.5, x axis = eta, y axis=pt   
                           if(minitree_lepton_leadingpt2->at(iel) >= 500. ) { ele1_pt = 499;}
                           if ( ele1_trigger_pt >= 200){ ele1_trigger_pt = 199;}
                           if ( ele1_trigger_pt <= 15){ ele1_trigger_pt = 15.1;}

                           if ( ele1_pt < 20)
                              {  
                                 Ele_Reco_SF1 = fh2DL2_ISO_SF1->GetBinContent(fh2DL2_ISO_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_ISO_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                                 Ele_Reco_SF1err = fh2DL2_ISO_SF1err->GetBinContent(fh2DL2_ISO_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_ISO_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                              }
                           else
                              {
                                 Ele_Reco_SF1 = fh2DL2_Reco_SF1->GetBinContent(fh2DL2_Reco_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_Reco_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                                 Ele_Reco_SF1err = fh2DL2_Reco_SF1err->GetBinContent(fh2DL2_Reco_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_Reco_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                              }
                           Ele_ID_SF1 = fh2DL2_ID_SF1->GetBinContent(fh2DL2_ID_SF1->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_ID_SF1->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                           Ele_ID_SF1err = fh2DL2_ID_SF1err->GetBinContent(fh2DL2_ID_SF1err->GetXaxis()->FindBin(fabs(minitree_lepton_leadingeta2->at(iel))),fh2DL2_ID_SF1err->GetYaxis()->FindBin(minitree_lepton_leadingpt2->at(iel)));
                        }
               }                       

         if (isMC && MuMu == 0) 
            {   

               bool skippy = false;

               if (ele1_trigger_pt < 15. || mu1_pt < 15. ||  ele1_trigger_pt > 200 || mu1_pt > 200){triggerSF=1.;triggerSFerr=0.;skippy = true;} // x and y axis 15 to 200
               if (!skippy)
                  {
                     // std::cout<<"here6A0: "<<std::endl;
                     triggerSF  = fhL1L2_TRG_SF->GetBinContent(fhL1L2_TRG_SF->GetXaxis()->FindBin(ele1_trigger_pt),fhL1L2_TRG_SF->GetYaxis()->FindBin(mu1_pt));
                     if (triggerSF == 0) {std::cout<<"ele1_pt : "<<ele1_pt<<std::endl;std::cout<<"mu1_pt : "<<mu1_pt<<std::endl;}
                     // std::cout<<"here6A1: "<<std::endl;
                     triggerSFerr = fhL1L2_TRG_SFerr->GetBinContent(fhL1L2_TRG_SFerr->GetXaxis()->FindBin(ele1_trigger_pt),fhL1L2_TRG_SFerr->GetYaxis()->FindBin(mu1_pt));
                     // std::cout<<"here6A2: "<<std::endl;
                  }
               //  std::cout<<"here6B: "<<std::endl;
            }
         else if (isMC && MuMu >0 ) 
            {
               triggerSF = fhL1L2_TRG_SF->GetBinContent(fhL1L2_TRG_SF->GetXaxis()->FindBin(fabs(mu1_eta)),fhL1L2_TRG_SF->GetYaxis()->FindBin(mu1_pt)); 
               triggerSFerr = fhL1L2_TRG_SFerr->GetBinContent(fhL1L2_TRG_SFerr->GetXaxis()->FindBin(fabs(mu1_eta)),fhL1L2_TRG_SFerr->GetYaxis()->FindBin(mu1_pt));
               if ((mu1_pt < 15.) ||  (mu1_pt > 120 )){triggerSF=1.;triggerSFerr=0.;}
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

      float ScaleWeight = 1;
      float PdfWeight = 1;

      std::vector<float>   PDFWeights;
      std::vector<float>   ScaleWeights;

      float top_pt_wt = 1.;
      float Prefweight = miniPrefweight->at(0);
      float PrefweightUp = 1.;
      float PrefweightDown = 1.;
      float PUweight = miniPUweight->at(0);
      float PUweightUp = 1.;
      float PUweightDown = 1.;
      float Roccor = 1;
      float RoccorUp = 1;
      float RoccorDown = 1;
      float TriggerSyst = triggerSF;
      Ele_SF =Ele_Reco_SF1*Ele_ID_SF1;
      Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * Mu_ISO_SF1;
      Mu_SF2 = Mu_Reco_SF2* Mu_ID_SF2 * Mu_ISO_SF2;
      if (Signal)
         {
            PUweightUp = miniPUweight_Up->at(0);
            PUweightDown = miniPUweight_Down->at(0);
         }

      NormFactorLumi =  XS*Lumi;  
      NormFactorLumiUp = XS*LumiUp;
      NormFactorLumiDown = XS*LumiDown;

      NormFactorXS =  XS*Lumi;  
      NormFactorXSUp = XS_up*Lumi;
      NormFactorXSDown = XS_down*Lumi;

      // rajouter l1 prefiring
      float NormFactorSYST = NormFactorLumi/norm;
if (isMC)
         {
         PUweight = miniPUweight->at(0);
         PUweightUp = miniPUweight_Up->at(0);
         PUweightDown = miniPUweight_Down->at(0);
         Prefweight = miniPrefweight->at(0);
         PrefweightUp = miniPrefweight_Up->at(0);
         PrefweightDown = miniPrefweight_Down->at(0);
         NormFactorSYST =  NormFactorLumi/norm;
}

         // It is supposed to run for several systematics sequentially, however the code oes not work when lloping on  different systematics ...
         // so we have to put the 0 instead of the index of the syst in the vector systlist

   if(thesample.Contains("TTTo2L2Nu") || thesample.Contains("TTToHadronic") ||  thesample.Contains("TTToSemiLeptonic")){
         top_pt_wt=minitree_genTop_Weight->at(0);
         // top_pt_wt=1;
      }
      else{
	      top_pt_wt=1;
      }


   if (!isMC)//<=> Data
      {
         NormFactorSYST = 1;                                                                                                                                                               
      }
   else{
     PUweight = miniPUweight->at(0);
         // !! regarer pour sf > 4 et fixer à 1
         PUweightUp = miniPUweight_Up->at(0);
         PUweightDown = miniPUweight_Down->at(0);
         Prefweight = miniPrefweight->at(0);
         PrefweightUp = miniPrefweight_Up->at(0);
         PrefweightDown = miniPrefweight_Down->at(0);
         //  std::cout<<"here8 : "<<std::endl;
         NormFactorSYST =  NormFactorLumi/norm;

         if (systlist.Contains("LumiUp"))
            {
               NormFactorSYST =  NormFactorLumiUp/norm;
            }
         else if (systlist.Contains("LumiDown"))
            {
               NormFactorSYST =  NormFactorLumiDown/norm;
            }
         else if (systlist.Contains("XSUp"))
            {
               NormFactorSYST =  NormFactorXSUp/norm;
            }
         else if (systlist.Contains("XSDown"))
            {
               NormFactorSYST =  NormFactorXSDown/norm;
            }
         else if (systlist.Contains("L1Up"))
            {
               Prefweight = PrefweightUp;
            }
         else if (systlist.Contains("L1Down"))
            {
               Prefweight = PrefweightDown;
            }
         else if (systlist.Contains("TriggerDown"))
            {
               TriggerSyst = triggerSF-triggerSFerr;
               
            }
         else if (systlist.Contains("TriggerUp"))
            {
               TriggerSyst = triggerSF+triggerSFerr;
            }
         else if (systlist.Contains("MuonIDUp"))
            {
               if (MuMu == 0)
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1+Mu_ID_SF1err) * Mu_ISO_SF1;
                     // Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1+Ele_ID_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1+Mu_ID_SF1err) * Mu_ISO_SF1;
                     Mu_SF2 = Mu_Reco_SF2* (Mu_ID_SF2+Mu_ID_SF2err) * Mu_ISO_SF2;
                  }
            }    
         else if (systlist.Contains("MuonIDDown"))
            {
               if (MuMu == 0)
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1-Mu_ID_SF1err) * Mu_ISO_SF1;
                     // Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1-Ele_ID_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1-Mu_ID_SF1err) * Mu_ISO_SF1;
                     Mu_SF2 = Mu_Reco_SF2* (Mu_ID_SF2-Mu_ID_SF2err) * Mu_ISO_SF2;
                  }
            }  
         else if (systlist.Contains("MuonISOUp"))
            {
               if (MuMu == 0)
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1+Mu_ISO_SF1err);
                     // Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1+Ele_Reco_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1+Mu_ISO_SF1err);
                     Mu_SF2 = Mu_Reco_SF2* Mu_ID_SF2 * (Mu_ISO_SF2+Mu_ISO_SF2err);
                  }
            } 
         else if (systlist.Contains("MuonISODown"))
            {
               if (MuMu == 0)
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1-Mu_ISO_SF1err);
                     // Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1-Ele_Reco_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1-Mu_ISO_SF1err);
                     Mu_SF2 = Mu_Reco_SF2* Mu_ID_SF2 * (Mu_ISO_SF2-Mu_ISO_SF2err);
                  }
            } 
                  else if (systlist.Contains("EleIDUp"))
            {
               if (MuMu == 0)
                  {
                     // Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1+abs(Mu_ID_SF1-Mu_ID_SF1err)) * Mu_ISO_SF1;
                     Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1+Ele_ID_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1+Mu_ID_SF1err) * Mu_ISO_SF1;
                     Mu_SF2 = Mu_Reco_SF2* (Mu_ID_SF2+Mu_ID_SF2err) * Mu_ISO_SF2;
                  }
            }    
         else if (systlist.Contains("EleIDDown"))
            {
               if (MuMu == 0)
                  {
                     // Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1- abs(Mu_ID_SF1-Mu_ID_SF1err)) * Mu_ISO_SF1;
                     Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1-Ele_ID_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* (Mu_ID_SF1-Mu_ID_SF1err) * Mu_ISO_SF1;
                     Mu_SF2 = Mu_Reco_SF2* (Mu_ID_SF2-Mu_ID_SF2err) * Mu_ISO_SF2;
                  }
            }  
         else if (systlist.Contains("EleISOUp"))
            {
               if (MuMu == 0)
                  {
                     // Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1+abs(Mu_ISO_SF1-Mu_ISO_SF1err));
                     Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1+Ele_Reco_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1+Mu_ISO_SF1err);
                     Mu_SF2 = Mu_Reco_SF2* Mu_ID_SF2 * (Mu_ISO_SF2+Mu_ISO_SF2err);
                  }
            } 
         else if (systlist.Contains("EleISODown"))
            {
               if (MuMu == 0)
                  {
                     // Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1-abs(Mu_ISO_SF1-Mu_ISO_SF1err));
                     Ele_SF=Ele_Reco_SF1*(Ele_ID_SF1-Ele_Reco_SF1err);
                  }
               else
                  {
                     Mu_SF = Mu_Reco_SF1* Mu_ID_SF1 * (Mu_ISO_SF1-Mu_ISO_SF1err);
                     Mu_SF2 = Mu_Reco_SF2* Mu_ID_SF2 * (Mu_ISO_SF2-Mu_ISO_SF2err);
                  }
            }
         else  if (systlist.Contains("PUUp"))
            {
               PUweight = PUweightUp;
            }
         else if (systlist.Contains("PUDown"))
            {
               PUweight = PUweightDown;
            }
         else if (systlist.Contains("SFEleDown"))
            {
               SFele_eta = SFele_eta - SFele_eta_err;
            }
         else if (systlist.Contains("SFEleUp"))
            {
               SFele_eta = SFele_eta + SFele_eta_err;
            }
         else if (systlist.Contains("TopPtUp"))
            {
               if(thesample.Contains("TTTo2L2Nu") || thesample.Contains("TTToHadronic") ||  thesample.Contains("TTToSemiLeptonic")){
                  top_pt_wt=minitree_genTop_Weight->at(0);
                  // top_pt_wt=1;
               }
               else{
                  top_pt_wt=1;
               }
            }
         else if (systlist.Contains("TopPtDown"))
            {
               top_pt_wt=1;
            }
         else if (systlist.Contains("PDFUp"))
            {  

               for (unsigned int i = 0; i < minitree_LHE_Weights->size(); i++)
                  {
                     if (i >=10 || i == 0) {PDFWeights.push_back(minitree_LHE_Weights->at(i));}
                     else  {ScaleWeights.push_back(minitree_LHE_Weights->at(i));}//if (i < 10 && i >=1)
                  }
               PDFVar PDFV(PDFWeights);
               ScaleVars ScaleV(ScaleWeights);
               PdfWeight = PDFV.GetPDFVarUp()/PDFV.GetPDFOriginal();
               GenWeight = GenWeight*PdfWeight;
               // std::cout<<"PdfWeight = "<<PdfWeight<<std::endl;
               // std::cout<<"PDFV.GetPDFVarUp() = "<<PDFV.GetPDFVarUp()<<std::endl;
               // std::cout<<"PDFV.GetPDFOriginal() = "<<PDFV.GetPDFOriginal()<<std::endl;

            }
         else if (systlist.Contains("PDFDown"))
            {
               for (unsigned int i = 0; i < minitree_LHE_Weights->size(); i++)
                  {
                     if (i >=10 || i == 0) {PDFWeights.push_back(minitree_LHE_Weights->at(i));}
                     else if (i < 10 && i >=1) {ScaleWeights.push_back(minitree_LHE_Weights->at(i));}
                  }
               
               PDFVar PDFV(PDFWeights);
               ScaleVars ScaleV(ScaleWeights);
               PdfWeight = PDFV.GetPDFVarDown()/PDFV.GetPDFOriginal();

               GenWeight = GenWeight*PdfWeight;
            }
         else if (systlist.Contains("ScaleUp"))
            {
               for (unsigned int i = 0; i < minitree_LHE_Weights->size(); i++)
                  {
                     if (i >=10 || i == 0) {PDFWeights.push_back(minitree_LHE_Weights->at(i));}
                     else if (i < 10 && i >=1) {ScaleWeights.push_back(minitree_LHE_Weights->at(i));}
                  }
               
               PDFVar PDFV(PDFWeights);
               ScaleVars ScaleV(ScaleWeights);
               ScaleWeight = ScaleV.GetScaleVarUp();
                if (std::isnan(ScaleWeight)) {ScaleWeight = 1;}
               GenWeight = GenWeight*ScaleWeight;
               // std::cout<<"ScaleWeight = "<<ScaleWeight<<std::endl;
            }
         else if (systlist.Contains("ScaleDown"))
            {
               for (unsigned int i = 0; i < minitree_LHE_Weights->size(); i++)
                  {
                     if (i >=10 || i == 0) {PDFWeights.push_back(minitree_LHE_Weights->at(i));}
                     else if (i < 10 && i >=1) {ScaleWeights.push_back(minitree_LHE_Weights->at(i));}
                  }
               
               PDFVar PDFV(PDFWeights);
               ScaleVars ScaleV(ScaleWeights);
               ScaleWeight = ScaleV.GetScaleVarDown();
               if (std::isnan(ScaleWeight)) {ScaleWeight = 1;}
               GenWeight = GenWeight*ScaleWeight;
            }
         else
            {
               NormFactorSYST =  NormFactorLumi/norm;
            }
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
            SFele_eta = 1;
            Pref_PU_gen_wt= (NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt) / FracEffEvent;
            NormFactor_mu =  (SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Mu_SF*TriggerSyst) / FracEffEvent;

            if (MuMu == 0)
               {
                  NormFactor =  (SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Mu_SF*Ele_SF*TriggerSyst) / FracEffEvent;
                  NormFactor_ele=  (SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Ele_SF*TriggerSyst) / FracEffEvent;
               }
            else
               {
                  NormFactor =  (SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Mu_SF*Mu_SF2*TriggerSyst) / FracEffEvent;
                  NormFactor_mu2 =  (SFele_eta*NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Mu_SF2*TriggerSyst) / FracEffEvent;
               }

            

         }
   if (Signal)
      {
            Pref_PU_gen_wt= NormFactorSYST*Prefweight*PUweight*top_pt_wt;
            NormFactor =  SFele_eta*NormFactorSYST*Prefweight*PUweight*top_pt_wt*Mu_SF*Mu_SF2*TriggerSyst;
            NormFactor_mu =  SFele_eta*NormFactorSYST*Prefweight*PUweight*top_pt_wt*Mu_SF*TriggerSyst;
            NormFactor_mu2 =  SFele_eta*NormFactorSYST*Prefweight*PUweight*top_pt_wt*Mu_SF2*TriggerSyst;
            NormFactor_ele=  SFele_eta*NormFactorSYST*Prefweight*PUweight*top_pt_wt*Ele_SF*TriggerSyst;
      }
      fillHisto("hData_Event_Weight","",samplename, NormFactor,1 );

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
      float Vtx_x = 0;
      float Vtx_y = 0;
      float Vtx_z= 0;
      float Vtx_r= 0;
      float Vtx_dR= 0;
      float Vtx_SumtrackWeight= 0;
      float Vtx_track_MeanDCA_d= 0;
      float Vtx_dist = 0;
      float Vtx_NChi = 0;
      float Vtx_ntrk10 = 0;
      float Vtx_step = 0;
      float highHpt = 100.;
      bool ping0 = false;
      bool ping1 = false;
      float dR0 = 0.;
      float dR1 = 0.;
      ////////////////// hemisphere pT
      float  hemi1_pt  = -1.;
      float  hemi2_pt  = -1.;
      ////////////////////////////////
      
      fillHisto("hData_Mmumu","",samplename, minitree_Mmumu->at(0),NormFactor );

      int FilterSample = -1;
      // if (showoutput) // std::cout<<"minitree_Mmumu->at(0) : "<<minitree_Mmumu->at(0)<<std::endl;
      if (MuMu == 2 || MuMu == 1)
         {
            if (DoubleMuon && minitree_trigger_doublelepton->at(0)  ) //&& !minitree_trigger_singlelepton->at(0)
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
      if (DoubleMuon && FilterSample != 2) continue;
      else if (!DoubleMuon && FilterSample != 1 && MuMu ==1) continue;
      else if (!DoubleMuon && FilterSample != 0 && MuMu ==0) continue;

         // std::cout<<" here 3 "<<std::endl;

      // !! --------------------------------------------------------//
      // !!--------------------------------------------------------//
      // !! --------------------------------------------------------//
      // !!                                                      // 
      // !!            ABCD regions code start here              //
      // !!                                                      //
      // !!                                                      //
      // !! -------------------------------------------------------//
      // !! --------------------------------------------------------//

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
      fillHisto("Hemisphere_leadingpt","", samplename, hemi_ptmax,NormFactor);
      fillHisto("Hemisphere_subleadingpt","", samplename, hemi_ptmin,NormFactor);

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
   float Vtx_x0 = minitree_Hemi_Vtx_x->at(0);
   float Vtx_x1 = minitree_Hemi_Vtx_x->at(1);
   float Vtx_y0 = minitree_Hemi_Vtx_y->at(0);
   float Vtx_y1 = minitree_Hemi_Vtx_y->at(1);
   float Vtx_z0 = minitree_Hemi_Vtx_z->at(0);
   float Vtx_z1 = minitree_Hemi_Vtx_z->at(1);
   float Vtx_r0 = minitree_Hemi_Vtx_r->at(0);
   float Vtx_r1 = minitree_Hemi_Vtx_r->at(1);
   float Vtx_dR0 = minitree_Hemi_Vtx_dR->at(0);
   float Vtx_dR1 = minitree_Hemi_Vtx_dR->at(1);
   float Vtx_SumtrackWeight0 = minitree_Hemi_Vtx_SumtrackWeight->at(0);
   float Vtx_SumtrackWeight1 = minitree_Hemi_Vtx_SumtrackWeight->at(1);
   float Vtx0_ntrk10 = minitree_Hemi_Vtx_BDT_ntrk10->at(0);
   float Vtx1_ntrk10 = minitree_Hemi_Vtx_BDT_ntrk10->at(1);

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
            Vtx_x0 = posx0;
            Vtx_y0 = posy0;
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
            Vtx_x1 = posx1;
            Vtx_y1 = posy1;
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
      float VtxChicut = 2.;
      float VtxChicut_low = 0.;
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
          Vtx_x = Vtx_x0;
          Vtx_y = Vtx_y0;
          Vtx_z = Vtx_z0;
          Vtx_r = Vtx_r0;
          Vtx_dR = Vtx_dR0;
          Vtx_SumtrackWeight = Vtx_SumtrackWeight0;
          Vtx_MeantrackWeight = Vtx_MeantrackWeight0;
          Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d0;
          Vtx_ntrk10 = Vtx0_ntrk10;
          Vtx_step = Vtx_step0;
         fillHisto("hData_Vtx_dist","Tight",samplename, Vtx_dist0,1);



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
          Vtx_x = Vtx_x0;
          Vtx_y = Vtx_y0;
          Vtx_z = Vtx_z0;
          Vtx_r = Vtx_r0;
          Vtx_dR = Vtx_dR0;
          Vtx_SumtrackWeight = Vtx_SumtrackWeight0;
          Vtx_MeantrackWeight = Vtx_MeantrackWeight0;
          Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d0;
          Vtx_ntrk10 = Vtx0_ntrk10;
          Vtx_step = Vtx_step0;
         fillHisto("hData_Vtx_dist","Loose",samplename, Vtx_dist0,1);

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
            if ( Vtx_x1 > Vtx_x )     Vtx_x           = Vtx_x1;
          if ( Vtx_y1 > Vtx_y )     Vtx_y           = Vtx_y1;
          if ( Vtx_r1 > Vtx_r )     Vtx_r           = Vtx_r1;
          if ( Vtx_z1 > Vtx_z )     Vtx_z           = Vtx_z1;
          if ( Vtx_dR1 > Vtx_dR )    Vtx_dR         = Vtx_dR1;
          if ( Vtx_SumtrackWeight1 > Vtx_SumtrackWeight )    Vtx_SumtrackWeight = Vtx_SumtrackWeight1;
          if ( Vtx_MeantrackWeight1 > Vtx_MeantrackWeight )    Vtx_MeantrackWeight = Vtx_MeantrackWeight1;
          if ( Vtx_track_MeanDCA_d1 < Vtx_track_MeanDCA_d ) Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d1;
          if ( Vtx1_ntrk10 > Vtx_ntrk10 ) Vtx_ntrk10 = Vtx1_ntrk10;
          if ( Vtx_step1 < Vtx_step && Vtx_step1 > 0) Vtx_step = Vtx_step1;
         fillHisto("hData_Vtx_dist","Tight",samplename, Vtx_dist1,1);

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
          if ( Vtx_x1 > Vtx_x )     Vtx_x           = Vtx_x1;
          if ( Vtx_y1 > Vtx_y )     Vtx_y           = Vtx_y1;
          if ( Vtx_r1 > Vtx_r )     Vtx_r           = Vtx_r1;
          if ( Vtx_z1 > Vtx_z )     Vtx_z           = Vtx_z1;
          if ( Vtx_dR1 > Vtx_dR )    Vtx_dR         = Vtx_dR1;
          if ( Vtx_SumtrackWeight1 > Vtx_SumtrackWeight )    Vtx_SumtrackWeight = Vtx_SumtrackWeight1;
          if ( Vtx_MeantrackWeight1 > Vtx_MeantrackWeight )    Vtx_MeantrackWeight = Vtx_MeantrackWeight1;
          if ( Vtx_track_MeanDCA_d1 < Vtx_track_MeanDCA_d ) Vtx_track_MeanDCA_d = Vtx_track_MeanDCA_d1;
          if ( Vtx1_ntrk10 > Vtx_ntrk10 ) Vtx_ntrk10 = Vtx1_ntrk10;
          if ( Vtx_step1 < Vtx_step && Vtx_step1 > 0) Vtx_step = Vtx_step1;
         fillHisto("hData_Vtx_dist","Loose",samplename, Vtx_dist1,1);

      }


      if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;// Tight TIght
      else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;// Tight and ???Loose or nothing

      if      ( isHemiVtx1Loose && isHemiVtx2Loose ) nVtxLoose = 2;// Loose and Loose
      else if ( isHemiVtx1Loose || isHemiVtx2Loose ) nVtxLoose = 1;// Loose and Tight or nothing

   nRecoTightVertex= nVtx;
   nRecoLooseVertex= nVtxLoose;
   if (nVtx+nVtxLoose == 2) {

         if (Vtx_SumtrackWeight0 <= VtxChicut) {

            fillHisto("VtxLowSTW_Hemipt","2Vtx",samplename,hemi1_pt,1);
         }

         if (Vtx_SumtrackWeight0 > VtxChicut) {

            fillHisto("VtxHighSTW_Hemipt","2Vtx",samplename,hemi1_pt,1);
         }
          
         if (Vtx_MeantrackWeight1 <= VtxChicut) {

            fillHisto("VtxLowSTW_Hemipt","2Vtx",samplename,hemi2_pt,1);
         }
         if (Vtx_MeantrackWeight1 > VtxChicut) {

            fillHisto("VtxHighSTW_Hemipt","2Vtx",samplename,hemi2_pt,1);
         }
   }
   float RATIOVTX = 0.;


   // !! New VTX BDT
   // mva_V_nTrks
   mva_V_chi =  Vtx_NChi0;
   mva_V_step =  Vtx_step0;
   mva_V_r = Vtx_r0;
   // mva_V_z = abs(Vtx_z0);
   mva_V_z = abs(Vtx_z0);
   mva_V_MTW = Vtx_SumtrackWeight0;
   mva_V_Mass = Vtx_Mass0;
   // mva_H_Mass);
   mva_V_dist = Vtx_dist0;
   // mva_V_ntrk10 = Vtx0_ntrk10;
   // mva_V_ntrk20);
   mva_V_MeanDCA = Vtx_track_MeanDCA_d0;
  double Vtx1_bdtVal = -10 ;
  double Vtx1_bdtVal_EVT = -10;
  if (isHemiVtx1 || isHemiVtx1Loose)
     {
      Vtx1_bdtVal = readerVtx->EvaluateMVA("BDTG");// values at -999 => thishappens if the Hemi_Mass and Hemi_Vtx_Mass are not definite

     }
   

   // mva_V_nTrks
   mva_V_chi =  Vtx_NChi1;
   mva_V_step = Vtx_step1;
   mva_V_r = Vtx_r1;
   // mva_V_z = abs(Vtx_z1);
   mva_V_z = abs(Vtx_z1);
   mva_V_MTW = Vtx_SumtrackWeight1;
   mva_V_Mass = Vtx_Mass1;
   // mva_H_Mass);
   mva_V_dist = Vtx_dist1;
   // mva_V_ntrk10 = Vtx1_ntrk10;
   // mva_V_ntrk20);
   mva_V_MeanDCA = Vtx_track_MeanDCA_d1;

   double Vtx2_bdtVal = -10;

  if (isHemiVtx2 || isHemiVtx2Loose)
     {
      Vtx2_bdtVal = readerVtx->EvaluateMVA("BDTG");// values at -999 => thishappens if the Hemi_Mass and Hemi_Vtx_Mass are not definite

     }
   
double Best_BDTVal = Vtx1_bdtVal;
   if (Vtx2_bdtVal > Vtx1_bdtVal)
      {
         Best_BDTVal = Vtx2_bdtVal;

      }

double SumBDTVal = -2;
double AveBDTVal = -2;
double SumSTW = -2;
double AveSTW = -2;

   if (nVtx+nVtxLoose == 2)
      {
         SumBDTVal = Vtx1_bdtVal + Vtx2_bdtVal;
         AveBDTVal = SumBDTVal/2.;
         SumSTW = Vtx_SumtrackWeight0 + Vtx_SumtrackWeight1;
         AveSTW = SumSTW/2.;
      }
   else if (nVtx+nVtxLoose == 1)
      {
         SumBDTVal = Best_BDTVal;
         AveBDTVal = Best_BDTVal;
      }
   else
      {
         SumBDTVal = -10;
         AveBDTVal = -10;
      }

   double Vtx_bdtVal = -10;
   double Vtx_bdtVal_NoSTW = -10;

      // mva_V_nTrks
   mva_V_chi =  Vtx_NChi;
   mva_V_step = Vtx_step;
   mva_V_r = Vtx_r;
   // mva_V_z = abs(Vtx_z);
   mva_V_z = abs(Vtx_z);

   mva_V_MTW = Vtx_SumtrackWeight;
   mva_V_Mass = VtxMass;
   // mva_H_Mass) = Vtx_HMass;
   mva_V_dist =  Vtx_dist;
   mva_V_ntrk10 = Vtx_ntrk10;
   // mva_V_ntrk20);
   mva_V_MeanDCA =  Vtx_track_MeanDCA_d ;

   if (nVtx+nVtxLoose == 1)
      {
         Vtx_bdtVal = readerVtx->EvaluateMVA("BDTG");
      }


  // !! New Vtx BDT at EVENT LEVEL taking as input the two vertices separately

   float VTX_BDTVal_EVT = -10;

   mva_V1_chi =  Vtx_NChi0;
   mva_V1_step =  Vtx_step0;
   mva_V1_r = Vtx_r0;
   mva_V1_z = abs(Vtx_z0);
   mva_V1_MTW = Vtx_SumtrackWeight0;
   mva_V1_Mass = Vtx_Mass0;
   // mva_V_dist = Vtx_dist0;
   mva_V1_MeanDCA = Vtx_track_MeanDCA_d0;

   mva_V2_chi =  Vtx_NChi1;
   mva_V2_step = Vtx_step1;
   mva_V2_r = Vtx_r1;
   mva_V2_z = abs(Vtx_z1);
   mva_V2_MTW = Vtx_SumtrackWeight1;
   mva_V2_Mass = Vtx_Mass1;
   // mva_V_dist = Vtx_dist1;
   mva_V2_MeanDCA = Vtx_track_MeanDCA_d1;

   if (nVtx+nVtxLoose == 2)
      {
         VTX_BDTVal_EVT = readerVtx_EVT->EvaluateMVA("BDTG");
      }

   // !! ---------------VtxRecoEffi ------------//
   if (Signal)
      {
         float rLLP1 = sqrt((minitree_Hemi_LLP_x->at(0)*minitree_Hemi_LLP_x->at(0))+(minitree_Hemi_LLP_y->at(0)*minitree_Hemi_LLP_y->at(0)));              
         float rLLP2 = sqrt((minitree_Hemi_LLP_x->at(1)*minitree_Hemi_LLP_x->at(1))+(minitree_Hemi_LLP_y->at(1)*minitree_Hemi_LLP_y->at(1)));              
         

         fillHisto("hSim_Hemi_Vtx_r","noSel",samplename,rLLP1,1);
         fillHisto("hSim_Hemi_Vtx_r","noSel",samplename,rLLP2,1);
         fillHisto("hSim_Hemi_Vtx_dist","noSel",samplename,minitree_Hemi_LLP_dist->at(0),1);
         fillHisto("hSim_Hemi_Vtx_dist","noSel",samplename,minitree_Hemi_LLP_dist->at(1),1);
      }

    if ( Vtx_NChi0>0 && Vtx_NChi0<10)//Reco Vtx criteria
      {//-->Goodrecovtx

         fillHisto("hData_Hemi_Vtx_r","Goodrecovtx",samplename,Vtx_r0,1);
         fillHisto("hData_Hemi_Vtx_dist","Goodrecovtx",samplename,Vtx_dist0,1);
        if (ping0)
          {

            fillHisto("hData_Hemi_Vtx_r","Ping",samplename,Vtx_r0,1);
            fillHisto("hData_Hemi_Vtx_dist","Ping",samplename,Vtx_dist0,1);
            //  Tightping
            if ( Vtx_step0 >=1 && Vtx_step0 <=2 ) fillHisto("hData_Hemi_Vtx_r","TightPing",samplename,Vtx_r0,1);
            if ( Vtx_step0 >=1 && Vtx_step0 <=2 ) fillHisto("hData_Hemi_Vtx_dist","TightPing",samplename,Vtx_dist0,1);
             // Looseping
            if ( Vtx_step0 >=3 && Vtx_step0 <=4 ) fillHisto("hData_Hemi_Vtx_r","LoosePing",samplename,Vtx_r0,1);
            if ( Vtx_step0 >=3 && Vtx_step0 <=4 ) fillHisto("hData_Hemi_Vtx_dist","LoosePing",samplename,Vtx_dist0,1);
          }
      }
   if ( Vtx_NChi1>0 && Vtx_NChi1<10)//Reco Vtx criteria
      {//-->Goodrecovtx
         fillHisto("hData_Hemi_Vtx_r","Goodrecovtx",samplename,Vtx_r1,1);
         fillHisto("hData_Hemi_Vtx_dist","Goodrecovtx",samplename,Vtx_dist1,1);
        if (ping1)
          {
            fillHisto("hData_Hemi_Vtx_r","Ping",samplename,Vtx_r1,1);
            fillHisto("hData_Hemi_Vtx_dist","Ping",samplename,Vtx_dist1,1);
            //  Tightping
            if ( Vtx_step1 >=1 && Vtx_step1 <=2 )fillHisto("hData_Hemi_Vtx_r","TightPing",samplename,Vtx_r1,1);
            // if ( step >=1 && step <=2 ) hData_Hemi_Vtx_eta_Tightping->Fill(eta);//to be changed
            if ( Vtx_step1>=1 && Vtx_step1 <=2 ) fillHisto("hData_Hemi_Vtx_dist","TightPing",samplename,Vtx_dist1,1);
             // Looseping
            if ( Vtx_step1 >=3 && Vtx_step1 <=4 )fillHisto("hData_Hemi_Vtx_r","LoosePing",samplename,Vtx_r1,1);
            // if ( step >=3 && step <=4 ) hData_Hemi_Vtx_eta_Looseping->Fill(eta);//to be changed
            if ( Vtx_step1 >=3 && Vtx_step1 <=4 ) fillHisto("hData_Hemi_Vtx_dist","LoosePing",samplename,Vtx_dist1,1);
          }
      }
   // !! ---------------VtxRecoEffi ------------//


   int nVtx1 = 0;
   int nVtx2 = 0;
   int nVtxLoose1 = 0;
   int nVtxLoose2 = 0;

   if (isHemiVtx1){nVtx1=1;}
   if (isHemiVtx2){nVtx2=1;}

   if (isHemiVtx1Loose){nVtxLoose1=1;}
   if (isHemiVtx2Loose){nVtxLoose2=1;}

   //-------------- //    
   float HemiAveragePt = (hemi_ptmax+hemi_ptmin)/2.;

   if ( (nRecoTightVertex + nRecoLooseVertex) == 2)
   {

      fillHisto("hData_VtxQualityTight_Hemipt","2Vtx",samplename,hemi_ptmax,nRecoTightVertex);
      fillHisto("hData_VtxQualityTight_Hemipt","2Vtx",samplename,hemi_ptmin,nRecoTightVertex);


      fillHisto("hData_VtxQualityTight_VtxBDT","2Vtx",samplename,Vtx1_bdtVal,nRecoTightVertex);
      fillHisto("hData_VtxQualityTight_VtxBDT","2Vtx",samplename,Vtx2_bdtVal,nRecoTightVertex);

      fillHisto("hData_VtxQualityLoose_VtxBDT","2Vtx",samplename,Vtx1_bdtVal,nRecoLooseVertex);
      fillHisto("hData_VtxQualityLoose_VtxBDT","2Vtx",samplename,Vtx2_bdtVal,nRecoLooseVertex);

      fillHisto2D("hData_Hemipt_VtxBDT","2Vtx",samplename,Vtx1_bdtVal,hemi1_pt,1);
      fillHisto2D("hData_Hemipt_VtxBDT","2Vtx",samplename,Vtx2_bdtVal,hemi2_pt,1);

      fillHisto("hData_VtxQualityLoose_Hemipt","2Vtx",samplename,hemi_ptmax,nRecoLooseVertex);
      fillHisto("hData_VtxQualityLoose_Hemipt","2Vtx",samplename,hemi_ptmin,nRecoLooseVertex);

   }

   if ((nVtx == 2 && nVtxLoose == 0)|| (nVtx == 0 && nVtxLoose == 2) || (nVtx == 1 && nVtxLoose == 1))
      {
          fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","NoSel",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,1);
          fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","NoSel",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,1);
          fillHisto("hData_Hemi_2VtxAll_NChi2","",samplename,Vtx_NChi0,1);
          fillHisto("hData_Hemi_2VtxAll_NChi2","",samplename,Vtx_NChi1,1);

      }

// !! $$ ££ ¤¤ NEW
      if (nVtx == 2)
         {
            fillHisto("hData_2Vtx_VTXMVA","Tight",samplename,Best_BDTVal,NormFactor );
            fillHisto("hData_2Vtx_VTXMVA_EVT","Tight",samplename,VTX_BDTVal_EVT,NormFactor );
         }
      else if (nVtx + nVtxLoose == 2)
         {
            fillHisto("hData_2Vtx_VTXMVA","TightLoose",samplename,Best_BDTVal,NormFactor );
            fillHisto("hData_2Vtx_VTXMVA_EVT","TightLoose",samplename,VTX_BDTVal_EVT,NormFactor );
         }
      else if (nVtxLoose == 2)
         {
            fillHisto("hData_2Vtx_VTXMVA","LooseLoose",samplename,Best_BDTVal,NormFactor );
            fillHisto("hData_2Vtx_VTXMVA_EVT","LooseLoose",samplename,VTX_BDTVal_EVT,NormFactor );   
         }
// !! $$ ££ ¤¤ NEW 2 

   if (Vtx_step0 >= 1 && Vtx_step0 <= 2 && Vtx_step1 >= 1 && Vtx_step1 <= 2 )
      {
            fillHisto("hData_2Vtx_VtxVTx_VTXMVA","TightTight",samplename,Best_BDTVal,NormFactor );
            fillHisto("hData_2Vtx_VtxVTx_VTXMVA_EVT","TightTight",samplename,VTX_BDTVal_EVT,NormFactor );  
      }
   else if (Vtx_step0 >= 1 && Vtx_step0 <= 2 && Vtx_step1 >= 3 && Vtx_step1 <= 4 )
      {
            fillHisto("hData_2Vtx_VtxVTx_VTXMVA","TightLoose",samplename,Best_BDTVal,NormFactor );
            fillHisto("hData_2Vtx_VtxVTx_VTXMVA_EVT","TightLoose",samplename,VTX_BDTVal_EVT,NormFactor );  
      }
   else if (Vtx_step0 >= 3 && Vtx_step0 <= 4 && Vtx_step1 >= 1 && Vtx_step1 <= 2 )
      {
            fillHisto("hData_2Vtx_VtxVTx_VTXMVA","LooseTight",samplename,Best_BDTVal,NormFactor );
            fillHisto("hData_2Vtx_VtxVTx_VTXMVA_EVT","LooseTight",samplename,VTX_BDTVal_EVT,NormFactor );  
      }
   else if (Vtx_step0 >= 3 && Vtx_step0 <= 4 && Vtx_step1 >= 3 && Vtx_step1 <= 4 )
      {
            fillHisto("hData_2Vtx_VtxVTx_VTXMVA","LooseLoose",samplename,Best_BDTVal,NormFactor );
            fillHisto("hData_2Vtx_VtxVTx_VTXMVA_EVT","LooseLoose",samplename,VTX_BDTVal_EVT,NormFactor );  
      }
// !! $$ ££ ¤¤ 
      //-----------------------------------------------------------//
      // ABCD using Hemipshere pt anf Tight+loose steps of vertexing 
      //-----------------------------------------------------------//

   isCutEvt = false;
    if (!BlindSR) {
         
      if (nVtx == 1 && nVtxLoose == 1 && isCutEvt) { // same as Loose since TL = LT
         fillHisto("hData_Hemi_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
         fillHisto("hData_Hemi_TLVtx_Mass","",samplename,VtxMass,NormFactor);

         fillHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,NormFactor);
         fillHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,NormFactor);
         fillHisto("hData_Hemi_TLVtxAll_Mass","",samplename,Vtx_Mass0,NormFactor);
         fillHisto("hData_Hemi_TLVtxAll_Mass","",samplename,Vtx_Mass1,NormFactor);
         fillHisto2D("hData_Hemi_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
         fillHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
         fillHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

      }
         
         if (nVtx == 2 && isCutEvt ) {
            fillHisto("hData_Hemi_2Vtx_Mmumu","",samplename,minitree_Mmumu->at(0),NormFactor);
            fillHisto("hData_Hemi_2Vtx_Mass","",samplename,VtxMass,NormFactor);
            fillHisto("hData_Hemi_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),NormFactor );
            fillHisto("hData_Hemi_2Vtx_MaxBDTvtx","",samplename, BDTvtx,NormFactor );
            fillHisto("hData_Hemi_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight,NormFactor );
            fillHisto2D("hData_Hemi_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);

            fillHisto("hData_Hemi_2VtxAll_Mass","",samplename,Vtx_Mass0,NormFactor);
            fillHisto("hData_Hemi_2VtxAll_Mass","",samplename,Vtx_Mass1,NormFactor);
            fillHisto("hData_Hemi_2VtxAll_BDTvtx","",samplename, BDTvtx1 ,NormFactor);
            fillHisto("hData_Hemi_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,NormFactor);
            fillHisto("hData_Hemi_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0,NormFactor );
            fillHisto("hData_Hemi_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1,NormFactor );
            fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
            fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

         }

      }
   else {

         if (hemi_ptmin > highHpt) isCutEvt = true;

      if (nVtx == 1 && nVtxLoose == 1 && isCutEvt ) { // same as Loose since TL = LT
         fillHisto("hData_Hemi_TLVtx_SumtrackWeight","",samplename,0,NormFactor);
         fillHisto("hData_Hemi_TLVtx_Mass","",samplename,0,NormFactor);

         fillHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",samplename,0,NormFactor);
         fillHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",samplename,0,NormFactor);
         fillHisto("hData_Hemi_TLVtxAll_Mass","",samplename,0,NormFactor);
         fillHisto("hData_Hemi_TLVtxAll_Mass","",samplename,0,NormFactor);
         fillHisto2D("hData_Hemi_TLVtx_STW_Ntrks","",samplename,0,0,NormFactor);
         fillHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",samplename,0,0,NormFactor);
         fillHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",samplename,0,0,NormFactor);

      }
         
         if (nVtx == 2 && isCutEvt ) {
            fillHisto("hData_Hemi_2Vtx_Mmumu","",samplename,minitree_Mmumu->at(0),NormFactor);
            fillHisto("hData_Hemi_2Vtx_Mass","",samplename,0,NormFactor);
            fillHisto("hData_Hemi_2Vtx_BDTevt","",samplename, 0,NormFactor );
            fillHisto("hData_Hemi_2Vtx_MaxBDTvtx","",samplename, 0,NormFactor );
            fillHisto("hData_Hemi_2Vtx_SumtrackWeight","",samplename, 0,NormFactor );
            fillHisto2D("hData_Hemi_2Vtx_STW_Ntrks","",samplename,0,Vtx_nTrks,NormFactor);

            fillHisto("hData_Hemi_2VtxAll_Mass","",samplename,0,NormFactor);
            fillHisto("hData_Hemi_2VtxAll_Mass","",samplename,0,NormFactor);
            fillHisto("hData_Hemi_2VtxAll_BDTvtx","",samplename, 0 ,NormFactor);
            fillHisto("hData_Hemi_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,NormFactor);
            fillHisto("hData_Hemi_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0,NormFactor );
            fillHisto("hData_Hemi_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1,NormFactor );
            fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
            fillHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);


         }

      }

      //--------CR : Low Pt --------//

      isCutEvt = false;


      if (Filter &&  ( (hemi1_pt >= lowHpt && hemi1_pt <highHpt && hemi2_pt >=highHpt) ||
             (hemi2_pt >= lowHpt && hemi2_pt <highHpt && hemi1_pt >=highHpt) ))//hemi_ptmin > lowHpt && hemi_ptmin <highHpt
         {

            if (nVtx == 1 && nVtxLoose == 1  ) { // same as Loose low pt since TL = LT
                 fillHisto("hData_CRlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                 fillHisto("hData_CRlowpt_TLVtx_Mass","",samplename,VtxMass,NormFactor);

                 fillHisto("hData_CRlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,NormFactor);
                 fillHisto("hData_CRlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,NormFactor);
                 fillHisto("hData_CRlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,NormFactor);
                 fillHisto("hData_CRlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,NormFactor);
                 fillHisto2D("hData_CRlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                 fillHisto2D("hData_CRlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                 fillHisto2D("hData_CRlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

            }

            if (nVtx == 2 )
               {
                  fillHisto("hData_CRlowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0) ,NormFactor);
                  fillHisto("hData_CRlowpt_2Vtx_Mass","",samplename,   VtxMass ,NormFactor);
                  fillHisto("hData_CRlowpt_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),NormFactor );
                  fillHisto("hData_CRlowpt_2Vtx_MaxBDTvtx","",samplename, BDTvtx ,NormFactor);
                  fillHisto("hData_CRlowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,NormFactor);
                  fillHisto2D("hData_CRlowpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  
                  fillHisto("hData_CRlowpt_2VtxAll_Mass","",samplename, Vtx_Mass0,NormFactor );
                  fillHisto("hData_CRlowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,NormFactor);
                  fillHisto("hData_CRlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,NormFactor );
                  fillHisto("hData_CRlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,NormFactor);
                  fillHisto("hData_CRlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,NormFactor);
                  fillHisto("hData_CRlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,NormFactor);         
                  fillHisto2D("hData_CRlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);


               }

         }

         //--------CR : Both vertices have Low Pt --------//
 
      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      // hemi1_pt = minitree_Hemi_pt->at(0);
      // hemi2_pt = minitree_Hemi_pt->at(1);
      if (Filter &&   hemi1_pt >= lowHpt && hemi1_pt <highHpt && lowHpt < hemi2_pt <highHpt )//hemi_ptmin > lowHpt && hemi_ptmin <highHpt
         {


            if (nVtx == 1 && nVtxLoose == 1  ) { // same as LooseLowLowpt since TL = LT
                  fillHisto("hData_CRlowlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                  fillHisto("hData_CRlowlowpt_TLVtx_Mass","",samplename,VtxMass,NormFactor);
                  fillHisto2D("hData_CRlowlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);

                  fillHisto("hData_CRlowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,NormFactor);
                  fillHisto("hData_CRlowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,NormFactor);
                  fillHisto2D("hData_CRlowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRlowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

            }

            if (nVtx == 2 )
               {
                  fillHisto("hData_CRlowlowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),NormFactor );
                  fillHisto("hData_CRlowlowpt_2Vtx_Mass","",samplename,   VtxMass ,NormFactor);
                  fillHisto("hData_CRlowlowpt_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),NormFactor );
                  fillHisto("hData_CRlowlowpt_2Vtx_MaxBDTvtx","",samplename, BDTvtx ,NormFactor);
                  fillHisto("hData_CRlowlowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,NormFactor);
                  fillHisto("hData_CRlowlowpt_2VtxAll_Mass","",samplename, Vtx_Mass0,NormFactor );
                  fillHisto("hData_CRlowlowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,NormFactor);
                  fillHisto("hData_CRlowlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,NormFactor );
                  fillHisto("hData_CRlowlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,NormFactor);
                  fillHisto("hData_CRlowlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,NormFactor);
                  fillHisto("hData_CRlowlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,NormFactor);

                  fillHisto2D("hData_CRlowlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRlowlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);


               }

         }
      //--------CR : Loose  --------//

      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if ( Filter ) 
         {

            if ( hemi_ptmin >highHpt ) isCutEvt = true;
            //$$ 


            if (nVtx == 1 && nVtxLoose == 1 && isCutEvt ) { 
                  fillHisto("hData_CRloose_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                  fillHisto("hData_CRloose_TLVtx_Mass","",samplename,VtxMass,NormFactor);
                  fillHisto2D("hData_CRloose_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto("hData_CRloose_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRloose_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRloose_TLVtxAll_Mass","",samplename,Vtx_Mass0,NormFactor);
                  fillHisto("hData_CRloose_TLVtxAll_Mass","",samplename,Vtx_Mass1,NormFactor);
                  fillHisto2D("hData_CRloose_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRloose_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);
            }

            if (nVtxLoose == 2 && isCutEvt)
               {
                  fillHisto("hData_CRloose_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0) ,NormFactor);
                  fillHisto("hData_CRloose_2Vtx_Mass","",samplename,   VtxMass,NormFactor );
                  fillHisto("hData_CRloose_2Vtx_BDTevt","",samplename, minitree_Evts_MVAval->at(0),NormFactor );
                  fillHisto("hData_CRloose_2Vtx_MaxBDTvtx","",samplename, BDTvtx,NormFactor );
                  fillHisto("hData_CRloose_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight,NormFactor );
                  fillHisto2D("hData_CRloose_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto2D("hData_CRloose_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRloose_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);
                  fillHisto("hData_CRloose_2VtxAll_Mass","",samplename, Vtx_Mass0 ,NormFactor);
                  fillHisto("hData_CRloose_2VtxAll_Mass","",samplename, Vtx_Mass1 ,NormFactor);
                  fillHisto("hData_CRloose_2VtxAll_BDTvtx","",samplename, BDTvtx1,NormFactor );
                  fillHisto("hData_CRloose_2VtxAll_BDTvtx","",samplename, BDTvtx2,NormFactor );           
                  fillHisto("hData_CRloose_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0,NormFactor );
                  fillHisto("hData_CRloose_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1,NormFactor );

               }  
         }

      //--------CR : Loose Lowpt  --------//

      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if (Filter && ( (hemi1_pt >= lowHpt && hemi1_pt <highHpt && hemi2_pt >=highHpt) ||
             (hemi2_pt >= lowHpt && hemi2_pt <highHpt && hemi1_pt >=highHpt) )  )//hemi_ptmin > lowHpt && hemi_ptmin <highHpt
         {

            if (nVtx == 1 && nVtxLoose == 1  ) {
                  fillHisto("hData_CRlooselowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                  fillHisto("hData_CRlooselowpt_TLVtx_Mass","",samplename,VtxMass,NormFactor);

                  fillHisto("hData_CRlooselowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlooselowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlooselowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,NormFactor);
                  fillHisto("hData_CRlooselowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,NormFactor);
                  fillHisto2D("hData_CRlooselowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto2D("hData_CRlooselowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRlooselowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

            }
            if (nVtxLoose == 2)
               {
                  fillHisto("hData_CRlooselowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),NormFactor );
                  fillHisto("hData_CRlooselowpt_2Vtx_Mass","",samplename,   VtxMass,NormFactor );

                  fillHisto("hData_CRlooselowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,NormFactor);

                  fillHisto2D("hData_CRlooselowpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto2D("hData_CRlooselowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRlooselowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

                  fillHisto("hData_CRlooselowpt_2VtxAll_Mass","",samplename, Vtx_Mass0 ,NormFactor);
                  fillHisto("hData_CRlooselowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,NormFactor);
                  fillHisto("hData_CRlooselowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,NormFactor );
                  fillHisto("hData_CRlooselowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,NormFactor);

                  fillHisto("hData_CRlooselowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,NormFactor);
                  fillHisto("hData_CRlooselowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,NormFactor);

               }

         }

               //--------CR : LooseLoose  LowLowpt  --------//

      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if (Filter &&  hemi1_pt >= lowHpt && hemi1_pt <highHpt &&  lowHpt <= hemi2_pt  &&  hemi2_pt <=highHpt   )//hemi_ptmin > lowHpt && hemi_ptmin <highHpt
         {


            if (nVtx == 1 && nVtxLoose == 1  ) {
                  fillHisto("hData_CRlooselowlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                  fillHisto("hData_CRlooselowlowpt_TLVtx_Mass","",samplename,VtxMass,NormFactor);


                  fillHisto2D("hData_CRlooselowlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto2D("hData_CRlooselowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRlooselowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

                  fillHisto("hData_CRlooselowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlooselowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlooselowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,NormFactor);
                  fillHisto("hData_CRlooselowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,NormFactor);
            }
            if (nVtxLoose == 2)
               {
                  fillHisto("hData_CRlooselowlowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),NormFactor );
                  fillHisto("hData_CRlooselowlowpt_2Vtx_Mass","",samplename,   VtxMass,NormFactor );

                  fillHisto("hData_CRlooselowlowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,NormFactor);

                  fillHisto2D("hData_CRlooselowlowpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto2D("hData_CRlooselowlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRlooselowlowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

                  fillHisto("hData_CRlooselowlowpt_2VtxAll_Mass","",samplename, Vtx_Mass0 ,NormFactor);
                  fillHisto("hData_CRlooselowlowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,NormFactor);
                  fillHisto("hData_CRlooselowlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,NormFactor );
                  fillHisto("hData_CRlooselowlowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,NormFactor);

                  fillHisto("hData_CRlooselowlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,NormFactor);
                  fillHisto("hData_CRlooselowlowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,NormFactor);


               }

         }


         // !! Tests 1 !! //
         // !! 
         // !! 
         // !! ---------------------------- !!//


         // !! rename la première catégorie !! //
         // Histograms for 2Vertices but TT and TL are combined together as well as HighLow and LowLow
         // Gives more stats in one box for the signal region
         // ends up in a new ABCD method and not ABCDEFGHI


      // this is new B, same as old C with LooseLooseLowLowpt
      if (Filter &&  hemi1_pt >= lowHpt && hemi1_pt <highHpt &&  lowHpt <= hemi2_pt  &&  hemi2_pt <=highHpt   )//hemi_ptmin > lowHpt && hemi_ptmin <highHpt
         {
            fillHisto("hData_CRlooselooselowlowpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),NormFactor );
            if ( nVtx == 0 && nVtxLoose == 2 ) {
                  // fillHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtx_Mass","",samplename,VtxMass,NormFactor);


                  fillHisto2D("hData_CRlooselooselowlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto2D("hData_CRlooselooselowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRlooselooselowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,NormFactor);

                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_MVA","",samplename,  minitree_Hemi_Vtx_MVAval_Tight->at(0),NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_MVA","",samplename,  minitree_Hemi_Vtx_MVAval_Tight->at(1),NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","",samplename,Vtx2_bdtVal,NormFactor );

                  // !! $$ ----------------------
                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","diffBin",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","6Bins",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","7Bins",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","8Bins",samplename,Best_BDTVal,NormFactor );

                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","Sum",samplename,SumBDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","Ave",samplename,AveBDTVal,NormFactor );

                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_STW","Sum",samplename,SumSTW,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_STW","Ave",samplename,AveSTW,NormFactor );

                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_STW","SumdiffBin",samplename,SumSTW,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2Vtx_STW","AvediffBin",samplename,AveSTW,NormFactor );

                  // !! $$ ----------------------


                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","6Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","6Bins",samplename,Vtx2_bdtVal,NormFactor );

                  // !! ----------------------- VtxBDTVariables

                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NChi2","",samplename,Vtx_NChi0,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NChi2","",samplename,Vtx_NChi1,NormFactor );

                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_z","",samplename,Vtx_z0,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_z","",samplename,Vtx_z1,NormFactor );

                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_dist","",samplename,Vtx_dist0,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_dist","",samplename,Vtx_dist1,NormFactor );

                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_ntrk10","",samplename,Vtx0_ntrk10,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_ntrk10","",samplename,Vtx1_ntrk10,NormFactor );

                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_MeanDCA","",samplename,Vtx_track_MeanDCA_d0,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_MeanDCA","",samplename,Vtx_track_MeanDCA_d1,NormFactor );

                  // !! --------------------

               // !! Control Plots

               fillHisto2D("TrackerMap","BestVtx_A",samplename,Vtx_x,Vtx_y,1);//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","rz_BestVtx_A",samplename,Vtx_z,Vtx_r,1);//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","BothVtx_A",samplename,Vtx_x0,Vtx_y0,1);//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","BothVtx_A",samplename,Vtx_x1,Vtx_y1,1 );//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","rz_BothVtx_A",samplename,Vtx_z0,Vtx_r0,1);//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","rz_BothVtx_A",samplename,Vtx_z1,Vtx_r1,1);//500,-25.,25.,500,-25.,25.

               fillHisto("Step","BestVtx_A",samplename,Vtx_step,NormFactor );
               fillHisto("r","BestVtx_A",samplename,Vtx_r,NormFactor );
               fillHisto("z","BestVtx_A",samplename,abs(Vtx_z),NormFactor );
               fillHisto("STW","BestVtx_A",samplename,Vtx_SumtrackWeight,NormFactor );
               fillHisto("Mass","BestVtx_A",samplename,VtxMass,NormFactor );
               fillHisto("ntrk10","BestVtx_A",samplename,Vtx_ntrk10,NormFactor );
               fillHisto("MeanDCA","BestVtx_A",samplename,Vtx_track_MeanDCA_d,NormFactor );
               // !! -------------

                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","7Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","7Bins",samplename,Vtx2_bdtVal,NormFactor );

                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","8Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","8Bins",samplename,Vtx2_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowlowpt_dist_VtxVtx","",samplename,  Vtx_Vtx_dist,NormFactor );


                  fillHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","diffBin",samplename,Vtx_SumtrackWeight,NormFactor);

                                    if (Vtx_SumtrackWeight0 > 8 ) { Vtx_SumtrackWeight0 = 8.5;}
                  if (Vtx_SumtrackWeight1 > 8 ) { Vtx_SumtrackWeight1 = 8.5;}
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight,NormFactor);
                                    if (Vtx_SumtrackWeight0 > 7 ) { Vtx_SumtrackWeight0 = 7.5;}
                  if (Vtx_SumtrackWeight1 > 7 ) { Vtx_SumtrackWeight1 = 7.5;}
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight,NormFactor);
                  if (Vtx_SumtrackWeight0 > 6 ) { Vtx_SumtrackWeight0 = 6.5;}
                  if (Vtx_SumtrackWeight1 > 6 ) { Vtx_SumtrackWeight1 = 6.5;}
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight,NormFactor);
            }
         }


//--------CR : TightLoose + ITghtTight for low lowpt --------//
// Old A and B combined to form new A
      isCutEvt = false;
      // BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if (Filter &&  hemi1_pt >= lowHpt && hemi1_pt <highHpt &&  lowHpt <= hemi2_pt  &&  hemi2_pt <=highHpt )//hemi_ptmin > lowHpt && hemi_ptmin <highHpt
         {

            if ( (nVtx == 1 && nVtxLoose == 1 ) || (nVtx == 2 && nVtxLoose == 0) ){
                  
                  fillHisto("hData_CRtightlowlowpt_TLVtx_Mass","",samplename,VtxMass,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","",samplename,Vtx_SumtrackWeight1,NormFactor);

                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass0,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_Mass","",samplename,Vtx_Mass1,NormFactor);
                  fillHisto2D("hData_CRtightlowlowpt_TLVtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto2D("hData_CRtightlowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRtightlowlowpt_TLVtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);


                  fillHisto("hData_CRtightlowlowpt_2VtxAll_MVA","",samplename,  minitree_Hemi_Vtx_MVAval_Tight->at(0),NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_MVA","",samplename,  minitree_Hemi_Vtx_MVAval_Tight->at(1),NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","",samplename,Vtx2_bdtVal,NormFactor );

                  // !! $$ ----------------------
                  fillHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","diffBin",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","6Bins",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","7Bins",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","8Bins",samplename,Best_BDTVal,NormFactor );

                  fillHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","Sum",samplename,SumBDTVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","Ave",samplename,AveBDTVal,NormFactor );

                  fillHisto("hData_CRtightlowlowpt_2Vtx_STW","Sum",samplename,SumSTW,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2Vtx_STW","Ave",samplename,AveSTW,NormFactor );

                  fillHisto("hData_CRtightlowlowpt_2Vtx_STW","SumdiffBin",samplename,SumSTW,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2Vtx_STW","AvediffBin",samplename,AveSTW,NormFactor );
                  // !! $$ ----------------------
                  // !! Control Plots


                  fillHisto2D("TrackerMap","BestVtx_B",samplename,Vtx_x,Vtx_y,1);//500,-25.,25.,500,-25.,25.
                  fillHisto2D("TrackerMap","rz_BestVtx_B",samplename,Vtx_z,Vtx_r,1);//500,-25.,25.,500,-25.,25.
                  fillHisto2D("TrackerMap","BothVtx_B",samplename,Vtx_x0,Vtx_y0,1);//500,-25.,25.,500,-25.,25.
                  fillHisto2D("TrackerMap","BothVtx_B",samplename,Vtx_x1,Vtx_y1,1 );//500,-25.,25.,500,-25.,25.

                  fillHisto2D("TrackerMap","rz_BothVtx_B",samplename,Vtx_z0,Vtx_r0,1);//500,-25.,25.,500,-25.,25.
                  fillHisto2D("TrackerMap","rz_BothVtx_B",samplename,Vtx_z1,Vtx_r1,1);//500,-25.,25.,500,-25.,25.

                  fillHisto("Step","BestVtx_B",samplename,Vtx_step,NormFactor );
                  fillHisto("r","BestVtx_B",samplename,Vtx_r,NormFactor );
                  fillHisto("z","BestVtx_B",samplename,abs(Vtx_z),NormFactor );
                  fillHisto("STW","BestVtx_B",samplename,Vtx_SumtrackWeight,NormFactor );
                  fillHisto("Mass","BestVtx_B",samplename,VtxMass,NormFactor );
                  fillHisto("ntrk10","BestVtx_B",samplename,Vtx_ntrk10,NormFactor );
                  fillHisto("MeanDCA","BestVtx_B",samplename,Vtx_track_MeanDCA_d,NormFactor );
         // !! ------------

                  // !! ----------------------- VtxBDT variables

                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_NChi2","",samplename,Vtx_NChi0,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_NChi2","",samplename,Vtx_NChi1,NormFactor );

                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_z","",samplename,Vtx_z0,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_z","",samplename,Vtx_z1,NormFactor );

                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_dist","",samplename,Vtx_dist0,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_dist","",samplename,Vtx_dist1,NormFactor );

                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_ntrk10","",samplename,Vtx0_ntrk10,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_ntrk10","",samplename,Vtx1_ntrk10,NormFactor );

                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_MeanDCA","",samplename,Vtx_track_MeanDCA_d0,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_MeanDCA","",samplename,Vtx_track_MeanDCA_d1,NormFactor );

                  // !! --------------------


                  fillHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","6Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","6Bins",samplename,Vtx2_bdtVal,NormFactor );

                  fillHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","7Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","7Bins",samplename,Vtx2_bdtVal,NormFactor );

                  fillHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","8Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","8Bins",samplename,Vtx2_bdtVal,NormFactor );
                  fillHisto("hData_CRtightlowlowpt_dist_VtxVtx","",samplename,  Vtx_Vtx_dist,NormFactor );


                  fillHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","diffBin",samplename,Vtx_SumtrackWeight,NormFactor);

                  if (Vtx_SumtrackWeight0 > 8 ) { Vtx_SumtrackWeight0 = 8.5;}
                  if (Vtx_SumtrackWeight1 > 8 ) { Vtx_SumtrackWeight1 = 8.5;}
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight,NormFactor);


                  if (Vtx_SumtrackWeight0 > 7 ) { Vtx_SumtrackWeight0 = 7.5;}
                  if (Vtx_SumtrackWeight1 > 7 ) { Vtx_SumtrackWeight1 = 7.5;}
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight,NormFactor);


                  if (Vtx_SumtrackWeight0 > 6 ) { Vtx_SumtrackWeight0 = 6.5;}
                  if (Vtx_SumtrackWeight1 > 6 ) { Vtx_SumtrackWeight1 = 6.5;}
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight,NormFactor);
            }
         }


//--------CR : LooseLoose for highhighpt or highlowpt --------//
// Old F and I combined to form new D


      if (Filter &&  ((hemi1_pt >= highHpt && hemi2_pt >= highHpt) || ( (hemi1_pt >= lowHpt && hemi1_pt <highHpt && hemi2_pt >=highHpt) ||
             (hemi2_pt >= lowHpt && hemi2_pt <highHpt && hemi1_pt >=highHpt) )  ))//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {

            if (nVtxLoose == 2 && nVtx == 0)
               {
                  fillHisto("hData_CRlooselooselowpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_Mass","",samplename,   VtxMass,NormFactor );

                  fillHisto("hData_CRlooselooselowpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,NormFactor);

                  fillHisto2D("hData_CRlooselooselowpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto2D("hData_CRlooselooselowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRlooselooselowpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_Mass","",samplename, Vtx_Mass0 ,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,NormFactor);

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,NormFactor);


                  fillHisto("hData_CRlooselooselowpt_2VtxAll_MVA","",samplename,  minitree_Hemi_Vtx_MVAval_Tight->at(0),NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_MVA","",samplename,  minitree_Hemi_Vtx_MVAval_Tight->at(1),NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","",samplename,Vtx2_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","6Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","6Bins",samplename,Vtx2_bdtVal,NormFactor );

                  // !! $$ ----------------------
                  fillHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","diffBin",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","6Bins",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","7Bins",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","8Bins",samplename,Best_BDTVal,NormFactor );

                  fillHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","Sum",samplename,SumBDTVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","Ave",samplename,AveBDTVal,NormFactor );

                  fillHisto("hData_CRlooselooselowpt_2Vtx_STW","Sum",samplename,SumSTW,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_STW","Ave",samplename,AveSTW,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_STW","SumdiffBin",samplename,SumSTW,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2Vtx_STW","AvediffBin",samplename,AveSTW,NormFactor );
                  // !! $$ ----------------------


                  // !! ----------------------- VtxBDT variables

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_NChi2","",samplename,Vtx_NChi0,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_NChi2","",samplename,Vtx_NChi1,NormFactor );

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_z","",samplename,Vtx_z0,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_z","",samplename,Vtx_z1,NormFactor );

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_dist","",samplename,Vtx_dist0,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_dist","",samplename,Vtx_dist1,NormFactor );

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_ntrk10","",samplename,Vtx0_ntrk10,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_ntrk10","",samplename,Vtx1_ntrk10,NormFactor );

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_MeanDCA","",samplename,Vtx_track_MeanDCA_d0,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_MeanDCA","",samplename,Vtx_track_MeanDCA_d1,NormFactor );

                  // !! --------------------
               // !! Control Plots

               fillHisto2D("TrackerMap","BestVtx_D",samplename,Vtx_x,Vtx_y,1);//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","rz_BestVtx_D",samplename,Vtx_z,Vtx_r,1);//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","BothVtx_D",samplename,Vtx_x0,Vtx_y0,1);//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","BothVtx_D",samplename,Vtx_x1,Vtx_y1,1 );//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","rz_BothVtx_D",samplename,Vtx_z0,Vtx_r0,1);//500,-25.,25.,500,-25.,25.
               fillHisto2D("TrackerMap","rz_BothVtx_D",samplename,Vtx_z1,Vtx_r1,1);//500,-25.,25.,500,-25.,25.


               fillHisto("Step","BestVtx_D",samplename,Vtx_step,NormFactor );
               fillHisto("r","BestVtx_D",samplename,Vtx_r,NormFactor );
               fillHisto("z","BestVtx_D",samplename,abs(Vtx_z),NormFactor );
               fillHisto("STW","BestVtx_D",samplename,Vtx_SumtrackWeight,NormFactor );
               fillHisto("Mass","BestVtx_D",samplename,VtxMass,NormFactor );
               fillHisto("ntrk10","BestVtx_D",samplename,Vtx_ntrk10,NormFactor );
               fillHisto("MeanDCA","BestVtx_D",samplename,Vtx_track_MeanDCA_d,NormFactor );
               // !! -------------


                  fillHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","7Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","7Bins",samplename,Vtx2_bdtVal,NormFactor );

                  fillHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","8Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","8Bins",samplename,Vtx2_bdtVal,NormFactor );

                  fillHisto("hData_CRlooselooselowpt_dist_VtxVtx","",samplename,  Vtx_Vtx_dist,NormFactor );


                  fillHisto("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","diffBin",samplename,Vtx_SumtrackWeight,NormFactor);


                  if (Vtx_SumtrackWeight0 > 8 ) { Vtx_SumtrackWeight0 = 8.5;}
                  if (Vtx_SumtrackWeight1 > 8 ) { Vtx_SumtrackWeight1 = 8.5;}
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight,NormFactor);


                  if (Vtx_SumtrackWeight0 > 7 ) { Vtx_SumtrackWeight0 = 7.5;}
                  if (Vtx_SumtrackWeight1 > 7 ) { Vtx_SumtrackWeight1 = 7.5;}
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight,NormFactor);


                  if (Vtx_SumtrackWeight0 > 6 ) { Vtx_SumtrackWeight0 = 6.5;}
                  if (Vtx_SumtrackWeight1 > 6 ) { Vtx_SumtrackWeight1 = 6.5;}
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight,NormFactor);
               }

         }

      //--------SR : TIghtLoose+ TIghtTight for highhighpt or highlowpt --------//
// Old F and I combined to form new C
if ( (nVtx +  nVtxLoose == 2))
   {
      fillHisto("StepEffi","",samplename, 5,1 );// at least 1 tight vtx
   }


if ((nVtx == 1 && nVtxLoose == 1) || (nVtx == 2 && nVtxLoose == 0))
   {
      fillHisto("StepEffi","",samplename, 6,1 );// at least 1 tight vtx
   }



      if (Filter &&  ((hemi1_pt >= highHpt && hemi2_pt >= highHpt) || ( (hemi1_pt >= lowHpt && hemi1_pt <highHpt && hemi2_pt >=highHpt) ||
             (hemi2_pt >= lowHpt && hemi2_pt <highHpt && hemi1_pt >=highHpt) )  ))//hemi_ptmin > lowHpt && hemi_ptmin < 80.
         {
            fillHisto("hData_CRtighthighpt_BDTevt","",samplename, minitree_Evts_MVAval->at(0),NormFactor );

            if (nVtx == 2 || (nVtx == 1 && nVtxLoose == 1 ) ){
   
                  fillHisto("hData_CRtighthighpt_2Vtx_Mmumu","",samplename,  minitree_Mmumu->at(0),NormFactor );
                  fillHisto("hData_CRtighthighpt_2Vtx_Mass","",samplename,   VtxMass,NormFactor );

                  fillHisto("hData_CRtighthighpt_2Vtx_SumtrackWeight","",samplename, Vtx_SumtrackWeight ,NormFactor);

                  fillHisto2D("hData_CRtighthighpt_2Vtx_STW_Ntrks","",samplename,Vtx_SumtrackWeight,Vtx_nTrks,NormFactor);
                  fillHisto2D("hData_CRtighthighpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight0,Vtx_nTrks0,NormFactor);
                  fillHisto2D("hData_CRtighthighpt_2VtxAll_STW_Ntrks","",samplename,Vtx_SumtrackWeight1,Vtx_nTrks1,NormFactor);

                  fillHisto("hData_CRtighthighpt_2VtxAll_Mass","",samplename, Vtx_Mass0 ,NormFactor);
                  fillHisto("hData_CRtighthighpt_2VtxAll_Mass","",samplename, Vtx_Mass1 ,NormFactor);
                  fillHisto("hData_CRtighthighpt_2VtxAll_BDTvtx","",samplename, BDTvtx1,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_BDTvtx","",samplename, BDTvtx2 ,NormFactor);

                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight0 ,NormFactor);
                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","",samplename, Vtx_SumtrackWeight1 ,NormFactor);

                  fillHisto("hData_CRtighthighpt_2VtxAll_MVA","",samplename,  minitree_Hemi_Vtx_MVAval_Tight->at(0),NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_MVA","",samplename,  minitree_Hemi_Vtx_MVAval_Tight->at(1),NormFactor );

                  // !! $$ ----------------------
                  fillHisto("hData_CRtighthighpt_2Vtx_NEWMVA","",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2Vtx_NEWMVA","diffBin",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2Vtx_NEWMVA","6Bins",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2Vtx_NEWMVA","7Bins",samplename,Best_BDTVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2Vtx_NEWMVA","8Bins",samplename,Best_BDTVal,NormFactor );

                  fillHisto("hData_CRtighthighpt_2Vtx_NEWMVA","Sum",samplename,SumBDTVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2Vtx_NEWMVA","Ave",samplename,AveBDTVal,NormFactor );

                  fillHisto("hData_CRtighthighpt_2Vtx_STW","Sum",samplename,SumSTW,NormFactor );
                  fillHisto("hData_CRtighthighpt_2Vtx_STW","Ave",samplename,AveSTW,NormFactor );

                  fillHisto("hData_CRtighthighpt_2Vtx_STW","SumdiffBin",samplename,SumSTW,NormFactor );
                  fillHisto("hData_CRtighthighpt_2Vtx_STW","AvediffBin",samplename,AveSTW,NormFactor );
                  // !! $$ ----------------------


                  // !! ----------------------- VtxBDT variables

                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_NChi2","",samplename,Vtx_NChi0,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_NChi2","",samplename,Vtx_NChi1,NormFactor );

                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_z","",samplename,Vtx_z0,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_z","",samplename,Vtx_z1,NormFactor );

                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_dist","",samplename,Vtx_dist0,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_dist","",samplename,Vtx_dist1,NormFactor );

                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_ntrk10","",samplename,Vtx0_ntrk10,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_ntrk10","",samplename,Vtx1_ntrk10,NormFactor );

                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_MeanDCA","",samplename,Vtx_track_MeanDCA_d0,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_MeanDCA","",samplename,Vtx_track_MeanDCA_d1,NormFactor );

                  // !! --------------------

                              // !! Control Plots
                  fillHisto2D("TrackerMap","BestVtx_C",samplename,Vtx_x,Vtx_y,1);//500,-25.,25.,500,-25.,25.
                  fillHisto2D("TrackerMap","rz_BestVtx_C",samplename,Vtx_z,Vtx_r,1);//500,-25.,25.,500,-25.,25.
                  fillHisto2D("TrackerMap","BothVtx_C",samplename,Vtx_x0,Vtx_y0,1);//500,-25.,25.,500,-25.,25.
                  fillHisto2D("TrackerMap","BothVtx_C",samplename,Vtx_x1,Vtx_y1,1 );//500,-25.,25.,500,-25.,25.
                  fillHisto2D("TrackerMap","rz_BothVtx_C",samplename,Vtx_z0,Vtx_r0,1);//500,-25.,25.,500,-25.,25.
                  fillHisto2D("TrackerMap","rz_BothVtx_C",samplename,Vtx_z1,Vtx_r1,1);//500,-25.,25.,500,-25.,25.

                  fillHisto("Step","BestVtx_C",samplename,Vtx_step,NormFactor );
                  fillHisto("r","BestVtx_C",samplename,Vtx_r,NormFactor );
                  fillHisto("z","BestVtx_C",samplename,abs(Vtx_z),NormFactor );
                  fillHisto("STW","BestVtx_C",samplename,Vtx_SumtrackWeight,NormFactor );
                  fillHisto("Mass","BestVtx_C",samplename,VtxMass,NormFactor );
                  fillHisto("ntrk10","BestVtx_C",samplename,Vtx_ntrk10,NormFactor );
                  fillHisto("MeanDCA","BestVtx_C",samplename,Vtx_track_MeanDCA_d,NormFactor );
                     // !! Control Plots

                  fillHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","",samplename,Vtx2_bdtVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","6Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","6Bins",samplename,Vtx2_bdtVal,NormFactor );

                  fillHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","7Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","7Bins",samplename,Vtx2_bdtVal,NormFactor );

                  fillHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","8Bins",samplename,Vtx1_bdtVal,NormFactor );
                  fillHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","8Bins",samplename,Vtx2_bdtVal,NormFactor );


                  fillHisto("hData_CRtighthighpt_2VtxAll_dist_VtxVtx","",samplename,  Vtx_Vtx_dist,NormFactor );


                  fillHisto("hData_CRtighthighpt_TLVtx_SumtrackWeight","",samplename,Vtx_SumtrackWeight,NormFactor);
                  fillHisto("hData_CRtighthighpt_TLVtx_SumtrackWeight","diffBin",samplename,Vtx_SumtrackWeight,NormFactor);


                  fillHisto("StepEffi","",samplename, 7,1 );// at least 1 tight vtx with at least 1 hemi witth pt > 100 GeV
                  if (Vtx_SumtrackWeight0 > 8 ) { Vtx_SumtrackWeight0 = 8.5;}
                  if (Vtx_SumtrackWeight1 > 8 ) { Vtx_SumtrackWeight1 = 8.5;}
                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRtighthighpt_TLVtx_SumtrackWeight","8Bins",samplename,Vtx_SumtrackWeight,NormFactor);


                                    if (Vtx_SumtrackWeight0 > 7 ) { Vtx_SumtrackWeight0 = 7.5;}
                  if (Vtx_SumtrackWeight1 > 7 ) { Vtx_SumtrackWeight1 = 7.5;}
                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRtighthighpt_TLVtx_SumtrackWeight","7Bins",samplename,Vtx_SumtrackWeight,NormFactor);

                  
                                    if (Vtx_SumtrackWeight0 > 6 ) { Vtx_SumtrackWeight0 = 6.5;}
                  if (Vtx_SumtrackWeight1 > 6 ) { Vtx_SumtrackWeight1 = 6.5;}
                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight0,NormFactor);
                  fillHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight1,NormFactor);
                  fillHisto("hData_CRtighthighpt_TLVtx_SumtrackWeight","6Bins",samplename,Vtx_SumtrackWeight,NormFactor);


               }

         }

         fillHisto("NewBDT_GlobalOutput","BestScore",samplename,Best_BDTVal,NormFactor );
         fillHisto("NewBDT_GlobalOutput","BothScore",samplename,Vtx1_bdtVal,NormFactor );
         fillHisto("NewBDT_GlobalOutput","BothScore",samplename,Vtx2_bdtVal,NormFactor );

         fillHisto("NewEVTBDT_GlobalOutput","Score",samplename,VTX_BDTVal_EVT,NormFactor );



        
   }// end of event loop

   //add end accolade
      NormFactor = NormFactor;
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


   fL1_Reco_SF->Close();
   fL1_ID_SF->Close();
   fL1_ISO_SF->Close();
   fL2_ISO_SF->Close();
   fL2_Reco_SF->Close();
   fL2_ID_SF->Close();
   fL1L2_TRG_SF->Close();
   fL2_Reco_SF2->Close(); 
   fL1_TRG_SF->Close();
   fL1L2_TRG_SFerr->Close();
   fEle_SF->Close();
   fEle_2DSF->Close();

   theoutputfile->Close();
   // gROOT->GetListOfFiles()->Print();
      // // std::cout<<"10"<<std::endl;

   // delete theoutputfile;


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

// !! $$ --

// addHistoDiffBin(TString var, TString selstep, TString sample, int nbins, float* binedges)
const int nbin4 = 4;
const int nbin5 = 5;
const int nbin6 = 6;
const int nbin7 = 7;
float binEdgesBestBDT[nbin6+1] = {-1,-0.7,-0.4,-0.1,0.2,0.5,1};
float binEdgesSTW[nbin5+1] = {1,2,3,4,5,20};
float binEdgesSumSTW[nbin5+1] = {2,3,4,5,6,20};
float binEdgesAveSTW[nbin4+1] = {1,2,3,4,20};
float binEdgesNewBDT[nbin5+1] = {-1,-0.8,-0.6,0.,0.8,1};
float binEdgesNewBDT_v2[nbin6+1] = {-1,-0.8,-0.6,0.,0.8,0.99,1};

// !! -- $$ ££ $$ ¤¤

addHisto("NewEVTBDT_GlobalOutput","Score",sample.Data(),20,-1,1);

// !! -- $$ ££ $$ ¤¤
addHisto("2Vtx_HighSTW_Highpt_BDT","C",sample.Data(),20,-1,1);
addHisto("2Vtx_LowSTW_Highpt_BDT","D",sample.Data(),20,-1,1 );
addHisto("2Vtx_LowSTW_Lowpt_BDT","B",sample.Data(),20,-1,1 );
addHisto("2Vtx_HighSTW_Lowpt_BDT","A",sample.Data(),20,-1,1);



addHistoDiffBin("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","diffBin",sample.Data(),nbin6,binEdgesBestBDT );
addHistoDiffBin("hData_CRlooselooselowlowpt_2Vtx_STW","SumdiffBin",sample.Data(),nbin5,binEdgesSumSTW );
addHistoDiffBin("hData_CRlooselooselowlowpt_2Vtx_STW","AvediffBin",sample.Data(),nbin4,binEdgesAveSTW );
addHistoDiffBin("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","diffBin",sample.Data(),nbin5,binEdgesSTW);

addHistoDiffBin("hData_CRtightlowlowpt_2Vtx_NEWMVA","diffBin",sample.Data(),nbin6,binEdgesBestBDT );
addHistoDiffBin("hData_CRtightlowlowpt_2Vtx_STW","SumdiffBin",sample.Data(),nbin5,binEdgesSumSTW );
addHistoDiffBin("hData_CRtightlowlowpt_2Vtx_STW","AvediffBin",sample.Data(),nbin4,binEdgesAveSTW );
addHistoDiffBin("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","diffBin",sample.Data(),nbin5,binEdgesSTW);

addHistoDiffBin("hData_CRlooselooselowpt_2Vtx_NEWMVA","diffBin",sample.Data(),nbin6,binEdgesBestBDT);
addHistoDiffBin("hData_CRlooselooselowpt_2Vtx_STW","SumdiffBin",sample.Data(),nbin5,binEdgesSumSTW );
addHistoDiffBin("hData_CRlooselooselowpt_2Vtx_STW","AvediffBin",sample.Data(),nbin4,binEdgesAveSTW );
addHistoDiffBin("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","diffBin",sample.Data(),nbin5,binEdgesSTW);

addHistoDiffBin("hData_CRtighthighpt_2Vtx_NEWMVA","diffBin",sample.Data(),nbin6,binEdgesBestBDT );
addHistoDiffBin("hData_CRtighthighpt_2Vtx_STW","SumdiffBin",sample.Data(),nbin5,binEdgesSumSTW );
addHistoDiffBin("hData_CRtighthighpt_2Vtx_STW","AvediffBin",sample.Data(),nbin4,binEdgesAveSTW );
addHistoDiffBin("hData_CRtighthighpt_TLVtx_SumtrackWeight","diffBin",sample.Data(),nbin5,binEdgesSTW);



// !! ££ Control Plots
addHisto("NewBDT_GlobalOutput","BestScore",sample.Data(),20,-1,1);
   addHisto("NewBDT_GlobalOutput","BothScore",sample.Data(),20,-1,1 );


addHisto2D("TrackerMap","BestVtx_C",sample.Data(),500,-25.,25.,500,-25.,25.);//500,-25.,25.,500,-25.,25.
addHisto2D("TrackerMap","BothVtx_C",sample.Data(),500,-25.,25.,500,-25.,25.);

addHisto2D("TrackerMap","BestVtx_B",sample.Data(),500,-25.,25.,500,-25.,25.);//500,-25.,25.,500,-25.,25.
addHisto2D("TrackerMap","BothVtx_B",sample.Data(),500,-25.,25.,500,-25.,25.);

addHisto2D("TrackerMap","BestVtx_D",sample.Data(),500,-25.,25.,500,-25.,25.);//500,-25.,25.,500,-25.,25.
addHisto2D("TrackerMap","BothVtx_D",sample.Data(),500,-25.,25.,500,-25.,25.);

addHisto2D("TrackerMap","BestVtx_A",sample.Data(),500,-25.,25.,500,-25.,25.);//500,-25.,25.,500,-25.,25.
addHisto2D("TrackerMap","BothVtx_A",sample.Data(),500,-25.,25.,500,-25.,25.);


addHisto2D("TrackerMap","rz_BestVtx_C",sample.Data(),1200,0.,120.,700,0.,70.);//1200,0.,120.,700,0.,70.
addHisto2D("TrackerMap","rz_BothVtx_C",sample.Data(),1200,0.,120.,700,0.,70.);

addHisto2D("TrackerMap","rz_BestVtx_B",sample.Data(),1200,0.,120.,700,0.,70.);//1200,0.,120.,700,0.,70.
addHisto2D("TrackerMap","rz_BothVtx_B",sample.Data(),1200,0.,120.,700,0.,70.);

addHisto2D("TrackerMap","rz_BestVtx_D",sample.Data(),1200,0.,120.,700,0.,70.);//1200,0.,120.,700,0.,70.
addHisto2D("TrackerMap","rz_BothVtx_D",sample.Data(),1200,0.,120.,700,0.,70.);

addHisto2D("TrackerMap","rz_BestVtx_A",sample.Data(),1200,0.,120.,700,0.,70.);//1200,0.,120.,700,0.,70.
addHisto2D("TrackerMap","rz_BothVtx_A",sample.Data(),1200,0.,120.,700,0.,70.);



   addHisto("VtxLowSTW_Hemipt","2Vtx",sample.Data(),25,0,500);
   addHisto("VtxHighSTW_Hemipt","2Vtx",sample.Data(),25,0,500);


   addHisto("Step","BestVtx_A",sample.Data(),4,0,4 );
   addHisto("r","BestVtx_A",sample.Data(),50,0,100);
   addHisto("z","BestVtx_A",sample.Data(),50,0,100);
   addHisto("STW","BestVtx_A",sample.Data(),50,0,5 );
   addHisto("Mass","BestVtx_A",sample.Data(),20,0,100 );
   addHisto("ntrk10","BestVtx_A",sample.Data(),20,0,20 );
   addHisto("MeanDCA","BestVtx_A",sample.Data(),100,0,10 );

   addHisto("Step","BestVtx_B",sample.Data(),4,0,4 );
   addHisto("r","BestVtx_B",sample.Data(),50,0,100);
   addHisto("z","BestVtx_B",sample.Data(),50,0,100);
   addHisto("STW","BestVtx_B",sample.Data(),50,0,5 );
   addHisto("Mass","BestVtx_B",sample.Data(),20,0,100);
   addHisto("ntrk10","BestVtx_B",sample.Data(),20,0,20 );
   addHisto("MeanDCA","BestVtx_B",sample.Data(),100,0,10);

   addHisto("Step","BestVtx_C",sample.Data(),4,0,4 );
   addHisto("r","BestVtx_C",sample.Data(),50,0,100);
   addHisto("z","BestVtx_C",sample.Data(),50,0,100 );
   addHisto("STW","BestVtx_C",sample.Data(),50,0,5 );
   addHisto("Mass","BestVtx_C",sample.Data(),20,0,100 );
   addHisto("ntrk10","BestVtx_C",sample.Data(),20,0,20 );
   addHisto("MeanDCA","BestVtx_C",sample.Data(),100,0,10 );

   addHisto("Step","BestVtx_D",sample.Data(),4,0,4 );
   addHisto("r","BestVtx_D",sample.Data(),50,0,100 );
   addHisto("z","BestVtx_D",sample.Data(),50,0,100 );
   addHisto("STW","BestVtx_D",sample.Data(),50,0,5 );
   addHisto("Mass","BestVtx_D",sample.Data(),20,0,100 );
   addHisto("ntrk10","BestVtx_D",sample.Data(),20,0,20 );
   addHisto("MeanDCA","BestVtx_D",sample.Data(),100,0,10 );
// !! ££


      addHisto("hData_Vtx_dist","Tight",sample.Data(), 30,0,150);
      addHisto("hData_Vtx_dist","Loose",sample.Data(), 30,0,150);


      addHisto("hData_VtxQualityTight_VtxBDT","2Vtx",sample.Data(),20,-1,1);
      addHisto("hData_VtxQualityLoose_VtxBDT","2Vtx",sample.Data(),20,-1,1);
      addHisto2D("hData_Hemipt_VtxBDT","2Vtx",sample.Data(),20,-1,1,20,30,230);


// !!
   addHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","",sample.Data(), 20,-1,1 );
   addHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","",sample.Data(), 20,-1,1 );
   addHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","",sample.Data(), 20,-1,1 );
   addHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","",sample.Data(), 20,-1,1 );

   addHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","6Bins",sample.Data(), 6,-1,1 );
   addHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","6Bins",sample.Data(), 6,-1,1 );
   addHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","6Bins",sample.Data(), 6,-1,1 );
   addHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","6Bins",sample.Data(), 6,-1,1 );


   addHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","7Bins",sample.Data(), 7,-1,1 );
   addHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","7Bins",sample.Data(), 7,-1,1 );
   addHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","7Bins",sample.Data(), 7,-1,1 );
   addHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","7Bins",sample.Data(), 7,-1,1 );

      addHisto("hData_CRlooselooselowlowpt_2VtxAll_NEWMVA","8Bins",sample.Data(), 8,-1,1 );
   addHisto("hData_CRtightlowlowpt_2VtxAll_NEWMVA","8Bins",sample.Data(), 8,-1,1 );
   addHisto("hData_CRlooselooselowpt_2VtxAll_NEWMVA","8Bins",sample.Data(), 8,-1,1 );
   addHisto("hData_CRtighthighpt_2VtxAll_NEWMVA","8Bins",sample.Data(), 8,-1,1 );


   // !! $$    -----------------
addHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","",sample.Data(),20,-1,1);
addHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","6Bins",sample.Data(),6,-1,1);
addHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","7Bins",sample.Data(),7,-1,1);
addHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","8Bins",sample.Data(),8,-1,1 );
addHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
addHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","8Bins",sample.Data(),8,1,9);
addHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","7Bins",sample.Data(),7,1,8);
addHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","6Bins",sample.Data(),6,1,7);

addHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","",sample.Data(),20,-1,1 );
addHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","6Bins",sample.Data(),6,-1,1);
addHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","7Bins",sample.Data(),7,-1,1 );
addHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","8Bins",sample.Data(),8,-1,1);
// addHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
addHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","8Bins",sample.Data(),8,1,9);
addHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","7Bins",sample.Data(),7,1,8);
addHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","6Bins",sample.Data(),6,1,7);
            
addHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","",sample.Data(),20,-1,1 );
addHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","6Bins",sample.Data(),6,-1,1 );
addHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","7Bins",sample.Data(),7,-1,1);
addHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","8Bins",sample.Data(),8,-1,1 );
addHisto("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
addHisto("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","8Bins",sample.Data(),8,1,9);
addHisto("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","7Bins",sample.Data(),7,1,8);
addHisto("hData_CRlooselooselowpt_TLVtx_SumtrackWeight","6Bins",sample.Data(),6,1,7);


addHisto("hData_CRtighthighpt_2Vtx_NEWMVA","",sample.Data(),20,-1,1);
addHisto("hData_CRtighthighpt_2Vtx_NEWMVA","6Bins",sample.Data(),6,-1,1 );
addHisto("hData_CRtighthighpt_2Vtx_NEWMVA","7Bins",sample.Data(),7,-1,1 );
addHisto("hData_CRtighthighpt_2Vtx_NEWMVA","8Bins",sample.Data(),8,-1,1 );
addHisto("hData_CRtighthighpt_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
addHisto("hData_CRtighthighpt_TLVtx_SumtrackWeight","8Bins",sample.Data(),8,1,9);
addHisto("hData_CRtighthighpt_TLVtx_SumtrackWeight","7Bins",sample.Data(),7,1,8);
addHisto("hData_CRtighthighpt_TLVtx_SumtrackWeight","6Bins",sample.Data(),6,1,7);
   // !! $$    -----------------


// !! $$$ --------------------

addHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","Sum",sample.Data(),20,-2,2);
addHisto("hData_CRlooselooselowlowpt_2Vtx_NEWMVA","Ave",sample.Data(),20,-1,1);

addHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","Sum",sample.Data(),20,-2,2);
addHisto("hData_CRtightlowlowpt_2Vtx_NEWMVA","Ave",sample.Data(),20,-1,1 );

addHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","Sum",sample.Data(),20,-2,2 );
addHisto("hData_CRlooselooselowpt_2Vtx_NEWMVA","Ave",sample.Data(),20,-1,1 );

addHisto("hData_CRtighthighpt_2Vtx_NEWMVA","Sum",sample.Data(),20,-2,2 );
addHisto("hData_CRtighthighpt_2Vtx_NEWMVA","Ave",sample.Data(),20,-1,1);


addHisto("hData_CRlooselooselowlowpt_2Vtx_STW","Sum",sample.Data(),19,1,20);
addHisto("hData_CRlooselooselowlowpt_2Vtx_STW","Ave",sample.Data(),19,1,20);

addHisto("hData_CRtightlowlowpt_2Vtx_STW","Sum",sample.Data(),19,1,20);
addHisto("hData_CRtightlowlowpt_2Vtx_STW","Ave",sample.Data(),19,1,20 );

addHisto("hData_CRlooselooselowpt_2Vtx_STW","Sum",sample.Data(),19,1,20 );
addHisto("hData_CRlooselooselowpt_2Vtx_STW","Ave",sample.Data(),19,1,20 );

addHisto("hData_CRtighthighpt_2Vtx_STW","Sum",sample.Data(),19,1,20 );
addHisto("hData_CRtighthighpt_2Vtx_STW","Ave",sample.Data(),19,1,20);
// !! $$$ --------------------
// !! 

   addHisto("hData_CRlooselooselowlowpt_2VtxAll_MVA","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlooselooselowlowpt_dist_VtxVtx","",sample.Data(),50,0,100);
   addHisto("hData_CRtightlowlowpt_2VtxAll_MVA","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRtightlowlowpt_dist_VtxVtx","",sample.Data(),50,0,100 );
   addHisto("hData_CRlooselooselowpt_2VtxAll_MVA","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlooselooselowpt_dist_VtxVtx","",sample.Data(),50,0,100 );
   addHisto("hData_CRtighthighpt_2VtxAll_MVA","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRtighthighpt_2VtxAll_dist_VtxVtx","",sample.Data(),50,0,100 );



                  // !! ----------------------- VtxBDTVariables

                  addHisto("hData_CRlooselooselowlowpt_2VtxAll_NChi2","",sample.Data(),100,0,10 );

                  addHisto("hData_CRlooselooselowlowpt_2VtxAll_z","",sample.Data(),20,-100,100 );

                  addHisto("hData_CRlooselooselowlowpt_2VtxAll_dist","",sample.Data(),30,0,150 );

                  addHisto("hData_CRlooselooselowlowpt_2VtxAll_ntrk10","",sample.Data(),30,0,30 );

                  addHisto("hData_CRlooselooselowlowpt_2VtxAll_MeanDCA","",sample.Data(),60,0,15 );


                  addHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_NChi2","",sample.Data(),100,0,10 );

                  addHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_z","",sample.Data(),20,-100,100 );

                  addHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_dist","",sample.Data(),30,0,150);

                  addHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_ntrk10","",sample.Data(),30,0,30);

                  addHisto("hData_CRtightlowlowpt_2VtxAll_2VtxAll_MeanDCA","",sample.Data(),60,0,15);


                  addHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_NChi2","",sample.Data(),100,0,10 );

                  addHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_z","",sample.Data(),20,-100,100  );

                  addHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_dist","",sample.Data(),30,0,150 );

                  addHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_ntrk10","",sample.Data(),30,0,30 );

                  addHisto("hData_CRlooselooselowpt_2VtxAll_2VtxAll_MeanDCA","",sample.Data(),60,0,15 );


                  addHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_NChi2","",sample.Data(),100,0,10 );

                  addHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_z","",sample.Data(),20,-100,100  );

                  addHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_dist","",sample.Data(),30,0,150);

                  addHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_ntrk10","",sample.Data(),30,0,30 );

                  addHisto("hData_CRtighthighpt_2VtxAll_2VtxAll_MeanDCA","",sample.Data(),60,0,15 );



addHisto("hSim_Hemi_Vtx_r","noSel",sample.Data(),50,0,100);

addHisto("hSim_Hemi_Vtx_dist","noSel",sample.Data(),50,0,100);

addHisto("hData_Hemi_Vtx_r","Goodrecovtx",sample.Data(),50,0,100);
addHisto("hData_Hemi_Vtx_dist","Goodrecovtx",sample.Data(),50,0,100);

addHisto("hData_Hemi_Vtx_r","Ping",sample.Data(),50,0,100);
addHisto("hData_Hemi_Vtx_dist","Ping",sample.Data(),50,0,100);

addHisto("hData_Hemi_Vtx_r","TightPing",sample.Data(),50,0,100);
addHisto("hData_Hemi_Vtx_dist","TightPing",sample.Data(),50,0,100);

addHisto("hData_Hemi_Vtx_r","LoosePing",sample.Data(),50,0,100);
addHisto("hData_Hemi_Vtx_dist","LoosePing",sample.Data(),50,0,100);



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

   addHisto("hData_VtxQualityTight_Hemipt","2Vtx",sample.Data(), 50,0,1000);
   addHisto("hData_VtxQualityLoose_Hemipt","2Vtx",sample.Data(), 50,0,1000);


   addHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","NoSel",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_Hemi_2VtxAll_NChi2","",sample.Data(),64,-1,15);


   addHisto("Hemisphere_leadingpt","", sample.Data(), 50,0,1000);
   addHisto("Hemisphere_subleadingpt","", sample.Data(), 50,0,1000);


   addHisto("hData_Hemi_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_Hemi_TLVtx_Mass","",sample.Data(),25,0.,100.);

   addHisto("hData_Hemi_TLVtxAll_SumtrackWeight","",sample.Data(),19,1,20);

   addHisto("hData_Hemi_TLVtxAll_Mass","",sample.Data(),25,0.,100.);

   addHisto2D("hData_Hemi_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_Hemi_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_Hemi_2Vtx_Mmumu","",sample.Data(),25,0.,500.);
   addHisto("hData_Hemi_2Vtx_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_Hemi_2Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_Hemi_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_Hemi_2Vtx_SumtrackWeight","",sample.Data(), 19,1,20 );
   addHisto2D("hData_Hemi_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_Hemi_2VtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_Hemi_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1);
   addHisto("hData_Hemi_2VtxAll_SumtrackWeight","",sample.Data(),19,1,20 );
   addHisto2D("hData_Hemi_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);

   addHisto("hData_CRlowpt_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_CRlowpt_TLVtxAll_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRlowpt_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlowpt_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlowpt_2Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlowpt_2Vtx_Mass","",sample.Data(),   25,0.,100.);
   addHisto("hData_CRlowpt_2Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowpt_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1);
   addHisto("hData_CRlowpt_2Vtx_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto2D("hData_CRlowpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlowpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRlowpt_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowpt_2VtxAll_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto2D("hData_CRlowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);

   addHisto("hData_CRlowlowpt_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRlowlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRlowlowpt_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlowlowpt_TLVtxAll_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRlowlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRlowlowpt_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlowlowpt_2Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlowlowpt_2Vtx_Mass","",sample.Data(),   25,0.,100.);
   addHisto("hData_CRlowlowpt_2Vtx_BDTevt","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlowlowpt_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1);
   addHisto("hData_CRlowlowpt_2Vtx_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto("hData_CRlowlowpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRlowlowpt_2VtxAll_BDTvtx","",sample.Data(),25,-1,1);
   addHisto("hData_CRlowlowpt_2VtxAll_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto2D("hData_CRlowlowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);

   addHisto("hData_CRloose_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRloose_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRloose_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRloose_TLVtxAll_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRloose_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRloose_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRloose_2Vtx_Mmumu","",sample.Data(),  25,0.,500.);
   addHisto("hData_CRloose_2Vtx_Mass","",sample.Data(),   25,0.,100. );
   addHisto("hData_CRloose_2Vtx_BDTevt","",sample.Data(), 25,-1,1);
   addHisto("hData_CRloose_2Vtx_MaxBDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRloose_2Vtx_SumtrackWeight","",sample.Data(), 19,1,20 );
   addHisto2D("hData_CRloose_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRloose_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRloose_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRloose_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1 );           
   addHisto("hData_CRloose_2VtxAll_SumtrackWeight","",sample.Data(), 19,1,20 );


   addHisto("hData_CRlooselowpt_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRlooselowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_CRlooselowpt_TLVtxAll_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRlooselowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRlooselowpt_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselowpt_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlooselowpt_2Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRlooselowpt_2Vtx_Mass","",sample.Data(),   25,0.,100. );
   addHisto("hData_CRlooselowpt_2Vtx_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto2D("hData_CRlooselowpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlooselowpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRlooselowpt_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlooselowpt_2VtxAll_SumtrackWeight","",sample.Data(), 19,1,20);


   addHisto("hData_CRlooselowlowpt_TLVtx_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRlooselowlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRlooselowlowpt_TLVtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselowlowpt_TLVtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlooselowlowpt_TLVtxAll_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto("hData_CRlooselowlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_CRlooselowlowpt_2Vtx_Mmumu","",sample.Data(), 25,0.,500. );
   addHisto("hData_CRlooselowlowpt_2Vtx_Mass","",sample.Data(),   25,0.,100. );
   addHisto("hData_CRlooselowlowpt_2Vtx_SumtrackWeight","",sample.Data(),19,1,20);
   addHisto2D("hData_CRlooselowlowpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselowlowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlooselowlowpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRlooselowlowpt_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1 );
   addHisto("hData_CRlooselowlowpt_2VtxAll_SumtrackWeight","",sample.Data(), 19,1,20);


   addHisto("hData_StepEff","",sample.Data(),9, 0, 9);
   addHisto("hData_StepEff","NonNormalized",sample.Data(),9, 0, 9);

   addHisto("StepEffi","",sample.Data(),8, 0, 8);

    // --------------- zzzzzzzzzzzzzzzzzzzzzzzzzzzz ------------------------------------//

//A
   addHisto("hData_CRtightlowlowpt_BDTevt","",sample.Data(),  25,-1,1  );
   addHisto("hData_CRtightlowlowpt_TLVtx_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto("hData_CRtightlowlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","6Bins",sample.Data(), 6,1,7);
   addHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","7Bins",sample.Data(), 7,1,8);
   addHisto("hData_CRtightlowlowpt_TLVtxAll_SumtrackWeight","8Bins",sample.Data(), 8,1,9);
   addHisto("hData_CRtightlowlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRtightlowlowpt_TLVtx_STW_Ntrks","",sample.Data(), 25,0,25,25,0,25);
   addHisto2D("hData_CRtightlowlowpt_TLVtxAll_STW_Ntrks","",sample.Data() ,25,0,25,25,0,25);

   //B
   addHisto("hData_CRlooselooselowlowpt_BDTevt","",sample.Data(),  25,-1,1  );
   // addHisto("hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto("hData_CRlooselooselowlowpt_TLVtx_Mass","",sample.Data(),25,0.,100.);
   addHisto2D("hData_CRlooselooselowlowpt_TLVtx_STW_Ntrks","",sample.Data() ,25,0,25,25,0,25);
   addHisto2D("hData_CRlooselooselowlowpt_TLVtxAll_STW_Ntrks","",sample.Data() ,25,0,25,25,0,25);
   addHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","6Bins",sample.Data(),6,1,7);
   addHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","7Bins",sample.Data(),7,1,8);
   addHisto("hData_CRlooselooselowlowpt_TLVtxAll_SumtrackWeight","8Bins",sample.Data(),8,1,9);
   addHisto("hData_CRlooselooselowlowpt_TLVtxAll_Mass","",sample.Data(),25,0.,100.);


//C
   addHisto("hData_CRtighthighpt_BDTevt","",sample.Data(),25,-1,1  );
   addHisto("hData_CRtighthighpt_2Vtx_Mmumu","",sample.Data(),  25,0.,500. );
   addHisto("hData_CRtighthighpt_2Vtx_Mass","",sample.Data(),   25,0.,100.);
   addHisto("hData_CRtighthighpt_2Vtx_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto2D("hData_CRtighthighpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRtighthighpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRtighthighpt_2VtxAll_Mass","",sample.Data(), 25,0.,100.);
   addHisto("hData_CRtighthighpt_2VtxAll_BDTvtx","",sample.Data(),  25,-1,1 );
   addHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","6Bins",sample.Data(), 6,1,7);
   addHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","7Bins",sample.Data(), 7,1,8);
   addHisto("hData_CRtighthighpt_2VtxAll_SumtrackWeight","8Bins",sample.Data(), 8,1,9);
   
//D
   addHisto("hData_CRlooselooselowpt_BDTevt","",sample.Data(), 25,-1,1  );
   addHisto("hData_CRlooselooselowpt_2Vtx_Mmumu","",sample.Data(),   25,0.,500.);
   addHisto("hData_CRlooselooselowpt_2Vtx_Mass","",sample.Data(),  25,0.,100. );
   addHisto("hData_CRlooselooselowpt_2Vtx_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto2D("hData_CRlooselooselowpt_2Vtx_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto2D("hData_CRlooselooselowpt_2VtxAll_STW_Ntrks","",sample.Data(),25,0,25,25,0,25);
   addHisto("hData_CRlooselooselowpt_2VtxAll_Mass","",sample.Data(),25,0.,100.);
   addHisto("hData_CRlooselooselowpt_2VtxAll_BDTvtx","",sample.Data(), 25,-1,1  );
   addHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","",sample.Data(), 19,1,20);
   addHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","6Bins",sample.Data(), 6,1,7);
   addHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","7Bins",sample.Data(), 7,1,8);
   addHisto("hData_CRlooselooselowpt_2VtxAll_SumtrackWeight","8Bins",sample.Data(), 8,1,9);


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

void TreeABCDReader::addHistoDiffBin(TString var, TString selstep, TString sample, int nbins, float* binedges){
 
   TString name =  sample+"_"+var+"_"+selstep;
  TH1F * thehisto = new TH1F(name,name,nbins,binedges);
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
