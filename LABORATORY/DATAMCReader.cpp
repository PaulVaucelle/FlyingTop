#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <TH1F.h>
#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/MCWeights.h"
#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/DATAMCREADER/EMU/PDFSCALEVar.h"
#include <TPad.h>
#include <TCanvas.h>


float get2DSFFromHisto(float pt, float eta,  TH2F* hist) {
    int bin_x = hist->GetXaxis()->FindFixBin(eta);
    int bin_y = hist->GetYaxis()->FindFixBin(pt);

    // Protection contre les dépassements de bin
    bin_x = std::max(1, std::min(bin_x, hist->GetNbinsX()));
    bin_y = std::max(1, std::min(bin_y, hist->GetNbinsY()));

    return hist->GetBinContent(bin_x, bin_y);
}
// FindBin


//  int  Year, bool isPostAPV,  int Channel ,bool DoubleMuon, TString thesystlist

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> <output_file.root>" << std::endl;
        return 1;
    }

    TString input_file = argv[1];
    TString output_file = argv[2];

    TString Prod = argv[3];
    int Year = atoi(argv[4]);
    bool isPostAPV = atoi(argv[5]);
   int Channel = atoi(argv[6]);
   bool DoubleMuon = atoi(argv[7]);
   // TString thesystlist = argv[8];

    float MeanGenW = 1;
    bool isMC = false;
    bool Signal = true;
   // !! Sample SF -- //
    float XS = 1.;
    float FracEffEvent = 1.;
    if (input_file.Contains("DYJetsToLL_M-10to50"))                     { XS = 22635; FracEffEvent = 1.;  MeanGenW = 1; isMC = true; Signal = false;}
    if (input_file.Contains("DYJetsToLL_M-50"))                         { XS = 6225.4;  FracEffEvent = 1.;  MeanGenW = 1; isMC = true;Signal = false;}
    if (input_file.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays")) { XS = 21.6;  FracEffEvent = 1.; MeanGenW = 32.5092; isMC = true;Signal = false;}
    if (input_file.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     { XS = 21.6;  FracEffEvent = 1.; MeanGenW = 32.4473; isMC = true;Signal = false;}
    if (input_file.Contains("TTJets_DiLept"))                          { XS = 53.07;  isMC = true;  Signal = false; }
    if (input_file.Contains("TTJets"))                          { XS = 831.76;  FracEffEvent = 0.362;  isMC = true; }
    if (input_file.Contains("TTTo2L2Nu"))                              { XS = 88.5;  FracEffEvent = 0.993;  MeanGenW = 72.6983 ; isMC = true;Signal = false;}
    if (input_file.Contains("TTToSemiLeptonic") )                      { XS = 366.3;  FracEffEvent = 0.993; MeanGenW = 303.358; isMC = true;Signal = false;}
    if (input_file.Contains("TTToHadronic"))                           { XS = 378.9;   isMC = true; Signal = false; }// TOP 20 -006 :377.6
    if (input_file.Contains("WWTo2L2Nu"))                              { XS = 11.09;  FracEffEvent = 0.997; MeanGenW = 11.1262;  isMC = true;Signal = false;}
    if (input_file.Contains("WZTo2Q2L_mllmin4p0"))                     { XS = 6.535;  FracEffEvent = 0.583; MeanGenW = 15.3793; isMC = true;Signal = false;}
    if (input_file.Contains("ZZTo2Q2L_mllmin4p0"))                     { XS = 3.676;  FracEffEvent = 0.612; MeanGenW =  8.51008 ;isMC = true;Signal = false; }
    if (input_file.Contains("ttWJetsToLNu_5f_EWK"))                    { XS = 0.290;  FracEffEvent = 0.551; MeanGenW =  0.0170462; isMC = true;Signal = false; } // not found on XSDB, no file on tier2...approximation
    if (input_file.Contains("TTZToLL_5f"))                             { XS = 0.253 ;FracEffEvent = 1.;   MeanGenW = 1.; isMC = true;Signal = false;}// !!  0.253 Not found on XSDB => used ana.py macro : 0.05188  +- 2.437e-04 pb
    if (input_file.Contains("TTWW"))                                   { XS = 0.006992; FracEffEvent = 1.;  MeanGenW = 1.;isMC = true;Signal = false;}//found on XSDB
    if (input_file.Contains("ST_t-channel_antitop_5f_InclusiveDecays")) {XS = 80.95; FracEffEvent = 0.995;  MeanGenW = 72.1425;isMC = true;Signal = false;}
    if (input_file.Contains("ST_t-channel_top_5f_InclusiveDecays")) {XS = 136.02;FracEffEvent = 0.995;    MeanGenW = 120.394;isMC = true;Signal = false;}
   if (input_file.Contains("DoubleMuon") || input_file.Contains("MuonEG")) {Signal = false; isMC = false; XS = 1.; FracEffEvent = 1.; MeanGenW = 1.;}

    float Lumi = 59800.;
    float LumiUp = 66000.0;
    float LumiDown = 55000.0;

   float ScaleWeight = 1;
   float PdfWeight = 1;

   std::vector<float>   PDFWeights;
   std::vector<float>   ScaleWeights;

   float PUweight = 1.;
   float PUweightUp = 1.;
   float PUweightDown = 1.;
   
   float Prefweight = 1.;
   float PrefweightUp = 1.;
   float PrefweightDown = 1.;

    float NormFactorLumi = 1;   
    float NormFactorLumiUp = 1.;  
    float NormFactorLumiDown = 1.;

    LumiWeights LUMI(Year);
    Lumi = LUMI.GetLumi();
    LumiUp = LUMI.GetLumiUp();
    LumiDown = LUMI.GetLumiDown();

    NormFactorLumi =  XS*Lumi;  
    NormFactorLumiUp = XS*LumiUp;
    NormFactorLumiDown = XS*LumiDown;

          //---------------------------------------------------------------------------//
      // !! ----------------------------- Scale Factors---------------------------------//
      //---------------------------------------------------------------------------//

   //      // Where L1 is always a muon and L2 is either a muon for the dimuon channel or a
   // // an electron for the Emu channel

   TFile*            fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
   TFile*            fL1_ID_SF= new TFile("../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root");
   TFile*           fL1_ISO_SF= new TFile("../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root");
   TFile*           fL2_ISO_SF= new TFile("../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root");

   TFile*           fL2_Reco_SF= new TFile("../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root");
   TFile*            fL2_Reco_SF2= new TFile("../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2018.root");
   TFile*            fL2_ID_SF= new TFile("../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_EGM2D_U18.root");
   TFile*            fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2018_ULv2.root");
   TFile*            fL1_TRG_SF = new TFile("../Scale_factors/Muon/DoubleMuon_2018/NUM_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_DEN_MiniIsoTight_and_TightID_abseta_pt.root");
   
   TFile*            fL1L2_TRG_SFerr = new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/Scale_factors/Trig_2018/trigSF_2Derr.root");

 // !! ----------------------------
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
               fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../Scale_factors/Muon/Run2018_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL1_ISO_SF= new TFile("../Scale_factors/Muon/Run2018_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

               fL2_Reco_SF= new TFile("../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2018.root");
               fL2_ISO_SF= new TFile("../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2018.root");
               fL2_ID_SF= new TFile("../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_EGM2D_U18.root");

               fL1L2_TRG_SF= new TFile("../Scale_factors/Trig_2018/trigSF_2D.root");
               fL1L2_TRG_SFerr = new TFile("../Scale_factors/Trig_2018/trigSF_2Derr.root");

               // cout<< " 2018 year  is here"<<endl;

            }
         if(Year == 2017)
            {
               //cout<<" year 2017"<<endl;
               fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../Scale_factors/Muon/Run2017_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL1_ISO_SF= new TFile("../Scale_factors/Muon/Run2017_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
               
               fL2_Reco_SF= new TFile("../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2017.root");
               fL2_ISO_SF= new TFile("../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2017.root");
               fL2_ID_SF= new TFile("../Scale_factors/Electron/egammaEffi.txt_EGM2D_Tight_UL17.root");

               // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2018_ULv2.root");
               fL1L2_TRG_SF= new TFile("../Scale_factors/Trig_2018/trigSF_2D.root");
               fL1L2_TRG_SFerr = new TFile("../Scale_factors/Trig_2017/trigSF_2Derr.root");
            }
         if(Year == 2016)
            {
               if (isPostAPV){
                  //cout<<" post 2016="<<endl;
                  // fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root");
                  // fL1_ID_SF= new TFile("../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root");
                  // fL1_ISO_SF= new TFile("../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root");

                  fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL1_ISO_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

                  fL2_Reco_SF= new TFile("../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root");
                  fL2_ISO_SF= new TFile("../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2016postVFP.root");
                  fL2_ID_SF= new TFile("../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_postVFP_EGM2D.root");

                  // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2016postVFP_ULv2.root");
                  fL1L2_TRG_SF= new TFile("../Scale_factors/Trig_2016/trigSF_2D.root");
                  fL1L2_TRG_SFerr = new TFile("../Scale_factors/Trig_2016/trigSF_2Derr.root");

               }
               else{
                  // cout<<" pre  2016="<<endl;
                  // fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root");
                  // fL1_ID_SF= new TFile("../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root");
                  // fL1_ISO_SF= new TFile("../Scale_factors/Muon/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root");

                  fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL1_ISO_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
                  
                  fL2_Reco_SF= new TFile("../Scale_factors/Electron/egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root");
                  fL2_ISO_SF= new TFile("../Scale_factors/Electron/egammaEffi_ptBelow20.txt_EGM2D_UL2016preVFP.root");
                  fL2_ID_SF= new TFile("../Scale_factors/Electron/egammaEffi.txt_Ele_Tight_preVFP_EGM2D.root");


                  // fL1L2_TRG_SF= new TFile("../Scale_factors/Electron/Top_trigger_group/TriggerSF_2016preVFP_ULv2.root");
                  fL1L2_TRG_SF= new TFile("../Scale_factors/Trig_2016_pre/trigSF_2D.root");
                  fL1L2_TRG_SFerr = new TFile("../Scale_factors/Trig_2016_pre/trigSF_2Derr.root");
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
               fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../Scale_factors/Muon/Run2018_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL1_ISO_SF= new TFile("../Scale_factors/Muon/Run2018_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
               
               fL2_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL2_ID_SF= new TFile("../Scale_factors/Muon/Run2018_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL2_ISO_SF= new TFile("../Scale_factors/Muon/Run2018_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

               fL1L2_TRG_SF = new TFile("../Scale_factors/Muon/Run2018_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
               fL1L2_TRG_SFerr = new TFile("../Scale_factors/Muon/Run2018_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");

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
               fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL1_ID_SF= new TFile("../Scale_factors/Muon/Run2017_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL1_ISO_SF= new TFile("../Scale_factors/Muon/Run2017_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
               
               fL2_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
               fL2_ID_SF= new TFile("../Scale_factors/Muon/Run2017_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
               fL2_ISO_SF= new TFile("../Scale_factors/Muon/Run2017_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

               fL1L2_TRG_SF = new TFile("../Scale_factors/Muon/Run2017_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
               fL1L2_TRG_SFerr = new TFile("../Scale_factors/Muon/Run2017_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
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
                  fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL1_ISO_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
                  
                  fL2_Reco_SF =  new TFile("../Scale_factors/Muon/Run2016_UL/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
                  fL2_ID_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL2_ISO_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

                  fL1L2_TRG_SF = new TFile("../Scale_factors/Muon/Run2016_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
                  fL1L2_TRG_SFerr = new TFile("../Scale_factors/Muon/Run2017_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");



               // histos associated to files

               }
               else{
                  // cout<<" pre  2016="<<endl;
                  fL1_Reco_SF =  new TFile("../Scale_factors/Muon/Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root");
                  fL1_ID_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL1_ISO_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");
                  
                  fL2_Reco_SF =  new TFile("../Scale_factors/Muon/Run2016_UL/Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root");
                  fL2_ID_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_TightID_DEN_TrackerMuons/NUM_TightID_DEN_TrackerMuons_abseta_pt.root");
                  fL2_ISO_SF= new TFile("../Scale_factors/Muon/Run2016_UL/NUM_MiniIsoTight_DEN_TightID/NUM_MiniIsoTight_DEN_TightID_abseta_pt.root");

                  fL1L2_TRG_SF = new TFile("../Scale_factors/Muon/Run2016_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
                  fL1L2_TRG_SFerr = new TFile("../Scale_factors/Muon/Run2017_UL/NUM_Trigger_DEN_MiniIsoTight/NUM_Trigger_DEN_MiniIsoTight_abseta_pt.root");
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



    fh2DL1_ID_SF1->SetDirectory(0);
    fh2DL2_ID_SF1->SetDirectory(0);
    fh2DL1_ID_SF1err->SetDirectory(0);
    fh2DL2_ID_SF1err->SetDirectory(0);
    fh2DL1_ISO_SF1->SetDirectory(0);
    fh2DL2_ISO_SF1->SetDirectory(0);
    fh2DL1_ISO_SF1err->SetDirectory(0);
    fh2DL2_ISO_SF1err->SetDirectory(0);
    fh2DL2_Reco_SF1err->SetDirectory(0);
    fh2DL1_Reco_SF1->SetDirectory(0);
    fh2DL2_Reco_SF1->SetDirectory(0);
    fhL1L2_TRG_SF->SetDirectory(0);
    fhL1L2_TRG_SFerr->SetDirectory(0); // Détache l'histogramme du fichier

    fL1_Reco_SF->Close();
    fL1_ID_SF->Close();
    fL1_ISO_SF->Close();
    fL2_ISO_SF->Close();
    fL2_Reco_SF->Close();
    fL2_Reco_SF2->Close();
    fL2_ID_SF->Close();
    fL1L2_TRG_SF->Close();
    fL1_TRG_SF->Close();
    fL1L2_TRG_SFerr->Close();


   // fhL1L2_TRG_SF = (TH2F*)(fL1L2_TRG_SF->Get("h2D_SF_emu_lepABpt_FullError")); // x- axiss pt and y axis pt  
   double norm   = 0.;
   if (!Signal && isMC)
      {
         TString hNorma = "hEvents_with_gen_wt";
         TFile* f1_DY= new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/"+input_file+".root");
         f1_DY->cd("");
         
         TDirectory* dir = f1_DY->GetDirectory("FlyingTop");
         if (dir) {
            dir->cd();  // Change to this directory if needed
            // Now you can access histograms, trees, etc. from the directory
               TH1D*  e1_DY = (TH1D*)gROOT->FindObject(hNorma);
            e1_DY->Sumw2();
            if  ( e1_DY->GetEntries() > 0 ) norm =  e1_DY->GetEntries();
         std::cout<<"norm : "<<norm<<std::endl;
         e1_DY->SetDirectory(0);

         }
         f1_DY->Close();
      }

    // !! Sample SF -- //

    // ROOT::EnableVerboseLogging(true);
    ROOT::EnableImplicitMT(8); // Multi-threading
    ROOT::RDataFrame df("ttree", "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniDATAMC_"+input_file+".root");
    TFile output("./EMU/DATAMC_"+output_file+".root", "RECREATE");

    float NormFactorSYST = NormFactorLumi/norm;
  
   // Appliquer le SF à chaque event via pt et eta
      auto df1 = df.Define("SF_PREF", [](const std::vector<double>& PREF) {
         if (PREF.empty() ) return 1.0; // Protection contre les events sans prefiring
         return PREF[0]; // Prend le premier muon
      }, {"miniPrefweight"});

   // Appliquer le SF à chaque event via pt et eta
      auto df2 = df1.Define("SF_PU", [](const std::vector<double>& PU) {
         if (PU.empty() ) return 1.0; // Protection contre les events sans PU
         return PU[0]; // Prend le premier muon
      }, {"miniPUweight"});

       // Appliquer le SF à chaque event via pt et eta
      auto df3 = df2.Define("SF_TopPT", [](const std::vector<double>& toppt) {
         if (toppt.empty() ) return 1.0; // Protection contre les events sans PU
         return toppt[0]; // Prend le premier muon
      }, {"minitree_genTop_Weight"});


             // Appliquer le SF à chaque event via pt et eta
      auto dfA = df3.Define("SF_GenWeight", [](const std::vector<double>& genweight) {
         if (genweight.empty() ) return 1.0; // Protection contre les events sans PU
         return genweight[0]; // Prend le premier muon
      }, {"minitree_only_gen_wt"});



      auto dfB = dfA.Define("SF_Norm", [NormFactorSYST]() {
         return NormFactorSYST;
      });

      auto dfC = dfB.Define("SF_MGW", [MeanGenW]() {
         return 1./MeanGenW;
      });

      auto dfD = dfC.Define("SF_Frac", [FracEffEvent]() {
         return FracEffEvent;
      });


   ROOT::RDF::RNode dfTopPTDown = dfD; // Valeur initiale au cas où aucune condition ne s'applique
   if (input_file.Contains("TTTo2L2Nu") || thesample.Contains("TTToHadronic") ||  thesample.Contains("TTToSemiLeptonic") )
      {
         auto dfTopPT = dfD.Define("SF_TopPT_Up", [](const std::vector<double>& toppt) {
            if (toppt.empty() ) return 1.0; // Protection contre les events sans PU
            return toppt[0]; 
         }, {"minitree_genTop_Weight"});
         dfTopPTDown = dfD.Define("SF_TopPT_Down", []() {
               return 1.;
            });
      }
   else
      {
         auto dfTopPT = dfD.Define("SF_TopPT_Up", []() {
               return 1.;
            });
         dfTopPTDown = dfD.Define("SF_TopPT_Down", []() {
               return 1.;
            });
      }

   // Appliquer le SF à chaque event via pt et eta
      auto dfE = dfTopPTDown.Define("SF_MUID", [fh2DL1_ID_SF1](const std::vector<float>& pt, const std::vector<float>& eta) {
         if (pt.empty() || eta.empty()) return 1.0f; // Protection contre les events sans muon
         return get2DSFFromHisto(pt[0], eta[0], fh2DL1_ID_SF1); // Prend le premier muon
      }, {"minitree_lepton_leadingpt", "minitree_lepton_leadingeta"});

      auto dfF = dfE.Define("SF_MUISO", [fh2DL1_ISO_SF1](const std::vector<float>& pt, const std::vector<float>& eta) {
         if (pt.empty() || eta.empty()) return 1.0f; // Protection contre les events sans muon
         return get2DSFFromHisto(pt[0], eta[0], fh2DL1_ISO_SF1); // Prend le premier muon
      }, {"minitree_lepton_leadingpt", "minitree_lepton_leadingeta"});

      auto dfG = dfF.Define("SF_MURECO", [fh2DL1_Reco_SF1](const std::vector<float>& pt, const std::vector<float>& eta) {
         if (pt.empty() || eta.empty()) return 1.0f; // Protection contre les events sans muon
         return get2DSFFromHisto(pt[0], eta[0], fh2DL1_Reco_SF1); // Prend le premier muon
      }, {"minitree_lepton_leadingpt", "minitree_lepton_leadingeta"});



   //  NormFactor =  (NormFactorSYST*GenWeight*Prefweight*PUweight*top_pt_wt*Mu_SF*Ele_SF*TriggerSyst) / FracEffEvent;

ROOT::RDF::RNode df9f = dfG; // Valeur initiale au cas où aucune condition ne s'applique

if (MuMu >= 1) {

   auto df6 = dfG.Define("SF2_MUID", [fh2DL2_ID_SF1](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      return get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1);
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});
   auto df7 = df6.Define("SF2_MUISO", [fh2DL2_ISO_SF1](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      return get2DSFFromHisto(pt[0], eta[0], fh2DL2_ISO_SF1); // <-- t'avais une petite coquille ici, tu appelais ID au lieu de ISO
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   auto df8 = df7.Define("SF2_MURECO", [fh2DL2_Reco_SF1](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      return get2DSFFromHisto(pt[0], eta[0], fh2DL2_Reco_SF1);
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   df9 = df8.Define("SF_TRIGGER", [fhL1L2_TRG_SF](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      return get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SF);
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   // !! Syst -- !!//


   auto df9a = df9.Define("SF2_MUID_Up", [fh2DL2_ID_SF1, fh2DL2_ID_SF1err](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1err);
      return f1+f2;
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   auto df9b = df9a.Define("SF2_MUID_Down", [fh2DL2_ID_SF1, fh2DL2_ID_SF1err](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1err);
      return f1+f2;
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   auto df9c = df9b.Define("SF2_MUISO_Up", [fh2DL2_ISO_SF1, fh2DL2_ISO_SF1err](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ISO_SF1);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ISO_SF1err);
      return f1+f2; // 
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   auto df9d = df9c.Define("SF2_MUISO_Down", [fh2DL2_ISO_SF1, fh2DL2_ISO_SF1err](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ISO_SF1);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ISO_SF1err);
      return f1-f2; // 
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   df9e = df9d.Define("SF_TRIGGER_Up", [fhL1L2_TRG_SF,fhL1L2_TRG_SFerr](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SF);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SFerr);
      return f1+f2;
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   df9f = df9e.Define("SF_TRIGGER_Down", [fhL1L2_TRG_SF,fhL1L2_TRG_SFerr](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SF);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SFerr);
      return f1-f2;
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});


   // !! Syst --  //


} else if (MuMu == 0) {
   auto df6 = dfG.Define("SF2_MUID", [fh2DL2_ID_SF1](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      return get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1);
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   auto df7 = df6.Define("SF2_MUISO", []() {
      return 1.0f; // Pas de ISO à appliquer
   });

   auto df8 = df7.Define("SF2_MURECO", [fh2DL2_Reco_SF1](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      return get2DSFFromHisto(pt[0], eta[0], fh2DL2_Reco_SF1);
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   df9 = df8.Define("SF_TRIGGER", [fhL1L2_TRG_SF](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f; 
      return get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SF);
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingpt"});

   // !! Syst
   
   auto df9a = df9.Define("SF2_MUID_Up", [fh2DL2_ID_SF1, fh2DL2_ID_SF1err](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1err);
      return f1+f2;
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   auto df9b = df9a.Define("SF2_MUID_Down", [fh2DL2_ID_SF1, fh2DL2_ID_SF1err](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_ID_SF1err);
      return f1+f2;
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   auto df9c = df9b.Define("SF2_MUISO_Up", [fh2DL2_Reco_SF1, fh2DL2_Reco_SF1err](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_Reco_SF1);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_Reco_SF1err);
      return f1+f2; // 
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   auto df9d = df9c.Define("SF2_MUISO_Down", [fh2DL2_Reco_SF1, fh2DL2_Reco_SF1err](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_Reco_SF1);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fh2DL2_Reco_SF1err);
      return f1-f2; // 
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   df9e = df9d.Define("SF_TRIGGER_Up", [fhL1L2_TRG_SF,fhL1L2_TRG_SFerr](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SF);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SFerr);
      return f1+f2;
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});

   df9f = df9e.Define("SF_TRIGGER_Down", [fhL1L2_TRG_SF,fhL1L2_TRG_SFerr](const std::vector<float>& pt, const std::vector<float>& eta) {
      if (pt.empty() || eta.empty()) return 1.0f;
      float f1 = get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SF);
      float f2 = get2DSFFromHisto(pt[0], eta[0], fhL1L2_TRG_SFerr);
      return f1-f2;
   }, {"minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2"});


   // !! -- !! //
}
   // !! Syst !!//
      auto df_10 = df9f.Define("SF_PUUp", [](const std::vector<double>& PU) {
         if (PU.empty() ) return 1.0; // Protection contre les events sans PU
         return-PU[0]; // Prend le premier muon
      }, {"miniPUweight_Up"});
      auto df_11 = df10.Define("SF_PUDown", [](const std::vector<double>& PU) {
         if (PU.empty() ) return 1.0; // Protection contre les events sans PU
         return PU[0]; // Prend le premier muon
      }, {"miniPUweight_Down"});

      auto df_11 = df10.Define("SF_PrefUp", [](const std::vector<double>& PU) {
         if (PU.empty() ) return 1.0; // Protection contre les events sans PU
         return PU[0]; // Prend le premier muon
      }, {"miniPrefweight_Up"});
      auto df_12 = df11.Define("SF_PrefDown", [](const std::vector<double>& PU) {
         if (PU.empty() ) return 1.0; // Protection contre les events sans PU
         return PU[0]; // Prend le premier muon
      }, {"miniPrefweight_Down"});

      auto df13 = df12.Define("SF_LumiUp", [LumiUp]() {
         return LumiUp;
      });
      auto df14 = df13.Define("SF_LumiDown", [LumiDown]() {
         return LumiDown;
      });


   // !! Syst !!//

      auto df_final = df14.Define("event_weight", [](double sf_pref, double sf_pu, double sf_toppt, double sf_gen,
                                                   float sf_norm, double sf_mgw, float sf_frac,
                                                   float sf_muid, float sf_muiso, float sf_mureco,
                                                   float sf2_muid, float sf2_muiso, float sf2_mureco,
                                                   float sf_trigger) {
         return sf_pref * sf_pu * sf_toppt * sf_gen *
               sf_norm * sf_mgw * sf_frac *
               sf_muid * sf_muiso * sf_mureco *
               sf2_muid * sf2_muiso * sf2_mureco *
               sf_trigger;
      }, {
         "SF_PREF", "SF_PU", "SF_TopPT", "SF_GenWeight",
         "SF_Norm", "SF_MGW", "SF_Frac",
         "SF_MUID", "SF_MUISO", "SF_MURECO",
         "SF2_MUID", "SF2_MUISO", "SF2_MURECO",
         "SF_TRIGGER"
      });

         
         // !! --- Syst --//
auto df_final_syst = df_final
.Vary("event_weight",{{"event_weight_PUUp","SF_PREF*miniPUweight_Up*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_PUDown","SF_PREF*miniPUweight_Down*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"PU_syst")
.Vary("event_weight",{{"event_weight_PrefUp","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_Prefown","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"Pref_syst")
.Vary("event_weight",{{"event_weight_LumiUp","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_LumiDown","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"Lumi_syst")
.Vary("event_weight",{{"event_weight_XSUp","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_XSDown","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"XS_syst")
.Vary("event_weight",{{"event_weight_TopPtUp","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_TopPtDown","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"TopPt_syst")
.Vary("event_weight",{{"event_weight_TriggerUp","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_TriggerDown","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"Trigger_syst")
.Vary("event_weight",{{"event_weight_MuonIDUp","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_MuonIDDown","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"MuonID_syst")
.Vary("event_weight",{{"event_weight_MuonISOUp","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_MuonISODown","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"MuonISO_syst")
.Vary("event_weight",{{"event_weight_EleIDUp","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_EleIDDown","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"EleID_syst")
.Vary("event_weight",{{"event_weight_EleRecoUp","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"},
                      {"event_weight_EleRecoDown","SF_PREF*SF_PU*SF_TopPT*SF_GenWeight*SF_Norm*SF_MGW*SF_Frac*SF_MUID*SF_MUISO*SF_MURECO*SF2_MUID*SF2_MUISO*SF2_MURECO*SF_TRIGGER"}}
                      ,"EleReco_syst")



// !! --- Syst --//
    //-------tree_Filter-------//
    // !! Becreful , everything needs to be on the same line, else it does not compile, thank you ROOT
    auto filtered_df = df_final.Filter("  ROOT::VecOps::All(minitree_Filter == true)  &&  ROOT::VecOps::All(minitree_njetNOmu > 0) && ROOT::VecOps::All(abs(minitree_Hemi_eta)<2.4) && ROOT::VecOps::All(minitree_Hemi_pt > 30) && ROOT::VecOps::All(minitree_lepton_leadingpt > 25) && ROOT::VecOps::All(minitree_lepton_leadingpt2 > 14) && ROOT::VecOps::All(minitree_Mmumu >20 )" , "Filtre global combine");

    auto Tree_filter = filtered_df.Histo1D({"Tree_filter", "Histogramme original;Pass;Counts", 2, 0, 2}, "minitree_Filter");
    Tree_filter->Write();

   // !! -- SF --//



   auto sf_PREF = filtered_df.Histo1D({"sf_PREF", "nosel;value;Counts", 200, 0, 2}, "SF_PREF");
   sf_PREF->Write();
   auto sf_PU = filtered_df.Histo1D({"sf_PU", "nosel;value;Counts", 200, 0, 2}, "SF_PU");
   sf_PU->Write();
   auto sf_TopPT = filtered_df.Histo1D({"sf_TopPT", "nosel;value;Counts", 200, 0, 2}, "SF_TopPT");
   sf_TopPT->Write();
   auto sf_GenWeight = filtered_df.Histo1D({"sf_GenWeight", "nosel;value;Counts", 500, 0, 500}, "SF_GenWeight");
   sf_GenWeight->Write();
   auto sf_Norm = filtered_df.Histo1D({"sf_Norm", "nosel;value;Counts", 100, 0, 100}, "SF_Norm");
   sf_Norm->Write();
   auto sf_MGW = filtered_df.Histo1D({"sf_MGW", "nosel;value;Counts", 200, 0, 2}, "SF_MGW");
   sf_MGW->Write();
   auto sf_Frac = filtered_df.Histo1D({"sf_Frac", "nosel;value;Counts", 200, 0, 2}, "SF_Frac");
   sf_Frac->Write();
   auto sf_MUID = filtered_df.Histo1D({"sf_MUID", "nosel;value;Counts", 200, 0, 2}, "SF_MUID");
   sf_MUID->Write();
   auto sf_MUISO = filtered_df.Histo1D({"sf_MUISO", "nosel;value;Counts", 200, 0, 2}, "SF_MUISO");
   sf_MUISO->Write();
   auto sf_MURECO = filtered_df.Histo1D({"sf_MURECO", "nosel;value;Counts", 200, 0, 2}, "SF_MURECO");
   sf_MURECO->Write();
   auto sf2_MUID = filtered_df.Histo1D({"sf2_MUID", "nosel;value;Counts", 200, 0, 2}, "SF2_MUID");
   sf2_MUID->Write();
   auto sf2_MUISO = filtered_df.Histo1D({"sf2_MUISO", "nosel;value;Counts", 200, 0, 2}, "SF2_MUISO");
   sf2_MUISO->Write();
   auto sf2_MURECO = filtered_df.Histo1D({"sf2_MURECO", "nosel;value;Counts", 200, 0, 2}, "SF2_MURECO");
   sf2_MURECO->Write();
   auto sf_TRIGGER = filtered_df.Histo1D({"sf_TRIGGER", "nosel;value;Counts", 200, 0, 2}, "SF_TRIGGER");
   sf_TRIGGER->Write();

   auto sf_event_weight= filtered_df.Histo1D({"sf_event_weight", "nosel;value;Counts", 200, 0, 20}, "event_weight");
   sf_event_weight->Write();

   // !! --------//

    auto Vertices_NoSel = filtered_df.Histo1D({"Vertices_NoSel", "h;nVtx;Counts", 100, 0, 100}, "minitree_nPV");
    Vertices_NoSel->Write();

    auto Vertices_filtercut_ = filtered_df.Histo1D({"Vertices_filtercut_", "h after filter;nVtx;Counts", 100, 0, 100}, "minitree_nPV");
    Vertices_filtercut_->Write();

    auto Tree_Mumu = filtered_df.Histo1D({"Tree_Mumu", "nosel;GeV;Counts",48,20.,500.}, "minitree_Mmumu", "event_weight");
    Tree_Mumu->Write();

    // !! ------------------

    auto Muon_pt_NoSel = filtered_df.Histo1D({"Muon_pt_NoSel", "h;GeV;Counts", 100, 0, 300}, "minitree_muon_pt");
    Muon_pt_NoSel->Write();

    auto Muon_eta_NoSel = filtered_df.Histo1D({"Muon_eta_NoSel", "h;eta;Counts", 25, -2.4, 2.4}, "minitree_muon_eta");
    Muon_eta_NoSel->Write();

    auto Electron_pt_NoSel = filtered_df.Histo1D({"Electron_pt_NoSel", "h;GeV;Counts", 100, 0, 300}, "minitree_electron_pt");
    Electron_pt_NoSel->Write();

    auto Electron_eta_NoSel = filtered_df.Histo1D({"Electron_eta_NoSel", "h;eta;Counts", 25, -2.4, 2.4}, "minitree_electron_eta");
    Electron_eta_NoSel->Write();

    // !! ------------------

    auto Leading_muon_pt_reco = filtered_df.Histo1D({"Leading_muon_pt_reco", "h;GeV;Counts", 100, 0, 400}, "minitree_lepton_leadingpt");
    Leading_muon_pt_reco->Write();

    auto Leading_lepton_pt_reco = filtered_df.Histo1D({"Leading_lepton_pt_reco", "h;GeV;Counts", 100, 0, 400}, "minitree_lepton_leadingpt2");
    Leading_lepton_pt_reco->Write();

    auto Leading_muon_eta_reco = filtered_df.Histo1D({"Leading_muon_eta_reco", "h;eta;Counts", 25, -2.4, 2.4}, "minitree_lepton_leadingeta");
    Leading_muon_eta_reco->Write();

    auto Leading_lepton_eta_reco = filtered_df.Histo1D({"Leading_lepton_eta_reco", "h;eta;Counts", 25, -2.4, 2.4}, "minitree_lepton_leadingeta2");
    Leading_lepton_eta_reco->Write();

    // !! ------------------

    auto Leading_muon_phi_reco = filtered_df.Histo1D({"Leading_muon_phi_reco", "h;phi;Counts", 25, -3.14, 3.14}, "minitree_lepton_leadingphi");
    Leading_muon_phi_reco->Write();

    auto Leading_lepton_phi_reco = filtered_df.Histo1D({"Leading_lepton_phi_reco", "h;phi;Counts", 25, -3.14, 3.14}, "minitree_lepton_leadingphi2");
    Leading_lepton_phi_reco->Write();

    auto Leading_lepton_pt_eta_reco = filtered_df.Histo2D({"Leading_lepton_pt_eta_reco", "h;pt;eta;Counts", 40, 0, 200, 25, -2.4, 2.4}, "minitree_lepton_leadingpt2", "minitree_lepton_leadingeta2");
    Leading_lepton_pt_eta_reco->Write();

    auto Leading_muon_dxy_reco = filtered_df.Histo1D({"Leading_muon_dxy_reco", "h;cm;Counts", 40, -0.2, 0.2}, "minitree_lepton_leadingdxy");
    Leading_muon_dxy_reco->Write();
  


    // !! ------------------

    auto Leading_lepton_dxy_reco = filtered_df.Histo1D({"Leading_lepton_dxy_reco", "h;cm;Counts", 200, -1, 1}, "minitree_lepton_leadingdxy2");
    Leading_lepton_dxy_reco->Write();

    auto Leading_muon_dz_reco = filtered_df.Histo1D({"Leading_muon_dz_reco", "h;cm;Counts", 50, -0.5, 0.5}, "minitree_lepton_leadingdz");
    Leading_muon_dz_reco->Write();

    auto Leading_lepton_dz_reco = filtered_df.Histo1D({"Leading_lepton_dz_reco", "h;cm;Counts", 200, -1, 1}, "minitree_lepton_leadingdz2");
    Leading_lepton_dz_reco->Write();

    auto Njet_NoSel = filtered_df.Histo1D({"Njet_NoSel", "h;Counts;Counts", 20, 0, 20}, "minitree_njet");
    Njet_NoSel->Write();
  
    // !! ------------------

    auto NjetNOmu_NoSel = filtered_df.Histo1D({"NjetNOmu_NoSel", "h;Counts;Counts", 20, 0, 20}, "minitree_njetNOmu");
    NjetNOmu_NoSel->Write();

    auto Hemisphere_leadingpt = filtered_df.Histo1D({"Hemisphere_leadingpt", "h;GeV;Counts", 100, 0, 500}, "minitree_Hemi_pt");// !! à changer ot leading 
    Hemisphere_leadingpt->Write();

    auto Hemisphere_subleadingpt = filtered_df.Histo1D({"Hemisphere_subleadingpt", "h;GeV;Counts", 100, 0, 500}, "minitree_Hemi_pt");// !! à changer to subleading
    Hemisphere_subleadingpt->Write();

    auto HData_jet_pt = filtered_df.Histo1D({"HData_jet_pt", "h;GeV;Counts", 100, 0, 500}, "minitree_jet_pt");
    HData_jet_pt->Write();
  
    // !! ------------------

    auto HData_jet_eta = filtered_df.Histo1D({"HData_jet_eta", "h;eta;Counts", 55, -5, 5}, "minitree_jet_eta");
    HData_jet_eta->Write();

    // auto HData_jet_btag_Deepjet = filtered_df.Histo1D({"HData_jet_btag_Deepjet", "h;Counts;Counts", 50, 0, 1}, "minitree_jet_btag_DeepJet");
    // HData_jet_btag_Deepjet->Write();

    // auto HData_jet_HadronFlavour = filtered_df.Histo1D({"HData_jet_HadronFlavour", "h;Counts;Counts", 7, -0.5, 6.5}, "minitree_jet_HadronFlavour");
    // HData_jet_HadronFlavour->Write();

    auto Leading_jet_pt = filtered_df.Histo1D({"Leading_jet_pt", "h;GeV;Counts", 100, 0, 500}, "minitree_jet_leadingpt");
    Leading_jet_pt->Write();

  
    // !! ------------------

    auto Leading_jet_eta = filtered_df.Histo1D({"Leading_jet_eta", "h;eta;Counts", 55, -2.5, 2.5}, "minitree_jet_leadingeta");
    Leading_jet_eta->Write();

    auto Subleading_jet_pt = filtered_df.Histo1D({"Subleading_jet_pt", "h;GeV;Counts", 100, 0, 500}, "minitree_jet_leadingpt2");
    Subleading_jet_pt->Write();

    auto Subleading_jet_eta = filtered_df.Histo1D({"Subleading_jet_eta", "h;eta;Counts", 55, -2.5, 2.5}, "minitree_jet_leadingeta2");
    Subleading_jet_eta->Write();

    auto Nmu = filtered_df.Histo1D({"Nmu", "h;Counts;Counts", 10, 0, 10}, "minitree_nmu");
    Nmu->Write();

  
    // !! ------------------


    auto Muon_pt = filtered_df.Histo1D({"Muon_pt", "h;GeV;Counts", 100, 0, 500}, "minitree_muon_pt");
    Muon_pt->Write();

    // auto Muon_PFIsoLoose = filtered_df.Histo1D({"Muon_PFIsoLoose", "h;Counts;Counts", 2, 0, 2}, "minitree_muon_PFIsoLoose");
    // Muon_PFIsoLoose->Write();

    // auto Muon_MiniIsoTight = filtered_df.Histo1D({"Muon_MiniIsoTight", "h;Counts;Counts", 2, 0, 2}, "minitree_muon_MiniIsoTight");
    // Muon_MiniIsoTight->Write();

    auto LeadingLeptons_dR = filtered_df.Histo1D({"LeadingLeptons_dR", "h;Counts;Counts", 50, 0, 5}, "minitree_lepton_lepton_dR");
    LeadingLeptons_dR->Write();

    // !! ------------------

    auto LeadingLeptons_dPhi = filtered_df.Histo1D({"LeadingLeptons_dPhi", "h;Counts;Counts", 35, 0, 3.5}, "minitree_lepton_lepton_dPhi");
    LeadingLeptons_dPhi->Write();

    auto LeadingJets_dR = filtered_df.Histo1D({"LeadingJets_dR", "h;Counts;Counts", 50, 0, 5}, "minitree_jet_jet_dR");
    LeadingJets_dR->Write();

    auto LeadingJets_dPhi = filtered_df.Histo1D({"LeadingJets_dPhi", "h;Counts;Counts", 35, 0, 3.5}, "minitree_jet_jet_dPhi");
    LeadingJets_dPhi->Write();

    auto LeadingLeptonJet_dRmax = filtered_df.Histo1D({"LeadingLeptonJet_dRmax", "h;Counts;Counts", 50, 0, 5}, "minitree_muon_jet_dRmax");
    LeadingLeptonJet_dRmax->Write();

    // !! ------------------

    auto LeadingLeptonJet_dRmin = filtered_df.Histo1D({"LeadingLeptonJet_dRmin", "h;Counts;Counts", 50, 0, 5}, "minitree_muon_jet_dRmin");
    LeadingLeptonJet_dRmin->Write();

    auto HemiAxis_Mu_dR = filtered_df.Histo1D({"HemiAxis_Mu_dR", "h;Counts;Counts", 50, 0, 5}, "minitree_HemiMu_dR");
    HemiAxis_Mu_dR->Write();

    auto HemiAxis_OpMu_dR = filtered_df.Histo1D({"HemiAxis_OpMu_dR", "h;Counts;Counts", 50, 0, 5}, "minitree_HemiMuOp_dR");
    HemiAxis_OpMu_dR->Write();

    auto HT = filtered_df.Histo1D({"HT", "h;GeV;Counts", 200, 0, 800}, "minitree_HT");
    HT->Write();

    // !! ------------------

    auto LT = filtered_df.Histo1D({"LT", "h;GeV;Counts", 100, 0, 500}, "minitree_LT");
    LT->Write();

    auto track_nTracks = filtered_df.Histo1D({"track_nTracks", "h;Counts;Counts", 100, 0, 100}, "minitree_nTracks");
    track_nTracks->Write();

    auto track_nLostTracks = filtered_df.Histo1D({"track_nLostTracks", "h;Counts;Counts", 25, 0, 25}, "minitree_nLostTracks");
    track_nLostTracks->Write();

    auto Hemi_pt = filtered_df.Histo1D({"Hemi_pt", "h;GeV;Counts", 100, 0, 600}, "minitree_Hemi_pt");
    Hemi_pt->Write();
  
    // !! ------------------


    auto Hemi_eta = filtered_df.Histo1D({"Hemi_eta", "h;eta;Counts", 55, -2.5, 2.5}, "minitree_Hemi_eta");
    Hemi_eta->Write();

    auto Hemi_nJet = filtered_df.Histo1D({"Hemi_nJet", "h;Counts;Counts", 10, 0, 10}, "minitree_Hemi_njet");
    Hemi_nJet->Write();

    auto Hemi_nJetNoMu = filtered_df.Histo1D({"Hemi_nJetNoMu", "h;Counts;Counts", 10, 0, 10}, "minitree_Hemi_njet_nomu");
    Hemi_nJetNoMu->Write();

    auto HemiMu_pt = filtered_df.Histo1D({"HemiMu_pt", "h;GeV;Counts", 100, 0, 600}, "minitree_HemiMu_pt");
    HemiMu_pt->Write();
  
    // !! ------------------

    auto HemiMu_Mass = filtered_df.Histo1D({"HemiMu_Mass", "h;GeV;Counts", 100, 0, 600}, "minitree_HemiMu_mass");
    HemiMu_Mass->Write();

    auto Hemi_nTrks = filtered_df.Histo1D({"Hemi_nTrks", "h;Counts;Counts", 25, 0, 25}, "minitree_Hemi_nTrks");
    Hemi_nTrks->Write();

    auto Hemi_Mass = filtered_df.Histo1D({"Hemi_Mass", "h;GeV;Counts", 50, 0, 300}, "minitree_Hemi_mass");
    Hemi_Mass->Write();

    auto track_TRACK_SIZE = filtered_df.Histo1D({"track_TRACK_SIZE", "h;Counts;Counts", 5001, -0.5, 5000.5}, "minitree_TRACK_SIZE");
    track_TRACK_SIZE->Write();
  
    // !! ------------------

    auto track_ipc = filtered_df.Histo1D({"track_ipc", "h;Counts;Counts", 50, 0, 5000}, "minitree_track_ipc");
    track_ipc->Write();

    auto track_dxyError = filtered_df.Histo1D({"track_dxyError", "h;Counts;Counts", 41, -20, 20}, "minitree_track_dxyError");
    track_dxyError->Write();
  
    // !! ------------------

    auto track_nHitPixel = filtered_df.Histo1D({"track_nHitPixel", "h;Counts;Counts", 15, 0, 15}, "minitree_track_nHitPixel");
    track_nHitPixel->Write();

    auto track_nHitTIB = filtered_df.Histo1D({"track_nHitTIB", "h;Counts;Counts", 15, 0, 15}, "minitree_track_nHitTIB");
    track_nHitTIB->Write();

    auto track_nHitTOB = filtered_df.Histo1D({"track_nHitTOB", "h;Counts;Counts", 15, 0, 15}, "minitree_track_nHitTOB");
    track_nHitTOB->Write();

    auto track_nHitTEC = filtered_df.Histo1D({"track_nHitTEC", "h;Counts;Counts", 15, 0, 15}, "minitree_track_nHitTEC");
    track_nHitTEC->Write();
  
    // !! ------------------

    auto track_nHitPXB = filtered_df.Histo1D({"track_nHitPXB", "h;Counts;Counts", 15, 0, 15}, "minitree_track_nHitPXB");
    track_nHitPXB->Write();

    auto track_nHitPXF = filtered_df.Histo1D({"track_nHitPXF", "h;Counts;Counts", 15, 0, 15}, "minitree_track_nHitPXF");
    track_nHitPXF->Write();

    auto track_isHitPixel = filtered_df.Histo1D({"track_isHitPixel", "h;Counts;Counts", 3000, 0, 1500}, "minitree_track_isHitPixel");
    track_isHitPixel->Write();

    auto track_nLayers = filtered_df.Histo1D({"track_nLayers", "h;Counts;Counts", 80, 0, 20}, "minitree_track_nLayers");
    track_nLayers->Write();
  
    // !! ------------------

    auto track_nLayersPixel = filtered_df.Histo1D({"track_nLayersPixel", "h;Counts;Counts", 80, 0, 20}, "minitree_track_nLayersPixel");
    track_nLayersPixel->Write();

    auto track_region = filtered_df.Histo1D({"track_region", "h;Counts;Counts", 6, -0.5, 5.5}, "minitree_track_region");
    track_region->Write();

    // auto track_btag = filtered_df.Histo1D({"track_btag", "h;Counts;Counts", 22, -0.5, 10.5}, "minitree_track_btag");
    // track_btag->Write();

    // auto track_Hemi = filtered_df.Histo1D({"track_Hemi", "h;Counts;Counts", 6, -0.5, 5.5}, "minitree_track_Hemi");
    // track_Hemi->Write();
  
    // // !! ------------------

    // auto track_lost = filtered_df.Histo1D({"track_lost", "h;Counts;Counts", 2, 0, 2}, "minitree_track_lost");
    // track_lost->Write();

    // auto track_dxy = filtered_df.Histo1D({"track_dxy", "h;Counts;Counts", 100, -50, 50}, "minitree_track_dxy");
    // track_dxy->Write();

    // auto track_dz = filtered_df.Histo1D({"track_dz", "h;Counts;Counts", 200, -100, 100}, "minitree_track_dz");
    // track_dz->Write();
  
    // // !! ------------------

    // auto track_pt = filtered_df.Histo1D({"track_pt", "h;Counts;Counts", 300, 0, 300}, "minitree_track_pt");
    // track_pt->Write();

    // auto track_eta = filtered_df.Histo1D({"track_eta", "h;Counts;Counts", 80, -4, 4}, "minitree_track_eta");
    // track_eta->Write();

    // auto track_NChi2 = filtered_df.Histo1D({"track_NChi2", "h;Counts;Counts", 5, 0, 5}, "minitree_track_NChi2");
    // track_NChi2->Write();

    // auto track_nhits = filtered_df.Histo1D({"track_nhits", "h;Counts;Counts", 40, 0, 40}, "minitree_track_nHit");
    // track_nhits->Write();
  
    // // !! ------------------

    // auto track_iJet = filtered_df.Histo1D({"track_iJet", "h;Counts;Counts", 22, -2, 20}, "minitree_track_iJet");
    // track_iJet->Write();

    // auto track_drSig = filtered_df.Histo1D({"track_drSig", "h;Counts;Counts", 1000, 0, 1000}, "minitree_track_drSig");
    // track_drSig->Write();

    // auto track_dzSig = filtered_df.Histo1D({"track_dzSig", "h;Counts;Counts", 5000, 0, 5000}, "minitree_track_dzSig");
    // track_dzSig->Write();

    // auto track_ntrk10 = filtered_df.Histo1D({"track_ntrk10", "h;Counts;Counts", 100, 0, 100}, "minitree_track_ntrk10");
    // track_ntrk10->Write();
  
    // // !! ------------------

    // auto track_ntrk20 = filtered_df.Histo1D({"track_ntrk20", "h;Counts;Counts", 100, 0, 100}, "minitree_track_ntrk20");
    // track_ntrk20->Write();

    // auto track_ntrk30 = filtered_df.Histo1D({"track_ntrk30", "h;Counts;Counts", 100, 0, 100}, "minitree_track_ntrk30");
    // track_ntrk30->Write();

    // auto track_ntrk40 = filtered_df.Histo1D({"track_ntrk40", "h;Counts;Counts", 100, 0, 100}, "minitree_track_ntrk40");
    // track_ntrk40->Write();

    // auto track_track_Hemi_dR = filtered_df.Histo1D({"track_track_Hemi_dR", "h;Counts;Counts", 50, 0, 5}, "minitree_track_Hemi_dR");
    // track_track_Hemi_dR->Write();
  
    // // !! ------------------

    // auto track_track_Hemi_dRmax = filtered_df.Histo1D({"track_track_Hemi_dRmax", "h;Counts;Counts", 50, 0, 5}, "minitree_track_Hemi_dRmax");
    // track_track_Hemi_dRmax->Write();

    // auto track_MVAVAL_TRK = filtered_df.Histo1D({"track_MVAVAL_TRK", "h;Counts;Counts", 100, -1, 1}, "minitree_track_MVAval");
    // track_MVAVAL_TRK->Write();

    // auto track_Track_firstHit = filtered_df.Histo1D({"track_Track_firstHit", "h;Counts;Counts", 4000, 0, 4000}, "minitree_track_firstHit");
    // track_Track_firstHit->Write();

    // auto track_Track_firstHit_x = filtered_df.Histo1D({"track_Track_firstHit_x", "h;Counts;Counts", 201, -100, 100}, "minitree_track_firstHit_x");
    // track_Track_firstHit_x->Write();

    // // !! ------------------

    // auto track_Track_firstHit_y = filtered_df.Histo1D({"track_Track_firstHit_y", "h;Counts;Counts", 201, -100, 100}, "minitree_track_firstHit_y");
    // track_Track_firstHit_y->Write();

    // // auto track_Track_firstHit_r = filtered_df.Histo1D({"track_Track_firstHit_r", "h;Counts;Counts", 600, 0, 200}, "minitree_track_firstHit_r");
    // // track_Track_firstHit_r->Write();

    // auto track_Track_firstHit_z = filtered_df.Histo1D({"track_Track_firstHit_z", "h;Counts;Counts", 801, -200, 200}, "minitree_track_firstHit_z");
    // track_Track_firstHit_z->Write();

    // auto track_Track_firstHit_X_Y = filtered_df.Histo2D({"track_Track_firstHit_X_Y", "h;Counts;Counts", 201, -100, 100, 201, -100, 100}, "minitree_track_firstHit_x", "minitree_track_firstHit_y");
    // track_Track_firstHit_X_Y->Write();

    // !! ------------------

    auto K0_mass = filtered_df.Histo1D({"K0_mass", "h;GeV;Counts", 202, 0.42, 0.58}, "minitree_K0_mass");
    K0_mass->Write();

    auto K0_pt = filtered_df.Histo1D({"K0_pt", "h;Counts;Counts", 50, 0, 50}, "minitree_K0_pt");
    K0_pt->Write();

    auto L0_mass = filtered_df.Histo1D({"L0_mass", "h;GeV;Counts", 202, 1.06, 1.18}, "minitree_L0_mass");
    L0_mass->Write();

    auto L0_pt = filtered_df.Histo1D({"L0_pt", "h;Counts;Counts", 500, 0, 50}, "minitree_L0_pt");
    L0_pt->Write();

    // !! ------------------


    // auto K0_df = filtered_df.Filter("(minitree_V0_reco_source ==1) ", "K0 Filter");
    // auto L0_df = filtered_df.Filter("(minitree_V0_reco_source ==2) ", "L0 Filter");

    // auto Reco_K0_mass = K0_df.Histo1D({"Reco_K0_mass", "h;GeV;Counts", 202, 0.42, 0.58}, "minitree_Reco_K0_mass");
    // Reco_K0_mass->Write();

    // auto Reco_K0_pt = K0_df.Histo1D({"Reco_K0_pt", "h;Counts;Counts", 50, 0, 50}, "minitree_Reco_K0_pt");
    // Reco_K0_pt->Write();

    // auto Reco_L0_mass = L0_df.Histo1D({"Reco_L0_mass", "h;GeV;Counts", 202, 1.06, 1.18}, "minitree_Reco_L0_mass");
    // Reco_L0_mass->Write();

    // auto Reco_L0_pt = L0_df.Histo1D({"Reco_L0_pt", "h;Counts;Counts", 50, 0, 50}, "minitree_Reco_L0_pt");
    // Reco_L0_pt->Write();

    // !! ------------------

    auto SecInt_mass = filtered_df.Histo1D({"SecInt_mass", "h;Counts;Counts", 20, 0, 2}, "minitree_SecInt_mass");
    SecInt_mass->Write();

    auto SecInt_drSig = filtered_df.Histo1D({"SecInt_drSig", "h;Counts;Counts", 200, 0, 2000}, "minitree_SecInt_drSig");
    SecInt_drSig->Write();

    auto SecInt_pt = filtered_df.Histo1D({"SecInt_pt", "h;Counts;Counts", 100, 0, 100}, "minitree_SecInt_pt");
    SecInt_pt->Write();

    auto SecInt_dzSig = filtered_df.Histo1D({"SecInt_dzSig", "h;Counts;Counts", 200, 0, 2000}, "minitree_SecInt_dzSig");
    SecInt_dzSig->Write();

    // !! ------------------

    auto SecInt_selec = filtered_df.Histo1D({"SecInt_selec", "h;Counts;Counts", 5, 0, 5}, "minitree_SecInt_selec");
    SecInt_selec->Write();

    auto SecInt_layer = filtered_df.Histo1D({"SecInt_layer", "h;Counts;Counts", 100, 0, 50}, "minitree_SecInt_layer");
    SecInt_layer->Write();

    auto SecInt_r = filtered_df.Histo1D({"SecInt_r", "h;Counts;Counts", 200, 0, 20}, "minitree_SecInt_r");
    SecInt_r->Write();

    auto SecInt_z = filtered_df.Histo1D({"SecInt_z", "h;Counts;Counts", 800, -200, 200}, "minitree_SecInt_z");
    SecInt_z->Write();

    // !! ------------------


    // auto SecInt_df = filtered_df.Filter("miniminitree_SecInt_selec != 0", "SecInt Selec Filter");

    // auto SecInt_mass_TrackerMatched = SecInt_df.Histo1D({"SecInt_mass_TrackerMatched", "h;Counts;Counts", 20, 0, 2}, "minitree_SecInt_mass");
    // SecInt_mass->Write();

    // auto SecInt_drSig_TrackerMatched = SecInt_df.Histo1D({"SecInt_drSig_TrackerMatched", "h;Counts;Counts", 200, 0, 2000}, "minitree_SecInt_drSig");
    // SecInt_drSig->Write();

    // auto SecInt_pt_TrackerMatched = SecInt_df.Histo1D({"SecInt_pt_TrackerMatched", "h;Counts;Counts", 100, 0, 100}, "minitree_SecInt_pt");
    // SecInt_pt->Write();

    // auto SecInt_dzSig_TrackerMatched = SecInt_df.Histo1D({"SecInt_dzSig_TrackerMatched", "h;Counts;Counts", 200, 0, 2000}, "minitree_SecInt_dzSig");
    // SecInt_dzSig->Write();

    // // !! ------------------
    // auto TMSecInt_df = SecInt_df.Filter("miniminitree_SecInt_layer != 0", "SecInt TM Filter");

    // auto SecInt_selec_TrackerMatched = TMSecInt_df.Histo1D({"SecInt_selec_TrackerMatched", "h;Counts;Counts", 5, 0, 5}, "minitree_SecInt_selec");
    // SecInt_selec->Write();

    // auto SecInt_layer_TrackerMatched = TMSecInt_df.Histo1D({"SecInt_layer_TrackerMatched", "h;Counts;Counts", 100, 0, 50}, "minitree_SecInt_layer");
    // SecInt_layer->Write();

    // auto SecInt_r_TrackerMatched = TMSecInt_df.Histo1D({"SecInt_r_TrackerMatched", "h;Counts;Counts", 200, 0, 20}, "minitree_SecInt_r");
    // SecInt_r->Write();

    // auto SecInt_z_TrackerMatched = TMSecInt_df.Histo1D({"SecInt_z_TrackerMatched", "h;Counts;Counts", 800, -200, 200}, "minitree_SecInt_z");
    // SecInt_z->Write();

    // !! ------------------

    auto SecVtx_NChi2 = filtered_df.Histo1D({"SecVtx_NChi2", "h;Counts;Counts", 15, 0, 15}, "minitree_Hemi_SecVtx_NChi2");
    SecVtx_NChi2->Write();

    auto SecVtx_nTrks = filtered_df.Histo1D({"SecVtx_nTrks", "h;Counts;Counts", 40, 0, 40}, "minitree_Hemi_SecVtx_nTrks");
    SecVtx_nTrks->Write();

    auto SecVtx_Mass = filtered_df.Histo1D({"SecVtx_Mass", "h;Counts;Counts", 150, 0, 1500}, "minitree_Hemi_SecVtx_Mass");
    SecVtx_Mass->Write();

    auto SecVtx_Dist = filtered_df.Histo1D({"SecVtx_Dist", "h;Counts;Counts", 100, 0, 100}, "minitree_Hemi_SecVtx_dist");
    SecVtx_Dist->Write();

    // !! ------------------

    auto Vtx_NChi2 = filtered_df.Histo1D({"Vtx_NChi2", "h;Counts;Counts", 15, 0, 15}, "minitree_Hemi_Vtx_NChi2");
    Vtx_NChi2->Write();

    auto Vtx_nTrks = filtered_df.Histo1D({"Vtx_nTrks", "h;Counts;Counts", 40, 0, 40}, "minitree_Hemi_Vtx_nTrks");
    Vtx_nTrks->Write();

    auto Vtx_Mass = filtered_df.Histo1D({"Vtx_Mass", "h;Counts;Counts", 10, 0, 100}, "minitree_Hemi_Vtx_Mass");
    Vtx_Mass->Write();

    auto Vtx_Dist = filtered_df.Histo1D({"Vtx_Dist", "h;Counts;Counts", 20, 0, 100}, "minitree_Hemi_Vtx_dist");
    Vtx_Dist->Write();

    // !! ------------------

    auto FinalVtx_nTrks = filtered_df.Histo1D({"FinalVtx_nTrks", "h;Counts;Counts", 40, 0, 40}, "minitree_Hemi_Vtx_BDT_nTrks");
    FinalVtx_nTrks->Write();

    auto FinalVtx_Step = filtered_df.Histo1D({"FinalVtx_Step", "h;Counts;Counts", 4, 1, 5}, "minitree_Hemi_Vtx_BDT_step");
    FinalVtx_Step->Write();

    auto FinalVtx_Mass = filtered_df.Histo1D({"FinalVtx_Mass", "h;Counts;Counts", 100, 0, 500}, "minitree_Hemi_Vtx_BDT_HMass");
    FinalVtx_Mass->Write();

    auto FinalVtx_HMass = filtered_df.Histo1D({"FinalVtx_HMass", "h;Counts;Counts", 100, 0, 500}, "minitree_Hemi_Vtx_BDT_HMass");
    FinalVtx_HMass->Write();

    // !! ------------------

    auto Vtx_Step = filtered_df.Histo1D({"Vtx_Step", "h;Counts;Counts", 4, 1, 5}, "minitree_Hemi_Vtx_step");
    Vtx_Step->Write();

    auto Vtx_dR = filtered_df.Histo1D({"Vtx_dR", "h;Counts;Counts", 50, 0, 5}, "minitree_Hemi_Vtx_dR");
    Vtx_dR->Write();

    auto SecVtx_Step = filtered_df.Histo1D({"SecVtx_Step", "h;Counts;Counts", 4, 1, 5}, "minitree_Hemi_SecVtx_step");
    SecVtx_Step->Write();

    auto SecVtx_dR = filtered_df.Histo1D({"SecVtx_dR", "h;Counts;Counts", 50, 0, 5}, "minitree_Hemi_SecVtx_dR");
    SecVtx_dR->Write();

    std::cout << "Analyse terminee. TTree sauvegarde dans " << output_file << std::endl;
    output.Close();

    return 0;
}

