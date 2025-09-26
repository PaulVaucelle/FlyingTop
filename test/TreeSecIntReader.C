#define TreeSecIntReader_cxx
#include "TreeSecIntReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>
 
#include "TTree.h"
#include "../HistogramManager.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TH1F.h"
using namespace std;

void TreeSecIntReader::Loop(TString Prod, TString sample)
{
//   In a ROOT session, you can do:
//      root> .L TreeSecIntReader.C
//      root> TreeSecIntReader t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//
    bool firstinit = true;
     //cout << " iiii " << i << endl;
     //cout << samplename << endl;

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TFile * theoutputfile = new TFile( ("./SECINT_"+sample+".root").Data() , "RECREATE");
   // std::ofstream ofs ("PROD_CSI_10_06_2024/Efficacity_"+thesample+".txt", std::ofstream::out);

   theoutputfile->cd();
        initializeHisto(sample, firstinit);
   std::cout << "Current directory: " << gDirectory->GetName() << std::endl;
   Long64_t nentries = fChain->GetEntries();
   std::cout<<"nentries = "<<nentries<<std::endl;

   float TOTALnSecIntFullSelec = 0;
   float TOTALnSecIntTrackerMatched = 0;


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {//nentries
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if ( jentry%10000 == 0 ) std::cout << "events : " << jentry << std::endl;

      int MNEU = 0;
      int MSMUON = 0;
      float CTAU = 0;
      if (minitree_smu_mass->size() > 0) //for signal samples basically
        {
            MNEU = minitree_neu_mass->at(0);
            MSMUON = minitree_smu_mass->at(0);
            CTAU = minitree_neu_ctau->at(0);
        }


      for (unsigned int i = 0; i < minitree_K0_reco_mass->size(); i++)
         {
            fillHisto("hData_reco_K0_mass","",  sample, minitree_K0_reco_mass->at(i),1.); 
         }
      for (unsigned int i = 0; i < minitree_L0_reco_mass->size(); i++)
         {
            fillHisto("hData_reco_L0_mass","",  sample, minitree_L0_reco_mass->at(i),1.);
         }

      for (unsigned int iV0 = 0 ; iV0 < minitree_V0_reco_source->size(); iV0++)
         {

            if ( abs(minitree_V0_reco_z->at(iV0)) <27 && abs(minitree_V0_reco_x->at(iV0))<25 && abs(minitree_V0_reco_y->at(iV0))<25 && abs(minitree_V0_reco_eta->at(iV0)) < 1.4)
               {
                  fillHisto2D("hData_reco_V0_xy","", sample, minitree_V0_reco_x->at(iV0),minitree_V0_reco_y->at(iV0),1.);
               }
               
            if (abs(minitree_V0_reco_z->at(iV0)) < 120  && abs(minitree_V0_reco_r->at(iV0))<70) 
               {
                  fillHisto2D("hData_reco_V0_rz","", sample, abs(minitree_V0_reco_z->at(iV0)), minitree_V0_reco_r->at(iV0),1.);
               }
         }


      for (unsigned int iK0 = 0 ; iK0 < minitree_K0_x->size(); iK0++)
         {
            if ( abs(minitree_K0_z->at(iK0)) <27 && abs(minitree_K0_x->at(iK0))<25 && abs(minitree_K0_y->at(iK0))<25 && abs(minitree_K0_eta->at(iK0)) < 1.4)
               {
                  fillHisto2D("hData_CMSSW_V0_xy","", sample, minitree_K0_x->at(iK0),minitree_K0_y->at(iK0),1.);
               }
               
            if (abs(minitree_K0_z->at(iK0)) < 120  && abs(minitree_K0_r->at(iK0))<70) 
               {
                  fillHisto2D("hData_CMSSW_V0_rz","", sample, abs(minitree_K0_z->at(iK0)), minitree_K0_r->at(iK0),1.);
               }
         }

      for (unsigned int iL0 = 0 ; iL0 < minitree_L0_x->size(); iL0++)
         {
            if ( abs(minitree_K0_z->at(iL0)) <27 && abs(minitree_K0_x->at(iL0))<25 && abs(minitree_K0_y->at(iL0))<25 && abs(minitree_K0_eta->at(iL0)) < 1.4)
               {
                  fillHisto2D("hData_CMSSW_V0_xy","", sample, minitree_K0_x->at(iL0),minitree_K0_y->at(iL0),1.);
               }
               
            if (abs(minitree_K0_z->at(iL0)) < 120  && abs(minitree_K0_r->at(iL0))<70) 
               {
                  fillHisto2D("hData_CMSSW_V0_rz","", sample, abs(minitree_K0_z->at(iL0)), minitree_K0_r->at(iL0),1.);
               }
         }

      for (unsigned int iYc = 0 ; iYc < minitree_Yc_x->size(); iYc++)
         {
            if ( abs(minitree_Yc_z->at(iYc)) <27 && abs(minitree_Yc_x->at(iYc))<25 && abs(minitree_Yc_y->at(iYc))<25 && abs(minitree_Yc_eta->at(iYc)) < 1.4)
               {
                  fillHisto2D("hData_CMSSW_Yc_xy","", sample, minitree_Yc_x->at(iYc),minitree_Yc_y->at(iYc),1.);
               }
               
            if (abs(minitree_Yc_z->at(iYc)) < 120  && abs(minitree_Yc_r->at(iYc))<70) 
               {
                  fillHisto2D("hData_CMSSW_Yc_rz","", sample, abs(minitree_Yc_z->at(iYc)), minitree_Yc_r->at(iYc),1.);
               }
         }
//   smalltree->Branch("minitree_Yc_layer",&minitree_Yc_layer);



      //*******************************
      //loop on Sec. Interactions
      //*******************************
      

      int nSecIntFullSelec = 0;
      int nSecIntTrackerMatched = 0;
      int nSecIntFullSelec_PU25 = 0;
      int nSecIntTrackerMatched_PU25 = 0;
      int nSecIntFullSelec_PU30 = 0;
      int nSecIntTrackerMatched_PU30 = 0;
      int nSecIntFullSelec_PU35 = 0;
      int nSecIntTrackerMatched_PU35 = 0;
      int nSecIntFullSelec_PU40 = 0;
      int nSecIntTrackerMatched_PU40 = 0;
      int nSecIntFullSelec_PU45 = 0;
      int nSecIntTrackerMatched_PU45 = 0;
      int nSecIntFullSelec_PU50 = 0;
      int nSecIntTrackerMatched_PU50 = 0;

      for (unsigned int iSecInt = 0; iSecInt <minitree_SecInt_mass->size(); iSecInt++)
         {

            if (minitree_SecInt_selec->at(iSecInt))
               {
                  fillHisto("hData_reco_SecInt_mass","FullSelec",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                  nSecIntFullSelec++;
                  TOTALnSecIntFullSelec++;
                  if(minitree_SecInt_layer->at(iSecInt)) 
                     {
                        fillHisto("hData_reco_SecInt_mass","FullTrackerMatched",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                        if (minitree_smu_mass->size() > 0)
                           {
                              fillHisto2D("hData_reco_3DSecInt","", sample, MNEU,CTAU,1.);
                           }
                        nSecIntTrackerMatched++;
                        TOTALnSecIntTrackerMatched++;
                     }

                  //to get a nice view of the inner tracker x vs y
                  //Selected
                  if ( abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4&&abs(minitree_SecInt_z->at(iSecInt))<27)
                    {
                        float R = sqrt(minitree_SecInt_x->at(iSecInt)*minitree_SecInt_x->at(iSecInt)+minitree_SecInt_y->at(iSecInt)*minitree_SecInt_y->at(iSecInt));
                        if (R > 2 && tree_SecInt_dca->at(i) < 1.)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","Selec", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","Selec",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","Selec", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);

                              fillHisto2D("hData_reco_SecInt_Sigxy_xy","Selec", sample, minitree_SecInt_drSig->at(iSecInt),minitree_SecInt_r->at(iSecInt),1.);
                              if (minitree_SecInt_drSig->at(iSecInt)<1000   )
                                 {
                                    fillHisto2D("hData_reco_SecInt_xy","Selec_Spe1000", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                                 }
                              if (minitree_SecInt_drSig->at(iSecInt)<9000   )
                                 {
                                    fillHisto2D("hData_reco_SecInt_xy","Selec_Spe9000", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                                 }
                           }

                    }
                  //TrackerMatched
                  if(minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4 && abs(minitree_SecInt_z->at(iSecInt))<27)
                    {
                      fillHisto2D("hData_reco_SecInt_xy","TrackerMatched", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                      fillHisto("hData_reco_SecInt_mass","TrackerMatched",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                      fillHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                    }

                  // to get a nice view of the tracker in the r vs z plane
                  //No PileUp selection
                  //Selected
                  if ( abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                    {
                      fillHisto2D("hData_reco_SecInt_rz","Selec", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                      fillHisto2D("hData_reco_SecInt_rz_Inner","Selec", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                    }
                  //
                  if ( minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                    {
                      fillHisto2D("hData_reco_SecInt_rz","TrackerMatched", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                      fillHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                    }

                  if (minitree_tree_nPV->at(0) < 25)
                     {
                        nSecIntFullSelec_PU25++;
                        if (minitree_SecInt_layer->at(iSecInt)!=0)
                           {
                              nSecIntTrackerMatched_PU25++;
                           }
                        if ( abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4&&abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","Selec_PU25", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","Selec_PU25",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU25", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                           }
                           //TrackerMatched
                        if(minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4 && abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU25", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","TrackerMatched_PU25",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU25", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                           }

                           // to get a nice view of the tracker in the r vs z plane
                           //No PileUp selection
                           //Selected
                        if ( abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","Selec_PU25", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU25", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                           }
                           //
                        if ( minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU25", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU25", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                           }
                     }

                  if (minitree_tree_nPV->at(0) < 30)
                     {
                        nSecIntFullSelec_PU30++;
                        if (minitree_SecInt_layer->at(iSecInt)!=0)
                           {
                              nSecIntTrackerMatched_PU30++;
                           }
                        if ( abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4&&abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","Selec_PU30", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","Selec_PU30",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU30", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                           }
                           //TrackerMatched
                        if(minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4 && abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU30", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","TrackerMatched_PU30",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU30", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                           }

                           // to get a nice view of the tracker in the r vs z plane
                           //No PileUp selection
                           //Selected
                        if ( abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","Selec_PU30", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU30", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                           }
                           //
                        if ( minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU30", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU30", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                           }
                     }

                  if (minitree_tree_nPV->at(0) < 35)
                     {
                        nSecIntFullSelec_PU35++;
                        if (minitree_SecInt_layer->at(iSecInt)!=0)
                           {
                              nSecIntTrackerMatched_PU35++;
                           }

                        if ( abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4&&abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","Selec_PU35", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","Selec_PU35",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU35", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);

                           }
                           //TrackerMatched
                        if(minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4 && abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU35", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","TrackerMatched_PU35",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU35", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);

                           }

                           // to get a nice view of the tracker in the r vs z plane
                           //No PileUp selection
                           //Selected
                        if ( abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","Selec_PU35", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU35", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                           }
                           //
                        if ( minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU35", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU35", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                           }
                     }

                  if (minitree_tree_nPV->at(0) < 40)
                     {
                        nSecIntFullSelec_PU40++;
                        if (minitree_SecInt_layer->at(iSecInt)!=0)
                           {
                              nSecIntTrackerMatched_PU40++;
                           }

                        if ( abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4&&abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","Selec_PU40", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","Selec_PU40",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU40", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                           }
                           //TrackerMatched
                        if(minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4 && abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU40", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","TrackerMatched_PU40",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU40", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                           }

                           // to get a nice view of the tracker in the r vs z plane
                           //No PileUp selection
                           //Selected
                        if ( abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","Selec_PU40", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU40", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                           }
                           //
                        if ( minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU40", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU40", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                           }
                     }

                  if (minitree_tree_nPV->at(0) < 45)
                     {
                        nSecIntFullSelec_PU45++;
                        if (minitree_SecInt_layer->at(iSecInt)!=0)
                           {
                              nSecIntTrackerMatched_PU45++;
                           }

                        if ( abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4&&abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","Selec_PU45", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","Selec_PU45",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU45", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);

                           }
                           //TrackerMatched
                        if(minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4 && abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU45", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","TrackerMatched_PU45",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU45", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                           }

                           // to get a nice view of the tracker in the r vs z plane
                           //No PileUp selection
                           //Selected
                        if ( abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","Selec_PU45", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU45", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                           }
                           //
                        if ( minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU45", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU45", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                           }
                     }
                  
                  if (minitree_tree_nPV->at(0) < 50)
                     {
                        nSecIntFullSelec_PU50++;
                        if (minitree_SecInt_layer->at(iSecInt)!=0)
                           {
                              nSecIntTrackerMatched_PU50++;
                           }

                        if ( abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4&&abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","Selec_PU50", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","Selec_PU50",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU50", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);

                           }
                           //TrackerMatched
                        if(minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4 && abs(minitree_SecInt_z->at(iSecInt))<27)
                           {
                              fillHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU50", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);
                              fillHisto("hData_reco_SecInt_mass","TrackerMatched_PU50",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU50", sample, minitree_SecInt_x->at(iSecInt),minitree_SecInt_y->at(iSecInt),1.);

                           }

                           // to get a nice view of the tracker in the r vs z plane
                           //No PileUp selection
                           //Selected
                        if ( abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","Selec_PU50", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU50", sample,abs(minitree_SecInt_z->at(iSecInt)), minitree_SecInt_r->at(iSecInt),1.);
                           }
                           //
                        if ( minitree_SecInt_layer->at(iSecInt)!=0 && abs(minitree_SecInt_r->at(iSecInt))<70 && abs(minitree_SecInt_z->at(iSecInt))<120)
                           {
                              fillHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU50", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                              fillHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU50", sample, abs(minitree_SecInt_z->at(iSecInt)),minitree_SecInt_r->at(iSecInt),1.);
                           }
                     }

                        //*******************************
                        //loop on Sec. Interactions tracks
                        //*******************************
                     float ddxy = 0.;
                     float ddz = 0.;
                     float firsthit_r = 0.;
                     if (minitree_SecInt_r->at(iSecInt) < 2.3 && minitree_SecInt_r->at(iSecInt)> 0.5)
                        {
                           for (unsigned int i = 0 ; i < minitree_V0_track_pt->size(); i++)
                              {
                                 if (minitree_V0_track_isFromSI->at(i) )
                                    {
                                       //fillhisto hitpattern
                                       fillHisto("hData_V0_track_hitpattern","Signal",  sample, minitree_V0_track_firstHit->at(i),1.);
                                       if (minitree_V0_track_firstHit->at(i) == 1160)
                                          {

                                                   fillHisto("hData_V0_track_pt","PF_Signal",  sample, minitree_V0_track_pt->at(i),1.);
                                                   fillHisto("hData_V0_track_eta","PF_Signal",  sample, minitree_V0_track_eta->at(i),1.);
                                                   fillHisto("hData_V0_track_phi","PF_Signal",  sample, minitree_V0_track_phi->at(i),1.);
                                                   fillHisto("hData_V0_track_charge","PF_Signal",  sample, minitree_V0_track_charge->at(i),1.);
                                                   fillHisto("hData_V0_track_NChi2","PF_Signal",  sample, minitree_V0_track_NChi2->at(i),1.);
                                                   fillHisto("hData_V0_track_dxy","PF_Signal",  sample, minitree_V0_track_dxy->at(i),1.);
                                                   if (minitree_V0_track_dxy->at(i) > 0) {ddxy = minitree_V0_track_drSig->at(i)/minitree_V0_track_dxy->at(i);}
                                                   fillHisto2D("hData_V0_track_ddxy","PF_Signal",  sample,minitree_V0_track_dxy->at(i),minitree_V0_track_drSig->at(i),1.);
                                                   fillHisto("hData_V0_track_dz","PF_Signal",  sample, minitree_V0_track_dz->at(i),1.);
                                                   if (minitree_V0_track_dz->at(i) > 0) {ddz = minitree_V0_track_dzSig->at(i)/minitree_V0_track_dz->at(i);}
                                                   fillHisto2D("hData_V0_track_ddz","PF_Signal",  sample, minitree_V0_track_dz->at(i),minitree_V0_track_dzSig->at(i),1.);

                                                   fillHisto("hData_V0_track_nHit","PF_Signal",  sample, minitree_V0_track_nHit->at(i),1.);
                                                   fillHisto("hData_V0_track_nHitPixel","PF_Signal",  sample, minitree_V0_track_nHitPixel->at(i),1.);
                                                   firsthit_r = sqrt(minitree_V0_track_firstHit_x->at(i)*minitree_V0_track_firstHit_x->at(i) + minitree_V0_track_firstHit_y->at(i)*minitree_V0_track_firstHit_y->at(i));
                                                   fillHisto("hData_V0_track_firstHit_x","PF_Signal",  sample, minitree_V0_track_firstHit_x->at(i),1.);
                                                   fillHisto("hData_V0_track_firstHit_y","PF_Signal",  sample, minitree_V0_track_firstHit_y->at(i),1.);
                                                   fillHisto2D("hData_V0_track_firstHit_xy","PF_Signal",  sample, minitree_V0_track_firstHit_x->at(i), minitree_V0_track_firstHit_y->at(i),1);
                                                   fillHisto("hData_V0_track_firstHit_r","PF_Signal",  sample, firsthit_r,1.);
                                                   fillHisto("hData_V0_track_firstHit_z","PF_Signal",  sample, minitree_V0_track_firstHit_z->at(i),1.);
                                                   fillHisto("hData_V0_track_iJet","PF_Signal",  sample, minitree_V0_track_iJet->at(i),1.);
                                                

                                          }
                                    }
                              }
                        }
                     else if (minitree_SecInt_r->at(iSecInt) < 3.35 && minitree_SecInt_r->at(iSecInt)> 2.3)
                        {
                           for (unsigned int i = 0 ; i < minitree_V0_track_pt->size(); i++)
                              {
                                 if (minitree_V0_track_isFromSI->at(i) )
                                    {
                                       fillHisto("hData_V0_track_hitpattern","BKG",  sample, minitree_V0_track_firstHit->at(i),1.);
                                       if (minitree_V0_track_firstHit->at(i) == 1160)
                                          {

                                                   fillHisto("hData_V0_track_pt","PF_BKG",  sample, minitree_V0_track_pt->at(i),1.);
                                                   fillHisto("hData_V0_track_eta","PF_BKG",  sample, minitree_V0_track_eta->at(i),1.);
                                                   fillHisto("hData_V0_track_phi","PF_BKG",  sample, minitree_V0_track_phi->at(i),1.);
                                                   fillHisto("hData_V0_track_charge","PF_BKG",  sample, minitree_V0_track_charge->at(i),1.);
                                                   fillHisto("hData_V0_track_NChi2","PF_BKG",  sample, minitree_V0_track_NChi2->at(i),1.);
                                                   fillHisto("hData_V0_track_dxy","PF_BKG",  sample, minitree_V0_track_dxy->at(i),1.);
                                                   if (minitree_V0_track_dxy->at(i) > 0) {ddxy = minitree_V0_track_drSig->at(i)/minitree_V0_track_dxy->at(i);}
                                                   fillHisto2D("hData_V0_track_ddxy","PF_BKG",  sample,minitree_V0_track_dxy->at(i),minitree_V0_track_drSig->at(i),1.);
                                                   fillHisto("hData_V0_track_dz","PF_BKG",  sample, minitree_V0_track_dz->at(i),1.);
                                                   if (minitree_V0_track_dz->at(i) > 0) {ddz = minitree_V0_track_dzSig->at(i)/minitree_V0_track_dz->at(i);}
                                                   fillHisto2D("hData_V0_track_ddz","PF_BKG",  sample, minitree_V0_track_dz->at(i),minitree_V0_track_dzSig->at(i),1.);

                                                   fillHisto("hData_V0_track_nHit","PF_BKG",  sample, minitree_V0_track_nHit->at(i),1.);
                                                   fillHisto("hData_V0_track_nHitPixel","PF_BKG",  sample, minitree_V0_track_nHitPixel->at(i),1.);
                                                   firsthit_r = sqrt(minitree_V0_track_firstHit_x->at(i)*minitree_V0_track_firstHit_x->at(i) + minitree_V0_track_firstHit_y->at(i)*minitree_V0_track_firstHit_y->at(i));
                                                   fillHisto("hData_V0_track_firstHit_x","PF_BKG",  sample, minitree_V0_track_firstHit_x->at(i),1.);
                                                   fillHisto("hData_V0_track_firstHit_y","PF_BKG",  sample, minitree_V0_track_firstHit_y->at(i),1.);
                                                   fillHisto2D("hData_V0_track_firstHit_xy","PF_BKG",  sample, minitree_V0_track_firstHit_x->at(i), minitree_V0_track_firstHit_y->at(i),1);
                                                   fillHisto("hData_V0_track_firstHit_r","PF_BKG",  sample, firsthit_r,1.);
                                                   fillHisto("hData_V0_track_firstHit_z","PF_BKG",  sample, minitree_V0_track_firstHit_z->at(i),1.);
                                                   fillHisto("hData_V0_track_iJet","PF_BKG",  sample, minitree_V0_track_iJet->at(i),1.);
                                                
                                          }
                                    }

                              }
                        }
                     else if (minitree_SecInt_x->at(iSecInt) < 0.25 && minitree_SecInt_x->at(iSecInt)> -0.04
                              && minitree_SecInt_y->at(iSecInt) < 5.0 && minitree_SecInt_y->at(iSecInt)> -5.0)
                        {
                           for (unsigned int i = 0 ; i < minitree_V0_track_pt->size(); i++)
                              {
                                 if (minitree_V0_track_isFromSI->at(i) )
                                    {
                                       fillHisto("hData_V0_track_hitpattern","PureSignal",  sample, minitree_V0_track_firstHit->at(i),1.);
                                       if (minitree_V0_track_firstHit->at(i) == 1160)
                                          {

                                                   fillHisto("hData_V0_track_pt","PF_PureSignal",  sample, minitree_V0_track_pt->at(i),1.);
                                                   fillHisto("hData_V0_track_eta","PF_PureSignal",  sample, minitree_V0_track_eta->at(i),1.);
                                                   fillHisto("hData_V0_track_phi","PF_PureSignal",  sample, minitree_V0_track_phi->at(i),1.);
                                                   fillHisto("hData_V0_track_charge","PF_PureSignal",  sample, minitree_V0_track_charge->at(i),1.);
                                                   fillHisto("hData_V0_track_NChi2","PF_PureSignal",  sample, minitree_V0_track_NChi2->at(i),1.);
                                                   fillHisto("hData_V0_track_dxy","PF_PureSignal",  sample, minitree_V0_track_dxy->at(i),1.);
                                                   if (minitree_V0_track_dxy->at(i) > 0) {ddxy = minitree_V0_track_drSig->at(i)/ minitree_V0_track_dxy->at(i);}
                                                   fillHisto2D("hData_V0_track_ddxy","PF_PureSignal",  sample, minitree_V0_track_dxy->at(i),minitree_V0_track_drSig->at(i),1.);
                                                   fillHisto("hData_V0_track_dz","PF_PureSignal",  sample, minitree_V0_track_dz->at(i),1.);
                                                   if (minitree_V0_track_dz->at(i) > 0) {ddz = minitree_V0_track_dzSig->at(i)/ minitree_V0_track_dz->at(i);}
                                                   fillHisto2D("hData_V0_track_ddz","PF_PureSignal",  sample, minitree_V0_track_dz->at(i),minitree_V0_track_dzSig->at(i),1.);

                                                   fillHisto("hData_V0_track_nHit","PF_PureSignal",  sample, minitree_V0_track_nHit->at(i),1.);
                                                   fillHisto("hData_V0_track_nHitPixel","PF_PureSignal",  sample, minitree_V0_track_nHitPixel->at(i),1.);
                                                   firsthit_r = sqrt(minitree_V0_track_firstHit_x->at(i)*minitree_V0_track_firstHit_x->at(i) + minitree_V0_track_firstHit_y->at(i)*minitree_V0_track_firstHit_y->at(i));
                                                   fillHisto("hData_V0_track_firstHit_x","PF_PureSignal",  sample, minitree_V0_track_firstHit_x->at(i),1.);
                                                   fillHisto("hData_V0_track_firstHit_y","PF_PureSignal",  sample, minitree_V0_track_firstHit_y->at(i),1);
                                                   fillHisto2D("hData_V0_track_firstHit_xy","PF_PureSignal",  sample, minitree_V0_track_firstHit_x->at(i), minitree_V0_track_firstHit_y->at(i),1);
                                                   fillHisto("hData_V0_track_firstHit_r","PF_PureSignal",  sample, firsthit_r,1.);
                                                   fillHisto("hData_V0_track_firstHit_z","PF_PureSignal",  sample, minitree_V0_track_firstHit_z->at(i),1.);
                                                   fillHisto("hData_V0_track_iJet","PF_PureSignal",  sample, minitree_V0_track_iJet->at(i),1.);
                                                
                                          }
                                    }

                              }
                        }


 
               }
         }//end loop on SecInt

      fillHisto("hData_reco_nSecInt","FullSelec",  sample, nSecIntFullSelec,1.);
      fillHisto("hData_reco_nSecInt","TrackerMatched",  sample, nSecIntFullSelec,1.);

      fillHisto("hData_reco_nSecInt","FullSelec_PU25",  sample, nSecIntFullSelec_PU25,1.);
      fillHisto("hData_reco_nSecInt","TrackerMatched_PU25",  sample, nSecIntFullSelec_PU25,1.);

      fillHisto("hData_reco_nSecInt","FullSelec_PU30",  sample, nSecIntFullSelec_PU30,1.);
      fillHisto("hData_reco_nSecInt","TrackerMatched_PU30",  sample, nSecIntFullSelec_PU30,1.);

      fillHisto("hData_reco_nSecInt","FullSelec_PU35",  sample, nSecIntFullSelec_PU35,1.);
      fillHisto("hData_reco_nSecInt","TrackerMatched_PU35",  sample, nSecIntFullSelec_PU35,1.);

      fillHisto("hData_reco_nSecInt","FullSelec_PU40",  sample, nSecIntFullSelec_PU40,1.);
      fillHisto("hData_reco_nSecInt","TrackerMatched_PU40",  sample, nSecIntFullSelec_PU40,1.);

      fillHisto("hData_reco_nSecInt","FullSelec_PU45",  sample, nSecIntFullSelec_PU45,1.);
      fillHisto("hData_reco_nSecInt","TrackerMatched_PU45",  sample, nSecIntFullSelec_PU45,1.);

      fillHisto("hData_reco_nSecInt","FullSelec_PU50",  sample, nSecIntFullSelec_PU50,1.);
      fillHisto("hData_reco_nSecInt","TrackerMatched_PU50",  sample, nSecIntFullSelec_PU50,1.);


   }

float NENTRIES = nentries;
float ratio = TOTALnSecIntTrackerMatched/NENTRIES;
fillHisto("hData_reco_nAveSecInt","TrackerMatched",  sample, ratio,1.);
   // GetHistDirectory();
   theoutputfile->Write();
   std::cout<<"end of loop : "<<theoutputfile->Write()<<std::endl;

   theoutputfile->Close();
   delete theoutputfile;
}

void TreeSecIntReader::initializeHisto(TString sample, bool isfirstset){


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
    
   addHisto("hData_reco_K0_mass","", sample.Data(), 140,0.43,0.57);
   addHisto("hData_reco_L0_mass","",  sample.Data(), 99,1.066,1.165);

   addHisto2D("hData_reco_V0_xy","", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_reco_V0_rz","", sample.Data(), 1200,0.,120.,700,0.,70.);
   addHisto2D("hData_CMSSW_V0_xy","", sample.Data(),500,-5,5,500,-5,5);
   addHisto2D("hData_CMSSW_V0_rz","", sample.Data(), 1200,0.,120.,700,0.,70.);
   addHisto2D("hData_CMSSW_Yc_xy","", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_CMSSW_Yc_rz","", sample.Data(), 1200,0.,120.,700,0.,70.);



   addHisto("hData_reco_SecInt_mass","FullSelec",  sample.Data(), 100,0,10);
   addHisto("hData_reco_SecInt_mass","FullTrackerMatched",  sample.Data(), 100,0,10);

   addHisto2D("hData_reco_SecInt_xy","Selec", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_xy","TrackerMatched", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_rz","Selec", sample.Data(),1200,0.,120.,700,0.,70.);
   addHisto2D("hData_reco_SecInt_rz","TrackerMatched", sample.Data(),1200,0.,120.,700,0.,70.);

   addHisto2D("hData_reco_SecInt_xy_Inner","Selec", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched", sample.Data(),500,-5,5,500,-5,5);
      addHisto2D("hData_reco_SecInt_rz_Inner","Selec", sample.Data(), 300,0,30,180,0,12);
   addHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched", sample.Data(),300,0,30,180,0,12);

   addHisto2D("hData_reco_SecInt_Sigxy_xy","Selec", sample.Data(), 36000,0,12000,250,0,5);
   addHisto2D("hData_reco_SecInt_xy","Selec_Spe1000", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_reco_SecInt_xy","Selec_Spe9000", sample.Data(), 500,-5,5,500,-5,5);

   addHisto("hData_reco_SecInt_mass","Selec",  sample.Data(), 100,0,10);
   addHisto("hData_reco_SecInt_mass","TrackerMatched",  sample.Data(), 100,0,10);

      addHisto("hData_reco_SecInt_mass","Selec_PU25",  sample.Data(), 100,0,10);
   addHisto("hData_reco_SecInt_mass","TrackerMatched_PU25",  sample.Data(), 100,0,10);
   addHisto2D("hData_reco_SecInt_xy","Selec_PU25", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU25", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_rz","Selec_PU25", sample.Data(),1200,0.,120.,700,0.,70.);
   addHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU25", sample.Data(), 1200,0.,120.,700,0.,70.);

      addHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU25", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU25", sample.Data(), 500,-5,5,500,-5,5);
         addHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU25", sample.Data(), 300,0,30,180,0,12);
   addHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU25", sample.Data(), 300,0,30,180,0,12);

      addHisto("hData_reco_SecInt_mass","Selec_PU30",  sample.Data(), 100,0,10);
   addHisto("hData_reco_SecInt_mass","TrackerMatched_PU30",  sample.Data(), 100,0,10);
   addHisto2D("hData_reco_SecInt_xy","Selec_PU30", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU30", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_rz","Selec_PU30", sample.Data(),1200,0.,120.,700,0.,70.);
   addHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU30", sample.Data(), 1200,0.,120.,700,0.,70.);

         addHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU30", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU30", sample.Data(), 500,-5,5,500,-5,5);
            addHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU30", sample.Data(), 300,0,30,180,0,12);
   addHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU30", sample.Data(), 300,0,30,180,0,12);


      addHisto("hData_reco_SecInt_mass","Selec_PU35",  sample.Data(), 100,0,10);
   addHisto("hData_reco_SecInt_mass","TrackerMatched_PU35",  sample.Data(), 100,0,10);
   addHisto2D("hData_reco_SecInt_xy","Selec_PU35", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU35", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_rz","Selec_PU35", sample.Data(),1200,0.,120.,700,0.,70.);
   addHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU35", sample.Data(), 1200,0.,120.,700,0.,70.);
      addHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU35", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU35", sample.Data(), 500,-5,5,500,-5,5);
         addHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU35", sample.Data(), 300,0,30,180,0,12);
   addHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU35", sample.Data(), 300,0,30,180,0,12);

   addHisto("hData_reco_SecInt_mass","Selec_PU40",  sample.Data(), 100,0,10);
   addHisto("hData_reco_SecInt_mass","TrackerMatched_PU40",  sample.Data(), 100,0,10);
   addHisto2D("hData_reco_SecInt_xy","Selec_PU40", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU40", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_rz","Selec_PU40", sample.Data(),1200,0.,120.,700,0.,70.);
   addHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU40", sample.Data(), 1200,0.,120.,700,0.,70.);
         addHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU40", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU40", sample.Data(), 500,-5,5,500,-5,5);
            addHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU40", sample.Data(), 300,0,30,180,0,12);
   addHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU40", sample.Data(), 300,0,30,180,0,12);


      addHisto("hData_reco_SecInt_mass","Selec_PU45",  sample.Data(), 100,0,10);
   addHisto("hData_reco_SecInt_mass","TrackerMatched_PU45",  sample.Data(), 100,0,10);
   addHisto2D("hData_reco_SecInt_xy","Selec_PU45", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU45", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_rz","Selec_PU45", sample.Data(),1200,0.,120.,700,0.,70.);
   addHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU45", sample.Data(), 1200,0.,120.,700,0.,70.);
         addHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU45", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU45", sample.Data(), 500,-5,5,500,-5,5);
            addHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU45", sample.Data(), 300,0,30,180,0,12);
   addHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU45", sample.Data(), 300,0,30,180,0,12);


   addHisto("hData_reco_SecInt_mass","Selec_PU50",  sample.Data(), 100,0,10);
   addHisto("hData_reco_SecInt_mass","TrackerMatched_PU50",  sample.Data(), 100,0,10);
   addHisto2D("hData_reco_SecInt_xy","Selec_PU50", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_xy","TrackerMatched_PU50", sample.Data(), 500,-25.,25.,500,-25.,25.);
   addHisto2D("hData_reco_SecInt_rz","Selec_PU50", sample.Data(),1200,0.,120.,700,0.,70.);
   addHisto2D("hData_reco_SecInt_rz","TrackerMatched_PU50", sample.Data(), 1200,0.,120.,700,0.,70.);
   addHisto2D("hData_reco_SecInt_xy_Inner","Selec_PU50", sample.Data(), 500,-5,5,500,-5,5);
   addHisto2D("hData_reco_SecInt_xy_Inner","TrackerMatched_PU50", sample.Data(), 500,-5,5,500,-5,5);
      addHisto2D("hData_reco_SecInt_rz_Inner","Selec_PU50", sample.Data(), 300,0,30,180,0,12);
   addHisto2D("hData_reco_SecInt_rz_Inner","TrackerMatched_PU50", sample.Data(), 300,0,30,180,0,12);


      //----------- tracks from SI  ------//

   addHisto("hData_V0_track_hitpattern","Signal",  sample.Data(),140 ,1150,1710);
   addHisto("hData_V0_track_hitpattern","BKG",  sample.Data(), 140 ,1150,1710);
   addHisto("hData_V0_track_hitpattern","PureSignal",  sample.Data(), 140 ,1150,1710);



   addHisto("hData_V0_track_pt","PF_Signal",  sample.Data(), 200,0,20);
   addHisto("hData_V0_track_eta","PF_Signal",  sample.Data(), 20,0,2);
   addHisto("hData_V0_track_phi","PF_Signal",  sample.Data(), 618,-3.14,3.14);
   addHisto("hData_V0_track_charge","PF_Signal",  sample.Data(), 4,-2,2);
   addHisto("hData_V0_track_NChi2","PF_Signal",  sample.Data(), 100,0,10);
   addHisto("hData_V0_track_dxy","PF_Signal",  sample.Data(), 100,-5,5);
   addHisto2D("hData_V0_track_ddxy","PF_Signal",  sample.Data(), 100,-5,5, 200 ,-100,100);
   addHisto("hData_V0_track_dz","PF_Signal",  sample.Data(),200 ,-10,10);
   addHisto2D("hData_V0_track_ddz","PF_Signal",  sample.Data(),200,-10,10, 200 ,-100,100);
   addHisto("hData_V0_track_nHit","PF_Signal",  sample.Data(), 30,0,30);
   addHisto("hData_V0_track_nHitPixel","PF_Signal",  sample.Data(),15,0,15);
   addHisto("hData_V0_track_firstHit_x","PF_Signal",  sample.Data(),  100,-5,5);
   addHisto("hData_V0_track_firstHit_y","PF_Signal",  sample.Data(),  100,-5,5);
   addHisto2D("hData_V0_track_firstHit_xy","PF_Signal",  sample.Data(), 100,-5,5,100,-5,5);
   addHisto("hData_V0_track_firstHit_r","PF_Signal",  sample.Data(),  50,0,5);
   addHisto("hData_V0_track_firstHit_z","PF_Signal",  sample.Data(), 400,-20,20);
   addHisto("hData_V0_track_iJet","PF_Signal",  sample.Data(),  30,0,30);


   addHisto("hData_V0_track_pt","PF_BKG",  sample.Data(), 200,0,20);
   addHisto("hData_V0_track_eta","PF_BKG",  sample.Data(), 20,0,2);
   addHisto("hData_V0_track_phi","PF_BKG",  sample.Data(), 618,-3.14,3.14);
   addHisto("hData_V0_track_charge","PF_BKG",  sample.Data(), 4,-2,2);
   addHisto("hData_V0_track_NChi2","PF_BKG",  sample.Data(), 100,0,10);
   addHisto("hData_V0_track_dxy","PF_BKG",  sample.Data(), 100,-5,5);
   addHisto2D("hData_V0_track_ddxy","PF_BKG",  sample.Data(), 100,-5,5, 200 ,-100,100);
   addHisto("hData_V0_track_dz","PF_BKG",  sample.Data(), 200 ,-10,10);
   addHisto2D("hData_V0_track_ddz","PF_BKG",  sample.Data(), 200,-10,10, 200 ,-100,100);
   addHisto("hData_V0_track_nHit","PF_BKG",  sample.Data(), 30,0,30);
   addHisto("hData_V0_track_nHitPixel","PF_BKG",  sample.Data(), 15,0,15);
   addHisto("hData_V0_track_firstHit_x","PF_BKG",  sample.Data(),  100,-5,5);
   addHisto("hData_V0_track_firstHit_y","PF_BKG",  sample.Data(),  100,-5,5);
   addHisto2D("hData_V0_track_firstHit_xy","PF_BKG",  sample.Data(), 100,-5,5,100,-5,5);
   addHisto("hData_V0_track_firstHit_r","PF_BKG",  sample.Data(),  50,0,5);
   addHisto("hData_V0_track_firstHit_z","PF_BKG",  sample.Data(), 400,-20,20);
   addHisto("hData_V0_track_iJet","PF_BKG",  sample.Data(),  30,0,30);


   addHisto("hData_V0_track_pt","PF_PureSignal",  sample.Data(), 200,0,20);
   addHisto("hData_V0_track_eta","PF_PureSignal",  sample.Data(), 20,0,2);
   addHisto("hData_V0_track_phi","PF_PureSignal",  sample.Data(), 618,-3.14,3.14);
   addHisto("hData_V0_track_charge","PF_PureSignal",  sample.Data(),4,-2,2);
   addHisto("hData_V0_track_NChi2","PF_PureSignal",  sample.Data(), 100,0,10);
   addHisto("hData_V0_track_dxy","PF_PureSignal",  sample.Data(), 100,-5,5);
   addHisto2D("hData_V0_track_ddxy","PF_PureSignal",  sample.Data(),100,-5,5, 200 ,-100,100);
   addHisto("hData_V0_track_dz","PF_PureSignal",  sample.Data(), 200 ,-10,10);
   addHisto2D("hData_V0_track_ddz","PF_PureSignal",  sample.Data(),200,-10,10, 200 ,-100,100);
   addHisto("hData_V0_track_nHit","PF_PureSignal",  sample.Data(), 30,0,30);
   addHisto("hData_V0_track_nHitPixel","PF_PureSignal",  sample.Data(), 15,0,15);
   addHisto("hData_V0_track_firstHit_x","PF_PureSignal",  sample.Data(),  100,-5,5);
   addHisto("hData_V0_track_firstHit_y","PF_PureSignal",  sample.Data(),  100,-5,5);
   addHisto2D("hData_V0_track_firstHit_xy","PF_PureSignal",  sample.Data(), 100,-5,5,100,-5,5);
   addHisto("hData_V0_track_firstHit_r","PF_PureSignal",  sample.Data(),  50,0,5);
   addHisto("hData_V0_track_firstHit_z","PF_PureSignal",  sample.Data(), 400,-20,20);
   addHisto("hData_V0_track_iJet","PF_PureSignal",  sample.Data(),  30,0,30);


   //-----------//
   addHisto("hData_reco_nSecInt","FullSelec",  sample.Data(), 20,0,20);
   addHisto("hData_reco_nSecInt","TrackerMatched",  sample.Data(),  20,0,20);
   addHisto("hData_reco_nAveSecInt","TrackerMatched", sample.Data(),100,0,10);

   addHisto("hData_reco_nSecInt","FullSelec_PU25",  sample.Data(),  20,0,20);
   addHisto("hData_reco_nSecInt","TrackerMatched_PU25",  sample.Data(),  20,0,20);

   addHisto("hData_reco_nSecInt","FullSelec_PU30",  sample.Data(),  20,0,20);
   addHisto("hData_reco_nSecInt","TrackerMatched_PU30",  sample.Data(),  20,0,20);

   addHisto("hData_reco_nSecInt","FullSelec_PU35",  sample.Data(), 20,0,20);
   addHisto("hData_reco_nSecInt","TrackerMatched_PU35",  sample.Data(),  20,0,20);

   addHisto("hData_reco_nSecInt","FullSelec_PU40",  sample.Data(),  20,0,20);
   addHisto("hData_reco_nSecInt","TrackerMatched_PU40",  sample.Data(),  20,0,20);

   addHisto("hData_reco_nSecInt","FullSelec_PU45",  sample.Data(),  20,0,20);
   addHisto("hData_reco_nSecInt","TrackerMatched_PU45",  sample.Data(), 20,0,20);

   addHisto("hData_reco_nSecInt","FullSelec_PU50",  sample.Data(),  20,0,20);
   addHisto("hData_reco_nSecInt","TrackerMatched_PU50",  sample.Data(),  20,0,20);

   addHisto2D("hData_reco_3DSecInt","", sample.Data(), 15,180,480,100,0,100);
}

//-------------------------------------------------------------
//instantiate and add
//first parameter is the variable name,
// second parameter is the selection step (like "afterleptsel")
//third parameter is the sample name (like "Z)
//others are TH1F binning
//creates one histograms per channel
//-------------------------------------------------------------
void TreeSecIntReader::addHisto(TString var, TString selstep, TString sample, int nbins, float min, float max){
 
   TString name =  sample+"_"+var+"_"+selstep;
  TH1F * thehisto = new TH1F(name,name,nbins,min,max);
  thehisto->Sumw2();
  thehisto->SetOption("HIST");
   // std::cout<<"adding histo with name : "<<name<<std::endl;
  histo_list_.push_back(thehisto);
  histo_map_[name.Data()] = numb_histo;
  numb_histo++;
}

void TreeSecIntReader::addHisto2D(TString var, TString selstep, TString sample, int nxbins, float xmin, float xmax, int nybins, float ymin, float ymax){
 
  TString name =  sample+"_"+var+"_"+selstep;
  TH2F * thehisto = new TH2F(name,name,nxbins,xmin,xmax,nybins,ymin,ymax);
  // thehisto->Sumw2();
  thehisto->SetOption("COL");
   // std::cout<<"adding histo with name : "<<name<<std::endl;
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
void TreeSecIntReader::fillHisto( TString var, TString selstep,TString sample, float val, float weight){
  TString name = sample+"_"+var+"_"+selstep;


  if(histo_map_[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histograms " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }else  {histo_list_[histo_map_[name.Data()]]->Fill(val, weight);}
  
}


void TreeSecIntReader::fillHisto2D( TString var, TString selstep, TString sample, float xval,float yval, float weight){
  TString name = sample+"_"+var+"_"+selstep;


  if(histo_map_2D_[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histograms " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }else  histo_list_2D_[histo_map_2D_[name.Data()]]->Fill(xval,yval, weight);
  
}


void TreeSecIntReader::GetHistDirectory(){
  for(unsigned int i=0; i<histo_list_.size(); i++){
    std::cout << "Histogram directory: " << histo_list_[i]->GetDirectory()->GetName() << std::endl;
  }
  for(unsigned int i=0; i<histo_list_2D_.size(); i++){
    std::cout << "Histogram directory: " << histo_list_2D_[i]->GetDirectory()->GetName() << std::endl;
  }
}


void TreeSecIntReader::deleteHisto(){
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
