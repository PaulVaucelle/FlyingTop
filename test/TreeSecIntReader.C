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

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {//nentries
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if ( jentry%10000 == 0 ) std::cout << "events : " << jentry << std::endl;

      for (unsigned int i = 0; i < minitree_K0_reco_mass->size(); i++)
         {
            fillHisto("hData_reco_K0_mass","",  sample, minitree_K0_reco_mass->at(i),1.);
         }
      for (unsigned int i = 0; i < minitree_L0_reco_mass->size(); i++)
         {
            fillHisto("hData_reco_L0_mass","",  sample, minitree_L0_reco_mass->at(i),1.);
         }

      //*******************************
      //loop on Sec. Interactions
      //*******************************

      for (unsigned int iSecInt = 0; iSecInt <minitree_SecInt_mass->size(); iSecInt++)
         {

            if (minitree_SecInt_selec->at(iSecInt))
               {
                  fillHisto("hData_reco_SecInt_mass","FullSelec",  sample, minitree_SecInt_mass->at(iSecInt),1.);

                  if(minitree_SecInt_layer->at(iSecInt)) 
                     {
                        fillHisto("hData_reco_SecInt_mass","FullTrackerMatched",  sample, minitree_SecInt_mass->at(iSecInt),1.);
                     }

                  //to get a nice view of the inner tracker x vs y
                  //Selected
                  if ( abs(minitree_SecInt_x->at(iSecInt))<25 && abs(minitree_SecInt_y->at(iSecInt))<25 && abs(minitree_SecInt_eta->at(iSecInt))<1.4&&abs(minitree_SecInt_z->at(iSecInt))<27)
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
 
               }
         }//end loop on SecInt

   }


   // GetHistDirectory();
   theoutputfile->Write();
   std::cout<<"end of loop : "<<theoutputfile->Write()<<std::endl;//files are not written somehow

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
