#define MiniSecIntNtuple_cxx
#include "MiniSecIntNtuple.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "HistogramManager.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "TTree.h"
using namespace std;



void MiniSecIntNtuple::Loop(TString sample , TString Production,bool Signal )
{

std::cout<<"//-------------------------//"<<std::endl;
//   In a ROOT session, you can do:
//      root> .L MiniSecIntNtuple.C
//      root> MiniSecIntNtuple t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

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
   TFile * myFile = new TFile( (Production+"/MiniSecInt_"+sample+".root").Data(), "recreate");
   TTree *smalltree = new TTree("ttree", "summary information");


   std::vector<int>    minirunNumber;
   std::vector<int>    minieventNumber;
   std::vector<int>    minilumiBlock;

  std::vector<float>      minitree_K0_reco_mass;
  std::vector<float>      minitree_L0_reco_mass;
  std::vector<int>       minitree_tree_nPV;
  std::vector<int>       minitree_nSecInt;
  std::vector<float>     minitree_SecInt_x;
  std::vector<float>     minitree_SecInt_y;
  std::vector<float>     minitree_SecInt_z;
  std::vector<float>     minitree_SecInt_r;
  std::vector<float>     minitree_SecInt_d;
  std::vector<float>     minitree_SecInt_drSig;
  std::vector<float>     minitree_SecInt_dzSig;
  std::vector<float>     minitree_SecInt_angleXY;
  std::vector<float>     minitree_SecInt_angleZ;
  std::vector<float>     minitree_SecInt_NChi2;
  std::vector<float>     minitree_SecInt_ndf;
  std::vector<float>     minitree_SecInt_mass;
  std::vector<float>     minitree_SecInt_pt;
  std::vector<float>     minitree_SecInt_eta;
  std::vector<float>     minitree_SecInt_phi;
  std::vector<int>       minitree_SecInt_charge;
  std::vector<bool>      minitree_SecInt_badTkHit;
  std::vector<float>     minitree_SecInt_dca;
  std::vector<bool>      minitree_SecInt_selec;
  std::vector<int>       minitree_SecInt_layer;
     
   // CMSSW collection
  vector<float>   minitree_K0_x;
  vector<float>   minitree_K0_y;
  vector<float>   minitree_K0_z;
  vector<float>   minitree_K0_r;
  vector<float>   minitree_K0_NChi2;
  vector<float>   minitree_K0_mass;
  vector<float>   minitree_K0_eta;

  // CMSSW collection
  vector<float>   minitree_L0_x;
  vector<float>   minitree_L0_y;
  vector<float>   minitree_L0_z;
  vector<float>   minitree_L0_r;
  vector<float>   minitree_L0_NChi2;
  vector<float>   minitree_L0_mass;
  vector<float>   minitree_L0_eta;

  // reco from us 
  vector<float>   minitree_V0_reco_x;
  vector<float>   minitree_V0_reco_y;
  vector<float>   minitree_V0_reco_z;
  vector<float>   minitree_V0_reco_r;
  vector<float>   minitree_V0_reco_NChi2;
  vector<float>   minitree_V0_reco_mass;
  vector<float>   minitree_V0_reco_eta;
  vector<int>     minitree_V0_reco_source;

   //CMSSW collection
  vector<float>   minitree_Yc_x;
  vector<float>   minitree_Yc_y;
  vector<float>   minitree_Yc_z;
  vector<float>   minitree_Yc_r;
  vector<int>     minitree_Yc_layer;
  vector<float>   minitree_Yc_NChi2;
  vector<float>   minitree_Yc_eta;
  vector<float>   minitree_Yc_mass;

//secint tracks
  vector<bool>    minitree_V0_track_isFromV0;
  vector<bool>    minitree_V0_track_isFromSI;
  vector<bool>    minitree_V0_track_lost;
  vector<float>   minitree_V0_track_pt;
  vector<float>   minitree_V0_track_eta;
  vector<float>   minitree_V0_track_phi;
  vector<int>     minitree_V0_track_charge;
  vector<float>   minitree_V0_track_NChi2;
  vector<float>   minitree_V0_track_dxy;
  vector<float>   minitree_V0_track_drSig;
  vector<float>   minitree_V0_track_dz;
  vector<float>   minitree_V0_track_dzSig;
  vector<int>     minitree_V0_track_nHit;
  vector<int>     minitree_V0_track_nHitPixel;
  vector<int>     minitree_V0_track_firstHit;
  vector<float>   minitree_V0_track_firstHit_x;
  vector<float>   minitree_V0_track_firstHit_y;
  vector<float>   minitree_V0_track_firstHit_z;
  vector<int>     minitree_V0_track_iJet;
  vector<float>   minitree_V0_track_ntrk10;
  vector<float>   minitree_V0_track_ntrk20;
  vector<float>   minitree_V0_track_ntrk30;
  vector<float>   minitree_V0_track_ntrk40;
  vector<int>     minitree_V0_track_Hemi;
  vector<float>   minitree_V0_track_Hemi_dR;
  vector<float>   minitree_V0_track_Hemi_dRmax;


  vector<int>           minitree_smu_mass;
  vector<int>           minitree_neu_mass;
  vector<float>         minitree_neu_ctau;

   smalltree->Branch("minirunNumber",&minirunNumber);
   smalltree->Branch("minieventNumber",&minieventNumber);
   smalltree->Branch("minilumiBlock",&minilumiBlock);


   smalltree->Branch("minitree_K0_reco_mass",&minitree_K0_reco_mass);
   smalltree->Branch("minitree_L0_reco_mass",&minitree_L0_reco_mass);
   smalltree->Branch("minitree_tree_nPV",&minitree_tree_nPV);
   smalltree->Branch("minitree_nSecInt",&minitree_nSecInt);
   smalltree->Branch("minitree_SecInt_x",&minitree_SecInt_x);
   smalltree->Branch("minitree_SecInt_y",&minitree_SecInt_y);
   smalltree->Branch("minitree_SecInt_z",&minitree_SecInt_z);
   smalltree->Branch("minitree_SecInt_r",&minitree_SecInt_r);
   smalltree->Branch("minitree_SecInt_d",&minitree_SecInt_d);
   smalltree->Branch("minitree_SecInt_drSig",&minitree_SecInt_drSig);
   smalltree->Branch("minitree_SecInt_dzSig",&minitree_SecInt_dzSig);
   smalltree->Branch("minitree_SecInt_angleXY",&minitree_SecInt_angleXY);
   smalltree->Branch("minitree_SecInt_angleZ",&minitree_SecInt_angleZ);
   smalltree->Branch("minitree_SecInt_NChi2",&minitree_SecInt_NChi2);
   smalltree->Branch("minitree_SecInt_ndf",&minitree_SecInt_ndf);
   smalltree->Branch("minitree_SecInt_mass",&minitree_SecInt_mass);
   smalltree->Branch("minitree_SecInt_pt",&minitree_SecInt_pt);
   smalltree->Branch("minitree_SecInt_eta",&minitree_SecInt_eta);
   smalltree->Branch("minitree_SecInt_phi",&minitree_SecInt_phi);
   smalltree->Branch("minitree_SecInt_charge",&minitree_SecInt_charge);
   smalltree->Branch("minitree_SecInt_badTkHit",&minitree_SecInt_badTkHit);
   smalltree->Branch("minitree_SecInt_dca",&minitree_SecInt_dca);
   smalltree->Branch("minitree_SecInt_selec",&minitree_SecInt_selec);
   smalltree->Branch("minitree_SecInt_layer",&minitree_SecInt_layer);

   // CMSSW collection
  smalltree->Branch("minitree_K0_x",&minitree_K0_x);
  smalltree->Branch("minitree_K0_y",&minitree_K0_y);
  smalltree->Branch("minitree_K0_z",&minitree_K0_z);
  smalltree->Branch("minitree_K0_r",&minitree_K0_r);
  smalltree->Branch("minitree_K0_NChi2",&minitree_K0_NChi2);
  smalltree->Branch("minitree_K0_mass",&minitree_K0_mass);
  smalltree->Branch("minitree_K0_eta",&minitree_K0_eta);

  // CMSSW collection
  smalltree->Branch("minitree_L0_x",&minitree_L0_x);
  smalltree->Branch("minitree_L0_y",&minitree_L0_y);
  smalltree->Branch("minitree_L0_z",&minitree_L0_z);
  smalltree->Branch("minitree_L0_r",&minitree_L0_r);
  smalltree->Branch("minitree_L0_NChi2",&minitree_L0_NChi2);
  smalltree->Branch("minitree_L0_mass",&minitree_L0_mass);
  smalltree->Branch("minitree_L0_eta",&minitree_L0_eta);

  // reco from us 
  smalltree->Branch("minitree_V0_reco_x",&minitree_V0_reco_x);
  smalltree->Branch("minitree_V0_reco_y",&minitree_V0_reco_y);
  smalltree->Branch("minitree_V0_reco_z",&minitree_V0_reco_z);
  smalltree->Branch("minitree_V0_reco_r",&minitree_V0_reco_r);
  smalltree->Branch("minitree_V0_reco_NChi2",&minitree_V0_reco_NChi2);
  smalltree->Branch("minitree_V0_reco_mass",&minitree_V0_reco_mass);
  smalltree->Branch("minitree_V0_reco_eta",&minitree_V0_reco_eta);
  smalltree->Branch("minitree_V0_reco_source",&minitree_V0_reco_source);

   //CMSSW collection
  smalltree->Branch("minitree_Yc_x",&minitree_Yc_x);
  smalltree->Branch("minitree_Yc_y",&minitree_Yc_y);
  smalltree->Branch("minitree_Yc_z",&minitree_Yc_z);
  smalltree->Branch("minitree_Yc_r",&minitree_Yc_r);
  smalltree->Branch("minitree_Yc_layer",&minitree_Yc_layer);
  smalltree->Branch("minitree_Yc_NChi2",&minitree_Yc_NChi2);
  smalltree->Branch("minitree_Yc_eta",&minitree_Yc_eta);
  smalltree->Branch("minitree_Yc_mass",&minitree_Yc_mass);

//secint tracks
  smalltree->Branch("minitree_V0_track_isFromV0",&minitree_V0_track_isFromV0);
  smalltree->Branch("minitree_V0_track_isFromSI",&minitree_V0_track_isFromSI);
  smalltree->Branch("minitree_V0_track_lost",&minitree_V0_track_lost);
  smalltree->Branch("minitree_V0_track_pt",&minitree_V0_track_pt);
  smalltree->Branch("minitree_V0_track_eta",&minitree_V0_track_eta);
  smalltree->Branch("minitree_V0_track_phi",&minitree_V0_track_phi);
  smalltree->Branch("minitree_V0_track_charge",&minitree_V0_track_charge);
  smalltree->Branch("minitree_V0_track_NChi2",&minitree_V0_track_NChi2);
  smalltree->Branch("minitree_V0_track_dxy",&minitree_V0_track_dxy);
  smalltree->Branch("minitree_V0_track_drSig",&minitree_V0_track_drSig);
  smalltree->Branch("minitree_V0_track_dz",&minitree_V0_track_dz);
  smalltree->Branch("minitree_V0_track_dzSig",&minitree_V0_track_dzSig);
  smalltree->Branch("minitree_V0_track_nHit",&minitree_V0_track_nHit);
  smalltree->Branch("minitree_V0_track_nHitPixel",&minitree_V0_track_nHitPixel);
  smalltree->Branch("minitree_V0_track_firstHit",&minitree_V0_track_firstHit);
  smalltree->Branch("minitree_V0_track_firstHit_x",&minitree_V0_track_firstHit_x);
  smalltree->Branch("minitree_V0_track_firstHit_y",&minitree_V0_track_firstHit_y);
  smalltree->Branch("minitree_V0_track_firstHit_z",&minitree_V0_track_firstHit_z);
  smalltree->Branch("minitree_V0_track_iJet",&minitree_V0_track_iJet);
  smalltree->Branch("minitree_V0_track_ntrk10",&minitree_V0_track_ntrk10);
  smalltree->Branch("minitree_V0_track_ntrk20",&minitree_V0_track_ntrk20);
  smalltree->Branch("minitree_V0_track_ntrk30",&minitree_V0_track_ntrk30);
  smalltree->Branch("minitree_V0_track_ntrk40",&minitree_V0_track_ntrk40);
  smalltree->Branch("minitree_V0_track_Hemi",&minitree_V0_track_Hemi);
  smalltree->Branch("minitree_V0_track_Hemi_dR",&minitree_V0_track_Hemi_dR);
  smalltree->Branch("minitree_V0_track_Hemi_dRmax",&minitree_V0_track_Hemi_dRmax);

   smalltree->Branch("minitree_smu_mass",&minitree_smu_mass);
   smalltree->Branch("minitree_neu_mass",&minitree_neu_mass);
   smalltree->Branch("minitree_neu_ctau",&minitree_neu_ctau);

   std::cout<<"//-------------------------//"<<std::endl;
   std::cout<<"// "<<Production<<" : "<<sample<<"  //"<<std::endl;
   std::cout<<"//-------------------------//"<<std::endl;

   bool debug = false;
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "Total Entries : " << nentries << std::endl;

   Long64_t nbytes = 0, nb = 0;

   int allevents = 0;

   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      allevents++;
      // if ( allevents%10000 == 0 ) std::cout << "events : " << allevents << std::endl;

      // if (Cut(ientry) < 0) continue;
      if ( jentry%10000 == 0 ) std::cout << "events : " << jentry << std::endl;

    // -------------- --------------------------------------------------------------------//
    // -----------------------------------------------------------------------------------//
    if ( !((tree_Filter || tree_FilterSameSign) && tree_njetNOmu>0) && !Signal ) continue;
    //------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------//

      minirunNumber.push_back(runNumber);
      minieventNumber.push_back(eventNumber);
      minilumiBlock.push_back(lumiBlock);
      
      minitree_tree_nPV.push_back(tree_nPV);
      minitree_nSecInt.push_back(tree_nSecInt);

      int neuMASS = 0;
      int smuMASS = 0;
      float neuCTAU = 0;
      if (Signal)
        {
          neuMASS = tree_neu_mass;
          smuMASS = tree_smu_mass;
          neuCTAU = tree_neu_ctau;
          minitree_smu_mass.push_back(smuMASS);
          minitree_neu_mass.push_back(neuMASS);
          minitree_neu_ctau.push_back(neuCTAU);
        }



      for (int i=0; i<tree_nV0_reco; i++)    // Loop on V0s
      {
        if ( tree_V0_reco_source->at(i) == 1 ) {
          minitree_K0_reco_mass.push_back( tree_V0_reco_mass->at(i) );
        }
        if ( tree_V0_reco_source->at(i) == 2 ) {
          minitree_L0_reco_mass.push_back( tree_V0_reco_mass->at(i) );
        }
      }

  for (int i = 0 ; i < tree_K0_x->size() ; i++)
  {
    minitree_K0_x.push_back(tree_K0_x->at(i));
    minitree_K0_y.push_back(tree_K0_y->at(i));
    minitree_K0_z.push_back(tree_K0_z->at(i));
    minitree_K0_r.push_back(tree_K0_r->at(i));
    minitree_K0_NChi2.push_back(tree_K0_NChi2->at(i));
    minitree_K0_mass.push_back(tree_K0_mass->at(i));
    minitree_K0_eta.push_back(tree_K0_eta->at(i));
  }

    for (int i = 0 ; i < tree_L0_x->size() ; i++)
  {
    minitree_K0_x.push_back(tree_L0_x->at(i));
    minitree_K0_y.push_back(tree_L0_y->at(i));
    minitree_K0_z.push_back(tree_L0_z->at(i));
    minitree_K0_r.push_back(tree_L0_r->at(i));
    minitree_K0_NChi2.push_back(tree_L0_NChi2->at(i));
    minitree_K0_mass.push_back(tree_L0_mass->at(i));
    minitree_K0_eta.push_back(tree_L0_eta->at(i));
  }
  
      for (int i = 0 ; i < tree_V0_reco_x->size() ; i++)
    {
      minitree_V0_reco_x.push_back(tree_V0_reco_x->at(i));
      minitree_V0_reco_y.push_back(tree_V0_reco_y->at(i));
      minitree_V0_reco_z.push_back(tree_V0_reco_z->at(i));
      minitree_V0_reco_r.push_back(tree_V0_reco_r->at(i));
      minitree_V0_reco_NChi2.push_back(tree_V0_reco_NChi2->at(i));
      minitree_V0_reco_mass.push_back(tree_V0_reco_mass->at(i));
      minitree_V0_reco_eta.push_back(tree_V0_reco_eta->at(i));
      minitree_V0_reco_source.push_back(tree_V0_reco_source->at(i));
    }
  
      for (int i = 0 ; i < tree_Yc_x->size() ; i++)
    {
      minitree_Yc_x.push_back(tree_Yc_x->at(i));
      minitree_Yc_y.push_back(tree_Yc_y->at(i));
      minitree_Yc_z.push_back(tree_Yc_z->at(i));
      minitree_Yc_r.push_back(tree_Yc_r->at(i));
      minitree_Yc_layer.push_back(tree_Yc_layer->at(i));
      minitree_Yc_NChi2.push_back(tree_Yc_NChi2->at(i));
      minitree_Yc_eta.push_back(tree_Yc_eta->at(i));
      minitree_Yc_mass.push_back(tree_Yc_mass->at(i));
    }
  
    for ( int i = 0 ; i < tree_nSecInt ; i++)
      {
        minitree_SecInt_x.push_back(tree_SecInt_x->at(i));
        minitree_SecInt_y.push_back(tree_SecInt_y->at(i));
        minitree_SecInt_z.push_back(tree_SecInt_z->at(i));
        minitree_SecInt_r.push_back(tree_SecInt_r->at(i));
        minitree_SecInt_d.push_back(tree_SecInt_d->at(i));
        minitree_SecInt_drSig.push_back(tree_SecInt_drSig->at(i));
        minitree_SecInt_dzSig.push_back(tree_SecInt_dzSig->at(i));
        minitree_SecInt_angleXY.push_back(tree_SecInt_angleXY->at(i));
        minitree_SecInt_angleZ.push_back(tree_SecInt_angleZ->at(i));
        minitree_SecInt_NChi2.push_back(tree_SecInt_NChi2->at(i));
        minitree_SecInt_ndf.push_back(tree_SecInt_ndf->at(i));
        minitree_SecInt_mass.push_back(tree_SecInt_mass->at(i));
        minitree_SecInt_pt.push_back(tree_SecInt_pt->at(i));
        minitree_SecInt_eta.push_back(tree_SecInt_eta->at(i));
        minitree_SecInt_phi.push_back(tree_SecInt_phi->at(i));
        minitree_SecInt_charge.push_back(tree_SecInt_charge->at(i));
        minitree_SecInt_badTkHit.push_back(tree_SecInt_badTkHit->at(i));
        minitree_SecInt_dca.push_back(tree_SecInt_dca->at(i));
        minitree_SecInt_selec.push_back(tree_SecInt_selec->at(i));
        minitree_SecInt_layer.push_back(tree_SecInt_layer->at(i));
      }  

    for (unsigned int i = 0 ; i < tree_V0_track_isFromSI->size() ; i++)
      {
        minitree_V0_track_isFromV0.push_back(tree_V0_track_isFromV0->at(i));
        minitree_V0_track_isFromSI.push_back(tree_V0_track_isFromSI->at(i));
        minitree_V0_track_lost.push_back(tree_V0_track_lost->at(i));
        minitree_V0_track_pt.push_back(tree_V0_track_pt->at(i));
        minitree_V0_track_eta.push_back(tree_V0_track_eta->at(i));
        minitree_V0_track_phi.push_back(tree_V0_track_phi->at(i));
        minitree_V0_track_charge.push_back(tree_V0_track_charge->at(i));
        minitree_V0_track_NChi2.push_back(tree_V0_track_NChi2->at(i));
        minitree_V0_track_dxy.push_back(tree_V0_track_dxy->at(i));
        minitree_V0_track_drSig.push_back(tree_V0_track_drSig->at(i));
        minitree_V0_track_dz.push_back(tree_V0_track_dz->at(i));
        minitree_V0_track_dzSig.push_back(tree_V0_track_dzSig->at(i));
        minitree_V0_track_nHit.push_back(tree_V0_track_nHit->at(i));
        minitree_V0_track_nHitPixel.push_back(tree_V0_track_nHitPixel->at(i));
        minitree_V0_track_firstHit.push_back(tree_V0_track_firstHit->at(i));
        minitree_V0_track_firstHit_x.push_back(tree_V0_track_firstHit_x->at(i));
        minitree_V0_track_firstHit_y.push_back(tree_V0_track_firstHit_y->at(i));
        minitree_V0_track_firstHit_z.push_back(tree_V0_track_firstHit_z->at(i));
        minitree_V0_track_iJet.push_back(tree_V0_track_iJet->at(i));
        minitree_V0_track_ntrk10.push_back(tree_V0_track_ntrk10->at(i));
        minitree_V0_track_ntrk20.push_back(tree_V0_track_ntrk20->at(i));
        minitree_V0_track_ntrk30.push_back(tree_V0_track_ntrk30->at(i));
        minitree_V0_track_ntrk40.push_back(tree_V0_track_ntrk40->at(i));
        minitree_V0_track_Hemi.push_back(tree_V0_track_Hemi->at(i));
        minitree_V0_track_Hemi_dR.push_back(tree_V0_track_Hemi_dR->at(i));
        minitree_V0_track_Hemi_dRmax.push_back(tree_V0_track_Hemi_dRmax->at(i));
      }



      //////////// FILLING //////////

      smalltree->Fill();

      //////////// CLEARING //////////

      minirunNumber.clear();
      minieventNumber.clear();
      minilumiBlock.clear();

      minitree_K0_reco_mass.clear();
      minitree_L0_reco_mass.clear();
      minitree_tree_nPV.clear();
      minitree_nSecInt.clear();
      minitree_SecInt_x.clear();
      minitree_SecInt_y.clear();
      minitree_SecInt_z.clear();
      minitree_SecInt_r.clear();
      minitree_SecInt_d.clear();
      minitree_SecInt_drSig.clear();
      minitree_SecInt_dzSig.clear();
      minitree_SecInt_angleXY.clear();
      minitree_SecInt_angleZ.clear();
      minitree_SecInt_NChi2.clear();
      minitree_SecInt_ndf.clear();
      minitree_SecInt_mass.clear();
      minitree_SecInt_pt.clear();
      minitree_SecInt_eta.clear();
      minitree_SecInt_phi.clear();
      minitree_SecInt_charge.clear();
      minitree_SecInt_badTkHit.clear();
      minitree_SecInt_dca.clear();
      minitree_SecInt_selec.clear();
      minitree_SecInt_layer.clear();

      minitree_K0_x.clear();
      minitree_K0_y.clear();
      minitree_K0_z.clear();
      minitree_K0_r.clear();
      minitree_K0_NChi2.clear();
      minitree_K0_mass.clear();
      minitree_K0_eta.clear();

      // CMSSW collection
      minitree_L0_x.clear();
      minitree_L0_y.clear();
      minitree_L0_z.clear();
      minitree_L0_r.clear();
      minitree_L0_NChi2.clear();
      minitree_L0_mass.clear();
      minitree_L0_eta.clear();

      // reco from us 
      minitree_V0_reco_x.clear();
      minitree_V0_reco_y.clear();
      minitree_V0_reco_z.clear();
      minitree_V0_reco_r.clear();
      minitree_V0_reco_NChi2.clear();
      minitree_V0_reco_mass.clear();
      minitree_V0_reco_eta.clear();
      minitree_V0_reco_source.clear();

      //CMSSW collection
      minitree_Yc_x.clear();
      minitree_Yc_y.clear();
      minitree_Yc_z.clear();
      minitree_Yc_r.clear();
      minitree_Yc_layer.clear();
      minitree_Yc_NChi2.clear();
      minitree_Yc_eta.clear();
      minitree_Yc_mass.clear();

      minitree_V0_track_isFromV0.clear();
      minitree_V0_track_isFromSI.clear();
      minitree_V0_track_lost.clear();
      minitree_V0_track_pt.clear();
      minitree_V0_track_eta.clear();
      minitree_V0_track_phi.clear();
      minitree_V0_track_charge.clear();
      minitree_V0_track_NChi2.clear();
      minitree_V0_track_dxy.clear();
      minitree_V0_track_drSig.clear();
      minitree_V0_track_dz.clear();
      minitree_V0_track_dzSig.clear();
      minitree_V0_track_nHit.clear();
      minitree_V0_track_nHitPixel.clear();
      minitree_V0_track_firstHit.clear();
      minitree_V0_track_firstHit_x.clear();
      minitree_V0_track_firstHit_y.clear();
      minitree_V0_track_firstHit_z.clear();
      minitree_V0_track_iJet.clear();
      minitree_V0_track_ntrk10.clear();
      minitree_V0_track_ntrk20.clear();
      minitree_V0_track_ntrk30.clear();
      minitree_V0_track_ntrk40.clear();
      minitree_V0_track_Hemi.clear();
      minitree_V0_track_Hemi_dR.clear();
      minitree_V0_track_Hemi_dRmax.clear();   

      minitree_smu_mass.clear();
      minitree_neu_mass.clear();
      minitree_neu_ctau.clear();  
   } // end loop on events 

   std::cout << "number of events  "<< allevents << std::endl;

   myFile->Write();
   delete myFile;

   // HistogramManager h ;
   // h.WriteAllHistogramsInFile((Production+"/Mini"+sample+".root").Data(),"recreate");
}

