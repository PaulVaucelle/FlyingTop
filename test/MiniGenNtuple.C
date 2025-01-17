#define MiniGenNtuple_cxx
#include "MiniGenNtuple.h"
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


void MiniGenNtuple::Loop(TString sample , TString Production,bool Signal )
{

   // int nthreads = 4;
   // ROOT::EnableImplicitMT(nthreads);

   TFile * myFile = new TFile( (Production+"/MiniGen_"+sample+".root").Data(), "recreate");
   TTree *smalltree = new TTree("ttree", "summary information");
  
    
   std::vector<int>    minitree_smu_mass; // used for signal only
   std::vector<int>    minitree_neu_mass; // used for signal only
   std::vector<float>  minitree_neu_ctau; // used for signal only
   
   std::vector<float>  minitree_genParticle_ct; //minitree_neu_ctau should be the same
   std::vector<float>  minitree_genParticle_ct0;
   std::vector<float>  minitree_genParticle_bg;

   std::vector<float>  minitree_Hemi_dR12;
   std::vector<float>  minitree_genAxis_dRneuneu;
   std::vector<float>  minitree_genAxis_dPhineuneu;
   std::vector<float>  minitree_genAxis_dEtaneuneu;
   std::vector<float>  minitree_GenAxes_Mass;
   std::vector<float>  minitree_GenAxes_CombinedHemiLeptonMass;
   std::vector<float>  minitree_GenAxis_Neu_dRmin;
   std::vector<float>  minitree_GenAxis_Neu_dRmax;
   std::vector<float>  minitree_GenAxis_RecoAxis_dRmin;
   std::vector<float>  minitree_GenAxis_RecoAxis_dRmax;

   //----- Vertex Efficiency variables ---//

   std::vector<float>  minitree_Hemi_LLP_r;
   std::vector<float>  minitree_Hemi_LLP_eta;
   std::vector<float>  minitree_Hemi_LLP_dist;
   std::vector<float>  minitree_Hemi_LLP_ping;

   std::vector<float>  minitree_Hemi_Vtx_NChi2;
   std::vector<int>  minitree_Hemi_Vtx_step;
   std::vector<float>  minitree_Hemi_Vtx_r;
   std::vector<float>  minitree_Hemi_Vtx_z;
   // std::vector<float>  minitree_Hemi_Vtx_eta;
   std::vector<float>  minitree_Hemi_Vtx_dist;
   std::vector<int>  minitree_Hemi_Vtx_nTrks;
   std::vector<float>  minitree_Hemi_Vtx_SumtrackWeight;
   std::vector<float>  minitree_Hemi;
   //--------------------------------------//

   smalltree->Branch("minitree_smu_mass",&minitree_smu_mass); // used for signal only
   smalltree->Branch("minitree_neu_mass",&minitree_neu_mass); // used for signal only
   smalltree->Branch("minitree_neu_ctau",&minitree_neu_ctau); // used for signal only

   smalltree->Branch("minitree_Hemi_dR12",&minitree_Hemi_dR12);
   smalltree->Branch("minitree_genParticle_ct",&minitree_genParticle_ct); // used for signal only
   smalltree->Branch("minitree_genParticle_ct0",&minitree_genParticle_ct0); // used for signal only
   smalltree->Branch("minitree_genParticle_bg",&minitree_genParticle_bg); // used for signal only

   smalltree->Branch("minitree_genAxis_dRneuneu",&minitree_genAxis_dRneuneu);
   smalltree->Branch("minitree_genAxis_dPhineuneu",&minitree_genAxis_dPhineuneu);
   smalltree->Branch("minitree_genAxis_dEtaneuneu",&minitree_genAxis_dEtaneuneu);
   smalltree->Branch("minitree_GenAxes_Mass",&minitree_GenAxes_Mass);
   smalltree->Branch("minitree_GenAxes_CombinedHemiLeptonMass",&minitree_GenAxes_CombinedHemiLeptonMass);
   smalltree->Branch("minitree_GenAxis_Neu_dRmin",&minitree_GenAxis_Neu_dRmin);
   smalltree->Branch("minitree_GenAxis_Neu_dRmax",&minitree_GenAxis_Neu_dRmax);
   smalltree->Branch("minitree_GenAxis_RecoAxis_dRmin",&minitree_GenAxis_RecoAxis_dRmin);
   smalltree->Branch("minitree_GenAxis_RecoAxis_dRmax",&minitree_GenAxis_RecoAxis_dRmax);

   smalltree->Branch("minitree_Hemi_LLP_r",&minitree_Hemi_LLP_r);
   smalltree->Branch("minitree_Hemi_LLP_eta",&minitree_Hemi_LLP_eta);
   smalltree->Branch("minitree_Hemi_LLP_dist",&minitree_Hemi_LLP_dist);
   smalltree->Branch("minitree_Hemi_LLP_ping",&minitree_Hemi_LLP_ping);

   smalltree->Branch("minitree_Hemi_Vtx_NChi2",&minitree_Hemi_Vtx_NChi2);
   smalltree->Branch("minitree_Hemi_Vtx_step",&minitree_Hemi_Vtx_step);
   smalltree->Branch("minitree_Hemi_Vtx_r",&minitree_Hemi_Vtx_r);
   smalltree->Branch("minitree_Hemi_Vtx_z",&minitree_Hemi_Vtx_z);
   // smalltree->Branch("minitree_Hemi_Vtx_eta",&minitree_Hemi_Vtx_eta);
   smalltree->Branch("minitree_Hemi_Vtx_dist",&minitree_Hemi_Vtx_dist);
   smalltree->Branch("minitree_Hemi_Vtx_nTrks",&minitree_Hemi_Vtx_nTrks);
   smalltree->Branch("minitree_Hemi_Vtx_SumtrackWeight",&minitree_Hemi_Vtx_SumtrackWeight);

   smalltree->Branch("minitree_Hemi",&minitree_Hemi);

   std::cout<<"//-------------------------//"<<std::endl;
   std::cout<<"// "<<Production<<" : "<<sample<<"  //"<<std::endl;
   std::cout<<"//-------------------------//"<<std::endl;

   bool debug = false;
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "Total Entries : " << nentries << std::endl;

   Long64_t nbytes = 0, nb = 0;

   int allevents = 0;
   int itest = 0;

   for (Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;//break
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      allevents++;
      // if ( allevents%10000 == 0 ) // std::cout << "events : " << allevents << std::endl;

      if (Cut(ientry) < 0) continue;
      // std::cout<<"jentry "<<jentry<<std::endl;
      if ( jentry%10000 == 0 )  std::cout << "events : " << jentry << std::endl;
      // if (jentry==0.1*nentries) {// std::cout<<"10/100 :"<<std::endl;}
      // if (jentry==0.2*nentries) {// std::cout<<"20/100 :"<<std::endl;}
      // if (jentry==0.3*nentries) {// std::cout<<"30/100 :"<<std::endl;}
      // if (jentry==0.4*nentries) {// std::cout<<"40/100 :"<<std::endl;}
      // if (jentry==0.5*nentries) {// std::cout<<"50/100 :"<<std::endl;}
      // if (jentry==0.6*nentries) {// std::cout<<"60/100 :"<<std::endl;}
      // if (jentry==0.7*nentries) {// std::cout<<"70/100 :"<<std::endl;}
      // if (jentry==0.8*nentries) {// std::cout<<"80/100 :"<<std::endl;}
      // if (jentry==0.9*nentries) {// std::cout<<"90/100 :"<<std::endl;}

//$$ 
      itest++;
//    if ( itest > 1000000 ) break;
//$$

      //----------------------//
      //        Event         //
      //----------------------//
    // -------------- --------------------------------------------------------------------//
    // -----------------------------------------------------------------------------------//
    if ( !((tree_Filter || tree_FilterSameSign) && tree_njetNOmu>0) && !Signal ) continue;
    //------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------//
//$$

      // miniAllPU_events_weight.push_back(AllPU_events_weight);

      if (Signal)
         {
            // if (debug)  { std::cout<<"Signal parameters "<<std::endl;}
            minitree_smu_mass.push_back(tree_smu_mass); // used for signal only
            minitree_neu_mass.push_back(tree_neu_mass); // used for signal only
            minitree_neu_ctau.push_back(tree_neu_ctau); // used for signal only

            for (unsigned int i = 0 ; i < tree_genParticle_ct->size(); i++)
               {
                  if ( tree_genParticle_ct->at(i)<=0 || tree_genParticle_ct0->at(i)<=0 ) continue;
                  minitree_genParticle_ct.push_back(tree_genParticle_ct->at(i));
                  minitree_genParticle_ct0.push_back(tree_genParticle_ct0->at(i));
                  float bg = tree_genParticle_ct->at(i)/tree_genParticle_ct0->at(i);
                  minitree_genParticle_bg.push_back(bg);
               }
            if (tree_nLLP>=2 && tree_genAxis_dRneuneu->size()>0)
               {  
                  // std::cout<<"tree_nLLP "<<tree_nLLP<<std::endl;
                  minitree_genAxis_dRneuneu.push_back(tree_genAxis_dRneuneu->at(0));
                  // std::cout<<"2 "<<std::endl;
                  minitree_genAxis_dPhineuneu.push_back(tree_genAxis_dPhineuneu->at(0));
                  // std::cout<<"3 "<<std::endl;
                  minitree_genAxis_dEtaneuneu.push_back(tree_genAxis_dEtaneuneu->at(0));
                  // std::cout<<"4 "<<std::endl;

               }
            if (tree_nLLP>=2 && tree_GenAxes_Mass->size()>1)
               {
                  minitree_GenAxes_Mass.push_back(tree_GenAxes_Mass->at(0));
                  // std::cout<<"5 "<<std::endl;
                  minitree_GenAxes_CombinedHemiLeptonMass.push_back(tree_GenAxes_CombinedHemiLeptonMass->at(0));

                  // std::cout<<"6 "<<std::endl;
                  minitree_GenAxes_Mass.push_back(tree_GenAxes_Mass->at(1));
                  // std::cout<<"7 "<<std::endl;
                  minitree_GenAxes_CombinedHemiLeptonMass.push_back(tree_GenAxes_CombinedHemiLeptonMass->at(1));


               }
            if (tree_nLLP>=2 && tree_GenAxis_Neu_dRmin->size()>1)
               {
                  // std::cout<<"8 "<<std::endl;
                  minitree_GenAxis_Neu_dRmin.push_back(tree_GenAxis_Neu_dRmin->at(0));
                  // std::cout<<"9 "<<std::endl;
                  minitree_GenAxis_Neu_dRmax.push_back(tree_GenAxis_Neu_dRmax->at(0));
                  // std::cout<<"10 "<<std::endl;
                  minitree_GenAxis_Neu_dRmin.push_back(tree_GenAxis_Neu_dRmin->at(1));
                  // std::cout<<"11 "<<std::endl;
                  minitree_GenAxis_Neu_dRmax.push_back(tree_GenAxis_Neu_dRmax->at(1));
               }
            if ( tree_nLLP>=2 && tree_GenAxis_RecoAxis_dRmin->size()>1)
               {
                  // std::cout<<"12 "<<std::endl;
                  minitree_GenAxis_RecoAxis_dRmin.push_back(tree_GenAxis_RecoAxis_dRmin->at(0));
                  // std::cout<<"13 "<<std::endl;
                  minitree_GenAxis_RecoAxis_dRmax.push_back(tree_GenAxis_RecoAxis_dRmax->at(0));
                  // std::cout<<"14 "<<std::endl;
                   minitree_GenAxis_RecoAxis_dRmin.push_back(tree_GenAxis_RecoAxis_dRmin->at(1));
                  // std::cout<<"15 "<<std::endl;
                  minitree_GenAxis_RecoAxis_dRmax.push_back(tree_GenAxis_RecoAxis_dRmax->at(1));
                  // std::cout<<"16 "<<std::endl;
               }
            for (unsigned int i = 0 ; i < tree_Hemi_dR12->size(); i++)
               {
                  minitree_Hemi_dR12.push_back(tree_Hemi_dR12->at(i));

                  float LLP_r = sqrt(tree_Hemi_LLP_x->at(i)*tree_Hemi_LLP_x->at(i) + tree_Hemi_LLP_y->at(i)*tree_Hemi_LLP_y->at(i));
                  minitree_Hemi_LLP_r.push_back(LLP_r);
                  minitree_Hemi_LLP_eta.push_back(tree_Hemi_LLP_eta->at(i));
                  minitree_Hemi_LLP_dist.push_back(tree_Hemi_LLP_dist->at(i));
                  minitree_Hemi_LLP_ping.push_back(tree_Hemi_LLP_ping->at(i));

                  minitree_Hemi_Vtx_NChi2.push_back(tree_Hemi_Vtx_NChi2->at(i));
                  minitree_Hemi_Vtx_step.push_back(tree_Hemi_Vtx_step->at(i));
                  minitree_Hemi_Vtx_r.push_back(tree_Hemi_Vtx_r->at(i));
                  minitree_Hemi_Vtx_z.push_back(tree_Hemi_Vtx_z->at(i));
                  // minitree_Hemi_Vtx_eta.push_back(tree_Hemi_Vtx_eta->at(i));
                  minitree_Hemi_Vtx_dist.push_back(tree_Hemi_Vtx_dist->at(i));
                  minitree_Hemi_Vtx_nTrks.push_back(tree_Hemi_Vtx_nTrks->at(i));
                  minitree_Hemi_Vtx_SumtrackWeight.push_back(tree_Hemi_Vtx_SumtrackWeight->at(i));

                  minitree_Hemi.push_back(tree_Hemi->at(i));
               }
         }
      else // for background, we cna add vtx information :/
         {
            // if (debug)  { std::cout<<"Background parameters "<<std::endl;}
            for (unsigned int i = 0 ; i < tree_Hemi->size(); i++)
               {
                  minitree_Hemi.push_back(tree_Hemi->at(i));
                  minitree_Hemi_Vtx_NChi2.push_back(tree_Hemi_Vtx_NChi2->at(i));
                  minitree_Hemi_Vtx_step.push_back(tree_Hemi_Vtx_step->at(i));
                  minitree_Hemi_Vtx_r.push_back(tree_Hemi_Vtx_r->at(i));
                  minitree_Hemi_Vtx_z.push_back(tree_Hemi_Vtx_z->at(i));
                  // minitree_Hemi_Vtx_eta.push_back(tree_Hemi_Vtx_eta->at(i));
                  minitree_Hemi_Vtx_dist.push_back(tree_Hemi_Vtx_dist->at(i));
                  minitree_Hemi_Vtx_nTrks.push_back(tree_Hemi_Vtx_nTrks->at(i));
                  minitree_Hemi_Vtx_SumtrackWeight.push_back(tree_Hemi_Vtx_SumtrackWeight->at(i));
               }
         }



      // std::cout<<"do we go here ? : "<<jentry<<std::endl;
      //////////// FILLING //////////

      smalltree->Fill();

      //////////// CLEARING //////////


      minitree_smu_mass.clear();
      minitree_neu_mass.clear();
      minitree_neu_ctau.clear();
      minitree_genParticle_ct.clear();
      minitree_genParticle_ct0.clear();
      minitree_genParticle_bg.clear();

      minitree_Hemi_dR12.clear();
      minitree_genAxis_dRneuneu.clear();
      minitree_genAxis_dPhineuneu.clear();
      minitree_genAxis_dEtaneuneu.clear();
      minitree_GenAxes_Mass.clear();
      minitree_GenAxes_CombinedHemiLeptonMass.clear();
      minitree_GenAxis_Neu_dRmin.clear();
      minitree_GenAxis_Neu_dRmax.clear();
      minitree_GenAxis_RecoAxis_dRmin.clear();
      minitree_GenAxis_RecoAxis_dRmax.clear();

      minitree_Hemi_dR12.clear();

      minitree_Hemi_LLP_r.clear();
      minitree_Hemi_LLP_eta.clear();
      minitree_Hemi_LLP_dist.clear();
      minitree_Hemi_LLP_ping.clear();

      minitree_Hemi_Vtx_NChi2.clear();
      minitree_Hemi_Vtx_step.clear();
      minitree_Hemi_Vtx_r.clear();
      minitree_Hemi_Vtx_z.clear();
      // minitree_Hemi_Vtx_eta.clear();
      minitree_Hemi_Vtx_dist.clear();
      minitree_Hemi_Vtx_nTrks.clear();
      minitree_Hemi_Vtx_SumtrackWeight.clear();

      minitree_Hemi.clear();
   } // end loop on events 

   std::cout << "number of events  "<< allevents << std::endl;

   myFile->Write();
   delete myFile;
}
