#define TreeReader_cxx
#include "TreeReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

void TreeReader::Loop(TString sample)//void TreeReader::Loop(TString sample)
{
//   In a ROOT session, you can do:
//      root> .L TreeReader.C
//      root> TreeReader t
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
   TString thesample  = sample;
   TFile * theoutputfile = new TFile( ("outputroot/histofile_"+thesample+".root").Data() , "recreate");
   std::ofstream ofs ("outputroot/Efficacity_"+thesample+".txt", std::ofstream::out);

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
   // ----------- Parameters --------------------//

   float TightWP = 0.85;
   float LooseWP = 0.;

   // -------- List of variables ---------------//

   float n_TightTrks = 0;
   float n_LooseTrks = 0;
   float n_TotalTrks = 0;
   float TightTrks_Eff = 0.;
   float LooseTrks_Eff = 0.;

   float nRecoVertex = 0;
   float MeanDistance = 0;
   float nRecoVertexTightWP = 0;
   float MeanDistanceTightWP = 0;
   float nRecoVertexLooseWP = 0;
   float MeanDistanceLooseWP = 0;
   float TotalnVertex = 0;
   float nEvts = 0;

//--------VtxSelection Variables-------

float MeanTWcut = 1.5;
float dTW = 0.001;
int nSteps = 1/dTW;
float nSelecVtx[1000] = {0};
float nSelecVtxStep1[1000] = {0};
//--------ntrk10------
int ntrk10Cut = 0;
int dntrk10 = 1;
int nStep_10 = (50-ntrk10Cut)/dntrk10;
float nEvts_ntrk10[50] = {0};
float nEvts_ntrk10_step1[50] = {0};
//-----------------------------
//--------ntrk20------
int ntrk20Cut = 0;
int dntrk20 = 1;
int nStep_20 = (50-ntrk20Cut)/dntrk20;
float nEvts_ntrk20[50] = {0};
float nEvts_ntrk20_step1[50] = {0};
//-----------------------------
//--------Vtx_dd------
int VtxddCut = 0;
int ddd = 1;
int nStep_dd = (250-VtxddCut)/ddd;
float nEvts_Vtx_dd[250] = {0};
float nEvts_Vtx_dd_step1[250] = {0};
//-----------------------------
//--------Vtx_nTrks------
int VtxnTrksCut = 2;
int dnTrks = 1;
int nStep_nTrks = (52-VtxddCut)/dnTrks;
float nEvts_Vtx_nTrks[50] = {0};
float nEvts_Vtx_nTrks_step1[50] = {0};
//-----------------------------
//--------Vtx_nVtx------
int VtxnVtxCut = 0;
int dnVtx = 1;
int nStep_nVtx = (4-VtxnVtxCut)/dnVtx;
float nEvts_Vtx_nVtx[4] = {0};
float nEvts_Vtx_nVtx_step1[4] = {0};
//-----------------------------------------------

//--------Vtx_Chi2------
float VtxChi2Cut = 0;
float dChi2 = 0.25;
int nStep_Chi2 = (10-VtxChi2Cut)/dChi2;
float nEvts_Vtx_Chi2[50] = {0};
float nEvts_Vtx_Chi2_step1[50] = {0};
//-----------------------------------------------

//--------Vtx_Step------
int VtxStepCut = 0;
int dStep = 1;
int nStep_Step = (4-VtxStepCut)/dStep;
float nEvts_Vtx_Step[4] = {0};
float nEvts_Vtx_Step_step1[4] = {0};
//-----------------------------------------------

//--------Vtx_InvMass------
int VtxMassCut = 0;
int dMass = 10;
int nStep_Mass = (14000-VtxMassCut)/dMass;
float nEvts_Vtx_Mass[1400] = {0};
float nEvts_Vtx_Mass_step1[1400] = {0};
//-----------------------------------------------

//--------Hemi_InvMass------
int HemiMassCut = 0;
int dHemiMass = 10;
int nStep_HemiMass = (1400-HemiMassCut)/dHemiMass;
float nEvts_Hemi_Mass[140] = {0};
float nEvts_Hemi_Mass_step1[140] = {0};
//-----------------------------------------------

//--------Vtx_DCA tracks------
float VtxDCACut = 0;
float dDCA = 1;
int nStep_DCA = (100-VtxDCACut)/dDCA;
float nEvts_Vtx_DCA[100] = {0};
float nEvts_Vtx_DCA_step1[100] = {0};
//-----------------------------------------------

//--------------------EvtsSelection Varaibles----
//--------Leadingjet pt------
int JetPTCut = 20;
int dpt = 10;
int nStep_jet_pt = (1020-JetPTCut)/dpt;//1000-20/dpt
float nEvts_jetpt[100] = {0};
//---------------------------
//--------subLeadingjet pt------
int Jet2PTCut = 20;
int dpt2 = 10;
int nStep_jet2_pt = (1020-Jet2PTCut)/dpt2;//1000-20/dpt
float nEvts_jetpt2[100] = {0};
//---------------------------
//--------Leadingmuon pt------
int MuonPTCut = 20;
int dpt_muon = 10;
int nStep_muon_pt = (520-MuonPTCut)/dpt_muon;//1000-20/dpt
float nEvts_muonpt[50] = {0};
//---------------------------
//--------subLeading muon pt------
int Muon2PTCut = 10;
int dpt_muon2 = 5;
int nStep_muon2_pt = (510-Muon2PTCut)/dpt_muon2;//1000-20/dpt
float nEvts_muon2pt[100] = {0};
//---------------------------
//--------HT------
int HTCut = 100;
int dpt_HT = 10;
int nStep_HT_pt = (2100-HTCut)/dpt_HT;//1000-20/dpt
float nEvts_HT[200] = {0};
//---------------------------
//--------ST------
int STCut = 30;
int dpt_ST = 10;
int nStep_ST_pt = (530-STCut)/dpt_ST;//1000-20/dpt
float nEvts_ST[50] = {0};
//---------------------------
//--------nJet------
int nJetCut = 0;
int dnJet_ST = 1;
int nStep_nJet = (20-nJetCut)/dnJet_ST;//1000-20/dpt
float nEvts_nJet[20] = {0};

//-------------nmu----------//
int nMuonCut = 0;
int dnMuon = 1;
int nStep_nMuon = (20-nMuonCut)/dnMuon;//1000-20/dpt
float nEvts_nMuon[20] = {0};
// //------------------------------

// //-------------muon_muon_dR----------//
float MuondRCut = 0;
float dMuondR = 0.125;
int nStep_MuondR = (5.-MuondRCut)/dMuondR;//1000-20/dpt
float nEvts_MuondR[40] = {0};
// //------------------------------

// //-------------muon_muon_dPhi----------//
float MuondPhiCut = 0.;
float dMuondPhi = 0.125;
int nStep_MuondPhi = (3.5-MuondPhiCut)/dMuondPhi;//1000-20/dpt
float nEvts_MuondPhi[28] = {0};
// //------------------------------

// //-------------muon_muon_dEta----------//
float MuondEtaCut = 0.;
float dMuondEta = 0.125;
int nStep_MuondEta = (4-MuondEtaCut)/dMuondEta;//1000-20/dpt
float nEvts_MuondEta[32] = {0};
// //------------------------------

// //-------------jet_jet_dR----------//
float JetdRCut = 0.4;
float dJetdR = 0.125;
int nStep_JetdR = (5.4-JetdRCut)/dJetdR;//1000-20/dpt
float nEvts_JetdR[40] = {0};
// //------------------------------

// //-------------jet_jet_dPhi----------//
float JetdPhiCut = 0.;
float dJetdPhi = 0.125;
int nStep_JetdPhi = (3.5-JetdPhiCut)/dJetdPhi;//1000-20/dpt
float nEvts_JetdPhi[28] = {0};
// //------------------------------

// //-------------jet_jet_dEta----------//
float JetdEtaCut = 0;
float dJetdEta = 0.125;
int nStep_JetdEta = (5-JetdEtaCut)/dJetdEta;//1000-20/dpt
float nEvts_JetdEta[40] = {0};
// //------------------------------

// //-------------muon_jet_dRmin----------//
float MuonJetdRminCut = 0;
float dMuonJetdRmin = 0.125;
int nStep_MuonJetdRmin = (5-MuonJetdRminCut)/dMuonJetdRmin;//1000-20/dpt
float nEvts_MuonJetdRmin[40] = {0};
// //------------------------------
// //-------------muon_jet_dRmax----------//
float MuonJetdRCutmax = 0;
float dMuonJetdRmax = 0.125;
int nStep_MuonMuonJetdRmax = (6-MuonJetdRCutmax)/dMuonJetdRmax;//1000-20/dpt
float nEvts_MuonJetdRmax[48] = {0};
//------------------------------
   // -------------End of list of variables ----------//

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

cout<< "Line : "  << __LINE__ << " " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      // std::cout<<"event : "<<jentry<<std::endl;
      fillHisto("Filter", "Offline+Online", thesample, tree_Filter, 1);

   if (tree_Filter)
      {
         nEvts++;
      //*******************************
      //loop on Muons
      //*******************************

      for(unsigned int iMuon = 0; iMuon <tree_muon_pt->size(); iMuon ++)
         {
            fillHisto("RecoMuo_pT", "noSel",  thesample,  tree_muon_pt->at(iMuon) , 1.);
            //add the selection
         }

      //*******************************
      //loop on V0 Candidates
      //*******************************

      for (unsigned int iV0 = 0; iV0 <tree_V0_reco_source->size(); iV0++)
         {
            if(tree_V0_reco_source->at(iV0)==1)//K0
               {
                  fillHisto("hData_reco_K0_mass","noSel",  thesample, tree_V0_reco_mass->at(iV0),1.);
                  fillHisto("hData_reco_K0_r","noSel",  thesample,    tree_V0_reco_r->at(iV0),1.);
                  fillHisto("hData_reco_K0_z","noSel",  thesample,    tree_V0_reco_z->at(iV0),1.);
               }

            if(tree_V0_reco_source->at(iV0)==2)//L0
               {

                  fillHisto("hData_reco_L0_mass","noSel",  thesample, tree_V0_reco_mass->at(iV0),1.);
                  fillHisto("hData_reco_L0_r","noSel",   thesample,   tree_V0_reco_r->at(iV0),1.);
                  fillHisto("hData_reco_L0_z","noSel",    thesample,  tree_V0_reco_z->at(iV0),1.);
               }

         }

      //*******************************
      //loop on Sec. Interactions

      for (unsigned int iSecInt = 0; iSecInt <tree_SecInt_mass->size(); iSecInt++)
         {
            fillHisto("hData_reco_SecInt_mass","noSel",        thesample,tree_SecInt_mass->at(iSecInt),1.);
            fillHisto("hData_reco_SecInt_x","noSel",           thesample, tree_SecInt_x->at(iSecInt),1.);
            fillHisto("hData_reco_SecInt_y","noSel",           thesample, tree_SecInt_y->at(iSecInt),1.);
            fillHisto("hData_reco_SecInt_r","noSel",           thesample,  tree_SecInt_r->at(iSecInt),1.);
            fillHisto("hData_reco_SecInt_z","noSel",           thesample, tree_SecInt_z->at(iSecInt),1.);

            if (tree_SecInt_selec->at(iSecInt))
               {
                  fillHisto("hData_reco_SecInt_mass","Selec",  thesample, tree_SecInt_mass->at(iSecInt),1.);
                  fillHisto("hData_reco_SecInt_x","Selec",     thesample, tree_SecInt_x->at(iSecInt),1.);
                  fillHisto("hData_reco_SecInt_y","Selec",     thesample, tree_SecInt_y->at(iSecInt),1.);
                  fillHisto("hData_reco_SecInt_r","Selec",     thesample, tree_SecInt_r->at(iSecInt),1.);
                  fillHisto("hData_reco_SecInt_z","Selec",     thesample,  tree_SecInt_z->at(iSecInt),1.);

                  if(tree_SecInt_layer->at(iSecInt)) 
                     {
                        fillHisto("hData_reco_SecInt_mass","TrackerMatched",  thesample, tree_SecInt_mass->at(iSecInt),1.);
                        fillHisto("hData_reco_SecInt_x","TrackerMatched",     thesample, tree_SecInt_x->at(iSecInt),1.);
                        fillHisto("hData_reco_SecInt_y","TrackerMatched",     thesample, tree_SecInt_y->at(iSecInt),1.);
                        fillHisto("hData_reco_SecInt_r","TrackerMatched",     thesample,tree_SecInt_r->at(iSecInt),1.);
                        fillHisto("hData_reco_SecInt_z","TrackerMatched",     thesample,tree_SecInt_z->at(iSecInt),1.);
                     }                
               }
         }

      //*******************************
      //loop on tracks
      //*******************************

      n_TotalTrks += tree_track_MVAval->size();

      for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++)
         {
          //--- Track pre-selected ---------//
          //get total number of tracks bfore any selection with tree_TRACK_SIZE

          //----------Tracks selected by the BDT
            fillHisto("hData_MVAVal","noSel", thesample,tree_track_MVAval->at(iTrk),1.);
            if ( tree_track_MVAval->at(iTrk) > TightWP )
               {
                  n_TightTrks++;
               }

            if ( tree_track_MVAval->at(iTrk) > LooseWP )
               {
                  n_LooseTrks++;
               }

         }

      TightTrks_Eff = n_TightTrks/n_TotalTrks;
      LooseTrks_Eff = n_LooseTrks/n_TotalTrks;

      //*******************************
      //loop on Axes
      //*******************************
      
      for(unsigned int iAxes = 0; iAxes < tree_Hemi_dR12->size(); iAxes ++)
         {
            fillHisto("hData_dR_RecoReco","noSel",  thesample, tree_Hemi_dR12->at(iAxes),1.);
         }

      //*******************************
      //loop on Reco Vertices
      //*******************************
 
      for(unsigned int iVtx = 0; iVtx <tree_Hemi->size(); iVtx ++)
         {
            fillHisto("hData_Hemi_Vtx_NChi2","noSel",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
            fillHisto("hData_Hemi_Vtx_nTrks","noSel",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
            TotalnVertex++;

            if (tree_Hemi_Vtx_NChi2->at(iVtx) != -10 )
               {
                  fillHisto("hData_Hemi_Vtx_NChi2","RecoVtx",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_nTrks","RecoVtx",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_dist","RecoVtx",   thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_step","RecoVtx",   thesample,tree_Hemi_Vtx_step->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_eta","RecoVtx",        thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_r","RecoVtx",   thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_z","RecoVtx",   thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_MVAval","RecoVtx",   thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);
                  

                  if ( tree_Hemi_Vtx_NChi2->at(iVtx)>0 && tree_Hemi_Vtx_NChi2->at(iVtx)<10)//Reco Vtx criteria
                     {
                        fillHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_dist","GoodRecoVtx",   thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_step","GoodRecoVtx",   thesample,tree_Hemi_Vtx_step->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_eta","GoodRecoVtx",    thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_r","GoodRecoVtx",   thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_z","GoodRecoVtx",   thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx",   thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);
                        nRecoVertex++;
                        MeanDistance +=  tree_Hemi_Vtx_dist->at(iVtx);              
                     }
                  

                  // -------------- TIGHT WP ------------------------//
                  if (tree_Hemi_Vtx_step->at(iVtx)==1 || tree_Hemi_Vtx_step->at(iVtx)==2)
                     {
                        fillHisto("hData_Hemi_Vtx_nTrks","RecoVtx_TightWP",  thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_NChi2","RecoVtx_TightWP",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_dist","RecoVtx_TightWP",   thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_eta","RecoVtx_TightWP",   thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_r","RecoVtx_TightWP",   thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_z","RecoVtx_TightWP",   thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_MVAval","RecoVtx_TightWP",   thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);
                       
                        
                        if ( tree_Hemi_Vtx_NChi2->at(iVtx)>0 && tree_Hemi_Vtx_NChi2->at(iVtx)<10)//Reco Vtx criteria
                           {
                              fillHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx_TightWP",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx_TightWP",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_dist","GoodRecoVtx_TightWP",   thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_eta","GoodRecoVtx_TightWP",   thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_r","GoodRecoVtx_TightWP",   thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_z","GoodRecoVtx_TightWP",   thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx_TightWP",   thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);
                              nRecoVertexTightWP++;
                              MeanDistanceTightWP +=  tree_Hemi_Vtx_dist->at(iVtx);              
                      
                           }
                     }// End Tight WP

                  // -------------- Loose WP ------------------------//
                  if (tree_Hemi_Vtx_step->at(iVtx)==3 || tree_Hemi_Vtx_step->at(iVtx)==4)
                     {
                        fillHisto("hData_Hemi_Vtx_nTrks","RecoVtx_LooseWP",  thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_NChi2","RecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_dist","RecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_eta","RecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_r","RecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_z","RecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_MVAval","RecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);

                        
                        if ( tree_Hemi_Vtx_NChi2->at(iVtx)>0 && tree_Hemi_Vtx_NChi2->at(iVtx)<10)//Reco Vtx criteria
                           {
                              fillHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_dist","GoodRecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_eta","GoodRecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_r","GoodRecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_z","GoodRecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);  
                              nRecoVertexLooseWP++;
                              MeanDistanceLooseWP +=  tree_Hemi_Vtx_dist->at(iVtx);              
                      
                           }
                     }//Loose WP

               }//RecoVtx
         }//End Loop on Vertices




         //---------------Vtx Selection Variables----------------//
  //------------------MWT CUT -------------------//
  for (int i =0 ;i<nSteps+1;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_MeantrackWeight->at(j)>(MeanTWcut+i*dTW) )
            {
              nSelecVtx[i]++;
              
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nSelecVtxStep1[i]++;
                }
            }
        }
    //-----------------------------------//
    }


//-----------------------------

//--------ntrk10------
//--------ntrk20------
// std::cout<<"debug 1"<<std::endl;
  for (int i =0 ;i<nStep_10;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_ntrk10->at(j)>=(ntrk10Cut+i*dntrk10))
            {
              nEvts_ntrk10[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nEvts_ntrk10_step1[i]++;
                }
            }
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_ntrk20->at(j)>=(ntrk20Cut+i*dntrk20))
            {
              nEvts_ntrk20[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nEvts_ntrk20_step1[i]++;
                }
            }
        }
    }
// std::cout<<"debug 2"<<std::endl;
// --------Vtx_dd------

  for (int i =0 ;i<nStep_dd;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_Vtx_dd->at(j)>=(VtxddCut+i*ddd))
            {
              nEvts_Vtx_dd[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2 )
                {
                  nEvts_Vtx_dd_step1[i]++;
                }
            }
        }
    }
// std::cout<<"debug 3"<<std::endl;
    //--------Vtx_nTrks------
//-----------------------------

  for (int i =0 ;i<nStep_nTrks;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_nTrks->at(j)>=(VtxnTrksCut+i*dnTrks))
            {
              nEvts_Vtx_nTrks[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nEvts_Vtx_nTrks_step1[i]++;
                }
            }
        }
    }
// std::cout<<"debug 4"<<std::endl;
//--------Vtx_nVtx------
  for (int i =0 ;i<nStep_nVtx;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi_Vtx_nVtx->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_nVtx->at(j)==(VtxnVtxCut+i*dnVtx))
            {
              nEvts_Vtx_nVtx[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2 )
                {
                  nEvts_Vtx_nVtx_step1[i]++;
                }
            }
        }
    }
//-----------------------------------------------

//--------Vtx_Chi2------

// float nEvts_Vtx_Chi2_step1[100] = {0};
  for (int i =0 ;i<nStep_Chi2;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_NChi2->at(j)>=(VtxChi2Cut+i*dChi2))
            {
              nEvts_Vtx_Chi2[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nEvts_Vtx_Chi2_step1[i]++;
                }
            }
        }
    }
//-----------------------------------------------

//--------Vtx_Step------

  for (int i =0 ;i<nStep_Step;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_step->at(j)<=(VtxStepCut+i*dStep))
            {
              nEvts_Vtx_Step[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nEvts_Vtx_Step_step1[i]++;
                }
            }
        }
    }
//-----------------------------------------------

//--------Vtx_InvMass------
  for (int i =0 ;i<nStep_Mass;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_Mass->at(j)>=(VtxMassCut+i*dMass))
            {
              nEvts_Vtx_Mass[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nEvts_Vtx_Mass_step1[i]++;
                }
            }
        }
    }
//-----------------------------------------------
//--------Hemi_InvMass------
  for (int i =0 ;i<nStep_HemiMass;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Mass->at(j)>=(HemiMassCut+i*dHemiMass))
            {
              nEvts_Hemi_Mass[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nEvts_Hemi_Mass_step1[i]++;
                }
            }
        }
    }

//--------Vtx_DCA tracks------

  for (int i =0 ;i<nStep_DCA;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_nTrks->at(j)>=(VtxDCACut+i*dDCA))
            {
              nEvts_Vtx_DCA[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nEvts_Vtx_DCA_step1[i]++;
                }
            }
        }
    }


  //-------------------------------------------------------------------//
  //--------------------------EVTS Selection Variables-----------------//
//--------Leadingjet pt------
for (int i = 0 ; i < nStep_jet_pt ; i++)
  {
    for (unsigned int j = 0 ; j < tree_jet_leadingpt->size() ; j++)
      {
        if (tree_jet_leadingpt->at(j) > (JetPTCut+i*dpt))
          {
            nEvts_jetpt[i]++;
          }
      }
  }

//--------subLeadingjet pt------
for (int i = 0 ; i < nStep_jet2_pt ; i++)
  {
    for (unsigned int j = 0 ; j < tree_jet_leadingpt2->size() ; j++)
      {
        if (tree_jet_leadingpt2->at(j) > (Jet2PTCut+i*dpt2))
          {
            nEvts_jetpt2[i]++;
          }
      }
  }

//--------Leadingmuon pt------
for (int i = 0 ; i < nStep_muon_pt ; i++)
  {
    for (unsigned int j = 0 ; j < tree_muon_leadingpt->size() ; j++)
      {
        if (tree_muon_leadingpt->at(j) > (MuonPTCut+i*dpt_muon))
          {
            nEvts_muonpt[i]++;
          }
      }
  }

//--------subLeading muon pt------
for (int i = 0 ; i < nStep_muon2_pt ; i++)
  {
    for (unsigned int j = 0 ; j < tree_muon_leadingpt2->size() ; j++)
      {
        if (tree_muon_leadingpt2->at(j) >= (Muon2PTCut+i*dpt_muon2))
          {
            nEvts_muon2pt[i]++;
          }
      }
  }

//--------HT------
for (int i = 0 ; i < nStep_HT_pt ; i++)
  {
    // for (unsigned int j = 0 ; j < tree_ST->size() ; j++)
    //   {
        if (tree_HT > (HTCut+i*dpt_HT))
          {
            nEvts_HT[i]++;
          }
      // }
  }

//--------ST------

for (int i = 0 ; i < nStep_ST_pt ; i++)
  {
    for (unsigned int j = 0 ; j < tree_ST->size() ; j++)
      {
        if (tree_ST->at(j) > (STCut+i*dpt_ST))
          {
            nEvts_ST[i]++;
          }
      }
  }

//--------nJets
for (int i = 0 ; i < nStep_nJet ; i++)
  {
    for (unsigned int j = 0 ; j < tree_ST->size() ; j++)
      {
        if (tree_njet > (nJetCut+i*dnJet_ST))
          {
            nEvts_nJet[i]++;
          }
      }
  }

  //-------------nmu----------//

for (int i = 0 ; i < nStep_nMuon ; i++)
  {
    for (unsigned int j = 0 ; j < tree_muon_nmu->size() ; j++)
      {
        if (tree_muon_nmu->at(j) > (nMuonCut+i*dnMuon))
          {
            nEvts_nMuon[i]++;
          }
      }
  }
// //------------------------------

// //-------------muon_muon_dR----------//

for (int i = 0 ; i < nStep_MuondR ; i++)
  {
    for (unsigned int j = 0 ; j < tree_muon_muon_dR->size() ; j++)
      {
        if (tree_muon_muon_dR->at(j) > (MuondRCut+i*dMuondR))
          {
            nEvts_MuondR[i]++;
          }
      }
  }
// //------------------------------

// //-------------muon_muon_dPhi----------//

for (int i = 0 ; i < nStep_MuondPhi ; i++)
  {
    for (unsigned int j = 0 ; j < tree_muon_muon_dPhi->size() ; j++)
      {
        if (tree_muon_muon_dPhi->at(j) > (MuondPhiCut+i*dMuondPhi))
          {
            nEvts_MuondPhi[i]++;
          }
      }
  }
// //------------------------------

// //-------------muon_muon_dEta----------//

for (int i = 0 ; i < nStep_MuondEta ; i++)
  {
    for (unsigned int j = 0 ; j < tree_muon_muon_dEta->size() ; j++)
      {
        if (tree_muon_muon_dEta->at(j) > (MuondEtaCut+i*dMuondEta ))
          {
            nEvts_MuondEta[i]++;
          }
      }
  }
// //------------------------------

// //-------------jet_jet_dR----------//

for (int i = 0 ; i < nStep_JetdR ; i++)
  {
    for (unsigned int j = 0 ; j < tree_jet_jet_dR->size() ; j++)
      {
        if (tree_jet_jet_dR->at(j) > (JetdRCut+i*dJetdR))
          {
            nEvts_JetdR[i]++;
          }
      }
  }
// //------------------------------

// //-------------jet_jet_dPhi----------//

for (int i = 0 ; i < nStep_JetdPhi ; i++)
  {
    for (unsigned int j = 0 ; j < tree_jet_jet_dPhi->size() ; j++)
      {
        if (tree_jet_jet_dPhi->at(j) > (JetdPhiCut+i*dJetdPhi))
          {
            nEvts_JetdPhi[i]++;
          }
      }
  }
// //------------------------------

// //-------------jet_jet_dEta----------//

for (int i = 0 ; i < nStep_JetdEta ; i++)
  {
    for (unsigned int j = 0 ; j < tree_jet_jet_dEta->size() ; j++)
      {
        if (tree_jet_jet_dEta->at(j) > (JetdEtaCut+i*dJetdEta))
          {
            nEvts_JetdEta[i]++;
          }
      }
  }
// //------------------------------

// //-------------muon_jet_dRmin----------//

for (int i = 0 ; i < nStep_MuonJetdRmin ; i++)
  {
    for (unsigned int j = 0 ; j < tree_muon_jet_dRmin->size() ; j++)
      {
        if (tree_muon_jet_dRmin->at(j) > (MuonJetdRminCut+i*dMuonJetdRmin))
          {
            nEvts_MuonJetdRmin[i]++;
          }
      }
  }
      // //------------------------------
      // //-------------muon_jet_dRmax----------//

      for (int i = 0 ; i < nStep_MuonMuonJetdRmax ; i++)
         {
            for (unsigned int j = 0 ; j < tree_muon_jet_dRmax->size() ; j++)
               {
                  if (tree_muon_jet_dRmax->at(j) > (MuonJetdRCutmax+i*dMuonJetdRmax))
                     {
                        nEvts_MuonJetdRmax[i]++;
                     }
               }
         }
      }// End of Tree_FIlter
   }// End Global Loop

// std::cout<<"bug starts here 1"<<std::endl;
//---------------------VTX Selections Variables----
  for (int i =0 ;i<nSteps;i++)
    {
      fillHisto("hData_Vtx_MeanTrackWeight","NoSel",  thesample,MeanTWcut+i*dTW ,nSelecVtx[i]/nRecoVertex);
      fillHisto("hData_Vtx_MeanTrackWeight","TightWP",  thesample,MeanTWcut+i*dTW ,nSelecVtxStep1[i]/nRecoVertexTightWP);         
    }
// std::cout<<"bug starts here 2"<<std::endl;
for (int h = 0 ; h < nStep_10; h++)
  {
      fillHisto("hData_Vtx_nTrk10","NoSel",  thesample,ntrk10Cut+h*dntrk10 ,nEvts_ntrk10[h]/nRecoVertex);
      fillHisto("hData_Vtx_nTrk10","TightWP",  thesample,ntrk10Cut+h*dntrk10 ,nEvts_ntrk10_step1[h]/nRecoVertexTightWP);
      fillHisto("hData_Vtx_nTrk20","NoSel",  thesample,ntrk20Cut+h*dntrk20 ,nEvts_ntrk20[h]/nRecoVertex);
      fillHisto("hData_Vtx_nTrk20","LooseWP",  thesample,ntrk20Cut+h*dntrk20 ,nEvts_ntrk20_step1[h]/(nRecoVertexLooseWP+nRecoVertexTightWP));
  }
// std::cout<<"bug starts here 3"<<std::endl;
for (int h = 0 ; h < nStep_dd; h++)
  {
      fillHisto("Data_Vtx_Vtx_dd","NoSel",  thesample,VtxddCut+h*ddd ,nEvts_Vtx_dd[h]/nRecoVertex);
      fillHisto("Data_Vtx_Vtx_dd","TightWP",  thesample,VtxddCut+h*ddd ,nEvts_Vtx_dd_step1[h]/nRecoVertexTightWP);
  }
std::cout<<"bug starts here 4"<<std::endl;
//   for (int h = 0 ; h < nStep_nTrks; h++)// !!Bug!!
//   {
//       fillHisto("hData_Vtx_nTrks","NoSel",  thesample,VtxnTrksCut+h*dnTrks ,nEvts_Vtx_nTrks[h]/nRecoVertex);
//       fillHisto("hData_Vtx_nTrks","TightWP",  thesample,VtxnTrksCut+h*dnTrks,nEvts_Vtx_nTrks_step1[h]/nRecoVertexTightWP);
//   }
std::cout<<"bug starts here 5"<<std::endl;
//--------------------------------//
  for (int h = 0 ; h < nStep_nVtx; h++)
  {
      fillHisto("hData_Vtx_nVtx","NoSel",  thesample,VtxnVtxCut+h*dnVtx ,2*nEvts_Vtx_nVtx[h]/nRecoVertex);
      fillHisto("hData_Vtx_nVtx","TightWP",  thesample,VtxnVtxCut+h*dnVtx ,2*nEvts_Vtx_nVtx_step1[h]/nRecoVertexTightWP);
  }
// std::cout<<"bug starts here 6"<<std::endl;

  for (int h = 0 ; h < nStep_Chi2; h++)
  {
      fillHisto("hData_Vtx_NChi2","NoSel",  thesample,VtxChi2Cut+h*dChi2 ,nEvts_Vtx_Chi2[h]/nRecoVertex);
      fillHisto("hData_Vtx_NChi2","TightWP",  thesample,VtxChi2Cut+h*dChi2 ,nEvts_Vtx_Chi2_step1[h]/nRecoVertexTightWP);
  }
std::cout<<"bug starts here 7"<<std::endl;
//--------Vtx_Step------
  for (int h = 0 ; h < nStep_Step; h++)
  {
      fillHisto("hData_Vtx_Step","NoSel",  thesample,VtxStepCut+h*dStep ,nEvts_Vtx_Step[h]/nRecoVertex);
      fillHisto("hData_Vtx_Step","TightWP",  thesample,VtxStepCut+h*dStep,nEvts_Vtx_Step_step1[h]/nRecoVertexTightWP);
  }
std::cout<<"bug starts here 8"<<std::endl;
//--------Vtx_InvMass------
  for (int h = 0 ; h < nStep_Mass; h++)
  {
      fillHisto("hData_Vtx_Mass","NoSel",  thesample,VtxMassCut+h*dMass,nEvts_Vtx_Mass[h]/nRecoVertex);
      fillHisto("hData_Vtx_Mass","TightWP",  thesample,VtxMassCut+h*dMass ,nEvts_Vtx_Mass_step1[h]/nRecoVertexTightWP);
  }
std::cout<<"bug starts here 9"<<std::endl;
//--------Hemi_InvMass------
  for (int h = 0 ; h < nStep_HemiMass; h++)
  {
      fillHisto("hData_Hemi_Mass","NoSel",  thesample,HemiMassCut+h*dHemiMass ,nEvts_Hemi_Mass[h]/nRecoVertex);
      fillHisto("hData_Hemi_Mass","TightWP",  thesample,HemiMassCut+h*dHemiMass ,nEvts_Hemi_Mass_step1[h]/nRecoVertexTightWP);
  }
//-----------------------------------------------
std::cout<<"bug starts here 10"<<std::endl;
//--------Vtx_DCA tracks------
  for (int h = 0 ; h < nStep_DCA; h++)
  {
      fillHisto("hData_Vtx_DCA","NoSel",  thesample,VtxDCACut+h*dDCA,nEvts_Vtx_DCA[h]/nRecoVertex);
      fillHisto("hData_Vtx_DCA","TightWP",  thesample,VtxDCACut+h*dDCA ,nEvts_Vtx_DCA_step1[h]/nRecoVertexTightWP);
  }

//--------------------------------------------
// std::cout<<"bug starts here 11"<<std::endl;
//------------Evts Selection Variables------------------------
for(int h =0 ; h < nStep_jet_pt;h++)
  {
    fillHisto("hData_Evts_LeadingJetPt","NoSel",  thesample,JetPTCut+h*dpt ,nEvts_jetpt[h]/nEvts);
  }
// std::cout<<"bug starts here 12"<<std::endl;
  for(int h =0 ; h < nStep_jet2_pt;h++)
  {
    fillHisto("hData_Evts_LeadingJet2Pt","NoSel",  thesample,Jet2PTCut+h*dpt,nEvts_jetpt2[h]/(2*nEvts));
  }
// std::cout<<"bug starts here 13"<<std::endl;
for(int h =0 ; h < nStep_muon_pt;h++)
  {
    fillHisto("hData_Evts_LeadingMuonPt","NoSel",  thesample,MuonPTCut+h*dpt_muon,nEvts_muonpt[h]/nEvts);
  }
// std::cout<<"bug starts here 14"<<std::endl;
for(int h =0 ; h < nStep_muon2_pt;h++)
  {
    fillHisto("hData_Evts_LeadingMuon2Pt","NoSel",  thesample,Muon2PTCut+h*dpt_muon2 ,nEvts_muon2pt[h]/nEvts);
  }
// std::cout<<"bug starts here 15"<<std::endl;
for(int h =0 ; h < nStep_HT_pt;h++)
  {
    fillHisto("hData_Evts_HT","NoSel",  thesample,HTCut+h*dpt_HT ,nEvts_HT[h]/nEvts);
  }
// std::cout<<"bug starts here 16"<<std::endl;
for(int h =0 ; h < nStep_ST_pt;h++)
  {
    fillHisto("hData_Evts_ST","NoSel",  thesample,STCut+h*dpt_ST ,nEvts_ST[h]/nEvts);
  }
  // std::cout<<"bug starts here 17"<<std::endl;
for(int h =0 ; h < nStep_nJet;h++)
  {
    fillHisto("hData_Evts_nJets","NoSel",  thesample,nJetCut+h*dnJet_ST ,nEvts_nJet[h]/nEvts);
  }
  // std::cout<<"bug starts here 18"<<std::endl;
//--------------------------------//
  for(int h =0 ; h < nStep_nMuon;h++)
  {
    fillHisto("hData_Evts_nMuon","NoSel",  thesample,nMuonCut+h*dnMuon ,nEvts_nMuon[h]/nEvts);
  }
// std::cout<<"bug starts here 19"<<std::endl;
// //------------------------------
  for(int h =0 ; h < nStep_MuondR;h++)
  {
    fillHisto("hData_Evts_MuondR","NoSel",  thesample,MuondRCut+h*dMuondR,nEvts_MuondR[h]/nEvts);
  }

// //------------------------------
// std::cout<<"bug starts here 20"<<std::endl;
// //-------------muon_muon_dPhi----------//
  for(int h =0 ; h < nStep_MuondPhi;h++)
  {
    fillHisto("hData_Evts_MuondPhi","NoSel",  thesample,MuondPhiCut+h*dMuondPhi,nEvts_MuondPhi[h]/nEvts);
  }

// //------------------------------
// std::cout<<"bug starts here 21"<<std::endl;
// //-------------muon_muon_dEta----------//
  for(int h =0 ; h < nStep_MuondEta;h++)
  {
    fillHisto("hData_Evts_MuondEta","NoSel",  thesample,MuondEtaCut+h*dMuondEta,nEvts_MuondEta[h]/nEvts);
  }

// //------------------------------
// std::cout<<"bug starts here 22"<<std::endl;
// //-------------jet_jet_dR----------//
  for(int h =0 ; h < nStep_JetdR;h++)
  {
    fillHisto("hData_Evts_jet_jet_dR","NoSel",  thesample,JetdRCut+h*dJetdR,nEvts_JetdR[h]/nEvts);
  }

// //------------------------------

// //-------------jet_jet_dPhi----------//
// std::cout<<"bug starts here 23"<<std::endl;
  for(int h =0 ; h < nStep_JetdPhi;h++)
  {
    fillHisto("hData_Evts_jet_jet_dPhi","NoSel",  thesample,JetdPhiCut+h*dJetdPhi,nEvts_JetdPhi[h]/nEvts);
  }
// //------------------------------
// std::cout<<"bug starts here 24"<<std::endl;
// //-------------jet_jet_dEta----------//
  for(int h =0 ; h < nStep_JetdEta;h++)
  {
    fillHisto("hData_Evts_jet_jet_dEta","NoSel",  thesample,JetdEtaCut+h*dJetdEta ,nEvts_JetdEta[h]/nEvts);
  }

// //------------------------------
// std::cout<<"bug starts here 25"<<std::endl;
// //-------------muon_jet_dRmin----------//
  for(int h =0 ; h < nStep_MuonJetdRmin;h++)
  {
    fillHisto("hData_Evts_muon_jet_dRmin","NoSel",  thesample,MuonJetdRminCut+h*dMuonJetdRmin,nEvts_MuonJetdRmin[h]/(2*nEvts));
  }
// std::cout<<"bug starts here 26"<<std::endl;
// //------------------------------
// //-------------muon_jet_dRmax----------//
  for(int h =0 ; h < nStep_MuonMuonJetdRmax;h++)
  {
    fillHisto("hData_Evts_muon_jet_dRmax","NoSel",  thesample, MuonJetdRCutmax+h*dMuonJetdRmax ,nEvts_MuonJetdRmax[h]/(2*nEvts));
  }

//------------------------------------------------

   float MeanD = MeanDistance/nRecoVertex;
   float MeanDTightWP = MeanDistanceTightWP/nRecoVertexTightWP;
   float MeanDLooseWP = MeanDistanceLooseWP/nRecoVertexLooseWP;
   float NRecoVertexEff = nRecoVertex/TotalnVertex;

   // Put in the file when the checks have been done
   // std::cout<<"Checks "<<std::endl;

   // std::cout<<"End of checks "<<std::endl;
   //---------------------------------------------
   float EvtSize = nentries;
   // Output Postscript
   ofs<<" //----------------- Efficacity for the Ntuple : " << thesample << " ---------------\\ "<<std::endl;
   ofs<<" //-Tracks selection efficiency by the Tight WP : " << TightTrks_Eff << " ---------------\\ "<<std::endl;
   ofs<<" //-Tracks selection efficiency by the Loose WP : " << LooseTrks_Eff << " ---------------\\ "<<std::endl;
   ofs<<"|| ----------------------------------------------------------------------------------"<<std::endl;
   ofs<<"|| Nvertex to be reco considering two vtx per evt :"<<nentries<<std::endl;
   ofs<<"|| NVertex To be reco : "<<TotalnVertex<<" which means "<<TotalnVertex/(2*nentries)<<" efficiency of the filter"<<std::endl;
   ofs<<"MeanDistance :"<<MeanD<<" cm "<<std::endl;
   ofs<<"MeanDistanceTightWP :"<<MeanDTightWP<<" cm "<<std::endl;
   ofs<<"MeanDistanceLooseWP :"<<MeanDLooseWP<<" cm "<<std::endl;
   ofs<<" Total NRecoVertexEff (vertex reco over vertex to be reconstructed): "<<NRecoVertexEff<<std::endl;
   ofs<<" NRecoVertexEff TightWp: "<<100*nRecoVertexTightWP/nRecoVertex<<" per cent of "<<NRecoVertexEff<<std::endl;
   ofs<<"  NRecoVertexEff  LooseWP: "<<100*nRecoVertexLooseWP/nRecoVertex<<" per cent of "<<NRecoVertexEff<<std::endl;
   ofs<<" //----------------- End of Table ---------------\\ "<<std::endl; 
   ofs.close();

   theoutputfile->Write();
   //deleteHisto();
   theoutputfile->Close();
   delete theoutputfile;

}



void TreeReader::initializeHisto(TString sample, bool isfirstset){


  cout << "#####################################" << endl;
  cout << "#####################################" << endl;
  cout << " initialize histograms               " << endl;
  cout << "#####################################" << endl;
  cout << "#####################################" << endl;
  

   if(isfirstset){
      numb_histo = 0;
      TH1F * first_emptyHisto = new TH1F("first_emptyHisto", "first_emptyHisto", 100, 0, 1000);
      histo_list_.push_back(first_emptyHisto);
  
      numb_histo++;
   }

   // ---------------------------- Filter (Online+Offline Selection)--------------
// Dimuon : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
// OR
// Single Muon : HLT_IsoMu24_v
// Offline selection of muons : Two Loose ID and PFIsoLoose muons required
// |dxy | < 0.1 cm & |dz | < 0.2 cm (prompt)
// pt1 > 25GeV & pt2 > 10GeV
// Mμμ > 10 GeV (remove low-resonances)

   addHisto("Filter", "Offline+Online", sample.Data(), 2,0,2);

   //----------------------------- Muons -----------------------------------------

   addHisto("RecoMuo_pT", "noSel", sample.Data(),  100, 0, 100);

   // ---------------------------- Electrons -------------------------------------

   // --------------------------  V0 Candidates ----------------------------------

   addHisto("hData_reco_K0_mass","noSel", sample.Data(),101,0.42,0.58);
   addHisto("hData_reco_K0_r","noSel",    sample.Data(),200,0,100);
   addHisto("hData_reco_K0_z","noSel",    sample.Data(),401,-200,200);
   addHisto("hData_reco_L0_mass","noSel", sample.Data(),101,1.06,1.18);
   addHisto("hData_reco_L0_r","noSel",    sample.Data(),200,0,100);
   addHisto("hData_reco_L0_z","noSel",    sample.Data(),401,-200,200);

   // ------------------------- Sec. Interactions --------------------------------

   addHisto("hData_reco_SecInt_mass","noSel", sample.Data(),100,0,10);
   addHisto("hData_reco_SecInt_x","noSel",    sample.Data(), 201,-100,100);
   addHisto("hData_reco_SecInt_y","noSel",    sample.Data(), 201,-100,100);
   addHisto("hData_reco_SecInt_r","noSel",    sample.Data(),200,0,100);
   addHisto("hData_reco_SecInt_z","noSel",    sample.Data(),401,-200,200);

   addHisto("hData_reco_SecInt_mass","Selec", sample.Data(),100,0,10);
   addHisto("hData_reco_SecInt_x","Selec",    sample.Data(),  201,-100,100);
   addHisto("hData_reco_SecInt_y","Selec",    sample.Data(),  201,-100,100);
   addHisto("hData_reco_SecInt_r","Selec",    sample.Data(),200,0,100);
   addHisto("hData_reco_SecInt_z","Selec",    sample.Data(),401,-200,200);

   addHisto("hData_reco_SecInt_mass","TrackerMatched", sample.Data(),100,0,10);
   addHisto("hData_reco_SecInt_x","TrackerMatched",    sample.Data(),  201,-100,100);
   addHisto("hData_reco_SecInt_y","TrackerMatched",    sample.Data(),  201,-100,100);
   addHisto("hData_reco_SecInt_r","TrackerMatched",    sample.Data(),200,0,100);
   addHisto("hData_reco_SecInt_z","TrackerMatched",    sample.Data(),401,-200,200);

   // ------------------------ Tracks -------------------------------------------

   addHisto("hData_MVAVal","noSel",                sample.Data(),101,-1,1);

   // ------------------------ Axes --------------------------------------------

  addHisto("hData_dR_RecoReco","noSel",            sample.Data(),201,0,6);

   //---------------------------Signal-Vertices---------------------------------------

   addHisto("hData_Hemi_Vtx_NChi2","noSel",        sample.Data(),70,-20.,50.);
   addHisto("hData_Hemi_Vtx_nTrks","noSel",        sample.Data(),40,0.5,40.5);

   addHisto("hData_Hemi_Vtx_NChi2","RecoVtx",      sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_nTrks","RecoVtx",      sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_dist","RecoVtx",       sample.Data(),50,0.,100.);
   addHisto("hData_Hemi_Vtx_step","RecoVtx",       sample.Data(),5,0,5);
   addHisto("hData_Hemi_Vtx_eta","RecoVtx",        sample.Data(),260,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","RecoVtx",          sample.Data(),200,0,100);
   addHisto("hData_Hemi_Vtx_z","RecoVtx",          sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","RecoVtx",     sample.Data(),101,-1,1);

   addHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx",  sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx",  sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_dist","GoodRecoVtx",   sample.Data(),50,0.,100.);
   addHisto("hData_Hemi_Vtx_step","GoodRecoVtx",   sample.Data(),5,0,5);
   addHisto("hData_Hemi_Vtx_eta","GoodRecoVtx",      sample.Data(),260,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","GoodRecoVtx",      sample.Data(),200,0,100);
   addHisto("hData_Hemi_Vtx_z","GoodRecoVtx",      sample.Data() ,401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx", sample.Data(),101,-1,1);  

   addHisto("hData_Hemi_Vtx_nTrks","RecoVtx_TightWP",    sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_NChi2","RecoVtx_TightWP",    sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_dist","RecoVtx_TightWP",     sample.Data(),50,0.,100.);
   addHisto("hData_Hemi_Vtx_eta","RecoVtx_TightWP",      sample.Data(),260,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","RecoVtx_TightWP",        sample.Data(),200,0,100);
   addHisto("hData_Hemi_Vtx_z","RecoVtx_TightWP",        sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","RecoVtx_TightWP",   sample.Data(),101,-1,1);

   addHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx_TightWP",   sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx_TightWP",   sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_dist","GoodRecoVtx_TightWP",    sample.Data(),50,0.,100.);
   addHisto("hData_Hemi_Vtx_eta","GoodRecoVtx_TightWP",    sample.Data(),260,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","GoodRecoVtx_TightWP",       sample.Data(),200,0,100);
   addHisto("hData_Hemi_Vtx_z","GoodRecoVtx_TightWP",       sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx_TightWP",  sample.Data(),101,-1,1);

   addHisto("hData_Hemi_Vtx_nTrks","RecoVtx_LooseWP",       sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_NChi2","RecoVtx_LooseWP",       sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_dist","RecoVtx_LooseWP",        sample.Data(),50,0.,100.);
   addHisto("hData_Hemi_Vtx_eta","RecoVtx_LooseWP",         sample.Data(),260,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","RecoVtx_LooseWP",           sample.Data(),200,0,100);
   addHisto("hData_Hemi_Vtx_z","RecoVtx_LooseWP",           sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","RecoVtx_LooseWP",      sample.Data(),101,-1,1);

   addHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx_LooseWP",   sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx_LooseWP",   sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_dist","GoodRecoVtx_LooseWP",    sample.Data(),50,0.,100.);
   addHisto("hData_Hemi_Vtx_eta","GoodRecoVtx_LooseWP",    sample.Data(),260,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","GoodRecoVtx_LooseWP",       sample.Data(),200,0,100);
   addHisto("hData_Hemi_Vtx_z","GoodRecoVtx_LooseWP",       sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx_LooseWP",  sample.Data(),101,-1,1);  


   //-------  Evt Level Variable ------------//

   addHisto("hData_Evts_LeadingJetPt","NoSel",             sample.Data(),100,30,1030);
   addHisto("hData_Evts_LeadingJet2Pt","NoSel",            sample.Data(),100,20,1020);
   addHisto("hData_Evts_LeadingMuonPt","NoSel",            sample.Data(),50,20,520);
   addHisto("hData_Evts_LeadingMuon2Pt","NoSel",           sample.Data(),100,10,510);
   addHisto("hData_Evts_HT","NoSel",                       sample.Data(),200,100,2100);
   addHisto("hData_Evts_ST","NoSel",                       sample.Data(),50,30,530);
   addHisto("hData_Evts_nJets","NoSel",                    sample.Data(),20,0,19);
   addHisto("hData_Evts_nMuon","NoSel",                    sample.Data(),20,0,20);
   addHisto("hData_Evts_MuondR","NoSel",                   sample.Data(),40,0,5);
   addHisto("hData_Evts_MuondPhi","NoSel",                 sample.Data(),28,0,3.5);
   addHisto("hData_Evts_MuondEta","NoSel",                 sample.Data(),32,0,4);
   addHisto("hData_Evts_jet_jet_dR","NoSel",               sample.Data(),40,0,5);
   addHisto("hData_Evts_jet_jet_dPhi","NoSel",             sample.Data(),28,0,3.5);
   addHisto("hData_Evts_jet_jet_dEta","NoSel",             sample.Data(),40,0,5);
   addHisto("hData_Evts_muon_jet_dRmin","NoSel",           sample.Data(),40,0,5);
   addHisto("hData_Evts_muon_jet_dRmax","NoSel",           sample.Data(),48,0,5);

   //-------- Vtx Level Variables -----------//

   addHisto("hData_Vtx_MeanTrackWeight","NoSel",           sample.Data(),1001,1.5,2.5);
   addHisto("hData_Vtx_MeanTrackWeight","TightWP",         sample.Data(),1001,1.5,2.5);         
   addHisto("hData_Vtx_nTrk10","NoSel",                    sample.Data(),50,0,50);
   addHisto("hData_Vtx_nTrk10","TightWP",                  sample.Data(),50,0,50);
   addHisto("hData_Vtx_nTrk20","NoSel",                    sample.Data(),50,0,50);
   addHisto("hData_Vtx_nTrk20","LooseWP",                  sample.Data(),50,0,50);
   addHisto("Data_Vtx_Vtx_dd","NoSel",                     sample.Data(),250,0,250);
   addHisto("Data_Vtx_Vtx_dd","TightWP",                   sample.Data(),250,0,250);
   addHisto("hData_Vtx_nTrks","NoSel",                     sample.Data(),50,2,52);
   addHisto("hData_Vtx_nTrks","TightWP",                   sample.Data(),50,2,52);
   addHisto("hData_Vtx_nVtx","NoSel",                      sample.Data(),4,0,4);
   addHisto("hData_Vtx_nVtx","TightWP",                    sample.Data(),4,0,4);
   addHisto("hData_Vtx_NChi2","NoSel",                     sample.Data(),40,0,10);
   addHisto("hData_Vtx_NChi2","TightWP",                   sample.Data(),40,0,10);
   addHisto("hData_Vtx_Step","NoSel",                      sample.Data(),5,0,5);
   addHisto("hData_Vtx_Step","TightWP",                    sample.Data(),5,0,5);
   addHisto("hData_Vtx_Mass","NoSel",                      sample.Data(),1400,0,14000);
   addHisto("hData_Vtx_Mass","TightWP",                    sample.Data(),1400,0,14000);
   addHisto("hData_Hemi_Mass","NoSel",                     sample.Data(),140,0,1400);
   addHisto("hData_Hemi_Mass","TightWP",                   sample.Data(),140,0,1400);
   addHisto("hData_Vtx_DCA","NoSel",                       sample.Data(),100,0,100);
   addHisto("hData_Vtx_DCA","TightWP",                     sample.Data(),100,0,100);

}


//-------------------------------------------------------------
//instantiate and add
//first parameter is the variable name,
// second parameter is the selection step (like "afterleptsel")
//third parameter is the sample name (like "Z)
//others are TH1F binning
//creates one histograms per channel
//-------------------------------------------------------------
void TreeReader::addHisto(TString var, TString selstep, TString sample, int nbins, float min, float max){
 
  TString name =  var+"_"+selstep+"__"+sample;
  TH1F * thehisto = new TH1F(name,name,nbins,min,max);
  thehisto->Sumw2();
  thehisto->SetOption("HIST");

  histo_list_.push_back(thehisto);
  histo_map_[name.Data()] = numb_histo;
  


  numb_histo++;
  
}


//-------------------------------------------------------------
//fill histograms
//first parameter is the channel,
//second parameter is the variable name,
//third parameter is the selection step (like "afterleptsel")
//forths parameter is the sample name (like "Z)
//others are value and weight
//-------------------------------------------------------------
void TreeReader::fillHisto( TString var, TString selstep, TString sample, float val, float weight){
  TString name = var+"_"+selstep+"__"+sample;


  if(histo_map_[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histograms " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }else  histo_list_[histo_map_[name.Data()]]->Fill(val, weight);
  
}


void TreeReader::deleteHisto(){
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
