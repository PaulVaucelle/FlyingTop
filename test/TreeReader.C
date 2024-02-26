#define TreeReader_cxx
#include "TreeReader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <ctime> 
using namespace std;

void TreeReader::Loop(TString sample, bool Signal)//void TreeReader::Loop(TString sample)
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
// TString thesample = sample.substr();
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

   float TightWP = 0.85;// for tracks => everytime Tight is mentioned, it is a refrence to this value.
   float LooseWP = 0.; // for tracks => same for loose

   float EVTSWP = 0;
   float VTXWP = 0.5;

   // -------- List of variables ---------------//

   float n_TightTrks = 0;
   float n_LooseTrks = 0;
   float n_TotalTrks = 0;
   float TightTrks_Eff = 0.;
   float LooseTrks_Eff = 0.;

   float nRecoVertex = 0;
   float nRecoVertex2Trks = 0;
   float MeanDistance = 0;
   float nRecoVertexTightWP = 0;
   float MeanDistanceTightWP = 0;
   float nRecoVertexLooseWP = 0;
   float MeanDistanceLooseWP = 0;
   float TotalnVertex = 0;
   float nEvts = 0;
   float nFilterEvt = 0;
   float nSelecBDTVtx = 0;
   float nSelecBDTVtxTight = 0;
   float nSelecBDTVtxTight_step12 = 0;


   float nEvts_w2Vtx = 0;
   float nEvts_w2TightVtx = 0;
   float nEvts_w2TightBDTVtx = 0;

  float nEvts_w2TightBDTVtx_step12 = 0;
  float nEvts_w1TightBDTVtx_step12 = 0;

  float nEvts_w1TightVtx = 0 ;
  float nEvts_w1TightBDTVtx = 0;
// Evts selectipn BDT-------------------

float MeanBDTcut = -1;
float dBDT = 0.02;
int nSteps_BDT = (1-MeanBDTcut)/dBDT;
float nSelecEvts[100] = {0};

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
int nStep_nTrks = (52-VtxnTrksCut)/dnTrks;
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

// //-------------HemiCombinedLeptonMass----------//
float Hemi_CombinedHemiLeptonMassCut = 0;
float dSmuonMass = 10;
int nStep_SmuonMass = (1400-Hemi_CombinedHemiLeptonMassCut)/dSmuonMass;//1000-20/dpt
float nEvts_Hemi_CombinedHemiLeptonMass[200] = {0};
//------------------------------


   // -------------End of list of variables ----------//

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

// Normalisation factor (XS)
  float NormFactor = 1. ;
  float XS = 1;

if (thesample.Contains("DYJetsToLL_M10to50"))                     { XS = 15910.0;   }
if (thesample.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays")) { XS = 10.8707;   }
if (thesample.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     { XS = 10.8908;   }
if (thesample.Contains("TTJets_DiLept"))                          { XS = 53.07;     }
if (thesample.Contains("WWTo2L2Nu_MLL_200To600"))                 { XS = 11.09;     }
if (thesample.Contains("WWTo2L2Nu_MLL_600To1200"))                { XS = 11.09;     }
if (thesample.Contains("WWTo2L2Nu"))                              { XS = 11.09;     }
if (thesample.Contains("WZTo2Q2L_mllmin4p0"))                     { XS = 6.535;     }
if (thesample.Contains("ZZTo2Q2L_mllmin4p0"))                     { XS = 3.676;     }
if (thesample.Contains("TTTo2L2Nu"))                              { XS = 88.3;      }
if (thesample.Contains("DYJetsToLL_M50"))                         { XS = 5379;      }
if (thesample.Contains("ttWJetsToLNu_5f_EWK"))                    { XS = 0.290;     } // not found on XSDB, no file on tier2...approximation
      //Took 0.868 pb (CMS-TOP-21-011)
     // as a starting point and then divided by 3 (lepton universality)
if (thesample.Contains("TTZToLL_5f"))                             { XS = 0.05188;   }//Not found on XSDB => used ana.py macro : 5.188e-02 +- 2.437e-04 pb
if (thesample.Contains("TTToHadronic"))                           { XS = 378.9;     }
if (thesample.Contains("TTWW"))                                   { XS = 0.006992;  }//found on XSDB
if (thesample.Contains("TTToSemiLeptonic")  )                      { XS = 365.34;    }


if (thesample.Contains("smu200"))
    {        
      XS = 0.01;//found on XSDB
    }
    if (thesample.Contains("smu250"))
    {        
      XS = 0.0045;//found on XSDB
    }
    if (thesample.Contains("smu300"))
    {        
      XS = 0.002;//found on XSDB
    }

    if (thesample.Contains("smu350"))
    {        
      XS = 0.001;//found on XSDB
    }

    if (thesample.Contains("smu400"))
    {        
      XS = 0.0006;//found on XSDB
    }

    if (thesample.Contains("smu450"))
    {        
      XS = 0.0004;//found on XSDB
    }

    if (thesample.Contains("smu500"))
    {        
      XS = 0.00025;//found on XSDB
    }
    if (thesample == "smu300_250_240")
    {        
      XS = 0.002;//found on XSDB
    }

  float lumiRun2 = 136000 ; //pb-1
  float Nevent = nentries;
  // if (nentries>1000000){nentries = 1000000;}
  
  NormFactor =  XS/nentries;

  cout<< "Line : "  << __LINE__ << " " << nentries << endl;
   Long64_t nbytes = 0, nb = 0;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      auto tim = std::chrono::system_clock::now();
      std::time_t start = std::chrono::system_clock::to_time_t(tim);
      if (jentry==0.1*nentries) {std::cout<<"10/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.2*nentries) {std::cout<<"20/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.3*nentries) {std::cout<<"30/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.4*nentries) {std::cout<<"40/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.5*nentries) {std::cout<<"50/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.6*nentries) {std::cout<<"60/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.7*nentries) {std::cout<<"70/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.8*nentries) {std::cout<<"80/100 :"<<std::ctime(&start)<<std::endl;}
      if (jentry==0.9*nentries) {std::cout<<"90/100 :"<<std::ctime(&start)<<std::endl;}
      
            // std::cout<<"event : "<<jentry<<std::endl;
      fillHisto("Filter", "Offline+Online", thesample, tree_Filter, 1);





   if (tree_Filter)
      {
        if (tree_njetNOmu == 0) continue;
      nFilterEvt++;

    // Daniel ----------------------//
    
    fillHisto("hData_njet","NoSel",  thesample,tree_njet ,1.);
    fillHisto("hData_njetNOmu","NoSel",  thesample,tree_njetNOmu ,1.);
    fillHisto("hData_Mmumu","NoSel",  thesample,tree_Mmumu ,1.);
    fillHisto("hData_BDTevt","NoSel",  thesample,tree_Evts_MVAval ,1.);
  
    int nHemi = tree_Hemi_njet->size();
    for (int i=0; i<nHemi; i++) 
      {   // Loop on hemispheres
        int njet = tree_Hemi_njet->at(i);
        fillHisto("hData_Hemi_njet","NoSel",  thesample,njet ,1.);

        if ( njet >= 2 ) fillHisto("hData_Hemi_mass","NoSel",  thesample,tree_Hemi_mass->at(i) ,1.);
        if ( njet >= 2 ) fillHisto("hData_Hemi_massTot","NoSel",  thesample,tree_HemiMu_mass->at(i) ,1.);
        if ( tree_Hemi_Vtx_NChi2->at(i) > 0. && tree_Hemi_Vtx_NChi2->at(i) < 10. ) 
          {
            fillHisto("hData_Hemi_Vtx_njet","NoSel",  thesample,njet ,1.);
            if ( njet >= 2 ) fillHisto("hData_Hemi_Vtx_mass","NoSel",  thesample,tree_Hemi_mass->at(i) ,1.);
            if ( njet >= 2 ) fillHisto("hData_Hemi_Vtx_massTot","NoSel",  thesample,tree_HemiMu_mass->at(i) ,1.);
          }
        if ( i == 0 ) 
          {
            fillHisto("hData_Hemi1_njet","NoSel",  thesample,njet ,1.);
            if ( njet >= 2 ) fillHisto("hData_Hemi1_mass","NoSel",  thesample,tree_Hemi_mass->at(i) ,1.);
            if ( njet >= 2 ) fillHisto("hData_Hemi1_massTot","NoSel",  thesample,tree_HemiMu_mass->at(i) ,1.);
            if ( tree_Hemi_Vtx_NChi2->at(i) > 0. && tree_Hemi_Vtx_NChi2->at(i) < 10. ) 
              {
                fillHisto("hData_Hemi1_Vtx_dist","NoSel",  thesample,tree_Hemi_Vtx_dist->at(i) ,1.);
                fillHisto("hData_Hemi1_Vtx_BDTvtx","NoSel",  thesample,tree_Hemi_Vtx_MVAval->at(i) ,1.);
                fillHisto("hData_Hemi1_Vtx_njet","NoSel",  thesample,njet ,1.);
                if ( njet >= 2 ) fillHisto("hData_Hemi1_Vtx_mass","NoSel",  thesample,tree_Hemi_mass->at(i) ,1.);
                if ( njet >= 2 ) fillHisto("hData_Hemi1_Vtx_massTot","NoSel",  thesample,tree_HemiMu_mass->at(i) ,1.);
              }
          }
        else 
          {
            fillHisto("hData_Hemi2_njet","NoSel",  thesample,njet ,1.);
            if ( njet >= 2 ) fillHisto("hData_Hemi2_mass","NoSel",  thesample,tree_Hemi_mass->at(i) ,1.);
            if ( njet >= 2 ) fillHisto("hData_Hemi2_massTot","NoSel",  thesample,tree_HemiMu_mass->at(i) ,1.);
            if ( tree_Hemi_Vtx_NChi2->at(i) > 0. && tree_Hemi_Vtx_NChi2->at(i) < 10. ) 
              {
                fillHisto("hData_Hemi2_Vtx_dist","NoSel",  thesample,tree_Hemi_Vtx_dist->at(i) ,1.);
                fillHisto("hData_Hemi2_Vtx_BDTvtx","NoSel",  thesample,tree_Hemi_Vtx_MVAval->at(i) ,1.);
                fillHisto("hData_Hemi2_Vtx_njet","NoSel",  thesample,njet ,1.);
                if ( njet >= 2 ) fillHisto("hData_Hemi2_Vtx_mass","NoSel",  thesample,tree_Hemi_mass->at(i) ,1.);
                if ( njet >= 2 ) fillHisto("hData_Hemi2_Vtx_massTot","NoSel",  thesample,tree_HemiMu_mass->at(i) ,1.);
              }
          }
      } // end loop on hemispheres

   if ( (tree_Hemi_Vtx_NChi2->at(0) > 0. && tree_Hemi_Vtx_NChi2->at(0) < 10.)
     || (tree_Hemi_Vtx_NChi2->at(1) > 0. && tree_Hemi_Vtx_NChi2->at(1) < 10.) ) 
      {
          fillHisto("hData_vtx1_Mmumu","NoSel",  thesample,tree_Mmumu ,1.);
          if ( tree_Evts_MVAval > 0. ) fillHisto("hData_BDTevt1_Mmumu","NoSel",  thesample,tree_Mmumu ,1.);
          if ( tree_Hemi_Vtx_MVAval->at(0) > 0.5 || tree_Hemi_Vtx_MVAval->at(1) > 0.5 ) 
            fillHisto("hData_BDTvtx1_Mmumu","NoSel",  thesample,tree_Mmumu ,1.);
          if ( tree_Evts_MVAval > 0. && 
        ( tree_Hemi_Vtx_MVAval->at(0) > 0.5 || tree_Hemi_Vtx_MVAval->at(1) > 0.5 ) ) 
            fillHisto("hData_BDTevtvtx1_Mmumu","NoSel",  thesample,tree_Mmumu ,1.);
      }

   if ( (tree_Hemi_Vtx_NChi2->at(0) > 0. && tree_Hemi_Vtx_NChi2->at(0) < 10.)
     && (tree_Hemi_Vtx_NChi2->at(1) > 0. && tree_Hemi_Vtx_NChi2->at(1) < 10.) ) {
     fillHisto("hData_vtx2_Mmumu","NoSel",  thesample,tree_Mmumu ,1.);
     if ( tree_Evts_MVAval > 0. )  fillHisto("hData_BDTevt2_Mmumu","NoSel",  thesample,tree_Mmumu ,1.);
     if ( tree_Hemi_Vtx_MVAval->at(0) > 0.5 && tree_Hemi_Vtx_MVAval->at(1) > 0.5 ) 
     fillHisto("hData_BDTvtx2_Mmumu","NoSel",  thesample,tree_Mmumu ,1.);
     if ( tree_Evts_MVAval > 0. && 
        ( tree_Hemi_Vtx_MVAval->at(0) > 0.5 && tree_Hemi_Vtx_MVAval->at(1) > 0.5 ) ) 
     fillHisto("hData_BDTevtvtx2_Mmumu","NoSel",  thesample,tree_Mmumu ,1.);
   }
    //end of  Daniel's part ----------------------//


      fillHisto("DiMuon_Mass","noSel", thesample,  tree_Mmumu, XS/nentries);
      if (abs(tree_Evts_MVAval)<2) fillHisto("Evts_MVAVal","noSel", thesample , tree_Evts_MVAval,1);
      

   if ( Signal)
    {        
      int Msmu = tree_smu_mass;
      int Mneu = tree_neu_mass;
      if ( Msmu > 180 && Msmu < 220 )      Msmu = 200;
      else if ( Msmu > 230 && Msmu < 270 ) Msmu = 250;
      else if ( Msmu > 280 && Msmu < 320 ) Msmu = 300;
      else if ( Msmu > 380 && Msmu < 420 ) Msmu = 400;
      else if ( Msmu > 480 && Msmu < 520 ) Msmu = 500;
      // else cout << " !!! smu mass out of range !!! " << Msmu;
      if ( Mneu > 170 && Mneu < 190 )      Mneu = 180;
      else if ( Mneu > 190 && Mneu < 210 ) Mneu = 200;
      else if ( Mneu > 240 && Mneu < 260 ) Mneu = 250;
      else if ( Mneu > 270 && Mneu < 290 ) Mneu = 280;
      else if ( Mneu > 290 && Mneu < 310 ) Mneu = 300;
      else if ( Mneu > 340 && Mneu < 360 ) Mneu = 350;
      else if ( Mneu > 370 && Mneu < 390 ) Mneu = 380;
      else if ( Mneu > 390 && Mneu < 410 ) Mneu = 400;
      else if ( Mneu > 440 && Mneu < 460 ) Mneu = 450;
      else if ( Mneu > 470 && Mneu < 490 ) Mneu = 480;
      // else cout << " !!! neu mass out of range !!! " << Mneu;


      
      if ((Msmu-Mneu) == 20){fillHisto("hSim_EVTBDT","dm20",  thesample,tree_Evts_MVAval ,1.);}
      if ((Msmu-Mneu) == 50){fillHisto("hSim_EVTBDT","dm50",  thesample,tree_Evts_MVAval ,1.);}
      if ((Msmu-Mneu) == 100){fillHisto("hSim_EVTBDT","dm100",  thesample,tree_Evts_MVAval ,1.);}
      if ((Msmu-Mneu) == 150){fillHisto("hSim_EVTBDT","dm150",  thesample,tree_Evts_MVAval ,1.);}
      if ((Msmu-Mneu) == 200){fillHisto("hSim_EVTBDT","dm200",  thesample,tree_Evts_MVAval ,1.);}
      if ((Msmu-Mneu) == 250){fillHisto("hSim_EVTBDT","dm250",  thesample,tree_Evts_MVAval ,1.);}
      if ((Msmu-Mneu) == 300){fillHisto("hSim_EVTBDT","dm300",  thesample,tree_Evts_MVAval ,1.);}
      if ((Msmu-Mneu) == 320){fillHisto("hSim_EVTBDT","dm320",  thesample,tree_Evts_MVAval ,1.);}
     

      int ngenpart =  tree_genParticle_pt->size();
      for (int i=0; i<ngenpart; i++)    // Loop on GenParticle
        {
          float pdgId = tree_genParticle_pdgId->at(i); 
          float mother_pdgId = tree_genParticle_mother_pdgId->at(i); 
          float ct0 = tree_genParticle_ct0->at(i);
          // top quark from neutralino
          if ( abs(pdgId) == 6 && abs(mother_pdgId) == 1000023 ) 
            {            
            

              float tree_jet_pt2 = 0;
              if (tree_jet_leadingpt2->size()>0){tree_jet_pt2 = tree_jet_leadingpt2->at(0);}
              if ((Msmu-Mneu) == 20){fillHisto("hSim_ct0","dm20",  thesample,ct0 ,1.);fillHisto("hBDT_MET","EVT_dm20", thesample,tree_PFMet_et,1.);fillHisto("hBDT_TRACK_SIZE","EVT_dm20", thesample,tree_TRACK_SIZE,1.);fillHisto("hBDT_muon_leadingpt","EVT_dm20", thesample,tree_lepton_leadingpt->at(0),1.);fillHisto("hBDT_muon_leadingpt2","EVT_dm20", thesample,tree_lepton_leadingpt2->at(0),1.);fillHisto("hBDT_jet_leadingpt","EVT_dm20", thesample,tree_jet_leadingpt->at(0),1.); fillHisto("hBDT_jet_leadingpt2","EVT_dm20", thesample,tree_jet_pt2 ,1.);fillHisto("hBDT_muon_muon_dR","EVT_dm20", thesample,tree_lepton_lepton_dR->at(0),1.);fillHisto("hBDT_muon_muon_dPhi","EVT_dm20", thesample,tree_lepton_lepton_dPhi->at(0),1.);fillHisto("hBDT_muon_muon_dEta","EVT_dm20", thesample,tree_lepton_lepton_dEta->at(0),1.);fillHisto("hBDT_jet_jet_dR","EVT_dm20", thesample,tree_jet_jet_dR->at(0),1.);fillHisto("hBDT_jet_jet_dPhi","EVT_dm20", thesample,tree_jet_jet_dPhi->at(0),1.);fillHisto("hBDT_jet_jet_dEta","EVT_dm20", thesample,tree_jet_jet_dEta->at(0),1.);fillHisto("hBDT_muon_jet_dRmin","EVT_dm20", thesample,tree_muon_jet_dRmin->at(0),1.);fillHisto("hBDT_muon_jet_dRmax","EVT_dm20", thesample,tree_muon_jet_dRmax->at(0),1.);fillHisto("hBDT_HT","EVT_dm20", thesample,tree_HT,1.);fillHisto("hBDT_ST","EVT_dm20", thesample,tree_LT,1.);fillHisto("hBDT_njet","EVT_dm20", thesample,tree_njet,1.);fillHisto("hBDT_muon_nmu","EVT_dm20", thesample,tree_nmu,1.);fillHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm20", thesample,tree_Hemi_LooseBTag_axes->at(0),1.); fillHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm20", thesample,tree_Hemi_MediumBTag_axes->at(0),1.);fillHisto("hBDT_Hemi_TightBTag_axes","EVT_dm20", thesample,tree_Hemi_TightBTag_axes->at(0),1.);}
              
              if ((Msmu-Mneu) == 50){fillHisto("hSim_ct0","dm50",  thesample,ct0 ,1.);fillHisto("hBDT_MET","EVT_dm50", thesample,tree_PFMet_et,1.);fillHisto("hBDT_TRACK_SIZE","EVT_dm50", thesample,tree_TRACK_SIZE,1.);fillHisto("hBDT_muon_leadingpt","EVT_dm50", thesample,tree_lepton_leadingpt->at(0),1.); fillHisto("hBDT_muon_leadingpt2","EVT_dm50", thesample,tree_lepton_leadingpt2->at(0),1.); fillHisto("hBDT_jet_leadingpt","EVT_dm50", thesample,tree_jet_leadingpt->at(0),1.);fillHisto("hBDT_jet_leadingpt2","EVT_dm50", thesample,tree_jet_pt2 ,1.);fillHisto("hBDT_muon_muon_dR","EVT_dm50", thesample,tree_lepton_lepton_dR->at(0),1.);fillHisto("hBDT_muon_muon_dPhi","EVT_dm50", thesample,tree_lepton_lepton_dPhi->at(0),1.);fillHisto("hBDT_muon_muon_dEta","EVT_dm50", thesample,tree_lepton_lepton_dEta->at(0),1.);fillHisto("hBDT_jet_jet_dR","EVT_dm50", thesample,tree_jet_jet_dR->at(0),1.);fillHisto("hBDT_jet_jet_dPhi","EVT_dm50", thesample,tree_jet_jet_dPhi->at(0),1.);fillHisto("hBDT_jet_jet_dEta","EVT_dm50", thesample,tree_jet_jet_dEta->at(0),1.);fillHisto("hBDT_muon_jet_dRmin","EVT_dm50", thesample,tree_muon_jet_dRmin->at(0),1.);fillHisto("hBDT_muon_jet_dRmax","EVT_dm50", thesample,tree_muon_jet_dRmax->at(0),1.);fillHisto("hBDT_HT","EVT_dm50", thesample,tree_HT,1.);fillHisto("hBDT_ST","EVT_dm50", thesample,tree_LT,1.);fillHisto("hBDT_njet","EVT_dm50", thesample,tree_njet,1.);fillHisto("hBDT_muon_nmu","EVT_dm50", thesample,tree_nmu,1.);fillHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm50", thesample,tree_Hemi_LooseBTag_axes->at(0),1.); fillHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm50", thesample,tree_Hemi_MediumBTag_axes->at(0),1.);fillHisto("hBDT_Hemi_TightBTag_axes","EVT_dm50", thesample,tree_Hemi_TightBTag_axes->at(0),1.);}
              if ((Msmu-Mneu) == 100){fillHisto("hSim_ct0","dm100",  thesample,ct0 ,1.);fillHisto("hBDT_MET","EVT_dm100", thesample,tree_PFMet_et,1.);fillHisto("hBDT_TRACK_SIZE","EVT_dm100", thesample,tree_TRACK_SIZE,1.);fillHisto("hBDT_muon_leadingpt","EVT_dm100", thesample,tree_lepton_leadingpt->at(0),1.); fillHisto("hBDT_muon_leadingpt2","EVT_dm100", thesample,tree_lepton_leadingpt2->at(0),1.); fillHisto("hBDT_jet_leadingpt","EVT_dm100", thesample,tree_jet_leadingpt->at(0),1.);fillHisto("hBDT_jet_leadingpt2","EVT_dm100", thesample,tree_jet_pt2 ,1.);fillHisto("hBDT_muon_muon_dR","EVT_dm100", thesample,tree_lepton_lepton_dR->at(0),1.);fillHisto("hBDT_muon_muon_dPhi","EVT_dm100", thesample,tree_lepton_lepton_dPhi->at(0),1.);fillHisto("hBDT_muon_muon_dEta","EVT_dm100", thesample,tree_lepton_lepton_dEta->at(0),1.);fillHisto("hBDT_jet_jet_dR","EVT_dm100", thesample,tree_jet_jet_dR->at(0),1.);fillHisto("hBDT_jet_jet_dPhi","EVT_dm100", thesample,tree_jet_jet_dPhi->at(0),1.);fillHisto("hBDT_jet_jet_dEta","EVT_dm100", thesample,tree_jet_jet_dEta->at(0),1.);fillHisto("hBDT_muon_jet_dRmin","EVT_dm100", thesample,tree_muon_jet_dRmin->at(0),1.);fillHisto("hBDT_muon_jet_dRmax","EVT_dm100", thesample,tree_muon_jet_dRmax->at(0),1.);fillHisto("hBDT_HT","EVT_dm100", thesample,tree_HT,1.);fillHisto("hBDT_ST","EVT_dm100", thesample,tree_LT,1.);fillHisto("hBDT_njet","EVT_dm100", thesample,tree_njet,1.);fillHisto("hBDT_muon_nmu","EVT_dm100", thesample,tree_nmu,1.);fillHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm100", thesample,tree_Hemi_LooseBTag_axes->at(0),1.); fillHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm100", thesample,tree_Hemi_MediumBTag_axes->at(0),1.);fillHisto("hBDT_Hemi_TightBTag_axes","EVT_dm100", thesample,tree_Hemi_TightBTag_axes->at(0),1.);}
              if ((Msmu-Mneu) == 150){fillHisto("hSim_ct0","dm150",  thesample,ct0 ,1.);fillHisto("hBDT_MET","EVT_dm150", thesample,tree_PFMet_et,1.);fillHisto("hBDT_TRACK_SIZE","EVT_dm150", thesample,tree_TRACK_SIZE,1.);fillHisto("hBDT_muon_leadingpt","EVT_dm150", thesample,tree_lepton_leadingpt->at(0),1.); fillHisto("hBDT_muon_leadingpt2","EVT_dm150", thesample,tree_lepton_leadingpt2->at(0),1.); fillHisto("hBDT_jet_leadingpt","EVT_dm150", thesample,tree_jet_leadingpt->at(0),1.);fillHisto("hBDT_jet_leadingpt2","EVT_dm150", thesample,tree_jet_pt2 ,1.);fillHisto("hBDT_muon_muon_dR","EVT_dm150", thesample,tree_lepton_lepton_dR->at(0),1.);fillHisto("hBDT_muon_muon_dPhi","EVT_dm150", thesample,tree_lepton_lepton_dPhi->at(0),1.);fillHisto("hBDT_muon_muon_dEta","EVT_dm150", thesample,tree_lepton_lepton_dEta->at(0),1.);fillHisto("hBDT_jet_jet_dR","EVT_dm150", thesample,tree_jet_jet_dR->at(0),1.);fillHisto("hBDT_jet_jet_dPhi","EVT_dm150", thesample,tree_jet_jet_dPhi->at(0),1.);fillHisto("hBDT_jet_jet_dEta","EVT_dm150", thesample,tree_jet_jet_dEta->at(0),1.);fillHisto("hBDT_muon_jet_dRmin","EVT_dm150", thesample,tree_muon_jet_dRmin->at(0),1.);fillHisto("hBDT_muon_jet_dRmax","EVT_dm150", thesample,tree_muon_jet_dRmax->at(0),1.);fillHisto("hBDT_HT","EVT_dm150", thesample,tree_HT,1.);fillHisto("hBDT_ST","EVT_dm150", thesample,tree_LT,1.);fillHisto("hBDT_njet","EVT_dm150", thesample,tree_njet,1.);fillHisto("hBDT_muon_nmu","EVT_dm150", thesample,tree_nmu,1.);fillHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm150", thesample,tree_Hemi_LooseBTag_axes->at(0),1.); fillHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm150", thesample,tree_Hemi_MediumBTag_axes->at(0),1.);fillHisto("hBDT_Hemi_TightBTag_axes","EVT_dm150", thesample,tree_Hemi_TightBTag_axes->at(0),1.);}
              if ((Msmu-Mneu) == 200){fillHisto("hSim_ct0","dm200",  thesample,ct0 ,1.);fillHisto("hBDT_MET","EVT_dm200", thesample,tree_PFMet_et,1.);fillHisto("hBDT_TRACK_SIZE","EVT_dm200", thesample,tree_TRACK_SIZE,1.);fillHisto("hBDT_muon_leadingpt","EVT_dm200", thesample,tree_lepton_leadingpt->at(0),1.); fillHisto("hBDT_muon_leadingpt2","EVT_dm200", thesample,tree_lepton_leadingpt2->at(0),1.); fillHisto("hBDT_jet_leadingpt","EVT_dm200", thesample,tree_jet_leadingpt->at(0),1.);fillHisto("hBDT_jet_leadingpt2","EVT_dm200", thesample,tree_jet_pt2 ,1.);fillHisto("hBDT_muon_muon_dR","EVT_dm200", thesample,tree_lepton_lepton_dR->at(0),1.);fillHisto("hBDT_muon_muon_dPhi","EVT_dm200", thesample,tree_lepton_lepton_dPhi->at(0),1.);fillHisto("hBDT_muon_muon_dEta","EVT_dm200", thesample,tree_lepton_lepton_dEta->at(0),1.);fillHisto("hBDT_jet_jet_dR","EVT_dm200", thesample,tree_jet_jet_dR->at(0),1.);fillHisto("hBDT_jet_jet_dPhi","EVT_dm200", thesample,tree_jet_jet_dPhi->at(0),1.);fillHisto("hBDT_jet_jet_dEta","EVT_dm200", thesample,tree_jet_jet_dEta->at(0),1.);fillHisto("hBDT_muon_jet_dRmin","EVT_dm200", thesample,tree_muon_jet_dRmin->at(0),1.);fillHisto("hBDT_muon_jet_dRmax","EVT_dm200", thesample,tree_muon_jet_dRmax->at(0),1.);fillHisto("hBDT_HT","EVT_dm200", thesample,tree_HT,1.);fillHisto("hBDT_ST","EVT_dm200", thesample,tree_LT,1.);fillHisto("hBDT_njet","EVT_dm200", thesample,tree_njet,1.);fillHisto("hBDT_muon_nmu","EVT_dm200", thesample,tree_nmu,1.);fillHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm200", thesample,tree_Hemi_LooseBTag_axes->at(0),1.); fillHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm200", thesample,tree_Hemi_MediumBTag_axes->at(0),1.);fillHisto("hBDT_Hemi_TightBTag_axes","EVT_dm200", thesample,tree_Hemi_TightBTag_axes->at(0),1.);}
              if ((Msmu-Mneu) == 250){fillHisto("hSim_ct0","dm250",  thesample,ct0 ,1.);fillHisto("hBDT_MET","EVT_dm250", thesample,tree_PFMet_et,1.);fillHisto("hBDT_TRACK_SIZE","EVT_dm250", thesample,tree_TRACK_SIZE,1.);fillHisto("hBDT_muon_leadingpt","EVT_dm250", thesample,tree_lepton_leadingpt->at(0),1.); fillHisto("hBDT_muon_leadingpt2","EVT_dm250", thesample,tree_lepton_leadingpt2->at(0),1.); fillHisto("hBDT_jet_leadingpt","EVT_dm250", thesample,tree_jet_leadingpt->at(0),1.);fillHisto("hBDT_jet_leadingpt2","EVT_dm250", thesample,tree_jet_pt2 ,1.);fillHisto("hBDT_muon_muon_dR","EVT_dm250", thesample,tree_lepton_lepton_dR->at(0),1.);fillHisto("hBDT_muon_muon_dPhi","EVT_dm250", thesample,tree_lepton_lepton_dPhi->at(0),1.);fillHisto("hBDT_muon_muon_dEta","EVT_dm250", thesample,tree_lepton_lepton_dEta->at(0),1.);fillHisto("hBDT_jet_jet_dR","EVT_dm250", thesample,tree_jet_jet_dR->at(0),1.);fillHisto("hBDT_jet_jet_dPhi","EVT_dm250", thesample,tree_jet_jet_dPhi->at(0),1.);fillHisto("hBDT_jet_jet_dEta","EVT_dm250", thesample,tree_jet_jet_dEta->at(0),1.);fillHisto("hBDT_muon_jet_dRmin","EVT_dm250", thesample,tree_muon_jet_dRmin->at(0),1.);fillHisto("hBDT_muon_jet_dRmax","EVT_dm250", thesample,tree_muon_jet_dRmax->at(0),1.);fillHisto("hBDT_HT","EVT_dm250", thesample,tree_HT,1.);fillHisto("hBDT_ST","EVT_dm250", thesample,tree_LT,1.);fillHisto("hBDT_njet","EVT_dm250", thesample,tree_njet,1.);fillHisto("hBDT_muon_nmu","EVT_dm250", thesample,tree_nmu,1.);fillHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm250", thesample,tree_Hemi_LooseBTag_axes->at(0),1.); fillHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm250", thesample,tree_Hemi_MediumBTag_axes->at(0),1.);fillHisto("hBDT_Hemi_TightBTag_axes","EVT_dm250", thesample,tree_Hemi_TightBTag_axes->at(0),1.);}
              if ((Msmu-Mneu) == 300){fillHisto("hSim_ct0","dm300",  thesample,ct0 ,1.);fillHisto("hBDT_MET","EVT_dm300", thesample,tree_PFMet_et,1.);fillHisto("hBDT_TRACK_SIZE","EVT_dm300", thesample,tree_TRACK_SIZE,1.);fillHisto("hBDT_muon_leadingpt","EVT_dm300", thesample,tree_lepton_leadingpt->at(0),1.); fillHisto("hBDT_muon_leadingpt2","EVT_dm300", thesample,tree_lepton_leadingpt2->at(0),1.); fillHisto("hBDT_jet_leadingpt","EVT_dm300", thesample,tree_jet_leadingpt->at(0),1.);fillHisto("hBDT_jet_leadingpt2","EVT_dm300", thesample,tree_jet_pt2 ,1.);fillHisto("hBDT_muon_muon_dR","EVT_dm300", thesample,tree_lepton_lepton_dR->at(0),1.);fillHisto("hBDT_muon_muon_dPhi","EVT_dm300", thesample,tree_lepton_lepton_dPhi->at(0),1.);fillHisto("hBDT_muon_muon_dEta","EVT_dm300", thesample,tree_lepton_lepton_dEta->at(0),1.);fillHisto("hBDT_jet_jet_dR","EVT_dm300", thesample,tree_jet_jet_dR->at(0),1.);fillHisto("hBDT_jet_jet_dPhi","EVT_dm300", thesample,tree_jet_jet_dPhi->at(0),1.);fillHisto("hBDT_jet_jet_dEta","EVT_dm300", thesample,tree_jet_jet_dEta->at(0),1.);fillHisto("hBDT_muon_jet_dRmin","EVT_dm300", thesample,tree_muon_jet_dRmin->at(0),1.);fillHisto("hBDT_muon_jet_dRmax","EVT_dm300", thesample,tree_muon_jet_dRmax->at(0),1.);fillHisto("hBDT_HT","EVT_dm300", thesample,tree_HT,1.);fillHisto("hBDT_ST","EVT_dm300", thesample,tree_LT,1.);fillHisto("hBDT_njet","EVT_dm300", thesample,tree_njet,1.);fillHisto("hBDT_muon_nmu","EVT_dm300", thesample,tree_nmu,1.);fillHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm300", thesample,tree_Hemi_LooseBTag_axes->at(0),1.); fillHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm300", thesample,tree_Hemi_MediumBTag_axes->at(0),1.);fillHisto("hBDT_Hemi_TightBTag_axes","EVT_dm300", thesample,tree_Hemi_TightBTag_axes->at(0),1.);}
              if ((Msmu-Mneu) == 320){fillHisto("hSim_ct0","dm320",  thesample,ct0 ,1.);fillHisto("hBDT_MET","EVT_dm320", thesample,tree_PFMet_et,1.);fillHisto("hBDT_TRACK_SIZE","EVT_dm320", thesample,tree_TRACK_SIZE,1.);fillHisto("hBDT_muon_leadingpt","EVT_dm320", thesample,tree_lepton_leadingpt->at(0),1.); fillHisto("hBDT_muon_leadingpt2","EVT_dm320", thesample,tree_lepton_leadingpt2->at(0),1.); fillHisto("hBDT_jet_leadingpt","EVT_dm320", thesample,tree_jet_leadingpt->at(0),1.);fillHisto("hBDT_jet_leadingpt2","EVT_dm320", thesample,tree_jet_pt2 ,1.);fillHisto("hBDT_muon_muon_dR","EVT_dm320", thesample,tree_lepton_lepton_dR->at(0),1.);fillHisto("hBDT_muon_muon_dPhi","EVT_dm320", thesample,tree_lepton_lepton_dPhi->at(0),1.);fillHisto("hBDT_muon_muon_dEta","EVT_dm320", thesample,tree_lepton_lepton_dEta->at(0),1.);fillHisto("hBDT_jet_jet_dR","EVT_dm320", thesample,tree_jet_jet_dR->at(0),1.);fillHisto("hBDT_jet_jet_dPhi","EVT_dm320", thesample,tree_jet_jet_dPhi->at(0),1.);fillHisto("hBDT_jet_jet_dEta","EVT_dm320", thesample,tree_jet_jet_dEta->at(0),1.);fillHisto("hBDT_muon_jet_dRmin","EVT_dm320", thesample,tree_muon_jet_dRmin->at(0),1.);fillHisto("hBDT_muon_jet_dRmax","EVT_dm320", thesample,tree_muon_jet_dRmax->at(0),1.);fillHisto("hBDT_HT","EVT_dm320", thesample,tree_HT,1.);fillHisto("hBDT_ST","EVT_dm320", thesample,tree_LT,1.);fillHisto("hBDT_njet","EVT_dm320", thesample,tree_njet,1.);fillHisto("hBDT_muon_nmu","EVT_dm320", thesample,tree_nmu,1.);fillHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm320", thesample,tree_Hemi_LooseBTag_axes->at(0),1.); fillHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm320", thesample,tree_Hemi_MediumBTag_axes->at(0),1.);fillHisto("hBDT_Hemi_TightBTag_axes","EVT_dm320", thesample,tree_Hemi_TightBTag_axes->at(0),1.);}
            
            } 

        }    // End Loop on GenParticle
        
      fillHisto("hGen_Msmu","noSel",  thesample, tree_smu_mass ,1.);
      fillHisto("hGen_Mneu","noSel",  thesample, tree_neu_mass ,1.);

    }
       
  // BDT EVTS SELCTION VARIABLES

  float tree_jet_pt2 = 0;
  if (tree_jet_leadingpt2->size()>0){tree_jet_pt2 = tree_jet_leadingpt2->at(0);}
    fillHisto("hBDT_MET","EVT", thesample,tree_PFMet_et,1.);
    fillHisto("hBDT_TRACK_SIZE","EVT", thesample,tree_TRACK_SIZE,1.);
    fillHisto("hBDT_muon_leadingpt","EVT", thesample,tree_lepton_leadingpt->at(0),1.);
    fillHisto("hBDT_muon_leadingpt2","EVT", thesample,tree_lepton_leadingpt2->at(0),1.);
    fillHisto("hBDT_jet_leadingpt","EVT", thesample,tree_jet_leadingpt->at(0),1.);
    fillHisto("hBDT_jet_leadingpt2","EVT", thesample,tree_jet_pt2 ,1.);
    fillHisto("hBDT_muon_muon_dR","EVT", thesample,tree_lepton_lepton_dR->at(0),1.);
    fillHisto("hBDT_muon_muon_dPhi","EVT", thesample,tree_lepton_lepton_dPhi->at(0),1.);
    fillHisto("hBDT_muon_muon_dEta","EVT", thesample,tree_lepton_lepton_dEta->at(0),1.);
    fillHisto("hBDT_jet_jet_dR","EVT", thesample,tree_jet_jet_dR->at(0),1.);
    fillHisto("hBDT_jet_jet_dPhi","EVT", thesample,tree_jet_jet_dPhi->at(0),1.);
    fillHisto("hBDT_jet_jet_dEta","EVT", thesample,tree_jet_jet_dEta->at(0),1.);
    fillHisto("hBDT_muon_jet_dRmin","EVT", thesample,tree_muon_jet_dRmin->at(0),1.);
    fillHisto("hBDT_muon_jet_dRmax","EVT", thesample,tree_muon_jet_dRmax->at(0),1.);
    fillHisto("hBDT_HT","EVT", thesample,tree_HT,1.);
    fillHisto("hBDT_ST","EVT", thesample,tree_LT,1.);
    fillHisto("hBDT_njet","EVT", thesample,tree_njet,1.);
    fillHisto("hBDT_muon_nmu","EVT", thesample,tree_muon_pt->size(),1.);
    fillHisto("hBDT_Hemi_LooseBTag_axes","EVT", thesample,tree_Hemi_LooseBTag_axes->at(0),1.);
    fillHisto("hBDT_Hemi_MediumBTag_axes","EVT", thesample,tree_Hemi_MediumBTag_axes->at(0),1.);
    fillHisto("hBDT_Hemi_TightBTag_axes","EVT", thesample,tree_Hemi_TightBTag_axes->at(0),1.);
    


      float HighMass = tree_HemiMu_mass->at(1);
      if(tree_HemiMu_mass->at(0)>tree_HemiMu_mass->at(1))
        {
          HighMass = tree_HemiMu_mass->at(0);
        }
      fillHisto("hData_HighSmuonMass","NoSel",thesample,HighMass,1);
      if (tree_Evts_MVAval <  EVTSWP){continue;}
      // if (HighMass < 200) continue;
      fillHisto("DiMuon_Mass","EVTSel", thesample,  tree_Mmumu, XS/nentries); //!!!!
      nEvts++;
 
      //*******************************
      //loop on Muons
      //*******************************

      for(unsigned int iMuon = 0; iMuon <tree_muon_pt->size(); iMuon ++)//test
         {
            fillHisto("RecoMuo_pT", "noSel",  thesample,  tree_muon_pt->at(iMuon) , 1.);
            
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
  //  std::cout<<"here"<<std::endl;
      //*******************************
      //loop on Sec. Interactions
      //*******************************

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

                  //to get a nice view of the inner tracker x vs y
                  if(tree_SecInt_selec->at(iSecInt)&&tree_SecInt_layer->at(iSecInt)!=0 && abs(tree_SecInt_x->at(iSecInt))<25 && abs(tree_SecInt_y->at(iSecInt))<25 && abs(tree_SecInt_eta->at(iSecInt))<1.4 && abs(tree_SecInt_z->at(iSecInt))<27)
                    {
                      fillHisto2D("hData_reco_SecInt_xy","TrackerMatched", thesample, tree_SecInt_x->at(iSecInt),tree_SecInt_y->at(iSecInt),1.);
                    }
                  if (tree_SecInt_selec->at(iSecInt) && abs(tree_SecInt_x->at(iSecInt))<25 && abs(tree_SecInt_y->at(iSecInt))<25 && abs(tree_SecInt_eta->at(iSecInt))<1.4&&abs(tree_SecInt_z->at(iSecInt))<27)
                    {
                      fillHisto2D("hData_reco_SecInt_xy","Selec", thesample, tree_SecInt_x->at(iSecInt),tree_SecInt_y->at(iSecInt),1.);
                    }

                  // to get a nice view of the tracker in the r vs z plane
                  if (tree_SecInt_selec->at(iSecInt) && tree_SecInt_layer->at(iSecInt)!=0 && abs(tree_SecInt_r->at(iSecInt))<70 && abs(tree_SecInt_z->at(iSecInt))<120)
                    {
                      fillHisto2D("hData_reco_SecInt_rz","TrackerMatched", thesample, abs(tree_SecInt_z->at(iSecInt)),tree_SecInt_r->at(iSecInt),1.);
                    }
                  if (tree_SecInt_selec->at(iSecInt) && tree_SecInt_layer->at(iSecInt)!=0 && abs(tree_SecInt_r->at(iSecInt))<70 && abs(tree_SecInt_z->at(iSecInt))<120)
                    {
                      fillHisto2D("hData_reco_SecInt_rz","Selec", thesample,abs(tree_SecInt_z->at(iSecInt)), tree_SecInt_r->at(iSecInt),1.);
                    }
               }
         }

      //*******************************
      //loop on jets
      //*******************************

    for (unsigned int k = 0 ; k < tree_jet_pt->size(); k++)
      {
        fillHisto("hData_jet_pt","",            thesample, tree_jet_pt->at(k),1.);
        fillHisto("hData_jet_eta","",           thesample, tree_jet_eta->at(k),1.);
        fillHisto("hData_jet_btag_Deepjet","",  thesample, tree_jet_btag_DeepJet->at(k),1.);
        fillHisto("hData_jet_HadronFlavour","", thesample, tree_jet_HadronFlavour->at(k),1.);
        
        if (abs(tree_jet_HadronFlavour->at(k))==5 && tree_jet_pileupID->at(k)>=0)//Probably an issue here since it is said that pileup ID has to be applied before any corrections 
                                                                                 // but corrections are already applied before running the code..... 
         {
            fillHisto2D("hData_BtagEff_Denom","",   thesample, tree_jet_pt->at(k),tree_jet_eta->at(k),1 );
            if (tree_jet_btag_DeepJet->at(k) >= 0.7100)//tight WP of Btag
               {
                  fillHisto2D("hData_BtagEff_Num","",     thesample, tree_jet_pt->at(k),tree_jet_eta->at(k),1 );
               }
         }
        
        
      }
    // std::cout<<"here2"<<std::endl;
      //*******************************
      //loop on tracks
      //*******************************

      n_TotalTrks += tree_track_MVAval->size();

    if ( Signal)
    {        
      
      int Msmu = tree_smu_mass;
      int Mneu = tree_neu_mass;
      if ( Msmu > 180 && Msmu < 220 )      Msmu = 200;
      else if ( Msmu > 230 && Msmu < 270 ) Msmu = 250;
      else if ( Msmu > 280 && Msmu < 320 ) Msmu = 300;
      else if ( Msmu > 380 && Msmu < 420 ) Msmu = 400;
      else if ( Msmu > 480 && Msmu < 520 ) Msmu = 500;
      // else cout << " !!! smu mass out of range !!! " << Msmu;
      if ( Mneu > 170 && Mneu < 190 )      Mneu = 180;
      else if ( Mneu > 190 && Mneu < 210 ) Mneu = 200;
      else if ( Mneu > 240 && Mneu < 260 ) Mneu = 250;
      else if ( Mneu > 270 && Mneu < 290 ) Mneu = 280;
      else if ( Mneu > 290 && Mneu < 310 ) Mneu = 300;
      else if ( Mneu > 340 && Mneu < 360 ) Mneu = 350;
      else if ( Mneu > 370 && Mneu < 390 ) Mneu = 380;
      else if ( Mneu > 390 && Mneu < 410 ) Mneu = 400;
      else if ( Mneu > 440 && Mneu < 460 ) Mneu = 450;
      else if ( Mneu > 470 && Mneu < 490 ) Mneu = 480;
      // else cout << " !!! neu mass out of range !!! " << Mneu;
    //  std::cout<<"Msmu and Mneu : "<<Msmu<<" and "<<Mneu<<std::endl;
      if ((Msmu-Mneu) == 20){ for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hSim_TRKBDT","dm20",  thesample,tree_track_MVAval->at(iTrk) ,1.);}}
      if ((Msmu-Mneu) == 50){ for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hSim_TRKBDT","dm50",  thesample,tree_track_MVAval->at(iTrk) ,1.);}}
      if ((Msmu-Mneu) == 100){ for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hSim_TRKBDT","dm100",  thesample,tree_track_MVAval->at(iTrk) ,1.);}}
      if ((Msmu-Mneu) == 150){ for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hSim_TRKBDT","dm150",  thesample,tree_track_MVAval->at(iTrk) ,1.);}}
      if ((Msmu-Mneu) == 200){ for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hSim_TRKBDT","dm200",  thesample,tree_track_MVAval->at(iTrk) ,1.);}}
      if ((Msmu-Mneu) == 250){ for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hSim_TRKBDT","dm250",  thesample,tree_track_MVAval->at(iTrk) ,1.);}}
      if ((Msmu-Mneu) == 300){ for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hSim_TRKBDT","dm300",  thesample,tree_track_MVAval->at(iTrk) ,1.);}}
      if ((Msmu-Mneu) == 320){ for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hSim_TRKBDT","dm320",  thesample,tree_track_MVAval->at(iTrk) ,1.);}}

      int ngenpart =  tree_genParticle_pt->size();
      for (int i=0; i<ngenpart; i++)    // Loop on GenParticle
        {
          float pdgId = tree_genParticle_pdgId->at(i); 
          float mother_pdgId = tree_genParticle_mother_pdgId->at(i); 
          float ct0 = tree_genParticle_ct0->at(i);

          // top quark from neutralino
          if ( abs(pdgId) == 6 && abs(mother_pdgId) == 1000023 ) 
            {            
              if ((Msmu-Mneu) == 20){fillHisto("hSim_ct0","dm20_TRKBDT",  thesample,ct0 ,1.);for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hBDT_lost","TRK_dm20", thesample,tree_track_lost->at(iTrk),1.);fillHisto("hBDT_dz","TRK_dm20", thesample,tree_track_dz->at(iTrk),1.);fillHisto("hBDT_dxy","TRK_dm20", thesample,tree_track_dxy->at(iTrk),1.);fillHisto("hBDT_pt","TRK_dm20", thesample,tree_track_pt->at(iTrk),1.);fillHisto("hBDT_eta","TRK_dm20", thesample,tree_track_eta->at(iTrk),1.);fillHisto("hBDT_NChi2","TRK_dm20", thesample,tree_track_NChi2->at(iTrk),1.);fillHisto("hBDT_nhits","TRK_dm20", thesample,tree_track_nHit->at(iTrk),1.);fillHisto("hBDT_iJet","TRK_dm20", thesample,tree_track_iJet->at(iTrk),1.);fillHisto("hBDT_drSig","TRK_dm20", thesample,tree_track_drSig->at(iTrk),1.);fillHisto("hBDT_dzSig","TRK_dm20", thesample,tree_track_dzSig->at(iTrk),1.);fillHisto("hBDT_ntrk10","TRK_dm20", thesample,tree_track_ntrk10->at(iTrk),1.);fillHisto("hBDT_ntrk20","TRK_dm20", thesample,tree_track_ntrk20->at(iTrk),1.);fillHisto("hBDT_ntrk30","TRK_dm20", thesample,tree_track_ntrk30->at(iTrk),1.);fillHisto("hBDT_ntrk40","TRK_dm20", thesample,tree_track_ntrk40->at(iTrk),1.);fillHisto("hBDT_Hemi_dR","TRK_dm20", thesample,tree_track_Hemi_dR->at(iTrk),1.);fillHisto("hBDT_Hemi_dRmax","TRK_dm20", thesample,tree_track_Hemi_dRmax->at(iTrk),1.);}}
              if ((Msmu-Mneu) == 50){fillHisto("hSim_ct0","dm50_TRKBDT",  thesample,ct0 ,1.);for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hBDT_lost","TRK_dm50", thesample,tree_track_lost->at(iTrk),1.);fillHisto("hBDT_dz","TRK_dm50", thesample,tree_track_dz->at(iTrk),1.);fillHisto("hBDT_dxy","TRK_dm50", thesample,tree_track_dxy->at(iTrk),1.);fillHisto("hBDT_pt","TRK_dm50", thesample,tree_track_pt->at(iTrk),1.);fillHisto("hBDT_eta","TRK_dm50", thesample,tree_track_eta->at(iTrk),1.);fillHisto("hBDT_NChi2","TRK_dm50", thesample,tree_track_NChi2->at(iTrk),1.);fillHisto("hBDT_nhits","TRK_dm50", thesample,tree_track_nHit->at(iTrk),1.);fillHisto("hBDT_iJet","TRK_dm50", thesample,tree_track_iJet->at(iTrk),1.);fillHisto("hBDT_drSig","TRK_dm50", thesample,tree_track_drSig->at(iTrk),1.);fillHisto("hBDT_dzSig","TRK_dm50", thesample,tree_track_dzSig->at(iTrk),1.);fillHisto("hBDT_ntrk10","TRK_dm50", thesample,tree_track_ntrk10->at(iTrk),1.);fillHisto("hBDT_ntrk20","TRK_dm50", thesample,tree_track_ntrk20->at(iTrk),1.);fillHisto("hBDT_ntrk30","TRK_dm50", thesample,tree_track_ntrk30->at(iTrk),1.);fillHisto("hBDT_ntrk40","TRK_dm50", thesample,tree_track_ntrk40->at(iTrk),1.);fillHisto("hBDT_Hemi_dR","TRK_dm50", thesample,tree_track_Hemi_dR->at(iTrk),1.);fillHisto("hBDT_Hemi_dRmax","TRK_dm50", thesample,tree_track_Hemi_dRmax->at(iTrk),1.);}}
              if ((Msmu-Mneu) == 100){fillHisto("hSim_ct0","dm100_TRKBDT",  thesample,ct0 ,1.);for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hBDT_lost","TRK_dm100", thesample,tree_track_lost->at(iTrk),1.);fillHisto("hBDT_dz","TRK_dm100", thesample,tree_track_dz->at(iTrk),1.);fillHisto("hBDT_dxy","TRK_dm100", thesample,tree_track_dxy->at(iTrk),1.);fillHisto("hBDT_pt","TRK_dm100", thesample,tree_track_pt->at(iTrk),1.);fillHisto("hBDT_eta","TRK_dm100", thesample,tree_track_eta->at(iTrk),1.);fillHisto("hBDT_NChi2","TRK_dm100", thesample,tree_track_NChi2->at(iTrk),1.);fillHisto("hBDT_nhits","TRK_dm100", thesample,tree_track_nHit->at(iTrk),1.);fillHisto("hBDT_iJet","TRK_dm100", thesample,tree_track_iJet->at(iTrk),1.);fillHisto("hBDT_drSig","TRK_dm100", thesample,tree_track_drSig->at(iTrk),1.);fillHisto("hBDT_dzSig","TRK_dm100", thesample,tree_track_dzSig->at(iTrk),1.);fillHisto("hBDT_ntrk10","TRK_dm100", thesample,tree_track_ntrk10->at(iTrk),1.);fillHisto("hBDT_ntrk20","TRK_dm100", thesample,tree_track_ntrk20->at(iTrk),1.);fillHisto("hBDT_ntrk30","TRK_dm100", thesample,tree_track_ntrk30->at(iTrk),1.);fillHisto("hBDT_ntrk40","TRK_dm100", thesample,tree_track_ntrk40->at(iTrk),1.);fillHisto("hBDT_Hemi_dR","TRK_dm100", thesample,tree_track_Hemi_dR->at(iTrk),1.);fillHisto("hBDT_Hemi_dRmax","TRK_dm100", thesample,tree_track_Hemi_dRmax->at(iTrk),1.);}}
              if ((Msmu-Mneu) == 150){fillHisto("hSim_ct0","dm150_TRKBDT",  thesample,ct0 ,1.);for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hBDT_lost","TRK_dm150", thesample,tree_track_lost->at(iTrk),1.);fillHisto("hBDT_dz","TRK_dm150", thesample,tree_track_dz->at(iTrk),1.);fillHisto("hBDT_dxy","TRK_dm150", thesample,tree_track_dxy->at(iTrk),1.);fillHisto("hBDT_pt","TRK_dm150", thesample,tree_track_pt->at(iTrk),1.);fillHisto("hBDT_eta","TRK_dm150", thesample,tree_track_eta->at(iTrk),1.);fillHisto("hBDT_NChi2","TRK_dm150", thesample,tree_track_NChi2->at(iTrk),1.);fillHisto("hBDT_nhits","TRK_dm150", thesample,tree_track_nHit->at(iTrk),1.);fillHisto("hBDT_iJet","TRK_dm150", thesample,tree_track_iJet->at(iTrk),1.);fillHisto("hBDT_drSig","TRK_dm150", thesample,tree_track_drSig->at(iTrk),1.);fillHisto("hBDT_dzSig","TRK_dm150", thesample,tree_track_dzSig->at(iTrk),1.);fillHisto("hBDT_ntrk10","TRK_dm150", thesample,tree_track_ntrk10->at(iTrk),1.);fillHisto("hBDT_ntrk20","TRK_dm150", thesample,tree_track_ntrk20->at(iTrk),1.);fillHisto("hBDT_ntrk30","TRK_dm150", thesample,tree_track_ntrk30->at(iTrk),1.);fillHisto("hBDT_ntrk40","TRK_dm150", thesample,tree_track_ntrk40->at(iTrk),1.);fillHisto("hBDT_Hemi_dR","TRK_dm150", thesample,tree_track_Hemi_dR->at(iTrk),1.);fillHisto("hBDT_Hemi_dRmax","TRK_dm150", thesample,tree_track_Hemi_dRmax->at(iTrk),1.);}}
              if ((Msmu-Mneu) == 200){fillHisto("hSim_ct0","dm200_TRKBDT",  thesample,ct0 ,1.);for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hBDT_lost","TRK_dm200", thesample,tree_track_lost->at(iTrk),1.);fillHisto("hBDT_dz","TRK_dm200", thesample,tree_track_dz->at(iTrk),1.);fillHisto("hBDT_dxy","TRK_dm200", thesample,tree_track_dxy->at(iTrk),1.);fillHisto("hBDT_pt","TRK_dm200", thesample,tree_track_pt->at(iTrk),1.);fillHisto("hBDT_eta","TRK_dm200", thesample,tree_track_eta->at(iTrk),1.);fillHisto("hBDT_NChi2","TRK_dm200", thesample,tree_track_NChi2->at(iTrk),1.);fillHisto("hBDT_nhits","TRK_dm200", thesample,tree_track_nHit->at(iTrk),1.);fillHisto("hBDT_iJet","TRK_dm200", thesample,tree_track_iJet->at(iTrk),1.);fillHisto("hBDT_drSig","TRK_dm200", thesample,tree_track_drSig->at(iTrk),1.);fillHisto("hBDT_dzSig","TRK_dm200", thesample,tree_track_dzSig->at(iTrk),1.);fillHisto("hBDT_ntrk10","TRK_dm200", thesample,tree_track_ntrk10->at(iTrk),1.);fillHisto("hBDT_ntrk20","TRK_dm200", thesample,tree_track_ntrk20->at(iTrk),1.);fillHisto("hBDT_ntrk30","TRK_dm200", thesample,tree_track_ntrk30->at(iTrk),1.);fillHisto("hBDT_ntrk40","TRK_dm200", thesample,tree_track_ntrk40->at(iTrk),1.);fillHisto("hBDT_Hemi_dR","TRK_dm200", thesample,tree_track_Hemi_dR->at(iTrk),1.);fillHisto("hBDT_Hemi_dRmax","TRK_dm200", thesample,tree_track_Hemi_dRmax->at(iTrk),1.);}}
              if ((Msmu-Mneu) == 250){fillHisto("hSim_ct0","dm250_TRKBDT",  thesample,ct0 ,1.);for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hBDT_lost","TRK_dm250", thesample,tree_track_lost->at(iTrk),1.);fillHisto("hBDT_dz","TRK_dm250", thesample,tree_track_dz->at(iTrk),1.);fillHisto("hBDT_dxy","TRK_dm250", thesample,tree_track_dxy->at(iTrk),1.);fillHisto("hBDT_pt","TRK_dm250", thesample,tree_track_pt->at(iTrk),1.);fillHisto("hBDT_eta","TRK_dm250", thesample,tree_track_eta->at(iTrk),1.);fillHisto("hBDT_NChi2","TRK_dm250", thesample,tree_track_NChi2->at(iTrk),1.);fillHisto("hBDT_nhits","TRK_dm250", thesample,tree_track_nHit->at(iTrk),1.);fillHisto("hBDT_iJet","TRK_dm250", thesample,tree_track_iJet->at(iTrk),1.);fillHisto("hBDT_drSig","TRK_dm250", thesample,tree_track_drSig->at(iTrk),1.);fillHisto("hBDT_dzSig","TRK_dm250", thesample,tree_track_dzSig->at(iTrk),1.);fillHisto("hBDT_ntrk10","TRK_dm250", thesample,tree_track_ntrk10->at(iTrk),1.);fillHisto("hBDT_ntrk20","TRK_dm250", thesample,tree_track_ntrk20->at(iTrk),1.);fillHisto("hBDT_ntrk30","TRK_dm250", thesample,tree_track_ntrk30->at(iTrk),1.);fillHisto("hBDT_ntrk40","TRK_dm250", thesample,tree_track_ntrk40->at(iTrk),1.);fillHisto("hBDT_Hemi_dR","TRK_dm250", thesample,tree_track_Hemi_dR->at(iTrk),1.);fillHisto("hBDT_Hemi_dRmax","TRK_dm250", thesample,tree_track_Hemi_dRmax->at(iTrk),1.);}}
              if ((Msmu-Mneu) == 300){fillHisto("hSim_ct0","dm300_TRKBDT",  thesample,ct0 ,1.);for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hBDT_lost","TRK_dm300", thesample,tree_track_lost->at(iTrk),1.);fillHisto("hBDT_dz","TRK_dm300", thesample,tree_track_dz->at(iTrk),1.);fillHisto("hBDT_dxy","TRK_dm300", thesample,tree_track_dxy->at(iTrk),1.);fillHisto("hBDT_pt","TRK_dm300", thesample,tree_track_pt->at(iTrk),1.);fillHisto("hBDT_eta","TRK_dm300", thesample,tree_track_eta->at(iTrk),1.);fillHisto("hBDT_NChi2","TRK_dm300", thesample,tree_track_NChi2->at(iTrk),1.);fillHisto("hBDT_nhits","TRK_dm300", thesample,tree_track_nHit->at(iTrk),1.);fillHisto("hBDT_iJet","TRK_dm300", thesample,tree_track_iJet->at(iTrk),1.);fillHisto("hBDT_drSig","TRK_dm300", thesample,tree_track_drSig->at(iTrk),1.);fillHisto("hBDT_dzSig","TRK_dm300", thesample,tree_track_dzSig->at(iTrk),1.);fillHisto("hBDT_ntrk10","TRK_dm300", thesample,tree_track_ntrk10->at(iTrk),1.);fillHisto("hBDT_ntrk20","TRK_dm300", thesample,tree_track_ntrk20->at(iTrk),1.);fillHisto("hBDT_ntrk30","TRK_dm300", thesample,tree_track_ntrk30->at(iTrk),1.);fillHisto("hBDT_ntrk40","TRK_dm300", thesample,tree_track_ntrk40->at(iTrk),1.);fillHisto("hBDT_Hemi_dR","TRK_dm300", thesample,tree_track_Hemi_dR->at(iTrk),1.);fillHisto("hBDT_Hemi_dRmax","TRK_dm300", thesample,tree_track_Hemi_dRmax->at(iTrk),1.);}}
              if ((Msmu-Mneu) == 320){fillHisto("hSim_ct0","dm320_TRKBDT",  thesample,ct0 ,1.);for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++){fillHisto("hBDT_lost","TRK_dm320", thesample,tree_track_lost->at(iTrk),1.);fillHisto("hBDT_dz","TRK_dm320", thesample,tree_track_dz->at(iTrk),1.);fillHisto("hBDT_dxy","TRK_dm320", thesample,tree_track_dxy->at(iTrk),1.);fillHisto("hBDT_pt","TRK_dm320", thesample,tree_track_pt->at(iTrk),1.);fillHisto("hBDT_eta","TRK_dm320", thesample,tree_track_eta->at(iTrk),1.);fillHisto("hBDT_NChi2","TRK_dm320", thesample,tree_track_NChi2->at(iTrk),1.);fillHisto("hBDT_nhits","TRK_dm320", thesample,tree_track_nHit->at(iTrk),1.);fillHisto("hBDT_iJet","TRK_dm320", thesample,tree_track_iJet->at(iTrk),1.);fillHisto("hBDT_drSig","TRK_dm320", thesample,tree_track_drSig->at(iTrk),1.);fillHisto("hBDT_dzSig","TRK_dm320", thesample,tree_track_dzSig->at(iTrk),1.);fillHisto("hBDT_ntrk10","TRK_dm320", thesample,tree_track_ntrk10->at(iTrk),1.);fillHisto("hBDT_ntrk20","TRK_dm320", thesample,tree_track_ntrk20->at(iTrk),1.);fillHisto("hBDT_ntrk30","TRK_dm320", thesample,tree_track_ntrk30->at(iTrk),1.);fillHisto("hBDT_ntrk40","TRK_dm320", thesample,tree_track_ntrk40->at(iTrk),1.);fillHisto("hBDT_Hemi_dR","TRK_dm320", thesample,tree_track_Hemi_dR->at(iTrk),1.);fillHisto("hBDT_Hemi_dRmax","TRK_dm320", thesample,tree_track_Hemi_dRmax->at(iTrk),1.);}}

            } 

        }    // End Loop on GenParticle  
    }

      for(unsigned int iTrk = 0; iTrk <tree_track_MVAval->size(); iTrk ++)
         {
          //--- Track pre-selected ---------//
          //get total number of tracks bfore any selection with tree_TRACK_SIZE
               
                fillHisto("hBDT_lost","TRK", thesample,tree_track_lost->at(iTrk),1.);
                fillHisto("hBDT_dz","TRK", thesample,tree_track_dz->at(iTrk),1.);
                fillHisto("hBDT_dxy","TRK", thesample,tree_track_dxy->at(iTrk),1.);
                fillHisto("hBDT_pt","TRK", thesample,tree_track_pt->at(iTrk),1.);
                fillHisto("hBDT_eta","TRK", thesample,tree_track_eta->at(iTrk),1.);
                fillHisto("hBDT_NChi2","TRK", thesample,tree_track_NChi2->at(iTrk),1.);
                fillHisto("hBDT_nhits","TRK", thesample,tree_track_nHit->at(iTrk),1.);
                fillHisto("hBDT_iJet","TRK", thesample,tree_track_iJet->at(iTrk),1.);
                fillHisto("hBDT_drSig","TRK", thesample,tree_track_drSig->at(iTrk),1.);
                fillHisto("hBDT_dzSig","TRK", thesample,tree_track_dzSig->at(iTrk),1.);
                fillHisto("hBDT_ntrk10","TRK", thesample,tree_track_ntrk10->at(iTrk),1.);
                fillHisto("hBDT_ntrk20","TRK", thesample,tree_track_ntrk20->at(iTrk),1.);
                fillHisto("hBDT_ntrk30","TRK", thesample,tree_track_ntrk30->at(iTrk),1.);
                fillHisto("hBDT_ntrk40","TRK", thesample,tree_track_ntrk40->at(iTrk),1.);
                fillHisto("hBDT_Hemi_dR","TRK", thesample,tree_track_Hemi_dR->at(iTrk),1.);
                fillHisto("hBDT_Hemi_dRmax","TRK", thesample,tree_track_Hemi_dRmax->at(iTrk),1.);

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
      if (Signal)
        {
          for(unsigned int iAxes = 0; iAxes < tree_Hemi_dR12->size(); iAxes ++)
            {
              fillHisto("hData_dR_RecoReco","noSel",  thesample, tree_Hemi_dR12->at(iAxes),1.);
              //Generated Info between the two axes (and also witht he reco axes)
              fillHisto("hSim_dR_GenGen","noSel", thesample, tree_Hemi_LLP_dR12->at(iAxes),1);
              fillHisto("hSim_dPhi_GenGen","noSel", thesample, tree_genAxis_dPhineuneu->at(iAxes),1);
              fillHisto("hSim_dEta_GenGen","noSel", thesample, tree_genAxis_dEtaneuneu->at(iAxes),1);
              fillHisto("hSim_dR_GenReco","noSel", thesample, tree_Hemi_LLP_dR->at(iAxes),1);
            }
        }


      //*******************************
      //loop on Reco Vertices
      //*******************************

      int nTightVertex = 0;
      float nSelecBDTVtxTight_evt = 0;
      float nSelecBDTVtxTight_step12_evt = 0;
      int countReco = 0;
      int countRecoGood = 0;
      int countRecoTight12 = 0;
      int countRecoTight12_BDT3 = 0;
      int countRecoTight12_BDT3bis = 0;
      for(unsigned int iVtx = 0; iVtx <tree_Hemi->size(); iVtx ++)
         {
            fillHisto("hData_Hemi_Vtx_NChi2","noSel",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
            fillHisto("hData_Hemi_Vtx_nTrks","noSel",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
            if(Signal)
              {
                float r = sqrt((tree_Hemi_LLP_x->at(iVtx)*tree_Hemi_LLP_x->at(iVtx))+(tree_Hemi_LLP_y->at(iVtx)*tree_Hemi_LLP_y->at(iVtx)));              
                fillHisto("hSim_Hemi_Vtx_r","noSel",      thesample,r,1); 
                fillHisto("hSim_Hemi_Vtx_eta","noSel",    thesample,tree_Hemi_LLP_eta->at(iVtx),1); 
                fillHisto("hSim_Hemi_Vtx_dist","noSel",   thesample,tree_Hemi_LLP_dist->at(iVtx),1);
              }
            TotalnVertex++;

            if (tree_Hemi_Vtx_NChi2->at(iVtx) != -10 )
               {
                  fillHisto("hData_Hemi_Vtx_NChi2","RecoVtx",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_nTrks","RecoVtx",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_dist","RecoVtx",    thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_step","RecoVtx",    thesample,tree_Hemi_Vtx_step->at(iVtx),1.);
                  // fillHisto("hData_Hemi_Vtx_eta","RecoVtx",     thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_r","RecoVtx",       thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_z","RecoVtx",       thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                  fillHisto("hData_Hemi_Vtx_MVAval","RecoVtx",  thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);
                  

                  if ( tree_Hemi_Vtx_NChi2->at(iVtx)>0 && tree_Hemi_Vtx_NChi2->at(iVtx)<10)//Reco Vtx criteria
                     {
                        fillHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_dist","GoodRecoVtx",    thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_step","GoodRecoVtx",    thesample,tree_Hemi_Vtx_step->at(iVtx),1.);
                        // fillHisto("hData_Hemi_Vtx_eta","GoodRecoVtx",     thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_r","GoodRecoVtx",       thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_z","GoodRecoVtx",       thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx",  thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);
                            fillHisto("hBDT_nTrks",     "VTX", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                            fillHisto("hBDT_NChi2",     "VTX", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                            fillHisto("hBDT_step",      "VTX", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);
                            fillHisto("hBDT_r",         "VTX", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                            fillHisto("hBDT_z",         "VTX", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                            fillHisto("hBDT_MWT",       "VTX", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);
                            fillHisto2D("MWTVsnTrks",   "VTX", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                            
                            fillHisto("hBDT_Vtx_Mass",  "VTX", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);
                            fillHisto("hBDT_Hemi_Mass", "VTX", thesample,tree_Hemi_mass->at(iVtx),1.);
                            fillHisto("hBDT_dist",      "VTX", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                            fillHisto("hBDT_ntrk10",    "VTX", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);
                            fillHisto("hBDT_ntrk20",    "VTX", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);
                            fillHisto("hBDT_MeanDCA",   "VTX", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);
                            countReco++;
                            
                        // if (tree_Hemi_Vtx_MVAval->size()!=tree_Hemi_Vtx_step->size())
                        //   {std::cout<<"trere is an issue"<<std::endl;}
                        if (tree_Hemi_Vtx_MVAval->at(iVtx)>VTXWP)
                          {
                            nSelecBDTVtx++;
                            countRecoGood++;
                            
                          }
                        if (Signal)
                          {
                            if (tree_Hemi_LLP_ping->at(iVtx))
                              {
                                fillHisto("hData_Hemi_Vtx_r","Ping",      thesample,tree_Hemi_Vtx_r->at(iVtx),1.); 
                                // fillHisto("hData_Hemi_Vtx_eta","Ping",    thesample,tree_Hemi_Vtx_eta->at(iVtx),1.); 
                                fillHisto("hData_Hemi_Vtx_dist","Ping",   thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);                   
                                float r = sqrt((tree_Hemi_LLP_x->at(iVtx)*tree_Hemi_LLP_x->at(iVtx))+(tree_Hemi_LLP_y->at(iVtx)*tree_Hemi_LLP_y->at(iVtx)));              
                                fillHisto("hSim_Hemi_Vtx_r","Ping",      thesample,r,1); 
                                fillHisto("hSim_Hemi_Vtx_eta","Ping",    thesample,tree_Hemi_LLP_eta->at(iVtx),1); 
                                fillHisto("hSim_Hemi_Vtx_dist","Ping",   thesample,tree_Hemi_LLP_dist->at(iVtx),1);                         
                              }
                          }
                        if (tree_Hemi_Vtx_nTrks->at(iVtx)==2)
                          {
                            nRecoVertex2Trks++;
                          }
                        
                        nRecoVertex++;
                        MeanDistance +=  tree_Hemi_Vtx_dist->at(iVtx);              
                     }
                  

                  // -------------- TIGHT WP ------------------------//
                  if (tree_Hemi_Vtx_step->at(iVtx)==1 || tree_Hemi_Vtx_step->at(iVtx)==2)
                     {
                        fillHisto("hData_Hemi_Vtx_nTrks","RecoVtx_TightWP",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_NChi2","RecoVtx_TightWP",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_dist","RecoVtx_TightWP",    thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                        // fillHisto("hData_Hemi_Vtx_eta","RecoVtx_TightWP",     thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_r","RecoVtx_TightWP",       thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_z","RecoVtx_TightWP",       thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_MVAval","RecoVtx_TightWP",  thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);
                       
                        
                        if ( tree_Hemi_Vtx_NChi2->at(iVtx)>0 && tree_Hemi_Vtx_NChi2->at(iVtx)<10)//Reco Vtx criteria
                           {
                              fillHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx_TightWP",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx_TightWP",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_dist","GoodRecoVtx_TightWP",    thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                              // fillHisto("hData_Hemi_Vtx_eta","GoodRecoVtx_TightWP",     thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_r","GoodRecoVtx_TightWP",       thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_z","GoodRecoVtx_TightWP",       thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx_TightWP",  thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);
                              nTightVertex++;
                              countRecoTight12++;


                              
                              if (tree_Hemi_Vtx_MVAval->at(iVtx)>VTXWP)
                                {
                                  nSelecBDTVtxTight++;
                                  nSelecBDTVtxTight_evt++;
                                  countRecoTight12_BDT3++;
                                  
                                }
                              if (tree_Hemi_Vtx_MVAval_Step1->at(iVtx)>VTXWP)
                                {
                                  nSelecBDTVtxTight_step12++;
                                  nSelecBDTVtxTight_step12_evt++;
                                  countRecoTight12_BDT3bis++;
                                  
                                }

                              if (Signal)
                                {
                                  if (tree_Hemi_LLP_ping->at(iVtx))
                                    {
                                      fillHisto("hData_Hemi_Vtx_r","Ping_TightWP",      thesample,tree_Hemi_Vtx_r->at(iVtx),1.); 
                                      // fillHisto("hData_Hemi_Vtx_eta","Ping_TightWP",    thesample,tree_Hemi_Vtx_eta->at(iVtx),1.); 
                                      fillHisto("hData_Hemi_Vtx_dist","Ping_TightWP",   thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);                   
                                      float r = sqrt((tree_Hemi_LLP_x->at(iVtx)*tree_Hemi_LLP_x->at(iVtx))+(tree_Hemi_LLP_y->at(iVtx)*tree_Hemi_LLP_y->at(iVtx)));               
                                      fillHisto("hSim_Hemi_Vtx_r","Ping_TightWP",      thesample,r,1); 
                                      fillHisto("hSim_Hemi_Vtx_eta","Ping_TightWP",    thesample,tree_Hemi_LLP_eta->at(iVtx),1); 
                                      fillHisto("hSim_Hemi_Vtx_dist","Ping_TightWP",   thesample,tree_Hemi_LLP_dist->at(iVtx),1);
                                    }
                                }

                              nRecoVertexTightWP++;
                              MeanDistanceTightWP +=  tree_Hemi_Vtx_dist->at(iVtx);    

                              if ( Signal)
                              {        
                                int Msmu = tree_smu_mass;
                                int Mneu = tree_neu_mass;
                                if ( Msmu > 180 && Msmu < 220 )      Msmu = 200;
                                else if ( Msmu > 230 && Msmu < 270 ) Msmu = 250;
                                else if ( Msmu > 280 && Msmu < 320 ) Msmu = 300;
                                else if ( Msmu > 380 && Msmu < 420 ) Msmu = 400;
                                else if ( Msmu > 480 && Msmu < 520 ) Msmu = 500;
                                // else cout << " !!! smu mass out of range !!! " << Msmu;
                                if ( Mneu > 170 && Mneu < 190 )      Mneu = 180;
                                else if ( Mneu > 190 && Mneu < 210 ) Mneu = 200;
                                else if ( Mneu > 240 && Mneu < 260 ) Mneu = 250;
                                else if ( Mneu > 270 && Mneu < 290 ) Mneu = 280;
                                else if ( Mneu > 290 && Mneu < 310 ) Mneu = 300;
                                else if ( Mneu > 340 && Mneu < 360 ) Mneu = 350;
                                else if ( Mneu > 370 && Mneu < 390 ) Mneu = 380;
                                else if ( Mneu > 390 && Mneu < 410 ) Mneu = 400;
                                else if ( Mneu > 440 && Mneu < 460 ) Mneu = 450;
                                else if ( Mneu > 470 && Mneu < 490 ) Mneu = 480;
                                // else cout << " !!! neu mass out of range !!! " << Mneu;
     
                                if ((Msmu-Mneu) == 20){ fillHisto("hSim_VTXBDT","dm20",  thesample,tree_Hemi_Vtx_MVAval_Step1->at(iVtx) ,1.);fillHisto("hBDT_nTrks",     "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);fillHisto("hBDT_NChi2",     "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);fillHisto("hBDT_step",      "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);fillHisto("hBDT_r",         "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);fillHisto("hBDT_z",         "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);fillHisto("hBDT_MWT",       "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);fillHisto("hBDT_Vtx_Mass",  "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);fillHisto("hBDT_Hemi_Mass", "VTX_TightWP_dm20", thesample,tree_Hemi_mass->at(iVtx),1.);fillHisto("hBDT_dist",      "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);fillHisto("hBDT_ntrk10",    "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);fillHisto("hBDT_ntrk20",    "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);fillHisto("hBDT_MeanDCA",   "VTX_TightWP_dm20", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);}
                                if ((Msmu-Mneu) == 50){ fillHisto("hSim_VTXBDT","dm50",  thesample,tree_Hemi_Vtx_MVAval_Step1->at(iVtx) ,1.);fillHisto("hBDT_nTrks",     "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);fillHisto("hBDT_NChi2",     "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);fillHisto("hBDT_step",      "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);fillHisto("hBDT_r",         "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);fillHisto("hBDT_z",         "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);fillHisto("hBDT_MWT",       "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);fillHisto("hBDT_Vtx_Mass",  "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);fillHisto("hBDT_Hemi_Mass", "VTX_TightWP_dm50", thesample,tree_Hemi_mass->at(iVtx),1.);fillHisto("hBDT_dist",      "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);fillHisto("hBDT_ntrk10",    "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);fillHisto("hBDT_ntrk20",    "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);fillHisto("hBDT_MeanDCA",   "VTX_TightWP_dm50", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);}
                                if ((Msmu-Mneu) == 100){ fillHisto("hSim_VTXBDT","dm100",  thesample,tree_Hemi_Vtx_MVAval_Step1->at(iVtx) ,1.);fillHisto("hBDT_nTrks",     "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);fillHisto("hBDT_NChi2",     "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);fillHisto("hBDT_step",      "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);fillHisto("hBDT_r",         "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);fillHisto("hBDT_z",         "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);fillHisto("hBDT_MWT",       "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);fillHisto("hBDT_Vtx_Mass",  "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);fillHisto("hBDT_Hemi_Mass", "VTX_TightWP_dm100", thesample,tree_Hemi_mass->at(iVtx),1.);fillHisto("hBDT_dist",      "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);fillHisto("hBDT_ntrk10",    "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);fillHisto("hBDT_ntrk20",    "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);fillHisto("hBDT_MeanDCA",   "VTX_TightWP_dm100", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);}
                                if ((Msmu-Mneu) == 150){ fillHisto("hSim_VTXBDT","dm150",  thesample,tree_Hemi_Vtx_MVAval_Step1->at(iVtx) ,1.);fillHisto("hBDT_nTrks",     "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);fillHisto("hBDT_NChi2",     "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);fillHisto("hBDT_step",      "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);fillHisto("hBDT_r",         "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);fillHisto("hBDT_z",         "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);fillHisto("hBDT_MWT",       "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);fillHisto("hBDT_Vtx_Mass",  "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);fillHisto("hBDT_Hemi_Mass", "VTX_TightWP_dm150", thesample,tree_Hemi_mass->at(iVtx),1.);fillHisto("hBDT_dist",      "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);fillHisto("hBDT_ntrk10",    "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);fillHisto("hBDT_ntrk20",    "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);fillHisto("hBDT_MeanDCA",   "VTX_TightWP_dm150", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);}
                                if ((Msmu-Mneu) == 200){ fillHisto("hSim_VTXBDT","dm200",  thesample,tree_Hemi_Vtx_MVAval_Step1->at(iVtx) ,1.);fillHisto("hBDT_nTrks",     "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);fillHisto("hBDT_NChi2",     "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);fillHisto("hBDT_step",      "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);fillHisto("hBDT_r",         "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);fillHisto("hBDT_z",         "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);fillHisto("hBDT_MWT",       "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);fillHisto("hBDT_Vtx_Mass",  "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);fillHisto("hBDT_Hemi_Mass", "VTX_TightWP_dm200", thesample,tree_Hemi_mass->at(iVtx),1.);fillHisto("hBDT_dist",      "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);fillHisto("hBDT_ntrk10",    "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);fillHisto("hBDT_ntrk20",    "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);fillHisto("hBDT_MeanDCA",   "VTX_TightWP_dm200", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);}
                                if ((Msmu-Mneu) == 250){ fillHisto("hSim_VTXBDT","dm250",  thesample,tree_Hemi_Vtx_MVAval_Step1->at(iVtx) ,1.);fillHisto("hBDT_nTrks",     "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);fillHisto("hBDT_NChi2",     "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);fillHisto("hBDT_step",      "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);fillHisto("hBDT_r",         "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);fillHisto("hBDT_z",         "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);fillHisto("hBDT_MWT",       "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);fillHisto("hBDT_Vtx_Mass",  "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);fillHisto("hBDT_Hemi_Mass", "VTX_TightWP_dm250", thesample,tree_Hemi_mass->at(iVtx),1.);fillHisto("hBDT_dist",      "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);fillHisto("hBDT_ntrk10",    "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);fillHisto("hBDT_ntrk20",    "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);fillHisto("hBDT_MeanDCA",   "VTX_TightWP_dm250", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);}
                                if ((Msmu-Mneu) == 300){ fillHisto("hSim_VTXBDT","dm300",  thesample,tree_Hemi_Vtx_MVAval_Step1->at(iVtx) ,1.);fillHisto("hBDT_nTrks",     "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);fillHisto("hBDT_NChi2",     "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);fillHisto("hBDT_step",      "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);fillHisto("hBDT_r",         "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);fillHisto("hBDT_z",         "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);fillHisto("hBDT_MWT",       "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);fillHisto("hBDT_Vtx_Mass",  "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);fillHisto("hBDT_Hemi_Mass", "VTX_TightWP_dm300", thesample,tree_Hemi_mass->at(iVtx),1.);fillHisto("hBDT_dist",      "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);fillHisto("hBDT_ntrk10",    "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);fillHisto("hBDT_ntrk20",    "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);fillHisto("hBDT_MeanDCA",   "VTX_TightWP_dm300", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);}
                                if ((Msmu-Mneu) == 320){ fillHisto("hSim_VTXBDT","dm320",  thesample,tree_Hemi_Vtx_MVAval_Step1->at(iVtx) ,1.);fillHisto("hBDT_nTrks",     "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);fillHisto("hBDT_NChi2",     "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);fillHisto("hBDT_step",      "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);fillHisto("hBDT_r",         "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);fillHisto("hBDT_z",         "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);fillHisto("hBDT_MWT",       "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);fillHisto("hBDT_Vtx_Mass",  "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);fillHisto("hBDT_Hemi_Mass", "VTX_TightWP_dm320", thesample,tree_Hemi_mass->at(iVtx),1.);fillHisto("hBDT_dist",      "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);fillHisto("hBDT_ntrk10",    "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);fillHisto("hBDT_ntrk20",    "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);fillHisto("hBDT_MeanDCA",   "VTX_TightWP_dm320", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);}

                                int ngenpart =  tree_genParticle_pt->size();
                                for (int i=0; i<ngenpart; i++)    // Loop on GenParticle
                                  {
                                    float pdgId = tree_genParticle_pdgId->at(i); 
                                    float mother_pdgId = tree_genParticle_mother_pdgId->at(i); 
                                    float ct0 = tree_genParticle_ct0->at(i);

                                    // top quark from neutralino
                                    if ( abs(pdgId) == 6 && abs(mother_pdgId) == 1000023 ) 
                                      {            
                                        if ((Msmu-Mneu) == 20){fillHisto("hSim_ct0","dm20_VTXBDT",  thesample,ct0 ,1.);}
                                        if ((Msmu-Mneu) == 50){fillHisto("hSim_ct0","dm50_VTXBDT",  thesample,ct0 ,1.);}
                                        if ((Msmu-Mneu) == 100){fillHisto("hSim_ct0","dm100_VTXBDT",  thesample,ct0 ,1.);}
                                        if ((Msmu-Mneu) == 150){fillHisto("hSim_ct0","dm150_VTXBDT",  thesample,ct0 ,1.);}
                                        if ((Msmu-Mneu) == 200){fillHisto("hSim_ct0","dm200_VTXBDT",  thesample,ct0 ,1.);}
                                        if ((Msmu-Mneu) == 250){fillHisto("hSim_ct0","dm250_VTXBDT",  thesample,ct0 ,1.);}
                                        if ((Msmu-Mneu) == 300){fillHisto("hSim_ct0","dm300_VTXBDT",  thesample,ct0 ,1.);}
                                        if ((Msmu-Mneu) == 320){fillHisto("hSim_ct0","dm320_VTXBDT",  thesample,ct0 ,1.);}
                                      } 

                                  }    // End Loop on GenParticle  
                              }
                            fillHisto("hBDT_nTrks",     "VTX_TightWP", thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                            fillHisto("hBDT_NChi2",     "VTX_TightWP", thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                            fillHisto("hBDT_step",      "VTX_TightWP", thesample,tree_Hemi_Vtx_step->at(iVtx),1.);
                            fillHisto("hBDT_r",         "VTX_TightWP", thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                            fillHisto("hBDT_z",         "VTX_TightWP", thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                            fillHisto("hBDT_MWT",       "VTX_TightWP", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),1.);
                            fillHisto2D("MWTVsnTrks",   "VTX_TightWP", thesample,tree_Hemi_Vtx_MeantrackWeight->at(iVtx),tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                            fillHisto("hBDT_Vtx_Mass",  "VTX_TightWP", thesample,tree_Hemi_Vtx_Mass->at(iVtx),1.);
                            fillHisto("hBDT_Hemi_Mass", "VTX_TightWP", thesample,tree_Hemi_mass->at(iVtx),1.);
                            fillHisto("hBDT_dist",      "VTX_TightWP", thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                            fillHisto("hBDT_ntrk10",    "VTX_TightWP", thesample,tree_Hemi_Vtx_ntrk10->at(iVtx),1.);
                            fillHisto("hBDT_ntrk20",    "VTX_TightWP", thesample,tree_Hemi_Vtx_ntrk20->at(iVtx),1.);
                            fillHisto("hBDT_MeanDCA",   "VTX_TightWP", thesample,tree_Hemi_Vtx_track_MeanDCA_d->at(iVtx),1.);
    
                           }
                          
                     }// End Tight WP

                  // -------------- Loose WP ------------------------//
                  if (tree_Hemi_Vtx_step->at(iVtx)==3 || tree_Hemi_Vtx_step->at(iVtx)==4)
                     {
                        fillHisto("hData_Hemi_Vtx_nTrks","RecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_NChi2","RecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_dist","RecoVtx_LooseWP",    thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                        // fillHisto("hData_Hemi_Vtx_eta","RecoVtx_LooseWP",     thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_r","RecoVtx_LooseWP",       thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_z","RecoVtx_LooseWP",       thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                        fillHisto("hData_Hemi_Vtx_MVAval","RecoVtx_LooseWP",  thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);

                        
                        if ( tree_Hemi_Vtx_NChi2->at(iVtx)>0 && tree_Hemi_Vtx_NChi2->at(iVtx)<10)//Reco Vtx criteria
                           {
                              fillHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_NChi2->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx_LooseWP",   thesample,tree_Hemi_Vtx_nTrks->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_dist","GoodRecoVtx_LooseWP",    thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);
                              // fillHisto("hData_Hemi_Vtx_eta","GoodRecoVtx_LooseWP",     thesample,tree_Hemi_Vtx_eta->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_r","GoodRecoVtx_LooseWP",       thesample,tree_Hemi_Vtx_r->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_z","GoodRecoVtx_LooseWP",       thesample,tree_Hemi_Vtx_z->at(iVtx),1.);
                              fillHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx_LooseWP",  thesample,tree_Hemi_Vtx_MVAval->at(iVtx),1.);
                              if (Signal)
                                {
                                  if (tree_Hemi_LLP_ping->at(iVtx))
                                    {
                                      fillHisto("hData_Hemi_Vtx_r","Ping_LooseWP",      thesample,tree_Hemi_Vtx_r->at(iVtx),1.); 
                                      // fillHisto("hData_Hemi_Vtx_eta","Ping_LooseWP",    thesample,tree_Hemi_Vtx_eta->at(iVtx),1.); 
                                      fillHisto("hData_Hemi_Vtx_dist","Ping_LooseWP",   thesample,tree_Hemi_Vtx_dist->at(iVtx),1.);                   
                                      float r = sqrt((tree_Hemi_LLP_x->at(iVtx)*tree_Hemi_LLP_x->at(iVtx))+(tree_Hemi_LLP_y->at(iVtx)*tree_Hemi_LLP_y->at(iVtx)));               
                                      fillHisto("hSim_Hemi_Vtx_r","Ping_LooseWP",      thesample,r,1); 
                                      fillHisto("hSim_Hemi_Vtx_eta","Ping_LooseWP",    thesample,tree_Hemi_LLP_eta->at(iVtx),1); 
                                      fillHisto("hSim_Hemi_Vtx_dist","Ping_LooseWP",   thesample,tree_Hemi_LLP_dist->at(iVtx),1);
                                    }
             
                                }
                                  nRecoVertexLooseWP++;
                                  MeanDistanceLooseWP +=  tree_Hemi_Vtx_dist->at(iVtx); 
                           }

                     }//Loose WP
              }//RecoVtx


         }//End Loop on Vertices

        if (nTightVertex==2){nEvts_w2TightVtx++;}
        if (nTightVertex==1){nEvts_w1TightVtx++;}
        if (nSelecBDTVtxTight_evt == 2){nEvts_w2TightBDTVtx++;}
        if (nSelecBDTVtxTight_evt == 1){nEvts_w1TightBDTVtx++;}
        if (nSelecBDTVtxTight_step12_evt == 2){ nEvts_w2TightBDTVtx_step12++;fillHisto("DiMuon_Mass","2TightVtx", thesample,  tree_Mmumu, XS/(nentries));} //!!!! 
        if (nSelecBDTVtxTight_step12_evt == 1){ nEvts_w1TightBDTVtx_step12++;fillHisto("DiMuon_Mass","1TightVtx", thesample,  tree_Mmumu, XS/(nentries));} //!!!!
        if (countReco > 0 ){fillHisto("DiMuon_Mass","RecoVtx", thesample,  tree_Mmumu, XS/(nentries));} //!!!!
        if (countRecoGood >0 ){fillHisto("DiMuon_Mass","RecoVtx_Good_BDT3", thesample,  tree_Mmumu, XS/(nentries));};//!!!!
        if (countRecoTight12 >0 ){fillHisto("DiMuon_Mass","RecoVtx_TightWP", thesample,  tree_Mmumu, XS/(nentries)); }
        if (countRecoTight12_BDT3>0){fillHisto("DiMuon_Mass","RecoVtx_TightWP_BDT3", thesample,  tree_Mmumu, XS/(nentries));}//!!!!
        if (countRecoTight12_BDT3bis>0){fillHisto("DiMuon_Mass","RecoVtx_TightWP_BDT3bis", thesample,  tree_Mmumu, XS/(nentries)); }//!!!!
        // std::cout<<"here3"<<std::endl;
      //*****************************************
      // Resolution on the reconstructed vertices
      //*****************************************

      	// in PXB:
	// ttree->Draw("tree_SecInt_LLP_dr","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>2.6&&tree_SecInt_r<20&&abs(tree_SecInt_z)<27&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dr>0.12")
	// in TIB:
	// ttree->Draw("tree_SecInt_LLP_dr","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>22&&tree_SecInt_r<52&&abs(tree_SecInt_z)<68&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dr>0.8")
	// in TOB:
	// ttree->Draw("tree_SecInt_LLP_dr","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>50&&tree_SecInt_r<60&&abs(tree_SecInt_z)<108&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dr>3")
	// in PXF:
	// ttree->Draw("tree_SecInt_LLP_dz","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>4.5&&tree_SecInt_r<16.5&&abs(tree_SecInt_z)>30&&abs(tree_SecInt_z)<50&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dz>1")
	// in TID:
	// ttree->Draw("tree_SecInt_LLP_dz","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>22&&tree_SecInt_r<52&&abs(tree_SecInt_z)>70&&abs(tree_SecInt_z)<110&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dz>2.5")
	// in TEC:
	// ttree->Draw("tree_SecInt_LLP_dz","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>20&&tree_SecInt_r<75&&abs(tree_SecInt_z)>120&&abs(tree_SecInt_z)<200&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dz>5")

  for (unsigned int i = 0; i < tree_SecInt_x->size(); i++)
    {
      // res in the transverse plane
      if (tree_SecInt_LLP->at(i)>0 && tree_SecInt_selec->at(i) && (tree_SecInt_LLP_dd->at(i)/tree_SecInt_d->at(i))<0.1)
        {
          fillHisto("TranseverseResolution","NoSel", thesample,  tree_SecInt_LLP_dr->at(i),1);
          fillHisto("TranseverseDistance","NoSel", thesample,  tree_SecInt_r->at(i),1);
          fillHisto2D("TransverseResVsDistance","NoSel",thesample,  tree_SecInt_r->at(i),tree_SecInt_LLP_dr->at(i),1);
        }
      // res along the z-axis
      if (tree_SecInt_LLP->at(i)>0 && tree_SecInt_selec->at(i) && (tree_SecInt_LLP_dd->at(i)/tree_SecInt_d->at(i))<0.1)
        {
          fillHisto("LongitudinalResolution","NoSel", thesample,  tree_SecInt_LLP_dz->at(i),1);
          fillHisto("LongitudinalDistance","NoSel", thesample,  tree_SecInt_z->at(i),1);
          fillHisto2D("LongitudinalResVsDistance","NoSel",thesample,  tree_SecInt_z->at(i),tree_SecInt_LLP_dz->at(i),1);

        }
            // res along the z-axis
      if (tree_SecInt_LLP->at(i)>0 && tree_SecInt_selec->at(i) && (tree_SecInt_LLP_dd->at(i)/tree_SecInt_d->at(i))<0.1)
        {
          fillHisto("3DResolution","NoSel", thesample,  tree_SecInt_LLP_dd->at(i),1);
          fillHisto("3DDistance","NoSel", thesample,  tree_SecInt_d->at(i),1);
          fillHisto2D("3DResVsDistance","NoSel",thesample,  tree_SecInt_d->at(i),tree_SecInt_LLP_dd->at(i),1);
        }
      
    }

         //---------------Vtx Selection Variables----------------//
  //------------------MWT CUT -------------------//
  for (int i =0 ;i<nSteps+1;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_MeantrackWeight->at(j)>(MeanTWcut+i*dTW) && tree_Hemi_Vtx_nTrks->at(j)==2 )
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

  // for (int i =0 ;i<nStep_dd;i++)
  //   {
  //     for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
  //       {
  //         if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_Vtx_Vtx_dd->at(j)>=(VtxddCut+i*ddd))
  //           {
  //             nEvts_Vtx_dd[i]++;
  //             if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2 )
  //               {
  //                 nEvts_Vtx_dd_step1[i]++;
  //               }
  //           }
  //       }
  //   }
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
      for (unsigned int j = 0 ; j< tree_event_nVtx->size();j++)
        {
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_event_nVtx->at(j)==(VtxnVtxCut+i*dnVtx))
            {
              nEvts_Vtx_nVtx[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2 )
                {
                  nEvts_Vtx_nVtx_step1[i]++;
                }
            }
        }
    }

      for (unsigned int j = 0 ; j< tree_event_nVtx->size();j++)
        {
          if (tree_event_nVtx->at(j)==2)
            {
              nEvts_w2Vtx++;
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
          if (tree_Hemi_Vtx_NChi2->at(j)>0 && tree_Hemi_Vtx_NChi2->at(j)<10 && tree_Hemi_mass->at(j)>=(HemiMassCut+i*dHemiMass))
            {
              nEvts_Hemi_Mass[i]++;
              if (tree_Hemi_Vtx_step->at(j)==1 || tree_Hemi_Vtx_step->at(j)==2)
                {
                  nEvts_Hemi_Mass_step1[i]++;
                }
            }
        }
    }

//--------Smuon_InvMass------
  for (int i =0 ;i<nStep_SmuonMass;i++)
    {
      for (unsigned int j = 0 ; j< tree_Hemi->size();j++)
        {
          if (tree_Hemi_mass->at(j)>=(Hemi_CombinedHemiLeptonMassCut+i*dSmuonMass))
            {
              nEvts_Hemi_CombinedHemiLeptonMass[i]++;
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
    for (unsigned int j = 0 ; j < tree_lepton_leadingpt->size() ; j++)
      {
        if (tree_lepton_leadingpt->at(j) > (MuonPTCut+i*dpt_muon))
          {
            nEvts_muonpt[i]++;
          }
      }
  }

//--------subLeading muon pt------
for (int i = 0 ; i < nStep_muon2_pt ; i++)
  {
    for (unsigned int j = 0 ; j < tree_lepton_leadingpt2->size() ; j++)
      {
        if (tree_lepton_leadingpt2->at(j) >= (Muon2PTCut+i*dpt_muon2))
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
    for (unsigned int j = 0 ; j < tree_lepton_leadingpt->size() ; j++)
      {
        if (tree_LT > (STCut+i*dpt_ST))
          {
            nEvts_ST[i]++;
          }
      }
  }

//--------nJets
for (int i = 0 ; i < nStep_nJet ; i++)
  {
    for (unsigned int j = 0 ; j < tree_lepton_leadingpt->size() ; j++)
      {
        if (tree_njet > (nJetCut+i*dnJet_ST))
          {
            nEvts_nJet[i]++;
          }
      }
  }

  //-------------nmu----------//

// for (int i = 0 ; i < nStep_nMuon ; i++)
//   {
//     for (unsigned int j = 0 ; j < tree_muon_pt->size() ; j++)
//       {
//         if (tree_muon_nmu->at(j) > (nMuonCut+i*dnMuon))
//           {
//             nEvts_nMuon[i]++;
//           }
//       }
//   }
// //------------------------------

// //-------------muon_muon_dR----------//

for (int i = 0 ; i < nStep_MuondR ; i++)
  {
    for (unsigned int j = 0 ; j < tree_lepton_lepton_dR->size() ; j++)
      {
        if (tree_lepton_lepton_dR->at(j) > (MuondRCut+i*dMuondR))
          {
            nEvts_MuondR[i]++;
          }
      }
  }
// //------------------------------

// //-------------muon_muon_dPhi----------//

for (int i = 0 ; i < nStep_MuondPhi ; i++)
  {
    for (unsigned int j = 0 ; j < tree_lepton_lepton_dPhi->size() ; j++)
      {
        if (tree_lepton_lepton_dPhi->at(j) > (MuondPhiCut+i*dMuondPhi))
          {
            nEvts_MuondPhi[i]++;
          }
      }
  }
// //------------------------------

// //-------------muon_muon_dEta----------//

for (int i = 0 ; i < nStep_MuondEta ; i++)
  {
    for (unsigned int j = 0 ; j < tree_lepton_lepton_dEta->size() ; j++)
      {
        if (tree_lepton_lepton_dEta->at(j) > (MuondEtaCut+i*dMuondEta ))
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


  for (int i =0 ;i<nSteps_BDT;i++)
    {

      for (unsigned int j = 0 ; j< tree_lepton_leadingpt->size();j++)
        {
          if ( tree_Evts_MVAval>(MeanBDTcut+i*dBDT) )
            {
              nSelecEvts[i]++;
            }
        }
    }
// std::cout<<"debugB"<<std::endl;

      //*******************************
      // Loop on Gen Particles
      //*******************************

////////                                                               ////////
////////                          GENERATION                           ////////
////////                                                               ////////

  if (Signal)
    {    
          
      int Msmu = tree_smu_mass;
      int Mneu = tree_neu_mass;
      if ( Msmu > 180 && Msmu < 220 )      Msmu = 200;
      else if ( Msmu > 230 && Msmu < 270 ) Msmu = 250;
      else if ( Msmu > 280 && Msmu < 320 ) Msmu = 300;
      else if ( Msmu > 380 && Msmu < 420 ) Msmu = 400;
      else if ( Msmu > 480 && Msmu < 520 ) Msmu = 500;
      // else cout << " !!! smu mass out of range !!! " << Msmu;
      if ( Mneu > 170 && Mneu < 190 )      Mneu = 180;
      else if ( Mneu > 190 && Mneu < 210 ) Mneu = 200;
      else if ( Mneu > 240 && Mneu < 260 ) Mneu = 250;
      else if ( Mneu > 270 && Mneu < 290 ) Mneu = 280;
      else if ( Mneu > 290 && Mneu < 310 ) Mneu = 300;
      else if ( Mneu > 340 && Mneu < 360 ) Mneu = 350;
      else if ( Mneu > 370 && Mneu < 390 ) Mneu = 380;
      else if ( Mneu > 390 && Mneu < 410 ) Mneu = 400;
      else if ( Mneu > 440 && Mneu < 460 ) Mneu = 450;
      else if ( Mneu > 470 && Mneu < 490 ) Mneu = 480;
      // else cout << " !!! neu mass out of range !!! " << Mneu;
        

      // fillHisto("hGen_Msmu","noSel",  thesample, tree_smu_mass ,1.);
      // fillHisto("hGen_Mneu","noSel",  thesample, tree_neu_mass ,1.);

      int ngenpart =  tree_genParticle_pt->size();
      for (int i=0; i<ngenpart; i++)    // Loop on GenParticle
        {
          float pdgId = tree_genParticle_pdgId->at(i); 
          float mother_pdgId = tree_genParticle_mother_pdgId->at(i); 
          float ct0 = tree_genParticle_ct0->at(i);

          // top quark from neutralino
          if ( abs(pdgId) == 6 && abs(mother_pdgId) == 1000023 ) 
            {
              fillHisto("hGen_ct0","noSel",  thesample, ct0 ,1.);
              
              if ( Msmu == 200 && Mneu == 180 ) fillHisto("hGen_ct0","smu200_neu180",  thesample, ct0 ,1.);
              if ( Msmu == 250 && Mneu == 200 ) fillHisto("hGen_ct0","smu250_neu200",  thesample, ct0 ,1.);
              if ( Msmu == 300 && Mneu == 180 ) fillHisto("hGen_ct0","smu300_neu180",  thesample, ct0 ,1.);
              if ( Msmu == 300 && Mneu == 200 ) fillHisto("hGen_ct0","smu300_neu200",  thesample, ct0 ,1.);
              if ( Msmu == 300 && Mneu == 250 ) fillHisto("hGen_ct0","smu300_neu250",  thesample, ct0 ,1.);
              if ( Msmu == 300 && Mneu == 280 ) fillHisto("hGen_ct0","smu300_neu280",  thesample, ct0 ,1.);
              if ( Msmu == 400 && Mneu == 180 ) fillHisto("hGen_ct0","smu400_neu180",  thesample, ct0 ,1.);
              if ( Msmu == 400 && Mneu == 200 ) fillHisto("hGen_ct0","smu400_neu200",  thesample, ct0 ,1.);
              if ( Msmu == 400 && Mneu == 250 ) fillHisto("hGen_ct0","smu400_neu250",  thesample, ct0 ,1.);
              if ( Msmu == 400 && Mneu == 300 ) fillHisto("hGen_ct0","smu400_neu300",  thesample, ct0 ,1.);
              if ( Msmu == 400 && Mneu == 350 ) fillHisto("hGen_ct0","smu400_neu350",  thesample, ct0 ,1.);
              if ( Msmu == 400 && Mneu == 380 ) fillHisto("hGen_ct0","smu400_neu380",  thesample, ct0 ,1.);
              if ( Msmu == 500 && Mneu == 180 ) fillHisto("hGen_ct0","smu500_neu180",  thesample, ct0 ,1.);
              if ( Msmu == 500 && Mneu == 200 ) fillHisto("hGen_ct0","smu500_neu200",  thesample, ct0 ,1.);
              if ( Msmu == 500 && Mneu == 250 ) fillHisto("hGen_ct0","smu500_neu250",  thesample, ct0 ,1.);
              if ( Msmu == 500 && Mneu == 300 ) fillHisto("hGen_ct0","smu500_neu300",  thesample, ct0 ,1.);
              if ( Msmu == 500 && Mneu == 350 ) fillHisto("hGen_ct0","smu500_neu350",  thesample, ct0 ,1.);
              if ( Msmu == 500 && Mneu == 400 ) fillHisto("hGen_ct0","smu500_neu400",  thesample, ct0 ,1.);
              if ( Msmu == 500 && Mneu == 450 ) fillHisto("hGen_ct0","smu500_neu450",  thesample, ct0 ,1.);
              if ( Msmu == 500 && Mneu == 480 ) fillHisto("hGen_ct0","smu500_neu480",  thesample, ct0 ,1.);
            } 

        }    // End Loop on GenParticle

////////                                                               ////////
////////                         PRESELECTION                          ////////
////////                                                               ////////

      bool Filter = tree_Filter; // Loose ID and Loose Iso for both muons
      //$$   if ( tree_GoodMu1 != 11 || tree_GoodMu2 != 11 ) Filter = false; // Tight ID and Tight Iso for both muons

      if ( Msmu == 200 && Mneu == 180 ) fillHisto("hSim_Filter","smu200_neu180",  thesample, Filter ,1.);
      if ( Msmu == 250 && Mneu == 200 ) fillHisto("hSim_Filter","smu250_neu200",  thesample, Filter ,1.);
      if ( Msmu == 300 && Mneu == 180 ) fillHisto("hSim_Filter","smu300_neu180",  thesample, Filter ,1.);
      if ( Msmu == 300 && Mneu == 200 ) fillHisto("hSim_Filter","smu300_neu200",  thesample, Filter ,1.);
      if ( Msmu == 300 && Mneu == 250 ) fillHisto("hSim_Filter","smu300_neu250",  thesample, Filter ,1.);
      if ( Msmu == 300 && Mneu == 280 ) fillHisto("hSim_Filter","smu300_neu280",  thesample, Filter ,1.);
      if ( Msmu == 400 && Mneu == 180 ) fillHisto("hSim_Filter","smu400_neu180",  thesample, Filter ,1.);
      if ( Msmu == 400 && Mneu == 200 ) fillHisto("hSim_Filter","smu400_neu200",  thesample, Filter ,1.);
      if ( Msmu == 400 && Mneu == 250 ) fillHisto("hSim_Filter","smu400_neu250",  thesample, Filter ,1.);
      if ( Msmu == 400 && Mneu == 300 ) fillHisto("hSim_Filter","smu400_neu300",  thesample, Filter ,1.);
      if ( Msmu == 400 && Mneu == 350 ) fillHisto("hSim_Filter","smu400_neu350",  thesample, Filter ,1.);
      if ( Msmu == 400 && Mneu == 380 ) fillHisto("hSim_Filter","smu400_neu380",  thesample, Filter ,1.);
      if ( Msmu == 500 && Mneu == 180 ) fillHisto("hSim_Filter","smu500_neu180",  thesample, Filter ,1.);
      if ( Msmu == 500 && Mneu == 200 ) fillHisto("hSim_Filter","smu500_neu200",  thesample, Filter ,1.);
      if ( Msmu == 500 && Mneu == 250 ) fillHisto("hSim_Filter","smu500_neu250",  thesample, Filter ,1.);
      if ( Msmu == 500 && Mneu == 300 ) fillHisto("hSim_Filter","smu500_neu300",  thesample, Filter ,1.);
      if ( Msmu == 500 && Mneu == 350 ) fillHisto("hSim_Filter","smu500_neu350",  thesample, Filter ,1.);
      if ( Msmu == 500 && Mneu == 400 ) fillHisto("hSim_Filter","smu500_neu400",  thesample, Filter ,1.);
      if ( Msmu == 500 && Mneu == 450 ) fillHisto("hSim_Filter","smu500_neu450",  thesample, Filter ,1.);
      if ( Msmu == 500 && Mneu == 480 ) fillHisto("hSim_Filter","smu500_neu480",  thesample, Filter ,1.);

      //$$
      if ( !Filter ) continue; 
      //$$


///////////////////////
// Delta R between hemisphere axis and closest neutralino

   int nHemi = tree_Hemi_LLP_dR->size();
   for (int i=0; i<nHemi; i++) {   // Loop on hemispheres
     float dR = tree_Hemi_LLP_dR->at(i);
     
     if ( Msmu == 200 && Mneu == 180 ) fillHisto("hSim_Hemi_dR","smu200_neu180",  thesample, dR ,1.);
     if ( Msmu == 250 && Mneu == 200 ) fillHisto("hSim_Hemi_dR","smu250_neu200",  thesample, dR ,1.);
     if ( Msmu == 300 && Mneu == 180 ) fillHisto("hSim_Hemi_dR","smu300_neu180",  thesample, dR ,1.);
     if ( Msmu == 300 && Mneu == 200 ) fillHisto("hSim_Hemi_dR","smu300_neu200",  thesample, dR ,1.);
     if ( Msmu == 300 && Mneu == 250 ) fillHisto("hSim_Hemi_dR","smu300_neu250",  thesample, dR ,1.);
     if ( Msmu == 300 && Mneu == 280 ) fillHisto("hSim_Hemi_dR","smu300_neu280",  thesample, dR ,1.);
     if ( Msmu == 400 && Mneu == 180 ) fillHisto("hSim_Hemi_dR","smu300_neu180",  thesample, dR ,1.);
     if ( Msmu == 400 && Mneu == 200 ) fillHisto("hSim_Hemi_dR","smu400_neu200",  thesample, dR ,1.);
     if ( Msmu == 400 && Mneu == 250 ) fillHisto("hSim_Hemi_dR","smu400_neu250",  thesample, dR ,1.);
     if ( Msmu == 400 && Mneu == 300 ) fillHisto("hSim_Hemi_dR","smu400_neu300",  thesample, dR ,1.);
     if ( Msmu == 400 && Mneu == 350 ) fillHisto("hSim_Hemi_dR","smu400_neu350",  thesample, dR ,1.);
     if ( Msmu == 400 && Mneu == 380 ) fillHisto("hSim_Hemi_dR","smu400_neu380",  thesample, dR ,1.);
     if ( Msmu == 500 && Mneu == 180 ) fillHisto("hSim_Hemi_dR","smu500_neu180",  thesample, dR ,1.);
     if ( Msmu == 500 && Mneu == 200 ) fillHisto("hSim_Hemi_dR","smu500_neu200",  thesample, dR ,1.);
     if ( Msmu == 500 && Mneu == 250 ) fillHisto("hSim_Hemi_dR","smu500_neu250",  thesample, dR ,1.);
     if ( Msmu == 500 && Mneu == 300 ) fillHisto("hSim_Hemi_dR","smu500_neu300",  thesample, dR ,1.);
     if ( Msmu == 500 && Mneu == 350 ) fillHisto("hSim_Hemi_dR","smu500_neu350",  thesample, dR ,1.);
     if ( Msmu == 500 && Mneu == 400 ) fillHisto("hSim_Hemi_dR","smu500_neu400",  thesample, dR ,1.);
     if ( Msmu == 500 && Mneu == 450 ) fillHisto("hSim_Hemi_dR","smu500_neu450",  thesample, dR ,1.);
     if ( Msmu == 500 && Mneu == 480 ) fillHisto("hSim_Hemi_dR","smu500_neu480",  thesample, dR ,1.);

     float dist = tree_Hemi_Vtx_dist->at(i);
     float NChi2 = tree_Hemi_Vtx_NChi2->at(i);
     
   if ( NChi2 > 0. && NChi2 < 10. ) {

      
     if ( Msmu == 200 && Mneu == 180 ) fillHisto("hData_Hemi_Vtx_dist","smu200_neu180",  thesample, dist ,1.);
     if ( Msmu == 250 && Mneu == 200 ) fillHisto("hData_Hemi_Vtx_dist","smu250_neu200",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 180 ) fillHisto("hData_Hemi_Vtx_dist","smu300_neu180",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 200 ) fillHisto("hData_Hemi_Vtx_dist","smu300_neu200",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 250 ) fillHisto("hData_Hemi_Vtx_dist","smu300_neu250",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 280 ) fillHisto("hData_Hemi_Vtx_dist","smu300_neu280",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 180 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu180",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 200 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu200",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 250 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu250",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 300 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu300",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 350 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu350",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 380 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu380",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 180 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu180",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 200 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu200",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 250 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu250",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 300 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu300",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 350 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu350",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 400 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu400",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 450 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu450",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 480 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu480",  thesample, dist ,1.);
   }
   else {
     if ( Msmu == 200 && Mneu == 180 ) fillHisto("hData_Hemi_Vtx_dist","smu200_neu180",  thesample, -1 ,1.);
     if ( Msmu == 250 && Mneu == 200 ) fillHisto("hData_Hemi_Vtx_dist","smu250_neu200",  thesample, -1 ,1.);
     if ( Msmu == 300 && Mneu == 180 ) fillHisto("hData_Hemi_Vtx_dist","smu300_neu180",  thesample, -1 ,1.);
     if ( Msmu == 300 && Mneu == 200 ) fillHisto("hData_Hemi_Vtx_dist","smu300_neu200",  thesample, -1 ,1.);
     if ( Msmu == 300 && Mneu == 250 ) fillHisto("hData_Hemi_Vtx_dist","smu300_neu250",  thesample, -1 ,1.);
     if ( Msmu == 300 && Mneu == 280 ) fillHisto("hData_Hemi_Vtx_dist","smu300_neu280",  thesample, -1 ,1.);
     if ( Msmu == 400 && Mneu == 180 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu180",  thesample, -1 ,1.);
     if ( Msmu == 400 && Mneu == 200 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu200",  thesample, -1 ,1.);
     if ( Msmu == 400 && Mneu == 250 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu250",  thesample, -1 ,1.);
     if ( Msmu == 400 && Mneu == 300 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu300",  thesample, -1 ,1.);
     if ( Msmu == 400 && Mneu == 350 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu350",  thesample, -1 ,1.);
     if ( Msmu == 400 && Mneu == 380 ) fillHisto("hData_Hemi_Vtx_dist","smu400_neu380",  thesample, -1 ,1.);
     if ( Msmu == 500 && Mneu == 180 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu180",  thesample, -1 ,1.);
     if ( Msmu == 500 && Mneu == 200 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu200",  thesample, -1 ,1.);
     if ( Msmu == 500 && Mneu == 250 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu250",  thesample, -1 ,1.);
     if ( Msmu == 500 && Mneu == 300 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu300",  thesample, -1 ,1.);
     if ( Msmu == 500 && Mneu == 350 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu350",  thesample, -1 ,1.);
     if ( Msmu == 500 && Mneu == 400 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu400",  thesample, -1 ,1.);
     if ( Msmu == 500 && Mneu == 450 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu450",  thesample, -1 ,1.);
     if ( Msmu == 500 && Mneu == 480 ) fillHisto("hData_Hemi_Vtx_dist","smu500_neu480",  thesample, -1 ,1.);
   }
    
   if ( NChi2 > 0. && NChi2 < 10. && tree_Hemi_LLP_ping->at(i) ) {
     if ( Msmu == 200 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu200_neu180",  thesample, dist ,1.);
     if ( Msmu == 250 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu250_neu200",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu180",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu200",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu250",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 280 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu280",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu180",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu200",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu250",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 300 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu300",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 350 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu350",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 380 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu380",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu180",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu200",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu250",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 300 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu300",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 350 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu350",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 400 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu400",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 450 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu450",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 480 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu480",  thesample, dist ,1.);
   }
   else {
     if ( Msmu == 200 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu200_neu180",  thesample, -1 ,1.);
     if ( Msmu == 250 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu250_neu200",  thesample, -1 ,1.);
     if ( Msmu == 300 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu180",  thesample, -1,1.);
     if ( Msmu == 300 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu200",  thesample, -1,1.);
     if ( Msmu == 300 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu250",  thesample, -1,1.);
     if ( Msmu == 300 && Mneu == 280 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu280",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu180",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu200",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu250",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 300 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu300",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 350 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu350",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 380 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu380",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu180",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu200",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu250",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 300 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu300",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 350 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu350",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 400 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu400",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 450 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu450",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 480 ) fillHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu480",  thesample, -1,1.);

   }
     
   if ( NChi2 > 0. && NChi2 < 10. && tree_Hemi_LLP_ping->at(i) && (tree_Hemi_Vtx_step->at(i) == 1 || tree_Hemi_Vtx_step->at(i) == 2) ) {
     if ( Msmu == 200 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu200_neu180",  thesample, dist ,1.);
     if ( Msmu == 250 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu250_neu200",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu180",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu200",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu250",  thesample, dist ,1.);
     if ( Msmu == 300 && Mneu == 280 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu280",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu180",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu200",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu250",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 300 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu300",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 350 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu350",  thesample, dist ,1.);
     if ( Msmu == 400 && Mneu == 380 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu380",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu180",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu200",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu250",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 300 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu300",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 350 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu350",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 400 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu400",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 450 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu450",  thesample, dist ,1.);
     if ( Msmu == 500 && Mneu == 480 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu480",  thesample, dist ,1.);

   }
   else {
     if ( Msmu == 200 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu200_neu180",  thesample, -1,1.);
     if ( Msmu == 250 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu250_neu200",  thesample, -1,1.);
     if ( Msmu == 300 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu180",  thesample, -1,1.);
     if ( Msmu == 300 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu200",  thesample, -1,1.);
     if ( Msmu == 300 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu250",  thesample, -1,1.);
     if ( Msmu == 300 && Mneu == 280 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu280",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu180",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu200",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu250",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 300 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu300",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 350 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu350",  thesample, -1,1.);
     if ( Msmu == 400 && Mneu == 380 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu380",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 180 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu180",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 200 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu200",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 250 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu250",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 300 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu300",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 350 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu350",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 400 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu400",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 450 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu450",  thesample, -1,1.);
     if ( Msmu == 500 && Mneu == 480 ) fillHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu480",  thesample, -1,1.);
   }

   }


    } // end loop on sample 200to500
  }// End of Tree_FIlter
}// End Global Loop

// std::cout<<"here5"<<std::endl;
    //-----------------------------------//
//------------Evts Bdt selection
  for (int i =0 ;i<nSteps_BDT;i++)
    {
      fillHisto("hData_Evts","BDT",  thesample,MeanBDTcut+i*dBDT ,nSelecEvts[i]/nEvts);
    }

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
      fillHisto("hData_Vtx_nTrk20","TightWP",  thesample,ntrk20Cut+h*dntrk20 ,nEvts_ntrk20_step1[h]/(nRecoVertexTightWP));
  }
// std::cout<<"bug starts here 3"<<std::endl;
for (int h = 0 ; h < nStep_dd; h++)
  {
      fillHisto("Data_Vtx_Vtx_dd","NoSel",  thesample,VtxddCut+h*ddd ,nEvts_Vtx_dd[h]/nEvts_w2Vtx);
      fillHisto("Data_Vtx_Vtx_dd","TightWP",  thesample,VtxddCut+h*ddd ,nEvts_Vtx_dd_step1[h]/nRecoVertexTightWP);
  }

  for (int h = 0 ; h < nStep_nTrks; h++)// !!Bug!!
  {
      fillHisto("hData_Vtx_nTrks","NoSel",  thesample,VtxnTrksCut+h*dnTrks ,nEvts_Vtx_nTrks[h]/nRecoVertex);
      fillHisto("hData_Vtx_nTrks","TightWP",  thesample,VtxnTrksCut+h*dnTrks,nEvts_Vtx_nTrks_step1[h]/nRecoVertexTightWP);
  }

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
// std::cout<<"bug starts here 7"<<std::endl;
//--------Vtx_Step------
  for (int h = 0 ; h < nStep_Step; h++)
  {
      fillHisto("hData_Vtx_Step","NoSel",  thesample,VtxStepCut+h*dStep ,nEvts_Vtx_Step[h]/nRecoVertex);
      fillHisto("hData_Vtx_Step","TightWP",  thesample,VtxStepCut+h*dStep,nEvts_Vtx_Step_step1[h]/nRecoVertexTightWP);
  }
// std::cout<<"bug starts here 8"<<std::endl;
//--------Vtx_InvMass------
  for (int h = 0 ; h < nStep_Mass; h++)
  {
      fillHisto("hData_Vtx_Mass","NoSel",  thesample,VtxMassCut+h*dMass,nEvts_Vtx_Mass[h]/nRecoVertex);
      fillHisto("hData_Vtx_Mass","TightWP",  thesample,VtxMassCut+h*dMass ,nEvts_Vtx_Mass_step1[h]/nRecoVertexTightWP);
  }
// std::cout<<"bug starts here 9"<<std::endl;
//--------Hemi_InvMass------
  for (int h = 0 ; h < nStep_HemiMass; h++)
  {
      fillHisto("hData_Hemi_Mass","NoSel",  thesample,HemiMassCut+h*dHemiMass ,nEvts_Hemi_Mass[h]/nRecoVertex);
      fillHisto("hData_Hemi_Mass","TightWP",  thesample,HemiMassCut+h*dHemiMass ,nEvts_Hemi_Mass_step1[h]/nRecoVertexTightWP);
  }

//Smuon mass
  for (int h = 0 ; h < nStep_SmuonMass; h++)
  {
      fillHisto("hData_Hemi_SmuonMass","NoSel",  thesample,Hemi_CombinedHemiLeptonMassCut+h*dSmuonMass ,nEvts_Hemi_CombinedHemiLeptonMass[h]/(2*nEvts));
  }

//-----------------------------------------------
// std::cout<<"bug starts here 10"<<std::endl;
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
    fillHisto("hData_Evts_LeadingJet2Pt","NoSel",  thesample,Jet2PTCut+h*dpt,nEvts_jetpt2[h]/(nEvts));
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

   //---------------------------------------------
   float EvtSize = nentries;
   // Output Postscript
   ofs<<" //----------------- Efficacity for the Ntuple : " << thesample << " ---------------\\ "<<std::endl;
   ofs<<" //-Tracks selection efficiency by the Tight WP : " << TightTrks_Eff << " ---------------\\ "<<std::endl;
   ofs<<" //-Tracks selection efficiency by the Loose WP : " << LooseTrks_Eff << " ---------------\\ "<<std::endl;
   ofs<<"|| ----------------------------------------------------------------------------------"<<std::endl;
   ofs<<"|| Initial number of event :"<<nentries<<std::endl;
   ofs<<"|| NVertex To be reco (twice the number of evts passing Tree_Filter): "<<TotalnVertex<<" which means "<<TotalnVertex/(2*nentries)<<" efficiency of the filter"<<std::endl;
   ofs<<"MeanDistance :"<<MeanD<<" cm "<<std::endl;
   ofs<<"MeanDistanceTightWP :"<<MeanDTightWP<<" cm "<<std::endl;
   ofs<<"MeanDistanceLooseWP :"<<MeanDLooseWP<<" cm "<<std::endl;
   ofs<<" Total NRecoVertexEff (vertex reco over vertex to be reconstructed): "<<NRecoVertexEff<<std::endl;
   ofs<<" NRecoVertexEff TightWp: "<<100*nRecoVertexTightWP/nRecoVertex<<" per cent of "<<NRecoVertexEff<<std::endl;
   ofs<<"  NRecoVertexEff  LooseWP: "<<100*nRecoVertexLooseWP/nRecoVertex<<" per cent of "<<NRecoVertexEff<<std::endl;
   ofs<<" //----------------- End of Table ---------------\\ "<<std::endl; 
      ofs<<" "<<std::endl;
   ofs<<"Event yields for each step : "<<std::endl;
   ofs<<"Nevents: "<<nentries*NormFactor*lumiRun2<<std::endl;
   ofs<<"Online+Offline Selection : "<<nFilterEvt*NormFactor*lumiRun2<<std::endl;
   ofs<<"EVTS BDT for WP "<<EVTSWP <<" : "<<nEvts*NormFactor*lumiRun2<<std::endl;
   ofs<<"Vertex Reco: "<<nRecoVertex*(NormFactor/2)*lumiRun2<<std::endl;
   ofs<<"Vtx BDT3 > cut "<<VTXWP <<": "<<nSelecBDTVtx*(NormFactor/2)*lumiRun2<<std::endl;
   ofs<<"BDT2-Tight Vertex only : "<<nRecoVertexTightWP*(NormFactor/2)*lumiRun2<<std::endl;
   ofs<<"BDT2-Tight BDT3 Vtx > cut "<<VTXWP <<": "<<nSelecBDTVtxTight*(NormFactor/2)*lumiRun2<<std::endl;
   ofs<<"Two BDT2-Tight vertices : "<<nEvts_w2TightVtx*(NormFactor)*lumiRun2<<std::endl;
   ofs<<"Two BDT2-Tight BDT3 Vtx > cut "<<VTXWP <<": "<<nEvts_w2TightBDTVtx*(NormFactor)*lumiRun2<<std::endl;
   ofs<<"Two BDT2-Tight BDT3bis Vtx > cut "<<VTXWP <<": "<<nEvts_w2TightBDTVtx_step12*(NormFactor)*lumiRun2<<std::endl;
   ofs<<"One BDT2-Tight Vertex: "<<nEvts_w1TightVtx*(NormFactor)*lumiRun2<<std::endl;
   ofs<<"One BDT2-Tight BDT3 Vtx > cut "<<VTXWP <<": "<<nEvts_w1TightBDTVtx*(NormFactor)*lumiRun2<<std::endl;
   ofs<<"One BDT2-Tight BDT3bis Vtx > cut "<<VTXWP <<": "<<nEvts_w1TightBDTVtx_step12*(NormFactor)*lumiRun2<<std::endl;
   ofs<<" //----------------- End of Event yields ---------------\\ "<<std::endl;
   ofs<<" "<<std::endl;
   ofs<<"Selection Efficiency of each step : "<<std::endl;
   ofs<<"Nevents: "<<(nentries/nentries) <<std::endl;
   ofs<<"Online+Offline Selection : "<<(nFilterEvt/nentries) <<std::endl;
   ofs<<"EVTS BDT for WP "<<EVTSWP <<" : "<<(nEvts/nFilterEvt) <<std::endl;
   ofs<<"Vertex Reco: "<<(nRecoVertex/(2*nEvts)) <<std::endl;
   ofs<<"Vtx BDT3 > cut "<<VTXWP <<": "<<(nSelecBDTVtx/(nRecoVertex)) <<std::endl;
    ofs<<"BDT2-Tight Vertex only : "<<(nRecoVertexTightWP/(nRecoVertex)) <<std::endl;
   ofs<<"BDT2-Tight BDT3 Vtx > cut "<<VTXWP <<": "<<(nSelecBDTVtxTight/(nRecoVertexTightWP)) <<std::endl;
     ofs<<"Two BDT2-Tight vertices : "<<(2*nEvts_w2TightVtx/nRecoVertexTightWP)*(1) <<std::endl;
  //  ofs<<"Two BDT2-Tight BDT3 Vtx > cut "<<VTXWP <<": "<<(2*nEvts_w2TightBDTVtx/nRecoVertexTightWP)*(1) <<std::endl;
   ofs<<"Two BDT2-Tight BDT3bis Vtx > cut "<<VTXWP <<": "<<(2*nEvts_w2TightBDTVtx_step12/nRecoVertexTightWP)*(1) <<std::endl;
   ofs<<"One BDT2-Tight Vertex: "<<(2*nEvts_w1TightVtx/nRecoVertexTightWP)*(1) <<std::endl;
  //  ofs<<"One BDT2-Tight BDT3 Vtx > cut "<<VTXWP <<": "<<(2*nEvts_w1TightBDTVtx/nRecoVertexTightWP)*(1) <<std::endl;
   ofs<<"One BDT2-Tight BDT3bis Vtx > cut "<<VTXWP <<": "<<(2*nEvts_w1TightBDTVtx_step12/nRecoVertexTightWP)*(1) <<std::endl;
   ofs<<" //----------------- End of selection efficiency ---------------\\ "<<std::endl;
   ofs<<"NormFactor : "<<NormFactor<<std::endl;
   ofs<<" "<<std::endl;
   ofs.close();
std::cout<<"here6"<<std::endl;
    fillHisto("hData_Trks_Selection","Step",thesample,0,LooseTrks_Eff);
    fillHisto("hData_Trks_Selection","Step",thesample,1,TightTrks_Eff);

      //first : Event Yields for Run 2 with 138 fb-1..
      // ... Normalized by Number of events to get the correct errors

      //then : Efficiency of each step

    // then : Normalized by Number of events * cross section
    
    // then : Normalized by Number of events

  for (int i = 0 ; i < nentries ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,0,NormFactor*lumiRun2);// We fill for the Number of event counted by the selectino (numerator) to have the statistical 
                                                                                  // instead of filling with a weight eaquals to theNumber of event selected/ nentries.....
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,0,NormFactor/XS);
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,0,NormFactor/XS);
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,0,1*NormFactor);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,0,NormFactor/XS);
    }
  for (int i = 0 ; i < nFilterEvt ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,1,NormFactor*lumiRun2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,1,NormFactor/XS);
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,1,NormFactor/XS);
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,1,1*NormFactor);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,1,NormFactor/XS);
    }

  for (int i = 0 ; i < nEvts ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,2,NormFactor*lumiRun2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,2,1/nFilterEvt);
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,2,NormFactor/XS);
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,2,1*NormFactor);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,2,NormFactor/XS);
    }
    

  for (int i = 0 ; i < nRecoVertex ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,3,lumiRun2*NormFactor/2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,3,1/(2*nEvts));
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,3,NormFactor/(2*XS));
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,3,1*NormFactor/2);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,3,NormFactor/(2*XS));
    }
    

  for (int i = 0 ; i < nSelecBDTVtx ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,4,lumiRun2*NormFactor/2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,4,1/(2*nRecoVertex));
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,4,NormFactor/(2*XS));
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,4,1*NormFactor/2);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,4,NormFactor/(2*XS));
    }

  for (int i = 0 ; i < nRecoVertexTightWP ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,5,lumiRun2*NormFactor/2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,5,1/(2*nSelecBDTVtx));
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,5,NormFactor/(2*XS));
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,5,1*NormFactor/2);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,5,NormFactor/(2*XS));
    }
    

  for (int i = 0 ; i < nSelecBDTVtxTight_step12 ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,6,lumiRun2*NormFactor/2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,6,1/(2*nRecoVertexTightWP));
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,6,NormFactor/(2*XS));
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,6,1*NormFactor/2);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,6,NormFactor/(2*XS));
    }
    

  for (int i = 0 ; i < nEvts_w2TightVtx ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,7,NormFactor*lumiRun2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,7,2*1/nRecoVertexTightWP);
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,7,NormFactor/XS);
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,7,1*NormFactor);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,7,NormFactor/(XS));
    }
    

  for (int i = 0 ; i < nEvts_w2TightBDTVtx_step12 ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,8,NormFactor*lumiRun2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,8,2*1/nRecoVertexTightWP);
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,8,NormFactor/XS);
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,8,1*NormFactor);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,8,NormFactor/(XS));
    }
    

  for (int i = 0 ; i < nEvts_w1TightVtx ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,9,NormFactor*lumiRun2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,9,2*1/nRecoVertexTightWP);
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,9,NormFactor/XS);
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,9,1*NormFactor);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,9,NormFactor/(XS));
    }
    

  for (int i = 0 ; i < nEvts_w1TightBDTVtx_step12 ; i++)
    {
      fillHisto("hData_EventYields","Step",thesample,10,NormFactor*lumiRun2);
      fillHisto("hData_Vtx_SelectionEff","Step",thesample,10,2*1/nRecoVertexTightWP);
      fillHisto("hData_Vtx_SelectionEff","Total",thesample,10,NormFactor/XS);
      fillHisto("hData_Vtx_SelectionNXSNormalized","Step",thesample,10,1*NormFactor);
      fillHisto("hData_Vtx_SelectionNNormalized","Step",thesample,10,NormFactor/(XS));
    }
    



  
//Non-normalized Event yields
    fillHisto("hData_Vtx_Selection","Step",thesample,0,nentries);
    fillHisto("hData_Vtx_Selection","Step",thesample,1,nFilterEvt);
    fillHisto("hData_Vtx_Selection","Step",thesample,2,nEvts);
    fillHisto("hData_Vtx_Selection","Step",thesample,3,nRecoVertex/2);
    fillHisto("hData_Vtx_Selection","Step",thesample,4,nSelecBDTVtx/2);
    fillHisto("hData_Vtx_Selection","Step",thesample,5,nRecoVertexTightWP/2);
    fillHisto("hData_Vtx_Selection","Step",thesample,6,nSelecBDTVtxTight_step12/2);
    fillHisto("hData_Vtx_Selection","Step",thesample,7,nEvts_w2TightVtx);
    fillHisto("hData_Vtx_Selection","Step",thesample,8,nEvts_w2TightBDTVtx_step12);
    fillHisto("hData_Vtx_Selection","Step",thesample,9,nEvts_w1TightVtx);
    fillHisto("hData_Vtx_Selection","Step",thesample,10,nEvts_w1TightBDTVtx_step12);


//Normalization histogram
    fillHisto("hData_Vtx_Total","Step",thesample,0,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,1,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,2,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,3,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,4,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,5,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,6,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,7,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,8,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,9,1/NormFactor);
    fillHisto("hData_Vtx_Total","Step",thesample,10,1/NormFactor);




   theoutputfile->Write();
   //deleteHisto();
   theoutputfile->Close();
   delete theoutputfile;
}



void TreeReader::initializeHisto(TString sample, bool isfirstset){


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

   // ---------------------------- Filter (Online+Offline Selection)--------------
// Dimuon : HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
// OR
// Single Muon : HLT_IsoMu24_v
// Offline selection of muons : Two Loose ID and PFIsoLoose muons required
// |dxy | < 0.1 cm & |dz | < 0.2 cm (prompt)
// pt1 > 25GeV & pt2 > 10GeV
// Mμμ > 10 GeV (remove low-resonances)

//-----------Daniel ------//
    addHisto("hData_njet","NoSel",            sample.Data(),21,-0.5,20.5);
    addHisto("hData_njetNOmu","NoSel",        sample.Data(),21,-0.5,20.5);
    addHisto("hData_Mmumu","NoSel",           sample.Data(),25,0.,500.);
    addHisto("hData_BDTevt","NoSel",          sample.Data(),50,-1.,1.);
    addHisto("hData_Hemi_njet","NoSel",       sample.Data(),21,-0.5,20.5);
    addHisto("hData_Hemi_mass","NoSel",       sample.Data(),30,0.,600.);
    addHisto("hData_Hemi_massTot","NoSel",    sample.Data(),25,0.,1000.);
    addHisto("hData_Hemi_Vtx_njet","NoSel",   sample.Data(),21,-0.5,20.5);
    addHisto("hData_Hemi_Vtx_mass","NoSel",   sample.Data(),30,0.,600.);
    addHisto("hData_Hemi_Vtx_massTot","NoSel",sample.Data(),25,0.,1000.);
    addHisto("hData_Hemi1_njet","NoSel",      sample.Data(),21,-0.5,20.5);
    addHisto("hData_Hemi1_mass","NoSel",      sample.Data(),30,0.,600.);
    addHisto("hData_Hemi1_massTot","NoSel",   sample.Data(),25,0.,1000.);
    addHisto("hData_Hemi1_Vtx_dist","NoSel",  sample.Data(),25,0.,100.);
    addHisto("hData_Hemi1_Vtx_BDTvtx","NoSel",sample.Data(),50,-1.,1.);
    addHisto("hData_Hemi1_Vtx_njet","NoSel",  sample.Data(),21,-0.5,20.5);
    addHisto("hData_Hemi1_Vtx_mass","NoSel",  sample.Data(),30,0.,600.);
    addHisto("hData_Hemi1_Vtx_massTot","NoSel",sample.Data(),25,0.,1000.);
    addHisto("hData_Hemi2_njet","NoSel",      sample.Data(),21,-0.5,20.5);
    addHisto("hData_Hemi2_mass","NoSel",      sample.Data(),30,0.,600.);
    addHisto("hData_Hemi2_massTot","NoSel",   sample.Data(),25,0.,1000.);
    addHisto("hData_Hemi2_Vtx_dist","NoSel",  sample.Data(),25,0.,100.);
    addHisto("hData_Hemi2_Vtx_BDTvtx","NoSel",sample.Data(),50,-1.,1.);
    addHisto("hData_Hemi2_Vtx_njet","NoSel",  sample.Data(),21,-0.5,20.5);
    addHisto("hData_Hemi2_Vtx_mass","NoSel",  sample.Data(),30,0.,600.);
    addHisto("hData_Hemi2_Vtx_massTot","NoSel",sample.Data(),25,0.,1000.);
    
    addHisto("hData_vtx1_Mmumu","NoSel",      sample.Data(),25,0.,500.);
    addHisto("hData_BDTevt1_Mmumu","NoSel",  sample.Data(),25,0.,500.);
    addHisto("hData_BDTvtx1_Mmumu","NoSel",   sample.Data(),25,0.,500.);
    addHisto("hData_BDTevtvtx1_Mmumu","NoSel",sample.Data(),25,0.,500.);
    addHisto("hData_vtx2_Mmumu","NoSel",      sample.Data(),25,0.,500.);
    addHisto("hData_BDTevt2_Mmumu","NoSel",  sample.Data(),25,0.,500.);
    addHisto("hData_BDTvtx2_Mmumu","NoSel",   sample.Data(),25,0.,500.);
    addHisto("hData_BDTevtvtx2_Mmumu","NoSel",sample.Data(),25,0.,500.);

  //------------------------Gen Infos --------------------------------------------

  addHisto("hGen_Msmu","noSel",sample.Data(),361,159.5,520.5);
  addHisto("hGen_Mneu","noSel",sample.Data(),361,159.5,520.5);
  addHisto("hGen_ct0","noSel",sample.Data(),200,0.,200.);

  //-------------------------- EVTS ----------------------------------------------
   addHisto("Filter", "Offline+Online", sample.Data(), 2,0,2);
   addHisto("hData_Evts","BDT",  sample.Data(),101,-1,1);
   //----------------------------- Muons -----------------------------------------

   addHisto("RecoMuo_pT", "noSel",                      sample.Data(),  100, 0, 100);

   addHisto("DiMuon_Mass","noSel",                      sample.Data(),  24,0,600); //!!!!
   addHisto("DiMuon_Mass","EVTSel",                     sample.Data(),  24,0,600); //!!!!
   addHisto("DiMuon_Mass","RecoVtx",                    sample.Data(),  24,0,600); //!!!!
   addHisto("DiMuon_Mass","RecoVtx_TightWP",            sample.Data(),  24,0,600); //!!!!
   addHisto("DiMuon_Mass","RecoVtx_Good_BDT3",          sample.Data(),  24,0,600); //!!!!
   addHisto("DiMuon_Mass","RecoVtx_TightWP_BDT3",       sample.Data(),  24,0,600); //!!!!
   addHisto("DiMuon_Mass","RecoVtx_TightWP_BDT3bis",    sample.Data(),  24,0,600); //!!!!
   addHisto("DiMuon_Mass","2TightVtx",                  sample.Data(),  24,0,600); //!!!! 
   addHisto("DiMuon_Mass","1TightVtx",                  sample.Data(),  24,0,600); //!!!!
   
   addHisto("Evts_MVAVal","noSel", sample.Data() , 101,-1,1);


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

   addHisto2D("hData_reco_SecInt_xy","TrackerMatched",  sample.Data(), 1000,-25,25, 1000,-25,25);
   addHisto2D("hData_reco_SecInt_xy","Selec",           sample.Data(), 1000,-25,25, 1000,-25,25);
   addHisto2D("hData_reco_SecInt_rz","TrackerMatched",  sample.Data(), 140,0,70, 240,0,120);
   addHisto2D("hData_reco_SecInt_rz","Selec",           sample.Data(), 140,0,70, 240,0,120);

   addHisto("TranseverseResolution","NoSel",  sample.Data(),    100 , -2 , 2);
   addHisto("TranseverseDistance","NoSel",    sample.Data(),    100 , 0 , 100);
   addHisto2D("TransverseResVsDistance","NoSel",sample.Data(),  100 , 0 , 100,100 , -2 , 2);

   addHisto("LongitudinalResolution","NoSel", sample.Data(),    100 , -2 , 2);
   addHisto("LongitudinalDistance","NoSel",   sample.Data(),    401 , -200 , 200);
   addHisto2D("LongitudinalResVsDistance","NoSel",sample.Data(),   401 , -200 , 200,100 , -2 , 2);

   addHisto("3DResolution","NoSel",           sample.Data(),    100 , 0 , 3);
   addHisto("3DDistance","NoSel",             sample.Data(),    200 , 0 , 200);
   addHisto2D("3DResVsDistance","NoSel",      sample.Data(),  100 , 0 , 3,200 , 0 , 200);

   // ------------------------ Jets   -------------------------------------------

   addHisto("hData_jet_pt","",            sample.Data(),  800,0,1600);
   addHisto("hData_jet_eta","",           sample.Data(),  26,-6.5,6.5);
   addHisto("hData_jet_btag_Deepjet","",  sample.Data(),  100,0,1);
   addHisto("hData_jet_HadronFlavour","", sample.Data(),  6,0,6);

   addHisto2D("hData_BtagEff_Denom","",   sample.Data(),100,0,1000,8,-4,4 );
   addHisto2D("hData_BtagEff_Num","",     sample.Data(),100,0,1000,26,-6.5,6.5 );
   // ------------------------ Tracks -------------------------------------------

   addHisto("hData_MVAVal","noSel",                sample.Data(),101,-1,1);

   // ------------------------ Axes --------------------------------------------

  addHisto("hData_dR_RecoReco","noSel",            sample.Data(),201,0,6);
   addHisto("hSim_dR_GenGen","noSel",    sample.Data(), 201, 0,6);
   addHisto("hSim_dPhi_GenGen","noSel",   sample.Data(), 201, 0,6);
   addHisto("hSim_dEta_GenGen","noSel",   sample.Data(), 201, 0,6);
   addHisto("hSim_dR_GenReco","noSel",    sample.Data(), 201, 0,6);
   //---------------------------Signal-Vertices---------------------------------------

   addHisto("hData_Hemi_Vtx_NChi2","noSel",        sample.Data(),70,-20.,50.);
   addHisto("hData_Hemi_Vtx_nTrks","noSel",        sample.Data(),40,0.5,40.5);

   addHisto("hData_Hemi_Vtx_NChi2","RecoVtx",      sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_nTrks","RecoVtx",      sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_dist","RecoVtx",       sample.Data(),20,0,100);
   addHisto("hData_Hemi_Vtx_step","RecoVtx",       sample.Data(),5,0,5);
   addHisto("hData_Hemi_Vtx_eta","RecoVtx",        sample.Data(),26,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","RecoVtx",          sample.Data(),50,0,100);
   addHisto("hData_Hemi_Vtx_z","RecoVtx",          sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","RecoVtx",     sample.Data(),101,-1,1);

   addHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx",  sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx",  sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_dist","GoodRecoVtx",   sample.Data(),20,0,100);
   addHisto("hData_Hemi_Vtx_step","GoodRecoVtx",   sample.Data(),5,0,5);
   addHisto("hData_Hemi_Vtx_eta","GoodRecoVtx",    sample.Data(),26,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","GoodRecoVtx",      sample.Data(),50,0,100);
   addHisto("hData_Hemi_Vtx_z","GoodRecoVtx",      sample.Data() ,401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx", sample.Data(),101,-1,1);

   addHisto("hData_Hemi_Vtx_nTrks","RecoVtx_TightWP",    sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_NChi2","RecoVtx_TightWP",    sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_dist","RecoVtx_TightWP",     sample.Data(),20,0,100);
   addHisto("hData_Hemi_Vtx_eta","RecoVtx_TightWP",      sample.Data(),26,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","RecoVtx_TightWP",        sample.Data(),50,0,100);
   addHisto("hData_Hemi_Vtx_z","RecoVtx_TightWP",        sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","RecoVtx_TightWP",   sample.Data(),101,-1,1);

   addHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx_TightWP",   sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx_TightWP",   sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_dist","GoodRecoVtx_TightWP",    sample.Data(),20,0,100);
   addHisto("hData_Hemi_Vtx_eta","GoodRecoVtx_TightWP",    sample.Data(),26,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","GoodRecoVtx_TightWP",       sample.Data(),50,0,100);
   addHisto("hData_Hemi_Vtx_z","GoodRecoVtx_TightWP",       sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx_TightWP",  sample.Data(),101,-1,1);
   
   addHisto("hData_Hemi_Vtx_nTrks","RecoVtx_LooseWP",       sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_NChi2","RecoVtx_LooseWP",       sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_dist","RecoVtx_LooseWP",        sample.Data(),20,0,100);
   addHisto("hData_Hemi_Vtx_eta","RecoVtx_LooseWP",         sample.Data(),26,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","RecoVtx_LooseWP",           sample.Data(),50,0,100);
   addHisto("hData_Hemi_Vtx_z","RecoVtx_LooseWP",           sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","RecoVtx_LooseWP",      sample.Data(),101,-1,1);

   addHisto("hData_Hemi_Vtx_NChi2","GoodRecoVtx_LooseWP",   sample.Data(),50,-10.,40.);
   addHisto("hData_Hemi_Vtx_nTrks","GoodRecoVtx_LooseWP",   sample.Data(),40,0.5,40.5);
   addHisto("hData_Hemi_Vtx_dist","GoodRecoVtx_LooseWP",    sample.Data(),20,0,100);
   addHisto("hData_Hemi_Vtx_eta","GoodRecoVtx_LooseWP",     sample.Data(),26,-6.5,6.5);
   addHisto("hData_Hemi_Vtx_r","GoodRecoVtx_LooseWP",       sample.Data(),50,0,100);
   addHisto("hData_Hemi_Vtx_z","GoodRecoVtx_LooseWP",       sample.Data(),401,-200,200);
   addHisto("hData_Hemi_Vtx_MVAval","GoodRecoVtx_LooseWP",  sample.Data(),101,-1,1);  

  //Gen Lvel info

   addHisto("hSim_Hemi_Vtx_r","noSel",      sample.Data(),50,0,100); 
   addHisto("hSim_Hemi_Vtx_eta","noSel",    sample.Data(),26,-6.5,6.5); 
   addHisto("hSim_Hemi_Vtx_dist","noSel",   sample.Data(),20,0,100);  


   addHisto("hData_Hemi_Vtx_r","Ping",      sample.Data(),50,0,100); 
   addHisto("hData_Hemi_Vtx_eta","Ping",    sample.Data(),26,-6.5,6.5); 
   addHisto("hData_Hemi_Vtx_dist","Ping",   sample.Data(),20,0,100);
   addHisto("hSim_Hemi_Vtx_r","Ping",      sample.Data(),50,0,100); 
   addHisto("hSim_Hemi_Vtx_eta","Ping",    sample.Data(),26,-6.5,6.5); 
   addHisto("hSim_Hemi_Vtx_dist","Ping",   sample.Data(),20,0,100);                     
 
   addHisto("hData_Hemi_Vtx_r","Ping_TightWP",      sample.Data(),50,0,100); 
   addHisto("hData_Hemi_Vtx_eta","Ping_TightWP",    sample.Data(),26,-6.5,6.5); 
   addHisto("hData_Hemi_Vtx_dist","Ping_TightWP",   sample.Data(),20,0,100);
   addHisto("hSim_Hemi_Vtx_r","Ping_TightWP",      sample.Data(),50,0,100); 
   addHisto("hSim_Hemi_Vtx_eta","Ping_TightWP",    sample.Data(),26,-6.5,6.5); 
   addHisto("hSim_Hemi_Vtx_dist","Ping_TightWP",   sample.Data(),20,0,100);                  

   addHisto("hData_Hemi_Vtx_r","Ping_LooseWP",      sample.Data(),50,0,100); 
   addHisto("hData_Hemi_Vtx_eta","Ping_LooseWP",    sample.Data(),26,-6.5,6.5); 
   addHisto("hData_Hemi_Vtx_dist","Ping_LooseWP",   sample.Data(),20,0,100);
   addHisto("hSim_Hemi_Vtx_r","Ping_LooseWP",      sample.Data(),50,0,100); 
   addHisto("hSim_Hemi_Vtx_eta","Ping_LooseWP",    sample.Data(),26,-6.5,6.5); 
   addHisto("hSim_Hemi_Vtx_dist","Ping_LooseWP",   sample.Data(),20,0,100);


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
   addHisto("hData_Vtx_nTrk20","TightWP",                  sample.Data(),50,0,50);
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
   addHisto("hData_Hemi_SmuonMass","NoSel",                sample.Data(),200,0,2000);
   addHisto("hData_HighSmuonMass","NoSel",                 sample.Data(),200,0,2000);

   //----------Tracks Selection -------------//

   addHisto("hData_Trks_Selection","Step",                 sample.Data(),3,0,3);

    //----------Vtx Selection -------------//
    addHisto("hData_EventYields","Step"   ,                sample.Data(),11,0,11);
    addHisto("hData_Vtx_Selection","Step",                 sample.Data(),11,0,11);
    addHisto("hData_Vtx_SelectionNXSNormalized","Step",    sample.Data(),11,0,11);

    addHisto("hData_Vtx_SelectionNNormalized","Step",       sample.Data(),11,0,11);
    addHisto("hData_Vtx_Total","Step",                      sample.Data(),11,0,11);
    addHisto("hData_Vtx_SelectionEff","Step",               sample.Data(),11,0,11);
    addHisto("hData_Vtx_SelectionEff","Total",              sample.Data(),11,0,11);

    //-------------Gen Infos---------------//

    
    addHisto("hGen_ct0","smu200_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu250_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu300_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu300_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu300_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu300_neu280",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu400_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu400_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu400_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu400_neu300",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu400_neu350",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu400_neu380",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu500_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu500_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu500_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu500_neu300",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu500_neu350",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu500_neu400",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu500_neu450",  sample.Data(), 200,0.,200.);
    addHisto("hGen_ct0","smu500_neu480",  sample.Data(), 200,0.,200.);

    addHisto("hSim_Filter","smu200_neu180",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu250_neu200",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu300_neu180",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu300_neu200",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu300_neu250",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu300_neu280",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu400_neu180",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu400_neu200",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu400_neu250",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu400_neu300",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu400_neu350",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu400_neu380",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu500_neu180",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu500_neu200",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu500_neu250",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu500_neu300",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu500_neu350",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu500_neu400",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu500_neu450",  sample.Data(),2,-0.5,1.5);
    addHisto("hSim_Filter","smu500_neu480",  sample.Data(),2,-0.5,1.5);

    addHisto("hSim_Hemi_dR","smu200_neu180",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu250_neu200",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu300_neu180",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu300_neu200",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu300_neu250",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu300_neu280",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu400_neu180",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu400_neu200",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu400_neu250",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu400_neu300",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu400_neu350",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu400_neu380",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu500_neu180",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu500_neu200",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu500_neu250",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu500_neu300",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu500_neu350",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu500_neu400",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu500_neu450",  sample.Data(),33,0.,3.3);
    addHisto("hSim_Hemi_dR","smu500_neu480",  sample.Data(),33,0.,3.3);

    addHisto("hData_Hemi_Vtx_dist","smu200_neu180",  sample.Data(),200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu250_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu300_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu300_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu300_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu300_neu280",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu400_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu400_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu400_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu400_neu300",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu400_neu350",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu400_neu380",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu500_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu500_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu500_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu500_neu300",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu500_neu350",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu500_neu400",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu500_neu450",  sample.Data(), 200,0.,200.);
    addHisto("hData_Hemi_Vtx_dist","smu500_neu480",  sample.Data(), 200,0.,200.);

    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu200_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu250_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu300_neu280",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu300",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu350",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu400_neu380",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu300",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu350",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu400",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu450",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_all","smu500_neu480",  sample.Data(), 200,0.,200.);

    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu200_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu250_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu300_neu280",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu300",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu350",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu400_neu380",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu180",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu200",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu250",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu300",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu350",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu400",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu450",  sample.Data(), 200,0.,200.);
    addHisto("hSim_Hemi_Vtx_dist_ping_step1","smu500_neu480",  sample.Data(), 200,0.,200.);


    addHisto("hSim_EVTBDT","dm20",   sample.Data(),101,-1,1);
    addHisto("hSim_EVTBDT","dm50",   sample.Data(),101,-1,1);
    addHisto("hSim_EVTBDT","dm100",   sample.Data(),101,-1,1);
    addHisto("hSim_EVTBDT","dm150",   sample.Data(),101,-1,1);
    addHisto("hSim_EVTBDT","dm200",   sample.Data(),101,-1,1);
    addHisto("hSim_EVTBDT","dm250",   sample.Data(),101,-1,1);
    addHisto("hSim_EVTBDT","dm300",   sample.Data(),101,-1,1);
    addHisto("hSim_EVTBDT","dm320",   sample.Data(),101,-1,1);

    addHisto("hSim_ct0","dm20",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm50",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm100",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm150",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm200",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm250",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm300",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm320",  sample.Data(),200,0.,200.);

    addHisto("hSim_TRKBDT","dm20",  sample.Data(),101,-1,1);
    addHisto("hSim_TRKBDT","dm50",  sample.Data(),101,-1,1);
    addHisto("hSim_TRKBDT","dm100",  sample.Data(),101,-1,1);
    addHisto("hSim_TRKBDT","dm150",  sample.Data(),101,-1,1);
    addHisto("hSim_TRKBDT","dm200",  sample.Data(),101,-1,1);
    addHisto("hSim_TRKBDT","dm250",  sample.Data(),101,-1,1);
    addHisto("hSim_TRKBDT","dm300",  sample.Data(),101,-1,1);
    addHisto("hSim_TRKBDT","dm320",  sample.Data(),101,-1,1);

    addHisto("hSim_ct0","dm20_TRKBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm50_TRKBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm100_TRKBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm150_TRKBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm200_TRKBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm250_TRKBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm300_TRKBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm320_TRKBDT",  sample.Data(),200,0.,200.);

    addHisto("hSim_VTXBDT","dm20",  sample.Data(),101,-1,1);
    addHisto("hSim_VTXBDT","dm50",  sample.Data(),101,-1,1);
    addHisto("hSim_VTXBDT","dm100",  sample.Data(),101,-1,1);
    addHisto("hSim_VTXBDT","dm150",  sample.Data(),101,-1,1);
    addHisto("hSim_VTXBDT","dm200",  sample.Data(),101,-1,1);
    addHisto("hSim_VTXBDT","dm250",  sample.Data(),101,-1,1);
    addHisto("hSim_VTXBDT","dm300",  sample.Data(),101,-1,1);
    addHisto("hSim_VTXBDT","dm320",  sample.Data(),101,-1,1);

    addHisto("hSim_ct0","dm20_VTXBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm50_VTXBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm100_VTXBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm150_VTXBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm200_VTXBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm250_VTXBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm300_VTXBDT",  sample.Data(),200,0.,200.);
    addHisto("hSim_ct0","dm320_VTXBDT",  sample.Data(),200,0.,200.);

    addHisto("hBDT_MET","EVT",                  sample.Data(),150,0,1500);
    addHisto("hBDT_TRACK_SIZE","EVT",           sample.Data(),4100,0,4100);
    addHisto("hBDT_muon_leadingpt","EVT",       sample.Data(),100,20,520);
    addHisto("hBDT_muon_leadingpt2","EVT",      sample.Data(),100,10,510);
    addHisto("hBDT_jet_leadingpt","EVT",        sample.Data(),100,30,1030);
    addHisto("hBDT_jet_leadingpt2","EVT",       sample.Data(),100,20,1020);
    addHisto("hBDT_muon_muon_dR","EVT",         sample.Data(),100,0,5);
    addHisto("hBDT_muon_muon_dPhi","EVT",       sample.Data(),70,0,3.5);
    addHisto("hBDT_muon_muon_dEta","EVT",       sample.Data(),80,0,4);
    addHisto("hBDT_jet_jet_dR","EVT",           sample.Data(),100,0,5);
    addHisto("hBDT_jet_jet_dPhi","EVT",         sample.Data(),70,0,3.5);
    addHisto("hBDT_jet_jet_dEta","EVT",         sample.Data(),80,0,4);
    addHisto("hBDT_muon_jet_dRmin","EVT",       sample.Data(),50,0,5);
    addHisto("hBDT_muon_jet_dRmax","EVT",       sample.Data(),50,0,5);
    addHisto("hBDT_HT","EVT",                   sample.Data(),200,100,2100);
    addHisto("hBDT_ST","EVT",                   sample.Data(),50,30,530);
    addHisto("hBDT_njet","EVT",                 sample.Data(),20,0,20);
    addHisto("hBDT_muon_nmu","EVT",             sample.Data(),20,0,20);
    addHisto("hBDT_Hemi_LooseBTag_axes","EVT",  sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_MediumBTag_axes","EVT", sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_TightBTag_axes","EVT",  sample.Data(),2,0,2);

    // eat sleep and repeat ---------------------//

    addHisto("hBDT_MET","EVT_dm20",                  sample.Data(),150,0,1500);
    addHisto("hBDT_TRACK_SIZE","EVT_dm20",           sample.Data(),4100,0,4100);
    addHisto("hBDT_muon_leadingpt","EVT_dm20",       sample.Data(),100,20,520);
    addHisto("hBDT_muon_leadingpt2","EVT_dm20",      sample.Data(),100,10,510);
    addHisto("hBDT_jet_leadingpt","EVT_dm20",        sample.Data(),100,30,1030);
    addHisto("hBDT_jet_leadingpt2","EVT_dm20",       sample.Data(),100,20,1020);
    addHisto("hBDT_muon_muon_dR","EVT_dm20",         sample.Data(),100,0,5);
    addHisto("hBDT_muon_muon_dPhi","EVT_dm20",       sample.Data(),70,0,3.5);
    addHisto("hBDT_muon_muon_dEta","EVT_dm20",       sample.Data(),80,0,4);
    addHisto("hBDT_jet_jet_dR","EVT_dm20",           sample.Data(),100,0,5);
    addHisto("hBDT_jet_jet_dPhi","EVT_dm20",         sample.Data(),70,0,3.5);
    addHisto("hBDT_jet_jet_dEta","EVT_dm20",         sample.Data(),80,0,4);
    addHisto("hBDT_muon_jet_dRmin","EVT_dm20",       sample.Data(),50,0,5);
    addHisto("hBDT_muon_jet_dRmax","EVT_dm20",       sample.Data(),50,0,5);
    addHisto("hBDT_HT","EVT_dm20",                   sample.Data(),200,100,2100);
    addHisto("hBDT_ST","EVT_dm20",                   sample.Data(),50,30,530);
    addHisto("hBDT_njet","EVT_dm20",                 sample.Data(),20,0,20);
    addHisto("hBDT_muon_nmu","EVT_dm20",             sample.Data(),20,0,20);
    addHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm20",  sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm20", sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_TightBTag_axes","EVT_dm20",  sample.Data(),2,0,2);

        addHisto("hBDT_MET","EVT_dm50",                  sample.Data(),150,0,1500);
    addHisto("hBDT_TRACK_SIZE","EVT_dm50",           sample.Data(),4100,0,4100);
    addHisto("hBDT_muon_leadingpt","EVT_dm50",       sample.Data(),100,20,520);
    addHisto("hBDT_muon_leadingpt2","EVT_dm50",      sample.Data(),100,10,510);
    addHisto("hBDT_jet_leadingpt","EVT_dm50",        sample.Data(),100,30,1030);
    addHisto("hBDT_jet_leadingpt2","EVT_dm50",       sample.Data(),100,20,1020);
    addHisto("hBDT_muon_muon_dR","EVT_dm50",         sample.Data(),100,0,5);
    addHisto("hBDT_muon_muon_dPhi","EVT_dm50",       sample.Data(),70,0,3.5);
    addHisto("hBDT_muon_muon_dEta","EVT_dm50",       sample.Data(),80,0,4);
    addHisto("hBDT_jet_jet_dR","EVT_dm50",           sample.Data(),100,0,5);
    addHisto("hBDT_jet_jet_dPhi","EVT_dm50",         sample.Data(),70,0,3.5);
    addHisto("hBDT_jet_jet_dEta","EVT_dm50",         sample.Data(),80,0,4);
    addHisto("hBDT_muon_jet_dRmin","EVT_dm50",       sample.Data(),50,0,5);
    addHisto("hBDT_muon_jet_dRmax","EVT_dm50",       sample.Data(),50,0,5);
    addHisto("hBDT_HT","EVT_dm50",                   sample.Data(),200,100,2100);
    addHisto("hBDT_ST","EVT_dm50",                   sample.Data(),50,30,530);
    addHisto("hBDT_njet","EVT_dm50",                 sample.Data(),20,0,20);
    addHisto("hBDT_muon_nmu","EVT_dm50",             sample.Data(),20,0,20);
    addHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm50",  sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm50", sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_TightBTag_axes","EVT_dm50",  sample.Data(),2,0,2);

        addHisto("hBDT_MET","EVT_dm100",                  sample.Data(),150,0,1500);
    addHisto("hBDT_TRACK_SIZE","EVT_dm100",           sample.Data(),4100,0,4100);
    addHisto("hBDT_muon_leadingpt","EVT_dm100",       sample.Data(),100,20,520);
    addHisto("hBDT_muon_leadingpt2","EVT_dm100",      sample.Data(),100,10,510);
    addHisto("hBDT_jet_leadingpt","EVT_dm100",        sample.Data(),100,30,1030);
    addHisto("hBDT_jet_leadingpt2","EVT_dm100",       sample.Data(),100,20,1020);
    addHisto("hBDT_muon_muon_dR","EVT_dm100",         sample.Data(),100,0,5);
    addHisto("hBDT_muon_muon_dPhi","EVT_dm100",       sample.Data(),70,0,3.5);
    addHisto("hBDT_muon_muon_dEta","EVT_dm100",       sample.Data(),80,0,4);
    addHisto("hBDT_jet_jet_dR","EVT_dm100",           sample.Data(),100,0,5);
    addHisto("hBDT_jet_jet_dPhi","EVT_dm100",         sample.Data(),70,0,3.5);
    addHisto("hBDT_jet_jet_dEta","EVT_dm100",         sample.Data(),80,0,4);
    addHisto("hBDT_muon_jet_dRmin","EVT_dm100",       sample.Data(),50,0,5);
    addHisto("hBDT_muon_jet_dRmax","EVT_dm100",       sample.Data(),50,0,5);
    addHisto("hBDT_HT","EVT_dm100",                   sample.Data(),200,100,2100);
    addHisto("hBDT_ST","EVT_dm100",                   sample.Data(),50,30,530);
    addHisto("hBDT_njet","EVT_dm100",                 sample.Data(),20,0,20);
    addHisto("hBDT_muon_nmu","EVT_dm100",             sample.Data(),20,0,20);
    addHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm100",  sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm100", sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_TightBTag_axes","EVT_dm100",  sample.Data(),2,0,2);

        addHisto("hBDT_MET","EVT_dm150",                  sample.Data(),150,0,1500);
    addHisto("hBDT_TRACK_SIZE","EVT_dm150",           sample.Data(),4100,0,4100);
    addHisto("hBDT_muon_leadingpt","EVT_dm150",       sample.Data(),100,20,520);
    addHisto("hBDT_muon_leadingpt2","EVT_dm150",      sample.Data(),100,10,510);
    addHisto("hBDT_jet_leadingpt","EVT_dm150",        sample.Data(),100,30,1030);
    addHisto("hBDT_jet_leadingpt2","EVT_dm150",       sample.Data(),100,20,1020);
    addHisto("hBDT_muon_muon_dR","EVT_dm150",         sample.Data(),100,0,5);
    addHisto("hBDT_muon_muon_dPhi","EVT_dm150",       sample.Data(),70,0,3.5);
    addHisto("hBDT_muon_muon_dEta","EVT_dm150",       sample.Data(),80,0,4);
    addHisto("hBDT_jet_jet_dR","EVT_dm150",           sample.Data(),100,0,5);
    addHisto("hBDT_jet_jet_dPhi","EVT_dm150",         sample.Data(),70,0,3.5);
    addHisto("hBDT_jet_jet_dEta","EVT_dm150",         sample.Data(),80,0,4);
    addHisto("hBDT_muon_jet_dRmin","EVT_dm150",       sample.Data(),50,0,5);
    addHisto("hBDT_muon_jet_dRmax","EVT_dm150",       sample.Data(),50,0,5);
    addHisto("hBDT_HT","EVT_dm150",                   sample.Data(),200,100,2100);
    addHisto("hBDT_ST","EVT_dm150",                   sample.Data(),50,30,530);
    addHisto("hBDT_njet","EVT_dm150",                 sample.Data(),20,0,20);
    addHisto("hBDT_muon_nmu","EVT_dm150",             sample.Data(),20,0,20);
    addHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm150",  sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm150", sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_TightBTag_axes","EVT_dm150",  sample.Data(),2,0,2);

        addHisto("hBDT_MET","EVT_dm200",                  sample.Data(),150,0,1500);
    addHisto("hBDT_TRACK_SIZE","EVT_dm200",           sample.Data(),4100,0,4100);
    addHisto("hBDT_muon_leadingpt","EVT_dm200",       sample.Data(),100,20,520);
    addHisto("hBDT_muon_leadingpt2","EVT_dm200",      sample.Data(),100,10,510);
    addHisto("hBDT_jet_leadingpt","EVT_dm200",        sample.Data(),100,30,1030);
    addHisto("hBDT_jet_leadingpt2","EVT_dm200",       sample.Data(),100,20,1020);
    addHisto("hBDT_muon_muon_dR","EVT_dm200",         sample.Data(),100,0,5);
    addHisto("hBDT_muon_muon_dPhi","EVT_dm200",       sample.Data(),70,0,3.5);
    addHisto("hBDT_muon_muon_dEta","EVT_dm200",       sample.Data(),80,0,4);
    addHisto("hBDT_jet_jet_dR","EVT_dm200",           sample.Data(),100,0,5);
    addHisto("hBDT_jet_jet_dPhi","EVT_dm200",         sample.Data(),70,0,3.5);
    addHisto("hBDT_jet_jet_dEta","EVT_dm200",         sample.Data(),80,0,4);
    addHisto("hBDT_muon_jet_dRmin","EVT_dm200",       sample.Data(),50,0,5);
    addHisto("hBDT_muon_jet_dRmax","EVT_dm200",       sample.Data(),50,0,5);
    addHisto("hBDT_HT","EVT_dm200",                   sample.Data(),200,100,2100);
    addHisto("hBDT_ST","EVT_dm200",                   sample.Data(),50,30,530);
    addHisto("hBDT_njet","EVT_dm200",                 sample.Data(),20,0,20);
    addHisto("hBDT_muon_nmu","EVT_dm200",             sample.Data(),20,0,20);
    addHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm200",  sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm200", sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_TightBTag_axes","EVT_dm200",  sample.Data(),2,0,2);

        addHisto("hBDT_MET","EVT_dm250",                  sample.Data(),150,0,1500);
    addHisto("hBDT_TRACK_SIZE","EVT_dm250",           sample.Data(),4100,0,4100);
    addHisto("hBDT_muon_leadingpt","EVT_dm250",       sample.Data(),100,20,520);
    addHisto("hBDT_muon_leadingpt2","EVT_dm250",      sample.Data(),100,10,510);
    addHisto("hBDT_jet_leadingpt","EVT_dm250",        sample.Data(),100,30,1030);
    addHisto("hBDT_jet_leadingpt2","EVT_dm250",       sample.Data(),100,20,1020);
    addHisto("hBDT_muon_muon_dR","EVT_dm250",         sample.Data(),100,0,5);
    addHisto("hBDT_muon_muon_dPhi","EVT_dm250",       sample.Data(),70,0,3.5);
    addHisto("hBDT_muon_muon_dEta","EVT_dm250",       sample.Data(),80,0,4);
    addHisto("hBDT_jet_jet_dR","EVT_dm250",           sample.Data(),100,0,5);
    addHisto("hBDT_jet_jet_dPhi","EVT_dm250",         sample.Data(),70,0,3.5);
    addHisto("hBDT_jet_jet_dEta","EVT_dm250",         sample.Data(),80,0,4);
    addHisto("hBDT_muon_jet_dRmin","EVT_dm250",       sample.Data(),50,0,5);
    addHisto("hBDT_muon_jet_dRmax","EVT_dm250",       sample.Data(),50,0,5);
    addHisto("hBDT_HT","EVT_dm250",                   sample.Data(),200,100,2100);
    addHisto("hBDT_ST","EVT_dm250",                   sample.Data(),50,30,530);
    addHisto("hBDT_njet","EVT_dm250",                 sample.Data(),20,0,20);
    addHisto("hBDT_muon_nmu","EVT_dm250",             sample.Data(),20,0,20);
    addHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm250",  sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm250", sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_TightBTag_axes","EVT_dm250",  sample.Data(),2,0,2);

        addHisto("hBDT_MET","EVT_dm300",                  sample.Data(),150,0,1500);
    addHisto("hBDT_TRACK_SIZE","EVT_dm300",           sample.Data(),4100,0,4100);
    addHisto("hBDT_muon_leadingpt","EVT_dm300",       sample.Data(),100,20,520);
    addHisto("hBDT_muon_leadingpt2","EVT_dm300",      sample.Data(),100,10,510);
    addHisto("hBDT_jet_leadingpt","EVT_dm300",        sample.Data(),100,30,1030);
    addHisto("hBDT_jet_leadingpt2","EVT_dm300",       sample.Data(),100,20,1020);
    addHisto("hBDT_muon_muon_dR","EVT_dm300",         sample.Data(),100,0,5);
    addHisto("hBDT_muon_muon_dPhi","EVT_dm300",       sample.Data(),70,0,3.5);
    addHisto("hBDT_muon_muon_dEta","EVT_dm300",       sample.Data(),80,0,4);
    addHisto("hBDT_jet_jet_dR","EVT_dm300",           sample.Data(),100,0,5);
    addHisto("hBDT_jet_jet_dPhi","EVT_dm300",         sample.Data(),70,0,3.5);
    addHisto("hBDT_jet_jet_dEta","EVT_dm300",         sample.Data(),80,0,4);
    addHisto("hBDT_muon_jet_dRmin","EVT_dm300",       sample.Data(),50,0,5);
    addHisto("hBDT_muon_jet_dRmax","EVT_dm300",       sample.Data(),50,0,5);
    addHisto("hBDT_HT","EVT_dm300",                   sample.Data(),200,100,2100);
    addHisto("hBDT_ST","EVT_dm300",                   sample.Data(),50,30,530);
    addHisto("hBDT_njet","EVT_dm300",                 sample.Data(),20,0,20);
    addHisto("hBDT_muon_nmu","EVT_dm300",             sample.Data(),20,0,20);
    addHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm300",  sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm300", sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_TightBTag_axes","EVT_dm300",  sample.Data(),2,0,2);

        addHisto("hBDT_MET","EVT_dm320",                  sample.Data(),150,0,1500);
    addHisto("hBDT_TRACK_SIZE","EVT_dm320",           sample.Data(),4100,0,4100);
    addHisto("hBDT_muon_leadingpt","EVT_dm320",       sample.Data(),100,20,520);
    addHisto("hBDT_muon_leadingpt2","EVT_dm320",      sample.Data(),100,10,510);
    addHisto("hBDT_jet_leadingpt","EVT_dm320",        sample.Data(),100,30,1030);
    addHisto("hBDT_jet_leadingpt2","EVT_dm320",       sample.Data(),100,20,1020);
    addHisto("hBDT_muon_muon_dR","EVT_dm320",         sample.Data(),100,0,5);
    addHisto("hBDT_muon_muon_dPhi","EVT_dm320",       sample.Data(),70,0,3.5);
    addHisto("hBDT_muon_muon_dEta","EVT_dm320",       sample.Data(),80,0,4);
    addHisto("hBDT_jet_jet_dR","EVT_dm320",           sample.Data(),100,0,5);
    addHisto("hBDT_jet_jet_dPhi","EVT_dm320",         sample.Data(),70,0,3.5);
    addHisto("hBDT_jet_jet_dEta","EVT_dm320",         sample.Data(),80,0,4);
    addHisto("hBDT_muon_jet_dRmin","EVT_dm320",       sample.Data(),50,0,5);
    addHisto("hBDT_muon_jet_dRmax","EVT_dm320",       sample.Data(),50,0,5);
    addHisto("hBDT_HT","EVT_dm320",                   sample.Data(),200,100,2100);
    addHisto("hBDT_ST","EVT_dm320",                   sample.Data(),50,30,530);
    addHisto("hBDT_njet","EVT_dm320",                 sample.Data(),20,0,20);
    addHisto("hBDT_muon_nmu","EVT_dm320",             sample.Data(),20,0,20);
    addHisto("hBDT_Hemi_LooseBTag_axes","EVT_dm320",  sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_MediumBTag_axes","EVT_dm320", sample.Data(),2,0,2);
    addHisto("hBDT_Hemi_TightBTag_axes","EVT_dm320",  sample.Data(),2,0,2);





    // ------------ zzzzzzzzzzzzzzzzzzzzzzzz --///////////////
    addHisto("hBDT_lost","TRK",                 sample.Data(),2,0,2);
    addHisto("hBDT_dz","TRK",                   sample.Data(),400,-200,200);
    addHisto("hBDT_dxy","TRK",                  sample.Data(),200,-100,100);
    addHisto("hBDT_pt","TRK",                   sample.Data(),100,0,200);
    addHisto("hBDT_eta","TRK",                  sample.Data(),80,-4,4);
    addHisto("hBDT_NChi2","TRK",                sample.Data(),5,0,5);
    addHisto("hBDT_nhits","TRK",                sample.Data(),35,0,35);
    addHisto("hBDT_iJet","TRK",                 sample.Data(),20,0,20);
    addHisto("hBDT_drSig","TRK",                sample.Data(),1000,0,2000);
    addHisto("hBDT_dzSig","TRK",                sample.Data(),1000,0,2000);
    addHisto("hBDT_ntrk10","TRK",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk20","TRK",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk30","TRK",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk40","TRK",               sample.Data(),100,0,100);
    addHisto("hBDT_Hemi_dR","TRK",              sample.Data(),100,0,6);
    addHisto("hBDT_Hemi_dRmax","TRK",           sample.Data(),100,0,6);

    // eat sleep and repeat -------------------------//

    addHisto("hBDT_lost","TRK_dm20",                 sample.Data(),2,0,2);
    addHisto("hBDT_dz","TRK_dm20",                   sample.Data(),400,-200,200);
    addHisto("hBDT_dxy","TRK_dm20",                  sample.Data(),200,-100,100);
    addHisto("hBDT_pt","TRK_dm20",                   sample.Data(),100,0,200);
    addHisto("hBDT_eta","TRK_dm20",                  sample.Data(),80,-4,4);
    addHisto("hBDT_NChi2","TRK_dm20",                sample.Data(),5,0,5);
    addHisto("hBDT_nhits","TRK_dm20",                sample.Data(),35,0,35);
    addHisto("hBDT_iJet","TRK_dm20",                 sample.Data(),20,0,20);
    addHisto("hBDT_drSig","TRK_dm20",                sample.Data(),1000,0,2000);
    addHisto("hBDT_dzSig","TRK_dm20",                sample.Data(),1000,0,2000);
    addHisto("hBDT_ntrk10","TRK_dm20",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk20","TRK_dm20",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk30","TRK_dm20",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk40","TRK_dm20",               sample.Data(),100,0,100);
    addHisto("hBDT_Hemi_dR","TRK_dm20",              sample.Data(),100,0,6);
    addHisto("hBDT_Hemi_dRmax","TRK_dm20",           sample.Data(),100,0,6);

        addHisto("hBDT_lost","TRK_dm50",                 sample.Data(),2,0,2);
    addHisto("hBDT_dz","TRK_dm50",                   sample.Data(),400,-200,200);
    addHisto("hBDT_dxy","TRK_dm50",                  sample.Data(),200,-100,100);
    addHisto("hBDT_pt","TRK_dm50",                   sample.Data(),100,0,200);
    addHisto("hBDT_eta","TRK_dm50",                  sample.Data(),80,-4,4);
    addHisto("hBDT_NChi2","TRK_dm50",                sample.Data(),5,0,5);
    addHisto("hBDT_nhits","TRK_dm50",                sample.Data(),35,0,35);
    addHisto("hBDT_iJet","TRK_dm50",                 sample.Data(),20,0,20);
    addHisto("hBDT_drSig","TRK_dm50",                sample.Data(),1000,0,2000);
    addHisto("hBDT_dzSig","TRK_dm50",                sample.Data(),1000,0,2000);
    addHisto("hBDT_ntrk10","TRK_dm50",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk20","TRK_dm50",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk30","TRK_dm50",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk40","TRK_dm50",               sample.Data(),100,0,100);
    addHisto("hBDT_Hemi_dR","TRK_dm50",              sample.Data(),100,0,6);
    addHisto("hBDT_Hemi_dRmax","TRK_dm50",           sample.Data(),100,0,6);


        addHisto("hBDT_lost","TRK_dm100",                 sample.Data(),2,0,2);
    addHisto("hBDT_dz","TRK_dm100",                   sample.Data(),400,-200,200);
    addHisto("hBDT_dxy","TRK_dm100",                  sample.Data(),200,-100,100);
    addHisto("hBDT_pt","TRK_dm100",                   sample.Data(),100,0,200);
    addHisto("hBDT_eta","TRK_dm100",                  sample.Data(),80,-4,4);
    addHisto("hBDT_NChi2","TRK_dm100",                sample.Data(),5,0,5);
    addHisto("hBDT_nhits","TRK_dm100",                sample.Data(),35,0,35);
    addHisto("hBDT_iJet","TRK_dm100",                 sample.Data(),20,0,20);
    addHisto("hBDT_drSig","TRK_dm100",                sample.Data(),1000,0,2000);
    addHisto("hBDT_dzSig","TRK_dm100",                sample.Data(),1000,0,2000);
    addHisto("hBDT_ntrk10","TRK_dm100",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk20","TRK_dm100",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk30","TRK_dm100",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk40","TRK_dm100",               sample.Data(),100,0,100);
    addHisto("hBDT_Hemi_dR","TRK_dm100",              sample.Data(),100,0,6);
    addHisto("hBDT_Hemi_dRmax","TRK_dm100",           sample.Data(),100,0,6);

        addHisto("hBDT_lost","TRK_dm150",                 sample.Data(),2,0,2);
    addHisto("hBDT_dz","TRK_dm150",                   sample.Data(),400,-200,200);
    addHisto("hBDT_dxy","TRK_dm150",                  sample.Data(),200,-100,100);
    addHisto("hBDT_pt","TRK_dm150",                   sample.Data(),100,0,200);
    addHisto("hBDT_eta","TRK_dm150",                  sample.Data(),80,-4,4);
    addHisto("hBDT_NChi2","TRK_dm150",                sample.Data(),5,0,5);
    addHisto("hBDT_nhits","TRK_dm150",                sample.Data(),35,0,35);
    addHisto("hBDT_iJet","TRK_dm150",                 sample.Data(),20,0,20);
    addHisto("hBDT_drSig","TRK_dm150",                sample.Data(),1000,0,2000);
    addHisto("hBDT_dzSig","TRK_dm150",                sample.Data(),1000,0,2000);
    addHisto("hBDT_ntrk10","TRK_dm150",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk20","TRK_dm150",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk30","TRK_dm150",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk40","TRK_dm150",               sample.Data(),100,0,100);
    addHisto("hBDT_Hemi_dR","TRK_dm150",              sample.Data(),100,0,6);
    addHisto("hBDT_Hemi_dRmax","TRK_dm150",           sample.Data(),100,0,6);

        addHisto("hBDT_lost","TRK_dm200",                 sample.Data(),2,0,2);
    addHisto("hBDT_dz","TRK_dm200",                   sample.Data(),400,-200,200);
    addHisto("hBDT_dxy","TRK_dm200",                  sample.Data(),200,-100,100);
    addHisto("hBDT_pt","TRK_dm200",                   sample.Data(),100,0,200);
    addHisto("hBDT_eta","TRK_dm200",                  sample.Data(),80,-4,4);
    addHisto("hBDT_NChi2","TRK_dm200",                sample.Data(),5,0,5);
    addHisto("hBDT_nhits","TRK_dm200",                sample.Data(),35,0,35);
    addHisto("hBDT_iJet","TRK_dm200",                 sample.Data(),20,0,20);
    addHisto("hBDT_drSig","TRK_dm200",                sample.Data(),1000,0,2000);
    addHisto("hBDT_dzSig","TRK_dm200",                sample.Data(),1000,0,2000);
    addHisto("hBDT_ntrk10","TRK_dm200",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk20","TRK_dm200",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk30","TRK_dm200",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk40","TRK_dm200",               sample.Data(),100,0,100);
    addHisto("hBDT_Hemi_dR","TRK_dm200",              sample.Data(),100,0,6);
    addHisto("hBDT_Hemi_dRmax","TRK_dm200",           sample.Data(),100,0,6);

        addHisto("hBDT_lost","TRK_dm250",                 sample.Data(),2,0,2);
    addHisto("hBDT_dz","TRK_dm250",                   sample.Data(),400,-200,200);
    addHisto("hBDT_dxy","TRK_dm250",                  sample.Data(),200,-100,100);
    addHisto("hBDT_pt","TRK_dm250",                   sample.Data(),100,0,200);
    addHisto("hBDT_eta","TRK_dm250",                  sample.Data(),80,-4,4);
    addHisto("hBDT_NChi2","TRK_dm250",                sample.Data(),5,0,5);
    addHisto("hBDT_nhits","TRK_dm250",                sample.Data(),35,0,35);
    addHisto("hBDT_iJet","TRK_dm250",                 sample.Data(),20,0,20);
    addHisto("hBDT_drSig","TRK_dm250",                sample.Data(),1000,0,2000);
    addHisto("hBDT_dzSig","TRK_dm250",                sample.Data(),1000,0,2000);
    addHisto("hBDT_ntrk10","TRK_dm250",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk20","TRK_dm250",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk30","TRK_dm250",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk40","TRK_dm250",               sample.Data(),100,0,100);
    addHisto("hBDT_Hemi_dR","TRK_dm250",              sample.Data(),100,0,6);
    addHisto("hBDT_Hemi_dRmax","TRK_dm250",           sample.Data(),100,0,6);

        addHisto("hBDT_lost","TRK_dm300",                 sample.Data(),2,0,2);
    addHisto("hBDT_dz","TRK_dm300",                   sample.Data(),400,-200,200);
    addHisto("hBDT_dxy","TRK_dm300",                  sample.Data(),200,-100,100);
    addHisto("hBDT_pt","TRK_dm300",                   sample.Data(),100,0,200);
    addHisto("hBDT_eta","TRK_dm300",                  sample.Data(),80,-4,4);
    addHisto("hBDT_NChi2","TRK_dm300",                sample.Data(),5,0,5);
    addHisto("hBDT_nhits","TRK_dm300",                sample.Data(),35,0,35);
    addHisto("hBDT_iJet","TRK_dm300",                 sample.Data(),20,0,20);
    addHisto("hBDT_drSig","TRK_dm300",                sample.Data(),1000,0,2000);
    addHisto("hBDT_dzSig","TRK_dm300",                sample.Data(),1000,0,2000);
    addHisto("hBDT_ntrk10","TRK_dm300",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk20","TRK_dm300",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk30","TRK_dm300",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk40","TRK_dm300",               sample.Data(),100,0,100);
    addHisto("hBDT_Hemi_dR","TRK_dm300",              sample.Data(),100,0,6);
    addHisto("hBDT_Hemi_dRmax","TRK_dm300",           sample.Data(),100,0,6);

            addHisto("hBDT_lost","TRK_dm320",                 sample.Data(),2,0,2);
    addHisto("hBDT_dz","TRK_dm320",                   sample.Data(),400,-200,200);
    addHisto("hBDT_dxy","TRK_dm320",                  sample.Data(),200,-100,100);
    addHisto("hBDT_pt","TRK_dm320",                   sample.Data(),100,0,200);
    addHisto("hBDT_eta","TRK_dm320",                  sample.Data(),80,-4,4);
    addHisto("hBDT_NChi2","TRK_dm320",                sample.Data(),5,0,5);
    addHisto("hBDT_nhits","TRK_dm320",                sample.Data(),35,0,35);
    addHisto("hBDT_iJet","TRK_dm320",                 sample.Data(),20,0,20);
    addHisto("hBDT_drSig","TRK_dm320",                sample.Data(),1000,0,2000);
    addHisto("hBDT_dzSig","TRK_dm320",                sample.Data(),1000,0,2000);
    addHisto("hBDT_ntrk10","TRK_dm320",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk20","TRK_dm320",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk30","TRK_dm320",               sample.Data(),100,0,100);
    addHisto("hBDT_ntrk40","TRK_dm320",               sample.Data(),100,0,100);
    addHisto("hBDT_Hemi_dR","TRK_dm320",              sample.Data(),100,0,6);
    addHisto("hBDT_Hemi_dRmax","TRK_dm320",           sample.Data(),100,0,6);


    // ---------- zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz--//

    addHisto("hBDT_nTrks","VTX",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX",                  sample.Data(),1001,1.5,2.5);
    addHisto("hBDT_Vtx_Mass","VTX",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX",              sample.Data(),100,0,1);

    addHisto("hBDT_nTrks","VTX_TightWP",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX_TightWP",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX_TightWP",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX_TightWP",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX_TightWP",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX_TightWP",                  sample.Data(),1001,1.5,2.5);
    addHisto2D("MWTVsnTrks","VTX_TightWP",              sample.Data(),50,0,50,50,0,50);
    addHisto2D("MWTVsnTrks","VTX",                      sample.Data(),50,0,50,50,0,50);
    addHisto("hBDT_Vtx_Mass","VTX_TightWP",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX_TightWP",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX_TightWP",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX_TightWP",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX_TightWP",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX_TightWP",              sample.Data(),100,0,1);

  // Eat sleep and repeat -----------------------------------------------------

      addHisto("hBDT_nTrks","VTX_TightWP_dm20",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX_TightWP_dm20",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX_TightWP_dm20",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX_TightWP_dm20",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX_TightWP_dm20",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX_TightWP_dm20",                  sample.Data(),1001,1.5,2.5);
    addHisto("hBDT_Vtx_Mass","VTX_TightWP_dm20",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX_TightWP_dm20",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX_TightWP_dm20",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX_TightWP_dm20",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX_TightWP_dm20",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX_TightWP_dm20",              sample.Data(),100,0,1);

    
    addHisto("hBDT_nTrks","VTX_TightWP_dm50",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX_TightWP_dm50",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX_TightWP_dm50",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX_TightWP_dm50",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX_TightWP_dm50",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX_TightWP_dm50",                  sample.Data(),1001,1.5,2.5);
    addHisto("hBDT_Vtx_Mass","VTX_TightWP_dm50",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX_TightWP_dm50",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX_TightWP_dm50",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX_TightWP_dm50",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX_TightWP_dm50",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX_TightWP_dm50",              sample.Data(),100,0,1);


    addHisto("hBDT_nTrks","VTX_TightWP_dm100",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX_TightWP_dm100",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX_TightWP_dm100",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX_TightWP_dm100",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX_TightWP_dm100",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX_TightWP_dm100",                  sample.Data(),1001,1.5,2.5);
    addHisto("hBDT_Vtx_Mass","VTX_TightWP_dm100",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX_TightWP_dm100",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX_TightWP_dm100",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX_TightWP_dm100",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX_TightWP_dm100",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX_TightWP_dm100",              sample.Data(),100,0,1);

    addHisto("hBDT_nTrks","VTX_TightWP_dm150",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX_TightWP_dm150",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX_TightWP_dm150",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX_TightWP_dm150",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX_TightWP_dm150",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX_TightWP_dm150",                  sample.Data(),1001,1.5,2.5);
    addHisto("hBDT_Vtx_Mass","VTX_TightWP_dm150",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX_TightWP_dm150",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX_TightWP_dm150",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX_TightWP_dm150",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX_TightWP_dm150",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX_TightWP_dm150",              sample.Data(),100,0,1);

    addHisto("hBDT_nTrks","VTX_TightWP_dm200",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX_TightWP_dm200",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX_TightWP_dm200",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX_TightWP_dm200",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX_TightWP_dm200",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX_TightWP_dm200",                  sample.Data(),1001,1.5,2.5);
    addHisto("hBDT_Vtx_Mass","VTX_TightWP_dm200",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX_TightWP_dm200",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX_TightWP_dm200",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX_TightWP_dm200",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX_TightWP_dm200",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX_TightWP_dm200",              sample.Data(),100,0,1);

    addHisto("hBDT_nTrks","VTX_TightWP_dm250",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX_TightWP_dm250",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX_TightWP_dm250",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX_TightWP_dm250",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX_TightWP_dm250",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX_TightWP_dm250",                  sample.Data(),1001,1.5,2.5);
    addHisto("hBDT_Vtx_Mass","VTX_TightWP_dm250",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX_TightWP_dm250",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX_TightWP_dm250",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX_TightWP_dm250",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX_TightWP_dm250",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX_TightWP_dm250",              sample.Data(),100,0,1);

    addHisto("hBDT_nTrks","VTX_TightWP_dm300",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX_TightWP_dm300",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX_TightWP_dm300",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX_TightWP_dm300",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX_TightWP_dm300",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX_TightWP_dm300",                  sample.Data(),1001,1.5,2.5);
    addHisto("hBDT_Vtx_Mass","VTX_TightWP_dm300",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX_TightWP_dm300",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX_TightWP_dm300",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX_TightWP_dm300",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX_TightWP_dm300",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX_TightWP_dm300",              sample.Data(),100,0,1);

    addHisto("hBDT_nTrks","VTX_TightWP_dm320",                sample.Data(),50,2,52);
    addHisto("hBDT_NChi2","VTX_TightWP_dm320",                sample.Data(),40,0,10);
    addHisto("hBDT_step","VTX_TightWP_dm320",                 sample.Data(),5,0,5);
    addHisto("hBDT_r","VTX_TightWP_dm320",                    sample.Data(),100,0,100);
    addHisto("hBDT_z","VTX_TightWP_dm320",                    sample.Data(),200,0,200);
    addHisto("hBDT_MWT","VTX_TightWP_dm320",                  sample.Data(),1001,1.5,2.5);
    addHisto("hBDT_Vtx_Mass","VTX_TightWP_dm320",             sample.Data(),1400,0,14000);
    addHisto("hBDT_Hemi_Mass","VTX_TightWP_dm320",            sample.Data(),140,0,1400);
    addHisto("hBDT_dist","VTX_TightWP_dm320",                 sample.Data(),200,0,200);
    addHisto("hBDT_ntrk10","VTX_TightWP_dm320",               sample.Data(),50,0,50);
    addHisto("hBDT_ntrk20","VTX_TightWP_dm320",               sample.Data(),50,0,50);
    addHisto("hBDT_MeanDCA","VTX_TightWP_dm320",              sample.Data(),100,0,1);

    // --------------- zzzzzzzzzzzzzzzzzzzzzzzzzzzz ------------------------------------//
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

void TreeReader::addHisto2D(TString var, TString selstep, TString sample, int nxbins, float xmin, float xmax, int nybins, float ymin, float ymax){
 
  TString name =  var+"_"+selstep+"__"+sample;
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
void TreeReader::fillHisto( TString var, TString selstep, TString sample, float val, float weight){
  TString name = var+"_"+selstep+"__"+sample;


  if(histo_map_[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histograms " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }else  histo_list_[histo_map_[name.Data()]]->Fill(val, weight);
  
}


void TreeReader::fillHisto2D( TString var, TString selstep, TString sample, float xval,float yval, float weight){
  TString name = var+"_"+selstep+"__"+sample;


  if(histo_map_2D_[name.Data()] == 0) {
    cout << "   WARNING trying to fill a non existing histograms " << endl;
    cout << "   please check the naming conventions " << endl;
    cout << "   histo name "  << name << endl;
  }else  histo_list_2D_[histo_map_2D_[name.Data()]]->Fill(xval,yval, weight);
  
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
