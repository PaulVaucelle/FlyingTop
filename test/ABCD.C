#define ABCD_cxx
#include "ABCD.h"
#include "TH1F.h"
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

void SetSebStyle()
{
 gStyle->SetTitleFillColor(42);
 gStyle->SetTitleFont(1);
 gStyle->SetStatColor(29);
 gStyle->SetCanvasColor(25);   
 gStyle->SetOptStat(1111111);
 gStyle->SetHistFillColor(5);
}

void SaveInFile(TH1* ahisto, TFile* afile)
{
 if (!ahisto) { std::cout <<     "!! no histo !!" << std::endl; return ;}
 TDirectory* current = gDirectory ;
 afile->cd();
 ahisto->Write();
 current->cd();
}

void ABCD::Loop(int aNN, float aTagCut, float aPtMin, float aPtMax, 
               float aEtaMin, float aEtaMax, float aFreeCut, int aIntCut, 
		         TString afilename, TString aweightFileMVA,
               TString sample, TString Production, bool Signal)
{
//   In a ROOT session, you can do:
//      root> .L ABCD.C
//      root> ABCD t
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

//$$
  int    NN           = aNN;      // 0/1 JP, 2/3 TC, 4/5 SSV, 6 CSV, 7 SL
  float  TagCut       = aTagCut;  // tag cut
  float  PtMin        = aPtMin;   // pt jet min
  float  PtMax        = aPtMax;   // pt jet max
  float  EtaMin       = aEtaMin;  // eta jet min
  float  EtaMax       = aEtaMax;  // eta jet max
  float  FreeCut      = aFreeCut; 
  int    IntCut       = aIntCut; 
  TString  filename   = afilename;
  TString  weightFile_ = aweightFileMVA;
//$$

//$$
  float HTcut =  aTagCut; 
//   float dRcut =  FreeCut; 
//$$
//**********************************
// Histograms
//**********************************
 TH1F* hData_Filter     = new TH1F("hData_Filter","",2,-0.5,1.5);

   //-----------------------------------------------------------//
   // ABCD using Evt and Tight+looseWP 
   //-----------------------------------------------------------//

   //------SR-1vtx----//
      TH1F* hData_EVT12_1Vtx_CutEvt_Mmumu    = new TH1F("hData_EVT12_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVT12_1Vtx_CutEvt_Mass     = new TH1F("hData_EVT12_1Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT12_1Vtx_BDTvtx          = new TH1F("hData_EVT12_1Vtx_BDTvtx","",100,-1,1);
   //------CR-----//
      TH1F* hData_NoEVT12_1Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVT12_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVT12_1Vtx_CutEvt_Mass   = new TH1F("hData_NoEVT12_1Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT12_1Vtx_BDTvtx        = new TH1F("hData_NoEVT12_1Vtx_BDTvtx","",100,-1,1);
   //------CR-----//
      TH1F* hData_EVT34_1Vtx_CutEvt_Mmumu    = new TH1F("hData_EVT34_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVT34_1Vtx_CutEvt_Mass     = new TH1F("hData_EVT34_1Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT34_1Vtx_BDTvtx          = new TH1F("hData_EVT34_1Vtx_BDTvtx","",100,-1,1);
   //------CR-----//
      TH1F* hData_NoEVT34_1Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVT34_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVT34_1Vtx_CutEvt_Mass   = new TH1F("hData_NoEVT34_1Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT34_1Vtx_BDTvtx        = new TH1F("hData_NoEVT34_1Vtx_BDTvtx","",100,-1,1);
   //------SR-2Vtx----//
      TH1F* hData_EVT12_2Vtx_CutEvt_Mmumu    = new TH1F("hData_EVT12_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVT12_2Vtx_CutEvt_Mass     = new TH1F("hData_EVT12_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT12_2Vtx_BDTvtx          = new TH1F("hData_EVT12_2Vtx_BDTvtx","",100,-1,1);
   //------CR-----//
      TH1F* hData_NoEVT12_2Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVT12_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVT12_2Vtx_CutEvt_Mass   = new TH1F("hData_NoEVT12_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT12_2Vtx_BDTvtx        = new TH1F("hData_NoEVT12_2Vtx_BDTvtx","",100,-1,1);
   //------CR-----//
      TH1F* hData_EVT34_2Vtx_CutEvt_Mmumu    = new TH1F("hData_EVT34_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVT34_2Vtx_CutEvt_Mass     = new TH1F("hData_EVT34_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVT34_2Vtx_BDTvtx          = new TH1F("hData_EVT34_2Vtx_BDTvtx","",100,-1,1);
   //------CR-----//
      TH1F* hData_NoEVT34_2Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVT34_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVT34_2Vtx_CutEvt_Mass   = new TH1F("hData_NoEVT34_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVT34_2Vtx_BDTvtx        = new TH1F("hData_NoEVT34_2Vtx_BDTvtx","",100,-1,1);


   //-----------------------------------------------------------//
   // ABCD using Evt and Vtx BDT
   //-----------------------------------------------------------//
      TH1F* hData_1Vtx_MVAval                = new TH1F("hData_1Vtx_MVAval","",,100,-1,1);
      TH1F* hData_2Vtx_MVAval                = new TH1F("hData_2Vtx_MVAval","",,100,-1,1);
      TH1F* hData_2VtxAll_MVAval             = new TH1F("hData_2VtxAll_MVAval","",100,-1,1);
      TH1F* hData_2VtxAll_MVAval             = new TH1F("hData_2VtxAll_MVAval","",100,-1,1);

   //------SR-1vtx----//
      TH1F* hData_EVTVtx_1Vtx_CutEvt_Mmumu   = new TH1F("hData_EVTVtx_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTVtx_1Vtx_CutEvt_Mass    = new TH1F("hData_EVTVtx_1Vtx_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_NoEVTVtx_1Vtx_CutEvt_Mmumu = new TH1F("hData_NoEVTVtx_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTVtx_1Vtx_CutEvt_Mass  = new TH1F("hData_NoEVTVtx_1Vtx_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_EVTNoVtx_1Vtx_CutEvt_Mmumu = new TH1F("hData_EVTNoVtx_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTNoVtx_1Vtx_CutEvt_Mass  = new TH1F("hData_EVTNoVtx_1Vtx_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_NoEVTNoVtx_1Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVTNoVtx_1Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTNoVtx_1Vtx_CutEvt_Mass   = new TH1F("hData_NoEVTNoVtx_1Vtx_CutEvt_Mass","",25,0.,100.);

   //------SR-2Vtx----//
      TH1F* hData_EVTVtx_2Vtx_CutEvt_Mmumu      = new TH1F("hData_EVTVtx_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTVtx_2Vtx_CutEvt_Mass       = new TH1F("hData_EVTVtx_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVTVtx_2VtxAll_CutEvt_Mmumu   = new TH1F("hData_EVTVtx_2VtxAll_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTVtx_2VtxAll_CutEvt_Mass    = new TH1F("hData_EVTVtx_2VtxAll_CutEvt_Mass","",25,0.,100.);
         
   //------CR-----//
      TH1F* hData_NoEVTVtx_2Vtx_CutEvt_Mmumu    = new TH1F("hData_NoEVTVtx_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTVtx_2Vtx_CutEvt_Mass     = new TH1F("hData_NoEVTVtx_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVTVtx_2VtxAll_CutEvt_Mmumu = new TH1F("hData_NoEVTVtx_2VtxAll_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTVtx_2VtxAll_CutEvt_Mass  = new TH1F("hData_NoEVTVtx_2VtxAll_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_EVTNoVtx_2Vtx_CutEvt_Mmumu    = new TH1F("hData_EVTNoVtx_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTNoVtx_2Vtx_CutEvt_Mass     = new TH1F("hData_EVTNoVtx_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_EVTNoVtx_2VtxAll_CutEvt_Mmumu = new TH1F("hData_EVTNoVtx_2VtxAll_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_EVTNoVtx_2VtxAll_CutEvt_Mass  = new TH1F("hData_EVTNoVtx_2VtxAll_CutEvt_Mass","",25,0.,100.);
   //------CR-----//
      TH1F* hData_NoEVTNoVtx_2Vtx_CutEvt_Mmumu  = new TH1F("hData_NoEVTNoVtx_2Vtx_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTNoVtx_2Vtx_CutEvt_Mass   = new TH1F("hData_NoEVTNoVtx_2Vtx_CutEvt_Mass","",25,0.,100.);
      TH1F* hData_NoEVTNoVtx_2VtxAll_CutEvt_Mmumu  = new TH1F("hData_NoEVTNoVtx_2VtxAll_CutEvt_Mmumu","",25,0.,500.);
      TH1F* hData_NoEVTNoVtx_2VtxAll_CutEvt_Mass   = new TH1F("hData_NoEVTNoVtx_2VtxAll_CutEvt_Mass","",25,0.,100.);
         
   //-----------------------------------------------------------//
   // ABCD using Hemisphere pt and Tight+looseWP
   //-----------------------------------------------------------//


//------SR-----//
   TH1F* hData_Hemi_0Vtx_CutEvt_Mmumu       = new TH1F("hData_Hemi_0Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_Hemi_1Vtx_CutEvt_Mass    = new TH1F("hData_Hemi_1Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_Hemi_1Vtx_CutEvt_Mmumu   = new TH1F("hData_Hemi_1Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_Hemi_2Vtx_CutEvt_Mass    = new TH1F("hData_Hemi_2Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_Hemi_2Vtx_CutEvt_Mmumu   = new TH1F("hData_Hemi_2Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_Hemi_2VtxAll_CutEvt_Mass    = new TH1F("hData_Hemi_2VtxAll_CutEvt_Mass","",25,0.,100.);
//------low pt-----//
   TH1F* hData_CRlowpt_0Vtx_CutEvt_Mmumu       = new TH1F("hData_CRlowpt_0Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRlowpt_1Vtx_CutEvt_Mass    = new TH1F("hData_CRlowpt_1Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRlowpt_1Vtx_CutEvt_Mmumu   = new TH1F("hData_CRlowpt_1Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRlowpt_2Vtx_CutEvt_Mass    = new TH1F("hData_CRlowpt_2Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRlowpt_2Vtx_CutEvt_Mmumu   = new TH1F("hData_CRlowpt_2Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRlowpt_2VtxAll_CutEvt_Mass    = new TH1F("hData_CRlowpt_2VtxAll_CutEvt_Mass","",25,0.,100.);
//------loose -----//
   TH1F* hData_CRloose_0Vtx_CutEvt_Mmumu       = new TH1F("hData_CRloose_0Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRloose_1Vtx_CutEvt_Mass    = new TH1F("hData_CRloose_1Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRloose_1Vtx_CutEvt_Mmumu   = new TH1F("hData_CRloose_1Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRloose_2Vtx_CutEvt_Mass    = new TH1F("hData_CRloose_2Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRloose_2Vtx_CutEvt_Mmumu   = new TH1F("hData_CRloose_2Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRloose_2VtxAll_CutEvt_Mass    = new TH1F("hData_CRloose_2VtxAll_CutEvt_Mass","",25,0.,100.);

//------loose low pt -----//
   TH1F* hData_CRlooselowpt_0Vtx_CutEvt_Mmumu       = new TH1F("hData_CRlooselowpt_0Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRlooselowpt_1Vtx_CutEvt_Mass    = new TH1F("hData_CRlooselowpt_1Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRlooselowpt_1Vtx_CutEvt_Mmumu   = new TH1F("hData_CRlooselowpt_1Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRlooselowpt_2Vtx_CutEvt_Mass    = new TH1F("hData_CRlooselowpt_2Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRlooselowpt_2Vtx_CutEvt_Mmumu   = new TH1F("hData_CRlooselowpt_2Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRlooselowpt_2VtxAll_CutEvt_Mass    = new TH1F("hData_CRlooselowpt_2VtxAll_CutEvt_Mass","",25,0.,100.);
//------ fwd -----//
   TH1F* hData_CRfwd_0Vtx_CutEvt_Mmumu       = new TH1F("hData_CRfwd_0Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRfwd_1Vtx_CutEvt_Mass    = new TH1F("hData_CRfwd_1Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRfwd_1Vtx_CutEvt_Mmumu   = new TH1F("hData_CRfwd_1Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRfwd_2Vtx_CutEvt_Mass    = new TH1F("hData_CRfwd_2Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRfwd_2Vtx_CutEvt_Mmumu   = new TH1F("hData_CRfwd_2Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRfwd_2VtxAll_CutEvt_Mass    = new TH1F("hData_CRfwd_2VtxAll_CutEvt_Mass","",25,0.,100.);
   //----- Samesign ----//  
   TH1F* hData_CRss_0Vtx_CutEvt_Mmumu       = new TH1F("hData_CRss_0Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRss_1Vtx_CutEvt_Mass    = new TH1F("hData_CRss_1Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRss_1Vtx_CutEvt_Mmumu   = new TH1F("hData_CRss_1Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRss_2Vtx_CutEvt_Mass    = new TH1F("hData_CRss_2Vtx_CutEvt_Mass","",25,0.,100.);
   TH1F* hData_CRss_2Vtx_CutEvt_Mmumu   = new TH1F("hData_CRss_2Vtx_CutEvt_Mmumu","",25,0.,500.);
   TH1F* hData_CRss_2VtxAll_CutEvt_Mass    = new TH1F("hData_CRss_2VtxAll_CutEvt_Mass","",25,0.,100.);

///////////////////////////////////////////////////////////////////


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

  float XS = 1;
   TString thesample  = sample;
   // XS are given in pb
   if (thesample.Contains("DYJetsToLL_M-10to50"))                    { XS = 15910.0;   }
   if (thesample.Contains("ST_tW_antitop_5f_NoFullyHadronicDecays")) { XS = 10.8707;   }
   if (thesample.Contains("ST_tW_top_5f_NoFullyHadronicDecays"))     { XS = 10.8908;   }
   if (thesample.Contains("TTJets_DiLept"))                          { XS = 53.07;     }
   if (thesample.Contains("WWTo2L2Nu"))                              { XS = 11.09;     }
   if (thesample.Contains("WZTo2Q2L_mllmin4p0"))                     { XS = 6.535;     }
   if (thesample.Contains("ZZTo2Q2L_mllmin4p0"))                     { XS = 3.676;     }
   if (thesample.Contains("TTTo2L2Nu"))                              { XS = 88.3;      }
   if (thesample.Contains("DYJetsToLL_M-50"))                        { XS = 5379;      }
   if (thesample.Contains("ttWJetsToLNu_5f_EWK"))                    { XS = 0.290;     } // not found on XSDB, no file on tier2...approximation
         //Took 0.868 pb (CMS-TOP-21-011)
      // as a starting point and then divided by 3 (lepton universality)
   if (thesample.Contains("TTZToLL_5f"))                             { XS = 0.05188;   }//Not found on XSDB => used ana.py macro : 5.188e-02 +- 2.437e-04 pb
   if (thesample.Contains("TTToHadronic"))                           { XS = 378.9;     }
   if (thesample.Contains("TTWW"))                                   { XS = 0.006992;  }//found on XSDB
   if (thesample.Contains("TTToSemiLeptonic") )                      { XS = 365.34;    }

   //Signal
   if (thesample.Contains("smu200")) { XS = 0.01;   }
   if (thesample.Contains("smu250")) { XS = 0.0045; }
   if (thesample.Contains("smu300")) { XS = 0.002;  }
   if (thesample.Contains("smu350")) { XS = 0.001;  }
   if (thesample.Contains("smu400")) { XS = 0.0006; }
   if (thesample.Contains("smu450")) { XS = 0.0004; }
   if (thesample.Contains("smu500")) { XS = 0.00025;}
   if (thesample.Contains("70_test")){ XS = 0.0035; } //this is for 14tev
   if (thesample.Contains("50_test")){ XS = 0.0025; } //this is for 14tev 
   if (thesample.Contains("30_test")){ XS = 0.002;  } //this is for 14tev 
   if (thesample.Contains("10_test")){ XS = 0.0035; } //this is for 14tev

   bool BlindSR = false;
   float Nevent = 0;
   bool signal = Signal;
   //   if (nentries>760000){nentries =760000;}

   cout<< "Line : "  << __LINE__ << " " << nentries << endl; 
   cout<< " XS : "<<XS<<endl;
   Long64_t nbytes = 0, nb = 0;
   int allevents = 0;
   // nentries = 100000;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      allevents++;
      if ( allevents%1000 == 0 ) std::cout << "events : " << allevents << std::endl;

       if ( signal &&  tree_nLLP != 2 ) continue; // protection against rare wrong signal events    // 
  
      //--------------------------------------------------------------//
      if (jentry==0.1*nentries) {std::cout<<"10/100 :"<<std::endl;}
      if (jentry==0.2*nentries) {std::cout<<"20/100 :"<<std::endl;}
      if (jentry==0.3*nentries) {std::cout<<"30/100 :"<<std::endl;}
      if (jentry==0.4*nentries) {std::cout<<"40/100 :"<<std::endl;}
      if (jentry==0.5*nentries) {std::cout<<"50/100 :"<<std::endl;}
      if (jentry==0.6*nentries) {std::cout<<"60/100 :"<<std::endl;}
      if (jentry==0.7*nentries) {std::cout<<"70/100 :"<<std::endl;}
      if (jentry==0.8*nentries) {std::cout<<"80/100 :"<<std::endl;}
      if (jentry==0.9*nentries) {std::cout<<"90/100 :"<<std::endl;}

      hData_Filter->Fill( tree_Filter );

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

      bool ping0 = false;
      bool ping1 = false;
      float dR0 = 0.;
      float dR1 = 0.;
      ////////////////// hemisphere pT
      float  hemi1_pt  = -1.;
      float  hemi2_pt  = -1.;
      ////////////////////////////////

      hData_Mmumu->Fill( tree_Mmumu );
      
      if ( !tree_Filter && !tree_FilterSameSign ) continue;//
      if (tree_njetNOmu < 1) continue;

      bool Filter = tree_Filter; // Tight ID and Mini IsoTight for both muons
      if ( tree_njetNOmu < 1 ) Filter = false; 
      hemi1_pt = tree_Hemi_pt->at(0);
      hemi2_pt = tree_Hemi_pt->at(1);

      float hemi_ptmin = hemi2_pt;
      if ( hemi2_pt > hemi1_pt ) hemi_ptmin = hemi1_pt;


   if (signal)
      {
         ping0 = tree_Hemi_LLP_ping->at(0);
         ping1 = tree_Hemi_LLP_ping->at(1);
         dR0 = tree_Hemi_LLP_dR->at(0);
         dR1 = tree_Hemi_LLP_dR->at(1);
      }

   int Vtx_step0 = tree_Hemi_Vtx_step->at(0);
   int Vtx_step1 = tree_Hemi_Vtx_step->at(1);
   float Vtx_NChi0 = tree_Hemi_Vtx_NChi2->at(0);
   float Vtx_NChi1 = tree_Hemi_Vtx_NChi2->at(1);
   float Vtx_Mass0 = tree_Hemi_Vtx_Mass->at(0);
   float Vtx_Mass1 = tree_Hemi_Vtx_Mass->at(1);
   float Vtx_dist0 = tree_Hemi_Vtx_dist->at(0);
   float Vtx_dist1 = tree_Hemi_Vtx_dist->at(1);

   float posx0 = tree_Hemi_Vtx_x->at(0);
   float posy0 = tree_Hemi_Vtx_y->at(0);
   float posz0 = tree_Hemi_Vtx_z->at(0);
   float posx1 = tree_Hemi_Vtx_x->at(1);
   float posy1 = tree_Hemi_Vtx_y->at(1);
   float posz1 = tree_Hemi_Vtx_z->at(1);
   float r0 = TMath::Sqrt( posx0*posx0 + posy0*posy0 );
   float z0 = TMath::Abs( posz0 );
   float r1 = TMath::Sqrt( posx1*posx1 + posy1*posy1 );
   float z1 = TMath::Abs( posz1 );
   float recX0 = posx0 - tree_PV_x;
   float recY0 = posy0 - tree_PV_y;
   float recZ0 = posz0 - tree_PV_z;
   float recX1 = posx1 - tree_PV_x;
   float recY1 = posy1 - tree_PV_y;
   float recZ1 = posz1 - tree_PV_z;

  float theta_Vtx0 = TMath::ATan2(sqrt(recX0*recX0+recY0*recY0),abs(recZ0)) ;
  float theta_Vtx1 = TMath::ATan2(sqrt(recX1*recX1+recY1*recY1),abs(recZ1)) ;

  float eta_Vtx0 = -TMath::Log(tan(theta_Vtx0/2));
  float eta_Vtx1 = -TMath::Log(tan(theta_Vtx1/2));
  if ( posz0 < 0 ) eta_Vtx0 = -eta_Vtx0;
  if ( posz1 < 0 ) eta_Vtx1 = -eta_Vtx1;

   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 1 && Vtx_step0 <= 2 ) {
     isHemiVtx1 = true;
   }
   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 3 && Vtx_step0 <= 4 ) {
     isHemiVtx1Loose = true;
   }
   if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 1 && Vtx_step1 <= 2 ) {
     isHemiVtx2 = true;
   }
   if ( Vtx_NChi1 > 0 && Vtx_NChi1< 10 && Vtx_step1 >= 3 && Vtx_step1 <= 4 ) {
     isHemiVtx2Loose = true;
   }
   if      ( isHemiVtx1 && isHemiVtx2 ) nVtxIni = 2;
   else if ( isHemiVtx1 || isHemiVtx2 ) nVtxIni = 1;

   if      ( isHemiVtx1Loose && isHemiVtx2Loose ) nVtxIniLoose = 2;
   else if ( isHemiVtx1Loose || isHemiVtx2Loose ) nVtxIniLoose = 1;

   bool Merging = true; // false if "no merge" or "no close vtx" output
   bool Protect = false;
   //$$
   if (signal)
      {
         if (tree_Hemi_SecLLP_ping->size() >= 1) Protect = true;
      }
   else
      {
         Protect = true;
      } 
    if ( Merging && Protect ) { 
         // protection again
     if ( tree_Hemi_SecVtx->size() >= 1 ) {
         if (signal)
         {ping0 = tree_Hemi_SecLLP_ping->at(0);}
       ping0 = false;
       Vtx_step0 = tree_Hemi_SecVtx_step->at(0);
       Vtx_NChi0 = tree_Hemi_SecVtx_NChi2->at(0);
       Vtx_Mass0 = tree_Hemi_SecVtx_Mass->at(0);
       Vtx_dist0 = tree_Hemi_SecVtx_dist->at(0);
       posx0 = tree_Hemi_SecVtx_x->at(0);
       posy0 = tree_Hemi_SecVtx_y->at(0);
       posz0 = tree_Hemi_SecVtx_z->at(0);
       r0 = tree_Hemi_SecVtx_r->at(0);
       Vtx_step1 = 0;
       Vtx_NChi1 = -1.;
       Vtx_Mass1 = 0.;
       Vtx_dist1 = 0.;
       ping1 = false;
       r1 = 0;
       float SecrecX0 = posx0 - tree_PV_x;
       float SecrecY0 = posy0 - tree_PV_y;
       float SecrecZ0 = posz0 - tree_PV_z;
       float theta_SecVtx0 = TMath::ATan2(sqrt(SecrecX0*SecrecX0+SecrecY0*SecrecY0),abs(SecrecZ0)) ;

       eta_Vtx0 = -TMath::Log(tan(theta_SecVtx0/2.));
       eta_Vtx1 = 0;
       if ( posz0 < 0 ) eta_Vtx0 = -eta_Vtx0;

     }
     if ( tree_Hemi_SecVtx->size() == 2 ) {
               if (signal)
         {ping1 = tree_Hemi_SecLLP_ping->at(1);}
       ping1 = false;
       Vtx_step1 = tree_Hemi_SecVtx_step->at(1);
       Vtx_NChi1 = tree_Hemi_SecVtx_NChi2->at(1);
       Vtx_Mass1 = tree_Hemi_SecVtx_Mass->at(1);
       Vtx_dist1 = tree_Hemi_SecVtx_dist->at(1);
       posx1 = tree_Hemi_SecVtx_x->at(1);
       posy1 = tree_Hemi_SecVtx_y->at(1);
       posz1 = tree_Hemi_SecVtx_z->at(1);
       r1 = tree_Hemi_SecVtx_r->at(1);
       float SecrecX1 = posx1 - tree_PV_x;
       float SecrecY1 = posy1 - tree_PV_y;
       float SecrecZ1 = posz1 - tree_PV_z;
       float theta_SecVtx1 = TMath::ATan2(sqrt(SecrecX1*SecrecX1+SecrecY1*SecrecY1),abs(SecrecZ1)) ;
       eta_Vtx1 = -TMath::Log(tan(theta_SecVtx1/2.));
       if ( posz1 < 0 ) eta_Vtx1 = -eta_Vtx1;
     }
   }
   isHemiVtx1 = false;
   isHemiVtx2 = false;

   isHemiVtx1Loose = false;
   isHemiVtx2Loose = false;        

   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 1 && Vtx_step0 <= 2 ) {
     isHemiVtx1 = true;
     VtxMass = Vtx_Mass0;
     BDTvtx1 = tree_Hemi_Vtx_MVAval_Step1->at(0);
     BDTvtx  = BDTvtx1;
   }
   if ( Vtx_NChi0 > 0 && Vtx_NChi0 < 10 && Vtx_step0 >= 3 && Vtx_step0 <= 4 ) {
     isHemiVtx1Loose = true;
     VtxMass = Vtx_Mass0;
     BDTvtx1 = tree_Hemi_Vtx_MVAval->at(0);
     BDTvtx  = BDTvtx1;
   }
   if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 1 && Vtx_step1 <= 2 ) {
     isHemiVtx2 = true;
     if ( Vtx_Mass1 > VtxMass ) VtxMass = Vtx_Mass1;
     BDTvtx2 = tree_Hemi_Vtx_MVAval_Step1->at(1);
     if ( BDTvtx2 > BDTvtx ) BDTvtx = BDTvtx2;
   }
   if ( Vtx_NChi1 > 0 && Vtx_NChi1 < 10 && Vtx_step1 >= 3 && Vtx_step1 <= 4 ) {
     isHemiVtx2Loose = true;
     if ( Vtx_Mass1 > VtxMass ) VtxMass = Vtx_Mass1;
     BDTvtx2 = tree_Hemi_Vtx_MVAval->at(1);
     if ( BDTvtx2 > BDTvtx ) BDTvtx = BDTvtx2;
   }
   if      ( isHemiVtx1 && isHemiVtx2 ) nVtx = 2;
   else if ( isHemiVtx1 || isHemiVtx2 ) nVtx = 1;

   if      ( isHemiVtx1Loose && isHemiVtx2Loose ) nVtxLoose = 2;
   else if ( isHemiVtx1Loose || isHemiVtx2Loose ) nVtxLoose = 1;

   int nHemi = tree_Hemi->size();//_LLP_dR
   for (int i = 0; i < nHemi; i++) 
      {   // Loop on hemispheres
         if ( i == 0 ) {
            dR    = dR0;
            dist  = Vtx_dist0;
            NChi2 = Vtx_NChi0;
            step  = Vtx_step0;
            ping  = ping0;
            r = r0;
            eta = eta_Vtx0;
         }
         else if ( i == 1 ) {
            dR    = dR1;
            dist  = Vtx_dist1;
            NChi2 = Vtx_NChi1;
            step  = Vtx_step1;
            ping  = ping1;
            r = r1;
            eta= eta_Vtx1;
         }
      }// end loop on Hemi

   //-----------------------------------------------------------//
   // ABCD using Evt and Tight+looseWP 
   //-----------------------------------------------------------//
   float EVTSWP = 0;
   float VTXWP = 0.5;
   hData_Evt_MVAval->Fill( tree_Evts_MVAval );
   //$$
   if (tree_Evts_MVAval > EVTSWP) isCutEvt = true;

   //$$
   //----------------------//
   if (nVtx == 1 && isCutEvt ) {
      hData_EVT12_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_EVT12_1Vtx_CutEvt_Mass->Fill(VtxMass);
      hData_EVT12_1Vtx_BDTvtx->Fill( BDTvtx );
   }

   if (nVtx == 1 && !isCutEvt) {
      hData_NoEVT12_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_NoEVT12_1Vtx_CutEvt_Mass->Fill(VtxMass);
      hData_NoEVT12_1Vtx_BDTvtx->Fill( BDTvtx );
   }

   if (nVtxLoose == 1 && isCutEvt) {
      hData_EVT34_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_EVT34_1Vtx_CutEvt_Mass->Fill(VtxMass);
      hData_EVT34_1Vtx_BDTvtx->Fill( BDTvtx );
   }

   if (nVtxLoose == 1 && !isCutEvt) {
      hData_NoEVT34_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_NoEVT34_1Vtx_CutEvt_Mass->Fill(VtxMass);
      hData_NoEVT34_1Vtx_BDTvtx->Fill( BDTvtx );
   }
   //----------------------//

   if (nVtx == 2 && isCutEvt) {
      hData_EVT12_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_EVT12_2Vtx_CutEvt_Mass->Fill(VtxMass);
      hData_EVT12_2Vtx_BDTvtx->Fill( BDTvtx );
   }

   if (nVtx == 2 && !isCutEvt) {
      hData_NoEVT12_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_NoEVT12_2Vtx_CutEvt_Mass->Fill(VtxMass);
      hData_NoEVT12_2Vtx_BDTvtx->Fill( BDTvtx );
   }

   if (nVtxLoose == 2 && isCutEvt) {
      hData_EVT34_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_EVT34_2Vtx_CutEvt_Mass->Fill(VtxMass);
      hData_EVT34_2Vtx_BDTvtx->Fill( BDTvtx );
   }

   if (nVtxLoose == 2 && !isCutEvt) {
      hData_NoEVT34_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_NoEVT34_2Vtx_CutEvt_Mass->Fill(VtxMass);
      hData_NoEVT34_2Vtx_BDTvtx->Fill( BDTvtx );
   }
    //----------------------//
   //$$

   //-----------------------------------------------------------//
   // ABCD using Evt and Vtx BDT
   //-----------------------------------------------------------//
   isCutEvt = false;
   isCutVtx = false;
   if (tree_Evts_MVAval > EVTSWP) isCutEvt = true;
   if (BDTvtx > VTXWP) isCutVtx = true;
   if (BDTvtx1 > VTXWP) isCutVtx1 = true;
   if (BDTvtx2 > VTXWP) isCutVtx2 = true;

   if (nVtx == 1 && isCutEvt)
      {
         hData_1Vtx_MVAval->Fill(BDTvtx);
      }
   if (nVtx == 2 && isCutEvt)
      {
         hData_2Vtx_MVAval->Fill(BDTvtx);
         hData_2VtxAll_MVAval->Fill(BDTvtx1);
         hData_2VtxAll_MVAval->Fill(BDTvtx2);
      }
   //----------------------//

   if (nVtx == 1 && isCutEvt  && isCutVtx) {
      hData_EVTVtx_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_EVTVtx_1Vtx_CutEvt_Mass->Fill(VtxMass);
   }

   if (nVtx == 1 && !isCutEvt  && isCutVtx) {
      hData_NoEVTVtx_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_NoEVTVtx_1Vtx_CutEvt_Mass->Fill(VtxMass);
   }

   if (nVtx == 1 && isCutEvt && !isCutVtx) {
      hData_EVTNoVtx_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_EVTNoVtx_1Vtx_CutEvt_Mass->Fill(VtxMass);
   }

   if (nVtx == 1 && !isCutEvt && !isCutVtx) {
      hData_NoEVTNoVtx_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_NoEVTNoVtx_1Vtx_CutEvt_Mass->Fill(VtxMass);
   }
   //----------------------//

   if (nVtx == 2 && isCutEvt && isCutVtx) {
      hData_EVTVtx_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_EVTVtx_2Vtx_CutEvt_Mass->Fill(VtxMass);
      if (isCutVtx1 && isCutVtx2)
         {
            hData_EVTVtx_2VtxAll_CutEvt_Mmumu->Fill(tree_Mmumu);
            hData_EVTVtx_2VtxAll_CutEvt_Mass->Fill(VtxMass);
         }
   }

   if (nVtx == 2 && !isCutEvt && isCutVtx) {
      hData_NoEVTVtx_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_NoEVTVtx_2Vtx_CutEvt_Mass->Fill(VtxMass);
      if (isCutVtx1 && isCutVtx2)
         {
            hData_NoEVTVtx_2VtxAll_CutEvt_Mmumu->Fill(tree_Mmumu);
            hData_NoEVTVtx_2VtxAll_CutEvt_Mass->Fill(VtxMass);
         }
   }

   if (nVtx == 2 && isCutEvt && !isCutVtx ) {
      hData_EVTNoVtx_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_EVTNoVtx_2Vtx_CutEvt_Mass->Fill(VtxMass);
      if (!isCutVtx1 && !isCutVtx2)
         {
            hData_EVTNoVtx_2VtxAll_CutEvt_Mmumu->Fill(tree_Mmumu);
            hData_EVTNoVtx_2VtxAll_CutEvt_Mass->Fill(VtxMass);
         }
   }

   if (nVtx == 2 && !isCutEvt && !isCutVtx ) {
      hData_NoEVTNoVtx_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
      hData_NoEVTNoVtx_2Vtx_CutEvt_Mass->Fill(VtxMass);
      if (!isCutVtx1 && !isCutVtx2)
         {
            hData_NoEVTNoVtx_2VtxAll_CutEvt_Mmumu->Fill(tree_Mmumu);
            hData_NoEVTNoVtx_2VtxAll_CutEvt_Mass->Fill(VtxMass);
         }
   }

   //-----------------------------------------------------------//
   // ABCD using Hemipshere pt anf Tight+loose steps of vertexing 
   //-----------------------------------------------------------//
   //--------SR--------//
   
      if (!BlindSR) {
         
         if (Filter && (abs(tree_Hemi_eta->at(0)) > 2.4 || abs(tree_Hemi_eta->at(1)) > 2.4)) continue;
         hData_Hemi_BDTevt->Fill( tree_Evts_MVAval );
         //$$
         if (VtxMass > 8.) isCutVtx = true;
         if (hemi_ptmin > 80.) isCutEvt = true;
         //$$
         
         if (nVtx == 0 && isCutEvt) {
            hData_Hemi_0Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
            hData_Hemi_0Vtx_BDTevt->Fill( tree_Evts_MVAval );
         }

         
         if (nVtx == 1 && isCutEvt) {
            hData_Hemi_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
            hData_Hemi_1Vtx_CutEvt_Mass->Fill(VtxMass);
            hData_Hemi_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
            hData_Hemi_1Vtx_BDTvtx->Fill( BDTvtx );
         }

         
         if (nVtx == 2 && isCutEvt) {
            hData_Hemi_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
            hData_Hemi_2Vtx_CutEvt_Mass->Fill(VtxMass);
            hData_Hemi_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
            hData_Hemi_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx );

            hData_Hemi_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass0);
            hData_Hemi_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass1);
            hData_Hemi_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1 );
            hData_Hemi_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2 );
         }

      }

      //--------CR : Low Pt --------//

      isCutVtx = false; isCutVtx1 = false; isCutVtx2 = false;
      isCutEvt = false;
      BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;

      if (Filter && hemi_ptmin > 40. && hemi_ptmin < 80.)
         {
             //$$
            if ( abs(tree_Hemi_eta->at(0)) > 2.4 || abs(tree_Hemi_eta->at(1)) > 2.4 ) continue;
            hData_CRlowpt_BDTevt->Fill( tree_Evts_MVAval );
            if ( VtxMass > 8. ) isCutVtx = true; 
            // if ( hemi_ptmin > 80. ) isCutEvt = true;   
            //$$
            if (nVtx == 0 )
               {
                  hData_CRlowpt_0Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_CRlowpt_0Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  
               }
            if (nVtx == 1 )
               {
                  hData_CRlowpt_1Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_CRlowpt_1Vtx_CutEvt_Mass->Fill(   VtxMass );
                  hData_CRlowpt_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  hData_CRlowpt_1Vtx_BDTvtx->Fill( BDTvtx );
               }
            if (nVtx == 2 )
               {
                  hData_CRlowpt_2Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_CRlowpt_2Vtx_CutEvt_Mass->Fill(   VtxMass );
                  hData_CRlowpt_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
                  hData_CRlowpt_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx );

                  hData_CRlowpt_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass0 );
                  hData_CRlowpt_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass1 );
                  hData_CRlowpt_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1 );
                  hData_CRlowpt_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2 );
               }

         }
      //--------CR : Loose  --------//

      isCutVtx = false; isCutVtx1 = false; isCutVtx2 = false;
      isCutEvt = false;
      BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if ( Filter && abs(tree_Hemi_eta->at(0)) < 2.4 && abs(tree_Hemi_eta->at(1)) < 2.4 ) 
         {
            hData_CRloose_BDTevt->Fill( tree_Evts_MVAval );
            //$$
            if ( VtxMass > 8. ) isCutVtx = true; 
            if ( hemi_ptmin > 80. ) isCutEvt = true;
            //$$ 
            if (nVtxLoose == 0 && isCutEvt)
               {
                  hData_CRloose_0Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_CRloose_0Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  
               }
            if (nVtxLoose == 1 && isCutEvt)
               {
                  hData_CRloose_1Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_CRloose_1Vtx_CutEvt_Mass->Fill(   VtxMass );

                  hData_CRloose_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  hData_CRloose_1Vtx_BDTvtx->Fill( BDTvtx );
               }
            if (nVtxLoose == 2 && isCutEvt)
               {
                  hData_CRloose_2Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_CRloose_2Vtx_CutEvt_Mass->Fill(   VtxMass );
                  hData_CRloose_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
                  hData_CRloose_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx );

                  hData_CRloose_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass0 );
                  hData_CRloose_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass1 );
                  hData_CRloose_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1 );
                  hData_CRloose_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2 );
               }  
         }

      //--------CR : Loose Lowpt  --------//
      isCutVtx = false; isCutVtx1 = false; isCutVtx2 = false;
      isCutEvt = false;
      BDTvtx = -2.; BDTvtx1 = -2.; BDTvtx2 = -2.;
      if (Filter && hemi_ptmin > 40. && hemi_ptmin < 80. )
         {
            hData_CRlooselowpt_BDTevt->Fill( tree_Evts_MVAval );
            //$$
            if ( VtxMass > 8. ) isCutVtx = true; 
            // if ( hemi_ptmin > 80. ) isCutEvt = true;   
            //$$
            if (nVtxLoose == 0 )
               {
                  hData_CRlooselowpt_0Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_CRlooselowpt_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  hData_CRlooselowpt_1Vtx_BDTvtx->Fill( BDTvtx );
                  
               }
            if (nVtxLoose == 1 )
               {
                  hData_CRlooselowpt_1Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_CRlooselowpt_1Vtx_CutEvt_Mass->Fill(   VtxMass );
                  hData_CRlooselowpt_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
                  hData_CRlooselowpt_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx );
               }
            if (nVtxLoose == 2)
               {
                  hData_CRlooselowpt_2Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                  hData_CRlooselowpt_2Vtx_CutEvt_Mass->Fill(   VtxMass );

                  hData_CRlooselowpt_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass0 );
                  hData_CRlooselowpt_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass1 );
                  hData_CRlooselowpt_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1 );
                  hData_CRlooselowpt_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2 );
               }

         }
      //--------CR : fwd  --------//
      if(Filter && ((abs(tree_Hemi_eta->at(0)) > 2.4 && abs(tree_Hemi_eta->at(0)) < 3.0) ||
         (abs(tree_Hemi_eta->at(1)) > 2.4 && abs(tree_Hemi_eta->at(1)) < 3.0)))
            {
               if ( VtxMass > 8. ) isCutVtx = true; 
               if ( hemi_ptmin > 80. ) isCutEvt = true;   
               //$$
               if (nVtx == 0 && isCutEvt)
                  {
                     hData_CRfwd_0Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                     hData_CRfwd_0Vtx_BDTevt->Fill( tree_Evts_MVAval );

                  }
               if (nVtx == 1 && isCutEvt)
                  {
                     hData_CRfwd_1Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                     hData_CRfwd_1Vtx_CutEvt_Mass->Fill(   VtxMass );
                     hData_CRfwd_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
                     hData_CRfwd_1Vtx_BDTvtx->Fill( BDTvtx );
                  }
               if (nVtx == 2 && isCutEvt)
                  {
                     hData_CRfwd_2Vtx_CutEvt_Mmumu->Fill(  tree_Mmumu );
                     hData_CRfwd_2Vtx_CutEvt_Mass->Fill(   VtxMass );
                     hData_CRfwd_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
                     hData_CRfwd_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx );

                     hData_CRfwd_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass0 );
                     hData_CRfwd_2VtxAll_CutEvt_Mass->Fill( Vtx_Mass1 );
                     hData_CRfwd_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1 );
                     hData_CRfwd_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2 );
                  }
            }
         //--------CR : SameSign  --------//
      if (tree_FilterSameSign && !Filter &&
            abs(tree_Hemi_eta->at(0)) < 2.4 && abs(tree_Hemi_eta->at(1)) < 2.4)
         {
               hData_SameSign_BDTevt->Fill( tree_Evts_MVAval );
               if ( VtxMass > 8. ) isCutVtx = true; 
               if ( hemi_ptmin > 80. ) isCutEvt = true;   
               //$$
               if (nVtx == 0 && isCutEvt) {
                  hData_CRss_0Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
                  hData_CRss_0Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  
               }

               if (nVtx == 1 && isCutEvt) {
                  hData_CRss_1Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
                  hData_CRss_1Vtx_CutEvt_Mass->Fill(VtxMass);
                  hData_CRss_1Vtx_BDTevt->Fill( tree_Evts_MVAval );
                  hData_CRss_1Vtx_BDTvtx->Fill( BDTvtx );
               }

               
               if (nVtx == 2 && isCutEvt) {
                  hData_CRss_2Vtx_CutEvt_Mmumu->Fill(tree_Mmumu);
                  hData_CRss_2Vtx_CutEvt_Mass->Fill(VtxMass);
                  hData_CRss_2Vtx_CutEvt_BDTevt->Fill( tree_Evts_MVAval );
                  hData_CRss_2Vtx_CutEvt_MaxBDTvtx->Fill( BDTvtx );

                  hData_CRss_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass0);
                  hData_CRss_2VtxAll_CutEvt_Mass->Fill(Vtx_Mass1);
                  hData_CRss_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx1 );
                  hData_CRss_2VtxAll_CutEvt_BDTvtx->Fill( BDTvtx2 );
               }
         }

   } // end LOOP on events

   HistogramManager h ;
   h.WriteAllHistogramsInFile((Production+"/ABCD_"+thesample+".root").Data(),"recreate");

} // end ABCD::Loop
