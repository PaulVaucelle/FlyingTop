//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep  6 07:46:26 2024 by ROOT version 6.14/09
// from TTree ttree/summary information
// found on file: MiniSecInt_DoubleMuon_UL2018_MiniAODv2_GT36-v1.root
//////////////////////////////////////////////////////////

#ifndef TreeSecIntReader_h
#define TreeSecIntReader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
class TreeSecIntReader {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<int>     *minirunNumber;
   vector<int>     *minieventNumber;
   vector<int>     *minilumiBlock;

   vector<float>   *minitree_K0_reco_mass;
   vector<float>   *minitree_L0_reco_mass;


   // CMSSW collection
  vector<float>   *minitree_K0_x;
  vector<float>   *minitree_K0_y;
  vector<float>   *minitree_K0_z;
  vector<float>   *minitree_K0_r;
  vector<float>   *minitree_K0_NChi2;
  vector<float>   *minitree_K0_mass;
  vector<float>   *minitree_K0_eta;

  // CMSSW collection
  vector<float>   *minitree_L0_x;
  vector<float>   *minitree_L0_y;
  vector<float>   *minitree_L0_z;
  vector<float>   *minitree_L0_r;
  vector<float>   *minitree_L0_NChi2;
  vector<float>   *minitree_L0_mass;
  vector<float>   *minitree_L0_eta;

  // reco from us 
  vector<float>   *minitree_V0_reco_x;
  vector<float>   *minitree_V0_reco_y;
  vector<float>   *minitree_V0_reco_z;
  vector<float>   *minitree_V0_reco_r;
  vector<float>   *minitree_V0_reco_NChi2;
  vector<float>   *minitree_V0_reco_mass;
  vector<float>   *minitree_V0_reco_eta;
  vector<int>     *minitree_V0_reco_source;

   //CMSSW collection
  vector<float>   *minitree_Yc_x;
  vector<float>   *minitree_Yc_y;
  vector<float>   *minitree_Yc_z;
  vector<float>   *minitree_Yc_r;
  vector<int>     *minitree_Yc_layer;
  vector<float>   *minitree_Yc_NChi2;
  vector<float>   *minitree_Yc_eta;
  vector<float>   *minitree_Yc_mass;


   vector<int>     *minitree_tree_nPV;
   vector<int>     *minitree_nSecInt;
   vector<float>   *minitree_SecInt_x;
   vector<float>   *minitree_SecInt_y;
   vector<float>   *minitree_SecInt_z;
   vector<float>   *minitree_SecInt_r;
   vector<float>   *minitree_SecInt_d;
   vector<float>   *minitree_SecInt_drSig;
   vector<float>   *minitree_SecInt_dzSig;
   vector<float>   *minitree_SecInt_angleXY;
   vector<float>   *minitree_SecInt_angleZ;
   vector<float>   *minitree_SecInt_NChi2;
   vector<float>   *minitree_SecInt_ndf;
   vector<float>   *minitree_SecInt_mass;
   vector<float>   *minitree_SecInt_pt;
   vector<float>   *minitree_SecInt_eta;
   vector<float>   *minitree_SecInt_phi;
   vector<int>     *minitree_SecInt_charge;
   vector<bool>    *minitree_SecInt_badTkHit;
   vector<float>   *minitree_SecInt_dca;
   vector<bool>    *minitree_SecInt_selec;
   vector<int>     *minitree_SecInt_layer;

   vector<bool>    *minitree_V0_track_isFromV0;
   vector<bool>    *minitree_V0_track_isFromSI;
   vector<bool>    *minitree_V0_track_lost;
   vector<float>   *minitree_V0_track_pt;
   vector<float>   *minitree_V0_track_eta;
   vector<float>   *minitree_V0_track_phi;
   vector<int>     *minitree_V0_track_charge;
   vector<float>   *minitree_V0_track_NChi2;
   vector<float>   *minitree_V0_track_dxy;
   vector<float>   *minitree_V0_track_drSig;
   vector<float>   *minitree_V0_track_dz;
   vector<float>   *minitree_V0_track_dzSig;
   vector<int>     *minitree_V0_track_nHit;
   vector<int>     *minitree_V0_track_nHitPixel;
   vector<int>     *minitree_V0_track_firstHit;
   vector<float>   *minitree_V0_track_firstHit_x;
   vector<float>   *minitree_V0_track_firstHit_y;
   vector<float>   *minitree_V0_track_firstHit_z;
   vector<int>     *minitree_V0_track_iJet;
   vector<float>   *minitree_V0_track_ntrk10;
   vector<float>   *minitree_V0_track_ntrk20;
   vector<float>   *minitree_V0_track_ntrk30;
   vector<float>   *minitree_V0_track_ntrk40;
   vector<int>     *minitree_V0_track_Hemi;
   vector<float>   *minitree_V0_track_Hemi_dR;
   vector<float>   *minitree_V0_track_Hemi_dRmax;

   vector<int>       *minitree_smu_mass;
   vector<int>       *minitree_neu_mass;
   vector<float>     *minitree_neu_ctau; 

   // List of branches
   TBranch        *b_minirunNumber;   //!
   TBranch        *b_minieventNumber;   //!
   TBranch        *b_minilumiBlock;   //!
   TBranch        *b_minitree_K0_reco_mass;
   TBranch        *b_minitree_L0_reco_mass;

      // CMSSW collection
   TBranch   *b_minitree_K0_x;
   TBranch   *b_minitree_K0_y;
   TBranch   *b_minitree_K0_z;
   TBranch   *b_minitree_K0_r;
   TBranch   *b_minitree_K0_NChi2;
   TBranch   *b_minitree_K0_mass;
   TBranch   *b_minitree_K0_eta;

   // CMSSW collection
   TBranch   *b_minitree_L0_x;
   TBranch   *b_minitree_L0_y;
   TBranch   *b_minitree_L0_z;
   TBranch   *b_minitree_L0_r;
   TBranch   *b_minitree_L0_NChi2;
   TBranch   *b_minitree_L0_mass;
   TBranch   *b_minitree_L0_eta;

   // reco from us 
   TBranch   *b_minitree_V0_reco_x;
   TBranch   *b_minitree_V0_reco_y;
   TBranch   *b_minitree_V0_reco_z;
   TBranch   *b_minitree_V0_reco_r;
   TBranch   *b_minitree_V0_reco_NChi2;
   TBranch   *b_minitree_V0_reco_mass;
   TBranch   *b_minitree_V0_reco_eta;
   TBranch   *b_minitree_V0_reco_source;

      //CMSSW collection
   TBranch   *b_minitree_Yc_x;
   TBranch   *b_minitree_Yc_y;
   TBranch   *b_minitree_Yc_z;
   TBranch   *b_minitree_Yc_r;
   TBranch   *b_minitree_Yc_layer;
   TBranch   *b_minitree_Yc_NChi2;
   TBranch   *b_minitree_Yc_eta;
   TBranch   *b_minitree_Yc_mass;

   TBranch        *b_minitree_tree_nPV;
   TBranch        *b_minitree_nSecInt;   //!
   TBranch        *b_minitree_SecInt_x;   //!
   TBranch        *b_minitree_SecInt_y;   //!
   TBranch        *b_minitree_SecInt_z;   //!
   TBranch        *b_minitree_SecInt_r;   //!
   TBranch        *b_minitree_SecInt_d;   //!
   TBranch        *b_minitree_SecInt_drSig;   //!
   TBranch        *b_minitree_SecInt_dzSig;   //!
   TBranch        *b_minitree_SecInt_angleXY;   //!
   TBranch        *b_minitree_SecInt_angleZ;   //!
   TBranch        *b_minitree_SecInt_NChi2;   //!
   TBranch        *b_minitree_SecInt_ndf;   //!
   TBranch        *b_minitree_SecInt_mass;   //!
   TBranch        *b_minitree_SecInt_pt;   //!
   TBranch        *b_minitree_SecInt_eta;   //!
   TBranch        *b_minitree_SecInt_phi;   //!
   TBranch        *b_minitree_SecInt_charge;   //!
   TBranch        *b_minitree_SecInt_badTkHit;   //!
   TBranch        *b_minitree_SecInt_dca;   //!
   TBranch        *b_minitree_SecInt_selec;   //!
   TBranch        *b_minitree_SecInt_layer;   //!


   TBranch     *b_minitree_V0_track_isFromV0;
   TBranch     *b_minitree_V0_track_isFromSI;
   TBranch     *b_minitree_V0_track_lost;
   TBranch     *b_minitree_V0_track_pt;
   TBranch     *b_minitree_V0_track_eta;
   TBranch     *b_minitree_V0_track_phi;
   TBranch     *b_minitree_V0_track_charge;
   TBranch     *b_minitree_V0_track_NChi2;
   TBranch     *b_minitree_V0_track_dxy;
   TBranch     *b_minitree_V0_track_drSig;
   TBranch     *b_minitree_V0_track_dz;
   TBranch     *b_minitree_V0_track_dzSig;
   TBranch     *b_minitree_V0_track_nHit;
   TBranch     *b_minitree_V0_track_nHitPixel;
   TBranch     *b_minitree_V0_track_firstHit;
   TBranch     *b_minitree_V0_track_firstHit_x;
   TBranch     *b_minitree_V0_track_firstHit_y;
   TBranch     *b_minitree_V0_track_firstHit_z;
   TBranch     *b_minitree_V0_track_iJet;
   TBranch     *b_minitree_V0_track_ntrk10;
   TBranch     *b_minitree_V0_track_ntrk20;
   TBranch     *b_minitree_V0_track_ntrk30;
   TBranch     *b_minitree_V0_track_ntrk40;
   TBranch     *b_minitree_V0_track_Hemi;
   TBranch     *b_minitree_V0_track_Hemi_dR;
   TBranch     *b_minitree_V0_track_Hemi_dRmax;

   TBranch     *b_minitree_smu_mass;
   TBranch     *b_minitree_neu_mass;
   TBranch     *b_minitree_neu_ctau; 

   TreeSecIntReader(TTree *tree=0,TString Prod ="", TString sample="");
   virtual ~TreeSecIntReader();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString Prod,TString sample);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

      void initializeHisto(TString sample, bool isfirstset);
   void addHisto( TString var, TString selstep, TString sample, int nbins, float min, float max);
   void fillHisto(TString var, TString selstep, TString sample, float val, float weight);
   void addHisto2D(TString var, TString selstep, TString sample, int nxbins, float xmin, float xmax, int nybins, float ymin, float ymax);
   void fillHisto2D( TString var, TString selstep, TString sample, float xval,float yval, float weight);
   void GetHistDirectory();

   std::vector<TH1F*> histo_list_;
   std::vector<TH2F*> histo_list_2D_;
   std::map<std::string, int> histo_map_;
   std::map<std::string, int> histo_map_2D_;


   int numb_histo;
   int numb_histo_2D_;
   void deleteHisto();
};

#endif

#ifdef TreeSecIntReader_cxx
TreeSecIntReader::TreeSecIntReader(TTree *tree, TString Prod , TString sample) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+sample+".root";

   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(Path);
      if (!f || !f->IsOpen()) {
         f = new TFile(Path);
      }
      f->GetObject("ttree",tree);

   }
   Init(tree);
}

TreeSecIntReader::~TreeSecIntReader()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TreeSecIntReader::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TreeSecIntReader::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TreeSecIntReader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   minirunNumber = 0;
   minieventNumber = 0;
   minilumiBlock = 0;
   minitree_tree_nPV = 0;
   
   minitree_K0_reco_mass=0;
   minitree_L0_reco_mass=0;

   
   // CMSSW collection
     minitree_K0_x=0;
     minitree_K0_y=0;
     minitree_K0_z=0;
     minitree_K0_r=0;
     minitree_K0_NChi2=0;
     minitree_K0_mass=0;
     minitree_K0_eta=0;

  // CMSSW collection
     minitree_L0_x=0;
     minitree_L0_y=0;
     minitree_L0_z=0;
     minitree_L0_r=0;
     minitree_L0_NChi2=0;
     minitree_L0_mass=0;
     minitree_L0_eta=0;

  // reco from us 
     minitree_V0_reco_x=0;
     minitree_V0_reco_y=0;
     minitree_V0_reco_z=0;
     minitree_V0_reco_r=0;
     minitree_V0_reco_NChi2=0;
     minitree_V0_reco_mass=0;
     minitree_V0_reco_eta=0;
     minitree_V0_reco_source=0;

   //CMSSW collection
     minitree_Yc_x=0;
     minitree_Yc_y=0;
     minitree_Yc_z=0;
     minitree_Yc_r=0;
     minitree_Yc_layer=0;
     minitree_Yc_NChi2=0;
     minitree_Yc_eta=0;
     minitree_Yc_mass=0;


   minitree_nSecInt = 0;
   minitree_SecInt_x = 0;
   minitree_SecInt_y = 0;
   minitree_SecInt_z = 0;
   minitree_SecInt_r = 0;
   minitree_SecInt_d = 0;
   minitree_SecInt_drSig = 0;
   minitree_SecInt_dzSig = 0;
   minitree_SecInt_angleXY = 0;
   minitree_SecInt_angleZ = 0;
   minitree_SecInt_NChi2 = 0;
   minitree_SecInt_ndf = 0;
   minitree_SecInt_mass = 0;
   minitree_SecInt_pt = 0;
   minitree_SecInt_eta = 0;
   minitree_SecInt_phi = 0;
   minitree_SecInt_charge = 0;
   minitree_SecInt_badTkHit = 0;
   minitree_SecInt_dca = 0;
   minitree_SecInt_selec = 0;
   minitree_SecInt_layer = 0;

   minitree_V0_track_isFromV0 = 0;
   minitree_V0_track_isFromSI= 0;
   minitree_V0_track_lost= 0;
   minitree_V0_track_pt= 0;
   minitree_V0_track_eta= 0;
   minitree_V0_track_phi= 0;
   minitree_V0_track_charge= 0;
   minitree_V0_track_NChi2= 0;
   minitree_V0_track_dxy= 0;
   minitree_V0_track_drSig= 0;
   minitree_V0_track_dz= 0;
   minitree_V0_track_dzSig= 0;
   minitree_V0_track_nHit= 0;
   minitree_V0_track_nHitPixel= 0;
   minitree_V0_track_firstHit= 0;
   minitree_V0_track_firstHit_x= 0;
   minitree_V0_track_firstHit_y= 0;
   minitree_V0_track_firstHit_z= 0;
   minitree_V0_track_iJet= 0;
   minitree_V0_track_ntrk10= 0;
   minitree_V0_track_ntrk20= 0;
   minitree_V0_track_ntrk30= 0;
   minitree_V0_track_ntrk40= 0;
   minitree_V0_track_Hemi= 0;
   minitree_V0_track_Hemi_dR= 0;
   minitree_V0_track_Hemi_dRmax= 0;

   minitree_smu_mass=0;
   minitree_neu_mass=0;
   minitree_neu_ctau=0; 
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("minirunNumber", &minirunNumber, &b_minirunNumber);
   fChain->SetBranchAddress("minieventNumber", &minieventNumber, &b_minieventNumber);
   fChain->SetBranchAddress("minilumiBlock", &minilumiBlock, &b_minilumiBlock);
   fChain->SetBranchAddress("minitree_K0_reco_mass", &minitree_K0_reco_mass, &b_minitree_K0_reco_mass);
   fChain->SetBranchAddress("minitree_L0_reco_mass",&minitree_L0_reco_mass, &b_minitree_L0_reco_mass);

      
   // CMSSW collection
     fChain->SetBranchAddress("minitree_K0_x",&minitree_K0_x,&b_minitree_K0_x);
     fChain->SetBranchAddress("minitree_K0_y",&minitree_K0_y,&b_minitree_K0_y);
     fChain->SetBranchAddress("minitree_K0_z",&minitree_K0_z,&b_minitree_K0_z);
     fChain->SetBranchAddress("minitree_K0_r",&minitree_K0_r,&b_minitree_K0_r);
     fChain->SetBranchAddress("minitree_K0_NChi2",&minitree_K0_NChi2,&b_minitree_K0_NChi2);
     fChain->SetBranchAddress("minitree_K0_mass",&minitree_K0_mass,&b_minitree_K0_mass);
     fChain->SetBranchAddress("minitree_K0_eta",&minitree_K0_eta,&b_minitree_K0_eta);

  // CMSSW collection
     fChain->SetBranchAddress("minitree_L0_x",&minitree_L0_x,&b_minitree_L0_x);
     fChain->SetBranchAddress("minitree_L0_y",&minitree_L0_y,&b_minitree_L0_y);
     fChain->SetBranchAddress("minitree_L0_z",&minitree_L0_z,&b_minitree_L0_z);
     fChain->SetBranchAddress("minitree_L0_r",&minitree_L0_r,&b_minitree_L0_r);
     fChain->SetBranchAddress("minitree_L0_NChi2",&minitree_L0_NChi2,&b_minitree_L0_NChi2);
     fChain->SetBranchAddress("minitree_L0_mass",&minitree_L0_mass,&b_minitree_L0_mass);
     fChain->SetBranchAddress("minitree_L0_eta",&minitree_L0_eta,&b_minitree_L0_eta);

  // reco from us 
     fChain->SetBranchAddress("minitree_V0_reco_x",&minitree_V0_reco_x,&b_minitree_V0_reco_x);
     fChain->SetBranchAddress("minitree_V0_reco_y",&minitree_V0_reco_y,&b_minitree_V0_reco_y);
     fChain->SetBranchAddress("minitree_V0_reco_z",&minitree_V0_reco_z,&b_minitree_V0_reco_z);
     fChain->SetBranchAddress("minitree_V0_reco_r",&minitree_V0_reco_r,&b_minitree_V0_reco_r);
     fChain->SetBranchAddress("minitree_V0_reco_NChi2",&minitree_V0_reco_NChi2,&b_minitree_V0_reco_NChi2);
     fChain->SetBranchAddress("minitree_V0_reco_mass",&minitree_V0_reco_mass,&b_minitree_V0_reco_mass);
     fChain->SetBranchAddress("minitree_V0_reco_eta",&minitree_V0_reco_eta,&b_minitree_V0_reco_eta);
     fChain->SetBranchAddress("minitree_V0_reco_source",&minitree_V0_reco_source,&b_minitree_V0_reco_source);

   //CMSSW collection
     fChain->SetBranchAddress("minitree_Yc_x",&minitree_Yc_x,&b_minitree_Yc_x);
     fChain->SetBranchAddress("minitree_Yc_y",&minitree_Yc_y,&b_minitree_Yc_y);
     fChain->SetBranchAddress("minitree_Yc_z",&minitree_Yc_z,&b_minitree_Yc_z);
     fChain->SetBranchAddress("minitree_Yc_r",&minitree_Yc_r,&b_minitree_Yc_r);
     fChain->SetBranchAddress("minitree_Yc_layer",&minitree_Yc_layer,&b_minitree_Yc_layer);
     fChain->SetBranchAddress("minitree_Yc_NChi2",&minitree_Yc_NChi2,&b_minitree_Yc_NChi2);
     fChain->SetBranchAddress("minitree_Yc_eta",&minitree_Yc_eta,&b_minitree_Yc_eta);
     fChain->SetBranchAddress("minitree_Yc_mass",&minitree_Yc_mass,&b_minitree_Yc_mass);


   fChain->SetBranchAddress("minitree_tree_nPV",&minitree_tree_nPV, &b_minitree_tree_nPV);
   fChain->SetBranchAddress("minitree_nSecInt", &minitree_nSecInt, &b_minitree_nSecInt);
   fChain->SetBranchAddress("minitree_SecInt_x", &minitree_SecInt_x, &b_minitree_SecInt_x);
   fChain->SetBranchAddress("minitree_SecInt_y", &minitree_SecInt_y, &b_minitree_SecInt_y);
   fChain->SetBranchAddress("minitree_SecInt_z", &minitree_SecInt_z, &b_minitree_SecInt_z);
   fChain->SetBranchAddress("minitree_SecInt_r", &minitree_SecInt_r, &b_minitree_SecInt_r);
   fChain->SetBranchAddress("minitree_SecInt_d", &minitree_SecInt_d, &b_minitree_SecInt_d);
   fChain->SetBranchAddress("minitree_SecInt_drSig", &minitree_SecInt_drSig, &b_minitree_SecInt_drSig);
   fChain->SetBranchAddress("minitree_SecInt_dzSig", &minitree_SecInt_dzSig, &b_minitree_SecInt_dzSig);
   fChain->SetBranchAddress("minitree_SecInt_angleXY", &minitree_SecInt_angleXY, &b_minitree_SecInt_angleXY);
   fChain->SetBranchAddress("minitree_SecInt_angleZ", &minitree_SecInt_angleZ, &b_minitree_SecInt_angleZ);
   fChain->SetBranchAddress("minitree_SecInt_NChi2", &minitree_SecInt_NChi2, &b_minitree_SecInt_NChi2);
   fChain->SetBranchAddress("minitree_SecInt_ndf", &minitree_SecInt_ndf, &b_minitree_SecInt_ndf);
   fChain->SetBranchAddress("minitree_SecInt_mass", &minitree_SecInt_mass, &b_minitree_SecInt_mass);
   fChain->SetBranchAddress("minitree_SecInt_pt", &minitree_SecInt_pt, &b_minitree_SecInt_pt);
   fChain->SetBranchAddress("minitree_SecInt_eta", &minitree_SecInt_eta, &b_minitree_SecInt_eta);
   fChain->SetBranchAddress("minitree_SecInt_phi", &minitree_SecInt_phi, &b_minitree_SecInt_phi);
   fChain->SetBranchAddress("minitree_SecInt_charge", &minitree_SecInt_charge, &b_minitree_SecInt_charge);
   fChain->SetBranchAddress("minitree_SecInt_badTkHit", &minitree_SecInt_badTkHit, &b_minitree_SecInt_badTkHit);
   fChain->SetBranchAddress("minitree_SecInt_dca", &minitree_SecInt_dca, &b_minitree_SecInt_dca);
   fChain->SetBranchAddress("minitree_SecInt_selec", &minitree_SecInt_selec, &b_minitree_SecInt_selec);
   fChain->SetBranchAddress("minitree_SecInt_layer", &minitree_SecInt_layer, &b_minitree_SecInt_layer);

   fChain->SetBranchAddress("minitree_V0_track_isFromV0",&minitree_V0_track_isFromV0,&b_minitree_V0_track_isFromV0);
   fChain->SetBranchAddress("minitree_V0_track_isFromSI",&minitree_V0_track_isFromSI,&b_minitree_V0_track_isFromSI);
   fChain->SetBranchAddress("minitree_V0_track_lost",&minitree_V0_track_lost,&b_minitree_V0_track_lost);
   fChain->SetBranchAddress("minitree_V0_track_pt",&minitree_V0_track_pt,&b_minitree_V0_track_pt);
   fChain->SetBranchAddress("minitree_V0_track_eta",&minitree_V0_track_eta,&b_minitree_V0_track_eta);
   fChain->SetBranchAddress("minitree_V0_track_phi",&minitree_V0_track_phi,&b_minitree_V0_track_phi);
   fChain->SetBranchAddress("minitree_V0_track_charge",&minitree_V0_track_charge,&b_minitree_V0_track_charge);
   fChain->SetBranchAddress("minitree_V0_track_NChi2",&minitree_V0_track_NChi2,&b_minitree_V0_track_NChi2);
   fChain->SetBranchAddress("minitree_V0_track_dxy",&minitree_V0_track_dxy,&b_minitree_V0_track_dxy);
   fChain->SetBranchAddress("minitree_V0_track_drSig",&minitree_V0_track_drSig,&b_minitree_V0_track_drSig);
   fChain->SetBranchAddress("minitree_V0_track_dz",&minitree_V0_track_dz,&b_minitree_V0_track_dz);
   fChain->SetBranchAddress("minitree_V0_track_dzSig",&minitree_V0_track_dzSig,&b_minitree_V0_track_dzSig);
   fChain->SetBranchAddress("minitree_V0_track_nHit",&minitree_V0_track_nHit,&b_minitree_V0_track_nHit);
   fChain->SetBranchAddress("minitree_V0_track_nHitPixel",&minitree_V0_track_nHitPixel,&b_minitree_V0_track_nHitPixel);
   fChain->SetBranchAddress("minitree_V0_track_firstHit",&minitree_V0_track_firstHit,&b_minitree_V0_track_firstHit);
   fChain->SetBranchAddress("minitree_V0_track_firstHit_x",&minitree_V0_track_firstHit_x,&b_minitree_V0_track_firstHit_x);
   fChain->SetBranchAddress("minitree_V0_track_firstHit_y",&minitree_V0_track_firstHit_y,&b_minitree_V0_track_firstHit_y);
   fChain->SetBranchAddress("minitree_V0_track_firstHit_z",&minitree_V0_track_firstHit_z,&b_minitree_V0_track_firstHit_z);
   fChain->SetBranchAddress("minitree_V0_track_iJet",&minitree_V0_track_iJet,&b_minitree_V0_track_iJet);
   fChain->SetBranchAddress("minitree_V0_track_ntrk10",&minitree_V0_track_ntrk10,&b_minitree_V0_track_ntrk10);
   fChain->SetBranchAddress("minitree_V0_track_ntrk20",&minitree_V0_track_ntrk20,&b_minitree_V0_track_ntrk20);
   fChain->SetBranchAddress("minitree_V0_track_ntrk30",&minitree_V0_track_ntrk30,&b_minitree_V0_track_ntrk30);
   fChain->SetBranchAddress("minitree_V0_track_ntrk40",&minitree_V0_track_ntrk40,&b_minitree_V0_track_ntrk40);
   fChain->SetBranchAddress("minitree_V0_track_Hemi",&minitree_V0_track_Hemi,&b_minitree_V0_track_Hemi);
   fChain->SetBranchAddress("minitree_V0_track_Hemi_dR",&minitree_V0_track_Hemi_dR,&b_minitree_V0_track_Hemi_dR);
   fChain->SetBranchAddress("minitree_V0_track_Hemi_dRmax",&minitree_V0_track_Hemi_dRmax,&b_minitree_V0_track_Hemi_dRmax);

   fChain->SetBranchAddress("minitree_smu_mass",&minitree_smu_mass,&b_minitree_smu_mass);
   fChain->SetBranchAddress("minitree_neu_mass",&minitree_neu_mass,&b_minitree_neu_mass);
   fChain->SetBranchAddress("minitree_neu_ctau",&minitree_neu_ctau,&b_minitree_neu_ctau);

   Notify();
}

Bool_t TreeSecIntReader::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TreeSecIntReader::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TreeSecIntReader::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//$$
// double DeltaR(double eta1, double phi1, double eta2, double phi2) {
//   double DeltaPhi = TMath::Abs(phi2 - phi1);
//   if (DeltaPhi > 3.141593 ) DeltaPhi = 2.*3.141593 - DeltaPhi;
//   return TMath::Sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
// }
// double DeltaPhi(double phi1, double phi2) {
//   double DeltaPhi = phi1 - phi2;
//   if (abs(DeltaPhi) > 3.141593 ) {
//     DeltaPhi = 2.*3.141593 - abs(DeltaPhi);
//     DeltaPhi = -DeltaPhi * (phi1 - phi2) / abs(phi1 - phi2);
//   }
//   return DeltaPhi;
// }
#endif // #ifdef TreeSecIntReader_cxx
