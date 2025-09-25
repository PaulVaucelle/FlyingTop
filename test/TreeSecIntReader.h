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

   // List of branches
   TBranch        *b_minirunNumber;   //!
   TBranch        *b_minieventNumber;   //!
   TBranch        *b_minilumiBlock;   //!
   TBranch        *b_minitree_K0_reco_mass;
   TBranch        *b_minitree_L0_reco_mass;
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
TString Path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_8_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/MiniSecInt_"+sample+".root";

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
