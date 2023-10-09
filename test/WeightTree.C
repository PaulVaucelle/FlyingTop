#include <iostream>
#include <TROOT.h>
#include "TVectorD.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TTree.h"

void WeightTree(TString thesample="")//std::vector<string> SampleList or string plot
{

    // if (SampleList.size()==0)
    //     {
    //         std::cout<<"Please, select samples to analyse and make sure they are all up to date"
    //     }
//  TFile *f1, *f2, *f3, *f4, *f5, *f6, *f7, *f8, *f9, *f10, *f11, *f12, *f13, *f14, *f15, *f16, *f17, *f18, *f19, *f20, *f21 , *f22 , *f23 , *f24 , *f25 , *f26, *f27, *f28, *f29, *f30 , *f31;

  //--------------------------------------------------------------------------------//
 float NormFactor = 1. ;
  float XS = 1;
if (thesample == "70_test"){XS = 0.0035;}
if (thesample == "50_test")  {  XS = 0.0025;}
if (thesample == "30_test") {  XS = 0.002; }
if (thesample == "10_test"){ XS = 0.0035;}
if (thesample == "DYJetsToLL_M10to50") {XS = 15910.0;}
if (thesample == "ST_tW_antitop_5f_NoFullyHadronicDecays"){  XS = 32.51;}
if (thesample == "ST_tW_top_5f_NoFullyHadronicDecays"){ XS = 32.45;}
if (thesample == "TTJets_DiLept"){  XS = 53.07;}
if (thesample == "WWTo2L2Nu_MLL_200To600"){ XS = 11.09;}
if (thesample == "WWTo2L2Nu_MLL_600To1200"){ XS = 11.09;}
if (thesample == "WWTo2L2Nu"){ XS = 11.09;}
if (thesample == "WZTo2Q2L_mllmin4p0"){ XS = 6.535;}
if (thesample == "ZZTo2Q2L_mllmin4p0"){  XS = 3.676;}
if (thesample == "TTTo2L2Nu"){XS = 88.3;}
if (thesample == "DYJetsToLL_M50"){ XS = 5379;}
if (thesample == "ttWJetsToLNu_5f_EWK"){        
    XS = 0.290;// not found on XSDB, no file on tier2...approximation
      //Took 0.868 pb (CMS-TOP-21-011)
     // as a starting point and then divided by 3 (lepton universality)
    }
if (thesample == "TTZToLL_5f"){XS = 0.05188;//Not found on XSDB => used ana.py macro : 5.188e-02 +- 2.437e-04 pb
}
if (thesample == "TTToHadronic"){XS = 687.1;}
if (thesample == "TTWW")
    {        
      XS = 0.006992;//found on XSDB
    }
if (thesample.Contains("smu200"))
    {        
      XS = 0.01;
    }
    if (thesample.Contains("smu250"))
    {        
      XS = 0.0045;
    }
    if (thesample.Contains("smu300"))
    {        
      XS = 0.002;
    }

    if (thesample.Contains("smu400"))
    {        
      XS = 0.0006;
    }

    if (thesample.Contains("smu500"))
    {        
      XS = 0.00025;
    }
    std::cout<<"------------Process Starting------------"<<std::endl;
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/Ntuple_"+thesample+".root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/Ntuple_"+thesample+".root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_20_FLY/src/FlyingTop/FlyingTop/test/Ntuple_"+thesample+".root:/FlyingTop");
      TString outfileName("Ntuple_"+thesample+"_weighted.root");
      TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
      TTree *TreeNW ;
      dir->GetObject("ttree",TreeNW);
      TTree *TreeW = TreeNW->CloneTree();
      float nentries = 1;
      nentries = TreeW->GetBranch("eventNumber")->GetEntries();
      NormFactor = XS/nentries;
TreeW->SetWeight(NormFactor);
std::cout<<" Processed sample : "<<thesample<<" with XS : "<<XS<<"pb and norm. factor (XS/Nevts): "<<NormFactor<<std::endl;
outputFile->cd();
TreeW->Write();
f->Close();
outputFile->Close();



}
