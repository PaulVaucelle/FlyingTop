#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

void plot_VtxStepEff() {
    // Open the root file
    TString Prod = "PROD_CSI_10_06_2024";
    TString Sample = "RPV_2018_ctauALL";
    TString filename = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/"+Prod+"/Mini"+Sample+".root";
    TFile* file = TFile::Open(filename);
    if (!file) {
        printf("Error opening root file: %s\n", filename);
        return;
    }
// MiniDYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.root  
    // Retrieve the tree from the file
    TTree* tree = dynamic_cast<TTree*>(file->Get("ttree"));
    if (!tree) {
        printf("Error retrieving tree from root file: %s\n", filename);
        file->Close();
        return;
    }

    // Retrieve the vector from the tree
    std::vector<float>* minitree_Hemi_Vtx_step = NULL;
    tree->SetBranchAddress("minitree_Hemi_Vtx_step", &minitree_Hemi_Vtx_step);

    // Create a histogram
    TH1F* histogram = new TH1F("VtxStepEffi", "VTx step Efficiency", 4, 1, 5);

    // Loop over the entries in the tree
    Long64_t numEntries = tree->GetEntries();
    for (Long64_t i = 0; i < numEntries; ++i) {
        tree->GetEntry(i);

        // Fill the histogram with the values from the vector
        for (unsigned int j = 0; j < minitree_Hemi_Vtx_step->size(); ++j) {
            histogram->Fill(minitree_Hemi_Vtx_step->at(j));
        }
    }
    histogram->Sumw2(); 
    if (histogram->Integral(1,5)>0)
        {
            histogram->Scale(1/histogram->Integral(1,5));
        
        }
    else    
        {
            printf("Error: histogram->Integral(0,4) = 0\n");
        }


    // Save the histogram to a root file
    TFile* outputFile = TFile::Open("./VtxStepEfficiency_"+Sample+".root", "RECREATE");
    if (outputFile) {
        histogram->Write();
        outputFile->Close();
    } else {
        printf("Error creating output root file\n");
    }
    // Clean up
    file->Close();
    // delete histogram;
}

