#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>

void computeTH2FRatio(TString Year) {
    // Ouvrir le premier fichier et récupérer l'histogramme

    int stati=0;
    bool fit= 1;
    bool logy=1;

    TString Yearcor = Year;
    if (Year == "2016PRE") Yearcor = "2016preVFP";
    if (Year == "2016POST") Yearcor = "2016";
    //MUMU
    TString extension = "";
    TFile* f1_Data_emu  = new TFile("../../DATA_EMU_"+Year+"_31_10_2024/DATAMC_MuonEG-UL2018_MiniAODv2_GT36-v1"+extension+".root");


    TFile* f1_DY  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8"+extension+".root");
    TFile* f2_DY  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"+extension+".root");
    TFile* f1_TT  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root");
    TFile* f2_TT  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8"+extension+".root");
    TFile* f1_ST  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"+extension+".root");
    TFile* f2_ST  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8"+extension+".root");
    TFile* f1_TTV = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8"+extension+".root");
    TFile* f2_TTV = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8"+extension+".root");
    TFile* f3_TTV = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_TTWW_TuneCP5_13TeV-madgraph-pythia8"+extension+".root");
    TFile* f1_VV  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8"+extension+".root");
    TFile* f2_VV  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8"+extension+".root");
    TFile* f3_VV  = new TFile("../../MC_EMU_"+Year+"_31_10_2024_v2/DATAMC_ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8"+extension+".root");
    
    TString DATAFILE[1] = {"MuonEG-UL2018_MiniAODv2_GT36-v1_"
    };
    if (Year == "2018") DATAFILE[0] = "MuonEG-UL2018_MiniAODv2_GT36-v1_";


    TString MCFILE[12] = {
                    "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_",
                    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_",
                    "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_",
                    "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_",
                    "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_",
                    "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8_",
                    "ttWJetsToLNu_5f_EWK_TuneCP5_13TeV_amcatnlo-pythia8_",
                    "TTZToLL_5f_TuneCP5_13TeV-madgraphMLM-pythia8_",
                    "TTWW_TuneCP5_13TeV-madgraph-pythia8_",
                    "WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_",
                    "WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_",
                    "ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8_"
    };

    TString ytitle = "#eta";
    TString xtitle = "pt"; 
    int nbinx = 40; 
    float xmin = 0;
    float xmax =  200;
    int nbiny = 25; 
    float ymin = -2.4;
    float ymax =  2.4;
    TString HeaderCMS = "CMS";
    TString htitleA = "leading_lepton_pt_eta_reco";
    TString HeaderA = "A";

    if (Year == "2016") HeaderCMS = "2016                            36.3 fb^{-1} (13 TeV)";
    if (Year == "2017") HeaderCMS = "2017                            41.5 fb^{-1} (13 TeV)";
    if (Year == "2018") HeaderCMS = "2018                            59.8 fb^{-1} (13 TeV)";


    TH2F* g1_DY = (TH2F*)gROOT->FindObject(MCFILE[0]+htitleA);//ok
    TH2F* g2_DY = (TH2F*)gROOT->FindObject(MCFILE[1]+htitleA);//ok
    TH2F*  h_DY = new TH2F("h_DY","",nbinx,xmin,xmax, nbiny,ymin,ymax);

    
    TH2F* g1_VV = (TH2F*)gROOT->FindObject(MCFILE[9]+htitleA);//ok
    TH2F* g2_VV = (TH2F*)gROOT->FindObject(MCFILE[10]+htitleA);//ok
    TH2F* g3_VV = (TH2F*)gROOT->FindObject(MCFILE[11]+htitleA);//ok
    TH2F*  h_VV = new TH2F("h_VV","",nbinx,xmin,xmax, nbiny,ymin,ymax);
    
    TH2F* g1_ST = (TH2F*)gROOT->FindObject(MCFILE[4]+htitleA);//ok
    TH2F* g2_ST = (TH2F*)gROOT->FindObject(MCFILE[5]+htitleA);//ok
    TH2F*  h_ST = new TH2F("h_ST","",nbinx,xmin,xmax, nbiny,ymin,ymax);

    TH2F* g1_TT = (TH2F*)gROOT->FindObject(MCFILE[2]+htitleA);//ok
    TH2F* g2_TT = (TH2F*)gROOT->FindObject(MCFILE[3]+htitleA);//ok
    TH2F*  h_TT = new TH2F("h_TT","",nbinx,xmin,xmax, nbiny,ymin,ymax);

    TH2F* g1_TTV = (TH2F*)gROOT->FindObject(MCFILE[6]+htitleA);//ok
    TH2F* g2_TTV = (TH2F*)gROOT->FindObject(MCFILE[7]+htitleA);//ok 
    TH2F* g3_TTV = (TH2F*)gROOT->FindObject(MCFILE[8]+htitleA);//ok
    TH2F*  h_TTV = new TH2F("h_TTV","",nbinx,xmin,xmax, nbiny,ymin,ymax);

    TH2F* htotMC  = new TH2F("htotMC","",nbinx,xmin,xmax, nbiny,ymin,ymax);
    TH2F* htotData  = new TH2F("htotData","",nbinx,xmin,xmax, nbiny,ymin,ymax);

    htotData->Sumw2();
    htotMC->Sumw2();

 
    f1_Data_emu->cd();
    TH2F* g1_Data_emu = (TH2F*)gROOT->FindObject(DATAFILE[0]+htitleA);
    htotData->Add(g1_Data_emu, htotData, 1, 0);


 f1_DY->cd();
 g1_DY = (TH2F*)gROOT->FindObject(MCFILE[0]+htitleA);
 g1_DY->Sumw2();
 h_DY = new TH2F("h_DY","",nbinx,xmin,xmax, nbiny,ymin,ymax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 g2_DY = (TH2F*)gROOT->FindObject(MCFILE[1]+htitleA);
 h_DY->Add(g2_DY, h_DY, 1,1);

 f1_VV->cd();
 g1_VV = (TH2F*)gROOT->FindObject(MCFILE[9]+htitleA);
 h_VV = new TH2F("h_VV","",nbinx,xmin,xmax, nbiny,ymin,ymax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH2F*)gROOT->FindObject(MCFILE[10]+htitleA);
 h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 g3_VV = (TH2F*)gROOT->FindObject(MCFILE[11]+htitleA);
 h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH2F*)gROOT->FindObject(MCFILE[6]+htitleA);
 h_TTV = new TH2F("h_TTV","",nbinx,xmin,xmax, nbiny,ymin,ymax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 g2_TTV = (TH2F*)gROOT->FindObject(MCFILE[7]+htitleA);
 h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH2F*)gROOT->FindObject(MCFILE[8]+htitleA);
 h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH2F("h_ST","",nbinx,xmin,xmax, nbiny,ymin,ymax);
 f1_ST->cd();
 g1_ST = (TH2F*)gROOT->FindObject(MCFILE[4]+htitleA);
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 g2_ST = (TH2F*)gROOT->FindObject(MCFILE[5]+htitleA);
 h_ST->Add(g2_ST, h_ST, 1, 1);

 f1_TT->cd();
g1_TT = (TH2F*)gROOT->FindObject(MCFILE[2]+htitleA);
 h_TT = new TH2F("h_TT","",nbinx,xmin,xmax, nbiny,ymin,ymax);
 h_TT->Add(g1_TT, h_TT, 0.5,0);

f2_TT->cd();
 g2_TT = (TH2F*)gROOT->FindObject(MCFILE[3]+htitleA);
 h_TT->Add(g2_TT, h_TT, 1,1);


 htotMC->Add(htotMC, h_DY, 1, 1);
 htotMC->Add(htotMC, h_VV, 1, 1);
 htotMC->Add(htotMC, h_TTV, 1, 1);
 htotMC->Add(htotMC, h_ST, 1, 1);
 htotMC->Add(htotMC, h_TT, 1, 1);

 h_VV->Add(h_VV, h_TTV, 1, 1);
 h_VV->Add(h_VV, h_ST, 1, 1);
 h_VV->Add(h_VV, h_TT, 1, 1);

 h_TTV->Add(h_TTV, h_ST, 1, 1);
 h_TTV->Add(h_TTV, h_TT, 1, 1);

 h_ST->Add(h_ST, h_TT, 1, 1);



    // Calcul du ratio htotData / h2_sum
    TH2F *h_ratio = (TH2F*)htotData->Clone("h_ratio");
    TH2F *h_err = (TH2F*)htotData->Clone("h_err");
    h_ratio->Reset();
    h_err->Reset();
    
    for (int x = 1; x <= htotData->GetNbinsX(); x++) {
        for (int y = 1; y <= htotData->GetNbinsY(); y++) {
            double num = htotData->GetBinContent(x, y);
            double den = htotMC->GetBinContent(x, y);
            double num_err = htotData->GetBinError(x, y);
            double den_err = htotMC->GetBinError(x, y);
            
            if (den != 0) {
                double ratio = num / den;
                double ratio_err = ratio * sqrt((num_err/num)*(num_err/num) + (den_err/den)*(den_err/den));
                h_ratio->SetBinContent(x, y, ratio);
                h_ratio->SetBinError(x, y, ratio_err);
                h_err->SetBinContent(x, y, ratio_err);
                h_err->SetBinError(x, y, 0);
            } else {
                h_ratio->SetBinContent(x, y, 0);
                h_ratio->SetBinError(x, y, 0);
                h_err->SetBinContent(x, y, 0);
                h_err->SetBinError(x, y, 0);
            }
        }
    }
    
    // !!
        // Style pour afficher la palette et les valeurs numÃ©riques
    gStyle->SetOptStat(0);   // DÃ©sactive la boÃ®te statistique
    gStyle->SetPaintTextFormat("4.0f"); // Format des valeurs numÃ©riques (ici entier)
    gStyle->SetPalette(kViridis);
    
    
    // Dessin de l'histogramme avec COLZ et texte
    h_err->Draw("COLZ TEXT");
    h_ratio->Draw("COLZ TEXT");
    // !!

    // Sauvegarder le résultat
    h_ratio->SaveAs("2DSFele.root");
    h_err->SaveAs("2DSFele_err.root");
    std::cout << "Ratio TH2F sauvegarde" << std::endl;
}
