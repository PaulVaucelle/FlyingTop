#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>    
#include <TLatex.h>
#include <TNamed.h>
#include <TString.h>
#include <TMath.h>


void ToysGen(TString Var, TString Channel, TString Year) {
 
  // HTL_VtxBDT_Ave_Corr 
  // HTL_EventBDT_Corr
  // HTL_STW_6Bins_Corr 
  // !! Parameters to change
    // !! §§ ---------------------
    // TString Var = "HTL_EventBDT_Corr"; // SI VTX ou EVT BDT borne : -1,1, SI STW borne : 1,8 
    // TString Channel = "MUMU"; // Channel name : "EMU", "MUMU"
    // TString Year= "2018"; // Year : "2016PRE", "2016POST", "2017", "2018"
    // !! §§ ---------------------


  TString outputfileName = "hsolveToys";
  TString outHistoName = "hToys";
 long nToys = 10000000;
  bool SpecificFit = false; // Si true, utilise un fit spécifique, sinon fit exponentiel décroissant par défaut
  float NormIntegral = 1.0;
  int nbin = 10; // Nombre de bins de l'histogramme à générer
  float xmin = -1.0;
  float xmax = 1.0;
  float BinEdgesSTW[7] = {1,2,3,4,5,6,7};
  float BinEdgesAveBDT[11] = {-1,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1};
  float binEdgesEvtBDT[11] = {-1,-0.8,-0.6,-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1};
  float* binEdges = nullptr;
  float par0 = 1.0;
  float par1 = 1.0;
  float par2 = 1.0;
  float par3 = 1.0;
  float par4 = 1.0;

  if (Var == "HTL_STW_6Bins_Corr") // Pas besoin
    {
        xmin = 1.0;
        xmax = 7;
        nbin = 6; // Nombre de bins pour STW
        binEdges = BinEdgesSTW; // Tableau des bords de bins pour STW
        outputfileName = "hsolveToys_"+Var+"_Data_"+Channel+"_"+Year; // Nom du fichier de sortie pour STW
        SpecificFit = false;

        // !! MuMu
        if (Year == "2018")
            {
                par0 = 5.67775e+02;
                par1 = 1.18707e+00;
                NormIntegral = 303.753;
                nToys = 1000000; // Nombre de toys pour STW
            }
        else if (Year == "2017")
            {
                par0 = 4.29017e+02;
                par1 = 1.03548e+00;
                NormIntegral = 163.616;  
                nToys = 1000000; // Nombre de toys pour STW  
            }
        else if (Year == "2016POST")
            {
                par0 = 2.13502e+02;
                par1 = 9.50614e-01;
                NormIntegral = 49.4446; 
                nToys = 1000000; // Nombre de toys pour STW   
            }
        else if (Year == "2016PRE")
            {
                par0 = 1.84999e+02;
                par1 = 9.15836e-01;
                NormIntegral = 52.6894; 
                nToys = 1000000; // Nombre de toys pour STW   
            }

    }
else if (Var == "HTL_VtxBDT_Ave_Corr")
    {
        xmin = -1.0;
        xmax = 1.0;
        nbin = 10; // Nombre de bins pour BDT
        binEdges = BinEdgesAveBDT; // Tableau des bords de bins pour VTX BDT
        outputfileName = "hsolveToys_"+Var+"_Data_"+Channel+"_"+Year;
        SpecificFit = true;
        // !! MuMu
        if (Year == "2018")
            {
                par0 = 5.36957e+00;
                par1 = 2.75681e-01;
                par2 = 1.00000e-03;
                par3 = 1.05850e-01;
                par4 = 2.80186e+01;
                NormIntegral = 319.593;
                nToys = 1000000; // Nombre de toys pour STW
            }
        else if (Year == "2017")
            {
                par0 = 2.70202e-01;
                par1 = 1.53807e-01;
                par2 = 1.00001e-02;
                par3 = 3.21315e-02;
                par4 = 9.71852e+00;
                NormIntegral = 164.066;  
                nToys = 1000000; // Nombre de toys pour STW  
            }
        else if (Year == "2016POST")
            {
                par0 = 9.27828e-01;
                par1 = 2.70480e-01;
                par2 = 1.00000e-02;
                par3 = -7.19460e-02;
                par4 = 3.82646e+00;
                NormIntegral = 52.8422; 
                nToys = 1000000; // Nombre de toys pour STW   
            }
        else if (Year == "2016PRE")
            {
                par0 = 8.19949e-01;
                par1 = 3.44227e-01;
                par2 = 1.00000e-02;
                par3 = -2.36741e-01;
                par4 = 4.04205e+01;
                NormIntegral = 80.32; 
                nToys = 1000000; // Nombre de toys pour STW   
            }
    }
else if (Var == "HTL_EventBDT_Corr")
    {
        xmin = -1.0;
        xmax = 1.0;
        nbin = 10; // Nombre de bins pour EVT BDT
        binEdges = binEdgesEvtBDT; // Tableau des bords de bins pour EVT BDT
        outputfileName = "hsolveToys_"+Var+"_Data_"+Channel+"_"+Year;
        // !! MuMu
        if (Year == "2018")
            {
                par0 = 1.00000e-03;
                par1 = 7.24093e-02;
                par2 = 6.69426e-03;
                par3 = -3.22336e-01;
                par4 = 2.54671e+01;

                NormIntegral = 315.52;

                nToys = 1000000; // Nombre de toys pour STW

                SpecificFit = true;
            }
        else if (Year == "2017")
            {
                par0 = 5.84920e-01;
                par1 = 1.61849e-01;
                NormIntegral = 158.123;

                nToys = 1000000; // Nombre de toys pour STW 

                SpecificFit = false;
            }
        else if (Year == "2016POST")
            {
                par0 = 4.86024e+00;
                par1 = 4.03230e-01;
                NormIntegral = 55.5108; 

                nToys = 1000000; // Nombre de toys pour STW  

                SpecificFit = false; 
            }
        else if (Year == "2016PRE")
            {
                par0 = 1.00001e-03;
                par1 = 8.84163e-02;
                par2 = 5.00000e-02;
                par3 = -3.00000e-01;
                par4 = 1.55290e+01;

                NormIntegral = 64.7945; 

                nToys = 1000000; // Nombre de toys pour STW 

                SpecificFit = true;  
            }
    }

  // !! §§ ---------------------

// !! -------------- !!
// !! -------------- !!
// !! -------------- !!
    // Définir la fonction (PDF)
    TF1* pdf = new TF1("pdf", "[0]*exp(-x/[1])", xmin,xmax); // ici une exponentielle entre 0 et 10
    pdf->SetParameter(0, par0); // lambda < 0 pour décroissance
    pdf->SetParameter(1, par1); //tau 
    if (SpecificFit)
        {
            pdf = new TF1("pdf", "[0]*exp(-x/[1])+ [4]*[2]/([2]+(x-[3])**2)", xmin, xmax); // Utiliser la fonction crystalBall
            pdf->SetParameter(0, par0); // lambda < 0 pour décroissance
            pdf->SetParameter(1, par1); //tau 
            pdf->SetParameter(2, par2); // LorentzNorm
            pdf->SetParameter(3, par3); //LorentzShift
            pdf->SetParameter(4, par4); // Scale
        }

// !! -------------- !!
// !! -------------- !!
// !! -------------- !!


    // Important : s'assurer que la fonction est positive
    if (pdf->Eval(pdf->GetXmin()) <= 0 || pdf->Eval(pdf->GetXmax()) <= 0) {
        std::cout << "Attention : la fonction n'est pas strictement positive sur le domaine !" << std::endl;
        return;
    }



    // Créer un histogramme pour les jouets
    TH1F* hToys = new TH1F("hToys", "Toys from TF1 PDF", nbin, binEdges); // Utiliser les bords de bins définis
    hToys->SetName(outHistoName); // Nom de l'histogramme
    

    // Générer les jouets
    for (int i = 0; i < nToys; ++i) {
        double x = pdf->GetRandom(); // tirage selon la forme de la fonction
        hToys->Fill(x);
    }

    hToys->Scale(1.0 / hToys->Integral()); // Normaliser l'histogramme à l'intégrale de 1
    hToys->Scale(NormIntegral); // Normaliser à l'intégrale souhaitée

    TH1F* hToys_NormUp = (TH1F*)hToys->Clone(); // Utiliser les bords de bins définis
    TH1F* hToys_NormDown =(TH1F*) hToys->Clone(); // Utiliser les bords de bins définis

    hToys_NormUp->Scale(1.5); // Normaliser l'histogramme
    hToys_NormDown->Scale(0.5); // Normaliser à l'intégrale souhaitée

    // hsolveToysNormUp 
    // hsolveToysNormDown 
    
        // Sauvegarde dans un fichier
    TFile* outFile = new TFile(outputfileName+".root", "RECREATE");
    hToys->Write();
    pdf->Write(); // si tu veux la récupérer plus tard
    outFile->Close();

            // Sauvegarde dans un fichier
    TFile* outFileUp = new TFile(outputfileName+"_NormUp.root", "RECREATE");
    hToys_NormUp->Write();
    pdf->Write(); // si tu veux la récupérer plus tard
    outFileUp->Close();


        // Sauvegarde dans un fichier
    TFile* outFileDown = new TFile(outputfileName+"_NormDown.root", "RECREATE");
    hToys_NormDown->Write();
    pdf->Write(); // si tu veux la récupérer plus tard
    outFileDown->Close();


    std::cout << " Toys générés et sauvegardés dans " << outputfileName << std::endl;
}
