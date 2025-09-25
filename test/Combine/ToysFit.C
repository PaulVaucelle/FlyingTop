
#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TCanvas.h>    
#include <TLatex.h>
#include <TNamed.h>
#include <TString.h>
#include <TMath.h>
#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/PlotCMS.h"

Double_t crystalBall(Double_t *x, Double_t *par)
{
  Double_t m = x[0];
  Double_t N = par[0];   // Normalisation
  Double_t mean = par[1];
  Double_t sigma = par[2];
  Double_t alpha = par[3];
  Double_t n = par[4];

  Double_t t = (m - mean) / sigma;
  if (alpha < 0) t = -t;

  Double_t absAlpha = fabs(alpha);
  Double_t A = pow(n / absAlpha, n) * exp(-0.5 * absAlpha * absAlpha);
  Double_t B = n / absAlpha - absAlpha;

  if (t > -absAlpha)
    return N * exp(-0.5 * t * t);
  else
    return N * A * pow(B - t, -n);
}


void ToysFit() 
{
  // Ouvre le fichier contenant l'histogramme
  // !! §§ Parameters to change
  // HTL_VtxBDT_Ave_Corr 
  // HTL_EventBDT_Corr
  // HTL_STW_6Bins_Corr 

  TString Var = "HTL_EventBDT_Corr"; // diffBin_EvtBDT ; HTL_VtxBDT_Ave_Corr
  TString Year= "2022B"; // Year : "2022A", "2022b", "2023A" , "2023B","2024"
  TString Channel = "MUMU"; // Channel name : "EMU", "MUMU"

  TString inputFileName = "ABCD__"+Var+"_Data_"+Year+".root"; // Nom du fichier d'entrée _"+Channel+"

  TString histoName = "hsolve";  
  TString outputFileName = "ToysFit_"+Channel+"_"+Var+"_Data_"+Year+".root"; // Nom du fichier de sortie
  TString outHistoName = "Fit";


  TString path = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/ABCDReader/"+Channel+"/"; // Chemin vers le fichier d'entrée
  bool SpecificFit = false; // Si true, utilise un fit spécifique, sinon fit exponentiel décroissant par défaut
  bool UltraSpecificFit = false; // Si true, utilise un fit spécifique avec unen lorentnzienne down et une lorentzienne up
  // !! §§ ---------------------

  TFile* inputFile = TFile::Open(path+inputFileName, "READ");
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Erreur: impossible d'ouvMrir " << inputFileName << std::endl;
    return;
  }

  // Récupère l'histogramme
  TH1* hist = dynamic_cast<TH1*>(inputFile->Get(histoName));
  if (!hist) {
    std::cerr << "Erreur: histogramme '" << histoName << "' introuvable." << std::endl;
    inputFile->Close();
    return;
  }
    std::cout<<" Hist integrale : "<<hist->Integral()<<std::endl;
  // Fit exponentiel décroissant : f(x) = A * exp(-x / tau)
  TF1* expoFit = new TF1("expoFit", "[0]*exp(-x/[1])", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());


  if (SpecificFit)
    {
        // expoFit = new TF1("expoFit", "[0]*exp(-x/[1])+ (1/[2])*exp(-((x-[3])/[4])**2)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
        // expoFit->SetParameters(18.4, 0.36,0.038,-0.1,0.05); // paramètres initiaux
        // expoFit->SetParNames("Amplitude", "Tau","GaussNorm","Gauss_Mean","Gauss1_Sigma"); // Noms des paramètres
        // expoFit->SetParLimits(1, 0.000001, 100); // tau
        // expoFit->SetParLimits(2, 0.02, 0.04); // GaussNorm
        // expoFit->SetParLimits(3, -0.15, 0.09); // Gauss Mean
        // expoFit->SetParLimits(4, 0.02, 0.1); // Gauss std dev

        // Lorentzian fit
        expoFit = new TF1("expoFit", "[0]*exp(-x/[1])+ [4]*[2]/([2]+(x-[3])**2)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
        expoFit->SetParameters(18.4, 0.36,0.038,-0.1,10); // paramètres initiaux
        expoFit->SetParNames("Amplitude", "Tau","LorentzNorm","LorentzShift","Scale"); // Noms des paramètres
        expoFit->SetParLimits(0, 0.001, 1); // amplitude 
        expoFit->SetParLimits(1, 0.0001, 1); // tau
        expoFit->SetParLimits(2, 0.05, 1); // LorentzNorm
        expoFit->SetParLimits(3, -0.7, -0.3); // LorentzShift
        expoFit->SetParLimits(4, 0.1, 100); //Scale


        // // crystall ball
        // // hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax() 
        // expoFit = new TF1("expoFit", "[0]*exp(-x/[1])+ ROOT::Math::crystalball_function(x,[2],[3],[4],[5]) ", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
        // expoFit->SetParameters(18.4, 0.36,0 ,5 , 0.15, -0.12); // paramètres initiaux
        // expoFit->SetParNames("Amplitude", "Tau","Crystall_Switch","Crystall_power","Crystall_stddev","Crystall_Mean"); // Noms des paramètres
        // expoFit->SetParLimits(0, 10, 25); // Amplitude
        // expoFit->SetParLimits(1, 0.0001, 1); // tau

        // expoFit->SetParLimits(2, -1, 1); // CRystall switch
        // expoFit->SetParLimits(3, -10, 100); // crystall_ power
        // expoFit->SetParLimits(4, 0.001, 0.5); // crystall std dev
        // expoFit->SetParLimits(5, -0.3, 0.); // crytsall mean


        // crystall ball
        // hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax() 
        // expoFit = new TF1("expoFit", "ROOT::Math::crystalball_function(x,[0],[1],[2],[3]) ", -0.4,0.3);
        // expoFit->SetParameters(0 ,5 , 0.15, -0.12); // paramètres initiaux
        // expoFit->SetParNames("Crystall_Switch","Crystall_power","Crystall_stddev","Crystall_Mean"); // Noms des paramètres
        // expoFit->SetParLimits(0, -1, 1); // CRystall switch
        // expoFit->SetParLimits(1, 1, 100); // crystall_ power
        // expoFit->SetParLimits(2, 0.001, 0.5); // crystall std dev
        // expoFit->SetParLimits(3, -0.3, 0.3); // crytsall mean


// double ROOT::Math::crystalball_function 	( 	double 	x,
// 		double 	alpha,
// 		double 	n,
// 		double 	sigma,
// 		double 	mean = 0 )
    }
else if (UltraSpecificFit)
  {
        // Lorentzian fit
        expoFit = new TF1("expoFit", "[0]*exp(-x/[1])+ [4]*[2]/([2]+(x-[3])**2)+[5]*[6]/([6]+(x-[7])**2)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
        expoFit->SetParameters(18.4, 0.36,0.038,0.1,10,  0.01,-0.3,15     ); // paramètres initiaux
        expoFit->SetParNames("Amplitude", "Tau","LorentzNorm2","LorentzShift2","Scale2","LorentzNorm","LorentzShift","Scale"); // Noms des paramètres
        expoFit->SetParLimits(0, 5, 25); // Amplitude
        expoFit->SetParLimits(1, 0.001, 1); // tau
        expoFit->SetParLimits(2, 0.001, 100); // LorentzNorm
        expoFit->SetParLimits(3, 0, .2); // LorentzShift
        // expoFit->SetParLimits(4, 0.1, 500); //Scale

        expoFit->SetParLimits(5, 0.001, 100); // LorentzNorm
        expoFit->SetParLimits(6, -0.6, -0.2); // LorentzShift
        // expoFit->SetParLimits(7, 0.1, 500); //Scale
  }
else
    {
        expoFit->SetParameters(hist->GetMaximum(), 1.0); // paramètres initiaux
        expoFit->SetParNames("Amplitude", "Tau"); // Noms des paramètres
        expoFit->SetParLimits(0, 0.0001, 50); // Amplitude >= 0
        expoFit->SetParLimits(1, 0.000001, 1); // tau >= 0
    }
    
  TFitResultPtr fitResult = hist->Fit(expoFit, "RS"); // "R" = range, "S" = retourne TFitResultPtr
  // Préparation pour sauvegarde
  TFile* outFile = TFile::Open(outputFileName, "RECREATE");
  if (!outFile || outFile->IsZombie()) {
    std::cerr << "Erreur: impossible de créer " << outputFileName << std::endl;
    inputFile->Close();
    return;
  }

  // Sauvegarde histogramme et fit
  hist->Write();
  expoFit->Write();

  // Sauvegarde des paramètres du fit dans un TNamed
  TString paramStr;
  paramStr += Form("Amplitude = %.4g\n", expoFit->GetParameter(0));
  paramStr += Form("Tau       = %.4g\n", expoFit->GetParameter(1));
  paramStr += Form("Chi2/NDF  = %.4g\n", expoFit->GetChisquare() / expoFit->GetNDF());
  std::cout << "expoFit->GetChisquare() / expoFit->GetNDF() : " << expoFit->GetChisquare() <<" / " << expoFit->GetNDF() << std::endl;
  TNamed* fitSummary = new TNamed("FitParameters", paramStr.Data());
  fitSummary->Write();

  // Création du canvas avec le fit
  TCanvas* c1 = new TCanvas("c1", "Fit Canvas", 800, 600);
  c1->SetFillColor(10);
  c1->SetFillStyle(4000);
  c1->SetBorderSize(2);
  c1->cd();
  TPad* pad1 = new TPad("pad1","This is pad1",0.01,0.03,0.99,0.99,21);
  pad1->SetFillColor(0);
  pad1->SetBorderMode(0);
  pad1->SetFrameFillColor(10);
  pad1->Draw();
  pad1->SetLogy(0);
  pad1->SetTopMargin(0.07);
  pad1->SetBottomMargin(0.15);
  pad1->SetRightMargin(0.04);
  pad1->SetLeftMargin(0.16);
  pad1->cd();


  hist->Draw("sameE2");
  hist->SetFillStyle(3354);
  hist->SetLineColor(kGray+1);
  hist->SetFillColor(kGray+1);
  hist->SetLineStyle(1);
  hist->SetLineWidth(1);

  hist->SetTickLength(0.03, "YZ");
  hist->SetTickLength(0.03,"X");
  hist->SetLabelOffset(0.015,"X");
  hist->SetLabelOffset(0.007,"Y");
  hist->SetLabelSize(0.045, "XYZ");
  hist->SetLabelFont(42, "XYZ"); 
  hist->SetTitleSize(0.055, "XYZ"); 
  hist->SetTitleFont(42, "XYZ");
  hist->SetTitleOffset(1.2,"X"); 
  hist->SetTitleOffset(1.3,"Y");
  hist->GetXaxis()->SetTitle("Event BDT Score");
  hist->GetXaxis()->SetTitleColor(1);
  hist->GetYaxis()->SetTitle("Events");
  hist->GetYaxis()->SetTitleColor(1);
  hist->SetNdivisions(509,"XYZ");


  expoFit->Draw("same");

  TLegend *  leg = new TLegend(0.6,0.6,0.85,0.8);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.04);
  leg->SetMargin(0.2);
  leg->AddEntry(hist," Prediction","FE4");
  leg->AddEntry(expoFit, "Fit", "L");
  leg->Draw();
  PlotCMSv4(pad1,Year,true);
  // Ajout du texte (chi2/ndf) sur le plot
  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.03);
  latex.DrawLatex(0.6, 0.85, Form("Chi2/NDF = %.2f", expoFit->GetChisquare() / expoFit->GetNDF()));

  c1->Write();
  c1->SaveAs("ToysFit_"+Channel+"_"+Var+"_Data_"+Year+".pdf");

  outFile->Close();
  inputFile->Close();

  std::cout << " Fit terminé. Résultats sauvegardés dans " << outputFileName << " + image fitResult.pdf" << std::endl;
}
