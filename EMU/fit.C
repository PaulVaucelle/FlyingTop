#include <TMath.h>

void fit() {
    // Ouvrir le fichier ROOT d'entrée
    TString inputFile = "ABCD_EMU__STW_2VtxAll.root";
    TString outputFile = "Fit_ABCD_STW_2VtxAll.root";
    TString histName = "hsolve" ;
    gROOT->SetBatch(kTRUE); // Désactiver l'affichage interactif
    TFile *file = TFile::Open(inputFile);
    if (!file || file->IsZombie()) {
        std::cout << "Erreur lors de l'ouverture du fichier " << inputFile << std::endl;
        return;
    }

    // Charger l'histogramme depuis le fichier
    TH1F *hist = (TH1F*)file->Get(histName);
    if (!hist) {
        std::cout << "L'histogramme " << histName << " n'a pas été trouvé dans le fichier." << std::endl;
        return;
    }

    // !! functions
    TF1 *gauss = new TF1("gauss", "gaus(0)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    TF1 *landau = new TF1("landau", "landau(0)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    TF1 *crystallball = new TF1("crystallball", "[0]*ROOT::Math::crystalball_function(x, [1], [2], [3], [4])", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    // !! end fit fucntions
    crystallball->SetParameters(1,1.5, -3, 16,1., 1);


    TF1 *exp = new TF1("exp", "[0]*exp(-x/[1]) + [2]", 3, 20);
    exp->SetParameters(1, 1, 0);

    TF1Convolution *convol = new TF1Convolution("gauss", 
        "ROOT::Math::crystalball_function", -10, 10, true);
    convol->SetNofPointsFFT(1000);
        // Créer la TF1 convoluée
    TF1 *fconv = new TF1("fconv", *convol, -10, 10, convol->GetNpar());
    
    // Initialiser les paramètres (optionnel mais souvent nécessaire)
    fconv->SetParameters(1, 1.5, 1,  -3, 16,1., 1);
    fconv->SetParNames("gaus_norm", "gaus_mean", "gaus_sigma", 
                       "cb_alpha", "cb_n", "cb_mean", "cb_sigma");



    TF1 *fitFunc = new TF1("fitFunc", "landau(x)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    fitFunc->SetParameters(200, hist->GetMean(), hist->GetRMS());
    // Initialiser les paramètres du fit

    fitFunc->SetParNames( "landau_norm", "landau_mean", "landau_width");//"gauss_norm", "gauss_mean", "gauss_sigma",

    // Fit l'histogramme
    // hist->Fit(crystallball, "R");
    hist->Fit(fitFunc, "R");

    // Afficher le fit et l'histogramme
    TCanvas *canvas = new TCanvas("canvas", "Fit", 800, 600);
    hist->Draw();
    hist->GetYaxis()->SetRangeUser(0, hist->GetMaximum() * 1.5);
    fitFunc->Draw("same");

    // Créer un fichier de sortie et enregistrer l'histogramme et le fit
    TFile *outFile = TFile::Open(outputFile, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cout << "Erreur lors de la création du fichier " << outputFile << std::endl;
        return;
    }

    // Sauvegarder l'histogramme et la fonction dans le fichier
    hist->Write();
    fitFunc->Write();
    
    // Fermer les fichiers
    outFile->Close();
    file->Close();

    std::cout << "Fit et histogramme enregistrés dans " << outputFile << std::endl;
}
