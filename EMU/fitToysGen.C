void genToys( ) {
    TString outputFile = "ToysGen.root";
    int nToys = 1000000; // Nombre de toys à générer
    // Définir la fonction (PDF)
    TF1* pdf = new TF1("pdf", "landau(x)", 0, 20); // ici une exponentielle entre 0 et 10
    pdf->SetParameters(2.31520e+03,1.72571,1.88699e-01); // lambda < 0 pour décroissance
    //1.15771e+01
    // Important : s'assurer que la fonction est positive
    // if (pdf->Eval(pdf->GetXmin()) <= 0 || pdf->Eval(pdf->GetXmax()) <= 0) {
    //     std::cout << "Attention : la fonction n'est pas strictement positive sur le domaine !" << std::endl;
    //     return;
    // }

    // Créer un histogramme pour les toys
    TH1F* hToys = new TH1F("hToys", "Toys from TF1 PDF", 20, 0, 20);

    // Générer les toys
    for (int i = 0; i < nToys; ++i) {
        double x = pdf->GetRandom(); // tirage selon la forme de la fonction
        hToys->Fill(x);
    }
    hToys->Scale(1./hToys->Integral(0,-1)); // Normaliser l'histogramme
    hToys->Scale(597.5806);
    // Sauvegarde dans un fichier
    TFile* outFile = new TFile(outputFile, "RECREATE");
    hToys->Write();
    pdf->Write(); // si tu veux la récupérer plus tard
    outFile->Close();

    std::cout << " Toys générés et sauvegardés dans " << outputFile << std::endl;
}
