float  correlation_ratio_Signal(TString Filename, TString HistonameTight, TString HistonameLoose) {
    // Nom du fichier ROOT et des histogrammes
    TString filename = Filename;
    TString histo1name = filename+"_"+HistonameTight;
    TString histo2name = filename+"_"+HistonameLoose;

    // Ouvrir le fichier ROOT
    TFile* f = TFile::Open("../../Signal_2018_L1/histofile_HT100_DM_OS_2p4_"+filename+".root");
    if (!f || f->IsZombie()) {
        std::cerr << "Erreur : impossible d'ouvrir " << filename << std::endl;
        return -10;
    }

    // Récupérer les histogrammes
    TH2* h1 = dynamic_cast<TH2*>(f->Get(histo1name));
    TH2* h2 = dynamic_cast<TH2*>(f->Get(histo2name));

    if (!h1 || !h2) {
        std::cerr << "Erreur : histogramme(s) manquant(s)." << std::endl;
        f->Close();
        return -10;
    }

    if (h1->GetNbinsX() != h2->GetNbinsX() || h1->GetNbinsY() != h2->GetNbinsY()) {
        std::cerr << "Erreur : histogrammes de dimensions incompatibles." << std::endl;
        f->Close();
        return -10;
    }

    // Créer un nouvel histogramme pour le ratio
    TH2D* h_ratio = new TH2D("h_ratio", "Ratio robuste", 
                             h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(),
                             h1->GetNbinsY(), h1->GetYaxis()->GetXmin(), h1->GetYaxis()->GetXmax());

    // Remplir manuellement le ratio avec vérification
    for (int ix = 1; ix <= h1->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h1->GetNbinsY(); ++iy) {
            double num = h1->GetBinContent(ix, iy);
            double den = h2->GetBinContent(ix, iy);

            if (den != 0 && std::isfinite(num) && std::isfinite(den)) {
                double ratio = num / den;
                h_ratio->SetBinContent(ix, iy, ratio);
            } else {
                h_ratio->SetBinContent(ix, iy, 0); // ou NaN si tu veux marquer l'absence
            }
        }
    }

    // // Calcul de la corrélation sur les bins valides uniquement
    // std::vector<double> x_vals;
    // std::vector<double> y_vals;

    // for (int ix = 1; ix <= h_ratio->GetNbinsX(); ++ix) {
    //     for (int iy = 1; iy <= h_ratio->GetNbinsY(); ++iy) {
    //         double val = h_ratio->GetBinContent(ix, iy);
    //         if (val != 0 && std::isfinite(val)) {
    //             x_vals.push_back(h_ratio->GetXaxis()->GetBinCenter(ix));
    //             y_vals.push_back(h_ratio->GetYaxis()->GetBinCenter(iy));
    //         }
    //     }
    // }

    // if (x_vals.empty()) {
    //     std::cerr << "Erreur : aucun bin valide trouvé pour la corrélation." << std::endl;
    //     f->Close();
    //     return -10;
    // }

    // Calcul via TMath::Correlation (vecteurs de même taille)
    // double corr = TMath::Correlation(&x_vals[0], &y_vals[0], x_vals.size());
    h_ratio->SaveAs(filename+"_"+histo1name+"_"+histo2name+".root");
    double corr = h_ratio->GetCorrelationFactor();
    std::cout << "Corrélation robuste (X,Y) sur le ratio " << histo1name << "/" << histo2name << " : " << corr << std::endl;

    f->Close();
    return corr;
}


float  correlation_ratio_BKG(TString Filename, TString HistonameTight, TString HistonameLoose) {
    // Nom du fichier ROOT et des histogrammes
    TString filename = Filename;
    TString histo1name = filename+"_"+HistonameTight;
    TString histo2name = filename+"_"+HistonameLoose;

    // Ouvrir le fichier ROOT
    TFile* f = TFile::Open("../../DATA_EMU_2018_03_02_2025/histofile_HT100_EM_OS_2p4_"+filename+".root");
    if (!f || f->IsZombie()) {
        std::cerr << "Erreur : impossible d'ouvrir " << filename << std::endl;
        return -10;
    }

    // Récupérer les histogrammes
    TH2* h1 = dynamic_cast<TH2*>(f->Get(histo1name));
    TH2* h2 = dynamic_cast<TH2*>(f->Get(histo2name));

    if (!h1 || !h2) {
        std::cerr << "Erreur : histogramme(s) manquant(s)." << std::endl;
        f->Close();
        return -10;
    }

    if (h1->GetNbinsX() != h2->GetNbinsX() || h1->GetNbinsY() != h2->GetNbinsY()) {
        std::cerr << "Erreur : histogrammes de dimensions incompatibles." << std::endl;
        f->Close();
        return -10;
    }

    // Créer un nouvel histogramme pour le ratio
    TH2D* h_ratio = new TH2D("h_ratio", "Ratio robuste", 
                             h1->GetNbinsX(), h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(),
                             h1->GetNbinsY(), h1->GetYaxis()->GetXmin(), h1->GetYaxis()->GetXmax());

    // Remplir manuellement le ratio avec vérification
    for (int ix = 1; ix <= h1->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h1->GetNbinsY(); ++iy) {
            double num = h1->GetBinContent(ix, iy);
            double den = h2->GetBinContent(ix, iy);

            if (den != 0 && std::isfinite(num) && std::isfinite(den)) {
                double ratio = num / den;
                h_ratio->SetBinContent(ix, iy, ratio);
            } else {
                h_ratio->SetBinContent(ix, iy, 0); // ou NaN si tu veux marquer l'absence
            }
        }
    }

    // // Calcul de la corrélation sur les bins valides uniquement
    // std::vector<double> x_vals;
    // std::vector<double> y_vals;

    // for (int ix = 1; ix <= h_ratio->GetNbinsX(); ++ix) {
    //     for (int iy = 1; iy <= h_ratio->GetNbinsY(); ++iy) {
    //         double val = h_ratio->GetBinContent(ix, iy);
    //         if (val != 0 && std::isfinite(val)) {
    //             x_vals.push_back(h_ratio->GetXaxis()->GetBinCenter(ix));
    //             y_vals.push_back(h_ratio->GetYaxis()->GetBinCenter(iy));
    //         }
    //     }
    // }

    // if (x_vals.empty()) {
    //     std::cerr << "Erreur : aucun bin valide trouvé pour la corrélation." << std::endl;
    //     f->Close();
    //     return -10;
    // }

    // Calcul via TMath::Correlation (vecteurs de même taille)
    // double corr = TMath::Correlation(&x_vals[0], &y_vals[0], x_vals.size());
    h_ratio->SaveAs("EMU_data_"+histo1name+"_"+histo2name+".root");
    double corr = h_ratio->GetCorrelationFactor();
    std::cout << "Corrélation robuste (X,Y) sur le ratio " << histo1name << "/" << histo2name << " : " << corr << std::endl;

    f->Close();
    return corr;
}


int main()
    {
        
        std::vector<TString> filenames;
        filenames.push_back("RPV_2018_smu200_neu180");
        filenames.push_back("RPV_2018_smu250_neu180");
        filenames.push_back("RPV_2018_smu250_neu200");
        filenames.push_back("RPV_2018_smu250_neu230");
        filenames.push_back("RPV_2018_smu300_neu180");
        filenames.push_back("RPV_2018_smu300_neu200");
        filenames.push_back("RPV_2018_smu300_neu250");
        filenames.push_back("RPV_2018_smu300_neu280");
        filenames.push_back("RPV_2018_smu350_neu180");
        filenames.push_back("RPV_2018_smu350_neu200");
        filenames.push_back("RPV_2018_smu350_neu250");
        filenames.push_back("RPV_2018_smu350_neu300");
        filenames.push_back("RPV_2018_smu350_neu330");
        filenames.push_back("RPV_2018_smu400_neu180");
        filenames.push_back("RPV_2018_smu400_neu200");
        filenames.push_back("RPV_2018_smu400_neu250");
        filenames.push_back("RPV_2018_smu400_neu300");
        filenames.push_back("RPV_2018_smu400_neu350");
        filenames.push_back("RPV_2018_smu400_neu380");
        filenames.push_back("RPV_2018_smu450_neu180");
        filenames.push_back("RPV_2018_smu450_neu200");
        filenames.push_back("RPV_2018_smu450_neu250");
        filenames.push_back("RPV_2018_smu450_neu300");
        filenames.push_back("RPV_2018_smu450_neu350");
        filenames.push_back("RPV_2018_smu450_neu400");
        filenames.push_back("RPV_2018_smu450_neu430");
        filenames.push_back("RPV_2018_smu500_neu180");
        filenames.push_back("RPV_2018_smu500_neu200");
        filenames.push_back("RPV_2018_smu500_neu250");
        filenames.push_back("RPV_2018_smu500_neu300");
        filenames.push_back("RPV_2018_smu500_neu350");
        filenames.push_back("RPV_2018_smu500_neu400");
        filenames.push_back("RPV_2018_smu500_neu450");
        filenames.push_back("RPV_2018_smu500_neu480");
        filenames.push_back("RPV_ctau001");
        filenames.push_back("RPV_ctau003");
        filenames.push_back("RPV_ctau010");
        filenames.push_back("RPV_ctau030");
        filenames.push_back("RPV_ctau100");
        filenames.push_back("RPV_ctau300");
        filenames.push_back("RPV_ctau1000");

        std::vector<TString> histonamesTight;
        histonamesTight.push_back("Tight_Hemileadingpt_2Vtx");
        histonamesTight.push_back("Tight_Hemisubleadingpt_2Vtx");
        histonamesTight.push_back("Tight_HemiAveragept_2Vtx");


        std::vector<TString> histonamesLoose;
        histonamesLoose.push_back("Loose_Hemileadingpt_2Vtx");
        histonamesLoose.push_back("Loose_Hemisubleadingpt_2Vtx");
        histonamesLoose.push_back("Loose_HemiAveragept_2Vtx");

        std::ofstream ofs ("./Correlationsv2.txt", std::ofstream::out);

        for(size_t i = 0; i < filenames.size(); ++i) {
            for (size_t j = 0; j < histonamesTight.size(); ++j) {
                // Call the correlation_calc function with the filename and histogram name
                ofs << "Correlation for Signal : " << filenames[i]<< " of histo : "<< histonamesTight[j] << " is : "<< correlation_ratio_Signal(filenames[i], histonamesTight[j], histonamesLoose[j])<<endl;
                
            }
        }

        for (size_t j = 0; j < histonamesTight.size(); ++j) {
            // Call the correlation_calc function with the filename and histogram name
            ofs << "Correlation for EMU Data of histo : "<< histonamesTight[j] << " is : "<< correlation_ratio_BKG("MuonEG-Run2018-UL2018_MiniAODv2_GT36-v1",  histonamesTight[j], histonamesLoose[j])<<endl;
            
        }
        ofs.close();
    }