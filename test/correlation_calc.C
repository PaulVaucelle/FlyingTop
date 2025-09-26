


float correlation_calc_Signal(TString Filename, TString Histoname) {
    // Nom du fichier ROOT et de l'histogramme
    TString filename = Filename;
    TString histoname = filename+"_"+Histoname;

    // Ouvrir le fichier
    TFile* f = TFile::Open("../../Signal_2018_L1/histofile_HT100_DM_OS_2p4_"+filename+".root");
    if (!f || f->IsZombie()) {
        std::cerr << "Erreur : impossible d'ouvrir " << filename << std::endl;
        return -10;
    }

    // Récupérer l'histogramme
    TH2* h2 = dynamic_cast<TH2*>(f->Get(histoname));
    if (!h2) {
        std::cerr << "Erreur : histogramme " << histoname << " introuvable." << std::endl;
        f->Close();
        return -10;
    }

    // Calculer le facteur de corrélation
    double corr = h2->GetCorrelationFactor();
    std::cout << "Corrélation entre X et Y dans " << histoname << " : " << corr << std::endl;

    f->Close();
    return corr;
}

float correlation_calc_BKG(TString Filename, TString Histoname) {
    // Nom du fichier ROOT et de l'histogramme
    TString filename = Filename;
    TString histoname = filename+"_"+Histoname;

    // Ouvrir le fichier
    TFile* f = TFile::Open("../../DATA_EMU_2018_03_02_2025/histofile_HT100_EM_OS_2p4_"+filename+"_BDT100.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Erreur : impossible d'ouvrir " << filename << std::endl;
        return -10;
    }

    // Récupérer l'histogramme
    TH2* h2 = dynamic_cast<TH2*>(f->Get(histoname));
    if (!h2) {
        std::cerr << "Erreur : histogramme " << histoname << " introuvable." << std::endl;
        f->Close();
        return -10;
    }

    // Calculer le facteur de corrélation
    double corr = h2->GetCorrelationFactor();
    std::cout << "Corrélation entre X et Y dans " << histoname << " : " << corr << std::endl;

    f->Close();
    return corr;
}


int main()
    {
        
        std::vector<TString> filenames;
        // filenames.push_back("RPV_2018_smu200_neu180");
        // filenames.push_back("RPV_2018_smu250_neu180");
        // filenames.push_back("RPV_2018_smu250_neu200");
        // filenames.push_back("RPV_2018_smu250_neu230");
        // filenames.push_back("RPV_2018_smu300_neu180");
        // filenames.push_back("RPV_2018_smu300_neu200");
        // filenames.push_back("RPV_2018_smu300_neu250");
        // filenames.push_back("RPV_2018_smu300_neu280");
        // filenames.push_back("RPV_2018_smu350_neu180");
        // filenames.push_back("RPV_2018_smu350_neu200");
        // filenames.push_back("RPV_2018_smu350_neu250");
        // filenames.push_back("RPV_2018_smu350_neu300");
        // filenames.push_back("RPV_2018_smu350_neu330");
        // filenames.push_back("RPV_2018_smu400_neu180");
        // filenames.push_back("RPV_2018_smu400_neu200");
        // filenames.push_back("RPV_2018_smu400_neu250");
        // filenames.push_back("RPV_2018_smu400_neu300");
        // filenames.push_back("RPV_2018_smu400_neu350");
        // filenames.push_back("RPV_2018_smu400_neu380");
        // filenames.push_back("RPV_2018_smu450_neu180");
        // filenames.push_back("RPV_2018_smu450_neu200");
        // filenames.push_back("RPV_2018_smu450_neu250");
        // filenames.push_back("RPV_2018_smu450_neu300");
        // filenames.push_back("RPV_2018_smu450_neu350");
        // filenames.push_back("RPV_2018_smu450_neu400");
        // filenames.push_back("RPV_2018_smu450_neu430");
        // filenames.push_back("RPV_2018_smu500_neu180");
        // filenames.push_back("RPV_2018_smu500_neu200");
        // filenames.push_back("RPV_2018_smu500_neu250");
        // filenames.push_back("RPV_2018_smu500_neu300");
        // filenames.push_back("RPV_2018_smu500_neu350");
        // filenames.push_back("RPV_2018_smu500_neu400");
        // filenames.push_back("RPV_2018_smu500_neu450");
        // filenames.push_back("RPV_2018_smu500_neu480");
        // filenames.push_back("RPV_ctau001");
        // filenames.push_back("RPV_ctau003");
        // filenames.push_back("RPV_ctau010");
        // filenames.push_back("RPV_ctau030");
        // filenames.push_back("RPV_ctau100");
        // filenames.push_back("RPV_ctau300");
        // filenames.push_back("RPV_ctau1000");

        std::vector<TString> histonames;
        // histonames.push_back("Ratio_Hemileadingpt_2Vtx");
        // histonames.push_back("Ratio_Hemisubleadingpt_2Vtx");
        // histonames.push_back("Ratio_HemiAveragept_2Vtx");
        histonames.push_back("hData_Hemipt_VtxBDT_2Vtx");

        std::ofstream ofs ("./Correlations.txt", std::ofstream::out);

        for(size_t i = 0; i < filenames.size(); ++i) {
            for (size_t j = 0; j < histonames.size(); ++j) {
                // Call the correlation_calc function with the filename and histogram name
                ofs << "Correlation for Signal : " << filenames[i]<< " of histo : "<< histonames[j] << " is : "<< correlation_calc_Signal(filenames[i], histonames[j])<<endl;
                
            }
        }

        for (size_t j = 0; j < histonames.size(); ++j) {
            // Call the correlation_calc function with the filename and histogram name
            ofs << "Correlation for EMU Data of histo : "<< histonames[j] << " is : "<< correlation_calc_BKG("MuonEG-Run2018-UL2018_MiniAODv2_GT36-v1", histonames[j])<<endl;
            
        }


        ofs.close();

    }
