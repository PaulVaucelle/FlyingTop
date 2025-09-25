// Include ROOT headers
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <iostream>
#include <vector>

// Function to gather histograms from different files into one

// We take all the signal samples with all the systematics and the data samples
// We take the histograms for the different regions
// We gather them into one file for each region =>  input for datacards of Combine

void GatherHistograms() {

    TString GlobalPath = "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/";
    // !! // Parameters to change
    TString Channel = "MUMU";// MUMU ou EMU
    TString Variable = "SumtrackWeight";// SumtrackWeight ou Mass
    TString nVtx = "2VtxAll";//2Vtx or 2VtxAll
    TString AltVtx = "TLVtxAll";
    TString Sign = "OS";
    TString MU = "DM";
    TString etamax = "2p4"; 
    TString ctau ="100";
    TString YEAR = "2018";
    TString DATANAME = "DoubleMuon_UL18";
    // !! // -------------------

    // !! Production directories 
    TString ProdSignal              = "Signal_"+YEAR;// SYST_CTAU100  ou SYST_EMU_CTAU100
    // TString ProdSignalLumiUp        = "Signal_"+YEAR+"_LumiUp";
    // TString ProdSignalLumiDown      = "Signal_"+YEAR+"_LumiDown";
    // TString ProdSignalXSUp          = "Signal_"+YEAR+"_XSUp";
    // TString ProdSignalXSDown        = "Signal_"+YEAR+"_XSDown";
    TString ProdSignalJECUp         = "Signal_"+YEAR+"_JECUp";
    TString ProdSignalJECDown       = "Signal_"+YEAR+"_JECDown";
    TString ProdSignalJERUp         = "Signal_"+YEAR+"_JERUp";
    TString ProdSignalJERDown       = "Signal_"+YEAR+"_JERDown";
    TString ProdData                = "DATA_MUMU_"+YEAR+"_19_08_2024";

    if (Channel == "EMU")
        {
            ProdSignal = "SYST_EMU_CTAU100";
            ProdData = "";
        }
    if (Channel == "MUMU") // !! Production direcotries relative to the year and the global path :D
        {
            ProdSignal              = "Signal_"+YEAR;
            // ProdSignalLumiUp        = "Signal_"+YEAR+"_LumiUp";
            // ProdSignalLumiDown      = "Signal_"+YEAR+"_LumiDown";
            // ProdSignalXSUp          = "Signal_"+YEAR+"_XSUp";
            // ProdSignalXSDown        = "Signal_"+YEAR+"_XSDown";
            ProdSignalJECUp         = "Signal_"+YEAR+"_JECUp";
            ProdSignalJECDown       = "Signal_"+YEAR+"_JECDown";
            ProdSignalJERUp         = "Signal_"+YEAR+"_JERUp";
            ProdSignalJERDown       = "Signal_"+YEAR+"_JERDown";
             ProdSignalRoccorDown       = "Signal_"+YEAR+"_RoccorDown";
            //MuonISOUp, MuonISODown, 
            //MuonIDUp, MuonIDDown, 
            //MuonTrigUp, MuonTrigDown

            ProdData = "DATA_MUMU_"+YEAR+"_19_08_2024";
        }
    std::vector<TString> ProdSignalSyst = {
        ProdSignal, 
        ProdSignal,ProdSignal, //lumi
        ProdSignal,ProdSignal, // L1
        ProdSignal,ProdSignal, //Trigger
        ProdSignal,ProdSignal, //Pu
        ProdSignal,ProdSignal, // SDEle
        ProdSignal,ProdSignal, // TopPt
        ProdSignal,ProdSignal, // PDF
        ProdSignal,ProdSignal, // Scale
        ProdSignalJECUp, ProdSignalJECDown,  // JEC
        ProdSignalJERUp, ProdSignalJERDown // JER
        ProdSignal,ProdSignalRoccorDown // Roccor
        ProdSignal, ProdSignal, //MuonISO
        ProdSignal, ProdSignal, //MuonID


        };


    // !! ----------------------------------------- 

    // !! ------------ Signal Samples ------------------ 
        TString SignalSet[34]={"RPV_"+YEAR+"_smu200_neu180_ctau"+ctau,"RPV_"+YEAR+"_smu250_neu180_ctau"+ctau,"RPV_"+YEAR+"_smu250_neu200_ctau"+ctau,
    "RPV_"+YEAR+"_smu250_neu230_ctau"+ctau,"RPV_"+YEAR+"_smu300_neu180_ctau"+ctau,"RPV_"+YEAR+"_smu300_neu200_ctau"+ctau,"RPV_"+YEAR+"_smu300_neu250_ctau"+ctau,"RPV_"+YEAR+"_smu300_neu280_ctau"+ctau,
    "RPV_"+YEAR+"_smu350_neu180_ctau"+ctau,"RPV_"+YEAR+"_smu350_neu200_ctau"+ctau,"RPV_"+YEAR+"_smu350_neu250_ctau"+ctau,"RPV_"+YEAR+"_smu350_neu300_ctau"+ctau,"RPV_"+YEAR+"_smu350_neu330_ctau"+ctau,
    "RPV_"+YEAR+"_smu400_neu180_ctau"+ctau,"RPV_"+YEAR+"_smu400_neu200_ctau"+ctau,"RPV_"+YEAR+"_smu400_neu250_ctau"+ctau,"RPV_"+YEAR+"_smu400_neu300_ctau"+ctau,"RPV_"+YEAR+"_smu400_neu350_ctau"+ctau,
    "RPV_"+YEAR+"_smu400_neu380_ctau"+ctau,"RPV_"+YEAR+"_smu450_neu180_ctau"+ctau,"RPV_"+YEAR+"_smu450_neu200_ctau"+ctau,"RPV_"+YEAR+"_smu450_neu250_ctau"+ctau,"RPV_"+YEAR+"_smu450_neu300_ctau"+ctau,
    "RPV_"+YEAR+"_smu450_neu350_ctau"+ctau,"RPV_"+YEAR+"_smu450_neu400_ctau"+ctau,"RPV_"+YEAR+"_smu450_neu430_ctau"+ctau,"RPV_"+YEAR+"_smu500_neu180_ctau"+ctau,"RPV_"+YEAR+"_smu500_neu200_ctau"+ctau,
    "RPV_"+YEAR+"_smu500_neu250_ctau"+ctau,"RPV_"+YEAR+"_smu500_neu300_ctau"+ctau,"RPV_"+YEAR+"_smu500_neu350_ctau"+ctau,"RPV_"+YEAR+"_smu500_neu400_ctau"+ctau,"RPV_"+YEAR+"_smu500_neu450_ctau"+ctau,
    "RPV_"+YEAR+"_smu500_neu480_ctau"+ctau};
    // !! ---------------------------------------------


    // !! Systematics 
    std::vector<TString> SYST= {
        "NOM",
        "LumiUp", "LumiDown",
        "L1Up", "L1Down",
        "TriggerUp", "TriggerDown",
        "PUUp", "PUDown",
        "SFEleUp", "SFEleDown",
        "TopPtUp", "TopPtDown",
        "PDFUp", "PDFDown",
        "ScaleUp", "ScaleDown",
        "JECUp", "JECDown",
        "JERUp", "JERDown",
        "RoccorUp", "RoccorDown"
        "MuonISOUp", "MuonISODown", 
        "MuonIDUp", "MuonIDDown"
    };

      TString DataSet[1]={ "DoubleMuon_UL2018_MiniAODv2_GT36-v1"
    };



    // !! à changer pour les données EMU
    if (Channel == "EMU")
        {
            DataSet[0] = "emu_2018A";

        }

    // => fileNames contient l'ensemble des paths des fichiers à lire

    // !! Name of the histograms to gather from the ABCDEFGHI regions (can be changed)
    TString htitleA = "hData_CRtightlowlowpt_"+AltVtx+"_"+Variable;
    TString htitleB = "hData_CRlooselooselowlowpt_"+AltVtx+"_"+Variable;
    TString htitleC = "hData_CRtighthighpt_"+nVtx+"_"+Variable;
    TString htitleD = "hData_CRlooselooselowpt_"+nVtx+"_"+Variable;

    // !!  Name of the output file (that will be an input of the datacards so do not change the name)
    std::vector<TString> REGIONS = {

        "Tight_LowLowPT_2Vtx_control_region",    //A
        "LooseLoose_LowLowPT_2Vtx_control_region",   //B
        "Tight_HighPT_2Vtx_signal_region",   //C
        "LooseLoose_LowPT_2Vtx_control_region"    //D

        };

    // !!  Loop over all regions to gather all histos into one file

        for (unsigned int i = 0 ; i < REGIONS.size() ; i++) // There are 9 regions to loop over for the 2 Vtx category : ABCDEFGHI ( 2Vtx and 2VtxAll )
            {
                TString outputFileName = REGIONS[i]+".root";
                TFile* outputFile = new TFile(outputFileName, "RECREATE");

                for (unsigned int j = 0 ; j < SYST.size(); j++) // loop over the systematics for signal samples
                    {
                        TString Prod = ProdSignalSyst[j];
                        for (unsigned int l = 0 ; l < 34 ; l++) // loop over the signal samples
                            {
                                TString Signal = SignalSet[l];
                                TString file = GlobalPath+Prod+"/histofile_"+MU+"_"+Sign+"_"+etamax+"_"+Signal+"_"+SYST[j]+".root";
                                
                                TString shortname = Signal+"_"+SYST[j];
                                std::cout << shortname << std::endl;
                                
                                // Open each file
                                TFile* inputFile = TFile::Open(file);

                                if (!inputFile || inputFile->IsZombie()) {
                                    std::cerr << "Error opening file: " << file << std::endl;
                                    continue;
                                }

                                // 
                                // inputFile->ls();


                                htitleA = shortname+"_hData_CRtightlowlowpt_"+AltVtx+"_"+Variable+"_";
                                htitleB = shortname+"_hData_CRlooselooselowlowpt_"+AltVtx+"_"+Variable+"_";
                                htitleC = shortname+"_hData_CRtighthighpt_"+nVtx+"_"+Variable+"_";
                                htitleD = shortname+"_hData_CRlooselooselowpt_"+nVtx+"_"+Variable+"_";

                            
                                std::vector<TString> HistoNames = {
                                    
                                        htitleA, 
                                        htitleB, 
                                        htitleC,
                                        htitleD                                    
                                    };
                                
                                // Get the histograms
                                // for (unsigned int m = 0 ; m < HistoNames.size() ; m++)
                                //     {
                                        // Retrieve the histogram
                                        inputFile->cd();
                                        TH1F* hist = (TH1F*)gROOT->FindObject(HistoNames[i]) ; //= (TH1F*)inputFile->Get(HistoNames[m]);
                                        // std::cout << "Getting  histogram: " << HistoNames[m] << " from file: " << file << std::endl;
                                        if (!hist) {
                                            std::cerr << "Error retrieving histogram: " << HistoNames[i] << " from file: " << file << std::endl;
                                            inputFile->Close();
                                            continue;//continue
                                        }

                                        // Optionally clone the histogram if you want to keep it after closing the file
                                        TH1F* histClone = (TH1F*)hist->Clone();
                                        histClone->SetDirectory(outputFile); // Attach to output file directory
                                        histClone->SetName(Signal+"_"+REGIONS[i]+"_"+SYST[j]);
                                        outputFile->cd();
                                        // Write the histogram to the output file
                                        histClone->Write();
                                    // } // loop over histos
                                    inputFile->Close();
                                    
                            }  // loop over signal samples

                    } /// loop over the systematics for signal samples
                        
                for (unsigned int j = 0 ; j < 1 ; j++) // loop over the data samples
                    {
                        TString DATA = DataSet[j];
                        TString file = GlobalPath+ProdData+"/histofile_"+MU+"_"+Sign+"_"+etamax+"_"+DATA+"_NOM.root";
                        
                        TString shortname = DATA+"_NOM";
                        std::cout << shortname << std::endl;

                        // Open each file
                        TFile* inputFile = TFile::Open(file);

                        if (!inputFile || inputFile->IsZombie()) {
                            std::cerr << "Error opening file: " << file << std::endl;
                            continue;
                        }
                       

                        htitleA = shortname+"_hData_CRtightlowlowpt_"+AltVtx+"_"+Variable+"_";
                        htitleB = shortname+"_hData_CRlooselooselowlowpt_"+AltVtx+"_"+Variable+"_";
                        htitleC = shortname+"_hData_CRtighthighpt_"+nVtx+"_"+Variable+"_";
                        htitleD = shortname+"_hData_CRlooselooselowpt_"+nVtx+"_"+Variable+"_";
                    
                        std::vector<TString> HistoNames = {
                            
                                htitleA, 
                                htitleB, 
                                htitleC,
                                htitleD
                            
                            };
                        
                        // Get the histograms
                        // for (unsigned int m = 0 ; m < HistoNames.size() ; m++)
                            // {
                                // Retrieve the histogram
                                inputFile->cd();
                                TH1F* hist = (TH1F*)gROOT->FindObject(HistoNames[i]) ; //= (TH1F*)inputFile->Get(HistoNames[m]);
                                // std::cout << "Getting  histogram: " << HistoNames[m] << " from file: " << file << std::endl;
                                if (!hist) {
                                    std::cerr << "Error retrieving histogram: " << HistoNames[i] << " from file: " << file << std::endl;
                                    inputFile->Close();
                                    continue;//continue
                                }

                                // Optionally clone the histogram if you want to keep it after closing the file
                                TH1F* histClone = (TH1F*)hist->Clone();
                                histClone->SetDirectory(outputFile); // Attach to output file directory
                                //Change name of the histograms to suit the datacards
                                histClone->SetName(DATA+"_"+REGIONS[i]+"_"+SYST[j]);

                                outputFile->cd();
                                // Write the histogram to the output file
                                histClone->Write();
                            // } // loop over histos
                            inputFile->Close();      
                    }// Loop over data samples



                outputFile->Close();
                std::cout << "Histograms have been gathered and saved in " << outputFileName << std::endl;
                
                int returnCode = gSystem->Exec("mv "+outputFileName+" /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/data/RPV_UDD_"+YEAR+"/ABCD_"+Channel+"_"+MU+"_"+Sign+"_"+etamax+"_"+nVtx+"/ ");
                // /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_11_3_4/src/HiggsAnalysis/CombinedLimit/data/RPV_UDD_2018/MUMU_DM_OS_2p4_2VtxAll/
                // Optionally, check the return code to see if the command was successful

                if (returnCode == 0) {
                    std::cout << "mv executed successfully!" << std::endl;
                } else {
                    std::cout << "mv failed with return code: " << returnCode << std::endl;
                }
            }// loop over regions

        // pour chaque histogram, hadd les fichiers correspondants
        // bouger le tout sur le setup combine

        
}// End of Gather Histograms
