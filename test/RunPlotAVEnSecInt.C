#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlotAVEnSecInt() {

bool First = false ;

if (First) 
    {
        // --- ctau fixed ---//

            gROOT->LoadMacro("plotAVEnSecInt_ctau.C");

            TString YEAR[1] = {"2018"};
            TString SELECTION[2] = {"FullSelec","TrackerMatched"};
            TString PU[7] = {"","_PU25","_PU30","_PU35","_PU40","_PU45","_PU50"};
            std::vector<TString> Sample;

            TString Msmu[7] = {"200","250","300","350","400","450","500"};
            TString Mneu[13] = {"180","200","230","250","280","300","330","350","380","400","430","450","480"};

            int INTMsmu[7] = {200,250,300,350,400,450,500};
            int INTMneu[13] = {180,200,230,250,280,300,330,350,380,400,430,450,480};

            TString msmu ;
            TString mneu ;
            int INTmsmu = 0;
            int INTmneu = 0;

            TString ctau[7] = {"001","003","010","030","100","300","1000"};
            TString NAMEctau[7] = {"0.1","0.3","1.0","3.0","10.0","30.0","100.0"};
            TString Fctau;
            TString NAMECTAU ;
            for (unsigned int m = 0 ; m < 7 ; m++)// loop on ctau
                {   
                    Fctau = ctau[m];
                    NAMECTAU = NAMEctau[m];
                    for (unsigned int j = 0 ; j < 7; j ++)
                        {
                            msmu = Msmu[j];
                            INTmsmu = INTMsmu[j];
                            for (unsigned int u = 0; u < 13 ; u++)
                                {
                                    mneu = Mneu[u];
                                    INTmneu = INTMneu[u];
                                    if ( (INTmneu < INTmsmu  && INTmneu == 180 ) )//% 50 to get all the samples // (INTmneu < INTmsmu && ((INTmsmu-INTmneu)==20 || INTmneu == 180 )) || (INTmneu < INTmsmu &&(INTmsmu-INTmneu)% 100 == 0)
                                        {
                                            TString NAME = "RPV_2018_smu"+msmu+"_neu"+mneu+"_ctau"+Fctau;
                                            Sample.push_back(NAME);
                                            std::cout<<" Sample : "<<NAME<<std::endl;
                                        }
                                    
                                }            
                        }

                
                    // //-- Core
                    for (unsigned int i = 0; i < 1; ++i) { // loop on the years (relevant only for the data)
                            for (unsigned int k = 1; k < 2; ++k) { // loop on the selection (Selec, TrackerMatched) 
                                for (unsigned int l = 0; l < 1; ++l) { // loop on the PU (empty, _PU25, _PU30, _PU35, _PU40, _PU45, _PU50)
                                    plot(Sample,YEAR[i], SELECTION[k], PU[l],Fctau, msmu, mneu, NAMECTAU);
                                }
                            }
                        }
                    // // -- End of Core
                    Sample.clear();
                }
    }    
   
else //works
    {
        // -mass depednance --// ---------------------------

            gROOT->LoadMacro("plotAVEnSecInt_mass.C");

            TString YEAR[1] = {"2018"};
            TString SELECTION[2] = {"FullSelec","TrackerMatched"};
            TString PU[7] = {"","_PU25","_PU30","_PU35","_PU40","_PU45","_PU50"};
            std::vector<TString> Sample;

            TString Msmu[7] = {"200","250","300","350","400","450","500"};
            TString Mneu[13] = {"180","200","230","250","280","300","330","350","380","400","430","450","480"};
            TString ctau[7] = {"001","003","010","030","100","300","1000"};

            int INTMsmu[7] = {200,250,300,350,400,450,500};
            int INTMneu[13] = {180,200,230,250,280,300,330,350,380,400,430,450,480};

            TString msmu ;
            TString mneu ;
            int INTmsmu = 0;
            int INTmneu = 0;
            TString Fctau;
                    
            for (unsigned int j = 0 ; j < 7; j ++)// Msmu
                {
                    msmu = Msmu[j];
                    INTmsmu = INTMsmu[j];
                    for (unsigned int u = 0; u < 13 ; u++) // Mneu
                        {
                            mneu = Mneu[u];
                            INTmneu = INTMneu[u];
                            for (unsigned int m = 0 ; m < 7 ; m++)//ctau
                                {
                                    Fctau = ctau[m];
                                    if ( (INTmneu < INTmsmu && ((INTmsmu-INTmneu)==20 || INTmneu == 180 )) || (INTmneu < INTmsmu &&(INTmsmu-INTmneu)% 50 == 0))
                                        {
                                            TString NAME = "RPV_2018_smu"+msmu+"_neu"+mneu+"_ctau"+Fctau;
                                            Sample.push_back(NAME);
                                            std::cout<<" Sample : "<<NAME<<std::endl;
                                        }

                                }         
                            //-- Core
                            if ( (INTmneu < INTmsmu && ((INTmsmu-INTmneu)==20 || INTmneu == 180 )) || (INTmneu < INTmsmu &&(INTmsmu-INTmneu)% 50 == 0))
                                        {
                                            for (unsigned int i = 0; i < 1; ++i) { // loop on the years (relevant only for the data)
                                                    for (unsigned int k = 1; k < 2; ++k) { // loop on the selection (Selec, TrackerMatched) 
                                                        for (unsigned int l = 0; l < 1; ++l) { // loop on the PU (empty, _PU25, _PU30, _PU35, _PU40, _PU45, _PU50)
                                                            plot(Sample,YEAR[i], SELECTION[k], PU[l],Fctau, msmu, mneu);
                                                        }
                                                    }
                                                }
                                            // -- End of Core
                                            Sample.clear(); 
                                        }  
                        }
                }
    }
}

