#include <iostream>
#include <TROOT.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TString.h>

void RunPlot3DnSecInt() {

gROOT->LoadMacro("plot3DSecInt.C");

TString YEAR[1] = {"2018"};
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


// for (unsigned int j = 0 ; j < 7; j ++)
//     {
//         msmu = Msmu[j];
//         INTmsmu = INTMsmu[j];
//         for (unsigned int u = 0; u < 13 ; u++)
//             {
//                 mneu = Mneu[u];
//                 INTmneu = INTMneu[u];
//                 if ( (INTmneu < INTmsmu  && INTmneu < INTmsmu &&(INTmsmu-INTmneu)% 50 == 0) )//% 50 to get all the samples // (INTmneu < INTmsmu && ((INTmsmu-INTmneu)==20 || INTmneu == 180 )) || (INTmneu < INTmsmu &&(INTmsmu-INTmneu)% 100 == 0)
//                     {
//                         TString NAME = "RPV_2018_smu"+msmu+"_neu"+mneu;
//                         Sample.push_back(NAME);
//                         std::cout<<" Sample : "<<NAME<<std::endl;
//                     }
                
//             }            
//     }

        Sample.push_back("RPV_2018_smu500_neu200");
        Sample.push_back("RPV_2018_smu500_neu250");
        Sample.push_back("RPV_2018_smu500_neu300");
        Sample.push_back("RPV_2018_smu500_neu350");
        Sample.push_back("RPV_2018_smu500_neu400");
        Sample.push_back("RPV_2018_smu500_neu450");
        // //-- Core
        for (unsigned int i = 0; i < 1; ++i) { // loop on the years (relevant only for the data)
                for (unsigned int l = 0; l < 1; ++l) { // loop on the PU (empty, _PU25, _PU30, _PU35, _PU40, _PU45, _PU50)
                    plot(Sample,YEAR[i], PU[l]);
                }
            }
        // // -- End of Core
        Sample.clear();
}

       


