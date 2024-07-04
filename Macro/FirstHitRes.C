#define FirstHitRes_cxx
#include "FirstHitRes.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "HistogramManager.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
void FirstHitRes::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L FirstHitRes.C
//      Root > FirstHitRes t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

//***********************
// Historam
//**********************

   TH1F *hData_res_r_1160 = new TH1F("hData_res_r_1160","hData_res_r_1160",100,0,10);
   TH1F *hData_res_d_1160 = new TH1F("hData_res_d_1160","hData_res_d_1160",100,0,10);
   TH1F *hData_res_z_1160 = new TH1F("hData_res_z_1160","hData_res_z_1160",100,0,10);

   TH1F *hData_res_r_1168 = new TH1F("hData_res_r_1168","hData_res_r_1168",100,0,10);
   TH1F *hData_res_d_1168 = new TH1F("hData_res_d_1168","hData_res_d_1168",100,0,10);
   TH1F *hData_res_z_1168 = new TH1F("hData_res_z_1168","hData_res_z_1168",100,0,10);

   TH1F *hData_res_r_1176 = new TH1F("hData_res_r_1176","hData_res_r_1176",100,0,10);
   TH1F *hData_res_d_1176 = new TH1F("hData_res_d_1176","hData_res_d_1176",100,0,10);
   TH1F *hData_res_z_1176 = new TH1F("hData_res_z_1176","hData_res_z_1176",100,0,10);

   TH1F *hData_res_r_1184 = new TH1F("hData_res_r_1184","hData_res_r_1184",100,0,10);
   TH1F *hData_res_d_1184 = new TH1F("hData_res_d_1184","hData_res_d_1184",100,0,10);
   TH1F *hData_res_z_1184 = new TH1F("hData_res_z_1184","hData_res_z_1184",100,0,10);

   TH1F *hData_res_r_1416 = new TH1F("hData_res_r_1416","hData_res_r_1416",100,0,10);
   TH1F *hData_res_d_1416 = new TH1F("hData_res_d_1416","hData_res_d_1416",100,0,10);
   TH1F *hData_res_z_1416 = new TH1F("hData_res_z_1416","hData_res_z_1416",100,0,10);

   TH1F *hData_res_r_1420 = new TH1F("hData_res_r_1420","hData_res_r_1420",100,0,10);
   TH1F *hData_res_d_1420 = new TH1F("hData_res_d_1420","hData_res_d_1420",100,0,10);
   TH1F *hData_res_z_1420 = new TH1F("hData_res_z_1420","hData_res_z_1420",100,0,10);

   TH1F *hData_res_r_1424 = new TH1F("hData_res_r_1424","hData_res_r_1424",100,0,10);
   TH1F *hData_res_d_1424 = new TH1F("hData_res_d_1424","hData_res_d_1424",100,0,10);
   TH1F *hData_res_z_1424 = new TH1F("hData_res_z_1424","hData_res_z_1424",100,0,10);

   TH1F *hData_res_r_1428 = new TH1F("hData_res_r_1428","hData_res_r_1428",100,0,10);
   TH1F *hData_res_d_1428 = new TH1F("hData_res_d_1428","hData_res_d_1428",100,0,10);
   TH1F *hData_res_z_1428 = new TH1F("hData_res_z_1428","hData_res_z_1428",100,0,10);

   TH1F *hData_res_r_1432 = new TH1F("hData_res_r_1432","hData_res_r_1432",100,0,10);
   TH1F *hData_res_d_1432 = new TH1F("hData_res_d_1432","hData_res_d_1432",100,0,10);
   TH1F *hData_res_z_1432 = new TH1F("hData_res_z_1432","hData_res_z_1432",100,0,10);

   TH1F *hData_res_r_1440 = new TH1F("hData_res_r_1440","hData_res_r_1440",100,0,10);
   TH1F *hData_res_d_1440 = new TH1F("hData_res_d_1440","hData_res_d_1440",100,0,10);
   TH1F *hData_res_z_1440 = new TH1F("hData_res_z_1440","hData_res_z_1440",100,0,10);

   //TOBL1
   TH1F *hData_res_r_1672 = new TH1F("hData_res_r_1672","hData_res_r_1672",100,0,10);
   TH1F *hData_res_d_1672 = new TH1F("hData_res_d_1672","hData_res_d_1672",100,0,10);
   TH1F *hData_res_z_1672 = new TH1F("hData_res_z_1672","hData_res_z_1672",100,0,10);

   //__________________________________________________________//
   //__________________________________________________________//
   //                   ENDCAP                                 //
   //__________________________________________________________//
   //__________________________________________________________//

   TH1F *hData_res_r_1288 = new TH1F("hData_res_r_1288","hData_res_r_1288",100,0,10);
   TH1F *hData_res_d_1288 = new TH1F("hData_res_d_1288","hData_res_d_1288",100,0,10);
   TH1F *hData_res_z_1288 = new TH1F("hData_res_z_1288","hData_res_z_1288",100,0,10);

   TH1F *hData_res_r_1296 = new TH1F("hData_res_r_1296","hData_res_r_1296",100,0,10);
   TH1F *hData_res_d_1296 = new TH1F("hData_res_d_1296","hData_res_d_1296",100,0,10);
   TH1F *hData_res_z_1296 = new TH1F("hData_res_z_1296","hData_res_z_1296",100,0,10);

   TH1F *hData_res_r_1304 = new TH1F("hData_res_r_1304","hData_res_r_1304",100,0,10);
   TH1F *hData_res_d_1304 = new TH1F("hData_res_d_1304","hData_res_d_1304",100,0,10);
   TH1F *hData_res_z_1304 = new TH1F("hData_res_z_1304","hData_res_z_1304",100,0,10);

   TH1F *hData_res_r_1544 = new TH1F("hData_res_r_1544","hData_res_r_1544",200,0,20);
   TH1F *hData_res_d_1544 = new TH1F("hData_res_d_1544","hData_res_d_1544",200,0,20);
   TH1F *hData_res_z_1544 = new TH1F("hData_res_z_1544","hData_res_z_1544",100,0,10);

   TH1F *hData_res_r_1548 = new TH1F("hData_res_r_1548","hData_res_r_1548",200,0,20);
   TH1F *hData_res_d_1548 = new TH1F("hData_res_d_1548","hData_res_d_1548",200,0,20);
   TH1F *hData_res_z_1548 = new TH1F("hData_res_z_1548","hData_res_z_1548",100,0,10);

   TH1F *hData_res_r_1552 = new TH1F("hData_res_r_1552","hData_res_r_1552",200,0,20);
   TH1F *hData_res_d_1552 = new TH1F("hData_res_d_1552","hData_res_d_1552",200,0,20);
   TH1F *hData_res_z_1552 = new TH1F("hData_res_z_1552","hData_res_z_1552",100,0,10);

   TH1F *hData_res_r_1556 = new TH1F("hData_res_r_1556","hData_res_r_1556",200,0,20);
   TH1F *hData_res_d_1556 = new TH1F("hData_res_d_1556","hData_res_d_1556",200,0,20);
   TH1F *hData_res_z_1556 = new TH1F("hData_res_z_1556","hData_res_z_1556",100,0,10);

   TH1F *hData_res_r_1560 = new TH1F("hData_res_r_1560","hData_res_r_1560",200,0,20);
   TH1F *hData_res_d_1560 = new TH1F("hData_res_d_1560","hData_res_d_1560",200,0,20);
   TH1F *hData_res_z_1560 = new TH1F("hData_res_z_1560","hData_res_z_1560",100,0,10);

   TH1F *hData_res_r_1564 = new TH1F("hData_res_r_1564","hData_res_r_1564",200,0,20);
   TH1F *hData_res_d_1564 = new TH1F("hData_res_d_1564","hData_res_d_1564",200,0,20);
   TH1F *hData_res_z_1564 = new TH1F("hData_res_z_1564","hData_res_z_1564",100,0,10);

   TH1F *hData_res_r_1800 = new TH1F("hData_res_r_1800","hData_res_r_1800",300,0,30);
   TH1F *hData_res_d_1800 = new TH1F("hData_res_d_1800","hData_res_d_1800",300,0,30);
   TH1F *hData_res_z_1800 = new TH1F("hData_res_z_1800","hData_res_z_1800",100,0,10);

   TH1F *hData_res_r_1804 = new TH1F("hData_res_r_1804","hData_res_r_1804",300,0,30);
   TH1F *hData_res_d_1804 = new TH1F("hData_res_d_1804","hData_res_d_1804",300,0,30);
   TH1F *hData_res_z_1804 = new TH1F("hData_res_z_1804","hData_res_z_1804",100,0,10);

   TH1F *hData_res_r_1808 = new TH1F("hData_res_r_1808","hData_res_r_1808",300,0,30);
   TH1F *hData_res_d_1808 = new TH1F("hData_res_d_1808","hData_res_d_1808",300,0,30);
   TH1F *hData_res_z_1808 = new TH1F("hData_res_z_1808","hData_res_z_1808",100,0,10);

   TH1F *hData_res_r_1812 = new TH1F("hData_res_r_1812","hData_res_r_1812",300,0,30);
   TH1F *hData_res_d_1812 = new TH1F("hData_res_d_1812","hData_res_d_1812",300,0,30);
   TH1F *hData_res_z_1812 = new TH1F("hData_res_z_1812","hData_res_z_1812",100,0,10);

   TH1F *hData_res_r_1816 = new TH1F("hData_res_r_1816","hData_res_r_1816",300,0,30);
   TH1F *hData_res_d_1816 = new TH1F("hData_res_d_1816","hData_res_d_1816",300,0,30);
   TH1F *hData_res_z_1816 = new TH1F("hData_res_z_1816","hData_res_z_1816",100,0,10);


   TH1F *hData_res_r_1820 = new TH1F("hData_res_r_1820","hData_res_r_1820",300,0,30);
   TH1F *hData_res_d_1820 = new TH1F("hData_res_d_1820","hData_res_d_1820",300,0,30);
   TH1F *hData_res_z_1820 = new TH1F("hData_res_z_1820","hData_res_z_1820",100,0,10);

   TH1F *hData_res_r_1824 = new TH1F("hData_res_r_1824","hData_res_r_1824",300,0,30);
   TH1F *hData_res_d_1824 = new TH1F("hData_res_d_1824","hData_res_d_1824",300,0,30);
   TH1F *hData_res_z_1824 = new TH1F("hData_res_z_1824","hData_res_z_1824",100,0,10);

   TH1F *hData_res_r_1828 = new TH1F("hData_res_r_1828","hData_res_r_1828",300,0,30);
   TH1F *hData_res_d_1828 = new TH1F("hData_res_d_1828","hData_res_d_1828",300,0,30);
   TH1F *hData_res_z_1828 = new TH1F("hData_res_z_1828","hData_res_z_1828",100,0,10);

   TH1F *hData_res_r_1832 = new TH1F("hData_res_r_1832","hData_res_r_1832",300,0,30);
   TH1F *hData_res_d_1832 = new TH1F("hData_res_d_1832","hData_res_d_1832",300,0,30);
   TH1F *hData_res_z_1832 = new TH1F("hData_res_z_1832","hData_res_z_1832",100,0,10);

   TH1F *hData_res_r_1836 = new TH1F("hData_res_r_1836","hData_res_r_1836",300,0,30);
   TH1F *hData_res_d_1836 = new TH1F("hData_res_d_1836","hData_res_d_1836",300,0,30);
   TH1F *hData_res_z_1836 = new TH1F("hData_res_z_1836","hData_res_z_1836",100,0,10);

   TH1F *hData_res_r_1840 = new TH1F("hData_res_r_1840","hData_res_r_1840",300,0,30);
   TH1F *hData_res_d_1840 = new TH1F("hData_res_d_1840","hData_res_d_1840",300,0,30);
   TH1F *hData_res_z_1840 = new TH1F("hData_res_z_1840","hData_res_z_1840",100,0,10);

   TH1F *hData_res_r_1844 = new TH1F("hData_res_r_1844","hData_res_r_1844",300,0,30);
   TH1F *hData_res_d_1844 = new TH1F("hData_res_d_1844","hData_res_d_1844",300,0,30);
   TH1F *hData_res_z_1844 = new TH1F("hData_res_z_1844","hData_res_z_1844",100,0,10);

   TH1F *hData_res_r_1848 = new TH1F("hData_res_r_1848","hData_res_r_1848",300,0,30);
   TH1F *hData_res_d_1848 = new TH1F("hData_res_d_1848","hData_res_d_1848",300,0,30);
   TH1F *hData_res_z_1848 = new TH1F("hData_res_z_1848","hData_res_z_1848",100,0,10);



   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      unsigned int nTrk = tree_track_pt->size();

      for (unsigned int i =0; i< nTrk ; i++)
         {
            int hitpattern = tree_track_hitpattern->at(i);
            if (hitpattern == 0) continue;
            float trk_reco_x = tree_track_firstHit_x->at(i);
            float trk_reco_y = tree_track_firstHit_y->at(i);
            float trk_reco_z = tree_track_firstHit_z->at(i);
            float trk_reco_r = sqrt(trk_reco_x*trk_reco_x + trk_reco_y*trk_reco_y);
            float trk_reco_d = sqrt(trk_reco_x*trk_reco_x + trk_reco_y*trk_reco_y + trk_reco_z*trk_reco_z);

            float trk_propa_x = tree_track_propagate_x->at(i);
            float trk_propa_y = tree_track_propagate_y->at(i);
            float trk_propa_z = tree_track_propagate_z->at(i);
            float trk_propa_r = sqrt(trk_propa_x*trk_propa_x + trk_propa_y*trk_propa_y);
            float trk_propa_d = sqrt(trk_propa_x*trk_propa_x + trk_propa_y*trk_propa_y + trk_propa_z*trk_propa_z);

            // float trk_dd = sqrt((trk_reco_x - trk_propa_x)*(trk_reco_x - trk_propa_x) + (trk_reco_y - trk_propa_y)*(trk_reco_y - trk_propa_y) + (trk_reco_z - trk_propa_z)*(trk_reco_z - trk_propa_z));
            // float trk_dr =     //sqrt((trk_reco_x - trk_propa_x)*(trk_reco_x - trk_propa_x) + (trk_reco_y - trk_propa_y)*(trk_reco_y - trk_propa_y)); 
            float trk_dz = trk_reco_z - trk_propa_z;
            float trk_dy = trk_reco_y - trk_propa_y;
            float trk_dx = trk_reco_x - trk_propa_x;
            float trk_dr = trk_propa_r-trk_reco_r;
            float trk_dd =  trk_propa_d - trk_reco_d;
            // float res_d = abs(trk_dd/trk_reco_d);
            // float res_r = abs(trk_dr/trk_reco_r);
            float res_z = trk_dz;
            float res_y = trk_dy;
            float res_x = trk_dx;
            // float res_r = trk_dr;
            // float res_d = trk_dd;

            float res_r = sqrt(res_x*res_x + res_y*res_y);
            float res_d = sqrt(res_x*res_x + res_y*res_y + res_z*res_z);

            if (hitpattern == 1160)
               {
                  hData_res_r_1160->Fill(res_r);
                  hData_res_d_1160->Fill(res_d);
                  hData_res_z_1160->Fill(res_z);
               }
                        if (hitpattern == 1168)
               {
                  hData_res_r_1168->Fill(res_r);
                  hData_res_d_1168->Fill(res_d);
                  hData_res_z_1168->Fill(res_z);
               }
                        if (hitpattern == 1176)
               {
                  hData_res_r_1176->Fill(res_r);
                  hData_res_d_1176->Fill(res_d);
                  hData_res_z_1176->Fill(res_z);
               }

                           if (hitpattern == 1184)
               {
                  hData_res_r_1184->Fill(res_r);
                  hData_res_d_1184->Fill(res_d);
                  hData_res_z_1184->Fill(res_z);
               }

                           if (hitpattern == 1416)
               {
                  hData_res_r_1416->Fill(res_r);
                  hData_res_d_1416->Fill(res_d);
                  hData_res_z_1416->Fill(res_z);
               }


                           if (hitpattern == 1420)
               {
                  hData_res_r_1420->Fill(res_r);
                  hData_res_d_1420->Fill(res_d);
                  hData_res_z_1420->Fill(res_z);
               }

                           if (hitpattern == 1424)
               {
                  hData_res_r_1424->Fill(res_r);
                  hData_res_d_1424->Fill(res_d);
                  hData_res_z_1424->Fill(res_z);
               }
                           if (hitpattern == 1428)
               {
                  hData_res_r_1428->Fill(res_r);
                  hData_res_d_1428->Fill(res_d);
                  hData_res_z_1428->Fill(res_z);
               }

                           if (hitpattern == 1432)
               {
                  hData_res_r_1432->Fill(res_r);
                  hData_res_d_1432->Fill(res_d);
                  hData_res_z_1432->Fill(res_z);
               }

                           if (hitpattern == 1440)
               {
                  hData_res_r_1440->Fill(res_r);
                  hData_res_d_1440->Fill(res_d);
                  hData_res_z_1440->Fill(res_z);
               }

               //TOBL1
                           if (hitpattern == 1672)
               {
                  hData_res_r_1672->Fill(res_r);
                  hData_res_d_1672->Fill(res_d);
                  hData_res_z_1672->Fill(res_z);
               }

               //__________________________________________________________//
               //__________________________________________________________//
               //                   ENDCAP                                 //
               //__________________________________________________________//
               //__________________________________________________________//


            if (hitpattern == 1288)
               {
                  hData_res_r_1288->Fill(res_r);
                  hData_res_d_1288->Fill(res_d);
                  hData_res_z_1288->Fill(res_z);
               }

            if (hitpattern == 1296)
               {
                  hData_res_r_1296->Fill(res_r);
                  hData_res_d_1296->Fill(res_d);
                  hData_res_z_1296->Fill(res_z);
               }

            if (hitpattern == 1304)
               {
                  hData_res_r_1304->Fill(res_r);
                  hData_res_d_1304->Fill(res_d);
                  hData_res_z_1304->Fill(res_z);
               }

            if (hitpattern == 1544)
               {
                  hData_res_r_1544->Fill(res_r);
                  hData_res_d_1544->Fill(res_d);
                  hData_res_z_1544->Fill(res_z);
               }

           if (hitpattern == 1548)
               {
                  hData_res_r_1548->Fill(res_r);
                  hData_res_d_1548->Fill(res_d);
                  hData_res_z_1548->Fill(res_z);
               }

           if (hitpattern == 1552)
               {
                  hData_res_r_1552->Fill(res_r);
                  hData_res_d_1552->Fill(res_d);
                  hData_res_z_1552->Fill(res_z);
               }

           if (hitpattern == 1556)
               {
                  hData_res_r_1556->Fill(res_r);
                  hData_res_d_1556->Fill(res_d);
                  hData_res_z_1556->Fill(res_z);
               }

           if (hitpattern == 1560)
               {
                  hData_res_r_1560->Fill(res_r);
                  hData_res_d_1560->Fill(res_d);
                  hData_res_z_1560->Fill(res_z);
               }

           if (hitpattern == 1564)
               {
                  hData_res_r_1564->Fill(res_r);
                  hData_res_d_1564->Fill(res_d);
                  hData_res_z_1564->Fill(res_z);
               }

           if (hitpattern == 1800)
               {
                  hData_res_r_1800->Fill(res_r);
                  hData_res_d_1800->Fill(res_d);
                  hData_res_z_1800->Fill(res_z);
               }

          if (hitpattern == 1804)
               {
                  hData_res_r_1804->Fill(res_r);
                  hData_res_d_1804->Fill(res_d);
                  hData_res_z_1804->Fill(res_z);
               }

         if (hitpattern == 1808)
               {
                  hData_res_r_1808->Fill(res_r);
                  hData_res_d_1808->Fill(res_d);
                  hData_res_z_1808->Fill(res_z);
               }

         if (hitpattern == 1812)
               {
                  hData_res_r_1812->Fill(res_r);
                  hData_res_d_1812->Fill(res_d);
                  hData_res_z_1812->Fill(res_z);
               }

         if (hitpattern == 1816)
               {
                  hData_res_r_1816->Fill(res_r);
                  hData_res_d_1816->Fill(res_d);
                  hData_res_z_1816->Fill(res_z);
               }

         if (hitpattern == 1820)
               {
                  hData_res_r_1820->Fill(res_r);
                  hData_res_d_1820->Fill(res_d);
                  hData_res_z_1820->Fill(res_z);
               }

         if (hitpattern == 1824)
               {
                  hData_res_r_1824->Fill(res_r);
                  hData_res_d_1824->Fill(res_d);
                  hData_res_z_1824->Fill(res_z);
               }

         if (hitpattern == 1828)
               {
                  hData_res_r_1828->Fill(res_r);
                  hData_res_d_1828->Fill(res_d);
                  hData_res_z_1828->Fill(res_z);
               }

         if (hitpattern == 1832)
               {
                  hData_res_r_1832->Fill(res_r);
                  hData_res_d_1832->Fill(res_d);
                  hData_res_z_1832->Fill(res_z);
               }

         if (hitpattern == 1836)
               {
                  hData_res_r_1836->Fill(res_r);
                  hData_res_d_1836->Fill(res_d);
                  hData_res_z_1836->Fill(res_z);
               }

                                          if (hitpattern == 1840)
               {
                  hData_res_r_1840->Fill(res_r);
                  hData_res_d_1840->Fill(res_d);
                  hData_res_z_1840->Fill(res_z);
               }

                                          if (hitpattern == 1844)
               {
                  hData_res_r_1844->Fill(res_r);
                  hData_res_d_1844->Fill(res_d);
                  hData_res_z_1844->Fill(res_z);
               }

                                          if (hitpattern == 1848)
               {
                  hData_res_r_1848->Fill(res_r);
                  hData_res_d_1848->Fill(res_d);
                  hData_res_z_1848->Fill(res_z);
               }

         }


   }
   HistogramManager h ;
   h.WriteAllHistogramsInFile("./FirstHitres.root","recreate");

}
