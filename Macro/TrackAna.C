#define TrackAna_cxx
#include "TrackAna.h"
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include "../HistogramManager.h"
#include <iostream>
#include <fstream>
#include <sstream>
void TrackAna::Loop(TString Prod, TString thesample)
{
//   In a ROOT session, you can do:
//      Root > .L TrackAna.C
//      Root > TrackAna t
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
//**********************************
// Histograms
//**********************************

   TH1F* hData_TRACK_SIZE = new TH1F("tree_TRACK_SIZE", "",     100, 0, 4000);
   TH1F* hData_nTracks = new TH1F("tree_nTracks",  "",           50, 0, 50); 
   TH1F* hData_nLostTracks = new TH1F("tree_nLostTracks", "",    50, 0, 50); 
   TH1F* hData_track_ipc = new TH1F("tree_track_ipc",  "",       70, 0, 3500);
   TH1F* hData_track_lost = new TH1F("tree_track_lost", "",      2, 0, 2);
   TH1F* hData_track_px = new TH1F("tree_track_px","",           50, 0, 500);
   TH1F* hData_track_py = new TH1F("tree_track_py","",           50, 0, 500);
   TH1F* hData_track_pz = new TH1F("tree_track_pz","",           50, 0, 500);
   TH1F* hData_track_pt = new TH1F("tree_track_pt","",           100, 0, 100);
   TH1F* hData_Strack_pt = new TH1F("tree_Strack_pt","",           100, 0, 100);
   TH1F* hData_Btrack_pt = new TH1F("tree_Btrack_pt","",           100, 0, 100);
   TH1F* hData_track_eta = new TH1F("tree_track_eta","",         60, -3, 3);
   TH1F* hData_track_phi = new TH1F("tree_track_phi","",         60, -3.15, 3.15);
   TH1F* hData_track_charge = new TH1F("tree_track_charge","",   3, -1, 1);
   TH1F* hData_track_NChi2 = new TH1F("tree_track_NChi2","",     25, 0, 25);
   TH1F* hData_Strack_NChi2 = new TH1F("tree_Strack_NChi2","",     25, 0, 25);
   TH1F* hData_Btrack_NChi2 = new TH1F("tree_Btrack_NChi2","",     25, 0, 25);
   TH1F* hData_isHighPurity = new TH1F("tree_track_isHighPurity","",  2, 0, 2);
   TH1F* hData_track_dxy = new TH1F("tree_track_dxy","",         160, -40, 40);
   TH1F* hData_track_dxyError = new TH1F("tree_track_dxyError","", 100, 0, 2);
   TH1F* hData_track_drSig = new TH1F("tree_track_drSig","",       100, 0, 500);
   TH1F* hData_Strack_drSig = new TH1F("tree_Strack_drSig","",       100, 0, 500);
   TH1F* hData_Btrack_drSig = new TH1F("tree_Btrack_drSig","",       100, 0, 500);
   TH1F* hData_track_dz = new TH1F("tree_track_dz","",            50, -100, 100);
   TH1F* hData_track_dzError = new TH1F("tree_track_dzError","",  50, 0, 10 );
   TH1F* hData_track_dzSig = new TH1F("tree_track_dzSig","",       50, 0, 500 );

   TH1F* hData_track_x = new TH1F("tree_track_x","",             60, -30, 30  );
   TH1F* hData_track_y = new TH1F("tree_track_y","",             60, -30, 30 );
   TH1F* hData_track_z = new TH1F("tree_track_z","",             200, -100, 100  );
   TH1F* hData_track_firstHit = new TH1F("tree_track_firstHit","",      200, 1100, 1900 );
   TH1F* hData_track_region = new TH1F("tree_track_region","",        2, 0, 2 );
   TH1F* hData_track_firstHit_x = new TH1F("tree_track_firstHit_x","",  60, -30, 30 );
   TH1F* hData_track_firstHit_y = new TH1F("tree_track_firstHit_y","",   60, -30, 30 );
   TH1F* hData_track_firstHit_z = new TH1F("tree_track_firstHit_z","",    200, -100, 100  );
   TH1F* hData_track_iJet = new TH1F("tree_track_iJet","",          20, 0, 20 );
   TH1F* hData_track_ntrk10 = new TH1F("tree_track_ntrk10","",       60, 0, 60 );
   TH1F* hData_track_ntrk20 = new TH1F("tree_track_ntrk20","",       60, 0, 60 );
   TH1F* hData_track_ntrk30 = new TH1F("tree_track_ntrk30","",       60, 0, 60);
   TH1F* hData_track_ntrk40 = new TH1F("tree_track_ntrk40","",       60, 0, 60 );
   TH1F* hData_track_MVAval = new TH1F("tree_track_MVAval","",        201, -1, 1 );
   TH1F* hData_track_btag = new TH1F("tree_track_btag","",          20, 0, 1 );
   TH1F* hData_track_energy = new TH1F("tree_track_energy","",        40, 0, 200 );

   TH1F*          hSignalEff = new TH1F("hSignalEff","hSignalEff",201,-1,1) ;
   TH1F*          hBkgEff = new TH1F("hBkgEff","hBkgEff",201,-1,1);
   TH2F*          hBkgEff_SigEff = new TH2F("hBkgEff_SigEff","hBkfEff_SigEff",201,0,1,201,0,1);
   TH1F*          h_sig = new TH1F("h_sig","h_sig",201,-1,1);

   TH1F* hData_track_MVAval_DM20 = new TH1F("hData_track_MVAval_DM20","",        201, -1, 1 );
   TH1F* hData_track_MVAval_DM50 = new TH1F("hData_track_MVAval_DM50","",        201, -1, 1 );
   TH1F* hData_track_MVAval_DM100 = new TH1F("hData_track_MVAval_DM100","",        201, -1, 1 );
   TH1F* hData_track_MVAval_DM150 = new TH1F("hData_track_MVAval_DM150","",        201, -1, 1 );
   TH1F* hData_track_MVAval_DM200 = new TH1F("hData_track_MVAval_DM200","",        201, -1, 1 );
   TH1F* hData_track_MVAval_DM250 = new TH1F("hData_track_MVAval_DM250","",        201, -1, 1 );
   TH1F* hData_track_MVAval_DM300 = new TH1F("hData_track_MVAval_DM300","",        201, -1, 1 );
   TH1F* hData_track_MVAval_DM320 = new TH1F("hData_track_MVAval_DM320","",        201, -1, 1 );

   TH1F* hData_reco_lepton_leadingpt= new TH1F( "hData_reco_lepton_leadingpt","", 100, 0, 500);
   TH1F* hData_reco_lepton_leadingpt2 = new TH1F("hData_reco_lepton_leadingpt2","", 100, 0, 500);
   TH1F* hData_jet_leadingpt = new TH1F("hData_jet_leadingpt","", 100, 0, 500);
   TH1F* hData_jet_leadingpt2 = new TH1F("hData_jet_leadingpt2","", 100, 0, 500);
   TH1F* hData_jet_jet_dR = new TH1F("hData_jet_jet_dR","", 100, 0, 5);

   TH1F* hData_lepton_lepton_dR = new TH1F("hData_lepton_lepton_dR","", 100, 0, 5);
   TH1F* hData_lepton_lepton_dPhi = new TH1F("hData_lepton_lepton_dPhi","", 100, 0, 5);
   TH1F* hData_lepton_lepton_dEta = new TH1F("hData_lepton_lepton_dEta","", 100, 0, 5);



   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
    int nSTracks = 0;
    int nDTracks = 0;
    float bdtcut = -1;
    float delta = 0.01;
    int nStep= 2* abs(bdtcut/delta);
    float nSelecTracks[200] = {0};
    float nDiscardTracks[200] =  {0};
    float sig[200] = {0};
   Long64_t nbytes = 0, nb = 0;
   nentries = 2000000;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if ( jentry%1000 == 0 ) std::cout << "events : " << jentry << std::endl;
      // if (Cut(ientry) < 0) continue;
      // Trk BDT Analysis + potential analysis about trks at the ntuple level (not minintuple)
      // Macro should take 5-6 hours per main MC sample (DY,TT)
      int Msmu = tree_smu_mass;
      int Mneu = tree_neu_mass;
      if ( Msmu > 180 && Msmu < 220 )      Msmu = 200;
      else if ( Msmu > 230 && Msmu < 270 ) Msmu = 250;
      else if ( Msmu > 280 && Msmu < 320 ) Msmu = 300;
      else if ( Msmu > 380 && Msmu < 420 ) Msmu = 400;
      else if ( Msmu > 480 && Msmu < 520 ) Msmu = 500;
      // else cout << " !!! smu mass out of range !!! " << Msmu;
      if ( Mneu > 170 && Mneu < 190 )      Mneu = 180;
      else if ( Mneu > 190 && Mneu < 210 ) Mneu = 200;
      else if ( Mneu > 240 && Mneu < 260 ) Mneu = 250;
      else if ( Mneu > 270 && Mneu < 290 ) Mneu = 280;
      else if ( Mneu > 290 && Mneu < 310 ) Mneu = 300;
      else if ( Mneu > 340 && Mneu < 360 ) Mneu = 350;
      else if ( Mneu > 370 && Mneu < 390 ) Mneu = 380;
      else if ( Mneu > 390 && Mneu < 410 ) Mneu = 400;
      else if ( Mneu > 440 && Mneu < 460 ) Mneu = 450;
      else if ( Mneu > 470 && Mneu < 490 ) Mneu = 480;
    // std::cout<<"here"<<std::endl;
    if ( (tree_Filter || tree_FilterSameSign))
      {
        // std::cout<<"here2"<<std::endl;
        if (  tree_Mmumu > 10. && tree_njetNOmu > 0 )
          {
            // std::cout<<"here3"<<std::endl;
            // if ( tree_lepton_leadingpt->size()>0)
            //   {
                // std::cout<<"here4"<<std::endl;
                hData_reco_lepton_leadingpt->Fill(tree_lepton_leadingpt->at(0));
                // std::cout<<"here5"<<std::endl;
                hData_reco_lepton_leadingpt2->Fill(tree_lepton_leadingpt2->at(0));
                // std::cout<<"here6"<<std::endl;
                hData_jet_leadingpt->Fill(tree_jet_leadingpt->at(0));
                hData_jet_leadingpt2->Fill(tree_jet_leadingpt2->at(0));
                hData_jet_jet_dR->Fill(tree_jet_jet_dR->at(0));

                hData_lepton_lepton_dR->Fill(tree_lepton_lepton_dR->at(0));
                hData_lepton_lepton_dPhi->Fill(tree_lepton_lepton_dPhi->at(0));
                hData_lepton_lepton_dEta->Fill(tree_lepton_lepton_dEta->at(0));
              // }

          }


      }


      hData_TRACK_SIZE->Fill(tree_TRACK_SIZE);
      hData_nTracks->Fill(tree_nTracks); 
      hData_nLostTracks->Fill(tree_nLostTracks); 
      int nTrk = tree_track_pt->size();
      for (int i = 0; i < nTrk; i++) 
        {

            hData_track_ipc->Fill(tree_track_ipc->at(i));
            hData_track_lost->Fill(tree_track_lost->at(i));
            hData_track_px->Fill(tree_track_px->at(i));
            hData_track_py->Fill(tree_track_py->at(i));
            hData_track_pz->Fill(tree_track_pz->at(i));
            hData_track_pt->Fill(tree_track_pt->at(i));
              if (tree_track_sim_LLP->at(i)>0)
                {
                  hData_Strack_pt->Fill(tree_track_pt->at(i));
                  hData_Strack_NChi2->Fill(tree_track_NChi2->at(i));
                  hData_Strack_drSig->Fill(tree_track_drSig->at(i));
                }
              else
                {
                  hData_Btrack_pt->Fill(tree_track_pt->at(i));
                  hData_Btrack_NChi2->Fill(tree_track_NChi2->at(i));
                  hData_Btrack_drSig->Fill(tree_track_drSig->at(i));
                }


            hData_track_eta->Fill(tree_track_eta->at(i));
            hData_track_phi->Fill(tree_track_phi ->at(i));
            hData_track_charge->Fill(tree_track_charge->at(i));
            hData_track_NChi2->Fill(tree_track_NChi2->at(i));
            hData_isHighPurity->Fill(tree_track_isHighPurity->at(i));
            hData_track_dxy->Fill(tree_track_dxy->at(i));
            hData_track_dxyError->Fill(tree_track_dxyError->at(i));
            hData_track_drSig->Fill(tree_track_drSig->at(i));
            hData_track_dz->Fill(tree_track_dz->at(i));
            hData_track_dzError->Fill(tree_track_dzError ->at(i));
            hData_track_dzSig->Fill(tree_track_dzSig->at(i));

            hData_track_x->Fill(tree_track_x ->at(i));
            hData_track_y->Fill(tree_track_y ->at(i));
            hData_track_z->Fill(tree_track_z->at(i));
            hData_track_firstHit->Fill(tree_track_firstHit->at(i));
            hData_track_region->Fill(tree_track_region->at(i));
            hData_track_firstHit_x->Fill(tree_track_firstHit_x->at(i));
            hData_track_firstHit_y->Fill(tree_track_firstHit_y->at(i));
            hData_track_firstHit_z->Fill(tree_track_firstHit_z->at(i));
            hData_track_iJet->Fill(tree_track_iJet->at(i));
            hData_track_ntrk10->Fill(tree_track_ntrk10->at(i));
            hData_track_ntrk20->Fill(tree_track_ntrk20->at(i));
            hData_track_ntrk30->Fill(tree_track_ntrk30->at(i));
            hData_track_ntrk40->Fill(tree_track_ntrk40->at(i));
            hData_track_MVAval->Fill(tree_track_MVAval->at(i));
            hData_track_btag->Fill(tree_track_btag->at(i));
            hData_track_energy->Fill(tree_track_energy->at(i));


            if ((Msmu-Mneu) == 20){hData_track_MVAval_DM20->Fill(tree_track_MVAval->at(i));}
            if ((Msmu-Mneu) == 50){hData_track_MVAval_DM50->Fill(tree_track_MVAval->at(i));}
            if ((Msmu-Mneu) == 100){hData_track_MVAval_DM100->Fill(tree_track_MVAval->at(i));}
            if ((Msmu-Mneu) == 150){hData_track_MVAval_DM150->Fill(tree_track_MVAval->at(i));}
            if ((Msmu-Mneu) == 200){hData_track_MVAval_DM200->Fill(tree_track_MVAval->at(i));}
            if ((Msmu-Mneu) == 250){hData_track_MVAval_DM250->Fill(tree_track_MVAval->at(i));}
            if ((Msmu-Mneu) == 300){hData_track_MVAval_DM300->Fill(tree_track_MVAval->at(i));}
            if ((Msmu-Mneu) == 320){hData_track_MVAval_DM320->Fill(tree_track_MVAval->at(i));}

        }// end loop on tracks


      for (unsigned int j = 0 ; j< tree_track_MVAval->size();j++)
        {
          if (tree_track_MVAval->at(j)>-1 && tree_track_sim_LLP->at(j)>0)
              {
                nSTracks++;
              }
            else if (tree_track_MVAval->at(j)>-1 && tree_track_sim_LLP->at(j)<0)
              {
                nDTracks++;
              }
        }

      for (int i =0 ;i<nStep;i++)
         {
            float nSelecTrack = 0;
            float nDiscardTrack = 0;
            for (unsigned int j = 0 ; j< tree_track_MVAval->size();j++)
                {
                if (tree_track_MVAval->at(j)>bdtcut+i*delta && tree_track_sim_LLP->at(j)>0)
                      {
                        nSelecTracks[i]++;
                        nSelecTrack++;
                      }
                else if (tree_track_MVAval->at(j)>bdtcut+i*delta && tree_track_sim_LLP->at(j)<0)
                      {
                        nDiscardTracks[i]++;
                        nDiscardTrack++;
                      }
                }

            float SignalEff = 0.;
            float BkgEff = 0.;
            if (nSTracks !=0)
              {
                SignalEff =  nSelecTrack/nSTracks;
              }
            if (nDTracks !=0)
              {
                BkgEff = nDiscardTrack/nDTracks;
              } 
         }


   }// end loop on events

    float sig_max;
    float bdtcutopti = 0;
    float inte = 0;
    for (int i =0; i<nStep;i++)
      {
        // std::cout<<bdtcut+i*delta<<": bdtcut+i*delta& "<<nSelecTracks[i]<<" nSelecTracks[i] "<<std::endl;
        hSignalEff->Fill(bdtcut+i*delta,nSelecTracks[i]/nSTracks);
        hBkgEff->Fill(bdtcut+i*delta,nDiscardTracks[i]/nDTracks);
        hBkgEff_SigEff->Fill(nSelecTracks[i]/nSTracks,1-nDiscardTracks[i]/nDTracks);
        sig[i]=nSelecTracks[i]/sqrt(nSelecTracks[i]+nDiscardTracks[i]);
        inte += sig[i];
        if(i>0)
          {
            sig_max=sig[0];
            if(sig[i]>sig_max)
              {
                sig_max=sig[i];
                bdtcutopti = bdtcut+i*delta;
              } 
          }
        h_sig->Fill(bdtcut+i*delta,sig[i]);
      }
    HistogramManager h ;
    h.WriteAllHistogramsInFile(("../"+Prod+"/TrackAna_"+thesample+".root").Data(),"recreate");
   //   inte= inte/nStep;
}// End Loop method
