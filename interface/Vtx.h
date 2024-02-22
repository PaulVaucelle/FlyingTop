/*----------INCLUDES-----------*/
// system include files
#include <vector>
#include <map>
#include <utility>
// user include files
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
/*---------------*/

      //----------------------------------------IAVF-----------------------------------------------//
      //                           Iterative Adaptive Vertex Fitter                                //
      //Input : Collections od displaced Tracks ordered by decreasing values of BDT => to have the //
      //        the best tracks at the beginning of the collection                                 //
      //                                                                                           //
      //Process : Within the collection of tracks, look for the first good seed made of 2 tracks,  //
      //          (good seed : vertex valid with a chi2 within a certain range (LowerLimit and     //
      //          UpperLimit)). From this seed, tracks are added one by one, and a vertex is built //
      //          at each step. The vertex is either valid or not. If valid in a certain range of  //
      //          Chi2, we keep going until there are no more tracks left                          //
      //                                                                                           //
      //Degrees of freedom: -LowerLimit (Chi2 restriction)                                         //
      //                    -UpperLimit (Chi2 restriction)                                         //
      //                    -BDtValue (Restriction on the BDT value of the tracks that formed  //
      //                     the vertex)                                                           //
      //                                                                                           //
      //Why? : Address the drop in efficiency due to MiniAOD datatier, turns out it improves the   //
      //       basic AVF implementation (even in RECO/AOD)                                         //
      //                                                                                           //
      //PS: Maximum efficiency is reached for MiniAOD when using the covariance matrix correction  //
      //-------------------------------------------------------------------------------------------//
using namespace reco;
class Vtx {
   public:

      //Constructor
      Vtx()
        {}
      
      //Destructor
      ~Vtx(){/*prolly wanna add something here*/
            Vtx_Weights.clear();
            Vtx_index.clear();}

      //Vertexing related methods


    void Vertexing(std::vector<reco::TransientTrack> VertexTracks, vector<std::pair<bool,TLorentzVector>> Track_FirstHit, bool ActivateStep = true, bool RequireGoodChi2Seed = false,bool RequireGoodChi2VertexIter = false, float Chi2down = 0., float Chi2up = 10., GlobalPoint *PV = nullptr )//return type to be chnged
        {
                static AdaptiveVertexFitter 
                theFitter_Vertex(
                 GeometricAnnealing ( sigmacut, Tini, ratio ), 
                 DefaultLinearizationPointFinder(),
                 KalmanVertexUpdator<5>(), 
                 KalmanVertexTrackCompatibilityEstimator<5>(), 
                 KalmanVertexSmoother() );
                 theFitter_Vertex.setParameters( maxshift, maxlpshift, maxstep, weightThreshold );

            // Values to be returned ----
            Vtx_ntk = VertexTracks.size();
            DCA_VTX_Meand = 0;
            int badtkhit_index = -1;
            bool success = false;
            MeanWeight=0;
            float tempMeanWeight=0;
            int ntracks = 0;
            std::vector<reco::TransientTrack> vTT;
            float tempchi2 = -10.;
            float tempx = -100.;
            float tempy = -100.;
            float tempz = -100.;
            Vtx_x = -100.;
            Vtx_y = -100.;
            Vtx_z = -100.;
            Vtx_chi = -10.;
            // Vtx_step = 0;
            Vtx_Weights.clear();
            Vtx_index.clear();
                 //-----------------//
            if (  Vtx_ntk > 1 && ActivateStep )
                {
                    for (int k = 0 ; k < Vtx_ntk-1; k++)
                        {
                            for (int p = k+1 ; p < Vtx_ntk ; p++)
                                {
                                  vTT.push_back(VertexTracks[p]);
                                  vTT.push_back(VertexTracks[k]);  
                                  ntracks = 2;
                                  TransientVertex TV = theFitter_Vertex.vertex(vTT); // We take the first "good-looking" seed to start
                                  if ( TV.isValid() )
                                    {
                                      if(RequireGoodChi2Seed && TV.normalisedChiSquared() < Chi2down && TV.normalisedChiSquared() > Chi2up) continue; 
                                      for (int m = 0; m < ntracks; m++) // we check that both tracks have their first hit "after" the vertex
                                        {
                                          if ( Track_FirstHit[m].first == true ) continue; //first hit of lost track is biased
                                          float PosFH = sqrt((Track_FirstHit[m].second.X()-PV->x())*(Track_FirstHit[m].second.X()-PV->x())+(Track_FirstHit[m].second.Y()-PV->y())*(Track_FirstHit[m].second.Y()-PV->y())+(Track_FirstHit[m].second.Z()-PV->z())*(Track_FirstHit[m].second.Z()-PV->z()));
                                          float PosVtx1 = sqrt((TV.position().x()-PV->x())*(TV.position().x()-PV->x())+(TV.position().y()-PV->y())*(TV.position().y()-PV->y())+(TV.position().z()-PV->z())*(TV.position().z()-PV->z()));
                                          if ( PosFH>PosVtx1 ) 
                                              {
                                                  success   = true; 
                                                  tempchi2  = TV.normalisedChiSquared();
                                                  tempx     = TV.position().x();
                                                  tempy     = TV.position().y();
                                                  tempz     = TV.position().z();
                                                  posError  = TV.positionError();
                                                  continue;
                                              }
                                          else badtkhit_index = m; // we keep in memory the index of the track that does not have a godd first hit
                                        }
                                      if (Vtx_ntk == 2 && !success ) break; // removing 1 track gives no other option than break
                                      else if ( success )
                                        {
                                          Vtx_index.clear();
                                          Vtx_index.push_back(k);
                                          Vtx_index.push_back(p);
                                          for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                                            {
                                              if (m == k || m == p) continue;
                                              ntracks++;
                                              vTT.push_back(VertexTracks[m]);
                                              TransientVertex updatedTV = theFitter_Vertex.vertex(vTT);
                                              if ( !updatedTV.isValid() ) 
                                                  {  
                                                    vTT.pop_back();
                                                    ntracks--;
                                                    updatedTV   = theFitter_Vertex.vertex(vTT);
                                                    tempchi2    = updatedTV.normalisedChiSquared();
                                                    tempx       = updatedTV.position().x();
                                                    tempy       = updatedTV.position().y();
                                                    tempz       = updatedTV.position().z();
                                                    posError    = updatedTV.positionError();
                                                    tempMeanWeight  = 0;
                                                    Vtx_Weights.clear();
                                                    for (int i = 0; i < ntracks; i++)
                                                      {
                                                        tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                                        Vtx_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                                        if (i == ntracks-1) MeanWeight = tempMeanWeight; 
                                                      }
                                                      continue;//not needed
                                                  }
                                              if ( updatedTV.isValid() ) 
                                                {

                                                  if(RequireGoodChi2VertexIter && updatedTV.normalisedChiSquared() < Chi2down && updatedTV.normalisedChiSquared() > Chi2up) success = false;
                                                  tempchi2  = updatedTV.normalisedChiSquared();
                                                  tempx     = updatedTV.position().x();
                                                  tempy     = updatedTV.position().y();
                                                  tempz     = updatedTV.position().z();
                                                  posError  = updatedTV.positionError();
                                                  float TPosFH = sqrt((Track_FirstHit[m].second.X()-PV->x())*(Track_FirstHit[m].second.X()-PV->x())+(Track_FirstHit[m].second.Y()-PV->y())*(Track_FirstHit[m].second.Y()-PV->y())+(Track_FirstHit[m].second.Z()-PV->z())*(Track_FirstHit[m].second.Z()-PV->z()));
                                                  float TPosVtx1 = sqrt((tempx-PV->x())*(tempx-PV->x())+(tempy-PV->y())*(tempy-PV->y())+(tempz-PV->z())*(tempz-PV->z()));
                                                  if (TPosFH>TPosVtx1 || Track_FirstHit[m].first == true) success = true; // continue not useful
                                                  else  
                                                    {
                                                      vTT.pop_back();
                                                      ntracks--;
                                                      updatedTV = theFitter_Vertex.vertex(vTT);
                                                      tempchi2  = updatedTV.normalisedChiSquared();
                                                      tempx     = updatedTV.position().x();
                                                      tempy     = updatedTV.position().y();
                                                      tempz     = updatedTV.position().z();
                                                      posError  = updatedTV.positionError();
                                                      tempMeanWeight=0;
                                                      Vtx_Weights.clear();
                                                      for(int i = 0; i < ntracks; i++)
                                                        {
                                                          tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                                          Vtx_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                                          if (i == ntracks-1) MeanWeight = tempMeanWeight; 
                                                        }
                                                          continue;
                                                    }
                                                  tempMeanWeight=0;
                                                  Vtx_Weights.clear();
                                                  Vtx_index.push_back(m);
                                                  for (int i = 0; i < ntracks; i++)
                                                    {
                                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                                      Vtx_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                                      if (i == ntracks-1) MeanWeight = tempMeanWeight; 
                                                    }
                                                }
                                            }
                                        }
                                        //----------------------------//
                                      else if (Vtx_ntk > 2 && !success)
                                        {
                                          Vtx_index.clear();
                                          if (badtkhit_index == k) 
                                            {
                                              vTT.erase(vTT.begin());
                                              ntracks--;
                                              Vtx_index.push_back(p);
                                            }
                                          else if (badtkhit_index == p) 
                                            {
                                              vTT.erase(vTT.end());
                                              ntracks--;
                                              Vtx_index.push_back(k);
                                            }
                                          else 
                                            {
                                              Vtx_index.push_back(k);
                                              Vtx_index.push_back(p);
                                            }
                                          for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                                            {
                                              if (m == k || m == p || m == badtkhit_index) continue; // we take care to not take into account the track with a wrong first hit
                                              ntracks++;//++
                                              tempMeanWeight=0;
                                              vTT.push_back(VertexTracks[m]);
                                              TransientVertex updatedTV = theFitter_Vertex.vertex(vTT);
                                              if ( !updatedTV.isValid() ) 
                                                {  
                                                  vTT.pop_back();
                                                  ntracks--;
                                                  if ( vTT.size() < 2 ) continue;
                                                  updatedTV   = theFitter_Vertex.vertex(vTT);
                                                  tempchi2    = updatedTV.normalisedChiSquared();
                                                  tempx       = updatedTV.position().x();
                                                  tempy       = updatedTV.position().y();
                                                  tempz       = updatedTV.position().z();
                                                  posError    = updatedTV.positionError();
                                                  Vtx_Weights.clear();
                                                  tempMeanWeight=0;
                                                  if ( ntracks < 2 )  success = false;
                                                  if ( ntracks >= 2 ) success = true;
                                                  for (int i=0; i<ntracks; i++)
                                                    {
                                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                                      Vtx_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                                      if (i == ntracks-1) MeanWeight = tempMeanWeight; 
                                                    }
                                                  continue;
                                                }
                                              if ( updatedTV.isValid() ) 
                                                {
                                                  if(RequireGoodChi2VertexIter && updatedTV.normalisedChiSquared() < Chi2down && updatedTV.normalisedChiSquared() > Chi2up) success = false;
                                                  tempchi2 = updatedTV.normalisedChiSquared();
                                                  tempx=updatedTV.position().x();
                                                  tempy=updatedTV.position().y();
                                                  tempz=updatedTV.position().z();
                                                  posError = updatedTV.positionError();
                                                  float TPosFH = sqrt((Track_FirstHit[m].second.X()-PV->x())*(Track_FirstHit[m].second.X()-PV->x())+(Track_FirstHit[m].second.Y()-PV->y())*(Track_FirstHit[m].second.Y()-PV->y())+(Track_FirstHit[m].second.Z()-PV->z())*(Track_FirstHit[m].second.Z()-PV->z()));
                                                  float TPosVtx1 = sqrt((tempx-PV->x())*(tempx-PV->x())+(tempy-PV->y())*(tempy-PV->y())+(tempz-PV->z())*(tempz-PV->z()));
                                                  if (TPosFH>TPosVtx1 || Track_FirstHit[m].first == true ) success = true; //continue not useful
                                                  else  
                                                    {
                                                      vTT.pop_back();
                                                      ntracks--;
                                                      if (vTT.size() < 2) continue;
                                                      updatedTV   = theFitter_Vertex.vertex(vTT);
                                                      tempchi2    = updatedTV.normalisedChiSquared();
                                                      tempx       = updatedTV.position().x();
                                                      tempy       = updatedTV.position().y();
                                                      tempz       = updatedTV.position().z();
                                                      posError    = updatedTV.positionError();
                                                      Vtx_Weights.clear();
                                                      tempMeanWeight=0;
                                                      if ( ntracks < 2 ) success = false;
                                                      else	       success = true;
                                                      for (int i=0; i<ntracks; i++)
                                                        {
                                                          tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                                          Vtx_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                                          if (i == ntracks-1){MeanWeight = tempMeanWeight;} 
                                                        }
                                                      continue;
                                                    }
                                                  Vtx_Weights.clear();
                                                  tempMeanWeight=0;
                                                  Vtx_index.push_back(m);
                                                  success = true;
                                                  for (int i=0; i<ntracks; i++)
                                                    {
                                                      tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                                      Vtx_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                                      if (i == ntracks-1) MeanWeight = tempMeanWeight; 
                                                    }
                                                }
                                              }
                                        }

                                      // We should have a Vertex after these conditions
                                      Vtx_ntk   = ntracks;
                                      Vtx_chi   = tempchi2;
                                      Vtx_x     = tempx;
                                      Vtx_y     = tempy;
                                      Vtx_z     = tempz;
                                      // Vtx_step  = 1;
                                      GlobalPoint RECOvtxPos(Vtx_x, Vtx_y, Vtx_z);
                                      DCA_VTX_Meand = 0;
                                      for (int k = 0; k< Vtx_ntk; k++)
                                        {
                                            TrajectoryStateClosestToPoint DCA_Vtx = vTT[k].trajectoryStateClosestToPoint(RECOvtxPos);
                                            if ( DCA_Vtx.isValid() ) // Be careful, all tracks are considered when looking at the DCA,
                                            //but one could look of the wieghts of the track at the same time
                                            { // The positions are given in the Global frame
                                                float pca_Vtx_x = DCA_Vtx.position().x();
                                                float pca_Vtx_y = DCA_Vtx.position().y();
                                                float pca_Vtx_z = DCA_Vtx.position().z();
                                                float refPoint_x = DCA_Vtx.referencePoint().x();
                                                float refPoint_y = DCA_Vtx.referencePoint().y();
                                                float refPoint_z = DCA_Vtx.referencePoint().z();
                                                float DCA_Vtx_x = refPoint_x-pca_Vtx_x;
                                                float DCA_Vtx_y = refPoint_y-pca_Vtx_y;
                                                float DCA_Vtx_z = refPoint_z-pca_Vtx_z;
                                                // float DCA_VTX_r = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y);
                                                float DCA_VTX_d = sqrt(DCA_Vtx_x*DCA_Vtx_x+DCA_Vtx_y*DCA_Vtx_y+DCA_Vtx_z*DCA_Vtx_z);
                                                DCA_VTX_Meand+=DCA_VTX_d;

                                                // tree_Hemi_Vtx_track_DCA_x.push_back(DCA_Vtx_x);
                                                // tree_Hemi_Vtx_track_DCA_y.push_back(DCA_Vtx_y);
                                                // tree_Hemi_Vtx_track_DCA_z.push_back(DCA_Vtx_z);
                                                // tree_Hemi_Vtx_track_DCA_r.push_back(DCA_VTX_r);
                                                // tree_Hemi_Vtx_track_DCA_d.push_back(DCA_VTX_d);
                                            }
                                        }
                                        DCA_VTX_Meand = DCA_VTX_Meand/(float)(Vtx_ntk);
                                        Tracks = vTT;
                                        if ( MeanWeight == 0 ) //<=>only two tracks in the valid vertex
                                          {
                                            Vtx_Weights.clear();
                                            for(int i = 0; i < ntracks; i++)
                                            {
                                              MeanWeight+=TV.trackWeight(vTT[i]);
                                              Vtx_Weights.push_back(TV.trackWeight(vTT[i]));
                                            }
                                          }
                                    }
                                    else
                                      {
                                        ntracks=0;
                                        vTT.clear();
                                      }
                                    if ( success ) break;
                                } // end loop on 2nd tracks
                            if ( success ) break;

                        }    // end loop on tracks
                }// end of TV.isValid()
              else 
                {
                  std::cout<<"Not enough tracks(<2) in input or step not activated :D "<<std::endl;
                }
        }   

      //Get data members

        //---AVF Parameters --/
      std::vector<double> GetAVFParameters(){return AVFParameters;}
      double Getmaxshift(){return maxshift;}
      double Getmaxstep(){return maxstep;} 
      double Getmaxlpshift(){return maxlpshift;}
      double GetweightThreshold(){return weightThreshold;}
      double Getsigmacut(){return sigmacut;}
      double GetTini(){return Tini;}
      double Getratio(){return ratio;}

      void Setmaxshift(double value ){maxshift = value;}
      void Setmaxstep(double value){maxstep = value;} 
      void Setmaxlpshift(double value){maxlpshift = value;}
      void SetweightThreshold(double value){weightThreshold = value ;}
      void Setsigmacut(double value){sigmacut = value;}
      void SetTini(double value){Tini = value;}
      void Setratio(double value){ratio = value;}
      // ---------------//

      //Vertex Information//
      float x(){return Vtx_x;}
      float y(){return Vtx_y;}
      float z(){return Vtx_z;}

      float chi(){return Vtx_chi;}
      int  step(){return Vtx_step;}
      int  nTrk(){return Vtx_ntk;}
      void SetStep(int nstep){Vtx_step = nstep;}
      std::vector<float> Vtx_TrkWeights(){return Vtx_Weights;}
      std::vector<unsigned int>  Vtx_TrkIndex(){return Vtx_index;}
      float MeanDCA(){return DCA_VTX_Meand;}
      float SumWeight(){return MeanWeight;}
      std::vector<TransientTrack> GetTTracks(){return Tracks;}
      GlobalError Vtx_PosErr(){return posError;}
      
   private:
         // ----------member data ---------------------------
             // parameters for the Adaptive Vertex Fitter (AVF)
             // These are the default parameters of the AVF from the twiki page, no need to change them, it does not do anything :)
        double maxshift        = 0.0001;
        double maxstep         = 30;
        double maxlpshift      = 0.1;
        double weightThreshold = 0.001;
        double sigmacut        = 3.;
        double Tini            = 256.;
        double ratio           = 0.25;
        std::vector<double>   AVFParameters = {maxshift,maxstep,maxlpshift,weightThreshold,sigmacut,Tini,ratio};
        //-------------//

        //--Output variables stored in the ntuple
        int   Vtx_ntk = 0;
        float DCA_VTX_Meand = 0;
        float MeanWeight=0;
        float Vtx_x = -100.;
        float Vtx_y = -100.;
        float Vtx_z = -100.;
        float Vtx_chi = -10.;
        int   Vtx_step = 0;
        GlobalError posError;
        std::vector<float> Vtx_Weights;
        std::vector<unsigned int> Vtx_index;
        std::vector<TransientTrack> Tracks;
        //-----------------//

        // vector<std::pair<float, TLorentzVector > > TrackInfo_Hemi1_mva to compute the mass
       //------------//

};



