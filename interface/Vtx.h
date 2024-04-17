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

using namespace reco;
class Vtx {
   public:

      //Constructor
      Vtx()
        {}
      
      //Destructor
      ~Vtx(){/*prolly wanna add something here*/
             }

    //Vertexing related methods
    //--------------------------------------------------------------------------------------------//

    /**
     * Performs vertexing using the Adaptive Vertex Fitter (AVF) algorithm.
     * 
     * @param VertexTracks The vector of transient tracks used for vertexing.
     * @return void
     */
    void AVFVertexing(std::vector<reco::TransientTrack> VertexTracks)//return type to be chnged
        {
          static AdaptiveVertexFitter 
          theFitter_Vertex(
                      GeometricAnnealing ( sigmacut, Tini, ratio ), 
                      DefaultLinearizationPointFinder(),
                      KalmanVertexUpdator<5>(), 
                      KalmanVertexTrackCompatibilityEstimator<5>(), 
                      KalmanVertexSmoother() );
          theFitter_Vertex.setParameters( maxshift, maxlpshift, maxstep, weightThreshold );

          TransientVertex TV_AVF = theFitter_Vertex.vertex(VertexTracks);
          Vtx_x = -100.;
          Vtx_y = -100.;
          Vtx_z = -100.;
          Vtx_chi = -10.;
          if (TV_AVF.isValid())
            {
              if (TV_AVF.normalisedChiSquared()<10 && TV_AVF.normalisedChiSquared()>0)
                {
                  Vtx_x = TV_AVF.position().x();
                  Vtx_y = TV_AVF.position().y();
                  Vtx_z = TV_AVF.position().z();
                  Vtx_chi = TV_AVF.normalisedChiSquared();
                }
            }
        }
    

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

    /**
     * Performs vertexing using the AdaptiveVertexFitter algorithm.
     * 
     * @param VertexTracks The vector of transient tracks representing the vertex tracks.
     * @param Track_FirstHit The vector of pairs indicating whether each track has its first hit after the vertex and the corresponding TLorentzVector.
     * @param ActivateStep Flag to activate the vertexing step.
     * @param RequireGoodChi2Seed Flag to require a good chi-squared value for the seed vertex.
     * @param RequireGoodChi2VertexIter Flag to require a good chi-squared value for each vertex iteration.
     * @param Chi2down The lower bound for the chi-squared value.
     * @param Chi2up The upper bound for the chi-squared value.
     * @param PV Pointer to the GlobalPoint representing the primary vertex.
     */
    void IAVFVertexing(std::vector<reco::TransientTrack> VertexTracks, vector<std::pair<bool,TLorentzVector>> Track_FirstHit, bool ActivateStep = true, bool RequireGoodChi2Seed = false,bool RequireGoodChi2VertexIter = false, float Chi2down = 0., float Chi2up = 10., GlobalPoint *PV = nullptr )//return type to be chnged
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
            // int badtkhit_index = -1;
            int badtkhit_index[2]= {-1,-1};
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
                                      if ( RequireGoodChi2Seed && (TV.normalisedChiSquared() < Chi2down || TV.normalisedChiSquared() > Chi2up) ) continue; 
                                      // if ( Track_FirstHit[k].first == true ) {success = true; }  //first hit of lost track is biased //
                                      float PosFH1 = sqrt((Track_FirstHit[k].second.X()-PV->x())*(Track_FirstHit[k].second.X()-PV->x())+(Track_FirstHit[k].second.Y()-PV->y())*(Track_FirstHit[k].second.Y()-PV->y())+(Track_FirstHit[k].second.Z()-PV->z())*(Track_FirstHit[k].second.Z()-PV->z()));
                                      float PosFH2 = sqrt((Track_FirstHit[p].second.X()-PV->x())*(Track_FirstHit[p].second.X()-PV->x())+(Track_FirstHit[p].second.Y()-PV->y())*(Track_FirstHit[p].second.Y()-PV->y())+(Track_FirstHit[p].second.Z()-PV->z())*(Track_FirstHit[p].second.Z()-PV->z()));
                                      float PosVtx1 = sqrt((TV.position().x()-PV->x())*(TV.position().x()-PV->x())+(TV.position().y()-PV->y())*(TV.position().y()-PV->y())+(TV.position().z()-PV->z())*(TV.position().z()-PV->z()));
                                      bool Pass1 = false;
                                      bool Pass2 = false;

                                      float diff1 = sqrt((Track_FirstHit[k].second.X()-TV.position().x())*(Track_FirstHit[k].second.X()-TV.position().x())+(Track_FirstHit[k].second.Y()-TV.position().y())*(Track_FirstHit[k].second.Y()-TV.position().y())+(Track_FirstHit[k].second.Z()-TV.position().z())*(Track_FirstHit[k].second.Z()-TV.position().z()));
                                      float diff2 = sqrt((Track_FirstHit[p].second.X()-TV.position().x())*(Track_FirstHit[p].second.X()-TV.position().x())+(Track_FirstHit[p].second.Y()-TV.position().y())*(Track_FirstHit[p].second.Y()-TV.position().y())+(Track_FirstHit[p].second.Z()-TV.position().z())*(Track_FirstHit[p].second.Z()-TV.position().z()));
                                      if ( PosFH1 > PosVtx1 || diff1 < 0.1 || diff1/PosVtx1 < 0.1 || Track_FirstHit[k].first == true) {Pass1 = true;}
                                      if ( PosFH2 > PosVtx1 || diff2 < 0.1 || diff2/PosVtx1 < 0.1 || Track_FirstHit[p].first == true) {Pass2 = true;}
                                      // if (Track_FirstHit[k].first == true;) {Pass1 = true; }
                                      if ( Pass1 && Pass2 ) 
                                        {
                                          success   = true; 
                                          tempchi2  = TV.normalisedChiSquared();
                                          tempx     = TV.position().x();
                                          tempy     = TV.position().y();
                                          tempz     = TV.position().z();
                                          posError  = TV.positionError();
                                          
                                        }
                                      else if (!Pass1) {badtkhit_index[0] = k;}
                                      else if (!Pass2) {badtkhit_index[1] = p;}
                                      // else if ( !(Pass1 && Pass2) ){success = false;}  useless since success is already set to false                             

                                      if (Vtx_ntk == 2 && !success ) break; // removing 1 track from a 2 track collection gives no other option than break => no Vtx
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
                                                    if (!updatedTV.isValid()) continue;
                                                    tempchi2    = updatedTV.normalisedChiSquared();
                                                    tempx       = updatedTV.position().x();
                                                    tempy       = updatedTV.position().y();
                                                    tempz       = updatedTV.position().z();
                                                    posError    = updatedTV.positionError();
                                                    tempMeanWeight  = 0;
                                                    Vtx_Weights.clear();
                                                    
                                                    ntracks = vTT.size();
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

                                                  if(RequireGoodChi2VertexIter && (updatedTV.normalisedChiSquared() < Chi2down || updatedTV.normalisedChiSquared() > Chi2up)) success = false;
                                                  tempchi2  = updatedTV.normalisedChiSquared();
                                                  tempx     = updatedTV.position().x();
                                                  tempy     = updatedTV.position().y();
                                                  tempz     = updatedTV.position().z();
                                                  posError  = updatedTV.positionError();
                                                  float TPosFH = sqrt((Track_FirstHit[m].second.X()-PV->x())*(Track_FirstHit[m].second.X()-PV->x())+(Track_FirstHit[m].second.Y()-PV->y())*(Track_FirstHit[m].second.Y()-PV->y())+(Track_FirstHit[m].second.Z()-PV->z())*(Track_FirstHit[m].second.Z()-PV->z()));
                                                  float TPosVtx1 = sqrt((tempx-PV->x())*(tempx-PV->x())+(tempy-PV->y())*(tempy-PV->y())+(tempz-PV->z())*(tempz-PV->z()));

                                                  float diff0 = sqrt((Track_FirstHit[m].second.X()-tempx)*(Track_FirstHit[m].second.X()-tempx)+(Track_FirstHit[m].second.Y()-tempy)*(Track_FirstHit[m].second.Y()-tempy)+(Track_FirstHit[m].second.Z()-tempz)*(Track_FirstHit[m].second.Z()-tempz));
                                                  if ( TPosFH > TPosVtx1 || diff0 < 0.1 || diff0/TPosVtx1 < 0.1 || Track_FirstHit[m].first == true) success = true;
                                                  else  
                                                    {
                                                      vTT.pop_back();
                                                      ntracks--;//not useful anymore
                                                      updatedTV = theFitter_Vertex.vertex(vTT);
                                                      if (!updatedTV.isValid()) continue;
                                                      tempchi2  = updatedTV.normalisedChiSquared();
                                                      tempx     = updatedTV.position().x();
                                                      tempy     = updatedTV.position().y();
                                                      tempz     = updatedTV.position().z();
                                                      posError  = updatedTV.positionError();
                                                      tempMeanWeight=0;
                                                      Vtx_Weights.clear();
                                                      ntracks = vTT.size();
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
                                                  // std::cout<<"ntracks : "<<ntracks<<" and Vtt.size : "<<vTT.size()<<std::endl;
                                                  ntracks = vTT.size();
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
                                          if (badtkhit_index[0] == k) 
                                            {
                                              vTT.erase(vTT.begin());
                                              ntracks--;
                                              Vtx_index.push_back(p);
                                            }
                                          else if (badtkhit_index[1] == p) 
                                            {
                                              vTT.erase(vTT.end());
                                              ntracks--;
                                              Vtx_index.push_back(k);
                                            }
                                          else 
                                            {
                                              //This should not happen
                                              std::cout<<"---------------> this should not happend"<<std::endl;
                                              Vtx_index.push_back(k);
                                              Vtx_index.push_back(p);
                                            }
                                          for (int m = 0; m < Vtx_ntk; m++) // We then add track by track to the vertex and check the validity of the vertex
                                            {
                                              if (m == k || m == p) continue; // we take care to not take into account the track with a wrong first hit 
                                              // even more security to avoid taking the same track twice
                                              ntracks++;
                                              tempMeanWeight=0;
                                              vTT.push_back(VertexTracks[m]);
                                              TransientVertex updatedTV = theFitter_Vertex.vertex(vTT);
                                              if ( !updatedTV.isValid() ) 
                                                {  
                                                  vTT.pop_back();
                                                  ntracks--;
                                                  if ( vTT.size() < 2 ) continue;
                                                  updatedTV   = theFitter_Vertex.vertex(vTT);
                                                  if (!updatedTV.isValid()) continue;
                                                  tempchi2    = updatedTV.normalisedChiSquared();
                                                  tempx       = updatedTV.position().x();
                                                  tempy       = updatedTV.position().y();
                                                  tempz       = updatedTV.position().z();
                                                  posError    = updatedTV.positionError();
                                                  Vtx_Weights.clear();
                                                  tempMeanWeight=0;
                                                  if ( ntracks < 2 )  success = false;
                                                  if ( ntracks >= 2 ) success = true;
                                                  ntracks = vTT.size();
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
                                                  if(RequireGoodChi2VertexIter && (updatedTV.normalisedChiSquared() < Chi2down || updatedTV.normalisedChiSquared() > Chi2up)) success = false;
                                                  tempchi2 = updatedTV.normalisedChiSquared();
                                                  tempx=updatedTV.position().x();
                                                  tempy=updatedTV.position().y();
                                                  tempz=updatedTV.position().z();
                                                  posError = updatedTV.positionError();
                                                  float TPosFH = sqrt((Track_FirstHit[m].second.X()-PV->x())*(Track_FirstHit[m].second.X()-PV->x())+(Track_FirstHit[m].second.Y()-PV->y())*(Track_FirstHit[m].second.Y()-PV->y())+(Track_FirstHit[m].second.Z()-PV->z())*(Track_FirstHit[m].second.Z()-PV->z()));
                                                  float TPosVtx1 = sqrt((tempx-PV->x())*(tempx-PV->x())+(tempy-PV->y())*(tempy-PV->y())+(tempz-PV->z())*(tempz-PV->z()));
                                                  //$$$$
                                                  /*
                                                  *         if (TPosFH>TPosVtx1 || Track_FirstHit[m].first == true ) success = true; // Track_FirstHit[m].first == true => remove lostTrack 
                                                  */
                                                  float diff0 = sqrt((Track_FirstHit[m].second.X()-tempx)*(Track_FirstHit[m].second.X()-tempx)+(Track_FirstHit[m].second.Y()-tempy)*(Track_FirstHit[m].second.Y()-tempy)+(Track_FirstHit[m].second.Z()-tempz)*(Track_FirstHit[m].second.Z()-tempz));
                                                  if ( TPosFH > TPosVtx1 || diff0 < 0.1 || diff0/TPosVtx1 < 0.1 || Track_FirstHit[m].first == true) success = true;
                                                              
                                                  else  
                                                    {
                                                      vTT.pop_back();
                                                      ntracks--;
                                                      if (vTT.size() < 2) continue;
                                                      updatedTV   = theFitter_Vertex.vertex(vTT);
                                                      if (!updatedTV.isValid()) continue; // It is not supposed to happen but somehow it happened once
                                                      tempchi2    = updatedTV.normalisedChiSquared();
                                                      // std::cout<<"here"<<std::endl;
                                                      tempx       = updatedTV.position().x();
                                                      tempy       = updatedTV.position().y();
                                                      tempz       = updatedTV.position().z();
                                                      posError    = updatedTV.positionError();
                                                      Vtx_Weights.clear();
                                                      tempMeanWeight=0;
                                                      if ( ntracks < 2 ) success = false;
                                                      else	       success = true;
                                                      ntracks = vTT.size();
                                                      ntracks = vTT.size();
                                                      for (int i=0; i<ntracks; i++)
                                                        {
                                                          tempMeanWeight+=updatedTV.trackWeight(vTT[i]);
                                                          Vtx_Weights.push_back(updatedTV.trackWeight(vTT[i]));
                                                          if (i ==ntracks-1){MeanWeight = tempMeanWeight;} 
                                                        }
                                                      continue;
                                                    }
                                                  Vtx_Weights.clear();
                                                  tempMeanWeight=0;
                                                  Vtx_index.push_back(m);
                                                  success = true;
                                                  ntracks = vTT.size();
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
                                      success = true;
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
                                            ntracks = vTT.size();
                                            for(int i = 0; i < ntracks; i++)
                                            {
                                              MeanWeight+=TV.trackWeight(vTT[i]);
                                              Vtx_Weights.push_back(TV.trackWeight(vTT[i]));
                                            }
                                          }
                                    }//end of tv is valid
                                  else
                                    {
                                      ntracks=0;
                                      vTT.clear();
                                    }
                                    if ( success ) break;
                                } // end loop on 2nd tracks
                            if ( success ) break;

                        }    // end 1st  loop on tracks
                }// end of ntk>1
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

      float chi2(){return Vtx_chi;}
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
// CMSSW_14_1_X_2024-02-22-2300/​RecoVertex/​LinearizationPointFinders/​interface/​FsmwLinearizationPointFinder.h
// 0005 
// 0006 /** A linearization point finder. It works the following way:
// 0007    *  1. Calculate in an optimal way 'n_pairs' different crossing points.
// 0008    *     Optimal in this context means the following:
// 0009    *     a. Try to use as many different tracks as possible;
// 0010    *        avoid using the same track all the time.
// 0011    *     b. Use the most energetic tracks.
// 0012    *     c. Try not to group the most energetic tracks together.
// 0013    *        Try to group more energetic tracks with less energetic tracks.
// 0014    *        We assume collimated bundles here, so this is why.
// 0015    *     d. Perform optimally. Do not sort more tracks (by total energy, see b)
// 0016    *        than necessary.
// 0017    *     e. If n_pairs >= (number of all possible combinations),
// 0018    *        do not leave any combinations out.
// 0019    *     ( a. and e. are almost but not entirely fulfilled in the current impl )



