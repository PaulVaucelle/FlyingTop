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
/*---------------*/

class VTTracks {
   public:

      //Constructor
      VTTracks(vector<reco::TransientTrack> VTT) : VectorTTracks (VTT){}
      
      //Destructor
      ~VTTracks(){/*prolly wanna add something here*/}

      //Get data members
      vector<reco::TransientTrack> GetVTTracks(){return VectorTTracks;}
      <reco::TransientTrack>  TTrack(int i){return VectorTTracks[i];}

      //Vector related methods
      unsigned int Size(){return VTTracks.size();}
      void PushBack(<reco::TransientTrack> TTrack){VTTracks.push_back(TTrack);}
      
      //3D position of the vertex//
      float x(int i){return VTTracks[i].innermostMeasurementState().globalPosition().x();}
      float y(int i){return VTTracks[i].innermostMeasurementState().globalPosition().y();}
      float z(int i){return VTTracks[i].innermostMeasurementState().globalPosition().z();}


      //Normalised chi2 of the vertex (following a Kalman filter atm//
      float Nchi2(int i){return ProtoVertex[i].second.normalisedChiSquared();}
      
   private:
         // ----------member data ---------------------------
      vector<reco::TransientTrack> VectorTTracks;
      
};


