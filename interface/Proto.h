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

class Proto {
   public:

      //Constructor
      Proto(std::vector< std::pair< std::vector<reco::TransientTrack>, TransientVertex> > prtVtx) : ProtoVertex (prtVtx){}
      
      //Destructor
      ~Proto(){/*prolly wanna add something here*/}

      //Get data members
      std::vector< std::pair< std::vector<reco::TransientTrack>, TransientVertex> > GetProto(){return ProtoVertex;}
      std::pair< std::vector<reco::TransientTrack>, TransientVertex> Pair(int i){return ProtoVertex[i];}
      std::vector<reco::TransientTrack> TTrack(int i){return ProtoVertex[i].first;}
      TransientVertex TVertex(int i){return ProtoVertex[i].second;}

      //Vector related methods
      unsigned int Size(){return GetProto().size();}
      unsigned int SizeTTracks(int i){return ProtoVertex[i].first.size(); }
      void PushBack(std::pair< std::vector<reco::TransientTrack>, TransientVertex> pair_iterators_vertex){ProtoVertex.push_back(pair_iterators_vertex);}
      
      //3D position of the vertex//
      float x(int i){return ProtoVertex[i].second.position().x();}
      float y(int i){return ProtoVertex[i].second.position().y();}
      float z(int i){return ProtoVertex[i].second.position().z();}


      //Normalised chi2 of the vertex (following a Kalman filter atm//
      float Nchi2(int i){return ProtoVertex[i].second.normalisedChiSquared();}
      
   private:
         // ----------member data ---------------------------
      std::vector< std::pair< std::vector<reco::TransientTrack>, TransientVertex> > ProtoVertex;
      vector<reco::TransientTrack> displacedTTracks;//Not used yet
};


