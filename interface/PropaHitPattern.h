/*----------INCLUDES-----------*/
// system include files
#include <vector>
#include <map>
#include <utility>
// user include files
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
/*---------------*/

class PropaHitPattern{
   public:

      //Constructor
      PropaHitPattern(){        
        Layer.push_back(make_pair(1160,2.959));//PIXBL1
        Layer.push_back(make_pair(1168,6.778));//PIXBL2
        Layer.push_back(make_pair(1176,10.89));//PIXBL3
        Layer.push_back(make_pair(1184,16));//PIXBL4
        Layer.push_back(make_pair(1416,23.83));//TIBL1
        Layer.push_back(make_pair(1420,27.02));//TIBL1stereo
        Layer.push_back(make_pair(1424,32.23));//TIBL2
        Layer.push_back(make_pair(1428,35.41));//TIBL2stereo
        Layer.push_back(make_pair(1432,41.75));//TIBL3
        Layer.push_back(make_pair(1440,49.71));//TIBL4
        Layer.push_back(make_pair(1672,60.43));//TOBL1

        Disk.push_back(make_pair(1288,32.35));//PXFdisk1
        Disk.push_back(make_pair(1296,39.41));//PXFdisk2
        Disk.push_back(make_pair(1304,48.96));//PXFdisk3
        Disk.push_back(make_pair(1544,77.9));//TIDWHeel1
        Disk.push_back(make_pair(1548,80.71));//TIDWHeel1stereo
        Disk.push_back(make_pair(1552,90.4));//TIDWHeel2
        Disk.push_back(make_pair(1556,94.02));//TIDWHeel2stereo
        Disk.push_back(make_pair(1560,102.7));//TIDWHeel3
        Disk.push_back(make_pair(1564,107.0));//TIDWHeel3stereo
        Disk.push_back(make_pair(1800,131.6));//TECWHeel1
        Disk.push_back(make_pair(1804,129.4));//TECWHeel1stereo
        Disk.push_back(make_pair(1808,145.5));//TECWHeel2
        Disk.push_back(make_pair(1812,142.8));//TECWHeel2stereo
        Disk.push_back(make_pair(1816,160.1));//TECWHeel3
        Disk.push_back(make_pair(1820,157.1));//TECWHeel3stereo
        Disk.push_back(make_pair(1824,174.2));//TECWHeel4
        Disk.push_back(make_pair(1828,172.7));//TECWHeel4stereo
        Disk.push_back(make_pair(1832,188.4));//TECWHeel5
        Disk.push_back(make_pair(1836,186.3));//TECWHeel5stereo
        Disk.push_back(make_pair(1840,203.2));//TECWHeel6
        Disk.push_back(make_pair(1844,203.8));//TECWHeel6stereo
        Disk.push_back(make_pair(1848,222.2));//TECWHeel7}
      }
      // Needed Initialisation of the Tracker DAtaBase. The radius and the uncertainties are computed from the first hit of the tracks in RECO dataTier.
      // The radius is the mean value as the distributinos for each layer are not well-defined (2-3-4 peaks, etc...)
      // The standard deviation is also available for new (better?) methods.
      // stereo means second layer of a given hitpattern
      //Destructor
      ~PropaHitPattern(){/*prolly wanna add something here*/}

      //-------Main Method--------//
      //Basic3DVector<float> is a frame independant object, however, both vectors have to be given in the same Frame (In this case, IT HAS TO BE GLOBAL)
      std::pair<int,GloballyPositioned<float>::PositionType> Main(uint16_t firsthit, AnalyticalPropagator* Prop,TrajectoryStateOnSurface tsos, float eta, float phi, float vz,Basic3DVector<float> PV,Basic3DVector<float> p )
        {
          std::pair<int,GloballyPositioned<float>::PositionType> FHPosition;
          if (firsthit==1288 || firsthit==1296 || firsthit==1304 || firsthit==1544 || firsthit==1548 || firsthit==1552 || firsthit==1556 || firsthit==1560 || firsthit==1564 || firsthit==1800 || firsthit==1804 || firsthit==1808 || firsthit==1812 || firsthit==1816 || firsthit==1820 || firsthit==1824 || firsthit==1828 || firsthit==1832 || firsthit==1836 || firsthit==1840 ||firsthit== 1844 || firsthit==1848)//supposed to be plane
            {
              FHPosition = make_pair(1,PropagateToDisk( firsthit, eta, phi, vz, PV, p ));
              return FHPosition;
            }
          else
            {
              FHPosition = make_pair(0,PropagateToCylinder( firsthit, Prop, tsos, eta, phi, vz ));
              return FHPosition;
            }
        }


      //-------Propagators---------//
      //The Propagate methods could be overloaded with a FreeTrajectoryState instead of a TSOS. The FTS can also be obtained from a Transient Track 
      //Disks
      GloballyPositioned<float>::PositionType PropagateToDisk(uint16_t firsthit,const float eta,const float phi, const float vz,Basic3DVector<float> PV,Basic3DVector<float> p )
        {
          float zlayers=0;
          std::pair<bool,Basic3DVector<float>> spairPlane;
          TkRotation<float> rot(1,0,0,0,1,0,0,0,1);//Cylinder/Plane are already well-orientated => along/normal to the z-axis
          float theta=2*atan(exp(-eta));//=> needed for propagation
          for (int i=0; i<22;i++)
            {
              if (Disk[i].first==firsthit)
                {
                  zlayers = Disk[i].second;
                  if (p.z()<0){zlayers=-zlayers;}//wih eta, THis is possibliy wrong. For 100 events, 1604 tracks in the disks: 2 wrong z
                  GloballyPositioned<float>::PositionType P3D_(0.,0.,zlayers);
                  Plane P(P3D_,rot);
                  StraightLinePlaneCrossing SLPC(PV,p);//Define the initial state
                  spairPlane = SLPC.position(P);//"Propagation" 
                  Basic3DVector<float> sPlane(-1000.,-1000.,-1000.);
                  if (spairPlane.first)//Check for validity
                    { 
                      sPlane = spairPlane.second;
                      return GloballyPositioned<float>::PositionType (sPlane.x(),sPlane.y(),sPlane.z());
                    }
                  else// It may happen that the propagation fails, so we use geometry
                    {
                      float R = (zlayers-vz)*tan(theta);
                      float x0 = R*cos(phi); 
                      float y0 = R*sin(phi);
                      return GloballyPositioned<float>::PositionType (x0,y0,zlayers);
                    }
                }
              else
                {
                  continue;
                }
            }
            //This return ,never happens, just needed for compilation
            return GloballyPositioned<float>::PositionType (-1000.,-1000.,-1000.);
        }

      //Cylinder
      GloballyPositioned<float>::PositionType  PropagateToCylinder(uint16_t firsthit, AnalyticalPropagator* Prop, TrajectoryStateOnSurface tsos,const float eta, const float phi,const float vz) //useful to use it for debug
       {
          float rad = 0;
          TrajectoryStateOnSurface PropTSOS;//TSOS for the Barrel
          TkRotation<float> rot(1,0,0,0,1,0,0,0,1);//Cylinder/Plane are already well-orientated => along/normal to the z-axis
          float theta=2*atan(exp(-eta));//=> needed for propagation
          for (int i=0; i<11;i++)
            {
              if (Layer[i].first==firsthit )
                {
                    rad = Layer[i].second;
                    Cylinder Cylind(rad);
                    PropTSOS = Prop->propagate(tsos,Cylind);
                    if (PropTSOS.isValid())//Propagator 
                      {
                        return GloballyPositioned<float>::PositionType (PropTSOS.globalPosition().x(),PropTSOS.globalPosition().y(),PropTSOS.globalPosition().z());
                      }
                    else //Propagator can fail => use geometry
                      {
                        float z0 = (rad+vz*tan(theta))/tan(theta);
                        float x0 = rad*cos(phi);
                        float y0 = rad*sin(phi); 
                        return GloballyPositioned<float>::PositionType (x0,y0,z0);
                      }
                }
              else 
                {
                  continue;
                }
            }
             
            //This return never happens, just needed for compilation
            return GloballyPositioned<float>::PositionType (-1000.,-1000.,-1000.);
        }

      //-----Access Data Members------//
      float* CylDB(){return Cyl;}
      float* diskDB(){return disk;} 

   private:
      // ----------member data ---------------------------
        float Cyl[11]={2.959,6.778,10.89,16,23.83,27.02,32.23,35.41,41.75,49.71,60.43};//radius
        float stdCyl[11]={0.2065,0.2,0.1856,0.1819,0.2647,0.2597,0.2747,0.2657,1.606,1.599,1.687};//stdvar
        float disk[22]={32.35,39.41,48.96,77.9,80.71,90.4,94.02,102.7,107.0,131.6,129.4,145.5,142.8,160.1,157.1,174.2,172.7,188.4,186.3,203.2,203.8,222.2};//z
        float stddisk[22]={1.038,1.211,1.217,1.915,1.613,2.209,1.468,2.131,1.465,3.851,3.585,3.557,3.33,3.519,3.42,3.38,3.508,3.501,3.563,0.618,3.545,0};//stdz
        std::vector<std::pair<uint16_t,float> > Layer;
        std::vector<std::pair<uint16_t,float> > Disk;
};


