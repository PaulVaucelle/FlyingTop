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
        Layer.push_back(make_pair(1184,16.));  //PIXBL4
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
        Disk.push_back(make_pair(1544,77.9)); //TIDWHeel1
        Disk.push_back(make_pair(1548,80.71));//TIDWHeel1stereo
        Disk.push_back(make_pair(1552,90.4)); //TIDWHeel2
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

        //Rminus - Rplus - abs(z) 
        // PXB
        //$$$$
//  *         DetailedLayer.push_back({2.26,3.5,27});// !!!!!! WARNING : Changes due to data/MC agreement being not that rgeat in the barrel region
        //  *         // !!! Putting an large hypothetical layer that gathers beam pipe, pxl first layer, and same for pxl layers , more patatoide things than
          //  *           // !! a really precise map of the tracker. This is due to th algignmeent of the tracker being different between data/MC
        //         DetailedLayer.push_back({2.74,2.855,27}); // 0.64
        //         DetailedLayer.push_back({3.055,3.165,27});
        //         DetailedLayer.push_back({3.275,3.375,27});
//$$$$
// corrected for PXB according to the observed reco layers in 2017-2018 data and MC
// see /ui2_data1/blochd/LLTopAna/SecIntAna.C and output/h_SecInt_2017_data.root, ...
        DetailedLayer.push_back({2.71,2.88,27}); // 0.64
        DetailedLayer.push_back({3.02,3.17,27});
        DetailedLayer.push_back({3.26,3.38,27});
//$$$$         DetailedLayer.push_back({6.4,7.3,27}); // !!! large layer
        //         DetailedLayer.push_back({6.578,6.627,27}); // 0.67
        //         DetailedLayer.push_back({6.935,6.982,27});
        //         DetailedLayer.push_back({7.205,7.25,27});
        DetailedLayer.push_back({6.57,6.64,27}); // 0.67
        DetailedLayer.push_back({6.92,7.01,27});
        DetailedLayer.push_back({7.19,7.26,27});
//$$$$         DetailedLayer.push_back({10.5,11.3,27}); // !!! large layer
        //         DetailedLayer.push_back({10.707,10.738,27}); // 0.63
        //         DetailedLayer.push_back({11.037,11.066,27});
        //         DetailedLayer.push_back({11.307,11.336,27});
        DetailedLayer.push_back({10.70,10.77,27}); // 0.63
        DetailedLayer.push_back({11.02,11.11,27});
        DetailedLayer.push_back({11.29,11.34,27});
//$$$$         DetailedLayer.push_back({15.6,16.40,27}); // !!! large layer
        //         DetailedLayer.push_back({15.815,15.836,27}); // 0.62
        //         DetailedLayer.push_back({16.146,16.166,27});
        //         DetailedLayer.push_back({16.416,16.436,27});
        DetailedLayer.push_back({15.80,15.87,27}); // 0.62
        DetailedLayer.push_back({16.13,16.22,27});
        DetailedLayer.push_back({16.40,16.47,27});
//$$$$

        //TIB
        DetailedLayer.push_back({23.5,24.45,66}); // 1.0
        DetailedLayer.push_back({26.7,27.85,66}); // 1.2
        DetailedLayer.push_back({31.85,32.8,66}); // 1.1
        DetailedLayer.push_back({35,36.2,66});    // 1.2
        DetailedLayer.push_back({40,40.86,66});   // 4.0
        DetailedLayer.push_back({43.15,44.05,66});
        DetailedLayer.push_back({47.85,48.78,66});// 4.0
        DetailedLayer.push_back({51.00,51.95,66}); 

        //TOB
        DetailedLayer.push_back({58.44,58.63,107});  // 4.5
        DetailedLayer.push_back({59.525,59.71,107});
        DetailedLayer.push_back({61.64,61.82,107});
        DetailedLayer.push_back({62.725,62.9,107});

        //abs(Zminus) - abs(Zplus) - Rminus - Rplus  .... 
	//PXF
        DetailedDisk.push_back({29.6,32.2,9.6,16.1}); // 5.4
        DetailedDisk.push_back({31,35,4.5,11});
        DetailedDisk.push_back({37,39.5,9.6,16.1});   // 5.5
        DetailedDisk.push_back({38.6,42.5,4.5,11});
        DetailedDisk.push_back({46.6,49.4,9.6,16.1}); // 5.4
        DetailedDisk.push_back({48,52,4.5,11});

	//TID
        DetailedDisk.push_back({74.38,76.82,38.5,50.2}); // 6.2
        DetailedDisk.push_back({77.65,80.55,22.5,34.5});
        DetailedDisk.push_back({81.18,83.78,31.5,42});   // 4.0
        DetailedDisk.push_back({87.3,89.8,38.5,50});     // 6.2
        DetailedDisk.push_back({90.6,93.5,22.5,35});
        DetailedDisk.push_back({94.1,97,31.5,42});       // 3.2
        DetailedDisk.push_back({100.27,102.77,38.5,50.5}); // 5.9
        DetailedDisk.push_back({103.57,106.18,22.5,35});
        DetailedDisk.push_back({107.6,109.68,31.5,42});  // 2.1

	//TEC
        DetailedDisk.push_back({126.42,127.45,60,70}); // 5.7
        DetailedDisk.push_back({126.42,127.45,22,32});
        DetailedDisk.push_back({127.15,127.7,40,50});
        DetailedDisk.push_back({129.65,130.2,50,62});
        DetailedDisk.push_back({129.87,130.57,32,41});
        DetailedDisk.push_back({133.6,134.6,60,68});
        DetailedDisk.push_back({136.8,137.4,50,62});
        DetailedDisk.push_back({134.3,134.7,40,50});
        DetailedDisk.push_back({137,137.7,33,40});
        DetailedDisk.push_back({134.3,134.6,22,32});

        DetailedDisk.push_back({140.42,141.45,60,70}); // 5.7
        DetailedDisk.push_back({140.42,141.45,22,32});
        DetailedDisk.push_back({141.15,141.7,40,50});
        DetailedDisk.push_back({143.65,144.2,50,62});
        DetailedDisk.push_back({143.87,144.57,32,41});
        DetailedDisk.push_back({147.6,148.6,60,68});
        DetailedDisk.push_back({150.8,151.4,50,62});
        DetailedDisk.push_back({148.3,148.7,40,50});
        DetailedDisk.push_back({151,151.7,33,40});
        DetailedDisk.push_back({148.3,148.6,22,32});

        DetailedDisk.push_back({154.42,155.45,60,70}); // 5.7
        DetailedDisk.push_back({154.42,155.45,22,32});
        DetailedDisk.push_back({155.15,155.7,40,50});
        DetailedDisk.push_back({157.65,158.2,50,62});
        DetailedDisk.push_back({157.87,158.57,32,41});
        DetailedDisk.push_back({161.6,162.6,60,68});
        DetailedDisk.push_back({164.8,165.4,50,62});
        DetailedDisk.push_back({162.3,162.7,40,50});
        DetailedDisk.push_back({165,165.7,33,40});
        DetailedDisk.push_back({162.3,162.6,22,32});

        DetailedDisk.push_back({168.42,169.45,60,70}); // 5.7
        DetailedDisk.push_back({169.15,169.7,40,50});
        DetailedDisk.push_back({171.65,172.2,50,62});
        DetailedDisk.push_back({171.87,172.57,32,41});
        DetailedDisk.push_back({175.6,176.6,60,68});
        DetailedDisk.push_back({178.8,179.4,50,62});
        DetailedDisk.push_back({176.3,176.7,40,50});
        DetailedDisk.push_back({179,179.7,33,40});

        DetailedDisk.push_back({182.42,183.45,60,70}); // 5.7
        DetailedDisk.push_back({183.15,183.7,40,50});
        DetailedDisk.push_back({185.65,186.2,50,62});
        DetailedDisk.push_back({185.87,186.57,32,41});
        DetailedDisk.push_back({189.6,190.6,60,68});
        DetailedDisk.push_back({192.8,193.4,50,62});
        DetailedDisk.push_back({190.3,190.7,40,50});
        DetailedDisk.push_back({193,193.7,33,40});

        DetailedDisk.push_back({200.42,201.45,60,70});
        DetailedDisk.push_back({201.15,201.7,40,50});
        DetailedDisk.push_back({203.65,204.2,50,62});
        DetailedDisk.push_back({203.87,204.57,32,41});
        DetailedDisk.push_back({207.6,208.6,60,68});
        DetailedDisk.push_back({210.8,211.4,50,62});
        DetailedDisk.push_back({208.3,208.7,40,50});
        DetailedDisk.push_back({211,211.7,33,40});
      }

      // Needed Initialisation of the Tracker DAtaBase. The radius and the uncertainties are computed from the first hit of the tracks in RECO dataTier.
      // The standard deviation is also available 
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
            //This return never happens, just needed for compilation
            return GloballyPositioned<float>::PositionType (-1000.,-1000.,-1000.);
        }

      //Cylinder
      GloballyPositioned<float>::PositionType  PropagateToCylinder(uint16_t firsthit, AnalyticalPropagator* Prop, TrajectoryStateOnSurface tsos,const float eta, const float phi,const float vz) 
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

      // The two following methods are using an approximation of the detector => the mean radius of each hitpattern assuming the real thickness of a layer + the resolution effects

      int VertexBelongsToBarrelLayer(float VTX_r, float VTX_z)
        {
          float* cyl = this->CylDB();
/*           float* stdcyl = this->stdCylDB();
 */
          float* CylThick = this->Cyl_ThickDB();
          float* stdcyl_Exp = this->stdCylDB_Exp();
          for (int j=0;j<11;j++)
            {
/*               float thickness = TMath::Sqrt(stdcyl[j]*stdcyl[j]+stdcyl_Exp[j]*stdcyl_Exp[j]);
 */ 
              float thickness = TMath::Sqrt(CylThick[j]*CylThick[j]/4+stdcyl_Exp[j]*stdcyl_Exp[j]); // CylThick[j]/2 to the square
              float down = cyl[j] - thickness;
              float up   = cyl[j] + thickness;
              if(VTX_r>=down && VTX_r<=up && abs(VTX_z)<=27 && VTX_r<20){return this->Layer[j].first;}//PXB
              if(VTX_r>=down && VTX_r<=up && abs(VTX_z)<=68 && VTX_r<55 && VTX_r>20){return this->Layer[j].first;}//TIB
              if(VTX_r>=down && VTX_r<=up && abs(VTX_z)<=108 && VTX_r>57){return this->Layer[j].first;}//TOB L1
              else if (j==10 && (VTX_r<down || VTX_r>up)){return 0;}
            }
            return 0;
        }

      int VertexBelongsToDiskLayer(float VTX_r, float VTX_z)
        {
          float* disk = this->diskDB();
/*           float* stdisk = this->stddiskDB();
 */
          float* DiskThick = this->disk_ThickDB();
          float* stdisk_Exp = this->stddiskDB_Exp();
          for (int j=0;j<22;j++)
            {
/*               float thickness = TMath::Sqrt(stdisk[j]*stdisk[j]+stdisk_Exp[j]*stdisk_Exp[j]);
 */
              float thickness = TMath::Sqrt(DiskThick[j]*DiskThick[j]/4+stdisk_Exp[j]*stdisk_Exp[j]);//DiskThick[j]/2 to the square
              float down = disk[j] - thickness;
              float up   = disk[j] + thickness;
              float result=this->Disk[j].first;
              result = -result;
                  if(VTX_z>=down && VTX_z<=up && VTX_r>=4.5 && VTX_r<= 16.1 && abs(VTX_z)<=55 && abs(VTX_z)>=20){return this->Disk[j].first;}//PXF
                  else if (-VTX_z<=-down && -VTX_z>=-up && VTX_r>=4.5 && VTX_r<=16.1&& abs(VTX_z)<=55 && abs(VTX_z)>=20){return result;}
                  else if(VTX_z>=down && VTX_z<=up && VTX_r>=22 && VTX_r<=52 && abs(VTX_z)<=110 && abs(VTX_z)>=70){return this->Disk[j].first;}//TID
                  else if (-VTX_z<=-down && -VTX_z>=-up && VTX_r>=22 && VTX_r<=52 && abs(VTX_z)<=110 && abs(VTX_z)>=70){return result;}
                  else if(VTX_z>=down && VTX_z<=up && VTX_r>=22 && VTX_r<=73 && abs(VTX_z)<=200 && abs(VTX_z)>=120){return this->Disk[j].first;}//TEC
                  else if (-VTX_z<=-down && -VTX_z>=-up && VTX_r>=22 && VTX_r<=73 && abs(VTX_z)<=200 && abs(VTX_z)>=120){return result;}
                  else if (j==10 && (VTX_z<down || VTX_z>up)){return 0;}
                  else if (j==10 && (-VTX_z>-down || -VTX_z<-up)){return 0;} 
            }
            return 0;
        }

//NEW 20/03/2023
      //The following method is a more detailed mapping of the cms tracker (not complete atm) => until TOB L1 and TEC Wheel5
      int VertexBelongsToTracker(float VTX_r, float VTX_z)
        {
          std::vector<std::vector<float> > disk = this->diskDDB();
          float* stdisk_Exp = this->stddiskDDB_Exp();
          std::vector<std::vector<float> > cyl = this->CylDDB();
          float* stdcyl_Exp = this->stdCylDDB_Exp();
          int IsInTracker = 0;
          
          for (int j=0; j<24; j++)
                      {
              float ModuleThickness_2 = (cyl[j][1]-cyl[j][0])/2; // (rmax-rmin)/2
              float thickness = TMath::Sqrt(ModuleThickness_2*ModuleThickness_2+stdcyl_Exp[j]*stdcyl_Exp[j]/1.64/1.64); // thickness/2 + 1 sigma resolution
              float radius = cyl[j][1]-ModuleThickness_2; // radius of the considered module
              float down = radius - thickness;
              float up   = radius + thickness;
              if (VTX_r >= down && VTX_r <= up && abs(VTX_z) <= cyl[j][2] ) {IsInTracker=1; return 1;} // && VTX_r<20
              // else if (j==22 && (VTX_r<down || VTX_r>up)){IsInTracker = 0;}
            }
          for (int j=0; j<69; j++)
            {
              float ModuleThickness_2 = (disk[j][1]-disk[j][0])/2; // (zmax-zmin)/2
              float thickness = TMath::Sqrt(ModuleThickness_2*ModuleThickness_2+stdisk_Exp[j]*stdisk_Exp[j]/1.64/1.64); // thickness/2 + 1 sigma resolution
              float z = disk[j][1]-ModuleThickness_2; // z of the considered module
              float down = z - thickness;
              float up   = z + thickness;
              if (abs(VTX_z) >= down && abs(VTX_z) <= up && VTX_r >= disk[j][2] && VTX_r <= disk[j][3] ) {IsInTracker=1; return 1;}//PXF
            }
          if ( IsInTracker == 0 ) return 0;
          else return 1;
        }

      //-----Access Data Members------//
      float* CylDB(){return Cyl;}
      float* stdCylDB(){return stdCyl;}
      float* diskDB(){return disk;}
      float* stddiskDB(){return stddisk;}
      float* stdCylDB_Exp(){return stdCyl_Exp;}
      float* stddiskDB_Exp(){return stddisk_Exp;}
      float* stdCylDDB_Exp(){return stdDCyl_Exp;}
      float* stddiskDDB_Exp(){return stdDdisk_Exp;}
      float* Cyl_ThickDB(){return Cyl_thick; }
      float* disk_ThickDB(){return Disk_thick; }

      std::vector<std::vector<float> > CylDDB(){return DetailedLayer; }//Detailed DataBase :)
      std::vector<std::vector<float> > diskDDB(){return DetailedDisk; }


   private:
      // ----------member data ---------------------------
        float Cyl[11]={2.959,6.778,10.89,16,23.83,27.02,32.23,35.41,41.75,49.71,60.43};//Mean radius for a given hitpattern: 
/*         float stdCyl[11]={0.2065,0.2,0.1856,0.1819,0.2647,0.2597,0.2747,0.2657,1.606,1.599,1.687};//stdvar
*/
        float stdCyl[11]={0.27,0.22,0.22,0.22,0.2647,0.2597,0.2747,0.2657,1.606,1.599,1.687};//the associated std dev to the radius != real thickness of the module
        float disk[22]={32.35,39.41,48.96,77.9,80.71,90.4,94.02,102.7,107.0,131.6,129.4,145.5,142.8,160.1,157.1,174.2,172.7,188.4,186.3,203.2,203.8,222.2};////Mean z for a given hitpattern: 
        float stddisk[22]={1.038,1.211,1.217,1.915,1.613,2.209,1.468,2.131,1.465,3.851,3.585,3.557,3.33,3.519,3.42,3.38,3.508,3.501,3.563,0.618,3.545,0};//the associated std dev to the z != real thickness of the module 
        // There is more ambiguity with stddisk since TEC module in a given Wheel are quite "far in z" from each other for a given hitpattern

        //#########################################################//
        //          Add errors due to reolution effects            //
        //#########################################################// 
        float stdCyl_Exp[11]={0.12,0.12,0.12,0.12,0.7,0.7,0.7,0.7,0.7,0.7,3.};
        float stddisk_Exp[22]={1.0,1.0,1.0,2.5,2.5,2.5,2.5,2.5,2.5,4.,4.,4.,4.,4.,4.,4.,4.,4.,4.,0.,0.,0.};

/*         // compute the experimental resolution as follows : 
	// use the LLP ct=50 cm signal MC, and for each tracker subdetector region consider the track pairs which both come from the same LLP,
	// plot the delta r (or z) between the reconstructed 2-LLP tracks SV and the nearest LLP decay point
	// then get the delta r (or z) cut which removes 90% of the events
	// (corresponding to a 1.64 sigma cut if gaussian).
	// Using the Ntuple, the selections are:
	// in PXB:
	// ttree->Draw("tree_SecInt_LLP_dr","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>2.6&&tree_SecInt_r<20&&abs(tree_SecInt_z)<27&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dr>0.12")
	// in TIB:
	// ttree->Draw("tree_SecInt_LLP_dr","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>22&&tree_SecInt_r<52&&abs(tree_SecInt_z)<68&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dr>0.8")
	// in TOB:
	// ttree->Draw("tree_SecInt_LLP_dr","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>50&&tree_SecInt_r<60&&abs(tree_SecInt_z)<108&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dr>3")
	// in PXF:
	// ttree->Draw("tree_SecInt_LLP_dz","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>4.5&&tree_SecInt_r<16.5&&abs(tree_SecInt_z)>30&&abs(tree_SecInt_z)<50&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dz>1")
	// in TID:
	// ttree->Draw("tree_SecInt_LLP_dz","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>22&&tree_SecInt_r<52&&abs(tree_SecInt_z)>70&&abs(tree_SecInt_z)<110&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dz>2.5")
	// in TEC:
	// ttree->Draw("tree_SecInt_LLP_dz","tree_SecInt_LLP>0&&tree_SecInt_selec&&tree_SecInt_r>20&&tree_SecInt_r<75&&abs(tree_SecInt_z)>120&&abs(tree_SecInt_z)<200&&tree_SecInt_LLP_dr/tree_SecInt_r<0.1&&tree_SecInt_LLP_dz>5")
        //
	// if dca < 0.1 (divided by 1.64):
        float stdDCyl_Exp[24]={0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.059,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.24,0.5,0.5,0.5,0.5}; // detailed
        float stdDdisk_Exp[69]={0.40,0.40,0.40,0.40,0.40,0.40,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.}; // detailed
	// if dca < 1 (NOT divided by 1.64):
  */        
        float stdDCyl_Exp[24]={0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.12,0.8,0.8,0.8,0.8,0.8,0.8,0.8,0.8,3.,3.,3.,3.}; // detailed
        float stdDdisk_Exp[69]={1.0,1.0,1.0,1.0,1.0,1.0,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,2.5,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.,5.}; // detailed

        std::vector<std::pair<uint16_t,float> > Layer;
        std::vector<std::pair<uint16_t,float> > Disk;

        std::vector<std::vector<float> > DetailedLayer;
        std::vector<std::vector<float> > DetailedDisk;
        

        //#########################################################//
        //                    Module thickness                     //
        //#########################################################// 
        // != from the std dev
        // Computed from MC reco ntuples using the hitpattern

        float Cyl_thick[11]={0.11,0.05,0.03,0.021,0.95,0.91,1.,1.,0.95,0.97,0.2};
        float Disk_thick[22]={3,3,3,2.4,3.5,2.4,3.5,2.4,3.5,0.5,0.7,0.5,0.7,0.5,0.7,0.5,0.7,0.5,0.7,0.5,0.7,0.5};
	//--
};
