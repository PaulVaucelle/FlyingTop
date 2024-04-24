{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L MiniEffi.C+g") ;
  
 if (gROOT->GetClass("MiniTree")==0) return;
 
 TChain c("T");

////////////////////////////////////////////////////////////////////////////////

c.Add("./24_03_2024/MiniRPV_2018_smu500_neu350_ctau100.root"); // BDTtk trained on all ctau100 samples and only TTTo2L2Nu: /opt/sbg/cms/ui2_data1/pvaucell/CMSSW_10_6_30_FLY/src/FlyingTop/FlyingTop/test/TestSample10cm/

////////////////////////////////////////////////////////////////////////////////

 MiniTree* t = new MiniTree(&c);

t->Loop(0,180.,20.,2999.,0.0,10.,0.,1,"./h_RPV_2018_smu200to500_ctau001_240422.root","");

////////////////////////////////////////////////////////////////////////////////

}