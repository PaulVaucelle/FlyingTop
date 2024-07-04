{
 gROOT->Reset() ; 

 // Compile user's analysis class //
   gROOT->ProcessLine(".L HistogramManager.C+g") ;
   gROOT->ProcessLine(".L FirstHitRes.C+g") ;
  
 if (gROOT->GetClass("FirstHitRes")==0) return;
 
 TChain c("trackingPerf/ttree");

      c.Add("UDD_bgctau50_smu275_snu225_hpsansalgontrk10.root");
      FirstHitRes* t = new FirstHitRes(&c);
      t->Loop();
}