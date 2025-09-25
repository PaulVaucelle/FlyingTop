#include <iostream>
#include <TROOT.h>
#include "TH1.h"
#include "TColor.h"
// #include "../MCWeights.h"
#include "/opt/sbg/cms/ui2_data1/pvaucell/CMSSW_14_0_20/src/FlyingTop/FlyingTop/test/PlotCMS.h"

TCanvas * plot(int method, TString Prod, TString Name, TString Year, TString Dmode,  TString Plots, bool Data, bool Mc, bool SIGNAL)
{
int stati=0;
bool fit= 1;
bool logy=0;

bool DATA = Data;
bool MC = Mc;
bool Signal = SIGNAL;
  float SumDataEvent = 0;
  float SumMCEvent = 0;
  float hmin = 0.5; // cannot be 0 for logy=1
  float hmax = 1E6;	         // for eta<2.4 pt>80 or CRlowlowpt
  float hmaxBD = 1E6;	         // for eta<2.4 pt>80 or CRlowlowpt
    float ReScaleXS = 1.;
 TString EXTRA = "HT100_";
TString HTcut = "100";
TString Yearcor = Year;
float scaleMC = 0.948;
TString ProdMC = "MC_EMU_2018_03_02_2025";
TString suffixDATA = "_Corr";//
TString suffixMC = "_BDT100";
TString SampleDATA = "MuonEG-Run2018-UL2018_MiniAODv2_GT36-v1";
TString SampleDATAExtra  = "MuonEG-Run2018-UL2018_MiniAODv2_GT36-v1"+suffixDATA;
TString ProdSignal = "RPV_2018";
 if (Year == "2022A") 
  {
    Yearcor = "2022A";
    scaleMC = 0.9142*0.3966*1;
    Prod = "DATA_EMU_2022_CDE_23_04_2025";
    ProdMC = "MC_EMU_2022_23_04_2025";
    SampleDATA = "MuonEG_Run2022-CDE-22Sep2023";
    SampleDATAExtra = "MuonEG_Run2022-CDE-22Sep2023"+suffixDATA;
    ProdSignal = "RPV_2022A";
  }
  
 if (Year == "2022B") 
  {
    Yearcor = "2022B";
    scaleMC = 0.9142*0.6004*1;
    Prod = "DATA_EMU_2022_FG_23_04_2025";
    ProdMC = "MC_EMU_2022_EFG_23_04_2025";
    SampleDATA = "MuonEG_Run2022-FG-22Sep2023";
    SampleDATAExtra = "MuonEG_Run2022-FG-22Sep2023"+suffixDATA;
    ProdSignal = "RPV_2022B";
  }
 if (Year == "2023A") 
  {
    Yearcor = "2023A";
    scaleMC = 0.9256*0.653*0.81875;
    Prod = "DATA_EMU_2023_C_23_04_2025";
    ProdMC = "MC_EMU_2023_C_23_04_2025";
    SampleDATA = "MuonEG_Run2023C-22Sep2023";
    SampleDATAExtra = "MuonEG_Run2023C-22Sep2023"+suffixDATA;
    ProdSignal = "RPV_2023A";
  }
 if (Year == "2023B") 
  {
    Yearcor = "2023B";
    scaleMC =  0.9256*0.347*1;
    Prod = "DATA_EMU_2023_D_23_04_2025";
    ProdMC = "MC_EMU_2023_D_23_04_2025";
    SampleDATA = "MuonEG_Run2023D-22Sep2023";
    SampleDATAExtra = "MuonEG_Run2023D-22Sep2023"+suffixDATA;
    ProdSignal = "RPV_2023B";
  }
 if (Year == "2024") 
  {
    Yearcor = "2024";
    scaleMC = 0.99225;
    Prod = "DATA_EMU_2024_23_04_2025";
    ProdMC = "MC_EMU_2024_23_04_2025";
    ProdSignal = "RPV_2024";
    SampleDATA = "MuonEG_Run2024";
    SampleDATAExtra = "MuonEG_Run2024"+suffixDATA;
  }
  


 // % of jobs that did not fail for emu data
// Dmu
TFile* f1_Data_emu  = new TFile("../../"+Prod+"/histofile_"+EXTRA+Dmode+"_OS_2p4_"+SampleDATAExtra+".root");//_BDT001, _BDT,010, _BDT100, BDT300
// TFile* f1_Data_emu  = new TFile("../../DATA_EMU_2018_M1_v2/histofile_"+EXTRA+Dmode+"_OS_2p4_MuonEG_Run2018-UL2018_MiniAODv2_GT36-v1_BDT100.root");//_BDT001, _BDT,010, _BDT100, BDT300

//emu


 TFile* f1_DY  = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.root");
 TFile* f2_DY  = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.root");

if (Year == "2024")
  {
    f1_DY  = new TFile("../../MC_MUMU_2023_D_23_04_2025/histofile_HT100_EM_OS_2p4_DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.root");
    f2_DY  = new TFile("../../MC_MUMU_2023_D_23_04_2025/histofile_HT100_EM_OS_2p4_DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.root");
  }


// if (Year == "2024")
//   {
//     f1_DY  = new TFile("../../MC_MUMU_"+YearCor+"_"+ERA_MC+Prod+"_DY/DATAMC_"+COMPARE+"DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8.root");
//     f2_DY  = new TFile("../../MC_MUMU_"+YearCor+"_"+ERA_MC+Prod+"_DY/DATAMC_"+COMPARE+"DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8.root");
//   }


// TFile* f3_DY  = new TFile("../../MC_MUMU_"+YearCor+"_"+ERA_MC+Prod+"/DATAMC_"+COMPARE+"DYto2Mu_"+CorBin+"MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8.root");
// TFile* f4_DY  = new TFile("../../MC_MUMU_"+YearCor+"_"+ERA_MC+Prod+"/DATAMC_"+COMPARE+"DYto2Mu_"+CorBin+"MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8.root");
// TFile* f5_DY  = new TFile("../../MC_MUMU_"+YearCor+"_"+ERA_MC+Prod+"/DATAMC_"+COMPARE+"DYto2Mu_"+CorBin+"MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8.root");

// if (Year == "2022A")
// {
//     f1_DY  = new TFile("../../MC_MUMU_2022_23_04_2025_DY/DATAMC_"+COMPARE+"DYto2Mu_MLL-10to50_TuneCP5_13p6TeV_powheg-pythia8.root");
//     f2_DY  = new TFile("../../MC_MUMU_2022_23_04_2025_DY/DATAMC_"+COMPARE+"DYto2Mu_MLL-50to120_TuneCP5_13p6TeV_powheg-pythia8.root");

//     f3_DY  = new TFile("../../MC_MUMU_2022_23_04_2025_DY/DATAMC_"+COMPARE+"DYto2Mu_MLL-120to200_TuneCP5_13p6TeV_powheg-pythia8.root");
//     f4_DY  = new TFile("../../MC_MUMU_2022_23_04_2025_DY/DATAMC_"+COMPARE+"DYto2Mu_MLL-200to400_TuneCP5_13p6TeV_powheg-pythia8.root");
//     f5_DY  = new TFile("../../MC_MUMU_2022_23_04_2025_DY/DATAMC_"+COMPARE+"DYto2Mu_MLL-400to800_TuneCP5_13p6TeV_powheg-pythia8.root");
// }

 TFile* f1_TT  = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8.root");
 TFile* f2_TT  = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8.root");

 TFile* f1_ST  = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8.root");
 TFile* f2_ST  = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8.root");

 TFile* f1_TTV = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8.root");
 TFile* f2_TTV = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8.root");
 TFile* f3_TTV = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8.root");

 TFile* f1_VV  = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8.root");
 TFile* f2_VV  = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8.root");
 TFile* f3_VV  = new TFile("../../"+ProdMC+"/histofile_"+EXTRA+Dmode+"_OS_2p4_ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8.root");
 

TString MSUON[3] = {"200","300","300"};
  TString MNEU[3] = {"180","180","200"};
  TString CTAU[3] = {"100","100","100"};

//signal
//  TFile* f1_LLP = new TFile("../../"+ProdSignal+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSUON[0]+"_neu"+MNEU[0]+"_ctau"+CTAU[0]+".root");
//  TFile* f2_LLP = new TFile("../../"+ProdSignal+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSUON[1]+"_neu"+MNEU[1]+"_ctau"+CTAU[1]+".root");
//  TFile* f3_LLP = new TFile("../../"+ProdSignal+"/histofile_"+EXTRA+Dmode+"_OS_2p4_RPV_"+Yearcor+"_smu"+MSUON[2]+"_neu"+MNEU[2]+"_ctau"+CTAU[2]+".root");

// !! ----------------- !! //
TString filename1 = "histofile_"+EXTRA+"DM_OS_2p4_RPV_"+Yearcor+"_Msmu-"+MSUON[2]+"_Mchi-"+MNEU[2]+"_ct-001.root";
TString filename2 = "histofile_"+EXTRA+"DM_OS_2p4_RPV_"+Yearcor+"_Msmu-"+MSUON[2]+"_Mchi-"+MNEU[2]+"_ct-003.root";
TString filename3 = "histofile_"+EXTRA+"DM_OS_2p4_RPV_"+Yearcor+"_Msmu-"+MSUON[2]+"_Mchi-"+MNEU[2]+"_ct-010.root";
if (Year == "2024") 
  {
    filename1 = "histofile_"+EXTRA+"DM_OS_2p4_RPV_"+Yearcor+"_Par-ct-001-MChi-"+MNEU[2]+"-MSmu-"+MSUON[2]+".root";
    filename2 = "histofile_"+EXTRA+"DM_OS_2p4_RPV_"+Yearcor+"_Par-ct-010-MChi-"+MNEU[2]+"-MSmu-"+MSUON[2]+".root";
    filename3 = "histofile_"+EXTRA+"DM_OS_2p4_RPV_"+Yearcor+"_Par-ct-100-MChi-"+MNEU[2]+"-MSmu-"+MSUON[2]+".root";
  }
 TFile* f1_LLP = new TFile("../../"+ProdSignal+"/"+filename1);
 TFile* f2_LLP = new TFile("../../"+ProdSignal+"/"+filename2);
 TFile* f3_LLP = new TFile("../../"+ProdSignal+"/"+filename3);


// !! ------------ !! //


 TString DATAFILE[1] = {SampleDATA+"_"
};
if (Year == "2018") DATAFILE[0] = SampleDATA+"_";


TString MCFILE[12] = {
  
    "DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_",
    "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_", 

    "TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_",
    "TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8_",

    "TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_",
    "TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_",

    "TTLL_MLL-4to50_TuneCP5_13p6TeV_amcatnlo-pythia8_",
    "TTLL_MLL-50_TuneCP5_13p6TeV_amcatnlo-pythia8_",
    "TTLNu-1Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_",

    "WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8_",
    "WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8_",
    "ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8_"


};

if (Year == "2024")
  {
    MCFILE[0] = "DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8_";
    MCFILE[1] = "DYto2Mu-4Jets_Bin-MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8_";
  }


TString  LLPFILE[3] = {
                  "RPV_"+Yearcor+"_Msmu-"+MSUON[2]+"_Mchi-"+MNEU[2]+"_ct-001_",
                  "RPV_"+Yearcor+"_Msmu-"+MSUON[2]+"_Mchi-"+MNEU[2]+"_ct-003_",
                  "RPV_"+Yearcor+"_Msmu-"+MSUON[2]+"_Mchi-"+MNEU[2]+"_ct-010_",
 };

if (Year == "2024")
  {
     LLPFILE[0] = "RPV_"+Yearcor+"_Par-ct-001-MChi-"+MNEU[2]+"-MSmu-"+MSUON[2]+"_";
     LLPFILE[1] = "RPV_"+Yearcor+"_Par-ct-003-MChi-"+MNEU[2]+"-MSmu-"+MSUON[2]+"_";
     LLPFILE[2] = "RPV_"+Yearcor+"_Par-ct-010-MChi-"+MNEU[2]+"-MSmu-"+MSUON[2]+"_";
  }

 TString ytitle = "Events"; 
 
    TString htitleA = "hData_EVT34_1Vtx_BDTvtx";
    TString htitleB = "hData_NoEVT34_1Vtx_BDTvtx";
    TString htitleC = "hData_EVT12_1Vtx_BDTvtx";
    TString htitleD = "hData_NoEVT12_1Vtx_BDTvtx";
    int nbin = 40; 
    float xmin = 0;
    float xmax =  40;
    TString HeaderA = "A";
    TString HeaderAbis = "Abis";
    TString HeaderB = "B";
    TString HeaderBbis = "Bbis";
    TString HeaderC = "C";
    TString HeaderCbis = "Cbis";
    TString HeaderD = "D";
    TString HeaderDbis = "Dbis";
    TString HeaderNVtx = "k Vtx";
    TString xtitle = "var";
    float rwTT = 1;//0.00923001376;//0.923001376;
  //-----------------------------------------------------------//
  // ABCD using Hemipt and Tight+looseWP 
  //-----------------------------------------------------------//
  int Method = method;
if (Method == 0)
  {
    htitleA = "hData_CRtightlowlowpt_TLVtx_SumtrackWeight_";//
    htitleB = "hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight_";//
    htitleC = "hData_CRtighthighpt_2Vtx_SumtrackWeight_";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_SumtrackWeight_";//

    nbin = 19; 
    xmin = 1;
    xmax = 20;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";


    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Sum of track weights at Vtx";
  }

  // !! &&
      //-----------------------------------------------------------//


    if (Method == 1)
  {
    htitleA = "hData_CRtightlowlowpt_2Vtx_NEWMVA_";//
    htitleB = "hData_CRlooselooselowlowpt_2Vtx_NEWMVA_";//
    htitleC = "hData_CRtighthighpt_2Vtx_NEWMVA_";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_NEWMVA_";//

    nbin = 10; 
    xmin = -1;
    xmax = 1;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Vtx BDT Score";
  }


    if (Method == 2)
  {
    htitleA = "hData_CRtightlowlowpt_2Vtx_NEWMVA_6Bins";//
    htitleB = "hData_CRlooselooselowlowpt_2Vtx_NEWMVA_6Bins";//
    htitleC = "hData_CRtighthighpt_2Vtx_NEWMVA_6Bins";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_NEWMVA_6Bins";//

    nbin = 6; 
    xmin = -1;
    xmax = 1;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Vtx BDT Score";
  }


    if (Method == 3)
  {
    htitleA = "hData_CRtightlowlowpt_2Vtx_NEWMVA_7Bins";//
    htitleB = "hData_CRlooselooselowlowpt_2Vtx_NEWMVA_7Bins";//
    htitleC = "hData_CRtighthighpt_2Vtx_NEWMVA_7Bins";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_NEWMVA_7Bins";//

    nbin = 7; 
    xmin = -1;
    xmax = 1;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Sum of track weights at Vtx";
  }

    if (Method == 4)
  {
    htitleA = "hData_CRtightlowlowpt_2Vtx_NEWMVA_8Bins";//
    htitleB = "hData_CRlooselooselowlowpt_2Vtx_NEWMVA_8Bins";//
    htitleC = "hData_CRtighthighpt_2Vtx_NEWMVA_8Bins";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_NEWMVA_8Bins";//

    nbin = 8; 
    xmin = -1;
    xmax = 1;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Sum of track weights at Vtx";
  }

    if (Method == 5)
  {
    htitleA = "hData_CRtightlowlowpt_TLVtx_SumtrackWeight_8Bins";//
    htitleB = "hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight_8Bins";//
    htitleC = "hData_CRtighthighpt_TLVtx_SumtrackWeight_8Bins";//
    htitleD = "hData_CRlooselooselowpt_TLVtx_SumtrackWeight_8Bins";//

    nbin = 8; 
    xmin = 1;
    xmax = 9;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Sum of track weights at Vtx";
  }


      if (Method == 6)
  {
    htitleA = "hData_CRtightlowlowpt_TLVtx_SumtrackWeight_7Bins";//
    htitleB = "hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight_7Bins";//
    htitleC = "hData_CRtighthighpt_TLVtx_SumtrackWeight_7Bins";//
    htitleD = "hData_CRlooselooselowpt_TLVtx_SumtrackWeight_7Bins";//

    nbin = 7; 
    xmin = 1;
    xmax = 8;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Sum of track weights at Vtx";
  }
      if (Method == 7)
  {
    htitleA = "hData_CRtightlowlowpt_TLVtx_SumtrackWeight_6Bins";//
    htitleB = "hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight_6Bins";//
    htitleC = "hData_CRtighthighpt_TLVtx_SumtrackWeight_6Bins";//
    htitleD = "hData_CRlooselooselowpt_TLVtx_SumtrackWeight_6Bins";//

    nbin = 6; 
    xmin = 1;
    xmax = 7;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Sum of track weights at Vtx";
  }


if (Method == 8)
  {
    htitleA = "hData_CRtightlowlowpt_2Vtx_NEWMVA_Sum";//
    htitleB = "hData_CRlooselooselowlowpt_2Vtx_NEWMVA_Sum";//
    htitleC = "hData_CRtighthighpt_2Vtx_NEWMVA_Sum";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_NEWMVA_Sum";//

    nbin = 10; 
    xmin = -2;
    xmax = 2;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Sum of BDT Scores";
  }
if (Method == 9)
  {
    htitleA = "hData_CRtightlowlowpt_2Vtx_NEWMVA_Ave";//
    htitleB = "hData_CRlooselooselowlowpt_2Vtx_NEWMVA_Ave";//
    htitleC = "hData_CRtighthighpt_2Vtx_NEWMVA_Ave";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_NEWMVA_Ave";//

    nbin = 10; 
    xmin = -1;
    xmax = 1;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Average of BDT Scores";
  }
if (Method == 10)
  {
    htitleA = "hData_CRtightlowlowpt_2Vtx_STW_Sum";//
    htitleB = "hData_CRlooselooselowlowpt_2Vtx_STW_Sum";//
    htitleC = "hData_CRtighthighpt_2Vtx_STW_Sum";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_STW_Sum";//

    nbin = 19; 
    xmin = 1;
    xmax = 20;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Sum of STW";
  }
if (Method == 11)
  {
    htitleA = "hData_CRtightlowlowpt_2Vtx_STW_Ave";//
    htitleB = "hData_CRlooselooselowlowpt_2Vtx_STW_Ave";//
    htitleC = "hData_CRtighthighpt_2Vtx_STW_Ave";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_STW_Ave";//

    nbin = 19; 
    xmin = 1;
    xmax = 20;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Average of STW";
  }
// !! New 
if (Method == 12)
  {
    htitleA = "hData_CRtightlowlowpt_TLVtx_SumtrackWeight_Closure";//
    htitleB = "hData_CRlooselooselowlowpt_TLVtx_SumtrackWeight_Closure";//
    htitleC = "hData_CRtighthighpt_TLVtx_SumtrackWeight_Closure";//
    htitleD = "hData_CRlooselooselowpt_TLVtx_SumtrackWeight_Closure";//

    nbin = 19; 
    xmin = 1;
    xmax = 20;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";



    HeaderAbis = "2 hemi., p_{T}< 60 GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > 60 GeV (<100 GeV)";
    HeaderBbis = "2 hemi., p_{T}< 60 GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > 60 GeV (<100 GeV)";

    HeaderNVtx = "2 Vertices";
    xtitle = "Sum of track weights at Vtx";
  }
if (Method == 13)
  {
    htitleA = "2Vtx_VtxVtx_EVT_MVA_tightlowlowpt_Closure";//
    htitleB = "2Vtx_VtxVtx_EVT_MVA_looselooselowlowpt_Closure";//
    htitleC = "2Vtx_VtxVtx_EVT_MVA_tighthighpt_Closure";//
    htitleD = "2Vtx_VtxVtx_EVT_MVA_looselooselowpt_Closure";//

    nbin = 20; 
    xmin = -1;
    xmax = 1;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

    HeaderAbis = "2 hemi., 30 < p_{T}< 60 GeV";
    HeaderCbis = ">= 1 hemi.,60  < p_{T} < 100 GeV ";
    HeaderBbis = "2 hemi., 30 < p_{T}< 60 GeV";
    HeaderDbis = ">= 1 hemi., 60 < p_{T} < 100 GeV ";

    HeaderNVtx = "2 Vertices";
    xtitle = "Evt BDT Score";
  }
  if (Method == 14)
  {
    htitleA = "hData_CRtightlowlowpt_2Vtx_NEWMVA_Ave_Closure";//
    htitleB = "hData_CRlooselooselowlowpt_2Vtx_NEWMVA_Ave_Closure";//
    htitleC = "hData_CRtighthighpt_2Vtx_NEWMVA_Ave_Closure";//
    htitleD = "hData_CRlooselooselowpt_2Vtx_NEWMVA_Ave_Closure";//

    nbin = 20; 
    xmin = -1;
    xmax = 1;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";

HeaderAbis = "2 hemi., 30 < p_{T}< 60 GeV";
    HeaderCbis = ">= 1 hemi.,60  < p_{T} < 100 GeV ";
    HeaderBbis = "2 hemi., 30 < p_{T}< 60 GeV";
    HeaderDbis = ">= 1 hemi., 60 < p_{T} < 100 GeV ";

    HeaderNVtx = "2 Vertices";
    xtitle = "Ave BDT Score";
  }
if (Method == 15)
  {
    htitleA = "2Vtx_VtxVtx_EVT_MVA_tightlowlowpt";//
    htitleB = "2Vtx_VtxVtx_EVT_MVA_looselooselowlowpt";//
    htitleC = "2Vtx_VtxVtx_EVT_MVA_tighthighpt";//
    htitleD = "2Vtx_VtxVtx_EVT_MVA_looselooselowpt";//

    nbin = 10; 
    xmin = -1;
    xmax = 1;
    HeaderA = ">= 1 tight";
    HeaderB = "2 loose";
    HeaderC = ">= 1 tight";
    HeaderD = "2 loose";


    HeaderAbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderCbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";
    HeaderBbis = "2 hemi., p_{T}< "+HTcut+" GeV";
    HeaderDbis = ">= 1 hemi., p_{T} > "+HTcut+" GeV";

    HeaderNVtx = "2 Vertices";
    xtitle = "Evt BDT Score";
  }



  TLegend* leg;
    
        // Couleurs pour les histogrammes
    Float_t r1 = 0.246;
Float_t g1 = 0.563;
Float_t b1 = 0.852;
TColor* color1 = new TColor(301,r1, g1, b1);
// color1->SetRGB(r1, g1, b1);
Int_t ColorBlue = color1->GetNumber();
// Int_t ColorBlue = TColor::GetColor(r1, g1, b1);

Float_t r2 = 1.000;
Float_t g2 = 0.661;
Float_t b2 = 0.055;
TColor* color2 = new TColor(302,r2, g2, b2);
// color2.SetRGB(r2, g2, b2);
Int_t ColorOrange = color2->GetNumber();

Float_t r3 = 0.739;
Float_t g3 = 0.122;
Float_t b3 = 0.004;
TColor* color3 = new TColor(303,r3, g3, b3);
Int_t ColorRed = color3->GetNumber();

Float_t r4 = 0.578;
Float_t g4 = 0.641;
Float_t b4 = 0.635;
TColor* color4 = new TColor(304,r4, g4, b4);
Int_t ColorGrey = color4->GetNumber();

Float_t r5 = 0.513;
Float_t g5 = 0.176;
Float_t b5 = 0.713;
TColor* color5 = new TColor(305,r5, g5, b5);
Int_t ColorDarkPurple = color5->GetNumber();

Float_t r6 = 0.661;
Float_t g6 = 0.418;
Float_t b6 = 0.348;
TColor* color6 = new TColor(306,r6, g6, b6);
Int_t ColorBrown = color6->GetNumber();

Float_t r7 = 0.905;
Float_t g7 = 0.387;
Float_t b7 = 0.000;
TColor* color7 = new  TColor(307,r7, g7, b7);
Int_t ColorDarkOrange = color7->GetNumber();

Float_t r8 = 0.723;
Float_t g8 = 0.672;
Float_t b8 = 0.438;
TColor* color8 = new TColor(308,r8, g8, b8);
Int_t ColorNeutral = color8->GetNumber();

Float_t r9 = 0.441;
Float_t g9 = 0.457;
Float_t b9 = 0.504;
TColor* color9 = new TColor(309,r9, g9, b9);
Int_t ColorDarkGrey = color9->GetNumber();

Float_t r10 = 0.571;
Float_t g10 = 0.852;
Float_t b10 = 0.867;
TColor* color10 = new TColor(310,r10, g10, b10);
// color10.SetRGB(r10, g10, b10);
Int_t ColorLightBlue = color10->GetNumber();
// *****************************************************************************

TCanvas *c1 = new TCanvas("c1", "plots",0,0,1300,1200);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

// !! -------------------
// !! Pads
float x1 = 0.03;
float x2 = 0.49;
float x3 = 0.51;
float x4 = 0.97;

float y1 = 0.62;
float y2 = 0.92;
float y3 = 0.17;
float y4 = 0.47;

// !! raps

float x5 = 0.03;
float x6 = 0.49;
float x7 = 0.51;
float x8 = 0.97;

float y5 = 0.48;
float y6 = 0.61;
float y7 = 0.03;
float y8 = 0.16;
// !! ------------------
if ((DATA && !MC ) || (!DATA && MC))
{

  y1 = 0.51;
  y2 = 0.92;
  y3 = 0.03;
  y4 = 0.44;

  y5 = 0.0;
  y6 = 0.005;
  y7 = 0.0;
  y8 = 0.0005;


}


  TPad* pad1 = new TPad("pad1","This is pad1",x1,y1,x2,y2,21);
TPad* pad2 = new TPad("pad2","This is pad2",x3,y1,x4,y2,21);
TPad* pad3 = new TPad("pad3","This is pad3",x1,y3,x2,y4,21);
TPad* pad4 = new TPad("pad4","This is pad4",x3,y3,x4,y4,21); 

  TPad* pad8 = new TPad("pad8","This is pad8",x5,y5,x6,y6,21);
  TPad* pad9 = new TPad("pad9","This is pad9",x7,y5,x8,y6,21);
  TPad* pad10 = new TPad("pad10","This is pad10",x5,y7,x6,y8,21);
  TPad* pad11 = new TPad("pad11","This is pad11",x7,y7,x8,y8,21);

float bottomMargin = 0.05;
if ((DATA && !MC) || ( !DATA && MC )) bottomMargin = 0.15;
pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogy(logy);
   pad1->SetTopMargin(0.07);
   pad1->SetBottomMargin(bottomMargin);
   pad1->SetRightMargin(0.04);
   pad1->SetLeftMargin(0.16);

pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogy(logy);
   pad2->SetTopMargin(0.07);
   pad2->SetBottomMargin(bottomMargin);
   pad2->SetRightMargin(0.04);
   pad2->SetLeftMargin(0.16);

pad3->SetFillColor(0);
pad3->SetBorderMode(0);
pad3->SetFrameFillColor(10);
pad3->Draw();
pad3->SetLogy(logy);
   pad3->SetTopMargin(0.07);
   pad3->SetBottomMargin(bottomMargin);
   pad3->SetRightMargin(0.04);
   pad3->SetLeftMargin(0.16);

pad4->SetFillColor(0);
pad4->SetBorderMode(0);
pad4->SetFrameFillColor(10);
pad4->Draw();
pad4->SetLogy(logy);
   pad4->SetTopMargin(0.07);
   pad4->SetBottomMargin(bottomMargin);
   pad4->SetRightMargin(0.04);
   pad4->SetLeftMargin(0.16);

pad8->SetFillColor(0);
pad8->SetBorderMode(0);
pad8->SetFrameFillColor(10);
pad8->Draw();
pad8->SetLogy(0);
   pad8->SetTopMargin(0.07);
   pad8->SetBottomMargin(0.35);
   pad8->SetRightMargin(0.04);
   pad8->SetLeftMargin(0.16);

pad9->SetFillColor(0);
pad9->SetBorderMode(0);
pad9->SetFrameFillColor(10);
pad9->Draw();
pad9->SetLogy(0);
   pad9->SetTopMargin(0.07);
   pad9->SetBottomMargin(0.35);
   pad9->SetRightMargin(0.04);
   pad9->SetLeftMargin(0.16);

pad10->SetFillColor(0);
pad10->SetBorderMode(0);
pad10->SetFrameFillColor(10);
pad10->Draw();
pad10->SetLogy(0);
   pad10->SetTopMargin(0.07);
   pad10->SetBottomMargin(0.35);
   pad10->SetRightMargin(0.04);
   pad10->SetLeftMargin(0.16);

pad11->SetFillColor(0);
pad11->SetBorderMode(0);
pad11->SetFrameFillColor(10);
pad11->Draw();
pad11->SetLogy(0);
   pad11->SetTopMargin(0.07);
   pad11->SetBottomMargin(0.35);
   pad11->SetRightMargin(0.04);
   pad11->SetLeftMargin(0.16);

gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleFont(62);
gStyle->SetTitleColor(1);
gStyle->SetTitleTextColor(1);
gStyle->SetTitleFillColor(10);
gStyle->SetTitleFontSize(0.05);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(0.09);
gStyle->SetOptStat(stati);
gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);
gROOT->SetBatch(kTRUE);
 TH1F* hsolve  = new TH1F("hsolve","",nbin,xmin,xmax);
 hsolve->Sumw2();

 TH1F* htotMC  = new TH1F("htotMC","",nbin,xmin,xmax);
TH1F* htotMCv2  = new TH1F("htotMCv2","",nbin,xmin,xmax);
 TH1F* htotData  = new TH1F("htotData","",nbin,xmin,xmax);
TH1F* htotDatav2  = new TH1F("htotData","",nbin,xmin,xmax);
TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 hsolve->Sumw2();
  htotData->Sumw2();
htotDatav2->Sumw2();
  htotMC->Sumw2();
htotMCv2->Sumw2();

TH1F* g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleA);

 TH1F* g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);//ok
 TH1F* g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);//ok
 TH1F*  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);


 TH1F* g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);//ok
 TH1F* g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);//ok
 TH1F* g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);//ok
 TH1F*  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);

TH1F* g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);//ok
 TH1F* g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);//ok
 TH1F*  h_ST = new TH1F("h_ST","",nbin,xmin,xmax);

 TH1F* g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);//ok
TH1F* g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);//ok
//   TH1F* g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);//ok
//  TH1F* g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);//ok
 TH1F*  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);

 TH1F* g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);//ok
 TH1F* g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);//ok
 TH1F* g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);//ok
 TH1F*  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);

 TH1F* g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleA);//ok
 TH1F* h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 TH1F* g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleA);//ok
 TH1F* h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 TH1F* g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleA);//ok
 TH1F* h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);

// *****************************************************************************

 pad1->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

 

if (MC)
  {
 f1_DY->cd();
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleA);
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleA);
  h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
  g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleA);
  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleA);
  h_VV->Add(g2_VV, h_VV,1, 1);

 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleA);
  h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleA);
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
  g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleA);
  h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleA);
  h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleA);
 
 h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleA);
  h_ST->Add(g2_ST, h_ST, 1, 1);

    //   f3_ST->cd();
    // g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleA);
    // h_ST->Add(g3_ST, h_ST, 1, 1);

    //   f4_ST->cd();
    // g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleA);
    // h_ST->Add(g4_ST, h_ST, 1, 1);

 f1_TT->cd();
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleA);
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, rwTT*1,0);

    f2_TT->cd();
    g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleA);
    h_TT->Add(g2_TT, h_TT, 1,1);

 //!! --------------------------
h_DY->Scale(1*scaleMC);
h_VV->Scale(1*scaleMC);
h_TTV->Scale(1*scaleMC);
h_ST->Scale(1*scaleMC);
h_TT->Scale(1*scaleMC);

// Daniel 
 h_VV->Add(h_DY, h_VV, 1, 1);
 h_ST->Add(h_ST, h_VV, 1, 1);
 h_TTV->Add(h_TTV, h_ST, 1, 1);
 htotMC->Add(h_TTV, htotMC, 1, 0);
 htotMC->Add(htotMC,h_TT, 1, 1);

if (MC && !DATA)
      {
        htotMC->Draw("PE1");
        htotMC->SetMarkerStyle(20);
        htotMC->SetMarkerSize(1);
        htotMC->SetMarkerColor(kBlack);
        htotMC->SetLineColor(kBlack);
        htotMC->SetLineWidth(1);
      }
    if (MC && DATA)
      {
 htotMC->Draw("HE"); 
  htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(ColorRed, 1);
 htotMC->SetLineColor(ColorRed);
htotMC->SetLineColorAlpha(ColorRed, 1);
}


 htotMC->SetLineStyle(1);
 htotMC->SetLineWidth(1);
 htotMC->SetTickLength(0.03, "YZ");
 htotMC->SetTickLength(0.03,"X");
 htotMC->SetLabelOffset(0.015,"X");
 htotMC->SetLabelOffset(0.007,"Y");
 htotMC->SetLabelSize(0.045, "XYZ");
 htotMC->SetLabelFont(42, "XYZ"); 
 htotMC->SetTitleSize(0.055, "XYZ"); 
 htotMC->SetTitleFont(42, "XYZ");
 htotMC->SetTitleOffset(1.2,"X"); 
 htotMC->SetTitleOffset(1.3,"Y");
 htotMC->GetXaxis()->SetTitle(xtitle);
 htotMC->GetXaxis()->SetTitleColor(1);
 htotMC->GetYaxis()->SetTitle(ytitle);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");

SumMCEvent += htotMC->Integral(0,nbin+1);

//  if (Method == 2 || Method == 3 || Method == 4 || Method == 5)
//    {
//       htotMC->GetXaxis()->SetRangeUser(0,20);
//    }


 //  htotMC->SetMinimum(hmin); 
 //  htotMC->SetMaximum(hmax); 
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);
 if (logy)
  {
    htotMC->SetMinimum(1); 
    htotMC->SetMaximum(htotMC->GetMaximum()*100); 
  }
else 
  {
    htotMC->SetMinimum(0); 
    htotMC->SetMaximum(htotMC->GetMaximum()*2); 
  }

}
if (DATA)
  {
    f1_Data_emu->cd();
    g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleA);
    htotData->Add(g1_Data_emu, htotData, 1, 0);

htotData->Draw("PE1same");
htotData->SetMarkerStyle(20);
htotData->SetMarkerSize(1);
htotData->SetMarkerColor(kBlack);
htotData->SetLineColor(kBlack);
htotData->SetLineWidth(1);

SumDataEvent += htotData->Integral(0,nbin+1);
if (DATA && !MC)
      {
        htotData->SetTickLength(0.03, "YZ");
        htotData->SetTickLength(0.03,"X");
        htotData->SetLabelOffset(0.015,"X");
        htotData->SetLabelOffset(0.007,"Y");
        htotData->SetLabelSize(0.045, "XYZ");
        htotData->SetLabelFont(42, "XYZ"); 
        htotData->SetTitleSize(0.055, "XYZ"); 
        htotData->SetTitleFont(42, "XYZ");
        htotData->SetTitleOffset(1.2,"X"); 
        htotData->SetTitleOffset(1.3,"Y");
            htotData->GetXaxis()->SetTitle(xtitle);
        htotData->GetXaxis()->SetTitleColor(1);
        htotData->GetYaxis()->SetTitle(ytitle);
        htotData->GetYaxis()->SetTitleColor(1);
        // htotData->GetXaxis()->SetRangeUser(0,20);
        if (logy)
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*100); 
          }
        else 
          {
            htotData->SetMinimum(0); 
            htotData->SetMaximum(htotData->GetMaximum()*2); 
          }
      }

  }
if (Signal)
  {
 f1_LLP->cd();
 
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleA);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

 f2_LLP->cd();
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleA);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

 f3_LLP->cd();
 
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleA);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 1,0);


 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(kRed-1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);

 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(kRed-2);
 h2_LLP->SetLineStyle(2);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(kRed-3);
 h3_LLP->SetLineStyle(3);
 h3_LLP->SetLineWidth(2);
}
 
 if (DATA )
  {
 hsolve->Add(hsolve, htotData, 0., 1.);
}
 else if (MC && !DATA)
  {
    hsolve->Add(hsolve, htotMC, 0., 1.);
  }

  // leg = new TLegend(0.17,0.94,0.50,0.98);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.06);
  // leg->SetHeader(HeaderCMS);
  // leg->Draw();

  float LEGY1 = 0.50;
  float LEGY2 = 0.89;
  float legsize = 0.06;
  float legsize2 = 0.04;
  if ( (DATA && !MC) || (!DATA && MC))
    {
      LEGY1 = 0.75;
      LEGY2 = 0.89;
    }
  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
if (DATA)
    {
  leg->AddEntry(htotData, " e#mu data","PE1");
}  
  if (MC && !DATA)
    {
  leg->AddEntry(htotMC, " Total MC","PE1");
  }
  if (MC && DATA)
    {
       leg->AddEntry(htotMC, " Total MC","F");
  }
  if (Signal)
    {
      leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
}


  leg->Draw();


  leg = new TLegend(0.23,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  
  leg = new TLegend(0.23,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderA);
  leg->Draw();

  if (Method == 2)
    {
      leg = new TLegend(0.23,0.69,0.35,0.74);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(legsize);
  leg->SetHeader(HeaderAbis);
  leg->Draw();
}
  else
    {
      leg = new TLegend(0.23,0.69,0.35,0.74);
      leg->SetBorderSize(0);
      leg->SetFillColor(kWhite);
      leg->SetTextFont(42);
      leg->SetTextSize(legsize);
      leg->SetHeader(HeaderAbis);
      leg->Draw();
    }

PlotCMSv3(pad1,Year,DATA);
  //Data/MC ---------------------
if (DATA  && MC)
    {
  pad8->cd();
    hDataMC->Divide(htotData,htotMC,1,1);
  hDataMC->Draw("PE1"); 
  hDataMC->SetFillStyle(1001);
  hDataMC->SetFillColorAlpha(kBlack, 1);
  hDataMC->SetLineColor(kBlack);
  hDataMC->SetLineStyle(1);
  hDataMC->SetLineWidth(1);
  hDataMC->SetTickLength(0.03, "YZ");
  hDataMC->SetTickLength(0.03,"X");
  hDataMC->SetLabelOffset(0.02,"X");
  hDataMC->SetLabelOffset(0.007,"Y");
  hDataMC->SetLabelSize(0.09, "XYZ");
  hDataMC->SetLabelFont(42, "XYZ"); 
  hDataMC->SetTitleSize(0.12, "XYZ"); 
  hDataMC->SetTitleFont(42, "XYZ");
  hDataMC->SetTitleOffset(1.1,"X"); 
  hDataMC->SetTitleOffset(0.5,"Y");
  hDataMC->GetXaxis()->SetTitle(xtitle);
  hDataMC->GetXaxis()->SetTitleColor(1);
  hDataMC->GetYaxis()->SetTitle("Data/MC");
  hDataMC->GetYaxis()->SetTitleColor(1);
  hDataMC->SetNdivisions(509,"XYZ");
  hDataMC->SetMinimum(0); 
  hDataMC->SetMaximum(2); 
// hDataMC->GetXaxis()->SetRangeUser(0,20);
}


// *****************************************************************************

 pad2->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);

hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 
if (MC)
  {
 f1_DY->cd();
g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleC);
 
   h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleC);
  h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleC);

  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleC);
 
  h_VV->Add(g2_VV, h_VV, 1, 1);
 f3_VV->cd();
 
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleC);
  h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleC);
 
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
  g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleC);
  h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
  g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleC);
  h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
  g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleC);
  h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
  g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleC);
  h_ST->Add(g2_ST, h_ST, 1, 1);

      //   f3_ST->cd();
      // g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleC);
      // h_ST->Add(g3_ST, h_ST, 1, 1);

      //   f4_ST->cd();
      // g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleC);
      // h_ST->Add(g4_ST, h_ST, 1, 1);


 f1_TT->cd();
 g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleC);
 
 h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, 1,0);

f2_TT->cd();
      g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleC);
      h_TT->Add(g2_TT, h_TT, 1,1);

//!! --------------------------
h_DY->Scale(1*scaleMC);
h_VV->Scale(1*scaleMC);
h_TTV->Scale(1*scaleMC);
h_ST->Scale(1*scaleMC);
h_TT->Scale(1*scaleMC);
// Daniel 
 h_VV->Add(h_DY, h_VV, 1, 1);
 h_ST->Add(h_ST, h_VV, 1, 1);
 h_TTV->Add(h_TTV, h_ST, 1, 1);
 htotMC->Add(h_TTV, htotMC, 1, 0);
 htotMC->Add(htotMC,h_TT, 1, 1);
// -----


if (MC && !DATA)
  {
// htotMCv2(htotMC);
        htotMCv2 = (TH1F*) htotMC->Clone("htotMC");  // Copie indépendante

    htotMC->Draw("PE1");
    htotMC->SetMarkerStyle(20);
    htotMC->SetMarkerSize(1);
    htotMC->SetMarkerColor(kBlack);
    htotMC->SetLineColor(kBlack);
    htotMC->SetLineWidth(1);

      }
    if (MC && DATA)
      {
        htotMCv2 = (TH1F*) htotMC->Clone("htotMC");  // Copie indépendante
        htotMC->Draw("HE"); 
          htotMC->SetFillStyle(1001);
        htotMC->SetFillColorAlpha(ColorRed, 1);
        htotMC->SetLineColor(ColorRed);
        htotMC->SetLineColorAlpha(ColorRed, 1);
              }
      if (logy)
        {
          htotMC->SetMinimum(1); 
          htotMC->SetMaximum(htotMC->GetMaximum()*100); 
        }
      else 
        {
          htotMC->SetMinimum(0); 
          htotMC->SetMaximum(htotMC->GetMaximum()*2); 
          }
 htotMC->SetLineStyle(1);
 htotMC->SetLineWidth(1);
 htotMC->SetTickLength(0.03, "YZ");
 htotMC->SetTickLength(0.03,"X");
 htotMC->SetLabelOffset(0.015,"X");
 htotMC->SetLabelOffset(0.007,"Y");
 htotMC->SetLabelSize(0.045, "XYZ");
 htotMC->SetLabelFont(42, "XYZ"); 
 htotMC->SetTitleSize(0.055, "XYZ"); 
 htotMC->SetTitleFont(42, "XYZ");
 htotMC->SetTitleOffset(1.2,"X"); 
 htotMC->SetTitleOffset(1.3,"Y");
 htotMC->GetXaxis()->SetTitle(xtitle);
 htotMC->GetXaxis()->SetTitleColor(1);
 htotMC->GetYaxis()->SetTitle(ytitle);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");

 SumMCEvent += htotMC->Integral(0,nbin+1);

 std::cout<<"EMu MC in SR "<<htotMC->Integral(0,nbin+1)<<std::endl;

//  if (Method == 2 || Method == 3 || Method == 4 || Method == 5)
//    {
//       htotMC->GetXaxis()->SetRangeUser(0,20);
//    }


}
 if (DATA)
  {
    f1_Data_emu->cd();
    g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleC);
    htotData->Add(g1_Data_emu, htotData, 1, 0);
htotDatav2 = (TH1F*) htotData->Clone("htotData");  // Copie indépendante
htotData->Draw("PE1same");
htotData->SetMarkerStyle(20);
htotData->SetMarkerSize(1);
htotData->SetMarkerColor(kBlack);
htotData->SetLineColor(kBlack);
htotData->SetLineWidth(1);
SumDataEvent += htotData->Integral(0,nbin+1);
std::cout<<"EMu data in SR "<<htotData->Integral(0,nbin+1)<<std::endl;
if (DATA && !MC)
      {
        htotData->SetTickLength(0.03, "YZ");
        htotData->SetTickLength(0.03,"X");
        htotData->SetLabelOffset(0.015,"X");
        htotData->SetLabelOffset(0.007,"Y");
        htotData->SetLabelSize(0.045, "XYZ");
        htotData->SetLabelFont(42, "XYZ"); 
        htotData->SetTitleSize(0.055, "XYZ"); 
        htotData->SetTitleFont(42, "XYZ");
        htotData->SetTitleOffset(1.2,"X"); 
        htotData->SetTitleOffset(1.3,"Y");
            htotData->GetXaxis()->SetTitle(xtitle);
        htotData->GetXaxis()->SetTitleColor(1);
        htotData->GetYaxis()->SetTitle(ytitle);
        htotData->GetYaxis()->SetTitleColor(1);
// htotData->GetXaxis()->SetRangeUser(0,20);
        if (logy)
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*100); 
          }
        else 
          {
            htotData->SetMinimum(0); 
            htotData->SetMaximum(htotData->GetMaximum()*2); 
          }
      }

  }


if(Signal)
  {
 f1_LLP->cd();
 
 
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleC);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

 f2_LLP->cd();
 
 
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleC);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

 f3_LLP->cd();
  
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleC);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(kRed-1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);

 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(kRed-2);
 h2_LLP->SetLineStyle(2);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(kRed-3);
 h3_LLP->SetLineStyle(3);
 h3_LLP->SetLineWidth(2);
}


  // leg = new TLegend(0.17,0.94,0.50,0.98);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.06);
  // leg->SetHeader(HeaderCMS);
  // leg->Draw();

  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
if(DATA)
    {
  leg->AddEntry(htotData, " e#mu data","PE1");
}
  if (MC && !DATA)
    {
  leg->AddEntry(htotMC, " Total MC","PE1");
  }
if (MC && DATA)
    {
       leg->AddEntry(htotMC, " Total MC","F");
    }
  // leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  // leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  // leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
leg->AddEntry(hsolve," Prediction","FE4");

//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();


  leg = new TLegend(0.23,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  
  leg = new TLegend(0.23,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderC);
  leg->Draw();

  
if (Method == 2 )
  {
    leg = new TLegend(0.23,0.69,0.35,0.74);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(legsize);
  leg->SetHeader(HeaderCbis);
  leg->Draw();
}
else
  {
        leg = new TLegend(0.23,0.69,0.35,0.74);
      leg->SetBorderSize(0);
      leg->SetFillColor(kWhite);
      leg->SetTextFont(42);
      leg->SetTextSize(legsize);
      leg->SetHeader(HeaderCbis);
      leg->Draw();
  }
PlotCMSv3(pad2,Year,DATA);
if(DATA && MC)
  {
  pad9->cd();
  
  // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
  hDataMC->Divide(htotData,htotMC,1,1);
  hDataMC->Draw("PE1"); 
  hDataMC->SetFillStyle(1001);
  hDataMC->SetFillColorAlpha(kBlack, 1);
  hDataMC->SetLineColor(kBlack);
  hDataMC->SetLineStyle(1);
  hDataMC->SetLineWidth(1);
  hDataMC->SetTickLength(0.03, "YZ");
  hDataMC->SetTickLength(0.03,"X");
  hDataMC->SetLabelOffset(0.015,"X");
  hDataMC->SetLabelOffset(0.007,"Y");
  hDataMC->SetLabelSize(0.09, "XYZ");
  hDataMC->SetLabelFont(42, "XYZ"); 
  hDataMC->SetTitleSize(0.12, "XYZ"); 
  hDataMC->SetTitleFont(42, "XYZ");
  hDataMC->SetTitleOffset(1.1,"X"); 
  hDataMC->SetTitleOffset(0.5,"Y");
  hDataMC->GetXaxis()->SetTitle(xtitle);
  hDataMC->GetXaxis()->SetTitleColor(1);
  hDataMC->GetYaxis()->SetTitle(ytitle);
  hDataMC->GetYaxis()->SetTitleColor(1);
  hDataMC->SetNdivisions(509,"XYZ");
  hDataMC->SetMinimum(0); 
  hDataMC->SetMaximum(2);
// hDataMC->GetXaxis()->SetRangeUser(0,20);
  
}

 
 /// !! --------------------------------------------------------------------------//
 /// !! --------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//
 ///--------------------------------------------------------------------------//


  TCanvas *c2 = new TCanvas("c2", "plots",200,0,700,700);
  c2->SetFillColor(10);
  c2->SetFillStyle(4000);
  c2->SetBorderSize(2);
  TPad* pad5 = new TPad("pad5","This is pad5",0.04,0.35,0.95,0.96,21);
  TPad* pad6 = new TPad("pad6","This is pad6",0.04,0.03,0.95,0.35,21);



  pad5->SetFillColor(0);
pad5->SetBorderMode(0);
pad5->SetFrameFillColor(10);
pad5->Draw();
pad5->SetLogy(0);
   pad5->SetTopMargin(0.07);
   pad5->SetBottomMargin(0.13);
   pad5->SetRightMargin(0.04);
   pad5->SetLeftMargin(0.16);

pad6->SetFillColor(0);
pad6->SetBorderMode(0);
pad6->SetFrameFillColor(10);
pad6->Draw();
pad6->SetLogy(0);
   pad6->SetTopMargin(0.07);
   pad6->SetBottomMargin(0.16);
   pad6->SetRightMargin(0.04);
   pad6->SetLeftMargin(0.16);

// pad7->SetFillColor(0);
// pad7->SetBorderMode(0);
// pad7->SetFrameFillColor(10);
// pad7->Draw();
// pad7->SetLogy(0);
//    pad7->SetTopMargin(0.07);
//    pad7->SetBottomMargin(0.13);
//    pad7->SetRightMargin(0.04);
//    pad7->SetLeftMargin(0.16);

/// !! --------------------------------------------------------------------------//
 /// !! --------------------------------------------------------------------------//
c2->cd();
  pad5->cd();
if (DATA)
  {
    if (DATA && !MC)
      {
          htotDatav2->Draw("E1"); 
          htotDatav2->SetFillColor(kBlack);
          htotDatav2->SetLineColor(kBlack);
          htotDatav2->SetLineStyle(1);
          htotDatav2->SetLineWidth(1);
          htotDatav2->SetMarkerStyle(20);
          htotDatav2->SetMarkerSize(1);
          htotDatav2->SetTickLength(0.03, "YZ");
          htotDatav2->SetTickLength(0.03,"X");
          htotDatav2->SetLabelOffset(0.015,"X");
          htotDatav2->SetLabelOffset(0.015,"Y");
          htotDatav2->SetLabelSize(0.045, "XYZ");
          htotDatav2->SetLabelFont(42, "XYZ"); 
          htotDatav2->SetTitleSize(0.045, "XYZ"); 
          htotDatav2->SetTitleFont(42, "XYZ");
          htotDatav2->SetTitleOffset(1.2,"X"); 
          htotDatav2->SetTitleOffset(1.3,"Y");
          htotDatav2->GetXaxis()->SetTitle(xtitle);
          htotDatav2->GetXaxis()->SetTitleColor(1);
          htotDatav2->GetYaxis()->SetTitle(ytitle);
          htotDatav2->GetYaxis()->SetTitleColor(1);
          htotDatav2->SetNdivisions(509,"XYZ");
            htotDatav2->SetMinimum(1); 
            htotDatav2->SetMaximum(htotDatav2->GetMaximum()*2); 
            std::cout<<" EMU SR  data integral :  "<<htotDatav2->Integral(0,nbin+1) << std::endl;
            //  if (Method == 2 || Method == 3 || Method == 4 || Method == 5)
            //     {
            //         htotDatav2->GetXaxis()->SetRangeUser(0,20);
            //     }
      }
    if (DATA && MC)
      {
          htotDatav2->Draw("E1"); 
          htotDatav2->SetFillColor(kBlack);
          htotDatav2->SetLineColor(kBlack);
          htotDatav2->SetLineStyle(1);
          htotDatav2->SetLineWidth(1);
          htotDatav2->SetMarkerStyle(20);
          htotDatav2->SetMarkerSize(1);
          htotDatav2->SetTickLength(0.03, "YZ");
          htotDatav2->SetTickLength(0.03,"X");
          htotDatav2->SetLabelOffset(0.015,"X");
          htotDatav2->SetLabelOffset(0.007,"Y");
          htotDatav2->SetLabelSize(0.045, "XYZ");
          htotDatav2->SetLabelFont(42, "XYZ"); 
          htotDatav2->SetTitleSize(0.045, "XYZ"); 
          htotDatav2->SetTitleFont(42, "XYZ");
          htotDatav2->SetTitleOffset(1.2,"X"); 
          htotDatav2->SetTitleOffset(1.3,"Y");
          htotDatav2->GetXaxis()->SetTitle(xtitle);
          htotDatav2->GetXaxis()->SetTitleColor(1);
          htotDatav2->GetYaxis()->SetTitle(ytitle);
          htotDatav2->GetYaxis()->SetTitleColor(1);
          htotDatav2->SetNdivisions(509,"XYZ");
            htotDatav2->SetMinimum(1); 
            htotDatav2->SetMaximum(htotDatav2->GetMaximum()*2);
                //          if (Method == 2 || Method == 3 || Method == 4 || Method == 5)
                // {
                //     htotDatav2->GetXaxis()->SetRangeUser(0,20);
                // }  
      }
  }  

if ( MC) { 
 
    if (MC && !DATA)
      {
        htotMCv2->Draw("PE1");
        htotMCv2->SetMarkerStyle(20);
        htotMCv2->SetMarkerSize(1);
        htotMCv2->SetMarkerColor(kBlack);
        htotMCv2->SetLineColor(kBlack);
        htotMCv2->SetLineWidth(1);

                   htotMCv2->SetTickLength(0.03, "YZ");
          htotMCv2->SetTickLength(0.03,"X");
          htotMCv2->SetLabelOffset(0.015,"X");
          htotMCv2->SetLabelOffset(0.007,"Y");
          htotMCv2->SetLabelSize(0.045, "XYZ");
          htotMCv2->SetLabelFont(42, "XYZ"); 
          htotMCv2->SetTitleSize(0.045, "XYZ"); 
          htotMCv2->SetTitleFont(42, "XYZ");
          htotMCv2->SetTitleOffset(1.2,"X"); 
          htotMCv2->SetTitleOffset(1.3,"Y");
          htotMCv2->GetXaxis()->SetTitle(xtitle);
          htotMCv2->GetXaxis()->SetTitleColor(1);
          htotMCv2->GetYaxis()->SetTitle(ytitle);
          htotMCv2->GetYaxis()->SetTitleColor(1);
          htotMCv2->SetNdivisions(509,"XYZ");


                    htotMCv2->SetMinimum(1); 
            htotMCv2->SetMaximum(htotMCv2->GetMaximum()*2); //htotMC->GetMaximum()*2
          
      }
    if (MC && DATA)
      {
        htotMCv2->Draw("HEsame"); 
 htotMCv2->SetFillStyle(1001);
 htotMCv2->SetFillColorAlpha(ColorRed, 1);
 htotMCv2->SetLineColor(ColorRed);
 htotMCv2->SetLineColorAlpha(ColorRed, 1);
    htotMCv2->SetMinimum(1); 
    htotMCv2->SetMaximum(htotMCv2->GetMaximum()*2); //htotMC->GetMaximum()*1
 }


              // if (Method == 2 || Method == 3 || Method == 4 || Method == 5)
              //   {
              //       htotMCv2->GetXaxis()->SetRangeUser(0,20);
              //   }
}

if (Signal)
  {
    f1_LLP->cd();
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleC);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

 f2_LLP->cd();
  g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleC);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

 f3_LLP->cd();
  g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleC);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(kRed-1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);

 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(kRed-2);
 h2_LLP->SetLineStyle(2);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(kRed-3);
 h3_LLP->SetLineStyle(3);
 h3_LLP->SetLineWidth(2);
}


  
  // leg = new TLegend(0.17,0.94,0.50,0.98);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.06);
  // leg->SetHeader(HeaderCMS);
  // leg->Draw();

 
  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
if (DATA)
    {
  leg->AddEntry(htotData, " e#mu data","PE1");
}
  if (MC && !DATA)
    {
      leg->AddEntry(htotMC, " Total MC","PE1");
      }
  if (MC  && DATA)
    {
      leg->AddEntry(htotMC, " Total MC","F");
    }
  // leg->AddEntry(htotData, " e#mu data","PE1");
leg->AddEntry(hsolve," Prediction","FE4");
//   leg->AddEntry(h4_LLP,"Signal, m_{#tilde{#mu}}= 500 GeV, m_{#tilde{#chi}^{0}}= 350 GeV","L");
  leg->Draw();

  leg = new TLegend(0.23,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  
  leg = new TLegend(0.23,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderC);
  leg->Draw();

leg = new TLegend(0.23,0.69,0.35,0.74);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(legsize);
  leg->SetHeader(HeaderCbis);
  leg->Draw();
PlotCMSv3(pad5,Year,DATA);




// pad7->cd();
// // // intégrale à droite
// TH1F* Inte  = new TH1F("ratio","",nbin,xmin,xmax);
// Inte->Sumw2();
// if (DATA)
//   {
// for (int i = 0; i< nbin-1; i++)
//   {
//     double sumPredi = 0;
//     double sumData = 0;
//     double Interatio = 0;
//     for (int j = i+1 ; j < nbin ; j++)
//       {
//           sumPredi += hsolve->GetBinContent(j);
//           sumData += htotData->GetBinContent(j);  
//       } 
//       if (sumPredi == 0){Inte->AddBinContent(i,0);} 
//       else {Inte->AddBinContent(i,sumData/sumPredi);}
 
// }
//   }
// if (!DATA)
//   {
//     for (int i = 0; i< nbin-1; i++)
//       {
//         double sumPredi = 0;
//         double sumData = 0;
//         double Interatio = 0;
//         for (int j = i+1 ; j < nbin ; j++)
//           {
//               sumPredi += hsolve->GetBinContent(j);
//               sumData += htotMC->GetBinContent(j);  
//           } 
//           if (sumPredi == 0){Inte->AddBinContent(i,0);} 
//           else {Inte->AddBinContent(i,sumData/sumPredi);}
    
//       }
//   }
//   Inte->Draw("PE1"); 
//   //  hsolve->SetFillStyle(3004);
// //  Inte->SetFillColor(kBlack);
//  Inte->SetLineColor(kBlack);
//  Inte->SetLineStyle(1);
//  Inte->SetLineWidth(1);
//  Inte->SetTickLength(0.03, "YZ");
//  Inte->SetTickLength(0.03,"X");
//  Inte->SetLabelOffset(0.015,"X");
//  Inte->SetLabelOffset(0.007,"Y");
//  Inte->SetLabelSize(0.045, "XYZ");
//  Inte->SetLabelFont(42, "XYZ"); 
//  Inte->SetTitleSize(0.045, "XYZ"); 
//  Inte->SetTitleFont(42, "XYZ");
//  Inte->SetTitleOffset(1.2,"X"); 
//  Inte->SetTitleOffset(1.3,"Y");
//  Inte->GetXaxis()->SetTitle(xtitle);
//  Inte->GetXaxis()->SetTitleColor(1);
//  Inte->GetYaxis()->SetTitle("Data/Prediction Integral ratio");
//  Inte->GetYaxis()->SetTitleColor(1);
//  Inte->SetNdivisions(509,"XYZ");
//  Inte->SetMinimum(0); 
//  Inte->SetMaximum(5); 
//  Inte->SetMarkerStyle(20);
//  Inte->SetMarkerSize(1);
//  Inte->GetXaxis()->SetRangeUser(0,20);
 /// !! --------------------------------------------------------------------------//
 /// !! --------------------------------------------------------------------------//

  c1->cd();
 pad3->cd();

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
 hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
 hmax = hmaxBD; 

 
 if(MC)
  {
 f1_DY->cd();
 
 g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleB);
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
  g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleB);
  h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
  g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleB);
  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
  g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleB);
  h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
  g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleB);
  h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
  g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleB);
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleB);
  h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleB);
  h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
  g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleB);
  h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleB);
  h_ST->Add(g2_ST, h_ST, 1, 1);

      //   f3_ST->cd();
      // g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleB);
      // h_ST->Add(g3_ST, h_ST, 1, 1);

      //   f4_ST->cd();
      // g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleB);
      // h_ST->Add(g4_ST, h_ST, 1, 1);

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleB);
 
  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, rwTT*1,0);

f2_TT->cd();
      g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleB);
      h_TT->Add(g2_TT, h_TT, 1,1);

//!! --------------------------
h_DY->Scale(1*scaleMC);
h_VV->Scale(1*scaleMC);
h_TTV->Scale(1*scaleMC);
h_ST->Scale(1*scaleMC);
h_TT->Scale(1*scaleMC);
// Daniel 
 h_VV->Add(h_DY, h_VV, 1, 1);
 h_ST->Add(h_ST, h_VV, 1, 1);
 h_TTV->Add(h_TTV, h_ST, 1, 1);
 htotMC->Add(h_TTV, htotMC, 1, 0);
 htotMC->Add(htotMC,h_TT, 1, 1);
// -----


    if (MC && !DATA)
      {
        htotMC->Draw("PE1");
        htotMC->SetMarkerStyle(20);
        htotMC->SetMarkerSize(1);
        htotMC->SetMarkerColor(kBlack);
        htotMC->SetLineColor(kBlack);
        htotMC->SetLineWidth(1);
      }
    if (MC && DATA)
      {
 htotMC->Draw("HE"); 
  htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(ColorRed, 1);
 htotMC->SetLineColor(ColorRed);
htotMC->SetLineColorAlpha(ColorRed, 1);
}
 htotMC->SetLineStyle(1);
 htotMC->SetLineWidth(1);
 htotMC->SetTickLength(0.03, "YZ");
 htotMC->SetTickLength(0.03,"X");
 htotMC->SetLabelOffset(0.015,"X");
 htotMC->SetLabelOffset(0.007,"Y");
 htotMC->SetLabelSize(0.045, "XYZ");
 htotMC->SetLabelFont(42, "XYZ"); 
 htotMC->SetTitleSize(0.055, "XYZ"); 
 htotMC->SetTitleFont(42, "XYZ");
 htotMC->SetTitleOffset(1.2,"X"); 
 htotMC->SetTitleOffset(1.3,"Y");
 htotMC->GetXaxis()->SetTitle(xtitle);
 htotMC->GetXaxis()->SetTitleColor(1);
 htotMC->GetYaxis()->SetTitle(ytitle);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
//  if (Method == 2 || Method == 3 || Method == 4 || Method == 5)
//    {
//       htotMC->GetXaxis()->SetRangeUser(0,20);
//    }

//  htotMC->SetMinimum(1); 
//  htotMC->SetMaximum(htotMC->GetMaximum()*2); 
 if (logy)
  {
    htotMC->SetMinimum(1); 
    htotMC->SetMaximum(htotMC->GetMaximum()*100); 
  }
else 
  {
    htotMC->SetMinimum(0); 
    htotMC->SetMaximum(htotMC->GetMaximum()*2); 
  }
//  htotMC->SetMarkerStyle(20);
//  htotMC->SetMarkerSize(1);

SumMCEvent += htotMC->Integral(0,nbin+1);
}
 


if (DATA)
  {
    f1_Data_emu->cd();
    g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleB);
    htotData->Add(g1_Data_emu, htotData, 1, 0);
htotData->Draw("PE1same");
htotData->SetMarkerStyle(20);
htotData->SetMarkerSize(1);
htotData->SetMarkerColor(kBlack);
htotData->SetLineColor(kBlack);
htotData->SetLineWidth(1);
SumDataEvent += htotData->Integral(0,nbin+1);
if (DATA && !MC)
      {
        htotData->SetTickLength(0.03, "YZ");
        htotData->SetTickLength(0.03,"X");
        htotData->SetLabelOffset(0.015,"X");
        htotData->SetLabelOffset(0.007,"Y");
        htotData->SetLabelSize(0.045, "XYZ");
        htotData->SetLabelFont(42, "XYZ"); 
        htotData->SetTitleSize(0.055, "XYZ"); 
        htotData->SetTitleFont(42, "XYZ");
        htotData->SetTitleOffset(1.2,"X"); 
        htotData->SetTitleOffset(1.3,"Y");
            htotData->GetXaxis()->SetTitle(xtitle);
        htotData->GetXaxis()->SetTitleColor(1);
        htotData->GetYaxis()->SetTitle(ytitle);
        htotData->GetYaxis()->SetTitleColor(1);
        // htotData->GetXaxis()->SetRangeUser(0,20);
                if (logy)
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*100); 
          }
        else 
          {
            htotData->SetMinimum(0); 
            htotData->SetMaximum(htotData->GetMaximum()*2); 
          }
      }
  }



if(Signal)
  {
 f1_LLP->cd();
  
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleB);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

 f2_LLP->cd();
 
 g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleB);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

 f3_LLP->cd();
 
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleB);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(kRed-1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);

 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(kRed-2);
 h2_LLP->SetLineStyle(2);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(kRed-3);
 h3_LLP->SetLineStyle(3);
 h3_LLP->SetLineWidth(2);
}


//  hsolve->Divide(hsolve, htotData, 1., 1.);
  if (DATA )
  {
 hsolve->Divide(hsolve, htotData, 1., 1.);
}
 else if (MC && !DATA)
  {
    hsolve->Divide(hsolve, htotMC, 1., 1.);
  }

  // leg = new TLegend(0.17,0.94,0.50,0.98);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.06);
  // leg->SetHeader(HeaderCMS);
  // leg->Draw();

  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
if (DATA)
    {
  leg->AddEntry(htotData, " e#mu data","PE1");
}
  if (MC && !DATA)
    {
  leg->AddEntry(htotMC, " Total MC","PE1");
  }
  if (MC && DATA)
    {
       leg->AddEntry(htotMC, " Total MC","F");
  }
if (Signal)
  {
    leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
}


  leg->Draw();

  leg = new TLegend(0.23,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  
  leg = new TLegend(0.23,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderB);
  leg->Draw();

    leg = new TLegend(0.23,0.69,0.35,0.74);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(legsize);
  leg->SetHeader(HeaderBbis);
  leg->Draw();

PlotCMSv3(pad3,Year,DATA);
if (DATA && MC)
  {
  pad10->cd();
  
  // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
  hDataMC->Divide(htotData,htotMC,1,1);
  hDataMC->Draw("PE1"); 
  hDataMC->SetFillStyle(1001);
  hDataMC->SetFillColorAlpha(kBlack, 1);
  hDataMC->SetLineColor(kBlack);
  hDataMC->SetLineStyle(1);
  hDataMC->SetLineWidth(1);
  hDataMC->SetTickLength(0.03, "YZ");
  hDataMC->SetTickLength(0.03,"X");
  hDataMC->SetLabelOffset(0.015,"X");
  hDataMC->SetLabelOffset(0.007,"Y");
  hDataMC->SetLabelSize(0.09, "XYZ");
  hDataMC->SetLabelFont(42, "XYZ"); 
  hDataMC->SetTitleSize(0.12, "XYZ"); 
  hDataMC->SetTitleFont(42, "XYZ");
  hDataMC->SetTitleOffset(0.8,"X"); 
  hDataMC->SetTitleOffset(0.5,"Y");
  hDataMC->GetXaxis()->SetTitle(xtitle);
  hDataMC->GetXaxis()->SetTitleColor(1);
  hDataMC->GetYaxis()->SetTitle(ytitle);
  hDataMC->GetYaxis()->SetTitleColor(1);
  hDataMC->SetNdivisions(509,"XYZ");
  hDataMC->SetMinimum(0); 
  hDataMC->SetMaximum(2);
// hDataMC->GetXaxis()->SetRangeUser(0,20);
}

// *****************************************************************************

 pad4->cd();
 hmax = hmaxBD; 

 htotData = new TH1F("htotData","",nbin,xmin,xmax);
 htotMC = new TH1F("htotMC","",nbin,xmin,xmax);
  hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);

 if (MC)
  {
 f1_DY->cd();
   g1_DY = (TH1F*)gROOT->FindObject(MCFILE[0]+htitleD);
  h_DY = new TH1F("h_DY","",nbin,xmin,xmax);
 h_DY->Add(g1_DY, h_DY, 1,0);

 f2_DY->cd();
 g2_DY = (TH1F*)gROOT->FindObject(MCFILE[1]+htitleD);
  h_DY->Add(g2_DY, h_DY, 1, 1);

 f1_VV->cd();
 g1_VV = (TH1F*)gROOT->FindObject(MCFILE[9]+htitleD);
  h_VV = new TH1F("h_VV","",nbin,xmin,xmax);
 h_VV->Add(g1_VV, h_VV, 1,0);

 f2_VV->cd();
 g2_VV = (TH1F*)gROOT->FindObject(MCFILE[10]+htitleD);
  h_VV->Add(g2_VV, h_VV, 1, 1);

 f3_VV->cd();
 g3_VV = (TH1F*)gROOT->FindObject(MCFILE[11]+htitleD);
  h_VV->Add(g3_VV, h_VV, 1, 1);

 f1_TTV->cd();
 g1_TTV = (TH1F*)gROOT->FindObject(MCFILE[6]+htitleD);
  h_TTV = new TH1F("h_TTV","",nbin,xmin,xmax);
 h_TTV->Add(g1_TTV, h_TTV, 1,0);

 f2_TTV->cd();
 g2_TTV = (TH1F*)gROOT->FindObject(MCFILE[7]+htitleD);
  h_TTV->Add(g2_TTV, h_TTV, 1, 1);

 f3_TTV->cd();
 g3_TTV = (TH1F*)gROOT->FindObject(MCFILE[8]+htitleD);
  h_TTV->Add(g3_TTV, h_TTV, 1, 1);

 h_ST = new TH1F("h_ST","",nbin,xmin,xmax);
 f1_ST->cd();
 g1_ST = (TH1F*)gROOT->FindObject(MCFILE[4]+htitleD);
  h_ST->Add(g1_ST, h_ST, 1,0);

 f2_ST->cd();
 g2_ST = (TH1F*)gROOT->FindObject(MCFILE[5]+htitleD);
  h_ST->Add(g2_ST, h_ST, 1, 1);

    //   f3_ST->cd();
    // g3_ST = (TH1F*)gROOT->FindObject(MCFILE[12]+htitleD);
    // h_ST->Add(g3_ST, h_ST, 1, 1);

    //   f4_ST->cd();
    // g4_ST = (TH1F*)gROOT->FindObject(MCFILE[13]+htitleD);
    // h_ST->Add(g4_ST, h_ST, 1, 1);

 f1_TT->cd();
g1_TT = (TH1F*)gROOT->FindObject(MCFILE[2]+htitleD);

  h_TT = new TH1F("h_TT","",nbin,xmin,xmax);
 h_TT->Add(g1_TT, h_TT, rwTT*1,0);

f2_TT->cd();
    g2_TT = (TH1F*)gROOT->FindObject(MCFILE[3]+htitleD);
    h_TT->Add(g2_TT, h_TT, 1,1);

//!! --------------------------
h_DY->Scale(1*scaleMC);
h_VV->Scale(1*scaleMC);
h_TTV->Scale(1*scaleMC);
h_ST->Scale(1*scaleMC);
h_TT->Scale(1*scaleMC);
// Daniel 
 h_VV->Add(h_DY, h_VV, 1, 1);
 h_ST->Add(h_ST, h_VV, 1, 1);
 h_TTV->Add(h_TTV, h_ST, 1, 1);
 htotMC->Add(h_TTV, htotMC, 1, 0);
 htotMC->Add(htotMC,h_TT, 1, 1);
// -----

if (MC && !DATA)
      {
        htotMC->Draw("PE1");
        htotMC->SetMarkerStyle(20);
        htotMC->SetMarkerSize(1);
        htotMC->SetMarkerColor(kBlack);
        htotMC->SetLineColor(kBlack);
        htotMC->SetLineWidth(1);
      }
    if (MC && DATA)
      {
 htotMC->Draw("HE"); 
  htotMC->SetFillStyle(1001);
 htotMC->SetFillColorAlpha(ColorRed, 1);
 htotMC->SetLineColor(ColorRed);
htotMC->SetLineColorAlpha(ColorRed, 1);
}
 htotMC->SetLineStyle(1);
 htotMC->SetLineWidth(1);
 htotMC->SetTickLength(0.03, "YZ");
 htotMC->SetTickLength(0.03,"X");
 htotMC->SetLabelOffset(0.015,"X");
 htotMC->SetLabelOffset(0.007,"Y");
 htotMC->SetLabelSize(0.045, "XYZ");
 htotMC->SetLabelFont(42, "XYZ"); 
 htotMC->SetTitleSize(0.055, "XYZ"); 
 htotMC->SetTitleFont(42, "XYZ");
 htotMC->SetTitleOffset(1.2,"X"); 
 htotMC->SetTitleOffset(1.3,"Y");
 htotMC->GetXaxis()->SetTitle(xtitle);
 htotMC->GetXaxis()->SetTitleColor(1);
 htotMC->GetYaxis()->SetTitle(ytitle);
 htotMC->GetYaxis()->SetTitleColor(1);
 htotMC->SetNdivisions(509,"XYZ");
//  if (Method == 2 || Method == 3 || Method == 4 || Method == 5)
//    {
//       htotMC->GetXaxis()->SetRangeUser(0,20);
//    }

 
 if (logy)
  {
    htotMC->SetMinimum(1); 
    htotMC->SetMaximum(htotMC->GetMaximum()*100); 
  }
else 
  {
    htotMC->SetMinimum(0); 
    htotMC->SetMaximum(htotMC->GetMaximum()*2); 
  }
SumMCEvent += htotMC->Integral(0,nbin+1);
}

if (DATA)
  {
    f1_Data_emu->cd();
    g1_Data_emu = (TH1F*)gROOT->FindObject(DATAFILE[0]+htitleD);
    htotData->Add(g1_Data_emu, htotData, 1, 0);

htotData->Draw("PE1same");
htotData->SetMarkerStyle(20);
htotData->SetMarkerSize(1);
htotData->SetMarkerColor(kBlack);
htotData->SetLineColor(kBlack);
htotData->SetLineWidth(1);
if (DATA && !MC)
      {
        SumDataEvent += htotData->Integral(0,nbin+1);
        htotData->SetTickLength(0.03, "YZ");
        htotData->SetTickLength(0.03,"X");
        htotData->SetLabelOffset(0.015,"X");
        htotData->SetLabelOffset(0.007,"Y");
        htotData->SetLabelSize(0.045, "XYZ");
        htotData->SetLabelFont(42, "XYZ"); 
        htotData->SetTitleSize(0.055, "XYZ"); 
        htotData->SetTitleFont(42, "XYZ");
        htotData->SetTitleOffset(1.2,"X"); 
        htotData->SetTitleOffset(1.3,"Y");
            htotData->GetXaxis()->SetTitle(xtitle);
        htotData->GetXaxis()->SetTitleColor(1);
        htotData->GetYaxis()->SetTitle(ytitle);
        htotData->GetYaxis()->SetTitleColor(1);
// htotData->GetXaxis()->SetRangeUser(0,20);
                if (logy)
          {
            htotData->SetMinimum(1); 
            htotData->SetMaximum(htotData->GetMaximum()*100); 
          }
        else 
          {
            htotData->SetMinimum(0); 
            htotData->SetMaximum(htotData->GetMaximum()*2); 
          }
      }
  }

if (Signal)
  {

 f1_LLP->cd();
 g1_LLP = (TH1F*)gROOT->FindObject(LLPFILE[0]+htitleD);
 g1_LLP->Sumw2();
 h1_LLP = new TH1F("h1_LLP","",nbin,xmin,xmax);
 h1_LLP->Add(g1_LLP, h1_LLP, 1,0);

 f2_LLP->cd();
  g2_LLP = (TH1F*)gROOT->FindObject(LLPFILE[1]+htitleD);
 g2_LLP->Sumw2();
 h2_LLP = new TH1F("h2_LLP","",nbin,xmin,xmax);
 h2_LLP->Add(g2_LLP, h2_LLP, 1,0);

 f3_LLP->cd();
 
 g3_LLP = (TH1F*)gROOT->FindObject(LLPFILE[2]+htitleD);
 g3_LLP->Sumw2();
 h3_LLP = new TH1F("h3_LLP","",nbin,xmin,xmax);
 h3_LLP->Add(g3_LLP, h3_LLP, 1,0);

 h1_LLP->Draw("HEsame"); 
 h1_LLP->SetLineColor(kRed-1);
 h1_LLP->SetLineStyle(1);
 h1_LLP->SetLineWidth(2);

 h2_LLP->Draw("HEsame"); 
 h2_LLP->SetLineColor(kRed-2);
 h2_LLP->SetLineStyle(2);
 h2_LLP->SetLineWidth(2);

 h3_LLP->Draw("HEsame");
 h3_LLP->SetLineColor(kRed-3);
 h3_LLP->SetLineStyle(3);
 h3_LLP->SetLineWidth(2);

}


//  hsolve->Multiply(hsolve, htotData, 1., 1.);
   if (DATA )
  {
 hsolve->Multiply(hsolve, htotData, 1., 1.);
}
 else if (MC && !DATA)
  {
    hsolve->Multiply(hsolve, htotMC, 1., 1.);
  }

  // leg = new TLegend(0.17,0.94,0.50,0.98);
  // leg->SetBorderSize(0);
  // leg->SetFillColor(kWhite);
  // leg->SetTextFont(42);
  // leg->SetTextSize(0.06);
  // leg->SetHeader(HeaderCMS);
  // leg->Draw();


  leg = new TLegend(0.69,LEGY1,0.89,LEGY2);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.045);
  leg->SetMargin(0.2);
if (DATA)
    {
  leg->AddEntry(htotData, " e#mu data","PE1");
}
  if (MC && !DATA)
    {
  leg->AddEntry(htotMC, " Total MC","PE1");
  }
  if (MC && DATA)
    {
       leg->AddEntry(htotMC, " Total MC","F");
  }
  if (Signal)
    {
      leg->AddEntry(h1_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 250 (200) GeV","L");
  leg->AddEntry(h2_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 300 (180) GeV","L");
  leg->AddEntry(h3_LLP," m_{#tilde{#mu} (#tilde{#chi}^{0})} = 400 (300) GeV","L");
}
  leg->Draw();

  leg = new TLegend(0.23,0.80,0.35,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderNVtx);
  leg->Draw();
  
  leg = new TLegend(0.23,0.75,0.35,0.79);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetHeader(HeaderD);
  leg->Draw();

leg = new TLegend(0.23,0.69,0.35,0.74);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(legsize);
  leg->SetHeader(HeaderDbis);
  leg->Draw();
PlotCMSv3(pad4,Year,DATA);


if (DATA && MC)
  {
  pad11->cd();
  
  // TH1F* hDataMC = new TH1F("hDataMC","",nbin,xmin,xmax);
  hDataMC->Divide(htotData,htotMC,1,1);
  hDataMC->Draw("PE1"); 
  hDataMC->SetFillStyle(1001);
  hDataMC->SetFillColorAlpha(kBlack, 1);
  hDataMC->SetLineColor(kBlack);
  hDataMC->SetLineStyle(1);
  hDataMC->SetLineWidth(1);
  hDataMC->SetTickLength(0.03, "YZ");
  hDataMC->SetTickLength(0.03,"X");
  hDataMC->SetLabelOffset(0.015,"X");
  hDataMC->SetLabelOffset(0.007,"Y");
  hDataMC->SetLabelSize(0.09, "XYZ");
  hDataMC->SetLabelFont(42, "XYZ"); 
  hDataMC->SetTitleSize(0.12, "XYZ"); 
  hDataMC->SetTitleFont(42, "XYZ");
  hDataMC->SetTitleOffset(0.8,"X"); 
  hDataMC->SetTitleOffset(0.5,"Y");
  hDataMC->GetXaxis()->SetTitle(xtitle);
  hDataMC->GetXaxis()->SetTitleColor(1);
  hDataMC->GetYaxis()->SetTitle(ytitle);
  hDataMC->GetYaxis()->SetTitleColor(1);
  hDataMC->SetNdivisions(509,"XYZ");
  hDataMC->SetMinimum(0); 
  hDataMC->SetMaximum(2);
// hDataMC->GetXaxis()->SetRangeUser(0,20);
}

// *****************************************************************************

//  pad2->cd();
//  hsolve->SetFillColor(kGray+1);
// const  TH1F*  hsolvetemp(hsolve);
// TGraphErrors TGsolve = new TGraphErrors(hsolvetemp);
//    TGsolve->SetFillColor(6);
//    TGsolve->SetFillStyle(3005);
//    TGsolve->Draw("a4same");
  
//  hsolve->SetLineColor(kBlack);
//  hsolve->SetFillColor(kBlack);
//  hsolve->SetFillStyle(3004);
// //  hsolve->SetLineStyle(1);
// //  hsolve->SetLineWidth(2);
//   hsolve->Draw("E4same");

  std::cout<<"EMU data 2 vertices : "<<SumDataEvent<<std::endl;
  std::cout<<"EMU MC 2 vertices : "<<SumMCEvent<<std::endl;
// *****************************************************************************

  c1->Update();

  //-------------------------------------------------//
//-------------------------------------------------//
//-------------------------------------------------//
//-------------------------------------------------//


  c2->cd();
  pad5->cd();

 if ( DATA && !MC ) {
   hsolve->SetLineColor(kGray+1);//kAzure+6
   hsolve->SetFillColor(kGray+1);
 }
 else {
   hsolve->SetLineColor(kGray+1);
   hsolve->SetFillColor(kGray+1);
 }
 hsolve->SetFillStyle(3354);
 hsolve->SetLineStyle(1);
 hsolve->SetLineWidth(1);
 hsolve->Draw("sameE2"); 

// hsolve->Draw("E4same");
  TString NAMELE =  Name+"_"+Plots;
  if (DATA && !MC) NAMELE += "_Data";
  if (MC && !DATA) NAMELE += "_MC";
  if (MC && DATA) NAMELE += "_DataMC";
  NAMELE += "_"+Year;
  hsolve->SaveAs(NAMELE+".root");
//-------------------------------------------------//
//-------------------------------------------------//
//-------------------------------------------------//
//-------------------------------------------------//

  pad6->cd();
//pulls between the predictions and data using the following formula:
//pull = (data - prediction) / sqrt(data+sigma^2_prediction)

// Poissionan errors are added due to low stats

 
TH1F* ratio  = new TH1F("ratio","",nbin,xmin,xmax);
ratio->Sumw2();
if (DATA )
  {
    for (int i = 0; i< nbin; i++)
      {
          double data = htotDatav2->GetBinContent(i);
          double pred = hsolve->GetBinContent(i);
          double sigma = hsolve->GetBinError(i);
          // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
          if (data == 0 || pred == 0) {ratio->AddBinContent(i,0);}    //when data == 0 and there is no prediction => you need this
          else
            {
              double pull = (data - pred) / sqrt(data + sigma*sigma);
              // double pull = data/pred;
              ratio->AddBinContent(i,pull);
            }
      }
  }
if ( !DATA)
  {
    for (int i = 0; i< nbin; i++)
      {
        double data = htotMCv2->GetBinContent(i);
        double pred = hsolve->GetBinContent(i);
        double sigma = hsolve->GetBinError(i);
        // std::cout<<"data = "<<data<<" pred = "<<pred<<" sigma = "<<sigma<<std::endl;
        if (data < 1 || pred < 1) {ratio->AddBinContent(i,0);}    //when data == 0 and there is no prediction => you need this
        else
          {
            double pull = (data - pred) / sqrt(data + sigma*sigma);
            // double pull = (data-pred)/pred;
            ratio->AddBinContent(i,pull);
          }
      }
 }
 ratio->Draw("PE1"); // PE11 ou E4
 ratio->SetFillColor(kRed);
 ratio->SetFillStyle(3001);
 ratio->SetLineColor(kBlack);
 ratio->SetLineStyle(1);
 ratio->SetLineWidth(1);
 ratio->SetTickLength(0.03, "YZ");
 ratio->SetTickLength(0.03,"X");
 ratio->SetLabelOffset(0.015,"X");
 ratio->SetLabelOffset(0.007,"Y");
 ratio->SetLabelSize(0.07, "XYZ");
 ratio->SetLabelFont(42, "XYZ"); 
 ratio->SetTitleSize(0.07, "XYZ"); 
 ratio->SetTitleFont(42, "XYZ");
 ratio->SetTitleOffset(1.2,"X"); 
 ratio->SetTitleOffset(0.8,"Y");
 ratio->GetXaxis()->SetTitle(xtitle);
 ratio->GetXaxis()->SetTitleColor(1);
 ratio->GetYaxis()->SetTitle("Pulls");
 ratio->GetYaxis()->SetTitleColor(1);
 ratio->SetNdivisions(509,"XYZ");
 ratio->SetMinimum(-3); 
 ratio->SetMaximum(3); 
 ratio->SetMarkerStyle(20);
 ratio->SetMarkerSize(1);
// ratio->GetXaxis()->SetRangeUser(0,20);

//-------------------------------------------------//
//-------------------------------------------------//
//-------------------------------------------------//
//-------------------------------------------------//


  c1->cd();
  pad2->cd();


  hsolve->SetFillStyle(3354);
 hsolve->SetLineStyle(1);
 hsolve->SetLineWidth(1);
 hsolve->Draw("sameE2");


  TString namele =  Name+"_"+Plots+"_v2";
  if (DATA && !MC) namele += "_Data";
  if (MC && !DATA) namele += "_MC";
  if (MC && DATA) namele += "_DataMC";
  namele += "_"+Year;
  c2->SaveAs(namele+".pdf");

  return c1;
}