{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Sat Jun 22 17:58:39 2024) by ROOT version5.34/36
      gROOT->LoadMacro("./tdrstyle.C");
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",1028,139,807,762);
   Canvas_1->Range(0.5,-0.1068211,5.5,0.9613902);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   // All signal then 001 003 010 ...
   TH1F *VtxStepEffi = new TH1F("VtxStepEffi","VTx step Efficiency",4,1,5);
   VtxStepEffi->SetBinContent(0,0.2883459);
   VtxStepEffi->SetBinContent(1,0.8133029);
   VtxStepEffi->SetBinContent(2,0.01044735);
   VtxStepEffi->SetBinContent(3,0.1601977);
   VtxStepEffi->SetBinContent(4,0.01605212);
   // VtxStepEffi->SetBinError(0,0.0003408201);
   // VtxStepEffi->SetBinError(1,0.0005723932);
   // VtxStepEffi->SetBinError(2,6.487411e-05);
   // VtxStepEffi->SetBinError(3,0.0002540367);
   // VtxStepEffi->SetBinError(4,8.041452e-05);
   // VtxStepEffi->SetEntries(3198128);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   VtxStepEffi->SetLineColor(ci);
   VtxStepEffi->GetXaxis()->SetLabelFont(42);
   VtxStepEffi->GetXaxis()->SetTitle("Vtx step");
   VtxStepEffi->GetYaxis()->SetTitle("Efficiency");
   VtxStepEffi->GetXaxis()->SetLabelSize(0.035);
   VtxStepEffi->GetXaxis()->SetTitleSize(0.035);
   VtxStepEffi->GetXaxis()->SetTitleFont(42);
   VtxStepEffi->GetYaxis()->SetLabelFont(42);
   VtxStepEffi->GetYaxis()->SetLabelSize(0.035);
   VtxStepEffi->GetYaxis()->SetTitleSize(0.035);
   VtxStepEffi->GetYaxis()->SetTitleFont(42);
   VtxStepEffi->GetZaxis()->SetLabelFont(42);
   VtxStepEffi->GetZaxis()->SetLabelSize(0.035);
   VtxStepEffi->GetZaxis()->SetTitleSize(0.035);
   VtxStepEffi->GetZaxis()->SetTitleFont(42);
   VtxStepEffi->SetLineWidth(2.0);
   VtxStepEffi->SetStats(0);
   VtxStepEffi->Draw("");
   

      VtxStepEffi2->SetLineColor(kRed);
   VtxStepEffi2->SetLineWidth(2.0);
   VtxStepEffi2->Draw("same");

   VtxStepEffi3->SetLineColor(kGreen+1);
   VtxStepEffi3->SetLineWidth(2.0);
   VtxStepEffi3->Draw("same");

   
   TPaveText *pt = new TPaveText(0.2676398,0.935,0.7323602,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("VTx step Efficiency");
   pt->Draw();


         setTDRStyle();
   //------Start of Copy Paste
TString cmsText     = "CMS";
float cmsTextFont   = 61;  // default is helvetic-bold

bool writeExtraText = true;
TString extraText   = "Private Work";
float extraTextFont = 52;  // default is helvetica-italics

// text sizes and text offsets with respect to the top frame
// in unit of the top margin size
float lumiTextSize     = 0.6;
float lumiTextOffset   = 0.2;
float cmsTextSize      = 0.75;
float cmsTextOffset    = 0.1;  // only used in outOfFrame version

float relPosX    = 0.045;
float relPosY    = 0.035;
float relExtraDY = 1.2;

// ratio of "CMS" and extra text size
float extraOverCmsTextSize  = 0.76;

TString lumi_13TeV = "";//137 fb^{-1}
TString lumi_sqrtS = " 13 TeV";
TString lumiText = lumi_13TeV+lumi_sqrtS;
  float H = Canvas_1->GetWh();
  float W = Canvas_1->GetWw();
  float l = Canvas_1->GetLeftMargin();
  float t = Canvas_1->GetTopMargin();
  float r = Canvas_1->GetRightMargin();
  float b = Canvas_1->GetBottomMargin();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    

  float extraTextSize = extraOverCmsTextSize*cmsTextSize;

  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);

      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);

float posX_=0;
  float posY_=0;
  int iPosX = 3;
  if( iPosX%10<=1 )
    {
      posX_ =   l + relPosX*(1-l-r);
    }
  else if( iPosX%10==2 )
    {
      posX_ =  l + 0.5*(1-l-r);
    }
  else if( iPosX%10==3 )
    {
      posX_ =  1-r - relPosX*(1-l-r);
    }
  float posY_ = 1-t - relPosY*(1-t-b);
  	  if( writeExtraText ) 
	    {
         posX_ =   l +  relPosX*(1-l-r);
         posY_ =   1-t+lumiTextOffset*t;
        int alignY_=3;
         int alignX_=2;
         if( iPosX/10==0 ) alignX_=1;
         if( iPosX==0    ) alignX_=1;
         if( iPosX==0    ) alignY_=1;
         if( iPosX/10==1 ) alignX_=1;
         if( iPosX/10==2 ) alignX_=2;
         if( iPosX/10==3 ) alignX_=3;
         //if( iPosX == 0  ) relPosX = 0.12;
         int align_ = 10*alignX_ + alignY_;
      latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(11);
      latex.DrawLatex(posX_+0.125, posY_, extraText);
	    }

   TLegend *legend = new TLegend(0.75,0.75,0.89,0.9);
   legend->SetHeader("MC Samples");
   legend->AddEntry(VtxStepEffi,"DYM50","l");
   legend->AddEntry(VtxStepEffi2,"TTTo2L2Nu","l");
   legend->AddEntry(VtxStepEffi3,"Signal Average","l");
   legend->Draw();

   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
}
