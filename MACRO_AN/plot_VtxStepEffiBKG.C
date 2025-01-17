{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Sat Jun 22 17:58:39 2024) by ROOT version5.34/36
      // gROOT->LoadMacro("./tdrstyle.C");
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",1028,139,807,762);
   Canvas_1->Range(0.5,-0.1068211,5.5,0.9613902);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   

Float_t r1 = 0.246;
Float_t g1 = 0.563;
Float_t b1 = 0.852;
TColor color1 = TColor(301,r1, g1, b1);
// color1.SetRGB(r1, g1, b1);
Int_t ColorBlue = color1.GetNumber();

Float_t r2 = 1.000;
Float_t g2 = 0.661;
Float_t b2 = 0.055;
TColor color2 = TColor(302,r2, g2, b2);
// color2.SetRGB(r2, g2, b2);
Int_t ColorOrange = color2.GetNumber();

Float_t r3 = 0.739;
Float_t g3 = 0.122;
Float_t b3 = 0.004;
TColor color3 = TColor(303,r3, g3, b3);
Int_t ColorRed = color3.GetNumber();

Float_t r4 = 0.578;
Float_t g4 = 0.641;
Float_t b4 = 0.635;
TColor color4 = TColor(304,r4, g4, b4);
Int_t ColorGrey = color4.GetNumber();

Float_t r5 = 0.513;
Float_t g5 = 0.176;
Float_t b5 = 0.713;
TColor color5 = TColor(305,r5, g5, b5);
Int_t ColorDarkPurple = color5.GetNumber();

Float_t r6 = 0.661;
Float_t g6 = 0.418;
Float_t b6 = 0.348;
TColor color6 = TColor(306,r6, g6, b6);
Int_t ColorBrown = color6.GetNumber();

Float_t r7 = 0.905;
Float_t g7 = 0.387;
Float_t b7 = 0.000;
TColor color7 = TColor(307,r7, g7, b7);
Int_t ColorDarkOrange = color7.GetNumber();

Float_t r8 = 0.723;
Float_t g8 = 0.672;
Float_t b8 = 0.438;
TColor color8 = TColor(308,r8, g8, b8);
Int_t ColorNeutral = color8.GetNumber();

Float_t r9 = 0.441;
Float_t g9 = 0.457;
Float_t b9 = 0.504;
TColor color9 = TColor(309,r9, g9, b9);
Int_t ColorDarkGrey = color9.GetNumber();

Float_t r10 = 0.571;
Float_t g10 = 0.852;
Float_t b10 = 0.867;
TColor color10 = TColor(310,r10, g10, b10);
// color10.SetRGB(r10, g10, b10);
Int_t ColorLightBlue = color10.GetNumber();

   // DYM50
   TH1F *VtxStepEffi = new TH1F("VtxStepEffi","",4,1,5);
   VtxStepEffi->SetBinContent(0,152.3681);
   VtxStepEffi->SetBinContent(1,0.004616726);
   VtxStepEffi->SetBinContent(2,0.0002001058);
   VtxStepEffi->SetBinContent(3,0.8837243);
   VtxStepEffi->SetBinContent(4,0.1114589);
   // VtxStepEffi->SetBinError(0,0.0003408201);
   // VtxStepEffi->SetBinError(1,0.0005723932);
   // VtxStepEffi->SetBinError(2,6.487411e-05);
   // VtxStepEffi->SetBinError(3,0.0002540367);
   // VtxStepEffi->SetBinError(4,8.041452e-05);
   // VtxStepEffi->SetEntries(3198128);

  //TTTo2L2Nu
  TH1F *VtxStepEffi2 = new TH1F("VtxStepEffi2","",4,1,5);
     VtxStepEffi2->SetBinContent(0,13.80965);
   VtxStepEffi2->SetBinContent(1,0.007654184);
   VtxStepEffi2->SetBinContent(2,0.0002913842);
   VtxStepEffi2->SetBinContent(3,0.860804);
   VtxStepEffi2->SetBinContent(4,0.1312505);
  //  VtxStepEffi2->SetEntries(2.06712e+07);

      TH1F *VtxStepEffi3 = new TH1F("VtxStepEffi3","",4,1,5);
   VtxStepEffi3->SetBinContent(0,0.2883459);
   VtxStepEffi3->SetBinContent(1,0.8133029);
   VtxStepEffi3->SetBinContent(2,0.01044735);
   VtxStepEffi3->SetBinContent(3,0.1601977);
   VtxStepEffi3->SetBinContent(4,0.01605212);   

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   VtxStepEffi->SetLineColor(ColorOrange);
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
   VtxStepEffi->SetMaximum(1.0);

      VtxStepEffi2->SetLineColor(ColorRed);
   VtxStepEffi2->SetLineWidth(2.0);
   VtxStepEffi2->Draw("same");

   VtxStepEffi3->SetLineColor(ColorLightBlue);
   VtxStepEffi3->SetLineWidth(2.0);
   VtxStepEffi3->Draw("same");

   
   TPaveText *pt = new TPaveText(0.2676398,0.935,0.7323602,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("");
   pt->Draw();


        //  setTDRStyle();
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
   posY_ = 1-t - relPosY*(1-t-b);
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

   TLegend *legend = new TLegend(0.7,0.75,0.89,0.87);
      legend ->SetBorderSize(0);
      legend->SetFillStyle(0);
      legend->SetFillColor(kWhite);
      legend->SetTextFont(42);
      legend->SetTextSize(0.025);
   legend->SetHeader("MC Samples");
   legend->AddEntry(VtxStepEffi,"DYM50","l");
   legend->AddEntry(VtxStepEffi2,"TTTo2L2Nu","l");
   legend->AddEntry(VtxStepEffi3,"Signal Average","l");
   legend->Draw();

   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
   Canvas_1->SaveAs("plot_VtxStepEffiBKG.pdf");
}
