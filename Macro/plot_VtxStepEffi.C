{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Sat Jun 22 18:16:16 2024) by ROOT version5.34/36
   gROOT->LoadMacro("./tdrstyle.C");

   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",1067,244,807,762);
   Canvas_1->Range(0.5,-0.1159888,5.5,1.043899);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
  //  //DYM50
  //  TH1F *VtxStepEffi = new TH1F("VtxStepEffi","VTx step Efficiency",4,1,5);
  //  VtxStepEffi->SetBinContent(0,152.3681);
  //  VtxStepEffi->SetBinContent(1,0.004616726);
  //  VtxStepEffi->SetBinContent(2,0.0002001058);
  //  VtxStepEffi->SetBinContent(3,0.8837243);
  //  VtxStepEffi->SetBinContent(4,0.1114589);
  //  // VtxStepEffi->SetEntries(1.073009e+07);
   
  //  //TTTo2L2Nu
  //     TH1F *VtxStepEffi2 = new TH1F("VtxStepEffi2","VTx step Efficiency",4,1,5);
  //  VtxStepEffi2->SetBinContent(0,13.80965);
  //  VtxStepEffi2->SetBinContent(1,0.007654184);
  //  VtxStepEffi2->SetBinContent(2,0.0002913842);
  //  VtxStepEffi2->SetBinContent(3,0.860804);
  //  VtxStepEffi2->SetBinContent(4,0.1312505);

   //AllSignal
      TH1F *VtxStepEffi = new TH1F("VtxStepEffi3","VTx step Efficiency",4,1,5);
   VtxStepEffi->SetBinContent(0,0.2883459);
   VtxStepEffi->SetBinContent(1,0.8133029);
   VtxStepEffi->SetBinContent(2,0.01044735);
   VtxStepEffi->SetBinContent(3,0.1601977);
   VtxStepEffi->SetBinContent(4,0.01605212);   

  //001
      TH1F *VtxStepEffi2 = new TH1F("VtxStepEffi2","VTx step Efficiency",4,1,5);
   VtxStepEffi2->SetBinContent(0,0.2526842);
   VtxStepEffi2->SetBinContent(1,0.7532045);
   VtxStepEffi2->SetBinContent(2,0.01236923);
   VtxStepEffi2->SetBinContent(3,0.2089717);
   VtxStepEffi2->SetBinContent(4,0.02545462);
   
   //003
   TH1F *VtxStepEffi3 = new TH1F("VtxStepEffi3","VTx step Efficiency",4,1,5);
   VtxStepEffi3->SetBinContent(0,0.1002829);
   VtxStepEffi3->SetBinContent(1,0.8736691);
   VtxStepEffi3->SetBinContent(2,0.009561902);
   VtxStepEffi3->SetBinContent(3,0.1050304);
   VtxStepEffi3->SetBinContent(4,0.01173854);
   
   //010
   TH1F *VtxStepEffi4 = new TH1F("VtxStepEffi4","VTx step Efficiency",4,1,5);
   VtxStepEffi4->SetBinContent(0,0.06434297);
   VtxStepEffi4->SetBinContent(1,0.9042034);
   VtxStepEffi4->SetBinContent(2,0.009283277);
   VtxStepEffi4->SetBinContent(3,0.07847837);
   VtxStepEffi4->SetBinContent(4,0.008034962);

  //030
   TH1F *VtxStepEffi5 = new TH1F("VtxStepEffi5","VTx step Efficiency",4,1,5);
   VtxStepEffi5->SetBinContent(0,0.06884967);
   VtxStepEffi5->SetBinContent(1,0.8875582);
   VtxStepEffi5->SetBinContent(2,0.009604014);
   VtxStepEffi5->SetBinContent(3,0.09354954);
   VtxStepEffi5->SetBinContent(4,0.00928821);
  //100
   TH1F *VtxStepEffi6 = new TH1F("VtxStepEffi6","VTx step Efficiency",4,1,5);
   VtxStepEffi6->SetBinContent(0,0.1585304);
   VtxStepEffi6->SetBinContent(1,0.8008596);
   VtxStepEffi6->SetBinContent(2,0.0115398);
   VtxStepEffi6->SetBinContent(3,0.1715609);
   VtxStepEffi6->SetBinContent(4,0.01603969);
  //300
   TH1F *VtxStepEffi7 = new TH1F("VtxStepEffi7","VTx step Efficiency",4,1,5);
   VtxStepEffi7->SetBinContent(0,0.544585);
   VtxStepEffi7->SetBinContent(1,0.691083);
   VtxStepEffi7->SetBinContent(2,0.01148293);
   VtxStepEffi7->SetBinContent(3,0.2733677);
   VtxStepEffi7->SetBinContent(4,0.02406638);
  //1000
     TH1F *VtxStepEffi8 = new TH1F("VtxStepEffi","VTx step Efficiency",4,1,5);
   VtxStepEffi8->SetBinContent(0,1.893163);
   VtxStepEffi8->SetBinContent(1,0.6004118);
   VtxStepEffi8->SetBinContent(2,0.009195838);
   VtxStepEffi8->SetBinContent(3,0.3594809);
   VtxStepEffi8->SetBinContent(4,0.03091146);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   VtxStepEffi->SetLineColor(kBlack);
   VtxStepEffi->GetXaxis()->SetTitle("Vtx step");
   VtxStepEffi->GetYaxis()->SetTitle("Efficiency");
   VtxStepEffi->GetXaxis()->SetLabelFont(42);
   VtxStepEffi->GetXaxis()->SetLabelSize(0.035);
   VtxStepEffi->GetXaxis()->SetTitleSize(0.035);
   VtxStepEffi->GetXaxis()->SetTitleFont(42);
   VtxStepEffi->GetYaxis()->SetLabelFont(42);
   VtxStepEffi->GetYaxis()->SetLabelSize(0.035);
   VtxStepEffi->GetYaxis()->SetTitleSize(0.035);
   VtxStepEffi->GetYaxis()->SetTitleOffset(1.5);
   VtxStepEffi->GetYaxis()->SetTitleFont(42);
   VtxStepEffi->GetZaxis()->SetLabelFont(42);
   VtxStepEffi->GetZaxis()->SetLabelSize(0.035);
   VtxStepEffi->GetZaxis()->SetTitleSize(0.035);
   VtxStepEffi->GetZaxis()->SetTitleFont(42);
   VtxStepEffi->SetLineWidth(2.0);
   VtxStepEffi->GetYaxis()->SetRangeUser(0,1.0);
   VtxStepEffi->SetStats(0);
   VtxStepEffi->Draw("");

   VtxStepEffi2->SetLineColor(kRed);
   VtxStepEffi2->SetLineWidth(2.0);
   VtxStepEffi2->Draw("same");

   VtxStepEffi3->SetLineColor(kGreen+1);
   VtxStepEffi3->SetLineWidth(2.0);
   VtxStepEffi3->Draw("same");


   VtxStepEffi4->SetLineColor(kBlue+2);
   VtxStepEffi4->SetLineWidth(2.0);
   VtxStepEffi4->Draw("same");


   VtxStepEffi5->SetLineColor(kYellow+1);
   VtxStepEffi5->SetLineWidth(2.0);
   VtxStepEffi5->Draw("same");


   VtxStepEffi6->SetLineColor(40);
   VtxStepEffi6->SetLineWidth(2.0);
   VtxStepEffi6->Draw("same");


   VtxStepEffi7->SetLineColor(41);
   VtxStepEffi7->SetLineWidth(2.0);
   VtxStepEffi7->Draw("same");


   VtxStepEffi8->SetLineColor(45);
   VtxStepEffi8->SetLineWidth(2.0);
   VtxStepEffi8->Draw("same");
   
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

   TLegend *legend = new TLegend(0.7,0.7,0.89,0.9);
   legend->SetHeader("MC Samples");
   legend->AddEntry(VtxStepEffi,"Signal Average","l");
   legend->AddEntry(VtxStepEffi2,"c#tau = 0.1 cm","l");
   legend->AddEntry(VtxStepEffi3,"c#tau = 0.3 cm","l");
   legend->AddEntry(VtxStepEffi4,"c#tau = 1.0 cm","l");
   legend->AddEntry(VtxStepEffi5,"c#tau = 3.0 cm","l");
   legend->AddEntry(VtxStepEffi6,"c#tau = 10.0 cm","l");
   legend->AddEntry(VtxStepEffi7,"c#tau = 30.0 cm","l");
   legend->AddEntry(VtxStepEffi8,"c#tau = 100.0 cm","l");
   legend->Draw();

   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
   Canvas_1->SaveAs("VtxStepEffi_Svsctau.pdf");

}
