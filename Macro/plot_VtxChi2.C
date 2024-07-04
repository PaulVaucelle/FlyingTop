{
//=========Macro generated from canvas: Canvas_1/Canvas_1
//=========  (Wed Jun 26 11:59:49 2024) by ROOT version5.34/36
   TCanvas *Canvas_1 = new TCanvas("Canvas_1", "Canvas_1",2795,152,792,660);
   Canvas_1->Range(-3,-540.75,17,4866.75);
   Canvas_1->SetFillColor(0);
   Canvas_1->SetBorderMode(0);
   Canvas_1->SetBorderSize(2);
   Canvas_1->SetFrameBorderMode(0);
   Canvas_1->SetFrameBorderMode(0);
   
   //300-180-100
   TH1F *hData_Hemi_Vtx_NChi2 = new TH1F("hData_Hemi_Vtx_NChi2","",64,-1,15);
   hData_Hemi_Vtx_NChi2->SetBinContent(5,4120);
   hData_Hemi_Vtx_NChi2->SetBinContent(6,1874);
   hData_Hemi_Vtx_NChi2->SetBinContent(7,1843);
   hData_Hemi_Vtx_NChi2->SetBinContent(8,1487);
   hData_Hemi_Vtx_NChi2->SetBinContent(9,1110);
   hData_Hemi_Vtx_NChi2->SetBinContent(10,735);
   hData_Hemi_Vtx_NChi2->SetBinContent(11,492);
   hData_Hemi_Vtx_NChi2->SetBinContent(12,332);
   hData_Hemi_Vtx_NChi2->SetBinContent(13,284);
   hData_Hemi_Vtx_NChi2->SetBinContent(14,225);
   hData_Hemi_Vtx_NChi2->SetBinContent(15,200);
   hData_Hemi_Vtx_NChi2->SetBinContent(16,126);
   hData_Hemi_Vtx_NChi2->SetBinContent(17,118);
   hData_Hemi_Vtx_NChi2->SetBinContent(18,79);
   hData_Hemi_Vtx_NChi2->SetBinContent(19,66);
   hData_Hemi_Vtx_NChi2->SetBinContent(20,46);
   hData_Hemi_Vtx_NChi2->SetBinContent(21,58);
   hData_Hemi_Vtx_NChi2->SetBinContent(22,31);
   hData_Hemi_Vtx_NChi2->SetBinContent(23,31);
   hData_Hemi_Vtx_NChi2->SetBinContent(24,33);
   hData_Hemi_Vtx_NChi2->SetBinContent(25,41);
   hData_Hemi_Vtx_NChi2->SetBinContent(26,29);
   hData_Hemi_Vtx_NChi2->SetBinContent(27,25);
   hData_Hemi_Vtx_NChi2->SetBinContent(28,22);
   hData_Hemi_Vtx_NChi2->SetBinContent(29,21);
   hData_Hemi_Vtx_NChi2->SetBinContent(30,17);
   hData_Hemi_Vtx_NChi2->SetBinContent(31,19);
   hData_Hemi_Vtx_NChi2->SetBinContent(32,14);
   hData_Hemi_Vtx_NChi2->SetBinContent(33,18);
   hData_Hemi_Vtx_NChi2->SetBinContent(34,13);
   hData_Hemi_Vtx_NChi2->SetBinContent(35,10);
   hData_Hemi_Vtx_NChi2->SetBinContent(36,10);
   hData_Hemi_Vtx_NChi2->SetBinContent(37,6);
   hData_Hemi_Vtx_NChi2->SetBinContent(38,3);
   hData_Hemi_Vtx_NChi2->SetBinContent(39,7);
   hData_Hemi_Vtx_NChi2->SetBinContent(40,11);
   hData_Hemi_Vtx_NChi2->SetBinContent(41,7);
   hData_Hemi_Vtx_NChi2->SetBinContent(42,6);
   hData_Hemi_Vtx_NChi2->SetBinContent(43,4);
   hData_Hemi_Vtx_NChi2->SetBinContent(44,9);
   hData_Hemi_Vtx_NChi2->SetEntries(13582);
   hData_Hemi_Vtx_NChi2->Sumw2();
   hData_Hemi_Vtx_NChi2->Scale(1/hData_Hemi_Vtx_NChi2->Integral(0,-1));
   hData_Hemi_Vtx_NChi2->GetYaxis()->SetRangeUser(0,1);
   hData_Hemi_Vtx_NChi2->SetStats(0);


   //TTTo2L2Nu
   TH1F *hData_Hemi_Vtx_NChi2_TT = new TH1F("hData_Hemi_Vtx_NChi2_TT","",64,-1,15);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(5,643926);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(6,76240);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(7,54009);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(8,42304);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(9,34780);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(10,29822);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(11,25319);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(12,22165);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(13,20899);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(14,21369);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(15,18737);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(16,17528);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(17,16433);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(18,11066);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(19,9307);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(20,7653);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(21,6781);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(22,6039);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(23,5435);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(24,5029);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(25,4569);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(26,4168);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(27,3913);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(28,3559);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(29,3420);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(30,3169);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(31,2874);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(32,2693);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(33,2517);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(34,2438);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(35,2272);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(36,2177);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(37,2134);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(38,1899);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(39,1837);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(40,1736);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(41,1689);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(42,1548);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(43,1463);
   hData_Hemi_Vtx_NChi2_TT->SetBinContent(44,1402);
   hData_Hemi_Vtx_NChi2_TT->SetEntries(1126318);
   
hData_Hemi_Vtx_NChi2_TT->Sumw2();
   hData_Hemi_Vtx_NChi2_TT->Scale(1/hData_Hemi_Vtx_NChi2_TT->Integral(0,-1));
   hData_Hemi_Vtx_NChi2_TT->SetStats(0);
   
   
   // TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   // ptstats->SetName("stats");
   // ptstats->SetBorderSize(1);
   // ptstats->SetFillColor(0);
   // ptstats->SetTextAlign(12);
   // ptstats->SetTextFont(42);
   // TText *text = ptstats->AddText("hData_Hemi_Vtx_NChi2");
   // text->SetTextSize(0.0368);
   // text = ptstats->AddText("Entries = 13582  ");
   // text = ptstats->AddText("Mean  =  0.942");
   // text = ptstats->AddText("RMS   =  1.231");
   // ptstats->SetOptStat(1111);
   // ptstats->SetOptFit(0);
   // ptstats->Draw();
   // hData_Hemi_Vtx_NChi2->GetListOfFunctions()->Add(ptstats);
   // ptstats->SetParent(hData_Hemi_Vtx_NChi2);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hData_Hemi_Vtx_NChi2->SetLineColor(ci);
   hData_Hemi_Vtx_NChi2->GetXaxis()->SetLabelFont(42);
   hData_Hemi_Vtx_NChi2->GetXaxis()->SetLabelSize(0.035);
   hData_Hemi_Vtx_NChi2->GetXaxis()->SetTitleSize(0.035);
   hData_Hemi_Vtx_NChi2->GetXaxis()->SetTitleFont(42);
   hData_Hemi_Vtx_NChi2->GetXaxis()->SetTitle("Vertex #chi^{2}");
   hData_Hemi_Vtx_NChi2->GetYaxis()->SetLabelFont(42);
   hData_Hemi_Vtx_NChi2->GetYaxis()->SetLabelSize(0.035);
   hData_Hemi_Vtx_NChi2->GetYaxis()->SetTitleSize(0.035);
   hData_Hemi_Vtx_NChi2->GetYaxis()->SetTitleFont(42);
   hData_Hemi_Vtx_NChi2->GetYaxis()->SetTitle("Entries");
   hData_Hemi_Vtx_NChi2->GetYaxis()->SetTitleOffset(1.5);
   hData_Hemi_Vtx_NChi2->GetZaxis()->SetLabelFont(42);
   hData_Hemi_Vtx_NChi2->GetZaxis()->SetLabelSize(0.035);
   hData_Hemi_Vtx_NChi2->GetZaxis()->SetTitleSize(0.035);
   hData_Hemi_Vtx_NChi2->GetZaxis()->SetTitleFont(42);
   hData_Hemi_Vtx_NChi2->Draw("L");


   hData_Hemi_Vtx_NChi2_TT->SetLineColor(kRed);
   hData_Hemi_Vtx_NChi2_TT->SetLineWidth(2);
   hData_Hemi_Vtx_NChi2_TT->GetYaxis()->SetRangeUser(0,1);
   hData_Hemi_Vtx_NChi2_TT->Draw("LSAME");

   TLegend *legend = new TLegend(0.75,0.75,0.9,0.9);
   legend->SetHeader("MC Samples");
   legend->AddEntry(hData_Hemi_Vtx_NChi2_TT,"TTTo2L2Nu","l");
   legend->AddEntry(hData_Hemi_Vtx_NChi2,"Signal","l");
   legend->Draw();

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
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t-0.01,lumiText);

      latex.SetTextFont(cmsTextFont);
      latex.SetTextAlign(11); 
      latex.SetTextSize(cmsTextSize*t);    
      latex.DrawLatex(l,1-t+lumiTextOffset*t-0.01,cmsText);

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
        latex.DrawLatex(posX_+0.11, posY_-0.01, extraText);
	    }


   Canvas_1->Modified();
   Canvas_1->cd();
   Canvas_1->SetSelected(Canvas_1);
   Canvas_1->SaveAs("plot_VtxChi.pdf");
}
