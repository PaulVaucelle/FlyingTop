void test_color_png() {
  TH1F* h = new TH1F("h",";x;N",50,-3,3);
  h->FillRandom("gaus",3000);

  Float_t r=0.246, g=0.563, b=0.852;
  Int_t myBlue = TColor::GetColor(r,g,b);

  TCanvas c("c","c",800,600);
  h->SetLineColor(myBlue);
  h->SetLineColorAlpha(myBlue,1);
  h->SetLineWidth(2);
  h->Draw("HIST");

  c.SaveAs("histo.pdf");  // doit sortir bleu
}