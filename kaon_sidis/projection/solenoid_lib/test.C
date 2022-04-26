void test(){
  gSystem->Load("libSole.so");
  
  TH1F *h1;
  TH1F *h2 = new TH1F("h2","h2",360,0.,360.);
  TH1F *h3 = new TH1F("h3","h3",360,0.,360.);
  
  ofstream outfile("1.log");
  Double_t x[360];
  TCanvas *c1 = new TCanvas("c1","c1",800,400);
  c1->SetFillColor(10);
  c1->Divide(2,1);
  TLatex *t = new TLatex();
  Double_t abc[4];
  for (Int_t a = 1;a!=3;a++){
    for (Int_t b = 1; b!=7;b++){
      for (Int_t c=1;c!=9;c++){
	TFile file(Form("./rootfiles/%d_%d_%d.root",a,b,c));
	for (Int_t d = 0; d!= file.GetNkeys();d++){
	  h1 = (TH1F*)file.Get(Form("h%d",d));
	  c1->cd(1);
	  h1->Draw();
	  for (Int_t i=0;i!=360;i++){
	    x[i] = h1->GetBinContent(i+1);
	  }
	  sole_inter sole;
	  sole.init(360,x);
	  sole.connect();
	  sole.get_new(x);
	  sole.get_limit(abc);
	  TLine l1(abc[0],0,abc[0],100000000);
	  TLine l2(abc[1],0,abc[1],100000000);
	  TLine l3(abc[2],0,abc[2],100000000);
	  TLine l4(abc[3],0,abc[3],100000000);
	  l1.Draw("same");
	  l2.Draw("same");
	  l3.Draw("same");
	  l4.Draw("same");
	  for (Int_t i=0;i!=360;i++){
	    h2->SetBinContent(i+1,x[i]);
	  }
	  h2->SetLineColor(2);
	  h2->Draw("same");
	  c1->cd(2);
	  Double_t sum=0;
	  Int_t count = 0;
	  for (Int_t i=0;i!=360;i++){
	    x[i] = sole.get_accep(i+0.5);
	    h3->SetBinContent(i+1,x[i]);
	    if (x[i]!=0){
	      sum+= x[i];
	      count++;
	    }
	  }
	  
	  h3->Draw();
	  t->DrawTextNDC(0.6,0.6,Form("%4.3f",sum/count));
	  c1->Update();
	  c1->Print(Form("./figures/%d_%d_%d_%d.gif",a,b,c,d));
	  h2->Reset();
	  h3->Reset();
	  if (fabs(sum/count-1)>0.02){
	    outfile << a << "\t" << b << "\t" << c << "\t" << d << "\t" << sum/count << endl;
	  }
	}
	file.Close();
      }
    }
  }
  outfile.close();
}
