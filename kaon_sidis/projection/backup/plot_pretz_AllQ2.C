double max1 = 8.1;
double min1 = -0.1;
double qmax1 = 0.0145;
double qmin1 = -0.0145;
double xmin1 = -0.01;
double xmax1 = 0.7;

/*void plot_neutron(target_flag,particle_flag){{{*/
void plot_neutron(Int_t target_flag = 3,Int_t i=1){
  //flag_t: 0->pretzelocity, 1->collins, 2->sivers
  //flag: 1->ploting Asym axis

  Float_t factor,shift;
  factor=(max1-min1)/(qmax1-qmin1);
  shift=max1/factor-qmax1;
  factor = 1; shift = 0.0;

  target_flag = 3;
  int particle_flag =i;

  /*Input{{{*/
  TString target = "X";
  if(target_flag==1)
    target ="p";
  else if(target_flag==2)
    target ="d2";
  else if(target_flag==3)
    target ="3he";
  else{
    cerr<<"I don't know this particle flag!"<<endl;
  }
  TString particle = "X";
  if(particle_flag==1){
    particle ="pip";
    qmin1 = -0.0025;
    qmax1 = 0.0065;
  }
  else if(particle_flag==2){
    particle ="pim";
    qmin1 = -0.0065;
    qmax1 = 0.0025;
  }
  else{
    cerr<<"I don't know this particle flag!"<<endl;
  }
  TString filename;
  filename.Form("./pretz_both_AllQ2/%s_%s.dat",target.Data(),particle.Data());
  ifstream infile(filename);

  /*}}}*/

  gStyle->SetOptStat(0);
  Int_t count1,count2;
  Double_t Q2[5000],x[5000],z[5000],pt[5000],y[5000],A[5000],Astat[5000],N[5000],coverage[5000],coef[3][5000];
  Double_t Q2[5000],x[5000];
  Double_t temp;
  Int_t pt_flag = -1;
  Int_t ncount=0;

  Int_t Q2_temp,x_temp;
  while(infile >> Q2_temp >> count2){
	  for (Int_t j1=0;j1<count2;j1++){
		  infile >> Q2_temp >> x_temp >> z[ncount] >> Q2[ncount] >> pt[ncount] >> x[ncount] >> y[ncount] >> Astat[ncount] >> N[ncount]>>
			  coverage[ncount] >> coef[0][ncount] >> coef[1][ncount] >> coef[2][ncount];

		  Astat[ncount] *= coef[2][ncount]; //coef[0]->Sivers,coef[1]->Collins, coef[2]->Pretzelosity

		  if (Astat[ncount]>0.&&Astat[ncount]<0.05){
			  cerr<<Form("---Q2=%f, x = %f, Asys= %f",Q2[ncount], x[ncount], Astat[ncount])<<endl;
			  A[ncount]=0.0;
			  Astat[ncount] *= factor; //coef[0]->Sivers,coef[1]->Collins, coef[2]->Pretzelosity
			  ncount ++;
		  }
	  }
  }

  infile.close();
  TString titlename;
  TLatex *t1 = new TLatex();
  t1->SetNDC();

  Float_t a[2]={0,1};
  Float_t b[2]={0,1};
  TGraph *g2 = new TGraph(2,a,b);
  g2->GetYaxis()->SetRangeUser(qmin1,qmax1);
  g2->GetXaxis()->SetRangeUser(xmin1,xmax1);
  g2->GetXaxis()->SetTitle("x");
  //g2->GetYaxis()->SetTitle("A_{UT}^{sin(3#phi-#phi_{S})}");
  g2->GetYaxis()->SetTitle("Pretzelocity Asymmetry");
  g2->GetXaxis()->SetNdivisions(506);
  g2->GetYaxis()->SetNdivisions(506);
  g2->GetXaxis()->CenterTitle(0);
  g2->GetXaxis()->SetTitleFont(22);
  g2->GetYaxis()->CenterTitle(0);
  g2->GetYaxis()->SetTitleFont(22);
  g2->GetXaxis()->SetTitleSize(0.06);
  g2->GetYaxis()->SetTitleSize(0.05);  
  g2->GetYaxis()->SetTitleOffset(1.4);
  g2->GetXaxis()->SetTitleOffset(0.8);

  g2->SetTitle(0);
  g2->Draw("AP");

  TGraphErrors *g1 = new TGraphErrors(ncount,x,A,0,Astat);
  g1->Draw("*same");
  g1->SetMarkerColor(2);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(2.0);
  g1->SetTitle();
  if(particle_flag==1){
    titlename.Form("#scale[1.2]{(#pi^{+}, neutron)}");
    t1->DrawLatex(0.20,0.85,titlename);
    titlename.Form("#scale[0.8]{<Q^{2}>=%2.1f GeV^{2}}", 2.5);
    t1->DrawLatex(0.65,0.20,titlename);
  }
  if(particle_flag==2){
    titlename.Form("#scale[1.2]{(#pi^{-}, neutron)}");
    t1->DrawLatex(0.20,0.85,titlename);
    titlename.Form("#scale[0.8]{<Q^{2}>=%2.1f GeV^{2}}", 2.5);
    t1->DrawLatex(0.65,0.20,titlename);
  }

  TLine *bb = new TLine(xmin1,0,xmax1,0);
  bb->Draw("same");
}
/*}}}*/

/*void plot_neutron_new(target_flag,particle_flag){{{*/
void plot_neutron_new(Int_t target_flag = 3,Int_t i=1){
  //flag_t: 0->pretzelocity, 1->collins, 2->sivers
  //flag: 1->ploting Asym axis

  Float_t factor,shift;
  factor=(max1-min1)/(qmax1-qmin1);
  shift=max1/factor-qmax1;
  factor = 1; shift = 0.0;

  target_flag = 3;
  int particle_flag =i;

  /*Input{{{*/
  TString target = "X";
  if(target_flag==1)
    target ="p";
  else if(target_flag==2)
    target ="d2";
  else if(target_flag==3)
    target ="3he";
  else{
    cerr<<"I don't know this particle flag!"<<endl;
  }
  TString particle = "X";
  if(particle_flag==1){
    particle ="pip";
    qmin1 = -0.0025;
    qmax1 = 0.0065;
  }
  else if(particle_flag==2){
    particle ="pim";
    qmin1 = -0.0065;
    qmax1 = 0.0025;
  }
  else{
    cerr<<"I don't know this particle flag!"<<endl;
  }
  TString filename;
  filename.Form("./pretz_both_AllQ2/%s_%s.dat",target.Data(),particle.Data());
  ifstream infile(filename);

  /*}}}*/

  gStyle->SetOptStat(0);
  Int_t count1,count2;
  Double_t Q2[5000],x[5000],z[5000],pt[5000],y[5000],A[5000],Astat[5000],N[5000],coverage[5000],coef[3][5000];
  Double_t Q2[5000],x[5000];
  Double_t temp;
  Int_t pt_flag = -1;
  Int_t ncount=0;

  Int_t Q2_temp,x_temp;
  while(infile >> Q2_temp >> count2){
	  for (Int_t j1=0;j1<count2;j1++){
		  infile >> Q2_temp >> x_temp >> z[ncount] >> Q2[ncount] >> pt[ncount] >> x[ncount] >> y[ncount] >> Astat[ncount] >> N[ncount]>>
			  coverage[ncount] >> coef[0][ncount] >> coef[1][ncount] >> coef[2][ncount];

		  Astat[ncount] *= coef[2][ncount]; //coef[0]->Sivers,coef[1]->Collins, coef[2]->Pretzelosity

		  if (Astat[ncount]>0.&&Astat[ncount]<0.05){
			  cerr<<Form("---Q2=%f, x = %f, Asys= %f",Q2[ncount], x[ncount], Astat[ncount])<<endl;
			  A[ncount]=0.0;
			  Astat[ncount] *= factor; //coef[0]->Sivers,coef[1]->Collins, coef[2]->Pretzelosity
			  ncount ++;
		  }
	  }
  }

  infile.close();
  TString titlename;
  TLatex *t1 = new TLatex();
  t1->SetNDC();

  Float_t a[2]={0,1};
  Float_t b[2]={0,1};
  TGraph *g2 = new TGraph(2,a,b);
  g2->GetYaxis()->SetRangeUser(qmin1,qmax1);
  g2->GetXaxis()->SetRangeUser(xmin1,xmax1);
  g2->GetXaxis()->SetTitle("x");
  //g2->GetYaxis()->SetTitle("A_{UT}^{sin(3#phi-#phi_{S})}");
  g2->GetYaxis()->SetTitle("Pretzelocity Asymmetry");
  g2->GetXaxis()->SetNdivisions(506);
  g2->GetYaxis()->SetNdivisions(506);
  g2->GetXaxis()->CenterTitle(0);
  g2->GetXaxis()->SetTitleFont(22);
  g2->GetYaxis()->CenterTitle(0);
  g2->GetYaxis()->SetTitleFont(22);
  g2->GetXaxis()->SetTitleSize(0.06);
  g2->GetYaxis()->SetTitleSize(0.05);  
  g2->GetYaxis()->SetTitleOffset(1.8);
  g2->GetXaxis()->SetTitleOffset(0.8);

  g2->SetTitle(0);
  g2->Draw("AP");

  TGraphErrors *g1 = new TGraphErrors(ncount,x,A,0,Astat);
  g1->Draw("*same");
  g1->SetMarkerColor(2);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(2.0);
  g1->SetTitle();
  if(particle_flag==1){
    titlename.Form("#scale[1.4]{(#pi^{+}, neutron)}");
    t1->DrawLatex(0.20,0.85,titlename);
 //   titlename.Form("#scale[1.0]{<Q^{2}>=%2.1f GeV^{2}}", 2.5);
 //   t1->DrawLatex(0.5,0.20,titlename);
  }
  if(particle_flag==2){
    titlename.Form("#scale[1.4]{(#pi^{-}, neutron)}");
    t1->DrawLatex(0.20,0.85,titlename);
 //   titlename.Form("#scale[1.0]{<Q^{2}>=%2.1f GeV^{2}}", 2.5);
 //   t1->DrawLatex(0.5,0.20,titlename);
  }

  TLine *bb = new TLine(xmin1,0,xmax1,0);
  bb->Draw("same");
}
/*}}}*/

/*void plot_neutron_sbs(target_flag,particle_flag){{{*/
void plot_neutron_sbs(Int_t target_flag = 3,Int_t i=1){
 
  target_flag = 3;
  int particle_flag =i;

  gStyle->SetOptStat(0);
  TString titlename;
  TLatex *t1 = new TLatex();
  t1->SetNDC();

  Float_t a[2]={0,1};
  Float_t b[2]={0,1};
  TGraph *g2 = new TGraph(2,a,b);
  g2->GetYaxis()->SetRangeUser(qmin1,qmax1);
  g2->GetXaxis()->SetRangeUser(xmin1,xmax1);
  g2->GetXaxis()->SetTitle("x");
  //g2->GetYaxis()->SetTitle("A_{UT}^{sin(3#phi-#phi_{S})}");
  g2->GetYaxis()->SetTitle("Pretzelocity Asymmetry");
  g2->GetXaxis()->SetNdivisions(506);
  g2->GetYaxis()->SetNdivisions(506);
  g2->GetXaxis()->CenterTitle(0);
  g2->GetXaxis()->SetTitleFont(22);
  g2->GetYaxis()->CenterTitle(0);
  g2->GetYaxis()->SetTitleFont(22);
  g2->GetXaxis()->SetTitleSize(0.06);
  g2->GetYaxis()->SetTitleSize(0.05);  
  g2->GetYaxis()->SetTitleOffset(1.8);
  g2->GetXaxis()->SetTitleOffset(0.8);

  g2->SetTitle(0);
  g2->Draw("AP");
  
  /*SBS Binning{{{*/
  double x_pip_sbs[5], Asys_pip_sbs[5];
  double x_pim_sbs[5], Asys_pim_sbs[5];
  double zero_sbs[5]={0.0,0.0,0.0,0.0,0.0}; 

  ifstream sbs_pip("SBS_SIDIS_x_pip.txt");  
  ifstream sbs_pim("SBS_SIDIS_x_pim.txt");  

  double xx[10], NN[10], dummy;
  TString Com;

  sbs_pip >> Com >> Com >> Com >> Com >> Com >> Com;
  for(int i=0;i<10;i++){
	  sbs_pip >> Com >> dummy >> xx[i] >> dummy >> NN[i] >> dummy;
	  cerr<<Form("pi+: #%d  x = %f, N = %f",i, xx[i],NN[i])<<endl;
  }
  sbs_pip.close();
  for(int i=0;i<5;i++){
	  if(fabs(xx[i]-xx[i+5])<0.001){
		  Asys_pip_sbs[i] = 1./sqrt(NN[i]+NN[i+5]);
		  x_pip_sbs[i] = xx[i];
		  cerr<<"--pi+: Found one, x = "<< x_pip_sbs[i]<<"  Asys = "<<Asys_pip_sbs[i]<<endl;
	  }
  }

  sbs_pim >> Com >> Com >> Com >> Com >> Com >> Com;
  for(int i=0;i<10;i++){
	  sbs_pim >> Com >> dummy >> xx[i] >> dummy >> NN[i] >> dummy;
	  cerr<<Form("pi-: #%d  x = %f, N = %f",i, xx[i],NN[i])<<endl;
  }
  sbs_pim.close();
  for(int i=0;i<5;i++){
	  if(fabs(xx[i]-xx[i+5])<0.001){
		  Asys_pim_sbs[i] = 1./sqrt(NN[i]+NN[i+5]);
		  x_pim_sbs[i] = xx[i];
		  cerr<<"--pi-: Found one, x = "<< x_pim_sbs[i]<<"  Asys = "<<Asys_pim_sbs[i]<<endl;
	  }
  }
  /*}}}*/
  
  /*SoLID Binning{{{*/
  ifstream solid_pip("./pretz_both_AllQ2/3he_pip_xfix.dat");  
  ifstream solid_pim("./pretz_both_AllQ2/3he_pim_xfix.dat");  

double zero_solid[6]={0.0,0.0,0.0,0.0,0.0,0.0};
double x_pip_solid[6]={0.1,0.2,0.3,0.4,0.5,0.6};
double x_pim_solid[6]={0.1,0.2,0.3,0.4,0.5,0.6};

  double Asys_pip_solid[6], Asys_pim_solid[6];
  int dum;

  cerr<<"--- SoLID Pip "<<endl;
  for(int i=0;i<6;i++){
	  solid_pip >> dum >> dum >> dummy >> dummy >> dummy >> dummy >> dummy >> Asys_pip_solid[i]
		  >> dummy >> dummy >> dummy >> dummy >> dummy;

	  x_pip_solid[i] -= 0.02;	
	  cerr<<Form("solid-pi+: #%d  Asys = %f",dum, Asys_pip_solid[i])<<endl;
  }
  solid_pip.close();

  cerr<<"--- SoLID Pim "<<endl;
  for(int i=0;i<6;i++){
	  solid_pim >> dum >> dum >> dummy >> dummy >> dummy >> dummy >> dummy >> Asys_pim_solid[i]
		  >> dummy >> dummy >> dummy >> dummy >> dummy;

	  x_pim_solid[i] -= 0.02;	
	  cerr<<Form("solid-pi-: #%d  Asys = %f",dum, Asys_pim_solid[i])<<endl;
  }
  solid_pim.close();
  /*}}}*/

  TGraphErrors *gr1 = new TGraphErrors(5, x_pip_sbs, zero_sbs, zero_sbs, Asys_pip_sbs);
  gr1->SetMarkerColor(2);
  gr1->SetMarkerStyle(21);
  gr1->SetLineColor(2);
  gr1->SetMarkerSize(1.5);
  gr1->SetTitle();

  TGraphErrors *gr11 = new TGraphErrors(6, x_pip_solid, zero_solid, zero_solid, Asys_pip_solid);
  gr11->SetMarkerColor(4);
  gr11->SetMarkerStyle(20);
  gr11->SetLineColor(4);
  gr11->SetMarkerSize(1.5);
  gr11->SetTitle();

  TGraphErrors *gr2 = new TGraphErrors(5, x_pim_sbs, zero_sbs, zero_sbs, Asys_pim_sbs);
  gr2->SetMarkerColor(2);
  gr2->SetMarkerStyle(21);
  gr2->SetLineColor(2);
  gr2->SetMarkerSize(1.5);
  gr2->SetTitle();

  TGraphErrors *gr22 = new TGraphErrors(6, x_pip_solid, zero_solid, zero_solid, Asys_pip_solid);
  gr22->SetMarkerColor(4);
  gr22->SetMarkerStyle(20);
  gr22->SetLineColor(4);
  gr22->SetMarkerSize(1.5);
  gr22->SetTitle();

  if(particle_flag==1){
     gr1->Draw("Psame");
     gr11->Draw("Psame");
  	  titlename.Form("#scale[1.4]{(#pi^{+}, neutron)}");
    t1->DrawLatex(0.20,0.85,titlename);
 //   titlename.Form("#scale[1.0]{<Q^{2}>=%2.1f GeV^{2}}", 2.5);
 //   t1->DrawLatex(0.5,0.20,titlename);
  }
  if(particle_flag==2){
     gr2->Draw("Psame");
     gr22->Draw("Psame");
     titlename.Form("#scale[1.4]{(#pi^{-}, neutron)}");
    t1->DrawLatex(0.20,0.85,titlename);
 //   titlename.Form("#scale[1.0]{<Q^{2}>=%2.1f GeV^{2}}", 2.5);
 //   t1->DrawLatex(0.5,0.20,titlename);
  }

  TLine *bb = new TLine(xmin1,0,xmax1,0);
  bb->Draw("same");
}
/*}}}*/

/*void plot_pretz_evolved(Int_t part=1,Int_t flag_t = 0){{{*/
void plot_pretz_evolved(Int_t targ, Int_t part){

  //targ: 1->proton, 2->neutron
  //part; 1->pip, 2->pim

  /*S_Wave{{{*/  
  ifstream infile1("./calculations/BARBARA/aut_sin_3phi_minus_phis_sdwaves_evolved.dat");
  Float_t x_s[100],ppip_s[100],ppim_s[100],ppi0_s[100];
  Float_t dpip_s[100],dpim_s[100],dpi0_s[100];
  Float_t npip_s[100],npim_s[100],npi0_s[100];

  Int_t i_full=0,i_p=0,i_s=0;
  while(!infile1.eof()){
	  infile1 >> x_s[i_s] 
		  >> ppip_s[i_s] >> ppim_s[i_s] >> ppi0_s[i_s] 
		  >> dpip_s[i_s] >> dpim_s[i_s] >> dpi0_s[i_s] 
		  >> npip_s[i_s] >> npim_s[i_s] >> npi0_s[i_s];

	  i_s++;
  }
  infile1.close();
  cerr<<Form("--- SD-Wave = %d Points", i_s-1)<<endl;

  TGraph *qqabc_ppip_s = new TGraph(i_s-1,x_s,ppip_s);
  TGraph *qqabc_ppim_s = new TGraph(i_s-1,x_s,ppim_s);
  TGraph *qqabc_dpip_s = new TGraph(i_s-1,x_s,dpip_s);
  TGraph *qqabc_dpim_s = new TGraph(i_s-1,x_s,dpim_s);
  TGraph *qqabc_npip_s = new TGraph(i_s-1,x_s,npip_s);
  TGraph *qqabc_npim_s = new TGraph(i_s-1,x_s,npim_s);

  qqabc_ppip_s->SetLineColor(2);
  qqabc_ppip_s->SetLineWidth(2.5);
  qqabc_ppim_s->SetLineColor(2);
  qqabc_ppim_s->SetLineWidth(2.5);
  qqabc_npip_s->SetLineColor(2);
  qqabc_npip_s->SetLineWidth(2.5);
  qqabc_npim_s->SetLineColor(2);
  qqabc_npim_s->SetLineWidth(2.5);
  qqabc_dpip_s->SetLineColor(2);
  qqabc_dpip_s->SetLineWidth(2.5);
  qqabc_dpim_s->SetLineColor(2);
  qqabc_dpim_s->SetLineWidth(2.5);
  /*}}}*/
  
  /*P_Wave{{{*/  
  ifstream infile2("./calculations/BARBARA/aut_sin_3phi_minus_phis_pwaves_evolved.dat");
  Float_t x_p[100],ppip_p[100],ppim_p[100],ppi0_p[100];
  Float_t dpip_p[100],dpim_p[100],dpi0_p[100];
  Float_t npip_p[100],npim_p[100],npi0_p[100];

  while(!infile2.eof()){
	  infile2 >> x_p[i_p] 
		  >> ppip_p[i_p] >> ppim_p[i_p] >> ppi0_p[i_p] 
		  >> dpip_p[i_p] >> dpim_p[i_p] >> dpi0_p[i_p] 
		  >> npip_p[i_p] >> npim_p[i_p] >> npi0_p[i_p];

	  i_p++;
  }
  infile2.close();
  cerr<<Form("--- P-Wave = %d Points", i_p-1)<<endl;

  TGraph *qqabc_ppip_p = new TGraph(i_p-1,x_p,ppip_p);
  TGraph *qqabc_ppim_p = new TGraph(i_p-1,x_p,ppim_p);
  TGraph *qqabc_dpip_p = new TGraph(i_p-1,x_p,dpip_p);
  TGraph *qqabc_dpim_p = new TGraph(i_p-1,x_p,dpim_p);
  TGraph *qqabc_npip_p = new TGraph(i_p-1,x_p,npip_p);
  TGraph *qqabc_npim_p = new TGraph(i_p-1,x_p,npim_p);

  qqabc_ppip_p->SetLineColor(4);
  qqabc_ppip_p->SetLineWidth(2.5);
  qqabc_ppim_p->SetLineColor(4);
  qqabc_ppim_p->SetLineWidth(2.5);
  qqabc_npip_p->SetLineColor(4);
  qqabc_npip_p->SetLineWidth(2.5);
  qqabc_npim_p->SetLineColor(4);
  qqabc_npim_p->SetLineWidth(2.5);
  qqabc_dpip_p->SetLineColor(4);
  qqabc_dpip_p->SetLineWidth(2.5);
  qqabc_dpim_p->SetLineColor(4);
  qqabc_dpim_p->SetLineWidth(2.5);
  /*}}}*/
  
  /*Total{{{*/  
  ifstream infile3("./calculations/BARBARA/aut_sin_3phi_minus_phis_evolved.dat");
  Float_t x_full[100],ppip_full[100],ppim_full[100],ppi0_full[100];
  Float_t dpip_full[100],dpim_full[100],dpi0_full[100];
  Float_t npip_full[100],npim_full[100],npi0_full[100];

  while(!infile3.eof()){
	  infile3 >> x_full[i_full] 
		  >> ppip_full[i_full] >> ppim_full[i_full] >> ppi0_full[i_full] 
		  >> dpip_full[i_full] >> dpim_full[i_full] >> dpi0_full[i_full] 
		  >> npip_full[i_full] >> npim_full[i_full] >> npi0_full[i_full];

	  i_full++;
  }
  infile3.close();

  TGraph *qqabc_ppip = new TGraph(i_full-1,x_full,ppip_full);
  TGraph *qqabc_ppim = new TGraph(i_full-1,x_full,ppim_full);
  TGraph *qqabc_dpip = new TGraph(i_full-1,x_full,dpip_full);
  TGraph *qqabc_dpim = new TGraph(i_full-1,x_full,dpim_full);
  TGraph *qqabc_npip = new TGraph(i_full-1,x_full,npip_full);
  TGraph *qqabc_npim = new TGraph(i_full-1,x_full,npim_full);

  qqabc_ppip->SetLineColor(1);
  qqabc_ppip->SetLineWidth(2.5);
  qqabc_ppim->SetLineColor(1);
  qqabc_ppim->SetLineWidth(2.5);
  qqabc_npip->SetLineColor(1);
  qqabc_npip->SetLineWidth(2.5);
  qqabc_npim->SetLineColor(1);
  qqabc_npim->SetLineWidth(2.5);
  qqabc_dpip->SetLineColor(1);
  qqabc_dpip->SetLineWidth(2.5);
  qqabc_dpim->SetLineColor(1);
  qqabc_dpim->SetLineWidth(2.5);
  /*}}}*/

  if (targ==1){
    if (part==1){
      qqabc_ppip->Draw("SL");
      qqabc_ppip_s->Draw("SL");
      qqabc_ppip_p->Draw("SL");
    }else{
      qqabc_ppim->Draw("SL");
      qqabc_ppim_s->Draw("SL");
      qqabc_ppim_p->Draw("SL");
    }
  }

  if (targ==2){
    if (part==1){
		qqabc_npip->Draw("SL");
		qqabc_npip_s->Draw("SL");
		qqabc_npip_p->Draw("SL");
	}else{
		qqabc_npim->Draw("SL");
		qqabc_npim_s->Draw("SL");
		qqabc_npim_p->Draw("SL");
	}
  }
/*
  TString titlename;
  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->SetTextColor(4);
  t1->SetTextFont(32);

  TLatex *t2 = new TLatex();
  t2->SetNDC();
  t2->SetTextColor(2);
  t2->SetTextFont(32);

  TLatex *t3 = new TLatex();
  t3->SetNDC();
  t3->SetTextColor(1);
  t3->SetTextFont(32);

  if(part==1){
	  titlename.Form("#scale[1.0]{P-P int.}");
	  t1->DrawLatex(0.725,0.38,titlename);
	  titlename.Form("#scale[1.0]{S-D int.}");
	  t2->DrawLatex(0.36,0.53,titlename);
	  titlename.Form("#scale[1.0]{TOT}");
	  t3->DrawLatex(0.33,0.64,titlename);
  }else{
	  titlename.Form("#scale[1.0]{P-P int.}");
	  t1->DrawLatex(0.665,0.67,titlename);
	  titlename.Form("#scale[1.0]{S-D int.}");
	  t2->DrawLatex(0.35,0.52,titlename);
	  titlename.Form("#scale[1.0]{TOT}");
	  t3->DrawLatex(0.35,0.38,titlename);
  }

*/

}
/*}}}*/

/*pretz(){{{*/
void pretz(){
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"Y");
  //gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"X");

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->SetFillColor(10);
  c1->cd();
  plot_neutron(2,1);//(target,particle,Q2)
  plot_pretz_evolved(2,1);//(target, particle)

  
  TLegend *l1 = new TLegend(0.5,0.80,0.89,0.94);
  Float_t x=0.,y=0.;
  
  TGraphErrors *g3 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g4 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g5 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g6 = new TGraphErrors(1,&x,&y,0,0);
  g3->SetLineColor(2);
  g3->SetLineWidth(2.5);

  g4->SetLineColor(4);
  g4->SetLineWidth(2.5);

  g5->SetLineColor(1);
  g5->SetLineWidth(2.5);

  g6->SetMarkerColor(2);
  g6->SetMarkerStyle(20);
  g6->SetMarkerSize(1.2);
  l1->AddEntry(g3," ","L");
  l1->AddEntry(g4,"Pasquini et. al.","L");
  l1->AddEntry(g5," ","L");
  l1->AddEntry(g6,"90 days SoLID","P");
  l1->SetTextSize(0.05);
  l1->SetTextFont(22);
  l1->Draw();
  
  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  c2->SetFillColor(10);
  //c2->Divide(1,2);
  c2->cd();
  plot_neutron(2,2);//(target,particle,Q2)
  plot_pretz_evolved(2,2);//(target, particle)
  l1->Draw();


  c1->Print("SoLID_neutron_pip_pretz_AllQ2.pdf");
  c1->Print("SoLID_neutron_pip_pretz_AllQ2.png");
  c2->Print("SoLID_neutron_pim_pretz_AllQ2.pdf");
  c2->Print("SoLID_neutron_pim_pretz_AllQ2.png");
}
/*}}}*/

/*pretz_new(){{{*/
void pretz_new(){
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"Y");
  //gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetPadBorderMode(0);

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->SetFillColor(10);
  c1->Divide(2,1);
  c1->cd(1);
//  gPad->SetLeftMargin(0.16);
//  gPad->SetRightMargin(0.002);
//  gPad->SetTopMargin(0.005);
//  gPad->SetBottomMargin(0.10);
  plot_neutron_new(2,1);//(target,particle,Q2)
  plot_pretz_evolved(2,1);//(target, particle)

  
  TLegend *l1 = new TLegend(0.6,0.80,0.98,0.94);
  Float_t x=0.,y=0.;
  
  TGraphErrors *g3 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g4 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g5 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g6 = new TGraphErrors(1,&x,&y,0,0);
  g3->SetLineColor(2);
  g3->SetLineWidth(2.5);

  g4->SetLineColor(4);
  g4->SetLineWidth(2.5);

  g5->SetLineColor(1);
  g5->SetLineWidth(2.5);

  g6->SetMarkerColor(2);
  g6->SetMarkerStyle(20);
  g6->SetMarkerSize(1.2);
  l1->AddEntry(g3,"S-D int., ","L");
  l1->AddEntry(g4,"P-P int.,Pasquini et. al.","L");
  l1->AddEntry(g5,"Total ","L");
  l1->AddEntry(g6,"90 days SoLID","P");
  l1->SetTextSize(0.04);
  l1->SetTextFont(22);
  l1->Draw();
  
 // TCanvas *c2 = new TCanvas("c2","c2",1000,800);
 // c2->SetFillColor(10);
  //c2->Divide(1,2);
  c1->cd(2);
  //gPad->SetLeftMargin(0.16);
  //gPad->SetRightMargin(0.002);
  //gPad->SetTopMargin(0.005);
  //gPad->SetBottomMargin(0.10);
   plot_neutron_new(2,2);//(target,particle,Q2)
  plot_pretz_evolved(2,2);//(target, particle)
  l1->Draw();


  c1->Print("SoLID_neutron_both_pretz_AllQ2.pdf");
  c1->Print("SoLID_neutron_both_pretz_AllQ2.png");
}
/*}}}*/

/*pretz_sbs(){{{*/
void pretz_sbs(){
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"Y");
  //gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"X");
  gStyle->SetPadBorderMode(0);

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  c1->SetFillColor(10);
  c1->Divide(2,1);
  c1->cd(1);
//  gPad->SetLeftMargin(0.16);
//  gPad->SetRightMargin(0.002);
//  gPad->SetTopMargin(0.005);
//  gPad->SetBottomMargin(0.10);
  plot_neutron_sbs(2,1);//(target,particle,Q2)
  plot_pretz_evolved(2,1);//(target, particle)

  
  TLegend *l1 = new TLegend(0.6,0.80,0.98,0.94);
  Float_t x=0.,y=0.;
  
  TGraphErrors *g3 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g4 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g5 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g6 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g7 = new TGraphErrors(1,&x,&y,0,0);
  g3->SetLineColor(2);
  g3->SetLineWidth(2.5);

  g4->SetLineColor(4);
  g4->SetLineWidth(2.5);

  g5->SetLineColor(1);
  g5->SetLineWidth(2.5);

  g6->SetMarkerColor(4);
  g6->SetMarkerStyle(20);
  g6->SetMarkerSize(1.2);

  g7->SetMarkerColor(2);
  g7->SetMarkerStyle(21);
  g7->SetMarkerSize(1.2);
  l1->AddEntry(g3,"S-D int. ","L");
  l1->AddEntry(g4,"P-P int. (Pasquini et. al.)","L");
  l1->AddEntry(g5,"Total","L");
l1->AddEntry(g6,"90 days SoLID","P");
  l1->AddEntry(g7,"60 days SBS","P");
  l1->SetTextSize(0.04);
  l1->SetTextFont(22);
  l1->Draw();
  
 // TCanvas *c2 = new TCanvas("c2","c2",1000,800);
 // c2->SetFillColor(10);
  //c2->Divide(1,2);
  c1->cd(2);
  //gPad->SetLeftMargin(0.16);
  //gPad->SetRightMargin(0.002);
  //gPad->SetTopMargin(0.005);
  //gPad->SetBottomMargin(0.10);
   plot_neutron_sbs(2,2);//(target,particle,Q2)
  plot_pretz_evolved(2,2);//(target, particle)
  l1->Draw();


  c1->Print("SoLID_neutron_both_pretz_SBS.pdf");
  c1->Print("SoLID_neutron_both_pretz_SBS.png");
}
/*}}}*/
