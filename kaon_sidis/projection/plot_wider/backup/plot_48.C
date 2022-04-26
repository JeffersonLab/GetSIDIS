void plot_48(int target_flag, int particle_flag){
   
	TCanvas *c1 = new TCanvas("c1","c1", 1000,1000);
	c1->Divide(8,6);
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
		return -1;
	}
	TString particle = "X";
	if(particle_flag==1)
		particle ="pip";
	else if(particle_flag==2)
		particle ="pim";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
		return -1;
	}
	/*}}}*/
	TString filename;
			
	TH2F *h1[8][6];
	TGraphErrors *g1[8][6];
	TLatex *t1[8][6];
    Int_t z_flag = -1, Q2_flag = -1;
	for(int i=0;i<8;i++){
		for(int j=0;j<1;j++){
			z_flag = i; Q2_flag = j;

			/*Get Z and Q2 Bin{{{*/
			double zmin =0., zmax = 0., Q2min=0., Q2max=0.;
			if (z_flag==1){
				zmin = 0.3; zmax = 0.35;
			}else if (z_flag==2){
				zmin = 0.35; zmax = 0.4;
			}else if (z_flag==3){
				zmin = 0.4; zmax = 0.45;
			}else if (z_flag==4){
				zmin = 0.45; zmax = 0.5;
			}else if (z_flag==5){
				zmin = 0.5; zmax = 0.55;
			}else if (z_flag==6){
				zmin = 0.55; zmax = 0.6;
			}else if (z_flag==7){
				zmin = 0.6; zmax = 0.65;
			}else if (z_flag==8){
				zmin = 0.65; zmax = 0.7;
			}
			if (Q2_flag==1){
				Q2min = 1.; Q2max = 2.0;
			}else if (Q2_flag==2){
				Q2min = 2.0; Q2max = 3.0;
			}else if (Q2_flag==3){
				Q2min = 3.0; Q2max = 4.0;
			}else if (Q2_flag==4){
				Q2min = 4.0; Q2max = 5.0;
			}else if (Q2_flag==5){
				Q2min = 5.0; Q2max = 6.0;
			}else if (Q2_flag==6){
				Q2min = 6.0; Q2max = 8.0;
			}else if (Q2_flag==7){//Not in used!!!
				Q2min = 8.0; Q2max = 10.0;
			}
			/*}}}*/

			filename.Form("./results/%s_%s_%d_%d.dat",target.Data(),particle.Data(),z_flag,Q2_flag);
			ifstream infile(filename);

			gStyle->SetOptStat(0);
			Int_t count1,count2;
			Double_t Q2[5000],x[5000],z[5000],pt[5000],y[5000],asy[5000],coverage[5000],coef[3][5000];
			Double_t Q2[5000],x[5000];
			Double_t temp;
			Int_t pt_flag = -1;
			Int_t ncount=0;

			infile >> count1;
			for (Int_t k=0;k!=count1;k++){
				infile >> pt_flag >> count2;
				if(pt_flag!=k)
					cerr<<"*** ERROR, something is wrong?"<<endl;

				for (Int_t l=0;l!=count2;l++){
					infile >> temp >> temp >> z[ncount] >> Q2[ncount] >> pt[ncount] >> x[ncount] >> y[ncount] >> asy[ncount] >> 
						coverage[ncount] >> coef[0][ncount] >> coef[1][ncount] >> coef[2][ncount];
					if (asy[ncount]>0.&&asy[ncount]<0.05&&coverage[ncount]>6.28/15.&&(coef[0][ncount]+coef[1][ncount]+coef[2][ncount])/3.<3.){ //&&asy[ncount]<0.05
						// 	asy[ncount] *= (coef[0][ncount]+coef[1][ncount]+coef[2][ncount])/3.; 
						asy[ncount] *= coef[0][ncount]; 
					}
					ncount ++;
				}
			}
			c1->cd((i+1)*(j+1));
			c1->SetGrid(1,1);
			h1[i][j] = new TH2F(Form("h1_%d_%d",i,j),"", 1000, 0.05,0.65, 1000, 0., 1.6);
			h1[i][j]->SetXTitle("x");
			h1[i][j]->GetXaxis()->CenterTitle(1);
			h1[i][j]->GetXaxis()->SetTitleFont(32);
			h1[i][j]->GetXaxis()->SetTitleSize(0.07);
			h1[i][j]->GetXaxis()->SetTitleOffset(0.7);
			h1[i][j]->GetXaxis()->SetNdivisions(510);
			h1[i][j]->SetYTitle("pt (GeV)");
			h1[i][j]->GetYaxis()->CenterTitle(1);
			h1[i][j]->GetYaxis()->SetTitleFont(32);
			h1[i][j]->GetYaxis()->SetTitleSize(0.06);
			h1[i][j]->GetYaxis()->SetTitleOffset(0.7);
			h1[i][j]->GetYaxis()->SetNdivisions(510);
			h1[i][j]->Draw();

			g1[i][j] = new TGraphErrors(ncount,x,pt,0,asy);
			g1[i][j]->SetMarkerStyle(20);
			g1[i][j]->SetMarkerColor(2);
			g1[i][j]->Draw("p");

			t1[i][j] = new TLatex();
			t1[i][j]->SetNDC();
			TString titlename;
			titlename.Form("%3.2f < z < %3.2f",zmin,zmax);
			t1[i][j]->DrawLatex(0.5,0.75,titlename);
			titlename.Form("%2.1f < Q^{2} < %2.1f (GeV^{2})",Q2min,Q2max);
			t1[i][j]->DrawLatex(0.5,0.85,titlename);
           
			c1->Update();	
			infile.close();

		}
	}

    c1->Print(Form("./figures/%s_%s_all48.png",target.Data(),particle.Data()));
    c1->Print(Form("./figures/%s_%s_all48.pdf",target.Data(),particle.Data()));
}


