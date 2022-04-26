#define max1 1.25;
#define min1 -0.1;
#define qmax1 0.35;
#define qmin1 -0.55;
#define xmin1 0.03
#define xmax1 0.51

/*void plot(Int_t z_flag =1, Int_t Q2_flag=1, Int_t flag_t=0, Int_t flag=0){{{*/
void plot(Int_t i=1, Int_t j=1, Int_t k =1, Int_t flag_t=0, Int_t flag=0){
	//flag_t: 1->collins, 2->sivers
	//flag: 1->ploting Asym axis

	int target_flag = 3;
	int particle_flag =i;
	int Q2_flag = j;
    int z_flag = k;

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
	if(particle_flag==1)
		particle ="pip";
	else if(particle_flag==2)
		particle ="pim";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}
	TString filename;
	filename.Form("./new_results/%s_%s_%d_%d.dat",target.Data(),particle.Data(),z_flag,Q2_flag);
	ifstream infile(filename);

	Float_t factor,shift;
//	factor=(max1-min1)/(qmax1-qmin1);
//	shift=max1/factor-qmax1;
	factor=(qmax1-qmin1)/(max1-min1);
	shift=qmax1/factor-max1;


	Int_t Q2min,Q2max;
	if (Q2_flag==1){
		Q2min = 1; Q2max =2;
	}else if (Q2_flag==2){
		Q2min = 2; Q2max =3;
	}else if (Q2_flag==3){
		Q2min = 3; Q2max =4;
	}else if (Q2_flag==4){
		Q2min = 4; Q2max =5;
	}else if (Q2_flag==5){
		Q2min = 5; Q2max =6;
	}else if (Q2_flag==6){
		Q2min = 6; Q2max =8;
	}

	Float_t zmin,zmax;
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
	
	/*}}}*/

	gStyle->SetOptStat(0);
	Int_t count1,count2;
	Double_t Q2[5000],x[5000],z[5000],pt[5000],y[5000],Astat[5000],N[5000],coverage[5000],coef[3][5000];
	Double_t Q2[5000],x[5000];
	Double_t temp;
	Int_t pt_flag = -1;
	Int_t ncount=0;
    TString pt_label[20];
    Float_t pt_pos[20];

	infile >> count1;
	Int_t pt_temp=0,x_temp=0;
	Float_t pt_min =0, pt_max =0;
	int pt_count = 0;
	while(infile >> pt_flag >> count2){
		//cerr<<" pt_flag="<<pt_flag<<", count="<<count2<<endl;

		for (Int_t j1=0;j1<count2;j1++){

			/*pt_flag{{{*/
			if (pt_flag==1){
				pt_min = 0.0; pt_max = 0.2;
			}else if (pt_flag==2){
				pt_min = 0.2; pt_max = 0.4;
			}else if (pt_flag==3){
				pt_min = 0.4; pt_max = 0.6;
			}else if (pt_flag==4){
				pt_min = 0.6; pt_max = 0.8;
			}else if (pt_flag==5){
				pt_min = 0.8; pt_max = 1.0;
			}else if (pt_flag==6){
				pt_min = 1.0; pt_max = 1.2;
			}else if (pt_flag==7){
				pt_min = 1.2; pt_max = 1.4;
			}else if (pt_flag==8){
				pt_min = 1.4; pt_max = 1.6;
			}
			/*}}}*/	
            if(pt_flag!=pt_temp){					
				pt_label[pt_count] = Form("<P_{T}>=%2.1f (GeV)", (pt_min+pt_max)*0.5);
				pt_pos[pt_count] = ((pt_min+pt_max)*0.5+shift)* factor; 
				pt_count++;	
			}

			infile >> pt_temp >> x_temp >> z[ncount] >> Q2[ncount] >> pt[ncount] >> x[ncount] >> y[ncount] >> Astat[ncount] >> N[ncount]>>
				coverage[ncount] >> coef[0][ncount] >> coef[1][ncount] >> coef[2][ncount];
			if(pt_temp==pt_flag && x_temp== j1){

				Astat[ncount] *= coef[2-flag_t][ncount]; //coef[0]->Sivers,coef[1]->Collins, coef[2]->Pretzelosity
				cerr<<Form("---z=%d, Q=%d pt = %f, x = %f, Asys= %f",z_flag,Q2_flag, pt[ncount],x[ncount], Astat[ncount])<<endl;

				//if (Astat[ncount]>0.&&Astat[ncount]<0.05&&coverage[ncount]>6.28/15.&&(coef[0][ncount]+coef[1][ncount]+coef[2][ncount])/3.<3.){ //
				if (Astat[ncount]>0.&&Astat[ncount]<0.05){
					// 	Astat[ncount] *= (coef[0][ncount]+coef[1][ncount]+coef[2][ncount])/3.;
				    pt[ncount] = (pt[ncount]+shift)* factor; 
					Astat[ncount] *= factor;
			 
					ncount ++;
				}
			}
		}
	}

		pt_label[pt_count] = Form("<P_{T}>=%2.1f (GeV)", (pt_min+pt_max)*0.5+0.2);
		pt_pos[pt_count] = ((pt_min+pt_max)*0.5+0.2+shift)* factor; 
		pt_count++;	


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
  if(particle_flag==1)
	  g2->GetYaxis()->SetTitle("A_{UT}^{Collins}(#pi^{+})");
  if(particle_flag==2)
	  g2->GetYaxis()->SetTitle("A_{UT}^{Collins}(#pi^{-})");
  g2->GetXaxis()->SetNdivisions(506);
  g2->GetYaxis()->SetNdivisions(506);
  g2->GetXaxis()->SetTitleOffset(0.8);
  g2->GetYaxis()->SetTitleOffset(0.8);
  g2->GetYaxis()->SetTitleSize(0.06);
  g2->GetXaxis()->SetTitleSize(0.07);
  g2->SetTitle(0);
  g2->Draw("AP");

  TGraphErrors *g1 = new TGraphErrors(ncount,x,pt,0,Astat);
  g1->Draw("*same");
  g1->SetMarkerColor(2);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(1.2);
  g1->SetTitle();

  titlename.Form("#scale[1.0]{%d < Q^{2} < %d}",Q2min,Q2max);
  t1->DrawLatex(0.75,0.35,titlename);
  titlename.Form("#scale[1.0]{%3.2f < z < %3.2f}",zmin,zmax);
  t1->DrawLatex(0.7,0.25,titlename);
  for(int i=0;i<pt_count;i++){
	  titlename.Form("#scale[0.8]{%s}",pt_label[i].Data());
	  t1->SetTextColor(1);
	  t1->SetTextFont(22);
	  t1->DrawLatex(0.18,0.28+0.322*(pt_pos[i]-pt_pos[0])/(pt_pos[pt_count]-pt_pos[0]),titlename);
  }

  TLine *bb = new TLine(xmin1,0,xmax1+0.006,0);
  bb->SetLineColor(12);
  bb->SetLineWidth(1.0);
  bb->Draw("same");

	}
	/*}}}*/

/*void plot_e06010(Int_t part=1, Int_t flag = 1, Int_t flag_6=0){{{*/
void plot_e06010(Int_t part=1, Int_t flag = 1, Int_t flag_6=0){
  // part = 1 for pi^+  2 for pi^-
  // flag = 1 for Collins 2 for Sivers
  
  Float_t factor,shift;
  //factor=(max1-min1)/(qmax1-qmin1);
  //shift=max1/factor-qmax1;
  factor = 1; shift =0;
  
  // Float_t x[4]={0.136689,     0.18919,        0.249406,       0.336276};
//   Float_t Collins_pip[4]={-0.0521008,        0.121314,       -0.0508186,     -0.279587};
//   Float_t Collins_pip_err[4]={0.113244,  0.107413,       0.103535,       0.102613};
//   Float_t Collins_pip_sys_err[4]={0.0640138,         0.0891116,      0.050367,       0.0855846};
  

//   Float_t Collins_pim[4]={-0.171281,         0.0155945,      0.042814,       0.0037289};
//   Float_t Collins_pim_err[4]={0.083695,  0.0783713,      0.072444,       0.0671274};
//   Float_t Collins_pim_sys_err[4]={0.0472002,         0.0468717,      0.0481302,      0.0399914};
  

//   Float_t Sivers_pip[4]={-0.010324,  -0.0867324,     -0.140877,      -0.11851};
//   Float_t Sivers_pip_err[4]={0.11159,    0.105726,       0.101663,       0.100828};
//   Float_t Sivers_pip_sys_err[4]={0.0867716,  0.143801,       0.0617867,      0.101248};
 

//   Float_t Sivers_pim[4]={-0.0755118,         0.0483624,      0.0658459,      0.0631764};
//   Float_t Sivers_pim_err[4]={ 0.0818826,  0.0768086,      0.0708994,      0.0656095};
//   Float_t Sivers_pim_sys_err[4]={0.0561873,  0.0649786,      0.0689067,      0.0452567};

  Float_t x[4]={(0.136926 + 0.136453) / 2+0.0204311, 
		 (0.189271 + 0.189109) / 2+0.018464,
		 (0.249588 + 0.249223) / 2+0.0165212, 
		 (0.33684 + 0.335713) / 2 +0.0140351};
  Float_t Collins_pip[4]={-0.0498879,0.146446,-0.0497506,-0.347633};
  Float_t Collins_pip_err[4]={0.128319,0.14281,0.0969608,0.15327};
  Float_t Collins_pip_sys_err[4]={0.138937,0.159033,0.106373,0.194805};
  

  Float_t Collins_pim[4]={-0.136839,0.0359655,0.0388251,0.0146006};
  Float_t Collins_pim_err[4]={0.0835516,0.1142,0.0667023,0.0761643};
  Float_t Collins_pim_sys_err[4]={0.0939522,0.128002,0.0748053,0.0853765};
  

  Float_t Sivers_pip[4]={-0.0111463,-0.108013,-0.129014,-0.162731};
  Float_t Sivers_pip_err[4]={0.125489,0.138209,0.0964687,0.13306};
  Float_t Sivers_pip_sys_err[4]={0.141468,0.158861,0.111911,0.158112};
 

  Float_t Sivers_pim[4]={-0.0589043,0.0796177,0.0586047,0.087065};
  Float_t Sivers_pim_err[4]={0.0820362,0.114184,0.0665064,0.0766304};
  Float_t Sivers_pim_sys_err[4]={0.0959498,0.134894,0.0785442,0.090822};
  
  
  Float_t xtemp[4],xtemp1[4],ytemp[4],ytemp1[4],yerrtemp[4],yerrtemp1[4],ytottemp[4],ytottemp1[4];
  
  xtemp[0]=x[0];
  xtemp[1]=x[1];
  xtemp1[0]=x[2];
  xtemp1[1]=x[3];
  
  if (part==1 && flag ==1){ //pi^+ Collins
    ytemp[0]=(Collins_pip[0]+shift)*factor;
    ytemp[1]=(Collins_pip[1]+shift)*factor;
    ytemp1[0]=(Collins_pip[2]+shift)*factor;
    ytemp1[1]=(Collins_pip[3]+shift)*factor;
    
    yerrtemp[0]=Collins_pip_err[0]*factor;
    yerrtemp[1]=Collins_pip_err[1]*factor;
    yerrtemp1[0]=Collins_pip_err[2]*factor;
    yerrtemp1[1]=Collins_pip_err[3]*factor;

    ytottemp[0]=sqrt(Collins_pip_err[0]*Collins_pip_err[0]*0.+Collins_pip_sys_err[0]*Collins_pip_sys_err[0])*factor;
    ytottemp[1]=sqrt(Collins_pip_err[1]*Collins_pip_err[1]*0.+Collins_pip_sys_err[1]*Collins_pip_sys_err[1])*factor;
    ytottemp1[0]=sqrt(Collins_pip_err[2]*Collins_pip_err[2]*0.+Collins_pip_sys_err[2]*Collins_pip_sys_err[2])*factor;
    ytottemp1[1]=sqrt(Collins_pip_err[3]*Collins_pip_err[3]*0.+Collins_pip_sys_err[3]*Collins_pip_sys_err[3])*factor;
  }else if (part==2 && flag ==1){ //pi^- Collins
    ytemp[0]=(Collins_pim[0]+shift)*factor;
    ytemp[1]=(Collins_pim[1]+shift)*factor;
    ytemp1[0]=(Collins_pim[2]+shift)*factor;
    ytemp1[1]=(Collins_pim[3]+shift)*factor;
    
    yerrtemp[0]=Collins_pim_err[0]*factor;
    yerrtemp[1]=Collins_pim_err[1]*factor;
    yerrtemp1[0]=Collins_pim_err[2]*factor;
    yerrtemp1[1]=Collins_pim_err[3]*factor;

    ytottemp[0]=sqrt(Collins_pim_err[0]*Collins_pim_err[0]*0.+Collins_pim_sys_err[0]*Collins_pim_sys_err[0])*factor;
    ytottemp[1]=sqrt(Collins_pim_err[1]*Collins_pim_err[1]*0.+Collins_pim_sys_err[1]*Collins_pim_sys_err[1])*factor;
    ytottemp1[0]=sqrt(Collins_pim_err[2]*Collins_pim_err[2]*0.+Collins_pim_sys_err[2]*Collins_pim_sys_err[2])*factor;
    ytottemp1[1]=sqrt(Collins_pim_err[3]*Collins_pim_err[3]*0.+Collins_pim_sys_err[3]*Collins_pim_sys_err[3])*factor;
  }else if (part==1 && flag ==2){ //pi^+ Sivers
    ytemp[0]=(Sivers_pip[0]+shift)*factor;
    ytemp[1]=(Sivers_pip[1]+shift)*factor;
    ytemp1[0]=(Sivers_pip[2]+shift)*factor;
    ytemp1[1]=(Sivers_pip[3]+shift)*factor;
    
    yerrtemp[0]=Sivers_pip_err[0]*factor;
    yerrtemp[1]=Sivers_pip_err[1]*factor;
    yerrtemp1[0]=Sivers_pip_err[2]*factor;
    yerrtemp1[1]=Sivers_pip_err[3]*factor;

    ytottemp[0]=sqrt(Sivers_pip_err[0]*Sivers_pip_err[0]*0.+Sivers_pip_sys_err[0]*Sivers_pip_sys_err[0])*factor;
    ytottemp[1]=sqrt(Sivers_pip_err[1]*Sivers_pip_err[1]*0.+Sivers_pip_sys_err[1]*Sivers_pip_sys_err[1])*factor;
    ytottemp1[0]=sqrt(Sivers_pip_err[2]*Sivers_pip_err[2]*0.+Sivers_pip_sys_err[2]*Sivers_pip_sys_err[2])*factor;
    ytottemp1[1]=sqrt(Sivers_pip_err[3]*Sivers_pip_err[3]*0.+Sivers_pip_sys_err[3]*Sivers_pip_sys_err[3])*factor;
  }else if (part==2 && flag ==2){ //pi^- Sivers
   ytemp[0]=(Sivers_pim[0]+shift)*factor;
    ytemp[1]=(Sivers_pim[1]+shift)*factor;
    ytemp1[0]=(Sivers_pim[2]+shift)*factor;
    ytemp1[1]=(Sivers_pim[3]+shift)*factor;
    
    yerrtemp[0]=Sivers_pim_err[0]*factor;
    yerrtemp[1]=Sivers_pim_err[1]*factor;
    yerrtemp1[0]=Sivers_pim_err[2]*factor;
    yerrtemp1[1]=Sivers_pim_err[3]*factor;

    ytottemp[0]=sqrt(Sivers_pim_err[0]*Sivers_pim_err[0]*0.+Sivers_pim_sys_err[0]*Sivers_pim_sys_err[0])*factor;
    ytottemp[1]=sqrt(Sivers_pim_err[1]*Sivers_pim_err[1]*0.+Sivers_pim_sys_err[1]*Sivers_pim_sys_err[1])*factor;
    ytottemp1[0]=sqrt(Sivers_pim_err[2]*Sivers_pim_err[2]*0.+Sivers_pim_sys_err[2]*Sivers_pim_sys_err[2])*factor;
    ytottemp1[1]=sqrt(Sivers_pim_err[3]*Sivers_pim_err[3]*0.+Sivers_pim_sys_err[3]*Sivers_pim_sys_err[3])*factor;
  }


  
  TGraphErrors *trans = new TGraphErrors(2,xtemp,ytemp,0,yerrtemp);
  trans->SetMarkerColor(1);
  trans->SetMarkerSize(1.5);  
  trans->SetMarkerStyle(21);
  TGraphErrors *trans2 = new TGraphErrors(2,xtemp,ytemp,0,ytottemp);
  trans2->SetMarkerColor(1);
  trans2->SetMarkerSize(1.5);  
  trans2->SetMarkerStyle(21);

  
  TGraphErrors *trans1 = new TGraphErrors(2,xtemp1,ytemp1,0,yerrtemp1);
  trans1->SetMarkerColor(1);
  trans1->SetMarkerSize(1.5);  
  trans1->SetMarkerStyle(21);
  TGraphErrors *trans3 = new TGraphErrors(2,xtemp1,ytemp1,0,ytottemp1);
  trans3->SetMarkerColor(1);
  trans3->SetMarkerSize(1.5);  
  trans3->SetMarkerStyle(21);
  
  if (flag_6==1){
    trans->Draw("Psame");
    trans2->Draw("Psame");
  }else if (flag_6==2){
    trans1->Draw("Psame");
    trans3->Draw("Psame");

  }else if (flag_6==3){
    trans->Draw("Psame");
    trans2->Draw("Psame");
    trans1->Draw("Psame");
    trans3->Draw("Psame");
  }
  
}
/*}}}*/

/*void plot_feng(Int_t part=1,Int_t flag_t = 0,Int_t flag_z = 0){{{*/
void plot_feng(Int_t part=1,Int_t flag_t = 0,Int_t flag_z = 0){
  
  Float_t factor,shift;
  //factor=(max1-min1)/(qmax1-qmin1);
  //shift=max1/factor-qmax1;
  factor = 1; shift =0;
  TString partname;
  TString FCollinsfile;
  TString FSiversfile;
  ifstream infile;

  if (part==1){
    partname = "pip_";
    FCollinsfile="./calculations/FENG/Collins_pip.dat";
    FSiversfile="./calculations/FENG/Sivers_pip.dat";
  }else if (part==2){
    partname = "pim_";
    FCollinsfile="./calculations/FENG/Collins_pim.dat";
    FSiversfile="./calculations/FENG/Sivers_pim.dat";
  }
  Float_t FCollinx[4][5], FCollinz[4][5], FCollinAn[4][5],FCollinAnh[4][5],FCollinAnl[4][5],FCollinAp[4][5],FCollinAph[4][5],FCollinApl[4][5];
    Float_t FSiverx[4][5], FSiverz[4][5], FSiverAn[4][5],FSiverAnh[4][5],FSiverAnl[4][5],FSiverAp[4][5],FSiverAph[4][5],FSiverApl[4][5];
    Float_t qqtemp[4][5];
    
    Float_t FCollinbandx[4][10],FCollinbandyp[4][10],FCollinbandyn[4][10];
    Float_t FSiverbandx[4][10],FSiverbandyp[4][10],FSiverbandyn[4][10];
    
    infile.open(FCollinsfile,ios::in);
    for (Int_t j=0;j!=5;j++){
      for (Int_t i=0;i!=4;i++){
	infile >> FCollinx[i][j] >> FCollinz[i][j] >> FCollinAn[i][j] >> FCollinAp[i][j] >> FCollinAnl[i][j] >> FCollinApl[i][j] >> FCollinAnh[i][j] >> FCollinAph[i][j];
	qqtemp[i][j]=0.;
	if (FCollinApl[i][j]>FCollinAp[i][j]){
	Float_t tmp;
	tmp=FCollinApl[i][j];
	FCollinApl[i][j]=FCollinAph[i][j];
	FCollinAph[i][j]=tmp;
      }
      if (FCollinAnl[i][j]>FCollinAn[i][j]){
	tmp=FCollinAnl[i][j];
	FCollinAnl[i][j]=FCollinAnh[i][j];
	FCollinAnh[i][j]=tmp;
      }
      
      FCollinbandx[i][j]=FCollinx[i][j];
      FCollinbandx[i][9-j]=FCollinx[i][j];
      FCollinbandyp[i][j]=(FCollinApl[i][j]+shift)*factor;
      FCollinbandyp[i][9-j]=(FCollinAph[i][j]+shift)*factor;
      FCollinbandyn[i][j]=(FCollinAnl[i][j]+shift)*factor;
      FCollinbandyn[i][9-j]=(FCollinAnh[i][j]+shift)*factor;
      


      FCollinApl[i][j]=(FCollinAp[i][j]-FCollinApl[i][j])*factor;
      FCollinAph[i][j]=(FCollinAph[i][j]-FCollinAp[i][j])*factor;
      FCollinAnl[i][j]=(FCollinAn[i][j]-FCollinAnl[i][j])*factor;
      FCollinAnh[i][j]=(FCollinAnh[i][j]-FCollinAn[i][j])*factor;


      FCollinAn[i][j]=(FCollinAn[i][j]+shift)*factor;
      }
    }
    infile.close();
    
    infile.open(FSiversfile,ios::in);
    for (Int_t j=0;j!=5;j++){
      for (Int_t i=0;i!=4;i++){
	
	infile >> FSiverx[i][j] >> FSiverz[i][j] >> FSiverAn[i][j] >> FSiverAp[i][j] >> FSiverAnl[i][j] >> FSiverApl[i][j] >> FSiverAnh[i][j] >> FSiverAph[i][j];
      if (FSiverApl[i][j]>FSiverAp[i][j]){
	Float_t tmp;
	tmp=FSiverApl[i][j];
	FSiverApl[i][j]=FSiverAph[i][j];
	FSiverAph[i][j]=tmp;
      }
      if (FSiverAnl[i][j]>FSiverAn[i][j]){
	tmp=FSiverAnl[i][j];
	FSiverAnl[i][j]=FSiverAnh[i][j];
	FSiverAnh[i][j]=tmp;
      }
      FSiverbandx[i][j]=FSiverx[i][j];
      FSiverbandx[i][9-j]=FSiverx[i][j];
      FSiverbandyp[i][j]=(FSiverApl[i][j]+shift)*factor;
      FSiverbandyp[i][9-j]=(FSiverAph[i][j]+shift)*factor;
      FSiverbandyn[i][j]=(FSiverAnl[i][j]+shift)*factor;
      FSiverbandyn[i][9-j]=(FSiverAnh[i][j]+shift)*factor;

      FSiverApl[i][j]=(FSiverAp[i][j]-FSiverApl[i][j])*factor;
      FSiverAph[i][j]=(FSiverAph[i][j]-FSiverAp[i][j])*factor;
      FSiverAnl[i][j]=(FSiverAn[i][j]-FSiverAnl[i][j])*factor;
      FSiverAnh[i][j]=(FSiverAnh[i][j]-FSiverAn[i][j])*factor;

      FSiverAn[i][j]=(FSiverAn[i][j]+shift)*factor;
      }
    }
    infile.close();
    TGraph** Fcollinasy=new TGraph*[8];
    TGraph** Fsiverasy=new TGraph*[8];
    TGraph** Fcollinbandn=new TGraph*[8];
    TGraph** Fsiverbandn=new TGraph*[8];
    for (Int_t i=0;i!=8;i++){
      Fcollinasy[i] = new TGraph(5,&FCollinx[Int_t(i/2.+0.1)][0],&FCollinAn[Int_t(i/2.+0.1)][0]);
      Fsiverasy[i] = new TGraph(5,&FSiverx[Int_t(i/2.+0.1)][0],&FSiverAn[Int_t(i/2.+0.1)][0]);
      
      Fcollinbandn[i] = new TGraph(10,&FCollinbandx[Int_t(i/2.+0.1)][0],&FCollinbandyn[Int_t(i/2.+0.1)][0]);
      Fsiverbandn[i] = new TGraph(10,&FSiverbandx[Int_t(i/2.+0.1)][0],&FSiverbandyn[Int_t(i/2.+0.1)][0]);
      
      Fcollinasy[i]->SetLineColor(3);
      Fsiverasy[i]->SetLineColor(3);
      Fcollinasy[i]->SetLineStyle(1);
      Fsiverasy[i]->SetLineStyle(2);
      Fcollinasy[i]->SetLineWidth(2.);
      Fsiverasy[i]->SetLineWidth(2.);
      
      Fcollinbandn[i]->SetFillColor(3);
      Fsiverbandn[i]->SetFillColor(3);
      Fcollinbandn[i]->SetFillStyle(0);
      Fsiverbandn[i]->SetFillStyle(0);      
    }

    if (flag_t==1){
      Fcollinasy[flag_z-1]->Draw("SL");
      Fcollinbandn[flag_z-1]->Draw("SLF");
    }else{
      Fsiverasy[flag_z-1]->Draw("SL");
      Fsiverbandn[flag_z-1]->Draw("SLF");
    }
    
    
}
/*}}}*/

/*void plot_ase(Int_t part=1,Int_t flag_t = 0){{{*/
void plot_ase(Int_t part=1,Int_t flag_t = 0){
  
  Float_t factor,shift;
  //factor=(max1-min1)/(qmax1-qmin1);
  //shift=max1/factor-qmax1;
  factor = 1; shift =0;
  Float_t alex_collins_pip_x[15]={0.15       ,0.2035714 ,0.2571429 ,0.3107143,0.3642857 ,0.4178571 ,0.4714286 ,0.525    ,0.5785714 ,0.6321429 ,0.6857143 ,0.7392857 ,0.7928571 ,0.8464286,0.9};
    Float_t alex_collins_pip_A[15]={0.004082549,-0.01631688,-0.03070491,-0.04104227,-0.04851713,-0.05346644,-0.05744219,-0.06021951,-0.06195057,-0.06273717,-0.06255476,-0.06130253,-0.05867733,-0.05404292,-0.04606274};
    
    Float_t alex_collins_pip_shadex[15]={0.15,0.2035714,0.2571429,0.3107143,0.3642857,0.4178571,0.4714286,0.4714286,0.4178571,0.3642857,0.3107143,0.2571429,0.2035714,0.15,0.15};
    Float_t alex_collins_pip_shadey[15]={0.00352671,-0.01836635,-0.03348654,-0.04391603,-0.05207097,-0.05840553,-0.06448088,-0.04508678,-0.04571461,-0.04387626,-0.03652231,-0.02715431,-0.01418629,0.004734356,0.00352671};

    Float_t alex_collins_pim_x[15]={0.15       ,0.2035714 ,0.2571429 ,0.3107143,0.3642857 ,0.4178571 ,0.4714286 ,0.525    ,0.5785714 ,0.6321429 ,0.6857143 ,0.7392857 ,0.7928571 ,0.8464286,0.9};
    Float_t alex_collins_pim_A[15]={-0.00468036,0.0193133,0.0371514,0.0502582,0.05945248,0.06414484,0.067623,0.06905009,0.06868634,0.06670207,0.06316258,0.05809289,0.05136418,0.04274524,0.03186207};

    Float_t alex_collins_pim_shadex[15]={0.15,0.2035714,0.2571429,0.3107143,0.3642857,0.4178571,0.4714286,0.4714286,0.4178571,0.3642857,0.3107143,0.2571429,0.2035714,0.15,0.15};
    Float_t alex_collins_pim_shadey[15]={-0.005230773,0.01714451,0.03309532,0.0452816,0.05320285,0.05303781,0.05114383,0.07745728,0.07162967,0.06496163,0.05407908,0.0404383,0.02109671,-0.004108014,-0.005230773};

    Float_t alex_sivers_pip_x[15]={0.01,0.0592857,0.108571,0.157857,0.207143,0.256429,0.305714,,0.355,0.404286,0.453571,0.502857,0.552143,0.601429,0.650714,0.7};
    Float_t alex_sivers_pip_A[15]={-0.0154218,-0.0914365,-0.150187,-0.186967,-0.204018,-0.205179,-0.194398,-0.175677,-0.151408,-0.125047,-0.0986611,-0.0740528,-0.0524703,-0.0346685,-0.0209667};

    Float_t alex_sivers_pip_shadex[21]={0.01,0.0592857,0.108571,0.157857,0.207143,0.256429,0.305714,0.355,0.404286,0.453571,0.453571,0.404286,0.355,0.305714,0.256429,0.207143,0.157857,0.108571,0.0592857,0.01,0.01};					
    Float_t alex_sivers_pip_shadey[21]={-0.0417382,-0.10967,-0.183874,-0.228258,-0.239553,-0.245266,-0.246228,-0.239004,-0.226254,-0.215257
					,-0.0240066,-0.0386121
					,-0.0558229,
					-0.0759046,
					-0.0966106,
					-0.113218,
					-0.111764,
					-0.0997433,
					-0.0638783,
					-0.00492083,
					-0.0417382};

    Float_t alex_sivers_pim_x[15]={0.01,0.0592857,0.108571,0.157857,0.207143,0.256429,0.305714,,0.355,0.404286,0.453571,0.502857,0.552143,0.601429,0.650714,0.7};
    Float_t alex_sivers_pim_A[15]={-0.0112223,-0.0425134,-0.0582766,-0.0627498,-0.0592402,-0.0509382,-0.0404207,-0.0302808,-0.0240403,-0.0166755,-0.0106578,-0.00605562,-0.00279292,-0.000718029,0.000368377};

    Float_t alex_sivers_pim_shadex[21]={0.01,0.0592857,0.108571,0.157857,0.207143,0.256429,0.305714,0.355,0.404286,0.453571,0.453571,0.404286,0.355,0.305714,0.256429,0.207143,0.157857,0.108571,0.0592857,0.01,0.01};
    Float_t alex_sivers_pim_shadey[21]={-0.0204835,-0.056709,-0.0763758,-0.0837231,-0.081482,-0.0735104,-0.0655461,-0.0562907,-0.0501532,-0.0407899,0.00702034,0.00321284,0.00086072,-0.00623796,-0.0142115,-0.0215203,-0.0269707,-0.027734,-0.0218104,-0.00381826,-0.0204835};
    

    for (Int_t qxqx=0;qxqx!=15;qxqx++){
      alex_collins_pip_A[qxqx] = (alex_collins_pip_A[qxqx]+shift)*factor;
      alex_collins_pim_A[qxqx] = (alex_collins_pim_A[qxqx]+shift)*factor;
      
      alex_sivers_pip_A[qxqx] = (alex_sivers_pip_A[qxqx]+shift)*factor;
      alex_sivers_pim_A[qxqx] = (alex_sivers_pim_A[qxqx]+shift)*factor;
      //cout << alex_collins_pip_A[i] << endl;
    }
    for (Int_t qxqx=0;qxqx!=15;qxqx++){
      alex_collins_pip_shadey[qxqx] = (alex_collins_pip_shadey[qxqx]+shift)*factor;
      alex_collins_pim_shadey[qxqx] = (alex_collins_pim_shadey[qxqx]+shift)*factor;
    }
    for (Int_t qxqx=0;qxqx!=21;qxqx++){
      alex_sivers_pip_shadey[qxqx] = (alex_sivers_pip_shadey[qxqx]+shift)*factor;
      alex_sivers_pim_shadey[qxqx] = (alex_sivers_pim_shadey[qxqx]+shift)*factor;
    }
    
    if (part==1){
      if (flag_t==1){
	TGraph *alex1 = new TGraph(7,alex_collins_pip_x,alex_collins_pip_A);
	alex1->SetLineColor(6);
	alex1->SetLineWidth(2.);
	alex1->SetLineStyle(1);
	alex1->Draw("SL");
      
	TGraph *alex1_band = new TGraph(15,alex_collins_pip_shadex,alex_collins_pip_shadey);
	alex1_band->SetFillColor(6);
	alex1_band->SetFillStyle(0);
	alex1_band->Draw("SLF");
      }else{
	TGraph *alex2 = new TGraph(10,alex_sivers_pip_x,alex_sivers_pip_A);
	alex2->SetLineColor(6);
	alex2->SetLineWidth(2.);
	alex2->SetLineStyle(2);
	alex2->Draw("SL");
	
	TGraph *alex2_band = new TGraph(21,alex_sivers_pip_shadex,alex_sivers_pip_shadey);
	alex2_band->SetFillColor(6);
	alex2_band->SetFillStyle(0);
	alex2_band->Draw("SLF");
      }    
    }else if(part==2){
      if (flag_t==1){
	TGraph *alex1 = new TGraph(7,alex_collins_pim_x,alex_collins_pim_A);
	alex1->SetLineColor(6);
	alex1->SetLineWidth(2.);
	alex1->SetLineStyle(1);
	alex1->Draw("SL");
	
	TGraph *alex1_band = new TGraph(15,alex_collins_pim_shadex,alex_collins_pim_shadey);
	alex1_band->SetFillColor(6);
	alex1_band->SetFillStyle(0);
	alex1_band->Draw("SLF");
      }else{
	TGraph *alex2 = new TGraph(10,alex_sivers_pim_x,alex_sivers_pim_A);
	alex2->SetLineColor(6);
	alex2->SetLineWidth(2.);
	alex2->SetLineStyle(2);
	alex2->Draw("SL");
	
	TGraph *alex2_band = new TGraph(21,alex_sivers_pim_shadex,alex_sivers_pim_shadey);
	alex2_band->SetFillColor(6);
	alex2_band->SetFillStyle(0);
	alex2_band->Draw("SLF");
      }
    }
}
/*}}}*/

/*void plot_pretz(Int_t part=1,Int_t flag_t = 0){{{*/
void plot_pretz(Int_t part=1,Int_t flag_t = 0){
  Float_t factor,shift;
  //factor=(max1-min1)/(qmax1-qmin1);
  //shift=max1/factor-qmax1;
  factor = 1; shift =0;
  
  ifstream infile3("./calculations/aut_sin_phi_plus_phis.dat");
  Float_t x_full[100],ppip_full[100],ppim_full[100],ppi0_full[100];
  Float_t dpip_full[100],dpim_full[100],dpi0_full[100];
  Float_t npip_full[100],npim_full[100],npi0_full[100];
  
  Int_t i_full=0,i_pwave=0,i_swave=0;
  while(!infile3.eof()){
    infile3 >> x_full[i_full] >> ppip_full[i_full] >> ppim_full[i_full] >> ppi0_full[i_full] >> dpip_full[i_full] >> dpim_full[i_full] >> dpi0_full[i_full] >> npip_full[i_full] >> npim_full[i_full] ;

    npip_full[i_full]=(npip_full[i_full]+shift)*factor;
    npim_full[i_full]=(npim_full[i_full]+shift)*factor;
    
    i_full++;
  }
  infile3.close();
  
  if (part==1){
    TGraph *qqabc = new TGraph(i_full-1,x_full,npip_full);
  }else{
    TGraph *qqabc = new TGraph(i_full-1,x_full,npim_full);
  }
  qqabc->SetLineColor(4);
  qqabc->SetLineWidth(2.5);

  ifstream infile4("./calculations/aut_sin_3phi_minus_phis.dat");
  Float_t x_pwave[100],ppip_pwave[100],ppim_pwave[100],ppi0_pwave[100];
  Float_t dpip_pwave[100],dpim_pwave[100],dpi0_pwave[100];
  Float_t npip_pwave[100],npim_pwave[100],npi0_pwave[100];
  while(!infile4.eof()){
    infile4 >> x_pwave[i_pwave] >> ppip_pwave[i_pwave] >> ppim_pwave[i_pwave] >> ppi0_pwave[i_pwave] >> dpip_pwave[i_pwave] >> dpim_pwave[i_pwave] >> dpi0_pwave[i_pwave] >> npip_pwave[i_pwave] >> npim_pwave[i_pwave] >> npi0_pwave[i_pwave] ;
    npip_pwave[i_pwave]=(npip_pwave[i_pwave]+shift)*factor;
    npim_pwave[i_pwave]=(npim_pwave[i_pwave]+shift)*factor;
    i_pwave++;
  }
  infile4.close();
  if (part==1){
    TGraph *qqabc1 = new TGraph(i_pwave-1,x_pwave,npip_pwave);
  }else{
    TGraph *qqabc1 = new TGraph(i_pwave-1,x_pwave,npim_pwave);
  }
  qqabc1->SetLineColor(7);
  qqabc1->SetLineWidth(2.5);

  if (flag_t==1){
    qqabc->Draw("SL");
  }else{
    qqabc1->Draw("SL");
  }
}
/*}}}*/

/*void plot_pretz_evolved(Int_t part=1,Int_t flag_t = 0){{{*/
void plot_pretz_evolved(Int_t part=1,Int_t flag_t = 0){
  Float_t factor,shift;
  //factor=(max1-min1)/(qmax1-qmin1);
  //shift=max1/factor-qmax1;
  factor = 1; shift =0;
  
  ifstream infile3("./calculations/BARBARA/aut_sin_phi_plus_phis_evolved.dat");
  Float_t x_full[100],ppip_full[100],ppim_full[100],ppi0_full[100];
  Float_t dpip_full[100],dpim_full[100],dpi0_full[100];
  Float_t npip_full[100],npim_full[100],npi0_full[100];
  
  Int_t i_full=0,i_pwave=0,i_swave=0;
  while(!infile3.eof()){
    infile3 >> x_full[i_full] >> ppip_full[i_full] >> ppim_full[i_full] >> ppi0_full[i_full] >> dpip_full[i_full] >> dpim_full[i_full] >> dpi0_full[i_full] >> npip_full[i_full] >> npim_full[i_full] >> npi0_full[i_full];

    npip_full[i_full]=(npip_full[i_full]+shift)*factor;
    npim_full[i_full]=(npim_full[i_full]+shift)*factor;
    
    i_full++;
  }
  infile3.close();
  
  if (part==1){
    TGraph *qqabc = new TGraph(i_full-1,x_full,npip_full);
  }else{
    TGraph *qqabc = new TGraph(i_full-1,x_full,npim_full);
  }
  qqabc->SetLineColor(4);
  qqabc->SetLineWidth(2.5);

  ifstream infile4("./calculations/aut_sin_3phi_minus_phis.dat");
  Float_t x_pwave[100],ppip_pwave[100],ppim_pwave[100],ppi0_pwave[100];
  Float_t dpip_pwave[100],dpim_pwave[100],dpi0_pwave[100];
  Float_t npip_pwave[100],npim_pwave[100],npi0_pwave[100];
  while(!infile4.eof()){
    infile4 >> x_pwave[i_pwave] >> ppip_pwave[i_pwave] >> ppim_pwave[i_pwave] >> ppi0_pwave[i_pwave] >> dpip_pwave[i_pwave] >> dpim_pwave[i_pwave] >> dpi0_pwave[i_pwave] >> npip_pwave[i_pwave] >> npim_pwave[i_pwave] >> npi0_pwave[i_pwave] ;
    npip_pwave[i_pwave]=(npip_pwave[i_pwave]+shift)*factor;
    npim_pwave[i_pwave]=(npim_pwave[i_pwave]+shift)*factor;
    i_pwave++;
  }
  infile4.close();
  if (part==1){
    TGraph *qqabc1 = new TGraph(i_pwave-1,x_pwave,npip_pwave);
  }else{
    TGraph *qqabc1 = new TGraph(i_pwave-1,x_pwave,npim_pwave);
  }
  qqabc1->SetLineColor(7);
  qqabc1->SetLineWidth(2.5);

  if (flag_t==1){
    qqabc->Draw("SL");
  }else{
    qqabc1->Draw("SL");
  }
}
/*}}}*/
