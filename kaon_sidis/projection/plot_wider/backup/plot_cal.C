void plot_cal(){

  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"Y");
  //gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"X"); 

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->SetFillColor(10);
  c1->Divide(2,2);
  c1->cd(1);
  plot_cal1(1);
 
  
  c1->cd(2);
  plot_cal1(2);
  
  c1->cd(3);
  plot_cal2(1);
  
  c1->cd(4);
  plot_cal2(2);

  
}

void plot_cal1(Int_t part=1){
    Float_t max1,min1,qmax1,qmin1;
  Float_t factor,shift;
  max1=0.7;min1=0.3;
  qmax1=0.4;qmin1=-0.6;
  factor=(max1-min1)/(qmax1-qmin1);
  shift=max1/factor-qmax1;
  
     Float_t alex_collins_pip_x[15]={0.01,0.08785714,0.1657143,0.2435714,0.3214286,0.3992857,0.4771429,0.555,0.6328571,0.7107143,0.7885714,0.8664286,0.9442857,1.022143,1.1};
    Float_t alex_collins_pip_A[15]={-0.0009071399,-0.007933148,-0.01478533,-0.02131613,-0.02738959,-0.0328861,-0.03770622,-0.04177358,-0.04503673,-0.04746977,-0.04907183,-0.04986565,-0.0498951,-0.04922211,-0.04792298};
    
    Float_t alex_collins_pip_shadex[30]={0.01,0.08785714,0.1657143,0.2435714,0.3214286,0.3992857,0.4771429,0.555,0.6328571,0.7107143,0.7885714,0.8664286,0.9442857,1.022143,1.1,1.1,1.022143,0.9442857,0.8664286,0.7885714,0.7107143,0.6328571,0.555,0.4771429,0.3992857,0.3214286,0.2435714,0.1657143,0.08785714,0.01};
    Float_t alex_collins_pip_shadey[30]={-0.001041205,-0.009095646,-0.01690406,-0.02425975,-0.03097618,-0.03690433,-0.04248016,-0.047281,-0.05124685,-0.05434208,-0.05705538,-0.05925904,-0.0608223,-0.06178749,-0.0620909,-0.03297334,-0.03577295,-0.03814901,-0.03994966,-0.04040578,-0.04006875,-0.03886992,-0.03614254,-0.03262004,-0.02844758,-0.02369118,-0.01843676,-0.01278761,-0.006861083,-0.0007840761};


    Float_t alex_collins_pim_x[15]={0.01,0.08785714,0.1657143,0.2435714,0.3214286,0.3992857,0.4771429,0.555,0.6328571,0.7107143,0.7885714,0.8664286,0.9442857,1.022143,1.1};
    Float_t alex_collins_pim_A[15]={0.001068012,0.009340011,0.01740736,0.02509632,0.03224686,0.03871812,0.04439303,0.0491817,0.05302354,0.05588805,0.05777422,0.05870881,0.05874349,0.05795115,0.05642163};

    Float_t alex_collins_pim_shadex[30]={0.01,0.08785714,0.1657143,0.2435714,0.3214286,0.3992857,0.4771429,0.555,0.6328571,0.7107143,0.7885714,0.8664286,0.9442857,1.022143,1.1,1.1,1.022143,0.9442857,0.8664286,0.7885714,0.7107143,0.6328571,0.555,0.4771429,0.3992857,0.3214286,0.2435714,0.1657143,0.08785714,0.01};
    Float_t alex_collins_pim_shadey[30]={0.0009300273,0.008137899,0.01518915,0.0219501,0.02829633,0.03411643,0.03931524,0.04373647,0.04715864,0.04891261,0.04924554,0.04795063,0.04578935,0.04293743,0.03957712,0.07085276,0.0704396,0.06933926,0.06745215,0.06498033,0.06183637,0.05817901,0.0536767,0.04871619,0.04296508,0.03611363,0.02830843,0.01973786,0.01062468,0.001216424};


    Float_t alex_sivers_pip_x[15]={0.01,0.0807143,0.151429,0.222143,0.292857,0.363571,0.434286,0.505,0.575714,0.646429,0.717143,0.787857,0.858571,0.929286,1};
    Float_t alex_sivers_pip_A[15]={-0.00535578,-0.0431111,-0.0803187,-0.116511,-0.15124,-0.184091,-0.214687,-0.242695,-0.267837,-0.289889,-0.308688,-0.324129,-0.336169,-0.344822,-0.350158};

    Float_t alex_sivers_pip_shadex[30]={0.01,0.0807143,0.151429,0.222143,0.292857,0.363571,0.434286,0.505,0.575714,0.646429,0.717143,0.787857,0.858571,0.929286,1,1,0.929286,0.858571,0.787857,0.717143,0.646429,0.575714,0.505,0.434286,0.363571,0.292857,0.222143,0.151429,0.0807143,0.01};					
    Float_t alex_sivers_pip_shadey[30]={-0.00617688,-0.0497386,-0.0927525,-0.134749,-0.175277,-0.213912,-0.250264,-0.283983,-0.314766,-0.342358,-0.36656,-0.387778,-0.406387,-0.421582,-0.433357,-0.199838,-0.196131,-0.190613,-0.183258,-0.17407,-0.16308,-0.150353,-0.135983,-0.120092,-0.102834,-0.0843867,-0.0649504,-0.0447455,-0.0240075,-0.00298202};

    Float_t alex_sivers_pim_x[15]={0.01,0.0807143,0.151429,0.222143,0.292857,0.363571,0.434286,0.505,0.575714,0.646429,0.717143,0.787857,0.858571,0.929286,1};
    Float_t alex_sivers_pim_A[15]={-0.00144856,-0.0116601,-0.0217235,-0.0315121,-0.0409053,-0.0497904,-0.0580654,-0.0656407,-0.0724408,-0.0784052,-0.0834896,-0.0876659,-0.0909222,-0.0932625,-0.0947057};

    Float_t alex_sivers_pim_shadex[30]={0.01,0.0807143,0.151429,0.222143,0.292857,0.363571,0.434286,0.505,0.575714,0.646429,0.717143,0.787857,0.858571,0.929286,1,1,0.929286,0.858571,0.787857,0.717143,0.646429,0.575714,0.505,0.434286,0.363571,0.292857,0.222143,0.151429,0.0807143,0.01};
    Float_t alex_sivers_pim_shadey[30]={-0.00202669,-0.0163148,-0.0304003,-0.04411,-0.0572784,-0.0697512,-0.0813879,-0.0920648,-0.101677,-0.110141,-0.117392,-0.123391,-0.129123,-0.133878,-0.137536,-0.0350713,-0.034252,-0.0331372,-0.0317254,-0.0300195,-0.0280268,-0.0257593,-0.0232332,-0.0204692,-0.0174921,-0.0143302,-0.0110152,-0.00758134,-0.00406523,-0.000504836};
    

    for (Int_t qxqx=0;qxqx!=15;qxqx++){
      alex_collins_pip_A[qxqx] = (alex_collins_pip_A[qxqx]+shift)*factor;
      alex_collins_pim_A[qxqx] = (alex_collins_pim_A[qxqx]+shift)*factor;
      
      alex_sivers_pip_A[qxqx] = (alex_sivers_pip_A[qxqx]+shift)*factor;
      alex_sivers_pim_A[qxqx] = (alex_sivers_pim_A[qxqx]+shift)*factor;
      //cout << alex_collins_pip_A[i] << endl;
    }
    for (Int_t qxqx=0;qxqx!=30;qxqx++){
      alex_collins_pip_shadey[qxqx] = (alex_collins_pip_shadey[qxqx]+shift)*factor;
      alex_collins_pim_shadey[qxqx] = (alex_collins_pim_shadey[qxqx]+shift)*factor;
    }
    for (Int_t qxqx=0;qxqx!=30;qxqx++){
      alex_sivers_pip_shadey[qxqx] = (alex_sivers_pip_shadey[qxqx]+shift)*factor;
      alex_sivers_pim_shadey[qxqx] = (alex_sivers_pim_shadey[qxqx]+shift)*factor;
    }

    TLine *a = new TLine(0.0,shift*factor,1.2,shift*factor);

    if (part==1){
      TGraph *alex1 = new TGraph(15,&alex_collins_pip_x[0],&alex_collins_pip_A[0]);
      alex1->SetLineColor(6);
      alex1->SetLineWidth(2.);
      alex1->SetLineStyle(1);
      alex1->GetYaxis()->SetRangeUser(min1,max1);
      alex1->GetXaxis()->SetTitle("P_{T} (GeV/c)");
      alex1->GetYaxis()->SetTitle("z");
      alex1->SetTitle("");
      
      alex1->Draw("AL");
      
      TGraph *alex1_band = new TGraph(30,alex_collins_pip_shadex,alex_collins_pip_shadey);
      alex1_band->SetFillColor(6);
      alex1_band->SetFillStyle(3001);
      alex1_band->Draw("SLF");

      
      TGraph *alex2 = new TGraph(15,&alex_sivers_pip_x[0],&alex_sivers_pip_A[0]);
      alex2->SetLineColor(6);
      alex2->SetLineWidth(2.);
      alex2->SetLineStyle(2);
      alex2->Draw("SL");
      
      TGraph *alex2_band = new TGraph(30,alex_sivers_pip_shadex,alex_sivers_pip_shadey);
      alex2_band->SetFillColor(6);
      alex2_band->SetFillStyle(3007);
      alex2_band->Draw("SLF");
    }else if(part==2){
      TGraph *alex1 = new TGraph(15,&alex_collins_pim_x[0],&alex_collins_pim_A[0]);
      alex1->SetLineColor(6);
      alex1->SetLineWidth(2.);
      alex1->SetLineStyle(1);
      alex1->GetYaxis()->SetRangeUser(min1,max1);
      alex1->Draw("AL");
      alex1->GetXaxis()->SetTitle("P_{T} (GeV/c)");
      alex1->GetYaxis()->SetTitle("z");
      alex1->SetTitle("");
      TGraph *alex1_band = new TGraph(30,alex_collins_pim_shadex,alex_collins_pim_shadey);
      alex1_band->SetFillColor(6);
      alex1_band->SetFillStyle(3001);
      alex1_band->Draw("SLF");
      
      TGraph *alex2 = new TGraph(15,&alex_sivers_pim_x[0],&alex_sivers_pim_A[0]);
      alex2->SetLineColor(6);
      alex2->SetLineWidth(2.);
      alex2->SetLineStyle(2);
      alex2->Draw("SL");
      
      TGraph *alex2_band = new TGraph(30,alex_sivers_pim_shadex,alex_sivers_pim_shadey);
      alex2_band->SetFillColor(6);
      alex2_band->SetFillStyle(3007);
      alex2_band->Draw("SLF");
    }

    a->Draw("same");
    TGaxis *axis = new TGaxis(alex1->GetXaxis()->GetXmax(),min1,alex1->GetXaxis()->GetXmax(),max1,qmin1,qmax1,510,"+L");
    if (part==1){
    axis->SetTitle("Asymmetry (#pi^{+})");
  }else if (part==2){
    axis->SetTitle("Asymmetry (#pi^{-})");
  }
    axis->Draw();
}


void plot_cal2(Int_t part=1){
  TString partname;
  TString FCollinsfile;
  TString FSiversfile;
  if (part==1){
    partname = "pip_";
    FCollinsfile="/w/halla/e03004/xin/12gev/TRANSVERSITY/calculations/FENG/Collins_pip.dat";
    FSiversfile="/w/halla/e03004/xin/12gev/TRANSVERSITY/calculations/FENG/Sivers_pip.dat";
  }else if (part==2){
    partname = "pim_";
    FCollinsfile="/w/halla/e03004/xin/12gev/TRANSVERSITY/calculations/FENG/Collins_pim.dat";
    FSiversfile="/w/halla/e03004/xin/12gev/TRANSVERSITY/calculations/FENG/Sivers_pim.dat";
  }
  
  ifstream infile;
  Float_t max1,min1,qmax1,qmin1;
  Float_t factor,shift;
  max1=1.2;min1=-0.2;
  qmax1=0.4;qmin1=-0.6;
  factor=(max1-min1)/(qmax1-qmin1);
  shift=max1/factor-qmax1;
  if (part==1||part==2){
    //Feng's calculations
    
    Float_t FCollinx[4][5], FCollinz[4][5], FCollinAn[4][5],FCollinAnh[4][5],FCollinAnl[4][5],FCollinAp[4][5],FCollinAph[4][5],FCollinApl[4][5];
    Float_t FSiverx[4][5], FSiverz[4][5], FSiverAn[4][5],FSiverAnh[4][5],FSiverAnl[4][5],FSiverAp[4][5],FSiverAph[4][5],FSiverApl[4][5];
    Float_t qqtemp[4][5];
    

    Float_t FCx[5][4],FCA[5][4],FSx[5][4],FSA[5][4];
    Float_t FCbx[5][8],FCbyp[5][8],FCbyn[5][8];
    Float_t FSbx[5][8],FSbyp[5][8],FSbyn[5][8];
    // Float_t FCollinbandx[4][10],FCollinbandyp[4][10],FCollinbandyn[4][10];
//     Float_t FSiverbandx[4][10],FSiverbandyp[4][10],FSiverbandyn[4][10];

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
      
      FCx[j][i] = FCollinz[i][j];
      FCA[j][i] = (FCollinAn[i][j]+shift)*factor;
      
      FCbx[j][i] = FCollinz[i][j];
      FCbx[j][7-i] = FCollinz[i][j];
      FCbyn[j][i] = (FCollinAnl[i][j]+shift)*factor;
      FCbyn[j][7-i] = (FCollinAnh[i][j]+shift)*factor;
      
 
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
     
      FSx[j][i] = FSiverz[i][j];
      FSA[j][i] = (FSiverAn[i][j]+shift)*factor;
      
      FSbx[j][i] = FSiverz[i][j];
      FSbx[j][7-i] = FSiverz[i][j];
      FSbyn[j][i] = (FSiverAnl[i][j]+shift)*factor;
      FSbyn[j][7-i] = (FSiverAnh[i][j]+shift)*factor;

      }
    }
    infile.close();
  }


  TGraph** Fcollinasy=new TGraph*[8];
  TGraph** Fsiverasy=new TGraph*[8];
  TGraph** Fcollinbandn=new TGraph*[8];
  TGraph** Fsiverbandn=new TGraph*[8];
  
  if (part==1||part==2){
    for (Int_t i=0;i!=8;i++){
      
      Fcollinasy[i] = new TGraph(4,&FCx[Int_t(i/2.+0.1)][0],&FCA[Int_t(i/2.+0.1)][0]);
      Fsiverasy[i] = new TGraph(4,&FSx[Int_t(i/2.+0.1)][0],&FSA[Int_t(i/2.+0.1)][0]);
      
      Fcollinbandn[i] = new TGraph(8,&FCbx[Int_t(i/2.+0.1)][0],&FCbyn[Int_t(i/2.+0.1)][0]);
      Fsiverbandn[i] = new TGraph(8,&FSbx[Int_t(i/2.+0.1)][0],&FSbyn[Int_t(i/2.+0.1)][0]);


      Fcollinasy[i]->SetLineColor(4);
      Fsiverasy[i]->SetLineColor(4);
      Fcollinasy[i]->SetLineStyle(1);
      Fsiverasy[i]->SetLineStyle(2);
      Fcollinasy[i]->SetLineWidth(2.);
      Fsiverasy[i]->SetLineWidth(2.);
      
      Fcollinbandn[i]->SetFillColor(4);
      Fsiverbandn[i]->SetFillColor(4);
      Fcollinbandn[i]->SetFillStyle(3001);
      Fsiverbandn[i]->SetFillStyle(3007);      


      //mgx[i]->Add(Fcollinasy[i]);
      //mgx[i]->Add(Fsiverasy[i]);
    }
    
  }
  

  TLine *a = new TLine(0.325,shift*factor,0.675,shift*factor);
  
  
    
    Fcollinasy[4]->Draw("AL");
    Fcollinasy[4]->GetYaxis()->SetRangeUser(min1,max1);
    Fcollinasy[4]->GetYaxis()->SetTitle("z");
    Fcollinasy[4]->GetXaxis()->SetTitle("x");
    Fcollinasy[4]->SetTitle("");
    
    Fsiverasy[4]->Draw("SL");
   

    Fcollinbandn[4]->Draw("SLF");
    Fsiverbandn[4]->Draw("SLF");
 
    TGaxis *axis = new TGaxis(Fcollinasy[4]->GetXaxis()->GetXmax(),min1,Fcollinasy[4]->GetXaxis()->GetXmax(),max1,qmin1,qmax1,510,"+L");
    if (part==1){
      axis->SetTitle("Asymmetry (#pi^{+})");
    }else if (part==2){
      axis->SetTitle("Asymmetry (#pi^{-})");
    }
    axis->Draw();
    a->Draw("same");


      Float_t alex_collins_pip_x[15]={0.1,0.1535714,0.2071429,0.2607143,0.3142857,0.3678571,0.4214286,0.475,0.5285714,0.5821429,0.6357143,0.6892857,0.7428571,0.7964286,0.85};
    Float_t alex_collins_pip_A[15]={-0.009797974,-0.01405856,-0.017792,-0.02113237,-0.02418771,-0.0270193,-0.02963044,-0.03198722,-0.03402365,-0.03569542,-0.03697892,-0.0378952,-0.03848495,-0.03880436,-0.03889044};
    
    Float_t alex_collins_pip_shadex[16]={0.3142857,0.3678571,0.4214286,0.475,0.5285714,0.5821429,0.6357143,0.6892857,0.6892857,0.6357143,0.5821429,0.5285714,0.475,0.4214286,0.3678571,0.3142857};
    Float_t alex_collins_pip_shadey[16]={-0.01105096,-0.01574964,-0.01999276,-0.02382328,-0.02733502,-0.03058051,-0.03356507,-0.03626354,-0.02769452,-0.02569908,-0.02348717,-0.02109419,-0.01851593,-0.01568564,-0.0124865,-0.00870852};

    Float_t alex_collins_pim_x[15]={0.1,0.1535714,0.2071429,0.2607143,0.3142857,0.3678571,0.4214286,0.475,0.5285714,0.5821429,0.6357143,0.6892857,0.7428571,0.7964286,0.85};
    Float_t alex_collins_pim_A[15]={0.01162085,0.01674229,0.02120608,0.02515162,0.02871348,0.03199167,0.03499749,0.03769569,0.03999106,0.04179464,0.04306359,0.0437946,0.04403291,0.04380168,0.04320618};

    Float_t alex_collins_pim_shadex[16]={0.3142857,0.3678571,0.4214286,0.475,0.5285714,0.5821429,0.6357143,0.6892857,0.6892857,0.6357143,0.5821429,0.5285714,0.475,0.4214286,0.3678571,0.3142857};
    Float_t alex_collins_pim_shadey[16]={0.009863652,0.01433008,0.01831446,0.02190201,0.02517484,0.0281946,0.03096816,0.03347408,0.04126722,0.03842056,0.03524082,0.03176774,0.02798429,0.02376431,0.0189341,0.01345915};

    Float_t alex_sivers_pip_x[15]={0.1,0.142857,0.185714,0.228571,0.271429,0.314286,0.357143,0.4,0.442857,0.485714,0.528571,0.571429,0.614286,0.657143,0.7};
    Float_t alex_sivers_pip_A[15]={-0.0332081,-0.0506665,-0.0695722,-0.0893083,-0.109296,-0.129017,-0.148151,-0.166584,-0.184394,-0.201888,-0.219373,-0.237083,-0.25516,-0.273582,-0.292074};

    Float_t alex_sivers_pip_shadex[20]={0.314286,0.357143,0.4,0.442857,0.485714,0.528571,0.571429,0.614286,0.657143,0.7,0.7,0.657143,0.614286,0.571429,0.528571,0.485714,0.442857,0.4,0.357143,0.314286};					
    Float_t alex_sivers_pip_shadey[20]={-0.152484,-0.174935,-0.196512,-0.217308,-0.237682,-0.257992,-0.278514,-0.299412,-0.32066,-0.34194,-0.162965,-0.152883,-0.142805,-0.13288,-0.123116,-0.113435,-0.103721,-0.0938117,-0.0835439,-0.0728779};

    Float_t alex_sivers_pim_x[15]={0.1,0.142857,0.185714,0.228571,0.271429,0.314286,0.357143,0.4,0.442857,0.485714,0.528571,0.571429,0.614286,0.657143,0.7};
    Float_t alex_sivers_pim_A[15]={-0.02631,-0.0339051,-0.0393613,-0.0431591,-0.0458633,-0.048037,-0.0501003,-0.0521866,-0.0541562,-0.0555858,-0.0559942,-0.0549027,-0.0520631,-0.0474668,-0.0413786};

    Float_t alex_sivers_pim_shadex[20]={0.314286,0.357143,0.4,0.442857,0.485714,0.528571,0.571429,0.614286,0.657143,0.7,0.7,0.657143,0.614286,0.571429,0.528571,0.485714,0.442857,0.4,0.357143,0.314286};
    Float_t alex_sivers_pim_shadey[20]={-0.0625295,-0.066583,-0.0705937,-0.0744239,-0.07766911,-0.0798686,-0.0805612,-0.0795127,-0.0767184,-0.0724283,-0.00271443,-0.00946943,-0.0143556,-0.0172728,-0.0191164,-0.020042,-0.0203554,-0.0203642,-0.0203391,-0.0203678};
    

    for (Int_t qxqx=0;qxqx!=15;qxqx++){
      alex_collins_pip_A[qxqx] = (alex_collins_pip_A[qxqx]+shift)*factor;
      alex_collins_pim_A[qxqx] = (alex_collins_pim_A[qxqx]+shift)*factor;
      
      alex_sivers_pip_A[qxqx] = (alex_sivers_pip_A[qxqx]+shift)*factor;
      alex_sivers_pim_A[qxqx] = (alex_sivers_pim_A[qxqx]+shift)*factor;
      //cout << alex_collins_pip_A[i] << endl;
    }
    for (Int_t qxqx=0;qxqx!=16;qxqx++){
      alex_collins_pip_shadey[qxqx] = (alex_collins_pip_shadey[qxqx]+shift)*factor;
      alex_collins_pim_shadey[qxqx] = (alex_collins_pim_shadey[qxqx]+shift)*factor;
    }
    for (Int_t qxqx=0;qxqx!=20;qxqx++){
      alex_sivers_pip_shadey[qxqx] = (alex_sivers_pip_shadey[qxqx]+shift)*factor;
      alex_sivers_pim_shadey[qxqx] = (alex_sivers_pim_shadey[qxqx]+shift)*factor;
    }
    

    if (part==1){
      TGraph *alex1 = new TGraph(8,&alex_collins_pip_x[4],&alex_collins_pip_A[4]);
      alex1->SetLineColor(6);
      alex1->SetLineWidth(2.);
      alex1->SetLineStyle(1);
      alex1->Draw("SL");
      
      TGraph *alex1_band = new TGraph(15,alex_collins_pip_shadex,alex_collins_pip_shadey);
      alex1_band->SetFillColor(6);
      alex1_band->SetFillStyle(3001);
      alex1_band->Draw("SLF");

      
      TGraph *alex2 = new TGraph(10,&alex_sivers_pip_x[5],&alex_sivers_pip_A[5]);
      alex2->SetLineColor(6);
      alex2->SetLineWidth(2.);
      alex2->SetLineStyle(2);
      alex2->Draw("SL");
      
      TGraph *alex2_band = new TGraph(20,alex_sivers_pip_shadex,alex_sivers_pip_shadey);
      alex2_band->SetFillColor(6);
      alex2_band->SetFillStyle(3007);
      alex2_band->Draw("SLF");
    }else if(part==2){
      TGraph *alex1 = new TGraph(8,&alex_collins_pim_x[4],&alex_collins_pim_A[4]);
      alex1->SetLineColor(6);
      alex1->SetLineWidth(2.);
      alex1->SetLineStyle(1);
      alex1->Draw("SL");
      
      TGraph *alex1_band = new TGraph(15,alex_collins_pim_shadex,alex_collins_pim_shadey);
      alex1_band->SetFillColor(6);
      alex1_band->SetFillStyle(3001);
      alex1_band->Draw("SLF");
      
      TGraph *alex2 = new TGraph(10,&alex_sivers_pim_x[5],&alex_sivers_pim_A[5]);
      alex2->SetLineColor(6);
      alex2->SetLineWidth(2.);
      alex2->SetLineStyle(2);
      alex2->Draw("SL");
      
      TGraph *alex2_band = new TGraph(20,alex_sivers_pim_shadex,alex_sivers_pim_shadey);
      alex2_band->SetFillColor(6);
      alex2_band->SetFillStyle(3007);
      alex2_band->Draw("SLF");

      TLegend *l1 = new TLegend(0.6,0.6,0.9,0.9);
      l1->SetHeader("calculation");
      l1->AddEntry(alex1,"Anselmino et al.(Collins)","L");
      l1->AddEntry(alex2,"Anselmino et al.(Sivers)","L");
      l1->AddEntry(Fcollinasy[4],"Vogelsang and Yuan (Collins)","L");
      l1->AddEntry(Fsiverasy[4],"Vogelsang and Yuan (Sivers)","L");
      l1->Draw();
    }


}
