

void plotcomparison(int processfiles = 0, int cutchoice = 0)
{
  gROOT->Reset();
  gStyle->SetOptStat(0);
  //gROOT->SetStyle("BABAR");
  TChain *T0 = new TChain("T");

  const Double_t EIC_Mom_Min_e = 0.5;
  const Double_t EIC_Mom_Max_e = 3.*10.0; //not in use
  const Double_t EIC_Mom_Min_h = 0.0;
  const Double_t EIC_Mom_Max_h = 10.0;

  const Double_t EIC_Th_Min_e = 0.0;
  const Double_t EIC_Th_Max_e  = 140.0;
  const Double_t EIC_Th_Min_h = 0.0;
  const Double_t EIC_Th_Max_h  = 180.0;

  const Double_t EIC_Ph_Min_e = 0.0;
  const Double_t EIC_Ph_Max_e  = 360.0;
  const Double_t EIC_Ph_Min_h = 0.0;
  const Double_t EIC_Ph_Max_h  = 360.0;

  const Double_t degtorad = TMath::Pi()/180.;

  //T0->Add(Form("./c12_pion/EIC_A12_pion_10_600_%d_%d.root", j, i));
  //EIC:  1 -> Q2<10, pt<1
  //               2 -> Q2<10, pt>1
  //           EIC:  3 -> Q2>10, pt<1
  //                 4 -> Q2>10, pt>1
 TString outputfileending = "-";

 Double_t histo_ptmin = 0;
 Double_t histo_ptmax = 5;


 cout << "Which files do you want to process? Choose 1 for (Q2<10, pt<1), 2 for (Q2<10, pt>1), 3 for (Q2>10, pt<1), 4 for (Q2>10, pt>1), " ;
 cout << "a combination of the two numbers in order for any two files (12,13,14,23,24,34) and anything else for all" << endl;
 if (processfiles == 0) {
   cout << "Are you sure you want to process all files? If yes choose 0, otherwise one of the above mentioned options " << endl;
   cin >> processfiles;
 }
 cout << "processfiles " << processfiles << endl;
 if (processfiles == 1)
 {
   cout << "File with Q2<10 and pt<1 is processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_1_2.root");
   outputfileending += "f1-" ;
   histo_ptmax = 1.1;
 }
 else if (processfiles == 2)
 {
   cout << "File with Q2<10 and pt>1 is processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_2_2.root");
   outputfileending += "f2-" ;
 }
 else if (processfiles == 3)
 {
   cout << "File with Q2>10 and pt<1 is processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_3_2.root");
   outputfileending += "f3-" ;
   histo_ptmax = 1.1;
 }
 else if (processfiles == 4)
 {
   cout << "File with Q2>10 and pt>1 is processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_4_2.root");
   outputfileending += "f4-" ;
 }
 else if (processfiles == 12)
 {
   cout << "Files with Q2<10 are processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_1_2.root");
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_2_2.root");
   outputfileending += "f12-" ;
 }
 else if (processfiles == 13)
 {
   cout << "Files with pt<1 are processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_1_2.root");
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_3_2.root");
   outputfileending += "f13-" ;
   histo_ptmax = 1.1;
 }
 else if (processfiles == 14)
 {
   cout << "Files with (Q2<10,pt<1) and (Q2>10,pt>1) are processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_1_2.root");
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_4_2.root");
   outputfileending += "f14-" ;
 }
 else if (processfiles == 23)
 {
   cout << "Files with (Q2<10,pt>1) and (Q2>10,pt<1) are processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_2_2.root");
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_3_2.root");
   outputfileending += "f23-" ;
 }
 else if (processfiles == 24)
 {
   cout << "Files with pt>1 are processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_2_2.root");
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_4_2.root");
   outputfileending += "f24-" ;
 }
 else if (processfiles == 34)
 {
   cout << "Files with Q2>10 are processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_3_2.root");
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_4_2.root");
   outputfileending += "f34-" ;
 }
 else //all files are processed
 {
   cout << "all files are processed" << endl;
   T0->Add("../../test/Runs_for_MCPStests-Jul17/EIC_A12_pion_10_600_*_2.root");
   outputfileending += "fall_" ;
 }

 cout << "Total events " << T0->GetEntries() << endl;

//Choice of cut for all histograms

cout << "Which cut? Choose 0 for no cut, choose 1 for no x cut and (weight_hp>0), choose 2 for (x>0.0 && x<=0.3) and (weight_hp>0), 3 for (x>0.05 && x<=0.1) and (weight_hp>0), " ;
cout << "choose 11 for no xcut and weight_hp (weighting), choose 12 for (x>0.0 && x<=0.3) and weight_hp (weighting), choose 13 for (x>0.05 && x<=0.1) and weight_hp (weighting)" << endl;
if (cutchoice == 0) {
  cout << "Are you sure you want to process without any cut? If yes choose 0, otherwise one of the above mentioned options " << endl;
  cin >> cutchoice;
}
cout << "choosen cut " << cutchoice << endl;
TString xbj_cut_wide = "(x>0.0 && x<=0.3)";
TString xbj_cut_small = "(x>0.05 && x<=0.1)";
TString weightcut = "(weight_hp>0)";
TString weighting = "(weight_hp)";
TString allcut;
if (cutchoice == 0) {
  allcut = "";
  outputfileending += "c_none" ;
}
else if (cutchoice == 1) {
  allcut = weightcut;
  outputfileending += "c_nox_w0" ;
}
else if (cutchoice == 2) {
  allcut = xbj_cut_wide+"*"+weightcut;
  outputfileending += "c_xwi_w0" ;
}
else if (cutchoice == 3) {
  allcut = xbj_cut_small+"*"+weightcut;
  outputfileending += "c_xsm_w0" ;
}
else if (cutchoice == 11) {
  allcut = weighting;
  outputfileending += "c_nox_w1" ;
}
else if (cutchoice == 12) {
  allcut = xbj_cut_wide+"*"+weighting;
  outputfileending += "c_xwi_w1" ;
}
else if (cutchoice == 13) {
  allcut = xbj_cut_small+"*"+weighting;
  outputfileending += "c_xsm_w1" ;
}
else {
  cout << "Choice of cut not defined!!!" << endl;
  throw std::exception();
}
 //allcut += "*(1000000/nsim)";
 outputfileending += ".pdf" ;
 cout << "file name extension: " << outputfileending << endl;

 //**Define**
 TString savestring; //to store file name of saved canvas
 //Definitions for determination of min and max values in labframe spectra
 Double_t ele_pmax = EIC_Mom_Max_e;
 Double_t had_pmax = EIC_Mom_Max_h;
 Double_t ele_themin = EIC_Th_Min_e;
 Double_t ele_themax = EIC_Th_Max_e;
 Double_t had_themax = EIC_Th_Max_h;
  //Definitions of Tree values
 Double_t Q2, W, Wp, x, y, z, pt, nu, s, epsilon,rapidity, jacoF;
 Double_t theta_gen= 0.0 , phi_gen = 0.0, mom_gen = 0.0;
 Double_t mom_gen_ele,mom_gen_had;
 Double_t theta_gen_ele,theta_gen_had;
 Double_t phi_gen_ele,phi_gen_had;
 Double_t theta_q, theta_s,phi_h,phi_s,mom_ele,mom_had,theta_ele, theta_had,phi_ele,phi_had;
 Double_t dxs_incl,dxs_hm,dxs_hp,dilute_hp,dilute_hm;
 Double_t px_ele, py_ele,pz_ele, px_had, py_had, pz_had, E_ele,E_had;
 Double_t u_pdf, d_pdf, s_pdf, g_pdf, ubar_pdf, dbar_pdf, sbar_pdf;
 Double_t weight_hp, weight_hm, weight_in;
 ULong64_t nsim = 0;
 //For Beam Position and Vertex info
 Double_t vx_ele, vy_ele, vz_ele, vx_had, vy_had, vz_had;
 Double_t D_fav, D_unfav, D_s, D_g;

 T0->SetBranchAddress("Q2",&Q2);
 T0->SetBranchAddress("W",&W);
 T0->SetBranchAddress("Wp",&Wp);
 T0->SetBranchAddress("x",&x );
 T0->SetBranchAddress("y",&y );
 T0->SetBranchAddress("z",&z );
 T0->SetBranchAddress("nu",&nu );
 T0->SetBranchAddress("s",&s );
 T0->SetBranchAddress("epsilon",&epsilon );
 T0->SetBranchAddress("pt",&pt );
 T0->SetBranchAddress("weight_hp",&weight_hp );
 T0->SetBranchAddress("weight_hm",&weight_hm );
 T0->SetBranchAddress("weight_in",&weight_in );
 T0->SetBranchAddress("rapidity",&rapidity );
 T0->SetBranchAddress("theta_q",&theta_q );
 T0->SetBranchAddress("theta_s",&theta_s );
 T0->SetBranchAddress("phi_h",&phi_h );
 T0->SetBranchAddress("phi_s",&phi_s );
 T0->SetBranchAddress("jacoF",&jacoF);
 T0->SetBranchAddress("dxs_hm",&dxs_hm);
 T0->SetBranchAddress("dxs_hp",&dxs_hp);
 T0->SetBranchAddress("dxs_incl",&dxs_incl);
 T0->SetBranchAddress("mom_ele",&mom_ele);
 //T0->SetBranchAddress("mom_gen_ele",&mom_gen_ele);
 T0->SetBranchAddress("mom_had",&mom_had);
 //T0->SetBranchAddress("mom_gen_had",&mom_gen_had);
 T0->SetBranchAddress("theta_ele",&theta_ele);
 //T0->SetBranchAddress("theta_gen_ele",&theta_gen_ele);
 T0->SetBranchAddress("theta_had",&theta_had);
 //T0->SetBranchAddress("theta_gen_had",&theta_gen_had);
 T0->SetBranchAddress("phi_ele",&phi_ele);
 //T0->SetBranchAddress("phi_gen_ele",&phi_gen_ele);
 T0->SetBranchAddress("phi_had",&phi_had);
 //T0->SetBranchAddress("phi_gen_had",&phi_gen_had);
 T0->SetBranchAddress("nsim",&nsim);
 T0->SetBranchAddress("dilute_p",&dilute_hp );
 T0->SetBranchAddress("dilute_m",&dilute_hm  );
 T0->SetBranchAddress("px_ele",&px_ele);
 T0->SetBranchAddress("py_ele",&py_ele);
 T0->SetBranchAddress("pz_ele",&pz_ele);
 T0->SetBranchAddress("E_ele",&E_ele);
 T0->SetBranchAddress("px_had",&px_had);
 T0->SetBranchAddress("py_had",&py_had);
 T0->SetBranchAddress("pz_had",&pz_had);
 T0->SetBranchAddress("E_had",&E_had);
 T0->SetBranchAddress("vx_ele",&vx_ele);
 T0->SetBranchAddress("vy_ele",&vy_ele);
 T0->SetBranchAddress("vz_ele",&vz_ele);
 T0->SetBranchAddress("vx_had",&vx_had);
 T0->SetBranchAddress("vy_had",&vy_had);
 T0->SetBranchAddress("vz_had",&vz_had);

 T0->SetBranchAddress("u_pdf", &u_pdf);
 T0->SetBranchAddress("d_pdf", &d_pdf);
 T0->SetBranchAddress("s_pdf", &s_pdf);
 T0->SetBranchAddress("g_pdf", &g_pdf);
 T0->SetBranchAddress("ubar_pdf", &ubar_pdf);
 T0->SetBranchAddress("dbar_pdf", &dbar_pdf);
 T0->SetBranchAddress("sbar_pdf", &sbar_pdf);

 T0->SetBranchAddress("D_g", &D_g);
 T0->SetBranchAddress("D_s", &D_s);
 T0->SetBranchAddress("D_fav", &D_fav);
 T0->SetBranchAddress("D_unfav", &D_unfav);
  int  number_of_hminusevents = 0;
  int  number_of_hplusevents = 0;
  for (int i = 0; i<T0->GetEntries() ; i++)
  {
    T0->GetEntry(i);
    if (x>0.0 && x<=0.3){
      if (weight_hp > 0){
        number_of_hplusevents++;
      }
      if (weight_hm > 0) {
        number_of_hminusevents++;
      }

    }
  }
  cout << "For (x>0.0 && x<=0.3) and weight > 0: pos hadrons events: " << number_of_hplusevents << " and neg. hadron events: " << number_of_hminusevents  << endl;
  TCanvas *c1 = new TCanvas("c1","Kinematics 2d Histos", 1200, 800);
  c1->Divide(3,2);
  TH2D *h1 = new TH2D("h1","^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; p_{T} (GeV); z (GeV)", 100, histo_ptmin, histo_ptmax,100, 0.15, 0.85);
  h1->GetXaxis()->CenterTitle(1);
  h1->GetYaxis()->CenterTitle(1);
  TH2D *h4 = new TH2D("h4","^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; p_{T} (GeV); log10(Q^{2}) ", 100, histo_ptmin, histo_ptmax,100, -0.1, 3.0);
  h4->GetXaxis()->CenterTitle(1);
  h4->GetYaxis()->CenterTitle(1);
  c1->cd(1);
//  gPad->SetLogz(1);
  T0->Draw("z:pt>>h1",TCut(allcut),"colz");
  c1->cd(4);
//  gPad->SetLogz(1);
  T0->Draw("log10(Q2):pt>>h4",TCut(allcut),"colz");
//    T0->Draw("z:pt>>h1","weight_hm*(weight_hp>0.0)","colz");

  TH2D *h2 = new TH2D("h2","^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; log10(x); z (GeV)", 100,-3.5, 0.0,100, 0.15, 0.85);
  h2->GetXaxis()->CenterTitle(1);
  h2->GetYaxis()->CenterTitle(1);
  TH2D *h5 = new TH2D("h5","^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; log10(x); log10(Q^{2}) ", 100,-3.5, 0.0,100, -0.1, 3.0);
  h5->GetXaxis()->CenterTitle(1);
  h5->GetYaxis()->CenterTitle(1);
  c1->cd(2);
  //gPad->SetLogz(1);
  T0->Draw("z:log10(x)>>h2",TCut(allcut),"colz");
  c1->cd(5);
//  gPad->SetLogz(1);
  T0->Draw("log10(Q2):log10(x)>>h5",TCut(allcut),"colz");

  TH2D *h3 = new TH2D("h3","^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; log10(Q^{2}) ; z (GeV)", 100, -0.1, 3.0,100, 0.15, 0.85);
  h3->GetXaxis()->CenterTitle(1);
  h3->GetYaxis()->CenterTitle(1);
  TH2D *h6 = new TH2D("h6","^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; log10(x); p_{T} (GeV)", 100,-3.5, 0.0, 100, histo_ptmin, histo_ptmax);
  h6->GetXaxis()->CenterTitle(1);
  h6->GetYaxis()->CenterTitle(1);
  c1->cd(3);
  //gPad->SetLogz(1);
  T0->Draw("z:log10(Q2)>>h3",TCut(allcut),"colz");
  c1->cd(6);
//  gPad->SetLogz(1);
  T0->Draw("pt:log10(x)>>h6",TCut(allcut),"colz");
  savestring.Clear();
  savestring = "c12_pip_kinvalues2d" + outputfileending;
  c1->Print(savestring);
//  c1->Print("c12_pip_kinvalues2d_xbcut.pdf");




//1-dim histograms of kinemtic values
TCanvas *c2 = new TCanvas("c2","Kinematics 1d Histos", 1024, 800);
c2->Divide(2,2);

TH1D *h1pt = new TH1D("h1pt","pt distribution, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; p_{T} (GeV)", 100,histo_ptmin,histo_ptmax);
h1pt->GetXaxis()->CenterTitle(1);
TH1D *h1z = new TH1D("h1z","z distribution, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; z (GeV)", 100,0,1);
h1z->GetXaxis()->CenterTitle(1);
TH1D *h1x = new TH1D("h1x","logx distribution, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; log10(x)", 100,-4,0);
h1x->GetXaxis()->CenterTitle(1);
TH1D *h1Q2 = new TH1D("h1Q2","LogQ2 distribution, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; log10(Q^{2})", 60,0,3);
h1Q2->GetXaxis()->CenterTitle(1);

c2->cd(1);
T0->Draw("pt>>h1pt",TCut(allcut),"");
c2->cd(2);
T0->Draw("z>>h1z",TCut(allcut),"");
c2->cd(3);
T0->Draw("log10(x)>>h1x",TCut(allcut),"");
c2->cd(4);
T0->Draw("log10(Q2)>>h1Q2",TCut(allcut),"");
savestring.Clear();
savestring = "c12_pip_kinvalues1d" + outputfileending;
c2->Print(savestring);


//1-dim histograms of generated lab values
TCanvas *c4 = new TCanvas("c4","Generated Lab values 1d Histos", 1024, 800);
c4->Divide(2,2);

TH1D *h1p_ele = new TH1D("h1p_ele","Electron Momentum Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; P_{e} (GeV/c)", int((EIC_Mom_Max_e+1)*4), 0.0, EIC_Mom_Max_e+1);
h1p_ele->GetXaxis()->CenterTitle(1);
TH1D *h1p_had = new TH1D("h1p_had","Hadron Momentum Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; P_{h} (GeV/c)", int((EIC_Mom_Max_h+1)*4), 0.0, EIC_Mom_Max_h+1);
h1p_had->GetXaxis()->CenterTitle(1);
TH1D *h1the_ele = new TH1D("h1the_ele","Electron Theta Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; #theta_{e} (Degree)", int(EIC_Th_Max_e/2), -0., EIC_Th_Max_e);
h1the_ele->GetXaxis()->CenterTitle(1);
TH1D *h1the_had = new TH1D("h1the_had","Hadron Theta Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; #theta_{h} (Degree)", int(EIC_Th_Max_h/2), -0., EIC_Th_Max_h);
h1the_had->GetXaxis()->CenterTitle(1);

c4->cd(1);
T0->Draw("mom_ele>>h1p_ele",TCut(allcut),"");
c4->cd(2);
T0->Draw("mom_had>>h1p_had",TCut(allcut),"");
c4->cd(3);
T0->Draw("theta_ele*180/3.14>>h1the_ele",TCut(allcut),"");
c4->cd(4);
T0->Draw("theta_had*180/3.14>>h1the_had",TCut(allcut),"");
savestring.Clear();
savestring = "c12_pip_labvalues1d" + outputfileending;
c4->Print(savestring);

 //Finding of limits for generated electron and hadron momentum and angle, adding extra 5bins uncertainty on the values and compared with original limits in generation
 ele_pmax = h1p_ele->GetBinCenter(h1p_ele->FindLastBinAbove(0,1)) + 5*h1p_ele->GetBinWidth(1);
 ele_themax = h1the_ele->GetBinCenter(h1the_ele->FindLastBinAbove(0,1)) + 5*h1the_ele->GetBinWidth(1);
 ele_themin = h1the_ele->GetBinCenter(h1the_ele->FindFirstBinAbove(0,1)) - 5*h1the_ele->GetBinWidth(1);
 had_pmax = h1p_had->GetBinCenter(h1p_had->FindLastBinAbove(0,1)) + 5*h1p_had->GetBinWidth(1);
 had_themax = h1the_had->GetBinCenter(h1the_had->FindLastBinAbove(0,1)) + 5*h1the_had->GetBinWidth(1);
// cout << "before if Determined Limits: mom ele max " << ele_pmax << " , theta ele " << ele_themin << " - " << ele_themax << " , mom had max " << had_pmax << " , had theta max " << had_themax << endl;

//Check if determined optimized limits are outside of generated limits, if so set the corresponding values to the generated values
 if (ele_pmax > EIC_Mom_Max_e)  ele_pmax = EIC_Mom_Max_e;
 if (ele_themax > EIC_Th_Max_e) ele_themax = EIC_Th_Max_e;
 if (ele_themin < EIC_Th_Min_e) ele_themin = EIC_Th_Min_e;
 if (had_pmax > EIC_Mom_Max_h)  had_pmax = EIC_Mom_Max_h;
 if (had_themax > EIC_Th_Max_h) had_themax = EIC_Th_Max_h;

 Double_t generated_electron_phase_space =(cos(EIC_Th_Min_e*degtorad) - cos(EIC_Th_Max_e*degtorad))*(EIC_Ph_Max_e*degtorad - EIC_Ph_Min_e*degtorad)*(EIC_Mom_Max_e - EIC_Mom_Min_e);
 Double_t generated_hadron_phase_space   =(cos(EIC_Th_Min_h*degtorad) - cos(EIC_Th_Max_h*degtorad))*(EIC_Ph_Max_h*degtorad - EIC_Ph_Min_h*degtorad)*(EIC_Mom_Max_h - EIC_Mom_Min_h);
 Double_t generated_phase_space=generated_electron_phase_space*generated_hadron_phase_space;
 Double_t useful_electron_phase_space  = (cos(ele_themin*degtorad) - cos(ele_themax*degtorad))*(EIC_Ph_Max_e*degtorad - EIC_Ph_Min_e*degtorad)*(ele_pmax - EIC_Mom_Min_e);
 Double_t useful_hadron_phase_space =  (cos(EIC_Th_Min_h*degtorad) - cos(had_themax*degtorad))*(EIC_Ph_Max_h*degtorad - EIC_Ph_Min_h*degtorad)*(had_pmax - EIC_Mom_Min_h);
 Double_t useful_phase_space =  useful_electron_phase_space * useful_hadron_phase_space;
 cout << "generated PS: " << generated_phase_space << ", useful_phase_space: " << useful_phase_space << ". Ratio (gen/useful) "  << generated_phase_space/useful_phase_space << endl;
 cout << "Determined Limits: mom ele max " << ele_pmax << " , theta ele " << ele_themin << " - " << ele_themax << " , mom had max " << had_pmax << " , had theta max " << had_themax << endl;

 //2-dim histograms of generated lab values
  TCanvas *c3 = new TCanvas("c3","Generated Lab values 2d Histos", 1024, 800);
  c3->Divide(2,2);

  TH2D *hTP_ele = new TH2D("hTP_ele","Electron Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; #theta_{e} (Degree); P_{e} (GeV/c)", int(EIC_Th_Max_e/2), -0., EIC_Th_Max_e, int((EIC_Mom_Max_e+1)*4), 0.0, EIC_Mom_Max_e+1);
  hTP_ele->GetXaxis()->CenterTitle(1);
  hTP_ele->GetYaxis()->CenterTitle(1);
  TH2D *hTP_had = new TH2D("hTP_had","Hadron Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; #theta_{h} (Degree); P_{h} (GeV/c)", int(EIC_Th_Max_h/2), -0., EIC_Th_Max_h, int((EIC_Mom_Max_h+1)*4), 0.0, EIC_Mom_Max_h+1);
  hTP_had->GetXaxis()->CenterTitle(1);
  hTP_had->GetYaxis()->CenterTitle(1);

  c3->cd(1); //gPad->SetLogz(1);
  T0->Draw("mom_ele:theta_ele*180/3.14>>hTP_ele",TCut(allcut),"colz");
  TBox *elerange1 = new TBox(ele_themin,EIC_Mom_Min_e,ele_themax,ele_pmax);
  elerange1->SetLineColor(2);
  elerange1->SetLineStyle(2);
  elerange1->SetLineWidth(2);
  elerange1->SetFillColorAlpha(2,0.1);
  elerange1->Draw("SAME");
  c3->cd(2); //gPad->SetLogz(1);
  T0->Draw("mom_had:theta_had*180/3.14>>hTP_had",TCut(allcut),"colz");
  TBox *hadrange1 = new TBox(EIC_Th_Min_h,EIC_Mom_Min_h,had_themax,had_pmax);
  hadrange1->SetLineColor(2);
  hadrange1->SetLineStyle(2);
  hadrange1->SetLineWidth(2);
  hadrange1->SetFillColorAlpha(2,0.1);
  hadrange1->Draw("SAME");

  TH2D *hThePhi_ele = new TH2D("hThePhi_ele","Electron Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; #theta_{e} (Degree); #phi_{e} (Degree)", int(EIC_Th_Max_e/4), -0., EIC_Th_Max_e, 90, -180, 180);
  hThePhi_ele->GetXaxis()->CenterTitle(1);
  hThePhi_ele->GetYaxis()->CenterTitle(1);
  TH2D *hThePhi_had = new TH2D("hThePhi_had","Hadron Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; #theta_{h} (Degree); #phi_{h} (Degree)", int(EIC_Th_Max_h/4), -0., EIC_Th_Max_h, 90, -180, 180);
  hThePhi_had->GetXaxis()->CenterTitle(1);
  hThePhi_had->GetYaxis()->CenterTitle(1);

  c3->cd(3); //gPad->SetLogz(1);
  T0->Draw("phi_ele*180/3.14:theta_ele*180/3.14>>hThePhi_ele",TCut(allcut),"colz");
  c3->cd(4);// gPad->SetLogz(1);
  T0->Draw("phi_had*180/3.14:theta_had*180/3.14>>hThePhi_had",TCut(allcut),"colz");
  savestring.Clear();
  savestring = "c12_pip_labvalues2d" + outputfileending;
  c3->Print(savestring);



//Further procedure for pt<1 and Q2<10, xcut and weighting -> file choice 1 and cutchoice either 12 or 13
if (processfiles == 1 && (cutchoice == 12 || cutchoice == 13)) {

  Double_t nof_nozbinQ2bincuts= (double) T0->GetEntries(TCut(allcut)) ;
  cout << "Number of events with Q2<10, pt<1 and choosen xcut (no zbin and Q2bin cuts): " << nof_nozbinQ2bincuts << endl;
  //double z_min = 0.0, z_max = 1.2;
  //const int zbin = 7;
  //const double z_cut[8] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};

  double z_binmin = 0.2, z_binmax = 0.8;
  const int zbin = 6;
  const double z_cut[7] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  //double Q2log_min = 0.0, Q2log_max = 1.8;
  //const int Q2bin=8;
  //const double Q2_log[9] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};

  double Q2binmin = 2.0, Q2binmax = 10;
  const int Q2bin=8;
  const double Q2_cut[9] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

  TH2D *test = new TH2D("test","Determined electron pmax ; z (GeV) ; Q2 (GeV^{2}) ",zbin,z_cut[0],z_cut[6],Q2bin,Q2_cut[0],Q2_cut[8]);
  TH2D *test2 = new TH2D("test2","Determined electron themax ; z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[6],Q2bin,Q2_cut[0],Q2_cut[8]);
  TH2D *test3 = new TH2D("test3","Determined electron themin ; z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[6],Q2bin,Q2_cut[0],Q2_cut[8]);
  TH2D *test4 = new TH2D("test3","Determined hadron pmax ; z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[6],Q2bin,Q2_cut[0],Q2_cut[8]);
  TH2D *pspart = new TH2D("pspart","Relative Counts per z and Q2 bin in percent ; z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[6],Q2bin,Q2_cut[0],Q2_cut[8]);
  test->GetXaxis()->CenterTitle(1);
  test->GetYaxis()->CenterTitle(1);
  test2->GetXaxis()->CenterTitle(1);
  test2->GetYaxis()->CenterTitle(1);
  test3->GetXaxis()->CenterTitle(1);
  test3->GetYaxis()->CenterTitle(1);
  test4->GetXaxis()->CenterTitle(1);
  test4->GetYaxis()->CenterTitle(1);
  pspart->GetXaxis()->CenterTitle(1);
  pspart->GetYaxis()->CenterTitle(1);

  TCanvas *canv_Q2bins_ele = new TCanvas("canv_Q2bins_ele","Generated 2d Lab values (e) for each zbin and one Q2bin ", 1200, 800);
  canv_Q2bins_ele->Divide(3,2);
  TCanvas *canv_Q2bins_had = new TCanvas("canv_Q2bins_had","Generated 2d Lab values (h) for each zbin and one Q2bin ", 1200, 800);
  canv_Q2bins_had->Divide(3,2);
  for (Int_t j=0;j<Q2bin;j++){
  //for (Int_t j=3;j<4;j++){//Test for one Q2 bin
    //canv_Q2bins_ele->Clear();
  //  canv_Q2bins_had->Clear();
    for (Int_t i=0;i<zbin;i++){
      TString cut = allcut+"*"+Form("(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
      double counts_in_bin = double(T0->GetEntries(TCut(cut)));
      cout << "Nof events for zbin " << i << " and Q2 bin " << j << ": "  << counts_in_bin << endl;
      pspart->SetBinContent(i+1,j+1,counts_in_bin/nof_nozbinQ2bincuts*100);

      TString histotitle = Form("%.1f #leq Q2 < %.1f and %.1f #leq z < %.1f ; #theta_{e} (Degree); P_{e} (GeV/c)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
      TH2D *htemp1 = new TH2D("htemp1",histotitle, 25, -0.5, 99.5, 24, -0.5, 11.5);
      htemp1->GetXaxis()->CenterTitle(1);
      htemp1->GetYaxis()->CenterTitle(1);
      histotitle.Clear();
      histotitle = Form("%.1f #leq Q2 < %.1f and %.1f #leq z < %.1f ; #theta_{h} (Degree); P_{h} (GeV/c)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
      TH2D *htemp2 = new TH2D("htemp2",histotitle, 47, -0.5, 187.5, 22, -0.5, 11.5);
      htemp2->GetXaxis()->CenterTitle(1);
      htemp2->GetYaxis()->CenterTitle(1);

      //Drawing of Electron Generated values (momentum and theta), determination of max and min range and drawing of box
      canv_Q2bins_ele->cd(i+1); //gPad->SetLogz(1);
      T0->Draw("mom_ele:theta_ele*180/3.14>>htemp1",TCut(cut),"colz");
      TH1D *htemp1x =  htemp1->ProjectionX("htemp1x",0);
      TH1D *htemp1y =  htemp1->ProjectionY("htemp1y",0);
      ele_pmax = htemp1y->GetBinCenter(htemp1y->FindLastBinAbove(0,1)) + 4*htemp1y->GetBinWidth(1);
      ele_themax = htemp1x->GetBinCenter(htemp1x->FindLastBinAbove(0,1)) + 4*htemp1x->GetBinWidth(1);
      ele_themin = htemp1x->GetBinCenter(htemp1x->FindFirstBinAbove(0,1)) - 4*htemp1x->GetBinWidth(1);

      if (ele_pmax > htemp1y->GetXaxis()->GetXmax())  ele_pmax = htemp1y->GetXaxis()->GetXmax();
      if (ele_themax > EIC_Th_Max_e) ele_themax = EIC_Th_Max_e;
      if (ele_themin < EIC_Th_Min_e) ele_themin = EIC_Th_Min_e;
      TBox *elerange2 = new TBox(ele_themin,EIC_Mom_Min_e,ele_themax,ele_pmax);
      elerange2->SetLineColor(2);
      elerange2->SetLineStyle(2);
      elerange2->SetLineWidth(2);
      elerange2->SetFillColorAlpha(2,0.1);
      elerange2->Draw("SAME");
      test->SetBinContent(i+1,j+1,ele_pmax);
      test2->SetBinContent(i+1,j+1,ele_themax);
      test3->SetBinContent(i+1,j+1,ele_themin);

    //  canv_Q2bins_ele->Update();
        //Drawing of Hadron Generated values (momentum and theta), determination of max and min range and drawing of box
      canv_Q2bins_had->cd(i+1); //gPad->SetLogz(1);
      T0->Draw("mom_had:theta_had*180/3.14>>htemp2",TCut(cut),"colz");
      TH1D *htemp2x =  htemp2->ProjectionX("htemp2x",0);
      TH1D *htemp2y =  htemp2->ProjectionY("htemp2y",0);
      had_pmax = htemp2y->GetBinCenter(htemp2y->FindLastBinAbove(0,1)) + 4*htemp2y->GetBinWidth(1);
      had_themax = htemp2x->GetBinCenter(htemp2x->FindLastBinAbove(0,1)) + 4*htemp2x->GetBinWidth(1);
      if (had_pmax > EIC_Mom_Max_h)  had_pmax = EIC_Mom_Max_h;
      if (had_themax > EIC_Th_Max_h) had_themax = EIC_Th_Max_h;
      TBox *hadrange2 = new TBox(EIC_Th_Min_h,EIC_Mom_Min_h,had_themax,had_pmax);
      hadrange2->SetLineColor(2);
      hadrange2->SetLineStyle(2);
      hadrange2->SetLineWidth(2);
      hadrange2->SetFillColorAlpha(2,0.1);
      hadrange2->Draw("SAME");
      test4->SetBinContent(i+1,j+1,had_pmax);
      cout << "Determined Limits: mom ele max " << ele_pmax << " , theta ele " << ele_themin << " - " << ele_themax << " , mom had max " << had_pmax << " , had theta max " << had_themax << endl;
    //  canv_Q2bins_had->Update();


    }
    savestring.Clear();
    savestring = Form("lab2d_ele_cut%i_Q2bin%i.pdf",cutchoice,j);
    canv_Q2bins_ele->Print(savestring);
    savestring.Clear();
    savestring = Form("lab2d_had_cut%i_Q2bin%i.pdf",cutchoice,j);
    canv_Q2bins_had->Print(savestring);
  }

 TCanvas *ct = new TCanvas("ct","testhist");
 ct->Divide(2,1);
 ct->cd(1);
 test->Draw("COLZ");
 ct->cd(2);
 test2->Draw("COLZ");
 //ct->cd(3);
 //test3->Draw("COLZ");
 //ct->cd(4);
 //test4->Draw("COLZ");
 ct->SaveAs("test.pdf");


 TCanvas *ct2 = new TCanvas("ct2","pspart");
 pspart->Draw("COLZ");
 savestring.Clear();
 savestring = "c12_pspart" + outputfileending;
 ct2->SaveAs(savestring);

}
delete T0;

/**
  TCanvas *c5 = new TCanvas("c5","c5", 1024, 800);
  c5->Divide(2,2);

  c5->cd(1); gPad->SetLogz(1);
  T0->Draw("mom_ele:theta_ele*180/3.14>>htemp1(110, -5.0, 105., 100, 0.0, 12)",TCut(cut),"colz");
  c5->cd(2); gPad->SetLogz(1);
  T0->Draw("mom_had:theta_had*180/3.14>>htemp2(80, -5.0, 185., 100, 0.0, 12.)",TCut(cut),"colz");
  c5->cd(3); gPad->SetLogz(1);
  T0->Draw("phi_ele*180/3.14:theta_ele*180/3.14>>htemp3(110, -5.0, 105., 90, -180, 180)",TCut(cut),"colz");
  c5->cd(4); gPad->SetLogz(1);
  T0->Draw("phi_had*180/3.14:theta_had*180/3.14>>htemp4(80, -5.0, 185., 90, -180, 180)",TCut(cut),"colz");
  c5->Update();
  TString savestring = "c12_pip_labvalues_Q2-";
  savestring += Form("%i_z-%i",j+1,i+1);
  savestring += ".pdf";
//  cout << savestring << endl;
  c5->Print(savestring);
  TCanvas *c6 = new TCanvas("c6","c6", 1024, 800);
  c6->Divide(2,2);

  c6->cd(1); gPad->SetLogz(1);
  T0->Draw("z:pt>>htemp11(100, 0.0, 5.0,100, 0.15, 0.85)",TCut(cut),"colz");
  c6->cd(4); gPad->SetLogz(1);
  T0->Draw("log10(Q2):log10(x)>>htemp21( 100,-3.5, 0.0,100, -0.1, 3.0)",TCut(cut),"colz");
  c6->cd(3); gPad->SetLogz(1);
  T0->Draw("log10(Q2):pt>>htemp31(100, 0.0, 5.0,100, -0.1, 3.0)",TCut(cut),"colz");
  c6->cd(2); gPad->SetLogz(1);
  T0->Draw("z:log10(x)>>htemp41(100,-3.5, 0.0,100, 0.15, 0.85)",TCut(cut),"colz");
  c6->Update();
  savestring.Clear();
  savestring = "c12_pip_kinvalues2d_Q2-";
  savestring += Form("%i_z-%i",j+1,i+1);
  savestring += ".pdf";
  c6->Print(savestring);
    }
  }**/
}
