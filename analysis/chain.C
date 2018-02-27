{
    TString target ="";  cout<<"--- Target (He3 or NH3) = "; cin>>target;
    TString hadron=""; cout<<"--- Hadron (pion or kaon) = "; cin>>hadron;
    Int_t Ebeam=""; cout<<"--- Beam Energy (11 or 8)= "; cin>>Ebeam;

    TChain* T0 = new TChain("T");
    T0->Add(Form("skim_%s_%s_E%d_z1_Qsq1.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z1_Qsq2.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z1_Qsq3.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z1_Qsq4.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z1_Qsq5.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z1_Qsq6.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z2_Qsq1.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z2_Qsq2.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z2_Qsq3.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z2_Qsq4.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z2_Qsq5.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z2_Qsq6.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z3_Qsq1.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z3_Qsq2.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z3_Qsq3.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z3_Qsq4.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z3_Qsq5.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z3_Qsq6.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z4_Qsq1.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z4_Qsq2.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z4_Qsq3.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z4_Qsq4.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z4_Qsq5.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z4_Qsq6.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z5_Qsq1.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z5_Qsq2.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z5_Qsq3.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z5_Qsq4.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z5_Qsq5.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z5_Qsq6.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z6_Qsq1.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z6_Qsq2.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z6_Qsq3.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z6_Qsq4.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z6_Qsq5.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z6_Qsq6.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z7_Qsq1.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z7_Qsq2.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z7_Qsq3.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z7_Qsq4.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z7_Qsq5.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z7_Qsq6.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z8_Qsq1.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z8_Qsq2.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z8_Qsq3.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z8_Qsq4.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z8_Qsq5.root",target.Data(), hadron.Data(), Ebeam));
    T0->Add(Form("skim_%s_%s_E%d_z8_Qsq6.root",target.Data(), hadron.Data(), Ebeam));

    TCanvas *c1 = new TCanvas("c1","c1", 800,1000);
    c1->Divide(2,3);
    c1->cd(1); gPad->SetLogz(1);
    T0->Draw("mom_ele:theta_ele>>h1(200, 0.,30, 200, 0.5, 9.0)","weight_hp*(acc_f_ele+acc_l_ele)*acc_f_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7)","colz");
    h1->SetXTitle("#theta_{ele}");
    h1->SetYTitle("P_{ele}");
    h1->SetTitle("");
    h1->GetXaxis()->CenterTitle(1);
    h1->GetYaxis()->CenterTitle(1);
    c1->cd(2); gPad->SetLogz(1);
    T0->Draw("mom_had:theta_had>>h2(200, 0.,30, 200, 0.5, 6.5)","weight_hp*(acc_f_ele+acc_l_ele)*acc_f_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7)","colz");
    h2->SetXTitle("#theta_{had}");
    h2->SetYTitle("P_{had}");
    h2->SetTitle("");
    h2->GetXaxis()->CenterTitle(1);
    h2->GetYaxis()->CenterTitle(1);
    c1->cd(3); gPad->SetLogz(1);
    T0->Draw("mom_ele:phi_ele>>h3(200, -181.,181, 200, 0.5, 9.0)","weight_hp*(acc_f_ele+acc_l_ele)*acc_f_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7)","colz");
    h3->SetXTitle("#phi_{ele}");
    h3->SetYTitle("P_{ele}");
    h3->SetTitle("");
    h3->GetXaxis()->CenterTitle(1);
    h3->GetYaxis()->CenterTitle(1);
    c1->cd(4); gPad->SetLogz(1);
    T0->Draw("mom_had:phi_had>>h4(200, -181.,181, 200, 0.5, 6.5)","weight_hp*(acc_f_ele+acc_l_ele)*acc_f_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7)","colz");
    h4->SetXTitle("#phi_{had}");
    h4->SetYTitle("P_{had}");
    h4->SetTitle("");
    h4->GetXaxis()->CenterTitle(1);
    h4->GetYaxis()->CenterTitle(1);
    c1->cd(5); gPad->SetLogz(1);
    T0->Draw("phi_ele:theta_ele>>h5(200,0.,30.,200, -181.,181)","weight_hp*(acc_f_ele+acc_l_ele)*acc_f_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7)","colz");
    h5->SetXTitle("#theta_{ele}");
    h5->SetYTitle("#phi_{ele}");
    h5->SetTitle("");
    h5->GetXaxis()->CenterTitle(1);
    h5->GetYaxis()->CenterTitle(1);
    c1->cd(6); gPad->SetLogz(1);
    T0->Draw("phi_had:theta_had>>h6(200,0.,30.,200, -181.,181)","weight_hp*(acc_f_ele+acc_l_ele)*acc_f_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7)","colz");
    h6->SetXTitle("#theta_{had}");
    h6->SetYTitle("#phi_{had}");
    h6->SetTitle("");
    h6->GetXaxis()->CenterTitle(1);
    h6->GetYaxis()->CenterTitle(1);

    c1.Print(Form("%s_acc_%s_%d.pdf",target.Data(), hadron.Data(), Ebeam));
    c1.Print(Form("%s_acc_%s_%d.png",target.Data(), hadron.Data(), Ebeam));

    TCanvas *c2 = new TCanvas("c2","c2", 800,1000);
    c2->Divide(2,3);
    c2->cd(1);
    T0->Draw("weight_hp*(acc_f_ele+acc_l_ele)*acc_f_hp>>kp1","mom_had<7.5&&isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7");
    c2->cd(2);
    T0->Draw("weight_hm*(acc_f_ele+acc_l_ele)*acc_f_hm>>km1","mom_had<7.5&&isphy_hm>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7");

    c2->cd(3);
    T0->Draw("weight_hp*(acc_f_ele+acc_l_ele)*acc_f_hp>>kp2","mom_had<5.0&&isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7");
    c2->cd(4);
    T0->Draw("weight_hm*(acc_f_ele+acc_l_ele)*acc_f_hm>>km2","mom_had<5.0&&isphy_hm>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7");

    c2->cd(5);
    T0->Draw("weight_hp*(acc_f_ele+acc_l_ele)*acc_f_hp>>kp3","mom_had<2.5&&isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7");
    c2->cd(6);
    T0->Draw("weight_hm*(acc_f_ele+acc_l_ele)*acc_f_hm>>km3","mom_had<2.5&&isphy_hm>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.3 && z<0.7");

    ofstream outf(Form("%s_rate_%s_%d.dat",target.Data(), hadron.Data(), Ebeam));

    cout<<"--- had+ (Ph<7.5GeV/c): "<< kp1->GetSum() * kp1->GetMean()<<endl;
    cout<<"--- had+ (Ph<5.0GeV/c): "<< kp2->GetSum() * kp2->GetMean()<<endl;
    cout<<"--- had+ (Ph<2.5GeV/c): "<< kp3->GetSum() * kp3->GetMean()<<endl;

    cout<<"--- had- (Ph<7.5GeV/c): "<< km1->GetSum() * km1->GetMean()<<endl;
    cout<<"--- had- (Ph<5.0GeV/c): "<< km2->GetSum() * km2->GetMean()<<endl;
    cout<<"--- had- (Ph<2.5GeV/c): "<< km3->GetSum() * km3->GetMean()<<endl;

    outf<<"--- had+ (Ph<7.5GeV/c): "<< kp1->GetSum() * kp1->GetMean()<<endl;
    outf<<"--- had+ (Ph<5.0GeV/c): "<< kp2->GetSum() * kp2->GetMean()<<endl;
    outf<<"--- had+ (Ph<2.5GeV/c): "<< kp3->GetSum() * kp3->GetMean()<<endl;

    outf<<"--- had- (Ph<7.5GeV/c): "<< km1->GetSum() * km1->GetMean()<<endl;
    outf<<"--- had- (Ph<5.0GeV/c): "<< km2->GetSum() * km2->GetMean()<<endl;
    outf<<"--- had- (Ph<2.5GeV/c): "<< km3->GetSum() * km3->GetMean()<<endl;
    outf.close();




}
