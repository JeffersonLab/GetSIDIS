{
    TString target="proton";
    TString hadron="pion";
    int Ebeam=5;
    int Ibeam=25;

    TChain* T0 = new TChain("T");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_1_0.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_1_1.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_1_2.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_2_0.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_2_1.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_2_2.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_3_0.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_3_1.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_3_2.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_4_0.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_4_1.root");
    T0->Add("./rootfiles/EIC_A1_pion_5_25_4_2.root");

    gStyle->SetOptStat(0);
    //1e34cm^-1*s^-1 for proton, 1e-33 is for nbar->cm
    double luminosity = 1e34 * 1e-33;// 
    double time = 100 * 24 * 60 * 60; //100 days, into seconds
    time /=3; //divided by 3 means that I add three sets of rootfiles together, equilibrium to run three-times longer

    //TCanvas *c1 = new TCanvas("c1","c1", 1000,600);
    //c1->Divide(2,1);
    //c1->cd(1); gPad->SetLogz(1);
    //T0->Draw("mom_ele:theta_ele*180/3.14>>h1(200, 10.,180, 200, 0.0, 20.)",Form("%f*%f*weight_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.2 && z<0.8)", luminosity, time),"colz");
    //h1->SetXTitle("#theta_{ele} (Degrees)");
    //h1->SetYTitle("P_{ele} (GeV/c)");
    //h1->SetTitle("");
    //h1->GetXaxis()->CenterTitle(1);
    //h1->GetYaxis()->CenterTitle(1);
    //c1->cd(2); gPad->SetLogz(1);
    //T0->Draw("mom_had:theta_had*180/3.14>>h2(200, 0.,180, 200, 0.0, 20.0)",Form("%f*%f*weight_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.2 && z<0.8)", luminosity, time),"colz");
    //h2->SetXTitle("#theta_{had} (Degrees)");
    //h2->SetYTitle("P_{had} (GeV/c)");
    //h2->SetTitle("");
    //h2->GetXaxis()->CenterTitle(1);
    //h2->GetYaxis()->CenterTitle(1);
    //c1.Print(Form("%s_acc_%s_%dGeV_%dGeV.pdf",target.Data(), hadron.Data(), Ebeam, Ibeam));
    //c1.Print(Form("%s_acc_%s_%dGeV_%dGeV.png",target.Data(), hadron.Data(), Ebeam, Ibeam));
    
    TCanvas *c2 = new TCanvas("c2","c2", 800,600);
    T0->Draw("Q2:x>>k6(200,0.,1.0,200, 0.5.,400.)",Form("%f*%f*weight_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.2 && z<0.8)", luminosity, time),"colz");
    k6->SetXTitle("x_{B}");
    k6->SetYTitle("Q^{2} (GeV^{2})");
    k6->SetTitle(Form("p(e,e'#pi^{+})X: P_{e}=%dGeV, P_{I}=%dGeV", Ebeam, Ibeam));
    k6->GetXaxis()->CenterTitle(1);
    k6->GetYaxis()->CenterTitle(1);
    k6->GetXaxis()->SetTitleSize(0.06);
    k6->GetYaxis()->SetTitleSize(0.06);
    k6->GetXaxis()->SetTitleOffset(0.8);
    k6->GetYaxis()->SetTitleOffset(0.8);
    gPad->SetLogz(1);
    c2.Print(Form("%s_xQ2_%s_%dGeV_%dGeV.pdf",target.Data(), hadron.Data(), Ebeam, Ibeam));
    c2.Print(Form("%s_xQ2_%s_%dGeV_%dGeV.png",target.Data(), hadron.Data(), Ebeam, Ibeam));
 

    TCanvas *c3 = new TCanvas("c3","c3", 800,600);
    T0->Draw("log10(Q2):log10(x)>>k5(200,-3.0.,0.0,200,-0.2,2.8)",Form("%f*%f*weight_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.2 && z<0.8)", luminosity, time),"colz");
    k5->SetXTitle("log10(x_{B})");
    k5->SetYTitle("log10(Q^{2})");
    k5->SetTitle(Form("p(e,e'#pi^{+})X: P_{e}=%dGeV, P_{I}=%dGeV", Ebeam, Ibeam));
    k5->GetXaxis()->CenterTitle(1);
    k5->GetYaxis()->CenterTitle(1);
    k5->GetXaxis()->SetTitleSize(0.06);
    k5->GetYaxis()->SetTitleSize(0.06);
    k5->GetXaxis()->SetTitleOffset(1.0);
    k5->GetYaxis()->SetTitleOffset(0.7);
    gPad->SetLogz(1);
    c3.Print(Form("%s_xQ2log_%s_%dGeV_%dGeV.pdf",target.Data(), hadron.Data(), Ebeam, Ibeam));
    c3.Print(Form("%s_xQ2log_%s_%dGeV_%dGeV.png",target.Data(), hadron.Data(), Ebeam, Ibeam));
 

    TCanvas *c4 = new TCanvas("c4","c4", 800,600);
    T0->Draw("pt:z>>k4(200,0.1,1.0,200, 0.0.,2.8)",Form("%f*%f*weight_hp*(isphy_hp>0&&Q2>1 && W>2.3 && Wp>1.6 && z>0.0 && z<1.2)", luminosity, time),"colz");
    k4->SetXTitle("z");
    k4->SetYTitle("p_{T} (GeV/c)");
    k4->SetTitle(Form("p(e,e'#pi^{+})X: P_{e}=%dGeV, P_{I}=%dGeV", Ebeam, Ibeam));
    k4->GetXaxis()->CenterTitle(1);
    k4->GetYaxis()->CenterTitle(1);
    k4->GetXaxis()->SetTitleSize(0.06);
    k4->GetYaxis()->SetTitleSize(0.06);
    k4->GetXaxis()->SetTitleOffset(0.8);
    k4->GetYaxis()->SetTitleOffset(0.7);
    gPad->SetLogz(1);
    c4.Print(Form("%s_zpt_%s_%dGeV_%dGeV.pdf",target.Data(), hadron.Data(), Ebeam, Ibeam));
    c4.Print(Form("%s_zpt_%s_%dGeV_%dGeV.png",target.Data(), hadron.Data(), Ebeam, Ibeam));
 

}
