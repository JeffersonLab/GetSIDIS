{
    gStyle->SetOptStat(0);
    TChain *T0 = new TChain("T");
    for(int i=0;i<100;i++){
        for(int j=1;j<=4;j++)
            //T0->Add(Form("./c12_pion/EIC_A12_pion_10_600_%d_%d.root", j, i));
            T0->Add(Form("./prot_pion/EIC_A1_pion_10_100_%d_%d.root", j, i));
    }

    TCanvas *c1 = new TCanvas("c1","c1", 800, 600);
    TH2D *h1 = new TH2D("h1","^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; p_{T} (GeV); z (GeV))", 100, 0.0, 5.0,100, 0.15, 0.85);
    h1->GetXaxis()->CenterTitle(1);
    h1->GetYaxis()->CenterTitle(1);
    gPad->SetLogz(1);
    T0->Draw("z:pt>>h1","weight_hm*(weight_hm>0.0)","colz");
    c1->Print("c12_pim_z_pt_A600.pdf");
    
    TCanvas *c2 = new TCanvas("c2","c2", 800, 600);
    gPad->SetLogz(1);
    TH2D *h2 = new TH2D("h2","^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; log10(x); log10(Q^{2}) (GeV^{2}))", 100,-3.5, 0.0,100, -0.1, 3.0);
    h2->GetXaxis()->CenterTitle(1);
    h2->GetYaxis()->CenterTitle(1);
    T0->Draw("log10(Q2):log10(x)>>h2","weight_hp","colz");
    c2->Print("c12_pip_Q2_x_log_A600.pdf");

    TCanvas *c3 = new TCanvas("c3","c3", 600, 800);
    c3->Divide(1,2);
    
    TH2D *hTP_ele = new TH2D("hTP_ele","Electron Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; #theta_{e} (Degree); P_{e} (GeV/c)", 150, -5.0, 145., 100, 0.0, 14);
    hTP_ele->GetXaxis()->CenterTitle(1);
    hTP_ele->GetYaxis()->CenterTitle(1);
    TH2D *hTP_had = new TH2D("hTP_had","Hadron Acceptance, ^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV; #theta_{h} (Degree); P_{h} (GeV/c)", 190, -5.0, 185., 100, 0.0, 12.);
    hTP_had->GetXaxis()->CenterTitle(1);
    hTP_had->GetYaxis()->CenterTitle(1);

    c3->cd(1); gPad->SetLogz(1);
    T0->Draw("mom_ele:theta_ele*180/3.14>>hTP_ele","weight_hp","colz");
    c3->cd(2); gPad->SetLogz(1);
    T0->Draw("mom_had:theta_had*180/3.14>>hTP_had","weight_hp","colz");
    c3->Print("c12_pip_acc_A600.pdf");

}

