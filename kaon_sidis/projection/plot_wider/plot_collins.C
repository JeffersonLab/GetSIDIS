#include "plot_collins1.h"
void plot_all(){
	gStyle->SetLabelSize(0.1,"Y");
	gStyle->SetTitleSize(0.1,"Y");
	gStyle->SetTitleOffset(0.8,"Y");
	gStyle->SetLabelSize(0.1,"X");
	gStyle->SetTitleSize(0.1,"X"); 

	gStyle->SetPadBorderMode(0);

    int target_flag=0; cout<<"-- Target (1->NH3, 3->He3)"; cin >> target_flag;
    int particle_flag=0; cout<<"-- Particle (1->pion, 2->kaon)"; cin >> particle_flag;
    int Ebeam=0; cout<<"-- EBeam (11, 8 (8.8GeV))"; cin >> Ebeam;

    TString target = "X";
	if(target_flag==1)
		target ="NH3";
	else if(target_flag==2)
		target ="D2";
	else if(target_flag==3)
		target ="He3";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}
	TString particle = "X";
	if(particle_flag==1)
		particle ="pion";
	else if(particle_flag==2)
		particle ="kaon";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}

	TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
	c1->SetFillColor(10);
	c1->Divide(8,6,0,0);
	for (Int_t i=0;i!=6;i++){
		for (Int_t j=0;j!=8;j++){
			c1->cd(i*8+j+1);
			if (i!=5){
				gPad->SetBottomMargin(0);
			}else{
				gPad->SetBottomMargin(0.2);
			}
			gPad->SetTopMargin(0);



			if (j!=0){
				gPad->SetLeftMargin(0);
			}else{
				gPad->SetLeftMargin(0.3);
			}
			if (j!=7){
				gPad->SetRightMargin(0);
			}else{
				gPad->SetRightMargin(0.3);
			}

			gPad->SetFillColor(10);
			if (j!=7){
				plot(Ebeam, target_flag, particle_flag,i+1,j+1,1,2);
			}else{
				plot(Ebeam, target_flag, particle_flag,i+1,j+1,1,1);
			}
			if (i==0){
				plot_feng(1,1,j+1);
				//plot_feng(1,2,j+1);

			}
			if (i==0&&j==0){
				plot_e06010(1,1,3);
				plot_ase(1,1);
				//plot_ase(1,2);
				plot_pretz(1,1);
				plot_pretz(1,2);
			}
		}
	}

	TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
	c2->SetFillColor(10);
	c2->Divide(8,6,0,0);
	for (Int_t i=0;i!=6;i++){
		for (Int_t j=0;j!=8;j++){
			c2->cd(i*8+j+1);
			gPad->SetFillColor(10);
			gPad->SetTopMargin(0);
			if (i!=5){
				gPad->SetBottomMargin(0);
			}else{
				gPad->SetBottomMargin(0.2);
			}

			if (j!=0){
				gPad->SetLeftMargin(0);
			}else{
				gPad->SetLeftMargin(0.30);
			}
			if (j!=7){
				gPad->SetRightMargin(0);
			}else{
				gPad->SetRightMargin(0.3);
			}


			if (j!=7){
				plot(Ebeam, target_flag, -1*particle_flag,i+1,j+1,1,2);
			}else{
				plot(Ebeam, target_flag, -1*particle_flag,i+1,j+1,1,1);
			}
			if (i==0){
				plot_feng(2,1,j+1);
				//plot_feng(2,2,j+1);

			}
			if (i==0&&j==0){
				plot_e06010(2,1,3);
				plot_ase(2,1);
				//plot_ase(2,2);
				plot_pretz(2,1);
				plot_pretz(2,2);

			}

		}
	}

	c1->Print(Form("project_%s_%s_hp_E%d_collins.pdf", target.Data(), particle.Data(), Ebeam));
	c1->Print(Form("project_%s_%s_hp_E%d_collins.png", target.Data(), particle.Data(), Ebeam));
	c2->Print(Form("project_%s_%s_hm_E%d_collins.pdf", target.Data(), particle.Data(), Ebeam));
	c2->Print(Form("project_%s_%s_hm_E%d_collins.png", target.Data(), particle.Data(), Ebeam));

}
