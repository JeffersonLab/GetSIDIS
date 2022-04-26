#include "plot_collins1.C"
void plot_all(){
	gStyle->SetLabelSize(0.1,"Y");
	gStyle->SetTitleSize(0.1,"Y");
	gStyle->SetTitleOffset(0.8,"Y");
	gStyle->SetLabelSize(0.1,"X");
	gStyle->SetTitleSize(0.1,"X"); 

	gStyle->SetPadBorderMode(0);


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
				gPad->SetLeftMargin(0.2);
			}
			if (j!=7){
				gPad->SetRightMargin(0);
			}else{
				gPad->SetRightMargin(0.2);
			}

			gPad->SetFillColor(10);
			if (j!=7){
				plot(1,i+1,j+1,1);
			}else{
				plot(1,i+1,j+1,1,1);
			}
			if (i==0){
				plot_feng(1,1,j+1);
				//plot_feng(1,2,j+1);

			}
			if (i==0&&j==0){
				plot_e06010(1,3);
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
				gPad->SetLeftMargin(0.2);
			}
			if (j!=7){
				gPad->SetRightMargin(0);
			}else{
				gPad->SetRightMargin(0.2);
			}


			if (j==7){
				plot(2,i+1,j+1,1,1);
			}else{
				plot(2,i+1,j+1,1);
			}
			if (i==0){
				plot_feng(2,1,j+1);
				//plot_feng(2,2,j+1);

			}
			if (i==0&&j==0){
				plot_e06010(2,3);
				plot_ase(2,1);
				//plot_ase(2,2);
				plot_pretz(2,1);
				plot_pretz(2,2);

			}

		}
	}

	   c1->Print("project_all_pip_collins.pdf");
	   c2->Print("project_all_pim_collins.pdf");

}
