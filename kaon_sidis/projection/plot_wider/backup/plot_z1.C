#include "plot.C"
void plot_z1(){
  TCanvas **c1 = new TCanvas*[8];
  for (Int_t j=0;j!=8;j++){
    c1[j] = new TCanvas(Form("c1_%d",j),Form("c1_%d",j),1200,900);
    c1[j]->SetFillColor(10);
    c1[j]->Divide(3,2);
    for (Int_t i=0;i!=6;i++){
      c1[j]->cd(i+1);
      plot(1,i+1,j+1);
    }
  }
  
  TCanvas **c2 = new TCanvas*[8];
  for (Int_t j=0;j!=8;j++){
    c2[j] = new TCanvas(Form("c2_%d",j),Form("c2_%d",j),1200,900);
    c2[j]->SetFillColor(10);
    c2[j]->Divide(3,2);
    for (Int_t i=0;i!=6;i++){
      c2[j]->cd(i+1);
      plot(2,i+1,j+1);
    }
  }
}
