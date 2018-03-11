//////////////////////////////////////////////////////
// SIDIS Events Generators for SoLID, CLAS12 or EIC //
//                                                  //
//Note: Basically the same as Xin Qian's "collider" //
//      but the model is coded in "SIDIS.h"         //
//  -- Zhihong Ye, 06/10/2014                       //
//////////////////////////////////////////////////////
#include "GetSIDIS.h"
#include "SIDIS_new.h"
//#include "SIDIS_Lite.h" //this version doesn't include LHAPDF
//#include "SIDIS_Lite_LO.h" //this version doesn't include LHAPDF, contributions from s, sbar and g, and only LO PDF

int main(Int_t argc, char *argv[]){
    cout<<endl;
    cout<<"oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo"<<endl;
    cout<<"oO0 SIDIS Events Generators for SoLID or CLAS12 or EIC or Spectrometer  0Oo//"<<endl;
    cout<<"oO0  with nPDF (EPS09) and free-PDF (LHDAPDF) implemented.    0Oo//"<<endl;
    cout<<"oO0  -- Zhihong Ye, updated in 08/12/2016                     0Oo//"<<endl;
    cout<<"oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo"<<endl;
    cout<<endl;

    gRandom->SetSeed(0);// uses time for seed

    /*Inputs&Output{{{*/
    //initialize
    TString inputfilename = argv[1];
    Init(inputfilename.Data());

    if(FileNo==0){
        if(argc == 3) FileNo = atoi(argv[2]);
        else{
            cout<<"&& Oops, Tell me FileNo (0, 1, ...) = "; cin >> FileNo;
        }
    }

    Double_t Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon,rapidity,physical, jacoF;

    if (config != "EIC" && config != "SoLID" && config != "CLAS12" && config != "SPECT"){
        cout << "not supported config = "<<config.Data() << endl;
        return -1;
    }    

    //define output file
    Int_t target_flag = A;
    TString prefix=config.Data();
    prefix += Form("_A%d", target_flag);

    Double_t charge_pos = 0, charge_neg = 0;
    int pid_pos = 0,pid_neg=0;
    Double_t mass_had = 0.0;
    if (particle_flag == 1){
        prefix += "_pion";
        charge_pos = 1; pid_pos = 211; mass_had = 139.57/1000.;//GeV
        charge_neg =-1; pid_neg =-211; 
    }else if (particle_flag == 2){
        prefix += "_kaon";
        charge_pos = 1; pid_pos = 321; mass_had = 493.68/1000.;//GeV
        charge_neg =-1; pid_neg =-321; 
    }else{
        cout << "particle_flag is wrong +-1 and +-2" << endl;
        return -1;
    }

    //Ion mass is the target mass in SoLID or fixed target experiment
    Double_t ion_mass = GetIonMass(A, Z);//GeV
    /*}}}*/

    /*Define{{{*/
    Double_t theta_gen= 0.0 , phi_gen = 0.0, mom_gen = 0.0;

    Int_t count[4] = {0,0,0,0};

    Double_t mom_gen_ele,mom_gen_had;
    Double_t theta_gen_ele,theta_gen_had;
    Double_t phi_gen_ele,phi_gen_had;
    Double_t theta_q, theta_s,phi_h,phi_s,mom_ele,mom_had,theta_ele, theta_had,phi_ele,phi_had;
    Double_t dxs_incl,dxs_hm,dxs_hp,dxs_hm_sidis,dxs_hp_sidis,dilute_hp,dilute_hm;
    Double_t px_ele, py_ele,pz_ele, px_had, py_had, pz_had, E_ele,E_had;
    Double_t u_pdf, d_pdf, s_pdf, g_pdf, ubar_pdf, dbar_pdf, sbar_pdf;
    Double_t D_fav, D_unfav, D_s, D_g;
    ULong64_t nsim = 0, Nsim1 = 0, Nsim2 = 0, Nsim3 = 0,Nsim4 = 0;
    //For Beam Position and Vertex info
    Double_t vx_ele, vy_ele, vz_ele, vx_had, vy_had, vz_had;
    double beamsize_x_ele=0.0, beamsize_y_ele=0.0;
    double vertex_length =0.0, vertex_center=0.0;
    int isphy_hp, isphy_hm;

    //The idea is to generate a phase-space which is slightly larger than the actual one
    Double_t Mom_Max_e = 0.0, Mom_Min_e = 0.0, Mom_Max_h = 0.0,Mom_Min_h = 0.0;
    Double_t Th_Max_e = 0.0,Th_Min_e = 0.0,Th_Max_h = 0.0, Th_Min_h = 0.0;
    Double_t Ph_Max_e = 0.0,Ph_Min_e = 0.0,Ph_Max_h = 0.0, Ph_Min_h = 0.0;
    if(config=="SoLID" ){/*{{{*/
        Mom_Min_e = SoLID_Mom_Min_e;  Mom_Max_e = momentum_ele; 
        Mom_Min_h = SoLID_Mom_Min_h;  Mom_Max_h = SoLID_Mom_Max_h;
        Th_Min_e = SoLID_Th_Min_e; Th_Max_e = SoLID_Th_Max_e; 
        Th_Min_h = SoLID_Th_Min_h; Th_Max_h = SoLID_Th_Max_h;
        Ph_Min_e = SoLID_Ph_Min_e; Ph_Max_e = SoLID_Ph_Max_e; 
        Ph_Min_h = SoLID_Ph_Min_h; Ph_Max_h = SoLID_Ph_Max_h;

        beamsize_x_ele = SoLID_BeamSizeX_ele;
        beamsize_y_ele = SoLID_BeamSizeY_ele;
        vertex_length = SoLID_Target_Length;
        vertex_center = SoLID_Target_Center;

    }/*}}}*/
    
    if(config=="CLAS12" ){/*{{{*/
        Mom_Min_e = CLAS12_Mom_Min_e;  Mom_Max_e = momentum_ele; 
        Mom_Min_h = CLAS12_Mom_Min_h;  Mom_Max_h = CLAS12_Mom_Max_h;
        Th_Min_e = CLAS12_Th_Min_e; Th_Max_e = CLAS12_Th_Max_e; 
        Th_Min_h = CLAS12_Th_Min_h; Th_Max_h = CLAS12_Th_Max_h;
        Ph_Min_e = CLAS12_Ph_Min_e; Ph_Max_e = CLAS12_Ph_Max_e; 
        Ph_Min_h = CLAS12_Ph_Min_h; Ph_Max_h = CLAS12_Ph_Max_h;

        beamsize_x_ele = CLAS12_BeamSizeX_ele;
        beamsize_y_ele = CLAS12_BeamSizeY_ele;
        vertex_length = CLAS12_Target_Length;
        vertex_center = CLAS12_Target_Center;

    }/*}}}*/
   
    //A rough guess but people claim EIC to be a full-acceptance device!
    else if(config=="EIC" ){/*{{{*/
        Mom_Min_e = EIC_Mom_Min_e;  Mom_Max_e =  abs(momentum_ele - momentum_ion);//allow the max mom of the ele to be the total momentum 
        Mom_Min_h = EIC_Mom_Min_h;  Mom_Max_h = EIC_Mom_Max_h;
        Th_Min_e = EIC_Th_Min_e; Th_Max_e = EIC_Th_Max_e; 
        Th_Min_h = EIC_Th_Min_h; Th_Max_h = EIC_Th_Max_h;    
        Ph_Min_e = EIC_Ph_Min_e; Ph_Max_e = EIC_Ph_Max_e; 
        Ph_Min_h = EIC_Ph_Min_h; Ph_Max_h = EIC_Ph_Max_h;    

        beamsize_x_ele = EIC_BeamSizeX_ele;
        beamsize_y_ele = EIC_BeamSizeY_ele;
        vertex_length = EIC_Vertex_Length;
        vertex_center = EIC_Vertex_Center;
    }/*}}}*/
   
    else if(config=="SPECT"){/*{{{*/
        Mom_Min_e = SPECT_Mom_Min_e;  Mom_Max_e = SPECT_Mom_Max_e; 
        Mom_Min_h = SPECT_Mom_Min_h;  Mom_Max_h = SPECT_Mom_Max_h;
        Th_Min_e = SPECT_Th_Min_e; Th_Max_e = SPECT_Th_Max_e; 
        Th_Min_h = SPECT_Th_Min_h; Th_Max_h = SPECT_Th_Max_h;    
        Ph_Min_e = SPECT_Ph_Min_e; Ph_Max_e = SPECT_Ph_Max_e; 
        Ph_Min_h = SPECT_Ph_Min_h; Ph_Max_h = SPECT_Ph_Max_h;    

        beamsize_x_ele = SPECT_BeamSizeX_ele;
        beamsize_y_ele = SPECT_BeamSizeY_ele;

        vertex_length = SPECT_Target_Length;
        vertex_center = SPECT_Target_Center;
    }/*}}}*/

    Double_t electron_phase_space =(cos(Th_Min_e/DEG) - cos(Th_Max_e/DEG))*(Ph_Max_e/DEG - Ph_Min_e/DEG)*(Mom_Max_e - Mom_Min_e);
    Double_t hadron_phase_space   =(cos(Th_Min_h/DEG) - cos(Th_Max_h/DEG))*(Ph_Max_h/DEG - Ph_Min_h/DEG)*(Mom_Max_h - Mom_Min_h);
    Double_t Phase_space=electron_phase_space*hadron_phase_space;           //electron*hadron phase space eg, for electron: delta_cos_theta*delta_phi*delta_energy
    cout<<" -- For Config="<<config<<" Phase_space: "<<electron_phase_space<<"	"<<hadron_phase_space<<"	"<<Phase_space<<endl;

    Double_t Q2_CutOff = 10.0; //A default setting
    if(config=="CLAS12")
        Q2_CutOff = 12.0; //A loose Q2 cut-off just for clas12 to see how far it can go

    /*}}}*/

    /*New ROOT files, Trees and  Branches{{{*/
    //create filename
    TString filename0 = Output_FileName;
    if(filename0=="NONE"){
        filename0.Form("_%d_%d",Int_t(momentum_ele),Int_t(momentum_ion));
        filename0 = prefix + filename0;
    }
    TString filename1 = Form("%s_1_%d.root",filename0.Data(), Int_t(FileNo));
    if(bXSMode)
        filename1 = Form("%s_pos_%d.root",filename0.Data(), Int_t(FileNo));
    TFile *file1 = new TFile(filename1,"RECREATE");
    TTree *t1 = new TTree("T","T");
    t1->SetDirectory(file1);

    t1->Branch("Q2",&Q2,"data/D");/*{{{*/
    t1->Branch("W",&W,"data/D");
    t1->Branch("Wp",&Wp,"data/D");
    t1->Branch("x",&x,"data/D");
    t1->Branch("y",&y,"data/D");
    t1->Branch("z",&z,"data/D");
    t1->Branch("nu",&nu,"data/D");
    t1->Branch("s",&s,"data/D");
    t1->Branch("pt",&pt,"data/D");
    t1->Branch("gamma",&gamma,"data/D");
    t1->Branch("epsilon", &epsilon,"data/D");
    t1->Branch("rapidity",&rapidity,"data/D");
    t1->Branch("physical",&physical,"data/D");
    t1->Branch("theta_q",&theta_q,"data/D");
    t1->Branch("theta_s",&theta_s,"data/D");
    t1->Branch("phi_h",&phi_h,"data/D");
    t1->Branch("phi_s",&phi_s,"data/D");
    t1->Branch("jacoF",&jacoF,"jacoF/D");
    t1->Branch("isphy_hm",&isphy_hm,"isphy_hm/I");
    t1->Branch("isphy_hp",&isphy_hp,"isphy_hp/I");
    t1->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t1->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
    t1->Branch("dxs_hm_sidis",&dxs_hm_sidis,"dxs_hm_sidis/D");
    t1->Branch("dxs_hp_sidis",&dxs_hp_sidis,"dxs_hp_sidis/D");
    t1->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
    t1->Branch("mom_ele",&mom_ele,"mom_ele/D");
    t1->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
    t1->Branch("mom_had",&mom_had,"mom_had/D");
    t1->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
    t1->Branch("theta_ele",&theta_ele,"theta_ele/D");
    t1->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
    t1->Branch("theta_had",&theta_had,"theta_had/D");
    t1->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
    t1->Branch("phi_ele",&phi_ele,"phi_ele/D");
    t1->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
    t1->Branch("phi_had",&phi_had,"phi_had/D");
    t1->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
    t1->Branch("nsim",&nsim,"nsim/l");
    t1->Branch("dilute_p",&dilute_hp,"data/D");
    t1->Branch("dilute_m",&dilute_hm ,"data/D");
    t1->Branch("px_ele",&px_ele, "px_ele/D");
    t1->Branch("py_ele",&py_ele, "py_ele/D");
    t1->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t1->Branch("E_ele",&E_ele, "E_ele/D");
    t1->Branch("px_had",&px_had, "px_had/D");
    t1->Branch("py_had",&py_had, "py_had/D");
    t1->Branch("pz_had",&pz_had, "pz_had/D");
    t1->Branch("E_had",&E_had, "E_had/D");
    t1->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t1->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t1->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t1->Branch("vx_had",&vx_had, "vx_had/D");
    t1->Branch("vy_had",&vy_had, "vy_had/D");
    t1->Branch("vz_had",&vz_had, "vz_had/D");
    
    t1->Branch("u_pdf", &u_pdf, "u_pdf/D");
    t1->Branch("d_pdf", &d_pdf, "d_pdf/D");
    t1->Branch("s_pdf", &s_pdf, "s_pdf/D");
    t1->Branch("g_pdf", &g_pdf, "g_pdf/D");
    t1->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
    t1->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
    t1->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
    
    t1->Branch("D_fav", &D_fav, "D_fav/D");
    t1->Branch("D_unfav", &D_unfav, "D_unfav/D");
    t1->Branch("D_s", &D_s, "D_s/D");
    t1->Branch("D_g", &D_g, "D_g/D");
    /*}}}*/

    TString filename2 = Form("%s_2_%d.root",filename0.Data(), Int_t(FileNo));
    if(bXSMode)
        filename2 = Form("%s_neg_%d.root",filename0.Data(), Int_t(FileNo));
    TFile *file2 = new TFile(filename2,"RECREATE");
    TTree *t2 = new TTree("T","T");
    t2->SetDirectory(file2);
    if(config=="SoLID" || config=="CLAS12" || config=="EIC" || bXSMode ){ //If it is in bXSModle, t1 save hp events and t2 save hm events
        t2->Branch("Q2",&Q2,"data/D");/*{{{*/
        t2->Branch("W",&W,"data/D");
        t2->Branch("Wp",&Wp,"data/D");
        t2->Branch("x",&x,"data/D");
        t2->Branch("y",&y,"data/D");
        t2->Branch("z",&z,"data/D");
        t2->Branch("nu",&nu,"data/D");
        t2->Branch("s",&s,"data/D");
        t2->Branch("pt",&pt,"data/D");
        t2->Branch("gamma",&gamma,"data/D");
        t2->Branch("epsilon", &epsilon,"data/D");
        t2->Branch("rapidity",&rapidity,"data/D");
        t2->Branch("physical",&physical,"data/D");
        t2->Branch("theta_q",&theta_q,"data/D");
        t2->Branch("theta_s",&theta_s,"data/D");
        t2->Branch("phi_h",&phi_h,"data/D");
        t2->Branch("phi_s",&phi_s,"data/D");
        t2->Branch("jacoF",&jacoF,"jacoF/D");
        t2->Branch("isphy_hm",&isphy_hm,"isphy_hm/I");
        t2->Branch("isphy_hp",&isphy_hp,"isphy_hp/I");
        t2->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
        t2->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
        t2->Branch("dxs_hm_sidis",&dxs_hm_sidis,"dxs_hm_sidis/D");
        t2->Branch("dxs_hp_sidis",&dxs_hp_sidis,"dxs_hp_sidis/D");
        t2->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
        t2->Branch("mom_ele",&mom_ele,"mom_ele/D");
        t2->Branch("mom_had",&mom_had,"mom_had/D");
        t2->Branch("theta_ele",&theta_ele,"theta_ele/D");
        t2->Branch("theta_had",&theta_had,"theta_had/D");
        t2->Branch("phi_ele",&phi_ele,"phi_ele/D");
        t2->Branch("phi_had",&phi_had,"phi_had/D");
        t2->Branch("nsim",&nsim,"nsim/l");
        t2->Branch("dilute_p",&dilute_hp ,"data/D");
        t2->Branch("dilute_m",&dilute_hm ,"data/D");
        t2->Branch("px_ele",&px_ele, "px_ele/D");
        t2->Branch("py_ele",&py_ele, "py_ele/D");
        t2->Branch("pz_ele",&pz_ele, "pz_ele/D");
        t2->Branch("E_ele",&E_ele, "E_ele/D");
        t2->Branch("px_had",&px_had, "px_had/D");
        t2->Branch("py_had",&py_had, "py_had/D");
        t2->Branch("pz_had",&pz_had, "pz_had/D");
        t2->Branch("E_had",&E_had, "E_had/D");
        t2->Branch("vx_ele",&vx_ele, "vx_ele/D");
        t2->Branch("vy_ele",&vy_ele, "vy_ele/D");
        t2->Branch("vz_ele",&vz_ele, "vz_ele/D");
        t2->Branch("vx_had",&vx_had, "vx_had/D");
        t2->Branch("vy_had",&vy_had, "vy_had/D");
        t2->Branch("vz_had",&vz_had, "vz_had/D");
        t2->Branch("u_pdf", &u_pdf, "u_pdf/D");
        t2->Branch("d_pdf", &d_pdf, "d_pdf/D");
        t2->Branch("s_pdf", &s_pdf, "s_pdf/D");
        t2->Branch("g_pdf", &g_pdf, "g_pdf/D");
        t2->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
        t2->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
        t2->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
        t2->Branch("D_fav", &D_fav, "D_fav/D");
        t2->Branch("D_unfav", &D_unfav, "D_unfav/D");
        t2->Branch("D_s", &D_s, "D_s/D");
        t2->Branch("D_g", &D_g, "D_g/D");
        /*}}}*/
    }

    TString filename3 = Form("%s_3_%d.root",filename0.Data(), Int_t(FileNo));
    TFile *file3 = new TFile(filename3,"RECREATE");
    TTree *t3 = new TTree("T","T");
    t3->SetDirectory(file3);

    TString filename4 = Form("%s_4_%d.root",filename0.Data(), Int_t(FileNo));
    TFile *file4 = new TFile(filename4,"RECREATE");
    TTree *t4 = new TTree("T","T");
    t4->SetDirectory(file4);

    if(config=="EIC" ){
        t3->Branch("Q2",&Q2,"data/D");/*{{{*/
        t3->Branch("W",&W,"data/D");
        t3->Branch("Wp",&Wp,"data/D");
        t3->Branch("x",&x,"data/D");
        t3->Branch("y",&y,"data/D");
        t3->Branch("z",&z,"data/D");
        t3->Branch("nu",&nu,"data/D");
        t3->Branch("s",&s,"data/D");
        t3->Branch("pt",&pt,"data/D");
        t3->Branch("gamma",&gamma,"data/D");
        t3->Branch("epsilon", &epsilon,"data/D");
        t3->Branch("rapidity",&rapidity,"data/D");
        t3->Branch("physical",&physical,"data/D");
        t3->Branch("theta_q",&theta_q,"data/D");
        t3->Branch("theta_s",&theta_s,"data/D");
        t3->Branch("phi_h",&phi_h,"data/D");
        t3->Branch("phi_s",&phi_s,"data/D");
        t3->Branch("jacoF",&jacoF,"jacoF/D");
        t3->Branch("isphy_hm",&isphy_hm,"isphy_hm/I");
        t3->Branch("isphy_hp",&isphy_hp,"isphy_hp/I");
        t3->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
        t3->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
        t3->Branch("dxs_hm_sidis",&dxs_hm_sidis,"dxs_hm_sidis/D");
        t3->Branch("dxs_hp_sidis",&dxs_hp_sidis,"dxs_hp_sidis/D");
        t3->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
        t3->Branch("mom_ele",&mom_ele,"mom_ele/D");
        t3->Branch("mom_had",&mom_had,"mom_had/D");
        t3->Branch("theta_ele",&theta_ele,"theta_ele/D");
        t3->Branch("theta_had",&theta_had,"theta_had/D");
        t3->Branch("phi_ele",&phi_ele,"phi_ele/D");
        t3->Branch("phi_had",&phi_had,"phi_had/D");
        t3->Branch("nsim",&nsim,"nsim/l");
        t3->Branch("dilute_p",&dilute_hp ,"data/D");
        t3->Branch("dilute_m",&dilute_hm ,"data/D");
        t3->Branch("px_ele",&px_ele, "px_ele/D");
        t3->Branch("py_ele",&py_ele, "py_ele/D");
        t3->Branch("pz_ele",&pz_ele, "pz_ele/D");
        t3->Branch("E_ele",&E_ele, "E_ele/D");
        t3->Branch("px_had",&px_had, "px_had/D");
        t3->Branch("py_had",&py_had, "py_had/D");
        t3->Branch("pz_had",&pz_had, "pz_had/D");
        t3->Branch("E_had",&E_had, "E_had/D");
        t3->Branch("vx_ele",&vx_ele, "vx_ele/D");
        t3->Branch("vy_ele",&vy_ele, "vy_ele/D");
        t3->Branch("vz_ele",&vz_ele, "vz_ele/D");
        t3->Branch("vx_had",&vx_had, "vx_had/D");
        t3->Branch("vy_had",&vy_had, "vy_had/D");
        t3->Branch("vz_had",&vz_had, "vz_had/D");
        t3->Branch("u_pdf", &u_pdf, "u_pdf/D");
        t3->Branch("d_pdf", &d_pdf, "d_pdf/D");
        t3->Branch("s_pdf", &s_pdf, "s_pdf/D");
        t3->Branch("g_pdf", &g_pdf, "g_pdf/D");
        t3->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
        t3->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
        t3->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
        t3->Branch("D_fav", &D_fav, "D_fav/D");
        t3->Branch("D_unfav", &D_unfav, "D_unfav/D");
        t3->Branch("D_s", &D_s, "D_s/D");
        t3->Branch("D_g", &D_g, "D_g/D");
        /*}}}*/

        t4->Branch("Q2",&Q2,"data/D");/*{{{*/
        t4->Branch("W",&W,"data/D");
        t4->Branch("Wp",&Wp,"data/D");
        t4->Branch("x",&x,"data/D");
        t4->Branch("y",&y,"data/D");
        t4->Branch("z",&z,"data/D");
        t4->Branch("nu",&nu,"data/D");
        t4->Branch("s",&s,"data/D");
        t4->Branch("pt",&pt,"data/D");
        t4->Branch("gamma",&gamma,"data/D");
        t4->Branch("epsilon", &epsilon,"data/D");
        t4->Branch("rapidity",&rapidity,"data/D");
        t4->Branch("physical",&physical,"data/D");
        t4->Branch("theta_q",&theta_q,"data/D");
        t4->Branch("theta_s",&theta_s,"data/D");
        t4->Branch("phi_h",&phi_h,"data/D");
        t4->Branch("phi_s",&phi_s,"data/D");
        t4->Branch("jacoF",&jacoF,"jacoF/D");
        t4->Branch("isphy_hm",&isphy_hm,"isphy_hm/I");
        t4->Branch("isphy_hp",&isphy_hp,"isphy_hp/I");
        t4->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
        t4->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
        t4->Branch("dxs_hm_sidis",&dxs_hm_sidis,"dxs_hm_sidis/D");
        t4->Branch("dxs_hp_sidis",&dxs_hp_sidis,"dxs_hp_sidis/D");
        t4->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
        t4->Branch("mom_ele",&mom_ele,"mom_ele/D");
        t4->Branch("mom_had",&mom_had,"mom_had/D");
        t4->Branch("theta_ele",&theta_ele,"theta_ele/D");
        t4->Branch("theta_had",&theta_had,"theta_had/D");
        t4->Branch("phi_ele",&phi_ele,"phi_ele/D");
        t4->Branch("phi_had",&phi_had,"phi_had/D");
        t4->Branch("nsim",&nsim,"nsim/l");
        t4->Branch("dilute_p",&dilute_hp ,"data/D");
        t4->Branch("dilute_m",&dilute_hm ,"data/D");
        t4->Branch("px_ele",&px_ele, "px_ele/D");
        t4->Branch("py_ele",&py_ele, "py_ele/D");
        t4->Branch("pz_ele",&pz_ele, "pz_ele/D");
        t4->Branch("E_ele",&E_ele, "E_ele/D");
        t4->Branch("px_had",&px_had, "px_had/D");
        t4->Branch("py_had",&py_had, "py_had/D");
        t4->Branch("pz_had",&pz_had, "pz_had/D");
        t4->Branch("E_had",&E_had, "E_had/D");
        t4->Branch("vx_ele",&vx_ele, "vx_ele/D");
        t4->Branch("vy_ele",&vy_ele, "vy_ele/D");
        t4->Branch("vz_ele",&vz_ele, "vz_ele/D");
        t4->Branch("vx_had",&vx_had, "vx_had/D");
        t4->Branch("vy_had",&vy_had, "vy_had/D");
        t4->Branch("vz_had",&vz_had, "vz_had/D");
        t4->Branch("u_pdf", &u_pdf, "u_pdf/D");
        t4->Branch("d_pdf", &d_pdf, "d_pdf/D");
        t4->Branch("s_pdf", &s_pdf, "s_pdf/D");
        t4->Branch("g_pdf", &g_pdf, "g_pdf/D");
        t4->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
        t4->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
        t4->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
       t4->Branch("D_fav", &D_fav, "D_fav/D");
       t4->Branch("D_unfav", &D_unfav, "D_unfav/D");
       t4->Branch("D_s", &D_s, "D_s/D");
       t4->Branch("D_g", &D_g, "D_g/D");
       /*}}}*/
    }

    /*}}}*/

    ofstream pos_gemc, neg_gemc;/*{{{*/
    TString filename_pos = Form("%s_%d_pos.LUND",filename0.Data(), Int_t(FileNo));
    TString filename_neg = Form("%s_%d_neg.LUND",filename0.Data(), Int_t(FileNo));
    if(bLUND){
        pos_gemc.open(filename_pos);
        neg_gemc.open(filename_neg);
    }/*}}}*/

    //Initialize XS Model here/*{{{*/
    //CTEQPDF or EPS09
    SIDIS *sidis = new SIDIS(model);
    sidis->Init(ion_mass, A, Z, particle_flag);

    ////////////////////////////////
    //*** If using EPS09, two CTEQ PDF sets are used for LO and NLO EPS09 sets
    //*** 1->LO need CTEQ6L1,  2->NLO need CTEQ6.1M, default is 2
    //////////
    //sidis->SetEPS09(2);
    ////////////////////////////////
   
    ////////////////////////////////
    /////////////
    //*** If using CTEQPDF, default is CTEQ6.1M, but you can specify the name of the PDF sets
    //*** For exmaple, 4--> CTEQ6L1 (LO), 200-->CTEQ6.1M (NLO)
    //*** See ./cteq-pdf-1.0.4/Cteq6Pdf-2008.txt for details
    /////////////
    //SetCTEQ( mode);
    ////////////////////////////////
       
    ////*** Or if using nCTEQ, add this to specify the associated PDF set of the target
    //sidis->SetLHAPDF(A, Z); 
    ////////////////////////////////
    /*}}}*/

    /*Start to generate events{{{*/
    bool exitcondition=true;	
    while(exitcondition){
        nsim ++;

        /*Generator{{{*/
        //For electron
        vx_ele = gRandom->Uniform(-beamsize_x_ele, beamsize_x_ele);
        vy_ele = gRandom->Uniform(-beamsize_y_ele, beamsize_y_ele);
        vz_ele = vertex_center + gRandom->Uniform(-vertex_length/2.0, vertex_length/2.0);

        phi_gen = gRandom->Uniform(Ph_Min_e/DEG,Ph_Max_e/DEG);
        theta_gen = acos(gRandom->Uniform(cos(Th_Max_e/DEG),cos(Th_Min_e/DEG)));
        mom_gen = gRandom->Uniform(Mom_Min_e, Mom_Max_e);

        mom_gen_ele = mom_gen; theta_gen_ele = theta_gen*DEG; phi_gen_ele = phi_gen*DEG;

        //For hadron, scattered electrons and hardon should come out from the same location
        vx_had = vx_ele;      vy_had = vy_ele;        vz_had = vz_ele;

        phi_gen = gRandom->Uniform(Ph_Min_h/DEG,Ph_Max_h/DEG);
        theta_gen = acos(gRandom->Uniform(cos(Th_Max_h/DEG),cos(Th_Min_h/DEG)));
        mom_gen = gRandom->Uniform(Mom_Min_h, Mom_Max_h);
        mom_gen_had = mom_gen; theta_gen_had = theta_gen*DEG; phi_gen_had = phi_gen*DEG;
        /*}}}*/

        sidis->SetKin(momentum_ele, momentum_ion,/*{{{*/
                mom_gen_ele, theta_gen_ele, phi_gen_ele,
                mom_gen_had, theta_gen_had, phi_gen_had);

        mom_ele = sidis->fMom_ele; theta_ele = sidis->fTheta_ele; phi_ele = sidis->fPhi_ele;
        mom_had = sidis->fMom_had; theta_had = sidis->fTheta_had; phi_had = sidis->fPhi_had;
        theta_q=sidis->fTheta_q;  theta_s=sidis->fTheta_s; phi_s= sidis->fPhi_s; phi_h = sidis->fPhi_h;
        px_ele = sidis->fPx_ele;	py_ele = sidis->fPy_ele;	pz_ele = sidis->fPz_ele; E_ele = sidis->fE_ele;
        px_had = sidis->fPx_had;	py_had = sidis->fPy_had;	pz_had = sidis->fPz_had; E_had = sidis->fE_had;

        x=sidis->fXb; y=sidis->fY; z=sidis->fZ_h; Q2=sidis->fQ2; W=sidis->fW; Wp=sidis->fWp;
        s=sidis->fS; nu=sidis->fNu; pt=sidis->fPt; gamma=sidis->fGamma; epsilon=sidis->fEpsilon;
        rapidity = sidis->fRapidity;
        physical = sidis->fPhysical;
        jacoF=sidis->fJacobF;/*}}}*/

        if(bXSMode){
            /*Generate Events based on XS and also for LUND output{{{*/
            if (x<0.0 || x>1.0 || Q2 <1.0 || W< 2.0) continue;
            //For EIC  
            if( (config=="EIC" && z>0.2&&z<0.9
                        &&((count[0]<number_of_events)|| (count[1]<number_of_events)))
                    ||((config=="SoLID" || config=="CLAS12") && z>0.3&&z<0.7 
                        &&((count[0]<number_of_events)|| (count[1]<number_of_events)))
                    ||(config=="SPECT" && z>0.2&&z<0.9 && count[0]<number_of_events)){

                sidis->CalcXS();/*{{{*/
                dxs_incl = sidis->GetXS_Inclusive();
                dxs_hp = sidis->GetXS_HP();
                dxs_hm = sidis->GetXS_HM();
                dxs_hp_sidis = sidis->GetXS_HP_SIDIS();
                dxs_hm_sidis = sidis->GetXS_HM_SIDIS();
                dilute_hp = sidis->GetDilute_HP();
                dilute_hm = sidis->GetDilute_HM();
                isphy_hp = sidis->IsPhy_HP();
                isphy_hm = sidis->IsPhy_HM();

                 u_pdf = sidis->get_uA();
                d_pdf = sidis->get_dA();
                s_pdf = sidis->get_s();
                g_pdf = sidis->get_g();
                ubar_pdf = sidis->get_ubar();
                dbar_pdf = sidis->get_dbar();
                sbar_pdf = sidis->get_sbar();
                
                D_fav = sidis->get_Dfav();
                D_unfav = sidis->get_Dunfav();
                D_s = sidis->get_Ds();
                D_g = sidis->get_Dg();
                /*}}}*/

                /*LUND For Positive Hadron{{{*/
                //These section will save events based on their XS distributions
                Double_t cdxs_max_rndm_hp = cdxs_max * gRandom->Uniform(0, 1);
                if(bLUND&&dxs_hp_sidis > cdxs_max_rndm_hp){
                    //Header:      1#part. 2#x 3#z 4#pt 5#Pol 6#Q2 7#W 8#cxs 9#phi_s 10#phi_h
                    pos_gemc<<Form("    %2d \t %10.4e \t %10.4e \t %10.4e \t %4.3f \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",/*{{{*/
                            2, //ele+had
                            x,
                            z,
                            pt,
                            1.0, //pol = 1.0 for now
                            Q2,
                            W,
                            phi_s,
                            phi_h,
                            dxs_hp_sidis						
                            )<<endl;

                    //electron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
                    pos_gemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
                            1, //index
                            -1.0,//charge
                            1, //=1 for active 
                            11,//pid
                            0,// parent pid, not in used now
                            0,// doughter for decay bookkeeping, not in used now
                            px_ele,
                            py_ele,
                            pz_ele,
                            E_ele,
                            0.0005, //mass not in used	
                            vx_ele, //vx
                            vy_ele, //vx
                            vz_ele  //vx
                            )<<endl;
                    //hadron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
                    pos_gemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
                            2, //index
                            charge_pos,//charge
                            1, //=1 for active 
                            pid_pos,//pid
                            0,// parent pid, not in used now
                            0,// doughter for decay bookkeeping, not in used now
                            px_had,
                            py_had,
                            pz_had,
                            E_had, 
                            mass_had, //mass not in used
                            vx_had, //vx
                            vy_had, //vx
                            vz_had  //vx
                            )<<endl;/*}}}*/
                }
                /*}}}*/

                /*LUND For Negative Hadron{{{*/
                //These section will save events based on their XS distributions
                Double_t cdxs_max_rndm_hm = cdxs_max * gRandom->Uniform(0, 1);
                if(bLUND&&dxs_hm_sidis > cdxs_max_rndm_hm){
                    //Header:      1#part. 2#x 3#z 4#pt 5#Pol 6#Q2 7#W 8#cxs 9#phi_s 10#phi_h
                    neg_gemc<<Form("    %2d \t %10.4e \t %10.4e \t %10.4e \t %4.3f \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",/*{{{*/
                            2, //ele+had
                            x,
                            z,
                            pt,
                            1.0, //pol = 1.0 for now
                            Q2,
                            W,
                            phi_s,
                            phi_h,
                            dxs_hm_sidis						
                            )<<endl;

                    //electron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
                    neg_gemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
                            1, //index
                            -1.0,//charge
                            1, //=1 for active 
                            11,//pid
                            0,// parent pid, not in used now
                            0,// doughter for decay bookkeeping, not in used now
                            px_ele,
                            py_ele,
                            pz_ele,
                            E_ele,
                            0.0005, //mass not in used	
                            vx_ele, //vx
                            vy_ele, //vx
                            vz_ele  //vx
                            )<<endl;
                    //hadron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
                    neg_gemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
                            2, //index
                            charge_neg,//charge
                            1, //=1 for active 
                            pid_neg,//pid
                            0,// parent pid, not in used now
                            0,// doughter for decay bookkeeping, not in used now
                            px_had,
                            py_had,
                            pz_had,
                            E_had, 
                            mass_had, //mass not in used
                            vx_had, //vx
                            vy_had, //vx
                            vz_had  //vx
                            )<<endl;/*}}}*/
                }
                /*}}}*/

                /*Save PI+ ROOT file based on XS distribution{{{*/
                if(dxs_hp_sidis > cdxs_max_rndm_hp){
                    t1->Fill();
                    //Just save events in one root files in this case
                    count[0] ++;//cout << 0 << " " << count[0] << endl;
                    cout << count[0] <<"\r";
                    Nsim1 = nsim;
                }
                /*}}}*/

                /*Save PI- ROOT file based on XS distribution{{{*/
                if(dxs_hm_sidis > cdxs_max_rndm_hm){
                    t2->Fill();
                    //Just save events in one root files in this case
                    count[1] ++;//cout << 0 << " " << count[0] << endl;
                    cout << count[1] <<"\r";
                    Nsim2 = nsim;
                }
                /*}}}*/
                cout << count[0] << "\t" << count[1] << "\t" << count[2] << "\t" << count[3] << "\r";
            }
            //judging exitcondition/*{{{*/
            if (config=="EIC"||config=="SoLID"||config=="CLAS12") {
                if (count[0] < number_of_events || count[1] < number_of_events 
                   ) exitcondition=true;
                else exitcondition=false;
            } 
            else if(config=="SPECT") {
                if (count[0] < number_of_events) exitcondition=true;
                else exitcondition=false;
            } /*}}}*/
            /*}}}*/
        }else{
            /*Generate Events Uniformly{{{*/
            if (x<0.0 || x>1.0 || Q2 <1.0 || W< 2.0) continue;
            //if (x<0.05 || x>0.3 || Q2 <1.0 || W< 2.0) continue;
            if ( (config=="EIC" && z>0.2&&z<0.9//&&y>0.05&&y<0.8
                        &&(   (count[0]<number_of_events&&pt<=1.0&&Q2<=Q2_CutOff)
                            ||(count[1]<number_of_events&&pt>1.0&&Q2<=Q2_CutOff)
                            ||(count[2]<number_of_events&&pt<=1.0&&Q2>Q2_CutOff)
                            ||(count[3]<number_of_events&&pt>1.0&&Q2>Q2_CutOff)))
                    ||((config=="SoLID"||config=="CLAS12") && z>0.3&&z<0.7 
                        &&(   (count[0]<number_of_events&&pt<=1.0)
                            ||(count[1]<number_of_events&&pt>1.0)))
                    ||(config=="SPECT" && z>0.2&&z<0.9 && count[0]<number_of_events))
            {

                sidis->CalcXS();/*{{{*/
                dxs_incl = sidis->GetXS_Inclusive();
                dxs_hp = sidis->GetXS_HP();
                dxs_hm = sidis->GetXS_HM();
                dxs_hp_sidis = sidis->GetXS_HP_SIDIS();
                dxs_hm_sidis = sidis->GetXS_HM_SIDIS();
                dilute_hp = sidis->GetDilute_HP();
                dilute_hm = sidis->GetDilute_HM();
                
                u_pdf = sidis->get_uA();
                d_pdf = sidis->get_dA();
                s_pdf = sidis->get_s();
                g_pdf = sidis->get_g();
                ubar_pdf = sidis->get_ubar();
                dbar_pdf = sidis->get_dbar();
                sbar_pdf = sidis->get_sbar();
                
                D_fav = sidis->get_Dfav();
                D_unfav = sidis->get_Dunfav();
                D_s = sidis->get_Ds();
                D_g = sidis->get_Dg();
                /*}}}*/

                /*LUND For Positive Hadron{{{*/
                //These section will save events based on their XS distributions
                Double_t cdxs_max_rndm_hp = cdxs_max * gRandom->Uniform(0, 1);
                if(bLUND&&dxs_hp_sidis > cdxs_max_rndm_hp){
                    //Header:      1#part. 2#x 3#z 4#pt 5#Pol 6#Q2 7#W 8#cxs 9#phi_s 10#phi_h
                    pos_gemc<<Form("    %2d \t %10.4e \t %10.4e \t %10.4e \t %4.3f \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",/*{{{*/
                            2, //ele+had
                            x,
                            z,
                            pt,
                            1.0, //pol = 1.0 for now
                            Q2,
                            W,
                            phi_s,
                            phi_h,
                            dxs_hp_sidis						
                            )<<endl;

                    //electron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
                    pos_gemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
                            1, //index
                            -1.0,//charge
                            1, //=1 for active 
                            11,//pid
                            0,// parent pid, not in used now
                            0,// doughter for decay bookkeeping, not in used now
                            px_ele,
                            py_ele,
                            pz_ele,
                            E_ele,
                            0.0005, //mass not in used	
                            vx_ele, //vx
                            vy_ele, //vx
                            vz_ele  //vx
                            )<<endl;
                    //hadron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
                    pos_gemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
                            2, //index
                            charge_pos,//charge
                            1, //=1 for active 
                            pid_pos,//pid
                            0,// parent pid, not in used now
                            0,// doughter for decay bookkeeping, not in used now
                            px_had,
                            py_had,
                            pz_had,
                            E_had, 
                            mass_had, //mass not in used
                            vx_had, //vx
                            vy_had, //vx
                            vz_had  //vx
                            )<<endl;/*}}}*/
                }
                /*}}}*/

                /*LUND For Negative Hadron{{{*/
                //These section will save events based on their XS distributions
                Double_t cdxs_max_rndm_hm = cdxs_max * gRandom->Uniform(0, 1);
                if(bLUND&&dxs_hm_sidis > cdxs_max_rndm_hm){
                    //Header:      1#part. 2#x 3#z 4#pt 5#Pol 6#Q2 7#W 8#cxs 9#phi_s 10#phi_h
                    neg_gemc<<Form("    %2d \t %10.4e \t %10.4e \t %10.4e \t %4.3f \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",/*{{{*/
                            2, //ele+had
                            x,
                            z,
                            pt,
                            1.0, //pol = 1.0 for now
                            Q2,
                            W,
                            phi_s,
                            phi_h,
                            dxs_hm_sidis						
                            )<<endl;

                    //electron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
                    neg_gemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
                            1, //index
                            -1.0,//charge
                            1, //=1 for active 
                            11,//pid
                            0,// parent pid, not in used now
                            0,// doughter for decay bookkeeping, not in used now
                            px_ele,
                            py_ele,
                            pz_ele,
                            E_ele,
                            0.0005, //mass not in used	
                            vx_ele, //vx
                            vy_ele, //vx
                            vz_ele  //vx
                            )<<endl;
                    //hadron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
                    neg_gemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
                            2, //index
                            charge_neg,//charge
                            1, //=1 for active 
                            pid_neg,//pid
                            0,// parent pid, not in used now
                            0,// doughter for decay bookkeeping, not in used now
                            px_had,
                            py_had,
                            pz_had,
                            E_had, 
                            mass_had, //mass not in used
                            vx_had, //vx
                            vy_had, //vx
                            vz_had  //vx
                            )<<endl;/*}}}*/
                }
                /*}}}*/

                /*Fill ROOT{{{*/
                if(config=="SoLID"||config=="CLAS12"||config=="EIC"){
                    if(Q2<=Q2_CutOff&&pt<=1.0){
                        t1->Fill();
                        count[0] ++;//cout << 0 << " " << count[0] << endl;
                        Nsim1 = nsim;
                    }
                    if (Q2<=Q2_CutOff&&pt>1.0){
                        t2->Fill();
                        count[1] ++;
                        Nsim2 = nsim;
                    }

                    if (config=="EIC"&&Q2>Q2_CutOff&&pt<=1.0){
                        t3->Fill();
                        count[2] ++;//cout << 2 << " " << count[2] << endl;
                        Nsim3 = nsim;
                    }
                    if (config=="EIC"&&Q2>Q2_CutOff&&pt>1.0){
                        t4->Fill();
                        count[3] ++;//cout << 3 << " " << count[3] <<  endl;
                        Nsim4 = nsim;
                    }
                }
                if(config=="SPECT"){
                    t1->Fill();
                    count[0] ++;//cout << 0 << " " << count[0] << endl;
                    Nsim1 = nsim;
                }
                /*}}}*/
                cout << count[0] << "\t" << count[1] << "\t" << count[2] << "\t" << count[3] << "\r";
            }

            //judging exitcondition/*{{{*/
            if (config=="EIC") {
                if (count[0] < number_of_events || count[1] < number_of_events 
                        || count[2] < number_of_events || count[3] < number_of_events) exitcondition=true;
                else exitcondition=false;
            } 
            else if (config=="SoLID"||config=="CLAS12") {
                if (count[0] < number_of_events || count[1] < number_of_events) exitcondition=true;
                else exitcondition=false;
            } 
            else if(config=="SPECT") {
                if (count[0] < number_of_events) exitcondition=true;
                else exitcondition=false;
            } /*}}}*/
            /*}}}*/ 
        }
    }
    /*}}}*/

    cout << count[0] << "\t" << count[1] << "\t" << count[2] << "\t" << count[3] << endl;

    file1->Write();/*{{{*/
    file1->Close();
    if(config=="SoLID" || config=="CLAS12" || config=="EIC" ){
        file2->Write();
        file2->Close();
    }
    if(config=="EIC" ){
        file3->Write();
        file3->Close();
        file4->Write();
        file4->Close();
    }
    delete sidis;

    if(bLUND){
        pos_gemc.close();
        neg_gemc.close();
    }
    /*}}}*/

    //Generate weights for uniform distributed events
    if(!bXSMode){/*{{{*/
        double weight_hp=1.0e-34;
        double weight_hm=1.0e-34;
        double weight_in=1.0e-34;
        /*Generate weights for file1{{{*/
        cout<<"--- Now insert weights to Root-file #1"<<filename1.Data()<<endl;
        TFile *f1 = new TFile(filename1.Data(), "update");
        TTree *T1=(TTree*) f1->Get("T");
        ULong64_t N1=T1->GetEntries();
        T1->SetBranchAddress("dxs_hp_sidis",&dxs_hp_sidis);
        T1->SetBranchAddress("dxs_hm_sidis",&dxs_hm_sidis);
        T1->SetBranchAddress("dxs_incl",&dxs_incl);
        ULong64_t nsim1; 
        T1->SetBranchAddress("nsim",&nsim1);
        T1->GetEntry(N1-1);          //get nsim for this rootfile
        ULong64_t Nsim10=nsim1;
        cout<<Form("--- Root#1: Nsim = %lld / %lld", Nsim1, Nsim10)<<endl;

        TBranch *branch_weight_in1=T1->Branch("weight_in",&weight_in,"weight_in/D");
        TBranch *branch_weight_hp1=T1->Branch("weight_hp",&weight_hp,"weight_hp/D");
        TBranch *branch_weight_hm1=T1->Branch("weight_hm",&weight_hm,"weight_hm/D");
        cout<<Form("---Filling weights foor ROOT#1, Nsim=%lld, Phase_space = %f", Nsim1, Phase_space)<<endl;
        for(ULong64_t i=0;i<N1;i++){
            T1->GetEntry(i);
            //warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
            weight_in=dxs_incl*electron_phase_space/Nsim1;   
            weight_hp=dxs_hp_sidis*Phase_space/Nsim1;   
            weight_hm=dxs_hm_sidis*Phase_space/Nsim1;

            if((weight_in)<1e-34) weight_in=1e-34;
            if((weight_hp)<1e-34) weight_hp=1e-34;
            if((weight_hm)<1e-34) weight_hm=1e-34;
            if(isnan(weight_in)) weight_in=1e-34;
            if(isnan(weight_hp)) weight_hp=1e-34;
            if(isnan(weight_hm)) weight_hm=1e-34;
            if(isinf(weight_in)) weight_in=1e-34;
            if(isinf(weight_hp)) weight_hp=1e-34;
            if(isinf(weight_hm)) weight_hm=1e-34;


            branch_weight_in1->Fill();
            branch_weight_hp1->Fill();
            branch_weight_hm1->Fill();
        }
        T1->Write("",TObject::kOverwrite);
        f1->Close();
        /*}}}*/

        if(config=="SoLID" || config=="CLAS12" || config=="EIC" ){
            /*Generate weights for file2{{{*/
            cout<<"--- Now insert weights to Root-file #2"<<filename2.Data()<<endl;
            TFile *f2 = new TFile(filename2.Data(), "update");
            TTree *T2=(TTree*) f2->Get("T");
            ULong64_t N2=T2->GetEntries();
            T2->SetBranchAddress("dxs_incl",&dxs_incl);
            T2->SetBranchAddress("dxs_hp_sidis",&dxs_hp_sidis);
            T2->SetBranchAddress("dxs_hm_sidis",&dxs_hm_sidis);
            ULong64_t nsim2; 
            T2->SetBranchAddress("nsim",&nsim2);
            T2->GetEntry(N2-1);          //get nsim for this rootfile
            ULong64_t Nsim20=nsim2;
            cout<<Form("--- Root#2: Nsim = %lld / %lld", Nsim2, Nsim20)<<endl;

            TBranch *branch_weight_in2=T2->Branch("weight_in",&weight_in,"weight_in/D");
            TBranch *branch_weight_hp2=T2->Branch("weight_hp",&weight_hp,"weight_hp/D");
            TBranch *branch_weight_hm2=T2->Branch("weight_hm",&weight_hm,"weight_hm/D");
            cout<<Form("---Filling weights foor ROOT#2, Nsim=%lld, Phase_space = %f", Nsim2, Phase_space)<<endl;
            for(ULong64_t i=0;i<N2;i++){
                T2->GetEntry(i);
                //warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
                weight_in=dxs_incl*electron_phase_space/Nsim2;   
                weight_hp=dxs_hp_sidis*Phase_space/Nsim2;   
                weight_hm=dxs_hm_sidis*Phase_space/Nsim2;

                if((weight_in)<1e-34) weight_in=1e-34;
                if((weight_hp)<1e-34) weight_hp=1e-34;
                if((weight_hm)<1e-34) weight_hm=1e-34;
                if(isnan(weight_in)) weight_in=1e-34;
                if(isnan(weight_hp)) weight_hp=1e-34;
                if(isnan(weight_hm)) weight_hm=1e-34;
                if(isinf(weight_in)) weight_in=1e-34;
                if(isinf(weight_hp)) weight_hp=1e-34;
                if(isinf(weight_hm)) weight_hm=1e-34;
                branch_weight_in2->Fill();
                branch_weight_hp2->Fill();
                branch_weight_hm2->Fill();
            }
            T2->Write("",TObject::kOverwrite);
            f2->Close();
            /*}}}*/
        }

        if(config=="EIC" ){
            /*Generate weights for file3{{{*/
            cout<<"--- Now insert weights to Root-file #3"<<filename3.Data()<<endl;
            TFile *f3 = new TFile(filename3.Data(), "update");
            TTree *T3=(TTree*) f3->Get("T");
            ULong64_t N3=T3->GetEntries();
            T3->SetBranchAddress("dxs_incl",&dxs_incl);
            T3->SetBranchAddress("dxs_hp_sidis",&dxs_hp_sidis);
            T3->SetBranchAddress("dxs_hm_sidis",&dxs_hm_sidis);
            ULong64_t nsim3; 
            T3->SetBranchAddress("nsim",&nsim3);
            T3->GetEntry(N3-1);          //get nsim for this rootfile
            ULong64_t Nsim30=nsim3;
            cout<<Form("--- Root#3: Nsim = %lld / %lld", Nsim3, Nsim30)<<endl;


            TBranch *branch_weight_in3=T3->Branch("weight_in",&weight_in,"weight_in/D");
            TBranch *branch_weight_hp3=T3->Branch("weight_hp",&weight_hp,"weight_hp/D");
            TBranch *branch_weight_hm3=T3->Branch("weight_hm",&weight_hm,"weight_hm/D");
            cout<<Form("---Filling weights foor ROOT#3, Nsim=%lld, Phase_space = %f", Nsim3, Phase_space)<<endl;
            for(ULong64_t i=0;i<N3;i++){
                T3->GetEntry(i);
                //warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
                weight_in=dxs_incl*electron_phase_space/Nsim3;   
                weight_hp=dxs_hp_sidis*Phase_space/Nsim3;   
                weight_hm=dxs_hm_sidis*Phase_space/Nsim3;

                if((weight_in)<1e-34) weight_in=1e-34;
                if((weight_hp)<1e-34) weight_hp=1e-34;
                if((weight_hm)<1e-34) weight_hm=1e-34;
                if(isnan(weight_in)) weight_in=1e-34;
                if(isnan(weight_hp)) weight_hp=1e-34;
                if(isnan(weight_hm)) weight_hm=1e-34;
                if(isinf(weight_in)) weight_in=1e-34;
                if(isinf(weight_hp)) weight_hp=1e-34;
                if(isinf(weight_hm)) weight_hm=1e-34;
                
                branch_weight_in3->Fill();
                branch_weight_hp3->Fill();
                branch_weight_hm3->Fill();
            }
            T3->Write("",TObject::kOverwrite);
            f3->Close();
            /*}}}*/

            /*Generate weights for file4{{{*/
            cout<<"--- Now insert weights to Root-file #4"<<filename4.Data()<<endl;
            TFile *f4 = new TFile(filename4.Data(), "update");
            TTree *T4=(TTree*) f4->Get("T");
            ULong64_t N4=T4->GetEntries();
            T4->SetBranchAddress("dxs_incl",&dxs_incl);
            T4->SetBranchAddress("dxs_hp_sidis",&dxs_hp_sidis);
            T4->SetBranchAddress("dxs_hm_sidis",&dxs_hm_sidis);
            ULong64_t nsim4; 
            T4->SetBranchAddress("nsim",&nsim4);
            T4->GetEntry(N4-1);          //get nsim for this rootfile
            ULong64_t Nsim40=nsim4;
            cout<<Form("--- Root#4: Nsim = %lld / %lld", Nsim4, Nsim40)<<endl;

            TBranch *branch_weight_in4=T4->Branch("weight_in",&weight_in,"weight_in/D");
            TBranch *branch_weight_hp4=T4->Branch("weight_hp",&weight_hp,"weight_hp/D");
            TBranch *branch_weight_hm4=T4->Branch("weight_hm",&weight_hm,"weight_hm/D");
            cout<<Form("---Filling weights foor ROOT#4, Nsim=%lld, Phase_space = %f", Nsim4, Phase_space)<<endl;
            for(ULong64_t i=0;i<N4;i++){
                T4->GetEntry(i);
                //warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
                weight_in=dxs_incl*electron_phase_space/Nsim4;   
                weight_hp=dxs_hp_sidis*Phase_space/Nsim4;   
                weight_hm=dxs_hm_sidis*Phase_space/Nsim4;

                if((weight_in)<1e-34) weight_in=1e-34;
                if((weight_hp)<1e-34) weight_hp=1e-34;
                if((weight_hm)<1e-34) weight_hm=1e-34;
                if(isnan(weight_in)) weight_in=1e-34;
                if(isnan(weight_hp)) weight_hp=1e-34;
                if(isnan(weight_hm)) weight_hm=1e-34;
                if(isinf(weight_in)) weight_in=1e-34;
                if(isinf(weight_hp)) weight_hp=1e-34;
                if(isinf(weight_hm)) weight_hm=1e-34;

                branch_weight_in4->Fill();
                branch_weight_hp4->Fill();
                branch_weight_hm4->Fill();
            }
            T4->Write("",TObject::kOverwrite);
            f4->Close();
            /*}}}*/
        }
    }
    /*}}}*/

    return 0;
    }

    /*Init(){{{*/
    void Init (const TString kInputFile){
        const Int_t CHAR_LEN = 1000;

        cout<<"============================================================"<<endl;
        cout<<"&& Initializing Parameters from "<<kInputFile<<" ..."<<endl;
        int i,j,k;
        vector<TString> inputdata;
        /*Read INPUTfile{{{*/
        FILE* INPUTfile;
        INPUTfile=fopen(kInputFile.Data(),"r");
        char buf[CHAR_LEN];
        char data[CHAR_LEN];
        while ( fgets(buf,CHAR_LEN,INPUTfile) )
        {
            i=0;
            while ( buf[i]==' '|| buf[i]=='\t' )
            {
                i++;
            }
            if ( buf[i]!='#' )
            {
                j=0;
                while ( buf[i]!='#' && buf[i]!='\0' && buf[i]!='\t' && buf[i]!='\n' )
                {
                    if( buf[i]!=' ' && buf[i]!='\t' && buf[i]!='\n' )
                        data[j]=buf[i];
                    i++; j++;
                }
                data[j]='\0';
                while ( data[--j]==' ' || data[j]=='\t'  || data[j]=='\n' )
                {
                    //remove space or tab at the end of data
                    data[j]='\0';
                }
                inputdata.push_back(data);
            }
            //else it's comment, skipped
        }

        fclose(INPUTfile);
        /*}}}*/

        /*Set Global Value{{{*/
        k=0;
        A=atoi(inputdata[k++]);
        Z=atoi(inputdata[k++]);
        particle_flag=atoi(inputdata[k++]);
        momentum_ele=atof(inputdata[k++]);
        momentum_ion=atof(inputdata[k++]);
        number_of_events=atoi(inputdata[k++]);
        FileNo=atoi(inputdata[k++]);
        config= inputdata[k++];
        model= inputdata[k++];
        cdxs_max=atof(inputdata[k++]);
        bLUND=atoi(inputdata[k++]);
        bXSMode=atoi(inputdata[k++]);
        Output_FileName= inputdata[k++];
        bDebug=atoi(inputdata[k++]);

        //A=atoi(inputdata[k++].c_str());
        //Z=atoi(inputdata[k++].c_str());
        //particle_flag=atoi(inputdata[k++].c_str());
        //momentum_ele=atof(inputdata[k++].c_str());
        //momentum_ion=atof(inputdata[k++].c_str());
        //FileNo=atoi(inputdata[k++].c_str());
        //number_of_events=atoi(inputdata[k++].c_str());
        //config= inputdata[k++].c_str();
        //model= inputdata[k++].c_str();
        //cdxs_max=atof(inputdata[k++].c_str());
        //bLUND=atoi(inputdata[k++].c_str());
        //bXSMode=atoi(inputdata[k++].c_str());
        //Output_FileName= inputdata[k++].c_str();
        //bDebug=atoi(inputdata[k++].c_str());
        /*}}}*/

        cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
        cout<<"---- A = "<<A<<endl;
        cout<<"---- Z = "<<Z<<endl;
        cout<<"---- Partile( 1-->Pion, 2-->Kaon) = ";
        if(particle_flag==1) cout<<"Pion"<<endl;
        else if(particle_flag==2) cout<<"Kaon"<<endl;
        else cout<<"*** ERROR, a wrong particle flag in your input-file!!!"<<endl;
        cout<<"---- P_e = "<<momentum_ele<<" GeV"<<endl;
        cout<<"---- P_A = "<<momentum_ion<<" GeV"<<endl;
        cout<<"---- #Events = "<<number_of_events<<endl;
        cout<<"---- File# = "<< FileNo;
        if(FileNo==0) cout<<" (from command line, e.g. ./GetSIDIS input.data FileNo)"<<endl;
        else cout<<endl;
        cout<<"---- Configs = "<<config.Data()<<endl;
        cout<<"---- Model = "<<model.Data()<<endl;
        cout<<"---- Save to LUND? = "<<bLUND<<endl;
        cout<<"---- Save in XSMode? = "<<bXSMode<<endl;
        cout<<"---- Rename files to (*.LUND, *_0.root)= "<<Output_FileName;
        if(Output_FileName=="NONE") cout <<" (use the default name)"<<endl;
        else cout<<endl;
        cout<<"---- Debug? = "<<bDebug<<endl;
        cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;

        inputdata.clear();
        cerr<<"&& Initialization done!"<<endl;
        cout<<"============================================================"<<endl;
    }
    /*}}}*/

    /*double GetIonMass(const int A, const int Z){{{*/
    double GetIonMass(const int A, const int Z){
        double ion_mass = 0.0;
        //if(A==1 && Z==1)//Hydrogen
            //ion_mass = mass_p;
        //else if(A==1 && Z==0)//Neutron
            //ion_mass = mass_n;
        //else if(A==2 && Z==1)//Deutron
            //ion_mass = mass_u*2.014102;
        //else if(A==3 && Z==1)//Tritium
            //ion_mass = mass_u*3.016049;
        //else if(A==3 && Z==2)//He3
            //ion_mass = mass_u*3.016029;
        //else if(A==4 && Z==2)//He4
            //ion_mass = mass_u*4.002602;
        //else if(A==7 && Z==3)//Lithium
            //ion_mass = mass_u*6.941000;
        //else if(A==9 && Z==4)//Be9
            //ion_mass = mass_u*9.012182;
        //else if(A==10 && Z==5)//Boron
            //ion_mass = mass_u*10.811000;
        //else if(A==12 && Z==6)//Carbon
            //ion_mass = mass_u*12.010700;
        //else if(A==14 && Z==7)//Nitrogen
            //ion_mass = mass_u*14.0067;
        //else if(A==16 && Z==8)//Oxygen
            //ion_mass = mass_u*15.9994;
        //else if(A==27 && Z==13)//Al
            //ion_mass = mass_u*26.981539;
        //else if(A==40 && Z==20)//C40
            //ion_mass = mass_u*40.07800;
        //else if(A==48 && Z==20)//C40
            //ion_mass = mass_u*47.952534;
        //else if(A==56 && Z==26)//Fe
            //ion_mass = mass_u*55.84500;
        //else if(A==64 && Z==29)//Cu
            //ion_mass = mass_u*63.546;
        //else if(A==197 && Z==79)//Cu
            //ion_mass = mass_u*196.966569;
////        if the target is not in the list above, simply add the total mass of protons and neutrons
////        but in general you can look at the PERIODIC TABLE for new target you want to add
        //else
            ion_mass = Z*mass_p+(A-Z)*mass_n;

        return ion_mass;
    }
    /*}}}*/
