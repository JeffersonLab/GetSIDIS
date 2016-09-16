//////////////////////////////////////////////////////
// SIDIS Events Generators for SoLID or EIC         //
//                                                  //
//Note: Basically the same as Xin Qian's "collider" //
//      but the model is coded in "SIDIS.h"         //
//  -- Zhihong Ye, 06/10/2014                       //
//////////////////////////////////////////////////////
#include "GetSIDIS.h"
//#include "SIDIS.h"
#include "SIDIS_Lite.h" //this version doesn't include LHAPDF

int main(Int_t argc, char *argv[]){
    cout<<endl;
    cout<<"oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo oO0Oo"<<endl;
    cout<<"oO0 SIDIS Events Generators for SoLID or EIC or Spectrometer  0Oo//"<<endl;
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

    Double_t Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon, jacoF;

    if (config != "EIC" && config != "SoLID" && config != "SPECT"){
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
    Double_t dxs_incl,dxs_hm,dxs_hp,dilute_hp,dilute_hm;
    Double_t px_ele, py_ele,pz_ele, px_had, py_had, pz_had, E0_ele,E0_had;
    Int_t nsim = 0;

    //For Beam Position and Vertex info
    Double_t vx_ele, vy_ele, vz_ele, vx_had, vy_had, vz_had;
    double beamsize_x_ele=0.0, beamsize_y_ele=0.0;
    double vertex_length =0.0, vertex_center=0.0;

    //The idea is to generate a phase-space which is slightly larger than the actual one
    Double_t Mom_Max_e = 0.0, Mom_Min_e = 0.0, Mom_Max_h = 0.0,Mom_Min_h = 0.0;
    Double_t Th_Max_e = 0.0,Th_Min_e = 0.0,Th_Max_h = 0.0, Th_Min_h = 0.0;
    Double_t Ph_Max_e = 0.0,Ph_Min_e = 0.0,Ph_Max_h = 0.0, Ph_Min_h = 0.0;
    if(config=="SoLID" ){
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

    }
    //A rough guess but people claim EIC to be a full-acceptance device!
    else if(config=="EIC" ){
        Mom_Min_e = EIC_Mom_Min_e;  Mom_Max_e =  momentum_ele * 3.0; 
        Mom_Min_h = EIC_Mom_Min_h;  Mom_Max_h = EIC_Mom_Max_h;
        Th_Min_e = EIC_Th_Min_e; Th_Max_e = EIC_Th_Max_e; 
        Th_Min_h = EIC_Th_Min_h; Th_Max_h = EIC_Th_Max_h;    
        Ph_Min_e = EIC_Ph_Min_e; Ph_Max_e = EIC_Ph_Max_e; 
        Ph_Min_h = EIC_Ph_Min_h; Ph_Max_h = EIC_Ph_Max_h;    

        beamsize_x_ele = EIC_BeamSizeX_ele;
        beamsize_y_ele = EIC_BeamSizeY_ele;
        vertex_length = EIC_Vertex_Length;
        vertex_center = EIC_Vertex_Center;
    }
    else if(config=="SPECT"){
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
    }

    Double_t electron_phase_space =(cos(Th_Min_e/DEG) - cos(Th_Max_e/DEG))*(Ph_Max_e/DEG - Ph_Min_e/DEG)*(Mom_Max_e - Mom_Min_e);
    Double_t hadron_phase_space   =(cos(Th_Min_h/DEG) - cos(Th_Max_h/DEG))*(Ph_Max_h/DEG - Ph_Min_h/DEG)*(Mom_Max_h - Mom_Min_h);
    Double_t Phase_space=electron_phase_space*hadron_phase_space;           //electron*hadron phase space eg, for electron: delta_cos_theta*delta_phi*delta_energy
    cout<<" -- For Config="<<config<<" Phase_space: "<<electron_phase_space<<"	"<<hadron_phase_space<<"	"<<Phase_space<<endl;

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
    t1->Branch("theta_q",&theta_q,"data/D");
    t1->Branch("theta_s",&theta_s,"data/D");
    t1->Branch("phi_h",&phi_h,"data/D");
    t1->Branch("phi_s",&phi_s,"data/D");
    t1->Branch("jacoF",&jacoF,"jacoF/D");
    t1->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t1->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
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
    t1->Branch("nsim",&nsim,"nsim/I");
    t1->Branch("dilute_p",&dilute_hp,"data/D");
    t1->Branch("dilute_m",&dilute_hm ,"data/D");
    t1->Branch("px_ele",&px_ele, "px_ele/D");
    t1->Branch("py_ele",&py_ele, "py_ele/D");
    t1->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t1->Branch("E0_ele",&E0_ele, "E0_ele/D");
    t1->Branch("px_had",&px_had, "px_had/D");
    t1->Branch("py_had",&py_had, "py_had/D");
    t1->Branch("pz_had",&pz_had, "pz_had/D");
    t1->Branch("E0_had",&E0_had, "E0_had/D");
    t1->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t1->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t1->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t1->Branch("vx_had",&vx_had, "vx_had/D");
    t1->Branch("vy_had",&vy_had, "vy_had/D");
    t1->Branch("vz_had",&vz_had, "vz_had/D");/*}}}*/

    TString filename2 = Form("%s_2_%d.root",filename0.Data(), Int_t(FileNo));
    if(bXSMode)
        filename2 = Form("%s_neg_%d.root",filename0.Data(), Int_t(FileNo));
    TFile *file2 = new TFile(filename2,"RECREATE");
    TTree *t2 = new TTree("T","T");
    t2->SetDirectory(file2);
    if(config=="SoLID" || config=="EIC" || bXSMode ){ //If it is in bXSModle, t1 save hp events and t2 save hm events
        t2->Branch("Q2",&Q2,"data/D");/*{{{*/
        t2->Branch("W",&W,"data/D");
        t2->Branch("Wp",&Wp,"data/D");
        t2->Branch("x",&x,"data/D");
        t2->Branch("y",&y,"data/D");
        t2->Branch("z",&z,"data/D");
        t2->Branch("nu",&nu,"data/D");
        t2->Branch("s",&s,"data/D");
        t2->Branch("pt",&pt,"data/D");
        t2->Branch("theta_q",&theta_q,"data/D");
        t2->Branch("theta_s",&theta_s,"data/D");
        t2->Branch("phi_h",&phi_h,"data/D");
        t2->Branch("phi_s",&phi_s,"data/D");
        t2->Branch("jacoF",&jacoF,"jacoF/D");
        t2->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
        t2->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
        t2->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
        t2->Branch("mom_ele",&mom_ele,"mom_ele/D");
        t2->Branch("mom_had",&mom_had,"mom_had/D");
        t2->Branch("theta_ele",&theta_ele,"theta_ele/D");
        t2->Branch("theta_had",&theta_had,"theta_had/D");
        t2->Branch("phi_ele",&phi_ele,"phi_ele/D");
        t2->Branch("phi_had",&phi_had,"phi_had/D");
        t2->Branch("nsim",&nsim,"nsim/I");
        t2->Branch("dilute_p",&dilute_hp ,"data/D");
        t2->Branch("dilute_m",&dilute_hm ,"data/D");
        t2->Branch("px_ele",&px_ele, "px_ele/D");
        t2->Branch("py_ele",&py_ele, "py_ele/D");
        t2->Branch("pz_ele",&pz_ele, "pz_ele/D");
        t2->Branch("E0_ele",&E0_ele, "E0_ele/D");
        t2->Branch("px_had",&px_had, "px_had/D");
        t2->Branch("py_had",&py_had, "py_had/D");
        t2->Branch("pz_had",&pz_had, "pz_had/D");
        t2->Branch("E0_had",&E0_had, "E0_had/D");
        t2->Branch("vx_ele",&vx_ele, "vx_ele/D");
        t2->Branch("vy_ele",&vy_ele, "vy_ele/D");
        t2->Branch("vz_ele",&vz_ele, "vz_ele/D");
        t2->Branch("vx_had",&vx_had, "vx_had/D");
        t2->Branch("vy_had",&vy_had, "vy_had/D");
        t2->Branch("vz_had",&vz_had, "vz_had/D");
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
        t3->Branch("theta_q",&theta_q,"data/D");
        t3->Branch("theta_s",&theta_s,"data/D");
        t3->Branch("phi_h",&phi_h,"data/D");
        t3->Branch("phi_s",&phi_s,"data/D");
        t3->Branch("jacoF",&jacoF,"jacoF/D");
        t3->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
        t3->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
        t3->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
        t3->Branch("mom_ele",&mom_ele,"mom_ele/D");
        t3->Branch("mom_had",&mom_had,"mom_had/D");
        t3->Branch("theta_ele",&theta_ele,"theta_ele/D");
        t3->Branch("theta_had",&theta_had,"theta_had/D");
        t3->Branch("phi_ele",&phi_ele,"phi_ele/D");
        t3->Branch("phi_had",&phi_had,"phi_had/D");
        t3->Branch("nsim",&nsim,"nsim/I");
        t3->Branch("dilute_p",&dilute_hp ,"data/D");
        t3->Branch("dilute_m",&dilute_hm ,"data/D");
        t3->Branch("px_ele",&px_ele, "px_ele/D");
        t3->Branch("py_ele",&py_ele, "py_ele/D");
        t3->Branch("pz_ele",&pz_ele, "pz_ele/D");
        t3->Branch("E0_ele",&E0_ele, "E0_ele/D");
        t3->Branch("px_had",&px_had, "px_had/D");
        t3->Branch("py_had",&py_had, "py_had/D");
        t3->Branch("pz_had",&pz_had, "pz_had/D");
        t3->Branch("E0_had",&E0_had, "E0_had/D");
        t3->Branch("vx_ele",&vx_ele, "vx_ele/D");
        t3->Branch("vy_ele",&vy_ele, "vy_ele/D");
        t3->Branch("vz_ele",&vz_ele, "vz_ele/D");
        t3->Branch("vx_had",&vx_had, "vx_had/D");
        t3->Branch("vy_had",&vy_had, "vy_had/D");
        t3->Branch("vz_had",&vz_had, "vz_had/D");
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
        t4->Branch("theta_q",&theta_q,"data/D");
        t4->Branch("theta_s",&theta_s,"data/D");
        t4->Branch("phi_h",&phi_h,"data/D");
        t4->Branch("phi_s",&phi_s,"data/D");
        t4->Branch("jacoF",&jacoF,"jacoF/D");
        t4->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
        t4->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
        t4->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
        t4->Branch("mom_ele",&mom_ele,"mom_ele/D");
        t4->Branch("mom_had",&mom_had,"mom_had/D");
        t4->Branch("theta_ele",&theta_ele,"theta_ele/D");
        t4->Branch("theta_had",&theta_had,"theta_had/D");
        t4->Branch("phi_ele",&phi_ele,"phi_ele/D");
        t4->Branch("phi_had",&phi_had,"phi_had/D");
        t4->Branch("nsim",&nsim,"nsim/I");
        t4->Branch("dilute_p",&dilute_hp ,"data/D");
        t4->Branch("dilute_m",&dilute_hm ,"data/D");
        t4->Branch("px_ele",&px_ele, "px_ele/D");
        t4->Branch("py_ele",&py_ele, "py_ele/D");
        t4->Branch("pz_ele",&pz_ele, "pz_ele/D");
        t4->Branch("E0_ele",&E0_ele, "E0_ele/D");
        t4->Branch("px_had",&px_had, "px_had/D");
        t4->Branch("py_had",&py_had, "py_had/D");
        t4->Branch("pz_had",&pz_had, "pz_had/D");
        t4->Branch("E0_had",&E0_had, "E0_had/D");
        t4->Branch("vx_ele",&vx_ele, "vx_ele/D");
        t4->Branch("vy_ele",&vy_ele, "vy_ele/D");
        t4->Branch("vz_ele",&vz_ele, "vz_ele/D");
        t4->Branch("vx_had",&vx_had, "vx_had/D");
        t4->Branch("vy_had",&vy_had, "vy_had/D");
        t4->Branch("vz_had",&vz_had, "vz_had/D");/*}}}*/
    }

    /*}}}*/

    ofstream pos_gemc, neg_gemc;/*{{{*/
    TString filename_pos = Form("%s_%d_pos.LUND",filename0.Data(), Int_t(FileNo));
    TString filename_neg = Form("%s_%d_neg.LUND",filename0.Data(), Int_t(FileNo));
    if(bLUND){
        pos_gemc.open(filename_pos);
        neg_gemc.open(filename_neg);
    }/*}}}*/

    //Only initialize once here/*{{{*/
    //LHAPDF, CTEQPDF or EPS09
    SIDIS *sidis = new SIDIS(model);

    //1->LO need CTEQ6L1,  2->NLO need CTEQ6.1M, default is 2
    //3->free L0 CTEQ6L1 PDF, 4->free NL0 CTEQ6.1M PDF, default is 4
    //sidis->SetCTEQOrder(2);/*}}}*/

    bool exitcondition=true;	
    while(exitcondition){/*{{{*/
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
                mom_gen_had, theta_gen_had, phi_gen_had,
                ion_mass, A, Z, particle_flag);

        mom_ele = sidis->mom_ele; theta_ele = sidis->theta_ele; phi_ele = sidis->phi_ele;
        mom_had = sidis->mom_had; theta_had = sidis->theta_had; phi_had = sidis->phi_had;
        theta_q=sidis->theta_q;  theta_s=sidis->theta_s; phi_s= sidis->phi_s; phi_h = sidis->phi_h;
        px_ele = sidis->px_ele;	py_ele = sidis->py_ele;	pz_ele = sidis->pz_ele; E0_ele = sidis->E0_ele;
        px_had = sidis->px_had;	py_had = sidis->py_had;	pz_had = sidis->pz_had; E0_had = sidis->E0_had;

        x=sidis->x; y=sidis->y; z=sidis->z; Q2=sidis->Q2; W=sidis->W; Wp=sidis->Wp;
        s=sidis->s; nu=sidis->nu; pt=sidis->pt; gamma=sidis->gamma; epsilon=sidis->epsilon;
        jacoF=sidis->jacoF;/*}}}*/

        /*Get XS{{{*/
        if (x>=0.0&&x<=1.0&&Q2 >=1.0 && W>= 2.0 //&&Wp>= 1.6
                &&( (config=="EIC" && z>0.2&&z<0.9//&&y>0.05&&y<0.8
                        && ((count[0]<number_of_events&&pt<=1.0&&Q2<=10.) 
                            || (count[1]<number_of_events&&pt>1.0&&Q2<=10.)
                            || (count[2]<number_of_events&&pt<=1.0&&Q2>10.)
                            || (count[3]<number_of_events&&pt>1.0&&Q2>10.)) 
                    )
                    ||(config=="SoLID" && z>0.3&&z<0.7 
                        && ((count[0]<number_of_events&&pt<=1.0&&Q2<=10.) 
                            || (count[1]<number_of_events&&pt>1.0&&Q2<=10.))
                      ) 
                    ||(config=="SPECT" && z>0.2&&z<0.9 
                        && count[0]<number_of_events 
                      )
                  )){
            sidis->CalcXS();/*{{{*/
            dxs_incl = sidis->GetXS_Inclusive();
            dxs_hp = sidis->GetXS_HP();
            dxs_hm = sidis->GetXS_HM();
            dilute_hp = sidis->GetDilute_HP();
            dilute_hm = sidis->GetDilute_HM();/*}}}*/

            //to avoid some wired behavior in log scale/*{{{*/
            if((dxs_incl)<1e-16) dxs_incl=1e-16;
            if((dxs_hp)<1e-16) dxs_hp=1e-16;
            if((dxs_hm)<1e-16) dxs_hm=1e-16;
            if((dilute_hp)<1e-16) dilute_hp=1e-16;
            if((dilute_hm)<1e-16) dilute_hm=1e-16;

            if(isnan(dxs_incl)) dxs_incl=1e-16;
            if(isnan(dxs_hp)) dxs_hp=1e-16;
            if(isnan(dxs_hm)) dxs_hm=1e-16;
            if(isnan(dilute_hp)) dilute_hp=1e-16;
            if(isnan(dilute_hm)) dilute_hm=1e-16;
            if(isinf(dxs_incl)) dxs_incl=1e-16;
            if(isinf(dxs_hp)) dxs_hp=1e-16;
            if(isinf(dxs_hm)) dxs_hm=1e-16;
            if(isinf(dilute_hp)) dilute_hp=1e-16;
            if(isinf(dilute_hm)) dilute_hm=1e-16;
            /*}}}*/

            /*LUND For Positive Hadron{{{*/
            //These section will save events based on their XS distributions
            Double_t cdxs_max_rndm_hp = cdxs_max * gRandom->Uniform(0, 1);
            if(bLUND&&dxs_hp < cdxs_max_rndm_hp){
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
                        dxs_hp						
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
                        E0_ele,
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
                        E0_had, 
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
            if(bLUND&&dxs_hm < cdxs_max_rndm_hm){
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
                        dxs_hm						
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
                        E0_ele,
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
                        E0_had, 
                        mass_had, //mass not in used
                        vx_had, //vx
                        vy_had, //vx
                        vz_had  //vx
                        )<<endl;/*}}}*/
            }
            /*}}}*/

            /*Save PI+ ROOT file based on XS distribution{{{*/
            if(bXSMode&&(dxs_hp < cdxs_max_rndm_hp)){
                t1->Fill();

                //Just save events in one root files in this case
                count[0] ++;//cout << 0 << " " << count[0] << endl;
                cout << count[0] <<"\r";
                count[2] = number_of_events;
                count[3] = number_of_events;
            }
            /*}}}*/

            /*Save PI- ROOT file based on XS distribution{{{*/
            if(bXSMode&&(dxs_hm < cdxs_max_rndm_hm)){
                t2->Fill();

                //Just save events in one root files in this case
                count[1] ++;//cout << 0 << " " << count[0] << endl;
                cout << count[1] <<"\r";
                count[2] = number_of_events;
                count[3] = number_of_events;
            }
            /*}}}*/

            if(!bXSMode){/*{{{*/
                if ((dxs_hp)!=0&&(dxs_hm)!=0){
                    if((config=="SPECT")||(Q2<=10.&&pt<=1.0&&(config=="SoLID"||config=="EIC"))){
                        t1->Fill();
                        count[0] ++;//cout << 0 << " " << count[0] << endl;
                    }

                    if(config=="SoLID"||config=="EIC"){
                        if (Q2<=10.&&pt>1.0){
                            t2->Fill();
                        }
                    }

                    if(config=="EIC"){
                        if (Q2>10.&&pt<=1.0){
                            t3->Fill();
                            count[2] ++;//cout << 2 << " " << count[2] << endl;
                        }
                        if (Q2>10.&&pt>1.0){
                            t4->Fill();
                            count[3] ++;//cout << 3 << " " << count[3] <<  endl;
                        }
                    } 
                }
                cout << count[0] << "\t" << count[1] << "\t" << count[2] << "\t" << count[3] << "\r";
                //cout << nsim << endl;
            }
            /*}}}*/
        }

        //judging exitcondition/*{{{*/
        if (config=="EIC") {
            if (count[0] < number_of_events || count[1] < number_of_events 
                    || count[2] < number_of_events || count[3] < number_of_events) exitcondition=true;
            else exitcondition=false;
        } 
        else if (config=="SoLID") {
            if (count[0] < number_of_events || count[1] < number_of_events) exitcondition=true;
            else exitcondition=false;
        } 
        else if(config=="SPECT") {
            if (count[0] < number_of_events) exitcondition=true;
            else exitcondition=false;
        } /*}}}*/
    }/*}}}*/

    cout << count[0] << "\t" << count[1] << "\t" << count[2] << "\t" << count[3] << endl;

    file1->Write();/*{{{*/
    file1->Close();
    if(config=="SoLID" || config=="EIC" ){
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
        double weight_hp=1.0e-16;
        double weight_hm=1.0e-16;
        double weight_in=1.0e-16;
        /*Generate weights for file1{{{*/
        cout<<"--- Now insert weights to Root-file #1"<<filename1.Data()<<endl;
        TFile *f1 = new TFile(filename1.Data(), "update");
        TTree *T1=(TTree*) f1->Get("T");
        Long64_t N1=T1->GetEntries();
        T1->SetBranchAddress("dxs_hp",&dxs_hp);
        T1->SetBranchAddress("dxs_hm",&dxs_hm);
        T1->SetBranchAddress("dxs_incl",&dxs_incl);
        T1->SetBranchAddress("nsim",&nsim);
        T1->GetEntry(N1-1);          //get nsim for this rootfile
        int Nsim1=nsim;

        TBranch *branch_weight_in1=T1->Branch("weight_in",&weight_in,"weight_in/D");
        TBranch *branch_weight_hp1=T1->Branch("weight_hp",&weight_hp,"weight_hp/D");
        TBranch *branch_weight_hm1=T1->Branch("weight_hm",&weight_hm,"weight_hm/D");
        cout<<Form("---Filling weights foor ROOT#1, Nsim=%d, Phase_space = %f", Nsim1, Phase_space)<<endl;
        for(Long64_t i=0;i<N1;i++){
            T1->GetEntry(i);
            //warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
            weight_in=dxs_incl*electron_phase_space/Nsim1;   
            weight_hp=dxs_hp*Phase_space/Nsim1;   
            weight_hm=dxs_hm*Phase_space/Nsim1;

            if((weight_in)<1e-16) weight_in=1e-16;
            if((weight_hp)<1e-16) weight_hp=1e-16;
            if((weight_hm)<1e-16) weight_hm=1e-16;
            if(isnan(weight_in)) weight_in=1e-16;
            if(isnan(weight_hp)) weight_hp=1e-16;
            if(isnan(weight_hm)) weight_hm=1e-16;
            if(isinf(weight_in)) weight_in=1e-16;
            if(isinf(weight_hp)) weight_hp=1e-16;
            if(isinf(weight_hm)) weight_hm=1e-16;


            branch_weight_in1->Fill();
            branch_weight_hp1->Fill();
            branch_weight_hm1->Fill();
        }
        T1->Write("",TObject::kOverwrite);
        f1->Close();
        /*}}}*/

        if(config=="SoLID" || config=="EIC" ){
            /*Generate weights for file2{{{*/
            cout<<"--- Now insert weights to Root-file #2"<<filename2.Data()<<endl;
            TFile *f2 = new TFile(filename2.Data(), "update");
            TTree *T2=(TTree*) f2->Get("T");
            Long64_t N2=T2->GetEntries();
            T2->SetBranchAddress("dxs_incl",&dxs_incl);
            T2->SetBranchAddress("dxs_hp",&dxs_hp);
            T2->SetBranchAddress("dxs_hm",&dxs_hm);
            T2->SetBranchAddress("nsim",&nsim);
            T2->GetEntry(N2-1);          //get nsim for this rootfile
            int Nsim2=nsim;

            TBranch *branch_weight_in2=T2->Branch("weight_in",&weight_in,"weight_in/D");
            TBranch *branch_weight_hp2=T2->Branch("weight_hp",&weight_hp,"weight_hp/D");
            TBranch *branch_weight_hm2=T2->Branch("weight_hm",&weight_hm,"weight_hm/D");
            cout<<Form("---Filling weights foor ROOT#2, Nsim=%d, Phase_space = %f", Nsim2, Phase_space)<<endl;
            for(Long64_t i=0;i<N2;i++){
                T2->GetEntry(i);
                //warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
                weight_in=dxs_incl*electron_phase_space/Nsim2;   
                weight_hp=dxs_hp*Phase_space/Nsim2;   
                weight_hm=dxs_hm*Phase_space/Nsim2;

                if((weight_in)<1e-16) weight_in=1e-16;
                if((weight_hp)<1e-16) weight_hp=1e-16;
                if((weight_hm)<1e-16) weight_hm=1e-16;
                if(isnan(weight_in)) weight_in=1e-16;
                if(isnan(weight_hp)) weight_hp=1e-16;
                if(isnan(weight_hm)) weight_hm=1e-16;
                if(isinf(weight_in)) weight_in=1e-16;
                if(isinf(weight_hp)) weight_hp=1e-16;
                if(isinf(weight_hm)) weight_hm=1e-16;
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
            Long64_t N3=T3->GetEntries();
            T3->SetBranchAddress("dxs_incl",&dxs_incl);
            T3->SetBranchAddress("dxs_hp",&dxs_hp);
            T3->SetBranchAddress("dxs_hm",&dxs_hm);
            T3->SetBranchAddress("nsim",&nsim);
            T3->GetEntry(N3-1);          //get nsim for this rootfile
            int Nsim3=nsim;

            TBranch *branch_weight_in3=T3->Branch("weight_in",&weight_in,"weight_in/D");
            TBranch *branch_weight_hp3=T3->Branch("weight_hp",&weight_hp,"weight_hp/D");
            TBranch *branch_weight_hm3=T3->Branch("weight_hm",&weight_hm,"weight_hm/D");
            cout<<Form("---Filling weights foor ROOT#3, Nsim=%d, Phase_space = %f", Nsim3, Phase_space)<<endl;
            for(Long64_t i=0;i<N3;i++){
                T3->GetEntry(i);
                //warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
                weight_in=dxs_incl*electron_phase_space/Nsim3;   
                weight_hp=dxs_hp*Phase_space/Nsim3;   
                weight_hm=dxs_hm*Phase_space/Nsim3;

                if((weight_in)<1e-16) weight_in=1e-16;
                if((weight_hp)<1e-16) weight_hp=1e-16;
                if((weight_hm)<1e-16) weight_hm=1e-16;
                if(isnan(weight_in)) weight_in=1e-16;
                if(isnan(weight_hp)) weight_hp=1e-16;
                if(isnan(weight_hm)) weight_hm=1e-16;
                if(isinf(weight_in)) weight_in=1e-16;
                if(isinf(weight_hp)) weight_hp=1e-16;
                if(isinf(weight_hm)) weight_hm=1e-16;
                
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
            Long64_t N4=T4->GetEntries();
            T4->SetBranchAddress("dxs_incl",&dxs_incl);
            T4->SetBranchAddress("dxs_hp",&dxs_hp);
            T4->SetBranchAddress("dxs_hm",&dxs_hm);
            T4->SetBranchAddress("nsim",&nsim);
            T4->GetEntry(N4-1);          //get nsim for this rootfile
            int Nsim4=nsim;

            TBranch *branch_weight_in4=T4->Branch("weight_in",&weight_in,"weight_in/D");
            TBranch *branch_weight_hp4=T4->Branch("weight_hp",&weight_hp,"weight_hp/D");
            TBranch *branch_weight_hm4=T4->Branch("weight_hm",&weight_hm,"weight_hm/D");
            cout<<Form("---Filling weights foor ROOT#4, Nsim=%d, Phase_space = %f", Nsim4, Phase_space)<<endl;
            for(Long64_t i=0;i<N4;i++){
                T4->GetEntry(i);
                //warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
                weight_in=dxs_incl*electron_phase_space/Nsim4;   
                weight_hp=dxs_hp*Phase_space/Nsim4;   
                weight_hm=dxs_hm*Phase_space/Nsim4;

                if((weight_in)<1e-16) weight_in=1e-16;
                if((weight_hp)<1e-16) weight_hp=1e-16;
                if((weight_hm)<1e-16) weight_hm=1e-16;
                if(isnan(weight_in)) weight_in=1e-16;
                if(isnan(weight_hp)) weight_hp=1e-16;
                if(isnan(weight_hm)) weight_hm=1e-16;
                if(isinf(weight_in)) weight_in=1e-16;
                if(isinf(weight_hp)) weight_hp=1e-16;
                if(isinf(weight_hm)) weight_hm=1e-16;

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
        if(A==1 && Z==1)//Hydrogen
            ion_mass = mass_p;
        else if(A==1 && Z==0)//Neutron
            ion_mass = mass_n;
        else if(A==2 && Z==1)//Deutron
            ion_mass = mass_u*2.014102;
        else if(A==3 && Z==1)//Tritium
            ion_mass = mass_u*3.016049;
        else if(A==3 && Z==2)//He3
            ion_mass = mass_u*3.016029;
        else if(A==4 && Z==2)//He4
            ion_mass = mass_u*4.002602;
        else if(A==7 && Z==3)//Lithium
            ion_mass = mass_u*6.941000;
        else if(A==9 && Z==4)//Be9
            ion_mass = mass_u*9.012182;
        else if(A==10 && Z==5)//Boron
            ion_mass = mass_u*10.811000;
        else if(A==12 && Z==6)//Carbon
            ion_mass = mass_u*12.010700;
        else if(A==14 && Z==7)//Nitrogen
            ion_mass = mass_u*14.0067;
        else if(A==16 && Z==8)//Oxygen
            ion_mass = mass_u*15.9994;
        else if(A==27 && Z==13)//Al
            ion_mass = mass_u*26.981539;
        else if(A==40 && Z==20)//C40
            ion_mass = mass_u*40.07800;
        else if(A==48 && Z==20)//C40
            ion_mass = mass_u*47.952534;
        else if(A==56 && Z==26)//Fe
            ion_mass = mass_u*55.84500;
        else if(A==64 && Z==29)//Cu
            ion_mass = mass_u*63.546;
        else if(A==197 && Z==79)//Cu
            ion_mass = mass_u*196.966569;
        //if the target is not in the list above, simply add the total mass of protons and neutrons
        //but in general you can look at the PERIODIC TABLE for new target you want to add
        else
            ion_mass = Z*mass_p+(A-Z)*mass_n;

        return ion_mass;
    }
    /*}}}*/
