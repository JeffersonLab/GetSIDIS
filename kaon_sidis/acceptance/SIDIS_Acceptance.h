class SIDIS_Acceptance{

    public:

        SIDIS_Acceptance(){ 
            //Init("He3"); 
        }

        virtual ~SIDIS_Acceptance(){
           file_ele->Close();
           file_kp->Close();
           file_km->Close();
           file_pip->Close();
           file_pim->Close();
           file_pho->Close();
        }

        void Init(TString kTarget_Name){/*{{{*/
            TString kDate="";
            if(kTarget_Name=="He3") kDate = "201701";
            if(kTarget_Name=="NH3") kDate = "201710";
            //CLEO
            const TString Acceptance_DIR = "/work/halla/solid/yez/sidis/acceptance";
            TString ele_filename = Form("%s/SIDIS_%s/acceptance_solid_SIDIS_%s_electron_%s_1e7_output_final.root", Acceptance_DIR.Data(),kTarget_Name.Data(),kTarget_Name.Data(),kDate.Data());
            TString pip_filename = Form("%s/SIDIS_%s/acceptance_solid_SIDIS_%s_pip_%s_1e7_output_final.root",      Acceptance_DIR.Data(),kTarget_Name.Data(),kTarget_Name.Data(),kDate.Data());
            TString pim_filename = Form("%s/SIDIS_%s/acceptance_solid_SIDIS_%s_pim_%s_1e7_output_final.root",      Acceptance_DIR.Data(),kTarget_Name.Data(),kTarget_Name.Data(),kDate.Data());
            TString kp_filename  = Form("%s/SIDIS_%s/acceptance_solid_SIDIS_%s_kp_%s_1e7_output_final.root",       Acceptance_DIR.Data(),kTarget_Name.Data(),kTarget_Name.Data(),kDate.Data());
            TString km_filename  = Form("%s/SIDIS_%s/acceptance_solid_SIDIS_%s_km_%s_1e7_output_final.root",       Acceptance_DIR.Data(),kTarget_Name.Data(),kTarget_Name.Data(),kDate.Data());
            TString pho_filename = Form("%s/SIDIS_%s/acceptance_solid_SIDIS_%s_photon_%s_1e7_output_final.root",   Acceptance_DIR.Data(),kTarget_Name.Data(),kTarget_Name.Data(),kDate.Data());

            cout<<endl<<"--- Loading the following SoLID-SIDIS Acceptance Profile: "<<endl;
            cout<<"     e-: "<<ele_filename.Data()<<endl;
            cout<<"    pi+: "<<pip_filename.Data()<<endl;
            cout<<"    pi-: "<<pim_filename.Data()<<endl;
            cout<<"     K+: "<<kp_filename.Data()<<endl;
            cout<<"     K-: "<<km_filename.Data()<<endl;
            cout<<"    pho: "<<pho_filename.Data()<<endl<<endl;

            TCanvas * dummycavans=new TCanvas();
            file_ele=new TFile(ele_filename.Data(),"r");
            accep_ele_forward=(TH2D*)file_ele->Get("acceptance_ThetaP_forwardangle"); 
            accep_ele_large=(TH2D*)file_ele->Get("acceptance_ThetaP_largeangle");    
            accep_ele_forward_3D=(TH3D*)file_ele->Get("acceptance_ThetaPhiP_forwardangle"); 
            accep_ele_large_3D=(TH3D*)file_ele->Get("acceptance_ThetaPhiP_largeangle");    


            file_kp=new TFile(kp_filename.Data(),"r");
            accep_kp_forward=(TH2D*)file_kp->Get("acceptance_ThetaP_forwardangle"); 
            accep_kp_large=(TH2D*)file_kp->Get("acceptance_ThetaP_largeangle");    
            accep_kp_forward_3D=(TH3D*)file_kp->Get("acceptance_ThetaPhiP_forwardangle"); 
            accep_kp_large_3D=(TH3D*)file_kp->Get("acceptance_ThetaPhiP_largeangle");    

            file_km=new TFile(km_filename.Data(),"r");
            accep_km_forward=(TH2D*)file_km->Get("acceptance_ThetaP_forwardangle");
            accep_km_large=(TH2D*)file_km->Get("acceptance_ThetaP_largeangle");   
            accep_km_forward_3D=(TH3D*)file_km->Get("acceptance_ThetaPhiP_forwardangle"); 
            accep_km_large_3D=(TH3D*)file_km->Get("acceptance_ThetaPhiP_largeangle");    
 
            file_pip=new TFile(pip_filename.Data(),"r");
            accep_pip_forward=(TH2D*)file_pip->Get("acceptance_ThetaP_forwardangle");
            accep_pip_large=(TH2D*)file_pip->Get("acceptance_ThetaP_largeangle");   
            accep_pip_forward_3D=(TH3D*)file_pip->Get("acceptance_ThetaPhiP_forwardangle"); 
            accep_pip_large_3D=(TH3D*)file_pip->Get("acceptance_ThetaPhiP_largeangle");    
    
            file_pim=new TFile(pim_filename.Data(),"r");
            accep_pim_forward=(TH2D*)file_pim->Get("acceptance_ThetaP_forwardangle");
            accep_pim_large=(TH2D*)file_pim->Get("acceptance_ThetaP_largeangle");   
            accep_pim_forward_3D=(TH3D*)file_pim->Get("acceptance_ThetaPhiP_forwardangle"); 
            accep_pim_large_3D=(TH3D*)file_pim->Get("acceptance_ThetaPhiP_largeangle");    
 
            file_pho=new TFile(pho_filename.Data(),"r");
            accep_pho_forward=(TH2D*)file_pho->Get("acceptance_ThetaP_forwardangle");
            accep_pho_large=(TH2D*)file_pho->Get("acceptance_ThetaP_largeangle");   
            accep_pho_forward_3D=(TH3D*)file_pho->Get("acceptance_ThetaPhiP_forwardangle"); 
            accep_pho_large_3D=(TH3D*)file_pho->Get("acceptance_ThetaPhiP_largeangle");    
        }/*}}}*/

        double GetAcc_OLD(TString kPID, TString kRegion, double kP, double kTheta){/*{{{*/
            //sidis_he3: theta=(0,40 degrees), phi = (-180, 180 degrees), p=(0,11 gev/c) 
            //sidis_nh3: theta=(0,50 degrees), phi = (-180, 180 degrees), p=(0,11 gev/c) 

            //old bin before 2017
            //int theta_bin=int(kTheta/0.2)+1;    //0.2 degree per bin
			//int p_bin=int(kP/0.05)+1;      //0.05 gev per bin for mom

            //new bin before 2017
            int theta_bin=int(kTheta/0.5)+1;    //0.5 degree per bin
			int p_bin=int(kP/0.01)+1;      //0.1 gev per bin for mom

            double acc=0.0;
            if(kRegion=="FA"||kRegion=="fa"||kRegion=="Forward"||kRegion=="forward"){
                if(kPID=="e-"||kPID=="ele"||kPID=="electron")
                    acc=accep_ele_forward->GetBinContent(theta_bin,p_bin);
                else if(kPID=="neutron"||kPID=="photon"||kPID=="pho"||kPID=="n"||kPID=="g")
                    acc=accep_pho_forward->GetBinContent(theta_bin,p_bin);
                else if(kPID=="pip"||kPID=="pionp"||kPID=="pi+")
                    acc=accep_pip_forward->GetBinContent(theta_bin,p_bin);
                else if(kPID=="pim"||kPID=="pionm"||kPID=="pi-")
                    acc=accep_pim_forward->GetBinContent(theta_bin,p_bin);
                else if(kPID=="kp"||kPID=="kaonp"||kPID=="K+")
                    acc=accep_kp_forward->GetBinContent(theta_bin,p_bin);
                else if(kPID=="kp"||kPID=="kaonm"||kPID=="K-")
                    acc=accep_km_forward->GetBinContent(theta_bin,p_bin);
                else{
                    cerr<<"*** ERROR, I don't know this particle = "<<kPID.Data()<<endl;
                    return 0.0;
                }
            }
            else if(kRegion=="LA"||kRegion=="la"||kRegion=="large"||kRegion=="Large"){
                if(kPID=="e-"||kPID=="ele"||kPID=="electron")
                    acc=accep_ele_large->GetBinContent(theta_bin,p_bin);
                else if(kPID=="neutron"||kPID=="photon"||kPID=="pho"||kPID=="n"||kPID=="g")
                    acc=accep_pho_large->GetBinContent(theta_bin,p_bin);
                else if(kPID=="pip"||kPID=="pionp"||kPID=="pi+")
                    acc=accep_pip_large->GetBinContent(theta_bin,p_bin);
                else if(kPID=="pim"||kPID=="pionm"||kPID=="pi-")
                    acc=accep_pim_large->GetBinContent(theta_bin,p_bin);
                else if(kPID=="kp"||kPID=="kaonp"||kPID=="K+")
                    acc=accep_kp_large->GetBinContent(theta_bin,p_bin);
                else if(kPID=="km"||kPID=="kaonm"||kPID=="K-")
                    acc=accep_km_large->GetBinContent(theta_bin,p_bin);
                else{
                    cerr<<"*** ERROR, I don't know this particle = "<<kPID.Data()<<endl;
                    return 0.0;
                }
            }
            else{
                cerr<<"*** ERROR, I don't know this region = "<<kRegion.Data()<< ", or particle = "<<kPID.Data()<<endl;
                return 0.0;
            }
            //cout<<Form("--- Region=%s, PID=%s, P=%f(%d), Theta=%f(%d), Acc = %f", kRegion.Data(), kPID.Data(),kP,p_bin, kTheta,theta_bin,acc)<<endl;

            return acc;
        }/*}}}*/

        double GetAcc(TString kPID, TString kRegion, double kP, double kTheta){/*{{{*/

            double acc=0.0;
            int bin =0, binx=0,biny=0,binz=0;
            if(kRegion=="FA"||kRegion=="fa"||kRegion=="Forward"||kRegion=="forward"){
                if(kPID=="e-"||kPID=="ele"||kPID=="electron")
                {
                    bin =accep_ele_forward->FindBin(kTheta, kP); 
                    accep_ele_forward->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_ele_forward->GetBinContent(binx, biny);
                }
                else if(kPID=="neutron"||kPID=="photon"||kPID=="pho"||kPID=="n"||kPID=="g")
                {
                    bin =accep_pho_forward->FindBin(kTheta, kP); 
                    accep_pho_forward->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pho_forward->GetBinContent(binx, biny);
                }
                else if(kPID=="pip"||kPID=="pionp"||kPID=="pi+")
                {
                    bin =accep_pip_forward->FindBin(kTheta, kP); 
                    accep_pip_forward->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pip_forward->GetBinContent(binx, biny);
                }
                else if(kPID=="pim"||kPID=="pionm"||kPID=="pi-")
                {
                    bin =accep_pim_forward->FindBin(kTheta, kP); 
                    accep_pim_forward->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pim_forward->GetBinContent(binx, biny);
                }
                else if(kPID=="kp"||kPID=="kaonp"||kPID=="K+")
                {
                    bin =accep_kp_forward->FindBin(kTheta, kP); 
                    accep_kp_forward->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_kp_forward->GetBinContent(binx, biny);
                }
                else if(kPID=="kp"||kPID=="kaonm"||kPID=="K-")
                {
                    bin =accep_kp_forward->FindBin(kTheta, kP); 
                    accep_km_forward->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_km_forward->GetBinContent(binx, biny);
                }
                else{
                    cerr<<"*** ERROR, I don't know this particle = "<<kPID.Data()<<endl;
                    return 0.0;
                }
            }
            else if(kRegion=="LA"||kRegion=="la"||kRegion=="large"||kRegion=="Large"){
                if(kPID=="e-"||kPID=="ele"||kPID=="electron")
                {
                    bin =accep_ele_large->FindBin(kTheta, kP); 
                    accep_ele_large->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_ele_large->GetBinContent(binx, biny);
                }
                else if(kPID=="neutron"||kPID=="photon"||kPID=="pho"||kPID=="n"||kPID=="g")
                {
                    bin =accep_pho_large->FindBin(kTheta, kP); 
                    accep_pho_large->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pho_large->GetBinContent(binx, biny);
                }
                else if(kPID=="pip"||kPID=="pionp"||kPID=="pi+")
                {
                    bin =accep_pip_large->FindBin(kTheta, kP); 
                    accep_pip_large->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pip_large->GetBinContent(binx, biny);
                }
                else if(kPID=="pim"||kPID=="pionm"||kPID=="pi-")
                {
                    bin =accep_pim_large->FindBin(kTheta, kP); 
                    accep_pim_large->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pim_large->GetBinContent(binx, biny);
                }
                else if(kPID=="kp"||kPID=="kaonp"||kPID=="K+")
                {
                    bin =accep_kp_large->FindBin(kTheta, kP); 
                    accep_kp_large->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_kp_large->GetBinContent(binx, biny);
                }
                else if(kPID=="km"||kPID=="kaonm"||kPID=="K-")
                {
                    bin =accep_km_large->FindBin(kTheta, kP); 
                    accep_km_large->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_km_large->GetBinContent(binx, biny);
                }
                else{
                    cerr<<"*** ERROR, I don't know this particle = "<<kPID.Data()<<endl;
                    return 0.0;
                }
            }
            else{
                cerr<<"*** ERROR, I don't know this region = "<<kRegion.Data()<< ", or particle = "<<kPID.Data()<<endl;
                return 0.0;
            }
            //cout<<Form("--- Region=%s, PID=%s, P=%f(%d), Theta=%f(%d), Acc = %f", kRegion.Data(), kPID.Data(),kP,p_bin, kTheta,theta_bin,acc)<<endl;

            return acc;
        }/*}}}*/
        
        double GetAcc3D(TString kPID, TString kRegion, double kP, double kTheta, double kPhi){/*{{{*/
            
            double acc=0.0;
            int bin =0, binx=0,biny=0,binz=0;
            if(kRegion=="FA"||kRegion=="fa"||kRegion=="Forward"||kRegion=="forward"){
                if(kPID=="e-"||kPID=="ele"||kPID=="electron")
                {
                    bin =accep_ele_forward_3D->FindBin(kTheta, kPhi, kP); 
                    accep_ele_forward_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_ele_forward_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="neutron"||kPID=="photon"||kPID=="pho"||kPID=="n"||kPID=="g")
                {
                    bin =accep_pho_forward_3D->FindBin(kTheta, kPhi, kP); 
                    accep_pho_forward_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pho_forward_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="pip"||kPID=="pionp"||kPID=="pi+")
                {
                    bin =accep_pip_forward_3D->FindBin(kTheta, kPhi, kP); 
                    accep_pip_forward_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pip_forward_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="pim"||kPID=="pionm"||kPID=="pi-")
                {
                    bin =accep_pim_forward_3D->FindBin(kTheta, kPhi, kP); 
                    accep_pim_forward_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pim_forward_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="kp"||kPID=="kaonp"||kPID=="K+")
                {
                    bin =accep_kp_forward_3D->FindBin(kTheta, kPhi, kP); 
                    accep_kp_forward_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_kp_forward_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="kp"||kPID=="kaonm"||kPID=="K-")
                {
                    bin =accep_km_forward_3D->FindBin(kTheta, kPhi, kP); 
                    accep_km_forward_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_km_forward_3D->GetBinContent(binx, biny, binz);
                }
                else{
                    cerr<<"*** ERROR, I don't know this particle = "<<kPID.Data()<<endl;
                    return 0.0;
                }
            }
            else if(kRegion=="LA"||kRegion=="la"||kRegion=="large"||kRegion=="Large"){
                if(kPID=="e-"||kPID=="ele"||kPID=="electron")
                {
                    bin =accep_ele_large_3D->FindBin(kTheta, kPhi, kP); 
                    accep_ele_large_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_ele_large_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="neutron"||kPID=="photon"||kPID=="pho"||kPID=="n"||kPID=="g")
                {
                    bin =accep_pho_large_3D->FindBin(kTheta, kPhi, kP); 
                    accep_pho_large_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pho_large_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="pip"||kPID=="pionp"||kPID=="pi+")
                {
                    bin =accep_pip_large_3D->FindBin(kTheta, kPhi, kP); 
                    accep_pip_large_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pip_large_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="pim"||kPID=="pionm"||kPID=="pi-")
                {
                    bin =accep_pim_large_3D->FindBin(kTheta, kPhi, kP); 
                    accep_pim_large_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_pim_large_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="kp"||kPID=="kaonp"||kPID=="K+")
                {
                    bin =accep_kp_large_3D->FindBin(kTheta, kPhi, kP); 
                    accep_kp_large_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_kp_large_3D->GetBinContent(binx, biny, binz);
                }
                else if(kPID=="km"||kPID=="kaonm"||kPID=="K-")
                {
                    bin =accep_km_large_3D->FindBin(kTheta, kPhi, kP); 
                    accep_km_large_3D->GetBinXYZ(bin, binx, biny, binz); 
                    acc=accep_km_large_3D->GetBinContent(binx, biny, binz);
                }
                else{
                    cerr<<"*** ERROR, I don't know this particle = "<<kPID.Data()<<endl;
                    return 0.0;
                }
            }
            else{
                cerr<<"*** ERROR, I don't know this region = "<<kRegion.Data()<< ", or particle = "<<kPID.Data()<<endl;
                return 0.0;
            }

            return acc;
        }/*}}}*/

    private:
        /*{{{*/
        TFile *file_ele;
        TH2D *accep_ele_forward;
        TH2D *accep_ele_large;  
        TH3D *accep_ele_forward_3D;
        TH3D *accep_ele_large_3D;  

        TFile *file_kp;
        TH2D *accep_kp_forward;
        TH2D *accep_kp_large;  
        TH3D *accep_kp_forward_3D;
        TH3D *accep_kp_large_3D;  

        TFile *file_km;
        TH2D *accep_km_forward;
        TH2D *accep_km_large; 
        TH3D *accep_km_forward_3D;
        TH3D *accep_km_large_3D;  

        TFile *file_pip;
        TH2D *accep_pip_forward;
        TH2D *accep_pip_large; 
        TH3D *accep_pip_forward_3D;
        TH3D *accep_pip_large_3D;  

        TFile *file_pim;
        TH2D *accep_pim_forward;
        TH2D *accep_pim_large; 
        TH3D *accep_pim_forward_3D;
        TH3D *accep_pim_large_3D;  

        TFile *file_pho;
        TH2D *accep_pho_forward;
        TH2D *accep_pho_large; 
        TH3D *accep_pho_forward_3D;
        TH3D *accep_pho_large_3D;  
        /*}}}*/
};
