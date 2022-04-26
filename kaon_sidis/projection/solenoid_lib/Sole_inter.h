#ifndef ROOT_sole_inter
#define ROOT_sole_inter

class sole_inter{

	public:
		sole_inter(){}
		virtual ~sole_inter(){}
		sole_inter(sole_inter const&);
		sole_inter& operator=(sole_inter const&){};

		//sole_inter(): Nbins(360),Ntotal(0.),Bsize(1.){};

		//initialize the number of bins and evnets of the histogram
		void init(Int_t nbins, Double_t* events){ //all phis(degree)
			Ntotal = 0.;
			Nbins = nbins;
			Bsize = 360./Nbins;
			if (Nbins>1000){
				cout << "Something wrong with Nbins !!" << endl;
				cout << "Force it to be 1000 !!" << endl;
				Nbins = 1000;
			}
			// if (Nbins%2!=0){
			//     cout << "Nbins must be 2*n!" << endl;
			//   }
			for (Int_t i=0;i<Nbins;i++){
				angle[i] = Bsize*i+Bsize/2.;
				x[i] = (*(events + i) + *(events+Nbins-1-i))/2.;
			}
		}

		//change the histogram, to remove all the 0 bins in the middle.
		void connect(){
			//current bin position
			Double_t current_angle;
			Int_t min_bin,max_bin;
			Int_t flag = 0;
			x[Int_t(Nbins/2.)] = (x[Int_t(Nbins/2.)]+x[Int_t(Nbins/2.)-1]+x[Int_t(Nbins/2.)+1])/3.;
			x[Int_t(Nbins/2.)-1] = x[Int_t(Nbins/2.)];
			x[Int_t(Nbins/2.)+1] = x[Int_t(Nbins/2.)];
			for (Int_t i=0;i!=Nbins;i++){
				current_angle = angle[i];
				flag = 0;
				if (x[i]==0){
					min_bin = findmin(i,current_angle);
					max_bin = findmax(i,current_angle);
					adjust(i,current_angle,min_bin,max_bin);
				}
			}

			for(Int_t i=0;i!=Nbins;i++){
				xp[i] = (x[i]+x[Nbins-1-i])/2.;
			}
			for(Int_t i=0;i!=Nbins;i++){
				x[i] = xp[i];
				xp[i] = 0.;
			}


			generate();
		}

		void get_new(Double_t *y){
			for(Int_t i=0;i!=Nbins;i++){
				*(y+i) = x[i];
			}
		}

		void get_limit(Double_t *y){
			for (Int_t i=0;i!=4;i++){
				*(y+i) = anglelimit[i];
			}
		}

		Double_t get_accep(Double_t x){ //x is angle(degree)

			Double_t acpt;
			Int_t count,count_min,count_max;
			count = Int_t(x/Bsize);
			Double_t mindis = 10000.;

			if (count ==0){
				count_min = 0;
				count_max = count+1;
			}else if (count==Nbins-1){
				count_min = count-1;
				count_max = Nbins-1;
			}else{
				count_min = count-1;
				count_max = count+1;
			}

			for (Int_t i=count_min;i<=count_max;i++){
				if (fabs(x-angle[i])<mindis){
					count = i;
					mindis = fabs(x-angle[i]);
				}
			}

			acpt = xp[count];
			return acpt;
		}

	private:
		void generate(){
			// count total number of events
			for(Int_t i=0;i<Nbins;i++){
				Ntotal += x[i];
			}
			find_edge();
			if (Ntotal!=0){
				Double_t sum = 0;
				Double_t count = 0;
				for(Int_t i=0;i<Nbins;i++){
					xp[i] = x[i]/( Ntotal*Bsize/(anglelimit[3]-anglelimit[2]+anglelimit[1]-anglelimit[0]));
					if (xp[i]!=0){
						sum+=xp[i];
						count++;
					}
				} 
				for(Int_t i=0;i<Nbins;i++){
					xp[i] = xp[i]/sum*count;
				}
			}else{
				for(Int_t i=0;i<Nbins;i++){
					xp[i] = x[i];
				}
			}
		}	

		// if one bin has 0 counts, find the bin before  it which is not 0
		Int_t findmin(Int_t current_bin,Double_t current_angle){
			Int_t minbin=current_bin;
			if (current_angle>180.){
				minbin = Int_t(Nbins/2.);
				for (Int_t i=current_bin;i>=Int_t(Nbins/2.);i--){
					if (x[i]!=0){
						minbin = i;
						break;
					}
				}
			}else{
				minbin = 0;
				for (Int_t i=current_bin;i>=0;i--){
					if (x[i]!=0){
						minbin = i;
						break;
					}
				}
			}
			return minbin;
		}

		// if one bin has 0 counts, find the bin after  it which is not 0
		Int_t findmax(Int_t current_bin, Double_t current_angle){
			Int_t maxbin=current_bin;
			if (current_angle>180.){
				maxbin = Nbins-1;
				for (Int_t i=current_bin;i<=Nbins-1;i++){
					if (x[i]!=0){
						maxbin = i;
						break;
					}
				}
			}else{
				maxbin = Int_t(Nbins/2.);
				for (Int_t i=current_bin;i<=Int_t(Nbins/2.);i++){
					if (x[i]!=0){
						maxbin = i;
						break;
					}
				}
			}
			return maxbin;
		}

		// fill in the events between the minbin and maxbin which is 0
		void adjust(Int_t current_bin,Double_t current_angle,Int_t minbin,Int_t maxbin){
			if ((current_angle>180.&&maxbin!=Nbins-1&&minbin!=Int_t(Nbins/2.))||(current_angle<180.&&minbin!=0&&maxbin!=Int_t(Nbins/2.))){
				Double_t a,b;
				Double_t sum = x[maxbin]+x[minbin];
				Double_t asy = (x[maxbin]-x[minbin])/sum;
				// cout << x[maxbin]+x[minbin] << " \t";
				if (asy==0){
					b = 0;
					a = sum/(maxbin-minbin+1);
				}else{
					b = sum/((maxbin+minbin)*(maxbin-minbin+1)/2.+(maxbin-minbin+1)*(asy*(maxbin+minbin)-(maxbin-minbin))/2./asy);
					a = b*(asy*(maxbin+minbin)-(maxbin-minbin))/2./asy;
				}
				Double_t temp_sum = 0;
				for (Int_t i=minbin;i!=maxbin+1;i++){
					x[i] = fabs(a+b*i);
					temp_sum+=x[i];
				}
			}
		}

		void find_edge(){
			binlimit[0] = 0; 
			binlimit[1] = Int_t(Nbins/2.);
			binlimit[2] = Int_t(Nbins/2.);
			binlimit[3] = Nbins-1;
			for (Int_t i=binlimit[0];i<=Int_t(Nbins/2.);i++){
				//cout << x[i] << endl;
				if (x[i]!=0) {
					binlimit[0] = i;
					break;
				}
			}
			for (Int_t i=binlimit[2];i<=Nbins-1;i++){
				if (x[i]!=0) {
					binlimit[2] = i;
					break;
				}
			}
			for (Int_t i=binlimit[1];i>=0;i--){
				if (x[i]!=0) {
					binlimit[1] = i;
					break;
				}
			}
			for (Int_t i=binlimit[3];i>=Int_t(Nbins/2.);i--){
				if (x[i]!=0) {
					binlimit[3] = i;
					break;
				}
			}

			anglelimit[0] = angle[binlimit[0]]-Bsize/2.;
			anglelimit[1] = angle[binlimit[1]]+Bsize/2.;
			anglelimit[2] = angle[binlimit[2]]-Bsize/2.;
			anglelimit[3] = angle[binlimit[3]]+Bsize/2.;

			if (anglelimit[2]<anglelimit[1]){
				anglelimit[1] = (anglelimit[2]+anglelimit[1])/2.;
				anglelimit[2] = anglelimit[1];
			}
		}


	private:
		Int_t binlimit[4];
		Double_t anglelimit[4];
		Double_t x[3600];
		Double_t xp[3600];
		Double_t angle[3600];
		Int_t Nbins; 
		Double_t Ntotal; 
		Double_t Bsize;
};
#endif

