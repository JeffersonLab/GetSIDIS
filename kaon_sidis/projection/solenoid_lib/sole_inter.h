#ifndef ROOT_sole_inter
#define ROOT_sole_inter

class sole_inter 
{
 public:
  sole_inter(): Nbins(360),Ntotal(0.),Bsize(1.){};
  virtual ~sole_inter();
  void init(Int_t, Double_t* );
  void connect();
  Double_t get_accep(Double_t xb);
  void get_new(Double_t* );
  void get_limit(Double_t*);

 protected:
  void generate();
  Int_t findmin(Int_t current_bin,Double_t current_angle);
  Int_t findmax(Int_t current_bin,Double_t current_angle);
  void adjust(Int_t current_bin,Double_t current_angle,Int_t minbin,Int_t maxbin);
  
  void find_edge();
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

