#include "iostream.h"
#include "fstream.h"
#include "LHAPDF/LHAPDF.h"
#include "CrossSection.h"
#include "math.h"
#include "TString.h"


using namespace std;
using namespace LHAPDF;

int main(int argc, char *argv[]){
  if (argc != 3) {
    cout << "usage: kine_name outfile_name" << endl;
    return 0;
  }else{
    const int SUBSET = 0;
    const string NAME = "MRST2004nlo";
    
    TString infile_name = argv[1];
    TString outfile_name = argv[2];
    
    setPDFPath("/home/xqian/lhapdf/bin/");
    cout << pdfsetsPath() << endl;
    
    LHAPDF::initPDFSet(NAME, LHAPDF::LHGRID, SUBSET);
    
    ifstream infile(infile_name);
    ofstream outfile(outfile_name);
    
    double logx,x,Q2,z,pt,s,temp;
    while(!infile.eof()){
      infile >> x >> Q2 >> z >> pt >> s;
//       x = pow(10,logx);
      
      outfile << x << " \t" << Cro_Sec_Unp(x,Q2,z,pt,s,1) << " \t" << 
	Cro_Sec_Tran(x,Q2,z,pt,s,1) << " \t" <<
	Cro_Sec_Siv(x,Q2,z,pt,s,1) << " \t" <<
	Cro_Sec_Pret(x,Q2,z,pt,s,1) << " \t" << 
	Cro_Sec_Unp(x,Q2,z,pt,s,-1) << " \t" << 
	Cro_Sec_Tran(x,Q2,z,pt,s,-1) << " \t" <<
	Cro_Sec_Siv(x,Q2,z,pt,s,-1) << " \t" <<
	Cro_Sec_Pret(x,Q2,z,pt,s,-1) << " \t" << endl;
    }
    infile.close();
    outfile.close();
  }
  return 0;
}
