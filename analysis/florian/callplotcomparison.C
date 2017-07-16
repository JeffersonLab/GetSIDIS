void callplotcomparison() {

      Int_t filechoices[7] = {0,1,2,3,4,12,13};
    //Int_t cutchoices[7] = {0,1,2,3,11,12,13};
      Int_t cutchoices[2] = {3,13};

    for (int i = 0; i<7; i++) {
     for (int j = 0; j<2; j++) {
    //  gROOT->ProcessLine(Form(".x plotcomparison.C(%i, %i)",filechoices[i],cutchoices[j]));
      cout << Form(".x plotcomparison.C(%i, %i)",filechoices[i],cutchoices[j]) << endl;
     }
    }
gROOT->ProcessLine(".x plotcomparison.C(1,12)");
gROOT->ProcessLine(".x plotcomparison.C(1,13)");
    // gROOT->ProcessLine(".x plotcomparison.C(1,1)");

}
